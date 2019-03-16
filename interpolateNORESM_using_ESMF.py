
from __future__ import print_function
import numpy as np
import matplotlib
matplotlib.use('Agg')
import string, os, sys
from netCDF4 import Dataset
import netCDF4
from datetime import datetime
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from pylab import *
import datetime
import glob
from shutil import copyfile

try:
    import ESMF
except:
    print("Unable to continue without ESMF: install using : conda install -c conda-forge esmpy")
    sys.exit()

__author__   = 'Trond Kristiansen'
__email__    = 'me (at) trondkristiansen.com'
__created__  = datetime.datetime(2019, 1, 21)
__modified__ = datetime.datetime(2019, 1, 25)
__version__  = "1.0"
__status__   = "Development, 21.1.2019, 25.01.2019, 16.03.2019"

# Prior to running this script you need to do the following.
#
# 1. Create one file containing all of the atmospheric variables from the NorESM model
# 2. In our case, CLDTOT came in a separate file which then had to be merged with the other atmospheric variables:
#    => cdo merge CLDTOT.cam2.hmlvl.2006-2100.nc atm.cam2.hmlvl.2006-2100.nc NORESM_ATM.cam2.hmlvl.2006-2100.nc
#
# 3. Create the detrended climatology from hindcast timeseries (ERA dataset):
#    => createClimatologyERA.sh
# 4. Create the detrended climatology and remove from teh timeseries to create residuals/deltas:
#    => createDeltasNorESM.sh
# 
# The results from these steps are stored in the results file and used in this script
#
# Copy to Fram:
# scp interpolateNORESM_using_ESMF.py trondk@fram.sigma2.no:/cluster/projects/nn9412k/A20/DELTA/.
#
# Login to fram: 
# module load Anaconda2/5.0.1
# source activate OpenDrift
# python interpolateNORESM_using_ESMF.py 
def createContourPlot(atmdata,lon,lat,basepath,varname,currentdate,filetype):
    plt.clf()
    mymap = Basemap(projection='npstere', boundinglat=45,lon_0=0)
    x, y = mymap(lon,lat)
    mymap.drawcoastlines()
    mymap.fillcontinents(color='grey',zorder=2)
    mymap.drawcountries()
    mymap.drawmapboundary()

    delta=(atmdata.max()-atmdata.min())/20
    levels=np.arange(atmdata.min(),atmdata.max(),delta)
    
    CS1 = mymap.contourf(x,y,atmdata,levels,cmap=plt.cm.jet, extend='both',alpha=1.0)
    plt.colorbar(CS1,orientation='vertical',extend='both', shrink=0.5)
    
    plt.title('Var:{} - date:{}'.format(varname,currentdate))
    
    if filetype=="interpolated":
        plotfile=basepath+'{}_{}_a20_delta.png'.format(varname,currentdate)
    if filetype=="original":
        plotfile=basepath+'{}_{}_noresm_delta.png'.format(varname,currentdate)
    if filetype=="final":
        plotfile=basepath+'{}_{}_a20_delta+clim.png'.format(varname,currentdate)
    
    plt.savefig(plotfile,dpi=300)

def setupESMF(noresmfile,romsgridfile):
    esmgrid = ESMF.Grid(filename=noresmfile, filetype=ESMF.FileFormat.GRIDSPEC,
                                          is_sphere=True, coord_names=['lat','lon'],
                                          add_mask=False)

    romsgrid = ESMF.Grid(filename=romsgridfile, filetype=ESMF.FileFormat.GRIDSPEC,
                                            is_sphere=True, coord_names=['lat_rho','lon_rho'], add_mask=False)

    fieldSrc = ESMF.Field(esmgrid, "fieldSrc", staggerloc=ESMF.StaggerLoc.CENTER)
    fieldDst = ESMF.Field(romsgrid, "fieldDst",staggerloc=ESMF.StaggerLoc.CENTER)

    regridSrc2Dst = ESMF.Regrid(fieldSrc, fieldDst,regrid_method=ESMF.RegridMethod.BILINEAR,
                                                        unmapped_action=ESMF.UnmappedAction.IGNORE)
    return regridSrc2Dst, fieldSrc, fieldDst

def convertToCelsius(val):
    return val-273.15

print("Starting logfile for ESMF")
manager = ESMF.Manager(debug=True)
debug=True

filebase="/cluster/projects/nn9412k/A20/DELTA/results/"
finalbase="/cluster/projects/nn9412k/A20/DELTA/results/final/"
romsgridfile="/cluster/projects/nn9412k/A20/Grid/Arctic-20km_grd.nc"
figbasepath="/cluster/projects/nn9412k/A20/DELTA/Figures/"

if not os.path.exists(figbasepath): os.makedirs(figbasepath)
if not os.path.exists(finalbase): os.makedirs(finalbase)

cdf=Dataset(romsgridfile)
lat=cdf.variables["lat_rho"][:]
lon=cdf.variables["lon_rho"][:]
angle=cdf.variables["angle"][:]
cdf.close()
 
# Not that the order in myvardescription and myvars must be mapped
NORESMvardescriptions=["Fraction cloud cover (0-1)",
"Downward longwave radiation (W/m2)",
"Downward shortwave radiation (W/m2)",
"Atmospheric pressure (Pa)",
"Specific humidity at 2m (kg/kg)",
"Total precipitation (kg/m2/s)",
"Air temperature 2m (K)",
"Xi-component of wind (m/s)",
"Eta-component of wind (m/s)"]

NORESMvars=["CLDTOT","FLDS","FSDS","PS","QREFHT","PRECT","TREFHT","U","V"]
ROMSvars=["cloud","lwrad","swrad","Pair","Qair","rain","Tair","Uwind","Vwind"]
ROMSvardescriptions=["Fraction (0-1)","W/m2","W/m2","Pa","kg/kg","kg/m2s" "degC", "m/s","m/s"]
first=True
plev=25

for NORESMvar, ROMSvar, NORESMvardescription,ROMSvardescription in zip(NORESMvars,ROMSvars,NORESMvardescriptions,ROMSvardescriptions):
        
    print("Working on NORESM variable {} (ROMS equiv. {}) = {} (Roms equiv. {})".format(NORESMvar,ROMSvar,NORESMvardescription,ROMSvardescription))

    # Make a copy of teh ROMS climatology file to be used as the new output file 
    # after interpolation and adding of deltas
    if (ROMSvar=="swrad"):
        inf="swrad_daymean"
        src="{}A20_{}_1980_2014_detrend_monclim.nc".format(filebase,inf)
    else:
        src="{}A20_{}_1980_2014_detrend_monclim.nc".format(filebase,ROMSvar)
    dst="{}A20_{}_noresm_era_biascorrected_projections.nc".format(finalbase,ROMSvar)
    if os.path.exists(dst): os.remove(dst)
    copyfile(src, dst)

    # Open the climatology ERA file already on ROMS grid. The deltas from NORESM will 
    # be added to these climatological values
    cdf=Dataset(src)

    if (ROMSvar=="lwrad" or ROMSvar=="Pair" or ROMSvar=="Uwind" or ROMSvar=="Vwind"): 
        if (ROMSvar=="Pair"):
            time_variable="pair_time"
        if (ROMSvar=="lwrad"):
            ROMSvar="lwrad_down"
            time_variable="lwrad_time"
        if (ROMSvar=="Uwind" or ROMSvar=="Vwind"):
            time_variable="wind_time"
    else:
        time_variable="{}_time".format(ROMSvar)

    ROMSclimatology=cdf.variables[ROMSvar][:]
    time=cdf.variables[time_variable][:]
    cdf.close()

    # Opening output file to write results to
    cdfout=Dataset(dst,"a")
    noresmfile_deltas=filebase+"NORESM_ATM.cam2.hmlvl.2006-2100_delta.nc"
    noresmfile_mondeltas=filebase+"NORESM_ATM.cam2.hmlvl.2006-2100_2006-2015_deltas_monthly.nc"

    if first:
        regridSrc2Dst, fieldSrc, fieldDst = setupESMF(noresmfile_deltas,romsgridfile)
        first=False

    # Opening the original NORESM delta file to be interpolated to ROMS grid
    cdf=Dataset(noresmfile_deltas)
    cdf_mondelta=Dataset(noresmfile_mondeltas)
    montas=cdf_mondelta.variables[NORESMvar][:]
    
    tas=cdf.variables[NORESMvar][:]
    lat_noresm=cdf.variables["lat"]
    lon_noresm=cdf.variables["lon"]
    lon_noresms,lat_noresms=np.meshgrid(lon_noresm,lat_noresm)
    timesteps=cdf.variables["time"][:]
    mycalendar = cdf.variables["time"].calendar
    myunits = cdf.variables["time"].units
    print("== Opended input file: {} containing {} timesteps".format(noresmfile_deltas,len(timesteps)))

    # Loop over all time steps, interpolate to ROMS grid and write to file
    if debug:
        yearstoplot=[2010,2020,2030,2040,2050]

    for t in range(len(timesteps)):
        
        currentdate = netCDF4.num2date(timesteps[t], myunits, mycalendar)
        monindex=int(currentdate.month-1)

        # The NorESM wind has a fourth component plev which we set to surface level (plev=0)
        if (ROMSvar=="Uwind" or ROMSvar=="Vwind"):
            inputdata=np.flipud(np.rot90(np.squeeze(tas[t,plev,:,:])))
        else:
            inputdata=np.flipud(np.rot90(np.squeeze(tas[t,:,:])))
       
        if ROMSvar=="Tair":
            inputdata=convertToCelsius(inputdata)
        
        if ROMSvar=="rain":
            inputdata_mondelta=np.flipud(np.rot90(np.squeeze(montas[monindex,:,:])))
            fieldSrc.data[:,:]=inputdata_mondelta
            fieldDst = regridSrc2Dst(fieldSrc, fieldDst)
            outputdata_mondeltas=np.flipud(np.rot90(fieldDst.data))

        fieldSrc.data[:,:]=inputdata
        fieldDst = regridSrc2Dst(fieldSrc, fieldDst)
        
        monclim=np.squeeze(ROMSclimatology[monindex])
        outputdata=np.flipud(np.rot90(fieldDst.data))
        if (ROMSvar=="ua"):
            scru=(outputdata*np.cos(angle)) + (outputdata*np.sin(angle))
            outputdata=scru
        if (ROMSvar=="va"):
            scrv=(outputdata*np.cos(angle)) - (outputdata*np.sin(angle))
            outputdata=scrv

        if (ROMSvar=="rain"):
            # Calculate the change as fractional
            finaldata=(outputdata/outputdata_mondeltas)*monclim
            finaldata=np.where(finaldata<0,0,finaldata)
            
        else:
            finaldata=outputdata+monclim
            
        if (ROMSvar=="cloud"):
            finaldata=np.where(finaldata<0,0,finaldata)
            finaldata=np.where(finaldata>1,1,finaldata)
        
        #  Calculate current timestamp using days since 1948
        
        JD=netCDF4.date2num(currentdate,"days since 1948-01-01 00:00:00",mycalendar)
        if (t==0 or t==len(timesteps)): print("First timestep {}: {} (julian day {})".format(ROMSvar,currentdate,JD))
       
        cdfout.variables[ROMSvar][t,:,:]=finaldata
        cdfout.variables[time_variable][t]=JD

        print("Finished doing interpolation for {} for date {}".format(ROMSvar,currentdate))
        
        if debug:
            if currentdate.year in yearstoplot and ROMSvar=="rain":
                # Routines that plots: 1. the original noresm data on the noresm grid, 2. the interpolated data on the A20 grid,
                # and 3. the final result of summing delta plus climatology on the a20 grid
                if (ROMSvar=="Uwind" or ROMSvar=="Vwind"):
                    createContourPlot(np.squeeze(tas[t,plev,:,:]),lon_noresms,lat_noresms,figbasepath,ROMSvar,currentdate,"original")
                else:
                    createContourPlot(np.squeeze(tas[t,:,:]),lon_noresms,lat_noresms,figbasepath,ROMSvar,currentdate,"original")
                
                createContourPlot(outputdata,lon,lat,figbasepath,ROMSvar,currentdate,"interpolated")
                createContourPlot(finaldata,lon,lat,figbasepath,ROMSvar,currentdate,"final")

    cdfout.close()




