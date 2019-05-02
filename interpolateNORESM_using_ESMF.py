
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
import xarray as xr
import cftime
import pandas as pd
try:
    from itertools import ifilter
except ImportError:  # Python 3
    ifilter = filter

try:
    import ESMF
except:
    print("Unable to continue without ESMF: install using : conda install -c conda-forge esmpy")
    sys.exit()

__author__   = 'Trond Kristiansen'
__email__    = 'me (at) trondkristiansen.com'
__created__  = datetime.datetime(2019, 1, 21)
__modified__ = datetime.datetime(2019, 4, 16)
__version__  = "1.0"
__status__   = "Development, 21.1.2019, 25.01.2019, 16.03.2019, 16.04.2019"

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
#
# srun --nodes=1 --time=00:30:00 --account=nn9297k --partition=bigmem --ntasks-per-node=1 --mem-per-cpu=128G --qos=devel --pty bash -i
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
debug=False

filebase="/cluster/projects/nn9412k/A20/DELTA/results/"
finalbase="/cluster/projects/nn9412k/A20/DELTA/results/final/"
romsgridfile="/cluster/projects/nn9412k/A20/Grid/Arctic-20km_grd.nc"
figbasepath="/cluster/projects/nn9412k/A20/DELTA/Figures/"

if debug:
    filebase="/Users/trondkr/Dropbox/NIVA/DownscaleA20/test/results/"
    finalbase="/Users/trondkr/Dropbox/NIVA/DownscaleA20/test/results/final/"
    romsgridfile="/Users/trondkr/Dropbox/NIVA/DownscaleA20/test/Arctic-20km_grd.nc"
    figbasepath="/Users/trondkr/Dropbox/NIVA/DownscaleA20/test/Figures/"

if not os.path.exists(figbasepath): os.makedirs(figbasepath)
if not os.path.exists(finalbase): os.makedirs(finalbase)

cdf=Dataset(romsgridfile)
lat=cdf.variables["lat_rho"][:]
lon=cdf.variables["lon_rho"][:]
angle=cdf.variables["angle"][:]
cdf.close()
 
# Not that the order in myvardescription and myvars must be mapped
NORESMvardescriptions=["ppm","Fraction cloud cover (0-1)",
"Downward longwave radiation (W/m2)",
"Downward shortwave radiation (W/m2)",
"Atmospheric pressure (Pa)",
"Specific humidity at 2m (kg/kg)",
"Total precipitation (kg/m2/s)",
"Air temperature 2m (K)",
"Xi-component of wind (m/s)",
"Eta-component of wind (m/s)"]

NORESMvars=["CO2","CLDTOT","FLDS","FSDS","PS","QREFHT","PRECT","TREFHT","U","V"]
ROMSvars=["xCO2atm","cloud","lwrad","swrad","Pair","Qair","rain","Tair","Uwind","Vwind"]
ROMSvardescriptions=["ppm","Fraction (0-1)","W/m2","W/m2","Pa","kg/kg","kg/m2s" "degC", "m/s","m/s"]

NORESMvars=["FLDS","FSDS","PS","QREFHT","PRECT","TREFHT","U","V"]
ROMSvars=["lwrad","swrad","Pair","Qair","rain","Tair","Uwind","Vwind"]
ROMSvardescriptions=["W/m2","W/m2","Pa","kg/kg","kg/m2s" "degC", "m/s","m/s"]

NORESMvardescriptions=["Downward longwave radiation (W/m2)",
"Downward shortwave radiation (W/m2)",
"Atmospheric pressure (Pa)",
"Specific humidity at 2m (kg/kg)",
"Total precipitation (kg/m2/s)",
"Air temperature 2m (K)",
"Xi-component of wind (m/s)",
"Eta-component of wind (m/s)"]

if debug:
    NORESMvars=["CLDTOT"]
    ROMSvars=["cloud"]
    ROMSvardescriptions=["Fraction (0-1)"]
    NORESMvardescriptions=["Fraction cloud cover (0-1)"]

first=True
plev=25

for NORESMvar, ROMSvar, NORESMvardescription,ROMSvardescription in zip(NORESMvars,ROMSvars,NORESMvardescriptions,ROMSvardescriptions):
        
    print("Working on NORESM variable {} (ROMS equiv. {}) = {} (Roms equiv. {})".format(NORESMvar,ROMSvar,NORESMvardescription,ROMSvardescription))

    # Make a copy of teh ROMS climatology file to be used as the new output file 
    # after interpolation and adding of deltas
    if (ROMSvar=="swrad"):
        inf="swrad_daymean"
        src="{}A20_{}_1980_2014_detrend_monclim.nc".format(filebase,inf)
        srcdaily="{}A20_{}_2006_2060_delta.nc".format(filebase,inf)
    elif (ROMSvar=="xCO2atm"):
        src="{}A20_{}_1979_2015_detrend_monclim.nc".format(filebase,ROMSvar)
    else:
        src="{}A20_{}_1980_2014_detrend_monclim.nc".format(filebase,ROMSvar)
        srcdaily="{}A20_{}_2006_2060_delta.nc".format(filebase,ROMSvar)
    dst="{}A20_{}_noresm_era_biascorrected_projections.nc".format(finalbase,ROMSvar)
    if os.path.exists(dst): os.remove(dst)
    copyfile(srcdaily, dst)
    print("Finished copying output file {}".format(dst))

    if (ROMSvar=="lwrad" or ROMSvar=="Pair" or ROMSvar=="Uwind" or ROMSvar=="Vwind" or ROMSvar=="xCO2atm"): 
        if (ROMSvar=="Pair"):
            time_variable="pair_time"
        if (ROMSvar=="lwrad"):
            ROMSvar="lwrad_down"
            time_variable="lwrad_time"
        if (ROMSvar=="Uwind" or ROMSvar=="Vwind"):
            time_variable="wind_time"
        if (ROMSvar=="xCO2atm"):
            time_variable="xCO2atm_time"
            plev=0
    else:
        time_variable="{}_time".format(ROMSvar)

    # Open the climatology ERA file already on ROMS grid. The deltas from NORESM will 
    # be added to these climatological values
    cdf=xr.open_dataset(src, chunks={time_variable: 365})
    cdfdaily=xr.open_dataset(srcdaily, chunks={time_variable: 365*6})

    ROMSclimatology=cdf[ROMSvar]
    time=cdf[time_variable]
    finaldates=cdfdaily[time_variable]

    print("Finished loading data into memory")
    # Opening output file to write results to
    cdfout=Dataset(dst,'a')
    print("Opened output file {}".format(dst))
    noresmfile_deltas=filebase+"NORESM_ATM.cam2.hmlvl.2006-2100_delta.nc"
    noresmfile_mondeltas=filebase+"NORESM_ATM.cam2.hmlvl.2006-2100_2006-2015_deltas_monthly.nc"
   
    if (ROMSvar=="xCO2atm"):
        noresmfile_deltas=filebase+"CO2.cam2.hmlvl.2006-2100_delta.nc"
        noresmfile_mondeltas=filebase+"CO2.cam2.hmlvl.2006-2100_2006-2015_deltas_monthly.nc"
    if first:
        regridSrc2Dst, fieldSrc, fieldDst = setupESMF(noresmfile_deltas,romsgridfile)
        first=False

    # Opening the original NORESM delta file to be interpolated to ROMS grid
    cdf=xr.open_dataset(noresmfile_deltas, chunks={"time": 365})
    cdf_mondelta=xr.open_dataset(noresmfile_mondeltas, chunks={"time": 365})
    montas=cdf_mondelta[NORESMvar]
    if ROMSvar=="xCO2atm":
        # Convert from [kg/kg] to [ppm (dry air)] (NorESM to ROMS)
        fact=1e6*29./(12.+2.*16.)
        montas=montas*fact

    lat_noresm=cdf["lat"]
    lon_noresm=cdf["lon"]
    lon_noresms,lat_noresms=np.meshgrid(lon_noresm,lat_noresm)
    timesteps=cdf["time"]
    
    print("== Opened input file: {} containing {} timesteps ({} to {})".format(noresmfile_deltas,len(timesteps),timesteps[0].dt.year.values,timesteps[-1].dt.year.values))
 
    # Loop over all time steps, interpolate to ROMS grid and write to file
    if debug:
        yearstoplot=[2010,2020,2030,2040,2050]

    for t,timeobject in enumerate(timesteps):

        if timeobject.dt.year < 2060:
            monindex=int(timeobject.dt.month-1)
            print("=> Working on NorESM timestamp: {}.{}.{}:{}:{}".format(timeobject.dt.year.values,
            timeobject.dt.month.values,
            timeobject.dt.day.values,
            timeobject.dt.hour.values,
            timeobject.dt.second.values))
            # First we read the monthly climatology from ERA sourced ROMS forcing files
            # Second we read the NorESM fields for the same variable
            # Third the NorESM value is interpolated to the ROMS grid and added to the climatology
            # Fourth, the climatology + delta from NorESM is added to the frequenct (6 hours) variability file
            # to create a bias corrected projection that uses the daily variability of the hindcast as modulation 
            # of the monthly mean clim+delta.

            # The NorESM wind has a fourth component plev which we set to surface level (plev=0)
            if (ROMSvar=="Uwind" or ROMSvar=="Vwind" or ROMSvar=="xCO2atm"):
                inputdata=np.flipud(np.rot90(np.squeeze(cdf[NORESMvar][t,plev,:,:])))
            else:
                inputdata=np.flipud(np.rot90(np.squeeze(cdf[NORESMvar][t,:,:])))
        
            if ROMSvar=="Tair":
                inputdata=convertToCelsius(inputdata)
            
            # Variables that we use the fractional approach on - need monthly hindcast residuals
            if ROMSvar=="rain":
                inputdata_mondelta=np.flipud(np.rot90(np.squeeze(montas[monindex,:,:])))
                fieldSrc.data[:,:]=inputdata_mondelta
                fieldDstRes = regridSrc2Dst(fieldSrc, fieldDst)
                outputdata_mondeltas=np.flipud(np.rot90(fieldDstRes.data))

            if ROMSvar=="xCO2atm":
                inputdata_mondelta=np.flipud(np.rot90(np.squeeze(montas[monindex,plev,:,:])))
            
                fieldSrc.data[:,:]=inputdata_mondelta
                fieldDstRes = regridSrc2Dst(fieldSrc, fieldDst)
                outputdata_mondeltas=np.flipud(np.rot90(fieldDstRes.data))
        
            if ROMSvar=="xCO2atm":
                # Convert from [kg/kg] to [ppm (dry air)] (NorESM to ROMS)
                fact=1e6*29./(12.+2.*16.)
                inputdata=inputdata*fact
            fieldSrc.data[:,:]=inputdata
        
            fieldDstRes = regridSrc2Dst(fieldSrc, fieldDst)
        
            monclim=np.squeeze(ROMSclimatology[monindex])
            outputdata=np.flipud(np.rot90(fieldDstRes.data))
        
            if (ROMSvar=="ua"):
                scru=(outputdata*np.cos(angle)) + (outputdata*np.sin(angle))
                outputdata=scru
            if (ROMSvar=="va"):
                scrv=(outputdata*np.cos(angle)) - (outputdata*np.sin(angle))
                outputdata=scrv

            if (ROMSvar=="rain" or ROMSvar=="xCO2atm"):
                finaldata=outputdata+monclim
                finaldata=np.where(finaldata<0,0,finaldata)
                #print("final co2", np.min(finaldata),np.max(finaldata))
            else:
                finaldata=outputdata+monclim
                
            if (ROMSvar=="cloud"):
                finaldata=np.where(finaldata<0,0,finaldata)
                finaldata=np.where(finaldata>1,1,finaldata)

            if (ROMSvar=="swrad" or ROMSvar=="lwrad"):
                finaldata=np.where(finaldata<0,0,finaldata)
               
            print("Finished with interpolation of NorESM field")
            # Now that we have monthly climatology and delta we need to loop over the more frequent variability file and
            # add the values
            month=timeobject.dt.month
            year=timeobject.dt.year
            dateobjects=np.asarray(list(ifilter(lambda d: (d[1].dt.month==month and d[1].dt.year==year), enumerate(finaldates))))
            dateindexes=np.asarray([x[0] for x in dateobjects])
            if len(dateindexes) > 0:

            #    if debug:
                for dd in dateobjects:
                    if dd == dateobjects[0] or dd==dateobjects[-1]:
                        print("   => Working on ERA timestamp: {}.{}.{}:{}:{}".format(dd[1].dt.year.values,
                        dd[1].dt.month.values,
                        dd[1].dt.day.values,
                        dd[1].dt.hour.values,
                        dd[1].dt.second.values))
                print("Extracted all indexes - now updating arrays at {} indexes".format(len(dateindexes)))
                variability=cdfdaily[ROMSvar][dateindexes,:,:].load()
                updatedValues=variability+finaldata
                cdfout[ROMSvar][dateindexes,:,:]=updatedValues
                print("Updated values in netcdf file")
            # JDs=[]
            # cdftime_jul = cftime.utime("days since 1948-01-01 00:00:00", calendar='julian')
            # for do in dateobjects:
            #     JDs.append(cdftime_jul.date2num(do[1].dt))
                    # if debug:
                    #    print("Editing JD {} to {}".format(do[1].dt,cdftime_jul.date2num(do[1].dt)))
                #cdfout[time_variable][dateindexes]=np.asarray(JDs)
                
            #    print("Writing value to file {} for date {}.{} to {}.{}".format(ROMSvar, dateobjects[0][1].dt.year,dateobjects[0][1].dt.month,dateobjects[-1][1].dt.year,dateobjects[-1][1].dt.month))
                
            if debug:
                if timeobject.dt.year in yearstoplot and ROMSvar=="xCO2am":
                    # Routines that plots: 1. the original noresm data on the noresm grid, 2. the interpolated data on the A20 grid,
                    # and 3. the final result of summing delta plus climatology on the a20 grid
                    if (ROMSvar=="Uwind" or ROMSvar=="Vwind"):
                        createContourPlot(np.squeeze(tas[t,plev,:,:]),lon_noresms,lat_noresms,figbasepath,ROMSvar,timeobject.dt,"original")
                    else:
                        createContourPlot(np.squeeze(tas[t,:,:]),lon_noresms,lat_noresms,figbasepath,ROMSvar,timeobject.dt,"original")
                    
                    createContourPlot(outputdata,lon,lat,figbasepath,ROMSvar,timeobject.dt,"interpolated")
                    createContourPlot(finaldata,lon,lat,figbasepath,ROMSvar,timeobject.dt,"final")
    cdfout.close()

   





