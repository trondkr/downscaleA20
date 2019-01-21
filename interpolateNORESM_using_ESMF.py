
from __future__ import print_function
import numpy as np
import matplotlib
matplotlib.use('Agg')
import string, os, sys
from netCDF4 import date2num, num2date, Dataset
import netCDF4
from datetime import datetime
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from pylab import *
import datetime
import glob


try:
    import ESMF
except:
    print("Unable to continue without ESMF: install using : conda install -c conda-forge esmpy")
    sys.exit()

__author__   = 'Trond Kristiansen'
__email__    = 'me (at) trondkristiansen.com'
__created__  = datetime.datetime(2019, 1, 21)
__modified__ = datetime.datetime(2019, 1, 21)
__version__  = "1.0"
__status__   = "Development, 21.1.2019"


def createContourPlot(atmdata,lon,lat,basepath,currentdate):
    mymap = Basemap(projection='npstere', boundinglat=45,lon_0=0)
    x, y = mymap(lon,lat)
    mymap.drawcoastlines()
    mymap.fillcontinents(color='grey',zorder=2)
    mymap.drawcountries()
    mymap.drawmapboundary()
    levels=np.arange(atmdata.min(),atmdata.max(),2)

    CS1 = mymap.contourf(x,y,atmdata,levels,cmap=plt.cm.jet, extend='both',alpha=1.0)
    plt.colorbar(CS1,orientation='vertical',extend='both', shrink=0.5)
    
    #plt.title('Var:%s - depth:%s - time:%s'%(myvar,mydepth,currentDate))
   
    plotfile=basepath+'tas_{}.png'.format(currentdate)

    plt.savefig(plotfile,dpi=300)
    print("Created plot")

print("Starting logfile for ESMF")
manager = ESMF.Manager(debug=True)

noresmfile="/cluster/projects/nn9412k/A20/DELTA/results/tas_Amon_NorESM1-ME_rcp85_r1i1p1_200601-204412_detrend_monclim.nc"
romsgridfile="/cluster/projects/nn9412k/A20/Grid/Arctic-20km_grd.nc"
figbasepath="/cluster/projects/nn9412k/A20/DELTA/Figures/"

if not os.path.exists(figbasepath): os.makedirs(figbasepath)

cdf=Dataset(romsgridfile)
lat=cdf.variables["lat_rho"][:]
lon=cdf.variables["lon_rho"][:]
cdf.close()

esmgrid = ESMF.Grid(filename=noresmfile, filetype=ESMF.FileFormat.GRIDSPEC,
                                          is_sphere=True, coord_names=['lat','lon'],
                                          add_mask=False)

romsgrid = ESMF.Grid(filename=romsgridfile, filetype=ESMF.FileFormat.GRIDSPEC,
                                          is_sphere=True, coord_names=['lat_rho','lon_rho'], add_mask=False)

fieldSrc = ESMF.Field(esmgrid, "fieldSrc", staggerloc=ESMF.StaggerLoc.CENTER)
fieldDst = ESMF.Field(romsgrid, "fieldDst",staggerloc=ESMF.StaggerLoc.CENTER)

regridSrc2Dst = ESMF.Regrid(fieldSrc, fieldDst,regrid_method=ESMF.RegridMethod.BILINEAR,
                                                    unmapped_action=ESMF.UnmappedAction.IGNORE)

cdf=Dataset(noresmfile)
tas=cdf.variables["tas"][:]
time=cdf.variables["time"][:]

mycalendar = cdf.variables["time"].calendar
myunits = cdf.variables["time"].units
    
for t in range(len(time)):
    inputdata=np.flipud(np.rot90(np.squeeze(tas[t,:,:])))
    fieldSrc.data[:,:]=inputdata
    fieldDst = regridSrc2Dst(fieldSrc, fieldDst)
    currentdate = netCDF4.num2date(time[t], myunits, mycalendar)
    outputdata=np.flipud(np.rot90(fieldDst.data))

    print("Finished doing interpolation for {}".format(currentdate))
    createContourPlot(outputdata,lon,lat,figbasepath,currentdate)




