
import numpy as np
import xarray as xr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import spatial
import KDTreeIndex as kd
import pandas as pd
from pandas import Series, DataFrame, Panel
import matplotlib.cm as cm
from netCDF4 import Dataset, num2date, date2num
import os
import collections
    
def createPlot(timeseries, locnames): #, currentdates):

    for tsindex, ts in enumerate(timeseries):
       
        if not os.path.exists('figures'): os.mkdir('figures')
        ax=plt.gca()
        colormap = plt.cm.plasma #nipy_spectral, Set1,Paired  
        colors = [colormap(i) for i in np.linspace(0, 0.9,len(timeseries))]  
 
        saveddata_A = ts.resample("5A").median()
        datatset = saveddata_A.resample("10080T").median() # 10080 is 7 days
        tsint = datatset.interpolate(method='cubic')
        ax = tsint.plot(color=colors[tsindex],lw=2, ax=ax, label=locnames[tsindex])            
        if tsindex==len(timeseries)-1:

            leg=ax.legend(loc='upper left')
            for ind,te in enumerate(locnames):
                leg.get_texts()[tsindex].set_text(te)

            plotfilename='figures/NorESM_A20_stations.png'
            print("Saving to file: {}".format(plotfilename))
            if os.path.exists(plotfilename):
                os.remove(plotfilename)
            plt.savefig(plotfilename,bbox='tight', dpi=300)
        
     
#gridfile="/Users/trondkr/Dropbox/NIVA/A20/Grid/A20niva_grd_v1.nc"
datafile = 'a20_bry_NORESM_20060115_to_20501215_biascorrected.nc'
da = xr.open_dataset(datafile)
da=da.assign_coords(lon=da.lon_rho)
da=da.assign_coords(lat=da.lat_rho)
timevarname='ocean_time'
varname1='temp_east'

# To get the index for a set of locations/coordinates we use the kdtree.py algorithms. This 
# requires input data to be 2D so we just extract teh first timestep 
tos=da[varname1][0,:,:]

#da.assign_coords(time = da.indexes['time'].to_datetimeindex()) 
#kdtree = kd.KDTreeIndex(tos)

#locations=[]; locnames=[]; loc_indexes=[]; timeseries=[]
# Locations to plot
#DoggerBank = (55.0, 2.0); locations.append(DoggerBank); locnames.append("Dogger Bank")
#SouthernBight = (52.0, 3.0); locations.append(SouthernBight); locnames.append("Southern Bight")

#da.close()
cdf=Dataset(datafile)

#for location, locname in zip(locations,locnames):
 #   loc_index = kdtree.query(location)
 #   loc_indexes.append(loc_index)
  #  to_ts=np.squeeze(cdf.variables[varname1][:,0,loc_index[0],loc_index[1]])
    
to_ts=np.squeeze(cdf.variables[varname1][:,:,:])
print(to_ts)
timeseries=[] 
ddict = collections.OrderedDict()
ddict['timestamp'] = da.ocean_time.values   
ddict['var'] = to_ts
df = pd.DataFrame(ddict)
df = df.set_index('timestamp')
timeseries.append(df) 
print(df)
locnames=[]
locnames.append("test")
createPlot(timeseries,locnames)
    
    