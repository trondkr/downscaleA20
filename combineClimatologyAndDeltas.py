from __future__ import print_function
import numpy as np
import string, os
from netCDF4 import date2num, num2date, Dataset
from datetime import datetime
from shutil import copyfile
from cftime import DatetimeNoLeap
import pandas as pd

# Trond Kristiansne, 21.01.2019

climfile="/Users/trondkr/Dropbox/NIVA/DownscaleA20/RESULTS/a20_bry_SODA_BCG_20060115_to_20151215_detrend_monclim.nc"
deltafile="/Users/trondkr/Dropbox/NIVA/DownscaleA20/RESULTS/a20_bry_NORESM_20060115_to_20491215_delta.nc"
resultfile="/Users/trondkr/Dropbox/NIVA/DownscaleA20/RESULTS/a20_bry_NORESM_20060115_to_20491215_biascorrected.nc"

if os.path.exists(resultfile): 
    os.remove(resultfile)
copyfile(deltafile, resultfile)

directions=["east","west","north","south"]
romsvariables=['temp','salt','zeta','u','v','ubar','vbar','O3_c','O3_TA','N1_p','N3_n','N5_s','O2_o']
cdf=Dataset(resultfile,"a")
cdfclim=Dataset(climfile)

times=cdf.variables["ocean_time"][:]

dates=num2date(times, units=cdf.variables["ocean_time"].units, calendar=cdf.variables["ocean_time"].calendar)

for romsvar in romsvariables:
    for direction in directions:

        currentvar = "{}_{}".format(romsvar,direction)
        print("Adding delta for {}".format(currentvar))
        climdata=cdfclim.variables[currentvar][:]

        for month in range(1,13,1):
            
            for i,d in enumerate(dates):

                if d.month==month:
                    if romsvar in ["vbar","ubar","zeta"]:
                        cdf.variables[currentvar][i,:]=np.squeeze(cdf.variables[currentvar][i,:])+np.squeeze(climdata[month-1,:])
                    else:
                        corrected=np.squeeze(cdf.variables[currentvar][i,:,:])+np.squeeze(climdata[month-1,:,:])

                        if romsvar in ['O3_c','O3_TA','N1_p','N3_n','N5_s','O2_o']:
                            corrected_clamped=np.where(corrected<0,0,corrected)
                        else:
                            corrected_clamped=corrected

                        cdf.variables[currentvar][i,:,:]=corrected_clamped
cdfclim.close()
            
