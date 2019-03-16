from __future__ import print_function
import numpy as np
import string, os
from netCDF4 import date2num, num2date, Dataset
from datetime import datetime
from shutil import copyfile
from cftime import DatetimeNoLeap
import pandas as pd

# Trond Kristiansen, 14.03.2019
#
# The new version of model2roms now produces BCG + physics into the same BRY file. But prior to that we had
# data into two separate files. This script combines those two files into one, based on the new output file used as template.
# This makes it easy to compare the historical hindcast with the projections and to caclulate climatology and 
# deltas for downscaling.

# Prior to running this script create files containing only the time period 2006-2015

 # cdo --verbose selyear,2006/2015 a20_bry_SODA3_19800115_to_20151215.nc a20_bry_SODA3_20060115_to_20151215.nc
 # cdo --verbose selyear,2006/2015 a20_bry_NORESMOCv1p2bc_19800115_to_20161215.nc a20_bry_NORESMOCv1p2bc_20060115_to_20151215.nc

 # cdo --verbose selyear,2006/2015 a20_bry_NORESM_20060115_to_20491215.nc a20_bry_template_20060115_to_20151215.nc

bcgfile="BRY/a20_bry_NORESMOCv1p2bc_20060115_to_20151215.nc"
sodafile="BRY/a20_bry_SODA3_20060115_to_20151215.nc"
templatefile="BRY/a20_bry_template_20060115_to_20151215.nc"
resultfile="a20_bry_SODA_BCG_20060115_to_20151215.nc"

if os.path.exists(resultfile): 
    os.remove(resultfile)
copyfile(templatefile, resultfile)

cdfresult=Dataset(resultfile,"a")
cdfsoda=Dataset(sodafile)
cdfbcg=Dataset(bcgfile)

romsvariables = ['temp','salt','zeta','u','v','ubar','vbar','O3_c','O3_TA','N1_p','N3_n','N5_s','O2_o']
directions=["east","west","north","south"]

times=cdfbcg.variables["ocean_time"][:]

dates=num2date(times, units=cdfbcg.variables["ocean_time"].units, calendar=cdfbcg.variables["ocean_time"].calendar)

for romsvar in romsvariables:
    for direction in directions:

        currentvar = "{}_{}".format(romsvar,direction)
        print("Adding data for : {}".format(currentvar))

        if romsvar in ['O3_c','O3_TA','N1_p','N3_n','N5_s','O2_o']:
            # No south direction in input data - using NORESM values
            if direction in ["east","west","north"]:
                cdfresult.variables[currentvar][:]=cdfbcg.variables[currentvar][:]
        else:
            cdfresult.variables[currentvar][:]=cdfsoda.variables[currentvar][:]
                   
cdfresult.close()
cdfsoda.close()
cdfbcg.close()
            
