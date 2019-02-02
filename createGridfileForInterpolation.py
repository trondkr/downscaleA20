from __future__ import print_function
import numpy as np
import string, os
from netCDF4 import date2num, num2date, Dataset
from datetime import datetime

# Trond Kristiansne, 21.01.2019
# Script written to create the rewquired grid info for using CDO to interpolate from one grid
# to another. In this case from the ESM grid to teh curvlinear grid of A20
#
# Run on fram using the Anaconda Python package:
# module load Anaconda2/5.0.1
# conda activate OpenDrift

#
# http://www.climate-cryosphere.org/wiki/index.php?title=Regridding_with_CDO
# TAPAS Project, NIVA 

print("Started grid generation script")
gridfile="/cluster/projects/nn9412k/A20/Grid/Arctic-20km_grd.nc"
outfilename="/cluster/projects/nn9412k/A20/Grid/Arctic-20km_grd_CDO.txt"

cdf=Dataset(gridfile)
print(cdf)
lat_rho=cdf.variables["lat_rho"][:]
lat_u=cdf.variables["lat_u"][:]
lat_v=cdf.variables["lat_v"][:]
lon_rho=cdf.variables["lon_rho"][:]
lon_u=cdf.variables["lon_u"][:]
lon_v=cdf.variables["lon_v"][:]
cdf.close()

# Create outputfile
if os.path.exists(outfilename): os.remove(outfilename)

outfile=open(outfilename,'a')
print("Output will be written to file: {}".format(outfilename))
print("Arctic 20 contains {} lat-lon points".format(len(lat_rho)*len(lon_rho)))
    
header="gridtype=curvlinear\ngridsize={}\nnvertex=4\n".format(len(lat_rho[:,0])*len(lon_rho[0,:]))
outfile.writelines(header)

#    ------- V(i,j+1)-------
#    |                     |
#  U(i,j)     RHO(i,j)   U(i+1,j)
#    |                     |
#    -------- V(i,j)--------

xvals="\nxvals ="; yvals="\nyvals ="; skip=10000; counter=0
xbounds="\nxbounds ="; yvals="\nybounds =";
sep=" "
counterBounds=0
for yi in range(len(lat_rho)):
     for xi in range(len(lon_rho)):

        xvals=xvals+sep+str(lon_rho[yi,xi])
        yvals=yvals+sep+str(lat_rho[yi,xi])

        boundsX=[xi+1,xi,xi,xi]
        boundsY=[yi,yi,yi,yi+1]

        xbounds=xbounds+sep+str(lon_v[yi,xi+1])+sep+str(lon_u[yi,xi])+sep+str(lon_v[yi,xi])+sep+str(lon_u[yi+1,xi])
        ybounds=ybounds+sep+str(lat_v[yi,xi+1])+sep+str(lat_u[yi,xi])+sep+str(lat_v[yi,xi])+sep+str(lat_u[yi+1,xi])
        counterBounds+=4

        counter+=1
        if counter==skip:
            print("Added {} points (bounds total {})".format(skip, counterBounds))
            skip=skip+10000
            
outfile.writelines(xvals)
outfile.writelines(xbounds)
outfile.writelines(yvals)
outfile.writelines(ybounds)

print("Finished writing grid file")
outfile.close()


