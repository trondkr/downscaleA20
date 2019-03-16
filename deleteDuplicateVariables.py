
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
__status__   = "Development, 21.1.2019, 25.01.2019"


def removeDuplicateVariables():
    print("Running remove duplicates")


infile="RESULTS/a20_bry_NORESM_20060115_to_20501215_delta.nc"

cdf=Dataset(infile,'a')
print(cdf.variables)

atts = cdf.ncattrs()
del cdf.s_rho_2


cdf.close()

    