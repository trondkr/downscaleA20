# Bias-correct NorESM climate data: create [ROMS](https://www.myroms.org/) forcing for downscaled projections

<img alt="GitHub" src="https://img.shields.io/github/license/trondkr/downscaleA20.svg"><img alt="GitHub last commit" src="https://img.shields.io/github/last-commit/trondkr/downscaleA20.svg">
## Ocean and biogeochemistry
The scripts used for this toolbox uses the Climate Data Operator toolbox heavily. Its a fantastic set of functions that enables working with large climate data much easier.

<b>Notes:</b>
When working here with smaller values of e.g. nitrate and silicate weneeded to ensure we used full precision of the dataset which is neabled with the options:
```bash 
cdo -b F64 .....
```

### Combine two BRY files
The new version of [model2roms](https://github.com/trondkr/model2roms) now produces BCG + PHYSICS into the same BRY file. But prior to that we had two separate BRY files. This script combines those two files into one, based on the new output file used as template. This makes it easy to compare the historical hindcast with the projections and to caclulate climatology and deltas for downscaling: (combineSODAAndBGCFileIntoOne.py)

### Step 1: createDeltas-NorESM-ocean.sh
```bash 
./createDeltas-NorESM-ocean.sh
```
Script that uses the [CDO toolbox](https://code.mpimet.mpg.de/projects/cdo/) to calculate statistics and trends from the NorESM files. The script calculates the detrended climatology and removes this climatology from the timeseries. The result is a file containig the residuals and trends inherent in the timeseries. These residuals are later added to teh hindcast climatology to create a biascorrected timeseries.

Load the module using 
```bash 
module load CDO/1.9.3-intel-2018a
```

The starting point is a collection of global NorESM files for longer time-periods (2006-2050). As the script performs the following calculations are created:
1. Extract data for period 2006-2015
2. Calculate the monthly climatology and remove from the time-series to calculate non-seasonal trends
3. Remove the trend from the original 2006-2015 dataset
4. Calculate the climatology for the de-trended dataset (2006-2015)
5. Remove the monthly climatology from the entire dataset (2006-2050)
6. The result is the deltas from the climatology which will later be added to the higher resolution ERA climatology to create downscaled projections of climate change

### Step 2: createDeltas-SODA-ocean.sh
```bash 
./createDeltas-SODA-ocean.sh
```
Creates the climatology of the hindcast timeseries. This timeseries is assumed to be reaalistic and in phase with observations. Here, we use the SODA (Simple Ocean Data Assimilation) dataset to create the hindcast using the model2roms toolbox. This toolbox handles SODA3 as input files to create the required hindcast BRY file.

### Step 3: combineClimatologyAndDeltas.py
This script reads the output from step 1 and 2 and combined the two BRY files. The deltas/residuals from step 1 is added to the climatology created in step 2. A new file, the bias-corrected BRY file, is created.

## Atmosphere
### interpolateNORESM_using_ESMF.py
Use this script to interpolate **from** NorESM grid **to** local ROMS grid. This script uses the fast and efficient interpolation of ESMF and requires the module to be installed to run. Install using Anaconda with command: 
```bash
conda install -c conda-forge esmpy
```

If you dont have Anaconda install follow these suggestions:

```bash
Create a new environment called OpenDrift:

conda config --add channels conda-forge
conda create --name OpenDrift python=3.7 basemap matplotlib gdal libnetcdf netCDF4 numpy scipy seaborn xarray

conda activate OpenDrift
conda install -c conda-forge basemap-data-hires
conda install -c conda-forge esmpy
conda install -c conda-forge seawater
conda install -c conda-forge geos proj4
conda install -c conda-forge cartopy
```

#### Prior to running this script you need to do the following.

 1. Create one file containing all of the atmospheric variables from the NorESM model
 2. In our case, CLDTOT came in a separate file which then had to be merged with the other atmospheric variables:
 ```bash
 cdo merge CLDTOT.cam2.hmlvl.2006-2100.nc atm.cam2.hmlvl.2006-2100.nc NORESM_ATM.cam2.hmlvl.2006-2100.nc
 ```
 
 3. Create the detrended climatology from hindcast timeseries [ERA dataset](https://www.ecmwf.int/en/forecasts/datasets/reanalysis-datasets/era-interim):
 
 ```bash 
 createClimatologyERA.sh
 ```
 
 4. Create the detrended climatology and remove from teh timeseries to create residuals/deltas:
 ```bash 
 createDeltasNorESM-atm.sh
 ```
    
### Step 3: createClimatologyERA.py
This script creates de-trended climatology for 2006-2015 of the ROMS grid interpolated forcing files created using the global AN ERA-INTERIM files. The atmospheric forcing files for the hindcast has already been created for the grid using the toolbox for atmosphere found in my [romstools](https://github.com/trondkr/romstools/tree/master/create_atmos_ROMS). This script then reads each individual file for each atmospheric variable and calulates the deltas and climatology. The variables used are :
```bash
declare -a ERAvars=( cloud lwrad swrad_daymean Pair Qair rain Tair Uwind Vwind )
```

### Step 4: createDeltasNorESM-atm.sh
This script creates the detrended climatology (2006-2015) and removes it from the timeseries (2006-2100) to create residuals/deltas. These deltas will be added to the climatology created in step 1.  

Finally, we run the script: 
```bash
interpolateNORESM_using_ESMF.py
```
This interpolates the ERA climatology and the NorESM deltas to the ROMS grid (A20 in this case). The interpolated deltas are added to the interpolated climatology for each grid point for each timestep. The result is written to a new NetCDF4 file. This result file is the bias-corrected atmospheric forcing file required to run ROMS.
