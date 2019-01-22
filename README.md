# downscaleA20

## createDeltasNorESM.sh
Script that uses the CDO toolbox to calculate statistics and trends from the NorESM files. The starting point is a collection of global NorESM files for longer time-periods (2006-2050). As the script performs the following calculations are created:
1. Extract data for period 2006-2015
2. Calculate the monthly climatology and remove from the time-series to calculate non-seasonal trends
3. Remove the trend from the original 2006-2015 dataset
4. Calculate the climatology for the de-trended dataset (2006-2015)
5. Remove the monthly climatology from the entire dataset (2006-2050)
6. The result is the deltas from the climatology which will later be added to the higher resolution ERA climatology to create downscaled projections of climate change

## interpolateNORESM_using_ESMF.py
Use this script to interpolate **from** NorESM grid **to** local ROMS grid. This script uses the fast and efficient interpolation of ESMF and requires the module to be installed to run. Install using Anaconda with command: 
```bash
conda install -c conda-forge esmpy
```
