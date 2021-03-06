
#!/bin/bash
# Load required CDO module
# module load CDO/1.9.3-intel-2018a

# Trond Kristiansen, 15-18 January 2019, San Francisco

export inputdir="/Users/trondkr/Dropbox/NIVA/DownscaleA20/BRY/"
export outdir="/Users/trondkr/Dropbox/NIVA/DownscaleA20/RESULTS/"
export rundir="/Users/trondkr/Dropbox/NIVA/DownscaleA20/"
export pattern="a20_bry_SODA_BCG_20060115_to_20151215.nc"

# Calculate global monthly, detrended values of atmospheric forcing from NorESM
for filename in $inputdir$pattern; do
    echo "Working with file:" $filename;
    
    # Get the filename without path and setup temporary files in work dir
    basefilename=$(basename $filename);
    export extractedyears=$outdir$"extractedyears.nc"
    export noseasonality=$outdir$"noseasonality.nc"
    export tmpfile3=$outdir$"tmp3.nc"
    export tmpfile4=$outdir$"tmp4.nc"
    export xstdh=$outdir$"Xstdh.nc"

    export detrended=$outdir${basefilename%.nc}"_detrend.nc"
    export detrendedclim=$outdir${basefilename%.nc}"_detrend_monclim.nc"
    export monthlymean=$outdir${basefilename%.nc}"_monclim.nc"

    # Extract the years of interest
    cdo --verbose -b F64 selyear,2006/2015 $filename $extractedyears
    
    # Calculate the monthly climatology
    cdo ymonmean $extractedyears $monthlymean

    # Remove the seasonality
    cdo ymonsub $extractedyears $monthlymean $noseasonality
    
    # Remove the linear trend from the timeseries with removed seasonality
    cdo trend $noseasonality $tmpfile3 $tmpfile4
    cdo subtrend $extractedyears $tmpfile3 $tmpfile4 $detrended

    # Calculate the monthly mean from the detrended timeseries 2006-2015
    cdo ymonmean $detrended $detrendedclim
    cdo ymonstd $detrended $xstdh
    
    # Calculate the deltas by subtracting the detrended, 
    # climatology from the full timeseries
    export deltafile=$outdir${basefilename%.nc}"_delta.nc"
    cdo ymonsub $filename $detrendedclim $deltafile
    cdo ymonsub $filename $detrendedclim $deltafile
    echo "Delta results saved to file:" $deltafile
    echo "Climatology results saved to file:" $detrendedclim

    rm $tmpfile3
    rm $tmpfile4
    rm $extractedyears
    rm $noseasonality
    rm $detrended
    rm $monthlymean
done
