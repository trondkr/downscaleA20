
#!/bin/bash
# Load required CDO module
module load CDO/1.9.3-intel-2018a

# Trond Kristiansen, 17th March 2019, San Francisco

export inputdir="/cluster/projects/nn9412k/A20/FORCING/HINDCAST/"
export outdir="/cluster/projects/nn9412k/A20/DELTA/results/"
export rundir="/cluster/projects/nn9412k/A20/DELTA/"

declare -a ERAvars=( xCO2atm )

for index in ${!ERAvars[@]}; do
    export filename=${inputdir}"a20_xCO2atm_1979_2015.nc"
    echo "Working with file:" $filename;
    
    # Get the filename without path and setup temporary files in work dir
    basefilename=$(basename $filename);
    export extractedyears=$outdir$"extractedyears.nc"
    export noseasonality=$outdir$"noseasonality.nc"
    export tmpfile3=$outdir$"tmp3.nc"
    export tmpfile4=$outdir$"tmp4.nc"
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
