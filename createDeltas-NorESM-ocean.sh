
#!/bin/bash
# Load required CDO module
module load CDO/1.9.5-intel-2018a

# Trond Kristiansen, 15-18 January 2019, San Francisco
#
# NOTE: Xstdh.nc is already calculated from historical data running createDeltas-SODA.sh

export inputdir="/Users/trondkr/Dropbox/NIVA/DownscaleA20/BRY/"
export outdir="/Users/trondkr/Dropbox/NIVA/DownscaleA20/RESULTS/"
export rundir="/Users/trondkr/Dropbox/NIVA/DownscaleA20/"
export pattern="a20_bry_NORESM_20060115_to_20491215.nc"

# Time-periods to scale up interannual variability
declare -a years1=( 2006 2015 2025 2035 )
declare -a years2=( 2015 2025 2035 2050 )

# Calculate global monthly, detrended values of atmospheric forcing from NorESM
for filename in $inputdir$pattern; do
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

    # Calculate the hindcast monthly std to use for scaling factor/inflation factor
    # cdo monstd $detrended Xstdh.nc

    for index in ${!years1[@]}; do 
        echo "Sub-extracting 10 year periods for removing trends" $years1[$index], $years2[$index];
        # Calculate the deltas by subtracting the detrended, 
        # climatology from the full timeseries
        export subextractedyears=$outdir$"subextractedyears.nc"
        export xstdp=$outdir$"Xstdp.nc"
        export xstdh=$outdir$"Xstdh.nc"
        export rfac=$outdir$"Rfac.nc"

        export deltafile=$outdir${basefilename%.nc}"_"${years1[$index]}"_"${years2[$index]}"_delta.nc"
        export deltafilecorrected=$outdir${basefilename%.nc}"_"${years1[$index]}"_"${years2[$index]}"_delta_stdcorr.nc"
        cdo --verbose -b F64 selyear,${years1[$index]}/${years2[$index]} $filename $subextractedyears

        # Calculate the projection monthly std to use for scaling factor/inflation factor for the variability
        cdo ymonstd $subextractedyears $xstdp
        cdo ymondiv $xstdh $xstdp $rfac

        # Remove the climatology from the timeseries
        #cdo ymonsub extracted.nc $detrendedclim $deltafile
        cdo ymonsub $subextractedyears $detrendedclim $deltafilecorrected

        # Inflate the variability using the ratio of std from hindcast to projections
        #cdo ymonmul -enlarge,$rfac $rfac $deltafile $deltafilecorrected

        echo "Delta results saved to file:" $deltafile
        echo "Climatology results saved to file:" $detrendedclim
        rm $subextractedyears
        rm $xstdh
        rm $xstdp 
        rm $rfac       
    done
        
    export deltafile=$outdir${basefilename%.nc}"_delta.nc"
    export deltafilecorrected=$outdir${basefilename%.nc}"*_delta_stdcorr.nc"
       
    cdo mergetime $deltafilecorrected $deltafile
    rm $tmpfile3
    rm $tmpfile4

    rm $deltafilecorrected
    rm $tmpfile3
    rm $subextractedyears
    rm $tmpfile4
    rm $extractedyears
    rm $noseasonality
    rm $detrended
    rm $monthlymean
done
