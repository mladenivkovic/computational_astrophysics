#!/bin/bash

#=============================================
# This script automatically runs the advection
# program with different parameters.
#=============================================



homedir=$PWD
advdir='/home/mivkov/UZH/Computational_Astrophysics/Projects/advection/'
param="advection-params"

export OMP_NUM_THREADS=3


for method in pwconst; do
   
    #-----------------------------------------
    # Create new directory for each method
    #-----------------------------------------

    mdir="$method" 
    if [[ ! -d "$mdir" ]]; then
        mkdir -p "$mdir"
    fi

    cd "$mdir"

    if [[ $method == pwconst ]]; then
        methodnr=0
    fi


    for nx in 100 1000 10000; do

        #-----------------------------------------
        # Create new directory for nx
        #-----------------------------------------

        nxdir="nx=""$nx"
        if [[ ! -d "$nxdir" ]]; then
            mkdir -p "$nxdir"
        fi

        cd "$nxdir"

        for profile in 0 1 2; do
            # Write param file
            echo "verbose = 1;" >> "$param"
            echo "courant_factor = 0.9;" >> "$param"
            echo "nx = ""$nx"";" >> "$param"
            echo "t_end = 100;" >> "$param"
            echo "density_profile = ""$profile"";" >> "$param"
            echo "method = ""$methodnr"";" >> "$param"


            # start run
            echo ""
            echo "============================================================"
            echo Started run for order method $method nx=$nx profile=$profile
            echo "============================================================"
            echo ""

            logfile=log_"$method"_"$nx""_profile=""$profile"
            unbuffer advection "$param" 2>&1 | tee "$logfile"

            # unbuffer: form "expect" package. The pipe doesn't wait until enough data to print has accumulated,
            # but directly passes.

            # remove param file to make room for next run
            rm "$param"

        done;


        # create plots
        plot_density.py

        cd ..
    done

    cd ..

done

