#!/bin/bash

#=============================================
#This script automatically runs the multipole
#nbody program with different parameters.
#=============================================



homedir=$PWD
nbdir='/home/mivkov/UZH/Computational_Astrophysics/Projects/Nbody/'
param="nbody-params"

export OMP_NUM_THREADS=3

for bucket in 01; do
# for bucket in 01 10 20; do
    bdir=bucket"$bucket"
    if [[ ! -d $bdir ]]; then
        mkdir -p "$bdir"
    fi


    for theta_max in 0.1 0.2 0.4 0.5 0.6 0.8 1; do

        thdir="$bdir"/"$theta_max"
        
        if [[ ! -d "$thdir" ]]; then
            mkdir -p "$thdir"
        fi
        
        cd "$thdir"

        for order in 0 2; do
            # Write param file
            echo "verbose = 1;" >> "$param"
            echo "direct_force = 0;" >> "$param"
            echo "multipole = 1;" >> "$param"
            echo "ncellpartmax = ""$bucket"";" >> "$param"
            echo "theta_max = ""$theta_max"";" >> "$param"
            echo "scale_cube = 0;" >> "$param"
            echo "multipole_order =""$order"";" >> "$param"
            echo "f_softening = 0.01;" >> "$param"
           
            # start run

            echo ""
            echo "============================================================"
            echo "Started run for order ""$order"" theta_max ""$theta_max"" bucket ""$bucket"
            echo "============================================================"
            echo ""

            logfile=log_"$order"-"$theta_max.log"
            unbuffer nbody "$param" "$nbdir"/files/data.ascii 2>&1 | tee "$logfile" 

            # unbuffer: form "expect" package. The pipe doesn't wait until enough data to print has accumulated,
            # but directly passes.

            # remove param file to make room for next run
            rm "$param"
        
        done 

        # create plots
        plot_multipole_forces.py
        # plot_deviations_multipole.py

        cd "$homedir"
    done 

done

