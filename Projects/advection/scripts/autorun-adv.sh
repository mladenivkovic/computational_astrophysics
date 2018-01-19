#!/bin/bash

#=============================================
# This script automatically runs the advection
# program with different parameters.
#=============================================



#------------------------
# SETUP
#------------------------



homedir=$PWD
advdir='/home/mivkov/UZH/Computational_Astrophysics/Projects/advection/'
param="advection-params"

export OMP_NUM_THREADS=3




#=============================
run_simulations(){
#=============================

    # Function to run the simulation.
    # Needs arguments:
    # run_simulations parameter-file nx methodnumber

    param=$1
    nx=$2
    methodnr=$3


    for profile in 0 1 2; do
        # Write param file
        echo "verbose = 1;" >> "$param"
        echo "courant_factor = 0.5;" >> "$param"
        echo "nx = ""$nx"";" >> "$param"
        echo "t_end = 100;" >> "$param"
        echo "density_profile = ""$profile"";" >> "$param"
        echo "method = ""$methodnr"";" >> "$param"


        # start run
        echo ""
        echo "============================================================"
        echo Started run for method $method nx=$nx profile=$profile
        echo "============================================================"
        echo ""

        logfile=log_"$method"_"$nx""_profile=""$profile"
        unbuffer advection "$param" 2>&1 | tee "$logfile"

        # unbuffer: form "expect" package. The pipe doesn't wait until enough data to print has accumulated,
        # but directly passes.

        # remove param file to make room for next run
        rm "$param"

    done;
}







#----------------------
# SIMULATION LOOP
#----------------------



# for method in pwconst pwlin ; do
for method in pwconst pwlin minmod VanLeer; do
   
    #-----------------------------------------
    # Create new directory for each method
    #-----------------------------------------

    mdir="$method" 
    if [[ ! -d "$mdir" ]]; then
        mkdir -p "$mdir"
    fi

    cd "$mdir"

    if [[ "$method" == "pwconst" ]]; then
        methodnr=0
    elif [[ "$method" == "pwlin" ]]; then
        methodnr=1
    elif [[ "$method" == "minmod" ]]; then
        methodnr=2
    elif [[ "$method" == "VanLeer" ]]; then
        methodnr=3
    fi


    # for nx in 100 1000 10000 500 2000 5000; do
    for nx in 200 ; do

        #-----------------------------------------
        # Create new directory for nx
        #-----------------------------------------

        nxdir="nx=""$nx"
        if [[ ! -d "$nxdir" ]]; then
            mkdir -p "$nxdir"
        fi

        cd "$nxdir"

        # run simulations
        run_simulations "$param" "$nx" "$methodnr"

        # create plots
        postprocess_advection1d.py

        cd ..
    done

    cd ..

done

