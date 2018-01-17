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
param="advection2d-params"

export OMP_NUM_THREADS=3




#=============================
run_simulations()
#=============================
{
    #-----------------------------------------------------
    # Function to run the simulation.
    # Needs arguments:
    # run_simulations parameter-file nx methodnumber u v
    #-----------------------------------------------------

    param=$1
    nx=$2
    methodnr=$3
    u=$4
    v=$5



    for profile in 0 1 2; do
        # Write param file
        echo "verbose = 1;" >> "$param"
        echo "use_2d = 1;" >> "$param"
        echo "courant_factor = 0.5;" >> "$param"
        echo "nx = ""$nx"";" >> "$param"
        echo "ny = ""$nx"";" >> "$param"
        echo "t_end = 50;" >> "$param"
        echo "density_profile = ""$profile"";" >> "$param"
        echo "method = ""$methodnr"";" >> "$param"
        echo "u = ""$u"";" >> $"param"
        echo "v = ""$v"";" >> $"param"


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


for u in 0 1 0.707106781; do

    #-----------------------------------------
    # Create new directory for each velocity
    #-----------------------------------------

    v=-1

    if [[ $u == 0 ]]; then
        v=1;
        veldir="u=0"
    elif [[ $u == 1 ]]; then
        v=0;
        veldir="u=1"
    elif [[ $u == 0.707106781 ]]; then
        v=$u
        veldir="u=0.5sqrt(2)"
    fi


    if [[ ! -d "$veldir" ]]; then
        mkdir -p "$veldir"
    fi

    cd "$veldir"



    # for method in pwconst pwlin minmod VanLeer; do
    for method in VanLeer; do
       
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


        for nx in 100 200 500; do

            #-----------------------------------------
            # Create new directory for nx
            #-----------------------------------------

            nxdir="nx=""$nx"
            if [[ ! -d "$nxdir" ]]; then
                mkdir -p "$nxdir"
            fi

            cd "$nxdir"

            # run simulations
            run_simulations "$param" "$nx" "$methodnr" "$u" "$v"

            # create plots
            postprocess_advection2d.py

            cd ..
        done

        cd ..

    done

    cd ..

done

