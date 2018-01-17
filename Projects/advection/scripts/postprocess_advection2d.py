#!/usr/bin/python3


#================================================
# Plots the density profile for the advection 
# problem in 2d.
#================================================




import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size
from os import getcwd, listdir
from sys import argv




#==========================
def get_filedata():
#==========================
    """
    Gets all relevant files, sorts out the data

    returns:
        filedata:   list of lists of filedata
        counts:     list of counts for plots
    """

    #-----------------------
    # get all output files
    #-----------------------

    files=listdir()

    srcfiles = []
    profiles = []
    methods = []
    ncells = []
    times = []

    for name in files:
        if (name[:7] == "output_") and (name[-7:]=="-2d.dat"):
            srcfiles.append(name)

    if (len(srcfiles)<1):
        print("No output files found. Are you in the right directory?")
        quit()
    else:
        # sort sourcefiles alphabetically
        srcfiles.sort()
    
    for name in srcfiles:
        short = name.replace("-2d.dat", "")
        short = short.replace("output_", "")
        split = short.split("_")

        method = split[-2]
        methods.append(method)

        profile = ""
        for i in range(len(split)-2):
            profile += split[i]+" "
        profiles.append(profile)

        numbers = split[-1]
        nc, time = numbers.split("-")
        nc = int(nc)
        time = float(time)
        ncells.append(nc)
        times.append(time)





    #---------------------------------------------
    # Count various methods and profiles read in
    #---------------------------------------------

    countprofiles = []
    countmethods = []
    counttimes = []
    countcells = []

    for i in range(len(methods)):
        if methods[i] not in countmethods:
            countmethods.append(methods[i])

        if profiles[i] not in countprofiles:
            countprofiles.append(profiles[i])

        if times[i] not in counttimes:
            counttimes.append(times[i])

        if ncells[i] not in countcells:
            countcells.append(ncells[i])



    filedata = [srcfiles, profiles, methods, ncells, times]
    counts = [countprofiles, countmethods, counttimes, countcells]

    return filedata, counts



#======================
def getdata(srcfile):
#======================
    """
    Reads data from file <filename>
    returns:
        rho, nx, ny, x, y
        t:      time of output
        rho:    density profile
        nx:     number of cells
        x:      cell centres
    """

    print("Reading data from", srcfile)
    
    rho = np.loadtxt(srcfile, dtype='float')
    
    nx = rho.shape[1]
    ny = rho.shape[0]
   
    print("Found nx", nx, "ny", ny)


    return rho, nx, ny










#=====================================
def gettitle(profile, method, nx, ny):
#=====================================
    """
    Constructs the title for the plot.
    Returns:
    string title
    """
    global profiles

    title = profile+"profile"
    if method == "pwconst":
        title += " piecewise constant method"
    elif method == "pwlin":
        title += " piecewise linear method"
    elif method == "minmod":
        title += " with minmod slope limiter"
    elif method == "VanLeer":
        title += " with VanLeer slope limiter"

    if (nx == ny):
        title+=" nx=ny="+str(nx)
    else:
        title += " nx="+str(nx)+" , ny="+str(ny)

    return title







#======================================
def get_theoretical_value(x, profile):
#======================================
    """
    Compute the theoretical value 
    returns: rho_theory 
    """
    
    rho_theory = 1.0*x

    if (profile == "linear step "):
        rho_theory[x<=0.3] = 1
        rho_theory[x>0.6] = 1
        rho_theory[(0.3<x) & (x<=0.6)] = 2 + 1.7*(x[(0.3<x) & (x<=0.6)]-0.3)
    
    elif (profile == "step function "):
        rho_theory[x<=0.3] = 1
        rho_theory[(x>0.3) & (x<=0.6)] = 2
        rho_theory[x>0.6] = 1

    elif (profile == "gauss "):
        rho_theory = 1 + np.exp(-(x-0.5)**2/0.1)


    return rho_theory





#==============================
def get_colorscale(profile):
#==============================
    """
    Returns the minimal/maximal value for the
    colorscale 
    """

    if profile == "linear step ":
        vmin = 1.8
        vmax = 5
    elif profile == "step function ":
        vmin = 1.8
        vmax = 5
    elif profile == "gauss ":
        vmin = 1.8
        vmax = 5

    return vmin, vmax







#===============================
if __name__ == "__main__":
#===============================


    #-----------
    # Setup
    #-----------

    print("===========================")
    print("Started advection plotting")
    print("===========================")

    
    filedata, counts = get_filedata()

    srcfiles = filedata[0] 
    profiles = filedata[1] 
    methods = filedata[2] 
    ncells = filedata[3]
    times = filedata[4]

    countprofiles = counts[0]
    countmethods = counts[1]
    counttimes = counts[2]
    countcells = counts[3]




    #-------------------
    # Plotting Loop
    #-------------------

    
    for i in range(len(srcfiles)):

        fig = plt.figure(figsize=(8,8), dpi=300)
        ax = fig.add_subplot(111)

                
        rho, nx, ny = getdata(srcfiles[i])

        #  vmin, vmax = get_colorscale(profiles[i])

        im=ax.imshow(rho, 
                interpolation='none', 
                cmap='jet', 
                origin='lower', 
                extent=[1/nx,1,1/ny,1])
                #  vmin=vmin,
                #  vmax=vmax,

        # Make colorbar same height as plot
        divider = make_axes_locatable(ax)
        cax1 = divider.append_axes("right", size="2%", pad=0.1)
        fig.colorbar(im, cax=cax1)

        title = gettitle(profiles[i], methods[i], nx, ny)
        ax.set_title(title,
            family='serif',
            size=16,
            position=(0.5, 1.05))

        ax.set_xlabel('x',
                labelpad=10,
                family='serif',
                size=12)
        ax.set_ylabel('y',
                labelpad=10,
                family='serif',
                size=12)


        plt.figtext(.02, .03, str("t ="+str(times[i])), family='serif', size=12)

        plt.tight_layout()

        #----------------------
        # Save figure
        #----------------------

        workdir= getcwd()
        profile_str = profiles[i].replace(" ","_")
        outputfilename = "plot_advection_2d_"+profile_str+methods[i]+"_nx="+str(nx)+"_ny="+str(ny)+"t="+str(times[i])
        fig_path = workdir+'/'+outputfilename+'.png'
        print("saving advection density plot as "+fig_path)
        plt.savefig(fig_path, format='png', facecolor=fig.get_facecolor(), transparent=False, dpi=300)#,bbox_inches='tight' )
        print("Done density profile plot\n")

        #  plt.close()







