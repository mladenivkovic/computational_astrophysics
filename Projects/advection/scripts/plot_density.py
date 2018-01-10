#!/usr/bin/python3


#===================================
# Plots the density profile for the
# advection problem.
#===================================




import numpy as np
import matplotlib.pyplot as plt
from os import getcwd
from sys import argv











if __name__ == "__main__":


    #-----------
    # Setup
    #-----------

    if (len(argv) != 2):
        print("Don't recognize your arguments.")
        print("Usage: plot_density.py srcfile.dat")
        quit(2)

    cwd = getcwd()
    srcfile = argv[1]

    print("Reading data from", srcfile)
    data = np.loadtxt(srcfile)
    t = data[0]
    rho = data[1:]
    
    nx = rho.shape[0]

    x = 0.0*rho 
    for i in range(nx):
        x[i] = i*nx

    
    


    #------------
    # Plot
    #------------

    fig = plt.figure(figsize=(16,12))

    ax = fig.add_subplot(111)

    ax.plot(x, rho)
    ax.set_xlabel('x',
            labelpad=10, 
            family='serif', 
            size=16)
    ax.set_ylabel('density',
            labelpad=10, 
            family='serif', 
            size=16)
            

    plt.figtext(.02, .03, str("t ="+str(t)), family='serif', size=12)

    plt.tight_layout()

    #----------------------
    # Save figure
    #----------------------

    workdir= str(getcwd())
    outputfilename = "density_plot_"+str(t)
    fig_path = workdir+'/'+outputfilename+'.png'
    print("saving density plot as "+fig_path)
    plt.savefig(fig_path, format='png', facecolor=fig.get_facecolor(), transparent=False, dpi=300)#,bbox_inches='tight' )
    print("Done density profile plot")

    plt.close()







