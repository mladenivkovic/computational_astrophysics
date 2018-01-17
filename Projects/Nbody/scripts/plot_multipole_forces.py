#!/usr/bin/python3


#=====================================================
# Creates a plot with both direct forces calculation
# and multipole method force calculation as function
# of radius.
# Doesn't need any cmdline args.
#=====================================================



from os import getcwd, listdir
from sys import argv
import numpy as np
import matplotlib.pyplot as plt




#==========================
def read_commons(datafile):
#==========================

    """ 

    reads in the data from given file.

    PARAMETERS:
        datafile:   string of filename


    RETURNS:
        radius, mass, force
        numpy arrays of read in data: radius and force

    """


    print("Reading in file")



    #===============
    # Read in data
    #===============

    radius, mass = np.loadtxt(datafile, skiprows=1, usecols=([3, 4]), unpack=True)


    return radius, mass




#======================================
def bin_particles(radius, mass):
#======================================
    """
    bins the particles according to their distance from origin
    calculates average force per bin


    PARAMETERS:
    radius: radii of particles
    mass:  masses of particles

    RETURNS:
    inds:           numpy array of bin index of each particle
    bins:           bin distance
    cum_mass:       cumulative mass profile
    parts_per_bin:  particles per bin
    """


    #============
    # define bins
    #============

    maxbin = radius.max()
    minbin = radius.min()

    bins = np.zeros((nbins))

    for i in range(nbins):
        bins[i] = minbin * (maxbin/minbin)**(i/nbins)
        #  bins[i] = (maxbin - minbin)/nbins * i

    bins = np.concatenate((bins, np.array([maxbin])))


    #=================
    # bin particles
    #=================

    inds = np.digitize(radius, bins, right=True)
    
    parts_per_bin = np.bincount(inds)

    #  av_force = np.bincount(inds, weights = force)
    #  av_force = av_force/parts_per_bin


    cum_mass = np.bincount(inds, weights=mass)

    for i in range(len(cum_mass)-1):
        cum_mass[i+1] += cum_mass[i]


    return inds, bins, cum_mass, parts_per_bin









#==============================
def find_a(bins, cum_mass):
#==============================
    """
    Find the parameter a, given by M(a) = M/4

    PARAMETERS:
    bins:       bin distances
    cum_mass:   cumulative mass profile

    RETURNS:
    a:          parameter a
    """

    a = None
    totmass = cum_mass[-1]

    for i in range(len(cum_mass)):
        if totmass/4 <= cum_mass[i]:
            a = bins[i]
            print("Found a =", a)
            break

    return a








#=========================================
def get_average_force(inds, ppb, datafile):
#=========================================
    """
    Computes the average force for softening f_soft

    PARAMETERS:
    inds:       particle bin indices
    ppb:        particles per bin
    f_soft:     softening parameter for filename


    RETURNS:
    av_force:   average force per bin

    """

   
    print("Computing average force ")


    force = np.loadtxt(datafile, skiprows=1, usecols=([8]), unpack=True)


    av_force = np.bincount(inds, weights = force)
    av_force[ppb>0] = av_force[ppb>0]/ppb[ppb>0]


    return av_force




#====================================================
def analytical_average_force(r1, r2, a, totmass, pm):
#====================================================
    
    """
    compute the analytical average force in the interval
    r1 -> r2.

    force = - grad phi;
    => average force = 1/(r2 - r1) * integral (force) dr r1 -> r2 
                     = 1/(r2 - r1) * [ - M * G /(r+a) ]_r1 ^r2

    Parameter:
    r1, r2:      radii that determine the interval
    a:           parameter a for hernquist model
    totmass:     total mass of system
    pm:          mass of 1 particle


    RETURNS:
    av_force:    average force in given interval
    """




    av_force = 1.0/(r2 - r1) * (-totmass) *( 1/(r2+a) - 1/(r1 + a) )
    av_force = av_force * pm


    return av_force








#=========================
if __name__=="__main__":
#=========================


    #-----------------------
    # get all filenames
    #-----------------------
    

    files=listdir()

    orders = []

    for name in files:
        if name[:17] == "output_multipole_":
            order = name
            order = order.replace("output_multipole_", "")
            order = order.replace(".dat", "")
            order, minus, thetamax = order.partition("-")
            orders.append(order)


    if (len(orders)<1):
        print("No multipole orders found. Are you in the right directory?")
        quit()



    orders.sort(key=float)
    print("found data for multipole orders:", orders)






    #--------------------------------------
    # Find data valid for all runs
    #--------------------------------------


    order = orders[0] #pick one

    datafile = "output_multipole_"+order+"-"+str(thetamax)+".dat"
    radius, mass = read_commons(datafile)

    nbins = 200
    inds, bins, cum_mass, ppb = bin_particles(radius, mass)

    a = find_a(bins, cum_mass)





    #-----------------------------
    # Prepare figure
    #-----------------------------

    fig = plt.figure(facecolor='white', figsize=(16, 9))
    ax1 = fig.add_subplot(111)




    #----------------------------------
    # Calculate and plot average force
    # for analytical solution
    #----------------------------------
   
    analytical_solution = analytical_average_force(bins[:-1], bins[1:], a, cum_mass[-1], mass[0])

    print("Plotting analytical solution")
    ax1.plot(bins[1:], analytical_solution,
            label='analytical average')





    #----------------------------------
    # Calculate and plot average force
    # for various multipoles
    #----------------------------------



    for o in orders:

        datafile = "output_multipole_"+o+"-"+str(thetamax)+".dat"
        print("Averaging and plotting ", datafile)
        av_force = get_average_force(inds, ppb, datafile)

        print("Plotting order =", o)
        #  ax1.errorbar(bins, av_force,yerr=np.sqrt(av_force),
        ax1.plot(bins, av_force,
                marker='o',
                linestyle='--',
                markersize=3,
                label=(r'order = ' + o))







    #-------------------------
    # Tweak plot
    #-------------------------

    ax1.legend()
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlabel(r'$r$ $[r_{max}=1]$', 
            labelpad=10, 
            family='serif', 
            size=16)
    ax1.set_ylabel(r'average force [code units]', 
            labelpad=10, 
            family='serif', 
            size=16)

    thetamax = thetamax.strip()
    ax1.set_title(r'average force profile for multipole method with $\theta =$'+str(thetamax), 
            family='serif', 
            size=18)
    
    ax1.grid()


    plt.tight_layout()

    
    plt.figtext(.02, .03, str(bins.shape[0]-1)+' bins used', family='serif', size=12)


    #-------------------
    # save figure
    #-------------------


    workdir= str(getcwd())
    outputfilename = 'multipole_forces_plot-'+str(thetamax)
    fig_path = workdir+'/'+outputfilename+'.png'
    print("saving average force profile plot as "+fig_path)
    plt.savefig(fig_path, format='png', facecolor=fig.get_facecolor(), transparent=False, dpi=300)#,bbox_inches='tight' )
    print("Done")

    plt.close()


