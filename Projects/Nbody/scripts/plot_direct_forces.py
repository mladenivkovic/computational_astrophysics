#!/usr/bin/python3



#=================================================
# Plots average (magnitude of) force of particles
# as a function of the particle distance from the
# origin, for various softening factors alpha,
# as well as the analytical profile.
# Doesn't need any cmdline args.
#=================================================


from os import getcwd, listdir
from sys import argv
import numpy as np
import matplotlib.pyplot as plt




#==========================
def read_commons(f_soft):
#==========================

    """ 

    reads in the data from given file.

    PARAMETERS:
        filename:   string of filename


    RETURNS:
        radius, mass, force
        numpy arrays of read in data: radius and force

    """


    print("Reading in file")


    #==================
    # read in header
    #==================

    #  infofile = "info_"+f_soft+'.txt'
    #
    #  f = open(infofile)
    #  header = f.readline() # skip first line
    #  values = f.readline()
    #  f.close()
    #
    #  values = values.split()
    #  scale_m = float(values[0])
    #  scale_l = float(values[1])
    #  scale_t = float(values[2])
    #  softening = float(values[3])




    #===============
    # Read in data
    #===============

    datafile = "output_direct_force_"+f_soft+".dat"
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

    nbins = 200
    maxbin = radius.max()
    minbin = radius.min()

    bins = np.zeros((nbins))

    for i in range(nbins):
        bins[i] = minbin * (maxbin/minbin)**(i/nbins)

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
def get_average_force(inds, ppb, f_soft):
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

   
    print("Computing average force for alpha =", s)


    datafile = "output_direct_force_"+f_soft+".dat"
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
    # get all softenings
    #-----------------------

    files=listdir()

    softenings = []
    for name in files:
        if name[:12] == "info_direct_":
            softening = name
            softening = softening.replace("info_direct_", "")
            softening = softening.replace(".txt","")
            softenings.append(softening)


    if (len(softenings)<1):
        print("No softenings found. Are you in the right directory?")
        quit()

    softenings.sort(key=float)
    print("found data for softenings:", softenings)



    #--------------------------------------
    # Find data valid for all softenings
    #--------------------------------------


    f_soft = softenings[0] #pick one

    radius, mass = read_commons(f_soft)

    inds, bins, cum_mass, ppb = bin_particles(radius, mass)

    a = find_a(bins, cum_mass)


    #-----------------------------
    # Prepare figure
    #-----------------------------

    fig = plt.figure(facecolor='white', figsize=(16, 9))
    ax1 = fig.add_subplot(111)




    #----------------------------------
    # Calculate and plot average force
    # for various softenings
    #----------------------------------


    for s in softenings:

        av_force = get_average_force(inds, ppb, s)

        print("Plotting alpha =", s)
        ax1.plot(bins[av_force>0], av_force[av_force>0],
                marker='o',
                linestyle='-',
                markersize=3,
                label=(r'$\alpha = $ ' + s))




    #----------------------------------
    # Calculate and plot average force
    # for analytical solution
    #----------------------------------
   
    analytical_solution = analytical_average_force(bins[:-1], bins[1:], a, cum_mass[-1], mass[0])

    print("Plotting analytical solution")
    ax1.plot(bins[1:], analytical_solution, 'k',
            label='analytical average')




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

    ax1.set_title(r'average force profile for varying softening', 
            family='serif', 
            size=18)
    
    ax1.grid()

    plt.tight_layout()

    
    plt.figtext(.02, .03, str(bins.shape[0]-1)+' bins used', family='serif', size=12)
    plt.figtext(.02, .01, r'$\epsilon = \alpha \cdot$ mean interparticle distance', family='serif', size=12)


    #-------------------
    # save figure
    #-------------------


    workdir= str(getcwd())
    outputfilename = 'average_force_plot'
    fig_path = workdir+'/'+outputfilename+'.png'
    print("saving average force profile plot as "+fig_path)
    plt.savefig(fig_path, format='png', facecolor=fig.get_facecolor(), transparent=False, dpi=300)#,bbox_inches='tight' )
    print("Done")

    plt.close()


