#!/usr/bin/python3


#=====================================================
# Compares time and accuracy of direct forces vs 
# multipole forces.
# Doesn't need any cmd line args, run in the
# Nbody/results directory.
#=====================================================



from os import getcwd, listdir
from sys import argv
import numpy as np
import matplotlib.pyplot as plt




#==========================
def read_commons(datafile):
#==========================

    """ 

    reads in the data common to all runs from given file.

    PARAMETERS:
        datafile:   string of filename


    RETURNS:
        radius, mass, force
        numpy arrays of read in data: radius and force

    """


    print("Reading in commons")



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



#=============================
def get_error(force_array):
#=============================
    """
    Compute the L1 error
    """

    global analytical_solution

    error = 0
    n = force_array.shape[0]
    for i in range(n-1):
        error += abs(force_array[i] - analytical_solution[i])
    error /= (n-1)
    return error



#==========================
def get_time(datafile):
#==========================
    """
    Get the time from datafile
    returns : time
    """

    time = np.loadtxt(datafile, skiprows=1, usecols=([4]))

    return time




#=========================
if __name__=="__main__":
#=========================


    #-----------------------
    # get all filenames
    #-----------------------
    
    src_mltp = 'multipole_forces/'
    src_direct = src_mltp+'direct_force/'

    dirforcefile = src_direct+'output_direct_force_0.01.dat' 
    dirinfofile = src_direct+'info_direct_0.01.txt' 

    src_mltp += 'bucket01/'

    monopole_files = []
    quadrupole_files = []
    mono_info_files = []
    quad_info_files = []
    theta = [0.1, 0.2, 0.4, 0.5, 0.6, 0.8, 1]

    for thetamax in theta:
        path=src_mltp+str(thetamax)+"/"
        monopole_files.append(path+"output_multipole_0-"+str(thetamax).rjust(3)+".dat")
        quadrupole_files.append(path+"output_multipole_2-"+str(thetamax).rjust(3)+".dat")
        mono_info_files.append(path+"info_multipole_0.txt")
        quad_info_files.append(path+"info_multipole_2.txt")
    





    radius, mass = read_commons(dirforcefile)

    nbins = 200
    inds, bins, cum_mass, ppb = bin_particles(radius, mass)

    a = find_a(bins, cum_mass)





    #-----------------------------
    # Prepare figure
    #-----------------------------

    fig = plt.figure(facecolor='white', figsize=(12, 6))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)




    #----------------------------------
    # Calculate and average force
    #----------------------------------

    analytical_solution = analytical_average_force(bins[:-1], bins[1:], a, cum_mass[-1], mass[0])
    direct_av_force = get_average_force(inds, ppb, dirforcefile)


    # get error on direct force
    error_direct = get_error(direct_av_force)
    time_direct = get_time(dirinfofile)

    error_monopole = []
    error_quadrupole = []

    time_monopole = []
    time_quadrupole = []




    for i in range(len(theta)):

        av_force_mono = get_average_force(inds, ppb,monopole_files[i])
        av_force_quad = get_average_force(inds, ppb,quadrupole_files[i])

        err_mono = get_error(av_force_mono)/error_direct
        error_monopole.append(err_mono)
        err_quad = get_error(av_force_quad)/error_direct
        error_quadrupole.append(err_quad)


        tm = get_time(mono_info_files[i])/time_direct
        time_monopole.append(tm)
        tq = get_time(quad_info_files[i])/time_direct
        time_quadrupole.append(tq)



    # plot errors
    ax1.plot(theta, error_monopole,
            marker='o',
            linestyle='-',
            markersize=3,
            label=('monopole'))
    ax1.plot(theta, error_quadrupole,
            marker='o',
            linestyle='-',
            markersize=3,
            label=('quadrupole'))

    ax2.plot(theta, time_monopole,
            marker='o',
            linestyle='-',
            markersize=3,
            label=('monopole'))
    ax2.plot(theta, time_quadrupole,
            marker='o',
            linestyle='-',
            markersize=3,
            label=('quadrupole'))





    #-------------------------
    # Tweak plot
    #-------------------------

    ax1.legend()
    #  ax1.set_xscale('log')
    #  ax1.set_yscale('log')
    ax1.set_xlabel(r'$\theta$',
            labelpad=10, 
            family='serif', 
            size=12)
    ax1.set_ylabel(r'Error[direct]/Error[multipole]',
            labelpad=10, 
            family='serif', 
            size=12)

    ax1.set_title(r'Error of multipole method in dependence of $\theta$', 
            family='serif', 
            size=14)
    
    ax1.grid()


    ax2.legend()
    #  ax1.set_xscale('log')
    #  ax1.set_yscale('log')
    ax2.set_xlabel(r'$\theta$',
            labelpad=10, 
            family='serif', 
            size=12)
    ax2.set_ylabel(r'time[multipole]/time[direct]',
            labelpad=10, 
            family='serif', 
            size=12)

    ax2.set_title(r'Time usage of multipole method in dependence of $\theta$', 
            family='serif', 
            size=14)
    
    ax2.grid()





    plt.tight_layout()

    
    plt.figtext(.02, .03, str(bins.shape[0]-1)+' bins used', family='serif', size=12)
    #  plt.figtext(.02, .01, r'$\epsilon = \alpha \cdot$ mean interparticle distance', family='serif', size=12)
    #

    #-------------------
    # save figure
    #-------------------


    workdir= str(getcwd())
    outputfilename = 'compare_direct_multipole'
    fig_path = workdir+'/'+outputfilename+'.png'
    print("saving comparison plot as "+fig_path)
    plt.savefig(fig_path, format='png', facecolor=fig.get_facecolor(), transparent=False, dpi=300)#,bbox_inches='tight' )
    print("Done")

    plt.close()


