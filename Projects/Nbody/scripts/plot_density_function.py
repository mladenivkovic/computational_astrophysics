#!/usr/bin/python3


#=====================================================
# Task 1 Step 0:
# Plots the density function rho(r)
# usage: plot_density_function.py data.ascii
# writes results in current workdir
#=====================================================




import numpy as np
from os import getcwd
from sys import argv







#==========================
def readdata(filename):
#==========================

    """

    reads in the data from given file.

    PARAMETERS:
        filename:   string of filename


    RETURNS:
        mass, x, y, z
        numpy arrays of read in data

    """
 

    print("Reading in file")



    #==================
    # read in header
    #==================

    f = open(filename)
    header = f.readline()
    f.close()

    header = header.split()
    npart = int(header[0])

    
    #==================
    # read in data
    #==================
    
    rawdata = np.loadtxt(filename, dtype='float', skiprows=1)
    
    mass = rawdata[0:npart]
    x = rawdata[npart:2*npart]
    y = rawdata[2*npart: 3*npart]
    z = rawdata[3*npart: 4*npart]
    # vx = rawdata[4*npart: 5*npart]
    # vy = rawdata[5*npart: 6*npart]
    # vz = rawdata[6*npart: 7*npart]
    # softening = rawdata[7*npart: 8*npart]
    # potential = rawdata[9*npart: 10*npart]

   



    return mass, x, y, z #, vx, vy, vz, softening, potential




#=====================================
def plot3d(x, y, z, outputfilename):
#=====================================
    """ 
    Creates a 3D scatterplot and saves it as 'outputfilename.png'

    PARAMETERS:
    x, y, z :       x, y, z coordinates to plot
    outputfilename: filename to save to
    
    
    RETURNS:
    nothing
    
    """


    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D


    print("Creating 3d plot")


    #==================
    # Create figure
    #==================

    fig = plt.figure(facecolor='white', figsize=(10,11))
    ax = fig.add_subplot(111, projection='3d')


    pointsize = 2
    mylabel = 'None'
    mylw = 0
    mymarker = 'o'
    pointalpha = 0.5

    ax.scatter(x, y, z, 
            s=pointsize, 
            label=mylabel, 
            lw=mylw, 
            marker=mymarker, 
            depthshade=True, 
            alpha=pointalpha)






    # set tick params (especially digit size)
    ax.tick_params(axis='both', which='major', labelsize=12,top=5)

    #label axes
    ax.set_xlabel(r'x')#, labelpad=20, family='serif',size=16)
    ax.set_ylabel(r'y')#, labelpad=20, family='serif',size=16)
    ax.set_zlabel(r'z')#, labelpad=20, family='serif',size=16)


    plt.tight_layout()




    #================
    # Save figure
    #================

    workdir= str(getcwd())
    fig_path = workdir+'/'+outputfilename+'.png'
    print("saving 3d plot as "+fig_path)
    plt.savefig(fig_path, format='png', facecolor=fig.get_facecolor(), transparent=False, dpi=100)#,bbox_inches='tight' )
    print("Done 3d plot")

    plt.close()

    return
   






#======================================
def get_profiles(x, y, z, mass):
#======================================
    """

    calculates the density profile.

    PARAMETERS:
    y, x, z, mass:  coordinates and mass of particles


    RETURNS:
    bins, density_profile, mass_profile    
    bins:           bin distance (particles sorted in bins by distance from origin)
    density_profile,
    mass_profile :  np.array of density/mass profile
    """

    print("Calculating density profile")


    # calculate radius from origin
    r = np.sqrt(x**2 + y**2 + z**2)

    # norm r and m
    r = r / r.max()
    mass = mass/sum(mass)*1e6


    nbins = 200
    maxbin = r.max()
    minbin = r.min()


    # give bin widths.
    bins = np.zeros((nbins))

    for i in range(nbins):
        # bins[i] = minbin + i * (maxbin-minbin)/nbins # for linear bins
        bins[i] = minbin*(maxbin/minbin)**(i/nbins)

    bins = np.concatenate( (bins, np.array([maxbin])) )


    # get mass profile
    inds = np.digitize(r, bins, right=True)

    mass_profile = np.bincount(inds, weights=mass)
    # mass_profile = mass_profile / np.sum(mass_profile) # norm mass profile to 1

    # get density profile
    density_profile = np.zeros((mass_profile.shape[0]))
    density_profile[1:] = mass_profile[1:] / (4 / 3 * np.pi * (bins[1:]**3-bins[:-1]**3))
    density_profile[0] = mass_profile[0]/ (4 / 3 * np.pi * bins[0]**3 )

    return bins, density_profile, mass_profile





#=========================================================================
def plot_profiles(bins, density_profile, mass_profile, outputfilename):
#=========================================================================
    """
    plots the density profile. 

    PARAMETERS:

        r:               distances from origin of each particle
        bins:            bin positions
        density_profile: density profile to plot
        mass_profile:    mass profile
        outputfilename:  output filename

    RETURNS:
        
        nothing
    """


    import matplotlib.pyplot as plt

    print("Creating density plot")

    #==========================
    # prepare analytical model
    #==========================

    #----------------------------------
    # find a: cum_mass(a) = tot_mass/4
    #----------------------------------
    cum_mass = mass_profile*1.0

    for i in range(len(mass_profile)-1):
        cum_mass[i+1] += cum_mass[i]

    a = None
    totmass = cum_mass[-1]


    for i in range(len(mass)):
        if totmass/4 <= cum_mass[i]:
            a = bins[i]
            print("Found a =", a)
            break


    #------------------------------
    # Define analytical functions
    #------------------------------

    def rho(r, a, totmass):
        return totmass/(2 * np.pi) * a / r * 1/(r + a)**3


    def cmp(r, a, totmass):
        return totmass * r**2/(r+a)**2


    mass_profile_theory = np.zeros((bins.shape[0]))
    mass_profile_theory[1:]=cmp(bins[1:], a, totmass )-cmp(bins[:-1], a, totmass)



    #==================
    # Plotting
    #==================

    fig = plt.figure(facecolor='white', figsize=(12, 6))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    #  ax3 = fig.add_subplot(133)

    #---------------
    # get errorbars:
    #---------------
    # in some cases, where density_profile < 1, the error sqrt(dens_prof) > dens_prof
    # => error bar would get negative, which is not permitted for logarithmic axes
    # => replace error bar length there with 0.99999*dens_prof
    derr_lower = np.sqrt(density_profile)
    derr_upper = derr_lower*1.0
    derr_lower[density_profile<=derr_lower]=density_profile[density_profile<=derr_lower]*0.999999

    cum_err_lower = np.sqrt(cum_mass)
    cum_err_upper = cum_err_lower*1.0
    cum_err_lower[cum_mass<=cum_err_lower]=cum_mass[cum_mass<=cum_err_lower]*0.999999

    mass_err_lower = np.sqrt(mass_profile)
    mass_err_upper = mass_err_lower*1.0
    mass_err_lower[mass_profile<=mass_err_lower]=mass_profile[mass_profile<=mass_err_lower]*0.999999




    #----------------------
    # Plot density profile
    #----------------------

    ax1.errorbar(bins, density_profile, 
            yerr=(derr_lower, derr_upper),
            fmt='o',
            markersize=1, 
            capsize=2,
            elinewidth=1,
            label='binned data')

    ax1.plot(bins, rho(bins, a, totmass), label='theoretical model')



    #-----------------------
    # Plot mass profile
    #-----------------------
    ax2.errorbar(bins, mass_profile, 
            yerr=(mass_err_lower, mass_err_upper),
            fmt='o',
            markersize=1, 
            capsize=2,
            elinewidth=1,
            label='binned data')
    ax2.plot(bins, mass_profile_theory, label='theoretical model')


    #  #--------------------------------
    #  # Plot cumulative mass profile
    #  #--------------------------------
    #  ax3.errorbar(bins,
    #          cum_mass,
    #          yerr=(cum_err_lower, cum_err_upper),
    #          fmt='o',
    #          markersize=1,
    #          capsize=2,
    #          elinewidth=1,
    #          label='binned data')
    #
    #  ax3.plot(bins, cmp(bins, a, totmass), label='theoretical model')
    #



    #----------------------
    # Tweak plots
    #----------------------

    # set tick params (especially digit size)
    ax1.tick_params(axis='both', which='major', labelsize=12, top=5)
    ax2.tick_params(axis='both', which='major', labelsize=12, top=5)
    #  ax3.tick_params(axis='both', which='major', labelsize=12, top=5)

    #label axes
    ax1.set_xlabel(r'$r$ $[r_{max} = 1]$', 
            labelpad=10, 
            family='serif', 
            size=16)
    ax1.set_ylabel(r'density', 
            labelpad=10, 
            family='serif', 
            size=16)
    ax2.set_xlabel(r'$r$ $[r_{max} = 1]$', 
            labelpad=10, 
            family='serif', 
            size=16)
    ax2.set_ylabel(r'mass $[M_{tot}=10^6]$', 
            labelpad=10, 
            family='serif', 
            size=16)
    #  ax3.set_xlabel(r'$r$ $[r_{max} = 1]$',
    #          labelpad=10,
    #          family='serif',
    #          size=16)
    #  ax3.set_ylabel(r'cumulative mass $[M_{tot}=10^6]$',
    #          labelpad=10,
    #          family='serif',
            #  size=16)

    #titles
    ax1.set_title(r'density profile', 
            family='serif', 
            size=18)
    ax2.set_title('mass profile', 
            family='serif', 
            size=18)
    #  ax3.set_title('cumulative mass profile',
    #          family='serif',
    #          size=18)


    # set log scale
    ax1.set_yscale('log')
    ax1.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_xscale('log')
    #  ax3.set_yscale('log')
    #  ax3.set_xscale('log')

    # set axis limits
    #  ax1.set_ylim([1e-2, 1e18])
    #  ax2.set_ylim([1e-3, 1e2])
    #  ax2.set_ylim([1e-3, 1e2])


    # grid
    ax1.grid()
    ax2.grid()
    #  ax3.grid()

    # legend
    ax1.legend()
    ax2.legend()
    #  ax3.legend()


    plt.tight_layout()

    # annotate figure with bin number
    plt.figtext(.02, .03, str(bins.shape[0]-1)+' bins used', family='serif', size=12)



    #----------------------
    # Save figure
    #----------------------

    workdir= str(getcwd())
    fig_path = workdir+'/'+outputfilename+'.png'
    print("saving density profile plot as "+fig_path)
    plt.savefig(fig_path, format='png', facecolor=fig.get_facecolor(), transparent=False, dpi=300)#,bbox_inches='tight' )
    print("Done density profile plot")

    plt.close()

    return
   



    








#=============================
if __name__ == "__main__":
#=============================


    if len(argv) != 2:
        print("I expect exactly 1 argument: the data file.")
        quit()


    filename = argv[1]

    # read data
    mass, x, y, z = readdata(filename)
   
    # get density and mass profile
    bins, density_profile, mass_profile = get_profiles(x, y, z, mass)



    # outputfilename = 'step0_3dplot'
    # plot3d(x, y, z, outputfilename)

    outputfilename = 'step0_density_plot'
    plot_profiles(bins, density_profile, mass_profile, outputfilename)
    
    
