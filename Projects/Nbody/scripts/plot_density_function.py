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
        mass, x, y, z, vx, vy, vz, softening, potential
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
    vx = rawdata[4*npart: 5*npart]
    vy = rawdata[5*npart: 6*npart]
    vz = rawdata[6*npart: 7*npart]
    softening = rawdata[7*npart: 8*npart]
    potential = rawdata[9*npart: 10*npart]

   



    return mass, x, y, z, vx, vy, vz, softening, potential




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
def get_density_profile(x, y, z, mass):
#======================================
    """

    calculates the density profile.

    PARAMETERS:
    y, x, z, mass:  coordinates and mass of particles


    RETURNS:
    r, bins, density_profile, mass_profile:    
    r:          np.array of distance from origin of each particle
    bins:       bin distance (particles sorted in bins by distance from origin)
    density_profile,
    mass_profile :  np.array of density/mass profile
    """

    print("Calculating density profile")


    # calculate radius from origin
    r = np.sqrt(x**2 + y**2 + z**2)

    # norm r
    r = r / r.max()


    nbins = 1000
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
    mass_profile = mass_profile / np.sum(mass_profile) # norm mass profile to 1

    # get density profile
    density_profile = np.zeros((mass_profile.shape[0]))
    density_profile[1:] = mass_profile[1:] / (4 / 3 * np.pi * (bins[1:]**3-bins[:-1]**3))
    density_profile[0] = mass_profile[0]/ (4 / 3 * np.pi * bins[0]**3 )

    return r, bins, density_profile, mass_profile





#=========================================================================
def plot_density(r, bins, density_profile, mass_profile, outputfilename):
#=========================================================================
    """
    plots the density profile. 

    PARAMETERS:

        bins:   bin positions
        density_profile: density profile to plot
        outputfilename: output filename

    RETURNS:
        
        nothing
    """


    import matplotlib.pyplot as plt

    print("Creating density plot")

    #==========================
    # prepare analytical model
    #==========================

    # find a: cum_mass(a) = tot_mass/4
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



    def rho(r, a, totmass):
        return totmass/(2 * np.pi) * a / r * 1/(r + a)**3


    def M(r, a, totmass):
        return totmass * r**2/(r+a)**2




    #==================
    # Plotting
    #==================

    fig = plt.figure(facecolor='white', figsize=(20,8))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)


    ax1.loglog(bins, density_profile, 'o', markersize=3)
    ax1.plot(bins, rho(bins, a, totmass))
    ax1.set_title(r'density profile', family='serif', size=24, pad=20)

    ax2.loglog(bins, cum_mass, 'o', markersize=3)
    ax2.plot(bins, M(bins, a, totmass))
    ax2.set_title('Cumulative mass profile', family='serif', size=24)



    # set tick params (especially digit size)
    ax1.tick_params(axis='both', which='major', labelsize=12, top=5)
    ax2.tick_params(axis='both', which='major', labelsize=12, top=5)

    #label axes
    ax1.set_xlabel(r'r', labelpad=20, family='serif', size=16)
    ax1.set_ylabel(r'density', labelpad=20, family='serif', size=16)
    ax2.set_xlabel(r'r', labelpad=20, family='serif', size=16)
    ax2.set_ylabel(r'cumulative mass', labelpad=20, family='serif', size=16)

    plt.tight_layout()



    # Save figure

    workdir= str(getcwd())
    fig_path = workdir+'/'+outputfilename+'.png'
    print("saving density profile plot as "+fig_path)
    plt.savefig(fig_path, format='png', facecolor=fig.get_facecolor(), transparent=False, dpi=100)#,bbox_inches='tight' )
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
    mass, x, y, z, vx, vy, vz, softening, potential = readdata(filename)
   



    r, bins, density_profile, mass_profile = get_density_profile(x, y, z, mass)



    # outputfilename = 'step0_3dplot'
    # plot3d(x, y, z, outputfilename)

    outputfilename = 'step0_density_plot'
    plot_density(r, bins, density_profile, mass_profile, outputfilename)
    
    
