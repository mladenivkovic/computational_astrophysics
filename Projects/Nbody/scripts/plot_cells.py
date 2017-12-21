#!/usr/bin/python3

#========================================================
# Plot the particles according to which cell they belong
#========================================================


from os import getcwd, listdir
from sys import argv
import numpy as np
import matplotlib.pyplot as plt



#==============================
def plot2d(x, y, z, cells):
#==============================

    ncells = int(cell.max())+1


    fig = plt.figure(facecolor = 'white', figsize = (20, 8))
    ax1 = fig.add_subplot(1,3,1, aspect='equal')
    ax2 = fig.add_subplot(1,3,2, aspect='equal')
    ax3 = fig.add_subplot(1,3,3, aspect='equal')


    tot = 0
    for i in range(ncells):

        xt = x[cell == i]
        yt = y[cell == i]
        zt = z[cell == i]

        ax1.scatter(xt, yt)
        ax2.scatter(yt, zt)
        ax3.scatter(xt, zt)


    ax1.set_xlabel('x', 
            labelpad=10, 
            family='serif', 
            size=16)

    ax1.set_ylabel('y', 
            labelpad=10, 
            family='serif', 
            size=16)

    ax2.set_xlabel('y', 
            labelpad=10, 
            family='serif', 
            size=16)

    ax2.set_ylabel('z', 
            labelpad=10, 
            family='serif', 
            size=16)

    ax3.set_xlabel('x', 
            labelpad=10, 
            family='serif', 
            size=16)

    ax3.set_ylabel('z', 
            labelpad=10, 
            family='serif', 
            size=16)




    ax1.grid()
    ax2.grid()
    ax3.grid()

    ax1.set_xlim((-1,1))
    ax1.set_ylim((-1,1))
    ax2.set_xlim((-1,1))
    ax2.set_ylim((-1,1))
    ax3.set_xlim((-1,1))
    ax3.set_ylim((-1,1))


    plt.tight_layout()


    workdir= str(getcwd())
    outputfilename = 'particles_by_cell_plot'
    fig_path = workdir+'/'+outputfilename+'.png'
    print("saving particle cell plot as "+fig_path)
    plt.savefig(fig_path, format='png', facecolor=fig.get_facecolor(), transparent=False, dpi=300)#,bbox_inches='tight' )
    print("Done")

    plt.close()


    return 


#==============================
def plot3d(x, y, z, cells):
#==============================

    from mpl_toolkits.mplot3d import Axes3D
    ncells = int(cell.max())+1


    fig = plt.figure(facecolor = 'white', figsize = (10, 10))
    ax1 = fig.add_subplot(1,1,1, aspect='equal', projection='3d')


    tot = 0
    for i in range(ncells):

        xt = x[cell == i]
        yt = y[cell == i]
        zt = z[cell == i]

        ax1.scatter(xt, yt, zt, depthshade=False, label='cell '+str(i))


    ax1.set_xlabel('x', 
            labelpad=10, 
            family='serif', 
            size=16)

    ax1.set_ylabel('y', 
            labelpad=10, 
            family='serif', 
            size=16)

    ax1.set_zlabel('z', 
            labelpad=10, 
            family='serif', 
            size=16)



    ax1.grid()
    ax1.legend()

    ax1.set_xlim((-1,1))
    ax1.set_ylim((-1,1))
    ax1.set_zlim((-1,1))


    plt.tight_layout()


    workdir= str(getcwd())
    outputfilename = 'particles_by_cell_plot-3D'
    fig_path = workdir+'/'+outputfilename+'.png'
    print("saving particle cell plot as "+fig_path)
    plt.savefig(fig_path, format='png', facecolor=fig.get_facecolor(), transparent=False, dpi=300)#,bbox_inches='tight' )
    print("Done")

    plt.close()


    return 


    
#=========================
if __name__=="__main__":
#=========================


    if (len(argv) != 2):
        print("Expecting exactly 1 argument: the output file.")
        quit()

    filename = argv[1]

    x, y, z = np.loadtxt(filename, skiprows=1, usecols=([0, 1, 2]), unpack=True)
    cell = np.loadtxt(filename, skiprows=1, dtype='int', usecols=([3])) 
    ncells = int(cell.max())+1
    print("Found", ncells, "cells")

    plot2d(x, y, z, cell)
    plot3d(x, y, z, cell)


