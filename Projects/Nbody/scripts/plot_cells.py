#!/usr/bin/python3

#========================================================
# Plot the particles according to which cell they belong
# Needs cellparticles.dat
#========================================================


from os import getcwd, listdir
from sys import argv
import numpy as np
import random
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D






#===============================================================
def plot2d(ax1, ax2, ax3, x, y, z, cell, cellcentres=True):
#===============================================================


    ncells = int(cell.max())+1

    colorcounter=0
    global mycolormap
    global indarray
    global indarray2

    if (cellcentres):
        ax1.scatter(x,y,c='black', s=10, label='cell centres')
        ax2.scatter(y,z,c='black', s=10, label='cell centres')
        ax3.scatter(x,z,c='black', s=10, label='cell centres')

    else:
        for i in range(ncells):

            xt = x[cell == i]
            yt = y[cell == i]
            zt = z[cell == i]

            if xt.shape[0]>0:
                
                ind = indarray[colorcounter]
                ind2 = indarray2[colorcounter]
                
                ax1.scatter(xt, yt, 
                        facecolor=mycolormap(ind), 
                        s = 15,
                        linewidth=1, 
                        edgecolor=mycolormap(ind2))
                ax2.scatter(yt, zt, 
                        facecolor=mycolormap(ind), 
                        s = 15,
                        linewidth=1, 
                        edgecolor=mycolormap(ind2))
                ax3.scatter(xt, zt, 
                        facecolor=mycolormap(ind), 
                        s = 15,
                        linewidth=1, 
                        edgecolor=mycolormap(ind2))
                colorcounter += 1


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

    ax1.legend()
    ax2.legend()
    ax3.legend()


    ax1.set_xlim((-1.1,1.1))
    ax1.set_ylim((-1.1,1.1))
    ax2.set_xlim((-1.1,1.1))
    ax2.set_ylim((-1.1,1.1))
    ax3.set_xlim((-1.1,1.1))
    ax3.set_ylim((-1.1,1.1))

    return 









#===============================================
def plot3d(ax1, x, y, z, cell, cellcentres):
#===============================================

    ncells = int(cell.max())+1

    colorcounter = 0
    global mycolormap
    global indarray
    global indarray2



    for i in range(ncells):

        if (cellcentres):
            ax1.scatter3D(x, y, z, 
                    depthshade=False, 
                    label='cell centres', 
                    c='black',
                    s=5)

            
        else:
            xt = x[cell == i]
            yt = y[cell == i]
            zt = z[cell == i]

            if (xt.shape[0]>0):
                ind = indarray[colorcounter]
                ind2 = indarray2[colorcounter]

                ax1.scatter3D(xt, yt, zt, 
                        depthshade=False, 
                        #  label='cell '+str(i),
                        facecolor=mycolormap(ind), 
                        s = 15,
                        linewidth=1, 
                        edgecolor=mycolormap(ind2))
                colorcounter += 1


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
    #  ax1.legend()

    ax1.set_xlim((-1.1,1.1))
    ax1.set_ylim((-1.1,1.1))
    ax1.set_zlim((-1.1,1.1))





    return 


   






#=========================
if __name__=="__main__":
#=========================


    #==========================
    # Preparation
    #==========================

    if (len(argv) != 3):
        print("Expecting exactly 2 arguments: \nthe output file of which cell each particle belongs to \nand cell centres.")
        quit()





    #===================================
    # Get particles by cell
    #===================================

    filename = argv[1]

    xp, yp, zp = np.loadtxt(filename, skiprows=1, usecols=([0, 1, 2]), unpack=True)
    cellp = np.loadtxt(filename, skiprows=1, dtype='int', usecols=([3])) 
    ncells = int(cellp.max()) - int(cellp.min()) + 1
    print("Found", ncells, "cells")







    #===================================
    # Get cell centres
    #===================================

    filename = argv[2]

    xc, yc, zc = np.loadtxt(filename, skiprows=1, usecols=([0, 1, 2]), unpack=True)
    cellc = np.loadtxt(filename, skiprows=1, dtype='int', usecols=([3])) 





    #=======================
    # Plotting prep
    #=======================


    # get default colormap
    mycolormap = plt.cm.get_cmap('tab20', ncells)
    
    # get random index to mix up colors
    random.seed(19)
    indarray = random.sample(range(ncells), ncells)
    indarray2 = random.sample(range(ncells), ncells)







    #=====================
    # Plot and save 2D
    #=====================

    fig2 = plt.figure(facecolor = 'white', figsize = (20, 8))
    ax1 = fig2.add_subplot(1,3,1, aspect='equal')
    ax2 = fig2.add_subplot(1,3,2, aspect='equal')
    ax3 = fig2.add_subplot(1,3,3, aspect='equal')

    plot2d(ax1, ax2, ax3, xp, yp, zp, cellp, cellcentres=False)
    plot2d(ax1, ax2, ax3, xc, yc, zc, cellc, cellcentres=True)

    plt.tight_layout()


    workdir= str(getcwd())
    outputfilename = 'cell_plot'
    fig_path = workdir+'/'+outputfilename+'.png'
    print("saving cellcentres/particle 2D plot as "+fig_path)
    plt.savefig(fig_path, format='png', facecolor=fig2.get_facecolor(), transparent=False, dpi=300)#,bbox_inches='tight' )
    print("Done")

    plt.close()








    #===================
    # Plot and save 3D
    #===================


    fig1 = plt.figure(facecolor = 'white', figsize = (10, 10))
    ax0 = fig1.add_subplot(1,1,1, aspect='equal', projection='3d')

    plot3d(ax0, xp, yp, zp, cellp, cellcentres=False)
    plot3d(ax0, xc, yc, zc, cellc, cellcentres=True)


    plt.tight_layout()
    wokdir= str(getcwd())


    outputfilename = 'cell_plot-3D'
    fig_path = workdir+'/'+outputfilename+'.png'
    print("saving cellcentres/particle 3D plot as "+fig_path)
    fig1.savefig(fig_path, format='png', facecolor=fig1.get_facecolor(), transparent=False, dpi=300)#,bbox_inches='tight' )
    print("Done")

    plt.close()







