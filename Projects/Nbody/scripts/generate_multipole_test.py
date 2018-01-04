#!/usr/bin/python3



#==================================================
# Generates simple test data in the same format
# as given data.
# First particle is test particle, others are in a
# rectangle to check multipoles.
#==================================================


import numpy as np
from sys import argv
import random



npart = 1000 # particles in rectangle


#=================
# Read in args
#=================

if (len(argv) > 1):
    try:
        npart = int(argv[1])

    except ValueError:
        print("Expecting integer value for number of particles in each dimension")
        quit(2)






#=================
# Generate cube
#=================


x = np.zeros(npart+1,dtype='float')
y = np.zeros(npart+1,dtype='float')
z = np.zeros(npart+1,dtype='float')

random.seed(1)
for i in range(0,npart):
    x[i+1] = random.uniform(-100,100)
    y[i+1] = random.uniform(-10,10)
    z[i+1] = random.uniform(10,20)


x[0] = -1
y[0] = -1
z[0] = -1



#===================================
# Generate other values for read-in
#===================================

counter = npart + 1

m = [1]*counter

vx = [1]*counter
vy = [1]*counter
vz = [1]*counter

softening = [1]*counter
potential = [1]*counter









#====================================
# Write results to file
#====================================
fname = "testdata_multipole.ascii"

f = open(fname, "w")

header = '{0:d} {1:d} {0:d}\n'.format(counter, 0)
f.write(header) 

data = [m, x, y, z, vx, vy, vz, softening, potential]

for d in data:
    for i in range(counter):
        f.write('{0:f}\n'.format(d[i]))




# plot if necessary
#======================
def plotme(x, y, z):
#======================
    import matplotlib.pyplot as plt 
    from mpl_toolkits.mplot3d import Axes3D

    
    fig = plt.figure(facecolor = 'white', figsize = (10, 10))
    ax1 = fig.add_subplot(1,1,1, aspect='equal', projection='3d')

    ax1.scatter3D(x,y,z)
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    ax1.set_zlabel('z')
    plt.show()


plotme(x,y,z)
