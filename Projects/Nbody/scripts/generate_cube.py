#!/usr/bin/python3



#==================================================
# Generates simple test data in the same format
# as given data.
# This gives particles in a simple cube.
#==================================================


import numpy as np
from sys import argv




npart = 2 # particles in each dimension
step = 1


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


x = np.zeros(npart**3,dtype='float')
y = np.zeros(npart**3,dtype='float')
z = np.zeros(npart**3,dtype='float')
counter = 0

for i in range(0,npart, step):
    for j in range(0, npart, step):
        for k in range(0, npart, step):
            x[counter] = i*1.0
            y[counter] = j*1.0
            z[counter] = k*1.0
            counter += 1

# norm values the same way you do in nbody code
radius = np.sqrt(x[-1]**2 + y[-1]**2 + z[-1]**2) * 0.5

x = x/radius
y = y/radius
z = z/radius

l = x[-1]-x[0]
x = x - l/2.0
y = y - l/2.0
z = z - l/2.0



print("Check 1:", x[0], y[0], z[0], np.sqrt(x[0]**2 + y[0]**2 + z[0]**2))
print("Check 2:", x[-1], y[-1], z[-1])





#===================================
# Generate other values for read-in
#===================================
m = [1]*counter

vx = [1]*counter
vy = [1]*counter
vz = [1]*counter

softening = [1]*counter
potential = [1]*counter









#====================================
# Write results to file
#====================================
fname = "testdata_cube"+str(npart**3)+".ascii"

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


#  plotme(x,y,z)
