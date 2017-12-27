#!/usr/bin/python3



#==================================================
# Generates simple test data in the same format
# as given data.
# This gives particles in a simple cube.
#==================================================




npart = 8 # particles in each dimension
step = 1


x = []
y = []
z = []
counter = 0

for i in range(0,npart, step):
    for j in range(0, npart, step):
        for k in range(0, npart, step):
            x.append(i - (npart-1)/2)
            y.append(j - (npart-1)/2)
            z.append(k - (npart-1)/2)
            counter += 1




m = [1]*counter

vx = [1]*counter
vy = [1]*counter
vz = [1]*counter

softening = [1]*counter
potential = [1]*counter



f = open("testdata_cube.ascii", "w")

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
