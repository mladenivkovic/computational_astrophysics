#!/usr/bin/python3



#==================================================
# Generates simple test data in the same format
#==================================================



import random 

npart = 16
minval = -10
maxval = 10

m = [abs(random.uniform(minval, maxval))*100 for i in range(npart)]
x = [random.uniform(minval, maxval) for i in range(npart)]
y = [random.uniform(minval, maxval) for i in range(npart)]
z = [random.uniform(minval, maxval) for i in range(npart)]

vx = [random.uniform(minval, maxval) for i in range(npart)]
vy = [random.uniform(minval, maxval) for i in range(npart)]
vz = [random.uniform(minval, maxval) for i in range(npart)]

softening = [random.uniform(minval, maxval) for i in range(npart)]
potential = [random.uniform(minval, maxval) for i in range(npart)]



f = open("testdata.ascii", "w")

header = '{0:d} {1:d} {0:d}\n'.format(npart, 0)
f.write(header) 

data = [m, x, y, z, vx, vy, vz, softening, potential]

for d in data:
    for i in range(npart):
        f.write('{0:f}\n'.format(d[i]))
