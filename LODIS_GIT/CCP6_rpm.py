# -*- coding: utf-8 -*-
"""
Created on: Sat Jun 17 11:04:59 2017
Final version: Sat Sep 9 20017

@author: NZ

Disclamer: Please note the code is not perfectly optimized and further improvements are possible.

# Variable are named in an intuitive way to help understand the code.

# Please note the code ONLY works for mono-metalic elements. 
It can be further developed to be able to support bi-metallic clusters. The work to be done added is not much.
If you are interested in having this developed and require guidance please get in touch.

***** I suggest adding a 'conversion' array that will assign the type of element to a give index [atom_type becomes array], this way a numerical value can be assigned to each type of element used, allowing for 'safety' checks.
I have created the CCC coefficient array in such a way to allow for easy implementation for the second set of coefficients as now we have dimensions [2,3,3,3] where the first index will allow for second element type.

# Double check the distribution of NN defined by code arithmetically with the one from PDF  #

# Tests confirm the accurate values, precaution is advised due to the importance of the cut-off distance

# Further optimization for cut-off distance is required to improve on the code

# Possible memory saving improvement and speed boost can be achieved by moving from a defaultdict to a specific size array

Contact:    Norbert.Zicher@physics.ox.ac.uk
"""
import os
import math
import numpy
import sys
from collections import defaultdict

"""
Variable list:
  Data Files:    f_init, f_final, g_rdf, g_pdf, substrate_pot_file
  Variables :    atom, atom_type, atomic_radius
                 i, j, k, l, m
                 N, NN, NN_cutoff, NN_dist
                 line, columns
                 ccc, ccA, ccB
                 metal_1, metal_2, aa, min_z, max_coord
                 t1, t2, t3, t4, t5, t6
                 ener_mgo_fixed, ener_mgo_min, mgo_min_dist
                 reposition_dist
                 
                 
  Position: init_x, init_y, init_z, final_z
            average_x, average_y, average_z
            dist_x, dist_y, dist_z, dist_tot     
    
"""

################################
#      Initiate files          #
################################
f_init = open(raw_input('Give the .xyz file name to be read: '),'r')      # Original Structure
substrate_pot_file = open (raw_input('Give the Mg0 pot file: '),'r')     # File containing Substrate parameters

outpath = raw_input('give the output directory path: ')
f_final  = open(os.path.join(outpath,raw_input('Give the output file name for the .xyz: ')),'w+')     # Final Structure

################################
#       Initialization         #
################################
print('===============')
print('Starting script')
print('===============')
print('Please ensure the correct path is provided to the xyz file of the structure')
print('Additionally, please make sure that the input file is of the required format')
print('-----   -----   -----   -----   -----   -----   -----   -----')
print()


N= int(f_init.readline())            # Number of atoms
f_elem = f_init.readline()        # Chemical species line
f_elem = f_elem.split()           #[elem1, elem2]

substrate_pot_file.readline() #text lines 
substrate_pot_file.readline()

line= substrate_pot_file.readline()
p_elem = line.split()


print (('%s, %s')%(f_elem, p_elem))
###           
if (f_elem[0] != f_elem[1]) or (p_elem[0] != p_elem[1]):
    print ('Type of elements different. Only mono atomic case supported')
    sys.exit

if (f_elem[0] != p_elem[0]):
    print ('Type of element in position file and potential file do not match')
    sys.exit
 

atom_type = str(f_elem[0])
atomic_radius = 0

if atom_type == 'Pt' :
    atomic_radius= 1.385
elif atom_type == 'Ag' :
    atomic_radius= 1.44
elif atom_type == 'Pd' :
    atomic_radius = 1.375
elif atom_type == 'Au':
	atomic_radius = 1.46
else:
    print('Unrecognized atom type,check input file and supported elements')
    sys.exit

print('Element type selected:',atom_type)
print('Atomic radius for selected species is:',atomic_radius)
print()
""" 
The above logical gate can be extended to incorporate new elements as the experimental data for the coefficients is published/determined.
This can be done in a rather straighforward mode by just adding the following lines before the last ### elif ### statement

########################
elif atom_type == 'XXXXXX' :
    atomic_radius= XXXXXX
########################

Please replace XXXXXX with the desired values. """

i=0
j=0
k=0
l=0

################################
# Defining the input vectors   #
################################

atom = defaultdict(str)   #Type of elements
NN = [0 for x in range(N)]     # Number of Nearest Neighbours

init_x = defaultdict(float)    # Initial positions in X
init_y = defaultdict(float)    # Initial positions in Y
init_z = defaultdict(float)    # Initial positions in Z

final_z = defaultdict(float)   # Final positions in Z
param = defaultdict(float)    #Input parameters from .pot file

print('Variable initialization completed')
print()

################################
#   Reading initial position   #
################################
j=0
for line in f_init:
    columns = line.split()   # used to split the input string data
    atom[j] = columns[0]     #initialing the input data from file
    init_x[j] =  float(columns[1])
    init_y[j] =  float(columns[2])
    init_z[j] =  float(columns[3])
    j=j+1

f_init.close()    

print('Initial positions of cluster have been read successfully')
print()

if N != (j) :
    print ('Number of atoms not equal with number of positions given')
    sys.exit
    

################################
#       NN calculation         #
################################

nn_cutoff = atomic_radius *2*math.sqrt(2)*0.8   #cut-off distance used. Atomic Radius *2*sqrt(2)*80%. This allows for some relaxation, but we do not consider 2nd order NN.
nn_dist= 0


for i in range(0,N-1):
    for j in range(i+1,N):
        nn_dist = math.sqrt((init_x[i]-init_x[j])**2 +(init_y[i]-init_y[j])**2+(init_z[i]-init_z[j])**2)
        if (nn_dist < nn_cutoff):
            NN[i] = NN[i]+1
            NN[j] = NN[j]+1

print('NN has been determined')
print()
################################
#     CCC coeff building       #
################################

ccc = [[[[0 for i in range(3)] for j in range(3)] for k in range(3)] for l in range(2)]

# I created it to have 2 dimensions although for the mono-metalic case I only use one as the two will be identical

ccA = [0 for i in range(3)]
ccB = [[0 for x in range(0,3)] for y in range(0,3)]
t3  = [0 for i in range(3)]
 
################################
#      CCC coeff reading       #
################################

substrate_pot_file.readline()
substrate_pot_file.readline()

for i in range(0,3):
    for j in range(0,3):
        line = substrate_pot_file.readline()
        columns = line.split()
        print(columns[0].replace('d0',''))
        ccc[0][i][j][0] = float(columns[0].replace('d0',''))
        ccc[0][i][j][1] = float(columns[1].replace('d0',''))
        ccc[0][i][j][2] = float(columns[2].replace('d0',''))

substrate_pot_file.readline()
substrate_pot_file.readline()            

for i in range(0,3):
    for j in range(0,3):
        line = substrate_pot_file.readline()
        columns = line.split()
        ccc[1][i][j][0] = float(columns[0].replace('d0',''))
        ccc[1][i][j][1] = float(columns[1].replace('d0',''))
        ccc[1][i][j][2] = float(columns[2].replace('d0',''))
        
substrate_pot_file.readline()
substrate_pot_file.readline()

aa_read = substrate_pot_file.readline()
aa_read = aa_read.split()

aa =float(aa_read[0].replace('d0','')) # mgo lattice size in A

print('Coefficients have been read successfully')
print('Lattice size is:',aa)
################################
#   Force field calculation    #
################################

#Finding atoms with lowest z

min_z = min(list(init_z.values()))
print ('min_z = %s'%(min_z))
max_coord= 1    # NN for the particle of interest

max_coord = max(NN)
m = NN.index(max_coord)     # index of the particle of interest
t1 = (math.cos(aa * init_x[m]) + math.cos(aa * init_y[m]))
t2 = (math.cos(aa * (init_x[m] + init_y[m]) ) + math.cos(aa * (init_x[m]- init_y[m])))
print('%s, %s, %s, %s'%(t1, ccc[0][0][0][0], ccc[0][0][0][1], ccc[0][0][0][2]))

for i in range(0,3):
    for j in range(0,3):
        ccB[i][j] = ccc[0][i][j][0] + ccc[0][i][j][1] *t1 + ccc[0][i][j][2]*t2            

t3[0] = math.exp(-max_coord/ccB[0][2])
t3[1] = math.exp(-max_coord/ccB[1][2])
t3[2] = math.exp(-max_coord/ccB[2][2])

ccA[0] = ccB[0][0] + ccB[0][1] * t3[0]   
ccA[1] = ccB[1][0] + ccB[1][1] * t3[0]  
ccA[2] = ccB[2][0] + ccB[2][1] * t3[0]            
            
ener_mgo_fixed = [0 for i in range(501)]    # store the values of the potential at different heights
ener_mgo_min = 0                                # the minimum value of the potential
mgo_min_dist = 0                                # distance where whe minimum is for max_coord (distance from the MgO substrate in [A])

for i in range(100,501):
        
    t4 = ( (float(i)/100) - ccA[2] )
    t5 = math.exp(-ccA[1]*t4)
    t6 = t5*t5

    ener_mgo_fixed[i] = ccA[0] * (t6 - 2. * t5)
    
ener_mgo_min = ener_mgo_fixed[101]

for i in range(101,501):
    if (ener_mgo_fixed[i] < ener_mgo_min):
         ener_mgo_min = ener_mgo_fixed[i]
         mgo_min_dist = (float(i)/100 )

# Now we have the mgo_min_dist in [A]
# All that is left to do is to reposition the substrate
print('Minima of Force Field has been calculated')
print('Minimum height from cluster has been determined')
print()          
################################
#        Repositioning         #
################################

reposition_dist = (mgo_min_dist - min_z)
reposition_dist = reposition_dist + 0.3     #

for i in range(N):
    final_z [i] = init_z[i] + reposition_dist

print('Cluster has been repositioned')
print()
#### Remember these are all in units of A #########

################################
#            Output            #
################################

f_final.write('         '+str(N) +' '+'\n')
f_final.write('%s %s \n'%(atom_type, atom_type))

for i in range(N):
    f_final.write( str(atom_type)+'  '+str('{0:.5f}'.format(init_x[i])).rjust(10) + '  '+str('{0:.5f}'.format(init_y[i])).rjust(10) +'  ' + str('{0:.5f}'.format(final_z[i])).rjust(10) + '  ' + str(NN[i]) +' '+'\n')

print('Script has finished running')
##########################
####   Finalization   ####
##########################
    

f_final.close()
substrate_pot_file.close()
