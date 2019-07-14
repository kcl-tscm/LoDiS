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

#######################################
#       Bimetallic procedure          #
#######################################
def bim_repos(p_elem, atom_rads, atom, f_init, init_x,
			  init_y, init_z, N, NN, aa, ccc,
			  final_z, f_final):

	print ('cluster_repos_func> Element 1: %s	radius: %s'%(p_elem[0], atom_rads[p_elem[0]]))
	print ('cluster_repos_func> Element 2: %s	radius: %s'%(p_elem[1], atom_rads[p_elem[1]]))
	###################################
	#     Reading initial positions    #
	###################################
	j=0
	
	for line in f_init:
		col = list(line.split())   # used to split the input string data
		atom[j] = str(col[0])
		init_x[j] =  float(col[1])
		init_y[j] =  float(col[2])
		init_z[j] =  float(col[3]) 
		j += 1

	f_init.close()    

	print('cluster_repos_func> Initial positions of cluster have been read successfully')
	print()

	if N != (j):
		print ('cluster_repos_func> Number of atoms not equal with number of positions given')
		print()
		sys.exit()
		
	################################
	#       NN calculation         #
	################################
	nn_cutoff = {}
	for a in atom_rads:
		nn_cutoff[a] = atom_rads[a] *2*math.sqrt(2)*0.8   #cut-off distances used. Atomic Radius *2*sqrt(2)*80%. This allows for some relaxation, but we do not consider 2nd order NN.
	
	nn_dist= 0
	NN_elem = {p_elem[0]:defaultdict(int), p_elem[1]:defaultdict(int)}     #dictionary of dictionaries for each element containing atom index and NN
	for i in range(0,N-1):
		for j in range(i+1,N):
			nn_dist = math.sqrt((init_x[i]-init_x[j])**2 +(init_y[i]-init_y[j])**2+(init_z[i]-init_z[j])**2)
			if atom[i] == atom[j]:                    
				if (nn_dist < nn_cutoff[atom[i]]):
					NN[i] = NN[i]+1
					NN[j] = NN[j]+1
					NN_elem[atom[i]][i] += 1               #Of the same element so fall into the same inner dictionary
					NN_elem[atom[i]][j] += 1
			else:                                                 #If atom i and j are of different elements
				if nn_cutoff[atom[i]] > nn_cutoff[atom[j]]:        #use the larger cut-off distance for NN
					if nn_dist < nn_cutoff[atom[i]]:
						NN[i] += 1
						NN[j] += 1
						NN_elem[atom[i]][i] += 1
						NN_elem[atom[j]][j] += 1
				else:
					if nn_dist < nn_cutoff[atom[j]]:
						NN[i] += 1
						NN[j] += 1
						NN_elem[atom[i]][i] += 1
						NN_elem[atom[j]][j] += 1				

		
		
	################################
	#   Force field calculation    #
	################################
	#Note p_elem[i] corresponds to ccc[i]
	
	ccA = [[0 for i in range(3)] for j in range(2)]     #[ccA for element 1, ccA for element 2]
	ccB = [[[0 for x in range(0,3)] for y in range(0,3)] for z in range(2)] #[ccB for element 1, ccB for element 2]
	t3  = [[0 for i in range(3)] for j in range(2)]    #[t3 for element 1, t3 for element 2]

	min_z = min(list(init_z.values())) #minimum z for cluster
	
	########## Element 1 ##########
	
	elem1 = p_elem[0]
	max_coord1 = 0              # max CNi for element 1
	m1 = 0                       #index of atom with largest CNi for element 1
	for key in NN_elem[elem1]:
		if max_coord1 < NN_elem[elem1][key]:
			max_coord1 = NN_elem[elem1][key]
			m1 = key
	
	t1_1 = (math.cos(aa * init_x[m1]) + math.cos(aa * init_y[m1]))
	t2_1 = (math.cos(aa * (init_x[m1] + init_y[m1]) ) + math.cos(aa * (init_x[m1]- init_y[m1])))

	for i in range(0,3):
		for j in range(0,3):
			ccB[0][i][j] = ccc[0][i][j][0] + ccc[0][i][j][1] *t1_1 + ccc[0][i][j][2]*t2_1            

	t3[0][0] = math.exp(-max_coord1/ccB[0][0][2])
	t3[0][1] = math.exp(-max_coord1/ccB[0][1][2])
	t3[0][2] = math.exp(-max_coord1/ccB[0][2][2])

	ccA[0][0] = ccB[0][0][0] + ccB[0][0][1] * t3[0][0]   
	ccA[0][1] = ccB[0][1][0] + ccB[0][1][1] * t3[0][0]  
	ccA[0][2] = ccB[0][2][0] + ccB[0][2][1] * t3[0][0]            
            
	ener_mgo_fixed1 = [0 for i in range(501)]    # store the values of the potential at different heights
	ener_mgo_min1 = 0                                # the minimum value of the potential
	mgo_min_dist1 = 0                                # distance where whe minimum is for max_coord (distance from the MgO substrate in [A])

	for i in range(100,501):
        
		t4_1 = ( (float(i)/100) - ccA[0][2] )
		t5_1 = math.exp(-ccA[0][1]*t4_1)
		t6_1 = t5_1*t5_1

		ener_mgo_fixed1[i] = ccA[0][0] * (t6_1 - 2. * t5_1)
    
	ener_mgo_min1 = ener_mgo_fixed1[101]

	for i in range(101,501):
		if (ener_mgo_fixed1[i] < ener_mgo_min1):
			ener_mgo_min1 = ener_mgo_fixed1[i]
			mgo_min_dist1 = (float(i)/100 )

	##########Element 2##########
	
	elem2 = p_elem[1]
	max_coord2 = 0              # max CNi for element 1
	m2 = 0                       #index of atom with largest CNi for element 1
	for key in NN_elem[elem2]:
		if max_coord2 < NN_elem[elem2][key]:
			max_coord2 = NN_elem[elem2][key]
			m2 = key
	
	t1_2 = (math.cos(aa * init_x[m2]) + math.cos(aa * init_y[m2]))
	t2_2 = (math.cos(aa * (init_x[m2] + init_y[m2]) ) + math.cos(aa * (init_x[m2]- init_y[m2])))


	for i in range(0,3):
		for j in range(0,3):
			ccB[1][i][j] = ccc[1][i][j][0] + ccc[1][i][j][1] *t1_2 + ccc[1][i][j][2]*t2_2            

	t3[1][0] = math.exp(-max_coord2/ccB[1][0][2])
	t3[1][1] = math.exp(-max_coord2/ccB[1][1][2])
	t3[1][2] = math.exp(-max_coord2/ccB[1][2][2])

	ccA[1][0] = ccB[1][0][0] + ccB[1][0][1] * t3[1][0]   
	ccA[1][1] = ccB[1][1][0] + ccB[1][1][1] * t3[1][0]  
	ccA[1][2] = ccB[1][2][0] + ccB[1][2][1] * t3[1][0]            
            
	ener_mgo_fixed2 = [0 for i in range(501)]    # store the values of the potential at different heights
	ener_mgo_min2 = 0                                # the minimum value of the potential
	mgo_min_dist2 = 0                                # distance where whe minimum is for max_coord (distance from the MgO substrate in [A])

	for i in range(100,501):
        
		t4_2 = ( (float(i)/100) - ccA[1][2] )
		t5_2 = math.exp(-ccA[1][1]*t4_2)
		t6_2 = t5_2*t5_2

		ener_mgo_fixed2[i] = ccA[1][0] * (t6_2 - 2. * t5_2)
    
	ener_mgo_min2 = ener_mgo_fixed2[101]

	for i in range(101,501):
		if (ener_mgo_fixed2[i] < ener_mgo_min2):
			ener_mgo_min2 = ener_mgo_fixed2[i]
			mgo_min_dist2 = (float(i)/100 )
	print()
	print('cluster_repos_func> minimum pot dist1: %s\n minimum pot dist2: %s\n minz: %s\n'%(mgo_min_dist1,mgo_min_dist2, min_z))

	if mgo_min_dist1 > mgo_min_dist2:
		reposition_dist = mgo_min_dist1 - min_z
		print('cluster_repos_func> minimum MgO potential distance of element 1 [%s] used'%(p_elem[0]))
	else:
		reposition_dist = mgo_min_dist2 - min_z
		print('cluster_repos_func> minimum MgO potential distance of element 2 [%s] used'%(p_elem[1]))
		
	# Now we have the mgo_min_dist in [A]
	# All that is left to do is to reposition the substrate

	################################
	#        Repositioning         #
	################################

	reposition_dist += 0.3     #Enable a soft-landing
	print('cluster_repos_func> repositioning distance: %s'%(reposition_dist))

	for i in range(N):
		final_z [i] = init_z[i] + reposition_dist

	#### Remember these are all in units of A #########

	################################
	#            Output            #
	################################

	f_final.write('         '+str(N) +' '+'\n')
	f_final.write('%s %s \n'%(p_elem[0], p_elem[1]))

	for i in range(N):
		f_final.write( str(atom[i])+'  '+str('{0:.5f}'.format(init_x[i])).rjust(10) + '  '+str('{0:.5f}'.format(init_y[i])).rjust(10) +'  ' + str('{0:.5f}'.format(final_z[i])).rjust(10) + '  ' + str(NN[i]) +' '+'\n')

	print('cluster_repos_func> Cluster has been repositioned')
	print()
	return
############################################################################################		
		
#######################################
#       Monometallic procedure        #
####################################### 
def mon_repos(p_elem, atom_rads, f_init, atom, init_x,
			  init_y, init_z, N, NN, aa, ccc, final_z, f_final):

	atomic_radius = atom_rads[p_elem[0]]

	print('cluster_repos_func> Element: %s, radius: %s'%(p_elem[0],atomic_radius))
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

	print('cluster_repos_func> Initial positions of cluster have been read successfully')
	print()

	if N != (j) :
		print ('cluster_repos_func> Number of atoms not equal with number of positions given')
		print()
		sys.exit()
    

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

 

	################################
	#   Force field calculation    #
	################################
	ccA = [0 for i in range(3)]
	ccB = [[0 for x in range(0,3)] for y in range(0,3)]
	t3  = [0 for i in range(3)]


	#Finding atoms with lowest z

	min_z = min(list(init_z.values()))
	print ('cluster_repos_func> min_z = %s'%(min_z))
	max_coord= 1    # NN for the particle of interest

	max_coord = max(NN)
	m = NN.index(max_coord)     # index of the particle of interest
	t1 = (math.cos(aa * init_x[m]) + math.cos(aa * init_y[m]))
	t2 = (math.cos(aa * (init_x[m] + init_y[m]) ) + math.cos(aa * (init_x[m]- init_y[m])))


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
			
	print('cluster_repos_func> minimum MgO pot dist: %s\n minimum z: %s'%(mgo_min_dist, min_z))
	reposition_dist = mgo_min_dist - min_z
	# Now we have the mgo_min_dist in [A]
	# All that is left to do is to reposition the substrate
	print()          

	################################
	#        Repositioning         #
	################################

	reposition_dist += 0.3     #Enable a soft-landing
	print('cluster_repos_func> repositioning distance: %s'%(reposition_dist))

	for i in range(N):
		final_z [i] = init_z[i] + reposition_dist


	#### Remember these are all in units of A #########

	################################
	#            Output            #
	################################

	f_final.write('         '+str(N) +' '+'\n')
	f_final.write('%s %s \n'%(p_elem[0], p_elem[0]))

	for i in range(N):
		f_final.write( str(p_elem[0])+'  '+str('{0:.5f}'.format(init_x[i])).rjust(10) + '  '+str('{0:.5f}'.format(init_y[i])).rjust(10) +'  ' + str('{0:.5f}'.format(final_z[i])).rjust(10) + '  ' + str(NN[i]) +' '+'\n')


	print('cluster_repos_func> Cluster has been repositioned')
	print()
	return
