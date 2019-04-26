#======================================
# Cluster Repositioning
#======================================
# This code is designed to translate the initial cluster to an optimal distance
# from the double square MgO substrate placed at z = 0.
# Run this code before starting a simulation with a substrate present.
# Works for both monometallic and bimetallic clusters.
#======================================
import os
import math
import numpy
import sys
from collections import defaultdict




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


###############################
# Checking file consistency   #
###############################
#Dictionary of atom radii of different elements
atom_rad = {'Pt': 1.385, 'Ag':1.44, 'Pd':1.375, 'Au':1.46}


N= int(f_init.readline())            # Number of atoms
f_elem = f_init.readline()        # Chemical species line
f_elem = f_elem.split()           #[elem1, elem2]

substrate_pot_file.readline() #text lines 
substrate_pot_file.readline()

line= substrate_pot_file.readline()
p_elem = line.split()


#print (('%s, %s')%(f_elem, p_elem))
###
if (f_elem[0] not in p_elem or f_elem[1] not in p_elem):           #Returns FALSE when f_elem == p_elem
    print ('Type of cluster in position file and potential file do not match')
    sys.exit

################################
# Defining the input vectors   #
################################

atom = defaultdict(str)   #Type of elements
NN = [0 for x in range(N)]     # Number of Nearest Neighbours
init_x = defaultdict(float)    # Initial positions in X
init_y = defaultdict(float)    # Initial positions in Y
init_z = defaultdict(float)    # Initial positions in Z

final_z = defaultdict(float)   # Final positions in Z


################################
#      CCC coeff reading       #
################################
ccc = [[[[0 for i in range(3)] for j in range(3)] for k in range(3)] for l in range(2)]
#ccc[:1][:2][:2][:2] is the upper bounds for each dimension, ccc[0] is the 27 parameters for the first metal, ccc[2] is for the second
# I created it to have 2 dimensions although for the mono-metallic case I only use one as the two will be identical

substrate_pot_file.readline()
substrate_pot_file.readline()

for i in range(0,3):
	for j in range(0,3):
		line = substrate_pot_file.readline()
		columns = line.split()
		#print(columns[0].replace('d0',''))
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

#######################################
#       Bimetallic procedure          #
#######################################    
if (f_elem[0] != f_elem[1]):  #Bimetallic cluster
	print ('Bimetallic cluster detected')
    
	        
	atomic_radius = {}
	for a in p_elem:
		if a in atom_rad:
			atomic_radius[a] = atom_rad[a]
		else:
			print('Unrecognised atom type, check input file and supported elements')
			sys.exit()
	print ('Element 1: %s	radius: %s'%(p_elem[0], atomic_radius[p_elem[0]]))
	print ('Element 2: %s	radius: %s'%(p_elem[1], atomic_radius[p_elem[1]]))	
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

	print('Initial positions of cluster have been read successfully')
	print()

	if N != (j):
		print ('Number of atoms not equal with number of positions given')
		sys.exit
		
	################################
	#       NN calculation         #
	################################
	nn_cutoff = {}
	for a in atomic_radius:
		nn_cutoff[a] = atomic_radius[a] *2*math.sqrt(2)*0.8   #cut-off distances used. Atomic Radius *2*sqrt(2)*80%. This allows for some relaxation, but we do not consider 2nd order NN.	
	
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

	print('NN has been determined')
	print()		
		
		
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
	print('First fitting parameter triplet for %s:'%(p_elem[0]))
	print('%s, %s, %s'%(ccc[0][0][0][0], ccc[0][0][0][1], ccc[0][0][0][2]))

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
	print('First fitting parameter triplet for %s'%(p_elem[1]))
	print('%s, %s, %s'%(ccc[1][0][0][0], ccc[1][0][0][1], ccc[1][0][0][2]))

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
	print('minimum pot dist1: %s\n minimum pot dist2: %s\n minz: %s\n'%(mgo_min_dist1,mgo_min_dist2, min_z))

	if mgo_min_dist1 > mgo_min_dist2:
		reposition_dist = mgo_min_dist1 - min_z
		print('minimum MgO potential distance of element 1 [%s] used'%(p_elem[0]))
	else:
		reposition_dist = mgo_min_dist2 - minz
		print('minimum MgO potential distance of element 2 [%s] used'%(p_elem[1]))
		
	# Now we have the mgo_min_dist in [A]
	# All that is left to do is to reposition the substrate
	print('Minima of Force Field has been calculated')
	print('Minimum height from cluster has been determined')
	print()  

	################################
	#        Repositioning         #
	################################

	reposition_dist += 0.3     #Enable a soft-landing
	print('repositioning distance: %s'%(reposition_dist))

	for i in range(N):
		final_z [i] = init_z[i] + reposition_dist

	print('Cluster has been repositioned')
	print()
	#### Remember these are all in units of A #########

	################################
	#            Output            #
	################################

	f_final.write('         '+str(N) +' '+'\n')
	f_final.write('%s %s \n'%(p_elem[0], p_elem[1]))

	for i in range(N):
	    f_final.write( str(atom[i])+'  '+str('{0:.5f}'.format(init_x[i])).rjust(10) + '  '+str('{0:.5f}'.format(init_y[i])).rjust(10) +'  ' + str('{0:.5f}'.format(final_z[i])).rjust(10) + '  ' + str(NN[i]) +' '+'\n')

	print('Script has finished running')        	
		
############################################################################################		
		
#######################################
#       Monometallic procedure        #
####################################### 
else:
	print ('monometallic cluster detected')
	atom_type = f_elem[0]
	atomic_radius = 0

	if atom_type in atom_rad:
		atomic_radius = atom_rad[atom_type]
	else:	
		print('Unrecognised atom type, check input file and supported elements')
		sys.exit()

	print('Element type selected:',atom_type)
	print('Atomic radius for selected species is:',atomic_radius)
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
	#   Force field calculation    #
	################################
	ccA = [0 for i in range(3)]
	ccB = [[0 for x in range(0,3)] for y in range(0,3)]
	t3  = [0 for i in range(3)]


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
			
	print('minimum MgO pot dist: %s\n minimum z: %s'%(mgo_min_dist, min_z))
	reposition_dist = mgo_min_dist - min_z
	# Now we have the mgo_min_dist in [A]
	# All that is left to do is to reposition the substrate
	print('Minima of Force Field has been calculated')
	print('Minimum height from cluster has been determined')
	print()          

	################################
	#        Repositioning         #
	################################

	reposition_dist += 0.3     #Enable a soft-landing
	print('repositioning distance: %s'%(reposition_dist))

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
