import os
import math
import numpy
import sys
from collections import defaultdict
import cluster_repos_func as crf

def read_intpot(filepot):

    f = open(filepot, 'r+')
    potread = f.readlines()

    charrep = [',', '\t', '\n']
    for l in range(len(potread)):
        for c in charrep:
            potread[l].replace(c, ' ')
        potread[l] = potread[l].split()
    in_elem = potread[2]
    atom_r = [float(r.replace('d0',' ')) for r in potread[12][:2]]
    atom_m = [float(r.replace('d0',' ')) for r in potread[13][:2]]
    atom_rads = {k:v for k,v in zip(in_elem, atom_r)}    #Contains 1 (Mono) or 2 (Bim) key/value pairs
    atom_mass = {k:v for k,v in zip(in_elem, atom_m)}
    f.close()
    return atom_rads, in_elem, atom_mass

def cluster_repos( posfile, mgofile, intfile, mgo_pres, run_process):               #mgofile will be 'None' when mgo_pres = False

    err_dict = {1: 'Elements in .pot and .xyz do not match',
                2: 'Unknown element(s) in .xyz file not accounted for in .pot file',
                3: 'Elements in .MgO.pot and .xyz do not match'}
    err_out = 0
    f_init = open(posfile,'r')      # Original Structure


    # Read elem1, elem2 and N
    charrep = ['\t', '\n', ',']            #Clean file of any characters

    f_read = f_init.readlines()
    N_r= f_read[0].split()            # Number of atoms
    N = int(N_r[0])
    elem_list = []
    for f in range(N+2):
        for c in charrep:
            f_read[f] = f_read[f].replace(c, ' ')
        f_read[f] = f_read[f].split()
        if f >= 2:
            elem_list.append(f_read[f][0])


    f_read = f_read[2:N+2]            #Remove the headers on the first two lines
    f_elem = []
    for e in elem_list:
        if e not in f_elem:
            f_elem.append(e)
    if len(f_elem) == 1:
        f_elem = f_elem*2

    print (f_elem)
    #read interaction .pot file to get atom radii and check consistency
    atom_rads, in_elem, atom_mass = read_intpot(intfile)

    f_elem.sort()               #Alphabetical order to compare to in_elem
    incheck = in_elem
    incheck.sort()
    if run_process not in ['Coalescence', 'Growth']:                           # .pot and .xyz must be the same
        if f_elem != incheck:
            err_out = [1, err_dict[1]]
            return f_elem[0], f_elem[1], N, atom_rads, err_out, atom_mass
    else:                                                                 #Growth and coalescence can have a monometallic cluster become bimetallic
        for f in f_elem:                                                  # .xyz must be in .pot but not the other way around
            if f not in incheck:
                err_out = [2, err_dict[2]]
                return f_elem[0], f_elem[1], N, atom_rads, err_out, atom_mass


    if mgo_pres:                          #If the substrate is present carry out repositioning
        print('---------------------------------------------------')
        print('cluster_repos_prep> Beginning repositioning process')
        print('---------------------------------------------------')

        if '_repos.xyz' not in posfile:                           #Not an already repositioned file
            f_final  = open(posfile.replace('.xyz', '_repos.xyz'),'w+')     # Final Structure
        else:
            f_final = open(posfile, 'w+')                         #Will overwrite existing _repos.xyz file

        substrate_pot_file = open (mgofile,'r')     # File containing Substrate parameters
        substrate_pot_file.readline() #text lines
        substrate_pot_file.readline()

        line= substrate_pot_file.readline()
        for i in charrep:
            line.replace(i, ' ')
        p_elem = line.split()

        p_elem.sort()
        if f_elem != p_elem:
            err_out = [3, err_dict[3]]
            return f_elem[0], f_elem[1], N, atom_rads, err_out, atom_mass

        # Defining the input vectors
        atom = defaultdict(str)   #Type of elements
        NN = [0 for x in range(N)]     # Number of Nearest Neighbours
        init_x = defaultdict(float)    # Initial positions in X
        init_y = defaultdict(float)    # Initial positions in Y
        init_z = defaultdict(float)    # Initial positions in Z

        final_z = defaultdict(float)   # Final positions in Z


        #      CCC coeff reading
        ccc = [[[[0 for i in range(3)] for j in range(3)] for k in range(3)] for l in range(2)]
        #ccc[:1][:2][:2][:2] is the upper bounds for each dimension, ccc[0] is the 27 parameters for the first metal, ccc[2] is for the second

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

        print('cluster_repos_prep> Coefficients have been read successfully')

        #       Bimetallic procedure
        if (f_elem[0] != f_elem[1]):  #Bimetallic cluster
            print ('cluster_repos_prep> Bimetallic cluster detected')
            crf.bim_repos(p_elem, atom_rads, atom, f_init, init_x,
                         init_y, init_z, N, NN, aa, ccc,
                         final_z, f_final)


        #       Monometallic procedure
        else:
            print ('cluster_repos_prep> monometallic cluster detected')
            crf.mon_repos(p_elem, atom_rads, f_init, atom, init_x,
                         init_y, init_z, N, NN, aa, ccc, final_z, f_final)

        substrate_pot_file.close()
        f_final.close()


    f_init.close()
    return f_elem[0], f_elem[1], N, atom_rads, err_out, atom_mass
