
# coding: utf-8

# # Pressure map of LoDiS movie.xyz
# 
# ## Requires:
#     - Input .POT file in the standard format
#     - The movie.xyz file


import numpy as np
from itertools import groupby
from collections import namedtuple



def readMovieFileXYZ(path_to_movie):
    """
    Reads a LoDiS movie.xyz file and fetches the coordinates
    for each atom for each frame.
    
    Input:
        path_to_movie: path to movie.xyz
    
    Returns:
        Named tuple read_movie:
            - read_movie.Frames: list of frames; each is an array of atoms each described by [Atom, x, y, z, Col]
            - read_movie.Headers: list of the movie frames headers    """
    
    read_file_chars = []

    with open(path_to_movie, 'r') as file:
        for line in file:
            read_file_chars.append(line)
    # 1. Delete line jump
    read_file_chars = [line[:-1] for line in read_file_chars]
    read_file_chars

    # 2. Separate line by line
    grouped_lines = [([list(group) for k, group in groupby(line,lambda x: x == " ") if not k]) for line in read_file_chars]

    # 3. Concatenate charaters
    joined_string = [[''.join(info_elem) for info_elem in grouped_line] for grouped_line in grouped_lines]

    # 4. Regroup into list of lists. Elements of outerlist are movie frames
    merged_frames = []
    current_frame = []

    for line in joined_string:
        if(line==joined_string[0]):

            if len(current_frame)!=0:
                merged_frames.append(current_frame)
            current_frame=[]
        else:
            current_frame.append(line)
    merged_frames.append(current_frame)

    # 5. Removing second line of header
    movie_headers_all_frames = [frame[0] for frame in merged_frames]
    merged_frames = [frame[1:] for frame in merged_frames]

    # 6. Converting coordinates and pressure to floats
    for frame in merged_frames:
        for line in frame:
            line[1] = float(line[1]) # x coord
            line[2] = float(line[2]) # y coord
            line[3] = float(line[3]) # z coord
            
    Movie = namedtuple('Movie', 'Frames Headers')
    read_movie = Movie(merged_frames, movie_headers_all_frames)

    return(read_movie) 

def readPotentialFile(path_to_potfile):
    """
    Reads the .pot file and extracts all parameters from it.
    Including the potential analytical continuation.
    
    Returns:
        potential (named tuple): contains all parameters extracted:
            'AtomTypes P Q A Qsi Cohesion Radius Mass Cutoff dik0 x3 x4 x5 a3 a4 a5'
    """
    read_file_chars = []

    with open(path_to_potfile, 'r') as file:
        for line in file:
            read_file_chars.append(line)
    # Delete line jump

    # Read all values
    read_file_chars = [line[:-1] for line in read_file_chars]

    atoms_type = [elem for elem in read_file_chars[2].split(' ') if elem != '']
    p_val = np.fromstring(read_file_chars[5], sep=' ')
    q_val = np.fromstring(read_file_chars[6], sep=' ')
    a_val = np.fromstring(read_file_chars[7], sep=' ')
    qsi_val = np.fromstring(read_file_chars[8], sep=' ')
    cohe_val = np.fromstring(read_file_chars[11], sep=' ')
    atom_rad_val = np.fromstring(read_file_chars[12], sep=' ')
    mass_val = np.fromstring(read_file_chars[13], sep=' ')
    cutoff_val = np.fromstring(read_file_chars[16], sep=' ')

    #Determines if system is bimetallic by comparing the p values

    if p_val[0]==p_val[1]:
        sys_bim = False
    else:
        sys_bim = True

    arete = [atom_rad_val[0]*np.sqrt(8)]
    if sys_bim:
        arete.append(atom_rad_val[1]*np.sqrt(8))
        arete.append((arete[0]+arete[1])/2)

    # Unit conversions to arete
    nn = arete/np.sqrt(2)
    dik0 = atom_rad_val[0]+atom_rad_val[1] # Minimal distance between atoms. Sum of atomic radii
    dist = [1/np.sqrt(2)]
    if sys_bim:
        dist.append(nn[1]/arete[0])
        dist.append(nn[2]/arete[0])


    # Converts cutoffs    
    cutoff_start = cutoff_val[0]#/arete[0]
    cutoff_end = cutoff_val[1]#/arete[0]


    # Analytical continuation of potential

    x3 = [0.0, 0.0, 0.0]
    x4 = [0.0, 0.0, 0.0]
    x5 = [0.0, 0.0, 0.0]

    a3 = [0.0, 0.0, 0.0]
    a4 = [0.0, 0.0, 0.0]
    a5 = [0.0, 0.0, 0.0]


    for i in range(3):
        # Old dik0, maybe to correct with units    dik0 = dist[min(i, len(dist)-1)]
        #dik0 = dist[min(i, len(dist)-1)]

        ar = -a_val[i]*np.exp(-p_val[i]*((cutoff_start/dik0)-1))/((cutoff_end-cutoff_start)**3)
        br = -(p_val[i]/dik0)*a_val[i]*np.exp(-p_val[i]*((cutoff_start/dik0)-1))/((cutoff_end-cutoff_start)**2)
        cr = -((p_val[i]/dik0)**2)*a_val[i]*np.exp(-p_val[i]*((cutoff_start/dik0)-1))/((cutoff_end-cutoff_start))

        ab = -qsi_val[i]*np.exp(-q_val[i]*((cutoff_start/dik0)-1))/((cutoff_end-cutoff_start)**3)
        bb = -(q_val[i]/dik0)*qsi_val[i]*np.exp(-q_val[i]*((cutoff_start/dik0)-1))/((cutoff_end-cutoff_start)**2)
        cb = -(((q_val[i]/dik0)**2)*qsi_val[i]*np.exp(-q_val[i]*((cutoff_start/dik0)-1))/((cutoff_end-cutoff_start)))

        x5[i] = (12*ab-6*bb+cb)/(2*((cutoff_end-cutoff_start)**2))
        x4[i] = (15*ab-7*bb+cb)/(((cutoff_end-cutoff_start)))
        x3[i] = (20*ab-8*bb+cb)/2

        a5[i] = (12*ar-6*br+cr)/(2*((cutoff_end-cutoff_start)**2))
        a4[i] = (15*ar-7*br+cr)/(((cutoff_end-cutoff_start)))
        a3[i] = (20*ar-8*br+cr)/2


    Potential = namedtuple('Potential','AtomTypes P Q A Qsi Cohesion Radius Mass CutoffStart CutoffEnd dik0 x3 x4 x5 a3 a4 a5')
    potential = Potential(atoms_type, p_val, q_val, a_val, qsi_val,cohe_val, atom_rad_val, mass_val, cutoff_val[0], cutoff_val[1], dik0, x3, x4, x5, a3, a4, a5)

    return(potential)

def getPressureTwoAtoms(atom_i, atom_j, potential):
    """
    For two atoms 1 and 2 -- given as arrays of form [Atom, x, y, z, Colour]--
    returns the bonding and repulsion pressures. The potential params are all
    stored in the potential named tuple.
    """
    
    ###### den_i should be set zero only wehn I change i but summing over all j
    
    # 1. Determining interaction type
    if atom_i[0] == potential.AtomTypes[0] and atom_j[0] == potential.AtomTypes[0]:
        interaction_type = 0 # Monometallic interaction between atoms of type 1 (cf Pot file)
    if atom_i[0] == potential.AtomTypes[1] and atom_j[0] == potential.AtomTypes[1]:
        interaction_type = 1 # Monometallic interaction between atoms of type 2 (cf Pot file)
    else:
        interaction_type = 2 # Bimetallic interaction
    
    # 2. Determining interatomic distance
    dist_ij = np.sqrt((atom_j[1]-atom_i[1])**2\
                          +(atom_j[2]-atom_i[2])**2\
                          +(atom_j[3]-atom_i[3])**2)
    
    # 3. Pressure calculations
    if dist_ij<=potential.CutoffStart and dist_ij>0: #Distances in A
        
    # 3.1 Pressure calculation for NN
        espo = (dist_ij/potential.dik0)-1

        pres_repul = potential.P[interaction_type]*potential.A[interaction_type]*\
        (np.exp(-potential.P[interaction_type]*espo))/potential.dik0

        pres_bond =-potential.Q[interaction_type]*(potential.Qsi[interaction_type]**2)*\
        (np.exp(-2*potential.Q[interaction_type]*espo))/potential.dik0

        #summing over all j-atoms contributing to force on i
        denom_i = -potential.Qsi[interaction_type]**2*np.exp(-2*potential.Q[interaction_type]*espo)
        # Units denom: eV**2


    # 3.2 Pressure calculation for AN
    elif dist_ij<=potential.CutoffEnd:
        dist_ij_m = dist_ij - potential.CutoffEnd
        
        pres_repul = (5*potential.a5[interaction_type]*(dist_ij_m**4))+\
        (4*potential.a4[interaction_type]*(dist_ij_m**3))+\
        (3*potential.a3[interaction_type]*dist_ij_m**2)
        
        pres_bond = (potential.x5[interaction_type]*(dist_ij_m**5)+\
                     potential.x4[interaction_type]*(dist_ij_m**4)+\
                     potential.x3[interaction_type]*(dist_ij_m**3))*\
        (5*potential.x5[interaction_type]*(dist_ij_m**4)+\
         4*potential.x4[interaction_type]*(dist_ij_m**3)+\
         3*potential.x3[interaction_type]*(dist_ij_m**2))
        
        denom_i = (potential.x5[interaction_type]*(dist_ij_m**5)+\
                   potential.x4[interaction_type]*(dist_ij_m**4)+\
                   potential.x3[interaction_type]*(dist_ij_m**3))**2
        # Units denom: eV**2
        
    # 3.3 No need to calculate pressure for far neighbours

    elif dist_ij>potential.CutoffEnd:
        pres_repul = 0
        pres_bond = 0
        denom_i = 0
    
    # 4. Calculate final pressure between atoms i and j, accounting for distance    
    return(np.array([pres_repul*dist_ij, pres_bond*dist_ij, denom_i]))


### Loop over all atoms to get pressure for each
def pressureMain():
    """
    Reads .pot and movie files and outputs to location
    the same movie with pressures added.
    """
    # 1. Read Files
    potential = readPotentialFile(PATH_TO_POT)
    movie = readMovieFileXYZ(PATH_TO_MOVIE)
    NATOM = int(movie.Headers[0][2])
    open(PATH_TO_NEW_MOVIE, 'w').close() #Clear old movie pressure file
    
    # 2. Startin the loop over all frames
    for frame_num, current_frame in enumerate(movie.Frames):
        print('Analyzing frame: {}/{}'.format(frame_num+1, len(movie.Frames)))

        # Pressure Calculation for that frame
        atom_pressures = []
        for i in range(len(current_frame)): # Loop over i
            pressure_repul_i = 0.0
            pressure_bond_i = 0.0
            summed_denom_i = 0.0
            for j in range(len(current_frame)): # Loop over j of i
                pressure_repul_i += getPressureTwoAtoms(current_frame[i], current_frame[j], potential)[0]
                pressure_bond_i += getPressureTwoAtoms(current_frame[i], current_frame[j], potential)[1]
                summed_denom_i +=getPressureTwoAtoms(current_frame[i], current_frame[j], potential)[2]
            
            pressure_bond_i = pressure_bond_i/np.sqrt(summed_denom_i)
            
            atom_pressures.append(pressure_repul_i+pressure_bond_i)

        # Check that we have one pressure value for each atom
        if (len(atom_pressures)!= NATOM):
            raise ValueError('CAREFUL: Not all atoms have one pressure value')

        #return(atom_pressures)
        # 3. Output Pressure to xyz file
        with open(PATH_TO_NEW_MOVIE, 'a+') as newmovie: # Mode chosen: append
            
            num_lines = sum(1 for line in open(PATH_TO_NEW_MOVIE))

            if (num_lines==0): # No newline for first line -- bugs Ovito if there is newline at beginning
                newmovie.write(str(NATOM)+'\n')
            else:
                newmovie.write('\n' + str(NATOM)+'\n')
                
            newmovie.write('\t'.join(str(item) for item in movie.Headers[frame_num]))
            
            for atom_index, atom_info in enumerate(current_frame):
                atom_info[-1] = atom_pressures[atom_index] # Adding pressure to tuple
                newmovie.write('\n')
                newmovie.write('  \t'.join(str(item) for item in atom_info))

                
for MOVIE_NAME in ['Pt147_sim7_md', 'Pt309_sim7_md']:#, 'Pt561_Sim5_md', 'Pt561_sim8_all_md', 'Pt1415_sim7_md']:
    PATH_TO_POT = './PtPt.pot'

    PATH_TO_MOVIE = './Movie_Files/{}.xyz'.format(MOVIE_NAME)
    PATH_TO_NEW_MOVIE = './Pressure_Results/{}_WithPres.xyz'.format(MOVIE_NAME)

    pressureMain()

