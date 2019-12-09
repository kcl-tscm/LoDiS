"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""##
#######################################################################
                                                                     ##
This script has the purpose of reading a given movie.xyz file        ##
and piecewise returning a given frame which may then be analysed.    ##
                                                                     ##
#######################################################################
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""##

from asap3 import FullNeighborList
from ase import Atoms
from ase.io import read
import numpy as np
import time




def read_trajectory(filename):
  
  """ Robert
  
  This function depends on the ase package and will extract all of the coordinates
  for every atom in every frame of a given trajectory.
  
  filename: This will typically be the name of the movie.xyz file which is pointed to
  by the "Manual_Input" function below. However; this could be any standard xyz file.
  """
  
  
  
    traj = read(filename, index = ':')
    all_positions = [atoms.get_positions() for atoms in traj]
    all_atoms = [atoms.get_chemical_symbols() for atoms in traj]
    return all_positions, all_atoms



def Manual_Input():
  
  """ Robert
  
  This function returns all of the external information files which may be salient
  for analysing the structure as a means of characterising and classifying it.
  
  At the moment, this is fed an energy.out (Could be named otherwise) file and a movie.xyz
  (could be named otherwise) file.
  This may be extended to include data pertaining to the velocity auto-correlation function
  for the purposes of calculating the vibrational density of states. This has not yet been implemented.
  """
  
  
  
    filepath=input("Please specify the absolute path to the directory containing "
                   " your trajectory and enery files.")

    print("Please state the name of the movie file to be analysed with the appropriate extension.")
    try:
        trajectory=input("If the name of your trajectory is "
                         " already 'movie.xyz', simply type "
                         " 'yes'. Otherwise, enter the full name.")
        if trajectory == 'yes':
            trajectory = 'movie.xyz'
    except NameError:
        print("That file does not exist.")


    print("Please state the name of the energy file to be analysed with the appropriate extension.")
    try:
        energy=input("If the name of your energy is "
                         " already 'energy.out', simply type "
                         " 'yes'. Otherwise, enter the full name.")
        if energy == 'yes':
            energy = 'energy.out'
    except NameError:
        print("That file does not exist.")
    tick = time.time()
    Energy=np.loadtxt(filepath+energy)
    Trajectory,Elements=read_trajectory(filepath+trajectory)
    tock = time.time()
    print("Time to read in files: %.2f [s]" %(tock-tick))
    return Energy, Trajectory, Elements, filepath
