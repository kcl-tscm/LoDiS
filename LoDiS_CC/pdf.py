from asap3 import FullNeighborList
from ase import Atoms
from ase.io import read
import numpy as np
from matplotlib import pyplot as plt
import time
from scipy.interpolate import interp1d

r_cut=10
""" (Armand)
Pdf function. Given a filename and a cutoff,
returns a list with all distances within the cutoff found appended.

It first finds the elements composition, haven't found if it does something with it for now
It loops over the frames after finding elements
and it finds the total neighbor list before getting the distances between neighbors
returns a list of the distances found between neighbors closer than 10A.
"""
def read_trajectory(filename, r_cut):
    traj = read(filename, index = ':')

###    This is not needed here, could be needed if PDFs for element-specific distances are needed in the future.

#     atom_number_list = [atoms.get_atomic_numbers() for atoms in traj]
#     flat_atom_number = np.concatenate(atom_number_list)
#     elements = np.unique(flat_atom_number, return_counts=False)
#     #checks the elemental composition, displays them in elements#

    all_distances = []
    for i, atoms in enumerate(traj): #loops over frames
        atoms.set_cell([[100, 0, 0], [0, 100, 0], [0, 0, 100]])
        nl = FullNeighborList(r_cut, atoms=atoms)#to note, the cut_off is HUGE here. 10A.
        for i in np.arange(len(atoms)):
            indices, positions, distances = nl.get_neighbors(i) #grabs the data from FullNeighborList
            all_distances.extend(np.sqrt(distances)) #extends the list with the new neighbor distances obtained.
    return all_distances

""" (Armand)
Splining function. Given the list of all distances,
bins all the distances and splines them to find a pdf function.
returns an array of the histogrammed values and the spline function
Useful for histogram presentation and to find the first minima later.
"""
def create_function(distances, r_cut):
    y, bin_edges = np.histogram(distances, bins = np.linspace(0, r_cut, r_cut*100+1))
    x = bin_edges[:-1] + 0.005
    y = np.array(y) / sum(np.array(y))
    spline = interp1d(x, y)
    values = np.linspace(1, 5, 100)
    x = np.asarray(x) #(Armand) not sure what's the point of this line
    return spline, values


""" (Matteo)
uses gaussian convolution to find minima. avoid loops and use arrays instead of lists.
"""
def get_cutoff_distance(distances, r_cut):
    y, bin_edges = np.histogram(distances, bins = np.linspace(0, r_cut, r_cut*100+1))
    x = bin_edges[:-1] + 0.005
    y = np.array(y) / sum(np.array(y))
    x = np.asarray(x)

    # (Claudio) adding a 
    # Tiny convolution of the the function y:
    y[1:-1] = y[2:]*0.2 + y[1:-1]*0.6+y[:-2]*0.2
    y[0] = y[0]*0.8 + y[1]*0.2
    y[-1] = y[-1]*0.8 + y[-2]*0.2
    
    # (Claudio) This is not efficient, and is also not correct as for i < 10 it
    # will look at values like y[-4] which is actually on the other side of the array.
    # Check that case please.
    for i in range(len(y)-10):
        if y[i] < y[i+10] and y[i] < y[i-10]:
            r_cut_cn = x[i]
            break # (Claudio) This stops the cycle so that you don't waste time looking for more minima
           
    return r_cut_cn

"""(Matteo)
recognise particular boundaries (422 planes).
Faster diffusion in proximity of which plane?
monitor dynamical behaviour of this planes. Pressure calculations?
get the versor of the planes.
"""
