# -*- coding: utf-8 -*-
"""
Claudio's Radial Distribution Function calcualtion averaged over all time frames.
r_cut_cn is the first big minimum of the plot, whihc is used as a cut off for the CN and GCN calculation. """
​
from asap3 import FullNeighborList
from ase import Atoms
from ase.io import read
import numpy as np
from matplotlib import pyplot as plt
import time
from scipy.interpolate import interp1d
import pandas as pd
​
r_cut=10
​
def read_trajectory(filename, r_cut):
    traj = read(filename, index = ':')
    atom_number_list = [atoms.get_atomic_numbers() for atoms in traj]
    flat_atom_number = np.concatenate(atom_number_list)
    elements = np.unique(flat_atom_number, return_counts=False)
    all_distances = []
    for i, atoms in enumerate(traj):
        atoms.set_cell([[100, 0, 0], [0, 100, 0], [0, 0, 100]])
        nl = FullNeighborList(r_cut, atoms=atoms)
        for i in np.arange(len(atoms)):
            indices, positions, distances = nl.get_neighbors(i)
            all_distances.extend(distances**0.5)
    return all_distances
​
def create_function(distances, r_cut):
    y, bin_edges = np.histogram(distances, bins = np.linspace(0, r_cut, r_cut*100+1))
    x = bin_edges[:-1] + 0.005
    y = np.array(y) / sum(np.array(y))
    spline = interp1d(x, y)
    values = np.linspace(1, 5, 100)
    x=np.asarray(x)
    minima=[]
    """ Finding the minima in a very general way, finding the 10th y-value (percentage of atoms within a certain radius) 
    value before and after a point on the x-axis """
    
    for i in range (len(y)-10):
        if y[i]<y[i+10] and y[i]<y[i-10]:               
            minima.append(x[i])
    global r_cut_cn
    """r_cut_cn is the first encountered minima, used for cn calculation"""
    r_cut_cn=minima[0]                        
    print("Radius for cutoff=", r_cut_cn)
    plt.plot(values, spline(values))
    plt.show()
    return spline, r_cut_cn  
​
if __name__ == '__main__':
    tic = time.time()
    distances  = read_trajectory("cu561/freezing_short_1200-400/i1/movie.xyz", r_cut)
    toc = time.time()
    print("Time to calculate distances: %.2f [s]" %(toc-tic))
​
    tic = time.time()
    create_function(distances, r_cut)
    toc = time.time()
    print("Time to create spline: %.2f [s]" %(toc-tic))
​
print(r_cut_cn)