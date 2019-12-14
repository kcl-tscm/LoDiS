import time
import numpy as np

from matplotlib import pyplot as plt

from asap3 import FullNeighborList
from ase import Atoms
from ase.io import read
from scipy.interpolate import interp1d

r_cut=10

"""
TODO(Armand): change function names to be more clear, most of them are
    too generic to be understandable at one glance.
TODO(Armand): Make the pdf NOT read every frames together!!!
TODO(Armand): Make a simple way to get the pdf of every single frame, so I can
    use it as clustering descriptor.
TODO(Armand): separate the PDF cutoff for each frame???.
    not sure if necessary for the cutoff.
TODO(Armand): spline variable name. Don't understand why spline is wanted,
    or if it makes a difference if a linear interpolation is used instead.
    linear interpolation is used for now, and only affects the graph of the PDF.
TODO(Armand): r_cut should be in the pp.py file instead of inside pdf.py.

"""

def read_trajectory(filename, r_cut):
    """ (Armand) Pair Distance Function (PDF).
        Reads a file, uses an arbitrary cutoff to find all pair distances.

    The function reads the file, and then loops over the frames.
    The cell is set to a 100^3 angstrom cell, the total neighbor list is found.
    then loops over each atom, finding the distances between all its neighbors.
    it finally extends a list with all distances found.

    Args:
        filename (str): string path to the filename of interest (xyz file).
        r_cut (float): arbitrary cutoff used for producing the PDF.
            (range:  around 10 angstrom)

    Returns:
        all_distances(:obj:`list` of :obj:`float`): list of all distances found
            from the FullNeighborList applied to every atom found within the
            xyz file.

    """
    traj = read(filename, index = ':')

### (Claudio) This is not needed here, could be needed if PDFs for
###     element-specific distances are needed in the future.
#
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


def create_function(distances, r_cut):
    """ (Armand) Linear Interpolation Function for plotting only.

    Builds a histogram of the all distances observed,
        making bins from 0 to r_cut, with r_cut*100+1 of them.
    The bins are offset by 0.005. (NO IDEA WHY. NEED TO ASK MATTEO ABOUT IT)
    Its last value is omitted. (I dont understand why, but the interp1d
        breaks without it)
    The PDF values are normalized to a probability distribution.
    After which a linear interpolation is performed using scipy.
    A set of values to put into the interpolated function is made, and returned.

    I should probably add... THIS ISNT A SPLINE. if you want a spline you have
        to use : spline = interp1d(x, y, kind=n), where n is the order of your
        spline, and only odd values accepted except for 2. I still dont
        understand why the function is called spline...

    Args:
        distances (:obj:`list` of :obj:`float`): list of all distances found
            using read_trajectory.
        r_cut (float): arbitrary cutoff used for binning all distances found
            (range:  around 10 angstrom)

    Returns:
        spline (:obj:`function`): linear interpolation of a built
            histogram of the Pair Distances.
        values (:obj:`array` of :obj:`float`): x coords for the spline function.

    """
    y, bin_edges = np.histogram(distances, bins = np.linspace(0, r_cut, r_cut*100+1))
    x = bin_edges[:-1] + 0.005
    y = np.array(y) / sum(np.array(y))
    spline = interp1d(x, y)
    values = np.linspace(1, 5, 100)
    #x = np.asarray(x) #(Armand) no purpose in making x an array, especially considering it already is an array...
    return spline, values


def get_cutoff_distance(distances, r_cut):
    """ (Matteo)
    uses gaussian convolution to find minima.
        avoid loops and use arrays instead of lists.

    (Armand) First Minimum finder function. I haven't found where the gaussian
        convolution is...

    Builds a histogram of the all distances observed,
        making bins from 0 to r_cut, with r_cut*100+1 of them.
    The bins are offset by 0.005. (NO IDEA WHY. NEED TO ASK MATTEO ABOUT IT)
    The PDF values are normalized to a probability distribution.
    The PDF values are then smoothed out by a factor of 60% of their current
        value, and 20% value of both of their neighbour values. The first and
        last value are smoothed by 80% and 20% for their neighbour.
    The first minimum is found by looking at values from the start, looking for
        any value where its 10th neighbours are both greater than itself.
    The first minimum distance is then returned.

    Args:
        distances (:obj:`list` of :obj:`float`): list of all distances found
            using read_trajectory.
        r_cut (float): arbitrary cutoff used for binning all distances found
            (range:  around 10 angstrom).

    Returns:
        r_cut_cn (float): distance found of the first minimum, used as a
            physical cutoff in the code.

    """
    y, bin_edges = np.histogram(distances, bins = np.linspace(0, r_cut, r_cut*100+1))
    x = bin_edges[:-1] + 0.005
    y = np.array(y) / sum(np.array(y))
    #x = np.asarray(x)#(Armand) no purpose in making x an array, especially considering it already is an array...

    """(Claudio) adding a tiny convolution of the function y:

    (Armand) This changes the r_cut_cn by about 0.040, still well within
    the first minimum, using the Au2223-FCC-Dh-Ih_test.xyz as benchmark.
    """
    
    y[1:-1] = y[2:]*0.2 + y[1:-1]*0.6+y[:-2]*0.2
    y[0] = y[0]*0.8 + y[1]*0.2
    y[-1] = y[-1]*0.8 + y[-2]*0.2
    
    """(Claudio) This is not efficient, and is also not correct as for i < 10 it
    will look at values like y[-4] which is actually on the other side of
    the array. Check that case please.

    (Armand) Changed the loop range to be after the first 10 steps.
    this should avoid looping backwards in the array.
    """
    for i in range(10, len(y)-10):
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
