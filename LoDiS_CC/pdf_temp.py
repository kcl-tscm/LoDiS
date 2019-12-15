import time
import numpy as np

from matplotlib import pyplot as plt

from asap3 import FullNeighborList
from ase import Atoms
from ase.io import read
from scipy.interpolate import interp1d

"""
#TODO(Armand):remove r_cut from this file, make it something fed by the main
    code file. Random values shouldn't be in module files.
    (Armand): removed it. put in pp.py for now to test with it, should be found
        in the main file anyway.
#TODO(Armand): change function names to be more clear, most of them are
    too generic to be understandable at one glance.
    (Armand): changed both first functions to names that describe their role
        clearer.
#TODO(Armand): Make the pdf NOT read every frames together!!!
    (Armand): God this took me so long. Made an array of arrays of each frame
        distances, due to different neighbour numbers... so couldn't make a
        single n by neighbour number array.
#TODO(Armand): Make a simple way to get the pdf of every single frame, so I can
    use it as clustering descriptor.
    (Armand): Done. Through using a larger r_cut, the whole pdf of every frame
        can be obtained. Might need to make a new def that only gets the pdf of
        a specific frame though. Would save a LOT of memory for the PDF array.
#TODO(Armand): separate the PDF cutoff for each frame???.
    not sure if necessary for the cutoff.
    (Armand): seperated the pdf of every frame to give the physical cutoff of
    every single frame inside of a list, and returns an array copy of that list
    to save memory. Adapting CN next to obtain the CNs with each frame's cutoff.
#TODO(Armand): spline variable name. Don't understand why spline is wanted,
    or if it makes a difference if a linear interpolation is used instead.
    linear interpolation is used for now, and only affects the graph of the PDF.
    (Armand): changed the variable name for now to pdf_interpolated,still need
    to ask if a spline is specifically wanted.
#TODO(Armand): Once every TODO is done, update the comments.
    (Armand): Update the args and returns, working on clearer function
        explanations now.
TODO(Armand): If necessary, code a PDF function that only finds the pdf of a
    specific frame, so that a huge r_cut can be used without waiting ages on
    movie files with a lot of frames. It isnt an issue when r_cut=10, but it
    would become one if a pdf descriptor requires r_cut=40, with 100 frames...
"""

def arbitrary_pair_distribution_function(filename, r_cut):
    """ (Armand) Pair Distance Function (PDF). Reads a file,
        uses an arbitrary cutoff to find all pair distances.

    The function reads the file, and then loops over the frames.
    The cell is set to a 100^3 angstrom cell, the total neighbor list is found.
    then loops over each atom, finding the distances between all its neighbors.
    To store the individual frame data, an array of arrays is made.
    Since we dont know how big the final array gets, to save memory,
    we first make a list which stores the distances found within one unique
    frame. After this, the list is transformed into an array and stored in a
    list, appended with each of the distances found for each frame.
    Once all the frames are looped over we finally transform this list into
    an array, giving us an array of arrays, of different length. To access it
    requires to use array[frame number][distance values]. Not simple to come up
    with, but fixes the problem of multiple frames with different neighbour
    numbers. This fix is applicable to almost every code used so far.

    Args:
        filename (str): string path to the filename of interest (xyz file).
        r_cut (float): arbitrary cutoff used for producing the PDF.
            (range:  around 10 angstrom)

    Returns:
        all_distances_array(:obj:`array` of :obj:`array` of float): array of
            arrays, where each array corresponds to their frame, and contains
            the distances found from the FullNeighborList applied to every atom
            of the frame.

    """
    traj = read(filename, index = ':')

### (Claudio) This is not needed here, could be needed if PDFs for
###     element-specific distances are needed in the future.
#
#     atom_number_list = [atoms.get_atomic_numbers() for atoms in traj]
#     flat_atom_number = np.concatenate(atom_number_list)
#     elements = np.unique(flat_atom_number, return_counts=False)
#     #checks the elemental composition, displays them in elements#

    all_distances_list=[] #list which stores the arrays of distances found.
    for i, atoms in enumerate(traj):
        all_distances_of_frame = [] #list which stores the distances of a single frame.
        atoms.set_cell([[100, 0, 0], [0, 100, 0], [0, 0, 100]])
        nl = FullNeighborList(r_cut, atoms=atoms)
        for i in np.arange(len(atoms)):
            indices, positions, distances = nl.get_neighbors(i)
            all_distances_of_frame.extend(np.sqrt(distances)) #filling the frame list
        all_distances_list.append(np.asarray(all_distances_of_frame))
        #line above: converting the frame distances list to an array,
            #and filling a list with them.
    all_distances_array=np.asarray(all_distances_list)
    #line above: finally, converting it all to an array of arrays of different length.
    return all_distances_array


def pdf_linear_interpolation(all_distances_array, r_cut):
    """ (Armand) Linear Interpolation Function for plotting only.

    Builds a histogram of the all distances observed, in ALL frames, done by
        flattening all_distances_array, and making bins from 0 to r_cut,
        with r_cut*100+1 of them.
    The bins are offset by 0.005. (NO IDEA WHY. NEED TO ASK MATTEO ABOUT IT)
    Its last value is omitted. (I dont understand why, but the interp1d
        breaks without it)
    The PDF values are normalized to a probability distribution.
    After which a linear interpolation is performed using scipy.
    A set of values to put into the interpolated function is made, and returned.

    I should probably add... THIS ISNT A SPLINE. if you want a spline you have
        to use : pdf_interpolated = interp1d(x, y, kind=n), where n is the order
        of your spline, and only odd values accepted except for 2.

    Args:
        all_distances_array(:obj:`array` of :obj:`array` of float):  array of
            arrays, where each array corresponds to their frame, and contains
            the distances found from the FullNeighborList applied to every atom
            of the frame.
        r_cut (float): arbitrary cutoff used for binning all distances found
            (range:  around 10 angstrom)

    Returns:
        pdf_interpolated (:obj:`function`): linear interpolation of a built
            histogram of the Pair Distances.
        values (:obj:`array` of float): x coords for pdf_interpolated.

    """
    flattened_distances=np.concatenate((all_distances_array[:][:]))
    y, bin_edges = np.histogram(flattened_distances, bins = np.linspace(0, r_cut, r_cut*100+1))
    x = bin_edges[:-1] + 0.005
    y = np.array(y) / sum(np.array(y))

    pdf_interpolated = interp1d(x, y)
    x_values = np.linspace(1, 5, 100)
    #x = np.asarray(x) #(Armand) not sure what's the point of this line
    return pdf_interpolated, x_values


def get_cutoff_distance(all_distances_array, r_cut):
    """(Armand) First Minimum finder function.

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

    (MAKE THIS COMMENT BETTER)
    Currently, this function is getting the first minimum of the pdf of ALL
    of the frames together. Im currently changing it so that r_cut_physical
    (r_cut_cn as it was before) becomes a list of the values of the first
    minimum for every single frame using the array of distances arrays. done.

    Args:
        all_distances_array(:obj:`array` of :obj:`array` of float):  array of
            arrays, where each array corresponds to their frame, and contains
            the distances found from the FullNeighborList applied to every atom
            of the frame.
        r_cut (float): arbitrary cutoff used for binning all distances found
            (range:  around 10 angstrom).

    Returns:
        r_cut_physical (:obj:`array` of float): array of the distances found of
        the first minimum, used as a physical cutoff in the code, where the
        first value corresponds to the first frame, etc.

    """
    r_cut_physical=[]
    for i in range(len(all_distances_array)):
        y, bin_edges = np.histogram(all_distances_array[i][:], bins = np.linspace(0, r_cut, r_cut*100+1))
        x = bin_edges[:-1] + 0.005
        y = np.array(y) / sum(np.array(y))
        #x = np.asarray(x) #(Armand) no purpose in making x an array, especially considering it already is an array...
        """(Claudio) adding a tiny convolution of the function y:

        (Armand) This changes the r_cut_cn by about 0.040, still well within
        the first minimum, using the Au2223-FCC-Dh-Ih_test.xyz as benchmark.
        """
        y[1:-1] = y[2:]*0.2 + y[1:-1]*0.6+y[:-2]*0.2
        y[0] = y[0]*0.8 + y[1]*0.2
        y[-1] = y[-1]*0.8 + y[-2]*0.2

        """(Claudio) This is not efficient, and is also not correct as for
        i < 10 it will look at values like y[-4] which is actually on the other
        side of the array. Check that case please.

        (Armand) Changed the loop range to be after the first 10 steps.
        this should avoid looping backwards in the array.
        """
        for j in range(10, len(y)-10):
            if y[j] < y[j+10] and y[j] < y[j-10]:
                r_cut_physical.append(x[j])
                break # (Claudio) This stops the cycle so that you don't waste time looking for more minima

    return np.asarray(r_cut_physical)
