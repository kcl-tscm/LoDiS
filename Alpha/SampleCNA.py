from ovito.io import import_file
from ovito.modifiers import CreateBondsModifier, CommonNeighborAnalysisModifier
from ovito.data import particles, BondsEnumerator

import sys
import numpy as np




def Frame_CNA(frame, R_Cut, Masterkey=None, filename=None):
    
    pipeline = import_file(filename)
    pipeline.modifiers.append(CreateBondsModifier(cutoff = R_Cut))

    pipeline.modifiers.append(CommonNeighborAnalysisModifier(
        mode = CommonNeighborAnalysisModifier.Mode.BondBased))
    data = pipeline.compute(frame)
    # The 'CNA Indices' bond property is a a two-dimensional array
    # containing the three CNA indices computed for each bond in the system.
    cna_indices = data.particles.bonds['CNA Indices']
    # This helper function takes a two-dimensional array and computes the frequency
    # histogram of the data rows using some NumPy magic.
    # It returns two arrays (of same length):
    #    1. The list of unique data rows from the input array
    #    2. The number of occurences of each unique row
    def row_histogram(a):
        ca = np.ascontiguousarray(a).view([('', a.dtype)] * a.shape[1])
        unique, indices, inverse = np.unique(ca, return_index=True, return_inverse=True)
        counts = np.bincount(inverse)
        return (a[indices], counts)
    # Used below for enumerating the bonds of each particle:
    bond_enumerator = BondsEnumerator(data.particles.bonds)
    # Loop over particles and print their CNA indices.
    all_cnas = {}
    for particle_index in range(data.particles.count):

        # Create local list with CNA indices of the bonds of the current particle.
        bond_index_list = list(bond_enumerator.bonds_of_particle(particle_index))
        local_cna_indices = cna_indices[bond_index_list]
        # Count how often each type of CNA triplet occurred.
        unique_triplets, triplet_counts = row_histogram(local_cna_indices)
        # Print list of triplets with their respective counts.
    
        for triplet, count in zip(unique_triplets, triplet_counts):
            try:
                all_cnas[tuple(triplet)] += count
                
            except KeyError:
                all_cnas[tuple(triplet)] = count
                if tuple(triplet) not in Masterkey:
                    Masterkey.append(tuple(triplet))
    return all_cnas, Masterkey


"""
Robert:
    This is here just for testing purposes for a given file.
"""

"""
if __name__ == '__main__':
    
    Masterkey= {}
    pipeline = import_file(filename)
    # Create bonds.
    pipeline.modifiers.append(CreateBondsModifier(cutoff = R_Cut))
    Masterkey = []
    
    Metadata={}
    for i in range(800):
        Metadata[i] = Frame_CNA(i)[0]
    filename="/home/k1899676/Documents/PhD/SampleTraj/movie.xyz"
    R_Cut=3.0
"""
            
            
def Get_Heights(Metadata, Masterkey, Norm = False):
    
    """ Jones
    
    Arguments:
        filename: The string name of your input xyz file
            Normally something like 'movie.xyz'
        
        Metadata: The dictionary containng the time ordered CNA signatures 
        and the number of observed occurances.
        
        
        MasterKey: The output from calling the Master function.
        This is to do pairwise comparrison for creating full 
        distributions without having to know what the craic is.
            
        Norm: Default - False
        Whether or not the user wishes to normalise the distribution of 
        CNA signatures for each frame in order to perform meaningful
        statistical analysis.
        
        
    Returns:
        
        Heights: np.array(Frames/Skip, len(MasterKey)) The array containing 
        the (if desired) normalised
        distribution of CNA signature occurances. 
        
    """
    

    Heights=np.zeros((len(Metadata),len(Masterkey)))
    

    for frame in range(len(Metadata)):
        Temp = Metadata[frame].keys()
        for x in Masterkey:
            if x not in Temp:
        
                Metadata[frame][x] = 0

            Heights[frame][Masterkey.index(x)] = Metadata[frame][x]
            
            if Norm == True:
                Heights[frame] = Heights[frame]/sum(Heights[frame])
            
    return Heights
