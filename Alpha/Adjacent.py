
import numpy as np
from ase.io import read
import wikiquote
import scipy.sparse as spa


if __name__ == '__main__':
    filename='/home/k1899676/Documents/JanusMeltMovie1.xyz'
    Cut=3 #The cut-off distance for an atom to be considered adjacent

#In theory, Cut should be evaluated framewise with respect to CNA



def Adjacency_Matrix(positions, R_Cut):
    
    """ Robert
        Args:
            Not yet fully implemented but in theory will take arguments of a single frame
            of an xyz trajectory to then be iterated over and the cut-off imposed to be the
            nearest nieghbour distances
            
        Returns:
            Distances:
                N x N array (N being number of atoms present) containing all pairwise distances 
                between atoms
                I.e., Entry[i][j] is the absolute distance between atoms i and j
            
            Adjacent: N x N array of 0s and 1s identifying the "Truth" of a given atom
            pairing being neighbours based on the criterion of R_Cut
            I.e., Dist<R_Cut returns a 1 for that entry as the condition has been met
            All diagonal elements are set to 0 as it is meaningless to be your own neighbour.
    """
    
    


    Distances=np.zeros((len(positions),len(positions)))
    for i in range(len(positions)):
        for j in range(len(positions)):
            Euc=(positions[i,0]-positions[j,0])**2+(positions[i,1]-positions[j,1])**2+(positions[i,2]-positions[j,2])**2
            
            Distances[i][j]=np.sqrt(Euc)
    Adjacent=(Distances<R_Cut).astype(int) #Evaluate if a pair are within R_Cut of eachother
    np.fill_diagonal(Adjacent,0)
    
    
    Adjacent = spa.csc_matrix(Adjacent)
            
            
    return Adjacent


#print(wikiquote.quotes(wikiquote.random_titles(max_titles=1)[0]), "\n")

#A,B=Adjacency_Matrix(filename, 0, 3.0)
