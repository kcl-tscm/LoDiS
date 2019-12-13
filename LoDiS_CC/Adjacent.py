import numpy as np
from ase.io import read
from Movie_Read import read_trajectory as read
import cna
import wikiquote


if __name__ == '__main__':
    filename='/home/k1899676/Documents/JanusMeltMovie1.xyz'
    Cut=3 #The cut-off distance for an atom to be considered adjacent

#In theory, Cut should be evaluated framewise with respect to CNA



def Adjacency_Matrix():
    
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
    
    

    A,B=read(filename)
    Atoms=[np.column_stack((B[i],A[i])) for i in range(len(B))][0]
    Distances=np.zeros((55,55), dtype=np.float); Temp=np.array(np.delete(Atoms, (0),axis=1), dtype=np.float64)
    for i in range(len(Atoms)):
        for j in range(len(Atoms)):
            Euc=(Temp[i,0]-Temp[j,0])**2+(Temp[i,1]-Temp[j,1])**2+(Temp[i,2]-Temp[j,2])**2
            
            Distances[i][j]=np.sqrt(Euc)
    Adjacent=(Distances<Cut).astype(int) #Evaluate if a pair are within R_Cut of eachother
    np.fill_diagonal(Adjacent,0)
    
    
            
            
    return Distances, Adjacent


print(wikiquote.quotes(wikiquote.random_titles(max_titles=1)[0]))

A,B=Adjacency_Matrix()