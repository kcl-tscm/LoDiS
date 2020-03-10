import numpy as np
import wikiquote
import scipy.sparse as spa

def Adjacency_Matrix(positions, distances, R_Cut):
    
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

    Adj = []; NumAdj = np.zeros((len(positions), len(positions)), dtype=np.float64)
    Tick = 0
    for i in range(1,len(positions)):
        Adj.append(distances[Tick:Tick+len(positions)-i])
        Tick += (len(positions)-i)
        
    for i in range(len(positions)):
        for j in range(len(positions)-(i+1)):
            NumAdj[i][j+i+1] = Adj[i][j]
            NumAdj[j+i+1][i] = Adj[i][j]
    
    Adjacent=(NumAdj<R_Cut).astype(int) #Evaluate if a pair are within R_Cut of eachother
    np.fill_diagonal(Adjacent,0)
        
    Adjacent = spa.csc_matrix(Adjacent)
                        
    return Adjacent

def get_coordination(Adj):
    Temp = spa.csr_matrix.todense(Adj)
    return [ Temp[i].sum() for i in range(len(Temp)) ]
    

def get_coordination_hetero(Adj, R_Cut):
    Adj = (Adj<R_Cut).astype(int)
    CoordA = [ Adj[i].sum() for i in range(len(Adj)) ]
    CoordB = [ Adj[:,j].sum() for j in range(len(Adj[0])) ]
    return (CoordA, CoordB)
    

def R(AdjT, AdjDeltaT):
    TempT = spa.csr_matrix.todense(AdjT)
    TempDeltaT = spa.csr_matrix.todense(AdjDeltaT)
    Temp = TempT-TempDeltaT
    return [ bool(x) for x in Temp.sum(1) ]

def Collectivity(R):
    return float(sum(R)/len(R))

def Concertedness(H1, H2):
    return abs(H2-H1) 
    
print(wikiquote.quotes(wikiquote.random_titles(max_titles=1)[0]), "\n")
