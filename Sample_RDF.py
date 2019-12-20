import numpy as np
import wikiquote
import matplotlib.pyplot as plt
from asap3.analysis.rdf import RadialDistributionFunction
from Movie_Read import read_trajectory 
from ase.io import read
import time

if __name__ == '__main__':

    filename='/home/k1899676/Documents/JanusMeltMovie1.xyz'
    A,B=read_trajectory(filename)
    print(wikiquote.quotes(wikiquote.random_titles(max_titles=1)[0]))

#Above are simply test parameters

def distance(a, b):
    
    """ Robert
    
    A simple distance function which takes arguments of
    
    a, b
        These are expected to be arrays containing three elements
        (x, y, z)
        Being the respective euclidean coordinates for atoms a and b
        at a given point in time.
        
    Reurns a single float being the euclidean distance between the atoms.
    
    """
    
    dx = abs(a[0] - b[0])
     
    dy = abs(a[1] - b[1])
     
    dz = abs(a[2] - b[2])
 
    return np.sqrt(dx**2 + dy**2 + dz**2)



def RDF(Resolution, Trajectory, R_Cut):
    
    """ Robert
    
    Args:
        Resolution: 
            int data type representing how finely you wish to make 
            the grid. Usually set in the order of 100
        
        Trajectory: 
            Single frame of xyz coordinates for a set of atoms
            Is expected to be iterated over and so will only take a single frame of xyz
        
        R_Cut: 
            Float type variable which indicates how far you wish to create
            the distribution for.
            Good practice is to set it to ~0.5 Diameter of the cluster
            Tested with 10 Angstroms

    Returns:
        Radii:
            A numpy array of all the radii the distribution has been computed over
            Will have length of "Resolution" and is to be used as the x axis on
            an RDF plot.
        
        G:
            A numpy array of the (unnormalised) calculated RDF values corresponding 
            to the respective radius in Radii. To be set on the y axis in a given
            RDF plot.
            
    Note bene:
        
        In the future, this function will be generalised to calculate 
            (full, homo, hetero)
        RDF plots. 
        Given that for the time being, we are mostly concerned with monometallic systems
        This is not a HUGE issue.
    """

    Resolution=100
    R_Cut=10
    dr = R_Cut / Resolution
    Radii = np.linspace(0,R_Cut,Resolution)
    Volumes=np.zeros(Resolution)
    G=np.zeros(Resolution)


    for i, atom1 in enumerate(Trajectory):
        for j in range(Resolution):
            r1 = j * dr #Inner radius for the spherical shell
            r2 = r1 + dr #Outer radius increased by increment dr
            v1 = 4.0 / 3.0 * np.pi * r1**3
            v2 = 4.0 / 3.0 * np.pi * r2**3
            Volumes[j] += v2 - v1 #Volume to consider when evaluating distribution

        for atom2 in Trajectory[i:]:
            Distance = distance(atom1, atom2)
            index = int(Distance / dr)
            if 0 < index < Resolution:
                G[index] += 2 #Identifies when there is an atom at this distance


    for i, value in enumerate(G):
        G[i] = value / Volumes[i] #Rescaling the distribution with respect to enclosing volume

    return Radii, G
if __name__ == '__main__':
    tick=time.time()
    Radii, G = RDF(100, A[0], 10.0)
    tock=time.time()
    print(tock-tick)

