import numpy as np

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



def get_CoM(positions):
    return (np.average(positions[:,0]), np.average(positions[:,1]), np.average(positions[:,2]) )

def get_subspecieslist(specie, elements, positions):
    Temp = np.column_stack((elements,positions))
    Temp = [x for x in Temp if x[0] == specie]
    return np.delete(Temp,0,1)

def CoM_Dist(positions, homo = False, specie = None, elements = None):
    
    if homo == False:
        CoM = get_CoM(positions)
        return [distance(x, CoM) for x in positions]
    elif homo == True:
        Temp = get_subspecieslist(specie, elements, positions)
        CoM = get_CoM(positions)
        return [distance(x, CoM) for x in Temp]
    else:
        raise TypeError("You weren't supposed to do that.")
        
    
                 
    


    

def Euc_Dist(i_frame, positions, homo = False, specie = None, elements = None):
    
    if homo == False:
        Distances=[]
        for i in range(len(positions)-1):
            for j in range(i+1,len(positions)):
                Euc = distance(positions[i],positions[j])

                Distances.append(Euc)
        return Distances
    
    elif homo == True:
        Temp = get_subspecieslist(specie, elements, positions)
        for i in range(len(Temp)-1):
            for j in range(i+1,len(Temp)):
                Euc = distance(Temp[i],Temp[j])

                Distances.append(Euc)
        return Distances
    else:
        raise TypeError("You weren't supposed to do that.")
    

        
def Hetero(i_frame, positions, specie, elements):
        
    """ Robert
    
    Note that no species need to be defined for this function as it is understood that LoDiS
    only has provision for mono/bimetallic systems (for the time being) although this
    function could be further generalised (albeit it a potential cost to computation time).
    """
    
    TempA = get_subspecieslist(specie, elements, positions)
    TempB = np.setdiff1d(TempA, positions)
    return [distance(a,b) for a in TempA for b in TempB]
    


    
def RDF(i_frame, positions, Res, R_Cut):
    
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
    
            
    dr = R_Cut / Res
    Radii = np.linspace(0, R_Cut, Res)
    Volumes=np.zeros(Res)
    G=np.zeros(Res)
    
    for i, atom1 in enumerate(positions):
        for j in range(Res):
            r1 = j * dr #Inner radius for the spherical shell
            r2 = r1 + dr #Outer radius increased by increment dr
            v1 = 4.0 / 3.0 * np.pi * r1**3
            v2 = 4.0 / 3.0 * np.pi * r2**3
            Volumes[j] += v2 - v1 #Volume to consider when evaluating distribution

        for atom2 in positions[i:]:
            Distance = distance(atom1, atom2)
            index = int(Distance / dr)
            if 0 < index < Res:
                G[index] += 2 #Identifies when there is an atom at this distance


    for i, value in enumerate(G):
        G[i] = value / Volumes[i] #Rescaling the distribution with respect to enclosing volume
    
    b = (np.diff(np.sign(np.diff(G))) > 0).nonzero()[0] + 1 # local min
    R=Radii[b[1]]
    return Radii, G, R

