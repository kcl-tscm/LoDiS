import numpy as np

import time



class Distance():

    def __init__(self, distances, species):
        
        self.distances = distances
        
        self.species = species
""" Robert
This class of functions is designed to extract information regarding atomic positions
for a whole trajectory and generate data regarding spatial distributions of all, homo, and
hetero combinations.
"""
        
    def Euc_Dist(Vector):
        
        """ Robert
        Vector: The tensor containing the xyz coordinates of all atoms in a given trajectory.
        This function will only act on one frame at any given time but can be iteratively 
        called so at to generate information regarding specific frames of a simulation.
        
        This is generally expected to be the output of "read_trajectory" in the 
        "Movie_Read.py" module which generates the xyz coordinates as well as the 
        numbered list of elements. These are initially kept separate although there
        does exist provision to stack these tensors to facilitate the homo/hereo 
        PDDFs.
        
        This function is indiscriminate regarding chemical species. It will calculate 
        pair-wise distances regardless of chemical species and return a 1D numpy array
        of all possible pairwise distances which is N!(N-1)!/2
        So bear in mind that this function scales poorly with system size!
        """
        
        
        Distances=[]; Temp=np.array(np.delete(Vector, (0),axis=1), dtype=np.float64)
        for i in range(len(Vector)-1):
            for j in range(i+1,len(Vector)):
                Euc=(Temp[i,0]-Temp[j,0])**2+(Temp[i,1]-Temp[j,1])**2+(Temp[i,2]-Temp[j,2])**2

                Distances.append(np.sqrt(Euc))
        return Distances




    
    def M1_M1(Vector,Species):
        
        """ Robert
        
        Vector: The tensor containing the xyz coordinates of all atoms in a given trajectory.
        This function will only act on one frame at any given time but can be iteratively 
        called so at to generate information regarding specific frames of a simulation.
        
        Species: The chemical symbol for the element to be considered. Parameter expected 
        to be input as such
        
        "M1_M1(Positions, Au)" if one were to analyse the xyz trajectory and was interested only 
        in gold.
        
        This function has near identical architecture to the one defined above, "Euc_Dist",
        with the exception that the function actively searches for similar chemical species
        from the input trajectory.
        """
        
        tick=time.time()
        Distances=[]; Temp=np.array(np.delete(Vector, (0),axis=1), dtype=np.float64)
        for i in range(len(Vector)-1):
            for j in range(i+1,len(Vector)):
                if Vector[i,0]==Species and Vector[j,0]==Species:
                    Euc=(Temp[i,0]-Temp[j,0])**2+(Temp[i,1]-Temp[j,1])**2+(Temp[i,2]-Temp[j,2])**2
                    Distances.append(np.sqrt(Euc))
        tock=time.time()
        print("Time to calculate distances: %.2f [s]" %(toc-tic))
        return Distances

    
    def M1_M2(Vector):
        
        """ Robert
        
        Vector: The tensor containing the xyz coordinates of all atoms in a given trajectory.
        This function will only act on one frame at any given time but can be iteratively 
        called so at to generate information regarding specific frames of a simulation.
        
        Note that no species need to be defined for this function as it is understood that LoDiS
        only has provision for mono/bimetallic systems (for the time being) although this
        function could be further generalised (albeit it a potential cost to computation time).
        """
        
        tick=time.time()
        Distances=[]; Temp=np.array(np.delete(Vector, (0),axis=1), dtype=np.float64)
        for i in range(len(Vector)-1):
            for j in range(i+1,len(Vector)):
                if Vector[i,0]!=Vector[j,0]:
                    Euc=(Temp[i,0]-Temp[j,0])**2+(Temp[i,1]-Temp[j,1])**2+(Temp[i,2]-Temp[j,2])**2
                    Distances.append(np.sqrt(Euc))
        tock=time.time()
        print("Time to calculate distances: %.2f [s]" %(tock-tick))
        return Distances
                
                

if __name__ == '__main__':
    import Movie_Read as Read

    Energy, Trajectory, Elements=Read.Manual_Input()
    Positions=[];Distances=[]

    for x in range(11000):
        Positions.append(np.column_stack((Elements[x],Trajectory[x])))
        Distances.append(Distance.Euc_Dist(Positions[x]))   #All possible pairings

