import numpy as np

import time



class Distance():

    def __init__(self, distances, species):
        
        self.distances = distances
        
        self.species = species


    def Euc_Dist(Vector):
        Distances=[]; Temp=np.array(np.delete(Vector, (0),axis=1), dtype=np.float64)
        for i in range(len(Vector)-1):
            for j in range(i+1,len(Vector)):
                Euc=(Temp[i,0]-Temp[j,0])**2+(Temp[i,1]-Temp[j,1])**2+(Temp[i,2]-Temp[j,2])**2

                Distances.append(np.sqrt(Euc))
        return Distances





    def M1_M1(Vector,Species):
        
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

