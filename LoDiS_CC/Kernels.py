import numpy as np
class Kernels():
    
    def __init__(self, Space, Density):
        self.Space = Space
        self.Density = Density



    def KDE_Gauss(Data,Band):
        Space=np.linspace(0,4.0,500);Density=[]
        for i in Space:
            P=0
            for j in range(len(Data)):
                X = (Data[j]-i)/Band
                Gauss = (1/np.sqrt(2*np.pi))*np.exp(-(X**2)/2)
                P+=Gauss
            Density.append(P/(len(Data)*Band))
        return Space, Density


    def KDE_Uniform(Data,Band):
        Space=np.linspace(0,4.0,500);Density=[]
        for i in Space:
            P=0
            for j in range(len(Data)):
                X = (Data[j]-i)/Band
                if -0.5*Band <= X <= 0.5*Band:
                    P+=1.0/Band
                else:
                    P+=0
            
            Density.append(P/(len(Data)*Band))
        return Space, Density
            
    
    def KDE_Epan(Data,Band):
        Space=np.linspace(0,4.0,500);Density=[]
        for i in Space:
            P=0
            for j in range(len(Data)):
                X = (Data[j]-i)/Band
                P+=0.75*max(1-X**2,0)
            
            Density.append(P/(len(Data)*Band))
        return Space, Density

    
    
    
def KB_Dist(P,Q):
    K=0
    for x in range(len(Q)):
        K-=P[x]*np.log(Q[x]/P[x])
    return K




if __name__ == '__main__':
    import Movie_Read as Read
    import Distances as Dist

    Energy, Trajectory, Elements=Read.Manual_Input()
    Positions=[];Distances=[];PDDF=[]

    for x in range(10):
        Positions.append(np.column_stack((Elements[x],Trajectory[x])))
        Distances.append(Dist.Distance.Euc_Dist(Positions[x]))   #All possible pairings
        PDDF.append(Kernels.KDE_Uniform(Distances[x],0.25))