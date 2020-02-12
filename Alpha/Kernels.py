import numpy as np
from joblib import Parallel, delayed
def minifunc(Data,Band,i):
    X = Data - i
    A1 = -0.5*Band <= X
    A2 = X <= 0.5*Band
    Temp = np.multiply(A1, A2)
    Temp = Temp/Band
    return np.sum(Temp)/(len(Data)*Band)
class Kernels():
    
    """ Robert
    
    This class contains three flavours of kernel desnity estimators.
    Functionally, these attempt to estimate a distribution function 
    for a data set (at this stage, only able to work in 1D but could
    be extended to N dimensions {in theory}).
    
    Essentially, these take a given data point and assign it a weight
    at each point on a pre-defined grid. This weight is determined by the
    type of function used and the proximity of the data point to the position 
    on the pre-defined grid.
    
    As this code is developed, it is likely that each of these functions will
    take additional arguments to determine the cut-off and the desnity of
    grid points.
    
    Note that given the way in which these functions are created and designed, 
    are already normalised.
    
    Note that additional kernels are likely to be implemented into this module
    in the fullnes of time as and when I get around to developing and testing them.
    However; at the moment, these three below appear to be sufficient in efficiency
    and accuracy for the task at hand. But variety IS the spice of life.
    """
    
    def __init__(self, Space, Density):
        self.Space = Space
        self.Density = Density



    def Gauss(Data,Band):
        
        """ Robert
        
        Data: At this stage, this expects the input of a 1D vector containing
        the data to be iterated over. This is expected to be the output of one
        of the functions from the "Distances.py" module. Although, in theory,
        this could be an arbitrary vector of data.
        
        Band: The bandwidth assigned to this kernel. Increasing the bandwidth results 
        in a smoother looking distribution although runs the risk of missing important
        features of the true distribution. Decreasing it tends towards the behaviour
        of each data point being a dirac-delta peak at its precise location in the limit
        of Band tends to zero.
        This would have the effect of re-paramterising the raw data.
        Good results have come from using band ~ 0.125 for estimating a PDDF. Although,
        one may wish to tune this parameter depending on how they wish to present their
        data.
        
        This particular function assigns the weight according to an underlying
        Gaussian distribution. I.e., the weight that a given data point has is 
        Gaussian distributed about the position on the grid under consideration.
        """
        
        Number = '5'
        
        Space=np.linspace(0,8.0,400);Density=[]
        for i in Space:
            P=0
            for j in range(len(Data)):
                X = (Data[j]-i)/Band
                Gauss = (1/np.sqrt(2*np.pi))*np.exp(-(X**2)/2)
                P+=Gauss
            Density.append(P/(len(Data)*Band))
        return Space, Density


    def minifunc(Data,Band,i):
        X = Data - i
        A1 = -0.5*Band <= X
        A2 = X <= 0.5*Band
        Temp = np.multiply(A1, A2)
        Temp = Temp/Band
        return np.sum(Temp)/(len(Data)*Band)

    def Uniform(Data, Band):
        
        """ Robert
        
        Data: At this stage, this expects the input of a 1D vector containing
        the data to be iterated over. This is expected to be the output of one
        of the functions from the "Distances.py" module. Although, in theory,
        this could be an arbitrary vector of data.
        
        Band: The bandwidth assigned to this kernel. Increasing the bandwidth results 
        in a smoother looking distribution although runs the risk of missing important
        features of the true distribution. Decreasing it tends towards the behaviour
        of each data point being a dirac-delta peak at its precise location in the limit
        of Band tends to zero.
        This would have the effect of re-paramterising the raw data.
        
        This assigns a uniform weight to a data point if it exists within an interval
        centred on the grid point under consideration. The weight and interval width are
        intrinsically linked for the purposes of distribution normalisation.
        
        Fine details of PDDFs (including peak splitting) has been best observed with 
        a bandwidth of ~ 0.25.
        """
        
        
        
        Space=np.linspace(0,6.0,300); Density = []
        for i in Space:
            Density.append(minifunc(Data, Band, i))
        return Space, Density
            
            
            
            
            
    
    def Epan(Data,Band):
        
        
        """ Robert
        
        Data: At this stage, this expects the input of a 1D vector containing
        the data to be iterated over. This is expected to be the output of one
        of the functions from the "Distances.py" module. Although, in theory,
        this could be an arbitrary vector of data.
        
        Band: The bandwidth assigned to this kernel. Increasing the bandwidth results 
        in a smoother looking distribution although runs the risk of missing important
        features of the true distribution. Decreasing it tends towards the behaviour
        of each data point being a dirac-delta peak at its precise location in the limit
        of Band tends to zero.
        This would have the effect of re-paramterising the raw data.
        
        This particular function utilises the Epanechnikov convention for assigning weights
        to each data point. In essence, this creates a small semi-circle of weight around 
        each grid point to weight the surroudning data points by.
        
        Testing has had good results for a bandwidth of 0.25 when analysing PDDFs.
        """
        
        Number = '7'
        
        
        Space=np.linspace(0,8.0,400);Density=[]
        for i in Space:
            P=0
            for j in range(len(Data)):
                X = (Data[j]-i)/Band
                P+=0.75*max(1-X**2,0)
            
            Density.append(P/(len(Data)*Band))
        return Space, Density

    
    
    
def KB_Dist(P,Q):
    
    """ Robert
    Calculates the Kullback-Liebler divergence between two distributions.
    
    P: The "initial" distribution against which one wishes to measure the mutual
    entropy of the distribution
    
    Q:
    
    At the moment, there is no actual provision to protect against zero division errors.
    One possible solution could be to define a local varaible, epsilon, which is added to 
    every point in P and prevents it from being zero at any point. 
    
    Note that these two distributions must have identical dimensions or the script
    will not run. 
    
    A reasonable work-around is to define both from an identical linspace.
    """
    
    
    K=0
    Epsilon=0.000001
    Q+=Epsilon
    P+=Epsilon
    for x in range(len(Q)):
        K-=P[x]*np.log(Q[x]/P[x])
    return K

def JSD(P,Q):

    K=0
    Epsilon=0.000001
    Q+=Epsilon
    P+=Epsilon
    for x in range(len(Q)):
        K-=0.5*(P[x]*np.log(2*Q[x]/(Q[x]+P[x]))+Q[x]*np.log(2*P[x]/(P[x]+Q[x])))
    return np.sqrt(K)


