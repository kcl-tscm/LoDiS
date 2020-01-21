from ase.io import read
import numpy as np
from scipy.stats import ks_2samp
from asap3.analysis import FullCNA
import Kernels


def get_all_cnas(filename, r_cut):
    """ (Claudio)
    Given a trajectory and a cutoff radius, returns a dictionary, sorted by value,
    with all the cnas that appear in the trajectory as keys and the number
    of times they appear as value.
    """
    
    
    """ Jones
    
    Arguments:
        filename: The string name of your input xyz file
            Normally something like 'movie.xyz'
        
        r_cut: The nearest neighbour spacing as a float
            This could be passed from a higher function
            or simply stated at the start of a run.
            
    Returns:
        sorted_cnas_dict: A dictionary whose keys are all of the identified CNA
        signtures appearing in the trajectory.
        The value associated with each key is the number of occurances.
        
    
    If one wishes to observe the entire distribution of CNA signatures over the 
    full production run. Simply call this function and plot its dictionary output.
    
    
    I personally use this to initially identify all of the signatures so that I may
    then identify frame-wise distributions for the cNA signatures.
    """
    
    all_cnas = {}
    traj=read(filename, index=':')
    for j, atoms in enumerate(traj):
        CNA = FullCNA(atoms, r_cut)
        atoms.set_cell([[100,0,0],[0,100,0],[0,0,100]])
        snapshot_cna = CNA.get_normal_cna()
        for i, atomic_cna in enumerate(snapshot_cna):
            for key in atomic_cna:
                try:
                    all_cnas[key] += atomic_cna[key]
                except KeyError:
                    all_cnas[key] = atomic_cna[key]

    sorted_cnas = sorted(all_cnas.items(), key=lambda kv: -kv[1])
    sorted_cnas_dict = {}
    for t in sorted_cnas:
        sorted_cnas_dict[t[0]] = t[1]

    return sorted_cnas_dict
 


def Master(filename, R_Cut):
    
    """ Jones
    
    Arguments: 
        filename: The string name of your input xyz file
            Normally something like 'movie.xyz'
        
        R_Cut: The nearest neighbour spacing as a float
            This could be passed from a higher function
            or simply stated at the start of a run.
           
    Returns:
        MasterKey:
            The sorted list containing all of the CNA signatures which appear in the 
            xyz file under consideration.
            
    """
        
            
    CNAS=get_all_cnas(filename,R_Cut)
    MasterKey=[]
    for keys in CNAS:
        MasterKey.append(keys)
    CNAS=0
    return MasterKey
    



def get_cnas(filename, R_Cut, j):
    """(Claudio)
    Given a trajectory and a cutoff radius, returns a dictionary, sorted by value,
    with all the cnas that appear in the trajectory as keys and the number
    of times they appear as value.
    """
    
    """ Jones
    
    Arguments:
        filename: The string name of your input xyz file
            Normally something like 'movie.xyz'
        
        R_Cut: The nearest neighbour spacing as a float
            This could be passed from a higher function
            or simply stated at the start of a run.
        
        j: Integer whicha specifies the frame of the trajectory file to be called
        
    Returns:
        (Key, Num) The tuple of the CNA signature alongside the number of times it has
        been observed in this given frame.
        
    In general, this function is to be called by the CNA_Sig_Frame and will proceed to be 
    processed by the follow-up guy.
    
    """
        
        
    
    all_cnas = {}
    traj=read(filename, index=j)
    traj.set_cell([[100,0,0],[0,100,0],[0,0,100]])
    CNA = FullCNA(traj, R_Cut)
    snapshot_cna = CNA.get_normal_cna()
    for i, atomic_cna in enumerate(snapshot_cna):
        for key in atomic_cna:
            try:
                all_cnas[key] += atomic_cna[key]
            except KeyError:
                all_cnas[key] = atomic_cna[key]

    sorted_cnas = sorted(all_cnas.items(), key=lambda kv: -kv[1])
    sorted_cnas_dict = {}
    for t in sorted_cnas:
        sorted_cnas_dict[t[0]] = t[1]

    Key=[]; Num=[]
    for keys in sorted_cnas_dict:
        Key.append(keys); Num.append(sorted_cnas_dict[keys])
    return (Key, Num)


def CNA_Sig_Frame(filename, MasterKey, R_Cut, Frames, Skip):
    
    """ Jones
    
    Arguments:
        filename: The string name of your input xyz file
            Normally something like 'movie.xyz'
        
        MasterKey: The output from calling the Master function.
        This is to do pairwise comparrison for creating full 
        distributions without having to know what the craic is.
            
        
        R_Cut: The nearest neighbour spacing as a float
            This could be passed from a higher function
            or simply stated at the start of a run.
        Frames: (int) The number of frames f your movie  file you wish to 
        consider up to
        
        Skip: (int) How many frames you wish to pass over before recording new data.
        
    Returns:
        
        Heights: np.array(Frames/Skip, len(MasterKey)) The array containing the normalised
        distribution of CNA signature occurances. 
        
    """
    
    FullList=[]
    FullSample=[]
    Heights=[]

    for frame in range(int(Frames/Skip)):
        FullList.append(get_cnas(filename,R_Cut,Skip*frame))
        Temp1=FullList[frame][0]
        for x in MasterKey:
            if x not in Temp1:
        
                FullList[frame][0].append(x)
                FullList[frame][1].append(0)
        Sample=[]
        for j in range(len(MasterKey)):
            Sample.append((FullList[frame][0][j],FullList[frame][1][j]))
        FullSample.append(Sample)
        FullSample[frame].sort()
        A,B=zip(*FullSample[frame])
        Heights.append(B/np.sum(B))
        
    return Heights
        
class Dist_Stats():
    
    """ Jones
    
    This class of functions is a group of statistical techniques
    that I am experimenting with as a means of identifying "significant"
    changes in distributions of random variables such as:
        
        Radial Distribution Function (RDF)
        
        Pair Distance Distribution Function (PDDF)
        
        NOT PDF as that means Probability Distribution Function 
        {Not to be confused with the Probability Mass Function}
        
        CNA signature distribution.
        
    Note that these tools do not a-priori require you to have normalised distributions.
    Where it is necessary that they are (PDDF, CNA Sigs), the functional form written
    ensures that they already are.
    
    Isn't life nice like that? :D
    
    In this realisation of the code, each analysis code is to be called for each time frame.
    See the example script.
    
    """
            
    def __init__(self, PStats, KL, JSD):
        self.PStats = PStats
        self.KL = KL
        self.JSD = JSD
        
        
    def PStat(Dist, frame):
        
        """ Jones
        
        Arguments:
            
            Dist: np.array() The Distribution to be analysed for a single frame.
            
            frame: (int) The Frame number under consideration.
            
        Returns:
            
            PStats: The Kolmogorov Smirnov Test statistic.
            Generally speaking, this being <0.05 is sufficing grounds to
            reject the null hypothesis that two sets of observations are drawn
            from the same distribution
            
            A fun wikiquoutes quote because I was bored and felt like learning while coding...
            
            """
            
            
        PStats = (ks_2samp(Dist[0],Dist[frame])[1]) #Performs KS testing and returns p statistic
        #if frame+1 == len(Dist):
           #print(wikiquote.quotes(wikiquote.random_titles(max_titles=1)[0]))
        return PStats
    
    def Kullback(Dist, frame):

        """ Jones
        
        Arguments:
            
            Dist: np.array() The Distribution to be analysed for a single frame.
            
            frame: (int) The Frame number under consideration.
            
        Returns:
            
            KL: The Kullback Liebler divergence:
                This is also known as the mutual information between two distributions.
                It may loosely (and dangerously) interpreted as the similarity between
                two distributions. 
                
                I care about plotting this as I suspect strong delineations in the growth
                of mutual entropy as the system undergoes a phase transition.
            
            A fun wikiquoutes quote because I was bored and felt like learning while coding...
            
            """
            
        KL = (Kernels.KB_Dist(Dist[0],Dist[frame]))
        #if frame+1 == len(Dist):
            #print(wikiquote.quotes(wikiquote.random_titles(max_titles=1)[0]))
        return KL
        
    def JSD(Dist, frame):
        
        """ Jones
        
        Arguments:
            
            Dist: np.array() The Distribution to be analysed for a single frame.
            
            frame: (int) The Frame number under consideration.
            
        Returns:
            
            J: Jenson-Shannon Distance which is a symmetric form the the KL distance above.
            I do not yet understand fully why this should be a superior function to KL but 
            it's another telling discriptor.
            
            A fun wikiquoutes quote because I was bored and felt like learning while coding...
            
            """
        
        J = (Kernels.JSD(Dist[0],Dist[frame]))
        #if frame+1 == len(Dist):
            #print(wikiquote.quotes(wikiquote.random_titles(max_titles=1)[0]))
        return J