from Adjacent import *
from CNA  import *
from Kernels import *
from Pressure import *
from RDF import *
from Energy import *
from ase.io import read

from DistFuncs import *

metadata={}
result_cache={}


"""
In the full version, the user will either be called for these arguments or may simply submit a script.
I'd like for both versions to be effective. But, for now, it shall remain hard-coded to faciliate 
debugging and support transparency.
"""


Start: int = 0
End: int = 50
Step: int = 10

Time=int((End-Start)/Step)

System = {
        'base_dir' : '/home/k1899676/Documents/PhD/Coding/Samples/',
        'movie_file_name' : 'JanusMeltMovie1.xyz',
        'energy_file_name' : 'JanusMeltEn1.out'}



settings: dict = None 
filename = System['base_dir']+System['movie_file_name']
Dataset = read(System['base_dir']+System['movie_file_name'], index = 0)
all_atoms = Dataset.get_chemical_symbols()
Masterkey=Master(System['base_dir']+System['movie_file_name'],3.0)


"""

Robert:
    
    I know that this quantities dictionary is a mess, the general idea is that this will
    be filled up with information from the class implementation 'automatically'. Because
    this has been thrown together in an afternoon with all of the calculators rejigged for
    ease of calling as opposed to efficiency, it looks ugly as sin! Sorry, team.
    
"""

Quantities = {
        'euc' : np.zeros((Time,int(len(all_atoms)*(len(all_atoms)-1)/2))), 'rdf' : np.zeros((Time,2,500)), 
        'pos' : np.zeros((Time,len(all_atoms),3)), 'cna' : np.zeros((Time,len(Masterkey))), 
        'adj' : np.zeros((Time, 2, len(all_atoms), len(all_atoms))), 'pdf' : np.zeros((Time, 2, 400)),
        'pdfhomo' : np.zeros((Time, 2, 400))}

for i in range(End-Start):
    First=True
    
    Dataset = read(System['base_dir']+System['movie_file_name'], index = i)
    all_positions = Dataset.get_positions()
    all_atoms = Dataset.get_chemical_symbols()
    
    used=set()
    Species = [x for x in all_atoms if x not in used and (used.add(x) or True)]
    
    NAtoms = len(all_atoms)
    
    
    metadata['Elements'] = all_atoms
    metadata['Species'] = Species
    metadata['NSpecies']=len(Species)
    metadata['NFrames'] = Time
    metadata['NAtoms'] = NAtoms
    
    
    """
    Robert:
        This little doozy under here just sets up the dictionary to add your observations to.
    """
    
    if i==Start:
        for key in Quantities:
            metadata[key] = Quantities[key]
    

    result_cache['pos'] = all_positions
    metadata['pos'][int(i/Step)] = result_cache['pos']

    result_cache['euc'] = Euc_Dist(i, result_cache['pos'])
    metadata['euc'] = result_cache['euc']
    
    result_cache['pdf'] = Kernels.Uniform(result_cache['euc'], 0.25)
    metadata['pdf'] = result_cache['pdf']
    
    
    
    """
    Robert:
        For the homo, it's a bit trickier for quick introduction stuff.
    """
    
    
    
    result_cache['homo'] = Homo(i, result_cache['pos'], metadata['Elements'], 'Au')
    metadata['homo'] = np.zeros((Time, len(result_cache['homo'])))
    metadata['homo'][int(i/Step)] = result_cache['homo']
    
    result_cache['pdfhomo'] = Kernels.Uniform(result_cache['homo'], 0.25)
    metadata['pdfhomo'] = result_cache['pdfhomo']
   
    
    result_cache['rdf'] = RDF(i, result_cache['pos'], 500, 10.0)
    metadata['rdf'][int(i/Step)] = result_cache['rdf']   
    
    result_cache['adj'] = Adjacency_Matrix(i, result_cache['pos'], R_Cut=3.0)
    metadata['adj'][int(i/Step)] = result_cache['adj']

    
    i+=Step*1


"""
Robert: 
    The way in which I've written this function makes it a pain to call at each frame.
    So instead, it can be called at the end, for now and then spruced up for the big release!
"""


metadata['cna'] = CNA_Sig_Frame(filename, Masterkey, R_Cut=3.0, Frames = End, Skip = Step)