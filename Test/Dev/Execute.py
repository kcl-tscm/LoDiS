from Adjacent import *
from CNA  import *
from Kernels import *
from Energy import *
from ase.io import read
import time
from DistFuncs import *
import pickle
metadata={}
result_cache={}


"""
In the full version, the user will either be called for these arguments or may simply submit a script.
I'd like for both versions to be effective. But, for now, it shall remain hard-coded to faciliate 
debugging and support transparency.




To accelerate the whole prrocess, we shall only be evaluating the PDF & Adjacency matrix every 200 frames
of the given step increment. The purpose of which is to identify means by which we may bring the 
analysis down to an approximate scaling of 0.01s/AtomFrame
"""


Start: int = 0
End: int = 800
Step: int = 1

Skip=100

"""
Robert:
    Change the arguments in the System dictionary to 
    reflect the directory and files you are working with.
"""

System = {
        'base_dir' : '/home/k1899676/Documents/PhD/SampleTraj/',
        'movie_file_name' : 'movie.xyz',
        'energy_file_name' : 'JanusMeltEn1.out',
        'Homo' : 'Cu',
        'Start' : 0, 'End' : 800, 'Step' : 1, 'Skip' : 100,
        'PdfStats' : True, 'HomoStats' : True, 'RdfStats' : True, 'CnaStats' : True}

Time=int((End-Start)/Step)

settings: dict = None 
filename = System['base_dir']+System['movie_file_name']
Dataset = read(filename, index = 0)
all_atoms = Dataset.get_chemical_symbols()
Masterkey=Master(filename,3.0)



"""

Robert:
    
    I know that this quantities dictionary is a mess, the general idea is that this will
    be filled up with information from the class implementation 'automatically'. Because
    this has been thrown together in an afternoon with all of the calculators rejigged for
    ease of calling as opposed to efficiency, it looks ugly as sin! Sorry, team.
    
"""

Quantities = {
        #'euc' : np.zeros((Time,int(len(all_atoms)*(len(all_atoms)-1)/2))),
        'rdf' : np.zeros((Time,2,500)), 
        'pos' : np.zeros((Time,len(all_atoms),3)), 'cna' : np.zeros((Time,len(Masterkey))), 
        'adj' : np.zeros((int(Time/Skip), len(all_atoms), len(all_atoms))), 'pdf' : np.zeros((int(Time/Skip), 2, 300))
        #'pdfhomo' : np.zeros((Time, 2, 300))
        }
FirstPDF = True
FirstRDF = True
First = True

for i in range(Start, End, Step):
    T0=time.time()
    
    
    Dataset = read(filename, index = i)
    all_positions = Dataset.get_positions()
    all_atoms = Dataset.get_chemical_symbols()
    
    used=set()
    Species = [x for x in all_atoms if x not in used and (used.add(x) or True)]
    
    NAtoms = len(all_atoms)
    
    if i == Start:
        metadata['Elements'] = all_atoms
        metadata['Species'] = Species
        metadata['NSpecies']=len(Species)
        metadata['NFrames'] = Time
        metadata['NAtoms'] = NAtoms
        metadata['Masterkey'] = Masterkey
        
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
    #metadata['euc'][int(i/Step)] = result_cache['euc']
    
    if i%Skip==0:      
        result_cache['pdf'] = Kernels.Uniform(result_cache['euc'], 0.25)
        metadata['pdf'][int(i/(Step*Skip))] = result_cache['pdf']
    
    
    
    if FirstPDF == True and System['PdfStats'] == True:
        metadata['pdfpstat'] = np.zeros((Time))
        metadata['pdfkld'] = np.zeros((Time))
        metadata['pdfjsd'] = np.zeros((Time))
        FirstPDF=False
        print("First PDF calculation, brah!")
    
    """
    if System['PdfStats'] == True:
        
        result_cache['pdfpstat'] = Dist_Stats.PStat(metadata['pdf'][:,1,:], i)
        result_cache['pdfkld'] = Dist_Stats.Kullback(metadata['pdf'][:,1,:], i)
        result_cache['pdfjsd'] = Dist_Stats.JSD(metadata['pdf'][:,1,:], i)
        
        metadata['pdfpstat'][int(i/Step)] = result_cache['pdfpstat']
        metadata['pdfkld'][int(i/Step)] = result_cache['pdfkld']
        metadata['pdfjsd'][int(i/Step)] = result_cache['pdfjsd']
    """
        
        
        
    
    """
    Robert:
        For the homo, it's a bit trickier for quick introduction stuff.
    """
    
    """
    Robert:
        Pro-tip. Be sure to change the final argument in Homo(frame, positions, elements, specie,)
        to whichever compound you want to work with.
        
    """
    
    
    """
    result_cache['homo'] = Homo(i, result_cache['pos'], metadata['Elements'], System['Homo'])
    
    if i == Start:
        metadata['homo'] = np.zeros((Time,  len(result_cache['homo'])))
        continue
    metadata['homo'][int(i/Step)] = result_cache['homo']
    
    result_cache['pdfhomo'] = Kernels.Uniform(result_cache['homo'], 0.25)
    metadata['pdfhomo'][int(i/Step)] = result_cache['pdfhomo']
    
    
    if i == Start and System['HomoStats'] == True:
        metadata['homopstat'] = np.zeros((Time))
        metadata['homokld'] = np.zeros((Time))
        metadata['homojsd'] = np.zeros((Time))
        
        continue
    if System['HomoStats'] == True:
        
        result_cache['homopstat'] = Dist_Stats.PStat(metadata['pdfhomo'][:,1,:], i)
        result_cache['homokld'] = Dist_Stats.Kullback(metadata['pdfhomo'][:,1,:], i)
        result_cache['homojsd'] = Dist_Stats.JSD(metadata['pdfhomo'][:,1,:], i)
        
        metadata['homopstat'][int(i/Step)] = result_cache['homopstat']
        metadata['homokld'][int(i/Step)] = result_cache['homokld']
        metadata['homojsd'][int(i/Step)] = result_cache['homojsd']
    """
        
        
        
    result_cache['rdf'] = RDF(i, result_cache['pos'], 500, 10.0)
    metadata['rdf'][int(i/Step)] = result_cache['rdf']   
    
    
    if FirstRDF == True and System['RdfStats'] == True:     
        metadata['rdfpstat'] = np.zeros((Time))
        metadata['rdfkld'] = np.zeros((Time))
        metadata['rdfjsd'] = np.zeros((Time))
        FirstRDF = False
        print("First RDF calculation, dudebro.")


    """"
    if System['RdfStats'] == True:
        
        result_cache['rdfpstat'] = Dist_Stats.PStat(metadata['rdf'][:,1,:], i)
        result_cache['rdfkld'] = Dist_Stats.Kullback(metadata['rdf'][:,1,:], i)
        result_cache['rdfjsd'] = Dist_Stats.JSD(metadata['rdf'][:,1,:], i)
        
        metadata['rdfpstat'][int(i/Step)] = result_cache['rdfpstat']
        metadata['rdfkld'][int(i/Step)] = result_cache['rdfkld']
        metadata['rdfjsd'][int(i/Step)] = result_cache['rdfjsd']    
    """
    
    if i%Skip==0:
        result_cache['adj'] = Adjacency_Matrix(i, result_cache['pos'], R_Cut=3.0)
        metadata['adj'][int(i/(Step*Skip))] = result_cache['adj']

    First = False
    i += Step
    T1=time.time()

    print("Step ", i," took ", T1-T0, " seconds.")



"""
Robert: 
    The way in which I've written this function makes it a pain to call at each frame.
    So instead, it can be called at the end, for now and then spruced up for the big release!
"""

tick = time.time()
metadata['cna'] = CNA_Sig_Frame(filename, Masterkey, R_Cut=3.0, Frames = End, Skip = Step)
tock = time.time()
print("CNA signatures for trajectory found in ", tock - tick, "seconds.")

if System['CnaStats'] == True:
    metadata['cnapstat'] = np.zeros((Time))
    metadata['cnakld'] = np.zeros((Time))
    metadata['cnajsd'] = np.zeros((Time))


"""
        
    result_cache['cnapstat'] = Dist_Stats.PStat(metadata['cna'][i], i)
    result_cache['cnakld'] = Dist_Stats.Kullback(metadata['cna'][i], i)
    result_cache['cnajsd'] = Dist_Stats.JSD(metadata['cna'][i], i)
        
    metadata['cnapstat'][int(i/Step)] = result_cache['cnapstat']
    metadata['cnakld'][int(i/Step)] = result_cache['cnakld']
    metadata['cnajsd'][int(i/Step)] = result_cache['cnajsd']
    
"""


with open("CuMetaTest.csv", "wb") as write:
    pickle.dump(metadata,write)