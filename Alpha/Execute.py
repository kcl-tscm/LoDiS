#Load the required modules for post-processing

from Adjacent import *
from Kernels import *
from ase.io import read
import time
from DistFuncs import *
from SampleCNA import *
from AGCN import *

import numpy as np
import csv

metadata={}
result_cache={}


"""
In the full version, the user will either be called for these arguments or may simply submit a script.
I'd like for both versions to be effective. But, for now, it shall remain hard-coded to faciliate 
debugging and support transparency.

"""
tick = time.time()

Start: int = 0
End: int = 800
Step: int = 50
Skip=100

R_Cut=3.0 #Don't get too used to this guy being here. I want to make him obsolete with the PDDF


"""
Robert:
    Change the arguments in the System dictionary to 
    reflect the directory and files you are working with.
"""

System = {
        'base_dir' : '/home/k1899676/Documents/PhD/SampleTraj/',
        'movie_file_name' : 'JanusMeltMovie1.xyz',
        'energy_file_name' : 'JanusMeltEn1.out',
        
        'Homo' : 'Cu',
        
        'Start' : 0, 'End' : 800, 'Step' : 50, 'Skip' : 50,
        'PdfStats' : True, 'HomoStats' : True, 'RdfStats' : True, 'CnaStats' : True,
        
        #Below are quantities related to the energy file
        
        'time': True, 'EPot': True, 'ETot' : True, 
        'EKin' : True, 'EDelta' : True, 'MeanETot' : True, 'Temp' : True
        }


settings: dict = None #Ignore this for the time being. It may yet be of use.

filename = System['base_dir']+System['movie_file_name']

Time=int((End-Start)/Step)

Dataset = read(filename, index = 0)

all_atoms = Dataset.get_chemical_symbols()

Quantities = {
        'euc' : np.zeros((Time,int(len(all_atoms)*(len(all_atoms)-1)/2))),
        'rdf' : np.empty((Time,), dtype=object), 
        'pos' : np.empty((Time,), dtype=object), 
        'cna' : np.empty((Time,), dtype=object), 
        'adj' : np.empty((int(Time),), dtype=object), 
        'pdf' : np.empty((int(Time*Step/Skip+1),), dtype=object),
        #'pdfhomo' : np.empty((Time,), dtype=object),
        'agcn' : np.empty((Time,), dtype=object),
        'nn' : np.empty((Time,), dtype=object),
        'time': None, 'EPot': None, 'ETot' : None, 
        'EKin' : None, 'EDelta' : None, 'MeanETot' : None, 'Temp' : None
        }



settings: dict = None #Ignore this for the time being. It may yet be of use.

filename = System['base_dir']+System['movie_file_name']

try: 
    System['energy_file_name']

    energy = np.loadtxt(System['base_dir']+System['energy_file_name'])
    print('Reading from the %s file.' %(System['energy_file_name']))
except KeyError:
    print("No energy file given, no quantities related to energy will be evaluated.", "\n")
    System['Time'] = False; System['EPot'] = False; System['ETot'] = False; System['EKin'] = False
    System['EDelta'] = False; System['MeanETot'] = False; System['Temp'] = False


Dataset = read(filename, index = 0)
all_atoms = Dataset.get_chemical_symbols()

print("Initialising system environment took %.3f seconds." %(time.time()-tick), "\n")


"""

Robert:
    
    I know that this quantities dictionary is a mess, the general idea is that this will
    be filled up with information from the class implementation 'automatically'. Because
    this has been thrown together in an afternoon with all of the calculators rejigged for
    ease of calling as opposed to efficiency, it looks ugly as sin! Sorry, team.
    
    
    Supported quantities as of this release:
        
        Euclidean distance: euc - Pairwise distance between all atoms
        
        RDF: rdf - radial distribution function
        
        Common Neighbour Analysis: cna - all signatures and the number of observed counts
        
        Adjacency matrix: adj - Sparse matrix of truth elements regarding whether or not two atoms are neighbours
        
        Pair distance distribution function: pdf - Kernel densiy estimator (uniform approximation) for the pdf.
        This will be updated in the near future to extrapolate a new R_Cut each time it is called. - NOT YET SUPPORTED -
        
        Homo atoms PDDF: pdfhomo: The PDDF for similar atoms only. - NOT YET FULLY TESTED. DO NOT USE -
        
        Atop generalised coordination number: agcn - ask Fra
        
        Nearest neighbours: nn - Number of nearest neighbours each atom has. 
    
"""

Supported=[
        'euc', 'rdf', 'cna', 'adj', 'pdf', 'pdfhomo', 'agcn', 'nn', 'pos',
        'SimTime', 'EPot', 'ETot', 'EKin', 'EDelta', 'MeanETot', 'Temp'
           ]

for x in Supported:
    try:
        Quantities[x]; print("Calculating the %s" %(x), "\n"); globals()[x] = True
    except KeyError:
        print("Will not calculate %s in this run." %(x), "\n"); globals()[x] = False
        
  
"""
try:
    Quantities['euc']; print('Euclidean distances will be calculated but not saved.', "\n"); Euc = True
except KeyError:
    print("No distances will be calculated. This will make you unable to calculate any PDDFs.", "\n" 
          "and R_Cut will be set to a default of 3.0 Angstrom.", "\n")

try: 
    Quantities['cna']; print('CNA signatures will be calculated.', "\n"); CNA = True
except KeyError:
    print('CNA signatures will not be calculated.', "\n"); CNA = False

try: 
    Quantities['pdf']
    if Euc is True:
        print('PDFs will be calculated.', "\n"); FirstPDF = True; PDF = True
    else:
        print("The PDF cannot be evaluated without calculating distances." 
              "No PDF will be evaluated this run. Setting R_Cut to 3.0 Angstrom by default.", "\n"); R_Cut = 3.0
except KeyError:
    print('PDFs will not be calculated.', "\n"); FirstPDF = False; PDF = False
    print('R_Cut will be set to a default of 3.0 Angstrom', '\n'); R_Cut = 3.0

try: 
    Quantities['rdf']; print('RDFs will be calculated.', "\n"); FirstRDF = True; RDF = True
except KeyError:
    print('RDFs will not be calculated.',"\n"); FirstRDF = False; RDF = False
"""


First = True; FirstPDF = True; FirstRDF = True



if bool(globals()['cna']) is True:
    tick = time.time()

    # Load the simulation dataset to be analyzed.
    pipeline = import_file(filename)
    Masterkey = []

    # Create bonds.

    """

    Robert:
        In principal, this line below should be called every time the R_Cut is updated
        so that more accurate signatures may be found. This, in theory, may be related 
        to pressures, thermal expansions etc...
    
    """

    pipeline.modifiers.append(CreateBondsModifier(cutoff = R_Cut))
    pipeline.modifiers.append(CommonNeighborAnalysisModifier(
        mode = CommonNeighborAnalysisModifier.Mode.BondBased))
    
    Masterkey = []
    
    print("Initialising CNA environment took %.3f seconds." %(time.time()-tick), "\n")

    # Compute CNA indices on the basis of the created bonds.






for i in range(Start, End, Step):
    T0=time.time()
    
    
    Dataset = read(filename, index = i)
    all_positions = Dataset.get_positions()
    all_atoms = Dataset.get_chemical_symbols()
    
    used=set()
    Species = [x for x in all_atoms if x not in used and (used.add(x) or True)]
    
    NAtoms = len(all_atoms)
    
    if i == Start:
        tick = time.time()
        metadata['Elements'] = all_atoms
        metadata['Species'] = Species
        metadata['NSpecies']=len(Species)
        metadata['NFrames'] = Time
        metadata['NAtoms'] = NAtoms
        
        for key in Quantities:
            metadata[key] = Quantities[key]
        print("Initialising Metadata took", time.time() - tick, "seconds.","\n")
        print("This system contains", NAtoms, "atoms","\n",
              "consisting of", Species,"\n",
              "across", Time, "frames of a trajectory.", "\n")
        
    if bool(globals()['pos']) is True:
        result_cache['pos'] = all_positions
        metadata['pos'][int(i/Step)] = result_cache['pos']
    
    

    
    if i%Skip==0 and bool(globals()['euc']*globals()['pdf']*globals()['pos']) == True:
        tick = time.time()
        result_cache['euc'] = Euc_Dist(i, result_cache['pos'])
        #metadata['euc'][int(i/Step)] = result_cache['euc']
        result_cache['pdf'] = Kernels.Uniform(result_cache['euc'], 0.25)
        metadata['pdf'][int(i/(Step*Skip))] = result_cache['pdf']
        print("PDF evaluated at frame", i,"and took %.3f seconds to calculate." %(time.time() - tick),"\n")
    
    
    
    if bool(FirstPDF*globals()['pdf']*globals()['pos']) is True:
        tick = time.time()
        metadata['pdfpstat'] = np.zeros((Time))
        metadata['pdfkld'] = np.zeros((Time))
        metadata['pdfjsd'] = np.zeros((Time))
        FirstPDF=False
        print("PDF statistics initialisation performed in %.3f seconds." %(time.time()-tick),"\n")
    
    
    if bool(globals()['cna']) is True:
        tick=time.time()
        result_cache['cna'] = Frame_CNA(i, pipeline, Masterkey)[0]
        metadata['cna'][int(i/Step)] = list(result_cache['cna'].items())
        metadata['cna'][int(i/Step)].sort()
        print("CNA signatures calculated in %.3f seconds." %(time.time()-tick),"\n")
        

    if bool(globals()['rdf']*globals()['pos']) is True:
        tick=time.time()
        result_cache['rdf'] = RDF(i, result_cache['pos'], 500, 10.0)
        metadata['rdf'][int(i/Step)] = result_cache['rdf']   
        print("RDF calculated in %.3f seconds." %(time.time()-tick),"\n")
    
    
    if bool(FirstRDF*globals()['rdf']*globals()['pos']) is True:     
        metadata['rdfpstat'] = np.zeros((Time))
        metadata['rdfkld'] = np.zeros((Time))
        metadata['rdfjsd'] = np.zeros((Time))
        FirstRDF = False
        print("RDF statistics initialisation performed in %.3f seconds." %(time.time()-tick),"\n")

    if bool(globals()['adj']*globals()['pos']) is True:
        tick=time.time()
        result_cache['adj'] = Adjacency_Matrix(result_cache['pos'], R_Cut=3.0)
        metadata['adj'][int(i/(Step))] = result_cache['adj']
        print("Adjacency matrix evaluated at frame", i,"in %.3f seconds." %(time.time()-tick), "\n")
    
    
    if bool(globals()['agcn']*globals()['nn']*globals()['adj']*globals()['pos']) is True:
        tick = time.time()
        Agcn, NN = agcn_generator(result_cache['adj'])
        metadata['agcn'][int(i/Step)] = Agcn; metadata['nn'][int(i/Step)] = NN
        print("AGCN and nearest neighbours found in %.3f seconds." %(time.time()-tick), "\n")
    elif bool(globals['agcn']*globals()['adj']*globals()['pos']) is True:
        tick = time.time()
        Agcn = agcn_generator(result_cache['adj'])[0]
        metadata['agcn'][int(i/Step)] = Agcn
        print("AGCN found in %.3f seconds." %(time.time()-tick), "\n")
    elif bool(globals['nn']*globals()['adj']*globals()['pos']) is True:
        tick = time.time()
        NN = agcn_generator(result_cache['adj'])[1]
        metadata['nn'][int(i/Step)] = NN
        print("Nearest neighbours found in %.3f seconds." %(time.time()-tick), "\n")
        

    """"
    if System['RdfStats'] == True:
        
        result_cache['rdfpstat'] = Dist_Stats.PStat(metadata['rdf'][:,1,:], i)
        result_cache['rdfkld'] = Dist_Stats.Kullback(metadata['rdf'][:,1,:], i)
        result_cache['rdfjsd'] = Dist_Stats.JSD(metadata['rdf'][:,1,:], i)
        
        metadata['rdfpstat'][int(i/Step)] = result_cache['rdfpstat']
        metadata['rdfkld'][int(i/Step)] = result_cache['rdfkld']
        metadata['rdfjsd'][int(i/Step)] = result_cache['rdfjsd']    

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

    if System['PdfStats'] == True:
        
        result_cache['pdfpstat'] = Dist_Stats.PStat(metadata['pdf'][:,1,:], i)
        result_cache['pdfkld'] = Dist_Stats.Kullback(metadata['pdf'][:,1,:], i)
        result_cache['pdfjsd'] = Dist_Stats.JSD(metadata['pdf'][:,1,:], i)
        
        metadata['pdfpstat'][int(i/Step)] = result_cache['pdfpstat']
        metadata['pdfkld'][int(i/Step)] = result_cache['pdfkld']
        metadata['pdfjsd'][int(i/Step)] = result_cache['pdfjsd']



    Robert:
        Pro-tip. Be sure to change the final argument in Homo(frame, positions, elements, specie,)
        to whichever compound you want to work with.
        
    """    
    
    
    

    First = False
    T1=time.time()

    print("Step", i,"took", T1-T0, "seconds.", "\n")
    i += Step


Masterkey.sort()
metadata['masterkey']=Masterkey

"""
Robert: 
    The way in which I've written this function makes it a pain to call at each frame.
    So instead, it can be called at the end, for now and then spruced up for the big release!



if System['CnaStats'] == True:
    metadata['cnapstat'] = np.zeros((Time))
    metadata['cnakld'] = np.zeros((Time))
    metadata['cnajsd'] = np.zeros((Time))



        
    result_cache['cnapstat'] = Dist_Stats.PStat(metadata['cna'][i], i)
    result_cache['cnakld'] = Dist_Stats.Kullback(metadata['cna'][i], i)
    result_cache['cnajsd'] = Dist_Stats.JSD(metadata['cna'][i], i)
        
    metadata['cnapstat'][int(i/Step)] = result_cache['cnapstat']
    metadata['cnakld'][int(i/Step)] = result_cache['cnakld']
    metadata['cnajsd'][int(i/Step)] = result_cache['cnajsd']
    
"""

"""

Robert:
    And now we check to see if the users wishes to evaluate any of the quantities
    from the energy file and add them to the metadata.

"""

'SimTime', 'EPot', 'ETot', 'EKin', 'EDelta', 'MeanETot', 'Temp'

if bool(globals()['SimTime']) is True:
    metadata['time'] = globals()['energy'][:,0]

if bool(globals()['EPot']) is True:
    metadata['EPot'] = globals()['energy'][:,1]
    
if bool(globals()['ETot']) is True:
    metadata['ETot'] = globals()['energy'][:,2]
    
if bool(globals()['EKin']) is True:
    metadata['EKin'] = globals()['energy'][:,3]
    
if bool(globals()['EDelta']) is True:
    metadata['EDelta'] = globals()['energy'][:,4]

if bool(globals()['MeanETot']) is True:
    metadata['MeanETot'] = globals()['energy'][:,5]
    
if bool(globals()['Temp']) is True:
    metadata['Temp'] = globals()['energy'][:,6]

with open((System['base_dir']+'Tested.csv'), 'w') as f:
    for key in metadata.keys():
        f.write("%s,%s\n"%(key,metadata[key]))