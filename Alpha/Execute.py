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


def Process(System = None, Quantities=None):

    tick = time.time()
    
    metadata={}
    
    result_cache={}
    
    filename = System['base_dir']+System['movie_file_name']

    print('Reading from the %s file.' %(filename), "\n")


    """
    In the full version, the user will either be called for these arguments or may simply submit a script.
    I'd like for both versions to be effective. But, for now, it shall remain hard-coded to faciliate 
    debugging and support transparency.
 
    """
    
    
    try:
        System['Start']
        if type(System['Start']) is not int:
            Start = 0
            print('Bad value set for initial frame. Start has been set to 0 by default. Please set an integer value in the future', "\n")
        else:

            Start = System['Start']
            print('Initial frame at %s.' %(Start), "\n")
        
    except KeyError:
        Start = 0
        print('No value set for initial frame. Start has been set to 0 by default.', "\n")
    
    metadata['Start'] = Start
    
    try:
        System['End']
        if type(System['End']) is not int:
            End  = len(read(filename, index= ':'))
            print('Bad value set for final frame. End has been set to %s, the final frame in this trajectory.'
                  'Please set an integer value in the future.' %(End), "\n")
        elif System['End']<Start:
            End  = len(read(filename, index= ':'))
            print('Bad value set for final frame. End has been set to %s, the final frame in this trajectory.'
                  'Please set a value greater than your start frame in the future.' %(End), "\n")
        
        else: 
            End = System['End']
            print('Final frame set to %s.' %(End), "\n")
        
    except KeyError:
        End  = len(read(filename, index= ':'))
        print('No value set for final frame. End has been set to %s, the final frame in this trajectory.'%(End), "\n")
        
    metadata['End'] = End
            
    try:
        System['Step']
        if type(System['Step']) is not int:
            Step = 1
            print('Bad value set for Step. This has been set to 1 by default. Please set an integer value in the future', "\n")
        else:
            Step = System['Step']
            print('Step set to %s.' %(Step), "\n")
    except KeyError:
        Step = 1
        print('No value set for Step. The default of 1 has been used instead.', "\n")
        
    metadata['Step'] = Step
    
    try:
        System['Skip']
        if type(System['Skip']) is not int:
            Skip = int(End-Start)/25.0
            print('Bad value set for Skip. This has been set to %s such that R_Cut will be evaluated roughly every 25 frames.'
                  'Be aware that this may slow down your processing considerably.' %(Skip), "\n")
        else:
            Skip = System['Skip']
            print('Skip has been set to %s.' %(Skip), "\n")
    except KeyError:
        Skip = int(End-Start)/25.0
        print('No value set for Skip. This has been set to %s such that R_Cut will be evaluated roughly every 25 frames.'
                  'Be aware that this may slow down your processing considerably.' %(Skip), "\n")
        
    metadata['Skip'] = Skip

    Time=int((End-Start)/Step)
    
    print("Reading trajectory from frames %s to %s with an increment of %s." %(Start, End, Step), "\n")
    print('The PDF and, by extension, R_Cut will be evaluated every %s frames.' %(Skip), "\n")

    Dataset = read(filename, index = 0)

    all_atoms = Dataset.get_chemical_symbols()


    try: 
        System['energy_file_name']

        energy = np.loadtxt(System['base_dir']+System['energy_file_name'])
        print('Reading from the %s file.' %(System['energy_file_name']), "\n")
    except KeyError:
        print("No energy file given, no quantities related to energy will be evaluated.", "\n")
        System['SimTime'] = False; System['EPot'] = False; System['ETot'] = False; System['EKin'] = False
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
            This function also sets a new R_Cut each time it is called and calculated.
        
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
            Quantities[x]; print("Calculating the %s." %(x), "\n"); globals()[x] = True
            Quantities[x] = np.empty((Time,), dtype=object)
            if x is 'pdf':
                Quantities[x] = np.empty((int((Time*Step)/(Skip)),), dtype=object)
        except KeyError:
            print("Will not calculate %s in this run." %(x), "\n"); globals()[x] = False
        


    First = True; FirstPDF = True; FirstRDF = True



    if bool(globals()['cna']) is True:
        tick = time.time()

        # Load the simulation dataset to be analyzed.
        Masterkey = []

        # Create bonds.

    
        Masterkey = []
    
        print("Initialising CNA environment took %.3f seconds." %(time.time()-tick), "\n")


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
            print("Initialising Metadata took %.3f seconds." %(time.time() - tick),"\n")
            print("This system contains", NAtoms, "atoms","\n",
                  "consisting of", Species, "as present atomic species.","\n")
        
        if bool(globals()['pos']) is True:
            result_cache['pos'] = all_positions
            metadata['pos'][int(i/Step)] = result_cache['pos']
    

        if bool(globals()['rdf']*globals()['pos']) is True:
            tick=time.time()
            result_cache['rdf'] = RDF(i, result_cache['pos'], 100, 10.0)
            metadata['rdf'][int(i/Step)] = result_cache['rdf'] 
            print("RDF calculated in %.3f seconds." %(time.time()-tick),"\n")
    
    

    
        if i%Skip==0 and bool(globals()['euc']*globals()['pdf']*globals()['pos']) == True:
            tick = time.time()
            result_cache['euc'] = Euc_Dist(i, result_cache['pos'])
            #metadata['euc'][int(i/Step)] = result_cache['euc']   #This is a very heavy list to carry around and the user is not advised to save it unless they DESPERATELY need it
            result_cache['pdf'] = Kernels.Uniform(result_cache['euc'], 0.25)
            metadata['pdf'][int(i/(Step*Skip))] = result_cache['pdf']
            R_Cut = result_cache['pdf'][-1]
            print("PDF evaluated at frame", i,"and took %.3f seconds to calculate." %(time.time() - tick),"\n")
            print("R_Cut is now set to %s." %(R_Cut))
    
    
    
        if bool(FirstPDF*globals()['pdf']*globals()['pos']) is True:
            tick = time.time()
            metadata['pdfpstat'] = np.zeros((Time))
            metadata['pdfkld'] = np.zeros((Time))
            metadata['pdfjsd'] = np.zeros((Time))
            FirstPDF=False
            print("PDF statistics initialisation performed in %.3f seconds." %(time.time()-tick),"\n")
    
    
        if bool(globals()['cna']) is True:
            tick=time.time()
            result_cache['cna'] = Frame_CNA(i, R_Cut, Masterkey, filename)[0]
            metadata['cna'][int(i/Step)] = list(result_cache['cna'].items())
            metadata['cna'][int(i/Step)].sort()
            print("CNA signatures calculated in %.3f seconds." %(time.time()-tick),"\n")
        

    
    
        if bool(FirstRDF*globals()['rdf']*globals()['pos']) is True:     
            metadata['rdfpstat'] = np.zeros((Time))
            metadata['rdfkld'] = np.zeros((Time))
            metadata['rdfjsd'] = np.zeros((Time))
            FirstRDF = False
            print("RDF statistics initialisation performed in %.3f seconds." %(time.time()-tick),"\n")

        if bool(globals()['adj']*globals()['pos']) is True:
            tick=time.time()
            result_cache['adj'] = Adjacency_Matrix(result_cache['pos'], R_Cut)
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
        
        Robert: 
            
            Below are many of the statistical analysis tools being called on various distributions.
            They are commented out at the moment as I feel it more useful for the user to use the
            functions already defined in a way that see fit for their data.
        
        
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
        
        result_cache['pdfhomo'] = Kernels.Uniform(result_cache['homo'], 0.25)globals(X[b[1]])
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
        metadata['cnajsd'] = np.zeros((Time))print("R_Cut is now set to %s." %(R_Cut))
    
    
    
            
        result_cache['cnapstat'] = Dist_Stats.PStat(metadata['cna'][i], i)
        result_cache['cnakld'] = Dist_Stats.Kullback(metadata['cna'][i], i)
        result_cache['cnajsd'] = Dist_Stats.JSD(metadata['cna'][i], i)
            
        metadata['cnapstat'][int(i/Step)] = result_cache['cnapstat']
        metadata['cnakld'][int(i/Step)] = result_cache['cnakld']
        metadata['cnajsd'][int(i/Step)] = result_cache['cnajsd']
        
        
    #############################################################################################
    
    Robert:
        And now we check to see if the users wishes to evaluate any of the quantities
        from the energy file and add them to the metadata.
    
    """


    if bool(globals()['SimTime']) is True:
        metadata['time'] = energy[:,0]

    if bool(globals()['EPot']) is True:
        metadata['EPot'] = energy[:,1]
    
    if bool(globals()['ETot']) is True:
        metadata['ETot'] = energy[:,2]
    
    if bool(globals()['EKin']) is True:
        metadata['EKin'] = energy[:,3]
    
    if bool(globals()['EDelta']) is True:
        metadata['EDelta'] = energy[:,4]

    if bool(globals()['MeanETot']) is True:
        metadata['MeanETot'] = energy[:,5]
    
    if bool(globals()['Temp']) is True:
        metadata['Temp'] = energy[:,6]

    return metadata
