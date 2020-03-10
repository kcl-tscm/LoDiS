#Load the required modules for post-processing
import Adjacent
import Kernels
from ase.io import read
import DistFuncs
import AGCN

import numpy as np
import time

def Process(System = None, Quantities=None):

    tick = time.time(); BigT = time.time()
    
    Supported=[
            'rdf', 'cna', 'adj', 'pdf', 'pdfhomo', 'agcn', 'nn', 'CoM',
            'SimTime', 'EPot', 'ETot', 'EKin', 'EDelta', 'MeanETot', 'Temp'
               ]
    
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
    
    
    try:
        System['UniformPDF']
        if System['UniformPDF'] is False:
            PDF = Kernels.Kernels.Gauss
            print('The set method for calculating the PDF is with a Gaussian kernel function. \n Be aware that this method'
                  'is slower than using a Uniform kernel. However; the distribution will be smoother.', "\n")
            metadata['pdftype'] = 'Gauss'
            try:
                System['Band']
                if bool(type(System['Band']) is float or int):
                    Band = System['Band']
                    print('Bandwidth for the Kernel Density Estimator set to %s.' %(Band), "\n")
                    metadata['Band'] = Band
                else:
                    Band = 0.05
                    print('Bad value set for the Kernel function bandwidth. \n Defaulting to % for the Gaussian Kernel Density Estimator.' %(Band), "\n")
                    metadata['Band'] = Band
            except KeyError:
                Band = 0.05
                print('Default setting for the Gaussian Kernel Density Estimator is set to %s.' %(Band), "\n")
                metadata['Band'] = Band
                
        else:
            PDF = Kernels.Kernels.Uniform
            print('The selected method for calculating the PDF is with a Uniform kernel function. \n Be aware that this method'
                  'may yield non-smooth distributions for certain structures. However; this is a much faster calculator.', "\n")
            metadata['pdftype'] = 'Uniform'
            try:
                System['Band']
                if bool(type(System['Band']) is float or int):
                    Band = System['Band']
                    print('Bandwidth for the Kernel Density Estimator set to %s.' %(Band), "\n")
                    metadata['Band'] = Band
                else:
                    Band = 0.25
                    print('Bad value set for the Kernel function bandwidth. \n Defaulting to % for the Uniform Kernel Density Estimator.' %(Band), "\n")
                    metadata['Band'] = Band
            except KeyError:
                Band = 0.25
                print('Default setting for the Uniform Kernel Density Estimator is set to %s.' %(Band), "\n")
                metadata['Band'] = Band
                
    except KeyError:
        PDF = Kernels.Kernels.Uniform
        print('The default method for calculating the PDF is with a Uniform kernel function. \n Be aware that this method'
              'may yield non-smooth distributions for certain structures. However; this is a much faster calculator.',"\n")
        metadata['pdftype'] = 'Uniform'
        try:
            System['Band']
            if bool(type(System['Band']) is float or int):
                Band = System['Band']
                print('Bandwidth for the Kernel Density Estimator set to %.' %(Band), "\n")
                metadata['Band'] = Band
            else:
                Band = 0.25
                print('Bad value set for the Kernel function bandwidth. \n Defaulting to % for the Uniform Kernel Density Estimator.' %(Band), "\n")
                metadata['Band'] = Band
        except KeyError:
            Band = 0.25
            print('Default setting for the Uniform Kernel Density Estimator is set to %.' %(Band), "\n")
            metadata['Band'] = Band

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
        
    for x in Supported:
        try:
            Quantities[x]; print("Calculating the %s." %(x), "\n"); globals()[x] = True
            Quantities[x] = np.empty((Time,), dtype=object)
            if x == 'pdf':
                Quantities[x] = np.empty((int((Time*Step)/(Skip)),), dtype=object)
            if x == 'rdf':
                Quantities[x] = np.empty((int((Time*Step)/(Skip)),), dtype=object)
        except KeyError:
            print("Will not calculate %s in this run." %(x), "\n"); globals()[x] = False


    Dataset = read(filename, index = 0)
    all_atoms = Dataset.get_chemical_symbols()

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
        
            Atop generalised coordination number: agcn - ask Fra
        
            Nearest neighbours: nn - Number of nearest neighbours each atom has. 
    
    """
            
    try: 
        System['HCStats']
        if bool(System['HCStats']) is not False:
            Quantities['h'] = np.empty((Time,), dtype=object); globals()['h'] = True
            Quantities['c'] = np.empty((Time,), dtype=object); globals()['c'] = True
            print("Will be calculating and evaluating collectednes and concertednes of cluster rearrangement.", "\n")
        else:
            print("Will not be calculating collectednes or concertednes of cluster rearrangements.", "\n")
    except KeyError:
        print("Will not be calculating collectednes or concertednes of cluster rearrangements.", "\n")
        
    print("Initialising system environment took %.3f seconds." %(time.time()-tick), "\n")

    if bool(globals()['cna']) is True:
        import SampleCNA
        tick = time.time()

        # Load the simulation dataset to be analyzed.
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
            
            print("Checking user input for calculating homo properties in this run.", "\n")
            try:
                System['Homo']
                
                if System['Homo'] is None:
                    try:
                        System['HomoQuants']
                        if System['HomoQuants'] is None:
                            print("No bimetallic properties for homo species will be calculated in this run.", "\n")
                        else:
                            System['Homo'] = metadata['Species']
                            print("No homo atom species requested, but you wish to calculate bimetallic homo properties." 
                                  "\n Instead we shall calculate homo properties for %s and hetero properties for the system." %(metadata['Species']), "\n")
                       
                            for x in System['HomoQuants']:
                                for y in System['Homo']:
                                    Quantities[x+y] = np.empty((Time,), dtype=object); globals()[x+y] = True
                                    print("Calculating %s as a homo property." %(x+y), "\n")
                                    if 'PDF' in x:
                                        Quantities[x+y] = np.empty((int((Time*Step)/(Skip)),), dtype=object)
                                    elif 'RDF' in x:
                                        Quantities[x+y] = np.empty((int((Time*Step)/(Skip)),), dtype=object)
                    except KeyError:
                        print("Will not be calculating any homo properties this run." , "\n") 
                        
                        
                elif False in [x not in metadata['Species'] for x in System['Homo']]:
                    print("Specie entered in homo not found in the system. The only observed species are %s and you have requested to observe %s." 
                          "\n Defaulting to atoms found in the system for evaluation." %(metadata['Species'], System['Homo']), "\n")
                    System['Homo'] = metadata['Species']
                    try:
                        System['HomoQuants']
                        if System['HomoQuants'] is None:
                            print("No homo properties will be calculated in this run.", "\n")
                        else:
                            for x in System['HomoQuants']:
                                for y in System['Homo']:
                                    Quantities[x+y] = np.empty((Time,), dtype=object); globals()[x+y] = True
                                    print("Calculating %s as a homo property." %(x+y), "\n")
                                    if 'PDF' in x:
                                        Quantities[x+y] = np.empty((int((Time*Step)/(Skip)),), dtype=object)
                                    elif 'RDF' in x:
                                        Quantities[x+y] = np.empty((int((Time*Step)/(Skip)),), dtype=object)
                    except KeyError:
                        print("Will not be calculating any homo properties this run as no qauntities have been given to calculate." , "\n") 
                        
                        
                else:
                    print("Homo atom properties will be caluclated for %s in this run." %(System['Homo']), "\n")
                    try:
                        System['HomoQuants']
                        if System['HomoQuants'] is None:
                            print("No bimetallic properties will be calculated in this run as none have been requested.", "\n")
                        else:
                            for x in System['HomoQuants']:
                                for y in System['Homo']:
                                    Quantities[x+y] = np.empty((Time,), dtype=object); globals()[x+y] = True
                                    print("Calculating %s as a homo property." %(x+y), "\n")
                                    if 'PDF' in x:
                                        Quantities[x+y] = np.empty((int((Time*Step)/(Skip)),), dtype=object)
                                    elif 'RDF' in x:
                                        Quantities[x+y] = np.empty((int((Time*Step)/(Skip)),), dtype=object)
                    except KeyError:
                        print("Will not be calculating any homo properties this run." , "\n")
                        
                        
            except KeyError:
                print("No homo atoms have been requested for calculation. Checking if bimetallic properties have been requested.", "\n")
                
                try:
                    System['HomoQuants']
                    if System['HomoQuants'] is None:
                        print("No homo properties have been requested, either. Continuing to calculate whole system properties, only.", "\n")
                    else:
                        print("You have requested to calculate %s while not asking for any atoms. Defaulting to considering all species identified in the system." %(System['HomoQuants']), "\n")
                        System['Homo'] = metadata['Species']
                        
                        for x in System['HomoQuants']:
                            for y in System['Homo']:
                                Quantities[x+y] = np.empty((Time,), dtype = object); globals()[x+y] = True
                                print("Calculating %s as a homo property." %(x+y), "\n")
                                if 'PDF' in x:
                                    Quantities[x+y] = np.empty((int((Time*Step)/(Skip)),), dtype = object)
                                elif 'RDF' in x:
                                    Quantities[x+y] = np.empty((int((Time*Step)/(Skip)),), dtype=object)
                                
                except KeyError:
                    print("No homo quantities have been requested, either.", "\n")
            
            print("Finished evaluating user input for homo atomic properties." , "\n")
            
            print("Checking user input for hetero atomic species.", "\n")
            
            try:
                System['Hetero']
                if System['Hetero'] is not True:
                    print("Bad input detected for the 'Hetero' argument'. \n Checking if the user has requested hetero quantities to calculate.", "\n")
                    try: 
                        System['HeteroQuants']
                        if System['HeteroQuants'] is None:
                            print("Bad input variable decalred for calculating hetero quantities. Nothing hetero will happen here, today!", "\n")
                        else:
                            print("User has requested hetero quantities without specifying the desire to do so. We shall assume that this is an error and calculate anyway.", "\n")
                            System['Hetero'] = True
                            for x in System['HeteroQuants']:
                                Quantities[x] = np.empty((Time,), dtype = object); globals()[x] = True
                                print("Calculating %s as a hetero property." %(x), "\n")
                                if 'PDF' in x:
                                    Quantities[x] = np.empty((int((Time*Step)/(Skip)),), dtype = object)
                                elif 'RDF' in x:
                                    Quantities[x] = np.empty((int((Time*Step)/(Skip)),), dtype =object)
                    except KeyError:
                        print("No hetero quantities requested and so none shall be calculated.", "\n")
                
                else:
                    print("Hetero quantities have been requested by the user.", "\n")
                    try:
                        System['HeteroQuants']
                        if System['HeteroQuants'] is None:
                            print("Bad input variable decalred for calculating hetero quantities. Nothing hetero will happen here, today!", "\n")
                        else:
                            print("User has requested hetero quantities.", "\n")
                            System['Hetero'] = True
                            for x in System['HeteroQuants']:
                                Quantities[x] = np.empty((Time,), dtype = object); globals()[x] = True
                                print("Calculating %s as a hetero property." %(x), "\n")
                                if 'PDF' in x:
                                    Quantities[x] = np.empty((int((Time*Step)/(Skip)),), dtype = object)
                                elif 'RDF' in x:
                                    Quantities[x] = np.empty((int((Time*Step)/(Skip)),), dtype =object)
                    except KeyError:
                        print("No hetero quantities requested and so none shall be calculated.", "\n")
            except KeyError:
                print("No input variable declared for 'Hetero' calculations. Checking if user has requested quantities without specifying the wish to calculate.", "\n")
                try:
                    System['HeteroQuants']
                    if System['HeteroQuants'] is None:
                        print("Bad input variable decalred for calculating hetero quantities. Nothing hetero will happen here, today!", "\n")
                    else:
                        print("User has requested hetero quantities.", "\n")
                        System['Hetero'] = True
                        for x in System['HeteroQuants']:
                            Quantities[x] = np.empty((Time,), dtype = object); globals()[x] = True
                            print("Calculating %s as a hetero property." %(x), "\n")
                            if 'PDF' in x:
                                Quantities[x] = np.empty((int((Time*Step)/(Skip)),), dtype =object)
                            elif 'RDF' in x:
                                Quantities[x] = np.empty((int((Time*Step)/(Skip)),), dtype =object)
                except KeyError:
                    print("No hetero quantities requested and so none shall be calculated.", "\n")
                    
                    
            print("Finished evaluating input arguments for homo/hetero calculations.", "\n")
                       
            #This block initialises the metadata
            for key in Quantities:
                metadata[key] = Quantities[key]
            print("Initialising Metadata took %.3f seconds." %(time.time() - tick),"\n")
            print("This system contains", NAtoms, "atoms","\n",
                  "consisting of", Species, "as present atomic species.","\n")

        temptime=time.time()
        result_cache['pos'] = all_positions
        result_cache['euc'] = DistFuncs.Euc_Dist(result_cache['pos'])
          
        
        #All RDF calculations performed in the following block
        if i%Skip==0 and bool(globals()['rdf']) is True: 
            result_cache['rdf'] = DistFuncs.RDF(result_cache['pos'], 100, 10.0)
            metadata['rdf'][int(i/(Step*Skip))] = result_cache['rdf'] 
            try:
                if bool(bool(System['Homo'])*bool('HoPDF' in System['HomoQuants'])) is True:
                    for x in System['Homo']:
                        result_cache['homopos'+x] = DistFuncs.get_subspecieslist(x, metadata['Elements'], result_cache['pos'])
                        metadata['HoRDF'+x][int(i/(Step*Skip))] = DistFuncs.RDF(result_cache['homopos'+x])
            except KeyError:
                pass
            try:
                if bool(bool(System['Hetero'])*globals()['HeRDF']) is True:
                    metadata['HeRDF'][int(i/(Step*Skip))] = DistFuncs.RDF(result_cache['pos'], Res=100, R_Cut=10.0, Hetero = True, 
                                                                          Species = metadata['Species'], Elements = metadata['Elements'])
            except KeyError:
                pass
    
        #All PDF calculations performed in the following block
        if i%Skip==0 and bool(globals()['pdf']) is True:
            result_cache['pdf'] = PDF(result_cache['euc'], Band)
            metadata['pdf'][int(i/(Step*Skip))] = result_cache['pdf']
            R_Cut = result_cache['pdf'][-1]
            print("R_Cut is now set to %s." %(R_Cut), "\n")
            try:
                if bool(bool(System['Homo'])*bool('HoPDF' in System['HomoQuants'])) is True:
                    for x in System['Homo']:
                        result_cache['homoed'+x] = DistFuncs.Euc_Dist(positions=result_cache['pos'], homo = True, specie = x, elements = metadata['Elements'])
                        if result_cache['homoed'+x] is not None:
                            metadata['HoPDF'+x][int(i/(Step*Skip))] = Kernels.Kernels.Uniform(result_cache['homoed'+x], Band, mon=True)
                        else:
                            pass
            except KeyError:
                pass
            try:
                if bool(System['Hetero']*globals()['HePDF']) is True:
                    result_cache['heteropos'] = DistFuncs.Hetero(result_cache['pos'], metadata['Species'][0], metadata['Elements'])
                    if result_cache['heteropos'] is not None:
                        Temp = np.concatenate(result_cache['heteropos']).ravel()
                        metadata['HePDF'] = Kernels.Kernels.Uniform(Temp, Band, mon=True)
                    else:
                        metadata['HePDF'] = None
                        print("There was an error with the heterogenous distance array. No PDF calculated for frame %s."%(i), "\n")
            except KeyError:
                pass
    
        #This block evaluates all of the CoM calculations
        if bool(globals()['CoM']) is True:
            metadata['CoM'] = DistFuncs.get_CoM(result_cache['pos'])
        try:
            if bool(bool(System['Homo'])*bool('CoM' in System['HomoQuants'])) is True:
                for x in System['Homo']:
                    metadata['CoM'+x][int(i/Step)] = DistFuncs.CoM_Dist(positions = result_cache['pos'], homo=True, specie = x, elements = metadata['Elements'])
        except KeyError:
            pass
 
        #This block calculates the CNA signatures for the whole system, only
        if bool(globals()['cna']) is True:
            result_cache['cna'] = SampleCNA.Frame_CNA(i, R_Cut, Masterkey, filename)[0]
            metadata['cna'][int(i/Step)] = list(result_cache['cna'].items())
            metadata['cna'][int(i/Step)].sort()
        
        
        #This block evaluates the adjacency matrices for the whole system, homo pair(s), & hetero atoms 
        if bool(globals()['adj']) is True:
            result_cache['adj'] = Adjacent.Adjacency_Matrix(result_cache['pos'], result_cache['euc'], R_Cut)
            metadata['adj'][int(i/(Step))] = result_cache['adj']
        try:
            if bool(bool(System['Homo'])*bool('HoAdj' in System['HomoQuants'])) is True:
                for x in System['Homo']:
                    result_cache['HomoED'+x] = DistFuncs.Euc_Dist(result_cache['pos'], homo = True, specie = x, elements = metadata['Elements'])
                    
                    metadata['HoAdj'+x] = Adjacent.get_coordination(Adjacent.Adjacency_Matrix(
                                                                                               DistFuncs.get_subspecieslist
                                                                                               (
                                                                                               x, metadata['Elements'], result_cache['pos']
                                                                                               ),
                                                                                               result_cache['HomoED'+x], R_Cut) )
        except KeyError:
            pass
        try:
            if bool(System['Hetero']*globals()['HeAdj']) is True:
                result_cache['HeDist'] = DistFuncs.Hetero(result_cache['pos'], metadata['Species'][0], metadata['Elements'])
                if result_cache['heteropos'] is not None:
                    metadata['HeAdj'][int(i/Step)] = Adjacent.get_coordination_hetero(result_cache['HeDist'], R_Cut)
                else:
                    metadata['HeAdj'] = None
                    print("There was an error with hetero positions, no respective adjacency matrix calculated for frame %s." %(i), "\n")
        except KeyError:
            pass
        
        
        #This  block calculates the concertedness and collectivity of atom rearrangements    
        if bool(System['HCStats']*i) is not False:
            result_cache['r'] = Adjacent.R(result_cache['adj'], metadata['adj'][int(i/(Step))-1])
            metadata['h'][int(i/(Step))-1] = Adjacent.Collectivity(result_cache['r'])
            if not(i<3):
                metadata['c'][int(i/(Step))-2] = Adjacent.Concertedness(metadata['h'][int(i/Step)-1], metadata['h'][int(i/(Step))-3])
    
    
        #This block evaluates the atop generalised coordination number for the whole system
        if bool(globals()['agcn']*globals()['nn']*globals()['adj']) is True:
            Agcn, NN = AGCN.agcn_generator(result_cache['adj'], NN = True)
            metadata['agcn'][int(i/Step)] = Agcn; metadata['nn'][int(i/Step)] = NN
        elif bool(globals['agcn']*globals()['adj']) is True:
            Agcn = AGCN.agcn_generator(result_cache['adj'])[0]
            metadata['agcn'][int(i/Step)] = Agcn
        elif bool(globals['nn']*globals()['adj']) is True:
            _,NN = AGCN.agcn_generator(result_cache['adj'], NN = True)
            metadata['nn'][int(i/Step)] = NN
    
        
        ##This is simply a progress updater which informs the user how every 5% is getting along.
        if i%int((End-Start)/20) == 0:
            Per = int(i/int((End-Start)/100))
            T1=time.time()
            print("Step %s computed in only %.4f nanocenturies." %(i, (T1-T0)/3.156), "\n")
            print("Currently performed %.3f%% of the calculation." %(Per), "\n")
            print('['+int(Per/5)*'##'+(20-int(Per/5))*'  '+']', "\n")
            if i >Start:
                print('Estimated time for completion is %s.' %(
                                                              time.strftime(
                                                                            "%H:%M:%S", 
                                                                             time.gmtime((100/Per)*(time.time()-BigT)))), "\n"                         
                                                                            )

        i += Step


    Masterkey.sort()
    metadata['masterkey']=Masterkey

    """
                
    #############################################################################################
    
    Robert:
        And now we check to see if the users wishes to evaluate any of the quantities
        from the energy file and add them to the metadata.
    
    """


    if bool(globals()['SimTime']) is True:
        metadata['SimTime'] = energy[:,0]

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
