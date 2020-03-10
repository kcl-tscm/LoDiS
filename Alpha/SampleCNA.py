import Process
import pickle

Supported=[
        'rdf', 'cna', 'adj', 'pdf', 'pdfhomo', 'agcn', 'nn',
        'SimTime', 'EPot', 'ETot', 'EKin', 'EDelta', 'MeanETot', 'Temp'
           ]

System = {
        'base_dir' : 'TestTraj/',
        'movie_file_name' : 'movie.xyz',
        'energy_file_name' : 'energy.out',
        
        'Homo' : ['Au', 'Pd'], 'HomoQuants' : [ 'HoPDF', 'HoRDF', 'CoM', 'HoAdj' ], 
        'Hetero' : True, 'HeteroQuants' : [ 'HePDF', 'HeRDF', 'HeAdj' ],
        
        'Start' : 0, 'End' : 1000, 'Step' : 1, 'Skip' : 50, 'UniformPDF' : False, 'Band' : 0.05,
        
        'HCStats' : True,
        
        'SimTime': True, 'EPot': True, 'ETot' : True, 
        'EKin' : True, 'EDelta' : True, 'MeanETot' : True, 'Temp' : True
        }

Quantities = {
        'euc' : None, 'rdf' : None, 'pos' : None, 'cna' : None, 
        'adj' : None, 'pdf' : None, 'agcn' : None, 'nn' : None, 'CoM' : None,
        'SimTime': None, 'EPot': None, 'ETot' : None, 
        'EKin' : None, 'EDelta' : None, 'MeanETot' : None, 'Temp' : None
        }

TestMetadata = Process.Process(System, Quantities)

with open('TestTraj/MetadataSample.csv', "wb") as f:
    pickle.dump(TestMetadata,f, pickle.HIGHEST_PROTOCOL)
