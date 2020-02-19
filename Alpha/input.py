from ExecuteTest import *
import pandas as pd


def Meta_Csv(Data, filename):
    df = pd.DataFrame.from_dict(Data, orient = "index")
    df.to_csv(System['base_dir']+filename+'.csv')

"""
Robert:
    Change the arguments in the System dictionary to 
    reflect the directory and files you are working with.
"""

"""

Robert:
    
    This list below is just to inform you as to what quantities you can calculate.
    
"""

Supported=[
        'euc', 'rdf', 'cna', 'adj', 'pdf', 'pdfhomo', 'agcn', 'nn', 'pos',
        'SimTime', 'EPot', 'ETot', 'EKin', 'EDelta', 'MeanETot', 'Temp'
           ]

System = {
        'base_dir' : '../../../Pd/1127/Melting100/',
        'movie_file_name' : 'PdMelt5-Movie.xyz',
        'energy_file_name' : 'PdMelt5-En.out',
        
        'Homo' : 'Cu', # Don't worry about this guy just yet, it will become relevant soon.
        
        'Start' : 0, 'End' : 100, 'Step' : 1, 'Skip' : 50, 'UniformPDF' : False, 'Band' : 0.05,
        'PdfStats' : True, 'HomoStats' : True, 'RdfStats' : True, 'CnaStats' : True,
        
        #Below are quantities related to the energy file
        
        'SimTime': True, 'EPot': True, 'ETot' : True, 
        'EKin' : True, 'EDelta' : True, 'MeanETot' : True, 'Temp' : True
        }



Quantities = {
        'euc' : None,
        'rdf' : None, 
        'pos' : None, 
        'cna' : None, 
        'adj' : None, 
        'pdf' : None,
        #'pdfhomo' : np.empty((Time,), dtype=object),
        'agcn' : None,
        'nn' : None,
        'SimTime': None, 'EPot': None, 'ETot' : None, 
        'EKin' : None, 'EDelta' : None, 'MeanETot' : None, 'Temp' : None
        }

TestMetadata = Process(System, Quantities)


Meta_Csv(TestMetadata, 'TestingMeta')
