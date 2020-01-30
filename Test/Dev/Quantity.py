class Quantity():
    name = ''
    display_name = ''
    dependencies = ''
    number = ''
    
    def __init__(self, params, System, Settings, Metadata):
        pass
    
    def not_implemented_string(self):
        return "Sorry, pal! We don't support " + self.name
    
    def calculate(self, i_frame, results_cache, Metadata):
        raise NotImplementedError(self.not_implemented_string())
        
    def get_dimensions(self, NAtoms, NSpecies):
        raise NotImplementedError(self.not_implemented_string())
        
    def cleanup(self):
        pass
    
class ChemQuantity(Quantity):
    def __init__(self, params, System, Settings, Metadata):
        if params is None:
            self.n_atoms_a = Metadata['NAtoms']  # Avalible
            self.filter = None
        else:
            self.n_atoms_a = Metadata['NAtoms'][params]
            self.filter = Metadata['Species'][params]
            self.display_name += ' (' + params + ')'