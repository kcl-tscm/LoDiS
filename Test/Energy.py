import numpy as np
from Quantity import Quantity
class Time(Quantity):
    name = 'time'
    display_name = 'Time'
    def calculate(self, i_frame, result_cache, metadata):
        return i_frame * metadata['time_step']
    def get_dimensions(self, n_atoms, n_elements):
        return 1, np.float
    
    
class EnergyFile(Quantity):
    display_name = 'Energy File'
    def __init__(self, params, System, Settings, Metadata):
        f_name = System['base_dir'] + System['energy_file_name']
        _, self.pot_energy, self.tot_energy, self.kin_energy, _, _, self.temperature = np.loadtxt(f_name, unpack=True)
        if self.pot_energy.shape[0] != Metadata['n_frames']:
            raise ValueError('Energy file and movie file have different number of frames')
    def calculate(self, i_frame, result_cache, metadata):
        return [self.pot_energy[i_frame], self.tot_energy[i_frame],
                self.kin_energy[i_frame], self.temperature[i_frame]]
        
        
class PotentialEnergy(Quantity):
    name = 'e_pot'
    display_name = 'Potential Energy'
    dependencies = [EnergyFile]
    def calculate(self, i_frame, result_cache, Metadata):
        return result_cache[EnergyFile][0]
    def get_dimensions(self, n_atoms, n_elements):
        return 1, np.float
    
    
class TotalEnergy(Quantity):
    name = 'e_tot'
    display_name = 'Total Energy'
    dependencies = [EnergyFile]
    def calculate(self, i_frame, result_cache, Metadata):
        return result_cache[EnergyFile][1]
    def get_dimensions(self, n_atoms, n_elements):
        return 1, np.float
    
    
class KineticEnergy(Quantity):
    name = 'e_kin'
    display_name = 'Kinetic Energy'
    dependencies = [EnergyFile]
    def calculate(self, i_frame, result_cache, Metadata):
        return result_cache[EnergyFile][2]
    def get_dimensions(self, n_atoms, n_elements):
        return 1, np.float
    
    
class Temperature(Quantity):
    name = 'temp'
    display_name = 'Temperature'
    dependencies = [EnergyFile]
    def calculate(self, i_frame, result_cache, Metadata):
        return result_cache[EnergyFile][3]
    def get_dimensions(self, n_atoms, n_elements):
        return 1, np.float