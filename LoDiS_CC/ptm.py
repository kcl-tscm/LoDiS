from asap3.analysis import PTM
from ase.visualize.primiplotter import *
from ase.io import read

#filename = "Au2223-FCC-Dh-Ih_test.xyz"

#Above was used for script testing purposes.

""" (Armand)
PTM function. Given a filename and a rmsd_max cutoff,
uses the asap3 implementation of the PTM,
returns a dictionary with NumPy arrays as values.
of interest:
'structure': 0 = none; 1 = FCC; 2 = HCP; 3 = BCC; 4 = Icosahedral; 5 = SC
             0 = none; 1 = SC; 2 = FCC; 3 = HCP; 4 = Icosahedral; 5 = BCC

about these two contraditory lines above: the asap3 documentation gives the first line
as the results identified. Not believing it, I went to look at the source code,
which gives us the second line. Requires to be tested, its on the to do list.

'rmsd': the RMSD error in the fitting of the template. No values will be accepted above rmsd_max,
        and the atom would instead be considered as unidentifiable.
        
        When you return to this. Please give a reasonable value for rmsd_max_i 
        - Rob -

this function still has the bug of only returning the last frame's PTM. Needs fixing.
"""

def get_PTM_array(filename,rmsd_max_i):
    traj = read(filename, index=slice(None), format='xyz')
    for frame, atoms in enumerate(traj):
        atoms.set_cell([[100, 0, 0], [0, 100, 0], [0, 0, 100]])
        ptmdata = PTM(atoms, rmsd_max=rmsd_max_i)
    return(ptmdata)
