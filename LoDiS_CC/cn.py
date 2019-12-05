from asap3 import FullNeighborList
from asap3.analysis import CoordinationNumbers
from ase import Atoms
from ase.io import read
import numpy as np

cnmax = 14
"""(Armand)
ACTUAL CNMAX FOUND IN THIS PARTICLE WAS 14, 12 was used previously,
not sure why.
"""
def get_coord_number(filename, r_cut_cn):

	traj = read(filename, index = ':', format='xyz')

	#################ELEMENT COMPOSITION CHECKER#################
	#############################################################
	atom_number_list = [atoms.get_atomic_numbers() for atoms in traj]
	flat_atom_number = np.concatenate(atom_number_list)
	elements = np.unique(flat_atom_number, return_counts=False)
	#checks the elemental composition, displays them in elements#

	for frame, atoms in enumerate(traj):
		atoms.set_cell([[100, 0, 0], [0, 100, 0], [0, 0, 100]])
		cn = []
		gcn =[]
		#atgcn = []
		nl = FullNeighborList(r_cut_cn, atoms = atoms)
		for atom in np.arange(len(atoms)):
			distance_atom = []
			indices, positions, dist = nl.get_neighbors(atom)
			distance_atom.extend(np.sqrt(dist))
			cn_atom = len(distance_atom)
			cn.append(cn_atom)
			gcn.append(sum(cn)/cnmax)
	return cn, gcn
"""(Matteo)
atop gcn: matrix with 0 and 1 (neigh. or no neigh.) times cn
bridge gcn: matrix with i+j. if i+j=2 convert to 1. Then multiply by cn. Use cosine function to normalize.
in pdf: use gaussian convolution.
"""


"""(Armand)
Running a cn using asap3 to compare results
"""
def get_coord_number_asap(filename, r_cut_cn):
	traj= read(filename, index = ':', format='xyz')
	for frame, atoms in enumerate(traj):
		atoms.set_cell([[100, 0, 0], [0, 100, 0], [0, 0, 100]])
		cn_atom = CoordinationNumbers(atoms,r_cut_cn)
	return cn_atom
