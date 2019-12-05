import numpy as np
from ase.io import read
from asap3.analysis import FullCNA


def transform_cna(cna, meaningful_cnas):
    """ (Claudio)
    Given a cna and a list containing the tuples defining the
    meaningful cnas, this method returns an array with length equal to the
    total number of cnas which contains the number of times each cna occurrs.
    """
    result = np.zeros(len(meaningful_cnas))
    for i, key in enumerate(meaningful_cnas):
        try:
            result[i] = cna[key]
        except KeyError:
            result[i] = 0
    return result

def get_atomic_cnas(traj, meaningful_cnas, r_cut):
    """ (Claudio)
    Given a trajectory, a list of tuples containing the CNAS the user is interested
    in, and a cutoff radius, this returns, for every atom for every snapshot in traj,
    an array which contains the count of cnas for that atom, arranged according to the order
    found in meaningful_cnas
    """
    transformed_cna = np.zeros((len(traj)*len(traj[0]), len(meaningful_cnas)))
    for j, atoms in enumerate(traj):
        cna = FullCNA(atoms, r_cut)
        atoms.set_cell([[100,0,0],[0,100,0],[0,0,100]])
        snapshot_cna = cna.get_normal_cna()
        for i, atomic_cna in enumerate(snapshot_cna):
            transformed_cna[j*len(traj[0])+i] = transform_cna(atomic_cna, meaningful_cnas)
    return transformed_cna

def get_all_cnas(traj, r_cut):
    """ (Claudio)
    Given a trajectory and a cutoff radius, returns a dictionary, sorted by value,
    with all the cnas that appear in the trajectory as keys and the number
    of times they appear as value.
    """
    all_cnas = {}
    for j, atoms in enumerate(traj):
        cna = FullCNA(atoms, r_cut)
        atoms.set_cell([[100,0,0],[0,100,0],[0,0,100]])
        snapshot_cna = cna.get_normal_cna()
        for i, atomic_cna in enumerate(snapshot_cna):
            for key in atomic_cna:
                try:
                    all_cnas[key] += atomic_cna[key]
                except KeyError:
                    all_cnas[key] = atomic_cna[key]

    sorted_cnas = sorted(all_cnas.items(), key=lambda kv: -kv[1])
    sorted_cnas_dict = {}
    for t in sorted_cnas:
        sorted_cnas_dict[t[0]] = t[1]

    return sorted_cnas_dict

def extract_cnas(traj, r_cut):
    """ (Claudio)
    Get all the cnas in the trajectory file, then extract the atomic CNA signatures for each CNA present.
    For each atom, the atomic cnas contains a row entry with dimensionality equal to the number of cnas present in the trajectory.
    The all_cnas variable contains a count of occurrance of each cna in the trajectory.
    """
    all_cnas = get_all_cnas(traj, r_cut)
    atomic_cnas = get_atomic_cnas(traj, all_cnas, r_cut)
    return atomic_cnas, all_cnas

def sample_uniform_cna(ntr, transformed_cnas):
    """ (Claudio)
    Sample from an array of transformed cnas a ntr number of indexes.
    For each cnas class, ntr//len(cna classes) atoms are selected which do
    contain at least one pair of that particular class.
    """
    tr_ind = []
    sampled_atoms = np.ones(len(transformed_cnas), dtype = 'bool')
    ntr = 500
    ntr_sampled = 0
    for i in range(transformed_cnas.shape[1]):
        indx_this_class = np.where(transformed_cnas[:,i][sampled_atoms] > 0)[0]
        ntr_this_class = min(len(indx_this_class), ntr//transformed_cnas.shape[1])
        sampled_inds = np.random.choice(indx_this_class, ntr_this_class, replace = False)
        sampled_atoms[sampled_inds] = False
        tr_ind.extend(sampled_inds)
        ntr_sampled += len(sampled_inds)

    if ntr_sampled < ntr:
        additional_inds = np.random.choice(np.arange(len(transformed_cnas))[sampled_atoms], ntr-ntr_sampled, replace = False)
        tr_ind.extend(additional_inds)
    return np.array(tr_ind)

def sample_cna(traj, r_cut, ntr, metric = "uniform"):
    """ (Claudio)
    From a trajectory file, calculate CNAS using r_cut as cutoff,
    order the classes and sample according to the sample_using_cna method
    """
    transformed_cnas, all_cnas = extract_cnas(traj,r_cut)
    print("CNA classes are: \n", all_cnas)
    if metric == "uniform":
        training_indexes = sample_uniform_cna(ntr, transformed_cnas)
    return training_indexes

#END OF DEFINITIONS
