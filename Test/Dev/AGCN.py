
"""Function that reads an .xyz file and returns a list or the 
coordination numbers (CNs) corresponding to each atoms, and their
atop generalized coordination numbers (aGCNs). 
​
The length of the two lists is equal to the number of atoms in the 
cluster. 
​
The CN is the number of nearest neighbours (neighbours within 
a distance of r_cut_cn, calculated through the PDF) of each atom.
​
The aGCN is the sum of all the CNs of the neighboring atoms, 
divided by the CN of the atoms in the bulk, which is equal to 12 
for the atop sites."""
​
def cn_generator (positions, r_cut_cn):
    """Creating two empty lists"""
    cn=[]
    agcn=[]
    for i, atoms in enumerate(positions):
        atoms.set_cell([[100, 0, 0], [0, 100, 0], [0, 0, 100]])
        """Creating an empty list for each atom, to which the 
        indices of all its nearest neighbours will be appended"""
        ind=[]
        for j in np.arange(len(atoms)):
            nl = FullNeighborList(r_cut_cn, atoms=atoms)
            indices, positions, distances = nl.get_neighbors(j)
            ind.append([int(k) for k in indices])
            distancej=[]    
            distancej.extend(distances**0.5)
            """The CN of each atom corresponds to the length of the 
            distance array calculated with asap3"""
            cnj=len(distancej)
            cn.append(cnj)
​
        for l in np.arange(len(atoms)):
            """List of all the indices of the neighbours of each atom"""
            cc=ind[l][:]
            list=[]
            """Appending the CN of each of these neighbors to a list and 
            calculating their sum divided by 12, giving the atop GCN"""
            for m in range(len(cc)):
                list.append(cn[ind[l][m]])
                sm=sum(list)/12
            agcn.append(sm)
    return (cn, agcn)