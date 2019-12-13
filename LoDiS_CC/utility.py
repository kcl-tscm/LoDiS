""" Robert

Creates a dictionary of possible classes to be calculated or iterated over

Items in "quantity_list" will be the class name for a given process
"""
from cn import *
from cna import *
from Distances import *
from ptm import *
from Kernels import *



quantity_list = [Euc_Dist, M1_M1, M1_M2, get_cutoff_distance,
                 Gauss, Uniform, Epan, KB_Dist, transform_cna,
                 get_atomic_cnas, get_all_cnas, extract_cnas,
                 sample_uniform_cna, sample_cna, get_coord_number,
                 get_coord_number_asap, cn_generator, get_PTM_array,]
                 

quantity_map = {q.name: q for q in quantity_list}

reverve_quantity_map = {v: k for k, v in quantity_map.items()}


def tuplize(q):
    return q if type(q) is tuple else (q, None)


def ordered_set(li: list):
    seen = []
    removed = []
    for el in li:
        if el not in seen:
            removed.append(el)
            seen.append(el)
    return removed


def input2key(inp):
    inp = tuplize(inp)
    return quantity_map[inp[0]], inp[1]


def key2input(key):
    key = reverve_quantity_map[key[0]], key[1]
    if key[1] is None:
        return key[0]
    else:
        return key
