import os
from typing import Any, Dict, List, Mapping, Optional, Tuple

import numpy as np

from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure

from featurization import get_seq_resid_rank

class Combine:
    def __init__(self,structures:List[Structure]):
        self.structures=structures

def borda_score(tuple:List[Tuple])->List[Tuple]:
    """
    Assume the last element of the tuple is the rank
    Args:
        tuple (list): list of tuples
    Returns:
        tuple (list): list of tuples
    """
    total_=len(tuple)
    new_tuples=[]
    for t in tuple:
        new_tuples.append((*t[0:-2],total_-t[-1]))
        #t[-1]=total_-t[-1]
    return new_tuples

def borda_count(structures:List[Structure], path=None):
    """
    Assume the structures only differ at bfactor column!!
    Args:
        structures (list): list of PDB structures
        path (pathlib.Path): saving path
    Returns:
        None
    """
    tuples=[get_seq_resid_rank(struct) for struct in structures]
    tuple_scores=[borda_score(tuple) for tuple in tuples]
    summed_=sum_tuple_score(tuple_scores)
    return summed_