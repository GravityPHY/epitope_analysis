import os

from typing import Any, Dict, List, Mapping, Optional, Tuple

import numpy as np

from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure

from ..common import utils


def majority_vote(structures:List[Structure], path):
    """
    The first structure in the structures list will be the reference structure
    Args:
        structures (list): list of PDB structures
        path (pathlib.Path): saving path
    Returns:
        None
    """

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

def sum_tuple_score(tuples:List[Tuple],index=-1):
    """
    Assume the tuple in tuples has the same length!
    Args:
        tuples (list): list of tuples
    Returns:
        summed_tuple (list): list of tuples
    """
    summed_tuple=[]
    for tuple in zip(*tuples):
        summed_score=sum(map(lambda x:x[index], tuple))
        summed_tuple.append((*tuple[0][0:-1],summed_score))
    return summed_tuple


def borda_count(structures:List[Structure], path=None):
    """
    Assume the structures only differ at bfactor column!!
    Args:
        structures (list): list of PDB structures
        path (pathlib.Path): saving path
    Returns:
        None
    """
    tuples=[utils.get_sequence_bfactor_rank(structure) for structure in structures]
    tuple_scores=[borda_score(tuple) for tuple in tuples]
    summed_=sum_tuple_score(tuple_scores)
    return summed_