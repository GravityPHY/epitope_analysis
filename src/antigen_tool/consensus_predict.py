import os
import numpy as np
import Bio.PDB.StructureBuilder
from Bio.PDB import PDBIO, PDBParser, PPBuilder, CaPPBuilder, Superimposer
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
from Bio import SeqIO, Align

import transform_utils


def majority_vote(structures, path):
    """
    The first structure in the structures list will be the reference structure
    Args:
        structures (list): list of PDB structures
        path (pathlib.Path): saving path
    Returns:
        None
    """

def borda_score(tuple):
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

def sum_tuple_score(tuples,index=-1):
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


def borda_count(structures, path=None):
    """
    Assume the structures only differ at bfactor column!!
    Args:
        structures (list): list of PDB structures
        path (pathlib.Path): saving path
    Returns:
        None
    """
    tuples=[transform_utils.get_sequence_bfactor_rank(structure) for structure in structures]
    tuple_scores=[borda_score(tuple) for tuple in tuples]
    summed_=sum_tuple_score(tuple_scores)
    return summed_



