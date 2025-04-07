import os
import csv
import argparse
import itertools

from typing import Any, Dict, List, Mapping, Optional, Tuple

from Bio import SeqIO, Align, pairwise2
from Bio.Seq import Seq
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB import PDBIO, PDBParser, MMCIFParser, PPBuilder, CaPPBuilder
from Bio.Data.SCOPData import protein_letters_3to1 as aa3to1


import warnings
warnings.filterwarnings('ignore')


def get_structure_cif(cif_file:str)->Structure:
    parser = MMCIFParser()
    structure = parser.get_structure('cif', cif_file)
    return structure


def get_structure_pdb(pdb_file:str)->Structure:
    parser = PDBParser()
    structure = parser.get_structure('PDB', pdb_file)
    return structure


def pdb_to_seq(pdb_file:str)-> Chain:
    """
    Give a path of .pdb with single chain, parse structure to get the sequence
    Args:
        pdb_file (str): path of pdb file, should end with .pdb
    Returns:
        PDB_chain_seq (Bio.Seq.Seq): amino acid sequence
    """
    parser = PDBParser()
    PDB_structure = parser.get_structure('PDB', pdb_file)
    for model in PDB_structure:
        for chain_name in model:
            PDB_chain = PDB_structure[0][chain_name.id]
    PDB_chain_seq = PPBuilder().build_peptides(PDB_chain, aa_only=False)[0].get_sequence()
    return PDB_chain_seq


def get_Seq(structure: Structure)->Seq:
    """
    Retrieve the AA sequence from a PDB structure as a Seq
    Args:
        structure (Bio.PDB.Structure.Structure): PDB structure
    Returns:
        PDB_chain_seq (Bio.Seq.Seq): amino acid sequence
    """
    _aainfo = lambda r: (r.id[1], aa3to1.get(r.resname, 'X'))
    seq = [_aainfo(r)[1] for r in structure.get_residues() if is_aa(r)]
    return Seq(''.join(seq))


def get_sequence_bfactor(structure: Structure, sum=False)->list:
    """
    Retrieve the AA sequence from a structure and the b-factor of the residue (CA atom)
    Args:
        structure (Bio.PDB.Structure.Structure): PDB structure
        sum (bool): the option of summing all b-factor values in a residue. True enable the summation, False disabled.
    Returns:
        seq_bfactor (list): list of (amino acid, bfactor) tuple. For example [('A',0.1),('G',0.8),('H',0)]
    """
    _aainfo = lambda r: (r.id[1], aa3to1.get(r.resname, 'X'))
    seq_bfactor = []
    for r in structure.get_residues():
        bfactor = 0
        counter = 0
        if sum:
            for atom in r.get_atoms():
                bfactor += atom.get_bfactor()
                counter += 1
            # if atom.get_id() == 'CA':
            #    bfactor = atom.get_bfactor()
            # seq_bfactor = [(_aainfo(r)[1], next(r.get_atoms()).get_bfactor()) for r in structure.get_residues() if is_aa(r)]
            seq_bfactor.append((_aainfo(r)[1], bfactor/counter))
        else:
            seq_bfactor = [(_aainfo(r)[1], next(r.get_atoms()).get_bfactor()) for r in structure.get_residues()
                           if is_aa(r)]
    return seq_bfactor


def get_ch_sequence_bfactor(structure: Structure, sum=False):

    # TODO: improve this function, with full residue name
    """
    Retrieve the AA sequence from a structure and the b-factor of the residue (CA atom)
    :param structure:
    :param sum: the option of summing all b-factor values in a residue. True enable the summation, False disabled.
    :return: dict of (amino acid, bfactor) tuple for each chain . For example {"A":[('A',0.1),('G',0.8),('H',0)]}
    """
    chains=[ch.id for ch in structure.get_chains()]
    _aainfo = lambda r: (r.id[1], aa3to1.get(r.resname, 'X'))
    ch_seq_bfactor={}
    for ch in structure.get_chains():
        ch_seq_bfactor[ch.id]= []
        for r in structure[0][ch.id].get_residues():
            bfactor = 0
            counter = 0
            if sum:
                for atom in r.get_atoms():
                    bfactor += atom.get_bfactor()
                    counter += 1
                # if atom.get_id() == 'CA':
                #    bfactor = atom.get_bfactor()
                # seq_bfactor = [(_aainfo(r)[1], next(r.get_atoms()).get_bfactor()) for r in structure.get_residues() if is_aa(r)]
                ch_seq_bfactor[ch.id].append((_aainfo(r)[1], bfactor))
            else:
                ch_seq_bfactor[ch.id] = [(_aainfo(r)[1], next(r.get_atoms()).get_bfactor()) for r in structure.get_residues()
                           if is_aa(r)]
    return combine_tuples(ch_seq_bfactor)


def combine_tuples(data:dict)-> list:
    assert isinstance(data,dict)
    combined=[]
    for key in sorted(data.keys(), key = lambda x: tuple(sorted(x))):
        combined.extend(data[key])
    return combined

def sum_atom_bfactor_residue(structure:Structure,
                             pdb_file:str)->None:
    """
    Only apply to AbEMap pdb at the moment, will save the structure in .pdb to pdb_file
    :param structure:
    :return:None
    """
    for r in structure.get_residues():
        bfactor = 0
        for atom in r.get_atoms():
            bfactor += atom.get_bfactor()
        for atom in r.get_atoms():
            atom.set_bfactor(bfactor)
    io = PDBIO()
    io.set_structure(structure)
    io.save(pdb_file)
## in transform_utils.py
##
def align_sequence_bfactor(aligned_seq_list:list,
                           seq_bfactor_list:list)->None:
    """
    update the seq_bfactor list with gap residue (represent in '0')
    :param aligned_seq_list:
    :param seq_bfactor_list: [('res', bfactor value)]
    :return:
    """
    index = 0
    for token in aligned_seq_list:
        if token == '-':
            seq_bfactor_list.insert(index, ('0', 0.0))
            index += 1
        else:
            index += 1


def align_seq_bfactor(structure_pred: Structure,
                      structure_true: Structure,
                      sum:bool):
    """
    Align the unbound sequence with bound sequence, and return the aligned (res, bfactor value) list
    :param structure_pred: structure with predict epitope labels
    :param structure_true: structure with true epitope labels
    :param sum: option to sum the atom bfactor or not. True - sum, False get the residue bfactor
    :return: both the prediction and true (res, bfactor value) lists
    """
    seq_pred = get_Seq(structure_pred)
    bfactor_pred = get_ch_sequence_bfactor(structure_pred, sum=sum)
    seq_true = get_Seq(structure_true)
    bfactor_true = get_ch_sequence_bfactor(structure_true)
    #aligner = Align.PairwiseAligner()
    #alignments = aligner.align(seq_pred, seq_true)
    #align_sequence_bfactor(alignments[0].query, bfactor_pred)
    #align_sequence_bfactor(alignments[0].target, bfactor_true)
    alignments = pairwise2.align.globalxx(seq_pred, seq_true)

    alignment = alignments[0]
    #print(alignment.seqA)
    align_sequence_bfactor(alignment.seqA, bfactor_pred)
    align_sequence_bfactor(alignment.seqB, bfactor_true)
    return bfactor_pred, bfactor_true


def combine_bfactor(bfactor_true, bfactor_pred1, bfactor_pred2):
    new = []
    for i, j, k in zip(bfactor_true, bfactor_pred1, bfactor_pred2):
        new.append((i[0], max(j[1], k[1])))
    return new