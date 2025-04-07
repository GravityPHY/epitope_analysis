import os
import csv
import argparse
import itertools

from Bio import SeqIO, Align
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.Seq import Seq
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB import PDBIO, PDBParser, MMCIFParser, PPBuilder, CaPPBuilder
from Bio.Data.SCOPData import protein_letters_3to1 as aa3to1

import warnings
warnings.filterwarnings('ignore')

def get_structure_cif(cif_file):
    parser = MMCIFParser()
    structure = parser.get_structure('cif', cif_file)
    return structure


def get_structure_pdb(pdb_file):
    parser = PDBParser()
    structure = parser.get_structure('PDB', pdb_file)
    return structure


def get_Seq(structure):
    """
    Retrieve the AA sequence from a PDB structure as a Seq
    Args:
        structure (PDB_structure): PDB structure
    Returns:
        Seq object defined in Biopython
    """
    _aainfo = lambda r: (r.id[1], aa3to1.get(r.resname, 'X'))
    seq = [_aainfo(r)[1] for r in structure.get_residues() if is_aa(r)]
    return Seq(''.join(seq))


def get_sequence_bfactor(structure, sum=False):
    """
    Retrieve the AA sequence from a structure, residue number and b-factor of the residue (CA atom)
    Args:
        structure (): PDB structure
        sum (bool): True enable the summation; False disable

    Returns:
        seq_resid_bfactor (list):
         list of (amino acid, res_id, bfactor) tuple. For example [('A',1,0.1),('G',2,0.8),('H',3, 0)]
    """
    _aainfo = lambda r: (r.id[1], aa3to1.get(r.resname, 'X'))
    seq_resid_bfactor = []
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
            seq_resid_bfactor.append((_aainfo(r)[1],
                                      r.id[1],
                                      bfactor))
        else:
            seq_resid_bfactor = [(_aainfo(r)[1],
                                  r.id[1],
                                  next(r.get_atoms()).get_bfactor()) for r in structure.get_residues()
                           if is_aa(r)]
    return seq_resid_bfactor


def get_ch_sequence_bfactor(structure, sum=False):
    """
    Retrieve the AA sequence from a structure and the b-factor of the residue (CA atom)
    :param structure:
    :param sum: the option of summing all atomic b-factor values in a residue. True enable the summation, False disabled.
    :return: dict of (amino acid, bfactor) tuple for each chain . For example {"A":[('A',0.1),('G',0.8),('H',0)]}
    """
    chains = [ch.id for ch in structure.get_chains()]
    _aainfo = lambda r: (r.id[1], aa3to1.get(r.resname, 'X'))
    ch_seq_bfactor = {}
    for ch in structure.get_chains():
        ch_seq_bfactor[ch.id] = []
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
                ch_seq_bfactor[ch.id] = [(_aainfo(r)[1], next(r.get_atoms()).get_bfactor()) for r in
                                         structure.get_residues()
                                         if is_aa(r)]
    return combine_tuples(ch_seq_bfactor)


def combine_tuples(data):
    assert isinstance(data, dict)
    combined = []
    for key in sorted(data.keys(), key=lambda x: tuple(sorted(x))):
        combined.extend(data[key])
    return combined


def get_sequence_bfactor_rank(structure, sum=False, tie_handling = "competition"):
    """
    Retrieve the AA sequence from a structure,
    residue number and b-factor of the residue (CA atom),
    ranking for b-factor values
    Args:
        structure (): PDB structure
        sum (bool): True enable the summation; False disable

    Returns:
        seq_resid_bfactor (list):
         list of (amino acid, res_id, bfactor) tuple. For example [('A',1,0.1),('G',2,0.8),('H',3, 0)]
    """
    seq_resid_bfactor = get_sequence_bfactor(structure, sum)
    indexed_values = [(tup[-1], i) for i, tup in enumerate(seq_resid_bfactor)]
    sorted_values = sorted(indexed_values,reverse=True)
    rank_map={}
    ranked_tuples = [None]*len(seq_resid_bfactor)
    if tie_handling == "competition":
        rank = 1
        for i, (value, original_index) in enumerate(sorted_values):
            if value not in rank_map:
                rank_map[value] = rank
            ranked_tuples[original_index] = seq_resid_bfactor[original_index] + (rank_map[value],)
            rank = i +2
    elif tie_handling == "dense":
        rank = 1
        prev_value = None
        for value, original_index in sorted_values:
            if value != prev_value:
                rank_map[value] = rank
                rank += 1
                prev_value = value
            ranked_tuples[original_index] = seq_resid_bfactor[original_index] + (rank_map[value],)
    else:
        raise ValueError("Invalid tie_handling method.")
    return ranked_tuples



def pdb_to_seq(pdb_file):
    """
    Give a path of .pdb with single chain, parse structure to get the sequence
    Args:
        param pdb_file (str): path of pdb file, should end with .pdb
    Returns:
        PDB_chain_seq (Seq): amino acid sequence in Seq object
    """

    PDB_structure = get_structure_pdb(pdb_file)
    for model in PDB_structure:
        for chain_name in model:
            PDB_chain = PDB_structure[0][chain_name.id]
    PDB_chain_seq = PPBuilder().build_peptides(PDB_chain, aa_only=False)[0].get_sequence()
    return PDB_chain_seq


def align_sequence_bfactor(aligned_seq_list, seq_bfactor_list):
    """
    update the seq_bfactor list with gap residue (represent in '0')
    Args:
    aligned_seq_list (list):
    seq_bfactor_list (list):
    Returns:
    None
    """
    index = 0
    for token in aligned_seq_list:
        if token == '-':
            seq_bfactor_list.insert(index, ('0', None,0.0))
            index += 1
        else:
            index += 1


def align_seq_bfactor(structure_pred, structure_true, sum):
    """
    Align the unbound sequence with bound sequence,
    and return the aligned (res, res_id, bfactor value) list
    Args:
        structure_pred (): structure with predict epitope labels
        structure_true (): structure with true epitope labels
        sum (bool): option to sum the atom bfactor or not. True - sum, False get the residue bfactor
    Returns:
        both the prediction and true lists of tuples (res, res_id, bfactor value)

    """
    seq_pred = get_Seq(structure_pred)
    bfactor_pred = get_ch_sequence_bfactor(structure_pred, sum=sum)
    seq_true = get_Seq(structure_true)
    bfactor_true = get_ch_sequence_bfactor(structure_true)
    #aligner = Align.PairwiseAligner()
    #alignments = aligner.align(seq_pred, seq_true)
    alignments = pairwise2.align.globalxx(str(seq_pred), str(seq_true))
    aligned_query, aligned_target, score, begin, end = alignments[0]
    align_sequence_bfactor(aligned_query, bfactor_pred)
    align_sequence_bfactor(aligned_target, bfactor_true)
    #align_sequence_bfactor(alignments[0].query, bfactor_pred)
    #align_sequence_bfactor(alignments[0].target, bfactor_true)
    return bfactor_pred, bfactor_true


def sum_atom_bfactor_residue(structure, pdb_file):
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


def map_to_structure(pdb_file, tuple,save_path):
    """
    Given the tuple with new bfactor,
    tuple in the format (res_name,res_id, bfactor)
    add to structure bfactor columns
    Args:
        structure (Structure):
        tuple (list):
    Returns:
        None
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", pdb_file)
    for model in structure:
        for chain in model:
            for residue,t in zip(chain,tuple):
                if residue.id[1] == t[1]:
                    for atom in residue:
                        atom.set_bfactor(t[2])
    io = PDBIO()
    io.set_structure(structure)
    io.save(save_path)
    print(f"updated: {save_path}")

