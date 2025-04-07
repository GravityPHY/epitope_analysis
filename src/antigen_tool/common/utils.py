import os

from typing import Any, Dict, List, Mapping, Optional, Tuple

from Bio.Seq import Seq
from Bio.PDB.Structure import Structure
from Bio.PDB import PDBIO, PDBParser, MMCIFParser, PPBuilder
from Bio.PDB.Polypeptide import is_aa
from Bio.Data.SCOPData import protein_letters_3to1 as aa3to1

def get_structure_cif(cif_file:str)->Structure:
    parser = MMCIFParser()
    structure = parser.get_structure('cif', cif_file)
    return structure

def get_structure_pdb(pdb_file:str)->Structure:
    parser = PDBParser()
    structure = parser.get_structure('pdb', pdb_file)
    return structure

def structure_to_seqs(structure: Structure) -> Dict[str, Seq]:
    """
    Given a structure, extract the amino acid sequences of all chains.

    Args:
        structure (Structure): a Biopython structure

    Returns:
        Dict[str, Seq]: Dictionary mapping chain IDs to amino acid sequences
    """
    model = next(structure.get_models())  # Usually just 1 model
    ppb = PPBuilder()
    chain_seqs = {}

    for chain in model:
        peptides = ppb.build_peptides(chain, aa_only=False)
        if peptides:
            seq = peptides[0].get_sequence()
            chain_seqs[chain.id] = seq
        else:
            chain_seqs[chain.id] = Seq("")  # Empty chain

    return chain_seqs


def get_sequence_bfactor(structure:Structure,
                         sum:bool=False)->List[Tuple[str, int,float]]:
    """
    Retrieve the AA sequence from a structure, residue number and b-factor of the residue (CA atom)
    Args:
        structure (Structure): a Biopython structure
        sum (bool): True enable the summation bfactor values of a residue; False disable

    Returns:
        seq_resid_bfactor (list):
         list of (amino acid, res_id, bfactor) tuple.
         For example [('A',1,0.1),('G',2,0.8),('H',3, 0)]
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
            seq_resid_bfactor.append((_aainfo(r)[1],
                                      r.id[1],
                                      bfactor))
        else:
            seq_resid_bfactor = [(_aainfo(r)[1],
                                  r.id[1],
                                  next(r.get_atoms()).get_bfactor()) for r in structure.get_residues()
                           if is_aa(r)]
    return seq_resid_bfactor

def get_sequence_bfactor_rank(structure: Structure,
                              sum:bool=False,
                              tie_handling:str = "competition")->List[Tuple[str, int,float]]:
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
