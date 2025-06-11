import os
from typing import Any, Dict, List, Mapping, Optional, Tuple

from Bio.Seq import Seq
from Bio.PDB.Structure import Structure

from Bio.PDB.SASA import ShrakeRupley
from Bio.SeqUtils.ProtParamData import kd  # Kyte–Doolittle scale

from structure import aa3to1


def structure_to_seqs(structure: Structure) -> Dict[str, Seq]:
    """
    Given a structure, extract the amino acid sequences from all chains in the structure.

    Args:
        structure (Structure): a Biopython structure

    Returns:
        Dict[str, Seq]: Dictionary mapping chain IDs to amino acid sequences
    """
    model = next(structure.get_models())  # Usually just 1 model
    chain_seqs = {}

    for chain in model:
        _aainfo = lambda r: (r.id[1], aa3to1.get(r.resname, 'X'))
        seq=""
        for r in chain.get_residues():
            seq+=aa3to1.get(r.resname, 'X')
        chain_seqs[chain.id] = Seq(seq) 
    return chain_seqs

def get_seq_resid_bfactor(structure: Structure,
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
        if sum:
            for atom in r.get_atoms():
                bfactor += atom.get_bfactor()
            seq_resid_bfactor.append((_aainfo(r)[1],
                                      r.id[1],
                                      bfactor))
        else:
            seq_resid_bfactor = [(_aainfo(r)[1],
                                  r.id[1],
                                  next(r.get_atoms()).get_bfactor()) for r in structure.get_residues()]
    return seq_resid_bfactor

def get_seq_resid_rank(structure: Structure,
                              sum:bool=False,
                              reversed:bool=False,
                              tie_handling:str = "competition")->List[Tuple[str, int,float]]:
    seq_resid_bfactor = get_seq_resid_bfactor(structure, sum)
    bfactor_values = [(tup[-1],i) for i, tup in enumerate(seq_resid_bfactor)]
    sorted_values = sorted(bfactor_values,reverse=reversed)
    rank_map={}
    ranked_tuples = [None]*len(seq_resid_bfactor)
    if tie_handling == "competition":
        rank = 1
        for i, (value, original_index) in enumerate(sorted_values):
            if value not in rank_map:
                rank_map[value]=rank
            ranked_tuples[original_index] = seq_resid_bfactor[original_index] + (rank_map[value],)
            rank = i + 2
    elif tie_handling == "dense":
        rank = 1
        prev_value = None
        for value, original_index in sorted_values:
            if value != prev_value:
                rank_map[value] = rank
                rank +=1
                prev_value = value
            ranked_tuples[original_index] = seq_resid_bfactor[original_index] + (rank_map[value],)
    else:
        raise ValueError("Invalid tie handling method.")
    return ranked_tuples

    
    

# hydrophobicity
def get_seq_resid_hyrophobicity(structure: Structure):
    ## TODO: hyrophobicity so far is zero
    sr = ShrakeRupley()
    sr.compute(structure, level='R')
    seq_resid_hydro = []
    for model in structure:
        for chain in model:
            for residue in chain:
                r=residue.get_resname()
                hydrophobicity = kd.get(r, 0.0)
                seq_resid_hydro.append((r,
                                        residue.id[1],
                                        hydrophobicity))
    return seq_resid_hydro
# solvent accessible surface area
def get_seq_resid_sasa(structure: Structure,
                       radius:float=5.0):
    sr = ShrakeRupley(probe_radius=radius)
    sr.compute(structure, level='R')
    seq_resid_sasa = []
    for model in structure:
        for chain in model:
            for residue in chain:
                seq_resid_sasa.append((residue.get_resname(),
                                          residue.id[1],
                                          residue.sasa))
    return seq_resid_sasa

# interface topology