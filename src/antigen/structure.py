from Bio.PDB import MMCIFParser
from Bio.PDB import PDBParser
from Bio.PDB import PDBIO
from Bio.Seq import Seq
from Bio.PDB.Chain import Chain
from Bio.PDB.Structure import Structure

import numpy as np
aa3to1={"ASH": "A",
    "ALA": "A",
    "CYX": "C",
    "CYS": "C",
    "ASP": "D",
    "GLU": "E",
    "PHE": "F",
    "GLY": "G",
    "HIS": "H",
    "HID": "H",
    "HIE": "H",
    "HIP": "H",
    "ILE": "I",
    "LYS": "K",
    "LEU": "L",
    "MET": "M",
    "MSE": "M",
    "ASN": "N",
    "PYL": "O",
    "HYP": "P",
    "PRO": "P",
    "GLN": "Q",
    "ARG": "R",
    "SER": "S",
    "THR": "T",
    "SEL": "U",
    "VAL": "V",
    "TRP": "W",
    "TYR": "Y",}

def get_structure_cif(cif_file:str,
                      name:str=None)->Structure:
    parser = MMCIFParser()
    if not name:
        structure = parser.get_structure('cif', cif_file)
    else:
        structure = parser.get_structure(name, cif_file)
    return structure

def get_structure_pdb(pdb_file:str,
                      name:str=None)->Structure:
    parser = PDBParser()
    if not name:
        structure = parser.get_structure('pdb', pdb_file)
    else:
        structure = parser.get_structure(name, pdb_file)
    return structure

def set_structure_color(pdb_file:str,
                        name:str=None,
                        output_file:str=None)->Structure:
    parser = PDBParser()
    if not name:
        structure = parser.get_structure('pdb', pdb_file)
    else:
        structure = parser.get_structure(name, pdb_file)
    new_structure=structure.copy()
    for model in new_structure:
        for chain in model:
            for residue in chain:
                for atom in residue.get_atoms():
                        atom.bfactor = np.mean([a.bfactor for a in residue.get_atoms()])
    io = PDBIO()
    io.set_structure(new_structure)
    io.save(output_file)
    return new_structure

def chain_to_sequence(chain: Chain) -> Seq:
    """
    Extract the amino acid sequence from a Biopython Chain object.

    Args:
        chain (Chain): Biopython Chain object.

    Returns:
        Seq: Amino acid sequence of the chain.
    """
    #ppb = PPBuilder()
    #peptides = ppb.build_peptides(chain, aa_only=False)
    _aainfo = lambda r: (r.id[1], aa3to1.get(r.resname, 'X'))
    seq = [_aainfo(r)[1] for r in chain.get_residues() if is_aa(r)]
    return Seq(''.join(seq))

