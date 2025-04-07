import os

from typing import Any, Dict, List, Mapping, Optional, Tuple

from Bio import SeqIO, Align, pairwise2
from Bio.Seq import Seq
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure
from Bio.PDB.Polypeptide import PPBuilder
from Bio.PDB.Polypeptide import is_aa
from Bio.Data.SCOPData import protein_letters_3to1 as aa3to1
from .utils import get_structure_pdb


def extract_chain_sequence(chain: Chain) -> Seq:
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


def get_ordered_residues(chain: Chain) -> List:
    """
    Get residues in a chain, excluding heteroatoms and waters.

    Args:
        chain (Chain): Biopython Chain object.

    Returns:
        List[Residue]: Ordered list of standard amino acid residues.
    """

    return [res for res in chain.get_residues() if res.id[0] == ' ']


def get_bfactor(residue, use_ca_only: bool) -> Optional[float]:
    """
    Get the B-factor from a residue.

    Args:
        residue (Residue): Biopython residue object.
        use_ca_only (bool): Whether to use only the CA atom or average all atoms.

    Returns:
        float or None: The B-factor value or None if CA atom is missing.
    """
    if use_ca_only:
        if 'CA' in residue:
            return residue['CA'].get_bfactor()
        else:
            return None
    else:
        atoms = list(residue.get_atoms())
        return sum(atom.get_bfactor() for atom in atoms) / len(atoms) if atoms else None


def align_and_compare_bfactors(
        chain1: Chain,
        chain2: Chain,
        use_ca_only1: bool = True,
        use_ca_only2: bool = True,
) -> Tuple[List, List]:
    """
    Align sequences from two chains and compare their B-factors residue-by-residue.

    Args:
        chain1 (Chain): First chain (e.g., from bound structure).
        chain2 (Chain): Second chain (e.g., from unbound structure).
        use_ca_only (bool): If True, compare only CA atom B-factors.
                            If False, use average B-factor over all atoms.

    Returns:
        List[Tuple[str, Optional[float], Optional[float]]]:
            A list of tuples: (residue_code, bfactor_chain1, bfactor_chain2).
            B-factors may be None if residue is missing or CA atom is not found.
    """
    seq1 = extract_chain_sequence(chain1)
    seq2 = extract_chain_sequence(chain2)

    alignment = pairwise2.align.globalxx(seq1, seq2)[0]
    aligned_seq1, aligned_seq2 = alignment.seqA, alignment.seqB

    residues1 = get_ordered_residues(chain1)
    residues2 = get_ordered_residues(chain2)

    i1 = i2 = 0
    bfactor_seq1 = []
    bfactor_seq2 = []

    for a1, a2 in zip(aligned_seq1, aligned_seq2):
        if a1 != '-' and a2 != '-':
            b1 = get_bfactor(residues1[i1], use_ca_only1)
            b2 = get_bfactor(residues2[i2], use_ca_only2)
            bfactor_seq1.append((i1, a1, b1))
            bfactor_seq2.append((i2, a2, b2))
            i1 += 1
            i2 += 1
        elif a1 != '-' and a2 == '-':
            b1 = get_bfactor(residues1[i1], use_ca_only1)
            bfactor_seq1.append((i1, a1, b1))
            bfactor_seq2.append((i2, a2, 0))
            i1 += 1
        elif a1 == '-' and a2 != '-':
            b2 = get_bfactor(residues2[i2], use_ca_only2)
            bfactor_seq1.append((i1, a1, 0))
            bfactor_seq2.append((i2, a2, b2))
            i2 += 1

    return bfactor_seq1, bfactor_seq2
