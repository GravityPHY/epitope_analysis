import os
import gc
import glob
import json
import shutil
import numpy as np
from pathlib import Path
from collections import defaultdict

from Bio.PDB import (PDBParser, MMCIFParser,
                     NeighborSearch, Selection, PDBIO)


class Mapping:
    def __init__(self, cutoff_distance: float = 5.0):
        self.cutoff_distance = cutoff_distance
        assert self.cutoff_distance >= 0.0, "cutoff_distance must be greater than 0."

    def get_receptor(self, receptor_path: str):
        parser = PDBParser()
        receptor_structure = parser.get_structure("receptor", receptor_path)
        return receptor_structure

    def get_ligand(self, ligand_path:str):
        parser = PDBParser()
        ligand_strucutre = parser.get_structure("ligand", ligand_path)
        return ligand_strucutre

    def get_ligand_with_rank(self, ligand_path:str):
        path=Path(ligand_path)
        parts=path.stem.split(".")
        try:
            number = int(parts[-1])
        except ValueError:
            raise ValueError
        parser = PDBParser()
        ligand_strucutre = parser.get_structure("ligand", ligand_path)
        return ligand_strucutre, number


    def get_ligand_interface_atoms(self, receptor_structure, ligand_structure):
        ligand_atoms = Selection.unfold_entities(ligand_structure, 'A')
        receptor_atoms = Selection.unfold_entities(receptor_structure, 'A')

        neighbor_search = NeighborSearch(receptor_atoms)
        ligand_surface_atoms = [(atom.get_serial_number(), atom) for atom in ligand_atoms if
                                neighbor_search.search(atom.coord, self.cutoff_distance)]
        return ligand_surface_atoms

    def get_atomic_contact_frequency(self, receptor, ligands:list):
        atom_numbers_list = {}
        for ligand in ligands:
            for t in self.get_ligand_interface_atoms(receptor, ligand):
                num, atom = t
                atom_numbers_list[num] = atom_numbers_list.get(num,0)+1
        return atom_numbers_list


    def project_interface_bfactors(self, structure, atom_numbers_dict, output_pdb):
        for atom in structure.get_atoms():
            atom.set_bfactor(0.0)
        for atom in structure.get_atoms():
            atomic_number = atom.get_serial_number()
            if atomic_number in atom_numbers_dict:
                atom.set_bfactor(atom_numbers_dict[atomic_number])
        io = PDBIO()
        io.set_structure(structure)
        io.save(output_pdb)
