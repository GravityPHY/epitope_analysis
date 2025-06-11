# 
# Date 2025/06/04
# Status [Developing]
# This script choosing the top scoring epitope residue and 
# define the neighboring residues as epitope patch
# 
# process one structure for now
# currently exploring
#   - neighbor radius (8 is better than 5)
#   - how to reset residue score

from Bio.PDB import PDBParser
from Bio.Seq import Seq
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure

from typing import Any, Dict, List, Mapping, Optional, Tuple

import copy
import random
import numpy as np

class Antigen:
    def __init__(self,name):
        self.name=name
        self.structure=None

    def set_structure(self, path:str,name:str=None)->None:
        """
        Args:
            path(str): assume it is a path to pdb file

        """
        parser = PDBParser()
        if not name:
            structure = parser.get_structure('pdb', path)
        else:
            structure = parser.get_structure(name, path)
        self.structure=structure

    def get_structure(self):
        """
        Args:
            path(string): assume it is a path to pdb file
        Returns:

        """
        if self.structure:
            return self.structure
        return None
    
    def show_structure(self,
                       chosen_residues:List[Residue]=None,
                       save_path=None,
                       set_random=False):
        chosen_residues_set = set(chosen_residues)
        structure_copy=self.structure.copy()
        for atom in structure_copy.get_atoms():
            parent_res=atom.get_parent()
            if parent_res in chosen_residues_set:
                if not set_random:
                    #atom.set_bfactor(parent_res['CA'].get_bfactor())
                    atom.set_bfactor(np.mean([a.bfactor for a in parent_res.get_atoms()]))
                else:
                    atom.set_bfactor(random.uniform(0, 100))
                #atom.set_bfactor(1.0)
            else:
                atom.set_bfactor(0.0)
        
        if save_path:
            from Bio.PDB import PDBIO
            io = PDBIO()
            io.set_structure(structure_copy)
            io.save(save_path)

    
    def get_surface_residues(self, 
                             probe_radius:float=5.0,
                             threshold:float=0.0,):
        """

        """
        from Bio.PDB.SASA import ShrakeRupley
        if self.structure is None:
            raise ValueError("Structure not set.")
        sr = ShrakeRupley(probe_radius=probe_radius)
        sr.compute(self.structure)
        surface_residues = []
        model=self.structure[0]
        for chain in model:
            for residue in chain:
                if residue.sasa>threshold:
                    surface_residues.append(residue)
        return surface_residues

    
    def get_top_patch_old(self, radius:float=25.0):
        """
        Find the surface residues within the highest 
        epitope likelihood score (value at b-factor column), 
        """
        from Bio.PDB import NeighborSearch
        from Bio.PDB.SASA import ShrakeRupley
        if self.structure is None:
            raise ValueError("Structure not set.")
        
        sr = ShrakeRupley(probe_radius=5.0)
        sr.compute(self.structure, level='R')
        model = self.structure[0]

        # find the top residues with the highest b-factors
        residues = []
        for chain in model:
            for residue in chain:
                bfac = residue['CA'].get_bfactor()
                residues.append((bfac, residue))
        residues.sort(reverse=True, key=lambda x: x[0])
        all_atoms = list(model.get_atoms())
        ns = NeighborSearch(all_atoms)

        def get_patch_residues(center_residue, exclude_residues=set()):
            center_ca = center_residue['CA']
            neighbor_atoms = ns.search(center_ca.get_coord(),
                                       radius,
                                       level='A')
            surface_residues=set()
            for atom in neighbor_atoms:
                parent_residue = atom.get_parent()
                if parent_residue.sasa>0 and parent_residue not in exclude_residues:
                    surface_residues.add(parent_residue)
            return surface_residues
        
        # Get the top patch
        top_residue = residues[0][1]
        top_patch_residues = get_patch_residues(top_residue)

        second_residue = None
        for bfac, residue in residues[1:]:
            if residue not in top_patch_residues:
                second_residue = residue
                break
        if second_residue:
            second_patch_residues = get_patch_residues(second_residue,
            exclude_residues=top_patch_residues)
        else:
            second_patch_residues = set()

        # Find the third highest residue not in the top or second patch
        third_residue = None
        for bfac, residue in residues[1:]:
            if residue not in top_patch_residues and residue not in second_patch_residues:
                third_residue = residue
                break
        if third_residue:
            third_patch_residues = get_patch_residues(third_residue,
            exclude_residues=top_patch_residues|second_patch_residues)
        else:
            third_patch_residues = set()
    
        return list(top_patch_residues), list(second_patch_residues), list(third_patch_residues)

    def get_top3_patch_residues(self,
                                radius:float=5,
                                threshold:float=0.0):
        from Bio.PDB.SASA import ShrakeRupley
        if self.structure is None:
            raise ValueError("Structure not set.")
        
        sr = ShrakeRupley(probe_radius=5.0)
        sr.compute(self.structure, level='R')
        model = self.structure[0]

        # find the top residues with the highest b-factors
        residues = []
        coords = []
        bfactors = []
        for chain in model:
            for residue in chain:
                atom = residue["CA"] if "CA" in residue else residue["CB"]
                coord = atom.coord
                b_avg = np.mean([a.bfactor for a in residue.get_atoms()])
                residues.append(residue)
                coords.append(coord)
                bfactors.append(b_avg)
        coords = np.array(coords)
        bfactors = np.array(bfactors)

        # Sort residues by B-factor descending
        sorted_indices = np.argsort(-bfactors)
        assigned = np.zeros(len(residues), dtype=bool)
        group_ids = -1 * np.ones(len(residues), dtype=int)

        group_counter = 0

    # Form spatially separated groups
        for idx in sorted_indices:
            if assigned[idx]:
                continue

        # Seed residue for this group
            seed_coord = coords[idx]
            dists = np.linalg.norm(coords - seed_coord, axis=1)
            group_members = np.where((dists <= radius) & (~assigned))[0]

            group_ids[group_members] = group_counter
            assigned[group_members] = True
            group_counter += 1

            # Stop after forming 3 groups
            if group_counter >= 3:
                break

     # Collect residues for the top 3 patches
        patch1 = [residues[i] for i in range(len(residues)) if group_ids[i] == 0]
        patch2 = [residues[i] for i in range(len(residues)) if group_ids[i] == 1]
        patch3 = [residues[i] for i in range(len(residues)) if group_ids[i] == 2]

        return patch1, patch2, patch3



    


if __name__=="__main__":
    test_antigen=Antigen("6oe4")
    test_antigen.set_structure("/Users/gravityphy/Downloads/epitope_mapping_tool/RSV_case/labeled/6OE4_A_lig_labeled.pdb")
    #print(test_antigen.get_top3_patch_residues())
    patch1,patch2,patch3=test_antigen.get_top3_patch_residues()#test_antigen.get_top_patch()
    test_antigen.show_structure(patch1,
                                save_path="./6oe4-1.pdb")
    test_antigen.show_structure(patch2,
                                save_path="./6oe4-2.pdb")

    test_antigen.show_structure(patch3,
                                save_path="./6oe4-3.pdb")
    surface=test_antigen.get_surface_residues()
    test_antigen.show_structure(surface,
                                save_path="./6oe4_surface.pdb",
                                set_random=True)
    
    test_antigen=Antigen("5toj")
    test_antigen.set_structure("/Users/gravityphy/Downloads/epitope_mapping_tool/RSV_case/labeled/5TOJ_A_lig_labeled.pdb")
    patch1,patch2,patch3=test_antigen.get_top_patch()
    test_antigen.show_structure(patch1,
                                save_path="./5toj-1.pdb")
    test_antigen.show_structure(patch2,
                                save_path="./5toj-2.pdb")

    test_antigen.show_structure(patch3,
                                save_path="./5toj-3.pdb")
    surface=test_antigen.get_surface_residues()
    test_antigen.show_structure(surface,
                                save_path="./5toj_surface.pdb",
                                set_random=True)
    
    test_antigen=Antigen("8ulk")
    test_antigen.set_structure("/Users/gravityphy/Downloads/epitope_mapping_tool/RSV_case/labeled/8ULK_CC_lig_labeled.pdb")
    patch1,patch2,patch3=test_antigen.get_top_patch()
    test_antigen.show_structure(patch1,
                                save_path="./8ulk-1.pdb")
    test_antigen.show_structure(patch2,
                                save_path="./8ulk-2.pdb")

    test_antigen.show_structure(patch3,
                                save_path="./8ulk-3.pdb")
    surface=test_antigen.get_surface_residues()
    test_antigen.show_structure(surface,
                                save_path="./8ulk_surface.pdb",
                                set_random=True)

    
    
        

