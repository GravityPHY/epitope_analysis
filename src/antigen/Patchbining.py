from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB.SASA import ShrakeRupley
import numpy as np
import pandas as pd

class ProteinPatchAnalyzer:
    def __init__(self, pdb_file: str, probe_radius: float = 5.0):
        self.pdb_file = pdb_file
        self.probe_radius = probe_radius
        self.structure = None
        self.res_coords = []
        self.res_ids = []
        self.bfactors = []
        self.sasas = []
        self.patches = []

    def load_structure(self):
        parser = PDBParser(QUIET=True)
        self.structure = parser.get_structure("protein", self.pdb_file)
        sr = ShrakeRupley(probe_radius=self.probe_radius)
        sr.compute(self.structure, level="R")

    def extract_residue_data(self):
        for model in self.structure:
            for chain in model:
                for residue in chain:
                    #if is_aa(residue) and "CA" in residue:
                    coord = residue["CA"].coord
                    b_avg = np.mean([atom.bfactor for atom in residue])
                    self.res_coords.append(coord)
                    self.bfactors.append(b_avg)
                    self.res_ids.append((chain.id, residue.id))
                    self.sasas.append(residue.sasa)

        self.res_coords = np.array(self.res_coords)
        self.bfactors = np.array(self.bfactors)
        self.sasas = np.array(self.sasas)

    def grow_patch(self, 
                   center_idx, 
                   selected, 
                   is_surface, 
                   nearby_radius,
                   min_residues_in_patch,
                   max_residues_in_patch,
                   bvalue_cutoff
                   ):
        if not is_surface[center_idx]:
            return []
        dists = np.linalg.norm(self.res_coords - self.res_coords[center_idx], axis=1)
        nearby = np.where((dists <= nearby_radius) & (~selected) & is_surface & (self.bfactors > 0))[0]

        if len(nearby) > max_residues_in_patch:
            return nearby[np.argsort(-self.bfactors[nearby])[:max_residues_in_patch]].tolist()

        patch = [center_idx]
        remaining = set(nearby) - {center_idx}
        while len(patch) < max_residues_in_patch and remaining:
            added = False
            for idx in list(remaining):
                if any(
                    np.linalg.norm(self.res_coords[idx] - self.res_coords[i]) <= 8.0  for i in patch
                ):
                    patch.append(idx)
                    remaining.remove(idx)
                    added = True
            if not added:
                break
        return patch if len(patch) >= min_residues_in_patch else []

    def build_patches(self,
                      nearby_radius=15.0,
                   min_residues_in_patch:int=10,
                   max_residues_in_patch:int=30,
                   bvalue_cutoff=0.01):
        selected = np.zeros(len(self.res_ids), dtype=bool)
        is_surface = self.sasas > 0
        print(self.pdb_file)
        while not np.all(selected | (self.bfactors <= 0)):
            unselected = np.where(~selected & is_surface & (self.bfactors > 0))[0]
            if len(unselected) == 0:
                break
            center_idx = unselected[np.argmax(self.bfactors[unselected])]
            patch = self.grow_patch(center_idx, selected, is_surface,
                                    nearby_radius,
                                    min_residues_in_patch,
                                    max_residues_in_patch,
                                    bvalue_cutoff)
            if patch:
                print(self.res_ids[center_idx],patch)
                selected[patch] = True
                self.patches.append(patch)
            else:
                selected[center_idx] = True

    def compile_results(self):
        patch_rows = []
        for patch_id, indices in enumerate(self.patches):
            for idx in indices:
                chain, res_id = self.res_ids[idx]
                patch_rows.append({
                    "patch_id": patch_id,
                    "chain": chain,
                    "residue_id": res_id,
                    "bfactor": self.bfactors[idx]
                })

        df = pd.DataFrame(patch_rows)
        return df

    def save_structure(self, output_file: str, df: pd.DataFrame):
        top3_patch_ids = (
            df.groupby("patch_id")["bfactor"]
            .mean()
            .sort_values(ascending=False)
            .index.tolist()
        )

        top3_df = df[df["patch_id"].isin(top3_patch_ids)].sort_values(by=["patch_id", "bfactor"], ascending=[True, False])
        res_to_patch = {
            (row["chain"], row["residue_id"]): row["patch_id"]
            for _, row in top3_df.iterrows()
        }

        unique_patch_ids = sorted(top3_df["patch_id"].unique(), reverse=True)
        patch_to_b = {pid: (i+1) * 20.0 for i, pid in enumerate(unique_patch_ids)}

        for model in self.structure:
            for chain in model:
                for residue in chain:
                    key = (chain.id, residue.id)
                    if key in res_to_patch:
                        patch_id = res_to_patch[key]
                        new_b = patch_to_b[patch_id]
                        for atom in residue:
                            atom.bfactor = new_b#np.mean([atom.bfactor for atom in residue])#new_b
                    else:
                        for atom in residue:
                            atom.bfactor = np.mean([atom.bfactor for atom in residue])

        io = PDBIO()
        io.set_structure(self.structure)
        io.save(output_file)

# Example usage

analyzer = ProteinPatchAnalyzer("/Users/gravityphy/Downloads/epitope_mapping_tool/RSV_case/AbEMap/8ULK_C_bound.pdb")
analyzer.load_structure()
analyzer.extract_residue_data()
analyzer.build_patches()
df = analyzer.compile_results()
analyzer.save_structure("8ulk_top5_patches_projected.pdb", df)

analyzer = ProteinPatchAnalyzer("/Users/gravityphy/Downloads/epitope_mapping_tool/RSV_case/AbEMap/6OE4_bound.pdb")
analyzer.load_structure()
analyzer.extract_residue_data()
analyzer.build_patches()
df = analyzer.compile_results()
analyzer.save_structure("6oe4_top5_patches_projected.pdb", df)

analyzer = ProteinPatchAnalyzer("/Users/gravityphy/Downloads/epitope_mapping_tool/RSV_case/AbEMap/5TOJ_bound.pdb")
analyzer.load_structure()
analyzer.extract_residue_data()
analyzer.build_patches()
df = analyzer.compile_results()
analyzer.save_structure("5toj_top5_patches_projected.pdb", df)