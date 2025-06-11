from Bio.PDB import PDBParser,PDBIO
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB.SASA import ShrakeRupley
from scipy.spatial import KDTree
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial.distance import pdist, squareform
import numpy as np
import pandas as pd



# Load structure
pdb_file = "/Users/gravityphy/Downloads/epitope_mapping_tool/RSV_case/AbEMap/8ULK_C_bound.pdb"  # Replace with your actual PDB path


parser = PDBParser(QUIET=True)
structure = parser.get_structure("protein", pdb_file)
sr = ShrakeRupley(probe_radius=5)
sr.compute(structure,level="R")

# Extract CA coordinates, B-factors, residue IDs
res_coords, res_ids, bfactors = [], [], []
sasas=[]

for model in structure:
    for chain in model:
        for residue in chain:
            if is_aa(residue) and "CA" in residue:
                coord = residue["CA"].coord
                b_avg = np.mean([atom.bfactor for atom in residue])
                res_coords.append(coord)
                bfactors.append(b_avg)
                res_ids.append((chain.id, residue.id))
                sasas.append(residue.sasa)

sasas= np.array(sasas)
res_coords = np.array(res_coords)
bfactors = np.array(bfactors)

# Determine surface residues using neighbor density
#tree = KDTree(res_coords)
#neighbor_counts = np.array([len(tree.query_ball_point(c, 10.0)) for c in res_coords])
is_surface = sasas >0

# Initialize patch tracking
selected = np.zeros(len(res_ids), dtype=bool)
patches = []

def grow_patch(center_idx):
    if not is_surface[center_idx]:
        return []
    dists = np.linalg.norm(res_coords - res_coords[center_idx], axis=1)
    nearby = np.where((dists <= 20.0) & (~selected) & is_surface & (bfactors>0))[0]

    if len(nearby) > 30:
        return nearby[np.argsort(-bfactors[nearby])[:30]].tolist()

    patch = [center_idx]
    remaining = set(nearby) - {center_idx}
    while len(patch) < 30 and remaining:
        added = False
        for idx in list(remaining):
            if any(
                np.linalg.norm(res_coords[idx] - res_coords[i]) <= 10.0 and
                abs(bfactors[idx] - bfactors[i]) <= 0.05 for i in patch
            ):
                patch.append(idx)
                remaining.remove(idx)
                added = True
        if not added:
            break
    return patch if len(patch) >= 5 else []

# Loop through all residues to build patches
while not np.all(selected | (bfactors <= 0)):
    unselected = np.where(~selected &  is_surface & (bfactors>0))[0]
    if len(unselected) == 0:
        break
    center_idx = unselected[np.argmax(bfactors[unselected])]
    
    patch = grow_patch(center_idx)
    
    if patch!=[]:
        print(res_ids[center_idx],patch)
        selected[patch] = True
        patches.append(patch)
    else:
        selected[center_idx] = True

# Compile result
patch_rows = []
for patch_id, indices in enumerate(patches):
    for idx in indices:
        chain, res_id = res_ids[idx]
        patch_rows.append({
            "patch_id": patch_id,
            "chain": chain,
            "residue_id": res_id,
            "bfactor": bfactors[idx]
        })


df = pd.DataFrame(patch_rows)

# Get top 3 patches by max B-factor
top3_patch_ids = (
    df.groupby("patch_id")["bfactor"]
    .mean()
    .sort_values(ascending=False)
    .index.tolist()
)

top3_df = df[df["patch_id"].isin(top3_patch_ids)].sort_values(by=["patch_id", "bfactor"], ascending=[True, False])
# Create residue key → patch_id lookup from top3_df
res_to_patch = {
    (row["chain"], row["residue_id"]): row["patch_id"]
    for _, row in top3_df.iterrows()
}

# Normalize patch IDs to B-factor values for coloring (e.g., 30.0, 60.0, 90.0)
unique_patch_ids = sorted(top3_df["patch_id"].unique(),reverse=True)
#unique_patch_ids = sorted(top3_df.groupby("patch_id")["bfactor"].mean().sort_values(ascending=False).index)
patch_to_b = {pid: (i+1) * 20.0 for i, pid in enumerate(unique_patch_ids)}

# Update B-factors in structure
for model in structure:
    for chain in model:
        for residue in chain:
            key = (chain.id, residue.id)
            if key in res_to_patch:
                patch_id = res_to_patch[key]
                new_b = patch_to_b[patch_id]
                for atom in residue:
                    atom.bfactor = new_b #residue['CA'].get_bfactor()  # assign color for visualization
            else:
                for atom in residue:
                    atom.bfactor = 0

# Save modified structure
io = PDBIO()
io.set_structure(structure)
io.save("8ulk_top5_patches_projected.pdb")