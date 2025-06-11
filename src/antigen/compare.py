import os
from typing import Any, Dict, List, Mapping, Optional, Tuple

from Bio import pairwise2
from Bio.Seq import Seq
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure

import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import (roc_auc_score,
                             average_precision_score,
                             auc,
                             roc_curve,
                             precision_recall_curve)

from structure import get_structure_pdb
from structure import get_structure_cif
from structure import aa3to1

from evaluation import get_performance


class Comparison:
    def __init__(self, 
                 task_name,
                 native:Structure,
                 model:Structure):
        """
        Assume native and model only have one chain for now
        """
        self.task_name=task_name
        self.native_structure=native[0]
        self.model_structure=model[0]
        self.native_chain=next(native.get_chains())
        self.model_chain=next(model.get_chains())

    def chain_to_sequence(self,
                          chain: Chain) -> Seq:
        _aainfo = lambda r: (r.id[1], aa3to1.get(r.resname, 'X'))
        seq = [_aainfo(r)[1] for r in chain.get_residues()]
        return Seq(''.join(seq))
    
    def get_ordered_residues(self, chain: Chain) -> List:
        """
        Get residues in a chain, excluding heteroatoms and waters.

        Args:
            chain (Chain): Biopython Chain object.

        Returns:
            List[Residue]: Ordered list of standard amino acid residues.
        """

        return [res for res in chain.get_residues() if res.id[0] == ' ']
    
    def get_bfactor(self,
                    residue:Residue,
                    use_ca_only:bool=True):
        if use_ca_only:
            if 'CA' in residue:
                return residue['CA'].get_bfactor()
            else:
                return None
        else:
            atoms = list(residue.get_atoms())
        return sum(atom.get_bfactor() for atom in atoms)

    
    def align_seq_property_list(self,
        chain1: Chain,
        chain2: Chain,
        use_ca_only1: bool = True,
        use_ca_only2: bool = True
            ) -> Tuple[List, List]:
        """
        Align chain1 and chain2 by their sequence 

        Args:
            chain1 (Chain): First chain (e.g., from native structure).
            chain2 (Chain): Second chain (e.g., from model structure).
            use_ca_only (bool): If True, compare only CA atom B-factors.
                            If False, use average B-factor over all atoms 
                            belong to the same residue.

        Returns:

        """
        seq1 = self.chain_to_sequence(chain1)
        seq2 = self.chain_to_sequence(chain2)

        alignment = pairwise2.align.globalxx(seq1,seq2)[0]
        aligned_seq1, aligned_seq2 = alignment.seqA, alignment.seqB

        residues1 = self.get_ordered_residues(chain1)
        residues2 = self.get_ordered_residues(chain2)
        i1 = i2 = 0
        bfactor_seq1 = []
        bfactor_seq2 = []
        for a1,a2 in zip(aligned_seq1, aligned_seq2):
            if a1!= '-' and a2!='-':
                b1=self.get_bfactor(residues1[i1], use_ca_only1)
                b2=self.get_bfactor(residues2[i2], use_ca_only2)
                bfactor_seq1.append((i1,a1,b1))
                bfactor_seq2.append((i2,a2,b2))
                i1+=1
                i2+=1
            elif a1!='-' and a2=='-':
                b1=self.get_bfactor(residues1[i1], use_ca_only1)
                bfactor_seq1.append((i1, a1, b1))
                bfactor_seq2.append((i2, a2, 0))
                i1+=1
            elif a1=='-' and a2!='-':
                b2=self.get_bfactor(residues2[i2], use_ca_only2)
                bfactor_seq1.append((i1,a1,0))
                bfactor_seq2.append((i2,a2,b2))
                i2+=1
        return bfactor_seq1, bfactor_seq2
    
    def get_metric(self,
                        native_chain,
                        model_chain,):
        bfactor_true,bfactor_pred=self.align_seq_property_list(native_chain,model_chain)
        y_pred ,y_true= [], []
        for pred_val, true_val in zip(bfactor_pred, bfactor_true):
            y_pred.append(pred_val[-1])
            y_true.append(true_val[-1])
        return get_performance(y_true, y_pred)
    
    def get_roc_figure(self,
                   native_chain,
                   model_chain,
                   name=None,
                   save_path=None):
        bfactor_true,bfactor_pred=self.align_seq_property_list(native_chain,model_chain)
        y_pred ,y_true= [], []
        for pred_val, true_val in zip(bfactor_pred, bfactor_true):
            y_pred.append(pred_val[-1])
            y_true.append(true_val[-1])
        fpr, tpr, thresholds = roc_curve(y_true, y_pred)
        roc_auc = np.round(roc_auc_score(y_true,y_pred),3)
        plt.xlim([-0.01,1.01])
        plt.ylim([-0.01,1.01])
        plt.ylabel("TPR")
        plt.xlabel("FPR")
        plt.plot(fpr,tpr,'b',alpha=0.5)
        plt.title(f"{name} - ROC-AUC {roc_auc}")
        if save_path:
            plt.savefig(save_path)
        plt.close()
    
    def get_pred_true(self,native_chain,
                   model_chain,
                   use_ca_only1=True,
                   use_ca_only2=True):
        bfactor_true,bfactor_pred=self.align_seq_property_list(native_chain,
                                                               model_chain,
                                                               use_ca_only1,
                                                               use_ca_only2)
        y_pred ,y_true= [], []
        for pred_val, true_val in zip(bfactor_pred, bfactor_true):
            if pred_val[0]==true_val[0]:
                y_pred.append(pred_val[-1])
                y_true.append(true_val[-1])
        return y_pred,y_true
        
    


if __name__=="__main__":
    native_sturcture=get_structure_pdb("/Users/gravityphy/Downloads/epitope_mapping_tool/RSV_case/labeled/8ULK_CC_lig_labeled.pdb")
    model_structure=get_structure_pdb("/Users/gravityphy/Downloads/epitope_mapping_tool/src/antigen/8ulk_top5_patches_projected.pdb")

    test_comparison=Comparison("test",
                               native_sturcture,
                               model_structure)
    print(test_comparison.get_metric(test_comparison.native_chain,
                                            test_comparison.model_chain))

    
    native_sturcture=get_structure_pdb("/Users/gravityphy/Downloads/epitope_mapping_tool/RSV_case/labeled/6OE4_A_lig_labeled.pdb")
    model_structure=get_structure_pdb("/Users/gravityphy/Downloads/epitope_mapping_tool/src/antigen/6oe4_top5_patches_projected.pdb")

    test_comparison=Comparison("test",
                               native_sturcture,
                               model_structure)
    print(test_comparison.get_metric(test_comparison.native_chain,
                                            test_comparison.model_chain))

    
    
    native_sturcture=get_structure_pdb("/Users/gravityphy/Downloads/epitope_mapping_tool/RSV_case/labeled/5TOJ_A_lig_labeled.pdb")
    model_structure=get_structure_pdb("/Users/gravityphy/Downloads/epitope_mapping_tool/src/antigen/5toj_top5_patches_projected.pdb")

    test_comparison=Comparison("test",
                               native_sturcture,
                               model_structure)
    print(test_comparison.get_metric(test_comparison.native_chain,
                                            test_comparison.model_chain))


