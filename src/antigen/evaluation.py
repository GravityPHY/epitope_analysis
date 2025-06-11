import os
import sys
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sklearn.metrics as metrics
from sklearn.metrics import (roc_auc_score,
                             average_precision_score,
                             auc,
                             balanced_accuracy_score,
                             roc_curve,
                             precision_recall_curve)


def F1(TP, FP, FN):
    return TP / (TP + 0.5 * (FP + FN))


def MCC(TP, FP, FN, TN):
    return (TP * TN - FP * FN) / np.sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))


def find_top_cutoff(tuple, top=10):
    """
    find the cutoff for the first top values in the list,
    if top was set to a value longer than the length of tuple, it will be the minimum in the list
    (in order to avoid out of index error)
    :param tuple:
    :param top:
    :return:
    """
    tuple_length = len(tuple)
    new_tuple = sorted(tuple, key=lambda x: x[-1], reverse=True)
    if top < tuple_length:
        return new_tuple[top][1]
    else:
        return new_tuple[-1][1]


def cal_confusion_matrix_scores(bfactor_pred, bfactor_true, top=None):
    try:
        assert len(bfactor_pred) == len(bfactor_true)
    except:
        print(len(bfactor_pred), '!=', len(bfactor_true))
    TP, TN, FP, FN = 0, 0, 0, 0
    tot = 0
    cutoff = 0.2
    if top is not None:
        cutoff = find_top_cutoff(bfactor_pred, top)
    for t1, t2 in zip(bfactor_pred, bfactor_true):
        if t2[0] != '0':
            # tot += 1
            if t1[-1] >= cutoff and int(t2[-1]) == 1:
                TP += 1
            elif t1[-1] >= cutoff and int(t2[-1]) == 0:
                FP += 1
            elif int(t2[-1]) == 0:
                TN += 1
            elif int(t2[-1]) == 1:
                FN += 1
    # assert tot == TP + TN + FP + FN
    return TP, FP, FN, TN


def get_optimal_threshold(y_true, y_pred):
    """Returns optimized threshold to maximize TPR/FPR"""
    # https://stackoverflow.com/questions/28719067/roc-curve-and-cut-off-point-python
    fpr_, tpr_, thresholds = metrics.roc_curve(y_true, y_pred)
    optimal_idx = np.argmax(tpr_ - fpr_)
    optimal_threshold = thresholds[optimal_idx]

    return optimal_threshold


def get_performance(y_true, y_pred, binary_threshold=None, verbose=True):
    """Print performance"""

    auc = np.round(metrics.roc_auc_score(y_true, y_pred), 3)
    pr_auc = np.round(metrics.average_precision_score(y_true, y_pred, average="weighted"), 3)

    if binary_threshold is None:
        binary_threshold = get_optimal_threshold(y_true, y_pred)

    # Convert binary
    y_pred_binary = y_pred >= binary_threshold

    # tpr = np.round(metrics.recall_score(y_true, y_pred_binary), 3)
    precision = np.round(metrics.precision_score(y_true, y_pred_binary), 3)
    mcc = np.round(metrics.matthews_corrcoef(y_true, y_pred_binary), 3)
    f1 = np.round(metrics.f1_score(y_true, y_pred_binary), 3)
    opt_t = np.round(binary_threshold, 3)

    # Confusion matrix
    conf = metrics.confusion_matrix(y_true, y_pred_binary)
    tn, fp, fn, tp = conf.ravel()
    tnr = np.round(tn / (tn + fp), 3)
    tpr = np.round(tp / (tp + fn), 3)
    # fpr = np.round(fp / (fp + tn), 4)

    if verbose:
        print(f"AUC {auc}, PR-AUC {pr_auc}, MCC {mcc}, F1 {f1}, Prec {precision}")
        print(f"TPR/Recall {tpr}, TNR/Specificity {tnr}, Optimal threshold {opt_t}")
        print(f"Epitope length {y_true.count(1)}, Antigen length {len(y_true)}")

    return {
        "auc": auc,
        "pr": pr_auc,
        "optimal_threshold": binary_threshold,
        "mcc": mcc,
        "precision": precision,
        "recall": tpr,
        "f1": f1,
        "conf": conf,
        "epitope": y_true.count(1),
        "residues": len(y_true)
    }


def cal_roc_auc(bfactor_pred, bfactor_true, value_index=-1):
    y_pred = []
    y_true = []
    for pred_val, true_val in zip(bfactor_pred, bfactor_true):
        y_pred.append(pred_val[value_index])
        y_true.append(true_val[value_index])
    return roc_auc_score(y_true, y_pred)


def cal_pr_auc(bfactor_pred, bfactor_true, value_index=-1):
    y_pred = []
    y_true = []
    for pred_val, true_val in zip(bfactor_pred, bfactor_true):
        y_pred.append(pred_val[value_index])
        y_true.append(true_val[value_index])
    precision, recall, thresholds = precision_recall_curve(y_true, y_pred)
    auc_precision_recall = auc(recall, precision)
    return auc_precision_recall


def optimal_metrics(bfactor_pred, bfactor_true, value_index=-1):
    """

    Args:
        bfactor_pred (Tuple):
        bfactor_true (Tuple):
        value_index (int):
    Returns:

    """
    y_pred = []
    y_true = []
    for pred_val, true_val in zip(bfactor_pred, bfactor_true):
        y_pred.append(pred_val[value_index])
        y_true.append(true_val[value_index])
    get_performance(y_true, y_pred, binary_threshold=None, verbose=True)


def make_prauc_plot(y_pred, y_true, name=None, saving_path=None):
    precision, recall, thresholds = precision_recall_curve(y_true, y_pred)
    avg_pr = np.round(average_precision_score(y_true, y_pred, average="weighted"), 3)
    pr_auc = np.round(auc(recall, precision), 3)
    plt.xlim([-0.01, 1.01])
    plt.ylim([-0.01, 1.01])
    plt.ylabel('Precision')
    plt.xlabel('Recall')
    #plt.show()
    plt.plot(recall, precision, 'b', alpha=0.5)
    plt.title(f"{name} - avg pr:{avg_pr},pr auc:{pr_auc}")
    if saving_path:
        plt.savefig(saving_path)
    plt.close()


def make_rocauc_plot(y_pred, y_true, name=None, saving_path=None):
    fpr, tpr, thresholds = roc_curve(y_true, y_pred)

    roc_auc = np.round(roc_auc_score(y_true, y_pred), 3)
    plt.xlim([-0.01, 1.01])
    plt.ylim([-0.01, 1.01])
    plt.ylabel('TPR')
    plt.xlabel('FPR')
    #plt.show()
    plt.plot(fpr, tpr, 'b', alpha=0.5)
    plt.title(f"{name} - ROC-AUC {roc_auc}")
    if saving_path:
        plt.savefig(saving_path)
    plt.close()
