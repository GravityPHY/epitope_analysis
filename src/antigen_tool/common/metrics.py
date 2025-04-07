import numpy as np
from sklearn.metrics import (roc_auc_score,
                             average_precision_score,
                             auc,
                             balanced_accuracy_score,
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


def cal_roc_auc(bfactor_pred, bfactor_true, value_index=-1):
    y_pred = []
    y_true = []
    for pred_val, true_val in zip(bfactor_pred,bfactor_true):
        y_pred.append(pred_val[value_index])
        y_true.append(true_val[-1])
    return roc_auc_score(y_true, y_pred)


def cal_pr_auc(bfactor_pred,bfactor_true, value_index=-1):
    y_pred = []
    y_true = []
    for pred_val, true_val in zip(bfactor_pred, bfactor_true):
        y_pred.append(pred_val[value_index])
        y_true.append(true_val[value_index])
    precision, recall, thresholds = precision_recall_curve(y_true, y_pred)
    auc_precision_recall = auc(recall, precision)
    return auc_precision_recall