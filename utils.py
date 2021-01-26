from bisect import bisect

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import (average_precision_score, precision_recall_curve,
                             roc_auc_score, roc_curve)
from tqdm import tqdm


def argsort_mat(NPP_cor_mat):
    # prep argsort:
    adj_argsort = np.zeros((NPP_cor_mat.shape[0], NPP_cor_mat.shape[0]))
    for i in tqdm(range(NPP_cor_mat.shape[0])):
        inds = np.argsort(-NPP_cor_mat[i, :])
        adj_argsort[i, inds] = np.arange(NPP_cor_mat.shape[0])
    return adj_argsort


class PairsIndex:
    """
    Pairs index caches the indices for each genes and contains several utilities to match a pair with a sorted all pairs index
    """

    def __init__(self, genes, sort=True):
        """
        Initializes a new instance.
        An instance contains the list of genes and a mapping from gene to place in list

        :param genes: list, all genes in index
        :param sort: whether to use a sorted index
        """
        if sort == True:
            self.genes = sorted(genes)
        else:
            self.genes = genes

        self.ind = {k: i for i, k in enumerate(self.genes)}
        self.rev_ind = {i: k for i, k in enumerate(self.genes)}
        self.n = len(genes)
        self.ind_size = (self.n * (self.n - 1)) // 2
        self.i0_inds = [((2 * self.n - i_0 - 1) * i_0) //
                        2 for i_0 in range(self.n)]

    def get_ind_solo(self, genes):
        """
        Returns the place on the genes in the list

        Assumes genes are found on the list!!

        :param genes: genes to query
        :return: indexes for these genes
        """
        if type(genes) == str:
            return self.ind[genes]
        else:
            return [self.ind[x] for x in genes]

    def _calc_pair_loc(self, genes):
        """
        Internal function, gets an ordered tuple and calculates its location
        :param genes: list, gene[0]<gene[1], also not the same..
        :return: int, location in pairs list
        """

        i_0, i_1 = self.get_ind_solo(genes)
        ind = self._calc_pair_loc_from_inds(i_0, i_1)
        return ind

    def _calc_pair_loc_from_inds(self, x, y):
        """
        Internal function, gets an ordered tuple and calculates its location
        :param genes: list, gene[0]<gene[1], also not the same..
        :return: int, location in pairs list
        """

        i_0, i_1 = sorted([x, y])
        ind = ((2 * self.n - i_0 - 1) * i_0) // 2 + (i_1 - i_0 - 1)
        return ind

    def get_genes_solo(self, inds):
        """
        Returns the gene in place ind

        Assumes genes are found on the list!!

        :param genes: genes to query
        :return: indexes for these genes
        """
        if type(inds) == int:
            return self.rev_ind[inds]
        else:
            return [self.rev_ind[x] for x in inds]

    def _calc_pair_genes(self, ind, return_inds=False):
        """
        Internal function, gets a location and calulate the indexes of the genes
        returns the genes as a tuple
        :param ind: location in pairs list
        :return: tuple, genes
        """
        i_0 = bisect(self.i0_inds, ind) - 1
        i_1 = ind - self.i0_inds[i_0] + i_0 + 1

        if return_inds:
            return (i_0, i_1)
        return (self.get_genes_solo(i_0), self.get_genes_solo(i_1))


def plot_roc_single(y_test, probas, target, title="ROC", label=None,
                    ax=None, lw=2, figsize=(10, 10), max_fpr=1.0):

    if type(y_test) == pd.core.frame.DataFrame:
        y_test = y_test.values

    if ax is None:
        fig, ax = plt.subplots(1, 1, constrained_layout=True, figsize=figsize)

    fpr, tpr, _ = roc_curve(y_test, probas)
    roc_auc = roc_auc_score(y_test, probas, max_fpr=max_fpr)
    if label is None:
        ax.plot(fpr, tpr, lw=lw,
                label='{0} (AUC = {1:0.3f})'.format(target, roc_auc))
    else:
        ax.plot(fpr, tpr, lw=lw, label=label)
    ax.plot([0, 1], [0, 1], 'k--', lw=2)
    ax.set_xlim([0.0, max_fpr])
    ax.set_ylim([0.0, 1.05])
    ax.set_xlabel('False Positive Rate')
    ax.set_ylabel('True Positive Rate')
    ax.set_title(title)
    ax.legend()
    return ax


def plot_PR_single(y_test, probas, target, label=None,
                   title="Precision - Recall",
                   ax=None, lw=2, figsize=(10, 10)):

    if type(y_test) == pd.core.frame.DataFrame:
        y_test = y_test.values

    if ax is None:
        fig, ax = plt.subplots(1, 1, constrained_layout=True, figsize=figsize)

    prec, rec, _ = precision_recall_curve(y_test, probas)
    AP = average_precision_score(y_test, probas)
    if label is None:
        ax.plot(rec, prec, lw=lw,
                label='{0} (AP = {1:0.2f})'.format(target, AP))
    else:
        ax.plot(rec, prec, lw=lw, label=label)
    ax.set_xlabel('Recall')
    ax.set_ylabel('Precision')
    ax.set_ylim([0.0, 1.05])
    ax.set_xlim([0.0, 1.0])
    ax.set_title(title)
    ax.legend()
    return ax
