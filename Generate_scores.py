import json
import joblib
import numpy as np
import pandas as pd
from numba import njit, prange
from scipy.spatial.distance import pdist
from tqdm import tqdm

from utils import argsort_mat

MIN_CLADE_SIZE = 20
clades_combs = json.load(open("clades_combs.json", "r"))

# %%
clades = json.load(open("clades_taxid.json", "r"))

with open("clades_taxid.json", "r") as f:
    clades = json.load(f)

    clades_dict = {}
    for k, v in clades.items():
        if v["size"] >= MIN_CLADE_SIZE:
            clades_dict[k] = v["taxid"]

NPP = pd.read_table("9606_NPP.tsv", delimiter="\t", index_col=0)
NPP_genes = {x: i for i, x in enumerate(NPP.index)}


def calculate_CladeOScope(calc="Min_Rank", combination="COMB5", fname=""):
    if combination == "COMB0":
        clades_keep = list(clades_dict.keys())
        print("all clades,", len(clades_keep))
    else:
        clades_keep = clades_combs[combination]

    BS_mat = pd.read_csv("9606_bitscore.tsv", sep="\t")
    BS_mat = BS_mat.iloc[:, 2:]
    BS_mat[BS_mat.isna()] = 0

    N = NPP.shape[0]

    calc = calc.lower()
    if calc == "min_rank":
        vals = np.full((N, N), 30000, dtype="float32")
    elif calc == "max_cor":
        vals = np.full((N, N), np.nan, dtype="float32")
    elif calc == "cummul_rank_log":
        vals = np.zeros((N, N), dtype="float32")
    elif calc == "npp_rank":
        clades_keep == ["Eukaryota"]

    for x in tqdm(clades_keep):
        # calc cor
        tmp = NPP.loc[:, clades_dict[x]].values.copy()
        tmp[np.quantile(BS_mat.loc[:, clades_dict[x]],
                        0.9, axis=1) <= 40, :] = np.nan
        tmp = np.corrcoef(tmp)
        if calc == "max_cor":
            vals = np.fmax(vals, tmp)
        else:
            # calc ranks
            tmp = argsort_mat(tmp).astype("float32") + 1
            tmp = np.sqrt(np.multiply(tmp, tmp.T))

            if calc == "npp_rank" and x == "Eukaryota":
                vals = tmp.copy()
            if calc == "min_rank":
                mask = tmp < vals
                vals[mask] = tmp[mask]
            elif calc == "cummul_rank_log":
                vals += np.log2(tmp) / len(clades_keep)
        del tmp

    vals = vals[np.triu_indices_from(vals, 1)]
    if calc == "max_cor":
        vals[np.isnan(vals)] = -1
    else:
        vals = -vals

    joblib.dump(vals, fname)
    del vals


def load_eval_mat():
    Eval_mat = pd.read_csv("9606_e-value.tsv", sep="\t")
    Eval_mat = Eval_mat.iloc[:, 2:].values
    Eval_mat = np.nan_to_num(Eval_mat)
    Eval_mat[Eval_mat <= 1e-3] = -1
    Eval_mat[Eval_mat > 1e-3] = 0
    Eval_mat = (Eval_mat * -1).astype('uint8')
    return Eval_mat


@njit(parallel=True)
def jaccard(X):
    M = X.shape[0]
    N = X.shape[1]
    D = np.zeros((M, M))
    sums = X.sum(axis=1)
    for i in prange(M - 1):
        if sums[i] == 0:
            continue
        for j in range(i + 1, M):
            if sums[j] == 0:
                continue
            u_sum = sums[i]
            v_sum = sums[j]
            dot = 0.0
            for k in range(N):
                dot += X[i, k] * X[j, k]
            if dot == 0:
                D[i, j] = 0
            else:
                D[i, j] = dot / (u_sum + v_sum - dot)
            D[j, i] = D[i, j]
    return D


def calculate_PPP(fname=""):
    Eval_mat = load_eval_mat()
    PPP_dist = jaccard(Eval_mat)
    PPP_dist_tmp = np.corrcoef(Eval_mat)

    PPP_dist = argsort_mat(PPP_dist)
    PPP_dist[PPP_dist_tmp < 0] = len(NPP_genes)

    PPP_dist = np.fmin(PPP_dist, PPP_dist.T)

    PPP_dist = -PPP_dist[np.triu_indices_from(PPP_dist, 1)]
    joblib.dump(PPP_dist, fname)
    del PPP_dist_tmp, PPP_dist, Eval_mat


def calcuate_BPP(fname=""):
    Eval_mat = load_eval_mat()

    # negative because of the implementation of the ROC expects higher -> better
    BPP_dist = -pdist(Eval_mat, "cityblock")
    BPP_dist = np.nan_to_num(BPP_dist)
    joblib.dump(BPP_dist, fname)
    print("BPP loaded and calculated")
    del BPP_dist, Eval_mat


def calculate_NPP_cor(fname=""):
    NPP_cors = np.corrcoef(NPP.iloc[:,1:].values)
    NPP_cors = NPP_cors[np.triu_indices_from(NPP_cors, 1)]
    joblib.dump(NPP_cors, fname)
    print("NPP loaded and calculated")
    del NPP_cors


map_to_func = {
    "BPP": calcuate_BPP,
    "PPP": calculate_PPP,
    "NPP_cors": calculate_NPP_cor,
    "CladeOScope": calculate_CladeOScope
}
if __name__ == "__main__":
    comps = json.load(open("comparisons.json", "r"))
    for k, v in comps.items():
        print(f"Calculating {k}")
        func = map_to_func[v['method']]

        func(fname=v['fname'], **v.get("kwargs", {}))
