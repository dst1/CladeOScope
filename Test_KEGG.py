# %%
import json

import joblib
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from tqdm import tqdm

from utils import PairsIndex, plot_PR_single, plot_roc_single

mpl.rcParams['legend.loc'] = 'lower right'

clades_combs = json.load(open("clades_combs.json", "r"))

# %%
suffix = ""
size = 7.5
sns.set_palette(sns.color_palette("tab10") + sns.color_palette("pastel"))
target = "KEGG"

# %%
NPP = pd.read_table("9606_NPP.tsv", delimiter="\t", index_col=0)
NPP_genes = {x: i for i, x in enumerate(NPP.index)}
ind = PairsIndex(list(NPP.index), sort=False)
del NPP

# %%
comps = json.load(open("comparisons.json", "r"))
comps = [(v['name'], v['fname']) for v in comps.values()]


# %%

df_pairs = pd.read_csv(f"{target}_pairs.tsv", sep="\t",
                       header=None, names=["gene_a", "gene_b", "status"])
print(df_pairs.shape)
df_pairs = df_pairs.loc[df_pairs.gene_a.isin(
    NPP_genes) & df_pairs.gene_b.isin(NPP_genes), :]
print(df_pairs.shape)

pos = df_pairs.loc[df_pairs['status'] == 1, :]
pos = [ind._calc_pair_loc((gene_a, gene_b))
       for gene_a, gene_b in zip(pos.gene_a, pos.gene_b)]
neg = df_pairs.loc[df_pairs['status'] == 0, :]
neg = [ind._calc_pair_loc((gene_a, gene_b))
       for gene_a, gene_b in zip(neg.gene_a, neg.gene_b)]

y_true = np.concatenate((np.ones_like(pos), np.zeros_like(neg)))

fig_roc, ax_roc = plt.subplots(figsize=(size, size))
fig_proc, ax_proc = plt.subplots(figsize=(size, size))
fig_pr, ax_pr = plt.subplots(figsize=(size, size))
for comp_name, comp_vals in tqdm(comps):
    comp_vals = joblib.load(comp_vals)
    comp_vals = np.nan_to_num(comp_vals)
    y_pred = np.concatenate((comp_vals[pos], comp_vals[neg]))

    ax_roc = plot_roc_single(y_true, y_pred, comp_name, title=f"{target} ROC",
                             ax=ax_roc)
    ax_proc = plot_roc_single(y_true, y_pred, comp_name, title=f"{target} pROC",
                              ax=ax_proc, max_fpr=0.1)
    ax_proc.set_ylim((0, 0.5))

    ax_pr = plot_PR_single(y_true, y_pred, target="", label=comp_name,
                           title=f"{target} PR", ax=ax_pr)
    ax_pr.legend(loc='upper right')
fig_pr.savefig(f"PR_{target}{suffix}.png", dpi=300)
fig_pr.savefig(f"PR_{target}{suffix}.svg")
fig_roc.savefig(f"ROC_{target}{suffix}.png", dpi=300)
fig_roc.savefig(f"ROC_{target}{suffix}.svg")
fig_proc.savefig(f"pROC_{target}{suffix}.png", dpi=300)
fig_proc.savefig(f"pROC_{target}{suffix}.svg")
