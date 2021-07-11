# CladeOScope

Paper: [![DOI:10.1093/nargab/lqab024](https://img.shields.io/badge/DOI-https%3A%2F%2Fdoi.org%2F10.1093%2Fnargab%2Flqab024-blue)](https://doi.org/10.1093/nargab/lqab024)  |  Data: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4464120.svg)](https://doi.org/10.5281/zenodo.4464120)  |  Code: [Github Repo](https://github.com/dst1/CladeOScope)

----------------

Code for reproducing the method described in:

### **CladeOScope: elucidating functional interactions via a clade co-evolution prism**    
Tomer Tsaban†, Doron Stupp†, Dana Sherill-Rofe, Idit Bloch, Elad Sharon, Ora Schueler-Furman, Reuven Wiener and Yuval Tabach*

†. These authors contributed equally     
\*. To whom correspondence should be addressed. YuvalTab@ekmd.huji.ac.il

Published in NAR Genomics and Bioinformatics, Apr. 2021

**Citation:**
```{bibtex}
@article{Tsaban_2021,
	doi = {10.1093/nargab/lqab024},
	url = {https://doi.org/10.1093%2Fnargab%2Flqab024},
	year = 2021,
	month = {apr},
	publisher = {Oxford University Press ({OUP})},
	volume = {3},
	number = {2},
	author = {Tomer Tsaban and Doron Stupp and Dana Sherill-Rofe and Idit Bloch and Elad Sharon and Ora Schueler-Furman and Reuven Wiener and Yuval Tabach},
	title = {{CladeOScope}: functional interactions through the prism of clade-wise co-evolution},
	journal = {{NAR} Genomics and Bioinformatics}
} 
```

----------------

## Description

CladeOScope is a bioinformatic method for the prediction of functional interactions between human genes based on their phylogenetic profile.
CladeOScope combines the phylogenetic profile similarity across 16 clades and All Eukaryotes to better predict these functional interactions.

The code in this repo reproduces figure 5 and supp. figures 1 and 2 from the paper - 

recreating the ROC, partial ROC and precision-recall curves comparing CladeOScope with various parameters and other PP methods.

## Installation:

### Conda Enviroment

One can install the required software to reproduce our results by using conda. Clone this repo into a directory, and using the conda command prompt (or terminal on linux) run:

```{bash}
conda env create -f PATH_TO_REPO/enviroment.yaml -n ENV_NAME
```

Where `PATH_TO_REPO` is the path to which you cloned the repo and `ENV_NAME` is the enviroment name for your choosing.

### Pip

Another option to install the required software to reproduce our results is by using pip. However, it is recommended to first create a virtual enviroment using `venv` or a similar tool. To install with pip, clone this repo into a directory, and using the terminal from inside that directory run:

```{bash}
pip install -r requirements.txt

```

Please take note that the scripts were ran with python version 3.7.5. Unlike conda pip doesn't enforce a python version and thus it is up to the user to make sure their python version is compatible. It remains untested whether the output is reproducible using other versions.

### Data files:

To run this example one needs to download the associated phylogenetic profiling matrices deposited on Zenodo at: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4464120.svg)](https://doi.org/10.5281/zenodo.4464120)

Specifically make sure you have the files:
- `9606_bitscore.tsv`
- `9606_e-value.tsv`
- `9606_NPP.tsv`

inside the folder of the cloned repo.

## Running:

Using the command line, run from inside the cloned repo:
```{bash}
conda activate ENV_NAME # if using conda
python Generate_scores.py
python Test_KEGG.py
```

The first line is not required if installed via pip

- `Generate_scores.py` creates several pickle files which are then used for comparison, each is about 1.5GB in size.
- `Test_KEGG.py` calculates the ROC, partial ROC and precision-recall curves for these comparisons against a list of functional interactions from KEGG (and random negative pairs) found in `KEGG_pairs.tsv`.

To change which methods are tested, change `comparisons.json`:
e.g. for a different combination:
```{json}
"CladeOScope_comb4": {
        "method": "CladeOScope",
        "fname": "CladeOScope_comb4.pkl",
        "name": "Min_Rank (COMB4)",
        "kwargs": {
            "calc": "min_rank",
            "combination": "COMB4"
        }
```

For calculating the rank in all Eukaryotes (`NPP (rank)` in the paper):
```{json}
"NPP_rank": {
        "method": "CladeOScope",
        "fname": "NPP_rank.pkl",
        "name": "NPP_rank",
        "kwargs": {
            "calc": "NPP_rank"
        }
```

## Citation

Please make sure to cite when using the method:
- When using the NPP table please cite: 
```{bibtex}
@article{Tabach2012,
  doi = {10.1038/nature11779},
  url = {https://doi.org/10.1038/nature11779},
  year = {2012},
  month = dec,
  publisher = {Springer Science and Business Media {LLC}},
  volume = {493},
  number = {7434},
  pages = {694--698},
  author = {Yuval Tabach and Allison C. Billi and Gabriel D. Hayes and Martin A. Newman and Or Zuk and Harrison Gabel and Ravi Kamath and Keren Yacoby and Brad Chapman and Susana M. Garcia and Mark Borowsky and John K. Kim and Gary Ruvkun},
  title = {Identification of small {RNA} pathway genes using patterns of phylogenetic conservation and divergence},
  journal = {Nature}
}
```
- When using the CladeOScope method please cite 
```{bibtex}
@article{Tsaban_2021,
	doi = {10.1093/nargab/lqab024},
	url = {https://doi.org/10.1093%2Fnargab%2Flqab024},
	year = 2021,
	month = {apr},
	publisher = {Oxford University Press ({OUP})},
	volume = {3},
	number = {2},
	author = {Tomer Tsaban and Doron Stupp and Dana Sherill-Rofe and Idit Bloch and Elad Sharon and Ora Schueler-Furman and Reuven Wiener and Yuval Tabach},
	title = {{CladeOScope}: functional interactions through the prism of clade-wise co-evolution},
	journal = {{NAR} Genomics and Bioinformatics}
} 
```
