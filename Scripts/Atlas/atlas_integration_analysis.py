import scanpy as sc
import numpy as np
import pandas as pd
import os
import seaborn as sns
import matplotlib as plt
from anndata import AnnData
from typing import Sequence, Optional

def flatten(l):
    return [item for sublist in l for item in sublist]

os.chdir("/Users/rosslaidlaw/R/cruzi_10x_files/")


datasets = ["early_mid_amast", "meta1", "meta2", "trypomast", "lifecycleMix1", "lifecycleMix2",  "epimast_LateAmast"]

features_common = []

for i in range(np.shape(datasets)[0]):
    adata = sc.read_h5ad("processed_files/"+ datasets[i] +".h5ad")

    features_common.append(list(adata.var_names[adata.var_names.str.startswith(("C", "M"))]))

common_genes = list(set.intersection(*map(set,features_common)))

object_list = []
object_list_raw = []

#Removing any genes which are only present in one sample. Might want to change this

#Cell cycle scoring - Chavez Cruzi

dir_chavez = "/Users/rosslaidlaw/R/cruzi/dm28cOrthologs_of_chavez_cruzi_cellCycle_genes_syntenic/"

chavez_G1 = pd.read_csv(dir_chavez + "G1.csv")
chavez_G2M = pd.read_csv(dir_chavez + "G2M.csv")
chavez_S = pd.read_csv(dir_chavez + "S.csv")

#Subset out genes where their are multiple input orthologs
chavez_G1_raw = chavez_G1.loc[chavez_G1["Input Ortholog(s)"].str.find(",") == -1 ,]
chavez_G2M_raw = chavez_G2M.loc[chavez_G2M["Input Ortholog(s)"].str.find(",") == -1 ,]
chavez_S_raw = chavez_S.loc[chavez_S["Input Ortholog(s)"].str.find(",") == -1 ,]

object_list = []
object_list_raw = []

for i in range(np.shape(datasets)[0]):
    adata = sc.read_h5ad("processed_files/"+ datasets[i] +".h5ad")

    # Doing a outer join concatenation of the objects doesn't work with the mt, human and cruzi information in var (due to using float64 nan, rather than float nan)
    # Adding then removing them for each dataset

    if datasets[i] in ["early_mid_amast", "lifecycleMix1", "lifecycleMix2", "epimast_LateAmast"]:
        del adata.var["cruzi"]
        del adata.var["human"]

    del adata.var["mt"]

    object_list.append(adata)

adata_concat = object_list[0].concatenate(object_list[1:len(object_list)], batch_categories = datasets, join = "outer")

adata_concat_cp = adata_concat.copy()

#Take the top 5000 genes, just so we can save out the counts without blowing the memory
sc.pp.highly_variable_genes(adata_concat_cp, n_top_genes=5000)
np.savetxt("epi_to_meta_expMtxRaw.csv", adata_concat_cp.layers["counts"].toarray(), delimiter=",")
adata_concat_cp = adata_concat_cp[:,adata_concat_cp.var.highly_variable]
del adata_concat_cp

sc.pp.highly_variable_genes(adata_concat, n_top_genes=2500)


sc.pl.highly_variable_genes(adata_concat)
genes = adata_concat.var.highly_variable

adata_concat.var.highly_variable["MT tcruzi"] = False
adata_concat.var.highly_variable["MT2tcruzi"] = False

#Subsetting out the variable genes gets rid of them in the counts layer.
#Need to keep those genes for tradeSeq analysis of metacyclogenesis

# adata_concat = adata_concat[:,adata_concat.var.highly_variable]

sc.pp.regress_out(adata_concat, ['gene_counts_cruzi'])


sc.pp.scale(adata_concat, max_value=10)

#By adding the use_highly_variable parameter we don't need to subset the matrix
sc.pp.pca(adata_concat, svd_solver = "arpack", use_highly_variable = True)
sc.pl.pca_variance_ratio(adata_concat, log=True, n_pcs=50, save='')

adata_concat_ridge = adata_concat.copy()

print(np.var(adata_concat.obsm["X_pca"], axis = 0))

sc.external.pp.bbknn(adata_concat, batch_key='batch', n_pcs = 17)  # running bbknn 1.3.6

sc.tl.umap(adata_concat)
sc.pl.umap(adata_concat, color = "sample")

sc.tl.leiden(adata_concat, resolution = 0.8)

sc.pl.umap(adata_concat, color = "leiden", legend_loc="on data")

new_cluster_names = [
    "trypomast_0", "late_amast_1", "epimast_2",
"meta_3", "trypomst_4", "early_mid_amast_5",
"epi_meta_trans_6", "epi_meta_trans_7", "meta_8",
"trypomast_9"]

adata_concat.rename_categories('leiden', new_cluster_names)

adata_concat.obs['stages'] = adata_concat.obs['leiden']

adata_concat.write("cruzi_bbknn_noAmast_cellCycle.h5ad")