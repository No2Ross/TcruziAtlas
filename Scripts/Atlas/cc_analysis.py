import scanpy as sc
import numpy as np
import pandas as pd
import os

import scipy
import seaborn as sns
import matplotlib as plt
from anndata import AnnData
import anndata as ad
from scanpy import get
from typing import Sequence, Optional

def flatten(l):
    return [item for sublist in l for item in sublist]

def score_genes_cell_cycle_custom(
    adata: AnnData,
    s_genes: Sequence[str],
    g2m_genes: Sequence[str],
    g1_genes: Sequence[str],
    obs_name: str,
    threshold: int = 0,
    copy: bool = False,
    **kwargs,
) -> Optional[AnnData]:

    adata = adata.copy() if copy else adata
    ctrl_size = min(len(s_genes), len(g2m_genes))
    # add s-score
    sc.tl.score_genes(
        adata, gene_list=s_genes, score_name='S_score', ctrl_size=ctrl_size, **kwargs
    )
    # add g2m-score
    sc.tl.score_genes(
        adata,
        gene_list=g2m_genes,
        score_name='G2M_score',
        ctrl_size=ctrl_size,
        **kwargs,
    )

    # add g1-score
    sc.tl.score_genes(
        adata,
        gene_list=g1_genes,
        score_name='G1_score',
        ctrl_size=ctrl_size,
        **kwargs,
    )

    scores = adata.obs[['S_score', 'G2M_score', 'G1_score']]

    score_index = np.array(['S_score', 'G2M_score', 'G1_score'])

    phase = pd.Series(score_index[np.argmax(np.array(scores), axis = 1)], index = scores.index)

    # default phase is S
    # phase = pd.Series('S', index=scores.index)

    # if G2M is higher than S, it's G2M
    # phase[scores.G2M_score > scores.S_score & scores.G2M_score > scores.G1_score] = 'G2M'
    #
    # phase[scores.G1_score > scores.G2M_score & scores.G1_score > scores.S_score] = "G1"

    # if all scores are negative, it's Unlabelled...
    phase[np.all(scores < threshold, axis=1)] = 'Unlabelled'

    adata.obs[obs_name] = phase

    return adata if copy else None

def assign_label(scores, first='S', second='G2M', third='G1', null='Unlabelled'):
    if all(scores < 1.05):
        return null
    else:
        max_indices = np.where(scores == max(scores))[0]
        if len(max_indices) > 1:
            return 'Undecided'
        else:
            return [first, second, third][max_indices[0]]


def MetaFeature(adata: AnnData,
    features: Sequence[str],
    obs_name: str,
    copy: bool = False,
    **kwargs,)-> Optional[AnnData]:

    newmat = adata.raw[:,features].X.toarray()

    genetotals = np.sum( newmat, axis = 0 )

    newmat = newmat / genetotals

    newdata = np.mean(newmat, axis = 1)

    adata.obs[obs_name] = newdata

    return adata if copy else None


os.chdir("/Users/rosslaidlaw/R/cruzi_paper/")

sc.set_figure_params(figsize = (10,10))

adata_concat = sc.read_h5ad("atlas/objects/cruzi_bbknn_noAmast_cellCycle.h5ad")


adata_concat.X = adata_concat.raw.X


adata_concat.obs['leiden'] = adata_concat.obs['leiden'].cat.reorder_categories(["early_mid_amast_5", "late_amast_1",  "trypomast_0", "trypomast_4", "trypomast_9",
 "epimast_2", "epi_meta_trans_6", "epi_meta_trans_7","meta_3", "meta_8"])

datasets = ["early_mid_amast", "lifecycleMix1", "lifecycleMix2", "epimast_LateAmast",
            "trypomast", "meta1", "meta2"]

object_list = [adata_concat[adata_concat.obs['batch'] == i] for i in datasets ]

dir_chavez = "/Users/rosslaidlaw/R/cruzi/dm28cOrthologs_of_chavez_cruzi_cellCycle_genes_syntenic/"

chavez_G1 = pd.read_csv(dir_chavez + "G1.csv")
chavez_G2M = pd.read_csv(dir_chavez + "G2M.csv")
chavez_S = pd.read_csv(dir_chavez + "S.csv")

#Subset out genes where their are multiple input orthologs
chavez_G1_raw = chavez_G1.loc[chavez_G1["Input Ortholog(s)"].str.find(",") == -1 ,]
chavez_G2M_raw = chavez_G2M.loc[chavez_G2M["Input Ortholog(s)"].str.find(",") == -1 ,]
chavez_S_raw = chavez_S.loc[chavez_S["Input Ortholog(s)"].str.find(",") == -1 ,]

output_list = []

for i in range(np.shape(datasets)[0]):
    adata = object_list[i]

    chavez_G1 = np.array([x for x in chavez_G1_raw["Gene ID"] if x in adata.raw.var_names])
    chavez_G2M = np.array([x for x in chavez_G2M_raw["Gene ID"] if x in adata.raw.var_names])
    chavez_S = np.array([x for x in chavez_S_raw["Gene ID"] if x in adata.raw.var_names])

    chavez_G1_tenP = np.array(flatten(
        np.sum(adata.raw[:, chavez_G1].X > 0, axis=0).tolist()))

    chavez_G1_tenP = chavez_G1[np.where(chavez_G1_tenP > len(adata.obs["n_genes_by_counts"]) * 0.1)]

    chavez_G2M_tenP = np.array(flatten(
        np.sum(adata.raw[:, chavez_G2M].X > 0, axis=0).tolist()))

    chavez_G2M_tenP = chavez_G2M[np.where(chavez_G2M_tenP > len(adata.obs["n_genes_by_counts"]) * 0.1)]

    chavez_S_tenP = np.array(flatten(
        np.sum(adata.raw[:, chavez_S].X > 0, axis=0).tolist()))

    chavez_S_tenP = chavez_S[np.where(chavez_S_tenP > len(adata.obs["n_genes_by_counts"]) * 0.1)]

    MetaFeature(adata, obs_name="S", features=chavez_S_tenP)

    adata.obs["chavez_S_ratio"] = adata.obs["S"] / np.mean(adata.obs["S"], axis=0)

    # sns.histplot(data=adata.obs, x="S_ratio")
    # plt.pyplot.show()

    MetaFeature(adata, obs_name="G2M", features=chavez_G2M_tenP)

    adata.obs["chavez_G2M_ratio"] = adata.obs["G2M"] / np.mean(adata.obs["G2M"], axis=0)

    # sns.histplot(data=adata.obs, x="G2M_ratio")
    # plt.pyplot.show()

    MetaFeature(adata, obs_name="G1", features=chavez_G1_tenP)

    adata.obs["chavez_G1_ratio"] = adata.obs["G1"] / np.mean(adata.obs["G1"], axis=0)

    # sns.histplot(data=adata.obs, x="G1_ratio")
    # plt.pyplot.show()

    df_pd = adata.obs[["chavez_S_ratio", "chavez_G2M_ratio", "chavez_G1_ratio"]]

    assignments = df_pd.apply(lambda row: assign_label(row), axis=1)

    adata.obs["chavez_phase"] = assignments

    output_list.append(adata)

print("why")

adata_concat = ad.concat([output_list[inds] for inds in range(np.shape(datasets)[0])], merge="first")

tmp = pd.crosstab(adata_concat.obs['chavez_phase'],adata_concat.obs['leiden'], normalize='columns').T.plot(kind='bar', stacked=True)
tmp.legend(title='chavez_phase', bbox_to_anchor=(1.26, 1.02),loc='upper right')
plt.pyplot.savefig("cellCycle/chavez_cellCycle_cluster_phase.pdf", bbox_inches="tight", format = "pdf")


sc.pl.umap(adata_concat, color = "chavez_phase", show = False, size = 10)
plt.pyplot.savefig("cellCycle/cell_cycle_chavez_umap_plot.pdf", bbox_inches="tight", format = "pdf")

print("why")

#Brucei orthologs
object_list = [adata_concat[adata_concat.obs['batch'] == i] for i in datasets ]

dir_brucei = "/Users/rosslaidlaw/R/cruzi_10x_files/brucei_cc_files/"

brucei_G1_pcf_raw = pd.read_csv(dir_brucei + "PCF_brucei_G1_markers_briggs2023_Dm28c_syntenic_clean.csv")
brucei_G2M_pcf_raw = pd.read_csv(dir_brucei + "PCF_brucei_G2M_markers_briggs2023_Dm28c_syntenic_clean.csv")
brucei_S_pcf_raw = pd.read_csv(dir_brucei + "PCF_brucei_S_markers_briggs2023_Dm28c_syntenic_clean.csv")

output_list = []

for i in range(np.shape(datasets)[0]):
    adata = object_list[i]

    brucei_G1 = np.array([x for x in brucei_G1_pcf_raw["Gene.ID"] if x in adata.raw.var_names])
    brucei_G2M = np.array([x for x in brucei_G2M_pcf_raw["Gene.ID"] if x in adata.raw.var_names])
    brucei_S = np.array([x for x in brucei_S_pcf_raw["Gene.ID"] if x in adata.raw.var_names])

    brucei_G1_tenP = np.array(flatten(
        np.sum(adata.raw[:, brucei_G1].X > 0, axis=0).tolist()))

    brucei_G1_tenP = brucei_G1[np.where(brucei_G1_tenP > len(adata.obs["n_genes_by_counts"]) * 0.1)]

    brucei_G2M_tenP = np.array(flatten(
        np.sum(adata.raw[:, brucei_G2M].X > 0, axis=0).tolist()))

    brucei_G2M_tenP = brucei_G2M[np.where(brucei_G2M_tenP > len(adata.obs["n_genes_by_counts"]) * 0.1)]

    brucei_S_tenP = np.array(flatten(
        np.sum(adata.raw[:, brucei_S].X > 0, axis=0).tolist()))

    brucei_S_tenP = brucei_S[np.where(brucei_S_tenP > len(adata.obs["n_genes_by_counts"]) * 0.1)]

    MetaFeature(adata, obs_name="S", features=brucei_S_tenP)

    adata.obs["PCF_S_ratio"] = adata.obs["S"] / np.mean(adata.obs["S"], axis=0)

    # sns.histplot(data=adata.obs, x="S_ratio")
    # plt.pyplot.show()

    MetaFeature(adata, obs_name="G2M", features=brucei_G2M_tenP)

    adata.obs["PCF_G2M_ratio"] = adata.obs["G2M"] / np.mean(adata.obs["G2M"], axis=0)

    # sns.histplot(data=adata.obs, x="G2M_ratio")
    # plt.pyplot.show()

    MetaFeature(adata, obs_name="G1", features=brucei_G1_tenP)

    adata.obs["PCF_G1_ratio"] = adata.obs["G1"] / np.mean(adata.obs["G1"], axis=0)

    # sns.histplot(data=adata.obs, x="G1_ratio")
    # plt.pyplot.show()

    df_pd = adata.obs[["PCF_S_ratio", "PCF_G2M_ratio", "PCF_G1_ratio"]]

    assignments = df_pd.apply(lambda row: assign_label(row), axis=1)

    adata.obs["bruceiPCF_phase"] = assignments

    output_list.append(adata)

print("why")

adata_concat = ad.concat([output_list[inds] for inds in range(np.shape(datasets)[0])], merge="first")

tmp = pd.crosstab(adata_concat.obs['bruceiPCF_phase'],adata_concat.obs['leiden'], normalize='columns').T.plot(kind='bar', stacked=True)
tmp.legend(title='bruceiPCF_phase', bbox_to_anchor=(1.26, 1.02),loc='upper right')
plt.pyplot.savefig("cellCycle/bruceiPCF_cellCycle_cluster_phase.pdf", bbox_inches="tight", format = "pdf")


sc.pl.umap(adata_concat, color = "bruceiPCF_phase", show = False, size = 10)
plt.pyplot.savefig("cellCycle/cell_cycle_bruceiPCF_umap_plot.pdf", bbox_inches="tight", format = "pdf")


#Brucei BSF
object_list = [adata_concat[adata_concat.obs['batch'] == i] for i in datasets ]

brucei_G1_bsf_raw = pd.read_csv(dir_brucei + "BSF_brucei_G1_markers_briggs2023_Dm28c_syntenic_clean.csv")
brucei_G2M_bsf_raw = pd.read_csv(dir_brucei + "BSF_brucei_G2M_markers_briggs2023_Dm28c_syntenic_clean.csv")
brucei_S_bsf_raw = pd.read_csv(dir_brucei + "BSF_brucei_S_markers_briggs2023_Dm28c_syntenic_clean.csv")

output_list = []

for i in range(np.shape(datasets)[0]):
    adata = object_list[i]

    brucei_G1 = np.array([x for x in brucei_G1_pcf_raw["Gene.ID"] if x in adata.raw.var_names])
    brucei_G2M = np.array([x for x in brucei_G2M_pcf_raw["Gene.ID"] if x in adata.raw.var_names])
    brucei_S = np.array([x for x in brucei_S_pcf_raw["Gene.ID"] if x in adata.raw.var_names])

    brucei_G1_tenP = np.array(flatten(
        np.sum(adata.raw[:, brucei_G1].X > 0, axis=0).tolist()))

    brucei_G1_tenP = brucei_G1[np.where(brucei_G1_tenP > len(adata.obs["n_genes_by_counts"]) * 0.1)]

    brucei_G2M_tenP = np.array(flatten(
        np.sum(adata.raw[:, brucei_G2M].X > 0, axis=0).tolist()))

    brucei_G2M_tenP = brucei_G2M[np.where(brucei_G2M_tenP > len(adata.obs["n_genes_by_counts"]) * 0.1)]

    brucei_S_tenP = np.array(flatten(
        np.sum(adata.raw[:, brucei_S].X > 0, axis=0).tolist()))

    brucei_S_tenP = brucei_S[np.where(brucei_S_tenP > len(adata.obs["n_genes_by_counts"]) * 0.1)]

    MetaFeature(adata, obs_name="S", features=brucei_S_tenP)

    adata.obs["BSF_S_ratio"] = adata.obs["S"] / np.mean(adata.obs["S"], axis=0)

    sns.histplot(data=adata.obs, x="S_ratio")
    plt.pyplot.show()

    MetaFeature(adata, obs_name="G2M", features=brucei_G2M_tenP)

    adata.obs["BSF_G2M_ratio"] = adata.obs["G2M"] / np.mean(adata.obs["G2M"], axis=0)

    sns.histplot(data=adata.obs, x="G2M_ratio")
    plt.pyplot.show()

    MetaFeature(adata, obs_name="G1", features=brucei_G1_tenP)

    adata.obs["BSF_G1_ratio"] = adata.obs["G1"] / np.mean(adata.obs["G1"], axis=0)

    sns.histplot(data=adata.obs, x="G1_ratio")
    plt.pyplot.show()

    df_pd = adata.obs[["BSF_S_ratio", "BSF_G2M_ratio", "BSF_G1_ratio"]]

    assignments = df_pd.apply(lambda row: assign_label(row), axis=1)

    adata.obs["bruceiBSF_phase"] = assignments

    output_list.append(adata)

print("why")

adata_concat = ad.concat([output_list[inds] for inds in range(np.shape(datasets)[0])], merge="first")

tmp = pd.crosstab(adata_concat.obs['bruceiBSF_phase'],adata_concat.obs['leiden'], normalize='columns').T.plot(kind='bar', stacked=True)
tmp.legend(title='bruceiBSF_phase', bbox_to_anchor=(1.26, 1.02),loc='upper right')
plt.pyplot.savefig("cellCycle/bruceiBSF_cellCycle_cluster_phase.pdf", bbox_inches="tight", format = "pdf")


sc.pl.umap(adata_concat, color = "bruceiBSF_phase", show = False, size = 10)
plt.pyplot.savefig("cellCycle/cell_cycle_bruceiBSF_umap_plot.pdf", bbox_inches="tight", format = "pdf")



#Combine datasets analysis
chavez_G1 = np.array([x for x in chavez_G1_raw["Gene ID"] if x in adata_concat.raw.var_names])
chavez_G2M = np.array([x for x in chavez_G2M_raw["Gene ID"] if x in adata_concat.raw.var_names])
chavez_S = np.array([x for x in chavez_S_raw["Gene ID"] if x in adata_concat.raw.var_names])

chavez_G1_tenP = np.array(flatten(
    np.sum(adata_concat.raw[:, chavez_G1].X > 0, axis=0).tolist()))

chavez_G1_tenP = chavez_G1[np.where(chavez_G1_tenP > len(adata_concat.obs["n_genes_by_counts"]) * 0.1)]

chavez_G2M_tenP = np.array(flatten(
     np.sum(adata_concat.raw[:, chavez_G2M].X > 0, axis=0).tolist()))

chavez_G2M_tenP = chavez_G2M[np.where(chavez_G2M_tenP > len(adata_concat.obs["n_genes_by_counts"]) * 0.1)]

chavez_S_tenP = np.array(flatten(
      np.sum(adata_concat.raw[:, chavez_S].X > 0, axis=0).tolist()))

chavez_S_tenP = chavez_S[np.where(chavez_S_tenP > len(adata_concat.obs["n_genes_by_counts"]) * 0.1)]

MetaFeature(adata_concat, obs_name="S", features=chavez_S_tenP)

adata_concat.obs["chavez_S_ratio"] = adata_concat.obs["S"] / np.mean(adata_concat.obs["S"], axis=0)

# sns.histplot(data=adata.obs, x="S_ratio")
# plt.pyplot.show()

MetaFeature(adata_concat, obs_name="G2M", features=chavez_G2M_tenP)

adata_concat.obs["chavez_G2M_ratio"] = adata_concat.obs["G2M"] / np.mean(adata_concat.obs["G2M"], axis=0)

# sns.histplot(data=adata.obs, x="G2M_ratio")
# plt.pyplot.show()

MetaFeature(adata_concat, obs_name="G1", features=chavez_G1_tenP)

adata_concat.obs["chavez_G1_ratio"] = adata_concat.obs["G1"] / np.mean(adata_concat.obs["G1"], axis=0)

# sns.histplot(data=adata.obs, x="G1_ratio")
# plt.pyplot.show()

df_pd = adata_concat.obs[["chavez_S_ratio", "chavez_G2M_ratio", "chavez_G1_ratio"]]

assignments = df_pd.apply(lambda row: assign_label(row), axis=1)

adata_concat.obs["chavez_phase"] = assignments

print("why")

tmp = pd.crosstab(adata_concat.obs['chavez_phase'],adata_concat.obs['leiden'], normalize='columns').T.plot(kind='bar', stacked=True)
tmp.legend(title='chavez_phase', bbox_to_anchor=(1.26, 1.02),loc='upper right')
plt.pyplot.savefig("cellCycle/ChavezPhase_concat_cellCycle_cluster_phase.pdf", bbox_inches="tight", format = "pdf")



#Combine datasets analysis
brucei_G1 = np.array([x for x in brucei_G1_pcf_raw["Gene.ID"] if x in adata.raw.var_names])
brucei_G2M = np.array([x for x in brucei_G2M_pcf_raw["Gene.ID"] if x in adata.raw.var_names])
brucei_S = np.array([x for x in brucei_S_pcf_raw["Gene.ID"] if x in adata.raw.var_names])

brucei_G1_tenP = np.array(flatten(
        np.sum(adata.raw[:, brucei_G1].X > 0, axis=0).tolist()))

brucei_G1_tenP = brucei_G1[np.where(brucei_G1_tenP > len(adata.obs["n_genes_by_counts"]) * 0.1)]

brucei_G2M_tenP = np.array(flatten(
        np.sum(adata.raw[:, brucei_G2M].X > 0, axis=0).tolist()))

brucei_G2M_tenP = brucei_G2M[np.where(brucei_G2M_tenP > len(adata.obs["n_genes_by_counts"]) * 0.1)]

brucei_S_tenP = np.array(flatten(
        np.sum(adata.raw[:, brucei_S].X > 0, axis=0).tolist()))

brucei_S_tenP = brucei_S[np.where(brucei_S_tenP > len(adata.obs["n_genes_by_counts"]) * 0.1)]

MetaFeature(adata, obs_name="S", features=brucei_S_tenP)

adata.obs["BSF_S_ratio"] = adata.obs["S"] / np.mean(adata.obs["S"], axis=0)

# sns.histplot(data=adata.obs, x="S_ratio")
# plt.pyplot.show()

MetaFeature(adata, obs_name="G2M", features=brucei_G2M_tenP)

adata.obs["BSF_G2M_ratio"] = adata.obs["G2M"] / np.mean(adata.obs["G2M"], axis=0)

# sns.histplot(data=adata.obs, x="G2M_ratio")
# plt.pyplot.show()

MetaFeature(adata, obs_name="G1", features=brucei_G1_tenP)

adata.obs["BSF_G1_ratio"] = adata.obs["G1"] / np.mean(adata.obs["G1"], axis=0)

 # sns.histplot(data=adata.obs, x="G1_ratio")
 # plt.pyplot.show()

df_pd = adata.obs[["BSF_S_ratio", "BSF_G2M_ratio", "BSF_G1_ratio"]]

assignments = df_pd.apply(lambda row: assign_label(row), axis=1)

adata.obs["bruceiBSF_phase"] = assignments

tmp = pd.crosstab(adata_concat.obs['bruceiBSF_phase'],adata_concat.obs['leiden'], normalize='columns').T.plot(kind='bar', stacked=True)
tmp.legend(title='bruceiBSF_phase', bbox_to_anchor=(1.26, 1.02),loc='upper right')
plt.pyplot.savefig("cellCycle/bruceiBSFPhase_concat_cellCycle_cluster_phase.pdf", bbox_inches="tight", format = "pdf")
