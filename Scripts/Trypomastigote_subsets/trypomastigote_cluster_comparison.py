import scanpy as sc
import numpy as np
import pandas as pd
import os
import seaborn as sns
import matplotlib as plt
from anndata import AnnData
from typing import Sequence, Optional
import scipy.stats as stats

def flatten(l):
    return [item for sublist in l for item in sublist]

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

def rank_genes_all_groups_df(adata: AnnData):

    if not adata.uns['rank_genes_groups']:
        raise Exception("sc.tl.rank_genes_groups not run")

    groups_key = adata.uns['rank_genes_groups']['params']['groupby']
    groups = adata.obs[groups_key].unique().tolist()

    output_df = None
    for group in groups:
        df = sc.get.rank_genes_groups_df(adata, group=group)

        df.insert(0, 'group', [group]*len(df.index))

        if output_df is None:
            output_df = df
            print("why")
        else:
            output_df = pd.concat([output_df, df], join="outer")

    return output_df

os.chdir("/Users/rosslaidlaw/R/cruzi_paper/")

adata_concat = sc.read_h5ad("atlas/objects/cruzi_bbknn_noAmast_cellCycle.h5ad")

adata_concat.obs['leiden'] = adata_concat.obs['leiden'].cat.reorder_categories(["early_mid_amast_5", "late_amast_1",  "trypomast_0", "trypomast_4", "trypomast_9",
 "epimast_2", "epi_meta_trans_6", "epi_meta_trans_7","meta_3", "meta_8"])

adata_concat.obs['lifecycleStage'] = adata_concat.obs['leiden'].copy()

adata_concat = adata_concat[adata_concat.obs['lifecycleStage'].isin(['trypomast_0','trypomast_4'])].copy()


sc.tl.rank_genes_groups(adata_concat, 'leiden', method='wilcoxon', corr_method = "bonferroni")
sc.pl.rank_genes_groups(adata_concat, n_genes=25, sharey=False)

groups = rank_genes_all_groups_df(adata_concat)

groups.to_csv("trypo_analysis/DEG_tables/top_markers_trypo0_vs_trypo4_no48hr_all.csv")

groups_filter = groups[groups["pvals_adj"] < 0.05]
groups_filter = groups_filter[groups_filter["logfoldchanges"] > 0.5]

groups_filter.to_csv("trypo_analysis/DEG_tables/top_markers_trypo0_vs_trypo4_no48hr_no48hr_customFiltered.csv")

#cluster 0
sc.pl.violin(adata_concat, keys = "C4B63_331g11", groupby= "lifecycleStage", show = False)
plt.pyplot.savefig("/Users/rosslaidlaw/R/cruzi_paper/trypo_analysis/plots/cluster0_marker_C4B63_331g11_violin.png", bbox_inches="tight")

sc.pl.violin(adata_concat, keys = "C4B63_2g163c", groupby= "lifecycleStage", show = False)
plt.pyplot.savefig("/Users/rosslaidlaw/R/cruzi_paper/trypo_analysis/plots/cluster0_marker_C4B63_2g163c_violin.png", bbox_inches="tight")

sc.pl.violin(adata_concat, keys = "C4B63_57g87", groupby= "lifecycleStage", show = False)
plt.pyplot.savefig("/Users/rosslaidlaw/R/cruzi_paper/trypo_analysis/plots/cluster0_marker_C4B63_57g87_violin.png", bbox_inches="tight")

sc.pl.violin(adata_concat, keys = "C4B63_19g183", groupby= "lifecycleStage", show = False)
plt.pyplot.savefig("/Users/rosslaidlaw/R/cruzi_paper/trypo_analysis/plots/cluster0_marker_C4B63_19g183_violin.png", bbox_inches="tight")

sc.pl.violin(adata_concat, keys = "C4B63_6g542", groupby= "lifecycleStage", show = False)
plt.pyplot.savefig("/Users/rosslaidlaw/R/cruzi_paper/trypo_analysis/plots/cluster0_marker_C4B63_6g542_violin.png", bbox_inches="tight")

sc.pl.violin(adata_concat, keys = "C4B63_42g129", groupby= "lifecycleStage", show = False)
plt.pyplot.savefig("/Users/rosslaidlaw/R/cruzi_paper/trypo_analysis/plots/cluster0_marker_C4B63_42g129_violin.png", bbox_inches="tight")


#cluster 4
sc.pl.violin(adata_concat, keys = "C4B63_47g92", groupby= "lifecycleStage", show = False)
plt.pyplot.savefig("/Users/rosslaidlaw/R/cruzi_paper/trypo_analysis/plots/cluster4_marker_C4B63_47g92_violin.png", bbox_inches="tight")

sc.pl.violin(adata_concat, keys = "C4B63_158g42", groupby= "lifecycleStage", show = False)
plt.pyplot.savefig("/Users/rosslaidlaw/R/cruzi_paper/trypo_analysis/plots/cluster4_marker_C4B63_158g42_violin.png", bbox_inches="tight")

sc.pl.violin(adata_concat, keys = "C4B63_16g33", groupby= "lifecycleStage", show = False)
plt.pyplot.savefig("/Users/rosslaidlaw/R/cruzi_paper/trypo_analysis/plots/cluster4_marker_C4B63_16g33_violin.png", bbox_inches="tight")

sc.pl.violin(adata_concat, keys = "C4B63_33g179", groupby= "lifecycleStage", show = False)
plt.pyplot.savefig("/Users/rosslaidlaw/R/cruzi_paper/trypo_analysis/plots/cluster4_marker_C4B63_33g179_violin.png", bbox_inches="tight")

sc.pl.violin(adata_concat, keys = "C4B63_419g2", groupby= "lifecycleStage", show = False)
plt.pyplot.savefig("/Users/rosslaidlaw/R/cruzi_paper/trypo_analysis/plots/cluster4_marker_C4B63_419g2_violin.png", bbox_inches="tight")

#RBP10 brucei ortholog
sc.pl.violin(adata_concat, keys = "C4B63_16g295", groupby= "lifecycleStage", show = False)
plt.pyplot.savefig("/Users/rosslaidlaw/R/cruzi_paper/trypo_analysis/plots/cluster4_marker_C4B63_16g295_violin.png", bbox_inches="tight")

#ZC3H32 brucei ortholog
sc.pl.violin(adata_concat, keys = "C4B63_13g223", groupby= "lifecycleStage", show = False)
plt.pyplot.savefig("/Users/rosslaidlaw/R/cruzi_paper/trypo_analysis/plots/cluster4_marker_C4B63_13g223_violin.png", bbox_inches="tight")

#Extract the trans-sialidase genes and make a meta feature

dir_cruzi = "/Users/rosslaidlaw/R/cruzi_paper/"

allGenes = pd.read_csv(dir_cruzi + "dm28c_2018_genes.csv")

trans_sialidase_genes = allGenes[allGenes['Product Description'].str.contains('trans-sial')]
GP63_genes = allGenes[allGenes['Product Description'].str.contains('GP63')]

trans_sialidase = np.array([x for x in trans_sialidase_genes["Gene ID"] if x in adata_concat.raw.var_names])
GP63 = np.array([x for x in GP63_genes["Gene ID"] if x in adata_concat.raw.var_names])


trans_sialidase_tenP = np.array(flatten(
    np.sum(adata_concat.raw[:, trans_sialidase].X > 0, axis=0).tolist()))

trans_sialidase_tenP = trans_sialidase[np.where(trans_sialidase_tenP > len(adata_concat.obs["n_genes_by_counts"]) * 0.05)]

MetaFeature(adata_concat, obs_name="Trans_Sialidase", features=trans_sialidase_tenP)

adata_concat.obs["Trans_Sialidase_ratio"] = adata_concat.obs["Trans_Sialidase"] / np.mean(adata_concat.obs["Trans_Sialidase"], axis=0)

# sns.histplot(data=adata.obs, x="S_ratio")
# plt.pyplot.show()

GP63_tenP = np.array(flatten(
    np.sum(adata_concat.raw[:, GP63].X > 0, axis=0).tolist()))

GP63_tenP = GP63[np.where(GP63_tenP > len(adata_concat.obs["n_genes_by_counts"]) * 0.05)]

MetaFeature(adata_concat, obs_name="GP63", features=GP63_tenP)

adata_concat.obs["GP63_ratio"] = adata_concat.obs["GP63"] / np.mean(adata_concat.obs["GP63"], axis=0)

sc.pl.violin(adata_concat, keys = "Trans_Sialidase_ratio", groupby="lifecycleStage", show=False)
plt.pyplot.savefig("/Users/rosslaidlaw/R/cruzi_paper/trypo_analysis/plots/transSialidase_metaFeature.pdf", bbox_inches="tight",
                   format = "pdf")

sc.pl.violin(adata_concat, keys = "GP63_ratio", groupby="lifecycleStage", show=False)
plt.pyplot.savefig("/Users/rosslaidlaw/R/cruzi_paper/trypo_analysis/plots/GP63_metaFeature.png", bbox_inches="tight")

trypo_0 = adata_concat[adata_concat.obs['stages'] == "trypomast_0"]
trypo_4 = adata_concat[adata_concat.obs['stages'] == "trypomast_4"]
stats.ttest_ind(a=trypo_0.obs['Trans_Sialidase'], b=trypo_4.obs['Trans_Sialidase'], equal_var=True)


