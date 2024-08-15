import scanpy as sc
import numpy as np
import pandas as pd
import os

import scipy
import seaborn as sns
import matplotlib as plt
from anndata import AnnData
from scanpy import get

def flatten(l):
    return [item for sublist in l for item in sublist]

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
# adata_concat_ind = sc.read_h5ad("cruzi_bbknn_scaleIndividual.h5ad")

adata_concat.obs['leiden'] = adata_concat.obs['leiden'].cat.reorder_categories(["early_mid_amast_5", "late_amast_1",  "trypomast_0", "trypomast_4", "trypomast_9",
 "epimast_2", "epi_meta_trans_6", "epi_meta_trans_7","meta_3", "meta_8"])

#Rename genes to be their broad lifecycle stage type, i.e. amastigote, trypomast, etc
#The transition ones will stay transition

old_to_new = dict(
    early_mid_amast_5='amasitgote',
    late_amast_1='amasitgote',
    trypomast_0='trypomastigote',
    trypomast_4='trypomastigote',
    trypomast_9='trypomastigote',
    epimast_2='epimastigote',
    epi_meta_trans_6='epi_meta_trans',
    epi_meta_trans_7='epi_meta_trans',
    meta_3='metacyclic',
    meta_8='metacyclic'
)
adata_concat.obs['broad_lifecycle'] = (
    adata_concat.obs['leiden']
    .map(old_to_new)
    .astype('category')
)

#The following section finds marker genes for the four different life cycle stages

sc.tl.rank_genes_groups(adata_concat, 'broad_lifecycle', method='wilcoxon', corr_method = "bonferroni")

sc.pl.rank_genes_groups(adata_concat, n_genes=25, sharey=False, show = False)
plt.pyplot.savefig("atlas/plots/stagesCoarse_rank_genes_groups_plot.pdf",format='pdf', bbox_inches="tight")


groups = rank_genes_all_groups_df(adata_concat)

groups.to_csv("atlas/DEG_tables/top_markers_stageCoarse_all.csv")

groups_filter = groups[groups["pvals_adj"] < 0.05]
groups_filter = groups_filter[groups_filter["logfoldchanges"] > 0.5]

groups_filter.to_csv("atlas/DEG_tables/top_markers_stagesCoarse_customfiltered.csv")


sc.tl.filter_rank_genes_groups(adata_concat, min_in_group_fraction = 0.5, max_out_group_fraction = 0.5)

filter_genes_pd = pd.DataFrame(adata_concat.uns['rank_genes_groups_filtered']['names'])

filter_genes_pd.to_csv("atlas/DEG_tables/top_markers_stagesCoarse_strict_minPCT_scanpyFiltered.csv")

sc.tl.filter_rank_genes_groups(adata_concat)

gene_names = []

for i in range(0,len(filter_genes_pd.keys())):

    current = filter_genes_pd[ filter_genes_pd.keys()[i] ].copy()

    current = current.dropna()

    gene_names.append(current.tolist())

gene_names = np.unique(flatten(gene_names))

filter_genes_pd.to_csv("atlas/DEG_tables/top_markers_stagesCoarse_scanpyFiltered.csv")

pd.DataFrame(adata_concat.uns['rank_genes_groups']['names']).to_csv("top_markers_stagesCoarse_raw.csv")

genes_plot = []
cluster_plot = []

for i in filter_genes_pd.columns.values:
    current = filter_genes_pd[i]

    current = current.dropna()

    current = current.tolist()

    print(i)

    if len(current) < 10:
        genes_plot.append(current)
        cluster_plot.append([i] * len(current))

    else:
        genes_plot.append(current[0:10])
        cluster_plot.append([i] * 10)


genes_plot = pd.Series(flatten(genes_plot))
cluster_plot = pd.Series(flatten(cluster_plot))

output_df = pd.concat({"genes":genes_plot, "cluster":cluster_plot}, axis = 1)

output_df.to_csv("atlas/DEG_tables/top_markers_stagesCoarse_filtered_heatmapPlot.csv")

while "MT tcruzi" in genes_plot:
    genes_plot.remove("MT tcruzi")

while "MT2tcruzi" in genes_plot:
    genes_plot.remove("MT2tcruzi")

sc.pl.heatmap(adata_concat,var_names = genes_plot, groupby='leiden', show = False)
plt.pyplot.savefig("atlas/plots/heatmap_rank_genes_stagesCoarse_filtered.pdf",format='pdf', bbox_inches="tight")


sc.set_figure_params(figsize = (20,60))

#Do dotplot of the marker genes for the different life cycle stages
markers = ['C4B63_341g12', 'C4B63_361g15', 'C4B63_32g230',"C4B63_34g140", 'C4B63_34g1526c', 'C4B63_12g175', 'C4B63_57g69', 'C4B63_212g9']
sc.pl.dotplot(adata_concat, markers, groupby='leiden', dendrogram=True, show = False)
plt.pyplot.savefig("atlas/plots/dotplot_marker_genes.pdf",format='pdf', bbox_inches="tight")



#The following code finds the marker genes of the atlas clusters
sc.tl.rank_genes_groups(adata_concat, 'leiden', method='wilcoxon', corr_method = "bonferroni")
sc.pl.rank_genes_groups(adata_concat, n_genes=25, sharey=False)

groups = rank_genes_all_groups_df(adata_concat)

groups.to_csv("atlas/DEG_tables/top_markers_clusters_all.csv")

groups_filter = groups[groups["pvals_adj"] < 0.05]
groups_filter = groups_filter[groups_filter["logfoldchanges"] > 0.5]

groups_filter.to_csv("atlas/DEG_tables/top_markers_clusters_customFiltered.csv")


sc.tl.filter_rank_genes_groups(adata_concat)

filter_genes_pd = pd.DataFrame(adata_concat.uns['rank_genes_groups_filtered']['names'])

pd.DataFrame(adata_concat.uns['rank_genes_groups']['names']).to_csv("atlas/DEG_tables/top_markers_clusters_raw.csv")

genes_plot = []
cluster_plot = []

for i in filter_genes_pd.columns.values:

    current = sc.get.rank_genes_groups_df(adata_concat, group=i)

    current = current[current["pvals_adj"] < 0.05]
    current = current[current["logfoldchanges"] > 0.5]
    current['stage'] = i

    print(i)

    if len(current['names']) < 20:
        genes_plot.append(current['names'])
        cluster_plot.append([i] * len(current['names']))

    else:
        genes_plot.append(current['names'][0:20])
        cluster_plot.append([i] * 20)

    current.to_csv("atlas/DEG_tables/clusterTopIndividual/top_markers_"+i+"_no48hr_filtered.csv")


genes_plot = pd.Series(flatten(genes_plot))
cluster_plot = pd.Series(flatten(cluster_plot))

while "MT tcruzi" in genes_plot:
    genes_plot.remove("MT tcruzi")

while "MT2tcruzi" in genes_plot:
    genes_plot.remove("MT2tcruzi")

sc.pl.heatmap(adata_concat,var_names = genes_plot, groupby='leiden', show = False)
plt.pyplot.savefig("atlas/plots/heatmap_rank_genes_filtered.pdf", format='pdf', bbox_inches="tight")



#The following section plots some visulisations of cluster distribution across UMAP and sample distribution across clusters and other such things

tmp = pd.crosstab(adata_concat.obs['sample'],adata_concat.obs['leiden'], normalize='columns').T.plot(kind='bar', stacked=True)
tmp.legend(title='Sample distribution', bbox_to_anchor=(1.26, 1.02),loc='upper right')
plt.pyplot.savefig("atlas/plots/stage_cluster_batch_barplot.pdf",format="pdf", bbox_inches="tight")

tmp = pd.crosstab(adata_concat.obs['sample'],adata_concat.obs['stages'], normalize='columns').T.plot(kind='bar', stacked=True)
tmp.legend(title='Sample distribution', bbox_to_anchor=(1.26, 1.02),loc='upper right')
plt.pyplot.savefig("atlas/plots/stage_cluster_unOrder_batch_barplot.pdf",format="pdf", bbox_inches="tight")

tmp = pd.crosstab(adata_concat[adata_concat.obs['sample'] != "lifecyclemix"].obs['sample'],
                  adata_concat[adata_concat.obs['sample'] != "lifecyclemix"].obs['leiden'],
                  normalize='columns').T.plot(kind='bar', stacked=True)
tmp.legend(title='Sample distribution', bbox_to_anchor=(1.26, 1.02),loc='upper right')
plt.pyplot.savefig("atlas/plots/noLifecycleMix_stage_cluster_batch_barplot.pdf",format="pdf", bbox_inches="tight")



tmp = pd.crosstab(adata_concat.obs['leiden'],adata_concat.obs['sample'], normalize='columns').T.plot(kind='bar', stacked=True)
tmp.legend(title='Sample distribution', bbox_to_anchor=(1.26, 1.02),loc='upper right')
plt.pyplot.savefig("atlas/plots/batch_stageCluster_barplot.pdf",format='pdf', bbox_inches="tight")

sc.set_figure_params()


sc.tl.embedding_density(adata_concat, groupby='batch', basis = "umap")
sc.pl.embedding_density(adata_concat, groupby='batch', show = False)
plt.pyplot.savefig("atlas/plots/UMAP_density_batch_plot.pdf", format='pdf', bbox_inches="tight")

sc.set_figure_params(figsize = (12,12))
sc.pl.umap(adata_concat, color = "leiden", legend_loc="on data", show = False)
plt.pyplot.savefig("atlas/plots/UMAP_clusters_plot.pdf", format='pdf', bbox_inches="tight")

sc.set_figure_params(figsize = (12,12))
sc.pl.umap(adata_concat, color = "sample", legend_loc="on data", show = False)
plt.pyplot.savefig("atlas/plots/UMAP_samples_plot.pdf", format='pdf', bbox_inches="tight")


sc.set_figure_params()

sc.pl.umap(adata_concat, color = "gene_counts_cruzi", legend_loc="on data", show = False)
plt.pyplot.savefig("atlas/plots/UMAP_UMI_plot.pdf", format='pdf', bbox_inches="tight")

sc.pl.umap(adata_concat, color = "pct_counts_cruzi", legend_loc="on data", show = False)
plt.pyplot.savefig("atlas/plots/UMAP_percent_cruzi_plot.pdf", format='pdf', bbox_inches="tight")


sc.pl.pca(adata_concat, color = "leiden", legend_loc="on data", show = False)
plt.pyplot.savefig("atlas/plots/PCA_clusters_plot.pdf", format='pdf', bbox_inches="tight")


