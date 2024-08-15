import scanpy as sc
import numpy as np
import pandas as pd
import os
import seaborn as sns
import matplotlib as plt
import upsetplot
from anndata import AnnData

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

#Removing 8 makes it look great

adata_concat = adata_concat[adata_concat.obs['leiden'].isin([
    "epimast_2", "meta_3",
 "epi_meta_trans_6", "epi_meta_trans_7"])]

#Recluster and go back to X being the normalized counts
adata_trans = adata_concat.raw.to_adata().copy()
adata_trans.layers["counts"] = adata_concat.layers["counts"].copy()
adata_trans.obs = adata_concat.obs.copy()

adata_trans.raw = adata_trans

sc.pp.highly_variable_genes(adata_trans, n_top_genes=1000)
sc.pl.highly_variable_genes(adata_trans)
genes = adata_trans.var.highly_variable

adata_trans.var.highly_variable["MT tcruzi"] = False
adata_trans.var.highly_variable["MT2tcruzi"] = False

sc.pp.regress_out(adata_trans, ['gene_counts_cruzi'])
sc.pp.scale(adata_trans, max_value=10)

sc.pp.pca(adata_trans, svd_solver = "arpack", use_highly_variable = True)
sc.pl.pca_variance_ratio(adata_trans, log=True, n_pcs=50)

print(np.var(adata_trans.obsm["X_pca"], axis = 0))

sc.external.pp.bbknn(adata_trans, batch_key='batch', n_pcs = 9)  # running bbknn 1.3.6

sc.pl.violin(adata_trans, ["total_counts", "C4B63_9g333"], multi_panel = True)

sc.tl.umap(adata_trans)
sc.pl.umap(adata_trans, color = "sample")
sc.pl.umap(adata_trans, color = "batch")

sc.pl.umap(adata_trans, color = "gene_counts_cruzi")


sc.pl.pca(adata_trans, color = "sample")

sc.tl.leiden(adata_trans, resolution = 0.3)

sc.pl.pca(adata_trans, color = "leiden")


sc.pl.umap(adata_trans, color = "leiden")

sc.pl.umap(adata_trans, color = "stages")

#Tubulin--tyrosine ligase-like protein 12
sc.pl.umap(adata_trans, color = "C4B63_2g456")

sc.tl.draw_graph(adata_trans)

sc.tl.paga(adata_trans, groups='leiden')

sc.pl.paga(adata_trans)

sc.tl.draw_graph(adata_trans, init_pos='paga')

sc.pl.draw_graph(adata_trans, color=['leiden'], legend_loc='on data', show = True)
sc.pl.draw_graph(adata_trans, color=['sample'], legend_loc='on data', show = True)
sc.pl.draw_graph(adata_trans, color=['stages'], legend_loc='on data', show = True)

sc.tl.draw_graph(adata_trans, init_pos='paga')

sc.pl.draw_graph(adata_trans, color=['leiden'], legend_loc='on data', show = False)

plt.pyplot.savefig("epi_to_meta/plots/epi_to_meta_paga_embedding_plot.pdf",format='pdf', bbox_inches="tight")

sc.pl.draw_graph(adata_trans, color=['sample'], legend_loc='on data', show = False)
plt.pyplot.savefig("epi_to_meta/plots/epi_to_meta_paga_embedding_sample_plot.pdf",format='pdf', bbox_inches="tight")


sc.tl.rank_genes_groups(adata_trans, 'leiden', method='t-test', use_raw = True)
sc.pl.rank_genes_groups(adata_trans, n_genes=25, sharey=False)

x = pd.DataFrame(adata_trans.uns['rank_genes_groups']['names'])

sc.tl.filter_rank_genes_groups(adata_trans)

filter_genes_pd = pd.DataFrame(adata_trans.uns['rank_genes_groups_filtered']['names'])

gene_names = ["C4B63_2g456", "C4B63_113g30", "C4B63_302g52c", "C4B63_6g311", "C4B63_51g451c",
              "C4B63_106g1", "C4B63_9g333", "C4B63_14g70", "C4B63_85g83", "C4B63_11g67",
            "C4B63_24g199", "C4B63_30g214", "C4B63_18g282", "C4B63_44g125", "C4B63_54g69"
              ]

gene_names = []

for i in range(0,len(filter_genes_pd.keys())):

    current = filter_genes_pd[ filter_genes_pd.keys()[i] ].copy()

    current = current.dropna()

    gene_names.append(current[0:10].tolist())


gene_names = flatten(gene_names)

adata_trans.obs['clusters'] = adata_trans.obs['leiden']  # just a cosmetic change
adata_trans.uns['clusters_colors'] = adata_trans.uns['leiden_colors']

adata_trans.uns['iroot'] = np.flatnonzero(adata_trans.obs['leiden']  == '0')[0]
sc.tl.dpt(adata_trans)

sc.pl.draw_graph(adata_trans, color = ["stages"])
sc.pl.draw_graph(adata_trans, color = ["dpt_pseudotime"])


adata_trans.obs['distance'] = adata_trans.obs['dpt_pseudotime']


new_cluster_names = ["epimastigote_sb", "metacyclic_sb","transitionary_sb"]

adata_trans.rename_categories('leiden', new_cluster_names)

#Change colours
adata_trans.uns["leiden_colors"] = ["#f8766d","#00ba38","#619cff"]

adata_trans.obs['stages_subset'] = adata_trans.obs['leiden']
adata_trans.uns['stages_subset_colors'] = adata_trans.uns['leiden_colors']

sc.pl.draw_graph(adata_trans[adata_trans.obs['distance']<0.65], color=['distance'], show = False)
plt.pyplot.savefig("epi_to_meta/plots/epi_to_meta_paga_embedding_distance_plot.pdf",format='pdf', bbox_inches="tight")

sc.pl.draw_graph(adata_trans[adata_trans.obs['distance']<0.65], color=['stages_subset'], legend_loc = 'on data', show = False)
plt.pyplot.savefig("epi_to_meta/plots/epi_to_meta_paga_embedding_newStages_plot.pdf",format='pdf', bbox_inches="tight")

sc.pl.paga(adata_trans, color=['stages_subset'], threshold = 0.1, show = False)
plt.pyplot.savefig("epi_to_meta/plots/epi_to_meta_paga_edge_plot.pdf",format='pdf', bbox_inches="tight")


sc.pl.rank_genes_groups_heatmap(adata_trans, groupby='stages_subset', n_genes = 10, show = False)
plt.pyplot.savefig("epi_to_meta/plots/heatmap_rank_genes_epi_meta.pdf",format='pdf', bbox_inches="tight")

x.to_csv("epi_to_meta/DEG_tables/top_markers_clusters_epi_to_meta.csv")


x = pd.DataFrame(adata_trans.uns['rank_genes_groups_filtered']['names'])
y = pd.DataFrame(adata_trans.uns['rank_genes_groups_filtered']['pvals_adj'])

x.to_csv("epi_to_meta/DEG_tables/top_markers_clusters_epi_to_meta_filtered.csv")
y.to_csv("epi_to_meta/DEG_tables/top_markers_clusters_epi_to_meta_filtered_adjPvals.csv")

np.savetxt("epi_to_meta/objects/epi_to_meta_expMtx.csv", adata_trans.raw.X.toarray(), delimiter=",")
np.savetxt("epi_to_meta/objects/epi_to_meta_expMtxRaw.csv", adata_trans.layers["counts"].toarray(), delimiter=",", header = "blah")
adata_trans.obs.to_csv("epi_to_meta/objects/epi_to_meta_metadata.csv")
pd.DataFrame(adata_trans.var_names).to_csv("epi_to_meta/objects/epi_to_meta_genes.csv")
