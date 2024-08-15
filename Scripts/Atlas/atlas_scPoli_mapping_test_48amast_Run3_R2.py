import scanpy as sc
import numpy as np
import pandas as pd
import os
import seaborn as sns
import matplotlib as plt
from anndata import AnnData
from typing import Sequence, Optional
from sklearn.metrics import classification_report
from scarches.models.scpoli import scPoli

os.chdir("/Users/rosslaidlaw/R/cruzi_paper/")

adata_concat = sc.read_h5ad("atlas/objects/cruzi_bbknn_noAmast_cellCycle.h5ad")
# adata_concat_ind = sc.read_h5ad("cruzi_bbknn_scaleIndividual.h5ad")

adata_concat.obs['leiden'] = adata_concat.obs['leiden'].cat.reorder_categories(["early_mid_amast_5", "late_amast_1",  "trypomast_0", "trypomast_4", "trypomast_9",
 "epimast_2", "epi_meta_trans_6", "epi_meta_trans_7","meta_3", "meta_8"])

adata_concat.obs['lifecycleStage'] = adata_concat.obs['leiden'].copy()

adata_concat.X = adata_concat.raw.X

adata_48 = sc.read_h5ad("/Users/rosslaidlaw/R/cruzi_10x_files/processed_files/late_amast_Run3_R2.h5ad")

adata_48.X = adata_48.layers["counts"].copy()

sc.pp.normalize_total(adata_48, target_sum=1e4)
sc.pp.log1p(adata_48)
adata_48.raw = adata_48

adata_48 = adata_48[:, adata_concat.var.highly_variable]

adata_concat = adata_concat[:, adata_concat.var.highly_variable]

adata_48.obs['batch'] = adata_48.obs['sample']


print("why")

early_stopping_kwargs = {
    "early_stopping_metric": "val_prototype_loss",
    "mode": "min",
    "threshold": 0,
    "patience": 20,
    "reduce_lr": True,
    "lr_patience": 13,
    "lr_factor": 0.1,
}

condition_key = ['batch', 'sample']
cell_type_key = 'lifecycleStage'


scpoli_model = scPoli(
    adata=adata_concat,
    condition_keys=condition_key,
    cell_type_keys=cell_type_key,
    embedding_dims=7,
    recon_loss='nb',
)
scpoli_model.train(
    n_epochs=50,
    pretraining_epochs=40,
    early_stopping_kwargs=early_stopping_kwargs,
    eta=5,
)


adata_48.obs['lifecycleStage'] = ""


scpoli_query = scPoli.load_query_data(
    adata=adata_48,
    reference_model=scpoli_model,
    labeled_indices=[]
)

scpoli_query.train(
    n_epochs=50,
    pretraining_epochs=40,
    eta=10
)


results_dict = scpoli_query.classify(adata_48, scale_uncertainties=True)

for i in range(len(cell_type_key)):
    preds = results_dict[cell_type_key]["preds"]
    results_dict[cell_type_key]["uncert"]
    classification_df = pd.DataFrame(
        classification_report(
            y_true=adata_48.obs[cell_type_key],
            y_pred=preds,
            output_dict=True,
        )
    ).transpose()
print(classification_df)



#get latent representation of reference data
scpoli_query.model.eval()
data_latent_source = scpoli_query.get_latent(
    adata_concat,
    mean=True
)

adata_latent_source = sc.AnnData(data_latent_source)
adata_latent_source.obs = adata_concat.obs.copy()

#get latent representation of query data
data_latent= scpoli_query.get_latent(
    adata_48,
    mean=True
)

#Problem with this section ->
adata_latent = sc.AnnData(data_latent)
adata_latent.obs = adata_48.obs.copy()

#get label annotations
adata_latent.obs['lifecycleStage_pred'] = results_dict['lifecycleStage']['preds'].tolist()
adata_latent.obs['lifecycleStage_uncert'] = results_dict['lifecycleStage']['uncert'].tolist()
adata_latent.obs['classifier_outcome'] = (
    adata_latent.obs['lifecycleStage_pred'] == adata_latent.obs['lifecycleStage']
)
# <- Problem with this section

#get prototypes
labeled_prototypes = scpoli_query.get_prototypes_info()
labeled_prototypes.obs['study'] = 'labeled prototype'
unlabeled_prototypes = scpoli_query.get_prototypes_info(prototype_set='unlabeled')
unlabeled_prototypes.obs['study'] = 'unlabeled prototype'



#join adatas
adata_latent_full = adata_latent_source.concatenate(
    [adata_latent, labeled_prototypes, unlabeled_prototypes],
    batch_key='query'
)
sc.pp.neighbors(adata_latent_full, n_neighbors=15)
# sc.external.pp.bbknn(adata_latent_full, batch_key='batch', n_pcs = 17)
sc.tl.umap(adata_latent_full)

#Save object
adata_latent_full.obs = adata_latent_full.obs.drop(columns = ['classifier_outcome'])
adata_latent_full.write("atlas/objects/cruziAtlas_scPoli_with_Run3_R2.h5ad")
adata_latent_full = sc.read_h5ad("atlas/objects/cruziAtlas_scPoli_with_Run3_R2.h5ad")

#get adata without prototypes
adata_no_prototypes = adata_latent_full[adata_latent_full.obs['query'].isin(['0', '1'])]
sc.pl.umap(
    adata_no_prototypes,
    color='lifecycleStage_pred',
    show=True,
    frameon=False,
)

adata_no_prototypes_query_coloured = adata_no_prototypes.copy()

sc.pl.umap(adata_no_prototypes_query_coloured, color = "lifecycleStage_pred", show = False)
plt.pyplot.savefig("atlas/plots/scPoli_amast48_lifecycleStage.png", bbox_inches="tight")

sc.pl.umap(adata_no_prototypes, color = "lifecycleStage_pred")

test = adata_no_prototypes.copy()
test_other = test[-test.obs['lifecycleStage'].isin(test.obs['lifecycleStage_pred'].cat.categories)].copy()
test.obs['lifecycleStage_pred'] = test.obs['lifecycleStage_pred'].cat.add_categories(test_other.obs['lifecycleStage'][test_other.obs['query'].isin(['0'])].cat.categories)
test.obs['lifecycleStage_pred'][test.obs['query'].isin(['0'])] =test.obs['lifecycleStage'][test.obs['query'].isin(['0'])].tolist()

sc.pl.umap(test, color = "lifecycleStage_pred", show = False)
plt.pyplot.savefig("atlas/plots/scPoli_amast48_lifecycleStage.png", bbox_inches="tight")

sc.pl.umap(test, color = "conditions_combined", show = False)
plt.pyplot.savefig("atlas/plots/scPoli_amast48_batch.png", bbox_inches="tight")


adata_48.obs['cell_type_pred'] = results_dict['lifecycleStage']['preds'].tolist()
adata_48.obs['cell_type_uncert'] = results_dict['lifecycleStage']['uncert'].tolist()
adata_48.obs['classifier_outcome'] = (
    adata_48.obs['cell_type_pred'] == adata_48.obs['lifecycleStage']
)

sc.pp.highly_variable_genes(adata_48)
genes = adata_48.var.highly_variable

adata_48.var.highly_variable["MT tcruzi"] = False
adata_48.var.highly_variable["MT2tcruzi"] = False

#Subsetting out the variable genes gets rid of them in the counts layer.
#Need to keep those genes for tradeSeq

# adata_concat = adata_concat[:,adata_concat.var.highly_variable]

sc.pp.regress_out(adata_48, ['gene_counts_cruzi'])


sc.pp.scale(adata_48, max_value=10)

#By adding the use_highly_variable parameter we don't need to subset the matrix
sc.pp.pca(adata_48, svd_solver = "arpack", use_highly_variable = True)
sc.pl.pca_variance_ratio(adata_48, log=True, n_pcs=50, save='')

sc.pp.neighbors(adata_48, n_neighbors=15)
sc.tl.umap(adata_48)
# sc.pl.umap(adata_48, color = "cell_type_pred")

adata_48.write("/Users/rosslaidlaw/R/cruzi_paper/mapping_test/amast_RunR3_scPoli_map.h5ad")

adata_48 = sc.read_h5ad("/Users/rosslaidlaw/R/cruzi_paper/mapping_test/amast_RunR3_scPoli_map.h5ad")

adata_48.obs['cell_type_pred'].value_counts()

sc.set_figure_params(figsize = (20,10))
sc.pl.violin(adata_48, groupby = 'cell_type_pred', keys = 'cell_type_uncert', show = False)
plt.pyplot.savefig("/Users/rosslaidlaw/R/cruzi_paper/atlas/plots/scPoli_AmastRun3_R2_pred_and_uncert_violin.pdf", bbox_inches="tight",
                   format = 'pdf')


df = pd.DataFrame({'count': adata_48.obs['cell_type_pred'].value_counts()})

#late_amast_1, epi_meta_trans_7, early_mid_amast_5, epimast_2, meta_3, trypomast_0, epi_meta_trans_6, meta_8
colors = ['#ff7f00', '#999999',
          '#377eb8', '#a65628',
          '#bcbd21', '#2aa02b',
          '#e377c2', '#19bdcf']

# Plotting
plt.pyplot.figure(figsize=(60, 40))
fig, ax = plt.pyplot.subplots()
ax.pie(df['count'], labels=df.index, colors = colors)

plt.pyplot.title('Pie chart of cell ID predictions')
plt.pyplot.savefig("/Users/rosslaidlaw/R/cruzi_paper/atlas/plots/scPoli_AmastRun3_R2_pieChart.pdf", bbox_inches="tight", format = 'pdf')

correct_guesses = df[df.index == "early_mid_amast_5"]['count'][0] +  df[df.index == "late_amast_1"]['count'][0]
incorrect_guesses = np.shape(adata_48)[0] - correct_guesses

tp = df[df.index == "early_mid_amast_5"]['count'][0] +  df[df.index == "late_amast_1"]['count'][0]
fp = np.shape(adata_48)[0] - tp
fn = 0
tn = 0

correct_percent = correct_guesses / np.shape(adata_48)[0]

for i in df['count']:

    print( (i / np.shape(adata_48)[0] ) * 100)

print("why")