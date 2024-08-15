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

test = sc.read_h5ad("/Users/rosslaidlaw/R/TrAGEDy_V2/figures/Immune_ALL_human.h5ad")

os.chdir("/Users/rosslaidlaw/R/cruzi_paper/")

adata_concat = sc.read_h5ad("atlas/objects/cruzi_bbknn_noAmast_cellCycle.h5ad")
# adata_concat_ind = sc.read_h5ad("cruzi_bbknn_scaleIndividual.h5ad")

adata_all = sc.read_h5ad("atlas/objects/cruzi_bbknn_noAmast_cellCycle.h5ad")


adata_concat.X = adata_concat.raw.X


adata_concat.obs['leiden'] = adata_concat.obs['leiden'].cat.reorder_categories(["early_mid_amast_5", "late_amast_1",  "trypomast_0", "trypomast_4", "trypomast_9",
 "epimast_2", "epi_meta_trans_6", "epi_meta_trans_7","meta_3", "meta_8"])

adata_concat.obs['lifecycleStage'] = adata_concat.obs['leiden'].copy()

adata_concat.obs["cell_id"] = adata_concat.obs.index

adata_concat = adata_concat[:, adata_concat.var.highly_variable]

adata_concat_cp = adata_concat.copy()

subsample_fractions = [0.05, 0.05,
                       0.25, 0.25,
                       0.5, 0.5,
                       0.75, 0.75,
                       0.9, 0.9]

dim_list = []

for i in range(0,10):

    adata_48 = adata_concat_cp.copy()
    sc.pp.subsample(adata_48, fraction=subsample_fractions[i])

    adata_48.write("atlas/objects/atlasSubsets/cruzi_bbknn_noAmast_cellCycle_raw_scPoliSubset_" + str(i) + "_" + str(subsample_fractions[i]) + ".h5ad")
    adata_48 = sc.read_h5ad("atlas/objects/atlasSubsets/cruzi_bbknn_noAmast_cellCycle_raw_scPoliSubset_" + str(i) + "_" + str(subsample_fractions[i]) + ".h5ad")

    adata_48.obs['lifecycleStage_original'] = adata_48.obs['lifecycleStage']

    adata_48.obs = adata_48.obs.drop(columns="lifecycleStage")

    print("pre-subset")
    print(adata_concat_cp.obs['lifecycleStage'].value_counts())

    adata_concat = adata_concat_cp[-adata_concat_cp.obs['cell_id'].isin(adata_48.obs.index)]

    adata_48.X = adata_48.layers["counts"].copy()

    sc.pp.normalize_total(adata_48, target_sum=1e4)
    sc.pp.log1p(adata_48)
    adata_48.raw = adata_48

    # adata_48.obs['batch'] = adata_48.obs['sample']

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

    print("why must you do this?")
    print(adata_48.obs['lifecycleStage_original'].value_counts())

    print("why must this be?")
    print(adata_concat.obs['lifecycleStage'].value_counts())

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

    for j in range(len(cell_type_key)):
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

    # get latent representation of reference data
    scpoli_query.model.eval()
    data_latent_source = scpoli_query.get_latent(
        adata_concat,
        mean=True
    )

    adata_latent_source = sc.AnnData(data_latent_source)
    adata_latent_source.obs = adata_concat.obs.copy()

    # get latent representation of query data
    data_latent = scpoli_query.get_latent(
        adata_48,
        mean=True
    )

    # Problem with this section ->
    adata_latent = sc.AnnData(data_latent)
    adata_latent.obs = adata_48.obs.copy()

    # get label annotations
    adata_latent.obs['lifecycleStage_pred'] = results_dict['lifecycleStage']['preds'].tolist()
    adata_latent.obs['lifecycleStage_uncert'] = results_dict['lifecycleStage']['uncert'].tolist()
    adata_latent.obs['classifier_outcome'] = (
            adata_latent.obs['lifecycleStage_pred'] == adata_latent.obs['lifecycleStage']
    )
    # <- Problem with this section

    # get prototypes
    labeled_prototypes = scpoli_query.get_prototypes_info()
    labeled_prototypes.obs['study'] = 'labeled prototype'
    unlabeled_prototypes = scpoli_query.get_prototypes_info(prototype_set='unlabeled')
    unlabeled_prototypes.obs['study'] = 'unlabeled prototype'

    # join adatas
    adata_latent_full = adata_latent_source.concatenate(
        [adata_latent, labeled_prototypes, unlabeled_prototypes],
        batch_key='query'
    )
    sc.pp.neighbors(adata_latent_full, n_neighbors=15)
    # sc.external.pp.bbknn(adata_latent_full, batch_key='batch', n_pcs = 17)
    sc.tl.umap(adata_latent_full)

    # Save object
    adata_latent_full.obs = adata_latent_full.obs.drop(columns=['classifier_outcome'])
    # adata_latent_full.write("atlas/objects/cruziAtlas_processed_scPoli_with_atlasSubset_",i,".h5ad")
    # adata_latent_full = sc.read_h5ad("atlas/objects/cruziAtlas_processed_scPoli_with_atlasSubset_",i,".h5ad")

    # get adata without prototypes
    adata_no_prototypes = adata_latent_full[adata_latent_full.obs['query'].isin(['0', '1'])]
    sc.pl.umap(
        adata_no_prototypes,
        color='lifecycleStage_pred',
        show=True,
        frameon=False,
    )

    adata_no_prototypes_query_coloured = adata_no_prototypes.copy()

    sc.pl.umap(adata_no_prototypes_query_coloured, color="lifecycleStage_pred", show=False)
    plt.pyplot.savefig("atlas/plots/atlasSubsets/scPoli_atlasSubset_" + str(i) + "_" + str(subsample_fractions[i]) + "_lifecycleStage.png", bbox_inches="tight")

    sc.pl.umap(adata_no_prototypes, color="lifecycleStage_pred")

    #The following is to plot the the embeddings of the query subset data in relation to the rest of the training dataset
    #It is not needed for the analysis
    test = adata_no_prototypes.copy()
    test_other = test[-test.obs['lifecycleStage'].isin(test.obs['lifecycleStage_pred'].cat.categories)].copy()
    test.obs['lifecycleStage_pred'] = test.obs['lifecycleStage_pred'].cat.add_categories(
        test_other.obs['lifecycleStage'][test_other.obs['query'].isin(['0'])].cat.categories)
    test.obs['lifecycleStage_pred'][test.obs['query'].isin(['0'])] = test.obs['lifecycleStage'][
        test.obs['query'].isin(['0'])].tolist()

    sc.pl.umap(test, color="lifecycleStage_pred", show=False)
    plt.pyplot.savefig("atlas/plots/atlasSubsets/scPoli_atlasSubset_" + str(i) + "_" + str(subsample_fractions[i]) + "_lifecycleStage.png", bbox_inches="tight")

    sc.pl.umap(test, color="conditions_combined", show=False)
    plt.pyplot.savefig("atlas/plots/atlasSubsets/scPoli_atlasSubset_" + str(i) + "_" + str(subsample_fractions[i]) + "_batch.png", bbox_inches="tight")

    sc.pl.umap(test, color="conditions_combined")

    adata_48.obs['cell_type_pred'] = results_dict['lifecycleStage']['preds'].tolist()
    adata_48.obs['cell_type_uncert'] = results_dict['lifecycleStage']['uncert'].tolist()
    adata_48.obs['classifier_outcome'] = (
            adata_48.obs['cell_type_pred'] == adata_48.obs['lifecycleStage']
    )



    adata_48.write("/Users/rosslaidlaw/R/cruzi_paper/atlas/objects/atlasSubsets/amast_atlasSubset_" + str(i) + "_" + str(subsample_fractions[i]) +"_processed_scPoli_map.h5ad")

    adata_48.obs.to_csv("/Users/rosslaidlaw/R/cruzi_paper/atlas/atlasSubset_metadata/atlas_subset"+ str(i) + "_" + str(subsample_fractions[i]) +"_scPoli_metadata.csv")

    adata_48 = sc.read_h5ad("/Users/rosslaidlaw/R/cruzi_paper/atlas/objects/atlasSubsets/amast_atlasSubset_" + str(i) + "_" + str(subsample_fractions[i]) + "_processed_scPoli_map.h5ad")

    sc.set_figure_params(figsize=(20, 10))
    sc.pl.violin(adata_48, groupby='cell_type_pred', keys='cell_type_uncert', show=False)
    plt.pyplot.savefig("/Users/rosslaidlaw/R/cruzi_paper/atlas/plots/atlasSubsets/scPoli_all_pred_and_uncert_violin.png",
                       bbox_inches="tight")

    df = pd.DataFrame({'count': adata_48.obs['cell_type_pred'].value_counts()})

    # Need 10 colours
    colors = ['red', 'blue', 'purple', 'orange', 'yellow', "green", "pink", "black", "cyan"]

    # Plotting
    plt.pyplot.figure(figsize=(8, 6))
    df.transpose().plot(kind='bar', stacked=True, color=colors)
    plt.pyplot.xlabel('Categories')
    plt.pyplot.ylabel('Count')
    plt.pyplot.title('Stacked Bar Plot')
    plt.pyplot.legend(title='Labels')
    plt.pyplot.savefig("/Users/rosslaidlaw/R/cruzi_paper/atlas/plots/atlasSubsets/scPoliSubset_" + str(i) + "_" + str(subsample_fractions[i]) + "_all_pred_and_uncert_stackedBarplot.png",
                       bbox_inches="tight")
    print(subsample_fractions[i])

    dim_list.append(np.shape(adata_48)[0])


print("why")