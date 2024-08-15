import scanpy as sc
import numpy as np
import pandas as pd
import os
import seaborn as sns
import matplotlib as plt
import bbknn

def flatten(l):
    return [item for sublist in l for item in sublist]

os.chdir("/Users/rosslaidlaw/R/cruzi_10x_files/")

datasets = ["amast","late_amast_Run3_R2","late_amast_Run3_R3","epimast", "lifecycleMix1", "lifecycleMix2", "meta1", "meta2", "trypomast"]

#Early mid amasts
adata = sc.read_10x_mtx("10x_files/amast/",
                        var_names='gene_symbols',
                        cache=True)


# sc.pp.filter_genes(adata, min_cells=10)

adata.var['mt'] = adata.var_names.str.startswith('TcruziDm28cUTR_MT_MT')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
# Human QC
adata.var['human'] = adata.var_names.str.startswith('Hsampien')
sc.pp.calculate_qc_metrics(adata, qc_vars=['human'], percent_top=None, log1p=False, inplace=True)

# Cruzi QC
adata.var['cruzi'] = adata.var_names.str.startswith('TcruziDm28cUTR')
sc.pp.calculate_qc_metrics(adata, qc_vars=['cruzi'], percent_top=None, log1p=False, inplace=True)


# nFeature_cruzi
adata.obs['gene_counts_cruzi'] = flatten(np.sum(adata.X[:, adata.var_names.str.startswith('TcruziDm28cUTR')] > 0, axis=1).tolist())

sc.pl.violin(adata, ['gene_counts_cruzi',
                     'pct_counts_mt',
                     'pct_counts_human'],jitter=0.4, multi_panel=True,
             show=False)
plt.pyplot.savefig("/Users/rosslaidlaw/R/cruzi_paper/QC_plots/early_mid_amast.pdf",format='pdf', bbox_inches="tight")



# subset out cruzi genes
adata = adata[:, adata.var_names.str.startswith('TcruziDm28cUTR')]

adata = adata[adata.obs.gene_counts_cruzi > 200, :]
adata = adata[adata.obs.gene_counts_cruzi < 750, :]
adata = adata[adata.obs.pct_counts_mt < 20, :]

adata = adata[adata.obs.pct_counts_human < 50, :]

adata.var_names = adata.var_names.str.replace("TcruziDm28cUTR_MT_", "")

#Remove MT genes
# adata = adata[ : , adata.var_names[ adata.var_names.str.startswith('C4B') ] ]

adata.obs['sample'] = "amast_6hr_24hr"

adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

adata.write("processed_files/early_mid_amast.h5ad")

#48 hour amast Run3_R2
adata = sc.read_10x_mtx("10x_files/late_amast_Run3_R2/",
                        var_names='gene_symbols',
                        cache=True)


# sc.pp.filter_genes(adata, min_cells=10)

adata.var['mt'] = adata.var_names.str.startswith('TcruziDm28cUTR_MT_MT')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
# Human QC
adata.var['human'] = adata.var_names.str.startswith('Hsampien')
sc.pp.calculate_qc_metrics(adata, qc_vars=['human'], percent_top=None, log1p=False, inplace=True)

# Cruzi QC
adata.var['cruzi'] = adata.var_names.str.startswith('TcruziDm28cUTR')
sc.pp.calculate_qc_metrics(adata, qc_vars=['cruzi'], percent_top=None, log1p=False, inplace=True)

# nFeature_cruzi
adata.obs['gene_counts_cruzi'] = flatten(np.sum(adata.X[:, adata.var_names.str.startswith('TcruziDm28cUTR')] > 0, axis=1).tolist())

# sc.pl.violin(adata, ['pct_counts_human'],jitter=0.4, return_fig = True)
# sc.pl.violin(adata, ['gene_counts_cruzi'],jitter=0.4, multi_panel=True)
# sc.pl.violin(adata, ['pct_counts_mt'],jitter=0.4, multi_panel=True)

# subset out cruzi genes
adata = adata[:, adata.var_names.str.startswith('TcruziDm28cUTR')]

adata = adata[adata.obs.gene_counts_cruzi > 200, :]
adata = adata[adata.obs.gene_counts_cruzi < 1500, :]
adata = adata[adata.obs.pct_counts_mt < 25, :]

adata = adata[adata.obs.pct_counts_human < 25, :]

adata.var_names = adata.var_names.str.replace("TcruziDm28cUTR_MT_", "")

#Remove MT genes
# adata = adata[ : , adata.var_names[ adata.var_names.str.startswith('C4B') ] ]
adata.obs['sample'] = "48hr_amast"

adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

adata.write("processed_files/late_amast_Run3_R2.h5ad")


#48 hour amast Run3_R3
adata = sc.read_10x_mtx("10x_files/late_amast_Run3_R3/",
                        var_names='gene_symbols',
                        cache=True)


# sc.pp.filter_genes(adata, min_cells=10)

adata.var['mt'] = adata.var_names.str.startswith('TcruziDm28cUTR_MT_MT')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
# Human QC
adata.var['human'] = adata.var_names.str.startswith('Hsampien')
sc.pp.calculate_qc_metrics(adata, qc_vars=['human'], percent_top=None, log1p=False, inplace=True)

# Cruzi QC
adata.var['cruzi'] = adata.var_names.str.startswith('TcruziDm28cUTR')
sc.pp.calculate_qc_metrics(adata, qc_vars=['cruzi'], percent_top=None, log1p=False, inplace=True)

# nFeature_cruzi
adata.obs['gene_counts_cruzi'] = flatten(np.sum(adata.X[:, adata.var_names.str.startswith('TcruziDm28cUTR')] > 0, axis=1).tolist())

# sc.pl.violin(adata, ['pct_counts_human'],jitter=0.4, return_fig = True)
# sc.pl.violin(adata, ['gene_counts_cruzi'],jitter=0.4, multi_panel=True)
# sc.pl.violin(adata, ['pct_counts_mt'],jitter=0.4, multi_panel=True)

# subset out cruzi genes
adata = adata[:, adata.var_names.str.startswith('TcruziDm28cUTR')]

adata = adata[adata.obs.gene_counts_cruzi > 200, :]
adata = adata[adata.obs.gene_counts_cruzi < 1750, :]
adata = adata[adata.obs.pct_counts_mt < 33, :]

adata = adata[adata.obs.pct_counts_human < 25, :]

adata.var_names = adata.var_names.str.replace("TcruziDm28cUTR_MT_", "")

#Remove MT genes
# adata = adata[ : , adata.var_names[ adata.var_names.str.startswith('C4B') ] ]
adata.obs['sample'] = "48hr_amast"

adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

adata.write("processed_files/late_amast_Run3_R3.h5ad")


#Epimast
adata = sc.read_10x_mtx("10x_files/epimast/",
                        var_names='gene_symbols',
                        cache=True)


# sc.pp.filter_genes(adata, min_cells=10)

adata.var['mt'] = adata.var_names.str.startswith('TcruziDm28cUTR_MT_MT')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
# Human QC
adata.var['human'] = adata.var_names.str.startswith('Hsampien')
sc.pp.calculate_qc_metrics(adata, qc_vars=['human'], percent_top=None, log1p=False, inplace=True)

# Cruzi QC
adata.var['cruzi'] = adata.var_names.str.startswith('TcruziDm28cUTR')
sc.pp.calculate_qc_metrics(adata, qc_vars=['cruzi'], percent_top=None, log1p=False, inplace=True)

# nFeature_cruzi
adata.obs['gene_counts_cruzi'] = flatten(np.sum(adata.X[:, adata.var_names.str.startswith('TcruziDm28cUTR')] > 0, axis=1).tolist())

# sc.pl.violin(adata, ['pct_counts_human'],jitter=0.4, return_fig = True)
# sc.pl.violin(adata, ['gene_counts_cruzi'],jitter=0.4, multi_panel=True)
# sc.pl.violin(adata, ['pct_counts_mt'],jitter=0.4, multi_panel=True)

sc.pl.violin(adata, ['gene_counts_cruzi',
                     'pct_counts_mt',
                     'pct_counts_human'],jitter=0.4, multi_panel=True,
             show=False)
plt.pyplot.savefig("/Users/rosslaidlaw/R/cruzi_paper/QC_plots/epi_lateAmast.pdf",format='pdf', bbox_inches="tight")

# subset out cruzi genes
adata = adata[:, adata.var_names.str.startswith('TcruziDm28cUTR')]

adata = adata[adata.obs.gene_counts_cruzi > 400, :]
adata = adata[adata.obs.gene_counts_cruzi < 2000, :]
adata = adata[adata.obs.pct_counts_mt < 30, :]

adata = adata[adata.obs.pct_counts_human < 10, :]

adata.var_names = adata.var_names.str.replace("TcruziDm28cUTR_MT_", "")

#Remove MT genes
# adata = adata[ : , adata.var_names[ adata.var_names.str.startswith('C4B') ] ]
adata.obs['sample'] = "epimast_amastDay5"

adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

adata.write("processed_files/epimast_LateAmast.h5ad")


#Lifecyclemix 1
#Increasing the percentage human higher than should be to see if we can capture more early/mid amastigotes

adata = sc.read_10x_mtx("10x_files/lifecycleMix1/",
                        var_names='gene_symbols',
                        cache=True)


# sc.pp.filter_genes(adata, min_cells=10)

adata.var['mt'] = adata.var_names.str.startswith('TcruziDm28cUTR_MT_MT')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
# Human QC
adata.var['human'] = adata.var_names.str.startswith('Hsampien')
sc.pp.calculate_qc_metrics(adata, qc_vars=['human'], percent_top=None, log1p=False, inplace=True)

# Cruzi QC
adata.var['cruzi'] = adata.var_names.str.startswith('TcruziDm28cUTR')
sc.pp.calculate_qc_metrics(adata, qc_vars=['cruzi'], percent_top=None, log1p=False, inplace=True)

# nFeature_cruzi
adata.obs['gene_counts_cruzi'] = flatten(np.sum(adata.X[:, adata.var_names.str.startswith('TcruziDm28cUTR')] > 0, axis=1).tolist())

# sc.pl.violin(adata, ['pct_counts_human'],jitter=0.4, return_fig = True)
# sc.pl.violin(adata, ['gene_counts_cruzi'],jitter=0.4, multi_panel=True)
# sc.pl.violin(adata, ['pct_counts_mt'],jitter=0.4, multi_panel=True)

sc.pl.violin(adata, ['gene_counts_cruzi',
                     'pct_counts_mt',
                     'pct_counts_human'],jitter=0.4, multi_panel=True,
             show=False)
plt.pyplot.savefig("/Users/rosslaidlaw/R/cruzi_paper/QC_plots/MIX1.pdf",format='pdf', bbox_inches="tight")

# subset out cruzi genes
adata = adata[:, adata.var_names.str.startswith('TcruziDm28cUTR')]

adata = adata[adata.obs.gene_counts_cruzi > 150, :]
adata = adata[adata.obs.gene_counts_cruzi < 1500, :]
adata = adata[adata.obs.pct_counts_mt < 28, :]

adata = adata[adata.obs.pct_counts_human < 30, :]

adata.var_names = adata.var_names.str.replace("TcruziDm28cUTR_MT_", "")

#Remove MT genes
# adata = adata[ : , adata.var_names[ adata.var_names.str.startswith('C4B') ] ]
adata.obs['sample'] = "lifecyclemix"

adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

adata.write("processed_files/lifecycleMix1.h5ad")



#Lifecycle mix 2
adata = sc.read_10x_mtx("10x_files/lifecycleMix2/",
                        var_names='gene_symbols',
                        cache=True)


# sc.pp.filter_genes(adata, min_cells=10)

adata.var['mt'] = adata.var_names.str.startswith('TcruziDm28cUTR_MT_MT')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
# Human QC
adata.var['human'] = adata.var_names.str.startswith('Hsampien')
sc.pp.calculate_qc_metrics(adata, qc_vars=['human'], percent_top=None, log1p=False, inplace=True)

# Cruzi QC
adata.var['cruzi'] = adata.var_names.str.startswith('TcruziDm28cUTR')
sc.pp.calculate_qc_metrics(adata, qc_vars=['cruzi'], percent_top=None, log1p=False, inplace=True)

# nFeature_cruzi
adata.obs['gene_counts_cruzi'] = flatten(np.sum(adata.X[:, adata.var_names.str.startswith('TcruziDm28cUTR')] > 0, axis=1).tolist())

# sc.pl.violin(adata, ['pct_counts_human'],jitter=0.4, return_fig = True)
# sc.pl.violin(adata, ['gene_counts_cruzi'],jitter=0.4, multi_panel=True)
# sc.pl.violin(adata, ['pct_counts_mt'],jitter=0.4, multi_panel=True)

sc.pl.violin(adata, ['gene_counts_cruzi',
                     'pct_counts_mt',
                     'pct_counts_human'],jitter=0.4, multi_panel=True,
             show=False)
plt.pyplot.savefig("/Users/rosslaidlaw/R/cruzi_paper/QC_plots/MIX2.pdf",format='pdf', bbox_inches="tight")

# subset out cruzi genes
adata = adata[:, adata.var_names.str.startswith('TcruziDm28cUTR')]

adata = adata[adata.obs.gene_counts_cruzi > 200, :]
adata = adata[adata.obs.gene_counts_cruzi < 1250, :]
adata = adata[adata.obs.pct_counts_mt < 30, :]

adata = adata[adata.obs.pct_counts_human < 30, :]

adata.var_names = adata.var_names.str.replace("TcruziDm28cUTR_MT_", "")


#Remove MT genes
# adata = adata[ : , adata.var_names[ adata.var_names.str.startswith('C4B') ] ]
adata.obs['sample'] = "lifecyclemix"

adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

adata.write("processed_files/lifecycleMix2.h5ad")



#Metacyclic 1
adata = sc.read_10x_mtx("10x_files/meta1/",
                        var_names='gene_symbols',
                        cache=True)


# sc.pp.filter_genes(adata, min_cells=10)

adata.var['mt'] = adata.var_names.str.startswith('MT')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
adata.obs['gene_counts_cruzi'] = adata.obs["n_genes_by_counts"]
adata.obs['pct_counts_human'] = 0

#Get rid of the massive outliers
adata = adata[adata.obs.n_genes_by_counts < 3000, :]


# sc.pl.violin(adata, ['gene_counts_cruzi'],jitter=0.4, multi_panel=True)
# sc.pl.violin(adata, ['pct_counts_mt'],jitter=0.4, multi_panel=True)

sc.pl.violin(adata, ['gene_counts_cruzi',
                     'pct_counts_mt'],jitter=0.4, multi_panel=True,
             show=False)
plt.pyplot.savefig("/Users/rosslaidlaw/R/cruzi_paper/QC_plots/meta1.pdf",format='pdf', bbox_inches="tight")

adata = adata[adata.obs.n_genes_by_counts > 150, :]
adata = adata[adata.obs.n_genes_by_counts < 1000, :]
adata = adata[adata.obs.pct_counts_mt < 40, :]

adata.obs['sample'] = "meta"

#Remove MT genes
# adata = adata[ : , adata.var_names[ adata.var_names.str.startswith('C4B') ] ]

adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

adata.write("processed_files/meta1.h5ad")



#Metacyclic 2
adata = sc.read_10x_mtx("10x_files/meta2/",
                        var_names='gene_symbols',
                        cache=True)


# sc.pp.filter_genes(adata, min_cells=10)

adata.var['mt'] = adata.var_names.str.startswith('MT')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
adata.obs['gene_counts_cruzi'] = adata.obs["n_genes_by_counts"]
adata.obs['pct_counts_human'] = 0

#Get rid of the massive outliers
adata = adata[adata.obs.n_genes_by_counts < 3000, :]


# sc.pl.violin(adata, ['gene_counts_cruzi'],jitter=0.4, multi_panel=True)
# sc.pl.violin(adata, ['pct_counts_mt'],jitter=0.4, multi_panel=True)

sc.pl.violin(adata, ['gene_counts_cruzi',
                     'pct_counts_mt'],jitter=0.4, multi_panel=True,
             show=False)
plt.pyplot.savefig("/Users/rosslaidlaw/R/cruzi_paper/QC_plots/meta2.pdf",format='pdf', bbox_inches="tight")

adata = adata[adata.obs.n_genes_by_counts > 200, :]
adata = adata[adata.obs.n_genes_by_counts < 1250, :]
adata = adata[adata.obs.pct_counts_mt < 40, :]
adata.obs['sample'] = "meta"

#Remove MT genes
# adata = adata[ : , adata.var_names[ adata.var_names.str.startswith('C4B') ] ]

adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

adata.write("processed_files/meta2.h5ad")

#Trypomast
adata = sc.read_10x_mtx("10x_files/trypomast/",
                        var_names='gene_symbols',
                        cache=True)


# sc.pp.filter_genes(adata, min_cells=10)

adata.var['mt'] = adata.var_names.str.startswith('MT')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
adata.obs['gene_counts_cruzi'] = adata.obs["n_genes_by_counts"]
adata.obs['pct_counts_human'] = 0


# sc.pl.violin(adata, ['gene_counts_cruzi'],jitter=0.4, multi_panel=True)
# sc.pl.violin(adata, ['pct_counts_mt'],jitter=0.4, multi_panel=True)

sc.pl.violin(adata, ['gene_counts_cruzi',
                     'pct_counts_mt'],jitter=0.4, multi_panel=True,
             show=False)
plt.pyplot.savefig("/Users/rosslaidlaw/R/cruzi_paper/QC_plots/Trypomast.pdf",format='pdf', bbox_inches="tight")

adata = adata[adata.obs.n_genes_by_counts > 200, :]
adata = adata[adata.obs.n_genes_by_counts < 1500, :]
adata = adata[adata.obs.pct_counts_mt < 25, :]

adata.obs['sample'] = "trypomast"

#Remove MT genes
# adata = adata[ : , adata.var_names[ adata.var_names.str.startswith('C4B') ] ]

adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

adata.write("processed_files/trypomast.h5ad")


