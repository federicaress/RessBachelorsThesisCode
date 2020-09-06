#scanpy HNSCC.py
import numpy as np
import pandas as pd
import scanpy as sc

sc.settings.verbosity = 3
sc.logging.print_versions()
results_file = '/home/ressf/Documenti/RessBachelorsThesisCode/Downstream_analysis/HNSCC/results_scanpy.h5ad'


adata = sc.read_text('/home/ressf/Documenti/RessBachelorsThesisCode/Downstream_analysis/HNSCC/hnscc_clean_trasp.txt', delimiter='\t', dtype ='float32')


adata.var_names_make_unique()
adata

#preprocessing
sc.pl.highest_expr_genes(adata, n_top=20, )

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)

sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts'],
             jitter=0.4, multi_panel=True)

sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')

adata = adata[adata.obs.n_genes_by_counts < 9000, :]

sc.pp.normalize_total(adata, target_sum=1e4)

sc.pp.log1p(adata)

sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

sc.pl.highly_variable_genes(adata)

adata.raw = adata


#filtering

adata = adata[:, adata.var.highly_variable]

sc.pp.regress_out(adata, ['total_counts'])
sc.pp.scale(adata, max_value=10)



#PCA

sc.tl.pca(adata, svd_solver='arpack')

sc.pl.pca(adata)

sc.pl.pca_variance_ratio(adata, log=True)
adata.write(results_file)
adata

#nn
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)
sc.pl.umap(adata)
sc.pl.umap(adata, use_raw=False)


#clustering
sc.tl.leiden(adata)
sc.pl.umap(adata, color=['leiden'])
adata.write(results_file)

#find marker genes
sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)

sc.settings.verbosity = 2 

sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')

sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)

adata.write(results_file)

sc.tl.rank_genes_groups(adata, 'leiden', method='logreg')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)

pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(10)
cluster =0
for j in range(0,31):
    for i in range(0,5869):
        if int(adata.obs.leiden[i]) == j:
            cluster = cluster+1
    print("Cluster number ")
    print(j)
    print("contains # cells = ")
    print(cluster)
    cluster=0

#mark the cell types
cnames = ["0 Basal cells",
"1 Fibroblasts",
"2 EC cells",
 "3 T cells",
"4 Basal cells",
"5 Fibroblasts",
 "6 Basal cells",
"7 Gamma delta T cells",
"8 Dendritic cells", "9 Natural killer cells",
"10 B cells", "11 EC cells",
"12 Basal cells",
"13 Basal cells",
"14 T cells", "15 T cells",
"16 Basal cells",
"17 Mast cells", "18 Monocytes", "19 Luminal epithelial cells",
"20 Fibroblasts",
"21 Smooth muscle cells",
"22 Basal cells",
"23 T cels",
"24 Basal cells",
"25 T memory cells",
"26 Dendritic cells",
"27 Basal cells",
"28 EC cells",
"29 Gamma delta T cells",
"30 Smooth muscle cells"]
adata.rename_categories('leiden', cnames)

sc.pl.umap(adata, color='leiden', legend_loc='on data', title='', frameon=False, save='.pdf')
