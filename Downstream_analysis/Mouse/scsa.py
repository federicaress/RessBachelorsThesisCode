

import numpy as np
import pandas as pd
import scanpy as sc

sc.settings.verbosity = 3
sc.logging.print_versions()
#sc.settings.set_figure_params(dpi=80, facecolor='white')
results_file = 'Mouse_results_scanpy.h5ad'
#
# # test= pd.DataFrame(data=pd.read_csv('Docu/Mouse/GSE121861_syngeneic_expression.csv', delimiter=',', dtype ='float32'), index=pd.read_csv('/Docu/Mouse/GSE121861_syngeneic_row_data.csv',dtype='string'))
# test= pd.DataFrame(data=pd.read_csv('~/Docu/Mouse/GSE121861_syngeneic_expression.csv', delimiter=',', dtype ='float32'), index=pd.read_csv('~/Docu/Mouse/GSE121861_syngeneic_row_data.csv',dtype='string'))
# #test= pd.DataFrame(data=pd.read_csv('Docu/Mouse/GSE121861_syngeneic_expression.csv', delimiter=',', dtype ='float32'), index=pd.read_csv('/Docu/Mouse/GSE121861_syngeneic_row_data.csv',dtype='string'), columns=pd.read_csv('/Docu/Mouse/GSE121861_syngeneic_column_data.csv',dtype='string'))
#
# tsynexp=pd.read_csv('~/Docu/Mouse/GSE121861_syngeneic_expression.csv',header = None, index_col=False)
# synrow=pd.read_csv('~/Docu/Mouse/GSE121861_syngeneic_column_data.csv')
# syncol=pd.read_csv('~/Docu/Mouse/GSE121861_syngeneic_row_data.csv')
# test=pd.DataFrame(data=synexp, index=synrow['Sample_ID'], columns=syncol['GeneSymbol'])


dataf = pd.read_csv("~/Docu/Mouse/GSE121861_syngeneic_expression.csv", delimiter=',', header=None, index_col=False)
cells = pd.read_csv('~/Docu/Mouse/GSE121861_syngeneic_column_data.csv', delimiter=',', index_col=False)
genes = pd.read_csv('~/Docu/Mouse/GSE121861_syngeneic_row_data.csv', delimiter=',', index_col=False)
dataf.index = cells['Sample_ID']
dataf.columns = genes['GeneSymbol']
dataf.to_csv("~/Docu/Mouse/mouseexp.csv")

# adata = sc.read('~/Docu/Mouse/mouseexp.csv', delimiter=',')
# adata = sc.read('/Docu/Mouse/GSE121861_syngeneic_expression.csv', delimiter=',', dtype ='float32')

adata = sc.read('mouseexp.csv', delimiter=',')

adata.var_names_make_unique()
adata
adata.obs_names_make_unique()

#preprocessing
sc.pl.highest_expr_genes(adata, n_top=20, )

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

#adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'

sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
# sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
#              jitter=0.4, multi_panel=True)

sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts'],jitter=0.4, multi_panel=True)

# sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')

adata = adata[adata.obs.n_genes < 3000, :]

# adata = adata[adata.obs.total_counts < 10000, :]
#adata = adata[adata.obs.pct_counts_mt < 5, :]

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

sc.settings.verbosity = 2  # reduce the verbosity

sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')

sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)

adata.write(results_file)

sc.tl.rank_genes_groups(adata, 'leiden', method='logreg')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)

#to do: list marker genes
marker_genes=[]


pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)

pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(10)

result = adata.uns['rank_genes_groups'].key()
result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names
dat = pd.DataFrame({group + '_' + key[:1]: result[key][group] for group in groups for key in ['names', 'logfoldchanges','scores','pvals']})
dat.to_csv("scsainput.csv")
