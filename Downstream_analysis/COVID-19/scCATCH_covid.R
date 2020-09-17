library(Seurat)
library(scCATCH)

# per installarlo:

# devtools::install_github('ZJUFanLab/scCATCH')

wilk <- readRDS("blish_covid.seu.rds")

clu_markers <- findmarkergenes(object=wilk, species="Human", cluster="All", match_CellMatch=TRUE, cancer=NULL, tissue="Peripheral blood", cell_min_pct=0.25, logfc=0.25, pvalue=0.05)

clu_ann <- scCATCH(object=clu_markers$clu_markers, species="Human", cancer=NULL, tissue="Peripheral blood")
