library(Seurat)
library(scCATCH)
library(limma)
library(dplyr)
library(Seurat)
library(patchwork)
#----------------- seurat part
library(limma)
library(dplyr)
library(Seurat)
library(patchwork)

hnscc_data <- read.table(file='GSE103322_HNSCC_all_data.txt' , row.names=1, skip=6, sep='\t')

#print(hnscc_data)

hnscc <- CreateSeuratObject(counts=hnscc_data)
print(hnscc)
VlnPlot(hnscc, features=c("nFeature_RNA","nCount_RNA"),ncol=2)
hnscc[["percent.mt"]] <- PercentageFeatureSet(hnscc, pattern = "^MT-")

VlnPlot(hnscc, features=c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3)

#normalization
hnscc_norm <- NormalizeData(hnscc, normalization.method="LogNormalize", scale.factor=10000)
#è già normalizzato il dataset!
#feature selection
hnscc_norm <- FindVariableFeatures(hnscc_norm, selection.method="vst", nfeatures= 5902)
top10 <- head(VariableFeatures(hnscc_norm), 10)
plot1 <-VariableFeaturePlot(hnscc_norm)
plot2 <-LabelPoints(plot= plot1, points=top10, repel =TRUE)
plot1+plot2
hnsccFeatures <- FindVariableFeatures(hnscc, selection.method="vst", nfeatures=5902)
# scaling the data
all.genes <-rownames(hnsccFeatures)

hnsccScaled<-ScaleData(hnscc_norm,features=all.genes)
#perform linear dimensional reduction
hnsccScaled<- RunPCA (hnsccScaled, features = VariableFeatures(object=hnsccScaled))

print(hnsccScaled[["pca"]],dims=1:5, nfeatures=5)
VizDimLoadings(hnsccScaled, dims=1:2, reduction="pca")
DimPlot(hnsccScaled, reduction="pca")
DimHeatmap(hnsccScaled, dims=1,balanced=TRUE)
DimHeatmap(hnsccScaled, dims=1:15,balanced=TRUE)

ElbowPlot(hnsccScaled)

#clustering
hnsccScaled <- FindNeighbors(hnsccScaled, dims = 1:10)
hnsccScaled <- FindClusters(hnsccScaled, resolution = 0.5)
head(Idents(hnsccScaled),5)

#run umap
hnsccScaled <- RunUMAP (hnsccScaled, dims=1:10)
DimPlot(hnsccScaled, reduction="umap")
saveRDS(hnsccScaled, file = "hnscc_final.rds")

# hnscc_data <- read.table(file="GSE103322_HNSCC_all_data.txt" , row.names=1, skip=6, sep='\t')

#print(hnscc_data)
#------------------------
#scCATCH
hnscc=hnsccScaled
clu_markers <- findmarkergenes(object=hnscc, species="Human", cluster="All", match_CellMatch=TRUE, cancer="Head and Neck Cancer", tissue="Oral cavity", cell_min_pct=0.25, logfc=0.25, pvalue=0.05)

hnscc_ann <- scCATCH(object=clu_markers$clu_markers, species="Human", cancer="Head and Neck Cancer", tissue="Oral cavity")
