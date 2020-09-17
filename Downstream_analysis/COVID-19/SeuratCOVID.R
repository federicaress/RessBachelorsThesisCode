library(limma)
library(dplyr)
library(Seurat)
library(patchwork)

covid_dataset <- readRDS(file='~/Docu/covid/blish_covid.seu.rds')

print(covid_data)
# pdf(“covidplots.pdf”)
VlnPlot(covid_data, features=c("nFeature_RNA","nCount_RNA"),ncol=2)
# #normalisation
# covid_norm <- NormalizeData(covid_data, normalization.method="LogNormalize", scale.factor=10000)
#
# #feature selection
# covid_norm <- FindVariableFeatures(covid_norm, selection.method="vst", nfeatures= 10473)
# top10 <- head(VariableFeatures(covid_norm), 10)
# plot1 <-VariableFeaturePlot(covid_norm)
# plot2 <-LabelPoints(plot= plot1, points=top10, repel =TRUE)
# plot1+plot2
covidFeatures <- FindVariableFeatures(covid_dataset, selection.method="vst", nfeatures=10473)
# scaling the data
all.genes <-rownames(covidFeatures)

# covidScaled<-ScaleData(covid_norm,features=all.genes)
#perform linear dimensional reduction
# covidScaled<- RunPCA (covidScaled, features = VariableFeatures(object=covidScaled))
#
# print(covidScaled[["pca"]],dims=1:5, nfeatures=5)
# VizDimLoadings(covidScaled, dims=1:2, reduction="pca")
DimPlot(covid_dataset, reduction="pca")
DimHeatmap(covid_dataset, dims=1,balanced=TRUE)
DimHeatmap(covid_dataset, dims=1:15,balanced=TRUE)
#Determine the dimensionality of the dataset
#non si riesce a fare (ucciso)
#covidScaled <- JackStraw (covidScaled, num.replicate = 100)
#covidScaled <- ScoreJackStraw (covidScaled, dims = 1:20)
#

#JackStrawPlot(covidScaled, dims = 1:15)
ElbowPlot(covid_dataset)

# #clustering
# covid_dataset <- FindNeighbors(covid_dataset, dims = 1:10)
# covid_dataset <- FindClusters(covid_dataset, resolution = 0.5)
# head(Idents(covid_dataset),5)
#
# #run umap
# covidScaled <- RunUMAP (covidScaled, dims=1:10)
# DimPlot(covidScaled, reduction="umap")
#
# saveRDS(covidScaled, file= '/home/ressf/Documenti/Tirocinio/singleCellProcessing/COVID/covidScaled.rds')
table(Idents(covid_dataset))

 #markers


# find markers for every cluster compared to all remaining cells, report only the positive ones
covid_dataset.markers <- FindAllMarkers(covid_dataset, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, )

covid_dataset.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
# pdf(“covidplotsumap.pdf”)
DimPlot(covid_dataset,reduction="umap")
# dev.off()
# pdf(“covidplotslabel.pdf”)

new.cluster.ids <- c( "0 Gamma delta T cells","1 T cells ","2 T memory cells", "3 Dendritic cells","4 Basal cells ","5 B cells","6 Dendritic cells","7 Fibroblasts" ,"8 Monocytes","9 Plasma cells","10 EC cells","11 Natural killer cells ","12 T cells ","13 Basal cells","14 Erythroid like cells ","15 Unknown/T cells ","16 Unknown/Inflammated or stressed cells ","17 Basal cells (subcluster)","18 Unknown/B cells (subcluster)","19 Unknown/T cells","20 Dendritic cells","21 Unknown/Basal cells","22 Unknown/T cells","23 Unknown/Basal cells","24 Unknown/T cells","25 Unknown/Dendritic cells","26 Unknown/Endothelial cells","27 Dendritic cells","28 Unknown/Dendritic cells","29 Unknown/Stressed cells " )

names(new.cluster.ids) <- levels(covid_dataset)
covid_dataset <- RenameIdents(covid_dataset, new.cluster.ids)
DimPlot(covid_dataset,reduction="umap")

saveRDS(covid_dataset,file='/home/ressf/Documenti/Tirocinio/singleCellProcessing/COVID/covid_final.rds' )
# dev.off()
