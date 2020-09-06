library(limma)
library(dplyr)
library(Seurat)
library(patchwork)

hnscc_data <- read.table(file='/home/ressf/Documenti/RessBachelorsThesisCode/Downstream_analysis/HNSCC/HNSCC_all_data.txt' , row.names=1, skip=6, sep='\t')

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
#Determine the dimensionality of the dataset
#non si riesce a fare (ucciso)
#hnsccScaled <- JackStraw (hnsccScaled, num.replicate = 100)
#hnsccScaled <- ScoreJackStraw (hnsccScaled, dims = 1:20)
#

#JackStrawPlot(hnsccScaled, dims = 1:15)
ElbowPlot(hnsccScaled)

#clustering
hnsccScaled <- FindNeighbors(hnsccScaled, dims = 1:10)
hnsccScaled <- FindClusters(hnsccScaled, resolution = 0.5)
head(Idents(hnsccScaled),5)

#run umap
hnsccScaled <- RunUMAP (hnsccScaled, dims=1:10)
DimPlot(hnsccScaled, reduction="umap")

saveRDS(hnsccScaled, file= '/home/ressf/Documenti/Tirocinio/singleCellProcessing/HNSCC/hnsccScaled.rds')
table(Idents(hnsccScaled))

 #markers
cluster1.markers <- FindMarkers(hnsccScaled, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)
cluster5.markers <- FindMarkers(hnsccScaled, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive ones
hnsccScaled.markers <- FindAllMarkers(hnsccScaled, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, )
hnsccScaled.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

# VlnPlot(hnsccScaled, features = c("S100A2", "S100A9"))
# VlnPlot(hnsccScaled, features = c("CXCR4", "PTPRC"), slot = "counts", log = TRUE)
# FeaturePlot(hnsccScaled, features=c("S100A2", "S100A9","CXCR4", "PTPRC","ACTA2","TAGLN","EPCAM","GSTA1","COL3A1","LUM"))
# top10 <- hnsccScaled.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
# DoHeatmap(hnsccScaled, features = top10$gene) + NoLegend()




plotmarker0<- FeaturePlot(hnsccScaled, features=c("S100A2", "S100A9", "KRT16", "ANXA8L1", "KRT14", "S100A8", "GJB2", "COL17A1", "KRT6A", "LGALS7"))

plotmarker1<- FeaturePlot(hnsccScaled, features= c("CXCR4", "PTPRC", "CD2", "CD3E", "LCK", "CD3D", "FYN", "SRGN", "GZMK", "CCR7"))

plotmarker2<-FeaturePlot(hnsccScaled, features=c("ACTA2", "TAGLN", "MUSTN1", "PLN", "A2M", "SPARCL1", "MYL9", "ADAMTS1", "MYH11", "FHL5"))

plotmarker3<-FeaturePlot(hnsccScaled, features=c("EPCAM", "GSTA1", "GSTM1", "CES1", "SPP1", "FGFBP2", "OSGIN1", "ALDH3A1", "GPX2", "UGT1A6"))

plotmarker4<-FeaturePlot(hnsccScaled, features=c("COL3A1", "COL6A3", "RARRES2", "MXRA8", "COL1A2", "LUM", "DCN", "CTSK", "COL1A1", "MMP2"))

plotmarker5<-FeaturePlot(hnsccScaled,features=c("CXCR6", "CD3D", "SIRPG", "CD2", "PRF1", "TIGIT", "ICOS", "LCK", "GNLY", "GZMB"))

plotmarker6<-FeaturePlot(hnsccScaled,features=c("CCNB1", "GJB6", "TK1", "PLEK2", "CDC20", "GJB2", "COL17A1", "FGFBP1", "S100A9", "S100A2"))

plotmarker7<-FeaturePlot(hnsccScaled,features=c("ELK2AP", "MZB1", "DERL3", "CD79A", "IGLL5", "IGJ", "NCF1", "BLNK", "SLAMF7", "FKBP11"))

plotmarker8<-FeaturePlot(hnsccScaled,features=c("PLVAP", "RAMP3", "RNASE1", "CLEC14A", "DARC", "AQP1", "VWF", "CDH5", "CLDN5", "HYAL2"))

plotmarker9<-FeaturePlot(hnsccScaled,features=c("CSN3", "FDCSP", "UBD", "TRIT1", "MGST1", "MYCL1", "ODAM", "SULT1E1", "CLDN1", "KRT19"))

plotmarker10<-FeaturePlot(hnsccScaled,features=c("APOD", "SFRP2", "DPT", "PODN", "ABCA8", "FBLN2", "MFAP4", "C3", "CXCL12", "CFD"))

plotmarker11<-FeaturePlot(hnsccScaled,features=c("GPX2", "AKR1C2", "AKR1C1", "ALDH3A1", "ADH7", "IPMK", "AKR1C3", "SAMD12", "SUZ12P1", "CYP4F11"))

plotmarker12<-FeaturePlot(hnsccScaled,features=c("AZGP1", "ZG16B", "CRISP3", "TCN1", "PIP", "STATH", "MUC7", "PRR4", "SLPI", "PIGR" ))

plotmarker13<-FeaturePlot(hnsccScaled,features=c("TPSAB1", "TPSB2", "CPA3", "TPSD1", "HDC", "CTSG", "KIT", "HPGDS", "IL1RL1", "MS4A2"))

plotmarker14<-FeaturePlot(hnsccScaled,features=c("C1QB", "AIF1", "C1QC", "MS4A6A", "C1QA", "FCER1G", "FCGR3A", "TYROBP", "IGSF6", "FCGR2A" ))


new.cluster.ids <-c("Basal cells", "T helper cells","Ectoderm","Basal cells, maybe cancer basal cells", "Fibroblast", "T cells", "Basal cells in phase S/M", "B cells", "Endothelial cells", "Mesoderm cells", "Fibroblast", "Aorta Endothelial cell", "Airway epithelial cell", "mast cell + neutrophils", "gamma delta T cells")
names(new.cluster.ids) <- levels(hnsccScaled)
hnsccScaled <- RenameIdents(hnsccScaled, new.cluster.ids)
DimPlot(hnsccScaled, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

saveRDS(hnsccScaled, file = "/home/ressf/Documenti/Tirocinio/singleCellProcessing/HNSCC/hnscc_final.rds")
