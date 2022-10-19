#clean evironment
remove(list = ls())
#required packages
#Bioconductor install install.packages("BiocManager")
#BiocManager::install("package name")
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(monocle3)
library(SeuratWrappers)
library(SingleCellExperiment)
library(tidyverse)
library(metap)
library(tibble)
library("devtools")
library(infercnv)
#########s2 clean########
s2.data = Read10X(data.dir = "/Users/joshuarivera/OneDrive - University of Massachusetts Boston/Shared_SPJR/scRNA-seq/raw data/s2_dfci/filtered_feature_bc_matrix/")
str(s2.data)

s2 = new("seurat", raw.data = s2.data)
str(s2)
class(s2)
s2 = CreateSeuratObject(counts = s2.data, min.cells= 5, min.features = 200, project = "midpoint_treated")


saveRDS(s2, file= "/Users/joshuarivera/s2_treated.rds")
readRDS(s2, file= "/Users/joshuarivera/s2_treated.rds")

s2[["percent.mt"]]= PercentageFeatureSet(s2, pattern = "^mt-")
s2[["percent.mt"]]

VlnPlot(s2, features = c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3)

s2.plot1 = FeatureScatter(s2, feature1="nCount_RNA", feature2 = "percent.mt",smooth = TRUE,span = TRUE,col="linen")
s2.plot1
s2.plot2 = FeatureScatter(s2, feature1="nCount_RNA",feature2= "nFeature_RNA",smooth=TRUE,span=TRUE, col="linen")
s2.plot2      
CombinePlots(list(s2.plot1,s2.plot2))  


s2 = subset(s2, subset = nFeature_RNA>200 & nFeature_RNA<3750 & percent.mt <5)
s2 = NormalizeData(s2, normalization.method = "LogNormalize", scale.factor = 1e4)


s2 = FindVariableFeatures(s2, selection.method = "vst", nfeatures = 3000)
s2.top10 = head(VariableFeatures(s2),10)
s2.top10
s2.top20 = head(VariableFeatures(s2),20)
s2.top20
s2.top50 = head(VariableFeatures(s2),50)
s2.top50
s2.plot3 = VariableFeaturePlot(s2)

s2.top10.p = LabelPoints(plot=s2.plot3, points = s2.top10, repel=TRUE)
s2.top20.p = LabelPoints(plot=s2.plot3, points = s2.top20, repel=TRUE)
s2.top50.p = LabelPoints(plot=s2.plot3, points = s2.top50, repel=TRUE)

plot(s2.top20.p)
CombinePlots(list(s2.top20.p,s1.top20.p))


####Scale and Principal component analysis
s2.all.genes = row.names(s2)
s2 = ScaleData(s2, features=s2.all.genes)
s2 = RunPCA(s2, features = VariableFeatures(s2), genes.print=10)

#####visualization of data####
VizDimLoadings(s2,dims=1:7, reduction="pca",col="rosybrown")
 DimPlot(s2, reduction="pca", dims = 3:4)

s2 = JackStraw(s2, num.replicate=100)
s2 - ScoreJackStraw(object=s2,dims = 1:8)
JackStrawPlot(s2, dims = 1:8)

ElbowPlot(s2)


###neighbor embedding#####
s2 = RunTSNE(s2, dims=1:8)
DimPlot(s2, reduction = "tsne")

s2 = FindNeighbors(s2, dims=1:8)
s2 = FindClusters(s2, resolution=.2)
s2 = RunUMAP(s2, dims=1:7)
DimPlot(s2, reduction = "umap", pt.size = .35, label = T)

#s2 = FindNeighbors(s2, dims=1:9)
#s2 = FindClusters(s2, resolution=.2)
#s2 = RunUMAP(s2, dims=1:8)
#DimPlot(s2, reduction = "umap", pt.size = .35, label = F)

current.cluster.ids = c(0,1,2,3,4,5,6,7,8,9)
new.cluster.ids =c("Fib","HSL","ML","B","LP","End","Tcells","Mac/DC","pericytes","TPCs")
names(new.cluster.ids) = levels(s2)     
s2 = RenameIdents(s2, new.cluster.ids)
DimPlot(s2, reduction="umap", label=T, pt.size=.05, cols = c("gray82","#8A9FD1","#89C75F","#90D5E4","#E6C2DC","gray82","gray82","gray82","gray82","#F37B7D"))

#HSL = "#8A9FD1", 
#Tum = "orange", 
# LP = "#E6C2DC", 
# B  = "#90D5E4",
# ML = "#89C75F", 
#TPCs= "#F37B7D"
#everything else = "gray82


s2.sub = subset(s2, idents = c("B","LP","TPCs"))

s2.subset.markers = FindAllMarkers(s2.sub,only.pos = T, min.pct= 0.25, thresh.use = 0.25)
s2.sub.top20= s2.subset.markers%>% group_by(cluster) %>% top_n(25, avg_log2FC)
DoHeatmap(object= s2.sub, features = s2.sub.top20$gene) + scale_fill_gradientn(colors = c("royal blue", "snow2", "orangered"))


