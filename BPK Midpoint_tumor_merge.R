#s2.fvb50 subset merge
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(monocle3)
library(RColorBrewer)
library(SingleCellExperiment)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(metap)
library(tibble)
library(AnnotationHub)
library(ensembldb)
library(DESeq2)
library(metap)
library(multtest)
library("devtools")
library(infercnv)
library(stringr)
library(pheatmap)
library(RColorBrewer)
library(scran)
library(scDblFinder)
library(harmony)
library(SeuratWrappers)



fvb50.data = Read10X(data.dir = "/Users/joshuarivera/OneDrive - University of Massachusetts Boston/Shared_SPJR/scRNA-seq/raw data/tumor_2018/filtered_feature_bc_matrix/")
s2.data = Read10X(data.dir = "/Users/joshuarivera/OneDrive - University of Massachusetts Boston/Shared_SPJR/scRNA-seq/raw data/s2_dfci/filtered_feature_bc_matrix/")


DimPlot(fvb50, reduction = "umap")
DimPlot(s2, reduction = "umap", label =T)
fvb50.subset = subset(fvb50, idents = c("tumor", "TPCs","LP"))
s2.subset = subset(s2, idents = c("HSL","B","ML","LP","TPCs"))

  p1 = DimPlot(fvb50.subset)
  p2 = DimPlot(s2.subset)
p1+p2  

s2.fvb50.subset = merge(s2.subset, y=fvb50.sub, add.cell.ids= c("midpoint","tumor"), project  = "TPCs")
  s2.fvb50.subset
s2.fvb50.subset[["percent.mt"]] = PercentageFeatureSet(s2.fvb50.subset, pattern = "^mt-")
  s2.fvb50.subset[["percent.mt"]]
  
VlnPlot(s2.fvb50.subset, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)

plot1 = FeatureScatter(s2.fvb50.subset, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 = FeatureScatter(s2.fvb50.subset, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1+plot2

s2.fvb50.subset = subset(s2.fvb50.subset, subset = nFeature_RNA>200 & nFeature_RNA<10000 & percent.mt <5)
s2.fvb50.subset = NormalizeData(s2.fvb50.subset, normalization.method = "LogNormalize", scale.factor = 1e4)

s2.fvb50.subset = FindVariableFeatures(s2.fvb50.subset, selection.method = "vst", nfeatures = 2000)
  s2.fvb50.subset.10 = head(VariableFeatures(s2.fvb50.subset),10)
  s2.fvb50.subset.plot = VariableFeaturePlot(s2.fvb50.subset)
  s2.fvb50.subset.top10.p = LabelPoints(plot = s2.fvb50.subset.plot, points = s2.fvb50.subset.10, repel = T, ynudge = 0, xnudge = 0)  

all.genes = rownames(s2.fvb50.subset)
s2.fvb50.subset = ScaleData(s2.fvb50.subset, features = all.genes)

s2.fvb50.subset = RunPCA(object = s2.fvb50.subset, features = VariableFeatures(object=s2.fvb50.subset), genes.print = 10)
ElbowPlot(s2.fvb50.subset)

s2.fvb50.subset = FindNeighbors(s2.fvb50.subset, dims = 1:15)
s2.fvb50.subset = FindClusters(s2.fvb50.subset, resolution = .4)
s2.fvb50.subset = RunUMAP(s2.fvb50.subset, dims = 1:10) # use 7
  p3 = DimPlot(s2.fvb50.subset, reduction = "umap", group.by = "orig.ident", pt.size = .15)
  p4 = DimPlot(s2.fvb50.subset, reduction = "umap", label=T)
  
p3+p4
#FeaturePlot(s2.fvb50.subset, features = c("Krt14","Krt8","Top2a","Csn3","Mfge8"))

DimPlot(s2.fvb50.subset, reduction = "umap", label = T)

###assigning cluster identities
s2.fvb50.subset.current.cluster.ids = c(0,1,2,3,4,5,6,7,8,9,10 )
s2.fvb50.subset.new.cluster.ids =c("HSL","LP","B","ML","tumor","tumor","RS-T","RS-LP","ML","dblt","HSL")
names(s2.fvb50.subset.new.cluster.ids) = levels(s2.fvb50.subset)     
s2.fvb50.subset = RenameIdents(s2.fvb50.subset, s2.fvb50.subset.new.cluster.ids)
DimPlot(s2.fvb50.subset, reduction="umap", label=F, pt.size=.15, label.size = 4, cols = c(HSL,LP,B,ML,Tum,RS_T,RS_LP,"gray82"))  

HSL = "#8A9FD1"
Tum = "orange"
LP = "#89C75F" 
B  = "#90D5E4"
ML = "#E6C2DC" 
RS_LP= "#F37B7D"
RS_T = "navyblue"
dblt = "#8A9FD1"
#everything else = "gray82
