######S4 vehicle########
s4.data = Read10X(data.dir = "/Volumes/JOSHUA/scRNA-seq/raw.data/s4.pbs_out/filtered_feature_bc_matrix/")
str(s4.data)

s4 = new("seurat", raw.data = s4.data)
str(s4)
class(s4)
s4 = CreateSeuratObject(counts = s4.data, min.cells= 5, min.features = 200, project = "s4_treated_het")

###mito qc
s4[["percent.mt"]]= PercentageFeatureSet(s4, pattern = "^mt-")
s4[["percent.mt"]]

VlnPlot(s4, features = c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3, cols = "skyblue3")

s4.plot1 = FeatureScatter(s4, feature1="nCount_RNA", feature2 = "percent.mt",smooth = TRUE,span = TRUE,col="skyblue3")
s4.plot1
s4.plot2 = FeatureScatter(s4, feature1="nCount_RNA",feature2= "nFeature_RNA",smooth=TRUE,span=TRUE, col="skyblue3")
s4.plot2      
CombinePlots(list(s4.plot1,s4.plot2))  


s4 = subset(s4, subset = nFeature_RNA>200 & nFeature_RNA<3750 & percent.mt <5)
s4 = NormalizeData(s4, normalization.method = "LogNormalize", scale.factor = 1e4)



#####variable gene expression
s4 = FindVariableFeatures(s4, selection.method = "vst", nfeatures = 2000)
s4.top10 = head(VariableFeatures(s4),10)
s4.top10
s4.top20 = head(VariableFeatures(s4),20)
s4.top20
s4.top50 = head(VariableFeatures(s4),50)
s4.top50
s4.plot3 = VariableFeaturePlot(s4)

s4.top10.p = LabelPoints(plot=s4.plot3, points = s4.top10, repel=TRUE)
s4.top20.p = LabelPoints(plot=s4.plot3, points = s4.top20, repel=TRUE)
s4.top50.p = LabelPoints(plot=s4.plot3, points = s4.top50, repel=TRUE)

plot(s4.top20.p)
CombinePlots(list(s4.top20.p,s1.top20.p))


####Scale and Principal component analysis#
s4.all.genes = row.names(s4)
s4 = ScaleData(s4, features=s4.all.genes)
s4 = RunPCA(s4, features = VariableFeatures(s4), genes.print=10)

#####visualization of data
VizDimLoadings(s4,dims=1:7, reduction="pca",col="rosybrown")
DimPlot(s4, reduction="pca")

s4 = JackStraw(s4, num.replicate=100)
s4 - ScoreJackStraw(object=s4,dims = 1:8)
JackStrawPlot(s4, dims = 1:8)

ElbowPlot(s4)


###neighbor embedding
s4 = RunTSNE(s4, dims=1)
DimPlot(s4, reduction = "tsne")

s4 = FindNeighbors(s4, dims=1:5)
s4 = FindClusters(s4, resolution=.2)
s4 = RunUMAP(s4, dims=1:5)
DimPlot(s4, reduction = "umap", pt.size = 1, label = TRUE, label.size = 6)

s4.current.cluster.ids = c(0,1,2,3,4,5,6,7,8,9)
s4.new.cluster.ids =c("B","ML","LP","HSL","E","F","Mac/DC","E","ML","E")
names(s4.new.cluster.ids)
levels(s4)
names(s4.new.cluster.ids) = levels(s4)     
s4 = RenameIdents(s4, s4.new.cluster.ids)
DimPlot(s4, reduction="umap", label=TRUE, pt.size=1.3)

p1 = DimPlot(s4, reduction = "umap", pt.size = .1, label = F, cols =c("#90D5E4","#89C75F","#E6C2DC","#8A9FD1","gray76","gray76","gray76","gray76","gray76","#89C75F","gray76"))
#HSL = "#8A9FD1", 
#Tum = "orange", 
# LP = "#E6C2DC", 
# B  = "#90D5E4",
# ML = "#89C75F", 
#TPCs= "#F37B7D"

  FeaturePlot(s4, features = c("Krt8","Krt14"), cols = c( "gray82","goldenrod1","indianred"))
  FeaturePlot(s4, features = c("Top2a","Wfdc18"), cols = c( "gray82","goldenrod1","indianred"))
  FeaturePlot(s4, features = c("Prlr","Ptprc"), cols = c( "gray82","goldenrod1","indianred"))
  FeaturePlot(s4, features = c("Pecam1","Fn1"), cols = c( "gray82","goldenrod1","indianred"))
  
