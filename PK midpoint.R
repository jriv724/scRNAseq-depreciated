###sample 3 Brca1wt/wt, Trp53flx/flx, K14cre postive injected for 13 injections #######
s3.data = Read10X(data.dir = "/Volumes/JOSHUA/scRNA-seq/raw.data/s3.treated_out/outs/filtered_feature_bc_matrix/")
str(s3.data)

s3 = new("seurat", raw.data = s3.data)
str(s3)
class(s3)
s3 = CreateSeuratObject(counts = s3.data, min.cells= 5, min.features = 200, project = "s3_treated_het")

###mito qc
s3[["percent.mt"]]= PercentageFeatureSet(s3, pattern = "^mt-")
s3[["percent.mt"]]

VlnPlot(s3, features = c("nFeature_RNA","nCount_RNA","percent.mt"))

s3.plot1 = FeatureScatter(s3, feature1="nCount_RNA", feature2 = "percent.mt",smooth = TRUE,span = TRUE,col="slategray3")
s3.plot1
s3.plot2 = FeatureScatter(s3, feature1="nCount_RNA",feature2= "nFeature_RNA",smooth=TRUE,span=TRUE, col="slategray3")
s3.plot2      
CombinePlots(list(s3.plot1,s3.plot2))  


#####variable gene expression
s3 = FindVariableFeatures(s3, selection.method = "vst", nfeatures = 2000)
s3.top10 = head(VariableFeatures(s3),10)
s3.top10
s3.top20 = head(VariableFeatures(s3),20)
s3.top20
s3.top50 = head(VariableFeatures(s3),50)
s3.top50
s3.plot3 = VariableFeaturePlot(s3)

s3.top10.p = LabelPoints(plot=s3.plot3, points = s3.top10, repel=TRUE)
s3.top20.p = LabelPoints(plot=s3.plot3, points = s3.top20, repel=TRUE)
s3.top50.p = LabelPoints(plot=s3.plot3, points = s3.top50, repel=TRUE)

plot(s3.top20.p)
CombinePlots(list(s3.top20.p,top20.p))


####Scale and Principal component analysis
s3.all.genes = row.names(s3)
s3 = ScaleData(s3, features=s3.all.genes)
s3 = RunPCA(s3, features = VariableFeatures(s3), genes.print=10)

#####visualization of data
VizDimLoadings(s3,dims=1:7, reduction="pca",col="darkslateblue")
DimPlot(s3, reduction="pca", cols = "darkslateblue")

s3 = JackStraw(s3, num.replicate=100)
s3 - ScoreJackStraw(object=s3,dims = 1:8)
JackStrawPlot(s3, dims = 1:8)

ElbowPlot(s3)


###neighbor embedding
s3 = RunTSNE(s3, dims=1:6)
DimPlot(s3, reduction = "tsne")


#actual
s3 = FindNeighbors(s3, dims=1:5)
s3 = FindClusters(s3, resolution=.3)
s3 = RunUMAP(s3, dims=1:5)

Idents(s3) = s3$RNA_snn_res.0.3
DimPlot(s3, reduction = "umap", pt.size = .35, label = T)



#####Assigning cell type identity to clusters
s3.current.cluster.ids = c(0,1,2,3,4,5,6,7,8,9)
s3.new.cluster.ids =c("B","End","HSL","ML","AV-LP","AV-LP","Fib","P","Im","HSL")
names(s3.new.cluster.ids) = levels(s3)     
s3 = RenameIdents(s3, s3.new.cluster.ids)
DimPlot(s3, reduction="umap", label=TRUE, pt.size=.75, label.size = .5)
p1 = DimPlot(s3, reduction = "umap", pt.size = .1, label = F, cols =c("#89C75F","#90D5E4","gray76","#8A9FD1","#E6C2DC","gray76","gray77"))


  FeaturePlot(s3, features = c("Krt8","Krt14"), cols = c( "gray82","goldenrod1","indianred"))
  FeaturePlot(s3, features = c("Top2a","Wfdc18"), cols = c( "gray82","goldenrod1","indianred"))
  FeaturePlot(s3, features = c("Prlr","Ptprc"), cols = c( "gray82","goldenrod1","indianred"))
  FeaturePlot(s3, features = c("Pecam1","Fn1"), cols = c( "gray82","goldenrod1","indianred"))

saveRDS(s3, file = "OneDrive - University of Massachusetts Boston/Shared_SPJR/scRNA-seq/RDS_obj/s3.rds")
  #HSL = "#8A9FD1", 
  #Tum = "orange", 
  # LP = "#E6C2DC", 
  # B  = "#90D5E4",
    # ML = "#89C75F", 
  #TPCs= "#F37B7D"
  #evthing else = "gray76"
