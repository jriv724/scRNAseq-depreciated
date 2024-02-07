library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)

###Data call on local host######
saveRDS(fvb50, file= "/Users/joshuarivera/FVB50_tumor.rds")
readRDS(fvb50, file= "/Users/joshuarivera/FVB50_tumor.rds")
UpdateSeuratObject(fvb50)
###tumor sample V2 commented, V3 active########
fvb50.data = Read10X(data.dir = "/Volumes/JOSHUA/scRNA-seq/raw.data/fvb50_out/filtered_gene_bc_matrices/mm10")
View(fvb50.data)

fvb50.dense = object.size(x = as.matrix(x= fvb50.data))  
fvb50.dense
fvb50.sparse = object.size(x = fvb50.data)  
fvb50.sparse
fvb50.dense / fvb50.sparse  

#dated# fvb51 = CreateSeuratObject(raw.data = fvb51.data, min.cells= 5, min.genes=200, project = "8.20.18_seq.fvb51")
fvb50 = new("seurat", raw.data = fvb50.data)
str(fvb50)
class(fvb50)
fvb50 = CreateSeuratObject(counts = fvb50.data, min.cells= 5, min.features = 200, project = "FVB50_tumor")

fvb50

str(fvb50)

#mitochondria QC##########
#fvb50.mito.genes = grep("^mt-*",rownames(fvb50@data), value= T)
fvb50[["percent.mt"]] = PercentageFeatureSet(fvb50, pattern = "^mt-")
fvb50[["percent.mt"]]

#fvb50.percent.mito = colSums(expm1(fvb50@data[fvb50.mito.genes,]))/colSums(expm1(fvb50@data))
# fvb50.percent.mito
#fvb50 = AddMetaData(fvb50, fvb50.percent.mito, "percent.mito.genes")
str(fvb50) 

VlnPlot(fvb50, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)

plot1 = FeatureScatter(fvb50, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot1
plot2 = FeatureScatter(fvb50, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
#par(mfrow=c(1,1))
#qc.vln = VlnPlot(fvb50, c("nGene", "nUMI", "percent.mito.genes"), nCol = 3, col = "indian red")
# qc.vln
#GenePlot is used to observe gene-gene relationship######
#gene.umi.cor =GenePlot(fvb50, "nUMI", "nGene",
#col = "indian red",
#pch=16,
#cex=.6,
#do.spline = TRUE,
#spline.span = .75)

#mito.qc = GenePlot(fvb50, "nUMI", "percent.mito.genes", col = "blue", pch=16, cex=.6)
fvb50 = subset(fvb50, subset = nFeature_RNA > 500 & nFeature_RNA < 8000 & percent.mt < 5)

#from here we decide which genes to omit (x<y<z) omit genes high and low cut off#
#Q: Where do we call cut off percent mito. There is a large subset of my genes at fall 
#fvb50 = FilterCells(object = fvb50, subset.names = c("nGene", "percent.mito.genes"),
#low.thresholds = c(500, -Inf),
#high.thresholds = c(7500,0.05))


fvb50 = NormalizeData(fvb50, normalization.method = "LogNormalize", scale.factor = 1e4)
#fvb50 = NormalizeData(fvb50)
####Seurat Dimensionality reduction####
#fvb50 = NormalizeData(object = fvb50, normalization.method = "LogNormalize", scale.factor = 1e4)

####Variable genes within clusters
fvb50 = FindVariableFeatures(fvb50, selection.method = "vst", nfeatures = 2000)
top10 = head(VariableFeatures(fvb50), 10)
top10
top20 = head(VariableFeatures(fvb50), 20)
top50 = head(VariableFeatures(fvb50), 50)
plot1 = VariableFeaturePlot(fvb50)

top10.p = LabelPoints(plot = plot1, points = top10, repel = TRUE)
top20.p = LabelPoints(plot = plot1, points = top20, repel = TRUE)
top50.p = LabelPoints(plot = plot1, points = top50, repel = TRUE, ynudge = 0, xnudge=0)

plot(top20.p)


#fvb50 = FindVariableGenes(object = fvb50, mean.function = ExpMean,
#dispersion.function = LogVMR,
#x.low.cutoff= 0.0125,
#x.high.cutoff= 3,
#y.cutoff=0.5)
#length = (x =(fvb50@var.genes))
#View(fvb50@var.genes)
#length

#scaling data, and exclusion based on UMI and Mitochondial gene outliers#
#fvb50 = ScaleData(object = fvb50, vars.to.regress =c("nUMI","percent.mito.genes"))
all.genes = row.names(fvb50)
fvb50 = ScaleData(fvb50, features=all.genes)
#dimentiality reduction after data normalization#
fvb50 = RunPCA(object = fvb50, features = VariableFeatures(object = fvb50), genes.print = 15)
print(fvb50[["pca"]], dims = 1:10, nfeatures = 10)
#Seurat PCA graphing methods######### whis is where the magic happens#
#PrintPCA(object = fvb50, pcs.print = 1:8, genes.print = 5, use.full = FALSE)
#VizPCA(object = fvb50, pcs.use = 1:4) #show genes that compose the most differentially expressed that defines clusters#
PCAPlot(object = fvb50, dim.1 = 1, dim.2=2, dim.3=3, do.hover = TRUE, data.hover = c("PC2"))
#PCHeatmap(object= fvb50, pc.use=1:4, cells.use = 500, do.balanced = TRUE, label.columns = TRUE)
#adding metadata of PCA developed
fvb50 = ProjectPCA(object= fvb50, do.print = TRUE) #Ask Aeidin about this###
VizDimLoadings(fvb50,dims=1:5, reduction = "pca")
DimPlot(fvb50, reduction = "pca")

DoHeatmap(object=fvb50,group.by = fvb50.markers)
DoHeatmap
###Jackstraw to determine statistically significant principle componenets#####
#fvb50 = JackStraw(object = fvb50, num.replicate= 100, display.progress = TRUE)
#JackStrawPlot(object= fvb50, PCs= 1:8)
fvb50 = JackStraw(fvb50, num.replicate = 100)
fvb50 = ScoreJackStraw(fvb50, dims=1:15)
JackStrawPlot(fvb50, dims = 1:15)

#another way to skin the cat - less computationally expensive#
PCAPlotelbowPlot(object= fvb50)
PCAPlot
ElbowPlot(fvb50)


##Find cluster: We will define tSNE cluster distribution based on the PC 
##significance from jackstraw and PCElbow. For this analysis we will go with 7 clusters#
#fvb50 = FindClusters(object= fvb50, reduction.type = "pca", dims.use=1:9,
#                    resolution = 0.5, print.output = 1, save.SNN= TRUE)
fvb50 = FindNeighbors(fvb50, dims=1:8)
fvb50 = FindClusters(fvb50, resolution = .3)

hist(fvb50$nFeature_RNA)
density(fvb50$nCount_RNA)
View(fvb50$nCount_RNA)
#define tSNE clusters print report
#PrintFindClustersParams(object=fvb50)
#cluster IDs of the first 5 cells
head(Idents(fvb50),25)
 
fvb50 = RunUMAP(fvb50, dims=1:8)
fvb50 = RunTSNE(fvb50, dims=1:8)

###non-linear dimensional reduction (GENERATION OF tSNE)######
fvb50 = RunTSNE(object= fvb50, dims.use = 1:8)
DimPlot(fvb50, reduction = "tsne")
DimPlot(fvb50, reduction = "umap", pt.size = 1,label = T)
DimPlot(fvb50, reduction = 'tsne')
#TSNEPlot(object = fvb50,
#         do.hover = TRUE,
#        data.hover = "PC1")
TSNEPlot(object = fvb50,
         pt.size = .75)

fvb50 = FindNeighbors(fvb50, dims=1:8)
fvb50 = FindClusters(fvb50, resolution = .3)
fvb50 = RunUMAP(fvb50, dims=1:5)
DimPlot(fvb50, reduction = "umap", pt.size = 1,label = T)

fvb50.sub = subset(fvb50, idents = c("tumor","TPCs","LP"))
###ALL ABOUT THE MARKERS ####
####marker lists###
{
fvb50.cluster0.markers = FindMarkers(object = fvb50, ident.1 = 0, min.pct = 0.25)
write.table(fvb50.cluster0.markers, file="/Users/joshuarivera/Desktop/V3.list/tumor2.xls", sep="\t", quote=FALSE,
            col.names= NA)


fvb50.cluster1.markers = FindMarkers(object = fvb50, ident.1 = 1, min.pct = 0.25)
print(cluster1.markers)
write.table(cluster1.markers, file="/Users/joshuarivera/Desktop/V3.list/tumor1.xls", sep="\t", quote=FALSE,
            col.names= NA)

fvb50.cluster2.markers = FindMarkers(object = fvb50, ident.1 = 2, min.pct=0.25)
print(cluster2.markers)
write.table(cluster2.markers, file="/Users/joshuarivera/Desktop/V3.list/myofibroblasts.xls", sep="\t", quote=FALSE,
            col.names= NA)

fvb50.cluster3.markers = FindMarkers(object = fvb50, ident.1 = 3, min.pct=0.25)
print(cluster3.markers)
write.table(cluster3.markers, file="/Users/joshuarivera/Desktop/V3.list/alveolar luminal.xls", sep="\t", quote=FALSE,
            col.names= NA)

fvb50.cluster4.markers = FindMarkers(object= fvb50, ident.1= 4, min.pct= 0.25)
print(cluster4.markers)
write.table(cluster4.markers, file="/Users/joshuarivera/Desktop/V3.list/mature luminal.xls", sep="\t", quote=FALSE,
            col.names= NA)

fvb50.cluster5.markers = FindMarkers(object= fvb50, ident.1= 5, min.pct= 0.25)
print(cluster5.markers)
write.table(cluster5.markers, file="/Users/joshuarivera/Desktop/V3.list/endothelial.xls", sep="\t", quote=FALSE,
            col.names= NA)

fvb50.cluster6.markers = FindMarkers(object= fvb50, ident.1= 6, min.pct= 0.25)
print(cluster6.markers)
write.table(cluster6.markers, file="/Users/joshuarivera/Desktop/V3.list/tumor macrophage.xls", sep="\t", quote=FALSE,
            col.names= NA)


fvb50.cluster8.markers = FindMarkers(fvb50, ident.1= 8, min.pct= 0.25)
print(cluster8.markers)
write.table(cluster7.markers, file="/Users/joshuarivera/Desktop/transitional population.xls", sep="\t", quote=FALSE,
            col.names= NA)
fvb50.cluster8.markers = FindMarkers(fvb50, ident.1= 8, min.pct = 0.25)}



#########
for(i in 0:7)
{
  cluster= paste('fvb50', 'cluster', i, 'markers')
  cluster = FindMarkers(fvb50, ident.1 = i, min.pct=0.25)
  filename = paste("/Users/joshuarivera/Desktop/fvb50_markers/fvb50.cluster",i,".xls")
  write.table(cluster, file=filename, sep="\t", quote=FALSE, col.names=NA)
}
####all markers#####
fvb50.markers = FindAllMarkers(object = fvb50, only.pos = TRUE, min.pct= 0.25, thresh.use = 0.25)
write.table(fvb50.markers, file="/Users/joshuarivera/Desktop/V3.list/FVB50.all.xls", sep="\t", quote=FALSE,
            col.names=NA)
fvb50.markers %>% group_by(cluster) %>% top_n(2,avg_logFC)

fvb50.markers
print(tails(fvb50.markers))

##roc curve for expectation maximization in each cluster classification##
cluster1.markers = FindMarkers(object= fvb50, ident.1 = 1, thresh.use = 0.25, test.use="roc", only.pos = FALSE)
print(cluster1.markers)
####Vln plot to show gene relationships based on clusters########
par(mfrow=c(1,2))
FeaturePlot(fvb50, features= c("Cd74"))
VlnPlot(fvb12, features = c("Il19","Krt79","Csf1","Csf1r","Gpx2"))
VlnPlot(fvb50, features = c("Brca1","Brca2","Trp53bp1", "Rad51"))
VlnPlot(object = fvb50, features.plot  =c("Zeb1","Itgb4","Snai3"))

CellPlot(fvb50,cell1 = 1, cell2 = 2, gene.ids = c("Krt14"),cols.use="black", nrpoints.use = Inf,
         pch=19, cex.lab=.5)
DotPlot(fvb50,genes.plot=c("Il19","Krt79","Csf1","Csf1r"),cols.use ="red",
        col.min= -2, col.max=4, dot.scale=6, scale.by="radius")
DotPlot()
DotPlot(fvb50,genes.plot=c("Krt8","Krt18","Krt14","Krt17"),cols.use ="red",
        col.min= -2, col.max=4, dot.scale=6, scale.by="radius")

RidgePlot(fvb50, features=c("Il19","Krt79","Csf1","Csf1r"))
##Vlnplot for for raw UMIs##
VlnPlot(fvb50, features =c("Krt8","Krt18","Krt14","Gpx2"))
VlnPlot(fvb50, features =c("Il19","Krt79","Csf1","Csf1r","Gpx2"), slot="counts", log=TRUE)
VlnPlot(fvb50, features= c("Brca1","Brca2","Trp53bp1", "Rad51"))

hsd = c("Cited1","Prlr","Esr1","Areg","Foxa1","Ly6a")
#hormone sensing progenitor
hsp = c("Aldh1a3","Kit","Cd14","Lypd3","Prlr")
#luminal progenitors
lp  = c("Aldh1a3","Kit","Cd14","Lypd3")
#alveolar progenitors
avp = c("Aldh1a3","Kit","Cd14","Lypd3","Rspo1")
#alveolar differentiated
avd = c("Fabp3","Thrsp","Wap","Glycam1","Olah")
#Basal epithelial cells
bsl = c("1500015O10Rik","Col7a1","Moxd1","Mia","Emid1","Igfbp3","Fst","Il17b","Emid1")
#myoepithelial cells
myo = c("Oxtr","Krt15","Igfbp6","Tnfvb50")
#proliferation markers
prlf= c("Mki67","Ccna2","Ccnd2","Rad51")
#general markers
gen.mkr = c("Aldh1a3","Krt14","Krt8","Wfdc18","Prlr","Csn3","Vim","Lyz2","Ccna2")

FeaturePlot(fvb50, features = hsd, label = TRUE)
FeaturePlot(fvb50, features = hsp, label = TRUE)
FeaturePlot(fvb50, features = lp,  label = TRUE)
FeaturePlot(fvb50, features = avp, label = TRUE)
FeaturePlot(fvb50, features = avd, label = TRUE)
FeaturePlot(fvb50, features = bsl, label = TRUE)
FeaturePlot(fvb50, features = c(myo,lp), label = TRUE)
FeaturePlot(fvb50, features = c(prlf), label = TRUE)
FeaturePlot(fvb50, features = c("Krt7","Krt14","Nrip2","Tagln"))

####PLOTS FOR GENE CLUSTER EXPRESSION####
luminal.progeniter = FeaturePlot(fvb50, features= c("Elf5","Cd14","Kit","Aldh1a1",
                                                    "Aldh1a3", "Foxc1","Foxc2","Itga2",
                                                    "Itgb3","Tspan8","Cd55"))
luminal.progeniter

project1 = FeaturePlot(fvb50, features = c("Il19","Krt79","Krt14"))
project1

basal.epithelial = FeaturePlot(fvb50, features =c("Krt14", "Krt17","Wnt10a","Krt5"))
basal.epithelial
tumor.markers = FeaturePlot(object=fvb50, features.plot=c("Cd44","Ccna2","Brca1","Trp53","Brca2","Rad51",
                                                          "Trp53bp1","Cdcp1","Cd80"),
                            cols.use=c('grey','red'),
                            reduction.use = "tsne",
                            pt.size = .5)


#cdcp1 - triple negative breast cancer prognosis marker  
par(mfrow=c(1,2))
#from alveolar 
#Mki67 - proliferation
#Cd68 - exhausted TAM
#Cd80 - binds Ctla4 
#

#####for paper#####
FeaturePlot(fvb50, features = c("Krt14","Krt8","Ccna2","Krt19", "Wfdc18","Prlr",
"Ccnd2","Sox9", "Prom1","Ly6a","Csn3","Mfge8"), cols = c('gray','red2'), label = T, pt.size =.5)

FeaturePlot(fvb50, feartures = "Nrip2")

mature.luminal.1 = FeaturePlot(fvb50, features = c("Krt7","Krt8","Krt18","Krt19",
                                                   "Epcam","Prlr","Areg","Sox9",
                                                   "Prom1","Ly6a","Csn3","Mfge8"));
mature.luminal.1

mature.luminal.2 = FeaturePlot(fvb50, features= c("Muc1","Trf","Cd24a","Gata3",
                                                  "Alcam","Ar","Cd200","Cdh1",
                                                  "Ceacam1","Cited1","Cxcl15"));
mature.luminal.2 

cd52 = FeaturePlot(fvb50, features=c("Cd52","Lyz2"))
  cd52

mature.luminal.3 = FeaturePlot(object=fvb50, features.plot = c("Egfr","Egr1","Ehf","Elf1",
                                                               "Elf2","Elf3","Elf4","Erbb2",
                                                               "Erbb3","Esr1","Foxa1"),
                               cols.use = c("grey","red"),
                               reduction.use="tsne")
mature.luminal.4 = FeaturePlot(object=fvb50, features.plot =c("Foxa2","Glycam1","Hefvb50","Hey1",
                                                              "Jag1","Jag2","Krt4","Lalba",  
                                                              "Notch1","Notch2","Notch3","Notch4"),
                               cols.use = c("grey","red"),
                               reduction.use="tsne")
mature.luminal.5 = FeaturePlot(object=fvb50, features.plot =c("Nrg4","Numb","Pgr","fvb5000a6",
                                                              "Sox10","Sox11","Sox12","Sox13",
                                                              "Sox15","Sox17","Sox18"),
                               cols.use = c("grey","red"),
                               reduction.use="tsne") 
mature.luminal.6 = FeaturePlot(object=fvb50, features.plot=c("Sox4","Sox5","Sox6","Sox7",
                                                             "Sox8","Stc1","Stc2","Wap",
                                                             "Wfdc18","Wnt4","Wnt5a","Ltf"),
                               cols.use=c('grey','red'),
                               reduction.use="tsne")

final.hormone.sensing.luminal = FeaturePlot(object=fvb50, features.plot= c("Prlr","Pgr","Esr1","Ly6a",
                                                                           "Prom1","Cited1","Esr1"),
                                            cols.use= c('grey','red'),
                                            reduction.use="tsne")

alveolar.luminal = FeaturePlot(fvb50, features=c("Csn3","Ltf","Mfge8","Muc1",
                                                             "Lalba","Trf","Wfdc18"))
 alveolar.luminal
 

myeoepithelial.1 = FeaturePlot(object=fvb50, features.plot=c("Moxd1","Ncam1","Nes","Ngfr",
                                                             "Nrg1","Nt5e","Oxtr","Pard3b",
                                                             "Pdpn","Procr","Pygo2","Sfrp1" ),
                               cols.use=c('grey','red'),
                               reduction.use="tsne")
myeoepithelial.2 = FeaturePlot(object=fvb50, features.plot=c("Snai1","Snai2","Snai3","Taz",
                                                             "Tgfbr1","Tgfbr2","Tgfbr3","Trp53",   
                                                             "Trp63","Twist1","Twist2","Wif1"),
                               cols.use=c('grey','red'),
                               reduction.use="tsne")

myeoepithelial.3 = FeaturePlot(object=fvb50, features.plot=c("Wnt10a","Wnt11","Wnt6","Zeb1",
                                                             "Zeb2","Krt5","Krt14","Krt17",
                                                             "Cd44","Vim","Acta2","Myh11"),
                               cols.use=c('grey','red'),
                               reduction.use="tsne")

myeoepithelial.4 = FeaturePlot(object=fvb50, features.plot=c("Myl9","Mylk","Cdh2","Cdh3",
                                                             "Mmp3","Mmp9","Igfbp4","Fn1",
                                                             "Col1a1","Col1a2","Col3a1","Sparc"),
                               cols.use=c('grey','red'),
                               reduction.use="tsne")

myeoepithelial.5 = FeaturePlot(object=fvb50, features.plot=c("Fzd7","Runx2","Cav1","Msn",
                                                             "Gng11","Dkk3","Bmp7","Lgalfvb50",
                                                             "Thy1","Inpp5d","Trp73" ),
                               cols.use=c('grey','red'),
                               reduction.use="tsne")

myeoepithelial.6 = FeaturePlot(object=fvb50, features.plot=c("Serpine1","Mme","Axin2","Itga6",
                                                             "Itgav","Lgr5","Lgr6","Lrp5",
                                                             "Lrp6","Icam1","Cd47","Id4" ),
                               cols.use=c('grey','red'),
                               reduction.use="tsne")

adipocyte = FeaturePlot(object= fvb50, features.plot = c("Pparg","Adipor1","Adipor2","Fabp4",
                                                         "Slc27a1","Slc27a4","Slc2a4","Ppargc1a"
                                                         ,"Slc27a6"),
                        cols.use=c('grey','red'),
                        reduction.use="tsne")
fibroblast = FeaturePlot(fvb50, features.plot=c("Pdgfra","Pdgfrb", "Acta2",
                                                "Fap","fvb5000a4", "Vim","Fn1",
                                                "Pparg",  "Cd47",   "Cd34"  ))
                         
 macrophage = FeaturePlot(fvb50, features=c("Cd14","Cd68","Cd163","Csf1r",
                                                                    "Adgre1","Itgam","Vim","Fn1",
                                                                    "Aif1","Spi1","Lyz2","Ptprc" ), label = T)
macrophage
                         
monocyte = FeaturePlot(object=fvb50, features.plot=c("Cd14","Fcgr1","Itgam","Csf1r", "Lyz2","Vmo1","Ccl2","fvb5000a9"), cols.use=c('grey','red'), reduction.use = "umap",
                       label = T)
                         
basophil = FeaturePlot(object=fvb50, features.plot=c("Il3ra","Gstm7","Scin",
                                                                              "Lmo4","Stx3","Ifitm1",
                                                                              "Ltb4r1","Cd63","Hdc","Ptprc" ),
                                                cols.use=c('grey','red'),
                                                reduction.use="tsne")
                         
tcell = FeaturePlot(object=fvb50, features.plot = c("Cd2","Cd3d","Cd3e","Cd3g",
                                                                             "Cd244","Cd8b1","Ccl5","Nkg7",
                                                                             "Cd4","Il7r","Ptprc"),
                                             cols.use = c('grey','red'),
                                             reduction.use="tsne")
bcell= FeaturePlot(object=fvb50, features.plot = c("Cd79a", "Cd79b", "Blnk",
                                                                            "Sfn","Ptprc","Cd274","Cd80",
                                                                            "Cd86","Cd70"),
                                            cols.use=c('grey','red'),
                                            reduction.use="tsne")
#Cd274 - Pd-L1
#Cd80 - Ctla4 ligand
endothelial = FeaturePlot(object=fvb50, features.plot=c("Pecam1", "Vwf","Cdh5","Sele",
                                                                                 "Kdr","Rasip1","Rgs5","Cav1",
                                                                                 "Msn","Cd34","Eng","fvb50pr1", 
                                                                                 "Emcn"  ),
                                                   cols.use = c('grey','red'),
                                                   reduction.use = "tsne",
                                                   pt.size=.5)

FeaturePlot(fvb50, features = c("Usp1","Krt14","Krt17","Cdh3","Nrg1","Ncam1"))
                         
#telomerase
telomerase = FeaturePlot(object=fvb50, features.plot=c("Tep1"),
                                                  cols.use=c('grey','red'),
                                                  reduction.use="tsne",
                                                  pt.size=.5)
######heatmap generation#######
fvb50.markers %>% group_by(cluster) %>% top_n(10, avg_logFC) -> top10
DoHeatmap(object = fvb50, genes.use = top10$gene, slim.col.label = TRUE, remove.key = FALSE,
                                   col.low = "royal blue",col.mid = "white",col.high = "red")
DoHeatmap
                         
                         
#####Assigning cell type identity to clusters#######
current.cluster.ids = c(0,1,2,3,4,5,6,7)
                        
new.cluster.ids =c("tumor","tumor","CAF","TBL-1","L","Im","endothelial","TBL-2")

names(new.cluster.ids) = levels(fvb50)     
fvb50 = RenameIdents(fvb50, new.cluster.ids)
DimPlot(fvb50, reduction="umap", label=F, pt.size=.75)

FeaturePlot(fvb50, features = c("Mki67","Ccna2","Ccnd2","Rad51"), label = TRUE)

####subsetting FVB50 to just epithelial cells#####
fvb50sub = subset(fvb50, idents = c("TPCs","basal","tumor","alveolar luminal","mature luminal"))
DimPlot(fvb50sub, reduction="umap", label=TRUE, pt.size=1, label.size = 4)

fvb50sub.markers = FindAllMarkers(fvb50sub, only.pos = T, min.pct= 0.25, thresh.use = 0.25)
fvb50sub.top50= fvb50sub.markers %>% group_by(cluster) %>% top_n(50, avg_logFC)
DoHeatmap(object= fvb50sub, features = fvb50sub.top30$gene)
