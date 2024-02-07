library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(monocle3)
library(RColorBrewer)
library(tibble)
####S2 - Brca1 wt.flx, Trp53 flx/flx, Krt14cre 12x 4NQO1####
s2.data = Read10X(data.dir = "/Volumes/JOSHUA/scRNA-seq/raw.data/s2.treated_out/filtered_feature_bc_matrix/")
str(s2.data)

s2 = new("seurat", raw.data = s2.data)
str(s2)
class(s2)
s2 = CreateSeuratObject(counts = s2.data, min.cells= 5, min.features = 200, project = "midpoint_treated")


saveRDS(s2, file= "/Users/joshuarivera/s2_treated.rds")
readRDS(s2, file= "/Users/joshuarivera/s2_treated.rds")

###mito qc####
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



#####variable gene expression#####
s2 = FindVariableFeatures(s2, selection.method = "vst", nfeatures = 2000)
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


####Scale and Principal component analysis#####
s2.all.genes = row.names(s2)
s2 = ScaleData(s2, features=s2.all.genes)
s2 = RunPCA(s2, features = VariableFeatures(s2), genes.print=10)

#####visualization of data####
VizDimLoadings(s2,dims=1:7, reduction="pca",col="rosybrown")
DimPlot(s2, reduction="pca", cols = "rosybrown")

s2 = JackStraw(s2, num.replicate=100)
s2 - ScoreJackStraw(object=s2,dims = 1:7)
JackStrawPlot(s2, dims = 1:8)

ElbowPlot(s2)


###neighbor embedding#####
s2 = RunTSNE(s2, dims=1:8)
DimPlot(s2, reduction = "tsne")

s2 = FindNeighbors(s2, dims=1:8)
s2 = FindClusters(s2, resolution=.2)
s2 = RunUMAP(s2, dims=1:8)
DimPlot(s2, reduction = "umap", pt.size = 1, label = F)


#alt for dot plot
s2 = FindNeighbors(s2, dims=1:7)
s2 = FindClusters(s2, resolution=.4)
s2 = RunUMAP(s2, dims=1:8)
DimPlot(s2, reduction = "umap", pt.size = 1, label = T)


FeaturePlot(s2, features= c("Wfdc18","Krt14","Krt5","Krt8","Mki67","Prlr","Cdh1","Csn3","Lgals7"),label = TRUE, label.size = .2, cols = c("grey","royal blue"))
FeaturePlot(s2, features= c("Wfdc18","Krt5","Aldh1a3","Krt18","Prlr","Cdh1","Csn3","Ccna2"),label = TRUE)
FeaturePlot(s2, features=c("Krt79","Itga2","Elf5","Kit","Mfge8","Vim","Lyz2"), pt.size=1,label = TRUE)
FeaturePlot(s2, features=c("Aldh1a3","Kit","Cd14","Lyz2"),label = TRUE)

FeaturePlot(s2, features=c("Il19","Krt79"), pt.size=1)

FeaturePlot(s2, features= c("Wfdc18","Csn3","Ltf","Mfge8","Muc1","Lalba","Trf"))

#epithelial stem cell gene
FeaturePlot(s2, features=c("Aldh1a1","Aldh1a3"))

###stem cell tf for luminal breast#
FeaturePlot(s2, features =c("Gata3","Ly6a","Cd14","Prlr","Csn2","Gng11","Lalba","Krt14"))


#updated epithelial markers
myl = c("Krt17","Krt14","Krt5","Acta2","Myl9","Mylk","Myh11")
lum = c("Krt19","Krt18","Krt8")
hs = c("Prlr","Cited1","Pgr","Prom1","Esr1")
av = c("Mfge8","Trf","Csn3","Wfdc18","Elf5","Ltf")
lp = c ("Kit","Aldh1a3","Cd14","Wap")
milk = c("Glycam1","Olah")

fibro = c("Col1a1","Col1a2","Col3a1","Fn1")
panVL = c("Pecam1","Cdh5","Eng")
ve = c("Sox17","Sele")
p = c("Rgs5","Des","Notch3")
LE = c("Mmrn1","Prox1","Flt4","Ccl21a")
Im = c("Ptprc")
My = c("Cd74","Lyz2")
DC = c("Traf1","Cd209a","Napsa","Flt3")
Mac = c("Csf1r","Fcgr3","Adgre1","Ms4a7")
Ma  = c("Mrc1","Cd209f","Cd163")
Mb  = c("Mmp12","Mmp13","Spic")
Ly  = c("Cd3g","Cd3d","Cd3e")
NK  = c("Gzma","Ncr1","Itgae")
Tcd8= c("Cd8a","Cd8b1")
Tcd4= c("Cd4")
B   = c("Cd79a","Cd79b")




FeaturePlot(s2, features = myl, label = TRUE)
FeaturePlot(s2, features = lum, label = TRUE)
FeaturePlot(s2, features = hs, label = TRUE)
FeaturePlot(s2, features = av, label = TRUE)
FeaturePlot(s2, features = lp, label = TRUE)
FeaturePlot(s2, features = milk, label = TRUE)

FeaturePlot(s2, features = fibro, label = TRUE)
FeaturePlot(s2, features = panVL, label = TRUE)
FeaturePlot(s2, features = ve, label = TRUE)
FeaturePlot(s2, features = p, label = TRUE)
FeaturePlot(s2, features = LE, label = TRUE)
FeaturePlot(s2, features = Im, label = TRUE)
FeaturePlot(s2, features = My, label = TRUE)
FeaturePlot(s2, features = DC, label = TRUE)
FeaturePlot(s2, features = Mac, label = TRUE)
FeaturePlot(s2, features = Ma, label = TRUE)
FeaturePlot(s2, features = Mb, label = TRUE)
FeaturePlot(s2, features = Ly, label = TRUE)
FeaturePlot(s2, features = NK, label = TRUE)
FeaturePlot(s2, features = Tcd8, label = TRUE)
FeaturePlot(s2, features = Tcd4, label = TRUE)
FeaturePlot(s2, features = B, label = TRUE)


#####pub clusters#####
#0 - Vascular endothelial/pericyte/Lymphatic endothelial
#1 - Luminal HS 
#2 - Myoepithelial 
#3 - Luminal
#4 - Luminal HS-AV
#5 - Luminal AV
#6 - T/B cells 
#7 - Macrophage (Ma)
#8 - fibroblasts
#9 - TPCs

#####DotPlot classifiers#####
#0 - Luminal HS
#1 - Myoepithleial
#2 - Lumial AV
#3 - Luminal HS-AV
#4 - Vascular endothelial
#5 - Lympatic enodothelial
#6 - Fibroblasts
#7 - Cd4+/Cd8+ T cells
#8 - Macrophage / DCs
#9 - TPCs
#10- Pericytes
#11- Luminal HS dblt
#12- Luminal HS


####break####
#mk1 of cluster names####
s2.current.cluster.ids = c(0,1,2,3,4,5,6,7,8,9)
#.tmp s2.new.cluster.ids =c("stroma","hormonal sensing","basal","alveolar/progenitor","mature luminal","endothelial","immune","immune","endothelial", "TPCs")
#s2.sub.ids = c("stroma","luminal","basal","luminal","luminal","endothelial","immune","immune","endothelial", "TPCs")
s2.sub.prog.ids = c("stroma","luminal","B","LP","luminal","endothelial","immune","immune","endothelial", "TPCs")
names(s2.sub.prog.ids) = levels(s2)     
s2 = RenameIdents(s2, s2.sub.prog.ids)
DimPlot(s2, reduction="umap", label=T, pt.size=1)
####updated####

#####Assigning cell type identity to clusters####### ##Super cluster UMAP IDs
current.cluster.ids = c(0,1,2,3,4,5,6,7,8,9)
new.cluster.ids =c("End","HSL","B","LP","ML","Fibro","Tcells","Mac/DC","pericytes","TPCs")
names(new.cluster.ids) = levels(s2)     
s2 = RenameIdents(s2, new.cluster.ids)
DimPlot(s2, reduction="umap", label=T, pt.size=.25)

#####Assigning cell type identity to clusters####### ##DotPlot classifiers 
current.cluster.ids = c(0,1,2,3,4,5,6,7,8,9,10,11,12)
new.cluster.ids =c("Luminal HS","Myoepithelial","Luminal AV","Luminal HS-AV","Vascular endothelia","Lymphatic endothelia","Fibroblasts","Cd4+/Cd8+ Tcells","Macrophage/ DCs", "TPCs", "Pericytes","Luminal HS dblt","Luminal HS")

names(new.cluster.ids) = levels(s2)     
s2 = RenameIdents(s2, new.cluster.ids)
DimPlot(s2, reduction="umap", label=F, pt.size=.75)


#######for ozgun analysis#########
current.cluster.ids = c(0,1,2,3,4,5,6,7,8,9)
new.cluster.ids =c("End","epi","epi","epi","epi","Fibro","Tcells","Mac/DC","pericytes","epi")
names(new.cluster.ids) = levels(s2)     
s2 = RenameIdents(s2, new.cluster.ids)
DimPlot(s2, reduction="umap", label=T, pt.size=.25)

s2.epi = FindMarkers(s2, ident.1 = "epi", min.pct = 0.25)
write.table(s2.epi, file = "/Users/joshuarivera/Desktop/cluster_list/4nqo_epithelial.xls", sep="\t", quote = FALSE,
            col.names=NA)


#s2.dotplot without subsetting and differential clustering
x = subset(s2, idents = c("HSL", "B","LP","ML","TPCs"))
s2.epitheliamk2 =  DotPlot(s2,idents = c("HSL", "B","LP","ML","TPCs"), cols = c("yellow","slateblue"), features = c(myl,lum,hs,av,lp))
s2.stromamk2 = DotPlot(s2, idents = c("TPCs","End","Fibro","Tcells", "Mac/DC","pericytes"),  cols = c("yellow","blue"), features = c(fibro, panVL, ve, p, LE, Im, My, DC, Mac, Ly, NK, Tcd8, Tcd4, B), cluster.idents = TRUE)


s2.epithelia = subset(s2, idents = c("HSL", "B","LP","ML","TPCs"))
  x = DotPlot(s2.epithelia, cols = c("light grey","firebrick"), features = c(myl,lum,hs,av,lp))


mylevels =   
s2.epithelia.1 = subset(s2, idents = c("TPCs", "B","LP","ML","HSL"))
 x =  DotPlot(s2.epithelia.1, cols = c("light grey","firebrick"), features = c(myl,lum,hs,av,lp), cluster.idents = TRUE)
  

  
s2.stroma = subset(s2, idents =c("TPCs","Vascular endothelia","Lymphatic endothelia","Fibroblasts","Cd4+/Cd8+ Tcells", "Macrophage/ DCs","Pericytes"))
 x = DotPlot(s2.stroma, cols = c("light grey","firebrick"), features = c(fibro, panVL, ve, p, LE, Im, My, DC, Mac, Ly, NK, Tcd8, Tcd4, B), cluster.idents = TRUE)



FeaturePlot(s2, features = hsd, label = TRUE)
FeaturePlot(s2, features = hsp, label = TRUE)
FeaturePlot(s2, features = lp,  label = TRUE)
FeaturePlot(s2, features = avp, label = TRUE)
FeaturePlot(s2, features = avd, label = TRUE)
FeaturePlot(s2, features = bsl, label = TRUE)
FeaturePlot(s2, features = c(myo,lp), label = TRUE)
FeaturePlot(s2, features = c(prlf,"Aldh1a3","Brca1"), label = TRUE)
FeaturePlot(s2, features = c(emt))
FeaturePlot(s2, features = c("Lyz2","Ptprc"), label = TRUE)
FeaturePlot(s2, features = c("Vim","Cd32"), label = TRUE)


####cluster ID #######

#####extra feature plots of random interest#####
#genes of interest from tumor - project1#
FeaturePlot(s2, features = c("Il19","Krt79"));

#basal.epithelial
FeaturePlot(fvb50, features =c("Krt14", "Krt17","Wnt10a","Krt5"))

#tumor markers
FeaturePlot(s2, features=c("Cd44","Ccna2","Brca1","Trp53","Brca2","Rad51"))

#initiation gene found in pioneer colon cancer cells#
FeaturePlot(s2, features = "Nrip2")

## nk cell activation found in cancer###
FeaturePlot(fvb50, features = "Cd14")

####stem cell MaSC markers####
FeaturePlot(s2, features=c("Aldh1a3","Ccna2","Ccnb2"))
FeaturePlot(s3, features=c("Aldh1a3","Ccna2","Ccnd2"))
FeaturePlot(fvb50, features=c("Aldh1a3","Ccna2","Ccnd2"))



#########carmans markers#######
#luminal.progeniter
FeaturePlot(s2, features = c("Elf5","Cd14","Kit","Aldh1a1",
                             "Aldh1a3", "Foxc1","Foxc2","Itga2",
                             "Itgb3","Tspan8","Cd55","Prlr"))

#mature.luminal.1 markers
FeaturePlot(s2, features = c("Krt7","Krt8","Krt18","Krt19",
                             "Epcam","Prlr","Areg","Sox9",
                             "Prom1","Ly6a","Csn3","Mfge8"));


#mature.luminal.2 markers
FeaturePlot(s2, features = c("Muc1","Trf","Cd24a","Gata3",
                             "Alcam","Ar","Cd200","Cdh1",
                             "Ceacam1","Cited1","Cxcl15"));


#mature.luminal.3 markers
FeaturePlot(s2, features = c("Egfr","Egr1","Ehf","Elf1",
                             "Elf2","Elf3","Elf4","Erbb2",
                             "Erbb3","Esr1","Foxa1"))

#mature.luminal.4 markers                             
FeaturePlot(s2, features = c("Foxa2","Glycam1","Hes2","Hey1",
                             "Jag1","Jag2","Krt4","Lalba",  
                             "Notch1","Notch2","Notch3","Notch4"))

#mature.luminal.5 markers
FeaturePlot(s2, features = c("Nrg4","Numb","Pgr","s200a6",
                             "Sox10","Sox11","Sox12","Sox13",
                             "Sox15","Sox17","Sox18"))

#mature.luminal.6 markers
FeaturePlot(s2, features = c("Sox4","Sox5","Sox6","Sox7",
                             "Sox8","Stc1","Stc2","Wap",
                             "Wfdc18","Wnt4","Wnt5a","Ltf"))

#final.hormone.sensing.luminal
FeaturePlot(s2, features = c("Prlr","Pgr","Esr1","Ly6a",
                             "Prom1","Cited1","Esr1"))

##cd52 marks tcells
FeaturePlot(s2, features = c("Il6","Lyz2"))


#alveolar.luminal markers
FeaturePlot(s2, features = c("Csn3","Ltf","Mfge8","Muc1",
                             "Lalba","Trf","Wfdc18"))

#myeoepithelial.1
FeaturePlot(s2, features = c("Moxd1","Ncam1","Nes","Ngfr",
                             "Nrg1","Nt5e","Oxtr","Pard3b",
                             "Pdpn","Procr","Pygo2","Sfrp1" ))

#myeoepithelial.2 
FeaturePlot(s2, features = c("Snai1","Snai2","Snai3","Taz",
                             "Tgfbr1","Tgfbr2","Tgfbr3","Trp53",   
                             "Trp63","Twist1","Twist2","Wif1"))
#myeoepithelial.3 
FeaturePlot(s2, features = c("Wnt10a","Wnt11","Wnt6","Zeb1",
                             "Zeb2","Krt5","Krt14","Krt17",
                             "Cd44","Vim","Acta2","Myh11"))

#myeoepithelial.4
FeaturePlot(s2, features = c("Myl9","Mylk","Cdh2","Cdh3",
                             "Mmp3","Mmp9","Igfbp4","Fn1",
                             "Col1a1","Col1a2","Col3a1","Sparc"))

#myeoepithelial.5 
FeaturePlot(s2, features = c("Fzd7","Runx2","Cav1","Msn",
                             "Gng11","Dkk3","Bmp7","Lgals2",
                             "Thy1","Inpp5d","Trp73" ))

#myeoepithelial.6 
FeaturePlot(s2, features = c("Serpine1","Mme","Axin2","Itga6",
                             "Itgav","Lgr5","Lgr6","Lrp5",
                             "Lrp6","Icam1","Cd47","Id4" ))

#adipocyte
FeaturePlot(s2, features = c("Pparg","Adipor1","Adipor2","Fabp4",
                             "Slc27a1","Slc27a4","Slc2a4","Ppargc1a",
                             "Slc27a6"))
#fibroblast 
FeaturePlot(s2, features = c("Pdgfra","Pdgfrb","Acta2",
                             "Fap","s200a4","Vim","Fn1",
                             "Pparg","Cd47","Cd34"))
#macrophage  
FeaturePlot(s2, features = c("Cd14","Cd68","Cd163","Csf1r",
                             "Adgre1","Itgam","Vim","Fn1",
                             "Aif1","Spi1","Lyz2","Ptprc","Inos"))

#monocyte
FeaturePlot(s2, features = c("Cd14","Fcgr1","Itgam","Csf1r",
                             "Lyz2","Vmo1","Ccl2","s200a9"))

#basophil 
FeaturePlot(s2, features = c("Il3ra","Gstm7","Scin",
                             "Lmo4","Stx3","Ifitm1",
                             "Ltb4r1","Cd63","Hdc","Ptprc" ))

#tcell 
FeaturePlot(s2, features = c("Cd2","Cd3d","Cd3e","Cd3g",
                             "Cd244","Cd8b1","Ccl5","Nkg7",
                             "Cd4","Il7r","Ptprc"))

#bcell 
FeaturePlot(s2, features = c("Cd79a", "Cd79b", "Blnk",
                             "Sfn","Ptprc","Cd274","Cd80",
                             "Cd86","Cd70"))

#endothelial 
FeaturePlot(s2, features = c("Pecam1", "Vwf","Cdh5","Sele",
                             "Kdr","Rasip1","Rgs5","Cav1",
                             "Msn","Cd34","Eng","s2pr1","Emcn"))

FeaturePlot(s2, features =c("Nrip2","Klf4"))

VlnPlot(s2, features = c("Lgr5","Krt5","Aldh1a3"))

#Cd274 - Pd-L1
#Cd80 - Ctla4 ligand

##### Cluster marker lists #####
s2.basal = FindMarkers(s2, ident.1 = "B", min.pct = 0.25)
print(s2.basal)
write.table(s2.basal, file = "/Users/joshuarivera/OneDrive - University of Massachusetts Boston/Shared_SPJR/Breast cancer paper/Supplimental/Cluster DGE in BPK midpoint/Basal.xls", sep="\t", quote = FALSE,
            col.names=NA)
s2.end = FindMarkers(s2, ident.1 = "End", min.pct = 0.25)
write.table(s2.end, file = "/Users/joshuarivera/OneDrive - University of Massachusetts Boston/Shared_SPJR/Breast cancer paper/Supplimental/Cluster DGE in BPK midpoint/Endothelial.xls", sep="\t", quote = FALSE,
            col.names=NA)
s2.fib = FindMarkers(s2, ident.1 = "Fib", min.pct = 0.25)
write.table(s2.fib, file = "/Users/joshuarivera/OneDrive - University of Massachusetts Boston/Shared_SPJR/Breast cancer paper/Supplimental/Cluster DGE in BPK midpoint/Fibroblasts.xls", sep="\t", quote = FALSE,
            col.names=NA)
s2.hsl = FindMarkers(s2, ident.1 = "HSL", min.pct = 0.25)
write.table(s2.hsl, file = "/Users/joshuarivera/OneDrive - University of Massachusetts Boston/Shared_SPJR/Breast cancer paper/Supplimental/Cluster DGE in BPK midpoint/HSL.xls", sep="\t", quote = FALSE,
            col.names=NA)
s2.ml = FindMarkers(s2, ident.1 = "ML", min.pct = 0.25)
write.table(s2.ml, file = "/Users/joshuarivera/OneDrive - University of Massachusetts Boston/Shared_SPJR/Breast cancer paper/Supplimental/Cluster DGE in BPK midpoint/ML.xls", sep="\t", quote = FALSE,
            col.names=NA)
s2.lp = FindMarkers(s2, ident.1 = "LP", min.pct = 0.25)
write.table(s2.lp, file = "/Users/joshuarivera/OneDrive - University of Massachusetts Boston/Shared_SPJR/Breast cancer paper/Supplimental/Cluster DGE in BPK midpoint/LP.xls", sep="\t", quote = FALSE,
            col.names=NA)
s2.TPCs = FindMarkers(s2, ident.1 = "TPCs", min.pct = 0.25)
write.table(s2.TPCs, file = "/Users/joshuarivera/OneDrive - University of Massachusetts Boston/Shared_SPJR/Breast cancer paper/Supplimental/Cluster DGE in BPK midpoint/RSC.xls", sep="\t", quote = FALSE,
            col.names=NA)
s2.tcells = FindMarkers(s2, ident.1 = "Tcells", min.pct = 0.25)
write.table(s2.tcells, file = "/Users/joshuarivera/OneDrive - University of Massachusetts Boston/Shared_SPJR/Breast cancer paper/Supplimental/Cluster DGE in BPK midpoint/lymphocytes.xls", sep="\t", quote = FALSE,
            col.names=NA)
s2.mac = FindMarkers(s2, ident.1 = "Mac/DC", min.pct = 0.25)
write.table(s2.mac, file = "/Users/joshuarivera/OneDrive - University of Massachusetts Boston/Shared_SPJR/Breast cancer paper/Supplimental/Cluster DGE in BPK midpoint/myeloid.xls", sep="\t", quote = FALSE,
            col.names=NA)
s2.peri = FindMarkers(s2, ident.1 = "pericytes", min.pct = 0.25)
write.table(s2.peri, file = "/Users/joshuarivera/OneDrive - University of Massachusetts Boston/Shared_SPJR/Breast cancer paper/Supplimental/Cluster DGE in BPK midpoint/pericytes.xls", sep="\t", quote = FALSE,
            col.names=NA)


s4.basal = FindMarkers(s4, ident.1 = "B", min.pct = 0.25)
write.table(s4.basal, file = "/Users/joshuarivera/OneDrive - University of Massachusetts Boston/Shared_SPJR/Breast cancer paper/Supplimental/DGE sheets/BPK vehicle/Basal.xls", sep="\t", quote = FALSE,
            col.names=NA)
s4.end = FindMarkers(s4, ident.1 = "E", min.pct = 0.25)
write.table(s4.end, file = "/Users/joshuarivera/OneDrive - University of Massachusetts Boston/Shared_SPJR/Breast cancer paper/Supplimental/DGE sheets/BPK vehicle/Endothelial.xls", sep="\t", quote = FALSE,
            col.names=NA)
s4.fib = FindMarkers(s4, ident.1 = "F", min.pct = 0.25)
write.table(s4.fib, file = "/Users/joshuarivera/OneDrive - University of Massachusetts Boston/Shared_SPJR/Breast cancer paper/Supplimental/DGE sheets/BPK vehicle/Fibroblasts.xls", sep="\t", quote = FALSE,
            col.names=NA)
s4.hsl = FindMarkers(s4, ident.1 = "HSL", min.pct = 0.25)
write.table(s4.hsl, file = "/Users/joshuarivera/OneDrive - University of Massachusetts Boston/Shared_SPJR/Breast cancer paper/Supplimental/DGE sheets/BPK vehicle/HSL.xls", sep="\t", quote = FALSE,
            col.names=NA)
s4.ml = FindMarkers(s4, ident.1 = "ML", min.pct = 0.25)
write.table(s4.ml, file = "/Users/joshuarivera/OneDrive - University of Massachusetts Boston/Shared_SPJR/Breast cancer paper/Supplimental/DGE sheets/BPK vehicle/LP.xls", sep="\t", quote = FALSE,
            col.names=NA)
s4.lp = FindMarkers(s4, ident.1 = "LP", min.pct = 0.25)
write.table(s4.lp, file = "/Users/joshuarivera/OneDrive - University of Massachusetts Boston/Shared_SPJR/Breast cancer paper/Supplimental/DGE sheets/BPK vehicle/ML.xls", sep="\t", quote = FALSE,
            col.names=NA)
s4.mac = FindMarkers(s4, ident.1 = "Mac/DC", min.pct = 0.25)
write.table(s4.mac, file = "/Users/joshuarivera/OneDrive - University of Massachusetts Boston/Shared_SPJR/Breast cancer paper/Supplimental/DGE sheets/BPK vehicle/myeloid.xls", sep="\t", quote = FALSE,
            col.names=NA)


s3.basal = FindMarkers(s3, ident.1 = "B", min.pct = 0.25)
write.table(s3.basal, file = "/Users/joshuarivera/OneDrive - University of Massachusetts Boston/Shared_SPJR/Breast cancer paper/Supplimental/DGE sheets/PK midpoint/Basal.xls", sep="\t", quote = FALSE,
            col.names=NA)
s3.end = FindMarkers(s3, ident.1 = "End", min.pct = 0.25)
write.table(s3.end, file = "/Users/joshuarivera/OneDrive - University of Massachusetts Boston/Shared_SPJR/Breast cancer paper/Supplimental/DGE sheets/PK midpoint/Endothelial.xls", sep="\t", quote = FALSE,
            col.names=NA)
s3.fib = FindMarkers(s3, ident.1 = "Str", min.pct = 0.25)
write.table(s3.fib, file = "/Users/joshuarivera/OneDrive - University of Massachusetts Boston/Shared_SPJR/Breast cancer paper/Supplimental/DGE sheets/PK midpoint/Fibroblasts.xls", sep="\t", quote = FALSE,
            col.names=NA)
s3.hsl = FindMarkers(s3, ident.1 = "HSL", min.pct = 0.25)
write.table(s3.hsl, file = "/Users/joshuarivera/OneDrive - University of Massachusetts Boston/Shared_SPJR/Breast cancer paper/Supplimental/DGE sheets/PK midpoint/HSL.xls", sep="\t", quote = FALSE,
            col.names=NA)
s3.ml = FindMarkers(s3, ident.1 = "ML", min.pct = 0.25)
write.table(s3.ml, file = "/Users/joshuarivera/OneDrive - University of Massachusetts Boston/Shared_SPJR/Breast cancer paper/Supplimental/DGE sheets/PK midpoint/ML.xls", sep="\t", quote = FALSE,
            col.names=NA)
s3.lp = FindMarkers(s3, ident.1 = "LP", min.pct = 0.25)
write.table(s3.lp, file = "/Users/joshuarivera/OneDrive - University of Massachusetts Boston/Shared_SPJR/Breast cancer paper/Supplimental/DGE sheets/PK midpoint/LP.xls", sep="\t", quote = FALSE,
            col.names=NA)
s3.mac = FindMarkers(s3, ident.1 = "Im", min.pct = 0.25)
write.table(s3.mac, file = "/Users/joshuarivera/OneDrive - University of Massachusetts Boston/Shared_SPJR/Breast cancer paper/Supplimental/DGE sheets/PK midpoint/myeloid.xls", sep="\t", quote = FALSE,
            col.names=NA)



for(i in 0:10)
{
  cluster= paste('s2', 'cluster', i, 'markers')
  cluster = FindMarkers(s2, ident.1 = i, min.pct=0.25)
  filename = paste("/Volumes/JOSHUA/scRNA-seq/sample.analysis/s2_het_treat/s2_markers/cluster",i,".xls")
  write.table(cluster, file=filename, sep="\t", quote=FALSE, col.names=NA)
}

s2@active.ident = s2$seurat_clusters

#####Assigning cell type identity to clusters#######
s2.current.cluster.ids = c(0,1,2,3,4,5,6,7,8,9,10)
s2.new.cluster.ids = c("End","HSL","ML","B","LP","Fib","Mac_DC","tcells","P","RSC","HSL")
names(s2.new.cluster.ids) = levels(s2)     
s2 = RenameIdents(s2, s2.new.cluster.ids)
DimPlot(s2, reduction="umap", label=T, pt.size=1)


FeaturePlot(s2, features=c("Ccna2","Top2a","Mki67","Ccnd2"))
VlnPlot(s2, features = c("Cdca8","Cdca3","Cdc25c","Tacc3") )
######interesting markers from this dataset 
#cd14- found in pioneer bladder cancer cells https://www.ncbi.nlm.nih.gov/pubmed/25825750
#mfg-e8 - enhances cell proliferation cyclin D1 - correlation to stem cells? 
#tgfb is important 
#Nrip2 in cluster 5 is shown to be involved in early cancers for colon cancer
#Galectin-1 https://en.wikipedia.org/wiki/Galectin-1 
#Cd34 cancer stem cells in living cancer https://www.ncbi.nlm.nih.gov/pubmed/25519836

#hypothesis does Brca1 happloinsufficiency promote transition to stem cell phenotype?#
s2.dense =(x = as.matrix(x= s2.data))  
s2.sparse=(x=s2.data)
print(s2.sparse)
str(a)

s2.sub = subset(s2, idents = c("TPCs","B","LP"))
DimPlot(s2.sub, reduction = "umap")
####heat map#######
#all luminal vs Basal vs TPCs
s2.sub.markers = FindAllMarkers(s2.sub,only.pos = T, min.pct= 0.25, thresh.use = 0.25)
s2.sub.top30= s2.sub.markers %>% group_by(cluster) %>% top_n(22, avg_log2FC)
DoHeatmap(object= s2.sub, features = s2.sub.top30$gene) + scale_fill_gradientn(colors = c("royalblue", "snow2", "orangered"))
#luminal progenitor vs Basal vs TPCs
s2.sub.prog = subset(s2, idents = c("TPCs","B","LP"))
DimPlot(s2.sub.prog, reduction = "umap")
s2.sub.prog.markers = FindAllMarkers(s2.sub.prog,only.pos = T, min.pct= 0.25, thresh.use = 0.25)
?FindAllMarkers
s2.sub.prog.top30= s2.sub.prog.markers %>% group_by(cluster) %>% top_n(30, avg_log2FC)
DoHeatmap(object= s2.sub.prog, features = s2.sub.prog.top30$gene, draw.lines = F, size = 8)+ scale_fill_gradientn(colors = c("royal blue", "snow2", "orangered"))+ theme(text = element_text(size = 15))

#15 markers
s2.sub.prog.top15= s2.sub.prog.markers %>% group_by(cluster) %>% top_n(15, avg_log2FC)
DoHeatmap(object= s2.sub.prog, features = s2.sub.prog.top15$gene, draw.lines = F, size = 8)+ scale_fill_gradientn(colors = c("royal blue", "snow2", "orangered"))+ theme(text = element_text(size = 25))


#20 markers
s2.sub.prog.top20= s2.sub.prog.markers %>% group_by(cluster) %>% top_n(20, avg_log2FC)
DoHeatmap(object= s2.sub.prog, features = s2.sub.prog.top20$gene, draw.lines = F, size = 8)+ scale_fill_gradientn(colors = c("royal blue", "snow2", "orangered"))+ theme(text = element_text(size = 25))

s2=readRDS(file = "OneDrive - University of Massachusetts Boston/Shared_SPJR/scRNA-seq/RDS_obj/s2.rds")
DoHeatmap(s2, features = s2.sub.prog.top20$gene, draw.lines = F, size = 8)+ scale_fill_gradientn(colors = c("royalblue", "snow2", "orangered"))+ theme(text = element_text(size = 25))

s2markers= FindAllMarkers(object = s2, min.pct = 0.25)
DoHeatmap(object=s2, group.by = s2markers, draw.lines = F, size = 8)+ scale_fill_gradientn(colors = c("royalblue", "snow2", "orangered"))+ theme(text = element_text(size = 25))

s2.top15 <- s2markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)
DoHeatmap(s2, features = s2.top15$gene, draw.lines = F, size = 10)+ 
  theme(text = element_text(size = 20))+
  scale_fill_gradientn(colors = c("royalblue", "snow2", "orangered"))+ theme(text = element_text(size = 25))





