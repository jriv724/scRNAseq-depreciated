library(tidyverse)
library(openxlsx)
library(Seurat)
library(ggplot2)
library(ggfortify)
library(jcolors)
require(gridExtra)
library(grid)
library(tidytext)
library(Seurat)




## For parallel calculations
num.cores <- detectCores(all.tests = FALSE, logical = TRUE)

getCurvePeakLoc <- function(t, y, prob = 0.8){
  
  ## Fitting the estimated kernel with smooth splines
  spline.fit <- smooth.spline(x = t, y = y)
  
  
  ## Compute the derivatives of the fitted splines
  s.0 <- predict(spline.fit, spline.fit$x, deriv=0)
  s.1 <- predict(spline.fit, spline.fit$x, deriv=1)
  s.derv <- data.frame(s0=s.0$y, s1=s.1$y)
  
  ## Get the location of the extrema
  locs <- rle(den.sign <- sign(s.derv$s1))
  
  
  ## Maxima
  inc.ind <- which(locs$values == 1)
  if(length(inc.ind) > 1){
    maxima.ind = {}
    for(i in inc.ind){
      maxima.ind = c(maxima.ind, sum(locs$lengths[1:i]))
    }
    ## Interpolate a point between the location where derivative changes sign
    maxima = (spline.fit$x[maxima.ind] + spline.fit$x[(maxima.ind + 1)]) / 2
    maxima = maxima[!is.na(maxima)]
    ## Get the maximum values
    maxval = predict(spline.fit, maxima)
    
    ## Get the outliers
    maxima.outliers = which(maxval$y >= quantile(maxval$y, prob = prob))
    
    ## Peaks for entities of interest
    entity.x = maxval$x[maxima.outliers]
    entity.y = maxval$y[maxima.outliers]
  }else{
    entity.x <- spline.fit$x[which.max(spline.fit$y)]
    entity.y <- spline.fit$y[which.max(spline.fit$y)]
  }
  
  return(entity.x)
}


## Plot specific splines
spline.fits <- readRDS('~/scMouseTumor/data/input/rds/spline_fits.rds')
trans.time <- readRDS( '~/scMouseTumor/data/input/rds/trans_time.rds')

## DEGs
#RS_LP.vs.LP = FindMarkers(s2.fvb50.sub, ident.1 = "RS-LP", ident.2 = "LP", min.pct = .15)
  #saveRDS(RS_LP.vs.LP, file = "~/scMouseTumor/data/input/rds/RS_LP.vs.LP_DEGs.rds")
#Tumor.vs.RS_LP = FindMarkers(s2.fvb50.sub, ident.1 = "tumor", ident.2 = "RS-LP", min.pct = .15)
  #saveRDS(Tumor.vs.RS_LP, file = "~/scMouseTumor/data/input/rds/Tumor.vs.RS_LP_DEGs.rds")
#RS_LP.vs.B = FindMarkers(s2.fvb50.sub, ident.1 = "RS-LP", ident.2 = "B", min.pct = .15)
  #saveRDS(RS_LP.vs.B, file = "~/scMouseTumor/data/input/rds/RS_LP.vs.B_DEGs.rds")


RS_LP.vs.LP <- readRDS('~/scMouseTumor/data/input/rds/RS_LP.vs.LP_DEGs.rds')
Tumor.vs.RS_LP <- readRDS('~/scMouseTumor/data/input/rds/Tumor.vs.RS_LP_DEGs.rds')
RS_LP.vs.B <- readRDS('~/scMouseTumor/data/input/rds/RS_LP.vs.B_DEGs.rds')


## Combining all
RS_LP.vs.LP$contrast <- 'RS-LP.vs.LP'
Tumor.vs.RS_LP$contrast <- 'Tumor.vs.RS-LP'
RS_LP.vs.B$contrast <- 'RS-LP.vs.B'

marker.genes <- bind_rows(RS_LP.vs.LP, Tumor.vs.RS_LP, RS_LP.vs.B)
marker.genes$GeneID = rownames(marker.genes)
colnames(marker.genes) = c("GeneID","p_val","avg_log2FC","pct.1","pct.2","p_val_adj","contrast")


## Turn the data into wide format (time by gene) and center & scale each gene
sc.rna.dtw.wide <- spline.fits %>% 
  pivot_wider(names_from = 'GeneID', values_from = 'y') %>% 
  mutate_at(vars(-1), scale) %>%
  as.data.frame()


sc.rna.mu.scale <- sc.rna.dtw.wide %>% 
  pivot_longer(-x, names_to = 'GeneID', values_to = 'expr')

sc.rna.mu.scale <- sc.rna.mu.scale[!is.na(sc.rna.mu.scale$expr), ]
sc.rna.peak.order <- sc.rna.mu.scale %>% group_by(GeneID) %>% summarise(peak.ord = getCurvePeakLoc(x, expr))
sc.rna.mu.scale <- left_join(sc.rna.mu.scale, sc.rna.peak.order, by = 'GeneID')

sc.rna.mu.scale$GeneID <- factor(sc.rna.mu.scale$GeneID, 
                                 levels = unique(sc.rna.mu.scale$GeneID[order(-sc.rna.mu.scale$peak.ord)]))

x = rownames(marker.genes)
sc.rna.mu.scale.sig <- sc.rna.mu.scale[sc.rna.mu.scale$GeneID %in% x, ]

p1 <- ggplot(sc.rna.mu.scale.sig, aes(x = x, y = GeneID, fill = expr)) + 
  geom_tile() + 
  #facet_grid(phase~., scales = "free",  space='free',
   #          labeller=label_wrap_gen(multi_line = TRUE))+
  ylab("Genes") + xlab("time/cells") + ggtitle('expression') + 
  #scale_fill_gradientn(colours = hm.palette(10)) +
  scale_fill_gradientn(colours = viridis::inferno(10)) +
  theme(panel.spacing = unit(0.1, "lines")) + 
  theme(
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y  = element_blank(),
    plot.title = element_text(size=14, face="bold"),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold"),
    legend.position = "none") 
    
theme(plot.margin=unit(c(2.5,2.5,2.5,2.5),"cm"))

plot(p1)


p2 = ggplot(sc.rna.mu.scale.sig, aes(x = x, y = GeneID, fill = expr)) + 
  geom_tile() + 
  #facet_grid(phase~., scales = "free",  space='free',
  #          labeller=label_wrap_gen(multi_line = TRUE))+
  ylab("Genes") + xlab("time/cells") + ggtitle('expression') + 
  #scale_fill_gradientn(colours = hm.palette(10)) +
  scale_fill_gradientn(colours = viridis::inferno(10)) +
  theme(panel.spacing = unit(0.1, "lines")) 

plot(p2)

ggsave(filename="~/scMouseTumor/heatmap.png",
       plot=p1,
       width = 4, height = 9,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)
