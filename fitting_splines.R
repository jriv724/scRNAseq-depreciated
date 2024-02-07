library(tidyverse)
library(openxlsx)
library(doParallel)
library(npreg)


xx = lp2t.cds@assays@data$logcounts #counts log transformed normalized counts from monocle3 
yy = lp2t.cds@principal_graph_aux$UMAP$pseudotime #pseudotime number scale / cell association
#xx <- readRDS('../../data/input/fromJosh/lognormCounts.rds')
#yy <- readRDS('../../data/input/fromJosh/pseudotimeData.rds')

yy <- data.frame(cells = names(yy), pt = c(yy))
cell.idents = as.data.frame(s2.fvb50.lp2t@active.ident)
cell.idents = rownames_to_column(cell.idents)
names(cell.idents)
colnames(cell.idents) = c("cells", "cellID")

cell.identsyy = left_join(yy, cell.idents, by = 'cells')

xx <- as.data.frame(as.matrix(xx))
xx$GeneName <- rownames(xx)
xx.long <- xx %>% pivot_longer(-GeneName, names_to = 'cells', values_to = 'expr')s2


xx.long <- left_join(xx.long, cell.identsyy, by = 'cells')
xx.long$st <- (xx.long$pt - min(xx.long$pt))/ (max(xx.long$pt) - min(xx.long$pt))
xx.long$cellType <- factor(xx.long$cellID, levels = c("LP", "RS-LP", "RS-T", "tumor"))
cell.to.time <- xx.long %>% dplyr::select(st, cellType) %>% distinct()
plot(cell.to.time$cellType, cell.to.time$st)

## Identify transition points
trans.time <- cell.to.time %>% group_by(cellType) %>% summarise(qq = quantile(st, probs = 0.75))
n.time <- length(trans.time$qq) 
trans.time$strt <- c(0,(trans.time$qq[1:(n.time - 1)] + trans.time$qq[2:(n.time)]) / 2)
trans.time$stp <- c((trans.time$qq[1:(n.time - 1)] + trans.time$qq[2:(n.time)]) / 2, 1)

saveRDS(trans.time, '~/scMouseTumor/data/input/rds/trans_time.rds')
## For parallel calculations

num.cores <- 16


genes <- unique(xx.long$GeneName)
i = which(genes == 'S100a1')
i = which(genes == 'Krt5')

lbx <- 2.0

## Expression graphs

spline.fits <- mclapply(1:length(genes), function(i){
  tmp <- xx.long %>% dplyr::filter(GeneName == genes[i]) %>%
    transmute(GeneID = GeneName, x = st, y = expr)
  
  
  y <- tmp$y
  t <- tmp$x
  
  rna.sp <- smooth.spline(t, y, lambda = 0.1)
  
  # w <- rep(1, length(y))
  # w[which(y == 0)] <- 1/3
  # 
  # rna.sp <- smooth.spline(t, y, cv = T, w = w)
  # sparx <- rna.sp$spar
  # 
  # if(sparx < lbx){
  #   
  #   rna.sp <- smooth.spline(t, y, spar = lbx, w = w)
  #   
  # }
  
  rna.sp <- predict(rna.sp, seq(0, 1, by = 0.01))
  
  ######
  #plot(t, y, col = 'red')
  #points(rna.sp$x, rna.sp$y, type = 'l', col = 'green', lwd = 2)
  mu <- data.frame(x = rna.sp$x, y = rna.sp$y)
  
  mu <- data.frame(GeneID = rep(tmp$GeneID[1], length(mu[,1])), x = mu[,1], y = mu[,2])
  return(mu)
}, mc.cores = num.cores)

spline.fits <- bind_rows(spline.fits)
saveRDS(spline.fits, '../../data/input/rds/spline_fits.rds')
