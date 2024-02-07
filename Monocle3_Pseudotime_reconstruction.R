########s2.fvb50.epi merge pseudotime#######
s2.fvb50.cds = as.cell_data_set(s2.fvb50.subset)
s2.fvb50.cds@rowRanges@elementMetadata@listData[["gene_short_name"]] = rownames(s2.fvb50.subset[["RNA"]])
s2.fvb50.cds = preprocess_cds(s2.fvb50.cds, num_dim = 100)
s2.fvb50.cds = cluster_cells(s2.fvb50.cds)
s2.fvb50.cds = learn_graph(s2.fvb50.cds)
s2.fvb50.cds = order_cells(s2.fvb50.cds, reduction_method = "UMAP")
s2.fvb50.cds_pr_test_res = graph_test(s2.fvb50.cds, neighbor_graph="principal_graph", cores=16)
pr_deg_ids = row.names(subset(s2.fvb50.cds_pr_test_res, q_value < 0.05))
gene_module_df = find_gene_modules(s2.fvb50.cds[pr_deg_ids,], resolution=1e-2)

plot_cells(
  cds = s2.fvb50.cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE
)

plot_cells(
  s2.fvb50.cds,
  label_groups_by_cluster = T,
  label_leaves = T,
  label_branch_points = T
)


p1 <- plot_cells(s2.fvb50.cds, show_trajectory_graph = TRUE)
p2 <- plot_cells(s2.fvb50.cds, color_cells_by = "partition", show_trajectory_graph = TRUE)
p1+p2  

plot_cells(
  s2.fvb50.cds,
  label_groups_by_cluster = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE
)

###identifying critical genes for pseudotime differentiation
#subset morans score genes for graphing
s2.fvb50.gt = graph_test(s2.fvb50.cds, neighbor_graph = "knn", cores = 8)
s2.fvb50.ids = row.names(subset(s2.fvb50.gt, q_value < 0.05))
  y = order(s2.fvb50.gt[s2.fvb50.ids,]$morans_I %in% !NA, decreasing = TRUE)
    z = s2.fvb50.gt[y,]
      a = order(z$morans_I, decreasing = T)
fin = s2.fvb50.gt[a,]
      #View(s2.fvb50.gt[a,])
  
    write.table(s2.fvb50.gt[s2.fvb50.ids,], file = "/Users/joshuarivera/Desktop/weighted_matrix_s2fvb50_sig.xls", sep = "\t", quote = FALSE, col.names = NA)

DotPlot(s2.fvb50.subset, features = c("Krt6a", rownames(fin[1:20,])), cols = c( "gold2","indianred"), group.by ="orig.ident")+ coord_flip()

#Finding modules of co-regulated genes
s2fvb50_module.df = gene_module_df
#The data frame gene_module_df contains a row for each gene and identifies the module it belongs to. 
#To see which modules are expressed in which clusters or partitions you can use two different approaches for visualization. T
#The first is just to make a simple table that shows the aggregate expression of all genes in each module across all the clusters. 
#Monocle provides a simple utility function called aggregate_gene_expression for this purpose:
cell_group.df = tibble::tibble(cell = row.names(colData(s2.fvb50.cds)),
                               cell_group = s2.fvb50.subset@active.ident)
agg_mat = aggregate_gene_expression(s2.fvb50.cds, s2fvb50_module.df, cell_group.df)
  row.names(agg_mat) = stringr::str_c("Module",row.names(agg_mat))
  colnames(agg_mat) = stringr::str_c("Partition", colnames(agg_mat))
pheatmap(agg_mat, cluster_rows = T, cluster_cols = T,
         scale="column", clustering_method = "ward.D2", fontsize = 6)



####plotting genes across pseudotime
plot_cells(s2.fvb50.cds, genes ,
           show_trajectory_graph = TRUE,
           label_cell_groups = F,
           label_leaves = T)
                                          
s2fvb50_module.df = find_gene_modules(s2.fvb50.cds[s2.fvb50.ids,], resolution=c(10^seq(-6,-1)))

#definind cluster internal "partitions"
cell_group.df = tibble::tibble(cell = row.names(colData(s2.fvb50.cds)),
                               cell_group = s2.fvb50.cds@colData$ident)

agg_mat = aggregate_gene_expression(s2.fvb50.cds, s2fvb50_module.df, cell_group.df)
  row.names(agg_mat) = stringr::str_c("Module",row.names(agg_mat))
 
pheatmap(agg_mat,
         scale="column", clustering_method = "ward.D2", fontsize = 10)

plot_cells(s2.fvb50.cds,
           genes = s2fvb50_module.df %>% dplyr::filter(module %in% c(4,1,7)),
           label_cell_groups = TRUE,
           show_trajectory_graph = TRUE)

plot_cells(s2.fvb50.cds,
           genes = s2fvb50_module.df %>% dplyr::filter(module %in% c(9,10,13)),
           label_cell_groups = TRUE,
           show_trajectory_graph = TRUE)

plot_cells(s2.fvb50.cds,
           genes = s2fvb50_module.df %>% dplyr::filter(module %in% c(8,12)),
           label_cell_groups = TRUE,
           show_trajectory_graph = TRUE)

mod1 = s2fvb50_module.df$module %in% 1
  module1_genes = s2fvb50_module.df[mod1,]
     write.table(module1_genes, file = "/Users/joshuarivera/OneDrive - University of Massachusetts Boston/Shared_SPJR/Pseudotime reconstruction/critical_gene_sets/module1_genes.xls", sep = "\t",
                 quote = F, col.names = NA)
mod4 = s2fvb50_module.df$module %in% 4
  module4_genes = s2fvb50_module.df[mod4,]
      write.table(module4_genes, file = "/Users/joshuarivera/OneDrive - University of Massachusetts Boston/Shared_SPJR/Pseudotime reconstruction/critical_gene_sets/module4_genes.xls", sep = "\t",
                 quote = F, col.names = NA)
mod7 = s2fvb50_module.df$module %in% 7
  module7_genes = s2fvb50_module.df[mod7,]
      write.table(module7_genes, file = "/Users/joshuarivera/OneDrive - University of Massachusetts Boston/Shared_SPJR/Pseudotime reconstruction/critical_gene_sets/module7_genes.xls", sep = "\t",
                 quote = F, col.names = NA)
mod9 = s2fvb50_module.df$module %in% 9
  module9_genes = s2fvb50_module.df[mod9,]
     write.table(module9_genes, file = "/Users/joshuarivera/OneDrive - University of Massachusetts Boston/Shared_SPJR/Pseudotime reconstruction/critical_gene_sets/module9_genes.xls", sep = "\t",
                 quote = F, col.names = NA)
mod2 = s2fvb50_module.df$module %in% 2
  module2_genes = s2fvb50_module.df[mod2,]
     write.table(module2_genes, file = "/Users/joshuarivera/OneDrive - University of Massachusetts Boston/Shared_SPJR/Pseudotime reconstruction/critical_gene_sets/module2_genes.xls", sep = "\t",
                 quote = F, col.names = NA)
mod16 = s2fvb50_module.df$module %in% 16
  module16_genes = s2fvb50_module.df[mod16,]
     write.table(module16_genes, file = "/Users/joshuarivera/OneDrive - University of Massachusetts Boston/Shared_SPJR/Pseudotime reconstruction/critical_gene_sets/module16_genes.xls", sep = "\t",
                 quote = F, col.names = NA)


#module8 is proliferation
#module 9 is HSL
  

plot_cells(s2.fvb50.cds,
           genes = s2fvb50_module.df %>% dplyr::filter(module %in% c(1,2,3,4)),
           label_cell_groups = TRUE,
           show_trajectory_graph = TRUE)

plot_cells(s2.fvb50.cds,
           genes = s2fvb50_module.df %>% dplyr::filter(module %in% c(5,6,7,8)),
           label_cell_groups = TRUE,
           show_trajectory_graph = TRUE)

plot_cells(s2.fvb50.cds,
           genes = s2fvb50_module.df %>% dplyr::filter(module %in% c(9,10,11,12)),
           label_cell_groups = TRUE,
           show_trajectory_graph = TRUE)

plot_cells(s2.fvb50.cds,
           genes = s2fvb50_module.df %>% dplyr::filter(module %in% c(13,14,15,16)),
           label_cell_groups = TRUE,
           show_trajectory_graph = TRUE)
