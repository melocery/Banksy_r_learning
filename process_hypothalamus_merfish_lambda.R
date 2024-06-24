# process_merfish_main.R
rm(list=ls())
graphics.off()
out.dir = 'out/'
check <- dir.exists(out.dir)
if (!check) dir.create(out.dir)

results.dir = 'out/lambda'
check <- dir.exists(results.dir)
if (!check) dir.create(results.dir)

USE_PROVIDED_BANKSY_OBJ = TRUE
number_of_animals = 2 # set this to 11 to use all naive animals. 

library(Banksy)
library(gridExtra) # grid.arrange
library(ggplot2) # facet_wrap 
# library(Seurat)
library(scales) # show_col
library(data.table) # fread
# library(irlba)
# library(Matrix)
# library(tidyverse)
library(plyr) # mapvalues
# library(scran)
library(tictoc)
library(ComplexHeatmap)
# library(dbscan)
# library(circlize)
library(peakRAM)
# library(qvalue) # check how many of these are needed.

results.dir = 'out/lambda/'
data.dir = 'data/'

k_geom = 15;lambda = c(0.2, 0.4, 0.6, 0.8, 1.0);npcs = 20;k_expr = 50;res = 0.5
list_of_animal_IDs = 1:number_of_animals 

# # Load data
all_mfish = fread('data/merfish_all_cells.csv') # see the readme file in the data dir
all_mfish <- all_mfish[,-c('Fos')]# remove Fos gene per Moffitt manuscript
all_mfish = cbind(cell_ids = paste0('cell_', 1:nrow(all_mfish)), all_mfish)
m.list = lapply(list_of_animal_IDs, function(x) all_mfish[all_mfish$Animal_ID==x,])
expr.list <- lapply(m.list, function(x){
  expr = t(as.matrix(x[,-c(1:10)]))
  colnames(expr) <- x$cell_ids
  expr})
locs.list <- lapply(m.list, function(x){
  df = as.data.frame(cbind(sdimx = x$Centroid_X, sdimy = x$Centroid_Y,
                           sdimz = 1000*x$Bregma))
  rownames(df) <- x$cell_ids
  df})
names(expr.list) <- names(locs.list) <- paste0('Animal_', list_of_animal_IDs)
  
# ------ start clustering block ---------
bank <- BanksyObject(own.expr = expr.list, cell.locs = locs.list)
all_mfish_by_animal_list = lapply(list_of_animal_IDs, function(x) all_mfish[Animal_ID %in% x,
                                                                            c('cell_ids',
                                                                              'Animal_ID',
                                                                              'Animal_sex',
                                                                              'Behavior',
                                                                              'Bregma',
                                                                              'Centroid_X',
                                                                              'Centroid_Y',
                                                                              'Cell_class',
                                                                              'Neuron_cluster_ID')])
bank@meta.data = cbind(bank@meta.data, do.call(rbind, all_mfish_by_animal_list))
  
# bank <- ComputeBanksy(bank, k_geom = k_geom)
ram_compute_banksy = peakRAM(bank <- ComputeBanksy(bank, k_geom = k_geom))

# bank <- ScaleBanksy(bank)
ram_scale_banksy = peakRAM(
  bank <- ScaleBanksy(bank)
)
  
lambdas<-c(0, lambda)
  
# bank <- Banksy:::RunBanksyPCA(bank, lambda = lambdas, npcs = npcs)
ram_pca_banksy = peakRAM(
  bank <- Banksy:::RunBanksyPCA(bank, lambda = lambdas, npcs = npcs)
)
  
# bank <- Banksy:::RunBanksyUMAP(bank, lambda = lambdas, npcs = npcs, nneighbors = k_expr)
ram_umap_banksy = peakRAM(
  bank <- Banksy:::RunBanksyUMAP(bank, lambda = lambdas, npcs = npcs, nneighbors = k_expr)
)
  
clust_major_numeric = as.numeric(factor(bank@meta.data$Cell_class))
clust_major_names = as.character(factor(bank@meta.data$Cell_class))
bank@meta.data = cbind(bank@meta.data, clust_major_num = clust_major_numeric)
set.seed(42)
ptm <- proc.time()
# bank <- ClusterBanksy(bank, lambda = lambdas, pca = TRUE, npcs = npcs,
#                       method = 'leiden', k.neighbors = k_expr, resolution = res)
ram_cluster_banksy = peakRAM(
  bank <- ClusterBanksy(bank, lambda = lambdas, pca = TRUE, npcs = npcs,
                        method = 'leiden', num.cores = 8,
                        k.neighbors = k_expr, resolution = res)
)
time_taken = proc.time() - ptm
print(time_taken)
saveRDS(bank, file = paste0(data.dir, 'banksyObj_naive_lambda_', number_of_animals, '.rds'))
# bank <- readRDS(file = paste0(data.dir, 'banksyObj_naive_run.rds'))
# ------ end clustering block ---------


reorder_genes <- function(bank){
  if (is.list(bank@own.expr)) {
    x<-bank@own.expr[[1]]
    d_gene <- dist(as.matrix(x))
    hc_gene <- hclust(d_gene)
    bank@own.expr <- lapply(bank@own.expr, function(x) x[hc_gene$order,])
    bank@nbr.expr <- lapply(bank@nbr.expr, function(x) x[hc_gene$order,])
    m1.list <- lapply(bank@harmonics, function(x) list(m1 = x$m1[hc_gene$order,]))
  } else {
    x<-bank@own.expr
    d_gene <- dist(as.matrix(x))
    hc_gene <- hclust(d_gene)
    bank@own.expr <- bank@own.expr[hc_gene$order,]    
    bank@nbr.expr <- bank@nbr.expr[hc_gene$order,] 
    bank@harmonics$m1 <- bank@harmonics$m1[hc_gene$order,]
  }
  return(bank)
}

non_ambig = bank@meta.data$cell_ID[which(bank@meta.data$Cell_class != 'Ambiguous')]
bank = SubsetBanksy(bank, cells = non_ambig)
nonspatial.main.run = 'clust_M0_lam0_k50_res0.5'
spatial.main.run = 'clust_M0_lam0.2_k50_res0.5'
moffitt.labels = 'clust_major_num'
bank = ConnectClusters(bank = bank, map.to = spatial.main.run)
num_clusters<-max(bank@meta.data[,clust.names(bank)])
hypo.cols<-Banksy:::getPalette(num_clusters)
names(hypo.cols)<-1:num_clusters
hypo.cols.old = hypo.cols
scales::show_col(hypo.cols)

bank<- reorder_genes(bank)
set.seed(42);
sampled.cells = bank@meta.data$cell_ID[sample.int(nrow(bank@meta.data),
                                                       min(80000, nrow(bank@meta.data)))]
bank.sampled = SubsetBanksy(bank, cells = sampled.cells)

hypo.cols['3'] = '#848482'
hypo.cols['5'] = '#008856'
hypo.cols['7'] = '#BE0032'
hypo.cols['8'] = '#F38400'
hypo.cols['10'] = '#F99379'
hypo.cols['11'] = '#5F6B6D' 
hypo.cols['13'] = '#0067A5'
hypo.cols['15'] = '#F6A6A0'

bank.runs = c('clust_M0_lam0_k50_res0.5', 
  'clust_M0_lam0.2_k50_res0.5','clust_M0_lam0.4_k50_res0.5', 
  'clust_M0_lam0.6_k50_res0.5','clust_M0_lam0.8_k50_res0.5',
  'clust_M0_lam1_k50_res0.5') 
bank.reductions<- paste0('umap_M0_lam', as.numeric(gsub('.*lam|_.*', '', bank.runs)))

png(paste0(results.dir, '/umap_banksy.png'), units = 'in',  res = 200,  
    height = 9*1.25, width = 40*1.25)
umapdims <- mapply(FUN = function(reduction, by, title) plotReduction(bank.sampled, 
                                                                      reduction = reduction, 
                                                                      main = title, legend = TRUE,
                                                                      by = by, col.discrete = hypo.cols,
                                                                      type = 'discrete', 
                                                                      pt.size = 0.2, main.size = 35), 
                   reduction = as.list(bank.reductions), 
                   by = as.list(c(bank.runs[2], bank.runs[2])), 
                   title = as.list(c('Non-spatial UMAP space', 'BANKSY UMAP space lam = 0.2', 
                                     'BANKSY UMAP space lam = 0.4','BANKSY UMAP space lam = 0.6',
                                     'BANKSY UMAP space lam = 0.8','BANKSY UMAP space lam = 1.0')),
                   SIMPLIFY = FALSE)
do.call("grid.arrange", c(umapdims, ncol = length(bank.runs)))
dev.off()

# spatial plot for animal 1. 
bank.animal1 = copy(bank)
animal_id = 'Animal_1'
bank.animal1@own.expr = bank.animal1@own.expr[[animal_id]]
bank.animal1@nbr.expr = bank.animal1@nbr.expr[[animal_id]]
bank.animal1@harmonics = bank.animal1@harmonics[[animal_id]]
bank.animal1@cell.locs = bank.animal1@cell.locs[[animal_id]]
bank.animal1@meta.data = bank.animal1@meta.data[bank.animal1@meta.data$dataset %in% c(animal_id), ]
bank.animal1@reduction$pca_M1_lam0$x <- bank.animal1@reduction$pca_M1_lam0$x[grep(paste0(animal_id, '_'), 
                                                                                  rownames(bank.animal1@reduction$pca_M1_lam0$x) ),
]
bank.animal1@reduction$pca_M1_lam0.2$x <- bank.animal1@reduction$pca_M1_lam0.2$x[grep(paste0(animal_id, '_'), 
                                                                                      rownames(bank.animal1@reduction$pca_M1_lam0.2$x) ),
]
bank.animal1@reduction$umap_M1_lam0 <- bank.animal1@reduction$umap_M1_lam0[grep(paste0(animal_id, '_'), 
                                                                                rownames(bank.animal1@reduction$umap_M1_lam0)),
]
bank.animal1@reduction$umap_M1_lam0.2 <- bank.animal1@reduction$umap_M1_lam0.2[grep(paste0(animal_id, '_'), 
                                                                                    rownames(bank.animal1@reduction$umap_M1_lam0.2)),
]


titles = c('Non-spatial UMAP space', 'BANKSY lam = 0.2', 
           'BANKSY lam = 0.4','BANKSY lam = 0.6',
           'BANKSY lam = 0.8','BANKSY lam = 1.0')
bank.runs = c('clust_M0_lam0_k50_res0.5', 
              'clust_M0_lam0.2_k50_res0.5','clust_M0_lam0.4_k50_res0.5', 
              'clust_M0_lam0.6_k50_res0.5','clust_M0_lam0.8_k50_res0.5',
              'clust_M0_lam1_k50_res0.5')

names(titles) = bank.runs

spatial_plots2 <- function(x, colorscheme = hypo.cols) {
  for (runid in bank.runs){
    png(paste0(results.dir, '/spatial_',runid,'.png'), height = 7, width = 7, 
        units = 'in', res = 300)
    layer_IDs <- unique(bank.animal1@meta.data$Bregma)
    # print(titles[runid])
    spatdims<-vector(mode = "list", length = length(layer_IDs))
    for (jj in 1:length(layer_IDs)){
      layer_id<-layer_IDs[jj]
      bank_layer<- SubsetBanksy(x,
                                cells = x@meta.data$cell_ID[x@meta.data$Bregma %in% layer_id]) # this works...
      spatdims[[jj]]<-plotSpatial(bank_layer, type = 'discrete',
                                  by = runid,
                                  col.discrete = colorscheme,
                                  pt.size = 0.2,legend = FALSE, 
                                  # main = paste0(layer_id, 'mm'),
                                  main.size = 15)+#facet_wrap(~feature, ncol = 5)+ 
        scale_x_discrete(labels = NULL, breaks = NULL) + labs(x = NULL) +
        scale_y_discrete(labels = NULL, breaks = NULL) + labs(y = NULL) +
        theme(axis.title.y=element_text(size = 50))
    }
    do.call("grid.arrange", c(spatdims, ncol = 4))
    dev.off()
  }
}
spatial_plots2(bank.animal1, hypo.cols)


