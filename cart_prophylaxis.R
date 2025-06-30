library(Rmagic)
library(plyr)
library(Seurat)
library(tidyverse)
library(cowplot)
library(ggplot2)
library(MAST)
library(DESeq2)
library(EnhancedVolcano)
library(limma)
library(scales)
library(metR)
library(ggpubr)
library(rstatix)
library(svglite)
library(viridis)
library(reshape)
library(harmony)
library(nichenetr)
library(RColorBrewer)
library(Libra)
library(Nebulosa)


#load data
{UPAR <- Read10X(data.dir = "/Users/vyom/data/Aging_cart_preventative/Amor_CAV13_UPAR1/outs/filtered_feature_bc_matrix/")
UT <- Read10X(data.dir = "/Users/vyom/data/Aging_cart_preventative/Amor_CAV13_UT2/outs/filtered_feature_bc_matrix/")

UPAR <- CreateSeuratObject(UPAR, project = "UPAR")
UT <- CreateSeuratObject(UT, project = "UT")}


#filter daata
{
d <- merge(UPAR, y = c(UT), add.cell.ids = c("UPAR","UT"), project = "Sen_preventative")
d[["percent.mt"]] <- PercentageFeatureSet(d, pattern = "mt-")
VlnPlot(d, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),pt.size = 0 ,ncol = 3)
d <- subset(d, subset = nCount_RNA > 1000 & nCount_RNA < 75000 & nFeature_RNA > 2500 & nFeature_RNA < 8000 & percent.mt < 12)

}

# Normilization and scaling
{
Data.list <- SplitObject(d, split.by = "ident")
Data.list <- Data.list[c("UPAR","UT")]
for (i in 1:length(Data.list)) {
  
  Data.list[[i]] <- SCTransform(Data.list[[i]], verbose = FALSE)
}

#select highly variable genes 
Data.features <- SelectIntegrationFeatures(object.list = Data.list, nfeatures = 2500)
Data.list <- PrepSCTIntegration(object.list = Data.list, anchor.features = Data.features, 
                                verbose = FALSE)
Data.anchors <- FindIntegrationAnchors(object.list = Data.list, normalization.method = "SCT", 
                                       anchor.features = Data.features, verbose = FALSE)

pre_obj <- IntegrateData(anchorset = Data.anchors, normalization.method = "SCT", 
                         verbose = TRUE)}

# Visulization and Clustering
{pre_obj <- RunPCA(pre_obj)

VizDimLoadings(pre_obj, dims = 1:2, reduction = "pca")

DimPlot(pre_obj, reduction = "pca")
ElbowPlot(pre_obj, ndims = 100, reduction = "pca")

pre_obj <- RunUMAP(pre_obj, dims = 1:25, n.epochs = 250, n.neighbors = 10)
pre_obj <- FindNeighbors(pre_obj, dims = 1:25)
pre_obj <- FindClusters(pre_obj, resolution = 1)
DimPlot(pre_obj, reduction = "umap", label = TRUE)
DimPlot(pre_obj, reduction = "umap",group.by = "Treatment", label = TRUE)
FeaturePlot(pre_obj, features = 'nCount_RNA')

new.cluster.ids <- c("Stem", "TA", "Ent.", "EP", "Ent.", "TA", "Ent.", "EP", "Gob.", "T", "Ent.", "Stem", "Stem", "T", "Gob.", "Tuft", "End.", "Ent.", "End.", "DC", "Pan.", "Gob.", "T")
pre_obj[["Cell_Type"]] <- Idents(pre_obj)
names(new.cluster.ids) <- levels(pre_obj)
pre_obj <- RenameIdents(pre_obj, new.cluster.ids)
pre_obj[["Cell_Type"]] <- Idents(pre_obj)
DimPlot(pre_obj, reduction = "umap", group.by= 'Cell_Type')}



#magic imputation
{
library(reticulate)
use_python("/Users/vyom/miniconda3/bin/python")
DefaultAssay(pre_obj) <- "SCT"
pre_obj <- magic(pre_obj, assay = 'SCT')
}

#score calculation
{All_Genes <- pre_obj@assays$RNA@data@Dimnames[[1]]

sen_mayo <- c("Acvr1b", "Ang", "Angpt1", "Angptl4", "Areg", "Axl", "Bex3", "Bmp2", "Bmp6",
              "C3", "Ccl1", "Ccl2", "Ccl20", "Ccl24", "Ccl26", "Ccl3", "Ccl4", "Ccl5", "Ccl7",
              "Ccl8", "Cd55", "Cd9", "Csf1", "Csf2", "Csf2rb", "Cst10", "Ctnnb1", "Ctsb",
              "Cxcl1", "Cxcl10", "Cxcl12", "Cxcl16", "Cxcl2", "Cxcl3", "Cxcr2", "Dkk1", "Edn1",
              "Egf", "Egfr", "Ereg", "Esm1", "Ets2", "Fas", "Fgf1", "Fgf2", "Fgf7", "Gdf15",
              "Gem", "Gmfg", "Hgf", "Hmgb1", "Icam1", "Icam5", "Igf1", "Igfbp1", "Igfbp2",
              "Igfbp3", "Igfbp4", "Igfbp5", "Igfbp6", "Igfbp7", "Il10", "Il13", "Il15", "Il18",
              "Il1a", "Il1b", "Il2", "Il6", "Il6st", "Il7", "Inha", "Iqgap2", "Itga2", "Itpka",
              "Jun", "Kitl", "Lcp1", "Mif", "Mmp13", "Mmp10", "Mmp12", "Mmp13", "Mmp14", "Mmp2",
              "Mmp3", "Mmp9", "Nap1l4", "Nrg1", "Pappa", "Pecam1", "Pgf", "Pigf", "Plat", "Plau",
              "Plaur", "Ptbp1", "Ptger2", "Ptges", "Rps6ka5", "Scamp4", "Selplg", "Sema3f",
              "Serpinb3a", "Serpine1", "Serpine2", "Spp1", "Spx", "Timp2", "Tnf", "Tnfrsf11b",
              "Tnfrsf1a", "Tnfrsf1b", "Tubgcp2", "Vegfa", "Vegfc", "Vgf", "Wnt16", "Wnt2"
)

All_Genes <- pre_obj@assays$SCT@data@Dimnames[[1]]
sen_mayo <- intersect(All_Genes, sen_mayo)

mean.exp <- zscore(colMeans(x = pre_obj@assays$MAGIC_SCT@data[sen_mayo, ], na.rm = TRUE), dist = 'norm')
if (all(names(x = mean.exp) == rownames(x = pre_obj@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  pre_obj@meta.data$sen_mayo <- mean.exp
  

Gene.list <- c('Lgr4','Myc', 'Sox9', 'Olfm4', 'Malat1', 'Hopx', 'Ccnd1' )
mean.exp <- zscore(colMeans(x = pre_obj@assays$MAGIC_SCT@data[Gene.list, ], na.rm = TRUE), dist ='norm')
if (all(names(x = mean.exp) == rownames(x = pre_obj@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  pre_obj@meta.data$Stem_Score <- mean.exp
}  
  
}
}

#Violin plots
for(i in gene.list) {
  selected_cells <- names(pre_obj$celltype)
  vln_data <- FetchData(pre_obj,
                        vars = c(i,"Treatment", "celltype"),
                        cells = selected_cells,
                        slot = "data")
  vln_data$Type
  vln_data <- melt(vln_data)
  
  All <- VlnPlot(pre_obj, split.by = "Treatment",group.by = 'celltype', features = i, pt.size = 0, assay = "MAGIC_SCT",  cols = c('#7CA1CC' ,'#FF4902'), log = TRUE, split.plot = TRUE) + 
    theme(legend.position = 'none') + 
    geom_boxplot(width=0.3,outlier.shape = NA, coef = 0, lwd= .2) + xlab('') + 
    theme(text = element_text(size=9), axis.text.x = element_text(size = 9), axis.text.y = element_text(size =9)) +
    stat_compare_means(data= vln_data, aes(x = celltype, y = value, fill = Treatment), method = "wilcox.test",  label = "p.signif", size = 2, hide.ns = TRUE) #+ scale_y_continuous(expand = expansion(mult = c(0, .1)))
  All$layers[[1]]$aes_params$size = .15
  ggsave(file = paste0('CART_Score_', i, '.pdf'), plot=All, width=3, height=3, units="in")
}

{
  genes_stem_prolif <- unique(c('Egf',	'Wnt9b',	'Wnt3',	'Bambi', 'Lgr5',"Ascl2",'Olfm4'))
  genes_antimicrobial <-  unique(c('DEFA5',	'DEFA6',	'DEFA27',	'DEFA30',	'DEFA26',	'DEFA20',	'MPTX1',	'MPTX2',	'LYZ1',	'LYZ2',	'KLF6',	'KLF4',	'KLF15',	'PLA2G2A',	'REG3B',	'REG3G',	'ITLN1',	'SPINK4',	'MT1',	'FABP1',	'FABP2',	'FABP4',	'MMP7',	'PLA2G12A', 'MUC2',	'CLCA1',	'PDIA4',	'PDIA5',	'ZG16',	'REG4',	'REG3B',	'ANG4', 'CASP6',	'CASP8',	'REG1',	'REG1',	'REG3A',	'REG3B',	'REG4'))%>% convert_human_to_mouse_symbols()         
  genes_hormone_secretion <-  unique(c('CHGA',	'GCG',	'GHRL',	'TRPA1',	'NTS',	'AFP',	'VEGFA',	'PCSK1',	'SCN3A',	'NEUROG3',	'NEUROD1',	'NKX2-2',	'FABP5',	'FABP1',	'FABP6',	'GIP',	'GLP1R',	'GLP2R',	'PYY',	'GAST',	'SCGN'))%>% convert_human_to_mouse_symbols()
  genes_nutrient_absorption <-  unique(c('SLC7A8',	'FABP1',	'CPS1',	'VDR',	'APOA4',	'RBP2',	'GSTA1',	'SLC5A1',	'SLC7A9',	'EHF',	'PPP1R14B',	'PP1R14D')) %>% convert_human_to_mouse_symbols()
  
  genes_stem_prolif <- genes_stem_prolif[!is.na(genes_stem_prolif)]
  genes_antimicrobial <- genes_antimicrobial[!is.na(genes_antimicrobial)]
  genes_hormone_secretion <- genes_hormone_secretion[!is.na(genes_hormone_secretion)]
  genes_nutrient_absorption <- genes_nutrient_absorption[!is.na(genes_nutrient_absorption)]
  
  All_Genes <- pre_obj@assays$SCT@data@Dimnames[[1]]
  genes_stem_prolif <- intersect(All_Genes, genes_stem_prolif)
  genes_antimicrobial <- intersect(All_Genes, genes_antimicrobial)   
  genes_hormone_secretion <- intersect(All_Genes, genes_hormone_secretion)  
  genes_nutrient_absorption <- intersect(All_Genes, genes_nutrient_absorption)
  
  pre_obj <- AddModuleScore(object = pre_obj, features = list(genes_stem_prolif), name = 'stem_prolif')
  pre_obj <- AddModuleScore(object = pre_obj, features = list(genes_antimicrobial), name = 'antimicrobial')
  pre_obj <- AddModuleScore(object = pre_obj, features = list(genes_hormone_secretion), name = 'hormone_secretion')
  pre_obj <- AddModuleScore(object = pre_obj, features = list(genes_nutrient_absorption), name = 'nutrient_absorption')
  
  Idents(pre_obj) <- pre_obj$Cell_Type
  Cell_Types <- c('Ent.', 'End.', 'Gob.', 'Pan.', 'Tuft')
  stem_prolif_FC <- list()
  nutrient_absorption_FC <- list()
  antimicrobial_FC <- list()
  hormone_secretion_FC <- list()
  
  for(i in Cell_Types){
    stem_prolif_FC[[i]] <- FoldChange(pre_obj,ident.1 = "UPAR", ident.2 = "UT",features = genes_stem_prolif, group.by = 'Treatment', subset.ident = i)
  }
  stem_prolif_FC_aggr <- c()
  for(i in Cell_Types){
    stem_prolif_FC_aggr <- c(stem_prolif_FC_aggr, mean(as.numeric(stem_prolif_FC[[i]]$avg_log2FC)))
  }
  
  for(i in Cell_Types){
    nutrient_absorption_FC[[i]] <- FoldChange(pre_obj,ident.1 = "UPAR", ident.2 = "UT",features = genes_nutrient_absorption, group.by = 'Treatment', subset.ident = i)
  }
  nutrient_absorption_FC_aggr <- c()
  for(i in Cell_Types){
    nutrient_absorption_FC_aggr <- c(nutrient_absorption_FC_aggr, mean(as.numeric(nutrient_absorption_FC[[i]]$avg_log2FC)))
  }
  
  for(i in Cell_Types){
    antimicrobial_FC[[i]] <- FoldChange(pre_obj,ident.1 = "UPAR", ident.2 = "UT",features = genes_antimicrobial, group.by = 'Treatment', subset.ident = i)
  }
  antimicrobial_FC_aggr <- c()
  for(i in Cell_Types){
    antimicrobial_FC_aggr <- c(antimicrobial_FC_aggr, mean(as.numeric(antimicrobial_FC[[i]]$avg_log2FC)))
  }
  
  for(i in Cell_Types){
    hormone_secretion_FC[[i]] <- FoldChange(pre_obj,ident.1 = "UPAR", ident.2 = "UT",features = genes_hormone_secretion, group.by = 'Treatment', subset.ident = i)
  }
  hormone_secretion_FC_aggr <- c()
  for(i in Cell_Types){
    hormone_secretion_FC_aggr <- c(hormone_secretion_FC_aggr, mean(as.numeric(hormone_secretion_FC[[i]]$avg_log2FC)))
  }
  
  
  Module_Heatmap <- data.frame(nutrient_absorption_FC_aggr, stem_prolif_FC_aggr, hormone_secretion_FC_aggr,antimicrobial_FC_aggr)
  rownames(Module_Heatmap) <- Cell_Types
  colnames(Module_Heatmap) <- c('nutrient_absorption','stem_prolif','hormone_secretion','antimicrobial')
  Module_Heatmap$Cell_Type <- rownames(Module_Heatmap)
  Module_Heatmap1 <-melt(Module_Heatmap)
  
  scores <- c('nutrient_absorption1','stem_prolif1','hormone_secretion1','antimicrobial1')
  Idents(pre_obj) <- pre_obj$Cell_Type
  pvals <- c()
  for(j in scores){
    for(i in Cell_Types){
      selected_cells <- WhichCells(object = pre_obj ,idents = c(i))
      vln_data <- FetchData(pre_obj,
                            vars = c(j,"Treatment", "Cell_Type"),
                            cells = selected_cells,
                            slot = "data")
      vln_data <- melt(vln_data)
      vln_data <- vln_data[c('Treatment', 'Cell_Type', 'value')]
      
      pvals <- c(pvals, t.test(vln_data[which(vln_data$Treatment == 'UPAR'),]$value, vln_data[which(vln_data$Treatment == 'UT'),]$value,  alternative = c("greater"))$p.value)
      print(i)
    }
  }
  
  Module_Heatmap1$pval <- -log(pvals) 
  
  my_levels <- c('Ent.', 'End.', 'Gob.', 'Pan.', 'Tuft')
  Module_Heatmap1$Cell_Type <- factor(Module_Heatmap1$Cell_Type, levels = my_levels)
  min(Module_Heatmap1$value)
  pathwayColorsDiff = rev(brewer.pal(9, "RdBu"))
  ggplot(Module_Heatmap1, aes(Cell_Type, variable, fill= value)) + geom_tile() + 
    scale_fill_gradientn(name = "Log2FC", colors = pathwayColorsDiff, limits = c(-1,1)) + coord_flip() +
    theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45), axis.text.y =   element_text( vjust = 0.1, hjust = 1), panel.background = element_blank(), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())
  
  ggsave(file="Preventative_CART_score_heatmap.pdf",  width=3, height=3, units="in")
  
  ggplot(Module_Heatmap1, aes(x = Cell_Type, y = variable, size = pval ,color= value)) + geom_point() + theme_vyom + 
    scale_color_gradientn(name = "Log2FC", colors = pathwayColorsDiff, limits = c(-1,1)) + coord_flip() + labs(size="-log(p)") +
    theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45), axis.text.y =   element_text( vjust = 0.1, hjust = 1), panel.background = element_blank(), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())
  
  ggsave(file="Preventative_CART_score_dotplot.pdf",  width=3, height=3, units="in")
}

#pseudotime analysis
{library(SeuratWrappers)
library(monocle3)
library(Signac)
library(org.Mm.eg.db)
gene_symbol <- as.list(org.Mm.egSYMBOL)
Idents(pre_obj) <- pre_obj$Cell_Type
pseudo_stem <- WhichCells(object = pre_obj, idents = 'Stem')


cds <- as.cell_data_set(pre_obj)
cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(pre_obj[["RNA"]])
cds <- cluster_cells(cds = cds, reduction_method = "UMAP")
cds <- learn_graph(cds, use_partition = FALSE,   close_loop = FALSE)
cds <- order_cells(cds, reduction_method = "UMAP")

plot_cells(cds,
           color_cells_by = "Cell_Type",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

pre_obj <- AddMetaData(
  object = pre_obj,
  metadata = cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "Pseudotime"
)

#pseudotime_gene
UPAR_rownames <- rownames(pre_obj@meta.data[pre_obj@meta.data$Treatment == 'UPAR',] )
Control_rownames <- rownames(pre_obj@meta.data[pre_obj@meta.data$Treatment == 'UT',] )
Cell_type_frame <- data.frame(pre_obj$Treatment,pre_obj$orig.ident )
all <-  cbind(type = Cell_type_frame, Pseudotime = pre_obj$Pseudotime)
colnames(all) <- c("Treatment",'Replicate',"Pseudotime")

all <- all %>% mutate(quantile = ntile(all$Pseudotime, 10))
quart <- data.frame(var1 = tapply(all$Pseudotime, all$quantile, mean), var2 = tapply(all$Pseudotime, all$quantile, max))
for(i in gene.list){
  colnames(pre_obj)
  DefaultAssay(pre_obj)
  
  score_iter = i
  #for gene
  #quart$maxi = max(pre_obj@assays$MAGIC_RNA@data[i, ])
  #for score
  quart$maxi = max(pre_obj$Stem_Score)
  
  UPAR <- subset(all, Treatment=='UPAR')
  UT <- subset(all, Treatment=='UT')
  
  #for gene
  #UPAR <-  cbind(UPAR, pre_obj[,UPAR_rownames]@assays$MAGIC_RNA@data[i, ])
  #UT <-  cbind(UT, pre_obj[,Control_rownames]@assays$MAGIC_RNA@data[i, ])
  
  #for score
  UPAR <-  cbind(UPAR, pre_obj[,UPAR_rownames]$Stem_Score)
  UT <-  cbind(UT, pre_obj[,Control_rownames]$Stem_Score)
  
  colnames(UPAR) <- c("Sample",'Replicate',"Pseudotime",'quantile' ,'score')
  colnames(UT) <- c("Sample",'Replicate',"Pseudotime",'quantile' ,'score')
  
  p.vals.prop = ''
  for(j in 1:10){
    UPARinqunat  <- subset(UPAR, quantile== j)
    Cinqunat  <- subset(UT, quantile== j)
    
    p.vals.prop = c(p.vals.prop, wilcox.test(UPARinqunat$score, Cinqunat$score, alternative = "two.sided")$p.value)
  }
  p.vals.prop = p.vals.prop[2:11]
  p.vals.prop <- as.numeric(p.vals.prop)
  p.vals.prop[p.vals.prop < 0.05] =  "*"
  p.vals.prop[p.vals.prop > 0.05] =  ""
  
  
  
  quart$p_asterik = p.vals.prop
  
  
  Pseudotimedot <- ggplot() +
    # Points
    #ggrastr::rasterise(geom_jitter(data=UT, aes(x=Pseudotime, y=score), size=0, alpha=0.05, width=0,stroke = 0, colour = "#7CA1CC" ), dpi = 400) + 
    #ggrastr::rasterise(geom_jitter(data=UPAR, aes(x=Pseudotime, y=score), size=0, alpha=0.05, width=0,stroke = 0, colour = "#FF4902"), dpi = 400) +
    # lines
    geom_smooth(data=UT, aes(x=Pseudotime, y=score), fill="black", colour="#7CA1CC", size=.5) + 
    geom_smooth(data=UPAR, aes(x=Pseudotime, y=score), fill="black", colour="#FF4902", size=.5) +
    theme( panel.background = element_rect(fill = "white", colour = "Black")) + ylab(score_iter) + theme_vyom +
    geom_vline(data = quart, aes(xintercept = var1), size = .25) + geom_text(data = quart, aes(var2 , maxi, label = p_asterik, fontface="bold"),size = 1.75 , position=position_dodge2(0.5), vjust=.75, inherit.aes = FALSE) + 
    theme(text = element_text(size=1),axis.text = element_text(size = 2) ,plot.title = element_text(size = 6), axis.text.x = element_text(size = 4), axis.text.y = element_text(size = 4), axis.title.x = element_text(size = 6), axis.title.y = element_text(size = 6)) +
    scale_x_continuous(expand = c(0, 0))
  Pseudotimedot
  ggsave(file= paste0(i,"_Prevention_UPAR_scrna_Pseudotime_expression.pdf"),plot = Pseudotimedot, width=3, height=2)
}

#pseudotime density plots
{
  dimredP = data.frame(pre_obj@reductions$umap@cell.embeddings)
  dimredP$Cell = rownames(dimredP)
  all.equal(names(pre_obj$seurat_clusters), dimredP$Cell)
  all.equal(names(pre_obj$Treatment), dimredP$Cell)
  
  dimredP$Cluster = pre_obj$Cell_Type # Overwrite clusters with cell types for current plots
  dimredP$Type = pre_obj$Treatment
  dimredP$CellType = dimredP$Cluster
  dimredP$Pseudotime = pre_obj$Pseudotime
  
  
  i = 'Stem'
  cell = i
  dimredP <- dimredP %>% mutate(quantile = ntile(Pseudotime, 10))
  dimredaub = dimredP %>% filter(CellType == i)
  quart <- data.frame(var1 = tapply(dimredP$Pseudotime, dimredP$quantile, mean), var2 = tapply(dimredP$Pseudotime, dimredP$quantile, max))
  p.vals.prop = ''
  dimredc <- dimredaub %>% filter(Type == 'UT')
  dimreda <- dimredaub %>% filter(Type == 'UPAR')
  
  for(i in 1:10){
    prop_table1 <- c(sum(dimreda$quantile == i), sum(dimredc$quantile == i))
    prop_table2 <- c(length(dimreda$quantile) - sum(dimreda$quantile == i), length(dimredc$quantile) - sum(dimredc$quantile == i)  )
    
    prop_table_fin <- cbind(prop_table1, prop_table2)
    p.vals.prop = c(p.vals.prop, fisher.test(prop_table_fin)$p.value)
    p.vals.prop = p.vals.prop[2:11]
  }
  p.vals.prop <- p.adjust(p.vals.prop, method = 'BH')
  quart$p_asterik = ''
  quart$p_asterik[p.vals.prop < 0.05] = "*"
  quart$p_asterik[p.vals.prop < 0.01] = "**"
  quart$p_asterik[p.vals.prop < 0.001] = "***"
  quart$p_asterik[is.na(quart$p_asterik)] <- ""
  quart$p_asterik
  
  plot11 <- dimredaub %>% group_by(Type) %>% 
    # calculate densities for each group over same range; store in list column
    summarise(d = list(density(Pseudotime, from = min(.$Pseudotime), to = max(.$Pseudotime)))) %>% 
    # make a new data.frame from two density objects
    do(data.frame(x = .$d[[1]]$x,    # grab one set of x values (which are the same)
                  y = .$d[[2]]$y - .$d[[1]]$y))
  quart$maxi = max(plot11$y)
  quart$maxi2 = .125
  
  Density_diff_stem_1 <- ggplot(plot11, aes(x, y, fill=y>=0)) + geom_col(position = 'dodge', width = .1) + ggtitle('Stem') + theme_vyom + scale_fill_manual(values = c( '#7CA1CC' ,'#FF4902'))+ theme(legend.position = "none") + ylab('Density') + xlab('Pseudotime') + geom_vline(data = quart, aes(xintercept = var2), size = .075) + geom_text(data = quart, aes(var1 , maxi, label = p_asterik), position=position_dodge2(0.2), size = 1.5, vjust=1, inherit.aes = FALSE) + theme(text = element_text(size=3),axis.text = element_text(size = 2) ,plot.title = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), plot.margin = unit(c(.2, 0, 0, 0), "cm"))
  Density_diff_stem_1
  ggsave(paste0("CART_UPAR_vs_UT_prevention_Pseduodensity_stem.pdf"), plot = Density_diff_stem_1, width = 3, height = 2, units="in")	
  
  #all cells
  
  dimredaub = dimredP
  quart <- data.frame(var1 = tapply(dimredP$Pseudotime, dimredP$quantile, mean), var2 = tapply(dimredP$Pseudotime, dimredP$quantile, max))
  p.vals.prop = ''
  dimredc <- dimredaub %>% filter(Type == 'UT')
  dimreda <- dimredaub %>% filter(Type == 'UPAR')
  
  for(i in 1:10){
    prop_table1 <- c(sum(dimreda$quantile == i), sum(dimredc$quantile == i))
    prop_table2 <- c(length(dimreda$quantile) - sum(dimreda$quantile == i), length(dimredc$quantile) - sum(dimredc$quantile == i)  )
    
    prop_table_fin <- cbind(prop_table1, prop_table2)
    p.vals.prop = c(p.vals.prop, fisher.test(prop_table_fin)$p.value)
    p.vals.prop = p.vals.prop[2:11]
  }
  p.vals.prop <- p.adjust(p.vals.prop, method = 'BH')
  quart$p_asterik = ''
  quart$p_asterik[p.vals.prop < 0.05] = "*"
  quart$p_asterik[p.vals.prop < 0.01] = "**"
  quart$p_asterik[p.vals.prop < 0.001] = "***"
  quart$p_asterik[is.na(quart$p_asterik)] <- ""
  quart$p_asterik
  dimredaub$Type
  plot11 <- dimredaub %>% group_by(Type) %>% 
    # calculate densities for each group over same range; store in list column
    summarise(d = list(density(Pseudotime, from = min(.$Pseudotime), to = max(.$Pseudotime)))) %>% 
    # make a new data.frame from two density objects
    do(data.frame(x = .$d[[1]]$x,    # grab one set of x values (which are the same)
                  y = .$d[[2]]$y - .$d[[1]]$y))
  quart$maxi = max(plot11$y)
  quart$maxi2 = .075
  
  Density_diff_all <- ggplot(plot11, aes(x, y, fill=y>=0)) + geom_col(position = 'dodge', width = .0978619) + ggtitle('All Cells') + theme_vyom + scale_fill_manual(values = c( '#7CA1CC' ,'#FF4902'))+ theme(legend.position = "none") + ylab('Density') + xlab('Pseudotime') + geom_vline(data = quart, aes(xintercept = var2), size = .075) + geom_text(data = quart, aes(var1 , maxi, label = p_asterik), position=position_dodge2(0.2), size = 1.5, vjust=1, inherit.aes = FALSE) + theme(text = element_text(size=3),axis.text = element_text(size = 2) ,plot.title = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), plot.margin = unit(c(.2, 0, 0, 0), "cm"))
  ggsave(paste0("CART_UPAR_vs_UT_prevention_Pseduodensity_All.pdf"), plot = Density_diff_all, width = 3, height = 2, units="in")	
  
}

}

