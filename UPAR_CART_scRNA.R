# Vyom Shah
# UPAR vs UT in old and young patients scRNA-seq

#load packages
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

#plotting themes
theme_linedraw2 = theme_linedraw() + theme(strip.background=element_rect(fill="grey80", colour="grey50", size=0.2), strip.text.x=element_text(colour="black"), strip.text.y=element_text(colour="black"))
theme_vyom = theme_linedraw2 + theme(legend.position="right", legend.title=element_text(size=15), legend.text=element_text(size=14), axis.text.x = element_text(size=12, angle=-90, hjust=0, vjust=0.5), axis.text.y=element_text(size=12), axis.title=element_text(size=15), axis.title.y=element_text(vjust=1), plot.title = element_text(size=18, vjust=1.5), strip.background = element_rect(fill="#EEEEEE"), strip.text = element_text(size = 11), panel.grid.major = element_line(colour = "grey98"), panel.grid.minor = element_blank())

# Load in the Data
uPAR_F_old <- Read10X(data.dir = "/Users/vyom/data/carT_Intestine/Amor_CAV05_C3NH_C3RH/outs/filtered_feature_bc_matrix/")
UT_F_old <- Read10X(data.dir = "/Users/vyom/data/carT_Intestine/Amor_CAV05_C32RH_C32LH/outs/filtered_feature_bc_matrix/")
uPAR_M_old <- Read10X(data.dir = "/Users/vyom/data/carT_Intestine/Amor_CAV06_C1RH_C1LH/outs/filtered_feature_bc_matrix/")
UT_M_old <- Read10X(data.dir = "/Users/vyom/data/carT_Intestine/Amor_CAV06_C12RH/outs/filtered_feature_bc_matrix/")

uPAR_F_young <- Read10X(data.dir = "/Users/vyom/data/carT_Intestine/cart_young/Amor_CAV07_C5NH_C5LH/outs/filtered_feature_bc_matrix/")
UT_F_young <- Read10X(data.dir = "/Users/vyom/data/carT_Intestine/cart_young/Amor_CAV07_C52RH_C52LH/outs/filtered_feature_bc_matrix/")
uPAR_M_young <- Read10X(data.dir = "/Users/vyom/data/carT_Intestine/cart_young/Amor_CAV08_C6NH/outs/filtered_feature_bc_matrix/")
UT_M_young <- Read10X(data.dir = "/Users/vyom/data/carT_Intestine/cart_young/Amor_CAV08_C62RH/outs/filtered_feature_bc_matrix/")


uPAR_F_old <- CreateSeuratObject(uPAR_F_old, project = "uPAR_F_old")
UT_F_old <- CreateSeuratObject(UT_F_old, project = "UT_F_old")
uPAR_M_old <- CreateSeuratObject(uPAR_M_old, project = "uPAR_M_old")
UT_M_old <- CreateSeuratObject(UT_M_old, project = "UT_M_old")
uPAR_F_young <- CreateSeuratObject(uPAR_F_young, project = "uPAR_F_young")
UT_F_young <- CreateSeuratObject(UT_F_young, project = "UT_F_young")
uPAR_M_young <- CreateSeuratObject(uPAR_M_young, project = "uPAR_M_young")
UT_M_young <- CreateSeuratObject(UT_M_young, project = "UT_M_young")

#merge and filter
d <- merge(uPAR_F_old, y = c(UT_F_old,uPAR_M_old,UT_M_old,uPAR_F_young,UT_F_young,uPAR_M_young,UT_M_young), add.cell.ids = c("uPAR_F_old", "UT_F_old","uPAR_M_old","UT_M_old","uPAR_F_young","UT_F_young","uPAR_M_young","UT_M_young"), project = "Immune")
d
d[["percent.mt"]] <- PercentageFeatureSet(d, pattern = "mt-")
VlnPlot(d, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),pt.size = 0 ,ncol = 3)
d <- subset(d, subset = nCount_RNA > 1000 & nCount_RNA < 75000 & nFeature_RNA > 2000 & nFeature_RNA < 8000 & percent.mt < 15)
d
d$orig.ident
Data.list <- SplitObject(d, split.by = "ident")
Data.list <- Data.list[c("uPAR_F_old", "UT_F_old","uPAR_M_old","UT_M_old","uPAR_F_young","UT_F_young","uPAR_M_young","UT_M_young")]
for (i in 1:length(Data.list)) {
  
  Data.list[[i]] <- SCTransform(Data.list[[i]], verbose = FALSE)
}

# Normilization
#select highly variable genes 
Data.features <- SelectIntegrationFeatures(object.list = Data.list, nfeatures = 2500)
Data.list <- PrepSCTIntegration(object.list = Data.list, anchor.features = Data.features, 
                                verbose = FALSE)
Data.anchors <- FindIntegrationAnchors(object.list = Data.list, normalization.method = "SCT", 
                                       anchor.features = Data.features, verbose = FALSE)
cart_obj <- IntegrateData(anchorset = Data.anchors, normalization.method = "SCT", 
                          verbose = TRUE)

# Visulization and Clustering
cart_obj <- RunPCA(cart_obj)

VizDimLoadings(cart_obj, dims = 1:2, reduction = "pca")

DimPlot(cart_obj, reduction = "pca")

#determine ideal PCs for umap
ElbowPlot(cart_obj, ndims = 100, reduction = "pca")

cart_obj <- RunUMAP(cart_obj, dims = 1:15, n.neighbors = 5, n.epochs = 500)
cart_obj <- FindNeighbors(cart_obj, dims = 1:15)

cart_obj <- FindClusters(cart_obj, resolution = 1)
DimPlot(cart_obj, reduction = "umap", label = TRUE)
allmarkers <- FindAllMarkers(cart_obj, min.pct = 0.5)

#annotate clusters
DotPlot_Sig <- c("Lgr5","Ascl2","Olfm4","Gkn3",'Cbx3','Larp1','Slc39a1','Hnf4a','Sox4','Mmp7','Dll1','Tff3',"Tubb5","Syce2","Stmn1","Fbxo5",'Cenpa','Ccna2','Ube2c','Cdkn3',"Apoa1","Apoa4","Fabp1","Adh6a","Chgb","Tac1","Tph1","Neurog3", "Muc2","Fcgbp","Atoh1","Agr2","Pou2f3","Avil","Tuba1a","Adh1","Lyz1","Defa17","Defa24","Ang4") 
DotPlot(cart_obj, features = DotPlot_Sig, assay = 'SCT') + labs(y= "Cell Type", x="") + scale_colour_distiller( palette ="RdYlBu") + scale_size(range = c(0, 1)) +
  theme(text = element_text(size=5), axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust= .01), axis.text.y = element_text(size = 6), axis.title.x = element_blank(), axis.title.y = element_blank())
DotPlot(cart_obj, features = DotPlot_Sig, dot.scale = 10, scale.max= 100, assay = 'SCT') + labs(y= "Cell type", x="") + scale_colour_distiller( palette ="RdYlBu") + scale_size(range = c(0, 1.3)) +
  theme(legend.key.size = unit(.1, "in"), legend.text= element_text(size=4), legend.title = element_text(size=4),axis.line = element_line(size = .3), axis.ticks = element_line(size = .3), text = element_text(size=4), axis.text.x = element_text(size = 4, angle = 90, hjust = 1, vjust= .01), axis.text.y = element_text(size = 4), axis.title.x = element_blank(), axis.title.y = element_blank())
Idents(cart_obj) <- cart_obj$integrated_snn_res.1

cart_obj <- FindClusters(cart_obj, resolution = 1)
new.cluster.ids <- c('Stem', 'Enterocyte (Proximal)', 'Stem', 'Enterocyte Progenitor', 'Enterocyte Progenitor', 'Enterocyte (Proximal)', 'Transit Amplifying', 'Transit Amplifying', 'Enterocyte (Distal)', 'Enterocyte (Distal)', 'Enterocyte (Proximal)', 'Enterocyte (Proximal)', 'Goblet', 'Goblet', 'Transit Amplifying', 'Stem', 'Tuft', 'Enterocyte (Proximal)', '', 'Enterocyte Progenitor', 'Enterocyte Progenitor', 'Enteroendocrine', 'Enterocyte (Proximal)', 'Goblet', 'Enteroendocrine', '', '', 'Goblet', 'Paneth', 'Goblet')     
cart_obj[["Cell_Type"]] <- Idents(cart_obj)
names(new.cluster.ids) <- levels(cart_obj)
cart_obj <- RenameIdents(cart_obj, new.cluster.ids)
cart_obj[["Cell_Type"]] <- Idents(cart_obj)
DimPlot(cart_obj, reduction = "umap", group.by= 'Cell_Type')

#refine immune clusters
Idents(cart_obj) <- cart_obj$Cell_Type
Cart_immune <- subset(cart_obj, idents = '')

Cart_immune <- RunPCA(Cart_immune, verbose = FALSE)

Cart_immune <- RunUMAP(Cart_immune, dims = 1:50, n.neighbors = 5, n.epochs = 500)

Cart_immune <- FindNeighbors(Cart_immune, dims = 1:50)
Cart_immune <- FindClusters(Cart_immune, resolution = .5)
DimPlot(Cart_immune, reduction = "umap", label = TRUE)


DotPlot_Sig <- unique(c('Ptprc',"S100a6","Ly6a","Anxa3", "Areg",'Cd3g','Cd3e','Cd8a','Cd4','Trac','Tcrg-C1','Lag3','Pdcd1','Havcr2','Tox','Tcf7','Gzmb','Tbx21','Ifng','Gata3','Il4','Rorc','Il17a','Il17f','Foxp3','Il10','Il2rb','Il2ra','Klrd1','Cd19','Ighm','Ighg1','Cd74','Ciita','Nrc1','Klre1','Itgam','Itgax','H2-Eb1','H2-Ab1','Arg1','Mrc1','Tgfbi','Ccr2','Vegfa','Prdx1','Clec4d','Ccl5','Cd83','Ccr7','Fcn1','Msrb1','Ly6g','Col3a1','Sparc'))
DotPlot(Cart_immune, features = DotPlot_Sig, dot.scale = 10, scale.max= 100, assay = 'SCT') + labs(y= "Cell type", x="") + scale_colour_distiller( palette ="RdYlBu") + scale_size(range = c(0, 1.3)) 

new.cluster.ids <- c('T Cell','Enterocyte (Distal)','T Cell','T Cell', 'Myeloid', 'B Cell', 'Myeloid', 'T Cell')

Cart_immune[["Cell_Type"]] <- Idents(Cart_immune)
names(new.cluster.ids) <- levels(Cart_immune)
Cart_immune <- RenameIdents(Cart_immune, new.cluster.ids)
Cart_immune[["Cell_Type"]] <- Idents(Cart_immune)
DimPlot(Cart_immune, reduction = "umap", label = TRUE)

cart_obj$Cell_type <- as.character(cart_obj$Cell_Type)
cart_obj$Cell_type[WhichCells(Cart_immune)] <- paste(Idents(Cart_immune))
Idents(cart_obj) <- cart_obj$Cell_type
DimPlot(cart_obj, reduction = "umap", group.by= 'Cell_type')

my_levels <- c('Stem', 'Transit Amplifying', 'Enterocyte Progenitor','Enterocyte (Proximal)','Enterocyte (Distal)', 'Enteroendocrine', 'Goblet', 'Paneth','Tuft', 'T Cell', 'B Cell', 'Myeloid')
cart_obj$Cell_Type <- factor(x = cart_obj$Cell_type, levels = my_levels)
DimPlot(cart_obj, reduction = "umap", group.by= 'Cell_Type')
ggsave(file = paste0('Cart_umap.pdf'), width=7, height=5, units="in")

#rename Idents
levels(factor(cart_obj@meta.data$orig.ident))
cart_obj[["Sex"]] <- Idents(cart_obj)
new.cluster.ids <- c("Female", "Female","Male","Male", "Female", "Female","Male","Male")
names(new.cluster.ids) <- levels(cart_obj)
cart_obj <- RenameIdents(cart_obj, new.cluster.ids)
cart_obj[["Sex"]] <- Idents(cart_obj)

Idents(cart_obj) <- factor(cart_obj$orig.ident)

levels(factor(cart_obj@meta.data$orig.ident))
cart_obj[["Treatment"]] <- Idents(cart_obj)
new.cluster.ids <- c("uPAR", "uPAR", "uPAR", "uPAR", "UT", "UT", "UT", "UT")
names(new.cluster.ids) <- levels(cart_obj)
cart_obj <- RenameIdents(cart_obj, new.cluster.ids)
cart_obj[["Treatment"]] <- Idents(cart_obj)

levels(factor(cart_obj@meta.data$orig.ident))
cart_obj[["Age"]] <- Idents(cart_obj)
new.cluster.ids <- c('Old','Young','Old','Young','Old','Young','Old','Young')
names(new.cluster.ids) <- levels(cart_obj)
cart_obj <- RenameIdents(cart_obj, new.cluster.ids)
cart_obj[["Age"]] <- Idents(cart_obj)

#Normalize and scale
DefaultAssay(cart_obj) <- "RNA"
all.genes <- rownames(cart_obj)
cart_obj <- NormalizeData(object = cart_obj, normalization.method = "LogNormalize", assay = "RNA")
cart_obj <- ScaleData(cart_obj, features = all.genes)
cart_obj <- RunPCA(cart_obj, features = all.genes)

#magic - imputation
cart_obj <- magic(cart_obj)

#Differential Expression analysis
Idents(cart_obj) <- cart_obj$Age
cart_o <- subset(cart_obj, idents = 'Old')
cart_y <- subset(cart_obj, idents = 'Young')

Idents(cart_y) <- cart_o$Treatment
cart_o <- PrepSCTFindMarkers(cart_o)
DE_treat_old <- FindMarkers(cart_o, ident.1 = "uPAR", ident.2 = "UT", test.use = "MAST", latent.vars ='Sex', logfc.threshold = .05, min.pct = .1, assay = 'RNA')
EnhancedVolcano(DE_treat_old, lab = rownames(DE_treat_old), x = 'avg_log2FC', y = 'p_val_adj', title = 'uPAR vs UT', pCutoff = 10e-1, FCcutoff = 0.25, xlim = c(-1.5,2.05), ylim = c(0,320), subtitle = 'All Cells' )
write.csv(DE_treat_old,'DE_Cart_old_upar_vs_ut.csv')

Idents(cart_y) <- cart_y$Treatment
cart_y <- PrepSCTFindMarkers(cart_y)
DE_treat_young <- FindMarkers(cart_y, ident.1 = "uPAR", ident.2 = "UT", test.use = "MAST", latent.vars ='Sex', logfc.threshold = .05, min.pct = .1, assay = 'RNA')
EnhancedVolcano(DE_treat_young, lab = rownames(DE_treat_young), x = 'avg_log2FC', y = 'p_val_adj', title = 'uPAR vs UT', pCutoff = 10e-1, FCcutoff = 0.25, xlim = c(-1,1.5), ylim = c(0,320), subtitle = 'All Cells' )
write.csv(DE_treat_young,'DE_Cart_young_upar_vs_ut.csv')


{
  cart_obj <- subset(cart_obj, idents = c('Enterocyte (Proximal)','Enterocyte (Distal)', 'Enteroendocrine', 'Goblet', 'Paneth'))  
  
  my_levels <- c('Enterocyte (Proximal)','Enterocyte (Distal)', 'Enteroendocrine', 'Goblet', 'Paneth')
  cart_obj$Cell_Type <- factor(x = cart_obj$Cell_Type, levels = my_levels)
  
  genes_stem_prolif <- unique(c('EGF',	'WNT9B',	'WNT3',	'BAMBI')) %>% convert_human_to_mouse_symbols()
  genes_antimicrobial <-  unique(c('DEFA5',	'DEFA6',	'DEFA27',	'DEFA30',	'DEFA26',	'DEFA20',	'MPTX1',	'MPTX2',	'LYZ1',	'LYZ2',	'KLF6',	'KLF4',	'KLF15',	'PLA2G2A',	'REG3B',	'REG3G',	'ITLN1',	'SPINK4',	'MT1',	'FABP1',	'FABP2',	'FABP4',	'MMP7',	'PLA2G12A', 'MUC2',	'CLCA1',	'PDIA4',	'PDIA5',	'ZG16',	'REG4',	'REG3B',	'ANG4', 'CASP6',	'CASP8',	'REG1',	'REG1',	'REG3A',	'REG3B',	'REG4'))%>% convert_human_to_mouse_symbols()         
  genes_hormone_secretion <-  unique(c('CHGA',	'GCG',	'GHRL',	'TRPA1',	'NTS',	'AFP',	'VEGFA',	'PCSK1',	'SCN3A',	'NEUROG3',	'NEUROD1',	'NKX2-2',	'FABP5',	'FABP1',	'FABP6',	'GIP',	'GLP1R',	'GLP2R',	'PYY',	'GAST',	'SCGN'))%>% convert_human_to_mouse_symbols()
  genes_nutrient_absorption <-  unique(c('SLC7A8',	'FABP1',	'CPS1',	'VDR',	'APOA4',	'RBP2',	'GSTA1',	'SLC5A1',	'SLC7A9',	'EHF',	'PPP1R14B',	'PP1R14D')) %>% convert_human_to_mouse_symbols()
  
  genes_stem_prolif <- genes_stem_prolif[!is.na(genes_stem_prolif)]
  genes_antimicrobial <- genes_antimicrobial[!is.na(genes_antimicrobial)]
  genes_hormone_secretion <- genes_hormone_secretion[!is.na(genes_hormone_secretion)]
  genes_nutrient_absorption <- genes_nutrient_absorption[!is.na(genes_nutrient_absorption)]
  
  All_Genes <- cart_obj@assays$RNA@data@Dimnames[[1]]
  genes_stem_prolif <- intersect(All_Genes, genes_stem_prolif)
  genes_antimicrobial <- intersect(All_Genes, genes_antimicrobial)   
  genes_hormone_secretion <- intersect(All_Genes, genes_hormone_secretion)  
  genes_nutrient_absorption <- intersect(All_Genes, genes_nutrient_absorption)
  
  cart_obj <- AddModuleScore(object = cart_obj, features = list(genes_stem_prolif), name = 'stem_prolif')
  cart_obj <- AddModuleScore(object = cart_obj, features = list(genes_antimicrobial), name = 'antimicrobial')
  cart_obj <- AddModuleScore(object = cart_obj, features = list(genes_hormone_secretion), name = 'hormone_secretion')
  cart_obj <- AddModuleScore(object = cart_obj, features = list(genes_nutrient_absorption), name = 'nutrient_absorption')
  
  Idents(cart_obj) <- cart_obj$Cell_Type
  Cell_Types <- levels(cart_obj$Cell_Type)
  stem_prolif_FC <- list()
  nutrient_absorption_FC <- list()
  antimicrobial_FC <- list()
  hormone_secretion_FC <- list()
  
  for(i in Cell_Types){
    stem_prolif_FC[[i]] <- FoldChange(cart_obj,ident.1 = "Old", ident.2 = "Young",features = genes_stem_prolif, group.by = 'Age', subset.ident = i)
  }
  stem_prolif_FC_aggr <- c()
  for(i in Cell_Types){
    stem_prolif_FC_aggr <- c(stem_prolif_FC_aggr, mean(as.numeric(stem_prolif_FC[[i]]$avg_log2FC)))
  }
  
  for(i in Cell_Types){
    nutrient_absorption_FC[[i]] <- FoldChange(cart_obj,ident.1 = "Old", ident.2 = "Young",features = genes_nutrient_absorption, group.by = 'Age', subset.ident = i)
  }
  nutrient_absorption_FC_aggr <- c()
  for(i in Cell_Types){
    nutrient_absorption_FC_aggr <- c(nutrient_absorption_FC_aggr, mean(as.numeric(nutrient_absorption_FC[[i]]$avg_log2FC)))
  }
  
  for(i in Cell_Types){
    antimicrobial_FC[[i]] <- FoldChange(cart_obj,ident.1 = "Old", ident.2 = "Young",features = genes_antimicrobial, group.by = 'Age', subset.ident = i)
  }
  antimicrobial_FC_aggr <- c()
  for(i in Cell_Types){
    antimicrobial_FC_aggr <- c(antimicrobial_FC_aggr, mean(as.numeric(antimicrobial_FC[[i]]$avg_log2FC)))
  }
  
  for(i in Cell_Types){
    hormone_secretion_FC[[i]] <- FoldChange(cart_obj,ident.1 = "Old", ident.2 = "Young",features = genes_hormone_secretion, group.by = 'Age', subset.ident = i)
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
  Idents(cart_obj) <- cart_obj$Cell_Type
  pvals <- c()
  for(j in scores){
    for(i in Cell_Types){
      selected_cells <- WhichCells(object = cart_obj ,idents = c(i))
      vln_data <- FetchData(cart_obj,
                            vars = c(j,"Age", "Cell_Type"),
                            cells = selected_cells,
                            slot = "data")
      vln_data <- melt(vln_data)
      vln_data <- vln_data[c('Age', 'Cell_Type', 'value')]
      
      pvals <- c(pvals, wilcox.test(vln_data[which(vln_data$Age == 'Old'),]$value, vln_data[which(vln_data$Age == 'Young'),]$value,  alternative = c("greater"))$p.value)
      print(i)
    }
  }
  pvals <- p.adjust(pvals, "BH", length(pvals))
  Module_Heatmap1$pval <- -log(pvals) 
  
  my_levels <- c('Enterocyte (Proximal)','Enterocyte (Distal)', 'Enteroendocrine', 'Goblet', 'Paneth')
  Module_Heatmap1$Cell_Type <- factor(Module_Heatmap1$Cell_Type, levels = my_levels)
  
  pathwayColorsDiff = rev(brewer.pal(9, "RdBu"))
  ggplot(Module_Heatmap1, aes(Cell_Type, variable, fill= value)) + geom_tile() + 
    scale_fill_gradientn(name = "Log2FC", colors = pathwayColorsDiff, limits = c(-1.11,1.11)) + coord_flip() +
    theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45), axis.text.y =   element_text( vjust = 0.1, hjust = 1), panel.background = element_blank(), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())
  
  ggsave(file="Control_CART_score_heatmap.pdf",  width=3, height=3, units="in")
  
  ggplot(Module_Heatmap1, aes(x = Cell_Type, y = variable, size = pval ,color= value)) + geom_point() + theme_vyom + 
    scale_color_gradientn(name = "Log2FC", colors = pathwayColorsDiff, limits = c(-1.11,1.11)) + coord_flip() + labs(size="-log(p)") +
    theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45), axis.text.y =   element_text( vjust = 0.1, hjust = 1), panel.background = element_blank(), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())
  
  ggsave(file="Control_CART_score_dotplot.pdf",  width=4, height=5, units="in")
}

sen_mayo <- c(
  "Acvr1b", "Ang", "Angpt1", "Angptl4", "Areg", "Axl", "Bex3", "Bmp2", "Bmp6",
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
Inflammatory_score_new <- c('Ccl4','Ccl5','Ccl9','Ccl20','Cd74','Cox4I1','Cox8A','Ctsb','Ctsd','H2-Ac','H2-Ab1','H2-Q7','Ifi47','Ifit1Bl1','Ifit2','Ifit3','Ifitm3','Irf7','Isg15','Stat2','Stat3','Tap1','Ucp2')

All_Genes <- cart_old@assays$RNA@data@Dimnames[[1]]
sen_mayo <- intersect(All_Genes, sen_mayo)
Inflammatory_score_new <- intersect(All_Genes, Inflammatory_score_new)


score.list <- list(sen_mayo, Inflammatory_score_new)
for(i in  names(score.list)){
  gene.list <- as.vector(score.list[[i]])
  mean.exp <- zscore(colMeans(x = cart_old@assays$RNA$data[gene.list, ], na.rm = TRUE), dist = 'norm')
  if (all(names(x = mean.exp) == rownames(x = cart_old@meta.data))) {
    cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
        "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
    cart_old@meta.data[[i]] <- mean.exp
  }
}


#sen score violin plots
Sen_score_list$SASP
library(readxl)
Sen_score_list <- read_excel("analysis/Senescence_and_SASP_genesets_2022_mouse_Corina.xlsx")

score_list <- colnames(Sen_score_list) 
for (i in score_list){

  selected_cells <- names(cart_old$Cell_Type)
  vln_data <- FetchData(cart_old,
                        vars = c(i,"Treatment", "Cell_Type"),
                        cells = selected_cells,
                        slot = "data")
  vln_data <- melt(vln_data)
  vln_data <- vln_data[c('Treatment', 'Cell_Type', 'value')]
  
  dodge <- position_dodge(width = .9)
  All <- VlnPlot(cart_old,group.by = 'Cell_Type' , split.by = "Treatment", features = i, pt.size = 0, assay = "MAGIC_RNA", cols = c('#7CA1CC' ,'#FF4902'), log = FALSE, split.plot = TRUE) + 
    theme(legend.position = 'none', axis.line = element_line(size = .3),axis.title.x = element_blank(), axis.ticks = element_line(size = .3), text = element_text(size=7), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6)) +
    geom_boxplot(width=0.3,outlier.shape = NA, coef = 0, lwd=.2, position = dodge) +
    stat_compare_means(data= vln_data, aes(x = Cell_Type, y = value, fill = Treatment), hide.ns = TRUE, method = "t.test",  label = "p.signif", size = 1)  
  All$layers[[1]]$aes_params$size = .15
  ggsave(file = paste0('CART_INTEST_OLD_scrna_vln_', i, '.pdf'), plot=All, width=3, height=2.5, units="in")
}


#pseudotime analysis
library(SeuratWrappers)
library(monocle3)
library(Signac)
library(org.Mm.eg.db)
gene_symbol <- as.list(org.Mm.egSYMBOL)
Idents(cart_o) <- cart_o$Cell_Type

cds <- as.cell_data_set(cart_o)
cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(cart_o[["RNA"]])
cds <- cluster_cells(cds = cds, reduction_method = "UMAP")
cds <- learn_graph(cds, use_partition = FALSE,   close_loop = FALSE)
cds <- order_cells(cds, reduction_method = "UMAP")


plot_cells(cds = cds, color_cells_by = "pseudotime", show_trajectory_graph = TRUE, cell_size = 2)
pseudo_umap <- plot_cells(cds, color_cells_by = "pseudotime", show_trajectory_graph = FALSE, label_cell_groups=FALSE, label_leaves=FALSE, label_branch_points=FALSE, trajectory_graph_color = "#A8A8A8", graph_label_size=1.5,  cell_size = 2) 
pseudo_umap + scale_fill_gradientn(colors = annColors$Pseudotime)  + geom_polygon(data = hulls2, aes(x = x, y = y, group=CellType), fill=NA, color="black", size=0.3, alpha = 0.3) + geom_shadowtext(data = clusterMedian, aes(x = UMAP_1,y= UMAP_2, label = CellType, group = as.factor(CellType)), size = 3.5, bg.colour="black") + labs(fill = "Pseudotime")
ggsave(file="cart_o_Pseudotime_Umap_no_trajectory.svg", width=10, height=10)
plot_cells(cds,
           color_cells_by = "Cell_Type",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

cart_o <- AddMetaData(
  object = cart_o,
  metadata = cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "Pseudotime"
)


Idents(cart_o) <- cart_o$Cell_Type
FeaturePlot(cart_o, c("Pseudotime"), pt.size = .0001, label = TRUE, repel = TRUE) +scale_fill_gradientn(colors = brewer.pal(11,'Spectral')) + scale_color_gradientn(colors =  brewer.pal(11,'Spectral'))
ggsave(file = 'cart_o_pseudotime_umap.pdf', width=5, height=5, units="in")

#pseudotime_gene
uPAR_rownames <- rownames(cart_o@meta.data[cart_o@meta.data$Treatment == 'uPAR',] )
Control_rownames <- rownames(cart_o@meta.data[cart_o@meta.data$Treatment == 'UT',] )
Cell_type_frame <- data.frame(cart_o$Treatment,cart_o$orig.ident )
all <-  cbind(type = Cell_type_frame, Pseudotime = cart_o$Pseudotime)
colnames(all) <- c("Treatment",'Replicate',"Pseudotime")

all <- all %>% mutate(quantile = ntile(all$Pseudotime, 10))
quart <- data.frame(var1 = tapply(all$Pseudotime, all$quantile, mean), var2 = tapply(all$Pseudotime, all$quantile, max))

gene.list <- c('H2-D1','Gzma','Gzmb','Cd74', "H2-Ab1", "H2-Aa",'Lgr5','Lgr4','Ascl2', 'Olfm4','Sox9','Malat1','Ccnd1')

for(i in gene.list){
  colnames(cart_o)
  DefaultAssay(cart_o)
  
  score_iter = i
  quart$maxi = max(cart_o@assays$MAGIC_RNA@data[i, ])
  uPAR <- subset(all, Treatment=='uPAR')
  UT <- subset(all, Treatment=='UT')
  
  uPAR <-  cbind(uPAR, cart_o[,uPAR_rownames]@assays$MAGIC_RNA@data[i, ])
  UT <-  cbind(UT, cart_o[,Control_rownames]@assays$MAGIC_RNA@data[i, ])
  
  colnames(uPAR) <- c("Sample",'Replicate',"Pseudotime",'quantile' ,'score')
  colnames(UT) <- c("Sample",'Replicate',"Pseudotime",'quantile' ,'score')
  
  p.vals.prop = ''
  for(j in 1:10){
    uPARinqunat  <- subset(uPAR, quantile== j)
    Cinqunat  <- subset(UT, quantile== j)
    
    p.vals.prop = c(p.vals.prop, wilcox.test(uPARinqunat$score, Cinqunat$score, alternative = "two.sided")$p.value)
  }
  p.vals.prop = p.vals.prop[2:11]
  p.vals.prop <- as.numeric(p.vals.prop)
  p.vals.prop[p.vals.prop < 0.05] =  "*"
  p.vals.prop[p.vals.prop > 0.05] =  ""
  
  
  
  quart$p_asterik = p.vals.prop
  
  
  Pseudotimedot <- ggplot() +
    # Points
    #ggrastr::rasterise(geom_jitter(data=UT, aes(x=Pseudotime, y=score), size=0, alpha=0.05, width=0,stroke = 0, colour = "#7CA1CC" ), dpi = 400) + 
    #ggrastr::rasterise(geom_jitter(data=uPAR, aes(x=Pseudotime, y=score), size=0, alpha=0.05, width=0,stroke = 0, colour = "#FF4902"), dpi = 400) +
    # lines
    geom_smooth(data=UT, aes(x=Pseudotime, y=score), fill="black", colour="#7CA1CC", size=.5) + 
    geom_smooth(data=uPAR, aes(x=Pseudotime, y=score), fill="black", colour="#FF4902", size=.5) +
    theme( panel.background = element_rect(fill = "white", colour = "Black")) + ylab(score_iter) + theme_vyom +
    geom_vline(data = quart, aes(xintercept = var1), size = .25) + geom_text(data = quart, aes(var2 , maxi, label = p_asterik, fontface="bold"),size = 1.75 , position=position_dodge2(0.5), vjust=.75, inherit.aes = FALSE) + 
    theme(text = element_text(size=1),axis.text = element_text(size = 2) ,plot.title = element_text(size = 6), axis.text.x = element_text(size = 4), axis.text.y = element_text(size = 4), axis.title.x = element_text(size = 6), axis.title.y = element_text(size = 6)) +
    scale_x_continuous(expand = c(0, 0))
  Pseudotimedot
  ggsave(file= paste0(i,"_uPAR_scrna_Pseudotime_expression.pdf"),plot = Pseudotimedot, width=3, height=2)
}

#pseudotime density plots
{
  dimredP = data.frame(cart_o@reductions$umap@cell.embeddings)
  dimredP$Cell = rownames(dimredP)
  all.equal(names(cart_o$seurat_clusters), dimredP$Cell)
  all.equal(names(cart_o$Treatment), dimredP$Cell)
  
  dimredP$Cluster = cart_o$Cell_Type # Overwrite clusters with cell types for current plots
  dimredP$Type = cart_o$Treatment
  dimredP$CellType = dimredP$Cluster
  dimredP$Pseudotime = cart_o$Pseudotime
  
  
  i = 'Stem'
  cell = i
  dimredP <- dimredP %>% mutate(quantile = ntile(Pseudotime, 10))
  dimredaub = dimredP %>% filter(CellType == i)
  quart <- data.frame(var1 = tapply(dimredP$Pseudotime, dimredP$quantile, mean), var2 = tapply(dimredP$Pseudotime, dimredP$quantile, max))
  p.vals.prop = ''
  dimredc <- dimredaub %>% filter(Type == 'UT')
  dimreda <- dimredaub %>% filter(Type == 'uPAR')
  
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
  ggsave(paste0("CART_UPAR_vs_UT_old_Pseduodensity_stem.pdf"), plot = Density_diff_stem_1, width = 3, height = 2, units="in")	
  
  #all cells
  
  dimredaub = dimredP
  quart <- data.frame(var1 = tapply(dimredP$Pseudotime, dimredP$quantile, mean), var2 = tapply(dimredP$Pseudotime, dimredP$quantile, max))
  p.vals.prop = ''
  dimredc <- dimredaub %>% filter(Type == 'UT')
  dimreda <- dimredaub %>% filter(Type == 'uPAR')
  
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
  ggsave(paste0("CART_UPAR_vs_UT_old_Pseduodensity_All.pdf"), plot = Density_diff_all, width = 3, height = 2, units="in")	
  
  Density_all <-ggplot(dimredP, aes(x=Pseudotime, fill=Type)) + geom_density(alpha=0.6, size = .1) + theme_vyom + scale_fill_manual(values=unlist(colorsType)) + theme(legend.position = "none") + ylab('Density') + xlab('Pseudotime') + geom_vline(data = quart, aes(xintercept = var2), size = .075) + geom_text(data = quart, aes(var1 , maxi2, label = p_asterik), position=position_dodge2(0.2), size = 1.5, vjust=.35, inherit.aes = FALSE) + theme(text = element_text(size=1),axis.text = element_text(size = 2) ,plot.title = element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(size = 6), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.ticks.x = element_blank(), plot.margin = unit(c(.1, .1, 0, 0), "cm"))
  Density_all <-ggplot(dimredP, aes(x=Pseudotime, fill=Type)) + geom_density(alpha=0.6, size = .1) + 
    theme_vyom + scale_fill_manual(values=unlist(colorsType)) + scale_y_continuous(expand = c(0.2,0), breaks = scales::breaks_extended(n = 5)) +
    theme(legend.position = "none") + ylab('Density') + xlab('Pseudotime') + geom_vline(data = quart, aes(xintercept = var2), size = .075) + geom_text(data = quart, aes(var1 , maxi2, label = p_asterik), position=position_dodge2(0.2), size = 1.5, vjust=.35, inherit.aes = FALSE) + theme(text = element_text(size=1),axis.text = element_text(size = 2) ,plot.title = element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(size = 6), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.ticks.x = element_blank(), plot.margin = unit(c(.2, .1, 0, .2), "cm"))
  
  
  plot_grid(Density_all, Density_diff_all, Density_diff_stem_1, Density_diff_stem_2, labels = c('','All', 'Stem 1', 'Stem 2'), label_size = 6, ncol = 1, align = "v", hjust = -.005) #'All Cells', 'Stem 1', 'Stem 2'
  ggsave(paste0("AA_pseudo_density.pdf"), width = 2.2, height = 3.1, units="in")	
  
  plot_grid(Density_all, Density_all_stem_1, Density_all_stem_2, labels = c('All', 'Stem 1', 'Stem 2'), label_size = 6, ncol = 1, align = "v", hjust = -.02) #'All Cells', 'Stem 1', 'Stem 2'
  ggsave(paste0("AA_pseudo_density_all.pdf"), width = 2.2, height = 3.1, units="in")	
}

#housekeeping
# saveRDS(cart_obj, file = "./data/Seurat_Objects/CART_Sobj.rds")
# cart_obj <- readRDS("./data/Seurat_Objects/CART_Sobj.rds", refhook = NULL)
cart_obj <- readRDS("./CSHL\ Dropbox\ Team\ Dropbox/Vyom\ Shah/Vyom\ Shah/Collabs/Ramesh\ EPCAM\ CD3/epcam_cd3_ALL.rds", refhook = NULL)
options(future.globals.maxSize = 4000 * 1024^10)
write.csv(cart_obj@meta.data, 'Treated_CAR_T_metadata.csv')

DimPlot(cart_obj)

