#UPAR+ and UPAR- scRNA seq analysis
# Vyom Shah 

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

theme_linedraw2 = theme_linedraw() + theme(strip.background=element_rect(fill="grey80", colour="grey50", size=0.2), strip.text.x=element_text(colour="black"), strip.text.y=element_text(colour="black"))
theme_vyom = theme_linedraw2 + theme(legend.position="right", legend.title=element_text(size=15), legend.text=element_text(size=14), axis.text.x = element_text(size=12, angle=-90, hjust=0, vjust=0.5), axis.text.y=element_text(size=12), axis.title=element_text(size=15), axis.title.y=element_text(vjust=1), plot.title = element_text(size=18, vjust=1.5), strip.background = element_rect(fill="#EEEEEE"), strip.text = element_text(size = 11), panel.grid.major = element_line(colour = "grey98"), panel.grid.minor = element_blank())

Vyom_color_Umap <- CustomPalette(low = "purple", high = "yellow", mid = "black", k = 100)

F_Neg <- Read10X(data.dir = "/Users/vyom/data/upar_scrna/count/Amor_CAV12_F_neg/outs/filtered_feature_bc_matrix/")
M_Neg <- Read10X(data.dir = "/Users/vyom/data/upar_scrna/count/Amor_CAV12_M_neg/outs/filtered_feature_bc_matrix/")
F_Plus <- Read10X(data.dir = "/Users/vyom/data/upar_scrna/count/Amor_CAV12_F_plus/outs/filtered_feature_bc_matrix/")
M_Plus <- Read10X(data.dir = "/Users/vyom/data/upar_scrna/count/Amor_CAV12_M_plus/outs/filtered_feature_bc_matrix/")

F_Neg <- CreateSeuratObject(F_Neg, project = "Female_Neg")
M_Neg <- CreateSeuratObject(M_Neg, project = "Male_Neg")
F_Plus <- CreateSeuratObject(F_Plus, project = "Female_Plus")
M_Plus <- CreateSeuratObject(M_Plus, project = "Male_Plus")

d <- merge(F_Neg, y = c(M_Neg, F_Plus,M_Plus), add.cell.ids = c("Female_Neg", "Male_Neg","Female_Plus", "Male_Plus"), project = "UPAR_sort")
d
d[["percent.mt"]] <- PercentageFeatureSet(d, pattern = "mt-")
VlnPlot(d, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),pt.size = 0 ,ncol = 3)

plot1 <- FeatureScatter(d, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(d, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
d <- subset(d, subset = percent.mt < 5)
d <- subset(d, subset = nCount_RNA > 250 & nCount_RNA < 30000)
d <- subset(d, subset = nFeature_RNA > 2000 )


levels(factor(d@meta.data$orig.ident))
Idents(d)<- d$orig.ident
new.cluster.ids <- c("UPAR-", "UPAR-", "UPAR+",  "UPAR+")     
d[["Treatment"]] <- Idents(d)
names(new.cluster.ids) <- levels(d)
d <- RenameIdents(d, new.cluster.ids)
d[["Treatment"]] <- Idents(d)

levels(factor(d@meta.data$orig.ident))
Idents(d)<- d$orig.ident
new.cluster.ids <- c("Female", "Male", "Female",  "Male")     
d[["Sex"]] <- Idents(d)
names(new.cluster.ids) <- levels(d)
d <- RenameIdents(d, new.cluster.ids)
d[["Sex"]] <- Idents(d)

Data.list <- SplitObject(d, split.by = "Sex")
Data.list <- Data.list[c("Female", "Male")]
for (i in 1:length(Data.list)) {
  
  Data.list[[i]] <- SCTransform(Data.list[[i]], verbose = FALSE)
}

# Normilization
#select highly variable genes 
Data.features <- SelectIntegrationFeatures(object.list = Data.list, nfeatures = 3000)
Data.list <- PrepSCTIntegration(object.list = Data.list, anchor.features = Data.features, 
                                verbose = FALSE)
Data.anchors <- FindIntegrationAnchors(object.list = Data.list, normalization.method = "SCT", 
                                       anchor.features = Data.features, verbose = FALSE)
UPAR_obj <- IntegrateData(anchorset = Data.anchors, normalization.method = "SCT", 
                          verbose = TRUE)

# Visulization and Clustering
UPAR_obj <- RunPCA(UPAR_obj)

UPAR_obj <- RunUMAP(UPAR_obj, dims = 1:50, n.neighbors = 5, n.epochs = 500)
DimPlot(UPAR_obj, reduction = "umap",group.by = 'orig.ident', label = TRUE)

UPAR_obj <- FindNeighbors(UPAR_obj, dims = 1:50, k.param = 30, compute.SNN = TRUE)
UPAR_obj <- FindClusters(UPAR_obj, resolution = 1)
DimPlot(UPAR_obj, reduction = "umap", label = TRUE)

meta_check <- UPAR_obj@meta.data

DimPlot(UPAR_obj, reduction = "umap", label = TRUE, group.by = 'Cell_Type')
prop.table(x = table(UPAR_obj$integrated_snn_res.1, UPAR_obj$Treatment), margin = 2)
DotPlot_Sig <- c("Lgr5","Ascl2","Olfm4","Gkn3","S100a6","Ly6a","Anxa3", "Areg","Tubb5","Syce2","Stmn1","Fbxo5",'Cenpa','Ccna2','Ube2c','Cdkn3',"Apoa1","Apoa4","Fabp1","Adh6a",'Tmigd1', 'Fabp6', 'Slc51b', 'Slc51a', "Chgb","Tac1","Tph1","Neurog3", "Muc2","Fcgbp","Atoh1","Agr2","Lyz1","Defa17","Defa24","Ang4","Pou2f3","Avil","Tuba1a","Adh1",'Ptprc','Cd8a','Cd4','Ighm') 
DotPlot_Sig <- unique(c('Cd3g','Cd3e','Cd8a','Cd4', 'Lag3','Pdcd1','Havcr2','Tox','Foxp3','Il10','Tcrg-C1','Gata3','Gzmb','Tbx21','Ifng','Il2rb','Klrd1','Cd19','Ighm','Ighg1','Cd74','Ciita','Nrc1','Klre1','Itgam','Itgax','H2-Eb1','H2-Ab1','Arg1','Mrc1','Tgfbi','Ccr2','Vegfa','Prdx1','Clec4d','Ccl5','Cd83','Ccr7','Fcn1','Msrb1','Ly6g','Col3a1','Sparc'))

DotPlot(UPAR_obj, features = DotPlot_Sig, dot.scale = 10, scale.max= 100, assay = 'SCT') + labs(y= "Cell type", x="") + scale_colour_distiller( palette ="RdYlBu") + scale_size(range = c(0, 1.3)) +
  theme(legend.key.size = unit(.1, "in"),axis.line = element_line(size = .3), axis.ticks = element_line(size = .3), axis.text.x = element_text( angle = 90, hjust = 1, vjust= .01),  axis.title.x = element_blank(), axis.title.y = element_blank())
ggsave(file="UPAR_Cluster_check.pdf", width=8, height=6)

# Create an empty vector of length 20
r_vector <- character(20)

# Assign the specified strings to the specified positions (adding 1 to each position)
r_vector[c(1)] <- "Stem"
r_vector[c(5, 6, 7)] <- "Transit Amplifying"
r_vector[c(3, 10, 13)] <- "Enterocyte Progenitor"
r_vector[c(2, 8, 9, 11)] <- "Enterocyte"
r_vector[c(15, 19)] <- "Enteroendocrine"
r_vector[c(4, 12)] <- "Goblet"
r_vector[c(20)] <- "Paneth"
r_vector[14] <- "Tuft"
r_vector[c(17, 18)] <- "T Cell"
r_vector[16] <- "Macrophage"

UPAR_obj[["Cell_Type"]] <- Idents(UPAR_obj)
names(r_vector) <- levels(UPAR_obj)
UPAR_obj <- RenameIdents(UPAR_obj, r_vector)
UPAR_obj[["Cell_Type"]] <- Idents(UPAR_obj)
DimPlot(UPAR_obj, reduction = "umap", group.by= 'Cell_Type')


as.factor(UPAR_obj$Treatment)
my_levels <- c('Stem', 'Transit Amplifying', 'Enterocyte Progenitor','Enterocyte', 'Enteroendocrine', 'Goblet', 'Paneth','Tuft', 'T Cell', 'Macrophage')
my_levels1 <- c("UPAR-", "UPAR+")
UPAR_obj$Cell_Type <- factor(x = UPAR_obj$Cell_Type, levels = my_levels)

DimPlot(UPAR_obj, reduction = "umap", group.by= 'Cell_Type', pt.size = .0001)
ggsave(file = paste0('UPAR_scrna_UMAP.pdf'), width=6, height=4, units="in")

DotPlot_Sig <- c("Lgr5","Ascl2","Olfm4","Gkn3","Tubb5","Syce2","Stmn1","Fbxo5",'Cenpa','Ccna2','Ube2c','Cdkn3',"Apoa1","Apoa4","Fabp1","Adh6a",'Tmigd1', 'Fabp6', 'Slc51b', 'Slc51a', "Chgb","Tac1","Tph1","Neurog3", "Muc2","Fcgbp","Atoh1","Agr2","Lyz1","Defa17","Defa24","Ang4","Pou2f3","Avil","Tuba1a","Adh1",'Ptprc','Cd8a','Cd4','Ighm') 
DotPlot(UPAR_obj, group.by = 'Cell_Type', features = DotPlot_Sig, dot.scale = 10, scale.max= 100, assay = 'SCT') + labs(y= "Cell type", x="") + scale_colour_distiller( palette ="RdYlBu") + scale_size(range = c(0, 1.3)) +
  theme(legend.key.size = unit(.1, "in"),axis.line = element_line(size = .3), axis.ticks = element_line(size = .3), axis.text.x = element_text( angle = 90, hjust = 1, vjust= .01),  axis.title.x = element_blank(), axis.title.y = element_blank())
DotPlot(UPAR_obj, group.by = 'Cell_Type', features = DotPlot_Sig, dot.scale = 10, scale.max= 100, assay = 'SCT') + labs(y= "Cell type", x="") + scale_colour_distiller( palette ="RdYlBu") + scale_size(range = c(0, 1.4)) +
  theme(axis.line = element_line(size = .3), axis.ticks = element_line(size = .3), text = element_text(size=4), axis.text.x = element_text(size = 4, angle = 90, hjust = 1, vjust= .01), axis.text.y = element_text(size = 4), axis.title.x = element_blank(), axis.title.y = element_blank())

ggsave(file="UPAR_Cluster_check.pdf", width=4, height=2)

Prop_table<- prop.table(x = table(UPAR_obj$Cell_Type, UPAR_obj$orig.ident), margin = 2)
Prop_Table <- as.data.frame(Prop_table, row.names = NULL, optional = FALSE,
                            make.names = TRUE, stringsAsFactors = default.stringsAsFactors())
levels(as.factor(UPAR_obj$orig.ident))
Prop_Table1 <- Prop_Table%>%mutate(Var2=recode(Var2,Female_Neg="UPAR-",Female_Plus="UPAR+", Male_Neg="UPAR-", Male_Plus="UPAR+"))

celltype_sample_norm3 = summarySE(Prop_Table1, measurevar="Freq", groupvars=c("Var1","Var2"))

my_levels <- c('Stem', 'Transit Amplifying', 'Enterocyte Progenitor','Enterocyte', 'Enteroendocrine', 'Goblet', 'Paneth','Tuft', 'T Cell', 'Macrophage')
celltype_sample_norm3$Var1 = factor(celltype_sample_norm3$Var1, levels = my_levels)
celltype_sample_norm3$Var2 = factor(celltype_sample_norm3$Var2, levels = c('UPAR-','UPAR+'))

ggplot(celltype_sample_norm3, aes(x = Var1, y = Freq, fill = Var2)) + geom_bar(stat="identity", position="dodge") + 
  geom_errorbar(aes(ymin=Freq-se, ymax=Freq+se),position=position_dodge(width = 0.85),width=0.3, size=0.25) + theme_vyom + 
  scale_fill_manual(values = colorsType) + expand_limits(y = c(0)) + 
  theme( axis.text.x = element_text(angle = 45, size =  6, vjust = 1, hjust = 1), panel.background = element_rect(fill = "white", colour = "Black"), axis.text.y = element_text(size = 6), axis.title.x = element_text(size = 7), axis.title.y = element_text(size = 7)) + 
  xlab("Cell Type") + ylab("Fraction of Cells") + labs(fill = "Age") + scale_y_continuous(expand = expansion(mult = c(0, .1)))
ggsave( "UPAR_Sort_prop_errbar.pdf", width=3.75, height=2, units="in")

#Normalize and scale
DefaultAssay(cart_obj) <- "RNA"
all.genes <- rownames(cart_obj)
cart_obj <- NormalizeData(object = cart_obj, normalization.method = "LogNormalize", assay = "RNA")
cart_obj <- ScaleData(cart_obj, features = all.genes)
cart_obj <- RunPCA(cart_obj, features = all.genes)

#run magic
DefaultAssay(UPAR_obj) <- "RNA"
UPAR_obj <- Rmagic::magic(UPAR_obj)

#differential expression:
Idents(UPAR_obj) <- UPAR_obj$Treatment
UPAR_obj <- PrepSCTFindMarkers(UPAR_obj)
DE_upar_sort <- FindMarkers(UPAR_obj, ident.1 = "UPAR+", ident.2 = "UPAR-", test.use = "MAST", latent.vars ='Sex', logfc.threshold = .05, min.pct = .1, assay = 'SCT')
EnhancedVolcano(DE_upar_sort, lab = rownames(DE_upar_sort), x = 'avg_log2FC', y = 'p_val_adj', title = 'uPAR+ vs uPAR-', pCutoff = 10e-1, FCcutoff = 0.25, xlim = c(-2.5,2), ylim = c(0,320), subtitle = 'All Cells' )
write.csv(DE_upar_sort,'DE_UPAR_Sort.csv')

#summarize clustering for some plotting functionality:
Idents(UPAR_obj) <- UPAR_obj$Cell_Type
current.cluster.ids <- c('Stem', 'Transit Amplifying', 'Enterocyte Progenitor','Enterocyte', 'Enteroendocrine', 'Goblet', 'Paneth','Tuft', 'T Cell', 'Macrophage')
new.cluster.ids <- c('Stem', 'Stem', 'Epithelial','Epithelial', 'Epithelial', 'Epithelial', 'Epithelial','Epithelial', 'Immune', 'Immune')
UPAR_obj@active.ident <- plyr::mapvalues(x = UPAR_obj@active.ident, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(UPAR_obj, reduction = 'umap')
UPAR_obj$Cell_Type_condense <- Idents(UPAR_obj)
DimPlot(UPAR_obj, group.by = "Cell_Type_condense",  reduction = 'umap')
DimPlot(UPAR_obj, group.by = "Cell_Type",  reduction = 'umap')

#load in scores of relevance:
senescence <- c('Psap','Timp2','Itm2b','Lgmn','Igf1','Arl6ip1','Ckb','Lsmem1','Igfbp4','Gsn','Zbtb20','Ccl2','Ccl7','Trio','Tulp4','Cxcl1')
stem_score <- c('Lgr4','Myc', 'Sox9', 'Olfm4', 'Malat1', 'Hopx', 'Ccnd1' )
MHC_score <- c('H2-D1','H2-K','H2-Q1','H2-Q2','H2-T23','H2-T3','Cd74','H2-Ac','H2-Ab1','H2-Q7')
Hormone <- c('Fabp1','Gast','Scgn', 'Pcsk1','Glp1r','Vegfa')
Absorption <- c('Cps1', 'Fabp1', 'Rbp2','Slc7a8', 'Vdr', 'Ppp1r14b')
antimicrobial <- c("Casp8","Itln1","Fabp2","Pdia4","Zg16","Muc2","Mmp7")
genes_stem_prolif <- unique(c('EGF',	'WNT9B',	'WNT3',	'BAMBI')) %>% convert_human_to_mouse_symbols()
genes_antimicrobial <-  unique(c('DEFA5',	'DEFA6',	'DEFA27',	'DEFA30',	'DEFA26',	'DEFA20',	'MPTX1',	'MPTX2',	'LYZ1',	'LYZ2',	'KLF6',	'KLF4',	'KLF15',	'PLA2G2A',	'REG3B',	'REG3G',	'ITLN1',	'SPINK4',	'MT1',	'FABP1',	'FABP2',	'FABP4',	'MMP7',	'PLA2G12A', 'MUC2',	'CLCA1',	'PDIA4',	'PDIA5',	'ZG16',	'REG4',	'REG3B',	'ANG4', 'CASP6',	'CASP8',	'REG1',	'REG1',	'REG3A',	'REG3B',	'REG4'))%>% convert_human_to_mouse_symbols()         
genes_hormone_secretion <-  unique(c('CHGA',	'GCG',	'GHRL',	'TRPA1',	'NTS',	'AFP',	'VEGFA',	'PCSK1',	'SCN3A',	'NEUROG3',	'NEUROD1',	'NKX2-2',	'FABP5',	'FABP1',	'FABP6',	'GIP',	'GLP1R',	'GLP2R',	'PYY',	'GAST',	'SCGN'))%>% convert_human_to_mouse_symbols()
genes_nutrient_absorption <-  unique(c('SLC7A8',	'FABP1',	'CPS1',	'VDR',	'APOA4',	'RBP2',	'GSTA1',	'SLC5A1',	'SLC7A9',	'EHF',	'PPP1R14B',	'PP1R14D')) %>% convert_human_to_mouse_symbols()
Inflammatory_score_new <- c('Ccl4','Ccl5','Ccl9','Ccl20','Cd74','Cox4I1','Cox8A','Ctsb','Ctsd','H2-Ac','H2-Ab1','H2-Q7','Ifi47','Ifit1Bl1','Ifit2','Ifit3','Ifitm3','Irf7','Isg15','Stat2','Stat3','Tap1','Ucp2')
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


genes_stem_prolif <- genes_stem_prolif[!is.na(genes_stem_prolif)]
genes_antimicrobial <- genes_antimicrobial[!is.na(genes_antimicrobial)]
genes_hormone_secretion <- genes_hormone_secretion[!is.na(genes_hormone_secretion)]
genes_nutrient_absorption <- genes_nutrient_absorption[!is.na(genes_nutrient_absorption)]
sen_BRD4 <- sen_BRD4[!is.na(sen_BRD4)]

All_Genes <- UPAR_obj@assays$RNA@data@Dimnames[[1]]
stem_score <- intersect(All_Genes, stem_score)
MHC_score <- intersect(All_Genes, MHC_score)
Hormone <- intersect(All_Genes, Hormone)
Absorption <- intersect(All_Genes, Absorption)
antimicrobial <- intersect(All_Genes, antimicrobial)
genes_stem_prolif <- intersect(All_Genes, genes_stem_prolif)
genes_antimicrobial <- intersect(All_Genes, genes_antimicrobial)
genes_hormone_secretion <- intersect(All_Genes, genes_hormone_secretion)
genes_nutrient_absorption <- intersect(All_Genes, genes_nutrient_absorption)
Inflammatory_score_new <- intersect(All_Genes, Inflammatory_score_new)
sen_mayo <- intersect(All_Genes, sen_mayo)

score.list <- list(sen_mayo, Inflammatory_score_new)
for(i in  names(score.list)){
  gene.list <- as.vector(score.list[[i]])
  mean.exp <- zscore(colMeans(x = UPAR_obj@assays$RNA$data[gene.list, ], na.rm = TRUE), dist = 'norm')
  if (all(names(x = mean.exp) == rownames(x = UPAR_obj@meta.data))) {
    cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
        "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
    UPAR_obj@meta.data[[i]] <- mean.exp
  }
}


UPAR_obj$Treatment
score_lsit <- c("senescence", "stem_score", "MHC_score", "Hormone" , "Absorption", "antimicrobial", "Inflammatory_score_new", "sen_mayo", "sen_BRD4" )
for(i in score_lsit) {
  selected_cells <- names(UPAR_obj$Cell_Type)
  # selected_cells <- rownames(UPAR_obj@meta.data[UPAR_obj$Treatment == c('Control','Arasco'),] )
  
  vln_data <- FetchData(UPAR_obj,
                        vars = c(i,"Treatment", "Cell_Type"),
                        cells = selected_cells,
                        slot = "data")
  vln_data <- melt(vln_data)
  vln_data <- vln_data[c('Treatment', 'Cell_Type', 'value')]
  
  All <- VlnPlot(UPAR_obj,group.by = 'Cell_Type' , split.by = "Treatment", features = i, pt.size = 0, assay = "MAGIC_RNA", cols = c('#7CA1CC' ,'#FF4902'), log = FALSE, split.plot = TRUE) + 
    theme(legend.position = 'none')  +
    theme(axis.line = element_line(size = .3),axis.title.x = element_blank(), axis.ticks = element_line(size = .3), text = element_text(size=7), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6)) #+ scale_y_continuous(expand = expansion(mult = c(0, .1)))
  dodge <- position_dodge(width = .9)
  All <- All + geom_boxplot(width=0.3,outlier.shape = NA, coef = 0, lwd=.2, position = dodge) +
    stat_compare_means(data= vln_data, aes(x = Cell_Type, y = value, fill = Treatment), hide.ns = TRUE, method = "t.test",  label = "p.signif", size = 1)  
  All$layers[[1]]$aes_params$size = .15
  ggsave(file = paste0('UPAR_scrna_vln_', i, '.pdf'), plot=All, width=3, height=2.5, units="in")
}

source("~/analysis/split_violin.R")
gene.list <- c('Lgr4', 'Lgr5', 'Sox9', 'Olfm4', 'Malat1', 'Hopx', 'Ccnd1' )

DefaultAssay(UPAR_obj) <- "MAGIC_RNA"
selected_cells <- names(UPAR_obj$Cell_Type[UPAR_obj$Cell_Type == c("Stem")])
vln_data <- FetchData(UPAR_obj,
                      vars = c(gene.list,"Treatment"),
                      cells = selected_cells,
                      slot = "data")

vln_data <- melt(vln_data)

All <- ggplot(vln_data, aes(x = variable, y = value, fill= Treatment)) + geom_split_violin(scale="width", trim = TRUE, size = .1) + theme_vyom + theme(legend.position = 'none') + 
  geom_boxplot(width=0.2,outlier.shape = NA, coef = 0, lwd=.2) + xlab('') + ylab('Log2 Expression Level')+ scale_fill_manual(values= c('#7CA1CC' ,'#FF4902')) +
  theme(text = element_text(size=7), axis.text.x = element_text(size = 6, angle = 45, hjust = 1, vjust = 1), axis.text.y = element_text(size = 6), axis.title.y = element_text(size = 7)) +
  stat_compare_means( method = "wilcox.test",  label = "p.signif", size = 2) + scale_y_continuous(trans='log2', expand = expansion(mult = c(0, 0.1)))
All
ggsave(file = 'UPAR_obj_Stem_score.pdf', plot=All, width=2.85, height=2, units="in")


#functional score plot:
{
  genes_stem_prolif <- unique(c('EGF',	'WNT9B',	'WNT3',	'BAMBI')) %>% convert_human_to_mouse_symbols()
  genes_antimicrobial <-  unique(c('DEFA5',	'DEFA6',	'DEFA27',	'DEFA30',	'DEFA26',	'DEFA20',	'MPTX1',	'MPTX2',	'LYZ1',	'LYZ2',	'KLF6',	'KLF4',	'KLF15',	'PLA2G2A',	'REG3B',	'REG3G',	'ITLN1',	'SPINK4',	'MT1',	'FABP1',	'FABP2',	'FABP4',	'MMP7',	'PLA2G12A', 'MUC2',	'CLCA1',	'PDIA4',	'PDIA5',	'ZG16',	'REG4',	'REG3B',	'ANG4', 'CASP6',	'CASP8',	'REG1',	'REG1',	'REG3A',	'REG3B',	'REG4'))%>% convert_human_to_mouse_symbols()         
  genes_hormone_secretion <-  unique(c('CHGA',	'GCG',	'GHRL',	'TRPA1',	'NTS',	'AFP',	'VEGFA',	'PCSK1',	'SCN3A',	'NEUROG3',	'NEUROD1',	'NKX2-2',	'FABP5',	'FABP1',	'FABP6',	'GIP',	'GLP1R',	'GLP2R',	'PYY',	'GAST',	'SCGN'))%>% convert_human_to_mouse_symbols()
  genes_nutrient_absorption <-  unique(c('SLC7A8',	'FABP1',	'CPS1',	'VDR',	'APOA4',	'RBP2',	'GSTA1',	'SLC5A1',	'SLC7A9',	'EHF',	'PPP1R14B',	'PP1R14D')) %>% convert_human_to_mouse_symbols()
  
  genes_stem_prolif <- genes_stem_prolif[!is.na(genes_stem_prolif)]
  genes_antimicrobial <- genes_antimicrobial[!is.na(genes_antimicrobial)]
  genes_hormone_secretion <- genes_hormone_secretion[!is.na(genes_hormone_secretion)]
  genes_nutrient_absorption <- genes_nutrient_absorption[!is.na(genes_nutrient_absorption)]
  
  All_Genes <- UPAR_obj@assays$RNA@data@Dimnames[[1]]
  genes_stem_prolif <- intersect(All_Genes, genes_stem_prolif)
  genes_antimicrobial <- intersect(All_Genes, genes_antimicrobial)   
  genes_hormone_secretion <- intersect(All_Genes, genes_hormone_secretion)  
  genes_nutrient_absorption <- intersect(All_Genes, genes_nutrient_absorption)
  
  UPAR_obj <- AddModuleScore(object = UPAR_obj, features = list(genes_stem_prolif), name = 'stem_prolif')
  UPAR_obj <- AddModuleScore(object = UPAR_obj, features = list(genes_antimicrobial), name = 'antimicrobial')
  UPAR_obj <- AddModuleScore(object = UPAR_obj, features = list(genes_hormone_secretion), name = 'hormone_secretion')
  UPAR_obj <- AddModuleScore(object = UPAR_obj, features = list(genes_nutrient_absorption), name = 'nutrient_absorption')
  
  Idents(UPAR_obj) <- UPAR_obj$Cell_Type
  Cell_Types <- levels(UPAR_obj$Cell_Type)
  stem_prolif_FC <- list()
  nutrient_absorption_FC <- list()
  antimicrobial_FC <- list()
  hormone_secretion_FC <- list()
  
  for(i in Cell_Types){
    stem_prolif_FC[[i]] <- FoldChange(UPAR_obj,ident.1 = "UPAR+", ident.2 = "UPAR-",features = genes_stem_prolif, group.by = 'Treatment', subset.ident = i)
  }
  stem_prolif_FC_aggr <- c()
  for(i in Cell_Types){
    stem_prolif_FC_aggr <- c(stem_prolif_FC_aggr, mean(as.numeric(stem_prolif_FC[[i]]$avg_log2FC)))
  }
  
  for(i in Cell_Types){
    nutrient_absorption_FC[[i]] <- FoldChange(UPAR_obj,ident.1 = "UPAR+", ident.2 = "UPAR-",features = genes_nutrient_absorption, group.by = 'Treatment', subset.ident = i)
  }
  nutrient_absorption_FC_aggr <- c()
  for(i in Cell_Types){
    nutrient_absorption_FC_aggr <- c(nutrient_absorption_FC_aggr, mean(as.numeric(nutrient_absorption_FC[[i]]$avg_log2FC)))
  }
  
  for(i in Cell_Types){
    antimicrobial_FC[[i]] <- FoldChange(UPAR_obj,ident.1 = "UPAR+", ident.2 = "UPAR-",features = genes_antimicrobial, group.by = 'Treatment', subset.ident = i)
  }
  antimicrobial_FC_aggr <- c()
  for(i in Cell_Types){
    antimicrobial_FC_aggr <- c(antimicrobial_FC_aggr, mean(as.numeric(antimicrobial_FC[[i]]$avg_log2FC)))
  }
  
  for(i in Cell_Types){
    hormone_secretion_FC[[i]] <- FoldChange(UPAR_obj,ident.1 = "UPAR+", ident.2 = "UPAR-",features = genes_hormone_secretion, group.by = 'Treatment', subset.ident = i)
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
  Idents(UPAR_obj) <- UPAR_obj$Cell_Type
  pvals <- c()
  for(j in scores){
    for(i in Cell_Types){
      selected_cells <- WhichCells(object = UPAR_obj ,idents = c(i))
      vln_data <- FetchData(UPAR_obj,
                            vars = c(j,"Treatment", "Cell_Type"),
                            cells = selected_cells,
                            slot = "data")
      vln_data <- melt(vln_data)
      vln_data <- vln_data[c('Treatment', 'Cell_Type', 'value')]
      
      pvals <- c(pvals, t.test(vln_data[which(vln_data$Treatment == 'UPAR+'),]$value, vln_data[which(vln_data$Treatment == 'UPAR-'),]$value,  alternative = c("greater"))$p.value)
      print(i)
    }
  }
  pvals <- p.adjust(pvals, "BH", length(pvals))
  Module_Heatmap1$pval <- -log(pvals) 
  
  my_levels <- Cell_Types#c('Enterocyte', 'Enteroendocrine', 'Goblet', 'Paneth')
  Module_Heatmap1$Cell_Type <- factor(Module_Heatmap1$Cell_Type, levels = my_levels)
  max(Module_Heatmap1$value)
  pathwayColorsDiff = rev(brewer.pal(9, "RdBu"))
  ggplot(Module_Heatmap1, aes(Cell_Type, variable, fill= value)) + geom_tile() + 
    scale_fill_gradientn(name = "Log2FC", colors = pathwayColorsDiff, limits = c(-.4,.4)) + coord_flip() +
    theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45), axis.text.y =   element_text( vjust = 0.1, hjust = 1), panel.background = element_blank(), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())
  
  ggsave(file="UPAR_sort_score_heatmap.pdf",  width=3, height=3, units="in")
  
  ggplot(Module_Heatmap1, aes(x = Cell_Type, y = variable, size = pval ,color= value)) + geom_point() + theme_vyom + 
    scale_color_gradientn(name = "Log2FC", colors = pathwayColorsDiff, limits = c(-.4,.4)) + coord_flip() + labs(size="-log(p)") +
    theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45), axis.text.y =   element_text( vjust = 0.1, hjust = 1), panel.background = element_blank(), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())
  
  ggsave(file="UPAR_sort_score_dotplot.pdf",  width=4, height=5, units="in")
}

{
  Gene_list <- c('H2-Aa', 'H2-Ab1', 'H2-DMb2', 'H2-Q7', 'Cd74', 'Tap1', 'Ifngr1', 'Ccl5', 'Ccl3', 'Il16', 'Il18')

  All_Genes <- UPAR_obj@assays$RNA@data@Dimnames[[1]]
  Gene_list <- intersect(All_Genes, Gene_list)
  
  Idents(UPAR_obj) <- UPAR_obj$Cell_Type
  Cell_Types <- levels(UPAR_obj$Cell_Type)
  UPAR_obj$Treatment
  UPAR_obj$celltype <- UPAR_obj$Cell_Type
  DE_tControlle_all = Idents(UPAR_obj) %>% levels() %>% intersect(Cell_Types) %>% lapply(get_lfc_celltype, seurat_obj = UPAR_obj, condition_colname = "Treatment", condition_oi = 'UPAR+', condition_reference = 'UPAR-', celltype_col = NULL, expression_pct = 0.01) %>% purrr::reduce(full_join)    
  DE_tControlle_all[is.na(DE_tControlle_all)] = 0
  
  DE_tControlle_all_frame <- data.frame(DE_tControlle_all)
  rownames(DE_tControlle_all_frame) <- DE_tControlle_all_frame$gene
  Gene_list1 <- c('H2.Aa', 'H2.Ab1', 'H2.DMb2', 'H2.Q7', 'Cd74', 'Tap1', 'Ifngr1', 'Ifitm3', 'Ifitm2', 'Ccl5', 'Ccl4', 'Ccl3', 'Il16', 'Il18')
  Gene_list1 <- c('Ccl4','Ccl5','Ccl9','Ccl20','Cd74','Cox4I1','Cox8A','Ctsb','Ctsd','H2.Ac','H2.Ab1','H2.Q7','Ifi47','Ifit1Bl1','Ifit2','Ifit3','Ifitm3','Irf7','Isg15','Stat2','Stat3','Tap1','Ucp2')
  
  DE_tControlle_all_filtered <- DE_tControlle_all_frame[Inflammatory_score_new,]
  DE_tControlle_all_filtered <- DE_tControlle_all_filtered[complete.cases(DE_tControlle_all_filtered), ]
  
  # make LFC heatmap
  lfc_matrix = DE_tControlle_all_filtered  %>% dplyr::select(-gene) %>% as.matrix() %>% magrittr::set_rownames(DE_tControlle_all_filtered$gene)
  rownames(lfc_matrix) = rownames(lfc_matrix) %>% make.names()
  
  vis_lfc = lfc_matrix
  
  colnames(vis_lfc) = vis_lfc %>% colnames() %>% make.names()
  vis_lfc <- data.frame(vis_lfc)
  vis_lfc$gene <- rownames(vis_lfc)
  
  vis_lfc1 <-melt(vis_lfc)
  vis_lfc2 <- vis_lfc1
  pathwayColorsDiff = rev(brewer.pal(9, "RdBu"))
  vis_lfc1$gene <- factor(vis_lfc1$gene,levels = Gene_list1)
  ggplot(vis_lfc1, aes(gene, variable, fill= value)) + geom_tile() + 
    scale_fill_gradientn(name = "Log2FC", colors = pathwayColorsDiff, limits = c(-.75,.8)) + coord_flip() +
    theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45), axis.text.y =   element_text( vjust = 0.1, hjust = 1), panel.background = element_blank(), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())
  
  Idents(UPAR_obj) <- UPAR_obj$Cell_Type
  pvals <- c()
  for(j in Gene_list){
    for(i in Cell_Types){
      selected_cells <- WhichCells(object = UPAR_obj ,idents = c(i))
      vln_data <- FetchData(UPAR_obj,
                            vars = c(j,"Treatment", "Cell_Type"),
                            cells = selected_cells,
                            slot = "data")
      vln_data <- melt(vln_data)
      vln_data <- vln_data[c('Treatment', 'Cell_Type', 'value')]
      
      pvals <- c(pvals, t.test(vln_data[which(vln_data$Treatment == 'UPAR+'),]$value, vln_data[which(vln_data$Treatment == 'UPAR-'),]$value,  alternative = c("greater"))$p.value)
      print(i)
    }
  }
  
  pvals <- p.adjust(pvals, "BH", length(pvals))
  vis_lfc2$pval <- -log(pvals) 
  
  my_levels <- c(Gene_list1)
  vis_lfc2$gene <- factor(vis_lfc2$gene, levels = my_levels)
  
  pathwayColorsDiff = rev(brewer.pal(9, "RdBu"))
  ggplot(vis_lfc2, aes(variable, gene, fill= value)) + geom_tile() + 
    scale_fill_gradientn(name = "Log2FC", colors = pathwayColorsDiff, limits = c(-.75,.8)) + coord_flip() +
    theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45), axis.text.y =   element_text( vjust = 0.1, hjust = 1), panel.background = element_blank(), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())
  
  ggsave(file="MHC_UPAR_heatmap.pdf",  width=3, height=3, units="in")
  
  ggplot(vis_lfc2, aes(x = gene, y = variable, size = pval ,color= value)) + geom_point() + theme_vyom + 
    scale_color_gradientn(name = "Log2FC", colors = pathwayColorsDiff, limits = c(-.75,.8)) + coord_flip() + labs(size="-log(p)") +
    theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45), axis.text.y =   element_text( vjust = 0.1, hjust = 1), panel.background = element_blank(), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())
  
  ggsave(file="MHC_UPAR_dotplot.pdf",  width=3.5, height=4, units="in")
}  


#Expression umaps for each score
score_lsit <- c("senescence", "stem_score", "MHC_score", "Hormone" , "Absorption", "antimicrobial", "Inflammatory_score_new", "sen_mayo", "sen_BRD4" )
for(i in score_lsit) {
  FeaturePlot(object = UPAR_obj, features = i, pt.size = .001, label = TRUE, label.size = 4, repel = TRUE) +
    theme(plot.title = element_blank(), text = element_text(size=6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), axis.title.x = element_blank(), axis.title.y = element_blank()) +
    scale_color_gradientn(colors=brewer.pal(n = 9, name = "YlGnBu"))
  ggsave(file = paste0("UPAR_sort_feature_Umap_",i,".pdf"), width=5, height=5, units="in")
}

#find senescent cells
UPAR_obj$sen_cells <- 'non'
UPAR_obj <- AddModuleScore(object = UPAR_obj, features = list(sen_mayo), name = 'sen_mayo')

# module score distribution
Senescence_module <- as.data.frame(UPAR_obj$sen_mayo1)
modulescores <- Senescence_module %>%
  rownames_to_column(var="id") %>%
  pivot_longer(-id, names_to="celltype", values_to="score")


p <- ggplot(modulescores)
#p <- ggplot(onescore)
inflection_plot <- p + geom_point(aes(x=fct_inorder(id), y=sort(score)))+ scale_y_continuous(breaks = seq(-.1, .3, by = .01)) +
  facet_wrap(~celltype) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
inflection_plot

sen_cells <-WhichCells(object = UPAR_obj, expression = sen_mayo1 > .02)
DimPlot(UPAR_obj, label=T,group.by = 'Treatment', cells.highlight= list(sen_cells),  cols.highlight = c("darkblue"),cols= "grey")

UPAR_obj$sen_cells[sen_cells] <- paste('senescent')

Prop_table<- prop.table(x = table(UPAR_obj$Treatment, UPAR_obj$sen_cells), margin = 2)
Prop_Table <- as.data.frame(Prop_table, row.names = NULL, optional = FALSE,
                            make.names = TRUE, stringsAsFactors = default.stringsAsFactors())
Prop_Table <- Prop_Table[Prop_Table$Var2 == 'senescent',]
melt(Prop_Table)
ggplot(Prop_Table, aes(fill=Var1, y=Freq, x=Var2)) + 
  geom_bar(position="fill", stat="identity")+theme_vyom +scale_fill_manual(values = c("#7CA1CC", "#FF4902"))
ggsave( "UPAR_sort_senescent_plot.pdf", width=2.5, height=3, units="in")

#subset to just senescent cells
Idents(UPAR_obj) <- UPAR_obj$sen_cells
senenscent_scrna <- subset(UPAR_obj,  idents = c('senescent'))

Prop_table<- prop.table(x = table(senenscent_scrna$Treatment, senenscent_scrna$Cell_Type ), margin = 2)
Prop_Table <- as.data.frame(Prop_table, row.names = NULL, optional = FALSE,
                            make.names = TRUE, stringsAsFactors = default.stringsAsFactors())
Prop_Table1 <- Prop_Table
my_levels <- c('Stem', 'Transit Amplifying', 'Enterocyte Progenitor','Enterocyte', 'Enteroendocrine', 'Goblet', 'Paneth','Tuft', 'T Cell', 'Macrophage')
my_levels1 <- c('UPAR-','UPAR+')

Prop_Table1$Var1 <- factor(Prop_Table1$Var1,levels = my_levels1)
Prop_Table1$Var2 <- factor(Prop_Table1$Var2,levels = my_levels)
plot <- ggplot(data = Prop_Table1, aes(Var2, Freq, fill=Var1)) + geom_bar(position="stack", stat="identity", na.rm = TRUE) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), panel.background = element_rect(fill = "white", colour = "Black")) +
  xlab("Cell Type") + ylab("Fraction of Cells") + labs(fill = "Treatment") + scale_fill_manual(values = c('#7CA1CC' ,'#FF4902'))
plot
ggsave(file = paste0('UPARsort_senescent_cell_proportion_percelltype.pdf'), width=4, height=3, units="in")

#housekeeping
# saveRDS(UPAR_obj, file = "./data/Seurat_Objects/UPAR_sort_SCRNA_Sobj.rds")
# UPAR_obj <- readRDS("./data/Seurat_Objects/UPAR_sort_SCRNA_Sobj.rds", refhook = NULL)

