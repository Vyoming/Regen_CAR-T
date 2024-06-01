#plotting functions

# Vyom Shah
# Cold Spring Harbor Laboratory
# Density Difference Expression Plot
library(plyr)
library(dplyr)
library(reshape2)
library(stringr)
library(RColorBrewer)
library(ggplot2)
library(matrixStats)
library(Seurat)
library(tidyverse)
library(ggrepel)
library(hdf5r)
library(scales)
library(future)
library(magrittr)
library(mgcv)
library(fgsea)

#load in themes/necessary functions/data
Seurat_Object
options(stringsAsFactors=FALSE)
colorsType = c(
  '25-30' = "#7CA1CC",
  '65-70' = "#FF4902"
)
colorsType = c(
  'PLAUR-' = "#7CA1CC",
  'PLAUR+' = "#FF4902"
)

theme_linedraw2 = theme_linedraw() + theme(strip.background=element_rect(fill="grey80", colour="grey50", size=0.2), strip.text.x=element_text(colour="black"), strip.text.y=element_text(colour="black"))
theme_vyom = theme_linedraw2 + theme(legend.position="right", legend.title=element_text(size=15), legend.text=element_text(size=14), axis.text.x = element_text(size=12, angle=-90, hjust=0, vjust=0.5), axis.text.y=element_text(size=12), axis.title=element_text(size=15), axis.title.y=element_text(vjust=1), plot.title = element_text(size=18, vjust=1.5), strip.background = element_rect(fill="#EEEEEE"), strip.text = element_text(size = 11), panel.grid.major = element_line(colour = "grey98"), panel.grid.minor = element_blank())

theme_linedraw2 = theme_linedraw() + theme(strip.background=element_rect(fill="grey80", colour="grey50", size=0.2), strip.text.x=element_text(colour="black"), strip.text.y=element_text(colour="black"))
theme_vyom = theme_linedraw2 + theme(legend.position="right", legend.title=element_text(size=15), legend.text=element_text(size=14), axis.text.x = element_text(size=12, angle=-90, hjust=0, vjust=0.5), axis.text.y=element_text(size=12), axis.title=element_text(size=15), axis.title.y=element_text(vjust=1), plot.title = element_text(size=18, vjust=1.5), strip.background = element_rect(fill="#EEEEEE"), strip.text = element_text(size = 11), panel.grid.major = element_line(colour = "grey98"), panel.grid.minor = element_blank())

findConcaveHull = function(df, X, Y, alpha, hullMargin=NA, extend=F, minPoint=3){
  library(igraph)
  library(alphahull)
  library(polyclip)
  
  if(nrow(df) < 3)
  {
    replicationMargin = 3*hullMargin/100
    df2 = rbind( data.frame(X=df[, X] + replicationMargin, Y=df[, Y] + replicationMargin), 
                 data.frame(X=df[, X] + replicationMargin, Y=df[, Y] - replicationMargin),
                 data.frame(X=df[, X] - replicationMargin, Y=df[, Y] + replicationMargin),
                 data.frame(X=df[, X] - replicationMargin, Y=df[, Y] - replicationMargin) )
  }else
  {
    df2 = data.frame(X=df[, X], Y=df[, Y])
  }
  
  if(minPoint <= 3) minPoint = 3
  
  if(is.na(hullMargin)) hullMargin = min( (max(df2$X, na.rm=T) - min(df2$X, na.rm=T) ), (max(df2$Y, na.rm=T) - min(df2$Y, na.rm=T)) ) * 0.05
  
  #shape = ashape(df2$X, df2$Y, alpha)
  
  ##plot(shape)
  #shapeGraph = graph.edgelist(cbind(as.character(shape$edges[, "ind1"]), as.character(shape$edges[, "ind2"])), directed = FALSE)
  ##plot(shapeGraph)
  
  shape2 = ahull(df2$X, df2$Y, alpha)	
  shapeArcs = shape2$arcs
  shapeArcs = shapeArcs[ shapeArcs[, "end1"] != shapeArcs[, "end2"], ]
  shapeGraph = graph.edgelist(cbind(as.character(shapeArcs[, "end1"]), as.character(shapeArcs[, "end2"])), directed = FALSE)
  
  graphComponents = decompose.graph(shapeGraph)
  edgesList = list()
  curIdx = 1
  for(curGraph in graphComponents)
  {
    if(length(V(curGraph)) < minPoint) next
    if(length(E(curGraph)) < minPoint) next
    #if (!is.connected(curGraph)) next
    #if (any(degree(curGraph) != 2)) next
    
    # find chain end points
    for(i in 1:length(E(curGraph)))
    {
      cutg = curGraph - E(curGraph)[i]
      ends = names(which(degree(cutg) == 1))
      if(length(ends) >= 2) break
    }
    
    if(length(ends) >= 2)
    {
      path = get.shortest.paths(cutg, ends[1], ends[2])[[1]]
      # this is an index into the points
      pathX = as.numeric(V(curGraph)[path[[1]]]$name)
      # join the ends
      pathX = c(pathX, pathX[1])
      edgesList[[curIdx]] = shape2$x[pathX, ]
      curIdx = curIdx + 1
    }else
    {
      # TODO: Edge case
      # Graph is shaped weirdly, do depth-first-traversal or similar way
      # dfs
    }
  }
  
  if(length(edgesList)==1)
  {
    curEdges = edgesList[[1]]
    allEdges = as.data.frame(curEdges)
    colnames(allEdges) = c("x", "y")
    allEdges$group = 1
    if(extend == T)
    {
      extended = polyoffset(list(x=allEdges[, 1], y=allEdges[, 2]), hullMargin, jointype="round", arctol=abs(hullMargin)/40)
      extended2 = data.frame(x=extended[[1]]$x, y=extended[[1]]$y, group=1)
      return(extended2)
    }else
    {
      return(allEdges)
    }
  }else if(length(edgesList) > 1)
  {
    mergedEdges = edgesList[[1]]
    mergedEdges = mergedEdges[1:(nrow(mergedEdges)-1), ]
    otherEdges = list()
    for(i in 2:length(edgesList))
    {
      curEdges = edgesList[[i]]
      curEdges = curEdges[1:(nrow(curEdges)-1), ]
      otherEdges[[i-1]] = list(x=curEdges[,1], y=curEdges[,2])
    }
    
    mergedShape = polyclip(list(x=mergedEdges[,1], y=mergedEdges[,2]), otherEdges, op="xor")
    
    mergedShape2 = data.frame()
    for(j in 1:length(mergedShape))
    {
      mergedShape[[j]]$x = c(mergedShape[[j]]$x, mergedShape[[j]]$x[1]) ## Extend the last point, otherwise ggplot freaks out with holes
      mergedShape[[j]]$y = c(mergedShape[[j]]$y, mergedShape[[j]]$y[1])
      if(extend == F)
      {
        mergedShape2 = rbind(mergedShape2, data.frame(x=mergedShape[[j]]$x, y=mergedShape[[j]]$y, group=j))
      }else
      {			
        mergedShapeEx = polyoffset(mergedShape[[j]], hullMargin/2, jointype="round", arctol=abs(hullMargin)/40)
        mergedShapeDF = data.frame(x=mergedShapeEx[[1]]$x, y=mergedShapeEx[[1]]$y, group=j)
        #mergedShapeDF = simplifyDF(mergedShapeDF, "x", "y", hullMargin, T)
        mergedShapeDF = data.frame(x=mergedShapeDF$x, y=mergedShapeDF$y, group=j)
        mergedShape2 = rbind(mergedShape2, mergedShapeDF)
      }
    }
    
    allEdges = as.data.frame(mergedShape2)
    allEdges = fixfeature(allEdges)
    #ggplot(allEdges, aes(x=x, y=y, group=group)) + geom_polygon()
    
    if(extend == T)
    {
      extended = polyoffset(list(x=allEdges$x, y=allEdges$y), hullMargin/2, jointype="round", arctol=abs(hullMargin)/40)
      extended = data.frame(x=extended[[1]]$x, y=extended[[1]]$y, group=1)
      return(extended)
    }else
    {
      return(data.frame(x=allEdges$x, y=allEdges$y, group=allEdges$group))
    }
  }else
  {
    return(data.frame(x=numeric(0), y=numeric(0), group=numeric(0)))
  }
}
concaveHull = function(df, X, Y, group, alpha, margin=NA, extend=F, minPoint = 3, marginMultiplier = 1){
  if(is.na(margin))  margin = min( (max(df[, X], na.rm=T) - min(df[, X], na.rm=T) ), (max(df[, Y], na.rm=T) - min(df[, Y], na.rm=T)) ) * 0.02 * marginMultiplier
  hulls = ddply(df, group, findConcaveHull, X, Y, alpha, margin, extend, minPoint)
  #hulls2 = ddply(hulls, group, extendHull, "X", "Y", margin)
  # hulls2 = ddply(hulls, group, smoothHull, "x", "y") ## if using bezier jointtype should be square
  
  return(hulls)
}
fixfeature <- function(df){
  ringstarts <- which(!duplicated(df$group))
  if(length(ringstarts) < 2) {
    return(df)
  } else {
    ringstarts <- c(ringstarts, nrow(df))
    indicies <- c(1:(ringstarts[2]-1), do.call(c, lapply(2:(length(ringstarts)-1), function(x) {
      c(1, ringstarts[x]:(ringstarts[x+1]-1))
    })), nrow(df))
    return(df[indicies,])
  }
}
smoothScore2d = function(score, x, y, type=NULL, numGrid = 100, knn = 100, m = 2, expand=0.05, xrng=NULL, yrng=NULL){
  library(ash)
  curDat = data.frame(score=score, x=x, y=y)
  if(is.null(xrng)) xrng = range(curDat$x)
  if(is.null(yrng)) yrng = range(curDat$y)
  xdiff = xrng[2] - xrng[1]
  ydiff = yrng[2] - yrng[1]
  xrng[1] = xrng[1] - xdiff*expand
  xrng[2] = xrng[2] + xdiff*expand
  yrng[1] = yrng[1] - ydiff*expand
  yrng[2] = yrng[2] + ydiff*expand
  bins = bin2(cbind(curDat$x, curDat$y), ab = rbind(xrng, yrng), nbin = c(numGrid, numGrid))
  binCounts = ash2(bins, m = c(m, m))
  gridDat = data.frame( expand.grid( x = binCounts$x, y = binCounts$y), density = melt(binCounts$z)[,3] )
  gridDat2 = gridDat[ gridDat$density > 0, ]
  
  if(is.null(type))
  {
    return( smoothScorePredict(curDat, gridDat2, knn = knn) )
  }else
  {
    curDat$type = type
    allPredicted = ddply(curDat, "type", function(x) { smoothScorePredict( x[, colnames(x) != "type"], gridDat2, knn = knn) } )
    return(allPredicted)
  }
}
smoothScorePredict = function(trainDat, newGrid, knn=30){
  scorePredicted = data.frame(score = as.numeric(NA), 
                              x = newGrid[, 1], 
                              y = newGrid[, 2])
  # run KKNN
  scoreKknn = kknn::kknn(score ~ ., 
                         train = trainDat, 
                         test = scorePredicted, 
                         kernel = "gaussian", 
                         k = knn)
  
  scorePredicted %<>% mutate(score = fitted(scoreKknn))
  return(scorePredicted)
}
sen_data1@assays$RNA@scale.data
calculatePathwayScores = function(seuratObj, dimred, pathwayList, method="seurat", minCells=10,  minGenes=2)
{
  curData = seuratObj@assays$RNA@scale.data
  curData = curData[ rowSums(curData > 0) >= minCells, ]
  curData = curData[ rowSds(curData) > 0, ]
  zscore = t(curData) # scale(t(curData))
  pathwayNames = unique( gsub("(_UP|_DN|_downreg|_upreg|_down_|_up_)", "\\*", names(pathwayList)) )
  
  pathwayScores = dimred
  for(curPathway in pathwayNames)
  {
    curPathwayName = gsub("\\*", "", curPathway)
    curUpSum = rep(0, nrow(zscore))
    curDnSum = rep(0, nrow(zscore))
    anyHits = F
    curUpGenesList = c()
    curDnGenesList = c()
    
    if(grepl("\\*", curPathway) & gsub("\\*", "_UP", curPathway) %in% names(pathwayList))
    {
      curUpGenesList = c(curUpGenesList, pathwayList[[ gsub("\\*", "_UP", curPathway) ]] )
    }
    if(grepl("\\*", curPathway) & gsub("\\*", "_upreg", curPathway) %in% names(pathwayList))
    {
      curUpGenesList = c(curUpGenesList, pathwayList[[ gsub("\\*", "_upreg", curPathway) ]] )
    }
    if(grepl("\\*", curPathway) & gsub("\\*", "_up_", curPathway) %in% names(pathwayList))
    {
      curUpGenesList = c(curUpGenesList, pathwayList[[ gsub("\\*", "_up_", curPathway) ]] )
    }
    
    if(grepl("\\*", curPathway) & gsub("\\*", "_DN", curPathway) %in% names(pathwayList))
    {
      curDnGenesList = c(curDnGenesList, pathwayList[[ gsub("\\*", "_DN", curPathway) ]] )
    }
    if(grepl("\\*", curPathway) & gsub("\\*", "_downreg", curPathway) %in% names(pathwayList))
    {
      curDnGenesList = c(curDnGenesList, pathwayList[[ gsub("\\*", "_downreg", curPathway) ]] )
    }
    if(grepl("\\*", curPathway) & gsub("\\*", "_down_", curPathway) %in% names(pathwayList))
    {
      curDnGenesList = c(curDnGenesList, pathwayList[[ gsub("\\*", "_down_", curPathway) ]] )
    }
    if( length(curUpGenesList) == 0 & length(curDnGenesList) == 0 & curPathway %in% names(pathwayList))
    {
      curUpGenesList = c(curUpGenesList, pathwayList[[ curPathway ]] )
    }
    curUpGenesList = intersect(curUpGenesList, colnames(zscore) )
    curDnGenesList = intersect(curDnGenesList, colnames(zscore) )
    
    if( length(curUpGenesList) + length(curDnGenesList) <= minGenes ) # min 10 genes
    {
      next
    }
    
    if(method == "seurat")
    {
      curUpGenesList = list( curUpGenesList )
      curDnGenesList = list( curDnGenesList )
      seuratObjUp = NULL
      seuratObjDn = NULL
      tryCatch({seuratObjUp = AddModuleScore(object = seuratObj, features = curUpGenesList, ctrl = 5, name = 'Score'); curUpSum = seuratObjUp$Score1; }, error = function(e) { print(e) })
      tryCatch({seuratObjDn = AddModuleScore(object = seuratObj, features = curDnGenesList, ctrl = 5, name = 'Score'); curDnSum = seuratObjDn$Score1; }, error = function(e) { print(e) })
      curScore = curUpSum - curDnSum 
      if(is.null(names(curScore))) next
      
    }else
    {
      if( length(curUpGenesList) > 0)
      {
        curUp = zscore[, colnames(zscore) %in% curUpGenesList, drop=F]
        curUpSum = rowSums(curUp)
      }
      if( length(curDnGenesList) > 0)
      {
        curDn = zscore[, colnames(zscore) %in% curDnGenesList, drop=F]
        curDnSum = rowSums(curDn)
      }
      
      curScore = curUpSum - curDnSum
      names(curScore) = rownames(zscore)
    }
    
    curScore = data.frame(Cell=names(curScore), Score=scale(curScore))
    colnames(curScore)[2] = curPathwayName
    pathwayScores = merge(pathwayScores, curScore, by="Cell")	
  }
  
  pathwayScores2 = melt(pathwayScores, colnames(dimred), variable.name="Pathway", value.name="Score")
  pathwayScores2$Pathway = as.character(pathwayScores2$Pathway)
  return(pathwayScores2)
}


dimred = data.frame(Seurat_Object@reductions$umap@cell.embeddings)
dimred$Cell = rownames(dimred)
all.equal(names(Seurat_Object$seurat_clusters), dimred$Cell)
all.equal(names(Seurat_Object$Treatment), dimred$Cell)
#Seurat_Object$Cell_Type <- Seurat_Object$CellType
dimred$Cluster = Seurat_Object$Cell_Type # Overwrite clusters with cell types for current plots
dimred$Treatment = Seurat_Object$Treatment
dimred$CellType = dimred$Cluster
assay = Matrix::t(Seurat_Object@assays$RNA$data)
dimred$Type = dimred$Treatment

getPathways = function(genesOrth = NULL)
{
  pathwaysM = c(gmtPathways("~/analysis/scRNA_Intestine/Beyaz_AA_final_organoid.gmt"))
  return(pathwaysM)
}
pathways = getPathways()
stem_score <- c('Lgr4','Myc', 'Sox9', 'Olfm4', 'Hopx', 'Ccnd1' )

clusterMedian = dimred %>%
  group_by(CellType) %>%
  summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2))

hulls2 = concaveHull(dimred, "UMAP_1", "UMAP_2", "CellType", alpha = 0.5, extend = T, minPoint = 20)

#use if calculating seperate genes
geneList = intersect(geneList, colnames(assay))
dimredG = merge(dimred, as.matrix(assay[, geneList, drop=F]), by.x="Cell", by.y="row.names")

#use if calculating score
getPathways = function(genesOrth = NULL)
{
  pathwaysM = c(gmtPathways("~/analysis/scRNA_Intestine/Beyaz_AA_final_organoid.gmt"))
  return(pathwaysM)
}
pathways = getPathways()
dimredG <- dimred
pathways$stem_score <- stem_score
pathways$MHC_score <- MHC_score
pathways$Inflammatory_score_new <- Inflammatory_score_new
pathways$sen_mayo <- sen_mayo

pathwayNames = gsub("(_UP|_DN|_downreg|_upreg|_down_|_up_)", "\\*", names(pathways))
pathwaysSelected = pathways[ pathwayNames %in% c("stem_score") ]

#pathwaysSelected = pathways

pathwayScoresZscore = calculatePathwayScores(Seurat_Object, dimred, pathwaysSelected, method="seurat")
All_Genes <- Seurat_Object@assays$RNA@data@Dimnames[[1]]
gene.list <- intersect(All_Genes, Creb_Targets_vs)

pathwayScores = pathwayScoresZscore

pathwayScoresW = dcast(pathwayScores, Cell+Cluster+Type+CellType+UMAP_1+UMAP_2~Pathway, value.var="Score")
#colnames(pathwayScoresW)[9:ncol(pathwayScoresW)] = paste0("Score_", colnames(pathwayScoresW)[9:ncol(pathwayScoresW)])

dimredG = pathwayScoresW
geneList <- gsub("(_UP|_DN|_downreg|_upreg|_down_|_up_)", "\\*", names(pathwaysSelected))

#colors
geneColors = rev(viridis::magma(10))#brewer.pal(9, "YlOrRd")[-c(1)] ## brewer.pal(9, "YlOrRd") ## rev(viridis::magma(10))
geneColorsDiff = rev(brewer.pal(9, "RdBu"))

curFilename = "CART"
curTitle = "UPAR vs UT"
Control_ident ="UT"
Experimental_ident = "UPAR"
for(curGene in geneList){
  rangeL = quantile(dimredG[, curGene], 0.01, na.rm = T)
  rangeH = quantile(dimredG[, curGene], 0.95, na.rm = T)
  
  if(rangeL == rangeH) { rangeL = min(dimredG[, curGene]); rangeH = max(dimredG[, curGene]); }
  
  sm = smoothScore2d(dimredG[, curGene], dimredG$UMAP_1, dimredG$UMAP_2, numGrid=100, knn=50, m=2)
  
  sm2 = ddply(dimredG, "Type", function(x) { smoothScore2d( x[,  curGene], x$UMAP_1, x$UMAP_2, numGrid=100, knn=50, m=2, xrng=range(dimredG$UMAP_1), yrng=range(dimredG$UMAP_2)) } )	
  
  sm3 = smoothScore2d(dimredG[, curGene], dimredG$UMAP_1, dimredG$UMAP_2, type=dimredG$Type, numGrid=100, knn=50, m=2)
  sm3W = dcast(sm3, x+y~type, value.var="score")
  sm3W$PMXS_vs_AB = sm3W[[Experimental_ident]] - sm3W[[Control_ident]]
  rangeDiff = max(c(abs(sm3W$PMXS_vs_AB)))
  
  ggplot(sm3W) + geom_tile(aes(x = x, y = y, fill = PMXS_vs_AB)) + theme_vyom +
    scale_fill_gradientn(name = "Log2fc", colors = geneColorsDiff, limits = c(-rangeDiff, rangeDiff)) +
    geom_polygon(data = hulls2, aes(x = x, y = y, group = CellType), alpha = 0.3, fill=NA, color="#666666", size=0.1) + ylab('UMAP_2') + xlab('UMAP_1') +
    geom_text_repel(data = clusterMedian, aes(x = UMAP_1, y = UMAP_2, label = CellType, group = as.factor(CellType)), size = 1.5, box.padding = .1, force = 5, color = "Black") +
    theme(legend.position = "none", text = element_text(size=5),  axis.text.x = element_text(size =  6), axis.text.y = element_text(size = 6), axis.title.x = element_text(size = 7), axis.title.y = element_text(size = 7))
  
  ggsave(file= paste0(curFilename, "_Expression_Desnity_Diff_", curGene, ".pdf"), width=2.5, height=2.5, units="in")
  
  ggplot(sm3W) + geom_tile(aes(x = x, y = y, fill = PMXS_vs_AB)) + theme_vyom +
    scale_fill_gradientn(name = "Log2fc", colors = geneColorsDiff, limits = c(-rangeDiff, rangeDiff)) +
    geom_polygon(data = hulls2, aes(x = x, y = y, group = CellType), alpha = 0.3, fill=NA, color="#666666", size=0.1) + ylab('UMAP_2') + xlab('UMAP_1') +
    geom_text_repel(data = clusterMedian, aes(x = UMAP_1, y = UMAP_2, label = CellType, group = as.factor(CellType)), size = 1.5, box.padding = .1, force = 5, color = "Black") +
    theme( text = element_text(size=5),  axis.text.x = element_text(size =  6), axis.text.y = element_text(size = 6), axis.title.x = element_text(size = 7), axis.title.y = element_text(size = 7))
  
  ggsave(file= paste0(curFilename, "_LABEL_Expression_Desnity_Diff_LABEL_", curGene, ".pdf"), width=2.5, height=2.5, units="in")
  
}

{
  library(MASS)
  library(reshape2)
  library(scales)
  library(ggrastr)
  # Calculate the common x and y range for geyser1 and geyser2
  xrng = range(dimred$UMAP_1)
  yrng = range(dimred$UMAP_2)
  
  extendRange = 0.06 # extend range by %6
  xDiff = (xrng[2] - xrng[1])
  yDiff = (yrng[2] - yrng[1])
  
  xrng[1] = xrng[1] - xDiff* extendRange
  xrng[2] = xrng[2] + xDiff * extendRange
  yrng[1] = yrng[1] - yDiff * extendRange
  yrng[2] = yrng[2] + yDiff * extendRange
  
  
  clusterMedian = dimred %>%
    group_by(CellType) %>%
    summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2))
  
  hulls2 = concaveHull(dimred, "UMAP_1", "UMAP_2", "CellType", alpha = 0.5, extend = T, minPoint = 20)
  
  
  
  for(caseType in unique(dimred$Type))
  {
    for(ctrlType in unique(dimred$Type))
    {
      if(caseType == ctrlType) next
      
      caseVsCtrlName = paste0(caseType, "_vs_", ctrlType)
      caseVsCtrlName <- gsub(" ", "_", caseVsCtrlName)
      colorLow = colorsType[ctrlType]
      colorMid = "white"
      colorHigh = colorsType[caseType]
      
      d_Case = kde2d( dimred$UMAP_1[dimred$Type == caseType], dimred$UMAP_2[dimred$Type == caseType ], lims = c(xrng, yrng), n = 500)
      d_Ctrl = kde2d( dimred$UMAP_1[dimred$Type == ctrlType  ], dimred$UMAP_2[dimred$Type == ctrlType ], lims = c(xrng, yrng), n = 500)
      
      # Confirm that the grid points for each density estimate are identical
      identical(d_Case$x, d_Ctrl$x) # TRUE
      identical(d_Case$y, d_Ctrl$y) # TRUE
      
      # Calculate the difference between the 2d density estimates
      diff_CaseVsCtrl = d_Ctrl 
      diff_CaseVsCtrl$z = d_Case$z - d_Ctrl$z
      diff_CaseVsCtrl$z = diff_CaseVsCtrl$z / max(diff_CaseVsCtrl$z)
      
      rownames(diff_CaseVsCtrl$z) = diff_CaseVsCtrl$x
      colnames(diff_CaseVsCtrl$z) = diff_CaseVsCtrl$y
      
      diff_CaseVsCtrlM = melt(diff_CaseVsCtrl$z, id.var = rownames(diff_CaseVsCtrl))
      names(diff_CaseVsCtrlM) = c("UMAP_1", "UMAP_2", 'caseVsCtrlName')
      
      ggplot(diff_CaseVsCtrlM, aes(x = UMAP_1, y = UMAP_2)) +
        ggrastr::rasterise(geom_tile(aes_string(fill = 'caseVsCtrlName'), alpha = 1), dpi = 200) +
        scale_fill_gradient2(low = colorLow, mid = colorMid, high = colorHigh, midpoint = 0) +
        coord_cartesian(xlim = xrng, ylim = yrng) +
        scale_color_manual(values = colorsType) +
        guides(colour = FALSE) + theme_vyom +
        geom_polygon(data = hulls2, aes(x = x, y = y, group = CellType), alpha = 0.1, size=0.1, fill=NA, color="Grey50") + 
        geom_text_repel(data = clusterMedian, aes(label = CellType, group = as.factor(CellType)), size = 2, box.padding = .1, force = 5, color = "Black") + 
        ggtitle(paste0("Density comparison of ", caseType, " vs ", ctrlType)) + theme(text = element_text(size=5), legend.key.size = unit(0.0, "cm"), legend.text=element_text(size=0), legend.title =element_text(size=0), axis.text.x = element_text(size = 6), plot.title = element_blank(), axis.text.y = element_text(size = 6), axis.title.x = element_text(size = 7), axis.title.y = element_text(size = 7))
      ggsave(paste0(curFilename, "_dimred_densityDiff_cellType", caseVsCtrlName, ".pdf"), width = 3, height = 3, units="in")
      
      ggplot(diff_CaseVsCtrlM, aes(x = UMAP_1, y = UMAP_2)) +
        ggrastr::rasterise(geom_tile(aes_string(fill = 'caseVsCtrlName'), alpha = 1), dpi = 200) +
        scale_fill_gradient2(low = colorLow, mid = colorMid, high = colorHigh, midpoint = 0) +
        coord_cartesian(xlim = xrng, ylim = yrng) +
        scale_color_manual(values = colorsType) +
        guides(colour = FALSE) + theme_vyom +
        geom_polygon(data = hulls2, aes(x = x, y = y, group = CellType), alpha = 0.1, size=0.1, fill=NA, color="Grey50") + 
        geom_text_repel(data = clusterMedian, aes(label = CellType, group = as.factor(CellType)), size = 2, box.padding = .1, force = 5, color = "Black") + 
        ggtitle(paste0("Density comparison of ", caseType, " vs ", ctrlType)) + theme(text = element_text(size=5),  axis.text.x = element_text(size = 6), plot.title = element_blank(), axis.text.y = element_text(size = 6), axis.title.x = element_text(size = 7), axis.title.y = element_text(size = 7))
      
      ggsave(paste0(curFilename, "_dimred_densityDiff_cellType_label", caseVsCtrlName, ".pdf"), width = 4, height = 4, units="in")
    }
  }
}


