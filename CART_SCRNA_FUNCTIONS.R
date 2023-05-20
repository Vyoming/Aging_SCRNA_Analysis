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
Cart_UT$Age
#load in themes/necessary functions/data
Seurat_Object <- Cart_UT
options(stringsAsFactors=FALSE)
colorsType = c(
  UT = "#7CA1CC",
  uPAR = "#FF4902"
)

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

dimred = data.frame(Seurat_Object@reductions$umap@cell.embeddings)
dimred$Cell = rownames(dimred)
all.equal(names(Seurat_Object$seurat_clusters), dimred$Cell)
all.equal(names(Seurat_Object$Treatment), dimred$Cell)

dimred$Cluster = Seurat_Object$Cell_Type # Overwrite clusters with cell types for current plots
dimred$Treatment = Seurat_Object$Age
dimred$CellType = dimred$Cluster
assay = Matrix::t(Seurat_Object@assays$RNA@data)
dimred$Type = dimred$Treatment
curFilename = "CART_all"
curTitle = "CART_all"
Seurat_Object$senescence
geneList = unique(c('senescence'))
geneList = unique(c( 'Cxcr3','Ctsb','Ccl2','Jun','Irf7','ppar_mouse'))

clusterMedian = dimred %>%
  group_by(CellType) %>%
  summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2))

hulls2 = concaveHull(dimred, "UMAP_1", "UMAP_2", "CellType", alpha = 0.5, extend = T, minPoint = 20)

geneList = intersect(geneList, colnames(assay))
dimredG = merge(dimred, as.matrix(assay[, geneList, drop=F]), by.x="Cell", by.y="row.names")

#if you are using score instead
dimredG <- dimred
dimredG$senescence <- Seurat_Object$senescence


geneColors = rev(viridis::magma(10))#brewer.pal(9, "YlOrRd")[-c(1)] ## brewer.pal(9, "YlOrRd") ## rev(viridis::magma(10))
geneColorsDiff = rev(brewer.pal(9, "RdBu"))
for(curGene in geneList)
{
  rangeL = quantile(dimredG[, curGene], 0.01, na.rm = T)
  rangeH = quantile(dimredG[, curGene], 0.95, na.rm = T)
  
  if(rangeL == rangeH) { rangeL = min(dimredG[, curGene]); rangeH = max(dimredG[, curGene]); }
  
  sm = smoothScore2d(dimredG[, curGene], dimredG$UMAP_1, dimredG$UMAP_2, numGrid=100, knn=50, m=2)
  
  sm2 = ddply(dimredG, "Treatment", function(x) { smoothScore2d( x[,  curGene], x$UMAP_1, x$UMAP_2, numGrid=100, knn=50, m=2, xrng=range(dimredG$UMAP_1), yrng=range(dimredG$UMAP_2)) } )	
  
  sm3 = smoothScore2d(dimredG[, curGene], dimredG$UMAP_1, dimredG$UMAP_2, type=dimredG$Treatment, numGrid=100, knn=50, m=2)
  sm3W = dcast(sm3, x+y~type, value.var="score")
  sm3W$PMXS_vs_AB = sm3W$Old - sm3W$Young
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
  
  ggsave(file= paste0(curFilename, "_Expression_Desnity_Diff_LEGEND_", curGene, ".pdf"), width=2.5, height=2.5, units="in")
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
      names(diff_CaseVsCtrlM) = c("UMAP_1", "UMAP_2", caseVsCtrlName)
      
      ggplot(diff_CaseVsCtrlM, aes(x = UMAP_1, y = UMAP_2)) +
        ggrastr::rasterise(geom_tile(aes_string(fill = caseVsCtrlName), alpha = 1), dpi = 200) +
        scale_fill_gradient2(low = colorLow, mid = colorMid, high = colorHigh, midpoint = 0) +
        coord_cartesian(xlim = xrng, ylim = yrng) +
        scale_color_manual(values = colorsType) +
        guides(colour = FALSE) + theme_vyom +
        geom_polygon(data = hulls2, aes(x = x, y = y, group = CellType), alpha = 0.1, size=0.1, fill=NA, color="Grey50") + 
        geom_text_repel(data = clusterMedian, aes(label = CellType, group = as.factor(CellType)), size = 2, box.padding = .1, force = 5, color = "Black") + 
        ggtitle(paste0("Density comparison of ", caseType, " vs ", ctrlType)) + theme(text = element_text(size=5), legend.key.size = unit(0.0, "cm"), legend.text=element_text(size=0), legend.title =element_text(size=0), axis.text.x = element_text(size = 6), plot.title = element_blank(), axis.text.y = element_text(size = 6), axis.title.x = element_text(size = 7), axis.title.y = element_text(size = 7))
      ggsave(paste0(curFilename, "_dimred_densityDiff_cellType", caseVsCtrlName, ".pdf"), width = 3, height = 3, units="in")
      
      ggplot(diff_CaseVsCtrlM, aes(x = UMAP_1, y = UMAP_2)) +
        ggrastr::rasterise(geom_tile(aes_string(fill = caseVsCtrlName), alpha = 1), dpi = 200) +
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


library(readxl)
Young_vs_Old_MHC_II_genes <- read_excel("Documents/Young vs Old MHC-II genes.xlsx")
View(Young_vs_Old_MHC_II_genes)
idents1 <- c('symbol','Young_ISC.norm','Old_ISC.norm')
Young_vs_Old_MHC_II_genes <- Young_vs_Old_MHC_II_genes[idents1]

rownames(Young_vs_Old_MHC_II_genes) <- Young_vs_Old_MHC_II_genes$symbol
fc <- c()
for(i in Young_vs_Old_MHC_II_genes$symbol){
  FC_ITER <- foldchange(Young_vs_Old_MHC_II_genes[i,]$Old_ISC.norm, Young_vs_Old_MHC_II_genes[i,]$Young_ISC.norm)
  fc <- c(fc,FC_ITER)
}
Young_vs_Old_MHC_II_genes$FC <- fc
Young_vs_Old_MHC_II_genes$Log2_FC <- log(fc, 2)
idents1 <- c('symbol','Log2_FC')
FC_data <- Young_vs_Old_MHC_II_genes[idents1]
library(reshape)
Heatmap_data <-melt(FC_data)
pathwayColorsDiff = rev(brewer.pal(9, "RdBu"))
ggplot(Heatmap_data, aes(symbol, variable, fill= value)) + geom_tile() + 
  scale_fill_gradientn(name = "Log2FC", colors = pathwayColorsDiff, limits = c(-3,3)) + coord_flip() +
  theme(axis.text.y =   element_text( vjust = 0.1, hjust = 1), panel.background = element_blank(), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())
ggsave(file="OldvsYoung_MHCII.pdf",  width=1.85, height=2, units="in")

Young_vs_Old_MHC_II_genes_mean
Young_vs_Old_MHC_II_genes
length(Young_vs_Old_MHC_II_genes_split$arasco[which(Young_vs_Old_MHC_II_genes_split$arasco$metabolite == i),]$measure)

zscores <- zscoreT(Young_vs_Old_MHC_II_genes$measure, df = length(Young_vs_Old_MHC_II_genes$measure))
aa_zscore <- Young_vs_Old_MHC_II_genes
aa_zscore$measure <- zscores

aa_zscore_mean <- aa_zscore %>%
  group_by(diet, metabolite) %>%
  summarise_at(vars(measure), list(name = mean))
aa_zscore_mean$diet

my_levels <- c('control','arasco')
aa_zscore_mean$diet <- factor(x = aa_zscore_mean$diet, levels = my_levels)

ggplot(Heatmap_data, aes(metabolite, variable, fill= value)) + geom_tile() + 
  scale_fill_gradientn(name = "Log2FC", colors = pathwayColorsDiff, limits = c(-1,3)) + coord_flip() +
  theme(axis.text.y =   element_text( vjust = 0.1, hjust = 1), panel.background = element_blank(), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())
ggsave(file="Ptger_Mouse_logfc_scale.pdf",  width=4, height=2, units="in")

ggplot(aa_zscore_mean, aes(metabolite, diet, fill= name)) + geom_tile() + 
  scale_fill_gradientn(name = "Log Counts", colors = pathwayColorsDiff, limits = c(0,7.5)) + coord_flip() +
  theme(axis.text.y =   element_text( vjust = 0.1, hjust = 1), panel.background = element_blank(), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())
ggsave(file="Ptger_Mouse_logfc_scale.pdf",  width=4, height=2, units="in")



