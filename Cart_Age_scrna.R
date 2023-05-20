#Vyom Shah
#February 2023
#CART Senscence project - Old vs Young Mice single-cell RNA seq analysis
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

# Load in the Data
Old1 <- Read10X(data.dir = "/Users/vyom/data/Old_vs_Young/SB26/filtered_feature_bc_matrix/")
Old2 <- Read10X(data.dir = "/Users/vyom/data/Old_vs_Young/SB27/filtered_feature_bc_matrix/")
Old3 <- Read10X(data.dir = "/Users/vyom/data/Old_vs_Young/SB31/filtered_feature_bc_matrix/")

Old1 <- CreateSeuratObject(Old1, project = "Old1")
Old2 <- CreateSeuratObject(Old2, project = "Old2")
Old3 <- CreateSeuratObject(Old3, project = "Old3")

Young1 <- Read10X_h5(file = "/Users/vyom/data/SB10_SB11_results/Individual_Samples/SB10/Beyaz_10_10xgex_control/filtered_feature_bc_matrix.h5")
Young2 <- Read10X_h5(file = "/Users/vyom/data/SB10_SB11_results/Individual_Samples/SB11/Beyaz_11_10xgex_control/filtered_feature_bc_matrix.h5")

Young1 <- CreateSeuratObject(Young1, project = "Young1")
Young2 <- CreateSeuratObject(Young2, project = "Young2")

d <- merge(Old1, y = c(Old2,Old3,Young1,Young2), add.cell.ids = c('Old1','Old2','Old3','Young1','Young2'), project = "Age")
d
d[["percent.mt"]] <- PercentageFeatureSet(d, pattern = "mt-")
VlnPlot(d, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),pt.size = 0 ,ncol = 3)
d <- subset(d, subset = nCount_RNA > 1000 & nCount_RNA < 20000 & nFeature_RNA > 100 & nFeature_RNA < 8000 & percent.mt < 15)
d
{
d$orig.ident
table(d$orig.ident)
s.genes <- cc.genes$s.genes %>% convert_human_to_mouse_symbols()
s.genes <- s.genes[!is.na(s.genes)]
g2m.genes <- cc.genes$g2m.genes %>% convert_human_to_mouse_symbols()
g2m.genes <- g2m.genes[!is.na(g2m.genes)]
d <- NormalizeData(d)
d <- FindVariableFeatures(d, selection.method = "vst")
d <- ScaleData(d, features = rownames(d))
d  <- RunPCA(d , features = VariableFeatures(d), ndims.print = 6:10, nfeatures.print = 10)
d <- CellCycleScoring(d, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# view cell cycle scores and phase assignments
head(d[[]])
d$CC.Difference <- d$S.Score - d$G2M.Score
age_obj <- SCTransform(
  d,
  vars.to.regress = c("S.Score", "G2M.Score", 'orig.ident'),
  verbose=TRUE
)
gc()
age_obj <- RunPCA(age_obj)
age_obj <- age_obj %>% 
  RunHarmony("orig.ident", assay.use="SCT", plot_convergence = TRUE)
rm(Data.list, d)

# UMAP and clustering with harmonized PCs
age_obj <- RunUMAP(age_obj, dims = 1:50, reduction='harmony', n.neighbors = 50, n.epochs = 500)
age_obj <- FindNeighbors(age_obj, reduction='harmony', dims = 1:50, k.param = 30, compute.SNN = TRUE)
age_obj <- FindClusters(age_obj, resolution = 1)
DimPlot(age_obj, reduction = "umap", label = TRUE)

DimPlot(age_obj, reduction = "umap", label = TRUE)
}
{
  Data.list <- SplitObject(d, split.by = "ident")
  Data.list <- Data.list[c('Old1','Old2','Old3','Young1','Young2')]
  for (i in 1:length(Data.list)) {
    
    Data.list[[i]] <- SCTransform(Data.list[[i]], verbose = FALSE)
  }
  
  # Normilization
  #select highly variable genes 
  Data.features <- SelectIntegrationFeatures(object.list = Data.list, nfeatures = 4000)
  options (future.globals.maxSize = 4000 * 1024^8)
  Data.list <- PrepSCTIntegration(object.list = Data.list, anchor.features = Data.features, 
                                  verbose = FALSE)
  Data.anchors <- FindIntegrationAnchors(object.list = Data.list, normalization.method = "SCT", 
                                         anchor.features = Data.features, verbose = FALSE)
  rm(Data.list, d)
  age_obj <- IntegrateData(anchorset = Data.anchors, normalization.method = "SCT", 
                           verbose = TRUE)
  rm(Data.anchors)
  rm(Old1,Old2,Old3,Young1,Young2)
}
# Visulization and Clustering
age_obj <- RunPCA(age_obj)

VizDimLoadings(age_obj, dims = 1:2, reduction = "pca")

DimPlot(age_obj, reduction = "pca")
ElbowPlot(age_obj, ndims = 100, reduction = "pca")




levels(factor(age_obj@meta.data$orig.ident))
age_obj[["Age"]] <- Idents(age_obj)
new.cluster.ids <- c("Old", "Old", "Old", "Young", 'Young')
names(new.cluster.ids) <- levels(age_obj)
age_obj <- RenameIdents(age_obj, new.cluster.ids)
age_obj[["Age"]] <- Idents(age_obj)

age_obj <- RunUMAP(age_obj, dims = 1:15, n.neighbors = 5, n.epochs = 500)
age_obj <- FindNeighbors(age_obj, dims = 1:50)
age_obj <- FindClusters(age_obj, resolution = .5)
DimPlot(age_obj, reduction = "umap", label = TRUE, group.by = 'orig.ident')
DimPlot(age_obj, reduction = "umap", label = TRUE)


DotPlot_Sig <- c("Lgr5","Ascl2","Olfm4","Gkn3","S100a6","Ly6a","Anxa3", "Areg","Tubb5","Syce2","Stmn1","Fbxo5",'Cenpa','Ccna2','Ube2c','Cdkn3',"Apoa1","Apoa4","Fabp1","Adh6a",'Tmigd1', 'Fabp6', 'Slc51b', 'Slc51a', "Chgb","Tac1","Tph1","Neurog3", "Muc2","Fcgbp","Atoh1","Agr2","Lyz1","Defa17","Defa24","Ang4","Pou2f3","Avil","Tuba1a","Adh1",'Ptprc','Cd8a','Cd4','Ighm') 
DotPlot(age_obj, features = DotPlot_Sig,dot.scale = 100, scale.max= 1000, assay = 'SCT') + labs(y= "Cell type", x="") + scale_colour_distiller( palette ="RdYlBu") + scale_size(range = c(0, 5)) +
  theme(legend.key.size = unit(.1, "in"), legend.text= element_text(size=4), legend.title = element_text(size=4),axis.line = element_line(size = .3), axis.ticks = element_line(size = .3), text = element_text(size=4), axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust= .01), axis.text.y = element_text(size = 10), axis.title.x = element_blank(), axis.title.y = element_blank())
ggsave(file="Cart_Cluster_check.pdf", width=4, height=2)


