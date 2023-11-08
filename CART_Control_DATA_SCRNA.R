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
OF2 <- Read10X(data.dir = "/Users/vyom/data/Cart_control_scrna/Amor_CAV09_OF2/outs/filtered_feature_bc_matrix/")
OM <- Read10X(data.dir = "/Users/vyom/data/Cart_control_scrna/Amor_CAV09_OM/outs/filtered_feature_bc_matrix/")
OF1 <- Read10X(data.dir = "/Users/vyom/data/Cart_control_scrna/Amor_CAV09_OF1/outs/filtered_feature_bc_matrix/")
YF2 <- Read10X(data.dir = "/Users/vyom/data/Cart_control_scrna/Amor_CAV09_YF2/outs/filtered_feature_bc_matrix/")
YF1 <- Read10X(data.dir = "/Users/vyom/data/Cart_control_scrna/Amor_CAV09_YF1/outs/filtered_feature_bc_matrix/")
YM <- Read10X(data.dir = "/Users/vyom/data/Cart_control_scrna/Amor_CAV09_YM/outs/filtered_feature_bc_matrix/")



OF2 <- CreateSeuratObject(OF2, project = "OF2")
OM <- CreateSeuratObject(OM, project = "OM")
OF1 <- CreateSeuratObject(OF1, project = "OF1")
YF2 <- CreateSeuratObject(YF2, project = "YF2")
YF1 <- CreateSeuratObject(YF1, project = "YF1")
YM <- CreateSeuratObject(YM, project = "YM")

d <- merge(OF2, y = c(OM,OF1,YF2,YF1,YM), add.cell.ids = c("OF2", "OM","OF1","YF2","YF1","YM"), project = "Cart_control")
d
d[["percent.mt"]] <- PercentageFeatureSet(d, pattern = "mt-")
VlnPlot(d, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),pt.size = 0,ncol = 3)
d <- subset(d, subset = nCount_RNA > 750 & nCount_RNA < 30000 & nFeature_RNA > 1000 & nFeature_RNA < 6000 & percent.mt < 15)
d
d$orig.ident
Data.list <- SplitObject(d, split.by = "ident")
Data.list <- Data.list[c("OF2", "OM","OF1","YF2","YF1","YM")]
for (i in 1:length(Data.list)) {
  
  Data.list[[i]] <- SCTransform(Data.list[[i]], verbose = FALSE)
}

# Normilization
#select highly variable genes 
Data.features <- SelectIntegrationFeatures(object.list = Data.list, nfeatures = 2500)
options (future.globals.maxSize = 4000 * 1024^8)
Data.list <- PrepSCTIntegration(object.list = Data.list, anchor.features = Data.features, 
                                verbose = FALSE)
Data.anchors <- FindIntegrationAnchors(object.list = Data.list, normalization.method = "SCT", 
                                       anchor.features = Data.features, verbose = FALSE)
rm(Data.list, d)
cart_control <- IntegrateData(anchorset = Data.anchors, normalization.method = "SCT", 
                          verbose = TRUE)
rm(Data.anchors)
rm(uPAR_F_young, uPAR_M_young, uPAR_F_old, uPAR_M_old,UT_F_young, UT_young,UT_F_old, UT_M_old)
# Visulization and Clustering
cart_control <- RunPCA(cart_control)

VizDimLoadings(cart_control, dims = 1:2, reduction = "pca")

DimPlot(cart_control, reduction = "pca")
ElbowPlot(cart_control, ndims = 100, reduction = "pca")

cart_control <- RunUMAP(cart_control, dims = 1:15, n.neighbors = 5, n.epochs = 500)



Idents(cart_control) <- factor(cart_control$orig.ident)

levels(factor(cart_control@meta.data$orig.ident))
cart_control[["Sex"]] <- Idents(cart_control)
new.cluster.ids <- c("Female", "Female","Male", "Female", "Female","Male")
names(new.cluster.ids) <- levels(cart_control)
cart_control <- RenameIdents(cart_control, new.cluster.ids)
cart_control[["Sex"]] <- Idents(cart_control)

Idents(cart_control) <- factor(cart_control$orig.ident)

levels(factor(cart_control@meta.data$orig.ident))
cart_control[["Age"]] <- Idents(cart_control)
new.cluster.ids <- c('Old','Old','Old','Young','Young','Young')
names(new.cluster.ids) <- levels(cart_control)
cart_control <- RenameIdents(cart_control, new.cluster.ids)
cart_control[["Age"]] <- Idents(cart_control)

levels(as.factor(cart_control$Sex))
levels(as.factor(cart_control$Age))
levels(as.factor(cart_control$orig.ident))
meta_check <- cart_control@meta.data
cart_control <- FindNeighbors(cart_control, dims = 1:30)
cart_control <- FindClusters(cart_control, resolution = 1)
DimPlot(cart_control, reduction = "umap", label = TRUE) + theme(legend.position = 0)

DimPlot(cart_control, reduction = "umap", group.by = 'Age')
DimPlot(cart_control, reduction = "umap", group.by = 'Sex')


DotPlot_Sig <- c('Ptprc',"Lgr5","Ascl2","Olfm4","Gkn3","S100a6","Ly6a","Anxa3", "Areg","Tubb5","Syce2","Stmn1","Fbxo5",'Cenpa','Ccna2','Ube2c','Cdkn3',"Apoa1","Apoa4","Fabp1","Adh6a",'Tmigd1', 'Fabp6', 'Slc51b', 'Slc51a', "Chgb","Tac1","Tph1","Neurog3", "Muc2","Fcgbp","Atoh1","Agr2","Lyz1","Defa17","Defa24","Ang4","Pou2f3","Avil","Tuba1a","Adh1",'Cd8a','Cd4','Ighm') 
DotPlot(cart_control, features = DotPlot_Sig, dot.scale = 10, scale.max= 100, assay = 'SCT') + labs(y= "Cell type", x="") + scale_colour_distiller( palette ="RdYlBu") + scale_size(range = c(0, 1.3)) +
  theme(legend.key.size = unit(.1, "in"), legend.text= element_text(size=4), legend.title = element_text(size=4),axis.line = element_line(size = .3), axis.ticks = element_line(size = .3), text = element_text(size=4), axis.text.x = element_text(size = 4, angle = 90, hjust = 1, vjust= .01), axis.text.y = element_text(size = 4), axis.title.x = element_blank(), axis.title.y = element_blank())

'Stem', 'Transit Amplifying', 'Enterocyte Progenitor','Enterocyte (Proximal)','Enterocyte (Distal)', 'Enteroendocrine', 'Goblet', 'Paneth','Tuft', 'T Cell', 'B Cell', 'Myeloid'


cart_control <- FindClusters(cart_control, resolution = 1)
new.cluster.ids <- c('Enterocyte Progenitor', 'Transit Amplifying', 'Transit Amplifying', 'Enterocyte (Proximal)', 'Stem', 'Transit Amplifying', 'Transit Amplifying', 'Goblet', 'Enterocyte (Proximal)', 'Stem', 'Enterocyte (Proximal)', 'Enterocyte (Proximal)', 'Enterocyte Progenitor', 'Goblet', 'Stem', 'Enterocyte (Proximal)', 'Enterocyte (Distal)', 'Tuft', 'Enterocyte (Distal)', 'Enterocyte Progenitor', '', 'Enterocyte (Distal)', 'Enterocyte (Proximal)', 'Enteroendocrine', 'Goblet', 'Goblet', 'Paneth', '', 'Paneth', 'Goblet')  
cart_control[["Cell_Type"]] <- Idents(cart_control)
names(new.cluster.ids) <- levels(cart_control)
cart_control <- RenameIdents(cart_control, new.cluster.ids)
cart_control[["Cell_Type"]] <- Idents(cart_control)
DimPlot(cart_control, reduction = "umap", group.by= 'Cell_Type')


Idents(cart_control) <- cart_control$Cell_Type
Cart_immune <- subset(cart_control, idents = '')

Cart_immune <- RunPCA(Cart_immune, verbose = FALSE)

Cart_immune <- RunUMAP(Cart_immune, dims = 1:50, n.neighbors = 5, n.epochs = 500)

Cart_immune <- FindNeighbors(Cart_immune, dims = 1:50)
Cart_immune <- FindClusters(Cart_immune, resolution = 1)
DimPlot(Cart_immune, reduction = "umap", label = TRUE)


DotPlot_Sig <- unique(c('Ptprc',"S100a6","Ly6a","Anxa3", "Areg",'Cd3g','Cd3e','Cd8a','Cd4','Trac','Tcrg-C1','Lag3','Pdcd1','Havcr2','Tox','Tcf7','Gzmb','Tbx21','Ifng','Gata3','Il4','Rorc','Il17a','Il17f','Foxp3','Il10','Il2rb','Il2ra','Klrd1','Cd19','Ighm','Ighg1','Cd74','Ciita','Nrc1','Klre1','Itgam','Itgax','H2-Eb1','H2-Ab1','Arg1','Mrc1','Tgfbi','Ccr2','Vegfa','Prdx1','Clec4d','Ccl5','Cd83','Ccr7','Fcn1','Msrb1','Ly6g','Col3a1','Sparc'))
DotPlot(Cart_immune, features = DotPlot_Sig, dot.scale = 10, scale.max= 100, assay = 'SCT') + labs(y= "Cell type", x="") + scale_colour_distiller( palette ="RdYlBu") + scale_size(range = c(0, 1.3)) +theme(axis.text.x = element_text(size = 7, angle = 90, hjust = 1, vjust= .01))

new.cluster.ids <- c('Goblet','T Cell','T Cell','T Cell','T Cell','Myeloid','B Cell','T Cell')

Cart_immune[["Cell_Type"]] <- Idents(Cart_immune)
names(new.cluster.ids) <- levels(Cart_immune)
Cart_immune <- RenameIdents(Cart_immune, new.cluster.ids)
Cart_immune[["Cell_Type"]] <- Idents(Cart_immune)
DimPlot(Cart_immune, reduction = "umap", label = TRUE)

cart_control$Cell_type <- as.character(cart_control$Cell_Type)
cart_control$Cell_type[WhichCells(Cart_immune)] <- paste(Idents(Cart_immune))
Idents(cart_control) <- cart_control$Cell_type
DimPlot(cart_control, reduction = "umap", group.by= 'Cell_type')

my_levels <- c('Stem', 'Transit Amplifying', 'Enterocyte Progenitor','Enterocyte (Proximal)','Enterocyte (Distal)', 'Enteroendocrine', 'Goblet', 'Paneth','Tuft', 'T Cell', 'B Cell', 'Myeloid')
cart_control$Cell_Type <- factor(x = cart_control$Cell_type, levels = my_levels)
DimPlot(cart_control, reduction = "umap", group.by= 'Cell_Type')
ggsave(file = paste0('Cart_umap.pdf'), width=7, height=5, units="in")


DefaultAssay(cart_control) <- "RNA"
all.genes <- rownames(cart_control)
cart_control <- NormalizeData(object = cart_control, normalization.method = "LogNormalize", assay = "RNA")
cart_control <- ScaleData(cart_control, features = all.genes)
cart_control <- RunPCA(cart_control, features = all.genes)
cart_control <- magic(cart_control)

#new
senescence <- c('Psap','Timp2','Itm2b','Lgmn','Igf1','Arl6ip1','Ckb','Lsmem1','Igfbp4','Gsn','Zbtb20','Ccl2','Ccl7','Trio','Tulp4','Cxcl1')
#old
senescence <- c('Acvr1b','Ang','Angpt1','Angptl4','Areg','Axl','Bex3','Bmp2','Bmp6','C3','Ccl1','Ccl2','Ccl20','Ccl24','Ccl26','Ccl3','Ccl4','Ccl5','Ccl7','Ccl8','Cd55','Cd9','Csf1','Csf2','Csf2rb','Cst10','Ctnnb1','Ctsb','Cxcl1','Cxcl10','Cxcl12','Cxcl16','Cxcl2','Cxcl3','Cxcr2','Dkk1','Edn1','Egf','Egfr','Ereg','Esm1','Ets2','Fas','Fgf1','Fgf2','Fgf7','Gdf15','Gem','Gmfg','Hgf','Hmgb1','Icam1','Icam5','Igf1','Igfbp1','Igfbp2','Igfbp3','Igfbp4','Igfbp5','Igfbp6','Igfbp7','Il10','Il13','Il15','Il18','Il1a','Il1b','Il2','Il6','Il6st','Il7','Inha','Iqgap2','Itga2','Itpka','Jun','Kitl','Lcp1','Mif','Mmp13','Mmp10','Mmp12','Mmp13','Mmp14','Mmp2','Mmp3','Mmp9','Nap1l4','Nrg1','Pappa','Pecam1','Pgf','Pigf','Plat','Plau','Plaur','Ptbp1','Ptger2','Ptges','Rps6ka5','Scamp4','Selplg','Sema3f','Serpinb3a','Serpine1','Serpine2','Spp1','Spx','Timp2','Tnf','Tnfrsf11b','Tnfrsf1a','Tnfrsf1b','Tubgcp2','Vegfa','Vegfc','Vgf','Wnt16','Wnt2')
#senmayo
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

All_Genes <- cart_control@assays$RNA@data@Dimnames[[1]]
senescence <- intersect(All_Genes, senescence)
senescence <- intersect(All_Genes, sen_mayo)
mean.exp <- zscore(colMeans(x = cart_control@assays$MAGIC_RNA@data[senescence, ], na.rm = TRUE), dist ='norm')
if (all(names(x = mean.exp) == rownames(x = cart_control@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  cart_control@meta.data$senescence <- mean.exp
}

Idents(cart_control) <- cart_control$Cell_Type
max(cart_control$senescence)
FeaturePlot(object = cart_control, features = 'senescence', pt.size = .001, label = TRUE, label.size = 4, repel = TRUE) +
  theme(plot.title = element_blank(), text = element_text(size=6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), axis.title.x = element_blank(), axis.title.y = element_blank()) +
  scale_color_gradientn(colors=brewer.pal(n = 11, name = "YlGnBu"), breaks = c(.0, .1,.2,.3),limits = c(.0,.3), trans = scales::boxcox_trans(.9))
ggsave(file = 'CART_control_Scenscence_BuGn.pdf', width=5, height=5, units="in")

All <- VlnPlot(cart_control, split.by = "Age",group.by = 'Cell_Type', features = 'senescence', pt.size = 0, assay = "RNA",  cols = c('#1b9e77','#d95f02'), log = TRUE, split.plot = TRUE) + 
  theme(legend.position = 'none') + 
  geom_boxplot(width=0.3,outlier.shape = NA, coef = 0, lwd= .2) + xlab('') + 
  theme(text = element_text(size=9), axis.text.x = element_text(size = 9), axis.text.y = element_text(size =9)) 
All$layers[[1]]$aes_params$size = .15
All
cart_control$Age

DotPlot_Sig <- c("Lgr5","Ascl2","Olfm4","Gkn3","Tubb5","Syce2","Stmn1","Fbxo5",'Cenpa','Ccna2','Ube2c','Cdkn3',"Apoa1","Apoa4","Fabp1","Adh6a",'Tmigd1', 'Fabp6', 'Slc51b', 'Slc51a', "Chgb","Tac1","Tph1","Neurog3", "Muc2","Fcgbp","Atoh1","Agr2","Lyz1","Defa17","Defa24","Ang4","Pou2f3","Avil","Tuba1a","Adh1",'Ptprc','Cd8a','Cd4','Cd19') 
DotPlot(cart_control, features = DotPlot_Sig,group.by = 'Cell_Type',dot.scale = 10, scale.max= 100, assay = 'SCT') + labs(y= "Cell type", x="") + scale_colour_distiller( palette ="RdYlBu") + scale_size(range = c(0, 1.3)) +
  theme(legend.key.size = unit(.1, "in"), legend.text= element_text(size=4), legend.title = element_text(size=4),axis.line = element_line(size = .3), axis.ticks = element_line(size = .3), text = element_text(size=4), axis.text.x = element_text(size = 4, angle = 90, hjust = 1, vjust= .01), axis.text.y = element_text(size = 4), axis.title.x = element_blank(), axis.title.y = element_blank())
ggsave(file="Cart_control_Cluster_check.pdf", width=4, height=2)

Prop_table<- prop.table(x = table(cart_control$Cell_Type, cart_control$orig.ident), margin = 2)
Prop_Table <- as.data.frame(Prop_table, row.names = NULL, optional = FALSE,
                            make.names = TRUE, stringsAsFactors = default.stringsAsFactors())
levels(as.factor(cart_control$orig.ident))
Prop_Table1 <- Prop_Table%>%mutate(Var2=recode(Var2,OF1="Old",OF2="Old",OM="Old",YF1="Young", YF2="Young", YM="Young"))

celltype_sample_norm3 = summarySE(Prop_Table1, measurevar="Freq", groupvars=c("Var1","Var2"))

my_levels <- c('Stem', 'Transit Amplifying', 'Enterocyte Progenitor','Enterocyte (Proximal)','Enterocyte (Distal)', 'Enteroendocrine', 'Goblet', 'Paneth','Tuft', 'T Cell', 'B Cell', 'Myeloid')
celltype_sample_norm3$Var1 = factor(celltype_sample_norm3$Var1, levels = my_levels)
celltype_sample_norm3$Var2 = factor(celltype_sample_norm3$Var2, levels = c('Young','Old'))

ggplot(celltype_sample_norm3, aes(x = Var1, y = Freq, fill = Var2)) + geom_bar(stat="identity", position="dodge") + 
  geom_errorbar(aes(ymin=Freq-se, ymax=Freq+se),position=position_dodge(width = 0.85),width=0.3, size=0.25) + theme_vyom + 
  scale_fill_manual(values = colorsType) + expand_limits(y = c(0)) + 
  theme( axis.text.x = element_text(angle = 45, size =  6, vjust = 1, hjust = 1), panel.background = element_rect(fill = "white", colour = "Black"), axis.text.y = element_text(size = 6), axis.title.x = element_text(size = 7), axis.title.y = element_text(size = 7)) + 
  xlab("Cell Type") + ylab("Fraction of Cells") + labs(fill = "Age") + scale_y_continuous(expand = expansion(mult = c(0, .1)))
ggsave( "CART_CONTROL_frac_of_cells.pdf", width=3.75, height=3, units="in")

colorsType = c(
  Young = "#7CA1CC",
  Old = "#FF4902"
)

All_Genes <- cart_control@assays$RNA@data@Dimnames[[1]]
Inflammatory_score_new <- c('Ccl4','Ccl5','Ccl9','Ccl20','Cd74','Cox4I1','Cox8A','Ctsb','Ctsd','H2-Ac','H2-Ab1','H2-Q7','Ifi47','Ifit1Bl1','Ifit2','Ifit3','Ifitm3','Irf7','Isg15','Stat2','Stat3','Tap1','Ucp2')
Inflammatory_Response <- intersect(All_Genes, Inflammatory_score_new)
MHC_score <- c('H2-D1','H2-K','H2-Q1','H2-Q2','H2-T23','H2-T3','Cd74','H2-Ac','H2-Ab1','H2-Q7')
MHC_score <- intersect(All_Genes, MHC_score)

mean.exp <- zscore(colMeans(x = cart_control@assays$MAGIC_RNA@data[Inflammatory_Response, ], na.rm = TRUE), dist = 'norm')
if (all(names(x = mean.exp) == rownames(x = cart_control@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  cart_control@meta.data$Inflammatory_Response <- mean.exp
}

mean.exp <- zscore(colMeans(x = cart_control@assays$MAGIC_RNA@data[MHC_score, ], na.rm = TRUE), dist = 'norm')
if (all(names(x = mean.exp) == rownames(x = cart_control@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  cart_control@meta.data$MHC_score <- mean.exp
}

cart_control$
i = 'Inflammatory_Response'
i = "MHC_score"
i = "MHC_score"
for(i in gene.list) {
  selected_cells <- names(cart_control$Cell_Type)
  vln_data <- FetchData(cart_control,
                        vars = c(i,"Age", "Cell_Type"),
                        cells = selected_cells,
                        slot = "data")
  vln_data <- melt(vln_data)
  vln_data <- vln_data[c('Age', 'Cell_Type', 'value')]
  
  #  my_comparisons <- list(c('WT_Control, WT_Arasco'),c('WT_Control, WT_Rev'))
  All <- VlnPlot(cart_control,group.by = 'Cell_Type', split.by = "Age", features = i, pt.size = 0, assay = "MAGIC_RNA", cols = c('#7CA1CC','#FF4902'), log = FALSE, split.plot = TRUE) + 
    theme(legend.position = 'none')  +
    theme(axis.line = element_line(size = .3),axis.title.x = element_blank(), axis.ticks = element_line(size = .3), text = element_text(size=7), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6))  +
    stat_compare_means(data= vln_data, aes(x = Cell_Type, y = value, fill = Age), hide.ns = TRUE, method = "t.test",  label = "p.signif", size = 1)  #+ scale_y_continuous(expand = expansion(mult = c(0, .1)))
  dodge <- position_dodge(width = .9)
  All <- All + geom_boxplot(width=0.3,outlier.shape = NA, coef = 0, lwd=.2, position = dodge) +
    stat_compare_means(data= vln_data, aes(x = Cell_Type, y = value, fill = Age), hide.ns = TRUE, method = "t.test",  label = "p.signif", size = 1)  
  All$layers[[1]]$aes_params$size = .15
  ggsave(file = paste0('Cart_Control_', i, '.pdf'), plot=All, width=2.5, height=2.5, units="in")
}
cart_control$Sen_cells
for(i in gene.list) {
  selected_cells <- names(cart_control$Cell_Type)
  vln_data <- FetchData(cart_control,
                        vars = c(i,"Age", "Sen_cells"),
                        cells = selected_cells,
                        slot = "data")
  vln_data <- melt(vln_data)
  vln_data <- vln_data[c('Age', 'Sen_cells', 'value')]
  
  #  my_comparisons <- list(c('WT_Control, WT_Arasco'),c('WT_Control, WT_Rev'))
  All <- VlnPlot(cart_control,group.by = 'Sen_cells', split.by = "Age", features = i, pt.size = 0, assay = "MAGIC_RNA", cols = c('#7CA1CC','#FF4902'), log = FALSE, split.plot = TRUE) + 
    theme(legend.position = 'none')  +
    theme(axis.line = element_line(size = .3),axis.title.x = element_blank(), axis.ticks = element_line(size = .3), text = element_text(size=7), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6))
  dodge <- position_dodge(width = .9)
  All <- All + geom_boxplot(width=0.3,outlier.shape = NA, coef = 0, lwd=.2, position = dodge) +
    stat_compare_means(data= vln_data, aes(x = Sen_cells, y = value, fill = Age), hide.ns = TRUE, method = "t.test", size = 2)  
  All$layers[[1]]$aes_params$size = .15
  All
  ggsave(file = paste0('Cart_Control_senescent_cells_', i, '.pdf'), plot=All, width=2.5, height=2.5, units="in")
}

IscI <- c('Pdgfa', 'Edn1', 'Rgcc', 'Lgr5', 'Cyp2e1', 'Jun', 'Agr3', 'Rnf43', 'Sorbs2', 'Fstl1', 'Filip1l', 'Nrn1', 'Sord', 'Nr1d2', 'Efna1', 'Rnf32', 'Fam13a', 'Arhgef26', 'Sypl', 'Tmem171', 'Slc12a2', 'Gm8979', 'Fgf1', 'Ogt', 'Fhdc1', 'Oat', 'Arid5b', 'Prss23', 'Snx10', 'Chp1', 'Sh3rf1', 'Lbh', 'Myo9a', 'Soat1', 'Npc2', 'Trim2', 'Eif3e', 'Slc14a1', 'Ivns1abp', 'Tprkb', 'Igbp1', 'Nsa2', 'Olfm4', 'Ces1d', 'Pik3r1', 'Kcnq1ot1', 'Tspan12', 'Aqp1', 'Eif4b', 'Prosc', 'Mtus1', 'Cox7a2l', 'Nme7')
IscII <- c('Rps27', 'Rpl26', 'Angptl4', 'H2-Ab1', 'Cd320', 'Ifitm3', 'Scn2b', 'Rpl34-ps1', 'Arglu1', 'Clca4b', 'Gm6654', 'Rn4.5s', 'Cd74', 'Zyx', 'Acot2', 'Nfkbia', 'H2-Eb1', 'Zcchc7', 'Rps15a-ps4', 'Dus2l', 'Rangrf', 'Gm11974', 'Myc', 'Rps28', 'Snx16', 'Per3', '4921508A21Rik', 'Car12', 'Clec2d', 'Zfp106', 'A230050P20Rik', 'Clca2', 'Sp140', '2410006H16Rik', 'Ddo', 'Rpl31-ps12', 'Nfic', 'Arhgef4', 'Zfp109', 'Pdrg1', 'Paics', 'Ndufa7', 'Tomm7', '1600029O15Rik', 'Rps19-ps3', 'C1ra', 'Gstm1', 'Snora43', 'Ankrd10', 'Relb', 'Slc15a2', 'Clca4', 'Rps17', 'Gm15772', 'Rps25', 'Rps18', 'Rpl37', 'Rps15a', 'Rps19', 'Rpl32', 'Gm10548', '6030458C11Rik', 'Rps13', 'Snhg8', 'Klhdc5', 'Ifitm2', 'Rps29', 'Rplp1', 'Rps12', 'Rpl36a', 'Rpl35', 'Rpl38', 'Dctd', 'Rps15a-ps6', 'Rpl37a', 'Rpl36', 'Urod', 'Shmt1', 'Rps26', '2210039B01Rik', 'Zbtb16', 'Noxa1', 'Rpl39', 'Rps23', 'Pex26', 'Cyba', 'Rps14', 'Shfm1', 'Gm12191', 'Rpl11', 'Mir703', 'Rps21', 'Rps10', 'Rpl35a', 'Rps24', 'Clec2g', 'Rpl31', 'Rpl23a', 'Rps16', 'Rpl18a') 
IscIII <- c('Hist1h2ao','Neil3','Top2a','Aurkb','Pbk','Ncaph','Hist2h3b','Cdca2','Ankle1',
            'Incenp','Fbxo5','Cdca5','Cenph','Prc1','Hist2h3c2','Cks1b','Cdca8','Bub1','Spag5',
            'Hist1h4d','Nusap1','Rhno1','Cdk1','Ccnf','Ndc80','Ncapd2','Cenpm','Esco2','Rad51ap1',
            'Hist1h1e','Psrc1','Haus5','Tubb5','Poc1a','Cit','Shcbp1','Mis18bp1','Asf1b','Cenpp',
            'Nrm','Nuf2','Trim59','Mis18a','Smc4','2810417H13Rik','Ska1','Kif20b','Kif20a','Rad51',
            'Cenpn','Mlf1ip','Haus8','Chek2','Rad54b','Rad18','Plk4','Kif18b','Ncapg','Eme1',
            'Kif15','Cdca4','Ttk','Ube2t','Kif11','Hist1h2ab','Mastl','Spc25','Tk1','Ctc1',
            'Kif2c','Hist1h1b','Tnfaip8l1','Rad54l','Fbln1','BC030867','Tuba1b','Ttf2','Rad51b',
            'Sgol1','Tinf2','Fancd2','Gtse1','Trip13','H2afx','H2-Q7','Ckap2l','Nup133','Ccdc34',
            'Cenpt','Ncapg2','Espl1','Cdc45','Aaas','Stil','Rfc5','Anln','Oip5','Ska3','Haus1',
            'Cyp39a1','Ect2','Smc2','2700094K13Rik','Zwilch','Hist1h2ad','Mns1','Suv39h1','Bub1b',
            'Ddx11','Hist1h2ae','Mxd3','Ezh2','Cenpk','Tyw1','Fanci','Hist1h2ak','6430706D22Rik',
            'Cmss1','Rsrc1','Cks2','Anapc15','Ncapd3','Sept10','Eri2','Rnf26','BC055324','Tmem194',
            'Sgol2','Rrm2','Ercc6l','Cep72','A730008H23Rik','Znhit3','Sass6','Dsn1','Ticrr','Mcm10',
            'Casc5','Cep57l1','Blm','Dbf4','AI450353','Fignl1','Prim2','Kif18a','Psmc3ip','Mad1l1',
            'Rbl1','Siah1b','Dnajc9','Melk','Cep110','Racgap1','Trim37','Pck2','Cenpc1','Nsl1','Ccdc77',
            'Stmn1','Brd8','Efcab11','Lmnb1','Exo1','Haus6','Smc1a','Med4','Rfc4','Haus4','Tacc3',
            'Nnt','4930579G24Rik','Hjurp','Traip','Pkmyt1','Vrk1','D030056L22Rik','Hirip3','Cenpl',
            'Rbmx2','Cenpq','Cdo1','Rfc3','Slc7a11','2700099C18Rik','Mgme1','Mki67','Kifc5b','Tyms',
            'Clspn','Nfyb','Rangap1','Tube1','Mis12','Pcbd2','Hist1h3b','1190002F15Rik','Cenpw',
            'Atad2','Hist1h2af','Cep128','4930523C07Rik','Pole','Hist1h4i','Cenpi','Haus3','Fancg',
            'Kntc1','Nelfe','Cep57','Cntrob','Spc24','Repin1','Prim1','Pask','Rbm15','Ccdc18',
            'Gmcl1','Nup37','Mms22l','Rrm1','Arhgap11a','Tmpo','Foxm1','Nmral1','Cpsf4','Bud13',
            'Miip','Sclt1','Exosc8','Iqgap3','Sun2','Xrcc3','Nsmce4a','Tmem107','Hist1h2bj','Gins3',
            'Gins4','Ppie','Cbx5','Naa38','9430015G10Rik','Chaf1a','Depdc1a','Pmf1','Cdk5rap2','Zranb3',
            'Xrcc2','Ipo9','Hmgb2','Brca2','E2f7','Wdr62','Kifc1','Spata24','Sephs1','Nudt1',
            'Cep192','Fen1','Rfwd3','Hist2h2ac','Asrgl1','Tiam1','Zfp367','Nup62','Clhc1','Troap',
            'Sin3a','Ube2s','Eefsec','Rad9a','Chchd6','Cep55','Mcph1','2310008H04Rik','Mybl2',
            'Mre11a','Rbbp8','Xrcc1','Hist1h3c','Cdc7','Parpbp','Pcnt','Nde1','Gmnn','Elof1','Brca1',
            'Kif22','Mnd1','Cenpf','Hnrnpul1','Itprip','Pola2','Med24','Aspm','Ino80e','Pms2',
            'Lrrc45','Krt15','Tcof1','Odf2','Lrr1','Mdc1','Palb2','Lage3','2810442I21Rik','Whsc1',
            'Ahctf1','Gtpbp10','Rasa3','Wdr8','Naa40','Ube2c','Xkr5','Cenpo','Nup214','Nucks1',
            'Parp2','Stra13','Tex10','Dclre1a','Pms1','C230052I12Rik','2810408I11Rik','Rpap1',
            'Topbp1','Dnmt1','Zfp41','Fxn','Ubap2','Dnph1','Bard1','E2f8','Fam84a','Polq',
            'Pold1','Ckap5','Mettl14','Hist1h2ai','Cox6b2','Atad5','5830418K08Rik','Ilf3',
            'Hyls1','Nup107','Fgfr1op','Tmco6','Hspa14','Spice1','2700029M09Rik','Kbtbd4','Pom121',
            'Phf6','Med16','Ephx1','Flywch1','2810428I15Rik','Dtymk','Nthl1','Ccp110','Fanca',
            'Topors','BC052040','Hist1h2bk','Rnaseh2b','Dgcr8','Nup160','Wdr31','Smarca5','G2e3',
            'Nos2','Hn1l','Sp1','Nup205','U2af2','Hist1h2ag','Ipo11','Gas2l3','Bcl2l12','Kif24',
            'Zfp1','Rad1','Vars','Skp2','Slc9a8','Nudt21','Mum1','Papd7','Lsm2','E4f1','Dek','Bckdk',
            'B3galtl','Zfp101','4932415G12Rik','Hist1h1a','Ncoa7','Dars2','Hist1h3f','Nt5dc2','Srrt',
            'Zfp828','Hist1h3e','Hist1h2bb','Slc20a2','Hemgn','Cse1l','Zdhhc15','1700063D05Rik',
            'Cdkn2c','Mtfr2','Nudc','Top3a','Pold2','Rtel1','Rcc2','Fam76b','Prdx4','Nfix',
            'Hist1h2bl','4930503L19Rik','Ogg1','Kif14','Anapc11','Suz12','Mtr','Notch2','Mbd4',
            'Ccdc15','Thada','Hmgn2','Pdik1l','Fam111a','Mcm8','Lrrc49','Fbxo48','Magoh','Pds5b',
            'Kat7','Polr2f','Gpr19','Rttn','Lonp1','Hist1h2bn','Cep250','Fbxl14','Hist1h3d','Asb3',
            'Tubgcp6','Cenpj','Ccdc14','Tbrg3','Pced1a','Prkd3','Tbc1d5','Eftud2','Cep135','Sae1',
            'Snrpd1','H2afv','Ssrp1','Diap3','Senp8','Stag1','4930558J18Rik','Depdc1b','BC053749',
            'Brip1','Fgr','Rps6ka6','Hist2h4','Alpk1','Ddx19a','Fam122b','Ecm2','Hist1h3i','Itga1',
            'Gen1','Zfp958','Hist1h2bg','Cep152','Dcp1a','Acsf2','Hist1h3a')
IscI <- intersect(All_Genes, IscI)
IscII <- intersect(All_Genes, IscII)
IscIII <- intersect(All_Genes, IscIII)

mean.exp <- zscore(colMeans(x = cart_control@assays$MAGIC_RNA@data[IscI, ], na.rm = TRUE), dist = 'norm')
if (all(names(x = mean.exp) == rownames(x = cart_control@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  cart_control@meta.data$IscI <- mean.exp
}
mean.exp <- zscore(colMeans(x = cart_control@assays$MAGIC_RNA@data[IscII, ], na.rm = TRUE), dist = 'norm')
if (all(names(x = mean.exp) == rownames(x = cart_control@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  cart_control@meta.data$IscII <- mean.exp
}
mean.exp <- zscore(colMeans(x = cart_control@assays$MAGIC_RNA@data[IscIII, ], na.rm = TRUE), dist = 'norm')
if (all(names(x = mean.exp) == rownames(x = cart_control@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  cart_control@meta.data$IscIII <- mean.exp
}

source("~/analysis/split_violin.R")
cart_control$Age
genes_score <- c('IscI','IscII','IscIII')
DefaultAssay(cart_control@assays$MAGIC_RNA) <- "MAGIC_RNA"
selected_cells <- names(cart_control$Cell_Type[cart_control$Cell_Type == c('Stem')])
vln_data <- FetchData(cart_control,
                      vars = c(genes_score,"Age"),
                      cells = selected_cells,
                      slot = "data")

vln_data <- melt(vln_data)
vln_data1 <- vln_data[,-3]
All <- ggplot(vln_data1, aes(x = variable, y = value, fill= Age)) + geom_split_violin(scale="width", trim = TRUE, size = .1) + theme_vyom + theme(legend.position = 'none') + 
  geom_boxplot(width=0.2,outlier.shape = NA, coef = 0, lwd=.2) + xlab('') + ylab('Expression Level') + scale_fill_manual(values= c('#7CA1CC' ,'#FF4902')) +
  theme(text = element_text(size=7), axis.text.x = element_text(size = 6, angle = 45, hjust = 1, vjust = 1), axis.text.y = element_text(size = 6), axis.title.y = element_text(size = 7)) +
  stat_compare_means( method = "wilcox.test",  label = "p.signif", size = 2) + scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
All
ggsave(file = 'Control_stem_ISC_differences.pdf', plot=All, width=2.5, height=2.5, units="in")

Idents(cart_control) <- cart_control$Cell_Type
cart_control[["celltype"]] <- Idents(cart_control)
new.cluster.ids <- c("Stem",  "Epithelial", "Epithelial", "Epithelial", "Epithelial", "Epithelial", "Epithelial", "Epithelial", "Epithelial", "Immune", "Immune", "Immune")
names(new.cluster.ids) <- levels(cart_control)
cart_control <- RenameIdents(cart_control, new.cluster.ids)
cart_control[["celltype"]] <- Idents(cart_control)

Idents(cart_control) <- cart_control$celltype
#Cart_epith <- subset(cart_control, idents = c('Stem','Epithelial'))

cart_control$senescence
i = 'senescence'
for(i in gene.list) {
  selected_cells <- names(cart_control$celltype)
  vln_data <- FetchData(cart_control,
                        vars = c(i,"Age", "celltype"),
                        cells = selected_cells,
                        slot = "data")
  vln_data <- melt(vln_data)
  vln_data <- vln_data[c('Age', 'celltype', 'value')]
  
  #  my_comparisons <- list(c('WT_Control, WT_Arasco'),c('WT_Control, WT_Rev'))
  All <- VlnPlot(cart_control,group.by = 'celltype', split.by = "Age", features = i, pt.size = 0, assay = "MAGIC_RNA", cols = c('#7CA1CC','#FF4902'), log = FALSE, split.plot = TRUE) + 
    theme(legend.position = 'none')  +
    theme(axis.line = element_line(size = .3),axis.title.x = element_blank(), axis.ticks = element_line(size = .3), text = element_text(size=7), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6))  +
    stat_compare_means(data= vln_data, aes(x = celltype, y = value, fill = Age), hide.ns = TRUE, method = "t.test",  label = "p.signif", size = 1)  #+ scale_y_continuous(expand = expansion(mult = c(0, .1)))
  dodge <- position_dodge(width = .9)
  All <- All + geom_boxplot(width=0.3,outlier.shape = NA, coef = 0, lwd=.2, position = dodge) +
    stat_compare_means(data= vln_data, aes(x = celltype, y = value, fill = Age), hide.ns = TRUE, method = "t.test",  label = "p.signif", size = 1)  
  All$layers[[1]]$aes_params$size = .15
  ggsave(file = paste0('Mouse_',i,'_control_vln.pdf'), plot=All, width=2.5, height=2.5, units="in")
}


Idents(Cart_epith) <- Cart_epith$Sen_cells
DE_treat <- FindMarkers(Cart_epith, ident.1 = "Senescent", ident.2 = "Non", test.use = "MAST", logfc.threshold = .1, min.pct = .05, assay = 'RNA')
EnhancedVolcano(DE_treat, lab = rownames(DE_treat), x = 'avg_log2FC', y = 'p_val_adj', title = 'uPAR vs UT', pCutoff = 10e-1, FCcutoff = 0.25, xlim = c(-1,1), ylim = c(0,320), subtitle = 'All Cells' )
DE_treat$gene_id <- convert_mouse_to_human_symbols(rownames(DE_treat))
library(fgsea)

res2 <- DE_treat %>% 
  dplyr::select(gene_id, avg_log2FC) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(gene_id) %>% 
  summarize(avg_log2FC=mean(avg_log2FC))
res2

ranks <- deframe(res2)
head(ranks, 20)
#c2.all.v2022.1.Hs.symbols.gmt - most comprehensive all
#c5.all.v2022.1.Hs.symbols.gmt - second most comprehensive
#h.all.v2022.1.Hs.symbols.gmt - overall summary

pathways.hallmark <- gmtPathways("Documents/c2.all.v2022.1.Hs.symbols.gmt")
pathways.hallmark %>% 
  head() %>% 
  lapply(head)

fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

# Show in a nice table:
fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -log2err) %>% 
  arrange(padj) %>% 
  DT::datatable()
fgseaResTidier <- fgseaResTidy[which(fgseaResTidy$padj < .05 ),]
fgseaResTidier %>% 
  dplyr::select(-leadingEdge, -ES, -log2err) %>% 
  arrange(padj) %>% 
  DT::datatable()
ggplot(fgseaResTidier, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=pval<0.01)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()
gsea_output <- as.data.frame(fgseaResTidier)
gsea_output$leadingEdge <- NULL
write.csv(gsea_output,'CarT_human_GSEA_sen_vs_non.csv')
colnames(gsea_output)
theTable$Position <- factor(theTable$Position, levels = c(...))
pathways <- c('FEVR_CTNNB1_TARGETS_UP', 'SANSOM_APC_TARGETS_DN', 'REACTOME_BIOLOGICAL_OXIDATIONS', 'KEGG_PPAR_SIGNALING_PATHWAY', 'REACTOME_METABOLISM_OF_LIPIDS', 'LEE_BMP2_TARGETS_DN', 'REACTOME_CELL_CYCLE', 'WANG_TUMOR_INVASIVENESS_UP', 'TIEN_INTESTINE_PROBIOTICS_24HR_DN', 'KEGG_SPLICEOSOME')
pathwayColors <- c( "#FDC926FF", "#FA9E3BFF", "#ED7953FF", "#D8576BFF", "#BD3786FF", "#9C179EFF", "#7301A8FF", "#47039FFF")
gsea_output1 <- gsea_output[order(-gsea_output$NES),]
gsea_output1 <- gsea_output1[ gsea_output1$pathway %in% pathways,]
gsea_output1$pathway <- factor(gsea_output1$pathway, levels = gsea_output1$pathway)
All <- ggplot(data = gsea_output1,aes(x=pathway,y=NES)) +  theme_vyom +
  theme(axis.text.y = element_text(size = 3) ,axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), panel.background = element_rect(fill = "white", colour = "Black", size = .75, linetype = 'solid')) +
  geom_segment( aes(x=pathway, xend=pathway, y=0, yend=NES, color = padj))+
  geom_point(color = 'black' ,size = 1.2 ) +
  geom_point(aes(color = padj) ,size = 1 ) + geom_vline(xintercept = 0) +
  scale_color_gradientn(colors=pathwayColors, breaks = c(.05, .01, .001, .00001),labels =  c('.05', '.01', '.001', '<.0001'),limits = c(0,.05), trans = scales::boxcox_trans(.25)) +
  coord_flip() + 
  scale_y_continuous(limits = c(-3,3.3),expand = expansion(mult = c(0, 0)), breaks = scales::breaks_extended(n = 3)) +
  labs(y= "Normalized Enrichment Score", x="Pathway") 
All + geom_vline(xintercept = 0)
ggsave(file = paste0('Fig 1 H_GSEA.pdf'), plot=All, width=5, height=3, units="in")

cell_type_list <- levels(factor(cart_control$Cell_Type))
for (i in cell_type_list){
  o_cells_de<- rownames(cart_control@meta.data[cart_control$Cell_Type == c(i) & cart_control$Age == c('Old'),] )
  y_cells_de<- rownames(cart_control@meta.data[cart_control$Cell_Type == c(i) & cart_control$Age == c('Young'),] )
  
  DE_subset <- FindMarkers(cart_control, ident.1 = o_cells_de, ident.2 = y_cells_de,  test.use = "MAST", logfc.threshold = .1, min.pct = .01, assay = 'SCT')
  write.csv(DE_subset, paste0(i, '_old_v_young_control.csv'))
}

#determine senescent cells new
sen_mayo
cart_control$sen_cells <- 'non'
cart_control <- AddModuleScore(object = cart_control, features = list(sen_mayo), name = 'sen_mayo')

# module score distribution
Senescence_module <- as.data.frame(cart_control$sen_mayo1)
modulescores <- Senescence_module %>%
  rownames_to_column(var="id") %>%
  pivot_longer(-id, names_to="celltype", values_to="score")


p <- ggplot(modulescores)
#p <- ggplot(onescore)
inflection_plot <- p + geom_point(aes(x=fct_inorder(id), y=sort(score)))+ scale_y_continuous(breaks = seq(-.1, .3, by = .01)) +
  facet_wrap(~celltype) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
inflection_plot

sen_cells <-WhichCells(object = cart_control, expression = sen_mayo1 > .05)
DimPlot(cart_control, label=T,group.by = 'celltype', cells.highlight= list(sen_cells),  cols.highlight = c("darkblue"),cols= "grey")

cart_control$sen_cells[sen_cells] <- paste('senescent')

Prop_table<- prop.table(x = table( cart_control$Age, cart_control$sen_cells), margin = 2)
Prop_Table <- as.data.frame(Prop_table, row.names = NULL, optional = FALSE,
                            make.names = TRUE, stringsAsFactors = default.stringsAsFactors())
Prop_Table <- Prop_Table[Prop_Table$Var2 == 'senescent',]
melt(Prop_Table)
ggplot(Prop_Table, aes(fill=Var1, y=Freq, x=Var2)) + 
  geom_bar(position="fill", stat="identity")+theme_vyom +scale_fill_manual(values = colorsType)
ggsave( "CART_CONTROL_senescent_plot.pdf", width=2.5, height=3, units="in")


levels(as.factor(cart_control$orig.ident))
Prop_Table1 <- Prop_Table%>%mutate(Var2=recode(Var2,OF1="Old",OF2="Old",OM="Old",YF1="Young", YF2="Young", YM="Young"))

celltype_sample_norm3 = summarySE(Prop_Table1, measurevar="Freq", groupvars=c("Var1","Var2"))

my_levels <- c('Stem', 'Transit Amplifying', 'Enterocyte Progenitor','Enterocyte (Proximal)','Enterocyte (Distal)', 'Enteroendocrine', 'Goblet', 'Paneth','Tuft', 'T Cell', 'B Cell', 'Myeloid')
celltype_sample_norm3$Var1 = factor(celltype_sample_norm3$Var1, levels = my_levels)
celltype_sample_norm3$Var2 = factor(celltype_sample_norm3$Var2, levels = c('Young','Old'))

ggplot(celltype_sample_norm3, aes(x = Var1, y = Freq, fill = Var2)) + geom_bar(stat="identity", position="dodge") + 
  geom_errorbar(aes(ymin=Freq-se, ymax=Freq+se),position=position_dodge(width = 0.85),width=0.3, size=0.25) + theme_vyom + 
  scale_fill_manual(values = colorsType) + expand_limits(y = c(0)) + 
  theme( axis.text.x = element_text(angle = 45, size =  6, vjust = 1, hjust = 1), panel.background = element_rect(fill = "white", colour = "Black"), axis.text.y = element_text(size = 6), axis.title.x = element_text(size = 7), axis.title.y = element_text(size = 7)) + 
  xlab("Cell Type") + ylab("Fraction of Cells") + labs(fill = "Age") + scale_y_continuous(expand = expansion(mult = c(0, .1)))
ggsave( "CART_CONTROL_frac_of_cells.pdf", width=3.75, height=3, units="in")
cart_control$Sen_cells



#housekeeping
#cart_control <- readRDS("./data/Seurat_Objects/CART_Control_Sobj.rds", refhook = NULL)
colorsType = c(
  Young = "#7CA1CC",
  Old = "#FF4902"
)

library(SeuratDisk)
SaveH5Seurat(cart_control, filename = "./data/Control_CAR_T.h5Seurat")
Convert("./data/Control_CAR_T.h5Seurat", dest = "h5ad", overwrite = TRUE)
write.csv(cart_control@meta.data, 'Control_CAR_T_metadata.csv')
