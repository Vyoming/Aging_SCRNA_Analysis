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
options (future.globals.maxSize = 4000 * 1024^8)
Data.list <- PrepSCTIntegration(object.list = Data.list, anchor.features = Data.features, 
                                verbose = FALSE)
Data.anchors <- FindIntegrationAnchors(object.list = Data.list, normalization.method = "SCT", 
                                       anchor.features = Data.features, verbose = FALSE)
rm(Data.list, d)
cart_obj <- IntegrateData(anchorset = Data.anchors, normalization.method = "SCT", 
                            verbose = TRUE)
rm(Data.anchors)
rm(uPAR_F_young, uPAR_M_young, uPAR_F_old, uPAR_M_old,UT_F_young, UT_young,UT_F_old, UT_M_old)
# Visulization and Clustering
cart_obj <- RunPCA(cart_obj)

VizDimLoadings(cart_obj, dims = 1:2, reduction = "pca")

DimPlot(cart_obj, reduction = "pca")
ElbowPlot(cart_obj, ndims = 100, reduction = "pca")

cart_obj <- RunUMAP(cart_obj, dims = 1:15, n.neighbors = 5, n.epochs = 500)


levels(factor(arasco_obj@meta.data$orig.ident))
arasco_obj[["Type"]] <- Idents(arasco_obj)
new.cluster.ids <- c("Arasco", "Arasco", "Control", "Control")
names(new.cluster.ids) <- levels(arasco_obj)
arasco_obj <- RenameIdents(arasco_obj, new.cluster.ids)


Idents(cart_obj) <- factor(cart_obj$orig.ident)

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

levels(as.factor(cart_obj$Sex))
levels(as.factor(cart_obj$Treatment))
levels(as.factor(cart_obj$orig.ident))
meta_check <- cart_obj@meta.data
cart_obj <- FindNeighbors(cart_obj, dims = 1:30)
cart_obj <- FindClusters(cart_obj, resolution = 1)
DimPlot(cart_obj, reduction = "umap", label = TRUE)

Idents(cart_obj) <- cart_obj$Cell_Type
markers <- FindAllMarkers(cart_obj, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 1)
markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
write.csv(markers,'Immutumor_Sigs_Per_Clust_anno.csv')

allmarkers <- FindAllMarkers(cart_obj, min.pct = 0.5)
#clustering
DotPlot_Sig <- c("Lgr5","Ascl2","Olfm4","Gkn3","S100a6","Ly6a","Anxa3", "Areg",'Cbx3','Larp1','Slc39a1','Hnf4a','Sox4','Mmp7','Dll1','Tff3',"Tubb5","Syce2","Stmn1","Fbxo5",'Cenpa','Ccna2','Ube2c','Cdkn3',"Apoa1","Apoa4","Fabp1","Adh6a","Chgb","Tac1","Tph1","Neurog3", "Muc2","Fcgbp","Atoh1","Agr2","Pou2f3","Avil","Tuba1a","Adh1","Lyz1","Defa17","Defa24","Ang4") 
DotPlot(cart_obj, features = DotPlot_Sig, assay = 'SCT') + labs(y= "Cell Type", x="") + scale_colour_distiller( palette ="RdYlBu") + scale_size(range = c(0, 1)) +
  theme(text = element_text(size=5), axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust= .01), axis.text.y = element_text(size = 6), axis.title.x = element_blank(), axis.title.y = element_blank())
DotPlot_Sig <- c('Ptprc',"Lgr5","Ascl2","Olfm4","Gkn3","S100a6","Ly6a","Anxa3", "Areg","Tubb5","Syce2","Stmn1","Fbxo5",'Cenpa','Ccna2','Ube2c','Cdkn3',"Apoa1","Apoa4","Fabp1","Adh6a",'Tmigd1', 'Fabp6', 'Slc51b', 'Slc51a', "Chgb","Tac1","Tph1","Neurog3", "Muc2","Fcgbp","Atoh1","Agr2","Lyz1","Defa17","Defa24","Ang4","Pou2f3","Avil","Tuba1a","Adh1",'Cd8a','Cd4','Ighm') 
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




Prop_table<- prop.table(x = table(cart_obj$Cell_Type, cart_obj$Sex), margin = 2)
Prop_Table <- as.data.frame(Prop_table, row.names = NULL, optional = FALSE,
                            make.names = TRUE, stringsAsFactors = default.stringsAsFactors())
Prop_Table1 <- Prop_Table
Prop_Table1$Var1 <- factor(Prop_Table1$Var1,levels = my_levels)

plot <- ggplot(data = Prop_Table1, aes(Var1, Freq, fill=Var2)) + geom_bar(position="dodge", stat="identity", na.rm = TRUE) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), panel.background = element_rect(fill = "white", colour = "Black")) + xlab("Cell Type") + ylab("Fraction of Cells") + labs(fill = "Sample") +scale_y_continuous(expand = expansion(mult = c(0, .1)))
plot
ggsave(file = paste0('Cart_Sex_Proportions.pdf'), width=4, height=3, units="in")

Prop_table<- prop.table(x = table(cart_obj$Cell_Type, cart_obj$Treatment), margin = 2)
Prop_Table <- as.data.frame(Prop_table, row.names = NULL, optional = FALSE,
                            make.names = TRUE, stringsAsFactors = default.stringsAsFactors())
Prop_Table1 <- Prop_Table
Prop_Table1$Var1 <- factor(Prop_Table1$Var1,levels = my_levels)

plot <- ggplot(data = Prop_Table1, aes(Var1, Freq, fill=Var2)) + geom_bar(position="dodge", stat="identity", na.rm = TRUE) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), panel.background = element_rect(fill = "white", colour = "Black")) + xlab("Cell Type") + ylab("Fraction of Cells") + labs(fill = "Sample") +scale_y_continuous(expand = expansion(mult = c(0, .1)))
plot
ggsave(file = paste0('Cart_treatment_Proportions.pdf'), width=4, height=3, units="in")

Prop_table<- prop.table(x = table(cart_obj$Cell_Type, cart_obj$Age), margin = 2)
Prop_Table <- as.data.frame(Prop_table, row.names = NULL, optional = FALSE,
                            make.names = TRUE, stringsAsFactors = default.stringsAsFactors())
Prop_Table1 <- Prop_Table
Prop_Table1$Var1 <- factor(Prop_Table1$Var1,levels = my_levels)

plot <- ggplot(data = Prop_Table1, aes(Var1, Freq, fill=Var2)) + geom_bar(position="dodge", stat="identity", na.rm = TRUE) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), panel.background = element_rect(fill = "white", colour = "Black")) + xlab("Cell Type") + ylab("Fraction of Cells") + labs(fill = "Sample") +scale_y_continuous(expand = expansion(mult = c(0, .1)))
plot
ggsave(file = paste0('Cart_age_Proportions.pdf'), width=4, height=3, units="in")

Prop_table<- prop.table(x = table(cart_obj$Cell_Type, cart_obj$Age), margin = 2)
Prop_Table <- as.data.frame(Prop_table, row.names = NULL, optional = FALSE,
                            make.names = TRUE, stringsAsFactors = default.stringsAsFactors())
Prop_Table1 <- Prop_Table
Prop_Table1$Var1 <- factor(Prop_Table1$Var1,levels = my_levels)

plot <- ggplot(data = Prop_Table1, aes(Var1, Freq, fill=Var2)) + geom_bar(position="dodge", stat="identity", na.rm = TRUE) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), panel.background = element_rect(fill = "white", colour = "Black")) + xlab("Cell Type") + ylab("Fraction of Cells") + labs(fill = "Sample") +scale_y_continuous(expand = expansion(mult = c(0, .1)))
plot
ggsave(file = paste0('Cart_Age_Proportions.pdf'), width=4, height=3, units="in")
#DE analysis

Idents(cart_obj) <- cart_obj$Treatment
DE_treat <- FindMarkers(cart_obj, ident.1 = "uPAR", ident.2 = "UT", test.use = "MAST", logfc.threshold = .15, min.pct = .1, assay = 'SCT')
EnhancedVolcano(DE_treat, lab = rownames(DE_sex), x = 'avg_log2FC', y = 'p_val_adj', title = 'uPAR vs UT', pCutoff = 10e-1, FCcutoff = 0.25, xlim = c(-1,1), ylim = c(0,320), subtitle = 'All Cells' )
write.csv(DE_treat,'DE_Cart_treatment.csv')

Idents(cart_obj) <- cart_obj$Sex
DE_sex1 <- FindMarkers(cart_obj, ident.1 = "Female", ident.2 = "Male", test.use = "MAST", logfc.threshold = .15, min.pct = .1, assay = 'SCT')
EnhancedVolcano(DE_sex1, lab = rownames(DE_sex1), x = 'avg_log2FC', y = 'p_val_adj', title = 'Female vs Male', pCutoff = 10e-1, FCcutoff = 0.25, xlim = c(-3.25,2.5), ylim = c(0,325), subtitle = 'All Cells' )
write.csv(DE_sex1,'DE_cart_sex.csv')


Idents(cart_obj) <- cart_obj$Age
DE_age <- FindMarkers(cart_obj, ident.1 = "Old", ident.2 = "Young", test.use = "MAST", logfc.threshold = .15, min.pct = .1, assay = 'SCT')
EnhancedVolcano(DE_sex1, lab = rownames(DE_sex1), x = 'avg_log2FC', y = 'p_val_adj', title = 'Female vs Male', pCutoff = 10e-1, FCcutoff = 0.25, xlim = c(-3.25,2.5), ylim = c(0,325), subtitle = 'All Cells' )
write.csv(DE_age,'DE_cart_obj_vs_young.csv')

DotPlot_Sig <- c("Lgr5","Ascl2","Olfm4","Gkn3","S100a6","Ly6a","Anxa3", "Areg","Tubb5","Syce2","Stmn1","Fbxo5",'Cenpa','Ccna2','Ube2c','Cdkn3',"Apoa1","Apoa4","Fabp1","Adh6a",'Tmigd1', 'Fabp6', 'Slc51b', 'Slc51a', "Chgb","Tac1","Tph1","Neurog3", "Muc2","Fcgbp","Atoh1","Agr2","Lyz1","Defa17","Defa24","Ang4","Pou2f3","Avil","Tuba1a","Adh1",'Ptprc','Cd8a','Cd4','Ighm') 
DotPlot(cart_obj, features = DotPlot_Sig,group.by = 'Cell_Type' ,dot.scale = 10, scale.max= 100, assay = 'SCT') + labs(y= "Cell type", x="") + scale_colour_distiller( palette ="RdYlBu") + scale_size(range = c(0, 1.3)) +
  theme(legend.key.size = unit(.1, "in"), legend.text= element_text(size=4), legend.title = element_text(size=4),axis.line = element_line(size = .3), axis.ticks = element_line(size = .3), text = element_text(size=4), axis.text.x = element_text(size = 4, angle = 90, hjust = 1, vjust= .01), axis.text.y = element_text(size = 4), axis.title.x = element_blank(), axis.title.y = element_blank())
ggsave(file="Cart_Cluster_check.pdf", width=4, height=2)

my_levels <- c('Stem', 'Transit Amplifying', 'Enterocyte Progenitor','Enterocyte (Proximal)','Enterocyte (Distal)', 'Enteroendocrine', 'Goblet', 'Paneth','Tuft', 'Immune')
Idents(cart_obj) <- cart_obj$Cell_Type

for(i in my_levels){
  subset_cell <- subset(cart_obj,  idents = i)
  Idents(subset_cell) <- subset_cell$Treatment
  DE_sex <- FindMarkers(cart_obj, ident.1 = "uPAR", ident.2 = "UT", test.use = "MAST", logfc.threshold = .25, min.pct = ., assay = 'SCT')
  write.csv(DE_subset, paste0(i,'_DE.csv')) }


cart_immune_obj <- subset(cart_obj,  idents = 'Immune')

cart_immune_obj <- RunUMAP(cart_immune_obj, dims = 1:10)
cart_immune_obj <- FindNeighbors(cart_immune_obj, dims = 1:10)
cart_immune_obj <- FindClusters(cart_immune_obj, resolution = 1)
DimPlot(cart_immune_obj, reduction = "umap", label = TRUE)


DotPlot_Sig <- c("Lgr5","Ascl2","Olfm4","Gkn3","S100a6","Ly6a","Anxa3", "Areg","Tubb5","Syce2","Stmn1","Fbxo5",'Cenpa','Ccna2','Ube2c','Cdkn3',"Apoa1","Apoa4","Fabp1","Adh6a",'Tmigd1', 'Fabp6', 'Slc51b', 'Slc51a', "Chgb","Tac1","Tph1","Neurog3", "Muc2","Fcgbp","Atoh1","Agr2","Lyz1","Defa17","Defa24","Ang4","Pou2f3","Avil","Tuba1a","Adh1",'Cd8a','Cd8b','Cd4','Ighm') 
DotPlot(cart_immune_obj, features = DotPlot_Sig,dot.scale = 10, scale.max= 100, assay = 'SCT') + labs(y= "Cell type", x="") + scale_colour_distiller( palette ="RdYlBu") + scale_size(range = c(0, 1.3)) +
  theme(legend.key.size = unit(.1, "in"), legend.text= element_text(size=4), legend.title = element_text(size=4),axis.line = element_line(size = .3), axis.ticks = element_line(size = .3), text = element_text(size=4), axis.text.x = element_text(size = 4, angle = 90, hjust = 1, vjust= .01), axis.text.y = element_text(size = 4), axis.title.x = element_blank(), axis.title.y = element_blank())

FeaturePlot(Cart_immune, features = 'Ptprc')
FeaturePlot(cart_obj, features = 'Plaur')
cart_obj$Treatment
VlnPlot(cart_obj, group.by = 'Cell_Type' , ncol = 1, split.by = "Treatment", features = c('Plaur'), pt.size = 0, assay = "SCT", cols = c('#1b9e77' ,'#d95f02')) + theme(legend.position = 'none') + geom_boxplot(width=0.5,outlier.shape = NA, coef = 0) + xlab('')

All_Genes <- cart_obj@assays$RNA@data@Dimnames[[1]]
library(readr)
HALLMARK_INFLAMMATORY_RESPONSE_v7_5_1 <- read_delim("/Users/vyom/documents/HALLMARK_INFLAMMATORY_RESPONSE.v7.5.1.gmt", 
                                                    delim = "\t", escape_double = FALSE, 
                                                    col_names = FALSE, trim_ws = TRUE)
library('nichenetr')
H_inflammatory <- as.character(HALLMARK_INFLAMMATORY_RESPONSE_v7_5_1[1,])
M_inflammatory = H_inflammatory %>% convert_human_to_mouse_symbols()
M_inflammatory <- unique(M_inflammatory)
Inflammatory_Response <- intersect(All_Genes, M_inflammatory)


mean.exp <- zscore(colMeans(x = cart_obj@assays$MAGIC_RNA@data[Inflammatory_Response, ], na.rm = TRUE), dist = 'norm')
if (all(names(x = mean.exp) == rownames(x = cart_obj@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  cart_obj@meta.data$Inflammatory_Response <- mean.exp
}
mean.exp <- zscore(colMeans(x = cart_obj@assays$MAGIC_RNA@data[Inflammatory_Response, ], na.rm = TRUE), dist = 'norm')
if (all(names(x = mean.exp) == rownames(x = cart_obj@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  cart_obj@meta.data$Inflammatory_Response <- mean.exp
}
mean.exp <- zscore(colMeans(x = cart_young@assays$MAGIC_RNA@data[Inflammatory_Response, ], na.rm = TRUE), dist = 'norm')
if (all(names(x = mean.exp) == rownames(x = cart_young@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  cart_young@meta.data$Inflammatory_Response <- mean.exp
}

#magic - imputation
DefaultAssay(cart_obj) <- "RNA"
all.genes <- rownames(cart_obj)
cart_obj <- NormalizeData(object = cart_obj, normalization.method = "LogNormalize", assay = "RNA")
cart_obj <- ScaleData(cart_obj, features = all.genes)
cart_obj <- RunPCA(cart_obj, features = all.genes)

cart_obj <- magic(cart_obj)

#Gene level analysis

gene.list <- c('Plaur','Myc','Sox9', 'Olfm4', 'Lgr5', 'S100a6', 'Hopx', 'Fabp1', 'Crebbp','Ep300')
gene.list <- c('Cd74','Ptprc','Cd74','H2-Ab1','H2-Aa','H2-Eb1','Ccl5','Gzmb','H2-D1','Ifi208', 'Rsad2', 'Ctsw', 'Cxcl10', 'Havcr2', 'Oxsr1', 'Il10', 'Ccl4', 'Ccl3')
i = 'Inflammatory_Response'
for(i in gene.list) {
  selected_cells <- names(cart_obj$Cell_Type)
  vln_data <- FetchData(cart_obj,
                        vars = c(i,"Treatment", "Cell_Type"),
                        cells = selected_cells,
                        slot = "data")
  vln_data <- melt(vln_data)
  vln_data <- vln_data[c('Treatment', 'Cell_Type', 'value')]
  
  #  my_comparisons <- list(c('WT_Control, WT_Arasco'),c('WT_Control, WT_Rev'))
  All <- VlnPlot(cart_obj,group.by = 'Cell_Type' , split.by = "Treatment", features = i, pt.size = 0, assay = "MAGIC_RNA", cols = c('#1b9e77' ,'#d95f02'), log = FALSE, split.plot = TRUE) + 
    theme(legend.position = 'none')  +
    theme(axis.line = element_line(size = .3),axis.title.x = element_blank(), axis.ticks = element_line(size = .3), text = element_text(size=7), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6))  +
    stat_compare_means(data= vln_data, aes(x = Cell_Type, y = value, fill = Treatment), hide.ns = TRUE, method = "t.test",  label = "p.signif", size = 1)  #+ scale_y_continuous(expand = expansion(mult = c(0, .1)))
  dodge <- position_dodge(width = .9)
  All <- All + geom_boxplot(width=0.3,outlier.shape = NA, coef = 0, lwd=.2, position = dodge) +
    stat_compare_means(data= vln_data, aes(x = Cell_Type, y = value, fill = Treatment), hide.ns = TRUE, method = "t.test",  label = "p.signif", size = 1)  
  All$layers[[1]]$aes_params$size = .15
  ggsave(file = paste0('Cart_stem_', i, '.pdf'), plot=All, width=2, height=2, units="in")
}

levels(cart_obj$Cell_Type)
subset_cell <- subset(cart_obj,  idents = i)
subset_cell <- subset(cart_obj,  idents = i)
cart_obj$Age
celltypes <- c('Stem', 'Transit Amplifying', 'Enterocyte Progenitor','Enterocyte (Proximal)','Enterocyte (Distal)', 'Enteroendocrine', 'Goblet', 'Paneth','Tuft', 'T Cell', 'B Cell', 'Myeloid')
celltypes <- c( 'T Cell', 'B Cell', 'Myeloid')
sexes <- c('Female', 'Male')
ages <- c('Old', 'Young')
Idents(cart_obj) <- cart_obj$Cell_Type
for(i in celltypes){
  subset_cell <- subset(cart_obj,  idents = i)
  Idents(subset_cell) <- subset_cell$Treatment
  DE_subset <- FindMarkers(subset_cell, ident.1 = "uPAR", ident.2 = "UT", test.use = "MAST", logfc.threshold = .25, min.pct = .05, assay = 'SCT')
  write.csv(DE_subset, paste0(i,'_DE.csv')) 
  for(j in sexes){
    Idents(subset_cell) <- subset_cell$Sex
    subset_cell1 <- subset(subset_cell,  idents = j)
    Idents(subset_cell1) <- subset_cell1$Treatment
    DE_subset <- FindMarkers(subset_cell1, ident.1 = "uPAR", ident.2 = "UT", test.use = "MAST", logfc.threshold = .25, min.pct = .05, assay = 'SCT')
    write.csv(DE_subset, paste0(j,'_',i,'_DE.csv')) }
  for(k in ages){
    Idents(subset_cell) <- subset_cell$Age
    subset_cell1 <- subset(subset_cell,  idents = k)
    Idents(subset_cell1) <- subset_cell1$Treatment
    DE_subset <- FindMarkers(subset_cell1, ident.1 = "uPAR", ident.2 = "UT", test.use = "MAST", logfc.threshold = .25, min.pct = .05, assay = 'SCT')
    write.csv(DE_subset, paste0(k,'_',i,'_DE.csv')) }
  }
k ='Young'
Idents(cart_obj) <- cart_obj$Age
subset_cell1 <- subset(cart_obj,  idents = k)
gene.list = 'Plaur'
for(i in gene.list) {
  selected_cells <- names(subset_cell1$Cell_Type)
  vln_data <- FetchData(subset_cell1,
                        vars = c(i,"Treatment", "Cell_Type"),
                        cells = selected_cells,
                        slot = "data")
  vln_data <- melt(vln_data)
  vln_data <- vln_data[c('Treatment', 'Cell_Type', 'value')]
  
  #  my_comparisons <- list(c('WT_Control, WT_Arasco'),c('WT_Control, WT_Rev'))
  All <- VlnPlot(subset_cell1,group.by = 'Cell_Type' , split.by = "Treatment", features = i, pt.size = 0, assay = "MAGIC_RNA", cols = c('#1b9e77' ,'#d95f02'), log = FALSE, split.plot = TRUE) + 
    theme(legend.position = 'none')  +
    theme(axis.line = element_line(size = .3),axis.title.x = element_blank(), axis.ticks = element_line(size = .3), text = element_text(size=7), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6))  +
    stat_compare_means(data= vln_data, aes(x = Cell_Type, y = value, fill = Treatment), hide.ns = TRUE, method = "t.test",  label = "p.signif", size = 1)  #+ scale_y_continuous(expand = expansion(mult = c(0, .1)))
  dodge <- position_dodge(width = .9)
  All <- All + geom_boxplot(width=0.3,outlier.shape = NA, coef = 0, lwd=.2, position = dodge) +
    stat_compare_means(data= vln_data, aes(x = Cell_Type, y = value, fill = Treatment), hide.ns = TRUE, method = "t.test",  label = "p.signif", size = 1)  
  All$layers[[1]]$aes_params$size = .15
  ggsave(file = paste0('Cart_age_young_', i, '.pdf'), plot=All, width=2, height=2, units="in")
}

DotPlot_Sig <- c("Lgr5","Ascl2","Olfm4","Gkn3","Tubb5","Syce2","Stmn1","Fbxo5",'Cenpa','Ccna2','Ube2c','Cdkn3',"Apoa1","Apoa4","Fabp1","Adh6a",'Tmigd1','Fabp6','Slc51b','Slc51a',"Chgb","Tac1","Tph1","Neurog3", "Muc2","Fcgbp","Atoh1","Agr2","Lyz1","Defa17","Defa24","Ang4","Pou2f3","Avil","Tuba1a","Adh1", 'Ptprc', 'Cd8a','Cd4','Cd19','Ighm','Cd74') 
DotPlot(cart_obj, features = DotPlot_Sig,group.by = 'Cell_Type', dot.scale = 10, scale.max= 100, assay = 'SCT') + labs(y= "Cell type", x="") + scale_colour_distiller( palette ="RdYlBu") + scale_size(range = c(0, 1.3)) +
  theme(legend.key.size = unit(.1, "in"), legend.text= element_text(size=4), legend.title = element_text(size=4),axis.line = element_line(size = .3), axis.ticks = element_line(size = .3), text = element_text(size=4), axis.text.x = element_text(size = 4, angle = 90, hjust = 1, vjust= .01), axis.text.y = element_text(size = 4), axis.title.x = element_blank(), axis.title.y = element_blank())
ggsave(file="Cart_Cluster_Check.pdf", width=5, height=3)

#old score
senescence <- c('Acvr1b',	'Ang',	'Angpt1',	'Angptl4',	'Areg',	'Axl',	'Bex3',	'Bmp2',	'Bmp6',	'C3',	'Ccl1',	'Ccl2',	'Ccl20',	'Ccl24',	'Ccl26',	'Ccl3',	'Ccl4',	'Ccl5',	'Ccl7',	'Ccl8',	'Cd55',	'Cd9',	'Csf1',	'Csf2',	'Csf2rb',	'Cst10',	'Ctnnb1',	'Ctsb',	'Cxcl1',	'Cxcl10',	'Cxcl12',	'Cxcl16',	'Cxcl2',	'Cxcl3',	'Cxcr2',	'Dkk1',	'Edn1',	'Egf',	'Egfr',	'Ereg',	'Esm1',	'Ets2',	'Fas',	'Fgf1',	'Fgf2',	'Fgf7',	'Gdf15',	'Gem',	'Gmfg',	'Hgf',	'Hmgb1',	'Icam1',	'Icam5',	'Igf1',	'Igfbp1',	'Igfbp2',	'Igfbp3',	'Igfbp4',	'Igfbp5',	'Igfbp6',	'Igfbp7',	'Il10',	'Il13',	'Il15',	'Il18',	'Il1a',	'Il1b',	'Il2',	'Il6',	'Il6st',	'Il7',	'Inha',	'Iqgap2',	'Itga2',	'Itpka',	'Jun',	'Kitl',	'Lcp1',	'Mif',	'Mmp13',	'Mmp10',	'Mmp12',	'Mmp13',	'Mmp14',	'Mmp2',	'Mmp3',	'Mmp9',	'Nap1l4',	'Nrg1',	'Pappa',	'Pecam1',	'Pgf',	'Pigf',	'Plat',	'Plau',	'Plaur',	'Ptbp1',	'Ptger2',	'Ptges',	'Rps6ka5',	'Scamp4',	'Selplg',	'Sema3f',	'Serpinb3a',	'Serpine1',	'Serpine2',	'Spp1',	'Spx',	'Timp2',	'Tnf',	'Tnfrsf11b',	'Tnfrsf1a',	'Tnfrsf1b',	'Tubgcp2',	'Vegfa',	'Vegfc',	'Vgf',	'Wnt16',	'Wnt2')
All_Genes <- cart_obj@assays$RNA@data@Dimnames[[1]]
senescence <- unique(senescence)
senescence <- intersect(All_Genes, senescence)

mean.exp <- zscore(colMeans(x = cart_obj@assays$MAGIC_RNA@data[senescence, ], na.rm = TRUE), dist = 'norm')
if (all(names(x = mean.exp) == rownames(x = cart_obj@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  cart_obj@meta.data$senescence <- mean.exp
}


Idents(cart_obj) <- cart_obj$Treatment
Cart_control <- subset(cart_obj, idents = 'UT')
gene.list = 'senescence'
DefaultAssay(Cart_control) <- 'MAGIC_RNA'

for(i in senescence) {
  selected_cells <- names(Cart_control$Cell_Type)
  vln_data <- FetchData(Cart_control,
                        vars = c(i,"Age", "Cell_Type"),
                        cells = selected_cells,
                        slot = "data")
  vln_data <- melt(vln_data)
  vln_data <- vln_data[c('Age', 'Cell_Type', 'value')]
  
  #  my_comparisons <- list(c('WT_Control, WT_Arasco'),c('WT_Control, WT_Rev'))
  All <- VlnPlot(Cart_control,group.by = 'Cell_Type' , split.by = "Age", features = i, pt.size = 0, assay = "MAGIC_RNA", cols = c('#7CA1CC' ,'#FF4902'), log = FALSE, split.plot = TRUE) + 
    theme(legend.position = 'none')  +
    theme(axis.line = element_line(size = .3),axis.title.x = element_blank(), axis.ticks = element_line(size = .3), text = element_text(size=7), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6))  +
    stat_compare_means(data= vln_data, aes(x = Cell_Type, y = value, fill = Age), hide.ns = TRUE, method = "t.test",  label = "p.signif", size = 1)  #+ scale_y_continuous(expand = expansion(mult = c(0, .1)))
  dodge <- position_dodge(width = .9)
  All <- All + geom_boxplot(width=0.3,outlier.shape = NA, coef = 0, lwd=.2, position = dodge) +
    stat_compare_means(data= vln_data, aes(x = Cell_Type, y = value, fill = Age), hide.ns = TRUE, method = "t.test",  label = "p.signif", size = 1)  
  All$layers[[1]]$aes_params$size = .15
  ggsave(file = paste0('Cart_age_UT_only_', i, '.pdf'), plot=All, width=2, height=2, units="in")
}
All <- VlnPlot(Cart_control, split.by = "Age", features = i, pt.size = 0, assay = "MAGIC_RNA", cols = c('#1b9e77' ,'#d95f02'), log = FALSE, split.plot = TRUE) + 
  theme(legend.position = 'none')  +
  theme(axis.line = element_line(size = .3),axis.title.x = element_blank(), axis.ticks = element_line(size = .3), text = element_text(size=7), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6))  
dodge <- position_dodge(width = .9)
All$layers[[1]]$aes_params$size = .15
All
ggsave(file = paste0('Cart_age_overall_', i, '.pdf'), plot=All, width=2, height=2, units="in")

Idents(cart_obj) <- cart_obj$Cell_Type
FeaturePlot(object = cart_obj, features = 'senescence', pt.size = .001, label = TRUE)  + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
  theme(plot.title = element_blank(), text = element_text(size=6), legend.position = 'none',axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), axis.title.x = element_blank(), axis.title.y = element_blank())
ggsave(file = 'CART_UT_Scenscence.pdf', width=5, height=5, units="in")

fp <- plot_density(cart_obj, reduction = "umap", 'senescence', combine=F, pal="inferno") 
ggsave(filename = paste0('CART_UT_',i,'_signatures_density_plots.pdf'), plot = fp, width=3.5, height=2.5, units="in" )


Idents(Cart_control) <- Cart_control$Age
DE_age <- FindMarkers(Cart_control, ident.1 = "Old", ident.2 = "Young", test.use = "MAST", logfc.threshold = .15, min.pct = .1, assay = 'SCT')
EnhancedVolcano(DE_age, lab = rownames(DE_age), x = 'avg_log2FC', y = 'p_val_adj', title = 'uPAR vs UT', pCutoff = 10e-1, FCcutoff = 0.25, xlim = c(-1,1), ylim = c(0,320), subtitle = 'All Cells' )
write.csv(DE_age,'DE_Cart_age.csv')

#gene set enrichment analysis
DE_age
library(biomaRt)
library('fgsea')
mart <- useDataset("hsapiens_gene_ensembl", mart=useMart("ensembl"))
bm <- getBM(attributes=c("ensembl_gene_id", "Gene name"), mart=mart) %>%
  distinct() %>%
  as_tibble() %>%
  na_if("") %>% 
  na.omit()
bm
DE_age$gene_id <- convert_mouse_to_human_symbols(rownames(DE_age))
DE_age1 <- DE_age[complete.cases(DE_age),]
res2 <- DE_age1 %>% 
  dplyr::select(gene_id, avg_log2FC) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(gene_id) %>% 
  summarize(avg_log2FC=mean(avg_log2FC))
res2

ranks <- deframe(res2)
head(ranks, 20)
pathways.hallmark <- gmtPathways("Documents/c5.all.v2022.1.Hs.symbols.gmt")
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
write.csv(gsea_output,'CarT_GSEA_Old_vs_Young.csv')
colnames(gsea_output)
theTable$Position <- factor(theTable$Position, levels = c(...))

pathwayColors <- c( "#FDC926FF", "#FA9E3BFF", "#ED7953FF", "#D8576BFF", "#BD3786FF", "#9C179EFF", "#7301A8FF", "#47039FFF")
gsea_output1 <- gsea_output[order(-gsea_output$NES),]
gsea_output1$pathway <- factor(gsea_output1$pathway, levels = gsea_output1$pathway)
All <- ggplot(data = gsea_output1,aes(x=pathway,y=NES)) +  theme_vyom +
  theme(axis.text.y = element_text(size = 3) ,axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), panel.background = element_rect(fill = "white", colour = "Black", size = .75, linetype = 'solid')) +
  geom_segment( aes(x=pathway, xend=pathway, y=0, yend=NES, color = padj))+
  geom_point(color = 'black' ,size = 1.2 ) +
  geom_point(aes(color = padj) ,size = 1 ) + geom_vline(xintercept = 0) +
  scale_color_gradientn(colors=pathwayColors, breaks = c(.05, .01, .001, .00001),labels =  c('.05', '.01', '.001', '<.0001'),limits = c(0,.05), trans = scales::boxcox_trans(.25)) +
  coord_flip() + 
  scale_y_continuous(limits = c(-3,3),expand = expansion(mult = c(0, 0)), breaks = scales::breaks_extended(n = 3)) +
  labs(y= "Normalized Enrichment Score", x="Pathway") 
All + geom_vline(xintercept = 0)
ggsave(file = paste0('Fig 1 C_GSEA.pdf'), plot=All, width=5, height=3, units="in")


CREB_Targets <- read_excel("analysis/Finalized_Signatures.xlsx", sheet = "Creb_Targets")
YAP_Targets <- read_excel("analysis/Finalized_Signatures.xlsx", sheet = "Yap_Targets")
BCatenin_Enhanced <- read_excel("analysis/Finalized_Signatures.xlsx", sheet = "bcatenin_enhanced")
PPAR_SCORE <- read.delim("~/Documents/KEGG_PPAR_SIGNALING_PATHWAY.v2022.1.Hs.gmt", header=FALSE)
PPAR_SCORE <- as.character(PPAR_SCORE[1,])
PPAR_SCORE = PPAR_SCORE %>% convert_human_to_mouse_symbols()
PPAR_targets <- unique(PPAR_SCORE)

All_Genes <- cart_obj@assays$RNA@data@Dimnames[[1]]
PPAR_targets <- intersect(All_Genes, PPAR_targets)
CREB_Targets_Genes <- intersect(All_Genes, CREB_Targets$Gene_Name)
BCatenin_Enhanced_Genes <- intersect(All_Genes, BCatenin_Enhanced$Gene_Name)
YAP_Targets_Genes <- intersect(All_Genes, YAP_Targets$Gene_Name)

mean.exp <- zscore(colMeans(x = cart_obj@assays$MAGIC_RNA@data[PPAR_targets, ], na.rm = TRUE), dist ='norm')
if (all(names(x = mean.exp) == rownames(x = cart_obj@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  cart_obj@meta.data$PPAR_Targets <- mean.exp
}

mean.exp <- zscore(colMeans(x = cart_obj@assays$MAGIC_RNA@data[CREB_Targets_Genes, ], na.rm = TRUE), dist ='norm')
if (all(names(x = mean.exp) == rownames(x = cart_obj@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  cart_obj@meta.data$CREB_Targets_Score <- mean.exp
}

mean.exp <- zscore(colMeans(x = cart_obj@assays$MAGIC_RNA@data[BCatenin_Enhanced_Genes, ], na.rm = TRUE), dist ='norm')
if (all(names(x = mean.exp) == rownames(x = cart_obj@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  cart_obj@meta.data$BCatenin_Enhanced_Score <- mean.exp
}

mean.exp <- zscore(colMeans(x = cart_obj@assays$MAGIC_RNA@data[YAP_Targets_Genes, ], na.rm = TRUE), dist ='norm')
if (all(names(x = mean.exp) == rownames(x = cart_obj@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  cart_obj@meta.data$YAP_Targets_Score <- mean.exp
}
Gene.list <- c('Lgr4','Myc', 'Sox9', 'Olfm4', 'Malat1', 'Hopx', 'Ccnd1' )
mean.exp <- zscore(colMeans(x = cart_obj@assays$MAGIC_RNA@data[Gene.list, ], na.rm = TRUE), dist ='norm')
if (all(names(x = mean.exp) == rownames(x = cart_obj@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  cart_obj@meta.data$Stem_Score <- mean.exp
}

cart_obj_stem <- AddModuleScore(object = cart_obj_stem, features = list('Ascl2'), name = 'Ascl2')

DefaultAssay(cart_obj) <- "MAGIC_RNA"
Gene.list <- c('YAP_Targets_Score','BCatenin_Enhanced_Score','CREB_Targets_Score','PPAR_Targets'   )
for(i in Gene.list) {
  selected_cells <- names(cart_obj$Cell_Type)
  vln_data <- FetchData(cart_obj,
                        vars = c(i,"Treatment", "Cell_Type"),
                        cells = selected_cells,
                        slot = "data")
  vln_data$Type
  vln_data <- melt(vln_data)
  
  All <- VlnPlot(cart_obj, split.by = "Treatment",group.by = 'Cell_Type', features = i, pt.size = 0, assay = "MAGIC_RNA",  cols = c('#1b9e77' ,'#d95f02', '#7570b3'), log = TRUE, split.plot = TRUE) + 
    theme(legend.position = 'none') + 
    geom_boxplot(width=0.3,outlier.shape = NA, coef = 0, lwd= .2) + xlab('') + 
    theme(text = element_text(size=9), axis.text.x = element_text(size = 9), axis.text.y = element_text(size =9)) +
    stat_compare_means(data= vln_data, aes(x = Cell_Type, y = value, fill = Treatment), method = "wilcox.test",  label = "p.signif", size = 2, hide.ns = TRUE) #+ scale_y_continuous(expand = expansion(mult = c(0, .1)))
  All$layers[[1]]$aes_params$size = .15
  All
  ggsave(file = paste0('CART_Score_', i, '.pdf'), plot=All, width=3, height=3, units="in")
}
Idents(cart_obj) <- cart_obj$Age
cart_obj <- subset(cart_obj, idents = 'Old')
Idents(cart_obj) <- cart_obj$Treatment
cart_obj_UT <- subset(cart_obj, idents = 'UT')
Cart_y <- subset(cart_obj, idents = 'Young')
Idents(Cart_y) <- Cart_y$Cell_Type
Cart_y_stem <- subset(Cart_y, idents = 'Stem')

DefaultAssay(cart_obj) <- "MAGIC_RNA"
Gene.list <- c('YAP_Targets_Score','BCatenin_Enhanced_Score','CREB_Targets_Score','PPAR_Targets'   )
Gene.list <- c('H2-Ab1','Inflammatory_Response'   )
Gene.list <- c(senescence)
for(i in Gene.list) {
  selected_cells <- names(cart_obj$Cell_Type)
  vln_data <- FetchData(cart_obj,
                        vars = c(i,"Treatment", "Cell_Type"),
                        cells = selected_cells,
                        slot = "data")
  vln_data$Type
  vln_data <- melt(vln_data)
  
  All <- VlnPlot(cart_obj, split.by = "Treatment",group.by = 'Cell_Type', features = i, pt.size = 0, assay = "MAGIC_RNA",  cols = c('#7CA1CC' ,'#FF4902'), log = TRUE, split.plot = TRUE) + 
    theme(legend.position = 'none') + 
    geom_boxplot(width=0.3,outlier.shape = NA, coef = 0, lwd= .2) + xlab('') + 
    theme(text = element_text(size=9), axis.text.x = element_text(size = 9), axis.text.y = element_text(size =9)) +
    stat_compare_means(data= vln_data, aes(x = Cell_Type, y = value, fill = Treatment), method = "wilcox.test",  label = "p.signif", size = 2, hide.ns = TRUE) #+ scale_y_continuous(expand = expansion(mult = c(0, .1)))
  All$layers[[1]]$aes_params$size = .15
  All
  ggsave(file = paste0('CART_Score_', i, '.pdf'), plot=All, width=3, height=3, units="in")
}

for(i in Gene.list) {
  selected_cells <- names(cart_obj$Cell_Type)
  vln_data <- FetchData(cart_obj,
                        vars = c(i,"Treatment", "Cell_Type"),
                        cells = selected_cells,
                        slot = "data")
  vln_data$Type
  vln_data <- melt(vln_data)
  
  All <- VlnPlot(cart_obj, split.by = "Treatment",group.by = 'Cell_Type', features = i, pt.size = 0, assay = "MAGIC_RNA",  cols = c('#7CA1CC' ,'#FF4902'), log = TRUE, split.plot = TRUE) + 
    theme(legend.position = 'none') + 
    geom_boxplot(width=0.3,outlier.shape = NA, coef = 0, lwd= .2) + xlab('') + 
    theme(text = element_text(size=9), axis.text.x = element_text(size = 9), axis.text.y = element_text(size =9)) +
    stat_compare_means(data= vln_data, aes(x = Cell_Type, y = value, fill = Treatment), method = "wilcox.test",  label = "p.signif", size = 2, hide.ns = TRUE) #+ scale_y_continuous(expand = expansion(mult = c(0, .1)))
  All$layers[[1]]$aes_params$size = .15
  All
  ggsave(file = paste0('CART_Score_', i, '.pdf'), plot=All, width=3, height=3, units="in")
}

cart_obj$Treatment

Idents(cart_obj) <- cart_obj$Cell_Type
cart_obj_stem <- subset(cart_obj, idents = 'Stem')
DefaultAssay(cart_obj_stem) <- "MAGIC_RNA"
Gene.list <- c('YAP_Targets_Score','BCatenin_Enhanced_Score','CREB_Targets_Score','PPAR_Targets'   )
Gene.list <- unique(c('Ifng', 'Il1b', 'Il2', 'Il4', 'Il6', 'Il10', 'Il15', 'Ccl2', 'Ccl3', 'Ccl4', 'Ccl5',  'Ccl12', 'Ccl20', 'Ccl17', 'Ifngr1',	'Irf7',	'Irf5',	'Stat2',	'Stat3',	'Isg15',	'Isg20',	'Ifi205',	'Ifi209',	'Ifi208',	'Ifi211',	'Ifi207',	'Ifi203',	'Ifi204',	'Ifi206',	'Ifi47',	'Ifi204',	'Ifi213',	'Ifi205',	'Ifi30',	'Ifit1',	'Ifit2',	'Ifit1bl1',	'Ifit3',	'Ifit3b',	'Ifitm2',	'Ifitm6',	'Ifitm3',	'H2-D1',	'H2-Q7',	'H2-T22',	'H2-Ab1',	'H2-Aa',	'H2-DMb2',	'Ctsb',	'Ctsd',	'Tap1',	'Cxcl10',	'Ccl9',	'Ccl6',	'Ccl3',	'Ccl2',	'Ccl7',	'Ccl12',	'Ccr1',	'Il1b',	'Il1r2',	'Il15',	'Junb',	'Jun',	'Jund',	'Fos',	'Cox8a',	'Cox4i1',	'Ucp2',	'Plin2'))
i = 'Cd74'
cart_obj$Treatment


Idents(cart_obj_stem) <- cart_obj_stem$Treatment
DE_cart <- FindMarkers(cart_obj_stem, ident.1 = "uPAR", ident.2 = "UT", test.use = "MAST", logfc.threshold = .01, min.pct = .01, assay = 'SCT')
EnhancedVolcano(DE_cart, lab = rownames(DE_cart), x = 'avg_log2FC', y = 'p_val_adj', title = 'uPAR vs UT', pCutoff = 10e-1, FCcutoff = 0.25, xlim = c(-1,1), ylim = c(0,320), subtitle = 'All Cells' )
write.csv(DE_cart,'DE_Cart_treatment.csv')

source("~/analysis/split_violin.R")
my_levels <- c('UT', 'uPAR' )
sense
Gene.list <- c('Psap','Timp2','Lgmn','Igf1','Arl6ip1','Lsmem1','Igfbp4','Gsn','Zbtb20','Ccl2','Ccl7','Trio','Tulp4','Cxcl1')
Gene.list <- c('Psap','Arl6ip1','Igfbp4','Gsn','Zbtb20','Tulp4')


DefaultAssay(cart_young) <- 'MAGIC_RNA'
cart_young$Treatment <- factor(x = cart_young$Treatment, levels = my_levels)
DefaultAssay(cart_young@assays$MAGIC_RNA) <- "MAGIC_RNA"
selected_cells <- names(cart_young$Cell_Type)

#selected_cells <- names(cart_young$Cell_Type[cart_young$Cell_Type == c("Stem")])
vln_data <- FetchData(cart_young,
                      vars = c(Gene.list,"Treatment"),
                      cells = selected_cells,
                      slot = "data")

vln_data <- melt(vln_data)

All <- ggplot(vln_data, aes(x = variable, y = value, fill= Treatment)) + geom_split_violin(scale="width", trim = TRUE, size = .1) + theme_vyom + theme(legend.position = 'none') + 
  geom_boxplot(width=0.2,outlier.shape = NA, coef = 0, lwd=.2) + xlab('') + ylab('Expression Level') + scale_fill_manual(values= c('#7CA1CC' ,'#FF4902')) +
  theme(text = element_text(size=7), axis.text.x = element_text(size = 6, angle = 45, hjust = 1, vjust = 1), axis.text.y = element_text(size = 6), axis.title.y = element_text(size = 7)) +
  stat_compare_means( method = "wilcox.test",  label = "p.signif", size = 2) + scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
All
ggsave(file = 'cart_young_Senescnce.pdf', plot=All, width=2, height=2, units="in")
DefaultAssay(cart_obj) <- "SCT"
cart_obj <- AddModuleScore(object = cart_obj, features = Gene.list, name = 'Stem_score', assay = 'SCT')

my_levels <- c('UT', 'uPAR' )
Cart_y$Treatment <- factor(x = Cart_y$Treatment, levels = my_levels)
DefaultAssay(Cart_y@assays$MAGIC_RNA) <- "MAGIC_RNA"
selected_cells <- names(Cart_y$Cell_Type[Cart_y$Cell_Type == c("Stem")])
vln_data <- FetchData(Cart_y,
                      vars = c(Gene.list,"Treatment"),
                      cells = selected_cells,
                      slot = "data")

vln_data <- melt(vln_data)

All <- ggplot(vln_data, aes(x = variable, y = value, fill= Treatment)) + geom_split_violin(scale="width", trim = TRUE, size = .1) + theme_vyom + theme(legend.position = 'none') + 
  geom_boxplot(width=0.2,outlier.shape = NA, coef = 0, lwd=.2) + xlab('') + ylab('Expression Level') + scale_fill_manual(values= c('#7CA1CC' ,'#FF4902')) +
  theme(text = element_text(size=7), axis.text.x = element_text(size = 6, angle = 45, hjust = 1, vjust = 1), axis.text.y = element_text(size = 6), axis.title.y = element_text(size = 7)) +
  stat_compare_means( method = "wilcox.test",  label = "p.signif", size = 2) + scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
All
ggsave(file = 'Fig D-E YOUNG_ONLY_STEM.pdf', plot=All, width=2, height=2, units="in")

Idents(cart_obj_stem) <- cart_obj_stem$Cell_Type
my_levels <- c("Stem")
cart_obj_stem$Cell_Type <- factor(x = cart_obj_stem$Cell_Type, levels = my_levels)
Prop_table <- prop.table(x = table(cart_obj_stem$Treatment, cart_obj_stem$Cell_Type))
Prop_Table <- as.data.frame(Prop_table, row.names = NULL, optional = FALSE,
                            make.names = TRUE, stringsAsFactors = default.stringsAsFactors())

Prop_Table$Var1 <- factor(Prop_Table$Var1,levels = c('UT','uPAR'))
plot <- ggplot(data = Prop_Table, aes(Var2, Freq, fill=Var1)) + geom_bar(position="stack", stat="identity", na.rm = TRUE) + scale_fill_manual(values = c('#7CA1CC' ,'#FF4902')) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size= 8), axis.text.y  = element_text(size = 8), axis.title.x = element_text(size = 9), axis.title.y = element_text(size = 9),  panel.background = element_rect(fill = "white", colour = "Black")) + xlab("Cell Type") + ylab("Fraction of Stem Cells") + labs(fill = "Sample") +scale_y_continuous(expand = expansion(mult = c(0, .1)))
plot
ggsave(file="Endo_stem_prop.pdf", width=2.5, height=2)
prop.table(x = table(cart_obj$Treatment, cart_obj$Cell_Type))

#old v young check senescence signatures
Idents(cart_obj) <- cart_obj$Age
cart_old <- subset(cart_obj,  idents = 'Old')
cart_young <- subset(cart_obj,  idents = 'Young')

Idents(cart_obj) <- cart_obj$Age
cart_old <- subset(cart_obj, idents = 'Old')
Idents(cart_obj) <- cart_obj$Treatment
cart_obj_UT <- subset(cart_obj, idents = 'UT')
Cart_y <- subset(cart_obj, idents = 'Young')
Idents(Cart_y) <- Cart_y$Cell_Type
Cart_y_stem <- subset(Cart_y, idents = 'Stem')

cart_obj_UT
#this is correct muscle signature
senescence <- c('Psap','Timp2','Itm2b','Lgmn','Igf1','Arl6ip1','Ckb','Lsmem1','Igfbp4','Gsn','Zbtb20','Ccl2','Ccl7','Trio','Tulp4','Cxcl1')
All_Genes <- cart_obj@assays$RNA@data@Dimnames[[1]]
senescence <- intersect(All_Genes, senescence)
mean.exp <- zscore(colMeans(x = cart_obj_UT@assays$MAGIC_RNA@data[senescence, ], na.rm = TRUE), dist ='norm')
if (all(names(x = mean.exp) == rownames(x = cart_obj_UT@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  cart_obj_UT@meta.data$senescence <- mean.exp
}
mean.exp <- zscore(colMeans(x = cart_obj@assays$MAGIC_RNA@data[senescence, ], na.rm = TRUE), dist ='norm')
if (all(names(x = mean.exp) == rownames(x = cart_obj@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  cart_obj@meta.data$senescence <- mean.exp
}
mean.exp <- zscore(colMeans(x = cart_young@assays$MAGIC_RNA@data[senescence, ], na.rm = TRUE), dist ='norm')
if (all(names(x = mean.exp) == rownames(x = cart_young@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  cart_young@meta.data$senescence <- mean.exp
}
mean.exp <- zscore(colMeans(x = cart_obj@assays$MAGIC_RNA@data[senescence, ], na.rm = TRUE), dist ='norm')
if (all(names(x = mean.exp) == rownames(x = cart_obj@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  cart_obj@meta.data$senescence <- mean.exp
}
Gene.list <- c('Lgr4','Myc', 'Sox9', 'Olfm4', 'Malat1', 'Hopx', 'Ccnd1' )
mean.exp <- zscore(colMeans(x = cart_young@assays$MAGIC_RNA@data[Gene.list, ], na.rm = TRUE), dist ='norm')
if (all(names(x = mean.exp) == rownames(x = cart_young@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  cart_young@meta.data$Stem_Score <- mean.exp
}

mean.exp <- zscore(colMeans(x = cart_young@assays$MAGIC_RNA@data[PPAR_targets, ], na.rm = TRUE), dist ='norm')
if (all(names(x = mean.exp) == rownames(x = cart_young@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  cart_young@meta.data$PPAR_Targets <- mean.exp
}

mean.exp <- zscore(colMeans(x = cart_young@assays$MAGIC_RNA@data[CREB_Targets_Genes, ], na.rm = TRUE), dist ='norm')
if (all(names(x = mean.exp) == rownames(x = cart_young@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  cart_young@meta.data$CREB_Targets_Score <- mean.exp
}

mean.exp <- zscore(colMeans(x = cart_young@assays$MAGIC_RNA@data[BCatenin_Enhanced_Genes, ], na.rm = TRUE), dist ='norm')
if (all(names(x = mean.exp) == rownames(x = cart_young@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  cart_young@meta.data$BCatenin_Enhanced_Score <- mean.exp
}

mean.exp <- zscore(colMeans(x = cart_young@assays$MAGIC_RNA@data[YAP_Targets_Genes, ], na.rm = TRUE), dist ='norm')
if (all(names(x = mean.exp) == rownames(x = cart_young@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  cart_young@meta.data$YAP_Targets_Score <- mean.exp
}

mean.exp <- zscore(colMeans(x = cart_obj@assays$MAGIC_RNA@data[PPAR_targets, ], na.rm = TRUE), dist ='norm')
if (all(names(x = mean.exp) == rownames(x = cart_obj@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  cart_obj@meta.data$PPAR_Targets <- mean.exp
}

mean.exp <- zscore(colMeans(x = cart_obj@assays$MAGIC_RNA@data[CREB_Targets_Genes, ], na.rm = TRUE), dist ='norm')
if (all(names(x = mean.exp) == rownames(x = cart_obj@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  cart_obj@meta.data$CREB_Targets_Score <- mean.exp
}

mean.exp <- zscore(colMeans(x = cart_obj@assays$MAGIC_RNA@data[BCatenin_Enhanced_Genes, ], na.rm = TRUE), dist ='norm')
if (all(names(x = mean.exp) == rownames(x = cart_obj@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  cart_obj@meta.data$BCatenin_Enhanced_Score <- mean.exp
}

mean.exp <- zscore(colMeans(x = cart_obj@assays$MAGIC_RNA@data[YAP_Targets_Genes, ], na.rm = TRUE), dist ='norm')
if (all(names(x = mean.exp) == rownames(x = cart_obj@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  cart_obj@meta.data$YAP_Targets_Score <- mean.exp
}

i = 'senescence'
Gene.list <- c('YAP_Targets_Score','BCatenin_Enhanced_Score','CREB_Targets_Score','PPAR_Targets'   )
for(i in Gene.list) {
selected_cells <- names(cart_young$Cell_Type)
vln_data <- FetchData(cart_young,
                      vars = c(i,"Treatment", "Cell_Type"),
                      cells = selected_cells,
                      slot = "data")
vln_data <- melt(vln_data)
All <- VlnPlot(cart_young, split.by = "Treatment",group.by = 'Cell_Type', features = i, pt.size = 0, assay = "MAGIC_RNA",  cols = c('#7CA1CC' ,'#FF4902'), log = TRUE, split.plot = TRUE) + 
  theme(legend.position = 'none') + 
  geom_boxplot(width=0.3,outlier.shape = NA, coef = 0, lwd= .2) + xlab('') + 
  theme(text = element_text(size=9), axis.text.x = element_text(size = 9), axis.text.y = element_text(size =9)) +
  stat_compare_means(data= vln_data, aes(x = Cell_Type, y = value, fill = Treatment), method = "wilcox.test",  label = "p.signif", size = 2, hide.ns = TRUE) #+ scale_y_continuous(expand = expansion(mult = c(0, .1)))
All$layers[[1]]$aes_params$size = .15
All
ggsave(file = paste0('CART_young_', i, '.pdf'), plot=All, width=3, height=3, units="in")
}

Idents(cart_obj) <- cart_obj$Treatment
Cart_UT <- subset(cart_obj, idents = 'UT')

for(i in senescence) {
  selected_cells <- names(Cart_UT$Cell_Type)
  vln_data <- FetchData(Cart_UT,
                        vars = c(i,"Age", "Cell_Type"),
                        cells = selected_cells,
                        slot = "data")
  vln_data$Type
  vln_data <- melt(vln_data)
  
  All <- VlnPlot(Cart_UT, split.by = "Age",group.by = 'Cell_Type', features = i, pt.size = 0, assay = "MAGIC_RNA",  cols = c('#7CA1CC' ,'#FF4902'), log = TRUE, split.plot = TRUE) + 
    theme(legend.position = 'none') + 
    geom_boxplot(width=0.3,outlier.shape = NA, coef = 0, lwd= .2) + xlab('') + 
    theme(text = element_text(size=9), axis.text.x = element_text(size = 9), axis.text.y = element_text(size =9)) +
    stat_compare_means(data= vln_data, aes(x = Cell_Type, y = value, fill = Age), method = "wilcox.test",  label = "p.signif", size = 2, hide.ns = TRUE) #+ scale_y_continuous(expand = expansion(mult = c(0, .1)))
  All$layers[[1]]$aes_params$size = .15
  All
  ggsave(file = paste0('CART_age_', i, '.pdf'), plot=All, width=3, height=3, units="in")
}
Idents(cart_obj) <- cart_obj$Cell_Type
FeaturePlot(object = cart_obj, features = 'senescence', pt.size = .001, label = TRUE)  + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")), trans = scales::boxcox_trans(.25)) +
  theme(plot.title = element_blank(), legend.position = 'none', text = element_text(size=6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), axis.title.x = element_blank(), axis.title.y = element_blank())
ggsave(file = 'CART_UT_Scenscence_label.pdf', width=5, height=5, units="in")

scale_color_gradientn(colors=pathwayColors, breaks = c(.05, .01, .001, .00001),labels =  c('.05', '.01', '.001', '<.0001'),limits = c(0,.05), trans = scales::boxcox_trans(.25))
Cart_UT$Age
Idents(cart_obj)
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
      
      pvals <- c(pvals, t.test(vln_data[which(vln_data$Age == 'Old'),]$value, vln_data[which(vln_data$Age == 'Young'),]$value,  alternative = c("greater"))$p.value)
      print(i)
    }
  }
  
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




#housekeeping
# saveRDS(cart_obj, file = "./data/Seurat_Objects/CART_Sobj.rds")
# cart_obj <- readRDS("./data/Seurat_Objects/CART_Sobj.rds", refhook = NULL)

my_levels <- c('UT', 'uPAR')
cart_obj$Treatment <- factor(x = cart_obj$Treatment, levels = my_levels)
options (future.globals.maxSize = 4000 * 1024^8)




