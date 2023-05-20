#senescence on human data
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
library(SeuratDisk)

#Convert("~/data/epi_log_counts02_v2.h5ad", dest = "h5seurat", overwrite = TRUE)
sen_data <- LoadH5Seurat("~/data/epi_log_counts02_v2.h5seurat")
sen_data1 <- sen_data
levels(sen_data1$Age)
Idents(sen_data1) <- sen_data1$Age
sen_data1 <- subset(sen_data1,  idents = c('25-30','70-75'))
my_levels <- c('25-30','70-75')
sen_data1$Age <- factor(x = sen_data1$Age, levels = my_levels)

DimPlot(sen_data1, group.by = 'Age')
DimPlot(sen_data1, group.by = 'annotation')

Prop_table<- prop.table(x = table(sen_data1$Age, sen_data1$annotation), margin = 2)

Idents(sen_data1) <- sen_data1$annotation
sen_data1 <- subset(sen_data1,  idents = c('Stem cells', 'TA', 'Enterocyte', 'Colonocyte', 'BEST4+ epithelial', 'BEST2+ Goblet cell', 'Goblet cell', 'Paneth', 'Tuft'))    
my_levels <- c('Stem cells', 'TA', 'Colonocyte', 'BEST4+ epithelial', 'BEST2+ Goblet cell', 'Goblet cell', 'Paneth', 'Tuft', 'Enterocyte')
sen_data1$annotation <- factor(x = sen_data1$annotation, levels = my_levels)
#new
senescence <- c('Psap','Timp2','Itm2b','Lgmn','Igf1','Arl6ip1','Ckb','Lsmem1','Igfbp4','Gsn','Zbtb20','Ccl2','Ccl7','Trio','Tulp4','Cxcl1') %>% convert_mouse_to_human_symbols()

All_Genes <- sen_data1@assays$RNA@data@Dimnames[[1]]
senescence <- intersect(All_Genes, senescence)
mean.exp <- zscore(colMeans(x = sen_data1@assays$RNA@data[senescence, ], na.rm = TRUE), dist ='norm')
if (all(names(x = mean.exp) == rownames(x = sen_data1@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  sen_data1@meta.data$senescence <- mean.exp
}

Idents(sen_data1) <- sen_data1$annotation
min(sen_data1$senescence)
FeaturePlot(object = sen_data1, features = 'senescence', pt.size = .001, label = TRUE, label.size = 4, repel = TRUE) +
  theme(plot.title = element_blank(), text = element_text(size=6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), axis.title.x = element_blank(), axis.title.y = element_blank()) +
  scale_color_gradientn(colors=brewer.pal(n = 9, name = "YlGnBu"), breaks = c(.0, .1,.2,.3,.5,1),limits = c(.0,1.1), trans = scales::boxcox_trans(1))
ggsave(file = 'Human_sen_umap_blue.pdf', width=5, height=5, units="in")

Idents(sen_data1) <- sen_data1$annotation
sen_data1[["celltype"]] <- Idents(sen_data1)
new.cluster.ids <- c("Stem",  "Epithelial", "Epithelial", "Epithelial", "Epithelial", "Epithelial", "Epithelial", "Epithelial", "Epithelial")
names(new.cluster.ids) <- levels(sen_data1)
sen_data1 <- RenameIdents(sen_data1, new.cluster.ids)
sen_data1[["celltype"]] <- Idents(sen_data1)
sen_data1$senescence1

i= 'senescence'
selected_cells <- names(sen_data1$celltype)
vln_data <- FetchData(sen_data1,
                      vars = c(i,"Age", "celltype"),
                      cells = selected_cells,
                      slot = "data")
vln_data <- melt(vln_data)
vln_data <- vln_data[c('Age', 'celltype', 'value')]

#  my_comparisons <- list(c('WT_Control, WT_Arasco'),c('WT_Control, WT_Rev'))
All <- VlnPlot(sen_data1,group.by = 'celltype' , split.by = "Age", features = i, pt.size = 0, assay = "RNA", cols = c('#7CA1CC' ,'#FF4902'), log = FALSE, split.plot = TRUE) + 
  theme(axis.line = element_line(size = .3),axis.title.x = element_blank(), axis.ticks = element_line(size = .3), text = element_text(size=7), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6))  +
  stat_compare_means(data= vln_data, aes(x = celltype, y = value, fill = Age), hide.ns = TRUE, method = "t.test",  label = "p.signif", size = 1)  #+ scale_y_continuous(expand = expansion(mult = c(0, .1)))
dodge <- position_dodge(width = .9)
All <- All + geom_boxplot(width=0.3,outlier.shape = NA, coef = 0, lwd=.2, position = dodge) +
  stat_compare_means(data= vln_data, aes(x = celltype, y = value, fill = Age), hide.ns = TRUE, method = "t.test",  label = "p.signif", size = 1)  
All$layers[[1]]$aes_params$size = .15
All
ggsave(file = paste0('FIG1D.pdf'), plot=All, width=3.5, height=2.5, units="in")

i= 'senescence'
  selected_cells <- names(sen_data1$annotation)
  vln_data <- FetchData(sen_data1,
                        vars = c(i,"Age", "annotation"),
                        cells = selected_cells,
                        slot = "data")
  vln_data <- melt(vln_data)
  vln_data <- vln_data[c('Age', 'annotation', 'value')]
  
  #  my_comparisons <- list(c('WT_Control, WT_Arasco'),c('WT_Control, WT_Rev'))
  All <- VlnPlot(sen_data1,group.by = 'annotation' , split.by = "Age", features = i, pt.size = 0, assay = "RNA", cols = c('#7CA1CC' ,'#FF4902'), log = FALSE, split.plot = TRUE) + 
    theme(axis.line = element_line(size = .3),axis.title.x = element_blank(), axis.ticks = element_line(size = .3), text = element_text(size=7), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6))  +
    stat_compare_means(data= vln_data, aes(x = annotation, y = value, fill = Age), hide.ns = TRUE, method = "t.test",  label = "p.signif", size = 1)  #+ scale_y_continuous(expand = expansion(mult = c(0, .1)))
  dodge <- position_dodge(width = .9)
  All <- All + geom_boxplot(width=0.3,outlier.shape = NA, coef = 0, lwd=.2, position = dodge) +
    stat_compare_means(data= vln_data, aes(x = annotation, y = value, fill = Age), hide.ns = TRUE, method = "t.test",  label = "p.signif", size = 1)  
  All$layers[[1]]$aes_params$size = .15
  All
  ggsave(file = paste0('Cart_HUMAN_', i, '.pdf'), plot=All, width=3.5, height=2.5, units="in")

sen_data1$senescence <- as.numeric(sen_data1$senescence)
sen_data1$sen_cells <- 'non'
sen_data1 <- AddModuleScore(object = sen_data1, features = list(senescence), name = 'senescence')

# module score distribution
Senescence_module <- as.data.frame(sen_data1$senescence1)
modulescores <- Senescence_module %>%
  rownames_to_column(var="id") %>%
  pivot_longer(-id, names_to="celltype", values_to="score")


p <- ggplot(modulescores)
#p <- ggplot(onescore)
p + geom_point(aes(x=fct_inorder(id), y=sort(score))) +
  facet_wrap(~celltype) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())


sen_cells <-WhichCells(object = sen_data1, expression = senescence1 > .2)
DimPlot(sen_data1, label=T,group.by = 'annotation', cells.highlight= list(sen_cells),  cols.highlight = c("darkblue"),cols= "grey")

sen_data1$sen_cells[sen_cells] <- paste('senescent')

All_Genes <- sen_data1@assays$RNA@data@Dimnames[[1]]
MHC_score <- c('H2-D1','H2-K','H2-Q1','H2-Q2','H2-T23','H2-T3','Cd74','H2-Ac','H2-Ab1','H2-Q7')
MHC_score <- MHC_score %>% convert_mouse_to_human_symbols()
MHC_score <- intersect(All_Genes, MHC_score)

mean.exp <- zscore(colMeans(x = sen_data1@assays$RNA@data[MHC_score, ], na.rm = TRUE), dist = 'norm')
if (all(names(x = mean.exp) == rownames(x = sen_data1@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  sen_data1@meta.data$MHC_score <- mean.exp
}
sen_data1$sen_cells
selected_cells <- names(sen_data1$Age)
vln_data <- FetchData(sen_data1,
                      vars = c(i,"Age", "sen_cells"),
                      cells = selected_cells,
                      slot = "data")
vln_data <- melt(vln_data)
vln_data <- vln_data[c('Age', 'sen_cells', 'value')]

#  my_comparisons <- list(c('WT_Control, WT_Arasco'),c('WT_Control, WT_Rev'))
All <- VlnPlot(sen_data1,group.by = 'Age', split.by = "sen_cells", features = i, pt.size = 0, assay = "RNA", cols = c('#7CA1CC','#FF4902'), log = FALSE, split.plot = TRUE) + 
  theme(legend.position = 'none')  +
  theme(axis.line = element_line(size = .3),axis.title.x = element_blank(), axis.ticks = element_line(size = .3), text = element_text(size=7), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6))
dodge <- position_dodge(width = .9)
All <- All + geom_boxplot(width=0.3,outlier.shape = NA, coef = 0, lwd=.2, position = dodge) +
  stat_compare_means(data= vln_data, aes(x = Age, y = value, fill = sen_cells), hide.ns = TRUE, method = "t.test", size = 2)  
All$layers[[1]]$aes_params$size = .15
All
ggsave(file = paste0('Fig 4 D human.pdf'), plot=All, width=2.5, height=2.5, units="in")

Idents(sen_data1) <- sen_data1$sen_cells
DE_treat <- FindMarkers(sen_data1, ident.1 = "senescent", ident.2 = "non", test.use = "MAST", logfc.threshold = .1, min.pct = .05, assay = 'RNA')
EnhancedVolcano(DE_treat, lab = rownames(DE_treat), x = 'avg_log2FC', y = 'p_val_adj', title = 'uPAR vs UT', pCutoff = 10e-1, FCcutoff = 0.25, xlim = c(-1,1), ylim = c(0,320), subtitle = 'All Cells' )
DE_treat$gene_id <- rownames(DE_treat)
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
write.csv(gsea_output,'CarT_human_GSEA_sen_vs_non.csv')
colnames(gsea_output)
theTable$Position <- factor(theTable$Position, levels = c(...))
pathways <- c('GOCC_MHC_CLASS_I_PROTEIN_COMPLEX', 'GOBP_REGULATED_EXOCYTOSIS', 'GOBP_CELL_ADHESION', 'GOBP_WOUND_HEALING', 'GOBP_TISSUE_DEVELOPMENT', 'GOBP_CHAPERONE_MEDIATED_PROTEIN_FOLDING', 'GOMF_DNA_POLYMERASE_BINDING', 'GOMF_FATTY_ACID_BINDING', 'GOBP_MRNA_PROCESSING', 'GOBP_ATP_METABOLIC_PROCESS')
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
  scale_y_continuous(limits = c(-3,3),expand = expansion(mult = c(0, 0)), breaks = scales::breaks_extended(n = 3)) +
  labs(y= "Normalized Enrichment Score", x="Pathway") 
All + geom_vline(xintercept = 0)
ggsave(file = paste0('Fig 1 I_GSEA.pdf'), plot=All, width=5, height=3, units="in")
  
  
