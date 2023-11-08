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

sen_data1 <- NormalizeData(sen_data1)
all.genes <- rownames(sen_data1)
sen_data1 <- ScaleData(sen_data1, features = all.genes)
sen_data1 <- magic(sen_data1)

Idents(sen_data1) <- sen_data1$annotation
sen_data1 <- subset(sen_data1,  idents = c('Stem cells', 'TA', 'Enterocyte', 'Colonocyte', 'BEST4+ epithelial', 'BEST2+ Goblet cell', 'Goblet cell', 'Paneth', 'Tuft'))    
my_levels <- c('Stem cells', 'TA', 'Colonocyte', 'BEST4+ epithelial', 'BEST2+ Goblet cell', 'Goblet cell', 'Paneth', 'Tuft', 'Enterocyte')
sen_data1$annotation <- factor(x = sen_data1$annotation, levels = my_levels)
#new
senescence <- c('Psap','Timp2','Itm2b','Lgmn','Igf1','Arl6ip1','Ckb','Lsmem1','Igfbp4','Gsn','Zbtb20','Ccl2','Ccl7','Trio','Tulp4','Cxcl1') %>% convert_mouse_to_human_symbols()
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
) %>% convert_mouse_to_human_symbols()
Inflammatory_score_new <- c('Ccl4','Ccl5','Ccl9','Ccl20','Cd74','Cox4I1','Cox8A','Ctsb','Ctsd','H2-Ac','H2-Ab1','H2-Q7','Ifi47','Ifit1Bl1','Ifit2','Ifit3','Ifitm3','Irf7','Isg15','Stat2','Stat3','Tap1','Ucp2') %>% convert_mouse_to_human_symbols()
HALLMARK_INFLAMMATORY_RESPONSE_v7_5_1 <- read_delim("/Users/vyom/documents/HALLMARK_INFLAMMATORY_RESPONSE.v7.5.1.gmt", 
                                                    delim = "\t", escape_double = FALSE, 
                                                    col_names = FALSE, trim_ws = TRUE)
H_inflammatory <- as.character(HALLMARK_INFLAMMATORY_RESPONSE_v7_5_1[1,])



All_Genes <- sen_data1@assays$RNA@data@Dimnames[[1]]
senescence <- intersect(All_Genes, senescence)
sen_mayo <- intersect(All_Genes, sen_mayo)
Inflammatory_Response <- intersect(All_Genes, H_inflammatory)
Inflammatory_score_new <- intersect(All_Genes, Inflammatory_score_new)

mean.exp <- zscore(colMeans(x = sen_data1@assays$RNA@data[Inflammatory_Response, ], na.rm = TRUE), dist ='norm')
if (all(names(x = mean.exp) == rownames(x = sen_data1@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  sen_data1@meta.data$Inflammatory_Response <- mean.exp
}

mean.exp <- zscore(colMeans(x = sen_data1@assays$RNA@data[Inflammatory_score_new, ], na.rm = TRUE), dist ='norm')
if (all(names(x = mean.exp) == rownames(x = sen_data1@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  sen_data1@meta.data$Inflammatory_score_new <- mean.exp
}

mean.exp <- zscore(colMeans(x = sen_data1@assays$RNA@data[senescence, ], na.rm = TRUE), dist ='norm')
if (all(names(x = mean.exp) == rownames(x = sen_data1@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  sen_data1@meta.data$senescence <- mean.exp
}
mean.exp <- zscore(colMeans(x = sen_data1@assays$RNA@data[sen_mayo, ], na.rm = TRUE), dist ='norm')
if (all(names(x = mean.exp) == rownames(x = sen_data1@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  sen_data1@meta.data$senescence <- mean.exp
}

Idents(sen_data1) <- sen_data1$annotation
min(sen_data1$senescence)
FeaturePlot(object = sen_data1, features = 'sen_mayo', pt.size = .001, label = TRUE, label.size = 4, repel = TRUE) +
  theme(plot.title = element_blank(), text = element_text(size=6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), axis.title.x = element_blank(), axis.title.y = element_blank()) +
  scale_color_gradientn(colors=brewer.pal(n = 9, name = "YlGnBu"))
ggsave(file = 'Human_sen_umap_blue.pdf', width=5, height=5, units="in")

Idents(sen_data1) <- sen_data1$annotation
sen_data1[["celltype"]] <- Idents(sen_data1)
new.cluster.ids <- c("Stem",  "Epithelial", "Epithelial", "Epithelial", "Epithelial", "Epithelial", "Epithelial", "Epithelial", "Epithelial")
names(new.cluster.ids) <- levels(sen_data1)
sen_data1 <- RenameIdents(sen_data1, new.cluster.ids)
sen_data1[["celltype"]] <- Idents(sen_data1)


sen_data1$Inflammatory_score_new
i= 'sen_mayo'
i= 'Inflammatory_Response'
i= 'Inflammatory_score_new'
i= 'PLAUR'

selected_cells <- names(sen_data1$celltype)
vln_data <- FetchData(sen_data1,
                      vars = c(i,"Age", "celltype"),
                      cells = selected_cells,
                      slot = "data")
vln_data <- melt(vln_data)
vln_data <- vln_data[c('Age', 'celltype', 'value')]

All <- VlnPlot(sen_data1,group.by = 'celltype' , split.by = "Age", features = i, pt.size = 0, assay = "MAGIC_RNA", cols = c('#7CA1CC' ,'#FF4902'), log = FALSE, split.plot = TRUE) + 
  theme(axis.line = element_line(size = .3),axis.title.x = element_blank(), axis.ticks = element_line(size = .3), text = element_text(size=7), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6))  +
  stat_compare_means(data= vln_data, aes(x = celltype, y = value, fill = Age), hide.ns = TRUE, method = "wilcox.test",  label = "p.signif", size = 1)  #+ scale_y_continuous(expand = expansion(mult = c(0, .1)))
dodge <- position_dodge(width = .9)
All <- All + geom_boxplot(width=0.3,outlier.shape = NA, coef = 0, lwd=.2, position = dodge) +
  stat_compare_means(data= vln_data, aes(x = celltype, y = value, fill = Age), hide.ns = TRUE, method = "wilcox.test",  label = "p.signif", size = 1)  
All$layers[[1]]$aes_params$size = .15
All
ggsave(file = paste0('Human_',i,'_violin.pdf'), plot=All, width=3.5, height=2.5, units="in")

source("~/analysis/split_violin.R")

gene.list <- c('LGR4','LGR5', 'OLFM4',  'LIG1', 'TERT', 'CCND1' )
DefaultAssay(sen_data1) <- "MAGIC_RNA"
selected_cells <- names(sen_data1$annotation)
selected_cells <- names(sen_data1$annotation[sen_data1$annotation == c("Stem cells")])
vln_data <- FetchData(sen_data1,
                      vars = c(gene.list,"Age"),
                      cells = selected_cells,
                      slot = "data")

vln_data <- melt(vln_data)
All <- ggplot(vln_data, aes(x = variable, y = value, fill= Age)) + geom_split_violin(scale="width", trim = TRUE, size = .1) + theme_vyom + theme(legend.position = 'none') + 
  geom_boxplot(width=0.2,outlier.shape = NA, coef = 0, lwd=.2) + xlab('') + ylab('Expression Level')+ scale_fill_manual(values= c('#7CA1CC' ,'#FF4902')) +
  theme(text = element_text(size=7), axis.text.x = element_text(size = 6, angle = 45, hjust = 1, vjust = 1), axis.text.y = element_text(size = 6), axis.title.y = element_text(size = 7)) +
  stat_compare_means( method = "wilcox.test",  label = "p.signif", size = 2) + scale_y_continuous( expand = expansion(mult = c(0, 0.1)))
All
ggsave(file = 'HUMAN_control_Stem_violin_select.pdf', plot=All, width=2.75, height=2, units="in")


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

  DefaultAssay(sen_data1) <- 'RNA'
sen_data1$senescence <- as.numeric(sen_data1$senescence)
sen_data1$sen_cells <- 'non'
sen_data1 <- AddModuleScore(object = sen_data1, features = list(sen_mayo), name = 'senescence')

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


sen_cells <-WhichCells(object = sen_data1, expression = senescence1 > .15)
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

sen_data1$Age
#identify senescent cells
sen_data1$sen_cells <- 'non'
sen_data1 <- AddModuleScore(object = sen_data1, features = list(sen_mayo), name = 'sen_mayo')

# module score distribution
Senescence_module <- as.data.frame(sen_data1$sen_mayo1)
modulescores <- Senescence_module %>%
  rownames_to_column(var="id") %>%
  pivot_longer(-id, names_to="celltype", values_to="score")


p <- ggplot(modulescores)
#p <- ggplot(onescore)
inflection_plot <- p + geom_point(aes(x=fct_inorder(id), y=sort(score)))+ scale_y_continuous(breaks = seq(-.1, .3, by = .01)) +
  facet_wrap(~celltype) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
inflection_plot

sen_cells <-WhichCells(object = sen_data1, expression = sen_mayo1 > .12)
DimPlot(sen_data1, label=T,group.by = 'Treatment', cells.highlight= list(sen_cells),  cols.highlight = c("darkblue"),cols= "grey")

sen_data1$sen_cells[sen_cells] <- paste('senescent')

Prop_table<- prop.table(x = table(sen_data1$Age, sen_data1$sen_cells), margin = 2)
Prop_Table <- as.data.frame(Prop_table, row.names = NULL, optional = FALSE,
                            make.names = TRUE, stringsAsFactors = default.stringsAsFactors())
Prop_Table <- Prop_Table[Prop_Table$Var2 == 'senescent',]
melt(Prop_Table)
ggplot(Prop_Table, aes(fill=Var1, y=Freq, x=Var2)) + 
  geom_bar(position="fill", stat="identity")+theme_vyom +scale_fill_manual(values = c("#7CA1CC", "#FF4902"))
ggsave( "HUMAN_senescent_cell_proportion.pdf", width=2.5, height=3, units="in")

#subset to just senescent cells
Idents(sen_data1) <- sen_data1$sen_cells
senenscent_scrna <- subset(sen_data1,  idents = c('senescent'))

Prop_table<- prop.table(x = table(senenscent_scrna$Age, senenscent_scrna$annotation ), margin = 2)
Prop_Table <- as.data.frame(Prop_table, row.names = NULL, optional = FALSE,
                            make.names = TRUE, stringsAsFactors = default.stringsAsFactors())
Prop_Table1 <- Prop_Table
my_levels <- c('Stem cells', 'TA', 'Colonocyte', 'BEST4+ epithelial', 'BEST2+ Goblet cell', 'Goblet cell', 'Paneth', 'Tuft', 'Enterocyte')
my_levels1 <- c('25-30','70-75')

Prop_Table1$Var1 <- factor(Prop_Table1$Var1,levels = my_levels1)
Prop_Table1$Var2 <- factor(Prop_Table1$Var2,levels = my_levels)
plot <- ggplot(data = Prop_Table1, aes(Var2, Freq, fill=Var1)) + geom_bar(position="stack", stat="identity", na.rm = TRUE) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), panel.background = element_rect(fill = "white", colour = "Black")) +
  xlab("Cell Type") + ylab("Fraction of Cells") + labs(fill = "Treatment") + scale_fill_manual(values = c('#7CA1CC' ,'#FF4902'))
plot
ggsave(file = paste0('HUMAN_senescent_cell_proportion_percelltype.pdf'), width=4, height=3, units="in")




#proportions per cell type
Prop_table<- prop.table(x = table(sen_data1$annotation, sen_data1$Age), margin = 2)
Prop_Table <- as.data.frame(Prop_table, row.names = NULL, optional = FALSE,
                            make.names = TRUE, stringsAsFactors = default.stringsAsFactors())
Prop_Table1 <- Prop_Table
my_levels <- c('Stem cells', 'TA', 'Colonocyte', 'BEST4+ epithelial', 'BEST2+ Goblet cell', 'Goblet cell', 'Paneth', 'Tuft', 'Enterocyte')
my_levels1 <- c('25-30','70-75')

Prop_Table1$Var1 <- factor(Prop_Table1$Var1,levels = my_levels)
Prop_Table1$Var2 <- factor(Prop_Table1$Var2,levels = my_levels1)
plot <- ggplot(data = Prop_Table1, aes(Var1, Freq, fill=Var2)) + geom_bar(position="dodge", stat="identity", na.rm = TRUE) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), panel.background = element_rect(fill = "white", colour = "Black")) +
  xlab("Cell Type") + ylab("Fraction of Cells") + labs(fill = "Treatment") + scale_fill_manual(values = c('#7CA1CC' ,'#FF4902'))
plot
ggsave(file = paste0('Human_control_scrna_Proportions.pdf'), width=4, height=3, units="in")



  
