#assorted plots
Idents(cart_obj) <- cart_obj@meta.data$orig.ident
levels(factor(cart_obj@meta.data$orig.ident))

cart_obj[["Treatment_age"]] <- Idents(cart_obj)
new.cluster.ids <- c("uPAR_old", "uPAR_young", "uPAR_old", "uPAR_young", "UT_old", "UT_young", "UT_old", "UT_young")
names(new.cluster.ids) <- levels(cart_obj)
cart_obj <- RenameIdents(cart_obj, new.cluster.ids)
cart_obj[["Treatment_age"]] <- Idents(cart_obj)
All_Genes <- cart_obj@assays$RNA@data@Dimnames[[1]]

Idents(cart_obj) <- cart_obj@meta.data$Cell_Type
levels(factor(cart_obj@meta.data$Cell_Type))

cart_obj[["celltype"]] <- Idents(cart_obj)
new.cluster.ids <- c("Stem", "Stem", "Epithelial", "Epithelial", "Epithelial", "Epithelial", "Epithelial", "Epithelial", "Epithelial", "Immune", "Immune", "Immune")
names(new.cluster.ids) <- levels(cart_obj)
cart_obj <- RenameIdents(cart_obj, new.cluster.ids)
cart_obj[["celltype"]] <- Idents(cart_obj)
All_Genes <- cart_obj@assays$RNA@data@Dimnames[[1]]


my_levels <- c("UT_young", "uPAR_young", "UT_old", "uPAR_old")
cart_obj$Treatment_age <- factor(x = cart_obj$Treatment_age, levels = my_levels)
my_levels <- c("Stem", "Epithelial", "Immune")
cart_obj$celltype <- factor(x = cart_obj$celltype, levels = my_levels)

Idents(cart_obj) <- cart_obj$Age
cart_old <- subset(cart_obj,  idents = 'Old')
cart_young <- subset(cart_obj,  idents = 'Young')

stem_gene <- c('Lgr4', 'Myc', 'Sox9', 'Olfm4', 'Malat1', 'Hopx', 'Ccnd1')
DefaultAssay(cart_obj@assays$MAGIC_RNA) <- "MAGIC_RNA"
selected_cells <- names(cart_obj$Cell_Type[Cart_old$Cell_Type == c('Stem')])
vln_data <- FetchData(cart_obj,
                      vars = c(stem_gene,"Treatment_age"),
                      cells = selected_cells,
                      slot = "data")

vln_data <- melt(vln_data)

All <- ggplot(vln_data, aes(x = variable, y = value, fill= Treatment_age)) + geom_violin(position=position_dodge(width = 0.9),scale="width", trim = TRUE, size = .1) + theme_vyom + theme(legend.position = 'none') + 
  geom_boxplot(position=position_dodge(width = 0.9), width=0.2,outlier.shape = NA, coef = 0, lwd=.2) + xlab('') + ylab('Expression Level')  +
  theme(text = element_text(size=7), axis.text.x = element_text(size = 6, angle = 45, hjust = 1, vjust = 1), axis.text.y = element_text(size = 6), axis.title.y = element_text(size = 7)) +
  stat_compare_means( method = "wilcox.test",  label = "p.signif", size = 2) + scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
All
ggsave(file = 'Figure_3_old_young_combined.pdf', plot=All, width=4, height=2, units="in")



Idents(cart_obj) <- cart_obj$Age
cart_old <- subset(cart_obj,  idents = 'Old')
cart_young <- subset(cart_obj,  idents = 'Young')
Idents(Cart_old) <- Cart_old$Cell_Type
Cart_old <- subset(cart_old, idents = c('Enterocyte (Proximal)','Enterocyte (Distal)', 'Enteroendocrine', 'Goblet', 'Paneth'))  

my_levels <- c('Enterocyte (Proximal)','Enterocyte (Distal)', 'Enteroendocrine', 'Goblet', 'Paneth')
Cart_old$Cell_Type <- factor(x = Cart_old$Cell_Type, levels = my_levels)

Idents(cart_young) <- cart_young$Cell_Type
cart_young <- subset(cart_young, idents = c('Enterocyte (Proximal)','Enterocyte (Distal)', 'Enteroendocrine', 'Goblet', 'Paneth'))  

my_levels <- c('Enterocyte (Proximal)','Enterocyte (Distal)', 'Enteroendocrine', 'Goblet', 'Paneth')
cart_young$Cell_Type <- factor(x = cart_young$Cell_Type, levels = my_levels)

source("~/analysis/split_violin.R")
#Refined
genes_stem_prolif 
genes_antimicrobial    

Hormone <- c('Fabp1','Gast','Scgn', 'Pcsk1','Glp1r','Vegfa')
Absorption <- c('Cps1', 'Fabp1', 'Rbp2','Slc7a8', 'Vdr', 'Ppp1r14b')
antimicrobial <- c("Casp8","Itln1","Fabp2","Pdia4","Zg16","Muc2","Mmp7")
my_levels <- c('UT', 'uPAR')
Cart_old$Treatment <- factor(x = Cart_old$Treatment, levels = my_levels)
DefaultAssay(Cart_old@assays$MAGIC_RNA) <- "MAGIC_RNA"
selected_cells <- names(Cart_old$Cell_Type[Cart_old$Cell_Type == c('Enterocyte (Proximal)','Enterocyte (Distal)')])
vln_data <- FetchData(Cart_old,
                      vars = c(genes_stem_prolif,"Treatment"),
                      cells = selected_cells,
                      slot = "data")

vln_data <- melt(vln_data)

All <- ggplot(vln_data, aes(x = variable, y = value, fill= Treatment)) + geom_split_violin(scale="width", trim = TRUE, size = .1) + theme_vyom + theme(legend.position = 'none') + 
  geom_boxplot(width=0.2,outlier.shape = NA, coef = 0, lwd=.2) + xlab('') + ylab('Expression Level') + scale_fill_manual(values= c('#7CA1CC' ,'#FF4902')) +
  theme(text = element_text(size=7), axis.text.x = element_text(size = 6, angle = 45, hjust = 1, vjust = 1), axis.text.y = element_text(size = 6), axis.title.y = element_text(size = 7)) +
  stat_compare_means( method = "wilcox.test",  label = "p.signif", size = 2) + scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
All
ggsave(file = 'old_ent_anti_microbial.pdf', plot=All, width=4, height=2, units="in")

my_levels <- c('UT', 'uPAR' )
cart_young$Treatment <- factor(x = cart_young$Treatment, levels = my_levels)
DefaultAssay(cart_young@assays$MAGIC_RNA) <- "MAGIC_RNA"
selected_cells <- names(cart_young$Cell_Type[cart_young$Cell_Type == c('Enterocyte (Proximal)','Enterocyte (Distal)')])
vln_data <- FetchData(cart_young,
                      vars = c(genes_antimicrobial,"Treatment"),
                      cells = selected_cells,
                      slot = "data")

vln_data <- melt(vln_data)

All <- ggplot(vln_data, aes(x = variable, y = value, fill= Treatment)) + geom_split_violin(scale="width", trim = TRUE, size = .1) + theme_vyom + theme(legend.position = 'none') + 
  geom_boxplot(width=0.2,outlier.shape = NA, coef = 0, lwd=.2) + xlab('') + ylab('Expression Level') + scale_fill_manual(values= c('#7CA1CC' ,'#FF4902')) +
  theme(text = element_text(size=7), axis.text.x = element_text(size = 6, angle = 45, hjust = 1, vjust = 1), axis.text.y = element_text(size = 6), axis.title.y = element_text(size = 7)) +
  stat_compare_means( method = "wilcox.test",  label = "p.signif", size = 2) + scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
All
ggsave(file = 'Young_ent_nutrient_absorption.pdf', plot=All, width=4, height=2, units="in")


{
  genes_stem_prolif <- unique(c('EGF',	'WNT9B',	'WNT3',	'BAMBI')) %>% convert_human_to_mouse_symbols()
  genes_antimicrobial <-  unique(c('DEFA5',	'DEFA6',	'DEFA27',	'DEFA30',	'DEFA26',	'DEFA20',	'MPTX1',	'MPTX2',	'LYZ1',	'LYZ2',	'KLF6',	'KLF4',	'KLF15',	'PLA2G2A',	'REG3B',	'REG3G',	'ITLN1',	'SPINK4',	'MT1',	'FABP1',	'FABP2',	'FABP4',	'MMP7',	'PLA2G12A', 'MUC2',	'CLCA1',	'PDIA4',	'PDIA5',	'ZG16',	'REG4',	'REG3B',	'ANG4', 'CASP6',	'CASP8',	'REG1',	'REG1',	'REG3A',	'REG3B',	'REG4'))%>% convert_human_to_mouse_symbols()         
  genes_hormone_secretion <-  unique(c('CHGA',	'GCG',	'GHRL',	'TRPA1',	'NTS',	'AFP',	'VEGFA',	'PCSK1',	'SCN3A',	'NEUROG3',	'NEUROD1',	'NKX2-2',	'FABP5',	'FABP1',	'FABP6',	'GIP',	'GLP1R',	'GLP2R',	'PYY',	'GAST',	'SCGN'))%>% convert_human_to_mouse_symbols()
  genes_nutrient_absorption <-  unique(c('SLC7A8',	'FABP1',	'CPS1',	'VDR',	'APOA4',	'RBP2',	'GSTA1',	'SLC5A1',	'SLC7A9',	'EHF',	'PPP1R14B',	'PP1R14D')) %>% convert_human_to_mouse_symbols()
  
  genes_stem_prolif <- genes_stem_prolif[!is.na(genes_stem_prolif)]
  genes_antimicrobial <- genes_antimicrobial[!is.na(genes_antimicrobial)]
  genes_hormone_secretion <- genes_hormone_secretion[!is.na(genes_hormone_secretion)]
  genes_nutrient_absorption <- genes_nutrient_absorption[!is.na(genes_nutrient_absorption)]
  
  All_Genes <- cart_old@assays$RNA@data@Dimnames[[1]]
  genes_stem_prolif <- intersect(All_Genes, genes_stem_prolif)
  genes_antimicrobial <- intersect(All_Genes, genes_antimicrobial)   
  genes_hormone_secretion <- intersect(All_Genes, genes_hormone_secretion)  
  genes_nutrient_absorption <- intersect(All_Genes, genes_nutrient_absorption)
  
  Cart_old <- AddModuleScore(object = Cart_old, features = list(genes_stem_prolif), name = 'stem_prolif')
  Cart_old <- AddModuleScore(object = Cart_old, features = list(genes_antimicrobial), name = 'antimicrobial')
  Cart_old <- AddModuleScore(object = Cart_old, features = list(genes_hormone_secretion), name = 'hormone_secretion')
  Cart_old <- AddModuleScore(object = Cart_old, features = list(genes_nutrient_absorption), name = 'nutrient_absorption')
  
  Idents(Cart_old) <- Cart_old$Cell_Type
  Cell_Types <- levels(Cart_old$Cell_Type)
  stem_prolif_FC <- list()
  nutrient_absorption_FC <- list()
  antimicrobial_FC <- list()
  hormone_secretion_FC <- list()
  
  for(i in Cell_Types){
    stem_prolif_FC[[i]] <- FoldChange(Cart_old,ident.1 = "uPAR", ident.2 = "UT",features = genes_stem_prolif, group.by = 'Treatment', subset.ident = i)
  }
  stem_prolif_FC_aggr <- c()
  for(i in Cell_Types){
    stem_prolif_FC_aggr <- c(stem_prolif_FC_aggr, mean(as.numeric(stem_prolif_FC[[i]]$avg_log2FC)))
  }
  
  for(i in Cell_Types){
    nutrient_absorption_FC[[i]] <- FoldChange(Cart_old,ident.1 = "uPAR", ident.2 = "UT",features = genes_nutrient_absorption, group.by = 'Treatment', subset.ident = i)
  }
  nutrient_absorption_FC_aggr <- c()
  for(i in Cell_Types){
    nutrient_absorption_FC_aggr <- c(nutrient_absorption_FC_aggr, mean(as.numeric(nutrient_absorption_FC[[i]]$avg_log2FC)))
  }
  
  for(i in Cell_Types){
    antimicrobial_FC[[i]] <- FoldChange(Cart_old,ident.1 = "uPAR", ident.2 = "UT",features = genes_antimicrobial, group.by = 'Treatment', subset.ident = i)
  }
  antimicrobial_FC_aggr <- c()
  for(i in Cell_Types){
    antimicrobial_FC_aggr <- c(antimicrobial_FC_aggr, mean(as.numeric(antimicrobial_FC[[i]]$avg_log2FC)))
  }
  
  for(i in Cell_Types){
    hormone_secretion_FC[[i]] <- FoldChange(Cart_old,ident.1 = "uPAR", ident.2 = "UT",features = genes_hormone_secretion, group.by = 'Treatment', subset.ident = i)
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
  Idents(Cart_old) <- Cart_old$Cell_Type
  pvals <- c()
  for(j in scores){
    for(i in Cell_Types){
      selected_cells <- WhichCells(object = Cart_old ,idents = c(i))
      vln_data <- FetchData(Cart_old,
                            vars = c(j,"Treatment", "Cell_Type"),
                            cells = selected_cells,
                            slot = "data")
      vln_data <- melt(vln_data)
      vln_data <- vln_data[c('Treatment', 'Cell_Type', 'value')]
      
      pvals <- c(pvals, t.test(vln_data[which(vln_data$Treatment == 'uPAR'),]$value, vln_data[which(vln_data$Treatment == 'UT'),]$value,  alternative = c("greater"))$p.value)
      print(i)
    }
  }
  
  Module_Heatmap1$pval <- -log(pvals) 
  
  my_levels <- c('Enterocyte (Proximal)','Enterocyte (Distal)', 'Enteroendocrine', 'Goblet', 'Paneth')
  Module_Heatmap1$Cell_Type <- factor(Module_Heatmap1$Cell_Type, levels = my_levels)
  
  pathwayColorsDiff = rev(brewer.pal(9, "RdBu"))
  ggplot(Module_Heatmap1, aes(Cell_Type, variable, fill= value)) + geom_tile() + 
    scale_fill_gradientn(name = "Log2FC", colors = pathwayColorsDiff, limits = c(-.25,.25)) + coord_flip() +
    theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45), axis.text.y =   element_text( vjust = 0.1, hjust = 1), panel.background = element_blank(), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())
  
  ggsave(file="Old_CART_score_heatmap.pdf",  width=3, height=3, units="in")
  
  ggplot(Module_Heatmap1, aes(x = Cell_Type, y = variable, size = pval ,color= value)) + geom_point() + theme_vyom + 
    scale_color_gradientn(name = "Log2FC", colors = pathwayColorsDiff, limits = c(-.25,.25)) + coord_flip() + labs(size="-log(p)") +
    theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45), axis.text.y =   element_text( vjust = 0.1, hjust = 1), panel.background = element_blank(), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())
  
  ggsave(file="Old_CART_score_dotplot.pdf",  width=4, height=5, units="in")
}


{
  genes_stem_prolif <- unique(c('EGF',	'WNT9B',	'WNT3',	'BAMBI')) %>% convert_human_to_mouse_symbols()
  genes_antimicrobial <-  unique(c('DEFA5',	'DEFA6',	'DEFA27',	'DEFA30',	'DEFA26',	'DEFA20',	'MPTX1',	'MPTX2',	'LYZ1',	'LYZ2',	'KLF6',	'KLF4',	'KLF15',	'PLA2G2A',	'REG3B',	'REG3G',	'ITLN1',	'SPINK4',	'MT1',	'FABP1',	'FABP2',	'FABP4',	'MMP7',	'PLA2G12A', 'MUC2',	'CLCA1',	'PDIA4',	'PDIA5',	'ZG16',	'REG4',	'REG3B',	'ANG4', 'CASP6',	'CASP8',	'REG1',	'REG1',	'REG3A',	'REG3B',	'REG4'))%>% convert_human_to_mouse_symbols()         
  genes_hormone_secretion <-  unique(c('CHGA',	'GCG',	'GHRL',	'TRPA1',	'NTS',	'AFP',	'VEGFA',	'PCSK1',	'SCN3A',	'NEUROG3',	'NEUROD1',	'NKX2-2',	'FABP5',	'FABP1',	'FABP6',	'GIP',	'GLP1R',	'GLP2R',	'PYY',	'GAST',	'SCGN'))%>% convert_human_to_mouse_symbols()
  genes_nutrient_absorption <-  unique(c('SLC7A8',	'FABP1',	'CPS1',	'VDR',	'APOA4',	'RBP2',	'GSTA1',	'SLC5A1',	'SLC7A9',	'EHF',	'PPP1R14B',	'PP1R14D')) %>% convert_human_to_mouse_symbols()
  
  genes_stem_prolif <- genes_stem_prolif[!is.na(genes_stem_prolif)]
  genes_antimicrobial <- genes_antimicrobial[!is.na(genes_antimicrobial)]
  genes_hormone_secretion <- genes_hormone_secretion[!is.na(genes_hormone_secretion)]
  genes_nutrient_absorption <- genes_nutrient_absorption[!is.na(genes_nutrient_absorption)]
  
  All_Genes <- cart_young@assays$RNA@data@Dimnames[[1]]
  genes_stem_prolif <- intersect(All_Genes, genes_stem_prolif)
  genes_antimicrobial <- intersect(All_Genes, genes_antimicrobial)   
  genes_hormone_secretion <- intersect(All_Genes, genes_hormone_secretion)  
  genes_nutrient_absorption <- intersect(All_Genes, genes_nutrient_absorption)
  
  cart_young <- AddModuleScore(object = cart_young, features = list(genes_stem_prolif), name = 'stem_prolif')
  cart_young <- AddModuleScore(object = cart_young, features = list(genes_antimicrobial), name = 'antimicrobial')
  cart_young <- AddModuleScore(object = cart_young, features = list(genes_hormone_secretion), name = 'hormone_secretion')
  cart_young <- AddModuleScore(object = cart_young, features = list(genes_nutrient_absorption), name = 'nutrient_absorption')
  
  Idents(cart_young) <- cart_young$Cell_Type
  Cell_Types <- levels(cart_young$Cell_Type)
  stem_prolif_FC <- list()
  nutrient_absorption_FC <- list()
  antimicrobial_FC <- list()
  hormone_secretion_FC <- list()
  
  for(i in Cell_Types){
    stem_prolif_FC[[i]] <- FoldChange(cart_young,ident.1 = "uPAR", ident.2 = "UT",features = genes_stem_prolif, group.by = 'Treatment', subset.ident = i)
  }
  stem_prolif_FC_aggr <- c()
  for(i in Cell_Types){
    stem_prolif_FC_aggr <- c(stem_prolif_FC_aggr, mean(as.numeric(stem_prolif_FC[[i]]$avg_log2FC)))
  }
  
  for(i in Cell_Types){
    nutrient_absorption_FC[[i]] <- FoldChange(cart_young,ident.1 = "uPAR", ident.2 = "UT",features = genes_nutrient_absorption, group.by = 'Treatment', subset.ident = i)
  }
  nutrient_absorption_FC_aggr <- c()
  for(i in Cell_Types){
    nutrient_absorption_FC_aggr <- c(nutrient_absorption_FC_aggr, mean(as.numeric(nutrient_absorption_FC[[i]]$avg_log2FC)))
  }
  
  for(i in Cell_Types){
    antimicrobial_FC[[i]] <- FoldChange(cart_young,ident.1 = "uPAR", ident.2 = "UT",features = genes_antimicrobial, group.by = 'Treatment', subset.ident = i)
  }
  antimicrobial_FC_aggr <- c()
  for(i in Cell_Types){
    antimicrobial_FC_aggr <- c(antimicrobial_FC_aggr, mean(as.numeric(antimicrobial_FC[[i]]$avg_log2FC)))
  }
  
  for(i in Cell_Types){
    hormone_secretion_FC[[i]] <- FoldChange(cart_young,ident.1 = "uPAR", ident.2 = "UT",features = genes_hormone_secretion, group.by = 'Treatment', subset.ident = i)
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
  Idents(cart_young) <- cart_young$Cell_Type
  pvals <- c()
  for(j in scores){
    for(i in Cell_Types){
      selected_cells <- WhichCells(object = cart_young ,idents = c(i))
      vln_data <- FetchData(cart_young,
                            vars = c(j,"Treatment", "Cell_Type"),
                            cells = selected_cells,
                            slot = "data")
      vln_data <- melt(vln_data)
      vln_data <- vln_data[c('Treatment', 'Cell_Type', 'value')]
      
      pvals <- c(pvals, t.test(vln_data[which(vln_data$Treatment == 'uPAR'),]$value, vln_data[which(vln_data$Treatment == 'UT'),]$value,  alternative = c("greater"))$p.value)
      print(i)
    }
  }
  
  Module_Heatmap1$pval <- -log(pvals) 
  
  my_levels <- c('Enterocyte (Proximal)','Enterocyte (Distal)', 'Enteroendocrine', 'Goblet', 'Paneth')
  Module_Heatmap1$Cell_Type <- factor(Module_Heatmap1$Cell_Type, levels = my_levels)
  
  pathwayColorsDiff = rev(brewer.pal(9, "RdBu"))
  ggplot(Module_Heatmap1, aes(Cell_Type, variable, fill= value)) + geom_tile() + 
    scale_fill_gradientn(name = "Log2FC", colors = pathwayColorsDiff, limits = c(-.25,.25)) + coord_flip() +
    theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45), axis.text.y =   element_text( vjust = 0.1, hjust = 1), panel.background = element_blank(), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())
  
  ggsave(file="young_CART_score_heatmap.pdf",  width=3, height=3, units="in")
  
  ggplot(Module_Heatmap1, aes(x = Cell_Type, y = variable, size = pval ,color= value)) + geom_point() + theme_vyom + 
    scale_color_gradientn(name = "Log2FC", colors = pathwayColorsDiff, limits = c(-.25,.25)) + coord_flip() + labs(size="-log(p)") +
    theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45), axis.text.y =   element_text( vjust = 0.1, hjust = 1), panel.background = element_blank(), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())
  
  ggsave(file="young_CART_score_dotplot.pdf",  width=4, height=5, units="in")
}

{
  Gene_list <- c('H2-Aa', 'H2-Ab1', 'H2-DMb2', 'H2-Q7', 'Cd74', 'Tap1', 'Ifngr1', 'Ifitm3', 'Ifitm2', 'Ccl5', 'Ccl4', 'Ccl3', 'Il16', 'Il18')
  
  All_Genes <- cart_old@assays$RNA@data@Dimnames[[1]]
  Gene_list <- intersect(All_Genes, Gene_list)
  
  Idents(cart_old) <- cart_old$celltype
  Cell_Types <- levels(cart_old$celltype)

  
  DE_tControlle_all = Idents(cart_old) %>% levels() %>% intersect(Cell_Types) %>% lapply(get_lfc_celltype, seurat_obj = cart_old, condition_colname = "Treatment", condition_oi = 'uPAR', condition_reference = 'UT', celltype_col = NULL, expression_pct = 0.01) %>% purrr::reduce(full_join)    
  DE_tControlle_all[is.na(DE_tControlle_all)] = 0
  
  DE_tControlle_all_frame <- data.frame(DE_tControlle_all)
  rownames(DE_tControlle_all_frame) <- DE_tControlle_all_frame$gene
  Gene_list1 <- c('H2.Aa', 'H2.Ab1', 'H2.DMb2', 'H2.Q7', 'Cd74', 'Tap1', 'Ifngr1', 'Ifitm3', 'Ifitm2', 'Ccl5', 'Ccl4', 'Ccl3', 'Il16', 'Il18')
  
  DE_tControlle_all_filtered <- DE_tControlle_all_frame[Gene_list,]
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
    scale_fill_gradientn(name = "Log2FC", colors = pathwayColorsDiff, limits = c(-.65,.65)) + coord_flip() +
    theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45), axis.text.y =   element_text( vjust = 0.1, hjust = 1), panel.background = element_blank(), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())
  
  Idents(cart_old) <- cart_old$celltype
  pvals <- c()
  for(j in Gene_list){
    for(i in Cell_Types){
      selected_cells <- WhichCells(object = cart_old ,idents = c(i))
      vln_data <- FetchData(cart_old,
                            vars = c(j,"Treatment", "celltype"),
                            cells = selected_cells,
                            slot = "data")
      vln_data <- melt(vln_data)
      vln_data <- vln_data[c('Treatment', 'celltype', 'value')]
      
      pvals <- c(pvals, t.test(vln_data[which(vln_data$Treatment == 'uPAR'),]$value, vln_data[which(vln_data$Treatment == 'UT'),]$value,  alternative = c("greater"))$p.value)
      print(i)
    }
  }
  
  vis_lfc2$pval <- -log(pvals) 
  
  my_levels <- c(Gene_list1)
  vis_lfc2$gene <- factor(vis_lfc2$gene, levels = my_levels)
  
  pathwayColorsDiff = rev(brewer.pal(9, "RdBu"))
  ggplot(vis_lfc2, aes(variable, gene, fill= value)) + geom_tile() + 
    scale_fill_gradientn(name = "Log2FC", colors = pathwayColorsDiff, limits = c(-.65,.65)) + coord_flip() +
    theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45), axis.text.y =   element_text( vjust = 0.1, hjust = 1), panel.background = element_blank(), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())
  
  ggsave(file="MHC_Old_heatmap.pdf",  width=3, height=3, units="in")
  
  ggplot(vis_lfc2, aes(x = gene, y = variable, size = pval ,color= value)) + geom_point() + theme_vyom + 
    scale_color_gradientn(name = "Log2FC", colors = pathwayColorsDiff, limits = c(-.65,.65)) + coord_flip() + labs(size="-log(p)") +
    theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45), axis.text.y =   element_text( vjust = 0.1, hjust = 1), panel.background = element_blank(), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())
  
  ggsave(file="MHC_Old_dotplot.pdf",  width=3.5, height=4, units="in")
}
{
  Gene_list <- c('H2-Aa', 'H2-Ab1', 'H2-DMb2', 'H2-Q7', 'Cd74', 'Tap1', 'Ifngr1', 'Ifitm3', 'Ifitm2', 'Ccl5', 'Ccl4', 'Ccl3', 'Il16', 'Il18')
  
  All_Genes <- cart_young@assays$RNA@data@Dimnames[[1]]
  Gene_list <- intersect(All_Genes, Gene_list)
  
  Idents(cart_young) <- cart_young$Cell_Type
  Cell_Types <- levels(cart_young$Cell_Type)
  

  DE_tControlle_all = Idents(cart_young) %>% levels() %>% intersect(Cell_Types) %>% lapply(get_lfc_celltype, seurat_obj = cart_young, condition_colname = "Treatment", condition_oi = 'uPAR', condition_reference = 'UT', celltype_col = NULL, expression_pct = 0.01) %>% purrr::reduce(full_join)    
  DE_tControlle_all[is.na(DE_tControlle_all)] = 0
  
  DE_tControlle_all_frame <- data.frame(DE_tControlle_all)
  rownames(DE_tControlle_all_frame) <- DE_tControlle_all_frame$gene
  Gene_list1 <- c('H2.Aa', 'H2.Ab1', 'H2.DMb2', 'H2.Q7', 'Cd74', 'Tap1', 'Ifngr1', 'Ifitm3', 'Ifitm2', 'Ccl5', 'Ccl4', 'Ccl3', 'Il16', 'Il18')
  
  DE_tControlle_all_filtered <- DE_tControlle_all_frame[Gene_list,]
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
    scale_fill_gradientn(name = "Log2FC", colors = pathwayColorsDiff, limits = c(-1,1)) + coord_flip() +
    theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45), axis.text.y =   element_text( vjust = 0.1, hjust = 1), panel.background = element_blank(), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())
  
  Idents(cart_young) <- cart_young$Cell_Type
  pvals <- c()
  for(j in Gene_list){
    for(i in Cell_Types){
      selected_cells <- WhichCells(object = cart_young ,idents = c(i))
      vln_data <- FetchData(cart_young,
                            vars = c(j,"Treatment", "Cell_Type"),
                            cells = selected_cells,
                            slot = "data")
      vln_data <- melt(vln_data)
      vln_data <- vln_data[c('Treatment', 'Cell_Type', 'value')]
      
      pvals <- c(pvals, t.test(vln_data[which(vln_data$Treatment == 'uPAR'),]$value, vln_data[which(vln_data$Treatment == 'UT'),]$value,  alternative = c("greater"))$p.value)
      print(i)
    }
  }
  
  vis_lfc2$pval <- -log(pvals) 
  
  my_levels <- c(Gene_list1)
  vis_lfc2$gene <- factor(vis_lfc2$gene, levels = my_levels)
  
  pathwayColorsDiff = rev(brewer.pal(9, "RdBu"))
  ggplot(vis_lfc2, aes(variable, gene, fill= value)) + geom_tile() + 
    scale_fill_gradientn(name = "Log2FC", colors = pathwayColorsDiff, limits = c(-1,1)) + coord_flip() +
    theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45), axis.text.y =   element_text( vjust = 0.1, hjust = 1), panel.background = element_blank(), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())
  
  ggsave(file="supp4 B heatmap.pdf",  width=4.5, height=3, units="in")
  
  ggplot(vis_lfc2, aes(x = gene, y = variable, size = pval ,color= value)) + geom_point() + theme_vyom + 
    scale_color_gradientn(name = "Log2FC", colors = pathwayColorsDiff, limits = c(-1,1)) + coord_flip() + labs(size="-log(p)") +
    theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45), axis.text.y =   element_text( vjust = 0.1, hjust = 1), panel.background = element_blank(), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())
  
  ggsave(file="supp4 B dotplot.pdf",  width=4.5, height=4, units="in")
}
#new inflammatory score
Inflammatory_score_new <- c('Ccl4',	'Ccl5',	'Ccl9',	'Ccl20',	'Cd74',	'Cox4I1',	'Cox8A',	'Ctsb',	'Ctsd',	'H2-Ac',	'H2-Ab1',	'H2-Q7',	'Ifi47',	'Ifit1Bl1',	'Ifit2',	'Ifit3',	'Ifitm3',	'Irf7',	'Isg15',	'Stat2',	'Stat3',	'Tap1',	'Ucp2')
Inflammatory_Response <- intersect(All_Genes, Inflammatory_score_new)

mean.exp <- zscore(colMeans(x = cart_old@assays$MAGIC_RNA@data[Inflammatory_Response, ], na.rm = TRUE), dist = 'norm')
if (all(names(x = mean.exp) == rownames(x = cart_old@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  cart_old@meta.data$Inflammatory_Response <- mean.exp
}

mean.exp <- zscore(colMeans(x = cart_young@assays$MAGIC_RNA@data[Inflammatory_Response, ], na.rm = TRUE), dist = 'norm')
if (all(names(x = mean.exp) == rownames(x = cart_young@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  cart_young@meta.data$Inflammatory_Response <- mean.exp
}

mean.exp <- zscore(colMeans(x = cart_obj@assays$MAGIC_RNA@data[Inflammatory_Response, ], na.rm = TRUE), dist = 'norm')
if (all(names(x = mean.exp) == rownames(x = cart_obj@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  cart_obj@meta.data$Inflammatory_Response <- mean.exp
}

senescence <- c('Psap','Timp2','Itm2b','Lgmn','Igf1','Arl6ip1','Ckb','Lsmem1','Igfbp4','Gsn','Zbtb20','Ccl2','Ccl7','Trio','Tulp4','Cxcl1') %>% convert_mouse_to_human_symbols()

senescence <- intersect(All_Genes, senescence)
mean.exp <- zscore(colMeans(x = cart_old@assays$MAGIC_RNA@data[senescence, ], na.rm = TRUE), dist = 'norm')
if (all(names(x = mean.exp) == rownames(x = cart_old@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  cart_old@meta.data$senescence <- mean.exp
}

mean.exp <- zscore(colMeans(x = cart_young@assays$MAGIC_RNA@data[senescence, ], na.rm = TRUE), dist = 'norm')
if (all(names(x = mean.exp) == rownames(x = cart_young@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  cart_young@meta.data$senescence <- mean.exp
}

mean.exp <- zscore(colMeans(x = cart_obj@assays$MAGIC_RNA@data[senescence, ], na.rm = TRUE), dist = 'norm')
if (all(names(x = mean.exp) == rownames(x = cart_obj@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  cart_obj@meta.data$senescence <- mean.exp
}

MHC_score <- c('H2-D1','H2-K','H2-Q1','H2-Q2','H2-T23','H2-T3','Cd74','H2-Ac','H2-Ab1','H2-Q7')
MHC_score <- intersect(All_Genes, MHC_score)

mean.exp <- zscore(colMeans(x = cart_old@assays$MAGIC_RNA@data[MHC_score, ], na.rm = TRUE), dist = 'norm')
if (all(names(x = mean.exp) == rownames(x = cart_old@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  cart_old@meta.data$MHC_score <- mean.exp
}

mean.exp <- zscore(colMeans(x = cart_young@assays$MAGIC_RNA@data[MHC_score, ], na.rm = TRUE), dist = 'norm')
if (all(names(x = mean.exp) == rownames(x = cart_young@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  cart_young@meta.data$MHC_score <- mean.exp
}

mean.exp <- zscore(colMeans(x = cart_obj@assays$MAGIC_RNA@data[MHC_score, ], na.rm = TRUE), dist = 'norm')
if (all(names(x = mean.exp) == rownames(x = cart_obj@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  cart_obj@meta.data$MHC_score <- mean.exp
}
cart_old$Inflammatory_Response
i = 'Inflammatory_Response'
for(i in gene.list) {
  selected_cells <- names(cart_old$celltype)
  vln_data <- FetchData(cart_old,
                        vars = c(i,"Treatment", "celltype"),
                        cells = selected_cells,
                        slot = "data")
  vln_data <- melt(vln_data)
  vln_data <- vln_data[c('Treatment', 'celltype', 'value')]
  
  All <- VlnPlot(cart_old,group.by = 'celltype' , split.by = "Treatment", features = i, pt.size = 0, assay = "MAGIC_RNA", cols = c('#7CA1CC' ,'#FF4902'), log = FALSE, split.plot = TRUE) + 
    theme(legend.position = 'none')  +
    theme(axis.line = element_line(size = .3),axis.title.x = element_blank(), axis.ticks = element_line(size = .3), text = element_text(size=7), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6))  +
    stat_compare_means(data= vln_data, aes(x = celltype, y = value, fill = Treatment), hide.ns = TRUE, method = "t.test",  label = "p.signif", size = 1)  #+ scale_y_continuous(expand = expansion(mult = c(0, .1)))
  dodge <- position_dodge(width = .9)
  All <- All + geom_boxplot(width=0.3,outlier.shape = NA, coef = 0, lwd=.2, position = dodge) +
    stat_compare_means(data= vln_data, aes(x = celltype, y = value, fill = Treatment), hide.ns = TRUE, method = "t.test",  label = "p.signif", size = 1)  
  All$layers[[1]]$aes_params$size = .15
  ggsave(file = paste0('Cart_Old_4 b', i, '.pdf'), plot=All, width=2.5, height=2.5, units="in")
}

for(i in gene.list) {
  selected_cells <- names(cart_young$celltype)
  vln_data <- FetchData(cart_young,
                        vars = c(i,"Treatment", "celltype"),
                        cells = selected_cells,
                        slot = "data")
  vln_data <- melt(vln_data)
  vln_data <- vln_data[c('Treatment', 'celltype', 'value')]
  
  #  my_comparisons <- list(c('WT_Control, WT_Arasco'),c('WT_Control, WT_Rev'))
  All <- VlnPlot(cart_young,group.by = 'celltype' , split.by = "Treatment", features = i, pt.size = 0, assay = "MAGIC_RNA", cols = c('#7CA1CC' ,'#FF4902'), log = FALSE, split.plot = TRUE) + 
    theme(legend.position = 'none')  +
    theme(axis.line = element_line(size = .3),axis.title.x = element_blank(), axis.ticks = element_line(size = .3), text = element_text(size=7), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6))  +
    stat_compare_means(data= vln_data, aes(x = celltype, y = value, fill = Treatment), hide.ns = TRUE, method = "t.test",  label = "p.signif", size = 1)  #+ scale_y_continuous(expand = expansion(mult = c(0, .1)))
  dodge <- position_dodge(width = .9)
  All <- All + geom_boxplot(width=0.3,outlier.shape = NA, coef = 0, lwd=.2, position = dodge) +
    stat_compare_means(data= vln_data, aes(x = celltype, y = value, fill = Treatment), hide.ns = TRUE, method = "t.test",  label = "p.signif", size = 1)  
  All$layers[[1]]$aes_params$size = .15
  ggsave(file = paste0('Cart_young_Supp 4_', i, '.pdf'), plot=All, width=2.5, height=2.5, units="in")
}
levels(cart_obj$Treatment_age)
  for(i in gene.list) {
    selected_cells <- names(cart_obj$Cell_Type)
    vln_data <- FetchData(cart_obj,
                          vars = c(i,"Treatment_age", "Cell_Type"),
                          cells = selected_cells,
                          slot = "data")
    vln_data <- melt(vln_data)
    vln_data <- vln_data[c('Treatment_age', 'Cell_Type', 'value')]
    
    #  my_comparisons <- list(c('WT_Control, WT_Arasco'),c('WT_Control, WT_Rev'))
    All <- VlnPlot(cart_obj,group.by = 'Cell_Type' , split.by = "Treatment_age", features = i, pt.size = 0, assay = "MAGIC_RNA",  log = FALSE, split.plot = TRUE) + 
      theme(legend.position = 'none')  +
      theme(axis.line = element_line(size = .3),axis.title.x = element_blank(), axis.ticks = element_line(size = .3), text = element_text(size=7), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6))  +
      stat_compare_means(data= vln_data, aes(x = Cell_Type, y = value, fill = Treatment_age), hide.ns = TRUE, method = "t.test",  label = "p.signif", size = 1)  #+ scale_y_continuous(expand = expansion(mult = c(0, .1)))
    dodge <- position_dodge(width = .9)
    All <- All + geom_boxplot(width=0.3,outlier.shape = NA, coef = 0, lwd=.2, position = dodge) +
      stat_compare_means(data= vln_data, aes(x = Cell_Type, y = value, fill = Treatment_age), hide.ns = TRUE, method = "t.test",  label = "p.signif", size = 1)  
    All$layers[[1]]$aes_params$size = .15
    ggsave(file = paste0('Cart_all_', i, '.pdf'), plot=All, width=3.5, height=2.5, units="in")
  }



All <- VlnPlot(cart_obj,group.by = 'Cell_Type' , split.by = "orig.ident", features = i, pt.size = 0, assay = "MAGIC_RNA", log = FALSE, split.plot = TRUE) + 
  theme(axis.line = element_line(size = .3),axis.title.x = element_blank(), axis.ticks = element_line(size = .3), text = element_text(size=7), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6))  +
  stat_compare_means(data= vln_data, aes(x = Cell_Type, y = value, fill = orig.ident), hide.ns = TRUE, method = "t.test",  label = "p.signif", size = 1)  #+ scale_y_continuous(expand = expansion(mult = c(0, .1)))
dodge <- position_dodge(width = .9)
All
ggsave(file = paste0('Cart_all_label', i, '.pdf'), plot=All, width=3, height=2.5, units="in")

cart_young
Prop_table<- prop.table(x = table(cart_young$Cell_Type, cart_young$orig.ident), margin = 2)
Prop_Table <- as.data.frame(Prop_table, row.names = NULL, optional = FALSE,
                            make.names = TRUE, stringsAsFactors = default.stringsAsFactors())
levels(as.factor(cart_young$orig.ident))
Prop_Table1 <- Prop_Table%>%mutate(Var2=recode(Var2,uPAR_F_young="uPAR",uPAR_M_young="uPAR", UT_F_young="UT", UT_M_young="UT"))

celltype_sample_norm3 = summarySE(Prop_Table1, measurevar="Freq", groupvars=c("Var1","Var2"))

my_levels <- c('Stem', 'Transit Amplifying', 'Enterocyte Progenitor','Enterocyte (Proximal)','Enterocyte (Distal)', 'Enteroendocrine', 'Goblet', 'Paneth','Tuft', 'T Cell', 'B Cell', 'Myeloid')
celltype_sample_norm3$Var1 = factor(celltype_sample_norm3$Var1, levels = my_levels)
celltype_sample_norm3$Var2 = factor(celltype_sample_norm3$Var2, levels = c('UT','uPAR'))

ggplot(celltype_sample_norm3, aes(x = Var1, y = Freq, fill = Var2)) + geom_bar(stat="identity", position="dodge") + 
  geom_errorbar(aes(ymin=Freq-se, ymax=Freq+se),position=position_dodge(width = 0.85),width=0.3, size=0.25) + theme_vyom + 
  scale_fill_manual(values = colorsType) + expand_limits(y = c(0)) + 
  theme( axis.text.x = element_text(angle = 45, size =  6, vjust = 1, hjust = 1), panel.background = element_rect(fill = "white", colour = "Black"), axis.text.y = element_text(size = 6), axis.title.x = element_text(size = 7), axis.title.y = element_text(size = 7)) + 
  xlab("Cell Type") + ylab("Fraction of Cells") + labs(fill = "Age") + scale_y_continuous(expand = expansion(mult = c(0, .1)))
ggsave( "Fig 3 supp.pdf", width=3.75, height=3, units="in")

heatmap_genes <- unique(c(genes_stem_prolif, genes_antimicrobial, genes_hormone_secretion, genes_nutrient_absorption))


Idents(Cart_old1) <- Cart_old1$Cell_Type
Cell_Types <- levels(Cart_old1$Cell_Type)#c('TReg','Proliferating TReg','NK Cell','M1 Macrophage','M2 Macrophage','Monocyte','DC','Neutrophil')##Cart_old1$Cell_Type
DE_tControlle_all = Idents(Cart_old1) %>% levels() %>% intersect(Cell_Types) %>% lapply(get_lfc_celltype, seurat_obj = Cart_old1, condition_colname = "Treatment", condition_oi = 'uPAR', condition_reference = 'UT', celltype_col = NULL, expression_pct = 0.01) %>% purrr::reduce(full_join)
DE_tControlle_all[is.na(DE_tControlle_all)] = 0

DE_tControlle_all_frame <- data.frame(DE_tControlle_all)
rownames(DE_tControlle_all_frame) <- DE_tControlle_all_frame$gene

gene.list <- heatmap_genes
DE_tControlle_all_filtered <- DE_tControlle_all_frame[gene.list,]
DE_tControlle_all_filtered <- DE_tControlle_all_filtered[complete.cases(DE_tControlle_all_filtered), ]

# make LFC heatmap
lfc_matrix = DE_tControlle_all_filtered  %>% dplyr::select(-gene) %>% as.matrix() %>% magrittr::set_rownames(DE_tControlle_all_filtered$gene)
rownames(lfc_matrix) = rownames(lfc_matrix) %>% make.names()

vis_lfc = lfc_matrix

colnames(vis_lfc) = vis_lfc %>% colnames() %>% make.names()
vis_lfc <- data.frame(vis_lfc)
vis_lfc$gene <- rownames(vis_lfc)

vis_lfc1 <- melt(vis_lfc)
pathwayColorsDiff = rev(brewer.pal(9, "RdBu"))
list.gene <- unique(vis_lfc1$gene)
vis_lfc1$gene <- factor(vis_lfc1$gene,levels = rev(list.gene))
ggplot(vis_lfc1, aes(gene, variable, fill= value)) + geom_tile() + 
  scale_fill_gradientn(name = "Log2FC", colors = pathwayColorsDiff, limits = c(-1,1)) + coord_flip() +
  theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45), axis.text.y =   element_text( vjust = 0.1, hjust = 1), panel.background = element_blank(), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())
cd8_de <- vis_lfc1
ggplot(vis_lfc1, aes(gene, variable, fill= value)) + geom_tile() + 
  scale_fill_gradientn(name = "Log2FC", colors = pathwayColorsDiff, limits = c(-2.5,2)) + coord_flip() +
  theme(legend.position = 'none', axis.text.x = element_text(size = 6,hjust = 1, vjust = 1, angle = 45), axis.text.y =   element_text(size = 6, vjust = 0.1, hjust = 1), panel.background = element_blank(), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())

ggsave(file="Figure 3 I.pdf",  width=1.5, height=4.5, units="in")

Idents(cart_young) <- cart_young$Treatment
DE_treat <- FindMarkers(cart_young, ident.1 = "uPAR", ident.2 = "UT", test.use = "MAST", logfc.threshold = .01, min.pct = .01, assay = 'SCT')
EnhancedVolcano(DE_treat, lab = rownames(DE_treat), x = 'avg_log2FC', y = 'p_val_adj', title = 'uPAR vs UT', pCutoff = 10e-1, FCcutoff = 0.25, xlim = c(-1,1), ylim = c(0,320), subtitle = 'All Cells' )
write.csv(DE_treat,'DE_Cart_treatment_young.csv')


de <- DE_treat
lfc <- .25
fdr <- .05
de$diffexpressed <- "NO"
de$diffexpressed[de$avg_log2FC > lfc & de$p_val_adj < fdr ] <- "UP"
de$diffexpressed[de$avg_log2FC < -lfc & de$p_val_adj < fdr ] <- "DOWN"
de$delabel <- NA
#specify genes
#de$delabel[de$gene_name %in% Interesting_gene_names] <- subset(de$gene_name, de$gene_name %in% Interesting_gene_names) 
#top genes
de$delabel[de$diffexpressed != "NO"] <- rownames(de[de$diffexpressed != "NO",])
secMin <- as.numeric(min( de$p_val_adj[de$p_val_adj!=min(de$p_val_adj)] ))
de["p_val_adj"][de["p_val_adj"] == 0] <- secMin
#gene_name_list <- Interesting_gene_names
gene_name_list <- de$delabel[de$diffexpressed != "NO"]

p1 <- ggplot(data=de, aes(x=avg_log2FC, y=-log10(p_val_adj), label=delabel)) +
  geom_point(aes(col=diffexpressed),size = 1, alpha = .5) + #ylim(c(0,600)) +
  theme(axis.text = element_text(size = 10, color="black"), axis.title = element_text(size = 10), legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_rect(colour = "#2E2E2E", fill=NA, size=0.5)) +
  geom_text_repel(size=3.5, fontface="bold", box.padding = 0.05, segment.size = 0.2, min.segment.length = 0.02, segment.color = '#585858', segment.curvature = 0, segment.ncp = 1, segment.angle = 20, arrow=arrow(angle = 20, length = unit(0.01, "inches"), ends = "last", type = "open")) +
  scale_color_manual(values=c("DOWN" = "#0101DF", "NO" = "#848484","UP" = "#FF0040")) +
  geom_vline(xintercept=c(-lfc, lfc), col="#585858", linetype="longdash", size=0.3) +
  geom_hline(yintercept=-log10(fdr), col="#585858", linetype="longdash", size=0.3) + ggtitle(paste0('uPAR vs UT')) + xlab("log2(fold-change)") + ylab("-log10(Pvalue)")

print(p1)

ggsave(file = 'young_UPARvsUT_Volcano.pdf',  width=5, height=4, units="in")


