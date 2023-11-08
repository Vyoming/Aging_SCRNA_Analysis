#Vyom Shah

#Import Dataset
Cvelo <- readRDS("./data/cancer_velocity/CSpat_132Tx_CancerOnly.Rds", refhook = NULL)

#cell reclustering 
Cvelo <- NormalizeData(Cvelo)
Cvelo <- FindVariableFeatures(Cvelo, selection.method = "vst", nfeatures = 2000, verbose = TRUE)
Cvelo <- ScaleData(Cvelo, verbose = TRUE, features = rownames(Cvelo) )
pcs <- 15
mres <- 0.3
DefaultAssay(Cvelo) <- "RNA"
Cvelo <- RunPCA(Cvelo, verbose = TRUE, npcs = pcs)
Cvelo <- FindNeighbors(Cvelo, reduction = "pca", dims = 1:pcs, verbose = TRUE)
Cvelo <- FindClusters(Cvelo, resolution = mres, algorithm=1, verbose = TRUE)
Cvelo <- RunUMAP(Cvelo, reduction = "pca", dims = 1:pcs, verbose = TRUE)

Cvelo$RNA_snn_res.0.3
DimPlot(Cvelo, group.by = 'RNA_snn_res.0.3')
DimPlot(Cvelo, group.by = 'PatientIDx')
Cvelo$Cell_types

#save important files
cvelo_meta <- Cvelo@meta.data
write.csv(cvelo_meta, './data/cancer_velocity/cvelo_meta.csv')

Cvelo[["umap"]]@cell.embeddings
write.csv(Cvelo[["umap"]]@cell.embeddings, './data/cancer_velocity/cvelo_umap.csv')

Cvelo.loom <- as.loom(Cvelo, filename = "./data/cancer_velocity/Cvelo_sobj.loom", assay = 'RNA', verbose = TRUE)
Cvelo.loom
Cvelo.loom$close_all()

library(SeuratDisk)
SaveH5Seurat(Cvelo, filename = "./data/cancer_velocity/Cvelo_sobj.h5Seurat")
Convert("./data/cancer_velocity/Cvelo_sobj.h5Seurat", dest = "h5ad", overwrite = TRUE)
