# Create a function to convert a Seurat object to an AnnData object
so <- readRDS("anotation_SFG.rds")
so_morabito <- so
so_morabito$tags <- Idents(so_morabito)
cells_rorb <- na.omit(colnames(so_morabito)[so_morabito$tags=="RORB+"])
so_morabito_RORB <- so_morabito[ ,cells_rorb]

atac_rorb <- readr::read_csv("atac_cellType_markers/RORB+_atac.csv") %>%
  atac_input()


library(Seurat)
library(SeuratData)
library(SeuratDisk)
SaveH5Seurat(seurat_obj, filename = "rorb.h5Seurat")
Convert("rorb.h5Seurat", dest = "h5ad")


