# Create a function to convert a Seurat object to an AnnData object
so <- readRDS("anotation_SFG.rds")

so_morabito$tags <- Idents(so_morabito)
cells_rorb <- na.omit(colnames(so_morabito)[so_morabito$tags=="RORB+"])
so_morabito_RORB <- so_morabito[ ,cells_rorb]

atac_rorb <- readr::read_csv("atac_cellType_markers/RORB+_atac.csv") %>%
  atac_input()


seurat_to_anndata <- function(seurat_obj) {
  # Extract the data from the Seurat object
  counts <- as.matrix(seurat_obj@assays$RNA@counts)
  obs <- seurat_obj@meta.data
  var <- rowData(seurat_obj)

  # Create an AnnData object
  adata <- anndata::AnnData(X = counts, obs = obs, var = var)

  return(adata)
}
