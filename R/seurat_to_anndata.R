# Create a function to convert a Seurat object to an AnnData object
so <- readRDS("anotation_SFG.rds")

seurat_to_anndata <- function(seurat_obj) {
  # Extract the data from the Seurat object
  counts <- as.matrix(seurat_obj@assays$RNA@counts)
  obs <- seurat_obj@meta.data
  var <- rowData(seurat_obj)

  # Create an AnnData object
  adata <- anndata::AnnData(X = counts, obs = obs, var = var)

  return(adata)
}
