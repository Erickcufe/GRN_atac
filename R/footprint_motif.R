
Seurat.object.NEW <- readRDS("motifs_ATAC.rds")
DefaultAssay(Seurat.object.NEW) <- "peaks"

Seurat.object <- CreateChromatinAssay(counts = Seurat.object.NEW@assays$peaks@counts,
                                                     genome = "hg38")
# extract position frequency matrices for the motifs
pwm <- getMatrixSet(
  x = JASPAR2022::JASPAR2022,
  opts = list(species = 9606, all_versions = FALSE)
)

# add motif information
Seurat.object <- AddMotifs(Seurat.object,
                           genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
                           pfm = pwm)

Seurat.object@assays$peaks@fragments <- nuc@assays$peak@fragments
# gather the footprinting information for sets of motifs
Seurat.object <- Footprint(
  object = Seurat.object,
  motif.name = c("GATA2","CEBPA","EBF1"),
  genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
  in.peaks = TRUE
)

# plot the footprint data for each group of cells
p2 <- PlotFootprint(Seurat.object.NEW, features = c("GATA2", "CEBPA", "EBF1"))
