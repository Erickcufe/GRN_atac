
Seurat.object.NEW <- readRDS("motifs_ATAC.rds")
DefaultAssay(Seurat.object.NEW) <- "peaks"

# gather the footprinting information for sets of motifs

Seurat.object <- Footprint(
  object = Seurat.object.NEW,
  motif.name = c("TBR1"),
  genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
  in.peaks = TRUE
)

# plot the footprint data for each group of cells
p2 <- PlotFootprint(Seurat.object, features = "TBR1")
Seurat.object$tags <- Idents(Seurat.object)

cells <- colnames(Seurat.object)[Seurat.object$tags %in% c("Pv", "RORB+","Ex", "Vip", "Non-Vip", "Sst")]
Seurat.object_pp <- Seurat.object[, cells]
p2 <- PlotFootprint(Seurat.object, features = "TBR1",
                    idents = c("Pv", "RORB+","Ex", "Vip", "Non-Vip", "Sst"),
                    split.by = c("disease")) +
  theme(text = element_text(size = 20))


jpeg("images/TBR1_motif.jpeg", units="in", width=10, height=10, res=300)
p2
dev.off()


DefaultAssay(mouse_brain) <- 'chromvar'

# look at the activity of Mef2c
p2 <- FeaturePlot(
  object = Seurat.object,
  features = "MA0497.1",
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1
)
p1 + p2
