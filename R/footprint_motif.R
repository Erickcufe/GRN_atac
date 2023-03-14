
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


p2 <- PlotFootprint(Seurat.object, features = "TBR1",
                    # idents = c("Pv", "RORB+","Ex", "Vip", "Non-Vip", "Sst"),
                    group.by = c("disease")) +
  theme(text = element_text(size = 20))


p3 <- PlotFootprint(Seurat.object, features = "TBR1",
                    idents = c("Pv", "RORB+","Ex", "Vip", "Non-Vip", "Sst")) +
  theme(text = element_text(size = 20))


jpeg("images/TBR1_disease_motif.jpeg", units="in", width=10, height=10, res=300)
p2
dev.off()


jpeg("images/TBR1_motif.jpeg", units="in", width=10, height=10, res=300)
p3
dev.off()

Seurat.object_pp <- Seurat.object
DefaultAssay(Seurat.object_pp) <- 'chromvar'

# look at the activity of Mef2c

colnames(Seurat.object_pp@reductions$umap@cell.embeddings) <- c("UMAP_1", "UMAP_2")

p4 <- Seurat::FeaturePlot(
  object = Seurat.object_pp,
  features = "MA0802.1",
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1,
  reduction = "umap", dims = c(1,2),cols = c("gray", "red"),
  raster = FALSE
) +
  theme(text = element_text(size = 20))
p4
