
Seurat.object.NEW <- readRDS("motifs_ATAC.rds")
DefaultAssay(Seurat.object.NEW) <- "peaks"

# gather the footprinting information for sets of motifs

Seurat.object <- Footprint(
  object = Seurat.object.NEW,
  motif.name = c("TBR1"),
  genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
  in.peaks = TRUE
)

Seurat.object <- Footprint(
  object = Seurat.object.NEW,
  motif.name = c("HSF4", "GLIS1", "MYC", "NFKB2", "KLF5", "FOXG1"),
  genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
  in.peaks = TRUE
)
# plot the footprint data for each group of cells

Seurat.object$tags <- Idents(Seurat.object)

saveRDS(Seurat.object, "footprint_TFs.rds")

## TBR1

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
  reduction = "umap", dims = c(1,2),cols = c("gray", "purple"),
  raster = FALSE
) +
  theme(text = element_text(size = 20))

jpeg("images/TBR1_activity.jpeg", units="in", width=10, height=10, res=300)
p4
dev.off()



## HSF4

p2 <- PlotFootprint(Seurat.object, features = "HSF4",
                    # idents = c("Pv", "RORB+","Ex", "Vip", "Non-Vip", "Sst"),
                    group.by = c("disease")) +
  theme(text = element_text(size = 20))


p3 <- PlotFootprint(Seurat.object, features = "HSF4",
                    idents = c("Pv", "RORB+","Ex", "Vip", "Non-Vip", "Sst")) +
  theme(text = element_text(size = 20))


jpeg("images/HSF4_disease_motif.jpeg", units="in", width=10, height=10, res=300)
p2
dev.off()


jpeg("images/HSF4_motif.jpeg", units="in", width=10, height=10, res=300)
p3
dev.off()

Seurat.object_pp <- Seurat.object
DefaultAssay(Seurat.object_pp) <- 'chromvar'

# look at the activity of Mef2c

colnames(Seurat.object_pp@reductions$umap@cell.embeddings) <- c("UMAP_1", "UMAP_2")
p4 <- Seurat::FeaturePlot(
  object = Seurat.object_pp,
  features = "MA0771.1",
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1,
  reduction = "umap", dims = c(1,2),cols = c("gray", "purple"),
  raster = FALSE
) +
  theme(text = element_text(size = 20))

jpeg("images/HSF4_activity.jpeg", units="in", width=10, height=10, res=300)
p4
dev.off()



## GLIS1

p2 <- PlotFootprint(Seurat.object, features = "GLIS1",
                    # idents = c("Pv", "RORB+","Ex", "Vip", "Non-Vip", "Sst"),
                    group.by = c("disease")) +
  theme(text = element_text(size = 20))


p3 <- PlotFootprint(Seurat.object, features = "GLIS1",
                    idents = c("Pv", "RORB+","Ex", "Vip", "Non-Vip", "Sst")) +
  theme(text = element_text(size = 20))


jpeg("images/GLIS1_disease_motif.jpeg", units="in", width=10, height=10, res=300)
p2
dev.off()


jpeg("images/GLIS1_motif.jpeg", units="in", width=10, height=10, res=300)
p3
dev.off()

Seurat.object_pp <- Seurat.object
DefaultAssay(Seurat.object_pp) <- 'chromvar'

# look at the activity of Mef2c

colnames(Seurat.object_pp@reductions$umap@cell.embeddings) <- c("UMAP_1", "UMAP_2")
p4 <- Seurat::FeaturePlot(
  object = Seurat.object_pp,
  features = "MA0735.1",
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1,
  reduction = "umap", dims = c(1,2),cols = c("gray", "purple"),
  raster = FALSE
) +
  theme(text = element_text(size = 20))

jpeg("images/GLIS1_activity.jpeg", units="in", width=10, height=10, res=300)
p4
dev.off()


## MYC

p2 <- PlotFootprint(Seurat.object, features = "MYC",
                    # idents = c("Pv", "RORB+","Ex", "Vip", "Non-Vip", "Sst"),
                    group.by = c("disease")) +
  theme(text = element_text(size = 20))


p3 <- PlotFootprint(Seurat.object, features = "MYC",
                    idents = c("Pv", "RORB+","Ex", "Vip", "Non-Vip", "Sst")) +
  theme(text = element_text(size = 20))


jpeg("images/MYC_disease_motif.jpeg", units="in", width=10, height=10, res=300)
p2
dev.off()


jpeg("images/MYC_motif.jpeg", units="in", width=10, height=10, res=300)
p3
dev.off()

Seurat.object_pp <- Seurat.object
DefaultAssay(Seurat.object_pp) <- 'chromvar'

# look at the activity of Mef2c

colnames(Seurat.object_pp@reductions$umap@cell.embeddings) <- c("UMAP_1", "UMAP_2")
p4 <- Seurat::FeaturePlot(
  object = Seurat.object_pp,
  features = "MA0147.3",
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1,
  reduction = "umap", dims = c(1,2),cols = c("gray", "purple"),
  raster = FALSE
) +
  theme(text = element_text(size = 20))

jpeg("images/MYC_activity.jpeg", units="in", width=10, height=10, res=300)
p4
dev.off()


## NFKB2

p2 <- PlotFootprint(Seurat.object, features = "NFKB2",
                    # idents = c("Pv", "RORB+","Ex", "Vip", "Non-Vip", "Sst"),
                    group.by = c("disease")) +
  theme(text = element_text(size = 20))


p3 <- PlotFootprint(Seurat.object, features = "NFKB2",
                    idents = c("Pv", "RORB+","Ex", "Vip", "Non-Vip", "Sst")) +
  theme(text = element_text(size = 20))


jpeg("images/NFKB2_disease_motif.jpeg", units="in", width=10, height=10, res=300)
p2
dev.off()


jpeg("images/NFKB2_motif.jpeg", units="in", width=10, height=10, res=300)
p3
dev.off()

Seurat.object_pp <- Seurat.object
DefaultAssay(Seurat.object_pp) <- 'chromvar'

# look at the activity of Mef2c

colnames(Seurat.object_pp@reductions$umap@cell.embeddings) <- c("UMAP_1", "UMAP_2")
p4 <- Seurat::FeaturePlot(
  object = Seurat.object_pp,
  features = "MA0778.1",
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1,
  reduction = "umap", dims = c(1,2),cols = c("gray", "purple"),
  raster = FALSE
) +
  theme(text = element_text(size = 20))

jpeg("images/NFKB2_activity.jpeg", units="in", width=10, height=10, res=300)
p4
dev.off()


## KLF5

p2 <- PlotFootprint(Seurat.object, features = "KLF5",
                    # idents = c("Pv", "RORB+","Ex", "Vip", "Non-Vip", "Sst"),
                    group.by = c("disease")) +
  theme(text = element_text(size = 20))


p3 <- PlotFootprint(Seurat.object, features = "KLF5",
                    idents = c("Pv", "RORB+","Ex", "Vip", "Non-Vip", "Sst")) +
  theme(text = element_text(size = 20))


jpeg("images/KLF5_disease_motif.jpeg", units="in", width=10, height=10, res=300)
p2
dev.off()


jpeg("images/KLF5_motif.jpeg", units="in", width=10, height=10, res=300)
p3
dev.off()

Seurat.object_pp <- Seurat.object
DefaultAssay(Seurat.object_pp) <- 'chromvar'

# look at the activity of Mef2c

colnames(Seurat.object_pp@reductions$umap@cell.embeddings) <- c("UMAP_1", "UMAP_2")
p4 <- Seurat::FeaturePlot(
  object = Seurat.object_pp,
  features = "MA0599.1",
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1,
  reduction = "umap", dims = c(1,2),cols = c("gray", "purple"),
  raster = FALSE
) +
  theme(text = element_text(size = 20))

jpeg("images/KLF5_activity.jpeg", units="in", width=10, height=10, res=300)
p4
dev.off()



## FOXG1

p2 <- PlotFootprint(Seurat.object, features = "FOXG1",
                    # idents = c("Pv", "RORB+","Ex", "Vip", "Non-Vip", "Sst"),
                    group.by = c("disease")) +
  theme(text = element_text(size = 20))


p3 <- PlotFootprint(Seurat.object, features = "FOXG1",
                    idents = c("Pv", "RORB+","Ex", "Vip", "Non-Vip", "Sst")) +
  theme(text = element_text(size = 20))


jpeg("images/FOXG1_disease_motif.jpeg", units="in", width=10, height=10, res=300)
p2
dev.off()


jpeg("images/FOXG1_motif.jpeg", units="in", width=10, height=10, res=300)
p3
dev.off()

Seurat.object_pp <- Seurat.object
DefaultAssay(Seurat.object_pp) <- 'chromvar'

# look at the activity of Mef2c

colnames(Seurat.object_pp@reductions$umap@cell.embeddings) <- c("UMAP_1", "UMAP_2")
p4 <- Seurat::FeaturePlot(
  object = Seurat.object_pp,
  features = "MA0613.1",
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1,
  reduction = "umap", dims = c(1,2),cols = c("gray", "purple"),
  raster = FALSE
) +
  theme(text = element_text(size = 20))

jpeg("images/FOXG1_activity.jpeg", units="in", width=10, height=10, res=300)
p4
dev.off()


## Other motifs

Seurat.object.NEW <- readRDS("footprint_TFs.rds")


Seurat.object <- Footprint(
  object = Seurat.object.NEW,
  motif.name = c("ZNF263", "FOXP1", "STAT1", "KLF5", "ELK1", "POU2F1"),
  genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
  in.peaks = TRUE
)
# plot the footprint data for each group of cells

Seurat.object$tags <- Idents(Seurat.object)

saveRDS(Seurat.object, "footprint_TFs.rds")



## ZNF263

p2 <- PlotFootprint(Seurat.object, features = "ZNF263",
                    # idents = c("Pv", "RORB+","Ex", "Vip", "Non-Vip", "Sst"),
                    group.by = c("disease")) +
  theme(text = element_text(size = 20))


p3 <- PlotFootprint(Seurat.object, features = "ZNF263",
                    idents = c("Pv", "RORB+","Ex", "Vip", "Non-Vip", "Sst")) +
  theme(text = element_text(size = 20))


jpeg("images/ZNF263_disease_motif.jpeg", units="in", width=10, height=10, res=300)
p2
dev.off()


jpeg("images/ZNF263_motif.jpeg", units="in", width=10, height=10, res=300)
p3
dev.off()

Seurat.object_pp <- Seurat.object
DefaultAssay(Seurat.object_pp) <- 'chromvar'

colnames(Seurat.object_pp@reductions$umap@cell.embeddings) <- c("UMAP_1", "UMAP_2")
p4 <- Seurat::FeaturePlot(
  object = Seurat.object_pp,
  features = "MA0528.2",
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1,
  reduction = "umap", dims = c(1,2),cols = c("gray", "purple"),
  raster = FALSE
) +
  theme(text = element_text(size = 20))

jpeg("images/ZNF263_activity.jpeg", units="in", width=10, height=10, res=300)
p4
dev.off()


## FOXP1

p2 <- PlotFootprint(Seurat.object, features = "FOXP1",
                    # idents = c("Pv", "RORB+","Ex", "Vip", "Non-Vip", "Sst"),
                    group.by = c("disease")) +
  theme(text = element_text(size = 20))


p3 <- PlotFootprint(Seurat.object, features = "FOXP1",
                    idents = c("Pv", "RORB+","Ex", "Vip", "Non-Vip", "Sst")) +
  theme(text = element_text(size = 20))


jpeg("images/FOXP1_disease_motif.jpeg", units="in", width=10, height=10, res=300)
p2
dev.off()


jpeg("images/FOXP1_motif.jpeg", units="in", width=10, height=10, res=300)
p3
dev.off()

Seurat.object_pp <- Seurat.object
DefaultAssay(Seurat.object_pp) <- 'chromvar'

colnames(Seurat.object_pp@reductions$umap@cell.embeddings) <- c("UMAP_1", "UMAP_2")
p4 <- Seurat::FeaturePlot(
  object = Seurat.object_pp,
  features = "MA0481.3",
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1,
  reduction = "umap", dims = c(1,2),cols = c("gray", "purple"),
  raster = FALSE
) +
  theme(text = element_text(size = 20))

jpeg("images/FOXP1_activity.jpeg", units="in", width=10, height=10, res=300)
p4
dev.off()



## STAT1

p2 <- PlotFootprint(Seurat.object, features = "STAT1",
                    # idents = c("Pv", "RORB+","Ex", "Vip", "Non-Vip", "Sst"),
                    group.by = c("disease")) +
  theme(text = element_text(size = 20))


p3 <- PlotFootprint(Seurat.object, features = "STAT1",
                    idents = c("Pv", "RORB+","Ex", "Vip", "Non-Vip", "Sst")) +
  theme(text = element_text(size = 20))


jpeg("images/STAT1_disease_motif.jpeg", units="in", width=10, height=10, res=300)
p2
dev.off()


jpeg("images/STAT1_motif.jpeg", units="in", width=10, height=10, res=300)
p3
dev.off()

Seurat.object_pp <- Seurat.object
DefaultAssay(Seurat.object_pp) <- 'chromvar'

colnames(Seurat.object_pp@reductions$umap@cell.embeddings) <- c("UMAP_1", "UMAP_2")
p4 <- Seurat::FeaturePlot(
  object = Seurat.object_pp,
  features = "MA0137.3",
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1,
  reduction = "umap", dims = c(1,2),cols = c("gray", "purple"),
  raster = FALSE
) +
  theme(text = element_text(size = 20))

jpeg("images/STAT1_activity.jpeg", units="in", width=10, height=10, res=300)
p4
dev.off()



## KLF5

p2 <- PlotFootprint(Seurat.object, features = "KLF5",
                    # idents = c("Pv", "RORB+","Ex", "Vip", "Non-Vip", "Sst"),
                    group.by = c("disease")) +
  theme(text = element_text(size = 20))


p3 <- PlotFootprint(Seurat.object, features = "KLF5",
                    idents = c("Pv", "RORB+","Ex", "Vip", "Non-Vip", "Sst")) +
  theme(text = element_text(size = 20))


jpeg("images/KLF5_disease_motif.jpeg", units="in", width=10, height=10, res=300)
p2
dev.off()


jpeg("images/KLF5_motif.jpeg", units="in", width=10, height=10, res=300)
p3
dev.off()

Seurat.object_pp <- Seurat.object
DefaultAssay(Seurat.object_pp) <- 'chromvar'

colnames(Seurat.object_pp@reductions$umap@cell.embeddings) <- c("UMAP_1", "UMAP_2")
p4 <- Seurat::FeaturePlot(
  object = Seurat.object_pp,
  features = "MA0599.1",
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1,
  reduction = "umap", dims = c(1,2),cols = c("gray", "purple"),
  raster = FALSE
) +
  theme(text = element_text(size = 20))

jpeg("images/KLF5_activity.jpeg", units="in", width=10, height=10, res=300)
p4
dev.off()




## ELK1

p2 <- PlotFootprint(Seurat.object, features = "ELK1",
                    # idents = c("Pv", "RORB+","Ex", "Vip", "Non-Vip", "Sst"),
                    group.by = c("disease")) +
  theme(text = element_text(size = 20))


p3 <- PlotFootprint(Seurat.object, features = "ELK1",
                    idents = c("Pv", "RORB+","Ex", "Vip", "Non-Vip", "Sst")) +
  theme(text = element_text(size = 20))


jpeg("images/ELK1_disease_motif.jpeg", units="in", width=10, height=10, res=300)
p2
dev.off()


jpeg("images/ELK1_motif.jpeg", units="in", width=10, height=10, res=300)
p3
dev.off()

Seurat.object_pp <- Seurat.object
DefaultAssay(Seurat.object_pp) <- 'chromvar'

colnames(Seurat.object_pp@reductions$umap@cell.embeddings) <- c("UMAP_1", "UMAP_2")
p4 <- Seurat::FeaturePlot(
  object = Seurat.object_pp,
  features = "MA0028.2",
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1,
  reduction = "umap", dims = c(1,2),cols = c("gray", "purple"),
  raster = FALSE
) +
  theme(text = element_text(size = 20))

jpeg("images/ELK1_activity.jpeg", units="in", width=10, height=10, res=300)
p4
dev.off()



## POU2F1

p2 <- PlotFootprint(Seurat.object, features = "POU2F1",
                    # idents = c("Pv", "RORB+","Ex", "Vip", "Non-Vip", "Sst"),
                    group.by = c("disease")) +
  theme(text = element_text(size = 20))


p3 <- PlotFootprint(Seurat.object, features = "POU2F1",
                    idents = c("Pv", "RORB+","Ex", "Vip", "Non-Vip", "Sst")) +
  theme(text = element_text(size = 20))


jpeg("images/POU2F1_disease_motif.jpeg", units="in", width=10, height=10, res=300)
p2
dev.off()


jpeg("images/POU2F1_motif.jpeg", units="in", width=10, height=10, res=300)
p3
dev.off()

Seurat.object_pp <- Seurat.object
DefaultAssay(Seurat.object_pp) <- 'chromvar'

colnames(Seurat.object_pp@reductions$umap@cell.embeddings) <- c("UMAP_1", "UMAP_2")
p4 <- Seurat::FeaturePlot(
  object = Seurat.object_pp,
  features = "MA0785.1",
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1,
  reduction = "umap", dims = c(1,2),cols = c("gray", "purple"),
  raster = FALSE
) +
  theme(text = element_text(size = 20))

jpeg("images/POU2F1_activity.jpeg", units="in", width=10, height=10, res=300)
p4
dev.off()

###################

 Seurat::FeaturePlot(
  object = Seurat.object_pp,
  features = "MA1650.1",
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1,
  reduction = "umap", dims = c(1,2),cols = c("gray", "purple"),
  raster = FALSE
) +
  theme(text = element_text(size = 20))



## Other motifs

Seurat.object.NEW <- readRDS("footprint_TFs.rds")

DefaultAssay(Seurat.object) <- "peaks"

Seurat.object <- Footprint(
  object = Seurat.object,
  motif.name = c("KLF16", "HINFP", "PATZ1", "TFDP1", "ZNF148", "FOS", "FOSL2", "MEF2C", "DBP"),
  genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
  in.peaks = TRUE
)
# plot the footprint data for each group of cells

Seurat.object$tags <- Idents(Seurat.object)

saveRDS(Seurat.object, "footprint_TFs.rds")

