library(ggplot2)
library(TFBSTools)
library(Seurat)
Seurat.object.NEW <- readRDS("motifs_ATAC.rds")

###############
##############
# RORB+ NEURONS
###############

differential.activity <- FindMarkers(
  object = Seurat.object.NEW,
  ident.1 = 'RORB+',
  ident.2 = c('Pv', 'Ex', 'Vip', 'Non-Vip', 'Sst'),
  only.pos = TRUE,
  test.use = 'LR',
  min.pct = 0.05,
  latent.vars = 'nCount_peaks'
)

differential.activity <- differential.activity[differential.activity$p_val_adj < 0.05,]
differential.activity$regions <- rownames(differential.activity)

readr::write_csv(differential.activity, "atac_cellType_markers/RORB+_atac.csv")

# get top differentially accessible peaks
top.da.peak <- rownames(differential.activity[differential.activity$p_val < 0.005, ])
enriched.motifs_RORB_Ex <- FindMotifs(
  object = Seurat.object.NEW,
  features = top.da.peak
)

readr::write_csv(enriched.motifs_RORB_Ex, "neuronal_motifs/RORB+_motifs.csv")

jpeg("images/RORB_atac_markers.jpeg", units="in", width=10, height=10, res=300)
MotifPlot(
  object = Seurat.object.NEW,
  motifs = head(rownames(enriched.motifs_RORB_Ex))
) + theme(text = element_text(size = 18))

dev.off()


###############
##############
# Pv NEURONS
###############


DefaultAssay(Seurat.object.NEW) <- "peaks"
differential.activity <- FindMarkers(
  object = Seurat.object.NEW,
  ident.1 = 'Pv',
  ident.2 = c('RORB+', 'Ex', 'Vip', 'Non-Vip', 'Sst'),
  only.pos = TRUE,
  test.use = 'LR',
  min.pct = 0.05,
  latent.vars = 'nCount_peaks'
)

differential.activity <- differential.activity[differential.activity$p_val_adj < 0.05,]
differential.activity$regions <- rownames(differential.activity)

readr::write_csv(differential.activity, "atac_cellType_markers/Pv_atac.csv")

# get top differentially accessible peaks
top.da.peak <- rownames(differential.activity[differential.activity$p_val < 0.005, ])
enriched.motifs_Pv_Ex <- FindMotifs(
  object = Seurat.object.NEW,
  features = top.da.peak
)

readr::write_csv(enriched.motifs_Pv_Ex, "neuronal_motifs/Pv_motifs.csv")

jpeg("images/Pv_atac_markers.jpeg", units="in", width=10, height=10, res=300)
MotifPlot(
  object = Seurat.object.NEW,
  motifs = head(rownames(enriched.motifs_Pv_Ex))
) + theme(text = element_text(size = 18))

dev.off()


###############
##############
# Ex NEURONS
###############

differential.activity <- FindMarkers(
  object = Seurat.object.NEW,
  ident.1 = 'Ex',
  ident.2 = c('RORB+', 'Pv', 'Vip', 'Non-Vip', 'Sst'),
  only.pos = TRUE,
  test.use = 'LR',
  min.pct = 0.05,
  latent.vars = 'nCount_peaks'
)

differential.activity <- differential.activity[differential.activity$p_val_adj < 0.05,]
differential.activity$regions <- rownames(differential.activity)

readr::write_csv(differential.activity, "atac_cellType_markers/Ex_atac.csv")

# get top differentially accessible peaks
top.da.peak <- rownames(differential.activity[differential.activity$p_val < 0.005, ])
enriched.motifs_Ex_Ex <- FindMotifs(
  object = Seurat.object.NEW,
  features = top.da.peak
)

readr::write_csv(enriched.motifs_Ex_Ex, "neuronal_motifs/Ex_motifs.csv")

jpeg("images/Ex_atac_markers.jpeg", units="in", width=10, height=10, res=300)
MotifPlot(
  object = Seurat.object.NEW,
  motifs = head(rownames(enriched.motifs_Ex_Ex))
) + theme(text = element_text(size = 18))

dev.off()


##############
# Vip NEURONS
###############

differential.activity <- FindMarkers(
  object = Seurat.object.NEW,
  ident.1 = 'Vip',
  ident.2 = c('RORB+', 'Pv', 'Ex', 'Non-Vip', 'Sst'),
  only.pos = TRUE,
  test.use = 'LR',
  min.pct = 0.05,
  latent.vars = 'nCount_peaks'
)

differential.activity <- differential.activity[differential.activity$p_val_adj < 0.05,]
differential.activity$regions <- rownames(differential.activity)

readr::write_csv(differential.activity, "atac_cellType_markers/Vip_atac.csv")

# get top differentially accessible peaks
top.da.peak <- rownames(differential.activity[differential.activity$p_val < 0.005, ])
enriched.motifs_Vip_Vip <- FindMotifs(
  object = Seurat.object.NEW,
  features = top.da.peak
)

readr::write_csv(enriched.motifs_Vip_Vip, "neuronal_motifs/Vip_motifs.csv")

jpeg("images/Vip_atac_markers.jpeg", units="in", width=10, height=10, res=300)
MotifPlot(
  object = Seurat.object.NEW,
  motifs = head(rownames(enriched.motifs_Vip_Vip))
) + theme(text = element_text(size = 18))

dev.off()


##############
# Non-Vip NEURONS
###############
DefaultAssay(Seurat.object.NEW)

differential.activity <- FindMarkers(
  object = Seurat.object.NEW,
  ident.1 = 'Non-Vip',
  ident.2 = c('RORB+', 'Pv', 'Ex', 'Vip', 'Sst'),
  only.pos = TRUE,
  test.use = 'LR',
  min.pct = 0.05,
  latent.vars = 'nCount_peaks'
)

differential.activity <- differential.activity[differential.activity$p_val_adj < 0.05,]
differential.activity$regions <- rownames(differential.activity)

readr::write_csv(differential.activity, "atac_cellType_markers/Non-Vip_atac.csv")

# get top differentially accessible peaks
top.da.peak <- rownames(differential.activity[differential.activity$p_val < 0.005, ])
enriched.motifs_Non_Vip <- FindMotifs(
  object = Seurat.object.NEW,
  features = top.da.peak
)

readr::write_csv(enriched.motifs_Non_Vip, "neuronal_motifs/Non-Vip_motifs.csv")

jpeg("images/Non-Vip_atac_markers.jpeg", units="in", width=10, height=10, res=300)
MotifPlot(
  object = Seurat.object.NEW,
  motifs = head(rownames(enriched.motifs_Non_Vip))
) + theme(text = element_text(size = 18))

dev.off()


##############
# Sst NEURONS
###############

differential.activity <- FindMarkers(
  object = Seurat.object.NEW,
  ident.1 = 'Sst',
  ident.2 = c('RORB+', 'Pv', 'Ex', 'Vip', 'Non-Vip'),
  only.pos = TRUE,
  test.use = 'LR',
  min.pct = 0.05,
  latent.vars = 'nCount_peaks'
)

differential.activity <- differential.activity[differential.activity$p_val_adj < 0.05,]
differential.activity$regions <- rownames(differential.activity)

readr::write_csv(differential.activity, "atac_cellType_markers/Sst_atac.csv")

# get top differentially accessible peaks
top.da.peak <- rownames(differential.activity[differential.activity$p_val < 0.005, ])
enriched.motifs_Sst <- FindMotifs(
  object = Seurat.object.NEW,
  features = top.da.peak
)

readr::write_csv(enriched.motifs_Sst, "neuronal_motifs/Sst_motifs.csv")

jpeg("images/Sst_atac_markers.jpeg", units="in", width=10, height=10, res=300)
MotifPlot(
  object = Seurat.object.NEW,
  motifs = head(rownames(enriched.motifs_Sst))
) + theme(text = element_text(size = 18))

dev.off()



# p1 <- DimPlot(Seurat.object.NEW, label = TRUE, pt.size = 0.1) + NoLegend()



############################
############################
########## CHROMVAR
DefaultAssay(Seurat.object.NEW) <- "chromvar"

## Pv

differential.activity <- FindMarkers(
  object = Seurat.object.NEW,
  ident.1 = 'Pv',
  ident.2 = c('RORB+', 'Ex', 'Vip', 'Non-Vip', 'Sst'),
  only.pos = TRUE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff"
)

differential.activity <- differential.activity[differential.activity$p_val_adj < 0.05,]
differential.activity$motif <- rownames(differential.activity)
readr::write_csv(differential.activity, "chromvar_motifs/pv_chromvar.csv")

jpeg("images/chromvar_Pv_atac_markers.jpeg", units="in", width=10, height=10, res=300)
MotifPlot(
  object = Seurat.object.NEW,
  motifs = head(rownames(differential.activity)),
  assay = 'peaks'
) + theme(text = element_text(size = 18))
dev.off()


## RORB+

differential.activity <- FindMarkers(
  object = Seurat.object.NEW,
  ident.1 = 'RORB+',
  ident.2 = c('Pv', 'Ex', 'Vip', 'Non-Vip', 'Sst'),
  only.pos = TRUE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff"
)

differential.activity <- differential.activity[differential.activity$p_val_adj < 0.05,]
differential.activity$motif <- rownames(differential.activity)
readr::write_csv(differential.activity, "chromvar_motifs/RORB_chromvar.csv")

jpeg("images/chromvar_RORB_atac_markers.jpeg", units="in", width=10, height=10, res=300)
MotifPlot(
  object = Seurat.object.NEW,
  motifs = head(rownames(differential.activity)),
  assay = 'peaks'
) + theme(text = element_text(size = 18))
dev.off()


## Vip

differential.activity <- FindMarkers(
  object = Seurat.object.NEW,
  ident.1 = 'Vip',
  ident.2 = c('Pv', 'Ex', 'RORB+', 'Non-Vip', 'Sst'),
  only.pos = TRUE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff"
)

differential.activity <- differential.activity[differential.activity$p_val_adj < 0.05,]
differential.activity$motif <- rownames(differential.activity)
readr::write_csv(differential.activity, "chromvar_motifs/Vip_chromvar.csv")

jpeg("images/chromvar_Vip_atac_markers.jpeg", units="in", width=10, height=10, res=300)
MotifPlot(
  object = Seurat.object.NEW,
  motifs = head(rownames(differential.activity)),
  assay = 'peaks'
) + theme(text = element_text(size = 18))
dev.off()


# ## Vip
#
# differential.activity <- FindMarkers(
#   object = Seurat.object.NEW,
#   ident.1 = 'Vip',
#   ident.2 = c('Pv', 'Ex', 'RORB+', 'Non-Vip', 'Sst'),
#   only.pos = TRUE,
#   mean.fxn = rowMeans,
#   fc.name = "avg_diff"
# )
#
# differential.activity <- differential.activity[differential.activity$p_val_adj < 0.05,]
# differential.activity$motif <- rownames(differential.activity)
# readr::write_csv(differential.activity, "chromvar_motifs/Vip_chromvar.csv")
#
# jpeg("images/chromvar_Vip_atac_markers.jpeg", units="in", width=10, height=10, res=300)
# MotifPlot(
#   object = Seurat.object.NEW,
#   motifs = head(rownames(differential.activity)),
#   assay = 'peaks'
# ) + theme(text = element_text(size = 18))
# dev.off()



## Non-Vip

differential.activity <- FindMarkers(
  object = Seurat.object.NEW,
  ident.1 = 'Non-Vip',
  ident.2 = c('Pv', 'Ex', 'RORB+', 'Vip', 'Sst'),
  only.pos = TRUE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff"
)

differential.activity <- differential.activity[differential.activity$p_val_adj < 0.05,]
differential.activity$motif <- rownames(differential.activity)
readr::write_csv(differential.activity, "chromvar_motifs/Non_Vip_chromvar.csv")

jpeg("images/chromvar_Non_Vip_atac_markers.jpeg", units="in", width=10, height=10, res=300)
MotifPlot(
  object = Seurat.object.NEW,
  motifs = head(rownames(differential.activity)),
  assay = 'peaks'
) + theme(text = element_text(size = 18))
dev.off()


## Sst

differential.activity <- FindMarkers(
  object = Seurat.object.NEW,
  ident.1 = 'Sst',
  ident.2 = c('Pv', 'Ex', 'RORB+', 'Vip', 'Non-Vip'),
  only.pos = TRUE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff"
)

differential.activity <- differential.activity[differential.activity$p_val_adj < 0.05,]
differential.activity$motif <- rownames(differential.activity)
readr::write_csv(differential.activity, "chromvar_motifs/Sst_chromvar.csv")

jpeg("images/chromvar_Sst_atac_markers.jpeg", units="in", width=10, height=10, res=300)
MotifPlot(
  object = Seurat.object.NEW,
  motifs = head(rownames(differential.activity)),
  assay = 'peaks'
) + theme(text = element_text(size = 18))
dev.off()


## Ex

differential.activity <- FindMarkers(
  object = Seurat.object.NEW,
  ident.1 = 'Ex',
  ident.2 = c('Pv', 'Sst', 'RORB+', 'Vip', 'Non-Vip'),
  only.pos = TRUE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff"
)

differential.activity <- differential.activity[differential.activity$p_val_adj < 0.05,]
differential.activity$motif <- rownames(differential.activity)
readr::write_csv(differential.activity, "chromvar_motifs/Ex_chromvar.csv")

jpeg("images/chromvar_Ex_atac_markers.jpeg", units="in", width=10, height=10, res=300)
MotifPlot(
  object = Seurat.object.NEW,
  motifs = head(rownames(differential.activity)),
  assay = 'peaks'
) + theme(text = element_text(size = 18))
dev.off()
