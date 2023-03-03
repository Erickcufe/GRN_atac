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
  ident.2 = NULL,
  only.pos = TRUE,
  test.use = 'LR',
  min.pct = 0.05,
  latent.vars = 'nCount_peaks'
)

differential.activity <- differential.activity[differential.activity$p_val < 0.05,]
differential.activity$regions <- rownames(differential.activity)

readr::write_csv(differential.activity, "atac_cellType_markers/RORB+_atac.csv")


###############
##############
# Pv NEURONS
###############


DefaultAssay(Seurat.object.NEW) <- "peaks"
differential.activity <- FindMarkers(
  object = Seurat.object.NEW,
  ident.1 = 'Pv',
  ident.2 = NULL,
  only.pos = TRUE,
  test.use = 'LR',
  min.pct = 0.05,
  latent.vars = 'nCount_peaks'
)

differential.activity <- differential.activity[differential.activity$p_val < 0.05,]
differential.activity$regions <- rownames(differential.activity)

readr::write_csv(differential.activity, "atac_cellType_markers/Pv_atac.csv")


###############
##############
# Ex NEURONS
###############

differential.activity <- FindMarkers(
  object = Seurat.object.NEW,
  ident.1 = 'Ex',
  ident.2 = NULL,
  only.pos = TRUE,
  test.use = 'LR',
  min.pct = 0.05,
  latent.vars = 'nCount_peaks'
)

differential.activity <- differential.activity[differential.activity$p_val < 0.05,]
differential.activity$regions <- rownames(differential.activity)

readr::write_csv(differential.activity, "atac_cellType_markers/Ex_atac.csv")

##############
# Vip NEURONS
###############

differential.activity <- FindMarkers(
  object = Seurat.object.NEW,
  ident.1 = 'Vip',
  ident.2 = NULL,
  only.pos = TRUE,
  test.use = 'LR',
  min.pct = 0.05,
  latent.vars = 'nCount_peaks'
)

differential.activity <- differential.activity[differential.activity$p_val < 0.05,]
differential.activity$regions <- rownames(differential.activity)

readr::write_csv(differential.activity, "atac_cellType_markers/Vip_atac.csv")


##############
# Non-Vip NEURONS
###############
DefaultAssay(Seurat.object.NEW)

differential.activity <- FindMarkers(
  object = Seurat.object.NEW,
  ident.1 = 'Non-Vip',
  ident.2 = NULL,
  only.pos = TRUE,
  test.use = 'LR',
  min.pct = 0.05,
  latent.vars = 'nCount_peaks'
)

differential.activity <- differential.activity[differential.activity$p_val_adj < 0.05,]
differential.activity$regions <- rownames(differential.activity)

readr::write_csv(differential.activity, "atac_cellType_markers/Non-Vip_atac.csv")


##############
# Sst NEURONS
###############

differential.activity <- FindMarkers(
  object = Seurat.object.NEW,
  ident.1 = 'Sst',
  ident.2 = NULL,
  only.pos = TRUE,
  test.use = 'LR',
  min.pct = 0.05,
  latent.vars = 'nCount_peaks'
)

differential.activity <- differential.activity[differential.activity$p_val < 0.05,]
differential.activity$regions <- rownames(differential.activity)

readr::write_csv(differential.activity, "atac_cellType_markers/Sst_atac.csv")
