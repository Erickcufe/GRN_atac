library(ggplot2)
library(TFBSTools)
library(Seurat)
Seurat.object.NEW <- readRDS("motifs_ATAC.rds")

differential.activity <- FindMarkers(
  object = Seurat.object.NEW,
  ident.1 = 'Ex',
  ident.2 = NULL,
  only.pos = TRUE,
  test.use = 'LR',
  min.pct = 0.05,
  latent.vars = 'nCount_peaks'
)

# get top differentially accessible peaks
top.da.peak <- rownames(differential.activity[differential.activity$p_val < 0.005, ])
enriched.motifs_RORB_Ex <- FindMotifs(
  object = Seurat.object.NEW,
  features = top.da.peak
)

differential.activity <- differential.activity[differential.activity$p_val_adj < 0.05,]
differential.activity$regions <- rownames(differential.activity)

readr::write_csv(differential.activity, "atac_cellType_markers/Ex_atac.csv")

#Plot de los motifs mas representativos
jpeg("images/Ex_atac_markers.jpeg", units="in", width=10, height=10, res=300)
MotifPlot(
  object = Seurat.object.NEW,
  motifs = head(rownames(enriched.motifs_RORB_Ex))
) + theme(text = element_text(size = 25))

dev.off()


## Pv
differential.activity <- FindMarkers(
  object = Seurat.object.NEW,
  ident.1 = 'Pv',
  ident.2 = NULL,
  only.pos = TRUE,
  test.use = 'LR',
  min.pct = 0.05,
  latent.vars = 'nCount_peaks'
)

# get top differentially accessible peaks
top.da.peak <- rownames(differential.activity[differential.activity$p_val < 0.005, ])
enriched.motifs_RORB_Ex <- FindMotifs(
  object = Seurat.object.NEW,
  features = top.da.peak
)

differential.activity <- differential.activity[differential.activity$p_val_adj < 0.05,]
differential.activity$regions <- rownames(differential.activity)

readr::write_csv(differential.activity, "atac_cellType_markers/Pv_atac.csv")

#Plot de los motifs mas representativos
jpeg("images/Pv_atac_markers.jpeg", units="in", width=10, height=10, res=300)
MotifPlot(
  object = Seurat.object.NEW,
  motifs = head(rownames(enriched.motifs_RORB_Ex))
) + theme(text = element_text(size = 25))

dev.off()


## Sst
differential.activity <- FindMarkers(
  object = Seurat.object.NEW,
  ident.1 = 'Sst',
  ident.2 = NULL,
  only.pos = TRUE,
  test.use = 'LR',
  min.pct = 0.05,
  latent.vars = 'nCount_peaks'
)

# get top differentially accessible peaks
top.da.peak <- rownames(differential.activity[differential.activity$p_val < 0.005, ])
enriched.motifs_RORB_Ex <- FindMotifs(
  object = Seurat.object.NEW,
  features = top.da.peak
)

differential.activity <- differential.activity[differential.activity$p_val_adj < 0.05,]
differential.activity$regions <- rownames(differential.activity)

readr::write_csv(differential.activity, "atac_cellType_markers/Sst_atac.csv")

#Plot de los motifs mas representativos
jpeg("images/Sst_atac_markers.jpeg", units="in", width=10, height=10, res=300)
MotifPlot(
  object = Seurat.object.NEW,
  motifs = head(rownames(enriched.motifs_RORB_Ex))
) + theme(text = element_text(size = 25))

dev.off()

## Vip
differential.activity <- FindMarkers(
  object = Seurat.object.NEW,
  ident.1 = 'Vip',
  ident.2 = NULL,
  only.pos = TRUE,
  test.use = 'LR',
  min.pct = 0.05,
  latent.vars = 'nCount_peaks'
)

# get top differentially accessible peaks
top.da.peak <- rownames(differential.activity[differential.activity$p_val < 0.005, ])
enriched.motifs_RORB_Ex <- FindMotifs(
  object = Seurat.object.NEW,
  features = top.da.peak
)

differential.activity <- differential.activity[differential.activity$p_val_adj < 0.05,]
differential.activity$regions <- rownames(differential.activity)

readr::write_csv(differential.activity, "atac_cellType_markers/Vip_atac.csv")

#Plot de los motifs mas representativos
jpeg("images/Vip_atac_markers.jpeg", units="in", width=10, height=10, res=300)
MotifPlot(
  object = Seurat.object.NEW,
  motifs = head(rownames(enriched.motifs_RORB_Ex))
) + theme(text = element_text(size = 25))

dev.off()


## Non-Vip
differential.activity <- FindMarkers(
  object = Seurat.object.NEW,
  ident.1 = 'Non-Vip',
  ident.2 = NULL,
  only.pos = TRUE,
  test.use = 'LR',
  min.pct = 0.05,
  latent.vars = 'nCount_peaks'
)

# get top differentially accessible peaks
top.da.peak <- rownames(differential.activity[differential.activity$p_val < 0.005, ])
enriched.motifs_RORB_Ex <- FindMotifs(
  object = Seurat.object.NEW,
  features = top.da.peak
)

differential.activity <- differential.activity[differential.activity$p_val_adj < 0.05,]
differential.activity$regions <- rownames(differential.activity)

readr::write_csv(differential.activity, "atac_cellType_markers/Non_vip_atac.csv")

#Plot de los motifs mas representativos
jpeg("images/Non_vip_atac_markers.jpeg", units="in", width=10, height=10, res=300)
MotifPlot(
  object = Seurat.object.NEW,
  motifs = head(rownames(enriched.motifs_RORB_Ex))
) + theme(text = element_text(size = 25))

dev.off()

