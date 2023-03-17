library(ggplot2)
library(dplyr)
search_motifBreak <- function(motif_name, region){

  suppressMessages(library(SNPlocs.Hsapiens.dbSNP155.GRCh38))
  suppressMessages(library(motifbreakR))
  suppressMessages(library(MotifDb))
  data(motifbreakR_motif)
  snps_ids <- SNPlocs.Hsapiens.dbSNP155.GRCh38::SNPlocs.Hsapiens.dbSNP155.GRCh38


  region <- stringr::str_remove(region, "chr")
  region <- snpsByOverlaps(snps_ids, GRanges(region))
  region <- data.frame(region)
  ids <- region$RefSNP_id

  snps.mb <- motifbreakR::snps.from.rsid(rsid = ids,
                                         dbSNP = SNPlocs.Hsapiens.dbSNP155.GRCh38::SNPlocs.Hsapiens.dbSNP155.GRCh38,
                                         search.genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)

  motif_list <- subset(MotifDb,
                       geneSymbol == motif_name, dataSource %in% c("jaspar2022"))

  motif_list_2 <- subset(motifbreakR_motif, geneSymbol == motif_name)
  results_1 <- motifbreakR::motifbreakR(snpList = snps.mb, filterp = TRUE,
                                      pwmList = motif_list,
                                      threshold = 1e-4,
                                      method = "ic",
                                      BPPARAM = BiocParallel::MulticoreParam(6))

  results_2 <- motifbreakR::motifbreakR(snpList = snps.mb, filterp = TRUE,
                                        pwmList = motif_list_2,
                                        threshold = 1e-4,
                                        method = "ic",
                                        BPPARAM = BiocParallel::MulticoreParam(6))

  return(results)
}

region_fosl2 <- RORB_net$promoter[RORB_net$TF=="FOSL2"]

FOSL2_snps <- search_motifBreak(motif_name = "FOSL2", region = region_fosl2)

Seurat.object <- readRDS("footprint_TFs.rds")

# Motifs: FOSL2, KLF5, KLF16, PATZ1, TBR1, ZNF148

##################################
## Posiciones claves de FOSL2 en RORB+
##################################

motif_list <- subset(motifbreakR_motif, geneSymbol == "FOSL2")

FOSL2 <- GetFootprintData(Seurat.object, features = "FOSL2", group.by = "disease")
motifs_IDs <- Seurat.object@assays$peaks@motifs@motif.names
motifs_list <- Seurat.object@assays$peaks@motifs@positions
scores_pos_FOSL2 <- motifs_list$MA0478.1 %>% data.frame()



PlotFootprint(Seurat.object, features = c("FOSL2"),
              idents = c("Pv", "RORB+","Ex", "Vip", "Non-Vip", "Sst")
)

