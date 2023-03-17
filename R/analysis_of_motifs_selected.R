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
                       geneSymbol == motif_name)

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
  results <- c(results_1, results_2)

  return(results)
}



Seurat.object <- readRDS("footprint_TFs.rds")

# Motifs: FOSL2, KLF5, KLF16, PATZ1, TBR1, ZNF148

##################################
## Posiciones claves de FOSL2 en RORB+
##################################
load("Ex/posciones_snps_motifbreak.rda")
load("Nets/CellType_Networks.RData")
load("RORB/posciones_snps_motifbreak.rda")
motifs_IDs <- Seurat.object@assays$peaks@motifs@motif.names
motifs_list <- Seurat.object@assays$peaks@motifs@positions
scores_pos_FOSL2 <- motifs_list$MA0478.1 %>% data.frame()

region_fosl2 <- RORB_net$promoter[RORB_net$TF=="FOSL2"]
FOSL2_snps <- search_motifBreak(motif_name = "FOSL2", region = region_fosl2)
saveRDS(FOSL2_snps, "RORB/FOSL2_snps.rds")
FOSL2_snps <- readRDS("RORB/FOSL2_snps.rds")

df_fosl2 <- data.frame(FOSL2_snps)
table(df_fosl2$SNP_id)

library(BSgenome.Hsapiens.UCSC.hg38)

plotMB(results = FOSL2_snps, rsid = "rs1021192569", effect = "strong")





PlotFootprint(Seurat.object, features = c("FOSL2"),
              idents = c("Pv", "RORB+","Ex", "Vip", "Non-Vip", "Sst")
)

