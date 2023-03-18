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


plot_SNPs <- function(FOSL2_snps, dir_path){

  FOSL2_snps <- subset(FOSL2_snps, effect == "strong")
  ids_snps <- unique(FOSL2_snps$SNP_id)
  suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
  suppressMessages(library(motifbreakR))


  for(i in 1:length(ids_snps)){
    a <- ids_snps[i]
    print(a)
    file_dir <- paste0(here::here(),"/",dir_path, "/", a, ".jpeg")
    jpeg(file_dir, units="in", width=10, height=10, res=300)
    plotMB(results = FOSL2_snps, rsid = a, effect = "strong")
    dev.off()
  }
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
saveRDS(FOSL2_snps, "RORB/FOSL2/FOSL2_snps.rds")
FOSL2_snps <- readRDS("RORB/FOSL2/FOSL2_snps.rds")
FOSL2_snps <- subset(FOSL2_snps, dataSource %in% c("jaspar2022", "HOCOMOCOv11-core-A"))
FOSL2_snps$Refpvalue <- NULL
FOSL2_snps$Altpvalue <- NULL

df_fosl2 <- data.frame(FOSL2_snps)
table(df_fosl2$SNP_id)
plot_SNPs(FOSL2_snps = FOSL2_snps, dir_path = "RORB/FOSL2")


##################################
## Posiciones claves de TRB1 en RORB+
##################################

region_tbr1 <- RORB_net$promoter[RORB_net$TF=="TBR1"]
TBR1_snps <- search_motifBreak(motif_name = "TBR1", region = region_tbr1)
saveRDS(TBR1_snps, "RORB/TBR1/TBR1_snps.rds")

TBR1_snps <- subset(TBR1_snps, dataSource %in% c("jaspar2022", "HOCOMOCOv11-core-A"))
TBR1_snps$Refpvalue <- NULL
TBR1_snps$Altpvalue <- NULL

plot_SNPs(FOSL2_snps = TBR1_snps, dir_path = "RORB/TBR1")
df <- data.frame(TBR1_snps)


##################################
## Posiciones claves de FOSL2 en Ex
##################################
region_fosl2 <- Ex_net$promoter[Ex_net$TF=="FOSL2"]
FOSL2_snps <- search_motifBreak(motif_name = "FOSL2", region = region_fosl2)
saveRDS(FOSL2_snps, "Ex/FOSL2/FOSL2_snps.rds")
FOSL2_snps <- subset(FOSL2_snps, dataSource %in% c("jaspar2022", "HOCOMOCOv11-core-A"))
FOSL2_snps$Refpvalue <- NULL
FOSL2_snps$Altpvalue <- NULL

plot_SNPs(FOSL2_snps = FOSL2_snps, dir_path = "Ex/FOSL2")



##################################
## Posiciones claves de TRB1 en Ex
##################################

region_tbr1 <- Ex_net$promoter[Ex_net$TF=="TBR1"]
TBR1_snps <- search_motifBreak(motif_name = "TBR1", region = region_tbr1)
saveRDS(TBR1_snps, "Ex/TBR1/TBR1_snps.rds")

TBR1_snps <- subset(TBR1_snps, dataSource %in% c("jaspar2022", "HOCOMOCOv11-core-A"))
TBR1_snps$Refpvalue <- NULL
TBR1_snps$Altpvalue <- NULL

plot_SNPs(FOSL2_snps = TBR1_snps, dir_path = "Ex/TBR1")
df <- data.frame(TBR1_snps)
