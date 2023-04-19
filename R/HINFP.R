
rorb_uniques <- readr::read_csv("rorb_uniques.csv")
region_HINFP <- rorb_uniques$promoter[rorb_uniques$TF=="HINFP"]
motif_name <- "HINFP"
region <- region_HINFP
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

results_1 <- motifbreakR::motifbreakR(snpList = snps.mb, filterp = TRUE,
                                      pwmList = motif_list,
                                      threshold = 1e-4,
                                      method = "ic",
                                      BPPARAM = BiocParallel::MulticoreParam(6))

HINFP_snps <- results_1


# HINFP_snps <- search_motifBreak(motif_name = "HINFP", region = region_HINFP)
dir.create("RORB/HINFP")
saveRDS(HINFP_snps, "RORB/HINFP/HINFP_snps.rds")
HINFP_snps <- readRDS("RORB/HINFP/HINFP_snps.rds")
HINFP_snps <- subset(HINFP_snps, dataSource %in% c("jaspar2022", "HOCOMOCOv11-core-A"))
HINFP_snps$Refpvalue <- NULL
HINFP_snps$Altpvalue <- NULL

df_HINFP <- data.frame(HINFP_snps)
df_HINFP <- df_HINFP[df_HINFP$effect=="strong",]
table(df_HINFP$SNP_id)
plot_SNPs(HINFP_snps, dir_path = "RORB/HINFP")
