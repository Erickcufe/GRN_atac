
rorb_uniques <- readr::read_csv("rorb_uniques.csv")
region_ZNF354A <- rorb_uniques$promoter[rorb_uniques$TF=="ZNF354A"]
motif_name <- "ZNF354A"
region <- region_ZNF354A
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


ZNF354A_snps <- search_motifBreak(motif_name = "ZNF354A", region = region_ZNF354A)
saveRDS(ZNF354A_snps, "RORB/ZNF354A/ZNF354A_snps.rds")
ZNF354A_snps <- readRDS("RORB/ZNF354A/ZNF354A_snps.rds")
ZNF354A_snps <- subset(ZNF354A_snps, dataSource %in% c("jaspar2022", "HOCOMOCOv11-core-A"))
ZNF354A_snps$Refpvalue <- NULL
ZNF354A_snps$Altpvalue <- NULL

df_ZNF354A <- data.frame(ZNF354A_snps)
table(df_ZNF354A$SNP_id)
plot_SNPs(ZNF354A_snps = ZNF354A_snps, dir_path = "RORB/ZNF354A")
