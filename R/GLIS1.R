
rorb_uniques <- readr::read_csv("rorb_uniques.csv")
region_GLIS1 <- rorb_uniques$promoter[rorb_uniques$TF=="GLIS1"]
motif_name <- "GLIS1"
region <- region_GLIS1
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

GLIS1_snps <- results_1


GLIS1_snps <- search_motifBreak(motif_name = "GLIS1", region = region_GLIS1)
dir.create("RORB/GLIS1")
saveRDS(GLIS1_snps, "RORB/GLIS1/GLIS1_snps.rds")
GLIS1_snps <- readRDS("RORB/GLIS1/GLIS1_snps.rds")
GLIS1_snps <- subset(GLIS1_snps, dataSource %in% c("jaspar2022", "HOCOMOCOv11-core-A"))
GLIS1_snps$Refpvalue <- NULL
GLIS1_snps$Altpvalue <- NULL

df_GLIS1 <- data.frame(GLIS1_snps)
df_GLIS1 <- df_GLIS1[df_GLIS1$effect=="strong",]
table(df_GLIS1$SNP_id)
plot_SNPs(GLIS1_snps, dir_path = "RORB/GLIS1")
