
rorb_uniques <- readr::read_csv("rorb_uniques.csv")
region_MAFF <- rorb_uniques$promoter[rorb_uniques$TF=="MAFF"]
motif_name <- "MAFF"
region <- region_MAFF
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

MAFF_snps <- results_1


# MAFF_snps <- search_motifBreak(motif_name = "MAFF", region = region_MAFF)
dir.create("RORB/MAFF")
saveRDS(MAFF_snps, "RORB/MAFF/MAFF_snps.rds")
MAFF_snps <- readRDS("RORB/MAFF/MAFF_snps.rds")
MAFF_snps <- subset(MAFF_snps, dataSource %in% c("jaspar2022", "HOCOMOCOv11-core-A"))
MAFF_snps$Refpvalue <- NULL
MAFF_snps$Altpvalue <- NULL

df_MAFF <- data.frame(MAFF_snps)
df_MAFF <- df_MAFF[df_MAFF$effect=="strong",]
table(df_MAFF$SNP_id)
plot_SNPs(MAFF_snps, dir_path = "RORB/MAFF")
