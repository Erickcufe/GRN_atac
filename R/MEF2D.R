
rorb_uniques <- readr::read_csv("rorb_uniques.csv")
region_MEF2D <- rorb_uniques$promoter[rorb_uniques$TF=="MEF2D"]
motif_name <- "MEF2D"
region <- region_MEF2D
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

MEF2D_snps <- results_1


# MEF2D_snps <- search_motifBreak(motif_name = "MEF2D", region = region_MEF2D)
dir.create("RORB/MEF2D")
saveRDS(MEF2D_snps, "RORB/MEF2D/MEF2D_snps.rds")
MEF2D_snps <- readRDS("RORB/MEF2D/MEF2D_snps.rds")
MEF2D_snps <- subset(MEF2D_snps, dataSource %in% c("jaspar2022", "HOCOMOCOv11-core-A"))
MEF2D_snps$Refpvalue <- NULL
MEF2D_snps$Altpvalue <- NULL

df_MEF2D <- data.frame(MEF2D_snps)
df_MEF2D <- df_MEF2D[df_MEF2D$effect=="strong",]
table(df_MEF2D$SNP_id)
plot_SNPs(MEF2D_snps, dir_path = "RORB/MEF2D")
