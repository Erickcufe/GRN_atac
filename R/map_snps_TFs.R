snps <- readRDS("all_snps_motifs_strong.rds")

search_overlap <- function(net){

  suppressMessages(library(GenomicRanges))
  suppressMessages(library(dplyr))
  snps <- readRDS("all_snps_motifs_strong.rds")

  s_ranges <- GRanges(seqnames = snps$seqnames,
                      IRanges(start = snps$start,
                              end = snps$start),
                      SNP = snps$SNP_id)

  chrs <- sapply(net$promoter, function(x){
    a <- stringr::str_split(x, ":")
    chrs <- a[[1]][[1]]
  }) %>% unname()

  seqranges <- sapply(net$promoter, function(x){
    a <- stringr::str_split(x, ":")
    seqranges <- a[[1]][[2]]
  }) %>% unname()

  seqranges <- stringr::str_split(seqranges, "-")

  starts <- sapply(seqranges, function(x){
    x[[1]]
  })

  ends <- sapply(seqranges, function(x){
    x[[2]]
  })

  p_ranges <- GRanges(seqnames = chrs,
                      IRanges(start = as.numeric(starts),
                              end = as.numeric(ends)),
                      TF = net$TF,
                      TG = net$TG)

  olap <- IRanges::findOverlaps(s_ranges, p_ranges)
  credpromoter <- s_ranges[S4Vectors::queryHits(olap)]

  rangos <-  data.frame(ranges(p_ranges[S4Vectors::subjectHits(olap)]))
  mcols(credpromoter) <- cbind(mcols(credpromoter),
                               mcols(p_ranges[S4Vectors::subjectHits(olap)]),
                               rangos)

}


## Posiciones claves de RORB+

chr4_174828789_174829804 <- snpsByOverlaps(snps_ids, GRanges("4:174828789-174829804"))
chr4_174828789_174829804 <- data.frame(chr4_174828789_174829804)

chr13_61414782_61415825 <- snpsByOverlaps(snps_ids, GRanges("13:61414782-61415825"))
chr13_61414782_61415825 <- data.frame(chr13_61414782_61415825)

chr5_114361595_114362564 <- snpsByOverlaps(snps_ids, GRanges("5:114361595-114362564"))
chr5_114361595_114362564 <- data.frame(chr5_114361595_114362564)

chr5_50440539_50441773 <- snpsByOverlaps(snps_ids, GRanges("5:50440539-50441773"))
chr5_50440539_50441773 <- data.frame(chr5_50440539_50441773)

chr5_49651964_49656963 <- snpsByOverlaps(snps_ids, GRanges("5:49651964-49656963"))
chr5_49651964_49656963 <- data.frame(chr5_49651964_49656963)

chr1_159141438_159142463 <- snpsByOverlaps(snps_ids, GRanges("1:159141438-159142463"))
chr1_159141438_159142463 <- data.frame(chr1_159141438_159142463)

library(motifbreakR)

ids <- prueba$RefSNP_id
snps.mb <- motifbreakR::snps.from.rsid(rsid = ids,
                                       dbSNP = SNPlocs.Hsapiens.dbSNP155.GRCh38::SNPlocs.Hsapiens.dbSNP155.GRCh38,
                                       search.genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)

data(motifbreakR_motif)
results <- motifbreakR::motifbreakR(snpList = snps.mb, filterp = TRUE,
                                    pwmList = motifbreakR_motif,
                                    threshold = 1e-4,
                                    method = "ic",
                                    BPPARAM = BiocParallel::MulticoreParam(6))





ids <- prueba$RefSNP_id
snps.mb <- motifbreakR::snps.from.rsid(rsid = ids,
                                       dbSNP = SNPlocs.Hsapiens.dbSNP155.GRCh38::SNPlocs.Hsapiens.dbSNP155.GRCh38,
                                       search.genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)

data(motifbreakR_motif)
results <- motifbreakR::motifbreakR(snpList = snps.mb, filterp = TRUE,
                                    pwmList = motifbreakR_motif,
                                    threshold = 1e-4,
                                    method = "ic",
                                    BPPARAM = BiocParallel::MulticoreParam(6))

r_strong <- results[results$effect=="strong",]
r_strong <- data.frame(r_strong)

rs1001158 <- results[names(results) %in% "rs1001158:A"]
rs1001158
rs1001158 <- calculatePvalue(rs1001158)
rs1001158[rs1001158$Altpvalue<=0.05]


saveRDS(results_pvalue, "results/motifs_all_AD_SNPs_PVALUE.rds")
rs13032148 <- results[names(results) %in% "rs1644757728"]

rs13032148 <- calculatePvalue(rs13032148)
rs13032148 <- data.frame(rs13032148)


library(BSgenome.Hsapiens.UCSC.hg38)

plotMB(results = results, rsid = "rs1644757728", effect = "strong")
