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

##################################
## Posiciones claves de RORB+
##################################

chr4_174828789_174829804 <- SNPlocs.Hsapiens.dbSNP155.GRCh38::snpsByOverlaps(snps_ids, GRanges("4:174828789-174829804"))
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

search_motifBreak <- function(ids){
  snps.mb <- motifbreakR::snps.from.rsid(rsid = ids,
                                         dbSNP = SNPlocs.Hsapiens.dbSNP155.GRCh38::SNPlocs.Hsapiens.dbSNP155.GRCh38,
                                         search.genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)

  data(motifbreakR_motif)
  results <- motifbreakR::motifbreakR(snpList = snps.mb, filterp = TRUE,
                                      pwmList = motifbreakR_motif,
                                      threshold = 1e-4,
                                      method = "ic",
                                      BPPARAM = BiocParallel::MulticoreParam(6))

  return(results)
}

list_positions <- list(chr1_159141438_159142463 = chr1_159141438_159142463$RefSNP_id,
                       chr5_49651964_49656963 = chr5_49651964_49656963$RefSNP_id,
                       chr5_50440539_50441773 = chr5_50440539_50441773$RefSNP_id,
                       chr5_114361595_114362564 = chr5_114361595_114362564$RefSNP_id,
                       chr13_61414782_61415825 = chr13_61414782_61415825$RefSNP_id,
                       chr4_174828789_174829804 = chr4_174828789_174829804$RefSNP_id)

results_snps_rorb <- purrr::map(list_positions, search_motifBreak)

save(results_snps_rorb, file = "RORB/posciones_snps_motifbreak.rda")


##################################
## Posiciones claves de Ex
##################################

chr1_53738678_53739782 <- snpsByOverlaps(snps_ids, GRanges("1:53738678-53739782"))
chr1_53738678_53739782 <- data.frame(chr1_53738678_53739782)

chr11_34639677_34644676 <- snpsByOverlaps(snps_ids, GRanges("11:34639677-34644676"))
chr11_34639677_34644676 <- data.frame(chr11_34639677_34644676)

chr11_7017847_7022846 <- snpsByOverlaps(snps_ids, GRanges("11:7017847-7022846"))
chr11_7017847_7022846 <- data.frame(chr11_7017847_7022846)

chr11_7017736_7022735 <- snpsByOverlaps(snps_ids, GRanges("11:7017736-7022735"))
chr11_7017736_7022735 <- data.frame(chr11_7017736_7022735)

chr3_25427787_25428711 <- snpsByOverlaps(snps_ids, GRanges("3:25427787-25428711"))
chr3_25427787_25428711 <- data.frame(chr3_25427787_25428711)

chr20_20366382_20368506 <- snpsByOverlaps(snps_ids, GRanges("20:20366382-20368506"))
chr20_20366382_20368506 <- data.frame(chr20_20366382_20368506)

chr12_20365604_20370603 <- snpsByOverlaps(snps_ids, GRanges("12:20365604-20370603"))
chr12_20365604_20370603 <- data.frame(chr12_20365604_20370603)

list_positions_ex <- list(chr12_20365604_20370603 = chr12_20365604_20370603$RefSNP_id,
                       chr20_20366382_20368506 = chr20_20366382_20368506$RefSNP_id,
                       chr3_25427787_25428711 = chr3_25427787_25428711$RefSNP_id,
                       chr11_7017736_7022735 = chr11_7017736_7022735$RefSNP_id,
                       chr11_7017847_7022846 = chr11_7017847_7022846$RefSNP_id,
                       chr11_34639677_34644676 = chr11_34639677_34644676$RefSNP_id,
                       chr1_53738678_53739782 = chr1_53738678_53739782$RefSNP_id)

results_snps_ex <- purrr::map(list_positions_ex, search_motifBreak)

save(results_snps_ex, file = "Ex/posciones_snps_motifbreak.rda")

########################
########################
########################

library(dplyr)
library(motifbreakR)
library(BSgenome.Hsapiens.UCSC.hg38)

load("Ex/posciones_snps_motifbreak.rda")
load("Nets/CellType_Networks.RData")
load("RORB/posciones_snps_motifbreak.rda")

chr4_174828789_174829804 <- results_snps_rorb[["chr4_174828789_174829804"]]$geneSymbol=="TBR1" %>%
  data.frame() %>% filter(effect == "strong")

results_snps_rorb[["chr4_174828789_174829804"]] <- na.omit(results_snps_rorb[["chr4_174828789_174829804"]])
results <- results_snps_rorb[["chr4_174828789_174829804"]]
results <- results[!is.na(elementMetadata(results)[,6])]
results <- results[results$geneSymbol %in% "TBR1"]
results <- results[2]

results <- calculatePvalue(results)
rs1286781862 <- results[names(results) %in% "rs1286781862"]
rs1286781862 <- calculatePvalue(rs1286781862)
results <- results[elementMetadata(results)]

jpeg("images/rs1286781862_TBR1.jpeg", units="in", width=10, height=10, res=300)
plotMB(results = results, rsid = "rs1286781862", effect = "strong")
dev.off()

library(clusterProfiler)
library(enrichplot)

library(org.Hs.eg.db)
AnnotationDbi::keytypes(org.Hs.eg.db)

gen_znf263 <- RORB_net[RORB_net$TF == "ZNF263",]
gen_znf263 <- gen_znf263$TG

# RORB
rorb_vs_norm_go <- clusterProfiler::enrichGO(RORB_net$TG,
                                             "org.Hs.eg.db",
                                             keyType = "SYMBOL",
                                             ont = "ALL",
                                             minGSSize = 50)

enr_go <- clusterProfiler::simplify(rorb_vs_norm_go)
View(enr_go@result)

jpeg("images/GO_TGs_RORB.jpeg", units="in", width=15, height=10, res=300)
enrichplot::emapplot(enrichplot::pairwise_termsim(enr_go),
                     showCategory = 30, cex_label_category = 0.7)
dev.off()



# Ex
ex_vs_norm_go <- clusterProfiler::enrichGO(Ex_net$TG,
                                             "org.Hs.eg.db",
                                             keyType = "SYMBOL",
                                             ont = "ALL",
                                             minGSSize = 50)

enr_go <- clusterProfiler::simplify(ex_vs_norm_go)
View(enr_go@result)

jpeg("images/GO_TGs_Ex.jpeg", units="in", width=15, height=10, res=300)
enrichplot::emapplot(enrichplot::pairwise_termsim(enr_go),
                     showCategory = 30, cex_label_category = 0.7)
dev.off()


########################
########################
########################

ids <- prueba$RefSNP_id




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

