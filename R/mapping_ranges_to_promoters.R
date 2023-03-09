mapping_ranges_to_promoters <- function(credranges){
  suppressMessages(library(GenomicRanges))
  suppressMessages(library(AnnotationHub))

  load("promoter_ranges.rda")

  p1 <- promoterranges[seqnames(promoterranges)==1,]
  p2 <- promoterranges[seqnames(promoterranges)==2,]
  p3 <- promoterranges[seqnames(promoterranges)==3,]
  p4 <- promoterranges[seqnames(promoterranges)==4,]
  p5 <- promoterranges[seqnames(promoterranges)==5,]
  p6 <- promoterranges[seqnames(promoterranges)==6,]
  p7 <- promoterranges[seqnames(promoterranges)==7,]
  p8 <- promoterranges[seqnames(promoterranges)==8,]
  p9 <- promoterranges[seqnames(promoterranges)==9,]
  p10 <- promoterranges[seqnames(promoterranges)==10,]
  p11 <- promoterranges[seqnames(promoterranges)==11,]
  p12 <- promoterranges[seqnames(promoterranges)==12,]
  p13 <- promoterranges[seqnames(promoterranges)==13,]
  p14 <- promoterranges[seqnames(promoterranges)==14,]
  p15 <- promoterranges[seqnames(promoterranges)==15,]
  p16 <- promoterranges[seqnames(promoterranges)==16,]
  p17 <- promoterranges[seqnames(promoterranges)==17,]
  p18 <- promoterranges[seqnames(promoterranges)==18,]
  p19 <- promoterranges[seqnames(promoterranges)==19,]
  p20 <- promoterranges[seqnames(promoterranges)==20,]
  p21 <- promoterranges[seqnames(promoterranges)==21,]
  p22 <- promoterranges[seqnames(promoterranges)==22,]
  px <- promoterranges[seqnames(promoterranges)=="X",]
  py <- promoterranges[seqnames(promoterranges)=="Y",]

  promoter_total <- c(p1, p2, p3, p4 ,p5,p6,p7,p8,p9,p10,p11,p12,p13,
                      p14,p15,p16,p17,p18,p19,p20,p21,p22,px,py)

  #change seq levels
  # seqlevels(promoter_total, pruning.mode="coarse") <- unique(seqnames(promoter_total))
  # seqlevels(promoter_total) <- paste0("chr", seqlevels(promoter_total))

  main.chroms <- standardChromosomes(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)
  seqlevels(promoterranges) <- paste0("chr", seqlevels(promoterranges))
  keep.peaks <- as.logical(seqnames(promoterranges) %in% main.chroms)
  promoterranges <- promoterranges[keep.peaks, ]
  seqlevels(promoterranges) <- main.chroms

  seqlevels(credranges) <- seqlevels(promoterranges)
  genome(credranges) <- genome(promoter_total)

  # Overlap credible SNPs with promoter regions.
  olap <- IRanges::findOverlaps(credranges, promoter_total)
  credpromoter <- credranges[S4Vectors::queryHits(olap)]

  rangos <-  data.frame(ranges(promoterranges[S4Vectors::subjectHits(olap)]))
  mcols(credpromoter) <- cbind(mcols(credpromoter),
                               mcols(promoterranges[S4Vectors::subjectHits(olap)]),
                               rangos)
  df <- data.frame(credpromoter)
  df <- na.omit(df)
  df <- df[!duplicated(df$start),]

  table(df$gene_name)
  summary(df$gene_name)
  genes <- df$gene_name
  gen_sel <- genes[genes != ""]
  gen_sel <- gen_sel[!duplicated(gen_sel)]
  message("This are the genes that match the Promoter region/ATAC region")
  print(gen_sel)
  return(credpromoter)

}


## El plan es mapear los fragmentos de atac seq con las regiones promotoras de los genes
credranges <- readRDS("ranges_ATAC_annoted.rds")


atac_rorb$hg38_Start <- as.numeric(atac_rorb$hg38_Start)
atac_rorb$hg38_Stop <- as.numeric(atac_rorb$hg38_Stop)
rorb_granges <- GRanges(seqnames = atac_rorb$hg38_Chromosome,
                        ranges = IRanges(start = atac_rorb$hg38_Start,
                                         end = atac_rorb$hg38_Stop))


# #change seq levels
# seqlevels(promoter_total, pruning.mode="coarse") <- unique(seqnames(promoter_total))
# seqlevels(promoter_total) <- paste0("chr", seqlevels(promoter_total))
#
# seqlevels(credranges) <- seqlevels(promoter_total)
# genome(credranges) <- genome(promoter_total)

# Overlap credible SNPs with promoter regions.
olap <- IRanges::findOverlaps(credranges, rorb_granges)
credpromoter <- credranges[S4Vectors::queryHits(olap)]

rangos <-  data.frame(ranges(promoterranges[S4Vectors::subjectHits(olap)]))
mcols(credpromoter) <- cbind(mcols(credpromoter),
                             mcols(promoterranges[S4Vectors::subjectHits(olap)]),
                             rangos)
df <- data.frame(credpromoter)
df <- na.omit(df)
df <- df[!duplicated(df$start),]

table(df$gene_name)
summary(df$gene_name)
genes <- df$gene_name
gen_sel <- genes[genes != ""]
gen_sel <- gen_sel[!duplicated(gen_sel)]
message("This are the genes that match the Promoter region/ATAC region")
print(gen_sel)
return(credpromoter)


mapping_ranges_to_promoters(rorb_granges)



results <- mapping_ranges_to_promoters(range_atac)
