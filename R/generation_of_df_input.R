load("promoter_ranges.rda")
credranges <- readRDS("ranges_ATAC_annoted.rds")
suppressMessages(library(GenomicRanges))
suppressMessages(library(AnnotationHub))

#change seq levels
# seqlevels(promoter_total, pruning.mode="coarse") <- unique(seqnames(promoter_total))
# seqlevels(promoter_total) <- paste0("chr", seqlevels(promoter_total))

main.chroms <- standardChromosomes(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)
seqlevels(promoterranges) <- paste0("chr", seqlevels(promoterranges))
keep.peaks <- as.logical(seqnames(promoterranges) %in% main.chroms)
promoterranges <- promoterranges[keep.peaks, ]
seqlevels(promoterranges) <- main.chroms

seqlevels(credranges) <- seqlevels(promoterranges)
genome(credranges) <- genome(promoterranges)

# Overlap credible SNPs with promoter regions.
olap <- IRanges::findOverlaps(credranges, promoterranges)
credpromoter <- credranges[S4Vectors::queryHits(olap)]

rangos <-  data.frame(ranges(promoterranges[S4Vectors::subjectHits(olap)]))
mcols(credpromoter) <- cbind(mcols(credpromoter),
                             mcols(promoterranges[S4Vectors::subjectHits(olap)]),
                             rangos)
df <- data.frame(credpromoter)
df <- na.omit(df)
df <- df[!duplicated(df$start),]

df_ready <- extract_TFs(df)


