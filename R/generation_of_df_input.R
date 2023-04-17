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
df_atac <- df


# Select per gene the merging range
df_gen <- data.frame()
for(i in levels(df_atac$seqnames)){
  df_tmp <- df_atac[df_atac$seqnames == i, ]

  nsel <- df_tmp[order(df_tmp$end),]
  nsel$gene_name <- factor(nsel$gene_name)

  for(k in levels(nsel$gene_name)){
    tmp <- nsel[nsel$gene_name == k, ]
    tmp <- tmp[!duplicated(tmp$start.1),]
    # min_r <- min(tmp$start.1)
    # max_r <- max(tmp$end.1)
    # extr_data <- tmp[1, c(-1, -2, -3, -4)]

    tgther <- tmp
    # tgther <- cbind(range_gen, extr_data)
    df_gen <- rbind(df_gen, tgther)
  }
}

df_gen_1 <- df_gen[df_gen$width >= 500, ]

# df_gen <- df_gen[-(which.max(df_gen$width)), ]

df_atac_ready <- data.frame(gene = df_gen_1$gene_name,
                        gene_chr = df_gen_1$seqnames,
                        promoter_start = df_gen_1$start.1,
                        promoter_end = df_gen_1$end.1)
saveRDS(df_atac_ready, "df_atac_ready.rds")



df_ready_1 <- extract_TFs(df_atac_ready[1:15860, ], num_cores = 6)
df_ready_2 <- extract_TFs(df_atac_ready[15860:31721, ], num_cores = 6)


saveRDS(df_ready, "INPUT_df.rds")

###########################################
###########################################
###########################################
###########################################


df_2 <- readRDS("INPUT_df.rds")
load("tfs_general_atac.rda")

# Mapping the ranges to a promoters
df$gene <- factor(df$gene)

df_modified <- data.frame()
for(i in levels(df$gene)){
  df_ranges <- df[df$gene==i,]
  split_coords <- sapply(df_ranges$promoter, function(x){
    a <- stringr::str_split(x, ":")
    b <- a[[1]][[2]]
    c <- stringr::str_split(b, "-")
    n1 <- as.numeric(c[[1]][[1]])
    n2 <- as.numeric(c[[1]][[2]])
    d <- n2 - n1
  })
  sel <- which(split_coords <= 5000)
  df_ranges <- df_ranges[sel,]
  df_modified <- rbind(df_modified, df_ranges)
}
df_ranges <- stringr::str_split(df$promoter, ":")

df_3 <- rbind(df, df_2)
df_3 <- df_3[!duplicated(df_3$promoter), ]

saveRDS(df_3, "TFs_input_df.rds")
