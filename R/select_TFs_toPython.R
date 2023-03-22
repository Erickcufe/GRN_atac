library(dplyr)
# In this script we incorporate the scATAC seq to select the TFs to the python pipile to GRN

# save(df, file = "tfs_general_atac.rda")
load("tfs_general_atac.rda")
load("image_for_networks.RData")

atac_rorb <- readr::read_csv("atac_cellType_markers/RORB+_atac.csv") %>%
  atac_input()

atac_ex <- readr::read_csv("atac_cellType_markers/Ex_atac.csv") %>%
  atac_input()

atac_vip <- readr::read_csv("atac_cellType_markers/Vip_atac.csv") %>%
  atac_input()

atac_pv <- readr::read_csv("atac_cellType_markers/Pv_atac.csv") %>%
  atac_input()

atac_sst <- readr::read_csv("atac_cellType_markers/Sst_atac.csv") %>%
  atac_input()

atac_non_vip <- readr::read_csv("atac_cellType_markers/Non-Vip_atac.csv") %>%
  atac_input()

atac_inh <- rbind(atac_non_vip, atac_pv, atac_sst, atac_pv)

select_TFs_by_ATAC <- function(open_chrom, df, extension_bps = 0){
  library(stringr)
  library(GenomicRanges)
  load("promoter_ranges.rda")
  open_chrom$hg38_Start <- as.numeric(open_chrom$hg38_Start)
  open_chrom$hg38_Stop <- as.numeric(open_chrom$hg38_Stop)

  atac_granges <- GRanges(seqnames = open_chrom$hg38_Chromosome,
                          ranges = IRanges(start = open_chrom$hg38_Start,
                                           end = open_chrom$hg38_Stop))

  main.chroms <- standardChromosomes(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)
  seqlevels(promoterranges) <- paste0("chr", seqlevels(promoterranges))
  keep.peaks <- as.logical(seqnames(promoterranges) %in% main.chroms)
  promoterranges <- promoterranges[keep.peaks, ]
  seqlevels(promoterranges) <- main.chroms
  olap <- IRanges::findOverlaps(atac_granges, promoterranges)
  credpromoter <- atac_granges[S4Vectors::queryHits(olap)]

  rangos <-  data.frame(ranges(promoterranges[S4Vectors::subjectHits(olap)]))
  mcols(credpromoter) <- cbind(mcols(credpromoter),
                               mcols(promoterranges[S4Vectors::subjectHits(olap)]),
                               rangos)

  df_ranges <- data.frame(credpromoter)
  df_ranges <- data.frame(gene = df_ranges$gene_name,
                          gene_chr = df_ranges$seqnames,
                          promoter_start = df_ranges$start,
                          promoter_end = df_ranges$end)
  df_ranges <- na.omit(df_ranges)
  df_ranges <- df_ranges[!duplicated(df_ranges$promoter_start),]
  df_ranges <- extract_TFs(df_ranges)


  df <- rbind(df, df_ranges)

  if(!is.null(open_chrom)){
    df$en_chrs = str_extract(df$promoter, "(?i)(?<=chr)\\d+")
    # Extract promoter start and end positions
    df$en_start = as.numeric(str_extract(df$promoter, "(?i)(\\d+){5}"))
    df$en_end = abs(as.numeric(str_extract(df$promoter, "(?i)-(\\d+){5}")))
    df$en_start_intend = df$en_start -extension_bps
    df$en_end_intend  = df$en_end + extension_bps
    df <- na.omit(df)
    colnames(open_chrom) = c('chrs', 'Start', 'Stop')
    open_chrom$chrs = gsub('chr','', open_chrom$chrs)
    open_chrom$Start <- as.numeric(open_chrom$Start)
    open_chrom$Stop <- as.numeric(open_chrom$Stop)

    library(GenomicRanges)
    en_scGRN = GRanges(seqnames = df$en_chrs,
                       IRanges(start = df$en_start_intend,
                               end = df$en_end_intend))

    en_openchrom = GRanges(seqnames = open_chrom$chrs,
                           IRanges(start = open_chrom$Start,
                                   end = open_chrom$Stop))

    overlap_en = as.data.frame(
      findOverlaps(
        query = en_scGRN,
        subject = en_openchrom,
        type = 'any',
        select = "all",
        ignore.strand = T
      )
    )

    unique_ens = unique(overlap_en$queryHits)
    df = df[unique_ens,]
  }

  df_tf <- unlist(df$promoter_TF)
  df_tf <- unique(df_tf)
  df_tf <- data.frame(tf_name = df_tf)

  return(df_tf)

}


rorb <- select_TFs_by_ATAC(open_chrom = atac_rorb, df = df, extension_bps = 0) %>%
  readr::write_csv("tf_names_rorb.csv")

ex <- select_TFs_by_ATAC(open_chrom = atac_ex, df = df, extension_bps = 0) %>%
  readr::write_csv("tf_names_ex.csv")

pv <- select_TFs_by_ATAC(open_chrom = atac_pv, df = df, extension_bps = 0) %>%
  readr::write_csv("tf_names_pv.csv")

sst <- select_TFs_by_ATAC(open_chrom = atac_sst, df = df, extension_bps = 0) %>%
  readr::write_csv("tf_names_sst.csv")

vip <- select_TFs_by_ATAC(open_chrom = atac_vip, df = df, extension_bps = 0) %>%
  readr::write_csv("tf_names_vip.csv")

non_vip <- select_TFs_by_ATAC(open_chrom = atac_non_vip, df = df, extension_bps = 0) %>%
  readr::write_csv("tf_names_non_vip.csv")

inh <- select_TFs_by_ATAC(open_chrom = atac_inh, df = df, extension_bps = 0) %>%
  readr::write_csv("tf_names_inh.csv")



