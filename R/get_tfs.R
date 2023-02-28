library(Seurat)
library(SingleCellExperiment)
library(dplyr)
nuc <- readRDS("SingleR_anotated_ATAC_NEW.rds")
Idents(nuc) <- nuc$Proplabels
DefaultAssay(nuc) <- "RNA"

sce_reference <- as.SingleCellExperiment(nuc, assay = "ATAC")
colData(sce_reference) <- as.data.frame(colData(sce_reference)) %>%
  mutate_if(is.character, as.factor) %>%
  DataFrame(row.names = colnames(sce_reference))

aa <- as.Seurat(ex_neu)
aa_aa <- subset(x = aa, subset =  nCount_ATAC > 3)


ex_neu <- sce_reference[,sce_reference$ident=="Ex"]

# new.data <- GetAssayData(object = nuc, slot = "counts", assay = "ATAC")
# nuc[["peaks"]] <- CreateAssayObject( new.data)
# nuc[,WhichCells(object = nuc, idents = "Ex")]
# subset(x = nuc, subset = Proplabels == "Ex")
DefaultAssay(nuc) <- "ATAC"
nuc_neuronEx <- nuc[,Idents(nuc)=="Ex"]
ranges_annoted <- nuc@assays$ATAC@annotation
ranges_annoted <- ranges_annoted[seqnames(ranges_annoted) != "chrX", ]
ranges_annoted <- ranges_annoted[seqnames(ranges_annoted) != "chrY", ]
ranges_annoted <- ranges_annoted[seqnames(ranges_annoted) != "chrM", ]

saveRDS(ranges_annoted, "ranges_ATAC_annoted.rds")

scGRNom_getTF <- function(df, database = JASPAR2022::JASPAR2022, species_type = 9606, min_score = 0.9,
                          pwm_type = 'prob',num_cores = 6){

  load("../CongresoNacionalBioquimica_Merida/data/GOAL_standard.rda")

  library(JASPAR2022)

  opts <- list()
  opts[["species"]] <- species_type
  opts[["all_versions"]] <- TRUE
  PFMatrixList <- TFBSTools::getMatrixSet(JASPAR2022::JASPAR2022, opts)
  pwmlist <- TFBSTools::toPWM(PFMatrixList, type = pwm_type)
  TF_names <- TFBSTools::name(pwmlist)
  names(TF_names) = NULL

  TF_names_splited = sapply(TF_names,data.table::tstrsplit,'::|\\(var.2\\)|\\(var.3\\)')

  df$promoter_id <- paste(df$gene_chr,':',df$promoter_start,'-',df$promoter_end,sep = '')
  df$enhancer_id <- paste(df$enh_chr,':',df$enh_start,'-',df$enh_end,sep = '')
  df = data.table::data.table(df)

  df_p <- dplyr::distinct(df[,c('gene_chr','promoter_start','promoter_end','promoter_id')])
  df_e <- dplyr::distinct(df[,c('enh_chr','enh_start','enh_end','enhancer_id')])

  df_p$gene_chr <- paste0("chr", df_p$gene_chr)
  df_e$enh_chr <- paste0("chr", df_e$enh_chr)
  suppressWarnings( G1 <- GenomicRanges::GRanges(seqnames = df_p$gene_chr,
                                                 IRanges::IRanges(start=df_p$promoter_start,
                                                                  end=df_p$promoter_end),
                                                 seqlengths = GenomeInfoDb::seqlengths(BSgenome.Hsapiens.UCSC.hg38::Hsapiens)[1:24]))
  G1 <- GenomicRanges::trim(G1)
  suppressWarnings( G2 <- GenomicRanges::GRanges(seqnames = df_e$enh_chr,
                                                 IRanges::IRanges(start=df_e$enh_start,
                                                                  end=df_e$enh_end),
                                                 seqlengths = GenomeInfoDb::seqlengths(BSgenome.Hsapiens.UCSC.hg38::Hsapiens)[1:24]))
  G2 <- GenomicRanges::trim(G2)


  cl <- parallel::makeCluster(num_cores) # not overload your computer
  doParallel::registerDoParallel(cl)

  `%dopar%` <- foreach::`%dopar%`
  df_p$promoter_TF <- foreach::foreach(i = 1:nrow(df_p), .combine = rbind,
                                       .packages = c('data.table','motifmatchr')) %dopar% {
                                         peak <- G1[i]
                                         motif_ix <- matchMotifs(pwmlist, peak,
                                                                 genome = "hg38",
                                                                 out = "scores"
                                         )
                                         result <- motifScores(motif_ix)[1,]
                                         curr_TF <- unique(unlist(TF_names_splited[result > quantile(result,min_score)]))
                                         if(length(curr_TF) == 0){
                                           curr_TF <- NA
                                         }
                                         data.table(promoter_TF = list(curr_TF))

                                       }
  parallel::stopCluster(cl)

  saveRDS(df_p, "promoters_TF.rds")

  cl <- parallel::makeCluster(num_cores)
  doParallel::registerDoParallel(cl)
  df_e$enhancer_TF <- foreach::foreach(i = 1:nrow(df_e), .combine = rbind,
                                       .packages = c('data.table','motifmatchr')) %dopar% {
                                         peak <- G2[i]
                                         motif_ix <- matchMotifs(pwmlist, peak,
                                                                 genome = "hg38",
                                                                 out = "scores")
                                         result <- motifScores(motif_ix)[1,]
                                         curr_TF <- unique(unlist(TF_names_splited[result > quantile(result,min_score)]))
                                         if(length(curr_TF) == 0){
                                           curr_TF <- NA
                                         }
                                         data.table(promoter_TF = list(curr_TF))
                                       }
  parallel::stopCluster(cl)


  df$promoter_TF <- df_p$promoter_TF[match(df$promoter_id, df_p$promoter_id)]
  df$enhancer_TF <- df_e$enhancer_TF[match(df$enhancer_id, df_e$enhancer_id)]

  df <- df[, c('gene','promoter_id','enhancer_id',
               'promoter_TF','enhancer_TF')]
  colnames(df) <- c('gene','promoter','enhancer',
                    'promoter_TF','enhancer_TF')
  return(df)


}
