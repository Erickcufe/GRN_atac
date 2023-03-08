atac_input <- function(atac){
  vect_split <- sapply(atac$regions, function(x){
    a <- stringr::str_split(x, "-")
  })
  df_atac <- data.frame()
  for(i in 1:length(vect_split)){
    df_tmp <- t(data.frame(unlist(vect_split[[i]])))
    colnames(df_tmp) <- c("hg38_Chromosome", "hg38_Start", "hg38_Stop")
    df_atac <- rbind(df_atac, df_tmp)
  }
  return(df_atac)
}
