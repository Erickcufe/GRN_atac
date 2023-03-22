

# save(df, file = "tfs_general_atac.rda")
load("tfs_general_atac.rda")

atac_rorb <- readr::read_csv("atac_cellType_markers/RORB+_atac.csv") %>%
  atac_input()

