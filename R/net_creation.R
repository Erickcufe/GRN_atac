####################################
############### Inputs  ############
####################################
df <- readRDS("TFs_input_df.rds")

save(df, atac_input, extract_TFs, GRNet_2TF, file = "image_for_networks.RData")

load("image_for_networks.RData")

# df <- readRDS("promoters_TF.rds")
# df <- df[!is.na(df$promoter_TF), ]
library(dplyr)
library(Seurat)
# snRNA
so_morabito <- readRDS("../merge_neurons/anotacion_Parcial_neuronas_neuronType.rds")
so_morabito$tags <- Idents(so_morabito)

## RORB
# cell_AD_rorb <- na.omit(colnames(so_morabito)[(so_morabito$disease=="AD" & so_morabito$tags=="RORB")])
# cell_ctrl_rorb <- na.omit(colnames(so_morabito)[(so_morabito$disease=="Control" & so_morabito$tags=="RORB")])
cells_rorb <- na.omit(colnames(so_morabito)[so_morabito$tags=="RORB"])

so_morabito_RORB <- so_morabito[ ,cells_rorb]
# so_morabito_RORB_AD <- so_morabito_RORB[,cell_AD_rorb]
# so_morabito_RORB_ctrl <- so_morabito_RORB[,cell_ctrl_rorb]

atac_rorb <- readr::read_csv("atac_cellType_markers/RORB+_atac.csv") %>%
  atac_input()

## Ex
# cell_AD_ex <- na.omit(colnames(so_morabito)[(so_morabito$disease=="AD" & so_morabito$tags=="Ex_1")])
# cell_ctrl_ex <- na.omit(colnames(so_morabito)[(so_morabito$disease=="Control" & so_morabito$tags=="Ex_1")])
cells_ex <- na.omit(colnames(so_morabito)[so_morabito$tags=="Ex_1"])

so_morabito_Ex <- so_morabito[ ,cells_ex]
# so_morabito_Ex_AD <- so_morabito_Ex[,cell_AD_ex]
# so_morabito_Ex_ctrl <- so_morabito_Ex[,cell_ctrl_ex]

atac_ex <- readr::read_csv("atac_cellType_markers/Ex_atac.csv") %>%
  atac_input()

## Pv
cell_AD_pv <- na.omit(colnames(so_morabito)[(so_morabito$disease=="AD" & so_morabito$tags=="Pv")])
cell_ctrl_pv <- na.omit(colnames(so_morabito)[(so_morabito$disease=="Control" & so_morabito$tags=="Pv")])
cells_Pv <- na.omit(colnames(so_morabito)[so_morabito$tags=="Pv"])

so_morabito_Pv <- so_morabito[ ,cells_Pv]
so_morabito_Pv_AD <- so_morabito_Pv[,cell_AD_pv]
so_morabito_Pv_ctrl <- so_morabito_Pv[,cell_ctrl_pv]

atac_pv <- readr::read_csv("atac_cellType_markers/Pv_atac.csv") %>%
  atac_input()

## Sst
cell_AD_Sst <- na.omit(colnames(so_morabito)[(so_morabito$disease=="AD" & so_morabito$tags=="Sst")])
cell_ctrl_Sst <- na.omit(colnames(so_morabito)[(so_morabito$disease=="Control" & so_morabito$tags=="Sst")])
cells_Sst <- na.omit(colnames(so_morabito)[so_morabito$tags=="Sst"])

so_morabito_Sst <- so_morabito[ ,cells_Sst]
so_morabito_Sst_AD <- so_morabito_Pv[, cell_AD_Sst]
so_morabito_Sst_ctrl <- so_morabito_Pv[,cell_ctrl_Sst]

atac_sst <- readr::read_csv("atac_cellType_markers/Sst_atac.csv") %>%
  atac_input()

## Vip
cell_AD_Vip <- na.omit(colnames(so_morabito)[(so_morabito$disease=="AD" & so_morabito$tags=="Vip")])
cell_ctrl_Vip <- na.omit(colnames(so_morabito)[(so_morabito$disease=="Control" & so_morabito$tags=="Vip")])
cells_Vip <- na.omit(colnames(so_morabito)[so_morabito$tags=="Vip"])

so_morabito_Vip <- so_morabito[ ,cells_Vip]
so_morabito_Vip_AD <- so_morabito_Pv[,cell_AD_Vip]
so_morabito_Vip_ctrl <- so_morabito_Pv[,cell_ctrl_Vip]

atac_vip <- readr::read_csv("atac_cellType_markers/Vip_atac.csv") %>%
  atac_input()

## Non-Vip
cell_AD_Non_Vip <- na.omit(colnames(so_morabito)[(so_morabito$disease=="AD" & so_morabito$tags=="Non-Vip")])
cell_ctrl_Non_Vip <- na.omit(colnames(so_morabito)[(so_morabito$disease=="Control" & so_morabito$tags=="Non-Vip")])
cells_Non_Vip <- na.omit(colnames(so_morabito)[so_morabito$tags=="Non-Vip"])

so_morabito_Non_Vip <- so_morabito[ ,cells_Non_Vip]
so_morabito_Non_Vip_AD <- so_morabito_Pv[,cell_AD_Non_Vip]
so_morabito_Non_Vip_ctrl <- so_morabito_Pv[,cell_ctrl_Non_Vip]

atac_non_vip <- readr::read_csv("atac_cellType_markers/Non-Vip_atac.csv") %>%
  atac_input()


## Inhibitory
cells_Inh <- na.omit(colnames(so_morabito)[so_morabito$tags=="Non-Vip" |  so_morabito$tags=="Vip" | so_morabito$tags=="Non-Vip" | so_morabito$tags=="Sst"])
so_morabito_Inh <- so_morabito[ ,cells_Inh]

atac_inh <- rbind(atac_non_vip, atac_pv, atac_sst, atac_pv)


####################################
############### Outputs  ###########
####################################

## RORB
RORB_net <- GRNet_2TF(df = df, gexpr = so_morabito_RORB@assays$RNA@counts, open_chrom = atac_rorb)
RORB_AD_net <- GRNet_2TF(df = df, gexpr = so_morabito_RORB_AD@assays$RNA@counts, open_chrom = atac_rorb)
RORB_Ctrl_net <- GRNet_2TF(df = df, gexpr = so_morabito_RORB_ctrl@assays$RNA@counts, open_chrom = atac_rorb)

## Ex
Ex_net <- GRNet_2TF(df = df, gexpr = so_morabito_Ex@assays$RNA@counts, open_chrom = atac_ex)
Ex_AD_net <- GRNet_2TF(df = df, gexpr = so_morabito_Ex_AD@assays$RNA@counts, open_chrom = atac_ex)
Ex_Ctrl_net <- GRNet_2TF(df = df, gexpr = so_morabito_Ex_ctrl@assays$RNA@counts, open_chrom = atac_ex)

## Pv
Pv_net <- GRNet_2TF(df = df, gexpr = so_morabito_Pv@assays$RNA@counts, open_chrom = atac_pv)
Pv_AD_net <- GRNet_2TF(df = df, gexpr = so_morabito_Pv_AD@assays$RNA@counts, open_chrom = atac_pv)
Pv_Ctrl_net <- GRNet_2TF(df = df, gexpr = so_morabito_Pv_ctrl@assays$RNA@counts, open_chrom = atac_pv)

## Sst
Sst_net <- GRNet_2TF(df = df, gexpr = so_morabito_Sst@assays$RNA@counts, open_chrom = atac_sst)
Sst_AD_net <- GRNet_2TF(df = df, gexpr = so_morabito_Sst_AD@assays$RNA@counts, open_chrom = atac_sst)
Sst_Ctrl_net <- GRNet_2TF(df = df, gexpr = so_morabito_Sst_ctrl@assays$RNA@counts, open_chrom = atac_sst)

## Vip
Vip_net <- GRNet_2TF(df = df, gexpr = so_morabito_Vip@assays$RNA@counts, open_chrom = atac_vip)
Vip_AD_net <- GRNet_2TF(df = df, gexpr = so_morabito_Vip_AD@assays$RNA@counts, open_chrom = atac_vip)
Vip_Ctrl_net <- GRNet_2TF(df = df, gexpr = so_morabito_Vip_ctrl@assays$RNA@counts, open_chrom = atac_vip)

## Non-Vip
Non_Vip_net <- GRNet_2TF(df = df, gexpr = so_morabito_Non_Vip@assays$RNA@counts, open_chrom = atac_non_vip)
Non_Vip_AD_net <- GRNet_2TF(df = df, gexpr = so_morabito_Non_Vip_AD@assays$RNA@counts, open_chrom = atac_non_vip)
Non_Vip_Ctrl_net <- GRNet_2TF(df = df, gexpr = so_morabito_Non_Vip_ctrl@assays$RNA@counts, open_chrom = atac_non_vip)


## Inhibitory

inh_net <- GRNet_2TF(df = df, gexpr = so_morabito_Inh@assays$RNA@counts, open_chrom = atac_inh)


save(RORB_net,Ex_net, Pv_net,
     Vip_net, Non_Vip_net, inh_net, file = "Nets/CellType_Networks.RData")

####################################
############### Igraph  ###########
####################################

library(igraph)
library(dplyr)

g_inh <- inh_net %>%
  select(TF, TG, coef, promoter, mse) %>%
  mutate(rmse = sqrt(mse)) %>%
  filter(rmse < mean(rmse)) %>%
  graph_from_data_frame(directed = TRUE)
write_graph(g_inh, file = "Nets/Graph_INH_SFG.graphml", format = "graphml")

inh_net <- inh_net %>%
  select(TF, TG, coef, promoter, mse) %>%
  mutate(rmse = sqrt(mse)) %>%
  filter(rmse < mean(rmse))
inh_net$edge_btw <- edge_betweenness(g_inh)

save(g_inh, inh_net, file = "Nets/CellTypes_Inh_Networks.RData")

g_rorb <- RORB_net %>%
  select(TF, TG, coef, promoter, mse) %>%
  mutate(rmse = sqrt(mse)) %>%
  filter(rmse < mean(rmse)) %>%
  graph_from_data_frame(directed = TRUE)
write_graph(g_rorb, file = "Nets/Graph_RORB_comparacion_SFG.graphml",format = "graphml")

RORB_net <- RORB_net %>%
  select(TF, TG, coef, promoter, mse) %>%
  mutate(rmse = sqrt(mse)) %>%
  filter(rmse < mean(rmse))
RORB_net$edge_btw <- edge_betweenness(g_rorb)

g_ex <- Ex_net %>%
  select(TF, TG, coef, promoter, mse) %>%
  mutate(rmse = sqrt(mse)) %>%
  filter(rmse <= mean(rmse)) %>%
  graph_from_data_frame(directed = TRUE)
write_graph(g_ex, file = "Nets/Graph_Ex_SFG.graphml",format = "graphml")

Ex_net <- Ex_net %>%
  select(TF, TG, coef, promoter, mse) %>%
  mutate(rmse = sqrt(mse)) %>%
  filter(rmse <= mean(rmse))
Ex_net$edge_btw <- edge_betweenness(g_ex)


g_pv <- Pv_net %>%
  select(TF, TG, coef, promoter, mse) %>%
  mutate(rmse = sqrt(mse)) %>%
  filter(rmse <= mean(rmse)) %>%
  graph_from_data_frame(directed = TRUE)
write_graph(g_pv, file = "Nets/Graph_Pv_contraste_SFG.graphml",format = "graphml")

Pv_net <- Pv_net %>%
  select(TF, TG, coef, promoter, mse) %>%
  mutate(rmse = sqrt(mse)) %>%
  filter(rmse <= mean(rmse))
Pv_net$edge_btw <- edge_betweenness(g_pv)


g_sst <- Sst_net %>%
  select(TF, TG, coef, promoter, mse) %>%
  mutate(rmse = sqrt(mse)) %>%
  filter(rmse <= mean(rmse)) %>%
  graph_from_data_frame(directed = TRUE)
write_graph(g_sst, file = "Nets/Graph_Sst_SFG.graphml",format = "graphml")

Sst_net <- Sst_net %>%
  select(TF, TG, coef, promoter, mse) %>%
  mutate(rmse = sqrt(mse)) %>%
  filter(rmse <= mean(rmse))
Sst_net$edge_btw <- edge_betweenness(g_sst)


g_vip <- Vip_net %>%
  select(TF, TG, coef, promoter, mse) %>%
  mutate(rmse = sqrt(mse)) %>%
  filter(rmse <= mean(rmse)) %>%
  graph_from_data_frame(directed = TRUE)
write_graph(g_vip, file = "Nets/Graph_Vip_SFG.graphml",format = "graphml")

Vip_net <- Vip_net %>%
  select(TF, TG, coef, promoter, mse) %>%
  mutate(rmse = sqrt(mse)) %>%
  filter(rmse <= mean(rmse))
Vip_net$edge_btw <- edge_betweenness(g_vip)


g_non_vip <- Non_Vip_net %>%
  select(TF, TG, coef, promoter, mse) %>%
  mutate(rmse = sqrt(mse)) %>%
  filter(rmse <= mean(rmse)) %>%
  graph_from_data_frame(directed = TRUE)
write_graph(g_non_vip, file = "Nets/Graph_Non_Vip_SFG.graphml",format = "graphml")

Non_Vip_net <- Non_Vip_net %>%
  select(TF, TG, coef, promoter, mse) %>%
  mutate(rmse = sqrt(mse)) %>%
  filter(rmse <= mean(rmse))
Non_Vip_net$edge_btw <- edge_betweenness(g_non_vip)

save(Non_Vip_net, g_non_vip, Vip_net, g_vip, Pv_net, g_pv,
     Ex_net, g_ex, RORB_net, g_rorb, file = "Nets/CellType_Networks.RData")
