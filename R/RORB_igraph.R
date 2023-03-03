library(igraph)
library(dplyr)
snps_motifs <- readRDS("all_snps_motifs_strong.rds")

grph_rorb <- readr::read_csv("Nets/RORB_net.csv")
grph_rorb_ad <- readr::read_csv("Nets/RORB_AD_net.csv")
grph_rorb_ct <- readr::read_csv("Nets/RORB_control_net.csv")


g <- graph_from_data_frame(grph_rorb, directed = TRUE)

plot(g)

write_graph(g, file = "Nets/Graph_RORB_SFG.graphml",format = "graphml")

degrees_ <- igraph::degree(g, mode = "in")
degrees_ <- degrees_[order(degrees_)]
