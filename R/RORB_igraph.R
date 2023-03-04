library(igraph)
library(dplyr)
snps_motifs <- readRDS("all_snps_motifs_strong.rds")

grph_rorb <- readr::read_csv("Nets/RORB_net.csv")
grph_rorb_ad <- readr::read_csv("Nets/RORB_AD_net.csv")
grph_rorb_ct <- readr::read_csv("Nets/RORB_control_net.csv")
grph_rorb_dbl <- readr::read_csv("Nets/double_net_RORB+.csv")

df_net <- data.frame(n1 = grph_rorb_dbl$TF,
                     n2 = grph_rorb_dbl$TG)

g <- graph_from_data_frame(df_net, directed = TRUE)

plot(g)

write_graph(g, file = "Nets/Graph_RORB_SFG.graphml",format = "graphml")

degrees_ <- igraph::degree(g, mode = "in")
degrees_ <- degrees_[order(degrees_)]

rb_attack <- brainGraph::robustness(g)
e_b <- igraph::edge.betweenness(g)

df_net$edge_bwt <- e_b
