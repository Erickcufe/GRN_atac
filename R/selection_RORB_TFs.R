TFs_uniques <- c("NFYC", "MEF2D", "GLIS1",
         "HMBOX1", "OVOL2", "ZNF354A",
         "ATF2", "ASCL1", "TAL1", "ONECUT2",
         "KLF15", "SOX15", "ZNF382",
         "FOXN3", "NFATC4", "CREB3",
         "MAFF", "HINFP", "PAX6",
         "ELK4", "IRF8", "TFAP2E",
         "TGIF2", "MXI1", "ZBTB14",
         "HNF4G")

Ex_RORB_tf <- c("YY1", "FOS", "NFATC3",
                "FOSL2", "MAFK", "TCF7L1",
                "REL", "POU2F1", "SOX12",
                "SOX8", "ETV1", "CEBPG",
                "BACH2", "TFE3", "NR2F1",
                "ETS2", "NFE2L1",
                "BHLHE41", "ETS1", "TCF12",
                "E2F3", "KLF7", "ETV3",
                "RARB", "ELF2", "BHLHE40",
                "GABPA", "TCF4", "TFEB",
                "NR2C2", "NFIB", "EGR1",
                "NFATC2", "FOSB", "NR2F2",
                "KLF16", "EMX1", "FOXG1",
                "NR2F6", "NFE2L2", "ARNT2",
                "KLF6", "E2F4", "ZBTB7A",
                "THRA", "RXRB")

RORB_all <- readr::read_csv("RORB_tfs_merge.csv")
TFs <- RORB_all$Tfs

rorb <- readRDS("Nets/Net_RORB_PAPER_SFG.rds")

rorb_sel <- rorb[rorb$TF %in% TFs,]
rorb_sel <- unique(rorb_sel)
rorb_sel$Is_marker <- 0

for(i in 1:nrow(rorb_sel)){
  for(j in 1:length(TFs_uniques)){
    if(rorb_sel[i,"TF"] == TFs_uniques[j]){
      rorb_sel$Is_marker[i] <- 1
    } else {
      next
    }
  }
}

library(igraph)
library(dplyr)

g_rorb <- rorb_sel %>%
  select(TF, TG, coef, promoter, mse, Is_marker) %>%
  graph_from_data_frame(directed = TRUE)

plot(g_rorb)
rorb_sel$edge_btw <- edge_betweenness(g_rorb)
degree_out <- data.frame(degree(g_rorb, mode = "out"))
degree_out$TF <- rownames(degree_out)
colnames(degree_out) <- c("degree_out", "TF")

degree_in <- data.frame(degree(g_rorb, mode = "in"))
degree_in$TF <- rownames(degree_in)
colnames(degree_in) <- c("degree_in", "TF")


write_graph(g_rorb, file = "Nets/Graph_RORB_SFG_TFs_common.graphml", format = "graphml")
readr::write_csv(rorb_sel, "Nets/Graph_rorb_PAPER_uniques.csv")

view_tfcommun_degree <- degree_out[degree_out$TF %in% TFs_uniques,]
view_tfcommun_degree <- view_tfcommun_degree[order(view_tfcommun_degree$degree_out),]

rorb_uniques <- rorb_sel[rorb_sel$TF %in% TFs_uniques,]
view_tfcommun_in <- degree_in[degree_in$TF %in% rorb_uniques$TG,]

jpeg("images/Degree_out_RORB.jpeg", units="in", width=15, height=10, res=300)
ggplot(view_tfcommun_degree, aes(x = reorder(TF, -degree_out), y = degree_out)) +
  geom_col(fill = "gray", color = "black") + theme_classic() +
  ylab("Degree out") + xlab("TF") +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 90))
dev.off()


jpeg("images/Degree_in_RORB.jpeg", units="in", width=18, height=10, res=300)
ggplot(view_tfcommun_in, aes(x = reorder(TF, -degree_in), y = degree_in)) +
  geom_col(fill = "gray", color = "black") + theme_classic() +
  ylab("Degree out") + xlab("Target gene") +
  theme(text = element_text(size = 16),
        axis.text.x = element_text(angle = 90))
dev.off()

rorb_uniques <- rorb_sel[rorb_sel$TF %in% TFs_uniques,]
readr::write_csv(rorb_uniques, "rorb_uniques.csv")
positions <- data.frame(table(rorb_uniques$promoter))
colnames(positions)[1] <- "promoter"

rorb_ <- right_join(rorb_uniques, positions, "promoter")
# Se uso la funcion search_overlap del archivo map_snps_TFs.R

snps_in_rorb <- search_overlap(rorb)
snps_in_rorb

# No hubo ninguna coincidencia con los SNPs asociados a Alzheimer en region promotora

###
###
###

#Ahora se realizará la busqueda de SNPs por motivo de unión

# se usó la funcion search_motifBreak del archivo analysis_of_motifs_selected

# ZNF354A

region_ZNF354A <- rorb_uniques$promoter[rorb_uniques$TF=="ZNF354A"]
ZNF354A_snps <- search_motifBreak(motif_name = "ZNF354A", region = region_ZNF354A)
saveRDS(ZNF354A_snps, "RORB/ZNF354A/ZNF354A_snps.rds")
ZNF354A_snps <- readRDS("RORB/ZNF354A/ZNF354A_snps.rds")
ZNF354A_snps <- subset(ZNF354A_snps, dataSource %in% c("jaspar2022", "HOCOMOCOv11-core-A"))
ZNF354A_snps$Refpvalue <- NULL
ZNF354A_snps$Altpvalue <- NULL

df_ZNF354A <- data.frame(ZNF354A_snps)
table(df_ZNF354A$SNP_id)
plot_SNPs(ZNF354A_snps = ZNF354A_snps, dir_path = "RORB/ZNF354A")

