library(ggplot2)
library(patchwork)

# load coloc data
irf5_eqtl <- read.csv("06_colocalisation/eqtl_catalogue/results/coloc_summary_irf5.csv")
tnpo3_eqtl <- read.csv("06_colocalisation/eqtl_catalogue/results/coloc_summary_tnpo3.csv")
gtex <- read.csv("06_colocalisation/gtex/results/coloc_summary_gtex_irf5.csv")

# add GTEx as a row in IRF5 dataset
gtex_row <- data.frame(dataset_id = "GTEx", dataset_name = "GTEx whole blood",
                       n_variants = gtex$n_overlapping, PP.H0 = gtex$PP.H0,
                       PP.H1 = gtex$PP.H1, PP.H2 = gtex$PP.H2,
                       PP.H3 = gtex$PP.H3, PP.H4 = gtex$PP.H4)
irf5_all <- rbind(irf5_eqtl, gtex_row)

# generate coloc summary heatmaps for each gene
irf5_all$label <- ifelse(irf5_all$PP.H4 > 0.5,
                         paste0(round(irf5_all$PP.H4, 2)),
                         paste0(round(irf5_all$PP.H4, 2)))
irf5_all$text_col <- ifelse(irf5_all$PP.H4 > 0.5, "white", "black")

tnpo3_eqtl$text_col <- ifelse(tnpo3_eqtl$PP.H4 > 0.5, "white", "black")
tnpo3_eqtl$label <- ifelse(tnpo3_eqtl$dataset_name == "GTEx whole blood",
                           paste0("PP.H4=", round(tnpo3_eqtl$PP.H4, 2),
                                  "\nPP.H3=", round(tnpo3_eqtl$PP.H3, 2)),
                           as.character(round(tnpo3_eqtl$PP.H4, 2)))

p_irf5 <- ggplot(irf5_all, aes(x = "IRF5", y = dataset_name, fill = PP.H4)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = label, colour = text_col), size = 2.8) +
  scale_fill_gradient(low = "grey90", high = "#8B0000", limits = c(0, 1), name = "PP.H4") +
  scale_colour_identity() +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 10) +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(face = "bold"))

p_tnpo3 <- ggplot(tnpo3_eqtl, aes(x = "TNPO3", y = dataset_name, fill = PP.H4)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = label, colour = text_col), size = 2.8) +
  scale_fill_gradient(low = "grey90", high = "#8B0000", limits = c(0, 1), name = "PP.H4") +
  scale_colour_identity() +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 10) +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(face = "bold"))

print(p_irf5)
print(p_tnpo3)

ggsave("08_visualisations/figures/coloc/coloc_heatmap_irf5.pdf", p_irf5, width = 4, height = 5.5)
ggsave("08_visualisations/figures/coloc/coloc_heatmap_irf5.png", p_irf5, width = 4, height = 5.5, dpi = 300)
ggsave("08_visualisations/figures/coloc/coloc_heatmap_tnpo3.pdf", p_tnpo3, width = 4, height = 2.5)
ggsave("08_visualisations/figures/coloc/coloc_heatmap_tnpo3.png", p_tnpo3, width = 4, height = 2.5, dpi = 300)

# generate GTEx IRF5 per-SNP PP.H4 plot
library(ggrepel)
library(locuszoomr)
library(EnsDb.Hsapiens.v86)

loc_genes <- locus(
  data = NULL,
  xrange = c(128930000, 129060000),
  seqname = "7",
  ens_db = "EnsDb.Hsapiens.v86"
)

snp_pp <- read.csv("06_colocalisation/gtex/results/coloc_snp_posteriors_gtex_irf5.csv")
snp_pp$pos_mb <- snp_pp$position / 1e6
snp_pp$in_promoter <- snp_pp$position >= 128930000 & snp_pp$position <= 128947119

lead_coloc <- snp_pp[which.max(snp_pp$SNP.PP.H4), ]
lead_coloc$label <- "rs3778754"

p_snp <- ggplot(snp_pp, aes(x = pos_mb, y = SNP.PP.H4, colour = in_promoter)) +
  geom_point(size = 2, alpha = 0.85) +
  geom_text_repel(data = lead_coloc,
                  aes(x = pos_mb, y = SNP.PP.H4, label = label),
                  size = 3.2, nudge_y = 0.03,
                  segment.size = 0.3, inherit.aes = FALSE) +
  scale_colour_manual(values = c("TRUE" = "#E31A1C", "FALSE" = "grey50"),
                      labels = c("TRUE" = "IRF5 promoter", "FALSE" = "Rest of region"),
                      name = NULL) +
  scale_x_continuous(limits = c(128.930, 129.060)) +
  scale_y_continuous(limits = c(0, 0.6)) +
  labs(x = NULL, y = "SNP PP.H4") +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "right")

p_genes_snp <- gg_genetracks(loc_genes,
                             filter_gene_name = c("IRF5", "TNPO3"),
                             italics = TRUE) +
  labs(x = "Position (Mb, GRCh38 chr7)")

fig <- p_snp / p_genes_snp + plot_layout(heights = c(4, 1))
print(fig)

ggsave("08_visualisations/figures/coloc/coloc_snp_pp.pdf", fig, width = 7, height = 5.5)
ggsave("08_visualisations/figures/coloc/coloc_snp_pp.png", fig, width = 7, height = 5.5, dpi = 300)