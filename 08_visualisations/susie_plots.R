library(ggplot2)
library(patchwork)
library(ggrepel)
library(locuszoomr)
library(EnsDb.Hsapiens.v86)

pip_all <- read.csv("04_susie_finemapping/results/all_variant_pip.csv")
pip_cs <- read.csv("04_susie_finemapping/results/credible_set_variants.csv")
pip_signals <- read.csv("04_susie_finemapping/results/conditional_signal_pip.csv")

pip_all$pos_mb <- pip_all$POS / 1e6
pip_all$in_cs <- pip_all$ID %in% pip_cs$ID

# gene track only, no association data needed
loc_genes <- locus(
  data = NULL,
  xrange = c(128930000, 129060000),
  seqname = "7",
  ens_db = "EnsDb.Hsapiens.v86"
)

pip_signals$label <- paste0(pip_signals$signal, "\n(", pip_signals$ID, ")")

p_pip <- ggplot(pip_all, aes(x = pos_mb, y = PIP)) +
  geom_point(data = subset(pip_all, !in_cs), shape = 1, size = 1.5, colour = "grey50") +
  geom_point(data = subset(pip_all, in_cs), shape = 16, size = 2, colour = "#E31A1C") +
  geom_hline(yintercept = 0.1, linetype = "dashed", linewidth = 0.4, colour = "black") +
  geom_text_repel(data = subset(pip_signals, signal != "signal2"),
                  aes(x = POS / 1e6, y = PIP, label = label),
                  size = 2.5, nudge_y = 0.02,
                  segment.size = 0.3, min.segment.length = 0,
                  inherit.aes = FALSE) +
  geom_text_repel(data = subset(pip_signals, signal == "signal2"),
                  aes(x = POS / 1e6, y = PIP, label = label),
                  size = 2.5, nudge_y = 0.05,
                  segment.size = 0.3, min.segment.length = 0,
                  inherit.aes = FALSE) +
  scale_x_continuous(limits = c(128.930, 129.060)) +
  scale_y_continuous(limits = c(0, 0.30)) +
  labs(x = NULL, y = "PIP") +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank())

p_genes <- gg_genetracks(loc_genes,
                         filter_gene_name = c("IRF5", "TNPO3"),
                         italics = TRUE) +
  labs(x = "Position (Mb, GRCh38 chr7)")

fig <- p_pip / p_genes + plot_layout(heights = c(4, 1))
print(fig)

ggsave("08_visualisations/figures/susie/susie_pip.pdf", fig, width = 7, height = 5.5)
ggsave("08_visualisations/figures/susie/susie_pip.png", fig, width = 7, height = 5.5, dpi = 300)