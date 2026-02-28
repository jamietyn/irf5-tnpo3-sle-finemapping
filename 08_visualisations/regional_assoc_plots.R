library(locuszoomr)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(patchwork)
library(ggrepel)

ld_vars <- read.table(
  "03_regional_association/results/ld_matrices/irf5_tnpo3_ld_r2.unphased.vcor2.vars",
  header = FALSE, col.names = "ID"
)
ld_mat <- as.matrix(read.table(
  "03_regional_association/results/ld_matrices/irf5_tnpo3_ld_r2.unphased.vcor2",
  header = FALSE
))

cond_signals <- read.table(
  "03_regional_association/results/conditional/cond_signals.txt",
  header = FALSE, col.names = "ID"
)
lead_r0 <- cond_signals$ID[1]
lead_r1 <- cond_signals$ID[2]
lead_r2 <- cond_signals$ID[3]

make_locus <- function(assoc_file, lead_id) {
  assoc <- read.table(assoc_file, header = TRUE, sep = "\t", comment.char = "")
  colnames(assoc)[1] <- "CHROM"
  assoc$ld <- ld_mat[match(assoc$ID, ld_vars$ID), which(ld_vars$ID == lead_id)]
  locus(
    data = assoc,
    xrange = c(128930000, 129060000),
    seqname = "7",
    ens_db = "EnsDb.Hsapiens.v86",
    chrom = "CHROM",
    pos = "POS",
    p = "P",
    labs = "ID",
    LD = "ld"
  )
}

make_fig <- function(loc, lead_id, conditioned_on = NULL,
                     pcutoff_line, ylim_top) {
  
  p_scatter <- gg_scatter(
    loc,
    index_snp = lead_id,
    labels = "index",
    pcutoff = Inf,
    size = 2,
    nudge_y = 0.3,
    ylim = c(0, ylim_top),
    xticks = FALSE,
    segment.alpha = 0
  ) +
    geom_hline(yintercept = -log10(pcutoff_line),
               linetype = "dashed", colour = "black", linewidth = 0.6) +
    theme(legend.position = "right")
  
  if (!is.null(conditioned_on)) {
    cond_data <- loc$data[loc$data$ID %in% conditioned_on, ]
    p_scatter <- p_scatter +
      geom_point(data = cond_data,
                 aes(x = POS / 1e6, y = -log10(P)),
                 shape = 4, size = 2.5, colour = "black",
                 inherit.aes = FALSE) +
      geom_text_repel(data = cond_data,
                      aes(x = POS / 1e6, y = -log10(P), label = ID),
                      size = 2.2, colour = "black", nudge_y = 0.3,
                      inherit.aes = FALSE)
  }
  
  p_genes <- gg_genetracks(
    loc,
    filter_gene_name = c("IRF5", "TNPO3"),
    italics = TRUE
  ) + labs(x = "Position (Mb, GRCh38 chr7)")
  
  p_scatter / p_genes + plot_layout(heights = c(5, 1))
}

loc_r0 <- make_locus(
  "03_regional_association/results/conditional/irf5_tnpo3_round0_sorted.txt",
  lead_r0
)
loc_r1 <- make_locus(
  "03_regional_association/results/conditional/irf5_tnpo3_round1_sorted.txt",
  lead_r1
)
loc_r2 <- make_locus(
  "03_regional_association/results/conditional/irf5_tnpo3_round2_sorted.txt",
  lead_r2
)

rd0 <- make_fig(loc_r0, lead_r0, NULL, 0.05 / 444, 7.0)
rd1 <- make_fig(loc_r1, lead_r1, lead_r0,0.01, 5.0)
rd2 <- make_fig(loc_r2, lead_r2, c(lead_r0, lead_r1), 0.01, 4.0)

print(rd0)
print(rd1)
print(rd2)

# save as pdf and png
ggsave("08_visualisations/figures/regional_assoc/fig_round0.pdf", rd0,  width = 7, height = 5.5)
ggsave("08_visualisations/figures/regional_assoc/fig_round0.png", rd0,  width = 7, height = 5.5, dpi = 300)
ggsave("08_visualisations/figures/regional_assoc/fig_round1.pdf", rd1, width = 7, height = 5.5)
ggsave("08_visualisations/figures/regional_assoc/fig_round1.png", rd1, width = 7, height = 5.5, dpi = 300)
ggsave("08_visualisations/figures/regional_assoc/fig_round2.pdf", rd2, width = 7, height = 5.5)
ggsave("08_visualisations/figures/regional_assoc/fig_round2.png", rd2, width = 7, height = 5.5, dpi = 300)