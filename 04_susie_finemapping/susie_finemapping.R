# SuSiE fine-mapping of IRF5-TNPO3 locus in UKB WGS data
# Window: chr7:128,930,429-129,058,173 (444 variants)

library(susieR)
library(dplyr)

source("../utils/variant_definitions.R")

# Load association results (round 0, unconditional)
assoc <- read.table("../03_regional_association/results/conditional/irf5_tnpo3_round0_sorted.txt",
                    header = TRUE, stringsAsFactors = FALSE, comment.char = "")
colnames(assoc)[1] <- "CHROM"
cat("Variants in window:", nrow(assoc), "\n")

# Load signed r LD matrix (required by susie_rss)
ld_r        <- as.matrix(read.table("../03_regional_association/results/ld/irf5_tnpo3_ld_r.unphased.vcor1", header = FALSE))
variant_ids <- readLines("../03_regional_association/results/ld/irf5_tnpo3_ld_r.unphased.vcor1.vars")
rownames(ld_r) <- variant_ids
colnames(ld_r) <- variant_ids

cat("LD matrix:", dim(ld_r)[1], "x", dim(ld_r)[2], "\n")
cat("Symmetric:", all(abs(ld_r - t(ld_r)) < 1e-10), "\n")
cat("Diagonal all 1s:", all(abs(diag(ld_r) - 1) < 1e-10), "\n")

# Compute z-scores: z = log(OR) / SE
# Valid for large-sample GWAS (n=405,448)
assoc <- assoc %>%
  mutate(
    beta = log(OR),
    z    = beta / LOG.OR._SE
  )

# Align z-scores with LD matrix column order
assoc_ordered <- assoc[match(colnames(ld_r), assoc$ID), ]
z <- assoc_ordered$z
names(z) <- colnames(ld_r)

cat("Z-score range:", range(z, na.rm = TRUE), "\n")
cat("Missing z-scores:", sum(is.na(z)), "\n")

# Run SuSiE
# n=405,448 total samples; L=10; s parameter not passed (susie_rss does not natively handle case-control)
set.seed(42)
susie_fit <- susie_rss(
  z                          = z,
  R                          = ld_r,
  n                          = 405448,
  L                          = 10,
  estimate_residual_variance = TRUE,
  verbose                    = TRUE
)

cat("Credible sets:", length(susie_fit$sets$cs), "\n")
cat("Converged:", susie_fit$converged, "\n")

# Extract PIPs
pip_results <- assoc_ordered %>%
  mutate(PIP = susie_fit$pip) %>%
  arrange(desc(PIP))

# Extract credible set variants and summary for each CS
if (length(susie_fit$sets$cs) > 0) {
  cs_variants <- data.frame()
  for (cs_name in names(susie_fit$sets$cs)) {
    cs_idx <- susie_fit$sets$cs[[cs_name]]
    cs_block <- assoc_ordered[cs_idx, ]
    cs_block$PIP <- susie_fit$pip[cs_idx]
    cs_block$credible_set <- cs_name
    cs_variants <- rbind(cs_variants, cs_block)
  }
  cs_variants <- cs_variants %>% arrange(credible_set, desc(PIP))

  cs_summary <- data.frame()
  for (cs_name in names(susie_fit$sets$cs)) {
    cs_idx <- susie_fit$sets$cs[[cs_name]]
    pos    <- assoc_ordered$POS[cs_idx]
    cs_row <- data.frame(
      credible_set = cs_name,
      n_variants   = length(cs_idx),
      coverage     = susie_fit$sets$coverage[cs_name],
      purity_min_r = susie_fit$sets$purity[cs_name, "min.abs.corr"],
      pos_min      = min(pos),
      pos_max      = max(pos),
      pos_span_kb  = (max(pos) - min(pos)) / 1000
    )
    cs_summary <- rbind(cs_summary, cs_row)
  }

  cat("Credible set sizes:", paste(cs_summary$n_variants, collapse = ", "), "\n")
  print(cs_summary)
} else {
  cat("No credible sets identified\n")
  cs_variants <- data.frame()
  cs_summary  <- data.frame()
}

# PIPs of literature variants (defined in utils/variant_definitions.R)
lit_pip <- pip_results %>%
  filter(ID %in% lit_var) %>%
  mutate(
    rsid            = names(lit_var)[match(ID, lit_var)],
    in_credible_set = ID %in% cs_variants$ID
  ) %>%
  select(rsid, ID, POS, OR, P, PIP, in_credible_set) %>%
  arrange(desc(PIP))

print(lit_pip, row.names = FALSE)

# PIPs of 3 conditional signals (defined in utils/variant_definitions.R)
signal_pip <- pip_results %>%
  filter(ID %in% signals) %>%
  mutate(
    signal          = names(signals)[match(ID, signals)],
    in_credible_set = ID %in% cs_variants$ID
  ) %>%
  select(signal, ID, POS, OR, P, PIP, in_credible_set) %>%
  arrange(desc(PIP))

print(signal_pip, row.names = FALSE)

# Save outputs
write.csv(pip_results,  "results/all_variant_pip.csv",           row.names = FALSE)
write.csv(cs_variants,  "results/credible_set_variants.csv",     row.names = FALSE)
write.csv(cs_summary,   "results/credible_set_summary.csv",      row.names = FALSE)
write.csv(lit_pip,      "results/lit_variant_pip.csv",           row.names = FALSE)
write.csv(signal_pip,   "results/conditional_signal_pip.csv",    row.names = FALSE)
saveRDS(susie_fit,      "results/susie_fit.rds")

# Save log Bayes factors for potential future coloc.susie use
lbf_df <- as.data.frame(t(susie_fit$lbf_variable))
colnames(lbf_df) <- paste0("L", 1:ncol(lbf_df))
lbf_df$variant_id <- colnames(ld_r)
write.csv(lbf_df, "results/log_bayes_factors.csv", row.names = FALSE)