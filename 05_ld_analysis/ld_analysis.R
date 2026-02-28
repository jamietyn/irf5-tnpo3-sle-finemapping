# LD analysis for IRF5-TNPO3 locus
# Pairwise r^2 between conditional signals, literature variants, and haplotype subgroups

source("../utils/variant_definitions.R")

# Load signed r LD matrix
ld_r_raw    <- as.matrix(read.table("../03_regional_association/results/ld_matrices/irf5_tnpo3_ld_r.unphased.vcor1", header = FALSE))
variant_ids <- readLines("../03_regional_association/results/ld_matrices/irf5_tnpo3_ld_r.unphased.vcor1.vars")
rownames(ld_r_raw) <- variant_ids
colnames(ld_r_raw) <- variant_ids

cat("Matrix symmetric:", all(abs(ld_r_raw - t(ld_r_raw)) < 1e-10), "\n")
cat("Diagonal all 1s:",  all(abs(diag(ld_r_raw) - 1) < 1e-10), "\n")
ld_r <- ld_r_raw
cat("LD matrix dimensions:", dim(ld_r), "\n")

# Check all variants of interest are present in LD matrix
all_vars <- c(signals, lit_var)
missing  <- all_vars[!all_vars %in% colnames(ld_r)]
if (length(missing) > 0) {
  cat("WARNING: Missing from LD matrix:\n")
  print(missing)
  stop("Cannot proceed with missing variants")
}
cat("All variants present in LD matrix\n")

# Pairwise r^2 between 3 conditional signals
# Confirms statistical independence of signals
signal_ld <- ld_r[signals, signals]^2
rownames(signal_ld) <- names(signals)
colnames(signal_ld) <- names(signals)
signal_ld_df <- as.data.frame(signal_ld)
signal_ld_df <- cbind(variant_id = unname(signals), signal = rownames(signal_ld_df), signal_ld_df)
rownames(signal_ld_df) <- NULL

cat("\nPairwise r² between conditional signals\n")
print(round(signal_ld_df[, -(1:2)], 3))
write.csv(signal_ld_df, "results/ld_signals_pairwise.csv", row.names = FALSE)

# r^2 between each signal and literature variants
# Maps signals onto published haplotypes; identifies potentially novel signals
signal_vs_lit <- ld_r[signals, lit_var]^2
rownames(signal_vs_lit) <- names(signals)
colnames(signal_vs_lit) <- names(lit_var)
signal_vs_lit_df <- as.data.frame(signal_vs_lit)
signal_vs_lit_df <- cbind(variant_id = unname(signals), signal = rownames(signal_vs_lit_df), signal_vs_lit_df)
rownames(signal_vs_lit_df) <- NULL

cat("\nr² between signals and literature variants\n")
print(round(signal_vs_lit_df[, -(1:2)], 3))
write.csv(signal_vs_lit_df, "results/ld_signals_vs_literature.csv", row.names = FALSE)

# IRF5 promoter haplotype LD structure
# Resolves whether rs142738614, rs4728142, rs3757387, rs2004640, rs10954213, rs35000415 tag the same signal
promoter_vars <- lit_var[c("rs142738614", "rs4728142", "rs3757387", "rs2004640", "rs10954213", "rs35000415")]

promoter_ld <- ld_r[promoter_vars, promoter_vars]^2
rownames(promoter_ld) <- names(promoter_vars)
colnames(promoter_ld) <- names(promoter_vars)
promoter_ld_df <- as.data.frame(promoter_ld)
promoter_ld_df <- cbind(variant_id = unname(promoter_vars), rsid = rownames(promoter_ld_df), promoter_ld_df)
rownames(promoter_ld_df) <- NULL

cat("\nIRF5 promoter haplotype LD structure\n")
print(round(promoter_ld_df[, -(1:2)], 3))
write.csv(promoter_ld_df, "results/ld_irf5_promoter_haplotype.csv", row.names = FALSE)

# TNPO3 haplotype LD structure
# Determines if signal1 tags the rs12534421/rs10488631/rs13239597 haplotype
tnpo3_vars <- c(
  signals["signal1"],
  lit_var[c("rs12534421", "rs10488631", "rs13239597")]
)

tnpo3_ld <- ld_r[tnpo3_vars, tnpo3_vars]^2
rownames(tnpo3_ld) <- names(tnpo3_vars)
colnames(tnpo3_ld) <- names(tnpo3_vars)
tnpo3_ld_df <- as.data.frame(tnpo3_ld)
tnpo3_ld_df <- cbind(variant_id = unname(tnpo3_vars), tnpo3_ld_df)
rownames(tnpo3_ld_df) <- NULL

cat("\nTNPO3 haplotype LD structure\n")
print(round(tnpo3_ld_df[, -1], 3))
write.csv(tnpo3_ld_df, "results/ld_tnpo3_haplotype.csv", row.names = FALSE)