# Post-colocalisation LD analysis: r^2 between coloc top variants, signals, and literature variants
# Inputs: coloc top variants from eQTL Catalogue and GTEx; LD matrix from regional association analysis

source("../utils/variant_definitions.R")

# Load LD matrix
ld_r_raw    <- as.matrix(read.table("../03_regional_association/results/ld_matrices/irf5_tnpo3_ld_r.unphased.vcor1",
                                    header = FALSE))
variant_ids <- readLines("../03_regional_association/results/ld_matrices/irf5_tnpo3_ld_r.unphased.vcor1.vars")
rownames(ld_r_raw) <- variant_ids
colnames(ld_r_raw) <- variant_ids

cat("Matrix symmetric:", all(abs(ld_r_raw - t(ld_r_raw)) < 1e-10), "\n")
cat("Diagonal all 1s:",  all(abs(diag(ld_r_raw) - 1) < 1e-10), "\n")
cat("LD matrix dimensions:", dim(ld_r_raw), "\n")
ld_r <- ld_r_raw

# Check all reference variants are present in LD matrix
all_refs <- c(signals, lit_var)
missing  <- all_refs[!all_refs %in% colnames(ld_r)]
if (length(missing) > 0) {
  cat("WARNING: Missing from LD matrix:\n")
  print(missing)
}

# Load coloc top variants from eQTL Catalogue (IRF5) and GTEx
eqtl_top <- read.csv("../06_colocalisation/eqtl_catalogue/results/coloc_top_variants_eqtlcat_irf5.csv",
                     stringsAsFactors = FALSE)
gtex_top  <- read.csv("../06_colocalisation/gtex/results/coloc_top_variants_gtex_irf5.csv",
                      stringsAsFactors = FALSE)

cat("eQTL Catalogue datasets:", nrow(eqtl_top), "\n")
cat("Unique coloc positions:", length(unique(eqtl_top$position)), "\n")

# Map coloc positions to variant IDs - coloc outputs store position; LD matrix uses chr:pos:ref:alt
assoc <- read.table("../03_regional_association/results/conditional/irf5_tnpo3_round0_sorted.txt",
                    header = TRUE, stringsAsFactors = FALSE, comment.char = "")

match_position_to_id <- function(positions, assoc_data) {
  ids <- character(length(positions))
  for (i in seq_along(positions)) {
    matches <- assoc_data$ID[assoc_data$POS == positions[i]]
    ids[i]  <- if (length(matches) > 0) matches[1] else NA
  }
  ids
}

eqtl_positions  <- unique(eqtl_top$position)
eqtl_ids        <- match_position_to_id(eqtl_positions, assoc)
names(eqtl_ids) <- paste0("eqtl_", eqtl_positions)

gtex_id        <- match_position_to_id(gtex_top$position, assoc)
names(gtex_id) <- paste0("gtex_", gtex_top$position)

# Filter to variants present in LD matrix; deduplicate across sources
eqtl_ids         <- eqtl_ids[!is.na(eqtl_ids) & eqtl_ids %in% colnames(ld_r)]
gtex_id          <- gtex_id[!is.na(gtex_id)   & gtex_id  %in% colnames(ld_r)]
unique_coloc_ids <- unique(c(eqtl_ids, gtex_id))

cat("Variants available for LD analysis:", length(unique_coloc_ids), "\n")

# Pairwise r^2 between coloc top variants - checks statistical independence of colocalised variants
if (length(unique_coloc_ids) > 0) {
  coloc_ld <- ld_r[unique_coloc_ids, unique_coloc_ids]^2
  rownames(coloc_ld) <- unique_coloc_ids
  colnames(coloc_ld) <- unique_coloc_ids
  signal_labels <- sapply(unique_coloc_ids, function(id) {
    sig_match <- names(signals)[signals == id]
    if (length(sig_match) > 0) return(sig_match)
    return(id)
  })
  output_matrix <- rbind(signal_labels, round(coloc_ld, 6))
  output_df <- as.data.frame(output_matrix, stringsAsFactors = FALSE)
  output_df <- cbind(variant_id = c("", rownames(coloc_ld)), output_df)
  write.table(output_df, "results/ld_coloc_top_variants_pairwise.csv", sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)
}

# r^2 between coloc top variants and conditional signals - determines which signal(s) each colocalized variant tags
coloc_vs_signals <- ld_r[unique_coloc_ids, signals]^2
rownames(coloc_vs_signals) <- unique_coloc_ids
colnames(coloc_vs_signals) <- names(signals)
signal_labels <- names(signals)
output_matrix <- rbind(signal_labels, round(coloc_vs_signals, 6))
output_df <- as.data.frame(output_matrix, stringsAsFactors = FALSE)
output_df <- cbind(variant_id = c("", rownames(coloc_vs_signals)), output_df)
write.table(output_df, "results/ld_coloc_vs_signals.csv", sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)

# r^2 between coloc top variants and known functional variants - maps colocalised variants onto published haplotypes
coloc_vs_lit <- ld_r[unique_coloc_ids, lit_var]^2
rownames(coloc_vs_lit) <- unique_coloc_ids
colnames(coloc_vs_lit) <- names(lit_var)
rsid_labels <- names(lit_var)
output_matrix <- rbind(rsid_labels, round(coloc_vs_lit, 6))
output_df <- as.data.frame(output_matrix, stringsAsFactors = FALSE)
output_df <- cbind(variant_id = c("", rownames(coloc_vs_lit)), output_df)
write.table(output_df, "results/ld_coloc_vs_literature.csv", sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)

# r^2 between coloc top variants and credible set variants
cs_variants <- read.csv("../04_susie_finemapping/results/credible_set_variants.csv",
                        stringsAsFactors = FALSE)
cs_ids_in_ld <- cs_variants$ID

coloc_vs_cs <- ld_r[unique_coloc_ids, cs_ids_in_ld]^2
rownames(coloc_vs_cs) <- unique_coloc_ids
colnames(coloc_vs_cs) <- cs_ids_in_ld
output_matrix <- rbind(cs_ids_in_ld, round(coloc_vs_cs, 6))
output_df <- as.data.frame(output_matrix, stringsAsFactors = FALSE)
output_df <- cbind(variant_id = c("", rownames(coloc_vs_cs)), output_df)
write.table(output_df, "results/ld_coloc_vs_credible_set.csv", sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)