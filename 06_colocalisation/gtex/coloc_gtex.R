# Colocalization analysis: SLE GWAS vs IRF5 eQTLs in GTEx v10 whole blood
# GTEx Portal API v2; query by genomic location (chr7:128,930,429-129,058,173)
# Handled separately from coloc_eqtl_catalogue.R: the GTEx Portal API queries by location
# rather than gene ID, returning all eGenes in the window and enabling regional breakdown.
# Standard error is not provided directly and must be approximated from normalised effect sizes and p-value.

library(coloc)
library(dplyr)
library(jsonlite)
library(httr)

# Load GWAS association results (unconditional, all 444 variants)
assoc <- read.table("../../03_regional_association/results/conditional/irf5_tnpo3_round0_sorted.txt",
                    header = TRUE, stringsAsFactors = FALSE, comment.char = "")

assoc_coloc <- assoc %>%
  mutate(
    beta     = log(OR),
    varbeta  = LOG.OR._SE^2,
    position = POS
  ) %>%
  select(ID, position, beta, varbeta, A1_FREQ, P)

cat("GWAS variants loaded:", nrow(assoc_coloc), "\n")
cat("Position range:", range(assoc_coloc$position), "\n")

# Query GTEx Portal API v2 by genomic location
# Returns all eGenes tested in this window; IRF5 is the primary eGene in GTEx whole blood
gtex_url    <- "https://gtexportal.org/api/v2/association/singleTissueEqtlByLocation"
gtex_params <- list(
  chromosome         = "chr7",
  start              = 128930429,
  end                = 129058173,
  tissueSiteDetailId = "Whole_Blood"
)

cat("\nQuerying GTEx portal for IRF5 eQTLs...\n")
gtex_response <- GET(gtex_url, query = gtex_params, accept_json())

if (status_code(gtex_response) == 200) {
  gtex_data_raw <- fromJSON(content(gtex_response, "text", encoding = "UTF-8"))
  gtex_irf5     <- gtex_data_raw$singleTissueEqtl
  cat("Total GTEx IRF5 eQTL variants:", nrow(gtex_irf5), "\n")
  cat("Position range:", range(gtex_irf5$pos), "\n")
} else {
  cat("GTEx API query failed with status:", status_code(gtex_response), "\n")
  stop("Cannot proceed without GTEx data")
}

# Prepare GTEx data for coloc
# SE approximated as |NES/z| where z = sign(NES) * |qnorm(p/2)| — standard approximation
# Variant ID constructed from GTEx variantId field (format: chr7_128933913_G_A_b38 → chr7:128933913:G:A)
gtex_irf5_coloc <- gtex_irf5 %>%
  mutate(
    z          = sign(nes) * abs(qnorm(pValue / 2)),
    se         = abs(nes / z),
    position   = pos,
    variant_id = gsub("_", ":", gsub("_b38$", "", variantId))
  ) %>%
  select(variant_id, position, beta = nes, se, pValue)

cat("\nGTEx data prepared for coloc\n")
cat("Unique variants:", nrow(gtex_irf5_coloc), "\n")
cat("Position range:", range(gtex_irf5_coloc$position), "\n")

# Match on variant ID (chr:pos:ref:alt) to avoid multiallelic position collisions
merged <- inner_join(
  assoc_coloc %>% select(ID, position, beta, varbeta, A1_FREQ),
  gtex_irf5_coloc %>% select(variant_id, beta, se),
  by     = c("ID" = "variant_id"),
  suffix = c("_gwas", "_eqtl")
)

cat("\nOverlapping variants:", nrow(merged), "\n")
cat("Position range:", range(merged$position), "\n")

if (nrow(merged) < 10) stop("Insufficient overlap for coloc analysis")

merged$snp <- merged$ID

# N=800: GTEx v10 whole blood with genotype (verified at gtexportal.org/home/tissue/Whole_Blood)
# sdY=1: GTEx expression data is standardised to unit variance
gwas_list <- list(
  beta    = merged$beta_gwas,
  varbeta = merged$varbeta,
  N       = 405448,
  type    = "cc",
  s       = 565/405448,
  MAF     = merged$A1_FREQ,
  snp     = merged$snp
)

eqtl_list <- list(
  beta    = merged$beta_eqtl,
  varbeta = merged$se^2,
  N       = 800,
  type    = "quant",
  sdY     = 1,
  snp     = merged$snp
)

cat("\nRunning colocalization...\n")
gtex_coloc_result <- coloc.abf(gwas_list, eqtl_list)

cat("\nColocalization results:\n")
print(gtex_coloc_result$summary)

# Save summary
gtex_summary <- data.frame(
  dataset       = "GTEx_whole_blood",
  gene          = "IRF5",
  n_overlapping = nrow(merged),
  PP.H0         = gtex_coloc_result$summary["PP.H0.abf"],
  PP.H1         = gtex_coloc_result$summary["PP.H1.abf"],
  PP.H2         = gtex_coloc_result$summary["PP.H2.abf"],
  PP.H3         = gtex_coloc_result$summary["PP.H3.abf"],
  PP.H4         = gtex_coloc_result$summary["PP.H4.abf"]
)
print(gtex_summary)
write.csv(gtex_summary, "results/coloc_summary_gtex_irf5.csv", row.names = FALSE)

# SNP-level posteriors (all overlapping variants)
gtex_snp_results <- gtex_coloc_result$results %>%
  mutate(position = merged$position) %>%
  arrange(desc(SNP.PP.H4))
write.csv(gtex_snp_results, "results/coloc_snp_posteriors_gtex_irf5.csv", row.names = FALSE)

# Top colocalised variant
top_variant <- gtex_snp_results %>% slice(1)
cat("\nTop colocalized variant:\n")
cat("Position:", top_variant$position, "\n")
cat("SNP.PP.H4:", round(top_variant$SNP.PP.H4, 4), "\n")

gtex_top_variant <- data.frame(
  dataset   = "GTEx_whole_blood",
  gene      = "IRF5",
  position  = top_variant$position,
  SNP.PP.H4 = top_variant$SNP.PP.H4
)
write.csv(gtex_top_variant, "results/coloc_top_variants_gtex_irf5.csv", row.names = FALSE)

# SNP.PP.H4 for the 3 conditional signals
# Note: signals may be absent from GTEx if not covered by Whole_Blood eQTL data
key_signals <- c(signal1 = 129013434, signal2 = 128984718, signal3 = 128930474)

signal_posteriors <- data.frame()
for (sig in names(key_signals)) {
  pos      <- key_signals[sig]
  in_gtex  <- pos %in% gtex_irf5$pos
  in_coloc <- pos %in% merged$position
  pp_h4    <- NA
  if (in_coloc) {
    match_row <- gtex_snp_results[gtex_snp_results$position == pos, ]
    if (nrow(match_row) > 0) pp_h4 <- match_row$SNP.PP.H4[1]
  }
  signal_posteriors <- rbind(signal_posteriors, data.frame(
    signal    = sig,
    position  = pos,
    in_gtex   = in_gtex,
    in_coloc  = in_coloc,
    SNP.PP.H4 = pp_h4
  ))
}
print(signal_posteriors, row.names = FALSE)
write.csv(signal_posteriors, "results/signal_posteriors_gtex.csv", row.names = FALSE)

# Regional breakdown: where is the colocalisation posterior mass - IRF5 or TNPO3 region?
# IRF5 gene body: chr7:128,930,179-128,955,260 (GRCh38, UCSC)
# TNPO3 gene body: chr7:128,954,180-129,058,173 (GRCh38, UCSC)
irf5_end    <- 128955260
tnpo3_start <- 128954180

irf5_stats  <- gtex_snp_results %>%
  filter(position < irf5_end) %>%
  summarise(n_variants = n(), max_PP.H4 = max(SNP.PP.H4), sum_PP.H4 = sum(SNP.PP.H4))

tnpo3_stats <- gtex_snp_results %>%
  filter(position >= tnpo3_start) %>%
  summarise(n_variants = n(), max_PP.H4 = max(SNP.PP.H4), sum_PP.H4 = sum(SNP.PP.H4))

cat("IRF5 region (<128,955,260):\n")
cat("  Variants:", irf5_stats$n_variants,
    " Max SNP.PP.H4:", round(irf5_stats$max_PP.H4, 4),
    " Sum:", round(irf5_stats$sum_PP.H4, 4), "\n")

cat("TNPO3 region (>=128,954,180):\n")
cat("  Variants:", tnpo3_stats$n_variants,
    " Max SNP.PP.H4:", round(tnpo3_stats$max_PP.H4, 4),
    " Sum:", round(tnpo3_stats$sum_PP.H4, 4), "\n")

regional_summary <- data.frame(
  region         = c("IRF5", "TNPO3"),
  position_range = c("<128955260", ">=128954180"),
  n_variants     = c(irf5_stats$n_variants,  tnpo3_stats$n_variants),
  max_PP.H4      = c(irf5_stats$max_PP.H4,   tnpo3_stats$max_PP.H4),
  sum_PP.H4      = c(irf5_stats$sum_PP.H4,   tnpo3_stats$sum_PP.H4)
)
write.csv(regional_summary, "results/regional_breakdown_gtex.csv", row.names = FALSE)