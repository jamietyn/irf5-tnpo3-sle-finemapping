# Colocalization analysis: SLE GWAS vs IRF5/TNPO3 eQTLs in immune cells
# Source: eQTL Catalogue v2 REST API (11 datasets)

library(coloc)
library(dplyr)
library(httr)
library(jsonlite)

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

# 11 immune cell datasets from eQTL Catalogue
# GTEx whole blood (QTD000356) handled separately in coloc_gtex.R
immune_datasets <- c(
  QTD000021 = "BLUEPRINT monocyte",
  QTD000026 = "BLUEPRINT neutrophil",
  QTD000379 = "Nedelec_2016 macrophage naive",
  QTD000384 = "Nedelec_2016 macrophage Listeria",
  QTD000389 = "Nedelec_2016 macrophage Salmonella",
  QTD000409 = "Quach_2016 monocyte naive",
  QTD000414 = "Quach_2016 monocyte LPS",
  QTD000419 = "Quach_2016 monocyte Pam3CSK4",
  QTD000424 = "Quach_2016 monocyte R848",
  QTD000429 = "Quach_2016 monocyte IAV",
  QTD000474 = "Schmiedel_2018 B cell naive"
)

irf5_gene_id  <- "ENSG00000128604"
tnpo3_gene_id <- "ENSG00000064419"

cat("Datasets to analyze:", length(immune_datasets), "\n")
cat("Genes: IRF5, TNPO3\n")

# Query eQTL Catalogue API with pagination, filtered to analysis window
fetch_eqtl_window <- function(dataset_id, gene_id,
                              window_start = 128930429,
                              window_end   = 129058173) {
  all_data  <- list()
  batch_size <- 1000
  start_pos  <- 0

  while (TRUE) {
    url      <- paste0("https://www.ebi.ac.uk/eqtl/api/v2/datasets/", dataset_id, "/associations")
    params   <- list(size = batch_size, start = start_pos, gene_id = gene_id)
    response <- GET(url, query = params, accept_json())

    if (status_code(response) != 200) break

    data <- fromJSON(content(response, "text", encoding = "UTF-8"))
    if (is.null(data) || nrow(data) == 0) break

    all_data[[length(all_data) + 1]] <- data
    if (nrow(data) < batch_size) break
    start_pos <- start_pos + batch_size
    Sys.sleep(0.2)
  }

  if (length(all_data) == 0) return(NULL)

  result <- bind_rows(all_data)
  result <- result %>% filter(position >= window_start & position <= window_end)
  return(result)
}

# Screen all datasets for IRF5 and TNPO3 eQTLs in the analysis window
cat("\nScreening for IRF5 eQTLs\n")
irf5_available <- list()

for (dataset_id in names(immune_datasets)) {
  cat("\nFetching", dataset_id, "-", immune_datasets[dataset_id], "...\n")
  eqtl_data <- fetch_eqtl_window(dataset_id, irf5_gene_id)

  if (!is.null(eqtl_data) && nrow(eqtl_data) > 0) {
    irf5_available[[dataset_id]] <- eqtl_data
    cat("  Found", nrow(eqtl_data), "IRF5 eQTL variants\n")
    cat("  Position range:", range(eqtl_data$position), "\n")
  } else {
    cat("  No IRF5 eQTLs in window\n")
  }
  Sys.sleep(0.5)
}

cat("Datasets with IRF5 eQTLs:", length(irf5_available), "\n")

irf5_screening <- data.frame()
for (ds in names(irf5_available)) {
  row <- data.frame(
    dataset_id   = ds,
    dataset_name = immune_datasets[ds],
    n_variants   = nrow(irf5_available[[ds]]),
    pos_min      = min(irf5_available[[ds]]$position),
    pos_max      = max(irf5_available[[ds]]$position)
  )
  irf5_screening <- rbind(irf5_screening, row)
}
write.csv(irf5_screening, "results/eqtl_screened_var_irf5.csv", row.names = FALSE)

cat("\nScreening for TNPO3 eQTLs\n")
tnpo3_available <- list()

for (dataset_id in names(immune_datasets)) {
  cat("\nFetching", dataset_id, "-", immune_datasets[dataset_id], "...\n")
  eqtl_data <- fetch_eqtl_window(dataset_id, tnpo3_gene_id)

  if (!is.null(eqtl_data) && nrow(eqtl_data) > 0) {
    tnpo3_available[[dataset_id]] <- eqtl_data
    cat("  Found", nrow(eqtl_data), "TNPO3 eQTL variants\n")
    cat("  Position range:", range(eqtl_data$position), "\n")
  } else {
    cat("  No TNPO3 eQTLs in window\n")
  }
  Sys.sleep(0.5)
}

cat("Datasets with TNPO3 eQTLs:", length(tnpo3_available), "\n")

tnpo3_screening <- data.frame()
for (ds in names(tnpo3_available)) {
  row <- data.frame(
    dataset_id   = ds,
    dataset_name = immune_datasets[ds],
    n_variants   = nrow(tnpo3_available[[ds]]),
    pos_min      = min(tnpo3_available[[ds]]$position),
    pos_max      = max(tnpo3_available[[ds]]$position)
  )
  tnpo3_screening <- rbind(tnpo3_screening, row)
}
write.csv(tnpo3_screening, "results/eqtl_screened_var_tnpo3.csv", row.names = FALSE)

# Run coloc for one dataset-gene pair
# Matches on variant ID (chr:pos:ref:alt) to avoid multiallelic position collisions
# sdY=1: eQTL Catalogue normalises expression to unit variance
run_coloc_analysis <- function(gwas_data, eqtl_data, n_samples) {

  # Construct variant ID from eQTL API fields
  eqtl_data <- eqtl_data %>%
    mutate(variant_id = paste0("chr", chromosome, ":", position, ":", ref, ":", alt))

  # Match on variant ID
  merged <- inner_join(
    gwas_data %>% select(ID, position, beta, varbeta, A1_FREQ),
    eqtl_data %>% select(variant_id, beta, se, maf),
    by     = c("ID" = "variant_id"),
    suffix = c("_gwas", "_eqtl")
  )

  if (nrow(merged) < 10) {
    cat("  Insufficient overlap (", nrow(merged), "variants)\n")
    return(NULL)
  }

  merged$snp <- merged$ID

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
    N       = n_samples,
    type    = "quant",
    sdY     = 1,
    MAF     = merged$maf,
    snp     = merged$snp
  )

  result             <- coloc.abf(gwas_list, eqtl_list)
  result$snp_results <- result$results %>% mutate(position = merged$position)
  return(result)
}

# Retrieve dataset sample size from eQTL Catalogue API
get_sample_size <- function(dataset_id) {
  url      <- paste0("https://www.ebi.ac.uk/eqtl/api/v2/datasets/", dataset_id)
  response <- GET(url, accept_json())
  if (status_code(response) == 200) {
    data <- fromJSON(content(response, "text", encoding = "UTF-8"))
    return(data$sample_size)
  }
  return(NA)
}

cat("\nRunning IRF5 colocalization\n")
irf5_coloc_results <- list()

for (dataset_id in names(irf5_available)) {
  cat("\nAnalyzing", dataset_id, "-", immune_datasets[dataset_id], "\n")
  n_samples <- get_sample_size(dataset_id)

  result <- run_coloc_analysis(
    gwas_data = assoc_coloc,
    eqtl_data = irf5_available[[dataset_id]],
    n_samples = n_samples
  )

  if (!is.null(result)) {
    irf5_coloc_results[[dataset_id]] <- result
    cat("  PP.H4:", round(result$summary["PP.H4.abf"], 4), "\n")
    cat("  Overlapping variants:", nrow(result$snp_results), "\n")
  }
  Sys.sleep(0.5)
}

cat("\nRunning TNPO3 colocalization\n")
tnpo3_coloc_results <- list()

for (dataset_id in names(tnpo3_available)) {
  cat("\nAnalyzing", dataset_id, "-", immune_datasets[dataset_id], "\n")
  n_samples <- get_sample_size(dataset_id)

  result <- run_coloc_analysis(
    gwas_data = assoc_coloc,
    eqtl_data = tnpo3_available[[dataset_id]],
    n_samples = n_samples
  )

  if (!is.null(result)) {
    tnpo3_coloc_results[[dataset_id]] <- result
    cat("  PP.H4:", round(result$summary["PP.H4.abf"], 4), "\n")
    cat("  Overlapping variants:", nrow(result$snp_results), "\n")
  }
  Sys.sleep(0.5)
}

cat("IRF5 results:", length(irf5_coloc_results), "\n")
cat("TNPO3 results:", length(tnpo3_coloc_results), "\n")

# Create and save coloc summary tables
irf5_coloc_summary <- data.frame()
for (dataset_id in names(irf5_coloc_results)) {
  row <- data.frame(
    dataset_id   = dataset_id,
    dataset_name = immune_datasets[dataset_id],
    n_variants   = nrow(irf5_coloc_results[[dataset_id]]$snp_results),
    PP.H0        = irf5_coloc_results[[dataset_id]]$summary["PP.H0.abf"],
    PP.H1        = irf5_coloc_results[[dataset_id]]$summary["PP.H1.abf"],
    PP.H2        = irf5_coloc_results[[dataset_id]]$summary["PP.H2.abf"],
    PP.H3        = irf5_coloc_results[[dataset_id]]$summary["PP.H3.abf"],
    PP.H4        = irf5_coloc_results[[dataset_id]]$summary["PP.H4.abf"]
  )
  irf5_coloc_summary <- rbind(irf5_coloc_summary, row)
}
print(irf5_coloc_summary)
write.csv(irf5_coloc_summary, "results/coloc_summary_irf5.csv", row.names = FALSE)

if (length(tnpo3_coloc_results) > 0) {
  tnpo3_coloc_summary <- data.frame()
  for (dataset_id in names(tnpo3_coloc_results)) {
    row <- data.frame(
      dataset_id   = dataset_id,
      dataset_name = immune_datasets[dataset_id],
      n_variants   = nrow(tnpo3_coloc_results[[dataset_id]]$snp_results),
      PP.H0        = tnpo3_coloc_results[[dataset_id]]$summary["PP.H0.abf"],
      PP.H1        = tnpo3_coloc_results[[dataset_id]]$summary["PP.H1.abf"],
      PP.H2        = tnpo3_coloc_results[[dataset_id]]$summary["PP.H2.abf"],
      PP.H3        = tnpo3_coloc_results[[dataset_id]]$summary["PP.H3.abf"],
      PP.H4        = tnpo3_coloc_results[[dataset_id]]$summary["PP.H4.abf"]
    )
    tnpo3_coloc_summary <- rbind(tnpo3_coloc_summary, row)
  }
  print(tnpo3_coloc_summary)
  write.csv(tnpo3_coloc_summary, "results/coloc_summary_tnpo3.csv", row.names = FALSE)
}

# Extract top colocalised variant per dataset
irf5_top_variants <- data.frame()
for (dataset_id in names(irf5_coloc_results)) {
  top_var <- irf5_coloc_results[[dataset_id]]$snp_results %>%
    arrange(desc(SNP.PP.H4)) %>%
    slice(1)
  row <- data.frame(
    dataset_id   = dataset_id,
    dataset_name = immune_datasets[dataset_id],
    position     = top_var$position,
    SNP.PP.H4    = top_var$SNP.PP.H4
  )
  irf5_top_variants <- rbind(irf5_top_variants, row)
  cat(dataset_id, "top variant:", top_var$position, "PP.H4 =", round(top_var$SNP.PP.H4, 4), "\n")
}
write.csv(irf5_top_variants, "results/coloc_top_variants_eqtlcat_irf5.csv", row.names = FALSE)

tnpo3_top_variants <- data.frame()
for (dataset_id in names(tnpo3_coloc_results)) {
  top_var <- tnpo3_coloc_results[[dataset_id]]$snp_results %>%
    arrange(desc(SNP.PP.H4)) %>%
    slice(1)
  row <- data.frame(
    dataset_id   = dataset_id,
    dataset_name = immune_datasets[dataset_id],
    position     = top_var$position,
    SNP.PP.H4    = top_var$SNP.PP.H4
  )
  tnpo3_top_variants <- rbind(tnpo3_top_variants, row)
  cat(dataset_id, "top variant:", top_var$position, "PP.H4 =", round(top_var$SNP.PP.H4, 4), "\n")
}
write.csv(tnpo3_top_variants, "results/coloc_top_variants_eqtlcat_tnpo3.csv", row.names = FALSE)

# Save all SNP-level posteriors stacked across datasets
irf5_snp_posteriors <- data.frame()
for (dataset_id in names(irf5_coloc_results)) {
  snp_df              <- irf5_coloc_results[[dataset_id]]$snp_results
  snp_df$dataset_id   <- dataset_id
  snp_df$dataset_name <- immune_datasets[dataset_id]
  irf5_snp_posteriors <- rbind(irf5_snp_posteriors, snp_df)
}
write.csv(irf5_snp_posteriors, "results/coloc_snp_posteriors_irf5.csv", row.names = FALSE)

if (length(tnpo3_coloc_results) > 0) {
  tnpo3_snp_posteriors <- data.frame()
  for (dataset_id in names(tnpo3_coloc_results)) {
    snp_df               <- tnpo3_coloc_results[[dataset_id]]$snp_results
    snp_df$dataset_id    <- dataset_id
    snp_df$dataset_name  <- immune_datasets[dataset_id]
    tnpo3_snp_posteriors <- rbind(tnpo3_snp_posteriors, snp_df)
  }
  write.csv(tnpo3_snp_posteriors, "results/coloc_snp_posteriors_tnpo3.csv", row.names = FALSE)
}