# Fine-Mapping the IRF5-TNPO3 Locus in Systemic Lupus Erythematosus
NUS undergraduate final year thesis project. Statistical fine-mapping of the IRF5-TNPO3 locus using UK Biobank whole-genome sequencing data to identify candidate causal variants for systemic lupus erythematosus (SLE).

## Dataset
UK Biobank DRAGEN WGS 500K release (GRCh38). Analysis window: chr7:128,930,429–129,058,173.

## Pipeline
| Stage | Description | Tools |
|---|---|---|
| 01 | Sample and phenotype QC | Python, Bash |
| 02 | Genome-wide variant QC | PLINK2 |
| 03 | Regional association and conditional analysis | PLINK2 |
| 04 | Bayesian fine-mapping | SuSiE (R) |
| 05 | LD analysis | PLINK2, R |
| 06 | Colocalization with eQTL datasets | coloc (R), eQTL Catalogue, GTEx |
| 07 | Functional annotation | VEP, HaploReg, RegulomeDB, UCSC API (Python) |
| 08 | Visualisations | ggplot2 (R) |

## Repository Structure
Scripts only. Results, large reference files, and raw genotype data are not tracked.

## Requirements
R (≥4.0): susieR, coloc, ggplot2, patchwork  
Python (≥3.8): requests  
PLINK2, bcftools, samtools