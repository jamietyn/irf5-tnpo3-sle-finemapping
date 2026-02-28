# variant_definitions.R
# Shared signal and literature variant definitions
# Sourced by: 04_susie_finemapping, 05_ld_analysis, 06_colocalization, 07_functional_annotation

# Conditional signals from stepwise Firth logistic regression
# Obtained from 03_regional_association/results/top_variants.txt
signals <- c(
  signal1 = "chr7:129013434:C:CAAAA",   # lead signal, TNPO3 intron
  signal2 = "chr7:128984718:T:C",       # independent signal, TNPO3 intron
  signal3 = "chr7:128930474:CA:C"       # independent signal, IRF5 upstream
)

# Key literature variants from prior SLE fine-mapping and functional studies
# Ordered by strength of functional/genetic evidence (rank 1 = strongest)
lit_var <- c(
  rs142738614 = "chr7:128937860:C:CGCGGG", # rank 1: IRF5 promoter CGGGG indel; highest Bayesian PP (0.71); SP1 EMSA; reporter assay (Sigurdsson 2008, Graham 2007)
  rs10488631  = "chr7:128954129:T:C",       # rank 2: TNPO3 haplotype anchor; Bayesian PP=1.0 (Sigurdsson 2008); Bentham 2015 P=9e-110; eQTL effect
  rs4728142   = "chr7:128933913:G:A",       # rank 3: IRF5 promoter; CRISPR triple validation; ZBTB3 allele-specific binding; PheWAS lead (Hou 2022)
  rs13239597  = "chr7:129055929:C:A",       # rank 4: TNPO3 promoter; CRISPR + luciferase + Hi-C + 3C; EVI1 binding; independent eQTL (Thynn 2020)
  rs12534421  = "chr7:128984019:C:A",       # rank 5: TNPO3 haplotype; trans-ancestral Bayesian fine-mapping; five-ancestry replication OR 1.66-2.44 (Kottyan 2014)
  rs2004640   = "chr7:128938247:T:G",       # rank 6: IRF5 exon 1B splice site; alt TSS discovery variant; independent conditional signal (Graham 2006)
  rs10954213  = "chr7:128949373:G:A",       # rank 7: IRF5 3'UTR polyadenylation site; conditional independence; canonical H2 haplotype (Graham 2007)
  rs3757387   = "chr7:128936032:T:C",       # rank 8: IRF5 promoter haplotype proxy; eQTL Catalogue coloc top variant (current study) (Thynn 2020)
  rs35000415  = "chr7:128945562:C:T"        # rank 9: IRF5 intron; Bentham 2015 P=1e-60, Langefeld 2018 P=1e-99; in SuSiE credible set; r2=0.109 with rs3757387
)
