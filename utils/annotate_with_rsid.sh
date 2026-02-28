#!/bin/bash
set -euo pipefail

# Annotate variant files (CSV or VCF) with rsIDs using bcftools
# Merges rsIDs back into the original input file as an additional column.
# NOTE: Specific to chromosome 7. Expects chr7 FASTA and pre-normalised dbSNP VCF.
#
# Usage:
#   bash annotate_with_rsid.sh <dbsnp_chr7_norm.vcf.gz> <chr7.fa> <input_file> <output_prefix>
#
# Arguments:
#   dbsnp_chr7_norm.vcf.gz  — Pre-normalised dbSNP VCF (already indexed with tabix)
#   chr7.fa                 — chr7 reference FASTA
#   input_file              — Input CSV (with header) or VCF; CSV must have ID column in
#                             format chr7:POS:REF:ALT
#   output_prefix           — Prefix for all intermediate and output files
#
# Output:
#   <output_prefix>_rsid_annotated.csv  — Original input CSV with rsID column appended
#                                         (rsID is "." if no match found in dbSNP)

DBSNP=$1
FASTA=$2
INPUT=$3
PREFIX=$4

SCRIPT_DIR="$(dirname "$0")"
VCF_RAW="${PREFIX}_raw.vcf"
VCF_NORM="${PREFIX}_norm.vcf.gz"
VCF_ANNOT="${PREFIX}_annot.vcf.gz"
RSID_MAP="${PREFIX}_rsid_map.txt"

# -------------------------
# Step 0: Normalise dbSNP if not already done
# Skipped if <dbsnp>.norm.vcf.gz already exists alongside the input dbSNP
# -------------------------
DBSNP_NORM="${DBSNP%.vcf.gz}.norm.vcf.gz"
if [ ! -f "$DBSNP_NORM" ]; then
  echo "Normalising dbSNP..."
  bcftools norm -f "$FASTA" -m -both "$DBSNP" -Oz -o "$DBSNP_NORM"
  tabix -p vcf "$DBSNP_NORM"
else
  echo "Using existing normalised dbSNP: $DBSNP_NORM"
fi
DBSNP="$DBSNP_NORM"

# -------------------------
# Step 1: Convert input to VCF
# -------------------------
if [[ "$INPUT" == *.csv ]]; then
  bash "${SCRIPT_DIR}/csv_to_vcf.sh" "$INPUT" "$VCF_RAW"
elif [[ "$INPUT" == *.vcf ]]; then
  cp "$INPUT" "$VCF_RAW"
else
  echo "Error: unsupported input format '${INPUT}'. Use .csv or .vcf."
  exit 1
fi

# -------------------------
# Step 2: Ensure VCF has valid header for bcftools norm
# VCF may use 'chr7' or '7' as contig name; inject correct contig header if missing.
# -------------------------

# Detect contig naming convention used in VCF
if grep -q '^chr7' "$VCF_RAW"; then
  CONTIG_ID="chr7"
  CONTIG_LINE="##contig=<ID=chr7,length=159345973>"
else
  CONTIG_ID="7"
  CONTIG_LINE="##contig=<ID=7,length=159345973>"
fi

VCF_FIXED="${PREFIX}_fixed.vcf"
if ! grep -q "##contig=<ID=${CONTIG_ID}" "$VCF_RAW"; then
  awk -v contig="$CONTIG_LINE" \
    '/^#CHROM/{print contig} {print}' "$VCF_RAW" > "$VCF_FIXED"
else
  cp "$VCF_RAW" "$VCF_FIXED"
fi

# -------------------------
# Step 3: Index FASTA if needed
# -------------------------
if [ ! -f "${FASTA}.fai" ]; then
  echo "Indexing FASTA..."
  samtools faidx "$FASTA"
fi

# -------------------------
# Step 4: Normalise user VCF against reference
# Left-aligns indels and splits multiallelic sites for consistent matching to dbSNP
# -------------------------
echo "Normalising variants..."
bcftools norm -f "$FASTA" -m -both "$VCF_FIXED" -Oz -o "${VCF_NORM}.tmp.vcf.gz"
bcftools sort -Oz -o "$VCF_NORM" "${VCF_NORM}.tmp.vcf.gz"
tabix -p vcf "$VCF_NORM"
rm -f "$VCF_FIXED" "${VCF_NORM}.tmp.vcf.gz"

# -------------------------
# Step 5: Annotate with dbSNP rsIDs
# Uses pre-normalised dbSNP VCF; matches by position + alleles after normalisation
# -------------------------
echo "Annotating with dbSNP rsIDs..."
bcftools annotate -a "$DBSNP" -c ID "$VCF_NORM" -Oz -o "$VCF_ANNOT"
tabix -p vcf "$VCF_ANNOT"

# -------------------------
# Step 6: Extract variant ID → rsID map
# Columns: original_id (chr:pos:ref:alt from INFO/original or reconstructed), rsID
# Output is tab-separated: CHROM:POS:REF:ALT <tab> rsID
# -------------------------
echo "Extracting rsID map..."
bcftools query -f '%CHROM:%POS:%REF:%ALT\t%ID\n' "$VCF_ANNOT" > "$RSID_MAP"

# -------------------------
# Step 7: Merge rsIDs back into original CSV
# Matches on the ID column (chr7:POS:REF:ALT) in the original CSV
# rsID is "." if no dbSNP match was found
# -------------------------
python3 - <<EOF
import csv

# Load rsID map keyed by POS:REF:ALT (chrom-independent)
# Only retain genuine rs-prefixed IDs
rsid_map = {}
with open("${RSID_MAP}") as f:
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) == 2:
            var_key, rsid = parts
            if rsid.startswith("rs"):
                fields = var_key.split(":")
                if len(fields) == 4:
                    rsid_map[":".join(fields[1:])] = rsid  # POS:REF:ALT

# Merge into CSV
rows = []
with open("${INPUT}") as fin:
    reader = csv.DictReader(fin)
    fieldnames = list(reader.fieldnames) + ["rsID"]
    for row in reader:
        key = f"{row['POS']}:{row['REF']}:{row['ALT']}"
        row["rsID"] = rsid_map.get(key, ".")
        rows.append(row)

with open("${PREFIX}_rsid_annotated.csv", "w", newline="") as fout:
    writer = csv.DictWriter(fout, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)

total      = len(rows)
annotated  = sum(1 for r in rows if r["rsID"] != ".")
unannotated = total - annotated
print(f"rsID annotation: {annotated}/{total} annotated, {unannotated}/{total} unannotated (.)")
print("Unannotated variants:")
for r in rows:
    if r["rsID"] == ".":
        print(f"  {r['ID']}")
EOF

# -------------------------
# Cleanup intermediates
# -------------------------
rm -f "$VCF_RAW" "$VCF_NORM" "${VCF_NORM}.tbi" "$VCF_ANNOT" "${VCF_ANNOT}.tbi" "$RSID_MAP"

echo "Done. Output: ${PREFIX}_rsid_annotated.csv"