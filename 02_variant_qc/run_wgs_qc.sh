#!/bin/bash
# WGS Variant QC: compile and submit bgens_qc.wdl on UK Biobank RAP
# Requires: sle_pqc.phe in /02.Phenotype_SampleQC/, bgens_qc_input.json in /03.Variant_QC/
# Output: final_WGS_snps_GRCh38_qc_pass.snplist in /03.Variant_QC/

set -e

echo "WGS Variant QC Workflow"

# Helper: check if a command is available
command_exists() { command -v "$1" >/dev/null 2>&1; }

# Install missing dependencies
echo "Checking dependencies..."

if ! command_exists python3; then
    sudo apt update && sudo apt install -y python3
fi
echo "Python3: OK"

if ! command_exists pip3; then
    sudo apt update && sudo apt install -y python3-pip
fi
echo "pip3: OK"

if ! command_exists dx; then
    pip3 install --user dxpy --break-system-packages
    export PATH="$HOME/.local/bin:$PATH"
    if ! command_exists dx; then
        echo "Error: DNAnexus CLI installation failed"; exit 1
    fi
fi
echo "DNAnexus CLI: OK"

# Java 17 required by dxCompiler 2.13.0
if ! command_exists java; then
    sudo apt update && sudo apt install -y openjdk-17-jre
fi
echo "Java: OK"

if ! command_exists wget; then
    sudo apt update && sudo apt install -y wget
fi
echo "wget: OK"

# DNAnexus login and project setup
if ! dx whoami >/dev/null 2>&1; then
    echo "DNAnexus login required:"
    dx login
fi
echo "Logged in as: $(dx whoami)"

read -p "Enter UK Biobank project ID: " PROJECT_ID
if [[ ! $PROJECT_ID =~ ^project- ]]; then
    echo "Error: Invalid project ID format (must start with 'project-')"; exit 1
fi
echo "Project: $PROJECT_ID"

# Download, compile, and submit WDL workflow
DXCOMPILER_VERSION="2.13.0"
DXCOMPILER_FILE="dxCompiler-${DXCOMPILER_VERSION}.jar"

echo "[1/5] Downloading dxCompiler ${DXCOMPILER_VERSION}..."
if [ ! -f "${DXCOMPILER_FILE}" ]; then
    wget -q "https://github.com/dnanexus/dxCompiler/releases/download/${DXCOMPILER_VERSION}/${DXCOMPILER_FILE}"
fi

echo "[2/5] Downloading WDL workflow..."
if [ ! -f "bgens_qc.wdl" ]; then
    wget -q https://raw.githubusercontent.com/dnanexus/UKB_RAP/main/end_to_end_gwas_phewas/bgens_qc/bgens_qc.wdl
fi

echo "[3/5] Downloading input configuration..."
dx download /03.Variant_QC/bgens_qc_input.json --overwrite

echo "[4/5] Compiling WDL workflow..."
java -jar ${DXCOMPILER_FILE} compile bgens_qc.wdl \
  -project ${PROJECT_ID} \
  -inputs bgens_qc_input.json \
  -archive \
  -folder /03.Variant_QC/

echo "[5/5] Submitting workflow..."
ANALYSIS_ID=$(dx run /03.Variant_QC/bgens_qc -f bgens_qc_input.dx.json --brief)

echo "Submitted: ${ANALYSIS_ID}"
echo "Monitor:   dx describe ${ANALYSIS_ID}"