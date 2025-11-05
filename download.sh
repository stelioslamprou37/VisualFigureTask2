#!/bin/bash

# ============================================================================
# Download Script: GSE48276 Data + Supplementary Survival Information
# ============================================================================

echo "============================================================================"
echo "Downloading GSE48276 Data + Supplementary Survival Information"
echo "============================================================================"
echo ""

mkdir -p input_data
cd input_data

# ============================================================================
# Step 1: Download expression data from GEO
# ============================================================================

echo "Step 1: Downloading GSE48276 series matrix from GEO..."
echo ""

wget -q --show-progress \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE48nnn/GSE48276/matrix/GSE48276_series_matrix.txt.gz"

if [ -f "GSE48276_series_matrix.txt.gz" ]; then
    echo "✓ Downloaded GSE48276_series_matrix.txt.gz"
    SIZE=$(ls -lh GSE48276_series_matrix.txt.gz | awk '{print $5}')
    echo "  File size: $SIZE"
    echo ""

    echo "Extracting..."
    gunzip -f GSE48276_series_matrix.txt.gz

    if [ -f "GSE48276_series_matrix.txt" ]; then
        echo "✓ Extracted GSE48276_series_matrix.txt"
        SIZE=$(ls -lh GSE48276_series_matrix.txt | awk '{print $5}')
        echo "  File size: $SIZE"
        echo ""
    else
        echo "✗ Error: Failed to extract"
        exit 1
    fi
else
    echo "✗ Error: Failed to download"
    exit 1
fi

# ============================================================================
# Step 2: Instructions for supplementary survival data
# ============================================================================

echo "============================================================================"
echo "Step 2: Download Supplementary Survival Data"
echo "============================================================================"
echo ""
echo "The survival data is embedded in the R script as 71 samples with outcomes"
echo ""
echo "Source: Supplementary File S1"
echo "Paper: Frontiers in Immunology 2025"
echo "DOI: 10.3389/fimmu.2024.1430583"
echo "URL: https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2024.1430583/full"
echo ""
echo "Manual download (optional reference):"
echo "  1. Go to paper URL above"
echo "  2. Look for 'Supplementary Materials' section"
echo "  3. Download 'Supplementary File S1' (Clinical data)"
echo "  4. Contains survival times (futime_OS_days) and event status for 71 samples"
echo ""

# ============================================================================
# Step 3: Data verification
# ============================================================================

echo "============================================================================"
echo "Data Files Summary"
echo "============================================================================"
echo ""
ls -lh
echo ""
echo "============================================================================"
echo "Ready to run analysis!"
echo "============================================================================"
echo ""
echo "For Docker:"
echo "  docker run -v \$(pwd):/workspace bladder-signature expert_solution.R"
echo ""
echo "For local R:"
echo "  Rscript ../expert_solution.R"
echo ""
