#!/bin/bash

# ============================================================================
# Download Script: GSE48276 Data + Supplementary Survival Information
# ============================================================================

echo "============================================================================"
echo "Downloading GSE48276 Data + Supplementary Survival Information"
echo "============================================================================"
echo ""

# ============================================================================
# Determine working directory structure
# ============================================================================

# Check if we're already in input_data or need to create it
if [[ "${PWD##*/}" == "input_data" ]]; then
    echo "Already in input_data directory"
    INPUT_DIR="."
    PARENT_DIR=".."
else
    echo "Creating input_data directory..."
    mkdir -p input_data
    INPUT_DIR="input_data"
    PARENT_DIR="."
fi

cd "$INPUT_DIR"
echo "Working directory: $(pwd)"
echo ""

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
echo "Step 2: Supplementary Survival Data Information"
echo "============================================================================"
echo ""
echo "The survival data is embedded in the analysis script as 71 samples with outcomes"
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

# ============================================================================
# Step 4: Detect analysis script and provide run instructions
# ============================================================================

echo "============================================================================"
echo "Ready to run analysis!"
echo "============================================================================"
echo ""

# Find R scripts in parent directory
ANALYSIS_SCRIPT=""
if [ -f "$PARENT_DIR/expert_solution.R" ]; then
    ANALYSIS_SCRIPT="expert_solution.R"
elif [ -f "$PARENT_DIR/analysis.R" ]; then
    ANALYSIS_SCRIPT="analysis.R"
elif [ -f "$PARENT_DIR/solution.R" ]; then
    ANALYSIS_SCRIPT="solution.R"
else
    # Look for any R script that might be the analysis script
    R_SCRIPTS=$(find "$PARENT_DIR" -maxdepth 1 -name "*.R" -type f 2>/dev/null)
    if [ -n "$R_SCRIPTS" ]; then
        ANALYSIS_SCRIPT=$(basename $(echo "$R_SCRIPTS" | head -n1))
    fi
fi

if [ -n "$ANALYSIS_SCRIPT" ]; then
    echo "Detected analysis script: $ANALYSIS_SCRIPT"
    echo ""
    echo "To run analysis:"
    echo ""
    echo "Option 1 - Docker:"
    echo "  docker run -v \$(pwd):/workspace bladder-signature $ANALYSIS_SCRIPT"
    echo ""
    echo "Option 2 - Local R (from parent directory):"
    echo "  cd .."
    echo "  Rscript $ANALYSIS_SCRIPT"
    echo ""
    echo "Option 3 - Local R (from this directory):"
    echo "  Rscript ../$ANALYSIS_SCRIPT"
    echo ""
else
    echo "No analysis script detected in parent directory."
    echo ""
    echo "Expected script names: expert_solution.R, analysis.R, solution.R"
    echo ""
    echo "To run your analysis script:"
    echo ""
    echo "Option 1 - Docker:"
    echo "  docker run -v \$(pwd):/workspace bladder-signature <your_script.R>"
    echo ""
    echo "Option 2 - Local R (from parent directory):"
    echo "  cd .."
    echo "  Rscript <your_script.R>"
    echo ""
    echo "Option 3 - Local R (from this directory):"
    echo "  Rscript ../<your_script.R>"
    echo ""
fi

echo "============================================================================"
echo "Download complete!"
echo "============================================================================"
