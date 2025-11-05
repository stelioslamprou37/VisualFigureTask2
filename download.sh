#!/bin/bash

# ============================================================================
# Download script for GSE48276 data
# ============================================================================

echo "============================================================================"
echo "Downloading GSE48276 data from GEO"
echo "============================================================================"

mkdir -p input_data
cd input_data

echo ""
echo "Step 1: Fetching GSE48276 series matrix file..."
wget -q --show-progress "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE48nnn/GSE48276/matrix/GSE48276_series_matrix.txt.gz"

if [ -f "GSE48276_series_matrix.txt.gz" ]; then
    echo "✓ Successfully downloaded GSE48276_series_matrix.txt.gz"

    SIZE=$(ls -lh GSE48276_series_matrix.txt.gz | awk '{print $5}')
    echo "  File size: $SIZE"

    echo ""
    echo "Step 2: Extracting compressed file..."
    gunzip -f GSE48276_series_matrix.txt.gz

    if [ -f "GSE48276_series_matrix.txt" ]; then
        echo "✓ Successfully extracted GSE48276_series_matrix.txt"
        SIZE=$(ls -lh GSE48276_series_matrix.txt | awk '{print $5}')
        echo "  File size: $SIZE"
    else
        echo "✗ Error: Failed to extract series matrix file"
        exit 1
    fi
else
    echo "✗ Error: Failed to download GSE48276 data"
    echo ""
    echo "Manual download: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE48276"
    exit 1
fi

echo ""
echo "============================================================================"
echo "Data download complete!"
echo "============================================================================"
ls -lh
echo "============================================================================"
