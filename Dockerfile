# ============================================================================
# Dockerfile: 18-Gene Glycolysis Signature Analysis
# Base: rocker/r-ver:4.3.0
# Purpose: Reproducible environment for bladder cancer ROC analysis
# ============================================================================

FROM rocker/r-ver:4.3.0

WORKDIR /workspace

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    wget \
    gzip \
    && rm -rf /var/lib/apt/lists/*

# Install BiocManager
RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org/')"

# Install Bioconductor packages
RUN R -e "BiocManager::install(c( \
    'GEOquery', \
    'preprocessCore' \
    ), ask = FALSE, update = FALSE)"

# Install CRAN packages
RUN R -e "install.packages(c( \
    'survival', \
    'survminer', \
    'timeROC', \
    'ggplot2' \
    ), repos='https://cloud.r-project.org/')"

# Create working directories
RUN mkdir -p /workspace/input_data /workspace/output

ENV DEBIAN_FRONTEND=noninteractive

# Set entrypoint for running R scripts
ENTRYPOINT ["/usr/local/bin/Rscript"]
CMD ["expert_solution.R"]

# ============================================================================
# Build and run:
#   docker build -t bladder-signature .
#   docker run -v $(pwd):/workspace bladder-signature expert_solution.R
# ============================================================================
