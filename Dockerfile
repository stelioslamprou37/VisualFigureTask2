# ============================================================================
# Dockerfile for 18-Gene Glycolysis Signature Validation
# Base image: rocker/r-ver with R 4.3.0
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

# Create directories
RUN mkdir -p /workspace/input_data

ENV DEBIAN_FRONTEND=noninteractive

CMD ["/bin/bash"]

# ============================================================================
# Build: docker build -t bladder-signature .
# Run: docker run -v $(pwd):/workspace -it bladder-signature Rscript expert_solution.R
# ============================================================================
