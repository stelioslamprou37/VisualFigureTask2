#!/usr/bin/env Rscript

# ============================================================================
# LOCAL VERSION: Clustered Heatmap of DEGs in Traumatic Brain Injury
# Figure 5 from PLOS ONE 2025 paper
# Data: GSE223245 (loaded from local file)
# ============================================================================

suppressPackageStartupMessages({
  library(GEOquery)
  library(limma)
  library(pheatmap)
  library(RColorBrewer)
})

cat("============================================================================\n")
cat("TBI Heatmap Analysis - GSE223245 (LOCAL VERSION)\n")
cat("============================================================================\n\n")

# ============================================================================
# Step 1: Load GSE223245 data from LOCAL file
# ============================================================================
cat("Step 1: Loading GSE223245 data from local file...\n")

# MODIFY THIS PATH to your local file location
# Windows example: "C:/Users/YourName/Desktop/project/input_data/GSE223245_series_matrix.txt"
# Mac/Linux example: "~/Desktop/project/input_data/GSE223245_series_matrix.txt"
data_file <- "./input_data/GSE223245_series_matrix.txt"

# Check if file exists
if(!file.exists(data_file)) {
  stop("ERROR: Cannot find ", data_file, "\n",
       "Please ensure GSE223245_series_matrix.txt is in the input_data folder.")
}

cat("  Loading:", data_file, "\n")
gset <- getGEO(filename = data_file, GSEMatrix = TRUE)

# Extract expression matrix and metadata
exprs_data <- exprs(gset)
pdata <- pData(gset)
fdata <- fData(gset)

cat("  ✓ Data loaded successfully\n")
cat("  Expression matrix dimensions:", dim(exprs_data), "\n")
cat("  Samples:", ncol(exprs_data), "\n")
cat("  Probes:", nrow(exprs_data), "\n\n")

# ============================================================================
# Step 2: Define sample groups
# ============================================================================
cat("Step 2: Defining sample groups...\n")

sample_names <- pdata$title
cat("  Sample names:\n")
print(sample_names)

# Define groups
groups <- ifelse(grepl("TBI|Patient", sample_names, ignore.case = TRUE), "TBI", "Control")

cat("\n  Group assignment:\n")
print(table(groups))
cat("\n")

# ============================================================================
# Step 3: Differential expression analysis with limma
# ============================================================================
cat("Step 3: Performing differential expression analysis...\n")

# Create design matrix
design <- model.matrix(~0 + factor(groups))
colnames(design) <- c("Control", "TBI")

# Fit linear model
fit <- lmFit(exprs_data, design)

# Make contrast
contrast_matrix <- makeContrasts(TBI_vs_Control = TBI - Control, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# Get results
results <- topTable(fit2, coef = 1, number = Inf, adjust.method = "BH")

cat("  Total genes tested:", nrow(results), "\n")

# Filter for significant DEGs
sig_degs <- results[results$adj.P.Val < 0.05 & abs(results$logFC) >= 1, ]
cat("  Significant DEGs (padj < 0.05, |log2FC| >= 1):", nrow(sig_degs), "\n\n")

# ============================================================================
# Step 4: Select top DEGs for heatmap
# ============================================================================
cat("Step 4: Selecting top DEGs for heatmap...\n")

# Select top 50 DEGs
if(nrow(sig_degs) >= 50) {
  top_degs <- sig_degs[1:50, ]
  n_genes <- 50
} else {
  top_degs <- sig_degs
  n_genes <- nrow(sig_degs)
}

cat("  Using top", n_genes, "DEGs for visualization\n\n")

# Get expression values
top_gene_ids <- rownames(top_degs)
heatmap_data <- exprs_data[top_gene_ids, ]

# Z-score normalization
heatmap_data_scaled <- t(scale(t(heatmap_data)))

# ============================================================================
# Step 5: Create annotation
# ============================================================================
cat("Step 5: Preparing heatmap annotations...\n")

annotation_col <- data.frame(
  Group = groups,
  row.names = colnames(heatmap_data_scaled)
)

annotation_colors <- list(
  Group = c(TBI = "#E41A1C", Control = "#377EB8")
)

# ============================================================================
# Step 6: Generate clustered heatmap
# ============================================================================
cat("Step 6: Generating clustered heatmap...\n")

png("output.png", width = 10, height = 12, units = "in", res = 300)

pheatmap(
  heatmap_data_scaled,
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  show_rownames = FALSE,
  show_colnames = TRUE,
  fontsize = 10,
  fontsize_row = 8,
  fontsize_col = 9,
  main = "Top DEGs in TBI vs Control",
  border_color = NA
)

dev.off()

cat("  ✓ Heatmap saved as output.png\n\n")

# ============================================================================
# Step 7: Export results to CSV
# ============================================================================
cat("Step 7: Exporting results to CSV...\n")

output_df <- data.frame(
  Gene_ID = rownames(top_degs),
  logFC = top_degs$logFC,
  AveExpr = top_degs$AveExpr,
  t_statistic = top_degs$t,
  P_Value = top_degs$P.Value,
  Adj_P_Value = top_degs$adj.P.Val,
  B_statistic = top_degs$B,
  stringsAsFactors = FALSE
)

# Add gene symbols if available
if("Gene.symbol" %in% colnames(fdata)) {
  output_df$Gene_Symbol <- fdata[rownames(top_degs), "Gene.symbol"]
} else if("GENE_SYMBOL" %in% colnames(fdata)) {
  output_df$Gene_Symbol <- fdata[rownames(top_degs), "GENE_SYMBOL"]
}

output_df <- output_df[order(output_df$Adj_P_Value), ]

write.csv(output_df, "output.csv", row.names = FALSE)

cat("  ✓ Results saved to output.csv\n")
cat("  Total DEGs exported:", nrow(output_df), "\n\n")

cat("  Summary:\n")
cat("    Upregulated (logFC > 1):", sum(output_df$logFC > 1), "\n")
cat("    Downregulated (logFC < -1):", sum(output_df$logFC < -1), "\n")
cat("    Mean |logFC|:", round(mean(abs(output_df$logFC)), 2), "\n\n")

cat("============================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("============================================================================\n")
cat("Outputs:\n")
cat("  - output.csv: Top DEGs with statistics\n")
cat("  - output.png: Clustered heatmap visualization\n")
cat("============================================================================\n")
