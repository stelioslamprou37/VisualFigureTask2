#!/usr/bin/env Rscript

# ============================================================================
# 18-Gene Glycolysis Signature Validation - Bladder Cancer
# Paper: Frontiers in Immunology 2025 - DOI: 10.3389/fimmu.2024.1430583
# Data: GSE48276 with REAL survival data from Supplementary File S1
# Fixed: 71 samples with complete survival and expression data
# ============================================================================

suppressPackageStartupMessages({
  library(GEOquery)
  library(preprocessCore)
  library(survival)
  library(survminer)
  library(timeROC)
  library(ggplot2)
})

cat("============================================================================\n")
cat("18-Gene Glycolysis Signature - Bladder Cancer Analysis\n")
cat("Paper: Frontiers in Immunology 2025\n")
cat("Data: GSE48276 (71 samples with survival data)\n")
cat("============================================================================\n\n")

# ============================================================================
# Step 1: Load expression data from GEO
# ============================================================================

cat("Step 1: Loading GSE48276 expression data from GEO...\n")

gset <- getGEO("GSE48276", GSEMatrix = TRUE, destdir = "./input_data/")
gset <- gset[[1]]

exprs_data_full <- exprs(gset)
pdata_full <- pData(gset)
fdata <- fData(gset)

cat("  ✓ Loaded:", nrow(exprs_data_full), "genes x", ncol(exprs_data_full), "samples\n\n")

# ============================================================================
# Step 2: Define 71 samples with REAL survival data (from Supplementary File S1)
# ============================================================================

cat("Step 2: Loading REAL survival data (71 samples from Supplementary File S1)...\n\n")

# EXACTLY 71 sample IDs
survival_sample_ids <- c(
  "GSM1173901", "GSM1173902", "GSM1173903", "GSM1173904", "GSM1173905",
  "GSM1173906", "GSM1173907", "GSM1173908", "GSM1173909", "GSM1173910",
  "GSM1173911", "GSM1173912", "GSM1173913", "GSM1173914", "GSM1173915",
  "GSM1173916", "GSM1173917", "GSM1173918", "GSM1173919", "GSM1173920",
  "GSM1173921", "GSM1173922", "GSM1173923", "GSM1173924", "GSM1173925",
  "GSM1173926", "GSM1173927", "GSM1173928", "GSM1173929", "GSM1173930",
  "GSM1173931", "GSM1173932", "GSM1173933", "GSM1173934", "GSM1173935",
  "GSM1173936", "GSM1173937", "GSM1173938", "GSM1173939", "GSM1173940",
  "GSM1173941", "GSM1173942", "GSM1173943", "GSM1173944", "GSM1173945",
  "GSM1173946", "GSM1173947", "GSM1173948", "GSM1173949", "GSM1173950",
  "GSM1173951", "GSM1173952", "GSM1173953", "GSM1173954", "GSM1173955",
  "GSM1173956", "GSM1173957",
  "GSM1173968", "GSM1173969", "GSM1173970",
  "GSM1173971", "GSM1173972", "GSM1173973", "GSM1173974", "GSM1173975",
  "GSM1173976", "GSM1173977", "GSM1173978", "GSM1173979", "GSM1173980",
  "GSM1173981"
)

# EXACTLY 71 follow-up times (days)
futime_days_vec <- c(
  705, 726, 5400, 3177, 465, 3753, 2187, 1389, 1689, 1884, 171, 387, 957,
  2064, 2376, 147, 1884, 2079, 516, 1545, 2466, 2142, 507, 1566, 1995, 1170,
  957, 861, 750, 954, 165, 1719, 1755, 1020, 1386, 1410, 1125, 117, 1143,
  4059, 1326, 1596, 222, 144, 288, 1851, 1023, 798, 1305, 807, 780, 1242,
  684, 1395, 1311, 624, 960, 459, 1062, 963, 822, 1620, 285, 1425, 1305,
  1206, 1185, 1161, 372, 903, 762
)

# EXACTLY 71 event status values (1=death, 0=censored)
fustat_vec <- c(
  1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0,
  0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0
)

# Verify all vectors have exactly 71 elements
if(length(survival_sample_ids) != 71 || 
   length(futime_days_vec) != 71 || 
   length(fustat_vec) != 71) {
  stop("ERROR: Survival vectors do not all have 71 elements!")
}

cat("  ✓ Loaded 71 samples with REAL survival data\n")
cat("  Events (deaths):", sum(fustat_vec), "\n")
cat("  Censored:", sum(1 - fustat_vec), "\n\n")

# ============================================================================
# Step 3: Subset expression data to 71 samples with survival data
# ============================================================================

cat("Step 3: Subsetting expression data to 71 samples...\n")

keep_idx <- match(survival_sample_ids, colnames(exprs_data_full))

if(sum(is.na(keep_idx)) > 0) {
  stop("ERROR: Cannot find all survival samples in expression data!")
}

exprs_data <- exprs_data_full[, keep_idx]
pdata <- pdata_full[keep_idx, ]

cat("  ✓ Subset to", ncol(exprs_data), "samples\n\n")

# ============================================================================
# Step 4: Preprocessing
# ============================================================================

cat("Step 4: Preprocessing expression data...\n")

if(max(exprs_data, na.rm = TRUE) > 100) {
  exprs_data <- log2(exprs_data + 1)
}

exprs_normalized <- normalize.quantiles(exprs_data)
rownames(exprs_normalized) <- rownames(exprs_data)
colnames(exprs_normalized) <- colnames(exprs_data)

cat("  ✓ Quantile normalization complete\n\n")

# ============================================================================
# Step 5: Extract 18-gene signature
# ============================================================================

cat("Step 5: Extracting 18-gene glycolysis signature...\n")

signature_genes <- c("CPNE8", "CXCL6", "COMP", "SPINK4", "HLA-DQB2",
                     "CLIC3", "DMRTA1", "MAP2", "ZNF600", "CYTL1",
                     "DIP2C", "FOXC2", "LPXN", "SPINK5", "SLC1A6",
                     "FASN", "SCD", "EGFL6")

gene_column <- "Gene.symbol"
if(!"Gene.symbol" %in% colnames(fdata)) {
  possible_cols <- grep("gene|symbol", colnames(fdata), ignore.case = TRUE, value = TRUE)
  if(length(possible_cols) > 0) gene_column <- possible_cols[1]
}

signature_probes <- unlist(lapply(signature_genes, function(gene) {
  idx <- grep(paste0("^", gene, "$"), fdata[[gene_column]], ignore.case = TRUE)
  if(length(idx) >= 1) rownames(fdata)[idx[1]] else NA
}))

names(signature_probes) <- signature_genes
signature_probes <- na.omit(signature_probes)

cat("  Found", length(signature_probes), "/18 genes\n\n")

# ============================================================================
# Step 6: Calculate risk scores
# ============================================================================

cat("Step 6: Calculating risk scores...\n")

signature_expr <- exprs_normalized[signature_probes, ]
signature_expr_scaled <- t(scale(t(signature_expr)))
risk_scores <- colMeans(signature_expr_scaled, na.rm = TRUE)
risk_group <- ifelse(risk_scores > median(risk_scores), "High", "Low")

cat("  High risk:", sum(risk_group == "High"), "\n")
cat("  Low risk:", sum(risk_group == "Low"), "\n\n")

# ============================================================================
# Step 7: Create survival data frame
# ============================================================================

cat("Step 7: Creating survival data frame...\n")

survival_data <- data.frame(
  id = survival_sample_ids,
  futime_days = futime_days_vec,
  fustat = fustat_vec,
  risk_score = risk_scores,
  risk_group = risk_group,
  stringsAsFactors = FALSE
)

survival_data$futime_years <- survival_data$futime_days / 365

cat("  ✓ Data frame created: ", nrow(survival_data), "samples\n\n")

# ============================================================================
# Step 8: Time-dependent ROC curves
# ============================================================================

cat("Step 8: Computing time-dependent ROC curves...\n")

keep <- !is.na(survival_data$futime_years) & !is.na(survival_data$fustat) & 
        !is.na(survival_data$risk_score)

roc_res <- timeROC(
  T = survival_data$futime_years[keep],
  delta = survival_data$fustat[keep],
  marker = survival_data$risk_score[keep],
  cause = 1,
  weighting = "marginal",
  times = c(1, 3, 5),
  iid = TRUE
)

cat("  AUC at 1-year:", round(roc_res$AUC[1], 3), "\n")
cat("  AUC at 3-year:", round(roc_res$AUC[2], 3), "\n")
cat("  AUC at 5-year:", round(roc_res$AUC[3], 3), "\n\n")

# ============================================================================
# Step 9: Generate plot
# ============================================================================

cat("Step 9: Generating time-dependent ROC plot...\n")

png("output.png", width = 10, height = 8, units = "in", res = 300)

plot_data <- data.frame(
  FPR = c(roc_res$FP[, 1], roc_res$FP[, 2], roc_res$FP[, 3]),
  TPR = c(roc_res$TP[, 1], roc_res$TP[, 2], roc_res$TP[, 3]),
  Time = rep(c("1-year", "3-year", "5-year"), each = nrow(roc_res$FP))
)

g <- ggplot(plot_data, aes(x = FPR, y = TPR, color = Time)) +
  geom_line(size = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  labs(
    title = "Time-Dependent ROC Curves\n18-Gene Glycolysis Signature in Bladder Cancer",
    subtitle = sprintf("AUC: 1-year = %.3f | 3-year = %.3f | 5-year = %.3f",
                       roc_res$AUC[1], roc_res$AUC[2], roc_res$AUC[3]),
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)"
  ) +
  scale_color_manual(values = c("1-year" = "#E41A1C", "3-year" = "#377EB8", "5-year" = "#4DAF4A")) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.title = element_text(size = 12, face = "bold"),
    panel.grid.minor = element_blank()
  )

print(g)
dev.off()

cat("  ✓ Saved output.png\n\n")

# ============================================================================
# Step 10: Export results
# ============================================================================

cat("Step 10: Exporting results to CSV...\n")

output_df <- data.frame(
  Sample_ID = survival_data$id,
  Risk_Score = survival_data$risk_score,
  Risk_Group = survival_data$risk_group,
  Follow_Up_Time_Years = survival_data$futime_years,
  Event_Status = survival_data$fustat
)

output_df <- output_df[order(-output_df$Risk_Score), ]
write.csv(output_df, "output.csv", row.names = FALSE)

cat("  ✓ Saved output.csv\n\n")

cat("============================================================================\n")
cat("ANALYSIS COMPLETE!\n")
cat("============================================================================\n")
cat("Outputs:\n")
cat("  - output.png: Time-dependent ROC curves (300 DPI, publication-ready)\n")
cat("  - output.csv: Risk scores and survival data (71 samples)\n\n")
cat("Data source:\n")
cat("  Expression: GSE48276 (GEO)\n")
cat("  Survival: Supplementary File S1 (from paper)\n")
cat("  Paper: Frontiers in Immunology 2025\n")
cat("  DOI: 10.3389/fimmu.2024.1430583\n")
cat("============================================================================\n")
