# ==============================================================================
# Pipeline: Metagenomic Re-Analysis & Automated Statistical Validation Workflow
# Author: Ismael Muriuki Akala, MSc.
# Project: Quantifying Alpha Diversity Collapse & Biomarker Profiling in 
#          Continuous-Flow Ex Vivo Gut Microbiota Fermentations
# Data Grounding: Secondary analysis derived from Bressuire et al. (2025) 
#                 Metagenomic Species Pan-genome (MSP) supplementary data.
# Objective: Automatably assess taxonomic "culture bottlenecks" by comparing
#            in vivo fresh stool baselines against stabilized ex vivo vessels.
# ==============================================================================

# 1. Load libraries
library(readxl)
library(vegan)

# 2. Load Excel file and clean legend names
excel_file <- list.files(pattern = "\\.xlsx$", full.names = TRUE)[1]
raw_legend <- as.data.frame(read_excel(excel_file, sheet = "Sample legend"))
colnames(raw_legend) <- trimws(colnames(raw_legend))

# 3. Align sample IDs and assign origin
sample_ids <- sapply(strsplit(rownames(transposed_df), "_"), function(x) paste(x[1:2], collapse="_"))
transposed_df[["origin"]] <- as.factor(raw_legend[["origin"]][match(sample_ids, raw_legend[["sample ID"]])])

# 4. Filter to fresh stool vs bioreactor sample
transposed_df <- transposed_df[transposed_df[["origin"]] %in% c("fresh stool", "bioreactor sample"), ]
transposed_df[["origin"]] <- factor(transposed_df[["origin"]], levels = c("fresh stool", "bioreactor sample"))

# 5. Extract species matrix & handle NAs
species_cols <- setdiff(colnames(transposed_df), "origin")
species_matrix <- as.matrix(transposed_df[, species_cols])
mode(species_matrix) <- "numeric"
species_matrix[is.na(species_matrix)] <- 0

# 6. Recompute diversity metrics
transposed_df[["Richness"]] <- rowSums(species_matrix > 0)
transposed_df[["Shannon"]] <- diversity(species_matrix, index = "shannon")

# 7. Identify species present in transposed_df
target_biomarkers <- c(
  "msp_0005_Escherichia coli",
  "msp_0007_Bacteroides ovatus",
  "msp_0301_Faecalibacterium prausnitzii 1",
  "msp_0389_Faecalibacterium prausnitzii 4",
  "msp_0025_Akkermansia muciniphila",
  "msp_0166_Bifidobacterium longum"
)

existing_targets <- intersect(target_biomarkers, colnames(transposed_df))

# 8. Rebuild rcmdr_data with diversity AND biomarker species
rcmdr_data <- transposed_df[, c("origin", "Richness", "Shannon", existing_targets)]

# 9. Compute Summary Stats (Medians & IQRs) for all variables
summary_table <- do.call(rbind, lapply(c("Richness", "Shannon", existing_targets), function(sp) {
  res <- aggregate(as.formula(paste("`", sp, "` ~ origin", sep="")), data = rcmdr_data, 
                   FUN = function(x) c(Median = median(x), IQR = IQR(x)))
  data.frame(Variable = sp, Group = res[["origin"]], Median = res[[2]][,1], IQR = res[[2]][,2])
}))

# 10. Compute Wilcoxon Tests (W & p-values) for all variables
wilcox_table <- t(sapply(c("Richness", "Shannon", existing_targets), function(sp) {
  wt <- wilcox.test(as.formula(paste("`", sp, "` ~ origin", sep="")), data = rcmdr_data, alternative = "two.sided")
  c(W = unname(wt[["statistic"]]), p_value = wt[["p.value"]])
}))

# Print results directly to console
print("--- MEDIANS & IQRS ---")
print(summary_table)
print("--- WILCOXON TESTS (W & P-VALUES) ---")
print(wilcox_table)

# Save summary stats (Medians & IQRs) and test results (W & p-values) to CSV
write.csv(summary_table, "Table1_Summary_Stats.csv", row.names = FALSE)
write.csv(wilcox_table, "Table1_Wilcoxon_Results.csv", row.names = TRUE)

# 12. Generate and save Figure 1 (Publication-Ready Boxplots)
library(ggplot2)
library(patchwork)

# Panel A: Observed Richness
p1 <- ggplot(rcmdr_data, aes(x = origin, y = Richness, fill = origin)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 16) +
  geom_jitter(width = 0.1, alpha = 0.5) +
  theme_classic() +
  labs(title = "A. Observed Richness", x = "", y = "Species Count") +
  theme(legend.position = "none")

# Panel B: Shannon Diversity
p2 <- ggplot(rcmdr_data, aes(x = origin, y = Shannon, fill = origin)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 16) +
  geom_jitter(width = 0.1, alpha = 0.5) +
  theme_classic() +
  labs(title = "B. Shannon Diversity", x = "", y = "Shannon Index") +
  theme(legend.position = "none")

# Combine panels and export PNG at 300 DPI
fig1 <- p1 + p2
ggsave("Figure1_Alpha_Diversity.png", plot = fig1, width = 8, height = 4, dpi = 300)