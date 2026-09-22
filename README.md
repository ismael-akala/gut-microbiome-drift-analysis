# Quantifying Taxonomic Drift and Alpha Diversity Collapse in Bioreactor-Derived Human Gut Microbiota
https://www.researchsquare.com/article/rs-11118628/v1
A reproducible statistical pipeline and data validation workflow tracking community restructuring, alpha diversity decay, and key biomarker abundance shifts during the host-to-vessel transition in ex vivo gut modeling.

## Data Grounding & Attribution
The raw data matrix analyzed in this pipeline originates from the public supplementary datasets published by **Bressuire et al. (2025)**: *Structural differences in prebiotic-resistant dextrins drive human gut microbiota and glycoside hydrolase profiles in ex vivo continuous fermentations*. 

## Project Objective
This independent re-analysis isolates the initial culture bottleneck that occurs when fresh human stool samples are stabilized in an ex vivo continuous fermentation environment. The pipeline automates data extraction, processes ecological distribution metrics, handles taxonomic missing values, and tracks the survival profile of specific core gut biomarkers.

## Pipeline Architecture & Directory
*   `microbiome_analysis.R`: The complete automated script executing data ingestion from Excel formats, string trimming, sample alignment, matrix numeric conversions, non-parametric hypothesis testing, and 300 DPI visualization generation.
*   `Table1_Summary_Stats.csv`: Automated export mapping the Medians and Interquartile Ranges (IQR) for community indices and target bacterial strains.
*   `Table1_Wilcoxon_Results.csv`: Statistical summaries containing the computed W-statistics and p-values from two-sided Wilcoxon rank-sum evaluations.
*   `Figure1_Alpha_Diversity.png`: High-resolution, multi-panel boxplot (Observed Richness vs. Shannon Diversity Index) generated via ggplot2 and patchwork.

## Tracked Target Biomarkers
In addition to global alpha diversity metrics, the script dynamically monitors the abundance profiles of a curated cohort of critical human gut commensals and generalist markers:
*   *Escherichia coli*
*   *Bacteroides ovatus*
*   *Faecalibacterium prausnitzii* (Strains 1 & 4)
*   *Akkermansia muciniphila*
*   *Bifidobacterium longum*

## Key Methodological Workflow
1.  **Metadata Alignment:** Extracts sample IDs and matches columns dynamically with the experimental sample legend using vector mapping.
2.  **Missing Data Strategy:** Zero-fills missing entries (`is.na`) following numeric matrix conversion to ensure complete cases for row total calculations.
3.  **Iterative Statistical Loops:** Uses automated formula mapping to run comparative statistical aggregates and tests across all indices and biological markers simultaneously, eliminating repetitive scripting.

## Data, Code, and Manuscript Availability Statement
The metagenomic species pangenome abundance matrix analyzed in this study was derived from Bressuire et al. (2025). All custom R scripts, statistical outputs, and summary tables generated during this analysis are publicly archived on GitHub at https://github.com/ismael-akala/gut-microbiome-drift-analysis. The definitive text of this technical note is hosted as an open-access preprint via Research Square at https://www.researchsquare.com/article/rs-11118628/v1, while the version-controlled computing environment and data assets are permanently archived via Zenodo at https://doi.org/10.5281/zenodo.22863070.



