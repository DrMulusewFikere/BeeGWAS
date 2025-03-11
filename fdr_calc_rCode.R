#---- Section 1. Objectives ----
# the following code calculates false discovery rate
# Using predefined nQTLs and GWAS summary statistics

#---- Section 2. Load required libraries ---- 
library(dplyr)

# ---- Section 3. Define the folder where GWAS files are stored  ----
gwasResult_folder <- "/gwasOutput/" 
output_file <- "fdr_results.csv"

# ---- Section 4. Read predefined QTLs  ----
predefined_QTLs <- read.csv("predefined_nQTL1.csv")$QTL

## ---- Section 4.1.Define significance threshold  ----
p_threshold <- 5e-6

# ---- Section 5. Read-in GWAS summary files ----
gwas_files <- list.files(gwasResult_folder, full.names = TRUE, pattern = "*.csv")

# ---- Section 6. Store results  ----
fdr_results <- data.frame(File = character(), TP = numeric(), FP = numeric(), FDR = numeric(), stringsAsFactors = FALSE)

## ---- Section 6.1. Loop through each GWAS summary statistics file ----
for (file in gwas_files) {
    gwas_results <- read.csv(file)
  
  # Identify significant SNPs
  significant_snps <- gwas_results %>%
    filter(P_value < p_threshold) %>%
    pull(SNP)
  # Count True Positives (TP): significant SNPs that are predefined QTLs
  TP <- sum(significant_snps %in% predefined_QTLs)
  # Count False Positives (FP): significant SNPs that are NOT predefined QTLs
  FP <- sum(!(significant_snps %in% predefined_QTLs))
 
## ---- Section 6.2. Calculate false discovery rate (FDR) ----
  FDR <- ifelse((FP + TP) > 0, FP / (FP + TP), NA)
  
## ---- Section 6.3. Store results ----
  fdr_results <- rbind(fdr_results, data.frame(File = basename(file), TP = TP, FP = FP, FDR = FDR))
}

# ---- Section 7. Save results to a CSV file
write.csv(fdr_results, output_file, row.names = FALSE)
