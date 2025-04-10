# ---- Section 1 - Load library ----
library(AlphaSimR)
library(SIMplyBee)
library(data.table)
library(parallel)
library(matrixcalc)
library(AGHmatrix)
library(sommer)
library(asreml)
library(snpReady)
library(MASS)
library(stringr)
library(dplyr)
library(Hmisc)
library(boot)

## ---- Section 1.1 Objectives ----
# i)   To simulate colony phenotype using a real genomic data \
# ii)  To dissect the queen and workers “contribution to the GWAS power” \
# iii) To determine sample size (n) to be sequenced to achieve an optimal GWAS power with a lower FDR \
# iv)  Revisiting p-value threshold to be used in honey bee GWAS 
# 
# # ---- Section 2. BeeGWAS Steps ----
# i)   Import phased haplotype genomic data and initializing phenotype simulation and \
# ii)  Make a cross \
# iii) PoolGenotype of worker \ 
# iv)  Then estimate colony level phenotype \
# v)  Export simulated phenotype and match genotype data \
# vi)  Perform GWAS for Queen and Worker and Queen + Worker using RRBLUP or GBLUP approach
# vii) Results from vi) has to be backsolved to get marker effect and GWAS summary stat  

# ---- Section 2 - Set workspace ----

setwd("/Users/mulusewfikere/Documents/managed AHB/data/gregor/")

# ---- Section 3 - Read genotype and map data ----
## ---- Section 3.1 - Read mapFile ----
ahbmngPos <- as.data.frame(fread("rawData/haplotypeMapFile.pos", h = FALSE))
colnames(ahbmngPos) <- c("chromosome", "position")
ahbmngPos$test <- sequence(rle(as.character(ahbmngPos$chromosome))$lengths)
ahbmngPos <- data.frame(markerName = paste(ahbmngPos$chromosome, ahbmngPos$test, sep = "_"), ahbmngPos)
ahbmngPos$test <- NULL
ahbmngPos$position <- ahbmngPos$position / 1000000
ahbmngPos$chromosome <- as.numeric(ahbmngPos$chromosome)

## ---- Section 3.2. - Read haplotype file ----
ahbmngGeno <- t(as.data.frame(fread("rawData/haplotypeGenoFile", h = FALSE)))
colnames(ahbmngGeno) <- ahbmngPos$markerName
rownames(ahbmngGeno) <- seq(1:nrow(ahbmngGeno))
ID <- rep(1:(nrow(ahbmngGeno) / 2), each = 2)
ahbmngGenoID <- as.data.frame(cbind(ID, ahbmngGeno))

## ---- Section 3.3. Simulation Parameters ----
# num_replications <- 10
# sampleSizes <- c(100,500, 1000, 1600)
# colonyH2s <- c(10, 45,150) # corresponds to "colonyH2s <- c(0.5, 0.2, 0.05)"
# nQtlsPerChr <- c(1, 10,100)
# nChips <- c(1,2,3)

num_replications <- 1
sampleSizes <- 500 #c(100,500, 1000, 1600)
colonyH2s <- 10 # c(10, 45,150) corresponds to "colonyH2s <- c(0.5, 0.2, 0.05)"
nQtlsPerChr <- 10 #c(1, 10,100)
nChips <- 1 # c(1,2,3)

nMelN           =   400
nWorkers        =   100
nDrones         =   30
nFathers        =   10
nDronesPerQueen =   20
nChip100 <- 100
nChip1000 <- 1000
nChip10000 <- 10000

set.seed(123456)

## ---- Section 3.4. Create a data frame to store results ----
Qgwas_fdr_results <- data.frame(Scenario = character(), TP = numeric(), FP = numeric(), FDR = numeric(), stringsAsFactors = FALSE)
Wgwas_fdr_results <- data.frame(Scenario = character(), TP = numeric(), FP = numeric(), FDR = numeric(), stringsAsFactors = FALSE)
QWgwas_Q_fdr_results <- data.frame(Scenario = character(), TP = numeric(), FP = numeric(), FDR = numeric(), stringsAsFactors = FALSE)
QWgwas_W_fdr_results <- data.frame(Scenario = character(), TP = numeric(), FP = numeric(), FDR = numeric(), stringsAsFactors = FALSE)

Qgwas_power_results <- data.frame(Scenario = character(), Power = numeric(), stringsAsFactors = FALSE)
Wgwas_power_results <- data.frame(Scenario = character(), Power = numeric(), stringsAsFactors = FALSE)
QWgwas_Q_power_results <- data.frame(Scenario = character(), Power = numeric(), stringsAsFactors = FALSE)
QWgwas_W_power_results <- data.frame(Scenario = character(), Power = numeric(), stringsAsFactors = FALSE)

intermidiate_scenario_QgwasGBLUP_summary <- data.frame()
intermidiate_scenario_WgwasGBLUP_W_summary <- data.frame()
intermidiate_scenario_QWgwasGBLUP_Q_summary <- data.frame()
intermidiate_scenario_QWgwasGBLUP_W_summary <- data.frame()

Qgwas_pmin_dist_results <- data.frame(
  GWAS_File = character(),
  Total_Markers = numeric(),
  Psig = numeric(),
  CI_Lower = numeric(),
  CI_Upper = numeric(),
  Neff = numeric(),
  Ratio = numeric(),
  stringsAsFactors = FALSE)

Wgwas_pmin_dist_results <- data.frame(
  GWAS_File = character(),
  Total_Markers = numeric(),
  Psig = numeric(),
  CI_Lower = numeric(),
  CI_Upper = numeric(),
  Neff = numeric(),
  Ratio = numeric(),
  stringsAsFactors = FALSE)

QWgwas_Q_pmin_dist_results <- data.frame(
  GWAS_File = character(),
  Total_Markers = numeric(),
  Psig = numeric(),
  CI_Lower = numeric(),
  CI_Upper = numeric(),
  Neff = numeric(),
  Ratio = numeric(),
  stringsAsFactors = FALSE)

QWgwas_W_pmin_dist_results <- data.frame(
  GWAS_File = character(),
  Total_Markers = numeric(),
  Psig = numeric(),
  CI_Lower = numeric(),
  CI_Upper = numeric(),
  Neff = numeric(),
  Ratio = numeric(),
  stringsAsFactors = FALSE)


# ---- Section 4. GWAS scenario: Simulation Loops ----
for (rep in 1:num_replications) {
  for (sampleSize in sampleSizes) {
    for (nChip in nChips) {
      for (colonyH2 in colonyH2s) {
       for (nQtlPerChr in nQtlsPerChr) {

      ## ---- Section 4.1. Initialize lists ----
      sampled_data <- list()
      ahbFounderGenomeHaplo <- list()
      SP <- list()
      
      ## ---- Section 4.2. Randomly select individuals ----
      sampled_ids <- sample(unique(ahbmngGenoID$ID), size = sampleSize)
      sampled_data <- ahbmngGenoID[ahbmngGenoID$ID %in% sampled_ids, -1]
      
      ## ---- Section 4.3. Import haplotype using SIMplyBee ----
      ahbFounderGenomeHaplo <- importHaplo(
        haplo = sampled_data,
        genMap = ahbmngPos,
        ploidy = 2L,
        ped = NULL)
      
      print("Status: Haplotypes are imported and founderGenome is generated")
      
      ## ---- Section 4.4.  Initialize SP ----
      SP <- SimParamBee$new(ahbFounderGenomeHaplo, csdChr = 3, nCsdAlleles = 100)
      SP$nWorkers <- nWorkers
      SP$nDrones <- nDrones
      SP$nFathers <- nFathers
      SP$swarmP <- 0.5
      SP$splitP <- 0.3
      SP$setTrackPed(TRUE)
      SP$setTrackRec(TRUE)
      
      ## ---- Section 4.5. Define variances ----
      meanP <- c(10, 10 / SP$nWorkers, 0, 0)
      varA <- c(1, 1 / SP$nWorkers, 1, 1 / SP$nWorkers)
      corA <- matrix(data = c(1.0, -0.5, 0.0, 0.0,
                              -0.5, 1.0, 0.0, 0.0,
                              0.0, 0.0, 1.0, -0.4,
                              0.0, 0.0, -0.4, 1.0), nrow = 4, byrow = TRUE)
      SP$addTraitA(nQtlPerChr = nQtlPerChr, mean = meanP, var = varA, corA = corA,
                          name = c("yieldQueenTrait", "yieldWorkersTrait", "calmQueenTrait", "calmWorkersTrait"))
      
      # Change varE values here to achieve colony level heritability
      # 15=0.5,45=0.2,150=0.05
      varE <- c(colonyH2, colonyH2 / SP$nWorkers, colonyH2, colonyH2 / SP$nWorkers) 
      varA / (varA + varE)
      corE <- matrix(data = c(1.0, 0.3, 0.0, 0.0,
                              0.3, 1.0, 0.0, 0.0,
                              0.0, 0.0, 1.0, 0.2,
                              0.0, 0.0, 0.2, 1.0), nrow = 4, byrow = TRUE)
      SP$setVarE(varE = varE, corE = corE)
      
      ## ---- Section 4.6. Add SNP chip ----
      SP$addSnpChip(nChip100)
      SP$addSnpChip(nChip1000)
      SP$addSnpChip(nChip10000)
      
      ## ---- Section 4.7. Save assigned QTLs  ----
      qtl <- colnames(pullQtlGeno(ahbFounderGenomeHaplo, simParam=SP))
      print("Status: SP global simulation completed")
      
      # Log details
      cat(" SampleSize:", sampleSize,"\n", 
          "Replication:", rep, "\n",
          "nQtlPerChr:", nQtlPerChr, "\n",
          "colonyH2s:", colonyH2, "\n",
          "nChips:", nChip100, "\n")
      
      # ---- Section 5. Create virgin queen and father ----
      
      virginQueens <- createVirginQueens(x = ahbFounderGenomeHaplo, n = nMelN)
      drones = createDrones(virginQueens, nInd = nDronesPerQueen)
      
      ## ---- Section 5.1. Make a cross ----
      
      colonies <- cross(x = virginQueens,
                        drones = drones,
                        crossPlan = "create",
                        spatial = FALSE,
                        nDrones = nFathersPoisson,
                        checkCross = "warning")
      
      ## ---- Section 5.2. Build ----
      
      age_0 <- createMultiColony(x = colonies)
      age_0 <- buildUp(age_0)
      nColony <- nColonies(age_0)
      
      ## ---- Section 5.3. Initialize data storage ----
      genoW <- vector("list", length = nColony)
      genoQ <- vector("list", length = nColony)
      phenoColony <- vector("list", length = nColony)

      ## ---- Section 5.4 - Calculate ColonyPhen ----
      calc_pheno_age_0 <- calcColonyPheno(age_0,queenTrait = c("yieldQueenTrait", "calmQueenTrait"),
                                          workersTrait = c("yieldWorkersTrait", "calmWorkersTrait"),
                                          traitName = c("yield", "calmness"),
                                          checkProduction = c(TRUE, FALSE)) |> as.data.frame()

      calc_gv_age_0 <- calcColonyGv(age_0,queenTrait = c("yieldQueenTrait", "calmQueenTrait"),
                                    workersTrait = c("yieldWorkersTrait", "calmWorkersTrait"),
                                    traitName = c("yield", "calmness")) |> as.data.frame()

      cor(calc_pheno_age_0$yield,calc_gv_age_0$yield)
      plot(calc_pheno_age_0$yield ~ calc_gv_age_0$yield, pch=25, col="maroon")
      
      ## ---- Section 5.5. Calculate colony level heritability ----

      diag(var(calc_gv_age_0)) / diag(var(calc_pheno_age_0))

      # ---- Section 6. Structure phenotype data and queen genotype data ----
      phenoColony <- calc_pheno_age_0
      phenoColony$QID <- seq(1:nrow(phenoColony))
      phenoColony$WID <- seq(nrow(phenoColony) + 1, length.out = nrow(phenoColony))

      # ---- Section 7. Structure queen genotype data ----      
      genoQ <- getSnpGeno(age_0, caste = "queen", snpChip = nChip, collapse = TRUE)
      rownames(genoQ) <-phenoColony$QID
      phenoColony$ID <- factor(phenoColony$QID, levels = rownames(genoQ))
      colnames(phenoColony)[1] <- "yield"
      phenoColony$QID <- as.factor(phenoColony$QID)
      phenoColony$WID <- as.factor(phenoColony$WID)  
      
      # ---- Section 8. Perform Queen GWAS ---- 
      
      # To be used for degrees of freedom
      n <- nrow(phenoColony)
      # To be used for degrees of freedom (number of levels in fixed effects)
      k <- 1
      
      ## ---- Section 8.1. Remove monomorphic loci, center M matrix and make GRM ---- 
  
     genoQ.ready <- raw.data(data = as.matrix(genoQ),
                       frame = "wide",
                       base = FALSE,
                       sweep.sample = 0.5,
                       call.rate = 0.95,
                       maf = 0.01,
                       imput = FALSE)

      dim(genoQ.ready$M.clean)
      genoQclean <- genoQ.ready$M.clean
      genoQ_centered <- scale(genoQclean, center = TRUE, scale = FALSE)
      grmQtmp <- Gmatrix(SNPmatrix=genoQ,
                      missingValue=-9,
                      maf=0.01,
                      method="VanRaden")
      grmQ <- grmQtmp + diag(nrow(grmQtmp))*0.000001
      grmQinv <- solve(grmQ)
      
      ## ---- Section 8.2. - Calculate allele freq for queen genotype data ----
      
      alleleFreq <- calcBeeAlleleFreq(x = genoQclean, sex = rep("F", nrow(genoQclean)))
      Qmaf <- ifelse(alleleFreq > 0.5, 1 - alleleFreq, alleleFreq)
      # Qpop.gen <- popgen(M = genoQ)
      # x <- Qpop.gen$whole$Markers$MAF
      plot(density(Qmaf), main = "Density Plot of MAF", 
           xlab = "MAF", col = "blue", lwd = 2)

      MTgrmQinv <- t(genoQ_centered)%*%grmQinv
      
      # ---- Section 8.3. fitting a queen model in ASREML-r
      Qmodel.asr <- asreml(fixed = yield ~ 1,
                           random = ~ vm(QID, grmQ),
                           data = phenoColony,
                           # ai.loadings = TRUE,
                           # ai.sing = TRUE,
                           workspace = "16gb")
      
      # ---- Section 8.4. Summary of the model ----
      summary(Qmodel.asr)$varcomp
      # Predict SNP effects using the model and Extract variance components and BLUP
      predictions <- predict(Qmodel.asr, classify = "QID", only = "vm(QID, grmQ)", pworkspace = "16gb", vcov = T)
      blup_individuals <- predictions$pvals
      pev_individuals <- as.matrix(predictions$vcov)
      # ---- Section 8.5. Start backsolve to compute SNP effects ----
      Z <- genoQ_centered
      G_inv <- solve(grmQ)
      snp_effects <- t(Z) %*% G_inv %*% blup_individuals$predicted.value
      # Assuming a residual variance (sigma2) from the model
      sigma2 <- Qmodel.asr$sigma^2
      # Standard error of SNP effects and Stabilize inversion
      se_snp_effects <- sqrt(diag(sigma2 * solve(t(Z) %*% G_inv %*% Z)))
      t_stats <- snp_effects / se_snp_effects
      chisq <- (snp_effects / se_snp_effects)^2
      lambda_gc <- median(chisq) / qchisq(0.5, df=1)
      QgwasGBLUP_pvalue <- 2 * pt(-abs(t_stats), df = nrow(phenoColony) - 1)
      QgwasGBLUP_pvalue_upd <- pchisq((snp_effects / se_snp_effects)^2 / lambda_gc, df=1, lower.tail=FALSE)
      
      
      # ---- Section 9. Fit a queen model Using mmer package ----
      QgwasGBLUP <- mmer(yield~1,
                     random=~vsr(QID, Gu=grmQ),
                     rcov=~units, nIters=3,
                     verbose = FALSE,
                     data=phenoColony)

      ## ---- Section 9.1. Back solve and calculate SNP effect, snpSE and Pvalues ----

      a.from.gQ <-MTgrmQinv%*%matrix(QgwasGBLUP$U$`u:QID`$yield,ncol=1)
      plot(genoQclean %*% a.from.gQ, QgwasGBLUP$U$`u:QID`$yield); abline(a=0,b=1)
      var.gQ <- kronecker(grmQ,QgwasGBLUP$sigma$`u:QID`)- QgwasGBLUP$PevU$`u:QID`$yield
      var.a.from.gQ <- t(genoQ_centered)%*%grmQinv%*% (var.gQ) %*% t(grmQinv)%*%genoQ_centered
      se.a.from.gQ <- sqrt(diag(var.a.from.gQ))
      t.stat.from.gQ <- a.from.gQ/se.a.from.gQ
      QpvalGBLUP <- dt(t.stat.from.gQ,df=n-k-1)
      
      ## ---- Section 9.2. Structure GWAS summary stat ----
      QgwasGBLUPsummary <- data.frame(markerName = colnames(genoQclean), 
                                      Beta = a.from.gQ, 
                                      SE = se.a.from.gQ, 
                                      P = QpvalGBLUP,
                                      MAF=Qmaf,
                                      N=sampleSize)
      
      ## ---- Section 9.3. Collate GWAS summary stat ----
      QgwasGBLUPsummary$index <- seq(1:ncol(genoQclean))
      snpchipMap <- merge(ahbmngPos, QgwasGBLUPsummary, by = "markerName")
      QgwasGBLUP_summary <- snpchipMap[order(as.integer(snpchipMap$index)), ]
      QgwasGBLUP_summary$position <- QgwasGBLUP_summary$position
     
      
      # ---- Section 10. Calculate FDR for Queen scenario ----
    
      ## ---- Section 10.1. Extract the predefined QTL from founderGenome and SP ----
      predefined_QTLlist <- as.data.frame(colnames(pullQtlGeno(ahbFounderGenomeHaplo, simParam=SP)))
      colnames(predefined_QTLlist) <- "markerName"
      predefined_QTLs <- as.data.frame(setDT(predefined_QTLlist)[, c("chromosome", "position") := tstrsplit(markerName, "_")])
      predefined_QTLs <- predefined_QTLs %>% mutate(position = as.numeric(position))
      QgwasGBLUP_summary <- QgwasGBLUP_summary %>% mutate(position = as.numeric(position))
      
      ## ---- Section 10.2. Define a thresholds ----
      p_threshold <- 10^(-6)
      qtl_window_size <- 5000 
      gwas_significant <- QgwasGBLUP_summary %>% filter(P < p_threshold)
      
      ## ---- Section 10.3. Check if significant SNPs fall within ±5k window of predefined QTL positions ----
      predefined_QTLs <- predefined_QTLs %>%
      rowwise() %>%
      mutate(is_significant = any(gwas_significant$chromosome == chromosome &
                                      gwas_significant$position >= (position - qtl_window_size) &
                                      gwas_significant$position <= (position + qtl_window_size))) %>%  ungroup() 
      
      ## ---- Section 10.4. Extract significant QTLs ----
      significant_qtl <- predefined_QTLs %>% filter(is_significant)
      ## ---- Section 10.5. Count True Positives (TP): significant SNPs that are predefined QTLs ----
      TP <- nrow(gwas_significant) + nrow(significant_qtl)
      ## ---- Section 10.6. Count False Positives (FP): significant SNPs that are NOT predefined QTLs ----
      FP <- abs(nrow(gwas_significant) - nrow(significant_qtl))
      ## ---- Section 10.7. Calculate false discovery rate (FDR) ----
      FDR <- ifelse((FP + TP) > 0, FP / (FP + TP), NA)
      ## ---- Section 10.8. Create a base file name ----
      QgwasGBLUP_summary_fdr <- paste0("QgwasSummary_",
                                         "nSample_", sampleSize,
                                         "_nQTL_",nQtlPerChr,
                                         "_h2_", colonyH2,
                                         "_nChip_",nChip,
                                         "_nRep_", rep)
        
      ## ---- Section 10.9. Save FDR results in all scenario ----
      Qgwas_fdr_results <- rbind(Qgwas_fdr_results, data.frame(Scenario = basename(QgwasGBLUP_summary_fdr), TP = TP, FP = FP, FDR = FDR))
      Qgwas_fdr_results[,1] <- gsub("h2_15", "h2_0.5", Qgwas_fdr_results[,1])
      Qgwas_fdr_results[,1] <- gsub("h2_45", "h2_0.2", Qgwas_fdr_results[,1])
      Qgwas_fdr_results[,1] <- gsub("h2_150", "h2_0.05", Qgwas_fdr_results[,1])
      
      # ---- Section 11. Calculate GWAS power for Queen scenario ----
    
      ## ---- Section 11.1. Set GWAS power parameters ----
      N <- QgwasGBLUP_summary$N[1]
      QgwasGBLUP_summary$Power <- apply(QgwasGBLUP_summary, 1, function(row) {
      MAF <- as.numeric(row["MAF"])
      Beta <- as.numeric(row["Beta"])
      SE <- as.numeric(row["SE"])
          
      ## ---- Section 11.2. Calculate Non-centrality Parameter (NCP) ----
      NCP <- (Beta / SE) * sqrt(N * 2 * MAF * (1 - MAF))
      Z_critical <- qnorm(1 - p_threshold / 2)
      power <- pnorm(NCP - Z_critical) + (1 - pnorm(NCP + Z_critical))
      return(power)
      })
    
      # Calculate mean GWAS power across all SNPs for this replication
      overall_power <- mean(QgwasGBLUP_summary$Power, na.rm = TRUE)
      se_power <- sd(QgwasGBLUP_summary$Power, na.rm = TRUE) / sqrt(sum(!is.na(QgwasGBLUP_summary$Power)))
    
      
      # Generate scenario file name
      QgwasGBLUP_summary_power <- paste0("QgwasSummary_",
                                         "nSample_", sampleSize,
                                         "_nQTL_", nQtlPerChr,
                                         "_h2_", colonyH2,
                                         "_nChip_", nChip,
                                         "_nRep_", rep)
          
      #  ---- Section 12. Save the power results for each scenario ----
      Qgwas_power_results <- rbind(Qgwas_power_results,data.frame(Scenario = basename(QgwasGBLUP_summary_power),
                                                      Power = overall_power, 
                                                      SE = se_power))
      
          Qgwas_power_results[,1] <- gsub("h2_15", "h2_0.5", Qgwas_power_results[,1])
          Qgwas_power_results[,1] <- gsub("h2_45", "h2_0.2", Qgwas_power_results[,1])
          Qgwas_power_results[,1] <- gsub("h2_150", "h2_0.05", Qgwas_power_results[,1])
          
      # ---- Section 12. Pmin distribution for queen model ----
          ## ---- Section 12.1. Pull the intermediate scenario of queen GWAS model ---- 
                        # Intermediate scenarios: h2=0.2=45; nQtlPerChr=10; nChip=3=10k (10kPerChr)
      
      # Scenario matches the target conditions
          if (colonyH2 == 45 & nQtlPerChr == 10 & nChip == 3) {
            
      # Use the intermediate QgwasGBLUP_summary for downstream analysis
        intermidiate_scenario_QgwasGBLUP_summary <- QgwasGBLUP_summary
        total_markers <- nrow(intermidiate_scenario_QgwasGBLUP_summary)
          
      # Extract and transform all p-values
        log_pvals <- -log10(intermidiate_scenario_QgwasGBLUP_summary$P)
          
      # Compute the 95th percentile using the Harrell–Davis estimator
        psig_log <- hdquantile(log_pvals, 0.95)
          
      # Convert back to p-value scale
        psig <- 10^(-psig_log)
          
      # Bootstrapping function for confidence intervals
        boot_fun <- function(data, indices) {
        return(hdquantile(data[indices], 0.95))
          }
          
      # Perform bootstrapping using the full vector of log-transformed p-values
        boot_results <- boot(data = log_pvals, statistic = boot_fun, R = 1000)
        ci <- boot.ci(boot_results, type = "perc")$percent[4:5]
          
      # Convert CI back to p-value scale
        ci_pval <- 10^(-ci)
          
      # Estimate effective number of independent variants
        alpha <- 0.05
        neff <- alpha / psig
        
      # Calculate the ratio of independent variants to total markers
        ratio <- neff / total_markers
          
      # Save results to CSV file
        pmin_dist_results <- data.frame(
          QgwasGBLUP_summary,
          Total_Markers = total_markers,
          Psig = psig,
          CI_Lower = ci_pval[1],
          CI_Upper = ci_pval[2],
          Neff = neff,
          Ratio = ratio)
                }
  # ----- Section 13. END OF QUEEN MODEL ----
  #           }
  #           
  #        }
  #       }
  #     }
  #   }
  # }

  # ----- Section 14. WORKER MODEL START HERE ----
        ## ---- Section 14.1. Pull worker genotype ----    
        getChipsName <- colnames(getSnpGeno(age_0, snpChip = nChip, caste = "workers",collapse = T))
        genoW <- do.call(rbind, lapply(X = getSegSiteGeno(x = age_0, caste = "workers"), 
                                       FUN = getPooledGeno, type = "mean"))[,getChipsName]
        # genoW <- round(poolWorkerChips, digits = 0)
        
        ## ---- Section 14.2. Set a threshold and manage monomorphic loci ---- 
        genoW.ready <- raw.data(data = as.matrix(genoW),
                                frame = "wide",
                                base = FALSE,
                                sweep.sample = 0.5,
                                call.rate = 0.95,
                                maf = 0.01,
                                imput = FALSE)
        genoWclean <- genoW.ready$M.clean
        
        genoW_centered <- scale(genoWclean, center = TRUE, scale = FALSE)
        grmWtmp <- Gmatrix(SNPmatrix=genoW,
                           missingValue=-9,
                           maf=0.01,
                           method="VanRaden")
        grmW <- grmWtmp + diag(nrow(grmWtmp))*0.000001
        grmWinv <- solve(grmW)
        
        ## ---- Section 14.3. - Calculate allele freq for Wueen genotype data ----
        
        alleleFreW <- calcBeeAlleleFreq(x = genoWclean, sex = rep("F", nrow(genoWclean)))
        Wmaf <- ifelse(alleleFreW > 0.5, 1 - alleleFreW, alleleFreW)
        # Wpop.gen <- popgen(M = genoW)
        # x <- Wpop.gen$whole$Markers$MAF
        plot(density(Wmaf), main = "Density Plot of MAF", 
             xlab = "MAF", col = "blue", lwd = 2)
        
        MTgrmWinv <- t(genoW_centered)%*%grmWinv
        
        ## ---- Section 14.4. Fit a Wueen model Using mmer package ----
        WgwasGBLUP <- mmer(yield~1,
                           random=~vsr(WID, Gu=grmW),
                           rcov=~units, nIters=3,
                           verbose = FALSE,
                           data=phenoColony)
        
        ## ---- Section 14.5. Back solve and calculate SNP effect, snpSE and Pvalues ----
        
        a.from.gW <-MTgrmWinv%*%matrix(WgwasGBLUP$U$`u:WID`$yield,ncol=1)
        plot(genoWclean %*% a.from.gW, WgwasGBLUP$U$`u:WID`$yield); abline(a=0,b=1)
        var.gW <- kronecker(grmW,WgwasGBLUP$sigma$`u:WID`)- WgwasGBLUP$PevU$`u:WID`$yield
        var.a.from.gW <- t(genoW_centered)%*%grmWinv%*% (var.gW) %*% t(grmWinv)%*%genoW_centered
        se.a.from.gW <- sWrt(diag(var.a.from.gW))
        t.stat.from.gW <- a.from.gW/se.a.from.gW
        WpvalGBLUP <- dt(t.stat.from.gW,df=n-k-1)
        
        ## ---- Section 14.6. Structure GWAS summary stat ----
        WgwasGBLUPsummary <- data.frame(markerName = colnames(genoWclean), 
                                        Beta = a.from.gW, 
                                        SE = se.a.from.gW, 
                                        P = WpvalGBLUP,
                                        MAF=Wmaf,
                                        N=sampleSize)
        
        ## ---- Section 14.7. Collate GWAS summary stat ----
        WgwasGBLUPsummary$index <- seq(1:ncol(genoWclean))
        snpchipMap <- merge(ahbmngPos, WgwasGBLUPsummary, by = "markerName")
        WgwasGBLUP_summary <- snpchipMap[order(as.integer(snpchipMap$index)), ]
        WgwasGBLUP_summary$position <- WgwasGBLUP_summary$position
        
        # ---- Section 15. Calculate FDR for worker scenario ----
        
        ## ---- Section 15.1. Extract the predefined QTL from founderGenome and SP ----
        WgwasGBLUP_summary <- WgwasGBLUP_summary %>% mutate(position = as.numeric(position))
        
        ## ---- Section 15.2. Define a thresholds ----
        p_threshold <- 10^(-6)
        qtl_window_size <- 5000 
        Wgwas_significant <- WgwasGBLUP_summary %>% filter(P < p_threshold)
        
        ## ---- Section 15.3. Check if significant SNPs fall within ±5k window of predefined QTL positions ----
        predefined_QTLs <- predefined_QTLs %>%
          rowwise() %>%
          mutate(is_significant = any(Wgwas_significant$chromosome == chromosome &
                                        Wgwas_significant$position >= (position - qtl_window_size) &
                                        Wgwas_significant$position <= (position + qtl_window_size))) %>%  ungroup() 
        
        ## ---- Section 15.4. Extract significant QTLs ----
        significant_qtl <- predefined_QTLs %>% filter(is_significant)
        ## ---- Section 15.5. Count True Positives (TP): significant SNPs that are predefined QTLs ----
        TP <- nrow(Wgwas_significant) + nrow(significant_qtl)
        ## ---- Section 15.6. Count False Positives (FP): significant SNPs that are NOT predefined QTLs ----
        FP <- abs(nrow(Wgwas_significant) - nrow(significant_qtl))
        ## ---- Section 15.7. Calculate false discovery rate (FDR) ----
        FDR <- ifelse((FP + TP) > 0, FP / (FP + TP), NA)
        ## ---- Section 15.8. Create a base file name ----
        WgwasGBLUP_summary_fdr <- paste0("WgwasSummary_",
                                         "nSample_", sampleSize,
                                         "_nQTL_",nQtlPerChr,
                                         "_h2_", colonyH2,
                                         "_nChip_",nChip,
                                         "_nRep_", rep)
        
        ## ---- Section 15.9. Save FDR results in all scenario and write a *.csv after the loop ----
        Wgwas_fdr_results <- rbind(Wgwas_fdr_results, data.frame(Scenario = basename(WgwasGBLUP_summary_fdr), TP = TP, FP = FP, FDR = FDR))
        Wgwas_fdr_results[,1] <- gsub("h2_15", "h2_0.5", Wgwas_fdr_results[,1])
        Wgwas_fdr_results[,1] <- gsub("h2_45", "h2_0.2", Wgwas_fdr_results[,1])
        Wgwas_fdr_results[,1] <- gsub("h2_150", "h2_0.05", Wgwas_fdr_results[,1])
      
        # ---- Section 16. Calculate GWAS power for Queen scenario ----
        
        ## ---- Section 16.1. Set GWAS power parameters ----
        N <- WgwasGBLUP_summary$N[1]
        WgwasGBLUP_summary$Power <- apply(WgwasGBLUP_summary, 1, function(row) {
          MAF <- as.numeric(row["MAF"])
          Beta <- as.numeric(row["Beta"])
          SE <- as.numeric(row["SE"])
          
          ## ---- Section 16.2. Calculate Non-centrality Parameter (NCP) ----
          NCP <- (Beta / SE) * sqrt(N * 2 * MAF * (1 - MAF))
          Z_critical <- qnorm(1 - p_threshold / 2)
          power <- pnorm(NCP - Z_critical) + (1 - pnorm(NCP + Z_critical))
          return(power)
        })
        
        # Calculate mean GWAS power across all SNPs for this replication
        overall_power <- mean(WgwasGBLUP_summary$Power, na.rm = TRUE)
        se_power <- sd(WgwasGBLUP_summary$Power, na.rm = TRUE) / sqrt(sum(!is.na(WgwasGBLUP_summary$Power)))
        
        
        # Generate scenario file name
        WgwasGBLUP_summary_power <- paste0("WgwasSummary_",
                                           "nSample_", sampleSize,
                                           "_nQTL_", nQtlPerChr,
                                           "_h2_", colonyH2,
                                           "_nChip_", nChip,
                                           "_nRep_", rep)
        
        #  ---- Section 17. Save the power results for each scenario ----
        Wgwas_power_results <- rbind(Wgwas_power_results,data.frame(Scenario = basename(WgwasGBLUP_summary_power),
                                                                    Power = overall_power, 
                                                                    SE = se_power))
        
        Wgwas_power_results[,1] <- gsub("h2_15", "h2_0.5", Wgwas_power_results[,1])
        Wgwas_power_results[,1] <- gsub("h2_45", "h2_0.2", Wgwas_power_results[,1])
        Wgwas_power_results[,1] <- gsub("h2_150", "h2_0.05", Wgwas_power_results[,1])
        
        # ---- Section 17. Pmin distribution for queen model ----
        ## ---- Section 17.1. Pull the intermediate scenario of queen GWAS model ---- 
        # Intermediate scenarios: h2=0.2=45; nQtlPerChr=10; nChip=3=10k (10kPerChr)
        
        # Scenario matches the target conditions
        if (colonyH2 == 45 & nQtlPerChr == 10 & nChip == 3) {
          
          # Use the intermediate WgwasGBLUP_summary for downstream analysis
          intermidiate_scenario_WgwasGBLUP_summary <- WgwasGBLUP_summary
          total_markers <- nrow(intermidiate_scenario_WgwasGBLUP_summary)
          
          # Extract and transform all p-values
          log_pvals <- -log10(intermidiate_scenario_WgwasGBLUP_summary$P)
          
          # Compute the 95th percentile using the Harrell–Davis estimator
          psig_log <- hdquantile(log_pvals, 0.95)
          
          # Convert back to p-value scale
          psig <- 10^(-psig_log)
          
          # Bootstrapping function for confidence intervals
          boot_fun <- function(data, indices) {
            return(hdquantile(data[indices], 0.95))
          }
          
          # Perform bootstrapping using the full vector of log-transformed p-values
          boot_results <- boot(data = log_pvals, statistic = boot_fun, R = 1000)
          ci <- boot.ci(boot_results, type = "perc")$percent[4:5]
          
          # Convert CI back to p-value scale
          ci_pval <- 10^(-ci)
          
          # Estimate effective number of independent variants
          alpha <- 0.05
          neff <- alpha / psig
          
          # Calculate the ratio of independent variants to total markers
          ratio <- neff / total_markers
          
          # Save results to CSV file
          Wgwas_pmin_dist_results <- data.frame(
            WgwasGBLUP_summary,
            Total_Markers = total_markers,
            Psig = psig,
            CI_Lower = ci_pval[1],
            CI_Upper = ci_pval[2],
            Neff = neff,
            Ratio = ratio)
      # ----- Section 18. END OF WORKER MODEL ----
          }
#           }
#       
#         }
#       }
#     }
#   }
# }
  
  # ----- Section 19. WORKER-QUEEN MODEL START HERE ----
          ## ----- Section 19.1. Structure the row levels -----
          phenoColony$idQ <- phenoColony$QID 
          phenoColony$QID <- factor(as.character(phenoColony$idQ, levels = rownames(grmQ)))
          phenoColony$WID <- factor(as.character(phenoColony$idW, levels = rownames(grmW)))
          
          ## ----- Section 19.2. Fit a model -----
          QWgwasGBLUP <- mmer(yield~1,
                              random=~vsr(QID, Gu=grmQ) + vsr(WID, Gu=grmW),
                              rcov=~units, nIters=3,
                              verbose = FALSE,
                              data=phenoColony)    
          
          ## ---- Section 19.3. Queen effect from Q + W model:Back solve ----
          a.from.QW.gQ <-MTgrmQinv%*%matrix(QWgwasGBLUP$U$`u:QID`$yield,ncol=1)
          plot(genoQclean %*% a.from.QW.gQ, QWgwasGBLUP$U$`u:QID`$yield); abline(a=0,b=1)
          var.QW.gQ <- kronecker(grmQ,QWgwasGBLUP$sigma$`u:QID`)- QWgwasGBLUP$PevU$`u:QID`$yield
          var.a.from.QW.gQ <- t(genoQ_centered)%*%grmQinv%*% (var.QW.gQ) %*% t(grmQinv)%*%genoQ_centered
          se.a.from.QW.gQ <- sqrt(diag(var.a.from.QW.gQ))
          t.stat.from.QW.gQ <- a.from.QW.gQ/se.a.from.QW.gQ
          QWpvalGBLUP.Q <- dt(t.stat.from.QW.gQ,df=n-k-1) 
          
          ## ---- Section 19.4. Collate GWAS for Queen effect -----
          QWgwas_Q_GBLUPsummary <- data.frame(markerName = colnames(genoQclean), 
                                              Beta = a.from.QW.gQ, 
                                              SE = se.a.from.QW.gQ, 
                                              P = QWpvalGBLUP.Q,
                                              MAF=Qmaf,
                                              N=sampleSize)
          
          ## ---- Section 19.5. Collate GWAS summary stat ----
          QWgwas_Q_GBLUPsummary$index <- seq(1:ncol(genoWclean))
          snpchipMap <- merge(ahbmngPos, QWgwas_Q_GBLUPsummary, by = "markerName")
          QWgwasGBLUP_summary <- snpchipMap[order(as.integer(snpchipMap$index)), ]
          QWgwasGBLUP_summary$position <- QWgwasGBLUP_summary$position
          
          # ---- Section 20. Calculate FDR for worker scenario ----
          
          ## ---- Section 20.1. Extract the predefined QTL from founderGenome and SP ----
          QWgwasGBLUP_summary <- QWgwasGBLUP_summary %>% mutate(position = as.numeric(position))
          
          ## ---- Section 20.2. Define a thresholds ----
          QWgwas_Q_significant <- QWgwasGBLUP_summary %>% filter(P < p_threshold)
          
          ## ---- Section 20.3. Check if significant SNPs fall within ±5k window of predefined QTL positions ----
          predefined_QTLs <- predefined_QTLs %>%
            rowwise() %>%
            mutate(is_significant = any(QWgwas_Q_significant$chromosome == chromosome &
                                          QWgwas_Q_significant$position >= (position - qtl_window_size) &
                                          QWgwas_Q_significant$position <= (position + qtl_window_size))) %>%  ungroup() 
          
          ## ---- Section 20.4. Extract significant QTLs ----
          significant_qtl <- predefined_QTLs %>% filter(is_significant)
          ## ---- Section 20.5. Count True Positives (TP): significant SNPs that are predefined QTLs ----
          TP <- nrow(QWgwas_Q_significant) + nrow(significant_qtl)
          ## ---- Section 20.6. Count False Positives (FP): significant SNPs that are NOT predefined QTLs ----
          FP <- abs(nrow(QWgwas_Q_significant) - nrow(significant_qtl))
          ## ---- Section 20.7. Calculate false discovery rate (FDR) ----
          FDR <- ifelse((FP + TP) > 0, FP / (FP + TP), NA)
          ## ---- Section 20.8. Create a base file name ----
          QWgwasGBLUP_Q_summary_fdr <- paste0("QWgwasSummary_Q_",
                                              "nSample_", sampleSize,
                                              "_nQTL_",nQtlPerChr,
                                              "_h2_", colonyH2,
                                              "_nChip_",nChip,
                                              "_nRep_", rep)
          
          ## ---- Section 20.9. Save FDR results in all scenario and write a *.csv after the loop ----
          QWgwas_Q_fdr_results <- rbind(QWgwas_Q_fdr_results, data.frame(Scenario = basename(QWgwasGBLUP_Q_summary_fdr), TP = TP, FP = FP, FDR = FDR))
          QWgwas_Q_fdr_results[,1] <- gsub("h2_15", "h2_0.5", QWgwas_Q_fdr_results[,1])
          QWgwas_Q_fdr_results[,1] <- gsub("h2_45", "h2_0.2", QWgwas_Q_fdr_results[,1])
          QWgwas_Q_fdr_results[,1] <- gsub("h2_150", "h2_0.05", QWgwas_Q_fdr_results[,1])
          # ---- Section 21. Calculate GWAS power for Queen scenario ----
          
          ## ---- Section 21.1. Set GWAS power parameters ----
          N <- QWgwasGBLUP_summary$N[1]
          QWgwasGBLUP_summary$Power <- apply(QWgwasGBLUP_summary, 1, function(row) {
            MAF <- as.numeric(row["MAF"])
            Beta <- as.numeric(row["Beta"])
            SE <- as.numeric(row["SE"])
            
            ## ---- Section 21.2. Calculate Non-centrality Parameter (NCP) ----
            NCP <- (Beta / SE) * sqrt(N * 2 * MAF * (1 - MAF))
            Z_critical <- qnorm(1 - p_threshold / 2)
            power <- pnorm(NCP - Z_critical) + (1 - pnorm(NCP + Z_critical))
            return(power)
          })
          
          # Calculate mean GWAS power across all SNPs for this replication
          overall_power <- mean(QWgwasGBLUP_summary$Power, na.rm = TRUE)
          se_power <- sd(QWgwasGBLUP_summary$Power, na.rm = TRUE) / sqrt(sum(!is.na(QWgwasGBLUP_summary$Power)))
          
          
          # Generate scenario file name
          QWgwasGBLUP_Q_summary_power <- paste0("QWgwasSummary_Q_",
                                                "nSample_", sampleSize,
                                                "_nQTL_", nQtlPerChr,
                                                "_h2_", colonyH2,
                                                "_nChip_", nChip,
                                                "_nRep_", rep)
          
          #  ---- Section 22. Save the power results for each scenario ----
          QWgwas_Q_power_results <- rbind(QWgwas_Q_power_results,data.frame(Scenario = basename(QWgwasGBLUP_Q_summary_power),
                                                                            Power = overall_power, 
                                                                            SE = se_power))
          
          QWgwas_Q_power_results[,1] <- gsub("h2_15", "h2_0.5", QWgwas_Q_power_results[,1])
          QWgwas_Q_power_results[,1] <- gsub("h2_45", "h2_0.2", QWgwas_Q_power_results[,1])
          QWgwas_Q_power_results[,1] <- gsub("h2_150", "h2_0.05", QWgwas_Q_power_results[,1])
          
          # ---- Section 23. Pmin distribution for queen model ----
          ## ---- Section 23.1. Pull the intermediate scenario of queen GWAS model ---- 
          # Intermediate scenarios: h2=0.2=45; nQtlPerChr=10; nChip=3=10k (10kPerChr)
          
          # Scenario matches the target conditions
          if (colonyH2 == 45 & nQtlPerChr == 10 & nChip == 3) {
            
            # Use the intermediate QWgwasGBLUP_summary for downstream analysis
            intermidiate_scenario_QWgwasGBLUP_Q_summary <- QWgwasGBLUP_summary
            total_markers <- nrow(intermidiate_scenario_QWgwasGBLUP_Q_summary)
            
            # Extract and transform all p-values
            log_pvals <- -log10(intermidiate_scenario_QWgwasGBLUP_Q_summary$P)
            
            # Compute the 95th percentile using the Harrell–Davis estimator
            psig_log <- hdquantile(log_pvals, 0.95)
            
            # Convert back to p-value scale
            psig <- 10^(-psig_log)
            
            # Bootstrapping function for confidence intervals
            boot_fun <- function(data, indices) {
              return(hdquantile(data[indices], 0.95))
            }
            
            # Perform bootstrapping using the full vector of log-transformed p-values
            boot_results <- boot(data = log_pvals, statistic = boot_fun, R = 1000)
            ci <- boot.ci(boot_results, type = "perc")$percent[4:5]
            
            # Convert CI back to p-value scale
            ci_pval <- 10^(-ci)
            
            # Estimate effective number of independent variants
            alpha <- 0.05
            neff <- alpha / psig
            
            # Calculate the ratio of independent variants to total markers
            ratio <- neff / total_markers
            
            # Save results to CSV file
            QWgwas_Q_pmin_dist_results <- data.frame(
              QWgwasGBLUP_summary,
              Total_Markers = total_markers,
              Psig = psig,
              CI_Lower = ci_pval[1],
              CI_Upper = ci_pval[2],
              Neff = neff,
              Ratio = ratio)
            # ----- Section 24. END OF QUEEN + WORKER (Q) MODEL ----
               }
            
 #           }
 #           
 #         }
 #      }
 #     }
 #   }
 # } 
          # ----- Section 25. QUEEN + WORKER (W) MODEL STARTS HERE ----
          ## ---- Section 25.1. Queen effect from Q + W model:Back solve ----
          a.from.QW.gW <-MTgrmQinv%*%matrix(QWgwasGBLUP$U$`u:WID`$yield,ncol=1)
          plot(genoWclean %*% a.from.QW.gW, QWgwasGBLUP$U$`u:WID`$yield); abline(a=0,b=1)
          var.QW.gW <- kronecker(grmQ,QWgwasGBLUP$sigma$`u:WID`)- QWgwasGBLUP$PevU$`u:WID`$yield
          var.a.from.QW.gW <- t(genoQ_centered)%*%grmQinv%*% (var.QW.gW) %*% t(grmQinv)%*%genoQ_centered
          se.a.from.QW.gW <- sqrt(diag(var.a.from.QW.gW))
          t.stat.from.QW.gW <- a.from.QW.gW/se.a.from.QW.gW
          QWpvalGBLUP.W <- dt(t.stat.from.QW.gW,df=n-k-1) 
          
          ## ---- Section 25.2. Collate GWAS for Queen effect -----
          QWgwas_W_GBLUPsummary <- data.frame(markerName = colnames(genoWclean), 
                                              Beta = a.from.QW.gW, 
                                              SE = se.a.from.QW.gW, 
                                              P = QWpvalGBLUP.W,
                                              MAF=Qmaf,
                                              N=sampleSize)
          
          ## ---- Section 25.3. Collate GWAS summary stat ----
          QWgwas_W_GBLUPsummary$index <- seq(1:ncol(genoWclean))
          snpchipMap <- merge(ahbmngPos, QWgwas_W_GBLUPsummary, by = "markerName")
          QWgwasGBLUP_summary <- snpchipMap[order(as.integer(snpchipMap$index)), ]
          QWgwasGBLUP_summary$position <- QWgwasGBLUP_summary$position
          
          # ---- Section 25.4. Calculate FDR for worker scenario ----
          
          ## ---- Section 25.6. Extract the predefined QTL from founderGenome and SP ----
          QWgwasGBLUP_summary <- QWgwasGBLUP_summary %>% mutate(position = as.numeric(position))
          
          ## ---- Section 25.7. Define a thresholds ----
          QWgwas_W_significant <- QWgwasGBLUP_summary %>% filter(P < p_threshold)
          
          ## ---- Section 25.8. Check if significant SNPs fall within ±5k window of predefined QTL positions ----
          predefined_QTLs <- predefined_QTLs %>%
            rowwise() %>%
            mutate(is_significant = any(QWgwas_W_significant$chromosome == chromosome &
                                          QWgwas_W_significant$position >= (position - qtl_window_size) &
                                          QWgwas_W_significant$position <= (position + qtl_window_size))) %>%  ungroup() 
          
          ## ---- Section 25.9. Extract significant QTLs ----
          significant_qtl <- predefined_QTLs %>% filter(is_significant)
          ## ---- Section 25.10. Count True Positives (TP): significant SNPs that are predefined QTLs ----
          TP <- nrow(QWgwas_W_significant) + nrow(significant_qtl)
          ## ---- Section 25.11. Count False Positives (FP): significant SNPs that are NOT predefined QTLs ----
          FP <- abs(nrow(QWgwas_W_significant) - nrow(significant_qtl))
          ## ---- Section 25.12. Calculate false discovery rate (FDR) ----
          FDR <- ifelse((FP + TP) > 0, FP / (FP + TP), NA)
          ## ---- Section 25.13. Create a base file name ----
          QWgwasGBLUP_W_summary_fdr <- paste0("QWgwasSummary_W_",
                                              "nSample_", sampleSize,
                                              "_nQTL_",nQtlPerChr,
                                              "_h2_", colonyH2,
                                              "_nChip_",nChip,
                                              "_nRep_", rep)
          
          ## ---- Section 25.14. Save FDR results in all scenario and write a *.csv after the loop ----
          QWgwas_W_fdr_results <- rbind(QWgwas_W_fdr_results, data.frame(Scenario = basename(QWgwasGBLUP_W_summary_fdr), TP = TP, FP = FP, FDR = FDR))
          QWgwas_W_fdr_results[,1] <- gsub("h2_15", "h2_0.5", QWgwas_W_fdr_results[,1])
          QWgwas_W_fdr_results[,1] <- gsub("h2_45", "h2_0.2", QWgwas_W_fdr_results[,1])
          QWgwas_W_fdr_results[,1] <- gsub("h2_150", "h2_0.05", QWgwas_W_fdr_results[,1])
            
          # ---- Section 26. Calculate GWAS power for Worker scenario ----
          
          ## ---- Section 26.1. Set GWAS power parameters ----
          N <- QWgwasGBLUP_summary$N[1]
          QWgwasGBLUP_summary$Power <- apply(QWgwasGBLUP_summary, 1, function(row) {
            MAF <- as.numeric(row["MAF"])
            Beta <- as.numeric(row["Beta"])
            SE <- as.numeric(row["SE"])
            
            ## ---- Section 26.2. Calculate Non-centrality Parameter (NCP) ----
            NCP <- (Beta / SE) * sqrt(N * 2 * MAF * (1 - MAF))
            Z_critical <- qnorm(1 - p_threshold / 2)
            power <- pnorm(NCP - Z_critical) + (1 - pnorm(NCP + Z_critical))
            return(power)
          })
          
          # Calculate mean GWAS power across all SNPs for this replication
          overall_power <- mean(QWgwasGBLUP_summary$Power, na.rm = TRUE)
          se_power <- sd(QWgwasGBLUP_summary$Power, na.rm = TRUE) / sqrt(sum(!is.na(QWgwasGBLUP_summary$Power)))
          
          
          # Generate scenario file name
          QWgwasGBLUP_W_summary_power <- paste0("QWgwasSummary_W_",
                                                "nSample_", sampleSize,
                                                "_nQTL_", nQtlPerChr,
                                                "_h2_", colonyH2,
                                                "_nChip_", nChip,
                                                "_nRep_", rep)
          
          #  ---- Section 27. Save the power results for each scenario ----
          QWgwas_W_power_results <- rbind(QWgwas_W_power_results,data.frame(Scenario = basename(QWgwasGBLUP_W_summary_power),
                                                                            Power = overall_power, 
                                                                            SE = se_power))
          
          QWgwas_W_power_results[,1] <- gsub("h2_15", "h2_0.5", QWgwas_W_power_results[,1])
          QWgwas_W_power_results[,1] <- gsub("h2_45", "h2_0.2", QWgwas_W_power_results[,1])
          QWgwas_W_power_results[,1] <- gsub("h2_150", "h2_0.05", QWgwas_W_power_results[,1])
          
          # ---- Section 28. Pmin distribution for queen model ----
          ## ---- Section 29.1. Pull the intermediate scenario of queen GWAS model ---- 
          # Intermediate scenarios: h2=0.2=45; nQtlPerChr=10; nChip=3=10k (10kPerChr)
          
          # Scenario matches the target conditions
          if (colonyH2 == 45 & nQtlPerChr == 10 & nChip == 3) {
            
            # Use the intermediate QWgwasGBLUP_summary for downstream analysis
            intermidiate_scenario_QWgwasGBLUP_W_summary <- QWgwasGBLUP_summary
            total_markers <- nrow(intermidiate_scenario_QWgwasGBLUP_W_summary)
            
            # Extract and transform all p-values
            log_pvals <- -log10(intermidiate_scenario_QWgwasGBLUP_W_summary$P)
            
            # Compute the 95th percentile using the Harrell–Davis estimator
            psig_log <- hdquantile(log_pvals, 0.95)
            
            # Convert back to p-value scale
            psig <- 10^(-psig_log)
            
            # Bootstrapping function for confidence intervals
            boot_fun <- function(data, indices) {
              return(hdquantile(data[indices], 0.95))
            }
            
            # Perform bootstrapping using the full vector of log-transformed p-values
            boot_results <- boot(data = log_pvals, statistic = boot_fun, R = 1000)
            ci <- boot.ci(boot_results, type = "perc")$percent[4:5]
            
            # Convert CI back to p-value scale
            ci_pval <- 10^(-ci)
            
            # Estimate effective number of independent variants
            alpha <- 0.05
            neff <- alpha / psig
            
            # Calculate the ratio of independent variants to total markers
            ratio <- neff / total_markers
            
            # Save results to CSV file
            QWgwas_W_pmin_dist_results <- data.frame(
              QWgwasGBLUP_summary,
              Total_Markers = total_markers,
              Psig = psig,
              CI_Lower = ci_pval[1],
              CI_Upper = ci_pval[2],
              Neff = neff,
              Ratio = ratio)
              }
          }
       }
      }
    }
  } 
  
  # ---- SECTION 30. END OF THE LOOP ----
  # ---- SECTION 31. Save result summary for each model and scenarios ----
  ##  ---- Section 31.1. Store pmin results ----
      write.csv(intermidiate_scenario_QgwasGBLUP_summary, file="Qgwas_Pmin_summary.csv", row.names = F)
      write.csv(intermidiate_scenario_WgwasGBLUP_summary, file="Qgwas_W_Pmin_summary.csv", row.names = F)
      write.csv(intermidiate_scenario_QWgwasGBLUP_Q_summary, file="QWgwas_Q_Pmin_summary.csv", row.names = F)
      write.csv(intermidiate_scenario_QWgwasGBLUP_W_summary, file="Qgwas_W_Pmin_summary.csv", row.names = F)

  #  ---- Section 31.2.Store fdr results ----
      write.csv(Qgwas_fdr_results, file="Qgwas_fdr_results.csv", row.names = F)
      write.csv(Wgwas_fdr_results, file="Wgwas_fdr_results.csv", row.names = F)
      write.csv(QWgwas_Q_fdr_results, file="QWgwas_Q_fdr_results.csv", row.names = F)
      write.csv(QWgwas_W_fdr_results, file="Qgwas_W_fdr_results.csv", row.names = F)
  
  #  ---- Section 31.3.Store power results ----
      write.csv(Qgwas_power_results, file="Qgwas_power_results.csv", row.names = F)
      write.csv(Wgwas_power_results, file="Wgwas_power_results.csv", row.names = F)
      write.csv(QWgwas_Q_power_results, file="QWgwas_Q_power_results.csv", row.names = F)
      write.csv(QWgwas_W_power_results, file="QWgwas_W_power_results.csv", row.names = F)

  #  ---- Section 31.4. Store pmin_dist_results ----
      write.csv(Qgwas_pmin_dist_results, "Qgwas_GWAS_Psig_Results.csv", row.names = FALSE)
      write.csv(Wgwas_pmin_dist_results, "Wgwas_GWAS_Psig_Results.csv", row.names = FALSE)
      write.csv(QWgwas_Q_pmin_dist_results, "QWgwas_Q_GWAS_Psig_Results.csv", row.names = FALSE)
      write.csv(QWgwas_W_pmin_dist_results, "QWgwas_W_GWAS_Psig_Results.csv", row.names = FALSE)

