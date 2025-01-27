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

# Do test here: on a subset of data
# ahbmngGenoIDD <- ahbmngGenoID[1:500, 1:1000]
# ahbmngPoss <- ahbmngPos[1:999, ]

# ---- Simulation Parameters ----
num_replications <- 1
sampleSizes <- 100 # c(100,500, 1000, 1600)
colonyH2s <- c(0.5, 0.2, 0.05)
nQtlsPerChr <- 10 #c(1, 10,100)
nChips <- c(1,2)

nMelN           =   400
nWorkers        =   100
nDrones         =   30
nFathers        =   10
nDronesPerQueen =   20
nChip100 <- 100
nChip1000 <- 1000
# nChip10000 <- 10000
colonyH2 <- 0.5

# ---- Section 4. GWAS scenario: Simulation Loops ----
for (rep in 1:num_replications) {
  for (sampleSize in sampleSizes) {
    for (nChip in nChips) {
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
      
      varE <- c(10, 10 / SP$nWorkers, 10, 10 / SP$nWorkers) # change values here for h2 scenario 10=0.5,
      varA / (varA + varE)
      corE <- matrix(data = c(1.0, 0.3, 0.0, 0.0,
                              0.3, 1.0, 0.0, 0.0,
                              0.0, 0.0, 1.0, 0.2,
                              0.0, 0.0, 0.2, 1.0), nrow = 4, byrow = TRUE)
      SP$setVarE(varE = varE, corE = corE)
      
      ## ---- Section 4.6. Add SNP chip ----
      SP$addSnpChip(nChip100)
      SP$addSnpChip(nChip1000)
      #SP$addSnpChip(nChip10000)
      
      ## ---- Section 4.7. Save assigned QTLs
      # qtl <- colnames(pullQtlGeno(ahbFounderGenomeHaplo, simParam=SP))
      
      print("Status: SP global simulation completed")
      
      # Log details
      cat("SampleSize:", sampleSize, "Replication:", rep, "\n")
      cat("Completed nQtlPerChr:", nQtlPerChr, "\n\n")
      
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

      ## ---- Section 5.5. Calculate colony heritability ----

      diag(var(calc_gv_age_0)) / diag(var(calc_pheno_age_0))

      # ---- Section 6. Structure phenotype data ----

      phenoColony <- calc_pheno_age_0
      phenoColony$QID <- seq(1:nrow(phenoColony))
      phenoColony$WID <- seq(nrow(phenoColony) + 1, length.out = nrow(phenoColony))

      # ---- Section 7. - Structure queenGenotype data ----      
      
      genoQ <- getSnpGeno(age_0, caste = "queen", snpChip = nChip, collapse = TRUE)
      rownames(genoQ) <-phenoColony$QID
      phenoColony$ID <- factor(phenoColony$QID, levels = rownames(genoQ))
      colnames(phenoColony)[1] <- "yield"
      phenoColony$QID <- as.factor(phenoColony$QID)
      phenoColony$WID <- as.factor(phenoColony$WID)  
      
      # ---- Section 8. Perform GWAS for Queen ---- 
      
      # To be used for degrees of freedom
      n <- nrow(phenoColony)
      # To be used for degrees of freedom (number of levels in fixed effects)
      k <- 1
      
      ## ---- Section 8.1. remove monomorphic loci, center M matrix and make GRM ---- 
      
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
      # grmQ <- calcBeeGRMIbs(x = genoQclean, sex = rep("F", nrow(genoQclean)))
      grmQTemp <- Gmatrix(SNPmatrix=genoQ, 
                      missingValue=-9, 
                      maf=0.01, 
                      method="VanRaden")
      
      # Make a semi-positive definite
      grmQ <- grmQTemp + diag(nrow(genoQ))*0.000001
      is.positive.definite(grmQ) # has to be TRUE
      grmQinv <- solve(grmQ + diag(1e-6, ncol(grmQ), ncol(grmQ)))
      MTgrmQinv<- t(genoQ_centered)%*%grmQinv
      
      # ---- Section 9. Fit a queen model----
      ## ---- Section 9.1. Structure pheno and geno data
      phenoColony$idQ <- phenoColony$QID 
      phenoColony$QID <- factor(as.character(phenoColony$idQ, levels = rownames(grmQ)))
      QgwasGBLUP <- mmer(yield~1,
                         random=~vsr(QID, Gu=grmQ),
                         rcov=~units, nIters=3,
                         verbose = FALSE,
                         data=phenoColony)
      
      ## ---- Section 9.2. Back solve and calculate SNP effect, snpSE and Pvalues ----
      
      a.from.gQ <-MTgrmQinv%*%matrix(QgwasGBLUP$U$`u:QID`$yield,ncol=1)
      plot(genoQclean %*% a.from.gQ, QgwasGBLUP$U$`u:QID`$yield); abline(a=0,b=1)
      var.gQ <- kronecker(grmQ,QgwasGBLUP$sigma$`u:QID`)- QgwasGBLUP$PevU$`u:QID`$yield
      var.a.from.gQ <- t(genoQ_centered)%*%grmQinv%*% (var.gQ) %*% t(grmQinv)%*%genoQ_centered
      se.a.from.gQ <- sqrt(diag(var.a.from.gQ))
      t.stat.from.gQ <- a.from.gQ/se.a.from.gQ
      QpvalGBLUP <- dt(t.stat.from.gQ,df=n-k-1) 
      
      ## ---- Section 9.3. Structure GWAS summary stat ----
      
      QgwasGBLUPsummary <- data.frame(markerName = colnames(genoQclean), 
                                      Beta = a.from.gQ, 
                                      t_stat = t.stat.from.gQ, 
                                      SE = se.a.from.gQ, 
                                      PValue = QpvalGBLUP)
      
      ## ---- Section 9.4. Collate GWAS summary stat ----
      
      QgwasGBLUPsummary$index <- seq(1:ncol(genoQclean))
      snpchipMap <- merge(ahbmngPos, QgwasGBLUPsummary, by = "markerName")
      snpchipMapOrd <- snpchipMap[order(as.integer(snpchipMap$index)), ]
      snpchipMapOrd$position <- snpchipMapOrd$position
      QgwasGBLUP_summary_tmp <- snpchipMapOrd[, c("markerName", "chromosome", "position", "PValue")]
      QgwasGBLUP_summary <- paste0("QgwasSummary_", 
                                   "nSample_", sampleSize, 
                                   "_nQTL_",nQtlPerChr,
                                   "_h2_", colonyH2,
                                   "_nChip_",nChip100, 
                                   "_nRep_", rep, ".csv")
      fwrite(QgwasGBLUP_summary_tmp, file = QgwasGBLUP_summary)
      
      ## ---- Section 9.5. Visualize GWAS man plot ----
      # CMplot(QgwasGBLUP_summary_tmp,type="p",
      #        plot.type="m",LOG10=TRUE,
      #        threshold=NULL,file="jpg",
      #        file.name="",dpi=300,
      #        file.output=FALSE,
      #        verbose=TRUE,
      #        width =14,height=6,
      #        chr.labels.angle=0,
      #        main = paste0("nSample_", sampleSize,
      #                      "_nQTL_",nQtlPerChr,
      #                      "_h2_", colonyH2,
      #                      "_nChip_",nChip,
      #                      "_nRep_", rep))
      
      ## ---- Section 9.6. Visualize GWAS QQ-plot ----
      # CMplot(QgwasGBLUP_summary_tmp,
      #        plot.type="q",
      #        box=FALSE,
      #        file="jpg",
      #        file.name="",
      #        dpi=300,
      #        conf.int=TRUE,
      #        conf.int.col=NULL,
      #        threshold.col="red",
      #        threshold.lty=2,
      #        file.output=FALSE,
      #        verbose=TRUE,
      #        width=5,
      #        height=5, 
      #        main = paste0("nSample_", sampleSize,
      #                     "_nQTL_",nQtlPerChr,
      #                     "_h2_", colonyH2,
      #                     "_nChip_",nChip,
      #                     "_nRep_", rep))  

# ---- Section 10. Fit a worker model ----     
## ---- Section 10.1. Pull worker genotype ----    
getChipsName <- colnames(getSnpGeno(age_0, snpChip = nChip, caste = "workers",collapse = T))
poolWorkerChips <- do.call(rbind, lapply(X = getSegSiteGeno(x = age_0, caste = "workers", nInd = 10), 
                                         FUN = getPooledGeno, type = "mean"))[,getChipsName]
genoW <- round(poolWorkerChips, digits = 0)

## ---- Section 10.2. Set a threshold and manage monomorphic loci ---- 
genoW.ready <- raw.data(data = as.matrix(genoW), 
                        frame = "wide", 
                        base = FALSE, 
                        sweep.sample = 0.5, 
                        call.rate = 0.95, 
                        maf = 0.01, 
                        imput = FALSE)

dim(genoW.ready$M.clean)
## ---- Section 10.3. Store a clean pooled worker data for down stream analysis ---- 
genoWclean <- genoW.ready$M.clean
genoW_centered <- scale(genoWclean, center = TRUE, scale = FALSE)
grmWTemp <- Gmatrix(SNPmatrix=genoW, 
                    missingValue=-9, 
                    maf=0.00, 
                    method="VanRaden")

# Make a semi-positive definite
grmW <- grmWTemp + diag(nrow(genoW))*0.000001
is.positive.definite(grmW) # has to be TRUE
grmWinv <- solve(grmW + diag(1e-6, ncol(grmW), ncol(grmW)))
MTgrmWinv <- t(genoW_centered)%*%grmWinv

## ---- Section 10.4. Fit a queen model
phenoColony$idW <- seq(nrow(phenoColony) + 1, length.out = nrow(phenoColony))
rownames(grmW) <- phenoColony$idW
colnames(grmW) <- phenoColony$idW
phenoColony$WID <- factor(as.character(phenoColony$idW, levels = rownames(grmW)))

WgwasGBLUP <- mmer(yield~1,
                   random=~vsr(WID, Gu=grmW),
                   rcov=~units, nIters=3,
                   verbose = FALSE,
                   data=phenoColony)

## ---- Section 10.5. Back solve and calculate SNP effect, snpSE and Pvalues ----
a.from.gW <- MTgrmWinv%*%matrix(WgwasGBLUP$U$`u:WID`$yield,ncol=1)
plot(genoWclean %*% a.from.gW, WgwasGBLUP$U$`u:WID`$yield); abline(a=0,b=1)
var.gW <- kronecker(grmW,WgwasGBLUP$sigma$`u:WID`)- WgwasGBLUP$PevU$`u:WID`$yield
var.a.from.gW <- t(genoW_centered)%*%grmWinv%*% (var.gW) %*% t(grmWinv)%*%genoW_centered
se.a.from.gW <- sqrt(diag(var.a.from.gW))
t.stat.from.gW <- a.from.gW/se.a.from.gW
WpvalGBLUP <- dt(t.stat.from.gW,df=n-k-1) 

## ---- Section 10.6. Structure worker GWAS summary stat ----
WgwasGBLUPsummary <- data.frame(markerName = colnames(genoWclean), 
                                Beta = a.from.gW, 
                                t_stat = t.stat.from.gW, 
                                SE = se.a.from.gW, 
                                PValue = WpvalGBLUP)

## ---- Section 10.7. Collate GWAS summary stat ----
WgwasGBLUPsummary$index <- seq(1:ncol(genoWclean))
snpchipMapW <- merge(ahbmngPos, WgwasGBLUPsummary, by = "markerName")
snpchipMapOrdW <- snpchipMapW[order(as.integer(snpchipMapW$index)), ]
snpchipMapOrdW$position <- snpchipMapOrdW$position
WgwasGBLUP_summary_tmp <- snpchipMapOrdW[, c("markerName", "chromosome", "position", "PValue")]
WgwasGBLUP_summary <- paste0("WgwasSummary_", 
                             "nSample_", sampleSize, 
                             "_nQTL_",nQtlPerChr,
                             "_h2_", colonyH2,
                             "_nChip_",nChip100, 
                             "_nRep_", rep, ".csv")
fwrite(WgwasGBLUP_summary_tmp, file = WgwasGBLUP_summary)

# ----- Section 11. Fit Queen + Worker model -----
## ----- Section 11.1. Structure the row levels -----
phenoColony$idQ <- phenoColony$QID 
phenoColony$QID <- factor(as.character(phenoColony$idQ, levels = rownames(grmQ)))
phenoColony$WID <- factor(as.character(phenoColony$idW, levels = rownames(grmW)))

## ----- Section 11.2. Fit a model -----
QWgwasGBLUP <- mmer(yield~1,
                    random=~vsr(QID, Gu=grmQ) + vsr(WID, Gu=grmW),
                    rcov=~units, nIters=3,
                    verbose = FALSE,
                    data=phenoColony)

## ---- Section 11.3. Queen effect from Q + W model: construct GRM and Ginverse and centered Z ----
grmQ <- grmQTemp + diag(nrow(genoQ))*0.000001
is.positive.definite(grmQ) # has to be TRUE
grmQinv <- solve(grmQ + diag(1e-6, ncol(grmQ), ncol(grmQ)))
MTgrmQinv <- t(genoQ_centered)%*%grmQinv

## ---- Section 11.4. Queen effect from Q + W model:Back solve and calculate SNP effect, snpSE and Pvalues ----
a.from.QW.gQ <-MTgrmQinv%*%matrix(QWgwasGBLUP$U$`u:QID`$yield,ncol=1)
plot(genoQclean %*% a.from.QW.gQ, QWgwasGBLUP$U$`u:QID`$yield); abline(a=0,b=1)
var.QW.gQ <- kronecker(grmQ,QWgwasGBLUP$sigma$`u:QID`)- QWgwasGBLUP$PevU$`u:QID`$yield
var.a.from.QW.gQ <- t(genoQ_centered)%*%grmQinv%*% (var.QW.gQ) %*% t(grmQinv)%*%genoQ_centered
se.a.from.QW.gQ <- sqrt(diag(var.a.from.QW.gQ))
t.stat.from.QW.gQ <- a.from.QW.gQ/se.a.from.QW.gQ
QWpvalGBLUP.Q <- dt(t.stat.from.QW.gQ,df=n-k-1) 

## ---- Section 11.5. Collate GWAS for Queen effect -----
QWgwasGBLUPsummaryQ <- data.frame(markerName = colnames(genoQclean), 
                                  Beta = a.from.QW.gQ, 
                                  t_stat = t.stat.from.QW.gQ, 
                                  SE = se.a.from.QW.gQ, 
                                  PValue = QWpvalGBLUP.Q)

QWgwasGBLUPsummaryQ$index <- seq(1:ncol(sdata_conv_v))
snpchipMapQ <- merge(ahbmngPos, QWgwasGBLUPsummaryQ, by = "markerName")
snpchipMapOrdQ <- snpchipMapW[order(as.integer(snpchipMapW$index)), ]
snpchipMapOrdQ$position <- snpchipMapOrdQ$position
QWgwasGBLUP_summary_Q <- snpchipMapOrdQ[, c("markerName", "chromosome", "position", "Pvalue")]

## ---- Section 11.6. Collate all GWAS summary stat ---- 
QWgwasGBLUPsummary_Q <- paste0("QWgwasSummary_Q", 
                               "nSample_", sampleSize, 
                               "_nQTL_",nQtlPerChr,
                               "_h2_", colonyH2,
                               "_nChip_",nChip100, 
                               "_nRep_", rep, ".csv")
fwrite(QWgwasGBLUP_summary_Q, file = QWgwasGBLUPsummary_Q)

## ---- Section 11.7. Worker effect from Q + W model: construct GRM and Ginverse and centered Z ----
grmW <- grmWTemp + diag(nrow(genoW))*0.000001
is.positive.definite(grmW) # has to be TRUE
grmWinv <- solve(grmW + diag(1e-6, ncol(grmW), ncol(grmW)))
MTgrmWinv <- t(genoW_centered)%*%grmWinv

## ---- Section 11.8. Worker effect from W + W model:Back solve and calculate SNP effect, snpSE and Pvalues ----
a.from.QW.gW <-MTgrmWinv%*%matrix(QWgwasGBLUP$U$`u:WID`$yield,ncol=1)
plot(genoWclean %*% a.from.QW.gW, QWgwasGBLUP$U$`u:WID`$yield); abline(a=0,b=1)
var.QW.gW <- kronecker(grmW,QWgwasGBLUP$sigma$`u:WID`)- QWgwasGBLUP$PevU$`u:WID`$yield
var.a.from.QW.gW <- t(genoW_centered)%*%grmWinv%*% (var.QW.gW) %*% t(grmWinv)%*%genoW_centered
se.a.from.QW.gW <- sWrt(diag(var.a.from.QW.gW))
t.stat.from.QW.gW <- a.from.QW.gW/se.a.from.QW.gW
QWpvalGBLUP.W <- dt(t.stat.from.QW.gW,df=n-k-1) 

## ---- Section 11.9. Collate GWAS for Worker effect -----
QWgwasGBLUPsummaryW <- data.frame(markerName = colnames(genoWclean), 
                                  Beta = a.from.QW.gW, 
                                  t_stat = t.stat.from.QW.gW, 
                                  SE = se.a.from.QW.gW, 
                                  PValue = QWpvalGBLUP.W)

QWgwasGBLUPsummaryW$index <- seW(1:ncol(sdata_conv_v))
snpchipMapW <- merge(ahbmngPos, QWgwasGBLUPsummaryW, by = "markerName")
snpchipMapOrdW <- snpchipMapW[order(as.integer(snpchipMapW$index)), ]
snpchipMapOrdW$position <- snpchipMapOrdW$position
QWgwasGBLUP_summary_W <- snpchipMapOrdW[, c("markerName", "chromosome", "position", "Pvalue")]

## ---- Section 11.10. Collate all GWAS summary stat ---- 
QWgwasGBLUPsummary_W <- paste0("QWgwasSummary_W", 
                                "nSample_", sampleSize, 
                                "_nQTL_",nQtlPerChr,
                                "_h2_", colonyH2,
                                "_nChip_",nChip100, 
                                "_nRep_", rep, ".csv")
fwrite(QWgwasGBLUP_summary_W, file = QWgwasGBLUPsummary_W)

      }
    }
  }
}
 
