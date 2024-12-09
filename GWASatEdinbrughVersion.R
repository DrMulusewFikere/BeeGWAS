##############################
# 0. Load library
##############################
library(AlphaSimR)
library(SIMplyBee)
library(data.table)
library(asreml)
library(qqman)
library(lme4)

# devtools::install_github(repo="HighlanderLab/SIMplyBee@main")
setwd("/Users/mulusewfikere/Documents/managed AHB/data/gregor/")

##############################
# 1. Read genotype and map data
##############################

      ####################
      # 1.1. Read genMap
      ####################
      ahbmngPos <- as.data.frame(fread("rawData/haplotypeMapFile.pos", h=F))
      colnames(ahbmngPos) <- c("chromosome", "position")
      ahbmngPos$test <- sequence(rle(as.character(ahbmngPos$chromosome))$lengths)
      ahbmngPos <- data.frame(markerName=paste(ahbmngPos$chromosome,ahbmngPos$test, sep="_"), ahbmngPos)
      ahbmngPos$test <- NULL
      ahbmngPos$position <- ahbmngPos$position/1000000 
      ahbmngPos$chromosome <- as.numeric(ahbmngPos$chromosome)

      ####################
      # 1.1. Read Geno file
      ####################
      ahbmngGeno <- t(as.data.frame(fread("rawData/haplotypeGenoFile", h=F)))
      colnames(ahbmngGeno) <- ahbmngPos$markerName
      rownames(ahbmngGeno) <- seq(1:nrow(ahbmngGeno))

##############################
# 2. Create a founder genome
##############################

      ####################
      # 2.1. Import haplotype after phasing VCF file
      ####################
    
      ahbFounderGenomeHaplo = importHaplo(haplo=ahbmngGeno,
                                          genMap=ahbmngPos,
                                          ploidy = 2L,
                                          ped=NULL)
      ahbFounderGenomeHaplo@nInd
      ahbFounderGenomeHaplo@nLoci
      ahbFounderGenomeHaplo@nChr
##############################
# 3. Initialize simulation
##############################
nMelN           =   400
nWorkers        =   100
nDrones         =   30
nFathers        =   10
nDronesPerQueen =   20
nRep = 2
nChip100 = 100
nChip1000 = 1000
nChip10000 = 10000
nQtl = 1 # 1,10, 100
h2= 0.5 # 0.5, 0.2, 0.05

      #################
      # 3.1. Create SP object and write in the global simulation/population parameters
      # Pre-define number of workers, drones, apiary, nQTLs
      #################
      SP <- SimParamBee$new(ahbFounderGenomeHaplo,csdChr = 3, nCsdAlleles = 100 )
      SP$nWorkers <- nWorkers
      SP$nDrones <- nDrones
      SP$nFathers <- nFathers
      SP$nVirginQueens <- nVirginQueens
      SP$swarmP <- 0.5
      SP$splitP <- 0.3
      # Track the pedigree
      SP$setTrackPed(TRUE)
      # Track the recombination
      SP$setTrackRec(TRUE)
      # define csd chromomsome
      csdChr <- SP$csdChr

      #################
      # 3.2. Pre-define variances
      #################
      # Start with h2=0.5
      meanP <- c(10, 10 / SP$nWorkers, 0, 0)
      varA <- c(1, 1 / SP$nWorkers, 1, 1 / SP$nWorkers)
      corA <- matrix(data = c( 1.0, -0.5,  0.0,  0.0,
                               -0.5,  1.0,  0.0,  0.0,
                               0.0,  0.0,  1.0, -0.4,
                               0.0,  0.0, -0.4,  1.0), nrow = 4, byrow = TRUE)
      SP$addTraitA(nQtlPerChr = nQtl, mean = meanP, var = varA, corA = corA,
                   name = c("yieldQueenTrait", "yieldWorkersTrait",
                            "calmQueenTrait", "calmWorkersTrait"))
      varE <- c(1, 1/ SP$nWorkers, 1, 1 / SP$nWorkers) # all 1=0.5, all 4=0.2, 7=0.05
      varA / (varA + varE)
      
      corE <- matrix(data = c(1.0, 0.3, 0.0, 0.0,
                              0.3, 1.0, 0.0, 0.0,
                              0.0, 0.0, 1.0, 0.2,
                              0.0, 0.0, 0.2, 1.0), nrow = 4, byrow = TRUE)
      SP$setVarE(varE = varE, corE = corE)

##############################
# 4. Assign SNPchip and extract QTL
##############################
SP$addSnpChip(nChip100)
SP$addSnpChip(nChip1000)
SP$addSnpChip(nChip10000)
qtl_1 <- colnames(pullQtlGeno(ahbFounderGenomeHaplo, simParam=SP))

##############################
# 5. Initialize data storage
##############################
# 
genoW <- vector("list", length = nColony)
genoQ <- vector("list", length = nColony)
phenoColony <- vector("list", length = nColony)

##############################
# 6. Create virgin queen and father
##############################
      
      #################
      # 6.1. Pre-define variances
      #################
      virginQueens <- createVirginQueens(x = ahbFounderGenomeHaplo, n = nMelN)
      drones = createDrones(virginQueens, nInd = nDronesPerQueen)
      #cross
      colonies <- cross(x = virginQueens,
                        drones = drones,
                        crossPlan = "create",
                        spatial = FALSE,
                        nDrones = nFathersPoisson,
                        checkCross = "warning")

      #################
      # 6.2. Build
      #################
      age_0 <- createMultiColony(x = colonies)
      age_0 <- buildUp(age_0)
      nColony <- nColonies(age_0)

      #################
      # 6.2. Calculate ColonyPheno
      #################
      calc_pheno_age_0 <- calcColonyPheno(age_0,queenTrait = c("yieldQueenTrait", "calmQueenTrait"),
                                    workersTrait = c("yieldWorkersTrait", "calmWorkersTrait"),
                                    traitName = c("yield05", "calmness05"),
                                    checkProduction = c(TRUE, FALSE)) |> as.data.frame()

                ##########
                # 6.2.1. Structure phenotype data
                ###########
                phenoColony <- calc_pheno_age_0
                phenoColony$QID <- seq(1:nrow(phenoColony))
                phenoColony$WID <- seq(nrow(phenoColony) + 1, length.out = nrow(phenoColony))

                ##########
                # 6.2.2. Structure queenGenotype
                ###########
                pooledGenoQ <- getSnpGeno(age_0, caste = "queen", snpChip = 2, collapse = TRUE)
                rownames(pooledGenoQ) <-phenoColony$QID
                
                phenoColony$ID <- factor(phenoColony$QID, levels = rownames(pooledGenoQ))
                colnames(phenoColony)[1] <- "yield05"
                phenoColony$QID <- as.factor(phenoColony$QID)
                phenoColony$WID <- as.factor(phenoColony$WID)  

##############################
# 7. Save RData
##############################
save.image("ahbFounderGenomeHaploAg0.RData")

##############################
# 8. Perform Queen GWAS
##############################

      #################
      # 8.1. Initialize vectors to store results
      #################

      beta_values <- numeric(ncol(pooledGenoQ))
      se_values <- numeric(ncol(pooledGenoQ))
      p_values <- numeric(ncol(pooledGenoQ))

      #################
      # 8.2. Fit model and loop through each SNP
      #################
      
      for (i in 1:ncol(pooledGenoQ)) {
        model <- lm(phenoColony[,1] ~ pooledGenoQ[,i])
        # Summary of the model
        summary_model <- summary(model)
        # Extract beta (coefficient) and standard error (SE)
        beta_values[i] <- summary_model$coefficients[2, 1] 
        se_values[i] <- summary_model$coefficients[2, 2]
        # p-value for the genotype effect
        p_values[i] <- summary_model$coefficients[2, 4]
      }
              ##########
              # 8.2.1. Collate result
              ###########
              Qgwas <- data.frame(markerName = colnames(pooledGenoQ), Beta = beta_values, SE = se_values, PValue = p_values)
              Qgwas$index <- rownames(Qgwas)
              # Get the SNP chip map and make GWAS summary stat
              snpchipMap <- merge(ahbmngPos, Qgwas, by = "markerName")
              snpchipMapOrd <- snpchipMap[order(as.integer(snpchipMap$index)), ]
              snpchipMapOrd$position <- snpchipMapOrd$position * 1e6
              snpchipMapOrdSel <- snpchipMapOrd[, c("markerName", "chromosome", "position", "PValue")]

##############################
# 9. Visualize
##############################
manhattan(
  snpchipMapOrdSel,
  chr = "chromosome",
  bp = "position",
  p = "PValue",
  snp = "markerName",
  main = "Queen GWAS",
  col = c("black", "gray"),
  suggestiveline = FALSE,
  genomewideline = FALSE)

      #################
      # 9.1. Add a text
      #################
 
      mtext(expression("nsnpChip = 1000"), side = 3, line = 3, cex = 0.8, adj = 1, col = "black")
      mtext(expression("nQTL = 1, h2=0.5"), side = 3, line = 2, cex = 0.8, adj = 1, col = "black")
      mtext(expression("nWorker = 100"), side = 3, line = 1, cex = 0.8, adj = 1, col = "black")

      #################
      # 9.2. QQ plot
      #################
      qq(snpchipMapOrdSel$PValue)

##############################
# 10. Pool worker GWAS
##############################

    
      
      
##############################
# 11. Queen and worker GWAS
##############################  
      
      
      
      
      
      
      
      
