##############################
# 0. Load library
##############################
library(AlphaSimR)
library(SIMplyBee)
library(data.table)
library(sommer)
library(CMplot)

# devtools::install_github(repo="HighlanderLab/SIMplyBee@main")
# setwd("/Users/mulusewfikere/Documents/managed\ AHB/data/gregor/")

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
      varE <- c(7, 7/ SP$nWorkers, 7, 7 / SP$nWorkers) # all 1=0.5, all 4=0.2, 7=0.05
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
      # 6.3. Initialize data storage
      #################
      genoW <- vector("list", length = nColony)
      genoQ <- vector("list", length = nColony)
      phenoColony <- vector("list", length = nColony)
      
      #################
      # 6.2. Calculate ColonyPheno
      #################
      calc_pheno_age_0 <- calcColonyPheno(age_0,queenTrait = c("yieldQueenTrait", "calmQueenTrait"),
                                    workersTrait = c("yieldWorkersTrait", "calmWorkersTrait"),
                                    traitName = c("yield05", "calmness05"),
                                    checkProduction = c(TRUE, FALSE)) |> as.data.frame()

      calc_gv_age_0 <- calcColonyGv(age_0,queenTrait = c("yieldQueenTrait", "calmQueenTrait"),
                                          workersTrait = c("yieldWorkersTrait", "calmWorkersTrait"),
                                          traitName = c("yield05", "calmness05")) |> as.data.frame()
      # Colony heritability
      diag(var(calc_gv_age_0)) / diag(var(calc_pheno_age_0)) # 0.6
      
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
save.image("ahbFounderGenomeHaploAg0_005.RData")

##############################
# 8. Perform GWAS for Queen
##############################
        #################
        # 8.1. Calculate ColonyPheno
        #################
        # To be used for degrees of freedom
        n <- nrow(phenoColony)
        # To be used for degrees of freedom (number of levels in fixed effects)
        k <- 1
   
        #################
        # 8.2. GWAS by RRBLUP approach
        #################             
        Z <- pooledGenoQ[as.character(phenoColony$QID),]
        QgwasRRBLUP <- mmer(yield05~1,
                            random=~vsr(Z),
                            rcov=~units, nIters=3,
                            verbose = FALSE,
                            data=phenoColony)
                    
              ##########
              # 8.2.1. Calculate SNP effect, SE and P values
              ########### 
              a <- QgwasRRBLUP$U$`u:Z`$yield05 # marker effects
              se.a <- sqrt(diag(kronecker(diag(ncol(Z)),QgwasRRBLUP$sigma$`u:Z`)- QgwasRRBLUP$PevU$`u:Z`$yield05)) # SE of
              t.stat <- a/se.a # t-statistic
              pvalRRBLUP <- dt(t.stat,df=n-k-1) # -log10(pval)
               
              ##########
              # 8.2.2. Calculate SNP effect, SE and P values
              ###########
              QgwasRRBLUP <- data.frame(markerName = colnames(pooledGenoQ), 
                                        Beta = a, 
                                        t_stat = t.stat, 
                                        SE = se.a, 
                                        PValue = pvalRRBLUP)
              
              QgwasRRBLUP$index <- seq(1:ncol(pooledGenoQ))
              # Get the SNP chip map and make GWAS summary stat
              snpchipMap <- merge(ahbmngPos, QgwasRRBLUP, by = "markerName")
              snpchipMapOrd <- snpchipMap[order(as.integer(snpchipMap$index)), ]
              snpchipMapOrd$position <- snpchipMapOrd$position
              QgwasRRBLUP_nQTL01_h2_005_nsnpChip1000 <- snpchipMapOrd[, c("markerName", "chromosome", "position", "PValue")]
              
        #################
        # 8.3. Visualize
        #################                 
        CMplot(QgwasRRBLUP_nQTL01_h2_005_nsnpChip1000,type="p",
               plot.type="m",LOG10=TRUE,
               threshold=NULL,file="jpg",
               file.name="",dpi=300,
               file.output=FALSE,
               verbose=TRUE,
               width =14,height=6,chr.labels.angle=0, 
               main = "QgwasRRBLUP")
              
              #################
              # 8.3.1. Add a custom note on the far right
              #################
              mtext(expression("nsnpChip = 1000"), side = 3, line = 3, cex = 0.8, adj = 1, col = "black")
              mtext(expression("nQTL = 1, intialH2=0.05"), side = 3, line = 2, cex = 0.8, adj = 1, col = "black")
              mtext(expression("nWorker = 100"), side = 3, line = 1, cex = 0.8, adj = 1, col = "black")
              
      #################
      # 8.4. GWAS by GBLUP approach
      #################   
      M <- pooledGenoQ
      Qgrm <- calcBeeGRMIbs(x = pooledGenoQ, sex = rep("F", nrow(pooledGenoQ)))
      Qgrminv<-solve(Qgrm + diag(1e-6, ncol(Qgrm), ncol(Qgrm))) ## inverse of MM'
      MTQgrminv<-t(M)%*%Qgrminv
    
              QgwasGBLUP <- mmer(yield05~1,
                                 random=~vsr(QID, Gu=Qgrm),
                                 rcov=~units, nIters=3,
                                 verbose = FALSE,
                                 data=phenoColony)
           ##########
           # 8.4.1. Calculate SNP effect, SE and P values
           ###########
           a.from.g <-MTQgrminv%*%matrix(QgwasGBLUP$U$`u:QID`$yield05,ncol=1)
           var.g <- kronecker(Qgrm,QgwasGBLUP$sigma$`u:QID`)- QgwasGBLUP$PevU$`u:QID`$yield05
           var.a.from.g <- t(M)%*%Qgrminv%*% (var.g) %*% t(Qgrminv)%*%M
           se.a.from.g <- sqrt(diag(var.a.from.g))
           t.stat.from.g <- a.from.g/se.a.from.g
           pvalGBLUP <- dt(t.stat.from.g,df=n-k-1) 
                
                 #######
                 # 8.4.1.1 Calculate SNP effect, SE and P values GBLUP
                 #######
                QgwasGBLUP <- data.frame(markerName = colnames(pooledGenoQ), 
                                         Beta = a.from.g, 
                                         t_stat = t.stat.from.g, 
                                         SE = se.a.from.g, 
                                         PValue = pvalGBLUP)
 
                QgwasGBLUP$index <- seq(1:ncol(pooledGenoQ))
                snpchipMap <- merge(ahbmngPos, QgwasGBLUP, by = "markerName")
                snpchipMapOrd <- snpchipMap[order(as.integer(snpchipMap$index)), ]
                snpchipMapOrd$position <- snpchipMapOrd$position
                QgwasGBLUP_nQTL01_h2_005_nsnpChip1000 <- snpchipMapOrd[, c("markerName", "chromosome", "position", "PValue")]
                
            ##########
            # 8.4.2. Visualize GWAS GBLUP
            ###########

            CMplot(QgwasGBLUP_nQTL01_h2_005_nsnpChip1000,type="p",
                   plot.type="m",LOG10=TRUE,
                   threshold=NULL,file="jpg",
                   file.name="",dpi=300,
                   file.output=FALSE,
                   verbose=TRUE,
                   width =14,height=6,chr.labels.angle=0, main = "QgwasGBLUP")
                
            ########
            # 8.4.3. Add a custom note on the far right
            ########
            mtext(expression("nsnpChip = 1000"), side = 3, line = 3, cex = 0.8, adj = 1, col = "black")
            mtext(expression("nQTL = 1, intialH2=0.05"), side = 3, line = 2, cex = 0.8, adj = 1, col = "black")
            mtext(expression("nWorker = 100"), side = 3, line = 1, cex = 0.8, adj = 1, col = "black")
            ##########
            # 8.4.4. QQ-plot GWAS GBLUP
            ###########
            CMplot(QgwasGBLUP_nQTL01_h2_005_nsnpChip1000,plot.type="q",box=FALSE,file="jpg",file.name="",dpi=300,
                     conf.int=TRUE,conf.int.col=NULL,threshold.col="red",threshold.lty=2,
                     file.output=FALSE,verbose=TRUE,width=5,height=5, main = "QgwasGBLUP")
              
##############################
# 9. Pool worker GWAS
##############################

      ################
      # 9.1 Pull worker genotype
      ################

        for (colony in 1:nColony) {
          pooledGenoW <- getPooledGeno(getSnpGeno(age_0, snpChip = 2, caste = "workers",collapse = T),type = "mean",
                                       sex = getCasteSex(age_0@colonies[[colony]], caste = "workers"))
          genoW[[colony]] <- pooledGenoW
      }
        pooledGenoW <- do.call(rbind, genoW)
        pooledGenoWrnd <- round(pooledGenoW, digits = 0)
        rownames(pooledGenoWrnd) <-phenoColony$WID
        
        # Pre-allocate the matrix for pooledGenoW
        pooledGenoW <- matrix(NA, nrow = length(phenoColony$WID), ncol = ncol(pooledGenoQ)) 
        
        # Loop through colonies
        for (colony in 1:nColony) {
          pooledGenoW[colony, ] <- getPooledGeno(
            getSnpGeno(age_0, snpChip = 2, caste = "workers", collapse = TRUE),
            type = "mean",
            sex = getCasteSex(age_0@colonies[[colony]], caste = "workers")
          )
        }
        
        # Round the pooled genotype
        pooledGenoWrnd <- round(pooledGenoW, digits = 0)
        rownames(pooledGenoWrnd) <- phenoColony$WID
        
##############################
# 10. Save RData
##############################
save.image("ahbFounderGenomeHaploAg0_QW.RData")
        
        
        
        
##############################
# 11. Queen and worker GWAS
##############################  
      
      
      
      
      
      
      
      
