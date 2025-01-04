# ---- Section 1 - Load library ----
library(AlphaSimR)
library(SIMplyBee)
library(data.table)
library(parallel)
library(matrixcalc)
library(AGHmatrix)
library(CMplot)

# ---- Section 2 - Set workspace ----
setwd("/Users/mulusewfikere/Documents/managed AHB/data/gregor/")

# ---- Section 3 - Read genotype and map data ----

## ---- Section 3.1 - Read mapFile ----
ahbmngPos <- as.data.frame(fread("rawData/haplotypeMapFile.pos", h=F))
colnames(ahbmngPos) <- c("chromosome", "position")
ahbmngPos$test <- sequence(rle(as.character(ahbmngPos$chromosome))$lengths)
ahbmngPos <- data.frame(markerName=paste(ahbmngPos$chromosome,ahbmngPos$test, sep="_"), ahbmngPos)
ahbmngPos$test <- NULL
ahbmngPos$position <- ahbmngPos$position/1000000 
ahbmngPos$chromosome <- as.numeric(ahbmngPos$chromosome)

## ---- Section 3.2 - Read genotype File ----
ahbmngGeno <- t(as.data.frame(fread("rawData/haplotypeGenoFile", h=F)))
colnames(ahbmngGeno) <- ahbmngPos$markerName
rownames(ahbmngGeno) <- seq(1:nrow(ahbmngGeno))

# ---- Section 4 - Create a founder genome ----

## ---- Section 4.1 - Import haplotype after phasing VCF fil ----

ahbFounderGenomeHaplo = importHaplo(haplo=ahbmngGeno,
                                    genMap=ahbmngPos,
                                    ploidy = 2L,
                                    ped=NULL)

### ---- Section 4.1.1 - Check Number of samples and SNP ----
ahbFounderGenomeHaplo@nInd
ahbFounderGenomeHaplo@nLoci
ahbFounderGenomeHaplo@nChr

# ---- Section 5 - Initialize simulation ----
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

## ---- Section 5.1 - SP global simulation ----

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

## ---- Section 5.2 - Pre-define variances ----
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

# ---- Section 6 - Assign SNPchip and extract QTL ----

SP$addSnpChip(nChip100)
SP$addSnpChip(nChip1000)
SP$addSnpChip(nChip10000)
qtl_1 <- colnames(pullQtlGeno(ahbFounderGenomeHaplo, simParam=SP))

# ---- Section 7 - Create virgin queen and father ----

virginQueens <- createVirginQueens(x = ahbFounderGenomeHaplo, n = nMelN)
drones = createDrones(virginQueens, nInd = nDronesPerQueen)

## ---- Section 7.1 - Make a cross ----

colonies <- cross(x = virginQueens,
                  drones = drones,
                  crossPlan = "create",
                  spatial = FALSE,
                  nDrones = nFathersPoisson,
                  checkCross = "warning")

## ---- Section 7.2 - Build ----

age_0 <- createMultiColony(x = colonies)
age_0 <- buildUp(age_0)
nColony <- nColonies(age_0)

## ---- Section 7.3 - Initialize data storage ----

genoW <- vector("list", length = nColony)
genoQ <- vector("list", length = nColony)
phenoColony <- vector("list", length = nColony)

## ---- Section 7.4 - Calculate ColonyPhen ----

calc_pheno_age_0 <- calcColonyPheno(age_0,queenTrait = c("yieldQueenTrait", "calmQueenTrait"),
                                    workersTrait = c("yieldWorkersTrait", "calmWorkersTrait"),
                                    traitName = c("yield05", "calmness05"),
                                    checkProduction = c(TRUE, FALSE)) |> as.data.frame()

calc_gv_age_0 <- calcColonyGv(age_0,queenTrait = c("yieldQueenTrait", "calmQueenTrait"),
                              workersTrait = c("yieldWorkersTrait", "calmWorkersTrait"),
                              traitName = c("yield05", "calmness05")) |> as.data.frame()

## ---- Section 7.5 - Colony heritability ----      

diag(var(calc_gv_age_0)) / diag(var(calc_pheno_age_0)) # 0.6

## ---- Section 7.6 - Structure phenotype data ----        

phenoColony <- calc_pheno_age_0
phenoColony$QID <- seq(1:nrow(phenoColony))
phenoColony$WID <- seq(nrow(phenoColony) + 1, length.out = nrow(phenoColony))

## ---- Section 7.7 - Structure queenGenotype data ----      

genoQ <- getSnpGeno(age_0, caste = "queen", snpChip = 1, collapse = TRUE)
rownames(genoQ) <-phenoColony$QID
phenoColony$ID <- factor(phenoColony$QID, levels = rownames(genoQ))
colnames(phenoColony)[1] <- "yield05"
phenoColony$QID <- as.factor(phenoColony$QID)
phenoColony$WID <- as.factor(phenoColony$WID)  

# ---- Section 8 - Save RData ---- 

save.image("ahbFounderGenomeHaploAg0_005.RData")

# ---- Section 9 - Perform GWAS for Queen ---- 

# To be used for degrees of freedom
n <- nrow(phenoColony)
# To be used for degrees of freedom (number of levels in fixed effects)
k <- 1


## ---- Section 9.1 - Center M matrix and make GRM ---- 
genoQ_centered <- scale(genoQ, center = TRUE, scale = FALSE)
# Qgrm <- calcBeeGRMIbs(x = genoQ, sex = rep("F", nrow(genoQ)))
Qgrm <- Gmatrix(SNPmatrix=genoQ, 
                missingValue=-9, 
                maf=0.05, 
                method="VanRaden")

Qgrminv <-solve(Qgrm + diag(1e-6, ncol(Qgrm), ncol(Qgrm)))
MTQgrminv<- t(genoQ_centered)%*%Qgrminv

## ---- Section 9.2 - GWAS by GBLUP approach: Fit a model ----

QgwasGBLUP <- mmer(yield05~1,
                   random=~vsr(QID, Gu=Qgrm),
                   rcov=~units, nIters=3,
                   verbose = FALSE,
                   data=phenoColony)

## ---- Section 9.3 - Calculate SNP effect, SE and P values ----

a.from.g <-MTQgrminv%*%matrix(QgwasGBLUP$U$`u:QID`$yield05,ncol=1)
plot(genoQ %*% a.from.g, QgwasGBLUP$U$`u:QID`$yield05); abline(a=0,b=1)
var.g <- kronecker(Qgrm,QgwasGBLUP$sigma$`u:QID`)- QgwasGBLUP$PevU$`u:QID`$yield05
var.a.from.g <- t(genoQ_centered)%*%Qgrminv%*% (var.g) %*% t(Qgrminv)%*%genoQ_centered
se.a.from.g <- sqrt(diag(var.a.from.g))
t.stat.from.g <- a.from.g/se.a.from.g
pvalGBLUP <- dt(t.stat.from.g,df=n-k-1) 

## ---- Section 9.4 - Calculate SNP effect, SE and P values GBLUP ----

QgwasGBLUP <- data.frame(markerName = colnames(genoQ), 
                         Beta = a.from.g, 
                         t_stat = t.stat.from.g, 
                         SE = se.a.from.g, 
                         PValue = pvalGBLUP)
 
## ---- Section 9.5 - Collate GWAS summary stat ----

QgwasGBLUP$index <- seq(1:ncol(genoQ))
snpchipMap <- merge(ahbmngPos, QgwasGBLUP, by = "markerName")
snpchipMapOrd <- snpchipMap[order(as.integer(snpchipMap$index)), ]
snpchipMapOrd$position <- snpchipMapOrd$position
QgwasGBLUP_nQTL01_h2_005_nsnpChip100 <- snpchipMapOrd[, c("markerName", "chromosome", "position", "PValue")]
                
# All nChip scenario
# QgwasGBLUP_nQTL01_h2_005_nsnpChip100 <- snpchipMapOrdSelGBLUP
# QgwasGBLUP_nQTL01_h2_005_nsnpChip1000 <- snpchipMapOrdSelGBLUP
# QgwasGBLUP_nQTL01_h2_005_nsnpChip10000 <- snpchipMapOrdSelGBLUP

## ---- Section 9.5 - Visualize GWAS GBLUP ----

CMplot(QgwasGBLUP_nQTL01_h2_005_nsnpChip1000,type="p",
       plot.type="m",LOG10=TRUE,
       threshold=NULL,file="jpg",
       file.name="",dpi=300,
       file.output=FALSE,
       verbose=TRUE,
       width =14,height=6,
       chr.labels.angle=0, 
       main = "QgwasGBLUP")
 
### ---- Section 9.5.1 - Add a custom note on the far right ----               
mtext(expression("nsnpChip = 1000"), 
      side = 3, 
      line = 3, 
      cex = 0.8, 
      adj = 1, 
      col = "black")
mtext(expression("nQTL = 1, intialH2=0.05"), 
      side = 3, 
      line = 2, 
      cex = 0.8, 
      adj = 1, 
      col = "black")
mtext(expression("nWorker = 100"), 
      side = 3, 
      line = 1, 
      cex = 0.8, 
      adj = 1, 
      col = "black")
 
## ---- Section 9.6 - QQ-plot GWAS GBLUP ----          

CMplot(QgwasGBLUP_nQTL01_h2_005_nsnpChip1000,
       plot.type="q",
       box=FALSE,
       file="jpg",
       file.name="",
       dpi=300,
       conf.int=TRUE,
       conf.int.col=NULL,
       threshold.col="red",
       threshold.lty=2,
       file.output=FALSE,
       verbose=TRUE,
       width=5,
       height=5, 
       main = "QgwasGBLUP")

# ---- Section 10 - Pool worker GWAS ----     

## ---- Section 10.1 - Pull worker genotype ----    

for (colony in 1:nColony) {
pooledGenoW <- getPooledGeno(getSnpGeno(age_0, snpChip = 1, 
                             caste = "workers",
                             collapse = T),
                             type = "mean",
                             sex = getCasteSex(age_0@colonies[[colony]], 
                             caste = "workers"))
              genoW[[colony]] <- pooledGenoW
      }

getChipsName <- colnames(getSnpGeno(age_0, snpChip = 1, caste = "workers",collapse = T))
poolWorkerChips_nQTL1_h05 <- do.call(rbind, lapply(X = getSegSiteGeno(x = age_0, caste = "workers", nInd = 1), 
                                         FUN = getPooledGeno, type = "mean"))[,getChipsName]

## ---- Section 10.2 - Collate Pull worker genotype ----    
pooledGenoW <- do.call(rbind, genoW)
pooledGenoWrnd <- round(pooledGenoW, digits = 0)
rownames(pooledGenoWrnd) <-phenoColony$WID

# ---- Section 11 - Fit worker model ----     

## ---- Section 11.1 - Center M matrix and make GRM ---- 
genoPoolWorker_centered <- scale(genoPoolWorker, center = TRUE, scale = FALSE)
Wgrm <- Gmatrix(SNPmatrix=genoPoolWorker, 
                missingValue=-9, 
                maf=0.05, 
                method="VanRaden")

Wgrminv <-solve(Wgrm + diag(1e-6, ncol(Wgrm), ncol(Wgrm)))
MTWgrminv<- t(genoPoolWorker_centered)%*%Wgrminv

## ---- Section 11.2 - GWAS by GBLUP approach: Fit a model ----

WgwasGBLUP <- mmer(yield05~1,
                   random=~vsr(WID, Gu=Wgrm),
                   rcov=~units, nIters=3,
                   verbose = FALSE,
                   data=phenoColony)

## ---- Section 11.3 - Calculate SNP effect, SE and P values ----

a.from.g <-MTWgrminv%*%matrix(WgwasGBLUP$U$`u:WID`$yield05,ncol=1)
plot(genoPoolWorker %*% a.from.g, WgwasGBLUP$U$`u:WID`$yield05); abline(a=0,b=1)
var.g <- kronecker(Wgrm,WgwasGBLUP$sigma$`u:WID`)- WgwasGBLUP$PevU$`u:WID`$yield05
var.a.from.g <- t(genoPoolWorker_centered)%*%Wgrminv%*% (var.g) %*% t(Wgrminv)%*%genoPoolWorker_centered
se.a.from.g <- sWrt(diag(var.a.from.g))
t.stat.from.g <- a.from.g/se.a.from.g
pvalGBLUP <- dt(t.stat.from.g,df=n-k-1) 

## ---- Section 11.4 - Calculate SNP effect, SE and P values GBLUP ----

WgwasGBLUP <- data.frame(markerName = colnames(genoPoolWorker), 
                         Beta = a.from.g, 
                         t_stat = t.stat.from.g, 
                         SE = se.a.from.g, 
                         PValue = pvalGBLUP)

## ---- Section 11.5 - Collate GWAS summary stat ----

WgwasGBLUP$index <- seW(1:ncol(genoPoolWorker))
snpchipMap <- merge(ahbmngPos, WgwasGBLUP, by = "markerName")
snpchipMapOrd <- snpchipMap[order(as.integer(snpchipMap$index)), ]
snpchipMapOrd$position <- snpchipMapOrd$position
WgwasGBLUP_nWTL01_h2_005_nsnpChip100 <- snpchipMapOrd[, c("markerName", "chromosome", "position", "PValue")]

## ---- Section 11.5 - Visualize GWAS GBLUP ----

CMplot(WgwasGBLUP_nWTL01_h2_005_nsnpChip100,type="p",
       plot.type="m",LOG10=TRUE,
       threshold=NULL,file="jpg",
       file.name="",dpi=300,
       file.output=FALSE,
       verbose=TRUE,
       width =14,height=6,
       chr.labels.angle=0, 
       main = "WgwasGBLUP")

### ---- Section 11.5.1 - Add a custom note on the far right ----               
mtext(expression("nsnpChip = 1000"), 
      side = 3, 
      line = 3, 
      cex = 0.8, 
      adj = 1, 
      col = "black")
mtext(expression("nWTL = 1, intialH2=0.05"), 
      side = 3, 
      line = 2, 
      cex = 0.8, 
      adj = 1, 
      col = "black")
mtext(expression("nWorker = 100"), 
      side = 3, 
      line = 1, 
      cex = 0.8, 
      adj = 1, 
      col = "black")

## ---- Section 11.6 - WW-plot GWAS GBLUP ----          

CMplot(WgwasGBLUP_nWTL01_h2_005_nsnpChip100,
       plot.type="W",
       box=FALSE,
       file="jpg",
       file.name="",
       dpi=300,
       conf.int=TRUE,
       conf.int.col=NULL,
       threshold.col="red",
       threshold.lty=2,
       file.output=FALSE,
       verbose=TRUE,
       width=5,
       height=5, 
       main = "WgwasGBLUP")      
        
# ---- Section 12 - Queen + worker GWAS Model ----           
## ---- Section 12.1. Combine the raw genotype ----
rownames(genoPoolWorker) <- phenoColony$WID
rownames(genoQ) <- phenoColony$QID
genoQW <- rbind(genoQ,genoPoolWorker)

## ---- Section 12.2 - Construct GRM ---- 
QWgrm <- Gmatrix(SNPmatrix=genoQW, 
                missingValue=-9, 
                maf=0.05, 
                method="VanRaden")

QWgrm.2 <- QWgrm + diag(nrow(genoQW))*0.000001
genoQW_centered <- scale(QWgrm.2, center = TRUE, scale = FALSE)

QWgrminv <-solve(QWgrm.2 + diag(1e-6, ncol(QWgrm.2), ncol(QWgrm.2)))
MTQWgrminv <- t(genoPoolWorker_centered)%*%QWgrminv

## ---- Section 12.3 - GWAS by GBLUP approach: Fit a model ----
phenoColony$QID <- factor(phenoColony$QID, levels = rownames(QWgrm.2))
phenoColony$WID <- factor(phenoColony$WID, levels = rownames(QWgrm.2))

asreml.options(dense = ~ vm(QID, QWgrm.2))
asreml.options(dense = ~ vm(WID, QWgrm.2))
QWmodel.asr <- asreml(fixed = yield05 ~ 1,
                      random = ~ vm(QID, QWgrm.2) + vm(WID, QWgrm.2),
                      data = phenoColony,
                      maxiter = 100,
                      # ai.loadings = TRUE,
                      # ai.sing = TRUE,
                      workspace = "32gb")
    
## ---- Section 12.4. QUEEN: Predict SNP effects using the model and Extract variance components ----
QW_Q_predictions <- predict(QWmodel.asr, classify = "QID", only = "vm(QID, Qgrm.2)", pworkspace = "16gb", vcov = T)
QW_Q_blup_individuals <- QW_Q_predictions$pvals
## ---- Section 12.5. WORKER: Predict SNP effects using the model and Extract variance components ----
QW_W_predictions <- predict(QWmodel.asr, classify = "WID", only = "vm(WID, Qgrm.2)", pworkspace = "16gb", vcov = T)
QW_W_blup_individuals <- QW_W_predictions$pvals
QW_W_pev_individuals <- as.matrix(predictions$vcov)

## ---- Section 12.6. Backsolve SNP effects Q + W model ----


