# ---- Section 1 - Population structure Q,W and original genotype data ----   
## ---- Section 1.1 - Eigen decomposition of the Qgrm ----
eigen_result <- eigen(Qgrm)
eigenvalues <- eigen_result$values
eigenvectors <- eigen_result$vectors
eigenvalues <- pmax(eigenvalues, 0)

## ---- Section 1.2 - Compute PCA scores ----
pca_scores <- eigenvectors %*% diag(sqrt(eigenvalues))
pca_df <- data.frame(
  PC1 = pca_scores[, 1],
  PC2 = pca_scores[, 2])

## ---- Section 1.3 - Plot the PCA results ----
ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point(size = 2, color = "blue") +
  # geom_text(vjust = -0.5, size = 3) +
  theme_minimal() +
  labs(title = "PCA of Queen Genotype nsnpChip-2",
       x = "PCA 1",
       y = "PCA 2")

# ----- Section 2. population structure -----
selected_genoFile_indices <- colnames(genoQ)
originalGeno <- pullSegSiteGeno(ahbFounderGenomeHaplo)[,selected_genoFile_indices]
sampleInfoQ <- paste(seq(1:nrow(genoQ)),"Q", sep = "_")
row.names(genoQ) <- sampleInfoQ
# worker
sampleInfoW <- paste(seq(1:nrow(poolWorkerChip2_nQTL1_h05_rwnd)),"W", sep = "_")
row.names(poolWorkerChip2_nQTL1_h05_rwnd) <- sampleInfoW

#original
sampleInfoOrg <- paste(seq(1:nrow(originalGeno)),"O", sep = "_")
row.names(originalGeno) <- sampleInfoOrg

Q <- as.data.frame(rep("Q", nrow(genoQ)))
colnames(Q) <- "sampleInfo"
W <- as.data.frame(rep("W", nrow(poolWorkerChip2_nQTL1_h05_rwnd)))
colnames(W) <- "sampleInfo"
initialGeno <- as.data.frame(rep("initialGeno", nrow(originalGeno)))
colnames(initialGeno) <- "sampleInfo"
sampleInfo <- rbind(Q,W,initialGeno)

# Combine all genofile
Q_W_initialGeno <- rbind(genoQ,poolWorkerChip2_nQTL1_h05_rwnd,originalGeno)

## per attribute
subAttribute <- as.matrix(sampleInfo$sampleInfo)
pop.gen_subAttribute <- popgen(M = Q_W_initialGeno, subgroups = subAttribute)
pop.gen_subAttribute$bygroup$F.stats$Genotypes
FstAttribute <- as.data.frame(pop.gen_subAttribute$bygroup$F.stats$Genotypes)
