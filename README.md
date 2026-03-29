# BeeGWAS

Genome-wide association study (GWAS) simulation and analysis pipeline for honey bee (*Apis mellifera*) colony traits.

This repository contains the full R workflow used to simulate colony-based genetic architectures and evaluate GWAS performance under different biological and statistical scenarios. 

The associated manuscript has been submitted to *Scientific Reports*.

## Overview

BeeGWAS is a simulation and analysis framework designed to study genome-wide associations in honey bees using realistic colony structure. The pipeline integrates queen, drone, and worker genetic contributions to model complex colony-level phenotypes.

It is designed to evaluate:

- GWAS power under different genetic architectures  
- False discovery rate (FDR) across simulation scenarios  
- Impact of sample size, heritability, and SNP density  
- Contribution of queen vs worker genomes  
- Joint modeling of colony genetic effects  

### Please review the workflow figure for the detail.

![BeeGWAS Workflow](https://raw.githubusercontent.com/DrMulusewFikere/BeeGWAS/main/workflowBeeGWAS.png)

## Biological Context

Honey bee colonies consist of:

- Queens (reproductive female, diploid)
- Workers (sterile females, colony contributors)
- Drones (haploid males)

Traits are influenced by both individual and colony-level genetics, making GWAS analysis more complex than in standard diploid populations.


## Workflow Summary

The pipeline includes the following steps:

### 1. Data Input
- Haplotyped genotype matrix
- Marker position file

### 2. Simulation Setup
- Define sample size
- Set heritability levels
- Configure number of QTLs
- Select SNP chip density

### 3. Colony Simulation
Using `SIMplyBee`:
- Generate virgin queens
- Produce drones
- Simulate mating
- Construct colonies
- Generate phenotypes

### 4. Trait Simulation
Traits include:
- Yield (production-related)
- Calmness (behavioral trait)

Phenotypes include:
- Queen genetic effects
- Worker genetic effects
- Environmental variation

## Fitting Models, GWAS statistical Power and FDR

Three GWAS models are implemented:

### Queen-only Model
Uses queen genotype with genomic relationship matrix (GRM):

```r
random = ~ vs(id, Gu = QGRM)
```
### Queen-only Model
Uses workers genotype with genomic relationship matrix (GRM):

```r
random = ~ vs(id, Gu = WGRM)
```
### Queen + Workers Model
Uses queen + workers genotype with genomic relationship matrix (GRM):

```r
random=~vsr(id, Gu=QGRM) + vsr(id, Gu=WGRM)
```

# GWAS Power and False Discovery Rate (FDR) Analysis

This section describes the evaluation of GWAS performance in terms of **statistical power** and **false discovery rate (FDR)** across different simulation scenarios for queen, worker, and combined colony models.

The analysis was performed across replicated simulations using predefined QTL positions from the founder genome.


## Overview

For each scenario, GWAS performance was evaluated using:

- **True Positive (TP)** detection of known simulated QTLs  
- **False Positive (FP)** detection of non-QTL signals  
- **False Discovery Rate (FDR)** estimation  
- **Statistical power** of SNP detection  
- **P-value distribution (Pmin analysis)** for threshold calibration  


# FDR Estimation

## Definition of significant SNPs

Genome-wide significant SNPs were defined using:

\[
P < 1 \times 10^{-6}
\]

A genomic window of:

\[
\pm 5,000 \text{ bp}
\]

was used to match significant SNPs with predefined causal QTL positions.

## True Positive (TP) definition

A QTL was considered detected if at least one significant SNP fell within the ±5 kb window of its position:

- **TP = detected predefined QTLs + overlapping significant SNP signals**


## False Positive (FP) definition

False positives were defined as significant SNPs that do not overlap any predefined QTL window:

- **FP = significant SNPs not matching known QTL regions**


## False Discovery Rate (FDR)

FDR was calculated as:

\[
FDR = \frac{FP}{TP + FP}
\]

Each simulation replicate contributes one FDR estimate per scenario.

## Output structure

For each scenario (sample size, heritability, SNP density, QTL number, replication), results are stored as:

| Scenario | TP | FP | FDR |
|----------|----|----|-----|
| simulated_condition | count | count | value |


#  GWAS Power Analysis

## SNP-wise power calculation

GWAS power was computed for each SNP using the **non-centrality parameter (NCP)**:

NCP = (β / SE) × √[ N × 2 × MAF × (1 − MAF) ]

where:

- \( beta \) = SNP effect size  
- \( SE \) = standard error  
- \( N \) = sample size  
- \( MAF \) = minor allele frequency  

## Power computation

Statistical power was calculated as:

\[
Power = P(Z > Z_{\alpha/2} | NCP)
\]

implemented as:

- Two-sided test at:
\[
\alpha = 10^{-6}
\]

- Using normal approximation of non-centrality distribution

## Scenario-level power

For each simulation replicate:

- SNP-level power was computed
- Mean GWAS power was calculated across all SNPs
- Standard error of power was estimated:

\[
SE = \frac{SD}{\sqrt{n}}
\]


## Output structure

Each scenario produces:

| Scenario | Power | SE |
|----------|-------|----|
| simulated_condition | mean power | standard error |


# Pmin Distribution and Effective Number of Tests

To characterize the empirical significance threshold, an intermediate simulation scenario was used:

- \( h^2 = 0.2 \)
- \( nQTL = 10 \)
- \( SNP density = 10k per chromosome \)


## Procedure

1. Transform p-values:

\[
-\log_{10}(P)
\]

2. Estimate empirical 95th percentile using **Harrell–Davis estimator**

3. Compute:

- Effective significance threshold:
\[
P_{sig}
\]

- Effective number of independent markers:
\[
N_{eff} = \frac{\alpha}{P_{sig}}
\]

- Marker independence ratio:
\[
\text{Ratio} = \frac{N_{eff}}{Total\ markers}
\]

---

## Output

| Metric | Description |
|--------|------------|
| Psig | Empirical significance threshold |
| CI | Bootstrap confidence interval |
| Neff | Effective number of independent variants |
| Ratio | Proportion of independent markers |

# Summary

This framework enables systematic evaluation of:

- GWAS detection power across biological scenarios  
- Accuracy of causal QTL recovery  
- False discovery behavior under complex colony genetics  
- Empirical calibration of genome-wide significance thresholds  

The pipeline is applied to **queen, worker, and combined colony genetic models** in honey bees.
