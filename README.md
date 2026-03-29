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


## Biological Context

Honey bee colonies consist of:

- Queens (reproductive female, diploid)
- Workers (sterile females, colony contributors)
- Drones (haploid males)

Traits are influenced by both individual and colony-level genetics, making GWAS analysis more complex than in standard diploid populations.


## ⚙️ Workflow Summary

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

## GWAS Models

Three GWAS models are implemented:

### Queen-only Model
Uses queen genotype with genomic relationship matrix (GRM):

```r
random = ~ vs(id, Gu = QGRM)

### Workers-only Model
Uses workers genotype with genomic relationship matrix (GRM):

```r
random = ~ vs(id, Gu = WGRM)

### Queen + Workers Model
Uses queen + workers genotype with genomic relationship matrix (GRM):

```r
random=~vsr(id, Gu=QGRM) + vsr(id, Gu=WGRM)
