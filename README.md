# Analysis of Pacific herring population structure

![herring-img](https://github.com/EleniLPetrou/pacific_herring_RADseq/blob/master/images/herring.jpg?raw=true)

This repository contains scripts used to analyze genomic data collected from wild spawning populations of Pacific herring.
These data were produced using restriction site-associated DNA (RAD) sequencing. 

## Study Overview

In brief, we collected tissue samples from  ~1,300 spawning adult herring from spawning sites across the Pacific Northwest Coast of North America.  Then, we sequenced each sample on an Illumina 2500 or 4000, using RAD sequencing.

Using the *stacks* pipeline, we assembled sequencing reads into loci and called genotypes. 
After filtering the data using *vcftools* , we conducted analyses on the spatial and temporal population structure of Pacific herring. These were mainly conducted using  R packages such as *bioconductor* and *tidyverse*. 

## Directory structure

Scripts are organized into folders, whose name describes the analysis conducted:
- R
-bayenv2
-bowtie2
-stacks
-vcftools




