# pacific_herring_RADseq

This repository contains the bash and R scripts used to analyze RAD sequencing data collected from wild spawning populations of Pacific herring. 

In brief, we collected tissue samples from  ~1,300 spawning adult herring from spawning sites across the Pacific Northwest Coast of North America.  Then, we sequenced each sample on an Illumina 2500 or 4000, using restriction site-associated DNA (RAD) sequencing.

Using the *stacks* pipeline, we assembled sequencing reads into loci and called genotypes. 
After filtering the data using *vcftools* , we conducted analyses on the spatial and temporal population structure of Pacific herring. These were mainly conducted using  R packages such as *bioconductor* and *tidyverse*. 

