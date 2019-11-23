# pacific_herring_RADseq

This repository contains the bash and R scripts that I used to analyze RAD sequencing data collected from wild spawning populations of Pacific herring. 

In a nutshell, we collected tissue samples from  ~1,300 spawning adult herring across the Pacific Northwest Coast of North America. That was an odyssey. Then, we sequenced each sample on an Illumina 2500 or 4000, using restriction site-associated DNA (RAD) sequencing, and following the protocol of Etter et al. 2011. 

Using the *stacks* pipeline, we assembled sequencing reads into loci and called genotypes. 
After filtering the data, we conducted analyses on the spatial and temporal population structure of Pacific herring. These were mainly conducted using  R packages, *bioconductor*,  and the fantastic *tidyverse*. 


