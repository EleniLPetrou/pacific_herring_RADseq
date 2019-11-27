# Calculate population statistics and export several output files.
#cd /mnt/hgfs/D/sequencing_data/Herring_Coastwide_PopulationStructure

# b= batch of data
# P = path to the stacks output files
# M = path to population map
# r = minimum percentage of individuals in a population required to process a locus for that population
# p = minimum number of populations a locus must be present to process that locus
# m = minimum stack depth required for individuals at a locus
# W= whitelist of loci to be used in analysis

populations -b 1 -P ./output_stacks -M ./output_stacks_populations/filtered_haplotypesANDsnps_1104indiv_7261loci/genome_alignment/whitelist_indiv.txt -m 10 -W ./output_stacks_populations/filtered_haplotypesANDsnps_1104indiv_7261loci/genome_alignment/whitelist_loci.txt -O ./output_stacks_populations/filtered_haplotypesANDsnps_1104indiv_7261loci/genome_alignment --min_maf 0.05 --vcf --genepop
