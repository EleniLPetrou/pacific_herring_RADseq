## Script for stacks version 1.46
# Move to your work directory
#cd /mnt/hgfs/D/sequencing_data/Herring_Coastwide_PopulationStructure/
# -P = Directory where all the stacks output files are found (created by pstacks)
# -M = file name where popmap for cstacks is found
# -g = use genomic position indicated from bowtie
# -p = number of threads to use

cstacks -P ./output_stacks -M ./popmaps/cstacks_popmap.txt -g -p 15