#!/bin/bash
#First line tells ubuntu what interpreter to use
#To run this script in ubuntu terminal, supply the terminal with the correct directory using cd command. Then, write ./bowtie2_alignment.sh


#Explanation of terms:
#bowtie2 -q -x <bt2-idx> -U <r> -S <sam>
#-x <bt2-idx> Indexed "reference genome" filename prefix (minus trailing .X.bt2).
#-U <r> Files with unpaired reads.
#-S <sam> File for SAM output (default: stdout)
#-f query input files are in fasta format
# you can use 2> to redirect stdout to a file (bowtie writes the summary log files to stdout)


# Set the source directory, where all of the files and the script is stored

src=/mnt/hgfs/D/sequencing_data/Herring_Coastwide_PopulationStructure/output_stacks_populations/filtered_haplotypesANDsnps_1104indiv_7261loci/genome_alignment


# Give the basename of the files that you want to align

files="
catalog_loci_7261
"

for file in $files
do
    bowtie2 -f -x Atlantic_herring -U ${file}.fasta -S ${file}.sam 2> ${file}_bowtie.log
done
