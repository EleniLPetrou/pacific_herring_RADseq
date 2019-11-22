#!/bin/bash
#First line tells ubuntu what interpreter to use
#To run this script in ubuntu terminal, supply the terminal with the correct directory using cd command. Then, write ./bowtie_bash_batch6.sh
#cd /mnt/hgfs/D/sequencing_data/Herring_Coastwide_PopulationStructure/data_IndividualData_Trimmed_90bp
src=/mnt/hgfs/D/sequencing_data/Herring_Coastwide_PopulationStructure/data_IndividualData_Trimmed_90bp

files="
Venn99_001
Venn99_009
Venn99_017
Venn99_025
Venn99_033
Venn99_041
Venn99_002
Venn99_010
Venn99_018
Venn99_026
Venn99_034
Venn99_042
Venn99_003
Venn99_011
Venn99_019
Venn99_027
Venn99_035
Venn99_043
Venn99_004
Venn99_012
Venn99_020
Venn99_028
Venn99_036
Venn99_044
Venn99_005
Venn99_013
Venn99_021
Venn99_029
Venn99_037
Venn99_045
Venn99_006
Venn99_014
Venn99_022
Venn99_030
Venn99_038
Venn99_046
Venn99_007
Venn99_015
Venn99_023
Venn99_031
Venn99_039
Venn99_047
Venn99_008
Venn99_016
Venn99_024
Venn99_032
Venn99_040
Skid99_001
Skid99_009
Skid99_017
Skid99_025
Skid99_033
Skid99_041
Skid99_002
Skid99_010
Skid99_018
Skid99_026
Skid99_034
Skid99_042
Skid99_003
Skid99_011
Skid99_019
Skid99_027
Skid99_035
Skid99_043
Skid99_004
Skid99_012
Skid99_020
Skid99_028
Skid99_036
Skid99_044
Skid99_005
Skid99_013
Skid99_021
Skid99_029
Skid99_037
Skid99_045
Skid99_006
Skid99_014
Skid99_022
Skid99_030
Skid99_038
Skid99_046
Skid99_007
Skid99_014_REP
Skid99_023
Skid99_031
Skid99_039
Skid99_047
Skid99_008
Skid99_016
Skid99_024
Skid99_032
Skid99_040
Skid99_048
"

# Run bowtie on the files; the i variable will be our ID for each sample we process.
i=1

for file in $files
do 
    bowtie -q -v 3 --norc --sam batch_5_ReferenceGenome /mnt/hgfs/D/sequencing_data/Herring_Coastwide_PopulationStructure/data_IndividualData_Trimmed_90bp/${file}.fq ${file}.sam
    let "i+=1";
done

