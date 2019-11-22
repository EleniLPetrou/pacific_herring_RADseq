# Eleni Petrou ran this script on October 17, 2016
# cd /mnt/hgfs/D/sequencing_data/Herring_PopulationStructure/stacks_scripts

#Cherry Point 2016
process_radtags \
-P -p /mnt/hgfs/D/sequencing_data/Herring_PopulationStructure/raw_data/PairedEnd/CherryPoint2016 \
-o /mnt/hgfs/D/sequencing_data/Herring_PopulationStructure/IndividualData_PairedEnd \
-b /mnt/hgfs/D/sequencing_data/Herring_PopulationStructure/barcodes/Barcodes_597_S2_L002.txt \
-i gzfastq \
-y fastq \
-e sbfI \
-r \
-c \
-q \
-E phred33 \
--filter_illumina


#Elliot Bay 2015
process_radtags \
-P -p /mnt/hgfs/D/sequencing_data/Herring_PopulationStructure/raw_data/PairedEnd/ElliotBay2015 \
-o /mnt/hgfs/D/sequencing_data/Herring_PopulationStructure/IndividualData_PairedEnd \
-b /mnt/hgfs/D/sequencing_data/Herring_PopulationStructure/barcodes/Barcodes_598_S3_L003.txt \
-i gzfastq \
-y fastq \
-e sbfI \
-r \
-c \
-q \
-E phred33 \
--filter_illumina