#cd /mnt/hgfs/D/sequencing_data/Herring_PopulationStructure/data_raw

######Paired-end samples

###################################################
#Cherry Point 2016 and Elliot Bay 2015_A
process_radtags \
-P -p /mnt/hgfs/D/sequencing_data/Herring_PopulationStructure/data_raw/PairedEnd/CherryPoint2016_ElliotBay2015_FirstLane \
-o /mnt/hgfs/D/sequencing_data/Herring_PopulationStructure/data_IndividualData_Trimmed_90bp \
-b /mnt/hgfs/D/sequencing_data/Herring_PopulationStructure/data_barcodes/Barcodes_597_S2_L002.txt \
-i gzfastq \
-y fastq \
-e sbfI \
-r \
-c \
-q \
-t 90 \
-E phred33 \
--filter_illumina


#Cherry Point 2016 and Elliot Bay 2015_B
process_radtags \
-P -p /mnt/hgfs/D/sequencing_data/Herring_PopulationStructure/data_raw/PairedEnd/CherryPoint2016_ElliotBay2015_SecondLane \
-o /mnt/hgfs/D/sequencing_data/Herring_PopulationStructure/data_IndividualData_Trimmed_90bp \
-b /mnt/hgfs/D/sequencing_data/Herring_PopulationStructure/data_barcodes/Barcodes_598_S3_L003.txt \
-i gzfastq \
-y fastq \
-e sbfI \
-r \
-c \
-q \
-t 90 \
-E phred33 \
--filter_illumina



##### Single-end samples


###########################################
#Cherry Point 2014
process_radtags \
-p /mnt/hgfs/D/sequencing_data/Herring_PopulationStructure/data_raw/SingleEnd/CherryPoint2014 \
-o /mnt/hgfs/D/sequencing_data/Herring_PopulationStructure/data_IndividualData_Trimmed_90bp \
-b /mnt/hgfs/D/sequencing_data/Herring_PopulationStructure/data_barcodes/Barcodes_490_S102_L007_R1_001.txt \
-i gzfastq \
-y fastq \
-e sbfI \
-r \
-c \
-q \
-t 90 \
-E phred33 \
--filter_illumina

#Port Gamble 2014
process_radtags \
-p /mnt/hgfs/D/sequencing_data/Herring_PopulationStructure/data_raw/SingleEnd/PortGamble2014 \
-o /mnt/hgfs/D/sequencing_data/Herring_PopulationStructure/data_IndividualData_Trimmed_90bp \
-b /mnt/hgfs/D/sequencing_data/Herring_PopulationStructure/data_barcodes/Barcodes_487_S99_L004_R1_001.txt \
-i gzfastq \
-y fastq \
-e sbfI \
-r \
-c \
-q \
-t 90 \
-E phred33 \
--filter_illumina


#Port Orchard 2014
process_radtags \
-p /mnt/hgfs/D/sequencing_data/Herring_PopulationStructure/data_raw/SingleEnd/PortOrchard2014 \
-o /mnt/hgfs/D/sequencing_data/Herring_PopulationStructure/data_IndividualData_Trimmed_90bp \
-b /mnt/hgfs/D/sequencing_data/Herring_PopulationStructure/data_barcodes/Barcodes_491_S103_L008_R1_001_V2.txt \
-i gzfastq \
-y fastq \
-e sbfI \
-r \
-c \
-q \
-t 90 \
-E phred33 \
--filter_illumina


#Quilcene Bay 2014
process_radtags \
-p /mnt/hgfs/D/sequencing_data/Herring_PopulationStructure/data_raw/SingleEnd/QuilceneBay2014 \
-o /mnt/hgfs/D/sequencing_data/Herring_PopulationStructure/data_IndividualData_Trimmed_90bp \
-b /mnt/hgfs/D/sequencing_data/Herring_PopulationStructure/data_barcodes/Barcodes_488_S100_L005_R1_001.txt \
-i gzfastq \
-y fastq \
-e sbfI \
-r \
-c \
-q \
-t 90 \
-E phred33 \
--filter_illumina


#Similk Bay 2015
process_radtags \
-p /mnt/hgfs/D/sequencing_data/Herring_PopulationStructure/data_raw/SingleEnd/SimilkBay2015 \
-o /mnt/hgfs/D/sequencing_data/Herring_PopulationStructure/data_IndividualData_Trimmed_90bp \
-b /mnt/hgfs/D/sequencing_data/Herring_PopulationStructure/data_barcodes/Barcodes_489_S101_L006_R1_001.txt \
-i gzfastq \
-y fastq \
-e sbfI \
-r \
-c \
-q \
-t 90 \
-E phred33 \
--filter_illumina

#Squaxin Pass 2014
process_radtags \
-p /mnt/hgfs/D/sequencing_data/Herring_PopulationStructure/data_raw/SingleEnd/SquaxinPass2014 \
-o /mnt/hgfs/D/sequencing_data/Herring_PopulationStructure/data_IndividualData_Trimmed_90bp \
-b /mnt/hgfs/D/sequencing_data/Herring_PopulationStructure/data_barcodes/Barcodes_492_S1_L001_R1_001_V2.txt \
-i gzfastq \
-y fastq \
-e sbfI \
-r \
-c \
-q \
-t 90 \
-E phred33 \
--filter_illumina

#Spiller Channel 2015 and Case Inlet 2007
#Raw data filename : 1090_S3_L003_R1_001.fastq.gz
process_radtags \
-p /mnt/hgfs/D/sequencing_data/Herring_PopulationStructure/data_raw/SingleEnd/SpillerChannel2015_CaseInlet2007 \
-o /mnt/hgfs/D/sequencing_data/Herring_PopulationStructure/data_IndividualData_Trimmed_90bp \
-b /mnt/hgfs/D/sequencing_data/Herring_PopulationStructure/data_barcodes/Barcodes_1090_S3_L003_R1_001.txt \
-i gzfastq \
-y fastq \
-e sbfI \
-r \
-c \
-q \
-t 90 \
-E phred33 \
--filter_illumina

#Harriet Island 2014 and Kwakume 2015
#Raw data filename : 1089_S2_L002_R1_001.fastq.gz
process_radtags \
-p /mnt/hgfs/D/sequencing_data/Herring_PopulationStructure/data_raw/SingleEnd/Harriet2014_Kwakume2015 \
-o /mnt/hgfs/D/sequencing_data/Herring_PopulationStructure/data_IndividualData_Trimmed_90bp \
-b /mnt/hgfs/D/sequencing_data/Herring_PopulationStructure/data_barcodes/Barcodes_1089_S2_L002_R1_001.txt \
-i gzfastq \
-y fastq \
-e sbfI \
-r \
-c \
-q \
-t 90 \
-E phred33 \
--filter_illumina

#Salisbury 2016 and Berner's Bay 2016
#Raw data filename : 1091_S4_L004_R1_001.fastq.gz 
process_radtags \
-p /mnt/hgfs/D/sequencing_data/Herring_PopulationStructure/data_raw/SingleEnd/BernersBay2016_Salisbury2016 \
-o /mnt/hgfs/D/sequencing_data/Herring_PopulationStructure/data_IndividualData_Trimmed_90bp \
-b /mnt/hgfs/D/sequencing_data/Herring_PopulationStructure/data_barcodes/Barcodes_1091_S4_L004_R1_001.txt \
-i gzfastq \
-y fastq \
-e sbfI \
-r \
-c \
-q \
-t 90 \
-E phred33 \
--filter_illumina

#################################################################
# Ellerslie 2015 and Gabriola 2015
#Raw data filename : 1350_S1_L004_R1_001.fastq.gz
 
process_radtags \
-p /mnt/hgfs/D/sequencing_data/Herring_PopulationStructure/data_raw/SingleEnd/Ellerslie2015_Gabriola2015 \
-o /mnt/hgfs/D/sequencing_data/Herring_PopulationStructure/data_IndividualData_Trimmed_90bp \
-b /mnt/hgfs/D/sequencing_data/Herring_PopulationStructure/data_barcodes/Barcodes_1350_S1_L004_R1_001.txt \
-i gzfastq \
-y fastq \
-e sbfI \
-r \
-c \
-q \
-t 90 \
-E phred33 \
--filter_illumina

# Masset 2016 and Newberry 2014
#Raw data filename : 1351_S2_L005_R1_001.fastq.gz

process_radtags \
-p /mnt/hgfs/D/sequencing_data/Herring_PopulationStructure/data_raw/SingleEnd/Masset2016_Newberry2014 \
-o /mnt/hgfs/D/sequencing_data/Herring_PopulationStructure/data_IndividualData_Trimmed_90bp \
-b /mnt/hgfs/D/sequencing_data/Herring_PopulationStructure/data_barcodes/Barcodes_1351_S2_L005_R1_001.txt \
-i gzfastq \
-y fastq \
-e sbfI \
-r \
-c \
-q \
-t 90 \
-E phred33 \
--filter_illumina


# Skidegate 2014 and Spiller Channel 2014
#Raw data filename : 1352_S3_L006_R1_001.fastq.gz #run on 11/5/2017

process_radtags \
-p /mnt/hgfs/D/sequencing_data/Herring_PopulationStructure/data_raw/SingleEnd/Skidegate2014_SpillerChannel2014 \
-o /mnt/hgfs/D/sequencing_data/Herring_PopulationStructure/data_IndividualData_Trimmed_90bp \
-b /mnt/hgfs/D/sequencing_data/Herring_PopulationStructure/data_barcodes/Barcodes_1352_S3_L006_R1_001.txt \
-i gzfastq \
-y fastq \
-e sbfI \
-r \
-c \
-q \
-t 90 \
-E phred33 \
--filter_illumina

# Bute 2004 and Metlakatla 2002
#Raw data filename : 1431_S99_L004_R1_001.fastq.gz #Run on 11/6/2017

process_radtags \
-p /mnt/hgfs/D/sequencing_data/Herring_PopulationStructure/data_raw/SingleEnd/Bute2004_Metlakatla2002 \
-o /mnt/hgfs/D/sequencing_data/Herring_PopulationStructure/data_IndividualData_Trimmed_90bp \
-b /mnt/hgfs/D/sequencing_data/Herring_PopulationStructure/data_barcodes/Barcodes_1431_S99_L004_R1_001.txt \
-i gzfastq \
-y fastq \
-e sbfI \
-r \
-c \
-q \
-t 90 \
-E phred33 \
--filter_illumina

# Sitka 2017 and Knight Inlet 2003
#Raw data filename : 1432_S100_L005_R1_001.fastq.gz #Run on 11/6/2017

process_radtags \
-p /mnt/hgfs/D/sequencing_data/Herring_PopulationStructure/data_raw/SingleEnd/Sitka2017_KnightInlet2003 \
-o /mnt/hgfs/D/sequencing_data/Herring_PopulationStructure/data_IndividualData_Trimmed_90bp \
-b /mnt/hgfs/D/sequencing_data/Herring_PopulationStructure/data_barcodes/Barcodes_1432_S100_L005_R1_001.txt \
-i gzfastq \
-y fastq \
-e sbfI \
-r \
-c \
-q \
-t 90 \
-E phred33 \
--filter_illumina

# Venn Passage 1999 and Skidegate 1999
#Raw data filename : 1433_S101_L006_R1_001.fastq.gz #Run on 11/7/2017


process_radtags \
-p /mnt/hgfs/D/sequencing_data/Herring_PopulationStructure/data_raw/SingleEnd/VennPassage1999_Skidegate1999 \
-o /mnt/hgfs/D/sequencing_data/Herring_PopulationStructure/data_IndividualData_Trimmed_90bp \
-b /mnt/hgfs/D/sequencing_data/Herring_PopulationStructure/data_barcodes/Barcodes_1433_S101_L006_R1_001.txt \
-i gzfastq \
-y fastq \
-e sbfI \
-r \
-c \
-q \
-t 90 \
-E phred33 \
--filter_illumina