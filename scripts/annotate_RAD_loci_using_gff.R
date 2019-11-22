# Part 3: Use a gff file to annotate  loci
######################################################################################

# Read a gff annotation file into R 
#Fields
#Fields must be tab-separated. Also, all but the final field in each feature line must contain a value; "empty" columns should be denoted with a '.'

#1. seqname - name of the chromosome or scaffold; chromosome names can be given with or without the 'chr' prefix. Important note: the seqname must be one used within Ensembl, i.e. a standard chromosome name or an Ensembl identifier such as a scaffold ID, without any additional content such as species or assembly. See the example GFF output below.
#2. source - name of the program that generated this feature, or the data source (database or project name)
#3. feature - feature type name, e.g. Gene, Variation, Similarity
#4. start - Start position of the feature, with sequence numbering starting at 1.
#5. end - End position of the feature, with sequence numbering starting at 1.
#6. score - A floating point value.
#7. strand - defined as + (forward) or - (reverse).
#8. frame - One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on..
#9. attribute - A semicolon-separated list of tag-value pairs, providing additional information about each feature.

###################################
# ok, this next part of the code is a pain in the ass to set up
# because it uses bioconductor packages.
# To load a bioconductor package, you need another package named biocLite
# This is an example of how to get biocLite:

#source("http://bioconductor.org/biocLite.R")
#biocLite("GenomicFeatures")


# bioconductor packages are strange and fickle creatures. 
# It took me a minute to download them, and in the process I temporarily fucked up a bunch of my other packages. 
# I apologize if this happens to you.

# Load the necessary bioconductor packages
library(GenomicRanges)
library(rtracklayer)
library(GenomicFeatures)
library(dplyr)
library(tidyr)


# set working directory
setwd("D:/sequencing_data/Herring_Coastwide_PopulationStructure/output_stacks_populations/filtered_haplotypesANDsnps_1104indiv_7261loci/annotation")


###########################################################
# Specify the file names for the input files
locus_file <- "Results_hierfstat_PerLocusFST.txt"

gff_file <- "D:/sequencing_data/Herring_Coastwide_PopulationStructure/output_stacks_populations/filtered_haplotypesANDsnps_1104indiv_7261loci/locus_metadata/Atlantic_herring.gff"

##########################################################
# Specify the file names for the output files

annoted_df <- "batch1_annotatedloci_ingenes.txt"


##########################################################
# Part 1 : Read in the data


# read in df of outlier loci

mydata <- read.table(locus_file, header = TRUE)
head(mydata)

# read in gff file
harengus_gff <- import.gff(gff_file)
harengus_gff
seqnames(harengus_gff)
ranges(harengus_gff)
elementMetadata(harengus_gff)
harengus_gff$Name


# Build a GRanges object from your df of outlier loci
# You can set the "width" around your sequence name, such that you only annotate the exact vicinity
# of your locus, 
# OR
# add some kind of buffer around the locus (say, look 5,000 bp upstream and downstream of locus)
# and collect all of the annotations in the "buffer neighborhood".

# To do this, build an IRanges object containing the ranges you are interested in
my_irange <- IRanges(start=mydata$pos, end=mydata$pos+90)

my_gr2 <- GRanges(seqname = mydata$rname, ranges = my_irange,
                  mcols = mydata)


seqnames(my_gr2)
ranges(my_gr2)
mcols(my_gr2)
my_gr2

# Does your GRanges object overlap with your gff?
## S4 method for signature 'GenomicRanges,GenomicRanges'
#findOverlaps(query, subject,
#maxgap=0L, minoverlap=1L,
#type=c("any", "start", "end", "within", "equal"),
# select=c("all", "first", "last", "arbitrary"),
#algorithm=c("nclist", "intervaltree"),
#ignore.strand=FALSE)


overlaps <- findOverlaps(my_gr2, harengus_gff,  type = "any", ignore.strand = TRUE)
overlaps # as far as I can tell, the overlaps are index positions on the original GRange objects

subjectHits(overlaps)
queryHits(overlaps)

harengus_gff[subjectHits(overlaps)]$Name

harengus_gff[subjectHits(overlaps)]$product
harengus_gff[subjectHits(overlaps)]$protein_id



subject_matches <- harengus_gff[subjectHits(overlaps)]
query_matches <- my_gr2[queryHits(overlaps)]


# Checking some stuff to make sure I did not fuck up any data
subject_matches@seqnames
length(subject_matches@seqnames)

query_matches@seqnames
length(query_matches@seqnames)


subject_matches@elementMetadata$gene
length(subject_matches@elementMetadata$gene)

subject_matches@elementMetadata$gene_biotype
length(subject_matches@elementMetadata$gene_biotype)

subject_matches@elementMetadata$type
length(subject_matches@elementMetadata$type)


subject_matches@elementMetadata$ID

subject_matches@elementMetadata$gene
length(subject_matches@elementMetadata$gene)
subject_matches@elementMetadata$product

query_matches@elementMetadata
length(query_matches@elementMetadata$mcols)
query_matches@ranges@start

# build a dope-ass annotation dataframe! Balllllliiinnnn!!!


#subject_matches@seqnames,
#query_matches@ranges@start,
#subject_matches@elementMetadata$gene_biotype,
#subject_matches@elementMetadata$type,
#subject_matches@elementMetadata$ID,
#subject_matches@elementMetadata$gene,
#subject_matches@elementMetadata$product

match_df <- data.frame(query_matches@elementMetadata,
                       subject_matches@seqnames,
                       query_matches@ranges@start,
                       subject_matches@elementMetadata$type,
                       subject_matches@elementMetadata$ID,
                       subject_matches@elementMetadata$gene,
                       subject_matches@elementMetadata$product)


head(match_df)

# For some odd reason, this gives me a "bloated" df, with rows repeated. 
# Basically, it lists each locus multiple times, if there is a hit to a region, gene, mRNA, etc.
# As a result, the df is like a Russian nesting doll, making it difficult to work with. 
# To get around this issue, I will first subsample the data frame, selecting for the most 
# general level, the region. This should retain information on all loci. 

match_region_df <- match_df %>%
  dplyr :: filter(subject_matches@elementMetadata$type == "region" ) %>%
  dplyr :: distinct(mcols.qseqid,mcols.rname, mcols.Sequence.Name, mcols.pos, mcols.Sequence.Length,
                    mcols.Fst, mcols.Ho, mcols.Fis)


# Now, create a second dataframe that contains some basic information (qseqid, elemet type, gene, product)
# for all loci that are in an mRNA.
# This data frame will contain a lot of duplicated rows (because of the SNPs within each RAD locus),
# so you also have to filter the dataframe by unique instances of the RAD locus name.

match_mRNA_df <- match_df %>%
  dplyr :: filter(subject_matches@elementMetadata$type == "mRNA" ) %>%
  dplyr :: select(mcols.qseqid, subject_matches.elementMetadata.type ,
                  subject_matches.elementMetadata.gene, subject_matches.elementMetadata.product ) %>%
  dplyr :: group_by(mcols.qseqid)  %>%
  dplyr :: top_n(1) %>%
  dplyr :: distinct()

# Combine both data frames together, to get a nice final data frame that retains some information
# (such as scaffold name, etc) about the loci that fo not have gene annotations. Yay! You did it!!

final_df <- left_join(match_region_df, match_mRNA_df, by = "mcols.qseqid")


write.table(final_df, file = annoted_df, 
            quote = FALSE, row.names = FALSE, sep = "\t" )


