## Running SnpEff version 5.0.0 to annotate variants

[SnpEff](http://snpeff.sourceforge.net/index.html) is a program that adds variant annotations to loci in a VCF file and also predicts the functional effect of a variant. The program requires:
- a VCF file: Contains variants
- a GFF file: Contains genome annotations
- Reference genome in FASTA format: Constitudes the reference for functional effect simulation



I installed snpeff using conda, into this directory: /home/ubuntu/miniconda2/share/snpeff-5.0-0

``` bash
conda install -c bioconda snpeff
```


I had downloaded the annotation file for the Atlantic herring genome assembly GCA_900700415 Ch_v2.0.2 from [here]https://www.ncbi.nlm.nih.gov/genome/15477?genome_assembly_id=495882 
This annotation file is named *GCF_900700415.1_Ch_v2.0.2_genomic.gff*

I also saved the "replicon info" table as a tab-delimited text file called *lookup_table_chromosome_names.txt*
The first couple of lines look like this:

```
Type	Name	RefSeq	INSDC	Size(Mb)	GC%	Protein	rRNA	tRNA	Other RNA	Gene	Pseudogene
Chr	1	NC_045152.1	LR535857.1	33.08	44.4	2,151	-	15	156	1,440	22
Chr	2	NC_045153.1	LR535858.1	33.01	43.9	1,862	-	11	152	1,174	30
Chr	3	NC_045154.1	LR535859.1	32.53	44.5	2,041	1	74	121	1,294	22
```

I did this because the chromosome names in the gff and .fna and .vcf files do not match. The gff file uses RefSeq names, while my other files name the chromosomes using the INSDC names. So, I need to rename the chromosomes in the .gff file (NC_...) so they match the chromosome names in the .vcf and .fna (genome) files (LR...). The only way I can think of doing this is by writing a python script that uses a dictionary to rename strings. 

``` bash
DIR=/mnt/hgfs/D/sequencing_data/Atlantic_herring_genome/GCA_900700415/refseq_gff
GFF_FILE=GCF_900700415.1_Ch_v2.0.2_genomic.gff
LOOKUP_TABLE=lookup_table_chromosome_names.txt
OUTFILE=GCF_900700415.1_Ch_v2.0.2_genomic.renamechroms.gff

cd $DIR
python rename_chroms.py $LOOKUP_TABLE $GFF_FILE $OUTFILE

```

The python script *rename_chroms.py* looks like this:

``` python
#rename_chroms.py

import sys
import re

#specify the file names
table_name = sys.argv[1]
gff_name=sys.argv[2]
out_name=sys.argv[3]

#open the files for reading & writing
table_file = open(table_name, "r")
gff_file = open(gff_name, "r")
out_file = open(out_name,"w")

# store specific values from the table_file in a dictionary:
my_dict = dict()

for line in table_file:
  if not line.startswith("#"): #ignore lines with comment char #
    line_list = line.strip().split("\t")
    #print line_list
    my_dict[line_list[2]] = line_list[3] #assign key-value pair for this file

table_file.close()

# print the keys to check that it worked
keys = my_dict.keys()
keys.sort()

for key in keys:
  print key + "\t" + my_dict[key]


# Open the gff file. :

for line in gff_file:
  for key in my_dict:
    if key in line:
      line = line.replace(key,my_dict[key])
  out_file.write(line) # write to output.txt

out_file.close()
```

Awesome. Now I hope I can follow Angela's instructions to build a snpeff database: 

## Build a new database with annotations from a GFF file
Detailed instructions on how to create a new database can be found [here](http://snpeff.sourceforge.net/SnpEff_manual.html#databases)

1. I directly modified the `snpEff.config` file by adding the herring genome to the non-standard databases field (line 127):

```
#---
# Non-standard Databases
#---

# Herring genome
Herring.genome : Herring
```

3. I saved a copy of the files required to create a SnpEff database for the herring genome in specific directories. First, I created a directory called `genomes` within the `data` directory of SnpEff (I used the same name as in the config file) `data/genomes/`. I copied to `data/genomes` the herring reference genome (`GCA_900700415.1_Ch_v2.0.2_genomic.fna`) with a modified name `Herring.fa`. Then, within `data` directory, I created a new folder called `Herring`. I copied to it the  GFF file with annotations for the herring genome that we made using step above, which is called `GCF_900700415.1_Ch_v2.0.2_genomic.renamechroms.gff`, and I changed its name for `genes.gff`.

Then, I used this command to create a new database:

``` bash

# Navigate to directory with snpeff
cd /home/ubuntu/miniconda2/share/snpeff-5.0-0

# run database building command
java -jar snpEff.jar build -gff3 -v Herring

```

4. Next, I used SnpEff to annotate my vcf file (run within the directory that has the `snpEff.jar` file):

```
#Specify file names and directories with data
VCF_DIR=/mnt/hgfs/D/sequencing_data/Herring_Coastwide_PopulationStructure/output_stacks_populations/filtered_haplotypesANDsnps_1104indiv_6718loci/vcf
VCF_FILE=batch_1_firstsnp_GCA900700415_mapq20.recode.sorted.vcf
OUT_FILE=batch_1_firstsnp_GCA900700415_mapq20.recode.sorted.annotations.vcf

# Navigate to directory with snpeff
cd /home/ubuntu/miniconda2/share/snpeff-5.0-0

java -Xmx100g -jar snpEff.jar -v Herring $VCF_DIR'/'$VCF_FILE > $VCF_DIR'/'$OUT_FILE

```
Oh my gosh. It actually worked! I am amazed.
 
# Install GATK and run

Install gatk by downloading .zip file from web and extraction to ~/Downloads folder on ubuntu machine. Directions can be found [here]https://gatk.broadinstitute.org/hc/en-us/articles/360036194592-Getting-started-with-


GATK syntax:
gatk [--java-options "jvm args like -Xmx4G go here"] ToolName [GATK args go here]

Navigate to directory with GATK jar files

``` bash 

# Specify some directories

GENOME_DIR=/mnt/hgfs/D/sequencing_data/Atlantic_herring_genome/GCA_900700415/fasta
GATK_DIR=/home/ubuntu/Downloads/gatk-4.1.9.0
VCF_DIR=/mnt/hgfs/D/sequencing_data/Herring_Coastwide_PopulationStructure/output_stacks_populations/filtered_haplotypesANDsnps_1104indiv_6718loci/vcf
VCF_INFILE=batch_1_firstsnp_GCA900700415_mapq20.recode.sorted.annotations.vcf
FASTA_INFILE=GCA_900700415.1_Ch_v2.0.2_genomic.fna
OUTFILE=batch_1_firstsnp_GCA900700415_mapq20.recode.sorted.annotations.ANN

###########################
# Go to directory where gatk files are
cd $GATK_DIR

#index the fasta file
samtools faidx $GENOME_DIR'/'$FASTA_INFILE

# create the gatk fasta dictionary sequence file
./gatk CreateSequenceDictionary -R $GENOME_DIR'/'$FASTA_INFILE

# run gatk VariantsToTable

./gatk VariantsToTable \
-R $GENOME_DIR'/'$FASTA_INFILE \
-V $VCF_DIR'/'$VCF_INFILE \
-O $VCF_DIR'/'$OUTFILE \
-F CHROM -F POS -F REF -F ALT -F ANN

```
The results are in a tab-delimited text file that is marginally more user-friendly than the annotated vcf.


