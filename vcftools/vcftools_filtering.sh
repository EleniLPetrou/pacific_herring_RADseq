# The purpose of this script is to filter a vcf file for things like missing data, minor allele frequency, etc.
# vcftools takes a vcf file as input

# Explanation of parameters:
# --vcf = name of vcf input file
# --max-missing <float> = Exclude sites on the basis of the proportion of missing data (defined to be between 0 and 1, where 0 allows sites that are completely missing and 1 indicates no missing data allowed).
# --maf <float> = Include only sites with a Minor Allele Frequency greater than or equal to the "--maf" value
# --minDP <float> = Exclude  all genotypes with a sequencing depth below that specified by this option
# --mac <integer> = Include only sites with Minor Allele Count greater than or equal to the "--mac" value


# First filtering step - rough cut, to remove loci with lots of missing genotypes, and a minimum sequencing depth of 10 reads
vcftools --vcf batch_1.vcf --max-missing 0.75 --mac 3  --minDP 10 --recode --recode-INFO-all --out batch_1.g8mac3dp10

# Assess individual levels of missing data
# --missing-indv = 
vcftools --vcf batch_1.g8mac3dp10.recode.vcf --missing-indv

# Use mawk to search for lines that do not contain the expression "IN". 
# Then, pipe those lines to the cut command, which removes sections from each line of files. 
# In this case, we will cut out the fifth field (or column), and save it to a file; this can be used for plotting.
mawk '!/IN/' out.imiss | cut -f5 > totalmissing

# Create a "badlist" of individuals with more than 20% missing data. 
mawk '$5 > 0.20' out.imiss | cut -f1 > lowDP.indv

# Now that there is a list of individuals to remove, we can feed that directly into VCFtools for filtering.
#--remove <filename>
#              Provide a file containing a list of individuals  to  exclude  in
#              subsequent  analysis.  Each individual ID (as defined in the VCF
#             headerline) should be included on a separate line. If  both  the
#              --keep and the --remove options are used, then the --keep option
#              is execute before the --remove option.

vcftools --vcf batch_1.g8mac3dp10.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out batch_1.g8mac3dplm

# Now that we have a filtered vcf file, we can extract the list of SNP names that we want to retain. 
# We will use this list as a "whitelist" in stacks-populations.
awk '{print $1,$3, $4, $5}' batch_1.g8mac3dplm.recode.vcf > whitelist.txt

# Using the filtered vcf file, calculate the mean sequencing depth for each SNP
# --site-mean-depth = Generates a file containing the mean depth per site averaged across all individuals. This output file has the suffix ".ldepth.mean".
vcftools --vcf batch_1.g8mac3dplm.recode.vcf --site-mean-depth

