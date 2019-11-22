# The first step for all bayenv 2 analyses is to estimate the covariance matrix. This should be done using a large number of markers with no or loose linkage between them. 
# The command line for estimating the matrix is:

# explanation of terms:
# -i : input file (SNPSFILE)
# -p : number of populations
# -k : number of iterations (Gunter and Coop used 100000 MCMC iterations in their 2013 Genetics paper)

# Place all input files to this directory, and navigate to:
cd ./Downloads/bayenv2/


./bayenv2 -i batch_1_firstsnp_bayenv.txt -p 28 -k 100000 > herring_matrix.out

# The rows and columns in the covariance matrix appear in the same population order as they appeared in the count file. 
# the MATRIXFILE used in the subsequent steps of baynv2 can be obtained by copying and pasting the last printed matrix from the output of the program.


# Calculate bayes factors for all SNPS of SNPSfile

# the archive of bayenv2 contains a small bash script (calc_bfs.sh) to calculate Bayes factors for all SNPs in SNPSfile. 
#to run it, type:

#./calc_bf.sh SNPSFILE ENVIRONFILE MATRIXFILE NUMPOPS NUMITER NUMENVIRON

./calc_bf.sh batch_1_firstsnp_bayenv.txt herring_environfile.txt herring_matrix_final.txt 28 1000 1


# It is also recommended that we do a non-parametric test, and see whether results from bayesian model and this test overlap.
# Sometimes the linear model underlying the Bayes Factors might not be correct, or outliers might misguide the model. 
# For such cases, Bayenv2 calculates the standardized allele frequencies (X), from which the covariance structure among population was removed. 
# Then it calculates the non-parametric Spearman's rank correlation coefficient (rho) in addition to Bayes factors, when the -c flag is set. 


# test for a single snp

./bayenv2 -i testSNP.txt -m herring_matrix_final.txt -e herring_environfile.txt -p 28 -k 10 -n 1 -t -r 429

# test for a single snp on their data

./bayenv2 -i testSNP.txt -m hgdp_matrix_1 -e PCs.env -p 28 -k 10 -n 1 -t -r 429