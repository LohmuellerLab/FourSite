# Foursite.R
README Foursite_shell.github.sh 

Written by Kirk Lohmueller, 2014-2015, UCLA.

This shell script calls an R script that will compute a maximum likelihood estimate of the per site heterozygosity (pi) and error rate across a number of sites within an individual genome.

This method is based on the underlying model of Lynch (2010 MBE) with modifications as described in Marsden et al.

As input, the program requires the number of sites where all 4 sequencing reads show the same base, the number of sites where 3 reads show one base and 1 read shows an alternative base, and finally, the number of sites where 2 reads show one base and 2 show the alternate base. These counts are provided in the input (see below) with 1 individual genome per line. See "outfile_test.txt" as an example of the input file for this script. 

If you wish to estimate heterozygosity for different functional categories, include the counts for each category on distinct lines in the input.

Usage:

./Foursite_shell.github.sh outfile_test.txt

outfile_test.txt is the input file.

Several parameters within the shell script can be modified:

pi_low This is the lower bound of pi to search over (per site)
pi_high This is the upper bound of pi to search over (per site)
error_low This is the lower bound of per read error rate (per site) to search over
error_high This is the upper bound of per read error rate (per site) to search over
STEP This is the step size in grid search

Note, pi_low and pi_hi should be set to sensible values for the organism of interest. Setting them too low or too high will yield odd behavior. If your MLE is at the boundary, you should increase the grid size. 



Output:

The script will generate a file with *.out.txt which will contain the output. The output will consist of the entire input with 3 additional columns:

error_est Estimated per read per site error rate
pi Estimated pi per site
llike Log-likelihood at the MLEs of the error rate and pi


Questions, bugs, email Kirk:
klohmueller@ucla.edu



