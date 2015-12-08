# CountReadFoursite

CountReadFoursite was developed to provide a method that estimates heterozygosity and performs equally well on low coverage and high coverage data. A minimum of 4 reads per site are required for a calculation to be made for a position. The code is written to calculate heterozygosity within specific regions of the genome, and is run on each chromosome seperately. We suggest that you run multiple iterations, as the calculations are based on randomly drawing of the reads observed at a site. The script can be modified to run on all sites within a chromosome though [more details on this will be provided soon].

#####usage 
    python CountReadFoursite.py number_iterations minimum_coverage maximum_coverage sites_file chromosome sampleid bamfile outfilename 

number_iterations # number of iterations to run
minimum_coverage # exclude sites with < this coverage. MUST be 4 or greater
maximum_coverage # exclude sites with > this coverage 
sites_file # name of the file with the sites of interest
chromosome # name of chromosome being assesed. This should match how chromosomes are labelled in your bam e.g. chr25
sampleid # name of sample being assessed
bamfile # name of bam infile
outfilename # name for results outfile

If any of the files specified are not in the current directory, you must provide the path

#### Dependancies
CountReadFoursite.py requires python 2.6 and the python pysam module.

#### Input file
A bam file and associated bai file is needed as input for CountReadFoursite.py. It is critical that your bai file is labelled .bam.bai i.e. if your bam is myfile.bam, your bai file should be myfile.bam.bai. See example in test dataset

Foursite.R takes the output from CountReadFoursite.py as the input file.

#### Sites file
The sites file is a text file with one site per line. 
You need a seperate sites file for each chromosome.
See example in the test dataset 

#### Pre-processing your bam file
In order for Foursite to accurately estimate heterozygosity, we have found it very necessary to pre-process your bam file to reduce error rates. We did the following steps.  1) Reads were first trimmed using sickle and scythe (https://github.com/ucdavis-bioinformatics) to remove adaptor contamination and low quality sequence at the end of the read (Q<30). 2) Base quality score recalibration (part of GaTK) was applied to the reads. 3) We used a python script to change any bases Q<30 within a read to N, so that are ignored downstream (https://github.com/cdmarsden/replace_lowqualitybases).[** In the next version of CountReadFoursite.py we hope to remove the need for this step by assessing the quality score of the base within the script]. 4) We removed reads with low overall quality (i.e. length < 40, or where > 20% of bases had Q<30 (https://github.com/cdmarsden/remove_lowqualityreads) 

It is noteworthy that within the script we count and output the number of 'triallelic sites'. i.e. sites where 3 bases are observed. In our experience, the presence of non-negligible numbers of triallelic sites is generally indicative of poor sequencing quality.

#### Test dataset
Provided are some test files
myinfile.bam and myinfile.bam.bai
chr25_sitesofinterest.txt
######To run:
    python CountReadFoursite.py 10 4 80 chr25_sitesofinterest.txt chr25 dog0920 myinfile.bam outfile_test.txt

