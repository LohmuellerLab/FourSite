# FourSite
Method to estimate heterozygosity directly from genomic sequencing reads

Two scripts are necessary to compute heterozygosity.

1) CountReadFoursite.py  This script iterates across positions in a bam file.  For each site, it samples four sequencing reads and records whether i) all four reads are the same base, ii) two reads are one base and two reads are a different base, or iii) one read is one base, and three reads are a different base. 

2) Foursite.R  Takes the output from CountReadFoursite.py and computes the likelihood of the heterozygosity and sequencing error rate as function of these counts across a particular functional category
