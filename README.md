# FourSite
Four site is a method to estimate heterozygosity directly from genomic sequencing reads. It was specifically designed for low coverage data, but can be applied to higher coverage data too. A minimum of four reads at any site is required. Comparisons with high coverage data (15x) show heterozygosity estimates computed with FourSite had excellent concordance with estimates derived from GaTK.

#### Two scripts are necessary to compute heterozygosity.

1) CountReadFoursite.py  This python script (python 2.6) iterates across positions in a bam file.  For each site, it samples four sequencing reads and records whether i) all four reads are the same base, ii) two reads are one base and two reads are a different base, or iii) one read is one base, and three reads are a different base. 

2) Foursite.R  Takes the output from CountReadFoursite.py and computes the likelihood of the heterozygosity and sequencing error rate as function of these counts across a particular functional category

Detailed README files for the above programmes are provided in specific folders
