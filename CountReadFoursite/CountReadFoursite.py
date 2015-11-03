#! /usr/bin/env python
## Written for python2.6
## to do change samefile = to bamfile = 
## update the usage error
import pysam
import random
import sys

if len(sys.argv) != 9:
	print "\n**** USAGE ERROR **** \n"
	print "python script.py number_iterations minimum_coverage maximum_coverage sites_file chromosome sampleid bamfile outfilename "
	sys.exit("exiting ERROR !!!!")

#########################
###### ***** User specified input
#########################
Num_iter=int(sys.argv[1])  # Number of iterations
Min_cov=int(sys.argv[2])  # Min coverage to consider a site
Max_cov=int(sys.argv[3])  # Max coverage to consider a site
Sites_file=sys.argv[4]  # This is a text file with just the position in it, one position per line.
chrom_sel = sys.argv[5] # The chromosome you are analysing. e.g. chr25
Sample = sys.argv[6] # ID of the sample in the bam file
samfile = pysam.Samfile(sys.argv[7], "rb") ## The bamfile you want to analyse
Outfilename = sys.argv[8]

if Min_cov < 4: 
	sys.exit(" ** Exiting error. You specified minimum coverage < 4. 4 is the minimum allowed")

### Outfile - open and write header
Outfile = open(Outfilename, 'w')
Outfile.write('#chromosome\titeration_number\tcount_excluded_triallelic_sites\tcount_exclued_error_sites\tcount_sites_excluded_for_lowhighcoverage\tcount_parsed_sites_from_sitesfile\ttype_of_site\tSampleID\tFourSite_count_parsed_sites_calcbasedon\tFoursite_count_4readsSame\tFoursite_count_3RdSame1Diff\tFoursite_count_2readsame2readdif\n')


#########################
###### ***** Get sites of interest
#########################
## Sites file has the list of sites you are interested in looking at. There is one file per chromosome. Read in the sites as a list then make it a set as this is more efficient

f = open(Sites_file, "rU")
listofsites = [int(x) for x in f]
sites = set(listofsites)

#########################
###### ***** Counters
#########################
count_iterations=0
count_screen_sites=0  # total num sites looked at - inc sites with too high coverage/N's
count_exc_cov=0 # count of sites excluded due to low or high coverage
count_error=0 # error check - if this is not 0, something went wrong
count_triallelic=0 # sites where number of bases > 2 after excluding N's

## Foursite counters
FS_4RdSame=0
FS_3RdSame1Diff=0
FS_2RdSame2Diff=0
FS_Ctsitesparse=0

#########################
###### ***** Functions
#########################
def site_homo(bases):
	"checks to see if the bases are identical or not"
	return all(x == bases[0] for x in bases)

#########################
###### ***** Main Script
#########################

## Check of position is site of interest - if not skip that site (may not use this if conduct computations only on sites of interest)
## For each site of interest, append all of the calls from the different reads to the list 'calls'
## Start counting sites
while count_iterations < Num_iter:
	count_iterations += 1
	for pileupcolumn in samfile.pileup(chrom_sel):	
		calls=[]
		if pileupcolumn.pos not in sites:
#		if pileupcolumn.pos == 0:  ## You would use this if you wanted the whole chromosome
			pass
		else:
			## get the bases for a position & store it in calls
			for pileupread in pileupcolumn.pileups:
				calls.append(pileupread.alignment.seq[pileupread.qpos])
			count_screen_sites += 1			
### 		List comprehension to remove N's from the list of potential calls to sample
			calls = [y for y in calls if y != 'N']
			num_bases=set(calls)

			## Skip sites where coverage is too low or too high - count the number of these sites 
			## Skip sites if they are triallelic and count these sites
			## Else randomly sample 4 reads from calls
			if len(calls) > Max_cov or len(calls) < Min_cov:
				count_exc_cov += 1
			elif len(num_bases) > 2:
				count_triallelic += 1
			elif len(num_bases) <= 0:
				count_error += 1
			else:
				fourreads=random.sample(calls,4)

			## Count if the call is the same, not the same, or if it isn't the same, different or an N (i.e. an error check)
	
				if len(fourreads) != 4:
					sys.exit("***ERROR exiting - script is meant to be written to sample four reads, but more were sampled. Something is wrong - bug?")	
				else:
					FS_Ctsitesparse+=1
					D_fourreads={}
					## This makes a dicitonary of what is in your four reads e.g. {'A': 3, 'C': 1}
					for item in fourreads:
						D_fourreads[item]=D_fourreads.get(item, 0) + 1

					## Extract out the dictionary values e.g. [3,1] for the above example
					val_fourreads=D_fourreads.values()
					key_fourreads=D_fourreads.keys()

					### Four reads the same [A,A,A,A], [4]
					if val_fourreads == [4]:
						FS_4RdSame+=1			
					### Two reads same : Two reads same e.g. [A,A,C,C] = [2,2]	
					elif val_fourreads == [2,2]:
						FS_2RdSame2Diff+=1
					### Three reads same : One diff e.g. [A,A,A,C] = [1,3] or [3,1]
					elif val_fourreads == [3,1] or val_fourreads == [1,3]:
						FS_3RdSame1Diff+=1
					else:
						print "Your read dictionary: ", D_fourread
						print "Your read values: ", val_fourreads
						sys.exit("***ERROR exiting - We extracted four reads, and counted whether they were all the same, 3 reads one base, 1 read another base OR 2 reads one base 2 reads another base.  You are none of these. Something is wrong - bug?")

	Outfile.write(str(chrom_sel) + "\t" + str(count_iterations) + "\t" + str(count_triallelic) + "\t" + str(count_error) + "\t" + str(count_exc_cov) + "\t" + str(count_screen_sites) + "\t" + str(Sites_file) + '\t' + str(Sample) + '\t' + str(FS_Ctsitesparse) + '\t' + str(FS_4RdSame) + '\t' + str(FS_3RdSame1Diff) + '\t' + str(FS_2RdSame2Diff) + '\n')

	
## reset the counters	
	count_screen_sites=0  # total number of sites looked at - including ones with too high coverage, or N's
	count_exc_cov=0 # count of sites excluded due to low or high coverage
	count_error=0 # error check - if this is not 0, something went wrong
	count_triallelic=0
	FS_4RdSame=0
	FS_3RdSame1Diff=0
	FS_2RdSame2Diff=0
	FS_Ctsitesparse=0

samfile.close()
Outfile.close()
