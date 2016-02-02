#!/bin/bash

#Here is the foursite code. This method will calcualate heterozygosity and an error rate from 4 sequencing reads across several sites within a genome.

#Code by Kirk Lohmueller, UCLA, November 2015
#Code is used in Marsden et al.

#useage: need to specify input filename on the command line
echo 'y'
pi_low=0.00005 #lower bound of pi to search over (per site)
pi_high=0.003 #upper bound of pi to search over (per site)
error_low=0.00001 #lower bound of per read error rate (per site) to search over
error_high=0.0001 #upper bound of per read error rate (per site) to search over
step=0.00005 #step size in grid search




R --quiet --no-save > /dev/null <<EOF

#this this an R script to compute heterozygosity from read counts

#


get_like_het<-function(het_prob) #This function will compute prob of read counts | het for a site
{	
	het_prob[1]<-dbinom(0,4,0.5)*2
	het_prob[2]<-dbinom(1,4,0.5)*2
	het_prob[3]<-dbinom(2,4,0.5)
	return(het_prob)
}


get_like_hom<-function(hom_prob,error) #This function will compute prob of read counts | hom for a site
{	
	hom_prob[1]<-(1-error)^4+error^4
	hom_prob[2]<-4*error*(1-error)^3+4*error^3*(1-error)
	hom_prob[3]<-choose(4,2)*error^2*(1-error)^2
	return(hom_prob)
}



get_like_for_point<-function(error,pi,het_prob,hom_prob,counts)
{
	out<-numeric()
	#first, compute the per site values:
	hom_prob<-get_like_hom(hom_prob,error)
	likelihood<-pi*het_prob+(1-pi)*hom_prob
	final_log_like<-sum(counts*log(likelihood))
	out<-cbind(error,pi,final_log_like)
	return(out)
}


search_grid<-function(error_low,error_high,pi_low,pi_high,het_prob,hom_prob,counts,step)
{
	final_out<-numeric()
	out<-numeric()
	for(i in seq(error_low,error_high,by=step))
	{
		for(j in seq(pi_low,pi_high,by=step))
		{
			out<-get_like_for_point(i,j,het_prob,hom_prob,counts)
			final_out<-rbind(final_out,out)
		}
	}
	
	
	return(final_out)
}




#read in data from a file:
data<-read.table("$1",header=T) #input datafile

MLE_vec<-numeric()

for (wl in c(1:length(data[,1])))
#for (wl in c(1:2))
{
	
	counts<-c(data[wl,10],data[wl,11],data[wl,12]) #these are the 3 colums with teh counts of 4 reads teh same, 3 reads 1 base (and 1 the other), and 2 reads one base (and 2 the other)
	#ok, here's the stuff for the acutal program:

	het_prob<-numeric() #matrix for calculations
	hom_prob<-numeric() #matrix for calcualtions

	het_prob<-get_like_het(het_prob) #this is a contsant

#now define the ranges of pi and error rates, these can be changed: 

#pi_low<-0.00005
#pi_high<-0.003
#error_low<-0.00001
#error_high<-0.0001
#step<-5e-5

final_out<-search_grid($error_low,$error_high,$pi_low,$pi_high,het_prob,hom_prob,counts,$step)

MLE<-subset(final_out,final_out[,3]==max(final_out[,3]))

MLE_vec<-rbind(MLE_vec,MLE)

}

out<-cbind(data,MLE_vec)

write.table(out, "$1.out.txt",quote=F,sep="\t", row.names=F,col.names=c(names(data),"error_est","pi","llike")) #use for finer grid 

EOF


