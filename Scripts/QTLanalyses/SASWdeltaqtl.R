#load packages
library(qtl2)
library(readxl)
library(Rcpp)
library(data.table)


#delta = (SA-SW)/((SA+SW)/2)
###############CREATING DELTA PHENO TABLE FROM SA AND SW PHENO TABLES####################################
#read in CSVs
setwd("/Users/marielelensink/Documents-local/SalicylicAcid/delta")
SAdf<-fread("../SWqtl/BayxSha_SWpheno.csv")
SWdf<-fread("../SAqtl/BayxSha_SApheno.csv")

rils<-SAdf[,2]
sam<-as.matrix(SAdf)
swm<-as.matrix(SWdf)
sam<-sam[,-c(1:2)]
swm<-swm[,-c(1:2)]
#calculate difference divided by mean
delta <- as.data.table((sam-swm)/((sam+swm)/2))
delta<-cbind(rils,delta)


#make this a CSV
write.csv(delta,file = "BayxSha_deltapheno.csv")

#########################################################################
#########QTL ANALYSIS########################
setwd("/Users/marielelensink/Documents-local/SalicylicAcid/delta/")
bayxsha_delta<-read_cross2("BayxSha_delta.yaml",quiet = FALSE)

#calculate conditional genotype probabilities given the marker data at each putative QTL position
##insert pseudomarkers into the genetic map
map_delta<-insert_pseudomarkers(bayxsha_delta$gmap,step=0)
#calculate genotype probabilities
pr_delta<- calc_genoprob(bayxsha_delta,map_delta,error_prob = 0.002)
#analyzing marker density
#determine grid of pseudomarkers
grid_delta <- calc_grid(bayxsha_delta$gmap, step=0)

####run the genome scan####

#scan1() takes as input the genotype probabilities, a matrix of phenotypes, 
#optional additive and interactive covariates, and the special X chromosome covariates.
out_delta <- scan1(pr_delta, bayxsha_delta$pheno)
#output is LOD scores.
#LOD stands for "logarithm of the odds." 
# positions X phenotypes

#find peaks
peaks_delta<-find_peaks(out_delta, map_delta, threshold=2, drop=1.5)
fwrite(peaks_delta,"Data/QTLAoutput/deltapeaks_save_final.txt")
peaks_delta<-fread("Data/QTLAoutput/deltapeaks_save_final.txt")
hist(table(peaks_delta$lodcolumn))

