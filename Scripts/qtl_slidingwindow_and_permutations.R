library(data.table)
library(ggplot2)
library(rlist)
library(reshape2)

peaks_SA<-fread("Data/QTLAoutput/SA_peakssave_final.txt")
peaks_SW<-fread("Data/QTLAoutput/SW_peaks_save_final.txt")
peaks_delta<-fread("Data/QTLAoutput/deltapeaks_save_final.txt")
#reordering peaks by chromosome then position on chromosome (not sure if necessary but it makes me feel better)
SA<-peaks_SA[order(peaks_SA$chr,peaks_SA$pos),]
SW<-peaks_SW[order(chr,pos),]
delta<-peaks_delta[order(chr,pos),]

#split by chromosome
SAlist<-split(SA,by='chr')
SWlist<-split(SW,by='chr')
deltalist<-split(delta,by='chr')

#BUILD SLIDING WINDOW FUNCTION#
###input list of qtl positions, window size, stepsize, and total length in centimorgans
centimorgan_slider<-function(positions,window,stepsize,cms) {
#set initial window ranges before loop 
  rangemin<-1
  rangemax<-rangemin+window-1
  peakcounts<-c()
while(rangemax <= cms){
    n<-sum(positions >= rangemin & positions < rangemax)
    peakcounts<-c(peakcounts,n)
    
    rangemin<- rangemin+stepsize
    rangemax<-rangemax+stepsize
    
  }
return(peakcounts)
}

#running the sliding window over list of tables split by chromosome
SAwindowcounts_split<-list.cbind(lapply(SAlist, function(i) data.table(centimorgan_slider(i[,pos],10,1,105))))
SWwindowcounts_split<-list.cbind(lapply(SWlist, function(i) data.table(centimorgan_slider(i[,pos],10,1,105))))
deltawindowcounts_split<-list.cbind(lapply(deltalist, function(i) data.table(centimorgan_slider(i[,pos],10,1,105))))

#reformatting tables for ggplot
index<-c(1:length(SAwindowcounts_split$`1.V1`))
colnames(SAwindowcounts_split)<-c('1','2','3','4','5')
SAcounts<-cbind(SAwindowcounts_split,index)
SAcounts$trt<-"SA"
colnames(SWwindowcounts_split)<-c('1','2','3','4','5')
SWcounts<-cbind(SWwindowcounts_split,index)
SWcounts$trt<-"SW"
colnames(deltawindowcounts_split)<-c('1','2','3','4','5')
deltacounts<-cbind(deltawindowcounts_split,index)
deltacounts$trt<-"delta"
#combine tables and finish formatting
peaks_all<-rbind(SAcounts,SWcounts,deltacounts)
peaks_all_melt<-melt(peaks_all,id.vars=c('trt','index'))
colnames(peaks_all_melt)<-c("Treatment",'Index','Chromosome','Count')
#make sliding window index represent cMs
tempind<-peaks_all_melt$Index
#tempind<-tempind*2.5-2.5
peaks_all_melt$Position<-tempind
#plot
png("Figures/QTLfrequency.png")
ggplot(peaks_all_melt,aes(x=Position,y=Count,color = Treatment))+
  geom_point(size = .5)+geom_line()+facet_wrap(peaks_all_melt$Chromosome,nrow = 1)+
  theme_bw()+
  ggtitle("Frequency of QTLs across the A. thaliana genome")
dev.off()

#adjust delta by dividing by 2.75 (normalize to SA and SW)
peakstable<-data.table(peaks_all_melt)
peaks_all_melt<-data.table(peaks_all_melt)
peakstable$Count <- as.double(peakstable$Count)
normpeaks<-peakstable[Treatment=="delta",Count:=Count/2.75]
normpeaks[Treatment=='SA']
peaks_all_melt[Treatment=='SA']

png("Figures/QTLfrequency_deltascaled.png")
ggplot(normpeaks,aes(x=Position,y=Count,color = Treatment))+
  geom_line()+facet_wrap(peaks_all_melt$Chromosome,nrow = 1)+
  theme_bw()+
  theme(text = element_text(size = 15),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Frequency of QTLs across the A. thaliana genome (scaled - delta/2.75")
dev.off()  

#####RANDOMIZE QTLs, 1,000 permutations######
#get number of qtls
deltaqtlcount<-length(delta$pos)
saqtlcount<-length(SA$pos)
swqtlcount<-length(SW$pos)

#get length of each chromosome, some to total length in centimorgans
totalcms<-max(delta[chr==1,pos])+max(delta[chr==2,pos])+max(delta[chr==3,pos])+
  max(delta[chr==4,pos])+max(delta[chr==5,pos])

#generate random qtls for each trt in 1000 permutations
deltaperms<-data.table()
SAperms<-data.table()
SWperms<-data.table()
for (i in 1:1000){
  deltatemp<-runif(deltaqtlcount,0,totalcms)
  deltanewcol<-as.data.table(deltatemp)
  deltaperms<-cbind(deltaperms,deltanewcol)
  satemp<-runif(saqtlcount,0,totalcms)
  sanewcol<-as.data.table(satemp)
  SAperms<-cbind(SAperms,sanewcol)
  swtemp<-runif(swqtlcount,0,totalcms)
  swnewcol<-as.data.table(swtemp)
  SWperms<-cbind(SWperms,swnewcol)
}

#run randomized sliding window on all of the permutations
deltaperm_windowcounts<-as.data.table(apply(deltaperms,2,function(x) centimorgan_slider(x,10,1,430))) 
SAperm_windowcounts<- as.data.table(apply(SAperms,2,function(x) centimorgan_slider(x,10,1,430))) 
SWperm_windowcounts<- as.data.table(apply(SWperms,2,function(x) centimorgan_slider(x,10,1,403))) 
#get 5% cutoffs for values in each group
deltacutoff<-quantile(as.vector(deltaperm_windowcounts),0.95)
SAcutoff<-quantile(as.vector(SAperm_windowcounts),0.95)
SWcutoff<-quantile(as.vector(SWperm_windowcounts),0.95)

delta5<- deltaperm_windowcounts>=291


peaks_all_melt<-data.table(peaks_all_melt)
countsum<-peaks_all_melt[,sum(Count),by = Treatment]
sum(countsum$V1)
