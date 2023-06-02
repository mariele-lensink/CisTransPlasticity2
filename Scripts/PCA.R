library(data.table)
library(ggplot2)
library(dplyr)

#setting up the input tables, rils and parents
rils<-as.data.table(read.csv("Data/rils_pheno_alltrts_AGI.csv"))
rils$X<-NULL
rils<-rils[,-2]
rils$ID<-as.character('Ril')
setcolorder(rils,c("ID","treatment",colnames(rils[,2:22747])))
parents<-as.data.table(read.csv("Data/parents_pheno_alltrts_AGI.csv"))
parents$X<-NULL
parents$rep<-NULL
colnames(parents)[1]<-"ID"
parentsavg<-parents[,lapply(.SD,mean),by = c("ID","treatment")]
#need to make the delta group for parents (SA+SW)/((SA+SW)/2)
dpbay<-((parentsavg[1,-c(1:2)]-parentsavg[3,-c(1:2)])/((parentsavg[1,-c(1:2)]+parentsavg[3,-c(1:2)])/2))
dpbay<-data.table(ID="Bay",treatment="delta",dpbay)
dpsha<-((parentsavg[2,-c(1:2)]-parentsavg[4,-c(1:2)])/((parentsavg[2,-c(1:2)]+parentsavg[4,-c(1:2)])/2))
dpsha<-data.table(ID="Sha",treatment="delta",dpsha)

#combined all rils and parents, all treatment groups into 1 table
dt<-rbind(parentsavg,dpbay,dpsha,rils)

#for each treatment group
subdt<-dt[treatment == 'delta',-c('treatment'),with=FALSE]
id<-subdt[,1]
subdt_pca_scaled<-prcomp(subdt[,-1],scale=T)

pcs<-data.frame(subdt_pca_scaled$x)
pcs$ID<-subdt$ID
pcs<-as.data.table(pcs)

pca_var<-subdt_pca_scaled$sdev^2
pca_var_per<-round(pca_var/sum(pca_var)*100,2)

png("Figures/pca_delta.png")
ggplot(pcs,aes(x=PC1,y=PC2,col=ID))+
  geom_point()+
  theme_bw()+
  geom_point(data=pcs[ID=="Bay"],aes(x=PC1,y=PC2,col=ID),
               pch=21,size=4,colour="red")+
  geom_point(data=pcs[ID=="Sha"],aes(x=PC1,y=PC2,col=ID),
             pch=21,size=4,colour="blue")+
  ggtitle("Delta (SA-SW)/((SA+SW)/2)")+
  xlab(paste0("PC1 - ",pca_var_per[1],"%"))+
  ylab(paste0("PC2 - ",pca_var_per[2],"%"))+
  theme(text = element_text(size=20))
dev.off()
