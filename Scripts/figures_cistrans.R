library(data.table)
library(ggplot2)
library(ggupset)

dt<-fread("~/Desktop/finalpeaks_withmarkernames_10cMwindow_LODinfo.txt")
cisdt<-dt[type =='cis',]
transdt<-dt[type =='trans',]

x<-data.table(table(dt[,type,trt]))
colnames(x)[3]<-'count'
ggplot(x, aes(fill = type, x = trt, y= count))+
  geom_col(position = 'fill')+
  scale_y_continuous()


cisdt<-cisdt[, transcript.marker := paste0(name,".",markername)]
cisdt<-cisdt[,.(transcript.marker,trt,lod)]
cisdt2<-cisdt[,.(alltrts=list(trt),alllods=list(lod)),by = transcript.marker]
cisdt3<-cisdt2[,.(transcript.marker,alltrts,LOD =unlist(alllods)),by = transcript.marker]
cisdt3<-cisdt3[,-1]

transdt<-transdt[, transcript.marker := paste0(name,".",markername)]
transdt<-transdt[,.(transcript.marker,trt,lod)]
transdt2<-transdt[,.(alltrts=list(trt),alllods=list(lod)),by = transcript.marker]
transdt3<-transdt2[,.(transcript.marker,alltrts,LOD =unlist(alllods)),by = transcript.marker]
transdt3<-transdt3[,-1]

transdt3$alltrts
rownum = 1
for (i in transdt3$alltrts){
  transdt3$alltrts[rownum] <-paste(unlist(i), collapse=" ")
  print(rownum)
  print(transdt3$alltrts[rownum])
  rownum = rownum + 1
}
unlist(transdt3$alltr)

test<-apply(transdt3,FUN = function(
    
  ))
#upset plot with transcript counts
ggplot(transdt2,aes(x=alltrts))+
  geom_bar()+
  scale_x_upset()+
  geom_text(stat='count', aes(label=after_stat(count)), vjust=-1)+
  labs(title = "Trans")+
  theme(text = element_text(size = 20)) 



#cis effect sizes over upset plot
ggplot(transdt3,aes(x=alltrts))+
  geom_jitter(aes(y = LOD,color = as.character(alltrts)))+
  scale_x_upset()+
  #geom_text(stat='count', aes(label=after_stat(count)), vjust=-1)+
  labs(title = "Trans - LOD scores")+
  theme(text = element_text(size = 15))+
  theme(legend.position = "none")



