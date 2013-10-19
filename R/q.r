# program to get the expression level of mRNA from qPCR data using 2deltadeltaT method
# by JessicaBan, 2013 04 24

#housekeeping gene
hk = "HPRT"
#target gene
tg = "DNMT3A"
np  = 7
#output file name
ofname<-paste(tg,'.csv')#,sep=',')

# panel data file cyc
for (i in 1:np){
delta<-0
dd <-0
el<-0
#input file name
ifname<-paste(i,'.txt',sep='')
a <- read.table(ifname,header=T)
#a$group <- a$Sample
a$group <- paste(a$Target,a$Sample,sep='_')
sample_stat <- aggregate(a$CT,list(a$group),function(x) c(mean(x),sd(x)))
stat_res <- data.frame(sample_stat[,2])
colnames(stat_res) <- c('mean','sd')
stat_res$group <- sample_stat[,1]
#dim(stat_res)
#head(stat_res)
#ID number
IDn<-dim(stat_res)[1]/2
#write ID into ID.txt
ID<- unique(sub("^[^_]+_","",stat_res[,3],perl=T))
#write(ID,"ID.txt")
# array CT
for (j in 1:IDn){
    #deltaT
    delta[j] <- stat_res[j,1]-stat_res[j+IDn,1]
    #deltadeltaT
    dd[j] <- delta[j]-delta[1]
    # expression level of each sample
    el[j] <- 2^(-dd[j])
}

    #integrate ID and its expression level
    fel <-data.frame(ID,el)
    #write the el to a csv file
    write.table(fel,ofname,append=TRUE,sep=",")



}
