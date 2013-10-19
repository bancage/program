# program to output basic statistics by JessicaBan,July 12th 2013.
a<-read.csv("./Pheno_Brown.csv",header=T)
pheno<-a[5:18]
Mean<-mean(pheno)
SD<-sd(pheno)
Max<-sapply(pheno,max)
Min<-sapply(pheno,min)
stat<-data.frame(Mean,SD,Max,Min)
write.csv(stat,"stat.csv")