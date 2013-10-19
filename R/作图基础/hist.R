autos_data <-read.table("autos.dat", header=T, sep="\t")
autos <-c(autos_data$cars,autos_data$trucks,autos_data$suvs) 
# Compute the largest y value used in the autos 
max_num<-max(autos)
hist(autos, col=heat.colors(max_num),breaks=max_num, xlim=c(0,max_num), right=F,main="Autos Histogram",las=1)
brk <-c(0,3,4,5,6,10,16)
hist(autos, col=heat.colors(length(brk)), breaks=brk,xlim=c(0,max_num), right=F, main="Probability Density", las=1, cex.axis=0.8, freq=F)