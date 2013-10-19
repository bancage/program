# Read values from tab-delimited autos.dat 
autos_data <-read.table("autos.dat", header=T, sep="\t") 
# Graph autos with adjacent bars using rainbow colors 
barplot(as.matrix(autos_data), main="Autos", ylab= "Total",beside=TRUE, col=rainbow(5)) 
# Place the legend at the top-left corner with no frame using rainbow colors 
legend("topleft", c("Mon","Tue","Wed","Thu","Fri"),cex=0.6, bty="n", fill=rainbow(5));