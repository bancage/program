f# Read values from tab-delimited autos.dat 
autos_data <-read.table("autos.dat", header=T, sep="\t")
# Expand right side of clipping rect to make room for the legend 
par(xpd=T, mar=par()$mar+c(0,0,0,4))
# Graph autos (transposing the matrix) using heat colors, 
# put 10% of the space between each bar, and make labels 
# smaller with horizontal y-axis labels 
barplot(t(autos_data), main="Autos", ylab="Total", col=heat.colors(3), space=0.1, 
cex.axis=0.8, las=1, names.arg=c("Mon","Tue","Wed","Thu","Fri"), cex=0.8) 
# Place the legend at (6,30) using heat colors 
legend(6, 30,names(autos_data), cex=0.8, fill=heat.colors(3));
# Restore default clipping rect 
par(mar=c(5, 4, 4, 2) + 0.1)
