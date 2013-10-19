cars<-c(1,3,6,4,9)
trucks<-c(2,5,4,5,12)
#g_range<-range(0,cars,trucks)
par(mar=c(6,6,3,6))
#plot(cars,type="o",col="blue",main="Autos",ylim=g_range,axes=F,ann=F)
plot(cars,type="o",col="blue",main="Autos",axes=F,ann=F)
axis(1,at=1:5,lab=c("Mon","Tue","Wed","Thu","Fri"))
axis(4,las=1,at=4*0:g_range[2])
box()
par(new=T)
plot(trucks,type="o",col="red",pch=22,lty=2,axes=F,ann=F)
axis(2)
legend("topleft",c("cars","trucks"),cex=0.8,col=c("blue","red"),pch=21:22,lty=1:2)
mtext("cars",side=4,srt=90,line=2)
#¼Ó×ó²à×Ý×ø±êÍ¼Àý
title(ylab="trucks")   


