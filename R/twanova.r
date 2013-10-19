# program to implementing two-way anova considering the interaction effect
# by JessicaBan，April 28,2013


ta<-read.table("heart.txt",header=T)

y<-ta$el    #the response variant
x1<-ta$b    #the breed variant
x2<-ta$e    #the environment variant
y.aov<-aov(y~x1+x2+x1:x2)
summary(y.aov)

