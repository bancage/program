#层次聚类及热点图
data(USArrests)
is.matrix(USArrests)
h<-hclust(dist(USArrests))
plot(h)
heatmap(as.matrix(USArrests))
