# program for analyzing quantitative pcr data by JessicaBan, May 9th 2103.
a<-read.table("heart.txt",header=T)
heart.aov<-aov(a$el~as.factor(a$f))
>summary(heart.aov)
               Df Sum Sq Mean Sq F value Pr(>F)  
as.factor(a$f)  3  5.783  1.9278   2.793 0.0684 .
Residuals      19 13.116  0.6903                 
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 
> bartlett.test(a$el~a$f)

	Bartlett test of homogeneity of variances

data:  a$el by a$f 
Bartlett's K-squared = 10.4773, df = 3, p-value = 0.01492

> pairwise.t.test(a$el,as.factor(a$f),p.adj="bonf")

	Pairwise comparisons using t tests with pooled SD 

data:  a$el and as.factor(a$f) 

  1    2    3   
2 0.07 -    -   
3 1.00 0.35 -   
4 1.00 0.81 1.00

P value adjustment method: bonferroni 
> fb<-as.factor(a$b)
> fe<-as.factor(a$e)
> heart.2aov<-aov(a$el~fb+fe+fb:fe)
> summary(heart.2aov)
            Df Sum Sq Mean Sq F value Pr(>F)  
fb           1  0.275   0.275   0.398 0.5356  
fe           1  3.493   3.493   5.060 0.0365 *
fb:fe        1  2.015   2.015   2.919 0.1038  
Residuals   19 13.116   0.690                 
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 
> bartlett.test(log(a$el)~a$f)

	Bartlett test of homogeneity of variances

data:  log(a$el) by a$f 
Bartlett's K-squared = 7.4669, df = 3, p-value = 0.05842

> heart.aov<-aov(log(a$el)~as.factor(a$f))
> summary(heart.aov)
               Df Sum Sq Mean Sq F value Pr(>F)  
as.factor(a$f)  3  4.007   1.336   2.497 0.0908 .
Residuals      19 10.165   0.535                 
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 
> heart.2aov<-aov(log(a$el)~fb+fe+fb:fe)
> summary(heart.2aov)
            Df Sum Sq Mean Sq F value Pr(>F)  
fb           1  0.217   0.217   0.406 0.5318  
fe           1  1.520   1.520   2.841 0.1082  
fb:fe        1  2.270   2.270   4.243 0.0534 .
Residuals   19 10.165   0.535                 
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 
> pairwise.t.test(a$el,a$f,p.adj="bonf")

	Pairwise comparisons using t tests with pooled SD 

data:  a$el and a$f 

  1    2    3   
2 0.07 -    -   
3 1.00 0.35 -   
4 1.00 0.81 1.00

P value adjustment method: bonferroni 
> pairwise.t.test(log(a$el),a$f,p.adj="bonf")

	Pairwise comparisons using t tests with pooled SD 

data:  log(a$el) and a$f 

  1     2     3    
2 0.094 -     -    
3 1.000 0.724 -    
4 1.000 0.501 1.000

P value adjustment method: bonferroni 

> stat <- aggregate(a$el,list(a$f),function(x) c(mean(x),sd(x)))
> stat
  Group.1       x.1       x.2
1       1 1.9684100 1.1660714
2       2 0.5652801 0.1922644
3       3 1.5296199 0.9644296
4       4 1.3141884 0.7539229


