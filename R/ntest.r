#function to perform normality of data, and equality of variance,
#using Shapiro-Wilk test and Bartlett's test respectively.
#x is the test data, v1 is the first factor, v2 is the second factor and so on.
#by JessicaBan, April 28,2013


###################################
##      function nev_test()      ##
###################################

"nev_test"<-function(x,v1,v2,....){
    #x<-data.frame(x)
    #normality test
    for (lev in unique(v1)){
        shapiro.test(x[v1==lev])
    }
    #equality of variance
    bartlett.test(x~v1)

    #normality test
    for (lev in unique(v2)){
        shapiro.test(x[v2==lev])
    }

    #equality of variance
    bartlett.test(x~v2)
    


}

