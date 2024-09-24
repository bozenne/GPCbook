### figureInference-1.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: sep 24 2024 (14:14) 
## Version: 
## Last-Updated: sep 24 2024 (14:35) 
##           By: Brice Ozenne
##     Update #: 5
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

library(BuyseTest)
library(xtable)
digits <- c(1,3)

## * simulate data
set.seed(10)
n.data <- 10
dtInference <- simBuyseTest(n.data)
dtInference[, score := round(score,1)]

##        id treatment  eventtime status toxicity score
##     <int>    <fctr>      <num>  <num>   <fctr> <num>
##  1:     1         C 0.12835260      0      yes  -1.2
##  2:     2         C 0.72453892      1      yes  -0.5
##  3:     3         C 1.28655942      1       no  -0.8
##  4:     4         C 0.04001154      1       no   0.3
##  5:     5         C 0.95792484      1       no   1.1
##  6:     6         C 0.09777004      1       no   1.2
##  7:     7         C 0.44152664      1      yes   0.7
##  8:     8         C 0.24508634      1       no  -0.5
##  9:     9         C 0.06606621      0      yes   0.6
## 10:    10         C 0.05038985      1       no  -1.2
## 11:    11         T 1.72276031      1       no  -0.6
## 12:    12         T 0.62656054      0       no  -2.2
## 13:    13         T 0.37864505      1      yes  -0.7
## 14:    14         T 0.46192620      1       no  -2.1
## 15:    15         T 0.09098929      1       no  -1.3
## 16:    16         T 0.13648098      1       no  -0.4
## 17:    17         T 0.27812384      1      yes  -0.7
## 18:    18         T 0.03946623      0      yes  -0.9
## 19:    19         T 1.23964915      1       no  -0.1
## 20:    20         T 0.75902556      1       no  -0.3

## * score matrix
M.score <- (outer(X = dtInference[treatment=="T",score], Y = dtInference[treatment=="C",score], FUN = ">")-0.5)*2
colnames(M.score) <- dtInference[treatment=="C",score]
rownames(M.score) <- dtInference[treatment=="T",score]
M.score
##      -1.2 -0.5 -0.8 0.3 1.1 1.2 0.7 -0.5 0.6 -1.2
## -0.6    1   -1    1  -1  -1  -1  -1   -1  -1    1
## -2.2   -1   -1   -1  -1  -1  -1  -1   -1  -1   -1
## -0.7    1   -1    1  -1  -1  -1  -1   -1  -1    1
## -2.1   -1   -1   -1  -1  -1  -1  -1   -1  -1   -1
## -1.3   -1   -1   -1  -1  -1  -1  -1   -1  -1   -1
## -0.4    1    1    1  -1  -1  -1  -1    1  -1    1
## -0.7    1   -1    1  -1  -1  -1  -1   -1  -1    1
## -0.9    1   -1   -1  -1  -1  -1  -1   -1  -1    1
## -0.1    1    1    1  -1  -1  -1  -1    1  -1    1
## -0.3    1    1    1  -1  -1  -1  -1    1  -1    1

## * generate table 1
GPC <- BuyseTest(treatment ~ cont(score), data = dtInference, trace = FALSE)
rowSums(M.score==1)
## -0.6 -2.2 -0.7 -2.1 -1.3 -0.4 -0.7 -0.9 -0.1 -0.3 
##    3    0    3    0    0    5    3    2    5    5 
sum(M.score==1) ## same as coef(GPC, statistic = "count.favorable")
## [1] 26
getIid(GPC, statistic = "favorable", scale = FALSE)[dtInference$treatment=="T"]
 ## [1]  0.04 -0.26  0.04 -0.26 -0.26  0.24  0.04 -0.06  0.24  0.24

##----------------------------------------------------------------------
### figureInference-1.R ends here
