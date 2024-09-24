### figureInference-3.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: sep 24 2024 (14:51) 
## Version: 
## Last-Updated: sep 24 2024 (14:54) 
##           By: Brice Ozenne
##     Update #: 4
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

## * score matrix
M.score <- (outer(X = dtInference[treatment=="T",score], Y = dtInference[treatment=="C",score], FUN = ">")-0.5)*2
colnames(M.score) <- dtInference[treatment=="C",score]
rownames(M.score) <- dtInference[treatment=="T",score]
M.score

rowSums(M.score==1)
## -0.6 -2.2 -0.7 -2.1 -1.3 -0.4 -0.7 -0.9 -0.1 -0.3 
##    3    0    3    0    0    5    3    2    5    5 
rowSums(M.score==-1)
## -0.6 -2.2 -0.7 -2.1 -1.3 -0.4 -0.7 -0.9 -0.1 -0.3 
##    7   10    7   10   10    5    7    8    5    5 
rowSums(M.score)
## -0.6 -2.2 -0.7 -2.1 -1.3 -0.4 -0.7 -0.9 -0.1 -0.3 
##   -4  -10   -4  -10  -10    0   -4   -6    0    0 
colSums(M.score)
## -1.2 -0.5 -0.8  0.3  1.1  1.2  0.7 -0.5  0.6 -1.2 
##    4   -4    2  -10  -10  -10  -10   -4  -10    4 

## * generate table 3
GPC <- BuyseTest(treatment ~ cont(score), data = dtInference, trace = FALSE, keep.pairScore = TRUE)

getPairScore(GPC)[,sum(favorable-unfavorable),by = index.T]
getPairScore(GPC)[,sum(favorable-unfavorable),by = index.C]

(rowMeans(M.score)-coef(GPC))^2
##   -0.6   -2.2   -0.7   -2.1   -1.3   -0.4   -0.7   -0.9   -0.1   -0.3 
## 0.0064 0.2704 0.0064 0.2704 0.2704 0.2304 0.0064 0.0144 0.2304 0.2304 
sum((rowMeans(M.score)-coef(GPC))^2)
## [1] 1.536
(colMeans(M.score)-coef(GPC))^2
##   -1.2   -0.5   -0.8    0.3    1.1    1.2    0.7   -0.5    0.6   -1.2 
## 0.7744 0.0064 0.4624 0.2704 0.2704 0.2704 0.2704 0.0064 0.2704 0.7744 
sum((colMeans(M.score)-coef(GPC))^2)
## [1] 3.376


##----------------------------------------------------------------------
### figureInference-3.R ends here
