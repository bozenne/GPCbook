### mainInference.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Oct  9 2023 (10:12) 
## Version: 
## Last-Updated: Oct  9 2023 (13:07) 
##           By: Brice Ozenne
##     Update #: 59
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

options(width=75)

## * Dependencies
library(BuyseTest)
library(data.table)
library(mvtnorm)

## * 1.3 Large sample distribution
## ** 1.3.1 Intuition

set.seed(10)
n.data <- 10
dtInference <- simBuyseTest(n.data)
dtInference[, score := round(score,1)]

## table 1 --> see tableInference-1.R

## ** 1.3.2 First order H-decomposition

BTinference.H1 <- BuyseTest(treatment ~ cont(score), data = dtInference, trace = FALSE)

#### Example (table 1.1)
coef(BTinference.H1, statistic = "favorable")
## [1] 0.26

#### Experimental group
dtInference[treatment=="T",score][1]
## [1] -0.6
dtInference[treatment=="T",score][1] > dtInference[treatment=="C",score]
## [1]  TRUE FALSE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE  TRUE
3/10-0.26
## [1] 0.04
getIid(BTinference.H1, statistic = "favorable", scale = FALSE, center = TRUE)[11]
## [1] 0.04

#### Control group
dtInference[treatment=="C",score][1] < dtInference[treatment=="T",score]
## [1]  TRUE FALSE  TRUE FALSE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE
7/10-0.26
## [1] 0.44
getIid(BTinference.H1, statistic = "favorable", scale = FALSE, center = TRUE)[1]
## [1] 0.44

## * 1.4 Comparison of inferential methods
## ** 1.4.1 Confidence intervals and p-values based on asymptotic approximation

#### Example (table 1.1)
score.T <- dtInference[treatment=="T",score]
score.C <- dtInference[treatment=="C",score]

Delta.i <- dtInference.T[,.(Delta = mean(score > score.C) - mean(score.C > score)),
                         by = "id"]
Delta.j <- dtInference.C[,.(Delta = mean(score.T > score) - mean(score > score.T)),
                         by = "id"]


range(c(Delta.i$Delta,Delta.j$Delta))
## [1] -1.0  0.4

#### No transformation
sigma.asym <- mean((Delta.i$Delta-coef(GPC))^2)/n.data + mean((Delta.j$Delta-coef(GPC))^2)/n.data
## [1] 0.04912
c(coef(GPC) - 1.96*sqrt(sigma.asym), coef(GPC) + 1.96*sqrt(sigma.asym))
## [1] -0.91439543 -0.04560457
2*(1-pnorm(abs(coef(GPC)/sqrt(sigma.asym))))
## [1] 0.03032887

#### Transformation
sigma.tasym <- sigma.asym / (1-coef(GPC)^2)^2
## [1] 0.08293317
tci <- c(atanh(coef(GPC)) - 1.96*sqrt(sigma.tasym), atanh(coef(GPC)) + 1.96*sqrt(sigma.tasym))
## [1] -1.08742741  0.04145885
tanh(tci)
## [1] -0.79593726  0.04143511
2*(1-pnorm(abs(atanh(coef(GPC))/sqrt(sigma.tasym))))
## [1] 0.06936481

## ** 1.4.2 Bootstrap confidence intervals and p-values

BTinference.boot <- BuyseTest(treatment ~ cont(score), data = dtInference, trace = FALSE,
                              seed = 10, method.inference  = "studentized bootstrap", strata.resampling = "treatment",
                              n.resampling = 15)

## *** confidence intervals
quantile(BTinference.boot@DeltaResampling[,"score","netBenefit"], c(0.025,0.975))
##   2.5%  97.5% 
## -0.753  0.066 

sigma.boot <- var(BTinference.boot@DeltaResampling[,"score","netBenefit"])
## [1] 0.06346667
c(coef(GPC) - 1.96*sqrt(sigma.boot), coef(GPC) + 1.96*sqrt(sigma.boot))
## [1] -0.97377479  0.01377479

t.boot <- (BTinference.boot@DeltaResampling[,"score","netBenefit"]-coef(BTinference.boot))/sqrt(BTinference.boot@covarianceResampling[,"score","netBenefit"])
qt.boot <- quantile(t.boot,c(0.025,0.975))
##      2.5%     97.5% 
## -1.811245  1.917857 
c(coef(GPC) + qt.boot[1]*sqrt(sigma.asym), coef(GPC) + qt.boot[2]*sqrt(sigma.asym))
##        2.5%       97.5% 
## -0.88142676 -0.05494471 

## *** p-value
quantile(BTinference.boot@DeltaResampling[,"score","netBenefit"],(13:15)/15) ## close to 0
##  86.66667%  93.33333%       100% 
## -0.2893333 -0.1906667  0.2200000 

sort(t.boot)[14:15]
##        10         4 
## 0.9486272 2.4397502 
coef(GPC) + sort(t.boot)[14:15]*sqrt(sigma.asym)
##          10           4 
## -0.26975545  0.06072262 
1/(15+1)
## [1] 0.0625

## *** transformation
BTinference.boot10000 <- BuyseTest(treatment ~ cont(score), data = dtInference, trace = FALSE,
                                   seed = 10, method.inference  = "studentized bootstrap", strata.resampling = "treatment", n.resampling = 1e4)
confint(BTinference.boot10000, method.ci = "studentized", transformation = FALSE)
##       estimate        se  lower.ci    upper.ci null p.value
## score    -0.48 0.2216303 -1.246266 -0.09410991    0  0.0194

confint(BTinference.boot10000, method.ci = "studentized", transformation = TRUE)
##       estimate        se   lower.ci   upper.ci null p.value
## score    -0.48 0.2216303 -0.7812602 0.02333005    0  0.0617

## ** 1.4.3 Permutation p-values
BTinference.perm <- BuyseTest(treatment ~ cont(score), data = dtInference, trace = FALSE,
                              seed = 10, method.inference  = "studentized permutation", n.resampling = 15)

abs(BTinference.perm@DeltaResampling[,"score","netBenefit"]) > abs(coef(BTinference.perm))
##     1     2     3     4     5     6     7     8     9    10    11    12    13 
## FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE 
##    14    15 
##  TRUE FALSE
BTinference.perm@DeltaResampling[,"score","netBenefit"][14]
##   14 
## 0.59 
2/16
## [1] 0.125

BTinference.perm10000 <- BuyseTest(treatment ~ cont(score), data = dtInference, trace = FALSE,
                                   seed = 10, method.inference  = "studentized permutation", n.resampling = 10000)

var(BTinference.perm10000@DeltaResampling[,"score","netBenefit"])
## [1] 0.0688505
var(BTinference.boot10000@DeltaResampling[,"score","netBenefit"])
## [1] 0.0514759

## *** exchangeability

set.seed(4)
X <- rnorm(10, sd = 1)
Y <- rnorm(100, sd = sqrt(0.01))
index <- sample.int(110,100,replace =FALSE)
Z1 <- c(X,Y)[index]
Z2 <- c(X,Y)[-index]
c(meanX = mean(X), meanY = mean(Y), sdX = sd(X), sdY = sd(Y), diffMean = mean(X)-mean(Y))
##       meanX       meanY         sdX         sdY    diffMean 
## 0.566529289 0.001919028 1.047353007 0.090627580 0.564610261 

## *** Studentization

(mean(X)-mean(Y))/sqrt(var(X)/10 + var(Y)/100)
## [1] 1.704092

coef(BTinference.perm)/confint(BTinference.perm, transformation = FALSE)$se
## [1] -2.165769

range(BTinference.perm@DeltaResampling[,"score","netBenefit"]/sqrt(BTinference.perm@covarianceResampling[,"score","netBenefit"]))
## [1] -1.220683  2.573011

sort(abs(BTinference.perm@DeltaResampling[,"score","netBenefit"]/sqrt(BTinference.perm@covarianceResampling[,"score","netBenefit"])))[14:15]
##       15       14 
## 1.220683 2.573011 

confint(BTinference.perm10000, transformation = FALSE)
##       estimate        se lower.ci upper.ci null    p.value
## score    -0.48 0.2216303       NA       NA    0 0.06069393

## ** 1.4.4 Empirical performance
allResS.tempo <- readRDS("results/aggregated-power.rds")
allResS.tempoW <- readRDS("results/aggregated-mismatch.rds")

## *** type 1 error
allResS.tempo[statistic == "netBenefit" & n == 200 & mu==0 & method.legend == "Permutation",.(n, rep, "type 1 error" = 100*power)]
##      n   rep type 1 error
## 8: 200 25000        5.836

(0.8/100)/sqrt(0.05*0.95/25000)
## [1] 5.80381

allResS.tempo[statistic == "netBenefit" & n==10 & mu==0 & method.legend == "Asymptotic",.(n, rep, "type 1 error" = 100*power)]
##      n   rep type 1 error
## 1:  10 25000        8.944

allResS.tempo[statistic == "netBenefit" & n == 10 & mu==0 & method.legend == "Percentile bootstrap",.(n, rep, "type 1 error" = 100*power)]
##      n   rep type 1 error
## 1:  10 25000        6.536

allResS.tempo[statistic == "netBenefit" & n == 10 & mu==0 & method.legend == "Asymptotic with transformation",.(n, rep, "type 1 error" = 100*power)]
##      n   rep type 1 error
## 1:  10 25000        4.140

## *** mismatch p-value / confidence interval
allResS.tempo[!is.na(mismatch) & mismatch >0 & statistic == "netBenefit" & method.legend == "Percentile bootstrap",
              .(mu,n,mismatch = 100*mismatch)]
##     mu   n mismatch
##  1:  0  10    0.620
##  2:  0  20    0.068
##  3:  0  50    0.012
##  4:  0 100    0.012
##  5:  0 150    0.020
##  6:  0 200    0.016
##  7:  1  10    1.556
##  8:  1  20    0.224
##  9:  1  50    0.080
## 10:  1 100    0.004
allRes.tempo[method=="boot-perc" & n==10 & mismatch>0 & mu == 1 & statistic=="netBenefit",
             unique(p.value)]
## [1] 0.0467 0.0311
allRes.tempo[method=="boot-perc" & n==10 & mismatch>0 & mu == 1 & statistic=="netBenefit",
             .(table(lower.ci))]
##    lower.ci   N
## 1:     -0.9   2
## 2:    -0.88   1
## 3:        0 291a
allResS.tempo[!is.na(mismatch) & mismatch >0 & statistic == "netBenefit" & method.legend == "Studentized bootstrap",
              .(mu,n,mismatch = 100*mismatch)]
##    mu   n mismatch
## 1:  0  10    0.012
## 2:  0  20    0.004
## 3:  0  75    0.004
## 4:  0 100    0.004
## 5:  0 150    0.016
## 6:  1  10    0.012
## 7:  1  20    0.016
## 8:  1  35    0.016
## 9:  1  50    0.008


## *** mismatch net benefit / win ratio
table.mismatch <- allResS.tempoW[mismatch>0]

unique(table.mismatch$method)
## [1] "Ustat"     "perm-perc" "perm-stud"

table.mismatch[method == "Ustat" & n %in% c(10,50,200) & mu == 0]
##      n mu sigma method   rep mismatch
## 1:  10  0     2  Ustat 25000  0.09400
## 2:  50  0     2  Ustat 25000  0.05164
## 3: 200  0     2  Ustat 25000  0.02568

table.mismatch[method != "Ustat" & n %in% c(10,50,200) & mu == 0]
##      n mu sigma    method   rep mismatch
## 1:  50  0     2 perm-perc 25000    2e-04
## 2: 200  0     2 perm-perc 25000    4e-05
## 3:  10  0     2 perm-stud 25000    2e-04

## *** time
dt.time <- dcast(allResS.tempo[,.(time=round(median(time),2)),by=c("method","n")],
                 method~n, value.var = "time")
dt.time
##         method   10   20   35   50   75   100   150   200
## 1:       Ustat 0.00 0.00 0.00 0.01 0.01  0.01  0.01  0.01
## 2: Ustat-trans 0.00 0.00 0.00 0.01 0.01  0.01  0.01  0.01
## 3:   boot-perc   NA   NA   NA   NA   NA    NA    NA    NA
## 4:   boot-stud 3.09 3.46 4.23 5.36 8.04 11.62 21.71 35.80
## 5:   perm-perc   NA   NA   NA   NA   NA    NA    NA    NA
## 6:   perm-stud 3.61 3.98 4.79 5.88 8.58 12.16 22.36 36.32

## * 1.5 Adjustment for multiple testing
## ** 1.5.2 Dunnett procedure for GPC

set.seed(10)
n.data <- 10
dtInference <- simBuyseTest(n.data)
dtInference[, score := round(score,1)]

score.T <- round(dtInference[treatment=="T",score],1)
score.C <- round(dtInference[treatment=="C",score],1)

## GPC
BTinference.H1 <- BuyseTest(treatment ~ cont(score), data = dtInference, trace = FALSE)
eSe0.BT <- sensitivity(BTinference.H1, threshold = 0:1, band = TRUE, adj.p.value = TRUE,
                       transformation = FALSE, trace = FALSE)
eSe0.BT[,c("score","estimate")]
##   score estimate
## 1     0    -0.48
## 2     1    -0.49

Delta0.i <- dtInference.T[,.(Delta = mean(score > score.C) - mean(score.C > score)),
                          by = "id"]
Delta0.j <- dtInference.C[,.(Delta = mean(score.T > score) - mean(score > score.T)),
                          by = "id"]

## check
range(attr(eSe0.BT,"iid")[,1]*10 - c(Delta0.j$Delta,Delta0.i$Delta))
## [1] 0.48 0.48

Delta1.i <- dtInference.T[,.(Delta = mean(score >= round(score.C+1,1)) - mean(score.C >= round(score+1,1))),
                          by = "id"]
Delta1.j <- dtInference.C[,.(Delta = mean(score.T >= round(score+1,1)) - mean(score >= round(score.T+1,1))),
                          by = "id"]

## check
range(attr(eSe0.BT,"iid")[,2]*10 - c(Delta1.j$Delta,Delta1.i$Delta))
## [1] 0.49 0.49

Delta1.i$Delta
 ## [1] -0.4 -1.0 -0.4 -0.8 -0.5 -0.4 -0.4 -0.5  0.0 -0.3
Delta1.j$Delta
 ## [1]  0.0 -0.2 -0.2 -0.4 -1.0 -1.0 -0.9 -0.2 -0.8  0.0
cor(c(Delta0.i$Delta,Delta0.j$Delta),c(Delta1.i$Delta,Delta1.j$Delta))
cor(c(Delta0.i$Delta,Delta0.j$Delta),c(Delta1.i$Delta,Delta1.j$Delta))
## [1] 0.8516011
cov2cor(crossprod(attr(eSe0.BT,"iid")))
##           [,1]      [,2]
## [1,] 1.0000000 0.8792872
## [2,] 0.8792872 1.0000000

attr(eSe0.BT,"iid")[,2] / c(Delta0.i$Delta + 0.48,Delta0.j$Delta + 0.48)

getIid(BTinference.H1, center = FALSE)*10 - c(Delta0.j$Delta,Delta0.i$Delta)
attr(eSe0.BT,"iid")[,2]*10 - c(Delta1.j$Delta,Delta1.i$Delta)


## ** 1.5.3 Finite sample performance


##----------------------------------------------------------------------
### mainInference.R ends here
