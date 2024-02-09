### mainInference.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Oct  9 2023 (10:12) 
## Version: 
## Last-Updated: feb  9 2024 (16:08) 
##           By: Brice Ozenne
##     Update #: 117
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
(coef(BTinference.H1, statistic = "favorable")-coef(BuyseTest(treatment ~ cont(score), data = dtInference[-11], trace = FALSE), statistic = "favorable"))*9
## [1] 0.04
9*(26/100 - 23/90)


#### Control group
dtInference[treatment=="C",score][1] < dtInference[treatment=="T",score]
## [1]  TRUE FALSE  TRUE FALSE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE
7/10-0.26
## [1] 0.44
getIid(BTinference.H1, statistic = "favorable", scale = FALSE, center = TRUE)[1]
## [1] 0.44

## first vs. second order
BuyseTest.options(order.Hprojection = 2)
BTinference.H2 <- BuyseTest(treatment ~ cont(score), data = dtInference, trace = FALSE)

confint(BTinference.H2, statistic = "netBenefit", order.Hprojection = 1)
##       estimate        se   lower.ci   upper.ci null    p.value
## score    -0.48 0.2216303 -0.7959335 0.04142476    0 0.06936481
confint(BTinference.H2, statistic = "netBenefit", order.Hprojection = 2)
##       estimate        se   lower.ci   upper.ci null    p.value
## score    -0.48 0.2278245 -0.8016426 0.05716099    0 0.07728498

BuyseTest.options(order.Hprojection = 1)

## * 1.4 Comparison of inferential methods
## ** 1.4.1 Confidence intervals and p-values based on asymptotic approximation

#### Example (table 1.1)
score.T <- dtInference[treatment=="T",score]
score.C <- dtInference[treatment=="C",score]

Delta.i <- dtInference[treatment=="T",.(PW = mean(score > score.C),
                                        PL = mean(score.C > score),
                                        Delta = mean(score > score.C) - mean(score.C > score)),
                       by = "id"]
Delta.j <- dtInference[treatment=="C",.(PW = mean(score.T > score),
                                        PL = mean(score > score.T),
                                        Delta = mean(score.T > score) - mean(score > score.T)),
                       by = "id"]


range(c(Delta.i$Delta,Delta.j$Delta))
## [1] -1.0  0.4

#### No transformation (first order)
Delta.Hat <- coef(BTinference.H1)

sigmaDelta.E <- mean((Delta.i$Delta-Delta.Hat)^2)
sigmaDelta.C <- mean((Delta.j$Delta-Delta.Hat)^2)

confint(BTinference.H2, statistic = "netBenefit", order.Hprojection = 1)$se^2
sigma.asym <- sigmaDelta.E/n.data + sigmaDelta.C/n.data
sigma.asym
## [1] 0.04912
c(Delta.Hat - 1.96*sqrt(sigma.asym), Delta.Hat + 1.96*sqrt(sigma.asym))
## [1] -0.91439543 -0.04560457
2*(1-pnorm(abs(Delta.Hat/sqrt(sigma.asym))))
## [1] 0.03032887

#### No transformation (second order)
PW.Hat <- coef(BTinference.H1, statistic = "favorable")
PL.Hat <- coef(BTinference.H1, statistic = "unfavorable")

sigmaPW.E <- mean((Delta.i$PW-PW.Hat)^2)
sigmaPW.C <- mean((Delta.j$PW-PW.Hat)^2)
sigmaPL.E <- mean((Delta.i$PL-PL.Hat)^2)
sigmaPL.C <- mean((Delta.j$PL-PL.Hat)^2)

sigmaPWL.E <- mean((Delta.i$PW-PW.Hat)*(Delta.i$PL-PL.Hat))
sigmaPWL.C <- mean((Delta.j$PW-PW.Hat)*(Delta.j$PL-PL.Hat))

confint(BTinference.H2, statistic = "favorable", order.Hprojection = 2)$se^2
((n.data-1)*sigmaPW.E + (n.data-1)*sigmaPW.C + PW.Hat * (1-PW.Hat))/n.data^2
## confint(BTinference.H2, statistic = "unfavorable", order.Hprojection = 2)$se^2
## (n.data-1)*sigmaPW.E/n.data^2 + (n.data-1)*sigmaPW.C/n.data^2 + PW.Hat * (1-PW.Hat)/n.data^2

confint(BTinference.H2, statistic = "netBenefit", order.Hprojection = 2)$se^2-confint(BTinference.H2, statistic = "netBenefit", order.Hprojection = 1)$se^2
(PW.Hat *(1-PW.Hat) + PL.Hat *(1-PL.Hat) + 2*PW.Hat*PL.Hat - sigmaPW.E - sigmaPW.C - sigmaPL.E - sigmaPL.C + 2*sigmaPWL.E + 2*sigmaPWL.C)/n.data^2
## [1] 0.002784

confint(BTinference.H2, statistic = "netBenefit", order.Hprojection = 2, transform = FALSE)$se^2
sigmaDelta.E/n.data + sigmaDelta.C/n.data + (PW.Hat *(1-PW.Hat) + PL.Hat *(1-PL.Hat) + 2*PW.Hat*PL.Hat - sigmaPW.E - sigmaPW.C - sigmaPL.E - sigmaPL.C + 2*sigmaPWL.E + 2*sigmaPWL.C)/n.data^2
((n.data-1) * sigmaDelta.E + (n.data-1) * sigmaDelta.C + PW.Hat *(1-PW.Hat) + PL.Hat *(1-PL.Hat) + 2*PW.Hat*PL.Hat)/n.data^2
## [1] 0.051904

#### Transformation
sigma.tasym <- sigma.asym / (1-Delta.Hat^2)^2
## [1] 0.08293317
tci <- c(atanh(Delta.Hat) - 1.96*sqrt(sigma.tasym), atanh(Delta.Hat) + 1.96*sqrt(sigma.tasym))
## [1] -1.08742741  0.04145885
tanh(tci)
## [1] -0.79593726  0.04143511
2*(1-pnorm(abs(atanh(Delta.Hat)/sqrt(sigma.tasym))))
## [1] 0.06936481

## ** 1.4.2 Bootstrap confidence intervals and p-values
BTinference.boot <- BuyseTest(treatment ~ cont(score), data = dtInference, trace = FALSE,
                              seed = 10, method.inference  = "studentized bootstrap", strata.resampling = "treatment",
                              n.resampling = 10^4)
eNBT.boot <- BTinference.boot@DeltaResampling[,"score","netBenefit"]
eNBTse.boot <- BTinference.boot@covarianceResampling[,"score","netBenefit"]
table(BTinference.boot@DeltaResampling[,"score","netBenefit"])["0"]

## *** confidence intervals (no transformation)
quantile(eNBT.boot, c(0.025,0.975))
##  2.5% 97.5% 
## -0.86  0.00 
quantile(eNBT.boot, probs = 1-c(0.0480,0.0481,0.0545,0.0596,0.0597)/2) ## close to 0
##+     97.6%   97.595%   97.275%    97.02%   97.015% 
##  0.000480  0.000000  0.000000  0.000000 -0.009403 
confint(BTinference.boot, method.ci.resampling = "percentile", transform = FALSE)
##       estimate       se lower.ci upper.ci null p.value
## score    -0.48 0.226883    -0.86        0    0  0.0545

sigma.boot <- var(eNBT.boot)
sigma.boot
## [1] 0.0514759
coef(BTinference.boot) + qnorm(c(0.025,0.975)) * sqrt(sigma.boot)
## [1] -0.92468254 -0.03531746
2*(1-pnorm(abs(coef(BTinference.boot)/sqrt(sigma.boot))))
## [1] 0.03437648
confint(BTinference.boot, method.ci.resampling = "gaussian", transform = FALSE)
##       estimate       se   lower.ci    upper.ci null    p.value
## score    -0.48 0.226883 -0.9246825 -0.03531746    0 0.03437648

t.boot <- (eNBT.boot-coef(BTinference.boot))/sqrt(eNBTse.boot)
qt.boot <- quantile(t.boot,c(0.025,0.975))
qt.boot
##      2.5%     97.5% 
## -3.457404  1.741143 
Delta.Hat + qt.boot*sqrt(sigma.asym)
##        2.5%       97.5% 
## -1.24626556 -0.09410991 
Delta.Hat + quantile(t.boot, 1-c(0.0192,0.0194,0.0195)/2) * sqrt(sigma.asym)
##       99.03% 
## 0.0003974146 
confint(BTinference.boot, method.ci.resampling = "studentized", transform = FALSE)
##       estimate        se  lower.ci    upper.ci null p.value
## score    -0.48 0.2216303 -1.246266 -0.09410991    0  0.0194


## *** confidence intervals (transformation)
eNBT.tboot <- atanh(eNBT.boot)
var(eNBT.tboot)
## [1] NaN
1.1*min(eNBT.tboot[!is.infinite(eNBT.tboot)])
## [1] -2.527316
eNBT.tboot[is.infinite(eNBT.tboot)] <- 1.1*min(eNBT.tboot[!is.infinite(eNBT.tboot)])
sigma.tboot <- var(eNBT.tboot)
sigma.tboot
## [1] 0.1162674
2*(1-pnorm(abs(atanh(coef(BTinference.boot))/sqrt(sigma.tboot))))
## [1] 0.1250868

confint(BTinference.boot, method.ci.resampling = "gaussian", transform = TRUE)
##       estimate       se   lower.ci upper.ci null   p.value
## score    -0.48 0.226883 -0.8309795  0.14431    0 0.1250868

t.tboot <- (atanh(eNBT.boot)-atanh(coef(BTinference.boot)))/sqrt(eNBTse.boot/(1-eNBT.boot^2)^2)
t.tboot[abs(eNBT.boot)==1] <- atanh(eNBT.boot)[abs(eNBT.boot)==1]
qt.tboot <- quantile(t.tboot,c(0.025,0.975), na.rm = TRUE)
qt.tboot
##      2.5%     97.5% 
## -1.825162  1.897063 

tanh(atanh(Delta.Hat) + qt.tboot*sqrt(sigma.asym/(1-Delta.Hat^2)^2))
##        2.5%       97.5% 
## -0.78126017  0.02333005 

confint(BTinference.boot, method.ci.resampling = "studentized", transform = TRUE)
##       estimate        se   lower.ci   upper.ci null p.value
## score    -0.48 0.2216303 -0.7812602 0.02333005    0  0.0617



## ** 1.4.3 Permutation p-values
BTinference.perm <- BuyseTest(treatment ~ cont(score), data = dtInference, trace = FALSE,
                              seed = 10, method.inference  = "studentized permutation", n.resampling = 10000)

confint(BTinference.perm, method = "percentile")
##       estimate        se lower.ci upper.ci null    p.value
## score    -0.48 0.2623938       NA       NA    0 0.06869313
confint(BTinference.perm, method = "studentized")
##       estimate        se lower.ci upper.ci null    p.value
## score    -0.48 0.2216303       NA       NA    0 0.05819418

var(BTinference.perm@DeltaResampling[,"score","netBenefit"])
## [1] 0.0688505
var(BTinference.boot@DeltaResampling[,"score","netBenefit"])
## [1] 0.0514759



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

## *** categorical vs continuous
res.continuousH0 <- allResS.tempo[outcome=="continuous" & mu == 0,]
res.categoricalH0 <- allResS.tempo[outcome=="categorical" & mu == 0,]
100*range(res.continuousH0$coverage-res.categoricalH0$coverage, na.rm = TRUE)
## [1] -0.1880075  0.8863975
100*range(res.continuousH0[n==200,coverage]-res.categoricalH0[n==200,coverage], na.rm = TRUE)
## [1] -0.084  0.188

res.continuousH1 <- allResS.tempo[outcome=="continuous" & mu == 1,]
res.categoricalH1 <- allResS.tempo[outcome=="categorical" & mu == 1,]
100*range(res.continuousH1$coverage-res.categoricalH1$coverage, na.rm = TRUE)
## [1] -0.368000  1.731282
100*range(res.continuousH1[n==200,coverage]-res.categoricalH1[n==200,coverage], na.rm = TRUE)
## [1] 0.208 1.432

## *** Win ratio vs Net Treatment Benefit
res.NTBH0 <- allResS.tempo[statistic=="netBenefit" & mu == 0,]
res.WRH0 <- allResS.tempo[statistic=="winRatio" & mu == 0,]
100*range(res.NTBH0$coverage-res.WRH0$coverage, na.rm = TRUE)
## [1] -0.1400000  0.9283707
100*range(res.NTBH0[n==200,coverage]-res.WRH0[n==200,coverage], na.rm = TRUE)
## [1] -0.140  0.084

res.NTBH1 <- allResS.tempo[statistic=="netBenefit" & mu == 1,]
res.WRH1 <- allResS.tempo[statistic=="winRatio" & mu == 1,]
100*range(res.NTBH1$coverage-res.WRH1$coverage, na.rm = TRUE)
## [1] -0.6244324  1.3040000
100*range(res.NTBH1[n==200,coverage]-res.WRH1[n==200,coverage], na.rm = TRUE)
## [1] -0.208  1.016

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

df.trial <- as.data.frame(list(group = c("A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B"),
                               status = c(1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,0L,0L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,0L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,0L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,0L,1L,0L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,0L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,0L,0L,1L,1L,1L,1L,1L,1L,1L,1L,1L,0L,1L,0L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,1L,0L),
                               time = c(13.9,15.7,7.5,19.8,13.8,6.6,21.4,11.5,21.4,14.3,18.2,18.6,10.3,30.3,32.4,25.8,36,9.4,16.2,51.1,21.8,17.2,33.1,8.6,15.9,9.5,21.9,11.1,59.5,13.4,18.1,29,17.8,12.6,30.1,34.2,4.6,22.2,16.2,6.7,8.6,23.7,23.1,6.9,31.1,10.9,91.5,9.6,27.4,47.2,38.1,20.3,11.8,10.5,12.4,8.1,33.3,19.6,30.8,14.5,10.7,10,8.8,26.2,8.3,6.6,16,6.6,8,31.9,22.9,4.7,6.1,6.6,46.6,81.7,14.5,27.1,30.7,15.9,7.6,9.1,4.6,9.7,23.6,24.7,7.8,14.4,15.5,9.1,4.7,17.7,32.9,21.3,35.6,17.8,36.4,5.3,26.1,17.6,13.6,19,18.4,32.5,21,39.8,25.2,36,26.8,29,6.9,31.8,20.1,17.4,58.2,33.9,54.2,50,11.4,25.5,34.6,50,22.9,4.1,19.6,10.8,48,8.4,11.3,30.2,52.8,30.6,30.5,14.2,21.9,30.6,10.6,4.9,5.5,29.6,56.4,49,42.8,60.4,16.9,17.7,10.5,36.7,37.1,44,18.7,12.7,22.8,27.5,29.3,31.4,44.6,33.4,35.7,35.8,30.8,37.8,43.7,34.9,16,17.9,22.3,40.4,5.7,5.1,24.4,11.3,20.5,37.7,93.3,68.4,16.8,21.4,6.4,4,14,30.4,9.2,6.6,98.1,4,56.7,20.1,16.2,8.2,28.7,39,17,30.8,14.1,9.7,11.2,25.7,41.7,19.8),
                               toxicity = c(0,0,1,0,0,0,0,0,0,0,1,0,1,0,0,0,0,1,1,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,1,1,0,1,1,0,0,1,1,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,1,0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,1,0,0,1,0,1,1,0,0,0,1,0,1,0,1,0,0,0,1,0,0,0,0,1,0,1,1,1,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0),
                               strata = structure(c(2L,2L,2L,1L,2L,1L,1L,1L,1L,2L,1L,1L,2L,1L,1L,1L,1L,1L,1L,1L,2L,2L,1L,1L,1L,1L,1L,1L,2L,1L,1L,1L,1L,2L,2L,1L,1L,1L,2L,1L,2L,1L,2L,2L,2L,2L,2L,2L,2L,2L,2L,1L,2L,1L,1L,2L,1L,2L,2L,2L,1L,1L,1L,2L,2L,1L,2L,1L,1L,1L,1L,1L,1L,1L,2L,1L,2L,1L,1L,2L,1L,1L,2L,2L,1L,1L,1L,2L,1L,2L,1L,2L,2L,1L,1L,1L,2L,1L,1L,2L,2L,2L,1L,1L,2L,1L,2L,2L,2L,2L,1L,2L,1L,1L,1L,2L,2L,2L,1L,2L,2L,2L,1L,1L,1L,2L,2L,2L,1L,2L,1L,1L,2L,2L,2L,1L,1L,1L,1L,2L,2L,1L,1L,2L,2L,1L,2L,2L,2L,2L,2L,2L,1L,1L,1L,2L,2L,1L,2L,1L,1L,1L,2L,1L,2L,1L,1L,1L,1L,1L,2L,1L,1L,2L,1L,1L,1L,1L,1L,1L,1L,2L,1L,1L,1L,1L,1L,1L,2L,2L,1L,2L,2L,2L,1L,1L,2L,1L,2L,2L),levels = c("1", "2"),class = "factor")))

e.mBT <- BuyseTest(group ~ cont(toxicity) + strata, data = df.trial)
confint(e.mBT, strata = TRUE)



## ** 1.5.2 Dunnett procedure for GPC

set.seed(10)
n.data <- 10
dtInference <- simBuyseTest(n.data)
dtInference[, score := round(score,1)]

score.T <- dtInference[treatment=="T",score]
score.C <- dtInference[treatment=="C",score]

## *** GPC
BTinference.H1 <- BuyseTest(treatment ~ cont(score), data = dtInference, trace = FALSE)
eSe0.BT <- sensitivity(BTinference.H1, threshold = 0:1, band = TRUE, adj.p.value = TRUE,
                       transformation = FALSE, trace = FALSE)
##   score estimate        se   lower.ci    upper.ci null      p.value
## 1     0    -0.48 0.2216303 -0.9143875 -0.04561255    0 0.0303288719
## 2     1    -0.49 0.1475805 -0.7792524 -0.20074756    0 0.0008994584
##   lower.band  upper.band adj.p.value
## 1 -0.9498640 -0.01013601  0.04475506
## 2 -0.8028758 -0.17712424  0.00147072

eSe0.BT[,c("score","estimate")]
##   score estimate
## 1     0    -0.48
## 2     1    -0.49

## *** correlation
Delta0.i <- dtInference[treatment=="T",.(Delta = mean(score > score.C) - mean(score.C > score)),
                          by = "id"]
Delta0.j <- dtInference[treatment=="C",.(Delta = mean(score.T > score) - mean(score > score.T)),
                          by = "id"]

## check
range(attr(eSe0.BT,"iid")[,1]*10 - c(Delta0.j$Delta,Delta0.i$Delta))
## [1] 0.48 0.48

## use round otherwise does not recognise that 0.3>=0.3 is TRUE probably due to rounding errors
Delta1.i <- dtInference[treatment=="T",.(Delta = mean(score >= round(score.C+1,1)) - mean(score.C >= round(score+1,1))),
                          by = "id"]
Delta1.j <- dtInference[treatment=="C",.(Delta = mean(score.T >= round(score+1,1)) - mean(score >= round(score.T+1,1))),
                          by = "id"]

## check
range(attr(eSe0.BT,"iid")[,2]*10 - c(Delta1.j$Delta,Delta1.i$Delta))
## [1] 0.49 0.49

Delta1.i$Delta
 ## [1] -0.4 -1.0 -0.4 -0.8 -0.5 -0.4 -0.4 -0.5  0.0 -0.3
Delta1.j$Delta
 ## [1]  0.0 -0.2 -0.2 -0.4 -1.0 -1.0 -0.9 -0.2 -0.8  0.0
cor(c(Delta0.i$Delta,Delta0.j$Delta),c(Delta1.i$Delta,Delta1.j$Delta))
## [1] 0.8792872
R <- cov2cor(crossprod(attr(eSe0.BT,"iid")))
R
##           [,1]      [,2]
## [1,] 1.0000000 0.8792872
## [2,] 0.8792872 1.0000000

## *** FWER
pmvnorm(lower = rep(-1.96,2), upper = rep(1.96,2), mean = c(0,0), sigma = R)
## [1] 0.9277588
## attr(,"error")
## [1] 1e-15
## attr(,"msg")
## [1] "Normal Completion"

## *** Adjusted p-values
eSe0.stats <- eSe0.BT$estimate / eSe0.BT$se
1-pmvnorm(lower = rep(eSe0.stats[1],2), upper = rep(-eSe0.stats[1],2), mean = c(0,0), sigma = R)
## [1] 0.04475506
## attr(,"error")
## [1] 1e-15
## attr(,"msg")
## [1] "Normal Completion"
1-pmvnorm(lower = rep(eSe0.stats[2],2), upper = rep(-eSe0.stats[2],2), mean = c(0,0), sigma = R)
## [1] 0.00147072
## attr(,"error")
## [1] 1e-15
## attr(,"msg")
## [1] "Normal Completion"

## *** critical threshold
qnorm.adj <- qmvnorm(0.95, mean = c(0,0), sigma = R, tail = "both.tails")
## $quantile
## [1] 2.120035

## $f.quantile
## [1] 1.471865e-05

## attr(,"message")
## [1] "Normal Completion"

cbind(lower = eSe0.BT$estimate - qnorm.adj$quantile * eSe0.BT$se,
      upper = eSe0.BT$estimate + qnorm.adj$quantile * eSe0.BT$se)
##           lower       upper
## [1,] -0.9498640 -0.01013601
## [2,] -0.8028758 -0.17712424

## *** transformation
eSe1.BT <- sensitivity(BTinference.H1, threshold = 0:1, band = TRUE, adj.p.value = TRUE,
                       transformation = TRUE, trace = FALSE)
##   score estimate        se   lower.ci    upper.ci null     p.value
## 1     0    -0.48 0.2216303 -0.7959335  0.04142476    0 0.069364812
## 2     1    -0.49 0.1475805 -0.7243353 -0.15417562    0 0.005776528
##   lower.band  upper.band adj.p.value
## 1 -0.8122186  0.08732288 0.098720663
## 2 -0.7387823 -0.12369087 0.009017036

eSe1.stats <- atanh(eSe0.BT$estimate) / (eSe0.BT$se/(1-eSe0.BT$estimate^2))
eSe1.stats
## [1] -1.816036 -2.760204
## 1-pmvnorm(lower = rep(eSe1.stats[1],2), upper = rep(-eSe1.stats[1],2), mean = c(0,0), sigma = R)
## [1] 0.09872066
## attr(,"error")
## [1] 1e-15
## attr(,"msg")
## [1] "Normal Completion"
1-pmvnorm(lower = rep(eSe1.stats[2],2), upper = rep(-eSe1.stats[2],2), mean = c(0,0), sigma = R)
## [1] 0.009017036
## attr(,"error")
## [1] 1e-15
## attr(,"msg")
## [1] "Normal Completion"
cbind(lower = tanh(atanh(eSe0.BT$estimate) - qnorm.adj$quantile * eSe0.BT$se/(1-eSe0.BT$estimate^2)),
      upper = tanh(atanh(eSe0.BT$estimate) + qnorm.adj$quantile * eSe0.BT$se/(1-eSe0.BT$estimate^2)))
##           lower       upper
## [1,] -0.8122186  0.08732288
## [2,] -0.7387823 -0.12369087

## ** 1.5.3 Finite sample performance
allResS.tempo2 <- readRDS("results/aggregated-FWER.rds")

## *** FWER and power
allResS.tempo2[mu==0, range(100*power)]
## [1] 8.832000 9.824393
allResS.tempo2[mu==0, range(100*power.bonf)]
## [1] 0.932000 1.212048
range(allResS.tempo2[mu==1, 100*(power.bonf-power)])
## [1] -35.604  -0.064
allResS.tempo2[mu==0, .(n,100*power.band)]
##      n       V2
## 1:  10 3.876155
## 2:  20 4.312000
## 3:  35 4.376000
## 4:  50 4.808000
## 5:  75 4.876000
## 6: 100 4.904000
## 7: 150 4.704000
## 8: 200 4.960000
range(allResS.tempo2[mu==1, 100*(power.band-power)])
## [1] -14.29372   0.00000

## *** coverage
allResS.tempo2[, range(100*coverage)]
## [1] 90.17561 95.16800
allResS.tempo2[mu==0, range(100*coverage.bonf)]
## [1] 98.78795 99.06800
allResS.tempo2[mu==0, range(100*coverage.band)]
## [1] 95.02400 96.13585

##----------------------------------------------------------------------
### mainInference.R ends here
