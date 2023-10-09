### tableInference-3.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Oct  9 2023 (10:38) 
## Version: 
## Last-Updated: Oct  9 2023 (11:40) 
##           By: Brice Ozenne
##     Update #: 8
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

## * generate table 3

## ** run GPC
GPC <- BuyseTest(treatment ~ cont(score), data = dtInference, trace = FALSE)
GPC.boot <- BuyseTest(treatment ~ cont(score), data = dtInference, trace = FALSE,
                      seed = 10, method.inference  = "studentized bootstrap", strata.resampling = "treatment",
                      n.resampling = 15)
GPC.perm <- BuyseTest(treatment ~ cont(score), data = dtInference, trace = FALSE,
                      seed = 10, method.inference  = "studentized permutation", n.resampling = 15)

## ** extract for table
Delta.i <- dtInference[treatment=="T",
                           .(Delta = mean(score > dtInference[treatment=="C",score]) - mean(score < dtInference[treatment=="C",score])),
                           by = "id"]
Delta.j <- dtInference[treatment=="C",
                           .(Delta = mean(dtInference[treatment=="T",score] > score) - mean(dtInference[treatment=="T",score] < score)),
                           by = "id"]

table3.inference <- data.frame(deltai = c(formatC(Delta.i$Delta, format = "f", digits = 1), rep(NA,5)),
                               deltaj = c(formatC(Delta.j$Delta, format = "f", digits = 1), rep(NA,5)),
                               s1 = NA,
                               delta.boot = GPC.boot@DeltaResampling[,"score","netBenefit"],
                               se.boot = round(GPC.boot@covarianceResampling[,"score","netBenefit"],3),
                               s2 = NA,
                               delta.perm = GPC.perm@DeltaResampling[,"score","netBenefit"],
                               se.perm = round(GPC.perm@covarianceResampling[,"score","netBenefit"],3))

## ** inference
table3.asymp <- confint(GPC, transformation = FALSE)
table3.tasymp <- confint(GPC, transformation = TRUE)
table3.tasymp$se <- ((atanh(table3.tasymp$upper.ci) - atanh(table3.tasymp$lower.ci))/(qnorm(0.975) - qnorm(0.025)))

table3.bootB <- confint(GPC.boot, method.ci = "gaussian", transformation = FALSE)
table3.bootT <- confint(GPC.boot, method.ci = "studentized", transformation = FALSE)
table3.bootT$q.upper <- (table3.bootT$upper.ci - table3.bootT$estimate)/table3.bootT$se
table3.bootT$q.lower <- (table3.bootT$lower.ci - table3.bootT$estimate)/table3.bootT$se
table3.bootP <- confint(GPC.boot, method.ci = "percentile", transformation = FALSE)

table3.perm <- confint(GPC.perm, method.ci = "percentile", transformation = FALSE)
table3.permT <- confint(GPC.perm, method.ci = "studentized", transformation = FALSE)

attr(table3.inference, "additional") <- list(
    "asymptotic" = data.frame(var = round(table3.asymp$se^2,4),
                              lower.ci = round(table3.asymp$lower.ci,3),
                              upper.ci = round(table3.asymp$upper.ci,3),
                              p = round(table3.asymp$p.value,3)),
    "transformed asymptotic" = data.frame(var = round(table3.tasymp$se^2,4),
                                          lower.ci = round(table3.tasymp$lower.ci,3),
                                          upper.ci = round(table3.tasymp$upper.ci,3),
                                          p = round(table3.tasymp$p.value,3)),
    "basic bootstrap" = data.frame(var = round(table3.bootB$se^2,4),
                                   lower.ci = round(table3.bootB$lower.ci,3),
                                   upper.ci = round(table3.bootB$upper.ci,3),
                                   p = round(table3.bootB$p.value,3)),
    "studentized bootstrap" = data.frame(q.upper = round(table3.bootT$q.upper,3),
                                         q.lower = round(table3.bootT$q.lower,3),
                                         lower.ci = round(table3.bootT$lower.ci,3),
                                         upper.ci = round(table3.bootT$upper.ci,3),
                                         p = round(table3.bootT$p.value,3)),
    "percentile bootstrap" = data.frame(lower.ci = round(table3.bootP$lower.ci,3),
                                        upper.ci = round(table3.bootP$upper.ci,3),
                                        p = round(table3.bootP$p.value,3)),
    "permutation" = data.frame(p = round(table3.perm$p.value,3)),
    "studentized permutation" = data.frame(p = round(table3.permT$p.value,3))
)

## * display table 3
print(xtable(table3.inference, digits = 3),include.rownames=FALSE)

## % latex table generated in R 4.2.0 by xtable 1.8-4 package
## % Mon Oct  9 10:45:01 2023
## \begin{table}[ht]
## \centering
## \begin{tabular}{lllrrlrr}
##   \hline
## deltai & deltaj & s1 & delta.boot & se.boot & s2 & delta.perm & se.perm \\ 
##   \hline
## -0.4 & 0.4 &  & -0.300 & 0.066 &  & -0.080 & 0.076 \\ 
##   -1.0 & -0.4 &  & -0.520 & 0.050 &  & 0.200 & 0.068 \\ 
##   -0.4 & 0.2 &  & -0.540 & 0.044 &  & -0.010 & 0.072 \\ 
##   -1.0 & -1.0 &  & 0.220 & 0.082 &  & -0.080 & 0.074 \\ 
##   -1.0 & -1.0 &  & -0.700 & 0.028 &  & 0.210 & 0.065 \\ 
##   0.0 & -1.0 &  & -0.640 & 0.039 &  & -0.100 & 0.068 \\ 
##   -0.4 & -1.0 &  & -0.420 & 0.070 &  & 0.190 & 0.067 \\ 
##   -0.6 & -0.4 &  & -0.360 & 0.062 &  & 0.050 & 0.069 \\ 
##   0.0 & -1.0 &  & -0.700 & 0.040 &  & -0.040 & 0.073 \\ 
##   0.0 & 0.4 &  & -0.220 & 0.075 &  & 0.300 & 0.068 \\ 
##    &  &  & -0.760 & 0.022 &  & 0.110 & 0.071 \\ 
##    &  &  & -0.740 & 0.024 &  & 0.090 & 0.067 \\ 
##    &  &  & -0.540 & 0.044 &  & -0.190 & 0.069 \\ 
##    &  &  & -0.420 & 0.055 &  & 0.590 & 0.053 \\ 
##    &  &  & -0.460 & 0.050 &  & -0.300 & 0.060 \\ 
##    \hline
## \end{tabular}
## \end{table}

attr(table3.inference, "additional")
## $asymptotic
##      var lower.ci upper.ci    p
## 1 0.0491   -0.914   -0.046 0.03

## $`transformed asymptotic`
##      var lower.ci upper.ci     p
## 1 0.0829   -0.796    0.041 0.069

## $`basic bootstrap`
##      var lower.ci upper.ci     p
## 1 0.0635   -0.974    0.014 0.057

## $`studentized bootstrap`
##   q.upper q.lower lower.ci upper.ci     p
## 1   1.918  -1.811   -0.881   -0.055 0.062

## $`percentile bootstrap`
##   lower.ci upper.ci     p
## 1   -0.753    0.066 0.067

## $permutation
##       p
## 1 0.125

## $`studentized permutation`
##       p
## 1 0.125

##----------------------------------------------------------------------
### tableInference-3.R ends here
