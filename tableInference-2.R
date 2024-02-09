### tableInference-2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Oct  9 2023 (10:38) 
## Version: 
## Last-Updated: feb  9 2024 (15:42) 
##           By: Brice Ozenne
##     Update #: 13
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
                      n.resampling = 1e4)
GPC.perm <- BuyseTest(treatment ~ cont(score), data = dtInference, trace = FALSE,
                      seed = 10, method.inference  = "studentized permutation", n.resampling = 1e4)


## ** inference
table3.asymp <- confint(GPC, transformation = FALSE)
table3.tasymp <- confint(GPC, transformation = TRUE)
table3.tasymp$se <- ((atanh(table3.tasymp$upper.ci) - atanh(table3.tasymp$lower.ci))/(qnorm(0.975) - qnorm(0.025)))
## table3.asymp$se/(1-table3.asymp$estimate^2)

table3.bootB <- confint(GPC.boot, method.ci = "gaussian", transformation = FALSE)
table3.bootT <- confint(GPC.boot, method.ci = "studentized", transformation = FALSE)
table3.bootT$q.upper <- (table3.bootT$upper.ci - table3.bootT$estimate)/table3.bootT$se
table3.bootT$q.lower <- (table3.bootT$lower.ci - table3.bootT$estimate)/table3.bootT$se
table3.bootP <- confint(GPC.boot, method.ci = "percentile", transformation = FALSE)

table3.tbootB <- confint(GPC.boot, method.ci = "gaussian", transformation = TRUE)
NTB.tbootB <- atanh(GPC.boot@DeltaResampling[,"score","netBenefit"])
NTB.tbootB[is.infinite(NTB.tbootB)] <- 1.1*min(NTB.tbootB[!is.infinite(NTB.tbootB)])
table3.tbootB$se <- sd(NTB.tbootB)
table3.tbootT <- confint(GPC.boot, method.ci = "studentized", transformation = TRUE)
table3.tbootT$q.upper <- (atanh(table3.tbootT$upper.ci) - atanh(table3.tbootT$estimate))/table3.tasymp$se
table3.tbootT$q.lower <- (atanh(table3.tbootT$lower.ci) - atanh(table3.tbootT$estimate))/table3.tasymp$se

table3.perm <- confint(GPC.perm, method.ci = "percentile")
table3.permT <- confint(GPC.perm, method.ci = "studentized", transformation = FALSE)
table3.tpermT <- confint(GPC.perm, method.ci = "studentized", transformation = TRUE)


table3 <- rbind(asymptotic.no = data.frame(var = signif(table3.asymp$se^2,2),
                                           CI = paste0("[",signif(table3.asymp$lower.ci,2), "; ",signif(table3.asymp$upper.ci,2),"]"),
                                           p.value = signif(table3.asymp$p.value,2)),
                asymptotic.yes = data.frame(var = signif(table3.tasymp$se^2,2),
                                            CI = paste0("[",signif(table3.tasymp$lower.ci,2), "; ",signif(table3.tasymp$upper.ci,2),"]"),
                                            p.value = signif(table3.tasymp$p.value,2)),
                basic.boot.no = data.frame(var = signif(table3.bootB$se^2,2),
                                           CI = paste0("[",signif(table3.bootB$lower.ci,3), "; ",signif(table3.bootB$upper.ci,2),"]"),
                                           p.value = signif(table3.bootB$p.value,2)),
                basic.boot.yes = data.frame(var = signif(table3.tbootB$se^2,3),
                                            CI = paste0("[",signif(table3.tbootB$lower.ci,3), "; ",signif(table3.tbootB$upper.ci,3),"]"),
                                            p.value = signif(table3.tbootB$p.value,3)),
                perc.boot = data.frame(var = signif(table3.bootP$se^2,2),
                                       CI = paste0("[",signif(table3.bootP$lower.ci,3), "; ",signif(table3.bootP$upper.ci,2),"]"),
                                       p.value = signif(table3.bootP$p.value,2)),
                stud.boot.no = data.frame(var = paste0("(",signif(table3.bootT$q.lower,3),";",signif(table3.bootT$q.upper,3),")"),
                                          CI = paste0("[",signif(table3.bootT$lower.ci,4), "; ",signif(table3.bootT$upper.ci,2),"]"),
                                          p.value = signif(table3.bootT$p.value,2)),
                stud.boot.yes = data.frame(var = paste0("(",signif(table3.tbootT$q.lower,3),";",signif(table3.tbootT$q.upper,3),")"),
                                           CI = paste0("[",signif(table3.tbootT$lower.ci,3), "; ",signif(table3.tbootT$upper.ci,2),"]"),
                                           p.value = signif(table3.tbootT$p.value,2)),
                basic.perm = data.frame(var = "",
                                        CI = "",
                                        p.value = signif(table3.perm$p.value,2)),
                stud.perm.no = data.frame(var = "",
                                       CI = "",
                                       p.value = signif(table3.permT$p.value,2)),
                stud.perm.yes = data.frame(var = "",
                                       CI = "",
                                       p.value = signif(table3.tpermT$p.value,2))
                )
                
                


## * display table 3
print(xtable(table3),include.rownames=FALSE)

## % latex table generated in R 4.2.0 by xtable 1.8-4 package
## % Fri Feb  9 15:32:52 2024
## \begin{table}[ht]
## \centering
## \begin{tabular}{llr}
##   \hline
## var & CI & p.value \\ 
##   \hline
## 0.049 & [-0.91; -0.046] & 0.03 \\ 
##   0.083 & [-0.8; 0.041] & 0.07 \\ 
##   0.051 & [-0.925; -0.035] & 0.03 \\ 
##   0.12 & [-0.831; 0.144] & 0.12 \\ 
##   0.051 & [-0.86; 0] & 0.05 \\ 
##   (-3.46;1.74) & [-1.246; -0.094] & 0.02 \\ 
##   (-1.83;1.9) & [-0.781; 0.023] & 0.06 \\ 
##    &  & 0.07 \\ 
##    &  & 0.06 \\ 
##    \hline
## \end{tabular}
## \end{table}

##----------------------------------------------------------------------
### tableInference-2.R ends here
