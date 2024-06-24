### tableInference-2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Oct  9 2023 (10:38) 
## Version: 
## Last-Updated: feb 21 2024 (17:27) 
##           By: Brice Ozenne
##     Update #: 15
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
table2.asymp <- confint(GPC, transformation = FALSE)
table2.tasymp <- confint(GPC, transformation = TRUE)
table2.tasymp$se <- ((atanh(table2.tasymp$upper.ci) - atanh(table2.tasymp$lower.ci))/(qnorm(0.975) - qnorm(0.025)))
## table2.asymp$se/(1-table2.asymp$estimate^2)

table2.bootB <- confint(GPC.boot, method.ci = "gaussian", transformation = FALSE)
table2.bootT <- confint(GPC.boot, method.ci = "studentized", transformation = FALSE)
table2.bootT$q.upper <- (table2.bootT$upper.ci - table2.bootT$estimate)/table2.bootT$se
table2.bootT$q.lower <- (table2.bootT$lower.ci - table2.bootT$estimate)/table2.bootT$se
table2.bootP <- confint(GPC.boot, method.ci = "percentile", transformation = FALSE)

table2.tbootB <- confint(GPC.boot, method.ci = "gaussian", transformation = TRUE)
NTB.tbootB <- atanh(GPC.boot@DeltaResampling[,"score","netBenefit"])
NTB.tbootB[is.infinite(NTB.tbootB)] <- 1.1*min(NTB.tbootB[!is.infinite(NTB.tbootB)])
table2.tbootB$se <- sd(NTB.tbootB)
table2.tbootT <- confint(GPC.boot, method.ci = "studentized", transformation = TRUE)
table2.tbootT$q.upper <- (atanh(table2.tbootT$upper.ci) - atanh(table2.tbootT$estimate))/table2.tasymp$se
table2.tbootT$q.lower <- (atanh(table2.tbootT$lower.ci) - atanh(table2.tbootT$estimate))/table2.tasymp$se

table2.perm <- confint(GPC.perm, method.ci = "percentile")
table2.permT <- confint(GPC.perm, method.ci = "studentized", transformation = FALSE)
table2.tpermT <- confint(GPC.perm, method.ci = "studentized", transformation = TRUE)


table2 <- rbind(asymptotic.no = data.frame(var = signif(table2.asymp$se^2,2),
                                           CI = paste0("[",signif(table2.asymp$lower.ci,2), "; ",signif(table2.asymp$upper.ci,2),"]"),
                                           p.value = signif(table2.asymp$p.value,2)),
                asymptotic.yes = data.frame(var = signif(table2.tasymp$se^2,2),
                                            CI = paste0("[",signif(table2.tasymp$lower.ci,2), "; ",signif(table2.tasymp$upper.ci,2),"]"),
                                            p.value = signif(table2.tasymp$p.value,2)),
                basic.boot.no = data.frame(var = signif(table2.bootB$se^2,2),
                                           CI = paste0("[",signif(table2.bootB$lower.ci,3), "; ",signif(table2.bootB$upper.ci,2),"]"),
                                           p.value = signif(table2.bootB$p.value,2)),
                basic.boot.yes = data.frame(var = signif(table2.tbootB$se^2,3),
                                            CI = paste0("[",signif(table2.tbootB$lower.ci,3), "; ",signif(table2.tbootB$upper.ci,3),"]"),
                                            p.value = signif(table2.tbootB$p.value,3)),
                perc.boot = data.frame(var = signif(table2.bootP$se^2,2),
                                       CI = paste0("[",signif(table2.bootP$lower.ci,3), "; ",signif(table2.bootP$upper.ci,2),"]"),
                                       p.value = signif(table2.bootP$p.value,2)),
                stud.boot.no = data.frame(var = paste0("(",signif(table2.bootT$q.lower,3),";",signif(table2.bootT$q.upper,3),")"),
                                          CI = paste0("[",signif(table2.bootT$lower.ci,4), "; ",signif(table2.bootT$upper.ci,2),"]"),
                                          p.value = signif(table2.bootT$p.value,2)),
                stud.boot.yes = data.frame(var = paste0("(",signif(table2.tbootT$q.lower,3),";",signif(table2.tbootT$q.upper,3),")"),
                                           CI = paste0("[",signif(table2.tbootT$lower.ci,3), "; ",signif(table2.tbootT$upper.ci,2),"]"),
                                           p.value = signif(table2.tbootT$p.value,2)),
                basic.perm = data.frame(var = "",
                                        CI = "",
                                        p.value = signif(table2.perm$p.value,2)),
                stud.perm.no = data.frame(var = "",
                                       CI = "",
                                       p.value = signif(table2.permT$p.value,2)),
                stud.perm.yes = data.frame(var = "",
                                       CI = "",
                                       p.value = signif(table2.tpermT$p.value,2))
                )
table2$p.value <- as.character(table2$p.value)                 
                


## * display table 3
print(xtable(table2),include.rownames=FALSE)

## \begin{table}[ht]
## \centering
## \begin{tabular}{lll}
##   \hline
## var & CI & p.value \\ 
##   \hline
## 0.049 & [-0.91; -0.046] & 0.03 \\ 
##   0.083 & [-0.8; 0.041] & 0.069 \\ 
##   0.051 & [-0.925; -0.035] & 0.034 \\ 
##   0.116 & [-0.831; 0.144] & 0.125 \\ 
##   0.051 & [-0.86; 0] & 0.054 \\ 
##   (-3.46;1.74) & [-1.246; -0.094] & 0.019 \\ 
##   (-1.83;1.9) & [-0.781; 0.023] & 0.062 \\ 
##    &  & 0.069 \\ 
##    &  & 0.061 \\ 
##    &  & 0.058 \\ 
##    \hline
## \end{tabular}
## \end{table}
##----------------------------------------------------------------------
### tableInference-2.R ends here
