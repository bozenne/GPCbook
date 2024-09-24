### tableInference-2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Oct  9 2023 (10:38) 
## Version: 
## Last-Updated: sep 24 2024 (14:12) 
##           By: Brice Ozenne
##     Update #: 19
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

## * warper to generate table
dataTOtable2 <- function(data){

    ## ** run GPC
    GPC <- BuyseTest(treatment ~ cont(score), data = data, trace = FALSE)
    GPC.boot <- BuyseTest(treatment ~ cont(score), data = data, trace = FALSE,
                          seed = 10, method.inference  = "studentized bootstrap", strata.resampling = "treatment",
                          n.resampling = 1e4)
    GPC.perm <- BuyseTest(treatment ~ cont(score), data = data, trace = FALSE,
                          seed = 10, method.inference  = "studentized permutation", n.resampling = 1e4)


    ## ** inference
    out.asymp <- confint(GPC, transformation = FALSE)
    out.tasymp <- confint(GPC, transformation = TRUE)
    out.tasymp$se <- ((atanh(out.tasymp$upper.ci) - atanh(out.tasymp$lower.ci))/(qnorm(0.975) - qnorm(0.025)))
    ## out.asymp$se/(1-out.asymp$estimate^2)

    out.bootB <- confint(GPC.boot, method.ci = "gaussian", transformation = FALSE)
    out.bootT <- confint(GPC.boot, method.ci = "studentized", transformation = FALSE)
    out.bootT$q.upper <- (out.bootT$upper.ci - out.bootT$estimate)/out.bootT$se
    out.bootT$q.lower <- (out.bootT$lower.ci - out.bootT$estimate)/out.bootT$se
    out.bootP <- confint(GPC.boot, method.ci = "percentile", transformation = FALSE)

    out.tbootB <- confint(GPC.boot, method.ci = "gaussian", transformation = TRUE)
    NTB.tbootB <- atanh(GPC.boot@DeltaResampling[,"score","netBenefit"])
    NTB.tbootB[is.infinite(NTB.tbootB)] <- 1.1*min(NTB.tbootB[!is.infinite(NTB.tbootB)])
    out.tbootB$se <- sd(NTB.tbootB)
    out.tbootT <- confint(GPC.boot, method.ci = "studentized", transformation = TRUE)
    out.tbootT$q.upper <- (atanh(out.tbootT$upper.ci) - atanh(out.tbootT$estimate))/out.tasymp$se
    out.tbootT$q.lower <- (atanh(out.tbootT$lower.ci) - atanh(out.tbootT$estimate))/out.tasymp$se

    out.perm <- confint(GPC.perm, method.ci = "percentile")
    out.permT <- confint(GPC.perm, method.ci = "studentized", transformation = FALSE)
    out.tpermT <- confint(GPC.perm, method.ci = "studentized", transformation = TRUE)


    out <- rbind(asymptotic.no = data.frame(var = signif(out.asymp$se^2,2),
                                               CI = paste0("[",signif(out.asymp$lower.ci,2), "; ",signif(out.asymp$upper.ci,2),"]"),
                                               p.value = signif(out.asymp$p.value,2)),
                    asymptotic.yes = data.frame(var = signif(out.tasymp$se^2,2),
                                                CI = paste0("[",signif(out.tasymp$lower.ci,2), "; ",signif(out.tasymp$upper.ci,2),"]"),
                                                p.value = signif(out.tasymp$p.value,2)),
                    basic.boot.no = data.frame(var = signif(out.bootB$se^2,2),
                                               CI = paste0("[",signif(out.bootB$lower.ci,3), "; ",signif(out.bootB$upper.ci,2),"]"),
                                               p.value = signif(out.bootB$p.value,2)),
                    basic.boot.yes = data.frame(var = signif(out.tbootB$se^2,3),
                                                CI = paste0("[",signif(out.tbootB$lower.ci,3), "; ",signif(out.tbootB$upper.ci,3),"]"),
                                                p.value = signif(out.tbootB$p.value,3)),
                    perc.boot = data.frame(var = signif(out.bootP$se^2,2),
                                           CI = paste0("[",signif(out.bootP$lower.ci,3), "; ",signif(out.bootP$upper.ci,2),"]"),
                                           p.value = signif(out.bootP$p.value,2)),
                    stud.boot.no = data.frame(var = paste0("(",signif(out.bootT$q.lower,3),";",signif(out.bootT$q.upper,3),")"),
                                              CI = paste0("[",signif(out.bootT$lower.ci,4), "; ",signif(out.bootT$upper.ci,2),"]"),
                                              p.value = signif(out.bootT$p.value,2)),
                    stud.boot.yes = data.frame(var = paste0("(",signif(out.tbootT$q.lower,3),";",signif(out.tbootT$q.upper,3),")"),
                                               CI = paste0("[",signif(out.tbootT$lower.ci,3), "; ",signif(out.tbootT$upper.ci,2),"]"),
                                               p.value = signif(out.tbootT$p.value,2)),
                    basic.perm = data.frame(var = "",
                                            CI = "",
                                            p.value = signif(out.perm$p.value,2)),
                    stud.perm.no = data.frame(var = "",
                                              CI = "",
                                              p.value = signif(out.permT$p.value,2)),
                    stud.perm.yes = data.frame(var = "",
                                               CI = "",
                                               p.value = signif(out.tpermT$p.value,2))
                    )
    out$p.value <- as.character(out$p.value)

    ## ** export
    return(out)

}

                
## * table 2 (in black)
set.seed(10)
n.data <- 10
dtInference10 <- simBuyseTest(n.data)
dtInference10[, score := round(score,1)]
table2.10 <- dataTOtable2(dtInference10)
##                         var               CI p.value
## asymptotic.no         0.049  [-0.91; -0.046]    0.03
## asymptotic.yes        0.083    [-0.8; 0.041]   0.069
## basic.boot.no         0.051 [-0.925; -0.035]   0.034
## basic.boot.yes        0.116  [-0.831; 0.144]   0.125
## perc.boot             0.051       [-0.86; 0]   0.054
## stud.boot.no   (-3.46;1.74) [-1.246; -0.094]   0.019
## stud.boot.yes   (-1.83;1.9)  [-0.781; 0.023]   0.062
## basic.perm                                     0.069
## stud.perm.no                                   0.061
## stud.perm.yes                                  0.058

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

## * table 2 (in grey)
set.seed(1)
n.data <- 100
dtInference100 <- simBuyseTest(n.data)
dtInference100[, score := round(score,2)]
table2.100 <- dataTOtable2(dtInference100)
##                         var               CI p.value
## asymptotic.no        0.0065    [-0.019; 0.3]   0.083
## asymptotic.yes       0.0068   [-0.021; 0.29]   0.088
## basic.boot.no        0.0066   [-0.0194; 0.3]   0.085
## basic.boot.yes      0.00697 [-0.0227; 0.295]  0.0914
## perc.boot            0.0066   [-0.0212; 0.3]   0.089
## stud.boot.no      (-1.97;2)  [-0.01896; 0.3]   0.084
## stud.boot.yes  (-1.98;1.92)  [-0.0221; 0.29]   0.091
## basic.perm                                     0.085
## stud.perm.no                                   0.085
## stud.perm.yes                                  0.085

##----------------------------------------------------------------------
### tableInference-2.R ends here
