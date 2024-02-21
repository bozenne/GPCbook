### BUILD.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: sep 26 2023 (09:18) 
## Version: 
## Last-Updated: feb 12 2024 (13:31) 
##           By: Brice Ozenne
##     Update #: 66
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:
## cd /projects/biostat01/people/hpl802/GPC/book/
library(reshape2)
library(data.table)
library(ggplot2)

path <- "x:/GPC/book/"
setwd(path)

calcTruth <- function(mu, sigma, outcome, statistic){
    if(outcome=="continuous"){
        favorable <- pnorm(mu/sqrt(1+sigma^2))
        unfavorable <- pnorm(-mu/sqrt(1+sigma^2))
        neutral <- 0
    }else if(outcome=="categorical"){
        breaks <- c(-Inf,(-3):3,Inf)
        n.break <- length(breaks)
        mygrid <- expand.grid(valueT = 2:n.break - 1, valueC = 2:n.break - 1)
        p.C <- pnorm(breaks[-1],mean=0,sd=1) - pnorm(breaks[-n.break],mean=0,sd=1)
        p.T <- pnorm(breaks[-1],mean=mu,sd=sigma) - pnorm(breaks[-n.break],mean=mu,sd=sigma)
        favorable <- sum( (mygrid$valueT > mygrid$valueC) * (p.T[mygrid$valueT] * p.C[mygrid$valueC]))
        unfavorable <- sum( (mygrid$valueT < mygrid$valueC) * (p.T[mygrid$valueT] * p.C[mygrid$valueC]))
        neutral <- sum( (mygrid$valueT == mygrid$valueC) * (p.T[mygrid$valueT] * p.C[mygrid$valueC]))
    }
    out <- switch(statistic,
                  "netBenefit" = favorable - unfavorable,
                  "winRatio" = (favorable+0.5*neutral) / (unfavorable+0.5*neutral))
    return(round(out,10))    
}
VcalcTruth <- Vectorize(calcTruth)
## VcalcTruth(mu = 0:1, sigma = 2:2, outcome = c("continuous","continuous"), statistic = c("netBenefit","netBenefit"))
## VcalcTruth(mu = 0:1, sigma = 2:2, outcome = c("categorical","categorical"), statistic = c("netBenefit","netBenefit"))
## VcalcTruth(mu = 0:1, sigma = 2:2, outcome = c("continuous","continuous"), statistic = c("winRatio","winRatio"))
## VcalcTruth(mu = 0:1, sigma = 2:2, outcome = c("categorical","categorical"), statistic = c("winRatio","winRatio"))

## * Single test
path1 <- file.path("Results","sim-type1-GPC")
allRes.tempo <- butils::sinkDirectory(path1, string.keep = "tempo")

allRes.tempo[, truth := calcTruth(mu = mu, sigma = sigma, outcome = outcome, statistic = statistic), by = c("mu","sigma","outcome","statistic")]
allRes.tempo[, null := ifelse(statistic=="netBenefit", 0, 1)]
allRes.tempo[, coverage := (lower.ci<=truth)*(truth<=upper.ci)]
allRes.tempo[, mismatch := (p.value<=0.05)!=(lower.ci>null|upper.ci<null)]
allRes.tempo[, mismatch2 := (round(p.value,5)<=0.05)!=(round(lower.ci,5)>null|round(upper.ci,5)<null)]

## ** Coverage/power
allResS.tempo <- allRes.tempo[,.(rep = .N, truth = truth[1], bias = mean(estimate-truth),
                                 power = mean(p.value<=0.05),
                                 coverage = mean(coverage,na.rm = TRUE),
                                 mismatch = mean(mismatch,na.rm = TRUE),
                                 time = mean(time, na.rm=TRUE)),
                              by = c("outcome","method","statistic","mu","sigma","n")]
convertion <- c("Ustat" = "Asymptotic",
                "Ustat-trans" = "Asymptotic with transformation",
                "boot-perc" = "Percentile bootstrap",
                "boot-stud" = "Studentized bootstrap",
                "perm-perc" = "Permutation",
                "perm-stud" = "Studentized permutation")
allResS.tempo[, method.legend := factor(method, levels = names(convertion), labels = convertion)]
allResS.tempo[, mu.legend := paste0("\u0394\u03bc=",mu)]
allResS.tempo[, statistic.legend := factor(statistic, levels = c("netBenefit","winRatio"), labels = c("Net Treatment Benefit","Win Ratio"))]

if(TRUE){
    saveRDS(allResS.tempo,file = "Results/aggregated-power.rds")
}

## ** Consistency between statistics
allRes.tempoW <- data.table::dcast(allRes.tempo[,.(i,iFile,outcome,method,statistic,mu,sigma,n,p.value)],
                                   i+iFile+mu+sigma+n+outcome+method~statistic, value.var = "p.value")
allResS.tempoW <- allRes.tempoW[, .(rep = .N,
                                    mismatch = mean((netBenefit<=0.05)!=(winRatio<=0.05), na.rm = TRUE),
                                    mismatch2 = mean((round(netBenefit,5)<=0.05)!=(round(winRatio,5)<=0.05), na.rm = TRUE)),
                                by = c("n","mu","sigma","method","outcome")]
allResS.tempoW[mismatch!=0]
allResS.tempoW[mismatch2!=0]

allRes.tempoW[method == "perm-perc"][(netBenefit<=0.05)!=(winRatio<=0.05)][n==50][1]

if(TRUE){
    saveRDS(allResS.tempoW,file = "Results/aggregated-mismatch.rds")
    ## allResS.tempoW <- readRDS("Results/aggregated-mismatch.rds")
}

## * Multiple tests
path2 <- file.path("Results","sim-FWER-GPC")
allRes.tempo2 <- butils::sinkDirectory(path2, string.keep = "tempo")

allRes.tempo2[, truth := pnorm((mu-threshold)/sqrt(1+sigma^2))-pnorm((-mu-threshold)/sqrt(1+sigma^2))]
allRes.tempo2[, null := 0]
allRes.tempo2[, coverage := (lower.ci<=truth)*(truth<=upper.ci)]
allRes.tempo2[, coverage.band := (lower.band<=truth)*(truth<=upper.band)]
allRes.tempo2[, coverage.bonf := (lower.bonf<=truth)*(truth<=upper.bonf)]

allResS.tempo2 <- allRes.tempo2[,.(rep = .N, bias = mean(estimate-truth),
                                   power = mean(p.value<=0.05),
                                   power.band = mean(adj.p.value<=0.05),
                                   power.bonf = mean(bonf.p.value<=0.05),
                                   coverage = mean(coverage,na.rm = TRUE),
                                   coverage.band = mean(coverage.band,na.rm = TRUE),
                                   coverage.bonf = mean(coverage.bonf,na.rm = TRUE)
                                   ),
                                by = c("mu","sigma","n")]
allResS.tempo2
##     mu sigma   n   rep          bias      power power.band power.bonf  coverage coverage.band coverage.bonf
##  1:  0     2  10 24999  0.0028313133 0.09824393 0.03876155 0.01212048 0.9017561     0.9613585     0.9878795
##  2:  0     2  20 25000  0.0007944000 0.09536000 0.04312000 0.00988000 0.9046400     0.9568400     0.9901200
##  3:  0     2  35 25000 -0.0002356898 0.08992000 0.04376000 0.00932000 0.9100800     0.9560400     0.9906800
##  4:  0     2  50 25000  0.0002427680 0.09208000 0.04808000 0.00936000 0.9079200     0.9518000     0.9906400
##  5:  0     2  75 25000 -0.0008693973 0.09184000 0.04876000 0.01076000 0.9081600     0.9513600     0.9892400
##  6:  0     2 100 25000  0.0006308280 0.09116000 0.04904000 0.00964000 0.9088400     0.9510000     0.9903600
##  7:  0     2 150 25000  0.0004775236 0.08832000 0.04704000 0.01032000 0.9116800     0.9531200     0.9896800
##  8:  0     2 200 25000 -0.0000112840 0.09072000 0.04960000 0.00956000 0.9092800     0.9502400     0.9904400
##  9:  1     2  10 24997  0.0698470314 0.33007961 0.18714246 0.07636916 0.9021083     0.9557547     0.9859059
## 10:  1     2  20 25000  0.0509444224 0.53800000 0.39692000 0.18196000 0.9256800     0.9612400     0.9904800
## 11:  1     2  35 25000  0.0356524078 0.74996000 0.64016000 0.39512000 0.9478800     0.9728800     0.9926800
## 12:  1     2  50 25000  0.0281556368 0.87872000 0.80432000 0.59240000 0.9462400     0.9726400     0.9941600
## 13:  1     2  75 25000  0.0195517427 0.96244000 0.93460000 0.81252000 0.9458800     0.9733200     0.9949600
## 14:  1     2 100 25000  0.0161516516 0.99108000 0.98256000 0.93000000 0.9484400     0.9723600     0.9948400
## 15:  1     2 150 25000  0.0122280123 0.99968000 0.99912000 0.99196000 0.9471600     0.9706800     0.9948800
## 16:  1     2 200 25000  0.0086929485 1.00000000 1.00000000 0.99936000 0.9516800     0.9748000     0.9948800
range(allResS.tempo2[mu==1,power] - allResS.tempo2[mu==1,power.bonf])
range(allResS.tempo2[mu==1,power] - allResS.tempo2[mu==1,power.band])
range(allResS.tempo2[mu==1,coverage])
range(allResS.tempo2[mu==1,coverage.band])
range(allResS.tempo2[mu==1,coverage.bonf])


allResSL.tempo2 <- data.table::melt(allResS.tempo2[,.(n,mu,power,power.band,power.bonf)], id.vars = c("n","mu"))
allResSL.tempo2[, mu.legend := paste0("\u0394\u03bc=",mu)]
allResSL.tempo2[, variable.legend := factor(variable,
                                            levels = c("power","power.band","power.bonf"),
                                            labels = c("None","Dunnett","Bonferroni"))]

if(FALSE){
    saveRDS(allResS.tempo2,file = "Results/aggregated-FWER.rds")
}
##----------------------------------------------------------------------
### BUILD.R ends here
