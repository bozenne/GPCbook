### BUILD.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: sep 26 2023 (09:18) 
## Version: 
## Last-Updated: okt  4 2023 (19:19) 
##           By: Brice Ozenne
##     Update #: 35
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

## * Single test
path1 <- file.path("Results","sim-type1-GPC")
allRes.tempo <- butils::sinkDirectory(path1, string.keep = "tempo")

allRes.tempo[statistic == "netBenefit", truth := pnorm(mu/sqrt(1+sigma^2))-pnorm(-mu/sqrt(1+sigma^2))]
allRes.tempo[, null := sapply(statistic, switch, netBenefit = 0, winRatio = 1)]
allRes.tempo[statistic == "winRatio", truth := pnorm(mu/sqrt(1+sigma^2))/pnorm(-mu/sqrt(1+sigma^2))]
allRes.tempo[, coverage := (lower.ci<=truth)*(truth<=upper.ci)]
allRes.tempo[, mismatch := (p.value<=0.05)!=(lower.ci>null|upper.ci<null)]

## error


## ** Coverage/power
allResS.tempo <- allRes.tempo[,.(rep = .N, truth = truth[1], bias = mean(estimate-truth),
                                 power = mean(p.value<=0.05),
                                 coverage = mean(coverage,na.rm = TRUE),
                                 mismatch = mean(mismatch,na.rm = TRUE),
                                 time = mean(time, na.rm=TRUE)),
                              by = c("method","statistic","mu","sigma","n")]
convertion <- c("Ustat" = "Asymptotic",
                "Ustat-trans" = "Asymptotic with transformation",
                "boot-perc" = "Percentile bootstrap",
                "boot-stud" = "Studentized bootstrap",
                "perm-perc" = "Permutation",
                "perm-stud" = "Studentized permutation")
allResS.tempo[, method.legend := factor(method, levels = names(convertion), labels = convertion)]
allResS.tempo[, mu.legend := paste0("\u0394\u03bc=",mu)]
allResS.tempo[, statistic.legend := factor(statistic, levels = c("netBenefit","winRatio"), labels = c("Net Treatment Benefit","Win Ratio"))]

if(FALSE){
    saveRDS(allResS.tempo,file = "Results/aggregated-power.rds")
}

## ** Consistency between statistics
allRes.tempoW <- dcast(allRes.tempo[,.(i,iFile,method,statistic,mu,sigma,n,p.value)],
                       i+iFile+mu+sigma+n+method~statistic, value.var = "p.value")
allResS.tempoW <- allRes.tempoW[, .(rep = .N, mismatch = mean((netBenefit<=0.05)!=(winRatio<=0.05), na.rm = TRUE)), by = c("n","mu","sigma","method")]
allResS.tempoW[mismatch!=0]

if(FALSE){
    saveRDS(allResS.tempoW,file = "Results/aggregated-mismatch.rds")
}

## ** Results

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
allResS.tempoW <- allRes.tempoW[, .(.N,
                                    powerNB = mean(netBenefit<=0.05),
                                    powerWR = mean(winRatio<=0.05),
                                    mismatch = 100*mean((netBenefit<=0.05)!=(winRatio<=0.05), na.rm = TRUE)),
              by = c("n","mu","method")]
allResS.tempoW[mismatch>0]
##       n mu    method     N powerNB powerWR    mismatch
##  1:  10  0     Ustat 25000 0.08944 0.09432  9.40000000
##  2:  10  0 perm-stud 25000 0.05128 0.05132  0.02000000
##  3:  20  0     Ustat 25000 0.06936 0.07456  7.39200000
##  4:  20  0 perm-perc 25000 0.05724 0.05792  0.06800000
##  5:  35  0     Ustat 25000 0.05740 0.06016  5.79600000
##  6:  35  0 perm-perc 25000 0.05568 0.05588  0.02000000
##  7:  50  0     Ustat 25000 0.05824 0.05628  5.16400000
##  8:  50  0 perm-perc 25000 0.05916 0.05936  0.02000000
##  9:  75  0     Ustat 25000 0.05668 0.05636  4.31200000
## 10:  75  0 perm-perc 25000 0.06120 0.06124  0.00400000
## 11: 100  0     Ustat 25000 0.05088 0.05028  3.59600000
## 12: 100  0 perm-perc 25000 0.05652 0.05656  0.00400000
## 13: 150  0     Ustat 25000 0.05364 0.05052  3.09600000
## 14: 200  0     Ustat 25000 0.05108 0.05092  2.56800000
## 15: 200  0 perm-perc 25000 0.05836 0.05840  0.00400000
## 16:  10  1     Ustat 25000 0.33300      NA 33.38536039
## 17:  10  1 perm-stud 25000 0.23868      NA  0.02401249
## 18:  20  1     Ustat 25000 0.50640 0.00044 50.66800000
## 19:  20  1 perm-perc 25000 0.47476 0.47492  0.01600000
## 20:  20  1 perm-stud 25000 0.45176 0.45180  0.00400000
## 21:  35  1     Ustat 25000 0.71244 0.22032 49.21200000
## 22:  35  1 perm-perc 25000 0.70836 0.70848  0.01200000
## 23:  50  1     Ustat 25000 0.84636 0.55680 28.95600000
## 24:  50  1 perm-perc 25000 0.84920 0.84972  0.05200000
## 25:  75  1     Ustat 25000 0.95264 0.85736  9.52800000
## 26:  75  1 perm-perc 25000 0.95528 0.95536  0.00800000
## 27: 100  1     Ustat 25000 0.98780 0.96216  2.56400000
## 28: 150  1     Ustat 25000 0.99940 0.99760  0.18000000
## 29: 200  1     Ustat 25000 0.99996 0.99988  0.00800000
##       n mu    method     N powerNB powerWR    mismatch

## *** coverage
allResS.tempo[n==200 & method == "perm-perc" & mu == 0]
##       method  statistic mu sigma   n   rep truth        bias   power coverage
## 1: perm-perc netBenefit  0     2 200 25000     0 0.000088040 0.05836      NaN
## 2: perm-perc   winRatio  0     2 200 25000     1 0.007422373 0.05840      NaN
allResS.tempo[n==10 & method != "perm-perc" & mu == 0 & statistic == "netBenefit"]
##         method  statistic mu sigma  n   rep truth       bias   power coverage
## 1:       Ustat netBenefit  0     2 10 25000     0 -0.0005704 0.08944  0.91056
## 2: Ustat-trans netBenefit  0     2 10 25000     0 -0.0005704 0.04140  0.95860
## 3:   perm-stud netBenefit  0     2 10 25000     0 -0.0005704 0.05128      NaN
## 4:   boot-perc netBenefit  0     2 10 25000     0 -0.0005704 0.06536  0.94084
## 5:   boot-stud netBenefit  0     2 10 25000     0 -0.0005704 0.03960  0.96052

gg.cov <- ggplot(allResS.tempo[allResS.tempo$method.legend %in% convertion[1:4]],
                 aes(x = n, y = coverage,
                     group = method.legend, shape = method.legend, color = method.legend, linetype = method.legend))
gg.cov <- gg.cov + geom_hline(yintercept = 0.95, linewidth = 2, color = "lightgray")
gg.cov <- gg.cov + geom_line(linewidth = 1.15) + geom_point(size = 3)
gg.cov <- gg.cov + facet_grid(mu.legend ~ statistic.legend) 
gg.cov <- gg.cov + scale_linetype_manual(values = c(1,1,2,2), breaks = convertion[1:4]) 
gg.cov <- gg.cov + scale_shape_manual(values = c(8,15,8,15), breaks = convertion[1:4]) 
gg.cov <- gg.cov + scale_color_manual(values = rep(c("darkgray","black"),2), breaks = convertion[1:4]) 
gg.cov <- gg.cov + coord_cartesian(y = c(0.9,0.96))
gg.cov <- gg.cov + labs(x = "Sample size (per group)", y = "Coverage", color = "Inferential method\n", shape = "Inferential method\n", linetype = "Inferential method\n")
gg.cov <- gg.cov + guides(shape = guide_legend(nrow = 2, byrow = TRUE))
gg.cov <- gg.cov + guides(color = guide_legend(nrow = 2, byrow = TRUE))
gg.cov <- gg.cov + guides(linetype = guide_legend(nrow = 2, byrow = TRUE))
gg.cov <- gg.cov + theme(text = element_text(size=15), 
                         axis.line = element_line(linewidth = 1.25),
                         axis.ticks = element_line(linewidth = 2),
                         axis.ticks.length=unit(.25, "cm"),
                         legend.key.width = unit(4,"line"),
                         legend.position = "bottom")
gg.cov

ggsave(gg.cov, filename = "c:/Users/hpl802/Documents/GitHub/BuyseTest/inst/book/figures/fig_inference_simulation-coverage.pdf", width = 11, height = 8,
        device = cairo_pdf)

## *** type 1 error
gg.type1 <- ggplot(allResS.tempo[mu==0], aes(x = n, y = power,
                                             group = method.legend, shape = method.legend, color = method.legend, linetype = method.legend))
gg.type1 <- gg.type1 + geom_hline(yintercept = 0.05, size = 2, color = "lightgray")
gg.type1 <- gg.type1 + geom_line(linewidth = 1.15) + geom_point(size = 3)
gg.type1 <- gg.type1 + facet_grid(mu.legend ~ statistic.legend)
gg.type1 <- gg.type1 + scale_linetype_manual(values = c(1,1,2,2,3,3), breaks = convertion) 
gg.type1 <- gg.type1 + scale_shape_manual(values = rep(c(8,15),3), breaks = convertion) 
gg.type1 <- gg.type1 + scale_color_manual(values = rep(c("darkgray","black"),3), breaks = convertion) 
gg.type1 <- gg.type1 + coord_cartesian(y = c(0.04,0.1))
gg.type1 <- gg.type1 + guides(color = guide_legend(nrow = 6, byrow = TRUE))
gg.type1 <- gg.type1 + labs(x = "Sample size (per group)", y = "Type 1 error", color = "Inferential method\n", shape = "Inferential method\n", linetype = "Inferential method\n")
gg.type1 <- gg.type1 + theme(text = element_text(size=15), 
                             axis.line = element_line(linewidth = 1.25),
                             axis.ticks = element_line(linewidth = 2),
                             axis.ticks.length=unit(.25, "cm"),
                             legend.text=element_text(size=11),
                             legend.key.width = unit(4,"line"))
gg.type1

ggsave(gg.type1, filename = "c:/Users/hpl802/Documents/GitHub/BuyseTest/inst/book/figures/fig_inference_simulation-type1.pdf", width = 9, height = 8,
       device = cairo_pdf)


library(ggpubr)
gg.type1cov <- ggarrange(gg.type1 + guides(color = guide_legend(nrow = 3, byrow = TRUE), linetype = guide_legend(nrow = 3, byrow = TRUE), shape = guide_legend(nrow = 3, byrow = TRUE)) + ggtitle("A"),
                         gg.cov + ggtitle("B"), common.legend = TRUE, legend = "bottom", ncol = 1,
                 heights = c(3, 4, 4), nrow = 2, align = "v")
print(gg.type1cov)

ggsave(gg.type1cov, filename = "c:/Users/hpl802/Documents/GitHub/BuyseTest/inst/book/figures/fig_inference_simulation-type1coverage.pdf", width = 9, height = 9,
       device = cairo_pdf)

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

## ** Results
## *** type 1 error
gg.FWER <- ggplot(allResSL.tempo2,aes(x = n, y = value))
gg.FWER <- gg.FWER + geom_hline(data = data.frame(x =  unique(allResSL.tempo2$n), y = 0.05, mu.legend = paste0("\u0394\u03bc=",0)), aes(yintercept = y) ,color = "gray", linewidth = 1.5)
gg.FWER <- gg.FWER + geom_line(aes(group = variable.legend, linetype = variable.legend), linewidth = 1.15) + geom_point(aes(group = variable, shape = variable.legend), size = 3)
gg.FWER <- gg.FWER + facet_wrap(~mu.legend)
gg.FWER <- gg.FWER + labs(x = "Sample size (per group)", y = "Rejection rate", linetype = "Multiplicity adjusment", shape = "Multiplicity adjusment")
gg.FWER <- gg.FWER + theme(text = element_text(size=15), 
                           axis.line = element_line(linewidth = 1.25),
                           axis.ticks = element_line(linewidth = 2),
                           axis.ticks.length=unit(.25, "cm"),
                           legend.text=element_text(size=11),
                           legend.key.width = unit(4,"line"),
                           legend.position = "bottom")
gg.FWER
ggsave(gg.FWER, filename = "c:/Users/hpl802/Documents/GitHub/BuyseTest/inst/book/figures/fig_inference_simulation-FWER.pdf", width = 9, height = 6,
       device = cairo_pdf)
##----------------------------------------------------------------------
### BUILD.R ends here
