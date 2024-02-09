### figureInference-4.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Oct  9 2023 (11:38) 
## Version: 
## Last-Updated: feb  9 2024 (15:11) 
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
library(data.table)
library(ggplot2)

## * simulate data
set.seed(10)
n.data <- 10
dtInference <- simBuyseTest(n.data)
dtInference[, score := round(score,1)]

## * generate figure 2

## ** run GPC
GPC.figure2 <- BuyseTest(treatment ~ cont(score), data = dtInference, trace = FALSE,
                         seed = 10, method.inference  = "studentized bootstrap", strata.resampling = "treatment", n.resampling = 1e4)

## ** gather results
NB.boot <- GPC.figure2@DeltaResampling[,,"netBenefit"]
seNB.boot <- GPC.figure2@covarianceResampling[,,"netBenefit"]

dt.boot <- rbind(data.table(estimate = NB.boot,
                            scale = "original scale", type = "bootstrap estimates"),
                 data.table(estimate = atanh(NB.boot),
                            scale = "atanh scale", type = "bootstrap estimates"),
                 data.table(estimate = (NB.boot-coef(GPC.figure2))/sqrt(seNB.boot),
                            scale = "original scale", type = "centered bootstrap statistics"),
                 data.table(estimate = (atanh(NB.boot)-atanh(coef(GPC.figure2)))/sqrt(seNB.boot/(1-NB.boot^2)^2),
                            scale = "atanh scale", type = "centered bootstrap statistics"))
dt.boot <- dt.boot[!is.na(estimate) & !is.infinite(estimate) & abs(estimate)<5]

dt.bootQ <- dt.boot[, .(Qlower = quantile(estimate, prob = 0.025),
                        Qupper = quantile(estimate, prob = 0.975)),
                    by = c("scale","type")]
dt.bootQ

dt.boot$type <- factor(dt.boot$type, levels = unique(dt.boot$type))
dt.boot$scale <- factor(dt.boot$scale, levels = unique(dt.boot$scale))
dt.bootQ$type <- factor(dt.bootQ$type, levels = unique(dt.boot$type))
dt.bootQ$scale <- factor(dt.bootQ$scale, levels = unique(dt.boot$scale))
dt.bootLine <- dt.boot[, .(estimate = seq(-5,5,length.out=1000),
                           density = dnorm(seq(-5,5,length.out=1000), mean = mean(estimate), sd = sd(estimate))), by = "type"][density > 0.001]


## ** generate figure
gg.histBoot <- ggplot()
gg.histBoot <- gg.histBoot + geom_histogram(data = dt.boot, mapping = aes(x=estimate, y=after_stat(density)), color = "black")
gg.histBoot <- gg.histBoot + geom_line(data = dt.bootLine, mapping = aes(x = estimate, y = density), color = "darkgrey", linewidth = 1.25)
gg.histBoot <- gg.histBoot + geom_point(data = dt.bootQ, mapping = aes(x=Qlower, y = -0.045, color = "quantile [2.5%;97.5%]", shape = "quantile [2.5%;97.5%]"), size = 4)
gg.histBoot <- gg.histBoot + geom_point(data = dt.bootQ, mapping = aes(x=Qupper, y = -0.045, color = "quantile [2.5%;97.5%]", shape = "quantile [2.5%;97.5%]"), size = 4)
gg.histBoot <- gg.histBoot + facet_grid(scale~type, scales="free")
gg.histBoot <- gg.histBoot + labs(y = "Relative frequency", x = NULL, shape = "", color = "")
gg.histBoot <- gg.histBoot + scale_colour_manual(values = c("grey")) + scale_shape_manual(values = 17)
gg.histBoot <- gg.histBoot + scale_x_continuous(n.breaks = 8)
gg.histBoot <- gg.histBoot + theme(text = element_text(size=15), 
                                   panel.grid.minor.x = element_blank(),
                                   panel.grid.major.x = element_blank(),
                                   axis.line = element_line(linewidth = 1.25),
                                   axis.ticks = element_line(linewidth = 1.5),
                                   axis.ticks.length=unit(.25, "cm"),
                                   legend.key.size = unit(3,"line"),
                                   legend.position = "bottom")
gg.histBoot

## gg.histBoot
ggsave(gg.histBoot, filename = "figures/fig_inference_bootstrap.pdf", width = 9, height = 6)

##----------------------------------------------------------------------
### figureInference-4.R ends here
