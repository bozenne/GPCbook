### figureInference-2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Oct  9 2023 (11:38) 
## Version: 
## Last-Updated: Oct  9 2023 (11:49) 
##           By: Brice Ozenne
##     Update #: 5
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
                            scale = "original scale", type = "bootstrap centered statistics"),
                 data.table(estimate = (atanh(NB.boot)-atanh(coef(GPC.figure2)))/sqrt(seNB.boot/(1-NB.boot^2)^2),
                            scale = "atanh scale", type = "bootstrap centered statistics"))

dt.bootQ <- dt.boot[, .(Qlower = quantile(estimate, prob = 0.025, na.rm = TRUE),
                        Qupper = quantile(estimate, prob = 0.975, na.rm = TRUE)),
                    by = c("scale","type")]
dt.bootQ

dt.boot$estimate[is.infinite(dt.boot$estimate)] <- NA
dt.boot$estimate[abs(dt.boot$estimate)>15] <- NA
dt.boot$type <- factor(dt.boot$type, levels = unique(dt.boot$type))
dt.boot$scale <- factor(dt.boot$scale, levels = unique(dt.boot$scale))
dt.bootQ$type <- factor(dt.bootQ$type, levels = unique(dt.boot$type))
dt.bootQ$scale <- factor(dt.bootQ$scale, levels = unique(dt.boot$scale))

## ** generate figure
gg.histBoot <- ggplot()
gg.histBoot <- gg.histBoot + geom_histogram(data = dt.boot, mapping = aes(x=estimate, y=after_stat(4 * count / sum(count))), color = "black")
gg.histBoot <- gg.histBoot + geom_vline(data = dt.bootQ, mapping = aes(xintercept=Qlower), color = "gray", linetype = 2, linewidth = 1.25)
gg.histBoot <- gg.histBoot + geom_vline(data = dt.bootQ, mapping = aes(xintercept=Qupper), color = "gray", linetype = 2, linewidth = 1.25)
gg.histBoot <- gg.histBoot + facet_grid(scale~type, scales="free")
gg.histBoot <- gg.histBoot + scale_y_continuous(labels = scales::percent)
gg.histBoot <- gg.histBoot + labs(y = "Relative frequency", x = NULL)
gg.histBoot <- gg.histBoot + theme(text = element_text(size=15), 
                                   axis.line = element_line(linewidth = 1.25),
                                   axis.ticks = element_line(linewidth = 2),
                                   axis.ticks.length=unit(.25, "cm"),
                                   legend.key.size = unit(3,"line"))
gg.histBoot

## gg.histBoot
ggsave(gg.histBoot, filename = "figures/fig_inference_bootstrap.pdf", width = 9, height = 6)

##----------------------------------------------------------------------
### figureInference-2.R ends here
