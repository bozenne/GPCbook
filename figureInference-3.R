### figureInference-3.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Oct  9 2023 (11:38) 
## Version: 
## Last-Updated: Oct  9 2023 (11:51) 
##           By: Brice Ozenne
##     Update #: 7
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

## * generate figure 3

## ** run GPC
GPC.figure3 <- BuyseTest(treatment ~ cont(score), data = dtInference, trace = FALSE,
                         seed = 10, method.inference  = "studentized permutation", n.resampling = 1e4)

## ** gather results
NB.perm <- GPC.figure3@DeltaResampling[,,"netBenefit"]
seNB.perm <- GPC.figure3@covarianceResampling[,,"netBenefit"]

dt.perm <- rbind(data.table(estimate = NB.perm,
                            type = "permutation estimates"),
                 data.table(estimate = NB.perm/sqrt(seNB.perm),
                            type = "permutation statistics")
                 )

dt.permQ <- dt.perm[, .(Qlower = quantile(estimate, prob = 0.025, na.rm = TRUE),
                        Qupper = quantile(estimate, prob = 0.975, na.rm = TRUE)),
                    by = "type"]
dt.permQ
dt.obsQ <- rbind(data.table(estimate = coef(GPC.figure3),
                            type = "permutation estimates"),
                 data.table(estimate = coef(GPC.figure3)/confint(GPC.figure3)$se,
                            type = "permutation statistics")
                 )

dt.perm[abs(estimate)>10,estimate := NA]
dt.perm$type <- factor(dt.perm$type, levels = unique(dt.perm$type))
dt.permQ$type <- factor(dt.permQ$type, levels = unique(dt.perm$type))


## ** generate figure

gg.histPerm <- ggplot()
gg.histPerm <- gg.histPerm + geom_histogram(data = dt.perm, mapping = aes(x=estimate, y=after_stat(2 * count / sum(count))), color = "black")
gg.histPerm <- gg.histPerm + geom_vline(data = dt.obsQ, mapping = aes(xintercept=estimate), color = "red", linetype = 1, size = 1.25)
gg.histPerm <- gg.histPerm + geom_vline(data = dt.permQ, mapping = aes(xintercept=Qlower), color = "gray", linetype = 2, size = 1.25)
gg.histPerm <- gg.histPerm + geom_vline(data = dt.permQ, mapping = aes(xintercept=Qupper), color = "gray", linetype = 2, size = 1.25)

gg.histPerm <- gg.histPerm + facet_grid(~type, scales="free")
gg.histPerm <- gg.histPerm + scale_y_continuous(labels = scales::percent)
gg.histPerm <- gg.histPerm + labs(y = "Relative frequency", x = NULL)
gg.histPerm <- gg.histPerm + theme(text = element_text(size=15), 
                                   axis.line = element_line(linewidth = 1.25),
                                   axis.ticks = element_line(linewidth = 2),
                                   axis.ticks.length=unit(.25, "cm"),
                                   legend.key.size = unit(3,"line"))

ggsave(gg.histPerm, filename = "figures/fig_inference_permutation.pdf", width = 9, height = 6)

##----------------------------------------------------------------------
### figureInference-3.R ends here
