### figureInference-5.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Oct  9 2023 (11:38) 
## Version: 
## Last-Updated: sep 24 2024 (17:44) 
##           By: Brice Ozenne
##     Update #: 14
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

## * generate figure 5

## ** run GPC
GPC.figure3 <- BuyseTest(treatment ~ cont(score), data = dtInference, trace = FALSE,
                         seed = 10, method.inference  = "studentized permutation", n.resampling = 1e4)

## ** gather results
NB.perm <- GPC.figure3@DeltaResampling[,,"netBenefit"]
seNB.perm <- GPC.figure3@covarianceResampling[,,"netBenefit"]

dt.perm <- rbind(data.table(estimate = NB.perm,
                            type = "Permuted NTB estimates"),
                 data.table(estimate = NB.perm/sqrt(seNB.perm),
                            type = "Studentized permuted NTB")
                 )[abs(estimate)<5]

dt.permQ <- dt.perm[, .(Qlower = quantile(estimate, prob = 0.025),
                        Qupper = quantile(estimate, prob = 0.975)),
                    by = "type"]

dt.obsQ <- rbind(data.table(estimate = coef(GPC.figure3),
                            type = "Permuted NTB estimates"),
                 data.table(estimate = coef(GPC.figure3)/confint(GPC.figure3)$se,
                            type = "Studentized permuted NTB")
                 )

dt.perm[,density := dnorm(estimate, mean = mean(estimate, na.rm = TRUE), sd = sd(estimate, na.rm = TRUE)), by = "type"]
dt.perm$type <- factor(dt.perm$type, levels = unique(dt.perm$type))
dt.permQ$type <- factor(dt.permQ$type, levels = unique(dt.perm$type))

## ** generate figure
theme_set(theme_bw())
gg.histPerm <- ggplot()
gg.histPerm <- gg.histPerm + geom_histogram(data = dt.perm, mapping = aes(x = estimate, y = after_stat(density)), color = "black")
gg.histPerm <- gg.histPerm + geom_line(data = dt.perm, mapping = aes(x = estimate, y = density), color = "darkgrey", linewidth = 1.25)
gg.histPerm <- gg.histPerm + geom_point(data = dt.permQ, mapping = aes(x=Qlower, y = -0.025, color = "quantile [2.5%;97.5%]", shape = "quantile [2.5%;97.5%]"), size = 5)
gg.histPerm <- gg.histPerm + geom_point(data = dt.permQ, mapping = aes(x=Qupper, y = -0.025, color = "quantile [2.5%;97.5%]", shape = "quantile [2.5%;97.5%]"), size = 5)
gg.histPerm <- gg.histPerm + geom_point(data = dt.obsQ, mapping = aes(x=estimate, y = -0.075, color = "observed estimate", shape = "observed estimate"), size = 5)
gg.histPerm <- gg.histPerm + facet_grid(~type, scales="free")
gg.histPerm <- gg.histPerm + coord_cartesian(ylim = c(-0.04,1.6))
gg.histPerm <- gg.histPerm + labs(y = "Relative frequency", x = NULL, shape = "", color = "")
gg.histPerm <- gg.histPerm + scale_colour_manual(values = c("black","grey"))
gg.histPerm <- gg.histPerm + theme(text = element_text(size=15),
                                   panel.grid.minor.x = element_blank(),
                                   panel.grid.major.x = element_blank(),
                                   axis.line = element_line(linewidth = 1.25),
                                   axis.ticks = element_line(linewidth = 2),
                                   axis.ticks.length=unit(.25, "cm"),
                                   legend.key.size = unit(3,"line"),
                                   legend.position = "bottom")
gg.histPerm

ggsave(gg.histPerm, filename = "figures/fig_inference_permutation.pdf", width = 9, height = 6.5)

##----------------------------------------------------------------------
### figureInference-5.R ends here
