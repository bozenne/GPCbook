### figureSoftware-6.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt  9 2023 (17:22) 
## Version: 
## Last-Updated: okt  9 2023 (18:09) 
##           By: Brice Ozenne
##     Update #: 10
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

library(BuyseTest)
library(ggplot2)

## * generate data
argsSurv <- list(name = c("OS","PFS"),
                 name.censoring = c("statusOS","statusPFS"),
                 scale.C = c(8.995655, 4.265128),
                 scale.T = c(13.76543, 7.884477),
                 shape.C = c(1.28993, 1.391015),
                 shape.T = c(1.275269, 1.327461),
                 scale.censoring.C = c(34.30562, 20.748712),
                 scale.censoring.T = c(27.88519, 17.484281),
                 shape.censoring.C = c(1.369449, 1.463876),
                 shape.censoring.T = c(1.490881, 1.835526))

argsTox <- list(name = "toxicity",
                p.C =  c(1.17, 2.92, 36.26, 39.18, 19.88, 0.59)/100,
                p.T = c(3.51, 4.09, 23.39, 47.37, 21.05, 0.59)/100,
                rho.T = 1, rho.C = 1)

set.seed(1)
dt.data <- simBuyseTest(n.T = 200, n.C = 200,
                        argsBin = argsTox,
                        argsCont = NULL,
                        argsTTE = argsSurv,
                        level.strata = c("M","F"), names.strata = "gender")

## * GPC
dt.data$toxicity.num <- as.numeric(dt.data$toxicity)

eRBB.BT <- BuyseTest(treatment ~ cont(toxicity.num, operator = "<0") + tte(OS, statusOS),
                     data=dt.data, trace = FALSE)
eRBB.Se <- sensitivity(eRBB.BT, threshold = list(1:5,c(0,5,10)),
                       band = TRUE, adj.p.value = TRUE, seed = 10, trace = FALSE)

## * generate figure 7
theme_set(theme_bw())
figure6 <- autoplot(eRBB.Se) + facet_wrap(~OS, labeller = label_both)
figure6 <- figure6 + scale_color_grey(start = 0.1, end = .6)
figure6 <- figure6 + scale_fill_grey(start = 0.1, end = .6)
figure6 <- figure6 + ylab("Net Treatment Benefit")
figure6 <- figure6 + theme(text = element_text(size=20), 
                           axis.line = element_line(linewidth = 1.25),
                           axis.ticks = element_line(linewidth = 1.25),
                           axis.ticks.length=unit(.25, "cm"),
                           legend.key.size = unit(2,"line"))

pdf("figures/fig_software_sensitivity.pdf", width = 12, height = 8)
figure6
dev.off()

##----------------------------------------------------------------------
### figureSoftware-6.R ends here
