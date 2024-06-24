### figureSoftware-5.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt  9 2023 (17:22) 
## Version: 
## Last-Updated: okt  9 2023 (17:40) 
##           By: Brice Ozenne
##     Update #: 6
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

eSH.BT <- BuyseTest(treatment ~ tte(OS, statusOS, threshold = 28)
                    + cont(toxicity.num, operator = "<0", threshold = 2)
                    + tte(OS, statusOS, threshold = 14)
                    + cont(toxicity.num, operator = "<0"),
                    data=dt.data, trace = FALSE)

eBRB.BT <- BuyseTest(treatment ~ tte(OS, statusOS) + cont(toxicity.num, operator = "<0"),
                     data=dt.data, trace = FALSE)

## * generate figure 5
theme_set(theme_bw())
colorG2R <- scales::seq_gradient_pal(low = rgb(0.9,0.9,0.9),
                                     high = rgb(0.1,0.1,0.1))


label5.A <- c("OS\n(any difference)","Toxicity\n(any difference)")
figure5.A <- plot(eBRB.BT, label.endpoint = label5.A)$plot + ggtitle("No threshold")
figure5.A <- figure5.A + scale_fill_grey(start = 0.2, end = .9)
figure5.A <- figure5.A + theme(text = element_text(size=20), 
                               axis.line = element_line(linewidth = 1.25),
                               axis.ticks = element_line(linewidth = 1.25),
                               axis.ticks.length=unit(.25, "cm"),
                               legend.key.size = unit(2,"line"))

label5.B <- c("OS\n(\U2265 28 days)","Toxicity\n(\U2265 2 grade)","OS\n(\U2265 14 days)","Toxicity\n(any difference)")
figure5.B <- plot(eSH.BT, label.endpoint = label5.B)$plot + ggtitle("With threshold")
figure5.B <- figure5.B + scale_fill_grey(start = 0.2, end = .9)
figure5.B <- figure5.B + theme(text = element_text(size=20), 
                               axis.line = element_line(linewidth = 1.25),
                               axis.ticks = element_line(linewidth = 1.25),
                               axis.ticks.length=unit(.25, "cm"),
                               legend.key.size = unit(2,"line"))


figure5 <- ggpubr::ggarrange(figure5.A, figure5.B,
                             common.legend = TRUE, legend = "bottom", widths = c(1,1.5))

pdf("figures/fig_software_hierarchical-threshold.pdf", width = 12, height = 8)
figure5
dev.off()

##----------------------------------------------------------------------
### figureSoftware-5.R ends here
