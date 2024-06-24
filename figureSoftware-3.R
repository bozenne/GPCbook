### figureSoftware-3.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt  9 2023 (15:51) 
## Version: 
## Last-Updated: okt  9 2023 (18:03) 
##           By: Brice Ozenne
##     Update #: 4
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

library(ggpubr)
library(ggplot2)
library(prodlim)
library(BuyseTest)

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

## * generate figure 4
dtPC.toxW <- prop.table(table(dt.data$treatment,
                              dt.data$toxicity))

dtPC.toxL <- as.data.frame(dtPC.toxW, responseName = "Probability")
names(dtPC.toxL)[1:2] <- c("treatment","grade")

## ** Toxicity
theme_set(theme_bw())
colorG2R <- scales::seq_gradient_pal(low = rgb(0.9,0.9,0.9),
                                     high = rgb(0.1,0.1,0.1))

figure3.A <- ggplot(dtPC.toxL, aes(x = treatment, fill = grade, y = Probability))
figure3.A <- figure3.A + geom_bar(position = position_fill(reverse = TRUE),
                            stat = "identity")
figure3.A <- figure3.A + scale_y_continuous(labels = scales::percent)
figure3.A <- figure3.A + scale_fill_manual("Worse\nadverse event",
                                     values = colorG2R(seq(0,1,length.out=6)))
figure3.A 


## ** assemble
pdf("figures/fig_software_hist-tox.pdf", width = 5, height = 5)
figure3.A + theme(text = element_text(size=15), 
                       axis.line = element_line(linewidth = 1.25),
                       axis.ticks = element_line(linewidth = 1.25),
                       axis.ticks.length=unit(.25, "cm"),
                       legend.key.size = unit(2,"line"))
dev.off()
pdf("figures/fig_software_KM-OS.pdf", width = 5, height = 5)
par(mar = rep(1,4))
plot(prodlim(Hist(OS,statusOS) ~ treatment, data = dt.data), col = c(rgb(0.1,0.1,0.1),rgb(0.7,0.7,0.7)),
     lty = c(4,1))
dev.off()
##----------------------------------------------------------------------
### figureSoftware-3.R ends here
