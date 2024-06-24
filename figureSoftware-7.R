### figureSoftware-8.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt  9 2023 (17:22) 
## Version: 
## Last-Updated: jun 13 2024 (11:05) 
##           By: Brice Ozenne
##     Update #: 12
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
library(lava)

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

eRBB.Hdecomp <- iid(eRBB.Se)
eRBB.cor <- cor(eRBB.Hdecomp)

## * generate figure 8

rownames(eRBB.cor) <- paste0("tox=",eRBB.Se$toxicity.num,";OS=",eRBB.Se$OS,"")
colnames(eRBB.cor) <- paste0("tox=",eRBB.Se$toxicity.num,";OS=",eRBB.Se$OS,"")

colTab <- scales::seq_gradient_pal(low = rgb(0.9,0.9,0.9),
                                   high = rgb(0.1,0.1,0.1))(seq(0,1,length.out = 30))

pdf("figures/fig_software_corIID.pdf", width = 8, height = 8)
par(mar  = c(6,6,2,2))
fields::image.plot(eRBB.cor[rev(colnames(eRBB.cor)),colnames(eRBB.cor)], axes = FALSE, col = colTab)
axis(1, at=seq(0,1,length.out=15), labels=rev(colnames(eRBB.cor)), las = 2)
axis(2, at=seq(0,1,length.out=15), labels=colnames(eRBB.cor), las = 2)
dev.off()

##----------------------------------------------------------------------
### figureSoftware-8.R ends here
