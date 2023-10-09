### mainSoftware.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Oct  9 2023 (10:13) 
## Version: 
## Last-Updated: okt  9 2023 (18:28) 
##           By: Brice Ozenne
##     Update #: 42
##----------------------------------------------------------------------
## 
### Commentary: 
## in sillico --> synthetic
## simplistic simulation since no correlation between endpoint --> good point
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

library(BuyseTest)

## * 3.2 Setting the stage

## ** 3.2.2 Generating synthetic data
set.seed(10) ## initialize the pseudo-random number generator 
dt0.data <- simBuyseTest(100)
dt0.data

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
dt.data

set.seed(10)
rbind(simBuyseTest(n.T = 100, n.C = 100,
                   argsBin = NULL,
                   argsCont = list(mu.C = 1, mu.T = 2),
                   argsTTE = NULL,
                   prefix.cluster = "M", level.strata = "M", names.strata = "gender"),
      simBuyseTest(n.T = 100, n.C = 100,
                   argsBin = NULL,
                   argsCont = list(mu.C = 10, mu.T = 20),
                   argsTTE = NULL,
                   prefix.cluster = "F", level.strata = "F", names.strata = "gender")
      )

dtPC.toxW <- prop.table(table(dt.data$treatment,
                              dt.data$toxicity))
dtPC.toxW * 100

library(ggplot2)
ggplot(dt.data, aes(x = toxicity, y = OS, fill = treatment)) + geom_boxplot()
ggplot(dt.data, aes(x = toxicity, y = PFS, fill = treatment)) + geom_boxplot()


## ** 3.2.3 Alternative software solutions
library(WR)

WRrec(ID = hfaction_cpx9$patid,
      time = hfaction_cpx9$time,
      status = hfaction_cpx9$status,
      trt = hfaction_cpx9$trt_ab,
      strata = hfaction_cpx9$age60,
      naive = TRUE)

## WRrec(ID = hfaction_cpx9$patid, time = hfaction_cpx9$time, status = hfaction_cpx9$status, 
##     trt = hfaction_cpx9$trt_ab, strata = hfaction_cpx9$age60, 
##     naive = TRUE)

##             N Rec. Event Death Med. Follow-up
## Control   221        571    57       28.62295
## Treatment 205        451    36       27.57377

## Analysis of last-event-assisted WR (LWR; recommended), first-event-assisted WR (FWR), and naive WR (NWR):
##     Win prob Loss prob WR (95% CI)*      p-value
## LWR 50.4%    38.2%     1.32 (1.05, 1.66) 0.0189 
## FWR 50.4%    38.3%     1.32 (1.04, 1.66) 0.0202 
## NWR 47%      35%       1.34 (1.05, 1.72) 0.0193 
## -----
## *Note: The scale of WR should be interpreted with caution as it depends on 
## censoring distribution without modeling assumptions.

## * 3.3 GPC with a single endpoint
dtPC.toxL <- as.data.frame(dtPC.toxW, responseName = "Probability")
names(dtPC.toxL)[1:2] <- c("treatment","grade")

## ** 3.3.2 GPC analysis for completely-observed data
dt.data$toxicity.num <- as.numeric(dt.data$toxicity)
wilcox.test(toxicity.num ~ treatment, data = dt.data)

eTox.BT <- BuyseTest(treatment ~ cont(toxicity.num, operator = "<0"),
                     data=dt.data)
summary(eTox.BT)
print(eTox.BT, percentage = FALSE)

## ** 3.3.3 Alternative GPC-based effects
confint(eTox.BT, statistic = "favorable")

eTox.BThalf <- BuyseTest(treatment ~ cont(toxicity.num, operator = "<0"),
                         add.halfNeutral = TRUE,
                         data = dt.data,
                         trace = FALSE)
confint(eTox.BThalf, statistic = "favorable")
confint(eTox.BThalf)

## ** 3.3.4 Additional options and comparison of the results with alternatitve software

BuyseTest.options(trace = 0)

library(asht)
dt.data$treatment2 <- relevel(dt.data$treatment,"T")
wmwTest(toxicity.num ~ treatment2, data = dt.data)

library(pim)
e.pim <- pim(toxicity.num ~ treatment2, data = dt.data)
summary(e.pim)

BuyseTest(treatment ~ cont(toxicity.num, operator = "<0"),
          add.halfNeutral = TRUE, method.inference = "permutation",
          data=dt.data, cpus = 5, n.resampling = 1e4, seed = 10)

## ** 3.3.5 Threshold of clinical relevance
eTox.BT2 <- BuyseTest(treatment ~ cont(toxicity.num, threshold = 2, operator = "<0"),
                     data=dt.data, keep.pairScore = TRUE, trace = FALSE)
print(eTox.BT2)

getPairScore(eTox.BT2)

dt.data[c(3:4,201),c("id","treatment","OS","statusOS","toxicity","gender")]

model.tables(eTox.BT, columns = "threshold")

## ** 3.3.6 Stratification
ffG <- treatment ~ cont(toxicity.num, operator = "<0") + strata(gender)
eTox.BTG <- BuyseTest(ffG, data=dt.data, keep.pairScore = TRUE, trace = FALSE)
summary(eTox.BTG)
summary(eTox.BTG, percentage = FALSE)

getPairScore(eTox.BTG)

confint(eTox.BTG, strata = TRUE)

e.pimS <- pim(toxicity.num ~ treatment + gender, data = dt.data,
              link = "identity")
summary(e.pimS)

## ** 3.3.7 Handling right-censoring when assessing efficacy

dt.data[,.(censoring=mean(statusOS==0)),by = "treatment"]

eEff.BT <- BuyseTest(treatment ~ tte(OS, statusOS), data=dt.data,
                     keep.pairScore = TRUE, trace = FALSE)

getPairScore(eEff.BT)[c(1,2,2623,8553),]

print(eEff.BT)

eEff.BT2 <- BuyseTest(treatment ~ tte(OS, statusOS), data=dt.data,
                      scoring.rule = "Gehan", keep.pairScore = TRUE, trace = FALSE)
print(eEff.BT2)

dt30.data <- data.table::copy(dt.data)
dt30.data[OS>30, c("OS", "statusOS") := .(30,0)]

eEff.BT30 <- BuyseTest(treatment ~ tte(OS, statusOS, restriction = 25), data=dt30.data,
                       keep.pairScore = TRUE, trace = FALSE)
print(eEff.BT30)

dt.data[c(44,211)]
getPairScore(eEff.BT30)[index.C==44 & index.T == 211,]
getPairScore(eEff.BT)[index.C==44 & index.T == 211,]

## * 3.4 Benefit-risk analysis using GPC

## ** 3.4.1 Prioritized & non-prioritized analyses
eBRB.BT <- BuyseTest(treatment ~ tte(OS, statusOS) + cont(toxicity.num, operator = "<0"),
                     data=dt.data, trace = FALSE)
print(eBRB.BT)

eRBB.BT <- BuyseTest(treatment ~ cont(toxicity.num, operator = "<0") + tte(OS, statusOS),
                     data=dt.data, trace = FALSE)
print(eRBB.BT)

eNH.BT <- BuyseTest(treatment ~ cont(toxicity.num, operator = "<0") + tte(OS, statusOS),
                    data=dt.data, hierarchical = FALSE, trace = FALSE)
print(eNH.BT)

library(ggplot2)
eRBB.plot <- plot(eRBB.BT)
eNH.plot <- plot(eNH.BT)
ggpubr::ggarrange(eRBB.plot$plot + ggtitle("Hierarchical"),
                  eNH.plot$plot + ggtitle("Non-hierarchical"),
                  common.legend = TRUE, legend = "bottom")

rbind("prioritized" = confint(eRBB.BT, transform = FALSE, endpoint = 1),
      "non-prioritized" = confint(eNH.BT, transform = FALSE, endpoint = 1))

## ** 3.4.2 Thresholds of clinical relevance

eSH.BT <- BuyseTest(treatment ~ tte(OS, statusOS, threshold = 28)
                              + cont(toxicity.num, operator = "<0", threshold = 2)
                              + tte(OS, statusOS, threshold = 14)
                              + cont(toxicity.num, operator = "<0"),
                    data=dt.data, trace = FALSE)
print(eSH.BT)

## ** 3.4.3 Encoding of the outcome

dt.data$OS2 <- dt.data$OS
dt.data$OS2[dt.data$statusOS==0] <- 150

print(BuyseTest(treatment ~ tte(OS2, statusOS), data=dt.data, trace = FALSE))

eD2.BT <- BuyseTest(treatment ~ bin(statusOS, operator = "<0") + tte(OS2, statusOS), data=dt.data, trace = FALSE)
print(eD2.BT)

dt.data$toxicity2 <- dt.data$toxicity.num
dt.data$toxicity2[dt.data$statusOS==1] <- -1

eBRB2.BT <- BuyseTest(treatment ~ bin(statusOS, operator = "<0") + cont(toxicity2, operator = "<0"), data=dt.data, trace = FALSE)
print(eBRB2.BT)

eR2.BT <- BuyseTest(treatment ~ cont(toxicity2, operator = "<0"),
                    data=dt.data[statusOS==0], trace = FALSE)
print(eR2.BT, percentage = FALSE)
(1947 - 2279)/40000

## ** 3.4.4 Sensitivity analyses

eRBB.Se <- sensitivity(eRBB.BT, threshold = list(1:5,c(0,5,10)),
                       band = TRUE, adj.p.value = TRUE, seed = 10, trace = FALSE)
eRBB.Se[c(1,2,6),]

autoplot(eRBB.Se) + facet_wrap(~OS, labeller = label_both)

library(lava)
eRBB.Hdecomp <- iid(eRBB.Se)
dim(eRBB.Hdecomp)

eRBB.cor <- cor(eRBB.Hdecomp)
range(eRBB.cor[lower.tri(eRBB.cor)])

range(eRBB.Se$adj.p.value/eRBB.Se$p.value)

e.MBT <- BuyseMultComp(list("OS-tox" = eBRB.BT, "tox-OS" = eRBB.BT, "threshold" = eSH.BT),
                       cluster = "id",
                       seed = 10)
e.MBT

## * 3.5 Power calculation for GPC

## ** 3.5.1 Data-generating mechanism
simFCT <- function(n.C, n.T){
     out <- rbind(data.frame(Y=stats::rt(n.C, df = 5), group=0),
                  data.frame(Y=stats::rt(n.T, df = 5) + 1, group=1))
     return(out)
}
set.seed(10)
simFCT(2,2)

simFCT2 <- function(n.T, n.C){
  out <- simBuyseTest(n.T, n.C,
                      argsBin = argsTox,
                      argsCont = NULL,
                      argsTTE = argsSurv,
                      level.strata = c("M","F"), names.strata = "gender")
  out$toxicity <- as.numeric(out$toxicity)
  return(out)
}
set.seed(10)
simFCT2(2,2)

## ** 3.5.2 Simulation-based power and sample size estimation
e.power <- powerBuyseTest(formula = treatment ~ tte(OS, statusOS, threshold = 5) + cont(toxicity, operator = "<0"),
                          sim = simFCT2, sample.size = c(10,50,100),
                          n.rep = 100, seed = 10)
summary(e.power)

e.power2 <- powerBuyseTest(formula = treatment ~ tte(OS, statusOS, threshold = 5) + cont(toxicity, operator = "<0"),
                           sim = simFCT2, sample.size = c(10,50,100),
                           n.rep = 1000, seed = 10, cpus = 5, export.cpus = c("argsTox", "argsSurv"))
print(e.power2, endpoint = "all")

e.nSearch <- powerBuyseTest(formula = treatment ~ tte(OS, statusOS, threshold = 5)
                            + cont(toxicity, operator = "<0"),
                            sim = simFCT2, power = 0.8, max.sample.size = 1000,
                            n.rep = c(1000,10), seed = 10, trace = 2, 
                            cpus = 5, export.cpus = c("argsTox", "argsSurv"))
print(e.nSearch)

##----------------------------------------------------------------------
### mainSoftware.R ends here
