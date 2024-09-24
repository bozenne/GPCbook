### mainSoftware.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Oct  9 2023 (10:13) 
## Version: 
## Last-Updated: sep 24 2024 (17:49) 
##           By: Brice Ozenne
##     Update #: 124
##----------------------------------------------------------------------
## 
### Commentary: 
## Sections with the flag [EXTRA] contains R code that is not mentionned nor explained in the book chapter.
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

library(BuyseTest)
library(asht)
library(pim)
library(WR)
library(data.table)
options(datatable.print.class = FALSE)

## * 16.2 Setting the stage

## ** 16.1 [Extra] Introduction
data(FEVData, package = "pim")

e.pim <- pim(FEV ~ Smoke , data = FEVData)
e.asht <- wmwTest(FEV ~ Smoke, data = FEVData)
e.BT <- BuyseTest(Smoke ~ cont(FEV), add.halfNeutral = TRUE, data =FEVData)

## e.WR <- WRrec(ID = 1:NROW(FEVData),
##               time = FEVData$FEV,
##               status = rep(1,NROW(FEVData)),
##               trt = FEVData$Smoke,
##               strata = rep(1,NROW(FEVData)),
##               naive = TRUE)

## *** point estimate
coef(e.BT, statistic = "favorable")
e.asht$estimate
1/(1+exp(-coef(e.pim)))

## coef(e.BT, statistic = "winRatio")
## exp(e.WR$log.WR) ## small difference

## *** statistical inference
confint(e.BT, statistic = "favorable")
1/(1+exp(-confint(e.pim)))
c(e.asht$conf.int, e.asht$p.value)
## e.WR$pval 

## ** 16.2.2 Generating synthetic data
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

## * 16.3 GPC with a single endpoint
## convert wide format to long format
dtPC.toxL <- as.data.frame(dtPC.toxW, responseName = "Probability")
names(dtPC.toxL)[1:2] <- c("treatment","grade")

## color scale
colorG2R <- scales::seq_gradient_pal(low = rgb(green=0.9,0,0),
                                     high = rgb(red=0.9,0,0))

## left panel
gg.tox <- ggplot(dtPC.toxL, aes(x = treatment, fill = grade, y = Probability))
gg.tox <- gg.tox + geom_bar(position = position_fill(reverse = TRUE),
                            stat = "identity")
gg.tox <- gg.tox + scale_y_continuous(labels = scales::percent)
gg.tox <- gg.tox + scale_fill_manual("Worse adverse event",
                                     values = colorG2R(seq(0,1,length.out=6)))
gg.tox

## right panel
library(prodlim)
plot(prodlim(Hist(OS,statusOS) ~ treatment, data = dt.data))                             

## ** 16.3.2 GPC analysis for completely-observed data
dt.data$toxicity.num <- as.numeric(dt.data$toxicity)

## traditional Wilcoxon test
wilcox.test(toxicity.num ~ treatment, data = dt.data)

## GPC approach
eTox.BT <- BuyseTest(treatment ~ cont(toxicity.num, operator = "<0"),
                     data = dt.data)
summary(eTox.BT)

print(eTox.BT, percentage = FALSE)

## equivalence with wilcoxon test
wilcox.test(toxicity.num ~ treatment, data = dt.data, correct = FALSE)$p.value

eTox.BTvarperm <- BuyseTest(treatment ~ cont(toxicity.num, operator = "<0"),
                 data = dt.data, method.inference = "varexact-permutation")
confint(eTox.BTvarperm)

## ** 16.3.3 Extracting treatment effect measures
confint(eTox.BT, statistic = "favorable", null = 0.5)

eTox.BThalf <- BuyseTest(treatment ~ cont(toxicity.num, operator = "<0"),
                         add.halfNeutral = TRUE,
                         data = dt.data,
                         trace = FALSE)
confint(eTox.BThalf, statistic = "favorable")

confint(eTox.BThalf)

BuyseTest.options(statistic = "winRatio", add.halfNeutral = FALSE)
BuyseTest.options(reinitialise = TRUE)


## ** 16.3.4 Methods for statistical inference
## retrieving bootstrap distribution
eTox.BTboot <- BuyseTest(treatment ~ cont(toxicity.num, operator = "<0"),
                         data = dt.data,
                         method.inference = "bootstrap", seed = 11)
confint(eTox.BTboot, method.ci.resampling = "percentile")
##              estimate         se   lower.ci   upper.ci null p.value
## toxicity.num  -0.0736 0.05690566 -0.1888924 0.02908674    0   0.177
confint(eTox.BTboot, method.ci.resampling = "gaussian", transformation = FALSE)
##              estimate         se  lower.ci   upper.ci null   p.value
## toxicity.num  -0.0736 0.05690566 -0.185133 0.03793304    0 0.1958836
confint(eTox.BT, transformation = FALSE)
##              estimate         se  lower.ci   upper.ci null   p.value
## toxicity.num  -0.0736 0.05617859 -0.183708 0.03650802    0 0.1901594

eTox.BTboot2 <- BuyseTest(treatment ~ cont(toxicity.num, operator = "<0"),
                          data = dt.data, cpus = 7, strata.resampling = "treatment",
                          method.inference = "bootstrap", n.resampling = 1e4, seed = 11)
confint(eTox.BTboot2, method.ci.resampling = "gaussian", transformation = FALSE)
##             estimate         se   lower.ci   upper.ci null   p.value
## toxicity.num  -0.0736 0.05620024 -0.1837504 0.03655045    0 0.1903302

## BuyseTest.options(order.Hprojection = 2)
## eTox.BT2 <- BuyseTest(treatment ~ cont(toxicity.num, operator = "<0"),
##                       data = dt.data, trace = FALSE)
## confint(eTox.BT2, transformation = FALSE)
## BuyseTest.options(order.Hprojection = 1)

## retrieving Wilcoxon test
eTox.BTperm <- BuyseTest(treatment ~ cont(toxicity.num, operator = "<0"),
                         data = dt.data,
                         method.inference = "permutation", n.resampling = 1e5, cpus = 7, seed = 11)
confint(eTox.BTperm)
##              estimate         se lower.ci upper.ci null   p.value
## toxicity.num  -0.0736 0.05593371       NA       NA    0 0.1897881

## ** 16.3.5 Threshold of clinical relevance
eTox.BT2 <- BuyseTest(treatment ~ cont(toxicity.num, threshold = 2, operator = "<0"),
                     data=dt.data, keep.pairScore = TRUE, trace = FALSE)
print(eTox.BT2)

getPairScore(eTox.BT2)

dt.data[c(3:4,201),c("id","treatment","OS","statusOS","toxicity","gender")]

model.tables(eTox.BT, columns = "threshold")

## ** 16.3.6 Stratification
ffG <- treatment ~ cont(toxicity.num, operator = "<0") + strata(gender)
eTox.BTG <- BuyseTest(ffG, data=dt.data, keep.pairScore = TRUE, trace = FALSE)
summary(eTox.BTG)
## summary(eTox.BTG, percentage = FALSE)

getPairScore(eTox.BTG)

confint(eTox.BTG, strata = TRUE)

## ** 16.3.7 Handling right-censoring when assessing efficacy

dt.data[,.(censoring=mean(statusOS==0)),by = "treatment"]

eEff.BT <- BuyseTest(treatment ~ tte(OS, statusOS), data=dt.data,
                     keep.pairScore = TRUE, trace = FALSE)

getPairScore(eEff.BT)[c(1,2,2623,8553),]

print(eEff.BT)

eEff.BT2 <- BuyseTest(treatment ~ tte(OS, statusOS), data=dt.data,
                      scoring.rule = "Gehan", keep.pairScore = TRUE, trace = FALSE)
print(eEff.BT2)

### getPairScore(eEff.BT2)[c(1,2,2623,8553),]

dt30.data <- data.table::copy(dt.data)
dt30.data[OS>30, c("OS", "statusOS") := .(30,0)]

eEff.BT30 <- BuyseTest(treatment ~ tte(OS, statusOS, restriction = 25), data=dt30.data,
                       keep.pairScore = TRUE, trace = FALSE)
print(eEff.BT30)

getPairScore(eEff.BT30)[index.C==44 & index.T == 211,]
dt.data[c(44,211)]

## ** 16.3 [Extra] Analyzing data from cross-over trial using stratification

## *** single comparison within patient
data("vasscoresL", package = "LMMstar")
vasscores <- subset(vasscoresL, group == "AB" & treatment != "C")
vasscores$treatment <- droplevels(vasscores$treatment)
vasscores <- vasscores[order(vasscores$id),]
head(vasscores)

eCO1.BT <- BuyseTest(treatment ~ cont(vas) + id, data = vasscores)
confint(eCO1.BT)
sqrt(var(coef(eCO1.BT, strata = TRUE))/10)

eCO1.BTperm <- BuyseTest(treatment ~ cont(vas) + id, data = vasscores,
                         method.inference = "permutation", strata.resampling = "id",
                         seed = 10)
confint(eCO1.BTperm)

vasB <- vasscores[vasscores$treatment == "B","vas"]
vasA <- vasscores[vasscores$treatment == "A","vas"]
prop.test(x = sum(vasB>vasA), n = length(vasB))$p.value

## *** multiple comparisons within patient
library(mvtnorm)
set.seed(2)
muCO2 <- rep(c(0,0.5),2)
rhoCO2 <- 0.75-0.5*as.matrix(dist(muCO2==0.5))+diag(0.25,4,4)
dfW.CO2 <- data.frame(id = 1:10, rmvnorm(10, mean = muCO2, sigma = rhoCO2))
dfL.CO2 <- reshape(dfW.CO2, direction = "long", idvar = "id",
                   times = rep(c("A","B"),2), timevar = "treatment",
                   varying = paste0("X",1:4), v.names = "vas")
dfL.CO2$vas <- (dfL.CO2$vas - min(dfL.CO2$vas))/diff(range(dfL.CO2$vas))
dfL.CO2$vas <- round(100*dfL.CO2$vas)
dfL.CO2[dfL.CO2$id==1,]

eCO2.BT <- BuyseTest(treatment ~ cont(vas) + id,  data = dfL.CO2)
confint(eCO2.BT)
## confintCO(eCO2.BT)

eCO2.Delta <- coef(eCO2.BT, strata = TRUE)
c(pool = mean(eCO2.Delta),eCO2.Delta)
with(dfW.CO2, 2*((X2>X1) + (X2>X3) + (X4>X1) + (X4>X3))/4 - 1)

sqrt((var(eCO2.Delta))/length(eCO2.Delta))

## *** simulation study - 1 sequence
confintCO <- function(object){ ## object <- eCO.BT
    iCI <- confint(object, strata = object@level.strata)
    iDelta <- mean(iCI$estimate)
    iSe <- sqrt((var(iCI$estimate))/NROW(iCI))
    iOut <- BuyseTest:::confint_Ustatistic(Delta = atanh(iDelta), Delta.se = iSe/(1-iDelta^2), null = 0,
                                           alternative = "two.sided", alpha = 0.05,
                                           endpoint = "Y", backtransform.delta = tanh, backtransform.se = function(x,se){se*(1-tanh(x)^2)})

    return(iOut)
}

warperCO51 <- function(i, n.arm, perm){ ## n.arm <- 50

    treat <- rep(c("A","B"),2)
    mu <- rep(0,4)
    Rho <- 0.75-0.5*as.matrix(dist(treat=="A")) + diag(0.25,4,4)
    dfW <- data.frame(id = 1:n.arm, mvtnorm::rmvnorm(n.arm, mean = mu, sigma = Rho))
    dfL <- reshape(dfW, direction = "long", idvar = "id",
                   times = rep(c("A","B"), 2), timevar = "treatment",
                   varying = 2:5, v.names = "vas")
    BT <- BuyseTest(treatment ~ cont(vas) + id,  data = dfL, trace = FALSE)
    BT.perm <- BuyseTest(treatment ~ cont(vas) + id,  data = dfL, trace = FALSE, method.inference = "permutation", strata.resampling = "id")

    CI.naive <- data.frame(type = "naive", confint(BT))
    CI.corrected <- data.frame(type = "corrected", confintCO(BT))
    if(perm){
        CI.perm <- data.frame(type = "perm", confint(BT.perm))
        return(cbind(seed = i, n = n.arm, rbind(CI.naive, CI.corrected, CI.perm)))
    }else{
        return(cbind(seed = i, n = n.arm, rbind(CI.naive, CI.corrected)))
    }
    
}

library(pbapply)
ls.resCO51 <- pblapply(1:1000, function(i){
    rbind(warperCO51(i, n.arm = 5, perm = FALSE),
          warperCO51(i, n.arm = 10, perm = TRUE),
          warperCO51(i, n.arm = 20, perm = FALSE),
          warperCO51(i, n.arm = 40, perm = FALSE),
          warperCO51(i, n.arm = 100, perm = FALSE))
}, cl = 100)

library(data.table)
dt.resCO51 <- as.data.table(do.call(rbind,ls.resCO51))
dtS.resCO51 <- dt.resCO51[, .(rep = .N, bias = mean(estimate), se = sd(estimate), seNA = sum(is.na(se)), sehat = mean(se, na.rm = TRUE), type1 = mean(p.value<=0.05, na.rm = TRUE)), by = c("n","type")]
dtS.resCO51
##         n      type   rep       bias         se  seNA      sehat      type1
##     <num>    <char> <int>      <num>      <num> <int>      <num>      <num>
##  1:     5     naive  1000  0.0025000 0.37808292     0 0.15117129 0.42500000
##  2:     5 corrected  1000  0.0025000 0.37808292     5 0.35591199 0.08844221
##  3:    10     naive  1000  0.0008000 0.26611121     0 0.11136246 0.39300000
##  4:    10 corrected  1000  0.0008000 0.26611121     0 0.25842224 0.04600000
##  5:    10      perm  1000  0.0008000 0.26611121     0 0.20371183 0.11100000
##  6:    20     naive  1000  0.0026000 0.18265844     0 0.08059622 0.38100000
##  7:    20 corrected  1000  0.0026000 0.18265844     0 0.18408208 0.04600000
##  8:    40     naive  1000 -0.0027125 0.13028487     0 0.05716630 0.39500000
##  9:    40 corrected  1000 -0.0027125 0.13028487     0 0.13044876 0.03800000
## 10:   100     naive  1000 -0.0006400 0.08022859     0 0.03650759 0.37900000
## 11:   100 corrected  1000 -0.0006400 0.08022859     0 0.08238939 0.04000000

## *** simulation study - 2 sequences
simCrossOver <- function(n.arm, n.time, delta, k.sigma, rho.common, rho.diff,
                         seed = NULL){
    if(!is.null(seed)){set.seed(seed)}

    treat <- rep(c("A","B"),5)
    muAB <- c(A = 0, B = delta)[treat]
    muBA <- c(A = 0, B = delta)[rev(treat)]

    Rho <- rho.common-(rho.common-rho.diff)*as.matrix(dist(treat=="A"))
    diag(Rho) <- 1
    OmegaAB <- tcrossprod(rep(c(1,k.sigma),n.time))*Rho
    OmegaBA <- tcrossprod(rep(c(k.sigma,1),n.time))*Rho

    dfW.AB <- data.frame(id = 1:n.arm, sequence = "AB", mvtnorm::rmvnorm(n.arm, mean = muAB, sigma = OmegaAB))
    dfL.AB <- reshape(dfW.AB, direction = "long", idvar = "id",
                      times = rep(c("A","B"),n.time), timevar = "treatment",
                      varying = 3:12, v.names = "Y")
    dfW.BA <- data.frame(id = (n.arm+1):(2*n.arm), sequence = "BA", mvtnorm::rmvnorm(n.arm, mean = muBA, sigma = OmegaBA))
    dfL.BA <- reshape(dfW.BA, direction = "long", idvar = "id",
                      times = rep(c("B","A"),n.time), timevar = "treatment",
                      varying = 3:12, v.names = "Y")
    return(rbind(dfL.AB,dfL.BA))
}


df.CO <- simCrossOver(n.arm = 10, n.time = 5, delta = 0.5, k.sigma = 1.5, rho.common = 0.75, rho.diff = 0.5, seed = 2)
eCO.BT <- BuyseTest(treatment ~ cont(Y) + id,  data = df.CO)
confint(eCO.BT)
confintCO(eCO.BT)

warperCO52 <- function(i, n.arm){ ## n.arm <- 50
    iData <- simCrossOver(n.arm = n.arm, n.time = 5, delta = 0, k.sigma = 1.5, rho.common = 0.75, rho.diff = 0.5)
    iBT <- BuyseTest(treatment ~ cont(Y) + id,  data = iData, trace = FALSE)
    iCI.naive <- data.frame(type = "naive", confint(iBT))
    iCI.corrected <- data.frame(type = "corrected", confintCO(iBT))
    iOut <- cbind(seed = i, n = n.arm*2, rbind(iCI.naive, iCI.corrected))
    return(iOut)      
}


library(pbapply)
ls.resCO52 <- pblapply(1:1000, function(i){
    rbind(warperCO52(i, n.arm = 5),
          warperCO52(i, n.arm = 10),
          warperCO52(i, n.arm = 20),
          warperCO52(i, n.arm = 40))
}, cl = 100)

library(data.table)
dt.resCO52 <- as.data.table(do.call(rbind,ls.resCO52))
dtS.resCO52 <- dt.resCO52[, .(rep = .N, bias = mean(estimate), se = sd(estimate), sehat = mean(se), type1 = mean(p.value<=0.05)), by = c("n","type")]
dtS.resCO52
##        n      type   rep      bias         se      sehat type1
##    <num>    <char> <int>     <num>      <num>      <num> <num>
## 1:    10     naive  1000  0.000248 0.21251767 0.08883136 0.436
## 2:    10 corrected  1000  0.000248 0.21251767 0.21062651 0.057
## 3:    20     naive  1000 -0.002488 0.15057789 0.06299482 0.429
## 4:    20 corrected  1000 -0.002488 0.15057789 0.14980091 0.048
## 5:    40     naive  1000  0.001736 0.10311579 0.04474594 0.397
## 6:    40 corrected  1000  0.001736 0.10311579 0.10612132 0.044
## 7:    80     naive  1000 -0.007290 0.07454888 0.03159631 0.409
## 8:    80 corrected  1000 -0.007290 0.07454888 0.07517525 0.047


## ** 16.3 [Extra] Multiple imputation
dt.data$toxicityNA <- ifelse(rbinom(NROW(dt.data), size = 1, prob = 0.25),NA,dt.data$toxicity)
table(dt.data$toxicityNA, useNA = "always")
dt.data$toxicityNA <- as.factor(dt.data$toxicityNA)

library(mice)
predictorMatrix <-  matrix(0, nrow = NCOL(dt.data), ncol = NCOL(dt.data),
                           dimnames = list(names(dt.data), names(dt.data)))
predictorMatrix["toxicityNA",c("gender","treatment")] <- 1
e.mice <- mice(dt.data, method = "polyreg", predictorMatrix = predictorMatrix, printFlag = FALSE,
               m = 10, maxit = 10, seed = 10)
df.mice <- complete(e.mice, action = "long")
table(df.mice$toxicityNA, useNA = "always")

eTox.MBT <- with(data = e.mice, BuyseTest(treatment ~ cont(toxicity.num, operator = "<0"), trace = FALSE,
                                          data = data.frame(treatment = treatment, toxicity.num = as.numeric(toxicityNA))))
Delta.mice <- do.call(rbind,lapply(eTox.MBT$analyses, confint))
mean(Delta.mice$estimate)
sqrt(mean(Delta.mice$se^2) +(1+1/10)*var(Delta.mice$estimate))

tidy.S4BuyseTest <- function(x, conf.int = FALSE, conf.level = 0.95, ...) {
    out <- confint(x)
    out$statistic <- out$estimate/out$se
    names(out) <- c("estimate","std.error","conf.low","conf.high","null","p.value","statistic")
    out
}
summary(pool(eTox.MBT))


## * 16.4 Benefit-risk analysis using GPC

## ** 16.4.1 Prioritized & non-prioritized analyses
eBRB.BT <- BuyseTest(treatment ~ tte(OS, statusOS) + cont(toxicity.num, operator = "<0"),
                     data=dt.data, trace = FALSE)
print(eBRB.BT)

eRBB.BT <- BuyseTest(treatment ~ cont(toxicity.num, operator = "<0") + tte(OS, statusOS),
                     data=dt.data, trace = FALSE)
print(eRBB.BT)

eNH.BT <- BuyseTest(treatment ~ cont(toxicity.num, operator = "<0") + tte(OS, statusOS),
                    data=dt.data, hierarchical = FALSE, trace = FALSE)
print(eNH.BT)

eRBB.plot <- plot(eRBB.BT)
eNH.plot <- plot(eNH.BT)
ggpubr::ggarrange(eRBB.plot$plot + ggtitle("Prioritized"),
                  eNH.plot$plot + ggtitle("Non-prioritized"),
                  common.legend = TRUE, legend = "bottom")

eRBB.BT2 <- BuyseTest(treatment ~ cont(toxicity.num, operator = "<0", weight = 0.8)
                      + tte(OS, statusOS, weight = 0.2),
                      data=dt.data, hierarchical = FALSE, trace = FALSE)
print(eRBB.BT2)

rbind("prioritized" = confint(eRBB.BT, transform = FALSE, endpoint = 1),
      "non-prioritized" = confint(eNH.BT, transform = FALSE, endpoint = 1))

## ** 16.4.2 Thresholds of clinical relevance

eSH.BT <- BuyseTest(treatment ~ tte(OS, statusOS, threshold = 28)
                              + cont(toxicity.num, operator = "<0", threshold = 2)
                              + tte(OS, statusOS, threshold = 14)
                              + cont(toxicity.num, operator = "<0"),
                    data=dt.data, trace = FALSE)
print(eSH.BT)

## ** 16.4.3 Encoding of the outcome

dt.data$OS2 <- dt.data$OS
dt.data$OS2[dt.data$statusOS==0] <- 150

print(BuyseTest(treatment ~ tte(OS2, statusOS), data=dt.data, trace = FALSE))

eD2.BT <- BuyseTest(treatment ~ bin(statusOS, operator = "<0") + tte(OS2, statusOS), data=dt.data, trace = FALSE)
print(eD2.BT)


## [EXTRA] ignoring some of the neutral pairs, e.g. only compare toxicity among survivors
dt.data$toxicity2 <- dt.data$toxicity.num
dt.data$toxicity2[dt.data$statusOS==1] <- -1

eBRB2.BT <- BuyseTest(treatment ~ bin(statusOS, operator = "<0") + cont(toxicity2, operator = "<0"), data=dt.data, trace = FALSE)
print(eBRB2.BT)
(1947 - 2279)/40000

eR2.BT <- BuyseTest(treatment ~ cont(toxicity2, operator = "<0"),
                    data=dt.data[statusOS==0], trace = FALSE)
print(eR2.BT, percentage = FALSE)

## ** 16.4.4 Sensitivity analyses
eRBB.Se <- sensitivity(eRBB.BT, threshold = list(1:5,c(0,5,10)),
                       band = TRUE, adj.p.value = TRUE, seed = 10, trace = FALSE)
eRBB.Se[c(1,2,6),]

autoplot(eRBB.Se) + facet_wrap(~OS, labeller = label_both)

library(lava)
eRBB.Hdecomp <- iid(eRBB.Se)
dim(eRBB.Hdecomp)

eRBB.cor <- cor(eRBB.Hdecomp)
eRBB.corlower <- eRBB.cor[lower.tri(eRBB.cor)]
range(eRBB.corlower)

range(eRBB.Se$adj.p.value/eRBB.Se$p.value)

e.MBT <- BuyseMultComp(list("OS-tox" = eBRB.BT, "tox-OS" = eRBB.BT, "threshold" = eSH.BT),
                       cluster = "id",
                       seed = 10)
e.MBT

## * 16.5 Power calculation for GPC

## ** 16.5.1 Data-generating mechanism
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

## ** 16.5.2 Simulation-based power and sample size estimation
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
