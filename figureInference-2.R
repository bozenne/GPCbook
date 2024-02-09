### FigureInference-2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Oct  9 2023 (10:17) 
## Version: 
## Last-Updated: feb  9 2024 (15:46) 
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
library(xtable)
digits <- c(1,3)

## * simulate data
set.seed(10)
n.data <- 10
dtInference <- simBuyseTest(n.data)
dtInference[, score := round(score,1)]

## * generate table 1

GPC <- BuyseTest(treatment ~ cont(score), data = dtInference, trace = FALSE)
GPCfav.iid <- getIid(GPC, statistic = "favorable", scale = FALSE)

GPCfav.iid[dtInference$treatment=="C"]
 ## [1]  0.44  0.04  0.34 -0.26 -0.26 -0.26 -0.26  0.04 -0.26  0.44
GPCfav.iid[dtInference$treatment=="T"]
 ## [1]  0.04 -0.26  0.04 -0.26 -0.26  0.24  0.04 -0.06  0.24  0.24

GPCfav.iid[dtInference$treatment=="C"]^2
 ## [1] 0.1936 0.0016 0.1156 0.0676 0.0676 0.0676 0.0676 0.0016 0.0676 0.1936
GPCfav.iid[dtInference$treatment=="T"]^2
 ## [1] 0.0016 0.0676 0.0016 0.0676 0.0676 0.0576 0.0016 0.0036 0.0576 0.0576

tapply(GPCfav.iid,
       dtInference$treatment=="C",
       crossprod)
## FALSE  TRUE 
## 0.384 0.844 


##----------------------------------------------------------------------
### FigureInference-2.R ends here
