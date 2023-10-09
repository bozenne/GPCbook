### tableInference-1.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Oct  9 2023 (10:17) 
## Version: 
## Last-Updated: Oct  9 2023 (11:40) 
##           By: Brice Ozenne
##     Update #: 11
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

djack <- sapply(1:20, function(i){
    coef(BuyseTest(treatment ~ cont(score), data = dtInference[-i], trace = FALSE), statistic = "favorable")
})

M.A <- cbind(score = dtInference$score[1:n.data],
             jack = (coef(GPC, statistic = "favorable") - jack)[1:n.data]*(n.data-1),
             h = getIid(GPC, statistic = "favorable", scale = FALSE)[1:n.data])
M.B <- cbind(score = dtInference$score[(n.data+1):(2*n.data)],
             jack = (coef(GPC, statistic = "favorable") - jack)[(n.data+1):(2*n.data)]*(n.data-1),
             h = getIid(GPC, statistic = "favorable", scale = FALSE)[(n.data+1):(2*n.data)])

table1.inference <- xtable(cbind(as.character(round(M.A[,1],digits[1])),
                                 as.character(round(M.A[,2],digits[2])),
                                 as.character(round(M.A[,3],digits[2])),
                                 NA,
                                 as.character(round(M.B[,1],digits[1])),
                                 as.character(round(M.B[,2],digits[2])),
                                 as.character(round(M.B[,3],digits[2]))))
          
attr(table1.inference,"additional") <- c(sigma_HE = mean(M.A[,"h"]^2)/n.data,
                                         sigma_HC = mean(M.B[,"h"]^2)/n.data,
                                         sigma_W = mean(M.A[,"h"]^2)/n.data+mean(M.B[,"h"]^2)/n.data,
                                         GS = confint(GPC, statistic = "favorable")$se^2)



## * display table 1
print(table1.inference, include.rownames=FALSE)
## % latex table generated in R 4.2.0 by xtable 1.8-4 package
## % Mon Oct  9 10:26:42 2023
## \begin{table}[ht]
## \centering
## \begin{tabular}{lllllll}
##   \hline
## 1 & 2 & 3 & 4 & 5 & 6 & 7 \\ 
##   \hline
## -1.2 & 0.44 & 0.44 &  & -0.6 & 0.04 & 0.04 \\ 
##   -0.5 & 0.04 & 0.04 &  & -2.2 & -0.26 & -0.26 \\ 
##   -0.8 & 0.34 & 0.34 &  & -0.7 & 0.04 & 0.04 \\ 
##   0.3 & -0.26 & -0.26 &  & -2.1 & -0.26 & -0.26 \\ 
##   1.1 & -0.26 & -0.26 &  & -1.3 & -0.26 & -0.26 \\ 
##   1.2 & -0.26 & -0.26 &  & -0.4 & 0.24 & 0.24 \\ 
##   0.7 & -0.26 & -0.26 &  & -0.7 & 0.04 & 0.04 \\ 
##   -0.5 & 0.04 & 0.04 &  & -0.9 & -0.06 & -0.06 \\ 
##   0.6 & -0.26 & -0.26 &  & -0.1 & 0.24 & 0.24 \\ 
##   -1.2 & 0.44 & 0.44 &  & -0.3 & 0.24 & 0.24 \\ 
##    \hline
## \end{tabular}
## \end{table}

print(attr(table1.inference, "additional"), include.rownames=FALSE)
## sigma_HE sigma_HC  sigma_W       GS 
##  0.00844  0.00384  0.01228  0.01228 

##----------------------------------------------------------------------
### tableInference-1.R ends here
