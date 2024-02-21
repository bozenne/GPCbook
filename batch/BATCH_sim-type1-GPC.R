## * Header 
## cd /projects/biostat01/people/hpl802/GPC/book/
## path <- "x:/GPC/book/"
## setwd(path)
## source("BATCH_sim-type1-GPC.R")
## sbatch -a 1-1 -J 'sim-type1-GPC' --output=/dev/null --error=/dev/null R CMD BATCH --vanilla BATCH_sim-type1-GPC.R /dev/null 

rm(list = ls())
gc()

## * seed
iter_sim <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
n.iter_sim <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_COUNT"))
if(is.na(iter_sim)){iter_sim <- 241}
if(is.na(n.iter_sim)){n.iter_sim <- 1000}
cat("iteration ",iter_sim," over ",n.iter_sim,"\n",sep="")

set.seed(1)
seqSeed <- sample(1:max(1e5,n.iter_sim),size=n.iter_sim,replace=FALSE)
iSeed <- seqSeed[iter_sim]
set.seed(iSeed)

cat("seed: ",iSeed,"\n")

## * prepare export
path <- "."
path.res <- file.path(path,"Results","sim-type1-GPC")
if(dir.exists(path.res)==FALSE){
    if(dir.exists(file.path(path,"Results"))==FALSE){
    dir.create(file.path(path,"Results"))
    }
    dir.create(path.res)
}
path.output <- file.path(path,"output","sim-type1-GPC")
if(dir.exists(path.output)==FALSE){
    if(dir.exists(file.path(path,"output"))==FALSE){
    dir.create(file.path(path,"output"))
    }
    dir.create(path.output)
}


## * libraries
library(BuyseTest)
BuyseTest.options(order.Hprojection = 2)
warper <- function(i, n, mu, sigma, n.resampling){
    if(length(n)==1){
        n <- rep(n,2)
    }

    data <- simBuyseTest(n.T = n[1], n.C = n[2],
                         argsCont = list(mu.C = 0, mu.T = mu, sigma.C = 1, sigma.T = sigma))
    data$category <- as.numeric(cut(data$score, breaks = c(-Inf,-3,-2,-1,-0,1,2,3,Inf)))

    ## continuous
    tps.CU <- system.time(
        BT.CU <- BuyseTest(treatment ~ cont(score), data = data, trace = FALSE,
                           method.inference = "u-statistic", add.halfNeutral = TRUE)
    )
    tps.Cboot <- system.time(
        BT.Cboot <- BuyseTest(treatment ~ cont(score), data = data, n.resampling = n.resampling[1], trace = FALSE,
                              method.inference = "studentized bootstrap", strata.resampling = "treatment", add.halfNeutral = TRUE)
    )
    tps.Cperm <- system.time(
        BT.Cperm <- BuyseTest(treatment ~ cont(score), data = data, n.resampling = n.resampling[1], trace = FALSE,
                              method.inference = "studentized permutation", add.halfNeutral = TRUE)
    )
  
    ## categorical
    tps.FU <- system.time(
        BT.FU <- BuyseTest(treatment ~ cont(category), data = data, trace = FALSE,
                           method.inference = "u-statistic", add.halfNeutral = TRUE)
    )
    tps.Fboot <- system.time(
        BT.Fboot <- BuyseTest(treatment ~ cont(category), data = data, n.resampling = n.resampling[1], trace = FALSE,
                              method.inference = "studentized bootstrap", strata.resampling = "treatment", add.halfNeutral = TRUE)
    )
    tps.Fperm <- system.time(
        BT.Fperm <- BuyseTest(treatment ~ cont(category), data = data, n.resampling = n.resampling[1], trace = FALSE,
                              method.inference = "studentized permutation", add.halfNeutral = TRUE)
    )

    out <- rbind(
        ## Net benefit (continuous outcome)
        cbind(outcome = "continuous", method = "Ustat", statistic = "netBenefit", confint(BT.CU, statistic = "netBenefit", transform = FALSE), time = tps.CU["elapsed"]),
        cbind(outcome = "continuous", method = "Ustat-trans", statistic = "netBenefit", confint(BT.CU, statistic = "netBenefit", transform = TRUE), time = tps.CU["elapsed"]),
        cbind(outcome = "continuous", method = "perm-perc", statistic = "netBenefit", confint(BT.Cperm, statistic = "netBenefit", method.ci.resampling = "percentile"), time = NA),
        cbind(outcome = "continuous", method = "perm-stud", statistic = "netBenefit", confint(BT.Cperm, statistic = "netBenefit", method.ci.resampling = "studentized"), time = tps.Cboot["elapsed"]),
        cbind(outcome = "continuous", method = "boot-perc", statistic = "netBenefit", confint(BT.Cboot, statistic = "netBenefit", method.ci.resampling = "percentile"), time = NA),
        cbind(outcome = "continuous", method = "boot-basic", statistic = "netBenefit", confint(BT.Cboot, statistic = "netBenefit", method.ci.resampling = "gaussian"), time = NA),
        cbind(outcome = "continuous", method = "boot-stud", statistic = "netBenefit", confint(BT.Cboot, statistic = "netBenefit", method.ci.resampling = "studentized"), time = tps.Cperm["elapsed"]),
        ## Win ratio (continuous outcome)
        cbind(outcome = "continuous", method = "Ustat", statistic = "winRatio", confint(BT.CU, statistic = "winRatio", transform = FALSE), time = tps.CU["sys.self"]),
        cbind(outcome = "continuous", method = "Ustat-trans", statistic = "winRatio", confint(BT.CU, statistic = "winRatio", transform = TRUE), time = tps.CU["sys.self"]),
        cbind(outcome = "continuous", method = "perm-perc", statistic = "winRatio", confint(BT.Cperm, statistic = "winRatio", method.ci.resampling = "percentile"), time = NA),
        cbind(outcome = "continuous", method = "perm-stud", statistic = "winRatio", confint(BT.Cperm, statistic = "winRatio", method.ci.resampling = "studentized"), time = tps.Cboot["elapsed"]),
        cbind(outcome = "continuous", method = "boot-perc", statistic = "winRatio", confint(BT.Cboot, statistic = "winRatio", method.ci.resampling = "percentile"), time = NA),
        cbind(outcome = "continuous", method = "boot-basic", statistic = "winRatio", confint(BT.Cboot, statistic = "winRatio", method.ci.resampling = "gaussian"), time = NA),
        cbind(outcome = "continuous", method = "boot-stud", statistic = "winRatio", confint(BT.Cboot, statistic = "winRatio", method.ci.resampling = "studentized"), time = tps.Cperm["elapsed"]),
        ## Net benefit (categorical outcome)
        cbind(outcome = "categorical", method = "Ustat", statistic = "netBenefit", confint(BT.FU, statistic = "netBenefit", transform = FALSE), time = tps.FU["elapsed"]),
        cbind(outcome = "categorical", method = "Ustat-trans", statistic = "netBenefit", confint(BT.FU, statistic = "netBenefit", transform = TRUE), time = tps.FU["elapsed"]),
        cbind(outcome = "categorical", method = "perm-perc", statistic = "netBenefit", confint(BT.Fperm, statistic = "netBenefit", method.ci.resampling = "percentile"), time = NA),
        cbind(outcome = "categorical", method = "perm-stud", statistic = "netBenefit", confint(BT.Fperm, statistic = "netBenefit", method.ci.resampling = "studentized"), time = tps.Fboot["elapsed"]),
        cbind(outcome = "categorical", method = "boot-perc", statistic = "netBenefit", confint(BT.Fboot, statistic = "netBenefit", method.ci.resampling = "percentile"), time = NA),
        cbind(outcome = "categorical", method = "boot-basic", statistic = "netBenefit", confint(BT.Fboot, statistic = "netBenefit", method.ci.resampling = "gaussian"), time = NA),
        cbind(outcome = "categorical", method = "boot-stud", statistic = "netBenefit", confint(BT.Fboot, statistic = "netBenefit", method.ci.resampling = "studentized"), time = tps.Fperm["elapsed"]),
        ## Win ratio (categorical outcome)
        cbind(outcome = "categorical", method = "Ustat", statistic = "winRatio", confint(BT.FU, statistic = "winRatio", transform = FALSE), time = tps.FU["sys.self"]),
        cbind(outcome = "categorical", method = "Ustat-trans", statistic = "winRatio", confint(BT.FU, statistic = "winRatio", transform = TRUE), time = tps.FU["sys.self"]),
        cbind(outcome = "categorical", method = "perm-perc", statistic = "winRatio", confint(BT.Fperm, statistic = "winRatio", method.ci.resampling = "percentile"), time = NA),
        cbind(outcome = "categorical", method = "perm-stud", statistic = "winRatio", confint(BT.Fperm, statistic = "winRatio", method.ci.resampling = "studentized"), time = tps.Fboot["elapsed"]),
        cbind(outcome = "categorical", method = "boot-perc", statistic = "winRatio", confint(BT.Fboot, statistic = "winRatio", method.ci.resampling = "percentile"), time = NA),
        cbind(outcome = "categorical", method = "boot-basic", statistic = "winRatio", confint(BT.Fboot, statistic = "winRatio", method.ci.resampling = "gaussian"), time = NA),
        cbind(outcome = "categorical", method = "boot-stud", statistic = "winRatio", confint(BT.Fboot, statistic = "winRatio", method.ci.resampling = "studentized"), time = tps.Fperm["elapsed"])
    )

    return(cbind(i = i, m = n[1], n = n[2], mu = mu, sigma  = sigma, out))
}


## * settings
rep.sim <- 25
n.resampling <- 1e4##c(1e4,1e3)

grid.sim <- expand.grid(n = c(10,20,35,50,75,100,150,200),
                        rho = 1,
                        mu = c(0,1),
                        sigma = 2)

## * function to execute
res <- NULL
for(iSim in 1:rep.sim){ ## iSim <- 1
    cat(iSim," ")
    for(iGrid in 1:NROW(grid.sim)){ ## iGrid <- 2
        iRes <- try(warper(iSim,
                           n = c(grid.sim[iGrid,"n"],grid.sim[iGrid,"rho"]*grid.sim[iGrid,"n"]),
                           mu = grid.sim[iGrid,"mu"],
                           sigma = grid.sim[iGrid,"sigma"],
                           n.resampling = n.resampling))
        if(!inherits(iRes,"try-error")){
            res <- rbind(res, iRes)
        }

        ## timing: (25*10*50)/3600
    }
    cat("\n")

    saveRDS(res, file = file.path(path.res,paste0("simul_",iter_sim,"(tempo).rds")))
}

## * export
saveRDS(res, file = file.path(path.res,paste0("simul_",iter_sim,".rds")))

## * R version
print(sessionInfo())

## system.time(
##     iRes <- warper(iSim,
##                    n = c(grid.sim[iGrid,"n"],grid.sim[iGrid,"rho"]*grid.sim[iGrid,"n"]),
##                    mu = grid.sim[iGrid,"mu"],
##                    sigma = grid.sim[iGrid,"sigma"],
##                    n.resampling = n.resampling)
## )

## * gather and process results
if(FALSE){




}

	
