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
if(is.na(iter_sim)){iter_sim <- 1}
if(is.na(n.iter_sim)){n.iter_sim <- 10}
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
warper <- function(i, n, mu, sigma, n.resampling){
    if(length(n)==1){
        n <- rep(n,2)
    }

    data <- simBuyseTest(n.T = n[1], n.C = n[2], argsCont = list(mu.C = 0, mu.T = mu, sigma.C = 1, sigma.T = sigma))
    tps.U <- system.time(
        BT.U <- BuyseTest(treatment ~ cont(score), data = data, trace = FALSE,
                          method.inference = "u-statistic")
    )
    tps.boot <- system.time(
        BT.boot <- BuyseTest(treatment ~ cont(score), data = data, n.resampling = n.resampling[1], trace = FALSE,
                             method.inference = "studentized bootstrap", strata.resampling = "treatment")
    )
    tps.perm <- system.time(
        BT.perm <- BuyseTest(treatment ~ cont(score), data = data, n.resampling = n.resampling[1], trace = FALSE,
                             method.inference = "studentized permutation")
    )
  
    out <- rbind(
        cbind(method = "Ustat", statistic = "netBenefit", confint(BT.U, statistic = "netBenefit", transform = FALSE), time = tps.U["elapsed"]),
        cbind(method = "Ustat-trans", statistic = "netBenefit", confint(BT.U, statistic = "netBenefit", transform = TRUE), time = tps.U["elapsed"]),
        cbind(method = "perm-perc", statistic = "netBenefit", confint(BT.perm, statistic = "netBenefit", method.ci.resampling = "percentile"), time = NA),
        cbind(method = "perm-stud", statistic = "netBenefit", confint(BT.perm, statistic = "netBenefit", method.ci.resampling = "studentized"), time = tps.boot["elapsed"]),
        cbind(method = "boot-perc", statistic = "netBenefit", confint(BT.boot, statistic = "netBenefit", method.ci.resampling = "percentile"), time = NA),
        cbind(method = "boot-stud", statistic = "netBenefit", confint(BT.boot, statistic = "netBenefit", method.ci.resampling = "studentized"), time = tps.perm["elapsed"]),
        cbind(method = "Ustat", statistic = "winRatio", confint(BT.U, statistic = "winRatio", transform = FALSE), time = tps.U["sys.self"]),
        cbind(method = "Ustat-trans", statistic = "winRatio", confint(BT.U, statistic = "winRatio", transform = TRUE), time = tps.U["sys.self"]),
        cbind(method = "perm-perc", statistic = "winRatio", confint(BT.perm, statistic = "winRatio", method.ci.resampling = "percentile"), time = NA),
        cbind(method = "perm-stud", statistic = "winRatio", confint(BT.perm, statistic = "winRatio", method.ci.resampling = "studentized"), time = tps.boot["elapsed"]),
        cbind(method = "boot-perc", statistic = "winRatio", confint(BT.boot, statistic = "winRatio", method.ci.resampling = "percentile"), time = NA),
        cbind(method = "boot-stud", statistic = "winRatio", confint(BT.boot, statistic = "winRatio", method.ci.resampling = "studentized"), time = tps.perm["elapsed"])
    )

    ## if(length(n.resampling)==2 && n.resampling[2]<n.resampling[1]){
    ##     BT.boot0 <- BT.boot
    ##     BT.boot0@DeltaResampling <- BT.boot@DeltaResampling[1:n.resampling[2],,,drop=FALSE]
    ##     BT.boot0@covarianceResampling <- BT.boot@covarianceResampling[1:n.resampling[2],,,drop=FALSE]
    ##     BT.boot0@weightStrataResampling <- BT.boot@weightStrataResampling[1:n.resampling[2],,drop=FALSE]
    ##     BT.boot0@n.resampling <- n.resampling[2]
    ##     out <- rbind(
    ##         out,
    ##         cbind(method = "boot-perc0", statistic = "netBenefit", confint(BT.boot0, statistic = "netBenefit", method.ci.resampling = "percentile"), time = NA),
    ##         cbind(method = "boot-stud0", statistic = "netBenefit", confint(BT.boot0, statistic = "netBenefit", method.ci.resampling = "studentized"), time = NA)
    ##     )
    ## }

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
    for(iGrid in 1:NROW(grid.sim)){ ## iGrid <- 1
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

	
