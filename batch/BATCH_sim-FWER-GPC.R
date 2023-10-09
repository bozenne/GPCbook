## * Header 
## cd /projects/biostat01/people/hpl802/GPC/book/
## path <- "x:/GPC/book/"
## setwd(path)
## source("BATCH_sim-FWER-GPC.R")
## sbatch -a 1-1 -J 'sim-FWER-GPC' --output=/dev/null --error=/dev/null R CMD BATCH --vanilla BATCH_sim-FWER-GPC.R /dev/null 

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
path.res <- file.path(path,"Results","sim-FWER-GPC")
if(dir.exists(path.res)==FALSE){
    if(dir.exists(file.path(path,"Results"))==FALSE){
    dir.create(file.path(path,"Results"))
    }
    dir.create(path.res)
}
path.output <- file.path(path,"output","sim-FWER-GPC")
if(dir.exists(path.output)==FALSE){
    if(dir.exists(file.path(path,"output"))==FALSE){
    dir.create(file.path(path,"output"))
    }
    dir.create(path.output)
}

## * libraries
library(BuyseTest)
warper <- function(i, n, mu, sigma, vec.threshold){
  if(length(n)==1){
    n <- rep(n,2)
  }
  data <- simBuyseTest(n.T = n[1], n.C = n[2], argsCont = list(mu.C = 0, mu.T = mu, sigma.C = 1, sigma.T = sigma))
  
  n.threshold <- length(vec.threshold)
  iBT <- BuyseTest(treatment ~ cont(score) + bin(toxicity), data = data, trace = FALSE,
                   method.inference = "u-statistic")
  iSe <- sensitivity(iBT, threshold = list(vec.threshold,0), band = TRUE, adj.p.value = TRUE, trace = FALSE)
  iSe$bonf.p.value <- pmin(1,iSe$p.value*n.threshold)
  iSe$lower.bonf <- tanh(atanh(iSe$estimate) + qnorm(0.025/n.threshold) * iSe$se/(1-iSe$estimate^2))
  iSe$upper.bonf <- tanh(atanh(iSe$estimate) + qnorm(1-0.025/n.threshold) * iSe$se/(1-iSe$estimate^2))
  colnames(iSe)[1] <- "threshold"
  out <- cbind(i = i, n = n[1], mu = mu, sigma = sigma, iSe[which.min(iSe$p.value)[1],,drop=FALSE])
  return(out)
}

## * settings
rep.sim <- 125
vec.threshold <- seq(0,2, length.out = 10)

grid.sim <- expand.grid(n = c(10,20,35,50,75,100,150,200),
                        mu = c(0,1),
                        sigma = 2)

## * function to execute
res <- NULL
for(iSim in 1:rep.sim){ ## iSim <- 1
    cat(iSim," ")
    for(iGrid in 1:NROW(grid.sim)){ ## iGrid <- 1
        iRes <- try(warper(iSim,
                           n = c(grid.sim[iGrid,"n"],grid.sim[iGrid,"n"]),
                           mu = grid.sim[iGrid,"mu"],
                           sigma = grid.sim[iGrid,"sigma"],
                           vec.threshold = vec.threshold))
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

## * gather and process results
if(FALSE){
path <- "x:/GPC/book/"
setwd(path)

path.sim-FWER-GPC <- file.path("Results","sim-FWER-GPC")
allRes.tempo <- butils::sinkDirectory(path.sim-FWER-GPC, string.keep = "tempo")
allRes.final <- butils::sinkDirectory(path.sim-FWER-GPC, string.exclude = "tempo")
}

	
