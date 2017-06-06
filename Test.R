#setwd("~/PolyGraph")
setwd("~/Dropbox/Research/PolyGraph")


#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
simidx <- args[1]

source("PolyGraphFunc.R")
source("PolyGraph.R")

runmode <- "alphas"
branchcandidate <- NaN
numsteps <- 1000000
numsample <- 1000


# Good combo for SimpleTree2 with branches of 0.02 length (spike and slab prior, i.e. ssfactor=25)
#selbranch <- c("B","q")
#selbpasted <- paste(selbranch,collapse="_")
#innerfreqs_proposize <- 0.3
#alpha_prior_stdev <- 0.1
#alpha_proposize <- 0.04
#ssfactor <- 25
#branchlength <- 0.02
#modegraph <- "simpletree2"


# Good combo for SimpleTree2 with branches of 0.05 length (spike and slab prior, i.e. ssfactor=25)
#selbranch <- c("B","q")
#selbpasted <- paste(selbranch,collapse="_")
#innerfreqs_proposize <- 0.3
#alpha_prior_stdev <- 0.1
#alpha_proposize <- 0.04
#ssfactor <- 25
#branchlength <- 0.05
#modegraph <- "simpletree2"




# Simple Graph - branch 0.02 (spike and slab prior, i.e. ssfactor=25)
selbranch <- c("v","q")
selbpasted <- paste(selbranch,collapse="_")
innerfreqs_proposize <- 0.1
alpha_prior_stdev <- 0.1
alpha_proposize <- 0.02
ssfactor <- 25
branchlength <- 0.02
modegraph <- "simplegraph"


# Simple Graph - branch 0.05 (spike and slab prior, i.e. ssfactor=25)
#selbranch <- c("v","q")
#selbpasted <- paste(selbranch,collapse="_")
#innerfreqs_proposize <- 0.1
#alpha_prior_stdev <- 0.1
#alpha_proposize <- 0.03
#ssfactor <- 25
#branchlength <- 0.05
#modegraph <- "simplegraph"




selalpha <- 0.1
numproc <- 5
qbfactor <- 3

simulate <- TRUE
  
if(modegraph == "simplegraph"){source("SimpleGraph.R")
} else if(modegraph == "simpletree2"){source("SimpleTree2.R")}
  
test <- mclapply(1:numproc, function(runnum) {
  tracefile <- paste("trace_file_alphas_",modegraph,"_branch",branchlength,"_",selbpasted,"_sim",simidx,"_run",runnum,".txt",sep="")
  print(tracefile)
  PolyGraph(tracefile,leaves_counts,neut_leaves_counts,effects,runmode,NaN,numsteps,numsample,innerfreqs_proposize,alpha_prior_stdev,alpha_proposize,ssfactor,qbfactor)
})


#selalpha <- 0.2
#allstats <- c()
#for(i in seq(1,50)){
#  simulate <- TRUE
#  source("SimpleGraph.R")
#  test <- ChiSquared(supergraph,leaves_freqs,effects,neut_leaves_freqs,total=FALSE)
#  print(i)
#  print(test)
#  allstats <- rbind(allstats,test)
#}

#plot(density(allstats[,1]),xlim=c(0,40))
#lines(seq(0.1,40,0.1),dchisq(seq(0.1,40,0.1),df=1),col="red")


#plot(hist(allstats[,3]))


