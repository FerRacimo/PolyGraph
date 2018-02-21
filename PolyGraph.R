PolyGraph <- function(tracefile,leaves_counts,neut_leaves_counts,effects,runmode,branchcandidate,numsteps,numsample,innerfreqs_proposize,alpha_prior_stdev,alpha_proposize,ssfactor,qbfactor,rootprior="unif"){

#print(branchcandidate) 
print("Running MCMC...")
  
# Run mode
#runmode <- "neutral"; branchcandidate <- NaN
#runmode <- "alphas"; branchcandidate <- NaN
#runmode <- "onealpha"; branchcandidate <- c("v","q")

# Set up prior parameters
innerfreqs_prior_lower <- 0
innerfreqs_prior_upper <- 1
alpha_prior_mean <- 0

# Input
input_graph <- supergraph
edgevalues <- supergraph[[2]]
input_leaf_counts <- leaves_counts
input_leaf_freqs <- ObtainFreqs(input_leaf_counts)
input_leaf_der <- as.matrix(as.data.frame(apply(input_leaf_counts,c(1,2),function(x){as.numeric(strsplit(x,",")[[1]][2])})))
input_leaf_tot <- as.matrix(as.data.frame(apply(input_leaf_counts,c(1,2),function(x){splitted <- strsplit(x,",")[[1]]; return(as.numeric(splitted[1])+as.numeric(splitted[2]))})))
numSNPs <- dim(input_leaf_freqs)[1]

# Neutral SNPs for F matrix
neut_leaf_freqs <- ObtainFreqs(neut_leaves_counts)

# Break input graph
broken_graph <- break_graph(input_graph[[1]])
num_inner_nodes <- length(broken_graph$inner_nodes)
num_nonmig_edges <- length(broken_graph$edges)
num_mig_edges <- length(broken_graph$admixtures)*2
num_total_edges <- num_nonmig_edges + num_mig_edges

# Starting values - inner freqs
start_innerfreqs <- jitter(cbind( input_leaf_freqs ,t(sapply(seq(1,numSNPs), function(x) {AncestralAverage(input_leaf_freqs[x,],input_graph)}))),0.5)
start_innerfreqs <- RemoveFixed(start_innerfreqs)
if(dim(start_innerfreqs)[1] == 1){start_innerfreqs <- t(start_innerfreqs)}

# Starting values - alphas
start_alphas <- cbind(as.character(edgevalues[,1]),as.character(edgevalues[,2]),rep(0,dim(input_graph[[2]])[1]))
colnames(start_alphas) <- c("child","parent","value")
num_alphas <- dim(start_alphas)[1]

# Spike-and-slab proportion (assuming only sparse selected branches)
propss_dev <- 0.1
minpropss <- 0.6
maxpropss <- 0.8
propss <- (minpropss + maxpropss)/2
old_propss <- propss
  
# Deconstruct graph
graphedges <- supergraph[[2]]
deconsgraph <- DeconstructGraph(input_graph)

# Compute chi-squared statistics
teststats <- ChiSquared(supergraph,input_leaf_freqs,effects,neut_leaf_freqs,total=FALSE)
Qstat <- teststats[,1]

# Initialize MCMC
old_innerfreqs <- start_innerfreqs
old_alphas <- start_alphas

# Set-up initial probabilities
log_prob_old_binom <- pbinom(input_leaf_der,input_leaf_tot,old_innerfreqs)
log_prob_old_graph <- pgraph(deconsgraph,old_innerfreqs,effects,old_alphas)

# Root frequency prior
if(rootprior == "unif"){ log_prob_old_innerfreqs <- 0
} else if(rootprior == "beta"){ log_prob_old_innerfreqs <- log_beta_prior(old_innerfreqs[,which(colnames(old_innerfreqs) == "r")],2,2)}
    
# Alpha prior
log_prob_old_alphas <- log_alphas_prior_spikeslab(as.numeric(old_alphas[,3]),alpha_prior_mean,alpha_prior_stdev,ssfactor,old_propss)

numaccepted_freqs <- 0
totalmoves_freqs <- 0
numaccepted_alphas <- 0
totalmoves_alphas <- 0
acceptance_rate_freqs <- NaN
acceptance_rate_alphas <- NaN
genvaluesvec <- c()

# Make header for trace file
header <- "step"
header <- c(header,as.vector(sapply(c(broken_graph$leaves,broken_graph$inner_nodes), function(x) {return( paste(x,seq(1,numSNPs,1),sep="") )})))
header <- c(header,sapply(colnames(old_innerfreqs),function(x){paste("nodevalue",x[1],sep="_")}))
header <- c(header,sapply(colnames(old_innerfreqs),function(x){paste("nodevar",x[1],sep="_")}))
header <- c(header,apply(old_alphas[,c(1,2)],1,function(x){paste("alpha",x[1],x[2],sep="_")}))
header <- c(header,"logprob","logprior")
header <- c(header,"propss")
header <- c(header,"acceptance_rate_freqs")
header <- c(header,"acceptance_rate_alphas")
header <- paste(as.vector(header),collapse="\t")
write(header,file=tracefile,append=FALSE)

# RUN MCMC
for( i in seq(1,numsteps)){

  # Propose new frequencies
  new_innerfreqs <- propose_innerfreqsB(old_innerfreqs,innerfreqs_proposize,numSNPs)
  colnames(new_innerfreqs) <- colnames(old_innerfreqs)
    
  freqtranslogdiff <- translogdiff_innerfreqs(old_innerfreqs,new_innerfreqs,innerfreqs_proposize,numSNPs)

  # Evaluate likelihood of new frequencies
  log_prob_new_binom <- pbinom(input_leaf_der,input_leaf_tot,new_innerfreqs)
  log_prob_new_graph <- pgraph(deconsgraph,new_innerfreqs,effects,old_alphas)
  
  # Root frequency prior
  if(rootprior == "unif"){ log_prob_new_innerfreqs <- 0
  } else if(rootprior == "beta"){ log_prob_new_innerfreqs <- log_beta_prior(new_innerfreqs[,which(colnames(new_innerfreqs) == "r")],2,2)}
    
  log_prob_trans_freq <- apply(freqtranslogdiff, 1, function(y){
    return(sum(y))
  })
  log_prob_new <- log_prob_new_binom + log_prob_new_graph + log_prob_new_innerfreqs
  log_prob_old <- log_prob_old_binom + log_prob_old_graph + log_prob_old_innerfreqs
    
  logpostdiff <-  log_prob_new - log_prob_old + log_prob_trans_freq
  logunifsample <- log(runif(numSNPs,0,1))

  replacevec <- (logpostdiff > logunifsample)

  # Update frequencies
  old_innerfreqs[which(replacevec),] <- new_innerfreqs[which(replacevec),]
  log_prob_old_binom[which(replacevec)] <- log_prob_new_binom[which(replacevec)]
  log_prob_old_graph[which(replacevec)] <- log_prob_new_graph[which(replacevec)]
  log_prob_old_innerfreqs <- log_prob_new_innerfreqs
  numaccepted_freqs <- numaccepted_freqs + sum(replacevec)
  totalmoves_freqs <- totalmoves_freqs + numSNPs

  acceptance_rate_freqs <- numaccepted_freqs / totalmoves_freqs
  
  
  if(runmode != "neutral"){
    
    # Propose new alphas
    if(runmode == "alphas"){
      new_alphas <- propose_alphasB(old_alphas,alpha_proposize,Qstat,qbfactor)
    } else if(runmode == "onealpha"){
      alphaidx <- which(old_alphas[,1] == branchcandidate[1] & old_alphas[,2] == branchcandidate[2])
      alphatochange <- old_alphas[alphaidx,]
      alphatochange <- matrix(alphatochange,ncol=3)
      new_candidate_alpha <- propose_alphas(alphatochange,alpha_proposize)
      new_alphas <- old_alphas
      new_alphas[alphaidx,] <- new_candidate_alpha
    }

    # Evaluate alpha move
    log_prob_new_graph <- pgraph(deconsgraph,old_innerfreqs,effects,new_alphas)
    log_prob_new_alphas <- log_alphas_prior_spikeslab(as.numeric(new_alphas[,3]),alpha_prior_mean,alpha_prior_stdev,ssfactor,old_propss)
    log_prob_new <- sum(log_prob_old_binom) + sum(log_prob_new_graph) + log_prob_new_alphas
    log_prob_old <- sum(log_prob_old_binom) + sum(log_prob_old_graph) + log_prob_old_alphas
    logpostdiff <- log_prob_new - log_prob_old
    unifsample <- runif(1,0,1)
    logunifsample <- log(runif(1,0,1))
    replace <- (logpostdiff > logunifsample)
    if(replace == TRUE){
      old_alphas <- new_alphas
      log_prob_old_graph <- log_prob_new_graph
      log_prob_old_alphas <- log_prob_new_alphas
      log_prob_old <- log_prob_new
    }
    numaccepted_alphas <- numaccepted_alphas + sum(replace)
    totalmoves_alphas <- totalmoves_alphas + 1
    acceptance_rate_alphas <- numaccepted_alphas / totalmoves_alphas
  }
  
  # Propose propss
  new_propss <- rtnorm(1,old_propss,propss_dev,lower=minpropss,upper=maxpropss)
  
  # Calculate propss transition log difference
  propsstranslogdiff <- log(dtnorm(old_propss,new_propss,propss_dev,lower=minpropss,upper=maxpropss)) - log(dtnorm(new_propss,old_propss,propss_dev,lower=minpropss,upper=maxpropss))
  
  #Evaluate propss move (hyperprior for propss is Unif[0,1], so no need to calculate hyperprior difference)
  log_prob_alphas_old_propss <- log_alphas_prior_spikeslab(as.numeric(old_alphas[,3]),alpha_prior_mean,alpha_prior_stdev,ssfactor,old_propss)
  log_prob_alphas_new_propss <- log_alphas_prior_spikeslab(as.numeric(old_alphas[,3]),alpha_prior_mean,alpha_prior_stdev,ssfactor,new_propss)
  logpropssdiff <- log_prob_alphas_new_propss  - log_prob_alphas_old_propss + propsstranslogdiff
  logunifsample <- log(runif(1,0,1))
  replace <- (logpropssdiff > logunifsample)
  if(replace == TRUE){
    old_propss <- new_propss
  }
  
  # Record current prior
  log_old_prior <- sum(log_prob_old_innerfreqs) + sum(log_prob_old_alphas)

  # Print out values
  genvalues <- apply(old_innerfreqs,2,function(x){sum(2*x*effects)})
  genvar <- apply(old_innerfreqs,2, function(x){sum(2*x*(1-x)*effects^2)})
  finalvec <- c(i,as.vector(old_innerfreqs),genvalues,genvar,as.numeric(old_alphas[,3]),sum(log_prob_old_graph)+sum(log_prob_old_binom),log_old_prior,old_propss,acceptance_rate_freqs,acceptance_rate_alphas)
  if( (i %% numsample) == 0){
    write(paste(finalvec,collapse="\t"),file=tracefile,append=TRUE)
    #print(finalvec)
  }
  #print(c(i,genvalues,genvar,as.numeric(old_alphas[,3]),sum(log_prob_old_graph)+sum(log_prob_old_binom),log_old_prior,old_propss,acceptance_rate_freqs,acceptance_rate_alphas))
}

}

