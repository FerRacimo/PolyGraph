# Plot graph
#plot(supergraph[[1]], show_inner_node_labels = TRUE, show_admixture_labels = TRUE)


if(exists("gwasfile")){
  
  # Load GWAS data  
  print("Loading data...")
  filename <- gwasfile
  data <- LoadCounts(filename, leaves)
  leaves_counts <- as.data.frame(data[,seq(5,dim(data)[2])])
  leaves_freqs <- ObtainFreqs(leaves_counts)
  effects <- as.numeric(data[,4])
  
  # Load Neutral data
  neutfilename <- neutfile
  neutdata <- LoadCounts(neutfilename, leaves)
  neut_leaves_counts <- as.data.frame(neutdata[,seq(5,dim(neutdata)[2])])
  neut_leaves_freqs <- ObtainFreqs(neut_leaves_counts)

  # Deconstruct graph
  graphedges <- supergraph[[2]]
  deconsgraph <- DeconstructGraph(supergraph)
  leaves <- supergraph[[1]]$leaves

  # Compute empirical covariance matrix
  checksegneut <- which( apply(neut_leaves_freqs,1,sum)/dim(neut_leaves_freqs)[2] < 0.99  & apply(neut_leaves_freqs,1,sum)/dim(neut_leaves_freqs)[2] > 0.01 )
  neut_leaves_freqs <- neut_leaves_freqs[checksegneut,]

  # Compute neutral hetero
  neut_leaves_freqs_means <- apply(neut_leaves_freqs, 1, mean)
  mean_hetero <- neut_leaves_freqs_means*(1-neut_leaves_freqs_means)
  numSNPs <- length(neut_leaves_freqs_means)

  # Compute F matrix
  print("Computing F matrix...")
  Fmat <- sapply(seq(1,dim(neut_leaves_freqs)[2]),function(x){
       sapply(seq(1,dim(neut_leaves_freqs)[2]),function(y){
	cov(neut_leaves_freqs[,x]/sqrt(mean_hetero), neut_leaves_freqs[,y]/sqrt(mean_hetero))
	})
       })
  colnames(Fmat) <- colnames(neut_leaves_freqs)
  rownames(Fmat) <- colnames(neut_leaves_freqs)

  # Compute contributions of each branch to each leaf
  contribmat <- c()
  contribmat_names <- c()
  for(branchidx in seq(1,dim(graphedges)[1])){
    targetbranch <- c(as.character(graphedges[branchidx,1]),as.character(graphedges[branchidx,2]))
    name <- paste(targetbranch,collapse="_")
    branchvec <- CollectBranches(supergraph,targetbranch,deconsgraph,leaves)
    contribmat <- rbind(contribmat,branchvec)
    contribmat_names <- rbind(contribmat_names,name)
  }
  rownames(contribmat) <- contribmat_names
  colnames(contribmat) <- leaves

  # Standardize vectors to be of unit length
  contribmat <- t(apply(contribmat,1, function(x) {
    x <- x - mean(x)
    return(x / norm_vec(x))
  }))


  # Calculate chi-squared statistics
  print("Computing Q_b statistics...")
  stats <- ChiSquaredReduced(graphedges,contribmat,Fmat,leaves_freqs,effects,total=FALSE,randomize=FALSE)
  qtab <- cbind(rownames(stats),round(stats[,1],3),round(stats[,2],3),stats[,3])
  colnames(qtab) <- c("branch","Q_B","q_B","Pval")
  totalstat <- ChiSquaredReduced(graphedges,contribmat,Fmat,leaves_freqs,effects,total=TRUE,randomize=FALSE)
  qtab <- rbind(qtab, cbind("Total",round(totalstat[1],3),round(totalstat[2],3),totalstat[3]))

  # For later use
  pvaltotal <- totalstat[3]
    
  # Calculate sign-randomized chi-squared statistics
  pseudorep <- 1000
  print(paste("Computing sign-randomized P-values, using ",pseudorep," pseudo-replicates...",sep=""))
  stats <- ChiSquaredReduced(graphedges,contribmat,Fmat,leaves_freqs,effects,total=FALSE,randomize=TRUE)
  totalstat <- ChiSquaredReduced(graphedges,contribmat,Fmat,leaves_freqs,effects,total=TRUE,randomize=TRUE)
  allstats <- c( round(stats[,1],3), round(totalstat[1],3) )
  names(allstats)[length(names(allstats))] <- "Total"
  randqtab <- allstats
  for(i in seq(2,pseudorep)){
    stats <- ChiSquaredReduced(graphedges,contribmat,Fmat,leaves_freqs,effects,total=FALSE,randomize=TRUE)
    totalstat <- ChiSquaredReduced(graphedges,contribmat,Fmat,leaves_freqs,effects,total=TRUE,randomize=TRUE)
    allstats <- c( round(stats[,1],3), round(totalstat[1],3) )
    randqtab <- cbind(randqtab, allstats)
  }
  
  # Calculate sign-randomized p-values
  randpvals <- sapply( seq(1, length(as.numeric(qtab[,2]))), function(i) { sum(randqtab[i,] > as.numeric(qtab[,2])[i]) / dim(randqtab)[2] } )
  qtab <- cbind(qtab,randpvals)
  colnames(qtab)[length(colnames(qtab))] <- "SignRandPval"
  
  # Write chi-squared statistics table
  write.table(qtab,file=qfile,quote=FALSE,col.names=TRUE,row.names=FALSE)

}
