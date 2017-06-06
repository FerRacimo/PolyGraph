library("admixturegraph")

# Define standard admixture graph
leaves <- c("YRI", "CEU", "CHB","PEL","CLM")
inner_nodes <- c("r","q","t","u","s","v")
edges <- parent_edges(c(edge("YRI","r"),
                        edge("q","r"),
                        edge("v","q"),
                        edge("CEU","v"),
                        edge("s","q"),
                        edge("CHB","s"),
                        edge("t","s"),
                        edge("PEL","t"),
                        edge("u","t"),
                        edge("CLM","u"),
                        edge("u","v"),
                        admixture_edge("u","v","t")))
admixtures <- admixture_proportions(c(
  admix_props("u", "v", "t", "gamma")
))
graph <- agraph(leaves, inner_nodes, edges, admixtures)

# Define admixture values
admvalues <- rbind(c("gamma",0.7655))

if(length(admvalues) > 0){
  colnames(admvalues) <- c("ratename","rate")
  admvalues <- as.data.frame(admvalues)
  admvalues[,2] <- as.numeric(as.character(admvalues[,2]))
}

# Define branch drift lengths
edgevalues <- rbind(c("YRI", "r",0.020052),
                    c("q", "r",0.2076),
                    c("v","q",0.028),
                    c("CEU","v",0.036),
                    c("s","q",0.049003),
                    c("CHB","s",0.090152),
                    c("t","s",0.069),
                    c("PEL","t",0.001),
                    c("u","t",0.001),
                    c("CLM","u",0.001),
                    c("u","v",0.001))


colnames(edgevalues) <- c("child","parent","value")
edgevalues <- as.data.frame(edgevalues)
edgevalues[,3] <- as.numeric(as.character(edgevalues[,3]))

# SUPERGRAPH OBJECTS HAS: 1) A GRAPH, 2) EDGE VALUES, 3) ADMIXTURE VALUES
supergraph <- list(graph,edgevalues,admvalues)

# Order edges topologically
supergraph <- OrderBranches(supergraph)

source("LoadFiles.R")