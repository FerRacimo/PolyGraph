#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

source("PolyGraphFunc.R")

file <- args[1]
qfile <- args[2]
phenoname <- args[3]
graphscript <- args[4]
output_boxplot <- args[5]
output_phenograph <- args[6]
output_phenograph_qfile <- args[7]


# Load admixture graph
source(graphscript)

# Load table
table <- read.table(file,header=TRUE)

# Remove burn-in
trace <- table[seq(10,dim(table)[1]-1,1),]

# Obtain alpha parameters
alphacols <- grep("alpha_",colnames(trace))

# Make pheno-graph from trace output
pdf(output_phenograph,width = 8, height = 8)
minsel <- -0.3
maxsel <- 0.3
MakeGraphPlot(file,edgevalues,"r",phenoname,minsel,maxsel)
dev.off()

# Make pheno-graph from qfile output
pdf(output_phenograph_qfile,width = 8, height = 8)
minsel <- -10
maxsel <- 10
MakeGraphPlotQfile(qfile,edgevalues,"r",phenoname,minsel,maxsel)
dev.off()


  
# Make boxplot
SNPs <- max(sapply(colnames(table),function(x){as.numeric(paste(unlist(strsplit(gsub("[^0-9]", "", x),"")),collapse=""))}),na.rm=TRUE)
titleboxplot <- paste("Posterior distributions of alpha parameters: ",phenoname," (num. SNPs: ",SNPs,")",sep="")
pdf(output_boxplot,width=15,height=5)
finalboxplot <- AlphaBoxPlot(trace,titleboxplot,FALSE,0)
grid.arrange(finalboxplot,ncol=1,nrow=1)
dev.off()
  