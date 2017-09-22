# PolyGraph

<img src="https://github.com/FerRacimo/PolyGraph/blob/master/HEIGHT_1KG_YRI_CEU_CHB_PEL_CLM.png" height="300">

# Required libraries

PolyGraph requires the following R libraries. Make sure they are installed before you run it:
- optparse
- admixturegraph
- msm
- reshape2
- pscl
- parallel
- ggplot2
- gridExtra

# Input files

To run, PolyGraph requires 3 input files:
- a GWAS input file
- a neutral input file
- a graph parameter file

# GWAS input file format

The GWAS input file should contain the SNPs that are known to be significantly associated with a particular trait. This file should contain 4+ columns. GWAS_HEIGHT_1000genomes_allpops.txt is an example GWAS input file. The first 3 columns correspond to the chromosome name, physical position and SNP id of each SNP that has been associated with a particular trait.

The fourth column contains the estimated effect size for the derived allele of that SNP. The estimated effect size (beta) can be computed as z * sqrt(v), where z is the signed z-score and v is the variance of the effect size, both of which can be obtained from a standard GWAS summary statistic file. It is important that the sign of the effect size is polarized with respect to the derived allele (rather than the minor or alternative allele, as may be traditionally done in some GWAS), so knowlege of the ancestral allele is required.

Finally the fifth columns and all columns after that one contain the number of ancestral and derived alleles in each population that is a leaf in our graph (in the format "[number of ancestral alleles],[number of derived alleles]"). All the leaf populations specificed in the graph parameter file should be included as columns in this file (in any order), with the same column name as their corresponding name in the parameter file. It is ok to include additional populations that are not specified in the parameter file (the program will just ignore these).

# Neutral input file format

The neutral input file should contain unlinked neutral SNPs that are not significantly associated with the trait in question, as this file will be used to compute the neutral population covariance matrix. This file has the same format as the GWAS file, except that the effect size column has no relevance and can be filled with any symbol. Neut_HEIGHT_1000genomes_allpops.txt is an example neutral input file (it is available in this repository but must be unzipped first).

# Graph parameter file

This file contains information about the admixture graph, which should have been previously inferred using a program like MixMapper (Lipson et al. 2013) or qpGraph (Patterson et al. 2012). SimpleGraph.R is an example graph parameter file. This file uses the R library admixturegraph (Leppälä et al. 2017) to build the graph.

# Running PolyGraph

The command line options for PolyGraph are as follows. The first 3 options (input file names) are required. We recommend tweaking the size of the proposal steps to obtain adequate MCMC acceptance probabilities.

    Rscript RunPhenotype.R \
    -w [GWAS input file name; default: NULL] \
    -e [Neutral input file name; default: NULL] \
    -r [Graph R file name; default: NULL] \
    -o [MCMC trace output file name; default: MCMC_trace.txt] \
    -q [Q_B statistic output file name; default: qfile.txt] \
    -n [Total number of steps in MCMC; default: 1000000] \
    -x [Sample (print) every X steps from the MCMC run; default: 1000] \
    -i [Size of std dev of proposal distribution for frequencies of inner nodes; default: 0.1] \
    -s [Standard deviation for alpha prior; default: 0.1] \
    -t [Size of std dev of proposal distribution for alpha parameters; default: 0.02] \
    -u [In spike-and-slab prior, number by which std dev of wide Normal dist. will be divived to obtain std dev of narrow Normal dist.; default: 25]
    -f [Number by which maximum Q_B score will be divided to obtain Q_B cutoff. If equal to 0, then switch to using the appropriate chi-squared significance cutoff instead; default: 3]


Here is an example of a full command line:

    Rscript RunPhenotype.R \
    -g GWAS_HEIGHT_1000genomes_allpops.txt \
    -e Neut_HEIGHT_1000genomes_allpops.txt \
    -r 1KG_YRI_CEU_CHB_PEL_CLM.R \
    -o trace_HEIGHT_1KG_YRI_CEU_CHB_PEL_CLM.txt \
    -q qfile_HEIGHT_1KG_YRI_CEU_CHB_PEL_CLM.txt \
    -n 1000000 \
    -x 1000 \
    -i 0.1 \
    -s 0.1 \
    -t 0.02 \
    -u 25 \
    -f 3


# Output files

PolyGraph should produce two output files:
- The list of test statistics for each branch in the graph (fast to compute). This contains:
    - The Q_B statistic for each branch (unsigned, chi-squared distributed with 1 degree of freedom under neutrality)
    - The q_B statistic for each branch(signed, Normal(0,1) distributed under neutrality)
    - The P-value for the Q_B statistic for each branch computed using the chi-squared distribution.
    - The P-value for the Q_B statistics for each branch computed using 1,000 pseudo-replicates of the trait-significant SNPs, by randomizing the sign of the estimated effect sizes.
    - The Q_X statistic (Berg et al. 2014) and its P-value (under row "Total"), which serves to test for significant evidence for selection in the entire set of populations
- The trace from the MCMC run (slower to compute; may take hours or days, depending on the complexity of the graph and the number of trait-affecting SNPs)

# Plotting output

To plot the output of an MCMC run in both a boxplot and a "poly-graph" (for both the alpha parameters of the trace and the q_b statistics), use the Plot_Trace.R script:

    Rscript Plot_Trace.R [trace_file.txt] [qfile.txt] [name_of_phenotype] [admixture_graph_file.R] [name_of_boxplot.pdf] [name_of_polygraph_from_mcmc.pdf] [name_of_polygraph_from_qb_statistics.pdf]
