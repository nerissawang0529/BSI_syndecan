install.packages("pcalg")
install.packages("rcausal")
# Install remotes if not already installed
install.packages("remotes")

# Install rcausal from GitHub
remotes::install_github("username/rcausal", auth_token = "your_personal_access_token")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Rgraphviz")
install.packages("igraph")
install.packages("devtools")
devtools::install_github("rje42/MixedGraphs")
BiocManager::install("RBGL")
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("bd2kccd/rcausal")

library(pcalg)
library(rcausal)

library(Rgraphviz)
library(igraph)
##BiocManager::install("Rgraphviz")
#devtools::install_github("rje42/MixedGraphs")
library(MixedGraphs)

##build in data
data("gmG")

##check which package is used for graph
class(gmG$g)

##variable names
nodes(gmG$g)

##underlying true graph
plot(gmG$g)

##Edge list
gmG$g@edgeL


##change edge list to matrix
gmg_r <-convert(gmG$g,cur_format='graphNEL')

gmg_r$edges$directed

gmg_r <- withAdjMatrix(gmg_r)

gmg_r$edges$directed

gmg_r$vnames

##check the true CPDAG
CPDAG <- ConstructPAG(convert(gmG$g,cur_format='graphNEL'))
CPDAG$vnames <- nodes(gmG$g)
CPDAG


#run PC algorithm

indepTest <- gaussCItest
alpha <- 0.95

pc_result <- pc(suffStat = list(C = cor(gmG$x), n = 8),
                indepTest = indepTest,
                alpha = alpha,
                labels = nodes(gmG$g),  # If colnames are missing, they’ll be auto-assigned
                verbose = TRUE)

plot(pc_result, main = "PC Algorithm Result")

CPDAG



# Run the GES algorithm

#define score function. Currently build-in only has only gaussian
score <- new("GaussL0penObsScore", gmG$x)

ges_result <- ges(score)

ig <- as(ges_result$essgraph, "graphNEL")
nodes(ig) <- nodes(gmG$g)
# Plot the estimated graph
plot(ig, main = "Estimated Graph from GES")


CPDAG



##using r-causal package

#list available algorithms
rcausal::tetradrunner.listAlgorithms()

#run gfci
rcausal_res<- tetradrunner(algoId = 'gfci',df = as.data.frame(gmG$x),
                           testId = 'fisher-z-test',scoreId = 'sem-bic-score',
                           dataType = 'continuous',maxDegree = -1,
                           maxPathLength = -1,completeRuleSetUsed = FALSE,
                           faithfulnessAssumed = TRUE, verbose = FALSE,numberResampling = 5,
                           resamplingEnsemble = 1, addOriginalDataset = TRUE)

#convert result graph's format to Mixedgraph package
rcausal_res_r <- convert(gfci_restoPAG(rcausal_res),cur_format = 'PAG')
rcausal_res_r$vnames <- nodes(gmG$g)
rcausal_res_r


GESMAG_res <- MAGsearch::GESMAG(gmG$x)
GESMAG_res$PAG_no_turn$vnames <- GESMAG_res$PAG$vnames <- nodes(gmG$g)
GESMAG_res
