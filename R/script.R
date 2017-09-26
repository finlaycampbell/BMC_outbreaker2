##===== Load functions, libraries and data =====#
source("R/internals.R")
load("data/ances.RData")
load("data/chains.RData")
libs('ape', 'outbreaker2', 'TransPhylo', 'o2mod.transphylo', 'viridis',
     'distcrete', 'phybreak', 'tidyr', 'dplyr', 'ggplot2')


##===== Run analysis yourself ======##
## 1. Simulate outbreak using phybreak
## 2. Reconstruct it using o2mod.TransPhylo
## 3. Reconstruct it using TransPhylo
## 4. Reconstruct it using outbreaker2
analysis <- run.analysis(runs = 1, size = 10)

## 5. Create a summary object called "store"
##       store$acc   - a dataframe with the accuracy of outbreak reconstruction
##       store$simil - a dataframe with the similarity of the inference to
##                     that made by TransPhylo
##       store$times - a dataframe with the time taken per reconstruction
store <- mk.summary(r = analysis)


##===== Or load stored results =====##
store <- create.store(dir = "data/runs/", load = T, dl = F)


##===== Visualise plots =====##
## Figure 1
vis.chains(chains$o2mod.res[[1]], chains$trans.res[[1]])

## Figure 2
vis.anc(ances$o2mod.res[[1]], ances$trans.res[[1]])
