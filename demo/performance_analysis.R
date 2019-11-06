library(senlm)
library(proftools)
library(Rgraphviz)

## --- Create a dataset --- #
x <- c(-2, -1, 1, 2)
y <- c(1, 1, 1, 1)
 
## --- Fit a model whilst recording profiling information --- #
Rprof(tmp <- tempfile(), line.profiling = TRUE)
fitted_model <- senlm(x, y, mean_fun="gaussian", err_dist="gaussian")
Rprof(NULL)

## --- Convert the profiling information into proftools format --- #
pd <- readProfileData(tmp)

## --- Plot Profiling Data --- #
plotProfileCallGraph(pd)
##calleeTreeMap(pd)
##flameGraph(pd)

## --- Annotate Sourcecode with Profiling Data -- #
annotateSource(pd)
##srcSummary(pd)

## --- Summarise the Profiling Data --- #
pathSummary(pd)
##funSummary(pd)
##callSummary(pd)
##hotPaths(pd)
##flatProfile(pd)
