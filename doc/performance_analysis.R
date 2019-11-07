## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ------------------------------------------------------------------------
library(senlm)
library(proftools)
library(Rgraphviz)

## --- Fit a model whilst recording profiling information --- #
Rprof(tmp <- tempfile(), line.profiling = TRUE)
fitted_model <- senlm (data=haul, xvar="x_depth", yvar="y_Sebastolobus.altivelis", 
              mean_fun="gaussian", err_dist="zip")
Rprof(NULL)

## --- Convert the profiling information into proftools format --- #
pd <- readProfileData(tmp)

## --- Plot Profiling Data --- #
plotProfileCallGraph(pd)
##calleeTreeMap(pd)
##flameGraph(pd)

## --- Summarise the Profiling Data --- #
pathSummary(pd)
##funSummary(pd)
##callSummary(pd)
##hotPaths(pd)
##flatProfile(pd)

## --- Annotate Sourcecode with Profiling Data -- #
# annotateSource is very useful, but won't work within a vignette.
#annotateSource(pd, show = TRUE)
##srcSummary(pd)

