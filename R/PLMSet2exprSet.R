#############################################################
##
## file: PLMset2exprSet
##
## aim: Convert a PLMset to ExpressionSet.
##
## Copyright (C) 2003   Ben Bolstad
##
## created by: B. M. Bolstad <bolstad@stat.berkeley.edu>
## created on: Sept 12, 2003 
##
##  Description: This is often useful
##      because people want to use exprSets and many
##      bioconductor functions work on exprSets
##
##      eventually this function will perform more specific
##      actions depending on what the specific model
##      that was actually fitted. However, for the initial
##      implementation we assume that the default model
##      (probe and sample effects) has been fit.
##     
##
## History
## Sept 12, 2003 - Initial version
## Nov 2, 2003 - fix error
##
##############################################################



PLMset2exprSet <- function(pset){
 new("ExpressionSet",
    exprs =  coefs(pset),
    se.exprs =  se(pset),
    phenoData = phenoData(pset),
    annotation = annotation(pset),
    experimentData = description(pset))
  ##FIXME: update this when notes is fixed
  # eset@notes <- notes(pset)
}


pset2eset <- function(pset){
  PLMset2exprSet(pset)
}
