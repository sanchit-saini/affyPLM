#############################################################
##
## file: normalize.scaling.R
##
## aim: scaling normalization functions
##
## Written by B. M. Bolstad <bolstad@stat.berkeley.edu>
##
## History
## Aug 23, 2003 - Initial version
##              - Added type argument to 
##
#############################################################


normalize.scaling <- function(X,trim=0.02,baseline=-1){
  .Call("R_normalize_scaling",X,trim,baseline)
}


normalize.AffyBatch.scaling <- function(abatch,type=c("together","pmonly","mmonly","separate"),trim=0.02,baseline=-1){

  type <- match.arg(type)

  if (type == "pmonly"){
    Index <- unlist(indexProbes(abatch,"pm"))
  } else if (type == "mmonly"){
    Index <- unlist(indexProbes(abatch,"mm"))
  } else if (type == "together"){
    Index <- unlist(indexProbes(abatch,"both"))
  } else if (type == "separate"){
    abatch <- normalize.AffyBatch.scaling(abatch,type="pmonly",trim=trim,baseline=baseline)
    Index <- unlist(indexProbes(abatch,"mm"))
  }
  
  
  col.names <- colnames(exprs(abatch))
  exprs(abatch)[Index,] <- normalize.scaling(exprs(abatch)[Index,],trim,baseline)
  colnames(exprs(abatch)) <- col.names
  return(abatch)
}


