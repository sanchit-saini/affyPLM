#############################################################
##
## file: normalize.scaling.R
##
## aim: scaling normalization functions
##
## Written by B. M. Bolstad <bmb@bmbolstad.com>
##
## History
## Aug 23, 2003 - Initial version
##              - Added type argument to normalize.AffyBatch.scaling
## Jul 10, 2006 - add in log.scalefactors (this implements the scale factors computed
##                using log2 scale data as suggested by Lu, Chao (2004) Improving the scaling
##                normalization for high-density oligonucleotide GeneChip expression microarrays
##                BMC Bioinformatics 5:103
##
#############################################################


normalize.scaling <- function(X,trim=0.02,baseline=-1,log.scalefactors=FALSE){
  .Call("R_normalize_scaling",X,trim,baseline,log.scalefactors,PACKAGE="affyPLM")
}


normalize.AffyBatch.scaling <- function(abatch,type=c("together","pmonly","mmonly","separate"),trim=0.02,baseline=-1,log.scalefactors=FALSE){

  type <- match.arg(type)

  if (type == "pmonly"){
    Index <- unlist(indexProbes(abatch,"pm"))
  } else if (type == "mmonly"){
    Index <- unlist(indexProbes(abatch,"mm"))
  } else if (type == "together"){
    Index <- unlist(indexProbes(abatch,"both"))
  } else if (type == "separate"){
    abatch <- normalize.AffyBatch.scaling(abatch,type="pmonly",trim=trim,baseline=baseline,log.scalefactors=log.scalefactors)
    Index <- unlist(indexProbes(abatch,"mm"))
  }
  
  
  col.names <- colnames(exprs(abatch))
  exprs(abatch)[Index,] <- normalize.scaling(exprs(abatch)[Index,],trim=trim,baseline=baseline,log.scalefactors=log.scalefactors)
  colnames(exprs(abatch)) <- col.names
  return(abatch)
}


