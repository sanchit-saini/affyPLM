####################################################################
#
# threestep - threestep interface to c code
#
# Copyright (C) 2003-2004    Ben Bolstad
#
# function by B. M. Bolstad <bolstad@stat.berkeley.edu>
#
# Originally based on rma.R from the affy package
# the three step method implemented in c code
#
# this code serves as interface to the c code.
#
#
# note this function does not leave the supplied
# AffyBatch unchanged if you select DESTRUCTIVE=TRUE. this is 
# for memory purposes but can be quite
# dangerous if you are not careful. Use destructive=FALSE if this is
# deemed likely to be a problem. NOTE DESTRUCTIVE item removed for now
#
# Feb 6, 2003 - additional summary methods added
# Mar 22, 2003 - methods so that can use LESN backgrounds
# Mar 24, 2003 - add in ability to use MAS style background
# Jul 23, 2003 - standard errors from three step methods
# Jul 24, 2003 - added scaling as an option
# Jul 26, 2003 - introduced normalization options parameter
#                converted background parameters in same manner
# Jul 27, 2003 - cleaned up parameter list
# Oct  5, 2003 - summary.param which controls options for summarization  added.
# Jan 18, 2004 - remove mapping functions to internalfunctions.R
# Feb 23, 2004 - "subset" parameter is now not ignored
# Aug 05, 2004 - some small changes to deal with new pre-process structure.
# Oct 10, 2006 - add verbosity.level argument to function (remove older and unused verbose argument
#
#####################################################################

threestep <- function(object,subset=NULL, normalize=TRUE,background=TRUE,background.method="RMA.2",normalize.method="quantile",summary.method="median.polish",background.param = list(),normalize.param=list(),summary.param=list(),verbosity.level=0){


  if (!is(object, "AffyBatch")) {
    stop(paste("argument is",class(object),"threestep requires AffyBatch"))
  }

  
  rows <- length(probeNames(object,subset))
  cols <- length(object)


  if (is.null(subset)){
    ngenes <- length(geneNames(object))
  } else {
    ngenes <- length(unique(subset))
  }
  
  
  modelparam <- verify.model.param(object,PM ~ -1 + samples)
  R.model <- PLM.designmatrix3(object,PM ~ -1 + samples,variable.type=c(default="factor"),constraint.type=c(default="contr.treatment"))
    
  background.param <- verify.bg.param(R.model, background.method,background.param = background.param)
  normalize.param <- verify.norm.param(R.model, normalize.method,normalize.param=normalize.param)
  
  s.param <- list(psi.type = "Huber", psi.k = NULL)
  s.param[names(summary.param)] <- summary.param
  if (is.null(s.param$psi.k)) {
    s.param$psi.k <- get.default.psi.k(s.param$psi.type)
  }
  s.param$psi.type <- get.psi.code(s.param$psi.type)

  
  # to avoid having to pass location information to the c code, we will just call the R code method
  if (is.element(background.method,c("MAS","MASIM")) & background){
    if (verbosity.level > 0){
      cat("Background Correcting\n")
    }
    object <- bg.correct.mas(object)
  }
  if (is.element(background.method,c("gcrma","GCRMA")) & background){
    if (verbosity.level > 0){
      cat("Background Correcting\n")
    }
    object <- bg.correct.gcrma(object)
  }
  results <- .Call("R_threestep_c",pm(object,subset), mm(object,subset), probeNames(object,subset), ngenes, normalize, background, background.method, normalize.method, get.summary.code(summary.method),background.param,normalize.param,s.param,verbosity.level, PACKAGE="affyPLM") 
  
  colnames(results[[1]]) <- sampleNames(object)
  colnames(results[[2]]) <- sampleNames(object)
  #se.exprs <- array(NA, dim(exprs)) # to be fixed later, besides which don't believe much in nominal se's with medianpolish
  
  phenodata <- phenoData(object)
  annotation <- annotation(object)
  description <- description(object) 
  notes <- notes(object)
  
  new("exprSet", exprs = results[[1]], se.exprs = results[[2]], phenoData = phenodata, 
       annotation = annotation, description = description, notes = notes)
}


