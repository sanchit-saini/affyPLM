#############################################################
##
## file: rmaPLM.R
##
## Aim: Implement traditional "RMA" (Background, QN and median
##      polish) as a PLM (probe level model).
##
## Copyright (C) 2003-2005     B. M. Bolstad 
##
## Created by: B. M. Bolstad <bolstad@stat.berkeley.edu>
##
## Created on: Sept 14, 2003
##
## History
## Sept 14, 2003 - Initial verison
## Dec 14, 2003 - Added model.description
## Jan 18, 2004 - removR to C code mapping functions to internalfunctions.R
## Feb 24, 2004 - added a "subset" parameter (and made it active)
## Aug 5, 2004 - some changes to deal with new pre-processing structure
##               Changes to deal with new PLMset structure
## Mar 12, 2005 - fix it so that probe effects are returned
## Oct 10, 2006 - add verbosity.level argument
## Jan 4, 2007 - change how const.coef etc are handled
##
#############################################################

rmaPLM <- function(object,subset=NULL,normalize=TRUE,background=TRUE,background.method="RMA.2",normalize.method="quantile",background.param = list(),normalize.param=list(),output.param=list(),model.param=list(),verbosity.level=0){

  if (!is(object, "AffyBatch")) {
    stop(paste("argument is",class(object),"fitPLM requires AffyBatch"))
  }

  rows <- length(probeNames(object,subset))
  cols <- length(object)


  if (is.null(subset)){
    ngenes <- length(geneNames(object))
  } else {
    ngenes <- length(unique(subset))
  }


  op.param <- list(weights=FALSE,residuals=TRUE, pseudo.SE=FALSE,resid.SE=FALSE)
  op.param[names(output.param)] <- output.param
  #output <- verify.output.param(output.param)
  modelparam <- verify.model.param(object,PM ~ -1 + probes + samples, model.param=model.param)
  R.model <- PLM.designmatrix3(object,PM ~ -1 + probes + samples, variable.type=list(default="factor"),constraint.type=list(default="contr.sum"))

  
  b.param <- verify.bg.param(R.model, background.method,background.param = background.param)
  n.param <- verify.norm.param(R.model, normalize.method,normalize.param=normalize.param)

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
  md.param <- list(psi.type = "Huber",psi.k=NULL)
  md.param[names(model.param)] <- model.param
  
  if (is.null(md.param$psi.k)){
    md.param$psi.k <- get.default.psi.k(md.param$psi.type)
  }
  md.param$psi.type <- get.psi.code(md.param$psi.type)
  
  fit.results <- .Call("R_rmaPLMset_c",pm(object,subset), mm(object,subset), probeNames(object,subset), ngenes, normalize, background, background.method, normalize.method,b.param,n.param,op.param,md.param,verbosity.level,PACKAGE="affyPLM") #, PACKAGE="affyPLM") 
  probenames <- rownames(pm(object,subset))
  colnames(fit.results[[1]]) <- sampleNames(object)
  colnames(fit.results[[4]]) <- sampleNames(object)
  if (op.param$weights){
    colnames(fit.results[[3]]) <- sampleNames(object)
    rownames(fit.results[[3]]) <- probenames
  }
  
  if (op.param$residuals){
    colnames(fit.results[[8]]) <- sampleNames(object)
    rownames(fit.results[[8]]) <- probenames 
  }

  colnames(fit.results[[2]]) <- "ProbeEffects"
  rownames(fit.results[[2]]) <- probenames
  colnames(fit.results[[5]]) <- "SEProbeEffects"
  rownames(fit.results[[5]]) <- probenames
  
  fit.results[[6]] <- matrix(0,0,0)
  fit.results[[7]] <- matrix(0,0,0)

  
  phenodata <- phenoData(object)
  annotation <- annotation(object)
  description <- description(object)
  notes <- notes(object)

  
  new("PLMset",
  chip.coefs=fit.results[[1]],
      probe.coefs= split(fit.results[[2]],probeNames(object,subset)),
      weights=list(PM.weights=fit.results[[3]],MM.weights=matrix(0,0,0)),
      se.chip.coefs=fit.results[[4]],
      se.probe.coefs=split(fit.results[[5]],probeNames(object,subset)),
      const.coefs=fit.results[[6]],
      se.const.coefs=fit.results[[7]],
      residuals=list(PM.resid=fit.results[[8]],MM.resid=matrix(0,0,0)),
      residualSE=fit.results[[9]],
      varcov = fit.results[[10]],
      phenoData = phenodata,
      annotation = annotation,
      experimentData= description,
      ##FIXME: remove # after notes is fixed
                                        # x@notes = notes,
      cdfName=object@cdfName,
      nrow=object@nrow,
      ncol=object@ncol,
      model.description = list(which.function="rmaPLM",
        preprocessing=list(bg.method=background.method,bg.param=b.param,
          background=background,norm.method=normalize.method,
          norm.param=n.param,normalize=normalize),
        modelsettings =list(model.param=md.param,summary.method=NULL,model=NULL,
          constraint.type=NULL,variable.type=NULL),
        outputsettings=op.param, 
        R.model=R.model)
 )

}
