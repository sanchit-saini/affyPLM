#############################################################
##
## file: threestepPLM.R
##
## Aim: Implement threestep expression summaries in PLM framework
##   
##
## Copyright (C) 2003-2004     B. M. Bolstad 
##
## Created by: B. M. Bolstad <bolstad@stat.berkeley.edu>
##
## Created on: Oct 9, 2003
##
## History
## Oct 9, 2003 - Initial verison
## Dec 14, 2003 - added model.description
## Jan 18, 2003 - remove internal functions
## Feb 24, 2003 - Add subset as a parameter (and make it active)
##
#############################################################




threestepPLM <- function(object,subset=NULL,normalize=TRUE,background=TRUE,background.method="RMA.2",normalize.method="quantile",summary.method="median.polish",background.param = list(),normalize.param=list(),output.param=list(),model.param=list()){

  rows <- length(probeNames(object,subset))
  cols <- length(object)

  if (is.null(subset)){
    ngenes <- length(geneNames(object))
  } else {
    ngenes <- length(unique(subset))
  }
  
  
  # background correction for RMA type backgrounds
  bg.dens <- function(x){density(x,kernel="epanechnikov",n=2^14)}


  # to avoid having to pass location information to the c code, we will just call the R code method
  if (is.element(background.method,c("MAS","MASIM")) & background){
    cat("Background Correcting\n")
    object <- bg.correct.mas(object)
  }
    
  LESN.param <-list(baseline=0.25,theta=4)
  LESN.param <- convert.LESN.param(LESN.param)

  b.param <- list(densfun =  body(bg.dens), rho = new.env(),lesnparam=LESN.param)
  b.param[names(background.param)] <- background.param
  
  n.param <- list(scaling.baseline=-4,scaling.trim=0.0,use.median=FALSE,use.log2=TRUE)
  n.param[names(normalize.param)] <- normalize.param
  
  
  op.param <- list(weights=FALSE,residuals=TRUE, pseudo.SE=FALSE,resid.SE=FALSE)
  op.param[names(output.param)] <- output.param

  md.param <- list(psi.type = "Huber",psi.k=NULL,summary.code=get.summary.code(summary.method))
  md.param[names(model.param)] <- model.param
  
  if (is.null(md.param$psi.k)){
    md.param$psi.k <- get.default.psi.k(md.param$psi.type)
  }
  md.param$psi.type <- get.psi.code(md.param$psi.type)
  
  fit.results <- .Call("R_threestepPLMset_c",pm(object,subset), mm(object,subset), probeNames(object,subset), ngenes, normalize, background, get.background.code(background.method), get.normalization.code(normalize.method),b.param,n.param,op.param,md.param, PACKAGE="affyPLM") 
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

  #colnames(fit.results[[2]]) <- "ProbeEffects"
  #rownames(fit.results[[2]]) <- probenames
  #colnames(fit.results[[5]]) <- "SEProbeEffects"
  #rownames(fit.results[[5]]) <- probenames
  
  fit.results[[6]] <- matrix(0,0,0)
  fit.results[[7]] <- matrix(0,0,0)

  
  phenodata <- phenoData(object)
  annotation <- annotation(object)
  description <- description(object)
  notes <- notes(object)

  
  x <- new("PLMset")
  x@chip.coefs=fit.results[[1]]
  x@probe.coefs= fit.results[[2]]
  x@weights=fit.results[[3]]
  x@se.chip.coefs=fit.results[[4]]
  x@se.probe.coefs=fit.results[[5]]
  x@exprs=fit.results[[6]]
  x@se.exprs=fit.results[[7]]
  x@residuals=fit.results[[8]]
  x@residualSE=fit.results[[9]]
  x@varcov = fit.results[[10]]
  x@phenoData = phenodata
  x@annotation = annotation
  x@description = description
  x@notes = notes
  x@cdfName=object@cdfName
  x@nrow=object@nrow
  x@ncol=object@ncol

  
  x@model.description = list(which.function="rmaPLM",preprocessing=list(bg.method=background.method,bg.param=b.param,background=background,norm.method=normalize.method,norm.param=n.param,normalize=normalize),modelsettings =list(model.param=md.param,summary.method=summary.method,model=NULL,constraint.type=NULL,variable.type=NULL),outputsettings=op.param)
  

  
  x

}
