#############################################################
##
## file: threestepPLM.R
##
## Aim: Implement threestep expression summaries in PLM framework
##
##
## Copyright (C) 2003     B. M. Bolstad
##
## Created by: B. M. Bolstad <bolstad@stat.berkeley.edu>
##
## Created on: Oct 9, 2003
##
## History
## Oct 9, 2003 - Initial verison
## Dec 14, 2003 - added model.description
##
#############################################################




threestepPLM <- function(object,normalize=TRUE,background=TRUE,background.method="RMA.2",normalize.method="quantile",summary.method="median.polish",background.param = list(),normalize.param=list(),output.param=list(),model.param=list()){

  get.background.code <- function(name) {
    background.names <- c("RMA.1", "RMA.2", "IdealMM","MAS","MASIM","LESN2","LESN1","LESN0","SB")
    if (!is.element(name, background.names)) {
      stop(paste(name, "is not a valid summary method. Please use one of:",
                 "RMA.1", "RMA.2", "IdealMM","LESN2","LESN1","LESN0","SB"))
    }
    code <- c(1, 2, 3, 4, 5, 6, 7, 8,9)[name == background.names]
    code
  }

  get.normalization.code <- function(name){
    normalization.names <- c("quantile","quantile.probeset","scaling")

    if (!is.element(name,normalization.names)){
      stop(paste(name,"is not a valid summary method. Please use one of:","quantile","quantile.probeset","scaling"))
    }
    code <- c(1,2,3)[name == normalization.names]
    code
  }

  convert.LESN.param <- function(param.list){
    defaults <- c(0.25,4)
    if (!is.null(param.list$baseline)){
      defaults[1] <- param.list$baseline
    }
    if (!is.null(param.list$theta)){
      defaults[2] <- param.list$theta
    }
    defaults
  }


  get.psi.code <- function(name){
    psi.names <- c("Huber","fair","Cauchy","Geman-Mclure","Welsch","Tukey","Andrews")
    if (!is.element(name, psi.names)) {
      stop(paste(name, "is not a valid Psi type. Please use one of:",
                 "Huber","fair","Cauchy","Geman-Mclure","Welsch","Tukey","Andrews"))
    }
    code <- c(0:6)[name == psi.names]
    code
  }

  get.default.psi.k <- function(name){
    psi.code <- get.psi.code(name)
    ## ** Huber - k = 1.345
    ## ** Fair - k = 1.3998
    ## ** Cauchy - k=2.3849
    ## ** Welsch - k = 2.9846
    ## ** Tukey Biweight - k = 4.6851
    ## ** Andrews Sine - K = 1.339
    if (psi.code == 0){
      psi.k <- 1.345
    } else if (psi.code == 1){
      psi.k <- 1.3998
    } else if (psi.code == 2){
      psi.k <- 2.3849
    } else if (psi.code == 4){
      psi.k <- 2.9846
    } else if (psi.code == 5){
       psi.k <- 4.6851
    } else if (psi.code == 6){
      psi.k <- 1.339
    } else {
      psi.k <- 1
    }
    psi.k
  }

   get.summary.code <- function(name){
    summary.names <- c("median.polish","tukey.biweight","average.log","rlm","log.average","log.median","median.log","log.2nd.largest","lm")

    if (!is.element(name,summary.names)){
      stop(paste(name,"is not a valid summary method. Please use one of:","median.polish","tukey.biweight","average.log","rlm","log.average","log.median","median.log","log.2nd.largest","lm"))
    }
    code <- c(1,2,3,4,5,6,7,8,9)[name ==summary.names]
    code
  }



  rows <- length(probeNames(object))
  cols <- length(object)

  ngenes <- length(geneNames(object))

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

  fit.results <- .Call("R_threestepPLMset_c",pm(object), mm(object), probeNames(object), ngenes, normalize, background, get.background.code(background.method), get.normalization.code(normalize.method),b.param,n.param,op.param,md.param, PACKAGE="affyPLM")
  probenames <- rownames(pm(object))
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
