####################################################################
#
# threestep - threestep interface to c code
#
# Copyright (C) 2003    Ben Bolstad
#
# function by B. M. Bolstad <bolstad@stat.berkeley.edu>
#
# based on rma.R from the affy package
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
#
#####################################################################

threestep <- function(object,subset=NULL, verbose=TRUE,normalize=TRUE,background=TRUE,background.method="RMA.2",normalize.method="quantile",summary.method="median.polish",background.param = list(),normalize.param=list(),summary.param=list()){

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


  get.summary.code <- function(name){
    summary.names <- c("median.polish","tukey.biweight","average.log","rlm","log.average","log.median","median.log","log.2nd.largest","lm")

    if (!is.element(name,summary.names)){
      stop(paste(name,"is not a valid summary method. Please use one of:","median.polish","tukey.biweight","average.log","rlm","log.average","log.median","median.log","log.2nd.largest","lm"))
    }
    code <- c(1,2,3,4,5,6,7,8,9)[name ==summary.names]
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

  s.param <- list(psi.type="Huber",psi.k=NULL)
  s.param[names(summary.param)] <- summary.param

  if (is.null(s.param$psi.k)){
    s.param$psi.k <- get.default.psi.k(s.param$psi.type)
  }

  s.param$psi.type <- get.psi.code(s.param$psi.type)



  results <- .Call("R_threestep_c",pm(object), mm(object), probeNames(object), ngenes, normalize, background, get.background.code(background.method), get.normalization.code(normalize.method), get.summary.code(summary.method),b.param,n.param,s.param, PACKAGE="AffyPLM")

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


