###########################################################
##
## file: fitPLM.R
##
## Copyright (C) 2003   Ben Bolstad
##
## created by: B. M. Bolstad <bolstad@stat.berkeley.edu>
## created on: Jan 15, 2003
##
## Last modified: Feb 22, 2003
##
## method for fitting robust linear models.
## interface to C code which does the heavy lifting
##
##
## by default model should be written in the following terms
##
## PM  : should appear on left hand side
## ~ : equivalent to  =
## -1 : remove the intercept term
## probes : include terms for probes in the model. a variable of type factor
## samples : defaults to be a factor, one for each chip or array
## other names for chip level factors or covariates. By default we will check the
## phenoData object to see if the named covariate or factor exists in the
## phenoData slot of the provided affybatch, then the parent enivronment is checked
## if the named variable does not exist there then procedure exits.
##
## an example model fitting terms for probes and samples
##
##    pm ~ -1 + probes + samples
##
## another example model
##
##    pm ~ -1 + probes + trt.group
##
## where trt.group is a factor in the phenoData object.
##
## A valid model should always include some chip effects
##
##
## Feb 1, 2003 - cache rownames(PM(object)) since this is slow operation
## Feb 15 - set model description to the function call: THIS was undone.
##         made it possible to fit models with an intercept term.
## Feb 17 -  start testing and fixing up column naming when model fit with no
##         intercept term, also make sure that se.exprs is set if there is
##         an intercept in the fitted model. re-order some parts of the function
##         work towards being able to fit models with no slope or probe-effects.
## Feb 22 - Add ability to fit models where parameters use a sum to zero constraint.
##          A similar mechanism to the variables.type parameter will be used.
##          note that R uses the convention for endpoint constraint that first
##          item is constrained to 0 and contr.sum that the last item is left unshown
##          we will try to follow this convention, except in the case constraints on
##          the probe coefficient.
##          Remember that the probe coefficient is handled by c code.
## Mar 22 - Add in ability to use LESN background correction.
## Mar 23 - Add ability to use MAS style background in fit PLM
## Jun 4-8 - Ability to specify other psi functions.
## Jul 22 - No longer need to specify psi.k, if left unspecified (ie null)
##          we set it to the default parameter. These defaults are the following:
## ** Huber - k = 1.345
## ** Fair - k = 1.3998
## ** Cauchy - k=2.3849
## ** Welsch - k = 2.9846
## ** Tukey Biweight - k = 4.6851
## ** Andrews Sine - K = 1.339
##          the other methods have no tuning parameter so set to 1.
## Jul 27 - Clean up parameters passed to function
## Sep 2  - Residuals are now stored in the outputed PLMset
## Sep 4  - may toggle what is actually outputed
##          in particular in respect to weights, residuals
##          var.cov and resid.SE
## Sep6 - Sept 8  -  More changes to how naming is done.
## Sep 12 - Fix how constraints are applied when factor variable
##          is defined in parent enivronment. Basically the constraint
##          was being ignored. Note that this fix was also applied to
##          the 0.5-14 Series package.
## Sep 12 - model.param was introduced as argument
##          se.type, psi.type, psi.k were placed into
##          this parameter
## Oct 10 - fix spelling of McClure
##
###########################################################

fitPLM <- function(object,model=PM ~ -1 + probes + samples,variable.type=c(default="factor"),constraint.type=c(default="contr.treatment"),background=TRUE,normalize=TRUE, background.method = "RMA.2",normalize.method = "quantile",background.param=list(),normalize.param=list(),output.param=list(),model.param=list()){


  get.background.code <- function(name) {
    background.names <- c("RMA.1", "RMA.2", "IdealMM","MAS","MASIM","LESN2","LESN1","LESN0")
    if (!is.element(name, background.names)) {
      stop(paste(name, "is not a valid summary method. Please use one of:",
                 "RMA.1", "RMA.2", "IdealMM","LESN2","LESN1","LESN0"))
    }
    code <- c(1, 2, 3, 4, 5, 6, 7, 8)[name == background.names]
    code
  }

  get.normalization.code <- function(name) {
    normalization.names <- c("quantile","quantile.probeset","scaling")
    if (!is.element(name, normalization.names)) {
      stop(paste(name, "is not a valid summary method. Please use one of:",
                 "quantile","quantile.probeset","scaling"))
    }
    code <- c(1,2,3)[name == normalization.names]
    code
  }

  get.psi.code <- function(name){
    psi.names <- c("Huber","fair","Cauchy","Geman-McClure","Welsch","Tukey","Andrews")
     if (!is.element(name, psi.names)) {
      stop(paste(name, "is not a valid Psi type. Please use one of:",
                 "Huber","fair","Cauchy","Geman-Mclure","Welsch","Tukey","Andrews"))
    }
    code <- c(0:6)[name == psi.names]
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






  if (class(object) != "AffyBatch"){
    stop(paste("argument is",class(object),"fitPLM requires AffyBatch"))
  }

#  oldconstraints <- options()$contrasts
#  if (constraint.type != "contr.treatment"){
#    stop("only endpoint constraint currently implemented")
#  }

#  options(contrasts=c(constraint.type,"contr.poly"))


  model.terms <- terms(model)
  mt.variables <- attr(model.terms,"variables")

  if ((mt.variables[[2]] != "PM") & (mt.variables[[2]] != "pm")){
    stop(paste("Response term in model should be 'PM' or 'pm'"))
  }

  mt.intercept <- attr(model.terms,"intercept")
  mt.termlabels <- attr(model.terms,"term.labels")


  length.parameters <- length(mt.termlabels)

  if (length.parameters < 1){
    stop("Insufficent parameters supplied in model")
  }

  # check to see if there is other chip level parameters and that they are valid
  mt.termlabels.abbrev <- mt.termlabels[mt.termlabels != "samples"]
  mt.termlabels.abbrev <- mt.termlabels.abbrev[mt.termlabels.abbrev != "probes"]

  has.probeeffects <- is.element("probes",mt.termlabels)
  has.chipeffects <- is.element("samples",mt.termlabels)




  #if (!has.probeeffects){
  #  stop("Fitting a model without probe-effects currently not supported")
  #}

  #check to see if anyone dumb enough to specify variable types for the "probes" and "samples" parameters

  if (!is.null(variable.type$probes)){
    cat("A variable type has been specified for probes parameter, this will be ignored. Assuming variable type is factor\n")
  }
  if (!is.null(variable.type$samples)){
    cat("A variable type has been specified for chips parameter, this will be ignored. Assuming variable type is factor\n")
  }

  #now go through the non default parameters seeing which type they should be treated as.

  #cat(mt.termlabels.abbrev,"\n")

  vt <- NULL

  #check to see that everything is either "factor" or "covariate"
  #figure out what the default type is
  if (is.na(variable.type["default"])){
    cat("No default type given. Assuming default variable type is factor\n")
    vt.default <- "factor"
  } else {
    if (is.element(variable.type["default"],c("factor","covariate"))){
      vt.default <- variable.type["default"]
    } else {
      stop("Incorrect default variable.type given")
    }
  }

  if (sum(!is.element(variable.type,c("factor","covariate")))){
    stop("An incorrect variable type provided")
  }


  if (length(mt.termlabels.abbrev) >=1){

    vt <- rep(vt.default,length(mt.termlabels.abbrev))
    names(vt) <- mt.termlabels.abbrev

    vt[names(variable.type)] <- variable.type
    #has.vt <- is.element(mt.termlabels.abbrev,names(variable.type))
  }
  #cat(has.vt,"\n");
  #print(vt)

  #figure out what the default constraint type is

  if (is.na(constraint.type["default"])){
    cat("No default type given. Assuming default constraint type is contr.treatment\n")
    ct.default <- "contr.treatment"
  } else {
    if (is.element(constraint.type["default"],c("contr.sum","contr.treatment"))){
      ct.default <- constraint.type["default"]
    } else {
      stop("Incorrect default constraint.type given")
    }
  }

  # check to see if there are any constraints on for the "probes" and "samples" parameters. Constraints on the first chip
  # level parameter will be ignored unless, the default constraint on probes will always "contr.sum", however
  # will allow the user to shoot themselves in the foot by setting it to something else (ie "contr.treatment") and there may be times when
  # this is useful.


  if (is.na(constraint.type["probes"])){
    ct.probe <- "contr.sum"
  } else {
    ct.probe <- constraint.type["probes"]
 }

  if (is.na(constraint.type["samples"])){
    ct.samples <- ct.default
  } else {
    ct.samples <- constraint.type["samples"]
  }


  # Constraint.types given for each variable.
  # note that if user puts a constraint on a covariate variable it will be ignored.

  if (length(mt.termlabels.abbrev) >=1){

    ct <- rep(ct.default,length(mt.termlabels.abbrev))
    names(ct) <- mt.termlabels.abbrev

    ct[names(constraint.type)] <- constraint.type
  #  print(ct)
  }
  #print(ct.probe)
  #print(ct.samples)

  if (mt.intercept == 0){
    if (has.probeeffects){
      if (ct.probe == "contr.sum"){
        model.type <- 0
      } else {
        model.type <- 10
      }
    } else {
      model.type <- 20
    }
  } else {
    #stop("models with intercept terms currently not supported")
    if (has.probeeffects){
      if (ct.probe == "contr.sum"){
        model.type <- 1
      } else {
        model.type <- 11
      }
    } else {
      model.type <- 21
    }
  }


  if (!has.chipeffects & (length(mt.termlabels.abbrev) < 1)){
    stop("Model does not have enough chip level parameters")
  }

  nsamples <- dim(pm(object))[2]

  # start making the chip.covariate matrix
  # when no intercept, we need no contraint on the chipeffects parameter
  #print var.type
  #print(nsamples)
  chip.param.names <- NULL


  if (has.chipeffects){
    our.samples <- 1:nsamples
    if (!mt.intercept){
      chip.effects <- model.matrix(~ -1 + as.factor(our.samples))[,1:nsamples]
      chip.param.names <- sampleNames(object)
    } else {
      chip.effects <- model.matrix( ~ C(as.factor(our.samples),ct.samples))[,2:nsamples]
      if (ct.samples == "contr.treatment"){
        chip.param.names <- sampleNames(object)[2:nsamples]
      } else {
        chip.param.names <- sampleNames(object)[1:(nsamples-1)]
      }
    }
  } else {
    chip.effects <- NULL
  }



  if (length(mt.termlabels.abbrev) >=1){

    in.pheno <- is.element(mt.termlabels.abbrev, names(pData(object)))
    in.parent.frame <- is.element(mt.termlabels.abbrev,ls(parent.frame()))

    #cat(mt.termlabels.abbrev)
    #cat(in.pheno,in.parent.frame,"\n")


    if (sum(in.pheno | in.parent.frame) != length(mt.termlabels.abbrev)){
      stop("Specified parameter does not exist in phenoData or parent frame")
    }

    chipeffect.names <-  NULL

    if (has.chipeffects){
      # chip.effects will handle intercept, use constraints on treatments
      # if chip.effects then the matrix will be singular, but we will go
      # through with the forming the chip.covariates matrix anyway.
      #options(contrasts=c("contr.treatment","contr.poly"))

      stop("Can not fit a model with an effect for every chip and additional parameters")

      #for (trt in mt.termlabels.abbrev){
      #   if (is.element(trt,  names(pData(object)))){
      #     trt.values <- pData(object)[,names(pData(object))==trt]
      #     if (vt[names(vt)==trt] == "covariate"){
      #       trt.effect <- model.matrix(~ -1 + trt.values)
      #     } else {
      #       trt.effect <- model.matrix(~ as.factor(trt.values))[,-1]
      #     }
      #   } else {
      #     trt.values <- eval(as.name(trt))
      #     if (length(trt.values) != nsamples){
      #       stop("Model parameter in parent environment is incorrect")
      #     }
      #     trt.effect <- model.matrix(~ as.factor(trt.values))[,-1]
      #   }
      #   chip.effects <- cbind(chip.effects,trt.effect)
      #   if (vt[names(vt)==trt] == "covariate"){
      #     chipeffect.names <- c(chipeffect.names,trt)
      #   } else {
      #     for (levs in levels(as.factor(trt.values))[-1]){
      #       chipeffect.names <- c(chipeffect.names,paste(trt,"_",levs,sep=""))
      #     }
      #   }
      # }
    } else {
      # no chipeffect, first factor treatment will be unconstrained if no intercept
      #
      #

      #options(contrasts=c("contr.treatment","contr.poly"))
      #print(ct)
      first.factor <- FALSE
      for (trt in mt.termlabels.abbrev){
                                        #     print(vt[names(vt)==trt])
        if (is.element(trt,  names(pData(object)))){
          trt.values <- pData(object)[,names(pData(object))==trt]
          if (vt[names(vt)==trt] == "covariate"){
                                        #         print("is covariate")
            trt.effect <- model.matrix(~ -1 + trt.values)
          } else {
            if (!first.factor){
                                        #first.factor <- TRUE
              if (!mt.intercept){
                trt.effect <- model.matrix(~ -1+ as.factor(trt.values))
              } else {
                trt.effect <- model.matrix(~ C(as.factor(trt.values),ct[trt]))[,-1]
              }
             } else {

               trt.effect <- model.matrix(~  C(as.factor(trt.values),ct[trt]))[,-1]
             }
          }
        } else {
          trt.values <- eval(as.name(trt))
          if (length(trt.values) != nsamples){
            stop("Model parameter in parent environment is incorrect")
          }
          if (vt[names(vt)==trt] == "covariate"){
            trt.effect <- model.matrix(~ -1 + trt.values)
          } else {
            if (!first.factor){
               #first.factor <- TRUE
              if (!mt.intercept){
                trt.effect <- model.matrix(~ -1+ as.factor(trt.values))
              } else {
                trt.effect <- model.matrix(~ C(as.factor(trt.values),ct[trt]))[,-1]
              }
            } else {
              trt.effect <- model.matrix(~ C(as.factor(trt.values),ct[trt]))[,-1]
            }
          }
        }
        chip.effects <- cbind(chip.effects,trt.effect)

        if (vt[names(vt)==trt] == "covariate"){
          chipeffect.names <- c(chipeffect.names,trt)
        } else {
          if (!first.factor){
            first.factor <- TRUE
            if (!mt.intercept){
              for (levs in levels(as.factor(trt.values))){
                chipeffect.names <- c(chipeffect.names,paste(trt,"_",levs,sep=""))
              }
            } else {
              if (ct[trt] == "contr.treatment"){
                for (levs in levels(as.factor(trt.values))[-1]){
                  chipeffect.names <- c(chipeffect.names,paste(trt,"_",levs,sep=""))
                }
              } else {
                for (levs in levels(as.factor(trt.values))[-length(levels(as.factor(trt.values)))]){
                  chipeffect.names <- c(chipeffect.names,paste(trt,"_",levs,sep=""))
                }
              }
            }
          } else {
            if (ct[trt] == "contr.treatment"){
                for (levs in levels(as.factor(trt.values))[-1]){
                  chipeffect.names <- c(chipeffect.names,paste(trt,"_",levs,sep=""))
                }
              } else {
                for (levs in levels(as.factor(trt.values))[-length(levels(as.factor(trt.values)))]){
                  chipeffect.names <- c(chipeffect.names,paste(trt,"_",levs,sep=""))
                }
            }
          }
        }
      }
    }
  }


  #print(chip.effects)
  #colnames(chip.effects) <- chipeffect.names
  chip.covariates <- chip.effects


  #print(chipeffect.names)
  #print(chip.effects)
  #print(colnames(chip.covariates))


  #check that chip covariates are not singular

  if (qr(chip.covariates)$rank < ncol(chip.covariates)){
    stop("chiplevel effects is singular: singular fits are not implemented in fitPLM")
  }

  #cat("Method is ",model.type,"\n")
  # now add other variables onto chip.covariates matrix

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


  op.param <- list(weights=TRUE,residuals=TRUE,varcov=c("none","chiplevel","all"),resid.SE=TRUE)
  op.param[names(output.param)] <- output.param

  op.param$varcov <- match.arg(op.param$varcov,c("none","chiplevel","all"))


  md.param <- list(se.type=4,psi.type="Huber",psi.k=NULL,max.its=20,init.method="ls")
  md.param[names(model.param)] <- model.param




  if (is.null(md.param$psi.k)){
    md.param$psi.k <- get.default.psi.k(md.param$psi.type)
  }

  md.param$psi.type <- get.psi.code(md.param$psi.type)



  # lets do the actual model fitting
  fit.results <- .Call("R_rlmPLMset_c", pm(object), mm(object),
                       probeNames(object), ngenes, normalize, background,
                       get.background.code(background.method),
                       get.normalization.code(normalize.method),
                       model.type,chip.covariates, b.param, n.param,
                       op.param, md.param, PACKAGE="affyPLM")



 #put names on matrices and return finished object

  #chip.param.names <- NULL

  if (has.chipeffects){
    #chip.param.names <- c(chip.param.names,sampleNames(object))
  } else {
    chip.param.names <- c(chip.param.names,chipeffect.names)
  }

  #options(contrasts=oldconstraints)


  #cache probenames rather than calling rownames(pm(object))
  probenames <- rownames(pm(object))
  #cat(chip.param.names,"\n")
  #colnames(fit.results[[1]]) <- sampleNames(object)
  colnames(fit.results[[1]]) <- chip.param.names
  if (op.param$weights){
    colnames(fit.results[[3]]) <- sampleNames(object)
    rownames(fit.results[[3]]) <- probenames
  }

  if (op.param$residuals){
   colnames(fit.results[[8]]) <- sampleNames(object)
   rownames(fit.results[[8]]) <- probenames
  }

  if (op.param$resid.SE){
     colnames(fit.results[[9]]) <- c("residSE","df")
     rownames(fit.results[[9]]) <- rownames(fit.results[[1]])
  }


  if (op.param$varcov != "none"){
    names(fit.results[[10]]) <- rownames(fit.results[[1]])

    if (op.param$varcov == "chiplevel"){
      tmp.colnames <- colnames(fit.results[[1]])

      name.matrix <- function(x,names){
        colnames(x) <- names
        rownames(x) <- names
        x
      }
      fit.results[[10]] <- lapply(fit.results[[10]],name.matrix,tmp.colnames)
    }
  }



  colnames(fit.results[[2]]) <- "ProbeEffects"
  rownames(fit.results[[2]]) <- probenames
  colnames(fit.results[[5]]) <- "SEProbeEffects"
  rownames(fit.results[[5]]) <- probenames
  colnames(fit.results[[4]]) <- chip.param.names
  if (mt.intercept){
      rownames(fit.results[[6]]) <- rownames(fit.results[[1]])
      rownames(fit.results[[7]]) <- rownames(fit.results[[1]])
       if (!has.probeeffects){
         fit.results[[2]] <- matrix(0,0,0)
         fit.results[[5]] <- matrix(0,0,0)
       }
  } else {
    if (has.probeeffects){
      fit.results[[6]] <- matrix(0,0,0)
      fit.results[[7]] <- matrix(0,0,0)
    } else {
      fit.results[[2]] <- matrix(0,0,0)
      fit.results[[5]] <- matrix(0,0,0)
      fit.results[[6]] <- matrix(0,0,0)
      fit.results[[7]] <- matrix(0,0,0)
    }
  }
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
  x@model.description = list(which.function="rmaPLM",preprocessing=list(bg.method=background.method,bg.param=b.param,background=background,norm.method=normalize.method,norm.param=n.param,normalize=normalize),modelsettings =list(model.param=md.param,summary.method=NULL,model=model,constraint.type=constraint.type,variable.type=variable.type),outputsettings=op.param)
  x
}

