###########################################################
##
## file: fitPLM.R
##
## Copyright (C) 2003-2004   Ben Bolstad
##
## created by: B. M. Bolstad <bolstad@stat.berkeley.edu>
## created on: Jan 15, 2003
##
## Last modified: See History below
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
## Jan 18, 2004 - move some of the sub functions to internalfunctions.R
## Feb 23, 2004 - subset paramter added
## Mar 1,  2004 - moved chip-level design matrix creation to its own function
## Mar 2, 2004 - rewrote design matrix function to handle interaction terms
## May 26, 2004 - allow use of MM as response term in fitPLM
## May 27, 2004 - if default model ie -1 + probes + samples flag that a optimized
##                algorithm should be used.
## Jun 28, 2004 - allow MM or PM as a covariate term in the model
##
###########################################################

fitPLM <- function(object,model=PM ~ -1 + probes + samples,variable.type=c(default="factor"),constraint.type=c(default="contr.treatment"),subset=NULL,background=TRUE,normalize=TRUE, background.method = "RMA.2",normalize.method = "quantile",background.param=list(),normalize.param=list(),output.param=list(),model.param=list()){


    if (! is(object, "AffyBatch") {
        stop(paste("argument is",class(object),"fitPLM requires AffyBatch"))
    }

#  oldconstraints <- options()$contrasts
#  if (constraint.type != "contr.treatment"){
#    stop("only endpoint constraint currently implemented")
#  }

#  options(contrasts=c(constraint.type,"contr.poly"))

  PLM.model.matrix <- PLM.designmatrix2(object,model,variable.type,constraint.type)

  #check that chip covariates are not singular

  if (qr(PLM.model.matrix$chip.covariates)$rank < ncol(PLM.model.matrix$chip.covariates)){
    stop("chiplevel effects is singular: singular fits are not implemented in fitPLM")
  }

  if (PLM.model.matrix$response.term == "MM"){
    tmp <- mm(object)
    mm(object) <- pm(object)
    pm(object) <- tmp
  }





  #cat("Method is ",model.type,"\n")
  # now add other variables onto chip.covariates matrix

  rows <- length(probeNames(object,subset))
  cols <- length(object)


  if (is.null(subset)){
    ngenes <- length(geneNames(object))
  } else {
    ngenes <- length(subset)
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


  op.param <- list(weights=TRUE,residuals=TRUE,varcov=c("none","chiplevel","all"),resid.SE=TRUE)
  op.param[names(output.param)] <- output.param

  op.param$varcov <- match.arg(op.param$varcov,c("none","chiplevel","all"))


  md.param <- list(se.type=4,psi.type="Huber",psi.k=NULL,max.its=20,init.method="ls",isdefaultmodel=PLM.model.matrix$isdefaultmodel,MMorPM.covariate=PLM.model.matrix$MMorPM.covariate)
  md.param[names(model.param)] <- model.param




  if (is.null(md.param$psi.k)){
    md.param$psi.k <- get.default.psi.k(md.param$psi.type)
  }

  md.param$psi.type <- get.psi.code(md.param$psi.type)



  # lets do the actual model fitting
  fit.results <- .Call("R_rlmPLMset_c",pm(object,subset),mm(object,subset),probeNames(object,subset),
        ngenes, normalize, background,
        get.background.code(background.method), get.normalization.code(normalize.method),PLM.model.matrix$model.type,PLM.model.matrix$chip.covariates,b.param,n.param,op.param,md.param,PACKAGE="affyPLM")



 #put names on matrices and return finished object

  #chip.param.names <- NULL

#  if (has.chipeffects){
    #chip.param.names <- c(chip.param.names,sampleNames(object))
#  } else {
    chip.param.names <- colnames(PLM.model.matrix$chip.covariates)
#  }

  #options(contrasts=oldconstraints)


  #cache probenames rather than calling rownames(pm(object))
  probenames <- rownames(pm(object,subset))
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
  if (PLM.model.matrix$mt.intercept){
      rownames(fit.results[[6]]) <- rownames(fit.results[[1]])
      rownames(fit.results[[7]]) <- rownames(fit.results[[1]])
       if (!PLM.model.matrix$has.probeeffects){
         fit.results[[2]] <- matrix(0,0,0)
         fit.results[[5]] <- matrix(0,0,0)
       }
  } else {
    if (PLM.model.matrix$has.probeeffects){
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
  x@model.description = list(which.function="fitPLM",preprocessing=list(bg.method=background.method,bg.param=b.param,background=background,norm.method=normalize.method,norm.param=n.param,normalize=normalize),modelsettings =list(model.param=md.param,summary.method=NULL,model=model,constraint.type=constraint.type,variable.type=variable.type),outputsettings=op.param)
  x
}


###
###
### OLD defunct function. Will be completely removed at a later date
###
### this function creates the chip-level part of the design matrix
###
###
### it's return value is the matrix
###

PLM.designmatrix <- function(object,model=PM ~ -1 + probes + samples,variable.type=c(default="factor"),constraint.type=c(default="contr.treatment")){

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
    in.parent.frame <- is.element(mt.termlabels.abbrev,c(ls(parent.frame(2)),ls(parent.frame())))

    #cat(mt.termlabels.abbrev)
    #cat(in.pheno,in.parent.frame,"\n")


    if (sum(in.pheno | in.parent.frame) != length(mt.termlabels.abbrev)){
      stop("Specified parameter does not exist in phenoData or parent frames")
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


  if (!has.chipeffects){
    colnames(chip.covariates) <- chipeffect.names
  } else {
    colnames(chip.covariates) <- chip.param.names
  }


  list(chip.covariates = chip.covariates,model.type=model.type,mt.intercept=mt.intercept,has.probeeffects=has.probeeffects)

}


###
###
### The redesigned function
###
###


PLM.designmatrix2 <- function(object,model=PM ~ -1 + probes + samples,variable.type=c(default="factor"),constraint.type=c(default="contr.treatment")){

  model.terms <- terms(model)
  mt.variables <- attr(model.terms,"variables")

  if (!is.element(as.character(mt.variables[[2]]),c("PM","pm","MM","mm"))){
    stop(paste("Response term in model should be 'PM' or 'pm' or 'MM' or 'mm'."))
  }

  if (is.element(as.character(mt.variables[[2]]),c("PM","pm"))){
    response.term <- "PM"
  } else {
    response.term <- "MM"
  }





  mt.intercept <- attr(model.terms,"intercept")
  mt.termlabels <- attr(model.terms,"term.labels")
  mt.order <- attr(model.terms,"order")

  length.parameters <- length(mt.termlabels)

  if (length.parameters < 1){
    stop("Insufficent parameters supplied in model")
  }


  # checking for probes and samples as model variables
  # Also check for  MM
  # check to see if there is other chip level parameters and that they are valid
  mt.termlabels.abbrev <- mt.termlabels[mt.termlabels != "samples"]
  mt.termlabels.abbrev <- mt.termlabels.abbrev[mt.termlabels.abbrev != "probes"]
  mt.order.abbrev <- mt.order[mt.termlabels != "samples" & mt.termlabels != "probes"]

  has.probeeffects <- is.element("probes",mt.termlabels)
  has.chipeffects <- is.element("samples",mt.termlabels)

  has.covariate.PMorMM <- is.element(c("MM","mm","PM","pm"),mt.termlabels)

  if(sum(has.covariate.PMorMM) > 1){
    stop("Can't have both PM and MM as covariates in model\n")
  }
  if (any(has.covariate.PMorMM)){
    mt.order.abbrev <- mt.order[mt.termlabels.abbrev != "PM" & mt.termlabels.abbrev != "pm" & mt.termlabels.abbrev != "mm" & mt.termlabels.abbrev != "MM" ]
    mt.termlabels.abbrev <- mt.termlabels.abbrev[mt.termlabels.abbrev != "PM" & mt.termlabels.abbrev != "pm" & mt.termlabels.abbrev != "mm" & mt.termlabels.abbrev != "MM" ]
    which <- is.element(c("MM","mm","PM","pm"),mt.termlabels)
    if (c("MM","mm","PM","pm")[which] == response.term){
      stop("Can not have response and covariate PM (or MM).")
    }
  }

  if((!is.null(variable.type$pm)) |  (!is.null(variable.type$mm)) |  (!is.null(variable.type$pm)) | (!is.null(variable.type$pm))){
    cat("Warning: A variable type has been specified for PM or MM variable, this will be ignored. Assuming variable type is covariate\n")
  }

  if (!is.null(variable.type$probes)){
    cat("Warning: A variable type has been specified for probes parameter, this will be ignored. Assuming variable type is factor\n")
  }
  if (!is.null(variable.type$samples)){
    cat("Warning: A variable type has been specified for chips parameter, this will be ignored. Assuming variable type is factor\n")
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


  vt <- rep(vt.default,length(mt.variables)-2)
  names(vt) <- as.character(mt.variables[3:length(mt.variables)])
  vt[names(variable.type)] <- variable.type




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

  nsamples <- dim(exprs(object))[2]

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

    if (has.chipeffects){
      stop("Cannot fit a model with an effect for every chip and additional parameters.\n")
    }

   # print(mt.termlabels.abbrev)
   # print(vt)
   # print(ct)

    if (mt.intercept){
      the.formula <- "~"
      firstfactor <- FALSE
    } else {
      the.formula <- "~ -1 +"
      firstfactor <- TRUE
    }
    the.frame <- NULL
    the.frame.names <- NULL
    terms.names <- NULL

    for (i in 1:length(mt.termlabels.abbrev)){
      if (mt.order.abbrev[i] == 1){
        if (is.element(mt.termlabels.abbrev[i],  names(pData(object)))){
          the.frame <- cbind(the.frame,pData(object)[,names(pData(object))==mt.termlabels.abbrev[i]])
          trt.values <- pData(object)[,names(pData(object))==mt.termlabels.abbrev[i]]
        } else {
          the.frame <- cbind(the.frame,eval(as.name(mt.termlabels.abbrev[i])))
          trt.values <- eval(as.name(mt.termlabels.abbrev[i]))
        }

        if (vt[mt.termlabels.abbrev[i]] == "covariate"){
          the.formula <- paste(the.formula,mt.termlabels.abbrev[i])
          terms.names <- c(terms.names,mt.termlabels.abbrev[i])
        } else {
          if (ct[mt.termlabels.abbrev[i]] == "contr.treatment"){
            the.formula <- paste(the.formula, "C(as.factor(",mt.termlabels.abbrev[i],"),","contr.treatment",")")
            if (!firstfactor){
              for (levs in levels(as.factor(trt.values))[-1]){
                terms.names <- c(terms.names,paste(mt.termlabels.abbrev[i],"_",levs,sep=""))
              }
            } else {
              for (levs in levels(as.factor(trt.values))){
                terms.names <- c(terms.names,paste(mt.termlabels.abbrev[i],"_",levs,sep=""))
              }
              firstfactor <- FALSE
            }
          } else {
            the.formula <- paste(the.formula, "C(as.factor(",mt.termlabels.abbrev[i],"),","contr.sum",")")
            if (!firstfactor){
              for (levs in levels(as.factor(trt.values))[-length(levels(as.factor(trt.values)))]){
                terms.names <- c(terms.names,paste(mt.termlabels.abbrev[i],"_",levs,sep=""))
              }
            } else {
              for (levs in levels(as.factor(trt.values))){
                terms.names <- c(terms.names,paste(mt.termlabels.abbrev[i],"_",levs,sep=""))
              }
              firstfactor <- FALSE
            }
          }
        }


        the.frame.names <- c(the.frame.names,mt.termlabels.abbrev[i])
      } else {
        # figure out which terms are in the higher order model
        in.term <- strsplit(mt.termlabels.abbrev[i],":")

        # first check to see if any of the corresponding lower order terms are going to be in the model
        # if not then need to add data to dataframe.
        # check that dealing with all factor variables


        for (current.term in in.term[[1]]){
          if (vt[current.term] != "factor"){
            stop("Can't have interaction terms involving covariate variables")
          }
        }


        for (current.term in in.term[[1]]){
          if (!(is.element(current.term,mt.termlabels.abbrev))){
            if (is.element(current.term,  names(pData(object)))){
              the.frame <- cbind(the.frame,pData(object)[,names(pData(object))==current.term])
            } else {
              the.frame <- cbind(the.frame,eval(as.name(current.term)))
            }
            the.frame.names <- c(the.frame.names,current.term)
          }
        }

        # now lets make the actual formula

        this.term <- "C("
        for (current.term in in.term[[1]]){
          this.term <- paste(this.term,"as.factor(",current.term,")")
          if (current.term != in.term[[1]][length(in.term[[1]])]){
            this.term <- paste(this.term,":")
          }
        }

        if (ct[mt.termlabels.abbrev[i]] == "contr.treatment"){
          this.term <- paste(this.term,",contr.treatment)")
        } else {
          this.term <- paste(this.term,",contr.sum)")
        }


        ### start off with first factor all levels then
        ### build up from that
        firstterm <- TRUE
        interaction.terms <- as.vector("")
        for (current.term in in.term[[1]]){
          if (is.element(current.term,  names(pData(object)))){
            trt.values <- pData(object)[,names(pData(object))==current.term]
          } else {
            trt.values <- eval(as.name(current.term))
          }

          levs <- levels(as.factor(trt.values))
          if (!firstterm){
            interaction.terms <-  as.vector(t(outer(interaction.terms,paste(":",current.term,"_",levs,sep=""),paste,sep="")))
          } else {
            interaction.terms <-  as.vector(t(outer(interaction.terms,paste(current.term,"_",levs,sep=""),paste,sep="")))
            firstterm <- FALSE
          }
          #print(interaction.terms)

        }
        if (!firstfactor){
          if (ct[mt.termlabels.abbrev[i]] == "contr.treatment"){
            interaction.terms <- interaction.terms[-1]
          } else {
            interaction.terms <- interaction.terms[-length(interaction.terms)]
          }
        } else {
          firstfactor <- FALSE
        }

        terms.names <- c(terms.names,interaction.terms)


        the.formula <- paste(the.formula,"+",this.term)
      }
      if ( i != length(mt.termlabels.abbrev) & mt.order.abbrev[i+1] == 1 ){
        the.formula <- paste(the.formula,"+")
      }
    }
    the.frame <- as.data.frame(the.frame)
    colnames(the.frame) <- the.frame.names  ###mt.termlabels.abbrev[mt.order.abbrev==1]

    #print(the.formula)
    #print(the.frame)
    #print(as.list(ct))
    #print(model.matrix(as.formula(the.formula),the.frame))
    chip.effects <- model.matrix(as.formula(the.formula),the.frame)

    if (mt.intercept){
      chip.effects <- as.matrix(chip.effects[,-1])
    }

  }


  #print(terms.names)


  #print(chip.effects)
  #colnames(chip.effects) <- chipeffect.names
  #if (mt.intercept){
  #  chip.covariates <- chip.effects[,-1]
  #} else {
    chip.covariates <- chip.effects
  #}


  if (!has.chipeffects){
    colnames(chip.covariates) <- terms.names

  } else {
      colnames(chip.covariates) <- chip.param.names
  }

  ## Check to see if it is the default model
  if ((model== (PM ~ -1 + probes + samples)) | (model== (PM ~ -1 + samples+probes))){
    if (ct.probe == "contr.sum"){
      isdefaultmodel <- TRUE
    } else {
      isdefaultmodel <- FALSE
    }
  } else {
    isdefaultmodel <- FALSE
  }



  list(chip.covariates = chip.covariates,model.type=model.type,mt.intercept=mt.intercept,has.probeeffects=has.probeeffects,response.term=response.term,isdefaultmodel=isdefaultmodel, MMorPM.covariate=any(has.covariate.PMorMM))

}
 #PLM.designmatrix2(Dilution,PM ~ -1 + probes+ liver*scanner)
#PLM.designmatrix2(Dilution,PM ~ probes+ liver + scanner + logliver,variable.type=c(logliver="covariate"))

#PLM.designmatrix2(Dilution,PM ~ -1 + samples + probes)
#PLM.designmatrix2(Dilution,PM ~ MM + -1 + samples + probes)
