###
###
### file: preprocess.R
###
### aim: implement background and normalization routines as a single function 
###
###


preprocess <- function(object,subset=NULL, normalize=TRUE,background=TRUE,background.method="RMA.2",normalize.method="quantile",background.param = list(),normalize.param=list(),verbosity.level=0){


  if (!is(object, "AffyBatch")) {
    stop(paste("argument is",class(object),"threestep requires AffyBatch"))
  }

  if (is.null(subset)){
    ngenes <- length(geneNames(object))
  } else {
    ngenes <- length(unique(subset))
  }
  

  ### it looks a little weird to be setting up a model here but it is because
  ### we are reusing some functions from the fitPLM pipeline 
  ### we also have to use a little trickery here to get the normalization
  ### method to handle together or separately.
  
  R.model <- PLM.designmatrix3(object,PM ~ -1 + samples,variable.type=c(default="factor"),constraint.type=c(default="contr.treatment"))


  background.param <- verify.bg.param(R.model, background.method,background.param = background.param)

  if (!is.null(normalize.param$type)){
    if ((normalize.param$type == "together") || (normalize.param$type == "separate")){
      R.model <- PLM.designmatrix3(object,PMMM ~ -1 + samples,variable.type=c(default="factor"),constraint.type=c(default="contr.treatment"))
    } else {
      R.model <- PLM.designmatrix3(object,PM ~ -1 + samples,variable.type=c(default="factor"),constraint.type=c(default="contr.treatment"))
    }
  }

  normalize.param <- verify.norm.param(R.model, normalize.method,normalize.param=normalize.param)
  
  
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
  
  print(background.param$type)
  print(normalize.param$type)

  if ((background.param$type == "pmonly") && (normalize.param$type == "pmonly")){
### SEXP pp_bothstages(SEXP PMmat, SEXP MMmat, SEXP ProbeNamesVec,SEXP N_probes,SEXP norm_flag, SEXP bg_flag, SEXP bg_type,SEXP norm_type, SEXP background_parameters,SEXP norm_parameters, SEXP verbosity)
    pm(object) <- .Call("pp_bothstages",pm(object,subset),
                        mm(object,subset),
                        probeNames(object,subset),
                        ngenes, normalize, background,
                        background.method, normalize.method,
                        background.param,normalize.param,
                        verbosity.level, PACKAGE="affyPLM")
  } else {
    x <- mm(object,subset)
    pm(object) <- .Call("pp_bothstages",pm(object,subset),
                        x,
                        probeNames(object,subset),
                        ngenes, normalize, background,
                        background.method, normalize.method,
                        background.param,normalize.param,
                        verbosity.level, PACKAGE="affyPLM")
    mm(object) <- x
  }



  
  object


}
