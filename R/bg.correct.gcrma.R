bg.correct.gcrma<- function(object,affinity.info=NULL,
                  type=c("fullmodel","affinities","mm","constant"),
                  k=6*fast+0.5*(1-fast),stretch=1.15*fast+1*(1-fast),
                  correction=1,rho=0.7,
                  optical.correct=TRUE,verbose=TRUE,fast=TRUE){

  type <- match.arg(type)


  pmonly <- (type=="affinities"|type=="constant")
  needaff <- (type=="fullmodel"|type=="affinities")
  if( needaff & is.null(affinity.info)){
    if(verbose) cat("Computing affinities")
    affinity.info <- compute.affinities(cdfName(object),verbose=verbose)
    if(verbose) cat("Done.\n")
  }

  if(optical.correct)
    object <- bg.adjust.optical(object,verbose=verbose)

  pm(object) <- gcrma.engine(pms=pm(object),
                             mms=mm(object),
                             pm.affinities=pm(affinity.info),
                             mm.affinities=mm(affinity.info),
                             type=type,k=k,
                             stretch=stretch,
                             correction=correction,rho=rho,
                             verbose=verbose,fast=fast)
  object
}
