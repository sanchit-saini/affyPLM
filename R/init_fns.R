############################################################
##
## file: init_fns.R
##
## Copyright (C) 2003   Ben Bolstad
##
## aim: implemement initialization functions for AffyExtensions
##
## Created by: B. M. Bolstad <bolstad@stat.berkeley.edu>
##
##
##
## Aug 22, 2003 - Added a initialization function registering
##                the a normalize method for exprSet objects
## Aug 23, 2003 - make sure to make scaling available via normalize
## Sep 11, 2003 - Added a boxplot function for exprSets
## Oct 29, 2003 - Port to R-1.8.0
## Mar 14, 2004 - added Mbox and MAplot functions for exprSets
##
############################################################


.initNormfunctions <- function(where){
  all.affy <- ls(where)
  start <- nchar("normalize.exprSet.")
  assign("normalize.exprSet.methods",
         substr(all.affy[grep("normalize\.exprSet\.*", all.affy)], start+1, 100),
         envir=as.environment(where))
 
  setMethod("normalize", signature(object="exprSet"),
            function(object, method=getOption("BioC")$affy$normalize.method, ...) {
              method <- match.arg(method, normalize.exprSet.methods)
              if (is.na(method))
                stop("unknown method")
              method <- paste("normalize.exprSet", method, sep=".")
              object <- do.call(method, alist(object, ...))
              return(object)
            })
  

}


.initExprSetFunctions <- function(where){

  if (!isGeneric("boxplot"))
    setGeneric("boxplot")
  
  setMethod("boxplot", signature(x="exprSet"),
            function(x,...){
              boxplot(data.frame(exprs(x)),...)
            })

  if (!isGeneric("Mbox"))
    setGeneric("Mbox",function(object,...)
               standardGeneric("Mbox"))
  
  
  setMethod("Mbox",signature("exprSet"),
            function(object,log=FALSE,...){
              if(log){
                x <- log2(exprs(object))
              } else {
                x <- exprs(object)
              }
              medianchip <- apply(x, 1, median)
              M <- sweep(x,1,medianchip,FUN='-')
              boxplot(data.frame(M),...)
            })
 
  if (!isGeneric("MAplot"))
    setGeneric("MAplot",function(object,...)
               standardGeneric("MAplot"))
  
  
  setMethod("MAplot",signature("exprSet"),
            function(object,log=FALSE,ref=NULL,...){
              if(log){
                x <- log2(exprs(object))
              } else {
                x <- exprs(object)
              }
              if (is.null(ref)){
                medianchip <- apply(x, 1, median)
              } else {
                medianchip <- x[,ref]
              }
              M <- sweep(x,1,medianchip,FUN='-')
              A <- 1/2*sweep(x,1,medianchip,FUN='+')
              if (is.null(ref)){
                for (i in 1:dim(x)[2]){
                  title <- paste(sampleNames(object)[i],"vs pseudo-median reference chip")
                  ma.plot(A[,i],M[,i],main=title,xlab="A",ylab="M",pch='.',...)
                }
              } else {
                for (i in (1:dim(x)[2])[-ref]){
                  title <- paste(sampleNames(object)[i],"vs",sampleNames(object)[ref])
                  ma.plot(A[,i],M[,i],main=title,xlab="A",ylab="M",pch='.',...)
                }
              }  
            })
  
}



.First.lib <- function(libname, pkgname) {
  s <- search() 
  
  require(affy,quietly = FALSE, warn.conflicts = FALSE)
  require(affydata,quietly = FALSE, warn.conflicts = FALSE)

  
  .initNormfunctions(match(paste("package:",pkgname,sep=""),search()))
  .initExprSetFunctions(match(paste("package:",pkgname,sep=""),search()))
  
  #if (length(search()) > length(s)) {
  #  detach("package:AffyExtensions")
  #  library(AffyExtensions,warn.conflicts=FALSE,verbose=FALSE)
  #} 	

  
  library.dynam("affyPLM",pkgname,libname,now=FALSE)
  
  current.normmethods <- get("normalize.AffyBatch.methods",envir=as.environment("package:affy"))
  
  assign("normalize.AffyBatch.methods",
         c(current.normmethods,"quantiles.probeset","scaling"),
         envir=as.environment(match("package:affy", search())))

  # load the Lapack library needed for some parts of fitPLM
  .C("Lapack_Init",PACKAGE="affyPLM")
  
}
