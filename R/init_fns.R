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
## Mar 14, 2004 - added Mbox and MAplot functions for exprSet
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
            function(object,log=FALSE,ref=NULL,subset=NULL,which=NULL,...){
              if(log){
                x <- log2(exprs(object))
              } else {
                x <- exprs(object)
              }
              if (is.null(which)){
                which <-1:dim(x)[2]
              }

              if (is.null(subset)){
                if (is.null(ref)){
                  medianchip <- apply(x, 1, median)
                } else {
                  medianchip <- x[,ref]
                }
                M <- sweep(x,1,medianchip,FUN='-')
                A <- 1/2*sweep(x,1,medianchip,FUN='+')
                if (is.null(ref)){
                  for (i in which){
                    title <- paste(sampleNames(object)[i],"vs pseudo-median reference chip")
                    ma.plot(A[,i],M[,i],main=title,xlab="A",ylab="M",pch='.',...)
                  }
                } else {
                  for (i in which){
                    if (i != ref){
                      title <- paste(sampleNames(object)[i],"vs",sampleNames(object)[ref])
                      ma.plot(A[,i],M[,i],main=title,xlab="A",ylab="M",pch='.',...)
                    }
                  }
                }
              } else {
                if (is.null(ref)){
                  medianchip <- apply(x[,subset], 1, median)
                } else {
                  if (is.element(ref,subset)){
                    medianchip <- x[,ref]
                  } else {
                    stop("Ref ",ref, "is not part of the subset")
                  }
                }
                if (!all(is.element(which,subset))){
                  stop("Specified arrays not part of subset")
                }
                M <- sweep(x,1,medianchip,FUN='-')
                A <- 1/2*sweep(x,1,medianchip,FUN='+')
                if (is.null(ref)){
                  for (i in which){
                    title <- paste(sampleNames(object)[i],"vs pseudo-median reference chip")
                    ma.plot(A[,i],M[,i],main=title,xlab="A",ylab="M",pch='.',...)
                  }
                } else {
                  for (i in which){
                    if (i != ref){
                      title <- paste(sampleNames(object)[i],"vs",sampleNames(object)[ref])
                      ma.plot(A[,i],M[,i],main=title,xlab="A",ylab="M",pch='.',...)
                    }
                  }
                }
                
              }
              
            })
  
}




.initAffyBatchFunctions <- function(where){


##  if( is.null(getGeneric("image.raw")))
##    setGeneric("image.raw",function(object,...)
##               standardGeneric("image.raw"))


###
### This is here because it is not worth fighting with
### affy authors about the correct orientation of these plots
### Note that artful use of axis() would allow tick marks
### to be drawn onto the image in such a way that the
### reversal would not be a problem.
###
  
###  setMethod("image.raw",signature("AffyBatch"),
###            function(object, which=0,transfo=log2, col=gray(c(0:256)/256),xlab="",ylab="", ...){
###             if (which == 0){
###               which <- 1:length(sampleNames(object))
###              }
###              for(i in which){
###                m <- object@exprs[,i]
###                if (is.function(transfo)) {
###                  m <- transfo(m)
###               }
###               m <-  matrix(m, nrow=nrow(object), ncol=ncol(object))
### m <- as.matrix(rev(as.data.frame(m)))
######               image(1:nrow(object), 1:ncol(object), m,
###                     col=col, main=sampleNames(object)[i],
###                     xlab=xlab, ylab=ylab, xaxt="n",yaxt="n", ...)
###             }
###           })
}











.First.lib <- function(libname, pkgname) {
  s <- search() 
  
  require(affy,quietly = FALSE, warn.conflicts = FALSE)
  require(affydata,quietly = FALSE, warn.conflicts = FALSE)
  require(gcrma,quietly = FALSE, warn.conflicts = FALSE)
  
  .initNormfunctions(match(paste("package:",pkgname,sep=""),search()))
  .initExprSetFunctions(match(paste("package:",pkgname,sep=""),search()))
  .initAffyBatchFunctions(match(paste("package:",pkgname,sep=""),search()))
  
  #if (length(search()) > length(s)) {
  #  detach("package:affyPLM")
  #  library(affyPLM,warn.conflicts=FALSE,verbose=FALSE)
  #} 	

  
  library.dynam("affyPLM",pkgname,libname,now=FALSE)
  
  current.normmethods <- get("normalize.AffyBatch.methods",envir=as.environment("package:affy"))
  
  assign("normalize.AffyBatch.methods",
         c(current.normmethods,"quantiles.probeset","scaling"),
         envir=as.environment(match("package:affy", search())))

  # load the Lapack library needed for some parts of fitPLM
  .C("Lapack_Init",PACKAGE="affyPLM")
  
}
