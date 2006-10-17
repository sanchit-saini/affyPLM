############################################################
##
## file: init_fns.R
##
## Copyright (C) 2003-2006   Ben Bolstad
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
## Sep 13, 2005 - add vignettes to windows menu
## Jun 22, 2006 - add pch to MAplot
## Jul 21, 2006 - MAplot now handles sampleName arguments. Removed the subset argument.
##                now ref and which can basically do the same things. Added pairs as a possible argument.
## Jul 26, 2006 - added hist() method for exprSet
## Oct 11, 2006 - make apply(x,1,median) a rowMedians() call
##
############################################################


.initNormfunctions <- function(where){
  all.affy <- ls(where)
  start <- nchar("normalize.exprSet.")
  assign("normalize.exprSet.methods",
         substr(all.affy[grep("normalize\\.exprSet\\.*", all.affy)], start+1, 100),
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
            function(x,range=0,...){
              boxplot(data.frame(exprs(x)),range=range,...)
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
              medianchip <- rowMedians(x)    ###apply(x, 1, median)
              M <- sweep(x,1,medianchip,FUN='-')
              boxplot(data.frame(M),...)
            })
 
  if (!isGeneric("MAplot"))
    setGeneric("MAplot",function(object,...)
               standardGeneric("MAplot"))
  
  
  setMethod("MAplot",signature("exprSet"),
            function(object,log=FALSE,groups=NULL,ref=NULL,which=NULL,pch=".",ref.fn=c("median","mean"),ref.title="vs pseudo-median reference chip",pairs=FALSE,...){



              if (is.null(groups)){
                if (is.character(ref)){
                  ref.indices <- match(ref,sampleNames(object))
                  if (all(is.na(ref.indices))){
                    stop("No known sampleNames in ref")
                  }
                  
                  if (any(is.na(ref.indices))){
                    warning(paste("Omitting the following from ref:",ref[is.na(ref.indices)], "because they can not be found."))
                  }
                  ref <- ref.indices[!is.na(ref.indices)]
                }
                
                if (is.character(which)){
                  which.indices <- match(which,sampleNames(object))
                  if (all(is.na(which.indices))){
                    stop("No known sampleNames in ref")
                  }
                  
                  if (any(is.na(which.indices))){
                    warning(paste("Omitting the following from which:",which[is.na(which.indices)], "because they can not be found."))
                  }
                  which <- which.indices[!is.na(which.indices)]
                }
                
                
                ref.fn <- match.arg(ref.fn)
                
                if(log){
                  x <- log2(exprs(object))
                } else {
                  x <- exprs(object)
                }
                
                
                if (!pairs){
                  if (is.null(which)){
                    which <-1:dim(x)[2]
                  }
                  
                  if (is.null(ref)){
                    if (ref.fn == "median"){
                      medianchip <- rowMedians(x)    ####apply(x, 1, median)
                    } else {
                      medianchip <- rowMeans(x)
                    }
                  } else if (length(ref) > 1){
                    if (ref.fn == "median"){
                      medianchip <- rowMedians(x[,ref])
                    } else {
                      medianchip <- rowMeans(x[,ref])
                    }
                  }  else {
                    medianchip <- x[,ref]
                  }
                  M <- sweep(x,1,medianchip,FUN='-')
                  A <- 1/2*sweep(x,1,medianchip,FUN='+')
                  if (is.null(ref)){
                    for (i in which){
                      title <- paste(sampleNames(object)[i],ref.title)
                      ma.plot(A[,i],M[,i],main=title,xlab="A",ylab="M",pch=pch,...)
                    }
                  } else {
                    for (i in which){
                      if (length(ref) == 1){
                        if (i != ref){
                          title <- paste(sampleNames(object)[i],"vs",sampleNames(object)[ref])
                          ma.plot(A[,i],M[,i],main=title,xlab="A",ylab="M",pch=pch,...)
                        }
                      } else {
                        title <- paste(sampleNames(object)[i],"vs",ref.title)
                        ma.plot(A[,i],M[,i],main=title,xlab="A",ylab="M",pch=pch,...)
                      }
                    }
                    
                  }
                } else {
                  if (!is.null(ref)) stop("Can't use pairs with non-null 'ref'")
                  if(is.null(which)) which <- 1:ncol(x)
                  mva.pairs(x[,which],log.it=FALSE,...)
                }
              } else {
                ## group labels have been given

                ## check that group variable is of same length as number of samples 
                
                if (dim(exprs(object))[2] != length(groups)){
                  stop("'groups' is of wrong length.")
                }
                

                ### group labels variable can be integer, character or factor variable.
                ### need to check that if any names supplied
                ### for ref or which can be found in group.labels
                
                if (!is.null(which)){
                  if (is.numeric(groups)){
                    if (!is.numeric(which)){
                      stop("'which' labels must also be found in 'groups'") 
                    } else {
                      if (!all(is.element(which,groups))){
                        stop("'which' labels must also be found in 'groups'") 
                      }
                    }
                  } else if (is.factor(groups)){
                    if (!is.character(which)){
                      stop("'which' should be character vector") 
                    } else {
                      if (!all(is.element(which,as.character(groups)))){
                        stop("'which' labels must also be found in 'groups'") 
                      }
                    }
                  } else if (is.character(groups)){
                    if (!is.character(which)){
                      stop("'which' should be character vector") 
                    } else {
                      if (!all(is.element(which,groups))){
                        stop("'which' labels must also be found in 'groups'") 
                      }
                    }
                  }
                }

                if (!is.null(ref)){
                  if (is.numeric(groups)){
                    if (!is.numeric(ref)){
                      stop("'ref' labels must also be found in 'groups'") 
                    } else {
                      if (!all(is.element(ref,groups))){
                        stop("'ref' labels must also be found in 'groups'") 
                      }
                    }
                  } else if (is.factor(groups)){
                    if (!is.character(ref)){
                      stop("'ref' should be character vector") 
                    } else {
                      if (!all(is.element(ref,as.character(groups)))){
                        stop("'ref' labels must also be found in 'groups'") 
                      }
                    }
                  } else if (is.character(groups)){
                    if (!is.character(ref)){
                      stop("'ref' should be character vector") 
                    } else {
                      if (!all(is.element(ref,groups))){
                        stop("'ref' labels must also be found in 'groups'") 
                      }
                    }
                  }
                }
                
                ref.fn <- match.arg(ref.fn)
                
                if(log){
                  x <- log2(exprs(object))
                } else {
                  x <- exprs(object)
                }
                groups.list <- split(1:dim(exprs(object))[2], as.factor(groups))


                grouped.data <- matrix(0,nrow(x),length(groups.list))
                colnames(grouped.data) <- names(groups.list)
                which.col <- 1
                for (group in groups.list){
                  grouped.data[,which.col] <- rowMeans(x[,group,drop=FALSE])
                  which.col <- which.col + 1
                }

                
                if (!pairs){
                  if (is.null(which)){
                    which <- names(groups.list)
                  }
                  
                  if (is.null(ref)){
                    if (ref.fn == "median"){
                      medianchip <- rowMedians(grouped.data)   ####apply(grouped.data, 1, median)
                    } else {
                      medianchip <- rowMeans(grouped.data)
                    }
                    
                  } else if (length(ref) == 1){
                    ref.name <- ref
                    ref <- match(ref,names(groups.list))
                    medianchip <- grouped.data[,ref]
                  } else {
                    ref <- match(ref,names(groups.list))
                    if (ref.fn == "median"){
                      medianchip <- rowMedians(grouped.data[,ref])
                    } else {
                      medianchip <- rowMeans(grouped.data[,ref])
                    }

                  }

                  M <- sweep(grouped.data,1,medianchip,FUN='-')
                  A <- 1/2*sweep(grouped.data,1,medianchip,FUN='+')
                  if (is.null(ref)){
                    for (i in which){
                      title <- paste(i,ref.title)
                      ma.plot(A[,i],M[,i],main=title,xlab="A",ylab="M",pch=pch,...)
                    }
                  } else {
                    for (i in which){
                      if (length(ref) == 1){
                        if (i != ref.name){
                          title <- paste(i,"vs",ref.name)
                          ma.plot(A[,i],M[,i],main=title,xlab="A",ylab="M",pch=pch,...)
                        }
                      } else {
                        title <- paste(i,ref.title)
                        ma.plot(A[,i],M[,i],main=title,xlab="A",ylab="M",pch=pch,...)
                      }
                    }
                    
                  }

                  

                } else {
                  if (!is.null(ref)) stop("Can't use pairs with non-null 'ref'")
                  if (is.null(which)){
                    which <- names(groups.list)
                  }
                  
                  mva.pairs(grouped.data[,which],log.it=FALSE,...)


                }
              }
            })


  plotDensity.exprSet <- function(x, col=1:6, log=FALSE,
                                  ylab="density",
                                  xlab=NULL,
                                  ...){
  
    x <- exprs(x)
    
    if(log){
      x <- log2(x)
      if(is.null(xlab)) xlab <- "log intensity"
    }
    else  if(is.null(xlab)) xlab <- "intensity"
    
    rv <- plotDensity(x, ylab=ylab, xlab=xlab, col=col, ...)
    
    invisible(rv)
  }

  

  if(!isGeneric("hist"))
    setGeneric("hist")
  
  setMethod("hist",signature(x="exprSet"), 
            function(x,...) plotDensity.exprSet(x,...))
  

  
  
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

  if(interactive() && .Platform$OS.type == "windows" &&
                  .Platform$GUI == "Rgui"){
             addVigs2WinMenu("affyPLM")
           }

  
}
