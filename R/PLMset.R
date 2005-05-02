###########################################################
##
## file: PLMset.R
##
## Copyright (C) 2003    Ben Bolstad
##
## created by: B. M. Bolstad <bolstad@stat.berkeley.edu>
## created on: Jan 14, 2003
##
##
## aim: define and implement the PLMset object class and
##      its methods
##
## The PLMset object should hold Probe Level Model fits
## in particular we will concentrate on fits by the
## robust linear model methodology.
##
## Will use some of the ideas from the exprSet class
## from Biobase.
##
## the PLMset object has slots to hold probe and chip coefficients
## their standard errors, and model weights, along with the
## usual phenoData, description, annotation and notes fields.
## the exprs and se slot of the parent exprSet will be used for
## storing the constant coef and its se. ie if you fit the model
##
## pm'_ij(k) = mu_(k) + probe_i(k) + chip_j(k) + chipcovariates_j(k) + \epsilon_ij
##
##  then mu(k) would be stored in the exprs slot, its standard error in the se slot
##  probe_i(k) is stored in the probe.coef slot, with its standard error in its respective slot
##  chip_j(k) and chipcovariates_j(k) would be stored in chip.coefs (and ses in se.chip.coefs)
## 
##
## Modification History
##
## Jan 14, 2003 - Initial version, weights, coef accessors
## Jan 16, 2003 - added slots for probe coefs and there standard errors. some people
##                might find these useful. Made PLMset extend exprSet. This
##                saves us having to redefine accessors to phenoData,
##                description, annotataion, notes.
## Feb 15, 2003 - Add in a "show" method. Add in an coefs.probe accessor function
## Apr 29, 2003 - Add a replacement function for se
## Sep 2, 2003 - added some new data members to the object.
##               in particular
##               residuals  - a matrix for storing residuals
##               residualSE - two column matrix residual SE and df
##               normVec - a vector that can be used to establish
##                  quantile normalization for data added at a
##                  later date.
##               model.description is now a list
##               accessors/replacement functions for above
##               image() now has options for display of residuals
## Sep 3, 2003 - image() will now draw a legend if requested
## Sep 5, 2003 - Variance Covariance matrix stored as list is added
##               as data member from object
## Sep 8, 2003 - accessor for resisualsSE and varcov.
##               made image check that weights or residual matrices exist.
## Sep 14, 2003 - fix up which parameter when PLMSet does not have weights
## Oct 10, 2003 - fix labeling on image when use.log =TRUE
## Oct 29, 2003 - port to R-1.8.0 (including some cross-porting from affyPLM in BioC1.3)
## Dec 8, 2003  - replace method for residuals
## Dec 9, 2003  - an indexing function to allow one to pull out appropriate
##                items from the weights, residuals (the accessor functions
##                have been modified to allow a genenames argument)
## Dec 10, 2003 - Residuals can now be given in standardized form
##                Summary function (simplistic)
## Dec 12, 2003 - model.description accessor
##                document the structure of the model description list
##                (see below for a description
##                of the list structure)
## Dec 14, 2003 - Adjust "show" to handle model.description
## Mar 14, 2004 - Added MAplot generic function
## June 23, 2004 - boxplot has type argument. Also the NUSE procedure attempts
##                 to construct a reasonable boxplot even if the default model
##                 has not been used.
## Aug 2, 2004 -  start making changes to the PLMset object.
##                in particular:
##                probe.coefs/se.probe.coefs are now lists
##                weights - list
##                residuals - list
## Feb 18, 2005 - remove ma.plot (it is now in affy)
## Mar 12, 2005 - Mbox() now includes a range arguement (with default value=0)
##                NUSE boxplot is also this way.
##                Added NUSE function (which gives either the NUSE boxplot or the values of NUSE)
##                Added RLE function
## Mar 14, 2005 - NUSE() and RLE() boxplots have default y-limits: ylim=c(0.92, 1.3) for NUSE, ylim=c(-0.75,0.75) for RLE
##                NUSE,RLE plots have horizontal lines
##                Speed up image() in certain situations
## Apr 6, 2005  - ability to change color maps on image()
##
###########################################################


  #creating the PLMset object

setClass("PLMset",
           representation(probe.coefs="list",
                          se.probe.coefs="list",
                          chip.coefs = "matrix",
                          se.chip.coefs = "matrix",
                          cdfName="character",
                          nrow="numeric",
                          ncol="numeric",
                          model.description="list",
                          model.call = "call",
                          weights="list",
                          residuals="list",
                          residualSE="matrix",
                          normVec="matrix", varcov="list"),
                         # phenoData="phenoData",
                         # description="characterORMIAME",
                         # annotation="character",
                         # notes="character"
           prototype=list(
             probe.coefs=list(),                           #matrix(nr=0,nc=0),
             se.probe.coefs=list(),                        #matrix(nr=0,nc=0),
             chip.coefs=matrix(nr=0,nc=0),
             se.chip.coefs=matrix(nr=0,nc=0),
             model.description=list(),
             weights=list(),                               #matrix(nr=0,nc=0),
             residuals =list(),                            #matrix(nr=0,nc=0),
             residualSE=matrix(nr=0,nc=0),
             normVec=matrix(nr=0,nc=0),
             varcov=list(),
             description=new("MIAME"),
             model.description=list(),
             annotation="",
             cdfName="",
             nrow=0, ncol=0,
             notes=""),contains="exprSet")


  #now some accessors.
  
if (is.null(getGeneric("cdfName")))
  setGeneric("cdfName", function(object)
             standardGeneric("cdfName"))

setMethod("cdfName", "PLMset", function(object)
          object@cdfName)


                                        #access weights
setMethod("weights",signature(object="PLMset"),
          function(object,genenames=NULL){ 
		if (is.null(genenames)){
			object@weights
		} else{
		 which <-indexProbesProcessed(object)[genenames]
		 which <- do.call("c",which)
                 if (object@model.description$R.model$response.variable == 0){
                   list(PM.weights=object@weights[[1]][which,],MM.weights=object@weights[[2]][which,])
                 } else if (object@model.description$R.model$response.variable == -1){
                   list(PM.weights=matrix(0,0,0),MM.weights=object@weights[[2]][which,])
                 } else if (object@model.description$R.model$response.variable == 1){
                   list(PM.weights=object@weights[[1]][which,],MM.weights=matrix(0,0,0))
                 }
		}	
	})



if (!isGeneric("weights<-"))
  setGeneric("weights<-",function(object,value)
             standardGeneric("weights<-"))


  #replace weights
setReplaceMethod("weights",signature(object="PLMset"),
                 function(object,value){
                   object@weights <- value
                   object
                 })



  #access parameter estimates (chip level coefficients)

if (!isGeneric("coefs"))
  setGeneric("coefs",function(object)
             standardGeneric("coefs"))
  
setMethod("coefs",signature(object="PLMset"),
            function(object) object@chip.coefs)

if (!isGeneric("coefs<-"))
  setGeneric("coefs<-",function(object,value)
             standardGeneric("coefs<-"))


                                        #replace coefs (chip level coefficients)
setReplaceMethod("coefs",signature(object="PLMset"),
                 function(object,value){
                   object@chip.coefs <- value
                   object
                 })


  #access the probe level coefficents
if (!isGeneric("coefs.probe"))
  setGeneric("coefs.probe",function(object)
             standardGeneric("coefs.probe"))

setMethod("coefs.probe",signature(object="PLMset"),
          function(object) object@probe.coefs)

if (!isGeneric("se"))
  setGeneric("se",function(object)
             standardGeneric("se"))
  
setMethod("se",signature(object="PLMset"),
          function(object) object@se.chip.coefs)

if (!isGeneric("se.probe"))
  setGeneric("se.probe",function(object)
             standardGeneric("se.probe"))
  
setMethod("se.probe",signature(object="PLMset"),
          function(object) object@se.probe.coefs)

if (!isGeneric("se<-"))
  setGeneric("se<-",function(object,value)
             standardGeneric("se<-"))


  #replace coefs (chip level coefficients)
setReplaceMethod("se",signature(object="PLMset"),
                 function(object,value){
                   object@se.chip.coefs <- value
                   object
                 })  

## indexProbes, similar to that used in the AffyBatch class
  ## use the cdfenv to get what we need.
  
if( !isGeneric("indexProbes") )
  setGeneric("indexProbes", function(object, which, ...)
             standardGeneric("indexProbes"))

setMethod("indexProbes", signature("PLMset", which="character"),
          function(object, which=c("pm", "mm","both"),
                   genenames=NULL, xy=FALSE) {
            
            which <- match.arg(which)
            
            i.probes <- match(which, c("pm", "mm", "both"))
            ## i.probes will know if "[,1]" or "[,2]"
            ## if both then [,c(1,2)]
            if(i.probes==3) i.probes=c(1,2)
            
            envir <- getCdfInfo(object)
            
            if(is.null(genenames)) 
              genenames <- ls(envir )
            
            ## shorter code, using the features of multiget
            ## (eventually more readable too)
            ## note: genenames could be confusing (the same gene can be
              ## found in several affyid (ex: the 3' and 5' controls)
            
            ans <-  mget(genenames, envir, ifnotfound=NA)
            
            ## this kind of thing could be included in 'multiget' as
            ## and extra feature. A function could be specified to
            ## process what is 'multiget' on the fly
            for (i in seq(along=ans)) {
              
              
                                        #this line needs to be changed for R 1.7.0
              if ( is.na(ans[[i]][1]) )
                next
              
              ##as.vector cause it might be a matrix if both
              tmp <- as.vector(ans[[i]][, i.probes])
              
              
              if (xy) {
                warning("flag 'xy' is deprecated")
                x <- tmp %% nrow(object)
                x[x == 0] <- nrow(object)
                y <- tmp %/% nrow(object) + 1
                tmp <- cbind(x, y)
              }
              
              ans[[i]] <- tmp
            }
            
            return(ans)
          })


if( !isGeneric("indexProbesProcessed") )
  setGeneric("indexProbesProcessed", function(object)
             standardGeneric("indexProbesProcessed"))

setMethod("indexProbesProcessed", signature("PLMset"),
	function(object){
		pmindex <-indexProbes(object,which="pm")	
		pmindex.length <- lapply(pmindex,length)

		cs <- cumsum(do.call("c",pmindex.length)) 
		cl  <- do.call("c",pmindex.length)
		for (i in 1:length(pmindex)){
			pmindex[[i]] <- cs[i] - (cl[i]:1)+1

		}
		return(pmindex)
	})





  

    
#  if( !isGeneric("image.weights") )
#    setGeneric("image.weights", function(x)
#               standardGeneric("image.weights"), where=where)

    
setMethod("image",signature(x="PLMset"),
          function(x,which=0,type=c("weights","resids","pos.resids","neg.resids","sign.resids"),use.log=TRUE,add.legend=FALSE,standardize=FALSE,col=NULL,...){


            if (is.null(col)){
              col.weights <- terrain.colors(25)
              col.resids <- pseudoPalette(low="blue",high="red",mid="white")
              col.pos.resids <- pseudoPalette(low="white",high="red")
              col.neg.resids <- pseudoPalette(low="blue",high="white")
            } else {
              col.weights <- col
              col.resids <- col
              col.pos.resids <-  col
              col.neg.resids <-  col
            }


       

            
            
            type <- match.arg(type)
            
            pm.index <- unique(unlist(indexProbes(x, "pm",row.names(coefs(x)))))
            rows <- x@nrow
            cols <- x@ncol
            pm.x.locs <- pm.index%%rows
            pm.x.locs[pm.x.locs == 0] <- rows
            pm.y.locs <- pm.index%/%rows + 1
            xycoor <- matrix(cbind(pm.x.locs,pm.y.locs),ncol=2)
            xycoor2 <- matrix(cbind(pm.x.locs,pm.y.locs+1),ncol=2)
            
            
            if (is.element(type,c("weights"))){
              if (any(dim(x@weights[[1]]) ==0) & any(dim(x@weights[[2]]) ==0)){
                stop("Sorry this PLMset does not appear to have weights\n");
              } 
              if (which == 0){
                which <- 1:max(dim(x@weights[[1]])[2], dim(x@weights[[2]])[2])
              }
            }
            
            if (is.element(type,c("resids","pos.resids","neg.resids","sign.resids"))){
              if (any(dim(x@residuals[[1]]) ==0) & any(dim(x@residuals[[2]]) ==0)){
                stop("Sorry this PLMset does not appear to have residuals\n");
              }
              if (which == 0){
                which <- 1:max(dim(x@residuals[[1]])[2],dim(x@residuals[[2]])[2])
              }
              if (standardize & type == "resids"){
                if (x@model.description$R.model$response.variable == 0){
                  resid.range <- c(-4,4)
                } else if (x@model.description$R.model$response.variable == -1){
                  resid.range <- range(resid(x,standardize)[[2]])
                } else if (x@model.description$R.model$response.variable == 1){
                  resid.range <- range(resid(x,standardize)[[1]])
                }
                
              } else {
                  if (x@model.description$R.model$response.variable == 0){
                    resid.range1 <- range(x@residuals[[1]])
                    resid.range2 <- range(x@residuals[[2]])
                    resid.range <- resid.range1
                    resid.range[1] <- min(resid.range1 , resid.range2)
                    resid.range[2] <- max(resid.range1 , resid.range2)
                  } else if (x@model.description$R.model$response.variable == -1){
                    resid.range <- range(x@residuals[[2]])
                  } else if (x@model.description$R.model$response.variable == 1){
                    resid.range <- range(x@residuals[[1]])
                  }
                }
            }
            
            
            
            for (i in which){
              if (type == "weights"){
                weightmatrix <-matrix(nrow=rows,ncol=cols)
                if (x@model.description$R.model$response.variable == 0){
                  weightmatrix[xycoor]<- x@weights[[1]][,i]
                  weightmatrix[xycoor2]<- x@weights[[2]][,i]
                } else if (x@model.description$R.model$response.variable == -1){
                  weightmatrix[xycoor]<- x@weights[[2]][,i]
                  weightmatrix[xycoor2]<- x@weights[[2]][,i]
                } else if (x@model.description$R.model$response.variable == 1){
                  weightmatrix[xycoor]<- x@weights[[1]][,i]
                  weightmatrix[xycoor2]<- x@weights[[1]][,i]
                }


                
                                        #this line flips the matrix around so it is correct
                weightmatrix <-as.matrix(rev(as.data.frame(weightmatrix)))
                if (add.legend){
                  layout(matrix(c(1, 2), 1, 2), width = c(9, 1))
                  par(mar = c(4, 4, 5, 3))
                }
                image(weightmatrix,col=col.weights,xaxt='n',
                      yaxt='n',main=sampleNames(x)[i],zlim=c(0,1))
                title(sampleNames(x)[i])
                if (add.legend){
                  par(mar = c(4, 0, 5, 3))
                  pseudoColorBar(seq(0,1,0.1), horizontal = FALSE, col = col.weights, main = "")
                  layout(1)
                  par(mar = c(5, 4, 4, 2) + 0.1)
                }
                
              }
              if (type == "resids"){
                residsmatrix <- matrix(nrow=rows,ncol=cols)
                if (standardize){
                  if (x@model.description$R.model$response.variable == 0){
                    residsmatrix[xycoor]<- resid(x,standardize)[[1]][,i]
                    residsmatrix[xycoor2]<- resid(x,standardize)[[2]][,i]
                  } else if  (x@model.description$R.model$response.variable == -1){
                    residsmatrix[xycoor]<- resid(x,standardize)[[2]][,i]
                    residsmatrix[xycoor2]<- resid(x,standardize)[[2]][,i]
                  } else if (x@model.description$R.model$response.variable == 1){
                    residsmatrix[xycoor]<- resid(x,standardize)[[1]][,i]
                    residsmatrix[xycoor2]<- resid(x,standardize)[[1]][,i]
                  }
                } else {
                  if (x@model.description$R.model$response.variable == 0){
                    residsmatrix[xycoor]<- x@residuals[[1]][,i]
                    residsmatrix[xycoor2]<- x@residuals[[2]][,i]
                  } else if (x@model.description$R.model$response.variable == -1){
                    residsmatrix[xycoor]<- x@residuals[[2]][,i]
                    residsmatrix[xycoor2]<- x@residuals[[2]][,i]
                  } else if (x@model.description$R.model$response.variable == 1){
                    residsmatrix[xycoor]<- x@residuals[[1]][,i]
                    residsmatrix[xycoor2]<- x@residuals[[1]][,i]
                  }
                    
                }
                                        #this line
                                        #flips the matrix around so it is correct
                residsmatrix<- as.matrix(rev(as.data.frame(residsmatrix)))
              ##  if (standardize){
              ##    if (x@model.description$R.model$response.variable == 0){
              ##      resid.range <- c(-4,4)
              ##    } else if (x@model.description$R.model$response.variable == -1){
              ##      resid.range <- range(resid(x,standardize)[[2]])
              ##    } else if (x@model.description$R.model$response.variable == 1){
              ##      resid.range <- range(resid(x,standardize)[[1]])
              ##    }
              ##      
              ##  } else {
              ##    if (x@model.description$R.model$response.variable == 0){
              ##      resid.range1 <- range(x@residuals[[1]])
              ##      resid.range2 <- range(x@residuals[[2]])
              ##      resid.range <- resid.range1
              ##      resid.range[1] <- min(resid.range1 , resid.range2)
              ##      resid.range[2] <- max(resid.range1 , resid.range2)
              ##    } else if (x@model.description$R.model$response.variable == -1){
              ##      resid.range <- range(x@residuals[[2]])
              ##    } else if (x@model.description$R.model$response.variable == 1){
              ##      resid.range <- range(x@residuals[[1]])
              ##    }
              ##  }
                if (use.log){
                  if (add.legend){
                    layout(matrix(c(1, 2), 1, 2), width = c(9, 1))
                    par(mar = c(4, 4, 5, 3))
                  }
                  residsmatrix <- sign(residsmatrix)*log2(abs(residsmatrix)+1)
                  image(residsmatrix,col=col.resids,xaxt='n',
                        yaxt='n',main=sampleNames(x)[i],zlim=c(-max(log2(abs(resid.range)+1)),max(log2(abs(resid.range)+1))))
                  if (add.legend){
                    par(mar = c(4, 0, 5, 3))
                    pseudoColorBar(seq(-max(log2(abs(resid.range)+1)),max(log2(abs(resid.range)+1)),0.1), horizontal = FALSE, col = col.resids, main = "",log.ticks=TRUE)
                    layout(1)
                    par(mar = c(5, 4, 4, 2) + 0.1)
                  } 
                  
                  
                  
                  
                } else {
                  
                  if (add.legend){
                    layout(matrix(c(1, 2), 1, 2), width = c(9, 1))
                    par(mar = c(4, 4, 5, 3))
                  }
                  image(residsmatrix,col=col.resids,xaxt='n',
                        yaxt='n',main=sampleNames(x)[i],zlim=c(-max(abs(resid.range)),max(abs(resid.range))))
                  if (add.legend){
                    par(mar = c(4, 0, 5, 3))
                    pseudoColorBar(seq(-max(abs(resid.range)),max(abs(resid.range)),0.1), horizontal = FALSE, col = col.resids, main = "")
                    layout(1)
                    par(mar = c(5, 4, 4, 2) + 0.1)
                  } 
                }
              }
              if (type == "pos.resids"){
                residsmatrix <- matrix(nrow=rows,ncol=cols)

                if (x@model.description$R.model$response.variable == 0){
                  residsmatrix[xycoor]<- pmax(x@residuals[[1]][,i],0)
                  residsmatrix[xycoor2]<- pmax(x@residuals[[2]][,i],0)
                } else if (x@model.description$R.model$response.variable == -1){
                  residsmatrix[xycoor]<- pmax(x@residuals[[2]][,i],0)
                  residsmatrix[xycoor2]<- pmax(x@residuals[[2]][,i],0)
                } else if (x@model.description$R.model$response.variable == 1){
                  residsmatrix[xycoor]<- pmax(x@residuals[[1]][,i],0)
                  residsmatrix[xycoor2]<- pmax(x@residuals[[1]][,i],0)
                }



                
                                        #this                 line flips the matrix around so it is correct
                residsmatrix <- as.matrix(rev(as.data.frame(residsmatrix)))


               ## if (x@model.description$R.model$response.variable == 0){
               ##   resid.range1 <- range(x@residuals[[1]])
               ##   resid.range2 <- range(x@residuals[[2]])
               ##   resid.range <- resid.range1
               ##   resid.range[1] <- min(resid.range1 , resid.range2)
               ##   resid.range[2] <- max(resid.range1 , resid.range2)
               ## } else if (x@model.description$R.model$response.variable == -1){
               ##   resid.range <- range(x@residuals[[2]])
               ## } else if (x@model.description$R.model$response.variable == 1){
               ##   resid.range <- range(x@residuals[[1]])
               ## }

                if (use.log){
                  if (add.legend){
                    layout(matrix(c(1, 2), 1, 2), width = c(9, 1))
                    par(mar = c(4, 4, 5, 3))
                  }
                  residsmatrix <- sign(residsmatrix)*log2(abs(residsmatrix) +1)
                  image(residsmatrix,col=col.pos.resids,xaxt='n',
                        yaxt='n',main=sampleNames(x)[i],zlim=c(0,max(log2(pmax(resid.range,0)+1))))
                  if (add.legend){
                    par(mar = c(4, 0, 5, 3))
                    pseudoColorBar(seq(0,max(log2(pmax(resid.range,0)+1)),0.1), horizontal = FALSE, col = col.pos.resids, main = "",log.ticks=TRUE)
                    layout(1)
                    par(mar = c(5, 4, 4, 2) + 0.1)
                  } 
                } else {
                  if (add.legend){
                    layout(matrix(c(1, 2), 1, 2), width = c(9, 1))
                    par(mar = c(4, 4, 5, 3))
                  }
                  image(residsmatrix,col=col.pos.resids,xaxt='n',
                        yaxt='n',main=sampleNames(x)[i],zlim=c(0,max(resid.range)))
                  if (add.legend){
                    par(mar = c(4, 0, 5, 3))
                    pseudoColorBar(seq(0,max(resid.range),0.1), horizontal = FALSE, col = col.pos.resids, main = "")
                    layout(1)
                    par(mar = c(5, 4, 4, 2) + 0.1)
                  } 
                }
              }
              if (type == "neg.resids"){
                residsmatrix <- matrix(nrow=rows,ncol=cols)
                if (x@model.description$R.model$response.variable == 0){
                  residsmatrix[xycoor]<- pmin(x@residuals[[1]][,i],0)
                  residsmatrix[xycoor2]<- pmin(x@residuals[[2]][,i],0)
                } else if (x@model.description$R.model$response.variable == -1){
                  residsmatrix[xycoor]<- pmin(x@residuals[[2]][,i],0)
                  residsmatrix[xycoor2]<- pmin(x@residuals[[2]][,i],0)
                } else if (x@model.description$R.model$response.variable == 1){
                  residsmatrix[xycoor]<- pmin(x@residuals[[1]][,i],0)
                  residsmatrix[xycoor2]<- pmin(x@residuals[[1]][,i],0)
                }


                  
                                        #this line flips the matrix around so it is correct
                residsmatrix <- as.matrix(rev(as.data.frame(residsmatrix)))

               ## if (x@model.description$R.model$response.variable == 0){
               ##   resid.range1 <- range(x@residuals[[1]])
               ##   resid.range2 <- range(x@residuals[[2]])
               ##   resid.range <- resid.range1
               ##   resid.range[1] <- min(resid.range1 , resid.range2)
               ##   resid.range[2] <- max(resid.range1 , resid.range2)
               ## } else if (x@model.description$R.model$response.variable == -1){
               ##   resid.range <- range(x@residuals[[2]])
               ## } else if (x@model.description$R.model$response.variable == 1){
               ##   resid.range <- range(x@residuals[[1]])
               ## }

                if(use.log){
                  if (add.legend){
                    layout(matrix(c(1, 2), 1, 2), width = c(9, 1))
                    par(mar = c(4, 4, 5, 3))
                  }
                  residsmatrix <- sign(residsmatrix)*log2(abs(residsmatrix) +1)
                  image(residsmatrix,col=col.neg.resids,xaxt='n',
                        yaxt='n',main=sampleNames(x)[i],zlim=c(-log2(abs(min(resid.range))+1),0))
                  if (add.legend){
                    par(mar = c(4, 0, 5, 3))
                    pseudoColorBar(seq(-max(log2(abs(pmin(resid.range,0))+1)),0,0.1), horizontal = FALSE, col = col.neg.resids, main = "",log.ticks=TRUE)
                    layout(1)
                    par(mar = c(5, 4, 4, 2) + 0.1)
                  } 
                  
                } else {
                  if (add.legend){
                    layout(matrix(c(1, 2), 1, 2), width = c(9, 1))
                    par(mar = c(4, 4, 5, 3))
                  }
                  image(residsmatrix,col=col.neg.resids,xaxt='n',
                        yaxt='n',main=sampleNames(x)[i],zlim=c(-abs(min(resid.range)),0))
                  if (add.legend){
                    par(mar = c(4, 0, 5, 3))
                    pseudoColorBar(seq(min(resid.range),0,0.1), horizontal = FALSE, col = col.neg.resids, main = "")
                    layout(1)
                    par(mar = c(5, 4, 4, 2) + 0.1)
                  } 
                }
                
              }
              if (type == "sign.resids"){

                residsmatrix <- matrix(nrow=rows,ncol=cols)
                if (x@model.description$R.model$response.variable == 0){
                  residsmatrix[xycoor]<- sign(x@residuals[[1]][,i])
                  residsmatrix[xycoor2]<- sign(x@residuals[[2]][,i])
                } else if (x@model.description$R.model$response.variable == -1){
                  residsmatrix[xycoor]<- sign(x@residuals[[2]][,i])
                  residsmatrix[xycoor2]<- sign(x@residuals[[2]][,i])
                } else if (x@model.description$R.model$response.variable == 1){
                  residsmatrix[xycoor]<- sign(x@residuals[[1]][,i])
                  residsmatrix[xycoor2]<- sign(x@residuals[[1]][,i])
                }

                                        #this line flips the matrix around so it is correct
                residsmatrix <- as.matrix(rev(as.data.frame(residsmatrix)))

                if (add.legend){
                  layout(matrix(c(1, 2), 1, 2), width = c(9, 1))
                  par(mar = c(4, 4, 5, 3))
                }
                image(residsmatrix,col=col.resids,xaxt='n',
                      yaxt='n',main=sampleNames(x)[i],zlim=c(-1,1))
                if (add.legend){
                  par(mar = c(4, 0, 5, 3))
                  pseudoColorBar(seq(-1,1,2), horizontal = FALSE, col = col.resids, main = "")
                  layout(1)
                  par(mar = c(5, 4, 4, 2) + 0.1)
                } 
                
              }
              

              
            }
          })


 

setMethod("boxplot",signature(x="PLMset"),
          function(x,type=c("NUSE","weights","residuals"),range=0,...){
           

            compute.nuse <- function(which){
              nuse <- apply(x@weights[which,],2,sum)
              1/sqrt(nuse)
            }
            
            
            type <- match.arg(type)
            model <- x@model.description$modelsettings$model
            if (type == "NUSE"){
              if (x@model.description$R.model$which.parameter.types[3] == 1 & x@model.description$R.model$which.parameter.types[1] == 0 ){
                grp.rma.se1.median <- apply(se(x), 1,median)
                grp.rma.rel.se1.mtx <- sweep(se(x),1,grp.rma.se1.median,FUN='/')
                boxplot(data.frame(grp.rma.rel.se1.mtx),range=range,...)
              } else {
                # not the default model try constructing them using weights.
                which <-indexProbesProcessed(x)
                ses <- matrix(0,length(which) ,4)
                
                for (i in 1:length(which))
                  ses[i,] <- compute.nuse(which[[i]])
                
                
                grp.rma.se1.median <- apply(ses, 1,median)
                grp.rma.rel.se1.mtx <- sweep(ses,1,grp.rma.se1.median,FUN='/')
                boxplot(data.frame(grp.rma.rel.se1.mtx),range=range,...)
              }
            } else if (type == "weights"){
              ow <- options("warn")
              options(warn=-1)
              if (x@model.description$R.model$response.variable == -1){
                boxplot(data.frame(x@weights[[2]]),...)
              } else if (x@model.description$R.model$response.variable == 1){
                boxplot(data.frame(x@weights[[1]]),...)
              } else {
                boxplot(data.frame(rbind(x@weights[[1]],x@weights[[2]])),...)
              }
              options(ow)
            } else if (type == "residuals"){
              ow <- options("warn")
              options(warn=-1)
              if (x@model.description$R.model$response.variable == -1){
                boxplot(data.frame(x@residuals[[2]]),...)
              } else if (x@model.description$R.model$response.variable == 1){
                boxplot(data.frame(x@wresiduals[[1]]),...)
              } else {
                boxplot(data.frame(rbind(x@residuals[[1]],x@residuals[[2]])),...)
              }
              options(ow)
            }
          })


setMethod("show", "PLMset",
          function(object) {
            
            cat("Probe level linear model (PLMset) object\n")
            cat("size of arrays=", object@nrow, "x", object@ncol,"\n",sep="")
            
            ## Location from cdf env
            try( cdf.env <- getCdfInfo(object) )
            if (! inherits(cdf.env, "try-error")) {
              num.ids <- length(ls(env=cdf.env))
            } else {
              warning("missing cdf environment !")
              num.ids <- "???"
            }
            
            cat("cdf=", object@cdfName,
                " (", num.ids, " probeset ids)\n",
                sep="")
            cat("number of samples=",dim(object@weights)[2],"\n",sep="")
            cat("number of probesets=", num.ids, "\n",sep="")
            cat("number of chip level parameters for each probeset=",dim(object@chip.coefs)[2],"\n")
            cat("annotation=",object@annotation,"\n",sep="")
            cat("notes=",object@notes,"\n\n",sep="")
            cat("PLMset settings\n")
            cat("Creating function:",object@model.description$which.function,"\n")
            cat("Preprocessing\n")
            cat("Background Correction=",object@model.description$preprocessing$background,sep="")
            if (object@model.description$preprocessing$background){
              cat(" Method=",object@model.description$preprocessing$bg.method)
            }
            cat("\n")
            
            cat("Normalization=",object@model.description$preprocessing$normalize,sep="")
            if (object@model.description$preprocessing$normalize){
              cat(" Method=",object@model.description$preprocessing$norm.method)
            }
            cat("\n")

            cat("\nModel/Summarization\n")
            print(object@model.description$modelsettings)
            cat("\n")
            cat("Output Settings\n")
            print(object@model.description$outputsettings)
            
            
          })

if (!isGeneric("coefs.const"))
  setGeneric("coefs.const",function(object)
             standardGeneric("coefs.const"))
  
  setMethod("coefs.const","PLMset",
            function(object){
              exprs(object)
            })


if (!isGeneric("se.const"))
  setGeneric("se.const",function(object)
             standardGeneric("se.const"))

setMethod("se.const","PLMset",
          function(object){
            se.exprs(object)
          })

#A summary method, to be cleaned up better at a later date.
 
setMethod("summary","PLMset",
          function(object,genenames=NULL){#

              if (is.null(genenames)){
                genenames <- rownames(object@chip.coefs)
              }
              cur.const.coef <-  NULL
              cur.const.se <- NULL

              allindexs <- indexProbesProcessed(object)
              for (probeset.names in genenames){
                if (all(dim( exprs(object)) != 0)){
                  cur.const.coef <- exprs(object)[grep(paste("^",probeset.names,sep=""),rownames(object@chip.coefs))]
                  cur.const.se <-  se.exprs(object)[grep(paste("^",probeset.names,sep=""),rownames(object@chip.coefs))]
                }
                inds <- allindexs[probeset.names]
                inds <- do.call("c",inds)
                cur.probe.coef <- object@probe.coefs[inds,]
                cur.se.probe.coef <- object@se.probe.coefs[inds,]
                cur.chip.coef <- object@chip.coefs[grep(paste("^",probeset.names,sep=""),rownames(object@chip.coefs)),]
                cur.chip.se <- object@se.chip.coefs[grep(paste("^",probeset.names,sep=""),rownames(object@se.chip.coefs)),]#

                
                cat("Probeset:", probeset.names,"\n")

                cat("Intercept Estimates\n")
                print(cbind(Coef=cur.const.coef,SE=cur.const.se))
                cat("\n")
                cat("Chip Effect Estimates\n")
                print(cbind(Coef=cur.chip.coef,SE=cur.chip.se))


                cat("\n")
                cat("Probe Effect Estimates\n")
                print(cbind(Coef=cur.probe.coef,SE=cur.se.probe.coef))

                cat("\nResiduals\n")
                print(object@residuals[inds,])
                
                 cat("\nWeights\n")
                print(object@weights[inds,])
                cat("\n\n")
              }
            })


if (!isGeneric("Mbox"))
  setGeneric("Mbox",function(object,...)
             standardGeneric("Mbox"))
  

  
setMethod("Mbox",signature("PLMset"),
          function(object,range=0,...){
            if (object@model.description$R.model$which.parameter.types[3] == 1){
              medianchip <- apply(coefs(object), 1, median)
              M <- sweep(coefs(object),1,medianchip,FUN='-')
              boxplot(data.frame(M),range=range,...)
            } else {
              stop("It doesn't appear that a model with sample effects was used.")
            }
          })



if (!isGeneric("resid<-"))
  setGeneric("resid<-",function(object,value)
             standardGeneric("resid<-"))


setReplaceMethod("resid",signature(object="PLMset"),
                 function(object,value){
                   object@residuals <- value
                   object
                 })




setMethod("resid",signature("PLMset"),
          function(object,genenames=NULL,standardize=FALSE){
	    if (!standardize){
              if (is.null(genenames)){	 	
                object@residuals
              } else {
                which <-indexProbesProcessed(object)[genenames]
                which <- do.call("c",which)
                if (object@model.description$R.model$response.variable == 0){
                  list(PM.resid=object@residuals[[1]][which,],MM.resid=object@residuals[[2]][which,])
                } else if (object@model.description$R.model$response.variable == -1){
                  list(PM.resid=matrix(0,0,0),MM.resid=object@residuals[[2]][which,])
                } else if (object@model.description$R.model$response.variable == 1){
                  list(PM.resid=object@residuals[[1]][which,],MM.resid=matrix(0,0,0))
                }
              }
	    } else {
              which <-indexProbesProcessed(object)
              if (!is.null(genenames)){
                which <- which[genenames]
              }
              if (object@model.description$R.model$response.variable == 0){
                results1 <- lapply(which,function(rowindex, x){
                  x[rowindex,]   
                },object@residuals[[1]])
                results2 <- lapply(which,function(rowindex, x){
                  x[rowindex,] 
                },object@residuals[[2]])
                for (i in 1:length(results1)){
                  cur.sd <- sd(c(as.vector(results1[[i]]),as.vector(results2[[i]])))
                  cur.mean <- mean(c(as.vector(results1[[i]]),as.vector(results2[[i]])))
                  results1[[i]] <- (results1[[i]]-cur.mean)/cur.sd
                  results2[[i]] <- (results2[[i]]-cur.mean)/cur.sd
                }
                return(list(PM.resid=do.call("rbind",results1),MM.resid=do.call("rbind",results2)))
              } else if (object@model.description$R.model$response.variable == -1){
                results <- lapply(which,function(rowindex, x){
                  (x[rowindex,]- mean(as.vector(x[rowindex,])))/sd(as.vector(x[rowindex,]))
                },object@residuals[[2]])
                return(list(PM.resid=matrix(0,0,0),MM.resid=do.call("rbind",results)))
              } else if (object@model.description$R.model$response.variable == 1){
                results <- lapply(which,function(rowindex, x){
                  (x[rowindex,]- mean(as.vector(x[rowindex,])))/sd(as.vector(x[rowindex,]))
                },object@residuals[[1]])
                return(list(PM.resid=do.call("rbind",results),MM.resid=matrix(0,0,0)))

              }
            }
          })


if (!isGeneric("residuals<-"))
  setGeneric("residuals<-",function(object,value)
             standardGeneric("residuals<-"))


setReplaceMethod("residuals",signature(object="PLMset"),
                 function(object,value){
                   object@residuals <- value
                   object
                 })


setMethod("residuals",signature("PLMset"),
            function(object,genenames=NULL,standardize=FALSE){
              resid(object,genenames,standardize)
	    })

if (!isGeneric("normvec"))
  setGeneric("normvec",function(object,...)
             standardGeneric("normvec"))
    
    
setMethod("normvec",signature("PLMset"),
          function(object){
            object@normVec
          })

if (!isGeneric("varcov"))
  setGeneric("varcov",function(object,...)
             standardGeneric("varcov"))

  
  setMethod("varcov",signature("PLMset"),
            function(object,...){
              object@varcov
            })
  
  

if (!isGeneric("residSE"))
  setGeneric("residSE",function(object,...)
             standardGeneric("residSE"))


  
setMethod("residSE",signature("PLMset"),
          function(object){
            return(object@residualSE)
          })



setMethod("sampleNames",signature("PLMset"),function(object){
  rownames(pData(object))
})



if (!isGeneric("sampleNames<-"))
  setGeneric("sampleNames<-",function(object,value)
             standardGeneric("sampleNames<-"))


setReplaceMethod("sampleNames",signature(object="PLMset"),
                 function(object,value){
                   rownames(pData(object)) <- value
		   if (!any(dim(object@weights) == 0)){
			colnames(object@weights) <- value			
		   }
		    if (!any(dim(object@residuals) == 0)){
			colnames(object@residuals) <- value			
		   }
		   object         
                 })




###################################################################
##
## model.description is a list
## $which.function - character string specifying name of function
##                   used to generate PLMset
## $preprocessing - a list of information for preprocessing
##             $bg.method - character.string
##             $bg.param - a list of settings relevant to bgc
##             $background - logical TRUE if background correction
##             $norm.method - character string
##             $norm.param - a list of settings relevant to normalization
##             $normalization -logical if normalization took places
## $modelsettings - a list of information related to summary/PLM model
##             $model.param - list of settings used
##             $summary.method - character string
##             $model - in the case of fitPLM, the model should be specified here, otherwise empty
##             $constraint.type - vector listing constraint's on terms in the model (fitPLM case)
##             $variable type - vector defining whether variables are factors or covariates (fitPLM case)
## $outputsettings - a list of output settings
##
##
##
##    fitPLM(object,model=PM ~ -1 + probes +samples,
##     variable.type=c(default="factor"),
##     constraint.type=c(default="contr.treatment"),
##     background=TRUE, normalize=TRUE, background.method = "RMA.2",normalize.method = "quantile",
##       background.param=list(),normalize.param=list(),output.param=list(),model.param=list())
##
##threestepPLM(object, normalize=TRUE,background=TRUE,background.method="RMA.2",normalize.method="quantile",summary.method="median.polish",background.param = list(),normalize.param=list(),output.param=list(), model.param=list())
##
##       rmaPLM(object,normalize=TRUE,background=TRUE,background.method="RMA.2",normalize.method="quantile",background.param = list(),normalize.param=list(),output.param=list(),model.param=list())
##
##
###################################################################



if (!isGeneric("model.description"))
  setGeneric("model.description",function(object,...)
             standardGeneric("model.description"))



setMethod("model.description", "PLMset", function(object)
          object@model.description)












pseudoPalette <-function (low = "white", high = c("green", "red"), mid = NULL,
                      k = 50)
{
    low <- col2rgb(low)/255
    high <- col2rgb(high)/255
    if (is.null(mid)) {
      r <- seq(low[1], high[1], len = k)
      g <- seq(low[2], high[2], len = k)
      b <- seq(low[3], high[3], len = k)
    }
    if (!is.null(mid)) {
        k2 <- round(k/2)
        mid <- col2rgb(mid)/255
        r <- c(seq(low[1], mid[1], len = k2), seq(mid[1], high[1],
            len = k2))
        g <- c(seq(low[2], mid[2], len = k2), seq(mid[2], high[2],
            len = k2))
        b <- c(seq(low[3], mid[3], len = k2), seq(mid[3], high[3],
            len = k2))
    }
    rgb(r, g, b)
  }

pseudoColorBar <- function (x, horizontal = TRUE, col = heat.colors(50), scale = 1:length(x),
    k = 11, log.ticks=FALSE,...)
{
    if (is.numeric(x)) {
        x <- x
        colmap <- col
    }
    else {
      colmap <- x
      low <- range(scale)[1]
      high <- range(scale)[2]
      x <- seq(low, high, length = length(x))
    }
    if (length(x) > k){
      x.small <- seq(x[1], x[length(x)], length = k)
      if (log.ticks){
        x.small <- sign(x.small)*(2^abs(x.small) -1)
        x <- sign(x)*(2^abs(x) -1)
      }
    }
    else{
      x.small <- x
      if (log.ticks){
        x.small <- sign(x.small)*(2^abs(x.small) -1)
        x <- sign(x)*(2^abs(x) -1)
      }
    }
    if (horizontal) {
        image(x, 1, matrix(x, length(x), 1), axes = FALSE, xlab = "",
            ylab = "", col = colmap, ...)
        axis(1, at = rev(x.small), labels = signif(rev(x.small),
            2), srt = 270)
    }
    if (!horizontal) {
      image(1, x, matrix(x, 1, length(x)), axes = FALSE, xlab = "",
            ylab = "", col = colmap, ...)
      par(las = 1)
      axis(4, at = rev(x.small), labels = signif(rev(x.small),2))
      par(las = 0)
    }
    box()
}



if (!isGeneric("MAplot"))
  setGeneric("MAplot",function(object,...)
             standardGeneric("MAplot"))


setMethod("MAplot",signature("PLMset"),
          function(object,ref=NULL,subset=NULL,which=NULL,...){
            x <- coefs(object)
            if (is.null(which)){
              which <- 1:dim(x)[2]
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



if (!isGeneric("nuse"))
  setGeneric("nuse",function(x,...)
             standardGeneric("nuse"))



setMethod("nuse",signature(x="PLMset"),
          function(x,type=c("plot","values","stats"),...){


            compute.nuse <- function(which){
              nuse <- apply(x@weights[which,],2,sum)
              1/sqrt(nuse)
            }

            type <- match.arg(type)
            model <- x@model.description$modelsettings$model
            if (type == "values" || type == "stats" ){
              
              if (x@model.description$R.model$which.parameter.types[3] == 1 & x@model.description$R.model$which.parameter.types[1] == 0 ){
                grp.rma.se1.median <- apply(se(x), 1,median)
                grp.rma.rel.se1.mtx <- sweep(se(x),1,grp.rma.se1.median,FUN='/')
              } else {
                                        # not the default model try constructing them using weights.
                which <-indexProbesProcessed(x)
                ses <- matrix(0,length(which) ,4)
                
                for (i in 1:length(which))
                  ses[i,] <- compute.nuse(which[[i]])
                
                
                grp.rma.se1.median <- apply(ses, 1,median)
                grp.rma.rel.se1.mtx <- sweep(ses,1,grp.rma.se1.median,FUN='/')
              }
              if (type == "values"){
                grp.rma.rel.se1.mtx
              } else {
                Medians <- apply(grp.rma.rel.se1.mtx,2,median)
                Quantiles <- apply(grp.rma.rel.se1.mtx,2,quantile,prob=c(0.25,0.75))
                nuse.stats <- rbind(Medians,Quantiles[2,] - Quantiles[1,])
                rownames(nuse.stats) <- c("median","IQR")
                nuse.stats
              }
            } else {
               boxplot(x,...)
            }
          })

if (!isGeneric("NUSE"))
  setGeneric("NUSE",function(x,...)
             standardGeneric("NUSE"))



setMethod("NUSE",signature(x="PLMset"),
          function(x,type=c("plot","values","stats"),ylim=c(0.9,1.2),add.line=TRUE,...){
           type <- match.arg(type)
            x <- nuse(x,type=type,ylim=ylim,...)
            if (add.line & (type == "plot")){
              abline(1,0)
            } else{
              x
            }
          })







            
if (!isGeneric("RLE"))
  setGeneric("RLE",function(x,...)
             standardGeneric("RLE"))




setMethod("RLE",signature(x="PLMset"),
            function(x,type=c("plot","values","stats"),ylim=c(-0.75,0.75),add.line=TRUE,...){


              type <- match.arg(type)
              model <- x@model.description$modelsettings$model
              if (type == "values" || type=="stats"){
                if (x@model.description$R.model$which.parameter.types[3] == 1){
                  medianchip <- apply(coefs(x), 1, median)
                  if (type == "values"){
                    sweep(coefs(x),1,medianchip,FUN='-')
                  } else {
                    RLE <- sweep(coefs(x),1,medianchip,FUN='-')
                    Medians <- apply(RLE,2,median)
                    Quantiles <- apply(RLE,2,quantile,prob=c(0.25,0.75))
                    RLE.stats <- rbind(Medians,Quantiles[2,] - Quantiles[1,])
                    rownames(RLE.stats) <- c("median","IQR")
                    RLE.stats
                  }
                } else {
                  stop("It doesn't appear that a model with sample effects was used.")
                }
              } else {
                Mbox(x,ylim=ylim,...)
                if (add.line){
                  abline(0,0)
                }
              }
            })
