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
## Sep 25, 2003 - Port to R-1.8
## Oct 7, 2003 - update getCDFInfo to reflect changes in affy package.
##
###########################################################

#creating the PLMset object

setClass("PLMset",
           representation(probe.coefs="matrix",
                          se.probe.coefs="matrix",
                          chip.coefs = "matrix",
                          se.chip.coefs = "matrix",
                          cdfName="character",
                          nrow="numeric",
                          ncol="numeric",
                          model.description="character",
                          model.call = "call",
                          weights="matrix"),
                         # phenoData="phenoData",
                         # description="characterORMIAME",
                         # annotation="character",
                         # notes="character"
           prototype=list(
             probe.coefs=matrix(nr=0,nc=0),
             se.probe.coefs=matrix(nr=0,nc=0),
             chip.coefs=matrix(nr=0,nc=0),
             se.chip.coefs=matrix(nr=0,nc=0),
             model.description="",
             weights=matrix(nr=0,nc=0),
             description=new("MIAME"),
             annotation="",
             cdfName="",
             nrow=0, ncol=0,
             notes=""),contains="exprSet")


  #now some accessors.

  #access weights
  setMethod("weights",signature(object="PLMset"),
            function(object) object@weights)



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
                genenames <- ls(envir)
              
              ## shorter code, using the features of multiget
              ## (eventually more readable too)
              ## note: genenames could be confusing (the same gene can be
              ## found in several affyid (ex: the 3' and 5' controls)
              
              ans <-  multiget(genenames, pos, envir, iffail=NA)

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
  

    
#  if( !isGeneric("image.weights") )
#    setGeneric("image.weights", function(x)
#               standardGeneric("image.weights"), where=where)

    
  setMethod("image",signature(x="PLMset"),
            function(x,which=0,...){
              pm.index <- unique(unlist(indexProbes(x, "pm")))
              rows <- x@nrow
              cols <- x@ncol
              pm.x.locs <- pm.index%%rows
              pm.x.locs[pm.x.locs == 0] <- rows
              pm.y.locs <- pm.index%/%rows + 1
              xycoor <- matrix(cbind(pm.x.locs,pm.y.locs),ncol=2)
              xycoor2 <- matrix(cbind(pm.x.locs,pm.y.locs+1),ncol=2)

              if (which == 0){
                which <- 1:dim(x@weights)[2]
              }
              for (i in which){
                weightmatrix <- matrix(nrow=rows,ncol=cols)
                weightmatrix[xycoor]<- x@weights[,i]
                weightmatrix[xycoor2]<- x@weights[,i]
                #this line flips the matrix around so it is correct
                weightmatrix <- as.matrix(rev(as.data.frame(weightmatrix)))
                image(weightmatrix,col=terrain.colors(12),xaxt='n', yaxt='n',main=sampleNames(x)[i])
                #title(sampleNames(x)[i])
              }
            }
            )
  
  setMethod("boxplot",signature(x="PLMset"),
            function(x,...){
              grp.rma.se1.median <- apply(se(x), 1, median)
              grp.rma.rel.se1.mtx <- sweep(se(x),1,grp.rma.se1.median,FUN='/')
              boxplot(data.frame(grp.rma.rel.se1.mtx),...)
            }
            )



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
              cat("notes=",object@notes,"\n",sep="")
              cat("The model fit for each probeset was",object@model.description,"\n")
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

# A summary method, to be cleaned up better at a later date.
# 
#  setMethod("summary","PLMset",
#          function(object,which=NULL){#
#
#              if (is.null(which)){
#                which <- rownames(object@chip.coefs)
#              }
#              cur.const.coef <-  NULL
#              cur.const.se <- NULL
#              for (probeset.names in which){
#                if (dim( exprs(object)) == 0){
#                  cur.const.coef <- exprs(object)[rownames(exprs(object)) == probeset.names]
#                  cur.const.se <-  se.exprs(object)[rownames(exprs(object)) == probeset.names]
#                }
#                inds <- grep(probeset.names,rownames(object@probe.coefs))
#                cur.probe.coef <- object@probe.coefs[inds,]
#                cur.se.probe.coef <- object@se.probe.coefs[inds,]
#                cur.chip.coef <- object@chip.coefs[grep(probeset.names,rownames(object@chip.coefs)),]
#                cur.chip.se <- object@se.chip.coefs[grep(probeset.names,rownames(object@se.chip.coefs)),]#
#
#                print(cbind(c(cur.const.coef, cur.chip.coef,cur.probe.coef),c(cur.const.se, cur.chip.se,cur.se.probe.coef)))                
#              }
#            },where=where)

  if (!isGeneric("Mbox"))
    setGeneric("Mbox",function(object,...)
               standardGeneric("Mbox"))
  

  
  setMethod("Mbox",signature("PLMset"),
            function(object,...){
              medianchip <- apply(coefs(object), 1, median)
              M <- sweep(coefs(object),1,medianchip,FUN='-')
              boxplot(data.frame(M),...)
            })






