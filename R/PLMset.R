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
  
  

  # getCdfInfo, using same code as in "affy" package.
  
  setMethod("getCdfInfo", signature("PLMset"), 
            function (object, how=getOption("BioC")$affy$probesloc) {
              ## "how" is a list. each element of the list must have an element
              ## tagged "what" and an element tagged "where"
              ## "what" can be "data", "package", "file" or "environment"
              ## "where" is where it can be found

              cdfname <- cleancdfname(object@cdfName)
              second.try <- FALSE

              if (debug.affy123)
                cat("Trying to get cdfenv for", cdfname, "\n")
              
              
              i <- 0
              while(i < length(how)) {
                i <- i+1
                
                if (debug.affy123)
                  cat(i, ":")
                
                what <- how[[i]]$what
                where <- how[[i]]$where

                if (debug.affy123) {
                  cat("what=", what, "where=")
                  print(where)
                }
                
                if (what == "data") {
                  ##if we can get it from data dir. otherwise load package
                  
                  if(cdfname %in% do.call("data", list(package=where))$results[, 3]) {
                    ##RI: package="affy" doesnt work it has to be package=affy
                    ##    fix if you can
                    ##LG: weird stuff with data... but I had a workaround...
                    
                    where.env <- pos.to.env(match(paste("package:", where, sep = ""), search()))
                    
                    ## check if the cdfenv is already loaded. If not load it *in* the environment
                    ## of the package (where.env)
                    if( ! exists(cdfname, where = where.env, inherits = FALSE)) {
                      path <- .path.package(where)
                      filename <- paste(cdfname, ".rda", sep="")
                      load(file.path(path, "data", filename) ,
                           envir = where.env)
                    }
                    cdfenv <- get(cdfname, envir=where.env)
                    return(cdfenv)
                  }
                  next
                }
                
                if (what == "package") {
                  loc <- .find.package(cdfname, lib.loc=where, quiet=TRUE)
                  
                  if (!second.try && identical(loc, character(0))) {
                    ## before jumping to the next option, check the possibility to
                    ## download the missing cdfenv pack
                    
                    if (how[[i]]$autoload) {
                      cat(paste("Environment",cdfname,"is not available.\n"))
                      cat("This environment contains needed probe location information.\n\n")
                      
                      cat(paste("We will try to download and install the",
                                cdfname,"package.\n\n"))
                      if (! "package:reposTools" %in% search()) {
                        on.exit(detach("package:reposTools"))
                      }
                      
                      if (! require(reposTools))
                        stop(paste("The package reposTools is required to download environments.",
                                   "Please download and install it.\n", sep="\n"))
                      
                      reposEntry <- getReposEntry(how[[i]]$repository)
                      
                      if (is.null(how[[i]]$installdir))
                        status.install <- install.packages2(cdfname, reposEntry)
                      else
                        status.install <- install.packages2(cdfname, reposEntry, how[[i]]$installdir)

                      if (length(statusList(status.install)) == 0) {
                        warning(paste("Data package", cdfname, "does not seem to exist", 
                                      "in the repository\n",
                                      how[[i]]$repository, "\n"))
                      } else {
                        ## rewind the iterator i and try again
                        i <- i-1
                        second.try <- TRUE
                      }
                    }
                    ## jump to next way to get the cdfenv
                    next
                  }
                  if (length(loc) > 1)
                    warning(paste("several packages with a matching name. Using the one at", loc[1]))
                  
                  existsnow<- .find.package(cdfname, lib.loc=where, quiet=TRUE)

                  if(!identical(loc, character(0))){
                    do.call("library", list(cdfname, lib.loc=dirname(loc[1])))
                    
                    return(get(cdfname, envir=as.environment(paste("package:", cdfname, sep=""))))
                  }
                  next
                }                
                
                if (what == "file") {
                  ##now this is an actual Affymetrix filename
                  cdfname <- paste(object@cdfName,".CDF",sep="")
                  cdf <- read.cdffile(file.path(path.expand(where), cdfname))
                  ## ---> extra paranoia <---
                  if (cdf@cdfName != object@cdfName)
                    warning(paste("The CDF file identifies as", cdf@cdfName,
                                  "while you probably want", object@cdfName))
                  ## ---> end of paranoia <---
                  return(getLocations.Cdf(cdf))
                }
                
                if (what == "environment") {
                  if(exists(object@cdfName,inherits=FALSE,where=where))
                    return(as.environment(get(object@cdfName,inherits=FALSE,envir=where)))
                  next
                }
              }
              
              warning(paste("\nWe could not find and/or install the necessary probe location information.\n",
                            "Here is a list of common problems and possible solutions:\n\n",
                            "Problem 1: You are not connected to the Internet.\n",
                            "Solution:  Try again once you connect to the Internet.\n\n",
                            "Problem 2: You do not have the necessary permissions to install packages.\n",
                            "Solution:  Ask your system administrator to install ",cleancdfname(object@cdfName), " package from:\n",
                            "           http://www.bioconductor.org/data/cdfenvs/cdfenvs.html\n\n",
                            "Problem 3: Necessary package not available from Bioconductor.\n",
                            "Solution:  Use makecdfenv package to create environment from CDF file.\n",
                            "           See the dealing_with_cdfenvs vignette for more details\n\n",
                            "NOTE: Once you install ",cleancdfname(object@cdfName)," you should not have this problem again.\n",
                            sep=""))
              warning(paste("To let you proceed for now, a dummy cdfenv", cleancdfname(object@cdfName),
                            "will be created..."))
              if (exists(object@cdfName, envir=.GlobalEnv)) {
                stop("Could not create dummy environment. Giving up.")
              }
              assign(object@cdfName,
                     new.env(parent=.GlobalEnv),
                     envir=.GlobalEnv)
              warning(paste("IMPORTANT: Depending on your settings, you might have to delete the object",
                            object@cdfName, " after you install right the package !\n(command 'rm(",
                            object@cdfName,")' )\n"))
              return(get(object@cdfName, envir=.GlobalEnv))
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






