##################################################
##
## file: normalize.exprSet.R
##
## aim: normalization routines as applied to
##      exprSets (ie normalize expression values)
##
## Methods implemented
## quantiles  (quantile normalization)
## loess  (cyclic loess)
## contrasts (contrast loess)
## qspline (quantile spline method)
## invariantset (Similar to method used in dChip)
##
## Most of this is just wrappers around the approrpiate
## routine from the affy package, just adapted to deal
## with exprSet
##
##
## created: Aug 22, 2003
## written by: B. M. Bolstad <bolstad@stat.berkeley.edu>
##
## History
## Aug 22, 2003 - Initial version
## Aug 23, 2003 - added nomalize.exprSet.scaling
## Aug 25, 2003 - Added parameters to control whether
##                a log/antilog transformation should be applied
##                this is most useful with RMA expression estimate
##                that is usually given in log2 scale
##                note that while normalization is carried
##                out in transformed scale, returned values will be in
##                original scale
##
##################################################

normalize.exprSet.quantiles <- function(eset,transfn=c("none","log","antilog")){
  transfn <- match.arg(transfn)
  
  col.names <- colnames(exprs(eset))
  row.names <- rownames(exprs(eset))
  if (transfn == "none"){
    exprs(eset) <- normalize.quantiles(exprs(eset))
  } else if (transfn == "antilog"){
    exprs(eset) <- log2(normalize.quantiles(2^exprs(eset)))
  } else {
    exprs(eset) <- 2^(normalize.quantiles(log2(exprs(eset))))
  }
  colnames(exprs(eset)) <- col.names
  rownames(exprs(eset)) <- row.names
  return(eset)
}


normalize.exprSet.loess <- function(eset,transfn=c("none","log","antilog"),...){
  transfn <- match.arg(transfn)

  if (transfn == "none"){
    exprs(eset) <- normalize.loess(exprs(eset),...)
  } else if (transfn == "antilog"){
    exprs(eset) <- log2(normalize.loess(2^exprs(eset),...))
  } else {
    exprs(eset) <- 2^(normalize.loess(log2(exprs(eset)),...))
  }
  return(eset)
}


normalize.exprSet.contrasts <- function(eset, span = 2/3, choose.subset = TRUE, subset.size = 5000, verbose = TRUE, family = "symmetric",transfn=c("none","log","antilog")){

  transfn <- match.arg(transfn)
  
  col.names <- colnames(exprs(eset))
  row.names <- rownames(exprs(eset))
  alldata <- exprs(eset)

  if (transfn=="antilog"){
    alldata <- 2^(alldata)
  }
  if (transfn=="log"){
    alldata <- log2(alldata)
  }
  
  if (choose.subset)
    subset1 <- maffy.subset(alldata, verbose = verbose, subset.size = subset.size)$subset
  else subset1 <- sample(1:dim(alldata)[1], subset.size)
  aux <- maffy.normalize(alldata, subset = subset1, verbose = verbose,
                         span = span, family = family)

  if (transfn=="antilog"){
    alldata <- log2(alldata)
  }

  if (transfn=="log"){
    alldata <- 2^(alldata)
  }
  
  exprs(eset) <- aux
  colnames(exprs(eset)) <- col.names
  rownames(exprs(eset)) <- row.names
  return(eset)
}

normalize.exprSet.qspline <- function(eset,transfn=c("none","log","antilog"),...){

  transfn <- match.arg(transfn)
  col.names <- colnames(exprs(eset))
  row.names <- rownames(exprs(eset))
  if (transfn == "none"){
    exprs(eset) <- normalize.qspline(exprs(eset),...)
  } else if (transfn == "antilog"){
     exprs(eset) <- log2(normalize.qspline(2^exprs(eset),...))
  } else {
    exprs(eset) <- 2^(normalize.qspline(log2(exprs(eset)),...))
  }
  colnames(exprs(eset)) <- col.names
  rownames(exprs(eset)) <- row.names
  return(eset)
}


normalize.exprSet.invariantset <- function(eset,prd.td = c(0.003, 0.007), verbose = FALSE,transfn=c("none","log","antilog")){
  
  transfn <- match.arg(transfn)

  require(modreg, quietly = TRUE)
  nc <- length(sampleNames(eset))
  m  <- vector("numeric", length = nc)


  alldata <- exprs(eset)
  if (transfn == "log"){
    alldata <- log2(alldata)
  }

  if (transfn == "antilog"){
    alldata <- 2^(alldata)
  }
   
  for (i in 1:nc) m[i] <- mean(alldata[, i])
  refindex <- trunc(median(rank(m)))
  if (verbose)
    cat("Data from", sampleNames(eset)[refindex], "used as baseline.\n")
  
  for (i in (1:nc)[-refindex]) {
    if (verbose)
      cat("normalizing array", sampleNames(eset)[i], "...")
    tmp <- normalize.invariantset(alldata[, i],
                                  alldata[, refindex], prd.td)
    tmp <- as.numeric(approx(tmp$n.curve$y, tmp$n.curve$x,
                             xout = alldata[, i], rule = 2)$y)
    alldata[, i] <- tmp
    if (verbose)
      cat("done.\n")
  }

  if (transfn == "log"){
    alldata <- 2^(alldata)
  }

  if (transfn == "antilog"){
    alldata <- log2(alldata)
  }

  exprs(eset) <- alldata
  return(eset)
}


normalize.exprSet.scaling <- function(eset,trim=0.02,baseline=-1,transfn=c("none","log","antilog")){

  transfn <- match.arg(transfn)
    
  col.names <- colnames(exprs(eset))
  row.names <- rownames(exprs(eset))
  if (transfn == "none"){
    exprs(eset) <- normalize.scaling(exprs(eset),trim,baseline)
  } else if (transfn =="antilog"){
    exprs(eset) <- log2(normalize.scaling(2^exprs(eset)))
  } else {
    exprs(eset) <- 2^(normalize.scaling(log2(exprs(eset))))
  }
  colnames(exprs(eset)) <- col.names
  rownames(exprs(eset)) <- row.names
  return(eset)
}


