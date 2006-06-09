normalize.AffyBatch.quantiles.probeset <- function(abatch,type=c("separate","pmonly","mmonly","together"),use.median=FALSE,use.log=TRUE) {

  type <- match.arg(type)
   
  
  rows <- length(probeNames(abatch))
  cols <- length(abatch)
  
  ngenes <- length(geneNames(abatch))

  if ((type == "pmonly")|(type == "separate")){
    pms <- unlist(pmindex(abatch))
    noNA <- apply(intensity(abatch)[pms,],1,function(x) all(!is.na(x)))
    pms <- pms[noNA]
    intensity(abatch)[pms,] <- .C("qnorm_probeset_R",as.double(intensity(abatch)[pms, ]),as.integer(rows), as.integer(cols),as.integer(ngenes),as.character(probeNames(abatch)),as.integer(use.median),as.integer(use.log),PACKAGE="affyPLM")[[1]]
  }

  if ((type == "mmonly")|(type == "separate")){
    mms <- unlist(mmindex(abatch))
    noNA <- apply(intensity(abatch)[mms, , drop = FALSE],
                  1, function(x) all(!is.na(x)))
    mms <- mms[noNA]
    intensity(abatch)[mms, ] <- .C("qnorm_probeset_R",as.double(intensity(abatch)[mms,]),as.integer(rows), as.integer(cols),as.integer(ngenes),as.character(probeNames(abatch)),as.integer(use.median),as.integer(use.log),PACKAGE="affyPLM")[[1]]
  }
  
  if (type == "together"){
    cat("together not current supported in quantiles.probeset\n")
  }
  

  ##this is MIAME we need to decide how to do this properly.
  ##for (i in 1:length(abatch)) {
  ##  history(abatch)[[i]]$name <- "normalized by quantiles"
  ##}
 
  return(abatch)
}






normalize.quantiles.determine.target <- function(x,target.length=NULL){

  if (!is.matrix(x)){
    stop("This function expects supplied argument to be matrix")
  }
  if (!is.numeric(x)){
    stop("Supplied argument should be a numeric matrix")
  }
  rows <- dim(x)[1]
  cols <- dim(x)[2]

  if (is.integer(x)){
    x <- matrix(as.double(x), rows, cols)
  }
  
  if (is.null(target.length)){
    target.length <- rows
  }
  
  if (target.length <= 0){
    stop("Need positive length for target.length")
  }


  return(.Call("R_qnorm_determine_target",x,target.length,PACKAGE="affyPLM"))
  

}



normalize.quantiles.use.target <- function(x,target,copy=TRUE){
  return(.Call("R_qnorm_using_target",x,target,copy,PACKAGE="affyPLM"))
}



normalize.quantiles.in.blocks <- function(x,blocks,copy=TRUE){

  rows <- dim(x)[1]
  cols <- dim(x)[2]

  if (rows != length(blocks)){
    stop("blocks is not vector of correct length")
  }

  if (is.factor(blocks)){
    blocks <- as.integer(blocks)
  }

  if (!is.numeric(blocks)){
    stop("non-numeric vector used for blocks")
  }


  return(.Call("R_qnorm_within_blocks",x,blocks,copy,PACKAGE="affyPLM"))



}






normalize.AffyBatch.quantiles.chromosome <- function(abatch,type=c("separate","pmonly","mmonly","together")){

  type <- match.arg(type)
  
  
  which.annot <- annotation(abatch)
  require(which.annot, character.only=TRUE)
  
  chromo.name <- paste(which.annot,"CHR",sep="")
  
  chromos <- do.call("c",as.list(eval(as.name(chromo.name))))
  
  which.chromos <- as.factor(chromos[probeNames(abatch)])
  
  
  if ((type == "pmonly")|(type == "separate")){
    pms <- unlist(pmindex(abatch))
    noNA <- apply(intensity(abatch)[pms,],1,function(x) all(!is.na(x)))
    pms <- pms[noNA]
    intensity(abatch)[pms,] <- normalize.quantiles.in.blocks(intensity(abatch)[pms, ],which.chromos,copy=FALSE)   
  }
  if ((type == "mmonly")|(type == "separate")){
    mms <- unlist(mmindex(abatch))
    noNA <- apply(intensity(abatch)[mms, , drop = FALSE],
                  1, function(x) all(!is.na(x)))
    mms <- mms[noNA]
    intensity(abatch)[pms,] <- normalize.quantiles.in.blocks(intensity(abatch)[pms, ],which.chromos,copy=FALSE) 
  }
  if (type == "together") {
    pms <- unlist(indexProbes(abatch, "both"))
    intensity(abatch)[pms, ] <- normalize.quantiles(intensity(abatch)[pms,, drop = FALSE],rep(which.chromos,2),copy = FALSE)
  }
  return(abatch)
}

