normalize.AffyBatch.quantiles.probeset <- function(abatch,type=c("separate","pmonly","mmonly","together"),use.median=FALSE,use.log=TRUE) {

  type <- match.arg(type)


  rows <- length(probeNames(abatch))
  cols <- length(abatch)

  ngenes <- length(geneNames(abatch))

  if ((type == "pmonly")|(type == "separate")){
    pms <- unlist(pmindex(abatch))
    noNA <- apply(intensity(abatch)[pms,],1,function(x) all(!is.na(x)))
    pms <- pms[noNA]
    intensity(abatch)[pms,] <-
    .C("qnorm_probeset_R",as.double(intensity(abatch)[pms,
    ]),as.integer(rows),
    as.integer(cols),as.integer(ngenes),as.character(probeNames(abatch)),as.integer(use.median),as.integer(use.log),
    PACKAGE="affyPLM")[[1]]
  }

  if ((type == "mmonly")|(type == "separate")){
    mms <- unlist(mmindex(abatch))
    noNA <- apply(intensity(abatch)[mms, , drop = FALSE],
                  1, function(x) all(!is.na(x)))
    mms <- mms[noNA]
    intensity(abatch)[mms, ] <-
    .C("qnorm_probeset_R",as.double(intensity(abatch)[mms,]),as.integer(rows),
    as.integer(cols),as.integer(ngenes),as.character(probeNames(abatch)),as.integer(use.median),as.integer(use.log),
    PACKAGE="affyPLM")[[1]]
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


