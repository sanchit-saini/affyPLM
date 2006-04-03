ReadRMAExpress <- function(filename,return.value=c("exprSet","matrix")){

  return.value <- match.arg(return.value)


  if (return.value == "matrix"){
    return(.Call("read_rmaexpress",filename,PACKAGE="affyPLM"))
  } else if (return.value == "exprSet"){
    intens <- .Call("read_rmaexpress",filename,PACKAGE="affyPLM")
    header <- .Call("read_rmaexpress_header",filename,PACKAGE="affyPLM")

    description <- new("MIAME")
    description@preprocessing$filenames <- header[[3]]
    description@preprocessing$rmaexpressversion <- paste("RMAExpress",header[[1]])

    samplenames <- sub("^/?([^/]*/)*", "", header[[3]],extended = TRUE)
    pdata <- data.frame(sample = 1:length(samplenames), row.names = samplenames)
    phenoData <- new("phenoData", pData = pdata, varLabels = list(sample = "arbitrary numbering"))

    return(new("exprSet",exprs=intens, se.exprs =array(NA, dim(intens)),annotation=cleancdfname(header[[2]], addcdf = FALSE),notes=paste("Processed using RMAExpress",header[[1]]),phenoData=phenoData,description=description))
  }
}
