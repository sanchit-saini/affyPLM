###
###
### May  1, 2006 - added support for user gzipped files (should read successfully)
### Feb 27, 2008 - removed support for exprSet objects
###



ReadRMAExpress <- function(filename,return.value=c("ExpressionSet","matrix")){

  return.value <- match.arg(return.value)

  file.type <- .Call("check_rmaexpress_format",filename,PACKAGE="affyPLM")

  if (file.type =="Uncompressed"){
    if (return.value == "matrix"){
      return(.Call("read_rmaexpress",filename,PACKAGE="affyPLM"))
    }
  } else if (file.type == "Compressed"){
     if (return.value == "matrix"){
      return(.Call("gz_read_rmaexpress",filename,PACKAGE="affyPLM"))
    } else if (return.value == "ExpressionSet"){
      intens <- .Call("gz_read_rmaexpress",filename,PACKAGE="affyPLM")
      header <- .Call("gz_read_rmaexpress_header",filename,PACKAGE="affyPLM")
      
      description <- new("MIAME")
      preproc(description)$filenames <- header[[3]]
      preproc(description)$rmaexpressversion <- paste("RMAExpress",header[[1]])
      
      samplenames <- sub("^/?([^/]*/)*", "", header[[3]],extended = TRUE)
      pdata <- data.frame(sample = 1:length(samplenames), row.names = samplenames)
      phenoData <- new("AnnotatedDataFrame", data = pdata, 
          varMetadata = data.frame(labelDescription= "arbitrary numbering"))
      
      return(new("ExpressionSet", 
             exprs = intens, 
             se.exprs = array(NA, dim(intens)),
             annotation = cleancdfname(header[[2]], addcdf = FALSE),
             ##FIXME: remove # below when notes is fixed
             #notes=paste("Processed using RMAExpress",header[[1]]),
             phenoData = phenoData,
             experimentData = description))
    }
   } else {
     stop("Don't recognize the format of this file.\n")
   }
}
