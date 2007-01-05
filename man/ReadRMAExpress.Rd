\name{ReadRMAExpress}
\alias{ReadRMAExpress}
\title{Read RMAExpress computed expression values}
\description{Read RMAExpress computed binary output files into a matrix
  or ExpressionSet
}
\usage{
ReadRMAExpress(filename,return.value=c("ExpressionSet","matrix"))
}
\arguments{
  \item{filename}{The name of the file containing RMAExpress output to
    be read in}
  \item{return.value}{should a \code{\link{matrix}} or an \code{\link[Biobase:ExpressionSet-class]{ExpressionSet}} be returned}
}
\value{
  returns an \code{\link[Biobase:ExpressionSet-class]{ExpressionSet}}
}
\references{
  \url{http://rmaexpress.bmbolstad.com}
}
\author{Ben Bolstad \email{bmb@bmbolstad.com}}
\keyword{manip}
