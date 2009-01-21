\name{PLMset2exprSet}
\alias{PLMset2exprSet}
\alias{pset2eset}
\title{Convert a PLMset to an ExpressionSet}
\description{
  This function converts a PLMset to an ExpressionSet. This is often useful
  since many Bioconductor functions operate on ExpressionSet objects.
}
\usage{
PLMset2exprSet(pset)
pset2eset(pset)
}
\arguments{
  \item{pset}{The \code{\link{PLMset}} to convert to \code{\link[Biobase]{ExpressionSet}}.}
}
\details{
  These functions convert PLMset objects to \code{ExpressionSet} objects.
  This is often useful since many Bioconductor functions operate on
  \code{ExpressionSet} objects. Note that the function \code{pset2eset}
  is a wrapper for \code{PLMset2exprSet}.
}
\value{
  returns a \code{\link[Biobase]{ExpressionSet}}
}
\author{Ben Bolstad \email{bmb@bmbolstad.com}}
\seealso{\code{\link[Biobase]{ExpressionSet}}}
\examples{
if (require(affydata)) {
  data(Dilution)
  Pset <- fitPLM(Dilution)
  eset <- pset2eset(Pset)
}
}
\keyword{manip}
