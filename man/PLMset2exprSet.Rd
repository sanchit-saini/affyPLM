\name{PLMset2exprSet}
\alias{PLMset2exprSet}
\alias{pset2eset}
\title{Convert a PLMset to an exprSet}
\description{
  This function converts a PLMset to an exprSet. This is often useful
  since many Bioconductor functions operate on exprSet objects.
}
\usage{
PLMset2exprSet(pset)
pset2eset(pset)
}

\arguments{
  \item{pset}{The \code{\link{PLMset}} to convert to \code{\link[Biobase]{exprSet}}}
}
\details{
  These functions convert  PLMset objects to exprSet objects. This is often useful
  since many Bioconductor functions operate on exprSet objects. Note
  that the finction \code{pset2eset} is a wrapper for \code{PLMset2exprSet}.
}

\value{
  returns a \code{\link[Biobase]{exprSet}}
  }



\author{Ben Bolstad \email{bolstad@stat.berkeley.edu}}

\seealso{\code{\link[Biobase]{exprSet}}}

\examples{
data(affybatch.example)
Pset <- fitPLM(affybatch.example)
eset <- pset2eset(Pset)
}
\keyword{manip}
