\name{normalize.scaling}
\alias{normalize.scaling}
\alias{normalize.AffyBatch.scaling}

\title{Scaling normalization}
\description{
  Allows the user to apply scaling normalization.
}
\usage{
normalize.scaling(X,trim=0.02, baseline=-1, log.scalefactors=FALSE)
normalize.AffyBatch.scaling(abatch,
    type=c("together","pmonly","mmonly","separate"),
    trim=0.02, baseline=-1, log.scalefactors=FALSE) 
}
\arguments{
  \item{X}{A matrix. The columns of which are to be normalized.}
  \item{abatch}{An \code{\link[affy:AffyBatch-class]{AffyBatch}}}
  \item{type}{A parameter controlling how normalization is applied to
    the Affybatch.}
  \item{trim}{How much to trim from the top and bottom before computing
    the mean when using the scaling normalization.}
  \item{baseline}{Index of array to use as baseline, negative values
    (-1,-2,-3,-4) control different baseline selection methods.}
  \item{log.scalefactors}{Compute the scale factors based on log2
    transformed data.}
  
}
\details{
  These function carries out scaling normalization of expression values.
}
\value{
  A normalized \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}}.
}
\examples{
if (require(affydata)) {
  data(Dilution)
  normalize.AffyBatch.scaling(Dilution)
}
}
\author{Ben Bolstad, \email{bmb@bmbolstad.com}}
\keyword{manip}
