\name{normalize.quantiles.probeset}
\alias{normalize.AffyBatch.quantiles.probeset}
\title{Quantile Normalization applied to probesets}
\description{
  Using a normalization based upon quantiles, this function
  normalizes a matrix of probe level intensities.
}
\usage{
  normalize.AffyBatch.quantiles.probeset(abatch,type=c("separate","pmonly","mmonly","together"),use.median=FALSE,use.log=TRUE)
}
\arguments{
  \item{abatch}{An \code{\link[affy:AffyBatch-class]{AffyBatch}}}
  \item{type}{how should MM and PM values be handled}
  \item{use.median}{use median rather than mean}
  \item{use.log}{take logarithms, then normalize}
}
\details{
  This function applies the \link[affy:normalize.quantiles]{quantile} method in a probeset specific manner.

  In particular a probeset summary is normalized using the quantile
  method and then the probes adjusted accordingly.

}

\value{
  A normalized \code{AffyBatch}.
}
\references{
  Bolstad, B (2001) \emph{Probe Level Quantile Normalization of High Density
    Oligonucleotide Array Data}. Unpublished manuscript
  \url{http://oz.berkeley.edu/~bolstad/stuff/qnorm.pdf}

  Bolstad, B. M., Irizarry R. A., Astrand, M, and Speed, T. P. (2003)
  \emph{A Comparison of Normalization Methods for High Density
    Oligonucleotide Array Data Based on Bias and Variance.}
   Bioinformatics 19(2) ,pp 185-193. \url{http://www.stat.berkeley.edu/~bolstad/normalize/normalize.html}
  }

\author{Ben Bolstad, \email{bmb@bmbolstad.com}}

\seealso{\code{\link[affy]{normalize}}, \code{\link[affy:normalize.quantiles]{normalize.quantiles}}}

\keyword{manip}
