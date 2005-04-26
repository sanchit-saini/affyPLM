\name{bg.correct.LESN}
\alias{bg.correct.LESN}
\title{LESN - Low End Signal is Noise Background corrections}
\description{
  This function background corrects PM probe data using LESN - Low End
  Signal is Noise concepts.
}
\usage{
bg.correct.LESN(object,method = 2,baseline = 0.25, theta=4)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{object}{an \code{\link[affy]{AffyBatch}}}
  \item{method}{an integer code specifying which method to use}
  \item{baseline}{A baseline value to use}
  \item{theta}{A parameter used in the background correction process}
}
\details{
  This method will be more formally documented at a later date.

  The basic concept is to consider that the lowest end of intensites is
  most likely just noise (and should be heavily corrected) and the
  highest end signals are most likely signal and should have little
  adjustment. Low end signals are made much smaller while high end
  signals get less adjustment relative adjustment.
}
\value{
An \code{\link[affy]{AffyBatch}}
}



\references{Bolstad, BM (2004) \emph{Low Level Analysis of High-density
    Oligonucleotide Array Data: Background, Normalization and
    Summarization}. PhD Dissertation. University of California, Berkeley.}


\author{Ben Bolstad \email{bolstad@stat.berkeley.edu}}

\examples{
data(affybatch.example)
affybatch.example.bgcorrect <- bg.correct.LESN(affybatch.example)
}
\keyword{manip}

