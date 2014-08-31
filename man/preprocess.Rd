\name{preprocess}
\alias{preprocess}
\title{Background correct and Normalize}
\description{
  This function pre-processes an \code{\link[affy:AffyBatch-class]{AffyBatch}}.
}
\usage{
preprocess(object, subset=NULL, normalize=TRUE, background=TRUE,
           background.method="RMA.2", normalize.method="quantile",
           background.param=list(), normalize.param=list(),
           verbosity.level=0)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{object}{an \code{\link[affy:AffyBatch-class]{AffyBatch}}}
  \item{subset}{a vector with the names of probesets to be used. If NULL
  then all probesets are used.}
  \item{normalize}{logical value. If \code{TRUE} normalize data using
    quantile normalization}
  \item{background}{logical value. If \code{TRUE} background correct
    using RMA background correction}
  \item{background.method}{name of background method to use.}
  \item{normalize.method}{name of normalization method to use.}
  \item{background.param}{list of parameters for background correction methods}
  \item{normalize.param}{list of parameters for normalization methods}
  \item{verbosity.level}{An integer specifying how much to print
    out. Higher values indicate more verbose. A value of 0 will print nothing}
	    
}
\details{
  This function carries out background correction and normalization
  pre-processing steps. It does not summarize to produce gene expression
  measures. All the same pre-processing methods supplied by
  \code{\link{threestep}} are supported by this function.
}
\value{
An \code{\link[affy:AffyBatch-class]{AffyBatch}}
}

\author{Ben Bolstad \email{bmb@bmbolstad.com}}

\references{Bolstad, BM (2004) \emph{Low Level Analysis of High-density
    Oligonucleotide Array Data: Background, Normalization and
    Summarization}. PhD Dissertation. University of California, Berkeley.}

\seealso{\code{\link[affy]{expresso}}, \code{\link[affy]{rma}}}

\examples{
if (require(affydata)) {
  data(Dilution)

  # should be equivalent to the bg and norm of rma()
  abatch.preprocessed <- preprocess(Dilution)
}
}
\keyword{manip}

