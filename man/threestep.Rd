\name{threestep}
\alias{threestep}
\title{Three Step expression measures}
\description{
  This function converts an \code{\link[affy]{AffyBatch}} into an
  \code{\link[Biobase]{exprSet}} using
 a three step expression measure.
}
\usage{
threestep(object,subset=NULL, normalize=TRUE,background=TRUE,background.method="RMA.2",normalize.method="quantile",summary.method="median.polish",background.param = list(),normalize.param=list(),summary.param=list(),verbosity.level=0)
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
  \item{summary.method}{name of summary method to use.}
  \item{background.param}{list of parameters for background correction methods}
  \item{normalize.param}{list of parameters for normalization methods}
  \item{summary.param}{list of parameters for summary methods}
  \item{verbosity.level}{An integer specifying how much to print
    out. Higher values indicate more verbose. A value of 0 will print nothing}
	    
}
\details{
  This function computes the expression measure using threestep
  methods. Greater details can be found in a vignette.}
}
\value{
An \code{\link[Biobase]{exprSet}}
}

\author{Ben Bolstad \email{bmb@bmbolstad.com}}

\references{Bolstad, BM (2004) \emph{Low Level Analysis of High-density
    Oligonucleotide Array Data: Background, Normalization and
    Summarization}. PhD Dissertation. University of California, Berkeley.}

\seealso{\code{\link[affy]{expresso}}, \code{\link[affy]{rma}}}

\examples{
data(affybatch.example)

# should be equivalent to rma()
eset <- threestep(affybatch.example)

# Using Tukey Biweight summarization
eset <- threestep(affybatch.example,summary.method="tukey.biweight")

# Using Average Log2 summarization
eset <- threestep(affybatch.example,summary.method="average.log")

# Using IdealMismatch background and Tukey Biweight and no normalization.
eset <- threestep(affybatch.example,normalize=FALSE,background.method="IdealMM",summary.method="tukey.biweight")

# Using average.log summarization and no background or normalization.
eset <- threestep(affybatch.example,background=FALSE,normalize=FALSE,background.method="IdealMM",summary.method="tukey.biweight")

# Use threestep methodology with the rlm model fit
eset <- threestep(affybatch.example,summary.method="rlm")

# Use threestep methodology with the log of the average
eset <- threestep(affybatch.example,summary.method="log.average")

# Use threestep methodology with log 2nd largest method
eset <- threestep(affybatch.example,summary.method="log.2nd.largest")

eset <- threestep(affybatch.example,background.method="LESN2")
}
\keyword{manip}

