\name{rmaPLM}
\alias{rmaPLM}
\title{Fit a RMA to Affymetrix Genechip Data as a PLMset}
\description{
  This function converts an \code{\link[affy]{AffyBatch}} into an
  \code{\link{PLMset}} by fitting a multichip model. In particular we
  concentrate on the RMA model.}
}
\usage{
rmaPLM(object,subset=NULL,normalize=TRUE,background=TRUE,background.method="RMA.2",normalize.method="quantile",background.param = list(),normalize.param=list(),output.param=list(),model.param=list())
}

\arguments{
  \item{object}{an \code{\link[affy]{AffyBatch}}}
  \item{subset}{a vector with the names of probesets to be used. If NULL
  then all probesets are used.}
  \item{normalize}{logical value. If \code{TRUE} normalize data using
    quantile normalization}
  \item{background}{logical value. If \code{TRUE} background correct
    using RMA background correction}
  \item{background.method}{name of background method to use.}
  \item{normalize.method}{name of normalization method to use.}
  \item{background.param}{A list of parameters for background routines}
  \item{normalize.param}{A list of parameters for normalization
    routines}
  \item{output.param}{A list of parameters controlling optional output
    from the routine.}
  \item{model.param}{A list of parameters controlling model procedure}
}
\details{
  This function fits the RMA as a Probe Level Linear models to all the probesets in
  an \code{\link[affy]{AffyBatch}}.
}
\value{
  An \code{\link{PLMset}}
}

\author{Ben Bolstad \email{bolstad@stat.berkeley.edu}}

\references{Bolstad, BM (2004) \emph{Low Level Analysis of High-density
    Oligonucleotide Array Data: Background, Normalization and
    Summarization}. PhD Dissertation. University of California,
  \cr \cr Irizarry RA, Bolstad BM, Collin F, Cope LM, Hobbs B and Speed
  TP (2003) \emph{Summaries of Affymetrix GeneChip probe level data}
  Nucleic Acids Research 31(4):e15
  \cr \cr Bolstad, BM, Irizarry RA, Astrand, M, and Speed, TP (2003)
  \emph{A Comparison of Normalization Methods for High Density
    Oligonucleotide Array Data Based on Bias and Variance.}
  Bioinformatics 19(2):185-193}

\seealso{\code{\link[affy]{expresso}},
  \code{\link[affy]{rma}}, \code{\link{threestep}},\code{\link{fitPLM}}, \code{\link{threestepPLM}}
}


\examples{
# A larger example testing weight image function
data(Dilution)
\dontrun{Pset <- rmaPLM(Dilution,output.param=list(weights=TRUE))}
\dontrun{image(Pset)}
}

\keyword{manip}
