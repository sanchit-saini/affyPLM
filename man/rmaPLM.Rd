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

\references{Under Preparation}
\seealso{\code{\link[affy]{expresso}},
  \code{\link[affy]{rma}}, \code{\link{threestep}},\code{\link{fitPLM}}, \code{\link{threestepPLM}}
}


\examples{
# A larger example testing weight image function
data(Dilution)
Pset <- rmaPLM(Dilution,output.param=list(weights=TRUE))
image(Pset)
}

\keyword{manip}
