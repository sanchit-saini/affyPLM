\name{fitPLM}
\alias{fitPLM}
\title{Fit a Probe Level Model to Affymetrix Genechip Data.}
\description{
  This function converts an \code{\link[affy]{AffyBatch}} into an
  \code{\link{PLMset}} by fitting a specified robust linear model to the
  probe level data}
}
\usage{
fitPLM(object,model=PM ~ -1 + probes +samples,
variable.type=c(default="factor"),
constraint.type=c(default="contr.treatment"),
background=TRUE, normalize=TRUE, background.method = "RMA.2",normalize.method = "quantile",se.type=4,psi.type="Huber",psi.k=NULL,background.param=list(),normalize.param=list())
}

\arguments{
  \item{object}{an \code{\link[affy]{AffyBatch}}}
  \item{model}{A formula describing the model to fit. This is slightly
    different from the standard method of specifying formulae in R. Read
    the description below}
  \item{variable.type}{a way to specify whether variables in the model
    are factors or standard variables}
  \item{constraint.type}{should factor variables sum to zero or have
    first variable set to zero (endpoint constraint)}
  \item{normalize}{logical value. If \code{TRUE} normalize data using
    quantile normalization}
  \item{background}{logical value. If \code{TRUE} background correct
    using RMA background correction}
  \item{background.method}{name of background method to use.}
  \item{normalize.method}{name of normalization method to use.}
  \item{se.type}{integer value indicating which of the four standard
    error methods to use}
  \item{psi.type}{A name of a weight function to use in M estimation}
  \item{psi.k}{A tuning parameter for \code{psi.type}}
  \item{background.param}{A list of parameters for background routines}
  \item{normalize.param}{A list of parameters for normalization routines}
  
}
\details{
  This function fits robust Probe Level Linear models to all the probesets in
  an \code{\link[affy]{AffyBatch}}.
}
\value{
  An \code{\link{PLMset}}
}

\author{Ben Bolstad \email{bolstad@stat.berkeley.edu}}

\references{Under Preparation}
\seealso{\code{\link[affy]{express}},\code{\link[affy]{expresso}},
  \code{\link[affy]{rma}}, \code{\link{threestep}}}

\examples{
data(affybatch.example)
Pset <- fitPLM(affybatch.example,model=PM ~ -1 + probes + samples)
se(Pset)[1:5,]


# A larger example testing weight image function
data(Dilution)
Pset <- fitPLM(Dilution,model=PM ~ -1 + probes + samples)
image(Pset)
boxplot(Pset)

#now lets try a wider class of models
Pset <- fitPLM(Dilution,model=PM ~ -1 + probes +liver,normalize=FALSE,background=FALSE)
coefs(Pset)[1:10,]

Pset <- fitPLM(Dilution,model=PM ~ -1 + probes + liver + scanner,normalize=FALSE,background=FALSE)
coefs(Pset)[1:10,]

#try liver as a covariate
logliver <- log2(c(20,20,10,10))
Pset <- fitPLM(Dilution,model=PM~-1+probes+logliver+scanner,normalize=FALSE,background=FALSE,variable.type=c(logliver="covariate"))
coefs(Pset)[1:10,]


#try a different se.type
Pset <- fitPLM(Dilution,model=PM~-1+probes+scanner,normalize=FALSE,background=FALSE,se.type=2)
se(Pset)[1:10,]
Pset <- fitPLM(Dilution,model=PM~-1+probes+scanner,normalize=FALSE,background=FALSE,se.type=3)
se(Pset)[1:10,]
}

\keyword{manip}

