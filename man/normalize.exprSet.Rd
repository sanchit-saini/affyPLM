\name{normalize.exprSet}
\alias{normalize.exprSet.quantiles}
\alias{normalize.exprSet.loess}
\alias{normalize.exprSet.contrasts}
\alias{normalize.exprSet.qspline}
\alias{normalize.exprSet.invariantset}
\alias{normalize.exprSet.scaling}
\alias{normalize.exprSet.methods}

\title{Normalization applied to exprSets}
\description{
  Allows the user to apply normalization routines to exprSets.
}
\usage{
  normalize.exprSet.quantiles(eset,transfn=c("none","log","antilog"))
  normalize.exprSet.loess(eset,transfn=c("none","log","antilog"),...)
  normalize.exprSet.contrasts(eset, span = 2/3, choose.subset = TRUE, subset.size = 5000, verbose = TRUE, family = "symmetric",transfn=c("none","log","antilog"))
  normalize.exprSet.qspline(eset,transfn=c("none","log","antilog"),...)
  normalize.exprSet.invariantset(eset,prd.td = c(0.003, 0.007),verbose=FALSE,transfn=c("none","log","antilog"))
normalize.exprSet.scaling(eset,trim=0.02,baseline=-1,transfn=c("none","log","antilog"))
}
\arguments{
  \item{eset}{An \code{\link[Biobase:exprSet-class]{exprSet}}}
  \item{span}{parameter to be passed to the function \code{\link[modreg]{loess}}.}
  \item{choose.subset}{}
  \item{subset.size}{}
  \item{verbose}{verbosity flag}
  \item{family}{parameter to be passed to the function
    \code{\link[modreg]{loess}}.}
  \item{prd.td}{cutoff parameter (details in the bibliographic
    reference)}
  \item{trim}{How much to trim from the top and bottom before computing
    the mean when using the scaling normalization}
  \item{baseline}{Index of array to use as baseline, negative values
    (-1,-2,-3,-4) control different baseline selection methods}
  \item{transfn}{Transform the exprSet before normalizing. Useful when
    dealing with expression values that are log-scale}
  \item{...}{Additional parameters that may be passed to the
    normalization routine}
  
}
\details{
  This function carries out normalization of expression values. In
  general you should either normalize at the probe level or at the
  expression value level, not both.

  Typing \code{normalize.exprSet.methods} should give you a list of
  methods that you may use. note that you can also use the
  \code{normalize} function on exprSets. Use \code{method} to select the
  normalization method.
}

\value{
  A normalized \code{\link[Biobase:exprSet-class]{exprSet}}.
}

\examples{
data(affybatch.example)
eset <- rma(affybatch.example,normalize=FALSE,background=FALSE)
normalize(eset)
}


\author{Ben Bolstad, \email{bolstad@stat.berkeley.edu}}

\seealso{\code{\link[affy]{normalize}}}

\keyword{manip}
