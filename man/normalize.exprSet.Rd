\name{normalize.ExpressionSet}
\alias{normalize.ExpressionSet.quantiles}
\alias{normalize.ExpressionSet.loess}
\alias{normalize.ExpressionSet.contrasts}
\alias{normalize.ExpressionSet.qspline}
\alias{normalize.ExpressionSet.invariantset}
\alias{normalize.ExpressionSet.scaling}
\alias{normalize.ExpressionSet.methods}

\title{Normalization applied to ExpressionSets}
\description{
  Allows the user to apply normalization routines to ExpressionSets.
}
\usage{
  normalize.ExpressionSet.quantiles(eset, transfn=c("none","log","antilog"))
  normalize.ExpressionSet.loess(eset, transfn=c("none","log","antilog"),...)
  normalize.ExpressionSet.contrasts(eset, span = 2/3,
      choose.subset=TRUE, subset.size=5000, verbose=TRUE, family="symmetric",
      transfn=c("none","log","antilog")) 
  normalize.ExpressionSet.qspline(eset, transfn=c("none","log","antilog"),...)
  normalize.ExpressionSet.invariantset(eset,prd.td=c(0.003, 0.007),
      verbose=FALSE, transfn=c("none","log","antilog"),
      baseline.type=c("mean","median","pseudo-mean","pseudo-median")) 
  normalize.ExpressionSet.scaling(eset, trim=0.02, baseline=-1,
      transfn=c("none","log","antilog"))
}
\arguments{
  \item{eset}{An \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}}}
  \item{span}{parameter to be passed to the function
    \code{\link[stats]{loess}}.} 
  \item{choose.subset}{use a subset of values to establish the
    normalization relationship}
  \item{subset.size}{number to use for subset}
  \item{verbose}{verbosity flag}
  \item{family}{parameter to be passed to the function
    \code{\link[stats]{loess}}.}
  \item{prd.td}{cutoff parameter (details in the bibliographic
    reference)}
  \item{trim}{How much to trim from the top and bottom before computing
    the mean when using the scaling normalization}
  \item{baseline}{Index of array to use as baseline, negative values
    (-1,-2,-3,-4) control different baseline selection methods}
  \item{transfn}{Transform the ExpressionSet before normalizing. Useful when
    dealing with expression values that are log-scale}
  \item{baseline.type}{A method of selecting the baseline array}
  \item{...}{Additional parameters that may be passed to the
    normalization routine}
  
}
\details{
  This function carries out normalization of expression values. In
  general you should either normalize at the probe level or at the
  expression value level, not both.

  Typing \code{normalize.ExpressionSet.methods} should give you a list of
  methods that you may use. note that you can also use the
  \code{normalize} function on ExpressionSets. Use \code{method} to select the
  normalization method.
}

\value{
  A normalized \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}}.
}

\examples{
if (require(affydata)) {
  data(Dilution)
  eset <- rma(Dilution, normalize=FALSE, background=FALSE)
  normalize(eset)
}
}


\author{Ben Bolstad, \email{bmb@bmbolstad.com}}

\references{Bolstad, BM (2004) \emph{Low Level Analysis of High-density
    Oligonucleotide Array Data: Background, Normalization and
    Summarization}. PhD Dissertation. University of California,
  Berkeley.}


\seealso{\code{\link[affy]{normalize}}}

\keyword{manip}
