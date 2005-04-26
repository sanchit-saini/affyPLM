\name{PLMset-class}
\docType{class}
\alias{PLMset-class}
\alias{PLMset}
\alias{weights<-}
\alias{coefs}
\alias{coefs<-}
\alias{se}
\alias{se<-}
\alias{coefs.probe}
\alias{coefs.const}
\alias{se.probe}
\alias{se.const}
\alias{Mbox}
\alias{normvec}
\alias{varcov}
\alias{residSE}
\alias{resid<-}
\alias{residuals<-}
\alias{indexProbesProcessed}
\alias{model.description}
\alias{sampleNames<-}
\alias{coefs,PLMset-method}
\alias{weights<-,PLMset-method}
\alias{weights,PLMset-method}
\alias{coefs<-,PLMset-method}
\alias{se,PLMset-method}
\alias{se<-,PLMset-method}
\alias{coefs.probe,PLMset-method}
\alias{coefs.const,PLMset-method}
\alias{se.probe,PLMset-method}
\alias{se.const,PLMset-method}
\alias{Mbox,PLMset-method}
\alias{image,PLMset-method}
\alias{show,PLMset-method}
\alias{boxplot,PLMset-method}
\alias{getCdfInfo,PLMset-method}
\alias{indexProbes,PLMset,character-method}
\alias{cdfName,PLMset-method}
\alias{normvec,PLMset-method}
\alias{resid,PLMset-method}
\alias{residuals,PLMset-method}
\alias{resid<-,PLMset-method}
\alias{residuals<-,PLMset-method}
\alias{residSE,PLMset-method}
\alias{varcov,PLMset-method}
\alias{indexProbesProcessed,PLMset-method}
\alias{model.description,PLMset-method}
\alias{summary,PLMset-method}
\alias{sampleNames,PLMset-method}
\alias{sampleNames<-,PLMset-method}
\alias{nuse}
\alias{NUSE}
\alias{RLE}
\alias{nuse,PLMset-method}
\alias{NUSE,PLMset-method}
\alias{RLE,PLMset-method}

\title{Class PLMset}
\description{This is a class representation for Probe level Linear
  Models fitted to Affymetrix GeneChip probe level data.}
\section{Objects from the Class}{
  Objects can be created using the function \code{\link{fitPLM}}}
\section{Slots}{
  \describe{
    \item{\code{probe.coefs}:}{Object of class "matrix". Contains model
      coefficients related to probe effects}
    \item{\code{se.probe.coefs}:}{Object of class "matrix". Contains
      standard error estimates for the probe coefficients}  
    \item{\code{chip.coefs}:}{Object of class "matrix". Contains model
      coefficients related to chip (or chip level) effects for each fit.}
    \item{\code{se.chip.coefs}:}{Object of class "matrix". Contains
      standard error estimates for the chip coefficients}
    \item{\code{model.description}:}{Object of class "character". This
      string describes the probe level model fitted.}
    \item{\code{weights}:}{Object of class "matrix". Contains probe
      weights for each fit. The matrix has columns for chips and rows
      are probes}
    \item{\code{phenoData}:}{Object of class "phenoData" This is an
      instance of class \code{phenoData} 
      containing the patient (or case) level data. The columns of the
      pData slot of this entity represent
      variables and the rows represent patients or cases. }
    \item{\code{annotation}}{A character string identifying the
      annotation that may be used for the \code{exprSet} instance.}
    \item{\code{description}:}{Object of class "MIAME". For
      compatibility with previous version of this class description can
      also be a "character". The class \code{characterOrMIAME} has been
      defined just for this.}
    \item{\code{cdfName}:}{A character string giving the name of the
      cdfFile}
    \item{\code{nrow}:}{Object of class "numeric". Number of rows in chip}
    \item{\code{ncol}:}{Object of class "numeric". Number of cols in chip}
    \item{\code{notes}:}{Object of class "character" Vector of
      explanatory text }
    \item{\code{varcov}:}{Object of class "list". A list of
      variance/covariance matrices}
    \item{\code{residualSE}:}{Object of class "matrix". Contains residual
      standard error and df}
    \item{\code{residuals}:}{Object of class "matrix". Contains
      residuals from model fit (if stored)}
  }
}
\section{Extends}{
Class \code{"exprSet"}, directly.
}
\section{Methods}{
  \describe{
    \item{weights<-}{\code{signature(object = "PLMset")}: replaces the weights }
    \item{weights}{\code{signature(object = "PLMset")}: extracts the
      model fit weights}
    \item{coefs<-}{\code{signature(object = "PLMset")}: replaces the
      chip coefs.  }
    \item{coefs}{\code{signature(object = "PLMset")}: extracts the
      chip coefs. }
    \item{se}{\code{signature(object = "PLMset")}: extracts the
      standard error estimates of the chip coefs. }
    \item{se<-}{\code{signature(object = "PLMset")}: replaces the
      standard error estimates of the chip coefs. }
    \item{coefs.probe}{\code{signature(object = "PLMset")}: extracts the
      probe coefs. }
    \item{se.probe}{\code{signature(object = "PLMset")}: extracts the standard error estimates of the probe coefs. }    
    \item{coefs.const}{\code{signature(object = "PLMset")}: extracts the
      intercept coefs. }
    \item{se.const}{\code{signature(object = "PLMset")}: extracts the standard error estimates of the intercept coefs. }
    \item{getCdfInfo}{\code{signature(object = "PLMset")}: retrieve
      the environment that defines the location of probes by probe set. }
    \item{image}{\code{signature(x = "PLMset")}: creates an
      image of the robust linear model fit weights for
      each sample.}
    \item{indexProbes}{\code{signature(object = "PLMset", which =
	"character")}: returns a list with locations of the probes in
      each probe  set. The list names defines the probe set
      names. \code{which} can be "pm", "mm", or "both". If "both" then
      perfect match locations are given followed by mismatch locations.}
    \item{Mbox}{\code{signature(object = "PLMset")}: gives a boxplot of
      M's for each chip. The M's are computed relative to a "median"
      chip.}
    \item{normvec}{\code{signature(x = "PLMset")}: will return the normalization vector (if it has been
      stored)}
    \item{residSE}{\code{signature(x = "PLMset")}: will return the residual SE (if it has been stored)}
    \item{boxplot}{\code{signature(x = "PLMset")}: Boxplot of Normalized
      Unscaled Standard Errors (NUSE)}
    \item{NUSE}{\code{signature(x = "PLMset")} : Boxplot of Normalized
      Unscaled Standard Errors (NUSE) or NUSE values}
    \item{RLE|}{\code{signature(x = "PLMset")} : Relative Log Expression
      boxplot or values}
 }
}
\note{This class is better described in the vignette.}

\author{B. M. Bolstad \email{bolstad@stat.berkeley.edu}}


\references{Bolstad, BM (2004) \emph{Low Level Analysis of High-density
    Oligonucleotide Array Data: Background, Normalization and
    Summarization}. PhD Dissertation. University of California,
  Berkeley.}

\examples{
}

\keyword{classes}
