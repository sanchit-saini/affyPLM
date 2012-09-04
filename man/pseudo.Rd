\name{pseudo.coloring}
\alias{pseudoColorBar}
\alias{pseudoPalette}
\title{Coloring pseudo chip images}
\description{
  These are routines used for coloring pseudo chip images.
}
\usage{
  pseudoPalette(low = "white", high = c("green", "red"), mid = NULL,k =50) 
  pseudoColorBar(x, horizontal = TRUE, col = heat.colors(50), scale = 1:length(x),k = 11, log.ticks=FALSE,...)
  
}
\arguments{
  \item{low}{color at low end of scale}
  \item{high}{color at high end of scale}
  \item{mid}{color at exact middle of scale}
  \item{k}{number of colors to have}
  \item{x}{A data series}
  \item{horizontal}{If \code{TRUE} then color bar is to be draw
    horizontally}
  \item{col}{colors for color bar}
  \item{scale}{tickmarks for \code{x} if \code{x} is not numeric}
  \item{log.ticks}{use a log type transformation to assign the colors}
  \item{...}{additional parameters to plotting routine}
}

\details{
  Adapted from similar tools in maPlots pacakge.
}



\author{Ben Bolstad \email{bmb@bmbolstad.com}}

\seealso{\code{\link[affy:AffyBatch-class]{AffyBatch}}, \code{\link[affy:read.affybatch]{read.affybatch}}}

\keyword{manip}
