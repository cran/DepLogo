% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DepLogoR.R
\name{getPWM}
\alias{getPWM}
\alias{getPWM.DLData}
\title{Position weight matrix from DLData object}
\usage{
getPWM(part)

\method{getPWM}{DLData}(part)
}
\arguments{
\item{part}{the DLData object}
}
\value{
the position weight matrix, where columns correspond to positions
  (columns of the DLData$data slot) and rows to symbols
}
\description{
Determines the position weight matrix from a DLData object as relative frequency of symbols
in each column of the data slot.
}
\examples{
data <- DLData(c("ACGT", "ATTA"))
getPWM(data)
}
\author{
Jan Grau <grau@informatik.uni-halle.de>
}
