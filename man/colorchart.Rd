% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DepLogoR.R
\name{colorchart}
\alias{colorchart}
\title{Plots a colorchart representation of a set of sequences}
\usage{
colorchart(part, yoff, ic.scale = TRUE)
}
\arguments{
\item{part}{the set of sequences as \link{DLData} object}

\item{yoff}{the offset in y-direction within the current plot}

\item{ic.scale}{ignored for colorcharts}
}
\value{
the vertical (y) offset after this plot
}
\description{
This function is a low-level plotting function (using \link[graphics]{image} with \code{add=TRUE}, internally).
}
\examples{
# read data and create DLData object
seqs <- read.table(system.file("extdata", "cjun.txt", package = "DepLogo"), 
    stringsAsFactors = FALSE)
data <- DLData(sequences = seqs[,1], weights = log1p(seqs[, 2]) )

# create high-level plot
plot(NULL, xlim = c(1, ncol(data$data) - 1), ylim = c(0, nrow(data$data)), 
    ylab = nrow(data$data), axes = FALSE)
# and add colorchart and axis
colorchart(data, yoff = nrow(data$data))
axis(1)
}
\author{
Jan Grau <grau@informatik.uni-halle.de>
}
