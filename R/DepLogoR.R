#' Determines the position weight matrix from a DLData object as relative frequency of symbols
#' in each column of the data slot.

#' @title Position weight matrix from DLData object
#'   
#' @param part the DLData object
#'   
#' @return the position weight matrix, where columns correspond to positions
#'   (columns of the DLData$data slot) and rows to symbols
#' @export
#' @author Jan Grau <grau@informatik.uni-halle.de>
#'   
#' @examples
#' data <- DLData(c("ACGT", "ATTA"))
#' getPWM(data)
getPWM <- function(part){
	UseMethod("getPWM", part)
}

getPWM.DLData <- function(part){
	alphabet <- part$alphabet$chars
	part <- part$data[, -ncol(x = part$data)]
	apply(part,2,function(a){
		t <- table(factor(x = a, levels = alphabet))
		t/sum(t)
	})
}

drawArc <- function(x, y, radius, col, ...){
	po <- seq(from = 0, to = pi, length = 100)
	a <- radius*cos(x = po)
	b <- radius*sin(x = po)
	al <- atan2(y = b, x = a)
	rad <- sqrt(x = a^2 + b^2)
	xp <- rad*cos(x = al) + x
	yp <- rad*sin(x = al) + y
	lines(x = xp, y = yp, col=col, ...)
}

getColor <- function(column, ic = 255, colors = c("green", "blue", "orange", "red")){
	vals <- col2rgb(col = colors)
	vals <- vals%*%column
	rgb(red = vals[1], green = vals[2], blue = vals[3], alpha = ic, maxColorValue = 255)
}

getICScale <- function(column){
	ic <- log2(x = length(x = column))
	max <- ic
	column <- column[column>0]
	ic <- ic + sum(column*log2(column))
	sqrt(x = ic/max)
}

#' Paritions data by most inter-dependent positions
#' 
#' Partitions \code{data} by the nucleotides at the most inter-dependent
#' positions as measures by pairwise mutual information. Paritioning is
#' performed recursively on the resulting subsets until i) the number of
#' sequences in a partition is less then \code{minElements}, ii) the average
#' pairwise dependency between the current position and \code{numBestForSorting}
#' other positions with the largest mutual information value drops below
#' \code{threshold}, or iii) \code{maxNum} recursive splits have already been
#' performed. If splitting results in smaller partitions than
#' \code{minElements}, these are added to the smallest partition with more than
#' \code{minElements} sequences.
#' 
#' @param data the data as \link{DLData} object
#' @param minElements the minimum number of elements to perform a further split
#' @param threshold the threshold on the average mutual information value
#' @param numBestForSorting the number of dependencies to other positions
#'   considered
#' @param maxNum the maximum number of recursive splits
#' @param sortByWeights if \code{TRUE}, partitions are ordered by their average
#'   weight value, if \code{false} by frequency of symbols at the partitioning
#'   position otherwise. If \code{NULL}, the \code{$sortByWeights} value of the
#'   \link{DLData} object is used
#' @param partition.by specify fixed positions to partition by
#'   
#' @return the partitions as list of \link{DLData} objects
#' @export
#' @author Jan Grau <grau@informatik.uni-halle.de>
#' 
#' @examples
#' # create DLData object
#' seqs <- read.table(system.file("extdata", "cjun.txt", package = "DepLogo"),
#'     stringsAsFactors = FALSE)
#' data <- DLData(sequences = seqs[, 1], weights = log1p(seqs[,2]) )
#' 
#' # partition data using default parameters
#' partitions <- partition(data)
#' 
#' # partition data using a threshold of 0.3 on the mutual 
#' # information value to the most dependent position, 
#' # sorting the resulting partitions by weight
#' partitions2 <- partition(data = data, threshold = 0.3, numBestForSorting = 1, sortByWeights = TRUE)
partition <- function(
	data,
	minElements = 10,
	threshold = 0.1,
	numBestForSorting = 3,
	maxNum = 6,
	sortByWeights = NULL,
	partition.by = NULL){
	UseMethod("partition", data)
}

partition.DLData<-function(
	data,
	minElements = 10,
	threshold = 0.1,
	numBestForSorting = 3,
	maxNum = 6,
	sortByWeights = NULL,
	partition.by = NULL){
	if( is.null(sortByWeights)){
		sortByWeights <- data$sortByWeights;
	}

	if(is.null(partition.by)){
		temp <- partitionRecursive(data = data$data,
								   minElements = minElements,
								   threshold = threshold,
								   numBestForSorting = numBestForSorting,
								   maxNum = maxNum,
								   sortByWeights = sortByWeights,
								   alphabet = data$alphabet$chars,
								   exclude = rep(FALSE, ncol(data$data) - 1))
	}else{
		temp <- split(x=data$data, f=data$data[, partition.by], drop=TRUE)
		temp <- joinSmall(partSort = temp, minElements = minElements, sortByWeights = FALSE)
		if(sortByWeights){
			ws <- sapply(temp,function(a){mean(a$weights)})
			o <- order(ws, decreasing = TRUE)
			temp <- temp[ o ]
		}
	}
	lapply(temp, function(a){
		li <- list(data = a, alphabet = data$alphabet, sortByWeights = sortByWeights, axis.labels = data$axis.labels)
		class(li) <- "DLData"
		li
	})
}

# data: data.frame with last column = weights
partitionRecursive <- function(data, 
							   minElements = 10, 
							   threshold = 0.1,
							   numBestForSorting = 3,
							   maxNum = 6,
							   sortByWeights = FALSE,
							   alphabet = c("A","C","G","T"),
							   exclude = rep(FALSE, ncol(data))){
	sortTemp <- data
	if(maxNum <= 0 | nrow(sortTemp) < minElements | sum(!exclude)<2){
		return( list(sortTemp) )
	}else{
		pair <- getInformation(data = sortTemp, numBestForSorting = numBestForSorting, exclude = exclude, alphabet=alphabet)
		inf <- pair$deps
		deps <- pair$mis
		best <- which.max(inf)

		if( inf[best]/nrow(sortTemp)/numBestForSorting < threshold){
			return(list(sortTemp))
		}else{
			exclude[best] <- TRUE
			
			partSort <- partition.2(sortTemp = sortTemp, curr = best, minElements = minElements, sortByWeights = sortByWeights)
			
			kls <- deps[best,]/nrow(sortTemp)
			o2 <- order(kls, decreasing = TRUE)
			secpos <- o2[ !exclude[o2] & kls[o2]>threshold ]
			
			if( length(secpos)>0 ){
				secpos <- secpos[1]

				li <- list();
				for(i in 1:length(partSort)){
					temp2 <- partition.2(sortTemp = partSort[[i]], curr = secpos, minElements = minElements, sortByWeights = sortByWeights)
					temp2 <- joinSmall(partSort = temp2, minElements = minElements, sortByWeights = sortByWeights)
					li <- c(li, temp2)
				}
				exclude[secpos] <- TRUE
				maxNum <- maxNum-1
				partSort <- li
			}
			
			partSort <- joinSmall(partSort = partSort, minElements = minElements, sortByWeights = FALSE)
			
			partitions <- list()
			for(i in 1:length(partSort)){
				
				part <- partitionRecursive(data = partSort[[i]],
										   minElements = minElements,
										   threshold = threshold,
										   numBestForSorting = numBestForSorting,
										   maxNum = maxNum - 1,
										   sortByWeights = sortByWeights,
										   alphabet = alphabet,
										   exclude = exclude)
				partitions <- c(partitions, part)	
			}

			return(partitions)
		}
	}
}

joinSmall <- function(partSort, minElements, sortByWeights){
	if(length(partSort)==1){
		return(partSort)
	}else{
		out <- rep(FALSE, length(partSort))
		nout <- 0
		minAbove <- Inf
		idx <- (-1)
		
		for(i in 1:length(partSort)){
			if(nrow(partSort[[i]]) < minElements){
				out[i] <- TRUE;
				nout <- nout + nrow(partSort[[i]])
			}else{
				if( nrow(partSort[[i]]) < minAbove ){
					minAbove <- nrow(partSort[[i]])
					idx <- i
				}
			}
		}
		
		if( nout == 0 & idx == -1 ){
			return(partSort)
		}else if( idx == -1 ){
			minAbove <- 0;
		}
		
		joined <- c();
		if(idx>-1){
			joined <- partSort[[idx]];
		}else{
			idx <- which(out)[1]
		}
		
		for(i in 1:length(partSort)){
			if(out[i]){
				joined <- rbind(joined, partSort[[i]])
			}
		}
		
		partSort[[idx]] <- joined
		out[idx] <- FALSE
		
		partSort <- partSort[!out]
		
		if(sortByWeights){
			meanw <- sapply(partSort, function(a){
				mean(a[, ncol(a)])
			})
			o <- order(meanw, decreasing = TRUE);
			partSort <- partSort[o];
		}
		
		return(partSort)
	}
}



#' Computes the dependencies (as measures by mutual information) between all
#' positions (columns) of discrete data. Specifically, it returns for each pair
#' of positions (i,j) the mutual information I(X_i,X_j) multiplied by the number
#' N of sequences (rows), which may also be used for testing the statistical
#' significance of mutual information values, as for large N, 2*N*I(X_i,X_j) is
#' approximately chi squared.
#' 
#' @title Compute dependencies between positions
#' @param data the data for computing mutual information. Either a DLData object
#'   or a data.frame; In the latter case, the symbols of the alphabet must be
#'   provided as a second parameter
#' @param ... the symbols of the alphabet as character vector, only if data is a
#'   data.frame
#'   
#' @return a matrix of the mutual information values, where the diagonal is
#'   fixed to zero
#' @export
#' @author Jan Grau <grau@informatik.uni-halle.de>
#' 
#' @examples
#' data <- DLData(c("ACGT", "ATTA"))
#' deps <- getDeps(data)
#' 
#' 
getDeps <- function(data, ...){
	UseMethod("getDeps",data)
}

getDeps.DLData <- function(data, ...){
	getDeps(data = data$data, alphabet = data$alphabet$chars)
}

getDeps.data.frame <- function(data, alphabet, ...){
	x <- data[, -ncol(data)]
	x <- data.frame(lapply(x, factor, levels=alphabet))
	sum<-log( nrow(x) )
	
	mis<-sapply(1:ncol(x),function(i){
		sapply(1:ncol(x),function(j){
			tab <- table( x[, c(i, j)] )
			rs <- rowSums(tab)
			cs <- colSums(tab)

			ref <- outer(X = rs, Y = cs, FUN = function(x,y){
				log(x) + log(y) - sum
			})
			
			mi <- tab*(log(tab) - ref);
			mi <- sum(mi[tab>0])
			
			mi
			
		})
	})
	mis[mis<0] <- 0
	diag(mis) <- 0
	mis
}



getInformation <- function(data, numBestForSorting, exclude, alphabet){
	mis <- getDeps(data = data, alphabet = alphabet)
	
	mis2 <- mis[!exclude,]
	vals2 <- apply(mis2, 1, function(a){
		so <- sort(x = a, decreasing = TRUE)
		sum(so[1:numBestForSorting])
	})
	vals <- rep(0, nrow(mis))
	vals[!exclude] <- vals2
	list(deps = vals, mis = mis)
}

partition.2 <- function(sortTemp, curr, minElements, sortByWeights){
	if(nrow(sortTemp)<minElements){
		return(list(sortTemp))
	}else{
		parts <- split( sortTemp, sortTemp[, curr] )
		freq <- sapply(parts, nrow)
		ws <- sapply(parts, function(a){ sum(a[, ncol(a)]) })
		if(sortByWeights){
			freq <- ws/freq;
		}
		ord <- order(freq, decreasing = TRUE);
		return(parts[ord])
	}
	
}

addLegend <- function(minp, maxp, minp.col, maxp.col, axis.at.bottom = TRUE, pvals = FALSE){
	bak <- par("mar")
	on.exit(par(bak))

	vals <- seq(minp, maxp, length = 20)
	step <- (maxp - minp)/19
	cols <- rgb(t(sapply(seq(0, 1, length = 20), function(a){a * (maxp.col - minp.col) + minp.col})), maxColorValue = 255)
	
	pmar <- par("mar")
	pmar[3] <- 2
	par(new = TRUE, mar = pmar)
	
	ylim <- exp(c(minp, maxp))
	if(pvals){
		ylim <- 10^-c(minp, maxp)
	}
	xlim <- c(minp, maxp + 3 * (maxp - minp))
	
	if(pvals){
		xlim <- xlim*-log(10)
		vals <- vals*-log(10)
		step <- step*-log(10)
	}
	
	if(axis.at.bottom){
		plot(NA, ylim = ylim, xlim = xlim, xlab = "", ylab = "", axes = FALSE, log = "y", yaxs = "i")
		ticks <- axTicks(2)
		axis(3, at = log(ticks), labels = ticks, las=0)
	}else{
		plot(NA, ylim = rev(ylim), xlim = xlim, xlab = "", ylab = "", axes = FALSE, log = "y", yaxs = "i")
		ticks <- axTicks(2)
		axis(1, at = log(ticks), labels = ticks, las = 0)
	}
	
	ybottom <- exp(maxp - (0.1) * (maxp - minp))
	ytop <- exp(maxp)
	if(pvals){
		ybottom <- 10^-(maxp - (0.1) * (maxp - minp))
		ytop <- 10^-maxp
	}
	for(i in 1:length(vals)){
		rect(xleft = vals[i] - step/2, xright = vals[i] + step/2, ybottom = ybottom, ytop = ytop, col = cols[i], border = cols[i])
	}
	
}



#' Plots a matrix representation of dependency values
#' 
#' Plots a representation of dependency values as a triangular matrix rotated by
#' 45 degrees. Internally, dependency values are computed using \link{getDeps}
#' on the \code{data} object.
#' 
#' @param data the \link{DLData} object containing the data
#' @param axis.at.bottom if \code{TRUE}, the x-axis is shown at the bottom
#'   (side=1) of the plot, and at the top (side=3) otherwise
#' @param add.legend if \code{TRUE} a legend of the color scale is added to the
#'   plot
#' @param show.pvals if \code{TRUE}, -log10 p-values (computed by
#'   \link[stats]{pchisq}) are shown instead of mutual information values
#' @param axis.labels the labels of the x-axis
#' @param threshold ignored
#'   
#' @export
#' @author Jan Grau <grau@informatik.uni-halle.de>
#'   
#' @examples
#' # create DLData object
#' seqs <- read.table(system.file("extdata", "cjun.txt", package = "DepLogo"), 
#'     stringsAsFactors = FALSE)
#' data <- DLData(sequences = seqs[, 1], weights = log1p(seqs[, 2]) )
#' 
#' # plot using default parameters
#' plotDepmatrix(data)
#' 
#' # plot with axis at top, without a legend (color scale), and using p-values
#' plotDepmatrix(data, axis.at.bottom = FALSE, add.legend = FALSE, show.pvals = TRUE)
plotDepmatrix <- function(data,
						  axis.at.bottom = TRUE,
						  add.legend = TRUE,
						  show.pvals = FALSE,
						  axis.labels = NULL,
						  threshold = 0.1){
	if(is.null(axis.labels)){
		axis.labels <- data$axis.labels
	}
	alphabet <- data$alphabet$chars
	if(show.pvals){
		stat <- getDeps(data) * 2
		mis <- (-log10(pchisq(stat, df = (length(alphabet) - 1)^2, lower.tail = FALSE)) )
		rang<-range(mis[!is.infinite(mis)], na.rm = TRUE);
		if(rang[2] > rang[1]){
			mis[is.infinite(mis)] <- rang[2]
		}else{
			mis[mis>100] <- 100
		}

		minp <- min(mis[upper.tri(mis)], na.rm = TRUE)
		maxp <- max(mis[upper.tri(mis)], na.rm = TRUE)
		if(maxp==minp){
			minp <- max(0, minp - 0.5)
			maxp <- maxp + 0.5
		}
	}else{
		mis <- log(getDeps(data) / nrow(data$data))
		
		minp <- log(0.01)
		maxp <- log(log(length(alphabet)))
	}
	
	if(is.null(axis.labels)){
		axis.labels <- 1:ncol(mis)
	}
	
	if(axis.at.bottom){
		plot(x = NA, xlim = c(0.5,ncol(mis) + 0.5), ylim = c(0, ncol(mis)/2), axes = FALSE, xlab = "", ylab = "", xaxs = "i")
		axis(side = 1, tcl = .5, line = 2, at = 1:ncol(mis), labels = NA)
		axis(side = 1, col=0, line = -.3 * max(nchar(as.character(axis.labels))), at = 1:ncol(mis), labels = axis.labels)
	}else{
		plot(x = NA, xlim = c(0.5, ncol(mis) + 0.5), ylim = c(ncol(mis)/2, 0), axes=FALSE, xlab = "", ylab = "", xaxs = "i")
		axis(side = 3, tcl = .5, line = 2, at = 1:ncol(mis), labels = NA)
		axis(side = 3, col = 0, line = -.3 * max(nchar(as.character(axis.labels)))/2, at = 1:ncol(mis), labels = axis.labels)
	}
	
	minp.col <- t(col2rgb("white"))
	maxp.col <- t(col2rgb("black"))
	for(i in 2:nrow(mis)){
		for(j in 1:(i-1)){
			mix <- max(0, (mis[i,j] - minp)/(maxp - minp))
			if(is.na(mix)){
				print(c(i, j, mix, mis[i,j], maxp, minp))
			}
			
			col <- rgb(mix * maxp.col + (1-mix) * minp.col, maxColorValue = 255)
			cent.x <- (i + j)/2
			cent.y <- (i - j)*0.5
			polygon(x = c(cent.x - 0.5, cent.x, cent.x + 0.5, cent.x), y = c(cent.y, cent.y + 0.5, cent.y, cent.y - 0.5), col = col, border = "grey")
		}
	}
	
	if(add.legend){	
		addLegend(minp, maxp, minp.col, maxp.col, axis.at.bottom, pvals = show.pvals)
	}
}

#' Plots a graph representation of dependency values
#' 
#' Plots a representation of dependency values as arcs between the sequence
#' positions. Internally, dependency values are computed using \link{getDeps} on
#' the \code{data} object.
#' 
#' @param data the \link{DLData} object containing the data
#' @param axis.at.bottom if \code{TRUE}, the x-axis is shown at the bottom
#'   (side=1) of the plot, and at the top (side=3) otherwise
#' @param add.legend if \code{TRUE} a legend of the color scale is added to the
#'   plot
#' @param show.pvals if \code{TRUE}, -log10 p-values (computed by
#'   \link[stats]{pchisq}) are shown instead of mutual information values
#' @param axis.labels the labels of the x-axis
#' @param threshold threshold in mutual information values, edges below this value are not shown; ignored in \code{show.pvals=TRUE}
#'   
#' @export
#' @author Jan Grau <grau@informatik.uni-halle.de>
#'   
#' @examples
#' # create DLData object
#' seqs <- read.table(system.file("extdata", "cjun.txt", package = "DepLogo"), 
#'     stringsAsFactors = FALSE)
#' data <- DLData(sequences = seqs[,1], weights = log1p(seqs[, 2]) )
#' 
#' # plot using default parameters
#' plotDeparcs(data)
#' 
#' # plot with axis at top, without a legend (color scale), and using p-values
#' plotDeparcs(data, axis.at.bottom = FALSE, add.legend = FALSE, show.pvals = TRUE)
plotDeparcs<-function(data,
					  axis.at.bottom = TRUE,
					  add.legend = TRUE,
					  show.pvals = FALSE,
					  axis.labels = NULL,
					  threshold = 0.1){
	if(is.null(axis.labels)){
		axis.labels <- data$axis.labels
	}
	alphabet <- data$alphabet$chars
	if(show.pvals){
		stat <- getDeps(data) * 2
		stat[stat<threshold * 2 * nrow(data$data)] <- NA
		
		mis<- (-log10(pchisq(stat,df = (length(alphabet)-1)^2, lower.tail = FALSE)) )
		suppressWarnings( rang <- range(mis[!is.infinite(mis)], na.rm = TRUE) )
		if(rang[2] > rang[1]){
			mis[is.infinite(mis)] <- rang[2]
		}else{
			mis[mis>100] <- 100
		}
		
		minp <- min(mis[upper.tri(mis)], na.rm = TRUE)
		maxp <- max(mis[upper.tri(mis)], na.rm = TRUE)
		if(maxp==minp){
			minp <- max(0,minp - 0.5)
			maxp <- maxp + 0.5
		}
	}else{
		mis <- log(getDeps(data) / nrow(data$data))
		
		minp <- log(threshold)
		maxp <- log(log(length(alphabet)))
	}
	
	maxd <- 0;
	temp <- which(mis>minp, arr.ind = TRUE)
	if(length(temp)>0){
		maxd <- max(abs(temp[, 1] - temp[, 2]))
	}
	
	max.y <- (maxd + 1)/2
	if(add.legend){
		max.y <- max.y * 1.5
	}
	
	if(is.null(axis.labels)){
		axis.labels <- 1:ncol(mis)
	}
	
	if(axis.at.bottom){
		plot(x = NA, xlim = c(0.5,ncol(mis) + 0.5), ylim = c(0, max.y), axes = F, xlab = "", ylab = "", xaxs = "i")
		axis(side = 1, tcl = .5, line = 2, at = 1:ncol(mis), labels = NA)
		axis(side = 1, col = 0, line = -.3*max(nchar(as.character(axis.labels))), at = 1:ncol(mis), labels = axis.labels)
	}else{
		plot(x = NA, xlim = c(0.5, ncol(mis) + 0.5), ylim = c(max.y, 0),axes = F, xlab = "", ylab = "", xaxs = "i")
		axis(side = 3, tcl = .5, line = 2, at=1:ncol(mis), labels = NA)
		axis(side = 3, col = 0, line = -.3 * max(nchar(as.character(axis.labels))), at = 1:ncol(mis), labels = axis.labels)
	}
	ypos <- 0
	
	for(i in 2:nrow(mis)){
		for(j in 1:(i-1)){
			if(!is.na(mis[i, j]) & mis[i, j]>minp){
				m <- (i + j)/2
				rad <- abs(i - j)/2
				
				drawArc(x = m, y = ypos, radius = rad, col = gray(0, min(1, (mis[i,j] - minp)/(maxp - minp) ) ), lwd = 2)
			}
		}
	}

	if(add.legend){	
		addLegend(minp, maxp, col2rgb("white"), col2rgb("black"), axis.at.bottom, pvals = show.pvals)
	}
	
}



#' Plots a colorchart representation of a set of sequences
#' 
#' This function is a low-level plotting function (using \link[graphics]{image} with \code{add=TRUE}, internally).
#'
#' @param part the set of sequences as \link{DLData} object
#' @param yoff the offset in y-direction within the current plot
#' @param ic.scale ignored for colorcharts
#'
#' @return the vertical (y) offset after this plot
#' @export
#' @author Jan Grau <grau@informatik.uni-halle.de>
#'
#' @examples
#' # read data and create DLData object
#' seqs <- read.table(system.file("extdata", "cjun.txt", package = "DepLogo"), 
#'     stringsAsFactors = FALSE)
#' data <- DLData(sequences = seqs[,1], weights = log1p(seqs[, 2]) )
#' 
#' # create high-level plot
#' plot(NULL, xlim = c(1, ncol(data$data) - 1), ylim = c(0, nrow(data$data)), 
#'     ylab = nrow(data$data), axes = FALSE)
#' # and add colorchart and axis
#' colorchart(data, yoff = nrow(data$data))
#' axis(1)
colorchart<-function(part, yoff, ic.scale = TRUE){

	mat <- as.matrix(part$data[, -ncol(part$data)])
	map <- 1:length(part$alphabet$chars)
	names(map) <- part$alphabet$chars
	mat <- matrix(map[mat], ncol = ncol(mat))

	image(x = 1:ncol(mat), y = (yoff - nrow(mat)):yoff, t(mat), col = part$alphabet$cols, add = TRUE)
	yoff - nrow(mat)
}

#' Plots a representation of a set of sequences by rectangles of (scaled) averaged color values of the symbols at each position
#' 
#' This function is a low-level plotting function (using \link[graphics]{rect}, internally).
#'
#' @title Rectangles of averaged colors
#' @param part the set of sequences as \link{DLData} object
#' @param yoff the offset in y-direction within the current plot
#' @param ic.scale if \code{TRUE}, alpha values of colors will be assigned based on "information content" of the distribution at each position
#'
#' @return the vertical (y) offset after this plot
#' @export
#' @author Jan Grau <grau@informatik.uni-halle.de>
#'
#' @examples
#' # read data and create DLData object
#' seqs <- read.table(system.file("extdata", "cjun.txt", package = "DepLogo"),  
#'    stringsAsFactors = FALSE)
#' data <- DLData(sequences = seqs[, 1],weights = log1p(seqs[, 2]) )
#' 
#' # create high-level plot
#' plot(NULL, xlim = c(1, ncol(data$data) - 1), ylim = c(0, nrow(data$data)), 
#'     ylab = nrow(data$data), axes = FALSE)
#' # and add deprects and axis
#' deprects(data, yoff = nrow(data$data))
#' axis(1)
deprects<-function(part, yoff, ic.scale = TRUE){
	pwm <- getPWM.DLData(part)
	size <- nrow(part$data)
	sapply(1:ncol(pwm), function(i){
		color <- getColor(pwm[, i], ifelse(ic.scale, getICScale(pwm[, i]) * 255, 255), part$alphabet$cols);
		rect(i - .5, yoff - size, i + .5, yoff, col = color, border=NA)
	})
	yoff - size
}

#' Plots a representation of a set of sequences as a sequence logo
#' 
#' This function is a low-level plotting function (using \link[graphics]{polygon}, internally).
#'
#' @title Sequence logo
#' @param part the set of sequences as \link{DLData} object
#' @param yoff the offset in y-direction within the current plot
#' @param ic.scale if \code{TRUE}, symbols are scaled by "information content" of the distribution at each position
#'
#' @return the vertical (y) offset after this plot
#' @export
#' @author Jan Grau <grau@informatik.uni-halle.de>
#'
#' @examples
#' # read data and create DLData object
#' seqs <- read.table(system.file("extdata", "cjun.txt", package = "DepLogo"), 
#'     stringsAsFactors = FALSE)
#' data <- DLData(sequences = seqs[, 1], weights = log1p(seqs[,2]) )
#' 
#' # create high-level plot
#' plot(NULL, xlim = c(1, ncol(data$data) - 1), ylim = c(0, nrow(data$data)),  
#'    ylab = nrow(data$data), axes = FALSE)
#' # and add sequence logo and axis
#' logo(data, yoff = nrow(data$data))
#' axis(1)
logo <- function(part, yoff, ic.scale = TRUE){
	pwm <- getPWM(part)
	size <- nrow(part$data)
	alphabet <- part$alphabet
	
	letters <- list(x = NULL, y = NULL, id = NULL, fill = NULL)
	npos <- ncol(pwm)
	wt <- 1
	x.pos <- 0.5
	eps <- 0
	heights <- c()
	ymins <- c()
	ymaxs <- c()
	for (j in 1:npos) {
		column <- pwm[, j]
		sh <- ifelse(ic.scale, getICScale(column)^2, 1)
		hts <- column * sh * size
		letterOrder <- order(abs(hts))
		ypos.pos <- yoff - size
		
		hts <- hts[letterOrder]
		chars <- alphabet$chars[letterOrder]
		cols <- alphabet$cols[letterOrder]
		
		for (i in 1:alphabet$size) {
			ht <- hts[i]
			y.pos <- ypos.pos
			ht <- ht - eps
			ypos.pos <- ypos.pos + ht + eps
			char <- chars[i]
			col <- cols[i]
			let <- getLetter(letterPolygons[[char]], x.pos, y.pos, ht, wt*0.99, col = col)
			polygon(let, col = let$col, border = NA)
		}
		
		x.pos <- x.pos + wt
	}
	yoff - size
}

#' Plots blocks of data
#' 
#' Plots the blocks of data in \code{data} by successive, vertically arranged sub-plots of the function provided as \code{block.fun}.
#' If \code{data} is a single \link{DLData} object, one block is plotted. Further arguments are provided to \code{block.fun}.
#'
#' @param data the data, a single \link{DLData} object or a list of \code{DLData} objects
#' @param show.number  if \code{true}, the number of sequences (in total) in data is displayed on the left side of the plot
#' @param block.fun the function called for each of the blocks
#' @param ic.scale if \code{TRUE}, output of \code{block.fun} may be scaled by "information content"
#' @param add if \code{TRUE}, the plot is added to an existing plot
#' @param ... if \code{add=FALSE} forwarded to the internal call to \link[graphics]{plot}
#'
#' @export
#' @author Jan Grau <grau@informatik.uni-halle.de>
#' @seealso \link{deprects}
#' @seealso \link{logo}
#' @seealso \link{colorchart}
#'
#' @examples
#' # read data and create DLData object
#' seqs <- read.table(system.file("extdata", "cjun.txt", package = "DepLogo"), 
#'     stringsAsFactors = FALSE)
#' data <- DLData(sequences = seqs[, 1], weights = log1p(seqs[, 2]) )
#' 
#' # plot all data
#' plotBlocks(data)
#' 
#' # partition data
#' partitions <- partition(data, threshold = 0.3)
#' # and plot partitions
#' plotBlocks(partitions)
#' 
#' # or plot partitions as sequence logos
#' plotBlocks(partitions, block.fun = logo)
plotBlocks <- function(data, show.number = TRUE, block.fun = deprects, ic.scale = TRUE, add = FALSE, ...){
	UseMethod("plotBlocks", data)
}

plotBlocks.DLData <- function(data, show.number = TRUE, block.fun = deprects, ic.scale = TRUE, add = FALSE, ...){
	plotBlocks.list(data = list(data), show.number = show.number, block.fun = block.fun, ic.scale = ic.scale, add = add, ...)	
}

plotBlocks.list <- function(data, show.number = TRUE, block.fun = deprects, ic.scale = TRUE, add=FALSE, ...){
	total <- sum(sapply(data, function(a){nrow(a$data)}))

	len <- ncol(data[[1]]$data) - 1
	if(!add){
		plot(NULL, xlim = c(0.5, len + 0.5), ylim = c(0, total + 2), axes = F, yaxs = 'i', xaxs = "i", xlab = "", ylab = "", ...)
	}
	if(show.number){
		mtext(side = 2, line = 1, text = paste("N =", total), cex = par("cex"))
	}

	t <- sapply(data, function(a){
		total <<- block.fun(a, total, ic.scale)
	})
	if(length(data)>1){
		lines(x = c(0.5,len + 0.5), y = c(total, total), col = 1, lwd = 2)
	}
}


#' Plots a dependency logo
#' 
#' The function \code{dep.fun} provided for plotting the representation of
#' dependencies is currently implemented in \link{plotDeparcs} and
#' \link{plotDepmatrix}. Custom implementations must have the same signature as
#' these functions and create a single plot without using
#' \link[graphics]{layout} (or similar).
#' 
#' The functions \code{block.fun} and \code{summary.fun} provided for plotting
#' the representation of individual partitions of the data generated in
#' dependency logos are currently implemented in \link{deprects},
#' \link{colorchart}, and \link{logo}. Custom implementations must have the same
#' signature as these functions and create a single plot without using
#' \link[graphics]{layout} (or similar).
#' 
#' The function \code{weight.fun} for plotting a representation of the
#' \code{weights} values of the sequences within one partition is currently
#' implemented in \link{subLines} and \link{subBoxes}. Custom implementations
#' must have the same signature as these functions and create a single plot
#' without using \link[graphics]{layout} (or similar).
#' 
#' @title Plot a dependency logo
#' @param data the data, currently implemented for \link{DLData} objects
#' @param dep.fun the function for plotting the representation of dependency
#'   values (as computed by \link{getDeps})
#' @param block.fun the function for plotting a representation of the individual
#'   partitions of the data generated in dependency logos.
#' @param summary.fun the function for plotting a representation of the summary
#'   plot for (one chunk of) the data
#' @param weight.fun the function for plotting a representation of the
#'   \code{weights} values of the sequences within one partition
#' @param chunks the size of chunks the data is split into. The sum of the chunk
#'   sizes must not be greater than the number of data points in data; The
#'   default value of NULL corresponds to one chunk containing all data points
#' @param chunk.height the (relative) height of the parts of the plot
#'   representing each of the chunks, one height for each chunk
#' @param summary.height the (relative) height of the block summaries in the
#'   plot
#' @param minPercent the minimum percentage of the (sub) data set that may
#'   constitute its own partition in the dependency logo
#' @param threshold the threshold on the dependency value for further splits
#' @param numBestForSorting the number of dependencies between position i and
#'   all other positions when computing the dependency value of position i
#' @param maxNum the maximum number of splits allowed
#' @param sortByWeights are partitions sorted by their average weight
#'   (descending)
#' @param dep.fun.legend if \code{TRUE}, a legend of the color scale used for
#'   plotting the dependency values in \code{dep.fun} is added to the plot
#' @param show.dependency.pvals is \code{TRUE}, p-values are used for plotting
#'   dependency values in \code{dep.fun} instead of mutual information values
#' @param axis.labels labels for the x-axis, vector of the same length as the
#'   individual sequences
#' @param weight.ratio the factor by which the plotting width for the main plot is larger than for \code{weight.fun}
#' @param partition.by specify fixed positions to partition by
#' @param ... forwarded to the high-level \code{plot} that contains the blocks
#'   plotted by \code{block.fun}
#'   
#' @return a list of \link{DLData} objects with the partitions created for the
#'   dependency logo
#' @export
#' @author Jan Grau <grau@informatik.uni-halle.de>
#'   
#' @examples
#' # read data and create DLData object
#' seqs <- read.table(system.file("extdata", "cjun.txt", package = "DepLogo"), 
#'     stringsAsFactors = FALSE)
#' data <- DLData(sequences = seqs[, 1],weights = log1p(seqs[, 2]) )
#' 
#' # plot default dependency logo
#' plotDeplogo(data)
#' 
#' # refine threshold for clearer picture
#' plotDeplogo(data, threshold = 0.3)
#' 
#' # customize different parts of the plot
#' plotDeplogo(data, threshold = 0.3, dep.fun = plotDepmatrix, block.fun = colorchart)
#' 
#' # add plots of the weights
#' plotDeplogo(data, weight.fun = subBoxes)
plotDeplogo<-function(data,
					  dep.fun = plotDeparcs,
					  block.fun = deprects,
					  summary.fun = logo,
					  weight.fun = NULL,
					  chunks = NULL,
					  chunk.height = 800,
					  summary.height = 100,
					  minPercent = 0.03,
					  threshold = 0.1,
					  numBestForSorting = 3,
					  maxNum = 6,
					  sortByWeights = NULL,
					  dep.fun.legend = TRUE,
					  show.dependency.pvals = FALSE,
					  axis.labels = NULL,
					  weight.ratio = 5,
					  partition.by = NULL,
					  ...){
	UseMethod("plotDeplogo", data)
}

plotDeplogo.DLData<-function(data,
							dep.fun = plotDeparcs,
							block.fun = deprects,
							summary.fun = logo,
							weight.fun = NULL,
							chunks = NULL,
							chunk.height = 800,
							summary.height = 100,
							minPercent = 0.03,
							threshold = 0.1,
							numBestForSorting = 3,
							maxNum = 6,
							sortByWeights = NULL,
							dep.fun.legend = TRUE,
							show.dependency.pvals = FALSE,
							axis.labels = NULL, 
							weight.ratio = 5,
							partition.by = NULL,
							...){
	
	if(is.null(chunks)){
		chunks <- nrow(data$data)
	}
	if(sum(chunks)>nrow(data$data)){
		stop("Chunk sizes larger than data.")
	}
	if(length(chunks) != length(chunk.height)){
		stop("Not the same number of chunks as number of chunk sizes")
	}
	
	if(is.null(sortByWeights)){
		sortByWeights <- data$sortByWeights;
	}
	
	if(sortByWeights && !data$sortByWeights){
		o <- order(data$data$weights, decreasing = TRUE);
		data$data <- data$data[o, ]
	}
	
	if(is.null(axis.labels)){
		axis.labels <- data$axis.labels
	}
	
	parts <- split(x = data$data, f = rep(1:length(chunks), chunks))
	parts <- lapply(parts, function(a){
		li <- list(data = a, alphabet = data$alphabet, sortByWeights = sortByWeights)
		class(li) <- "DLData"
		li
	})
	
	total.height <- sum(chunk.height) + summary.height * (length(chunks) - 1) + summary.height * 3.5
	numPlots <- length(chunks) * 2 + 1 + 1
	height <- c(summary.height * 3, as.vector(rbind(chunk.height, rep(summary.height, length(chunk.height)))), summary.height)/total.height

	bak <- par(no.readonly = TRUE)
	on.exit(par(bak))
	par(mar = c(2, 2.5, 0, ifelse(!is.null(weight.fun), 0, 1)))
	
	if(!is.null(weight.fun)){
		layout(mat = matrix(1:(2 * numPlots), ncol=2), widths = c(ncol(data$data) - 1, (ncol(data$data) - 1)/weight.ratio), heights = height)
	}else{
		layout(mat = matrix(1:numPlots, ncol = 1), widths = c(1), heights = height)
	}
	
	dep.fun( data = data, add.legend = dep.fun.legend, show.pvals = show.dependency.pvals, axis.labels = axis.labels, threshold = threshold)
	par(mar = c(0, 2.5, 0, ifelse(!is.null(weight.fun), 0, 1)))

	w.li <- list()
	sapply(1:length(parts), function(i){
		sub.parts <- partition.DLData(data = parts[[i]], minElements = nrow(parts[[i]]$data) * minPercent,
							 threshold = threshold, numBestForSorting = numBestForSorting, maxNum = maxNum, partition.by = partition.by)

		plotBlocks(data = sub.parts, show.number = TRUE, block.fun = block.fun, ic.scale = TRUE, add = FALSE, ...)
		plotBlocks(data = parts[[i]], show.number = FALSE, block.fun = summary.fun, ic.scale = TRUE, add = FALSE, ... )
	
		for(j in 1:length(sub.parts)){
			sub.parts[[j]]$data <- sub.parts[[j]]$data[order(sub.parts[[j]]$data$weights, decreasing = T), ]
		}
		
		w.li <<- c(w.li, list(sub.parts))
	})
	
	par(mar = c(1.5, 2.5, 0, ifelse(!is.null(weight.fun), 0, 1)))
	plot(x = NA, xlim = c(0.5, ncol(data$data) - 0.5), ylim = c(0, 1),
		 axes = F, frame.plot = F, xlab = "", ylab = "", yaxs = 'i', xaxs = "i")
	axis(side = 1, at = 1:(ncol(data$data) - 1), labels = axis.labels, cex = par("cex"), pos = 1)
	
	if(!is.null(weight.fun)){
		par(mar = c(0, .5, 0, 1))
		plot.new()
		range <- range(data$data$weights)
		
		for(i in 1:length(w.li)){
			el <- w.li[[i]]
			weight.fun(el, range, axis.above = (i==1), axis.below = TRUE)
			
			plot.new()
		}	
	}

	layout(1)
	invisible(w.li)
}



#' Plots weights as lines
#' 
#' Plots a representation of the weights of a list of \link{DLData} objects.
#' Each entry of the list is shown as an independent line with the median value
#' shown as a red vertical line. Plots of list entries are separated by
#' horizontal grey lines.
#' 
#' @param sub.parts a list of \link{DLData} objects
#' @param range the range of values shown in the plot (i.e., the \code{xlim}
#'   value of the call to \link[graphics]{plot})
#' @param axis.above if \code{TRUE}, an axis at the top of the plot (side=3) is shown
#' @param axis.below if \code{TRUE}, an axis at the bottom of the plot (side=1) is shown
#' 
#' @export
#' @author Jan Grau <grau@informatik.uni-halle.de>
#' 
#' @examples
#' # read data and create DLData object
#' seqs <- read.table(system.file("extdata", "nrsf.txt", package = "DepLogo"), 
#'     stringsAsFactors = FALSE)
#' data <- DLData(sequences = seqs[, 1], weights = log1p(seqs[, 2]) )
#' 
#' # create dependency logo with plotted weights
#' plotDeplogo(data, threshold = 0.03, weight.fun = subLines)
subLines <- function(sub.parts, range, axis.above = TRUE, axis.below = TRUE){
	total <- sum(sapply(sub.parts, function(a){nrow(a$data)}))
	plot(NULL, xlim = range, ylim = c(total, 1), axes = FALSE, yaxs = 'i')
	if(axis.below) axis(1)
	if(axis.above) axis(3)
	cum <- 1
	for(part in sub.parts){
		if(cum>1){
			abline(h = cum, lwd = 1, col = "grey")
		}
		lines(x = part$data$weights, y = cum:(cum + nrow(part$data) - 1))
		lines(x = c(median(part$data$weights), median(part$data$weights)), y = c(cum, cum + nrow(part$data) - 1), col = 2)
		cum <- cum + nrow(part$data)
	}
}

#' Plots weights as boxplots
#' 
#' Plots a representation of the weights of a list of \link{DLData} objects.
#' Each entry of the list is shown as an independent boxplot.
#' 
#' @param sub.parts a list of \link{DLData} objects
#' @param range the range of values shown in the plot (i.e., the \code{xlim}
#'   value of the call to \link[graphics]{plot})
#' @param axis.above if \code{TRUE}, an axis at the top of the plot (side=3) is shown
#' @param axis.below if \code{TRUE}, an axis at the bottom of the plot (side=1) is shown
#' 
#' @export
#' @author Jan Grau <grau@informatik.uni-halle.de>
#' 
#' @examples
#' # read data and create DLData object
#' seqs <- read.table(system.file("extdata", "nrsf.txt", package = "DepLogo"), 
#'     stringsAsFactors = FALSE)
#' data <- DLData(sequences = seqs[, 1], weights = log1p(seqs[, 2]) )
#' 
#' # create dependency logo with plotted weights
#' plotDeplogo(data, threshold = 0.03, weight.fun = subBoxes)
subBoxes <- function(sub.parts, range, axis.above = TRUE, axis.below = TRUE){
	total <- sum(sapply(sub.parts, function(a){nrow(a$data)}))
	plot(NULL, xlim = range, ylim = c(total, 1), axes = FALSE, yaxs = 'i')
	if(axis.below) axis(1)
	if(axis.above) axis(3)
	cum <- 1
	for(part in sub.parts){
		stats <- boxplot.stats(part$data$weights)
		lines(x = c(stats$stats[1], stats$stats[5]), y = c(cum + nrow(part$data)/2, cum + nrow(part$data)/2), lty = "dashed")
		lines(x = c(stats$stats[1], stats$stats[1]), y = c(cum + nrow(part$data) * 0.35, cum + nrow(part$data) * 0.65))
		lines(x = c(stats$stats[5], stats$stats[5]), y = c(cum + nrow(part$data) * 0.35, cum + nrow(part$data) * 0.65))
		rect(xleft = stats$stats[2], ybottom = cum + nrow(part$data) * 0.2, xright = stats$stats[4], ytop = cum + nrow(part$data) * 0.8, col = "white")
		lines(x = c(stats$stats[3], stats$stats[3]), y = cum + nrow(part$data) * c(0.2, 0.8),lwd = 2, col = 4)
		points(x = stats$out, y = rep(cum + nrow(part$data)/2, length(stats$out)), col = rgb(0,0,0,0.2))
		cum<-cum + nrow(part$data)
	}
}


length.DLData <- function(x){
	nrow(x$data)
}

dim.DLData <- function(x){
	dim(x$data)
}



#' Creates a new \code{DLData} object from a set of input sequences.
#' 
#' Sequences may either be provided as a \code{character} vector or as a
#' \code{data.frame}. All symbols occurring in these sequences need to be
#' defined and assigned to colors, which are used for plotting later. Colors do
#' not need to be unique, but symbols with identical colors may become
#' indistinguishable in subsequent plots (which might even be desired, for
#' instance, when visualizing protein properties instead of amino acids). 
#' Sequences may have an associated weight, which is used to order sequences,
#' e.g., for creating chunks/blocks of sequences in subsequent plots (see
#' \code{chunks} parameter of \code{plotDeplogo}).
#' 
#' @title Create \code{DLData} object
#' @param sequences the input sequences, may be provided as i) \code{character}
#'   vector or ii) a \code{data.frame} with sequences organized in rows and one 
#'   symbol per column
#' @param weights weights associated with the sequences, numeric vector of the 
#'   same length as \code{sequences} has sequences
#' @param symbols the symbols (alphabet) over which the sequences are defined
#' @param colors colors for each of the \code{symbols}, not necessarily unique
#' @param delim delimiter between the symbols in the input sequences, ignored if
#'   \code{sequences} as a \code{data.frame}
#' @param sortByWeights if \code{TRUE}, \code{sequences} will be ordered by 
#'   their \code{weight} in decreasing order
#' @param axis.labels the labels of the individual sequence positions; 
#'   if \code{NULL}, indexes from 1 to to total number of positions will be used
#'   
#' @return the \code{DLData} object
#' @export
#' @exportClass DLData
#' @seealso \link{plotDeplogo}
#' @author Jan Grau <grau@informatik.uni-halle.de>
#'   
#' @examples 
#' # creating a DLData object using default (DNA) alphabet and colors
#' # from a character vector with two entries
#' data <- DLData(c("ACGT", "ATTA"))
#' 
#' # creating a DLData object using a custom, binary alphabet and custom colors
#' data2 <- DLData(c("A,B,B,A,B", "A,B,B,A,A", "A,B,A,A,B"),
#'     symbols = c("A", "B"), colors = c("red","green"), delim = ",")
#' 
#' # creating a DLData object from a data frame 
#' # (created from a character vector, in this case)
#' vec <- c("A,B,B,A,B", "A,B,B,A,A", "A,B,A,A,B")
#' df <- as.data.frame(t(sapply(vec, function(a){strsplit(a, ",")[[1]]})))
#' data.df <- DLData(df, symbols = c("A", "B"), colors = c("red", "green"))
#' 
#' # creating a DLData object from sequences and weights, read from a tabular file
#' seqs <- read.table(system.file("extdata", "cjun.txt", package = "DepLogo"), 
#'     stringsAsFactors = FALSE)
#' data3 <- DLData(sequences = seqs[, 1], weights = log1p(seqs[, 2]) )
DLData<-function(sequences,
				 weights = NULL,
				 symbols = alphabet.dna$alphabet,
				 colors = alphabet.dna$colors,
				 delim = "",
				 sortByWeights = !is.null(weights),
				 axis.labels = NULL){
	sortByWeights <- sortByWeights
	if(class(sequences)=="character"){
		if( length( unique(sapply(sequences,nchar)) ) != 1 ){
			stop("All sequences must have the same length.");
		}
		
		seqs <- t(sapply(sequences, function(a){strsplit(a, delim)[[1]]}))
		rownames(seqs) <- NULL
	}else if(class(sequences)=="data.frame"){
		seqs <- as.matrix(sequences)
	}else{
		stop("Only data.frames or character vectors allowed.")
	}
	if( is.null(weights) ){
		weights <- rep(1, nrow(seqs))
	}
	syms <- unique(as.vector(seqs))
	if( sum(!syms%in%symbols)>0 ){
		stop("Provided symbols do not match those of the data.")
	}
	if(length(colors) != length(symbols)){
		stop("Number of colors not equal to number of symbols.")
	}
	
	df <- data.frame(seqs, weights)
	
	if(sortByWeights){
		o <- order(weights, decreasing = TRUE)
		df <- df[o, ]
	}
	
	li <- list(data = df, alphabet = Alphabet(symbols, colors), sortByWeights = sortByWeights, axis.labels = axis.labels)
	class(li) <- "DLData"
	li
}



summary.list <- function(object, delete.gaps = FALSE, ...){
	if(length(object)>0 & length(unique(sapply(object, class))) ==1 & class(object[[1]]) == "DLData"){
		df<-c();
		for(el in object){
			df<-rbind(df,data.frame(summary(el, delete.gaps = delete.gaps, ...), stringsAsFactors = FALSE), stringsAsFactors = FALSE)
		}
		df
	}else{
		NextMethod("summary")
	}
}


#' Summarizing DLData objects
#' 
#' \code{summary} method for class "DLData".
#' The summary includes the number of sequences, the consensus sequence
#' and the number of sequences in \code{object} that match the consensus.
#'
#' @param object an object of class "DLData"
#' @param delete.gaps if gaps should be removed from the consensus
#' @param ... further arguments passed to or from other methods
#'
#' @return a \code{list} with elements \code{members} containing the number of sequences,
#'   \code{consensus} containing the consensus sequences, and \code{equal.consensus} containing the
#'   number of sequences in \code{object} that are identical to \code{consensus}
#' @export
#' @author Jens Keilwagen, Jan Grau <grau@informatik.uni-halle.de>
#'
#' @examples
#' seqs <- read.table(system.file("extdata", "cjun.txt", package = "DepLogo"), 
#'     stringsAsFactors = FALSE)
#' data <- DLData(sequences = seqs[, 1], weights = log1p(seqs[, 2]) )
#' summary(data)
summary.DLData <- function(object, delete.gaps = FALSE, ...){
	m <- object$data
	rn <- names(m)
	m <- m[, rn!="weights"]
	
	consensus <- apply(m, 2, function(x){ 
		tab <- table(x)
		names(tab)[which.max(tab)]
	} )
	
	num <- nrow(m);
	
	idx <- 1:num;
	i <- 1;
	r <- ncol(m);
	while( length(idx)>0 && i <= r ) {
		idx <- idx[ m[idx, i]==consensus[i] ]
		i <- i+1;
	}
	
	if( delete.gaps ) {
		consensus <- consensus[ consensus != "-" ]
	}             
	seq <- paste(consensus, collapse="")
	
	return( list( members = num, consensus = seq, equals.consensus = length(idx) ) )
}




#' Filters columns (sequence positions) by gaps
#'
#' @param percent.gap the maximum fraction of gaps allowed to retain a column
#'
#' @return function that, given a \link{DLData} object, returns \code{TRUE} 
#'   for every column that does not exceed the specified number of gaps
#' @export
#' @author Jan Grau <grau@informatik.uni-halle.de>
#'
#' @examples
#' fun <- filter.by.gaps(percent.gap = 0.1)
filter.by.gaps <- function(percent.gap){
	function(data){
		rats <- apply(data$data[, -ncol(data$data)], 2, function(a){sum(a=="-")})/nrow(data$data)
		sel <- (rats<percent.gap)
		list(selected = sel, range = range(rats))
	}
}


#' Filters columns (sequence positions) by conservation
#' 
#' Filters columns based on the relative information content of each column
#' which is the standard information content normalized to the interval [0,1],
#' where 0 corresponds to uniform distribution and 1 to perfect conservation
#' of one nucleotide or amino acid, respectively.
#'
#' @param relative.ic the maximum relative information content allowed
#'   to retain a position
#'
#' @return function that, given a \link{DLData} object, returns \code{TRUE} 
#'   for every column that does not exceed the specified relative information content
#' @export
#' @author Jan Grau <grau@informatik.uni-halle.de>
#'
#' @examples
#' fun <- filter.by.conservation(relative.ic = 0.9)
filter.by.conservation <- function(relative.ic){
	function(data){
		pwm <- getPWM(data)
		ics <- apply(pwm, 2, function(a){getICScale(a)^2})
		sel <- (ics<relative.ic)
		list(selected = sel, range = range(ics, na.rm = TRUE))
	}
}

#' Filters columns (sequence positions) by dependency
#' 
#' Filters columns based on the average or maximum mutual information of a column
#' to all other columns. Mutual information is normalized to to interval
#' [0,1], where 0 corresponds to independence and 1 to perfect dependence.
#'
#' @param mi.threshold the minimum average or maximum mutual information required
#' @param use.max if \code{TRUE}, the maximum and otherwise the average mutual
#'   information will be considered
#'
#' @return function that, given a \link{DLData} object, returns \code{TRUE} 
#'   for every column that does exceed the specified average mutual information
#' @export
#' @author Jan Grau <grau@informatik.uni-halle.de>
#'
#' @examples
#' fun <- filter.by.dependencies(mi.threshold = 0.3)
filter.by.dependencies <- function(mi.threshold, use.max = FALSE){
	function(data){
		deps <- getDeps(data)/nrow(data$data)/log(length(data$alphabet))
		if(use.max){
			deps <- apply(deps, 1, max, na.rm = TRUE)	
		}else{
			deps <- rowSums(deps)/(ncol(data$data) - 1)
		}
		sel <- deps>=mi.threshold
		list(selected = sel, range = range(deps, na.rm = TRUE))
	}
}


#' Filters data columns by some filter function
#'
#' Filters the columns of the input data, i.e., positions of input sequences,
#' by a filter function that, given a \link{DLData} object, returns a list containing
#' i) as element \code{$selected} a vector with entries \code{TRUE} for every column 
#' that should be retained in the filtered data and ii) as element \code{$range} the
#' range of values obtained for the filtering criterion.
#'
#' @param data the data as \link{DLData} object
#' @param filter.fun the filter function
#'
#' @return a \link{DLData} object containing the filtered columns and the indexes of the remaining in its \code{axis.labels} field
#' @export
#' @author Jan Grau <grau@informatik.uni-halle.de>
#' @seealso \link{filter.by.gaps}
#' @seealso \link{filter.by.dependencies}
#' @seealso \link{filter.by.conservation}
#'
#' @examples
#' # read data and create DLData object
#' seqs <- read.table(system.file("extdata", "cjun.txt", package = "DepLogo"), 
#'     stringsAsFactors = FALSE)
#' data <- DLData(sequences = seqs[, 1], weights = log1p(seqs[, 2]) )
#' 
#' # create a filter function based on the percentage of gap symbols (at most 10%)
#' fun <- filter.by.gaps(percent.gap = 0.1)
#' data2 <- filterColumns(data, fun)
filterColumns <- function(data, filter.fun){
	UseMethod("filterColumns", data)
}

filterColumns.DLData<-function(data, filter.fun){
	ax.lab <- data$axis.labels
	if(is.null(ax.lab)){
		ax.lab <- 1:(ncol(data$data) - 1)
	}
	
	temp <- filter.fun(data)
	sel <- temp$selected
	if(sum(sel)==0){
		stop(paste0("No column matched the filter criteria. Range of values: [", temp$range[1], ", ", temp$range[2], "]"))
	}
	dat2 <- data
	dat2$data <- dat2$data[, c(sel, TRUE)]
	dat2$axis.labels <- ax.lab[sel]
	
	dat2
}





#' Suggests colors for symbols
#' 
#' Suggests colors for the symbols in \code{data} based on the co-occurrence of
#' symbols at common positions, weighted by the dependency values at those positions.
#' The idea is to assign similar colors only to symbols that either mostly occur at 
#' different positions or that are present at positions with low inter-dependencies 
#' to other positions.
#'
#' @param data the data
#'
#' @return the colors
#' @export
#' @author Jan Grau <grau@informatik.uni-halle.de>
#' @seealso \link{replaceColors}
#'
#' @examples
#' # read data and create DLData object
#' seqs <- read.table(system.file("extdata", "cjun.txt", package = "DepLogo"), 
#'     stringsAsFactors = FALSE)
#' data <- DLData(sequences = seqs[, 1] ,weights = log1p(seqs[, 2]) )
#' 
#' suggestColors(data)
suggestColors <- function(data){
	UseMethod("suggestColors", data)
}

suggestColors.DLData <- function(data){
	syms <- data$alphabet$chars
	deps <- getDeps(data)
	maxs <- apply(deps, 1, max)
	maxs <- maxs/sum(maxs)
	freqs <- apply(data$data[, -ncol(data$data)], 2, function(a){table(factor(a, levels = syms))})/nrow(data$data)
	
	d <- matrix(NA, nrow = length(syms), ncol = length(syms))
	rownames(d) <- colnames(d) <- syms
	for(a in syms){
		for(b in syms){
			d[a, b] <- sum( apply(freqs[c(a, b), ], 2, prod) * maxs )
		}
	}
	d <- matrix(rank(d), ncol = ncol(d))
	rownames(d) <- colnames(d) <- syms
	diag(d) <- 0

	scale <- cmdscale(d, k = 1, add = TRUE)$points[, 1]
	names(scale) <- syms
	ran <- rank(scale, ties.method = "min")
	colors <- rainbow(max(ran))[ran]
	if("-"%in%syms){
		colors[syms=="-"] <- "black"
	}
	colors
}


#' Replaces colors in \link{DLData} object
#'
#' @param data the data
#' @param colors the new colors
#'
#' @return the modified \link{DLData} object
#' @export
#' @author Jan Grau <grau@informatik.uni-halle.de>
#' @seealso \link{replaceColors}
#'
#' @examples
#' # read data and create DLData object
#' seqs <- read.table(system.file("extdata", "cjun.txt", package = "DepLogo"), 
#'     stringsAsFactors = FALSE)
#' data <- DLData(sequences = seqs[, 1], weights = log1p(seqs[, 2]) )
#' 
#' replaceColors(data, c("red", "green", "blue", "yellow"))
replaceColors <- function(data, colors){
	UseMethod("replaceColors", data)
}


replaceColors.DLData <- function(data, colors){
	if(length(colors) != length(data$alphabet$cols)){
		stop("Number of colors does not match the number of symbols.")
	}
	data$alphabet$cols <- colors
	data
}


#' Reverse complement
#'
#' Determine the reverse complementary \code{DLData} object. Only works for DNA or RNA. Data may include gap symbols.
#'
#' @param data the data
#'
#' @return the reverse complement
#' @export
#' @author Jan Grau <grau@informatik.uni-halle.de>
#'
#' @examples
#' data <- DLData(c("ACGT", "ATTA"))
#' revcom(data)
revcom<-function(data){
	UseMethod("revcom",data)
}

revcom.DLData <- function(data){
	alphabets<-list(
		dna = c("A" = "T", "C" = "G", "G" = "C", "T" = "A"),
		dnag = c("A" = "T", "C" = "G", "G" = "C", "T" = "A", "-" = "-"),
		rna = c("A" = "U", "C" = "G", "G" = "C", "U" = "A"),
		rnag = c("A" = "U", "C" = "G", "G" = "C", "U" = "A", "-" = "-"))
	
	data.alph <- data$alphabet$chars
	
	sel <- sapply(1:length(alphabets), function(i){
		ref <- names(alphabets[[i]])
		if(length(data.alph)!=length(ref)){
			FALSE
		}else{
			a <- sort(intersect(data.alph, ref))
			b <- sort(union(data.alph, ref))
			if(length(a)==length(b)){
				sum(a!=b)==0
			}else{
				FALSE
			}
		}
	})
	if(sum(sel)==0){
		stop("Alphabet is not DNA or RNA")
	}
	used <- alphabets[sel][[1]]
	
	df <- data$data[, -ncol(data$data)]
	df2 <- apply(df, 2, function(a){used[as.character(a)]})
	df2 <- df2[, ncol(df2):1]
	df2 <- data.frame(df2, weights=data$data[, ncol(data$data)])
	data$data <- df2
	
	data
}


