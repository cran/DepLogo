#' DNA alphabet
alphabet.dna <- list(alphabet=c("A", "C", "G", "T"), colors = c("green", "blue", "orange", "red"))
#' DNA alphabet with gaps
alphabet.dna.gap <- list(alphabet = c("A", "C", "G", "T", "-"), colors = c("green", "blue", "orange", "red", "black"))
#' RNA alphabet
alphabet.rna <- list(alphabet = c("A", "C", "G", "U"), colors = c("green", "blue", "orange", "red"))
#' RNA alphabet with gaps
alphabet.rna.gap <- list(alphabet = c("A", "C", "G", "U", "-"), colors = c("green", "blue", "orange", "red", "black"))
#' Amino acid alphabet
alphabet.protein <- list(alphabet = c("I", "L", "V", "F", "M", "C", "A", "G", "P", "T", "S", "Y", "W", "Q", "N", "H", "E", "D", "K", "R", "B", "Z", "X"), colors = rainbow(23))
#' Amino acid alphabet with gaps
alphabet.protein.gap <- list(alphabet = c("I", "L", "V", "F", "M", "C", "A", "G", "P", "T", "S", "Y", "W", "Q", "N", "H", "E", "D", "K", "R", "B", "Z", "X","-"), colors = rainbow(24))

# appends the letter which to the object letters
addLetter <- function (letters, letterPolygon, x.pos, y.pos, ht, wt, col = "black") 
{
	x <- x.pos + wt * letterPolygon$x
	y < y.pos + ht * letterPolygon$y
	polygons <- sum(is.na(x)) + 1  # a letter can consist of more then one polygon
	letters$x <- c(letters$x, NA, x)
	letters$y <- c(letters$y, NA, y)
	letters$col <- c(letters$col, rep(col, polygons))
	letters
}


# appends the letter which to the object letters
getLetter <- function (letterPolygon, x.pos, y.pos, ht, wt, col = "black") 
{
	x <- x.pos + wt * letterPolygon$x
	y <- y.pos + ht * letterPolygon$y
	polygons <- sum(is.na(x)) + 1  # a letter can consist of more then one polygon

	list(x = x, y = y, id = NULL, fill = NULL, col = rep(col, polygons))
}


##########
# Class to represent letterPolygons. A Letter consists of a x and a y vector
Letter <- function(x, y) {
	if(length(x) != length(y)) {
		stop("The length of vector x is different from length of vector y")
	}
	pts <- list(x = x, y = y)
	class(pts) <- "Letter"
	return(pts)
}

##' builts an object of class Alphabet from the given set of symbols and colors
##'
##' @title built alphabet
##' @param chars set of symbols
##' @param cols set of colors; one for each symbol
##' @return the Alphabet object
##' @export
##' @exportClass Alphabet
##' @author Martin Nettling
##' @examples
##' DNA <- Alphabet(c("A", "C", "G", "T"), c("green4", "blue", "orange", "red"))
Alphabet <- function(chars, cols) {
	obj <- list(chars = chars, cols = cols, size = length(chars))
	class(obj) <- "Alphabet"
	return(obj)
}

#################### Letter Polygons must be created
letterPolygons <- list()
class(letterPolygons) <- "LetterPolygons"
############## A
letterPolygons$A = Letter(
	c(0,  4,  6, 10,  8,  5, 2, 0, NA,2.2,2.6,7.4,7.8,2.2) * 0.1,
	c(0, 10, 10,  0,  0,7.5, 0, 0, NA,  3,  4,  4,  3,  3) * 0.1
)
############## P
tmpx.1=sin(seq(0, 1*pi, length = 80))*5 +5.0
tmpy.1=cos(seq(0, 1*pi, length = 80))*2.875 +7.125
tmpx.2=rev(sin(seq(0, 1*pi, length = 80))*3 +5.0)
tmpy.2=rev(cos(seq(0, 1*pi, length = 80))*1.375 +7.125)
letterPolygons$P = Letter(
	c(0, 2, 2,0,NA,1, 5, tmpx.1,5,1,1,5,tmpx.2,5,1)*0.1,
	c(10,10,0,0,NA,10,10,tmpy.1,4.25,4.25,5.75,5.75,tmpy.2,8.5,8.5)*0.1
)
############## B
letterPolygons$B = Letter(
	c(0, 2, 2,0,NA,1, 5, 5,  1,  NA,1,   5,   5,   1,   NA,1,5,5,  1,  NA,tmpx.1,tmpx.2,NA,tmpx.1,tmpx.2)*0.1,
	c(10,10,0,0,NA,10,10,8.5,8.5,NA,5.75,5.75,4.25,4.25,NA,0,0,1.5,1.5,NA,tmpy.1,tmpy.2,NA,tmpy.1-4.25,tmpy.2-4.25)*0.1
)
############## R
letterPolygons$R = Letter(
	c(0, 2, 2,0,NA,1, 5, 5,  1,  NA,1,   5,   5,   1,   NA,5, 10, 8, 3,NA,tmpx.1,tmpx.2)*0.1,
	c(10,10,0,0,NA,10,10,8.5,8.5,NA,5.75,5.75,4.25,4.25,NA,5, 0 , 0, 5,NA,tmpy.1,tmpy.2)*0.1
)
############## C is copied from SeqLogo
angle1 = seq(0.3 + pi/2, pi, length = 40)
angle2 = seq(pi, 1.5 * pi, length = 40)
x.l1 = 0.5 + 0.5 * sin(angle1)
y.l1 = 0.5 + 0.5 * cos(angle1)
x.l2 = 0.5 + 0.5 * sin(angle2)
y.l2 = 0.5 + 0.5 * cos(angle2)
x.l = c(x.l1, x.l2)
y.l = c(y.l1, y.l2)
x = c(x.l, rev(x.l))
y = c(y.l, 1 - rev(y.l))
x.i1 = 0.5 + 0.35 * sin(angle1)
y.i1 = 0.5 + 0.35 * cos(angle1)
x.i1 = x.i1[y.i1 <= max(y.l1)]
y.i1 = y.i1[y.i1 <= max(y.l1)]
y.i1[1] = max(y.l1)
x.i2 = 0.5 + 0.35 * sin(angle2)
y.i2 = 0.5 + 0.35 * cos(angle2)
x.i = c(x.i1, x.i2)
y.i = c(y.i1, y.i2)
letterPolygons$C = Letter(
	c(x, rev(c(x.i, rev(x.i))))/max(x),
	c(y, rev(c(y.i, 1 - rev(y.i))))
)
############## G is copied from SeqLogo
letterPolygons$G = Letter(
	c(letterPolygons$C$x, NA, 1,  0.5,0.5,0.8,0.8,1,1),
	c(letterPolygons$C$y, NA, 0.4,0.4,0.3,0.3,0,  0,0.4)
)
############## D
tmpx.1=sin(seq(0, 1*pi, length = 80))*5 +5.0
tmpy.1=cos(seq(0, 1*pi, length = 80))*5 +5
tmpx.2=rev(sin(seq(0, 1*pi, length = 40))*3.5 +5.0)
tmpy.2=rev(cos(seq(0, 1*pi, length = 40))*3.5 +5.0)
letterPolygons$D = Letter(
	c(0, 2, 2,0,NA,1, 5, 5,  1,  NA,1,5,5, 1,   NA,tmpx.1,tmpx.2)*0.1,
	c(10,10,0,0,NA,10,10,8.5,8.5,NA,0,0,1.5,1.5,NA,tmpy.1,tmpy.2)*0.1
)
############## E
letterPolygons$E = Letter(
	c(0,  10, 10,  2,   2, 0, 0,  NA, 1, 9, 9,   1,   NA, 1,10,10  ,1  ) * 0.1,
	c(10, 10, 8.5, 8.5, 0, 0, 10, NA, 4, 4, 5.5, 5.5, NA, 0,0,1.5,1.5) * 0.1
)
############## F
letterPolygons$F = Letter(
	c(0,  10, 10,  2,   2, 0, 0,  NA, 1, 8, 8, 1 ) * 0.1,
	c(10, 10, 8.5, 8.5, 0, 0, 10, NA, 4, 4, 5.5, 5.5 ) * 0.1
)
############## H
letterPolygons$H = Letter(
	c(0,  2,  2, 0, NA, 8, 10, 10, 8, NA,1,9,9,1) * 0.1,
	c(10, 10, 0, 0, NA, 10, 10, 0, 0, NA,4,4,6,6) * 0.1
)
############## H
letterPolygons$I = Letter(
	c(3,  7,  7, 3) * 0.1,
	c(10, 10, 0, 0) * 0.1
)
############## J
tmpx.1=sin(seq(0.5*pi, 1.5*pi, length = 40))*3 +5
tmpy.1=cos(seq(0.5*pi, 1.5*pi, length = 40))*2.5 +2.5
tmpx.2=rev(sin(seq(.5*pi, 1.5*pi, length = 40))*1 +5)
tmpy.2=rev(cos(seq(.5*pi, 1.5*pi, length = 40))*1 +2.5)
letterPolygons$J = Letter(
	c(6,  8 , 8, 6,NA,tmpx.1,tmpx.2) * 0.1,
	c(10, 10, 2.5, 2.5,NA,tmpy.1,tmpy.2) * 0.1
)
############## K
letterPolygons$K = Letter(
	c(0,  2,  2, 0, NA, 0.2,8, 10, 0.2, NA, 5, 10, 8, 3) * 0.1,
	c(10, 10, 0, 0, NA, 4,10,10,2, NA, 6, 0 , 0, 6) * 0.1
)
############## L
letterPolygons$L = Letter(
	c(0,  2,  2,   10,  10, 0) * 0.1,
	c(10, 10, 1.5, 1.5, 0,  0 ) * 0.1
)

############## gap
letterPolygons$`-` = Letter(
	c(0,  10,  10,  0) * 0.1,
	c(10, 10,   0,  0) * 0.1
)


############## M
letterPolygons$M = Letter(
	c(0,  2,  2, 0, NA, 8, 10, 10, 8, NA, 1.5,4,6.0,3.5 , NA, 8.5,6,4.,6.5) * 0.1,
	c(10, 10, 0, 0, NA, 10, 10, 0, 0, NA, 10 ,1,1,  10, NA, 10, 1,  1,10) * 0.1
)
############## N
letterPolygons$N = Letter(
	c(0,  2,  2, 0, NA, 8, 10, 10, 8, NA, .5,2.5, 9.5,7.5 ) * 0.1,
	c(10, 10, 0, 0, NA, 10, 10, 0, 0, NA, 10, 10,  0,  0) * 0.1
)
############## O
tmpx=sin(seq(0, 2*pi, length = 80))/2+0.5
tmpy=cos(seq(0, 2*pi, length = 80))/2+0.5
letterPolygons$O = Letter(
	c(tmpx,rev(tmpx)*0.6+0.2),
	c(tmpy,tmpy*0.7+0.15)
)
############## Q
letterPolygons$Q = Letter(
	c(letterPolygons$O$x,NA,.5,.7,1.0,.8),
	c(letterPolygons$O$y,NA,.3,.3,.0,.0)
)
############## S
tmpx.1=(sin(seq(0.0*pi, 1.5*pi, length = 100))*5 +5)
tmpy.1=cos(seq(0.0*pi, 1.5*pi, length = 100))*2.825 +2.825
tmpx.2=rev(sin(seq(.0*pi, 1.5*pi, length = 100)))*3 +5
tmpy.2=rev(cos(seq(.0*pi, 1.5*pi, length = 100))*1.5 +2.825)

letterPolygons$S = Letter(
	c(tmpx.1,tmpx.2,NA,-tmpx.1+10,-tmpx.2+10) * 0.1,
	c(tmpy.1,tmpy.2,NA,-tmpy.1+10,-tmpy.2+10) * 0.1
)

############## T
letterPolygons$T = Letter(
	c(0, 10, 10, 6, 6, 4, 4, 0) * 0.1,
	c(10, 10, 9, 9, 0, 0, 9, 9) * 0.1
)
############## U
tmpx.1=rev(sin(seq(0.5*pi, 1.5*pi, length = 80))*5 +5)
tmpy.1=cos(seq(0.5*pi, 1.5*pi, length = 80))*3 +3
tmpx.2=(sin(seq(.5*pi, 1.5*pi, length = 80))*3 +5)
tmpy.2=rev(cos(seq(.5*pi, 1.5*pi, length = 80))*1.75 +3)
letterPolygons$U = Letter(
	c(0,  0,tmpx.1,10, 10,8 ,8,tmpx.2,2,2 ) * 0.1,
	c(10, 3,tmpy.1,3, 10,10,3,tmpy.2,3,10) * 0.1
)
############## V
letterPolygons$V = Letter(
	c(0,4,6,10,8,5,2) * 0.1,
	c(10,0,0,10,10,2,10) * 0.1
)
############## W
letterPolygons$W = Letter(
	c(0, 2,4,5.5, 4, 3.0,1.5, NA,4.5, 6,8,10, 8.5,7,6) * 0.1,
	c(10,0,0,10, 10,  2,  10, NA,10,0,0,10,10,2,10) * 0.1
)
############## X
letterPolygons$X = Letter(
	c(0,2,10,8,NA,0,2,10,8) * 0.1,
	c(10,10,0,0,NA,0,0,10,10) * 0.1
)
############## Y
letterPolygons$Y = Letter(
	c(0,2,6,4,NA,4,6,10,8,NA,4,6,6,4) * 0.1,
	c(10,10,4.5,4.5,NA,4.5,4.5,10,10,NA,5,5,0,0) * 0.1
)
############## Z
letterPolygons$Z = Letter(
	c(0,2.5,10,7.5,NA,0,10,10,0,NA,0,10,10,0) * 0.1,
	c(1.5,1.5,8.5,8.5,NA,10,10,8.5,8.5,NA,0,0,1.5,1.5) * 0.1
)


informationContent<-function (p) 
{
	obj = list()
	ic = log2(length(p)) + sum(p * log2(p), na.rm = T)
	obj$height = ic
	obj$ylab = "Information Content [bits]"
	return(obj)
}



plot.minlogo<-function(pwm,sparse,alphabet){
	letters = list(x = NULL, y = NULL, id = NULL, fill = NULL)
	npos = ncol(pwm)
	ylim.negMax = 0
	ylim.posMax = 0
	wt = 1
	x.pos = 0.5
	eps = 0
	heights = c()
	ymins = c()
	ymaxs = c()
	for (j in 1:npos) {
		column = pwm[, j]
		sh = informationContent(column)
		hts = column * sh$height
		letterOrder = order(abs(hts))
		yneg.pos = 0
		ypos.pos = 0
		for (i in 1:alphabet$size) {
			ht = hts[letterOrder[i]]
			y.pos = ypos.pos
			ht = ht - eps
			ypos.pos = ypos.pos + ht + eps
			char = alphabet$chars[letterOrder[i]]
			col = alphabet$cols[letterOrder[i]]
			letters = addLetter(letters, letterPolygons[[char]], 
										   x.pos, y.pos, ht, wt * 0.99, col = col)
		}
		x.pos = x.pos + wt
	}
	if(sparse){
		plot(NA,xlim = c(0.5, x.pos), ylim = c(0, log2(alphabet$size)),
			 axes=F,frame.plot=F,xlab="",ylab="",yaxs='i',xaxs="i")
		axis(2,at = 0:2,labels=0:2,pos=0.5,cex=par("cex"))
	}else{
		plot(NA,xlim = c(0.5, x.pos), ylim = c(0, log2(alphabet$size)),
			 axes=F,frame.plot=F,xlab="",ylab="",xaxs="i")
		axis(2,at = 0:2,labels = 0:2,pos=0.5,cex=par("cex"))
		axis(1,at=1:ncol(pwm),cex=par("cex"))
	}
	
	polygon(letters, col = letters$col, border = NA)
}
