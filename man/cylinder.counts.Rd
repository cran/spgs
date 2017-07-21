%This file is part of the source code for
%SPGS: an R package for identifying statistical patterns in genomic sequences.
%Copyright (C) 2015  Universidad de Chile and INRIA-Chile
%
%This program is free software; you can redistribute it and/or modify
%it under the terms of the GNU General Public License as published by
%the Free Software Foundation; either version 2 of the License, or
%(at your option) any later version.
%
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.
%
%A copy of Version 2 of the GNU Public License is available in the 
%share/licenses/gpl-2 file in the R installation directory or from 
%http://www.R-project.org/Licenses/GPL-2.

\name{cylinder.counts}
\alias{cylinder.counts}
\title{
Count Cylinders (Fixed-Offset Patterns) in Character Vectors
}
\description{
Count fixed tuples of not necessarily adjacent symbols/elements in a character vector.
}
\usage{
cylinder.counts(x, cylinder, case=c("lower", "upper", "as is"), circular=TRUE)
}
\arguments{
  \item{x}{
  a character vector or an object that can be coersed to a character vector.
  }
  \item{cylinder}{
  A vector of indices specifying the form of cylinders to count.  See \sQuote{Details}.
  }
\item{case}{
  determines how labels for the array should be generated: in
  lowercase, in uppercase or left as is, in which case labels such as \dQuote{b}
  and \dQuote{B} will be seen as distinct symbols and counted separately.
}
\item{circular}{
Determines if the vector should be treated as circular or not.  The default is 
\code{TRUE}, meaning that the start and end of the sequence will be joined 
together for the purpose of counting.
}
}
\details{
\code{cylinder} represents a set of symbol patterns that one wishes to count in 
the sequence \code{x}. For example, if \code{cylinder} is \code{c(1,3,5)}, then 
this function will count occurrences of all patterns of the form \samp{u.v.w}, 
where \samp{u}, \samp{v} and \samp{w} can be any symbol present in \code{x} and 
\code{.} stands for a symbol whose value is not relevant to the pattern.    

Suppose that \code{x} is a sequence of the nucleotides \code{a}, \code{c}, 
\code{g} and \code{t}.  Then, \code{cylinder=1:2} will count the occurrences of 
all 16 dinucleotides: \code{aa}, \code{ac}, \code{ag}, \code{at}, \code{ca}, 
\code{cc}, \dots.  In contrast, \code{cylinder=c(1,3)} will counts 16 sets of 
trinucleotides:  \code{a.a}, \code{a.c}, \code{a.g}, \code{a.t}, \code{c.a}, 
\code{c.c}, \code{c.g}, \dots.  the dot \dQuote{\code{.}} stands for any 
nucleotide, so that \code{a.c} represents the set {\code{aac}, \code{acc}, 
\code{agc}, \code{atg}}.  In both of these examples, a \eqn{4\times 4}{4 X 4} 
array of counts will be produced, but in the first case the array will 
represent counts of dinucleotides, while in the second case it will represent 
counts of groups of trinucleotides.

If \code{circular} is \code{TRUE}, the vector \code{x} is treated as circular so that the 
some of all the counts in the resulting array is equal to the length of the 
vector and the sums across all dimentions of the array are equivalent, that is: 
writing\cr
\code{counts <- cylinder.counts(x, cylinder=c(1,3,5))} \cr
for some character sequence x, then \cr
\code{apply(counts,1,sum)}, \code{apply(counts,2,sum)} and \code{apply(counts,3,sum)}\cr
will all be identical.

On the other hand, if \code{circular} is \code{FALSE}, the sum of all the 
entries in the counts array will be less than the length of the vector and 
there will be a discrepancy between the sums over the various dimensions.
}
\value{
An \eqn{n}-dimensional array of counts, where \eqn{n} is the length of 
\code{cylinder}.
}
\note{
\code{\link{table}}is more efficient (by almost a factor of 2) at computing the 
counts of cylinders of length 1, whereas \code{cylinder.counts} is faster and 
uses less memory than \code{} for cylinders of length greater than 1.  
}
\author{
Andrew Hart and Servet Martínez
}
\seealso{
\code{\link{pair.counts}}, \code{\link{triple.counts}}, 
\code{\link{quadruple.counts}}, 
\code{\link{array2vector}}, \code{\link{table2vector}}
}
\examples{
#Generate an IID uniform DNA sequence
seq <- simulateMarkovChain(5000, matrix(0.25, 4, 4), states=c("a","c","g","t"))
cylinder.counts(seq, 1) #essentially the same as unclass(table(seq))
cylinder.counts(seq, 1:5) #counts of all 5-mers in the sequence

 #counts of all patterns of the form a.b where a and b represent
 #specific symbols and . denotes an arbitrary symbol.
 pat <- cylinder.counts(seq, c(1, 3))
#For example, pat["a","c"] gives the number of times that any of 
#the following 4 words appears in the sequence:  aac, acc, agc, atc.
identical(cylinder.counts(seq, c(1,3)), apply(cylinder.counts(seq, 1:3), c(1, 3), sum))

##some relationships between cylinder.counts and other functionns
identical(cylinder.counts(seq, 1:2), pair.counts(seq))
identical(cylinder.counts(seq, 1:3), triple.counts(seq))
identical(cylinder.counts(seq, 1:4), quadruple.counts(seq))

#The following relationship means that counts on circular sequences are 
#invariant under translationn
identical(cylinder.counts(seq, 1:6), cylinder.counts(seq, 10:15))

#Treating seq as non circular, most of the preceding relationships continue to hold
identical(cylinder.counts(seq, 1:2, circular=FALSE), 
  pair.counts(seq, circular=FALSE))
identical(cylinder.counts(seq, 1:3, circular=FALSE), 
triple.counts(seq, circular=FALSE))
identical(cylinder.counts(seq, 1:4, circular=FALSE), 
  quadruple.counts(seq, circular=FALSE))
#The following relationship no longer holds; that is, non-circular counts
#are not invariant under translation.
identical(cylinder.counts(seq, 1:6, circular=FALSE), 
  cylinder.counts(seq, 10:15, circular=FALSE))
}
\keyword{array}
