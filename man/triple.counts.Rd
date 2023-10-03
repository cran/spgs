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

\name{triple.counts}
\alias{triple.counts}
\title{
Count Triplets in Character Vectors
}
\description{
Count triples of adjacent symbols/elements in a character vector.
}
\usage{
triple.counts(x, case=c("lower", "upper", "as is"), circular=TRUE)
}
\arguments{
  \item{x}{
  a character vector or an object that can be coersed to a character vector.
  }
  \item{case}{
  determines how labels for the array should be generated: in
  'lower'  case, in ' upper'  case or 'as is', in which case labels such as 'b'
  and 'B' will be counted as distinct elements and counted separately.
  }
\item{circular}{
Determines if the vector should be treated as circular or not.  The default is 
\code{TRUE}, meaning that the start and end of the sequence will be joined 
together for the purpose of counting.
}
}
\details{
If \code{circular} is \code{TRUE}, the vector is treated as circular so that the 
some of all the counts in the resulting array is equal to the length of the 
vector and the sums across all dimentions of the array are equivalent, that is: 
if we write\cr
\code{t <- triple.counts(x)} \cr
for some character sequence x, then \cr
\code{apply(t,1,sum)}, \code{apply(t,2,sum)} and \code{apply(t,3,sum)} \cr
are all identical.

On the other hand, if \code{circular} is \code{FALSE}, the sum of all the 
entries in the counts array will be two less than the length of the vector and 
there will be a discrepancy between the sums over the various dimensions.
}
\value{
A 3-dimensional array of counts.  The labels of the \eqn{i}-th dimension correspond to
the \eqn{i}-th element of each triple, where \eqn{i} is either 1, 2 or 3.
}
\author{
Andrew Hart and Servet Martínez
}
\seealso{
\code{\link{pair.counts}}, \code{\link{quadruple.counts}},
\code{\link{cylinder.counts}}, 
\code{\link{array2vector}}, \code{\link{table2vector}}
}
\keyword{array}
