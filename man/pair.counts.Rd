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

\name{pair.counts}
\alias{pair.counts}
\title{
Count Pairs in Character Vectors
}
\description{
Count pairs of adjacent symbols/elements in a character vector.
}
\usage{
pair.counts(x, case=c("lower", "upper", "as is"), circular=TRUE)
}
\arguments{
  \item{x}{
  a character vector or an object that can be coersed to a character vector.
}
\item{case}{
  determines how labels for the array should be generated: in
  'lower'  case, in ' upper'  case or 'as is', in which case labels such as 'b'
  and 'B' will be considered as distinct elements and counted separately.
}
\item{circular}{
Determines if the vector should be treated as circular or not.  The default is 
\code{TRUE}, meaning that the start and end of the sequence will be joined 
together for the purpose of counting.
}
}
\details{
When \code{circular} is \code{TRUE}, the vector is treated as circular so that 
the some of all the counts in the resulting matrix is equal to the length of the 
vector and the row and column sums are equivalent.  When \code{circular} is 
\code{FALSE}, the sum of all the entries in the counts matrix will be one less 
than the length of the vector and there will be a discrepancy between the row 
and column sums.
}
\value{
A matrix of counts.  The row and column labels correspond to
the first and second element of each pair, respectively.
}
\author{
Andrew Hart and Servet Martínez
}
\seealso{
\code{\link{triple.counts}}, \code{\link{quadruple.counts}},
\code{\link{cylinder.counts}}, 
\code{\link{array2vector}}, \code{\link{table2vector}}
}
\keyword{array}
