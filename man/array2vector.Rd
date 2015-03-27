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

\name{array2vector}
\alias{array2vector}
\alias{table2vector}
\title{
Convert Arrays and Tables to Vectors with Named Elements
}
\description{
Convert an array/table into an integer vector, preserving the names 
corresponding to each element in a sensible way.  These functions differ from 
\code{\link{as.vector}} in that they name each element of the resulting vector 
by combining appropriate names from the various dimensions that together 
uniquely identify the position of each element in the original array/table.
}
\usage{
array2vector(x, sep=".", sort=FALSE, rev=FALSE, ...)
table2vector(x, sep=".", sort=FALSE, rev=FALSE, ...)
}
\arguments{
  \item{x}{
  an array or table.  
 }
  \item{sep}{
    a character string to be used to separate the name corresponding to each 
  dimension when constructing the element names for the vector.  The default value is \code{"."}.
    }
\item{sort}{
Should the elements in the resulting vector be sorted in lexicographic order 
according to the names they are assigned?  The default is \code{FALSE}.
}
\item{rev}{
For the purposes of sorting, should the names of the vector's elements be read 
from right-to-left, i.e. in reverse order?  The default is \code{FALSE}.
}
\item{\dots}{
Arguments to be passed from or to other functions.
}
}
\details{
\code{table2vector} is merely a convenience alias for \code{array2vector}, which 
converts a multi-dimensional array or table to a vector using 
\code{\link{as.vector}}, but names each of the elements in the resulting vector 
according to the names contained in its \code{dimnames} attribute.

the name of each element is constructed by concatenating names (one from each 
\code{dimnames} member) seperated by the value specified in \code{sep}.

Note that dimensions of \code{x} which lack a corresponding vector of names in 
the \code{dimnames} attribute will be assigned a names vector of the form 
\code{1:d} where \code{d} is the dimension size specified in the corresponding 
entry of the \code{dim} attribute.
}
\value{
An integer vector with names set as described in \sQuote{Details}.
}
\author{
Andrew Hart and Servet Martínez
}
\seealso{
\code{\link{pair.counts}}, \code{\link{triple.counts}}, 
\code{\link{quadruple.counts}}, \code{\link{cylinder.counts}}
}
\examples{
a <- array(1:8, dim=rep(2,3), dimnames=list(c("a","b"), c("x","p"), c("v","u")))
array2vector(a)
array2vector(a, sep="")
array2vector(a, sep="", sort=TRUE)
array2vector(a, sep="", sort=TRUE, rev=TRUE)
array2vector(a, sep="", sort=TRUE, decreasing=TRUE)
}
\keyword{array}
