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

\name{rstochmat}
\alias{rstochmat}
\title{
Random Generation of Stochastic Matrices
}
\description{
Randomly generates stochastic matrices.
}
\usage{
rstochmat(n, labels)
}
\arguments{
  \item{n}{the dimension of the matrix.  If n is not specified, it is inferred from the lenth of \samp{labels}.}
  \item{labels}{
  a vector of labels for the rows and columns of the matrix.  If \samp{labels}
  is not specified, \samp{n} must be specified and the value
  \code{as.character(1:n)} is assumed.
}
}
\details{
Stochastic matrices are non-negative matrices whose rows all sum to unity.  This
function uniformly generates samples from the set of \eqn{n\times n}{n X n} stochastic matrices.

At least one of the arguments must be specified.  The missing argument is infered from the other.
}
\value{
An \eqn{n\times n}{n X n} stochastic matrix with rows and columns labelled according to \samp{labels}.
}
\author{
Andrew Hart and Servet Martínez
}
\seealso{
\code{\link{rcspr2mat}}, \code{\link{estimateMarkovChain}}, \code{\link{simulateMarkovChain}}
}
\examples{
rstochmat(4)
rstochmat(3, c("a", "b", "c"))
rstochmat(labels=c("r", "R"))
}
\keyword{array}
\keyword{distribution}
\keyword{datagen}

