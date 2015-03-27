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

\name{rstochvec}
\alias{rstochvec}
\title{
Random Generation of Stochastic (Probability) Vectors
}
\description{
Randomly generate probability vectors, that is, non-negative vectors whose 
elements sum to unity.
}
\usage{
rstochvec(n, labels)
}
\arguments{
  \item{n}{the length of the vector.  If n is not specified, it is inferred from the lenth of \samp{labels}.}
  \item{labels}{
  a vector of labels for the elements of the vector.  If \samp{labels} is not 
  specified, n must be specified and the value \code{as.character(1:n)} is assumed.
  }
}
\details{
Stochastic (or probability) vectors are non-negative vectors that sum to unity.
This function uniformly generates samples from the set of probability vector sof
length \eqn{n}.

At least one of the arguments must be specified.  The missing argument is infered from the other.
}
\value{
A probability vector of length \eqn{n} with elements named according to \samp{labels}.
}
\author{
Andrew Hart and Servet Martínez
}
%\seealso{
%}
\examples{
rstochvec(4)
rstochvec(3, c("a", "b", "c"))
rstochvec(labels=c("r", "R"))
}
\keyword{distribution}
\keyword{datagen}

