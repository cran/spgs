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

\name{estimateMarkovChain}
\alias{estimateMarkovChain}
\title{
Fit a first-Order Markov Chain to a Sequence of Finite Symbols
}
\description{
Estimates the transition matrix and stationary distribution of  a first-order
Markov chain from an observed sequence of symbols.
}
\usage{
estimateMarkovChain(x, circular=TRUE)
}
\arguments{
  \item{x}{The sequence of observed symbols as a character vector.}
  purposes  \item{circular}{Should the sequence be treated as circular for the purpose of estimation?  The default is \samp{TRUE}.}
}
\value{
A list with class \samp{FiniteStateMarkovChain} having the following components:

\item{trans.mat}{The stochastic transition matrix estimated from x.}
\item{stat.dist}{The stationary distribution estimated from x.}
\item{states}{the state labels}
}
\author{
Andrew Hart and Servet Martínez
}
\seealso{
\code{\link{markov.test}}, \code{\link{markov.disturbance}}, \code{\link{simulateMarkovChain}}
}
\examples{
#Obtain a random 3 x 3 stochastic matrix with rows and columns labelled "A", "B", "C"
mat <- rstochmat(3, labels=c("A", "B", "C"))
mat

#Simulate a Markov chain of length 500 using mat as the transition matrix
seq <- simulateMarkovChain(500, mat)

#Estimate mat and the stationary distribution for the Markov chain which generated seq
estimateMarkovChain(seq)
}
\keyword{models}
