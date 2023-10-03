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

\name{simulateMarkovChain}
\alias{simulateMarkovChain}
\title{
Simulate a first-Order Markov Chain
}
\description{
Simulates a first-order Markov chain.
}
\usage{
simulateMarkovChain(n, trans.mat, init.dist=NULL, states=colnames(trans.mat))
}
\arguments{
  \item{n}{the length of the sample path to simulate.}
  \item{trans.mat}{The transition matrix of the Markov chain to simulate.}
\item{init.dist}{
The initial distribution to use for starting the simulation.  If it is not
specified, the stationary distribution of the Markov chain will be computed from
trans.mat and used to start the simulation in the steady state.
}
\item{states}{This argument can be used to override the labels on the transition matrix (if any) and name the states output on the sample path.}
}
\details{
\samp{trans.mat} must be a stochastic matrix.  It must either have both row and
column names, in which case they must agree, or no row and column names at all.
The row/column names will be used to label the states visited by the Markov
chain in the sample path simulated.  If \samp{states} is specified, it will be
used to label the states of the Markov chain instead of the row/column names of
\samp{trans.mat}, in which the length of \samp{states} must agree with the
dimension of \samp{trans.mat}. If \samp{trans.mat} has no row/column names and
\samp{states} is not specified, then the states of the Markov chain will be
labelled \eqn{1,\ldots,n}, where \eqn{n} is the dimension of \samp{trans.mat}.
}
\value{
A vector of length n containing a realisation of the specified Markov chain.
}
\author{
Andrew Hart and Servet Martínez
}
\seealso{
\code{\link{estimateMarkovChain}}, \code{\link{rstochmat}}, \code{\link{rcspr2mat}}
}
\examples{
simulateMarkovChain(50, matrix(c(.8, .2, .2, .8), ncol=2))
simulateMarkovChain(50, rstochmat(3), states=c("yes", "no", "maybe"))
}
\keyword{models}
\keyword{ts}
\keyword{distribution}
\keyword{datagen}
