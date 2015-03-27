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

\name{lb.test}
\alias{lb.test}
\title{
The Ljung-Box Test for Uncorrelated Data
}
\description{
This function is a convenient wrapper for using \code{Box.test} to perform the Ljung-
Box Q test of uncorrelated data without having to specify \samp{type}.  In other
words, \code{lb.test(x, ...)} is equivalent to 
\code{Box.test(x, type="Ljung-Box", ...)}.
}
\usage{
lb.test(x, ...)
}
\arguments{
  \item{x}{
  a numeric vector or univariate time series.
}
\item{\dots}{
parameters to pass to \code{Box.test}.
}
}
\value{
A list with class "htest" containing the following components:

\item{statistic}{the value of the test statistic.}
\item{parameter}{the degrees of freedom of the approximate chi-squared distribution of the test statistic (taking \code{fitdf} into account).}
\item{p.value}{the p-value of the test.}
\item{method}{a character string indicating what type of test was performed.}
\item{data.name}{a character string giving the name of the data.}
}
\author{
Andrew Hart and Servet Martínez
}
\seealso{
\code{\link{Box.test}}, \code{\link{markov.test}}, \code{\link{diid.test}}
\code{\link{diffsign.test}}, \code{\link{turningpoint.test}}, \code{\link{rank.test}}
}
\examples{
#Generate an IID standard normal sequence
n <- rnorm(1000)
lb.test(n)
}
\keyword{htest}
