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

\name{ks.unif.test}
\alias{ks.unif.test}
\title{
Using \code{ks.test} to test for Uniformity on the Unit Interval
}
\description{
Uses \code{ks.test} to test that a data vector is uniform on the unit interval.
\code{ks.unif.test(x)} is merely convenient shorthand for 
\code{ks.test(x,
punif)}.
}
\usage{
ks.unif.test(x)
}
\arguments{
  \item{x}{
  a numeric vector or univariate time series.
}
}
\value{
A list with class "htest" containing the following components:

\item{statistic}{the value of the test statistic.}
\item{p.value}{the p-value of the test.}
\item{method}{a character string indicating what type of test was performed.}
\item{data.name}{a character string giving the name of the data.}
}
\author{
Andrew Hart and Servet Martínez
}
\seealso{
\code{\link{chisq.unif.test}}, \code{\link{markov.test}}, \code{\link{diid.test}}
}
\examples{
#Generate an IID uniform(0,1) sequence
u <- runif(1000)
ks.unif.test(u)
}
\keyword{htest}
