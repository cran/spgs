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

\name{rank.test}
\alias{rank.test}
\title{
The Rank Test of Statistical Independence
}
\description{
Test for a trend in a data series by comparing the number of increasing pairs in the series with the number expected in an i.i.d. series.
}
\usage{
rank.test(x)
}
\arguments{
  \item{x}{
a numeric vector or univariate time series.
}
}
\details{
Perform a test for trend based on the number of increasing ordered pairs in a data series.  
Consider pairs of the form (x(i), x(j)), where i<j.  An increasing pair is 
any such pair for which x_i<x_j.  This function counts the number of 
increasing pairs in the data, standardises it to have mean 0 and variance 
1 and asymptotically tests it against a standard normal distribution.
the test statistic is:

R = (pairs-mu)/sigma, where \cr
pairs is the number of increasing pairs in the data, \cr
mu = n*(n-1)/4, \cr
  sigma = sqrt(n*(n-1)*(2*n+5)/72) and \cr
  n is the number of data points in the series.

The test is set up as follows:

\eqn{H_0}{H0}:  the data series is i.i.d. (not trending) \cr
\eqn{H_1}{H1}:  the data series is not i.i.d. (trending)
}
\value{
A list with class "htest" containing the following components:

\item{statistic}{the value of the test statistic.}
\item{p.value}{the p-value of the test.}
\item{method}{a character string indicating what type of test was performed.}
\item{data.name}{a character string giving the name of the data.}
\item{pairs}{the number of increasing pairs counted in the data series.}
\item{n}{the number of points in the data series.}
\item{mu}{The expected number of increasing pairs that would be seen in an i.i.d. series.}
\item{sigma}{The standard deviation of the number of increasing pairs that would be seen in an i.i.d. series.}
}
\note{
Missing values are not handled. 

Points followed by a point having the exact same value are removed from the data series before computing the test statistic.

This test is useful for detecting linear trends in  data series.
}
\section{Warning}{
If the spgs shared object was successfully compiled with support for a 
64-bit unsigned integer type, then the following line should yield the value 0:

\kbd{rank.test(1:92683)$pairs-2^32-55607}

if not, then the package is only using 32-bit integer arithmetic for computing 
the rank test statistic and this will restrict \code{rank.test} to analysing 
series whose length is at most 92682.  In this case, attempting to apply 
\code{rank.test} to a series longer than 92682 will result in a warning about an 
integer overflow having occurred and the results of the test should \strong{not} 
be trusted.
}
\references{
Brockwell, Peter J., Davis, Richard A. (2002) \emph{Introduction to Time Series and Forecasting}.
Springer Texts in Statistics, Springer-Verlag, New York.
}
\author{
Andrew Hart and Servet Martínez
}
\seealso{
\code{\link{diffsign.test}}, \code{\link{turningpoint.test}}, \code{\link{lb.test}},
\code{\link{markov.test}}, \code{\link{diid.test}}
}
\examples{
#Generate an IID standard normal sequence
n <- rnorm(1000)
rank.test(n)
}
\keyword{htest}
