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

\name{diffsign.test}
\alias{diffsign.test}
\title{
the Differents-Sign Test of Statistical Independence
}
\description{
Tests for a trend in a data series by comparing the number of positive differences between successive elements in the series to the number expected in an i.i.d. series.
}
\usage{
diffsign.test(x)
}
\arguments{
  \item{x}{
  a numeric vector or univariate time series.
}
}
\details{
Perform a test for trend based on the signs of successive differences in a data series.
#this function counts the number of positive successive differences in the data, standardises #it to have mean 0 and variance 1 and asymptotically tests it against a 
standard normal distribution.  the test statistic is:

D = (pd - mu)/sigma, where \cr
pd is the number of positive differences in the data series, \cr
mu = (n-1)/2, \cr
sigma = sqrt((n+1)/12) and \cr
n is the number of points in the data series.

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
\item{n}{the number of points in the data series.}
\item{mu}{The expected number of positive differences that would be seen in an i.i.d. series.}
\item{sigma}{The standard deviation of the number of positive differences that would be seen in an i.i.d. series.}
}
\note{
Missing values are not handled. 

Points followed by a point having the exact same value are removed from the data series before computing the test statistic.

This test is useful for detecting linear trends in  data series.
}
\references{
Brockwell, Peter J., Davis, Richard A. (2002) \emph{Introduction to Time Series and Forecasting}.
Springer Texts in Statistics, Springer-Verlag, New York.
}
\author{
Andrew Hart and Servet Martínez
}
\seealso{
\code{\link{turningpoint.test}}, \code{\link{rank.test}}, \code{\link{lb.test}}
\code{\link{markov.test}}, \code{\link{diid.test}},
}
\examples{
#Generate an IID standard normal sequence
n <- rnorm(1000)
diffsign.test(n)
}
\keyword{htest}
