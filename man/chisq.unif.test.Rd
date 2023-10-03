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

\name{chisq.unif.test}
\alias{chisq.unif.test}
\title{
Test of Uniformity Based on Pearsons's Chi-Squared test
}
\description{
Tests if a set of data points is uniformly distributed over a specified interval [a,b].
}
\usage{
chisq.unif.test(x, bins=NULL, interval=c(0,1), min.bin.size=10, all.inside=TRUE, 
rightmost.closed=TRUE, ...)
}
\arguments{
  \item{x}{A numeric vector of data values.}
  \item{bins}{If specified, the number of bins to use to discretise the
  interval.  Otherwise, the number of bins will be chosen automatically.}
  \item{interval}{A two-element vector giving the support of the uniform distribution.  The
  default is \code{c(0,1)}.}
  \item{min.bin.size}{The minimum number of data points to have in each bin.  If
  bins cannot be chosen without violating this constraint, an error is
  generated.  The default is 10.  This parameter is ignored if \samp{bins} 
  is specified.}
  \item{all.inside}{Determines if data points outside the interval should be counted as belonging to the extremal bins.  The default is \samp{TRUE}.}
  \item{rightmost.closed}{Determines if data points that coinside
  with\code{interval[2]} are counted as belonging to the last bin.  This
  parameter only has an effect if \code{all.inside} is set to \code{FALSE}.  The
  default is \code{TRUE}.}
  \item{\dots}{Additional parameters to be passed to \code{\link{chisq.test}}.}
}
\details{

This function tests the fit of a set of data points to a uniform distribution on
\eqn{[a,b]} by partitioning \eqn{[a,b]} into \samp{bins} bins, counting how many
points fall in each bin and then testing that the points are equally distributed
among the bins using Pearson's chi-squared test.

When \samp{bins} is not specified, its value is selected using the following
heuristic.  Let \eqn{n} be the number of data points.  If \eqn{n>200}, then
\samp{bins} is set to 20.  Otherwise, \samp{bins} is set to \eqn{n/10}.
Next,while there is a bin containing fewer than \samp{min.bin.size} points in
the resulting partition, \samp{bins} is decremented by one.  This process stops
either when \samp{bins} is equal to 1 or every bin contains at least
\samp{min.bin.size} points.
}
\section{Warning}{
The arguments \samp{all.inside} and \samp{rightmost.closed} are included 
for experimentation and should be altered with caution.
}
\value{
A list with class \dQuote{htest} containing the following components:

\item{statistic}{the value of the test statistic.}
\item{parameter}{A vector containing the degrees of freedom  of the chi-squared
test (df), the left end of the interval (a) and the right end of the interval
(b).}
\item{p.value}{the p-value of the test.}
\item{method}{a character string indicating what type of test was performed.}
\item{data.name}{a character string giving the name of the data.}
\item{bins}{The number of bins used for the test.}
\item{min.bin.size}{The minimum bin size.}
\item{interval}{The interval used for the test.}
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Andrew Hart and Servet Martínez
}
\seealso{
\code{\link{chisq.test}}, \code{\link{findInterval}}, \code{\link{ks.unif.test}}
}
\examples{
#Generate an IID uniform(0,1) sequence
u <- runif(1000)
chisq.unif.test(u)
}
\keyword{htest}
