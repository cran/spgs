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

\name{chargaff1.test}
\alias{chargaff1.test}
\title{Matrix Test of CSPR for Mononucleotides}
\description{
Performs a test of Chargaff's second parity rule (CSPR) for mononucleotides on a genomic sequence using a 4X4 stochastic matrix 
estimated from the sequence.  The test was proposed by
Hart and Martínez (2011).
}
 \usage{
chargaff1.test(x, alg=c("table", "simulate", "upper"), n, no.p.value=FALSE)
}
\arguments{
\item{x}{either a vector containing the relative frequencies of each of the 4 nucleotides 
A, C, G, T, a character vector representing a DNA sequence in which each element contains a single nucleotide, or
a DNA sequence stored using the SeqFastadna class from the \pkg{\link[seqinr]{seqinr}} package.}
\item{alg}{the algorithm for computing the p-value.  If set to \dQuote{\samp{simulate}}, the
p-value is obtained via Monte Carlo simulation.  If set to \dQuote{\samp{upper}}, an analytic upper
bound on the p-value is computed.  \dQuote{\samp{upper}} are based on formulae
in Hart and Martínez (2011).  If
\samp{type} is specified as \dQuote{\samp{table}} (the default value),the p-value for the test is
obtained from a linear interpolation of a look-up table.  See the note below for further details.}
  \item{n}{The number of replications to use for Monte Carlo simulation.  
  If computationally feasible, a value >= 10000000 is recommended.}
\item{no.p.value}{If \samp{TRUE}, do not compute the p-value.  
The default is \samp{FALSE}.}
}
\details{
The first argument may be a character vector representing a DNA sequence, a DNA
sequence represented using the SeqFastadna class from the \pkg{seqinr} package,
or a vector containing the relative frequencies of the A, C, G and T nucleic
acids.

This function performs a test of Chargaff's second parity rule for
mononucleotides based on a 4X4 stochastic matrix \eqn{P} estimated from the
empirical dinucleotide distribution of a genomic sequence .  The \eqn{a,b)}
entry of \eqn{P} gives the empirical probability (relative frequency) that a
nucleotide \eqn{a} is followed by a nucleotide \eqn{b} in the sequence. The test
is set up as follows:

\eqn{H_0}{H0}:  the sequence (or matrix \eqn{P}) does not comply with CSPR for mononucleotides \cr
\eqn{H_1}{H1}:  the sequence (or matrix \eqn{P}) complies with CSPR for mononucleotides
}
\value{
A list with class "htest.ext" containing the following components:

\item{statistic}{the value of the test statistic.}
\item{p.value}{the p-value of the test. Only included if no.p.value is FALSE.}
\item{method}{a character string indicating what type of test was performed.}
\item{data.name}{a character string giving the name of the data.}
\item{f}{the 2-element vector used in calculating the test statistic.}
\item{estimate}{the stochastic matrix \eqn{P} used to derive the test statistic.}
\item{stat.desc}{a brief description of the test statistic.}
\item{null}{the null hypothesis (\eqn{H_0}{H0}) of the test.}
\item{alternative}{the alternative hypothesis (\eqn{H_1}{H1}) of the test.}
}
\note{
Currently, the look-up table that is employed when \samp{alg} is set to
\dQuote{\samp{table}} does not provide an accurate p-value when the statistic is
smaller than 0.00806.  Care should be taken when adjusting p-values for multiple testing.

The algebraically computed upper bound on the p-value obtained when \samp{alg}
is set to \dQuote{\samp{upper}} is not very tight and not suitable for real-
world applications.

\samp{no.p.value} suppresses computation of the p-value when it is set to
\samp{TRUE}.  This may be useful wen using this function to help simulate the
test statistic.
}
\references{
Hart, A.G. and Martínez, S. (2011)
Statistical testing of Chargaff's second parity rule in bacterial genome sequences.
\emph{Stoch. Models} \bold{27(2)}, 1--46.
}
\author{
Andrew Hart and Servet Martínez
}
\seealso{
\code{\link{chargaff0.test}}, \code{\link{chargaff2.test}}, 
\code{\link{agct.test}}, \code{\link{ag.test}},
\code{\link{chargaff.gibbs.test}}
}
\examples{
#Demonstration on real bacterial sequence
data(nanoarchaeum)
chargaff1.test(nanoarchaeum)

#Simulate synthetic DNA sequence that does not satisfy Chargaff's second parity rule
trans.mat <- matrix(c(.4, .1, .4, .1, .2, .1, .6, .1, .4, .1, .3, .2, .1, .2, .4, .3), 
ncol=4, byrow=TRUE)
seq <- simulateMarkovChain(500000, trans.mat, states=c("a", "c", "g", "t"))
chargaff1.test(seq)
}
\keyword{htest}
