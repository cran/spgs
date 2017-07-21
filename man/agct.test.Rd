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

\name{agct.test}
\alias{agct.test}
\title{Test of Purine-Pyrimidine Parity Based on Euclidean distance}
\description{
Performs a test proposed by Hart and Martínez (2011) for the equivalence of the
relative frequencies of purines (\eqn{A+G}) and pyrimidines (\eqn{C+T}) in DNA
sequences.  It does this by checking whether or not the mononucleotide
frequencies of a DNA sequence satisfy the relationship A+G=C+T.
}
 \usage{
agct.test(x, alg=c("exact", "simulate", "lower", "Lower", "upper"), n)
}
\arguments{
\item{x}{either a vector containing the relative frequencies of each of the 4 nucleotides 
A, C, G, T, a character vector representing a DNA sequence in which each element contains a single nucleotide, or
a DNA sequence stored using the SeqFastadna class from the \pkg{\link[seqinr]{seqinr}} package.}
\item{alg}{the algorithm for computing the p-value.  If set to \dQuote{\samp{simulate}}, the
p-value is obtained via Monte Carlo simulation.  If set to \dQuote{\samp{lower}}, an analytic
lower bound on the p-value is computed.  If set to \dQuote{\samp{upper}}, an analytic upper
bound on the p-value is computed.  \dQuote{\samp{lower}} and \dQuote{\samp{upper}} are based on formulae
in Hart and Martínez (2011).  a Tighter (though unpublished) lower bound on the p-value  may be obtained by specifying \dQuote{\samp{Lower}}. If
\samp{alg} is specified as \dQuote{\samp{exact}} (the default value), the p-value for the test is
computed exactly.}
  \item{n}{The number of replications to use for Monte Carlo simulation.  
  If computationally feasible, a value >= 10000000 is recommended.}
}
\details{
The first argument may be a character vector representing a DNA sequence, a DNA
sequence represented using the SeqFastadna class from the \pkg{seqinr} package,
or a vector containing the relative frequencies of the A, C, G and T nucleic
acids.

Let A, C, G and T denote the relative frequencies of the nucleotide bases
appearing in  a DNA sequence. This function carries out a statistical hypothesis
test that the relative frequencies satisfy the relation \eqn{A+G=C+T}, or that
purines \eqn{\{A, G\}}{{A,G}} occur equally as often as pyrimidines \eqn{\{C,T\}}{{C,T}} in a DNA sequence.
The relationship can be rewritten as \eqn{A-T=C-G}, from which it is easy to see
that the property being tested is a generalisation of Chargaff's second parity
rule for mononucleotides, which states that \eqn{A=T} and \eqn{C=G}. The test is
set up as follows:

\eqn{H_0}{H0}: \eqn{A+G \neq C+T}{A+G != C+T} \cr
\eqn{H_1}{H1}: \eqn{A+G = C+T}

The vector \eqn{(A,C,G,T)} is assumed to come from a Dirichlet(1,1,1,1)
distribution on the 3-simplex under the null hypothesis.

The test statistic \eqn{\eta_V}{etaV} is the Euclidean distance from the
relative frequency vector \eqn{(A,C,G,T)} to the closest point in the square set
\eqn{\theta_V=\{(x,y,1/2-x,1/2-y) : 0 <= x,y <= 1/2\}}{thetaV = {(x,y,1/2-x,1/2-
y) : 0 <= x,y <= 1/2}}, which divides the 3-simplex into two equal parts.
\eqn{\eta_V}{etaV} lies in the range \eqn{[0,\sqrt{3/8}]}{[0,sqrt(3/8)]}.
}
\note{
\code{agct.test(x, alg="upper")} is equivalent to \code{ag.test(x,
alg="simplex")} except that the p-value computed using the formula for
\samp{alg="upper"} is exact for the test statistic \eqn{\eta_V^*}{etaV*} used in
\code{\link{ag.test}}, whereas it is merely an upper bound on the p-value for
\eqn{\eta_V}{etaV}.
}
\value{
A list with class "htest.ext" containing the following components:

\item{statistic}{the value of the test statistic.}
\item{p.value}{the p-value of the test.}
\item{method}{a character string indicating what type of test was performed.}
\item{data.name}{a character string giving the name of the data.}
\item{estimate}{the probability vector used to derive the test statistic.}
\item{stat.desc}{a brief description of the test statistic.}
\item{null}{the null hypothesis (\eqn{H_0}{H0}) of the test.}
\item{alternative}{the alternative hypothesis (\eqn{H_1}{H1}) of the test.}
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
\code{\link{chargaff0.test}}, \code{\link{chargaff1.test}}, 
\code{\link{chargaff2.test}}, \code{\link{ag.test}},
\code{\link{chargaff.gibbs.test}}
}
\examples{
#Demonstration on real viral sequence
data(pieris)
agct.test(pieris)

#Simulate synthetic DNA sequence that does not exhibit Purine-Pyrimidine parity
trans.mat <- matrix(c(.4, .1, .4, .1, .2, .1, .6, .1, .4, .1, .3, .2, .1, .2, .4, .3), 
ncol=4, byrow=TRUE)
seq <- simulateMarkovChain(500000, trans.mat, states=c("a", "c", "g", "t"))
agct.test(seq)
}
\keyword{htest}

