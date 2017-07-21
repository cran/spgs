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

\name{chargaff0.test}
\alias{chargaff0.test}
\title{Vector Test of Chargaff's Second Parity Rule (CSPR) for Mononucleotides}
\description{
Performs the vector test of Chargaff's second parity rule (CSPR) for mononucleotides proposed in Hart and Martínez (2001).
}
 \usage{
chargaff0.test(x, alg=c("exact", "simulate", "lower", "upper", "Lower", "Upper"), n, 
no.p.value=FALSE)
}
\arguments{
\item{x}{either a vector containing the relative frequencies of each of the 4 nucleotides 
A, C, G, T, a character vector representing a DNA sequence in which each element contains a single nucleotide, or
a DNA sequence stored using the SeqFastadna class from the \pkg{\link[seqinr]{seqinr}} package.}
\item{alg}{the algorithm for computing the p-value.  If set to \dQuote{\samp{simulate}}, the
p-value is obtained via Monte Carlo simulation.  If set to \dQuote{\samp{lower}}, an analytic
lower bound on the p-value is computed.  If set to \dQuote{\samp{upper}}, an analytic upper
bound on the p-value is computed.  \dQuote{\samp{lower}} and \dQuote{\samp{upper}} are based on formulae
in Hart and Martínez (2011).  a Tighter (though unpublished) lower
/upper bound on the p-value  may be obtained by specifying \dQuote{\samp{Lower}}/\dQuote{\samp{Upper}}. If
\samp{type} is specified as \dQuote{\samp{exact}} (the default value),the p-value for the test is
computed exactly for small values of the test statistic and crudely approximated for large values.  See the note below for further details.}
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

Letting A, C, G and T denote the relative frequencies of their corresponding nucleic acids, this function performs the following hypothesis test:

\eqn{H_0}{H0}: \eqn{A\neq T}{A!=T} or \eqn{C\neq G}{C!=G} \cr
\eqn{H_1}{H1}: \eqn{A=T} and \eqn{C=G}

The vector \eqn{(A,C,G,T)} is assumed to come from a Dirichlet(1,1,1,1)
distribution on the 3-simplex under the null hypothesis.

The test statistic is \eqn{\eta_0 = \sqrt{(A-T)^2/2+(C-G)^2/2}}{eta0 = sqrt((A-T)^2/2+(C-G)^2/2)}.
}
\note{
Currently, regardless of the algorithm (\code{alg}) selected, the p-value or
bound is only computed correctly when the test statistic is smaller than or
equal to \eqn{\sqrt{2}/4}{sqrt(2)/4}.  A value of 1 is returned when 
the test statistic is greater than \eqn{\sqrt{2}/4}{sqrt(2)/4}.  This is not
accurate, but shouldn't matter as it is well within the acceptance region of the
null hypothesis.

The algebraically computed bounds on the p-value obtained when \samp{alg} is set
to either \dQuote{\samp{lower}}or  \dQuote{\samp{upper}} are not as tight as
those corresponding to \dQuote{\samp{Lower}} and \dQuote{\samp{Upper}}, which  should be
generally preferred.  However, \dQuote{\samp{exact}} or \dQuote{\samp{simulate}}
should be employed for real- world analysis.

\samp{no.p.value} suppresses computation of the p-value when it is set to
\samp{TRUE}.  This may be useful wen using this function to help simulate the
test statistic.
}
\value{
A list with class "htest.ext" containing the following components:

\item{statistic}{the value of the test statistic.}
\item{p.value}{the p-value of the test. Only included if no.p.value is FALSE.}
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
\code{\link{chargaff1.test}}, \code{\link{chargaff2.test}}, 
\code{\link{agct.test}}, \code{\link{ag.test}}, 
\code{\link{chargaff.gibbs.test}}
}
\examples{
#Demonstration on real bacterial sequence
data(nanoarchaeum)
chargaff0.test(nanoarchaeum)

#Simulate synthetic DNA sequence that does not satisfy Chargaff's second parity rule
trans.mat <- matrix(c(.4, .1, .4, .1, .2, .1, .6, .1, .4, .1, .3, .2, .1, .2, .4, .3), 
ncol=4, byrow=TRUE)
seq <- simulateMarkovChain(500000, trans.mat, states=c("a", "c", "g", "t"))
chargaff0.test(seq)
}
\keyword{htest}
