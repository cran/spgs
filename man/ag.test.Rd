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

\name{ag.test}
\alias{ag.test}
\title{Test of Purine-Pyrimidine Parity Based on Purine Count}
\description{
Performs a test proposed by Hart and Martínez (2011) for the equivalence of the
relative frequencies of purines (\eqn{A+G}) and pyrimidines (\eqn{C+T}) in DNA
sequences.  It does this by checking whether or not the mononucleotide
frequencies of a DNA sequence satisfy the relationship A+G=C+T.
}
 \usage{
ag.test(x, type=c("interval", "simplex"))
}
\arguments{
\item{x}{either a vector containing the relative frequencies of each of the 4 nucleotides 
A, C, G, T, a character vector representing a DNA sequence in which each element contains a single nucleotide, or
a DNA sequence stored using the SeqFastadna class from the \pkg{\link[seqinr]{seqinr}} package.}
\item{type}{Specifies one of two possible tests to perform, both of which are
based on the same test statistic, but assuming different forms of the Dirichlet
distribution under the null. \dQuote{\samp{simplex}} assumes a Dirichlet(1,1,1,1) distribution
on the 3-simplex while \dQuote{\samp{interval}} assumes a Dirichlet(1,1) (uniform)
distribution on the unit interval.  The default is \dQuote{\samp{interval}}.}
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

If \samp{type} is set to \dQuote{\samp{simplex}}, the vector \eqn{(A,C,G,T)} is
assumed to come from a Dirichlet(1,1,1,1) distribution on the 3-simplex under
the null hypothesis.  Otherwise, if \samp{type} is set to
\dQuote{\samp{interval}}, it is assumed under the null hypothesis that
\eqn{(A+G,C+T)} ~ Dirichlet(1,1) or, in other words, \eqn{A+G} and \eqn{C+T} are
uniformly distributed on the unit interval and satisfy \eqn{A+G+C+T=1}.

In both cases, the test statistic is \eqn{\eta_V^* = |A+G-0.5|}{etaV* = |A+G-0.5|}.
}
\value{
A list with class "htest.ext" containing the following components:

\item{statistic}{the value of the test statistic.}
\item{p.value}{the p-value of the test. Only included if \samp{no.p.value} is \samp{FALSE}.}
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
\code{\link{chargaff2.test}}, \code{\link{agct.test}},
\code{\link{chargaff.gibbs.test}}
}
\examples{
#Demonstration on real viral sequence
data(pieris)
ag.test(pieris)
ag.test(pieris, type="simplex")

#Simulate synthetic DNA sequence that does not exhibit Purine-Pyrimidine parity
trans.mat <- matrix(c(.4, .1, .4, .1, .2, .1, .6, .1, .4, .1, .3, .2, .1, .2, .4, .3), 
ncol=4, byrow=TRUE)
seq <- simulateMarkovChain(500000, trans.mat, states=c("a", "c", "g", "t"))
ag.test(seq)
}
\keyword{htest}

