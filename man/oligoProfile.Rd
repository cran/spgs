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

\name{oligoProfile}
\alias{oligoProfile}
\alias{plot.OligoProfile}
\alias{print.OligoProfile}
\alias{oligoCorr}
\title{
Oligo Profiles and Oligo Profile Correlation Plots of Nucleotide Sequences
}
\description{
Construct a k-mer oligo profile of a nucleotide sequence and print such a 
profile or its reverse complement.  There is also a plot function for producing 
plots of the profile or its reverse complement and for comparing primary and 
complementary strand profiles.
}
\usage{
oligoProfile(x, k, content=c("dna", "rna"), 
case=c("lower", "upper", "as is"), circular=TRUE, disambiguate=TRUE, 
plot=TRUE, ...)
\S3method{plot}{OligoProfile}(x, which=1L, units=c("percentage", "count", "proportion"),
 main=NULL, xlab=NULL, ylab=NULL, ...)
\S3method{print}{OligoProfile}(x, which=1L, units=c("percentage", "count", "proportion"), 
digits=switch(units, percentage=3L, count=NULL, proportion=3L), ...)
%oligoCorr(x, content=c("dna", "rna"), 
%case=c("lower", "upper", "as is"), circular=TRUE, disambiguate=TRUE, ...)
}
\arguments{
  \item{x}{
  a character vector or an object that can be coersed to a character vector.
  }
  \item{k}{
  the k-mer profile to produce.
  }
  \item{content}{
    The content type (\dQuote{\code{dna}} or \dQuote{\code{rna}}) of the input 
  sequence.  \code{oligoProfile} can often detect this automatically based on 
  the presence/absence of \code{t}'s or \code{u}'s, but if neither is present, 
  the \code{content} argument is consulted.  The default value is 
  \dQuote{\code{dna}}.
  }
\item{case}{
  determines how labels for the array should be generated: in
  lowercase, in uppercase or left as is, in which case labels such as \dQuote{b}
  and \dQuote{B} will be seen as distinct symbols and counted separately.
}
\item{circular}{
Determines if the vector should be treated as circular or not.  The default is 
\code{TRUE}, meaning that the start and end of the sequence will be joined 
together for the purpose of counting.
}
\item{disambiguate}{
if set to the default of \code{true}, makes the input sequence unambiguous 
before generating the profile.  Otherwise, ambiguous symbols are treated like 
any other symbols and k-mer counts including them will be computed.
}
\item{plot}{
should a plot of the profile be produced?  The default is \code{TRUE}.
}
\item{which}{
For \code{print}, specifies whether to display the profile for the sequence used 
to generate the OligoProfile object (\code{1}) or the profile of its reverse 
complement (\code{2}).

For the \code{plot} method, \code{which} determines what should be plotted. 
Values \code{1} and \code{2} cause the profile for the original sequence 
(primary strand) or its reverse complement (complementary strand) to be plotted, 
respectively.  Specifying \code{which=3} will plot a comparison of the two 
profiles which can be used to assess compliance with Chargaff's second parity 
rule.

the \code{which} argument may also be specified when calling 
\code{oligoProfile}, in which case it will be passed on to the \code{plot} 
method if the \code{plot} argument is set to \code{TRUE}. 
}
\item{units}{
The oligo profiles can be scaled according to three different units for 
presentation on plots: \dQuote{\code{percentage}}, \dQuote{\code{count}} or 
\dQuote{\code{proportion}}. The default is \dQuote{\code{percentage}}.
}
\item{main}{
The title of the plot.  See \code{\link{plot.default}}.  If not specified, an 
appropriate title is automatically generated.
}
\item{xlab}{
a label for the x-axis of the plot.  See \code{\link{plot.default}}.  If not specified, an 
appropriate label is automatically generated.
}
\item{ylab}{
a label for the y-axis of the plot.  See \code{\link{plot.default}}.  If not specified, an 
appropriate label is automatically generated.
}
\item{digits}{
The number of significant digits to print.  The default is \code{0L} when \code{units} is set to \dQuote{\code{count}} and \code{3L} otherwise.
}
\item{\dots}{
arguments to be passed from or to other functions
%for example, arguments such 
%as \code{which}, \code{units} or other graphical parameters (see 
%\code{\link{par}}) may be passed to the plot function.
}
}
\details{

This function returns the oligo profile for a sequence in an \code{OligoProfile} 
object, which is printed on screen if the \code{plot} parameter is \code{FALSE}. 
An oligo profile is simply the counts of all \code{k}-mers in a sequence for 
some specified value of \code{k}.  

By default, \code{oligoProfile} produces a plot of the oligo profile expressed 
in terms of percentages.  The \code{plot} argument determines if the plot 
should be generated or not and plotting parameters such as \code{main}, 
\code{sub}, etc., may be passed as arguments to the function when \code{plot} is 
\code{TRUE}.  

The \code{plot} method, either called directly or indirectly via the 
\code{oligoProfile} function, can produce either the oligo profile of \code{x} 
(\code{which = 1}), the oligo profile of its reverse complement (\code{which = 
2}), or an interstrand k-mer correlation plot comparing the k-oligo profile 
of\code{x} with that of its reverse complement (\code{which = 3)}.  Such 

Correlation plots effectively show the relationship between k-mers on the 
primary and complementary strands in a DNA duplex and can be used to assess 
compliance with CSPR.  More precisely, one would conclude that a genomic 
sequence complies with CSPR if all the plotted points lie on a diagonal line 
running from the bottom-left corner to the top-right corner of the graph.
}
\value{
A list with class \dQuote{OligoProfile} containing the following components:

\item{name}{a name to identify the source of the profile.}
\item{wordLength}{the value of k used to derive the k-mer profile.}
\item{content}{indicates if the profile pertains to a DNA or RNA sequence.}
\item{case}{indicates how the case of letters was processed before 
producing the profile.}
\item{circular}{indicates whether or not the sequence was considered circular for the purpose of producing the profile.}
\item{disambiguate}{indicates if the sequence was made unambiguous before 
producing the profile.}
\item{profile}{a vector containing the raw counts (frequencies) of all k-mers.}
}
\references{
Albrecht-Buehler, G. (2006)  Asymptotically increasing compliance of genomes 
with Chargaff's second parity rules through inversions and inverted 
transpositions. \emph{PNAS} \bold{103(47)}, 17828--17833.
}
\author{
Andrew Hart and Servet Martínez
}
\seealso{
\code{\link{pair.counts}}, \code{\link{triple.counts}}, 
\code{\link{quadruple.counts}}, \code{\link{cylinder.counts}},
\code{\link{array2vector}}, \code{\link{table2vector}}, \code{\link{disambiguate}}
}
\examples{
data(nanoarchaeum)
#Get the 3-oligo profile of Nanoarchaeum without plotting it
nano.prof <- oligoProfile(nanoarchaeum, 3, plot=FALSE)
nano.prof #print oligo profile as percentages
print(nano.prof, units="count") #print oligo profile as counts
plot(nano.prof) #oligo profile plotted as percentages
plot(nano.prof, units="count") #plot it as counts

#plot the 2-oligo profile of Nanoarchaeum as proportions
oligoProfile(nanoarchaeum, k=3, units="proportion")
}
\keyword{array}
