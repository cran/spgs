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

\name{Pieris Rapae Granulovirus Genome}
\docType{data}
\alias{pieris}
\title{DNA sequence for the Pieris Rapae Granulovirus Genome}
\description{
This data set contains the DNA sequence for the \emph{Pieris rapae} granulovirus genome.  The Accession number for this sequence is 
NC_013797.1.
  }
\usage{pieris}
\format{a \code{\link[seqinr]{SeqFastadna}} object.}
\source{The NCBI ftp server at ftp://ftp.ncbi.nlm.nih.gov in the 
/genomes/viruses directory.}
\seealso{
\code{\link{nanoarchaeum}}
}
\examples{
data(pieris)
pair.counts(pieris)
}
\keyword{datasets}

