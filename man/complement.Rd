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

\name{complement}
\alias{complement}
\alias{complement.default}
\alias{complement.SeqFastadna}
\alias{complement.list}
\title{
Complement of a DNA/RNA Sequence
}
\description{
Compute the complement of a DNA or RNA sequence.
}
\usage{
\method{complement}{default}(x, content=c("dna", "rna"), 
case=c("lower", "upper", "as is"), \dots)
\method{complement}{SeqFastadna}(x, \dots)
\method{complement}{list}(x, \dots)
}
\arguments{
  \item{x}{
A character vector, an object that can be coersed to a character vector or a 
  list of objects that canbe be converted to character vectors.  this argument 
  can also be a \code{SeqFastadna} object provided by the \pkg{seqinr} package.
}
  \item{content}{
  The content type of sequence(s).  At present, supported types include 
  \dQuote{\code{dna}} and \dQuote{\code{rna}}.  the default type is 
  \dQuote{\code{dna}}. 
    }
\item{case}{
Determines how symbols in \code{x} should be treated before translating them 
into their complements. \dQuote{\code{lower}}, the default behaviour, converts 
all symbols to lowercase while \dQuote{\code{upper}} converts them to 
uppercase.  \dQuote{\code{as is}} allows the symbols to pass unchanged so that the case 
of each output symbol matches that of the corresponding input symbol.
}
\item{\dots}{
Arguments to be passed from or to other functions.
}
}
\details{
If \code{x} is a SeqFastadna object or a character vector in which each element 
is a single nucleobase, then it represents a single sequence and its 
complementary sequence will be returned in the same form.

On the other hand, if \code{x} is a vector of character strings, each of which 
represents a nucleic sequence, then the result will bea a character vector in 
which each element contains the complement of the corresponding element 
in \code{x} as a character string.
}
\value{
According to the input \code{x}, a character vector, SeqFastadna object or list 
containing the complement(s) of the sequence(s) in \code{x}.
}
\author{
Andrew Hart and Servet Martínez
}
\seealso{
\code{\link{reverseComplement}}
}
\examples{
complement("actg")
complement(c("t", "g", "a"))

#List of sequences
some.dna <- list("atgcgtcgttaa", c("g", "t", "g", "a", "a", "a"))
complement(some.dna)

#RNA sequence example
complement(c("a", "u", "g"), content="rna")

#Examples of lowercase, uppercase and as-is conversion
mixed.case <- c("t", "G", "g", "C", "a")
complement(mixed.case)
complement(mixed.case, case="upper")
complement(mixed.case, case="as is")
}
\keyword{datagen}
