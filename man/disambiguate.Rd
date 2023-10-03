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

\name{disambiguate}
\alias{disambiguate}
\alias{disambiguate.default}
\alias{disambiguate.SeqFastadna}
\alias{disambiguate.list}
\title{
Disambiguate a Nucleic Sequence
}
\description{

Make a DNA/RNA sequence unambiguous by stripping out all symbols that do 
\bold{not} uniquely specify nucleic acids.  In other words, remove all symbols 
other than \code{a}'s, \code{c}'s, \code{g}'s, \code{t}'s or \code{u}'s from the 
sequence.
}
\usage{
\method{disambiguate}{default}(x, case=c("lower", "upper", "as is"), \dots)
\method{disambiguate}{SeqFastadna}(x, \dots)
\method{disambiguate}{list}(x, \dots)
}
\arguments{
  \item{x}{
A character vector, an object that can be coersed to a character vector or a 
  list of objects that canbe be converted to character vectors.  this argument 
  can also be a \code{SeqFastadna} object provided by the \pkg{seqinr} package.
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
is a single nucleobase, then it represents a single sequence. It will be made 
unambiguous and returned in the same form.

On the other hand, if \code{x} is a vector of character strings, each of which 
represents a nucleic sequence, then the result will bea a character vector in 
which each element contains the unambiguous sequence corresponding to the 
element in \code{x} as a character string.
}
\value{
According to the input \code{x}, a character vector, SeqFastadna object or list 
containing the  completely unambiguous sequence(s) in \code{x}.
}
\author{
Andrew Hart and Servet Martínez
}
%\seealso{
%}
\keyword{datagen}
%\keyword{ts}
