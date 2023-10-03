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

\name{diid.disturbance}
\alias{diid.disturbance}
\title{
Construct feasible Random Noise Generating a Bernoulli Process
}
\description{
Produces a sequence of random noise which would generate an observed sequence of
finite symbols provided that the sequence of symbols results from a Bernoulli 
process.
}
\usage{
diid.disturbance(x, random = TRUE,  estimates = FALSE)
}
\arguments{
  \item{x}{
A sequence of finite symbols represented  as a  character vector.
}
    \item{random}{This can be a logical value or a number in the range 0-1. If \samp{TRUE}, random
noise will be generated.  If \samp{FALSE}, the constant value 0.5 will be used
as the noise source.  If a value in the range 0-1 is specified, that value will
be used as a constant noise source.  the default value is \samp{TRUE}.
}
  \item{estimates}{
A logical value specifying if the distribution estimated for the Bernoulli process should be included in
the return.
}
}
\value{
If \samp{estimates} is \samp{TRUE}, returns a list containing the following components:

\item{disturbance}{the sequence of random noise as a numeric vector.}
\item{stat.dist}{The stationary distribution estimated from x.}

Otherwise, if \samp{estimate} is \samp{FALSE}, returns the sequence of random
noise as a numeric vector.
}
\author{
Andrew Hart and Servet Martínez
}
\seealso{
\code{\link{markov.test}}, \code{\link{diid.test}}, \code{\link{markov.disturbance}}
}
\keyword{htest}
\keyword{datagen}
\keyword{models}
