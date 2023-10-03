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

\name{rcspr2mat}
\alias{rcspr2mat}
\title{
Random Stochastic Matrices Complying with Chargaff's Second Parity Rule for Dinucleotides
}
\description{
Randomly generate a 4 X 4 stochastic matrix that satisfies Chargaff's second parity rule for dinucleotides.
}
\usage{
rcspr2mat(labels=c("a", "c", "g", "t"))
}
\arguments{
  \item{labels}{
    a vector of labels for the rows and columns of the matrix.  By
  default, this is set to the set of four nucleotides a, c, g and t.
  }
}
\details{
This function randomly generates Stochastic matrices of the form 

\deqn{
\left(\begin{array}{cccc}
b_1&b_2&b_3&1-(b_1+b_2+b_3) \\
zb_6&b_4 & 1-(zb_6+b_4+zb_3) & zb_3 \\
zb_5 &1-(zb_5+b_4+zb_2) & b_4&zb_2 \\
1-(b_5+b_6+b_1)&b_5&b_6&b_1
\end{array}\right)
}{
          b1             b2             b3 1-(b1+b2+b3)
         zb6             b4 1-(zb6+b4+zb3)          zb3
         zb5 1-(zb5+b4+zb2)             b4          zb2
1-(b5+b6+b1)             b5             b6           b1}

where \eqn{b_1}{b1}, \ldots, \eqn{b_6}{b6} are values in the interval (0,1) and
\eqn{z} is a positive number.

Such matrices characterize sequences of DNA that comply with Chargaff's second
parity rule for dinucleotides. See the reference for further information.
}
\value{
A 4 X 4 stochastic matrix satisfying Chargaff's second parity rule.  The rows and columns are labelled according to labels.
}
\references{
Hart, A.G. and Martínez, S. (2011)
Statistical testing of Chargaff's second parity rule in bacterial genome sequences.
\emph{Stoch. Models} \bold{27(2)}, 1--46.
}
\note{
This function is only intended for obtaining samples of matrices complying with
CSPR. It doe snot sample uniformly from the set of all such matrices and hence
is not appropriate for simulation experiments requiring uniformly drawn samples.
}
\author{
Andrew Hart and Servet Martínez
}
\seealso{
\code{\link{rstochmat}}
}
\keyword{array}
\keyword{datagen}
