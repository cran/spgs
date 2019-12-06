\name{NEWS}
\title{News for Package \pkg{spgs}}

\section{Changes in version 1.0-3}{
\itemize{
\item Fixed a bit of bad programming style that breaks chargaff.gibbs.test in R 
4.0. More precisely, I now use \code{inherits(x, "try-error")} instead of 
\code{class(x)=="try-error"}.
}
}

\section{Changes in version 1.0-2}{
  \itemize{
    \item Tweaks were made to spgs-package.Rd in response to CRAN maintainer
    since some commented lines were being detected by CRAN checks.
    The file now takes advantage of some package meta data.
    There is no observable difference to the user.
  }
}

\section{Changes in version 1.0-1}{
  \itemize{
    \item Added a couple more cross-references to the \dQuote{Also see} section of a few help files.
    \item Fixed minor internal bugs, e.g., line endings that should be lf, etc.
    \item Refacctored the interface to the C code to use explicit registration in R and have the R code call the registered symbols directly.
  }
}
