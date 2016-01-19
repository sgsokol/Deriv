\name{Deriv-package}
\alias{Deriv-package}
\docType{package}
\title{
Symbolic Differentiation
}
\description{
R already contains two differentiation functions: D and deriv.

R's existing functions have several limitations:

\itemize{
 \item the derivatives table can't be modified at runtime, and is only available
in C.
 \item function cannot substitute function calls.  eg:
      
      \code{f <- function(x, y) x + y; deriv(~f(x, x^2), "x")}
}

The advantages of this package include:

\itemize{
 \item It is entirely written in R, so would be easier to maintain.
 \item Can differentiate function calls:
   \itemize{
      \item if the function is in the derivative table, then the chain rule
      is applied.
      \item if the function is not in the derivative table (or it is anonymous),
      then the function body is substituted in.
      \item these two methods can be mixed.  An entry in the derivative table
      need not be self-contained -- you don't need to provide an infinite
      chain of derivatives.
   }
 \item It's easy to add custom entries to the derivatives table, e.g.
      
      \code{drule[["cos"]] <- alist(x=-sin(x))}
 \item The output can be an executable function, which makes it suitable for use in
optimization problems.
}
}

\details{
\tabular{ll}{
Package: \tab Deriv\cr
Type: \tab Package\cr
Version: \tab 3.6.1\cr
Date: \tab 2016-01-19\cr
License: \tab GPL (>= 3)\cr
}
Two main functions are Deriv() for differentiating
and Simplify() for simplifying symbolically.
}

\author{
Andrew Clausen, Serguei Sokol

Maintainer: Serguei Sokol (sokol at insa-toulouse.fr)
}
\references{
\url{https://andrewclausen.net/computing/deriv.html}
}
\keyword{ package }
\seealso{
\code{\link{D}}, \code{\link{deriv}}, packages Ryacas, rSymPy
}
\examples{
\dontrun{f <- function(x) x^2}
\dontrun{Deriv(f)}
# function (x)
# 2 * x
}

