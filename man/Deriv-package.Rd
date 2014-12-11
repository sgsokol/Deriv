\name{Deriv-package}
\alias{Deriv-package}
\docType{package}
\title{
Symbolic differentiation
}
\description{
R already contains two differentiation functions: D and deriv.  D does
simple univariate differentiation.  "deriv" uses D do to multivariate
differentiation.  The output of "D" is an expression, whereas the output of
"deriv" is an executable function.

R's existing functions have several limitations.  They can probably be fixed,
but since they are written in C, this would probably require a lot of work.
Limitations include:
\itemize{
 \item The derivatives table can't be modified at runtime, and is only available
in C.
 \item The output of "deriv" can not be differentiated again.
 \item Neither function can substitute function calls.  eg:
      f <- function(x, y) x + y; deriv(f(x, x^2), "x")
}
So, here are the advantages and disadvantages of this implementation:

GOOD POINTS:
\itemize{
 \item It is entirely written in R, so would be easier to maintain.
 \item Can do multi-variate differentiation.
 \item Can differentiate function calls:
      - if the function is in the derivative table, then the chain rule
      is applied.  For example, if you declared that the derivative of
      sin is cos, then it would figure out how to call cos correctly.
      - if the function is not in the derivative table (or it is anonymous),
      then the function body is substituted in.
      - these two methods can be mixed.  An entry in the derivative table
      need not be self-contained -- you don't need to provide an infinite
      chain of derivatives.
 \item It's easy to add custom entries to the derivatives table.  It could be
easier though... it would be nice if something like
      add.deriv("cos", "-sin(expr)")
worked, rather than the clunky function definitions.  (This is purely a
cosmetic issue, though... everything works as is.)
 \item The output is an executable function, which makes it suitable for use in
optimization problems.
}

BAD POINTS:
\itemize{
 \item Differentiating vector-valued functions doesn't work properly, since
the multiplication code doesn't know when to use scalar vs matrix
multiplication.  Unfortunately, solving this is a hard problem because
we would need to know if an arbitrary expression is a vector or not.
We would have to add extra metadata to do this.  Bottom line: can compute
gradients but not Jacobians or Hessians.
 \item Gives useless error messages when it gets stuck.  This could be fixed.
}
}
\details{
\tabular{ll}{
Package: \tab Deriv\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-12-10\cr
License: \tab GPL (>= 3)\cr
}
The package contains two useful functions: Deriv() for differentiating
expressions, and Deriv.function() for differentiating functions symbolically.
It uses a series of Simplify.something() functions.
}
\author{
Andrew Clausen

Maintainer: Serguei Sokol <sokol@insa-toulouse.fr>
}
\references{
\url{https://andrewclausen.net/computing/deriv.html}
}
\keyword{ package }
\seealso{
\code{\link{D}}, \code{\link{deriv}}, \code{\link[Ryacas]{yacas}}
}
\examples{
\dontrun{f <- function(x) x^2}
\dontrun{Deriv.function(f)}
# function (x) 
# 2 * x
}
