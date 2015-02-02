Deriv
=====

Symbollic differentiation
-------------------------
This software is written in R by Andrew Clausen (clausen at econ.upenn.edu) in 2007.

Thanks to Mark Reid (mark.reid at anu.edu.au) for a patch, applied 21/2/2009.

In 2014, Andrew has passed the maintenance to Serguei Sokol (sokol at insa-toulouse.fr).

Installation
------------

    > devtools::install_github("sgsokol/Deriv")

Usage
-----
In R session do:

    > library(Deriv)
    > f <- function(x, n=2) x^n+sin(n*x)     # user defined function to diffierentiate
    > (df <- Deriv(f))                       # -> c(x = n * x^(n - 1) + n * cos(n * x), n = log(x) * x^n + x * cos(n * x))
    > df(2, 3)                               # ->         x         n
                                             # -> 14.880511  7.465518
    
    > Deriv(expression(f(y, 3)), "y")        # -> expression(3 * y^2 + 3 * cos(3 * y))
    > Deriv(~ f(y, 3), "y")                  # -> 3 * y^2 + 3 * cos(3 * y)
    > y <- 2; eval(Deriv(~ f(y, 3), "y"))    # -> 14.88051

For more information and examples:

    > ?Deriv
