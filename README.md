Deriv
=====

Symbollic differentiation
-------------------------
This software is written in R by Andrew Clausen <clausen at econ.upenn.edu> in 2007.

Thanks to Mark Reid <mark.reid at anu.edu.au> for a patch, applied 21/2/2009.

In 2014, Andrew has passed the maintenance to Serguei Sokol <sokol at insa-toulouse.fr>.

Installation
------------
devtools::install_github("sgsokol/Deriv")

Usage
-----
In R session do:
> f <- function(x) x^2

> Deriv.function(f)

For more information and examples:
> ?Deriv
