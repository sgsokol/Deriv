Sys.setenv("R_TESTS" = "")

library(testthat)
library(Deriv)

test_check("Deriv")
