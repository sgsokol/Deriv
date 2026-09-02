library(testthat)
library(Deriv)

test_that("Deriv_cpp handles basic derivatives and numeric folding", {
  # Constant derivatives
  expect_equal(Deriv_cpp(quote(42), "x"), 0)
  expect_equal(Deriv_cpp(quote(x), "y"), 0)
  expect_equal(Deriv_cpp(quote(x), "x"), 1)
  
  # Linear & polynomial
  expect_equal(Deriv_cpp(quote(2 * x + 5), "x"), 2)
  expect_equal(Deriv_cpp(quote(x^2), "x"), quote(2 * x))
  expect_equal(Deriv_cpp(quote(x^3 + 3 * x), "x"), quote(3 * (x * x) + 3))
})

test_that("Deriv_cpp simplifies unary minus and sign hoisting", {
  # Double negation: -(-x) => x
  expect_equal(Simplify_cpp(quote(-(-x))), quote(x))
  
  # Multiplication sign hoisting: (-a) * b or a * (-b) => -(a * b)
  expect_equal(Deriv_cpp(quote((-x) * y), "x"), quote(-y))
  expect_equal(Deriv_cpp(quote(x * (-y)), "x"), quote(-y))
  expect_equal(Deriv_cpp(quote((-x) * (-y)), "x"), quote(y))
  
  # Addition with negation: x + (-y) => x - y, (-x) + y => y - x
  expect_equal(Simplify_cpp(quote(x + (-y))), quote(x - y))
  expect_equal(Simplify_cpp(quote((-x) + y)), quote(y - x))
  
  # Division sign hoisting: (-a) / b or a / (-b) => -(a / b)
  expect_equal(Simplify_cpp(quote((-x) / y)), quote(-(x / y)))
  expect_equal(Simplify_cpp(quote(x / (-y))), quote(-(x / y)))
})

test_that("Deriv_cpp handles indexed variables and vector targeting", {
  # Named target vector mapping c(x = "1") => targets x[1]
  expect_equal(Deriv_cpp(quote(x[1] * 2 + y), c(x = "1")), 2)
  expect_equal(Deriv_cpp(quote(x[1] * y), c(x = "1")), quote(y))
  expect_equal(Deriv_cpp(quote(x[1]^2 + x[2]), c(x = "1")), quote(2 * x[1]))
  
  # Multiple indexed targets in a single call
  res <- Deriv_cpp(quote(x[1]^2 + x[2]^3), c(x = "1", x = "2"))
  expect_equal(res, quote(c(x_1 = 2 * x[1], x_2 = 3 * (x[2] * x[2]))))
})

test_that("Deriv_cpp supports multiple standard differentiation targets", {
  res <- Deriv_cpp(quote(x^2 + y^3 + z), c("x", "y", "z"))
  expect_equal(res, quote(c(x = 2 * x, y = 3 * (y * y), z = 1)))
})

test_that("Deriv_cpp handles user function expansion in dedicated env", {
  # Create custom environment parented to baseenv()
  test_env <- new.env(parent = baseenv())
  test_env$my_fun <- function(a, b) a^2 + b
  expect_equal(Deriv_cpp(quote(my_fun(x, y)), "x", env = test_env), quote(2 * x))
  expect_equal(Deriv_cpp(quote(my_fun(x, y)), "y", env = test_env), 1)
})
test_that("Deriv_cpp throws an error when function is in neither drule nor env", {
  test_env <- new.env(parent = baseenv())

  expect_error(
    Deriv_cpp(quote(non_existent_fn(x, y)), "x", env = test_env),
    regexp = "non_existent_fn" # Matches function name in error message
  )
})
test_that("Deriv_cpp substitute default arguments if absent", {
  expect_equal(Deriv_cpp(quote(log(x)), "x"), quote(1 / x)) # Default base = exp(1)
  expect_equal(Deriv_cpp(quote(log(x, base=10)), "x"), quote(1 / (2.30258509299405 * x)))
})
