context("Symbolic differentiation rules")
f=function(x) {} # empty place holder

expect_equal_deriv <- function(t, r, nmvar="x") {
   test=substitute(t)
   ref=substitute(r)
   # compare as language
   ans=Deriv(test, nmvar, cache.exp=FALSE)
   #print(deparse(ans))
   eval(bquote(expect_equal(format1(quote(.(ans))), format1(quote(.(ref))))))
   # compare as string
   ans=Deriv(format1(test), nmvar, cache.exp=FALSE)
   #print(ans)
   eval(bquote(expect_equal(.(ans), format1(quote(.(ref))))))
   # compare as formula
   ans=Deriv(call("~", test), nmvar, cache.exp=FALSE)
   #print(deparse(ans))
   eval(bquote(expect_equal(format1(quote(.(ans))), format1(quote(.(ref))))))
   # compare as expression
   ans=Deriv(as.expression(test), nmvar, cache.exp=FALSE)
   #print(deparse(ans))
   eval(bquote(expect_equal(format1(.(ans)), format1(expression(.(ref))))))
   # compare as function
   body(f)=test
   ans=Deriv(f, nmvar, cache.exp=FALSE)
   #print(deparse(ans))
   body(f)=ref
   eval(bquote(expect_equal(quote(.(ans)), quote(.(f)))))
   # compare with central differences
   x=seq(0.1, 1, len=10)
   h=1.e-7
   suppressWarnings(f1 <- try(sapply(x-h, function(val) eval(test, list(x=val))), silent=TRUE))
   suppressWarnings(f2 <- try(sapply(x+h, function(val) eval(test, list(x=val))), silent=TRUE))
   if (!inherits(f1, "try-error") && !inherits(f2, "try-error")) {
      numder=(f2-f1)/h/2
      refder=sapply(x, function(val) eval(ref, list(x=val)))
      i=is.finite(refder) & is.finite(numder)
      expect_more_than(sum(i), 0, label=sprintf("length of central diff for %s", format1(test)))
      expect_equal(numder[i], refder[i], tolerance=5.e-8, label=sprintf("Central diff. of '%s'", format1(test)), expected.label=sprintf("'%s'", format1(ref)))
   }
}
expect_equal_format1 <- function(t, r) {
   eval(bquote(expect_equal(format1(.(t)), format1(.(r)))))
}
test_that("elementary functions", {
   expect_equal(Deriv("x", "x"), "1")
   expect_equal(Deriv(quote(x), "x"), 1)
   expect_equal(Deriv(quote((x)), "x"), 1)
   expect_equal_deriv(x**2, 2*x)
   expect_equal_deriv(x**n, n*x^(n-1))
   expect_equal_deriv(sin(x), cos(x))
   expect_equal_deriv(cos(x), -sin(x))
   expect_equal_deriv(tan(x), 1/cos(x)^2)
   expect_equal_deriv(asin(x), 1/sqrt(1 - x^2))
   expect_equal_deriv(acos(x), -(1/sqrt(1 - x^2)))
   expect_equal_deriv(atan(x), 1/(1+x^2))
   expect_equal_deriv(atan2(x, y), y/(x^2+y^2))
   expect_equal_deriv(exp(x), exp(x))
   expect_equal_deriv(expm1(x), exp(x))
   expect_equal_deriv(log(x), 1/x)
   expect_equal_deriv(log1p(x), 1/(x+1))
   expect_equal_deriv(abs(x), sign(x))
   expect_equal_deriv(sign(x), 0)
   expect_equal_deriv(sinpi(x), pi*cospi(x))
   expect_equal_deriv(cospi(x), -(pi*sinpi(x)))
   expect_equal_deriv(tanpi(x), pi/cospi(x)**2)
   expect_equal_deriv(sinh(x), cosh(x))
   expect_equal_deriv(cosh(x), sinh(x))
   expect_equal_deriv(tanh(x), 1-tanh(x)^2)
})
test_that("special functions", {
   expect_equal_deriv(beta(x, y), beta(x, y) * (digamma(x) - digamma(x + y)))
   expect_equal_deriv(beta(x, y), beta(x, y) * (digamma(y) - digamma(x + y)), "y")
   expect_equal_deriv(besselI(x, 0), besselI(x, 1))
   expect_equal_deriv(besselI(x, 0, FALSE), besselI(x, 1, FALSE))
   expect_equal_deriv(besselI(x, 0, TRUE), besselI(x, 1, TRUE)-besselI(x, 0, TRUE))
   expect_equal_deriv(besselI(x, 1), 0.5 * (besselI(x, 0) + besselI(x, 2)))
   expect_equal_deriv(besselI(x, 1, FALSE), 0.5 * (besselI(x, 0, FALSE) + besselI(x, 2, FALSE)))
   expect_equal_deriv(besselI(x, 1, TRUE), 0.5 * (besselI(x, 0, TRUE) + besselI(x, 2, TRUE))-besselI(x, 1, TRUE))
   expect_equal_deriv(besselI(x, n), if (n == 0) besselI(x, 1) else 0.5 * (besselI(x, n - 1) + besselI(x, n + 1)))
   expect_equal_deriv(besselI(x, n, TRUE), (if (n == 0) besselI(x, 1, TRUE) else 0.5 * (besselI(x, n - 1, TRUE) + besselI(x, n + 1, TRUE)))-besselI(x, n, TRUE))
   expect_equal_deriv(besselK(x, 0), -besselK(x, 1))
   expect_equal_deriv(besselK(x, 0, FALSE), -besselK(x, 1, FALSE))
   expect_equal_deriv(besselK(x, 0, TRUE), besselK(x, 0, TRUE)-besselK(x, 1, TRUE))
   expect_equal_deriv(besselK(x, 1), -(0.5 * (besselK(x, 0) + besselK(x, 2))))
   expect_equal_deriv(besselK(x, 1, FALSE), -(0.5 * (besselK(x, 0, FALSE) + besselK(x, 2, FALSE))))
   expect_equal_deriv(besselK(x, 1, TRUE), besselK(x, 1, TRUE)-0.5 * (besselK(x, 0, TRUE) + besselK(x, 2, TRUE)))
   expect_equal_deriv(besselK(x, n), if (n == 0) -besselK(x, 1) else -(0.5 * (besselK(x, n - 1) + besselK(x, n + 1))))
   expect_equal_deriv(besselK(x, n, FALSE), if (n == 0) -besselK(x, 1, FALSE) else -(0.5 * (besselK(x, n - 1, FALSE) + besselK(x, n + 1, FALSE))))
   expect_equal_deriv(besselK(x, n, TRUE), (if (n == 0) -besselK(x, 1, TRUE) else -(0.5 * (besselK(x, n - 1, TRUE) + besselK(x, n + 1, TRUE))))+besselK(x, n, TRUE))
   expect_equal_deriv(besselJ(x, 0), -besselJ(x, 1))
   expect_equal_deriv(besselJ(x, 1), 0.5 * (besselJ(x, 0) - besselJ(x, 2)))
   expect_equal_deriv(besselJ(x, n), if (n == 0) -besselJ(x, 1) else 0.5 * (besselJ(x, n - 1) - besselJ(x, n + 1)))
   expect_equal_deriv(besselY(x, 0), -besselY(x, 1))
   expect_equal_deriv(besselY(x, 1), 0.5 * (besselY(x, 0) - besselY(x, 2)))
   expect_equal_deriv(besselY(x, n), if (n == 0) -besselY(x, 1) else 0.5 * (besselY(x, n - 1) - besselY(x, n + 1)))
   expect_equal_deriv(gamma(x), digamma(x) * gamma(x))
   expect_equal_deriv(lgamma(x), digamma(x))
   expect_equal_deriv(digamma(x), trigamma(x))
   expect_equal_deriv(trigamma(x), psigamma(x, 2L))
   expect_equal_deriv(psigamma(x), psigamma(x, 1L))
   expect_equal_deriv(psigamma(x, n), psigamma(x, n+1L))
   expect_equal_deriv(beta(x, y), beta(x, y) * (digamma(x) - digamma(x + y)))
   expect_equal_deriv(beta(x, y), beta(x, y) * (digamma(y) - digamma(x + y)), "y")
   expect_equal_deriv(lbeta(x, y), digamma(x) - digamma(x + y))
   expect_equal_deriv(lbeta(x, y), digamma(y) - digamma(x + y), "y")
})
test_that("chain rule: multiply by a const", {
   expect_equal_deriv(a*x, a)
   expect_equal_deriv((a*x)**2, 2*(a^2*x))
   expect_equal_deriv((a*x)**n, a*n*(a*x)^(n-1))
   expect_equal_deriv(sin(a*x), a*cos(a*x))
   expect_equal_deriv(cos(a*x), -(a*sin(a*x)))
   expect_equal_deriv(tan(a*x), a/cos(a*x)^2)
   expect_equal_deriv(exp(a*x), a*exp(a*x))
   expect_equal_deriv(log(a*x), 1/x)
})
test_that("special cases", {
   expect_equal_deriv(log(x, x), 0)
})

# test AD and caching
# gaussian function
g <- function(x, m=0, s=1) exp(-0.5*(x-m)^2/s^2)/s/sqrt(2*pi)
g1c <- Deriv(g, "x") # cache enabled by default
g1n <- Deriv(g, "x", cache.exp=FALSE) # cache disabled
g2c <- Deriv(g1c, "x") # cache enabled by default
g2n <- Deriv(g1n, "x", cache.exp=FALSE) # cache disabled
m <- 0.5
s <- 3.
x=seq(-2, 2, len=11)
test_that("expression cache test", {
   expect_equal_deriv(exp(-0.5*(x-m)^2/s^2)/s/sqrt(2*pi), -(exp(-(0.5 * ((x - m)^2/s^2))) * (x - m)/(s^3 * sqrt(2 * pi))))
   expect_equal(g2n(x, m, s), g2c(x, m, s))
})


# doc examples
fsq <- function(x) x^2
fsc <- function(x, y) sin(x) * cos(y)
f_ <- Deriv(fsc)
fc <- function(x, h=0.1) if (abs(x) < h) 0.5*h*(x/h)**2 else abs(x)-0.5*h
test_that("doc examples", {
   expect_equal_format1(Deriv(fsq), function (x) 2 * x)
   expect_equal_format1(Deriv(fsc), function (x, y) c(x = cos(x) * cos(y), y = -(sin(x) * sin(y))))
   expect_equal(f_(3, 4), c(x=0.6471023, y=0.1068000), tolerance = 1.e-7)
   expect_equal(Deriv(~ fsc(x, y^2), "y"), quote(-(2 * (y * sin(x) * sin(y^2)))))
   expect_equal(Deriv(quote(fsc(x, y^2)), c("x", "y"), cache.exp=FALSE), quote(c(x = cos(x) * cos(y^2), y = -(2 * (y * sin(x) * sin(y^2))))))
   expect_equal(Deriv(expression(sin(x^2) * y), "x"), expression(2 * (x * y * cos(x^2))))
   expect_equal(Deriv("sin(x^2) * y", "x"), "2 * (x * y * cos(x^2))")
   expect_equal(Deriv(fc, "x", cache=FALSE), function(x, h=0.1) if (abs(x) < h) x/h else sign(x))
})
