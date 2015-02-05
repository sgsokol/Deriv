context("Symbolic derivation rules")
f=function(x) {} # empty place holder

expect_equal_deriv <- function(t, r) {
   test=substitute(t)
   ref=substitute(r)
   # compare as language
   ans=Deriv(test, "x")
   #print(deparse(ans))
   eval(bquote(expect_equal(format1(quote(.(ans))), format1(quote(.(ref))))))
   # compare as string
   ans=Deriv(format1(test), "x")
   #print(ans)
   eval(bquote(expect_equal(.(ans), format1(quote(.(ref))))))
   # compare as formula
   ans=Deriv(call("~", test), "x")
   #print(deparse(ans))
   eval(bquote(expect_equal(format1(quote(.(ans))), format1(quote(.(ref))))))
   # compare as expression
   ans=Deriv(as.expression(test), "x")
   #print(deparse(ans))
   eval(bquote(expect_equal(format1(.(ans)), format1(expression(.(ref))))))
   # compare as function
   body(f)=test
   ans=Deriv(f, "x")
   #print(deparse(ans))
   body(f)=ref
   eval(bquote(expect_equal(quote(.(ans)), quote(.(f)))))
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
   expect_equal_deriv(exp(x), exp(x))
   expect_equal_deriv(log(x), 1/x)
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

fsq <- function(x) x^2
fsc <- function(x, y) sin(x) * cos(y)
f_ <- Deriv(fsc)
test_that("doc examples", {
   expect_equal_format1(Deriv(fsq), function (x) 2 * x)
   expect_equal_format1(Deriv(fsc), function (x, y) c(x = cos(x) * cos(y), y = -(sin(x) * sin(y))))
   expect_equal(f_(3, 4), c(x=0.6471023, y=0.1068000), tolerance = 1.e-7)
   expect_equal(Deriv(~ fsc(x, y^2), "y"), quote(-(2 * (sin(x) * y * sin(y^2)))))
   expect_equal(Deriv(quote(fsc(x, y^2)), c("x", "y")), quote(c(x = cos(x) * cos(y^2), y = -(2 * (sin(x) * y * sin(y^2))))))
   expect_equal(Deriv(expression(sin(x^2) * y), "x"), expression(2 * (x * cos(x^2) * y)))
   expect_equal(Deriv("sin(x^2) * y", "x"), "2 * (x * cos(x^2) * y)")
})
