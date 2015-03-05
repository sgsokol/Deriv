context("Symbolic simplifications")

expect_equal_lang=function(t, r) {
   test=substitute(t)
   right=substitute(r)
   ans=Simplify(as.expression(test))[[1]]
   #print(deparse(ans))
   eval(bquote(expect_equal(format1(quote(.(ans))), format1(quote(.(right))))))
}
test_that("rational simplifications", {
   expect_equal_lang(a*b, a*b) # no change must occur
   expect_equal_lang(a/2, a/2) # no change must occur
   expect_equal_lang(a*2, 2*a) # numeric first
   expect_equal_lang(a/a, 1) # complete single simplification
   expect_equal_lang(2/2, 1) # complete single numeric simplification
   expect_equal_lang(a*b/(b*a), 1) # complete multiple simplification
   expect_equal_lang(a/(b*a), 1/b) # single simplification in denominator
   expect_equal_lang(-a/(b*a), -(1/b)) # single negative simplification in denominator
   expect_equal_lang(2/(b*2), 1/b) # single numeric simplification in denominator
   expect_equal_lang(-2/(b*2), -(1/b)) # single negative numeric simplification in denominator
   expect_equal_lang((a*b)/b, a) # single simplification in numerator
   expect_equal_lang((a*b)/-b, -a) # single simplification in numerator
   expect_equal_lang((a*2)/2, a) # single numeric simplification in numerator
   expect_equal_lang((a*2)/-2, -a) # single negative numeric simplification in numerator
   expect_equal_lang(a*c/(c*b*a), 1/b) # multiple simplification (denominator)
   expect_equal_lang(a*-c/(c*b*a), -(1/b)) # multiple negative simplification (denominator)
   expect_equal_lang((a*c*b)/(c*a), b) # multiple simplification (numerator)
   expect_equal_lang((-a*c*b)/(c*a), -b) # multiple negative simplification (numerator)
})
test_that("log simplifications", {
   expect_equal_lang(log(a), log(a)) # no change must occur
   expect_equal_lang(log(a*b), log(a)+log(b))
   expect_equal_lang(log(exp(a)), a)
   expect_equal_lang(log(a^n), n*log(a))
   expect_equal_lang(log(sqrt(a)), 0.5*log(a))
})
test_that("sqrt simplifications", {
   expect_equal_lang(sqrt(a), sqrt(a)) # no change must occur
   expect_equal_lang(sqrt(a*a), abs(a))
   expect_equal_lang(sqrt(a^4), a^2)
   expect_equal_lang(sqrt(exp(a)), exp(a/2))
   expect_equal_lang(sqrt(a^n), abs(a)^(n/2))
   expect_equal_lang(sqrt(sqrt(a)),a^0.25)
})
test_that("abs simplifications", {
   expect_equal_lang(abs(a), abs(a)) # no change must occur
   expect_equal_lang(abs(a*a), a^2)
})
test_that("factorizations", {
   expect_equal_lang(a+b, a+b) # no change must occur
   expect_equal_lang(a*a+b*a, a*(a+b))
   expect_equal_lang(a^2+b*a^3, a^2*(1+a*b))
   expect_equal_lang(a^2/c**5+b*a^3/d/c**3, a^2*(1/c^2+a*b/d)/c^3)
})
