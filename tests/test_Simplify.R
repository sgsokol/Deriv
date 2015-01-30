library(testthat)
library(Deriv)

context("Symbolic simplifications")
f=function(x) {}

expect_equal_lang=function(t, r) {
   test=substitute(t)
   right=substitute(r)
   ans=Simplify(as.expression(test))[[1]]
   #print(deparse(ans))
   eval(bquote(expect_equal(format1(quote(.(ans))), format1(quote(.(right))))))
}
test_that("rational simplifications", {
   expect_equal_lang(a*b, a*b) # no change must occur
   expect_equal_lang(a/a, 1) # complete single simplification
   expect_equal_lang(2/2, 1) # complete single numeric simplification
   expect_equal_lang(a*b/(b*a), 1) # complete multiple simplification
   expect_equal_lang(a/(b*a), 1/b) # non complete single simplification in denominator
   expect_equal_lang(-a/(b*a), -1/b) # non complete single negative simplification in denominator
   expect_equal_lang(2/(b*2), 1/b) # non complete single numeric simplification in denominator
   expect_equal_lang(-2/(b*2), -1/b) # non complete single negative numeric simplification in denominator
   expect_equal_lang((a*b)/b, a) # non complete single simplification in numerator
   expect_equal_lang((a*b)/-b, -a) # non complete single simplification in numerator
   expect_equal_lang((a*2)/2, a) # non complete single numeric simplification in numerator
   expect_equal_lang((a*2)/-2, -a) # non complete single negative numeric simplification in numerator
   expect_equal_lang(a*c/(c*b*a), 1/b) # non complete multiple simplification (denominator)
   expect_equal_lang(a*-c/(c*b*a), -1/b) # non complete multiple negative simplification (denominator)
   expect_equal_lang((a*c*b)/(c*a), b) # non complete multiple simplification (numerator)
   expect_equal_lang((-a*c*b)/(c*a), -b) # non complete multiple negative simplification (numerator)
})
