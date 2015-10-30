context("Symbolic simplifications")
lc_orig=Sys.getlocale(category = "LC_COLLATE")
Sys.setlocale(category = "LC_COLLATE", locale = "C")

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
   expect_equal_lang(a/(b*a^x), a^(1 - x)/b) # diff numeric - symbol
   expect_equal_lang(a^x/(b*a), a^(x - 1)/b) # diff symbol - numeric
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
   expect_equal_lang(log(1+x), log1p(x))
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
   expect_equal_lang(3-5*x, 3-5*x) # no change must occur
   expect_equal_lang(a*a+b*a, a*(a+b))
   expect_equal_lang(a^2+b*a^3, a^2*(1+a*b))
   expect_equal_lang(a^2/c**5+b*a^3/d/c**3, a^2*(1/c^2+a*b/d)/c^3)
   expect_equal_lang(a+1-a, 1)
   expect_equal_lang(1+x-1, x)
})
test_that("term order", {
   expect_equal_lang(a(1+x)+a(x-1), a(1+x)+a(x-1)) # no change must occur
   expect_equal_lang(a(x-1)+a(x+1), a(1+x)+a(x-1)) # "+" is before "-" in C collate. It is inverse in fr_FR.UTF8
   expect_equal_lang(1+a, 1+a) # no change must occur
   expect_equal_lang(a+1, 1+a)
})
test_that("{...; const}", {
   expect_equal_lang({a=x^2; 0}, 0)
})

context("expression caching")

test_that("Cache test", {
   expect_equal(Cache(quote(c(a+b, sin(a+b)))), quote({
    .e1 <- a + b
    c(.e1, sin(.e1))
    }
    ))
   expect_equal(Cache(quote({t=a+b; c(a+b, sin(a+b))})), quote({
    t = a + b
    c(t, sin(t))
    }
    ))
   expect_equal(Cache(Simplify(deCache(quote({t=a+b; c(a+b, sin(a+b))})))), quote({
    .e1 <- a + b
    c(.e1, sin(.e1))
    }
    ))
   expect_equal(Cache(quote({a=x^2; b=x^2})), quote({
    a = x^2
    b = a
    }
    ))
})
Sys.setlocale(category = "LC_COLLATE", locale = lc_orig)
