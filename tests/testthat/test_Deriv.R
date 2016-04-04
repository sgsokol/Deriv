context(paste("Symbolic differentiation rules v", packageVersion("Deriv"), sep=""))
lc_orig=Sys.getlocale(category = "LC_COLLATE")
Sys.setlocale(category = "LC_COLLATE", locale = "C")

num_test_deriv <- function(fun, larg, narg=1, h=1.e-5, tolerance=2000*h^2) {
   # test the first derivative of a function fun() (given as a character
   # string) by Deriv() and central difference.
   # larg is a named list of parameters to pass to fun
   # narg indicates by which of fun's arguments the differentiation must be made
   # h is the small perturbation in the central differentiation: x-h and x+h 
   # Parameter tolerance is used in comparison test.
   if (length(names(larg)) == 0)
      stop(sprintf("No argument for function %s() to differentiate. There must be at leat one argument.", fun))
   if (h <= 0)
      stop("Parameter h must be positive")
   larg_ph=larg_mh=larg
   larg_ph[[narg]]=larg_ph[[narg]]+h
   larg_mh[[narg]]=larg_mh[[narg]]-h
   f_ph=do.call(fun, larg_ph)
   f_mh=do.call(fun, larg_mh)
   dnum=(f_ph-f_mh)/(2*h)
   sym_larg=larg
   nm_x=names(larg)[narg]
   sym_larg[[narg]]=as.symbol(nm_x)
   flang=as.symbol(fun)
   dsym=do.call(as.function(c(sym_larg, Deriv(as.call(c(flang, sym_larg)), nm_x))), larg, quote=TRUE)
#cat(sprintf("comparing %s by %s\n", format1(as.call(c(flang, larg))), nm_x))
   expect_equal(dnum, dsym, tolerance=tolerance, info=sprintf("%s by %s", format1(as.call(c(flang, larg))), nm_x))
}

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
   body(f)=ref
#cat("\nf deriv=", format1(ans), "\n", sep="")
#cat("\nsimplify=", format1(Simplify(ans)), "\n", sep="")
#cat("f ref=", format1(f), "\n", sep="")
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
   expect_equal_deriv(2**x, 0.693147180559945 * 2^x)
   expect_equal_deriv(sin(x), cos(x))
   expect_equal_deriv(cos(x), -sin(x))
   expect_equal_deriv(tan(x), 1/cos(x)^2)
   expect_equal_deriv(asin(x), 1/sqrt(1 - x^2))
   expect_equal_deriv(acos(x), -(1/sqrt(1 - x^2)))
   expect_equal_deriv(atan(x), 1/(1+x^2))
   expect_equal_deriv(atan2(x, y), y/(x^2+y^2))
   expect_equal_deriv(atan2(0.5, x), -(0.5/(0.25 + x^2)))
   expect_equal_deriv(exp(x), exp(x))
   expect_equal_deriv(expm1(x), exp(x))
   expect_equal_deriv(log(x), 1/x)
   expect_equal_deriv(log1p(x), 1/(1+x))
   expect_equal_deriv(abs(x), sign(x))
   expect_equal_deriv(sign(x), 0)
   expect_equal_deriv(sinh(x), cosh(x))
   expect_equal_deriv(cosh(x), sinh(x))
   expect_equal_deriv(tanh(x), 1-tanh(x)^2)
})
if (getRversion() >= "3.1.0") {
   test_that("trigonometric functions with pi", {
      expect_equal_deriv(sinpi(x), pi*cospi(x))
      expect_equal_deriv(cospi(x), -(pi*sinpi(x)))
      expect_equal_deriv(tanpi(x), pi/cospi(x)**2)
   })
}
test_that("special functions", {
   expect_equal_deriv(beta(x, y), beta(x, y) * (digamma(x) - digamma(x + y)))
   expect_equal_deriv(beta(x, y), beta(x, y) * (digamma(y) - digamma(x + y)), "y")
   expect_equal_deriv(besselI(x, 0), besselI(x, 1))
   expect_equal_deriv(besselI(x, 0, FALSE), besselI(x, 1))
   expect_equal_deriv(besselI(x, 0, TRUE), besselI(x, 1, TRUE)-besselI(x, 0, TRUE))
   expect_equal_deriv(besselI(x, 1), 0.5 * (besselI(x, 0) + besselI(x, 2)))
   expect_equal_deriv(besselI(x, 1, FALSE), 0.5 * (besselI(x, 0) + besselI(x, 2)))
   expect_equal_deriv(besselI(x, 1, TRUE), 0.5 * (besselI(x, 0, TRUE) + besselI(x, 2, TRUE))-besselI(x, 1, TRUE))
   expect_equal_deriv(besselI(x, n), if (n == 0) besselI(x, 1) else 0.5 * (besselI(x, 1 + n) + besselI(x, n - 1)))
   expect_equal_deriv(besselI(x, n, TRUE), (if (n == 0) besselI(x, 1, TRUE) else 0.5 * (besselI(x, 1 + n, TRUE) + besselI(x, n - 1, TRUE)))-besselI(x, n, TRUE))
   expect_equal_deriv(besselK(x, 0), -besselK(x, 1))
   expect_equal_deriv(besselK(x, 0, FALSE), -besselK(x, 1))
   expect_equal_deriv(besselK(x, 0, TRUE), besselK(x, 0, TRUE)-besselK(x, 1, TRUE))
   expect_equal_deriv(besselK(x, 1), -(0.5 * (besselK(x, 0) + besselK(x, 2))))
   expect_equal_deriv(besselK(x, 1, FALSE), -(0.5 * (besselK(x, 0) + besselK(x, 2))))
   expect_equal_deriv(besselK(x, 1, TRUE), besselK(x, 1, TRUE)-0.5 * (besselK(x, 0, TRUE) + besselK(x, 2, TRUE)))
   expect_equal_deriv(besselK(x, n), if (n == 0) -besselK(x, 1) else -(0.5 * (besselK(x, 1 + n) + besselK(x, n - 1))))
   expect_equal_deriv(besselK(x, n, FALSE), if (n == 0) -besselK(x, 1) else -(0.5 * (besselK(x, 1 + n) + besselK(x, n - 1))))
   expect_equal_deriv(besselK(x, n, TRUE), besselK(x, n, TRUE)+if (n == 0) -besselK(x, 1, TRUE) else -(0.5 * (besselK(x, 1 + n, TRUE) + besselK(x, n - 1, TRUE))))
   expect_equal_deriv(besselJ(x, 0), -besselJ(x, 1))
   expect_equal_deriv(besselJ(x, 1), 0.5 * (besselJ(x, 0) - besselJ(x, 2)))
   expect_equal_deriv(besselJ(x, n), if (n == 0) -besselJ(x, 1) else 0.5 * (besselJ(x, n - 1) - besselJ(x, 1 + n)))
   expect_equal_deriv(besselY(x, 0), -besselY(x, 1))
   expect_equal_deriv(besselY(x, 1), 0.5 * (besselY(x, 0) - besselY(x, 2)))
   expect_equal_deriv(besselY(x, n), if (n == 0) -besselY(x, 1) else 0.5 * (besselY(x, n - 1) - besselY(x, 1 + n)))
   expect_equal_deriv(gamma(x), digamma(x) * gamma(x))
   expect_equal_deriv(lgamma(x), digamma(x))
   expect_equal_deriv(digamma(x), trigamma(x))
   expect_equal_deriv(trigamma(x), psigamma(x, 2L))
   expect_equal_deriv(psigamma(x), psigamma(x, 1L))
   expect_equal_deriv(psigamma(x, n), psigamma(x, 1L+n))
   expect_equal_deriv(beta(x, y), beta(x, y) * (digamma(x) - digamma(x + y)))
   expect_equal_deriv(beta(x, y), beta(x, y) * (digamma(y) - digamma(x + y)), "y")
   expect_equal_deriv(lbeta(x, y), digamma(x) - digamma(x + y))
   expect_equal_deriv(lbeta(x, y), digamma(y) - digamma(x + y), "y")
})
test_that("probability densities", {
   expect_equal_deriv(dbinom(5,3,x), 3 * ((3 - 5 * x) * dbinom(5, 2, x)/(1 - x)^2))
   expect_equal_deriv(dnorm(x, m=0.5), -(dnorm(x, 0.5, 1) * (x - 0.5)))
})
test_that("chain rule: multiply by a const", {
   expect_equal_deriv(a*x, a)
   expect_equal_deriv(a[1]*x, a[1])
   expect_equal_deriv(a[[1]]*x, a[[1]])
   expect_equal_deriv(a$b*x, a$b)
   expect_equal_deriv((a*x)**2, 2*(a^2*x))
   expect_equal_deriv((a*x)**n, a*n*(a*x)^(n-1))
   expect_equal_deriv(sin(a*x), a*cos(a*x))
   expect_equal_deriv(cos(a*x), -(a*sin(a*x)))
   expect_equal_deriv(tan(a*x), a/cos(a*x)^2)
   expect_equal_deriv(exp(a*x), a*exp(a*x))
   expect_equal_deriv(log(a*x), 1/x)
})
test_that("particular cases", {
   expect_equal_deriv(log(x, x), 0)
   expect_equal_deriv(x^n+sin(n*x), n * (cos(n * x) + x^(n - 1)))
   expect_equal_deriv(x*(1-x), 1-2*x)
   expect_equal_deriv(x^x, x^x+x^x*log(x))
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
f <- function(a) (1+a)^(1/a)
f1c <- Deriv(f)
f2c <- Deriv(f1c)
f3c <- Deriv(f2c)
f1 <- Deriv(f, cache.exp=FALSE)
f2 <- Deriv(f1, cache.exp=FALSE)
f3 <- Deriv(f2, cache.exp=FALSE)
a=seq(0.01, 2, len=11)

test_that("expression cache test", {
   expect_equal_deriv(exp(-0.5*(x-m)^2/s^2)/s/sqrt(2*pi), -(exp(-(0.5 * ((x - m)^2/s^2))) * (x - m)/(s^3 * sqrt(2 * pi))))
   expect_equal(g2n(x, m, s), g2c(x, m, s))
   expect_equal(f3(a), f3c(a))
})

# composite function differentiation/caching (issue #6)
f<-function(x){ t<-x^2; log(t) }
g<-function(x) cos(f(x))
test_that("composite function", {
   expect_equal(Deriv(g,"x"), function (x) -(2 * (sin(f(x))/x)))
})

# user function with non diff arguments
ifel<-ifelse
drule[["ifel"]]<-alist(test=NULL, yes=(test)*1, no=(!test)*1)
suppressWarnings(rm(t))
expect_equal(Deriv(~ifel(abs(t)<0.1, t**2, abs(t)), "t"), quote({
    .e2 <- abs(t) < 0.1
    (!.e2) * sign(t) + 2 * (t * .e2)
}))
drule[["ifel"]]<-NULL


# test error reporting
test_that("error reporting", {
   expect_error(Deriv(rnorm), "is not in derivative table", fixed=TRUE)
   expect_error(Deriv(~rnorm(x), "x"), "is not in derivative table", fixed=TRUE)
   expect_error(Deriv(~x+rnorm(x), "x"), "is not in derivative table", fixed=TRUE)
})

# systematic central difference tests
set.seed(7)
test_that("central differences", {
   for (nm_f in ls(drule)) {
      rule <- drule[[nm_f]]
      larg <- rule
      narg <- length(larg)
      larg[] <- runif(narg)
      # possible logical parameters are swithed on/off
      fargs=formals(nm_f)
      ilo=sapply(fargs, is.logical)
      if (any(ilo))
         logrid=do.call(expand.grid, rep(list(c(TRUE, FALSE)), sum(ilo)))
      for (iarg in seq_len(narg)) {
         if (is.null(rule[[iarg]]))
            next
         if (is.null(fargs) || !any(ilo)) {
            suppressWarnings(num_test_deriv(nm_f, larg, narg=iarg))
         } else {
            apply(logrid, 1, function(lv) {
               lolarg=larg
               lolarg[ilo]=lv
               suppressWarnings(num_test_deriv(nm_f, lolarg, narg=iarg))
            })
         }
      }
   }
})

tmp <- Deriv(Deriv(quote(dnorm(x ** 2 - x)), "x"), "x")
test_that("dsym cleaning after nested call", {
   expect_identical(Deriv(quote(.e1*x), "x"), quote(.e1)) # was issue #2
})

# doc examples
fsq <- function(x) x^2
fsc <- function(x, y) sin(x) * cos(y)
f_ <- Deriv(fsc)
fc <- function(x, h=0.1) if (abs(x) < h) 0.5*h*(x/h)**2 else abs(x)-0.5*h
myfun <- function(x, y=TRUE) NULL # do something usefull
dmyfun <- function(x, y=TRUE) NULL # myfun derivative by x.
drule[["myfun"]] <- alist(x=dmyfun(x, y), y=NULL) # y is just a logical
#cat("Deriv(myfun)=", format1(Deriv(myfun)), "\n")
theta <- list(m=0.1, sd=2.)
x <- names(theta)
names(x)=rep("theta", length(theta))

test_that("doc examples", {
   expect_equal_format1(Deriv(fsq), function (x) 2 * x)
   expect_equal_format1(Deriv(fsc), function (x, y) c(x = cos(x) * cos(y), y = -(sin(x) * sin(y))))
   expect_equal(f_(3, 4), c(x=0.6471023, y=0.1068000), tolerance = 1.e-7)
   expect_equal(Deriv(~ fsc(x, y^2), "y"), quote(-(2 * (y * sin(x) * sin(y^2)))))
   expect_equal(Deriv(quote(fsc(x, y^2)), c("x", "y"), cache.exp=FALSE), quote(c(x = cos(x) * cos(y^2), y = -(2 * (y * sin(x) * sin(y^2))))))
   expect_equal(Deriv(expression(sin(x^2) * y), "x"), expression(2 * (x * y * cos(x^2))))
   expect_equal(Deriv("sin(x^2) * y", "x"), "2 * (x * y * cos(x^2))")
   expect_equal(Deriv(fc, "x", cache=FALSE), function(x, h=0.1) if (abs(x) < h) x/h else sign(x))
   expect_equal(Deriv(myfun(z^2, FALSE), "z"), quote(2 * (z * dmyfun(z^2, FALSE))))
   expect_equal(Deriv(~exp(-(x-theta$m)**2/(2*theta$sd)), x, cache.exp=FALSE),
    quote(c(theta_m = exp(-((x - theta$m)^2/(2 * theta$sd))) * (x - theta$m)/theta$sd, 
    theta_sd = 2 * (exp(-((x - theta$m)^2/(2 * theta$sd))) * 
        (x - theta$m)^2/(2 * theta$sd)^2))))
})
drule[["myfun"]] <- NULL
Sys.setlocale(category = "LC_COLLATE", locale = lc_orig)
