context(paste("Symbolic differentiation rules v", packageVersion("Deriv"), sep=""))
lc_orig=Sys.getlocale(category = "LC_COLLATE")
Sys.setlocale(category = "LC_COLLATE", locale = "C")

num_test_deriv <- function(fun, larg, narg, h=1.e-5, tolerance=2000*h^2) {
   # test the first derivative of a function fun() (given as a character
   # string) by Deriv() and central difference.
   # larg is a named list of parameters to pass to fun
   # narg indicates the argument name by which the differentiation must be made
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

   sym_larg[[narg]]=as.symbol(narg)
   flang=as.symbol(fun)
   dsym=try(do.call(as.function(c(sym_larg, list(Deriv(as.call(c(flang, sym_larg)), narg)))), larg, quote=TRUE))
   if (inherits(dsym, "try-error")) {
      stop(sprintf("failed to calculate symbolic derivative of '%s'", format1(as.call(c(flang, sym_larg)))))
   }
#cat(sprintf("comparing %s by %s\n", format1(as.call(c(flang, larg))), nm_x))
   expect_equal(as.vector(dnum), as.vector(dsym), tolerance=tolerance, info=sprintf("%s by %s", format1(as.call(c(flang, larg))), narg))
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
   eval(bquote(expect_equal(quote(.(ans)), quote(.(f)), check.environment=FALSE)))
   # compare with central differences
   x=seq(0.1, 1, len=10)
   h=1.e-7
   suppressWarnings(f1 <- try(sapply(x-h, function(val) eval(test, list(x=val))), silent=TRUE))
   suppressWarnings(f2 <- try(sapply(x+h, function(val) eval(test, list(x=val))), silent=TRUE))
   if (!inherits(f1, "try-error") && !inherits(f2, "try-error")) {
      numder=(f2-f1)/h/2
      refder=sapply(x, function(val) eval(ref, list(x=val)))
      i=is.finite(refder) & is.finite(numder)
      expect_gt(sum(i), 0, label=sprintf("length of central diff for %s", format1(test)))
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
   expect_equal_deriv(dbinom(1,3,x), (1 - 3 * x) * dbinom(1, 3, x)/(x * (1 - x)))
   expect_equal_deriv(dnorm(x, m=0.5), -(dnorm(x, 0.5, 1) * (x - 0.5)))
})
test_that("normal quantile", {
   expect_equal_deriv(qnorm(x, mu, lower.tail=FALSE), -(1/dnorm(qnorm(x, mean = mu, sd = 1, lower.tail = FALSE, log.p = FALSE), mean = mu, sd = 1)))
   expect_equal_deriv(qnorm(x, mu, lower.tail=TRUE), 1/dnorm(qnorm(x, mean = mu, sd = 1, lower.tail = TRUE, log.p = FALSE), mean = mu, sd = 1))
   expect_equal_deriv(qnorm(x, mu, log.p=TRUE), exp(x)/dnorm(qnorm(x, mean = mu, sd = 1, lower.tail = TRUE, log.p = TRUE), mean = mu, sd = 1))
   expect_equal_deriv(qnorm(x, mu, log.p=FALSE), 1/dnorm(qnorm(x, mean = mu, sd = 1, lower.tail = TRUE, log.p = FALSE), mean = mu, sd = 1))
})
a=0.1
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
   expect_equal_deriv(x^x, x^x*(1+log(x)))
})
test_that("indexing", {
   expect_equal_deriv(a[['b']], 0)
})
test_that("matrix calculus", {
   expect_equal_deriv(solve(matrix(c(1, x, x**2, x**3), nrow=2, ncol=2)), -solve(matrix(c(1, x, x^2, x^3), nrow = 2, ncol = 2)) %*% matrix(c(0, 1, 2 * x, 3 * x^2), nrow = 2, ncol = 2, byrow = , dimnames = ) %*% solve(matrix(c(1, x, x^2, x^3), nrow = 2, ncol = 2)))
})

test_that("language constructs", {
   expect_equal_deriv(ifelse(x>0, x^2, x^3), ifelse(test=x>0, yes=2*x, no=3*x^2))
   expect_equal_deriv(with(list(c=2), x^c), with(list(c = 2), c * x^(c - 1)))
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
test_that("reused variables", { # (issue #12)
   expect_equal(Deriv(~{sum=x; sum=sum*(1+x); sum=sum*y}, c("x", "y")), quote(c(x = y * (1 + 2 * x), y = x * (1 + x))))
})

# composite function differentiation/caching (issue #6)
f <- function(x){ t<-x^2; log(t) }
g <- function(x) cos(f(x))
test_that("composite function", {
   expect_equal(Deriv(g,"x"), function (x) -(2 * (sin(f(x))/x)), check.environment=FALSE)
})

# user function with non diff arguments
ifel <- ifelse
drule[["ifel"]]<-alist(test=NULL, yes=(test)*1, no=(!test)*1)
suppressWarnings(rm(t))
expect_equal(Deriv(~ifel(abs(t)<0.1, t**2, abs(t)), "t"), quote({
    .e2 <- abs(t) < 0.1
    (!.e2) * sign(t) + 2 * (t * .e2)
}))
rm("ifel", envir=drule)

# long function name (issu #26)
eedddddddddddddddddddddddlog=function(x) log(x)
expect_error(Deriv(function(x) eedddddddddddddddddddddddlog(x)^(1-sig)*exp(x)*h, "x"), NA)

# test error reporting
test_that("error reporting", {
   expect_error(Deriv(rnorm), "is not in derivative table", fixed=TRUE)
   expect_error(Deriv(~rnorm(x), "x"), "is not in derivative table", fixed=TRUE)
   expect_error(Deriv(~x+rnorm(x), "x"), "is not in derivative table", fixed=TRUE)
})

# systematic central difference tests
set.seed(7)
test_that("central differences", {
#browser()
   for (nm_f in ls(drule)) {
      fargs=head(as.list(args(nm_f)), -1L)
      fargs[["..."]]=NULL
      ilo=sapply(fargs, isTRUE) | sapply(fargs, isFALSE)
      rule <- drule[[nm_f]]
      larg <- fargs
      narg <- length(larg)
      if (nm_f == "rep.int") {
         larg["x"]=pi
         larg["times"]=2
      } else if (nm_f == "rep.int") {
         larg["x"]=pi
         larg["length.out"]=2
      } else {
         larg[] <- runif(narg)
      }
      if (nm_f == "det") {
         larg[["x"]]=as.matrix(larg[["x"]])
      } else if (nm_f == "acosh") {
         larg[["x"]]=1+larg[["x"]]
      } else if (nm_f == "diag" || nm_f == "matrix") {
         larg[["nrow"]]=larg[["ncol"]]=1L
         if (nm_f == "matrix") {
#browser()
            larg[["dimnames"]]=NULL
            ilo=ilo[-which(names(ilo) %in% "dimnames")]
         }
      }
      # possible logical parameters are swithed on/off
      if (any(ilo))
         logrid=do.call(expand.grid, rep(list(c(TRUE, FALSE)), sum(ilo)))
      for (arg in names(rule)) {
         if (is.null(rule[[arg]]) || arg == "_missing")
            next
         if (is.null(fargs) || !any(ilo)) {
            
            tmp=try(num_test_deriv(nm_f, larg, narg=arg), silent=TRUE)
            if (inherits(tmp, "try-error")) {
               stop(sprintf("Failed num. deriv test on '%s(%s)'", nm_f, paste(names(larg), larg, sep="=", collapse=", ")))
            }
         } else {
            apply(logrid, 1, function(lv) {
#browser()
               lolarg=larg
               lolarg[ilo]=lv
               if (nm_f == "qnorm" && isTRUE(lolarg[["log.p"]])) {
                  lolarg[["p"]]=log(lolarg[["p"]])
               }
               suppressWarnings(num_test_deriv(nm_f, lolarg, narg=arg))
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

# GMM example
set.seed(777)
ncomp=2
a=runif(ncomp)
a=a/sum(a) # amplitude or weight of each component
m=rnorm(ncomp) # mean
s=runif(ncomp) # sd
# two column matrix of probabilities: one row per x value, one column per component
pn=function(x, a, m, s, log=FALSE) {
  n=length(a)
  structure(vapply(seq(n), function(i) a[i]*dnorm(x, m[i], s[i], log),
    double(length(x))), dim=c(length(x), n))
}
p=function(x, a, m, s) rowSums(pn(x, a, m, s)) # overall probability
dp=Deriv(p, "x")



test_that("doc examples", {
   expect_equal_format1(Deriv(fsq), function (x) 2 * x)
   expect_equal_format1(Deriv(fsc), function (x, y) c(x = cos(x) * cos(y), y = -(sin(x) * sin(y))))
   expect_equal(f_(3, 4), c(x=0.6471023, y=0.1068000), tolerance = 1.e-7)
   expect_equal(Deriv(~ fsc(x, y^2), "y"), quote(-(2 * (y * sin(x) * sin(y^2)))))
   expect_equal(Deriv(quote(fsc(x, y^2)), c("x", "y"), cache.exp=FALSE), quote(c(x = cos(x) * cos(y^2), y = -(2 * (y * sin(x) * sin(y^2))))))
   expect_equal(Deriv(expression(sin(x^2) * y), "x"), expression(2 * (x * y * cos(x^2))))
   expect_equal(Deriv("sin(x^2) * y", "x"), "2 * (x * y * cos(x^2))")
   expect_equal(Deriv(fc, "x", cache=FALSE), function(x, h=0.1) if (abs(x) < h) x/h else sign(x), check.environment=FALSE)
   expect_equal(Deriv(~myfun(z^2, FALSE), "z"), quote(2 * (z * dmyfun(z^2, FALSE))))
   expect_equal(Deriv(~exp(-(x-theta$m)**2/(2*theta$sd)), x, cache.exp=FALSE),
    quote(c(theta_m = exp(-((x - theta$m)^2/(2 * theta$sd))) * (x - theta$m)/theta$sd, 
    theta_sd = 2 * (exp(-((x - theta$m)^2/(2 * theta$sd))) * 
        (x - theta$m)^2/(2 * theta$sd)^2))))
   expect_equal(Deriv(~with(theta, exp(-(x-m)**2/(2*sd))), x, cache.exp=FALSE),
    quote(c(theta_m = with(theta, exp(-((x - m)^2/(2 * sd))) * (x - m)/sd), 
    theta_sd = with(theta, 2 * (exp(-((x - m)^2/(2 * sd))) * 
        (x - m)^2/(2 * sd)^2)))))
   expect_equal(dp(0, a, m, s), -0.9547048, tolerance=1.e-6)
})
drule[["myfun"]] <- NULL

# test renaming primitive function (issue #10)
f=cos
g = function(f) Deriv(f)
test_that("renaming primitive", {
   expect_identical(g(f), Deriv(cos))
})

# test returning a constant vector of length > 1 and c()-argument  (issues #14 and #15)
f <- function(x, y) x + y
res=c(x=1, y=1)
fd=as.function(alist(x=, y=, res))
body(fd)=res
f2=function(x, y) c(x, y)^2
test_that("multivar diff", {
   expect_identical(Deriv(f), fd)
   expect_equal(Deriv(f2, cache=FALSE), function (x, y) c(x = c(2, 0) * c(x, y), y = c(0, 2) * c(x, y)), check.environment=FALSE)
})

Sys.setlocale(category = "LC_COLLATE", locale = lc_orig)
