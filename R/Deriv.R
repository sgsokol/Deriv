#' @name Deriv
#' @title Symbollic differentiation of an expression or function
#' @aliases Deriv drule
#' @concept symbollic differentiation
# \usage{
# Deriv(f, x=if (is.function(f)) names(formals(f)) else all.vars(if (is.character(f)) parse(text=f) else f), env=if (is.function(f)) environment(f) else parent.frame(), use.D=FALSE, cache.exp=TRUE)
# }
#' 
#' 
#' @param f An expression or function to be differentiated.
#'  f can be \itemize{
#'   \item a user defined function: \code{function(x) x**n}
#'   \item a string: \code{"x**n"}
#'   \item an expression: \code{expression(x**n)}
#'   \item a call: \code{call("^", quote(x), quote(n))}
#'   \item a language: \code{quote(x**n)}
#'   \item a right hand side of a formula: \code{~ x**n} or \code{y ~ x**n}
#'  }
#' @param x An optiona character vector with variable name(s) with resptect to which
#'  \code{f} must be differentiated. If not provided, x is guessed from
#'  \code{names(formals(f))}, if \code{f} is a function, or from all variables in f
#'  in other cases. If f is a primitive
#'  function, x is set to \code{names(formals(args(f)))}
#' @param env An environment where the symbols and functions are searched for.
#'  Defaults to \code{parent.frame()} for \code{f} expression and to
#'  \code{environment(f)} if \code{f} is a function. For primitive function,
#'  it is set by default to .GlobalEnv
#' @param use.D An optional logical (default FALSE), indicates if base::D()
#'  must be used for differentiation of basic expressions.
#' @param cache.exp An optional logical (default TRUE), indicates if
#'  final expression must be optimized with cached subexpressions.
#'  If enabled, repeated calculations are made only once and their
#'  results stored in cache variables which are then reused.
#' 
#' @return \itemize{
#'  \item a function if \code{f} is a function
#'  \item an expression if \code{f} is an expression
#'  \item a character string if \code{f} is a character string
#'  \item a language (usually a so called 'call' but may be also a symbol or just a numeric) for other types of \code{f}
#' }
#'
#' @details
#' R already contains two differentiation functions: D and deriv. D does
#' simple univariate differentiation.  "deriv" uses D do to multivariate
#' differentiation.  The output of "D" is an expression, whereas the output of
#' "deriv" can be an executable function.
#' 
#' R's existing functions have several limitations.  They can probably be fixed,
#' but since they are written in C, this would probably require a lot of work.
#' Limitations include:
#' \itemize{
#'  \item The derivatives table can't be modified at runtime, and is only available
#' in C.
#'  \item Function cannot substitute function calls.  eg:
#'	f <- function(x, y) x + y; deriv(~f(x, x^2), "x")
#' }
#'
#' So, here are the advantages of this implementation:
#' 
#' \itemize{
#'  \item It is entirely written in R, so would be easier to maintain.
#'  \item Can do multi-variate differentiation.
#'  \item Can differentiate function calls:
#'  \itemize{
#'	   \item if the function is in the derivative table, then the chain rule
#'	is applied.  For example, if you declared that the derivative of
#'	sin is cos, then it would figure out how to call cos correctly.
#'	   \item if the function is not in the derivative table (or it is anonymous),
#'	then the function body is substituted in.
#'	   \item these two methods can be mixed.  An entry in the derivative table
#'	need not be self-contained -- you don't need to provide an infinite
#'	chain of derivatives.
#'  }
#'  \item It's easy to add custom entries to the derivatives table, e.g.
#'   
#'   \code{drule[["cos"]] <- list(-._d1*sin(._1))}
#'  \item The output is an executable function, which makes it suitable
#'      for use in optimization problems.
#'  \item Compound functions (i.e. piece-wise functions based on if-else operator) can
#'      be differentiated (cf. examples section).
#' }
#' 
#' Two working environments \code{drule} and \code{simplifications} are
#' exported in the package namescape.
#' As their names indicate, they contain tables of derivative and
#' simplification rules.
#' To see the list of defined rules do \code{ls(drule)}.
#' To add your own derivative rule for a function called say \code{sinpi(x)} calculating sin(pi*x), do \code{drule[["sinpi"]] <- list(quote(pi*._d1*cospi(._1)))}.
#' Here, "._1" stands for the "first arguments", "._d1" for the first derivative of the first arguments. For a function that might have more than one argument,
#' e.g. \code{log(x, base=exp(1))}, the drule entry must be a list with one rule
#' per number of possible arguments. See \code{drule$log} for an example to follow.
#' After adding \code{sinpi} you can differentiate expressions like \code{Deriv(~ sinpi(x^2), "x")}
#' 
#' NB. In \code{abs()} and \code{sign()} function, singularity treatment at point 0 is left to user's care.
#' 
#' NB2. In Bessel functions, derivatives are calculated only by the first argument,
#'      not by the \code{nu} argument which is supposed to be constant.
#' @examples
#'
#' \dontrun{f <- function(x) x^2}
#' \dontrun{Deriv(f)}
#' # function (x)
#' # 2 * x
#' 
#' \dontrun{f <- function(x, y) sin(x) * cos(y)}
#' \dontrun{Deriv(f)}
#' # function (x, y)
#' # c(x = cos(x) * cos(y), y = -(sin(x) * sin(y)))
#'
#' \dontrun{f_ <- Deriv(f)}
#' \dontrun{f_(3, 4)}
#' #              x         y
#' # [1,] 0.6471023 0.1068000
#' 
#' \dontrun{Deriv(~ f(x, y^2), "y")}
#' # -(2 * (y * sin(x) * sin(y^2)))
#' 
#' \dontrun{Deriv(quote(f(x, y^2)), c("x", "y"))}
#' # c(x = cos(x) * cos(y^2), y = -(2 * (y * sin(x) * sin(y^2))))
#' 
#' \dontrun{Deriv(expression(sin(x^2) * y), "x")}
#' # expression(2*(x*y*cos(x^2)))
#' 
#' Deriv("sin(x^2) * y", "x") # differentiate only by x
#' "2 * (x * y * cos(x^2))"
#' 
#' Deriv("sin(x^2) * y", cache.exp=FALSE) # differentiate by all variables (here by x and y)
#' "c(x = 2 * (x * y * cos(x^2)), y = sin(x^2))"
#' 
#' # Compound function example (here abs(x) smoothed near 0)
#' fc <- function(x, h=0.1) if (abs(x) < h) 0.5*h*(x/h)**2 else abs(x)-0.5*h
#' Deriv("fc(x)", "x", cache.exp=FALSE)
#' "if (abs(x) < h) x/h else sign(x)"
#' 
#' # Example of a first argument that cannot be evaluated in the current environment:
#' \dontrun{
#'   suppressWarnings(rm("xx", "yy"))
#'   Deriv(xx^2+yy^2)
#' }
#' # c(xx = 2 * xx, yy = 2 * yy)
#' 
#' # Automatic differentiation (AD), note itermediate variable \code{d} assignment
#' \dontrun{Deriv(~{d <- ((x-m)/s)^2; exp(-0.5*d)}, "x")}
#' #{
#' #   d <- ((x - m)/s)^2
#' #   .d_x <- 2 * ((x - m)/s^2)
#' #   -(0.5 * (.d_x * exp(-(0.5 * d))))
#' #}
#' 

# This wrapper of Deriv_ returns an appropriate expression (wrapped up "properly") depending on the type of argument to differentiate
Deriv <- function(f, x=if (is.function(f)) names(formals(f)) else all.vars(if (is.character(f)) parse(text=f) else f), env=if (is.function(f)) environment(f) else parent.frame(), use.D=FALSE, cache.exp=TRUE) {
	tf <- try(f, silent=TRUE)
	if (inherits(tf, "try-error")) {
		f <- substitute(f)
	}
	# clean dsym env
	rm(list=ls(dsym), envir=dsym)
	if (is.null(x)) {
		# primitive function
		fch <- deparse(substitute(f))
		af <- formals(args(f))
		x <- names(af)
		if ("..." %in% x) {
			stop(sprintf("Undefined list of arguments for %s()", fch))
		}
		if (is.null(env))
			env <- .GlobalEnv
		rule <- drule[[fch]]
		if (is.null(rule) && !fch %in% dlin) {
			stop(sprintf("Undefined rule for '%s()' differentiation", fch))
		}
		res <- Deriv_(as.call(c(as.symbol(fch), lapply(x, as.symbol))), x, env, use.D)
		if (cache.exp)
			res <- Cache(res)
		return(as.function(c(af, res), envir=env))
	}
	x <- as.character(x)
	if (any(nchar(x) == 0)) {
		stop("Names in the second argument must not be empty")
	}
	if (is.character(f)) {
		# f is to parse
		res <- Deriv_(parse(text=f)[[1]], x, env, use.D)
		if (cache.exp)
			res <- Cache(res)
		format1(res)
	} else if (is.function(f)) {
		b <- body(f)
		if (is.call(b) && b[[1]] == as.symbol(".Internal") || b[[1]] == as.symbol(".External")) {
			fch <- deparse(substitute(f))
			if (fch %in% dlin || !is.null(drule[[fch]])) {
				arg <- lapply(names(formals(args(f))), as.symbol)
				acall <- as.call(c(as.symbol(fch), arg))
				res <- Deriv_(acall, x, env, use.D)
				if (cache.exp)
					res <- Cache(res)
				as.function(c(formals(f), res), envir=env)
			} else {
				stop(sprintf("Internal or external function '%s()' is not in derivative table.", fch))
			}
		} else {
			res <- Deriv_(b, x, env, use.D)
			if (cache.exp)
				res <- Cache(res)
			as.function(c(formals(f), res), envir=env)
		}
	} else if (is.expression(f)) {
		res <- Deriv_(f[[1]], x, env, use.D)
		if (cache.exp)
			res <- Cache(res)
		as.expression(res)
	} else if (is.language(f)) {
		if (is.call(f) && f[[1]] == as.symbol("~")) {
			# rhs of the formula
			res <- Deriv_(f[[length(f)]], x, env, use.D)
			if (cache.exp)
				res <- Cache(res)
			res
		} else {
			# plain call derivation
			res <- Deriv_(f, x, env, use.D)
			if (cache.exp)
				res <- Cache(res)
			res
		}
	} else {
		f <- substitute(f)
		if (length(x) == 0)
			x <- all.vars(f)
		res <- Deriv_(f, x, env, use.D)
		if (cache.exp)
			res <- Cache(res)
		res
		#stop("Invalid type of 'f' for differentiation")
	}
}

# workhorse function doing the main work of differentiation
Deriv_ <- function(st, x, env, use.D) {
	# Make x scalar and wrap results in a c() call if length(x) > 1
	if (length(x) > 1) {
		# many variables => recursive call on single name
		res <- lapply(x, function(xi) Deriv_(st, xi, env, use.D))
		names(res) <- x;
		return(as.call(c(as.symbol("c"), res)))
	}
	# differentiate R statement 'st' (a call, or a symbol or numeric) by a name in 'x'
	if (is.unumeric(st)) {
		return(0)
	} else if (is.symbol(st)) {
		stch <- as.character(st)
		if (stch == x) {
			return(1)
		} else if (is.null(dsym[[stch]])) {
			return(0)
		} else {
			return(dsym[[stch]])
		}
	} else if ((is.uminus(st) || is.uplus(st)) && is.symbol(st[[2]])) {
		if (as.character(st[[2]]) == x) {
			return(if (is.uminus(st)) -1 else 1)
		} else {
			return(0)
		}
	} else if (is.call(st)) {
#browser()
		stch <- as.character(st[[1]])
		args <- as.list(st)[-1]
		if (stch %in% dlin) {
			# linear case
			# differentiate all arguments then pass them to the function
			dargs <- lapply(args, Deriv_, x, env, use.D)
			return(Simplify_(as.call(c(st[[1]], dargs))))
		}
		nb_args=length(st)-1
		# special cases: out of rule table or args(stch) -> NULL
		if (stch == "{") {
#browser()
			# AD differentiation
			res=list(st[[1]])
			for (a in args) {
				if (is.call(a) && (a[[1]] == as.symbol("<-") || a[[1]] == as.symbol("="))) {
					if (!is.symbol(a[[2]]))
						stop(sprintf("In AD mode, don't know to deal with a non symbol '%s' at lhs", format1(a[[2]])))
					res <- append(res, a)
					de_a <- Deriv_(a[[3]], x, env, use.D)
					ach <- as.character(a[[2]])
					if (de_a == 0) {
						next
					} else if (!is.call(de_a)) {
						dsym[[ach]] <- de_a
						next
					}
					d_a <- as.symbol(paste(".", ach, "_", x, sep=""))
					dsym[[ach]] <- d_a
					res <- append(res, call("<-", d_a, de_a))
				} else {
					res <- append(res, Deriv_(a, x, env, use.D))
				}
			}
			return(as.call(res))
		} else if (is.uminus(st)) {
			return(Simplify(call("-", Deriv_(st[[2]], x, env, use.D))))
		} else if (stch == "(") {
			return(Simplify(Deriv_(st[[2]], x, env, use.D)))
		} else if (stch == "if") {
			return(if (nb_args == 2)
				Simplify(call("if", st[[2]], Deriv_(st[[3]], x, env, use.D))) else
				Simplify(call("if", st[[2]], Deriv_(st[[3]], x, env, use.D),
					Deriv_(st[[4]], x, env, use.D))))
		}
		rule <- drule[[stch]]
		if (is.null(rule)) {
			# no derivative rule for this function
			# try to get the body and differentiate it
			ff <- get(stch, mode="function", envir=env)
			bf <- body(ff)
			if (is.null(bf)) {
				stop(sprintf("Could not retrieve body of '%s()'", stch))
			}
			mc <- match.call(ff, st)
			st <- Simplify_(do.call("substitute", list(bf, as.list(mc)[-1])))
			return(Deriv_(st, x, env, use.D))
		}
		# there is a rule!
		if (use.D) {
			return(Simplify(D(st, x)))
		}
		# prepare replacement list
		da <- args(stch)
		mc <- as.list(match.call(def=da, call=st))[-1]
		da <- as.list(da)
		da <- da[-length(da)] # all declared arguments with default values
		aa <- modifyList(da, mc) # all arguments with actual values
		# actualize the rule with actual arguments
		rule <- lapply(rule, function(r) do.call("substitute", list(r, aa)))
#browser()		
		# which arguments have to be differentiated?
		iad <- which(!sapply(rule, is.null))
		rule <- rule[iad]
		if (!any(names(which(mc==x)) == names(rule))) {
			warning(sprintf("A call %s cannot be differentiated by the argument '%s'", format1(st), x))
			return(NULL)
		}
		# dargs are ordered by rule names
		dargs <- lapply(mc[names(iad)], Deriv_, x, env, use.D)
		ize <- sapply(dargs, `==`, 0)
		dargs <- dargs[!ize]
		rule <- rule[!ize]
		if (length(rule) == 0) {
			return(0)
		}
		
		# apply chain rule where needed
		ione <- sapply(dargs, `==`, 1)
		imone <- sapply(dargs, `==`, -1)
		for (i in seq_along(rule)[!(ione|imone)]) {
			rule[[i]] <- call("*", dargs[[i]], rule[[i]])
		}
		for (i in seq_along(rule)[imone]) {
			rule[[i]] <- call("-", rule[[i]])
		}
		return(Simplify(li2sum(rule)))
	} else if (is.function(st)) {
		# differentiate its body if can get it
		args <- as.list(st)[-1]
		names(args)=names(formals(ff))
		if (is.null(names(args))) {
			stop(sprintf("Could not retrieve arguments of '%s()'", stch))
		}
		st <- do.call("substitute", list(body(ff), args))
		Deriv_(st, x, env, use.D)
	} else {
		stop("Invalid type of 'st' argument. It must be numeric, symbol or a call.")
	}
}

sl <- function(...) {
	# substitute arguments and return the list
	mc <- match.call()
	as.list(mc)[-1]
}
drule <- new.env()
dsym <- new.env()

# linear functions, i.e. d(f(x))/dx == f(d(arg)/dx)
dlin=c("+", "-", "c", "t", "sum", "cbind", "rbind")


drule[["*"]] <- sl(e1=e2, e2=e1)
drule[["^"]] <- sl(e1=e2*e1^(e2-1), e2=e1^e2/log(e1))
drule[["/"]] <- sl(e1=1/e2, e2=-e1/e2^2)
drule[["sqrt"]] <- sl(x=0.5/sqrt(x))
drule[["log"]] <- sl(x=1/(x*log(base)), base=-log(base, x)/(base*log(base)))
drule[["logb"]] <- drule[["log"]]
drule[["log2"]] <- sl(x=1/(x*log(2)))
drule[["log10"]] <- sl(x=1/(x*log(10)))
drule[["log1p"]] <- sl(x=1/(x+1))
drule[["exp"]] <- sl(x=exp(x))
drule[["expm1"]] <- sl(x=exp(x))
# trigonometric
drule[["sin"]] <- sl(x=cos(x))
drule[["cos"]] <- sl(x=-sin(x))
drule[["tan"]] <- sl(x=1/cos(x)^2)
drule[["asin"]] <- sl(x=1/sqrt(1-x^2))
drule[["acos"]] <- sl(x=-1/sqrt(1-x^2))
drule[["atan"]] <- sl(x=1/(1+x^2))
drule[["atan2"]] <- sl(y=x/(x^2+y^2), x=-y/(x^2+y^2))
drule[["sinpi"]] <- sl(x=pi*cospi(x))
drule[["cospi"]] <- sl(x=-pi*sinpi(x))
drule[["tanpi"]] <- sl(x=pi/cospi(x)^2)
# hyperbolic
drule[["sinh"]] <- sl(x=cosh(x))
drule[["cosh"]] <- sl(x=sinh(x))
drule[["tanh"]] <- sl(x=(1-tanh(x)^2))
drule[["asinh"]] <- sl(x=1/sqrt(x^2+1))
drule[["acosh"]] <- sl(x=1/sqrt(x^2-1))
drule[["atanh"]] <- sl(x=1/(1-x^2))
# sign depending functions
drule[["abs"]] <- sl(x=sign(x))
drule[["sign"]] <- sl(x=0)
# special functions
drule[["besselI"]] <- sl(x=if (nu == 0) besselI(x, 1, expon.scaled) else 0.5*(besselI(x, nu-1, expon.scaled) + besselI(x, nu+1, expon.scaled))-if (expon.scaled) besselI(x, nu, TRUE) else 0, nu=NULL, expon.scaled=NULL)
drule[["besselK"]] <- sl(x=if (nu == 0) -besselK(x, 1, expon.scaled) else -0.5*(besselK(x, nu-1, expon.scaled) + besselK(x, nu+1, expon.scaled))+if (expon.scaled) besselK(x, nu, TRUE) else 0)
drule[["besselJ"]] <- sl(x=if (nu == 0) -besselJ(x, 1) else 0.5*(besselJ(x, nu-1) - besselJ(x, nu+1)), nu=NULL)
drule[["besselY"]] <- sl(x=if (nu == 0) -besselY(x, 1) else 0.5*(besselY(x, nu-1) - besselY(x, nu+1)), nu=NULL)
drule[["gamma"]] <- sl(x=gamma(x)*digamma(x))
drule[["lgamma"]] <- sl(x=digamma(x))
drule[["digamma"]] <- sl(x=trigamma(x))
drule[["trigamma"]] <- sl(x=psigamma(x, 2L))
drule[["psigamma"]] <- sl(x=psigamma(x, deriv+1L), deriv=NULL)
drule[["beta"]] <- sl(a=beta(a, b)*(digamma(a)-digamma(a+b)), b=beta(a, b)*(digamma(b)-digamma(a+b)))
drule[["lbeta"]] <- sl(a=digamma(a)-digamma(a+b), b=digamma(b)-digamma(a+b))
# probability densities
drule[["dbinom"]] <- sl(x=NULL, size=NULL, log=NULL, prob=if (size == 0) -x*(1-prob)^(x-1) else if (x == size) size*prob^(size-1) else (size-x*prob)*(x-size+1)*dbinom(x, size-1, prob)/(1-prob)^2/(if (log) dbinom(x, size, prob) else 1))
drule[["dnorm"]] <- sl(x=-(x-mean)/sd^2*if (log) 1 else dnorm(x, mean, sd),
	mean=(x-mean)/sd^2*if (log) 1 else dnorm(x, mean, sd),
	sd=(((x - mean)/sd)^2 - 1)/sd * if (log) 1 else dnorm(x, mean, sd),
	log=NULL)
