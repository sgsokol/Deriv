#' @name Deriv
#' @title Symbollic differentiation of an expression or function
#' @aliases Deriv drule qlist
#' @concept symbollic differentiation
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
#' @param ... (in \code{qlist()}) is a suite of named unevaluated expressions.
#'  It is used to add derivative rules to \code{drule} environment.
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
#'   \code{drule[["cos"]] <- qlist(x=-sin(x))}
#'   
#'   The chain rule will be automatically applied if needed.
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
#' To add your own derivative rule for a function called say \code{sinpi(x)} calculating sin(pi*x), do \code{drule[["sinpi"]] <- qlist(x=pi*cospi(x))}.
#' Here, "x" stands for the first and unique argument in \code{sinpi()} definition. For a function that might have more than one argument,
#' e.g. \code{log(x, base=exp(1))}, the drule entry must be a list with a named rule
#' per argument. See \code{drule$log} for an example to follow.
#' After adding \code{sinpi} you can differentiate expressions like \code{Deriv(~ sinpi(x^2), "x")}. The chain rule will automatically apply.
#' 
#' NB. In \code{abs()} and \code{sign()} function, singularity treatment at point 0 is left to user's care.
#' 
#' NB2. In Bessel functions, derivatives are calculated only by the first argument,
#'      not by the \code{nu} argument which is supposed to be constant.
#' @author Andrew Clausen (original version) and Serguei Sokol (maintainer)
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
#' \dontrun{Deriv(quote(f(x, y^2)), c("x", "y"), cache.exp=FALSE)}
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
#' # Custom derivation rule
#' \dontrun{
#'   myfun <- function(x, y=TRUE) NULL # do something usefull
#'   dmyfun <- function(x, y=TRUE) NULL # myfun derivative by x.
#'   drule[["myfun"]] <- qlist(x=dmyfun(x, y), y=NULL) # y is just a logical
#'   Deriv(myfun(z^2, FALSE), "z")
#'   # 2 * (z * dmyfun(z^2, FALSE))
#' }

Deriv <- function(f, x=if (is.function(f)) names(formals(f)) else all.vars(if (is.character(f)) parse(text=f) else f), env=if (is.function(f)) environment(f) else parent.frame(), use.D=FALSE, cache.exp=TRUE) {
	tf <- try(f, silent=TRUE)
	if (inherits(tf, "try-error")) {
		f <- substitute(f)
	}
	# clean dsym env
	rm(list=ls(dsym), envir=dsym)
	if (is.null(env))
		env <- .GlobalEnv
	if (is.null(x)) {
		# primitive function
		fch <- deparse(substitute(f))
		af <- formals(args(f))
		x <- names(af)
		if ("..." %in% x) {
			stop(sprintf("Undefined list of arguments for %s()", fch))
		}
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
		if ((is.call(b) && (b[[1]] == as.symbol(".Internal") || b[[1]] == as.symbol(".External"))) || (is.null(b) && is.primitive(f))) {
			fch <- deparse(substitute(f))
			if (fch %in% dlin || !is.null(drule[[fch]])) {
				arg <- lapply(names(formals(args(f))), as.symbol)
				acall <- as.call(c(as.symbol(fch), arg))
				res <- Deriv_(acall, x, env, use.D)
				if (cache.exp)
					res <- Cache(res)
				as.function(c(formals(args(f)), res), envir=env)
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
		res <- lapply(x, function(xi) {rm(list=ls(dsym), envir=dsym); Deriv_(st, xi, env, use.D)})
		names(res) <- x;
		return(as.call(c(as.symbol("c"), res)))
	}
	# differentiate R statement 'st' (a call, or a symbol or numeric) by a name in 'x'
	if (is.conuloch(st)) {
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
						stop(sprintf("In AD mode, don't know how to deal with a non symbol '%s' at lhs", format1(a[[2]])))
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
		mc <- as.list(match.call(definition=da, call=st))[-1]
		da <- as.list(da)
		da <- da[-length(da)] # all declared arguments with default values
		aa <- modifyList(da, mc) # all arguments with actual values
		# actualize the rule with actual arguments
		rule <- lapply(rule, function(r) do.call("substitute", list(r, aa)))
#browser()		
		# which arguments have to be differentiated?
		iad <- which(!sapply(rule, is.null))
		rule <- rule[iad]
		lsy <- ls(dsym, all.names=TRUE)
		if (!any(names(which(sapply(mc, function(it) {av <- all.vars(it); any(x == av) || any(av %in% lsy)}))) == names(rule))) {
			#warning(sprintf("A call %s cannot be differentiated by the argument '%s'", format1(st), x))
			return(0)
		}
		# rules and dargs are ordered by mc names
		rule <- rule[names(mc)]
		dargs <- lapply(mc, Deriv_, x, env, use.D)
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
		stop("Invalid type of 'st' argument. It must be constant, symbol or a call.")
	}
}

#' @rdname Deriv
qlist <- function(...) {
	# substitute arguments and return the list
	mc <- match.call()
	as.list(mc)[-1]
}
drule <- new.env()
dsym <- new.env()

# linear functions, i.e. d(f(x))/dx == f(d(arg)/dx)
dlin=c("+", "-", "c", "t", "sum", "cbind", "rbind")


drule[["*"]] <- qlist(e1=e2, e2=e1)
drule[["^"]] <- qlist(e1=e2*e1^(e2-1), e2=e1^e2*log(e1))
drule[["/"]] <- qlist(e1=1/e2, e2=-e1/e2^2)
drule[["sqrt"]] <- qlist(x=0.5/sqrt(x))
drule[["log"]] <- qlist(x=1/(x*log(base)), base=-log(base, x)/(base*log(base)))
drule[["logb"]] <- drule[["log"]]
drule[["log2"]] <- qlist(x=1/(x*log(2)))
drule[["log10"]] <- qlist(x=1/(x*log(10)))
drule[["log1p"]] <- qlist(x=1/(x+1))
drule[["exp"]] <- qlist(x=exp(x))
drule[["expm1"]] <- qlist(x=exp(x))
# trigonometric
drule[["sin"]] <- qlist(x=cos(x))
drule[["cos"]] <- qlist(x=-sin(x))
drule[["tan"]] <- qlist(x=1/cos(x)^2)
drule[["asin"]] <- qlist(x=1/sqrt(1-x^2))
drule[["acos"]] <- qlist(x=-1/sqrt(1-x^2))
drule[["atan"]] <- qlist(x=1/(1+x^2))
drule[["atan2"]] <- qlist(y=x/(x^2+y^2), x=-y/(x^2+y^2))
drule[["sinpi"]] <- qlist(x=pi*cospi(x))
drule[["cospi"]] <- qlist(x=-pi*sinpi(x))
drule[["tanpi"]] <- qlist(x=pi/cospi(x)^2)
# hyperbolic
drule[["sinh"]] <- qlist(x=cosh(x))
drule[["cosh"]] <- qlist(x=sinh(x))
drule[["tanh"]] <- qlist(x=(1-tanh(x)^2))
drule[["asinh"]] <- qlist(x=1/sqrt(x^2+1))
drule[["acosh"]] <- qlist(x=1/sqrt(x^2-1))
drule[["atanh"]] <- qlist(x=1/(1-x^2))
# sign depending functions
drule[["abs"]] <- qlist(x=sign(x))
drule[["sign"]] <- qlist(x=0)
# special functions
drule[["besselI"]] <- qlist(x=(if (nu == 0) besselI(x, 1, expon.scaled) else 0.5*(besselI(x, nu-1, expon.scaled) + besselI(x, nu+1, expon.scaled)))-if (expon.scaled) besselI(x, nu, TRUE) else 0, nu=NULL, expon.scaled=NULL)
drule[["besselK"]] <- qlist(x=(if (nu == 0) -besselK(x, 1, expon.scaled) else -0.5*(besselK(x, nu-1, expon.scaled) + besselK(x, nu+1, expon.scaled)))+if (expon.scaled) besselK(x, nu, TRUE) else 0)
drule[["besselJ"]] <- qlist(x=if (nu == 0) -besselJ(x, 1) else 0.5*(besselJ(x, nu-1) - besselJ(x, nu+1)), nu=NULL)
drule[["besselY"]] <- qlist(x=if (nu == 0) -besselY(x, 1) else 0.5*(besselY(x, nu-1) - besselY(x, nu+1)), nu=NULL)
drule[["gamma"]] <- qlist(x=gamma(x)*digamma(x))
drule[["lgamma"]] <- qlist(x=digamma(x))
drule[["digamma"]] <- qlist(x=trigamma(x))
drule[["trigamma"]] <- qlist(x=psigamma(x, 2L))
drule[["psigamma"]] <- qlist(x=psigamma(x, deriv+1L), deriv=NULL)
drule[["beta"]] <- qlist(a=beta(a, b)*(digamma(a)-digamma(a+b)), b=beta(a, b)*(digamma(b)-digamma(a+b)))
drule[["lbeta"]] <- qlist(a=digamma(a)-digamma(a+b), b=digamma(b)-digamma(a+b))
# probability densities
drule[["dbinom"]] <- qlist(x=NULL, size=NULL, log=NULL, prob=if (size == 0) -x*(1-prob)^(x-1) else if (x == size) size*prob^(size-1) else (size-x*prob)*(x-size+1)*dbinom(x, size-1, prob)/(1-prob)^2/(if (log) dbinom(x, size, prob) else 1))
drule[["dnorm"]] <- qlist(x=-(x-mean)/sd^2*if (log) 1 else dnorm(x, mean, sd),
	mean=(x-mean)/sd^2*if (log) 1 else dnorm(x, mean, sd),
	sd=(((x - mean)/sd)^2 - 1)/sd * if (log) 1 else dnorm(x, mean, sd),
	log=NULL)
