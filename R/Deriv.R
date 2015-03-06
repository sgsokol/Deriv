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
		if (is.call(b) && b[[1]] == as.symbol(".Internal")) {
			fch <- deparse(substitute(f))
			if (fch %in% dlin || !is.null(drule[[fch]])) {
				res <- Deriv_(as.call(c(substitute(f), lapply(x, as.symbol))), x, env, use.D)
				if (cache.exp)
					res <- Cache(res)
				as.function(c(formals(f), res), envir=env)
			} else {
				stop(sprintf("Internal function '%s()' is not in derivative table.", fch))
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
			dargs <- lapply(as.list(st)[-1], function(a) Deriv_(a, x, env, use.D))
			return(Simplify_(as.call(c(st[[1]], dargs))))
		}
		nb_args=length(st)-1
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
		}
		if (is.null(drule[[stch]])) {
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
		} else if (is.null(drule[[stch]][[nb_args]]) && !use.D) {
			stop(sprintf("Don't know how to differentiate function or operator '%s' when it is called with %s arguments", stch, nb_args))
		}
		# there is a rule!
		if (use.D) {
			return(Simplify(D(st, x)))
		}
		# prepare replacement list ._1 -> first argument, ._d1 -> derivative of the first argument and so on
		args <- as.list(st)[-1]
		names(args) <- paste("._", seq_len(nb_args), sep="")
		# which arguments have to be differentiated?
		dstr <- format1(drule[[stch]][[nb_args]])
		ma <- gregexpr("\\._d[1-9][0-9]*", dstr)
		dgrep <- substring(dstr, ma[[1]], ma[[1]]+attr(ma[[1]], "match.length")-1)
		if (any(nchar(dgrep) != 0)) {
			dnum <- as.integer(sub("._d", "", unique(dgrep)))
			dargs <- lapply(args[dnum], Deriv_, x, env, use.D)
#cat("st=", format1(st), "\n", sep="")
#cat("stch=", stch, "\n", sep="")
#cat("dgrep=", dgrep, "\n", sep=", ")
#cat("dstr=", dstr, "\n", sep="")
#cat("dnum=", dnum, "\n", sep=", ")
			names(dargs) <- paste("._d", dnum, sep="")
		} else {
			dargs <- NULL
		}
		return(Simplify_(do.call("substitute", list(drule[[stch]][[nb_args]], c(args, dargs)))))
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

drule <- new.env()
dsym <- new.env()

# linear functions, i.e. d(f(x))/dx == f(d(arg)/dx)
dlin=c("+", "-", "c", "t", "sum", "cbind", "rbind")

# first item in the list correspond to a call with one argument
# second (if any) for two, third for three. NULL means that with
# this number of argument this function can not be called
drule[["("]] <- list(quote(._d1)) # (._1) => omit paranthesis
# linear arithmetics are already handled by dlin
# exception is maid for unitary plus (it is just omitted)
#drule[["+"]] <- list(quote(._d1), quote(._d1+._d2)) # +._1, ._1+._2
#drule[["-"]] <- list(quote(-._d1), quote(._d1-._d2)) # -._1, ._1-._2
# arithmetic non linear rules
drule[["*"]] <- list(NULL, quote(._d1*._2+._1*._d2)) # ._1*._2
drule[["/"]] <- list(NULL, quote(._d1/._2-._1*._d2/._2^2)) # ._1*._2
# power functions
drule[["^"]] <- list(NULL, quote(._d1*._2*._1^(._2-1)+._d2*log(._1)*._1^._2)) # ._1^._2
# example of recursive call
#drule[["sqrt"]] <- list(Deriv(call("^", as.symbol("._1"), 0.5), "._1", NULL, use.D)) # sqrt(._1)
# but we prefer a sqrt() formula
drule[["sqrt"]] <- list(quote(0.5*._d1/sqrt(._1)))
drule[["log"]] <- list(quote(._d1/._1), quote(._d1/(._1*log(._2))- ._d2*log(._2, ._1)/(._2*log(._2)))) # log(._1), log(._1, b)
drule[["logb"]] <- drule[["log"]]
drule[["log2"]] <- list(quote(._d1/(._1*log(2))))
drule[["log10"]] <- list(quote(._d1/(._1*log(10))))
drule[["log1p"]] <- list(quote(._d1/(._1+1)))
drule[["exp"]] <- list(quote(._d1*exp(._1)))
drule[["expm1"]] <- list(quote(._d1*exp(._1)))
# trigonometric
drule[["sin"]] <- list(quote(._d1*cos(._1)))
drule[["cos"]] <- list(quote(-._d1*sin(._1)))
drule[["tan"]] <- list(quote(._d1/cos(._1)^2))
drule[["asin"]] <- list(quote(._d1/sqrt(1-._1^2)))
drule[["acos"]] <- list(quote(-._d1/sqrt(1-._1^2)))
drule[["atan"]] <- list(quote(._d1/(1+._1^2)))
drule[["atan2"]] <- list(NULL, quote((._d1*._2-._d2*._1)/(._1^2+._2^2)))
drule[["sinpi"]] <- list(quote(pi*._d1*cospi(._1)))
drule[["cospi"]] <- list(quote(-pi*._d1*sinpi(._1)))
drule[["tanpi"]] <- list(quote(pi*._d1/cospi(._1)^2))
# hyperbolic
drule[["sinh"]] <- list(quote(._d1*cosh(._1)))
drule[["cosh"]] <- list(quote(._d1*sinh(._1)))
drule[["tanh"]] <- list(quote(._d1*(1-tanh(._1)^2)))
drule[["asinh"]] <- list(quote(._d1/sqrt(._1^2+1)))
drule[["acosh"]] <- list(quote(._d1/sqrt(._1^2-1)))
drule[["atanh"]] <- list(quote(._d1/(1-._1^2)))
# code control
drule[["if"]] <- list(NULL, quote(if (._1) ._d2), quote(if (._1) ._d2 else ._d3))
# sign depending functions
drule[["abs"]] <- list(quote(._d1*sign(._1)))
drule[["sign"]] <- list(0)
# special functions
drule[["besselI"]] <- list(NULL,
	quote(if (._2 == 0) ._d1*besselI(._1, 1) else 0.5*._d1*(besselI(._1, ._2-1) + besselI(._1, ._2+1))),
	quote((if (._2 == 0) ._d1*besselI(._1, 1, ._3) else 0.5*._d1*(besselI(._1, ._2-1, ._3) + besselI(._1, ._2+1, ._3)))-if (._3) besselI(._1, ._2, TRUE) else 0)
)
drule[["besselK"]] <- list(NULL,
	quote(if (._2 == 0) -._d1*besselK(._1, 1) else -0.5*._d1*(besselK(._1, ._2-1) + besselK(._1, ._2+1))),
	quote((if (._2 == 0) -._d1*besselK(._1, 1, ._3) else -0.5*._d1*(besselK(._1, ._2-1, ._3) + besselK(._1, ._2+1, ._3)))+if (._3) besselK(._1, ._2, TRUE) else 0)
)
drule[["besselJ"]] <- list(NULL,
	quote(if (._2 == 0) -._d1*besselJ(._1, 1) else 0.5*._d1*(besselJ(._1, ._2-1) - besselJ(._1, ._2+1)))
)
drule[["besselY"]] <- list(NULL,
	quote(if (._2 == 0) -._d1*besselY(._1, 1) else 0.5*._d1*(besselY(._1, ._2-1) - besselY(._1, ._2+1)))
)
drule[["gamma"]] <- list(quote(._d1*gamma(._1)*digamma(._1)))
drule[["lgamma"]] <- list(quote(._d1*digamma(._1)))
drule[["digamma"]] <- list(quote(._d1*trigamma(._1)))
drule[["trigamma"]] <- list(quote(._d1*psigamma(._1, 2L)))
drule[["psigamma"]] <- list(quote(._d1*psigamma(._1, 1L)), quote(._d1*psigamma(._1, ._2+1L)))
drule[["beta"]] <- list(NULL, quote(beta(._1, ._2)*(._d1*digamma(._1)+._d2*digamma(._2)-digamma(._1+._2)*(._d1+._d2))))
drule[["lbeta"]] <- list(NULL, quote(._d1*digamma(._1)+._d2*digamma(._2)-digamma(._1+._2)*(._d1+._d2)))
