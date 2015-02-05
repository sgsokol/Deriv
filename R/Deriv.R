#' @name Deriv
#' @title Symbollic differentiation of an expression or function
#' @aliases Deriv drule
#' @concept symbollic differentiation
# \usage{
# Deriv(f, x=if (is.function(f)) names(formals(f)) else "", env=if (is.function(f)) environment(f) else parent.frame())
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
#' @param x A character string with variable name(s) with resptect to which
#'  \code{f} must be differentiated. If \code{f} is a function \code{x} is
#'  optional andis set to \code{names(formals(f))}. If f is a primitive
#'  function, x is set to c("x" [, "y"[, "z"]]) depending on the maximal
#'  number of allowed arguments in drule table. For functions with
#'  unlimited number of arguments, like sum(), three names are used:
#'  "x", "y" and "z".
#' @param env An environment where the symbols and functions are searched for.
#'  Defaults to \code{parent.frame()} for \code{f} expression and to
#'  \code{environment(f)} if \code{f} is a function.
#' @return \itemize{
#'  \item a function if \code{f} is a function
#'  \item an expression if \code{f} is an expression
#'  \item a character string if \code{f} is a character string
#'  \item a language (usually a so called 'call' but may be also a symbol or just a numeric) for other types of \code{f}
#' }
#'
#' @details
#' R already contains two differentiation functions: D and deriv.  D does
#' simple univariate differentiation.  "deriv" uses D do to multivariate
#' differentiation.  The output of "D" is an expression, whereas the output of
#' "deriv" is an executable function.
#
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
#	- if the function is in the derivative table, then the chain rule
#	is applied.  For example, if you declared that the derivative of
#	sin is cos, then it would figure out how to call cos correctly.
#	- if the function is not in the derivative table (or it is anonymous),
#	then the function body is substituted in.
#	- these two methods can be mixed.  An entry in the derivative table
#	need not be self-contained -- you don't need to provide an infinite
#	chain of derivatives.
#'  \item It's easy to add custom entries to the derivatives table, e.g.
#	drule[["cos"]] <- list(-._d1*sin(._1))
#'  \item The output is an executable function, which makes it suitable
#' for use in optimization problems.
#' }
#' 
#' Two working environments \code{drule} and \code{simplifications} are
#' exported in the package namescape.
#' As their names indicates, they contain tables of derivative and
#' simplification rules.
#' To see the list of defined rules do \code{ls(drule)}.
#' To add your own derivative rule for a function called say \code{sinpi(x)} calculating sin(pi*x), do \code{drule[["sinpi"]] <- list(quote(pi*._d1*cospi(._1)))}.
#' Here, "._1" stands for the "first arguments", "._d1" for the first derivative of the first arguments. For a function that might have more than one argument,
#' e.g. log(x, base=exp(1)), the drule entry must be a list with one rule
#' per number of possible arguments. See \code{drule$log} for an example to follow.
#' After adding \code{sinpi} you can derivate expressions like \code{Deriv(~ sinpi(x^2), "x")}
#' 
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
#' # -(2 * (sin(x) * y * sin(y^2)))
#' 
#' \dontrun{Deriv(quote(f(x, y^2)), c("x", "y"))}
#' # c(x = cos(x) * cos(y^2), y = -(2 * (sin(x) * y * sin(y^2))))
#' 
#' \dontrun{Deriv(expression(sin(x^2) * y), "x")}
#' # expression(cos(x^2) * (2 * x) * y)
#' 
#' Deriv("sin(x^2) * y", "x")
#' "2 * (x * cos(x^2) * y)"

# This wrapper of Deriv_ returns an appropriate expression (wrapped up "properly") depending on the type of argument to differentiate
Deriv <- function(f, x=if (is.function(f)) names(formals(f)) else NA, env=if (is.function(f)) environment(f) else parent.frame()) {
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
		return(as.function(c(af, Deriv_(as.call(c(as.symbol(fch), lapply(x, as.symbol))), x, env)), envir=env))
	} else if (is.na(x[1])) {
		stop("First argument is not a function, so variable name(s) must be supplied in the second argument")
	}
	x <- as.character(x)
	if (any(nchar(x) == 0)) {
		stop("Names in the second argument must not be empty")
	}
	if (is.character(f)) {
		# f is to parse
		format1(Deriv_(parse(text=f)[[1]], x, env))
	} else if (is.function(f)) {
		if (is.primitive(f)) {
			fch <- deparse(substitute(f))
		} else {
			as.function(c(formals(f), Deriv_(body(f), x, env)), envir=env)
		}
	} else if (is.expression(f)) {
		as.expression(Deriv_(f[[1]], x, env))
	} else if (is.language(f)) {
		if (is.call(f) && f[[1]] == as.symbol("~")) {
			# rhs of the formula
			Deriv_(f[[length(f)]], x, env)
		} else {
			# plain call derivation
			Deriv_(f, x, env)
		}
	} else {
		stop("Invalid type of 'f' for differentiation")
	}
}

# workhorse function doing the main work of differentiation
Deriv_ <- function(st, x, env) {
	# Make x scalar and wrap results in a c() call if length(x) > 1
	if (length(x) > 1) {
		# many variables => recursive call on single name
		res <- lapply(x, function(xi) Deriv_(st, xi, env))
		names(res) <- x;
		return(as.call(c(as.symbol("c"), res)))
	}
	# differentiate R statement 'st' (a call, or a symbol or numeric) by a name in 'x'
	if (is.unumeric(st)) {
		return(0)
	} else if (is.symbol(st)) {
		if (as.character(st) == x) {
			return(1)
		} else {
			return(0)
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
		if (stch %in% dlin) {
			# linear case
			# differentiate all arguments then pass them to the function
			dargs <- lapply(st[-1], function(a) Simplify_(Deriv_(a, x, env)))
			return(as.call(c(st[[1]], dargs)))
		}
		nb_args=length(st)-1
		if (is.null(drule[[stch]])) {
			# no derivative rule for this function
			# try to get the body and differentiate it
			ff <- get(stch, mode="function", envir=env)
			bf <- body(ff)
			if (is.null(bf)) {
				stop(sprintf("Could not retrieve body of '%s()'", stch))
			}
			args <- as.list(st[-1])
			mc <- match.call(ff, st)
			st <- Simplify_(eval(call("substitute", bf, as.list(mc[-1]))))
			return(Deriv_(st, x, env))
		} else if (is.null(drule[[stch]][[nb_args]])) {
			stop(sprintf("Don't know how to differentiate function or operator '%s' when it is called with %s arguments", stch, nb_args))
		}
		# there is a rule!
		# prepare replacement list ._1 -> first argument, ._d1 -> derivative of the first argument and so on
		args <- as.list(st[-1])
		names(args) <- paste("._", seq_len(nb_args), sep="")
		# which arguments have to be differentiated?
		dstr <- format1(drule[[stch]][[nb_args]])
		ma <- gregexpr("\\._d[1-9][0-9]*", dstr)
		dgrep <- substring(dstr, ma[[1]], ma[[1]]+attr(ma[[1]], "match.length")-1)
		dnum <- as.integer(sub("._d", "", unique(dgrep)))
		dargs <- lapply(args[dnum], Deriv_, x, env)
#cat("st=", format1(st), "\n", sep="")
#cat("stch=", stch, "\n", sep="")
#cat("dgrep=", dgrep, "\n", sep=", ")
#cat("dstr=", dstr, "\n", sep="")
#cat("dnum=", dnum, "\n", sep=", ")
		names(dargs) <- paste("._d", dnum, sep="")
		return(Simplify_(eval(call("substitute", drule[[stch]][[nb_args]], c(args, dargs)))))
	} else if (is.function(st)) {
		# differentiate its body if can get it
		args <- st[-1]
		names(args)=names(formals(ff))
		if (is.null(names(args))) {
			stop(sprintf("Could not retrieve arguments of '%s()'", stch))
		}
		st <- eval(call("substitute", body(ff), args))
		Deriv_(st, x, env)
	} else {
		stop("Invalid type of 'st' argument. It must be numeric, symbol or a call.")
	}
}

drule <- new.env()

# linear functions, i.e. d(f(x))/dx == f(d(arg)/dx)
dlin=c("-", "c", "t", "sum")

# first item in the list correspond to a call with one argument
# second (if any) for two, third for three. NULL means that with
# this number of argument this function can not be called
drule[["("]] <- list(quote(._d1)) # (._1) => omit paranthesis
# linear arithmetics are already handled by dlin
# exception is maid for unitary plus (it is just omitted)
drule[["+"]] <- list(quote(._d1), quote(._d1+._d2)) # +._1, ._1+._2
#drule[["-"]] <- list(quote(-._d1), quote(._d1-._d2)) # -._1, ._1-._2
# arithmetic non linear rules
drule[["*"]] <- list(NULL, quote(._d1*._2+._1*._d2)) # ._1*._2
drule[["/"]] <- list(NULL, quote(._d1/._2-._1*._d2/._2^2)) # ._1*._2
# power functions
drule[["^"]] <- list(NULL, quote(._d1*._2*._1^(._2-1)+._d2*log(._1)*._1^._2)) # ._1^._2
# example of recursive call
#drule[["sqrt"]] <- list(Deriv_(call("^", as.symbol("._1"), 0.5), "._1", NULL)) # sqrt(._1)
# but we prefer a sqrt() formula
drule[["sqrt"]] <- list(quote(0.5*._d1/sqrt(._1)))
drule[["log"]] <- list(quote(._d1/._1), quote(._d1/(._1*log(._2))- ._d2*log(._2, ._1)/(._2*log(._2)))) # log(._1), log(._1, b)
drule[["logb"]] <- drule[["log"]]
drule[["log2"]] <- list(quote(._d1/(._1*log(2))))
drule[["log10"]] <- list(quote(._d1/(._1*log(10))))
drule[["exp"]] <- list(quote(._d1*exp(._1)))
# trigonometric
drule[["sin"]] <- list(quote(._d1*cos(._1)))
drule[["cos"]] <- list(quote(-._d1*sin(._1)))
drule[["tan"]] <- list(quote(._d1/cos(._1)^2))
drule[["asin"]] <- list(quote(._d1/sqrt(1-._1^2)))
drule[["acos"]] <- list(quote(-._d1/sqrt(1-._1^2)))
drule[["atan"]] <- list(quote(._d1/(1+._1^2)))
# hyperbolic
drule[["sinh"]] <- list(quote(._d1*cosh(._1)))
drule[["cosh"]] <- list(quote(._d1*sinh(._1)))
drule[["tanh"]] <- list(quote(._d1*(1-tanh(._1)^2)))
drule[["asinh"]] <- list(quote(._d1/sqrt(._1^2+1)))
drule[["acosh"]] <- list(quote(._d1/sqrt(._1^2-1)))
drule[["atanh"]] <- list(quote(._d1/(1-._1^2)))
# code control
drule[["if"]] <- list(NULL, quote(if (._1) ._d2), quote(if (._1) ._d2 else ._d3))
