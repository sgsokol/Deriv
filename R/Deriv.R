#' @name Deriv
#' @title Symbollic differentiation of an expression or function
#' @aliases Deriv derivatives simplifications
#' @concept symbollic derivation
# \usage{
# Deriv(f, x = names(formals(f)), env = environment(f))
# }
#' 
#' 
#' @param f An expression or function to be differentiated.
#'  f can be \itemize{
#'   \item a user defined function \code{function(x) x**n}
#'   \item a string \code{"x**n"}
#'   \item an expression \code{expression(x**n)}
#'   \item a call \code{call("^", quote(x), quote(n))}
#'   \item a language \code{quote(x**n)}
#'   \item a right hand side of a formula \code{~ x**n} or \code{y ~ x**n}
#'  }
#' @param x A character string with variable name(s) with resptect to which
#'  \code{f} must be differentiated. If \code{f} is a function \code{x} is
#'  optional and defaults to \code{names(formals(f))}
#' @param env An environment where the symbols and functions are searched for.
#'  Defaults to \code{parent.frame()} for \code{f} expression and to
#'  \code{environment(f)} if \code{f} is a function.
#' @return \itemize{
#'  \item a function if \code{f} is a function
#'  \item an expression if \code{f} is an expression
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
#'  \item The output of "deriv" can not be differentiated again.
#'  \item Neither function can substitute function calls.  eg:
#'	f <- function(x, y) x + y; deriv(f(x, x^2), "x")
#' }
#'
#' So, here are the advantages and disadvantages of this implementation:
#
#' GOOD POINTS:
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
#'  \item It's easy to add custom entries to the derivatives table.  It could be
#' easier though... it would be nice if something like
#	add.deriv("cos", "-sin(expr)")
#' worked, rather than the clunky function definitions.  (This is purely a
#' cosmetic issue, though... everything works as is.)
#'  \item The output is an executable function, which makes it suitable for use in
#' optimization problems.
#' }
#'
#' BAD POINTS:
#' \itemize{
#'  \item Differentiating vector-valued functions doesn't work properly, since
#' the multiplication code doesn't know when to use scalar vs matrix
#' multiplication.  Unfortunately, solving this is a hard problem because
#' we would need to know if an arbitrary expression is a vector or not.
#' We would have to add extra metadata to do this.  Bottom line: can compute
#' gradients but not Jacobians or Hessians.
#'  \item Gives useless error messages when it gets stuck.  This could be fixed.
#' }
#' Two working environments derivatives and simplifications are created in the package
#' namescape. As their names indicates, they contain tables of derivatives and
#' simplification rules. A priori, user does not had to manipulate them directly.
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
#' # sin(x) * -sin(y) * t(c(0, 1)) + cos(x) * t(c(1, 0)) *  cos(y)
#'
#' \dontrun{f_ <- Deriv(f)}
#' \dontrun{f_(3, 4)}
#' #              x         y
#' # [1,] 0.6471023 0.1068
#
#' \dontrun{Deriv(~ f(x, y^2)), "y")}
#' # 2 * (sin(x) * -sin(y^2) * y)
#'
#' \dontrun{Deriv(quote(f(x, y^2)), c("x", "y"))}
#' # 2 * (sin(x) * -sin(y^2) * y * t(c(0, 1))) + cos(x) * 
#' #     t(c(1, 0)) * cos(y^2))
#'
#' \dontrun{Deriv(expression(sin(x^2) * y), "x")}
#' # expression(cos(x^2) * (2 * x) * y)

# This wrapper of Deriv_ returns an expression (wrapped up "properly")
# or call or function depending on type of 'f'
Deriv <- function(f, x=if (length(find(deparse(substitute(f)))) && is.function(f)) names(formals(f)) else stop("Argument 'f' is not a function, so variable name(s) must be supplied in 'x' argument"), env=if (is.function(f)) environment(f) else parent.frame()) {
	x # referense x here for a possible error message
	if (is.character(f)) {
		# f is to parse
		format1(mderst(parse(t=f)[[1]], x, env))
	} else if (is.function(f)) {
		as.function(c(formals(f), mderst(body(f), x, env)), envir=env)
	} else if (is.expression(f)) {
		as.expression(mderst(f[[1]], x, env))
	} else if (is.language(f)) {
		if (is.call(f) && f[[1]] == as.symbol("~")) {
			# rhs of the formula
			mderst(f[[length(f)]], x, env)
		} else {
			# plain call derivation
			mderst(f, x, env)
		}
	} else {
		stop("Invalid type of 'f' for differentiation")
	}
}

# Wrapper for derst(). Can get unknown functions (substitute body) and
# differentiate by many variables. In the latter case a named vector is returned
# (i.e. the gradient vector expression)
mderst <- function(f, x, env) {
	if (is.numeric(f)) {
		result <- numeric(length(x))
		if (length(x) > 1) {
			names(result) <- x
		}
		return(result)
	}
	if (length(x) > 1) {
		# many variables
		res <- lapply(x, function(xi) derst(f, xi, env))
		return(as.call(c(as.symbol("c"), res)))
	} else {
		# only one variable to differentiate by
		return(derst(f, x, env))
	}
}

Deriv.ifelse <- function(f, x, env)
{
	f[[3]] <- Deriv_(f[[3]], x, env)
	f[[4]] <- Deriv_(f[[4]], x, env)
	f
}

derst <- function(st, x, env) {
#browser()
	# differentiate R statement 'st' (a call, or a symbol or numeric) by a name in 'x'
	if (is.numeric(st) || ((is.uminus(st) || is.uplus(st)) && is.numeric(st[[2]]))) {
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
			dargs <- lapply(st[-1], function(a) Simplify_(derst(a, x, env)))
			return(as.call(c(st[[1]], dargs)))
		}
		nb_args=length(st)-1
		if (is.null(drule[[stch]])) {
			# no derivative rule for this function
			# try to get the body and differentiate it
			ff <- get(stch, mode="function", envir=env)
			args <- st[-1]
			names(args) <- names(formals(ff))
			if (is.null(names(args))) {
				stop(sprintf("Could not retrieve arguments of '%s()'", stch))
			}
		} else if (is.null(drule[[stch]][[nb_args]])) {
			stop(sprintf("Don't know how to differentiate function or operator '%s' when it is called with %s arguments", stch, nb_args))
		}
		if (nb_args <= 3) {
			lrepl=list(
				a=st[[2]],
				d_a=derst(st[[2]], x, env)
			)
		}
		if (nb_args > 1 && nb_args <= 3) {
			lrepl$b <- st[[3]]
			lrepl$d_b <- derst(st[[3]], x, env)
		}
		if (nb_args > 2 && nb_args <= 3) {
			lrepl$c <- st[[4]]
			lrepl$d_c <- derst(st[[4]], x, env)
		}
		
		return(Simplify_(eval(call("substitute", drule[[stch]][[nb_args]], lrepl))))
	} else if (is.function(st)) {
		# differentiate its body if can get it
		args <- st[-1]
		names(args)=names(formals(ff))
		if (is.null(names(args))) {
			stop(sprintf("Could not retrieve arguments of '%s()'", stch))
		}
		st <- eval(call("substitute", body(ff), args))
	} else {
		stop("Invalid type of 'st' argument. It must be numeric, symbol or a call.")
	}
}

.onLoad <- function(libname, pkgname) {
   assign("simplifications", new.env(), envir=environment(Deriv))
   assign("derivatives", new.env(), envir=environment(Deriv))
   assign("drule", new.env(), envir=environment(Deriv))
   
   # linear functions, i.e. d(f(x))/dx == f(d(arg)/dx)
   dlin=c("-", "c", "t", "sum")
   assign("dlin", dlin, envir=environment(Deriv))
   
   # first item in the list correspond to a call with one argument
   # second (if any) for two, third for three. NULL means that with
   # this number of argument this function can not be called
   drule[["("]] <- list(quote(d_a)) # (a) => omit paranthesis
   # linear arithmetics are already handled by dlin
   # exception is maid for unitary plus (it is just omitted)
   drule[["+"]] <- list(quote(d_a), quote(d_a+d_b)) # +a, a+b
   #drule[["-"]] <- list(quote(-d_a), quote(d_a-d_b)) # -a, a-b
   # arithmetic non linear rules
   drule[["*"]] <- list(NULL, quote(d_a*b+a*d_b)) # a*b
   drule[["/"]] <- list(NULL, quote(d_a/b-a*d_b/b^2)) # a*b
   # power functions
   drule[["^"]] <- list(NULL, quote(d_a*b*a^(b-1)+d_b*log(a)*a^b)) # a^b
   # example of recursive call
   #drule[["sqrt"]] <- list(derst(call("^", as.symbol("a"), 0.5), "a", NULL)) # sqrt(a)
   # but we prefer a sqrt() formula
   drule[["sqrt"]] <- list(quote(0.5*d_a/sqrt(a)))
   drule[["log"]] <- list(quote(d_a/a), quote(d_a/(a*log(b)))) # log(a), log(a, b)
   drule[["logb"]] <- drule[["log"]]
   drule[["log2"]] <- list(quote(d_a/(a*log(2))))
   drule[["log10"]] <- list(quote(d_a/(a*log(10))))
   drule[["exp"]] <- list(quote(d_a*exp(a)))
   # trigonometric
   drule[["sin"]] <- list(quote(d_a*cos(a)))
   drule[["cos"]] <- list(quote(-d_a*sin(a)))
   drule[["tan"]] <- list(quote(d_a/cos(a)^2))
   drule[["asin"]] <- list(quote(d_a/sqrt(1-a^2)))
   drule[["acos"]] <- list(quote(-d_a/sqrt(1-a^2)))
   drule[["atan"]] <- list(quote(d_a/(1+a^2)))
   # hyperbolic
   drule[["sinh"]] <- list(quote(d_a*cosh(a)))
   drule[["cosh"]] <- list(quote(d_a*sinh(a)))
   drule[["tanh"]] <- list(quote(d_a*(1-tanh(a)^2)))
   drule[["asinh"]] <- list(quote(d_a/sqrt(a^2+1)))
   drule[["acosh"]] <- list(quote(d_a/sqrt(a^2-1)))
   drule[["atanh"]] <- list(quote(d_a/(1-a^2)))
   # stats
   drule[["dnorm"]] <- list(NULL, quote(d_a) )

   assign("+", `Simplify.+`, envir=simplifications)
   assign("-", `Simplify.-`, envir=simplifications)
   assign("*", `Simplify.*`, envir=simplifications)
   assign("/", `Simplify./`, envir=simplifications)
   assign("(", `Simplify.(`, envir=simplifications)
   assign("^", `Simplify.^`, envir=simplifications)

   #assign("ifelse", Deriv.ifelse, envir=derivatives)
}
