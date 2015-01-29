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
		Deriv(parse(t=f), x, env)
	} else if (is.function(f)) {
		as.function(c(as.list(formals(f)), Deriv_(body(f), x, env)),  envir=env)
	} else if (is.expression(f)) {
		as.expression(Deriv_(f[[1]], x, env))
	} else if (is.language(f)) {
		if (is.call(f) && f[[1]] == as.symbol("~")) {
			# rhs of a formula
			Deriv_(f[[length(f)]], x, env)
		} else {
			# plain call derivation
			Deriv_(f, x, env)
		}
	} else {
		stop("Unrecognized type of 'f' for derivation")
	}
}

# This is the main function.  It takes an expression (with the first layer
# stripped -- I don't know the terminology) f as an argument, and
# differentiates with respect to the variables listed in x.  The environment
# parameter, env specifies where to search for functions invoked in the
# expression

Deriv_ <- function(f, x, env)
{
	if (is.numeric(f)) {
		result <- numeric(length(x))
		names(result) <- x
		result
	} else if (is.symbol(f)) {
		# if f is the x in "2x + 3", and we're differentiating with
		# respect to x only, then the derivative of the "x" bit is
		# just "1".  If we're differentiating wrt "y", then it is 0.
		# If we're differentiating with respect to both "x" and "y",
		# then it is t(c(1, 0)).
		vector <- mapply(function(var) 1 * (as.character(f) == var), x)
		if (length(x) == 1) {
			vector
		} else {
			result <- expression(t(x))[[1]]
			result[[2]] <- vector
			result
		}
	} else if (is.language(f) && is.symbol(f[[1]]) ) {
		# is there a rule in the table?
		sym.name <- as.character(f[[1]])
		if (class(try(Deriv.rule <- get(sym.name, envir=derivatives,
						inherits=FALSE), silent=TRUE))
				!= "try-error") {
			# There is a rule... apply it.
			Simplify_(Deriv.rule(f, x, env))
		} else {
			# There is no rule... substitute the function in,
			# and differentiate the result.
			sym <- get(sym.name, envir=env)
			stopifnot(is.function(sym))
			args <- formals(sym)
			if (is.null(args)) {
				stop(sprintf("Could not retrieve arguments of '%s()'", sym.name))
			}
			for (arg in 1:length(args))
				args[[arg]] <- f[[1 + arg]]
			subst.fun <- eval(call("substitute", body(sym), args))
			Simplify_(Deriv_(subst.fun, x, env))
		}
	}
}

# Wrapper for derst(). Can get unknown functions (substitute body) and
# differentiate by many variables. In the latter case a named vector is returned
# (i.e. the gradient vector expression)
mderst <- function(f, x, env) {
	if (is.numeric(f)) {
		result <- numeric(length(x))
		names(result) <- x
		return(result)
	}
	if (length(x) > 1) {
		# many variables
		res <- lapply(x, function(xi) derst(st, xi))
		return(as.call(c(as.symbol("c"), res)))
	} else {
		# only one variable to differentiate by
		return(derst(st, xi))
	}
}

subop <- function(expr) function(a, b)
{
	expr[[1]][[2]] <- a
	expr[[1]][[3]] <- b
	expr[[1]]
}

`expr.+` <- subop(expression(a + b))
`expr.-` <- subop(expression(a - b))
`expr.*` <- subop(expression(a * b))
`expr.%*%` <- subop(expression(a %*% b))
`expr.^` <- subop(expression(a ^ b))

Deriv.constant <- function(f, x, env)
	0

`Deriv.+` <- function(f, x, env)
	`expr.+`(Deriv_(f[[2]], x, env), Deriv_(f[[3]], x, env))

`Deriv.-` <- function(f, x, env)
{
	if(length(f) == 3) {
		`expr.-`(Deriv_(f[[2]], x, env), Deriv_(f[[3]], x, env))
	} else {
		result <- expression(-a)[[1]]
		result[[2]] <- Deriv_(f[[2]], x, env)
		result
	}
}

# FIXME: don't know when both arguments (or their derivatives) are vectors
`Deriv.*` <- function(f, x, env)
	`expr.+`(
		`expr.*`(f[[2]], Deriv_(f[[3]], x, env)),
		`expr.*`(Deriv_(f[[2]], x, env), f[[3]]))

`Deriv.%*%` <- function(f, x, env)
	`expr.+`(
		`expr.%*%`(f[[2]], Deriv_(f[[3]], x, env)),
		`expr.%*%`(Deriv_(f[[2]], x, env), f[[3]]))


`Deriv./` <- function(f, x, env)
{
	expr <- expression(u * reciprocal(v))[[1]]
	expr[[2]] <- f[[2]]
	expr[[3]][[2]] <- f[[3]]
	`Deriv.*`(expr, x, env)
}

reciprocal <- function(x)
{
	if (is.matrix(x)) {
		solve(x)
	} else {
		1/x
	}
}

reciprocal.deriv <- function(x) -1/x^2

`Deriv.(` <- function(f, x, env)
	Deriv_(f[[2]], x, env)

`Deriv.^` <- function(f, x, env)
{
	b <- f[[3]]
	if (is.numeric(b) && b == 0) {
		0
	} else {
		b=Simplify_(b)
		`expr.*`(
			b,
			`expr.*`(`expr.^`(f[[2]], substitute(b - 1)),
				 Deriv_(f[[2]], x, env)))
	}
}

do.rbind <- function(rows)
{
	n <- length(rows)
	result <- expression(rbind(rows))[[1]]
	for (i in 1:n)
		result[i+1] <- rows[i]
	result
}

do.t <- function(cols)
{
	n <- length(cols)
	result <- expression(t(cols))[[1]]
	for (i in 1:n)
		result[i+1] <- cols[i]
	result
}

Deriv.c <- function(f, x, env)
	do.rbind(lapply(f[-1], function (arg) Deriv_(arg, x, env)))

Deriv.t <- function(f, x, env)
	do.t(lapply(f[-1], function (arg) Deriv_(arg, x, env)))

# This function applies the chain rule.  For example,
#	chain.rule("cos")
# would be a valid function for differentiating "sin".
#
# f.deriv gives the derivative of the function.
# The chain rule evaluates f'(g(x)) * g'(x).
# The "u" term corresponds to f'(g(x)).
# The "v" term corresponds to g'(x).
chain.rule <- function(f.deriv, vec.valued=FALSE) function(expr, x, env)
{
	u <- expr
	u[[1]] <- as.symbol(f.deriv)

	n <- length(expr)
	v <- lapply(expr[2:n], function(arg) Deriv_(arg, x, env))

	if (vec.valued) {
		`expr.%*%`(u, do.rbind(v))
	} else {
		`expr.*`(u, v[[1]])
	}
}

Deriv.sum <- function(f, x, env)
{
	result <- expression(sum(terms))[[1]]
	result[[2]] <- Deriv_(f[[2]], x, env)
	result
}

sqrt.deriv <- function(x) 1/(2*sqrt(x))
tan.deriv <- function(x) 1/cos(x)**2

neg.sin <- function(x) -sin(x)

Deriv.ifelse <- function(f, x, env)
{
	f[[3]] <- Deriv_(f[[3]], x, env)
	f[[4]] <- Deriv_(f[[4]], x, env)
	f
}
dnorm.slow <- function(x, mean, sd)
	1/(sqrt(2*pi)*sd) * exp(-(x - mean)^2/(2*sd^2))
dnorm.deriv <- function(x, mean, sd) {
	f <- function(x, mean, sd){}
	body(f) <- Deriv_(body(dnorm.slow), "x")
	return(f)
}

repl_ab <- function(st, lrepl) {
	# In a statement st, replace symbols a, d_a etc by the content
	# the correponding item in the list lrepl
	if (is.numeric(st)) {
		return(st)
	} else if (is.symbol(st)) {
		stch <- as.character(st)
		repl <- lrepl[[stch]]
		return(if (is.null(repl)) st else repl)
	} else if (is.call(st)) {
		# recursive call
		for (i in 2:length(st)) {
			# replace a,b, in all arguments
			st[[i]] <- repl_ab(st[[i]], lrepl)
		}
		return(st)
	}
}
derst <- function(st, x, env) {
	# differentiate R statement 'st' (a call, or a symbol or numeric) by a name in 'x'
	if (is.numeric(st)) {
		return(0)
	} else if (is.symbol(st)) {
		if (as.character(st) == x) {
			return(1)
		} else {
			return(0)
		}
	} else if (is.call(st)) {
		stch <- as.character(st[[1]])
		# prepare expression to differentiate
		if (is.null(drule[[stch]])) {
			# no known rule for that function
			# differentiate its body if can get it
			ff <- get(stch, envir=env, mode="function")
			args <- st[-1]
			names(args)=names(formals(ff))
			if (is.null(names(args))) {
				stop(sprintf("Could not retrieve arguments of '%s()'", stch))
			}
			st <- eval(call("substitute", body(ff), args))
		}
		nb_args=length(st)-1
		if (is.null(drule[[stch]][[nb_args]])) {
			stop(sprintf("Don't know how to differentiate function or operator '%s' when it is called with %s arguments", stch, nb_args))
		}
		lrepl=list(
			a=st[[2]],
			d_a=derst(st[[2]], x, env)
		)
		if (nb_args == 2) {
			lrepl$b <- st[[3]]
			lrepl$d_b <- derst(st[[3]], x, env)
		}
		
		return(Simplify_(eval(call("substitute", drule[[stch]][[nb_args]], lrepl))))
	} else {
		stop("Invalid type of 'st' argument. It must be numeric, symbol or a call.")
	}
}

.onLoad <- function(libname, pkgname) {
   assign("simplifications", new.env(), envir=environment(Deriv))
   assign("derivatives", new.env(), envir=environment(Deriv))
   assign("drule", new.env(), envir=environment(Deriv))
   
   # arithmetic rules
   # first item in the list correspond to a call with one argument
   # second (if any) for two. NULL means that with this number of argument
   # a function can not be called
   drule[["+"]] <- list(quote(d_a), quote(d_a+d_b)) # +a, a+b
   drule[["-"]] <- list(quote(-d_a), quote(d_a-d_b)) # -a, a-b
   drule[["*"]] <- list(NULL, quote(d_a*b+a*d_b)) # a*b
   drule[["/"]] <- list(NULL, quote((d_a*b-a*d_b)/b^2)) # a*b
   # power functions
   drule[["^"]] <- list(NULL, quote(d_a*b*a^(b-1)+d_b*log(a)*a^b)) # a^b
   # example of recursive call
   #drule[["sqrt"]] <- list(derst(call("^", as.symbol("a"), 0.5), "a")) # sqrt(a)
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

   assign("+", `Simplify.+`, envir=simplifications)
   assign("-", `Simplify.-`, envir=simplifications)
   assign("*", `Simplify.*`, envir=simplifications)
   assign("/", `Simplify./`, envir=simplifications)
   assign("(", `Simplify.(`, envir=simplifications)
   assign("c", `Simplify.c`, envir=simplifications)
   assign("^", `Simplify.^`, envir=simplifications)
   assign("reciprocal", `Simplify.reciprocal`, envir=simplifications)
   assign("neg.sin", `Simplify.neg.sin`, envir=simplifications)
   assign("reciprocal.deriv", `Simplify.reciprocal.deriv`, envir=simplifications)
   assign("sqrt.deriv", `Simplify.sqrt.deriv`, envir=simplifications)
   assign("tan.deriv", `Simplify.tan.deriv`, envir=simplifications)

   assign("+", `Deriv.+`, envir=derivatives)
   assign("-", `Deriv.-`, envir=derivatives)
   assign("*", `Deriv.*`, envir=derivatives)
   assign("%*%", `Deriv.%*%`, envir=derivatives)
   assign("/", `Deriv./`, envir=derivatives)
   assign("reciprocal", chain.rule("reciprocal.deriv"), envir=derivatives)
   assign("^", `Deriv.^`, envir=derivatives)
   assign("(", `Deriv.(`, envir=derivatives)
   #assign("c", Deriv.c, envir=derivatives)
   #assign("t", Deriv.t, envir=derivatives)
   assign("sin", chain.rule("cos"), envir=derivatives)
   assign("cos", chain.rule("neg.sin"), envir=derivatives)
   assign("tan", chain.rule("tan.deriv"), envir=derivatives)
   assign("exp", chain.rule("exp"), envir=derivatives)
   assign("log", chain.rule("reciprocal"), envir=derivatives)
   assign("length", Deriv.constant, envir=derivatives)
   assign("sum", Deriv.sum, envir=derivatives)
   assign("sqrt", chain.rule("sqrt.deriv"), envir=derivatives)
   assign("ifelse", Deriv.ifelse, envir=derivatives)
   assign("dnorm", chain.rule("dnorm.deriv", TRUE), envir=derivatives)
}
