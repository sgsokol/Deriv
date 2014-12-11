#' @name Deriv
#' @title Symbollic differentiation of an expression or function
#' @aliases Deriv.function Deriv derivatives simplifications
#' @concept symbollic derivation
# \usage{
# Deriv.function(f, x = names(formals(f)), env = environment(f))
# }
#' 
#' 
#' @param f An expression or function to be differentiated
#' @param x A character string with variable name with resptect to which
#'  \code{f} must be differentiated. If \code{f} is function \code{x} is
#'  optional and defaults to \code{names(formals(f))}
#' @param env An environment where the symbols and functions are searched for.
#'  Defaults to \code{parent.frame()} for \code{f} expression and to
#'  \code{environment(f)} if \code{f} is a function.
#' @return An expression (for \code{Deriv}) or a function (for
#'  \code{Deriv.function}) with the derivative of \code{f}.
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
#' \dontrun{Deriv.function(f)}
#' # function (x) 
#' # 2 * x
#'
#' \dontrun{f <- function(x, y) sin(x) * cos(y)}
#' \dontrun{Deriv.function(f)}
#' # function (x, y) 
#' # sin(x) * (sin(y) * c(0, 1)) + cos(y) * (cos(x) * c(1, 0))
#'
#' \dontrun{f_ <- Deriv.function(f)}
#' \dontrun{f_(3, 4)}
#' #              x         y
#' # [1,] 0.6471023 0.1068000
#
#' \dontrun{Deriv(expression(f(x, y^2)), "y")}
#' # expression(sin(x) * (neg.sin(y^2) * (2 * y)))
#'
#' \dontrun{Deriv(expression(f(x, y^2)), c("x", "y"))}
#' # expression(sin(x) * (neg.sin(y^2) * (2 * (y * t(c(0, 1))))) + 
#' #     cos(x) * t(c(1, 0)) * cos(y^2))
#'
#' \dontrun{Deriv(expression(sin(x^2) * y), "x")}
#' # expression(cos(x^2) * (2 * x) * y)

# This wrapper of Deriv_ returns an expression (wrapped up "properly")
Deriv <- function(f, x, env=parent.frame())
	as.expression(Deriv_(f[[1]], x, env))

# This wrapper of Deriv_ takes a function rather than an expression.
# By default, it automatically figures out which variables to differentiate
# with respect to (all of them).  It returns an executable function.

Deriv.function <- function(f, x=names(formals(f)), env=environment(f))
{
	stopifnot(is.function(f))
	as.function(c(as.list(formals(f)),
			      Deriv_(body(f), x, env)),
		    envir=env)
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
			for (arg in 1:length(args))
				args[[arg]] <- f[[1 + arg]]
			subst.fun <- eval(call("substitute", body(sym), args))
			Simplify_(Deriv_(subst.fun, x, env))
		}
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

# FIXME: only works if the exponent is a constant.
`Deriv.^` <- function(f, x, env)
{
	stopifnot(is.numeric(f[[3]]))
	b <- f[[3]]
	if (b == 0) {
		0
	} else {
		`expr.*`(
			b,
			`expr.*`(`expr.^`(f[[2]], b - 1),
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

neg.sin <- function(x) -sin(x)

Deriv.ifelse <- function(f, x, env)
{
	f[[3]] <- Deriv_(f[[3]], x, env)
	f[[4]] <- Deriv_(f[[4]], x, env)
	f
}

.onLoad <- function(libname, pkgname) {
   assign("simplifications", new.env(), envir=environment(Deriv))

   assign("+", `Simplify.+`, envir=simplifications)
   assign("-", `Simplify.-`, envir=simplifications)
   assign("*", `Simplify.*`, envir=simplifications)
   assign("(", `Simplify.(`, envir=simplifications)
   assign("c", `Simplify.c`, envir=simplifications)
   assign("^", `Simplify.^`, envir=simplifications)

   assign("derivatives", new.env(), envir=environment(Deriv))

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
   assign("exp", chain.rule("exp"), envir=derivatives)
   assign("log", chain.rule("reciprocal"), envir=derivatives)
   assign("length", Deriv.constant, envir=derivatives)
   assign("sum", Deriv.sum, envir=derivatives)
   assign("sqrt", chain.rule("sqrt.deriv"), envir=derivatives)
   assign("ifelse", Deriv.ifelse, envir=derivatives)

   dnorm.slow <- function(x, mean, sd)
	         1/(sqrt(2*pi)*sd) * exp(-(x - mean)^2/(2*sd^2))
   dnorm.deriv <- Deriv.function(dnorm.slow)
   assign("dnorm", chain.rule("dnorm.deriv", TRUE), envir=derivatives)
}
