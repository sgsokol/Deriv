# Simplify.R -- symbollic simplification
# written by Andrew Clausen <clausen@econ.upenn.edu> in 2007
# thanks to a bug fix from Mark Reid <mark.reid@anu.edu.au> in 21/2/2009
#
# This isn't a serious attempt at simplification code.  It just does some
# obvious things like 0 + x => x.  It was written to support Deriv.R.

Simplify_ <- function(expr)
{
	if (is.symbol(expr)) {
		expr
	} else if (is.language(expr) && is.symbol(expr[[1]])) {
		# is there a rule in the table?
		sym.name <- as.character(expr[[1]])
		if (class(try(Simplify.rule <-
				get(sym.name, envir=simplifications,
				    inherits=FALSE), silent=TRUE))
				!= "try-error")
			return(Simplify.rule(expr))
	}
	expr
}

Simplify <- function(expr)
	as.expression(Simplify_(expr[[1]]))

Simplify.function <- function(f, x=names(formals(f)), env=parent.frame())
{
	stopifnot(is.function(f))
	as.function(c(as.list(formals(f)),
			      Simplify_(body(f))),
		    envir=env)
}

`Simplify.+` <- function(expr)
{
	a <- Simplify_(expr[[2]])
	b <- Simplify_(expr[[3]])

	if (is.numeric(a) && all(a == 0)) {
		b
	} else if (is.numeric(b) && all(b == 0)) {
		a
	} else if (is.numeric(a) && is.numeric(b)) {
		a + b
	} else {
		expr[[2]] <- a
		expr[[3]] <- b
		expr
	}
}

`Simplify.-` <- function(expr)
{
	if (length(expr) == 2)
	{
		if (is.numeric(expr[[2]]))
			return(-expr[[2]])
		return(expr)
	}

	a <- Simplify_(expr[[2]])
	b <- Simplify_(expr[[3]])

	if (is.numeric(a) && all(a == 0)) {
		if (is.numeric(b)) -b else substitute(-b)
	} else if (is.numeric(b) && all(b == 0)) {
		a
	} else if (is.numeric(a) && is.numeric(b)) {
		a - b
	} else {
		expr[[2]] <- a
		expr[[3]] <- b
		expr
	}
}

`Simplify.(` <- function(expr)
	expr[[2]]

`Simplify.*` <- function(expr)
{
	a <- Simplify_(expr[[2]])
	b <- Simplify_(expr[[3]])

	if (is.numeric(a) && all(a == 0)) {
		0
	} else if (is.numeric(b) && all(b == 0)) {
		0
	} else if (is.numeric(a) && all(a == 1)) {
		b
	} else if (is.numeric(b) && all(b == 1)) {
		a
	} else if (is.numeric(a) && is.numeric(b)) {
		a * b
	} else {
		expr[[2]] <- a
		expr[[3]] <- b
		expr
	}
}

`Simplify.^` <- function(expr)
{
	a <- Simplify_(expr[[2]])
	b <- Simplify_(expr[[3]])

	if (is.numeric(a) && all(a == 0)) {
		0
	} else if (is.numeric(b) && all(b == 0)) {
		1
	} else if (is.numeric(a) && all(a == 1)) {
		1
	} else if (is.numeric(b) && all(b == 1)) {
		a
	} else if (is.numeric(a) && is.numeric(b)) {
		a ^ b
	} else {
		expr[[2]] <- a
		expr[[3]] <- b
		expr
	}
}

`Simplify.c` <- function(expr)
{
	args <- expr[-1]
	args.simplified <- lapply(args, Simplify_)
	if (all(lapply(args.simplified, is.numeric))) {
		as.numeric(args.simplified)
	} else {
		for (i in 1:length(args))
			expr[[i + 1]] <- args.simplified[[i]]
		expr
	}
}
