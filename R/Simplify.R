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
#browser()
		# get numerator and denumerator for a and b than combine them
		nd_a=Numden(a)
		nd_b=Numden(b)
		nd=list(num=c(nd_a$num, nd_b$num), den=c(nd_a$den, nd_b$den))
		# reduce numerics to only one factor
		fa=list()
		for (na in c("num", "den")) {
			inu=which(sapply(nd[[na]], is.numeric))
			if (length(inu)) {
				fa[[na]]=prod(unlist(nd[[na]][inu]))
				# remove numerics
				nd[[na]]=nd[[na]][-inu]
			}
		}
		fa$num=if (length(fa$num)) fa$num else 1
		# simplify identical terms in num and denum
		nd_eq=outer(sapply(nd$den, format1), sapply(nd$num, format1), `==`)
		ipair=matrix(0, nrow=2, ncol=0)
		for (inum in seq(len=ncol(nd_eq))) {
			iden=which(nd_eq[,inum])
			iden=iden[!iden %in% ipair[2,]]
			if (length(iden)) {
				# remove this pair
				ipair=cbind(ipair, c(inum, iden[1]))
			}
		}
		if (ncol(ipair) > 0) {
			nd$num=nd$num[-ipair[1,]]
			nd$den=nd$den[-ipair[2,]]
		}
		# form symbolic products
		eprod=list()
		for (na in c("num", "den")) {
			if (length(nd[[na]]) == 0) next
			eprod[[na]]=nd[[na]][[1]]
			for (term in nd[[na]][-1]) {
				eprod[[na]]=call("*", eprod[[na]], term)
			}
			if (fa[[na]] != 1) {
				eprod[[na]]=if (is.null(eprod[[na]])) fa[[na]] else call("*", fa[[na]], eprod[[na]])
			}
		}
		eprod$num=if (is.null(eprod$num)) 1 else eprod$num
		if (is.null(eprod$den)) {
			# we have no denominator
			eprod$num
		} else {
			call("/", eprod$num, eprod$den)
		}
	}
}
`Simplify./` <- function(expr)
{
	a <- Simplify_(expr[[2]])
	b <- Simplify_(expr[[3]])

	if (is.numeric(a) && all(a == 0)) {
		0
	} else if (is.numeric(b) && all(b == 0)) {
		Inf
	} else if (is.numeric(b) && all(b == 1)) {
		a
	} else if (is.numeric(a) && is.numeric(b)) {
		a/b
	} else if (a==1) {
		expr
	} else {
		`Simplify.*`(call("*", a, substitute(1/b)))
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
`Simplify.reciprocal` <- function(expr)
{
	a <- Simplify_(expr[[2]])
	if (is.numeric(a)) {
		1/a
	} else {
		substitute(1/a)
	}
}
`Simplify.neg.sin` <- function(expr)
{
	a <- Simplify_(expr[[2]])
	if (is.numeric(a)) {
		-sin(a)
	} else {
		substitute(-sin(a))
	}
}
`Simplify.reciprocal.deriv` <- function(expr)
{
	a <- Simplify_(expr[[2]])
	if (is.numeric(a)) {
		-1/(a*a)
	} else {
		substitute(-1/(a*a))
	}
}
`Simplify.sqrt.deriv` <- function(expr)
{
	a <- Simplify_(expr[[2]])
	if (is.numeric(a)) {
		0.5/sqrt(a)
	} else {
		substitute(0.5/sqrt(a))
	}
}
`Simplify.tan.deriv` <- function(expr)
{
	a <- Simplify_(expr[[2]])
	if (is.numeric(a)) {
		1/cos(a)**2
	} else {
		substitute(1/cos(a)**2)
	}
}
Numden <- function(expr) {
	# Return a list with "num" as numerator and and "den" as denominator sublists.
	# Each sublist regroups the language expressions which are not products neither
	# divisions
	if (is.numeric(expr) || is.symbol(expr)) {
		list(num=list(expr), den=list(1))
	} else if (expr[[1]] == as.symbol("*")) {
		# recursive call
		a=Numden(expr[[2]])
		b=Numden(expr[[3]])
		list(num=c(a$num, b$num), den=c(a$den, b$den))
	} else if (expr[[1]] == as.symbol("/")) {
		# recursive call
		a=Numden(expr[[2]])
		b=Numden(expr[[3]])
		list(num=c(a$num, b$den), den=c(a$den, b$num))
	} else {
		list(num=list(expr), den=list(1))
	}
}
format1 <- function(expr) {
	res <- if (is.symbol(expr)) as.character(expr) else format(expr)
	if (length(res) > 1) {
		res=paste(res, collapse="")
	}
	return(res)
}
