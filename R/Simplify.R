# Simplify.R -- symbollic simplification
# written by Andrew Clausen <clausen@econ.upenn.edu> in 2007
# thanks to a bug fix from Mark Reid <mark.reid@anu.edu.au> in 21/2/2009
#
# This isn't a serious attempt at simplification code.  It just does some
# obvious things like 0 + x => x.  It was written to support Deriv.R.

Simplify_ <- function(expr)
{
	if (is.unumeric(expr)) {
		eval(expr)
	} else if (is.call(expr)) {
		args <- lapply(expr[-1], Simplify_)
		if (all(sapply(args, is.numeric))) {
			# if all arguments are numeric, evaluate them
			return(eval(expr))
		} else {
			# is there a rule in the table?
			sym.name <- as.character(expr[[1]])
			if (class(try(Simplify.rule <-
					get(sym.name, envir=simplifications,
					inherits=FALSE), silent=TRUE))
					!= "try-error") {
				expr[-1]=args
				return(Simplify.rule(expr))
			} else {
				expr[-1]=args
				return(expr)
			}
		}
	} else {
		expr
	}
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

# in what follows no need to Simplify_ args neither to check if
# all arguments are unumeric. It is done in upper Simplify_()
`Simplify.(` <- function(expr)
{
	expr[[2]]
}
`Simplify.+` <- function(expr)
{
	if (length(expr) == 2)
	{
		return(expr[[2]])
	}
	a <- expr[[2]]
	b <- expr[[3]]

	if (a == 0) {
		b
	} else if (b == 0) {
		a
	} else if (is.uminus(b)) {
		call("-", a, b[[2]])
	} else {
		expr
	}
}

`Simplify.-` <- function(expr)
{
	if (length(expr) == 2)
	{
		if (is.uminus(expr[[2]])) {
			return(Simplify_(expr[[2]][[2]]))
		} else if (is.uplus(expr[[2]])) {
			return(Simplify_(substitute(-expr[[2]][[2]])))
		}
		a <- 0
		b <- expr[[2]]
	} else {
		a <- expr[[2]]
		b <- expr[[3]]
	}

	if (a == 0) {
		if (is.uminus(b)) b[[2]] else substitute(-b)
	} else if (b == 0) {
		a
	} else if (is.uminus(b)) {
		call("+", a, b[[2]])
	} else {
		expr
	}
}

`Simplify.*` <- function(expr, div=FALSE)
{
	a <- expr[[2]]
	b <- expr[[3]]
	if (is.uminus(a)) {
		sminus <- TRUE
		a <- a[[2]]
	} else {
		sminus <- FALSE
	}
	if (is.uminus(b)) {
		sminus <- !sminus
		b <- b[[2]]
	}
#browser()
	if (a == 0 || (b == 0 && !div)) {
		0
	} else if (a == 1 && !div) {
		if (sminus) substitute(-b) else b
	} else if (b == 1) {
		if (sminus) substitute(-a) else a
	} else {
#browser()
		# get numerator and denominator for a and b than combine them
		nd_a <- Numden(a)
		nd_b <- Numden(b)
		if (div) {
			nd <- list(num=c(nd_a$num, nd_b$den), den=c(nd_a$den, nd_b$num))
		} else {
			nd <- list(num=c(nd_a$num, nd_b$num), den=c(nd_a$den, nd_b$den))
		}
		# reduce numerics to only one factor
		fa=list()
		for (na in c("num", "den")) {
			inu=if (length(nd[[na]])) which(sapply(nd[[na]], is.numeric)) else integer(0)
			if (length(inu)) {
				fa[[na]] <- prod(unlist(nd[[na]][inu]))
				# remove numerics
				nd[[na]] <- nd[[na]][-inu]
			} else {
				fa[[na]] <- 1
			}
			# make factors positive
			if (fa[[na]] < 0) {
				fa[[na]] <- -fa[[na]]
				sminus <- !sminus
			}
		}
		if (fa$num == 0) {
			return(0)
		}
		if ((as.integer(fa$num) != fa$num ||
				as.integer(fa$den) != fa$den) ||
				fa$num == fa$den || fa$num == -fa$den) {
			fa$num <- fa$num/fa$den
			fa$den <- 1
		}
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
			if (length(nd[[na]]) == 0 && fa[[na]] == 1)
				next
			eprod[[na]]=if (length(nd[[na]])) nd[[na]][[1]] else fa[[na]]
			for (term in nd[[na]][-1]) {
				eprod[[na]] <- call("*", eprod[[na]], term)
			}
			if (length(nd[[na]]) && fa[[na]] != 1) {
				eprod[[na]] <- if (fa[[na]] == -1) call("-", eprod[[na]]) else call("*", fa[[na]], eprod[[na]])
			}
		}
		eprod$num=if (is.null(eprod$num)) fa$num else eprod$num
		if (is.null(eprod$den)) {
			# we have no denominator
			expr <- eprod$num
		} else {
			expr <- call("/", eprod$num, eprod$den)
		}
		return(if (sminus) substitute(-expr) else expr)
	}
}
`Simplify./` <- function(expr)
{
	`Simplify.*`(expr, div=TRUE)
}
`Simplify.^` <- function(expr)
{
	a <- expr[[2]]
	b <- expr[[3]]

	if (a == 0) {
		0
	} else if (b == 0 || a == 1) {
		1
	} else if (b == 1) {
		a
	} else if (b == 0.5) {
		substitute(sqrt(a))
	} else if (b == -0.5) {
		substitute(1/sqrt(a))
	} else if (is.call(a)) {
		if (as.character(a[[1]]) == "^") {
			# product of exponents
			b <- Simplify_(call("*", a[[3]], b))
			a <- a[[2]]
		} else if (as.character(a[[1]]) == "sqrt") {
			# divide by 2
			b <- Simplify_(call("/", b, 2))
			a <- a[[2]]
		}
		expr[[2]] <- a
		expr[[3]] <- b
		expr
	} else {
		expr
	}
}

Numden <- function(expr) {
	# Return a list with "num" as numerator and "den" as denominator sublists.
	# Each sublist regroups the language expressions which are not products neither
	# divisions
	if (is.uminus(expr)) {
		if (is.numeric(expr[[2]])) {
			list(num=list(-expr[[2]]))
		} else {
			a=Numden(expr[[2]])
			a$num=c(-1, a$num)
			a
		}
	} else if (is.uplus(expr)) {
		Numden(expr[[2]])
	} else if (is.symbol(expr) || is.numeric(expr)) {
		list(num=list(expr))
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
		list(num=list(expr))
	}
}
format1 <- function(expr) {
	res <- if (is.symbol(expr)) as.character(expr) else format(expr)
	if (length(res) > 1) {
		res=paste(res, collapse="")
	}
	return(res)
}
is.uminus <- function(e) {
	# detect if e is unitary minus, e.g. "-a"
	return(is.call(e) && length(e) == 2 && e[[1]] == as.symbol("-"))
}
is.uplus <- function(e) {
	# detect if e is unitary plus, e.g. "+a"
	return(is.call(e) && length(e) == 2 && e[[1]] == as.symbol("+"))
}
is.unumeric <- function(e) {
	# detect if numeric with optional unitary sign(s)
	return(is.numeric(e) || ((is.uminus(e) || is.uplus(e)) && is.unumeric(e[[2]])))
}
