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
			nd <- list(
				num=list(b=c(nd_a$num$b, nd_b$den$b),
				p=c(nd_a$num$p, nd_b$den$p)),
				den=list(b=c(nd_a$den$b, nd_b$num$b),
				p=c(nd_a$den$p, nd_b$num$p))
			)
		} else {
			nd <- list(
				num=list(b=c(nd_a$num$b, nd_b$num$b),
				p=c(nd_a$num$p, nd_b$num$p)),
				den=list(b=c(nd_a$den$b, nd_b$den$b),
				p=c(nd_a$den$p, nd_b$den$p))
			)
		}
		# reduce numerics to only one factor
		fa=list()
		for (na in c("num", "den")) {
			inu=if (length(nd[[na]]$b)) which(sapply(nd[[na]]$b, is.numeric) &
				sapply(nd[[na]]$p, is.numeric)) else integer(0)
			if (length(inu)) {
				# the power of numerics are always 1 after Simplify_()
				fa[[na]] <- prod(unlist(nd[[na]]$b[inu]))
				# remove numerics
				nd[[na]]$b <- nd[[na]]$b[-inu]
				nd[[na]]$p <- nd[[na]]$p[-inu]
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
		# group identical bases by adding their powers
#browser()
		for (na in c("num", "den")) {
			if (length(nd[[na]]$b) <= 1)
				next
			nd_eq <- outer(sapply(nd[[na]]$b, format1), sapply(nd[[na]]$b, format1), `==`)
			for (inum in seq(len=ncol(nd_eq))) {
				isim <- which(nd_eq[,inum])
				isim <- isim[isim > inum & sapply(nd[[na]]$p[isim], `!=`, 0)]
				if (length(isim)) {
					# add powers for this base
					p <- nd[[na]]$p[[inum]]
					for (i in isim) {
						p <- call("+", p, nd[[na]]$p[[i]])
					}
					nd[[na]]$p[[inum]] <- Simplify_(p)
					# set grouped powers to 0
					nd[[na]]$p[isim] <- 0
				}
			}
			# remove power==0 terms
			ize=which(sapply(nd[[na]]$p, `==`, 0))
			if (length(ize)) {
				nd[[na]]$b <- nd[[na]]$b[-ize]
				nd[[na]]$p <- nd[[na]]$p[-ize]
			}
		}
		# simplify identical terms in num and denum by subtracting powers
		nd_eq <- outer(sapply(nd$den$b, format1), sapply(nd$num$b, format1), `==`)
		ipair <- matrix(0, nrow=2, ncol=0)
		for (inum in seq(len=ncol(nd_eq))) {
			iden <- which(nd_eq[,inum]) # of length at most 1 as terms are already grouped
			iden <- iden[!iden %in% ipair[2,]]
			if (length(iden)) {
				# simplify power for this pair
				ipair <- cbind(ipair, c(inum, iden))
				nd$num$p[[inum]] <- Simplify_(call("-", nd$num$p[[inum]], nd$den$p[[iden]]))
				nd$den$p[[iden]] <- 0
			}
		}
		if (ncol(ipair) > 0) {
			# remove power==0 terms
			for (na in c("num", "den")) {
				ize=which(sapply(nd[[na]]$p, `==`, 0))
				if (length(ize)) {
					nd[[na]]$b <- nd[[na]]$b[-ize]
					nd[[na]]$p <- nd[[na]]$p[-ize]
				}
			}
		}
		# form symbolic products
		eprod <- list()
		for (na in c("num", "den")) {
			if (length(nd[[na]]$b) == 0 && fa[[na]] == 1)
				next
			# remove power==0 terms
			ize=sapply(nd[[na]]$p, `==`, 0)
			nd[[na]]$b <- nd[[na]]$b[!ize]
			nd[[na]]$p <- nd[[na]]$p[!ize]
			if (length(nd[[na]]$b) == 0)
				next
			eprod[[na]] <- if (length(nd[[na]])) Simplify_(call("^", nd[[na]]$b[[1]], nd[[na]]$p[[1]])) else fa[[na]]
			for (i in 1+seq_along(nd[[na]]$b[-1])) {
				term <- Simplify_(call("^", nd[[na]]$b[[i]], nd[[na]]$p[[i]]))
				eprod[[na]] <- call("*", eprod[[na]], term)
			}
			if (length(nd[[na]]$b) && fa[[na]] != 1) {
				eprod[[na]] <- if (fa[[na]] == -1) call("-", eprod[[na]]) else call("*", fa[[na]], eprod[[na]])
			}
		}
		eprod$num <- if (is.null(eprod$num)) fa$num else eprod$num
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
	# divisions. The terms are decomposed in b^p sublists
	if (is.uminus(expr)) {
		a=Numden(expr[[2]])
		a$num$b=c(-1, a$num$b)
		a$num$p=c(1, a$num$p)
		a
	} else if (is.uplus(expr)) {
		Numden(expr[[2]])
	} else if (is.symbol(expr) || is.numeric(expr)) {
		list(num=list(b=expr, p=1))
	} else if (expr[[1]] == as.symbol("*")) {
		# recursive call
		a=Numden(expr[[2]])
		b=Numden(expr[[3]])
		list(num=list(b=c(a$num$b, b$num$b), p=c(a$num$p, b$num$p)),
			den=list(b=c(a$den$b, b$den$b), p=c(a$den$p, b$den$p)))
	} else if (expr[[1]] == as.symbol("/")) {
		# recursive call
		a=Numden(expr[[2]])
		b=Numden(expr[[3]])
		list(num=list(b=c(a$num$b, b$den$b), p=c(a$num$p, b$den$p)),
			den=list(b=c(a$den$b, b$num$b), p=c(a$den$p, b$num$p)))
	} else if (is.call(expr) && expr[[1]] == as.symbol("^")) {
		if (expr[[3]] < 0) {
			# make the power look positive
			list(den=list(b=expr[[2]], p=if (is.numeric(expr[[3]])) -expr[[3]] else expr[[3]][[2]]))
		} else {
			list(num=list(b=expr[[2]], p=expr[[3]]))
		}
	} else {
		list(num=list(b=expr, p=1))
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
