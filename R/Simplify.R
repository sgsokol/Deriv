# Simplify.R -- symbollic simplification
# written by Andrew Clausen <clausen@econ.upenn.edu> in 2007
# thanks to a bug fix from Mark Reid <mark.reid@anu.edu.au> in 21/2/2009
#
# This isn't a serious attempt at simplification code.  It just does some
# obvious things like 0 + x => x.  It was written to support Deriv.R.

Simplify_ <- function(expr)
{
	if (is.symbol(expr) || is.numeric(expr)) {
		expr
	} else if (is.language(expr) && is.symbol(expr[[1]])) {
		# is there a rule in the table?
		sym.name <- as.character(expr[[1]])
		if (class(try(Simplify.rule <-
				get(sym.name, envir=simplifications,
					inherits=FALSE), silent=TRUE))
					!= "try-error")
			return(Simplify.rule(expr))
		else if (all(sapply(expr[-1], function(arg) is.numeric(Simplify_(arg)))))
			return(eval(expr))
		# if all arguments are numeric, evaluate them
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
	if (length(expr) == 2)
	{
		if (is.numeric(expr[[2]]))
			return(expr[[2]])
		else {
			return(Simplify_(expr[[2]]))
		}
	}
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
		if (is.uminus(expr[[2]])) {
			return(Simplify_(expr[[2]][[2]]))
		} else if (is.uplus(expr[[2]])) {
			return(Simplify_(substitute(-expr[[2]][[2]])))
		}
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
	Simplify_(expr[[2]])

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
	if (is.numeric(a) && a == 0) {
		0
	} else if (is.numeric(b) && b == 0 && !div) {
		0
	} else if (is.numeric(a) && a == 1 && !div) {
		if (sminus) Simplify_(substitute(-b)) else Simplify_(b)
	} else if (is.numeric(b) && b == 1) {
		if (sminus) Simplify_(substitute(-a)) else Simplify_(a)
	} else if (is.numeric(a) && is.numeric(b)) {
		res <- if (div) a/b else a * b
		if (sminus) -res else res
	} else {
#browser()
		# get numerator and denominator for a and b than combine them
		nd_a <- Numden(Simplify_(a))
		nd_b <- Numden(Simplify_(b))
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
		fa$num=if (sminus) -fa$num else fa$num
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
			if (length(nd[[na]]) == 0 && fa[[na]] == 1) next
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
			return(eprod$num)
		} else {
			return(call("/", eprod$num, eprod$den))
		}
	}
}
`Simplify./` <- function(expr)
{
	`Simplify.*`(expr, div=TRUE)
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
		if (is.call(a)) {
			if (as.character(a[[1]]) == "^") {
				# product of exponents
				b <- Simplify_(call("*", a[[3]], b))
				a <- a[[2]]
			} else if (as.character(a[[1]]) == "sqrt") {
				# divide by 2
				b <- Simplify_(call("/", b, 2))
				a <- a[[2]]
			}
		}
		if (is.numeric(b)) {
			if (all(b == 1)) {
				return(a)
			} else if (all(b == 0.5)) {
				return(substitute(sqrt(a)))
			} else if (all(b == -0.5)) {
				return(substitute(1/sqrt(a)))
			}
		}
		expr[[2]] <- a
		expr[[3]] <- b
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
