#' @name Simplify
#' @title Symbollic simplification of an expression or function
#' @aliases Simplify simplifications
#' @concept symbolic simplification
# \usage{
# Simplify(expr, env=parent.frame())
# }
#' 
#' 
#' @param expr An expression to be simplified, expr can be
#' \itemize{
#'    \item an expression: \code{as.expression(x+x)}
#'    \item an string: \code{"x+x"}
#'    \item a function: \code{function(x) x+x}
#'    \item a right hand side of a formula: \code{~x+x}
#'    \item a language: \code{quote(x+x)}
#' }
#' @param env An environment in wich a simplified function is created
#'  if \code{expr} is a function. This argument is ignored is all other cases.
#' @return A simplified expression. The result is of the same type as
#'  \code{expr} except for formula, where a language is returned.
#' @details An environment \code{simplifications} containing simplification rules, is exported in the user namespace.
Simplify <- function(expr, env=parent.frame())
	if (is.expression(expr)) {
		as.expression(Simplify_(expr[[1]]))
	} else if (is.function(expr)) {
		as.function(c(as.list(formals(expr)),
			Simplify_(body(expr))),
			envir=env)
	} else if (is.call(expr) && expr[[1]] == as.symbol("~")) {
		Simplify_(expr[[length(expr)]])
	} else if (is.character(expr)) {
		format1(Simplify_(parse(text=expr)[[1]]))
	} else {
		Simplify_(expr)
	}
#' @name format1
#' @title Wrapper for base::format() function
# \usage{
# format1(expr)
# }
#' 
#' 
#' @param expr An expression or symbol or language to be converted to a string.
#' @return A character vector of length 1 contrary to base::format() which
#'  can split its output over several lines.
format1 <- function(expr) {
	res <- if (is.symbol(expr)) as.character(expr) else format(expr)
	if (length(res) > 1) {
		res=paste(res, collapse="")
	}
	return(res)
}

Simplify_ <- function(expr)
{
	if (is.unumeric(expr)) {
		eval(expr)
	} else if (is.call(expr)) {
		args <- lapply(expr[-1], Simplify_)
		if (all(sapply(args, is.numeric))) {
			# if all arguments are numeric, evaluate them
			return(eval(as.call(c(expr[[1]], args))))
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


# in what follows no need to Simplify_ args neither to check if
# all arguments are unumeric. It is done in upper Simplify_()
`Simplify.(` <- function(expr)
{
	expr[[2]]
}
`Simplify.+` <- function(expr, add=TRUE)
{
	if (length(expr) == 2)
	{
		if (add)
			return(expr[[2]])
		else if (is.uminus(expr[[2]]))
			return(expr[[2]][[2]])
		else if (is.uplus(expr[[2]]))
			return(call("-", expr[[2]][[2]]))
		else
			return(expr)
	}
	a <- expr[[2]]
	b <- expr[[3]]
	sminus <- FALSE
	
	if (a == 0) {
		return(if (add) b else call("-", b))
	} else if (b == 0) {
		return(a)
	} else if (format1(a) == format1(b) && !add) {
		return(0)
	} else if (is.uminus(b)) {
		b <- b[[2]]
		add <- !add
		if (!add && is.uminus(a)) {
			a <- a[[2]]
			expr <- call("-", call("+", a, b))
			sminus <- TRUE
			add <- TRUE
		} else {
			expr <- call(if (add) "+" else "-", a, b)
		}
	}
	# group and simplify identical terms
	alc <- Lincomb(a)
	blc <- Lincomb(b)
	if (add)
		lc <- list(co=c(alc$co, blc$co), it=c(alc$it, blc$it))
	else
		lc <- list(co=c(alc$co, -blc$co), it=c(alc$it, blc$it))
	if (sminus)
		lc$co <- -lc$co # sminus is no more used here after
	# group all numerics
	inum <- sapply(lc$it, `==`, 1)
	nu <- sum(lc$co[inum])
	# remove numerics from it list
	lc$it <- lc$it[!inum]
	lc$co <- lc$co[!inum]
	# group identical terms
	itch <- sapply(lc$it, format1)
	ta <- table(itch) # count of unique entries
	if (sum(inum) <= 1 && all(ta == 1)) {
		# nothing to simplify
		return(expr)
	}
	tsim <- outer(itch, names(ta), `==`)
	isim <- lapply(seq_len(ncol(tsim)), function(i) which(tsim[,i]))
	# keep only unique terms
	irm <- c()
	for (i in isim) {
		if (length(i) == 0)
			next
		lc$co[i[1]] <- sum(lc$co[i])
		irm <- c(irm, i[-1])
	}
	# remove grouped terms but first
	if (length(irm)) {
		lc$co <- lc$co[-irm]
		lc$it <- lc$it[-irm]
		itch <- itch[irm]
	}
	
	# remove co==0
	inul <- lc$co == 0
	lc$co <- lc$co[!inul]
	lc$it <- lc$it[!inul]
	itch <- itch[!inul]
	if (length(lc$it) == 0) {
		return(0)
	}
#browser()
	ipn=list(pos=lc$co > 0, neg=lc$co < 0)
	for (pn in c("pos", "neg")) {
		# where abs(co) != 1, replace term by nb_repeat*term
		i <- ipn[[pn]] & abs(lc$co[ipn[[pn]]]) != 1
		lc$it[i] <- lapply(which(i), function(ii) {
			e <- lc$it[[ii]]
			if (is.call(e) && e[[1]] == as.symbol("/") && e[[2]] == 1)
				e[[2]] <- abs(lc$co[ii])
			else
				e <- Simplify_(call("*", abs(lc$co[ii]), e))
			return(e)
		})
	}
	if (nu != 0) {
		lc$it <- c(lc$it, abs(nu))
		lc$co <- c(lc$co, sign(nu))
		itch <- c(itch, format1(abs(nu)))
	}
	# form final symbolic expression
	iord <- order(lc$co < 0, itch) # positives first
	expr <- lc$it[[iord[1]]]
	sminus <- lc$co[iord[1]] < 0 # all negatives
	
	for (i in iord[-1]) {
		expr <- call(if (sminus || lc$co[i] > 0) "+" else "-", expr, lc$it[[i]])
	}
	if (sminus)
		return(call("-", expr))
	else
		return(expr)
}

`Simplify.-` <- function(expr)
{
	`Simplify.+`(expr, add=FALSE)
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
	} else if (div && format1(a) == format1(b)) {
		if (sminus) -1 else 1
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
			sminus=xor(sminus, xor(nd_a$sminus, nd_b$sminus))
		} else {
			nd <- list(
				num=list(b=c(nd_a$num$b, nd_b$num$b),
				p=c(nd_a$num$p, nd_b$num$p)),
				den=list(b=c(nd_a$den$b, nd_b$den$b),
				p=c(nd_a$den$p, nd_b$den$p))
			)
			sminus=xor(sminus, xor(nd_a$sminus, nd_b$sminus))
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
			if (fa[[na]] < 0) {
				sminus = !sminus
				fa[[na]] <- -fa[[na]]
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
#browser()
		# form symbolic products
		eprod <- list()
		for (na in c("num", "den")) {
			if (length(nd[[na]]$b) == 0)
				next
			# remove power==0 terms
			ize=sapply(nd[[na]]$p, `==`, 0)
			nd[[na]]$b <- nd[[na]]$b[!ize]
			nd[[na]]$p <- nd[[na]]$p[!ize]
			if (length(nd[[na]]$b) == 0)
				next
			for (i in seq_along(nd[[na]]$b)) {
				p <- nd[[na]]$p[[i]]
				term <- if (p == 1) nd[[na]]$b[[i]] else Simplify_(call("^", nd[[na]]$b[[i]], p))
				if (is.null(eprod[[na]]))
					eprod[[na]] <- term # start the sequence
				else
					eprod[[na]] <- call("*", eprod[[na]], term)
			}
		}
		expr <- if (is.null(eprod$num)) 1 else eprod$num
		if (!is.null(eprod$den)) {
			expr <- call("/", expr, eprod$den)
		}
		# put numeric factor at first place
		if (fa$num != 1 && fa$den != 1) {
			# add to both num. and denom.
			if (!is.null(eprod$den)) {
				expr[[2]] <- call("*", fa$num, expr[[2]])
				expr[[3]] <- call("*", fa$den, expr[[3]])
			} else {
				expr <- call("/", call("*", fa$num, expr), fa$den)
			}
		} else if (fa$num != 1) {
			if (is.call(expr) && expr[[1]] == as.symbol("/") && expr[[2]] == 1)
				expr[[2]] <- fa$num
			else
				expr <- call("*", fa$num, expr)
		} else if (fa$den != 1) {
			if (is.call(expr) && expr[[1]] == as.symbol("/"))
				expr[[3]] <- call("*", fa$den, expr[[3]])
			else
				expr <- call("/", expr, fa$den)
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
		if (a[[1]] == as.symbol("^")) {
			# product of exponents
			b <- Simplify_(call("*", a[[3]], b))
			a <- a[[2]]
		} else if (a[[1]] == as.symbol("sqrt")) {
			# divide by 2
			b <- Simplify_(call("/", b, 2))
			a <- a[[2]]
		} else if (a[[1]] == as.symbol("abs") && is.numeric(b) && b%%2 == 0) {
			# remove abs() for even power
			a <- a[[2]]
		}
		expr[[2]] <- a
		expr[[3]] <- b
		expr
	} else {
		expr
	}
}
Simplify.log <- function(expr) {
	if (is.call(expr[[2]])) {
		# the argument of log is a function
		if (expr[[2]][[1]] == as.symbol("^")) {
			p <- expr[[2]][[3]]
			expr[[2]] <- expr[[2]][[2]]
			expr <- Simplify_(call("*", p, expr))
		} else if (expr[[2]][[1]] == as.symbol("exp")) {
			if (length(expr) == 2)
				expr <- expr[[2]][[2]]
			else
				expr <- Simplify_(call("/", expr[[2]][[2]], call("log", expr[[3]])))
		} else if (expr[[2]][[1]] == as.symbol("sqrt")) {
			expr[[2]] <- expr[[2]][[2]]
			expr <- Simplify_(call("*", 0.5, expr))
		} else if (expr[[2]][[1]] == as.symbol("*")) {
			a <- expr
			a[[2]] <- expr[[2]][[2]]
			expr[[2]] <- expr[[2]][[3]] # unitary "+" cannot appear here
			expr <- Simplify_(call("+", a, expr))
		} else if (expr[[2]][[1]] == as.symbol("/")) {
			a <- expr
			a[[2]] <- expr[[2]][[2]]
			expr[[2]] <- expr[[2]][[3]] # unitary "+" cannot appear here
			expr <- Simplify_(call("-", a, expr))
		} else {
			expr
		}
	}
	if (is.call(expr) && expr[[1]] == as.symbol("log") && length(expr) == 3 && format1(expr[[2]]) == format1(expr[[3]])) {
		1
	} else {
		expr
	}
}
Simplify.sqrt <- function(expr) {
	if (is.call(expr[[2]])) {
		# the argument of sqrt is a function
		if (expr[[2]][[1]] == as.symbol("^")) {
			p <- expr[[2]][[3]]
			Simplify_(call("^",  call("abs", expr[[2]][[2]]), call("/", p, 2)))
		} else if (expr[[2]][[1]] == as.symbol("exp")) {
			expr[[2]][[2]] <- Simplify_(call("/", expr[[2]][[2]], 2))
			expr[[2]]
		} else if (expr[[2]][[1]] == as.symbol("sqrt")) {
			Simplify_(call("^", expr[[2]][[2]], 0.25))
		} else if (expr[[2]][[1]] == as.symbol("*") && format1(expr[[2]][[2]]) == format1(expr[[2]][[3]])) {
			Simplify_(call("abs", expr[[2]][[2]]))
		} else {
			expr
		}
	} else {
		expr
	}
}
Simplify.abs <- function(expr) {
	if (is.uminus(expr[[2]])) {
		expr[[2]] <- expr[[2]][[2]]
	} else if (is.call(expr[[2]])) {
		if (expr[[2]][[1]] == as.symbol("^")) {
			p <- expr[[2]][[3]]
			if (is.numeric(p) && p%%2 == 0)
				expr <- expr[[2]]
		} else if (expr[[2]][[1]] == as.symbol("exp") || expr[[2]][[1]] == as.symbol("sqrt")) {
			expr <- expr[[2]]
		}
	}
	expr
}

Numden <- function(expr) {
	# Return a list with "num" as numerator and "den" as denominator sublists.
	# Each sublist regroups the language expressions which are not products neither
	# divisions. The terms are decomposed in b^p sublists
	if (is.uminus(expr)) {
		a=Numden(expr[[2]])
		a$sminus <- !a$sminus
		a
	} else if (is.uplus(expr)) {
		Numden(expr[[2]])
	} else if (is.symbol(expr)) {
		if (is.numeric(expr)) {
			sminus <- expr < 0
			expr <- if (sminus) -expr else expr
		} else {
			sminus <- FALSE
		}
		list(num=list(b=expr, p=1), sminus=sminus)
	} else if (is.call(expr)) {
		if (expr[[1]] == as.symbol("*")) {
			# recursive call
			a=Numden(expr[[2]])
			b=Numden(expr[[3]])
			list(num=list(b=c(a$num$b, b$num$b), p=c(a$num$p, b$num$p)),
				den=list(b=c(a$den$b, b$den$b), p=c(a$den$p, b$den$p)),
				sminus=xor(a$sminus, b$sminus))
		} else if (expr[[1]] == as.symbol("/")) {
			# recursive call
			a=Numden(expr[[2]])
			b=Numden(expr[[3]])
			list(num=list(b=c(a$num$b, b$den$b), p=c(a$num$p, b$den$p)),
				den=list(b=c(a$den$b, b$num$b), p=c(a$den$p, b$num$p)),
				sminus=xor(a$sminus, b$sminus))
		} else if (expr[[1]] == as.symbol("^")) {
			if (expr[[3]] < 0) {
				# make the power look positive
				list(den=list(b=expr[[2]], p=if (is.numeric(expr[[3]])) -expr[[3]] else expr[[3]][[2]]), sminus=FALSE)
			} else {
				list(num=list(b=expr[[2]], p=expr[[3]]), sminus=FALSE)
			}
		} else {
			list(num=list(b=expr, p=1), sminus=FALSE)
		}
	} else {
		list(num=list(b=expr, p=1), sminus=FALSE)
	}
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
Lincomb <- function(expr) {
	# decompose sum and diff in linear combibnation of num.coeff*item paires.
	# For numerics, item is 1
	# numerical factor is supposed to be at first place in products
	if (is.uminus(expr)) {
		a <- Lincomb(expr[[2]])
		a$co=-a$co
		a
	} else if (is.uplus(expr)) {
		Lincomb(expr[[2]])
	} else if (is.symbol(expr)) {
		list(it=list(expr), co=1)
	} else if (is.numeric(expr)) {
		list(co=expr, it=1)
	} else if (is.call(expr)) {
		if (expr[[1]] == as.symbol("+")) {
			# recursive call
			a <- Lincomb(expr[[2]])
			b <- Lincomb(expr[[3]])
			list(co=c(a$co, b$co), it=c(a$it, b$it))
		} else if (expr[[1]] == as.symbol("-")) {
			# recursive call
			a <- Lincomb(expr[[2]])
			b <- Lincomb(expr[[3]])
			list(co=c(a$co, -b$co), it=c(a$it, b$it))
		} else if (expr[[1]] == as.symbol("*") && is.numeric(expr[[2]])) {
			list(co=expr[[2]], it=expr[[3]])
		} else if (expr[[1]] == as.symbol("/")) {
			if (is.numeric(expr[[2]]))
				list(it=list(call("/", 1, expr[[3]])), co=expr[[2]])
			else if (is.call(expr[[2]]) && expr[[2]][[1]] == as.symbol("*") && is.numeric(expr[[2]][[2]]))
				list(it=list(call("/", expr[[2]][[3]]), expr[[3]]), co=expr[[2]][[2]])
			else
				list(co=1, it=list(expr))
		} else {
			list(co=1, it=list(expr))
		}
	} else {
		list(co=1, it=list(expr))
	}
}

simplifications <- new.env()

assign("+", `Simplify.+`, envir=simplifications)
assign("-", `Simplify.-`, envir=simplifications)
assign("*", `Simplify.*`, envir=simplifications)
assign("/", `Simplify./`, envir=simplifications)
assign("(", `Simplify.(`, envir=simplifications)
assign("^", `Simplify.^`, envir=simplifications)
assign("log", `Simplify.log`, envir=simplifications)
assign("logb", `Simplify.log`, envir=simplifications)
assign("sqrt", `Simplify.sqrt`, envir=simplifications)
assign("abs", `Simplify.abs`, envir=simplifications)
