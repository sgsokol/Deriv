#' @name Simplify
#' @title Symbollic simplification of an expression or function
#' @description Symbollic simplification of an expression or function
#' @aliases Simplify simplifications Cache deCache
#' @concept symbolic simplification
# \usage{
# Simplify(expr, env=parent.frame(), scache=new.env())
# }
#' 
#' 
#' @param expr An expression to be simplified, expr can be
#' \itemize{
#'    \item an expression: \code{expression(x+x)}
#'    \item a string: \code{"x+x"}
#'    \item a function: \code{function(x) x+x}
#'    \item a right hand side of a formula: \code{~x+x}
#'    \item a language: \code{quote(x+x)}
#' }
#' @param env An environment in which a simplified function is created
#'  if \code{expr} is a function. This argument is ignored in all other cases.
#' @param scache An environment where there is a list in which simplified expression are cached
#' @param st A language expression to be cached
#' @param prefix A string to start the names of the cache variables
#' @return A simplified expression. The result is of the same type as
#'  \code{expr} except for formula, where a language is returned.
#' @details An environment \code{simplifications} containing simplification rules, is exported in the namespace accessible by the user.
#'  Cache() is used to remove redundunt calculations by storing them in
#'  cache variables. Default parameters to Cache() does not have to be provided
#'  by user. deCache() makes the inverse job -- a series of assignements
#'  are replaced by only one big expression without assignement.
#'  Sometimes it is usefull to
#'  apply deChache() and only then pass its result to Cache().
Simplify <- function(expr, env=parent.frame(), scache=new.env()) {
	if (is.null(scache$l))
		scache$l <- list() # for stand alone use of Simplify
	if (is.expression(expr)) {
		as.expression(Simplify_(expr[[1]], scache))
	} else if (is.function(expr)) {
		as.function(c(as.list(formals(expr)),
			Simplify_(body(expr), scache)),
			envir=env)
	} else if (is.call(expr) && expr[[1]] == as.symbol("~")) {
		Simplify_(expr[[length(expr)]], scache)
	} else if (is.character(expr)) {
		format1(Simplify_(parse(text=expr)[[1]], scache))
	} else {
		Simplify_(expr, scache)
	}
}

#' @name format1
#' @title Wrapper for base::format() function
#' @description Wrapper for base::format() function
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
	n <- length(res)
	if (n > 1) {
		if (res[1] == "{" && res[n] == "}" && n > 3) {
			b <- paste(res[-c(1,n)], collapse="; ")
			res <- paste("{", b, "}", collapse="")
		} else {
			res <- paste(res, collapse="")
		}
	}
	return(res)
}

Simplify_ <- function(expr, scache) {
	if (is.call(expr)) {
		che <- format1(expr)
		res <- scache$l[[che]]
		if (!is.null(res)) {
			if (typeof(res) == "logical" && is.na(res)) {
				# recursive infinite call
				scache$l[[che]] <- expr
				return(expr)
			} else {
				return(res)
			}
		}
		scache$l[[che]] <- NA # token holder
#cat("simp expr=", format1(expr), "\n", sep="")
		args <- lapply(as.list(expr)[-1], Simplify_, scache)
		expr[-1]=args
		if (all(sapply(args, is.conuloch))) {
			# if all arguments are like numeric, evaluate them
			res <- eval(expr)
			scache$l[[che]] <- res
			return(res)
		} else {
			# is there a rule in the table?
			sym.name <- as.character(expr[[1]])
			Simplify.rule <- simplifications[[sym.name]]
			res <- if (!is.null(Simplify.rule)) Simplify.rule(expr, scache=scache) else expr
			scache$l[[che]] <- res
			return(res)
		}
	} else {
		expr
	}
}

# in what follows no need to Simplify_ args neither to check if
# all arguments are numeric. It is done in the upper Simplify_()
`Simplify.(` <- function(expr, scache=NULL) {
	expr[[2]]
}
`Simplify.+` <- function(expr, add=TRUE, scache=NULL) {
	if (length(expr) == 2) {
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
	
	if (a == 0) {
		return(if (add) b else call("-", b))
	} else if (b == 0) {
		return(a)
	} else if (add && is.uminus(a) && !is.uminus(b)) {
		a <- b
		b <- expr[[2]][[2]]
		add <- FALSE
		expr <- call("-", a, b)
	} else if (identical(a, b)) {
		return(if (add) Simplify_(call("*", 2, a), scache) else 0)
	} else if (!is.call(a) && !is.call(b)) {
		if (add) {
			# just reorder
			expr[-1] <- expr[1+order(sapply(expr[-1], as.character))]
		}
		return(expr)
	}
	# factorise most repeated terms
	alc <- Lincomb(a)
	blc <- Lincomb(b)
	if (add) {
		lc <- c(alc, blc)
	} else {
		# inverse sminus in b
		blc <- lapply(blc, function(it) {it$sminus <- !it$sminus; it})
		lc <- c(alc, blc)
	}
#browser()
	# sum purely numeric terms
	inum <- which(sapply(lc, function(it) length(it$num)==0 && length(it$den)==0))
	if (length(inum) > 1) {
		term <- sum(sapply(lc[inum], function(it) (if (it$sminus) -1 else 1)*it$fa$num/it$fa$den))
		lc[[inum[1]]] <- list(fa=list(num=abs(term), den=1), sminus=term<0)
		lc <- lc[-inum[-1]]
	}
	bch <- ta <- tsim <- po <- ilc <- ind <- list()
	for (cnd in c("num", "den")) {
		# character bases in num/den
		bch[[cnd]] <- unlist(lapply(lc, function(it) {lapply(it[[cnd]]$b, format1)}))
		# powers
		po[[cnd]] <- do.call(c, lapply(lc, function(it) it[[cnd]]$p), quote=TRUE)
		# index of the lc term for each bnch
		ta[[cnd]] <- table(bch[[cnd]])
		ta[[cnd]] <- ta[[cnd]][ta[[cnd]] > 1] # keep only repeated bases
		tsim[[cnd]] <- outer(bch[[cnd]], names(ta[[cnd]]), `==`)
		ilc[[cnd]] <- unlist(lapply(seq_along(lc), function(i) {rep(i, length(lc[[i]][[cnd]]$b))}))
		# index of the base in a given term (nd) for each bnch
		ind[[cnd]] <- unlist(lapply(seq_along(lc), function(i) {seq_along(lc[[i]][[cnd]]$b)}))
	}
#browser()
	# fnd will be the name "num" or "den" where the first factor
	# will be taken. ond is the "other" name (if fnd=="num", then ond == "den")
	# we select the cadidate which is most repeated provided that it
	# has at least one numeric power occurance.
	taa <- unlist(ta)
	ota <- order(taa, decreasing=TRUE)
	ntan <- length(ta$num)
	fnd <- NA
	for (i in ota) {
		cnd <- if (i > ntan) "den" else "num"
		ita <- i - if (i > ntan) ntan else 0
		ib <- bch[[cnd]] == names(ta[[cnd]])[ita]
		if (any(sapply(po[[cnd]][ib], is.numeric))) {
			fnd <- cnd
			iit <- which(ib) # the bases equal to factor
			p_fa <- min(sapply(po[[cnd]][ib], function(p) if (is.numeric(p)) p else NA), na.rm=TRUE)
			i_lc <- ilc[[cnd]][iit]
			i_nd <- ind[[cnd]][iit]
			break
		}
	}
#browser()
	if (is.na(fnd))
		return(lc2expr(lc, scache)) # nothing to factorize, just order terms
	ond <- if (fnd == "num") "den" else "num"
	# create nd with the first factor
	fa_nd <- list(num=list(b=list(), p=list()),
		den=list(b=list(), p=list()),
		sminus=FALSE, fa=list(num=1, den=1))
	fa_nd[[fnd]]$b <- lc[[i_lc[1]]][[fnd]]$b[i_nd[1]]
	fa_nd[[fnd]]$p <- list(p_fa)
	# decrease p in the lc terms
	for (i in seq_along(i_lc)) {
		lc[[i_lc[i]]][[fnd]]$p[[i_nd[i]]] <- Simplify_(call("-", lc[[i_lc[i]]][[fnd]]$p[[i_nd[i]]], p_fa), scache)
	}
	
	for (cnd in c(fnd, ond)) {
		# see if other side can provide factors
		for (i in seq_along(ta[[cnd]])) {
			if ((cnd == fnd && i == ita) || ta[[fnd]][ita] != ta[[cnd]][i] || any(ilc[[cnd]][tsim[[cnd]][,i]] != i_lc)) {
				next # no common layout with the factor
			}
			ib <- bch[[cnd]] == names(ta[[cnd]])[i]
			# see if it has numeric power
			if (!any(sapply(po[[cnd]][ib], is.numeric))) {
				next
			}
			iit <- which(ib) # the bases equal to factor
			p_fa <- min(sapply(po[[cnd]][ib], function(p) if (is.numeric(p)) p else NA), na.rm=TRUE)
			i_lc <- ilc[[cnd]][iit]
			i_nd <- ind[[cnd]][iit]
			fa_nd[[cnd]]$b <- append(fa_nd[[cnd]]$b, lc[[i_lc[1]]][[cnd]]$b[i_nd[1]])
			fa_nd[[cnd]]$p <- append(fa_nd[[cnd]]$p,
p_fa)
			# decrease p in the lc terms
			for (i in seq_along(i_lc)) {
				lc[[i_lc[i]]][[cnd]]$p[[i_nd[i]]] <- Simplify_(call("-", lc[[i_lc[i]]][[cnd]]$p[[i_nd[i]]], p_fa), scache)
			}
		}
	}
#browser()
	# form final symbolic expression
	# replace all i_lc by one product of fa_nd and lincomb of the reduced nds
	rest <- Simplify_(lc2expr(lc[i_lc], scache), scache)
	if (is.neg.expr(rest)) {
		rest <- negate.expr(rest)
		fa_nd$sminus <- !fa_nd$sminus
	}
	fa_nd$num$b <- append(fa_nd$num$b, rest)
	fa_nd$num$p <- append(fa_nd$num$p, 1)
	lc <- c(list(fa_nd), lc[-i_lc])
	return(lc2expr(lc, scache))
}

`Simplify.-` <- function(expr, scache=NULL)
{
	`Simplify.+`(expr, add=FALSE, scache=scache)
}

`Simplify.*` <- function(expr, div=FALSE, scache=NULL)
{
#print(expr)
#browser()
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
		if (sminus) call("-", b) else b
	} else if (b == 1) {
		if (sminus) call("-", a) else a
	} else if (div && identical(a, b)) {
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
		if (div) {
			fa$num <- nd_a$fa$num*nd_b$fa$den
			fa$den <- nd_a$fa$den*nd_b$fa$num
		} else {
			fa$num <- nd_a$fa$num*nd_b$fa$num
			fa$den <- nd_a$fa$den*nd_b$fa$den
		}
		res <- fa$num/fa$den
		if (as.integer(res) == res) {
			fa$num <- res
			fa$den <- 1
		} else if (fa$den != 1) {
			res <- fa$den/fa$num
			if (as.integer(res) == res) {
				fa$num <- 1
				fa$den <- res
			}
		}
		# group identical bases by adding their powers
#browser()
		bch=list()
		for (na in c("num", "den")) {
			bch[[na]] <- sapply(nd[[na]]$b, format1)
			if (length(nd[[na]]$b) <= 1)
				next
			ta <- table(bch[[na]])
			ta <- ta[ta > 1]
			if (length(ta) == 0)
				next
			nd_eq <- outer(bch[[na]], names(ta), `==`)
			for (inum in seq(len=ncol(nd_eq))) {
				isim <- which(nd_eq[,inum])
				if (length(isim)) {
					# add powers for this base
					nd[[na]]$p[[isim[1]]] <- Simplify_(li2sum(nd[[na]]$p[isim]), scache)
					# set grouped powers to 0
					nd[[na]]$p[isim[-1]] <- 0
				}
			}
			# remove power==0 terms
			ize <- isim[-1]
			if (length(ize)) {
				nd[[na]]$b <- nd[[na]]$b[-ize]
				nd[[na]]$p <- nd[[na]]$p[-ize]
				bch[[na]] <- bch[[na]][-ize]
			}
		}
		# simplify identical terms in num and denum by subtracting powers
		nd_eq <- outer(bch$den, bch$num, `==`)
		ipair <- matrix(0, nrow=2, ncol=0)
		for (inum in seq(len=ncol(nd_eq))) {
			iden <- which(nd_eq[,inum]) # of length at most 1 as terms are already grouped
			if (length(iden)) {
				# simplify power for this pair
				ipair <- cbind(ipair, c(inum, iden))
				res <- Simplify_(call("-", nd$num$p[[inum]], nd$den$p[[iden]]), scache)
				if (is.neg.expr(res)) {
					nd$num$p[[inum]] <- 0
					nd$den$p[[iden]] <- negate.expr(res)
				} else {
					nd$num$p[[inum]] <- res
					nd$den$p[[iden]] <- 0
				}
			}
		}
#browser()
		# remove power==0 terms
		for (na in c("num", "den")) {
			if (length(nd[[na]]$b) == 0)
				next
			ize=sapply(nd[[na]]$p, `==`, 0)
			nd[[na]]$b <- nd[[na]]$b[!ize]
			nd[[na]]$p <- nd[[na]]$p[!ize]
		}
		nd[["fa"]] <- fa
		nd[["sminus"]] <- sminus
		expr <- nd2expr(nd, scache)
		expr
	}
}
`Simplify./` <- function(expr, scache=NULL)
{
	`Simplify.*`(expr, div=TRUE, scache=scache)
}
`Simplify.^` <- function(expr, scache=NULL)
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
			b <- Simplify_(call("*", a[[3]], b), scache)
			a <- a[[2]]
		} else if (a[[1]] == as.symbol("sqrt")) {
			# divide by 2
			b <- Simplify_(call("/", b, 2), scache)
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
Simplify.log <- function(expr, scache=NULL) {
	if (is.call(expr[[2]])) {
		# the argument of log is a function
		subf <- as.character(expr[[2]][[1]])
		if (subf == "^") {
			p <- expr[[2]][[3]]
			expr[[2]] <- expr[[2]][[2]]
			expr <- Simplify_(call("*", p, expr), scache)
		} else if (subf == "exp") {
			if (length(expr) == 2)
				expr <- expr[[2]][[2]]
			else
				expr <- Simplify_(call("/", expr[[2]][[2]], call("log", expr[[3]])), scache)
		} else if (subf == "sqrt") {
			expr[[2]] <- expr[[2]][[2]]
			expr <- Simplify_(call("*", 0.5, expr), scache)
		} else if (subf == "*") {
			a <- expr
			a[[2]] <- expr[[2]][[2]]
			expr[[2]] <- expr[[2]][[3]] # unitary "+" cannot appear here
			expr <- Simplify_(call("+", a, expr), scache)
		} else if (subf == "/") {
			a <- expr
			a[[2]] <- expr[[2]][[2]]
			expr[[2]] <- expr[[2]][[3]] # unitary "+" cannot appear here
			expr <- Simplify_(call("-", a, expr), scache)
		} else if (subf == "+") {
			# replace log(1+x) by log1p(x)
			if (expr[[2]][[2]] == 1) {
				expr <- call("log1p", expr[[2]][[3]])
			} else if (expr[[2]][[3]] == 1) {
				expr <- call("log1p", expr[[2]][[2]])
			}
		}
	}
	if (length(expr) == 3 && identical(expr[[2]], expr[[3]])) {
		1
	} else {
		expr
	}
}
Simplify.sqrt <- function(expr, scache=NULL) {
	if (is.call(expr[[2]])) {
		# the argument of sqrt is a function
		subf <- as.character(expr[[2]][[1]])
		if (subf == "^") {
			p <- expr[[2]][[3]]
			Simplify_(call("^",  call("abs", expr[[2]][[2]]), call("/", p, 2)), scache)
		} else if (subf == "exp") {
			expr[[2]][[2]] <- Simplify_(call("/", expr[[2]][[2]], 2), scache)
			expr[[2]]
		} else if (subf == "sqrt") {
			Simplify_(call("^", expr[[2]][[2]], 0.25), scache)
		} else if (subf == "*" && identical(expr[[2]][[2]], expr[[2]][[3]])) {
			Simplify_(call("abs", expr[[2]][[2]]), scache)
		} else {
			expr
		}
	} else {
		expr
	}
}
Simplify.abs <- function(expr, scache=NULL) {
	if (is.uminus(expr[[2]])) {
		expr[[2]] <- expr[[2]][[2]]
	} else if (is.call(expr[[2]])) {
		subf <- as.character(expr[[2]][[1]])
		if (subf == "^") {
			p <- expr[[2]][[3]]
			if (is.numeric(p) && p%%2 == 0)
				expr <- expr[[2]]
		} else if (subf == "exp" || subf == "sqrt") {
			expr <- expr[[2]]
		}
	}
	expr
}
Simplify.sign <- function(expr, scache=NULL) {
	if (is.uminus(expr[[2]])) {
		expr[[2]] <- expr[[2]][[2]]
		expr <- call("-", expr)
	} else if (is.call(expr[[2]])) {
		subf <- as.character(expr[[2]][[1]])
		if (subf == "^") {
			p <- expr[[2]][[3]]
			if (is.numeric(p) && p%%2 == 0)
				expr <- 1
		} else if (subf == "exp" || subf == "sqrt") {
			expr <- 1
		}
	}
	expr
}
Simplify.if <- function(expr, scache=NULL) {
	cond <- expr[[2]]
	if ((is.logical(cond) || is.numeric(cond)) && isTRUE(!!cond)) {
		expr <- expr[[3]]
	} else if (length(expr) == 4) {
		if ((is.logical(cond) || is.numeric(cond)) && isTRUE(!cond)) {
			expr <- expr[[4]]
		} else if (identical(expr[[3]], expr[[4]])) {
			expr <- expr[[3]]
		}
	}
	expr
}
Simplify.bessel <- function(expr, scache=NULL) {
	if (length(expr) < 4)
		return(expr)
	cond <- expr[[4]]
	if ((is.logical(cond) || is.numeric(cond)) && isTRUE(!cond)) {
		expr[[4]] <- NULL
	}
	expr
}
`Simplify.=` <- function(expr, scache=NULL) {
	# just strore the rhs in the scache
	if (is.symbol(expr[[2]]) && is.call(expr[[3]])) {
		scache$l[[format1(expr[[3]])]] <- expr[[2]]
	}
	expr
}
`Simplify.{` <- function(expr, scache=NULL) {
	# if the last expression is a constant just return it
	n <- length(expr)
	la <- expr[[n]]
	if (is.conuloch(la)) {
		expr <- la
	}
	expr
}

Numden <- function(expr) {
	# Return a list with "num" as numerator and "den" as denominator sublists.
	# "fa" field is for numeric factors in "num" and "den" subfields.
	# "sminus" is logical for applying or not "-" to the whole expression
	# Each sublist regroups the language expressions which are not products neither
	# divisions. The terms are decomposed in b^p sublists
	if (is.uminus(expr)) {
		a=Numden(expr[[2]])
		a$sminus <- !a$sminus
		a
	} else if (is.uplus(expr)) {
		Numden(expr[[2]])
	} else if (is.symbol(expr)) {
		list(num=list(b=list(expr), p=list(1)),
			sminus=FALSE,
			fa=list(num=1, den=1))
	} else if (is.numeric(expr)) {
		sminus <- expr < 0
		list(fa=list(num=if (sminus) -expr else expr, den=1),
			sminus=sminus)
	} else if (is.call(expr)) {
		if (expr[[1]] == as.symbol("*")) {
			# recursive call
			a=Numden(expr[[2]])
			b=Numden(expr[[3]])
			list(num=list(b=c(a$num$b, b$num$b), p=c(a$num$p, b$num$p)),
				den=list(b=c(a$den$b, b$den$b), p=c(a$den$p, b$den$p)),
				sminus=xor(a$sminus, b$sminus),
				fa=list(num=a$fa$num*b$fa$num, den=a$fa$den*b$fa$den))
		} else if (expr[[1]] == as.symbol("/")) {
			# recursive call
			a=Numden(expr[[2]])
			b=Numden(expr[[3]])
			list(num=list(b=c(a$num$b, b$den$b), p=c(a$num$p, b$den$p)),
				den=list(b=c(a$den$b, b$num$b), p=c(a$den$p, b$num$p)),
				sminus=xor(a$sminus, b$sminus),
				fa=list(num=a$fa$num*b$fa$den, den=a$fa$den*b$fa$num))
		} else if (expr[[1]] == as.symbol("^")) {
			if (is.neg.expr(expr[[3]])) {
				# make the power look positive
				list(den=list(b=list(expr[[2]]), p=list(negate.expr(expr[[3]]))),
					sminus=FALSE,
					fa=list(num=1, den=1)
				)
			} else {
				list(num=list(b=list(expr[[2]]), p=list(expr[[3]])),
					sminus=FALSE,
					fa=list(num=1, den=1)
				)
			}
		} else {
			list(num=list(b=list(expr), p=list(1)),
				sminus=FALSE,
				fa=list(num=1, den=1))
		}
	} else {
		list(num=list(b=list(expr), p=list(1)),
			sminus=FALSE,
			fa=list(num=1, den=1))
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
is.conuloch <- function(e) {
	# detect if e is complex, numeric, logical or character
	return(is.numeric(e) || is.logical(e) || is.complex(e) || is.character(e))
}
is.neg.expr <- function(e) {
	# detect if e is a negative expression, i.e. is one of:
	#  - negative real number
	#  - unitary minus (-a)
	return((is.numeric(e) && e < 0) || is.uminus(e))
}
negate.expr <- function(e) {
	# make negative expression looking positive or inverse the difference
	if (is.numeric(e)) 
		-e
	else # e is supposed to be a unitary minus
		e[[2]]
}
is.assign <- function(e) {
	# detect if it is an assignment operator
	is.call(e) && (e[[1]] == as.symbol("<-") || e[[1]] == as.symbol("="))
}
is.subindex <- function(e) {
	# is e a simple subindex expression?
	is.call(e) && any(as.character(e[[1]]) == c("$", "[", "[[")) && (is.symbol(e[[2]]) && (is.symbol(e[[3]]) || is.conuloch(e[[3]])))
}
Lincomb <- function(expr) {
	# decompose expr in a list of product terms (cf Numden)
	# the sign of each term is determined by the nd$sminus logical item.
	if (is.call(expr) && length(expr) == 3) {
		if (expr[[1]] == as.symbol("+")) {
			# recursive call
			c(Lincomb(expr[[2]]), Lincomb(expr[[3]]))
		} else if (expr[[1]] == as.symbol("-")) {
			# recursive call
			a <- Lincomb(expr[[2]])
			b <- Lincomb(expr[[3]])
			# inverse the sign in b terms
			b <- lapply(b, function(it) {it$sminus <- !it$sminus; it})
			c(a, b)
		} else {
			list(Numden(expr))
		}
	} else {
		list(Numden(expr))
	}
}

# return an environement in wich stored subexpressions with
# an index giving the position of each subexpression in the
# whole statement st ("rhs" entry). Index is given as a string i1.i2.i3...
# where the integeres iN refer to st[[i2]][[i3]][[...]]
# "lhs" is index to char mapping (what is defined where)
# "def" is a mapping of lhs (char) to rhs (char)
# "{" where accolade operators are
Leaves <- function(st, ind="1", res=new.env()) {
	if (is.null(res$rhs)) {
		res$rhs <- list()
		res$lhs <- list()
		res$def <- list() # store definitions by asignments
		res[["{"]] <- list()
	}
	if (is.call(st)) {
		if (st[[1]] != as.symbol("<-") && st[[1]] != as.symbol("=")) {
			res$rhs[[ind]] <- format1(st)
			if (st[[1]] == as.symbol("{")) {
				res[["{"]] <- append(res[["{"]], ind)
			}
		} else {
			if (!is.null(res$lhs[[ind]]))
				stop("Re-assignment is not supported yet in caching.")
			if (is.call(st[[2]]))
				stop("Cannot handle yet indexing in left values.")
			lhs <- as.character(st[[2]])
			res$lhs[[ind]] <- lhs # we cannot handle yet `[`, `$` etc.
			res$def[[lhs]] <- format1(st[[3]])
			# exclude this assignement from replacements if .eX
			#if (regexpr("\\.+e[0-9]+", lhs) > 0)
			#	return(res)
		}
		args <- as.list(st)[-1]
		l <- lapply(seq_along(args), function(i) Leaves(args[[i]], paste(ind, i+1, sep="."), res))
	}
	return(res)
}

# convert index calculated by Leaves() to a call like st[[i2]][[i3]]...
# the first two chars "1." are striped out
ind2call <- function(ind, st="st")
	if (ind == "1") as.symbol(st) else parse(text=sprintf("%s[[%s]]", st, gsub("\\.", "]][[", substring(ind, 3))))[[1]]

# replace repeated subexpressions by cached values
# prefix is used to form auxiliary variable
##' @rdname Simplify
Cache <- function(st, env=Leaves(st), prefix="") {
	stch <- if (is.call(st)) as.character(st[[1]]) else ""
	env$lhs <- unlist(env$lhs)
	#if (stch == "<-" || stch == "=") {
	#	return(call("<-", st[[2]], Cache(st[[3]], env=env, prefix=paste(".", st[[2]], sep=""))))
	#} else if (stch == "{" || stch == "c") {
	#	return(as.call(c(list(st[[1]]), lapply(as.list(st)[-1], Cache, env=env))))
	#}
	alva <- all.vars(st)
	p <- grep(sprintf("^%s.e[0-9]+", prefix), alva, value=T)
	if (nchar(prefix) == 0 && length(p) > 0) {
		prefix <- max(p)
	}
	ve <- unlist(env$rhs)
	defs <- unlist(env$def)
	tdef <- outer(ve, defs, "==")
#browser()
	# if the subexpression is in defs, replace it with the symbol in the lhs
	for (ic in seq_len(ncol(tdef))) {
		v <- tdef[,ic]
		nme <- colnames(tdef)[ic]
		idef <- names(which(env$lhs==nme))
		for (i in which(v)) {
			ind <- names(v)[i] # subexpression index in st
			# skip self assignment
			ispl <- strsplit(ind, ".", fixed=TRUE)[[1]]
			indup <- paste(ispl[-length(ispl)], collapse=".")
			stup <- eval(ind2call(indup))
			if ((is.assign(stup) && (as.character(stup[[2]]) == nme || natcompare(indup, idef) < 0)))
				next
			ve[i] <- NA
			do.call(`<-`, list(ind2call(ind), quote(as.symbol(nme))))
		}
	}

	suppressWarnings(ve <- ve[!is.na(ve)])
	# skip simple subindex
	isi <- sapply(ve, function(e) is.subindex(parse(text=e)[[1]]))
	ve <- ve[!isi]
	
	ta <- table(ve)
	ta <- ta[ta > 1]
	if (length(ta) == 0)
		return(st)
	e <- list() # will store the result code
	alva <- list()
	for (sub in names(sort(ta, decreasing=TRUE))) {
		# get st indexes for this subexpression
		isubs <- names(which(ve == sub))
		for (i in seq_along(isubs)) {
			isub <- isubs[i]
			subst <- ind2call(isub)
			if (i == 1) {
				esubst <- try(eval(subst), silent=TRUE)
				if (inherits(esubst, "try-error"))
					break # was already cached
				# add subexpression to the final code
				ie=length(e)+1
				estr <- sprintf("%s.e%d", prefix, ie)
				esub <- as.symbol(estr)
				e[[ie]] <- call("<-", esub, esubst)
				alva[[estr]] <- all.vars(esubst)
			}
			# replace subexpression in st by .eX
			do.call(`<-`, list(subst, as.symbol("esub")))
		}
	}
	alva[["end"]] <- all.vars(st)
	# where .eX are used? If only once, develop, replace and remove it
	wh <- lapply(seq_along(e), function(i) {
		it=sprintf("%s.e%d", prefix, i)
		which(sapply(alva, function(v) any(it == v)))
	})
	dere <- sapply(wh, function(it) if (length(it) == 1 && names(it) != "end") it[[1]] else 0)
	for (i in which(dere != 0)) {
		idest <- dere[i]
		li <- list()
		li[[sprintf("%s.e%d", prefix, i)]] <- e[[i]][[3]]
		e[[idest]][[3]] <- do.call("substitute", c(e[[idest]][[3]], list(li)))
	}
	e <- e[which(!dere)]
#browser()
	# place auxiliary vars after the definition of the used vars
	if (stch != "{") {
		l <- c(as.symbol("{"), e, st)
		st <- as.call(l)
	} else {
		n <- length(st)
		res <- c(e, as.list(st)[-c(1, n)])
		alva <- lapply(res, all.vars)
		i <- toporder(alva)
		res <- c(as.symbol("{"), res[i], st[[n]])
		st <- as.call(res)
	}
	return(st)
}

##' @rdname Simplify
deCache <- function(st) {
	# do the job inverse to Cache(), i.e. substitute all auxiliary expressions
	# in the final one
	# NB side effect: all assignement not used in the last operation in {...} are
	# just lost.
	if (!is.call(st)) {
		return(st)
	}
	stch <- as.character(st[[1]])
	stl <- as.list(st)
	if (stch == "{") {
		repl <- list()
		for (op in stl[-1]) {
			# gather substitutions
			if (is.assign(op)) {
				repl[[as.character(op[[2]])]] <- do.call("substitute", list(deCache(op[[3]]), repl))
			}
		}
		# the last operation subst
		la <- stl[[length(stl)]]
		if (is.assign(la)) {
			st <- repl[[length(repl)]]
		} else {
			st <- do.call("substitute", list(deCache(la), repl))
		}
	} else {
		# recurrsive call to deCache on all arguments of the statement
		stl <- lapply(stl, deCache)
		st <- as.call(stl)
	}
	return(st)
}

nd2expr <- function(nd, scache, sminus=NULL) {
	# form symbolic products
	# if sminus is not null, use it instead of the nd's one
	if (length(nd) == 0)
		return(0)
	eprod <- list()
	sminus <- (!is.null(sminus) && sminus) || (is.null(sminus) && nd$sminus)
	for (na in c("num", "den")) {
		if (length(nd[[na]]$b) == 0)
			next
		# alphabetic order for bases, symbols first, then calls
		for (i in order(sapply(nd[[na]]$b, is.call), sapply(nd[[na]]$b, format1))) {
			p <- nd[[na]]$p[[i]]
			if (p == 0)
				next
			term <- if (p == 1) nd[[na]]$b[[i]] else Simplify_(call("^", nd[[na]]$b[[i]], p), scache)
			if (term == 0)
				return(if (na == "num") 0 else if (sminus) -Inf else Inf)
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
	fa=nd$fa
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
		else if (expr == 1)
			expr <- fa$num
		else
			expr <- call("*", fa$num, expr)
	} else if (fa$den != 1) {
		if (is.call(expr) && expr[[1]] == as.symbol("/"))
			expr[[3]] <- call("*", fa$den, expr[[3]])
		else
			expr <- call("/", expr, fa$den)
	}
	expr <- if (sminus) call("-", expr) else expr
#print(sprintf("nd->%s", format1(expr)))
	return(expr)
}

lc2expr <- function(lc, scache) {
	# form symbolic sum and diff form a list of nds
	# separate in positive and negative
	smin <- sapply(lc, "[[", "sminus")
	epos <- lapply(lc[which(!smin)], nd2expr, scache)
	if (length(epos) > 1) {
#cat("epos orig=", sapply(epos, format1), sep="\n")
#cat("epos order=", order(sapply(epos, format1)), sep="\n")
#cat("order - +=", order(c("-", "+")), sep="\n")
		epos <- epos[order(sapply(epos, format1), decreasing = FALSE)]
#cat("epos=", sapply(epos, format1), sep="\n")
	}
	eneg <- lapply(lc[which(smin)], nd2expr, scache, sminus=FALSE)
	if (length(eneg) > 1) {
		eneg <- eneg[order(sapply(eneg, format1))]
	}
	if (length(epos) == 0)
		return(if (length(eneg) == 0) 0 else call("-", li2sum(eneg)))
	else
		return(if (length(eneg) == 0) li2sum(epos) else Simplify_(call("-", li2sum(epos), li2sum(eneg)), scache))
}

li2sum <- function(li) {
	# form a long sum of expressions from the list li
	len <- length(li)
	if (len == 0)
		0
	else if (len == 1)
		li[[1]]
	else if (len == 2)
		if (li[[1]] == 0)
			li[[2]]
		else if (li[[2]] == 0)
			li[[1]]
		else
			call("+", li[[1]], li[[2]])
	else
		call("+", li2sum(li[-len]), li[[len]])
}
toporder <- function(l, ind=seq_along(l), vars=sapply(l, `[[`, 1)) {
	# Topological ordering of assignement operators
	# l is a list whose memebers are resulted from all.vars(op)
	# ind is a subindexing vector for l (for recursive call)
	# vars is a vector of variable which are assigned in ops[ind]
	# return a vector of indexes like in order()
	
	# find independet assignements, i.e. whose rhs vars are not in vars
#cat("ind=", ind, "\n")
	if (length(ind) <= 1) {
		return(ind)
	}
	rhsvar <- lapply(l[ind], `[`, -1)
	indep <- which(!sapply(rhsvar, function(v) any(v %in% vars)))
#cat("indep=", ind[indep], "\n")
	return(c(ind[indep], toporder(l, ind[-indep], vars[-indep])))
}
natcompare <- function(s1, s2, sep="[^0-9]+") {
	# Compare two strings in natural ordering,
	# i.e. natlower("1.12", "1.2") returns 1 (i.e s1 is greater than s2)
	# while plain "1.12" < "1.2" returns TRUE
	# sep is separator for string splitting
	# By default any non number chain of characters
	# is used as a single separator and thus is exlculed
	# from comparison.
	# The fields after string splitting are compared as numerics
	# Empty string or NA are considered as -Inf, i.e. they are less
	# than any other finite number.
	# Return -1 if s1 is lower s2, 0 if s1 is equal to s2 and 1 otherwise
	# 
	v1 <- as.numeric(strsplit(s1, sep[1])[[1]])
	v1[is.na(v1)] <- -Inf
	v2 <- as.numeric(strsplit(s2, sep[1])[[1]])
	v2[is.na(v2)] <- -Inf
	l1 <- length(v1)
	l2 <- length(v2)
	lmin <- min(l1, l2)
	# complete the shortest vector by -Inf
	v1 <- c(v1, rep(-Inf, l2-lmin))
	v2 <- c(v2, rep(-Inf, l1-lmin))
	m1 <- v1 < v2
	eq <- v1 == v2
	p1 <- v1 > v2
	if (all(m1) || (any(m1) && all(!p1)) || any(head(which(m1), 1) < head(which(p1), 1))) {
		-1
	} else if (all(eq)) {
		0
	} else {
		1
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
assign("sign", `Simplify.sign`, envir=simplifications)
assign("if", `Simplify.if`, envir=simplifications)
assign("besselI", `Simplify.bessel`, envir=simplifications)
assign("besselK", `Simplify.bessel`, envir=simplifications)
#assign("<-", `Simplify.=`, envir=simplifications)
#assign("=", `Simplify.=`, envir=simplifications)
assign("{", `Simplify.{`, envir=simplifications)
