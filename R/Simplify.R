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
#'    \item an expression: \code{expression(x+x)}
#'    \item a string: \code{"x+x"}
#'    \item a function: \code{function(x) x+x}
#'    \item a right hand side of a formula: \code{~x+x}
#'    \item a language: \code{quote(x+x)}
#' }
#' @param env An environment in wich a simplified function is created
#'  if \code{expr} is a function. This argument is ignored in all other cases.
#' @return A simplified expression. The result is of the same type as
#'  \code{expr} except for formula, where a language is returned.
#' @details An environment \code{simplifications} containing simplification rules, is exported in the user namespace.
Simplify <- function(expr, env=parent.frame()) {
	te <- try(expr, silent=TRUE)
	if (inherits(te, "try-error")) {
		expr <- substitute(expr)
	}
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
	if (is.call(expr)) {
		args <- lapply(as.list(expr)[-1], Simplify_)
		expr[-1]=args
		if (all(sapply(args, is.numeric))) {
			# if all arguments are numeric, evaluate them
			return(eval(expr))
		} else {
			# is there a rule in the table?
			sym.name <- as.character(expr[[1]])
			Simplify.rule <- simplifications[[sym.name]]
			if (!is.null(Simplify.rule)) {
				return(Simplify.rule(expr))
			} else {
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
	
	if (a == 0) {
		return(if (add) b else call("-", b))
	} else if (b == 0) {
		return(a)
	} else if (add && is.uminus(a) && !is.uminus(b)) {
		a <- b
		b <- expr[[2]][[2]]
		expr <- call("-", a, b)
	} else if (format1(a) == format1(b)) {
		return(if (add) call("*", 2, a) else 0)
	} else if (!is.call(a) && !is.call(b)) {
		return(expr) # nothing to simplify
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
		if (any(sapply(po[[cnd]], is.numeric))) {
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
		return(expr) # nothing to factorize
	ond <- if (fnd == "num") "den" else "num"
	# create nd with the first factor
	fa_nd <- list(num=list(b=list(), p=list()),
		den=list(b=list(), p=list()),
		sminus=FALSE, fa=list(num=1, den=1))
	fa_nd[[fnd]]$b <- lc[[i_lc[1]]][[fnd]]$b[i_nd[1]]
	fa_nd[[fnd]]$p <- list(p_fa)
	# decrease p in the lc terms
	for (i in seq_along(i_lc)) {
		lc[[i_lc[i]]][[fnd]]$p[[i_nd[i]]] <- Simplify_(call("-", lc[[i_lc[i]]][[fnd]]$p[[i_nd[i]]], p_fa))
	}
	
	for (cnd in c(fnd, ond)) {
		# see if other side can provide factors
		for (i in seq_along(ta[[cnd]])) {
			if ((cnd == fnd && i == ita) || ta[[fnd]][ita] != ta[[cnd]][i] || any(ilc[[cnd]][tsim[[cnd]][,i]] != i_lc)) {
				next # no common layout with factor
			}
			ib <- bch[[cnd]] == names(ta[[cnd]])[i]
			# see if it has numeric power
			if (!any(sapply(po[[cnd]], is.numeric))) {
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
				lc[[i_lc[i]]][[cnd]]$p[[i_nd[i]]] <- Simplify_(call("-", lc[[i_lc[i]]][[cnd]]$p[[i_nd[i]]], p_fa))
			}
		}
	}
#browser()
	# form final symbolic expression
	# replace all i_lc by one product of fa_nd and lincomb of the reduced nds
	fa_nd$num$b <- append(fa_nd$num$b, lc2expr(lc[i_lc]))
	fa_nd$num$p <- append(fa_nd$num$p, 1)
	lc <- c(list(fa_nd), lc[-i_lc])
	return(lc2expr(lc))
}

`Simplify.-` <- function(expr)
{
	`Simplify.+`(expr, add=FALSE)
}

`Simplify.*` <- function(expr, div=FALSE)
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
				res <- Simplify_(call("-", nd$num$p[[inum]], nd$den$p[[iden]]))
				if (res > 0) {
					nd$num$p[[inum]] <- res
					nd$den$p[[iden]] <- 0
				} else {
					nd$num$p[[inum]] <- 0
					nd$den$p[[iden]] <- Simplify_(substitute(-res))
				}
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
		expr <- nd2expr(nd)
		expr
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
Simplify.sign <- function(expr) {
	if (is.uminus(expr[[2]])) {
		expr[[2]] <- expr[[2]][[2]]
		expr <- call("-", expr)
	} else if (is.call(expr[[2]])) {
		if (expr[[2]][[1]] == as.symbol("^")) {
			p <- expr[[2]][[3]]
			if (is.numeric(p) && p%%2 == 0)
				expr <- 1
		} else if (expr[[2]][[1]] == as.symbol("exp") || expr[[2]][[1]] == as.symbol("sqrt")) {
			expr <- 1
		}
	}
	expr
}
Simplify.if <- function(expr) {
	cond <- expr[[2]]
	if ((is.logical(cond) || is.numeric(cond)) && isTRUE(!!cond)) {
		expr <- expr[[3]]
	} else if (length(expr) == 4) {
		if ((is.logical(cond) || is.numeric(cond)) && isTRUE(!cond)) {
			expr <- expr[[4]]
		} else if (format1(expr[[3]]) == format1(expr[[4]])) {
			expr <- expr[[3]]
		}
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
		list(num=list(b=list(expr), p=1),
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
			if (expr[[3]] < 0) {
				# make the power look positive
				list(den=list(b=list(expr[[2]]), p=if (is.numeric(expr[[3]])) -expr[[3]] else expr[[3]][[2]]),
					sminus=FALSE,
					fa=list(num=1, den=1))
			} else {
				list(num=list(b=list(expr[[2]]), p=expr[[3]]),
					sminus=FALSE,
					fa=list(num=1, den=1))
			}
		} else {
			list(num=list(b=list(expr), p=1),
				sminus=FALSE,
				fa=list(num=1, den=1))
		}
	} else {
		list(num=list(b=list(expr), p=1),
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
# whole statement st
Leaves <- function(st, ind="1", res=new.env()) {
	if (is.call(st)) {
		res[[ind]] <- format1(st)
		args <- as.list(st)[-1]
		l <- lapply(seq_along(args), function(i) Leaves(args[[i]], paste(ind, i+1, sep="."), res))
	}
	return(res)
}

# replace repeated subexpressions by cached values
Cache <- function(st, env=Leaves(st), prefix="") {
	stch <- if (is.call(st)) as.character(st[[1]]) else ""
	if (stch == "<-" || stch == "=") {
		return(call("<-", st[[2]], Cache(st[[3]], prefix=paste(".", st[[1]], sep=""))))
	} else if (stch == "{") {
		return(as.call(c(list(st[[1]]), lapply(as.list(st)[-1], Cache))))
	}
	alva <- all.vars(st)
	p <- grep(sprintf("^%s.e[0-9]+", prefix), alva, value=T)
	if (length(p) > 0) {
		prefix <- max(p)
	}
	ve <- unlist(as.list(env))
	ta <- table(ve)
	ta <- ta[ta > 1]
	if (length(ta) == 0)
		return(st)
	e <- list() # will store the result code
	alva <- list()
	for (sub in names(sort(ta, decreasing=TRUE))) {
		# get indexes for this subexpression
		isubs <- names(which(ve == sub))
		for (i in seq_along(isubs)) {
			isub <- isubs[i]
			subst <- parse(text=sprintf("st[[%s]]", gsub("\\.", "]][[", substring(isub, 3))))[[1]]
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
#browser()
	alva[["end"]] <- all.vars(st)
	# where .eX are used? If only once, develop, replace and remove it
	wh <- lapply(seq_along(e), function(i) {
		it=sprintf("%s.e%d", prefix, i)
		which(sapply(alva, function(v) any(it == v)))
	})
	# the final touch
	e[[ie+1]] <- st
	dere <- sapply(wh, function(it) if (length(it) == 1 && names(it) != "end") it[[1]] else 0)
	for (i in which(dere != 0)) {
		idest <- dere[i]
		li <- list()
		li[[sprintf("%s.e%d", prefix, i)]] <- e[[i]][[3]]
		e[[idest]][[3]] <- do.call("substitute", c(e[[idest]][[3]], list(li)))
	}
	e <- c(list(as.symbol("{")), e[which(!dere)], e[length(e)])
	return(as.call(e))
}
nd2expr <- function(nd, sminus=NULL) {
	# form symbolic products
	# if sminus is not null, use it instead of the nd's one
	if (length(nd) == 0)
		return(0)
	eprod <- list()
	for (na in c("num", "den")) {
		if (length(nd[[na]]$b) == 0)
			next
		# alphabetic order for bases, symbols first, then calls
		for (i in order(sapply(nd[[na]]$b, is.call), sapply(nd[[na]]$b, format1))) {
			p <- nd[[na]]$p[[i]]
			if (p == 0)
				next
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
		else
			expr <- call("*", fa$num, expr)
	} else if (fa$den != 1) {
		if (is.call(expr) && expr[[1]] == as.symbol("/"))
			expr[[3]] <- call("*", fa$den, expr[[3]])
		else
			expr <- call("/", expr, fa$den)
	}
	expr <- if ((!is.null(sminus) && sminus) || (is.null(sminus) && nd$sminus)) substitute(-expr) else expr
#print(sprintf("nd->%s", format1(expr)))
	return(expr)
}
lc2expr <- function(lc) {
	# form symbolic sum and diff form a list of nds
	# separate in positive and negative
	smin <- sapply(lc, "[[", "sminus")
	epos <- lapply(lc[which(!smin)], nd2expr)
	eneg <- lapply(lc[which(smin)], nd2expr, sminus=FALSE)
	if (length(epos) == 0)
		return(if (length(eneg) == 0) 0 else substitute(-li2sum(eneg)))
	else
		return(if (length(eneg) == 0) li2sum(epos) else call("-", li2sum(epos), li2sum(eneg)))
}
li2sum <- function(li) {
	# form a long sum of expressions from the list li
	if (length(li) == 0)
		0
	else if (length(li) == 1)
		li[[1]]
	else if (length(li) == 2)
		call("+", li[[1]], li[[2]])
	else
		call("+", li[[1]], li2sum(li[-1]))
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
