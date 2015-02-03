# Simplify.R -- symbolic simplification
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
	} else if (format1(a) == format1(b) && !add) {
		return(0)
	}
	# group and simplify identical terms
	apn <- Sumdiff(a)
	bpn <- Sumdiff(b)
	if (add)
		pn <- list(pos=c(apn$pos, bpn$pos), neg=c(apn$neg, bpn$neg))
	else
		pn <- list(pos=c(apn$pos, bpn$neg), neg=c(apn$neg, bpn$pos))
	# group all numerics
	nu <- list()
	for (na in c("pos", "neg")) {
		inum <- sapply(pn[[na]], is.numeric)
		nu[[na]] <- sum(unlist(pn[[na]][inum]))
		# remove numerics from main term lists
		pn[[na]] <- pn[[na]][!inum]
	}
	nu <- nu$pos-nu$neg
	if (nu > 0) {
		pn$pos <- c(nu, pn$pos)
	} else if (nu < 0) {
		pn$neg <- c(-nu, pn$neg)
	}
	# group identical terms
	ta <- list()
	pnch <- list()
	for (na in c("pos", "neg")) {
		if (length(pn[[na]]) == 0)
			next
		pnch[[na]] <- sapply(pn[[na]], format1)
		ta[[na]] <- table(pnch[[na]]) # count of unique entries
		isim <- apply(outer(pnch[[na]], names(ta[[na]]), `==`), 2, function(v) which(v)[1])
		# keep only unique terms
		pn[[na]] <- pn[[na]][isim]
		pnch[[na]] <- pnch[[na]][isim]
	}
	# simplify ta for identical terms in pos and neg
	ijsim <- outer(pnch$neg, pnch$pos, `==`) # each column has at most one TRUE
	irm=c()
	for (ipos in seq_along(pnch$pos)) {
		ineg <- which(ijsim[,ipos])
		if (length(ineg) == 0)
			next
		ta$pos[ipos] <- ta$pos[ipos] - ta$neg[ineg]
		irm <- c(irm, ineg)
	}
	# remove simplified terms from neg
	if (length(irm)) {
		pn$neg <- pn$neg[-irm]
		ta$neg <- ta$neg[-irm]
	}
	# move negative coefs from pos to neg
	ineg <- which(ta$pos < 0)
	if (length(ineg)) {
		pn$neg <- c(pn$neg, pn$pos[ineg])
		pn$pos <- pn$pos[-ineg]
		ta$neg <- c(ta$neg, -ta$pos)
		ta$pos <- ta$pos[-ineg]
		pnch$neg <- c(pnch$neg, pnch$pos[ineg])
		pnch$pos <- pnch$pos[-ineg]
	}
	# remove ta==0 in pos
	inul <- ta$pos == 0
	pn$pos <- pn$pos[!inul]
	ta$pos <- ta$pos[!inul]
	pnch$pos <- pnch$pos[!inul]
	
	for (na in c("pos", "neg")) {
		# where ta > 1, replace term by nb_repeat*term
		pn[[na]][ta[[na]] > 1] <- lapply(which(ta[[na]] > 1), function(i) Simplify_(call("*", ta[[na]][i], pn[[na]][i])))
	}
	# form final symbolic expression
	if (length(pn$pos) == 0 && length(pn$neg) == 0) {
			return(0)
	}
	res=list()
	for (na in c("pos", "neg")) {
		if (length(pn[[na]]) == 0)
			next
		iord <- order(pnch[[na]])
		expr <- pn[[na]][[iord[1]]]
		for (i in iord[-1]) {
			expr <- call("+", expr, pn[[na]][[i]])
		}
		res[[na]] <- expr
	}
	if (is.null(res$pos))
		return(call("-", res$neg))
	else if (is.null(res$neg))
		return(res$pos)
	else
		return(call("-", res$pos, res$neg))
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
			if (length(nd[[na]]$b) == 0 && fa[[na]] == 1)
				next
			# remove power==0 terms
			ize=sapply(nd[[na]]$p, `==`, 0)
			nd[[na]]$b <- nd[[na]]$b[!ize]
			nd[[na]]$p <- nd[[na]]$p[!ize]
			if (length(nd[[na]]$b) == 0 && fa[[na]] == 1)
				next
			eprod[[na]] <- if (fa[[na]] != 1) fa[[na]] else NULL
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
	} else if (expr[[1]] == as.symbol("*")) {
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
	} else if (is.call(expr) && expr[[1]] == as.symbol("^")) {
		if (expr[[3]] < 0) {
			# make the power look positive
			list(den=list(b=expr[[2]], p=if (is.numeric(expr[[3]])) -expr[[3]] else expr[[3]][[2]]), sminus=FALSE)
		} else {
			list(num=list(b=expr[[2]], p=expr[[3]]), sminus=FALSE)
		}
	} else {
		list(num=list(b=expr, p=1), sminus=FALSE)
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
Sumdiff <- function(expr) {
	# Return a list with "pos" as a list of positive terms in a sum
	# and "neg" as a list of negative terms in a sum
	if (is.uminus(expr)) {
		a <- Sumdiff(expr[[2]])
		list(pos=a$neg, neg=a$pos)
	} else if (is.uplus(expr)) {
		Sumdiff(expr[[2]])
	} else if (is.symbol(expr) || is.numeric(expr)) {
		if (expr > 0) list(pos=expr) else list(neg=-expr)
	} else if (is.call(expr) && expr[[1]] == as.symbol("+")) {
		# recursive call
		a <- Sumdiff(expr[[2]])
		b <- Sumdiff(expr[[3]])
		list(pos=c(a$pos, b$pos), neg=c(a$neg, b$neg))
	} else if (is.call(expr) && expr[[1]] == as.symbol("-")) {
		# recursive call
		a <- Sumdiff(expr[[2]])
		b <- Sumdiff(expr[[3]])
		list(pos=c(a$pos, b$neg), neg=c(a$neg, b$pos))
	} else {
		list(pos=list(expr))
	}
}
Simplify.log <- function(expr) {
	if (is.call(expr[[2]])) {
		if (expr[[2]][[1]] == as.symbol("^")) {
			p <- expr[[2]][[3]]
			expr[[2]] <- expr[[2]][[2]]
			Simplify_(call("*", p, expr))
		} else if (expr[[2]][[1]] == as.symbol("exp")) {
			if (length(expr) == 2)
				expr[[2]][[2]]
			else
				Simplify_(call("/", expr[[2]][[2]], call("log", expr[[3]])))
		} else if (expr[[2]][[1]] == as.symbol("*")) {
			a <- expr
			a[[2]] <- expr[[2]][[2]]
			expr[[2]] <- expr[[2]][[3]] # unitary "+" cannot appear here
			Simplify_(call("+", Simplify_(a), Simplify_(expr)))
		} else if (expr[[2]][[1]] == as.symbol("/")) {
			a <- expr
			a[[2]] <- expr[[2]][[2]]
			expr[[2]] <- expr[[2]][[3]] # unitary "+" cannot appear here
			Simplify_(call("-", a, expr))
		} else {
			expr
		}
	} else if (length(expr) == 3 && format1(expr[[2]]) == format1(expr[[3]])) {
		1
	} else {
		expr
	}
}
