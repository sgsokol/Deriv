#' @name Deriv
#' @title Symbollic differentiation of an expression or function
#' @description Symbollic differentiation of an expression or function
#' @aliases Deriv drule
#' @concept symbollic differentiation
#' 
#' @param f An expression or function to be differentiated.
#'  f can be \itemize{
#'   \item a user defined function: \code{function(x) x**n}
#'   \item a string: \code{"x**n"}
#'   \item an expression: \code{expression(x**n)}
#'   \item a call: \code{call("^", quote(x), quote(n))}
#'   \item a language: \code{quote(x**n)}
#'   \item a right hand side of a formula: \code{~ x**n} or \code{y ~ x**n}
#'  }
#' @param x An optional character vector with variable name(s) with resptect to which
#'  \code{f} must be differentiated. If not provided (i.e. x=NULL), x is
#'  guessed either from\ code{names(formals(f))} (if \code{f} is a function)
#'  or from all variables in f in other cases.
#'  To differentiate expressions including components of lists or vectors, i.e. by expressions like
#'  \code{p[1]}, \code{theta[["alpha"]]} or \code{theta$beta}, the vector of
#'  variables \code{x}
#'  must be a named vector. For the cited examples, \code{x} must be given
#'  as follows \code{c(p="1", theta="alpha", theta="beta")}. Note the repeated name \code{theta} which must be provided for every component of the list \code{theta} by which a
#'  differerentiation is required.
#' @param env An environment where the symbols and functions are searched for.
#'  Defaults to \code{parent.frame()} for \code{f} expression and to
#'  \code{environment(f)} if \code{f} is a function. For primitive function,
#'  it is set by default to .GlobalEnv
#' @param use.D An optional logical (default FALSE), indicates if base::D()
#'  must be used for differentiation of basic expressions.
#' @param cache.exp An optional logical (default TRUE), indicates if
#'  final expression must be optimized with cached subexpressions.
#'  If enabled, repeated calculations are made only once and their
#'  results stored in cache variables which are then reused.
#' @param nderiv An optional integer vector of derivative orders to calculate.
#'  Default NULL value correspond to one differentiation. If length(nderiv)>1,
#'  the resulting expression is a list where each component corresponds to derivative order
#'  given in nderiv. Value 0 corresponds to the original function or expression  non
#'  differentiated. All values must be non negative. If the entries in nderiv
#'  are named, their names are used as names in the returned list. Otherwise
#'  the value of nderiv component is used as a name in the resulting list.
#' 
#' @return \itemize{
#'  \item a function if \code{f} is a function
#'  \item an expression if \code{f} is an expression
#'  \item a character string if \code{f} is a character string
#'  \item a language (usually a so called 'call' but may be also a symbol or just a numeric) for other types of \code{f}
#' }
#'
#' @details
#' R already contains two differentiation functions: D and deriv. D does
#' simple univariate differentiation.  "deriv" uses D do to multivariate
#' differentiation.  The output of "D" is an expression, whereas the output of
#' "deriv" can be an executable function.
#' 
#' R's existing functions have several limitations.  They can probably be fixed,
#' but since they are written in C, this would probably require a lot of work.
#' Limitations include:
#' \itemize{
#'  \item The derivatives table can't be modified at runtime, and is only available
#' in C.
#'  \item Function cannot substitute function calls.  eg:
#'	f <- function(x, y) x + y; deriv(~f(x, x^2), "x")
#' }
#'
#' So, here are the advantages of this implementation:
#' 
#' \itemize{
#'  \item It is entirely written in R, so would be easier to maintain.
#'  \item Can do multi-variate differentiation.
#'  \item Can differentiate function calls:
#'  \itemize{
#'	   \item if the function is in the derivative table, then the chain rule
#'	is applied.  For example, if you declared that the derivative of
#'	sin is cos, then it would figure out how to call cos correctly.
#'	   \item if the function is not in the derivative table (or it is anonymous),
#'	then the function body is substituted in.
#'	   \item these two methods can be mixed.  An entry in the derivative table
#'	need not be self-contained -- you don't need to provide an infinite
#'	chain of derivatives.
#'  }
#'  \item It's easy to add custom entries to the derivatives table, e.g.
#'   
#'   \code{drule[["cos"]] <- alist(x=-sin(x))}
#'   
#'   The chain rule will be automatically applied if needed.
#'  \item The output is an executable function, which makes it suitable
#'      for use in optimization problems.
#'  \item Compound functions (i.e. piece-wise functions based on if-else operator) can
#'      be differentiated (cf. examples section).
#'  \item in case of multiple derivatives (e.g. gradient and hessian calculation),
#'      caching can make calculation economies for both
#' }
#' 
#' Two work environments \code{drule} and \code{simplifications} are
#' exported in the package namescape.
#' As their names indicate, they contain tables of derivative and
#' simplification rules.
#' To see the list of defined rules do \code{ls(drule)}.
#' To add your own derivative rule for a function called say \code{sinpi(x)} calculating sin(pi*x), do \code{drule[["sinpi"]] <- alist(x=pi*cospi(x))}.
#' Here, "x" stands for the first and unique argument in \code{sinpi()} definition. For a function that might have more than one argument,
#' e.g. \code{log(x, base=exp(1))}, the drule entry must be a list with a named rule
#' per argument. See \code{drule$log} for an example to follow.
#' After adding \code{sinpi} you can differentiate expressions like \code{Deriv(~ sinpi(x^2), "x")}. The chain rule will automatically apply.
#' 
#' NB. In \code{abs()} and \code{sign()} function, singularity treatment at point 0 is left to user's care.
#' 
#' NB2. In Bessel functions, derivatives are calculated only by the first argument,
#'      not by the \code{nu} argument which is supposed to be constant.
#' 
#' @author Andrew Clausen (original version) and Serguei Sokol (maintainer)
#' @examples
#'
#' \dontrun{f <- function(x) x^2}
#' \dontrun{Deriv(f)}
#' # function (x)
#' # 2 * x
#' 
#' \dontrun{f <- function(x, y) sin(x) * cos(y)}
#' \dontrun{Deriv(f)}
#' # function (x, y)
#' # c(x = cos(x) * cos(y), y = -(sin(x) * sin(y)))
#'
#' \dontrun{f_ <- Deriv(f)}
#' \dontrun{f_(3, 4)}
#' #              x         y
#' # [1,] 0.6471023 0.1068000
#' 
#' \dontrun{Deriv(~ f(x, y^2), "y")}
#' # -(2 * (y * sin(x) * sin(y^2)))
#' 
#' \dontrun{Deriv(quote(f(x, y^2)), c("x", "y"), cache.exp=FALSE)}
#' # c(x = cos(x) * cos(y^2), y = -(2 * (y * sin(x) * sin(y^2))))
#' 
#' \dontrun{Deriv(expression(sin(x^2) * y), "x")}
#' # expression(2*(x*y*cos(x^2)))
#' 
#' Deriv("sin(x^2) * y", "x") # differentiate only by x
#' "2 * (x * y * cos(x^2))"
#' 
#' Deriv("sin(x^2) * y", cache.exp=FALSE) # differentiate by all variables (here by x and y)
#' "c(x = 2 * (x * y * cos(x^2)), y = sin(x^2))"
#' 
#' # Compound function example (here abs(x) smoothed near 0)
#' fc <- function(x, h=0.1) if (abs(x) < h) 0.5*h*(x/h)**2 else abs(x)-0.5*h
#' Deriv("fc(x)", "x", cache.exp=FALSE)
#' "if (abs(x) < h) x/h else sign(x)"
#' 
#' # Example of a first argument that cannot be evaluated in the current environment:
#' \dontrun{
#'   suppressWarnings(rm("xx", "yy"))
#'   Deriv(xx^2+yy^2)
#' }
#' # c(xx = 2 * xx, yy = 2 * yy)
#' 
#' # Automatic differentiation (AD), note itermediate variable 'd' assignment
#' \dontrun{Deriv(~{d <- ((x-m)/s)^2; exp(-0.5*d)}, "x")}
#' #{
#' #   d <- ((x - m)/s)^2
#' #   .d_x <- 2 * ((x - m)/s^2)
#' #   -(0.5 * (.d_x * exp(-(0.5 * d))))
#' #}
#' 
#' # Custom derivation rule
#' \dontrun{
#'   myfun <- function(x, y=TRUE) NULL # do something usefull
#'   dmyfun <- function(x, y=TRUE) NULL # myfun derivative by x.
#'   drule[["myfun"]] <- alist(x=dmyfun(x, y), y=NULL) # y is just a logical
#'   Deriv(myfun(z^2, FALSE), "z")
#'   # 2 * (z * dmyfun(z^2, FALSE))
#' }
#' # Differentiantion by list components
#' \dontrun{
#'   theta <- list(m=0.1, sd=2.)
#'   x <- names(theta)
#'   names(x)=rep("theta", length(theta))
#'   Deriv(~exp(-(x-theta$m)**2/(2*theta$sd)), x, cache.exp=FALSE)
#' # c(theta_m = exp(-((x - theta$m)^2/(2 * theta$sd))) *
#' #  (x - theta$m)/theta$sd, theta_sd = 2 * (exp(-((x - theta$m)^2/
#' #  (2 * theta$sd))) * (x - theta$m)^2/(2 * theta$sd)^2))
#' }

Deriv <- function(f, x=if (is.function(f)) NULL else all.vars(if (is.character(f)) parse(text=f) else f), env=if (is.function(f)) environment(f) else parent.frame(), use.D=FALSE, cache.exp=TRUE, nderiv=NULL) {
	tf <- try(f, silent=TRUE)
	fch <- deparse(substitute(f))
	if (inherits(tf, "try-error")) {
		f <- substitute(f)
	}
	# create dsym and scache local envs (to keep clean nested calls)
	dsym <- new.env()
	dsym$l <- list()
	scache <- new.env()
	scache$l <- list()
	
	if (is.null(env))
		env <- .GlobalEnv
	if (is.null(x)) {
		# primitive function
		af <- formals(args(f))
		x <- names(af)
		rule <- drule[[fch]]
		if (!is.null(rule)) {
			# exclude arguments by which we cannot not differentiate from x
			x=as.list(x)
			x[sapply(rule, is.null)] <- NULL
			if (length(x) == 0) {
				stop(sprintf("There is no differentiable argument in the function %s", fch))
			}
			x=unlist(x)
		}
		fd <- as.call(c(as.symbol(fch), lapply(x, as.symbol)))
		pack_res <- as.call(alist(as.function, c(af, res), envir=env))
	} else {
		x[] <- as.character(x)
		if (any(nchar(x) == 0)) {
			stop("Names in the second argument must not be empty")
		}
		fd <- NULL
	}
	# prepare fd (a call to differentiate)
	# and pack_res (a cal to evaluate and return as result)
	if (!is.null(fd)) {
		; # we are already set
	} else if (is.character(f)) {
		# f is to parse
		fd <- parse(text=f)[[1]]
		pack_res <- as.call(alist(format1, res))
	} else if (is.function(f)) {
#browser()
		b <- body(f)
		if ((is.call(b) && (b[[1]] == as.symbol(".Internal") || b[[1]] == as.symbol(".External") || b[[1]] == as.symbol(".Call"))) || (is.null(b) && (is.primitive(f)) || !is.null(drule[[fch]]))) {
			if (fch %in% dlin || !is.null(drule[[fch]])) {
				arg <- lapply(names(formals(args(f))), as.symbol)
				fd <- as.call(c(as.symbol(fch), arg))
				pack_res <- as.call(alist(as.function, c(formals(args(f)), res), envir=env))
			} else {
				stop(sprintf("Internal or external function '%s()' is not in derivative table.", fch))
			}
		} else {
			fd <- b
			pack_res <- as.call(alist(as.function, c(formals(f), res), envir=env))
		}
	} else if (is.expression(f)) {
		fd <- f[[1]]
		pack_res <- as.call(alist(as.expression, res))
	} else if (is.language(f)) {
		if (is.call(f) && f[[1]] == as.symbol("~")) {
			# rhs of the formula
			fd <- f[[length(f)]]
			pack_res <- quote(res)
		} else {
			# plain call derivation
			fd <- f
			pack_res <- quote(res)
		}
	} else {
		fd <- substitute(f)
		pack_res <- quote(res)
		#stop("Invalid type of 'f' for differentiation")
	}
	res <- Deriv_(fd, x, env, use.D, dsym, scache)
	if (!is.null(nderiv)) {
		# multiple derivatives
		# prepare their names
		if (any(nderiv < 0)) {
			stop("All entries in nderiv must be non negative")
		}
		nm_deriv <- names(nderiv)
		nderiv <- as.integer(nderiv)
		if (is.null(nm_deriv))
			nm_deriv <- nderiv
		iempt <- nchar(nm_deriv)==0
		nm_deriv[iempt] <- seq_along(nderiv)[iempt]
		# prepare list of repeated derivatives
		lrep <- as.list(nderiv)
		names(lrep) <- nm_deriv
		# check if 0 is nderiv
		iz <- nderiv==0
		lrep[iz] <- list(fd)
		# set first derivative
		i <- nderiv==1
		lrep[i] <- list(res)
		
		maxd <- max(nderiv)
		for (ider in seq_len(maxd)) {
			if (ider < 2)
				next
			res <- Deriv_(res, x, env, use.D, dsym, scache)
			i <- ider == nderiv
			lrep[i] <- list(res)
		}
		if (length(lrep) == 1) {
			res <- lrep[[1]]
		} else {
			res <- as.call(c(quote(list), lrep))
		}
	}
#browser()
	if (cache.exp)
		res <- Cache(Simplify(deCache(res), scache=scache))
	eval(pack_res)
}

# workhorse function doing the main work of differentiation
Deriv_ <- function(st, x, env, use.D, dsym, scache) {
	stch <- as.character(if (is.call(st)) st[[1]] else st)
	# Make x scalar and wrap results in a c() call if length(x) > 1
	iel=which("..." == x)
	if (length(iel) > 0) {
		# remove '...' from derivable arguments
		x=as.list(x)
		x[iel]=NULL
		x=unlist(x)
	}
	nm_x <- names(x)
	if (!is.null(nm_x))
		nm_x[is.na(nm_x)] <- ""
	else
		nm_x <- rep("", length(x))
#browser()
	if (length(x) > 1 && stch != "{") {
		# many variables => recursive call on single name
		res <- lapply(seq_along(x), function(ix) Deriv_(st, x[ix], env, use.D, dsym, scache))
		names(res) <- if (is.null(nm_x)) x else ifelse(is.na(nm_x) | nchar(nm_x) == 0, x, paste(nm_x, x, sep="_"));
		return(as.call(c(as.symbol("c"), res)))
	}
	# differentiate R statement 'st' (a call, or a symbol or numeric) by a name in 'x'
	get_sub_x <- !(is.null(nm_x) | nchar(nm_x) == 0 | is.na(nm_x))
	is_index_expr <- is.call(st) && any(as.character(st[[1]]) == c("$", "[", "[["))
	is_sub_x <- is_index_expr &&
				format1(st[[2]]) == nm_x && format1(st[[3]]) == x
	if (is.conuloch(st) || (is_index_expr && !is_sub_x)) {
		return(0)
	} else if (is.symbol(st) || (get_sub_x && is_index_expr)) {
#browser()
		stch <- format1(st)
		if ((stch == x && !get_sub_x) || (get_sub_x && is_sub_x)) {
			return(1)
		} else if ((get_sub_x && is_index_expr && !is_sub_x) ||
				(if (get_sub_x) is.null(dsym$l[[nm_x]][[x]][[stch]]) else
				is.null(dsym$l[[x]][[stch]]))) {
			return(0)
		} else {
			return(if (get_sub_x) dsym$l[[nm_x]][[x]][[stch]] else dsym$l[[x]][[stch]])
		}
	} else if (is.call(st)) {
#browser()
		stch <- as.character(st[[1]])
		args <- as.list(st)[-1]
		if (stch %in% dlin) {
			# linear case
			# differentiate all arguments then pass them to the function
			dargs <- lapply(args, Deriv_, x, env, use.D, dsym, scache)
			return(Simplify_(as.call(c(st[[1]], dargs)), scache))
		}
		nb_args=length(st)-1
		# special cases: out of rule table or args(stch) -> NULL
		if (stch == "{") {
#browser()
			# AD differentiation (may be with many x)
			res=list(st[[1]])
			# initiate dsym[[x[ix]]] or dsym[[nm_x[ix]}}[[x[ix]]]
			for (ix in seq_along(x)) {
				if (get_sub_x[ix]) {
					if (is.null(dsym$l[[nm_x[ix]]][[x[ix]]]))
						dsym$l[[nm_x[ix]]][[x[ix]]] <- list()
				} else {
					if (is.null(dsym$l[[x[ix]]]))
						dsym$l[[x[ix]]] <- list()
				}
			}
			# collect defined var names (to avoid redifferentiation)
			defs <- sapply(as.list(st)[-1], function(e) if (is.assign(e)) as.character(e[[2]]) else "")
			alva=list()
			for (iarg in seq_along(args)) {
#browser()
				a <- args[[iarg]]
				if (is.assign(a)) {
					if (!is.symbol(a[[2]]))
						stop(sprintf("In AD mode, don't know how to deal with a non symbol '%s' at lhs", format1(a[[2]])))
					# put in scache the assignement
					Simplify_(a, scache)
					res <- append(res, a)
					alva <- append(alva, list(all.vars(a)))
					ach <- as.character(a[[2]])
					for (ix in seq_along(x)) {
						d_ach <- paste(".", ach, "_", x[ix], sep="")
						d_a <- as.symbol(d_ach)
						if (any(d_a == defs)) {
							# already differentiated in previous calls
							if (get_sub_x[ix])
								dsym$l[[nm_x[ix]]][[x[ix]]][[ach]] <- d_a
							else
								dsym$l[[x[ix]]][[ach]] <- d_a
							next
						}
						de_a <- Deriv_(a[[3]], x[ix], env, use.D, dsym, scache)
						if (de_a == 0) {
							if (iarg < length(args))
								next
						} else if (!is.call(de_a)) {
							if (get_sub_x[ix])
								dsym$l[[nm_x[ix]]][[x[ix]]][[ach]] <- de_a
							else
								dsym$l[[x[ix]]][[ach]] <- de_a
							if (iarg < length(args))
								next
						}
						if (get_sub_x[ix])
							dsym$l[[nm_x[ix]]][[x[ix]]][[ach]] <- d_a
						else
							dsym$l[[x[ix]]][[ach]] <- d_a
						res <- append(res, call("<-", d_a, de_a))
						alva <- append(alva, list(c(d_ach, all.vars(de_a))))
						# store it in scache too
						#scache$l[[format1(de_a)]] <- as.symbol(d_a)
					}
				} else {
					de_a <- lapply(seq_along(x), function(ix) Deriv_(a, x[ix], env, use.D, dsym, scache))
					if (length(x) > 1) {
						names(de_a) <- ifelse(get_sub_x, paste(nm_x, x, sep="_"), x)
						res <- append(res, as.call(c(as.symbol("c"), de_a)))
					} else {
						res <- append(res, de_a)
					}
				}
			}
#browser()
			if (length(alva) == length(res)) {
				i <- toporder(alva[-length(alva)]) # the last expression must stay the last
			} else {
				i <- toporder(alva)
			}
			res[-c(1, length(res))] <- res[-c(1, length(res))][i]
			return(Simplify(as.call(res)))
		} else if (is.uminus(st)) {
			return(Simplify(call("-", Deriv_(st[[2]], x, env, use.D, dsym, scache)), scache=scache))
		} else if (stch == "(") {
#browser()
			return(Simplify(Deriv_(st[[2]], x, env, use.D, dsym, scache), scache=scache))
		} else if(stch == "ifelse") {
			return(Simplify(call("ifelse", st[[2]], Deriv_(st[[3]], x, env, use.D, dsym, scache),
				Deriv_(st[[4]], x, env, use.D, dsym, scache)), scache=scache))
		} else if(stch == "rep") {
#browser()
			# 'x' argument is named or positional?
			i=if ("x" %in% names(st)) "x" else 2
			dst=st
			dst[[i]]=Simplify(Deriv_(st[[i]], x, env, use.D, dsym, scache), scache=scache)
			return(dst)
		} else if (stch == "if") {
			return(if (nb_args == 2)
				Simplify(call("if", st[[2]], Deriv_(st[[3]], x, env, use.D, dsym, scache)), scache=scache) else
				Simplify(call("if", st[[2]], Deriv_(st[[3]], x, env, use.D, dsym, scache),
					Deriv_(st[[4]], x, env, use.D, dsym, scache)), scache=scache))
		}
		rule <- drule[[stch]]
		if (is.null(rule)) {
#browser()
			# no derivative rule for this function
			# see if its arguments depend on x. If not, just send 0
			dargs <- lapply(args, Deriv_, x, env, use.D, dsym, scache)
			if (all(sapply(dargs, is.numeric)) && all(dargs == 0)) {
				return(0)
			}
			# otherwise try to get the body and differentiate it
			ff <- get(stch, mode="function", envir=env)
			bf <- body(ff)
			if (is.null(bf)) {
				stop(sprintf("Could not retrieve body of '%s()'", stch))
			}
			if (is.call(bf) && (bf[[1]] == as.symbol(".External") || bf[[1]] == as.symbol(".Internal") || bf[[1]] == as.symbol(".Call"))) {
#cat("aha\n")
				stop(sprintf("Function '%s()' is not in derivative table", stch))
			}
			mc <- match.call(ff, st)
			st <- Simplify_(do.call("substitute", list(bf, as.list(mc)[-1])), scache)
			return(Deriv_(st, x, env, use.D, dsym, scache))
		}
		# there is a rule!
		if (use.D) {
			return(Simplify(D(st, x), scache=scache))
		}
#if (stch == "myfun")
#browser()
		# prepare replacement list
		da <- try(args(stch), silent=TRUE)
		if (inherits(da, "try-error")) {
			# last chance to get unknown function definition
			# may be it is a user defined one?
			da <- args(get(stch, mode="function", envir=env))
		}
		mc <- as.list(match.call(definition=da, call=st))[-1]
		da <- as.list(da)
		da <- da[-length(da)] # all declared arguments with default values
		aa <- modifyList(da, mc) # all arguments with actual values
		# actualize the rule with actual arguments
		rule <- lapply(rule, function(r) do.call("substitute", list(r, aa)))
#browser()
		# which arguments can be differentiated?
		iad <- which(!sapply(rule, is.null))
		rule <- rule[iad]
		lsy <- unlist(lapply(dsym$l, function(it) if (get_sub_x && is.list(it)) unlist(lapply(it, ls, all.names=TRUE)) else ls(it, all.names=TRUE)))
		if (!any(names(which(sapply(mc, function(it) {av <- all.vars(it); (if (get_sub_x) any(nm_x == av) else any(x == av)) || any(av %in% lsy)}))) %in% names(rule))) {
			#warning(sprintf("A call %s cannot be differentiated by the argument '%s'", format1(st), x))
			return(0)
		}
		dargs <- lapply(names(rule), function(nm_a) if (is.null(mc[[nm_a]])) 0 else Deriv_(mc[[nm_a]], x, env, use.D, dsym, scache))
		ize <- sapply(dargs, `==`, 0)
		dargs <- dargs[!ize]
		rule <- rule[!ize]
		if (length(rule) == 0) {
			return(0)
		}
		
		# apply chain rule where needed
		ione <- sapply(dargs, `==`, 1)
		imone <- sapply(dargs, `==`, -1)
		for (i in seq_along(rule)[!(ione|imone)]) {
			rule[[i]] <- Simplify(call("*", dargs[[i]], rule[[i]]), scache=scache)
		}
		for (i in seq_along(rule)[imone]) {
			rule[[i]] <- Simplify(call("-", rule[[i]]), scache=scache)
		}
		return(Simplify(li2sum(rule), scache=scache))
	} else if (is.function(st)) {
#browser()
		# differentiate its body if can get it
		args <- as.list(st)[-1]
		names(args)=names(formals(ff))
		if (is.null(names(args))) {
			stop(sprintf("Could not retrieve arguments of '%s()'", stch))
		}
		st <- do.call("substitute", list(body(ff), args))
		Deriv_(st, x, env, use.D, dsym, scache)
	} else {
		stop("Invalid type of 'st' argument. It must be constant, symbol or a call.")
	}
}

drule <- new.env()

# linear functions, i.e. d(f(x))/dx == f(d(arg)/dx)
dlin=c("+", "-", "c", "t", "sum", "cbind", "rbind")

# rule table
# arithmetics
drule[["*"]] <- alist(e1=e2, e2=e1)
drule[["^"]] <- alist(e1=e2*e1^(e2-1), e2=e1^e2*log(e1))
drule[["/"]] <- alist(e1=1/e2, e2=-e1/e2^2)
# log, exp, sqrt
drule[["sqrt"]] <- alist(x=0.5/sqrt(x))
drule[["log"]] <- alist(x=1/(x*log(base)), base=-log(x, base)/(base*log(base)))
drule[["logb"]] <- drule[["log"]]
drule[["log2"]] <- alist(x=1/(x*log(2)))
drule[["log10"]] <- alist(x=1/(x*log(10)))
drule[["log1p"]] <- alist(x=1/(x+1))
drule[["exp"]] <- alist(x=exp(x))
drule[["expm1"]] <- alist(x=exp(x))
# trigonometric
drule[["sin"]] <- alist(x=cos(x))
drule[["cos"]] <- alist(x=-sin(x))
drule[["tan"]] <- alist(x=1/cos(x)^2)
drule[["asin"]] <- alist(x=1/sqrt(1-x^2))
drule[["acos"]] <- alist(x=-1/sqrt(1-x^2))
drule[["atan"]] <- alist(x=1/(1+x^2))
drule[["atan2"]] <- alist(y=x/(x^2+y^2), x=-y/(x^2+y^2))
if (getRversion() >= "3.1.0") {
	drule[["sinpi"]] <- alist(x=pi*cospi(x))
	drule[["cospi"]] <- alist(x=-pi*sinpi(x))
	drule[["tanpi"]] <- alist(x=pi/cospi(x)^2)
}
# hyperbolic
drule[["sinh"]] <- alist(x=cosh(x))
drule[["cosh"]] <- alist(x=sinh(x))
drule[["tanh"]] <- alist(x=(1-tanh(x)^2))
drule[["asinh"]] <- alist(x=1/sqrt(x^2+1))
drule[["acosh"]] <- alist(x=1/sqrt(x^2-1))
drule[["atanh"]] <- alist(x=1/(1-x^2))
# sign depending functions
drule[["abs"]] <- alist(x=sign(x))
drule[["sign"]] <- alist(x=0)
# special functions
drule[["besselI"]] <- alist(x=(if (nu == 0) besselI(x, 1, expon.scaled) else 0.5*(besselI(x, nu-1, expon.scaled) + besselI(x, nu+1, expon.scaled)))-if (expon.scaled) besselI(x, nu, TRUE) else 0, nu=NULL, expon.scaled=NULL)
drule[["besselK"]] <- alist(x=(if (nu == 0) -besselK(x, 1, expon.scaled) else -0.5*(besselK(x, nu-1, expon.scaled) + besselK(x, nu+1, expon.scaled)))+if (expon.scaled) besselK(x, nu, TRUE) else 0, nu=NULL, expon.scaled=NULL)
drule[["besselJ"]] <- alist(x=if (nu == 0) -besselJ(x, 1) else 0.5*(besselJ(x, nu-1) - besselJ(x, nu+1)), nu=NULL)
drule[["besselY"]] <- alist(x=if (nu == 0) -besselY(x, 1) else 0.5*(besselY(x, nu-1) - besselY(x, nu+1)), nu=NULL)
drule[["gamma"]] <- alist(x=gamma(x)*digamma(x))
drule[["lgamma"]] <- alist(x=digamma(x))
drule[["digamma"]] <- alist(x=trigamma(x))
drule[["trigamma"]] <- alist(x=psigamma(x, 2L))
drule[["psigamma"]] <- alist(x=psigamma(x, deriv+1L), deriv=NULL)
drule[["beta"]] <- alist(a=beta(a, b)*(digamma(a)-digamma(a+b)), b=beta(a, b)*(digamma(b)-digamma(a+b)))
drule[["lbeta"]] <- alist(a=digamma(a)-digamma(a+b), b=digamma(b)-digamma(a+b))
# probability densities
drule[["dbinom"]] <- alist(x=NULL, size=NULL, prob=if (size == 0) -x*(1-prob)^(x-1) else if (x == size) size*prob^(size-1) else (size-x*prob)*(x-size+1)*dbinom(x, size-1, prob)/(1-prob)^2/(if (log) dbinom(x, size, prob) else 1), log=NULL)
drule[["dnorm"]] <- alist(x=-(x-mean)/sd^2*if (log) 1 else dnorm(x, mean, sd),
	mean=(x-mean)/sd^2*if (log) 1 else dnorm(x, mean, sd),
	sd=(((x - mean)/sd)^2 - 1)/sd * if (log) 1 else dnorm(x, mean, sd),
	log=NULL)
drule[["pnorm"]] <- alist(q=dnorm(q, mean, sd)*(if (lower.tail) 1 else -1)/(if (log.p) pnorm(q, mean, sd, lower.tail) else 1), mean=dnorm(q, mean, sd)*(if (lower.tail) -1 else 1)/(if (log.p) pnorm(q, mean, sd, lower.tail) else 1), sd=dnorm(q, mean, sd)*(mean-q)/sd*(if (lower.tail) 1 else -1)/(if (log.p) pnorm(q, mean, sd, lower.tail) else 1), lower.tail=NULL, log.p=NULL)
