library(testthat)


devtools::install()


expect_identical_expression = function(a, b) {
  # testthat compares result of deparse which doesn't always match even though the expressions are identical
  a = .reparse_expression(a)
  b = .reparse_expression(b)
  testthat::expect_identical(a, b)
}


.reparse_expression = function(expression) {
  result = parse(text=deparse(expression), keep.source=F)
  return(result)
}


## matrix

matrix_derivatives = list(data=quote(matrix(.derivative, nrow=nrow, ncol=ncol, byrow=byrow, dimnames=dimnames)))
# assumes
# - nrow and ncol are always specified

assign('matrix', matrix_derivatives, envir=Deriv::drule)

test_that('symbolic derivatives of matrix', {
  expect_identical_expression(Deriv::Deriv(quote(
    matrix(2)
  ), 'x'),
    0
  )
  
  expect_identical_expression(Deriv::Deriv(quote(
    matrix(2*x)
  ), 'x'), quote(
    matrix(2, nrow=1, ncol=1, byrow=FALSE, dimnames=NULL)
  ))
  
  expect_identical_expression(Deriv::Deriv(quote(
    matrix(c(2*x, 3), nrow=2)
  ), 'x'), quote(
    matrix(c(2, 0), nrow=2, ncol=1, byrow=FALSE, dimnames=NULL)
  ))
})


## matrix multiplication

assign('%*%', list(x=quote(y), y=quote(x)), envir=Deriv::drule)

test_that('symbolic derivatives of %*%', {
  expect_identical_expression(Deriv::Deriv(quote(
    matrix(2, nrow=1, ncol=1)%*%matrix(3, nrow=1, ncol=1)
  ), 'x'),
    0
  )
  
  expect_identical_expression(Deriv::Deriv(quote(
    matrix(2*x, nrow=1, ncol=1)%*%matrix(3*y, nrow=1, ncol=1)
  ), 'x'), quote(
    matrix(2, nrow=1, ncol=1, byrow=FALSE, dimnames=NULL)%*%matrix(3*y, nrow=1, ncol=1)
  ))
  
  expect_identical_expression(Deriv::Deriv(quote(
    matrix(2*y, nrow=1, ncol=1)%*%matrix(3*x, nrow=1, ncol=1)
  ), 'x'), quote(
    matrix(2*y, nrow=1, ncol=1)%*%matrix(3, nrow=1, ncol=1, byrow=FALSE, dimnames=NULL)
  ))
  
  expect_identical_expression(Deriv::Deriv(quote(
    matrix(c(2*x, 3), nrow=1, ncol=2)%*%matrix(c(4*x, 5), nrow=2, ncol=1)
  ), 'x'), quote(
    matrix(c(2*x, 3), nrow=1, ncol=2)%*%matrix(c(4, 0), nrow=2, ncol=1, byrow=FALSE, dimnames=NULL)
    + matrix(c(2, 0), nrow=1, ncol=2, byrow=FALSE, dimnames=NULL)%*%matrix(c(4*x, 5), nrow=2, ncol=1)
  ))
})


## det

det_derivatives = list(x=quote(det(x)*sum(diagonal(solve(x)%*%(.derivative)))))  # TODO only works for invertible matrices
# assumes
# - matrix is invertible
# - function diagonal exists (see below)

assign('det', det_derivatives, envir=Deriv::drule)

test_that('symbolic derivatives of det', {
  expect_identical_expression(Deriv::Deriv(quote(
    det(matrix(2, nrow=1, ncol=1))
  ), 'x'),
    0
  )
  
  expect_identical_expression(Deriv::deCache(Deriv::Deriv(quote(
    det(matrix(2*x, nrow=1, ncol=1))
  ), 'x')), quote(
    det(matrix(2*x, nrow=1, ncol=1))*sum(diagonal(
      solve(matrix(2*x, nrow=1, ncol=1))%*%matrix(2, nrow=1, ncol=1, byrow=FALSE, dimnames=NULL)
    ))
  ))
})


## solve

solve_derivatives = list(a=quote(-solve(a)%*%(.derivative)%*%solve(a)),
                         b=quote(not_implemented))
# assumes
# - solve is called with a single argument

assign('solve', solve_derivatives, envir=Deriv::drule)

test_that('symbolic derivatives of solve', {
  expect_identical_expression(Deriv::Deriv(quote(
    solve(matrix(2, nrow=1, ncol=1))
  ), 'x'),
    0
  )
  
  expect_identical_expression(Deriv::deCache(Deriv::Deriv(quote(
    solve(matrix(2*x, nrow=1, ncol=1))
  ), 'x')), quote(
    -solve(matrix(2*x, nrow=1, ncol=1))%*%matrix(2, nrow=1, ncol=1, byrow=FALSE, dimnames=NULL)
    %*%solve(matrix(2*x, nrow=1, ncol=1))
  ))
  
  expect_identical_expression(Deriv::deCache(Deriv::Deriv(quote(
    solve(matrix(c(x^3, x^2, x, 1), nrow=2, ncol=2))
  ), 'x')), quote(
    -solve(matrix(c(x^3, x^2, x, 1), nrow=2, ncol=2))
    %*%matrix(c(3*x^2, 2*x, 1, 0), nrow=2, ncol=2, byrow=FALSE, dimnames=NULL)
    %*%solve(matrix(c(x^3, x^2, x, 1), nrow=2, ncol=2))
  ))
})


## diag

diag_derivatives = list(x=quote(if (!is.matrix(x) && length(x) == 1) matrix(0, nrow=x, ncol=x) else diag(x=.derivative, nrow, ncol, names=names)))

assign('diag', diag_derivatives, envir=Deriv::drule)

test_that('symbolic derivatives of diag', {
  expect_identical_expression(Deriv::Deriv(quote(
    diag(2)
  ), 'x'),
    0
  )
  
  expect_identical_expression(Deriv::deCache(Deriv::Deriv(quote(
    diag(2*x)
  ), 'x')), quote(
    if (!is.matrix(2*x) && length(2*x) == 1) matrix(0, nrow=2*x, ncol=2*x)
    else diag(x=2, , , names=TRUE)
  ))
})


## diag alternatives

# diagonal

diagonal = function(x) diag(x)

assign('diagonal', list(x=quote(diagonal(.derivative))), envir=Deriv::drule)

test_that('symbolic derivatives of diagonal', {
  expect_identical_expression(Deriv::Deriv(quote(
    diagonal(matrix(2))
  ), 'x'),
  0
  )
  
  expect_identical_expression(Deriv::Deriv(quote(
    diagonal(matrix(2*x))
  ), 'x'), quote(
    diagonal(matrix(2, nrow=1, ncol=1, byrow=FALSE, dimnames=NULL))
  ))
})

# identity matrix

identity_matrix = function(x) diag(x)

assign('identity_matrix', list(), envir=Deriv::drule)  # empty list means derivative is 0 and body of function is not considered

test_that('symbolic derivatives of identity_matrix', {
  expect_identical_expression(Deriv::Deriv(quote(
    identity_matrix(2)
  ), 'x'),
    0
  )
  
  expect_identical_expression(Deriv::Deriv(quote(
    identity_matrix(2*x)
  ), 'x'),
    0
  )
})


## other

# length

assign('length', list(x=quote(0)), envir=Deriv::drule)

test_that('symbolic derivatives of length', {
  expect_identical_expression(Deriv::Deriv(quote(
    length(2)
  ), 'x'),
    0
  )
  
  expect_identical_expression(Deriv::Deriv(quote(
    length(2*x)
  ), 'x'),
    0
  )
  
  expect_identical_expression(Deriv::Deriv(quote(
    2*x*length(3*x)
  ), 'x'), quote(
    2*length(3*x)
  ))
})


# $

test_that('symbolic derivatives of $', {
  expect_identical_expression(Deriv::Deriv(quote(
    list(a=2)$a
  ), 'x'), quote(
    list(a=0)$a
  ))
  
  expect_identical_expression(Deriv::Deriv(quote(
    list(a=2*x)$a
  ), 'x'), quote(
    list(a=2)$a
  ))
  
  expect_identical_expression(Deriv::Deriv(quote(
    list(x=2)$x
  ), 'x'), quote(
    list(x=0)$x
  ))
  
  expect_identical_expression(Deriv::Deriv(quote(
    list(x=2*x)$x
  ), 'x'), quote(
    list(x=2)$x
  ))
})


# [

test_that('symbolic derivatives of [', {
  expect_identical_expression(Deriv::Deriv(quote(
    2[1]
  ), 'x'),
    0
  )
  
  expect_identical_expression(Deriv::Deriv(quote(
    (2*x)[1]
  ), 'x'), quote(
    2
  ))
  
  expect_identical_expression(Deriv::Deriv(quote(
    (2*x)[3*x]
  ), 'x'), quote(
    2[3*x]
  ))
  
  expect_identical_expression(Deriv::Deriv(quote(
    2[3*x]
  ), 'x'), quote(
    0[3*x]
  ))
})


# [[

test_that('symbolic derivatives of [[', {
  expect_identical_expression(Deriv::Deriv(quote(
    2[[1]]
  ), 'x'),
    0
  )
  
  expect_identical_expression(Deriv::Deriv(quote(
    (2*x)[[1]]
  ), 'x'), quote(
    2
  ))
  
  expect_identical_expression(Deriv::Deriv(quote(
    (2*x)[[3*x]]
  ), 'x'), quote(
    2[[3*x]]
  ))
  
  expect_identical_expression(Deriv::Deriv(quote(
    2[[3*x]]
  ), 'x'), quote(
    0[[3*x]]
  ))
})


# qnorm

qnorm_derivatives = list(p=quote(1/dnorm(qnorm(p, mean=mean, sd=sd, lower.tail=lower.tail, log.p=log.p), mean=mean, sd=sd)),
                         mean=quote(1),
                         sd=quote((qnorm(p, mean=mean, sd=sd, lower.tail=lower.tail, log.p=log.p) - mean)/sd))

assign('qnorm', qnorm_derivatives, envir=Deriv::drule)

test_that('symbolic derivatives of qnorm', {
  expect_identical_expression(Deriv::Deriv(quote(
    qnorm(2)
  ), 'x'),
    0
  )
  
  expect_identical_expression(Deriv::Deriv(quote(
    qnorm(2*x)
  ), 'x'), quote(
    2/dnorm(qnorm(2*x, mean=0, sd=1, lower.tail=TRUE, log.p=FALSE), mean=0, sd=1)
  ))
  
  expect_identical_expression(Deriv::Deriv(quote(
    qnorm(2*x, mean=3*mu)
  ), 'mu'),
    3
  )
  
  expect_identical_expression(Deriv::deCache(Deriv::Deriv(quote(
    qnorm(2*x, mean=3*mu, sd=4*sd)
  ), 'sd')), quote(
    (qnorm(2*x, mean=3*mu, sd=4*sd, lower.tail=TRUE, log.p=FALSE) - 3*mu)/sd
  ))
})


## real example

x = c(1, 2); x
mu = c(1 + 1/3, 2 - 1/4); mu
Sigma = matrix(c(1, 1/2, 1/2, 3), nrow=2, ncol=2); Sigma

u = pnorm(x, mean=mu, sd=sqrt(diag(Sigma))); u
rho = Sigma[1, 2] / prod(sqrt(diag(Sigma))); rho
Rho = matrix(c(1, rho, rho, 1), nrow=2, ncol=2); Rho


# normal copula density

normal_copula_density = quote(det(Rho)^(-1/2)*exp(-1/2*(t(qnorm(u))%*%(solve(Rho) - identity_matrix(length(u)))%*%qnorm(u))[1]))

eval(normal_copula_density)*prod(dnorm(x, mean=mu, sd=sqrt(diag(Sigma))))


# derivatives

. = Deriv::deCache( Deriv::Deriv(substituteDirect(normal_copula_density, list(u=quote(c(u1, u2)))), 'u1') ); .
eval(., list(u1=u[1], u2=u[2], Rho=Rho))  # 0.2052829

. = Deriv::deCache( Deriv::Deriv(., 'u1') ); .
eval(., list(u1=u[1], u2=u[2], Rho=Rho))  # -0.7928689

. = Deriv::deCache( Deriv::Deriv(substituteDirect(normal_copula_density, list(Rho=quote(matrix(c(1, rho, rho, 1), nrow=2, ncol=2)))), 'rho') ); .
eval(.)  # 0.2122524
