# Second order linear ODEs with 2F1 hypergeometric solutions, after
# sympy/solvers/ode/hypergeometric.py.

# hypergeometric.match_2nd_hypergeometric: list(A, B) of the normalized
# equation f'' + A f' + B f = 0, or NULL
.ode_match_2nd_hypergeometric <- function(eq, func) {
  x <- func$args[[1L]]
  df <- .ode_d(func, x)
  d2 <- .ode_d(func, x, 2L)
  ex <- list(func, df, d2)
  a3 <- Wild("a3", exclude = ex)
  b3 <- Wild("b3", exclude = ex)
  c3 <- Wild("c3", exclude = ex)
  deq <- a3 * d2 + b3 * df + c3 * func
  r <- sym_match(collect(eq, list(d2, df, func)), deq)
  if (length(r)) {
    if (!all(vapply(r, function(v) isTRUE(is_polynomial(v)), TRUE))) {
      nd <- .as_numer_denom(eq)
      eq <- expand(nd[[1L]])
      r <- sym_match(collect(eq, list(d2, df, func)), deq)
    }
  }
  if (length(r) && !identical((r$a3 %||% S_Zero)$key, S_Zero$key)) {
    A <- cancel(.div(r$b3 %||% S_Zero, r$a3))
    B <- cancel(.div(r$c3 %||% S_Zero, r$a3))
    return(list(A, B))
  }
  NULL
}

# the powers of x in an expression (the _power_counting helper)
.ode_power_counting <- function(vals, x) {
  pw <- list(S_Zero)
  addp <- function(p) if (!(p$key %in% vapply(pw, .k, ""))) pw[[length(pw) + 1L]] <<- p
  for (val in vals) {
    if (!.has1(val, x)) next
    if (val$type == "Pow" && identical(.base_exp(val)[[1L]]$key, x$key)) addp(.base_exp(val)[[2L]])
    else if (identical(val$key, x$key)) addp(S_One)
    else for (p in .ode_power_counting(.sympy_args(val), x)) addp(p)
  }
  pw
}

# gcd of rational numbers over QQ
.ode_gcd_rationals <- function(vals) {
  g <- S_Zero
  for (v in vals) {
    if (is_zero_num(g)) {
      g <- Abs(v)
      next
    }
    if (is_zero_num(v)) next
    a <- .as_numer_denom(g)
    b <- .as_numer_denom(v)
    num <- gcd(a[[1L]], b[[1L]])
    den <- lcm(a[[2L]], b[[2L]])
    g <- Abs(.div(num, den))
  }
  g
}

# hypergeometric.equivalence
.ode_hyper_equivalence <- function(max_num_pow, dem_pow) {
  dp <- vapply(dem_pow, .ode_num, 0)
  same <- function(opts) any(vapply(opts, function(o) length(o) == length(dp) && all(o == dp), TRUE))
  if (max_num_pow == 2) {
    if (same(list(c(2, 2), c(2, 2, 2)))) return("2F1")
  } else if (max_num_pow == 1) {
    if (same(list(c(1, 2, 2), c(2, 2, 2), c(1, 2), c(2, 2)))) return("2F1")
  } else if (max_num_pow == 0) {
    if (same(list(c(1, 1, 2), c(2, 2), c(1, 2, 2), c(1, 1), 2, c(1, 2), c(2, 2)))) return("2F1")
  }
  NULL
}

# hypergeometric.equivalence_hypergeometric
.ode_equivalence_hypergeometric <- function(A, B, func) {
  x <- func$args[[1L]]
  I1 <- factor(cancel(.add(list(.div(.diff1(A, x), S_Two), .div(.pow(A, S_Two), .int(4)), .neg(B)))))
  J1 <- factor(cancel(.add(list(.mul(list(.pow(x, S_Two), I1)), .rat(1L, 4L)))))
  nd <- .as_numer_denom(J1)
  num <- powdenest(expand(nd[[1L]]))
  dem <- powdenest(expand(nd[[2L]]))
  pw <- .ode_power_counting(list(num), x)
  pd <- .ode_power_counting(list(dem), x)
  for (p in pw) if (!(p$key %in% vapply(pd, .k, ""))) pd[[length(pd) + 1L]] <- p
  k <- .ode_gcd_rationals(pd)
  I0 <- powdenest(simplify(factor(.div(.sub(.div(J1, .pow(k, S_Two)), .rat(1L, 4L)), .pow(.pow(x, k), S_Two)))), force = TRUE)
  I0 <- factor(cancel(powdenest(subs(I0, x, .pow(x, .inv(k))), force = TRUE)))
  if (!.ode_is_rational_function(I0, x)) return(NULL)
  nd <- .as_numer_denom(I0)
  num <- nd[[1L]]
  dem <- nd[[2L]]
  max_num_pow <- max(vapply(.ode_power_counting(list(num), x), .ode_num, 0))
  sing_point <- list()
  dem_pow <- list()
  for (arg in .sympy_args(dem)) {
    if (!.has1(arg, x)) next
    be <- .base_exp(arg)
    dem_pow[[length(dem_pow) + 1L]] <- be[[2L]]
    base <- if (arg$type == "Pow") be[[1L]] else arg
    rs <- roots(base, x)
    sing_point[[length(sing_point) + 1L]] <- rs[[1L]]
  }
  if (length(dem_pow)) dem_pow <- dem_pow[order(vapply(dem_pow, .ode_num, 0))]
  if (identical(.ode_hyper_equivalence(max_num_pow, dem_pow), "2F1"))
    return(list(I0 = I0, k = k, sing_point = sing_point, type = "2F1"))
  NULL
}

# hypergeometric.match_2nd_2F1_hypergeometric
.ode_match_2nd_2F1_hypergeometric <- function(I, k, sing_point, func) {
  x <- func$args[[1L]]
  c <- Wild("c")
  s <- Wild("s")
  r <- Wild("r")
  alpha <- Wild("alpha")
  beta <- Wild("beta")
  gamma <- Wild("gamma")
  delta <- Wild("delta")
  a <- Wild("a")
  b <- Wild("b")
  I0 <- .div(.add(list(.mul(list(.add(list(a, .neg(b), S_One)), .add(list(a, .neg(b), S_NegativeOne)), .pow(x, S_Two))),
                       .mul(list(S_Two, .add(list(.mul(list(.sub(.sub(S_One, a), b), c)), .mul(list(S_Two, a, b)))), x)),
                       .mul(list(c, .sub(c, S_Two))))),
             .mul(list(.int(4), .pow(x, S_Two), .pow(.sub(x, S_One), S_Two))))
  is01 <- length(sing_point) == 2L && is_zero_num(sing_point[[1L]]) && identical(sing_point[[2L]]$key, S_One$key)
  if (!is01) {
    if (length(sing_point) < 2L) .ode_stop("IndexError", "list index out of range")
    .beta <- .neg(.mul(list(alpha, sing_point[[1L]])))
    .delta <- .neg(.mul(list(gamma, sing_point[[2L]])))
    .gamma <- alpha
    if (length(sing_point) == 3L)
      .gamma <- .div(.add(list(.beta, .mul(list(sing_point[[3L]], alpha)))), .sub(sing_point[[3L]], sing_point[[2L]]))
    mob <- .div(.add(list(.mul(list(alpha, x)), beta)), .add(list(.mul(list(gamma, x)), delta)))
    mob <- subs(mob, beta, .beta)
    mob <- subs(mob, delta, .delta)
    mob <- subs(mob, gamma, .gamma)
    mob <- cancel(mob)
    t <- .div(.sub(beta, .mul(list(delta, x))), .sub(.mul(list(gamma, x)), alpha))
    t <- cancel(subs(subs(subs(t, beta, .beta), delta, .delta), gamma, .gamma))
  } else {
    mob <- x
    t <- x
  }
  I <- subs(I, x, t)
  I <- .mul(list(I, .pow(.diff1(t, x), S_Two)))
  I <- factor(I)
  I0_dem <- .as_numer_denom(I0)[[2L]]
  x2 <- .pow(x, S_Two)
  cd <- collect(expand(cancel(.mul(list(I, I0_dem)))), list(x2, x), evaluate = FALSE)
  getc <- function(key) {
    v <- cd[[.str(key)]]
    if (is.null(v)) S_Zero else v
  }
  lhs <- list(getc(x2), getc(x), getc(S_One))
  rhs <- list(.sub(.pow(s, S_Two), S_One),
              .add(list(.mul(list(S_Two, .sub(S_One, r), c)), .mul(list(.add(list(r, s)), .sub(r, s))))),
              .mul(list(c, .sub(c, S_Two))))
  .c <- .sub(S_One, factor(.sqrt(.add(list(S_One, lhs[[3L]])))))
  if (!length(atoms(.c, c("Symbol", "Dummy", "Wild")))) {
    cr <- roots(.sub(lhs[[3L]], rhs[[3L]]), c)
    cv <- vapply(cr, function(z) as.complex(.neval(z, list())), 0i)
    if (any(abs(Im(cv)) > 0)) .ode_stop("TypeError", "Invalid comparison of non-real")
    .c <- cr[[which.min(Re(cv))]]
  }
  .s <- factor(.sqrt(.add(list(lhs[[1L]], S_One))))
  .r <- .sub(.c, factor(.sqrt(.add(list(.pow(.c, S_Two), .pow(.s, S_Two), lhs[[2L]], .neg(.mul(list(S_Two, .c))))))))
  .a <- .div(.add(list(.r, .s)), S_Two)
  .b <- .div(.sub(.r, .s), S_Two)
  list(a = simplify(.a), b = simplify(.b), c = simplify(.c), k = k, mobius = mob, type = "2F1")
}

# hypergeometric.get_sol_2F1_hypergeometric
.ode_get_sol_2F1_hypergeometric <- function(eq, func, match_object) {
  x <- func$args[[1L]]
  cs <- .ode_get_numbered_constants(eq, num = 2L)
  C0 <- cs[[1L]]
  C1 <- cs[[2L]]
  a <- match_object$a
  b <- match_object$b
  c <- match_object$c
  A <- match_object$A
  sol <- NULL
  if (isFALSE(ask_fact(c, "integer"))) {
    sol <- .add(list(.mul(list(C0, hyper(list(a, b), list(c), x))),
                     .mul(list(C1, hyper(list(.add(list(a, .neg(c), S_One)), .add(list(b, .neg(c), S_One))), list(.sub(S_Two, c)), x),
                               .pow(x, .sub(S_One, c))))))
  } else if (identical(c$key, S_One$key)) {
    h <- hyper(list(a, b), list(c), x)
    y2 <- .mul(list(Integral(.div(fcall("exp", list(Integral(.div(.add(list(.mul(list(.neg(.add(list(a, b, S_One))), x)), c)),
                                                                    .sub(.pow(x, S_Two), x)), x))),
                                  .pow(hyperexpand(h), S_Two)), x), h))
    sol <- .add(list(.mul(list(C0, h)), .mul(list(C1, y2))))
  } else if (isFALSE(ask_fact(.sub(.sub(c, a), b), "integer"))) {
    omx <- .sub(S_One, x)
    sol <- .add(list(.mul(list(C0, hyper(list(a, b), list(.add(list(S_One, a, b, .neg(c)))), omx))),
                     .mul(list(C1, hyper(list(.sub(c, a), .sub(c, b)), list(.add(list(S_One, c, .neg(a), .neg(b)))), omx),
                               .pow(omx, .sub(.sub(c, a), b))))))
  }
  if (is.null(sol) || !.ode_true(sol)) return(sol)
  mob <- match_object$mobius
  dtdx <- simplify(.inv(.diff1(mob, x)))
  .B <- .mul(list(subs(.sub(.mul(list(.add(list(a, b, S_One)), x)), c), x, mob), dtdx))
  .B <- factor(.add(list(.B, .mul(list(subs(.sub(.pow(x, S_Two), x), x, mob), .diff1(dtdx, x), dtdx)))))
  .A <- factor(.mul(list(subs(.sub(.pow(x, S_Two), x), x, mob), .pow(dtdx, S_Two))))
  e <- fcall("exp", list(logcombine(Integral(cancel(.div(.B, .mul(list(S_Two, .A)))), x), force = TRUE)))
  k <- match_object$k
  sol <- subs(sol, x, mob)
  sol <- subs(sol, x, .pow(x, k))
  e <- subs(e, x, .pow(x, k))
  xk <- .pow(x, .div(.add(list(.neg(k), S_One)), S_Two))
  if (!isTRUE(.sp_fact(A, "zero"))) {
    e1 <- fcall("exp", list(logcombine(Integral(.div(A, S_Two), x), force = TRUE)))
    return(Eq(func, .mul(list(cancel(.mul(list(.div(e, e1), xk))), sol))))
  }
  Eq(func, .mul(list(cancel(.mul(list(e, xk))), sol)))
}
