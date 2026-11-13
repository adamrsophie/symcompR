# Rational particular and general solutions of Riccati ODEs, after
# sympy/solvers/ode/riccati.py (Kovacic's algorithm for the normal form
# f'(x) + f(x)^2 = a(x)). The polynomials of the original are kept as
# expressions in x here, their coefficients brought to the canonical form of
# an algebraic number field element where SymPy converts them.

# an expression from a Poly or an expression
.ode_rx <- function(p) if (inherits(p, "Poly")) as_expr(p) else as_sym(p)

# Poly(e, x).degree(), -Inf for the zero polynomial
.ode_rdeg <- function(e, x) {
  e <- expand(.ode_rx(e))
  if (is_zero_num(e)) return(-Inf)
  as.numeric(degree(e, x))
}

# the coefficients of a polynomial in x, highest degree first
.ode_rcoeffs <- function(e, x) {
  e <- expand(.ode_rx(e))
  if (is_zero_num(e)) return(list(S_Zero))
  lapply(all_coeffs(Poly(e, x)), as_sym)
}
.ode_rlc <- function(e, x) .ode_rcoeffs(e, x)[[1L]]

# an element of QQ<a> as AlgebraicField.to_sympy writes it
.ode_anp <- function(e) expand(radsimp(e))

# Poly.cancel(g, include=True) over the coefficient field: the reduced
# numerator and denominator, the denominator with a positive leading
# coefficient
.ode_rcancel <- function(n, d, x) {
  f <- expand(.ode_rx(n))
  g <- expand(.ode_rx(d))
  r <- cancel(.div(f, g), extension = TRUE)
  nd <- .as_numer_denom(r)
  p <- expand(nd[[1L]])
  q <- expand(nd[[2L]])
  if (.ode_alg_domain(c(.ode_rcoeffs(f, x), .ode_rcoeffs(g, x)))) {
    # over QQ<a> the gcd is monic, so q keeps the leading coefficient of g
    s <- .div(.ode_rlc(g, x), .ode_rlc(q, x))
    p <- .ode_anp_poly(.mul(list(p, s)), x)
    q <- .ode_anp_poly(.mul(list(q, s)), x)
  }
  if (.ode_anp_negative(.ode_rlc(q, x))) {
    p <- expand(.neg(p))
    q <- expand(.neg(q))
  }
  list(p, q)
}

# whether construct_domain(coeffs, extension=True) gives an algebraic field
.ode_alg_domain <- function(cs) {
  any(vapply(cs, function(c) !(c$type %in% c("Integer", "Rational")) && is.null(.cd_pure_complex(c)) &&
               !length(free_symbols(c)) && isTRUE(ask_fact(c, "algebraic")), TRUE))
}

# AlgebraicField.is_negative: the sign of the leading coefficient of the
# element as a polynomial in the generator
.ode_anp_negative <- function(c) {
  e <- expand(c)
  if (e$type %in% c("Integer", "Rational")) return(.ode_num(e) < 0)
  rad <- Filter(function(p) p$args[[2L]]$type == "Rational", atoms(e, "Pow"))
  if (length(rad)) {
    b <- coeff(e, rad[[1L]])
    if (b$type %in% c("Integer", "Rational") && !is_zero_num(b)) return(.ode_num(b) < 0)
  }
  isTRUE(ask_fact(e, "negative"))
}

# riccati.riccati_normal
.ode_riccati_normal <- function(w, x, b1, b2) {
  .sub(.sub(.mul(list(.neg(b2), w)), .div(.diff1(b2, x), .mul(list(S_Two, b2)))), .div(b1, S_Two))
}

# riccati.riccati_inverse_normal
.ode_riccati_inverse_normal <- function(y, x, b1, b2, bp = NULL) {
  if (is.null(bp))
    bp <- .sub(.div(.neg(.diff1(b2, x)), .mul(list(S_Two, .pow(b2, S_Two)))), .div(b1, .mul(list(S_Two, b2))))
  .add(list(.div(.neg(y), b2), bp))
}

.ode_riccati_a <- function(x, b0, b1, b2) {
  d1 <- .diff1(b1, x)
  d2 <- .diff1(b2, x)
  .add(list(.mul(list(.neg(b0), b2)), .div(.pow(b1, S_Two), .int(4)), .neg(.div(d1, S_Two)),
            .div(.mul(list(.int(3), .pow(d2, S_Two))), .mul(list(.int(4), .pow(b2, S_Two)))),
            .div(.mul(list(b1, d2)), .mul(list(S_Two, b2))),
            .neg(.div(.diff_spec_eval(b2, list(list(x, 2L))), .mul(list(S_Two, b2))))))
}

# riccati.riccati_reduced
.ode_riccati_reduced <- function(eq, f, x) {
  mf <- .ode_match_riccati(eq, f, x)
  if (!isTRUE(mf[[1L]])) return(FALSE)
  b <- mf[[2L]]
  a <- .ode_riccati_a(x, b[[1L]], b[[2L]], b[[3L]])
  fx <- if (is.function(f)) f(x) else f
  .add(list(.ode_d(fx, x), .pow(fx, S_Two), .neg(a)))
}

# Expr.is_rational_function(x)
.ode_is_rational_function <- function(e, x) {
  if (e$type %in% c("NaN", "Infinity", "NegativeInfinity", "ComplexInfinity")) return(FALSE)
  if (!.has1(e, x)) return(TRUE)
  .is_rational_function_in(e, list(x))
}

# riccati.match_riccati
.ode_match_riccati <- function(eq, f, x) {
  fx <- if (is.function(f)) f(x) else f
  no <- .pytuple(list(FALSE, list()))
  if (eq$type == "Equality") eq <- .sub(eq$args[[1L]], eq$args[[2L]])
  eq <- collect(expand(eq), list(fx))
  df <- .ode_d(fx, x)
  cf <- coeff(eq, df)
  if (!is_zero_num(cf) && eq$type == "Add") {
    eq <- collect(.add(lapply(eq$args, function(t) cancel(.div(t, cf)))), list(fx))
    b1 <- .neg(coeff(eq, fx))
    b2 <- .neg(coeff(eq, .pow(fx, S_Two)))
    b0 <- expand(.add(list(df, .neg(.mul(list(b1, fx))), .neg(.mul(list(b2, .pow(fx, S_Two)))), .neg(eq))))
    funcs <- list(b0, b1, b2)
    if (any(vapply(funcs, function(z) length(atoms(z, c("Symbol", "Dummy"))) > 1L || length(atoms(z, "Float")) > 0L, TRUE)))
      return(no)
    hasf <- length(Filter(function(a) identical(a$name, fx$name), atoms(b0, "Function"))) > 0L
    if (hasf || identical(b2$key, S_Zero$key) || !all(vapply(funcs, .ode_is_rational_function, TRUE, x = x)))
      return(no)
    return(.pytuple(list(TRUE, funcs)))
  }
  no
}

# riccati.val_at_inf
.ode_val_at_inf <- function(num, den, x) {
  v <- .ode_rdeg(den, x) - .ode_rdeg(num, x)
  if (is.finite(v)) as.integer(v) else v
}

.ode_num <- function(v) if (is_sym(v)) as.numeric(.neval(v, list())) else as.numeric(v)

# riccati.check_necessary_conds
.ode_check_necessary_conds <- function(val_inf, muls) {
  v <- .ode_num(val_inf)
  ms <- vapply(muls, .ode_num, 0)
  (v >= 2 || (v <= 0 && v %% 2 == 0)) && all(ms == 1 | (ms %% 2 == 0 & ms >= 2))
}

# x**n * f(1/x) for a polynomial f of degree n (Poly.transform(1, x))
.ode_reverse_poly <- function(e, x) {
  cs <- .ode_rcoeffs(e, x)
  n <- length(cs) - 1L
  .add(lapply(seq_along(cs), function(i) .mul(list(cs[[i]], .pow(x, .int(i - 1L))))))
}

.ode_inverse_transform <- function(num, den, x) {
  num <- .ode_rx(num)
  den <- .ode_rx(den)
  pwr <- .ode_val_at_inf(num, den, x)
  if (pwr >= 0) {
    if (!is_zero_num(num)) {
      num <- expand(.mul(list(.ode_reverse_poly(num, x), .pow(x, .int(pwr)))))
      den <- .ode_reverse_poly(den, x)
    }
  } else {
    num <- .ode_reverse_poly(num, x)
    den <- expand(.mul(list(.ode_reverse_poly(den, x), .pow(x, .int(-pwr)))))
  }
  .ode_rcancel(num, den, x)
}

# riccati.inverse_transform_poly
.ode_inverse_transform_poly <- function(num, den, x) {
  nd <- .ode_inverse_transform(num, den, x)
  .pytuple(list(Poly(nd[[1L]], x), Poly(nd[[2L]], x)))
}

# riccati.limit_at_inf
.ode_limit_at_inf <- function(num, den, x) {
  pwr <- -.ode_val_at_inf(num, den, x)
  if (pwr > 0) return(.mul(list(S_Infinity, fcall("sign", list(.div(.ode_rlc(num, x), .ode_rlc(den, x)))))))
  if (pwr == 0) return(.div(.ode_rlc(num, x), .ode_rlc(den, x)))
  S_Zero
}

# riccati.construct_c_case_1
.ode_construct_c_case_1 <- function(num, den, x, pole) {
  nd <- .ode_rcancel(.mul(list(.ode_rx(num), expand(.pow(.sub(x, pole), S_Two)))), den, x)
  r <- .div(.ode_anp(subs(nd[[1L]], x, pole)), .ode_anp(subs(nd[[2L]], x, pole)))
  if (!identical(r$key, .rat(-1L, 4L)$key)) {
    s <- .sqrt(.add(list(S_One, .mul(list(.int(4), r)))))
    return(list(list(.div(.add(list(S_One, s)), S_Two)), list(.div(.sub(S_One, s), S_Two))))
  }
  list(list(S_Half))
}

# a Laurent series dictionary: a pydict, or a list of (key, value) pairs
.ode_ser_get <- function(ser, k) {
  if (inherits(ser, "pydict")) {
    ks <- vapply(ser$keys, function(z) .ode_num(z), 0)
    return(ser$values[[match(k, ks)]])
  }
  for (p in ser) if (.ode_num(p[[1L]]) == k) return(p[[2L]])
  .ode_stop("KeyError", as.character(k))
}

.ode_same_list <- function(a, b) length(a) == length(b) && all(vapply(seq_along(a), function(i) identical(a[[i]]$key, b[[i]]$key), TRUE))

# riccati.construct_c_case_2
.ode_construct_c_case_2 <- function(num, den, x, pole, mul) {
  mul <- .ode_num(mul)
  ri <- mul %/% 2
  ser <- .ode_rational_laurent_series(num, den, x, pole, mul, 6L)
  sg <- function(k) .ode_ser_get(ser, k)
  cplus <- rep(list(S_Zero), ri)
  cplus[[ri]] <- .sqrt(sg(2 * ri))
  s <- ri - 1
  sm <- S_Zero
  if (ri - 1 >= 1) {
    for (s in seq(ri - 1, 1)) {
      sm <- S_Zero
      if (s + 1 <= ri - 1) for (j in seq(s + 1, ri - 1)) sm <- .add(list(sm, .mul(list(cplus[[j]], cplus[[ri + s - j]]))))
      if (s != 1) cplus[[s]] <- .div(.sub(sg(ri + s), sm), .mul(list(S_Two, cplus[[ri]])))
    }
  }
  cminus <- lapply(cplus, .neg)
  cplus[[1L]] <- .div(.add(list(sg(ri + s), .neg(sm), .neg(.mul(list(.int(ri), cplus[[ri]]))))), .mul(list(S_Two, cplus[[ri]])))
  cminus[[1L]] <- .div(.add(list(sg(ri + s), .neg(sm), .neg(.mul(list(.int(ri), cminus[[ri]]))))), .mul(list(S_Two, cminus[[ri]])))
  if (!.ode_same_list(cplus, cminus)) return(list(cplus, cminus))
  cplus
}

# riccati.construct_c_case_3
.ode_construct_c_case_3 <- function() list(list(S_One))

# riccati.construct_c
.ode_construct_c <- function(num, den, x, poles, muls) {
  c <- list()
  for (i in seq_along(poles)) {
    mul <- .ode_num(muls[[i]])
    c[[i]] <- if (mul == 1) .ode_construct_c_case_3()
              else if (mul == 2) .ode_construct_c_case_1(num, den, x, poles[[i]])
              else .ode_construct_c_case_2(num, den, x, poles[[i]], mul)
  }
  c
}

# riccati.construct_d_case_4; the list index i + 1 holds d_i and the last
# entry d_(-1)
.ode_construct_d_case_4 <- function(ser, N) {
  N <- .ode_num(N)
  sg <- function(k) .ode_ser_get(ser, k)
  dplus <- rep(list(S_Zero), N + 2)
  dplus[[N + 1]] <- .sqrt(sg(2 * N))
  for (s in seq(N - 1, -1)) {
    sm <- S_Zero
    if (s + 1 <= N - 1) for (j in seq(s + 1, N - 1)) sm <- .add(list(sm, .mul(list(dplus[[j + 1]], dplus[[N + s - j + 1]]))))
    if (s != -1) dplus[[s + 1]] <- .div(.sub(sg(N + s), sm), .mul(list(S_Two, dplus[[N + 1]])))
  }
  dminus <- lapply(dplus, .neg)
  last <- N + 2
  dplus[[last]] <- .div(.add(list(sg(N + s), .neg(.mul(list(.int(N), dplus[[N + 1]]))), .neg(sm))), .mul(list(S_Two, dplus[[N + 1]])))
  dminus[[last]] <- .div(.add(list(sg(N + s), .neg(.mul(list(.int(N), dminus[[N + 1]]))), .neg(sm))), .mul(list(S_Two, dminus[[N + 1]])))
  if (!.ode_same_list(dplus, dminus)) return(list(dplus, dminus))
  dplus
}

# riccati.construct_d_case_5
.ode_construct_d_case_5 <- function(ser) {
  d0 <- .sqrt(.ode_ser_get(ser, 0))
  dplus <- list(d0, .div(.ode_ser_get(ser, -1), .mul(list(S_Two, d0))))
  dminus <- lapply(dplus, .neg)
  if (!.ode_same_list(dplus, dminus)) return(list(dplus, dminus))
  dplus
}

# riccati.construct_d_case_6
.ode_construct_d_case_6 <- function(num, den, x) {
  s_inf <- .ode_limit_at_inf(expand(.mul(list(.pow(x, S_Two), .ode_rx(num)))), den, x)
  if (!identical(s_inf$key, .rat(-1L, 4L)$key)) {
    s <- .sqrt(.add(list(S_One, .mul(list(.int(4), s_inf)))))
    return(list(list(.div(.add(list(S_One, s)), S_Two)), list(.div(.sub(S_One, s), S_Two))))
  }
  list(list(S_Half))
}

# riccati.construct_d
.ode_construct_d <- function(num, den, x, val_inf) {
  N <- (-val_inf) %/% 2
  mul <- if (val_inf < 0) -val_inf else 0
  ser <- .ode_rational_laurent_series(num, den, x, S_Infinity, mul, 1L)
  if (val_inf < 0) .ode_construct_d_case_4(ser, N)
  else if (val_inf == 0) .ode_construct_d_case_5(ser)
  else .ode_construct_d_case_6(num, den, x)
}

# riccati.rational_laurent_series: the coefficients of the Laurent series of
# num/den about r, keyed by the power of (x - r)
.ode_rational_laurent_series <- function(num, den, x, r, m, n) {
  num <- .ode_rx(num)
  den <- .ode_rx(den)
  r <- as_sym(r)
  m <- .ode_num(m)
  n <- .ode_num(n)
  if (identical(r$key, S_Infinity$key)) {
    nd <- .ode_inverse_transform(num, den, x)
    num <- nd[[1L]]
    den <- nd[[2L]]
    r <- S_Zero
  }
  if (!is_zero_num(r)) {
    num <- .ode_anp_poly(subs(num, x, .add(list(x, r))), x)
    den <- .ode_anp_poly(subs(den, x, .add(list(x, r))), x)
  }
  nd <- .ode_rcancel(.mul(list(num, .pow(x, .int(m)))), den, x)
  num <- .ode_anp_poly(nd[[1L]], x)
  den <- .ode_anp_poly(nd[[2L]], x)
  maxdegree <- 1 + max(.ode_rdeg(num, x), .ode_rdeg(den, x))
  ncs <- rev(.ode_rcoeffs(num, x))
  dcs <- rev(.ode_rcoeffs(den, x))
  nk <- function(k) if (k + 1 <= length(ncs)) ncs[[k + 1]] else S_Zero
  dk <- function(k) if (k + 1 <= length(dcs)) dcs[[k + 1]] else S_Zero
  # the triangular system num_k = sum(den_(k-i)*a_i), solved by linsolve
  coeffs <- list()
  for (k in seq_len(maxdegree) - 1L) {
    acc <- nk(k)
    if (k >= 1) for (i in seq_len(k) - 1L) acc <- .sub(acc, .mul(list(dk(k - i), coeffs[[i + 1L]])))
    coeffs[[k + 1L]] <- .ode_anp(.div(acc, dk(0)))
  }
  div <- dcs[[1L]]
  rec_rhs <- dcs[-1L]
  series <- coeffs
  while (length(series) < n) {
    L <- length(series)
    terms <- lapply(seq_along(rec_rhs), function(j) .mul(list(rec_rhs[[j]], series[[L - j + 1L]])))
    series[[L + 1L]] <- .neg(.div(.add(terms), div))
  }
  .pydict(lapply(seq_along(series) - 1L, function(i) .int(m - i)), series)
}

# a polynomial in x with its coefficients in canonical form
.ode_anp_poly <- function(e, x) {
  cs <- rev(.ode_rcoeffs(e, x))
  .add(lapply(seq_along(cs), function(i) .mul(list(.ode_anp(cs[[i]]), .pow(x, .int(i - 1L))))))
}

# riccati.compute_m_ybar
.ode_compute_m_ybar <- function(x, poles, choice, N) {
  last <- choice[[length(choice)]]
  m <- last[[length(last)]]
  dybar <- list()
  for (i in seq_along(poles)) {
    for (j in seq_along(choice[[i]])) dybar[[length(dybar) + 1L]] <- .div(choice[[i]][[j]], .pow(.sub(x, poles[[i]]), .int(j)))
    m <- .sub(m, choice[[i]][[1L]])
  }
  ybar <- .add(dybar)
  if (N >= 0) for (i in 0:N) ybar <- .add(list(ybar, .mul(list(last[[i + 1L]], .pow(x, .int(i))))))
  list(.ode_anp(m), ybar)
}

# riccati.solve_aux_eq
.ode_solve_aux_eq <- function(numa, dena, numy, deny, x, m) {
  psyms <- lapply(seq_len(m) - 1L, function(i) Dummy(paste0("C", i)))
  psol <- .pow(x, .int(m))
  for (k in seq_along(psyms)) psol <- .add(list(psol, .mul(list(psyms[[k]], .pow(x, .int(m - k))))))
  dx <- function(e) .diff1(e, x)
  auxeq <- .mul(list(.sub(.mul(list(dena, .add(list(.mul(list(dx(numy), deny)), .neg(.mul(list(numy, dx(deny)))), .pow(numy, S_Two))))),
                          .mul(list(numa, .pow(deny, S_Two)))), psol))
  if (m >= 1) {
    px <- dx(psol)
    auxeq <- .add(list(auxeq, .mul(list(px, S_Two, numy, deny, dena))))
  }
  if (m >= 2) auxeq <- .add(list(auxeq, .mul(list(dx(px), .pow(deny, S_Two), dena))))
  auxeq <- expand(auxeq)
  if (m != 0) {
    eqs <- Filter(function(e) !is_zero_num(e), .ode_rcoeffs(auxeq, x))
    coeffs <- NULL
    if (!length(eqs)) coeffs <- lapply(psyms, function(s) list(s, s))
    else {
      sl <- tryCatch(.sp_solve(eqs, psyms, list(dict = TRUE, `_raw` = TRUE)), error = function(e) list())
      if (length(sl)) {
        d <- sl[[1L]]
        coeffs <- lapply(psyms, function(s) {
          j <- match(s$key, vapply(d$k, .k, ""))
          list(s, if (is.na(j)) s else d$v[[j]])
        })
      }
    }
    if (is.null(coeffs)) return(list(psol, list(), TRUE, FALSE))
    pk <- vapply(psyms, .k, "")
    is_general <- any(vapply(coeffs, function(p) any(vapply(free_symbols(p[[2L]]), .k, "") %in% pk), TRUE))
    return(list(psol, coeffs, TRUE, is_general))
  }
  list(S_One, auxeq, is_zero_num(auxeq), FALSE)
}

# riccati.remove_redundant_sols
.ode_remove_redundant_sols <- function(sol1, sol2, x) {
  nsym <- function(e) {
    a <- atoms(e, c("Symbol", "Dummy"))
    length(unique(c(vapply(a, .k, ""), x$key)))
  }
  nd1 <- .as_numer_denom(together(sol1))
  nd2 <- .as_numer_denom(together(sol2))
  e <- expand(.sub(.mul(list(nd1[[1L]], nd2[[2L]])), .mul(list(nd1[[2L]], nd2[[1L]]))))
  cs <- Filter(function(z) !is_zero_num(z), .ode_rcoeffs(e, x))
  syms <- Filter(function(s) !identical(s$key, x$key), atoms(e, c("Symbol", "Dummy")))
  ok <- !length(cs) ||
    (length(syms) && length(tryCatch(.sp_solve(cs, syms, list(dict = TRUE, `_raw` = TRUE)), error = function(err) list())) > 0L)
  if (ok) {
    n1 <- nsym(sol1)
    n2 <- nsym(sol2)
    if (n1 > n2) return(sol2)
    return(sol1)
  }
  NULL
}

# riccati.get_gen_sol_from_part_sol
.ode_get_gen_sol_from_part_sol <- function(part_sols, a, x) {
  if (!length(part_sols)) return(list())
  if (length(part_sols) == 1L) {
    y1 <- part_sols[[1L]]
    i <- fcall("exp", list(Integral(.mul(list(S_Two, y1)), x)))
    z <- .mul(list(i, .add(list(Integral(.inv(i), x), Dummy("C1")))))
    z <- doit(z)
    if (is_zero_num(z)) return(y1)
    return(.add(list(y1, .inv(z))))
  }
  if (length(part_sols) == 2L) {
    y1 <- part_sols[[1L]]
    y2 <- part_sols[[2L]]
    C1 <- Dummy("C1")
    u <- .mul(list(C1, doit(fcall("exp", list(Integral(.sub(y2, y1), x))))))
    if (identical(u$key, S_One$key)) return(y2)
    return(.div(.sub(.mul(list(y2, u)), y1), .sub(u, S_One)))
  }
  y1 <- part_sols[[1L]]
  y2 <- part_sols[[2L]]
  y3 <- part_sols[[3L]]
  C1 <- Dummy("C1")
  .div(.mul(list(.add(list(C1, S_One)), y2, .sub(y1, y3))),
       .add(list(.mul(list(C1, y1)), y2, .neg(.mul(list(.add(list(C1, S_One)), y3))))))
}

# itertools.product over a list of option lists
.ode_product <- function(opts) {
  out <- list(list())
  for (o in opts) {
    nxt <- list()
    for (p in out) for (v in o) nxt[[length(nxt) + 1L]] <- c(p, list(v))
    out <- nxt
  }
  out
}

# riccati.solve_riccati
.ode_solve_riccati <- function(fx, x, b0, b1, b2, gensol = FALSE) {
  a <- .ode_riccati_a(x, b0, b1, b2)
  a_t <- together(a)
  nd0 <- .as_numer_denom(a_t)
  nd <- .ode_rcancel(nd0[[1L]], nd0[[2L]], x)
  num <- nd[[1L]]
  den <- nd[[2L]]
  presol <- list()
  if (is_zero_num(num)) {
    presol[[1L]] <- list(.inv(.add(list(x, Dummy("C1")))), TRUE)
  } else if (.ode_rdeg(num, x) <= 0 && .ode_rdeg(den, x) <= 0) {
    presol <- list(list(.sqrt(a), FALSE), list(.neg(.sqrt(a)), FALSE))
  }
  pr <- roots(den, x)
  mult <- attr(pr, "multiplicity")
  poles <- unclass(pr)
  attributes(poles) <- NULL
  val_inf <- .ode_val_at_inf(num, den, x)
  if (length(poles)) {
    if (!.ode_check_necessary_conds(val_inf, as.list(mult))) .ode_stop("ValueError", "Rational Solution doesn't exist")
    c <- .ode_construct_c(num, den, x, poles, as.list(mult))
    d <- .ode_construct_d(num, den, x, val_inf)
    c[[length(c) + 1L]] <- d
    for (choice in .ode_product(c)) {
      mb <- .ode_compute_m_ybar(x, poles, choice, (-val_inf) %/% 2)
      m <- mb[[1L]]
      ybar <- mb[[2L]]
      ndy <- .as_numer_denom(together(ybar))
      numy <- expand(ndy[[1L]])
      deny <- expand(ndy[[2L]])
      if (isTRUE(ask_fact(m, "nonnegative")) && isTRUE(ask_fact(m, "integer"))) {
        ax <- .ode_solve_aux_eq(num, den, numy, deny, x, as.integer(.ode_num(m)))
        psol <- ax[[1L]]
        coeffs <- ax[[2L]]
        if (isTRUE(ax[[3L]])) {
          if (identical(psol$key, S_One$key) && is_sym(coeffs) && is_zero_num(coeffs)) {
            presol[[length(presol) + 1L]] <- list(ybar, FALSE)
          } else if (is.list(coeffs) && !is_sym(coeffs) && length(coeffs)) {
            psol <- xreplace(psol, coeffs)
            sol <- .add(list(ybar, .div(.diff1(psol, x), psol)))
            presol[[length(presol) + 1L]] <- list(sol, ax[[4L]])
          }
        }
      }
    }
  }
  remove <- character()
  if (length(presol) >= 2L) {
    for (i in seq_len(length(presol) - 1L)) for (j in seq(i + 1L, length(presol))) {
      rem <- .ode_remove_redundant_sols(presol[[i]][[1L]], presol[[j]][[1L]], x)
      if (!is.null(rem)) remove <- c(remove, rem$key)
    }
  }
  sols <- Filter(function(s) !(s[[1L]]$key %in% remove), presol)
  bp <- .sub(.div(.neg(.diff1(b2, x)), .mul(list(S_Two, .pow(b2, S_Two)))), .div(b1, .mul(list(S_Two, b2))))
  if (gensol) {
    gen <- NULL
    for (s in sols) if (isTRUE(s[[2L]])) {
      gen <- list(s[[1L]])
      break
    }
    sols <- if (!is.null(gen)) gen else list(.ode_get_gen_sol_from_part_sol(lapply(sols, `[[`, 1L), a, x))
  } else sols <- lapply(sols, `[[`, 1L)
  lapply(sols, function(y) Eq(fx, cancel(.ode_riccati_inverse_normal(y, x, b1, b2, bp), extension = TRUE)))
}
