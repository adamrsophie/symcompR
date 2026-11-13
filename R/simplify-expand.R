# Expansion, after Expr.expand() in sympy/core/expr.py and sympy/core/function.py.

#' Expand an expression
#'
#' Distributes products over sums and expands powers of sums, like SymPy's
#' `expand()`. Hints switch individual rules: `mul` (distribute), `multinomial`
#' (powers of sums), `power_exp` (`x^(a+b) = x^a*x^b`), `power_base`
#' (`(x*y)^a = x^a*y^a` when valid), `log` (`log(x*y) = log(x) + log(y)` when
#' valid), `trig` (angle addition), `complex` (split into real and imaginary
#' parts) and `func` (special function identities). `force = TRUE` ignores
#' the validity assumptions of `power_base` and `log`.
#'
#' The helpers `expand_mul()`, `expand_multinomial()`, `expand_log()`,
#' `expand_power_exp()`, `expand_power_base()`, `expand_trig()`,
#' `expand_complex()` and `expand_func()` apply one rule.
#'
#' @param e A `sym` expression, matrix or list.
#' @param deep Expand inside function arguments too.
#' @param mul,multinomial,power_exp,power_base,log,trig,complex,func Hints.
#' @param force Apply `log` and `power_base` regardless of assumptions.
#' @param basic The `basic` hint, as in SymPy.
#' @param factor For the `log` hint: factor integer arguments of logarithms.
#' @param ... Ignored.
#' @return A `sym` expression.
#' @examples
#' x <- Symbol("x"); y <- Symbol("y")
#' expand((x + y)^3)
#' expand((x + 1) * (x - 1))
#' expand_trig(sin(x + y))
#' @export
expand <- function(e, deep = TRUE, mul = TRUE, multinomial = TRUE,
                   power_exp = TRUE, power_base = TRUE, log = TRUE,
                   trig = FALSE, complex = FALSE, func = FALSE,
                   force = FALSE, basic = TRUE, factor = FALSE, ...) {
  if (inherits(e, "TensExpr")) return(.tn_expand(e))
  if (inherits(e, "symmatrix")) return(.mat_map(e, function(z)
    expand(z, deep, mul, multinomial, power_exp, power_base, log, trig, complex, func, force)))
  if (is.list(e) && !inherits(e, "sym")) return(lapply(e, expand, deep = deep, mul = mul,
                                                        multinomial = multinomial, power_exp = power_exp,
                                                        power_base = power_base, log = log, trig = trig,
                                                        complex = complex, func = func, force = force))
  e <- as_sym(e)
  # Expectation, Variance and Covariance expand by linearity (sympy.stats)
  if (e$type %in% c("Expectation", "Variance", "Covariance")) return(.stats_expand(e))
  if (.is_bd(e)) return(.bd_expand(e, deep = deep, mul = mul, multinomial = multinomial, power_exp = power_exp,
                                  power_base = power_base, log = log, trig = trig, complex = complex, func = func,
                                  force = force, basic = basic, factor = factor))
  h <- list(mul = mul, multinomial = multinomial, power_exp = power_exp, power_base = power_base, log = log,
            trig = trig, complex = complex, func = func, force = force, basic = basic, factor = factor)
  # Expr.expand: one pass per hint in sorted order (mul after multinomial),
  # then multinomial, mul and log to a fixed point
  for (hint in c("basic", "complex", "func", "log", "multinomial", "mul", "power_base", "power_exp", "trig"))
    if (isTRUE(h[[hint]])) e <- .expand_hint(e, hint, deep, h)
  repeat {
    was <- e$key
    if (multinomial) e <- .expand_hint(e, "multinomial", deep, h)
    if (mul) e <- .expand_hint(e, "mul", deep, h)
    if (log) e <- .expand_hint(e, "log", deep, h)
    if (identical(e$key, was)) break
  }
  e
}

# Expr._expand_hint: the arguments first (always deep below the top), then
# the node's _eval_expand_<hint>
.expand_hint <- function(e, hint, deep = TRUE, h) {
  if (deep && length(e$args)) {
    args <- lapply(e$args, .expand_hint, hint = hint, deep = TRUE, h = h)
    if (!identical(vapply(args, .k, ""), vapply(e$args, .k, ""))) e <- sym_rebuild(e, args)
  }
  r <- .expand_node(e, hint, h)
  if (is.null(r)) e else r
}

.all_nonneg_or_nonpos <- function(a) {
  args <- if (a$type == "Add") a$args else list(a)
  all(vapply(args, function(z) isTRUE(ask_fact(z, "nonnegative")), NA)) || all(vapply(args, function(z) isTRUE(ask_fact(z, "nonpositive")), NA))
}

.expand_node <- function(e, hint, h) {
  t <- e$type
  switch(hint,
    basic = {
      # AddWithLimits._eval_expand_basic
      if (t %in% c("Integral", "Sum", "Product")) {
        summand <- do.call(expand, c(list(e$args[[1L]]), h[c("mul", "multinomial", "power_exp", "power_base", "log", "trig", "complex", "func", "force")]))
        fin <- !any(vapply(e$args[-1L], function(l) any(vapply(l$args[-1L], function(z) isTRUE(ask_fact(z, "infinite")), NA)), NA))
        if (summand$type == "Add" && (h$force || (.is_comm(summand) && fin)))
          return(.add(lapply(summand$args, function(a) new_sym(t, c(list(a), e$args[-1L])))))
        if (!identical(summand$key, e$args[[1L]]$key)) return(new_sym(t, c(list(summand), e$args[-1L])))
      }
      NULL
    },
    complex = if (t %in% c("Symbol", "Dummy", "Add", "Mul", "Pow", "Function")) .expand_complex1(e) else NULL,
    func = if (t == "Function") { fh <- .expand_hooks[[e$name]]; if (is.null(fh)) NULL else fh(e) } else NULL,
    log = if (t == "Function" && e$name == "log" && length(e$args) == 1L) .expand_log1(e$args[[1L]], h$force, h$factor) else NULL,
    multinomial = if (t == "Pow") .expand_multinomial1(e, h) else NULL,
    mul = if (t == "Mul") .expand_mul1(e) else NULL,
    power_base = if (t == "Pow" && e$args[[1L]]$type == "Mul") .eval_expand_power_base(e, h$force) else NULL,
    power_exp = {
      if (t == "Function" && e$name == "exp" && e$args[[1L]]$type == "Add" && .is_comm(e$args[[1L]]))
        return(.mul(lapply(e$args[[1L]]$args, function(q) fcall("exp", list(q)))))
      if (t == "Pow") {
        b <- e$args[[1L]]; p <- e$args[[2L]]
        if (p$type == "Add" && (h$force || identical(ask_fact(b, "zero"), FALSE) || .all_nonneg_or_nonpos(p))) {
          if (.is_comm(p)) return(.mul(lapply(p$args, function(q) .pow(b, q))))
          if (.is_comm(b)) {
            cm <- vapply(p$args, .is_comm, NA)
            if (any(cm)) return(.mul(c(lapply(p$args[cm], function(q) .pow(b, q)), list(.pow(b, new_sym("Add", p$args[!cm]))))))
          }
        }
      }
      NULL
    },
    trig = if (t == "Function" && e$name %in% c("sin", "cos", "tan", "cot", "sec", "csc", "sinh", "cosh", "tanh", "coth", "sech", "csch")) .expand_trig1(e, h) else NULL,
    NULL)
}

# Mul._eval_expand_mul (hints carry no `deep`, so new products of sums are
# left to the fixed point loop of expand)
.expand_mul1 <- function(e) {
  nd <- fraction(e)
  n <- nd[[1L]]; d <- nd[[2L]]
  if (d$type == "Mul") {
    if (n$type == "Mul") n <- .expand_mul1(n)
    d <- .expand_mul1(d)
  }
  e <- .div(n, d)
  if (e$type != "Mul") return(e)
  fs <- e$args
  is_sum <- vapply(fs, function(f) f$type == "Add", NA)
  if (!any(is_sum) && all(vapply(fs, .is_comm, NA))) return(e)
  plain <- .mul(fs[!is_sum & vapply(fs, .is_comm, NA)])
  sums <- fs[is_sum | !vapply(fs, .is_comm, NA)]
  if (!any(is_sum)) return(e)
  # Mul._expandsums
  terms <- list(list())
  for (f in sums) {
    parts <- if (f$type == "Add") f$args else list(f)
    terms <- unlist(lapply(terms, function(tl) lapply(parts, function(p) c(tl, list(p)))), recursive = FALSE)
  }
  .add(lapply(terms, function(tl) if (all(vapply(tl, .is_comm, NA))) .mul(c(list(plain), tl)) else new_sym("Mul", c(list(plain), tl))))
}

# Pow._eval_expand_multinomial
.expand_multinomial1 <- function(e, h) {
  b <- e$args[[1L]]; p <- e$args[[2L]]
  rat <- p$type %in% c("Integer", "Rational")
  if (rat && p$value > 0 && b$type == "Add") {
    if (p$type == "Rational") {
      q <- gmp::as.bigq(p$value)
      n <- gmp::as.bigz(gmp::numerator(q) %/% gmp::denominator(q))
      if (n == 0) return(NULL)
      radical <- .pow(b, make_rational(q - gmp::as.bigq(n)))
      ebn <- .pow(b, .zs(n))
      if (ebn$type == "Pow") { r <- .expand_multinomial1(ebn, h); if (!is.null(r)) ebn <- r }
      return(.add(lapply(if (ebn$type == "Add") ebn$args else list(ebn), function(tm) .mul(list(tm, radical)))))
    }
    return(.multinomial_expand(b, as.numeric(p$value)))
  }
  if (rat && p$value < 0 && b$type == "Add") {
    q <- gmp::as.bigq(p$value)
    if (abs(gmp::numerator(q)) > gmp::denominator(q)) {
      inner <- .pow(b, .neg(p))
      r <- if (inner$type == "Pow") .expand_multinomial1(inner, h) else NULL
      return(.inv(if (is.null(r)) inner else r))
    }
    return(NULL)
  }
  if (p$type == "Add" && is_number_type(b) && (h$force || identical(ask_fact(b, "zero"), FALSE) || .all_nonneg_or_nonpos(p))) {
    num <- vapply(p$args, is_number_type, NA)
    return(.mul(c(lapply(p$args[num], function(z) .pow(b, z)), list(.pow(b, .add(p$args[!num]))))))
  }
  NULL
}

.expand_complex1 <- function(e) {
  ri <- .as_real_imag(e)
  .add(list(ri[[1L]], .mul(list(S_I, ri[[2L]]))))
}

# Distribute a product over the sums among its factors.
.expand_mul_dist <- function(e, h = NULL) {
  if (e$type == "Add") {
    return(.add(lapply(e$args, function(t) if (t$type == "Mul") .expand_mul_dist(t, h) else t)))
  }
  if (e$type != "Mul") return(e)
  fs <- e$args
  is_sum <- vapply(fs, function(f) f$type == "Add", TRUE)
  if (!any(is_sum)) return(e)
  comm <- vapply(fs, .is_comm, TRUE)
  if (!all(comm)) {
    terms <- list(list())
    for (f in fs) {
      parts <- if (f$type == "Add") f$args else list(f)
      terms <- unlist(lapply(terms, function(tl) lapply(parts, function(p) c(tl, list(p)))),
                      recursive = FALSE)
    }
    return(.add(lapply(terms, .mul)))
  }
  rest <- .mul(fs[!is_sum])
  acc <- list(rest)
  for (f in fs[is_sum]) {
    nxt <- vector("list", length(acc) * length(f$args))
    k <- 0L
    for (a in acc) for (b in f$args) {
      k <- k + 1L
      nxt[[k]] <- .mul(list(a, b))
    }
    s <- .add(nxt)
    acc <- if (s$type == "Add") s$args else list(s)
  }
  res <- .add(lapply(acc, function(t) {
    if (t$type == "Mul" && any(vapply(t$args, function(f) f$type == "Add", TRUE)))
      .expand_mul_dist(t, h) else t
  }))
  res
}

# (t1 + ... + tk)^n by the multinomial theorem.
.multinomial_expand <- function(b, n) {
  terms <- b$args
  k <- length(terms)
  if (n == 1) return(b)
  if (!all(vapply(terms, .is_comm, TRUE))) {
    # noncommutative terms: distribute ordered products
    # ponytail: k^n products, fine for the small powers this is used with
    r <- b
    for (i in seq_len(n - 1)) r <- .add(unlist(lapply(if (r$type == "Add") r$args else list(r), function(u)
      lapply(terms, function(t) .mul(list(u, t)))), recursive = FALSE))
    return(r)
  }
  if (n > 2000) stop("exponent too large to expand")
  fact <- gmp::factorialZ(0:n)
  out <- list()
  rec <- function(i, left, exps) {
    if (i == k) {
      exps <- c(exps, left)
      coef <- fact[n + 1]
      for (ev in exps) coef <- gmp::divq.bigz(coef, fact[ev + 1])
      parts <- list(make_rational(gmp::as.bigq(coef)))
      for (j in seq_len(k)) if (exps[j] > 0) parts[[length(parts) + 1L]] <- .pow(terms[[j]], .int(exps[j]))
      out[[length(out) + 1L]] <<- .mul(parts)
      return(invisible())
    }
    for (e in left:0) rec(i + 1L, left - e, c(exps, e))
  }
  rec(1L, n, integer())
  .add(out)
}

# log._eval_expand_log, iterated to a fixed point as Expr.expand does; NULL
# when log(a) does not change
.expand_log1 <- function(a, force, factor = FALSE) {
  full <- function(z) {
    lz <- fcall("log", list(z))
    if (!(lz$type == "Function" && lz$name == "log")) return(lz)
    r <- .expand_log1(z, force, factor)
    if (is.null(r)) lz else r
  }
  if (a$type == "Integer") {
    pp <- if (abs(a$value) > 1) .perfect_power(gmp::as.bigz(a$value)) else FALSE
    arg <- a; coeff <- S_One; logarg <- NULL
    if (!isFALSE(pp)) { arg <- .int(pp[[1L]]); coeff <- .int(pp[[2L]]); logarg <- fcall("log", list(arg)) }
    if (isTRUE(factor) && abs(arg$value) > 1) {
      fi <- factorint(arg)
      ks <- names(fi)
      if (!(as.character(arg$value) %in% ks))
        logarg <- .add(lapply(seq_along(fi), function(i) .mul(list(.int(fi[[i]]), fcall("log", list(as_sym(gmp::as.bigz(ks[i]))))))))
    }
    if (!is.null(logarg)) return(.mul(list(coeff, logarg)))
    return(NULL)
  }
  if (a$type == "Rational") {
    q <- gmp::as.bigq(a$value)
    return(.sub(full(.zs(gmp::numerator(q))), full(.zs(gmp::denominator(q)))))
  }
  if (a$type == "Mul") {
    ex <- list(); nonpos <- list()
    for (f in .sympy_args(a)) {
      if (force || isTRUE(ask_fact(f, "positive")) || isTRUE(.is_polar(f))) ex[[length(ex) + 1L]] <- full(f)
      else if (isTRUE(ask_fact(f, "negative"))) { ex[[length(ex) + 1L]] <- fcall("log", list(.neg(f))); nonpos[[length(nonpos) + 1L]] <- S_NegativeOne }
      else nonpos[[length(nonpos) + 1L]] <- f
    }
    r <- .add(c(ex, list(fcall("log", list(.mul(nonpos))))))
    if (identical(r$key, fcall("log", list(a))$key)) return(NULL)
    return(r)
  }
  isexp <- a$type == "Function" && a$name == "exp"
  if (a$type == "Pow" || isexp) {
    b <- if (isexp) S_E else a$args[[1L]]
    p <- if (isexp) a$args[[1L]] else a$args[[2L]]
    if (force || (isTRUE(ask_fact(p, "extended_real")) && (isTRUE(ask_fact(b, "positive")) ||
          (isTRUE(ask_fact(.add(list(p, S_One)), "positive")) && isTRUE(ask_fact(.sub(p, S_One), "nonpositive"))))) ||
        isTRUE(.is_polar(b))) {
      return(.mul(list(unpolarify(p), full(b))))
    }
  }
  NULL
}

# _eval_expand_trig of the trigonometric and hyperbolic functions; h are
# the hints of the expand() call (sinh and cosh expand their argument)
.expand_trig1 <- function(e, h = NULL) {
  a <- e$args[[1L]]
  nm <- e$name
  cm <- function(z) {
    if (z$type == "Mul" && z$args[[1L]]$type %in% c("Integer", "Rational")) return(list(z$args[[1L]], .mul(z$args[-1L])))
    if (z$type %in% c("Integer", "Rational")) return(list(z, S_One))
    list(S_One, z)
  }
  two <- function(z) { ar <- .sympy_args(z); list(ar[[1L]], if (length(ar) == 2L) ar[[2L]] else .add(ar[-1L])) }
  et <- function(name, z) .expand_trig1(fnode(name, list(z)), h)
  F <- function(name, ...) fcall(name, list(...))
  switch(nm,
    sin = {
      if (a$type == "Add") {
        xy <- two(a)
        return(.add(list(.mul(list(et("sin", xy[[1L]]), et("cos", xy[[2L]]))), .mul(list(et("sin", xy[[2L]]), et("cos", xy[[1L]]))))))
      }
      if (a$type == "Mul") {
        nx <- cm(a); n <- nx[[1L]]; x <- nx[[2L]]
        if (n$type == "Integer") {
          ni <- as.numeric(n$value)
          if (ni %% 2 != 0) return(.mul(list(.pow(S_NegativeOne, .int((ni - 1) / 2)), F("chebyshevt", n, sin(x)))))
          return(expand_mul(.mul(list(.pow(S_NegativeOne, .int(ni / 2 - 1)), cos(x), F("chebyshevu", .int(ni - 1), sin(x)))), deep = FALSE))
        }
      }
      sin(a)
    },
    cos = {
      if (a$type == "Add") {
        xy <- two(a)
        return(.sub(.mul(list(et("cos", xy[[1L]]), et("cos", xy[[2L]]))), .mul(list(et("sin", xy[[1L]]), et("sin", xy[[2L]])))))
      }
      if (a$type == "Mul") {
        nx <- cm(a)
        if (nx[[1L]]$type == "Integer") return(F("chebyshevt", nx[[1L]], cos(nx[[2L]])))
      }
      cos(a)
    },
    tan = , cot = {
      if (a$type == "Add") {
        n <- length(a$args)
        TX <- lapply(.sympy_args(a), function(x) et(nm, x))
        Y <- lapply(seq_len(n) - 1L, function(i) Symbol(paste0("Y", i)))
        p <- list(S_Zero, S_Zero)
        if (nm == "tan") for (i in 0:n) {
          k <- 2L - i %% 2L
          p[[k]] <- .add(list(p[[k]], .mul(list(symmetric_poly(i, Y), .int((-1)^((i %% 4L) %/% 2L))))))
        } else for (i in n:0) {
          k <- (n - i) %% 2L + 1L
          p[[k]] <- .add(list(p[[k]], .mul(list(symmetric_poly(i, Y), .int((-1)^(((n - i) %% 4L) %/% 2L))))))
        }
        r <- .div(p[[1L]], p[[2L]])
        for (i in seq_len(n)) r <- subs(r, Y[[i]], TX[[i]])
        return(r)
      }
      if (a$type == "Mul") {
        ct <- cm(a)
        if (ct[[1L]]$type == "Integer" && ct[[1L]]$value > 1) {
          z <- Symbol("dummy", real = TRUE)
          P <- expand(.pow(if (nm == "tan") .add(list(S_One, .mul(list(S_I, z)))) else .add(list(z, S_I)), ct[[1L]]))
          r <- if (nm == "tan") .div(F("im", P), F("re", P)) else .div(F("re", P), F("im", P))
          return(subs(r, z, F(nm, ct[[2L]])))
        }
      }
      F(nm, a)
    },
    sec = , csc = {
      t <- F(if (nm == "sec") "cos" else "sin", a)
      if (!(t$type == "Function" && t$name %in% c("cos", "sin"))) return(NULL)
      .inv(.expand_trig1(t, h))
    },
    sinh = , cosh = {
      arg <- if (is.null(h)) expand(a) else do.call(expand, c(list(a), h[c("mul", "multinomial", "power_exp", "power_base", "log", "trig", "complex", "func", "force")]))
      x <- NULL
      if (arg$type == "Add") { xy <- two(arg); x <- xy[[1L]]; y <- xy[[2L]] }
      else {
        ct <- cm(arg)
        if (!is_one_num(ct[[1L]]) && ct[[1L]]$type == "Integer" && !is_one_num(ct[[2L]])) { x <- ct[[2L]]; y <- .mul(list(.sub(ct[[1L]], S_One), x)) }
      }
      if (!is.null(x)) {
        r <- if (nm == "sinh") .add(list(.mul(list(sinh(x), cosh(y))), .mul(list(sinh(y), cosh(x)))))
             else .add(list(.mul(list(cosh(x), cosh(y))), .mul(list(sinh(x), sinh(y)))))
        return(expand(r, trig = TRUE))
      }
      F(nm, arg)
    },
    tanh = {
      if (a$type == "Add") {
        n <- length(a$args)
        TX <- lapply(.sympy_args(a), function(x) et("tanh", x))
        p <- list(S_Zero, S_Zero)
        for (i in 0:n) p[[i %% 2L + 1L]] <- .add(list(p[[i %% 2L + 1L]], symmetric_poly(i, TX)))
        return(.div(p[[2L]], p[[1L]]))
      }
      if (a$type == "Mul") {
        ct <- if (is_number_type(a$args[[1L]])) list(a$args[[1L]], .mul(a$args[-1L])) else list(S_One, a)
        if (ct[[1L]]$type == "Integer" && ct[[1L]]$value > 1) {
          cn <- as.integer(.int_val(ct[[1L]])); T <- fcall("tanh", list(ct[[2L]]))
          nn <- lapply(seq(1L, cn, by = 2L), function(k) .mul(list(.int(choose(cn, k)), .pow(T, .int(k)))))
          dd <- lapply(seq(0L, cn, by = 2L), function(k) .mul(list(.int(choose(cn, k)), .pow(T, .int(k)))))
          return(.div(.add(nn), .add(dd)))
        }
      }
      fcall("tanh", list(a))
    },
    coth = {
      if (a$type == "Add") {
        CX <- lapply(.sympy_args(a), function(x) et("coth", x))
        n <- length(a$args); p <- list(list(), list())
        for (i in n:0) p[[(n - i) %% 2L + 1L]][[length(p[[(n - i) %% 2L + 1L]]) + 1L]] <- symmetric_poly(i, CX)
        return(.div(.add(p[[1L]]), .add(p[[2L]])))
      }
      if (a$type == "Mul") {
        ct <- cm(a)
        if (ct[[1L]]$type == "Integer" && ct[[1L]]$value > 1) {
          cn <- as.integer(.int_val(ct[[1L]])); c <- fnode("coth", list(ct[[2L]])); p <- list(list(), list())
          for (i in cn:0) p[[(cn - i) %% 2L + 1L]][[length(p[[(cn - i) %% 2L + 1L]]) + 1L]] <- .mul(list(binomial(.int(cn), .int(i)), .pow(c, .int(i))))
          return(.div(.add(p[[1L]]), .add(p[[2L]])))
        }
      }
      fcall("coth", list(a))
    },
    sech = , csch = {
      t <- F(if (nm == "sech") "cosh" else "sinh", a)
      if (!(t$type == "Function" && t$name %in% c("cosh", "sinh"))) return(NULL)
      .inv(.expand_trig1(t, h))
    },
    NULL)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# Pow._eval_expand_power_base
.eval_expand_power_base <- function(e, force = FALSE) {
  b <- e$args[[1L]]; ex <- e$args[[2L]]
  if (b$type != "Mul") return(e)
  comm <- vapply(b$args, .is_comm, TRUE)
  cargs <- b$args[comm]; nc <- b$args[!comm]
  if (length(nc)) {
    nc <- lapply(nc, function(i) if (i$type == "Pow") .eval_expand_power_base(i, force) else i)
    if (ex$type == "Integer") {
      n <- as.integer(.int_val(ex))
      rv <- if (n > 0) new_sym("Mul", rep(nc, n)) else new_sym("Mul", rep(lapply(rev(nc), function(i) .pow(i, S_NegativeOne)), -n))
      if (length(cargs)) rv <- .mul(list(rv, .pow(.mul(cargs), ex)))
      return(rv)
    }
    if (!length(cargs)) return(new_sym("Pow", list(.mul(nc), ex)))
    nc <- list(.mul(nc))
  }
  other <- list(); maybe_real <- list()
  for (x in cargs) {
    if (identical(ask_fact(x, "extended_real"), FALSE) && !isTRUE(ask_fact(x, "polar"))) other[[length(other) + 1L]] <- x
    else maybe_real[[length(maybe_real) + 1L]] <- x
  }
  nonneg <- list(); neg <- list(); imag <- list()
  for (x in maybe_real) {
    if (is_const(x, "I")) { imag[[length(imag) + 1L]] <- x; next }
    if (isTRUE(ask_fact(x, "polar"))) { nonneg[[length(nonneg) + 1L]] <- x; next }
    nn <- ask_fact(x, "extended_nonnegative")
    if (isTRUE(nn)) nonneg[[length(nonneg) + 1L]] <- x
    else if (identical(nn, FALSE)) neg[[length(neg) + 1L]] <- x
    else other[[length(other) + 1L]] <- x
  }
  if (length(imag)) {
    i <- length(imag) %% 4L
    if (i == 1L) other[[length(other) + 1L]] <- S_I
    else if (i >= 2L) {
      if (length(neg)) {
        nonn <- .neg(neg[[length(neg)]]); neg <- neg[-length(neg)]
        if (!is_one_num(nonn)) nonneg[[length(nonneg) + 1L]] <- nonn
      } else neg[[length(neg) + 1L]] <- S_NegativeOne
      if (i == 3L) other[[length(other) + 1L]] <- S_I
    }
  }
  if (force || isTRUE(ask_fact(ex, "integer"))) {
    cargs <- c(nonneg, neg, other)
    other <- nc
  } else {
    if (length(neg) > 1L) {
      o <- S_One
      if (!length(other) && is_number_type(neg[[1L]])) { o <- .mul(list(o, neg[[1L]])); neg <- neg[-1L] }
      if (length(neg) %% 2L == 1L) o <- .neg(o)
      for (n in neg) nonneg[[length(nonneg) + 1L]] <- .neg(n)
      if (!is_one_num(o)) other[[length(other) + 1L]] <- o
    } else if (length(neg) && length(other)) {
      if (is_number_type(neg[[1L]]) && !is_negone_num(neg[[1L]])) {
        other[[length(other) + 1L]] <- S_NegativeOne
        nonneg[[length(nonneg) + 1L]] <- .neg(neg[[1L]])
      } else other <- c(other, neg)
    } else other <- c(other, neg)
    cargs <- nonneg
    other <- c(other, nc)
  }
  rv <- S_One
  if (length(cargs)) {
    if (ex$type %in% c("Integer", "Rational")) {
      isnp <- vapply(cargs, function(x) x$type == "Pow" && x$args[[2L]]$type %in% c("Integer", "Rational") && !length(free_symbols(x$args[[1L]])), TRUE)
      rv <- .mul(lapply(cargs[isnp], function(x) .pow(sym_rebuild(x, x$args), ex)))
      cargs <- cargs[!isnp]
    }
    rv <- .mul(list(rv, .mul(lapply(cargs, function(x) .pow(x, ex)))))
  }
  if (length(other)) rv <- .mul(list(rv, .pow(.mul(other), ex)))
  rv
}

.hint_expand <- function(on) {
  force(on)
  function(e, deep = TRUE, force = FALSE) {
    args <- list(e, deep = deep, mul = FALSE, multinomial = FALSE, power_exp = FALSE, power_base = FALSE,
                 log = FALSE, trig = FALSE, complex = FALSE, func = FALSE, force = force, basic = FALSE)
    args[on] <- TRUE
    do.call(expand, args)
  }
}

#' @rdname expand
#' @export
expand_mul <- .hint_expand("mul")
#' @rdname expand
#' @export
expand_multinomial <- .hint_expand("multinomial")
#' @rdname expand
#' @export
expand_log <- function(e, deep = TRUE, force = FALSE, factor = FALSE) {
  if (inherits(e, "symmatrix")) return(.mat_map(e, function(z) expand_log(z, deep, force, factor)))
  e <- as_sym(e)
  if (!factor) {
    islog_int <- function(z) z$type == "Function" && z$name == "log" && z$args[[1L]]$type == "Integer"
    handle_mul <- function(x) {
      nd <- fraction(x)
      nl <- Filter(islog_int, .uniq_syms(.all_subexprs(nd[[1L]])))
      dl <- Filter(islog_int, .uniq_syms(.all_subexprs(nd[[2L]])))
      if (length(nl) == 1L && length(dl) == 1L) {
        n <- nl[[1L]]; d <- dl[[1L]]
        m <- as.integer(multiplicity(d$args[[1L]], n$args[[1L]]))
        if (m) {
          q <- .zs(gmp::as.bigz(n$args[[1L]]$value) %/% gmp::as.bigz(d$args[[1L]]$value)^m)
          r <- .add(list(.int(m), .div(fcall("log", list(q)), d)))
          x <- subs(x, n, .mul(list(d, r)))
        }
      }
      x1 <- expand_mul(expand_log(x, deep = deep, force = force, factor = TRUE))
      if (.mi_count_fn(x1, "log") <= .mi_count_fn(x, "log")) x1 else x
    }
    rat_log <- function(j) any(vapply(if (j$type == "Add") j$args else list(j), function(i) i$type == "Function" && i$name == "log" && i$args[[1L]]$type %in% c("Integer", "Rational"), NA))
    e <- .map_bottom_up(e, function(z) if (z$type == "Mul" && all(vapply(.as_numer_denom(z), rat_log, NA))) handle_mul(z) else z)
  }
  expand(e, deep = deep, log = TRUE, mul = FALSE, power_exp = FALSE, power_base = FALSE, multinomial = FALSE,
         basic = FALSE, force = force, factor = factor)
}
#' @rdname expand
#' @export
expand_power_exp <- .hint_expand("power_exp")
#' @rdname expand
#' @export
expand_power_base <- .hint_expand("power_base")
#' @rdname expand
#' @export
expand_trig <- .hint_expand("trig")
#' @rdname expand
#' @export
expand_complex <- .hint_expand("complex")
#' @rdname expand
#' @export
expand_func <- .hint_expand("func")
