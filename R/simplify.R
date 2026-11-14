# Simplification, after sympy/simplify.

.measure <- function(e) .cops(e)

.shortest <- function(cands) {
  cands <- Filter(Negate(is.null), cands)
  ops <- vapply(cands, .measure, 0L)
  len <- vapply(cands, function(z) nchar(z$key), 0L)
  cands[[order(ops, len)[1L]]]
}


# ---------------------------------------------------------------------------
# Powers and logarithms


# ---------------------------------------------------------------------------
# Radicals


# ---------------------------------------------------------------------------
# sqrtdenest: sympy/simplify/sqrtdenest.py

.sdn_is_sqrt <- function(e) e$type == "Pow" && e$args[[2L]]$type == "Rational" &&
  gmp::denominator(e$args[[2L]]$value) == 2 && abs(gmp::numerator(e$args[[2L]]$value)) == 1
.sdn_depth <- function(p) {
  if (is_const(p, "I")) return(1L)
  if (is_atom(p)) return(0L)
  if (p$type %in% c("Add", "Mul")) return(max(vapply(p$args, .sdn_depth, 0L)))
  if (.sdn_is_sqrt(p)) return(.sdn_depth(p$args[[1L]]) + 1L)
  0L
}
.sdn_is_algebraic <- function(p) {
  if (p$type %in% c("Integer", "Rational")) return(TRUE)
  if (is_atom(p)) return(FALSE)
  if (.sdn_is_sqrt(p) || (p$type == "Pow" && p$args[[2L]]$type == "Integer")) return(.sdn_is_algebraic(p$args[[1L]]))
  if (p$type %in% c("Add", "Mul")) return(all(vapply(p$args, .sdn_is_algebraic, NA)))
  FALSE
}
.sdn_subsets <- function(n) {
  if (n == 1L) return(list(1L))
  if (n == 2L) return(list(c(1L, 0L), c(0L, 1L), c(1L, 1L)))
  if (n == 3L) return(list(c(1L, 0L, 0L), c(0L, 1L, 0L), c(1L, 1L, 0L), c(0L, 0L, 1L), c(1L, 0L, 1L), c(0L, 1L, 1L), c(1L, 1L, 1L)))
  b <- .sdn_subsets(n - 1L)
  c(lapply(b, function(x) c(x, 0L)), list(c(rep(0L, n - 1L), 1L)), lapply(b, function(x) c(x, 1L)))
}
# Python truth of a numeric relational; undecided is FALSE
.sdn_true <- function(r) is_sym(r) && r$type == "BooleanAtom" && r$name == "true"
.sdn_sorted <- function(args) .key_sort(lapply(args, .sort_key), args)
.sdn_is_rat <- function(e) e$type %in% c("Integer", "Rational")
.sdn_stop <- function() stop(structure(class = c("SqrtdenestStopIteration", "error", "condition"),
                                       list(message = "sqrtdenest stop", call = NULL)))

#' @rdname radsimp
#' @param max_iter For `sqrtdenest()`: the largest number of denesting passes.
#' @export
sqrtdenest <- function(e, max_iter = 3L) {
  e <- expand_mul(as_sym(e))
  for (i in seq_len(max_iter)) {
    z <- .sdn0(e)
    if (identical(e$key, z$key)) return(e)
    e <- z
  }
  e
}

.sdn_match <- function(p) {
  p <- .mexpand(p)
  if (is_number_type(p)) return(list(p, S_Zero, S_Zero))
  if (p$type == "Add") {
    pargs <- .sdn_sorted(p$args)
    if (all(vapply(pargs, function(x) { s <- .pow(x, S_Two); .sdn_is_rat(s) && isTRUE(ask_fact(s, "positive")) }, NA))) {
      sp <- split_surds(p)
      return(list(sp[[3L]], sp[[2L]], sp[[1L]]))
    }
    depths <- vapply(pargs, .sdn_depth, 0L)
    # max over (depth, x, i) by default_sort_key: the deepest, then the
    # largest term by sort key
    cand <- which(depths == max(depths))
    i <- cand[[length(cand)]]
    if (length(cand) > 1L) {
      keys <- lapply(pargs[cand], .sort_key)
      i <- cand[[match(.key_sort(keys, pargs[cand])[[length(cand)]]$key, vapply(pargs[cand], `[[`, "", "key"))]]
    }
    depth <- depths[[i]]
    if (depth == 0L) return(list())
    r <- pargs[[i]]
    rest <- pargs[-i]; rdep <- depths[-i]
    b <- S_One
    if (r$type == "Mul") {
      ra <- r$args
      low <- vapply(ra, function(x) .sdn_depth(x) < depth, NA)
      b <- .mul(ra[low]); r <- .mul(ra[!low])
    }
    a1 <- list(); b1 <- list(b)
    for (j in seq_along(rest)) {
      x1 <- rest[[j]]
      if (rdep[[j]] < depth) a1[[length(a1) + 1L]] <- x1
      else if (identical(x1$key, r$key)) b1[[length(b1) + 1L]] <- S_One
      else if (x1$type == "Mul" && r$key %in% vapply(x1$args, `[[`, "", "key")) {
        xa <- x1$args
        xa <- xa[-match(r$key, vapply(xa, `[[`, "", "key"))]
        b1[[length(b1) + 1L]] <- .mul(xa)
      } else a1[[length(a1) + 1L]] <- x1
    }
    return(list(.add(a1), .add(b1), .pow(r, S_Two)))
  }
  cm <- as_coeff_mul_split(p)
  if (.sdn_is_sqrt(cm[[2L]])) list(S_Zero, cm[[1L]], .pow(cm[[2L]], S_Two)) else list()
}

.sdn0 <- function(e) {
  if (.sdn_is_sqrt(e)) {
    nd <- .as_numer_denom(e); n <- nd[[1L]]; d <- nd[[2L]]
    if (is_one_num(d)) {
      if (n$args[[1L]]$type == "Add") {
        args <- .sdn_sorted(n$args[[1L]]$args)
        if (length(args) > 2L && all(vapply(args, function(x) .pow(x, S_Two)$type == "Integer", NA))) {
          r <- tryCatch(.sdn_rec(n), SqrtdenestStopIteration = function(c) NULL)
          if (!is.null(r)) return(r)
        }
        e <- .sqrt(.mexpand(.add(lapply(args, .sdn0))))
      }
      return(.sdn1(e))
    }
    return(.div(.sdn0(n), .sdn0(d)))
  }
  if (e$type == "Add") {
    cm <- lapply(.sympy_args(e), as_coeff_mul_split)
    cs <- lapply(cm, `[[`, 1L); args <- lapply(cm, `[[`, 2L)
    if (all(vapply(cs, .sdn_is_rat, NA)) && all(vapply(args, .sdn_is_sqrt, NA))) return(.sdn_ratcomb(cs, args))
  }
  if (length(e$args)) return(sym_rebuild(e, lapply(e$args, .sdn0)))
  e
}

.sdn_rec <- function(e) {
  if (e$type != "Pow") return(sqrtdenest(e))
  base <- e$args[[1L]]
  if (.sdn_true(Lt(base, S_Zero))) return(.mul(list(S_I, .sdn_rec(.sqrt(.neg(base))))))
  sp <- split_surds(base)
  a <- .mul(list(sp[[2L]], .sqrt(sp[[1L]]))); b <- sp[[3L]]
  if (.sdn_true(Lt(a, b))) { t <- a; a <- b; b <- t }
  c2 <- .mexpand(.sub(.pow(a, S_Two), .pow(b, S_Two)))
  if (length(c2$args) > 2L) {
    sp1 <- split_surds(c2)
    a1 <- .mul(list(sp1[[2L]], .sqrt(sp1[[1L]]))); b1 <- sp1[[3L]]
    if (.sdn_true(Lt(a1, b1))) { t <- a1; a1 <- b1; b1 <- t }
    c2_1 <- .mexpand(.sub(.pow(a1, S_Two), .pow(b1, S_Two)))
    c_1 <- .sdn_rec(.sqrt(c2_1))
    d_1 <- .sdn_rec(.sqrt(.add(list(a1, c_1))))
    nd <- rad_rationalize(b1, d_1)
    c <- .mexpand(.add(list(.div(d_1, .sqrt(S_Two)), .div(nd[[1L]], .mul(list(nd[[2L]], .sqrt(S_Two)))))))
  } else c <- .sdn1(.sqrt(c2))
  if (.sdn_depth(c) > 1L) .sdn_stop()
  ac <- .add(list(a, c))
  if (length(ac$args) >= length(e$args) && count_ops(ac) >= count_ops(base)) .sdn_stop()
  d <- sqrtdenest(.sqrt(ac))
  if (.sdn_depth(d) > 1L) .sdn_stop()
  nd <- rad_rationalize(b, d)
  r <- .add(list(.div(d, .sqrt(S_Two)), .div(nd[[1L]], .mul(list(nd[[2L]], .sqrt(S_Two))))))
  .mexpand(radsimp(r))
}

.sdn1 <- function(e, denester = TRUE) {
  if (!.sdn_is_sqrt(e)) return(e)
  a <- e$args[[1L]]
  if (is_atom(a)) return(e)
  val <- .sdn_match(a)
  if (!length(val)) return(e)
  a <- val[[1L]]; b <- val[[2L]]; r <- val[[3L]]
  d2 <- .mexpand(.sub(.pow(a, S_Two), .mul(list(.pow(b, S_Two), r))))
  if (.sdn_is_rat(d2)) {
    if (d2$value > 0) {
      z <- .sdn_numeric(a, b, r, d2)
      if (!is.null(z)) return(z)
    } else {
      dr2 <- .mexpand(.neg(.mul(list(d2, r))))
      dr <- .sqrt(dr2)
      if (.sdn_is_rat(dr)) {
        z <- .sdn_numeric(.mexpand(.mul(list(b, r))), a, r, dr2)
        if (!is.null(z)) return(.div(z, .pow(r, .rat(1, 4))))
      }
    }
  } else {
    z <- .sdn_symbolic(a, b, r)
    if (!is.null(z)) return(z)
  }
  if (!denester || !.sdn_is_algebraic(e)) return(e)
  res <- .sdn_biquadratic(e, a, b, r, d2)
  if (!is.null(res)) return(res)
  av0 <- new.env()
  av0$v <- list(a, b, r, d2)
  z <- .sdn_denester(list(radsimp(.pow(e, S_Two))), av0, 0L, .sdn_depth(e))[[1L]]
  if (is.null(av0$v[[2L]])) return(e)
  if (!is.null(z)) {
    if (.sdn_depth(z) == .sdn_depth(e) && count_ops(z) > count_ops(e)) return(e)
    return(z)
  }
  e
}

.sdn_symbolic <- function(a, b, r) {
  rval <- .sdn_match(r)
  if (!length(rval)) return(NULL)
  ra <- rval[[1L]]; rb <- rval[[2L]]; rr <- rval[[3L]]
  if (is_zero_num(rb)) return(NULL)
  y <- Dummy("y", positive = TRUE)
  newa <- tryCatch(Poly(subs(a, .sqrt(rr), .div(.sub(.pow(y, S_Two), ra), rb)), y), error = function(err) NULL)
  if (is.null(newa) || degree(newa) != 2) return(NULL)
  cf <- lapply(all_coeffs(newa), as_sym)
  ca <- cf[[1L]]; cb <- .add(list(cf[[2L]], b)); cc <- cf[[3L]]
  if (!is_zero_num(simplify(.mexpand(.sub(.pow(cb, S_Two), .mul(list(.int(4), ca, cc))))))) return(NULL)
  z <- .sqrt(.mul(list(ca, .pow(.add(list(.sqrt(r), .div(cb, .mul(list(S_Two, ca))))), S_Two))))
  if (.is_number(z)) {
    cp <- as_content_primitive(z)
    z <- .mexpand(.mul(list(cp[[1L]], cp[[2L]])))
  }
  z
}

.sdn_numeric <- function(a, b, r, d2) {
  d <- .sqrt(d2)
  s <- .add(list(a, d))
  if (.sdn_depth(s) < .sdn_depth(r) + 1L || .sdn_is_rat(.pow(s, S_Two))) {
    s1 <- fcall("sign", list(s)); s2 <- fcall("sign", list(b))
    if (is_number_type(s1) && is_number_type(s2) && .int_val(s1) == -1 && .int_val(s2) == -1) { s1 <- S_One; s2 <- S_One }
    res <- .div(.mul(list(.add(list(.mul(list(s1, .sqrt(.add(list(a, d))))), .mul(list(s2, .sqrt(.sub(a, d)))))), .sqrt(S_Two))), S_Two)
    return(expand(res))
  }
  NULL
}

.sdn_biquadratic <- function(e, a, b, r, d2) {
  if (!.sdn_true(Gt(r, S_Zero)) || .sdn_true(Lt(d2, S_Zero)) || is_zero_num(b) || .sdn_depth(e$args[[1L]]) < 2L) return(NULL)
  for (x in list(a, b, r)) for (y in x$args) {
    y2 <- .pow(y, S_Two)
    if (y2$type != "Integer" || y2$value <= 0) return(NULL)
  }
  sqd <- .mexpand(sqrtdenest(.sqrt(radsimp(d2))))
  if (.sdn_depth(sqd) > 1L) return(NULL)
  for (x in list(.add(list(.div(a, S_Two), .div(sqd, S_Two))), .sub(.div(a, S_Two), .div(sqd, S_Two)))) {
    A <- sqrtdenest(.sqrt(x))
    if (.sdn_depth(A) > 1L) next
    nd <- rad_rationalize(b, .mexpand(.mul(list(S_Two, A))))
    z <- .add(list(A, .mul(list(.div(nd[[1L]], nd[[2L]]), .sqrt(r)))))
    if (.sdn_true(Lt(z, S_Zero))) z <- .neg(z)
    return(.mexpand(z))
  }
  NULL
}

# _denester; av0 is an environment holding the list av0 that SymPy mutates
.sdn_denester <- function(nested, av0, h, max_depth_level) {
  if (h > max_depth_level) return(list(NULL, NULL))
  if (is.null(av0$v[[2L]])) return(list(NULL, NULL))
  n <- length(nested)
  if (is.null(av0$v[[1L]]) && all(vapply(nested, is_number_type, NA))) {
    for (f in .sdn_subsets(n)) {
      p <- .mexpand(.mul(nested[f == 1L]))
      if (sum(f) > 1L && f[[length(f)]] == 1L) p <- .neg(p)
      sqp <- .sqrt(p)
      if (.sdn_is_rat(sqp)) return(list(sqp, f))
    }
    return(list(.sqrt(nested[[n]]), rep(0L, n)))
  }
  R <- NULL
  if (!is.null(av0$v[[1L]])) {
    values <- list(av0$v[1:2])
    R <- av0$v[[3L]]
    nested2 <- list(av0$v[[4L]], R)
    av0$v[1L] <- list(NULL)
  } else {
    values <- Filter(length, lapply(nested, .sdn_match))
    for (v in values) if (!is_zero_num(v[[3L]])) {
      if (!is.null(R)) {
        if (!identical(R$key, v[[3L]]$key)) { av0$v[2L] <- list(NULL); return(list(NULL, NULL)) }
      } else R <- v[[3L]]
    }
    if (is.null(R)) return(list(.sqrt(nested[[n]]), rep(0L, n)))
    nested2 <- c(lapply(values, function(v) .sub(.mexpand(.pow(v[[1L]], S_Two)), .mexpand(.mul(list(R, .pow(v[[2L]], S_Two)))))), list(R))
  }
  df <- .sdn_denester(nested2, av0, h + 1L, max_depth_level)
  d <- df[[1L]]; f <- df[[2L]]
  if (!length(f)) return(list(NULL, NULL))
  if (!any(f[seq_len(n)] == 1L)) {
    v <- values[[length(values)]]
    return(list(.sqrt(.add(list(v[[1L]], .mexpand(.mul(list(v[[2L]], d)))))), f))
  }
  p <- .mul(nested[which(f[seq_len(n)] == 1L)])
  v <- .sdn_match(p)
  if (1L %in% f && match(1L, f) < n && f[[n]] == 1L) { v[[1L]] <- .neg(v[[1L]]); v[[2L]] <- .neg(v[[2L]]) }
  if (f[[n + 1L]] == 0L) {
    vad <- .mexpand(.add(list(v[[1L]], d)))
    if (.sdn_true(Le(vad, S_Zero))) return(list(.sqrt(nested[[n]]), rep(0L, n)))
    if (!(.sdn_depth(vad) <= .sdn_depth(R) + 1L || is_number_type(.pow(vad, S_Two)))) { av0$v[2L] <- list(NULL); return(list(NULL, NULL)) }
    sqvad <- .sdn1(.sqrt(vad), denester = FALSE)
    if (!(.sdn_depth(sqvad) <= .sdn_depth(R) + 1L)) { av0$v[2L] <- list(NULL); return(list(NULL, NULL)) }
    sqvad1 <- radsimp(.inv(sqvad))
    res <- .mexpand(.add(list(.div(sqvad, .sqrt(S_Two)), .div(.mul(list(v[[2L]], .sqrt(R), sqvad1)), .sqrt(S_Two)))))
    return(list(res, f))
  }
  s2 <- .add(list(.mexpand(.mul(list(v[[2L]], R))), d))
  if (.sdn_true(Le(s2, S_Zero))) return(list(.sqrt(nested[[n]]), rep(0L, n)))
  FR <- .pow(.mexpand(R), .rat(1, 4)); s <- .sqrt(s2)
  list(.mexpand(.add(list(.div(s, .mul(list(.sqrt(S_Two), FR))), .div(.mul(list(v[[1L]], FR)), .mul(list(.sqrt(S_Two), s)))))), f)
}

.sdn_ratcomb <- function(cs, args) {
  idx <- NULL
  n <- length(args)
  if (n > 1L) for (i in seq_len(n - 1L)) {
    for (j in (i + 1L):n) {
      p <- .mexpand(.mul(list(args[[i]]$args[[1L]], args[[j]]$args[[1L]])))
      s <- sqrtdenest(.sqrt(p))
      if (!identical(s$key, .sqrt(p)$key)) { idx <- list(s, i, j); break }
    }
    if (!is.null(idx)) break
  }
  if (is.null(idx)) return(.add(Map(function(c, a) .mul(list(c, a)), cs, args)))
  s <- idx[[1L]]; i1 <- idx[[2L]]; i2 <- idx[[3L]]
  c2 <- cs[[i2]]; cs <- cs[-i2]; args <- args[-i2]
  cs[[i1]] <- .add(list(cs[[i1]], radsimp(.div(.mul(list(c2, s)), args[[i1]]$args[[1L]]))))
  .sdn_ratcomb(cs, args)
}

is_const_num_exp <- function(z, q) z$args[[2L]]$type == "Rational" && z$args[[2L]]$value == q

# ---------------------------------------------------------------------------
# Collection and rational simplification


#' Separate variables
#'
#' `separatevars(e)` factors `e` into a product of factors that each depend
#' on fewer symbols, when possible.
#'
#' @param e An expression.
#' @param syms Optional symbols; with `dict = TRUE` returns a named list.
#' @param dict Return a list of factors per symbol.
#' @return An expression or a list.
#' @export
separatevars <- function(e, syms = NULL, dict = FALSE) {
  e <- as_sym(e)
  f <- .separatevars(e)
  if (!dict) return(f)
  syms <- if (is.null(syms)) free_symbols(e) else .as_gens(list(syms))
  facs <- if (f$type == "Mul") f$args else list(f)
  out <- list(coeff = S_One)
  for (s in syms) out[[s$name]] <- S_One
  for (fc in facs) {
    dep <- Filter(function(s) .has1(fc, s), syms)
    if (length(dep) > 1L) return(NULL)
    if (!length(dep)) out$coeff <- .mul(list(out$coeff, fc))
    else out[[dep[[1L]]$name]] <- .mul(list(out[[dep[[1L]]$name]], fc))
  }
  out
}

# SymPy _separatevars: keep a Mul and separate its factors, separate the
# base of a power and expand it, else factor.
.separatevars <- function(e) {
  if (length(free_symbols(e)) < 2L) return(e)
  if (e$type == "Mul") return(.mul(lapply(e$args, .separatevars)))
  if (e$type == "Pow" && !is_const(e$args[[1L]], "E"))
    e <- expand_power_base(.pow(.separatevars(e$args[[1L]]), e$args[[2L]]))
  factor_sym(e)
}

.nsimplify_num <- function(v, constants, tol, rational) {
  if (v == 0) return(S_Zero)
  try_rat <- function(w, maxden = 10000) {
    r <- .best_rational(w, maxden)
    if (!is.null(r) && abs(as.double(r) - w) <= tol * max(1, abs(w))) return(make_rational(r))
    NULL
  }
  r <- try_rat(v)
  if (!is.null(r) || rational) return(if (is.null(r)) make_rational(gmp::as.bigq(v)) else r)
  for (cst in c(constants, list(S_Pi, S_E))) {
    cv <- .neval(as_sym(cst), list())
    r <- try_rat(v / cv, 1000)
    if (!is.null(r)) return(.mul(list(r, as_sym(cst))))
  }
  r2 <- try_rat(v^2, 1000)
  if (!is.null(r2)) return(.mul(list(if (v < 0) S_NegativeOne else S_One, .sqrt(r2))))
  for (k in 2:12) {
    r <- try_rat(v - sqrt(k), 1000)
    if (!is.null(r)) return(.add(list(r, .sqrt(.int(k)))))
    r <- try_rat(v + sqrt(k), 1000)
    if (!is.null(r)) return(.sub(r, .sqrt(.int(k))))
  }
  Float(v)
}

.best_rational <- function(w, maxden) {
  # continued fraction convergents
  a <- floor(w)
  h0 <- 1; h1 <- a
  k0 <- 0; k1 <- 1
  frac <- w - a
  best <- gmp::as.bigq(a)
  for (i in 1:40) {
    if (abs(frac) < 1e-15) break
    w <- 1 / frac
    a <- floor(w)
    frac <- w - a
    h2 <- a * h1 + h0
    k2 <- a * k1 + k0
    if (k2 > maxden) break
    best <- gmp::as.bigq(h2, k2)
    h0 <- h1; h1 <- h2
    k0 <- k1; k1 <- k2
  }
  best
}

#' Simplify signs
#'
#' `signsimp()` takes the sign out of sums inside products, as SymPy's
#' `signsimp()` (`-x*(-x + y) + y` becomes `x*(x - y) + y`), and `posify()`
#' replaces symbols with positive ones.
#'
#' @param e An expression.
#' @return An expression; `posify()` returns `list(expr, reps)`.
#' @export
signsimp <- function(e) {
  e <- as_sym(e)
  if (!length(e$args)) return(e)
  # sub_pre: a sum from which -1 can be extracted becomes -1*(-sum), and
  # re-evaluating the products folds that -1 into their coefficients
  adds <- Filter(could_extract_minus_sign, atoms(e, "Add"))
  reps <- list()
  for (a in adds) {
    na <- .neg(a)
    if (na$type == "Mul") next
    reps[[length(reps) + 1L]] <- list(a, new_sym("Mul", list(S_NegativeOne, na)))
  }
  if (length(reps)) e <- xreplace(e, reps)
  e <- .map_bottom_up(e, function(z) if (z$type == "Mul") .mul(z$args) else z)
  if (e$type == "Add") return(.add(lapply(e$args, signsimp)))
  e
}

#' @rdname signsimp
#' @export
posify <- function(e) {
  e <- as_sym(e)
  fs <- Filter(function(s) is.na(ask_fact(s, "positive")), free_symbols(e))
  reps <- lapply(fs, function(s) list(s, Dummy(s$name, positive = TRUE)))
  list(expr = xreplace(e, reps), reps = lapply(reps, rev))
}

