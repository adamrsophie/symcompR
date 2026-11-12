# The Fu et al. trigonometric simplification, after sympy/simplify/fu.py:
# the transformation rules TR0 to TR22, TR2i, TR10i, TR12i, TRmorrie,
# TRpower, and the fu() driver that combines them greedily.

# ---------------------------------------------------------------------------
# Helpers standing in for SymPy primitives

.fu_trig <- c("sin", "cos", "tan", "cot", "sec", "csc")
.fu_is <- function(z, nm) is_sym(z) && z$type == "Function" && z$name %in% nm && length(z$args) == 1L
.fu_arg <- function(z) z$args[[1L]]
.fu_f <- function(nm, a) fcall(nm, list(a))
.fu_raw <- function(nm, a) fnode(nm, list(a))
# SymPy's bottom_up: f applies to every node that has arguments
.fu_bu <- function(e, f) .map_bottom_up(as_sym(e), function(z) if (length(z$args)) f(z) else z)
.fu_ordered <- function(xs) .key_sort(lapply(xs, .sort_key), xs)
.fu_margs <- function(z) if (z$type == "Mul") z$args else list(z)
.fu_nd <- function(z) .as_numer_denom(z)

# Expr.as_numer_denom: a sum is put over a common denominator (after its
# rational content is taken out), a product combines its factors' parts,
# and a power splits its base when the sign of the denominator allows.
.as_numer_denom <- function(e) {
  e <- as_sym(e)
  switch(e$type,
    Rational = list(.zs(gmp::numerator(e$value)), .zs(gmp::denominator(e$value))),
    Mul = {
      if (!.is_comm(e)) return(list(e, S_One))
      nds <- lapply(e$args, .as_numer_denom)
      list(.mul(lapply(nds, `[[`, 1L)), .mul(lapply(nds, `[[`, 2L)))
    },
    Pow = {
      if (!.is_comm(e)) return(list(e, S_One))
      b <- e$args[[1L]]; x <- e$args[[2L]]
      nd <- .as_numer_denom(b); n <- nd[[1L]]; d <- nd[[2L]]
      neg <- isTRUE(ask_fact(x, "negative"))
      if (x$type == "Mul" && !neg && !isTRUE(ask_fact(x, "positive"))) neg <- could_extract_minus_sign(x)
      ie <- isTRUE(ask_fact(x, "integer"))
      if (!(isTRUE(ask_fact(d, "real")) || ie)) { n <- b; d <- S_One }
      dnp <- ask_fact(d, "nonpositive")
      if (isTRUE(dnp)) { n <- .neg(n); d <- .neg(d) }
      else if (is.na(dnp) && !ie) { n <- b; d <- S_One }
      if (neg) { x <- .neg(x); t <- n; n <- d; d <- t }
      if (is_inf_const(x)) {
        if (is_one_num(n) && !is_one_num(d)) return(list(n, .pow(d, x)))
        if (!is_one_num(n) && is_one_num(d)) return(list(.pow(n, x), d))
      }
      list(.pow(n, x), .pow(d, x))
    },
    # ExpBase.as_numer_denom
    Function = {
      if (e$name %in% c("exp", "exp_polar") && .is_comm(e)) {
        x <- e$args[[1L]]
        ng <- isTRUE(ask_fact(x, "negative"))
        if (!ng && !isTRUE(ask_fact(.neg(x), "negative"))) ng <- could_extract_minus_sign(x)
        if (ng) return(list(S_One, fcall(e$name, list(.neg(x)))))
      }
      list(e, S_One)
    },
    Add = {
      # Add.primitive: the rational content of the terms
      cs <- lapply(e$args, function(t) as_coeff_mul_split(t)[[1L]])
      rat <- all(vapply(cs, function(c) c$type %in% c("Integer", "Rational"), TRUE))
      con <- S_One
      if (rat) {
        qs <- lapply(cs, function(c) gmp::as.bigq(c$value))
        g <- Reduce(gmp::gcd.bigz, lapply(qs, function(q) abs(gmp::numerator(q))))
        l <- Reduce(function(a, b) a * b / gmp::gcd.bigz(a, b), lapply(qs, gmp::denominator))
        con <- make_rational(gmp::as.bigq(g, l))
      }
      terms <- if (is_one_num(con)) e$args else lapply(e$args, function(t) .div(t, con))
      ncon <- .zs(gmp::numerator(gmp::as.bigq(con$value))); dcon <- .zs(gmp::denominator(gmp::as.bigq(con$value)))
      groups <- list(); dens <- list()
      for (t in terms) {
        nd <- .as_numer_denom(t)
        k <- nd[[2L]]$key
        groups[[k]] <- c(groups[[k]], list(nd[[1L]])); dens[[k]] <- nd[[2L]]
      }
      if (length(groups) == 1L)
        return(list(.mul(list(ncon, .add(groups[[1L]]))), .mul(list(dcon, dens[[1L]]))))
      nums <- lapply(groups, .add)
      ks <- names(groups)
      n <- .add(lapply(seq_along(ks), function(i) .mul(c(unname(dens[ks[-i]]), list(nums[[i]])))))
      list(.mul(list(ncon, n)), .mul(list(dcon, .mul(unname(dens)))))
    },
    list(e, S_One))
}
.fu_has <- function(e, nms) any(vapply(atoms(e, "Function"), function(f) f$name %in% nms, TRUE))
.fu_eq <- function(a, b) identical(a$key, b$key)
.fu_div <- function(a, b) .mul(list(a, .pow(b, S_NegativeOne)))
.fu_is_int <- function(z) isTRUE(ask_fact(z, "integer"))
.fu_is_pos <- function(z) isTRUE(ask_fact(z, "positive"))

# Mul.as_powers_dict: base key -> list(base, exponent)
.fu_powers <- function(z) {
  out <- list()
  for (t in .fu_margs(z)) {
    be <- .base_exp(t)
    k <- be[[1L]]$key
    out[[k]] <- if (is.null(out[[k]])) be else list(be[[1L]], .add(list(out[[k]][[2L]], be[[2L]])))
  }
  out
}
.fu_from_powers <- function(d) .mul(lapply(Filter(Negate(is.null), d), function(p) .pow(p[[1L]], p[[2L]])))

# exprtools.Factors: normal() and gcd() of two products
.fac_expr <- function(m) .mul(lapply(m, function(p) .pow(p[[1L]], p[[2L]])))
.fac_normal <- function(a, b) {
  fa <- .factors_map(a); fb <- .factors_map(b)
  sa <- fa; sb <- fb
  for (k in names(fa)) {
    if (is.null(fb[[k]])) next
    delta <- .sub(fa[[k]][[2L]], fb[[k]][[2L]])
    if (isTRUE(ask_fact(delta, "nonnegative"))) {
      sb[[k]] <- NULL
      if (is_zero_num(delta)) sa[[k]] <- NULL else sa[[k]][[2L]] <- delta
    } else if (isTRUE(ask_fact(delta, "nonpositive"))) {
      sa[[k]] <- NULL
      sb[[k]][[2L]] <- .neg(delta)
    }
  }
  list(sa, sb)
}
.fac_gcd <- function(a, b) {
  fa <- .factors_map(a); fb <- .factors_map(b)
  g <- list()
  for (k in names(fa)) {
    if (is.null(fb[[k]])) next
    lt <- ask_fact(.sub(fa[[k]][[2L]], fb[[k]][[2L]]), "negative")
    if (isTRUE(lt)) g[[k]] <- fa[[k]] else if (isFALSE(lt)) g[[k]] <- fb[[k]]
  }
  .fac_expr(g)
}

# The measure fu() minimizes: the number of trigonometric functions, then
# the operation count.
.fu_L <- function(e) {
  n <- 0L
  walk <- function(z) { if (.fu_is(z, .fu_trig)) n <<- n + 1L; for (a in z$args) walk(a) }
  walk(e)
  n
}
.fu_measure <- function(e) c(.fu_L(e), count_ops(e))
.fu_less <- function(a, b) { for (i in seq_along(a)) if (a[i] != b[i]) return(a[i] < b[i]); FALSE }
.fu_min <- function(cands, measure) {
  best <- cands[[1L]]; mb <- measure(best)
  for (c in cands[-1L]) { mc <- measure(c); if (.fu_less(mc, mb)) { best <- c; mb <- mc } }
  best
}
# strategies.tree.greedy: a "chain" applies its steps in order, a "choice"
# applies each branch and keeps the best by the measure.
.fu_chain <- function(...) structure(list(...), class = "fu_chain")
.fu_choice <- function(...) structure(list(...), class = "fu_choice")
.fu_apply <- function(tree, e, measure) {
  if (is.function(tree)) return(tree(e))
  if (inherits(tree, "fu_chain")) { for (s in tree) e <- .fu_apply(s, e, measure); return(e) }
  .fu_min(lapply(tree, .fu_apply, e = e, measure = measure), measure)
}

# ---------------------------------------------------------------------------
# The transformation rules

#' Fu trigonometric simplification
#'
#' `fu(rv)` simplifies trigonometric expressions with the algorithm of Fu,
#' Zhong and Zeng, as SymPy's `fu()`: it applies the transformation rules
#' below in greedy combinations and keeps the result with the fewest
#' trigonometric functions, then the fewest operations. Each rule is also
#' exported: `TR0` (normal, factor, expand), `TR1` (sec, csc to cos, sin),
#' `TR2` (tan, cot to sin, cos), `TR2i` (the reverse), `TR3` (induced
#' formulas), `TR4` (special angles), `TR5` and `TR6` (sin^2 and cos^2 to
#' the other), `TR7` (lower the degree of cos^2), `TR8` (products to sums),
#' `TR9` (sums to products), `TR10` and `TR10i` (expand and contract sums
#' of angles), `TR11` and `TR11` with a base (double angles), `TR12` and
#' `TR12i` (tangent of sums), `TR13` (products of tangents), `TRmorrie`
#' (Morrie's law), `TR14` (factored powers of sin and cos), `TR15`, `TR16`
#' and `TR22` (powers in terms of cot, tan, sec and csc), `TR111` (negative
#' powers to the reciprocal functions) and `TRpower` (powers to multiple
#' angles).
#'
#' @param rv An expression.
#' @param measure The function to minimize, returning a numeric vector
#'   compared lexicographically.
#' @param half For `TR2i()`, also recognize half angle forms.
#' @param first Internal; leave the default.
#' @param base For `TR11()`, the angle to write in terms of.
#' @param max,pow For the power rules: the largest exponent, and whether
#'   perfect powers are handled.
#' @return An expression.
#' @examples
#' x <- Symbol("x"); y <- Symbol("y")
#' fu(sin(50)^2 + cos(50)^2 + sin(pi/6))
#' fu(sqrt(6)*cos(x) + sqrt(2)*sin(x))
#' TR10i(cos(1)*cos(3) + sin(1)*sin(3))
#' TR8(cos(x)*sin(y))
#' @name fu
NULL

#' @rdname fu
#' @export
TR0 <- function(rv) {
  nd <- .fu_nd(as_sym(rv))
  expand(sym_factor(.fu_div(nd[[1L]], nd[[2L]])))
}

#' @rdname fu
#' @export
TR1 <- function(rv) .fu_bu(rv, function(z) {
  if (.fu_is(z, "sec")) return(.pow(.fu_f("cos", .fu_arg(z)), S_NegativeOne))
  if (.fu_is(z, "csc")) return(.pow(.fu_f("sin", .fu_arg(z)), S_NegativeOne))
  z
})

#' @rdname fu
#' @export
TR2 <- function(rv) .fu_bu(rv, function(z) {
  if (.fu_is(z, "tan")) return(.fu_div(.fu_f("sin", .fu_arg(z)), .fu_f("cos", .fu_arg(z))))
  if (.fu_is(z, "cot")) return(.fu_div(.fu_f("cos", .fu_arg(z)), .fu_f("sin", .fu_arg(z))))
  z
})

#' @rdname fu
#' @export
TR2i <- function(rv, half = FALSE) .fu_bu(rv, function(rv) {
  if (rv$type != "Mul") return(rv)
  nd <- .fu_nd(rv)
  if (!length(nd[[1L]]$args) || !length(nd[[2L]]$args)) return(rv)
  ok <- function(k, e) {
    (.fu_is_int(e) || .fu_is_pos(k)) &&
      (.fu_is(k, c("sin", "cos")) ||
         (half && k$type == "Add" && length(k$args) >= 2L &&
            any(vapply(k$args, function(a) any(vapply(.fu_margs(a), function(ai) .fu_is(ai, "cos") ||
              (ai$type == "Pow" && .fu_is(ai$args[[1L]], "cos")), TRUE)), TRUE))))
  }
  split <- function(z) {
    d <- .fu_powers(z)
    done <- list()
    for (k in names(d)) if (!ok(d[[k]][[1L]], d[[k]][[2L]])) { done[[length(done) + 1L]] <- d[[k]]; d[[k]] <- NULL }
    list(d, done)
  }
  sn <- split(nd[[1L]]); n <- sn[[1L]]; ndone <- sn[[2L]]
  if (!length(n)) return(rv)
  sd <- split(nd[[2L]]); d <- sd[[1L]]; ddone <- sd[[2L]]
  if (!length(d)) return(rv)
  factorize <- function(d, done) {
    newk <- list()
    for (k in names(d)) {
      b <- d[[k]][[1L]]
      if (b$type == "Add" && length(b$args) > 1L) {
        knew <- if (half) sym_factor(b) else factor_terms(b)
        if (!.fu_eq(knew, b)) newk[[k]] <- list(b, knew, d[[k]][[2L]])
      }
    }
    if (length(newk)) {
      for (k in names(newk)) d[[k]] <- NULL
      pd <- .fu_powers(.mul(lapply(newk, function(t) .pow(t[[2L]], t[[3L]]))))
      for (k in names(pd)) {
        v <- if (is.null(d[[k]])) pd[[k]][[2L]] else .add(list(d[[k]][[2L]], pd[[k]][[2L]]))
        if (ok(pd[[k]][[1L]], v)) d[[k]] <- list(pd[[k]][[1L]], v)
        else { d[[k]] <- NULL; done[[length(done) + 1L]] <- list(pd[[k]][[1L]], v) }
      }
    }
    list(d, done)
  }
  r <- factorize(n, ndone); n <- r[[1L]]; ndone <- r[[2L]]
  r <- factorize(d, ddone); d <- r[[1L]]; ddone <- r[[2L]]
  t <- list()
  for (k in names(n)) {
    if (is.null(n[[k]])) next
    b <- n[[k]][[1L]]; e <- n[[k]][[2L]]
    if (.fu_is(b, "sin")) {
      a <- .fu_raw("cos", .fu_arg(b))
      if (!is.null(d[[a$key]]) && .fu_eq(d[[a$key]][[2L]], e)) {
        t[[length(t) + 1L]] <- .pow(.fu_f("tan", .fu_arg(b)), e); n[[k]] <- NULL; d[[a$key]] <- NULL
      } else if (half) {
        a1 <- .add(list(S_One, a))
        if (!is.null(d[[a1$key]]) && .fu_eq(d[[a1$key]][[2L]], e)) {
          t[[length(t) + 1L]] <- .pow(.fu_f("tan", .div(.fu_arg(b), S_Two)), e); n[[k]] <- NULL; d[[a1$key]] <- NULL
        }
      }
    } else if (.fu_is(b, "cos")) {
      a <- .fu_raw("sin", .fu_arg(b))
      if (!is.null(d[[a$key]]) && .fu_eq(d[[a$key]][[2L]], e)) {
        t[[length(t) + 1L]] <- .pow(.fu_f("tan", .fu_arg(b)), .neg(e)); n[[k]] <- NULL; d[[a$key]] <- NULL
      }
    } else if (half && b$type == "Add" && is_one_num(b$args[[1L]]) && length(b$args) == 2L && .fu_is(b$args[[2L]], "cos")) {
      a <- .fu_raw("sin", .fu_arg(b$args[[2L]]))
      if (!is.null(d[[a$key]]) && .fu_eq(d[[a$key]][[2L]], e) && (.fu_is_int(d[[a$key]][[2L]]) || .fu_is_pos(a))) {
        t[[length(t) + 1L]] <- .pow(.fu_f("tan", .div(.fu_arg(a), S_Two)), .neg(e)); n[[k]] <- NULL; d[[a$key]] <- NULL
      }
    }
  }
  if (length(t)) {
    rv <- .fu_div(.mul(c(t, list(.fu_from_powers(n)))), .fu_from_powers(d))
    rv <- .mul(list(rv, .fu_div(.fu_from_powers(ndone), .fu_from_powers(ddone))))
  }
  rv
})

#' @rdname fu
#' @export
TR3 <- function(rv) {
  # numeric products inside trigonometric arguments are rebuilt first
  rv <- .fu_bu(rv, function(z) z)
  .fu_bu(rv, function(z) {
    if (!.fu_is(z, .fu_trig)) return(z)
    z <- .fu_f(z$name, signsimp(.fu_arg(z)))
    if (!.fu_is(z, .fu_trig)) return(z)
    a <- .fu_arg(z)
    if (isTRUE(ask_fact(.sub(a, .div(S_Pi, .int(4))), "positive")) && isTRUE(ask_fact(.sub(.div(S_Pi, S_Two), a), "positive"))) {
      fmap <- c(cos = "sin", sin = "cos", tan = "cot", cot = "tan", sec = "csc", csc = "sec")
      z <- .fu_f(fmap[[z$name]], .sub(.div(S_Pi, S_Two), a))
    }
    z
  })
}

#' @rdname fu
#' @export
TR4 <- function(rv) .fu_bu(rv, function(z) {
  if (!.fu_is(z, .fu_trig)) return(z)
  r <- .div(.fu_arg(z), S_Pi)
  if (r$type %in% c("Integer", "Rational") && as.integer(as.character(gmp::denominator(gmp::as.bigq(r$value)))) %in% c(1L, 2L, 3L, 4L, 6L))
    return(.fu_f(z$name, .fu_arg(z)))
  z
})

.TR56 <- function(rv, f, g, h, max, pow) .fu_bu(rv, function(z) {
  if (!(z$type == "Pow" && .fu_is(z$args[[1L]], f))) return(z)
  e <- z$args[[2L]]
  if (!isTRUE(ask_fact(e, "real"))) return(z)
  if (.is_true(Lt(e, S_Zero))) return(z)
  if (.is_true(Gt(e, as_sym(max)))) return(z)
  if (is_one_num(e)) return(z)
  x <- .fu_arg(z$args[[1L]])
  hg <- function() h(.pow(.fu_f(g, x), S_Two))
  if (e$type == "Integer" && .int_val(e) == 2) return(hg())
  if (e$type != "Integer") return(z)
  ev <- .Z(e)
  if (ev %% 2 == 1) return(.mul(list(.fu_f(f, x), .pow(hg(), .zs(ev %/% 2)))))
  if (ev == 4) k <- gmp::as.bigz(2)
  else if (!pow) { if (ev %% 2 != 0) return(z); k <- ev %/% 2 }
  else { if (isFALSE(perfect_power(ev))) return(z); k <- ev %/% 2 }
  .pow(hg(), .zs(k))
})

#' @rdname fu
#' @export
TR5 <- function(rv, max = 4, pow = FALSE) .TR56(rv, "sin", "cos", function(x) .sub(S_One, x), max, pow)
#' @rdname fu
#' @export
TR6 <- function(rv, max = 4, pow = FALSE) .TR56(rv, "cos", "sin", function(x) .sub(S_One, x), max, pow)

#' @rdname fu
#' @export
TR7 <- function(rv) .fu_bu(rv, function(z) {
  if (!(z$type == "Pow" && .fu_is(z$args[[1L]], "cos") && identical(z$args[[2L]]$key, S_Two$key))) return(z)
  .div(.add(list(S_One, .fu_f("cos", .mul(list(S_Two, .fu_arg(z$args[[1L]])))))), S_Two)
})

#' @rdname fu
#' @export
TR8 <- function(rv, first = TRUE) .fu_bu(rv, function(z) {
  if (!(z$type == "Mul" || (z$type == "Pow" && .fu_is(z$args[[1L]], c("cos", "sin")) &&
                             (.fu_is_int(z$args[[2L]]) || .fu_is_pos(z$args[[1L]]))))) return(z)
  if (first) {
    nd <- lapply(.fu_nd(z), expand_mul)
    nn <- TR8(nd[[1L]], first = FALSE); nd2 <- TR8(nd[[2L]], first = FALSE)
    if (!.fu_eq(nn, nd[[1L]]) || !.fu_eq(nd2, nd[[2L]])) {
      z <- gcd_terms(.fu_div(nn, nd2))
      if (z$type == "Mul" && length(z$args) == 2L && z$args[[1L]]$type %in% c("Integer", "Rational") && z$args[[2L]]$type == "Add")
        z <- .mul(list(z$args[[1L]], z$args[[2L]]))
    }
    return(z)
  }
  cs <- list(); sn <- list(); other <- list()
  for (a in .fu_margs(z)) {
    if (.fu_is(a, "cos")) cs[[length(cs) + 1L]] <- .fu_arg(a)
    else if (.fu_is(a, "sin")) sn[[length(sn) + 1L]] <- .fu_arg(a)
    else if (a$type == "Pow" && a$args[[2L]]$type == "Integer" && .int_val(a$args[[2L]]) > 0 && .fu_is(a$args[[1L]], c("cos", "sin"))) {
      reps <- rep(list(.fu_arg(a$args[[1L]])), .int_val(a$args[[2L]]))
      if (a$args[[1L]]$name == "cos") cs <- c(cs, reps) else sn <- c(sn, reps)
    } else other[[length(other) + 1L]] <- a
  }
  if (!((length(cs) && length(sn)) || length(cs) > 1L || length(sn) > 1L)) return(z)
  pop <- function(v) { x <- v[[length(v)]]; list(x, v[-length(v)]) }
  for (i in seq_len(min(length(cs), length(sn)))) {
    p1 <- pop(sn); a1 <- p1[[1L]]; sn <- p1[[2L]]
    p2 <- pop(cs); a2 <- p2[[1L]]; cs <- p2[[2L]]
    other[[length(other) + 1L]] <- .div(.add(list(.fu_f("sin", .add(list(a1, a2))), .fu_f("sin", .sub(a1, a2)))), S_Two)
  }
  while (length(cs) > 1L) {
    p1 <- pop(cs); a1 <- p1[[1L]]; cs <- p1[[2L]]; p2 <- pop(cs); a2 <- p2[[1L]]; cs <- p2[[2L]]
    other[[length(other) + 1L]] <- .div(.add(list(.fu_f("cos", .add(list(a1, a2))), .fu_f("cos", .sub(a1, a2)))), S_Two)
  }
  if (length(cs)) other[[length(other) + 1L]] <- .fu_f("cos", cs[[1L]])
  while (length(sn) > 1L) {
    p1 <- pop(sn); a1 <- p1[[1L]]; sn <- p1[[2L]]; p2 <- pop(sn); a2 <- p2[[1L]]; sn <- p2[[2L]]
    other[[length(other) + 1L]] <- .div(.add(list(.neg(.fu_f("cos", .add(list(a1, a2)))), .fu_f("cos", .sub(a1, a2)))), S_Two)
  }
  if (length(sn)) other[[length(other) + 1L]] <- .fu_f("sin", sn[[1L]])
  TR8(expand_mul(.mul(other)))
})

# The pairwise search over the terms of a sum shared by TR9 and TR10i.
.fu_pair_do <- function(rv, two_fn) {
  do <- function(rv) {
    if (rv$type != "Add") return(rv)
    args <- .fu_ordered(rv$args)
    if (length(args) != 2L) {
      hit <- FALSE
      for (i in seq_along(args)) {
        ai <- args[[i]]
        if (is.null(ai)) next
        for (j in seq_along(args)) {
          if (j <= i) next
          aj <- args[[j]]
          if (is.null(aj)) next
          was <- .add(list(ai, aj))
          new <- do(was)
          if (!.fu_eq(new, was)) { args[i] <- list(new); args[j] <- list(NULL); hit <- TRUE; break }
        }
      }
      if (hit) {
        rv <- .add(Filter(Negate(is.null), args))
        if (rv$type == "Add") rv <- do(rv)
      }
      return(rv)
    }
    r <- two_fn(args[[1L]], args[[2L]])
    if (is.null(r)) rv else r
  }
  do
}

#' @rdname fu
#' @export
TR9 <- function(rv) .fu_bu(rv, function(rv) {
  if (rv$type != "Add") return(rv)
  do <- .fu_pair_do(rv, function(u, v) {
    sp <- .trig_split(u, v)
    if (is.null(sp)) return(NULL)
    gcd <- sp$gcd; n1 <- sp$n1; n2 <- sp$n2; a <- sp$a; b <- sp$b
    if (sp$flag) {
      if (n1 == n2) return(.mul(list(gcd, .int(n1), S_Two, .fu_f("cos", .div(.add(list(a, b)), S_Two)), .fu_f("cos", .div(.sub(a, b), S_Two)))))
      if (n1 < 0) { t <- a; a <- b; b <- t }
      return(.mul(list(.int(-2), gcd, .fu_f("sin", .div(.add(list(a, b)), S_Two)), .fu_f("sin", .div(.sub(a, b), S_Two)))))
    }
    if (n1 == n2) return(.mul(list(gcd, .int(n1), S_Two, .fu_f("sin", .div(.add(list(a, b)), S_Two)), .fu_f("cos", .div(.sub(a, b), S_Two)))))
    if (n1 < 0) { t <- a; a <- b; b <- t }
    .mul(list(S_Two, gcd, .fu_f("cos", .div(.add(list(a, b)), S_Two)), .fu_f("sin", .div(.sub(a, b), S_Two))))
  })
  .process_common_addends(rv, do)
})

#' @rdname fu
#' @export
TR10 <- function(rv, first = TRUE) .fu_bu(rv, function(z) {
  if (!.fu_is(z, c("cos", "sin"))) return(z)
  f <- z$name; arg <- .fu_arg(z)
  if (arg$type != "Add") return(z)
  args <- if (first) .fu_ordered(arg$args) else arg$args
  a <- args[[length(args)]]
  b <- .add(args[-length(args)])
  if (b$type == "Add") {
    if (f == "sin") return(.add(list(.mul(list(.fu_f("sin", a), TR10(.fu_f("cos", b), first = FALSE))),
                                    .mul(list(.fu_f("cos", a), TR10(.fu_f("sin", b), first = FALSE))))))
    return(.sub(.mul(list(.fu_f("cos", a), TR10(.fu_f("cos", b), first = FALSE))),
                .mul(list(.fu_f("sin", a), TR10(.fu_f("sin", b), first = FALSE)))))
  }
  if (f == "sin") return(.add(list(.mul(list(.fu_f("sin", a), .fu_f("cos", b))), .mul(list(.fu_f("cos", a), .fu_f("sin", b))))))
  .sub(.mul(list(.fu_f("cos", a), .fu_f("cos", b))), .mul(list(.fu_f("sin", a), .fu_f("sin", b))))
})

#' @rdname fu
#' @export
TR10i <- function(rv) .fu_bu(rv, function(rv) {
  if (rv$type != "Add") return(rv)
  do <- .fu_pair_do(rv, function(u, v) {
    sp <- .trig_split(u, v, two = TRUE)
    if (is.null(sp)) return(NULL)
    gcd <- .mul(list(.int(sp$n1), sp$gcd)); a <- sp$a; b <- sp$b
    if (sp$flag) return(.mul(list(gcd, .fu_f("cos", if (sp$n1 == sp$n2) .sub(a, b) else .add(list(a, b))))))
    .mul(list(gcd, .fu_f("sin", if (sp$n1 == sp$n2) .add(list(a, b)) else .sub(b, a))))
  })
  rv <- .process_common_addends(rv, do, function(z) paste(vapply(.fu_ordered(free_symbols(z)), .k, ""), collapse = ","))
  r3 <- .sqrt(.int(3)); ir3 <- .pow(r3, S_NegativeOne)
  while (rv$type == "Add") {
    byrad <- list(); keys <- list()
    for (a in rv$args) {
      kk <- S_One
      if (a$type == "Mul") for (ai in a$args)
        if (ai$type == "Pow" && identical(ai$args[[2L]]$key, S_Half$key) && ai$args[[1L]]$type == "Integer") { kk <- ai; break }
      byrad[[kk$key]] <- c(byrad[[kk$key]], list(a)); keys[[kk$key]] <- kk
    }
    args <- list()
    for (ak in names(byrad)) {
      for (bb in list(.mul(list(r3, keys[[ak]])), ir3)) {
        bk <- bb$key
        if (is.null(byrad[[bk]])) next
        for (i in seq_along(byrad[[ak]])) {
          if (is.null(byrad[[ak]][[i]])) next
          for (j in seq_along(byrad[[bk]])) {
            if (is.null(byrad[[bk]][[j]])) next
            was <- .add(list(byrad[[ak]][[i]], byrad[[bk]][[j]]))
            new <- do(was)
            if (!.fu_eq(new, was)) {
              args[[length(args) + 1L]] <- new
              byrad[[ak]][i] <- list(NULL); byrad[[bk]][j] <- list(NULL)
              break
            }
          }
        }
      }
    }
    if (length(args)) rv <- .add(c(args, lapply(byrad, function(v) .add(Filter(Negate(is.null), v)))))
    else { rv <- do(rv); break }
  }
  rv
})

#' @rdname fu
#' @export
TR11 <- function(rv, base = NULL) .fu_bu(rv, function(z) {
  if (!.fu_is(z, c("cos", "sin"))) return(z)
  if (!is.null(base)) {
    base <- as_sym(base)
    t <- .fu_f(z$name, .mul(list(base, S_Two)))
    co <- S_One
    if (t$type == "Mul") { sp <- as_coeff_mul_split(t); co <- sp[[1L]]; t <- sp[[2L]] }
    if (!.fu_is(t, c("cos", "sin"))) return(z)
    if (.fu_eq(.fu_arg(z), .fu_arg(t))) {
      c <- .fu_f("cos", base); s <- .fu_f("sin", base)
      if (z$name == "cos") return(.fu_div(.sub(.pow(c, S_Two), .pow(s, S_Two)), co))
      return(.fu_div(.mul(list(S_Two, c, s)), co))
    }
    return(z)
  }
  if (is_number_type(.fu_arg(z))) return(z)
  sp <- as_coeff_mul_split(.fu_arg(z))
  c <- sp[[1L]]; m <- sp[[2L]]
  if (!(c$type %in% c("Integer", "Rational"))) return(z)
  p <- gmp::numerator(gmp::as.bigq(c$value)); q <- gmp::denominator(gmp::as.bigq(c$value))
  if (p %% 2 == 0) {
    arg <- .div(.mul(list(.zs(p %/% 2), m)), .zs(q))
    cc <- TR11(.fu_f("cos", arg)); ss <- TR11(.fu_f("sin", arg))
    z <- if (z$name == "sin") .mul(list(S_Two, ss, cc)) else .sub(.pow(cc, S_Two), .pow(ss, S_Two))
  }
  z
})

#' @rdname fu
#' @export
TR12 <- function(rv, first = TRUE) .fu_bu(rv, function(z) {
  if (!.fu_is(z, "tan")) return(z)
  arg <- .fu_arg(z)
  if (arg$type != "Add") return(z)
  args <- if (first) .fu_ordered(arg$args) else arg$args
  a <- args[[length(args)]]
  b <- .add(args[-length(args)])
  tb <- if (b$type == "Add") TR12(.fu_f("tan", b), first = FALSE) else .fu_f("tan", b)
  ta <- .fu_f("tan", a)
  .fu_div(.add(list(ta, tb)), .sub(S_One, .mul(list(ta, tb))))
})

#' @rdname fu
#' @export
TR13 <- function(rv) .fu_bu(rv, function(z) {
  if (z$type != "Mul") return(z)
  tn <- list(); ct <- list(); other <- list()
  for (a in .fu_margs(z)) {
    if (.fu_is(a, "tan")) tn[[length(tn) + 1L]] <- .fu_arg(a)
    else if (.fu_is(a, "cot")) ct[[length(ct) + 1L]] <- .fu_arg(a)
    else other[[length(other) + 1L]] <- a
  }
  if (length(tn) < 2L && length(ct) < 2L) return(z)
  while (length(tn) > 1L) {
    t1 <- tn[[length(tn)]]; t2 <- tn[[length(tn) - 1L]]; tn <- tn[seq_len(length(tn) - 2L)]
    ts <- .fu_f("tan", .add(list(t1, t2)))
    other[[length(other) + 1L]] <- .sub(S_One, .add(list(.fu_div(.fu_f("tan", t1), ts), .fu_div(.fu_f("tan", t2), ts))))
  }
  if (length(tn)) other[[length(other) + 1L]] <- .fu_f("tan", tn[[1L]])
  while (length(ct) > 1L) {
    t1 <- ct[[length(ct)]]; t2 <- ct[[length(ct) - 1L]]; ct <- ct[seq_len(length(ct) - 2L)]
    cs <- .fu_f("cot", .add(list(t1, t2)))
    other[[length(other) + 1L]] <- .add(list(S_One, .mul(list(.fu_f("cot", t1), cs)), .mul(list(.fu_f("cot", t2), cs))))
  }
  if (length(ct)) other[[length(other) + 1L]] <- .fu_f("cot", ct[[1L]])
  .mul(other)
})

#' @rdname fu
#' @export
TRmorrie <- function(rv) {
  f <- function(z, first = TRUE) {
    if (z$type != "Mul") return(z)
    if (first) { nd <- .fu_nd(z); return(.fu_div(f(nd[[1L]], FALSE), f(nd[[2L]], FALSE))) }
    args <- list(); argk <- list(); coss <- list(); other <- list()
    for (c in z$args) {
      be <- .base_exp(c)
      if (be[[2L]]$type == "Integer" && .fu_is(be[[1L]], "cos")) {
        sp <- as_coeff_mul_split(.fu_arg(be[[1L]]))
        k <- sp[[2L]]$key
        argk[[k]] <- sp[[2L]]
        args[[k]] <- c(args[[k]], list(sp[[1L]]))
        coss[[be[[1L]]$key]] <- be[[2L]]
      } else other[[length(other) + 1L]] <- c
    }
    new <- list()
    for (ak in names(args)) {
      a <- argk[[ak]]
      cl <- args[[ak]]
      cl <- cl[order(vapply(cl, function(v) num_value_dbl(v), 0))]
      while (length(cl)) {
        k <- 0L
        ci <- cl[[1L]]; cc <- ci
        inl <- function(v) any(vapply(cl, function(w) .fu_eq(w, v), TRUE))
        while (inl(cc)) { k <- k + 1L; cc <- .mul(list(cc, S_Two)) }
        if (k > 1L) {
          newarg <- .fu_div(.fu_div(.fu_f("sin", .mul(list(.pow(S_Two, .int(k)), ci, a))), .pow(S_Two, .int(k))), .fu_f("sin", .mul(list(ci, a))))
          take <- NULL; ccs <- list()
          for (i in seq_len(k)) {
            cc <- .div(cc, S_Two)
            key <- .fu_raw("cos", .mul(list(a, cc)))$key
            ccs[[length(ccs) + 1L]] <- cc
            take <- if (is.null(take)) coss[[key]] else if (.is_true(Lt(coss[[key]], take))) coss[[key]] else take
          }
          for (i in seq_len(k)) {
            cc <- ccs[[length(ccs)]]; ccs <- ccs[-length(ccs)]
            key <- .fu_raw("cos", .mul(list(a, cc)))$key
            coss[[key]] <- .sub(coss[[key]], take)
            if (is_zero_num(coss[[key]])) cl <- Filter(function(w) !.fu_eq(w, cc), cl)
          }
          new[[length(new) + 1L]] <- .pow(newarg, take)
        } else {
          b <- .fu_f("cos", .mul(list(cl[[1L]], a)))
          cl <- cl[-1L]
          other[[length(other) + 1L]] <- .pow(b, if (is.null(coss[[b$key]])) S_One else coss[[b$key]])
        }
      }
      args[[ak]] <- cl
    }
    if (length(new)) {
      rest <- unlist(lapply(names(args), function(ak) lapply(args[[ak]], function(k) .fu_raw("cos", .mul(list(k, argk[[ak]]))))), recursive = FALSE)
      z <- .mul(c(new, other, rest))
    }
    z
  }
  .fu_bu(rv, f)
}

# SymPy's as_f_sign_1: e = g*(f + s) with s = 1 or -1, or NULL.
.as_f_sign_1 <- function(e) {
  if (e$type != "Add" || length(e$args) != 2L) return(NULL)
  a <- e$args[[1L]]; b <- e$args[[2L]]
  if (.fu_eq(a, S_NegativeOne) || .fu_eq(a, S_One)) {
    g <- S_One
    if (b$type == "Mul" && is_number_type(b$args[[1L]]) && num_sign(b$args[[1L]]) < 0) { a <- .neg(a); b <- .neg(b); g <- .neg(g) }
    return(list(g, b, a))
  }
  nm <- .fac_normal(a, b)
  gcd <- .fac_gcd(a, b)
  ua <- nm[[1L]]; ub <- nm[[2L]]
  m1 <- S_NegativeOne$key
  if (!is.null(ua[[m1]])) { ua[[m1]] <- NULL; n1 <- -1L; n2 <- 1L }
  else if (!is.null(ub[[m1]])) { ub[[m1]] <- NULL; n1 <- 1L; n2 <- -1L }
  else { n1 <- 1L; n2 <- 1L }
  a <- .fac_expr(ua); b <- .fac_expr(ub)
  if (is_one_num(a)) { t <- a; a <- b; b <- t; t <- n1; n1 <- n2; n2 <- t }
  if (n1 == -1L) { gcd <- .neg(gcd); n2 <- -n2 }
  if (is_one_num(b)) return(list(gcd, a, .int(n2)))
  NULL
}

#' @rdname fu
#' @export
TR12i <- function(rv) .fu_bu(rv, function(rv) {
  if (!(rv$type %in% c("Add", "Mul", "Pow"))) return(rv)
  nd <- .fu_nd(rv); n <- nd[[1L]]; d <- nd[[2L]]
  if (!length(d$args) || !length(n$args)) return(rv)
  dok <- list(); dkey <- list()
  okd <- function(di) {
    m <- .as_f_sign_1(di)
    if (!is.null(m)) {
      g <- m[[1L]]; f <- m[[2L]]; s <- m[[3L]]
      if (.fu_eq(s, S_NegativeOne) && f$type == "Mul" && length(f$args) == 2L && all(vapply(f$args, .fu_is, TRUE, "tan")))
        return(list(g, f))
    }
    NULL
  }
  d_args <- .fu_margs(d)
  i <- 1L
  while (i <= length(d_args)) {
    di <- d_args[[i]]
    m <- okd(di)
    if (!is.null(m)) {
      s <- .add(lapply(m[[2L]]$args, .fu_arg))
      dok[[s$key]] <- S_One; dkey[[s$key]] <- s
      d_args[[i]] <- m[[1L]]
    } else if (di$type == "Add") {
      di <- sym_factor(di)
      if (di$type == "Mul") { d_args <- c(d_args, di$args); d_args[[i]] <- S_One }
    } else if (di$type == "Pow" && (.fu_is_int(di$args[[2L]]) || .fu_is_pos(di$args[[1L]]))) {
      m <- okd(di$args[[1L]])
      if (!is.null(m)) {
        s <- .add(lapply(m[[2L]]$args, .fu_arg))
        dok[[s$key]] <- di$args[[2L]]; dkey[[s$key]] <- s
        d_args[[i]] <- .pow(m[[1L]], di$args[[2L]])
      } else {
        di <- sym_factor(di)
        if (di$type == "Mul") { d_args <- c(d_args, di$args); d_args[[i]] <- S_One }
      }
    }
    i <- i + 1L
  }
  if (!length(dok)) return(rv)
  okn <- function(ni) {
    if (ni$type == "Add" && length(ni$args) == 2L && .fu_is(ni$args[[1L]], "tan") && .fu_is(ni$args[[2L]], "tan")) return(ni$args)
    NULL
  }
  n_args <- .fu_margs(factor_terms(n))
  hit <- FALSE
  i <- 1L
  while (i <= length(n_args)) {
    ni <- n_args[[i]]
    m <- okn(ni)
    skip <- FALSE
    if (is.null(m)) {
      m <- okn(.neg(ni))
      if (!is.null(m)) n_args[[i]] <- S_NegativeOne
      else if (ni$type == "Add") {
        ni <- sym_factor(ni)
        if (ni$type == "Mul") { n_args <- c(n_args, ni$args); n_args[[i]] <- S_One }
        skip <- TRUE
      } else if (ni$type == "Pow" && (.fu_is_int(ni$args[[2L]]) || .fu_is_pos(ni$args[[1L]]))) {
        m <- okn(ni$args[[1L]])
        if (!is.null(m)) n_args[[i]] <- S_One
        else {
          ni <- sym_factor(ni)
          if (ni$type == "Mul") { n_args <- c(n_args, ni$args); n_args[[i]] <- S_One }
          skip <- TRUE
        }
      } else skip <- TRUE
    } else n_args[[i]] <- S_One
    if (!skip) {
      s <- .add(lapply(m, .fu_arg))
      if (!is.null(dok[[s$key]])) {
        hit <- TRUE
        ed <- dok[[s$key]]
        newed <- .sub(ed, S_One)
        if (is_zero_num(newed)) dok[[s$key]] <- NULL else dok[[s$key]] <- newed
        n_args[[i]] <- .mul(list(n_args[[i]], .neg(.fu_f("tan", s))))
      }
    }
    i <- i + 1L
  }
  if (hit) {
    rest <- lapply(names(dok), function(k) .pow(.sub(.add(lapply(dkey[[k]]$args, function(a) .fu_f("tan", a))), S_One), dok[[k]]))
    rv <- .fu_div(.fu_div(.mul(n_args), .mul(d_args)), .mul(rest))
  }
  rv
})

#' @rdname fu
#' @export
TR14 <- function(rv, first = TRUE) .fu_bu(rv, function(z) {
  if (z$type != "Mul") return(z)
  if (first) {
    nd <- .fu_nd(z)
    if (!is_one_num(nd[[2L]])) {
      nn <- TR14(nd[[1L]], first = FALSE); dd <- TR14(nd[[2L]], first = FALSE)
      if (!.fu_eq(nn, nd[[1L]]) || !.fu_eq(dd, nd[[2L]])) z <- .fu_div(nn, dd)
      return(z)
    }
  }
  other <- list(); process <- list()
  for (a in z$args) {
    e <- S_One
    if (a$type == "Pow") {
      b <- a$args[[1L]]; e <- a$args[[2L]]
      if (!(.fu_is_int(e) || .fu_is_pos(b))) { other[[length(other) + 1L]] <- a; next }
      a <- b
    }
    m <- .as_f_sign_1(a)
    if (is.null(m) || !.fu_is(m[[2L]], c("cos", "sin"))) {
      other[[length(other) + 1L]] <- if (is_one_num(e)) a else .pow(a, e)
      next
    }
    process[[length(process) + 1L]] <- list(g = m[[1L]], num = is_number_type(e), e = e, f = m[[2L]], si = m[[3L]], a = a)
  }
  process <- .fu_ordered_tuples(process)
  nother <- length(other)
  trig2 <- function(A, B, take) {
    t <- if (A$f$name == "cos") "sin" else "cos"
    .pow(.mul(list(.neg(A$g), B$g, .pow(.fu_f(t, .fu_arg(A$f)), S_Two))), take)
  }
  while (length(process)) {
    A <- process[[1L]]; process <- process[-1L]
    if (length(process)) {
      B <- process[[1L]]
      if (A$num && B$num) {
        if (.fu_eq(A$f, B$f) && !.fu_eq(A$si, B$si)) {
          process <- process[-1L]
          take <- if (.is_true(Lt(B$e, A$e))) B$e else A$e
          if (!.fu_eq(B$e, take)) { rem <- B; rem$e <- .sub(rem$e, take); process <- c(list(rem), process) }
          else if (!.fu_eq(A$e, take)) { rem <- A; rem$e <- .sub(rem$e, take); process <- c(list(rem), process) }
          other[[length(other) + 1L]] <- trig2(A, B, take)
          next
        }
      } else if (.fu_eq(A$e, B$e)) {
        if (.fu_eq(A$f, B$f) && !.fu_eq(A$si, B$si)) {
          process <- process[-1L]
          other[[length(other) + 1L]] <- trig2(A, B, A$e)
          next
        }
      }
    }
    other[[length(other) + 1L]] <- .pow(A$a, A$e)
  }
  if (length(other) != nother) z <- .mul(other)
  z
})

# ordered() of the TR14 tuples (g, e.is_Number, e, f, si, a), compared in turn
.fu_ordered_tuples <- function(ps) {
  if (length(ps) < 2L) return(ps)
  key <- function(p) list(.sort_key(p$g), .sort_key(if (p$num) S_true else S_false), .sort_key(p$e), .sort_key(p$f), .sort_key(p$si), .sort_key(p$a))
  cmp <- function(x, y) { for (i in seq_along(x)) { c <- .key_cmp(x[[i]], y[[i]]); if (c != 0) return(c) }; 0 }
  ks <- lapply(ps, key)
  o <- integer()
  for (i in seq_along(ks)) {
    j <- length(o)
    while (j > 0L && cmp(ks[[o[j]]], ks[[i]]) > 0) j <- j - 1L
    o <- append(o, i, after = j)
  }
  ps[o]
}

#' @rdname fu
#' @export
TR15 <- function(rv, max = 4, pow = FALSE) .fu_bu(rv, function(z) {
  if (!(z$type == "Pow" && .fu_is(z$args[[1L]], "sin"))) return(z)
  e <- z$args[[2L]]
  if (e$type == "Integer" && .Z(e) %% 2 == 1) return(.fu_div(TR15(.pow(z$args[[1L]], .add(list(e, S_One)))), z$args[[1L]]))
  ia <- .pow(z, S_NegativeOne)
  a <- .TR56(ia, "sin", "cot", function(x) .add(list(S_One, x)), max, pow)
  if (!.fu_eq(a, ia)) z <- a
  z
})

#' @rdname fu
#' @export
TR16 <- function(rv, max = 4, pow = FALSE) .fu_bu(rv, function(z) {
  if (!(z$type == "Pow" && .fu_is(z$args[[1L]], "cos"))) return(z)
  e <- z$args[[2L]]
  if (e$type == "Integer" && .Z(e) %% 2 == 1) return(.fu_div(TR15(.pow(z$args[[1L]], .add(list(e, S_One)))), z$args[[1L]]))
  ia <- .pow(z, S_NegativeOne)
  a <- .TR56(ia, "cos", "tan", function(x) .add(list(S_One, x)), max, pow)
  if (!.fu_eq(a, ia)) z <- a
  z
})

#' @rdname fu
#' @export
TR111 <- function(rv) .fu_bu(rv, function(z) {
  if (!(z$type == "Pow" && (.fu_is_pos(z$args[[1L]]) || (.fu_is_int(z$args[[2L]]) && isTRUE(ask_fact(z$args[[2L]], "negative")))))) return(z)
  b <- z$args[[1L]]; e <- z$args[[2L]]
  if (.fu_is(b, "tan")) return(.pow(.fu_f("cot", .fu_arg(b)), .neg(e)))
  if (.fu_is(b, "sin")) return(.pow(.fu_f("csc", .fu_arg(b)), .neg(e)))
  if (.fu_is(b, "cos")) return(.pow(.fu_f("sec", .fu_arg(b)), .neg(e)))
  z
})

#' @rdname fu
#' @export
TR22 <- function(rv, max = 4, pow = FALSE) .fu_bu(rv, function(z) {
  if (!(z$type == "Pow" && .fu_is(z$args[[1L]], c("cot", "tan")))) return(z)
  z <- .TR56(z, "tan", "sec", function(x) .sub(x, S_One), max, pow)
  .TR56(z, "cot", "csc", function(x) .sub(x, S_One), max, pow)
})

#' @rdname fu
#' @export
TRpower <- function(rv) .fu_bu(rv, function(z) {
  if (!(z$type == "Pow" && .fu_is(z$args[[1L]], c("sin", "cos")))) return(z)
  b <- z$args[[1L]]; n <- z$args[[2L]]
  if (!(n$type == "Integer" && .int_val(n) > 0)) return(z)
  x <- .fu_arg(b); nv <- .int_val(n)
  bin <- function(k) sym_binomial(.int(nv), .int(k))
  pre <- .pow(S_Two, .int(1L - nv))
  if (nv %% 2 == 1 && b$name == "cos")
    z <- .mul(list(pre, .add(lapply(seq_len((nv + 1L) %/% 2L) - 1L, function(k) .mul(list(bin(k), .fu_f("cos", .mul(list(.int(nv - 2L*k), x)))))))))
  else if (nv %% 2 == 1)
    z <- .mul(list(pre, .pow(S_NegativeOne, .int((nv - 1L) %/% 2L)), .add(lapply(seq_len((nv + 1L) %/% 2L) - 1L, function(k)
      .mul(list(bin(k), .pow(S_NegativeOne, .int(k)), .fu_f("sin", .mul(list(.int(nv - 2L*k), x)))))))))
  else if (b$name == "cos")
    z <- .mul(list(pre, .add(lapply(seq_len(nv %/% 2L) - 1L, function(k) .mul(list(bin(k), .fu_f("cos", .mul(list(.int(nv - 2L*k), x)))))))))
  else
    z <- .mul(list(pre, .pow(S_NegativeOne, .int(nv %/% 2L)), .add(lapply(seq_len(nv %/% 2L) - 1L, function(k)
      .mul(list(bin(k), .pow(S_NegativeOne, .int(k)), .fu_f("cos", .mul(list(.int(nv - 2L*k), x)))))))))
  if (nv %% 2 == 0) z <- .add(list(z, .mul(list(.pow(S_Two, .int(-nv)), sym_binomial(.int(nv), .int(nv %/% 2L))))))
  z
})


# SymPy's process_common_addends: group the terms by their absolute
# coefficient (and key2), apply do to each group, and add the results.
.process_common_addends <- function(rv, do, key2 = NULL) {
  groups <- list(); gc <- list()
  for (a in rv$args) {
    sp <- as_coeff_mul_split(a)
    c <- sp[[1L]]; t <- sp[[2L]]
    if (num_sign(c) < 0) { c <- .neg(c); t <- .neg(t) }
    k <- paste(c$key, if (!is.null(key2)) key2(t) else "1", sep = "|")
    groups[[k]] <- c(groups[[k]], list(t)); gc[[k]] <- c
  }
  args <- list(); hit <- FALSE
  for (k in names(groups)) {
    v <- groups[[k]]
    if (length(v) > 1L) {
      e <- .add(v)
      new <- do(e)
      if (!.fu_eq(new, e)) { e <- new; hit <- TRUE }
      args[[length(args) + 1L]] <- .mul(list(gc[[k]], e))
    } else args[[length(args) + 1L]] <- .mul(list(gc[[k]], v[[1L]]))
  }
  if (hit) .add(args) else rv
}

# SymPy's trig_split(a, b, two): the pieces of a + b when it is a sum or
# difference formula, as list(gcd, n1, n2, a, b, flag), or NULL.
.trig_split <- function(a, b, two = FALSE) {
  nm <- .fac_normal(a, b)
  gcd <- .fac_gcd(a, b)
  ua <- nm[[1L]]; ub <- nm[[2L]]
  n1 <- 1L; n2 <- 1L
  m1 <- S_NegativeOne$key
  if (!is.null(ua[[m1]])) { ua[[m1]] <- NULL; n1 <- -n1 }
  else if (!is.null(ub[[m1]])) { ub[[m1]] <- NULL; n2 <- -n2 }
  a <- .fac_expr(ua); b <- .fac_expr(ub)
  pcs <- function(a) {
    c <- NULL; s <- NULL; co <- S_One
    if (a$type == "Mul") {
      sp <- as_coeff_mul_split(a); co <- sp[[1L]]; a <- sp[[2L]]
      if ((a$type == "Mul" && length(a$args) > 2L) || !two) return(NULL)
      args <- .fu_margs(a)
      a <- args[[1L]]; args <- args[-1L]
      if (.fu_is(a, "cos")) c <- a
      else if (.fu_is(a, "sin")) s <- a
      else if (a$type == "Pow" && identical(a$args[[2L]]$key, S_Half$key)) co <- .mul(list(co, a))
      else return(NULL)
      if (length(args)) {
        b <- args[[1L]]
        if (.fu_is(b, "cos")) { if (!is.null(c)) s <- b else c <- b }
        else if (.fu_is(b, "sin")) { if (!is.null(s)) c <- b else s <- b }
        else if (b$type == "Pow" && identical(b$args[[2L]]$key, S_Half$key)) co <- .mul(list(co, b))
        else return(NULL)
      }
      return(list(co = if (is_one_num(co)) NULL else co, c = c, s = s))
    }
    if (.fu_is(a, "cos")) c <- a else if (.fu_is(a, "sin")) s <- a
    if (is.null(c) && is.null(s)) return(NULL)
    list(co = NULL, c = c, s = s)
  }
  ma <- pcs(a); if (is.null(ma)) return(NULL)
  mb <- pcs(b); if (is.null(mb)) return(NULL)
  coa <- ma$co; ca <- ma$c; sa <- ma$s; cob <- mb$co; cb <- mb$c; sb <- mb$s
  if ((is.null(ca) && !is.null(cb)) || (!is.null(ca) && .fu_is(ca, "sin"))) {
    t <- list(coa, ca, sa); coa <- cob; ca <- cb; sa <- sb; cob <- t[[1L]]; cb <- t[[2L]]; sb <- t[[3L]]
    t <- n1; n1 <- n2; n2 <- t
  }
  if (!two) {
    c <- if (!is.null(ca)) ca else sa
    s <- if (!is.null(cb)) cb else sb
    if (is.null(c) || is.null(s) || c$name != s$name) return(NULL)
    return(list(gcd = gcd, n1 = n1, n2 = n2, a = .fu_arg(c), b = .fu_arg(s), flag = c$name == "cos"))
  }
  if (is.null(coa) && is.null(cob)) {
    if (!is.null(ca) && !is.null(cb) && !is.null(sa) && !is.null(sb)) {
      if ((ca$name == sa$name) != (cb$name == sb$name)) return(NULL)
      ak <- c(.fu_arg(ca)$key, .fu_arg(sa)$key)
      if (!all(c(.fu_arg(cb)$key, .fu_arg(sb)$key) %in% ak)) return(NULL)
      return(list(gcd = gcd, n1 = n1, n2 = n2, a = .fu_arg(ca), b = .fu_arg(sa), flag = ca$name == sa$name))
    }
  }
  if ((!is.null(ca) && !is.null(sa)) || (!is.null(cb) && !is.null(sb)) ||
      (two && ((is.null(ca) && is.null(sa)) || (is.null(cb) && is.null(sb))))) return(NULL)
  c <- if (!is.null(ca)) ca else sa
  s <- if (!is.null(cb)) cb else sb
  if (!.fu_eq(.fu_arg(c), .fu_arg(s))) return(NULL)
  if (is.null(coa)) coa <- S_One
  if (is.null(cob)) cob <- S_One
  if (.fu_eq(coa, cob)) return(list(gcd = .mul(list(gcd, .sqrt(S_Two))), n1 = n1, n2 = n2, a = .fu_arg(c), b = .div(S_Pi, .int(4)), flag = FALSE))
  r <- .fu_div(coa, cob)
  if (.fu_eq(r, .sqrt(.int(3)))) return(list(gcd = .mul(list(gcd, S_Two, cob)), n1 = n1, n2 = n2, a = .fu_arg(c), b = .div(S_Pi, .int(3)), flag = FALSE))
  if (.fu_eq(r, .pow(.sqrt(.int(3)), S_NegativeOne))) return(list(gcd = .mul(list(gcd, S_Two, coa)), n1 = n1, n2 = n2, a = .fu_arg(c), b = .div(S_Pi, .int(6)), flag = FALSE))
  NULL
}

#' @rdname fu
#' @export
fu <- function(rv, measure = .fu_measure) {
  was <- as_sym(rv)
  rv <- TR1(was)
  RL1 <- .fu_chain(TR4, TR3, TR4, TR12, TR4, TR13, TR4, TR0)
  CTR1 <- .fu_choice(.fu_chain(TR5, TR0), .fu_chain(TR6, TR0), identity)
  CTR2 <- .fu_chain(TR11, .fu_choice(.fu_chain(TR5, TR0), .fu_chain(TR6, TR0), TR0))
  CTR3 <- .fu_choice(.fu_chain(TRmorrie, TR8, TR0), .fu_chain(TRmorrie, TR8, TR10i, TR0), identity)
  CTR4 <- .fu_choice(.fu_chain(TR4, TR10i), identity)
  RL2 <- .fu_choice(.fu_chain(TR4, TR3, TR10, TR4, TR3, TR11), .fu_chain(TR5, TR7, TR11, TR4),
                    .fu_chain(CTR3, CTR1, TR9, CTR2, TR4, TR9, TR9, CTR4), identity)
  if (.fu_has(rv, c("tan", "cot"))) {
    rv1 <- .fu_apply(RL1, rv, measure)
    if (.fu_less(measure(rv1), measure(rv))) rv <- rv1
    if (.fu_has(rv, c("tan", "cot"))) rv <- TR2(rv)
  }
  if (.fu_has(rv, c("sin", "cos"))) {
    rv1 <- .fu_apply(RL2, rv, measure)
    rv2 <- TR8(TRmorrie(rv1))
    rv <- .fu_min(list(was, rv, rv1, rv2), measure)
  }
  .fu_min(list(TR2i(rv), rv), measure)
}
