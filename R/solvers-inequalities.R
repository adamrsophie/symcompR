# Polynomial, rational and absolute value inequalities after
# sympy/solvers/inequalities.py: solve_poly_inequality,
# solve_poly_inequalities, solve_rational_inequalities,
# reduce_rational_inequalities, reduce_abs_inequality(ies) and
# reduce_inequalities.

.ineq_poly_parts <- function(p) {
  if (inherits(p, "Poly")) return(list(expr = as_expr(p), gen = p$gens[[1L]]))
  e <- as_sym(p)
  fs <- free_symbols(e)
  list(expr = e, gen = if (length(fs)) fs[[1L]] else NULL)
}

.ineq_rel <- function(a, b, rel) {
  switch(rel, "==" = Eq(a, b), "!=" = Ne(a, b), "<" = Lt(a, b), "<=" = Le(a, b), ">" = Gt(a, b), ">=" = Ge(a, b),
         .sp_stop("ValueError", paste0("'", rel, "' is not a valid relation")))
}

# Poly.real_roots(multiple=False): the real roots in increasing order with
# their multiplicities, radicals for linear, quadratic and binomial factors
# and CRootOf for the others (CRootOf.real_roots(radicals=True)).
.ineq_real_roots <- function(f, x) {
  fl <- factor_list(f, x)
  facs <- fl[[2L]]
  out <- list()
  for (pr in facs) {
    g <- pr[[1L]]
    m <- as.integer(pr[[2L]])
    d <- degree(g, x)
    rs <- NULL
    if (d == 1) {
      rs <- .roots_sympy(g, x)
    } else if (d == 2 || length(.aargs(expand(g))) == 2L) {
      rr <- .roots_sympy(g, x)
      rs <- rr[seq_along(rr)]
    }
    if (!is.null(rs)) {
      rs <- Filter(function(r) isTRUE(.sp_fact(r, "real")), rs[seq_along(rs)])
    } else {
      nreal <- count_roots(g, x)
      rs <- lapply(seq_len(nreal) - 1L, function(k) CRootOf(g, x, k))
    }
    for (r in rs) out[[length(out) + 1L]] <- list(r, m)
  }
  if (length(out) > 1L) {
    # exact comparison for rational roots, numeric otherwise
    less <- function(a, b) {
      d <- .sub(a, b)
      if (is_exact_number(d)) return(num_sign(d) < 0)
      as.numeric(Re(as.complex(.neval(d, list())))) < 0
    }
    o <- seq_along(out)
    for (i in seq_along(o)[-1L]) {
      j <- i
      while (j > 1L && less(out[[o[j]]][[1L]], out[[o[j - 1L]]][[1L]])) {
        tmp <- o[j]; o[j] <- o[j - 1L]; o[j - 1L] <- tmp
        j <- j - 1L
      }
    }
    out <- out[o]
  }
  out
}

#' Inequalities in one real variable
#'
#' These follow SymPy's `sympy.solvers.inequalities`.
#' `solve_poly_inequality(poly, rel)` solves `poly rel 0` for a [Poly()]
#' and returns a list of intervals. `solve_poly_inequalities(polys)` joins
#' the solutions of a list of `list(poly, rel)` pairs.
#' `solve_rational_inequalities(eqs)` solves `eqs`, a list of lists whose
#' entries `list(list(numer, denom), rel)` mean `numer/denom rel 0`; inner
#' lists are joined by `And`, the outer list by `Or`.
#' `reduce_rational_inequalities(exprs, gen)` does the same for relationals
#' or `list(expr, rel)` pairs. `reduce_abs_inequality(expr, rel, gen)` and
#' `reduce_abs_inequalities(exprs, gen)` handle nested absolute values, and
#' `reduce_inequalities(inequalities, symbols)` reduces a system of
#' inequalities.
#'
#' @param poly A `Poly` (or a polynomial expression).
#' @param rel A relation: `"<"`, `"<="`, `">"`, `">="`, `"=="` or `"!="`.
#' @param polys A list of `list(poly, rel)` pairs.
#' @param eqs A nested list of inequalities.
#' @param exprs A nested list (rational) or a list of `list(expr, rel)` pairs
#'   (absolute values).
#' @param gen The real variable.
#' @param relational Return a relational instead of a set.
#' @param expr An expression.
#' @param inequalities A relational or a list of them.
#' @param symbols The symbols to solve for.
#' @return A list of sets, a set or a relational.
#' @examples
#' x <- Symbol("x")
#' solve_poly_inequality(Poly(x^2 - 1, x), ">")
#' reduce_rational_inequalities(list(list(x^2 < 4)), x)
#' reduce_abs_inequality(Abs(x - 5) - 3, "<", x)
#' reduce_inequalities(list(x^2 < 4, x > 0), x)
#' @name inequality-solvers
NULL

#' @rdname inequality-solvers
#' @export
solve_poly_inequality <- function(poly, rel) {
  pp <- .ineq_poly_parts(poly)
  e <- pp$expr
  x <- pp$gen
  if (.sp_is_number(e)) {
    t <- .ineq_rel(e, S_Zero, rel)
    if (.is_true(t)) return(list(Reals))
    if (.is_false(t)) return(list(S_EmptySet))
    .sp_stop("NotImplementedError", paste0("could not determine truth value of ", .sp_py_str(t)))
  }
  reals <- .ineq_real_roots(e, x)
  intervals <- list()
  if (rel == "==") {
    for (r in reals) intervals[[length(intervals) + 1L]] <- Interval(r[[1L]], r[[1L]])
  } else if (rel == "!=") {
    left <- S_NegativeInfinity
    for (r in c(reals, list(list(S_Infinity, 1L)))) {
      intervals[[length(intervals) + 1L]] <- Interval(left, r[[1L]], TRUE, TRUE)
      left <- r[[1L]]
    }
  } else {
    sign <- if (.is_true(Gt(LC(e, x), S_Zero))) 1L else -1L
    eq_sign <- switch(rel, ">" = 1L, "<" = -1L, ">=" = 1L, "<=" = -1L,
                      .sp_stop("ValueError", paste0("'", rel, "' is not a valid relation")))
    equal <- rel %in% c(">=", "<=")
    right <- S_Infinity
    right_open <- TRUE
    for (r in rev(reals)) {
      left <- r[[1L]]
      if (r[[2L]] %% 2L) {
        if (sign == eq_sign) intervals <- c(list(Interval(left, right, !equal, right_open)), intervals)
        sign <- -sign
        right <- left
        right_open <- !equal
      } else {
        if (sign == eq_sign && !equal) {
          intervals <- c(list(Interval(left, right, TRUE, right_open)), intervals)
          right <- left
          right_open <- TRUE
        } else if (sign != eq_sign && equal) {
          intervals <- c(list(Interval(left, left)), intervals)
        }
      }
    }
    if (sign == eq_sign) intervals <- c(list(Interval(S_NegativeInfinity, right, TRUE, right_open)), intervals)
  }
  intervals
}

#' @rdname inequality-solvers
#' @export
solve_poly_inequalities <- function(polys) {
  sets <- list()
  for (p in polys) sets <- c(sets, solve_poly_inequality(p[[1L]], p[[2L]]))
  do.call(Union, sets)
}

#' @rdname inequality-solvers
#' @export
solve_rational_inequalities <- function(eqs) {
  result <- S_EmptySet
  for (grp in eqs) {
    if (!length(grp)) next
    global <- list(Interval(S_NegativeInfinity, S_Infinity))
    for (q in grp) {
      np <- .ineq_poly_parts(q[[1L]][[1L]])
      dp <- .ineq_poly_parts(q[[1L]][[2L]])
      rel <- q[[2L]]
      x <- if (!is.null(np$gen)) np$gen else if (!is.null(dp$gen)) dp$gen else Symbol("x")
      nd <- expand(.mul(list(np$expr, dp$expr)))
      numer_intervals <- solve_poly_inequality(Poly(nd, x), rel)
      denom_intervals <- solve_poly_inequality(Poly(dp$expr, x), "==")
      intervals <- list()
      for (ni in numer_intervals) for (gi in global) {
        iv <- Intersection(ni, gi)
        if (!identical(iv$key, S_EmptySet$key)) intervals[[length(intervals) + 1L]] <- iv
      }
      global <- intervals
      intervals <- list()
      for (gi in global) {
        for (di in denom_intervals) gi <- Complement(gi, di)
        if (!identical(gi$key, S_EmptySet$key)) intervals[[length(intervals) + 1L]] <- gi
      }
      global <- intervals
      if (!length(global)) break
    }
    for (iv in global) result <- Union(result, iv)
  }
  result
}

# Set.evalf(): the numbers in a set evaluated, the infinities kept
.ineq_set_evalf <- function(s) {
  ev <- function(z) if (is_inf_const(z)) z else N(z)
  switch(s$type,
         FiniteSet = do.call(FiniteSet, lapply(s$args, ev)),
         Interval = Interval(ev(s$args[[1L]]), ev(s$args[[2L]]), .is_true(s$args[[3L]]), .is_true(s$args[[4L]])),
         Union = do.call(Union, lapply(s$args, .ineq_set_evalf)),
         s)
}

# Poly.to_exact(): each Float as the simplest rational within its precision,
# taken from its shortest decimal form
.ineq_to_exact <- function(e) {
  .map_bottom_up(e, function(z) {
    if (z$type != "Float") return(z)
    s <- .fmt_float_sym(z, strip = TRUE)
    m <- regmatches(s, regexec("^(-?)([0-9]*)\\.?([0-9]*)(e([-+]?[0-9]+))?$", s))[[1L]]
    if (!length(m)) return(nsimplify(z, rational = TRUE))
    digits <- paste0(m[3L], m[4L])
    k <- nchar(m[4L]) - (if (nzchar(m[6L])) as.integer(m[6L]) else 0L)
    q <- gmp::as.bigq(gmp::as.bigz(digits))
    q <- if (k >= 0L) q / gmp::as.bigz(10)^k else q * gmp::as.bigz(10)^(-k)
    make_rational(if (nzchar(m[2L])) -q else q)
  })
}

.ineq_poly_error <- function()
  stop(.sympy_condition("PolynomialError", .filldedent("only polynomials and rational functions are supported in this context.")))

# coefficients of e as a polynomial in gen, or an error
.ineq_coeffs <- function(e, gen) {
  if (!isTRUE(is_polynomial(e, gen))) .ineq_poly_error()
  p <- .sp_as_poly(e, list(gen))
  if (is.null(p)) .ineq_poly_error()
  lapply(p$terms, `[[`, "coeff")
}

.ineq_lost_factors <- function(a, b, gen) {
  a <- .ineq_to_exact(a); b <- .ineq_to_exact(b)
  ln <- .as_numer_denom(a); rn <- .as_numer_denom(b)
  for (e in c(ln, rn)) .ineq_coeffs(e, gen)
  ld <- ln[[2L]]; rd <- rn[[2L]]
  if (!.has_free(ld, list(gen)) || !.has_free(rd, list(gen))) return(S_One)
  common <- sqf_part(gcd(ld, rd, gen), gen)
  if (!.has_free(common, list(gen))) return(S_One)
  dd <- .as_numer_denom(.sub(a, b))[[2L]]
  remaining <- sqf_part(dd, gen)
  exquo(common, gcd(common, remaining, gen), gen)
}

#' @rdname inequality-solvers
#' @export
reduce_rational_inequalities <- function(exprs, gen, relational = TRUE) {
  gen <- as_sym(gen)
  exact <- TRUE
  solution <- S_EmptySet
  for (ex in exprs) {
    if (!length(ex)) next
    eqs <- list()
    sol <- Reals
    for (item in ex) {
      rel <- NULL
      if (is.list(item) && !is_sym(item)) {
        expr <- as_sym(item[[1L]])
        rel <- item[[2L]]
        nd <- .as_numer_denom(expr)
      } else {
        expr <- as_sym(item)
        if (is_relational(expr) && expr$type != "BooleanAtom") {
          lost <- .ineq_lost_factors(expr$args[[1L]], expr$args[[2L]], gen)
          if (!is_one_num(lost)) {
            # the domain of the lost factors is inexact when a denominator has a Float
            if (.has_type(lost, "Float") || any(vapply(expr$args, function(z) .has_type(.as_numer_denom(z)[[2L]], "Float"), TRUE))) {
              lost <- .ineq_to_exact(lost)
              exact <- FALSE
            }
            eqs[[length(eqs) + 1L]] <- list(list(lost, S_One), "!=")
          }
          rel <- .rel_op(expr)
          nd <- .as_numer_denom(.sub(expr$args[[1L]], expr$args[[2L]]))
        }
      }
      if (is.null(rel)) {
        rel <- "=="
        nd <- if (.is_true(expr)) list(S_Zero, S_One) else if (.is_false(expr)) list(S_One, S_One) else .as_numer_denom(expr)
      }
      numer <- nd[[1L]]; den <- nd[[2L]]
      cf <- c(.ineq_coeffs(numer, gen), .ineq_coeffs(den, gen))
      if (any(vapply(cf, function(c) .has_type(c, "Float"), TRUE))) {
        numer <- .ineq_to_exact(numer)
        den <- .ineq_to_exact(den)
        cf <- c(.ineq_coeffs(numer, gen), .ineq_coeffs(den, gen))
        exact <- FALSE
      }
      if (!all(vapply(cf, is_exact_number, TRUE))) {
        sol <- Intersection(sol, solve_univariate_inequality(.ineq_rel(.div(numer, den), S_Zero, rel), gen,
                                                             relational = FALSE))
      } else eqs[[length(eqs) + 1L]] <- list(list(numer, den), rel)
    }
    if (length(eqs)) sol <- Intersection(sol, solve_rational_inequalities(list(lapply(eqs, function(q)
      list(list(Poly(q[[1L]][[1L]], gen), Poly(q[[1L]][[2L]], gen)), q[[2L]])))))
    solution <- Union(solution, sol)
  }
  if (!exact && !identical(solution$key, S_EmptySet$key)) solution <- .ineq_set_evalf(solution)
  if (relational) solution <- as_relational(solution, gen)
  solution
}

#' @rdname inequality-solvers
#' @export
reduce_abs_inequality <- function(expr, rel, gen) {
  gen <- as_sym(gen)
  if (isFALSE(.sp_fact(gen, "extended_real")))
    .sp_stop("TypeError", .filldedent("Cannot solve inequalities with absolute values containing non-real variables."))
  scan <- function(e) {
    out <- list()
    if (e$type %in% c("Add", "Mul")) {
      for (a in e$args) {
        sub <- scan(a)
        if (!length(out)) {
          out <- sub
        } else {
          nxt <- list()
          for (p in out) for (q in sub)
            nxt[[length(nxt) + 1L]] <- list(if (e$type == "Add") .add(list(p[[1L]], q[[1L]])) else .mul(list(p[[1L]], q[[1L]])),
                                            c(p[[2L]], q[[2L]]))
          out <- nxt
        }
      }
    } else if (e$type == "Pow") {
      n <- e$args[[2L]]
      if (n$type != "Integer") .sp_stop("ValueError", "Only Integer Powers are allowed on Abs.")
      out <- lapply(scan(e$args[[1L]]), function(p) list(.pow(p[[1L]], n), p[[2L]]))
    } else if (e$type == "Function" && e$name == "Abs") {
      for (p in scan(e$args[[1L]])) {
        out[[length(out) + 1L]] <- list(p[[1L]], c(p[[2L]], list(list(p[[1L]], ">="))))
        out[[length(out) + 1L]] <- list(.neg(p[[1L]]), c(p[[2L]], list(list(p[[1L]], "<"))))
      }
    } else out <- list(list(e, list()))
    out
  }
  mapping <- c("<" = ">", "<=" = ">=")
  inequalities <- list()
  for (p in scan(as_sym(expr))) {
    e <- if (!(rel %in% names(mapping))) list(p[[1L]], rel) else list(.neg(p[[1L]]), mapping[[rel]])
    inequalities[[length(inequalities) + 1L]] <- c(list(e), p[[2L]])
  }
  reduce_rational_inequalities(inequalities, gen)
}

#' @rdname inequality-solvers
#' @export
reduce_abs_inequalities <- function(exprs, gen) {
  do.call(And, lapply(exprs, function(q) reduce_abs_inequality(q[[1L]], q[[2L]], gen)))
}

.ineq_symbol_atoms <- function(e) Filter(function(z) z$type %in% c("Symbol", "Dummy"), atoms(e))

.reduce_inequalities <- function(inequalities, symbols) {
  poly_part <- list(); poly_gens <- list()
  abs_part <- list(); abs_gens <- list()
  other <- list()
  add_to <- function(lst, gens, gen, item) {
    k <- gen$key
    if (is.null(lst[[k]])) { lst[[k]] <- list(); gens[[k]] <- gen }
    lst[[k]][[length(lst[[k]]) + 1L]] <- item
    list(lst, gens)
  }
  for (ie in inequalities) {
    expr <- .sub(ie$args[[1L]], ie$args[[2L]])
    rel <- .rel_op(ie)
    gens <- .ineq_symbol_atoms(ie)
    if (length(gens) == 1L) {
      gen <- gens[[1L]]
    } else {
      common <- Filter(function(s) .sp_in(s, symbols), free_symbols(ie))
      if (length(common) == 1L) {
        other[[length(other) + 1L]] <- .solve_inequality(ie, common[[1L]])
        next
      }
      .sp_stop("NotImplementedError", .filldedent("inequality has more than one symbol of interest."))
    }
    bden <- .sp_uniq(c(denoms(ie$args[[1L]], gen), denoms(ie$args[[2L]], gen)))
    eden <- denoms(expr, gen)
    if (any(vapply(bden, function(d) !.sp_in(d, eden), TRUE))) {
      other[[length(other) + 1L]] <- .solve_inequality(ie, gen)
      next
    }
    if (isTRUE(is_polynomial(expr, gen))) {
      r <- add_to(poly_part, poly_gens, gen, list(expr, rel)); poly_part <- r[[1L]]; poly_gens <- r[[2L]]
    } else {
      components <- .sp_find(expr, function(u) .has1(u, gen) &&
                               (u$type == "Function" || (u$type == "Pow" && u$args[[2L]]$type != "Integer")))
      if (length(components) && all(vapply(components, function(i) i$type == "Function" && i$name == "Abs", TRUE))) {
        r <- add_to(abs_part, abs_gens, gen, list(expr, rel)); abs_part <- r[[1L]]; abs_gens <- r[[2L]]
      } else {
        other[[length(other) + 1L]] <- .solve_inequality(ie, gen)
      }
    }
  }
  poly_reduced <- lapply(names(poly_part), function(k) reduce_rational_inequalities(list(poly_part[[k]]), poly_gens[[k]]))
  abs_reduced <- lapply(names(abs_part), function(k) reduce_abs_inequalities(abs_part[[k]], abs_gens[[k]]))
  do.call(And, c(poly_reduced, abs_reduced, other))
}

#' @rdname inequality-solvers
#' @export
reduce_inequalities <- function(inequalities, symbols = list()) {
  if (!is.list(inequalities) || is_sym(inequalities)) inequalities <- list(inequalities)
  inequalities <- lapply(inequalities, as_sym)
  gens <- unique_syms(lapply(inequalities, free_symbols))
  if (!is.list(symbols) || is_sym(symbols)) symbols <- list(symbols)
  symbols <- lapply(symbols, as_sym)
  symbols <- if (length(symbols)) Filter(function(s) .sp_in(s, gens), symbols) else gens
  if (any(vapply(symbols, function(i) isFALSE(.sp_fact(i, "extended_real")), TRUE)))
    .sp_stop("TypeError", .filldedent("inequalities cannot contain symbols that are not real."))
  recast <- list()
  for (i in gens) if (is.na(.sp_fact(i, "extended_real"))) recast[[length(recast) + 1L]] <- list(i, Dummy(i$name, extended_real = TRUE))
  if (length(recast)) {
    # xreplace rebuilds each relational, which evaluates it again
    reval <- function(i) if (is_relational(i)) .ineq_rel(i$args[[1L]], i$args[[2L]], .rel_op(i)) else i
    inequalities <- lapply(inequalities, function(i) reval(xreplace(i, recast)))
    symbols <- lapply(symbols, function(i) xreplace(i, recast))
  }
  keep <- list()
  for (i in inequalities) {
    if (!is_relational(i) && i$type != "BooleanAtom") i <- Eq(i, S_Zero)
    if (.is_true(i)) next
    if (.is_false(i)) return(S_false)
    if (!length(free_symbols(i)) && .sp_is_number(.sub(i$args[[1L]], i$args[[2L]])))
      .sp_stop("NotImplementedError", paste0("could not determine truth value of ", .sp_py_str(i)))
    keep[[length(keep) + 1L]] <- i
  }
  rv <- .reduce_inequalities(keep, symbols)
  if (length(recast)) rv <- xreplace(rv, lapply(recast, function(p) list(p[[2L]], p[[1L]])))
  rv
}
