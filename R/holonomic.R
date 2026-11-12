# Holonomic functions, after sympy/holonomic (holonomic.py, recurrence.py,
# numerical.py). The base ring QQ[x] holds polynomials as expanded
# expressions; its fraction field holds cancelled rational functions. The
# linear systems are solved through their reduced row echelon form, which
# is unique, so particular solutions and null spaces match DomainMatrix's.

#' Holonomic functions
#'
#' `DifferentialOperators(base, generator)` returns the algebra and `Dx`;
#' the base is `holonomic_ring(x)` (SymPy's `QQ.old_poly_ring(x)`).
#' Operators combine with `+`, `-`, `*` and `^`. `HolonomicFunction(ann, x,
#' x0, y0)` is a holonomic function; `+`, `*`, `^`, `integrate_holonomic()`,
#' `diff_holonomic()`, `composition()`, `to_sequence()`, `series_holonomic()`,
#' `to_hyper()`, `to_expr_holonomic()`, `to_meijerg()`, `change_ics()`,
#' `shift_x()`, `change_x()`, `evalf_holonomic()` and `holonomic_indicial()`
#' are SymPy's methods. `expr_to_holonomic()`, `from_hyper()` and
#' `from_meijerg()` convert expressions. `RecurrenceOperators()` and
#' `HolonomicSequence()` are the recurrence counterparts.
#'
#' @param base A base ring from `holonomic_ring()`.
#' @param x,z,n Symbols.
#' @param dom The coefficient domain name.
#' @param generator The generator name.
#' @param annihilator A differential operator.
#' @param x0 The point of the initial conditions.
#' @param y0 Initial conditions: a list, or a named list keyed by exponents
#'   for singular conditions.
#' @param f A holonomic function.
#' @param recurrence A recurrence operator.
#' @param u0 Initial values.
#' @param expr,func An expression.
#' @param limits `x` or `list(x, a, b)`.
#' @param initcond Compute initial conditions.
#' @param lb Return the lower bound too.
#' @param coefficient Return coefficients.
#' @param order Add an Order term.
#' @param as_list Return a list of terms.
#' @param b,a Points.
#' @param lenics The number of initial conditions.
#' @param points Points for numerical evaluation.
#' @param method `"RK4"` or `"Euler"`.
#' @param step The step size.
#' @param derivatives Return derivatives too.
#' @param evalf Evaluate initial conditions numerically.
#' @param domain The coefficient domain.
#' @param ... Initial point and conditions for `composition()`.
#' @return An operator, holonomic function, sequence or expression.
#' @examples
#' x <- Symbol("x")
#' R <- DifferentialOperators(holonomic_ring(x), "Dx"); Dx <- R$Dx
#' p <- HolonomicFunction(Dx - 1, x, 0, list(1))
#' q <- HolonomicFunction(Dx^2 + 1, x, 0, list(0, 1))
#' p + q
#' p * q
#' to_expr_holonomic(q * q)
#' @name holonomic
NULL

#' @rdname holonomic
#' @export
holonomic_ring <- function(x, dom = "QQ") structure(list(x = as_sym(x), dom = dom), class = "holonomic_ring")
.hr_eq <- function(a, b) identical(a$x$key, b$x$key) && identical(a$dom, b$dom)
.hp <- function(e, R) expand(cancel(as_sym(e)))            # an element of the base ring
.hk <- function(e) cancel(as_sym(e))                         # an element of the fraction field
.hk_zero <- function(e) is_zero_num(cancel(e))
.hnum <- function(e) fraction(cancel(e))[[1L]]
.hden <- function(e) fraction(cancel(e))[[2L]]
.hmonic <- function(p, x) {
  p <- expand(p)
  if (is_zero_num(p)) return(p)
  expand(.div(p, .lc_in(p, x)))
}
.hgcd <- function(a, b, x) {
  if (is_zero_num(expand(a))) return(.hmonic(b, x))
  if (is_zero_num(expand(b))) return(.hmonic(a, x))
  .hmonic(gcd(a, b), x)
}
.hlcm <- function(a, b, x) .hmonic(lcm(a, b), x)
.hdeg <- function(p, x) { p <- expand(p); if (is_zero_num(p)) -Inf else .deg_in(p, x) }

# ---------------------------------------------------------------------------
# differential operators

#' @rdname holonomic
#' @export
DifferentialOperators <- function(base, generator = "Dx") {
  alg <- structure(list(base = base, gen_symbol = Symbol(generator, commutative = FALSE)), class = "DifferentialOperatorAlgebra")
  list(algebra = alg, Dx = .diffop(list(S_Zero, S_One), alg))
}

#' @export
print.DifferentialOperatorAlgebra <- function(x, ...) {
  cat("Univariate Differential Operator Algebra in intermediate ", .str(x$gen_symbol), " over the base ring ",
      x$base$dom, "[", .str(x$base$x), "]\n", sep = "")
  invisible(x)
}

.diffop <- function(polys, parent) {
  polys <- lapply(polys, function(p) .hp(p, parent$base))
  structure(list(polys = polys, parent = parent, order = length(polys) - 1L), class = "DifferentialOperator")
}
.hadd_lists <- function(a, b) {
  n <- max(length(a), length(b))
  lapply(seq_len(n), function(i) {
    if (i > length(a)) b[[i]] else if (i > length(b)) a[[i]] else .add(list(a[[i]], b[[i]]))
  })
}

.do_mul <- function(self, other) {
  x <- self$parent$base$x
  lo <- if (inherits(other, "DifferentialOperator")) other$polys else list(as_sym(other))
  ls <- self$polys
  sol <- lapply(lo, function(i) .mul(list(i, ls[[1L]])))
  for (i in seq_along(ls)[-1L]) {
    sol1 <- c(list(S_Zero), lo)
    sol2 <- lapply(lo, function(z) .diff1(z, x))
    lo <- .hadd_lists(sol1, sol2)
    sol <- .hadd_lists(sol, lapply(lo, function(z) .mul(list(z, ls[[i]]))))
  }
  .diffop(sol, self$parent)
}

.do_pow <- function(self, n) {
  if (n == 1) return(self)
  result <- .diffop(list(S_One), self$parent)
  if (n == 0) return(result)
  if (.do_equal(self, .diffop(list(S_Zero, S_One), self$parent)))
    return(.diffop(c(rep(list(S_Zero), n), list(S_One)), self$parent))
  xx <- self
  repeat {
    if (n %% 2) result <- .do_mul(result, xx)
    n <- n %/% 2
    if (!n) break
    xx <- .do_mul(xx, xx)
  }
  result
}

.do_keys <- function(d) vapply(d$polys, function(p) expand(p)$key, "")
.do_equal <- function(a, b) {
  if (inherits(b, "DifferentialOperator")) return(identical(.do_keys(a), .do_keys(b)))
  identical(expand(a$polys[[1L]])$key, as_sym(b)$key) && all(vapply(a$polys[-1L], is_zero_num, NA))
}

#' @export
chooseOpsMethod.DifferentialOperator <- function(x, y, mx, my, cl, reverse) TRUE
#' @export
Ops.DifferentialOperator <- function(e1, e2) {
  if (missing(e2)) {
    if (.Generic == "-") return(.do_mul(.diffop(list(S_NegativeOne), e1$parent), e1))
    return(e1)
  }
  op1 <- inherits(e1, "DifferentialOperator"); op2 <- inherits(e2, "DifferentialOperator")
  par <- if (op1) e1$parent else e2$parent
  as_op <- function(z) if (inherits(z, "DifferentialOperator")) z else .diffop(list(as_sym(z)), par)
  switch(.Generic,
    "*" = if (op1) .do_mul(e1, e2) else .diffop(lapply(e2$polys, function(p) .mul(list(as_sym(e1), p))), par),
    "+" = .diffop(.hadd_lists(as_op(e1)$polys, as_op(e2)$polys), par),
    "-" = .diffop(.hadd_lists(as_op(e1)$polys, .do_mul(.diffop(list(S_NegativeOne), par), as_op(e2))$polys), par),
    "/" = .do_mul(e1, .inv(as_sym(e2))),
    "^" = .do_pow(e1, as.integer(e2)),
    "==" = .do_equal(e1, e2),
    stop("unsupported operator ", .Generic, call. = FALSE))
}

.do_str <- function(d) {
  out <- ""
  g <- .str(d$parent$gen_symbol)
  for (i in seq_along(d$polys)) {
    j <- d$polys[[i]]
    if (is_zero_num(j)) next
    s <- .str(j)
    if (i == 1L) { out <- paste0(out, "(", s, ")"); next }
    if (nzchar(out)) out <- paste0(out, " + ")
    out <- paste0(out, "(", s, ")*", g, if (i > 2L) paste0("^", i - 1L) else "")
  }
  out
}
#' @export
print.DifferentialOperator <- function(x, ...) { cat(.do_str(x), "\n"); invisible(x) }
#' @export
format.DifferentialOperator <- function(x, ...) .do_str(x)

.do_is_singular <- function(d, x0) {
  rs <- roots(d$polys[[length(d$polys)]], d$parent$base$x)
  any(vapply(rs, function(r) identical(r$key, as_sym(x0)$key), NA))
}

# ---------------------------------------------------------------------------
# linear algebra over the fraction field: rref, particular solution and
# null space of A x = b, and _find_nonzero_solution

.hrref <- function(M) {
  nr <- length(M); nc <- if (nr) length(M[[1L]]) else 0L
  piv <- integer(); r <- 1L
  for (cidx in seq_len(nc)) {
    if (r > nr) break
    pr <- NA_integer_
    for (i in r:nr) if (!.hk_zero(M[[i]][[cidx]])) { pr <- i; break }
    if (is.na(pr)) next
    if (pr != r) { tmp <- M[[r]]; M[[r]] <- M[[pr]]; M[[pr]] <- tmp }
    pv <- M[[r]][[cidx]]
    M[[r]] <- lapply(M[[r]], function(z) .hk(.div(z, pv)))
    for (i in seq_len(nr)) if (i != r && !.hk_zero(M[[i]][[cidx]])) {
      f <- M[[i]][[cidx]]
      M[[i]] <- Map(function(a, b) .hk(.sub(a, .mul(list(f, b)))), M[[i]], M[[r]])
    }
    piv <- c(piv, cidx); r <- r + 1L
  }
  list(M = M, piv = piv)
}

# particular solution (free variables zero) and null space rows
.hsolve <- function(A, b) {
  nr <- length(A); nc <- length(A[[1L]])
  aug <- Map(function(row, bi) c(row, list(bi)), A, b)
  rr <- .hrref(aug)
  M <- rr$M; piv <- rr$piv[rr$piv <= nc]
  part <- rep(list(S_Zero), nc)
  for (k in seq_along(piv)) part[[piv[k]]] <- M[[k]][[nc + 1L]]
  free <- setdiff(seq_len(nc), piv)
  ns <- lapply(free, function(j) {
    v <- rep(list(S_Zero), nc); v[[j]] <- S_One
    for (k in seq_along(piv)) v[[piv[k]]] <- .hk(.neg(M[[k]][[j]]))
    v
  })
  list(particular = part, nullspace = ns)
}
.find_nonzero_solution <- function(A) {
  nr <- length(A)
  s <- .hsolve(A, rep(list(S_Zero), nr))
  sol <- s$particular
  for (v in s$nullspace) sol <- Map(function(a, b) .hk(.add(list(a, b))), sol, v)
  sol
}
.hzero_vec <- function(v) all(vapply(v, .hk_zero, NA))

# _normalize converts each coefficient with K.from_sympy: the numerator and the
# denominator go through dict_from_basic and QQ.from_sympy, which raises
# CoercionFailed on an irrational number (sympy/polys/domains/rationalfield.py:82).
# "QQ" here also stands for SymPy's QQ(a, ...) of expr_to_holonomic, whose
# FracField.from_expr refuses a coefficient that is not a rational function of
# a, ... (sympy/polys/fields.py:315); its symbols are taken in name order.
.hqq_check <- function(cs, x) {
  ratfun <- function(v) switch(v$type, Integer = , Rational = , Float = , Symbol = , Dummy = TRUE,
    Add = , Mul = all(vapply(v$args, ratfun, NA)), Pow = v$args[[2L]]$type == "Integer" && ratfun(v$args[[1L]]), FALSE)
  syms <- setdiff(sort(unique(unlist(lapply(cs, function(j) vapply(free_symbols(j), .str, ""))))), .str(x))
  for (j in cs) for (p in fraction(cancel(j))) {
    p <- expand(p)
    cf <- list()
    for (t in if (p$type == "Add") .sympy_args(p) else list(p)) {
      ce <- .so_as_coeff_exponent(t, x)
      k <- paste0("m", ce[[2L]]$key)
      cf[[k]] <- if (is.null(cf[[k]])) ce[[1L]] else .add(list(cf[[k]], ce[[1L]]))
    }
    for (v in cf) {
      if (!length(free_symbols(v))) {
        if (!(v$type %in% c("Integer", "Rational", "Float"))) stop(CoercionFailed(paste0("expected `Rational` object, got ", .str(v))))
      } else if (!ratfun(v)) stop(.sympy_condition("ValueError", paste0("expected an expression convertible to a rational function in Rational function field in ",
                                   paste(syms, collapse = ", "), " over QQ with lex order, got ", .str(v))))
    }
  }
}

# _normalize: clear the denominators and divide by the gcd of the numerators
.hnormalize <- function(lst, parent, negative = TRUE) {
  x <- parent$base$x
  cs <- lapply(lst, .hk)
  nums <- lapply(cs, .hnum); dens <- lapply(cs, .hden)
  L <- S_One
  for (d in dens) L <- .hlcm(d, L, x)
  if (negative) L <- .neg(L)
  cs <- lapply(cs, function(z) .hk(.mul(list(z, L))))
  G <- expand(cs[[length(cs)]])
  for (nm in nums) G <- .hgcd(nm, G, x)
  .diffop(lapply(cs, function(z) expand(cancel(.div(z, G)))), parent)
}

# ---------------------------------------------------------------------------
# holonomic functions

#' @rdname holonomic
#' @export
HolonomicFunction <- function(annihilator, x, x0 = 0, y0 = NULL)
  structure(list(annihilator = annihilator, x = as_sym(x), x0 = as_sym(x0),
                 y0 = if (is.null(y0)) NULL else .hy0(y0)), class = "HolonomicFunction")
.hy0 <- function(y0) {
  if (!is.list(y0)) y0 <- as.list(y0)
  if (!is.null(names(y0))) return(lapply(y0, function(v) lapply(as.list(v), as_sym)))
  lapply(y0, as_sym)
}
.is_singularics <- function(h) !is.null(h$y0) && !is.null(names(h$y0))
.have_ics <- function(h) length(h$y0) > 0L
.y0_str <- function(y0) {
  if (!is.null(names(y0))) return(paste0("{", paste(paste0(names(y0), ": [", vapply(y0, function(v) paste(vapply(v, .str, ""), collapse = ", "), ""), "]"), collapse = ", "), "}"))
  paste0("[", paste(vapply(y0, .str, ""), collapse = ", "), "]")
}
#' @export
format.HolonomicFunction <- function(x, ...) {
  if (.have_ics(x)) paste0("HolonomicFunction(", .do_str(x$annihilator), ", ", .str(x$x), ", ", .str(x$x0), ", ", .y0_str(x$y0), ")")
  else paste0("HolonomicFunction(", .do_str(x$annihilator), ", ", .str(x$x), ")")
}
#' @export
print.HolonomicFunction <- function(x, ...) { cat(format(x), "\n"); invisible(x) }

.hextend_y0 <- function(h, n) {
  if (.do_is_singular(h$annihilator, h$x0) || .is_singularics(h)) return(h$y0)
  ann <- h$annihilator; a <- ann$order; y0 <- h$y0
  x <- h$x
  if (length(y0) < a || n <= length(y0)) return(y0)
  lp <- lapply(ann$polys, .hk)
  red <- lapply(seq_len(a), function(i) .hk(.neg(.div(lp[[i]], lp[[a + 1L]]))))
  y1 <- y0[seq_len(min(length(y0), a))]
  for (k in seq_len(n - a)) {
    sol <- S_Zero
    for (i in seq_along(red)) {
      if (i > length(y1)) break
      r <- subs(red[[i]], x, h$x0)
      if (isFALSE(ask_fact(r, "finite")) || is_nan(r)) return(y0)
      sol <- .add(list(sol, .mul(list(y1[[i]], r))))
    }
    y1[[length(y1) + 1L]] <- sol
    red <- .hderivate(red, x)
  }
  c(y0, y1[seq_along(y1) > length(y0)])
}
.hderivate <- function(lp, x) {
  a <- length(lp)
  sol <- list(.hk(.diff1(lp[[1L]], x)))
  for (i in seq_along(lp)[-1L]) sol[[i]] <- .hk(.add(list(.diff1(lp[[i]], x), lp[[i - 1L]])))
  sol[[a + 1L]] <- lp[[a]]
  sol
}

.hrows <- function(ops, dim) lapply(ops, function(e) lapply(seq_len(dim + 1L), function(i) if (i <= length(e$polys)) .hk(e$polys[[i]]) else S_Zero))
.htranspose <- function(M) lapply(seq_along(M[[1L]]), function(j) lapply(M, `[[`, j))

.h_add <- function(self, other) {
  d1 <- self$annihilator$order; d2 <- other$annihilator$order
  dim <- max(d1, d2)
  gen <- self$annihilator$parent; Dx <- .diffop(list(S_Zero, S_One), gen)
  rs <- list(self$annihilator); ro <- list(other$annihilator)
  for (i in seq_len(dim - d1)) rs[[length(rs) + 1L]] <- .do_mul(Dx, rs[[length(rs)]])
  for (i in seq_len(dim - d2)) ro[[length(ro) + 1L]] <- .do_mul(Dx, ro[[length(ro)]])
  repeat {
    r <- .htranspose(.hrows(c(rs, ro), dim))
    sol <- .find_nonzero_solution(r)
    if (!.hzero_vec(sol)) break
    dim <- dim + 1L
    rs[[length(rs) + 1L]] <- .do_mul(Dx, rs[[length(rs)]])
    ro[[length(ro) + 1L]] <- .do_mul(Dx, ro[[length(ro)]])
  }
  sol <- sol[seq_len(dim + 1L - d1)]
  s1 <- .hnormalize(sol, gen)
  s2 <- .do_mul(s1, self$annihilator)
  sol <- .hnormalize(s2$polys, gen, negative = FALSE)
  x <- self$x
  if (!(.have_ics(self) && .have_ics(other))) return(HolonomicFunction(sol, x))
  if (!.is_singularics(self) && !.is_singularics(other)) {
    if (identical(self$x0$key, other$x0$key)) {
      y1 <- .hextend_y0(self, sol$order); y2 <- .hextend_y0(other, sol$order)
      k <- min(length(y1), length(y2))
      return(HolonomicFunction(sol, x, self$x0, Map(function(a, b) .add(list(a, b)), y1[seq_len(k)], y2[seq_len(k)])))
    }
    s0 <- .do_is_singular(self$annihilator, 0); o0 <- .do_is_singular(other$annihilator, 0)
    if (is_zero_num(self$x0) && !s0 && !o0) return(.h_add(self, change_ics(other, 0)))
    if (is_zero_num(other$x0) && !s0 && !o0) return(.h_add(change_ics(self, 0), other))
    sx <- .do_is_singular(self$annihilator, self$x0); ox <- .do_is_singular(other$annihilator, self$x0)
    if (!sx && !ox) return(.h_add(self, change_ics(other, self$x0)))
    return(.h_add(change_ics(self, other$x0), other))
  }
  if (!identical(self$x0$key, other$x0$key)) return(HolonomicFunction(sol, x))
  tosing <- function(h) list(`0` = lapply(seq_along(h$y0), function(i) .div(h$y0[[i]], fcall("factorial", list(.int(i - 1L))))))
  y1 <- if (.is_singularics(self)) self$y0 else tosing(self)
  y2 <- if (.is_singularics(other)) other$y0 else tosing(other)
  y0 <- list()
  for (i in names(y1)) y0[[i]] <- if (!is.null(y2[[i]])) Map(function(a, b) .add(list(a, b)), y1[[i]][seq_len(min(length(y1[[i]]), length(y2[[i]])))], y2[[i]]) else y1[[i]]
  for (i in names(y2)) if (is.null(y1[[i]])) y0[[i]] <- y2[[i]]
  HolonomicFunction(sol, x, self$x0, y0)
}

.h_mul <- function(self, other) {
  if (!inherits(other, "HolonomicFunction")) {
    other <- as_sym(other)
    if (.has1(other, self$x)) stop("Can't multiply a HolonomicFunction and expressions/functions.", call. = FALSE)
    if (!.have_ics(self)) return(self)
    y0 <- .hextend_y0(self, self$annihilator$order)
    return(HolonomicFunction(self$annihilator, self$x, self$x0, lapply(y0, function(j) expand(.mul(list(j, other))))))
  }
  as_ <- self$annihilator; ao <- other$annihilator
  a <- as_$order; b <- ao$order
  x <- self$x
  ls <- lapply(as_$polys, .hk); lo <- lapply(ao$polys, .hk)
  sred <- lapply(seq_len(a), function(i) .hk(.neg(.div(ls[[i]], ls[[a + 1L]]))))
  ored <- lapply(seq_len(b), function(i) .hk(.neg(.div(lo[[i]], lo[[b + 1L]]))))
  cm <- lapply(seq_len(a + 1L), function(i) rep(list(S_Zero), b + 1L))
  cm[[1L]][[1L]] <- S_One
  flat <- function() unlist(lapply(seq_len(a), function(i) lapply(seq_len(b), function(j) cm[[i]][[j]])), recursive = FALSE)
  els <- list(flat())
  sol <- .find_nonzero_solution(.htranspose(els))
  while (.hzero_vec(sol)) {
    for (i in rev(seq_len(a))) for (j in rev(seq_len(b))) {
      cm[[i]][[j + 1L]] <- .hk(.add(list(cm[[i]][[j + 1L]], cm[[i]][[j]])))
      cm[[i + 1L]][[j]] <- .hk(.add(list(cm[[i + 1L]][[j]], cm[[i]][[j]])))
      cm[[i]][[j]] <- .hk(.diff1(cm[[i]][[j]], x))
    }
    for (i in seq_len(a + 1L)) {
      if (.hk_zero(cm[[i]][[b + 1L]])) next
      for (j in seq_len(b)) cm[[i]][[j]] <- .hk(.add(list(cm[[i]][[j]], .mul(list(ored[[j]], cm[[i]][[b + 1L]])))))
      cm[[i]][[b + 1L]] <- S_Zero
    }
    for (j in seq_len(b)) {
      if (.hk_zero(cm[[a + 1L]][[j]])) next
      for (i in seq_len(a)) cm[[i]][[j]] <- .hk(.add(list(cm[[i]][[j]], .mul(list(sred[[i]], cm[[a + 1L]][[j]])))))
      cm[[a + 1L]][[j]] <- S_Zero
    }
    els[[length(els) + 1L]] <- flat()
    sol <- .find_nonzero_solution(.htranspose(els))
  }
  ann <- .hnormalize(sol, as_$parent, negative = FALSE)
  if (!(.have_ics(self) && .have_ics(other))) return(HolonomicFunction(ann, x))
  if (!.is_singularics(self) && !.is_singularics(other)) {
    if (identical(self$x0$key, other$x0$key)) {
      ys <- .hextend_y0(self, ann$order); yo <- .hextend_y0(other, ann$order)
      y0 <- list(.mul(list(ys[[1L]], yo[[1L]])))
      kk <- min(length(ys), length(yo))
      if (kk > 1L) for (i in 1:(kk - 1L)) {
        s <- S_Zero
        for (j in 0:i) s <- .add(list(s, .mul(list(sym_binomial(i, j), ys[[j + 1L]], yo[[i - j + 1L]]))))
        y0[[length(y0) + 1L]] <- s
      }
      return(HolonomicFunction(ann, x, self$x0, y0))
    }
    s0 <- .do_is_singular(as_, 0); o0 <- .do_is_singular(ao, 0)
    if (is_zero_num(self$x0) && !s0 && !o0) return(.h_mul(self, change_ics(other, 0)))
    if (is_zero_num(other$x0) && !s0 && !o0) return(.h_mul(change_ics(self, 0), other))
    sx <- .do_is_singular(as_, self$x0); ox <- .do_is_singular(ao, self$x0)
    if (!sx && !ox) return(.h_mul(self, change_ics(other, self$x0)))
    return(.h_mul(change_ics(self, other$x0), other))
  }
  if (!identical(self$x0$key, other$x0$key)) return(HolonomicFunction(ann, x))
  tosing <- function(h) list(`0` = lapply(seq_along(h$y0), function(i) .div(h$y0[[i]], fcall("factorial", list(.int(i - 1L))))))
  y1 <- if (.is_singularics(self)) self$y0 else tosing(self)
  y2 <- if (.is_singularics(other)) other$y0 else tosing(other)
  y0 <- list()
  for (i in names(y1)) for (j in names(y2)) {
    k <- min(length(y1[[i]]), length(y2[[j]]))
    cc <- lapply(seq_len(k) - 1L, function(aa) .add(lapply(0:aa, function(bb) .mul(list(y1[[i]][[bb + 1L]], y2[[j]][[aa - bb + 1L]])))))
    key <- .str(.add(list(sympify(i), sympify(j))))
    y0[[key]] <- if (is.null(y0[[key]])) cc else Map(function(u, v) .add(list(u, v)), cc, y0[[key]])
  }
  HolonomicFunction(ann, x, self$x0, y0)
}

.h_pow <- function(self, n) {
  ann <- self$annihilator
  if (ann$order <= 1) {
    y0 <- if (is.null(self$y0)) NULL else list(.pow(self$y0[[1L]], as_sym(n)))
    sol <- list(expand(.mul(list(ann$polys[[1L]], as_sym(n)))), ann$polys[[2L]])
    return(HolonomicFunction(.diffop(sol, ann$parent), self$x, self$x0, y0))
  }
  if (n < 0) stop("Negative Power on a Holonomic Function", call. = FALSE)
  Dx <- .diffop(list(S_Zero, S_One), ann$parent)
  result <- HolonomicFunction(Dx, self$x, 0, list(S_One))
  if (n == 0) return(result)
  xx <- self
  repeat {
    if (n %% 2) result <- .h_mul(result, xx)
    n <- n %/% 2
    if (!n) break
    xx <- .h_mul(xx, xx)
  }
  result
}

#' @export
chooseOpsMethod.HolonomicFunction <- function(x, y, mx, my, cl, reverse) TRUE
#' @export
Ops.HolonomicFunction <- function(e1, e2) {
  if (missing(e2)) return(if (.Generic == "-") .h_mul(e1, -1) else e1)
  h1 <- inherits(e1, "HolonomicFunction")
  switch(.Generic,
    "+" = .h_add(e1, e2),
    "-" = if (h1) .h_add(e1, .h_mul(e2, -1)) else .h_add(.h_mul(e2, -1), e1),
    "*" = if (h1) .h_mul(e1, e2) else .h_mul(e2, e1),
    "/" = .h_mul(e1, .inv(as_sym(e2))),
    "^" = .h_pow(e1, as.integer(e2)),
    "==" = .h_equal(e1, e2),
    stop("unsupported operator ", .Generic, call. = FALSE))
}
.h_equal <- function(a, b) {
  if (!.do_equal(a$annihilator, b$annihilator) || !identical(a$x$key, b$x$key)) return(FALSE)
  if (.have_ics(a) && .have_ics(b)) return(identical(a$x0$key, b$x0$key) && identical(.y0_str(a$y0), .y0_str(b$y0)))
  TRUE
}

#' @rdname holonomic
#' @export
integrate_holonomic <- function(f, limits, initcond = FALSE) {
  D <- .diffop(list(S_Zero, S_One), f$annihilator$parent)
  annD <- .do_mul(f$annihilator, D)
  x <- f$x
  if (.is_singularics(f)) {
    r <- .singularics_to_ord(f)
    if (!is.null(r)) return(integrate_holonomic(r, limits, initcond))
    y0 <- list()
    for (i in names(f$y0)) {
      cc <- f$y0[[i]]; iv <- sympify(i)
      y0[[.str(.add(list(iv, S_One)))]] <- lapply(seq_along(cc), function(j) {
        if (is_zero_num(cc[[j]])) return(S_Zero)
        den <- .add(list(iv, .int(j)))
        if (is_zero_num(den)) stop("logarithmic terms in the series are not supported", call. = FALSE)
        .div(cc[[j]], den)
      })
    }
    if (is.list(limits)) stop("Definite integration for singular initial conditions", call. = FALSE)
    return(HolonomicFunction(annD, x, f$x0, y0))
  }
  if (!.have_ics(f)) return(if (initcond) HolonomicFunction(annD, x, f$x0, list(S_Zero)) else HolonomicFunction(annD, x))
  definite <- FALSE
  if (is.list(limits) && length(limits) == 3L && identical(as_sym(limits[[1L]])$key, x$key)) {
    x0 <- f$x0; a <- as_sym(limits[[2L]]); b <- as_sym(limits[[3L]]); definite <- TRUE
  }
  y0 <- c(list(S_Zero), f$y0)
  ind <- HolonomicFunction(annD, x, f$x0, y0)
  if (!definite) return(ind)
  if (!identical(x0$key, a$key)) {
    ie <- tryCatch(to_expr_holonomic(ind), error = function(e) NULL)
    lower <- if (!is.null(ie)) { v <- subs(ie, x, a); if (is_nan(v)) limit(ie, x, a) else v } else evalf_holonomic(ind, a)
    if (identical(b$key, x$key)) { y0[[1L]] <- .sub(y0[[1L]], lower); return(HolonomicFunction(annD, x, x0, y0)) }
    if (is_number_type(b)) {
      upper <- if (!is.null(ie)) { v <- subs(ie, x, b); if (is_nan(v)) limit(ie, x, b) else v } else evalf_holonomic(ind, b)
      return(.sub(upper, lower))
    }
  }
  if (identical(b$key, x$key)) return(HolonomicFunction(annD, x, a, y0))
  if (is_number_type(b)) {
    s <- tryCatch(to_expr_holonomic(HolonomicFunction(annD, x, a, y0)), error = function(e) NULL)
    if (is.null(s)) return(evalf_holonomic(HolonomicFunction(annD, x, a, y0), b))
    v <- subs(s, x, b)
    return(if (!is_nan(v)) v else limit(s, x, b))
  }
  HolonomicFunction(annD, x)
}

.singularics_to_ord <- function(h) {
  a <- names(h$y0)[1L]; b <- h$y0[[1L]]; av <- suppressWarnings(as.numeric(a))
  if (length(h$y0) == 1L && !is.na(av) && av == floor(av) && av > 0) {
    y0 <- c(rep(list(S_Zero), av), lapply(seq_along(b), function(i) .mul(list(b[[i]], fcall("factorial", list(.int(av + i - 1L)))))))
    return(HolonomicFunction(h$annihilator, h$x, h$x0, y0))
  }
  NULL
}

#' @rdname holonomic
#' @export
diff_holonomic <- function(f, x = NULL, n = 1L) {
  if (!is.null(x) && !identical(as_sym(x)$key, f$x$key)) return(S_Zero)
  if (n > 1L) { for (k in seq_len(n)) f <- diff_holonomic(f); return(f) }
  ann <- f$annihilator
  if (is_zero_num(ann$polys[[1L]]) && ann$order == 1L) return(S_Zero)
  if (is_zero_num(ann$polys[[1L]])) {
    sol <- .diffop(ann$polys[-1L], ann$parent)
    if (.have_ics(f) && !.is_singularics(f)) return(HolonomicFunction(sol, f$x, f$x0, f$y0[-1L]))
    return(HolonomicFunction(sol, f$x))
  }
  x <- f$x
  seq <- lapply(ann$polys, .hk)
  rhs <- c(list(S_Zero), lapply(seq[-1L], function(i) .hk(.div(i, seq[[1L]]))))
  sol <- .hderivate(rhs, x)
  sol <- .hadd_lists(sol, list(S_Zero, S_One))
  sol <- .hnormalize(sol[-1L], ann$parent, negative = FALSE)
  if (!.have_ics(f) || .is_singularics(f)) return(HolonomicFunction(sol, x))
  y0 <- .hextend_y0(f, sol$order + 1L)[-1L]
  HolonomicFunction(sol, x, f$x0, y0)
}

#' @rdname holonomic
#' @export
composition <- function(f, expr, ...) {
  extra <- list(...)
  x <- f$x; expr <- as_sym(expr)
  a <- f$annihilator$order
  dd <- .diff1(expr, x)
  lp <- f$annihilator$polys
  r <- subs(lp[[a + 1L]], x, expr)
  sbs <- lapply(seq_len(a), function(i) .neg(.div(subs(lp[[i]], x, expr), r)))
  coeffs <- c(list(S_One), rep(list(S_Zero), a - 1L))
  system <- list(coeffs)
  repeat {
    cn <- lapply(coeffs, function(p) .diff1(p, x))
    if (a > 1L) for (i in seq_len(a - 1L)) cn[[i + 1L]] <- .add(list(cn[[i + 1L]], .mul(list(coeffs[[i]], dd))))
    for (i in seq_len(a)) cn[[i]] <- .add(list(cn[[i]], .mul(list(coeffs[[a]], sbs[[i]], dd))))
    coeffs <- cn
    system[[length(system) + 1L]] <- coeffs
    # gauss_jordan_solve of the transposed system against zero, the free
    # parameter set to 1
    M <- lapply(seq_len(a), function(i) lapply(system, function(cv) .hk(cv[[i]])))
    s <- .hsolve(M, rep(list(S_Zero), a))
    if (length(s$nullspace)) break
  }
  sol <- s$particular
  # gauss_jordan_solve parametrizes with the first free variable; tau = 1
  for (v in s$nullspace[1L]) sol <- Map(function(p, q) .hk(.add(list(p, q))), sol, v)
  if (identical(f$annihilator$parent$base$dom, "QQ")) .hqq_check(sol, x)
  sol <- .hnormalize(sol, f$annihilator$parent, negative = FALSE)
  if (length(extra)) return(HolonomicFunction(sol, x, extra[[1L]], extra[[2L]]))
  HolonomicFunction(sol, x)
}

#' @rdname holonomic
#' @export
holonomic_degree <- function(f) max(vapply(f$annihilator$polys, function(p) .hdeg(p, f$x), 0))

# ---------------------------------------------------------------------------
# recurrences

#' @rdname holonomic
#' @export
RecurrenceOperators <- function(base, generator = "Sn") {
  alg <- structure(list(base = base, gen_symbol = Symbol(generator, commutative = FALSE)), class = "RecurrenceOperatorAlgebra")
  list(algebra = alg, Sn = .recop(list(S_Zero, S_One), alg))
}
.recop <- function(polys, parent) structure(list(polys = lapply(polys, function(p) expand(as_sym(p))), parent = parent,
                                                 order = length(polys) - 1L), class = "RecurrenceOperator")
.ro_str <- function(d) {
  out <- ""
  for (i in seq_along(d$polys)) {
    j <- d$polys[[i]]
    if (is_zero_num(j)) next
    s <- .str(j)
    if (i == 1L) { out <- paste0(out, "(", s, ")"); next }
    if (nzchar(out)) out <- paste0(out, " + ")
    out <- paste0(out, "(", s, ")Sn", if (i > 2L) paste0("^", i - 1L) else "")
  }
  out
}
#' @export
print.RecurrenceOperator <- function(x, ...) { cat(.ro_str(x), "\n"); invisible(x) }
.ro_mul <- function(self, other) {
  n <- self$parent$base$x
  lo <- if (inherits(other, "RecurrenceOperator")) other$polys else list(as_sym(other))
  ls <- self$polys
  sol <- lapply(lo, function(i) .mul(list(i, ls[[1L]])))
  for (i in seq_along(ls)[-1L]) {
    lo <- c(list(S_Zero), lapply(lo, function(z) expand(subs(z, n, .add(list(n, S_One))))))
    sol <- .hadd_lists(sol, lapply(lo, function(z) .mul(list(z, ls[[i]]))))
  }
  .recop(sol, self$parent)
}
#' @export
chooseOpsMethod.RecurrenceOperator <- function(x, y, mx, my, cl, reverse) TRUE
#' @export
Ops.RecurrenceOperator <- function(e1, e2) {
  r1 <- inherits(e1, "RecurrenceOperator")
  par <- if (r1) e1$parent else e2$parent
  as_op <- function(z) if (inherits(z, "RecurrenceOperator")) z else .recop(list(as_sym(z)), par)
  switch(.Generic,
    "*" = if (r1) .ro_mul(e1, e2) else .recop(lapply(e2$polys, function(p) .mul(list(as_sym(e1), p))), par),
    "+" = .recop(.hadd_lists(as_op(e1)$polys, as_op(e2)$polys), par),
    "-" = .recop(.hadd_lists(as_op(e1)$polys, .ro_mul(.recop(list(S_NegativeOne), par), as_op(e2))$polys), par),
    "^" = { n <- as.integer(e2); res <- .recop(list(S_One), par); for (k in seq_len(n)) res <- .ro_mul(res, e1); res },
    stop("unsupported operator ", .Generic, call. = FALSE))
}

#' @rdname holonomic
#' @export
HolonomicSequence <- function(recurrence, u0 = list()) {
  if (!is.list(u0)) u0 <- list(u0)
  structure(list(recurrence = recurrence, u0 = lapply(u0, as_sym), n = recurrence$parent$base$x), class = "HolonomicSequence")
}
#' @export
format.HolonomicSequence <- function(x, ...) {
  s <- paste0("HolonomicSequence(", .ro_str(x$recurrence), ", ", .str(x$n), ")")
  if (!length(x$u0)) return(s)
  paste0(s, paste0(", u(", seq_along(x$u0) - 1L, ") = ", vapply(x$u0, .str, ""), collapse = ""))
}
#' @export
print.HolonomicSequence <- function(x, ...) { cat(format(x), "\n"); invisible(x) }

.zroots <- function(p, n) {
  rs <- Filter(function(r) r$type == "Integer", roots(p, n, filter = "Z"))
  if (!length(rs)) return(integer())
  vapply(rs, function(r) as.numeric(as.character(r$value)), 0)
}

#' @rdname holonomic
#' @export
to_sequence <- function(f, lb = TRUE) {
  if (!is_zero_num(f$x0)) return(to_sequence(shift_x(f, f$x0)))
  if (.do_is_singular(f$annihilator, f$x0)) return(.frobenius(f, lb))
  n <- Symbol("n", integer = TRUE)
  R <- RecurrenceOperators(holonomic_ring(n, f$annihilator$parent$base$dom), "Sn")$algebra
  x <- f$x
  d1 <- list(); dkeys <- list()
  addd <- function(key, v) {
    k <- paste(key, collapse = ",")
    if (is.null(d1[[k]])) { d1[[k]] <<- v; dkeys[[k]] <<- key } else d1[[k]] <<- .add(list(d1[[k]], v))
  }
  for (i in seq_along(f$annihilator$polys)) {
    cl <- rev(.coeff_list(f$annihilator$polys[[i]], x))
    for (k in seq_along(cl)) {
      cf <- cl[[k]]
      if (is_zero_num(cf)) next
      ii <- i - 1L; kk <- k - 1L
      addd(c(ii - kk, kk), .mul(list(cf, fcall("RisingFactorial", list(.add(list(n, .int(-kk + 1L))), .int(ii))))))
    }
  }
  keylist <- vapply(dkeys, `[[`, 0, 1L)
  lower <- min(keylist); upper <- max(keylist)
  degree <- holonomic_degree(f)
  smallest_n <- lower + degree
  sol <- lapply(lower:upper, function(j) {
    if (j %in% keylist) .add(lapply(names(d1)[keylist == j], function(k) subs(d1[[k]], n, .sub(n, .int(lower))))) else S_Zero
  })
  sol <- .recop(sol, R)
  order <- sol$order
  ar <- .zroots(sol$polys[[length(sol$polys)]], n)
  if (length(ar)) smallest_n <- max(max(ar) + 1, smallest_n)
  order <- order + smallest_n
  y0 <- .hextend_y0(f, order)
  u0 <- lapply(seq_along(y0), function(i) .div(y0[[i]], fcall("factorial", list(.int(i - 1L)))))
  if (length(u0) < order) {
    dummys <- list(); unknowns <- list(); eqs <- list()
    for (i in seq_len(degree) - 1L) {
      eq <- S_Zero
      for (k in names(d1)) {
        jj <- dkeys[[k]]; idx <- i + jj[1L]
        key <- as.character(idx)
        if (idx < 0) dummys[[key]] <- S_Zero
        else if (idx < length(u0)) dummys[[key]] <- u0[[idx + 1L]]
        else if (is.null(dummys[[key]])) { dummys[[key]] <- Symbol(paste0("C_", idx)); unknowns[[length(unknowns) + 1L]] <- dummys[[key]] }
        if (jj[2L] <= i) eq <- .add(list(eq, .mul(list(subs(d1[[k]], n, .int(i)), dummys[[key]]))))
      }
      eqs[[length(eqs) + 1L]] <- eq
    }
    soleqs <- if (length(unknowns)) solve(eqs, unknowns, dict = TRUE) else list()
    sd <- if (length(soleqs)) soleqs[[1L]] else list()
    for (i in seq_len(order - length(u0)) + length(u0) - 1L) {
      key <- as.character(i)
      if (is.null(dummys[[key]])) dummys[[key]] <- Symbol(paste0("C_", i))
      dn <- dummys[[key]]$name
      u0[[length(u0) + 1L]] <- if (!is.null(dummys[[key]]$name) && !is.null(sd[[dn]])) sd[[dn]] else dummys[[key]]
    }
  }
  if (lb) return(list(list(HolonomicSequence(sol, u0), smallest_n)))
  list(HolonomicSequence(sol, u0))
}

# HolonomicFunction._frobenius: series solutions at a regular singular point
.frobenius <- function(f, lb = TRUE) {
  ir <- holonomic_indicial(f)
  mult <- attr(ir, "multiplicity")
  ordk <- unlist(.key_sort(lapply(ir, .sort_key), seq_along(ir)))
  reals <- list(); compl <- list()
  for (k in ordk) {
    i <- ir[[k]]
    if (isTRUE(ask_fact(i, "real"))) reals <- c(reals, rep(list(i), mult[k]))
    else { ri <- .re_im(i); compl <- c(compl, rep(list(list(i, ri[[1L]], ri[[2L]])), mult[k])) }
  }
  if (length(compl)) {
    o <- order(vapply(compl, function(t) num_value_dbl(t[[2L]]), 0)); compl <- compl[o]
    o <- order(vapply(compl, function(t) num_value_dbl(t[[3L]]), 0)); compl <- compl[o]
  }
  if (length(reals)) reals <- reals[order(vapply(reals, function(r) num_value_dbl(r), 0))]
  isint <- function(v) isTRUE(ask_fact(v, "integer")) || (is_number_type(v) && v$type != "Float" && gmp::denominator(gmp::as.bigq(v$value)) == 1)
  grp <- list()
  for (i in reals) {
    placed <- FALSE
    for (g in seq_along(grp)) if (isint(.sub(grp[[g]][[1L]], i))) { grp[[g]] <- c(grp[[g]], list(i)); placed <- TRUE; break }
    if (!placed) grp[[length(grp) + 1L]] <- list(i)
  }
  independent <- all(lengths(grp) == 1L)
  allpos <- all(vapply(reals, function(i) num_value_dbl(i) >= 0, NA))
  allint <- all(vapply(reals, isint, NA))
  if (.is_singularics(f)) {
    rts <- list()
    for (i in names(f$y0)[order(vapply(names(f$y0), function(k) num_value_dbl(sympify(k)), 0))])
      for (k in ordk) if (isTRUE(.is_true(Eq(ir[[k]], sympify(i))))) rts[[length(rts) + 1L]] <- sympify(i)
  } else if (allpos && allint) {
    rts <- list(reals[[1L]])
  } else if (independent) {
    rts <- c(lapply(grp, `[[`, 1L), lapply(compl, `[[`, 1L))
  } else if (!allint) {
    rts <- Filter(function(i) !isint(i), reals)
  } else {
    y00 <- if (.have_ics(f)) f$y0[[1L]] else NULL
    if (!.have_ics(f) || isFALSE(ask_fact(if (is.list(y00)) y00[[1L]] else y00, "finite"))) rts <- list(reals[[1L]])
    else rts <- list(Filter(function(i) num_value_dbl(i) >= 0, reals)[[1L]])
  }
  n <- Symbol("n", integer = TRUE)
  R <- RecurrenceOperators(holonomic_ring(n, f$annihilator$parent$base$dom), "Sn")$algebra
  x <- f$x
  finalsol <- list()
  ch <- utf8ToInt("C")
  for (p in rts) {
    d1 <- list(); dkeys <- list()
    addd <- function(key, v) {
      k <- paste(key, collapse = ",")
      if (is.null(d1[[k]])) { d1[[k]] <<- v; dkeys[[k]] <<- key } else d1[[k]] <<- .add(list(d1[[k]], v))
    }
    for (i in seq_along(f$annihilator$polys)) {
      cl <- rev(.coeff_list(f$annihilator$polys[[i]], x))
      for (k in seq_along(cl)) {
        cf <- cl[[k]]
        if (is_zero_num(cf)) next
        ii <- i - 1L; kk <- k - 1L
        addd(c(ii - kk, kk - ii), .mul(list(cf, fcall("RisingFactorial", list(.add(list(n, .int(-kk + 1L), p)), .int(ii))))))
      }
    }
    keylist <- vapply(dkeys, `[[`, 0, 1L)
    lower <- min(keylist); upper <- max(keylist)
    degree <- max(vapply(dkeys, `[[`, 0, 2L)); degree2 <- min(vapply(dkeys, `[[`, 0, 2L))
    smallest_n <- lower + degree
    sol <- lapply(lower:upper, function(j) {
      if (j %in% keylist) .add(lapply(names(d1)[keylist == j], function(k) subs(d1[[k]], n, .sub(n, .int(lower))))) else S_Zero
    })
    sol <- .recop(sol, R)
    order <- sol$order
    ar <- .zroots(sol$polys[[length(sol$polys)]], n)
    if (length(ar)) smallest_n <- max(max(ar) + 1, smallest_n)
    order <- order + smallest_n
    u0 <- list()
    pint <- isint(p) && num_value_dbl(p) >= 0
    if (.is_singularics(f)) u0 <- f$y0[[.str(p)]]
    else if (pint && length(rts) == 1L) {
      pv <- .int_val(p)
      y0 <- .hextend_y0(f, order + pv)
      if (length(y0) > pv) u0 <- lapply((pv + 1L):length(y0), function(i) .div(y0[[i]], fcall("factorial", list(.int(i - 1L)))))
    }
    if (length(u0) < order) {
      dummys <- list(); unknowns <- list(); eqs <- list()
      if (degree > degree2) for (i in degree2:(degree - 1)) {
        eq <- S_Zero
        for (k in names(d1)) {
          jj <- dkeys[[k]]; idx <- i + jj[1L]; key <- as.character(idx)
          if (idx < 0) dummys[[key]] <- S_Zero
          else if (idx < length(u0)) dummys[[key]] <- u0[[idx + 1L]]
          else if (is.null(dummys[[key]])) { dummys[[key]] <- Symbol(paste0(intToUtf8(ch), "_", idx)); unknowns[[length(unknowns) + 1L]] <- dummys[[key]] }
          if (jj[2L] <= i) eq <- .add(list(eq, .mul(list(subs(d1[[k]], n, .int(i)), dummys[[key]]))))
        }
        eqs[[length(eqs) + 1L]] <- eq
      }
      soleqs <- if (length(unknowns)) solve(eqs, unknowns, dict = TRUE) else list()
      sd <- if (length(soleqs)) soleqs[[1L]] else list()
      for (i in seq_len(order - length(u0)) + length(u0) - 1L) {
        key <- as.character(i)
        if (is.null(dummys[[key]])) dummys[[key]] <- Symbol(paste0(intToUtf8(ch), "_", i))
        dn <- dummys[[key]]$name
        u0[[length(u0) + 1L]] <- if (!is.null(dn) && !is.null(sd[[dn]])) sd[[dn]] else dummys[[key]]
      }
    }
    finalsol[[length(finalsol) + 1L]] <- if (lb) list(HolonomicSequence(sol, u0), p, smallest_n) else list(HolonomicSequence(sol, u0), p)
    ch <- ch + 1L
  }
  finalsol
}

#' @rdname holonomic
#' @export
shift_x <- function(f, a) {
  x <- f$x
  sol <- .diffop(lapply(f$annihilator$polys, function(p) subs(p, x, .add(list(x, as_sym(a))))), f$annihilator$parent)
  x0 <- .sub(f$x0, as_sym(a))
  if (!.have_ics(f)) return(HolonomicFunction(sol, x))
  HolonomicFunction(sol, x, x0, f$y0)
}

#' @rdname holonomic
#' @export
change_x <- function(f, z) {
  z <- as_sym(z)
  R <- holonomic_ring(z, f$annihilator$parent$base$dom)
  par <- DifferentialOperators(R, "Dx")$algebra
  sol <- .diffop(lapply(f$annihilator$polys, function(p) subs(p, f$x, z)), par)
  HolonomicFunction(sol, z, f$x0, f$y0)
}

.unpack_recur <- function(rec) {
  if (length(rec) == 1L && is.list(rec[[1L]]) && !inherits(rec[[1L]], "HolonomicSequence")) rec <- rec[[1L]]
  if (length(rec) == 2L) return(list(seq = rec[[1L]], smallest = rec[[2L]], cp = 0))
  list(seq = rec[[1L]], cp = rec[[2L]], smallest = rec[[3L]])
}

#' @rdname holonomic
#' @export
series_holonomic <- function(f, n = 6, coefficient = FALSE, order = TRUE) {
  ur <- .unpack_recur(to_sequence(f))
  rec <- ur$seq; cp <- ur$cp
  n <- n - as.integer(cp)
  l <- length(rec$u0) - 1L; k <- rec$recurrence$order
  x <- f$x; nn <- rec$n
  seq <- lapply(rec$recurrence$polys, .hk)
  sub <- lapply(seq_len(k), function(i) .hk(.neg(.div(seq[[i]], seq[[k + 1L]]))))
  sol <- rec$u0
  if (l + 1L < n) for (i in (l + 1L - k):(n - k - 1L)) {
    terms <- lapply(seq_len(k) - 1L, function(j) if (i + j >= 0) .mul(list(subs(sub[[j + 1L]], nn, .int(i)), sol[[i + j + 1L]])) else S_Zero)
    sol[[length(sol) + 1L]] <- .add(terms)
  }
  if (coefficient) return(sol)
  ser <- .add(lapply(seq_along(sol), function(i) .mul(list(.pow(x, .add(list(.int(i - 1L), as_sym(cp)))), sol[[i]]))))
  if (order) ser <- .add(list(ser, O(.pow(x, .int(n + as.integer(cp))), x)))
  if (!is_zero_num(f$x0)) ser <- subs(ser, x, .sub(x, f$x0))
  ser
}

#' @rdname holonomic
#' @export
holonomic_indicial <- function(f) {
  if (!is_zero_num(f$x0)) return(holonomic_indicial(shift_x(f, f$x0)))
  lc <- f$annihilator$polys; x <- f$x
  pole_degree <- function(p) {
    rs <- roots(p, x, filter = "Z"); m <- attr(rs, "multiplicity")
    w <- which(vapply(rs, is_zero_num, NA))
    if (length(w)) m[w[1L]] else 0
  }
  degree <- max(vapply(lc, function(p) .hdeg(p, x), 0))
  inf <- 10 * (max(1, degree) + max(1, f$annihilator$order))
  deg <- function(q) if (is_zero_num(q)) inf else pole_degree(q)
  b <- min(vapply(seq_along(lc), function(j) deg(lc[[j]]) - (j - 1L), 0))
  s <- S_Zero; y <- S_One
  for (i in seq_along(lc)) {
    cl <- .coeff_list(lc[[i]], x); dg <- length(cl) - 1L
    ii <- i - 1L
    if (0 <= ii + b && ii + b <= dg) s <- .add(list(s, .mul(list(cl[[dg - ii - b + 1L]], y))))
    y <- expand(.mul(list(y, .sub(x, .int(ii)))))
  }
  roots(s, x)
}

#' @rdname holonomic
#' @export
to_hyper <- function(f, as_list = FALSE) {
  ur <- .unpack_recur(to_sequence(f))
  rec <- ur$seq; smallest_n <- ur$smallest; cp <- as_sym(ur$cp)
  u0 <- rec$u0; r <- rec$recurrence; x <- f$x; x0 <- f$x0; nn <- rec$n
  m <- r$order
  if (m == 0L) {
    nz <- roots(r$polys[[1L]], nn, filter = "R")
    sol <- S_Zero
    for (j in seq_along(nz)) {
      i <- nz[[j]]
      if (!isTRUE(ask_fact(i, "integer")) || num_sign(i) < 0) next
      iv <- .int_val(i)
      sol <- .add(list(sol, .mul(list(if (iv < length(u0)) u0[[iv + 1L]] else Symbol(paste0("C_", j - 1L)), .pow(x, i)))))
    }
    sol <- .mul(list(sol, .pow(x, cp)))
    if (as_list) return(list(list(if (!is_zero_num(x0)) subs(sol, x, .sub(x, x0)) else sol)))
    return(if (!is_zero_num(x0)) subs(sol, x, .sub(x, x0)) else sol)
  }
  if (smallest_n + m > length(u0)) stop(.sympy_condition("NotImplementedError", "Can't compute sufficient Initial Conditions"))
  # NotHyperSeriesError(self, self.x0) (sympy/holonomic/holonomic.py:1904)
  if (m > 1L && any(!vapply(r$polys[2:m], is_zero_num, NA)))
    stop(.sympy_condition("NotHyperSeriesError", paste0("Power series expansion of ", format(f), " about ", .str(x0), " is not hypergeometric")))
  a <- r$polys[[1L]]; b <- r$polys[[length(r$polys)]]
  da <- .hdeg(a, nn); db <- .hdeg(b, nn)
  cc <- .neg(.div(.mul(list(.lc_in(a, nn), .pow(.int(m), .int(da)))), .mul(list(.lc_in(b, nn), .pow(.int(m), .int(db))))))
  arg1 <- roots(a, nn); arg2 <- roots(b, nn)
  m1 <- attr(arg1, "multiplicity"); m2 <- attr(arg2, "multiplicity")
  ord1 <- .key_sort(lapply(arg1, .sort_key), seq_along(arg1)); ord2 <- .key_sort(lapply(arg2, .sort_key), seq_along(arg2))
  sol <- S_Zero; lst <- list()
  for (i in seq_len(smallest_n + m) - 1L) {
    if (i < smallest_n) {
      t <- .mul(list(u0[[i + 1L]], .pow(x, .int(i))))
      if (as_list) lst[[length(lst) + 1L]] <- list(subs(.mul(list(u0[[i + 1L]], .pow(x, .add(list(.int(i), cp))))), x, .sub(x, x0)))
      else sol <- .add(list(sol, t))
      next
    }
    if (is_zero_num(u0[[i + 1L]])) next
    ap <- list(); bq <- list()
    for (k in unlist(ord1)) ap <- c(ap, rep(list(nsimplify(.div(.sub(.int(i), arg1[[k]]), .int(m)))), m1[k]))
    for (k in unlist(ord2)) bq <- c(bq, rep(list(nsimplify(.div(.sub(.int(i), arg2[[k]]), .int(m)))), m2[k]))
    w <- which(vapply(bq, is_one_num, NA))
    if (length(w)) bq <- bq[-w[1L]] else ap <- c(ap, list(S_One))
    h <- hyper(ap, bq, .mul(list(cc, .pow(x, .int(m)))))
    if (as_list) lst[[length(lst) + 1L]] <- list(subs(.mul(list(u0[[i + 1L]], .pow(x, .add(list(.int(i), cp))))), x, .sub(x, x0)), subs(h, x, .sub(x, x0)))
    else sol <- .add(list(sol, .mul(list(u0[[i + 1L]], h, .pow(x, .int(i))))))
  }
  if (as_list) return(lst)
  sol <- .mul(list(sol, .pow(x, cp)))
  if (!is_zero_num(x0)) return(subs(sol, x, .sub(x, x0)))
  sol
}

#' @rdname holonomic
#' @export
to_expr_holonomic <- function(f) simplify(hyperexpand(to_hyper(f)))

#' @rdname holonomic
#' @export
change_ics <- function(f, b, lenics = NULL) {
  if (is.null(lenics) && length(f$y0) > f$annihilator$order) lenics <- length(f$y0)
  sol <- tryCatch(expr_to_holonomic(to_expr_holonomic(f), x = f$x, x0 = b, lenics = lenics), error = function(e) NULL)
  if (!is.null(sol) && identical(sol$x0$key, as_sym(b)$key)) return(sol)
  y0 <- evalf_holonomic(f, b, derivatives = TRUE)
  HolonomicFunction(f$annihilator, f$x, b, y0)
}

.hyper_to_meijerg <- function(func) {
  ap <- .hy_ap(func); bq <- .hy_bq(func)
  if (any(vapply(ap, function(i) isTRUE(ask_fact(i, "integer")) && !isTRUE(ask_fact(i, "positive")), NA))) return(hyperexpand(func))
  z <- func$args[[3L]]
  k <- .div(.mul(lapply(bq, function(i) fcall("gamma", list(i)))), .mul(lapply(ap, function(i) fcall("gamma", list(i)))))
  .mul(list(k, meijerg(lapply(ap, function(i) .sub(S_One, i)), list(), list(S_Zero), lapply(bq, function(i) .sub(S_One, i)), .neg(z))))
}

#' @rdname holonomic
#' @export
to_meijerg <- function(f) {
  rep <- to_hyper(f, as_list = TRUE)
  .add(lapply(rep, function(i) if (length(i) == 1L) i[[1L]] else .mul(list(i[[1L]], .hyper_to_meijerg(i[[2L]])))))
}

# ---------------------------------------------------------------------------
# numerical evaluation (numerical.py)

#' @rdname holonomic
#' @export
evalf_holonomic <- function(f, points, method = "RK4", step = 0.05, derivatives = FALSE) {
  lp <- FALSE
  h <- step
  if (!is.list(points) && length(points) == 1L) {
    lp <- TRUE
    b <- as_sym(points)
    if (identical(f$x0$key, b$key)) return(.h_evalf(f, list(b), method, derivatives)[[1L]])
    a <- num_value_dbl(f$x0); bv <- num_value_dbl(b)
    if (a > bv) h <- -h
    nst <- as.integer((bv - a) / h)
    pts <- list(a + h)
    for (i in seq_len(nst - 1L)) pts[[length(pts) + 1L]] <- pts[[length(pts)]] + h
    points <- pts
  }
  for (r in roots(f$annihilator$polys[[length(f$annihilator$polys)]], f$x)) {
    rv <- tryCatch(num_value_dbl(r), error = function(e) NA)
    if (identical(r$key, f$x0$key) || (!is.na(rv) && any(vapply(points, function(p) isTRUE(all.equal(as.numeric(p), rv)), NA))))
      stop("Holonomic function has a singularity at ", .str(r), call. = FALSE)
  }
  res <- .h_evalf(f, points, method, derivatives)
  if (lp) res[[length(res)]] else res
}

.h_evalf <- function(f, points, method, derivatives) {
  ann <- f$annihilator; a <- ann$order; x <- f$x
  lp <- lapply(ann$polys, function(p) p)
  red <- lapply(seq_len(a), function(i) .div(.neg(lp[[i]]), lp[[a + 1L]]))
  y0 <- vapply(f$y0, function(v) num_value_dbl(N(v)), 0)
  if (length(y0) < a) stop("Not Enough Initial Conditions", call. = FALSE)
  fx <- lapply(red, function(r) { fn <- function(xv) num_value_dbl(N(subs(r, x, Float(xv)))); fn })
  dfun <- function(xv, y) c(y[-1L], sum(vapply(seq_len(a), function(i) fx[[i]](xv) * y[i], 0)))
  x0 <- num_value_dbl(N(f$x0)); y <- y0[seq_len(a)]
  out <- list()
  for (p in points) {
    pv <- as.numeric(if (is_sym(p)) num_value_dbl(N(p)) else p)
    h <- pv - x0
    if (method == "Euler") y <- y + h * dfun(x0, y)
    else {
      k1 <- h * dfun(x0, y); k2 <- h * dfun(x0 + h / 2, y + k1 / 2)
      k3 <- h * dfun(x0 + h / 2, y + k2 / 2); k4 <- h * dfun(x0 + h, y + k3)
      y <- y + (k1 + 2 * k2 + 2 * k3 + k4) / 6
    }
    x0 <- pv
    out[[length(out) + 1L]] <- if (derivatives) lapply(y, Float) else Float(y[1L])
  }
  out
}

# ---------------------------------------------------------------------------
# conversion from expressions

.find_conditions <- function(func, x, x0, order, evalf = FALSE, use_limit = TRUE) {
  y0 <- list()
  for (k in seq_len(order)) {
    val <- subs(func, x, x0)
    if (evalf) val <- N(val)
    if (use_limit && is_nan(val)) val <- limit(func, x, x0)
    if (isFALSE(ask_fact(val, "finite")) || is_nan(val)) return(NULL)
    y0[[length(y0) + 1L]] <- val
    func <- .diff1(func, x)
  }
  y0
}

#' @rdname holonomic
#' @export
from_hyper <- function(func, x0 = 0, evalf = FALSE) {
  a <- .hy_ap(func); b <- .hy_bq(func); z <- func$args[[3L]]
  x <- free_symbols(z)[[1L]]
  D <- DifferentialOperators(holonomic_ring(x), "Dx"); Dx <- D$Dx
  xDx <- x * Dx
  r1 <- .diffop(list(S_One), D$algebra)
  for (ai in a) r1 <- r1 * (xDx + ai)
  xDx1 <- xDx - 1
  r2 <- Dx
  for (bi in b) r2 <- r2 * (xDx1 + bi)
  sol <- r1 - r2
  simp <- hyperexpand(func)
  if (is_inf_const(simp)) return(composition(HolonomicFunction(sol, x), z))
  x0 <- as_sym(x0)
  if (!(simp$type == "Function" && simp$name == "hyper")) {
    y0 <- .find_conditions(simp, x, x0, sol$order, use_limit = FALSE)
    while (is.null(y0)) { x0 <- .add(list(x0, S_One)); y0 <- .find_conditions(simp, x, x0, sol$order, use_limit = FALSE) }
    return(composition(HolonomicFunction(sol, x), z, x0, y0))
  }
  x0 <- S_One
  y0 <- .find_conditions(simp, x, x0, sol$order, evalf, use_limit = FALSE)
  while (is.null(y0)) { x0 <- .add(list(x0, S_One)); y0 <- .find_conditions(simp, x, x0, sol$order, evalf, use_limit = FALSE) }
  composition(HolonomicFunction(sol, x), z, x0, y0)
}

#' @rdname holonomic
#' @export
from_meijerg <- function(func, x0 = 0, evalf = FALSE, initcond = TRUE, domain = "QQ") {
  a <- .mg_ap(func); b <- .mg_bq(func)
  n <- length(.mg_an(func)); m <- length(.mg_bm(func)); p <- length(a)
  z <- func$args[[length(func$args)]]
  x <- free_symbols(z)[[1L]]
  D <- DifferentialOperators(holonomic_ring(x, domain), "Dx"); Dx <- D$Dx
  xDx <- x * Dx; xDx1 <- xDx + 1
  r1 <- .diffop(list(.mul(list(x, .pow(S_NegativeOne, .int(m + n - p))))), D$algebra)
  for (ai in a) r1 <- r1 * (xDx1 - ai)
  r2 <- .diffop(list(S_One), D$algebra)
  for (bi in b) r2 <- r2 * (xDx - bi)
  sol <- r1 - r2
  if (!initcond) return(composition(HolonomicFunction(sol, x), z))
  simp <- hyperexpand(func)
  if (is_inf_const(simp)) return(composition(HolonomicFunction(sol, x), z))
  x0 <- as_sym(x0)
  if (!(simp$type == "Function" && simp$name == "meijerg")) {
    y0 <- .find_conditions(simp, x, x0, sol$order, use_limit = FALSE)
    while (is.null(y0)) { x0 <- .add(list(x0, S_One)); y0 <- .find_conditions(simp, x, x0, sol$order, use_limit = FALSE) }
    return(composition(HolonomicFunction(sol, x), z, x0, y0))
  }
  x0 <- S_One
  y0 <- .find_conditions(simp, x, x0, sol$order, evalf, use_limit = FALSE)
  while (is.null(y0)) { x0 <- .add(list(x0, S_One)); y0 <- .find_conditions(simp, x, x0, sol$order, evalf, use_limit = FALSE) }
  composition(HolonomicFunction(sol, x), z, x0, y0)
}

.h_table <- function(x1) {
  D <- DifferentialOperators(holonomic_ring(x1), "Dx"); Dx <- D$Dx
  sp <- sqrt(sym_pi)
  list(
    sin = HolonomicFunction(Dx^2 + 1, x1, 0, list(0, 1)),
    cos = HolonomicFunction(Dx^2 + 1, x1, 0, list(1, 0)),
    exp = HolonomicFunction(Dx - 1, x1, 0, list(1)),
    log = HolonomicFunction(Dx + x1 * Dx^2, x1, 1, list(0, 1)),
    erf = HolonomicFunction(2 * x1 * Dx + Dx^2, x1, 0, list(0, .div(S_Two, sp))),
    erfc = HolonomicFunction(2 * x1 * Dx + Dx^2, x1, 0, list(1, .neg(.div(S_Two, sp)))),
    erfi = HolonomicFunction(-2 * x1 * Dx + Dx^2, x1, 0, list(0, .div(S_Two, sp))),
    sinh = HolonomicFunction(Dx^2 - 1, x1, 0, list(0, 1)),
    cosh = HolonomicFunction(Dx^2 - 1, x1, 0, list(1, 0)),
    sinc = HolonomicFunction(x1 + 2 * Dx + x1 * Dx^2, x1),
    Si = HolonomicFunction(x1 * Dx + 2 * Dx^2 + x1 * Dx^3, x1),
    Ci = HolonomicFunction(x1 * Dx + 2 * Dx^2 + x1 * Dx^3, x1),
    Shi = HolonomicFunction(-x1 * Dx + 2 * Dx^2 + x1 * Dx^3, x1))
}

.convert_poly_rat_alg <- function(func, x, x0 = 0, y0 = NULL, lenics = NULL, domain = "QQ", initcond = TRUE) {
  ispoly <- isTRUE(is_polynomial(func, x))
  israt <- ispoly || .is_ratfun(func, x)
  is_alg <- ispoly || israt
  if (!is_alg) {
    be <- .abe(func); basepoly <- be[[1L]]; ratexp <- be[[2L]]
    if (isTRUE(is_polynomial(basepoly, x)) && is_number_type(ratexp)) {
      if (ratexp$type == "Float") ratexp <- nsimplify(ratexp)
      m <- gmp::numerator(gmp::as.bigq(ratexp$value)); nq <- gmp::denominator(gmp::as.bigq(ratexp$value))
      is_alg <- TRUE
    }
  }
  if (!(ispoly || israt || is_alg)) return(NULL)
  D <- DifferentialOperators(holonomic_ring(x, domain), "Dx"); Dx <- D$Dx
  x0 <- as_sym(x0)
  if (!.has1(func, x)) return(HolonomicFunction(Dx, x, 0, list(func)))
  if (ispoly) {
    sol <- func * Dx - .diff1(func, x)
    sol <- .hnormalize(sol$polys, sol$parent, negative = FALSE)
    sing <- .do_is_singular(sol, x0)
    if (is.null(y0) && is_zero_num(x0) && sing) {
      cl <- rev(.coeff_list(func, x))
      i <- which(!vapply(cl, is_zero_num, NA))[1L]
      y0 <- stats::setNames(list(cl[i:length(cl)]), as.character(i - 1L))
    }
  } else if (israt) {
    pq <- .as_numer_denom(func); p <- pq[[1L]]; q <- pq[[2L]]
    sol <- .mul(list(p, q)) * Dx + .sub(.mul(list(p, .diff1(q, x))), .mul(list(q, .diff1(p, x))))
    sol <- .hnormalize(sol$polys, sol$parent, negative = FALSE)
  } else {
    sol <- .int(nq) * (x / make_rational(gmp::as.bigq(m))) * Dx - 1
    sol <- composition(HolonomicFunction(sol, x), basepoly)$annihilator
    sing <- .do_is_singular(sol, x0)
    if (is.null(y0) && is_zero_num(x0) && sing && (is.null(lenics) || lenics <= 1)) {
      cl <- rev(.coeff_list(basepoly, x))
      i <- which(!vapply(cl, is_zero_num, NA))[1L]
      cf <- .pow(cl[[i]], ratexp)
      y0 <- stats::setNames(list(list(cf)), .str(.mul(list(.int(i - 1L), ratexp))))
    }
  }
  if (!is.null(y0) || !initcond) return(HolonomicFunction(sol, x, x0, y0))
  if (is.null(lenics)) lenics <- sol$order
  if (.do_is_singular(sol, x0)) {
    r <- holonomic_indicial(HolonomicFunction(sol, x, x0))
    mult <- attr(r, "multiplicity")
    if (length(r) == 1L && mult[1L] == 1L) {
      g <- .div(func, .pow(.sub(x, x0), r[[1L]]))
      si <- .find_conditions(g, x, x0, lenics)
      si <- lapply(seq_along(si), function(i) .div(si[[i]], fcall("factorial", list(.int(i - 1L)))))
      return(HolonomicFunction(sol, x, x0, stats::setNames(list(si), .str(r[[1L]]))))
    }
  }
  y0 <- .find_conditions(func, x, x0, lenics)
  while (is.null(y0)) { x0 <- .add(list(x0, S_One)); y0 <- .find_conditions(func, x, x0, lenics) }
  HolonomicFunction(sol, x, x0, y0)
}

#' @rdname holonomic
#' @export
expr_to_holonomic <- function(func, x = NULL, x0 = 0, y0 = NULL, lenics = NULL, domain = NULL, initcond = TRUE) {
  func <- as_sym(func)
  syms <- free_symbols(func)
  if (is.null(x)) {
    if (length(syms) == 1L) x <- syms[[1L]] else stop("Specify the variable for the function", call. = FALSE)
  }
  x <- as_sym(x); x0 <- as_sym(x0)
  if (is.null(domain)) domain <- if (.has_type(func, "Float")) "RR" else "QQ"
  sp <- .convert_poly_rat_alg(func, x, x0 = x0, y0 = y0, lenics = lenics, domain = domain, initcond = initcond)
  if (!is.null(sp)) return(sp)
  if (func$type == "Function") {
    x1 <- Dummy("x_1")
    tab <- .h_table(x1)
    ent <- tab[[func$name]]
    if (is.null(ent)) {
      sol <- .convert_meijerint(func, x, initcond = FALSE, domain = domain)
      if (is.null(sol)) stop("NotImplementedError", call. = FALSE)
      if (!is.null(y0)) sol$y0 <- .hy0(y0)
      if (!is.null(y0) || !initcond) { sol$x0 <- x0; return(sol) }
      if (is.null(lenics)) lenics <- sol$annihilator$order
      yy <- .find_conditions(func, x, x0, lenics)
      while (is.null(yy)) { x0 <- .add(list(x0, S_One)); yy <- .find_conditions(func, x, x0, lenics) }
      return(HolonomicFunction(sol$annihilator, x, x0, yy))
    }
    sol <- change_x(ent, x)
    if (!is.null(y0) || !initcond) {
      sol <- composition(sol, func$args[[1L]])
      if (!is.null(y0)) sol$y0 <- .hy0(y0)
      sol$x0 <- x0
      return(sol)
    }
    if (is.null(lenics)) lenics <- sol$annihilator$order
    yy <- .find_conditions(func, x, x0, lenics)
    while (is.null(yy)) { x0 <- .add(list(x0, S_One)); yy <- .find_conditions(func, x, x0, lenics) }
    return(composition(sol, func$args[[1L]], x0, yy))
  }
  args <- if (func$type %in% c("Add", "Mul")) .sympy_args(func) else func$args
  sol <- expr_to_holonomic(args[[1L]], x = x, initcond = FALSE, domain = domain)
  if (func$type == "Add") for (i in seq_along(args)[-1L]) sol <- sol + expr_to_holonomic(args[[i]], x = x, initcond = FALSE, domain = domain)
  else if (func$type == "Mul") for (i in seq_along(args)[-1L]) sol <- sol * expr_to_holonomic(args[[i]], x = x, initcond = FALSE, domain = domain)
  else if (func$type == "Pow") sol <- .h_pow(sol, .int_val(args[[2L]]))
  sol$x0 <- x0
  if (!is.null(y0)) sol$y0 <- .hy0(y0)
  if (!is.null(y0) || !initcond) return(sol)
  if (.have_ics(sol)) return(sol)
  if (is.null(lenics)) lenics <- sol$annihilator$order
  if (.do_is_singular(sol$annihilator, x0)) {
    r <- holonomic_indicial(sol); mult <- attr(r, "multiplicity")
    if (length(r) == 1L && mult[1L] == 1L) {
      g <- .div(func, .pow(.sub(x, x0), r[[1L]]))
      si <- .find_conditions(g, x, x0, lenics)
      si <- lapply(seq_along(si), function(i) .div(si[[i]], fcall("factorial", list(.int(i - 1L)))))
      return(HolonomicFunction(sol$annihilator, x, x0, stats::setNames(list(si), .str(r[[1L]]))))
    }
  }
  yy <- .find_conditions(func, x, x0, lenics)
  while (is.null(yy)) { x0 <- .add(list(x0, S_One)); yy <- .find_conditions(func, x, x0, lenics) }
  HolonomicFunction(sol$annihilator, x, x0, yy)
}

# _convert_meijerint: a sum of Meijer G functions from meijerint._rewrite1
.convert_meijerint <- function(func, x, initcond = TRUE, domain = "QQ") {
  args <- .mi_rewrite1(func, x)
  if (is.null(args)) return(NULL)
  fac <- args[[1L]]; po <- args[[2L]]; g <- args[[3L]]
  t <- if (po$type == "Pow") list(po$args[[1L]], po$args[[2L]]) else list(po, S_One)
  s0 <- if (identical(t[[1L]]$key, x$key)) t[[2L]] else S_Zero
  shift <- function(G, s) {
    z <- G$args[[3L]]
    if (.contains_key(z, S_I$key)) z <- .map_bottom_up(z, function(w) if (w$type == "Function" && w$name == "exp_polar") fcall("exp", w$args) else w)
    ab <- .mi_coeff_exp(z, x)
    r <- .div(s, ab[[2L]])
    tr <- function(l) lapply(l, function(i) .add(list(i, r)))
    list(.pow(ab[[1L]], .neg(r)), meijerg(list(tr(.mg_an(G)), tr(.mg_aother(G))), list(tr(.mg_bm(G)), tr(.mg_bother(G))), z))
  }
  sol <- NULL
  for (tm in g) {
    cm <- shift(tm[[3L]], .add(list(s0, tm[[2L]])))
    h <- .mul(list(fac, tm[[1L]], cm[[1L]])) * from_meijerg(cm[[2L]], initcond = initcond, domain = domain)
    sol <- if (is.null(sol)) h else sol + h
  }
  sol
}
