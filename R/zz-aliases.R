# SymPy names that map onto existing functionality, and small utilities
# from sympy.core, sympy.utilities, sympy.calculus, sympy.assumptions and
# other modules.

# ---------------------------------------------------------------------------
# Core constructors and helpers

#' Core constructors and helpers
#'
#' `Add()`, `Mul()` and `Pow()` build canonical sums, products and powers.
#' `Equality()`, `Unequality()`, `StrictLessThan()`, `LessThan()`,
#' `StrictGreaterThan()` and `GreaterThan()` are the long names of [Eq()] and
#' friends, and `Rel(a, b, op)` builds a relation from its operator string.
#' `Subs()` is an unevaluated substitution (evaluate with [doit()]).
#' `nfloat()` evaluates the numbers in an expression and keeps symbols.
#' `bottom_up()` applies a function to every node, leaves first; atoms are
#' passed to the function only with `atoms = TRUE`.
#' `default_sort_key()` returns SymPy's `default_sort_key` encoded as a
#' string, so that `order()` of these strings is SymPy's sorted order, and
#' `sym_ordered()` sorts expressions.
#' `comp()` compares two numbers within a tolerance. `check_assumptions()`
#' and `failing_assumptions()` test facts. `sstrrepr()` is [sstr()],
#' `true` and `false` are the boolean atoms.
#'
#' @param ... Expressions, or assumption facts as `name = value`.
#' @param a,b,e Expressions.
#' @param op One of `"=="`, `"!="`, `"<"`, `"<="`, `">"`, `">="`.
#' @param old,new Old and new values for `Subs()`.
#' @param f A function.
#' @param x A list of expressions.
#' @param tol Relative tolerance.
#' @param atoms If `TRUE`, `bottom_up()` also applies `f` to atoms.
#' @param n Significant digits for `nfloat()`.
#' @param exponent If `TRUE`, `nfloat()` also converts exponents.
#' @param dkeys If `TRUE`, `nfloat()` also converts the names of a list.
#' @param evaluate If `FALSE`, `Equality()` keeps the relation unevaluated.
#' @return An expression, list, logical or string.
#' @examples
#' x <- Symbol("x")
#' Add(x, 1, x)
#' doit(Subs(x^2, x, 3))
#' check_assumptions(x^2 + 1, positive = TRUE)
#' @name core-helpers
NULL

#' @rdname core-helpers
#' @export
Add <- function(...) {
  if (any(vapply(list(...), .me_is_matrix, NA))) return(.me_postprocess("Add", list(...)))
  a <- lapply(list(...), as_sym)
  bd <- Filter(.is_bd, a)
  if (length(bd)) return(.bd_add(a, .bd_fam(bd[[1L]])))
  .add(a)
}
#' @rdname core-helpers
#' @export
Mul <- function(...) {
  if (any(vapply(list(...), .me_is_matrix, NA))) return(.me_postprocess("Mul", list(...)))
  a <- lapply(list(...), as_sym)
  bd <- Filter(.is_bd, a)
  if (length(bd)) {
    # Mul.flatten collects equal factors into a power first
    if (anyDuplicated(vapply(bd, .k, ""))) .vec_stop("TypeError", "Power operation is not supported for vectors")
    return(.bd_mul(a, .bd_fam(bd[[1L]])))
  }
  .mul(a)
}
#' @rdname core-helpers
#' @export
Pow <- function(a, b) {
  if (.is_vec(as_sym(a)) || .is_vec(as_sym(b))) .vec_stop("TypeError", "Power operation is not supported for vectors")
  .pow(as_sym(a), as_sym(b))
}
#' @rdname core-helpers
#' @export
Equality <- function(a, b, evaluate = TRUE) Eq(a, b, evaluate)
#' @rdname core-helpers
#' @export
Unequality <- function(a, b) Ne(a, b)
#' @rdname core-helpers
#' @export
StrictLessThan <- function(a, b) Lt(a, b)
#' @rdname core-helpers
#' @export
LessThan <- function(a, b) Le(a, b)
#' @rdname core-helpers
#' @export
StrictGreaterThan <- function(a, b) Gt(a, b)
#' @rdname core-helpers
#' @export
GreaterThan <- function(a, b) Ge(a, b)
#' @rdname core-helpers
#' @export
Rel <- function(a, b, op = "==") {
  switch(op, "==" = Eq(a, b), "!=" = Ne(a, b), "<" = Lt(a, b), "<=" = Le(a, b),
         ">" = Gt(a, b), ">=" = Ge(a, b), stop("unknown relational operator ", op))
}

#' @rdname core-helpers
#' @export
Subs <- function(e, old, new) {
  e <- as_sym(e); old <- as_sym(old); new <- as_sym(new)
  vars <- if (old$type == "Tuple") old$args else list(old)
  pts <- if (new$type == "Tuple") new$args else list(new)
  if (length(vars) != length(pts)) stop("Number of point values must be the same as the number of variables.", call. = FALSE)
  if (!length(pts)) return(e)
  # denest
  if (e$type == "Subs") {
    vp <- .subs_parts(e)
    vars <- c(vp[[1L]], vars); pts <- c(vp[[2L]], pts); e <- e$args[[1L]]
    old <- do.call(Tuple, vars); new <- do.call(Tuple, pts)
  }
  # Subs._hashable_content: the variables replaced by symbols named after
  # their points, so Subs objects differing only in their dummies are equal
  # The placeholders are plain commutative symbols with the default facts.
  reps <- lapply(seq_along(vars), function(i) list(vars[[i]], new_sym("Symbol", name = paste0("_", .str(pts[[i]])),
    extra = list(facts = .symbol_facts(logical()), given = logical()), key = paste0("SubsPt(", pts[[i]]$key, ")"))))
  ce <- xreplace(e, reps)
  extra <- vapply(seq_along(vars), function(i) if (.has1(e, vars[[i]])) "" else paste0(vars[[i]]$key, "=", pts[[i]]$key), "")
  key <- paste0("Subs{", ce$key, "|", paste(sort(extra[nzchar(extra)]), collapse = ","), "}")
  new_sym("Subs", list(e, old, new), key = key)
}
register_rebuilder("Subs", function(e, args) Subs(args[[1L]], args[[2L]], args[[3L]]))
.subs_parts <- function(e) {
  V <- e$args[[2L]]
  P <- e$args[[3L]]
  if (V$type == "Tuple") list(V$args, P$args) else list(list(V), list(P))
}
# Subs.doit: self mappings drop; a Derivative differentiates with respect to
# the variables its expression contains, then substitutes the remaining
# pairs in turn and takes the other derivatives in the order given; any other
# expression is evaluated and substituted. deep = TRUE evaluates the result.
register_doit("Subs", function(e, deep = TRUE) {
  ex <- e$args[[1L]]
  vp <- .subs_parts(e)
  keep <- !vapply(seq_along(vp[[1L]]), function(i) identical(vp[[1L]][[i]]$key, vp[[2L]][[i]]$key), TRUE)
  undone <- Map(list, vp[[1L]][keep], vp[[2L]][keep])
  if (!length(undone)) return(ex)
  if (ex$type == "Derivative") {
    D <- Dummy()
    arg <- ex$args[[1L]]
    undone2 <- list()
    for (u in undone) {
      if (is.null(.free_sym_list(xreplace(ex, list(list(u[[1L]], D))))[[D$key]])) {
        if (.has1(arg, u[[1L]])) ex <- subs(ex, u[[1L]], u[[2L]])
      } else undone2[[length(undone2) + 1L]] <- u
    }
    undone <- undone2
    wrt <- list()
    D <- Dummy()
    expr <- ex$args[[1L]]
    free <- names(.free_sym_list(expr))
    for (tp in ex$args[-1L]) {
      vi <- tp$args[[1L]]
      if ((vi$type %in% c("Symbol", "Dummy") && vi$key %in% free) ||
          !is.null(.free_sym_list(subs(expr, vi, D))[[D$key]]))
        expr <- diff(expr, list(vi, tp$args[[2L]]))
      else wrt[[length(wrt) + 1L]] <- tp
    }
    rv <- if (length(undone)) subs(expr, undone) else expr
    for (tp in wrt) rv <- diff(rv, list(tp$args[[1L]], tp$args[[2L]]))
  } else {
    rv <- subs(doit(ex, deep = deep), undone)
  }
  if (deep && !identical(rv$key, e$key)) rv <- doit(rv)
  rv
})
register_diff_hook("Subs", function(e, s) {
  f <- e$args[[1L]]
  vp <- .subs_parts(e)
  terms <- list()
  for (k in seq_along(vp[[1L]])) {
    dp <- .diff1(vp[[2L]][[k]], s)
    if (is_zero_num(dp)) next
    terms[[length(terms) + 1L]] <- .mul(list(dp, doit(Subs(.diff1(f, vp[[1L]][[k]]), e$args[[2L]], e$args[[3L]]))))
  }
  vk <- vapply(vp[[1L]], .k, "")
  free <- setdiff(vapply(free_symbols(f), .k, ""), vk)
  if (length(intersect(free, vapply(free_symbols(s), .k, ""))))
    terms[[length(terms) + 1L]] <- doit(Subs(.diff1(f, s), e$args[[2L]], e$args[[3L]]))
  .add(terms)
})
register_str_printer("Subs", function(e) paste0("Subs(", .str(e$args[[1L]]), ", ", .str(e$args[[2L]]), ", ",
                                                .str(e$args[[3L]]), ")"))
register_bound_vars("Subs", function(e) e$args[[2L]]$key)

#' @rdname core-helpers
#' @export
nfloat <- function(e, n = 15, exponent = FALSE, dkeys = FALSE) {
  # SymPy's nfloat: Rationals become Floats except in exponents (unless
  # exponent = TRUE) and in undefined functions
  rec <- function(z) nfloat(z, n = n, exponent = exponent, dkeys = dkeys)
  if (inherits(e, "symmatrix")) return(.mat_map(e, rec))
  if (is.list(e) && !is_sym(e)) {
    nm <- names(e)
    out <- lapply(e, rec)
    if (dkeys && !is.null(nm)) names(out) <- vapply(nm, function(k) .str(rec(sympify(k))), "")
    return(out)
  }
  rv <- as_sym(e)
  if (is_number_type(rv)) return(Float(rv, dps = n))
  if (!length(free_symbols(rv)) && .is_expr_number(rv)) {
    rv <- N(rv, n)
    if (is_number_type(rv)) rv <- Float(N(rv, n), dps = n)
    return(rv)
  }
  if (!length(rv$args)) return(rv)
  if (is_relational(rv)) return(sym_rebuild(rv, lapply(rv$args, rec)))
  ros <- atoms(rv, "CRootOf")
  if (length(ros)) rv <- xreplace(rv, lapply(ros, function(ro) list(ro, N(ro, n))))
  if (!exponent) {
    pows <- atoms(rv, "Pow")
    reps <- lapply(pows, function(p) list(p, .pow(p$args[[1L]], Dummy())))
    rv <- xreplace(rv, reps)
  }
  rv <- N(rv, n)
  if (!exponent) {
    rv <- xreplace(rv, lapply(seq_along(reps), function(i) list(reps[[i]][[2L]]$args[[2L]], pows[[i]]$args[[2L]])))
  } else {
    rv <- .xreplace_transform(rv, function(x) x$type == "Pow" && x$args[[2L]]$type == "Integer",
                              function(x) .pow(x$args[[1L]], Float(x$args[[2L]], dps = n)))
  }
  .xreplace_transform(rv, function(x) x$type == "Function" && !.is_undefined_application(x),
                      function(x) sym_rebuild(x, lapply(x$args, function(a) nfloat(a, n = n, exponent = exponent))))
}

# xreplace with a Transform: the first node (top down) that passes `keep` is
# replaced by f(node) and not searched further.
.xreplace_transform <- function(e, keep, f) {
  if (keep(e)) return(f(e))
  if (!length(e$args)) return(e)
  sym_rebuild(e, lapply(e$args, function(a) if (is_sym(a)) .xreplace_transform(a, keep, f) else a))
}

# Expr.is_number: no free symbols and not a Boolean, set or other non-Expr.
.is_expr_number <- function(e) {
  !(e$type %in% c(.BOOL_TYPES, .REL_TYPES)) && !is_set(e)
}

#' @rdname core-helpers
#' @export
bottom_up <- function(e, f, atoms = FALSE) {
  # sympy.simplify.simplify.bottom_up: atoms are passed to f only on request
  rec <- function(z) {
    if (!length(z$args)) return(if (atoms) as_sym(f(z)) else z)
    args <- lapply(z$args, rec)
    if (!identical(vapply(args, .k, ""), vapply(z$args, .k, ""))) z <- sym_rebuild(z, args)
    as_sym(f(z))
  }
  rec(as_sym(e))
}

#' @rdname core-helpers
#' @export
default_sort_key <- function(e) .sort_key_string(.sort_key(as_sym(e)))

# SymPy's sort_key tuple as a string over [0-9a-z] whose order (in any
# collation) is the tuple order: each item is tagged (1 number, 2 string,
# 3 tuple) and strings and tuples end with a marker below every item.
.sort_key_string <- function(k) {
  if (is.list(k)) return(paste0("3", paste(vapply(k, .sort_key_string, ""), collapse = ""), "0"))
  if (is.character(k)) {
    cp <- utf8ToInt(k)
    return(paste0("2", paste(vapply(cp, function(i) .base36(i, 4L), ""), collapse = ""), "0000"))
  }
  if (!inherits(k, "bigq") && is.infinite(k)) return(if (k > 0) "14" else "10")
  q <- gmp::as.bigq(k)
  if (q == 0) return("12")
  if (q > 0) return(paste0("13", .sort_key_pos(q)))
  d <- strsplit(.sort_key_pos(-q), "")[[1L]]
  paste0("11", paste(9L - as.integer(d), collapse = ""))
}

# a positive rational as a fixed width decimal exponent and 30 digits
.sort_key_pos <- function(q) {
  e <- floor(log10(as.double(q)))
  p10 <- function(n) if (n >= 0) gmp::as.bigq(gmp::as.bigz(10)^n) else 1 / gmp::as.bigq(gmp::as.bigz(10)^(-n))
  while (p10(e) > q) e <- e - 1
  while (p10(e + 1) <= q) e <- e + 1
  m <- gmp::as.bigz(floor(q / p10(e - 29)))
  paste0(formatC(e + 50000, width = 6, flag = "0", format = "d"), as.character(m))
}

.base36 <- function(i, width) {
  d <- character(width)
  for (j in width:1) {
    d[j] <- c(0:9, letters)[i %% 36 + 1]
    i <- i %/% 36
  }
  paste(d, collapse = "")
}

#' @noRd
ordered <- function(x, ...) {
  # SymPy's ordered(): node count, then default_sort_key
  if (is.list(x) && length(x) && all(vapply(x, is_sym, TRUE))) return(.bool_print_order(x))
  base::ordered(x, ...)
}

#' @rdname core-helpers
#' @export
comp <- function(a, b, tol = 1e-12) {
  a <- as_sym(a); b <- as_sym(b)
  if (is_exact_number(a) && is_exact_number(b)) {
    d <- abs(a$value - b$value)
    return(d <= gmp::as.bigq(tol) * max(gmp::as.bigq(1), abs(b$value)))
  }
  va <- .neval(a, list()); vb <- .neval(b, list())
  base::Mod(va - vb) <= tol * max(1, base::Mod(vb))
}

#' @rdname core-helpers
#' @export
check_assumptions <- function(e, ...) {
  want <- list(...)
  f <- .facts(as_sym(e))
  res <- vapply(names(want), function(nm) {
    v <- f[[nm]]
    if (is.na(v)) NA else v == isTRUE(want[[nm]])
  }, NA)
  if (any(res %in% FALSE)) return(FALSE)
  if (anyNA(res)) return(NA)
  TRUE
}

#' @rdname core-helpers
#' @export
failing_assumptions <- function(e, ...) {
  want <- list(...)
  f <- .facts(as_sym(e))
  out <- list()
  for (nm in names(want)) {
    v <- f[[nm]]
    if (is.na(v) || v != isTRUE(want[[nm]])) out[[nm]] <- v
  }
  out
}

#' @rdname core-helpers
#' @export
sstrrepr <- function(e) sstr(e)

#' @rdname core-helpers
#' @export
true <- S_true

#' @rdname core-helpers
#' @export
false <- S_false

# ---------------------------------------------------------------------------
# Polynomial aliases and constructions

#' Polynomial aliases and constructions
#'
#' `sym_poly()` is [Poly()], `PurePoly()` is a [Poly()] whose generators
#' are the placeholders `_0`, `_1`, ... so that polynomials differing only in
#' their generators compare equal, `GroebnerBasis()` is
#' [groebner()], and `rootof()`, `RootOf()` and `ComplexRootOf()` are
#' [CRootOf()]. `monic()` divides by the leading coefficient, `cofactors()`
#' returns the gcd and both cofactors, `half_gcdex()` the Bezout coefficient
#' and gcd, `pexquo()` exact pseudo division, `all_roots()` all roots with
#' multiplicity, solving linear and quadratic factors and binomials by
#' radicals and giving `CRootOf` for the other factors. `cyclotomic_poly()`,
#' `symmetric_poly()`, `interpolating_poly()` (in Lagrange form), `random_poly()`,
#' `itermonomials()` and `viete()` build polynomials and formulas;
#' `legendre_poly()`, `chebyshevt_poly()`, `chebyshevu_poly()`,
#' `hermite_poly()`, `hermite_prob_poly()`, `laguerre_poly()` and
#' `jacobi_poly()` return the explicit orthogonal polynomials. `lex`,
#' `grlex` and `grevlex` name the monomial orders.
#'
#' @param x,f,g Expressions.
#' @param ... Generators.
#' @param n Degree, or number of variables.
#' @param k Order of the symmetric polynomial.
#' @param vars A list of symbols.
#' @param X,Y Interpolation nodes and values.
#' @param inf,sup Coefficient range for `random_poly()`.
#' @param max_degree,min_degree Degree bounds.
#' @param roots A list of root symbols for `viete()`.
#' @param a,b Jacobi parameters.
#' @return An expression, list or `Poly`.
#' @examples
#' x <- Symbol("x")
#' cyclotomic_poly(6, x)
#' cofactors(x^2 - 1, x - 1)
#' viete(x^2 + Symbol("p") * x + Symbol("q"), x)
#' @name poly-aliases
NULL

#' @noRd
poly <- function(x, ...) if (is_sym(x)) Poly(x, ...) else stats::poly(x, ...)
#' @rdname poly-aliases
#' @export
PurePoly <- function(f, ...) .purepoly(f, ...)
#' @rdname poly-aliases
#' @export
GroebnerBasis <- function(f, ...) groebner(f, ...)
#' @rdname poly-aliases
#' @export
rootof <- function(f, x, k) CRootOf(f, x, k)
#' @rdname poly-aliases
#' @export
RootOf <- rootof
#' @rdname poly-aliases
#' @export
ComplexRootOf <- rootof
#' @rdname poly-aliases
#' @export
lex <- "lex"
#' @rdname poly-aliases
#' @export
grlex <- "grlex"
#' @rdname poly-aliases
#' @export
grevlex <- "grevlex"

#' @rdname poly-aliases
#' @export
monic <- function(f, ...) {
  f <- as_sym(f)
  if (is_zero_num(f)) return(S_Zero)
  lc <- LC(f, ...)
  expand(.div(f, lc))
}

#' @rdname poly-aliases
#' @export
cofactors <- function(f, g, ...) {
  h <- gcd(f, g, ...)
  list(h = h, cff = cancel(.div(as_sym(f), h)), cfg = cancel(.div(as_sym(g), h)))
}

#' @rdname poly-aliases
#' @export
half_gcdex <- function(f, g, ...) {
  r <- gcdex(f, g, ...)
  list(s = r$s, h = r$h)
}

#' @rdname poly-aliases
#' @export
pexquo <- function(f, g, ...) {
  d <- pdiv(f, g, ...)
  r <- if (inherits(d$r, "Poly")) as_expr(d$r) else d$r
  if (!is_zero_num(r)) stop("pseudo division is not exact")
  d$q
}

#' @rdname poly-aliases
#' @export
all_roots <- function(f, x = NULL) {
  if (inherits(f, "RootSystem")) return(.rs_all_roots(f))
  if (!is.null(r <- .all_roots_rootof(f, x))) return(r)
  r <- roots(f, x, multiple = TRUE, rootof = TRUE)
  if (is.null(r)) list() else r
}

#' @rdname poly-aliases
#' @export
cyclotomic_poly <- function(n, x = Symbol("x")) {
  x <- as_sym(x)
  n <- as.integer(n)
  cache <- list()
  phi <- function(m) {
    key <- as.character(m)
    if (!is.null(cache[[key]])) return(cache[[key]])
    p <- .sub(.pow(x, .int(m)), S_One)
    for (d in seq_len(m - 1L)) if (m %% d == 0L) p <- quo(p, phi(d), x)
    p <- expand(p)
    cache[[key]] <<- p
    p
  }
  phi(n)
}

#' @rdname poly-aliases
#' @export
symmetric_poly <- function(k, vars) {
  vars <- lapply(vars, as_sym)
  if (k == 0L) return(S_One)
  .add(lapply(combinations(seq_along(vars), k), function(ix) .mul(vars[ix])))
}

#' @rdname poly-aliases
#' @export
interpolating_poly <- function(n, x, X = NULL, Y = NULL) {
  x <- as_sym(x)
  if (is.null(X)) X <- lapply(seq_len(n) - 1L, function(i) Symbol(paste0("x", i)))
  if (is.null(Y)) Y <- lapply(seq_len(n) - 1L, function(i) Symbol(paste0("y", i)))
  .lagrange_poly(n, x, lapply(X, as_sym), lapply(Y, as_sym))
}

#' @rdname poly-aliases
#' @export
random_poly <- function(x, n, inf = -10, sup = 10) {
  x <- as_sym(x)
  pool <- inf:sup
  cf <- pool[sample.int(length(pool), n + 1L, replace = TRUE)]
  .add(lapply(0:n, function(k) .mul(list(.int(cf[n - k + 1L]), .pow(x, .int(k))))))
}

#' @rdname poly-aliases
#' @export
itermonomials <- function(vars, max_degree, min_degree = 0) {
  vars <- lapply(vars, as_sym)
  out <- list()
  rec <- function(i, left, cur) {
    if (i > length(vars)) {
      d <- max_degree - left
      if (d >= min_degree) out[[length(out) + 1L]] <<- cur
      return(invisible())
    }
    for (e in 0:left) rec(i + 1L, left - e, .mul(list(cur, .pow(vars[[i]], .int(e)))))
  }
  rec(1L, max_degree, S_One)
  .sort_by_key(out[!duplicated(vapply(out, .k, ""))])
}

#' @rdname poly-aliases
#' @export
viete <- function(f, x = NULL, roots = NULL) {
  f <- as_sym(f)
  if (is.null(x)) x <- free_symbols(f)[[1L]]
  cl <- .coeff_list(f, as_sym(x))
  n <- length(cl) - 1L
  if (is.null(roots)) roots <- lapply(seq_len(n), function(i) Symbol(paste0("r", i)))
  lapply(seq_len(n), function(k) list(symmetric_poly(k, roots), cancel(.mul(list(.int((-1)^k), cl[[k + 1L]], .inv(cl[[1L]]))))))
}

#' @rdname poly-aliases
#' @export
legendre_poly <- function(n, x = Symbol("x")) legendre(n, x)
#' @rdname poly-aliases
#' @export
chebyshevt_poly <- function(n, x = Symbol("x")) chebyshevt(n, x)
#' @rdname poly-aliases
#' @export
chebyshevu_poly <- function(n, x = Symbol("x")) chebyshevu(n, x)
#' @rdname poly-aliases
#' @export
hermite_poly <- function(n, x = Symbol("x")) hermite(n, x)
#' @rdname poly-aliases
#' @export
hermite_prob_poly <- function(n, x = Symbol("x")) hermite_prob(n, x)
#' @rdname poly-aliases
#' @export
laguerre_poly <- function(n, x = Symbol("x")) laguerre(n, x)
#' @rdname poly-aliases
#' @export
jacobi_poly <- function(n, a, b, x = Symbol("x")) expand(jacobi(n, a, b, x))

# ---------------------------------------------------------------------------
# Function aliases

#' Function aliases and small special functions
#'
#' `sym_rf()` and `ff()` are [RisingFactorial()] and [FallingFactorial()],
#' `ln()` is the natural logarithm, `E1(z)` the exponential integral
#' `expint(1, z)`, `Li(x)` the offset logarithmic integral `li(x) - li(2)`,
#' `erf2(x, y)` is `erf(y) - erf(x)`, `Id(x)` the identity, and
#' `LeviCivita()` (alias `Eijk()`) the permutation symbol of integer
#' arguments, extended as in SymPy by
#' `prod((a[j] - a[i]) / sym_factorial(i))`. `chebyshevt_root(n, k)` and `chebyshevu_root(n, k)` give the
#' `k`-th root (from 0) of the Chebyshev polynomials. `transpose()` and
#' `adjoint()` act on matrices.
#'
#' @param x,y,z,k,n Arguments.
#' @param ... Integer arguments.
#' @param M A matrix.
#' @return An expression or matrix.
#' @examples
#' sym_rf(Symbol("x"), 3)
#' LeviCivita(1, 2, 3)
#' chebyshevt_root(3, 0)
#' @name function-aliases
NULL

#' @noRd
rf <- function(x, k, ...) {
  if (missing(x) || missing(k) || ...length()) return(.base_call(quote(stats::rf)))
  fcall("RisingFactorial", list(as_sym(x), as_sym(k)))
}
#' @rdname function-aliases
#' @export
ff <- function(x, k) fcall("FallingFactorial", list(as_sym(x), as_sym(k)))
#' @rdname function-aliases
#' @export
ln <- function(x) fcall("log", list(as_sym(x)))
#' @rdname function-aliases
#' @export
E1 <- function(z) fcall("expint", list(S_One, as_sym(z)))
#' @rdname function-aliases
#' @export
Li <- function(x) fcall("Li", list(as_sym(x)))
#' @rdname function-aliases
#' @export
erf2 <- function(x, y) fcall("erf2", list(as_sym(x), as_sym(y)))

# Li.eval and Li.fdiff; values through li(z) - li(2)
def_function("Li",
  eval = function(args) {
    z <- args[[1L]]
    if (is_const(z, "oo")) return(S_Infinity)
    if (identical(z$key, S_Two$key)) return(S_Zero)
    NULL
  },
  diff = function(args, i) .inv(fcall("log", args)),
  numeric = function(z) .FUNS$li$numeric(z) - .FUNS$li$numeric(2))

# erf2.eval and erf2.fdiff
def_function("erf2", nargs = 2L,
  eval = function(args) {
    x <- args[[1L]]; y <- args[[2L]]
    d <- function() .sub(fcall("erf", list(y)), fcall("erf", list(x)))
    if (is_nan(x) || is_nan(y)) return(S_NaN)
    if (identical(x$key, y$key)) return(S_Zero)
    chk <- c(S_Infinity$key, S_NegativeInfinity$key, S_Zero$key)
    if (x$key %in% chk || y$key %in% chk) return(d())
    if (y$type == "Function" && y$name == "erf2inv" && identical(y$args[[1L]]$key, x$key)) return(y$args[[2L]])
    rinf <- function(v) isTRUE(ask_fact(v, "extended_real")) && isTRUE(ask_fact(v, "infinite"))
    if (isTRUE(ask_fact(x, "zero")) || isTRUE(ask_fact(y, "zero")) || rinf(x) || rinf(y)) return(d())
    sx <- could_extract_minus_sign(x); sy <- could_extract_minus_sign(y)
    if (sx && sy) return(.neg(fcall("erf2", list(.neg(x), .neg(y)))))
    if (sx || sy) return(d())
    NULL
  },
  diff = function(args, i) {
    v <- args[[i]]
    .div(.mul(list(if (i == 1L) .int(-2) else S_Two, fcall("exp", list(.neg(.sq(v)))))), .sqrt(S_Pi))
  },
  numeric = function(x, y) .FUNS$erf$numeric(y) - .FUNS$erf$numeric(x))
.expand_hooks$erf2 <- function(e) .sub(fcall("erf", list(e$args[[2L]])), fcall("erf", list(e$args[[1L]])))
#' @rdname function-aliases
#' @export
Id <- function(x) as_sym(x)
#' @rdname function-aliases
#' @export
LeviCivita <- function(...) {
  a <- lapply(list(...), as_sym)
  if (!length(a)) return(S_One)
  if (has_dups(a)) return(S_Zero)
  if (!all(vapply(a, is_int_num, TRUE))) return(fnode("LeviCivita", a))
  v <- lapply(a, function(z) gmp::as.bigz(z$value))
  n <- length(v)
  num <- gmp::as.bigz(1L)
  for (i in seq_len(n - 1L)) for (j in (i + 1L):n) num <- num * (v[[j]] - v[[i]])
  make_rational(gmp::as.bigq(num, prod(gmp::factorialZ(seq_len(n) - 1L))))
}
#' @rdname function-aliases
#' @export
Eijk <- LeviCivita
#' @rdname function-aliases
#' @export
chebyshevt_root <- function(n, k) fcall("cos", list(.mul(list(S_Pi, .rat(2L * k + 1L, 2L * n)))))
#' @rdname function-aliases
#' @export
chebyshevu_root <- function(n, k) fcall("cos", list(.mul(list(S_Pi, .rat(k + 1L, n + 1L)))))
#' @rdname function-aliases
#' @export
transpose <- function(M) t(Matrix(M))
#' @rdname function-aliases
#' @export
adjoint <- function(M) {
  if (is_sym(M) && .is_matexpr(M)) return(.madjoint(M))
  t(.mat_map(Matrix(M), function(z) fcall("conjugate", list(z))))
}

# ---------------------------------------------------------------------------
# Series, simplification and transform helpers

#' Series and summation helpers
#'
#' `limit_seq(e, n)` the limit as `n`
#' tends to infinity, `difference_delta(e, n, step)` the forward
#' difference, and `residue(e, x, x0)` the residue at a pole (the
#' coefficient of `1/(x - x0)` in the Laurent series). `hypersimp(f, k)`
#' returns the rational term ratio `f(k + 1)/f(k)` of a hypergeometric term
#' (or `NULL`), and `hypersimilar(f, g, k)` tests whether `f/g` is rational.
#'
#' @param e,f,g Expressions.
#' @param x,n,k Variables.
#' @param x0 Point.
#' @param step Step size.
#' @return An expression, `NULL` or a logical.
#' @examples
#' x <- Symbol("x")
#' residue(1 / (x^2 + 1), x, sym_I)
#' hypersimp(sym_factorial(Symbol("k")), Symbol("k"))
#' @name series-helpers
NULL

#' @rdname series-helpers
#' @export
difference_delta <- function(e, n = NULL, step = 1) {
  e <- as_sym(e); step <- as_sym(step)
  if (is.null(n)) {
    fs <- free_symbols(e)
    if (!length(fs)) return(S_Zero)
    if (length(fs) > 1L) stop("Since there is more than one variable in the expression, a variable must be supplied to take the difference of ", .str(e), call. = FALSE)
    n <- fs[[1L]]
  }
  if (identical(ask_fact(step, "finite"), FALSE) || length(free_symbols(step)))
    stop("Step should be a finite number.", call. = FALSE)
  .difference_delta(e, as_sym(n), step)
}
#' @rdname series-helpers
#' @export
residue <- function(e, x, x0 = 0) {
  # sympy/series/residues.py
  e <- as_sym(e); x <- as_sym(x); x0 <- as_sym(x0)
  if (!is_zero_num(x0)) e <- subs(e, x, .add(list(x, x0)))
  for (n in c(0L, 1L, 2L, 4L, 8L, 16L, 32L)) {
    s <- nseries(e, x, n = n)
    if (!.has_type(s, "Order") || .is_true(Ge(.so_getn(s), 0))) break
  }
  s <- collect(.so_removeO(s), x)
  terms <- if (s$type == "Add") s$args else list(s)
  res <- S_Zero
  for (a in terms) {
    cm <- .as_coeff_mul_dep(a, x)
    m <- .mul(cm[[2L]])
    if (!(is_one_num(m) || identical(m$key, x$key) || (m$type == "Pow" && m$args[[2L]]$type == "Integer")))
      .notimpl_error(paste("term of unexpected form:", .str(m)))
    if (identical(m$key, .inv(x)$key)) res <- .add(list(res, cm[[1L]]))
  }
  res
}
#' @rdname series-helpers
#' @export
hypersimp <- function(f, k) .hypersimp(as_sym(f), as_sym(k))
#' @rdname series-helpers
#' @export
hypersimilar <- function(f, g, k) {
  # SymPy's hypersimilar: (f/g) rewritten with gamma, expand_func, rational in k
  h <- expand_func(.rewrite_gamma_sympy(.div(as_sym(f), as_sym(g))))
  .is_ratfun(h, as_sym(k))
}



# ---------------------------------------------------------------------------
# Calculus utilities (sympy.calculus)

#' Calculus utilities
#'
#' `singularities(f, x, domain)` returns the points of `domain` where `f` is
#' undefined (zeros of denominators and of logarithm arguments); the default
#' domain is `Reals` for a real `x` and `Complexes` otherwise. `stationary_points()`
#' solves `f'(x) = 0` on a domain; `minimum()` and `maximum()` return the
#' extreme values on a closed interval or the real line, including
#' one-sided limits at poles.
#' `is_increasing()`, `is_strictly_increasing()`, `is_decreasing()`,
#' `is_strictly_decreasing()`, `is_monotonic()` and `is_convex()` test the
#' sign of derivatives on an interval. `periodicity(f, x)` returns the
#' period of trigonometric expressions. `euler_equations(L, funcs, vars)`
#' gives the list of Euler-Lagrange equations of a Lagrangian, for
#' derivatives of any order. `finite_diff_weights(order, x_list, x0)` returns
#' Fornberg weights `w[[m + 1]][[k]][[j]]`: derivative order `m`, first `k`
#' grid points, point `j`. `apply_finite_diff()` applies the full-grid
#' weights.
#'
#' @param f,L Expressions.
#' @param x The variable.
#' @param domain An interval, `Reals` or `Complexes`; for `singularities()`,
#'   `NULL` picks from the assumptions on `x`.
#' @param interval An interval, or `Reals`.
#' @param funcs Function applications such as `u(x)`, a list of them, or
#'   `NULL` for all undefined functions in `L`.
#' @param vars Independent variables, or `NULL` for the arguments of the
#'   first function.
#' @param order Derivative order.
#' @param x_list,y_list Grid points and values.
#' @param x0 Evaluation point.
#' @return A set, expression, logical or list.
#' @examples
#' x <- Symbol("x")
#' singularities(1 / (x^2 - 1), x)
#' minimum(x^2 - 2*x, x)
#' is_increasing(exp(x), Reals, x)
#' finite_diff_weights(1, list(-1, 0, 1), 0)[[2]][[3]]
#' @name calculus-utils
NULL

#' @rdname calculus-utils
#' @export
singularities <- function(f, x, domain = NULL) {
  # sympy.calculus.singularities.singularities
  f <- as_sym(f); x <- as_sym(x)
  if (is.null(domain)) domain <- if (isTRUE(ask_fact(x, "real"))) Reals else Complexes
  domain <- as_sym(domain)
  trig <- c("tan", "cot", "sec", "csc", "tanh", "coth", "sech", "csch")
  has_trig <- any(vapply(atoms(f, "Function"), function(g) g$name %in% trig && .has1(g, x), TRUE))
  if (!has_trig) {
    pts <- FiniteSet(.sort_set_elems(.real_singularities(f, x, real = FALSE)))
    return(if (identical(domain$key, Complexes$key)) pts else Intersection(pts, domain))
  }
  # the poles of the trigonometric and hyperbolic functions come from the
  # negative powers of their cos and cosh rewrites
  e <- rewrite(rewrite(f, "cos"), "cosh")
  notimpl <- function() .sset_err("NotImplementedError", .filldedent(
    "Methods for determining the singularities of this function have not been developed."))
  sings <- S_EmptySet
  for (i in Filter(function(a) a$type == "Pow", preorder_traversal(e))) {
    b <- i$args[[1L]]; ex <- i$args[[2L]]
    if (is_zero_num(b)) { sings <- Union(sings, solveset(Lt(ex, 0), x, domain)); next }
    if (isTRUE(ask_fact(ex, "infinite"))) notimpl()
    if (isTRUE(ask_fact(ex, "negative"))) sings <- Union(sings, solveset(b, x, domain))
    else if (.has1(ex, x)) sings <- Union(sings, Intersection(solveset(b, x, domain), solveset(Lt(ex, 0), x, domain)))
  }
  fa <- Filter(function(a) is_sym(a) && a$type == "Function", preorder_traversal(f))
  for (i in Filter(function(a) a$name %in% c("log", "asech", "acsch"), fa)) sings <- Union(sings, solveset(i$args[[1L]], x, domain))
  for (i in Filter(function(a) a$name %in% c("atanh", "acoth"), fa)) {
    sings <- Union(sings, solveset(.sub(i$args[[1L]], S_One), x, domain))
    sings <- Union(sings, solveset(.add(list(i$args[[1L]], S_One)), x, domain))
  }
  sings
}

#' @rdname calculus-utils
#' @export
stationary_points <- function(f, x, domain = Reals) solveset(.diff1(as_sym(f), as_sym(x)), as_sym(x), domain)

.extreme <- function(f, x, domain, sgn) {
  f <- as_sym(f); x <- as_sym(x); domain <- as_sym(domain)
  pick <- function(vals) {
    b <- vals[[1L]]
    for (v in vals[-1L]) if (.sgn_cmp(v, b) == sgn) b <- v
    b
  }
  if (domain$type == "FiniteSet") return(pick(lapply(domain$args, function(p) simplify(subs(f, x, p)))))
  if (domain$type == "EmptySet") stop("empty domain")
  if (!(domain$type == "Interval" || identical(domain$key, Reals$key))) stop("unsupported domain for extrema")
  if (!.has1(f, x)) return(f)
  # function_range: a periodic function is examined over one period
  per <- periodicity(f, x)
  if (!is.null(per) && !is_zero_num(per) && (identical(domain$key, Reals$key) ||
      (domain$type == "Interval" && (is_inf_const(domain$args[[1L]]) || is_inf_const(domain$args[[2L]])))))
    domain <- Interval(S_Zero, per)
  cands <- list()
  st <- stationary_points(f, x, domain)
  if (st$type == "FiniteSet") {
    cands <- st$args
  } else if (st$type != "EmptySet") {
    per <- periodicity(f, x)
    bounded <- domain$type == "Interval" && !is_inf_const(domain$args[[1L]]) && !is_inf_const(domain$args[[2L]])
    if (is.null(per) || !bounded) stop("could not find all stationary points")
    base <- .solve_univariate(.diff1(f, x), x, check = FALSE)
    lo <- as.numeric(N(domain$args[[1L]])); hi <- as.numeric(N(domain$args[[2L]])); Tn <- as.numeric(N(per))
    for (b in base) {
      bn <- as.numeric(N(b))
      for (k in seq(floor((lo - bn) / Tn) - 1, ceiling((hi - bn) / Tn) + 1)) {
        p <- simplify(.add(list(b, .mul(list(.int(k), per)))))
        pn <- as.numeric(N(p))
        if (pn >= lo - 1e-12 && pn <= hi + 1e-12) cands[[length(cands) + 1L]] <- p
      }
    }
  }
  ends <- if (domain$type == "Interval") list(domain$args[[1L]], domain$args[[2L]]) else list(S_NegativeInfinity, S_Infinity)
  lo <- as.numeric(N(ends[[1L]])); hi <- as.numeric(N(ends[[2L]]))
  vals <- lapply(cands, function(p) simplify(subs(f, x, p)))
  # Approach poles and endpoints from inside the domain.
  side <- function(p, dir) {
    if (is_inf_const(p)) return(limit(f, x, p))
    v <- simplify(subs(f, x, p))
    if ((v$type == "Constant" && v$name == "nan") || is_inf_const(v)) limit(f, x, p, dir) else v
  }
  vals <- c(vals, list(side(ends[[1L]], "+"), side(ends[[2L]], "-")))
  for (p in .real_singularities(f, x)) {
    pn <- as.numeric(N(p))
    if (pn > lo && pn <= hi) vals[[length(vals) + 1L]] <- limit(f, x, p, "-")
    if (pn >= lo && pn < hi) vals[[length(vals) + 1L]] <- limit(f, x, p, "+")
  }
  ok <- vapply(vals, function(v) !is.na(tryCatch(as.numeric(N(v)), error = function(e) NA)), TRUE)
  if (!all(ok)) stop("could not compare candidate values")
  pick(vals)
}

#' @rdname calculus-utils
#' @export
minimum <- function(f, x, domain = Reals) .extreme(f, x, domain, -1)
#' @rdname calculus-utils
#' @export
maximum <- function(f, x, domain = Reals) .extreme(f, x, domain, 1)

.first_free <- function(f) {
  fs <- free_symbols(f)
  if (length(fs)) fs[[1L]] else Dummy("x")
}

.sign_on <- function(g, x, interval, rel) {
  s <- solve_univariate_inequality(.make_rel(rel, g, S_Zero, evaluate = FALSE), x, relational = FALSE,
                                   domain = as_sym(interval))
  identical(Complement(as_sym(interval), s)$type, "EmptySet") ||
    identical(Complement(as_sym(interval), s)$key, S_EmptySet$key)
}

#' @rdname calculus-utils
#' @export
is_increasing <- function(f, interval = Reals, x = NULL) {
  f <- as_sym(f); if (is.null(x)) x <- .first_free(f)
  .sign_on(.diff1(f, as_sym(x)), as_sym(x), interval, "GreaterThan")
}
#' @rdname calculus-utils
#' @export
is_strictly_increasing <- function(f, interval = Reals, x = NULL) {
  f <- as_sym(f); if (is.null(x)) x <- .first_free(f)
  .sign_on(.diff1(f, as_sym(x)), as_sym(x), interval, "StrictGreaterThan")
}
#' @rdname calculus-utils
#' @export
is_decreasing <- function(f, interval = Reals, x = NULL) {
  f <- as_sym(f); if (is.null(x)) x <- .first_free(f)
  .sign_on(.diff1(f, as_sym(x)), as_sym(x), interval, "LessThan")
}
#' @rdname calculus-utils
#' @export
is_strictly_decreasing <- function(f, interval = Reals, x = NULL) {
  f <- as_sym(f); if (is.null(x)) x <- .first_free(f)
  .sign_on(.diff1(f, as_sym(x)), as_sym(x), interval, "StrictLessThan")
}
#' @rdname calculus-utils
#' @export
is_monotonic <- function(f, interval = Reals, x = NULL) {
  # sympy.calculus.singularities.is_monotonic: no interior singularity and
  # no turning point (a zero of the derivative) in the interval
  f <- as_sym(f)
  interval <- as_sym(interval)
  fs <- free_symbols(f)
  if (is.null(x) && length(fs) > 1L)
    stop(.sympy_condition("NotImplementedError", "is_monotonic has not yet been implemented for all multivariate expressions."))
  x <- if (!is.null(x)) as_sym(x) else if (length(fs)) fs[[1L]] else Symbol("x")
  sings <- tryCatch(singularities(f, x, interval), error = function(e) NULL)
  if (!is.null(sings) && isTRUE(is_subset(interval, Reals)) && interval$type %in% c("Interval", "SpecialSet")) {
    iv <- .iv(interval)
    inner <- if (interval$type == "Interval") Interval(iv$a, iv$b, TRUE, TRUE) else interval
    if (!identical(Intersection(inner, sings)$key, S_EmptySet$key)) return(FALSE)
  }
  tp <- solveset(diff(f, x), x, interval)
  identical(Intersection(interval, tp)$key, S_EmptySet$key)
}
#' @rdname calculus-utils
#' @export
is_convex <- function(f, x = NULL, domain = Reals) {
  f <- as_sym(f); if (is.null(x)) x <- .first_free(f)
  x <- as_sym(x)
  .sign_on(.diff1(.diff1(f, x), x), x, domain, "GreaterThan")
}

#' @rdname calculus-utils
#' @export
periodicity <- function(f, x) {
  f <- as_sym(f); x <- as_sym(x)
  # periodicity: a function free of the symbol has period 0
  if (!.has1(f, x)) return(S_Zero)
  p <- tryCatch(.ineq_period(f, x), error = function(e) NULL)
  if (is.null(p)) return(NULL)
  p <- simplify(p)
  # the candidate must leave f unchanged: x + sin(x) has no period
  d <- simplify(.sub(subs(f, x, .add(list(x, p))), f))
  if (is_zero_num(d)) return(p)
  pts <- c(0.37, 1.91, -2.63)
  v <- tryCatch(vapply(pts, function(t) base::Mod(.neval(d, stats::setNames(list(t), x$name))), 0), error = function(e) NA)
  if (all(is.finite(v)) && all(v < 1e-9)) p else NULL
}

#' @rdname calculus-utils
#' @export
euler_equations <- function(L, funcs = NULL, vars = NULL) {
  L <- as_sym(L)
  as_list <- function(v) if (is_sym(v)) list(v) else lapply(v, as_sym)
  funcs <- if (is.null(funcs)) Filter(.is_undefined_application, atoms(L, "Function")) else as_list(funcs)
  if (!length(funcs)) stop("no function found in L")
  vars <- if (is.null(vars)) funcs[[1L]]$args else as_list(vars)
  fkeys <- vapply(funcs, `[[`, "", "key")
  ders <- Filter(function(d) d$args[[1L]]$key %in% fkeys, atoms(L, "Derivative"))
  order <- max(c(0, vapply(ders, function(d) sum(vapply(d$args[-1L], function(t) as.numeric(t$args[[2L]]), 0)), 0)))
  eqs <- lapply(funcs, function(f) {
    eq <- .diff1(L, f)
    for (i in seq_len(order)) for (p in combinations_with_replacement(vars, i)) {
      t <- .diff1(L, Reduce(.diff1, p, f))
      for (v in p) t <- .diff1(t, v)
      eq <- .add(list(eq, if (i %% 2L) .neg(t) else t))
    }
    Eq(eq, S_Zero)
  })
  Filter(function(e) e$type == "Equality", eqs)
}

#' @rdname calculus-utils
#' @export
finite_diff_weights <- function(order, x_list, x0 = 1) {
  xs <- lapply(x_list, as_sym)
  x0 <- as_sym(x0)
  n <- length(xs)
  M <- order
  delta <- array(list(S_Zero), c(M + 1L, n, n))
  delta[[1L, 1L, 1L]] <- S_One
  c1 <- S_One
  for (nn in seq_len(n - 1L)) {
    c2 <- S_One
    for (nu in 0:(nn - 1L)) {
      c3 <- .sub(xs[[nn + 1L]], xs[[nu + 1L]])
      c2 <- .mul(list(c2, c3))
      for (m in 0:min(nn, M)) {
        prev <- delta[[m + 1L, nn, nu + 1L]]
        lower <- if (m > 0L) delta[[m, nn, nu + 1L]] else S_Zero
        delta[[m + 1L, nn + 1L, nu + 1L]] <- cancel(.div(.sub(.mul(list(.sub(xs[[nn + 1L]], x0), prev)),
                                                              .mul(list(.int(m), lower))), c3))
      }
    }
    for (m in 0:min(nn, M)) {
      prev <- delta[[m + 1L, nn, nn]]
      lower <- if (m > 0L) delta[[m, nn, nn]] else S_Zero
      delta[[m + 1L, nn + 1L, nn + 1L]] <- cancel(.mul(list(.div(c1, c2), .sub(.mul(list(.int(m), lower)),
                                                                              .mul(list(.sub(xs[[nn]], x0), prev))))))
    }
    c1 <- c2
  }
  lapply(0:M, function(m) lapply(seq_len(n), function(k) lapply(seq_len(n), function(j) delta[[m + 1L, k, j]])))
}

#' @rdname calculus-utils
#' @export
apply_finite_diff <- function(order, x_list, y_list, x0 = 0) {
  w <- finite_diff_weights(order, x_list, x0)[[order + 1L]][[length(x_list)]]
  # SymPy sums weight times value without simplifying the sum
  .add(Map(function(a, b) .mul(list(simplify(a), as_sym(b))), w, y_list))
}

# ---------------------------------------------------------------------------
# Solver helpers

#' Solver helpers
#'
#' `solve_linear_system(M, syms, ...)` solves an augmented matrix; further
#' symbols may follow `syms`. It returns `NULL` when there is no solution.
#' `solve_linear(lhs, rhs, symbols, exclude)` solves one equation for one of
#' `symbols` when it is linear in it, and returns the pair `(symbol, value)`,
#' `(0, 1)` when nothing is free, or `(numerator, denominator)` otherwise.
#' `checksol(f, symbol, sol)` tests a solution by substitution and returns
#' `TRUE`, `FALSE` or `NA` when undecided.
#' `solve_undetermined_coeffs(eq, coeffs, ...)` makes `eq` vanish identically
#' in the symbols `...`, as SymPy: a named list for one solution, a list of
#' them for several, and `NULL` when no solution exists.
#' `solve_poly_system(eqs, ...)` solves a zero dimensional polynomial system
#' in the generators `...` and returns a list of solution tuples; `det_quick()` and `inv_quick()` are
#' [sym_det()] and [inverse()].
#'
#' @param eqs A list of equations or expressions.
#' @param syms,coeffs Lists of symbols.
#' @param ... Further symbols for `solve_linear_system()`,
#'   `solve_undetermined_coeffs()`, or flags such as `numerical = FALSE` for
#'   `checksol()`.
#' @param dict,set Output forms of `solve_undetermined_coeffs()`, as for
#'   [solve()].
#' @param M An augmented matrix.
#' @param lhs,rhs Sides of an equation.
#' @param x A variable; for `solve_linear()` an alias of `symbols`.
#' @param symbols,exclude Symbols to solve for and symbols to skip.
#' @param strict Kept for SymPy compatibility; roots that cannot be found are
#'   left out.
#' @param symbol A symbol, or a named list (or list of `list(symbol, value)`
#'   pairs) giving the whole solution.
#' @param f,eq An equation or expression.
#' @param sol The value of `symbol`.
#' @return A list, expression, logical or set.
#' @examples
#' x <- Symbol("x"); y <- Symbol("y")
#' checksol(x^2 - 4, list(x = 2))
#' @name solver-helpers
NULL

#' @rdname solver-helpers
#' @export
solve_linear_system <- function(M, syms, ...) {
  M <- Matrix(M)
  syms <- .as_gens(c(list(syms), list(...)))
  n <- length(syms)
  if (M$ncol != n + 1L) stop("the matrix needs one column per symbol plus one", call. = FALSE)
  # solve_lin_sys on system * Matrix(symbols + (-1,)), as SymPy
  eqs <- lapply(.mrows(M), function(r) .add(c(Map(function(c, s) .mul(list(c, s)), r[seq_len(n)], syms),
                                              list(.neg(r[[n + 1L]])))))
  sol <- .sp_solve_lin_sys(lapply(eqs, function(e) .sp_as_poly(e, syms)), syms)
  if (is.null(sol)) return(NULL)
  keep <- !vapply(seq_along(sol$k), function(i) identical(sol$k[[i]]$key, sol$v[[i]]$key), TRUE)
  .sd_to_r(list(k = sol$k[keep], v = sol$v[keep]))
}

#' @rdname solver-helpers
#' @export
solve_linear <- function(lhs, rhs = 0, symbols = list(), exclude = list(), x = NULL) {
  if (!is.null(x)) symbols <- x
  if (is_sym(symbols)) symbols <- list(symbols)
  if (is_sym(exclude)) exclude <- list(exclude)
  .pytuple(.sp_solve_linear(lhs, rhs, lapply(symbols, as_sym), lapply(exclude, as_sym)))
}

#' @rdname solver-helpers
#' @export
checksol <- function(f, symbol, sol = NULL, ...) {
  flags <- list(...)
  d <- .sd_new()
  if (!is.null(sol)) {
    d <- .sd_set(d, as_sym(symbol), as_sym(sol))
  } else if (is.list(symbol) && !is_sym(symbol)) {
    nm <- names(symbol)
    if (length(symbol) && !is.null(nm) && all(nzchar(nm))) {
      # a named list keyed by symbol names: the matching free symbols of f
      fs <- unique_syms(lapply(if (is.list(f) && !is_sym(f)) f else list(f), function(z) free_symbols(as_sym(z))))
      for (k in nm) {
        hit <- Filter(function(z) identical(z$name, k), fs)
        d <- .sd_set(d, if (length(hit)) hit[[1L]] else Symbol(k), as_sym(symbol[[k]]))
      }
    } else for (p in symbol) d <- .sd_set(d, as_sym(p[[1L]]), as_sym(p[[2L]]))
  } else {
    stop(.sympy_condition("ValueError", paste0("Expecting (sym, val) or ({sym: val}, None) but got (",
                                               .str(as_sym(symbol)), ", None)")))
  }
  .sp_checksol(f, d, flags)
}

#' @rdname solver-helpers
#' @export
solve_undetermined_coeffs <- function(eq, coeffs, ..., dict = FALSE, set = FALSE) {
  flags <- list()
  if (isTRUE(dict)) flags$dict <- TRUE
  if (isTRUE(set)) flags$set <- TRUE
  coeffs <- if (is_sym(coeffs)) list(coeffs) else lapply(coeffs, as_sym)
  syms <- lapply(list(...), function(s) if (is_sym(s)) s else lapply(s, as_sym))
  .sp_solve_undetermined_coeffs(as_sym(eq), coeffs, syms, flags)
}

# The coefficients of eq in x set to zero and solved linearly for coeffs: a
# named list, or NULL when inconsistent (used by the recurrence solvers).
.solve_undetermined_linear <- function(eq, coeffs, x) {
  e <- expand(.eq_to_expr(eq))
  x <- as_sym(x)
  coeffs <- .as_gens(list(coeffs))
  .solve_linear_system(.coeff_list(e, x), coeffs)
}

#' @rdname solver-helpers
#' @export
solve_poly_system <- function(eqs, ..., strict = FALSE) {
  if (!is.list(eqs) || is_sym(eqs) || inherits(eqs, "Poly")) eqs <- list(eqs)
  gens <- list(...)
  gens <- if (length(gens)) .as_gens(gens) else list()
  exprs <- lapply(eqs, function(e) if (inherits(e, "Poly")) as_expr(e) else as_sym(e))
  if (!length(gens)) {
    gens <- if (all(vapply(eqs, inherits, TRUE, "Poly")))
      unlist(lapply(eqs, function(e) e$gens[seq_len(.poly_ngens(e))]), recursive = FALSE)
    else unique_syms(lapply(exprs, free_symbols))
    gens <- .sp_sorted(.sp_uniq(gens))
  }
  polys <- lapply(exprs, function(e) .sp_as_poly(e, gens))
  if (any(vapply(polys, is.null, TRUE)))
    stop(.sympy_condition("ComputationFailed", paste0("solve_poly_system(", length(eqs), ") failed without generators")))
  res <- .sp_solve_poly_system(polys, gens)
  if (is.null(res)) NULL else lapply(res, .pytuple)
}
#' @rdname solver-helpers
#' @export
det_quick <- function(M) det(Matrix(M))
#' @rdname solver-helpers
#' @export
inv_quick <- function(M) inverse(Matrix(M))

# ---------------------------------------------------------------------------
# Matrix constructions

#' Matrix constructions
#'
#' `rot_axis1()`, `rot_axis2()` and `rot_axis3()` are the rotation matrices
#' used by SymPy (clockwise about the axis when looking down it);
#' `rot_ccw_axis1()` and friends rotate counterclockwise. `jordan_cell(e, n)`
#' builds a Jordan block, `casoratian(seqs, n)` the Casoratian determinant,
#' `randMatrix(r, c, min, max)` a random integer matrix, `symarray(prefix,
#' shape)` a matrix of symbols, `permanent()` (alias `per()`) the permanent,
#' and `matrix_multiply_elementwise()` the Hadamard product.
#' `ImmutableMatrix()`, `MutableMatrix()`, `MutableDenseMatrix()` and
#' `ImmutableDenseMatrix()` are [Matrix()].
#'
#' @param theta Angle.
#' @param e Diagonal value.
#' @param n Size, or the index variable for `casoratian()`.
#' @param seqs A list of sequences in `n`.
#' @param zero Evaluate at `n = 0, 1, ...` (SymPy's default) instead of `n, n + 1, ...`.
#' @param r,c Dimensions.
#' @param min,max Entry range.
#' @param prefix Symbol prefix.
#' @param shape Dimensions.
#' @param M,A,B Matrices.
#' @param ... Arguments for [Matrix()].
#' @return A matrix or expression.
#' @examples
#' rot_axis3(sym_pi / 2)
#' jordan_cell(Symbol("a"), 3)
#' permanent(Matrix(list(list(1, 2), list(3, 4))))
#' @name matrix-constructions
NULL

#' @rdname matrix-constructions
#' @export
rot_axis1 <- function(theta) {
  ct <- fcall("cos", list(as_sym(theta))); st <- fcall("sin", list(as_sym(theta)))
  Matrix(list(list(1, 0, 0), list(0, ct, st), list(0, .neg(st), ct)))
}
#' @rdname matrix-constructions
#' @export
rot_axis2 <- function(theta) {
  ct <- fcall("cos", list(as_sym(theta))); st <- fcall("sin", list(as_sym(theta)))
  Matrix(list(list(ct, 0, .neg(st)), list(0, 1, 0), list(st, 0, ct)))
}
#' @rdname matrix-constructions
#' @export
rot_axis3 <- function(theta) {
  ct <- fcall("cos", list(as_sym(theta))); st <- fcall("sin", list(as_sym(theta)))
  Matrix(list(list(ct, st, 0), list(.neg(st), ct, 0), list(0, 0, 1)))
}
#' @rdname matrix-constructions
#' @export
rot_ccw_axis1 <- function(theta) t(rot_axis1(theta))
#' @rdname matrix-constructions
#' @export
rot_ccw_axis2 <- function(theta) t(rot_axis2(theta))
#' @rdname matrix-constructions
#' @export
rot_ccw_axis3 <- function(theta) t(rot_axis3(theta))
#' @rdname matrix-constructions
#' @export
jordan_cell <- function(e, n) Matrix(n, n, function(i, j) if (i == j) e else if (j == i + 1L) 1 else 0)
#' @rdname matrix-constructions
#' @export
casoratian <- function(seqs, n, zero = TRUE) {
  n <- as_sym(n)
  seqs <- lapply(seqs, as_sym)
  k <- length(seqs)
  at <- function(i) if (zero) .int(i) else .add(list(n, .int(i)))
  if (!k) return(S_One)
  simplify(det(Matrix(lapply(seq_len(k) - 1L, function(i) lapply(seqs, function(s) subs(s, n, at(i)))))))
}
#' @rdname matrix-constructions
#' @export
randMatrix <- function(r, c = r, min = 0, max = 99) {
  pool <- min:max
  Matrix(r, c, function(i, j) pool[sample.int(length(pool), 1L)])
}
#' @rdname matrix-constructions
#' @export
symarray <- function(prefix, shape) {
  shape <- as.integer(shape)
  if (length(shape) == 1L) return(Matrix(lapply(seq_len(shape) - 1L, function(i) Symbol(paste0(prefix, "_", i)))))
  if (length(shape) == 2L)
    return(Matrix(shape[1L], shape[2L], function(i, j) Symbol(paste0(prefix, "_", i - 1L, "_", j - 1L))))
  # higher ranks: an Array in row-major order
  idx <- expand.grid(rev(lapply(shape, function(n) seq_len(n) - 1L)))[, rev(seq_along(shape)), drop = FALSE]
  Array(lapply(seq_len(nrow(idx)), function(r) Symbol(paste0(prefix, "_", paste(unlist(idx[r, ]), collapse = "_")))), shape)
}
#' @rdname matrix-constructions
#' @export
permanent <- function(M) {
  M <- Matrix(M)
  if (M$nrow > M$ncol) M <- t(M)
  n <- M$nrow
  if (M$ncol > 9L) stop("permanent: matrix too large")
  # rectangular: sum over injections of the rows into the columns
  .add(lapply(permutations(seq_len(M$ncol), n), function(p) .mul(lapply(seq_len(n), function(i) .mget(M, i, p[i])))))
}
#' @rdname matrix-constructions
#' @export
per <- permanent
#' @rdname matrix-constructions
#' @export
matrix_multiply_elementwise <- function(A, B) hadamard_product(A, B)
#' @rdname matrix-constructions
#' @export
ImmutableMatrix <- function(...) Matrix(...)
#' @rdname matrix-constructions
#' @export
MutableMatrix <- ImmutableMatrix
#' @rdname matrix-constructions
#' @export
MutableDenseMatrix <- ImmutableMatrix
#' @rdname matrix-constructions
#' @export
ImmutableDenseMatrix <- ImmutableMatrix

# ---------------------------------------------------------------------------
# Arrays, indexed objects, discrete transforms

#' Array aliases and indexed objects
#'
#' `NDimArray()`, `MutableDenseNDimArray()`, `ImmutableDenseNDimArray()`
#' and `DenseNDimArray()` are [Array()], and `shape()` its dimensions.
#' Indexed objects are documented in [Indexed()].
#'
#' @param ... Arguments for [Array()].
#' @param a An array.
#' @return An array, function or `sym` expression.
#' @examples
#' shape(Array(list(1, 2, 3)))
#' @name array-aliases
NULL

#' @rdname array-aliases
#' @export
NDimArray <- function(...) Array(...)
#' @rdname array-aliases
#' @export
MutableDenseNDimArray <- NDimArray
#' @rdname array-aliases
#' @export
ImmutableDenseNDimArray <- NDimArray
#' @rdname array-aliases
#' @export
DenseNDimArray <- NDimArray
#' @rdname array-aliases
#' @export
SparseNDimArray <- NDimArray
#' @rdname array-aliases
#' @export
MutableSparseNDimArray <- NDimArray
#' @rdname array-aliases
#' @export
ImmutableSparseNDimArray <- NDimArray
#' @rdname array-aliases
#' @export
shape <- function(a) if (inherits(a, "SymArray")) a$shape else dim(a)

#' Discrete transform aliases
#'
#' `ifft()` is [sym_ifft()] for lists of symbolic or exact values. `fwht()`
#' and `ifwht()` are the Walsh-Hadamard transform and its inverse.
#' `covering_product()` and `intersecting_product()` multiply sequences
#' indexed by bitmasks, combining over unions and intersections of the
#' index sets.
#'
#' @param seq,a,b Sequences.
#' @return A list of expressions.
#' @examples
#' covering_product(list(1, 2), list(3, 4))
#' @name discrete-aliases
NULL

#' @noRd
fft <- function(z, inverse = FALSE) {
  if (is.numeric(z) || is.complex(z)) return(stats::fft(z, inverse = inverse))
  if (inverse) sym_ifft(z) else sym_fft(z)
}
#' @rdname discrete-aliases
#' @export
ifft <- function(seq) sym_ifft(seq)
#' @rdname discrete-aliases
#' @export
fwht <- function(seq) walsh_hadamard_transform(seq)
#' @rdname discrete-aliases
#' @export
ifwht <- function(seq) walsh_hadamard_transform(seq, inverse = TRUE)
#' @rdname discrete-aliases
#' @export
covering_product <- function(a, b) {
  if (!length(a) || !length(b)) return(list())
  a <- .dc_seq(a); b <- .dc_seq(b)
  n <- max(length(a), length(b))
  A <- mobius_transform(c(as.list(a), rep(list(0), n - length(a))))
  B <- mobius_transform(c(as.list(b), rep(list(0), n - length(b))))
  lapply(inverse_mobius_transform(Map(function(u, v) .mul(list(u, v)), A, B)), expand)
}
#' @rdname discrete-aliases
#' @export
intersecting_product <- function(a, b) {
  if (!length(a) || !length(b)) return(list())
  a <- .dc_seq(a); b <- .dc_seq(b)
  n <- max(length(a), length(b))
  A <- mobius_transform(c(as.list(a), rep(list(0), n - length(a))), subset = FALSE)
  B <- mobius_transform(c(as.list(b), rep(list(0), n - length(b))), subset = FALSE)
  lapply(inverse_mobius_transform(Map(function(u, v) .mul(list(u, v)), A, B), subset = FALSE), expand)
}

# ---------------------------------------------------------------------------
# Sets, geometry, printing, utilities

#' Set aliases
#'
#' `imageset()` is the simplified image of sets under a map, as SymPy's
#' `imageset()`: the map is a `Lambda`, a function name such as `"cos"`, or
#' a symbol (or list of symbols) followed by an expression, and the sets
#' come after it. `Contains(x, s)` is the membership condition, and
#' `PowerSet(s)` the set of all subsets of a finite set.
#'
#' @param ... The map, then the sets.
#' @param x,s An element and a set.
#' @return A set or boolean.
#' @examples
#' x <- Symbol("x")
#' imageset(x, 2 * x, Interval(0, 1))
#' imageset(Lambda(x, x^2), Interval(-2, 1))
#' @export
imageset <- function(...) {
  # sympy.sets.sets.imageset
  args <- list(...)
  if (length(args) < 2L) .sset_err("ValueError", sprintf("imageset expects at least 2 args, got: %d", length(args)))
  a1 <- args[[1L]]
  if ((is_sym(a1) && a1$type == "Symbol") || (is.list(a1) && !is_sym(a1)) || (is_sym(a1) && a1$type == "Tuple")) {
    if (length(args) > 2L) {
      body <- args[[2L]]
      if (is.list(body) && !is_sym(body)) body <- do.call(Tuple, lapply(body, as_sym))
      f <- Lambda(if (is_sym(a1) && a1$type == "Tuple") a1$args else a1, body)
      set_list <- args[-(1:2)]
    } else { f <- a1; set_list <- args[-1L] }
  } else { f <- a1; set_list <- args[-1L] }
  if (is.character(f) && length(f) == 1L) f <- Lambda(Symbol("x"), fcall(f, list(Symbol("x"))))
  else if (is.function(f)) f <- Lambda(Symbol("x"), as_sym(f(Symbol("x"))))
  else if (!(is_sym(f) && f$type == "Lambda"))
    .sset_err("TypeError", .filldedent(sprintf("expecting lambda, Lambda, or FunctionClass, not '%s'.",
                                               if (is_sym(f)) .sympy_class(f) else class(f)[1L])))
  set_list <- lapply(set_list, function(s) if (is.list(s) && !is_sym(s)) s else as_sym(s))
  if (!all(vapply(set_list, function(s) is_sym(s) && is_set(s), NA))) {
    nm <- vapply(set_list, function(s) if (!is_sym(s)) "list" else if (s$type == "SpecialSet") s$name else .sympy_class(s), "")
    .sset_err("ValueError", sprintf("arguments after mapping should be sets, not [%s]", paste0("'", nm, "'", collapse = ", ")))
  }
  if (length(set_list) == 1L) return(.imageset(f, set_list[[1L]]))
  do.call(ImageSet, c(list(f), set_list))
}

#' @rdname imageset
#' @export
Contains <- function(x, s) {
  x <- if (inherits(x, "symmatrix")) .as_block(x) else as_sym(x)
  s <- as_sym(s)
  # Contains.__new__: only a true or false _contains result evaluates
  r <- .set_contains_raw(s, x)
  if (!is.null(r) && r$type == "BooleanAtom") r else new_sym("Contains", list(x, s))
}
#' @rdname imageset
#' @export
PowerSet <- function(s) {
  s <- as_sym(s)
  # PowerSet.__new__ keeps the set unevaluated; rewrite(, "FiniteSet") lists it
  if (!is_set(s)) stop(structure(class = c("ValueError", "error", "condition"),
                                 list(message = paste0(.str(s), " must be a set."), call = NULL)))
  new_sym("PowerSet", list(s))
}
register_str_printer("PowerSet", function(e) paste0("PowerSet(", .str(e$args[[1L]]), ")"))
# PowerSet._eval_rewrite_as_FiniteSet: arg.powerset() of a finite set
.powerset_finite <- function(s) {
  if (s$type == "EmptySet") return(FiniteSet(S_EmptySet))
  if (s$type != "FiniteSet") return(NULL)
  FiniteSet(lapply(subsets(s$args), function(sub) if (length(sub)) FiniteSet(sub) else S_EmptySet))
}

#' Geometry aliases and helpers
#'
#' `Point2D()`, `Line2D()` and `Segment2D()` are the 2D constructors.
#' `rad()` and `deg()` convert angles. `RegularPolygon(center, r, n)` gives
#' the vertices of a regular polygon as a [Polygon()]. `convex_hull()`
#' returns the hull of points (Andrew's monotone chain), and
#' `closest_points()` and `farthest_points()` the pairs at minimal and
#' maximal distance.
#'
#' @param ... Coordinates or points.
#' @param p,q Points.
#' @param x An angle.
#' @param center,r,n Regular polygon data.
#' @return A geometric entity, expression or list.
#' @examples
#' convex_hull(Point(0, 0), Point(1, 1), Point(2, 0), Point(1, 0.5), Point(1, -1))
#' @name geometry-aliases
NULL

#' @rdname geometry-aliases
#' @export
Point2D <- function(...) Point(...)
#' @rdname geometry-aliases
#' @export
Line2D <- function(p, q) Line(p, q)
#' @rdname geometry-aliases
#' @export
Segment2D <- function(p, q) Segment(p, q)
#' @rdname geometry-aliases
#' @export
rad <- function(x) .mul(list(as_sym(x), S_Pi, .rat(1L, 180L)))
#' @rdname geometry-aliases
#' @export
deg <- function(x) .mul(list(as_sym(x), .int(180), .inv(S_Pi)))
#' @rdname geometry-aliases
#' @export
RegularPolygon <- function(center, r, n) {
  c <- center$coords
  do.call(Polygon, lapply(seq_len(n) - 1L, function(k) {
    ang <- .mul(list(.rat(2L * k, n), S_Pi))
    Point(.add(list(c[[1L]], .mul(list(as_sym(r), fcall("cos", list(ang)))))),
          .add(list(c[[2L]], .mul(list(as_sym(r), fcall("sin", list(ang)))))))
  }))
}

# Sign of a - b, exactly when possible.
.sgn_cmp <- function(a, b) {
  r <- .cmp(as_sym(a), as_sym(b))
  if (!is.na(r)) return(r)
  sign(as.numeric(N(.sub(as_sym(a), as_sym(b)))))
}

#' @rdname geometry-aliases
#' @export
convex_hull <- function(...) {
  pts <- list(...)
  pts <- pts[!duplicated(vapply(pts, function(p) paste(vapply(p$coords, .k, ""), collapse = ","), ""))]
  if (length(pts) == 1L) return(pts[[1L]])
  before <- function(p, q) {
    c1 <- .sgn_cmp(p$coords[[1L]], q$coords[[1L]])
    c1 < 0 || (c1 == 0 && .sgn_cmp(p$coords[[2L]], q$coords[[2L]]) < 0)
  }
  # ponytail: insertion sort with exact comparisons; fine for small point sets
  o <- list()
  for (p in pts) {
    k <- length(o)
    while (k > 0 && before(p, o[[k]])) k <- k - 1L
    o <- append(o, list(p), after = k)
  }
  pts <- o
  cross2 <- function(o, a, b) {
    P <- pts[[o]]$coords; A <- pts[[a]]$coords; B <- pts[[b]]$coords
    .sgn_cmp(.mul(list(.sub(A[[1L]], P[[1L]]), .sub(B[[2L]], P[[2L]]))),
             .mul(list(.sub(A[[2L]], P[[2L]]), .sub(B[[1L]], P[[1L]]))))
  }
  lower <- integer(); upper <- integer()
  for (i in seq_along(pts)) {
    while (length(lower) >= 2L && cross2(lower[length(lower) - 1L], lower[length(lower)], i) <= 0) lower <- lower[-length(lower)]
    lower <- c(lower, i)
  }
  for (i in rev(seq_along(pts))) {
    while (length(upper) >= 2L && cross2(upper[length(upper) - 1L], upper[length(upper)], i) <= 0) upper <- upper[-length(upper)]
    upper <- c(upper, i)
  }
  idx <- c(lower[-length(lower)], upper[-length(upper)])
  if (length(idx) == 2L) return(Segment(pts[[idx[1L]]], pts[[idx[2L]]]))
  do.call(Polygon, pts[idx])
}
.pair_search <- function(pts, best) {
  pts <- pts[!duplicated(vapply(pts, function(p) paste(vapply(p$coords, .k, ""), collapse = ","), ""))]
  out <- list()
  bv <- NULL
  for (i in seq_len(length(pts) - 1L)) for (j in (i + 1L):length(pts)) {
    d <- distance(pts[[i]], pts[[j]])
    c1 <- if (is.null(bv)) best else .sgn_cmp(d, bv)
    if (c1 == best) {
      bv <- d
      out <- list(list(pts[[i]], pts[[j]]))
    } else if (c1 == 0) {
      out[[length(out) + 1L]] <- list(pts[[i]], pts[[j]])
    }
  }
  out
}
#' @rdname geometry-aliases
#' @export
closest_points <- function(...) .pair_search(list(...), -1)
#' @rdname geometry-aliases
#' @export
farthest_points <- function(...) .pair_search(list(...), 1)

#' Printing helpers
#'
#' `print_latex()`, `print_ccode()`, `print_fcode()`, `print_rcode()` and
#' `print_jscode()` print the corresponding strings. `julia_code()`,
#' `cxxcode()`, `mathematica_code()` and `maple_code()` are further code
#' printers. `print_tree()` shows the expression tree.
#'
#' @param e An expression.
#' @param assign_to Optional variable name.
#' @param ... Printer settings, as in SymPy (for example `standard`).
#' @return The string, invisibly for the `print_` functions.
#' @examples
#' x <- Symbol("x")
#' mathematica_code(sin(x)^2 + exp(x))
#' print_tree(x^2 + 1)
#' @name printing-helpers
NULL

.code_styles$julia <- list(plain_frac = TRUE, pow = function(b, e) paste0(b, " .^ ", e), sqrt = "sqrt", pi = "pi", e = "e",
                           num = function(q) paste0(gmp::numerator(q), " // ", gmp::denominator(q)), abs = "abs",
                           ln = "log", ceiling = "ceil", I = "im", inf = "Inf")
.code_styles$cxx <- list(pow = function(b, e) paste0("std::pow(", b, ", ", e, ")"), sqrt = "std::sqrt", pi = "M_PI", e = "M_E",
                         num = function(q) paste0(gmp::numerator(q), ".0/", gmp::denominator(q), ".0"), abs = "std::abs",
                         ln = "std::log", ceiling = "std::ceil", I = "std::complex<double>(0, 1)", inf = "INFINITY",
                         prefix = "std::")
.code_styles$maple <- list(pow = function(b, e) paste0(b, "^", e), sqrt = "sqrt", pi = "Pi", e = "exp(1)",
                           num = function(q) paste0(gmp::numerator(q), "/", gmp::denominator(q)), abs = "abs",
                           ln = "ln", ceiling = "ceil", I = "I", inf = "infinity")

#' @rdname printing-helpers
#' @export
julia_code <- .mk_printer("julia", "%s = %s")
#' @rdname printing-helpers
#' @export
cxxcode <- .mk_printer("cxx", "%s = %s;")
#' @rdname printing-helpers
#' @export
maple_code <- .mk_printer("maple", "%s := %s")

.MMA_NAMES <- c(exp = "Exp", log = "Log", sin = "Sin", cos = "Cos", tan = "Tan", cot = "Cot", sec = "Sec",
                csc = "Csc", asin = "ArcSin", acos = "ArcCos", atan = "ArcTan", sinh = "Sinh", cosh = "Cosh",
                tanh = "Tanh", asinh = "ArcSinh", acosh = "ArcCosh", atanh = "ArcTanh", gamma = "Gamma",
                Abs = "Abs", erf = "Erf", zeta = "Zeta", factorial = "Factorial", binomial = "Binomial",
                floor = "Floor", ceiling = "Ceiling", sign = "Sign", besselj = "BesselJ", bessely = "BesselY")

.mma <- function(e) {
  t <- e$type
  paren <- function(z) { s <- .mma(z); if (.prec(z) < PREC$Mul) paste0("(", s, ")") else s }
  switch(t,
         Integer = as.character(gmp::numerator(e$value)),
         Rational = paste0(gmp::numerator(e$value), "/", gmp::denominator(e$value)),
         Float = sprintf("%.15g", e$value),
         Symbol = , Dummy = e$name,
         Constant = switch(EXPR = e$name, pi = "Pi", E = "E", I = "I", oo = "Infinity", `-oo` = "-Infinity",
                           EulerGamma = "EulerGamma", GoldenRatio = "GoldenRatio", Catalan = "Catalan",
                           zoo = "ComplexInfinity", nan = "Indeterminate", e$name),
         Add = {
           terms <- .add_print_order(e$args)
           out <- vapply(seq_along(terms), function(i) {
             z <- terms[[i]]
             if (i > 1L && could_extract_minus_sign(z)) paste0(" - ", paren(.neg(z)))
             else paste0(if (i > 1L) " + ", .mma(z))
           }, "")
           paste(out, collapse = "")
         },
         Mul = {
           if (could_extract_minus_sign(e)) return(paste0("-", paren(.neg(e))))
           paste(vapply(.mul_print_order(e$args), paren, ""), collapse = "*")
         },
         Pow = {
           b <- e$args[[1L]]; x <- e$args[[2L]]
           bs <- .mma(b); if (.prec(b) < PREC$Atom) bs <- paste0("(", bs, ")")
           xs <- .mma(x); if (.prec(x) < PREC$Atom || x$type == "Rational") xs <- paste0("(", xs, ")")
           paste0(bs, "^", xs)
         },
         Function = {
           nm <- .MMA_NAMES[e$name]
           if (is.na(nm)) nm <- paste0(toupper(substr(e$name, 1L, 1L)), substring(e$name, 2L))
           paste0(nm, "[", paste(vapply(e$args, .mma, ""), collapse = ", "), "]")
         },
         Integral = {
           lims <- e$args[-1L]
           args <- if (length(lims) == 1L && length(lims[[1L]]$args) == 1L) list(e$args[[1L]], lims[[1L]]$args[[1L]]) else e$args
           paste0("Hold[Integrate[", paste(vapply(args, .mma, ""), collapse = ", "), "]]")
         },
         Sum = paste0("Hold[Sum[", paste(vapply(e$args, .mma, ""), collapse = ", "), "]]"),
         Tuple = paste0("{", paste(vapply(e$args, .mma, ""), collapse = ", "), "}"),
         stop("mathematica_code: unsupported ", t))
}

#' @rdname printing-helpers
#' @export
mathematica_code <- function(e) .mma(as_sym(e))

#' @rdname printing-helpers
#' @export
print_latex <- function(e) { s <- latex(e); cat(s, "\n"); invisible(s) }
#' @rdname printing-helpers
#' @export
print_ccode <- function(e, assign_to = NULL) { s <- ccode(e, assign_to); cat(s, "\n"); invisible(s) }
#' @rdname printing-helpers
#' @export
print_fcode <- function(e, assign_to = NULL) { s <- fcode(e, assign_to); cat(s, "\n"); invisible(s) }
#' @rdname printing-helpers
#' @export
print_rcode <- function(e, assign_to = NULL) { s <- rcode(e, assign_to); cat(s, "\n"); invisible(s) }
#' @rdname printing-helpers
#' @export
print_jscode <- function(e, assign_to = NULL) { s <- jscode(e, assign_to); cat(s, "\n"); invisible(s) }
#' @rdname printing-helpers
#' @export
print_tree <- function(e) {
  lines <- character()
  walk <- function(z, prefix, last) {
    head <- if (length(z$args)) z$type else paste0(z$type, ": ", .str(z))
    if (z$type == "Function") head <- paste0(z$name, "()")
    lines <<- c(lines, paste0(prefix, if (nchar(prefix)) (if (last) "+-" else "|-") else "", head))
    np <- paste0(prefix, if (nchar(prefix)) (if (last) "  " else "| ") else "  ")
    for (i in seq_along(z$args)) walk(z$args[[i]], np, i == length(z$args))
  }
  walk(as_sym(e), "", TRUE)
  cat(lines, sep = "\n")
  invisible(lines)
}

#' Text plot
#'
#' `textplot(expr, a, b, W, H)` prints an ASCII plot of a function of one
#' variable on `[a, b]`, with the layout of SymPy's `textplot()`.
#'
#' @param expr An expression in one variable.
#' @param a,b Range; numbers or numeric expressions such as `2*pi`.
#' @param W,H Width and height in characters.
#' @return The character rows, invisibly.
#' @examples
#' textplot(Symbol("x")^2, -2, 2, W = 20, H = 6)
#' @export
textplot <- function(expr, a, b, W = 55, H = 21) {
  # sympy.plotting.textplot.textplot_str
  e <- as_sym(expr)
  fs <- free_symbols(e)
  if (length(fs) > 1L) stop(.sympy_condition("ValueError", sprintf("The expression must have a single variable. (Got {%s})",
                                                                   paste(vapply(fs, .str, ""), collapse = ", "))))
  W <- as.integer(W); H <- as.integer(H)
  a <- as.numeric(N(as_sym(a))); b <- as.numeric(N(as_sym(b)))
  xs <- vapply(0:(W - 1L), function(k) a + (b - a) * k / (W - 1L), 0)
  ys <- if (length(fs)) .eval_grid(e, list(fs[[1L]]), list(xs)) else rep(as.numeric(N(e)), W)
  ys <- suppressWarnings(as.numeric(ys))
  ok <- is.finite(ys)
  if (any(ok)) {
    ma <- max(ys[ok]); mi <- min(ys[ok])
    if (ma == mi) {
      if (ma != 0) { r <- sort(c(0, 2 * ma)); mi <- r[1L]; ma <- r[2L] } else { mi <- -1; ma <- 1 }
    }
  } else {
    mi <- -1; ma <- 1
  }
  prec <- -(floor(log10(ma - mi)) - 1)
  mi <- round(mi, prec); ma <- round(ma, prec)
  norm <- ma - mi; off <- (ma + mi) / 2
  yr <- rep(NA_integer_, W)
  for (i in which(ok)) {
    nv <- (ys[i] - off) / norm
    # Float(v).round(): half to even on the 15 significant digit decimal value
    if (is.finite(nv)) yr[i] <- as.integer(round(signif((nv * H + H / 2) * (H - 1L) / H, 15L)))
  }
  ybins <- vapply(0:(H - 1L), function(k) mi + (ma - mi) * k / (H - 1L), 0)
  eqh <- function(i, h) !is.na(yr[i]) && yr[i] == h
  margin <- 7L
  out <- character()
  for (h in (H - 1L):0L) {
    s <- rep(" ", W)
    for (i in seq_len(W)) {
      if (eqh(i, h)) {
        s[i] <- if ((i == 1L || eqh(i - 1L, h - 1L)) && (i == W || eqh(i + 1L, h + 1L))) "/"
          else if ((i == 1L || eqh(i - 1L, h + 1L)) && (i == W || eqh(i + 1L, h - 1L))) "\\"
          else "."
      }
    }
    if (h == 0L) s <- rep("_", W)
    prefix <- if (h %in% c(0L, H %/% 2L, H - 1L)) substr(sprintf("%7s", sprintf("%g", ybins[h + 1L])), 1L, margin) else strrep(" ", margin)
    s <- paste(s, collapse = "")
    if (h == H %/% 2L) s <- gsub(" ", "-", s, fixed = TRUE)
    out <- c(out, paste0(prefix, " |", s))
  }
  half <- W %/% 2L
  bottom <- paste0(strrep(" ", margin + 2L), sprintf("%-*s", half, sprintf("%g", xs[1L])),
                   sprintf("%-*s", if (W %% 2L == 1L) half else half - 1L, sprintf("%g", xs[half + 1L])),
                   sprintf("%g", xs[W]))
  out <- c(out, bottom)
  cat(out, sep = "\n")
  invisible(out)
}

#' Iteration utilities
#'
#' Small helpers from `sympy.utilities.iterables`: `flatten()`, `unflatten()`,
#' `group()` (runs of equal elements), `take()`, `variations()`, `sift()`,
#' `has_dups()`, `has_variety()`, `rotations()`, `prefixes()`, `postfixes()`,
#' `numbered_symbols()`, `dict_merge()` and `topological_sort()`.
#'
#' @param repetition Allow repeated elements.
#' @param f A key function.
#' @param prefix Symbol prefix.
#' @param start First index.
#' @param ... Lists to merge.
#' @param vertices,edges A graph, with edges as pairs.
#' @param key For `topological_sort()`, a function giving the sort key of a
#'   vertex (the vertex itself by default).
#' @param multiple For `group()`, return runs instead of value and count.
#' @return A list or logical.
#' @examples
#' flatten(list(1, list(2, list(3))))
#' group(c(1, 1, 2, 3, 3, 3))
#' topological_sort(1:4, list(c(1, 2), c(2, 3), c(1, 4)))
#' @name iterables
NULL

#' @rdname iterables
#' @export
flatten <- function(x) {
  out <- list()
  walk <- function(z) {
    if (is.list(z) && !inherits(z, "sym")) for (u in z) walk(u)
    else if (!inherits(z, "sym") && length(z) > 1L) for (i in seq_along(z)) walk(z[[i]])
    else out[[length(out) + 1L]] <<- z
  }
  walk(x)
  out
}
#' @rdname iterables
#' @export
unflatten <- function(x, n = 2) {
  if (!length(x)) return(list())
  lapply(seq(1, length(x), by = n), function(i) x[i:(i + n - 1L)])
}
#' @rdname iterables
#' @export
group <- function(x, multiple = TRUE) {
  key <- function(v) if (is_sym(v)) v$key else paste(deparse(v), collapse = "")
  out <- list()
  for (i in seq_along(x)) {
    v <- x[[i]]
    n <- length(out)
    if (n && identical(key(out[[n]]$v), key(v))) out[[n]]$l <- out[[n]]$l + 1L
    else out[[n + 1L]] <- list(v = v, l = 1L)
  }
  if (multiple) return(lapply(out, function(g) if (is.list(x)) rep(list(g$v), g$l) else rep(g$v, g$l)))
  lapply(out, function(g) list(g$v, g$l))
}
#' @rdname iterables
#' @export
take <- function(x, n) {
  # list(islice(iter, n)) on a generator such as series(n = NULL)
  if (inherits(x, "sym_generator")) {
    out <- list()
    while (length(out) < n && !is.null(v <- x$nxt())) out[[length(out) + 1L]] <- v
    return(out)
  }
  x[seq_len(min(n, length(x)))]
}
#' @rdname iterables
#' @export
variations <- function(x, n, repetition = FALSE) {
  if (!repetition) return(permutations(x, n))
  if (n == 0L) return(list(x[0]))
  g <- expand.grid(rep(list(seq_along(x)), n))[, rev(seq_len(n)), drop = FALSE]
  lapply(seq_len(nrow(g)), function(i) x[unlist(g[i, ])])
}
#' @rdname iterables
#' @export
sift <- function(x, f) split(x, vapply(x, function(z) as.character(f(z)), ""))
#' @rdname iterables
#' @export
has_dups <- function(x) anyDuplicated(if (is.list(x)) vapply(x, function(z) if (is_sym(z)) z$key else paste(z), "") else x) > 0L
#' @rdname iterables
#' @export
has_variety <- function(x) length(unique(x)) > 1L
#' @rdname iterables
#' @export
rotations <- function(x) lapply(seq_along(x) - 1L, function(k) x[(seq_along(x) + k - 1L) %% length(x) + 1L])
#' @rdname iterables
#' @export
prefixes <- function(x) lapply(seq_along(x), function(k) x[seq_len(k)])
#' @rdname iterables
#' @export
postfixes <- function(x) lapply(rev(seq_along(x)), function(k) x[k:length(x)])
#' @rdname iterables
#' @export
numbered_symbols <- function(prefix = "x", n = 10, start = 0) lapply(start + seq_len(n) - 1L, function(i) Symbol(paste0(prefix, i)))
#' @rdname iterables
#' @export
dict_merge <- function(...) {
  out <- list()
  for (d in list(...)) out[names(d)] <- d
  out
}
#' @rdname iterables
#' @export
topological_sort <- function(vertices, edges, key = NULL) {
  # sympy.utilities.iterables.topological_sort: S kept in descending key
  # order and popped from its end, freed vertices inserted in place
  if (is.null(key)) key <- function(z) z
  gt <- function(a, b) {
    a <- key(a); b <- key(b)
    if (is.character(a) && is.character(b)) {
      x <- utf8ToInt(a); y <- utf8ToInt(b); n <- min(length(x), length(y))
      d <- which(x[seq_len(n)] != y[seq_len(n)])
      return(if (length(d)) x[d[1L]] > y[d[1L]] else length(x) > length(y))
    }
    if (is_sym(a) || is_sym(b)) return(as.numeric(num_value_dbl(as_sym(a))) > as.numeric(num_value_dbl(as_sym(b))))
    a > b
  }
  same <- function(a, b) identical(paste(a), paste(b))
  E <- lapply(edges, function(e) as.list(e))
  V <- as.list(vertices)
  S <- list()
  for (z in V) if (!any(vapply(E, function(e) same(e[[2L]], z), NA)) && !any(vapply(S, same, NA, z))) S[[length(S) + 1L]] <- z
  # sorted(S, key=key, reverse=True): insertion keeps the sort stable
  srt <- list()
  for (z in S) {
    i <- 1L
    while (i <= length(srt) && !gt(z, srt[[i]])) i <- i + 1L
    srt <- append(srt, list(z), after = i - 1L)
  }
  S <- srt
  L <- list()
  while (length(S)) {
    node <- S[[length(S)]]
    S <- S[-length(S)]
    L[[length(L) + 1L]] <- node
    for (e in E) {
      if (!same(e[[1L]], node)) next
      k <- which(vapply(E, function(f) same(f[[1L]], e[[1L]]) && same(f[[2L]], e[[2L]]), NA))[1L]
      E <- E[-k]
      v <- e[[2L]]
      if (any(vapply(E, function(f) same(f[[2L]], v), NA))) next
      pos <- which(vapply(S, function(s) gt(v, s), NA))
      S <- if (length(pos)) append(S, list(v), after = pos[1L] - 1L) else c(S, list(v))
    }
  }
  if (length(E)) stop(structure(class = c("ValueError", "error", "condition"), list(message = "cycle detected", call = NULL)))
  if (all(vapply(L, function(z) is.atomic(z) && length(z) == 1L, NA))) unlist(L) else L
}
