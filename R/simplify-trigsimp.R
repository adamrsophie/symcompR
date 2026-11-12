# trigsimp after sympy/simplify/trigsimp.py: futrig (the default
# "matching" method, a greedy search over the fu rules), "fu", "groebner",
# "combined" and "old", and hyper_as_trig of sympy/simplify/fu.py.

.TRIG_FUNS <- c("sin", "cos", "tan", "cot", "sec", "csc")
.HYP_FUNS <- c("sinh", "cosh", "tanh", "coth", "sech", "csch")

# SymPy's simplify.bottom_up: F on every node with arguments, after the
# arguments; atoms are left alone.
.bottom_up_args <- function(e, F) {
  if (!is_sym(e) || !length(e$args)) return(e)
  args <- lapply(e$args, function(a) .bottom_up_args(a, F))
  if (!identical(vapply(args, .k, ""), vapply(e$args, .k, ""))) e <- sym_rebuild(e, args)
  F(e)
}

# _eapply: func on e when e is an expression (all of ours are)
.eapply <- function(func, e, cond = NULL) {
  if (!is_sym(e)) return(e)
  if (is_relational(e) || e$type %in% c("And", "Or", "Not", "Xor", "Implies", "Equivalent"))
    return(sym_rebuild(e, lapply(e$args, function(a) if (is.null(cond) || cond(a)) .eapply(func, a) else a)))
  func(e)
}

.as_independent_funs <- function(e, names) {
  fs <- if (e$type == "Mul") e$args else list(e)
  dep <- vapply(fs, function(f) .has_function(f, names), TRUE)
  list(.mul(fs[!dep]), .mul(fs[dep]))
}

# _split_add: the terms of a sum that take part in trigonometric
# simplification, and the rest
.split_add_funs <- function(e, names) {
  terms <- Filter(function(a) .has_function(a, names), e$args)
  coeffs <- lapply(terms, function(a) .as_independent_funs(a, names)[[1L]])
  fsyms <- unique(unlist(lapply(terms, function(a) lapply(Filter(function(f) f$name %in% names, atoms(a, "Function")),
                                                           function(f) vapply(free_symbols(f$args[[1L]]), .k, "")))))
  related_exp <- function(a) {
    for (f in Filter(function(f) f$name == "exp", atoms(a, "Function"))) {
      arg <- f$args[[1L]]
      if (length(intersect(vapply(free_symbols(arg), .k, ""), fsyms)) &&
          (identical(names, .HYP_FUNS) || .contains_key(arg, S_I$key))) return(TRUE)
    }
    FALSE
  }
  ind <- list(); act <- list()
  for (a in e$args) {
    isnum <- !length(free_symbols(a))
    if (.has_function(a, names) || isnum || related_exp(a) ||
        any(vapply(coeffs, function(cc) !length(free_symbols(.div(a, cc))), TRUE))) act[[length(act) + 1L]] <- a
    else ind[[length(ind) + 1L]] <- a
  }
  list(.add(ind), .add(act))
}

.expr_normal <- function(e) {
  nd <- .as_numer_denom(e)
  if (is_one_num(nd[[2L]])) return(nd[[1L]])
  # Expr.normal: _unevaluated_Mul(n, 1/d) for a numeric d
  if (is_number_type(nd[[2L]])) return(.umul(nd[[1L]], .inv(nd[[2L]])))
  .div(nd[[1L]], nd[[2L]])
}

.mexpand <- function(e) expand_mul(expand_multinomial(e))

# _TR11: double angles in the numerator against single angles in the
# denominator
.TR11_ <- function(rv) .bottom_up_args(rv, function(rv) {
  sincos_args <- function(flat) {
    out <- list(sin = list(), cos = list())
    for (fi in (if (flat$type == "Mul") flat$args else list(flat))) {
      be <- .base_exp(fi)
      if (be[[2L]]$type == "Integer" && num_sign(be[[2L]]) > 0 && be[[1L]]$type == "Function" && be[[1L]]$name %in% c("sin", "cos")) {
        a <- be[[1L]]$args[[1L]]
        out[[be[[1L]]$name]][[a$key]] <- a
      }
    }
    out
  }
  nd <- .as_numer_denom(rv)
  na <- sincos_args(nd[[1L]]); da <- sincos_args(nd[[2L]])
  handle <- function(rv, num_args, den_args) {
    for (narg in num_args$sin) {
      half <- .div(narg, S_Two)
      if (!is.null(den_args$cos[[half$key]])) fn <- "cos" else if (!is.null(den_args$sin[[half$key]])) fn <- "sin" else next
      rv <- TR11(rv, half)
      den_args[[fn]][[half$key]] <- NULL
    }
    rv
  }
  rv <- handle(rv, na, da)
  handle(rv, da, na)
})

.fu_count_trig <- function(e) {
  n <- 0L
  walk <- function(z) { if (z$type == "Function" && z$name %in% .TRIG_FUNS) n <<- n + 1L; for (a in z$args) walk(a) }
  walk(e); n
}

# _futrig: the greedy rule tree of SymPy
.futrig1 <- function(e) {
  if (!.has_function(e, .TRIG_FUNS)) return(e)
  combine <- NULL; coeff <- NULL
  if (e$type == "Mul") {
    sp <- .as_independent_funs(e, .TRIG_FUNS); coeff <- sp[[1L]]; e <- sp[[2L]]; combine <- function(a, b) .mul(list(a, b))
  } else if (e$type == "Add" && .is_comm(e)) {
    sp <- .split_add_funs(e, .TRIG_FUNS); coeff <- sp[[1L]]; e <- sp[[2L]]; combine <- function(a, b) .add(list(a, b))
  }
  Lops <- function(x) c(.fu_count_trig(x), count_ops(x), .nodes(x), length(x$args), as.numeric(x$type == "Add"))
  trigs <- function(x) .has_function(x, .TRIG_FUNS)
  id <- strat_identity
  tree <- list(id, strat_seq(
    TR3, TR1, TR12,
    function(x) .eapply(factor_sym, x, trigs),
    TR2,
    list(id, function(x) .eapply(.mexpand, x, trigs)),
    TR2i,
    function(x) .eapply(function(i) factor_sym(.expr_normal(i)), x, trigs),
    TR14, TR5, TR10, TR11, .TR11_, TR6,
    function(x) .eapply(factor_sym, x, trigs),
    TR14,
    list(id, function(x) .eapply(.mexpand, x, trigs)),
    TR10i, TRmorrie,
    list(id, TR8),
    list(id, function(x) TR2i(TR2(x))),
    list(function(x) .eapply(expand_mul, TR5(x), trigs), function(x) .eapply(expand_mul, TR15(x), trigs)),
    list(function(x) .eapply(expand_mul, TR6(x), trigs), function(x) .eapply(expand_mul, TR16(x), trigs)),
    TR111,
    list(id, TR2i),
    list(id, function(x) .eapply(expand_mul, TR22(x), trigs)),
    TR1, TR2, TR2i,
    list(id, function(x) .eapply(factor_terms, TR12(x), trigs))))
  e <- greedy(tree, objective = Lops)(e)
  if (!is.null(coeff)) e <- combine(coeff, e)
  e
}

# _osborne and _osbornei: hyperbolic functions as trigonometric ones of
# d*x and back
.osborne <- function(e, d) .bottom_up_args(e, function(rv) {
  if (!(rv$type == "Function" && rv$name %in% .HYP_FUNS)) return(rv)
  a <- rv$args[[1L]]
  a <- if (a$type != "Add") .mul(list(a, d)) else .add(lapply(a$args, function(i) .mul(list(i, d))))
  switch(rv$name,
    sinh = .mul(list(S_I, fcall("sin", list(a)))), cosh = fcall("cos", list(a)),
    tanh = .mul(list(S_I, fcall("tan", list(a)))), coth = .div(fcall("cot", list(a)), S_I),
    sech = fcall("sec", list(a)), csch = .div(fcall("csc", list(a)), S_I))
})
.osbornei <- function(e, d) .bottom_up_args(e, function(rv) {
  if (!(rv$type == "Function" && rv$name %in% c(.TRIG_FUNS, .HYP_FUNS))) return(rv)
  if (rv$name %in% .HYP_FUNS) return(fcall(rv$name, list(xreplace(rv$args[[1L]], list(list(d, S_I))))))
  arg <- rv$args[[1L]]
  parts <- if (arg$type == "Add") arg$args else list(arg)
  dep <- vapply(parts, function(p) .has1(p, d), TRUE)
  const <- .add(parts[!dep]); xx <- .add(parts[dep])
  a <- .add(list(xreplace(xx, list(list(d, S_One))), .mul(list(const, S_I))))
  switch(rv$name,
    sin = .div(fcall("sinh", list(a)), S_I), cos = fcall("cosh", list(a)),
    tan = .div(fcall("tanh", list(a)), S_I), cot = .mul(list(fcall("coth", list(a)), S_I)),
    sec = fcall("sech", list(a)), csc = .mul(list(fcall("csch", list(a)), S_I)))
})

#' @rdname fu
#' @export
hyper_as_trig <- function(rv) {
  rv <- as_sym(rv)
  trigs <- Filter(function(f) f$name %in% .TRIG_FUNS, unique(atoms(rv, "Function")))
  reps <- lapply(trigs, function(t) list(t, Dummy("t")))
  masked <- if (length(reps)) xreplace(rv, reps) else rv
  back <- lapply(reps, function(p) list(p[[2L]], p[[1L]]))
  d <- Dummy("d")
  list(.osborne(masked, d), function(x) {
    y <- .osbornei(x, d)
    if (length(back)) y <- xreplace(y, back)
    collect(signsimp(y), S_I)
  })
}

#' @rdname trigsimp
#' @export
futrig <- function(e, hyper = TRUE) {
  e <- as_sym(e)
  if (!length(e$args)) return(e)
  old <- e
  e <- .bottom_up_args(e, .futrig1)
  if (hyper && .has_function(e, .HYP_FUNS)) {
    coeff <- NULL
    if (e$type == "Add" && .is_comm(e)) { sp <- .split_add_funs(e, .HYP_FUNS); coeff <- sp[[1L]]; e <- sp[[2L]] }
    ht <- hyper_as_trig(e)
    e <- ht[[2L]](.bottom_up_args(ht[[1L]], .futrig1))
    if (!is.null(coeff)) e <- .add(list(e, coeff))
  }
  if (!identical(e$key, old$key) && e$type == "Mul" && e$args[[1L]]$type %in% c("Integer", "Rational")) {
    sp <- as_coeff_mul_split(e)
    e <- .mul(list(sp[[1L]], sp[[2L]]))
  }
  e
}

# _trigsimp_inverse: inverse trigonometric functions of their inverses
.trigsimp_inverse <- function(rv) {
  inv <- c(asin = "sin", acos = "cos", atan = "tan", acot = "cot", asec = "sec", acsc = "csc")
  f <- function(rv) {
    if (rv$type == "Function" && rv$name %in% names(inv) && rv$args[[1L]]$type == "Function" && rv$args[[1L]]$name == inv[[rv$name]])
      return(rv$args[[1L]]$args[[1L]])
    if (rv$type == "Function" && rv$name == "atan2") {
      y <- rv$args[[1L]]; x <- rv$args[[2L]]
      neg <- function(z) { sp <- as_coeff_mul_split(z); is_number_type(sp[[1L]]) && num_sign(sp[[1L]]) < 0 }
      if (neg(y)) return(.neg(f(fcall("atan2", list(.neg(y), x)))))
      if (neg(x)) return(.sub(S_Pi, f(fcall("atan2", list(y, .neg(x))))))
      if (length(x$args) && length(y$args) && identical(x$args[[1L]]$key, y$args[[1L]]$key)) {
        if (y$type == "Function" && y$name == "sin" && x$type == "Function" && x$name == "cos") return(x$args[[1L]])
        if (y$type == "Function" && y$name == "cos" && x$type == "Function" && x$name == "sin") return(.sub(.div(S_Pi, S_Two), x$args[[1L]]))
      }
    }
    rv
  }
  .bottom_up_args(rv, f)
}

#' Trigonometric simplification
#'
#' `trigsimp()` simplifies trigonometric and hyperbolic expressions as
#' SymPy's `trigsimp()`. The default `method = "matching"` runs `futrig()`,
#' a greedy search over the Fu rules. `"fu"` runs `fu()`, `"groebner"`
#' reduces against the ideal of trigonometric identities with
#' `trigsimp_groebner()`, `"combined"` runs the Groebner method with tangent
#' hints and then `futrig()`, and `"old"` (or `old = TRUE`) runs the pattern
#' based algorithm of earlier SymPy versions. `inverse = TRUE` also
#' simplifies inverse trigonometric functions of their inverses.
#' `exptrigsimp()` rewrites exponentials into trigonometric or hyperbolic
#' functions where shorter. `futrig()` is the "matching" method itself.
#'
#' @param e,expr An expression, or a matrix of expressions.
#' @param inverse Also simplify `asin(sin(x))` and the like.
#' @param method One of `"matching"`, `"fu"`, `"groebner"`, `"combined"`
#'   or `"old"`.
#' @param old Use the old method.
#' @param ... Options of the method: `measure` for `"fu"`; `hints`,
#'   `quick`, `order` and `polynomial` for `"groebner"`; `deep`,
#'   `recursive` and `first` for `"old"`.
#' @param hints,quick,order,polynomial See `trigsimp_groebner()` in SymPy:
#'   generator hints (integers and function names), a quick search, the
#'   monomial order and whether to keep the result polynomial.
#' @param hyper Also simplify hyperbolic functions through their
#'   trigonometric counterparts.
#' @return An expression.
#' @examples
#' x <- Symbol("x")
#' trigsimp(sin(x)^4 - cos(x)^4)
#' trigsimp(sin(x) / cos(x))
#' trigsimp(2 * sin(x)^2 + 2 * cos(x)^2, method = "groebner")
#' @export
trigsimp <- function(e, inverse = FALSE, method = "matching", old = FALSE, ...) {
  opts <- list(...)
  if (inherits(e, "symmatrix")) return(.mat_map(e, function(z) trigsimp(z, inverse = inverse, method = method, old = old, ...)))
  e <- as_sym(e)
  if (.is_bd(e)) return(.bd_trigsimp(e, inverse = inverse, method = method, old = old, ...))
  if (old) { opts$method <- method; method <- "old" } else opts$deep <- opts$recursive <- NULL
  gsimp <- function(ex, o) {
    traverse <- function(z) {
      if (!length(z$args)) return(z)
      args <- lapply(z$args, traverse)
      if (z$type %in% c("Function", "Pow")) args <- lapply(args, function(a) do.call(trigsimp_groebner, c(list(a), o)))
      sym_rebuild(z, args)
    }
    new <- traverse(ex)
    if (is_relational(new) || new$type %in% c("And", "Or", "Not", "Xor", "Implies", "Equivalent", "BooleanAtom")) return(new)
    do.call(trigsimp_groebner, c(list(new), o))
  }
  r <- switch(match.arg(method, c("matching", "fu", "groebner", "combined", "old")),
    matching = futrig(e),
    fu = do.call(fu, c(list(e), opts)),
    groebner = gsimp(e, opts),
    combined = futrig(gsimp(e, list(polynomial = TRUE, hints = list(2L, "tan")))),
    old = do.call(trigsimp_old, c(list(e), opts)))
  if (inverse) r <- .trigsimp_inverse(r)
  r
}

.TRIGS_ALL <- c(.TRIG_FUNS, .HYP_FUNS)

#' @rdname trigsimp
#' @param first,deep,recursive Options of the old method: split the
#'   expression into separable factors first, simplify inside every
#'   argument, and simplify common subexpressions one at a time.
#' @export
trigsimp_old <- function(expr, first = TRUE, method = "matching", deep = FALSE, recursive = FALSE, ...) {
  opts <- list(...)
  expr <- as_sym(expr)
  old <- expr
  sub <- function(v) do.call(trigsimp, c(list(v, old = TRUE, method = method, deep = deep, recursive = recursive, first = FALSE), opts))
  if (first) {
    if (!.has_function(expr, .TRIGS_ALL)) return(expr)
    trigsyms <- unique(unlist(lapply(Filter(function(f) f$name %in% .TRIGS_ALL, atoms(expr, "Function")),
                                     function(f) lapply(free_symbols(f), function(s) s)), recursive = FALSE))
    trigsyms <- trigsyms[!duplicated(vapply(trigsyms, .k, ""))]
    if (length(trigsyms) > 1L) {
      d <- separatevars(expr)
      if (d$type == "Mul") { dd <- separatevars(d, dict = TRUE); if (is.list(dd) && !is_sym(dd)) d <- dd }
      if (is.list(d) && !is_sym(d)) {
        expr <- S_One
        for (v in d) {
          was <- v
          v <- expand_mul(v)
          vnew <- sub(v)
          if (identical(vnew$key, v$key)) vnew <- was
          expr <- .mul(list(expr, vnew))
        }
        old <- expr
      } else if (d$type == "Add") {
        for (s in trigsyms) {
          ri <- .as_independent_add(expr, s)
          if (!is_zero_num(ri[[1L]])) {
            expr <- .add(list(ri[[1L]], sub(ri[[2L]])))
            if (expr$type != "Add") break
          }
        }
        old <- expr
      }
    }
  }
  gsimp <- function(ex, d, o) {
    traverse <- function(z) {
      if (!length(z$args)) return(z)
      args <- lapply(z$args, traverse)
      if (z$type %in% c("Function", "Pow")) args <- lapply(args, function(a) do.call(trigsimp_groebner, c(list(a), o)))
      sym_rebuild(z, args)
    }
    if (d) ex <- traverse(ex)
    do.call(trigsimp_groebner, c(list(ex), o))
  }
  tfun <- switch(method,
    matching = function(x, d) .trigsimp_old1(x, d),
    groebner = function(x, d) gsimp(x, d, opts[setdiff(names(opts), "compare")]),
    combined = function(x, d) .trigsimp_old1(gsimp(x, d, list(polynomial = TRUE, hints = list(2L, "tan"))), d),
    stop("unknown method: ", method))
  if (recursive) {
    cs <- cse(expr)
    g <- tfun(cs$reduced[[1L]], deep)
    for (r in rev(cs$replacements)) g <- tfun(subs(g, r[[1L]], r[[2L]]), deep)
    result <- g
  } else result <- tfun(expr, deep)
  if (isTRUE(opts$compare)) {
    f <- futrig(old)
    if (!identical(f$key, result$key)) message("\tfutrig: ", sstr(f))
  }
  result
}

.as_independent_add <- function(e, s) {
  ts <- if (e$type == "Add") e$args else list(e)
  dep <- vapply(ts, function(t) .has1(t, s), TRUE)
  list(.add(ts[!dep]), .add(ts[dep]))
}

.trigsimp_old1 <- function(expr, deep = FALSE) if (.has_function(expr, .TRIGS_ALL)) .trigsimp_old2(expr, deep) else expr

# _dotrig: same head and both with trigonometric or both with hyperbolic
# functions
.dotrig <- function(a, b) identical(a$type, b$type) &&
  ((.has_function(a, .TRIG_FUNS) && .has_function(b, .TRIG_FUNS)) || (.has_function(a, .HYP_FUNS) && .has_function(b, .HYP_FUNS)))

# _replace_mul_fpowxgpow: f(x)^e*g(x)^rexp(e) to h(x)^rexph(e)
.replace_mul_fpowxgpow <- function(expr, f, g, rexp, h, rexph) {
  fargs <- list(); gargs <- list(); keys <- list(); args <- list()
  for (x in expr$args) {
    if (x$type == "Pow" || (x$type == "Function" && x$name %in% c(f, g))) {
      be <- .base_exp(x); b <- be[[1L]]; ee <- be[[2L]]
      if (isTRUE(ask_fact(b, "positive")) || isTRUE(ask_fact(ee, "integer"))) {
        if (b$type == "Function" && b$name == f) { k <- b$args[[1L]]$key; keys[[k]] <- b$args[[1L]]; fargs[[k]] <- .add(list(if (is.null(fargs[[k]])) S_Zero else fargs[[k]], ee)); next }
        if (b$type == "Function" && b$name == g) { k <- b$args[[1L]]$key; keys[[k]] <- b$args[[1L]]; gargs[[k]] <- .add(list(if (is.null(gargs[[k]])) S_Zero else gargs[[k]], ee)); next }
      }
    }
    args[[length(args) + 1L]] <- x
  }
  hit <- FALSE
  for (k in intersect(names(fargs), names(gargs))) {
    fe <- fargs[[k]]; ge <- gargs[[k]]
    if (identical(fe$key, rexp(ge)$key)) {
      hb <- if (is.null(h)) S_One else fcall(h, list(keys[[k]]))
      args[[length(args) + 1L]] <- .pow(hb, rexph(fe))
      fargs[[k]] <- NULL; gargs[[k]] <- NULL
      hit <- TRUE
    }
  }
  if (!hit) return(expr)
  for (k in names(fargs)) args[[length(args) + 1L]] <- .pow(fcall(f, list(keys[[k]])), fargs[[k]])
  for (k in names(gargs)) args[[length(args) + 1L]] <- .pow(fcall(g, list(keys[[k]])), gargs[[k]])
  .mul(args)
}

.match_div_rewrite <- function(expr, i) {
  idn <- function(x) x; midn <- function(x) .neg(x)
  spec <- switch(as.character(i),
    "0" = list("sin", "cos", midn, "tan", idn), "1" = list("tan", "cos", idn, "sin", idn),
    "2" = list("cot", "sin", idn, "cos", idn), "3" = list("tan", "sin", midn, "cos", midn),
    "4" = list("cot", "cos", midn, "sin", midn), "5" = list("cot", "tan", idn, NULL, idn),
    "8" = list("sinh", "cosh", midn, "tanh", idn), "9" = list("tanh", "cosh", idn, "sinh", idn),
    "10" = list("coth", "sinh", idn, "cosh", idn), "11" = list("tanh", "sinh", midn, "cosh", midn),
    "12" = list("coth", "cosh", midn, "sinh", midn), "13" = list("coth", "tanh", idn, NULL, idn),
    NULL)
  if (is.null(spec)) return(NULL)
  .replace_mul_fpowxgpow(expr, spec[[1L]], spec[[2L]], spec[[3L]], spec[[4L]], spec[[5L]])
}

.trigpat_cache <- new.env()
.trigpats <- function() {
  if (!is.null(.trigpat_cache$p)) return(.trigpat_cache$p)
  a <- Wild("a"); b <- Wild("b"); c <- Wild("c"); d <- Wild("d")
  f <- function(nm, x) fcall(nm, list(x))
  P <- function(x, n) .pow(x, n)
  M <- function(...) .mul(list(...)); A <- function(...) .add(list(...))
  nc <- .neg(c)
  div <- list(
    list(M(a, P(f("sin", b), c), P(f("cos", b), nc)), M(a, P(f("tan", b), c)), f("sin", b), f("cos", b)),
    list(M(a, P(f("tan", b), c), P(f("cos", b), c)), M(a, P(f("sin", b), c)), f("sin", b), f("cos", b)),
    list(M(a, P(f("cot", b), c), P(f("sin", b), c)), M(a, P(f("cos", b), c)), f("sin", b), f("cos", b)),
    list(M(a, P(f("tan", b), c), P(f("sin", b), nc)), M(a, P(f("cos", b), nc)), f("sin", b), f("cos", b)),
    list(M(a, P(f("cot", b), c), P(f("cos", b), nc)), M(a, P(f("sin", b), nc)), f("sin", b), f("cos", b)),
    list(M(a, P(f("cot", b), c), P(f("tan", b), c)), a, f("sin", b), f("cos", b)),
    list(M(a, P(A(f("cos", b), S_One), c), P(A(f("cos", b), S_NegativeOne), c)), M(a, P(.neg(P(f("sin", b), S_Two)), c)),
         A(f("cos", b), S_One), A(f("cos", b), S_NegativeOne)),
    list(M(a, P(A(f("sin", b), S_One), c), P(A(f("sin", b), S_NegativeOne), c)), M(a, P(.neg(P(f("cos", b), S_Two)), c)),
         A(f("sin", b), S_One), A(f("sin", b), S_NegativeOne)),
    list(M(a, P(f("sinh", b), c), P(f("cosh", b), nc)), M(a, P(f("tanh", b), c)), S_One, S_One),
    list(M(a, P(f("tanh", b), c), P(f("cosh", b), c)), M(a, P(f("sinh", b), c)), S_One, S_One),
    list(M(a, P(f("coth", b), c), P(f("sinh", b), c)), M(a, P(f("cosh", b), c)), S_One, S_One),
    list(M(a, P(f("tanh", b), c), P(f("sinh", b), nc)), M(a, P(f("cosh", b), nc)), S_One, S_One),
    list(M(a, P(f("coth", b), c), P(f("cosh", b), nc)), M(a, P(f("sinh", b), nc)), S_One, S_One),
    list(M(a, P(f("coth", b), c), P(f("tanh", b), c)), a, S_One, S_One),
    list(M(c, A(f("tanh", a), f("tanh", b)), P(A(S_One, M(f("tanh", a), f("tanh", b))), S_NegativeOne)),
         M(f("tanh", A(a, b)), c), S_One, S_One))
  add <- list(
    list(A(M(c, f("sin", a), f("cos", b)), M(c, f("cos", a), f("sin", b)), d), A(M(f("sin", A(a, b)), c), d)),
    list(A(M(c, f("cos", a), f("cos", b)), .neg(M(c, f("sin", a), f("sin", b))), d), A(M(f("cos", A(a, b)), c), d)),
    list(A(M(c, f("sin", a), f("cos", b)), .neg(M(c, f("cos", a), f("sin", b))), d), A(M(f("sin", .sub(a, b)), c), d)),
    list(A(M(c, f("cos", a), f("cos", b)), M(c, f("sin", a), f("sin", b)), d), A(M(f("cos", .sub(a, b)), c), d)),
    list(A(M(c, f("sinh", a), f("cosh", b)), M(c, f("sinh", b), f("cosh", a)), d), A(M(f("sinh", A(a, b)), c), d)),
    list(A(M(c, f("cosh", a), f("cosh", b)), M(c, f("sinh", a), f("sinh", b)), d), A(M(f("cosh", A(a, b)), c), d)))
  sq <- function(x) P(x, S_Two)
  ident <- list(
    list(M(a, sq(f("sin", b))), .sub(a, M(a, sq(f("cos", b))))),
    list(M(a, sq(f("tan", b))), .sub(M(a, sq(.inv(f("cos", b)))), a)),
    list(M(a, sq(f("cot", b))), .sub(M(a, sq(.inv(f("sin", b)))), a)),
    list(M(a, f("sin", A(b, c))), M(a, A(M(f("sin", b), f("cos", c)), M(f("sin", c), f("cos", b))))),
    list(M(a, f("cos", A(b, c))), M(a, .sub(M(f("cos", b), f("cos", c)), M(f("sin", b), f("sin", c))))),
    list(M(a, f("tan", A(b, c))), M(a, .div(A(f("tan", b), f("tan", c)), .sub(S_One, M(f("tan", b), f("tan", c)))))),
    list(M(a, sq(f("sinh", b))), .sub(M(a, sq(f("cosh", b))), a)),
    list(M(a, sq(f("tanh", b))), .sub(a, M(a, sq(.inv(f("cosh", b)))))),
    list(M(a, sq(f("coth", b))), A(a, M(a, sq(.inv(f("sinh", b)))))),
    list(M(a, f("sinh", A(b, c))), M(a, A(M(f("sinh", b), f("cosh", c)), M(f("sinh", c), f("cosh", b))))),
    list(M(a, f("cosh", A(b, c))), M(a, A(M(f("cosh", b), f("cosh", c)), M(f("sinh", b), f("sinh", c))))),
    list(M(a, f("tanh", A(b, c))), M(a, .div(A(f("tanh", b), f("tanh", c)), A(S_One, M(f("tanh", b), f("tanh", c)))))))
  art <- function(base, res, ex) list(base, res, ex)
  ad <- M(a, d)
  arts <- list(
    art(A(a, .neg(M(a, sq(f("cos", b)))), c), A(M(a, sq(f("sin", b))), c), "cos"),
    art(A(a, .neg(M(a, sq(.inv(f("cos", b))))), c), A(.neg(M(a, sq(f("tan", b)))), c), "cos"),
    art(A(a, .neg(M(a, sq(.inv(f("sin", b))))), c), A(.neg(M(a, sq(f("cot", b)))), c), "sin"),
    art(A(a, .neg(M(a, sq(f("cosh", b)))), c), A(.neg(M(a, sq(f("sinh", b)))), c), "cosh"),
    art(A(a, .neg(M(a, sq(.inv(f("cosh", b))))), c), A(M(a, sq(f("tanh", b))), c), "cosh"),
    art(A(a, M(a, sq(.inv(f("sinh", b)))), c), A(M(a, sq(f("coth", b))), c), "sinh"),
    art(A(ad, .neg(M(ad, sq(f("cos", b)))), c), A(M(ad, sq(f("sin", b))), c), "cos"),
    art(A(ad, .neg(M(ad, sq(.inv(f("cos", b))))), c), A(.neg(M(ad, sq(f("tan", b)))), c), "cos"),
    art(A(ad, .neg(M(ad, sq(.inv(f("sin", b))))), c), A(.neg(M(ad, sq(f("cot", b)))), c), "sin"),
    art(A(ad, .neg(M(ad, sq(f("cosh", b)))), c), A(.neg(M(ad, sq(f("sinh", b)))), c), "cosh"),
    art(A(ad, .neg(M(ad, sq(.inv(f("cosh", b))))), c), A(M(ad, sq(f("tanh", b))), c), "cosh"),
    art(A(ad, M(ad, sq(.inv(f("sinh", b)))), c), A(M(ad, sq(f("coth", b))), c), "sinh"))
  .trigpat_cache$p <- list(a = a, b = b, c = c, d = d, div = div, add = add, ident = ident, arts = arts)
  .trigpat_cache$p
}

.wsubs <- function(e, m) if (length(m)) xreplace(e, unname(lapply(names(m), function(k) list(m[[k]]$w, m[[k]]$v)))) else e

# sym_match result as a list keyed by wild key holding the wild and value
.wmatch <- function(e, pattern, wilds) {
  r <- tryCatch(sym_match(e, pattern), error = function(err) NULL)
  if (is.null(r)) return(NULL)
  out <- list()
  for (w in wilds) { v <- r[[w$name]]; if (!is.null(v)) out[[w$key]] <- list(w = w, v = v) }
  out
}
.wget <- function(m, w) if (is.null(m[[w$key]])) NULL else m[[w$key]]$v

.trig_arg_in <- function(e, val) any(vapply(Filter(function(f) f$name %in% .TRIGS_ALL, atoms(e, "Function")),
                                            function(f) identical(f$args[[1L]]$key, val$key), TRUE))

# __trigsimp
.trigsimp_old2 <- function(expr, deep = FALSE) {
  tp <- .trigpats(); a <- tp$a; b <- tp$b; c <- tp$c; d <- tp$d; W <- list(a, b, c, d)
  if (expr$type == "Mul") {
    for (i in seq_along(tp$div) - 1L) {
      pt <- tp$div[[i + 1L]]
      if (!.dotrig(expr, pt[[1L]])) next
      newexpr <- .match_div_rewrite(expr, i)
      if (!is.null(newexpr)) {
        if (!identical(newexpr$key, expr$key)) { expr <- newexpr; break } else next
      }
      res <- .wmatch(expr, pt[[1L]], W)
      cv <- .wget(res, c)
      if (!is.null(res) && !is.null(cv) && !is_zero_num(cv)) {
        if (!isTRUE(ask_fact(cv, "integer"))) {
          if (!isTRUE(ask_fact(.wsubs(pt[[3L]], res), "positive"))) next
          if (!isTRUE(ask_fact(.wsubs(pt[[4L]], res), "positive"))) next
        }
        av <- .wget(res, a); bv <- .wget(res, b)
        if (!is.null(av) && !is.null(bv) && .trig_arg_in(av, bv)) next
        expr <- .wsubs(pt[[2L]], res)
        break
      }
    }
  }
  if (expr$type == "Add") {
    args <- lapply(expr$args, function(term) {
      term <- .trigsimp_old1(term, deep)
      for (pr in tp$ident) {
        res <- .wmatch(term, pr[[1L]], W)
        if (!is.null(res)) { term <- .wsubs(pr[[2L]], res); break }
      }
      term
    })
    if (!identical(vapply(args, .k, ""), vapply(expr$args, .k, ""))) {
      expr <- .add(args)
      ex <- expand(expr)
      if (count_ops(ex) < count_ops(expr)) expr <- ex
    }
    if (expr$type == "Add") {
      for (pr in tp$add) {
        if (!.dotrig(expr, pr[[1L]])) next
        expr <- TR10i(expr)
        if (.has_function(expr, .HYP_FUNS)) {
          res <- .wmatch(expr, pr[[1L]], W)
          dv <- .wget(res, d)
          if (is.null(res) || is.null(.wget(res, a)) || is.null(.wget(res, b)) ||
              (!is.null(dv) && (.trig_arg_in(dv, .wget(res, a)) || .trig_arg_in(dv, .wget(res, b))))) next
          expr <- .wsubs(pr[[2L]], res)
          break
        }
      }
    }
    for (pr in tp$arts) {
      if (!.dotrig(expr, pr[[1L]])) next
      exf <- pr[[3L]]
      a_t <- Wild("a_t", properties = list(function(v) !.has_function(v, exf)))
      pattern <- xreplace(pr[[1L]], list(list(a, a_t))); result <- xreplace(pr[[2L]], list(list(a, a_t)))
      W2 <- list(a_t, b, c, d)
      m <- .wmatch(expr, pattern, W2)
      was <- NULL
      while (!is.null(m) && !identical(was, expr$key)) {
        was <- expr$key
        at <- .wget(m, a_t); cv <- .wget(m, c); if (is.null(cv)) cv <- S_Zero
        if (is.null(at) || is_zero_num(at) || any(vapply(if (cv$type == "Add") cv$args else list(), function(t) identical(t$key, .neg(at)$key), TRUE)) ||
            is_zero_num(.add(list(at, cv)))) break
        dv <- .wget(m, d)
        if (!is.null(dv) && is_zero_num(.add(list(.mul(list(at, dv)), cv)))) break
        if (is.null(m[[c$key]])) m[[c$key]] <- list(w = c, v = S_Zero)
        expr <- .wsubs(result, m)
        m <- .wmatch(expr, pattern, W2)
      }
    }
  } else if (expr$type %in% c("Mul", "Pow") || (deep && length(expr$args))) {
    expr <- sym_rebuild(expr, lapply(expr$args, function(z) .trigsimp_old1(z, deep)))
  }
  # rewrite in exponentials when that adds no new exponential
  if (.has_function(expr, .TRIGS_ALL)) {
    e0 <- vapply(Filter(function(f) f$name == "exp", atoms(expr, "Function")), .k, "")
    new <- if (deep) rewrite(expr, "exp") else .rewrite_top(expr, "exp")
    fnew <- factor_sym(new)
    if (!identical(fnew$key, new$key) && count_ops(fnew) < count_ops(new)) new <- fnew
    e1 <- vapply(Filter(function(f) f$name == "exp", atoms(new, "Function")), .k, "")
    if (!length(setdiff(e1, e0))) expr <- new
  }
  expr
}

# rewrite of the top node only (deep = False)
.rewrite_top <- function(e, target) {
  rules <- .rewrite_rules[[target]]
  if (e$type == "Function" && !is.null(rules[[e$name]])) {
    r <- do.call(rules[[e$name]], e$args)
    if (!is.null(r)) return(r)
  }
  e
}

# rewrite of the functions named in `names` only (SymPy rewrite(pattern, rule))
.rewrite_funs <- function(e, target, names) {
  rules <- .rewrite_rules[[target]]
  .map_bottom_up(e, function(z) {
    if (z$type == "Function" && z$name %in% names && !is.null(rules[[z$name]])) {
      r <- do.call(rules[[z$name]], z$args)
      if (!is.null(r)) return(r)
    }
    z
  })
}

.exptrig_choice <- function(e) {
  choices <- list(e)
  if (.has_function(e, .TRIGS_ALL)) {
    op <- .rewrite_funs(e, "exp", .TRIGS_ALL)
    choices[[2L]] <- if (e$type == "Add") factor_terms(op) else op
    choices[[3L]] <- rewrite(e, "sin")
  }
  choices[[length(choices) + 1L]] <- rewrite(e, "cos")
  ops <- vapply(choices, count_ops, 0)
  choices[[which.min(ops)]]
}

# the Mul pass of exptrigsimp: 1 +- exp(x) factors against exp(x/2) and
# each other
.exptrig_mul <- function(rv) {
  if (rv$type != "Mul") return(rv)
  bases <- list(); ex <- list()
  put <- function(b, e) {
    k <- b$key
    if (is.null(bases[[k]])) { bases[[k]] <<- b; ex[[k]] <<- e } else ex[[k]] <<- .add(list(ex[[k]], e))
  }
  for (t in rv$args) {
    if (t$type == "Function" && t$name == "exp") put(S_E, t$args[[1L]])
    else { be <- .base_exp(t); put(be[[1L]], be[[2L]]) }
  }
  nb <- bases; nd <- ex
  get <- function(b) if (is.null(nd[[b$key]])) S_Zero else nd[[b$key]]
  inc <- function(b, v) { if (is.null(nb[[b$key]])) nb[[b$key]] <<- b; nd[[b$key]] <<- .add(list(get(b), v)) }
  signlog <- function(z, sign = 1L) {
    if (is_const(z, "E")) return(list(sign, S_One))
    if (z$type == "Function" && z$name == "exp") return(list(sign, z$args[[1L]]))
    if (sign == 1L) return(signlog(.neg(z), -1L))
    NULL
  }
  ee <- if (is.null(ex[[S_E$key]])) S_Zero else ex[[S_E$key]]
  Ex <- function(x) fcall("exp", list(x))
  for (k in names(bases)) {
    kb <- bases[[k]]
    if (kb$type != "Add" || length(kb$args) != 2L) next
    c <- kb$args[[1L]]
    sl <- signlog(.div(kb$args[[2L]], c))
    if (is.null(sl)) next
    sg <- sl[[1L]]; x <- sl[[2L]]
    m <- ex[[k]]
    inc(kb, .neg(m))
    if (identical(ee$key, .div(.neg(.mul(list(x, m))), S_Two)$key)) {
      inc(S_E, .neg(ee)); ee <- S_Zero
      if (sg == 1L) inc(.mul(list(S_Two, c, fcall("cosh", list(.div(x, S_Two))))), m)
      else inc(.mul(list(S_NegativeOne, S_Two, c, fcall("sinh", list(.div(x, S_Two))))), m)
    } else {
      other <- .sub(S_One, .mul(list(.int(sg), Ex(x))))
      if (identical(get(other)$key, .neg(m)$key)) {
        nb[[other$key]] <- NULL; nd[[other$key]] <- NULL
        if (sg == 1L) inc(.neg(.div(c, fcall("tanh", list(.div(x, S_Two))))), m)
        else inc(.neg(.mul(list(c, fcall("tanh", list(.div(x, S_Two)))))), m)
      } else {
        inc(.add(list(S_One, .mul(list(.int(sg), Ex(x))))), m)
        inc(c, m)
      }
    }
  }
  .mul(lapply(names(nb), function(k) {
    b <- nb[[k]]; v <- nd[[k]]
    if (is_const(b, "E")) Ex(v) else .pow(b, v)
  }))
}

#' @rdname trigsimp
#' @export
exptrigsimp <- function(expr) {
  expr <- as_sym(expr)
  newexpr <- .bottom_up_args(expr, .exptrig_choice)
  newexpr <- .bottom_up_args(newexpr, .exptrig_mul)
  if (.has_function(newexpr, .HYP_FUNS)) {
    ht <- hyper_as_trig(newexpr)
    newexpr <- ht[[2L]](TR2i(ht[[1L]]))
  }
  if (.has_function(newexpr, .TRIG_FUNS)) newexpr <- TR2i(newexpr)
  if (!(.contains_key(newexpr, S_I$key) && !.contains_key(expr, S_I$key))) expr <- newexpr
  expr
}
