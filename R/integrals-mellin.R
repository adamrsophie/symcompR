# Mellin transforms, after sympy/integrals/transforms.py: IntegralTransform.doit
# for the Mellin pair, _mellin_transform, _rewrite_sin, _rewrite_gamma and
# _inverse_mellin_transform.

.it_error <- function(transform, msg) {
  stop(structure(class = c("IntegralTransformError", "NotImplementedError", "error", "condition"),
                 list(message = paste0(transform, " Transform could not be computed: ", msg, "."), call = NULL)))
}
.strip_error <- function(msg) {
  stop(structure(class = c("MellinTransformStripError", "error", "condition"), list(message = msg, call = NULL)))
}
.notimpl_error <- function(msg) {
  stop(structure(class = c("NotImplementedError", "error", "condition"), list(message = msg, call = NULL)))
}

.disjuncts <- function(e) if (e$type == "Or") e$args else list(e)
.conjuncts <- function(e) if (e$type == "And") e$args else list(e)

# the lesser and greater sides of an inequality

.tr_simplify <- function(expr, doit) if (doit) simplify(powdenest(piecewise_fold(expr), polar = TRUE)) else expr

.default_integrator <- function(f, x) integrate(f, list(x, S_Zero, S_Infinity))

# _mellin_transform: list(F, list(a, b), aux)
.mellin_core <- function(f, x, s_, integrator = .default_integrator, simplify = TRUE) {
  s <- .mi_dummy("s", "mellin-transform", f)
  F <- integrator(.mul(list(.pow(x, .sub(s, S_One)), f)), x)
  if (!.has_type(F, "Integral"))
    return(list(.tr_simplify(subs(F, s, s_), simplify), list(S_NegativeInfinity, S_Infinity), S_true))
  if (F$type != "Piecewise") .it_error("Mellin", "could not compute integral")
  pc <- F$args[[1L]]
  F <- pc$args[[1L]]; cond <- pc$args[[2L]]
  if (.has_type(F, "Integral")) .it_error("Mellin", "integral in unexpected form")
  process_conds <- function(cond) {
    a <- S_NegativeInfinity; b <- S_Infinity; aux <- S_true
    t <- Dummy("t", real = TRUE)
    res <- fcall("re", list(s))
    for (c in .conjuncts(to_cnf(cond))) {
      a_ <- S_Infinity; b_ <- S_NegativeInfinity; aux_ <- list()
      for (d in .disjuncts(c)) {
        d_ <- .map_bottom_up(d, function(z) if (z$type == "Function" && z$name == "re") .as_real_imag(z$args[[1L]])[[1L]] else z)
        d_ <- subs(d_, res, t)
        if (!is_relational(d) || d$type %in% c("Equality", "Unequality") || .has1(d_, s) || !.has1(d_, t)) {
          aux_[[length(aux_) + 1L]] <- d; next
        }
        soln <- .solve_inequality(d_, t)
        if (!is_relational(soln) || soln$type %in% c("Equality", "Unequality")) {
          aux_[[length(aux_) + 1L]] <- d; next
        }
        if (identical(.rel_lts(soln)$key, t$key)) b_ <- Max(.rel_gts(soln), b_)
        else a_ <- Min(.rel_lts(soln), a_)
      }
      if (!is_const(a_, "oo") && !identical(a_$key, b$key)) a <- Max(a_, a)
      else if (!is_const(b_, "-oo") && !identical(b_$key, a$key)) b <- Min(b_, b)
      else aux <- And(aux, .mi_OR(aux_))
    }
    list(a, b, aux)
  }
  conds <- lapply(.disjuncts(cond), process_conds)
  conds <- Filter(function(z) !identical(z[[3L]]$key, S_false$key), conds)
  if (!length(conds)) .it_error("Mellin", "no convergence found")
  if (length(conds) > 1L) {
    # list.sort with key (a - b, count_ops(aux)); SymPy comparisons that
    # cannot be decided raise TypeError
    key <- lapply(conds, function(z) list(.sub(z[[1L]], z[[2L]]), count_ops(z[[3L]])))
    less <- function(u, v) {
      if (!identical(u[[1L]]$key, v[[1L]]$key)) {
        if (.is_true(Lt(u[[1L]], v[[1L]]))) return(TRUE)
        if (.is_true(Ge(u[[1L]], v[[1L]]))) return(FALSE)
        stop("cannot determine truth value of Relational")
      }
      u[[2L]] < v[[2L]]
    }
    ord <- seq_along(conds)
    for (i in seq_along(ord)[-1L]) {
      j <- i
      while (j > 1L && less(key[[ord[j]]], key[[ord[j - 1L]]])) { tmp <- ord[j]; ord[j] <- ord[j - 1L]; ord[j - 1L] <- tmp; j <- j - 1L }
    }
    conds <- conds[ord]
  }
  list(.tr_simplify(subs(F, s, s_), simplify), list(conds[[1L]][[1L]], conds[[1L]][[2L]]), conds[[1L]][[3L]])
}

.applied_undef_in <- function(f, x) any(vapply(atoms(f, "Function"), function(g) .is_undefined_application(g) && .has1(g, x), NA))

# IntegralTransform.doit for MellinTransform
.mellin_doit <- function(f, x, s, integrator, simplify, needeval, noconds) {
  T <- NULL
  if (!.applied_undef_in(f, x))
    T <- tryCatch(.mellin_core(f, x, s, integrator, simplify), IntegralTransformError = function(e) NULL)
  fn <- if (f$type == "Add") f else expand_mul(f)
  if (!is.null(T)) return(if (noconds) T[[1L]] else T)
  if (fn$type == "Add") {
    res <- lapply(.sympy_args(fn), function(a) .mellin_doit(a, x, s, integrator, simplify, needeval, noconds))
    ress <- list(); extra <- list()
    for (r in res) {
      if (is_sym(r)) { ress[[length(ress) + 1L]] <- r; next }
      ress[[length(ress) + 1L]] <- r[[1L]]
      extra[[length(extra) + 1L]] <- r[-1L]
    }
    tot <- if (isTRUE(simplify)) simplify(.add(ress)) else .add(ress)
    if (!length(extra)) return(tot)
    col <- tryCatch({
      a <- lapply(extra, function(z) z[[1L]][[1L]]); b <- lapply(extra, function(z) z[[1L]][[2L]])
      cc <- do.call(And, lapply(extra, `[[`, 2L))
      A <- do.call(Max, a); B <- do.call(Min, b)
      if (.is_true(Ge(A, B)) || identical(cc$key, S_false$key)) .it_error("Mellin", "no combined convergence.")
      list(list(A, B), cc)
    }, IntegralTransformError = function(e) NULL)
    if (!is.null(col)) return(c(list(tot), col))
  }
  if (needeval) .it_error("Mellin", "needeval")
  cm <- .as_coeff_mul_dep(fn, x)
  .mul(list(cm[[1L]], MellinTransform(.mul(cm[[2L]]), x, s)))
}

#' Mellin transforms
#'
#' Ports of SymPy's `mellin_transform` and `inverse_mellin_transform`
#' (sympy/integrals/transforms.py). `mellin_transform(f, x, s)` returns
#' `list(F, list(a, b), cond)`: the transform, its fundamental strip
#' `a < re(s) < b` and an auxiliary condition; with `noconds = TRUE` only `F`.
#' `inverse_mellin_transform(F, s, x, strip)` inverts `F` over the strip by
#' rewriting it as a Meijer G function of `x`.
#'
#' @param f,F Expressions.
#' @param x,s Variables.
#' @param strip `list(a, b)`; either may be `NULL`.
#' @param integrator A function of `(f, x)` computing the integral over
#'   `(0, oo)`; by default [integrate()].
#' @param simplify Simplify the result.
#' @param needeval Signal an error instead of returning an unevaluated
#'   transform.
#' @param noconds Return only the transform.
#' @param as_meijerg Leave the inverse as a Meijer G function.
#' @return The transform.
#' @examples
#' x <- Symbol("x"); s <- Symbol("s")
#' mellin_transform(exp(-x), x, s)
#' inverse_mellin_transform(gamma(s), s, x, list(0, oo))
#' @export
mellin_transform <- function(f, x, s, integrator = NULL, simplify = TRUE, needeval = FALSE, noconds = FALSE) {
  f <- as_sym(f); x <- as_sym(x); s <- as_sym(s)
  if (is.null(integrator)) integrator <- .default_integrator
  .mellin_doit(f, x, s, integrator, simplify, needeval, noconds)
}

# _rewrite_sin
.rewrite_sin <- function(m_n, s, a, b) {
  m <- expand_mul(.div(m_n[[1L]], S_Pi)); n <- expand_mul(.div(m_n[[2L]], S_Pi))
  r <- fcall("ceiling", list(.sub(.neg(.mul(list(m, a))), .as_real_imag(n)[[1L]])))
  list(gamma(.add(list(.mul(list(m, s)), n, r))), gamma(.sub(.sub(.sub(S_One, n), r), .mul(list(m, s)))),
       .mul(list(.pow(S_NegativeOne, r), S_Pi)))
}

.rg_rational_pq <- function(v) { q <- gmp::as.bigq(v$value); list(p = gmp::numerator(q), q = gmp::denominator(q)) }

# _rewrite_gamma: list(list(an, ap), list(bm, bq), arg, exponent, fac)
.mellin_rewrite_gamma <- function(f, s, a, b) {
  a_ <- if (is.null(a)) NULL else as_sym(a); b_ <- if (is.null(b)) NULL else as_sym(b)
  left <- function(c, is_numer) {
    c <- expand(fcall("re", list(c)))
    if (is.null(a_) && !is.null(b_) && is_const(b_, "oo")) return(TRUE)
    if (is.null(a_)) return(.mi_bool(Lt(c, b_)))
    if (is.null(b_)) return(.mi_bool(Le(c, a_)))
    if (.is_true(Ge(c, b_))) return(FALSE)
    if (.is_true(Le(c, a_))) return(TRUE)
    if (is_numer) return(NA)
    if (length(free_symbols(a_)) || length(free_symbols(b_)) || length(free_symbols(c))) return(NA)
    .strip_error("Pole inside critical strip?")
  }
  mults <- list()
  for (g in atoms(f, "Function")) {
    if (!(g$name %in% c("gamma", "sin", "cos", "tan", "cot")) || !.has1(g, s)) next
    ar <- g$args[[1L]]
    if (ar$type == "Add") ar <- .add(Filter(function(z) .has1(z, s), ar$args))
    cf <- .as_coeff_mul_dep(ar, s)[[1L]]
    mults[[length(mults) + 1L]] <- if (g$name == "gamma") cf else .div(cf, S_Pi)
  }
  mults <- lapply(mults, function(z) if (isTRUE(ask_fact(z, "extended_real"))) Abs(z) else z)
  common <- S_One
  for (z in mults) if (!(z$type %in% c("Integer", "Rational"))) { common <- z; break }
  mults <- lapply(mults, function(z) .div(z, common))
  if (!(all(vapply(mults, function(z) z$type %in% c("Integer", "Rational"), NA)) && isTRUE(ask_fact(common, "extended_real"))))
    .it_error("Gamma", "Nonrational multiplier")
  lq <- Reduce(gmp::lcm.bigz, lapply(mults, function(z) .rg_rational_pq(z)$q), gmp::as.bigz(1))
  s_mult <- .div(common, .zs(lq))
  if (identical(s_mult$key, common$key)) {
    if (!length(mults)) s_mult <- common
    else s_mult <- .mul(list(common, .zs(Reduce(gmp::gcd.bigz, lapply(mults, function(z) .rg_rational_pq(z)$p)))))
  }
  f <- subs(f, s, .div(s, s_mult))
  fac <- .inv(s_mult); exponent <- .inv(s_mult)
  if (!is.null(a_)) a_ <- .mul(list(a_, s_mult))
  if (!is.null(b_)) b_ <- .mul(list(b_, s_mult))
  nd <- .as_numer_denom(f)
  margs <- function(e) if (e$type == "Mul") .sympy_args(e) else list(e)
  args <- c(lapply(margs(nd[[1L]]), function(z) list(z, TRUE)), lapply(margs(nd[[2L]]), function(z) list(z, FALSE)))
  facs <- list(); dfacs <- list(); numer_gammas <- list(); denom_gammas <- list(); exponentials <- list()
  exception <- function(fact) .it_error("Inverse Mellin", paste0("Unrecognised form '", .str(fact), "'"))
  linear_arg <- function(arg, fact) {
    if (!isTRUE(is_polynomial(arg, s))) exception(fact)
    cl <- .coeff_list(expand(arg), s)
    if (length(cl) != 2L) exception(fact)
    cl
  }
  pushg <- function(is_numer, upper, v) {
    if (xor(is_numer, upper)) denom_gammas[[length(denom_gammas) + 1L]] <<- v
    else numer_gammas[[length(numer_gammas) + 1L]] <<- v
  }
  pushf <- function(is_numer, v) {
    if (is_numer) facs[[length(facs) + 1L]] <<- v else dfacs[[length(dfacs) + 1L]] <<- v
  }
  while (length(args)) {
    it <- args[[length(args)]]; args[[length(args)]] <- NULL
    fact <- it[[1L]]; is_numer <- it[[2L]]
    if (!.has1(fact, s)) { pushf(is_numer, fact); next }
    if (fact$type == "Pow" || (fact$type == "Function" && fact$name == "exp")) {
      if (fact$type == "Pow") { base <- fact$args[[1L]]; ex <- fact$args[[2L]] }
      else { base <- fcall("exp_polar", list(S_One)); ex <- fact$args[[1L]] }
      if (ex$type == "Integer") {
        cond <- is_numer
        if (.int_val(ex) < 0) cond <- !cond
        args <- c(args, rep(list(list(base, cond)), abs(.int_val(ex))))
        next
      } else if (!.has1(base, s)) {
        ab <- linear_arg(ex, fact)
        if (!is_numer) base <- .inv(base)
        exponentials[[length(exponentials) + 1L]] <- .pow(base, ab[[1L]])
        facs[[length(facs) + 1L]] <- .pow(base, ab[[2L]])
      } else exception(fact)
    } else if (isTRUE(is_polynomial(fact, s))) {
      cl <- .coeff_list(expand(fact), s)
      if (length(cl) - 1L != 1L) {
        coeff <- cl[[1L]]
        rts <- .rg_roots(fact, s)
        pushf(is_numer, coeff)
        args <- c(args, lapply(rts, function(cc) list(.sub(s, cc), is_numer)))
        next
      }
      A <- cl[[1L]]; c <- cl[[2L]]
      pushf(is_numer, A)
      c <- .div(c, .neg(A))
      lv <- left(c, is_numer)
      if (isTRUE(lv)) {
        pushg(is_numer, TRUE, list(S_One, .add(list(.neg(c), S_One))))
        pushg(is_numer, FALSE, list(S_One, .neg(c)))
      } else {
        pushf(is_numer, S_NegativeOne)
        pushg(is_numer, TRUE, list(S_NegativeOne, .add(list(c, S_One))))
        pushg(is_numer, FALSE, list(S_NegativeOne, c))
      }
    } else if (fact$type == "Function" && fact$name == "gamma") {
      ab <- linear_arg(fact$args[[1L]], fact); A <- ab[[1L]]; B <- ab[[2L]]
      if (is_numer) {
        lv <- left(.div(.neg(B), A), is_numer)
        if ((.is_true(Gt(A, 0)) && identical(lv, FALSE)) || (.is_true(Lt(A, 0)) && identical(lv, TRUE)))
          .notimpl_error("Gammas partially over the strip.")
      }
      pushg(is_numer, TRUE, list(A, B))
    } else if (fact$type == "Function" && fact$name == "sin") {
      a0 <- fact$args[[1L]]
      if (is_numer) { g1 <- gamma(.div(a0, S_Pi)); g2 <- gamma(.sub(S_One, .div(a0, S_Pi))); fac_ <- S_Pi }
      else { rs <- .rewrite_sin(linear_arg(a0, fact), s, a_, b_); g1 <- rs[[1L]]; g2 <- rs[[2L]]; fac_ <- rs[[3L]] }
      args <- c(args, list(list(g1, !is_numer), list(g2, !is_numer)))
      pushf(is_numer, fac_)
    } else if (fact$type == "Function" && fact$name == "tan") {
      a0 <- fact$args[[1L]]
      args <- c(args, list(list(fnode("sin", list(a0)), is_numer), list(fnode("sin", list(.sub(.div(S_Pi, S_Two), a0))), !is_numer)))
    } else if (fact$type == "Function" && fact$name == "cos") {
      a0 <- fact$args[[1L]]
      args <- c(args, list(list(fnode("sin", list(.sub(.div(S_Pi, S_Two), a0))), is_numer)))
    } else if (fact$type == "Function" && fact$name == "cot") {
      a0 <- fact$args[[1L]]
      args <- c(args, list(list(fnode("sin", list(.sub(.div(S_Pi, S_Two), a0))), is_numer), list(fnode("sin", list(a0)), !is_numer)))
    } else exception(fact)
  }
  fac <- .mul(list(fac, .div(.mul(facs), .mul(dfacs))))
  an <- list(); ap <- list(); bm <- list(); bq <- list()
  for (is_numer in c(TRUE, FALSE)) {
    gammas <- if (is_numer) numer_gammas else denom_gammas
    while (length(gammas)) {
      ac <- gammas[[length(gammas)]]; gammas[[length(gammas)]] <- NULL
      A <- ac[[1L]]; c <- ac[[2L]]
      if (!is_one_num(A) && !identical(A$key, S_NegativeOne$key)) {
        if (A$type != "Integer") stop("a is not an integer")
        p <- Abs(A); newa <- .div(A, p); newc <- .div(c, p)
        for (k in seq_len(.int_val(p)) - 1L) gammas[[length(gammas) + 1L]] <- list(newa, .add(list(newc, .div(.int(k), p))))
        gf <- .mul(list(.pow(.mul(list(S_Two, S_Pi)), .div(.sub(S_One, p), S_Two)), .pow(p, .sub(c, S_Half))))
        if (is_numer) { fac <- .mul(list(fac, gf)); exponentials[[length(exponentials) + 1L]] <- .pow(p, A) }
        else { fac <- .div(fac, gf); exponentials[[length(exponentials) + 1L]] <- .pow(p, .neg(A)) }
        next
      }
      if (is_one_num(A)) { if (is_numer) an[[length(an) + 1L]] <- .sub(S_One, c) else bq[[length(bq) + 1L]] <- .sub(S_One, c) }
      else { if (is_numer) bm[[length(bm) + 1L]] <- c else ap[[length(ap) + 1L]] <- c }
    }
  }
  srt <- function(l) if (length(l) > 1L) .key_sort(lapply(l, .sort_key), l) else l
  list(list(srt(an), srt(ap)), list(srt(bm), srt(bq)), .mul(exponentials), exponent, fac)
}

# roots of a polynomial in s with multiplicity: roots(), else CRootOf.all_roots
.rg_roots <- function(p, s) {
  rs <- roots(p, s)
  if (length(rs) != degree(p, s)) return(all_roots(p, s))
  unclass(rs)[seq_along(rs)]
}

# _inverse_mellin_transform: list(value, cond)
.inverse_mellin_core <- function(F, s, x_, strip, as_meijerg = FALSE) {
  x <- .mi_dummy("t", "inverse-mellin-transform", F, positive = TRUE)
  F <- rewrite(F, "gamma")
  for (g in list(sym_factor(F), expand_mul(F), expand(F))) {
    if (g$type == "Add") {
      ress <- lapply(.sympy_args(g), function(G) .inverse_mellin_core(G, s, x, strip, as_meijerg))
      res <- .add(lapply(ress, `[[`, 1L))
      if (!as_meijerg) res <- sym_factor(res)
      return(list(subs(res, x, x_), do.call(And, lapply(ress, `[[`, 2L))))
    }
    rg <- tryCatch(.mellin_rewrite_gamma(g, s, strip[[1L]], strip[[2L]]), IntegralTransformError = function(e) NULL)
    if (is.null(rg)) next
    G <- tryCatch(meijerg(rg[[1L]], rg[[2L]], .div(rg[[3L]], .pow(x, rg[[4L]]))), error = function(e) NULL)
    if (is.null(G)) next
    if (as_meijerg) h <- G
    else {
      h <- hyperexpand(G)
      if (h$type == "Piecewise" && length(h$args) == 3L) {
        C <- rg[[3L]]
        h <- .add(list(.mul(list(Heaviside(.sub(x, Abs(C))), h$args[[1L]]$args[[1L]])),
                       .mul(list(Heaviside(.sub(Abs(C), x)), h$args[[2L]]$args[[1L]]))))
      }
    }
    za <- G$args[[3L]]
    lhs <- Abs(fcall("arg", list(za)))
    cond <- Or(Lt(lhs, .mul(list(.mg_delta(G), S_Pi))),
               And(Or(.mi_tb(length(.mg_ap(G)) != length(.mg_bq(G))), Ge(0, .add(list(fcall("re", list(.mg_nu(G))), S_One)))),
                   Eq(lhs, .mul(list(.mg_delta(G), S_Pi)))))
    if (identical(.mi_bool(cond), FALSE)) .it_error("Inverse Mellin", "does not converge")
    return(list(subs(.mul(list(h, rg[[5L]])), x, x_), cond))
  }
  .it_error("Inverse Mellin", "")
}

.IMT_ALLOWED <- c("exp", "gamma", "sin", "cos", "tan", "cot", "cosh", "sinh", "tanh", "coth", "factorial", "RisingFactorial")

.inverse_mellin_doit <- function(F, s, x, strip, as_meijerg, needeval, noconds, simplify) {
  T <- NULL
  if (!.applied_undef_in(F, s)) {
    T <- tryCatch({
      for (g in atoms(F, "Function")) if (.has1(g, s) && !(g$name %in% .IMT_ALLOWED))
        .it_error("Inverse Mellin", paste0("Component ", .str(g), " not recognised"))
      r <- .inverse_mellin_core(F, s, x, strip, as_meijerg)
      if (noconds) r[[1L]] else r
    }, IntegralTransformError = function(e) NULL)
  }
  fn <- if (F$type == "Add") F else expand_mul(F)
  if (!is.null(T)) return(T)
  if (fn$type == "Add") {
    res <- lapply(.sympy_args(fn), function(a) .inverse_mellin_doit(a, s, x, strip, as_meijerg, needeval, noconds, simplify))
    ress <- list(); extra <- list()
    for (r in res) {
      if (is_sym(r)) { ress[[length(ress) + 1L]] <- r; next }
      ress[[length(ress) + 1L]] <- r[[1L]]
      extra[[length(extra) + 1L]] <- r[[2L]]
    }
    tot <- if (isTRUE(simplify)) simplify(.add(ress)) else .add(ress)
    if (!length(extra)) return(tot)
    cc <- do.call(And, extra)
    if (!identical(cc$key, S_false$key)) return(list(tot, cc))
  }
  if (needeval) .it_error("Inverse Mellin", "needeval")
  cm <- .as_coeff_mul_dep(fn, s)
  none <- .imt_none()
  .mul(list(cm[[1L]], InverseMellinTransform(.mul(cm[[2L]]), s, x, if (is.null(strip[[1L]])) none else strip[[1L]],
                                               if (is.null(strip[[2L]])) none else strip[[2L]])))
}
.imt_none <- function() { if (is.null(.mi_state$imt_none)) .mi_state$imt_none <- Dummy("None"); .mi_state$imt_none }

#' @rdname mellin_transform
#' @export
inverse_mellin_transform <- function(F, s, x, strip, as_meijerg = FALSE, needeval = FALSE, noconds = TRUE, simplify = TRUE) {
  F <- as_sym(F); s <- as_sym(s); x <- as_sym(x)
  strip <- lapply(strip, function(z) if (is.null(z)) NULL else if (is.numeric(z) && is.infinite(z)) (if (z > 0) S_Infinity else S_NegativeInfinity) else as_sym(z))
  if (length(strip) < 2L) strip <- c(strip, rep(list(NULL), 2L - length(strip)))
  .inverse_mellin_doit(F, s, x, strip, as_meijerg, needeval, noconds, simplify)
}
