# Laplace transforms, after sympy/integrals/laplace.py: the rule tables,
# the programmatic rules, integration with condition processing, the
# inverse transform with its rational, irrational and shift rules.

.lt_m <- function(e, pat) {
  m <- tryCatch(sym_match(e, pat), error = function(err) NULL)
  if (is.null(m) || !length(m)) NULL else m
}
.lt_pairs <- function(m) {
  if (is.null(m)) return(list())
  lapply(names(m), function(nm) list(attr(m, "wilds")[[nm]] %||% Wild(nm), m[[nm]]))
}
# xreplace of the wilds of a match result
.lt_xw <- function(e, m, wilds) {
  prs <- lapply(names(wilds), function(nm) if (!is.null(m[[nm]])) list(wilds[[nm]], m[[nm]]) else NULL)
  prs <- Filter(Negate(is.null), prs)
  if (!length(prs)) return(e)
  xreplace(e, prs)
}
.lt_pos <- function(e) isTRUE(ask_fact(e, "positive"))
.lt_neg <- function(e) isTRUE(ask_fact(e, "negative"))
.lt_re <- function(e) fcall("re", list(e))
.lt_im <- function(e) fcall("im", list(e))
.lt_H <- function(e) fcall("Heaviside", list(e))
.lt_D <- function(e, n = NULL) fcall("DiracDelta", if (is.null(n)) list(e) else list(e, n))
.lt_exp <- function(e) fcall("exp", list(e))
.lt_max <- function(l) if (length(l) == 1L) l[[1L]] else do.call(Max, l)
.lt_and <- function(l) if (!length(l)) S_true else if (length(l) == 1L) l[[1L]] else do.call(And, l)

# _canonical: every relational in canonical form
.lt_canonical <- function(cond) .map_bottom_up(cond, function(z) if (z$type %in% .REL_TYPES) .rel_canonical(z) else z)

# _simplifyconds
.lt_simplifyconds <- function(expr, s, a) {
  power <- function(ex) {
    if (identical(ex$key, s$key)) return(S_One)
    if (ex$type == "Pow" && identical(ex$args[[1L]]$key, s$key)) return(ex$args[[2L]])
    NULL
  }
  bigger <- function(ex1, ex2) {
    if (.has1(ex1, s) && .has1(ex2, s)) return(NA)
    if (ex1$type == "Function" && ex1$name == "Abs") ex1 <- ex1$args[[1L]]
    if (ex2$type == "Function" && ex2$name == "Abs") ex2 <- ex2$args[[1L]]
    if (.has1(ex1, s)) return(bigger(.inv(ex2), .inv(ex1)))
    n <- power(ex2)
    if (is.null(n)) return(NA)
    if (.is_true(Gt(n, 0)) && .is_true(Le(Abs(ex1), .pow(Abs(a), n)))) return(FALSE)
    if (.is_true(Lt(n, 0)) && .is_true(Ge(Abs(ex1), .pow(Abs(a), n)))) return(TRUE)
    NA
  }
  isabs <- function(z) z$type == "Function" && z$name == "Abs"
  replie <- function(x, y) {
    if (!(.lt_pos(x) || isabs(x)) || !(.lt_pos(y) || isabs(y))) return(Lt(x, y))
    r <- bigger(x, y)
    if (!is.na(r)) return(if (r) S_false else S_true)
    Lt(x, y)
  }
  replue <- function(x, y) {
    b <- bigger(x, y)
    if (!is.na(b)) return(S_true)
    Ne(x, y)
  }
  if (identical(expr$key, S_true$key) || identical(expr$key, S_false$key)) return(expr)
  expr <- collect_abs(expr)
  expr <- .map_bottom_up(expr, function(z) if (z$type == "StrictLessThan") replie(z$args[[1L]], z$args[[2L]]) else z)
  expr <- .map_bottom_up(expr, function(z) if (z$type == "StrictGreaterThan") replie(z$args[[2L]], z$args[[1L]]) else z)
  .map_bottom_up(expr, function(z) if (z$type == "Unequality") replue(z$args[[1L]], z$args[[2L]]) else z)
}

#' @rdname laplace_transform
#' @export
expand_dirac_delta <- function(expr) {
  expr <- as_sym(expr)
  ds <- .uniq_syms(Filter(function(z) z$type == "Function" && z$name == "DiracDelta", .all_subexprs(expr)))
  e <- expand(expr)
  coeffs <- list(); rest <- list()
  for (term in if (e$type == "Add") e$args else list(e)) {
    hit <- NULL
    for (d in ds) {
      q <- cancel(.div(term, d))
      if (!.has1(q, d) && !any(vapply(ds, function(z) .has1(q, z), NA))) { hit <- list(d, q); break }
    }
    if (is.null(hit)) rest[[length(rest) + 1L]] <- term
    else { k <- hit[[1L]]$key; coeffs[[k]] <- list(hit[[1L]], if (is.null(coeffs[[k]])) hit[[2L]] else .add(list(coeffs[[k]][[2L]], hit[[2L]]))) }
  }
  list(.add(rest), coeffs)
}

# _laplace_transform_integration
.lt_integration <- function(f, t, s_, simplify) {
  s <- Dummy("s")
  if (.has_function(f, "DiracDelta")) return(NULL)
  F <- integrate(.mul(list(f, .lt_exp(.neg(.mul(list(s, t)))))), list(t, S_Zero, S_Infinity))
  if (!.has_type(F, "Integral")) return(list(.tr_simplify(subs(F, s, s_), simplify), S_NegativeInfinity, S_true))
  if (F$type != "Piecewise") return(NULL)
  cond <- F$args[[1L]]$args[[2L]]; F <- F$args[[1L]]$args[[1L]]
  if (.has_type(F, "Integral")) return(NULL)
  process_conds <- function(conds) {
    a <- S_NegativeInfinity; aux <- S_true
    conds <- .conjuncts(to_cnf(conds))
    W <- function(nm) Wild(nm, exclude = list(s))
    p <- W("p"); q <- W("q"); w1 <- W("w1"); w2 <- W("w2"); w3 <- W("w3"); w4 <- W("w4"); w5 <- W("w5")
    fA <- function(z) fcall("Abs", list(z))
    pa <- function(z, per) fcall("periodic_argument", list(z, per))
    patterns <- list(Lt(p * fA(fcall("arg", list((s + w3) * q))), w2), Le(p * fA(fcall("arg", list((s + w3) * q))), w2),
                     Lt(fA(pa((s + w3)^p * q, w1)), w2), Le(fA(pa((s + w3)^p * q, w1)), w2),
                     Lt(fA(pa(fcall("polar_lift", list(s + w3))^p * q, w1)), w2), Le(fA(pa(fcall("polar_lift", list(s + w3))^p * q, w1)), w2))
    for (c in conds) {
      a_ <- S_Infinity; aux_ <- list()
      for (d in .disjuncts(c)) {
        isrel <- d$type %in% .REL_TYPES
        if (isrel && .has1(d$args[[2L]], s)) d <- .rel_reversed(d)
        if (d$type %in% c("GreaterThan", "StrictGreaterThan")) d <- .rel_reversedsign(d)
        m <- NULL
        for (pat in patterns) { m <- .lt_m(d, pat); if (!is.null(m)) break }
        if (!is.null(m) && .lt_pos(m$q) && .mi_eq(.div(m$w2, m$p), .div(S_Pi, S_Two)))
          d <- Lt(.neg(.lt_re(.add(list(s, m$w3)))), 0)
        m <- .lt_m(d, Lt(p - cos(w1 * fA(fcall("arg", list(s * w5))) * w2) * fA(s^w3)^w4, 0))
        if (is.null(m)) m <- .lt_m(d, Lt(cos(p - fA(pa(s^w1 * w5, q)) * w2) * fA(s^w3)^w4, 0))
        if (is.null(m)) m <- .lt_m(d, Lt(p - cos(fA(pa(fcall("polar_lift", list(s))^w1 * w5, q)) * w2) * fA(s^w3)^w4, 0))
        if (!is.null(m) && all(vapply(c("w1", "w2", "w3", "w4", "w5"), function(w) .lt_pos(m[[w]]), NA))) d <- Gt(.lt_re(s), m$p)
        d_ <- subs(.map_bottom_up(d, function(z) if (z$type == "Function" && z$name == "re") .as_real_imag(expand(z$args[[1L]]))[[1L]] else z), .lt_re(s), t)
        isrel <- d$type %in% .REL_TYPES
        if (!isrel || d$type %in% c("Equality", "Unequality") || .has1(d_, s) || !.has1(d_, t)) { aux_[[length(aux_) + 1L]] <- d; next }
        soln <- tryCatch(.solve_inequality(d_, t), error = function(e) NULL)
        if (is.null(soln) || !(soln$type %in% .REL_TYPES) || soln$type %in% c("Equality", "Unequality")) { aux_[[length(aux_) + 1L]] <- d; next }
        if (identical(.rel_lts(soln)$key, t$key)) return(NULL)
        a_ <- Min(.rel_lts(soln), a_)
      }
      if (!identical(a_$key, S_Infinity$key)) a <- Max(a_, a)
      else aux <- And(aux, if (length(aux_) == 1L) aux_[[1L]] else do.call(Or, aux_))
    }
    list(a, if (aux$type %in% .REL_TYPES) .rel_canonical(aux) else aux)
  }
  conds <- lapply(.disjuncts(cond), process_conds)
  if (any(vapply(conds, is.null, NA))) return(NULL)
  conds2 <- Filter(function(x) !identical(x[[2L]]$key, S_false$key) && !identical(x[[1L]]$key, S_NegativeInfinity$key), conds)
  if (!length(conds2)) conds2 <- Filter(function(x) !identical(x[[2L]]$key, S_false$key), conds)
  if (!length(conds2)) return(NULL)
  cnt <- function(e) if (e$type == "BooleanAtom") 0 else as.numeric(count_ops(e))
  # list(ordered(conds2)) then a stable sort by (-a, count_ops(aux))
  conds2 <- conds2[order(vapply(conds2, function(x) paste(.str(x[[1L]]), .str(x[[2L]])), ""))]
  key1 <- vapply(conds2, function(x) { v <- tryCatch(num_value_dbl(x[[1L]]), error = function(e) NA_real_); if (is.na(v)) 0 else -v }, 0)
  key2 <- vapply(conds2, function(x) cnt(x[[2L]]), 0)
  conds2 <- conds2[order(key1, key2)]
  a <- conds2[[1L]][[1L]]; aux <- conds2[[1L]][[2L]]
  if (simplify) { F <- .lt_simplifyconds(F, s, a); aux <- .lt_simplifyconds(aux, s, a) }
  list(.tr_simplify(subs(F, s, s_), simplify), subs(a, s, s_), .lt_canonical(subs(aux, s, s_)))
}

# _laplace_deep_collect
.lt_deep_collect <- function(f, t) {
  if (!is_sym(f) || f$type %in% c(.REL_TYPES, "And", "Or", "Not", "BooleanAtom")) return(f)
  p <- P_or_null(f, t)
  if (!is.null(p)) return(P_expr(p))
  if (!length(f$args)) return(f)
  sym_rebuild(f, lapply(f$args, .lt_deep_collect, t = t))
}

.lt_state <- new.env()
.lt_build_rules <- function() {
  if (!is.null(.lt_state$rules)) return(.lt_state$rules)
  t <- Dummy("t"); s <- Dummy("s")
  W <- function(nm) Wild(nm, exclude = list(t))
  a <- W("a"); b <- W("b"); n <- W("n"); tau <- W("tau"); omega <- W("omega")
  wilds <- list(a = a, b = b, n = n, tau = tau, omega = omega)
  E <- function(z) fcall("exp", list(as_sym(z))); f1 <- function(nm, ...) fcall(nm, lapply(list(...), as_sym))
  A <- function(z) fcall("Abs", list(z)); R <- .lt_re; I_ <- .lt_im; arg <- function(z) fcall("arg", list(z))
  Ei <- function(z) f1("Ei", z); erfc_ <- function(z) f1("erfc", z); sq <- function(z) sqrt(z)
  gam <- function(z) f1("gamma", z); EG <- S_EulerGamma
  r <- function(td, sd, chk, pl) list(td, sd, chk, pl)
  pi_ <- S_Pi; half <- S_Half
  rules <- list(
    r(a, a / s, S_true, S_Zero),
    r(.lt_D(a * t - b), E(-s * b / a) / A(a), Or(And(Gt(a, 0), Ge(b, 0)), And(Lt(a, 0), Le(b, 0))), S_NegativeInfinity),
    r(.lt_D(a * t - b), S_Zero, Or(And(Lt(a, 0), Ge(b, 0)), And(Gt(a, 0), Le(b, 0))), S_NegativeInfinity),
    r(.lt_H(a * t - b), E(-s * b / a) / s, And(Gt(a, 0), Gt(b, 0)), S_Zero),
    r(.lt_H(a * t - b), (1 - E(-s * b / a)) / s, And(Lt(a, 0), Lt(b, 0)), S_Zero),
    r(.lt_H(a * t - b), 1 / s, And(Gt(a, 0), Le(b, 0)), S_Zero),
    r(.lt_H(a * t - b), S_Zero, And(Lt(a, 0), Gt(b, 0)), S_Zero),
    r(t, 1 / s^2, S_true, S_Zero),
    r(1 / (a * t + b), -E(-b / a * s) * Ei(-b / a * s) / a, Lt(A(arg(b / a)), pi_), S_Zero),
    r(1 / sq(a * t + b), sq(a * pi_ / s) * E(b / a * s) * erfc_(sq(b / a * s)) / a, Lt(A(arg(b / a)), pi_), S_Zero),
    r((a * t + b)^Rational(-3, 2), 2 * b^Rational(-1, 2) - 2 * (pi_ * s / a)^half * E(b / a * s) * erfc_(sq(b / a * s)) / a, Lt(A(arg(b / a)), pi_), S_Zero),
    r(sq(t) / (t + b), sq(pi_ / s) - pi_ * sq(b) * E(b * s) * erfc_(sq(b * s)), Lt(A(arg(b)), pi_), S_Zero),
    r(1 / (a * sq(t) + t^Float(1.5)), pi_ * a^half * E(a * s) * erfc_(sq(a * s)), S_true, S_Zero),
    r(t^n, gam(n + 1) / s^(n + 1), Gt(n, -1), S_Zero),
    r((a * t + b)^n, f1("uppergamma", n + 1, b / a * s) * E(b / a * s) / s^(n + 1) / a, And(Gt(n, -1), Lt(A(arg(b / a)), pi_)), S_Zero),
    r(t^n / (t + a), a^n * gam(n + 1) * f1("uppergamma", -n, a * s), And(Gt(n, -1), Lt(A(arg(a)), pi_)), S_Zero),
    r(E(a * t - tau), E(-tau) / (s - a), S_true, R(a)),
    r(t * E(a * t - tau), E(-tau) / (s - a)^2, S_true, R(a)),
    r(t^n * E(a * t), gam(n + 1) / (s - a)^(n + 1), Gt(R(n), -1), R(a)),
    r(E(-a * t^2), sq(pi_ / 4 / a) * E(s^2 / 4 / a) * erfc_(s / sq(4 * a)), Gt(R(a), 0), S_Zero),
    r(t * E(-a * t^2), 1 / (2 * a) - 2 / sq(pi_) / (4 * a)^Rational(3, 2) * s * erfc_(s / sq(4 * a)), Gt(R(a), 0), S_Zero),
    r(E(-a / t), 2 * sq(a / s) * f1("besselk", 1, 2 * sq(a * s)), Ge(R(a), 0), S_Zero),
    r(sq(t) * E(-a / t), half * sq(pi_ / s^3) * (1 + 2 * sq(a * s)) * E(-2 * sq(a * s)), Ge(R(a), 0), S_Zero),
    r(E(-a / t) / sq(t), sq(pi_ / s) * E(-2 * sq(a * s)), Ge(R(a), 0), S_Zero),
    r(E(-a / t) / (t * sq(t)), sq(pi_ / a) * E(-2 * sq(a * s)), Gt(R(a), 0), S_Zero),
    r(t^n * E(-a / t), 2 * (a / s)^((n + 1) / 2) * f1("besselk", n + 1, 2 * sq(a * s)), Gt(R(a), 0), S_Zero),
    r(E(-a * E(-t)), a^(-s) * f1("lowergamma", s, a), S_true, S_Zero),
    r(E(-a * E(t)), a^s * f1("uppergamma", -s, a), Gt(R(a), 0), S_Zero),
    r(log(a * t), -log(E(EG) * s / a) / s, Gt(a, 0), S_Zero),
    r(log(1 + a * t), -E(s / a) / s * Ei(-s / a), Lt(A(arg(a)), pi_), S_Zero),
    r(log(a * t + b), (log(b) - E(s / b / a) / s * a * Ei(-s / b)) / s * a, And(Gt(a, 0), Lt(A(arg(b)), pi_)), S_Zero),
    r(log(t) / sq(t), -sq(pi_ / s) * log(4 * s * E(EG)), S_true, S_Zero),
    r(t^n * log(t), gam(n + 1) * s^(-n - 1) * (f1("polygamma", 0, n + 1) - log(s)), Gt(R(n), -1), S_Zero),
    r(log(a * t)^2, (log(E(EG) * s / a)^2 + pi_^2 / 6) / s, Gt(a, 0), S_Zero),
    r(sin(omega * t), omega / (s^2 + omega^2), S_true, A(I_(omega))),
    r(A(sin(omega * t)), omega / (s^2 + omega^2) * f1("coth", pi_ * s / 2 / omega), Gt(omega, 0), S_Zero),
    r(sin(omega * t) / t, atan(omega / s), S_true, A(I_(omega))),
    r(sin(omega * t)^2 / t, log(1 + 4 * omega^2 / s^2) / 4, S_true, 2 * A(I_(omega))),
    r(sin(omega * t)^2 / t^2, omega * atan(2 * omega / s) - s * log(1 + 4 * omega^2 / s^2) / 4, S_true, 2 * A(I_(omega))),
    r(cos(omega * t), s / (s^2 + omega^2), S_true, A(I_(omega))),
    r(cos(omega * t)^2, (s^2 + 2 * omega^2) / (s^2 + 4 * omega^2) / s, S_true, 2 * A(I_(omega))),
    r(sin(a * t) * sin(b * t), 2 * a * b * s / (s^2 + (a + b)^2) / (s^2 + (a - b)^2), S_true, A(I_(a)) + A(I_(b))),
    r(cos(a * t) * sin(b * t), b * (s^2 - a^2 + b^2) / (s^2 + (a + b)^2) / (s^2 + (a - b)^2), S_true, A(I_(a)) + A(I_(b))),
    r(cos(a * t) * cos(b * t), s * (s^2 + a^2 + b^2) / (s^2 + (a + b)^2) / (s^2 + (a - b)^2), S_true, A(I_(a)) + A(I_(b))),
    r(sinh(a * t), a / (s^2 - a^2), S_true, A(R(a))),
    r(cosh(a * t), s / (s^2 - a^2), S_true, A(R(a))),
    r(sinh(a * t)^2, 2 * a^2 / (s^3 - 4 * a^2 * s), S_true, 2 * A(R(a))),
    r(cosh(a * t)^2, (s^2 - 2 * a^2) / (s^3 - 4 * a^2 * s), S_true, 2 * A(R(a))),
    r(sinh(a * t) / t, log((s + a) / (s - a)) / 2, S_true, A(R(a))),
    r(t^n * sinh(a * t), gam(n + 1) / 2 * ((s - a)^(-n - 1) - (s + a)^(-n - 1)), Gt(n, -2), A(a)),
    r(t^n * cosh(a * t), gam(n + 1) / 2 * ((s - a)^(-n - 1) + (s + a)^(-n - 1)), Gt(n, -1), A(a)),
    r(f1("erf", a * t), E(s^2 / (2 * a)^2) * erfc_(s / (2 * a)) / s, Lt(4 * A(arg(a)), pi_), S_Zero),
    r(f1("besselj", n, a * t), a^n / (sq(s^2 + a^2) * (s + sq(s^2 + a^2))^n), Gt(R(n), -1), A(I_(a))),
    r(t^b * f1("besselj", n, a * t), 2^n / sq(pi_) * gam(n + half) * a^n * (s^2 + a^2)^(-n - half), And(Gt(R(n), -half), Eq(b, n)), A(I_(a))),
    r(t^b * f1("besselj", n, a * t), 2^(n + 1) / sq(pi_) * gam(n + Rational(3, 2)) * a^n * s * (s^2 + a^2)^(-n - Rational(3, 2)), And(Gt(R(n), -1), Eq(b, n + 1)), A(I_(a))),
    r(f1("besselj", 0, a * sq(t^2 + b * t)), E(b * s - b * sq(s^2 + a^2)) / sq(s^2 + a^2), Lt(A(arg(b)), pi_), A(I_(a))),
    r(f1("besseli", n, a * t), a^n / (sq(s^2 - a^2) * (s + sq(s^2 - a^2))^n), Gt(R(n), -1), A(R(a))),
    r(t^b * f1("besseli", n, a * t), 2^n / sq(pi_) * gam(n + half) * a^n * (s^2 - a^2)^(-n - half), And(Gt(R(n), -half), Eq(b, n)), A(R(a))),
    r(t^b * f1("besseli", n, a * t), 2^(n + 1) / sq(pi_) * gam(n + Rational(3, 2)) * a^n * s * (s^2 - a^2)^(-n - Rational(3, 2)), And(Gt(R(n), -1), Eq(b, n + 1)), A(R(a))),
    r(f1("bessely", 0, a * t), -2 / pi_ * asinh(s / a) / sq(s^2 + a^2), S_true, A(I_(a))),
    r(f1("besselk", 0, a * t), log((s + sq(s^2 - a^2)) / a) / sq(s^2 - a^2), S_true, -R(a)))
  .lt_state$rules <- list(rules = rules, t = t, s = s, wilds = wilds)
  .lt_state$rules
}

.lt_apply_simple_rules <- function(f, t, s) {
  rr <- .lt_build_rules()
  prep_f <- .lt_deep_collect(subs(f, t, rr$t), rr$t)
  for (rule in rr$rules) {
    ma <- .lt_m(prep_f, rule[[1L]])
    if (is.null(ma)) next
    c <- tryCatch(.lt_xw(rule[[3L]], ma, rr$wilds), error = function(e) NULL)
    if (is.null(c)) next
    if (.is_true(c)) return(list(subs(.lt_xw(rule[[2L]], ma, rr$wilds), rr$s, s), .lt_xw(rule[[4L]], ma, rr$wilds), S_true))
  }
  NULL
}

.lt_rule_timescale <- function(f, t, s) {
  a <- Wild("a", exclude = list(t))
  if (f$type == "Function" && length(f$args) == 1L) {
    arg <- collect(f$args[[1L]], t)
    ma2 <- .lt_m(arg, a * t)
    if (!is.null(ma2) && .lt_pos(ma2$a) && !is_one_num(ma2$a)) {
      r <- .lt_transform(.mul(list(.inv(ma2$a), fcall(f$name, list(t)))), t, .div(s, ma2$a), simplify = FALSE)
      return(r)
    }
  }
  NULL
}

.lt_rule_heaviside <- function(f, t, s) {
  a <- Wild("a", exclude = list(t)); y <- Wild("y"); g <- Wild("g")
  ma1 <- .lt_m(f, .lt_H(y) * g)
  if (!is.null(ma1)) {
    ma2 <- .lt_m(ma1$y, t - a)
    if (!is.null(ma2)) {
      if (.lt_pos(ma2$a)) {
        r <- .lt_transform(subs(ma1$g, t, .add(list(t, ma2$a))), t, s, simplify = FALSE)
        return(list(.mul(list(.lt_exp(.neg(.mul(list(ma2$a, s)))), r[[1L]])), r[[2L]], r[[3L]]))
      }
      if (.lt_neg(ma2$a)) return(.lt_transform(ma1$g, t, s, simplify = FALSE))
    }
    ma2 <- .lt_m(ma1$y, a - t)
    if (!is.null(ma2)) {
      if (.lt_pos(ma2$a)) return(.lt_transform(.mul(list(.sub(S_One, .lt_H(.sub(t, ma2$a))), ma1$g)), t, s, simplify = FALSE))
      if (.lt_neg(ma2$a)) return(list(S_Zero, S_Zero, S_true))
    }
  }
  NULL
}

.lt_rule_exp <- function(f, t, s) {
  a <- Wild("a", exclude = list(t)); y <- Wild("y"); z <- Wild("z")
  ma1 <- .lt_m(f, fnode("exp", list(y)) * z)
  if (!is.null(ma1)) {
    ma2 <- .lt_m(collect(ma1$y, t), a * t)
    if (!is.null(ma2)) {
      r <- .lt_transform(ma1$z, t, .sub(s, ma2$a), simplify = FALSE)
      return(list(r[[1L]], .add(list(r[[2L]], .lt_re(ma2$a))), r[[3L]]))
    }
  }
  NULL
}

.lt_rule_delta <- function(f, t, s) {
  a <- Wild("a", exclude = list(t)); b <- Wild("b", exclude = list(t)); y <- Wild("y"); z <- Wild("z")
  ma1 <- .lt_m(f, fnode("DiracDelta", list(y)) * z)
  if (!is.null(ma1) && !.has_function(ma1$z, "DiracDelta")) {
    ma2 <- .lt_m(collect(ma1$y, t), b * t - a)
    if (!is.null(ma2)) {
      loc <- .div(ma2$a, ma2$b)
      if (.is_true(Ge(.lt_re(loc), 0)) && is_zero_num(.lt_im(loc))) {
        fn <- .mul(list(.lt_exp(.neg(.mul(list(loc, s)))), ma1$z))
        if (.has_function(fn, c("sin", "cos"))) fn <- ratsimp(rewrite(fn, "sinc"))
        nd <- lapply(.as_numer_denom(fn), function(x) subs(x, t, loc))
        if (!is_zero_num(nd[[2L]])) return(list(.div(.div(nd[[1L]], nd[[2L]]), ma2$b), S_NegativeInfinity, S_true))
        return(NULL)
      }
      return(list(S_Zero, S_NegativeInfinity, S_true))
    }
    if (isTRUE(is_polynomial(ma1$y, t))) {
      ro <- tryCatch(roots(ma1$y, t), error = function(e) NULL)
      if (!is.null(ro) && length(ro) && all(attr(ro, "multiplicity") == 1L)) {
        slope <- .diff1(ma1$y, t)
        terms <- list()
        for (x in ro) if (is_zero_num(.lt_im(x)) && .is_true(Ge(.lt_re(x), 0)))
          terms[[length(terms) + 1L]] <- .div(.mul(list(.lt_exp(.neg(.mul(list(x, s)))), subs(ma1$z, t, s))), subs(slope, t, x))
        return(list(.add(terms), S_NegativeInfinity, S_true))
      }
    }
  }
  NULL
}

.lt_trig_split <- function(fn) {
  trigs <- list(S_One); other <- list(S_One)
  for (term in .margs(fn)) {
    if (.has_function(term, c("sin", "cos", "sinh", "cosh", "exp"))) trigs[[length(trigs) + 1L]] <- term
    else other[[length(other) + 1L]] <- term
  }
  list(.mul(trigs), .mul(other))
}

.lt_trig_expsum <- function(f, t) {
  c1 <- Wild("c1", exclude = list(t)); c0 <- Wild("c0", exclude = list(t)); p <- Wild("p", exclude = list(t))
  xm <- list(); xn <- list()
  x1 <- expand(rewrite(f, "exp"))
  for (term in if (x1$type == "Add") x1$args else list(x1)) {
    if (!.has1(term, t)) { xm[[length(xm) + 1L]] <- list(k = term, a = S_Zero, re = S_Zero, im = S_Zero); next }
    term <- .lt_deep_collect(powsimp(term, combine = "exp"), t)
    r <- .lt_m(term, p * fnode("exp", list(c1 * t + c0)))
    if (!is.null(r)) xm[[length(xm) + 1L]] <- list(k = .mul(list(r$p, .lt_exp(r$c0))), a = r$c1, re = .lt_re(r$c1), im = .lt_im(r$c1))
    else xn[[length(xn) + 1L]] <- term
  }
  list(xm, xn)
}

.lt_trig_ltex <- function(xm, t, s) {
  results <- list(); planes <- list()
  simpc <- function(coeffs) lapply(coeffs, function(c) {
    ri <- .as_real_imag(c)
    if (.has_function(ri[[1L]], "im")) rewrite(c, "cos") else rewrite(.add(list(ri[[1L]], .mul(list(S_I, ri[[2L]])))), "cos")
  })
  poly_s <- function(cs) { n <- length(cs); .add(lapply(seq_len(n), function(i) .mul(list(cs[[i]], .pow(s, .int(n - i)))))) }
  I <- S_I
  quadpole <- function(t1, k1, k2, k3) {
    a <- t1$a; k0 <- t1$k; ar <- t1$re; ai <- t1$im
    nc <- list(k0 + k1 + k2 + k3,
               a * (k0 + k1 - k2 - k3) - 2 * I * ai * k1 + 2 * I * ai * k2,
               a^2 * (-k0 - k1 - k2 - k3) + a * (4 * I * ai * k0 + 4 * I * ai * k3) + 4 * ai^2 * k0 + 4 * ai^2 * k3,
               a^3 * (-k0 - k1 + k2 + k3) + a^2 * (4 * I * ai * k0 + 2 * I * ai * k1 - 2 * I * ai * k2 - 4 * I * ai * k3) + a * (4 * ai^2 * k0 - 4 * ai^2 * k3))
    dc <- list(S_One, S_Zero, 2 * ai^2 - 2 * ar^2, S_Zero, ai^4 + 2 * ai^2 * ar^2 + ar^4)
    .div(poly_s(simpc(nc)), poly_s(dc))
  }
  ccpole <- function(t1, k1) {
    a <- t1$a; k0 <- t1$k; ar <- t1$re; ai <- t1$im
    .div(poly_s(simpc(list(k0 + k1, -a * k0 - a * k1 + 2 * I * ai * k0))), poly_s(list(S_One, -2 * ar, ai^2 + ar^2)))
  }
  rspole <- function(t1, k2) {
    a <- t1$a; k0 <- t1$k; ar <- t1$re; ai <- t1$im
    .div(poly_s(simpc(list(k0 + k2, a * k0 - a * k2 - 2 * I * ai * k0))), poly_s(list(S_One, -2 * I * ai, -ai^2 - ar^2)))
  }
  sypole <- function(t1, k3) {
    a <- t1$a; k0 <- t1$k
    .div(poly_s(simpc(list(k0 + k3, a * (k0 - k3)))), poly_s(list(S_One, S_Zero, -a^2)))
  }
  eqk <- function(u, v) .mi_eq(u, v)
  while (length(xm)) {
    t1 <- xm[[length(xm)]]; xm <- xm[-length(xm)]
    ii <- NULL; ir <- NULL; ip <- NULL
    for (i in seq_along(xm)) {
      real_eq <- eqk(t1$re, xm[[i]]$re); realsym <- eqk(t1$re, .neg(xm[[i]]$re))
      imag_eq <- eqk(t1$im, xm[[i]]$im); imagsym <- eqk(t1$im, .neg(xm[[i]]$im))
      if (realsym && imagsym && !is_zero_num(t1$re) && !is_zero_num(t1$im)) ip <- i
      else if (realsym && imag_eq && !is_zero_num(t1$re)) ir <- i
      else if (real_eq && imagsym && !is_zero_num(t1$im)) ii <- i
    }
    if (!is.null(ii) && !is.null(ir) && !is.null(ip)) {
      results[[length(results) + 1L]] <- quadpole(t1, xm[[ii]]$k, xm[[ir]]$k, xm[[ip]]$k)
      planes[[length(planes) + 1L]] <- Abs(.lt_re(t1$a))
      xm <- xm[-c(ii, ir, ip)]
    } else if (!is.null(ii)) {
      results[[length(results) + 1L]] <- ccpole(t1, xm[[ii]]$k); planes[[length(planes) + 1L]] <- t1$re; xm <- xm[-ii]
    } else if (!is.null(ir)) {
      results[[length(results) + 1L]] <- rspole(t1, xm[[ir]]$k); planes[[length(planes) + 1L]] <- Abs(t1$re); xm <- xm[-ir]
    } else if (!is.null(ip)) {
      results[[length(results) + 1L]] <- sypole(t1, xm[[ip]]$k); planes[[length(planes) + 1L]] <- Abs(t1$re); xm <- xm[-ip]
    } else {
      results[[length(results) + 1L]] <- .div(t1$k, .sub(s, t1$a)); planes[[length(planes) + 1L]] <- t1$re
    }
  }
  list(.add(results), .lt_max(planes))
}

.lt_rule_trig <- function(fn, t_, s) {
  t <- Dummy("t", real = TRUE)
  if (!.has_function(fn, c("sin", "cos", "sinh", "cosh"))) return(NULL)
  fg <- .lt_trig_split(subs(fn, t_, t))
  es <- .lt_trig_expsum(fg[[1L]], t)
  xm <- es[[1L]]; xn <- es[[2L]]
  if (length(xn)) return(NULL)
  g <- fg[[2L]]
  if (!.has1(g, t)) {
    rp <- .lt_trig_ltex(xm, t, s)
    return(list(.mul(list(g, rp[[1L]])), rp[[2L]], S_true))
  }
  G <- .lt_transform(g, t, s, simplify = FALSE)
  results <- lapply(xm, function(x1) .mul(list(x1$k, subs(G[[1L]], s, .sub(s, x1$a)))))
  planes <- lapply(xm, function(x1) .add(list(G[[2L]], .lt_re(x1$a))))
  list(subs(.add(results), t, t_), .lt_max(planes), G[[3L]])
}

.lt_rule_diff <- function(f, t, s) {
  a <- Wild("a", exclude = list(t)); n <- Wild("n", exclude = list(t)); g <- WildFunction("g")
  if (!.has_type(f, "Derivative")) return(NULL)
  ma1 <- .lt_m(f, a * new_sym("Derivative", list(g, Tuple(t, n))))
  if (!is.null(ma1) && isTRUE(ask_fact(ma1$n, "integer"))) {
    m <- vapply(ma1$g$args, function(z) .has1(z, t), NA)
    if (sum(m) == 1L) {
      nn <- .int_val(ma1$n)
      d <- lapply(seq_len(nn) - 1L, function(k) {
        y <- if (k == 0L) subs(ma1$g, t, S_Zero) else subs(Derivative(ma1$g, list(t, k)), t, S_Zero)
        .mul(list(.pow(s, .int(nn - k - 1L)), y))
      })
      r <- .lt_transform(ma1$g, t, s, simplify = FALSE)
      return(list(.mul(list(ma1$a, .sub(.mul(list(.pow(s, ma1$n), r[[1L]])), .add(d)))), r[[2L]], r[[3L]]))
    }
  }
  NULL
}

.lt_rule_sdiff <- function(f, t, s) {
  if (f$type == "Mul") {
    pfac <- list(); ofac <- list()
    for (fac in .margs(f)) if (.mi_is_poly(fac, t)) pfac[[length(pfac) + 1L]] <- fac else ofac[[length(ofac) + 1L]] <- fac
    if (length(pfac)) {
      pc <- P_all_coeffs(.P(.mul(pfac), t))
      N <- length(pc)
      if (N > 1L) {
        rr <- .lt_transform(.mul(ofac), t, s, simplify = FALSE)
        r_ <- rr[[1L]]
        deri <- list(r_)
        d1 <- tryCatch(.neg(.diff1(r_, s)), error = function(e) NULL)
        if (.has_type(r_, "LaplaceTransform")) {
          for (k in seq_len(N - 1L)) deri[[length(deri) + 1L]] <- .mul(list(.int((-1L)^k), Derivative(r_, list(s, k))))
        } else if (!is.null(d1) && !is_zero_num(d1)) {
          deri[[length(deri) + 1L]] <- d1
          if (N > 2L) for (k in seq_len(N - 2L)) deri[[length(deri) + 1L]] <- .neg(.diff1(deri[[length(deri)]], s))
        }
        if (!is.null(d1) && !is_zero_num(d1))
          return(list(.add(lapply(seq_len(N) - 1L, function(n) .mul(list(pc[[N - n]], deri[[n + 1L]])))), rr[[2L]], rr[[3L]]))
      }
    }
  }
  n <- Wild("n", exclude = list(t)); g <- Wild("g")
  ma1 <- .lt_m(f, t^n * g)
  if (!is.null(ma1) && isTRUE(ask_fact(ma1$n, "integer")) && .lt_pos(ma1$n)) {
    rr <- .lt_transform(ma1$g, t, s, simplify = FALSE)
    dr <- if (ma1$n$type == "Integer") diff(rr[[1L]], s, .int_val(ma1$n)) else new_sym("Derivative", list(rr[[1L]], Tuple(s, ma1$n)))
    return(list(.mul(list(.pow(S_NegativeOne, ma1$n), dr)), rr[[2L]], rr[[3L]]))
  }
  NULL
}

.lt_expand <- function(f, t, s) {
  r <- expand(f, deep = FALSE)
  if (r$type == "Add") return(.lt_transform(r, t, s, simplify = FALSE))
  r <- expand_mul(f)
  if (r$type == "Add") return(.lt_transform(r, t, s, simplify = FALSE))
  r <- expand(f)
  if (r$type == "Add") return(.lt_transform(r, t, s, simplify = FALSE))
  if (!.mi_eq(r, f)) return(.lt_transform(r, t, s, simplify = FALSE))
  r <- expand(expand_trig(f))
  if (r$type == "Add") return(.lt_transform(r, t, s, simplify = FALSE))
  NULL
}

.lt_apply_prog_rules <- function(f, t, s) {
  for (rule in list(.lt_rule_heaviside, .lt_rule_delta, .lt_rule_timescale, .lt_rule_exp, .lt_rule_trig, .lt_rule_diff, .lt_rule_sdiff)) {
    L <- rule(f, t, s)
    if (!is.null(L)) return(L)
  }
  NULL
}

.lt_piecewise_to_heaviside <- function(f, t) {
  if (!isTRUE(ask_fact(t, "real"))) {
    r <- Dummy("r", real = TRUE)
    return(xreplace(.lt_piecewise_to_heaviside(xreplace(f, list(list(t, r))), r), list(list(r, t))))
  }
  x <- piecewise_exclusive(f)
  out <- list()
  isrel <- function(c) c$type %in% .REL_TYPES
  for (pr in x$args) {
    fn <- pr$args[[1L]]; cond <- pr$args[[2L]]
    if (isrel(cond) && any(vapply(cond$args, function(z) .mi_eq(z, t), NA))) {
      if (cond$type %in% c("Equality", "Unequality")) return(f)
      out[[length(out) + 1L]] <- .mul(list(.lt_H(.sub(.rel_gts(cond), .rel_lts(cond))), fn))
    } else if (cond$type == "Or" && length(cond$args) == 2L) {
      for (c2 in cond$args) {
        if (isrel(c2) && .mi_eq(c2$args[[1L]], t)) out[[length(out) + 1L]] <- .mul(list(.lt_H(.sub(.rel_gts(c2), .rel_lts(c2))), fn))
        else return(f)
      }
    } else if (cond$type == "And" && length(cond$args) == 2L) {
      c0 <- cond$args[[1L]]; c1 <- cond$args[[2L]]
      if (isrel(c0) && isrel(c1) && .mi_eq(c0$args[[1L]], t) && .mi_eq(c1$args[[1L]], t)) {
        if (grepl(">", .rel_op(c0), fixed = TRUE)) { tmp <- c0; c0 <- c1; c1 <- tmp }
        out[[length(out) + 1L]] <- .mul(list(.sub(.lt_H(.sub(.rel_gts(c1), .rel_lts(c1))), .lt_H(.sub(.rel_lts(c0), .rel_gts(c0)))), fn))
      } else return(f)
    } else return(f)
  }
  .add(out)
}

.lt_applied_undef_has <- function(f, t) any(vapply(atoms(f, "Function"), function(g) .is_undefined_application(g) && .has1(g, t), NA))

# _laplace_transform
.lt_transform <- function(fn, t_, s_, simplify) {
  terms <- list()
  for (ff in .mi_make_args_add(fn)) {
    sp <- .as_coeff_mul_dep(ff, t_); k <- sp[[1L]]; ft <- .mul(sp[[2L]])
    if (.has_function(ft, "SingularityFunction")) {
      for (tm in .mi_make_args_add(rewrite(ft, "Heaviside"))) { sp1 <- .as_coeff_mul_dep(tm, t_); terms[[length(terms) + 1L]] <- list(.mul(list(k, sp1[[1L]])), .mul(sp1[[2L]])) }
    } else if (ft$type == "Piecewise" && !.has1(ft, .lt_D(t_))) {
      for (tm in .mi_make_args_add(.lt_piecewise_to_heaviside(ft, t_))) { sp1 <- .as_coeff_mul_dep(tm, t_); terms[[length(terms) + 1L]] <- list(.mul(list(k, sp1[[1L]])), .mul(sp1[[2L]])) }
    } else terms[[length(terms) + 1L]] <- list(k, ft)
  }
  terms_s <- list(); planes <- list(); conditions <- list()
  for (kt in terms) {
    k <- kt[[1L]]; ft <- kt[[2L]]
    if (.has_function(ft, "SingularityFunction")) r <- list(LaplaceTransform(ft, t_, s_), S_NegativeInfinity, S_true)
    else {
      if (.has1(ft, .lt_H(t_)) && !.has1(ft, .lt_D(t_))) ft <- subs(ft, .lt_H(t_), S_One)
      r <- .lt_apply_simple_rules(ft, t_, s_)
      if (is.null(r)) r <- .lt_apply_prog_rules(ft, t_, s_)
      if (is.null(r)) r <- .lt_expand(ft, t_, s_)
      if (is.null(r)) {
        if (.lt_applied_undef_has(ft, t_)) r <- list(LaplaceTransform(ft, t_, s_), S_NegativeInfinity, S_true)
        else {
          r <- .lt_integration(ft, t_, s_, simplify)
          if (is.null(r)) r <- list(LaplaceTransform(ft, t_, s_), S_NegativeInfinity, S_true)
        }
      }
    }
    terms_s[[length(terms_s) + 1L]] <- .mul(list(k, r[[1L]]))
    planes[[length(planes) + 1L]] <- r[[2L]]
    conditions[[length(conditions) + 1L]] <- r[[3L]]
  }
  result <- .add(terms_s)
  if (simplify) result <- simplify(result, doit = FALSE)
  list(result, .lt_max(planes), .lt_and(conditions))
}

#' Laplace transforms
#'
#' Ports of SymPy's `laplace_transform` and `inverse_laplace_transform`
#' (sympy/integrals/laplace.py). `laplace_transform(f, t, s)` returns
#' `list(F, a, cond)`: the transform, its half plane of convergence
#' `re(s) > a` and an auxiliary condition; with `noconds = TRUE` only `F`.
#' `inverse_laplace_transform(F, s, t)` returns the inverse (by default
#' without conditions). `expand_dirac_delta()` splits an expression linear in
#' DiracDelta terms into the remainder and the coefficients.
#'
#' @param f,F,expr Expressions.
#' @param t,s Variables.
#' @param plane The half plane of convergence, or `NULL`.
#' @param noconds Return only the transform.
#' @param simplify Simplify the result.
#' @return The transform, or a list with the conditions.
#' @examples
#' t <- Symbol("t", positive = TRUE); s <- Symbol("s")
#' laplace_transform(exp(-2*t), t, s)
#' inverse_laplace_transform(1/(s + 2), s, t)
#' @export
laplace_transform <- function(f, t, s, noconds = FALSE, simplify = FALSE) {
  if (inherits(f, "symmatrix")) {
    els <- .mat_map(f, function(z) laplace_transform(z, t, s, noconds = TRUE, simplify = simplify))
    if (noconds) return(els)
    rs <- lapply(f$data, function(z) laplace_transform(z, t, s, noconds = FALSE, simplify = simplify))
    return(list(els, .lt_max(lapply(rs, `[[`, 2L)), .lt_and(lapply(rs, `[[`, 3L))))
  }
  f <- as_sym(f); t <- as_sym(t); s <- as_sym(s)
  r <- .lt_transform(f, t, s, simplify)
  if (noconds) r[[1L]] else r
}

# ---------------------------------------------------------------------------
# inverse

.ilt_integration <- function(F, s, t_, plane, simplify) {
  t <- Dummy("t", real = TRUE)
  pw_simp <- function(pw) {
    if (length(pw$args) != 3L) return(pw)
    arg <- pw$args[[3L]]$args[[1L]]$args[[1L]]$args[[1L]]
    ce <- .mi_coeff_exp(arg, t)
    e1 <- pw$args[[1L]]$args[[1L]]; e2 <- pw$args[[2L]]$args[[1L]]
    .add(list(.mul(list(.lt_H(.sub(.inv(Abs(ce[[1L]])), .pow(t, ce[[2L]]))), e1)), .mul(list(.lt_H(.sub(.pow(t, ce[[2L]]), .inv(Abs(ce[[1L]])))), e2))))
  }
  if (.mi_is_ratfun(F, s)) F <- apart(F, s)
  if (F$type == "Add") {
    f <- .add(lapply(F$args, function(X) { r <- .ilt_integration(X, s, t, plane, simplify); if (is.null(r)) stop("no inverse") else r[[1L]] }))
    return(list(.tr_simplify(subs(f, t, t_), simplify), S_true))
  }
  r <- tryCatch(inverse_mellin_transform(F, s, .lt_exp(.neg(t)), list(NULL, S_Infinity), needeval = TRUE, noconds = FALSE),
                IntegralTransformError = function(e) NULL)
  if (!is.null(r)) { f <- r[[1L]]; cond <- r[[2L]] }
  else {
    f <- meijerint_inversion(F, s, t)
    if (is.null(f)) return(NULL)
    if (f$type == "Piecewise") {
      cond <- f$args[[1L]]$args[[2L]]; f <- f$args[[1L]]$args[[1L]]
      if (.has_type(f, "Integral")) return(NULL)
    } else cond <- S_true
    f <- .map_bottom_up(f, function(z) if (z$type == "Piecewise") pw_simp(z) else z)
  }
  if (f$type == "Piecewise") return(list(subs(f, t, t_), cond))
  u <- Dummy("u")
  simp_heaviside <- function(z) {
    arg <- z$args[[1L]]
    a <- subs(arg, .lt_exp(.neg(t)), u)
    if (.has1(a, t)) return(z)
    rel <- .solve_inequality(Gt(a, 0), u)
    if (.mi_eq(.rel_lts(rel), u)) return(fcall("Heaviside", c(list(.add(list(t, log(.rel_gts(rel))))), z$args[-1L])))
    fcall("Heaviside", c(list(.neg(.add(list(t, log(.rel_lts(rel)))))), z$args[-1L]))
  }
  f <- .map_bottom_up(f, function(z) if (z$type == "Function" && z$name == "Heaviside") simp_heaviside(z) else z)
  f <- .map_bottom_up(f, function(z) if (z$type == "Function" && z$name == "exp") expand_complex(z) else z)
  list(.tr_simplify(subs(f, t, t_), simplify), cond)
}

.ilt_build_rules <- function() {
  if (!is.null(.lt_state$irules)) return(.lt_state$irules)
  s <- Dummy("s"); t <- Dummy("t")
  a <- Wild("a", exclude = list(s)); b <- Wild("b", exclude = list(s)); c <- Wild("c", exclude = list(s))
  wilds <- list(a = a, b = b, c = c)
  rules <- list(
    list(a / s, a),
    list(b * (s + a)^(-c), t^(c - 1) * fcall("exp", list(-a * t)) / fcall("gamma", list(c))),
    list(1 / (s^2 + a^2)^2, (sin(a * t) - a * t * cos(a * t)) / (2 * a^3)),
    list(1 / s^b, t^(b - 1) / fcall("gamma", list(b))),
    list(1 / (s * (s + a)^b), fcall("lowergamma", list(b, a * t)) / (a^b * fcall("gamma", list(b)))))
  .lt_state$irules <- list(rules = rules, s = s, t = t, wilds = wilds)
  .lt_state$irules
}

.ilt_apply_simple_rules <- function(f, s, t) {
  if (is_one_num(f)) return(list(.lt_D(t), S_true))
  rr <- .ilt_build_rules()
  F <- subs(f, s, rr$s)
  for (rule in rr$rules) {
    ma <- .lt_m(F, rule[[1L]])
    if (!is.null(ma)) return(list(.mul(list(.lt_H(t), subs(.lt_xw(rule[[2L]], ma, rr$wilds), rr$t, t))), S_true))
  }
  NULL
}

.ilt_diff <- function(f, s, t, plane) {
  a <- Wild("a", exclude = list(s)); n <- Wild("n", exclude = list(s)); g <- Wild("g")
  if (!.has_type(f, "Derivative")) return(NULL)
  ma <- .lt_m(f, a * new_sym("Derivative", list(g, Tuple(s, n))))
  if (!is.null(ma) && isTRUE(ask_fact(ma$n, "integer"))) {
    r <- .ilt_transform(ma$g, s, t, plane, simplify = FALSE, dorational = FALSE)
    return(list(.mul(list(.pow(.neg(t), ma$n), r[[1L]])), r[[2L]]))
  }
  NULL
}

.ilt_time_shift <- function(F, s, t, plane) {
  a <- Wild("a", exclude = list(s)); g <- Wild("g")
  if (!.has1(F, s)) return(list(.mul(list(F, .lt_D(t))), S_true))
  if (!.has_function(F, "exp")) return(NULL)
  ma1 <- .lt_m(F, fnode("exp", list(a * s)))
  if (!is.null(ma1)) {
    if (.lt_neg(ma1$a)) return(list(.lt_D(.add(list(t, ma1$a))), S_true))
    return(list(InverseLaplaceTransform(F, s, t, plane), S_true))
  }
  ma1 <- .lt_m(F, fnode("exp", list(a * s)) * g)
  if (!is.null(ma1)) {
    if (.lt_neg(ma1$a)) return(.ilt_transform(ma1$g, s, .add(list(t, ma1$a)), plane, simplify = FALSE, dorational = TRUE))
    return(list(InverseLaplaceTransform(F, s, t, plane), S_true))
  }
  NULL
}

.ilt_freq_shift <- function(F, s, t, plane) {
  if (!.has1(F, s)) return(list(.mul(list(F, .lt_D(t))), S_true))
  if (length(F$args) == 1L) {
    a <- Wild("a", exclude = list(s))
    ma <- .lt_m(F$args[[1L]], s - a)
    if (!is.null(ma) && isTRUE(ask_fact(.lt_re(ma$a), "positive")))
      return(list(.mul(list(.lt_exp(.mul(list(ma$a, t))), InverseLaplaceTransform(sym_rebuild(F, list(s)), s, t, plane))), S_true))
  }
  NULL
}

.ilt_time_diff <- function(F, s, t, plane) {
  n <- Wild("n", exclude = list(s)); g <- Wild("g")
  ma1 <- .lt_m(F, s^n * g)
  if (!is.null(ma1) && isTRUE(ask_fact(ma1$n, "integer")) && .lt_pos(ma1$n)) {
    r <- .ilt_transform(ma1$g, s, t, plane, simplify = FALSE, dorational = TRUE)
    rr <- subs(r[[1L]], .lt_H(t), S_One)
    if (.has_type(rr, "InverseLaplaceTransform")) return(list(diff(rr, t, .int_val(ma1$n)), r[[2L]]))
    return(list(.mul(list(.lt_H(t), diff(rr, t, .int_val(ma1$n)))), r[[2L]]))
  }
  NULL
}

.ilt_irrational <- function(fn, s, t, plane) {
  a <- Wild("a", exclude = list(s)); b <- Wild("b", exclude = list(s)); m <- Wild("m", exclude = list(s)); n <- Wild("n", exclude = list(s))
  fa <- .mul_ordered_factors(fn)
  ma <- lapply(fa, function(x) .lt_m(x, (a * s^m + b)^n))
  if (any(vapply(ma, is.null, NA))) return(NULL)
  constants <- S_One; zeros <- list(); poles <- list(); rest <- list()
  for (i in seq_along(ma)) {
    term <- ma[[i]]
    if (is_zero_num(term$a)) constants <- .mul(list(constants, fa[[i]]))
    else if (.lt_pos(term$n)) zeros[[length(zeros) + 1L]] <- term
    else if (.lt_neg(term$n)) poles[[length(poles) + 1L]] <- term
    else rest[[length(rest) + 1L]] <- term
  }
  srt <- function(l) {
    if (length(l) < 2L) return(l)
    k1 <- vapply(l, function(x) num_value_dbl(x$n), 0)
    k2 <- vapply(l, function(x) !is_zero_num(x$b), NA)
    k3 <- vapply(l, function(x) { v <- tryCatch(num_value_dbl(x$b), error = function(e) NA_real_); if (is.na(v)) 0 else v }, 0)
    l[order(k1, k2, k3)]
  }
  poles <- srt(poles); zeros <- srt(zeros)
  if (length(rest)) return(NULL)
  E <- function(z) .lt_exp(z); sq <- .sqrt; erfc_ <- function(z) fcall("erfc", list(z)); erf_ <- function(z) fcall("erf", list(z))
  eqv <- function(x, v) .mi_eq(x, as_sym(v))
  P <- function(i, k) poles[[i]][[k]]; Z <- function(i, k) zeros[[i]][[k]]
  pi_ <- S_Pi; h <- S_Half
  result <- NULL
  np <- length(poles); nz <- length(zeros)
  if (np == 1L && nz == 0L) {
    if (eqv(P(1, "n"), -1) && eqv(P(1, "m"), h)) {
      a_ <- P(1, "b") / P(1, "a"); k_ <- 1 / P(1, "a") * constants
      if (.lt_pos(a_)) result <- k_ / sq(pi_) / sq(t) - k_ * a_ * E(a_^2 * t) * erfc_(a_ * sq(t))
    } else if (eqv(P(1, "n"), -2) && eqv(P(1, "m"), h)) {
      a_sq <- P(1, "b") / P(1, "a"); a_ <- a_sq^2; k_ <- 1 / P(1, "a")^2 * constants
      if (.lt_pos(a_sq)) result <- k_ * (1 - 2 / sq(pi_) * sq(a_) * sq(t) + (1 - 2 * a_ * t) * E(a_ * t) * (erf_(sq(a_) * sq(t)) - 1))
    } else if (eqv(P(1, "n"), -3) && eqv(P(1, "m"), h)) {
      a_ <- P(1, "b") / P(1, "a"); k_ <- 1 / P(1, "a")^3 * constants
      if (.lt_pos(a_)) result <- k_ * (2 / sq(pi_) * (a_^2 * t + 1) * sq(t) - a_ * t * E(a_^2 * t) * (2 * a_^2 * t + 3) * erfc_(a_ * sq(t)))
    } else if (eqv(P(1, "n"), -4) && eqv(P(1, "m"), h)) {
      a_ <- P(1, "b") / P(1, "a"); k_ <- 1 / P(1, "a")^4 * constants / 3
      if (.lt_pos(a_)) result <- k_ * (t * (4 * a_^4 * t^2 + 12 * a_^2 * t + 3) * E(a_^2 * t) * erfc_(a_ * sq(t)) - 2 / sq(pi_) * a_^3 * t^Rational(5, 2) * (2 * a_^2 * t + 5))
    } else if (eqv(P(1, "n"), -h) && eqv(P(1, "m"), 2)) {
      a_ <- sq(P(1, "b") / P(1, "a")); k_ <- 1 / sq(P(1, "a")) * constants
      result <- k_ * fcall("besselj", list(S_Zero, a_ * t))
    }
  } else if (np == 1L && nz == 1L) {
    if (eqv(P(1, "n"), -3) && eqv(P(1, "m"), h) && eqv(Z(1, "n"), h) && is_zero_num(Z(1, "b"))) {
      a_ <- P(1, "b"); k_ <- sq(Z(1, "a")) / P(1, "a") * constants
      result <- k_ * (2 * a_^4 * t^2 + 5 * a_^2 * t + 1) * E(a_^2 * t) * erfc_(a_ * sq(t)) - 2 / sq(pi_) * a_ * (a_^2 * t + 2) * sq(t)
    }
    if (eqv(P(1, "n"), -1) && eqv(P(1, "m"), 1) && eqv(Z(1, "n"), h) && eqv(Z(1, "m"), 1)) {
      a_ <- Z(1, "b") / Z(1, "a"); b_ <- P(1, "b") / P(1, "a"); k_ <- sq(Z(1, "a")) / P(1, "a") * constants
      result <- k_ * (E(-a_ * t) / sq(t) / sq(pi_) + sq(a_ - b_) * E(-b_ * t) * erf_(sq(a_ - b_) * sq(t)))
    }
  } else if (np == 2L && nz == 0L) {
    if (eqv(P(1, "n"), -1) && eqv(P(1, "m"), 1) && eqv(P(2, "n"), -h) && eqv(P(2, "m"), 1) && is_zero_num(P(2, "b"))) {
      a_ <- -P(1, "b") / P(1, "a"); k_ <- 1 / sq(P(2, "a")) / P(1, "a") * constants
      if (.lt_pos(a_)) result <- k_ / sq(a_) * E(a_ * t) * erf_(sq(a_) * sq(t))
    } else if (eqv(P(1, "n"), -1) && eqv(P(1, "m"), 1) && is_zero_num(P(1, "b")) && eqv(P(2, "n"), -1) && eqv(P(2, "m"), h)) {
      a_ <- P(2, "b") / P(2, "a"); k_ <- 1 / P(1, "a") / P(2, "a") / a_ * constants
      if (.lt_pos(a_)) result <- k_ * (1 - E(a_^2 * t) * erfc_(a_ * sq(t)))
    } else if (eqv(P(1, "n"), -1) && eqv(P(1, "m"), h) && eqv(P(2, "n"), -h) && eqv(P(2, "m"), 1) && is_zero_num(P(2, "b"))) {
      a_ <- P(1, "b") / P(1, "a"); k_ <- 1 / (P(1, "a") * sq(P(2, "a"))) * constants
      if (.lt_pos(a_)) result <- k_ * E(a_^2 * t) * erfc_(a_ * sq(t))
    } else if (eqv(P(1, "n"), Rational(-3, 2)) && eqv(P(1, "m"), 1) && is_zero_num(P(1, "b")) && eqv(P(2, "n"), -1) && eqv(P(2, "m"), h)) {
      a_ <- P(2, "b") / P(2, "a"); k_ <- 1 / (P(1, "a")^Rational(3, 2) * P(2, "a")) / a_^2 * constants
      if (.lt_pos(a_)) result <- k_ * (2 / sq(pi_) * a_ * sq(t) + E(a_^2 * t) * erfc_(a_ * sq(t)) - 1)
    } else if (eqv(P(1, "n"), -2) && eqv(P(1, "m"), h) && eqv(P(2, "n"), -1) && eqv(P(2, "m"), 1) && is_zero_num(P(2, "b"))) {
      a_sq <- P(1, "b") / P(1, "a"); a_ <- a_sq^2; k_ <- 1 / P(1, "a")^2 / P(2, "a") * constants
      if (.lt_pos(a_sq)) result <- k_ * (1 / a_ + (2 * t - 1 / a_) * E(a_ * t) * erfc_(sq(a_) * sq(t)) - 2 / sq(pi_) / sq(a_) * sq(t))
    } else if (eqv(P(1, "n"), -2) && eqv(P(1, "m"), h) && eqv(P(2, "n"), -h) && eqv(P(2, "m"), 1) && is_zero_num(P(2, "b"))) {
      a_ <- P(1, "b") / P(1, "a"); k_ <- 1 / P(1, "a")^2 / sq(P(2, "a")) * constants
      if (.lt_pos(a_)) result <- k_ * (2 / sq(pi_) * sq(t) - 2 * a_ * t * E(a_^2 * t) * erfc_(a_ * sq(t)))
    } else if (eqv(P(1, "n"), -3) && eqv(P(1, "m"), h) && eqv(P(2, "n"), -h) && eqv(P(2, "m"), 1) && is_zero_num(P(2, "b"))) {
      a_ <- P(1, "b"); k_ <- constants / sq(P(2, "a")) / P(1, "a")
      result <- k_ * ((2 * a_^2 * t + 1) * t * E(a_^2 * t) * erfc_(a_ * sq(t)) - 2 / sq(pi_) * a_ * t^Rational(3, 2))
    } else if (eqv(P(1, "n"), -1) && eqv(P(1, "m"), 1) && eqv(P(2, "n"), -h) && eqv(P(2, "m"), 1)) {
      a_ <- P(1, "b") / P(1, "a"); b_ <- P(2, "b") / P(2, "a"); k_ <- constants / sq(P(2, "a")) / P(1, "a")
      result <- k_ * (1 / sq(b_ - a_) * E(-a_ * t) * erf_(sq(b_ - a_) * sq(t)))
    }
  } else if (np == 2L && nz == 1L) {
    if (eqv(P(1, "n"), -1) && eqv(P(1, "m"), 1) && eqv(P(2, "n"), -1) && eqv(P(2, "m"), h) && eqv(Z(1, "n"), h) && eqv(Z(1, "m"), 1) && is_zero_num(Z(1, "b"))) {
      a_sq <- P(2, "b") / P(2, "a"); a_ <- a_sq^2; b_ <- -P(1, "b") / P(1, "a")
      k_ <- sq(Z(1, "a")) / P(1, "a") / P(2, "a") / (a_ - b_) * constants
      if (.lt_pos(a_sq) && .lt_pos(b_)) result <- k_ * (a_ * E(a_ * t) * erfc_(sq(a_) * sq(t)) + sq(a_) * sq(b_) * E(b_ * t) * erfc_(sq(b_) * sq(t)) - b_ * E(b_ * t))
    } else if (eqv(P(1, "n"), -1) && eqv(P(1, "m"), 1) && is_zero_num(P(1, "b")) && eqv(P(2, "n"), -1) && eqv(P(2, "m"), h) && eqv(Z(1, "n"), 1) && eqv(Z(1, "m"), h)) {
      a_num <- Z(1, "b") / Z(1, "a"); a_ <- P(2, "b") / P(2, "a")
      if (is_zero_num(a_ + a_num)) { k_ <- Z(1, "a") / P(1, "a") / P(2, "a") * constants; result <- k_ * (2 * E(a_^2 * t) * erfc_(a_ * sq(t)) - 1) }
    } else if (eqv(P(2, "n"), -1) && eqv(P(2, "m"), 1) && is_zero_num(P(2, "b")) && eqv(P(1, "n"), -2) && eqv(P(1, "m"), h) && eqv(Z(1, "n"), 2) && eqv(Z(1, "m"), h)) {
      a_num <- Z(1, "b") / Z(1, "a"); a_ <- P(1, "b") / P(1, "a")
      if (is_zero_num(a_ + a_num)) { k_ <- Z(1, "a")^2 / P(2, "a") / P(1, "a")^2 * constants; result <- k_ * (1 + 8 * a_^2 * t * E(a_^2 * t) * erfc_(a_ * sq(t)) - 8 / sq(pi_) * a_ * sq(t)) }
    } else if (eqv(P(2, "n"), -1) && eqv(P(2, "m"), 1) && is_zero_num(P(2, "b")) && eqv(P(1, "n"), -3) && eqv(P(1, "m"), h) && eqv(Z(1, "n"), 3) && eqv(Z(1, "m"), h)) {
      a_num <- Z(1, "b") / Z(1, "a"); a_ <- P(1, "b") / P(1, "a")
      if (is_zero_num(a_ + a_num)) { k_ <- Z(1, "a")^3 / P(2, "a") / P(1, "a")^3 * constants; result <- k_ * (2 * (8 * a_^4 * t^2 + 8 * a_^2 * t + 1) * E(a_^2 * t) * erfc_(a_ * sq(t)) - 8 / sq(pi_) * a_ * sq(t) * (2 * a_^2 * t + 1) - 1) }
    }
  } else if (np == 3L && nz == 0L) {
    if (eqv(P(1, "n"), -1) && is_zero_num(P(1, "b")) && eqv(P(1, "m"), 1) && eqv(P(2, "n"), -1) && eqv(P(2, "m"), 1) && eqv(P(3, "n"), -h) && eqv(P(3, "m"), 1)) {
      a_ <- -P(2, "b") / P(2, "a"); k_ <- 1 / P(1, "a") / P(2, "a") / sq(P(3, "a")) * constants
      if (.lt_pos(a_)) result <- k_ * (a_^Rational(-3, 2) * E(a_ * t) * erf_(sq(a_) * sq(t)) - 2 / a_ / sq(pi_) * sq(t))
    } else if (eqv(P(1, "n"), -1) && eqv(P(1, "m"), 1) && eqv(P(2, "n"), -1) && eqv(P(2, "m"), h) && eqv(P(3, "n"), -h) && eqv(P(3, "m"), 1) && is_zero_num(P(3, "b"))) {
      a_sq <- P(2, "b") / P(2, "a"); a_ <- a_sq^2; b_ <- -P(1, "b") / P(1, "a")
      k_ <- 1 / P(1, "a") / P(2, "a") / sq(P(3, "a")) / (sq(b_) * (a_ - b_))
      if (.lt_pos(a_sq) && .lt_pos(b_)) result <- k_ * (sq(b_) * E(a_ * t) * erfc_(sq(a_) * sq(t)) + sq(a_) * E(b_ * t) * erf_(sq(b_) * sq(t)) - sq(b_) * E(b_ * t))
    }
  }
  if (is.null(result)) return(NULL)
  list(.mul(list(.lt_H(t), result)), S_true)
}
# Mul.as_ordered_factors: the factors in printing order
.mul_ordered_factors <- function(e) if (e$type == "Mul") .sympy_args(e) else list(e)

.ilt_laplace <- function(fn, s, t, plane) {
  if (fn$type != "LaplaceTransform") return(NULL)
  if (.mi_eq(fn$args[[3L]], s)) return(list(.mul(list(subs(fn$args[[1L]], fn$args[[2L]], t), .lt_H(t))), S_true))
  NULL
}
.ilt_early_prog_rules <- function(F, s, t, plane) {
  for (rule in list(.ilt_laplace, .ilt_irrational)) { r <- rule(F, s, t, plane); if (!is.null(r)) return(r) }
  NULL
}
.ilt_apply_prog_rules <- function(F, s, t, plane) {
  for (rule in list(.ilt_time_shift, .ilt_freq_shift, .ilt_time_diff, .ilt_diff, .ilt_irrational)) { r <- rule(F, s, t, plane); if (!is.null(r)) return(r) }
  NULL
}
.ilt_expand <- function(fn, s, t, plane) {
  if (fn$type == "Add") return(NULL)
  r <- expand(fn, deep = FALSE)
  if (r$type == "Add") return(.ilt_transform(r, s, t, plane, simplify = FALSE, dorational = TRUE))
  r <- expand_mul(fn)
  if (r$type == "Add") return(.ilt_transform(r, s, t, plane, simplify = FALSE, dorational = TRUE))
  r <- expand(fn)
  if (r$type == "Add") return(.ilt_transform(r, s, t, plane, simplify = FALSE, dorational = TRUE))
  if (.mi_is_ratfun(fn, s)) r <- doit(apart(fn, s))
  if (r$type == "Add") return(.ilt_transform(r, s, t, plane, simplify = FALSE, dorational = TRUE))
  NULL
}

.ilt_rational <- function(fn, s, t, plane, simplify) {
  x_ <- Symbol("x_")
  f <- apart(fn, s)
  terms_t <- list(); conditions <- list(S_true)
  for (term in .mi_make_args_add(f)) {
    nd <- .as_numer_denom(term)
    dc <- P_all_coeffs(.P(nd[[2L]], s))
    lead <- dc[[1L]]
    dc <- lapply(dc, function(x) .div(x, lead))
    nc <- lapply(P_all_coeffs(.P(nd[[1L]], s)), function(x) .div(x, lead))
    if (length(dc) == 1L) {
      N <- length(nc) - 1L
      for (i in seq_along(nc)) terms_t[[length(terms_t) + 1L]] <- .mul(list(nc[[i]], .lt_D(t, .int(N - (i - 1L)))))
    } else if (length(dc) == 2L) {
      terms_t[[length(terms_t) + 1L]] <- .mul(list(.lt_H(t), nc[[1L]], .lt_exp(.neg(.mul(list(dc[[2L]], t))))))
    } else if (length(dc) == 3L) {
      a <- .div(dc[[2L]], S_Two)
      b <- sym_factor(.sub(dc[[3L]], .sq(a)))
      if (length(nc) == 1L) nc <- c(list(S_Zero), nc)
      l <- nc[[1L]]; m <- nc[[2L]]
      ea <- .lt_exp(.neg(.mul(list(a, t))))
      if (is_zero_num(b)) r <- .mul(list(.add(list(.mul(list(m, t)), .mul(list(l, .sub(S_One, .mul(list(a, t))))))), ea))
      else {
        hyp <- FALSE
        if (.lt_neg(b)) { b <- .neg(b); hyp <- TRUE }
        b2 <- roots(.sub(.sq(x_), b), x_)[[1L]]
        bs <- simplify(.sqrt(b))
        if (hyp) r <- .add(list(.mul(list(l, ea, cosh(.mul(list(b2, t))))), .mul(list(.div(.sub(m, .mul(list(a, l))), bs), ea, sinh(.mul(list(bs, t)))))))
        else r <- .add(list(.mul(list(l, ea, cos(.mul(list(b2, t))))), .mul(list(.div(.sub(m, .mul(list(a, l))), bs), ea, sin(.mul(list(bs, t)))))))
      }
      terms_t[[length(terms_t) + 1L]] <- .mul(list(.lt_H(t), r))
    } else {
      rr <- .ilt_transform(term, s, t, plane, simplify = simplify, dorational = FALSE)
      terms_t[[length(terms_t) + 1L]] <- rr[[1L]]
      conditions[[length(conditions) + 1L]] <- rr[[2L]]
    }
  }
  result <- .add(terms_t)
  if (simplify) result <- simplify(result, doit = FALSE)
  list(result, .lt_and(conditions))
}

.ilt_transform <- function(fn, s_, t_, plane, simplify, dorational) {
  terms_t <- list(); conditions <- list()
  for (term in .mi_make_args_add(fn)) {
    if (.has_function(term, "exp")) term <- subs(together(subs(term, s_, .neg(s_))), s_, .neg(s_))
    sp <- .as_coeff_mul_dep(term, s_); k <- sp[[1L]]; f <- .mul(sp[[2L]])
    r <- NULL
    if (dorational && .mi_is_ratfun(term, s_)) r <- .ilt_rational(f, s_, t_, plane, simplify)
    if (is.null(r)) r <- .ilt_apply_simple_rules(f, s_, t_)
    if (is.null(r)) r <- .ilt_early_prog_rules(f, s_, t_, plane)
    if (is.null(r)) r <- .ilt_expand(f, s_, t_, plane)
    if (is.null(r)) r <- .ilt_apply_prog_rules(f, s_, t_, plane)
    if (is.null(r)) {
      if (.lt_applied_undef_has(f, s_)) r <- list(InverseLaplaceTransform(f, s_, t_, plane), S_true)
      else {
        r <- tryCatch(.ilt_integration(f, s_, t_, plane, simplify), error = function(e) NULL)
        if (is.null(r)) r <- list(InverseLaplaceTransform(f, s_, t_, plane), S_true)
      }
    }
    terms_t[[length(terms_t) + 1L]] <- .mul(list(k, r[[1L]]))
    conditions[[length(conditions) + 1L]] <- r[[2L]]
  }
  result <- .add(terms_t)
  if (simplify) result <- simplify(result, doit = FALSE)
  list(result, .lt_and(conditions))
}

#' @rdname laplace_transform
#' @export
inverse_laplace_transform <- function(F, s, t, plane = NULL, noconds = TRUE, simplify = FALSE) {
  if (inherits(F, "symmatrix")) return(.mat_map(F, function(z) inverse_laplace_transform(z, s, t, plane, noconds = noconds, simplify = simplify)))
  F <- as_sym(F); s <- as_sym(s); t <- as_sym(t)
  if (!is.null(plane)) plane <- as_sym(plane)
  r <- .ilt_transform(F, s, t, plane, simplify = simplify, dorational = TRUE)
  if (noconds) r[[1L]] else r
}
