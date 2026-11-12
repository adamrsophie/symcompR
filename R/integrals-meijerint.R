# Integration by Meijer G functions, after sympy/integrals/meijerint.py
# (part 1: the lookup table, helpers, the condition simplifier, the Saxena
# rewrites and the convergence antecedents).

.mi_state <- new.env()
.mi_z <- function() { if (is.null(.mi_state$z)) .mi_state$z <- Dummy("z"); .mi_state$z }

# _mytype: the sorted classes of the functions of x in f
.mi_mytype <- function(f, x) {
  if (!.has1(f, x)) return(character())
  if (f$type == "Function") return(f$name)
  ts <- unlist(lapply(f$args, .mi_mytype, x = x))
  if (!length(ts)) return(character())
  # class_key order: functions by their class index, then name
  ks <- lapply(ts, function(t) .class_key(new_sym("Function", list(), name = t)))
  unlist(.key_sort(ks, as.list(ts)))
}
.mi_typekey <- function(t) paste0("T:", paste(t, collapse = ","))

.mi_mg <- function(an, ap, bm, bq, arg) meijerg(an, ap, bm, bq, arg)

# _create_lookup_table
.mi_table <- function() {
  if (!is.null(.mi_state$table)) return(.mi_state$table)
  z <- .mi_z()
  wild <- function(n) Wild(n, exclude = list(z))
  p <- wild("p"); q <- wild("q"); a <- wild("a"); b <- wild("b"); c <- wild("c")
  n <- Wild("n", properties = list(function(v) v$type == "Integer" && v$value > 0))
  t <- p * z^q
  .mi_state$wilds <- list(p = p, q = q, a = a, b = b, c = c, n = n)
  tab <- list()
  push <- function(formula, inst, cond, hint) {
    k <- .mi_typekey(.mi_mytype(formula, z))
    tab[[k]] <<- c(tab[[k]], list(list(formula = formula, inst = inst, cond = cond, hint = hint)))
  }
  add <- function(formula, an, ap, bm, bq, arg = t, fac = S_One, cond = S_true, hint = S_true)
    push(formula, list(list(as_sym(fac), .mi_mg(an, ap, bm, bq, arg))), as_sym(cond), as_sym(hint))
  addi <- function(formula, inst, cond, hint = S_true) push(formula, inst, as_sym(cond), as_sym(hint))
  constant <- function(aa) list(list(aa, .mi_mg(list(1), list(), list(), list(0), z)), list(aa, .mi_mg(list(), list(1), list(0), list(), z)))
  tab[[.mi_typekey(character())]] <- list(list(formula = a, inst = constant(a), cond = S_true, hint = S_true))
  pi <- S_Pi; h <- S_Half
  add(Heaviside(t - b) * (t - b)^(a - 1), list(a), list(), list(), list(0), t / b, gamma(a) * b^(a - 1), Gt(b, 0))
  add(Heaviside(b - t) * (b - t)^(a - 1), list(), list(a), list(0), list(), t / b, gamma(a) * b^(a - 1), Gt(b, 0))
  add(Heaviside(z - (b / p)^(1 / q)) * (t - b)^(a - 1), list(a), list(), list(), list(0), t / b, gamma(a) * b^(a - 1), Gt(b, 0))
  add(Heaviside((b / p)^(1 / q) - z) * (b - t)^(a - 1), list(), list(a), list(0), list(), t / b, gamma(a) * b^(a - 1), Gt(b, 0))
  add((b + t)^(-a), list(1 - a), list(), list(0), list(), t / b, b^(-a) / gamma(a), hint = Not(fcall("IsNonPositiveInteger", list(a))))
  add(Abs(b - t)^(-a), list(1 - a), list((1 - a) / 2), list(0), list((1 - a) / 2), t / b,
      2 * sin(pi * a / 2) * gamma(1 - a) * Abs(b)^(-a), Lt(re(a), 1))
  add((t^a - b^a) / (t - b), list(0, a), list(), list(0, a), list(), t / b, b^(a - 1) * sin(a * pi) / pi)
  A1 <- function(r, sgn, nu) pi^Rational(-1, 2) * (-sgn * nu / 2)^(1 - 2 * r)
  for (rs in list(list(S_Zero, 1), list(S_Zero, -1), list(h, 1), list(h, -1))) {
    r <- rs[[1L]]; sg <- rs[[2L]]
    add((sqrt(a^2 + t) + sg * a)^b / (a^2 + t)^r, list((1 + b) / 2, 1 - 2 * r + b / 2), list(),
        list((b - sg * b) / 2), list((b + sg * b) / 2), t / a^2, a^(b - 2 * r) * A1(r, sg, b))
  }
  for (rs in list(list(S_Zero, 1), list(S_Zero, -1), list(h, 1), list(h, -1))) {
    r <- rs[[1L]]; sg <- rs[[2L]]
    add((sqrt(a + p * z^q) + sg * sqrt(p) * z^(q / 2))^b / (a + p * z^q)^r, list(1 - r + sg * b / 2), list(1 - r - sg * b / 2),
        list(0, h), list(), p * z^q / a, a^(b / 2 - r) * A1(r, sg, b))
  }
  add(exp(polar_lift(-1) * t), list(), list(), list(0), list())
  add(sinh(t), list(), list(1), list(h), list(1, 0), t^2 / 4, pi^Rational(3, 2))
  add(cosh(t), list(), list(h), list(0), list(h, h), t^2 / 4, pi^Rational(3, 2))
  add(sin(t), list(), list(), list(h), list(0), t^2 / 4, sqrt(pi))
  add(cos(t), list(), list(), list(0), list(h), t^2 / 4, sqrt(pi))
  add(fcall("sinc", list(t)), list(), list(), list(0), list(Rational(-1, 2)), t^2 / 4, sqrt(pi) / 2)
  make_log1 <- function(subs) { N <- .int_val(subs[["n"]])
    list(list((-1)^N * factorial(Integer(N)), .mi_mg(list(), rep(list(1), N + 1), rep(list(0), N + 1), list(), t))) }
  make_log2 <- function(subs) { N <- .int_val(subs[["n"]])
    list(list(factorial(Integer(N)), .mi_mg(rep(list(1), N + 1), list(), list(), rep(list(0), N + 1), t))) }
  addi(log(t)^n * Heaviside(1 - t), make_log1, S_true)
  addi(log(t)^n * Heaviside(t - 1), make_log2, S_true)
  addi(log(t)^n, function(subs) c(make_log1(subs), make_log2(subs)), S_true)
  addi(log(t + a), c(constant(log(a)), list(list(S_One, .mi_mg(list(1, 1), list(), list(1), list(0), t / a)))), S_true)
  addi(log(Abs(t - a)), c(constant(log(Abs(a))), list(list(pi, .mi_mg(list(1, 1), list(h), list(1), list(0, h), t / a)))), S_true)
  addi(Ei(t), c(constant(-sym_I * pi), list(list(S_NegativeOne, .mi_mg(list(), list(1), list(0, 0), list(), t * polar_lift(-1))))), S_true)
  add(Si(t), list(1), list(), list(h), list(0, 0), t^2 / 4, sqrt(pi) / 2)
  add(Ci(t), list(), list(1), list(0, 0), list(h), t^2 / 4, -sqrt(pi) / 2)
  add(Shi(t), list(h), list(), list(0), list(Rational(-1, 2), Rational(-1, 2)), polar_lift(-1) * t^2 / 4, t * sqrt(pi) / 4)
  add(Chi(t), list(), list(h, 1), list(0, 0), list(h, h), t^2 / 4, -pi^Rational(3, 2) / 2)
  add(expint(a, t), list(), list(a), list(a - 1, 0), list(), t)
  add(erf(t), list(1), list(), list(h), list(0), t^2, 1 / sqrt(pi))
  add(erfc(t), list(), list(1), list(0, h), list(), t^2, 1 / sqrt(pi))
  add(erfi(t), list(h), list(), list(0), list(Rational(-1, 2)), -t^2, t / sqrt(pi))
  add(fresnels(t), list(1), list(), list(Rational(3, 4)), list(0, Rational(1, 4)), pi^2 * t^4 / 16, h)
  add(fresnelc(t), list(1), list(), list(Rational(1, 4)), list(0, Rational(3, 4)), pi^2 * t^4 / 16, h)
  add(besselj(a, t), list(), list(), list(a / 2), list(-a / 2), t^2 / 4)
  add(bessely(a, t), list(), list(-(a + 1) / 2), list(a / 2, -a / 2), list(-(a + 1) / 2), t^2 / 4)
  add(besseli(a, t), list(), list((1 + a) / 2), list(a / 2), list(-a / 2, (1 + a) / 2), t^2 / 4, pi)
  add(besselk(a, t), list(), list(), list(a / 2, -a / 2), list(), t^2 / 4, h)
  add(elliptic_k(t), list(h, h), list(), list(0), list(0), -t, h)
  add(elliptic_e(t), list(h, 3 * h), list(), list(0), list(0), -t, Rational(-1, 2) / 2)
  .mi_state$table <- tab
  tab
}

# IsNonPositiveInteger of the table: True or False for integers
def_function("IsNonPositiveInteger", eval = function(args) {
  a <- unpolarify(args[[1L]])
  if (a$type == "Integer") return(if (a$value <= 0) S_true else S_false)
  NULL
})

.mi_err <- function(msg) stop(structure(class = c("mi_coeff_exp_error", "error", "condition"), list(message = msg, call = NULL)))

# _get_coeff_exp
.mi_coeff_exp <- function(expr, x) {
  e <- expand_power_base(powsimp(as_sym(expr)))
  cm <- .as_coeff_mul_dep(e, x)
  c <- cm[[1L]]; m <- cm[[2L]]
  if (!length(m)) return(list(c, S_Zero))
  if (length(m) != 1L) .mi_err("expr not of form a*x**b")
  m <- m[[1L]]
  if (m$type == "Pow") {
    if (!identical(m$args[[1L]]$key, x$key)) .mi_err("expr not of form a*x**b")
    return(list(c, m$args[[2L]]))
  }
  if (identical(m$key, x$key)) return(list(c, S_One))
  .mi_err(paste0("expr not of form a*x**b: ", .str(expr)))
}

.mi_exponents <- function(expr, x) {
  res <- list()
  rec <- function(e) {
    if (identical(e$key, x$key)) { res[["1"]] <<- S_One; return(invisible()) }
    if (e$type == "Pow" && identical(e$args[[1L]]$key, x$key)) { res[[e$args[[2L]]$key]] <<- e$args[[2L]]; return(invisible()) }
    for (a in e$args) rec(a)
  }
  rec(expr)
  unname(res)
}
.mi_functions <- function(expr, x) unique(vapply(Filter(function(g) .has1(g, x), atoms(expr, "Function")), `[[`, "", "name"))

.mi_find_splitting_points <- function(expr, x) {
  p <- Wild("p", exclude = list(x)); q <- Wild("q", exclude = list(x))
  res <- list(); stack <- list(expr)
  while (length(stack)) {
    e <- stack[[length(stack)]]; stack[[length(stack)]] <- NULL
    if (!.is_expr_node(e)) next
    m <- tryCatch(sym_match(e, p * x + q), error = function(err) NULL)
    if (!is.null(m) && !is_zero_num(m[["p"]])) { v <- .neg(.div(m[["q"]], m[["p"]])); res[[v$key]] <- v }
    else if (length(e$args)) stack <- c(stack, e$args)
  }
  unname(res)
}

# _split_mul: f = fac * po * g with fac free of x and po a power of x
.mi_split_mul <- function(f, x) {
  fac <- S_One; po <- S_One; g <- S_One
  f <- expand_power_base(f)
  for (a in (if (f$type == "Mul") .sympy_args(f) else list(f))) {
    if (identical(a$key, x$key)) { po <- .mul(list(po, x)); next }
    if (!.has1(a, x)) { fac <- .mul(list(fac, a)); next }
    if (a$type == "Pow" && !.has1(a$args[[2L]], x)) {
      cm <- .as_coeff_mul_dep(a$args[[1L]], x); c <- cm[[1L]]; tt <- cm[[2L]]
      if (!(length(tt) == 1L && identical(tt[[1L]]$key, x$key))) {
        cm <- .as_coeff_mul_dep(expand_mul(a$args[[1L]]), x); c <- cm[[1L]]; tt <- cm[[2L]]
      }
      if (length(tt) == 1L && identical(tt[[1L]]$key, x$key)) {
        po <- .mul(list(po, .pow(x, a$args[[2L]])))
        fac <- .mul(list(fac, unpolarify(polarify(.pow(c, a$args[[2L]]), subs = FALSE))))
        next
      }
    }
    g <- .mul(list(g, a))
  }
  list(fac, po, g)
}

.mi_mul_args <- function(f) {
  gs <- list()
  for (g in (if (f$type == "Mul") .sympy_args(f) else list(f))) {
    if (g$type == "Pow" && g$args[[2L]]$type == "Integer") {
      n <- .int_val(g$args[[2L]]); base <- g$args[[1L]]
      if (n < 0) { n <- -n; base <- .inv(base) }
      gs <- c(gs, rep(list(base), n))
    } else gs[[length(gs) + 1L]] <- g
  }
  gs
}

.mi_mul_as_two_parts <- function(f) {
  gs <- .mi_mul_args(f)
  if (length(gs) < 2L) return(NULL)
  if (length(gs) == 2L) return(list(gs))
  lapply(multiset_partitions(gs, 2), function(pp) list(.mul(pp[[1L]]), .mul(pp[[2L]])))
}

.mi_inflate_g <- function(g, n) {
  inflate <- function(params) unlist(lapply(params, function(a) lapply(seq_len(n) - 1L, function(i) .div(.add(list(a, .int(i))), .int(n)))), recursive = FALSE)
  v <- .int(length(.mg_ap(g)) - length(.mg_bq(g)))
  C <- .pow(.int(n), .add(list(S_One, .mg_nu(g), .div(v, S_Two))))
  C <- .div(C, .pow(.mul(list(S_Two, S_Pi)), .mul(list(.int(n - 1L), .mg_delta(g)))))
  list(C, .mi_mg(inflate(.mg_an(g)), inflate(.mg_aother(g)), inflate(.mg_bm(g)), inflate(.mg_bother(g)),
                 .mul(list(.pow(g$args[[3L]], .int(n)), .pow(.int(n), .mul(list(.int(n), v)))))))
}

.mi_flip_g <- function(g) {
  tr <- function(l) lapply(l, function(a) .sub(S_One, a))
  .mi_mg(tr(.mg_bm(g)), tr(.mg_bother(g)), tr(.mg_an(g)), tr(.mg_aother(g)), .inv(g$args[[3L]]))
}

.mi_inflate_fox_h <- function(g, a) {
  if (num_sign(a) < 0) return(.mi_inflate_fox_h(.mi_flip_g(g), .neg(a)))
  P <- gmp::numerator(gmp::as.bigq(a$value)); Q <- gmp::denominator(gmp::as.bigq(a$value))
  p <- .zs(P); q <- as.integer(as.character(Q))
  r <- .mi_inflate_g(g, q); D <- r[[1L]]; g <- r[[2L]]
  zz <- g$args[[3L]]
  D <- .div(D, .mul(list(.pow(.mul(list(S_Two, S_Pi)), .div(.sub(S_One, p), S_Two)), .pow(p, Rational(-1, 2)))))
  zz <- .div(zz, .pow(p, p))
  pi_ <- as.integer(as.character(P))
  bs <- lapply(seq_len(pi_) - 1L, function(nn) .div(.int(nn + 1L), p))
  list(D, .mi_mg(.mg_an(g), .mg_aother(g), .mg_bm(g), c(.mg_bother(g), bs), zz))
}

.mi_dummy_ <- function(name, token, ...) {
  k <- paste(name, token)
  if (is.null(.mi_state$dummies[[k]])) .mi_state$dummies[[k]] <- Dummy(name, ...)
  .mi_state$dummies[[k]]
}

.mi_dummy <- function(name, token, expr, ...) {
  k <- paste(name, token)
  if (is.null(.mi_state$dummies[[k]])) .mi_state$dummies[[k]] <- Dummy(name, ...)
  d <- .mi_state$dummies[[k]]
  if (.has1(expr, d)) return(Dummy(name, ...))
  d
}

.mi_is_analytic <- function(f, x) !any(vapply(atoms(f, "Function"), function(e) e$name %in% c("Heaviside", "Abs") && .has1(e, x), NA))

# ---------------------------------------------------------------------------
# _condsimp

.mi_cond_rules <- function() {
  if (!is.null(.mi_state$rules)) return(.mi_state$rules)
  p <- Wild("p"); q <- Wild("q"); r <- Wild("r"); pi <- S_Pi; I <- sym_I
  ar <- function(e) fcall("arg", list(e)); ua <- function(e) unbranched_argument(e)
  rules <- list(
    list(Or(Lt(p, q), Eq(p, q)), Le(p, q)),
    list(And(Le(Abs(ar(p)), pi), Le(Abs(ar(p) - 2 * pi), pi)), Eq(ar(p) - pi, 0)),
    list(And(Le(Abs(2 * ar(p) + pi), pi), Le(Abs(2 * ar(p) - pi), pi)), Eq(ar(p), 0)),
    list(And(Lt(Abs(2 * ar(p) + pi), pi), Le(Abs(2 * ar(p) - pi), pi)), S_false),
    list(And(Le(Abs(ar(p) - pi / 2), pi / 2), Le(Abs(ar(p) + pi / 2), pi / 2)), Eq(ar(p), 0)),
    list(And(Le(Abs(ar(p) - pi / 2), pi / 2), Lt(Abs(ar(p) + pi / 2), pi / 2)), S_false),
    list(And(Lt(Abs(ar(p^2 / 2 + 1)), pi), Ne(Abs(ar(p^2 / 2 + 1)), pi)), S_true),
    list(Or(Lt(Abs(ar(p^2 / 2 + 1)), pi), Ne(1 / (p^2 / 2 + 1), 0)), S_true),
    list(And(Le(Abs(ua(p)), pi), Le(Abs(ua(exp_polar(-2 * pi * I) * p)), pi)), Eq(ua(exp_polar(-I * pi) * p), 0)),
    list(And(Le(Abs(ua(p)), pi / 2), Le(Abs(ua(exp_polar(-pi * I) * p)), pi / 2)), Eq(ua(exp_polar(-I * pi / 2) * p), 0)),
    list(Or(Le(p, q), And(Lt(p, q), r)), Le(p, q)),
    list(And(Ne(p^2, 1), Gt(p^2, 1)), Gt(p^2, 1)),
    list(And(Ne(1 / p, 1), Gt(cos(Abs(ar(p))) * Abs(p), 1)), Gt(Abs(p), 1)),
    list(And(Ne(p, 2), Gt(cos(Abs(ar(p))) * Abs(p), 2)), Gt(Abs(p), 2)),
    list(And(Lt(Abs(ar(p)), pi / 2), Gt(cos(Abs(ar(p))) * sqrt(Abs(p^2)), 1)), Gt(p^2, 1)))
  .mi_state$rules <- rules
  rules
}

.mi_is_bool_fn <- function(e) is_sym(e) && e$type %in% c("And", "Or", "Not", "Xor", "Implies", "Equivalent", "ITE")
.mi_bool_rebuild <- function(type, args) switch(type, And = do.call(And, args), Or = do.call(Or, args), Not = Not(args[[1L]]),
                                                 sym_rebuild(new_sym(type, args), args))

.mi_condsimp <- function(cond, first = TRUE) {
  cond <- as_sym(cond)
  if (first) { cond <- .map_bottom_up(cond, function(z) if (is_relational(z)) .canonical_coeff(z) else z); first <- FALSE }
  if (!.mi_is_bool_fn(cond)) return(cond)
  rules <- .mi_cond_rules()
  rkey <- Wild("r")$key
  cond <- .mi_bool_rebuild(cond$type, lapply(cond$args, .mi_condsimp, first = first))
  change <- TRUE
  while (change) {
    change <- FALSE
    if (!.mi_is_bool_fn(cond)) break
    for (rule in rules) {
      fro <- rule[[1L]]; to <- rule[[2L]]
      if (!identical(fro$type, cond$type)) next
      applied <- FALSE
      for (nn in seq_along(cond$args)) {
        arg1 <- cond$args[[nn]]
        if (.contains_key(fro$args[[1L]], rkey)) { m <- tryCatch(sym_match(arg1, fro$args[[2L]]), error = function(e) NULL); num <- 2L }
        else { num <- 1L; m <- tryCatch(sym_match(arg1, fro$args[[1L]]), error = function(e) NULL) }
        if (is.null(m) || !length(m)) next
        reps <- lapply(names(m), function(k) list(Wild(k), m[[k]]))
        others <- lapply(fro$args[-num], function(z) subs(z, reps))
        olist <- nn
        for (a2 in others) for (k in seq_along(cond$args)) {
          if (k %in% olist) next
          a3 <- cond$args[[k]]
          if (identical(a2$key, a3$key)) { olist <- c(olist, k); break }
          if (a3$type == "And" && a2$type == "And" && identical(a2$args[[2L]]$key, rkey) && any(vapply(a3$args, function(z) identical(z$key, a2$args[[1L]]$key), NA))) { olist <- c(olist, k); break }
          if (a3$type == "And" && a2$type == "And" && identical(a2$args[[1L]]$key, rkey) && any(vapply(a3$args, function(z) identical(z$key, a2$args[[2L]]$key), NA))) { olist <- c(olist, k); break }
        }
        if (length(olist) != length(others) + 1L) next
        newargs <- c(cond$args[-olist], list(subs(to, reps)))
        cond <- .mi_bool_rebuild(cond$type, newargs)
        change <- TRUE; applied <- TRUE
        break
      }
      if (applied) break
    }
  }
  touch <- function(rel) {
    if (rel$type != "Equality" || !is_zero_num(rel$args[[2L]])) return(rel)
    L <- rel$args[[1L]]; p <- Wild("p"); q <- Wild("q")
    m <- tryCatch(sym_match(L, fcall("arg", list(p))^q), error = function(e) NULL)
    if (is.null(m) || !length(m)) m <- tryCatch(sym_match(L, unbranched_argument(polar_lift(p)^q)), error = function(e) NULL)
    if (is.null(m) || !length(m)) {
      if (L$type == "Function" && L$name == "periodic_argument" && !isTRUE(.is_polar(L$args[[1L]])) && is_const(L$args[[2L]], "oo"))
        return(Gt(L$args[[1L]], 0))
      return(rel)
    }
    Gt(m[["p"]], 0)
  }
  .map_bottom_up(cond, function(z) if (is_relational(z)) touch(z) else z)
}

.mi_eval_cond <- function(cond) {
  if (is.logical(cond)) return(cond)
  .mi_condsimp(doit(cond))
}
.mi_bool <- function(c) { if (is.logical(c)) return(c); if (identical(c$key, S_true$key)) TRUE else if (identical(c$key, S_false$key)) FALSE else NA }

.mi_principal_branch <- function(expr, period, full_pb = FALSE) {
  res <- principal_branch(expr, period)
  if (!full_pb) res <- .map_bottom_up(res, function(z) if (z$type == "Function" && z$name == "principal_branch") z$args[[1L]] else z)
  res
}

.mi_rewrite_saxena_1 <- function(fac, po, g, x) {
  s <- .mi_coeff_exp(po, x)[[2L]]
  ab <- .mi_coeff_exp(g$args[[3L]], x); a <- ab[[1L]]; b <- ab[[2L]]
  a <- .mi_principal_branch(a, .mg_get_period(g))
  C <- .div(fac, .mul(list(Abs(b), .pow(a, .sub(.div(.add(list(s, S_One)), b), S_One)))))
  tr <- function(l) lapply(l, function(aa) .sub(.add(list(aa, .div(.add(list(S_One, s)), b))), S_One))
  list(C, .mi_mg(tr(.mg_an(g)), tr(.mg_aother(g)), tr(.mg_bm(g)), tr(.mg_bother(g)), .mul(list(a, x))))
}

.mi_AND <- function(l) do.call(And, lapply(l, as_sym))
.mi_OR <- function(l) do.call(Or, lapply(l, as_sym))
.mi_tb <- function(v) if (isTRUE(v)) S_true else S_false

.mi_check_antecedents_1 <- function(g, x, helper = FALSE) {
  delta <- .mg_delta(g)
  eta <- .mi_coeff_exp(g$args[[3L]], x)[[1L]]
  m <- length(.mg_bm(g)); n <- length(.mg_an(g)); p <- length(.mg_ap(g)); q <- length(.mg_bq(g))
  if (p > q) {
    tr <- function(l) lapply(l, function(a) .sub(S_One, a))
    return(.mi_check_antecedents_1(.mi_mg(tr(.mg_bm(g)), tr(.mg_bother(g)), tr(.mg_an(g)), tr(.mg_aother(g)), .div(x, eta)), x))
  }
  I <- function(v) .int(v)
  tmp <- c(lapply(.mg_bm(g), function(b) Lt(-re(b), 1)), lapply(.mg_an(g), function(a) Lt(1, 1 - re(a))))
  cond_3 <- .mi_AND(tmp)
  tmp <- c(tmp, lapply(.mg_bother(g), function(b) Lt(-re(b), 1)), lapply(.mg_aother(g), function(a) Lt(1, 1 - re(a))))
  cond_3_star <- .mi_AND(tmp)
  cond_4 <- Gt(-re(.mg_nu(g)) + Rational(q + 1 - p, 2), q - p)
  ua <- function(e) unbranched_argument(e)
  conds <- list()
  tmp1 <- list(Le(1, n), Lt(p, q), Le(1, m))
  tmp2 <- list(Le(1, p), Le(1, m), Eq(I(q), p + 1), Not(And(Eq(I(n), 0), Eq(I(m), p + 1))))
  tmp3 <- list(Le(1, p), Eq(I(q), p))
  for (k in 0:.int_val(fcall("ceiling", list(.div(delta, S_Two))))) tmp3[[length(tmp3) + 1L]] <- Ne(Abs(ua(eta)), (delta - 2 * k) * S_Pi)
  tt <- list(Gt(delta, 0), Lt(Abs(ua(eta)), delta * S_Pi))
  extra <- if (helper) list() else list(Ne(eta, 0), cond_3)
  for (t in list(tmp1, tmp2, tmp3)) conds[[length(conds) + 1L]] <- .mi_AND(c(t, tt, extra))
  extra <- if (helper) list() else list(cond_3)
  conds[[length(conds) + 1L]] <- .mi_AND(c(list(Eq(I(n), 0), Le(p + 1, m), Le(m, q), Lt(Abs(ua(eta)), delta * S_Pi)), extra))
  extra <- if (helper) list() else list(cond_3, cond_4)
  conds[[length(conds) + 1L]] <- .mi_AND(c(list(Lt(p, q), Le(1, m), Gt(delta, 0), Eq(Abs(ua(eta)), delta * S_Pi)), extra))
  conds[[length(conds) + 1L]] <- .mi_AND(c(list(Le(p, q - 2), Eq(delta, 0), Eq(Abs(ua(eta)), 0)), extra))
  ce <- list(Eq(I(p), q), Eq(delta, 0), Eq(ua(eta), 0), Ne(eta, 0))
  if (!helper) ce <- c(ce, list(cond_3))
  s <- Map(function(a, b) .sub(b, a), .mg_ap(g), .mg_bq(g))
  ce <- c(ce, list(Lt(re(.add(s)), 0)))
  conds[[length(conds) + 1L]] <- .mi_AND(ce)
  ce2 <- list(Gt(delta, 0), Lt(Abs(ua(eta)), delta * S_Pi))
  if (!helper) ce2 <- c(ce2, list(cond_3))
  conds[[length(conds) + 1L]] <- .mi_AND(ce2)
  .mi_OR(conds)
}

.mi_int0oo_1 <- function(g, x) {
  eta <- .mi_coeff_exp(g$args[[3L]], x)[[1L]]
  res <- .inv(eta)
  for (b in .mg_bm(g)) res <- .mul(list(res, gamma(.add(list(b, S_One)))))
  for (a in .mg_an(g)) res <- .mul(list(res, gamma(.sub(.sub(S_One, a), S_One))))
  for (b in .mg_bother(g)) res <- .div(res, gamma(.sub(.sub(S_One, b), S_One)))
  for (a in .mg_aother(g)) res <- .div(res, gamma(.add(list(a, S_One))))
  gammasimp(unpolarify(res))
}

.mi_rewrite_saxena <- function(fac, po, g1, g2, x, full_pb = FALSE) {
  pb <- function(g) {
    ab <- .mi_coeff_exp(g$args[[3L]], x)
    .mi_mg(.mg_an(g), .mg_aother(g), .mg_bm(g), .mg_bother(g), .mul(list(.mi_principal_branch(ab[[1L]], .mg_get_period(g), full_pb), .pow(x, ab[[2L]]))))
  }
  s <- .mi_coeff_exp(po, x)[[2L]]
  b1 <- .mi_coeff_exp(g1$args[[3L]], x)[[2L]]; b2 <- .mi_coeff_exp(g2$args[[3L]], x)[[2L]]
  if (.is_true(Lt(b1, 0))) { b1 <- .neg(b1); g1 <- .mi_flip_g(g1) }
  if (.is_true(Lt(b2, 0))) { b2 <- .neg(b2); g2 <- .mi_flip_g(g2) }
  if (!(b1$type %in% c("Integer", "Rational")) || !(b2$type %in% c("Integer", "Rational"))) return(NULL)
  q1 <- gmp::as.bigq(b1$value); q2 <- gmp::as.bigq(b2$value)
  m1 <- gmp::numerator(q1); n1 <- gmp::denominator(q1); m2 <- gmp::numerator(q2); n2 <- gmp::denominator(q2)
  tau <- gmp::lcm.bigz(m1 * n2, m2 * n1)
  r1 <- as.integer(as.character(tau %/% (m1 * n2))); r2 <- as.integer(as.character(tau %/% (m2 * n1)))
  i1 <- .mi_inflate_g(g1, r1); i2 <- .mi_inflate_g(g2, r2)
  g1 <- pb(i1[[2L]]); g2 <- pb(i2[[2L]])
  fac <- .mul(list(fac, i1[[1L]], i2[[1L]]))
  ab1 <- .mi_coeff_exp(g1$args[[3L]], x); a1 <- ab1[[1L]]; b <- ab1[[2L]]
  a2 <- .mi_coeff_exp(g2$args[[3L]], x)[[1L]]
  ex <- .sub(.div(.add(list(s, S_One)), b), S_One)
  fac <- .div(fac, .mul(list(Abs(b), .pow(a1, ex))))
  tr <- function(l) lapply(l, function(a) .add(list(a, ex)))
  g1 <- .mi_mg(tr(.mg_an(g1)), tr(.mg_aother(g1)), tr(.mg_bm(g1)), tr(.mg_bother(g1)), .mul(list(a1, x)))
  g2 <- .mi_mg(.mg_an(g2), .mg_aother(g2), .mg_bm(g2), .mg_bother(g2), .mul(list(a2, x)))
  list(powdenest(fac, polar = TRUE), g1, g2)
}

# _check_antecedents: the conditions of Prudnikov et al. for the integral of
# a product of two G functions over (0, oo)
.mi_check_antecedents <- function(g1, g2, x) {
  sigma <- .mi_coeff_exp(g1$args[[3L]], x)[[1L]]; omega <- .mi_coeff_exp(g2$args[[3L]], x)[[1L]]
  s <- length(.mg_bm(g1)); t <- length(.mg_an(g1)); u <- length(.mg_ap(g1)); v <- length(.mg_bq(g1))
  m <- length(.mg_bm(g2)); n <- length(.mg_an(g2)); p <- length(.mg_ap(g2)); q <- length(.mg_bq(g2))
  I <- function(k) .int(k); PI <- S_Pi
  bstar <- I(s + t) - Rational(u + v, 2); cstar <- I(m + n) - Rational(p + q, 2)
  rho <- .mg_nu(g1) + Rational(u - v, 2) + 1; mu <- .mg_nu(g2) + Rational(p - q, 2) + 1
  phi <- I(q - p - (v - u)); eta <- 1 - (v - u) - mu - rho
  ua <- function(e) unbranched_argument(e)
  psi <- (PI * (q - m - n) + Abs(ua(omega))) / (q - p)
  theta <- (PI * (v - s - t) + Abs(ua(sigma))) / (v - u)
  pos <- function(e) .mi_tb(isTRUE(ask_fact(e, "positive")))
  neg <- function(e) .mi_tb(isTRUE(ask_fact(e, "negative")))
  c1 <- {
    ok <- TRUE
    for (g in list(g1, g2)) for (i in .mg_an(g)) for (j in .mg_bm(g)) {
      d <- .sub(i, j)
      if (isTRUE(ask_fact(d, "integer")) && isTRUE(ask_fact(d, "positive"))) ok <- FALSE
    }
    .mi_tb(ok)
  }
  c2 <- .mi_AND(unlist(lapply(.mg_bm(g1), function(i) lapply(.mg_bm(g2), function(j) Gt(re(1 + i + j), 0))), recursive = FALSE))
  c3 <- .mi_AND(unlist(lapply(.mg_an(g1), function(i) lapply(.mg_an(g2), function(j) Lt(re(1 + i + j), 1 + 1))), recursive = FALSE))
  c4 <- .mi_AND(lapply(.mg_an(g1), function(i) Gt((p - q) * re(1 + i - 1) - re(mu), Rational(-3, 2))))
  c5 <- .mi_AND(lapply(.mg_bm(g1), function(i) Gt((p - q) * re(1 + i) - re(mu), Rational(-3, 2))))
  c6 <- .mi_AND(lapply(.mg_an(g2), function(i) Gt((u - v) * re(1 + i - 1) - re(rho), Rational(-3, 2))))
  c7 <- .mi_AND(lapply(.mg_bm(g2), function(i) Gt((u - v) * re(1 + i) - re(rho), Rational(-3, 2))))
  w <- re((rho - 1) * (q - p) + (v - u) * (q - p) + (mu - 1) * (v - u))
  c8 <- Gt(Abs(phi) + 2 * w, 0); c9 <- Gt(Abs(phi) - 2 * w, 0)
  c10 <- Lt(Abs(ua(sigma)), bstar * PI); c11 <- Eq(Abs(ua(sigma)), bstar * PI)
  c12 <- Lt(Abs(ua(omega)), cstar * PI); c13 <- Eq(Abs(ua(omega)), cstar * PI)
  z0 <- exp(-(bstar + cstar) * PI * sym_I)
  zos <- unpolarify(z0 * omega / sigma); zso <- unpolarify(z0 * sigma / omega)
  if (identical(zos$key, .inv(zso)$key)) {
    c14 <- And(Eq(phi, 0), Le(bstar + cstar, 1), Or(Ne(zos, 1), Lt(re(mu + rho + v - u), 1), Lt(re(mu + rho + q - p), 1)))
  } else {
    cnd <- function(zz) if (is_one_num(zz)) S_false else Lt(Abs(fcall("arg", list(1 - zz))), PI)
    c14 <- And(Eq(phi, 0), Le(bstar - 1 + cstar, 0), Or(And(Ne(zos, 1), cnd(zos)), And(Lt(re(mu + rho + v - u), 1), Eq(zos, 1))))
    c14a <- And(Eq(phi, 0), Le(cstar - 1 + bstar, 0), Or(And(Ne(zso, 1), cnd(zso)), And(Lt(re(mu + rho + q - p), 1), Eq(zso, 1))))
    c14 <- Or(c14, c14a)
  }
  c15 <- tryCatch({
    lc <- (q - p) * Abs(omega)^Rational(1, q - p) * cos(psi) + (v - u) * Abs(sigma)^Rational(1, v - u) * cos(theta)
    if (!identical(.mi_bool(.mi_eval_cond(Gt(lc, 0))), FALSE)) Gt(lc, 0)
    else {
      ls0 <- function(a, b) a * (q - p) * Abs(omega)^Rational(1, q - p) * sin(psi) + b * (v - u) * Abs(sigma)^Rational(1, v - u) * sin(theta)
      so <- sign(ua(omega)); ss <- sign(ua(sigma))
      lsx <- Piecewise(list(ls0(1, 1) * ls0(-1, -1), And(Eq(ua(sigma), 0), Eq(ua(omega), 0))),
                       list(ls0(so, 1) * ls0(so, -1), And(Eq(ua(sigma), 0), Ne(ua(omega), 0))),
                       list(ls0(1, ss) * ls0(-1, ss), And(Ne(ua(sigma), 0), Eq(ua(omega), 0))),
                       list(ls0(so, ss), S_true))
      Or(Gt(lc, 0), And(Eq(lc, 0), Ne(lsx, 0), Gt(re(eta), -1)), And(Eq(lc, 0), Eq(lsx, 0), Gt(re(eta), 0)))
    }
  }, error = function(e) S_false)
  A <- function(...) .mi_AND(list(...))
  conds <- list(
    A(.mi_tb(m * n * s * t != 0), pos(bstar), pos(cstar), c1, c2, c3, c10, c12),
    A(Eq(I(u), v), Eq(bstar, 0), pos(cstar), pos(sigma), Lt(re(rho), 1), c1, c2, c3, c12),
    A(Eq(I(p), q), Eq(cstar, 0), pos(bstar), pos(omega), Lt(re(mu), 1), c1, c2, c3, c10),
    A(Eq(I(p), q), Eq(I(u), v), Eq(bstar, 0), Eq(cstar, 0), pos(sigma), pos(omega), Lt(re(mu), 1), Lt(re(rho), 1), Ne(sigma, omega), c1, c2, c3),
    A(Eq(I(p), q), Eq(I(u), v), Eq(bstar, 0), Eq(cstar, 0), pos(sigma), pos(omega), Lt(re(mu + rho), 1), Ne(omega, sigma), c1, c2, c3),
    A(Gt(p, q), pos(I(s)), pos(bstar), Ge(cstar, 0), c1, c2, c3, c5, c10, c13),
    A(Lt(p, q), pos(I(t)), pos(bstar), Ge(cstar, 0), c1, c2, c3, c4, c10, c13),
    A(Gt(u, v), pos(I(m)), pos(cstar), Ge(bstar, 0), c1, c2, c3, c7, c11, c12),
    A(Lt(u, v), pos(I(n)), pos(cstar), Ge(bstar, 0), c1, c2, c3, c6, c11, c12),
    A(Gt(p, q), Eq(I(u), v), Eq(bstar, 0), Ge(cstar, 0), pos(sigma), Lt(re(rho), 1), c1, c2, c3, c5, c13),
    A(Lt(p, q), Eq(I(u), v), Eq(bstar, 0), Ge(cstar, 0), pos(sigma), Lt(re(rho), 1), c1, c2, c3, c4, c13),
    A(Eq(I(p), q), Gt(u, v), Ge(bstar, 0), Eq(cstar, 0), pos(omega), Lt(re(mu), 1), c1, c2, c3, c7, c11),
    A(Eq(I(p), q), Lt(u, v), Ge(bstar, 0), Eq(cstar, 0), pos(omega), Lt(re(mu), 1), c1, c2, c3, c6, c11),
    A(Lt(p, q), Gt(u, v), Ge(bstar, 0), Ge(cstar, 0), c1, c2, c3, c4, c7, c11, c13),
    A(Gt(p, q), Lt(u, v), Ge(bstar, 0), Ge(cstar, 0), c1, c2, c3, c5, c6, c11, c13),
    A(Gt(p, q), Gt(u, v), Ge(bstar, 0), Ge(cstar, 0), c1, c2, c3, c5, c7, c8, c11, c13, c14),
    A(Lt(p, q), Lt(u, v), Ge(bstar, 0), Ge(cstar, 0), c1, c2, c3, c4, c6, c9, c11, c13, c14),
    A(Eq(I(t), 0), pos(I(s)), pos(bstar), pos(phi), c1, c2, c10),
    A(Eq(I(s), 0), pos(I(t)), pos(bstar), neg(phi), c1, c3, c10),
    A(Eq(I(n), 0), pos(I(m)), pos(cstar), neg(phi), c1, c2, c12),
    A(Eq(I(m), 0), pos(I(n)), pos(cstar), pos(phi), c1, c3, c12),
    A(Eq(I(s * t), 0), pos(bstar), pos(cstar), c1, c2, c3, c10, c12),
    A(Eq(I(m * n), 0), pos(bstar), pos(cstar), c1, c2, c3, c10, c12))
  mt1 <- .mi_check_antecedents_1(g1, x, helper = TRUE); mt2 <- .mi_check_antecedents_1(g2, x, helper = TRUE)
  conds <- c(conds, list(
    A(mt2, Eq(I(t), 0), Lt(u, s), pos(bstar), c10, c1, c2, c3),
    A(mt2, Eq(I(s), 0), Lt(v, t), pos(bstar), c10, c1, c2, c3),
    A(mt1, Eq(I(n), 0), Lt(p, m), pos(cstar), c12, c1, c2, c3),
    A(mt1, Eq(I(m), 0), Lt(q, n), pos(cstar), c12, c1, c2, c3)))
  r <- .mi_OR(conds)
  if (!identical(.mi_bool(.mi_eval_cond(r)), FALSE)) return(r)
  conds <- c(conds, list(
    A(Gt(m + n, p), Eq(I(t), 0), Eq(phi, 0), pos(I(s)), pos(bstar), neg(cstar), Lt(Abs(ua(omega)), (m + n - p + 1) * PI), c1, c2, c10, c14, c15),
    A(Gt(m + n, q), Eq(I(s), 0), Eq(phi, 0), pos(I(t)), pos(bstar), neg(cstar), Lt(Abs(ua(omega)), (m + n - q + 1) * PI), c1, c3, c10, c14, c15),
    A(Eq(I(p), q - 1), Eq(I(t), 0), Eq(phi, 0), pos(I(s)), pos(bstar), Ge(cstar, 0), Lt(cstar * PI, Abs(ua(omega))), c1, c2, c10, c14, c15),
    A(Eq(I(p), q + 1), Eq(I(s), 0), Eq(phi, 0), pos(I(t)), pos(bstar), Ge(cstar, 0), Lt(cstar * PI, Abs(ua(omega))), c1, c3, c10, c14, c15),
    A(Lt(p, q - 1), Eq(I(t), 0), Eq(phi, 0), pos(I(s)), pos(bstar), Ge(cstar, 0), Lt(cstar * PI, Abs(ua(omega))), Lt(Abs(ua(omega)), (m + n - p + 1) * PI), c1, c2, c10, c14, c15),
    A(Gt(p, q + 1), Eq(I(s), 0), Eq(phi, 0), pos(I(t)), pos(bstar), Ge(cstar, 0), Lt(cstar * PI, Abs(ua(omega))), Lt(Abs(ua(omega)), (m + n - q + 1) * PI), c1, c3, c10, c14, c15),
    A(Eq(I(n), 0), Eq(phi, 0), Gt(s + t, 0), pos(I(m)), pos(cstar), neg(bstar), Lt(Abs(ua(sigma)), (s + t - u + 1) * PI), c1, c2, c12, c14, c15),
    A(Eq(I(m), 0), Eq(phi, 0), Gt(s + t, v), pos(I(n)), pos(cstar), neg(bstar), Lt(Abs(ua(sigma)), (s + t - v + 1) * PI), c1, c3, c12, c14, c15),
    A(Eq(I(n), 0), Eq(phi, 0), Eq(I(u), v - 1), pos(I(m)), pos(cstar), Ge(bstar, 0), Lt(bstar * PI, Abs(ua(sigma))), Lt(Abs(ua(sigma)), (bstar + 1) * PI), c1, c2, c12, c14, c15),
    A(Eq(I(m), 0), Eq(phi, 0), Eq(I(u), v + 1), pos(I(n)), pos(cstar), Ge(bstar, 0), Lt(bstar * PI, Abs(ua(sigma))), Lt(Abs(ua(sigma)), (bstar + 1) * PI), c1, c3, c12, c14, c15),
    A(Eq(I(n), 0), Eq(phi, 0), Lt(u, v - 1), pos(I(m)), pos(cstar), Ge(bstar, 0), Lt(bstar * PI, Abs(ua(sigma))), Lt(Abs(ua(sigma)), (s + t - u + 1) * PI), c1, c2, c12, c14, c15),
    A(Eq(I(m), 0), Eq(phi, 0), Gt(u, v + 1), pos(I(n)), pos(cstar), Ge(bstar, 0), Lt(bstar * PI, Abs(ua(sigma))), Lt(Abs(ua(sigma)), (s + t - v + 1) * PI), c1, c3, c12, c14, c15)))
  .mi_OR(conds)
}

.mi_int0oo <- function(g1, g2, x) {
  eta <- .mi_coeff_exp(g1$args[[3L]], x)[[1L]]; omega <- .mi_coeff_exp(g2$args[[3L]], x)[[1L]]
  ng <- function(l) lapply(l, .neg)
  a1 <- c(ng(.mg_bm(g1)), .mg_an(g2)); a2 <- c(.mg_aother(g2), ng(.mg_bother(g1)))
  b1 <- c(ng(.mg_an(g1)), .mg_bm(g2)); b2 <- c(.mg_bother(g2), ng(.mg_aother(g1)))
  .div(.mi_mg(a1, a2, b1, b2, .div(omega, eta)), eta)
}

.mi_rewrite_inversion <- function(fac, po, g, x) {
  s <- .mi_coeff_exp(po, x)[[2L]]
  ab <- .mi_coeff_exp(g$args[[3L]], x); a <- ab[[1L]]; b <- ab[[2L]]
  tr <- function(l) lapply(l, function(t) .add(list(t, .div(s, b))))
  list(powdenest(.div(fac, .pow(a, .div(s, b))), polar = TRUE),
       .mi_mg(tr(.mg_an(g)), tr(.mg_aother(g)), tr(.mg_bm(g)), tr(.mg_bother(g)), g$args[[3L]]))
}

.mi_check_antecedents_inversion <- function(g, x) {
  z <- g$args[[3L]]
  e <- .mi_coeff_exp(z, x)[[2L]]
  if (.is_true(Lt(e, 0))) return(.mi_check_antecedents_inversion(.mi_flip_g(g), x))
  PI <- S_Pi; I <- sym_I
  statement_half <- function(a, b, c, zz, plus) {
    ce <- .mi_coeff_exp(zz, x); coeff <- ce[[1L]]; ex <- ce[[2L]]
    a <- a * ex; b <- b * coeff^c; c <- c * ex
    w <- if (plus) b * exp(I * re(c) * PI / 2) else b * exp(-I * re(c) * PI / 2)
    Or(And(Or(Eq(b, 0), Le(re(c), 0)), Le(re(a), -1)),
       And(Ne(b, 0), Eq(im(c), 0), Gt(re(c), 0), Lt(re(w), 0)),
       And(Ne(b, 0), Eq(im(c), 0), Gt(re(c), 0), Le(re(w), 0), Le(re(a), -1)))
  }
  statement <- function(a, b, c, zz) And(statement_half(a, b, c, zz, TRUE), statement_half(a, b, c, zz, FALSE))
  m <- length(.mg_bm(g)); n <- length(.mg_an(g)); p <- length(.mg_ap(g)); q <- length(.mg_bq(g))
  tau <- m + n - p; nu <- q - m - n; rho <- Rational(tau - nu, 2); sigma <- q - p
  epsilon <- if (sigma == 1) S_Half else if (sigma > 1) S_One else S_NaN
  theta <- (Rational(1 - sigma, 2) + .add(.mg_bq(g)) - .add(.mg_ap(g))) / sigma
  delta <- .mg_delta(g)
  if (!(.is_true(Ge(delta, e / 2)) || (p >= 1 && p >= q))) return(S_false)
  for (a in .mg_an(g)) for (b in .mg_bm(g)) if (isTRUE(ask_fact(.sub(a, b), "integer")) && .is_true(Gt(a, b))) return(S_false)
  E <- function(zz) .mi_AND(lapply(.mg_an(g), function(a) statement(a - 1, 0, 0, zz)))
  if (p >= q) return(E(z))
  H <- function(zz) statement(theta, -sigma, 1 / Integer(sigma), zz)
  Hp <- function(zz) statement_half(theta, -sigma, 1 / Integer(sigma), zz, TRUE)
  Hm <- function(zz) statement_half(theta, -sigma, 1 / Integer(sigma), zz, FALSE)
  conds <- list(
    And(Le(1, n), Le(1, m), Ge(rho * PI - delta, PI / 2), Gt(delta, 0), E(z * exp(I * PI * (nu + 1)))),
    And(Le(p + 1, m), Le(m + 1, q), Gt(delta, 0), Lt(delta, PI / 2), .mi_tb(n == 0), Ge((m - p + 1) * PI - delta, PI / 2),
        Hp(z * exp(I * PI * (q - m))), Hm(z * exp(-I * PI * (q - m)))),
    And(.mi_tb(m == q), .mi_tb(n == 0), Gt(delta, 0), Ge((sigma + epsilon) * PI - delta, PI / 2), H(z)),
    And(Or(And(Le(p, q - 2), Le(1, tau), Le(tau, Rational(sigma, 2))), And(Le(p + 1, m + n), Le(m + n, Rational(p + q, 2)))),
        Gt(delta, 0), Lt(delta, PI / 2), Ge((tau + 1) * PI - delta, PI / 2), Hp(z * exp(I * PI * nu)), Hm(z * exp(-I * PI * nu))),
    And(Le(1, m), Gt(rho, 0), Gt(delta, 0), Lt(delta + rho * PI, PI / 2), Ge((tau + epsilon) * PI - delta, PI / 2),
        Hp(z * exp(I * PI * nu)), Hm(z * exp(-I * PI * nu))),
    .mi_tb(m == 0))
  .mi_OR(conds)
}

.mi_int_inversion <- function(g, x, t) {
  ba <- .mi_coeff_exp(g$args[[3L]], x); b <- ba[[1L]]; a <- ba[[2L]]
  r <- .mi_inflate_fox_h(.mi_mg(.mg_an(g), .mg_aother(g), .mg_bm(g), .mg_bother(g), .div(b, .pow(t, a))), .neg(a))
  .mul(list(.div(r[[1L]], t), r[[2L]]))
}

# _rewrite_single: f as a sum of c * x^s * G(a * x^b), with a condition
.mi_rewrite_single <- function(f, x, recursive = TRUE) {
  if (f$type == "Function" && f$name == "meijerg") {
    cm <- .as_coeff_mul_dep(sym_factor(f$args[[3L]]), x)
    coeff <- cm[[1L]]; m <- cm[[2L]]
    if (length(m) > 1L) return(NULL)
    m <- m[[1L]]
    if (m$type == "Pow") { if (!identical(m$args[[1L]]$key, x$key) || !(m$args[[2L]]$type %in% c("Integer", "Rational"))) return(NULL) }
    else if (!identical(m$key, x$key)) return(NULL)
    return(list(list(list(S_One, S_Zero, .mi_mg(.mg_an(f), .mg_aother(f), .mg_bm(f), .mg_bother(f), .mul(list(coeff, m))))), S_true))
  }
  z <- .mi_z()
  f_ <- f
  f <- subs(f, x, z)
  tab <- .mi_table()
  ent <- tab[[.mi_typekey(.mi_mytype(f, z))]]
  for (row in ent) {
    sb <- tryCatch(sym_match(f, row$formula, old = TRUE), error = function(e) NULL)
    if (is.null(sb) || !length(sb)) next
    sb <- lapply(sb, function(to) unpolarify(polarify(to, lift = TRUE), exponents_only = TRUE))
    reps <- lapply(names(sb), function(k) list(.mi_state$wilds[[k]], sb[[k]]))
    hint <- subs(row$hint, reps)
    if (identical(hint$key, S_false$key)) next
    cond <- row$cond
    if (!(cond$type == "BooleanAtom")) cond <- unpolarify(subs(cond, reps))
    if (identical(.mi_bool(.mi_eval_cond(cond)), FALSE)) next
    terms <- if (is.function(row$inst)) row$inst(sb) else row$inst
    res <- list()
    for (tm in terms) {
      r1 <- tryCatch(.mi_coeff_exp(unpolarify(subs(subs(tm[[1L]], reps), z, x), exponents_only = TRUE), x), error = function(e) NULL)
      if (is.null(r1)) next
      g <- tryCatch(subs(subs(tm[[2L]], reps), z, x), error = function(e) NULL)
      if (is.null(g)) next
      bad <- c(S_Infinity$key, S_ComplexInfinity$key, S_NegativeInfinity$key)
      if (any(vapply(list(r1[[1L]], r1[[2L]], g), function(e) any(vapply(bad, function(k) .contains_key(e, k), NA)), NA))) next
      if (!(g$type == "Function" && g$name == "meijerg")) next
      g <- .mi_mg(.mg_an(g), .mg_aother(g), .mg_bm(g), .mg_bother(g), unpolarify(g$args[[3L]], exponents_only = TRUE))
      res[[length(res) + 1L]] <- list(r1[[1L]], r1[[2L]], g)
    }
    if (length(res)) return(list(res, cond))
  }
  if (!recursive) return(NULL)
  my_imt <- function(F, s, x, strip) {
    tryCatch(inverse_mellin_transform(F, s, x, strip, as_meijerg = TRUE, needeval = TRUE),
             MellinTransformStripError = function(e)
               inverse_mellin_transform(simplify(cancel(expand(F))), s, x, strip, as_meijerg = TRUE, needeval = TRUE))
  }
  f <- f_
  s <- .mi_dummy("s", "rewrite-single", f)
  my_integrator <- function(f, x) {
    r <- .mi_definite_4(f, x, only_double = TRUE)
    if (!is.null(r)) {
      res <- .mi_my_unpolarify(hyperexpand(r[[1L]], rewrite = "nonrepsmall"))
      return(Piecewise(list(res, r[[2L]]), list(Integral(f, list(x, S_Zero, S_Infinity)), S_true)))
    }
    Integral(f, list(x, S_Zero, S_Infinity))
  }
  g <- tryCatch({
    mt <- mellin_transform(f, x, s, integrator = my_integrator, simplify = FALSE, needeval = TRUE)
    my_imt(mt[[1L]], s, x, mt[[2L]])
  }, IntegralTransformError = function(e) NULL)
  if (is.null(g)) {
    a <- .mi_dummy_("a", "rewrite-single")
    if (!.has1(f, a) && .mi_is_analytic(f, x)) {
      g <- tryCatch({
        mt <- mellin_transform(subs(f, x, .mul(list(a, x))), x, s, integrator = my_integrator, needeval = TRUE, simplify = FALSE)
        subs(my_imt(mt[[1L]], s, x, mt[[2L]]), a, S_One)
      }, IntegralTransformError = function(e) NULL)
    }
  }
  if (is.null(g) || any(vapply(c(S_Infinity$key, S_NaN$key, S_ComplexInfinity$key), function(k) .contains_key(g, k), NA))) return(NULL)
  res <- list()
  for (fz in (if (g$type == "Add") .sympy_args(g) else list(g))) {
    cm <- .as_coeff_mul_dep(fz, x)
    if (length(cm[[2L]]) > 1L) .notimpl_error("Unexpected form...")
    gg <- cm[[2L]][[1L]]
    ab <- .mi_coeff_exp(gg$args[[3L]], x)
    res[[length(res) + 1L]] <- list(cm[[1L]], S_Zero, .mi_mg(.mg_an(gg), .mg_aother(gg), .mg_bm(gg), .mg_bother(gg),
      .mul(list(unpolarify(polarify(ab[[1L]], lift = TRUE), exponents_only = TRUE), .pow(x, ab[[2L]])))))
  }
  list(res, S_true)
}

.mi_rewrite1 <- function(f, x, recursive = TRUE) {
  sp <- .mi_split_mul(f, x)
  g <- .mi_rewrite_single(sp[[3L]], x, recursive)
  if (!is.null(g)) list(sp[[1L]], sp[[2L]], g[[1L]], g[[2L]]) else NULL
}

.mi_rewrite2 <- function(f, x) {
  sp <- .mi_split_mul(f, x); fac <- sp[[1L]]; po <- sp[[2L]]; g <- sp[[3L]]
  if (any(vapply(.mi_mul_args(g), function(e) is.null(.mi_rewrite_single(e, x, FALSE)), NA))) return(NULL)
  l <- .mi_mul_as_two_parts(g)
  if (!length(l)) return(NULL)
  k1 <- vapply(l, function(pp) max(length(.mi_exponents(pp[[1L]], x)), length(.mi_exponents(pp[[2L]], x))), 0)
  k2 <- vapply(l, function(pp) max(length(.mi_functions(pp[[1L]], x)), length(.mi_functions(pp[[2L]], x))), 0)
  k3 <- vapply(l, function(pp) max(length(.mi_find_splitting_points(pp[[1L]], x)), length(.mi_find_splitting_points(pp[[2L]], x))), 0)
  # ordered(l, keys): the keys, then default_sort_key of the pairs
  nk <- vapply(l, function(pp) 1 + .nodes(pp[[1L]]) + .nodes(pp[[2L]]), 0)
  dk <- lapply(l, function(pp) list(.sort_key(pp[[1L]]), .sort_key(pp[[2L]])))
  ord <- .key_sort(Map(function(a, b, c, n, d) list(a, b, c, n, d), k1, k2, k3, nk, dk), seq_along(l))
  l <- l[unlist(ord)]
  for (rec in c(FALSE, TRUE)) for (pp in l) {
    g1 <- .mi_rewrite_single(pp[[1L]], x, rec); g2 <- .mi_rewrite_single(pp[[2L]], x, rec)
    if (!is.null(g1) && !is.null(g2)) {
      cond <- And(g1[[2L]], g2[[2L]])
      if (!identical(cond$key, S_false$key)) return(list(fac, po, g1[[1L]], g2[[1L]], cond))
    }
  }
  NULL
}

.mi_has <- function(res, names) {
  res <- piecewise_fold(res)
  if (res$type == "Piecewise") return(all(vapply(res$args, function(pc) .has_function(pc, names), NA)))
  .has_function(res, names)
}

.mi_first_ordered <- function(results) {
  keys <- lapply(results, function(r) {
    if (is.list(r) && !is_sym(r)) list(sum(vapply(r, function(z) .nodes(as_sym(z)), 0)) + 1, lapply(r, function(z) .sort_key(as_sym(z))))
    else list(.nodes(r), .sort_key(r))
  })
  .key_sort(keys, results)[[1L]]
}

.mi_my_unpolarify <- function(f) .mi_eval_cond(unpolarify(f))

.mi_rewrite_hyperbolics_as_exp <- function(e) .map_bottom_up(e, function(z)
  if (z$type == "Function" && z$name %in% .HYP_FUNS) rewrite(z, "exp") else z)

#' Integration by Meijer G functions
#'
#' Ports of SymPy's `meijerint_indefinite`, `meijerint_definite` and
#' `meijerint_inversion` (sympy/integrals/meijerint.py). The integrand is
#' rewritten as a product of Meijer G functions, which are integrated with
#' the formulas of Prudnikov, Brychkov and Marichev and then expanded back
#' to elementary and special functions with [hyperexpand()].
#'
#' `meijerint_indefinite(f, x)` returns an antiderivative (possibly a
#' Piecewise) or `NULL`. `meijerint_definite(f, x, a, b)` returns a list of
#' the value and its convergence condition, or `NULL`.
#' `meijerint_inversion(f, x, t)` computes the inverse Laplace transform of
#' `f`, or `NULL`.
#'
#' @param f The integrand.
#' @param x The integration variable.
#' @param a,b The limits.
#' @param t The variable of the inverse Laplace transform.
#' @return An expression, a list `(value, condition)`, or `NULL`.
#' @examples
#' x <- Symbol("x")
#' meijerint_indefinite(sin(x), x)
#' meijerint_definite(exp(-x^2), x, 0, oo)
#' @export
meijerint_indefinite <- function(f, x) {
  f <- as_sym(f); x <- as_sym(x)
  pts <- c(.mi_find_splitting_points(f, x), list(S_Zero))
  pts <- pts[!duplicated(vapply(pts, `[[`, "", "key"))]
  pts <- .key_sort(lapply(pts, .sort_key), pts)
  results <- list()
  for (a in pts) {
    res <- .mi_indefinite_1(subs(f, x, .add(list(x, a))), x)
    if (is.null(res)) next
    res <- subs(res, x, .sub(x, a))
    if (.mi_has(res, c("hyper", "meijerg"))) results[[length(results) + 1L]] <- res
    else return(res)
  }
  if (.has_function(f, .HYP_FUNS)) {
    rv <- meijerint_indefinite(.mi_rewrite_hyperbolics_as_exp(f), x)
    # collect(factor_terms(rv), rv.atoms(exp)) (sympy/integrals/meijerint.py:1684)
    if (!is.null(rv)) return(collect(factor_terms(rv), Filter(function(f) f$name == "exp", atoms(rv, "Function"))))
  }
  if (length(results)) .mi_first_ordered(results) else NULL
}

.mi_clean <- function(res, x) {
  res <- expand_mul(cancel(res), deep = FALSE)
  terms <- if (res$type == "Add") res$args else list(res)
  .add(Filter(function(t) .has1(t, x), terms))
}

.mi_indefinite_1 <- function(f, x) {
  gs <- .mi_rewrite1(f, x)
  if (is.null(gs)) return(NULL)
  fac <- gs[[1L]]; po <- gs[[2L]]; gl <- gs[[3L]]; cond <- gs[[4L]]
  res <- S_Zero
  t <- .mi_dummy("t", "meijerint-indefinite", S_One)
  for (tm in gl) {
    C <- tm[[1L]]; s <- tm[[2L]]; g <- tm[[3L]]
    ab <- .mi_coeff_exp(g$args[[3L]], x); a <- ab[[1L]]; b <- ab[[2L]]
    c <- .add(list(.mi_coeff_exp(po, x)[[2L]], s))
    fac_ <- .div(.mul(list(fac, C, .pow(x, .add(list(S_One, c))))), b)
    rho <- .div(.add(list(c, S_One)), b)
    trbm <- lapply(.mg_bm(g), function(z) .add(list(z, rho)))
    if (any(vapply(trbm, function(z) isTRUE(ask_fact(z, "integer")) && .is_true(Le(z, 0)), NA)))
      r <- .neg(.mi_mg(.mg_an(g), c(.mg_aother(g), list(.sub(S_One, rho))), c(.mg_bm(g), list(.neg(rho))), .mg_bother(g), t))
    else
      r <- .mi_mg(c(.mg_an(g), list(.sub(S_One, rho))), .mg_aother(g), .mg_bm(g), c(.mg_bother(g), list(.neg(rho))), t)
    f0 <- subs(f, x, S_Zero)
    place <- if (isTRUE(ask_fact(b, "extended_nonnegative")) && !.contains_key(f0, S_NaN$key) && !.contains_key(f0, S_ComplexInfinity$key)) 0 else NULL
    r <- hyperexpand(subs(r, t, .mul(list(a, .pow(x, b)))), place = place)
    res <- .add(list(res, powdenest(.mul(list(fac_, r)), polar = TRUE)))
  }
  res <- piecewise_fold(res, evaluate = NULL)
  if (res$type == "Piecewise") {
    res <- .piecewise_noeval(lapply(res$args, function(pc) list(.mi_my_unpolarify(.mi_clean(pc$args[[1L]], x)), pc$args[[2L]])))
  } else res <- .mi_my_unpolarify(.mi_clean(res, x))
  Piecewise(list(res, .mi_my_unpolarify(cond)), list(Integral(f, x), S_true))
}

#' @rdname meijerint_indefinite
#' @export
meijerint_definite <- function(f, x, a, b) {
  f <- as_sym(f); x <- as_sym(x); a <- as_sym(a); b <- as_sym(b)
  if (.has_function(f, c("DiracDelta", "SingularityFunction"))) return(NULL)
  f_ <- f; x_ <- x; a_ <- a; b_ <- b
  d <- Dummy("x")
  f <- subs(f, x, d); x <- d
  if (identical(a$key, b$key)) return(list(S_Zero, S_true))
  results <- list()
  isinf <- function(e, n) is_const(e, n)
  if (isinf(a, "-oo") && !isinf(b, "oo")) return(meijerint_definite(subs(f, x, .neg(x)), x, .neg(b), .neg(a)))
  if (isinf(a, "-oo")) {
    inner <- .mi_find_splitting_points(f, x)
    inner <- rev(.key_sort(lapply(inner, .sort_key), inner))
    for (c in c(inner, list(S_Zero))) {
      if (!isTRUE(ask_fact(c, "extended_real"))) next
      r1 <- .mi_definite_2(subs(f, x, .add(list(x, c))), x)
      if (is.null(r1)) next
      r2 <- .mi_definite_2(subs(f, x, .sub(c, x)), x)
      if (is.null(r2)) next
      cond <- .mi_condsimp(And(r1[[2L]], r2[[2L]]))
      if (identical(.mi_bool(cond), FALSE)) next
      return(list(.add(list(r1[[1L]], r2[[1L]])), cond))
    }
  } else if (isinf(a, "oo")) {
    res <- meijerint_definite(f, x, b, S_Infinity)
    if (is.null(res)) return(NULL)
    return(list(.neg(res[[1L]]), res[[2L]]))
  } else if (is_zero_num(a) && isinf(b, "oo")) {
    res <- .mi_definite_2(f, x)
    if (!is.null(res)) {
      if (.mi_has(res[[1L]], "meijerg")) results[[length(results) + 1L]] <- res else return(res)
    }
  } else {
    if (isinf(b, "oo")) for (sp in .mi_find_splitting_points(f, x)) {
      if (.is_true(Ge(.sub(a, sp), 0))) {
        res <- .mi_definite_2(.mul(list(subs(f, x, .add(list(x, sp))), Heaviside(.sub(.add(list(x, sp)), a)))), x)
        if (!is.null(res)) {
          if (.mi_has(res[[1L]], "meijerg")) results[[length(results) + 1L]] <- res else return(res)
        }
      }
    }
    f <- subs(f, x, .add(list(x, a)))
    b <- .sub(b, a); a <- S_Zero
    if (!isinf(b, "oo")) {
      phi <- exp(sym_I * fcall("arg", list(b)))
      b <- Abs(b)
      f <- subs(f, x, .mul(list(phi, x)))
      f <- .mul(list(f, Heaviside(.sub(b, x)), phi))
      b <- S_Infinity
    }
    res <- .mi_definite_2(f, x)
    if (!is.null(res)) {
      if (.mi_has(res[[1L]], "meijerg")) results[[length(results) + 1L]] <- res else return(res)
    }
  }
  if (.has_function(f_, .HYP_FUNS)) {
    rv <- meijerint_definite(.mi_rewrite_hyperbolics_as_exp(f_), x_, a_, b_)
    # (sympy/integrals/meijerint.py:1912)
    if (!is.null(rv)) return(list(collect(factor_terms(rv[[1L]]), Filter(function(f) f$name == "exp", atoms(rv[[1L]], "Function"))), rv[[2L]]))
  }
  if (length(results)) .mi_first_ordered(results) else NULL
}

.mi_guess_expansion <- function(f, x) {
  res <- list(f); saw <- f$key
  addc <- function(e) if (!(e$key %in% saw)) { res[[length(res) + 1L]] <<- e; saw <<- c(saw, e$key) }
  addc(expand_mul(f))
  addc(expand(f))
  if (.has_function(f, c("sin", "cos", "tan", "cot", "sec", "csc", .HYP_FUNS))) addc(expand_mul(expand_trig(f)))
  if (.has_function(f, c("sin", "cos"))) addc(TR8(expand_mul(TRpower(f))))
  res
}

.mi_definite_2 <- function(f, x) {
  dummy <- .mi_dummy("x", "meijerint-definite2", f, positive = TRUE)
  f <- subs(f, x, dummy); x <- dummy
  if (is_zero_num(f)) return(list(S_Zero, S_true))
  for (g in .mi_guess_expansion(f, x)) {
    res <- .mi_definite_3(g, x)
    if (!is.null(res)) return(res)
  }
  NULL
}

.mi_definite_3 <- function(f, x) {
  res <- .mi_definite_4(f, x)
  if (!is.null(res) && !identical(.mi_bool(res[[2L]]), FALSE)) return(res)
  if (f$type == "Add") {
    ress <- lapply(.sympy_args(f), .mi_definite_4, x = x)
    if (all(!vapply(ress, is.null, NA))) {
      c <- do.call(And, lapply(ress, `[[`, 2L))
      if (!identical(.mi_bool(c), FALSE)) return(list(.add(lapply(ress, `[[`, 1L)), c))
    }
  }
  NULL
}

.mi_definite_4 <- function(f, x, only_double = FALSE) {
  if (!only_double) {
    gs <- .mi_rewrite1(f, x, recursive = FALSE)
    if (!is.null(gs)) {
      fac <- gs[[1L]]; po <- gs[[2L]]; cond <- gs[[4L]]
      res <- S_Zero
      for (tm in gs[[3L]]) {
        C <- tm[[1L]]
        if (is_zero_num(C)) next
        r <- .mi_rewrite_saxena_1(.mul(list(fac, C)), .mul(list(po, .pow(x, tm[[2L]]))), tm[[3L]], x)
        res <- .add(list(res, .mul(list(r[[1L]], .mi_int0oo_1(r[[2L]], x)))))
        cond <- And(cond, .mi_check_antecedents_1(r[[2L]], x))
        if (identical(.mi_bool(cond), FALSE)) break
      }
      cond <- .mi_my_unpolarify(cond)
      if (!identical(.mi_bool(cond), FALSE)) return(list(.mi_my_unpolarify(hyperexpand(res)), cond))
    }
  }
  gs <- .mi_rewrite2(f, x)
  if (is.null(gs)) return(NULL)
  for (full_pb in c(FALSE, TRUE)) {
    fac <- gs[[1L]]; po <- gs[[2L]]; g1 <- gs[[3L]]; g2 <- gs[[4L]]; cond <- gs[[5L]]
    res <- S_Zero; stop_all <- FALSE
    for (t1 in g1) {
      for (t2 in g2) {
        r <- .mi_rewrite_saxena(.mul(list(fac, t1[[1L]], t2[[1L]])), .mul(list(po, .pow(x, .add(list(t1[[2L]], t2[[2L]]))))), t1[[3L]], t2[[3L]], x, full_pb)
        if (is.null(r)) return(NULL)
        cond <- And(cond, .mi_check_antecedents(r[[2L]], r[[3L]], x))
        if (identical(.mi_bool(cond), FALSE)) { stop_all <- TRUE; break }
        res <- .add(list(res, .mul(list(r[[1L]], .mi_int0oo(r[[2L]], r[[3L]], x)))))
      }
      if (stop_all) break
    }
    cond <- .mi_my_unpolarify(cond)
    if (!identical(.mi_bool(cond), FALSE)) {
      if (only_double) return(list(res, cond))
      return(list(.mi_my_unpolarify(hyperexpand(res)), cond))
    }
  }
  NULL
}

#' @rdname meijerint_indefinite
#' @export
meijerint_inversion <- function(f, x, t) {
  f <- as_sym(f); x <- as_sym(x); t_ <- as_sym(t)
  f_ <- f
  t <- Dummy("t", polar = TRUE)
  f <- subs(f, t_, t)
  if (!.mi_is_analytic(f, x)) return(NULL)
  shift <- S_Zero
  args <- if (f$type == "Mul") .sympy_args(f) else if (f$type == "Function" && f$name == "exp") list(f) else NULL
  if (!is.null(args)) {
    newargs <- list(); expo <- list()
    while (length(args)) {
      arg <- args[[length(args)]]; args[[length(args)]] <- NULL
      if (arg$type == "Function" && arg$name == "exp") {
        arg2 <- expand(arg)
        if (arg2$type == "Mul") { args <- c(args, .sympy_args(arg2)); next }
        ab <- tryCatch(.mi_coeff_exp(arg$args[[1L]], x), mi_coeff_exp_error = function(e) list(NULL, S_Zero))
        if (is_one_num(ab[[2L]])) expo[[length(expo) + 1L]] <- ab[[1L]] else newargs[[length(newargs) + 1L]] <- arg
      } else if (arg$type == "Pow") {
        arg2 <- expand(arg)
        if (arg2$type == "Mul") { args <- c(args, .sympy_args(arg2)); next }
        if (!.has1(arg$args[[1L]], x)) {
          ab <- tryCatch(.mi_coeff_exp(arg$args[[2L]], x), mi_coeff_exp_error = function(e) list(NULL, S_Zero))
          if (is_one_num(ab[[2L]])) expo[[length(expo) + 1L]] <- .mul(list(ab[[1L]], log(arg$args[[1L]])))
        }
        newargs[[length(newargs) + 1L]] <- arg
      } else newargs[[length(newargs) + 1L]] <- arg
    }
    shift <- .add(expo)
    f <- .mul(newargs)
  }
  if (!.has1(f, x)) {
    cond <- Eq(im(shift), 0)
    if (identical(.mi_bool(cond), FALSE)) return(NULL)
    res <- .mul(list(f, DiracDelta(.add(list(t, shift)))))
    return(Piecewise(list(subs(res, t, t_), cond)))
  }
  gs <- .mi_rewrite1(f, x)
  if (is.null(gs)) return(NULL)
  fac <- gs[[1L]]; po <- gs[[2L]]; cond <- gs[[4L]]
  res <- S_Zero
  for (tm in gs[[3L]]) {
    r <- .mi_rewrite_inversion(.mul(list(fac, tm[[1L]])), .mul(list(po, .pow(x, tm[[2L]]))), tm[[3L]], x)
    res <- .add(list(res, .mul(list(r[[1L]], .mi_int_inversion(r[[2L]], x, t)))))
    cond <- And(cond, .mi_check_antecedents_inversion(r[[2L]], x))
    if (identical(.mi_bool(cond), FALSE)) break
  }
  cond <- .mi_my_unpolarify(cond)
  if (identical(.mi_bool(cond), FALSE)) return(NULL)
  res <- .mi_my_unpolarify(hyperexpand(res))
  if (!.has_function(res, "Heaviside")) res <- .mul(list(res, Heaviside(t)))
  res <- subs(res, t, .add(list(t, shift)))
  if (!(cond$type == "BooleanAtom")) cond <- subs(cond, t, .add(list(t, shift)))
  Piecewise(list(subs(res, t, t_), cond), list(InverseLaplaceTransform(subs(f_, t, t_), x, t_, NULL), S_true))
}
