# Inverse trigonometric functions of trigonometric values (the
# isinstance(arg, cos) branches of acos.eval, atan.eval and acot.eval in
# sympy/functions/elementary/trigonometric.py) and the Fibonacci and Lucas
# numbers and polynomials of sympy/functions/combinatorial/numbers.py.

# a > b for comparable numbers, as the Relational truth value SymPy uses
.fs_gt <- function(u, v) isTRUE(Re(.neval(.sub(u, v), list())) > 0)

# acos.eval: acos(cos(ang)), acos(-cos(ang)), acos(sin(ang)), acos(-sin(ang))
.fs_acos_trig <- function(a) {
  sa <- .sympy_args(a)
  if (a$type == "Mul" && length(sa) == 2L && is_negone_num(sa[[1L]])) {
    narg <- sa[[2L]]
    minus <- TRUE
  } else {
    narg <- a
    minus <- FALSE
  }
  if (narg$type != "Function" || !(narg$name %in% c("cos", "sin"))) return(NULL)
  ang <- narg$args[[1L]]
  if (!.is_comparable(ang)) return(NULL)
  hp <- .mul(list(S_Half, S_Pi))
  if (narg$name == "cos") {
    if (minus) ang <- .sub(S_Pi, ang)
    tp <- .mul(list(S_Two, S_Pi))
    ang <- Mod(ang, tp)
    if (.fs_gt(ang, S_Pi)) ang <- .sub(tp, ang)
    return(ang)
  }
  if (minus) return(.add(list(hp, fcall("asin", list(narg)))))
  .sub(hp, fcall("asin", list(narg)))
}

# atan.eval (name "atan", matching tan) and acot.eval (name "acot",
# matching cot): the angle reduced into (-pi/2, pi/2]
.fs_atan_trig <- function(name, a) {
  if (a$type != "Function" || !(a$name %in% c("tan", "cot"))) return(NULL)
  ang <- a$args[[1L]]
  if (!.is_comparable(ang)) return(NULL)
  hp <- .mul(list(S_Half, S_Pi))
  same <- if (name == "atan") "tan" else "cot"
  if (a$name == same) {
    ang <- Mod(ang, S_Pi)
  } else {
    other <- if (name == "atan") "acot" else "atan"
    ang <- .sub(hp, fcall(other, list(a)))
  }
  if (.fs_gt(ang, hp)) ang <- .sub(ang, S_Pi)
  ang
}

# fibonacci._fibpoly: F1 = 1, F2 = x, Fn = (F(n-2) + x*F(n-1)).expand()
.fs_fibpoly <- function(n, x) {
  prev <- list(S_One, x)
  if (n <= 2) return(prev[[n]])
  for (i in 3:n) prev <- list(prev[[2L]], expand(.add(list(prev[[1L]], .mul(list(x, prev[[2L]]))))))
  prev[[2L]]
}

.fs_fib_int <- function(n) {
  a <- gmp::as.bigz(0); b <- gmp::as.bigz(1)
  if (n == 0) return(S_Zero)
  for (i in seq_len(n - 1)) { t <- a + b; a <- b; b <- t }
  make_rational(gmp::as.bigq(b))
}

# fibonacci.eval(n, sym=None)
.fs_fibonacci_eval <- function(args) {
  n <- args[[1L]]
  if (is_const(n, "oo")) return(S_Infinity)
  if (n$type != "Integer") return(NULL)
  k <- .int_val(n)
  if (length(args) == 1L) {
    if (abs(k) >= 5000) return(NULL)
    if (k < 0) return(.mul(list(.pow(S_NegativeOne, .int(k + 1)), .fs_fib_int(-k))))
    return(.fs_fib_int(k))
  }
  if (k < 1) stop(.sympy_condition("ValueError", "Fibonacci polynomials are defined only for positive integer indices."))
  d <- Dummy("x")
  subs(.fs_fibpoly(k, d), d, args[[2L]])
}

# lucas.eval: fibonacci(n + 1) + fibonacci(n - 1)
.fs_lucas_eval <- function(args) {
  n <- args[[1L]]
  if (is_const(n, "oo")) return(S_Infinity)
  if (n$type != "Integer" || abs(.int_val(n)) >= 5000) return(NULL)
  k <- .int_val(n)
  .add(list(.fs_fibonacci_eval(list(.int(k + 1))), .fs_fibonacci_eval(list(.int(k - 1)))))
}

# Mul.flatten's neg1e: the imaginary unit counts as (-1)^(1/2) in the total
# exponent of -1 collected from positive rational powers of -1; an odd
# integer part becomes the factor -1, a half becomes I, and any other
# remainder stays (-1)^(p/q) for the join with a positive base.
.fs_neg1e <- function(fs) {
  ii <- which(vapply(fs, function(f) identical(f$type, "Constant") && identical(f$name, "I"), NA))
  if (!length(ii)) return(fs)
  neg <- which(vapply(fs, function(f) f$type == "Pow" && f$args[[1L]]$type %in% c("Integer", "Rational") &&
                        f$args[[1L]]$value < 0 && f$args[[2L]]$type == "Rational" &&
                        f$args[[2L]]$value > 0, TRUE))
  if (!length(ii) || !length(neg)) return(fs)
  e <- gmp::as.bigq(length(ii), 2)
  # a negative base b gives (-1)^e and (-b)^e
  pos <- list()
  for (k in neg) {
    e <- e + fs[[k]]$args[[2L]]$value
    b <- fs[[k]]$args[[1L]]
    if (!is_negone_num(b)) pos[[length(pos) + 1L]] <- .pow(make_rational(-b$value), fs[[k]]$args[[2L]])
  }
  fs <- c(fs[-c(ii, neg)], pos)
  p <- gmp::numerator(e)
  q <- gmp::denominator(e)
  r <- p %% q
  n <- (p - r) %/% q
  if (n %% 2 != 0) fs[[length(fs) + 1L]] <- S_NegativeOne
  if (q == 2) {
    fs[[length(fs) + 1L]] <- S_I
  } else if (r != 0) {
    fs[[length(fs) + 1L]] <- new_sym("Pow", list(S_NegativeOne, make_rational(gmp::as.bigq(r, q))))
  }
  fs
}

# DiracDelta._eval_expand_diracdelta: DiracDelta(g(x)) as the sum of
# DiracDelta(x - r)/|g'(r)| over the simple real roots r of g, with respect
# to the hint wrt (or the only free symbol)
.fs_expand_diracdelta <- function(e, h = list()) {
  if (!(e$type == "Function" && e$name == "DiracDelta")) return(NULL)
  wrt <- h$wrt
  if (is.null(wrt)) {
    free <- free_symbols(e)
    if (length(free) != 1L)
      stop(.sympy_condition("TypeError", paste0("\nWhen there is more than 1 free symbol or variable in the expression,\n",
                                                "the 'wrt' keyword is required as a hint to expand when using the\nDiracDelta hint.")))
    wrt <- free[[1L]]
  }
  wrt <- as_sym(wrt)
  a <- e$args[[1L]]
  if (!.has1(a, wrt) || (length(e$args) > 1L && !is_zero_num(e$args[[2L]]))) return(e)
  rts <- tryCatch(roots(a, wrt), error = function(err) NULL)
  if (is.null(rts)) return(e)
  mult <- attr(rts, "multiplicity")
  darg <- fcall("Abs", list(diff(a, wrt)))
  terms <- list()
  for (k in seq_along(rts)) {
    r <- rts[[k]]
    if (identical(ask_fact(r, "real"), FALSE) || mult[k] != 1L) return(e)
    terms[[k]] <- .div(fcall("DiracDelta", list(.sub(wrt, r))), subs(darg, wrt, r))
  }
  .add(terms)
}
.expand_custom$diracdelta <- .fs_expand_diracdelta

# heurisch_wrapper: solve([Eq(k, v) for every case], dict=True, exclude=(x,))
.fs_heur_solve_system <- function(slns, x) {
  eqs <- unlist(lapply(slns, function(sd) lapply(sd, function(kv) Eq(kv[[1L]], kv[[2L]]))), recursive = FALSE)
  eqs <- Filter(function(q) !.is_true(q), eqs)
  if (!length(eqs)) return(list())
  syms <- .poly_sort_gens(Filter(function(s) !identical(s$key, x$key), free_symbols(do.call(Tuple, eqs))))
  if (!length(syms)) return(list())
  sol <- tryCatch(solve(eqs, syms, dict = TRUE), error = function(e) list())
  if (!is.list(sol) || !length(sol)) return(list())
  if (!is.null(names(sol)) && all(nzchar(names(sol)))) sol <- list(sol)
  bynm <- stats::setNames(syms, vapply(syms, .str, ""))
  lapply(sol, function(d) lapply(names(d), function(nm) list(bynm[[nm]], as_sym(d[[nm]]))))
}

# heurisch solves for the coefficients over the field of the numbers in the
# integrand (construct_domain(non_syms, field=True)), whose elements print
# as a polynomial in the generator; a constant value from the rref over
# expressions is brought to that form by rationalizing the denominator
.fs_heur_canon <- function(v) {
  if (v$type %in% c("Integer", "Rational") || length(free_symbols(v)) || !.has_type(v, "Pow")) return(v)
  r <- tryCatch(expand(radsimp(v)), error = function(e) v)
  if (length(free_symbols(r))) v else r
}
