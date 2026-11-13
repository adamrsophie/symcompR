# Appell sequences, after sympy/polys/appellseqs.py. Dense coefficient
# vectors (highest degree first) are gmp rationals.

.dup_integrate1 <- function(p) {
  n <- length(p)
  if (!n) return(p)
  c(p / gmp::as.bigq(n:1), gmp::as.bigq(0))
}

.as_dup_eval <- function(p, a) {
  r <- gmp::as.bigq(0)
  for (i in seq_along(p)) r <- r * a + p[i]
  r
}

.dup_sub_ground <- function(p, c) {
  p[length(p)] <- p[length(p)] - c
  p
}

.dup_bernoulli <- function(n) {
  if (n < 1L) return(gmp::as.bigq(1))
  p <- gmp::as.bigq(c(1, -1), c(1, 2))
  if (n >= 2L) for (i in 2:n) {
    p <- .dup_integrate1(p * i)
    if (i %% 2L == 0L)
      p <- .dup_sub_ground(p, .as_dup_eval(p, gmp::as.bigq(1, 2)) *
                             gmp::as.bigq(gmp::as.bigz(2)^(i - 1L), gmp::as.bigz(2)^i - 1L))
  }
  p
}

.dup_bernoulli_c <- function(n) {
  p <- gmp::as.bigq(1)
  if (n >= 1L) for (i in 1:n) {
    p <- .dup_integrate1(p * i)
    if (i %% 2L == 0L)
      p <- .dup_sub_ground(p, .as_dup_eval(p, gmp::as.bigq(1)) *
                             gmp::as.bigq(gmp::as.bigz(2)^(i - 1L) - 1L, gmp::as.bigz(2)^i - 1L))
  }
  p
}

.dup_genocchi <- function(n) {
  if (n < 1L) return(gmp::as.bigq(0))
  p <- gmp::as.bigq(-1)
  if (n >= 2L) for (i in 2:n) {
    p <- .dup_integrate1(p * i)
    if (i %% 2L == 0L) {
      v <- gmp::numerator(.as_dup_eval(p, gmp::as.bigq(1)))
      p <- .dup_sub_ground(p, gmp::as.bigq(v %/% 2L))
    }
  }
  p
}

.dup_euler <- function(n) .dup_genocchi(n + 1L) / gmp::as.bigq(-n - 1L)

.dup_andre <- function(n) {
  p <- gmp::as.bigq(1)
  if (n >= 1L) for (i in 1:n) {
    p <- .dup_integrate1(p * i)
    if (i %% 2L == 0L) p <- .dup_sub_ground(p, .as_dup_eval(p, gmp::as.bigq(1)))
  }
  p
}

# sympy.polys.polytools.named_poly
.named_poly <- function(n, f, domain, name, x, polys) {
  if (is_sym(n)) n <- num_value_dbl(n)
  if (n < 0) stop(.sympy_condition("ValueError", sprintf("Cannot generate %s of index %s", name, format(n))))
  p <- f(as.integer(n))
  while (length(p) > 1L && p[1L] == 0) p <- p[-1L]
  if (length(p) == 1L && p[1L] == 0) p <- p[0L]
  # SymPy's PurePoly writes its generator as the symbol _0
  gen <- if (is.null(x)) Symbol("_0") else as_sym(x)
  k <- length(p)
  e <- if (!k) S_Zero else .add(lapply(seq_len(k), function(i)
    .mul(list(make_rational(p[i]), .pow(gen, .int(k - i))))))
  if (polys) Poly(e, gen, domain = domain) else e
}

#' Appell sequences
#'
#' The Bernoulli, central Bernoulli, Genocchi, Euler and Andre polynomials,
#' as SymPy's `bernoulli_poly()`, `bernoulli_c_poly()`, `genocchi_poly()`,
#' `euler_poly()` and `andre_poly()`. With `x = NULL` the variable is the
#' symbol `_0`, as SymPy writes it; with `polys = TRUE` the result is a
#' [Poly()].
#'
#' @param n The index, a nonnegative integer.
#' @param x The variable, or `NULL`.
#' @param polys Return a `Poly`.
#' @return An expression or a `Poly`.
#' @examples
#' x <- Symbol("x")
#' bernoulli_poly(4, x)
#' bernoulli_c_poly(4, x)
#' genocchi_poly(5, x)
#' euler_poly(3, x)
#' andre_poly(4, x)
#' @export
bernoulli_poly <- function(n, x = NULL, polys = FALSE)
  .named_poly(n, .dup_bernoulli, "QQ", "Bernoulli polynomial", x, polys)

#' @rdname bernoulli_poly
#' @export
bernoulli_c_poly <- function(n, x = NULL, polys = FALSE)
  .named_poly(n, .dup_bernoulli_c, "QQ", "central Bernoulli polynomial", x, polys)

#' @rdname bernoulli_poly
#' @export
genocchi_poly <- function(n, x = NULL, polys = FALSE)
  .named_poly(n, .dup_genocchi, "ZZ", "Genocchi polynomial", x, polys)

#' @rdname bernoulli_poly
#' @export
euler_poly <- function(n, x = NULL, polys = FALSE)
  .named_poly(n, .dup_euler, "QQ", "Euler polynomial", x, polys)

#' @rdname bernoulli_poly
#' @export
andre_poly <- function(n, x = NULL, polys = FALSE)
  .named_poly(n, .dup_andre, "ZZ", "Andre polynomial", x, polys)

#' Steps of the extended Euclidean algorithm
#'
#' `gcdex_steps(f, g, ...)` lists the intermediate triples `(s, t, r)` with
#' `f*s + g*t = r` of the extended Euclidean algorithm, as SymPy's
#' `gcdex_steps()`; the last `r` is a greatest common divisor of `f` and `g`.
#' Multivariate polynomials need the generator, given in `...` or as `gens`.
#' Integers use the algorithm on the integers.
#'
#' @param f,g Expressions.
#' @param ... The generator.
#' @param gens The generator.
#' @return A list of lists `list(s, t, r)`.
#' @examples
#' x <- Symbol("x")
#' gcdex_steps(x^4 - 2*x^3 - 6*x^2 + 12*x + 15, x^3 + x^2 - 4*x - 4)
#' @export
gcdex_steps <- function(f, g, ..., gens = NULL) {
  gl <- c(list(...), if (is.null(gens)) list() else if (is_sym(gens)) list(gens) else as.list(gens))
  if (inherits(f, "Poly")) { if (!length(gl)) gl <- f$gens; f <- as_expr(f) }
  if (inherits(g, "Poly")) g <- as_expr(g)
  f <- as_sym(f); g <- as_sym(g)
  if (!length(gl)) gl <- Filter(function(z) !identical(z$key, S_I$key), .find_gens(list(expand(f), expand(g))))
  if (!length(gl)) {
    # domain elements: the integers or the rationals
    if (!all(vapply(list(f, g), .is_rat, TRUE)))
      stop(.sympy_condition("ComputationFailed", "gcdex_steps only works on numbers in a Euclidean domain"))
    zz <- f$type == "Integer" && g$type == "Integer"
    a <- gmp::as.bigq(f$value); b <- gmp::as.bigq(g$value)
    s1 <- gmp::as.bigq(1); s2 <- gmp::as.bigq(0); t1 <- gmp::as.bigq(0); t2 <- gmp::as.bigq(1)
    out <- list()
    while (b != 0) {
      out[[length(out) + 1L]] <- lapply(list(s2, t2, b), make_rational)
      if (zz) { q <- gmp::as.bigq(gmp::numerator(a) %/% gmp::numerator(b)); r <- a - q * b }
      else { q <- a / b; r <- gmp::as.bigq(0) }
      a <- b; b <- r
      tmp <- s1 - q * s2; s1 <- s2; s2 <- tmp
      tmp <- t1 - q * t2; t1 <- t2; t2 <- tmp
    }
    return(out)
  }
  if (length(gl) > 1L) stop(.sympy_condition("ValueError", "univariate polynomial expected"))
  x <- as_sym(gl[[1L]])
  r1 <- .ex_trim(.coeff_list(f, x)); r2 <- .ex_trim(.coeff_list(g, x))
  one <- list(S_One)
  s1 <- one; s2 <- list(); t1 <- list(); t2 <- one
  df <- length(r1) - 1L; dg <- length(r2) - 1L
  if (df < dg) {
    tmp <- s1; s1 <- t1; t1 <- tmp
    tmp <- s2; s2 <- t2; t2 <- tmp
    tmp <- r1; r1 <- r2; r2 <- tmp
  }
  out <- list()
  for (k in seq_len(max(df, dg) + 1L)) {
    out[[length(out) + 1L]] <- lapply(list(s2, t2, r2), .ex_to_expr, x = x)
    d <- .ex_divmod(r1, r2)
    if (!length(d$r)) break
    r1 <- r2; r2 <- d$r
    tmp <- .ex_sub(s1, .ex_mul(d$q, s2)); s1 <- s2; s2 <- tmp
    tmp <- .ex_sub(t1, .ex_mul(d$q, t2)); t1 <- t2; t2 <- tmp
  }
  out
}
