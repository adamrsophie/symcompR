# Galois groups of polynomials of degree at most 6, after
# sympy/polys/numberfields/galoisgroups.py, galois_resolvents.py and
# sympy/combinatorics/galois.py. The resolvent coefficient table of
# resolvent_lookup.py lives in galois-resolvent-lookup.R.
#
# Polynomials over ZZ are handled as dense bigq vectors with integer values,
# leading coefficient first (the up_* helpers of polys-core.R).

# ---------------------------------------------------------------------------
# Exceptions: GaloisGroupException and its subclasses

.gg_error <- function(name, message = "") {
  parents <- if (name == "GaloisGroupException") character() else "GaloisGroupException"
  structure(list(message = message, call = NULL),
            class = c(name, parents, "sympy_error", "error", "condition"))
}

# ---------------------------------------------------------------------------
# Transitive subgroups of S1, ..., S6 (sympy/combinatorics/galois.py)

.gg_enum <- function(cls, names, values = names) {
  stats::setNames(Map(function(n, v) structure(list(cls = cls, name = n, value = v), class = "TransitiveSubgroupName"),
                      names, values), names)
}

#' Names of the transitive subgroups of S1 to S6
#'
#' Each of these lists holds the names of the transitive subgroups of the
#' symmetric group of the given degree, up to conjugacy, as in Cohen's
#' *A Course in Computational Algebraic Number Theory* and SymPy's
#' `sympy.combinatorics.galois`. A name prints like SymPy's enum member
#' (`S4TransitiveSubgroups.V`); [get_perm_group()] returns the group itself.
#' [galois_group()] returns these names when `by_name = TRUE`.
#'
#' @format Named lists of `TransitiveSubgroupName` objects.
#' @examples
#' S4TransitiveSubgroups$V
#' names(S6TransitiveSubgroups)
#' @name TransitiveSubgroups
#' @export
S1TransitiveSubgroups <- .gg_enum("S1TransitiveSubgroups", "S1")

#' @rdname TransitiveSubgroups
#' @export
S2TransitiveSubgroups <- .gg_enum("S2TransitiveSubgroups", "S2")

#' @rdname TransitiveSubgroups
#' @export
S3TransitiveSubgroups <- .gg_enum("S3TransitiveSubgroups", c("A3", "S3"))

#' @rdname TransitiveSubgroups
#' @export
S4TransitiveSubgroups <- .gg_enum("S4TransitiveSubgroups", c("C4", "V", "D4", "A4", "S4"))

#' @rdname TransitiveSubgroups
#' @export
S5TransitiveSubgroups <- .gg_enum("S5TransitiveSubgroups", c("C5", "D5", "M20", "A5", "S5"))

#' @rdname TransitiveSubgroups
#' @export
S6TransitiveSubgroups <- .gg_enum(
  "S6TransitiveSubgroups",
  c("C6", "S3", "D6", "A4", "G18", "A4xC2", "S4m", "S4p", "G36m", "G36p", "S4xC2", "PSL2F5", "G72",
    "PGL2F5", "A6", "S6"),
  c("C6", "S3", "D6", "A4", "G18", "A4 x C2", "S4-", "S4+", "G36-", "G36+", "S4 x C2", "PSL2(F5)", "G72",
    "PGL2(F5)", "A6", "S6"))

#' @export
format.TransitiveSubgroupName <- function(x, ...) paste0(x$cls, ".", x$name)

#' @export
print.TransitiveSubgroupName <- function(x, ...) {
  cat(format(x), "\n", sep = "")
  invisible(x)
}

# repr() of the enum member, as it appears inside a tuple
.gg_name_repr <- function(x) paste0("<", x$cls, ".", x$name, ": '", x$value, "'>")

#' Permutation group of a named transitive subgroup
#'
#' `get_perm_group(name)` returns the permutation group for a name from
#' [S1TransitiveSubgroups] to [S6TransitiveSubgroups], as the method of the
#' same name in SymPy. The constructors `four_group()`, `M20()`,
#' `S3_in_S6()`, `A4_in_S6()`, `S4m()`, `S4p()`, `A4xC2()`, `S4xC2()`,
#' `G18()`, `G36m()`, `G36p()`, `G72()`, `PSL2F5()` and `PGL2F5()` build the
#' groups that are not symmetric, alternating, cyclic or dihedral, with
#' SymPy's generators.
#'
#' @param name A `TransitiveSubgroupName`.
#' @return A `PermutationGroup`.
#' @examples
#' G <- get_perm_group(S4TransitiveSubgroups$V)
#' group_order(G)
#' group_order(M20())
#' @export
get_perm_group <- function(name) {
  switch(format(name),
         S1TransitiveSubgroups.S1 = SymmetricGroup(1),
         S2TransitiveSubgroups.S2 = SymmetricGroup(2),
         S3TransitiveSubgroups.A3 = AlternatingGroup(3),
         S3TransitiveSubgroups.S3 = SymmetricGroup(3),
         S4TransitiveSubgroups.C4 = CyclicGroup(4),
         S4TransitiveSubgroups.V = four_group(),
         S4TransitiveSubgroups.D4 = DihedralGroup(4),
         S4TransitiveSubgroups.A4 = AlternatingGroup(4),
         S4TransitiveSubgroups.S4 = SymmetricGroup(4),
         S5TransitiveSubgroups.C5 = CyclicGroup(5),
         S5TransitiveSubgroups.D5 = DihedralGroup(5),
         S5TransitiveSubgroups.M20 = M20(),
         S5TransitiveSubgroups.A5 = AlternatingGroup(5),
         S5TransitiveSubgroups.S5 = SymmetricGroup(5),
         S6TransitiveSubgroups.C6 = CyclicGroup(6),
         S6TransitiveSubgroups.S3 = S3_in_S6(),
         S6TransitiveSubgroups.D6 = DihedralGroup(6),
         S6TransitiveSubgroups.A4 = A4_in_S6(),
         S6TransitiveSubgroups.G18 = G18(),
         S6TransitiveSubgroups.A4xC2 = A4xC2(),
         S6TransitiveSubgroups.S4m = S4m(),
         S6TransitiveSubgroups.S4p = S4p(),
         S6TransitiveSubgroups.G36m = G36m(),
         S6TransitiveSubgroups.G36p = G36p(),
         S6TransitiveSubgroups.S4xC2 = S4xC2(),
         S6TransitiveSubgroups.PSL2F5 = PSL2F5(),
         S6TransitiveSubgroups.G72 = G72(),
         S6TransitiveSubgroups.PGL2F5 = PGL2F5(),
         S6TransitiveSubgroups.A6 = AlternatingGroup(6),
         S6TransitiveSubgroups.S6 = SymmetricGroup(6))
}

# Permutation(a, b, ...)(c, d, ...) of SymPy, padded to `size`
.gg_perm <- function(..., size = NULL) Permutation(list(...), size = size)

#' @rdname get_perm_group
#' @export
four_group <- function() PermutationGroup(list(.gg_perm(c(0, 1), c(2, 3)), .gg_perm(c(0, 2), c(1, 3))))

#' @rdname get_perm_group
#' @export
M20 <- function() {
  G <- PermutationGroup(list(.gg_perm(c(0, 1, 2, 3, 4)), .gg_perm(c(1, 2, 4, 3))))
  G$degree <- 5L
  G$order <- 20
  G$is_transitive <- TRUE
  G$is_sym <- FALSE
  G$is_alt <- FALSE
  G$is_cyclic <- FALSE
  G$is_dihedral <- FALSE
  G
}

#' @rdname get_perm_group
#' @export
S3_in_S6 <- function() {
  # the symmetries of a triangular prism
  G <- PermutationGroup(list(.gg_perm(c(0, 1, 2), c(3, 4, 5)), .gg_perm(c(0, 3), c(2, 4), c(1, 5))))
  .set_symmetric_group_properties(G, 3, 6)
  G
}

#' @rdname get_perm_group
#' @export
A4_in_S6 <- function() {
  G <- PermutationGroup(list(.gg_perm(c(0, 4, 5), c(1, 3, 2)), .gg_perm(c(0, 1, 2), c(3, 5, 4))))
  .set_alternating_group_properties(G, 4, 6)
  G
}

#' @rdname get_perm_group
#' @export
S4m <- function() {
  G <- PermutationGroup(list(.gg_perm(c(1, 4, 5, 3)), .gg_perm(c(0, 4), c(1, 5), c(2, 3))))
  .set_symmetric_group_properties(G, 4, 6)
  G
}

#' @rdname get_perm_group
#' @export
S4p <- function() {
  G <- PermutationGroup(list(.gg_perm(c(0, 2, 4, 1), c(3, 5)), .gg_perm(c(0, 3), c(4, 5))))
  .set_symmetric_group_properties(G, 4, 6)
  G
}

#' @rdname get_perm_group
#' @export
A4xC2 <- function() {
  PermutationGroup(list(.gg_perm(c(0, 4, 5), c(1, 3, 2)), .gg_perm(c(0, 1, 2), c(3, 5, 4)),
                        .gg_perm(c(2, 4), size = 6)))
}

#' @rdname get_perm_group
#' @export
S4xC2 <- function() {
  PermutationGroup(list(.gg_perm(c(1, 4, 5, 3)), .gg_perm(c(0, 4), c(1, 5), c(2, 3)), .gg_perm(c(1, 4), c(3, 5))))
}

#' @rdname get_perm_group
#' @export
G18 <- function() {
  PermutationGroup(list(.gg_perm(c(0, 1, 2), size = 6), .gg_perm(c(3, 4, 5)), .gg_perm(c(0, 4), c(1, 5), c(2, 3))))
}

#' @rdname get_perm_group
#' @export
G36m <- function() {
  PermutationGroup(list(.gg_perm(c(0, 1, 2), size = 6), .gg_perm(c(3, 4, 5)), .gg_perm(c(1, 2), c(3, 5)),
                        .gg_perm(c(0, 4), c(1, 5), c(2, 3))))
}

#' @rdname get_perm_group
#' @export
G36p <- function() {
  PermutationGroup(list(.gg_perm(c(0, 1, 2), size = 6), .gg_perm(c(3, 4, 5)), .gg_perm(c(0, 5, 2, 3), c(1, 4))))
}

#' @rdname get_perm_group
#' @export
G72 <- function() {
  PermutationGroup(list(.gg_perm(c(0, 1, 2), size = 6), .gg_perm(c(0, 4, 1, 3), c(2, 5)),
                        .gg_perm(c(0, 3), c(1, 4), c(2, 5))))
}

#' @rdname get_perm_group
#' @export
PSL2F5 <- function() {
  G <- PermutationGroup(list(.gg_perm(c(0, 4, 5), c(1, 3, 2)), .gg_perm(c(0, 4, 3, 1, 5))))
  .set_alternating_group_properties(G, 5, 6)
  G
}

#' @rdname get_perm_group
#' @export
PGL2F5 <- function() {
  G <- PermutationGroup(list(.gg_perm(c(0, 1, 2, 3, 4)), .gg_perm(c(0, 5), c(1, 2), c(3, 4))))
  .set_symmetric_group_properties(G, 5, 6)
  G
}

# ---------------------------------------------------------------------------
# Dense univariate polynomials over ZZ (bigq vectors with integer entries)

.gg_key <- function(f) paste(as.character(f), collapse = ",")

# dup_resultant over QQ by the Euclidean remainder sequence:
# Res(f, g) = (-1)^(deg f deg g) lc(g)^(deg f - deg r) Res(g, r), r = f mod g
.gg_resultant <- function(f, g) {
  f <- up_trim(f)
  g <- up_trim(g)
  if (!length(f) || !length(g)) return(gmp::as.bigq(0))
  res <- gmp::as.bigq(1)
  repeat {
    n <- up_deg(f)
    m <- up_deg(g)
    if (m == 0L) return(res * g[1L]^n)
    r <- up_divmod(f, g)$r
    if (!length(r)) return(gmp::as.bigq(0))
    if ((n * m) %% 2L == 1L) res <- -res
    res <- res * g[1L]^(n - up_deg(r))
    f <- g
    g <- r
  }
}

# dup_discriminant
.gg_discriminant <- function(f) {
  d <- up_deg(f)
  if (d <= 0L) return(gmp::as.bigq(0))
  s <- if (((d * (d - 1L)) %/% 2L) %% 2L == 1L) -1L else 1L
  .gg_resultant(f, up_diff(f)) / (f[1L] * s)
}

.gg_is_square <- function(d) .is_square_int(gmp::numerator(gmp::as.bigq(d)))

# has_square_disc
.gg_has_square_disc <- function(f) .gg_is_square(.gg_discriminant(.gg_coeffs(f)))

# dup_sqf_p: f has no repeated factor (the zero polynomial counts as square free)
.gg_sqf_p <- function(f) length(up_gcd(f, up_diff(f))) <= 1L

# dup_factor_list over ZZ: list of list(factor, multiplicity)
.gg_factor_list <- function(f) {
  fl <- qp_factor_list(up_to_qpoly(f))
  lapply(fl$factors, function(pr) list(up_from_qpoly(pr[[1L]], 1L), pr[[2L]]))
}

# dup_irreducible_p: no factor, or one factor of multiplicity one
.gg_irreducible_p <- function(f) {
  fl <- .gg_factor_list(f)
  !length(fl) || (length(fl) == 1L && fl[[1L]][[2L]] == 1L)
}

# coefficient vector of a Poly, an expression, or a vector
.gg_coeffs <- function(T) {
  if (inherits(T, "Poly")) return(up_from_qpoly(T$rep, 1L))
  if (inherits(T, "bigq")) return(T)
  if (is.numeric(T) || inherits(T, "bigz")) return(gmp::as.bigq(T))
  up_from_qpoly(Poly(T)$rep, 1L)
}

# ---------------------------------------------------------------------------
# Resolvent lookup (get_resolvent_by_lookup)

.gg_eval_terms <- function(tm, s) {
  n <- length(s)
  E <- matrix(as.integer(tm$e), ncol = n, byrow = TRUE)
  v <- gmp::as.bigz(tm$c)
  for (i in seq_len(n)) if (any(E[, i] > 0L)) v <- v * gmp::pow.bigz(s[i], E[, i])
  sum(v)
}

# The resolvent `number` of the monic integer polynomial T, as a dense
# coefficient vector. As in SymPy, the table's arguments s1, ..., sn receive
# the coefficients of T after the leading one.
.gg_resolvent_by_lookup <- function(T, number) {
  Tc <- .gg_coeffs(T)
  L <- .GG_RESOLVENT_COEFF_TERMS[[paste0(up_deg(Tc), "_", number)]]
  s <- gmp::numerator(Tc[-1L])
  gmp::as.bigq(c(gmp::as.bigz(1), do.call(c, lapply(L, .gg_eval_terms, s = s))))
}

# ---------------------------------------------------------------------------
# Tschirnhausen transformations

# coeff_search(m, R) of numberfields/utilities.py, as a stateful generator:
# each call returns the next coefficient vector.
.gg_coeff_search <- function(m, R) {
  R0 <- R
  cf <- rep(R, m)
  started <- FALSE
  function() {
    repeat {
      if (started) {
        j <- m
        while (cf[j] == -R) j <- j - 1L
        cf[j] <<- cf[j] - 1L
        if (j < m) cf[(j + 1L):m] <<- R
        if (all(cf == 0L)) {
          R <<- R + 1L
          cf <<- rep(R, m)
        }
      }
      started <<- TRUE
      if (R == R0 || R %in% cf || -R %in% cf) return(cf)
    }
  }
}

# history of tschirnhausen_transformation: an environment used as a set of
# coefficient vectors
.gg_hist_add <- function(h, f) assign(.gg_key(f), TRUE, envir = h)
.gg_hist_has <- function(h, f) exists(.gg_key(f), envir = h, inherits = FALSE)

# random.randint(a, b) and dup_random(n, a, b, ZZ) on Python's generator
.gg_randint <- function(a, b) as.numeric(.py_randrange(a, b + 1))
.gg_dup_random <- function(n, a, b) {
  f <- vapply(seq_len(n + 1L), function(i) .gg_randint(a, b), 0)
  while (f[1L] == 0) f[1L] <- .gg_randint(a, b)
  gmp::as.bigq(f)
}

# U = Res_x(T(x), X - A(x)) = prod (X - A(r)) over the roots r of the monic T:
# the characteristic polynomial of A(theta), from the power sums
# Tr(A(theta)^k) and Newton's identities.
.gg_newton_coeffs <- function(p) {
  n <- length(p)
  u <- gmp::as.bigq(rep(0, n))
  for (k in seq_len(n)) {
    acc <- p[k]
    for (i in seq_len(k - 1L)) acc <- acc + u[i] * p[k - i]
    u[k] <- -acc / k
  }
  c(gmp::as.bigq(1), u)
}

# power sums P_k = sum r^k over the roots r of T, k = 0..K (Newton's identities)
.gg_power_sums <- function(Tc, K) {
  n <- up_deg(Tc)
  cT <- Tc[-1L] / Tc[1L]
  P <- gmp::as.bigq(rep(0, K + 1L))
  P[1L] <- n
  for (k in seq_len(K)) {
    acc <- if (k <= n) k * cT[k] else gmp::as.bigq(0)
    for (i in seq_len(min(k - 1L, n))) acc <- acc + cT[i] * P[k - i + 1L]
    P[k + 1L] <- -acc
  }
  P
}

.gg_tsch_U <- function(Tc, a) {
  n <- up_deg(Tc)
  P <- .gg_power_sums(Tc, n - 1L)
  trace <- function(h) {
    h <- up_trim(h)
    if (!length(h)) return(gmp::as.bigq(0))
    d <- up_deg(h)
    sum(h * P[(d:0) + 1L])
  }
  A <- up_divmod(a, Tc)$r
  B <- gmp::as.bigq(1)
  p <- gmp::as.bigq(rep(0, n))
  for (k in seq_len(n)) {
    B <- up_divmod(up_mul(B, A), Tc)$r
    p[k] <- trace(B)
  }
  .gg_newton_coeffs(p)
}

# tschirnhausen_transformation on coefficient vectors: list(A, U)
.gg_tsch <- function(Tc, max_coeff = 10, max_tries = 30, history = NULL, fixed_order = TRUE) {
  n <- up_deg(Tc)
  if (is.null(history)) history <- new.env()
  .gg_hist_add(history, Tc)
  if (fixed_order) {
    coeff_generators <- list()
    deg_coeff_sum <- 3L
    current_degree <- 2L
  }
  get_coeff_generator <- function(degree) {
    key <- as.character(degree)
    if (is.null(coeff_generators[[key]])) coeff_generators[[key]] <<- .gg_coeff_search(degree, 1L)
    coeff_generators[[key]]
  }
  for (i in seq_len(max_tries) - 1L) {
    # Linear A(x) are never used: a fixed linear change of the roots only
    # multiplies the discriminant by a square.
    if (fixed_order) {
      # move through the (degree, max coeff) space along lines of constant sum;
      # a draw that overshoots the sum is consumed and dropped, as in SymPy
      coeffs <- get_coeff_generator(current_degree)()
      m <- max(abs(coeffs))
      if (current_degree + m > deg_coeff_sum) {
        if (current_degree == 2L) {
          deg_coeff_sum <- deg_coeff_sum + 1L
          current_degree <- deg_coeff_sum - 1L
        } else {
          current_degree <- current_degree - 1L
        }
        coeffs <- get_coeff_generator(current_degree)()
      }
      a <- gmp::as.bigq(c(1L, coeffs))
    } else {
      # progressive coefficient bound, five tries per bound
      C <- min(i %/% 5L + 1L, max_coeff)
      d <- .gg_randint(2, n - 1L)
      a <- .gg_dup_random(d, -C, C)
    }
    U <- .gg_tsch_U(Tc, a)
    if (!.gg_hist_has(history, U) && .gg_sqf_p(U)) return(list(A = a, U = U))
  }
  stop(.gg_error("MaxTriesException"))
}

#' Tschirnhausen transformation of a polynomial
#'
#' Given a monic irreducible polynomial `T` over the integers, find another
#' such polynomial `U` defining the same number field, together with the
#' polynomial `A` that maps the roots of `T` to the roots of `U`. This is
#' Algorithm 6.3.4 of Cohen, as in SymPy's
#' `sympy.polys.numberfields.galoisgroups.tschirnhausen_transformation`.
#'
#' With `fixed_order = TRUE` the candidates `A` are tried in SymPy's fixed
#' order, from small degrees and coefficients upward; otherwise they are drawn
#' at random from the Python compatible generator.
#'
#' @param T A monic `Poly` over the integers.
#' @param max_coeff Bound on the coefficients of random candidates.
#' @param max_tries Number of candidates to try.
#' @param history `NULL` or an environment used as a set of coefficient
#'   vectors that `U` must avoid; `T` is added to it.
#' @param fixed_order Use the deterministic order.
#' @return A list `(A, U)` of `Poly` objects; `U` is in a dummy variable `X`.
#' @examples
#' x <- Symbol("x")
#' tschirnhausen_transformation(Poly(x^2 - 2, x))
#' @export
tschirnhausen_transformation <- function(T, max_coeff = 10, max_tries = 30, history = NULL, fixed_order = TRUE) {
  x <- T$gens[[1L]]
  r <- .gg_tsch(.gg_coeffs(T), max_coeff, max_tries, history, fixed_order)
  X <- Dummy("X")
  list(Poly(up_to_expr(r$A, x), x), Poly(up_to_expr(r$U, X), X))
}

# ---------------------------------------------------------------------------
# Resolvents computed from approximate roots (class Resolvent)

.gg_cx <- function(re, im = gmp::as.bigq(0)) list(re = gmp::as.bigq(re), im = gmp::as.bigq(im))
.gg_cx_add <- function(a, b) list(re = a$re + b$re, im = a$im + b$im)
.gg_cx_sub <- function(a, b) list(re = a$re - b$re, im = a$im - b$im)
.gg_cx_mul <- function(a, b) list(re = a$re * b$re - a$im * b$im, im = a$re * b$im + a$im * b$re)
.gg_cx_div <- function(a, b) {
  d <- b$re^2 + b$im^2
  list(re = (a$re * b$re + a$im * b$im) / d, im = (a$im * b$re - a$re * b$im) / d)
}
.gg_cx_pow <- function(a, k) {
  r <- .gg_cx(1)
  for (i in seq_len(k)) r <- .gg_cx_mul(r, a)
  r
}
.gg_cx_horner <- function(f, z) {
  r <- .gg_cx(0)
  for (c in as.list(f)) r <- .gg_cx_add(.gg_cx_mul(r, z), .gg_cx(c))
  r
}

# nearest integer of a rational (round_mpf)
.gg_nint <- function(q) {
  q <- q + gmp::as.bigq(1, 2)
  gmp::numerator(q) %/% gmp::denominator(q)
}

# q rounded to a multiple of 2^-k
.gg_dyadic <- function(q, k) {
  s <- gmp::pow.bigz(2, k)
  gmp::as.bigq(.gg_nint(q * s), s)
}

#' @noRd
.gg_Resolvent <- function(F, X, s) {
  rep <- do.call(Poly, c(list(F), X))$rep
  m <- length(s)
  t <- NULL
  r <- 0L
  for (k in seq_along(rep$c)) {
    if (abs(rep$c[k]) != 1) stop(.gg_error("ResolventException", "Resolvent class expects forms with unit coeffs"))
    tk <- sum(rep$e[k, ])
    if (!is.null(t) && tk != t) stop(.gg_error("ResolventException", "Resolvent class expects homogeneous forms"))
    t <- tk
    r <- r + 1L
  }
  if (!(m < r * 2^t)) stop(.gg_error("ResolventException", "Resolvent class expects m < r*2^t"))
  R <- list(F = F, X = X, s = s, m = m, t = t, r = r, e = rep$e, c = rep$c)
  # The conjugate (sigma F) is F evaluated at X_i = (argument sigma(i)).
  R$root_lambdas <- lapply(s, function(sigma) {
    force(sigma)
    function(...) .gg_form_eval(R, sigma, list(...))
  })
  structure(R, class = "Resolvent")
}

# F at X_i = vals[sigma(i)], symbolically or in complex rational arithmetic
.gg_form_eval <- function(R, sigma, vals) {
  idx <- sigma$array + 1L
  if (is_sym(vals[[1L]])) {
    terms <- lapply(seq_along(R$c), function(k) {
      fs <- lapply(which(R$e[k, ] > 0L), function(i) .pow(vals[[idx[i]]], .int(R$e[k, i])))
      .mul(c(list(make_rational(R$c[k])), fs))
    })
    return(.add(terms))
  }
  out <- .gg_cx(0)
  for (k in seq_along(R$c)) {
    term <- .gg_cx(R$c[k])
    for (i in which(R$e[k, ] > 0L)) term <- .gg_cx_mul(term, .gg_cx_pow(vals[[idx[i]]], R$e[k, i]))
    out <- .gg_cx_add(out, term)
  }
  out
}

# get_prec: bits of accuracy needed in the roots of T, for roots bounded by M
.gg_res_get_prec <- function(R, M, target = "coeffs") {
  M <- max(M, 2)
  m <- R$m; t <- R$t; r <- R$r
  lg <- if (target == "coeffs") m * log2(r) + log2(m * t + 1) + (m * t - 1) * log2(M)
        else log2(r) + log2(t + 1) + (t - 1) * log2(M)
  as.integer(ceiling(1 + lg)) + 1L
}

# approximate_roots_of_poly: the roots of T in CRootOf order (real roots
# increasing, then the roots of positive imaginary part by real and imaginary
# part, each preceded by its conjugate), refined by Newton's method in exact
# arithmetic. Only the roots in the upper half plane are refined; their
# conjugates are taken exactly, so a conjugate pair is never split.
# ponytail: double precision starting points from polyroot() and the order
# taken from the root values, not from SymPy's isolating rectangles; a
# refinement that does not converge raises ResolventException rather than
# rounding.
.gg_res_approx_roots <- function(R, Tc, target = "coeffs") {
  x <- Symbol("x")
  v <- vapply(.numeric_roots(up_to_expr(Tc, x), x), as.complex, 0i)
  reals <- sort(Re(v[Im(v) == 0]))
  upper <- v[Im(v) > 0]
  upper <- upper[order(Re(upper), Im(upper))]
  if (length(reals) + 2L * length(upper) != up_deg(Tc)) {
    stop(.gg_error("ResolventException", "could not separate the real and complex roots"))
  }
  M <- max(Mod(v)) + 1
  m <- .gg_res_get_prec(R, M, target)
  bits <- m + as.integer(ceiling(log2(M))) + 9L
  dT <- up_diff(Tc)
  tol2 <- gmp::as.bigq(1, gmp::pow.bigz(2, 2L * bits))
  refine <- function(z0) {
    z <- .gg_cx(Re(z0), Im(z0))
    for (it in seq_len(60L)) {
      st <- .gg_cx_div(.gg_cx_horner(Tc, z), .gg_cx_horner(dT, z))
      z <- .gg_cx_sub(z, st)
      z <- list(re = .gg_dyadic(z$re, bits + 8L), im = .gg_dyadic(z$im, bits + 8L))
      if (st$re^2 + st$im^2 < tol2) return(z)
    }
    stop(.gg_error("ResolventException", "root approximation did not converge"))
  }
  out <- lapply(reals, function(r) refine(complex(real = r)))
  for (u in upper) {
    z <- refine(u)
    out <- c(out, list(list(re = z$re, im = -z$im), z))
  }
  out
}

# round_roots_to_integers_for_poly: named list, names are the (0 based)
# indices of the conjugates whose value rounds to a rational integer
.gg_res_round_roots <- function(R, T) {
  rt <- .gg_res_approx_roots(R, .gg_coeffs(T), "roots")
  out <- list()
  for (j in seq_len(R$m)) {
    v <- .gg_form_eval(R, R$s[[j]], rt)
    if (.gg_nint(v$im) == 0) out[[as.character(j - 1L)]] <- .gg_nint(v$re)
  }
  out
}

# eval_for_poly: list(R as coefficient vector, integer root or NULL, its 0
# based index or NULL)
.gg_res_eval_for_poly <- function(R, T, find_integer_root = FALSE) {
  rt <- .gg_res_approx_roots(R, .gg_coeffs(T), "coeffs")
  vals <- lapply(R$s, function(sigma) .gg_form_eval(R, sigma, rt))
  cf <- list(.gg_cx(1))
  for (v in vals) {
    nxt <- c(cf, list(.gg_cx(0)))
    for (k in seq_along(cf)) nxt[[k + 1L]] <- .gg_cx_sub(nxt[[k + 1L]], .gg_cx_mul(cf[[k]], v))
    cf <- nxt
  }
  Rd <- gmp::as.bigq(do.call(c, lapply(cf, function(c) {
    if (.gg_nint(c$im) != 0) {
      stop(.gg_error("ResolventException", paste0("Got non-integer coeff for resolvent: ",
                                                  format(as.numeric(c$re)), "+", format(as.numeric(c$im)), "j")))
    }
    .gg_nint(c$re)
  })))
  a0 <- NULL
  i0 <- NULL
  if (find_integer_root) {
    for (i in seq_along(vals)) {
      if (.gg_nint(vals[[i]]$im) != 0) next
      a <- .gg_nint(vals[[i]]$re)
      if (up_eval(Rd, a) == 0) {
        a0 <- a
        i0 <- i - 1L
        break
      }
    }
  }
  list(Rd, a0, i0)
}

# ---------------------------------------------------------------------------
# define_resolvents: the forms, variables and coset representatives behind
# the lookup table

.gg_X <- function(n) lapply(seq_len(n) - 1L, function(i) Symbol(paste0("X", i)))
.gg_idperm <- function(n) Permutation(n - 1L)
.gg_transp <- function(i, j, n) Permutation(list(c(i, j)), size = n)

# F(X_sigma(0), ..., X_sigma(n-1)), the form f() taking a list of variables
.gg_permuted <- function(f, X, sigma) f(X[sigma$array + 1L])

.gg_F51 <- function(X) {
  X[[1]]^2 * (X[[2]] * X[[5]] + X[[3]] * X[[4]]) + X[[2]]^2 * (X[[3]] * X[[1]] + X[[4]] * X[[5]]) +
    X[[3]]^2 * (X[[4]] * X[[2]] + X[[5]] * X[[1]]) + X[[4]]^2 * (X[[5]] * X[[3]] + X[[1]] * X[[2]]) +
    X[[5]]^2 * (X[[1]] * X[[4]] + X[[2]] * X[[3]])
}
.gg_s51 <- function() c(list(.gg_idperm(5)), lapply(1:4, function(j) .gg_transp(0, j, 5)), list(.gg_transp(1, 4, 5)))

.gg_define_resolvents <- function() {
  X4 <- .gg_X(4)
  X <- X4
  F40 <- expand(X[[1]] * X[[2]]^2 + X[[2]] * X[[3]]^2 + X[[3]] * X[[4]]^2 + X[[4]] * X[[1]]^2)
  s40 <- list(.gg_idperm(4), .gg_transp(0, 1, 4), .gg_transp(0, 2, 4), .gg_transp(0, 3, 4),
              .gg_transp(1, 2, 4), .gg_transp(2, 3, 4))
  F41 <- X[[1]] * X[[3]] + X[[2]] * X[[4]]
  s41 <- list(.gg_idperm(4), .gg_transp(0, 1, 4), .gg_transp(0, 3, 4))
  X5 <- .gg_X(5)
  F51 <- expand(.gg_F51(X5))
  s51 <- .gg_s51()
  X6 <- .gg_X(6)
  term0 <- function(X) X[[1]]^2 * X[[6]]^2 * (X[[2]] * X[[5]] + X[[3]] * X[[4]])
  terms <- list()
  for (sigma in group_elements(PGL2F5())) {
    tm <- expand(.gg_permuted(term0, X6, sigma))
    terms[[tm$key]] <- tm
  }
  F61 <- .add(unname(terms))
  s61 <- c(list(.gg_idperm(6)), lapply(1:5, function(j) .gg_transp(0, j, 6)))
  F62 <- X6[[1]] * X6[[2]] * X6[[3]] + X6[[4]] * X6[[5]] * X6[[6]]
  s62 <- c(list(.gg_idperm(6)), unlist(lapply(0:2, function(i) lapply(0:2, function(j) .gg_transp(i, j + 3L, 6))),
                                      recursive = FALSE))
  list("4_0" = list(F40, X4, s40), "4_1" = list(F41, X4, s41), "5_1" = list(F51, X5, s51),
       "6_1" = list(F61, X6, s61), "6_2" = list(F62, X6, s62))
}

# ---------------------------------------------------------------------------
# Galois groups by degree. Each returns list(name, is_alternating).

.galois_group_degree_3 <- function(T, max_tries = 30, randomize = FALSE) {
  S3 <- S3TransitiveSubgroups
  if (.gg_has_square_disc(.gg_coeffs(T))) list(S3$A3, TRUE) else list(S3$S3, FALSE)
}

# Algorithm 6.3.7 of Cohen, by root approximation only
.galois_group_degree_4_root_approx <- function(T, max_tries = 30, randomize = FALSE) {
  S4 <- S4TransitiveSubgroups
  Tc <- .gg_coeffs(T)
  X <- .gg_X(4)
  # first resolvent: F = X0*X2 + X1*X3, G = S4, stabilizer D4
  s1 <- list(.gg_idperm(4), .gg_transp(0, 1, 4), .gg_transp(0, 3, 4))
  R1 <- .gg_Resolvent(X[[1]] * X[[3]] + X[[2]] * X[[4]], X, s1)
  F2_pre <- function(X) X[[1]] * X[[2]]^2 + X[[2]] * X[[3]]^2 + X[[3]] * X[[4]]^2 + X[[4]] * X[[1]]^2
  s2_pre <- list(.gg_idperm(4), .gg_transp(0, 2, 4))
  history <- new.env()
  for (i in seq_len(max_tries) - 1L) {
    if (i > 0L) Tc <- .gg_tsch(Tc, max_tries = max_tries, history = history, fixed_order = !randomize)$U
    ev <- .gg_res_eval_for_poly(R1, Tc, find_integer_root = TRUE)
    if (!.gg_sqf_p(ev[[1L]])) next
    # Gal(T) is in A4 iff disc(T) is a square (Prop 6.3.1)
    sq_disc <- .gg_has_square_disc(Tc)
    i0 <- ev[[3L]]
    if (is.null(i0)) return(if (sq_disc) list(S4$A4, TRUE) else list(S4$S4, FALSE))
    if (sq_disc) return(list(S4$V, TRUE))
    # Gal(T) is C4 or D4, in the conjugate of D4 given by the integer root
    sigma <- s1[[i0 + 1L]]
    F2 <- .gg_permuted(F2_pre, X, sigma)
    s2 <- lapply(s2_pre, function(tau) sigma * tau * sigma)
    R2 <- .gg_Resolvent(F2, X, s2)
    d <- .gg_discriminant(.gg_res_eval_for_poly(R2, Tc)[[1L]])
    if (d == 0) next
    return(if (.gg_is_square(d)) list(S4$C4, FALSE) else list(S4$D4, FALSE))
  }
  stop(.gg_error("MaxTriesException"))
}

# Algorithm 6.3.6 of Cohen, with resolvent coefficient lookup
.galois_group_degree_4_lookup <- function(T, max_tries = 30, randomize = FALSE) {
  S4 <- S4TransitiveSubgroups
  Tc <- .gg_coeffs(T)
  history <- new.env()
  ok <- FALSE
  for (i in seq_len(max_tries)) {
    R_dup <- .gg_resolvent_by_lookup(Tc, 0)
    if (.gg_sqf_p(R_dup)) {
      ok <- TRUE
      break
    }
    Tc <- .gg_tsch(Tc, max_tries = max_tries, history = history, fixed_order = !randomize)$U
  }
  if (!ok) stop(.gg_error("MaxTriesException"))
  # degrees of the irreducible factors of R, in increasing order
  L <- sort(unlist(lapply(.gg_factor_list(R_dup), function(pr) rep(up_deg(pr[[1L]]), pr[[2L]]))))
  if (identical(L, 6L)) return(if (.gg_has_square_disc(Tc)) list(S4$A4, TRUE) else list(S4$S4, FALSE))
  if (identical(L, c(1L, 1L, 4L))) return(list(S4$C4, FALSE))
  if (identical(L, c(2L, 2L, 2L))) return(list(S4$V, TRUE))
  stopifnot(identical(L, c(2L, 4L)))
  list(S4$D4, FALSE)
}

# Algorithm 6.3.9 of Cohen: resolvent coefficient lookup, then root
# approximation for C5 against D5
.galois_group_degree_5_hybrid <- function(T, max_tries = 30, randomize = FALSE) {
  S5 <- S5TransitiveSubgroups
  Tc <- .gg_coeffs(T)
  X <- .gg_X(5)
  s51 <- .gg_s51()
  R51 <- .gg_Resolvent(expand(.gg_F51(X)), X, s51)
  F2_pre <- function(X) {
    X[[1]] * X[[2]]^2 + X[[2]] * X[[3]]^2 + X[[3]] * X[[4]]^2 + X[[4]] * X[[5]]^2 + X[[5]] * X[[1]]^2
  }
  s2_pre <- list(.gg_idperm(5), Permutation(list(c(0, 1), c(2, 4)), size = 5))
  history <- new.env()
  reached_second_stage <- FALSE
  for (i in seq_len(max_tries) - 1L) {
    if (i > 0L) Tc <- .gg_tsch(Tc, max_tries = max_tries, history = history, fixed_order = !randomize)$U
    R51_dup <- .gg_resolvent_by_lookup(Tc, 1)
    if (!.gg_sqf_p(R51_dup)) next
    if (!reached_second_stage) {
      sq_disc <- .gg_has_square_disc(Tc)
      if (.gg_irreducible_p(R51_dup)) return(if (sq_disc) list(S5$A5, TRUE) else list(S5$S5, FALSE))
      if (!sq_disc) return(list(S5$M20, FALSE))
    }
    reached_second_stage <- TRUE
    # the conjugate of F51 that is an integer root of R51 picks the conjugate
    # of the second resolvent; if no candidate is a root, the last one stays
    rounded_roots <- .gg_res_round_roots(R51, Tc)
    permutation_index <- NULL
    for (nm in names(rounded_roots)) {
      permutation_index <- as.integer(nm)
      if (up_eval(R51_dup, rounded_roots[[nm]]) == 0) break
    }
    sigma <- s51[[permutation_index + 1L]]
    F2 <- .gg_permuted(F2_pre, X, sigma)
    s2 <- lapply(s2_pre, function(tau) sigma * tau * sigma)
    R2 <- .gg_Resolvent(F2, X, s2)
    d <- .gg_discriminant(.gg_res_eval_for_poly(R2, Tc)[[1L]])
    if (d == 0) next
    return(if (.gg_is_square(d)) list(S5$C5, TRUE) else list(S5$D5, TRUE))
  }
  stop(.gg_error("MaxTriesException"))
}

# Algorithm 6.3.9 of Cohen, with resolvent coefficient lookup and
# factorization over the number field of T
.galois_group_degree_5_lookup_ext_factor <- function(T, max_tries = 30, randomize = FALSE) {
  S5 <- S5TransitiveSubgroups
  Tc <- .gg_coeffs(T)
  T0 <- Tc
  history <- new.env()
  ok <- FALSE
  for (i in seq_len(max_tries)) {
    R_dup <- .gg_resolvent_by_lookup(Tc, 1)
    if (.gg_sqf_p(R_dup)) {
      ok <- TRUE
      break
    }
    Tc <- .gg_tsch(Tc, max_tries = max_tries, history = history, fixed_order = !randomize)$U
  }
  if (!ok) stop(.gg_error("MaxTriesException"))
  sq_disc <- .gg_has_square_disc(Tc)
  if (.gg_irreducible_p(R_dup)) return(if (sq_disc) list(S5$A5, TRUE) else list(S5$S5, FALSE))
  if (!sq_disc) return(list(S5$M20, FALSE))
  # Gal(T) is D5 or C5; it has order 5 iff the original T splits completely
  # over the field generated by one of its roots.
  if (.gg_ext_factor_count(T0) == 5L) list(S5$C5, TRUE) else list(S5$D5, TRUE)
}

# Number of irreducible factors of the monic, square free T over Q(theta),
# theta a root of T: dup_ext_factor (Trager) finds the first shift s = 0, 1,
# ... whose norm N_s(x) = Res_y(T(y), T(x - s*y)) = prod T(x - s*theta_i) is
# square free; the factors of T over Q(theta) then correspond one to one to
# the irreducible factors of N_s over QQ. N_s is built from its power sums
# sum_{i,j} (theta_j + s*theta_i)^k, so no symbolic resultant is needed.
.gg_ext_factor_count <- function(Tc) {
  n <- up_deg(Tc)
  N <- n * n
  P <- .gg_power_sums(Tc, N)
  s <- 0L
  repeat {
    p <- gmp::as.bigq(rep(0, N))
    for (k in seq_len(N)) {
      l <- 0:k
      p[k] <- sum(gmp::chooseZ(k, l) * gmp::pow.bigz(s, l) * P[l + 1L] * P[k - l + 1L])
    }
    r <- .gg_newton_coeffs(p)
    if (.gg_sqf_p(r)) break
    s <- s + 1L
  }
  length(.gg_factor_list(r))
}

# Algorithm 6.3.10 of Cohen, with resolvent coefficient lookup
.galois_group_degree_6_lookup <- function(T, max_tries = 30, randomize = FALSE) {
  S6 <- S6TransitiveSubgroups
  Tc <- .gg_coeffs(T)
  # first resolvent
  history <- new.env()
  ok <- FALSE
  for (i in seq_len(max_tries)) {
    R_dup <- .gg_resolvent_by_lookup(Tc, 1)
    if (.gg_sqf_p(R_dup)) {
      ok <- TRUE
      break
    }
    Tc <- .gg_tsch(Tc, max_tries = max_tries, history = history, fixed_order = !randomize)$U
  }
  if (!ok) stop(.gg_error("MaxTriesException"))
  # the factors grouped by degree
  factors_by_deg <- list()
  for (pr in .gg_factor_list(R_dup)) {
    key <- as.character(up_deg(pr[[1L]]))
    factors_by_deg[[key]] <- c(factors_by_deg[[key]], list(pr[[1L]]))
  }
  L <- sort(unlist(lapply(names(factors_by_deg), function(d) rep(as.integer(d), length(factors_by_deg[[d]])))))
  T_has_sq_disc <- .gg_has_square_disc(Tc)
  if (identical(L, c(1L, 2L, 3L))) {
    f1 <- factors_by_deg[["3"]][[1L]]
    return(if (.gg_has_square_disc(f1)) list(S6$C6, FALSE) else list(S6$D6, FALSE))
  }
  if (identical(L, c(3L, 3L))) {
    f <- factors_by_deg[["3"]]
    any_square <- .gg_has_square_disc(f[[1L]]) || .gg_has_square_disc(f[[2L]])
    return(if (any_square) list(S6$G18, FALSE) else list(S6$G36m, FALSE))
  }
  if (identical(L, c(2L, 4L))) {
    if (T_has_sq_disc) return(list(S6$S4p, TRUE))
    f1 <- factors_by_deg[["4"]][[1L]]
    return(if (.gg_has_square_disc(f1)) list(S6$A4xC2, FALSE) else list(S6$S4xC2, FALSE))
  }
  if (identical(L, c(1L, 1L, 4L))) return(if (T_has_sq_disc) list(S6$A4, TRUE) else list(S6$S4m, FALSE))
  if (identical(L, c(1L, 5L))) return(if (T_has_sq_disc) list(S6$PSL2F5, TRUE) else list(S6$PGL2F5, FALSE))
  if (identical(L, c(1L, 1L, 1L, 3L))) return(list(S6$S3, FALSE))
  stopifnot(identical(L, 6L))
  # second resolvent
  history <- new.env()
  ok <- FALSE
  for (i in seq_len(max_tries)) {
    R_dup <- .gg_resolvent_by_lookup(Tc, 2)
    if (.gg_sqf_p(R_dup)) {
      ok <- TRUE
      break
    }
    Tc <- .gg_tsch(Tc, max_tries = max_tries, history = history, fixed_order = !randomize)$U
  }
  if (!ok) stop(.gg_error("MaxTriesException"))
  T_has_sq_disc <- .gg_has_square_disc(Tc)
  if (.gg_irreducible_p(R_dup)) {
    if (T_has_sq_disc) list(S6$A6, TRUE) else list(S6$S6, FALSE)
  } else {
    if (T_has_sq_disc) list(S6$G36p, TRUE) else list(S6$G72, FALSE)
  }
}

# Poly.make_monic_over_integers_by_scaling_roots on a coefficient vector:
# list(g, c) with g(x) = c^n f(x/c) / lc(f) monic over ZZ
.gg_make_monic_over_integers <- function(f) {
  if (f[1L] == 1 && all(gmp::denominator(f) == 1)) return(list(f, gmp::as.bigz(1)))
  fm <- f / f[1L]
  c <- Reduce(gmp::lcm.bigz, as.list(gmp::denominator(fm)))
  list(fm * gmp::as.bigq(gmp::pow.bigz(c, seq_along(fm) - 1L)), c)
}

#' Galois group of a polynomial
#'
#' Computes the Galois group of an irreducible polynomial over the integers
#' or rationals of degree at most 6, as SymPy's `galois_group()`, with the
#' algorithms of Cohen's *A Course in Computational Algebraic Number Theory*
#' (Section 6.3) and SymPy's lookup tables of resolvent coefficients. When a
#' resolvent has a repeated root the polynomial is replaced by a
#' Tschirnhausen transformation, tried in SymPy's fixed order unless
#' `randomize = TRUE`.
#'
#' @param f An irreducible polynomial expression or `Poly`.
#' @param gens Optional generator (or list of generators) for converting `f`.
#' @param by_name Return the group by name, an element of
#'   [S1TransitiveSubgroups] to [S6TransitiveSubgroups], instead of a
#'   `PermutationGroup`.
#' @param max_tries Attempts allowed in the steps that generate
#'   Tschirnhausen transformations.
#' @param randomize Use random coefficients (from the Python compatible
#'   generator) for the Tschirnhausen transformations.
#' @return A `galois_group_result`: a list `(G, alt)` printed like SymPy's
#'   tuple, where `G` is the group or its name and `alt` says whether the
#'   group is contained in the alternating group of the degree of `f`.
#' @examples
#' x <- Symbol("x")
#' galois_group(x^4 + 1)
#' galois_group(x^4 - 2, by_name = TRUE)
#' group_order(galois_group(x^3 + 2)[[1]])
#' @export
galois_group <- function(f, gens = NULL, by_name = FALSE, max_tries = 30, randomize = FALSE) {
  # poly_from_expr. `odd` are the coefficients outside QQ: SymPy makes
  # number atoms such as sqrt(2) generators when no generator is given, and
  # puts them in the domain (EX for algebraic numbers, a polynomial ring
  # otherwise, RR for floats) when the generator is given.
  bad <- "Polynomial must be irreducible and univariate over ZZ or QQ."
  if (inherits(f, "Poly")) {
    P <- f
    if (.poly_ngens(P) != 1L) stop(.sympy_condition("ValueError", bad))
    if (length(P$gens) > 1L) stop(.sympy_condition("ValueError", bad))
  } else {
    e <- expand(as_sym(f))
    g <- .gens_of(list(e), if (is_sym(gens)) list(gens) else gens)
    if (!length(g)) {
      stop(.sympy_condition("ComputationFailed", paste0("galois_group(", sstr(as_sym(f)), ") failed without generators")))
    }
    if (length(g) != 1L) stop(.sympy_condition("ValueError", bad))
    odd <- Filter(function(c) !is_exact_number(c), .coeff_list(e, g[[1L]]))
    if (length(odd)) {
      ex <- !is.null(gens) && any(vapply(odd, function(c) {
        !length(free_symbols(c)) && !.has_type(c, "Float") && isTRUE(ask_fact(c, "algebraic"))
      }, TRUE))
      if (ex) stop(.sympy_condition("DomainError", "factorization not supported over EX"))
      stop(.sympy_condition("ValueError", bad))
    }
    P <- Poly(e, g[[1L]])
  }
  fl <- qp_factor_list(P$rep)$factors
  irreducible <- !length(fl) || (length(fl) == 1L && fl[[1L]][[2L]] == 1L)
  if (!irreducible) stop(.sympy_condition("ValueError", bad))
  f_c <- up_from_qpoly(P$rep, 1L)
  n <- up_deg(f_c)
  if (n > 6L) stop(.sympy_condition("ValueError", "Only polynomials up to degree 6 are supported."))
  if (n < 1L) stop(.sympy_condition("ValueError", "Constant polynomial has no Galois group."))
  res <- if (n == 1L) {
    list(S1TransitiveSubgroups$S1, TRUE)
  } else if (n == 2L) {
    list(S2TransitiveSubgroups$S2, FALSE)
  } else {
    g <- .gg_make_monic_over_integers(f_c)[[1L]]
    fun <- switch(as.character(n),
                  "3" = .galois_group_degree_3,
                  "4" = .galois_group_degree_4_lookup,
                  "5" = .galois_group_degree_5_lookup_ext_factor,
                  "6" = .galois_group_degree_6_lookup)
    fun(g, max_tries = max_tries, randomize = randomize)
  }
  G <- if (by_name) res[[1L]] else get_perm_group(res[[1L]])
  structure(list(G, res[[2L]]), class = "galois_group_result")
}

.gg_repr <- function(v) {
  if (inherits(v, "TransitiveSubgroupName")) return(.gg_name_repr(v))
  if (is.logical(v)) return(if (isTRUE(v)) "True" else "False")
  format(v)
}

#' @export
format.galois_group_result <- function(x, ...) paste0("(", .gg_repr(x[[1L]]), ", ", .gg_repr(x[[2L]]), ")")

#' @export
print.galois_group_result <- function(x, ...) {
  cat(format(x), "\n", sep = "")
  invisible(x)
}
