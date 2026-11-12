# Diophantine equations, after sympy/solvers/diophantine/diophantine.py.
# Integers are gmp bigz throughout; Python floor division is %/% and the
# Python remainder is .pmod().

.Z <- function(x) {
  if (inherits(x, "bigz")) return(x)
  if (is_sym(x)) {
    if (x$type != "Integer") stop("expected an integer, got ", .str(x), call. = FALSE)
    return(gmp::as.bigz(gmp::numerator(gmp::as.bigq(x$value))))
  }
  gmp::as.bigz(x)
}
.zs <- function(z) make_rational(gmp::as.bigq(z))
.pmod <- function(a, b) a - b * (a %/% b)
.zsign <- function(a) { a <- .Z(a); if (a > 0) 1L else if (a < 0) -1L else 0L }
.zabs <- function(a) { a <- .Z(a); if (a < 0) -a else a }
.zgcd <- function(...) {
  v <- lapply(list(...), .Z)
  g <- gmp::as.bigz(0)
  for (a in v) g <- gmp::gcd.bigz(g, a)
  g
}
.zlcm <- function(...) {
  v <- lapply(list(...), .Z)
  l <- gmp::as.bigz(1)
  for (a in v) l <- if (a == 0) gmp::as.bigz(0) else gmp::lcm.bigz(l, a)
  .zabs(l)
}

.isqrtz <- function(n) {
  n <- .Z(n)
  if (n < 0) stop("n must be nonnegative", call. = FALSE)
  if (n < gmp::as.bigz(2)^50) {
    r <- gmp::as.bigz(floor(sqrt(as.numeric(n))))
  } else {
    r <- gmp::as.bigz(2)^((gmp::sizeinbase(n, 2) + 1L) %/% 2L)
    repeat {
      r1 <- (r + n %/% r) %/% 2
      if (r1 >= r) break
      r <- r1
    }
  }
  while (r * r > n) r <- r - 1
  while ((r + 1) * (r + 1) <= n) r <- r + 1
  r
}
.is_sqz <- function(n) { n <- .Z(n); n >= 0 && .isqrtz(n)^2 == n }
# integer_nthroot(n, 2): floor root and exactness; ValueError for n < 0.
.zroot2 <- function(n) { r <- .isqrtz(n); list(r, r * r == .Z(n)) }

# SymPy gcdext (pure Python ground types): returns (g, x, y).
.gcdext <- function(a, b) {
  a <- .Z(a); b <- .Z(b)
  if (a == 0 || b == 0) {
    g <- .zabs(a) + .zabs(b)
    if (g == 0) return(list(g, g, g))
    return(list(g, a %/% g, b %/% g))
  }
  xs <- .zsign(a); a <- .zabs(a)
  ys <- .zsign(b); b <- .zabs(b)
  x <- gmp::as.bigz(1); r <- gmp::as.bigz(0)
  y <- gmp::as.bigz(0); s <- gmp::as.bigz(1)
  while (b != 0) {
    q <- a %/% b; c <- a - q * b
    a <- b; b <- c
    t <- x - q * r; x <- r; r <- t
    t <- y - q * s; y <- s; s <- t
  }
  list(a, x * xs, y * ys)
}

.remove_gcd <- function(...) {
  x <- list(...)
  if (length(x) == 1L && is.list(x[[1L]])) x <- x[[1L]]
  if (all(vapply(x, function(v) inherits(v, "bigz") || is.numeric(v) || (is_sym(v) && v$type == "Integer"), TRUE))) {
    z <- lapply(x, .Z)
    g <- do.call(.zgcd, z)
    if (g == 1 || g == 0) return(z)
    return(lapply(z, function(v) v %/% g))
  }
  # symbolic entries: divide by the gcd of the numeric contents
  ex <- lapply(x, as_sym)
  nz <- Filter(function(v) !is_zero_num(v), ex)
  if (length(nz) < 2L) return(ex)
  cs <- lapply(nz, function(v) gmp::as.bigz(gmp::numerator(.q_content(v))))
  g <- do.call(.zgcd, cs)
  if (g == 1) return(ex)
  lapply(ex, function(v) expand(.div(v, .zs(g))))
}

.rational_pq <- function(a, b) .remove_gcd(.zsign(b) * .Z(a), .zabs(b))
.nint_or_floor <- function(p, q) {
  w <- p %/% q
  r <- p - q * w
  if (.zabs(r) <= .zabs(q) %/% 2) w else w + 1
}

# ---------------------------------------------------------------------------
# Equation analysis

.dio_parse <- function(eq, vars = NULL) {
  eq <- expand(as_sym(eq))
  if (is.null(vars)) vars <- free_symbols(eq)
  if (!length(vars)) stop("equation should have 1 or more free symbols", call. = FALSE)
  p <- expr_to_qpoly(eq, vars, strict = TRUE)
  if (is.null(p)) stop("Equation should be a polynomial with Rational coefficients.", call. = FALSE)
  if (any(gmp::denominator(p$c) != 1)) stop("Coefficients should be Integers", call. = FALSE)
  E <- p$e
  if (!is.matrix(E)) E <- matrix(E, ncol = length(vars))
  cz <- gmp::as.bigz(gmp::numerator(p$c))
  deg <- if (nrow(E)) rowSums(E) else integer()
  list(eq = eq, vars = vars, E = E, c = cz, n = length(vars),
       total_degree = if (length(deg)) max(deg) else 0L,
       homogeneous = !any(deg == 0L),
       homogeneous_order = !any(deg == 1L),
       has_mul = any(rowSums(E > 0L) >= 2L))
}

# Coefficient of the monomial with exponent vector ex (0 if absent).
.dio_c <- function(info, ex) {
  if (!nrow(info$E)) return(gmp::as.bigz(0))
  hit <- which(apply(info$E, 1L, function(r) all(r == ex)))
  if (length(hit)) info$c[hit[1L]] else gmp::as.bigz(0)
}
.dio_mono <- function(info, ...) {
  ex <- integer(info$n)
  for (i in c(...)) ex[i] <- ex[i] + 1L
  .dio_c(info, ex)
}
.dio_const <- function(info) .dio_c(info, integer(info$n))

.dio_type <- function(info) {
  td <- info$total_degree
  n <- info$n
  E <- info$E
  nonconst <- if (nrow(E)) rowSums(E) > 0L else logical()
  if (td == 1L) return("linear")
  if (n == 1L) return("univariate")
  if (td == 2L && n == 2L) return("binary_quadratic")
  if (td == 2L && n == 3L && info$homogeneous) {
    if (!info$homogeneous_order) return("inhomogeneous_ternary_quadratic")
    sq <- all(vapply(seq_len(n), function(i) { ex <- integer(n); ex[i] <- 2L; .dio_c(info, ex) != 0 }, TRUE))
    nz <- sum(info$c != 0)
    if (nz == 3L && sq) return("homogeneous_ternary_quadratic_normal")
    return("homogeneous_ternary_quadratic")
  }
  if (td == 2L && n >= 3L) {
    if (!info$homogeneous_order) return("inhomogeneous_general_quadratic")
    if (info$has_mul) return(if (info$homogeneous) "homogeneous_general_quadratic" else "inhomogeneous_general_quadratic")
    if (all(info$c[nonconst] == 1)) return("general_sum_of_squares")
    if (all(vapply(seq_along(info$c), function(i) .is_sqz(.zabs(info$c[i])), TRUE)) &&
        abs(sum(vapply(seq_along(info$c), function(i) .zsign(info$c[i]), 0L))) == n - 2L)
      return("general_pythagorean")
  }
  if (td == 3L && n == 2L) return("cubic_thue")
  if (td > 3L && td %% 2L == 0L) {
    pw <- nrow(E) && all(vapply(which(nonconst), function(i) sum(E[i, ] > 0L) == 1L && max(E[i, ]) == td, TRUE))
    if (pw && all(info$c[nonconst] == 1)) return("general_sum_of_even_powers")
  }
  NULL
}

# SymPy DiophantineSolutionSet.add: no -p for a parameter p not present.
.dio_canon <- function(sol) {
  sol <- lapply(sol, as_sym)
  keys <- vapply(sol, .k, "")
  for (i in seq_along(sol)) {
    nx <- .neg(sol[[i]])
    if (is_symbol_like(nx) && !(nx$key %in% keys)) {
      sol <- lapply(sol, function(s) subs(s, nx, sol[[i]]))
      keys <- vapply(sol, .k, "")
    }
  }
  sol
}
.dio_set_add <- function(set, sol) {
  sol <- .dio_canon(sol)
  k <- paste(vapply(sol, .k, ""), collapse = "|")
  set[[k]] <- sol
  set
}

.int_param <- function(name) Symbol(name, integer = TRUE)

# ---------------------------------------------------------------------------
# Linear

.base_solution_linear <- function(c, a, b, t = NULL) {
  abc <- .remove_gcd(.Z(a), .Z(b), .Z(c))
  a <- abc[[1L]]; b <- abc[[2L]]; c <- abc[[3L]]
  if (c == 0) {
    if (is.null(t)) return(list(S_Zero, S_Zero))
    if (b < 0) t <- .neg(t)
    return(list(.mul(list(.zs(b), t)), .mul(list(.zs(-a), t))))
  }
  g <- .gcdext(.zabs(a), .zabs(b))
  x0 <- g[[2L]] * .zsign(a)
  y0 <- g[[3L]] * .zsign(b)
  d <- g[[1L]]
  if (.pmod(c, d) != 0) return(list(NULL, NULL))
  if (is.null(t)) return(list(.zs(c * x0), .zs(c * y0)))
  if (b < 0) t <- .neg(t)
  list(.add(list(.zs(c * x0), .mul(list(.zs(b), t)))), .sub(.zs(c * y0), .mul(list(.zs(a), t))))
}

.dio_linear <- function(info, params) {
  var <- info$vars
  n <- info$n
  c <- .zs(-.dio_const(info))
  A <- lapply(seq_len(n), function(i) .dio_mono(info, i))
  if (n == 1L) {
    q <- .Z(c) %/% A[[1L]]
    r <- .Z(c) - q * A[[1L]]
    return(if (r == 0) list(list(.zs(q))) else list())
  }
  B <- list()
  if (n > 2L) {
    B[[1L]] <- .zgcd(A[[n - 1L]], A[[n]])
    A[[n - 1L]] <- A[[n - 1L]] %/% B[[1L]]
    A[[n]] <- A[[n]] %/% B[[1L]]
    if (n - 2L >= 2L) for (i in (n - 2L):2L) {
      g <- .zgcd(B[[1L]], A[[i]])
      B[[1L]] <- B[[1L]] %/% g
      A[[i]] <- A[[i]] %/% g
      B <- c(list(g), B)
    }
  }
  B[[length(B) + 1L]] <- A[[n]]
  pk <- vapply(params, .k, "")
  solutions <- list()
  for (j in seq_along(B)) {
    Ai <- A[[j]]; Bi <- B[[j]]
    tot_x <- list(); tot_y <- list()
    terms <- if (c$type == "Add") c$args else list(c)
    for (arg in terms) {
      if (arg$type == "Integer") {
        k <- arg; p <- S_One; pnew <- params[[1L]]
      } else {
        sp <- as_coeff_mul_split(arg)
        k <- sp[[1L]]; p <- sp[[2L]]
        pnew <- params[[match(p$key, pk) + 1L]]
      }
      sol <- .base_solution_linear(k, Ai, Bi, pnew)
      if (identical(p$key, S_One$key)) {
        if (is.null(sol[[1L]]) || is.null(sol[[2L]])) return(list())
      } else {
        fix <- function(s) if (s$type == "Add") .add(list(.mul(list(s$args[[1L]], p)), s$args[[2L]])) else s
        sol <- lapply(sol, fix)
      }
      tot_x[[length(tot_x) + 1L]] <- sol[[1L]]
      tot_y[[length(tot_y) + 1L]] <- sol[[2L]]
    }
    solutions[[length(solutions) + 1L]] <- .add(tot_x)
    c <- .add(tot_y)
  }
  solutions[[length(solutions) + 1L]] <- c
  list(solutions)
}

# ---------------------------------------------------------------------------
# Pell type equations x^2 - D*y^2 = N

.square_factor <- function(a) {
  a <- .zabs(a)
  if (a <= 1) return(gmp::as.bigz(1))
  f <- .zfactor(a)
  r <- gmp::as.bigz(1)
  for (i in seq_along(f$p)) r <- r * f$p[i]^(f$e[i] %/% 2L)
  r
}

# Prime factorization as list(p = bigz primes, e = exponents).
.zfactor <- function(a) {
  a <- .zabs(a)
  fs <- gmp::factorize(a)
  if (!length(fs) || (length(fs) == 1L && fs[1L] == 1)) return(list(p = gmp::as.bigz(integer()), e = integer()))
  u <- unique(as.character(fs))
  list(p = gmp::as.bigz(u), e = vapply(u, function(s) sum(as.character(fs) == s), 0L, USE.NAMES = FALSE))
}

.zdivisors <- function(n) {
  n <- .zabs(n)
  if (n == 0) return(gmp::as.bigz(integer()))
  f <- .zfactor(n)
  d <- gmp::as.bigz(1)
  for (i in seq_along(f$p)) {
    pw <- f$p[i]^(0:f$e[i])
    d <- gmp::as.bigz(as.vector(outer(as.character(d), as.character(pw), function(a, b) as.character(gmp::as.bigz(a) * gmp::as.bigz(b)))))
  }
  sort(d)
}

.zsqrt_mod <- function(a, m, all = FALSE) {
  m <- .zabs(m)
  a <- .pmod(.Z(a), m)
  if (m == 1) return(gmp::as.bigz(0))
  r <- sqrt_mod(a, m, all_roots = all)
  if (is.null(r) || !length(r)) return(NULL)
  r <- gmp::as.bigz(r)
  if (all) sort(r) else r[1L]
}

.cornacchia <- function(a, b, m) {
  a <- .Z(a); b <- .Z(b); m <- .Z(m)
  sols <- list()
  addsol <- function(r, s) sols[[paste(r, s)]] <<- list(r, s)
  if (a + b > m) {
    if (a == 1) {
      sr <- .zroot2(m %/% a)
      if (sr[[2L]]) addsol(sr[[1L]], gmp::as.bigz(0))
      if (a == b) return(sols)
    }
    if (.pmod(m, b) == 0) {
      sr <- .zroot2(m %/% b)
      if (sr[[2L]]) addsol(gmp::as.bigz(0), sr[[1L]])
    }
    return(sols)
  }
  inv <- gmp::inv.bigz(.pmod(a, m), m)
  for (t in as.list(.zsqrt_mod(-b * inv, m, all = TRUE))) {
    if (t < m %/% 2) next
    u <- m; r <- t
    while ((m1 <- m - a * r^2) <= 0) {
      tmp <- r; r <- .pmod(u, r); u <- tmp
    }
    q <- m1 %/% b
    if (m1 - q * b != 0) next
    sr <- .zroot2(q)
    if (sr[[2L]]) {
      s <- sr[[1L]]
      if (a == b && r < s) { tmp <- r; r <- s; s <- tmp }
      addsol(r, s)
    }
  }
  sols
}

# Terms of SymPy's PQa generator, computed on demand.
.pqa <- function(P0, Q0, D) {
  sqD <- .isqrtz(D)
  st <- new.env()
  st$A2 <- gmp::as.bigz(0); st$B1 <- gmp::as.bigz(0)
  st$A1 <- gmp::as.bigz(1); st$B2 <- gmp::as.bigz(1)
  st$G1 <- .Z(Q0); st$G2 <- -.Z(P0)
  st$P <- .Z(P0); st$Q <- .Z(Q0)
  function() {
    a <- (st$P + sqD) %/% st$Q
    t <- a * st$A1 + st$A2; st$A2 <- st$A1; st$A1 <- t
    t <- a * st$B1 + st$B2; st$B2 <- st$B1; st$B1 <- t
    t <- a * st$G1 + st$G2; st$G2 <- st$G1; st$G1 <- t
    out <- list(P = st$P, Q = st$Q, a = a, A = st$A1, B = st$B1, G = st$G1)
    st$P <- a * st$Q - st$P
    st$Q <- (D - st$P^2) %/% st$Q
    out
  }
}

# Number of terms (non repeating plus period) of the continued fraction of
# (P + s*sqrt(D))/Q, as SymPy's length() over continued_fraction_periodic().
.cf_length <- function(P, Q, D, s = 1L) {
  P <- .Z(P); Q <- .Z(Q); D <- .Z(D)
  r <- .isqrtz(D)
  if (r * r == D) {
    # rational: count Euclid steps
    a <- P + s * r; b <- Q
    k <- 0L
    while (b != 0) { q <- a %/% b; t <- a - q * b; a <- b; b <- t; k <- k + 1L }
    return(k)
  }
  if (Q < 0) { P <- -P; Q <- -Q; s <- -s }
  # floor of (P + s*sqrt(D))/q for an irrational numerator
  fl <- function(P, q, D) {
    r <- .isqrtz(D)
    L <- P + (if (s > 0) r else -r - 1)
    if (q > 0) L %/% q else (-L - 1) %/% (-q)
  }
  if (s > 0 && P < 0 && P * P > D || s < 0 && (P <= 0 || P * P < D)) {
    # negative value: SymPy expands 1 - frac(-n), which has the same length
    w <- (-fl(P, Q, D)) - 1
    return(.cf_length(P + (w + 1) * Q, Q, D, s))
  }
  if (.pmod(D - P^2, Q) != 0) {
    D <- D * Q^2; P <- P * Q; Q <- Q * Q
  }
  seen <- character()
  repeat {
    k <- paste(P, Q)
    if (k %in% seen) return(length(seen))
    seen <- c(seen, k)
    a <- fl(P, Q, D)
    P <- a * Q - P
    Q <- (D - P^2) %/% Q
  }
}

.special_diop_DN <- function(D, N) {
  sqrtD <- .isqrtz(D)
  Fm <- list()
  for (f in as.list(.zdivisors(.square_factor(.zabs(N))))) Fm[[as.character(N %/% f^2)]] <- f
  P <- gmp::as.bigz(0); Q <- gmp::as.bigz(1)
  G0 <- gmp::as.bigz(0); G1 <- gmp::as.bigz(1)
  B0 <- gmp::as.bigz(1); B1 <- gmp::as.bigz(0)
  sols <- list()
  repeat {
    for (i in 1:2) {
      a <- (P + sqrtD) %/% Q
      P <- a * Q - P
      Q <- (D - P^2) %/% Q
      t <- a * G1 + G0; G0 <- G1; G1 <- t
      t <- a * B1 + B0; B0 <- B1; B1 <- t
      s <- as.character(G1^2 - D * B1^2)
      if (!is.null(Fm[[s]])) sols[[length(sols) + 1L]] <- list(Fm[[s]] * G1, Fm[[s]] * B1)
    }
    if (Q == 1) break
  }
  sols
}

# Solutions of x^2 - D*y^2 = N, as a list of pairs (bigz, or sym when a
# parameter t appears).
.diop_DN <- function(D, N, t = .int_param("t")) {
  D <- .Z(D); N <- .Z(N)
  if (D < 0) {
    if (N == 0) return(list(list(gmp::as.bigz(0), gmp::as.bigz(0))))
    if (N < 0) return(list())
    sol <- list()
    for (d in as.list(.zdivisors(.square_factor(N)))) {
      for (xy in .cornacchia(1, -D, N %/% d^2)) {
        sol[[length(sol) + 1L]] <- list(d * xy[[1L]], d * xy[[2L]])
        if (D == -1) sol[[length(sol) + 1L]] <- list(d * xy[[2L]], d * xy[[1L]])
      }
    }
    return(sol)
  }
  if (D == 0) {
    if (N < 0) return(list())
    if (N == 0) return(list(list(S_Zero, t)))
    sr <- .zroot2(N)
    return(if (sr[[2L]]) list(list(.zs(sr[[1L]]), t)) else list())
  }
  sr <- .zroot2(D)
  sD <- sr[[1L]]
  if (sr[[2L]]) {
    if (N == 0) return(list(list(.mul(list(.zs(sD), t)), t)))
    sol <- list()
    hi <- gmp::as.bigz(floor(gmp::as.bigq(.zsign(N) * (N - 1), 2 * sD)))
    y <- gmp::as.bigz(0)
    while (y <= hi) {
      v <- D * y^2 + N
      if (v >= 0) {
        q <- .zroot2(v)
        if (q[[2L]]) sol[[length(sol) + 1L]] <- list(q[[1L]], y)
      }
      y <- y + 1
    }
    return(sol)
  }
  if (1 < N^2 && N^2 < D) return(.special_diop_DN(D, N))
  if (N == 0) return(list(list(gmp::as.bigz(0), gmp::as.bigz(0))))
  sol <- list()
  if (.zabs(N) == 1) {
    nx <- .pqa(0, 1, D)
    first <- nx()
    prev_B <- first$B; prev_G <- first$G
    j <- -1L
    repeat {
      j <- j + 1L
      tm <- nx()
      if (tm$a == 2 * sD) break
      prev_B <- tm$B; prev_G <- tm$G
    }
    if (j %% 2L) {
      if (N == 1) sol[[1L]] <- list(prev_G, prev_B)
      return(sol)
    }
    if (N == -1) return(list(list(prev_G, prev_B)))
    for (i in seq_len(j)) tm <- nx()
    return(list(list(tm$G, tm$B)))
  }
  for (f in as.list(.zdivisors(.square_factor(N)))) {
    m <- N %/% f^2
    am <- .zabs(m)
    roots <- .zsqrt_mod(D, am, all = TRUE)
    for (sqm in as.list(roots)) {
      z <- .Z(symmetric_residue(sqm, am))
      nx <- .pqa(z, am, D)
      first <- nx()
      prev_B <- first$B; prev_G <- first$G
      L <- .cf_length(z, am, D)
      if (L - 1L >= 1L) for (i in seq_len(L - 1L)) {
        tm <- nx()
        if (.zabs(tm$Q) == 1) {
          if (prev_G^2 - D * prev_B^2 == m) {
            sol[[length(sol) + 1L]] <- list(f * prev_G, f * prev_B)
          } else {
            a <- .diop_DN(D, -1)
            if (length(a)) {
              a0 <- a[[1L]][[1L]]; a1 <- a[[1L]][[2L]]
              sol[[length(sol) + 1L]] <- list(f * (prev_G * a0 + prev_B * D * a1), f * (prev_G * a1 + prev_B * a0))
            }
          }
          break
        }
        prev_B <- tm$B; prev_G <- tm$G
      }
    }
  }
  sol
}

# ---------------------------------------------------------------------------
# Binary quadratic

.dio_qm <- function(...) lapply(list(...), function(v) gmp::as.bigq(v))

# Transformation (P, Q) with (x, y) = P*(X, Y) + Q to the form X^2 - D*Y^2 = N.
.transformation_to_DN <- function(cf) {
  g <- .remove_gcd(cf$a, cf$b, cf$c, cf$d, cf$e, cf$f)
  a <- g[[1L]]; b <- g[[2L]]; c <- g[[3L]]; d <- g[[4L]]; e <- g[[5L]]; f <- g[[6L]]
  mmul <- function(M, N) list(M[[1L]] * N[[1L]] + M[[2L]] * N[[3L]], M[[1L]] * N[[2L]] + M[[2L]] * N[[4L]],
                              M[[3L]] * N[[1L]] + M[[4L]] * N[[3L]], M[[3L]] * N[[2L]] + M[[4L]] * N[[4L]])
  mvec <- function(M, v) list(M[[1L]] * v[[1L]] + M[[2L]] * v[[2L]], M[[3L]] * v[[1L]] + M[[4L]] * v[[2L]])
  if (b != 0) {
    BC <- .rational_pq(2 * a, b); B <- BC[[1L]]; C <- BC[[2L]]
    AT <- .rational_pq(a, B^2); A <- AT[[1L]]; T <- AT[[2L]]
    r <- .transformation_to_DN(list(a = A * B, b = 0, c = B * (c * T - A * C^2), d = d * T,
                                    e = B * e * T - d * T * C, f = f * T * B))
    M <- .dio_qm(gmp::as.bigq(1, B), gmp::as.bigq(-C, B), 0, 1)
    return(list(mmul(M, r[[1L]]), mvec(M, r[[2L]])))
  }
  if (d != 0) {
    BC <- .rational_pq(2 * a, d); B <- BC[[1L]]; C <- BC[[2L]]
    AT <- .rational_pq(a, B^2); A <- AT[[1L]]; T <- AT[[2L]]
    r <- .transformation_to_DN(list(a = A, b = 0, c = c * T, d = 0, e = e * T, f = f * T - A * C^2))
    M <- .dio_qm(gmp::as.bigq(1, B), 0, 0, 1)
    v <- mvec(M, r[[2L]])
    return(list(mmul(M, r[[1L]]), list(v[[1L]] - gmp::as.bigq(C, B), v[[2L]])))
  }
  if (e != 0) {
    BC <- .rational_pq(2 * c, e); B <- BC[[1L]]; C <- BC[[2L]]
    AT <- .rational_pq(c, B^2); A <- AT[[1L]]; T <- AT[[2L]]
    r <- .transformation_to_DN(list(a = a * T, b = 0, c = A, d = 0, e = 0, f = f * T - A * C^2))
    M <- .dio_qm(1, 0, 0, gmp::as.bigq(1, B))
    v <- mvec(M, r[[2L]])
    return(list(mmul(M, r[[1L]]), list(v[[1L]], v[[2L]] - gmp::as.bigq(C, B))))
  }
  list(.dio_qm(gmp::as.bigq(1, a), 0, 0, 1), .dio_qm(0, 0))
}

.find_DN <- function(cf) {
  PQ <- .transformation_to_DN(cf)
  P <- PQ[[1L]]; Q <- PQ[[2L]]
  # substitute x = P11 X + P12 Y + Q1, y = P21 X + P22 Y + Q2 and read off
  # the coefficients of X^2, Y^2 and 1
  u <- list(P[[1L]], P[[2L]], Q[[1L]])
  v <- list(P[[3L]], P[[4L]], Q[[2L]])
  A <- gmp::as.bigq(cf$a); B <- gmp::as.bigq(cf$b); C <- gmp::as.bigq(cf$c)
  Dd <- gmp::as.bigq(cf$d); E <- gmp::as.bigq(cf$e); F <- gmp::as.bigq(cf$f)
  cX2 <- A * u[[1L]]^2 + B * u[[1L]] * v[[1L]] + C * v[[1L]]^2
  cY2 <- A * u[[2L]]^2 + B * u[[2L]] * v[[2L]] + C * v[[2L]]^2
  c1 <- A * u[[3L]]^2 + B * u[[3L]] * v[[3L]] + C * v[[3L]]^2 + Dd * u[[3L]] + E * v[[3L]] + F
  list(-cY2 / cX2, -c1 / cX2)
}

.bq_coeffs <- function(info) {
  list(a = .dio_mono(info, 1L, 1L), b = .dio_mono(info, 1L, 2L), c = .dio_mono(info, 2L, 2L),
       d = .dio_mono(info, 1L), e = .dio_mono(info, 2L), f = .dio_const(info))
}

.is_int_expr <- function(e) as_sym(e)$type == "Integer"

.dio_binary_quadratic <- function(info, params) {
  var <- info$vars
  x <- var[[1L]]; y <- var[[2L]]
  cf <- .bq_coeffs(info)
  g <- .remove_gcd(cf$a, cf$b, cf$c, cf$d, cf$e, cf$f)
  A <- g[[1L]]; B <- g[[2L]]; C <- g[[3L]]; D <- g[[4L]]; E <- g[[5L]]; F <- g[[6L]]
  t <- params[[1L]]; u <- params[[2L]]
  out <- list()
  add <- function(s) out <<- .dio_set_add(out, s)
  discr <- B^2 - 4 * A * C
  if (A == 0 && C == 0 && B != 0) {
    if (D * E - B * F == 0) {
      q <- E %/% B
      if (E - q * B == 0) add(list(.zs(-q), t))
      q <- D %/% B
      if (D - q * B == 0) add(list(t, .zs(-q)))
    } else {
      dv <- .zdivisors(D * E - B * F)
      dv <- c(as.list(dv), lapply(as.list(dv), function(z) -z))
      for (d in dv) {
        x0 <- (d - E) %/% B
        if ((d - E) - x0 * B != 0) next
        q <- (D * E - B * F) %/% d
        if ((D * E - B * F) - q * d != 0) next
        y0 <- (q - D) %/% B
        if ((q - D) - y0 * B == 0) add(list(.zs(x0), .zs(y0)))
      }
    }
  } else if (discr == 0) {
    if (A == 0) {
      sw <- .dio_parse(info$eq, list(y, x))
      for (s in .dio_binary_quadratic(sw, list(t, u))) add(list(s[[2L]], s[[1L]]))
    } else {
      gg <- .zsign(A) * .zgcd(A, C)
      a <- A %/% gg; c <- C %/% gg
      e <- .zsign(B) * .zsign(A)
      sqa <- .isqrtz(a); sqc <- .isqrtz(c)
      cc <- e * sqc * D - sqa * E
      if (cc == 0) {
        z <- Symbol("z", real = TRUE)
        zeq <- .add(list(.mul(list(.zs(sqa * gg), .pow(z, S_Two))), .mul(list(.zs(D), z)), .zs(sqa * F)))
        rts <- Filter(function(r) r$type == "Integer", .unique_syms(roots(zeq, z)))
        for (root in rts) {
          ans <- .dio_solve_eq(.sub(.add(list(.mul(list(.zs(sqa), x)), .mul(list(.zs(e * sqc), y)))), root), .int_param("t"))
          add(list(ans[[1L]], ans[[2L]]))
        }
      } else {
        for (z0 in seq_len(as.integer(.zabs(cc))) - 1L) {
          z0 <- gmp::as.bigz(z0)
          if (.pmod(sqa * gg * z0^2 + D * z0 + sqa * F, cc) == 0 &&
              .pmod(e * sqc * gg * z0^2 + E * z0 + e * sqc * F, cc) == 0) {
            sx <- .add(list(.mul(list(.zs(-e * sqc * gg * cc), .pow(t, S_Two))),
                            .mul(list(.zs(-(E + 2 * e * sqc * gg * z0)), t)),
                            .zs(-((e * sqc * gg * z0^2 + E * z0 + e * sqc * F) %/% cc))))
            sy <- .add(list(.mul(list(.zs(sqa * gg * cc), .pow(t, S_Two))),
                            .mul(list(.zs(D + 2 * sqa * gg * z0), t)),
                            .zs((sqa * gg * z0^2 + D * z0 + sqa * F) %/% cc)))
            add(list(sx, sy))
          }
        }
      }
    }
  } else if (.is_sqz(discr)) {
    if (A != 0) {
      r <- .isqrtz(discr)
      uu <- .int_param("u"); vv <- .int_param("v")
      eq2 <- expand(.add(list(.mul(list(.zs(4 * A * r), uu, vv)),
                              .mul(list(.zs(4 * A * D), .add(list(.mul(list(.zs(B), vv)), .mul(list(.zs(r), uu)),
                                                                  .mul(list(.zs(r), vv)), .mul(list(.zs(-B), uu)))))),
                              .mul(list(.zs(8 * A^2 * E), .sub(uu, vv))),
                              .zs(16 * A^2 * r * F))))
      for (st in .dio_solve_eq(eq2, t)) {
        s0 <- st[[1L]]; t0 <- st[[2L]]
        num <- expand(.add(list(.mul(list(.zs(B), t0)), .mul(list(.zs(r), s0)), .mul(list(.zs(r), t0)), .mul(list(.zs(-B), s0)))))
        x_0 <- .div(num, .zs(4 * A * r))
        y_0 <- .div(.sub(s0, t0), .zs(2 * r))
        if (length(free_symbols(s0)) || length(free_symbols(t0))) {
          for (s in .check_param(x_0, y_0, 4 * A * r, params)) add(s)
        } else if (.is_int_expr(x_0) && .is_int_expr(y_0)) {
          if (is_zero_num(expand(subs(info$eq, list(list(x, x_0), list(y, y_0)))))) add(list(x_0, y_0))
        }
      }
    } else {
      sw <- .dio_parse(info$eq, list(y, x))
      for (s in .dio_binary_quadratic(sw, list(t, u))) add(rev(s))
    }
  } else {
    PQ <- .transformation_to_DN(cf)
    P <- PQ[[1L]]; Q <- PQ[[2L]]
    DN <- .find_DN(cf)
    Dq <- DN[[1L]]; Nq <- DN[[2L]]
    Dz <- gmp::as.bigz(gmp::numerator(Dq)); Nz <- gmp::as.bigz(gmp::numerator(Nq))
    pell <- .diop_DN(Dz, Nz)
    apply_PQ <- function(X, Y) list(
      .add(list(.mul(list(make_rational(P[[1L]]), X)), .mul(list(make_rational(P[[2L]]), Y)), make_rational(Q[[1L]]))),
      .add(list(.mul(list(make_rational(P[[3L]]), X)), .mul(list(make_rational(P[[4L]]), Y)), make_rational(Q[[2L]]))))
    if (Dz < 0) {
      for (p in pell) for (X in list(-p[[1L]], p[[1L]])) for (Y in list(-p[[2L]], p[[2L]])) {
        s <- apply_PQ(.zs(X), .zs(Y))
        if (all(vapply(s, .is_int_expr, TRUE))) add(s)
      }
    } else {
      keys <- character()
      ps <- list()
      for (p in c(pell, lapply(pell, function(p) list(-p[[1L]], -p[[2L]])))) {
        k <- paste(p[[1L]], p[[2L]])
        if (!(k %in% keys)) { keys <- c(keys, k); ps[[length(ps) + 1L]] <- p }
      }
      a <- .diop_DN(Dz, 1)
      T <- a[[1L]][[1L]]; U <- a[[1L]][[2L]]
      sD <- .sqrt(.zs(Dz))
      allint <- all(vapply(c(P, Q), function(q) gmp::denominator(q) == 1, TRUE))
      pell_t <- function(X, Y, Tk, Uk, ex = TRUE) {
        aa <- .mul(list(.add(list(.zs(X), .mul(list(.zs(Y), sD)))), .pow(.add(list(.zs(Tk), .mul(list(.zs(Uk), sD)))), t)))
        bb <- .mul(list(.sub(.zs(X), .mul(list(.zs(Y), sD))), .pow(.sub(.zs(Tk), .mul(list(.zs(Uk), sD))), t)))
        r <- list(.div(.add(list(aa, bb)), S_Two), .div(.sub(aa, bb), .mul(list(S_Two, sD))))
        if (ex) lapply(r, expand) else r
      }
      if (allint) {
        for (p in ps) {
          xy <- pell_t(p[[1L]], p[[2L]], T, U)
          add(lapply(apply_PQ(xy[[1L]], xy[[2L]]), expand))
        }
      } else {
        L <- do.call(.zlcm, lapply(c(P, Q), function(q) gmp::as.bigz(gmp::denominator(q))))
        k <- 1L; Tk <- T; Uk <- U
        while (.pmod(Tk - 1, L) != 0 || .pmod(Uk, L) != 0) {
          t1 <- Tk * T + Dz * Uk * U
          Uk <- Tk * U + Uk * T
          Tk <- t1
          k <- k + 1L
        }
        for (p in ps) {
          X <- p[[1L]]; Y <- p[[2L]]
          for (i in seq_len(k)) {
            s0 <- apply_PQ(.zs(X), .zs(Y))
            if (all(vapply(s0, .is_int_expr, TRUE))) {
              xy <- pell_t(X, Y, Tk, Uk, ex = FALSE)
              add(apply_PQ(xy[[1L]], xy[[2L]]))
            }
            t1 <- X * T + Dz * U * Y
            Y <- X * U + Y * T
            X <- t1
          }
        }
      }
    }
  }
  unname(out)
}

.unique_syms <- function(lst) lst[!duplicated(vapply(lst, .k, ""))]

# Positive rational content: gcd of the numerators over the lcm of the
# denominators of the numeric term coefficients.
.q_content <- function(e) {
  e <- expand(as_sym(e))
  terms <- if (e$type == "Add") e$args else list(e)
  q <- gmp::as.bigz(1); pn <- gmp::as.bigz(0)
  for (tm in terms) {
    c <- as_coeff_mul_split(tm)[[1L]]
    v <- if (is_exact_number(c)) gmp::as.bigq(c$value) else gmp::as.bigq(1)
    q <- gmp::lcm.bigz(q, gmp::denominator(v))
    pn <- gmp::gcd.bigz(pn, gmp::numerator(v))
  }
  if (pn == 0) gmp::as.bigq(1) else gmp::as.bigq(pn, q)
}

# SymPy clear_coefficients(expr, rhs)[1]: divide out the rational content
# and move the rational constant of expr, doing the same to rhs.
.clear_coefficients <- function(expr, rhs) {
  expr <- expand(as_sym(expr))
  rhs <- as_sym(rhs)
  if (is_exact_number(expr)) return(expand(.sub(rhs, expr)))
  was <- ""
  while (!is_zero_num(expr) && !identical(was, expr$key)) {
    was <- expr$key
    m <- make_rational(.q_content(expr))
    expr <- expand(.div(expr, m))
    rhs <- .div(rhs, m)
    terms <- if (expr$type == "Add") expr$args else list(expr)
    cst <- .add(Filter(is_exact_number, terms))
    if (!is_zero_num(cst)) {
      expr <- expand(.sub(expr, cst))
      rhs <- .sub(rhs, cst)
    }
  }
  if (could_extract_minus_sign(expr)) rhs <- .neg(rhs)
  expand(rhs)
}

# SymPy check_param: integer points of a rational parametric solution. The
# unknowns m and n stand for the values of x and y.
.check_param <- function(x, y, a, params) {
  x <- as_sym(x); y <- as_sym(y)
  if (!length(free_symbols(x)) && x$type != "Integer") return(list())
  if (!length(free_symbols(y)) && y$type != "Integer") return(list())
  m <- .int_param("m"); n <- .int_param("n")
  c <- .q_content(.add(list(.mul(list(m, x)), .mul(list(n, y)))))
  if (.pmod(.Z(a), gmp::as.bigz(gmp::denominator(c))) != 0) return(list())
  eq <- expand(.sub(.clear_coefficients(x, m), .clear_coefficients(y, n)))
  eq <- expand(.div(eq, make_rational(.q_content(eq))))
  info <- .dio_parse(eq)
  ty <- .dio_type(info)
  if (is.null(ty)) return(list())
  .dio_run(info, ty, params)
}

# ---------------------------------------------------------------------------
# Ternary quadratics

.sqf_normal <- function(a, b, c) {
  ABC <- .remove_gcd(.Z(a), .Z(b), .Z(c))
  sq <- lapply(ABC, .square_factor)
  sqf <- Map(function(i, j) i %/% j^2, ABC, sq)
  A <- gmp::as.bigq(sqf[[1L]]); B <- gmp::as.bigq(sqf[[2L]]); C <- gmp::as.bigq(sqf[[3L]])
  zq <- function(v) gmp::as.bigz(gmp::numerator(v))
  pc <- .zgcd(zq(A), zq(B)); A <- A / pc; B <- B / pc
  pa <- .zgcd(zq(B), zq(C)); B <- B / pa; C <- C / pa
  pb <- .zgcd(zq(A), zq(C)); A <- A / pb; B <- B / pb
  A <- A * pa; B <- B * pb; C <- C * pc
  list(sq, sqf, list(zq(A), zq(B), zq(C)))
}

.reconstruct <- function(A, B, z) {
  f <- .zfactor(.zgcd(A, B))
  for (i in seq_along(f$p)) {
    if (f$e[i] != 1L) stop("a and b should be square-free", call. = FALSE)
    z <- z * f$p[i]
  }
  z
}

.gaussian_reduce <- function(w, a, b) {
  a <- .zabs(a)
  dot <- function(u, v) u[[1L]] * v[[1L]] + a * u[[2L]] * v[[2L]]
  u <- list(.Z(b), gmp::as.bigz(0))
  v <- if (.Z(b) * .Z(w) >= 0) list(.Z(w), gmp::as.bigz(1)) else list(-.Z(w), gmp::as.bigz(-1))
  if (.Z(b)^2 < .Z(w)^2 + a) { tmp <- u; u <- v; v <- tmp }
  while (dot(u, u) > (dv <- dot(v, v))) {
    k <- dot(u, v) %/% dv
    nv <- list(u[[1L]] - k * v[[1L]], u[[2L]] - k * v[[2L]])
    u <- v; v <- nv
  }
  cc <- list(v[[1L]] - u[[1L]], v[[2L]] - u[[2L]])
  if (dot(cc, cc) <= dot(u, u) && dot(u, u) <= 2 * dot(u, v)) return(cc)
  u
}

.descent <- function(A, B) {
  A <- .Z(A); B <- .Z(B)
  if (.zabs(A) > .zabs(B)) {
    r <- .descent(B, A)
    return(list(r[[1L]], r[[3L]], r[[2L]]))
  }
  if (B == 1) return(.qz(1, 0, 1))
  if (A == 1) return(.qz(1, 1, 0))
  if (B == -A) return(.qz(0, 1, 1))
  if (B == A) {
    r <- .descent(-1, A)
    return(list(A * r[[3L]], r[[2L]], r[[1L]]))
  }
  w <- .zsqrt_mod(A, B)
  xz <- .gaussian_reduce(w, A, B)
  x_0 <- xz[[1L]]; z_0 <- xz[[2L]]
  t <- (x_0^2 - A * z_0^2) %/% B
  t_2 <- .square_factor(t)
  t_1 <- t %/% t_2^2
  r <- .descent(A, t_1)
  x_1 <- r[[1L]]; z_1 <- r[[2L]]; y_1 <- r[[3L]]
  .remove_gcd(x_0 * x_1 + A * z_0 * z_1, z_0 * x_1 + x_0 * z_1, t_1 * t_2 * y_1)
}
.qz <- function(...) lapply(list(...), gmp::as.bigz)

.holzer <- function(x, y, z, a, b, c) {
  k <- if (c %% 2 != 0) 2 * c else c %/% 2
  small <- a * b * c
  step <- 0L
  repeat {
    t1 <- a * x^2; t2 <- b * y^2; t3 <- c * z^2
    if (t1 + t2 != t3) {
      if (step == 0L) stop("bad starting solution", call. = FALSE)
      break
    }
    x_0 <- x; y_0 <- y; z_0 <- z
    if (max(c(t1, t2, t3)) <= small) break
    uv <- .base_solution_linear(k, y_0, -x_0)
    if (is.null(uv[[1L]]) || is.null(uv[[2L]])) break
    u <- .Z(uv[[1L]]); v <- .Z(uv[[2L]])
    p <- -(a * u * x_0 + b * v * y_0); q <- c * z_0
    if (c %% 2 == 0) {
      w <- .nint_or_floor(p, q)
    } else {
      w <- p %/% q
      if ((a * u + b * v + c * w) %% 2 != 0) w <- w + 1
    }
    A <- a * u^2 + b * v^2 + c * w^2
    B <- a * u * x_0 + b * v * y_0 + c * w * z_0
    x <- (x_0 * A - 2 * u * B) %/% k
    y <- (y_0 * A - 2 * v * B) %/% k
    z <- (z_0 * A - 2 * w * B) %/% k
    step <- step + 1L
  }
  list(x_0, y_0, z_0)
}

.tern_c <- function(cf, i, j) {
  k <- paste(sort(c(i, j)), collapse = "")
  v <- cf[[k]]
  if (is.null(v)) gmp::as.bigz(0) else v
}
# Ternary coefficients by variable index: "11" for x^2, "12" for x*y, ...
.tern_coeffs <- function(info) {
  out <- list()
  for (i in 1:3) for (j in i:3) out[[paste0(i, j)]] <- .dio_mono(info, i, j)
  out
}

.tern_normal <- function(cf) {
  a <- .tern_c(cf, 1, 1); b <- .tern_c(cf, 2, 2); c <- .tern_c(cf, 3, 3)
  sn <- .sqf_normal(a, b, c)
  sqf_a <- sn[[1L]][[1L]]; sqf_b <- sn[[1L]][[2L]]; sqf_c <- sn[[1L]][[3L]]
  a_1 <- sn[[2L]][[1L]]; b_1 <- sn[[2L]][[2L]]; c_1 <- sn[[2L]][[3L]]
  a_2 <- sn[[3L]][[1L]]; b_2 <- sn[[3L]][[2L]]; c_2 <- sn[[3L]][[3L]]
  A <- -a_2 * c_2
  B <- -b_2 * c_2
  if (A < 0 && B < 0) return(NULL)
  sm <- function(v, m) is.null(.zsqrt_mod(v, m))
  if (sm(-b_2 * c_2, a_2) || sm(-c_2 * a_2, b_2) || sm(-a_2 * b_2, c_2)) return(NULL)
  d <- .descent(A, B)
  z_0 <- d[[1L]]; x_0 <- d[[2L]]; y_0 <- d[[3L]]
  zq <- .rational_pq(z_0, .zabs(c_2))
  z_0 <- zq[[1L]]; q <- zq[[2L]]
  x_0 <- x_0 * q; y_0 <- y_0 * q
  r <- .remove_gcd(x_0, y_0, z_0)
  x_0 <- r[[1L]]; y_0 <- r[[2L]]; z_0 <- r[[3L]]
  if (.zsign(a) == .zsign(b)) {
    h <- .holzer(x_0, y_0, z_0, .zabs(a_2), .zabs(b_2), .zabs(c_2))
    x_0 <- h[[1L]]; y_0 <- h[[2L]]; z_0 <- h[[3L]]
  } else if (.zsign(a) == .zsign(c)) {
    h <- .holzer(x_0, z_0, y_0, .zabs(a_2), .zabs(c_2), .zabs(b_2))
    x_0 <- h[[1L]]; z_0 <- h[[2L]]; y_0 <- h[[3L]]
  } else {
    h <- .holzer(y_0, z_0, x_0, .zabs(b_2), .zabs(c_2), .zabs(a_2))
    y_0 <- h[[1L]]; z_0 <- h[[2L]]; x_0 <- h[[3L]]
  }
  x_0 <- .reconstruct(b_1, c_1, x_0)
  y_0 <- .reconstruct(a_1, c_1, y_0)
  z_0 <- .reconstruct(a_1, b_1, z_0)
  L <- .zlcm(sqf_a, sqf_b, sqf_c)
  .remove_gcd(.zabs(x_0 * L %/% sqf_a), .zabs(y_0 * L %/% sqf_b), .zabs(z_0 * L %/% sqf_c))
}

# Permute the variables of a ternary coefficient table.
.tern_perm <- function(cf, p) {
  out <- list()
  for (i in 1:3) for (j in i:3) out[[paste(sort(c(i, j)), collapse = "")]] <- .tern_c(cf, p[i], p[j])
  out
}

.is_tern_normal <- function(cf) {
  nz <- sum(vapply(cf, function(v) v != 0, TRUE))
  nz == 3L && all(vapply(1:3, function(i) .tern_c(cf, i, i) != 0, TRUE))
}

# A particular solution of a homogeneous ternary quadratic, or NULL.
.tern_solve <- function(cf) {
  if (.is_tern_normal(cf)) return(.tern_normal(cf))
  sw <- function(p) {
    r <- .tern_solve(.tern_perm(cf, p))
    if (is.null(r)) return(NULL)
    out <- vector("list", 3L)
    out[p] <- r
    out
  }
  sq <- vapply(1:3, function(i) .tern_c(cf, i, i) != 0, TRUE)
  if (!any(sq)) {
    if (.tern_c(cf, 1, 3) != 0) {
      x <- .int_param("x"); z <- .int_param("z")
      sols <- .diophantine_set(.sub(.add(list(.mul(list(.zs(.tern_c(cf, 1, 2)), x)), .mul(list(.zs(.tern_c(cf, 2, 3)), z)))),
                                    .mul(list(x, z))), list(x, z))
      sc <- vapply(sols, function(s) as.numeric(abs(.Z(s[[1L]])) + abs(.Z(s[[2L]]))), 0)
      s <- sols[[which.min(sc)]]
      return(.remove_gcd(.Z(s[[1L]]), -.tern_c(cf, 1, 3), .Z(s[[2L]])))
    }
    return(sw(c(2L, 1L, 3L)))
  }
  if (.tern_c(cf, 1, 1) == 0) {
    if (.tern_c(cf, 2, 2) == 0) return(sw(c(3L, 2L, 1L)))
    return(sw(c(2L, 1L, 3L)))
  }
  if (.tern_c(cf, 1, 2) != 0 || .tern_c(cf, 1, 3) != 0) {
    A <- .tern_c(cf, 1, 1); B <- .tern_c(cf, 1, 2); C <- .tern_c(cf, 1, 3)
    D <- .tern_c(cf, 2, 2); E <- .tern_c(cf, 2, 3); F <- .tern_c(cf, 3, 3)
    nc <- list("11" = 4 * A^2, "22" = 4 * A * D - B^2, "33" = 4 * A * F - C^2,
               "23" = 4 * A * E - 2 * B * C, "12" = gmp::as.bigz(0), "13" = gmp::as.bigz(0))
    r <- .tern_solve(nc)
    if (is.null(r)) return(NULL)
    pq <- .rational_pq(B * r[[2L]] + C * r[[3L]], 2 * A)
    return(.remove_gcd(r[[1L]] * pq[[2L]] - pq[[1L]], r[[2L]] * pq[[2L]], r[[3L]] * pq[[2L]]))
  }
  if (.tern_c(cf, 2, 3) != 0) {
    if (.tern_c(cf, 2, 2) == 0) {
      if (.tern_c(cf, 3, 3) == 0) {
        ba <- .rational_pq(-.tern_c(cf, 2, 3), .tern_c(cf, 1, 1))
        return(.remove_gcd(ba[[1L]], ba[[2L]], ba[[1L]]))
      }
      return(sw(c(3L, 2L, 1L)))
    }
    return(sw(c(2L, 1L, 3L)))
  }
  .tern_normal(cf)
}

.parametrize_ternary <- function(sol, info) {
  if (is.null(sol)) return(list(NULL, NULL, NULL))
  if (sum(vapply(sol, function(v) v == 0, TRUE)) >= 2L) return(list(NULL, NULL, NULL))
  v <- info$vars
  if (sol[[1L]] == 0) {
    info2 <- info
    info2$vars <- v[c(2L, 1L, 3L)]
    r <- .parametrize_ternary(sol[c(2L, 1L, 3L)], info2)
    return(r[c(2L, 1L, 3L)])
  }
  x0 <- .zs(sol[[1L]]); y0 <- .zs(sol[[2L]]); z0 <- .zs(sol[[3L]])
  r <- .int_param("r"); p <- .int_param("p"); q <- .int_param("q")
  e1 <- expand(xreplace(info$eq, list(list(v[[1L]], .mul(list(r, x0))), list(v[[2L]], .add(list(.mul(list(r, y0)), p))),
                                      list(v[[3L]], .add(list(.mul(list(r, z0)), q))))))
  ind <- .as_independent(e1, r)
  A <- ind[[1L]]; B <- ind[[2L]]
  .remove_gcd(expand(.mul(list(A, x0))), expand(.sub(.mul(list(A, y0)), expand(.div(.mul(list(B, p)), r)))),
              expand(.sub(.mul(list(A, z0)), expand(.div(.mul(list(B, q)), r)))))
}

# ---------------------------------------------------------------------------
# Sums of squares and powers

.prime_as_sum_of_two_squares <- function(p) {
  p <- .Z(p)
  if (p %% 4 != 1) return(NULL)
  b <- if (p %% 8 == 5) 2 else if (p %% 12 == 5) 3 else if ((p %% 5) %in% c(2, 3)) 5 else {
    bb <- 7
    while (.dio_jacobi(bb, p) == 1) bb <- as.numeric(.Z(nextprime(bb)))
    bb
  }
  b <- gmp::powm(gmp::as.bigz(b), p %/% 4, p)
  a <- p
  while (b^2 > p) { t <- .pmod(a, b); a <- b; b <- t }
  list(.pmod(a, b), b)
}

.dio_jacobi <- function(a, n) {
  a <- .pmod(.Z(a), .Z(n)); n <- .Z(n)
  r <- 1L
  while (a != 0) {
    while (a %% 2 == 0) {
      a <- a %/% 2
      if ((n %% 8) %in% c(3, 5)) r <- -r
    }
    t <- a; a <- n; n <- t
    if (a %% 4 == 3 && n %% 4 == 3) r <- -r
    a <- .pmod(a, n)
  }
  if (n == 1) r else 0L
}

.remove4 <- function(n) {
  v <- 0L
  while (n != 0 && n %% 4 == 0) { n <- n %/% 4; v <- v + 1L }
  list(n, v)
}

.sum_of_three_squares <- function(n) {
  special <- list("1" = c(0, 0, 1), "2" = c(0, 1, 1), "3" = c(1, 1, 1), "10" = c(0, 1, 3), "34" = c(3, 3, 4),
                  "58" = c(0, 3, 7), "85" = c(0, 6, 7), "130" = c(0, 3, 11), "214" = c(3, 6, 13),
                  "226" = c(8, 9, 9), "370" = c(8, 9, 15), "526" = c(6, 7, 21), "706" = c(15, 15, 16),
                  "730" = c(0, 1, 27), "1414" = c(6, 17, 33), "1906" = c(13, 21, 36), "2986" = c(21, 32, 39),
                  "9634" = c(56, 57, 57))
  n <- .Z(n)
  if (n < 0) stop("n should be a non-negative integer", call. = FALSE)
  if (n == 0) return(.qz(0, 0, 0))
  rv <- .remove4(n); n <- rv[[1L]]; v <- gmp::as.bigz(2)^rv[[2L]]
  if (n %% 8 == 7) return(NULL)
  sp <- special[[as.character(n)]]
  if (!is.null(sp)) return(lapply(sp, function(i) v * i))
  sr <- .zroot2(n)
  s <- sr[[1L]]
  if (sr[[2L]]) return(list(gmp::as.bigz(0), gmp::as.bigz(0), v * s))
  srt <- function(l) l[order(vapply(l, function(z) as.numeric(z), 0))]
  if (n %% 8 == 3) {
    if (s %% 2 == 0) s <- s - 1
    x <- s
    while (x >= 0) {
      N <- (n - x^2) %/% 2
      if (isprime(N)) {
        yz <- .prime_as_sum_of_two_squares(N)
        return(srt(list(v * x, v * (yz[[1L]] + yz[[2L]]), v * .zabs(yz[[1L]] - yz[[2L]]))))
      }
      x <- x - 2
    }
    stop("sum_of_three_squares failed", call. = FALSE)
  }
  if (!xor(n %% 2 == 1, s %% 2 == 1)) s <- s - 1
  x <- s
  while (x >= 0) {
    N <- n - x^2
    if (isprime(N)) {
      yz <- .prime_as_sum_of_two_squares(N)
      return(srt(list(v * x, v * yz[[1L]], v * yz[[2L]])))
    }
    x <- x - 2
  }
  stop("sum_of_three_squares failed", call. = FALSE)
}

.can_do_sum_of_squares <- function(n, k) {
  if (k < 1L || n < 0) return(FALSE)
  if (n == 0) return(TRUE)
  if (k == 1L) return(.is_sqz(n))
  if (k == 2L) {
    if (n %in% c(1, 2)) return(TRUE)
    if (isprime(n)) return(if (n %% 4 == 1) 1L else FALSE)
    f <- .zfactor(n)
    return(all(vapply(seq_along(f$p), function(i) f$p[i] %% 4 != 3 || f$e[i] %% 2L == 0L, TRUE)))
  }
  if (k == 3L) return(.remove4(n)[[1L]] %% 8 != 7)
  TRUE
}

# Integer p-th root floor and exactness.
.zroot <- function(n, p) {
  n <- .Z(n)
  if (p == 2L) return(.zroot2(n))
  unname(.iroot(n, p))
}

# All representations of n as k p-th powers, in SymPy's generation order;
# stops after limit results.
.power_representation <- function(n, p, k, zeros = FALSE, limit = Inf) {
  n <- .Z(n); p <- as.integer(p); k <- as.integer(k)
  out <- list()
  emit <- function(t) {
    out[[length(out) + 1L]] <<- t
    length(out) >= limit
  }
  if (n < 0) {
    if (p %% 2L) for (t in .power_representation(-n, p, k, zeros, limit)) if (emit(lapply(t, function(i) -i))) break
    return(out)
  }
  if (p < 1L || k < 1L) stop("Expecting positive integers for `(p, k)`", call. = FALSE)
  if (n == 0) {
    if (zeros) emit(rep(list(gmp::as.bigz(0)), k))
    return(out)
  }
  if (k == 1L) {
    if (p == 1L) emit(list(n))
    else if (n == 1) emit(list(gmp::as.bigz(1)))
    else {
      be <- perfect_power(n)
      if (!isFALSE(be) && !is.null(be) && length(be) == 2L) {
        b <- .Z(be[[1L]]); e <- as.integer(be[[2L]])
        if (e %% p == 0L) emit(list(b^(e %/% p)))
      }
    }
    return(out)
  }
  if (p == 1L) {
    parts <- if (!zeros) .ordered_partitions(as.integer(as.character(n)), k)
             else do.call(c, lapply(seq_len(k), function(m) lapply(.ordered_partitions(as.integer(as.character(n)), m),
                                                                     function(a) c(rep(0L, k - m), a))))
    for (a in parts) if (emit(lapply(a, gmp::as.bigz))) break
    return(out)
  }
  if (p == 2L) {
    if (k == 3L) {
      rv <- .remove4(n)
      if (rv[[2L]] > 0L) {
        v <- gmp::as.bigz(2)^rv[[2L]]
        for (t in .power_representation(rv[[1L]], p, k, zeros, limit)) if (emit(lapply(t, function(i) i * v))) break
        return(out)
      }
    }
    feasible <- .can_do_sum_of_squares(n, k)
    if (isFALSE(feasible)) return(out)
    if (!zeros) {
      if (n > 33 && k >= 5L && k <= n && (n - k) %in% c(13, 10, 7, 5, 4, 2, 1)) return(out)
      if (k == 4L && (n %in% c(1, 3, 5, 9, 11, 17, 29, 41) || .remove4(n)[[1L]] %in% c(2, 6, 14))) return(out)
      if (k == 3L && n %in% c(1, 2, 5, 10, 13, 25, 37, 58, 85, 130)) return(out)
    }
    if (!isTRUE(feasible)) {
      emit(.prime_as_sum_of_two_squares(n))
      return(out)
    }
  }
  if (k == 2L && p > 2L) {
    be <- perfect_power(n)
    if (!isFALSE(be) && !is.null(be) && length(be) == 2L && as.integer(be[[2L]]) %% p == 0L) return(out)
  }
  done <- FALSE
  rec <- function(n_i, kk, rem, terms) {
    if (done) return(invisible())
    if (n_i <= 0 || kk <= 0L) return(invisible())
    if (rem < kk) return(invisible())
    if (kk * n_i^p < rem) return(invisible())
    if (kk == 1L) {
      r <- .zroot(rem, p)
      if (r[[2L]] && r[[1L]] <= n_i) {
        if (emit(rev(c(terms, list(r[[1L]]))))) done <<- TRUE
      }
      return(invisible())
    }
    nt <- gmp::as.bigz(1)
    while (nt <= n_i) {
      res <- rem - nt^p
      if (res < 0) break
      rec(nt, kk - 1L, res, c(terms, list(nt)))
      if (done) return(invisible())
      nt <- nt + 1
    }
  }
  if (n >= k) rec(.zroot(n - (k - 1L), p)[[1L]], k, n, list())
  if (zeros && !done) {
    a <- .zroot(n, p)[[1L]]
    for (i in seq_len(k - 1L)) {
      zpad <- rep(list(gmp::as.bigz(0)), k - i)
      emit0 <- emit
      emit <- function(t) emit0(c(zpad, t))
      rec(a, i, n, list())
      emit <- emit0
      if (done) break
    }
  }
  out
}

# SymPy ordered_partitions(n, m) with sort = TRUE: partitions of n into m
# parts in nondecreasing order, in SymPy's generation order.
.ordered_partitions <- function(n, m) {
  if (n < 1L || m < 1L) return(list(integer()))
  if (m == 1L) return(list(n))
  if (n == m) return(list(rep(1L, n)))
  out <- list()
  for (b in seq_len(n %/% m)) {
    a <- rep(b, m)
    x <- n - b * m
    if (!x) {
      out[[length(out) + 1L]] <- a
      next
    }
    for (mi in seq_len(m - 1L)) for (ax in .ordered_partitions(x, mi)) {
      if (!length(ax)) next
      a2 <- a
      a2[(m - mi + 1L):m] <- ax + b
      out[[length(out) + 1L]] <- a2
    }
  }
  out
}

.permute_signs <- function(t) {
  nz <- which(vapply(t, function(v) !is_zero_num(as_sym(v)), TRUE))
  out <- list()
  for (m in 0:(2^length(nz) - 1)) {
    s <- t
    for (j in seq_along(nz)) if (bitwAnd(m, 2^(j - 1)) != 0) s[[nz[j]]] <- .neg(as_sym(s[[nz[j]]]))
    out[[length(out) + 1L]] <- lapply(s, as_sym)
  }
  out
}

.signed_permutations <- function(t) {
  n <- length(t)
  perms <- if (n == 1L) list(1L) else {
    rec <- function(v) if (length(v) <= 1L) list(v) else do.call(c, lapply(seq_along(v), function(i) lapply(rec(v[-i]), function(r) c(v[i], r))))
    rec(seq_len(n))
  }
  out <- list()
  for (pm in perms) out <- c(out, .permute_signs(t[pm]))
  out[!duplicated(vapply(out, function(s) paste(vapply(s, .k, ""), collapse = "|"), ""))]
}

# ---------------------------------------------------------------------------
# Solving by type

.dio_run <- function(info, ty, params = NULL, limit = Inf) {
  n <- info$n
  var <- info$vars
  if (ty == "linear") {
    if (is.null(params)) params <- lapply(seq_len(n) - 1L, function(i) .int_param(paste0("t_", i)))
    return(.dio_linear(info, params))
  }
  if (ty == "univariate") {
    rs <- Filter(function(r) r$type == "Integer", .unique_syms(roots(info$eq, var[[1L]])))
    return(lapply(rs, list))
  }
  if (ty == "binary_quadratic") {
    if (is.null(params)) params <- list(.int_param("t_0"), .int_param("t_1"))
    return(.dio_binary_quadratic(info, params))
  }
  if (ty %in% c("homogeneous_ternary_quadratic", "homogeneous_ternary_quadratic_normal")) {
    r <- .tern_solve(.tern_coeffs(info))
    return(if (is.null(r)) list() else list(lapply(r, .zs)))
  }
  if (ty == "general_sum_of_squares") {
    k <- -.dio_const(info)
    if (k < 0 || limit < 1) return(list())
    signs <- vapply(var, function(v) if (isTRUE(.facts(v)[["nonpositive"]])) -1L else 1L, 0L)
    out <- list()
    for (t in .power_representation(k, 2L, n, zeros = TRUE, limit = limit))
      out <- .dio_set_add(out, Map(function(s, j) .zs(s * j), signs, t))
    return(unname(out))
  }
  if (ty == "general_pythagorean") {
    cf <- lapply(seq_len(n), function(i) { ex <- integer(n); ex[i] <- 2L; .dio_c(info, ex) })
    if (sum(vapply(cf[1:3], .zsign, 0L)) < 0) cf <- lapply(cf, function(v) -v)
    if (is.null(params)) params <- lapply(seq_len(n - 1L) - 1L, function(i) .int_param(paste0("t_", i)))
    m <- params
    index <- 0L
    for (i in seq_len(n)) if (.zsign(cf[[i]]) == -1L) index <- i - 1L
    ith <- .add(lapply(m, function(mi) .pow(mi, S_Two)))
    L <- c(list(.sub(ith, .mul(list(S_Two, .pow(m[[n - 1L]], S_Two))))),
           lapply(seq_len(n - 2L), function(i) .mul(list(S_Two, m[[i]], m[[n - 1L]]))))
    sol <- c(L[seq_len(index)], list(ith), L[setdiff(seq_along(L), seq_len(index))])
    lcm <- gmp::as.bigz(1)
    for (i in seq_len(n)) {
      i0 <- i - 1L
      if (i0 == index || (index > 0L && i0 == 0L) || (index == 0L && i0 == 1L)) {
        lcm <- .zlcm(lcm, .isqrtz(.zabs(cf[[i]])))
      } else {
        s <- .isqrtz(cf[[i]])
        lcm <- .zlcm(lcm, if (s %% 2 != 0) s else s %/% 2)
      }
    }
    sol <- lapply(seq_len(n), function(i) expand(.div(.mul(list(.zs(lcm), sol[[i]])), .zs(.isqrtz(.zabs(cf[[i]]))))))
    return(list(sol))
  }
  if (ty == "general_sum_of_even_powers") {
    p <- info$total_degree
    nn <- -.dio_const(info)
    if (nn < 0 || limit < 1) return(list())
    signs <- vapply(var, function(v) if (isTRUE(.facts(v)[["nonpositive"]])) -1L else 1L, 0L)
    out <- list()
    for (t in .power_representation(nn, p, n, limit = limit))
      out <- .dio_set_add(out, Map(function(s, j) .zs(s * j), signs, t))
    return(unname(out))
  }
  stop("No solver has been written for ", ty, ".", call. = FALSE)
}

# classify_diop raises (sympy/solvers/diophantine/diophantine.py:1658)
.dio_unrecognized <- function() stop(.sympy_condition("NotImplementedError", filldedent("
        This equation is not yet recognized or else has not been
        simplified sufficiently to put it in a form recognized by
        diop_classify().")))

# diop_solve(): solutions of one irreducible equation in the parametrized
# form SymPy's diophantine() merges.
.dio_solve_eq <- function(eq, param = .int_param("t")) {
  info <- .dio_parse(eq)
  ty <- .dio_type(info)
  if (is.null(ty)) .dio_unrecognized()
  n <- info$n
  pn <- if (is_sym(param)) param$name else NULL
  switch(ty,
    linear = {
      params <- if (is.null(pn)) NULL else lapply(seq_len(n) - 1L, function(i) .int_param(paste0(pn, "_", i)))
      r <- .dio_run(info, ty, params)
      if (is.null(pn) && length(r)) r <- lapply(r, function(s) lapply(s, function(v) xreplace(v, lapply(free_symbols(v), function(z) list(z, S_Zero)))))
      if (length(r)) r[[1L]] else rep(list(NULL), n)
    },
    binary_quadratic = .dio_run(info, ty, list(if (is.null(pn)) .int_param("t") else param, .int_param("u"))),
    homogeneous_ternary_quadratic = , homogeneous_ternary_quadratic_normal = {
      .parametrize_ternary(.tern_solve(.tern_coeffs(info)), info)
    },
    general_pythagorean = {
      params <- if (is.null(pn)) NULL else lapply(seq_len(n - 1L), function(i) .int_param(paste0(pn, i)))
      .dio_run(info, ty, params)[[1L]]
    },
    univariate = .dio_run(info, ty),
    general_sum_of_squares = , general_sum_of_even_powers = .dio_run(info, ty, limit = Inf),
    stop("No solver has been written for ", ty, ".", call. = FALSE))
}

.merge_solution <- function(var, var_t, solution) {
  if (any(vapply(solution, is.null, TRUE))) return(NULL)
  tk <- vapply(var_t, .k, "")
  sol <- list()
  j <- 0L
  np <- 0L
  for (v in var) {
    if (v$key %in% tk) {
      j <- j + 1L
      sol[[length(sol) + 1L]] <- as_sym(solution[[j]])
    } else {
      np <- np + 1L
      sol[[length(sol) + 1L]] <- .int_param(paste0("n", np))
    }
  }
  for (i in seq_along(var)) if (identical(.fits_assumptions(sol[[i]], var[[i]]), FALSE)) return(NULL)
  sol
}

.diophantine_set <- function(eq, syms = NULL, param = .int_param("t"), permute = FALSE) {
  eq <- .eq_to_expr(eq)
  var <- free_symbols(expand(eq))
  if (!is.null(syms)) {
    syms <- Filter(function(s) any(vapply(var, function(v) identical(v$key, s$key), TRUE)), .as_gens(list(syms)))
    if (!identical(vapply(syms, .k, ""), vapply(var, .k, ""))) {
      # SymPy picks t[index of v in syms] for each v in var
      idx <- match(vapply(var, .k, ""), vapply(syms, .k, ""))
      return(lapply(.diophantine_set(eq, NULL, param, permute), function(t) t[idx]))
    }
  }
  nd <- fraction(together(eq))
  if (!length(free_symbols(nd[[1L]]))) return(list())
  if (length(free_symbols(nd[[2L]]))) {
    dsol <- .diophantine_set(nd[[2L]], NULL, param)
    dk <- vapply(dsol, function(s) paste(vapply(s, .k, ""), collapse = "|"), "")
    good <- Filter(function(s) !(paste(vapply(s, .k, ""), collapse = "|") %in% dk), .diophantine_set(nd[[1L]], NULL, param))
    return(Filter(function(s) !is_zero_num(expand(xreplace(nd[[2L]], Map(list, var, s)))), good))
  }
  e <- factor_terms(nd[[1L]])
  facs <- if (e$type == "Mul") e$args else list(e)
  e <- .mul(Filter(function(f) length(free_symbols(f)) > 0, facs))
  e <- expand(e)
  if (!is_polynomial(e)) stop("Equation should be a polynomial with Rational coefficients.", call. = FALSE)
  # clear rational coefficients, as SymPy's as_numer_denom() does
  qc <- .q_content(e)
  if (qc != 1) e <- expand(.div(e, make_rational(qc)))
  info <- .dio_parse(e)
  ty <- .dio_type(info)
  flags <- list(signs = FALSE, few = FALSE, var = FALSE)
  if (!is.null(ty) && permute) {
    if (ty %in% c("general_sum_of_squares", "general_sum_of_even_powers")) flags$var <- TRUE
    else if (ty %in% c("homogeneous_ternary_quadratic", "homogeneous_ternary_quadratic_normal", "binary_quadratic")) {
      nv <- info$n
      if (nv %in% 2:3) {
        prs <- utils::combn(nv, 2L, simplify = FALSE)
        xy_coeff <- all(vapply(prs, function(p) .dio_mono(info, p[1L], p[2L]) != 0, TRUE))
        x_coeff <- all(vapply(seq_len(nv), function(i) .dio_mono(info, i) != 0, TRUE))
        if (!xy_coeff && !x_coeff) flags$signs <- TRUE
        else if (!x_coeff) flags$few <- TRUE
      }
    }
  }
  if (!is.null(ty) && ty == "general_sum_of_squares") {
    terms <- list(e)
  } else {
    fl <- factor_list(e)
    if (!identical(fl$coeff$key, S_One$key))
      return(.diophantine_set(expand(.div(e, fl$coeff)), syms, param, permute))
    terms <- lapply(fl$factors, `[[`, 1L)
  }
  sols <- list()
  addsol <- function(s) {
    if (is.null(s)) return(invisible())
    sols[[paste(vapply(s, .k, ""), collapse = "|")]] <<- s
  }
  for (base in terms) {
    binfo <- .dio_parse(base)
    bty <- .dio_type(binfo)
    if (could_extract_minus_sign(base)) base <- .neg(base)
    sp <- as_coeff_mul_split(base)
    base <- sp[[2L]]
    solution <- .dio_solve_eq(base, param)
    if (is.null(bty)) stop("unhandled type", call. = FALSE)
    if (bty %in% c("linear", "homogeneous_ternary_quadratic", "homogeneous_ternary_quadratic_normal", "general_pythagorean")) {
      addsol(.merge_solution(var, binfo$vars, solution))
    } else if (bty %in% c("binary_quadratic", "general_sum_of_squares", "general_sum_of_even_powers", "univariate")) {
      for (s in solution) addsol(.merge_solution(var, binfo$vars, s))
    } else stop("unhandled type: ", bty, call. = FALSE)
  }
  null <- rep(list(S_Zero), length(var))
  if (!length(sols) && is_zero_num(expand(xreplace(e, Map(list, var, null))))) {
    if (all(vapply(var, function(v) !identical(.fits_assumptions(S_Zero, v), FALSE), TRUE))) addsol(null)
  }
  final <- list()
  addf <- function(s) final[[paste(vapply(s, .k, ""), collapse = "|")]] <<- s
  for (s in sols) {
    if (all(vapply(s, function(v) v$type == "Integer", TRUE))) {
      if (flags$signs) for (q in .permute_signs(s)) addf(q)
      else if (flags$few) {
        prod12 <- .mul(list(s[[1L]], s[[2L]]))$key
        for (q in .permute_signs(s)) if (identical(.mul(list(q[[1L]], q[[2L]]))$key, prod12)) addf(q)
      } else if (flags$var) for (q in .signed_permutations(s)) addf(q)
      else addf(s)
    } else addf(s)
  }
  unname(final)
}

#' Diophantine equations
#'
#' `diophantine(eq)` finds the integer solutions of a polynomial equation
#' with integer coefficients, as a `FiniteSet` of `Tuple`s ordered like the
#' sorted free symbols (or `syms`). Parametric families use integer
#' parameters such as `t_0`, `t_1`. Linear, binary quadratic (including
#' Pell equations), homogeneous ternary quadratic, general Pythagorean,
#' sums of squares and sums of even powers are supported, and products of
#' such factors. `permute = TRUE` adds the sign and order permutations
#' SymPy adds. `classify_diop(eq)` returns the variables, coefficients and
#' type; `diop_solve(eq)` solves one equation of a known type.
#'
#' @param eq An equation or an expression equal to zero.
#' @param param The base name of the parameters.
#' @param syms Optional order of the variables.
#' @param permute Add sign and order permutations of the solutions.
#' @return A set of tuples; `classify_diop()` returns a list.
#' @examples
#' x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
#' diophantine(2*x + 3*y - 5)
#' diophantine(x^2 - 2*y^2 - 1)
#' diophantine(x^2 + y^2 - z^2)
#' classify_diop(x^2 + y^2 - x*y + x + 5)
#' @export
diophantine <- function(eq, param = Symbol("t", integer = TRUE), syms = NULL, permute = FALSE) {
  sols <- .diophantine_set(eq, syms, if (is.null(param)) NULL else as_sym(param), permute)
  if (!length(sols)) return(S_EmptySet)
  FiniteSet(lapply(sols, function(s) do.call(Tuple, s)))
}

#' @rdname diophantine
#' @export
classify_diop <- function(eq) {
  info <- .dio_parse(.eq_to_expr(eq))
  ty <- .dio_type(info)
  if (is.null(ty)) .dio_unrecognized()
  mons <- lapply(seq_len(nrow(info$E)), function(i) .mul(Map(function(v, k) .pow(v, .int(k)), info$vars, info$E[i, ])))
  cf <- stats::setNames(lapply(info$c, .zs), vapply(mons, .str, ""))
  list(info$vars, cf, ty)
}

#' @rdname diophantine
#' @export
diop_solve <- function(eq, param = Symbol("t", integer = TRUE)) {
  r <- .dio_solve_eq(.eq_to_expr(eq), if (is.null(param)) NULL else as_sym(param))
  if (!length(r)) return(S_EmptySet)
  if (is.list(r[[1L]]) && !is_sym(r[[1L]])) FiniteSet(lapply(r, function(s) do.call(Tuple, s)))
  else if (any(vapply(r, is.null, TRUE))) r else do.call(Tuple, r)
}
