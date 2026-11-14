# nsolve through mpmath's findroot (mpmath/calculus/optimization.py): the
# Secant solver for one equation and MDNewton for systems, run at the
# working precision prec + 20 bits, with findroot's default tolerance and
# verification. The root is rounded to the requested precision, as
# sympify() of an mpf.

.fs_mpc_eval <- function(e, env, wp) {
  # ponytail: a value below 2^-(2^20) enters as 0, since the function
  # evaluators allocate on the size of the exponent; exact where it only
  # meets sums of normal size
  env <- lapply(env, function(z) .mpc(.fs_flush(z$re), .fs_flush(z$im)))
  v <- .mpeval(e, wp, env)
  if (is.null(v)) stop(.sympy_condition("TypeError", "cannot evaluate the function numerically"))
  .mpc(.mp_round(v$re$m, v$re$e, wp), .mp_round(v$im$m, v$im$e, wp))
}
.fs_flush <- function(a) if (a$m != 0 && .mp_mag(a) < -2^20) .mp_zero else a
# an mpc value passed in the env of .mpeval as the value of a symbol
.fs_is_mpc <- function(v) is.list(v) && !is_sym(v) && identical(names(v), c("re", "im"))
.fs_mpc_zero <- function(z) z$re$m == 0 && z$im$m == 0
# normalized mpf values are equal when mantissa and exponent are
.fs_mp_eq <- function(a, b) a$m == b$m && (a$m == 0 || a$e == b$e)
.fs_mpc_eq <- function(a, b) .fs_mp_eq(a$re, b$re) && .fs_mp_eq(a$im, b$im)
# the sign of a - b; rounding to a few bits never changes it, and this
# avoids .mp_cmp's exact subtraction for exponents far apart
.fs_mp_cmp <- function(a, b) .mp_sign(.mp_sub(a, b, 8L))
.fs_mp_lt <- function(a, b) .fs_mp_cmp(a, b) < 0
.fs_mp_str <- function(a, wp) .mp_to_str(a, .prec_to_dps(wp), strip_zeros = TRUE)

.fs_mpc_from <- function(v, wp) {
  v <- as_sym(v)
  if (is_exact_number(v)) return(.mpc(.mp_from_q(v$value, wp)))
  r <- .mpeval(v, wp, list())
  if (is.null(r)) stop(.sympy_condition("TypeError", paste0("cannot create mpf from ", .str(v))))
  r
}

.fs_findroot_fail <- function(res2, tol, wp)
  stop(.sympy_condition("ValueError", sprintf(
    "Could not find root within given tolerance. (%s > %s)\nTry another starting point or tweak arguments.",
    .fs_mp_str(res2, wp), .fs_mp_str(tol, wp))))

.fs_nsolve1 <- function(fe, x, x0, prec, tol, maxsteps) {
  p <- .dps_to_prec(prec)
  wp <- p + 20L
  tol <- if (is.null(tol)) .mpf(1, 11L - wp) else .fs_tol_mp(tol, wp)
  f <- function(z) .fs_mpc_eval(fe, stats::setNames(list(z), x$name), wp)
  a0 <- .fs_mpc_from(x0, wp)
  fx0 <- f(a0)
  if (.fs_mpc_zero(fx0)) return(.mpc_sym(.mpc(.mp_round(a0$re$m, a0$re$e, p), .mp_round(a0$im$m, a0$im$e, p)), p))
  if (is.null(maxsteps)) maxsteps <- 30L
  # Secant: x1 = x0 + 0.25
  x0m <- a0
  x1m <- .mpc_add(a0, .mpc(.mpf(1, -2L)), wp)
  f0 <- fx0
  i <- 0L
  xr <- NULL
  repeat {
    f1 <- f(x1m)
    l <- .mpc_sub(x1m, x0m, wp)
    if (.fs_mpc_zero(l)) break
    s <- .mpc_div(.mpc_sub(f1, f0, wp), l, wp)
    if (.fs_mpc_zero(s)) break
    x0m <- x1m
    x1m <- .mpc_sub(x1m, .mpc_div(f1, s, wp), wp)
    f0 <- f1
    xr <- x1m
    err <- .mpc_abs(l, wp)
    i <- i + 1L
    nx <- .mpc_abs(xr, wp)
    bound <- .mp_mul(tol, if (.fs_mp_lt(nx, .mp_one)) .mp_one else nx, wp)
    if (.fs_mp_lt(err, bound) || i >= maxsteps) break
  }
  if (i == 0L)
    stop(.sympy_condition("ValueError", "Could not find root using the given solver.\nTry another starting point or tweak arguments."))
  r <- .mpc_abs(f(xr), wp)
  r2 <- .mp_mul(r, r, wp)
  if (.fs_mp_cmp(r2, tol) > 0) .fs_findroot_fail(r2, tol, wp)
  .mpc_sym(.mpc(.mp_round(xr$re$m, xr$re$e, p), .mp_round(xr$im$m, xr$im$e, p)), p)
}

.fs_tol_mp <- function(tol, wp) .mpc_abs(.fs_mpc_from(tol, wp), wp)

# mpmath's lu_solve at wp + 10 bits: LU with partial pivoting by modulus
.fs_lu_solve <- function(A, b, wp) {
  w <- wp + 10L
  n <- length(b)
  for (k in seq_len(n)) {
    mods <- vapply(k:n, function(i) .mp_to_dbl(.mpc_abs(A[[i]][[k]], 53L)), 0)
    pk <- k - 1L + which.max(mods)
    if (mods[pk - k + 1L] == 0) stop(.sympy_condition("ZeroDivisionError", "matrix is numerically singular"))
    if (pk != k) {
      tA <- A[[k]]; A[[k]] <- A[[pk]]; A[[pk]] <- tA
      tb <- b[[k]]; b[[k]] <- b[[pk]]; b[[pk]] <- tb
    }
    if (k < n) for (i in (k + 1L):n) {
      fct <- .mpc_div(A[[i]][[k]], A[[k]][[k]], w)
      for (j in k:n) A[[i]][[j]] <- .mpc_sub(A[[i]][[j]], .mpc_mul(fct, A[[k]][[j]], w), w)
      b[[i]] <- .mpc_sub(b[[i]], .mpc_mul(fct, b[[k]], w), w)
    }
  }
  xs <- vector("list", n)
  for (i in rev(seq_len(n))) {
    s <- b[[i]]
    if (i < n) for (j in (i + 1L):n) s <- .mpc_sub(s, .mpc_mul(A[[i]][[j]], xs[[j]], w), w)
    xs[[i]] <- .mpc_div(s, A[[i]][[i]], w)
  }
  lapply(xs, function(z) .mpc(.mp_round(z$re$m, z$re$e, wp), .mp_round(z$im$m, z$im$e, wp)))
}

.fs_nsolve_system <- function(fs, xs, x0, prec, tol, maxsteps) {
  p <- .dps_to_prec(prec)
  wp <- p + 20L
  tol <- if (is.null(tol)) .mpf(1, 11L - wp) else .fs_tol_mp(tol, wp)
  nm <- vapply(xs, function(v) v$name, "")
  J <- lapply(fs, function(fi) lapply(xs, function(v) .diff1(fi, v)))
  envof <- function(xv) stats::setNames(xv, nm)
  fvec <- function(xv) { env <- envof(xv); lapply(fs, function(fi) .fs_mpc_eval(fi, env, wp)) }
  jmat <- function(xv) { env <- envof(xv); lapply(J, function(row) lapply(row, function(d) .fs_mpc_eval(d, env, wp))) }
  norm <- function(v) Reduce(function(a, b) if (.fs_mp_lt(a, b)) b else a, lapply(v, .mpc_abs, p = wp))
  x0m <- lapply(x0, .fs_mpc_from, wp = wp)
  fx <- fvec(x0m)
  fxnorm <- norm(fx)
  fin <- function(xv) stats::setNames(lapply(xv, function(z) .mpc_sym(.mpc(.mp_round(z$re$m, z$re$e, p), .mp_round(z$im$m, z$im$e, p)), p)), nm)
  if (fxnorm$m == 0) return(fin(x0m))
  if (is.null(maxsteps)) maxsteps <- 10L
  i <- 0L
  cancel <- FALSE
  while (!cancel) {
    s <- .fs_lu_solve(jmat(x0m), lapply(fx, .mpc_neg), wp)
    l <- .mpc(.mp_one)
    x1 <- Map(function(a, d) .mpc_add(a, d, wp), x0m, s)
    repeat {
      if (all(unlist(Map(.fs_mpc_eq, x1, x0m)))) {
        cancel <- TRUE
        break
      }
      fx <- fvec(x1)
      newnorm <- norm(fx)
      if (.fs_mp_lt(newnorm, fxnorm)) {
        fxnorm <- newnorm
        x0m <- x1
        break
      }
      l <- .mpc_mul(l, .mpc(.mpf(1, -1L)), wp)
      x1 <- Map(function(a, d) .mpc_add(a, .mpc_mul(l, d, wp), wp), x0m, s)
    }
    i <- i + 1L
    nx <- norm(x0m)
    bound <- .mp_mul(tol, if (.fs_mp_lt(nx, .mp_one)) .mp_one else nx, wp)
    if (.fs_mp_lt(fxnorm, bound) || i >= maxsteps) break
  }
  r <- norm(fvec(x0m))
  r2 <- .mp_mul(r, r, wp)
  if (.fs_mp_cmp(r2, tol) > 0) .fs_findroot_fail(r2, tol, wp)
  fin(x0m)
}
