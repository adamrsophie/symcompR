# Wigner, Clebsch-Gordan, Racah and Gaunt coefficients and the Wigner D
# matrix, after sympy/physics/wigner.py. The coefficients are computed
# exactly with big integers (math.comb in SymPy, gmp::chooseZ here).

# Python's repr() of a float, for error messages.
.py_float_repr <- function(x) {
  if (is.na(x)) return("nan")
  if (is.infinite(x)) return(if (x > 0) "inf" else "-inf")
  d <- 1L
  while (d < 17L && as.numeric(sprintf("%.*e", d - 1L, x)) != x) d <- d + 1L
  s <- sprintf("%.*e", d - 1L, x)
  mant <- sub("e.*$", "", s)
  ex <- as.integer(sub("^.*e", "", s))
  digits <- gsub("[-.]", "", mant)
  neg <- startsWith(mant, "-")
  if (ex < -4L || ex >= 16L) {
    m <- if (nchar(digits) > 1L) paste0(substr(digits, 1L, 1L), ".", substring(digits, 2L)) else digits
    return(paste0(if (neg) "-", m, "e", if (ex < 0) "-" else "+", sprintf("%02d", abs(ex))))
  }
  if (ex >= 0L) {
    if (nchar(digits) <= ex + 1L) {
      r <- paste0(digits, strrep("0", ex + 1L - nchar(digits)), ".0")
    } else {
      r <- paste0(substr(digits, 1L, ex + 1L), ".", substring(digits, ex + 2L))
    }
  } else {
    r <- paste0("0.", strrep("0", -ex - 1L), digits)
  }
  paste0(if (neg) "-", r)
}

# str() of a Python number or SymPy object, for error messages.
.py_value_str <- function(v) {
  if (is_sym(v)) return(sstr(v))
  if (is.double(v) && length(v) == 1L) return(.py_float_repr(v))
  as.character(v)
}

# as_int(value, strict=False)
.wig_as_int <- function(v) {
  bad <- function() .physics_error("ValueError", paste0(.py_value_str(v), " is not an integer"))
  if (is_sym(v)) {
    if (v$type == "Integer") return(as.numeric(v$value))
    if (v$type == "Rational") bad()
    if (v$type == "Float") {
      x <- as.numeric(v)
      if (x != trunc(x)) bad()
      return(trunc(x))
    }
    bad()
  }
  if (is.logical(v)) return(as.numeric(v))
  if (!is.numeric(v) || length(v) != 1L) bad()
  if (v != trunc(v)) bad()
  as.numeric(v)
}

# twice a value that must be an integer or a half integer
.wig_doubled_int <- function(v) {
  w <- if (is_sym(v)) .mul(list(v, S_Two)) else v * 2
  r <- tryCatch(.wig_as_int(w), ValueError = function(e) NULL)
  if (is.null(r)) .physics_error("ValueError", paste0("expecting integer or half-integer, got ", .py_value_str(v)))
  r
}

.wig_check_couple <- function(dj1, dj2, dj3) {
  if ((dj1 + dj2 + dj3) %% 2 != 0) return(FALSE)
  if (dj1 > dj2 + dj3) return(FALSE)
  if (dj2 > dj1 + dj3) return(FALSE)
  if (dj3 > dj1 + dj2) return(FALSE)
  TRUE
}

.wig_comb <- function(n, k) {
  if (k < 0 || n < 0) .physics_error("ValueError", "k must be a non-negative integer")
  gmp::chooseZ(n, k)
}

.wig_sign <- function(k) if (k %% 2 == 0) 1 else -1

.wig_rat <- function(num, den) make_rational(gmp::as.bigq(num, den))

.wig_prec <- function(res, prec) if (!is.null(prec) && !isFALSE(prec) && !identical(prec, 0)) N(res, prec) else res

#' Wigner 3j, 6j and 9j symbols and related coefficients
#'
#' Exact values of the Wigner 3j, 6j and 9j symbols, the Clebsch-Gordan
#' coefficients, the Racah symbol, the Gaunt coefficient and its real
#' counterpart, after `sympy.physics.wigner`. Arguments are integers or
#' half integers (numbers or exact rationals); with `prec` the result is
#' evaluated numerically to that many digits.
#'
#' @param j_1,j_2,j_3,j_4,j_5,j_6,j_7,j_8,j_9 Angular momenta.
#' @param m_1,m_2,m_3 Projections.
#' @param aa,bb,cc,dd,ee,ff Arguments of the Racah symbol.
#' @param l_1,l_2,l_3 Degrees of the spherical harmonics.
#' @param mu_1,mu_2,mu_3 Orders of the real spherical harmonics.
#' @param prec Optional number of digits for a numerical result.
#' @return A `sym` number.
#' @examples
#' wigner_3j(2, 6, 4, 0, 0, 0)
#' clebsch_gordan(Rational(3, 2), S_Half, 2, Rational(3, 2), S_Half, 2)
#' wigner_6j(3, 3, 3, 3, 3, 3)
#' wigner_9j(1, 1, 1, 1, 1, 1, 1, 1, 0, prec = 20)
#' racah(3, 3, 3, 3, 3, 3)
#' gaunt(1, 0, 1, 1, 0, -1)
#' real_gaunt(2, 2, 2, 0, 2, 2)
#' @export
wigner_3j <- function(j_1, j_2, j_3, m_1, m_2, m_3, prec = NULL) {
  d <- vapply(list(j_1, j_2, j_3, m_1, m_2, m_3), .wig_doubled_int, 0)
  dj1 <- d[1]; dj2 <- d[2]; dj3 <- d[3]; dm1 <- d[4]; dm2 <- d[5]; dm3 <- d[6]
  if (dm1 + dm2 + dm3 != 0) return(S_Zero)
  if (!.wig_check_couple(dj1, dj2, dj3)) return(S_Zero)
  if (abs(dm1) > dj1 || abs(dm2) > dj2 || abs(dm3) > dj3) return(S_Zero)
  if (!((dj1 - dm1) %% 2 == 0 && (dj2 - dm2) %% 2 == 0 && (dj3 - dm3) %% 2 == 0)) return(S_Zero)
  sumj <- (dj1 + dj2 + dj3) %/% 2
  jm1 <- sumj - dj1; jm2 <- sumj - dj2; jm3 <- sumj - dj3
  j1mm1 <- (dj1 - dm1) %/% 2; j2mm2 <- (dj2 - dm2) %/% 2; j3mm3 <- (dj3 - dm3) %/% 2
  j1pm1 <- (dj1 + dm1) %/% 2
  imin <- max(0, j1pm1 - jm2, j2mm2 - jm1)
  imax <- min(jm3, j1pm1, j2mm2)
  sumres <- gmp::as.bigz(0)
  for (ii in seq_len(max(0, imax - imin + 1)) + imin - 1) {
    ti <- .wig_comb(jm3, ii) * .wig_comb(jm2, j1pm1 - ii) * .wig_comb(jm1, j2mm2 - ii)
    sumres <- ti - sumres
  }
  if (sumres == 0) return(S_Zero)
  res_num <- .wig_comb(dj1, jm2) * .wig_comb(dj2, jm1)
  res_den <- .wig_comb(sumj, jm3) * .wig_comb(dj1, j1mm1) * .wig_comb(dj2, j2mm2) * .wig_comb(dj3, j3mm3) * (sumj + 1)
  ressqrt <- .sqrt(.wig_rat(res_num, res_den))
  phase <- .wig_sign(dj1 + (dj3 + dm3) %/% 2 + imax)
  res <- .mul(list(.int(phase), ressqrt, make_rational(gmp::as.bigq(sumres))))
  .wig_prec(res, prec)
}

#' @rdname wigner_3j
#' @export
clebsch_gordan <- function(j_1, j_2, j_3, m_1, m_2, m_3, prec = NULL) {
  w <- wigner_3j(j_1, j_2, j_3, m_1, m_2, .neg_any(m_3), prec = prec)
  ex <- .add(list(as_sym(j_1), .neg(as_sym(j_2)), as_sym(m_3)))
  .mul(list(.pow(S_NegativeOne, ex), .sqrt(Integer(.wig_as_int(.add(list(.mul(list(S_Two, as_sym(j_3))), S_One))))), w))
}

.neg_any <- function(v) if (is_sym(v)) .neg(v) else -v

#' @rdname wigner_3j
#' @export
wigner_6j <- function(j_1, j_2, j_3, j_4, j_5, j_6, prec = NULL) {
  d <- vapply(list(j_1, j_2, j_3, j_4, j_5, j_6), .wig_doubled_int, 0)
  dj1 <- d[1]; dj2 <- d[2]; dj3 <- d[3]; dj4 <- d[4]; dj5 <- d[5]; dj6 <- d[6]
  if (!.wig_check_couple(dj1, dj2, dj3)) return(S_Zero)
  if (!.wig_check_couple(dj1, dj5, dj6)) return(S_Zero)
  if (!.wig_check_couple(dj4, dj2, dj6)) return(S_Zero)
  if (!.wig_check_couple(dj4, dj5, dj3)) return(S_Zero)
  j123 <- (dj1 + dj2 + dj3) %/% 2; j156 <- (dj1 + dj5 + dj6) %/% 2
  j426 <- (dj4 + dj2 + dj6) %/% 2; j453 <- (dj4 + dj5 + dj3) %/% 2
  jpm123 <- (dj1 + dj2 - dj3) %/% 2; jpm132 <- (dj1 + dj3 - dj2) %/% 2
  jpm231 <- (dj2 + dj3 - dj1) %/% 2; jpm156 <- (dj1 + dj5 - dj6) %/% 2
  jpm426 <- (dj4 + dj2 - dj6) %/% 2; jpm453 <- (dj4 + dj5 - dj3) %/% 2
  sumres <- gmp::as.bigz(0)
  imin <- max(j123, j453, j426, j156)
  imax <- max(jpm123 + j453, jpm132 + j426, jpm231 + j156)
  for (ii in seq_len(max(0, imax - imin + 1)) + imin - 1) {
    ti <- .wig_comb(ii + 1, j123 + 1) * .wig_comb(jpm123, ii - j453) *
      .wig_comb(jpm132, ii - j426) * .wig_comb(jpm231, ii - j156)
    sumres <- ti - sumres
  }
  if (sumres == 0) return(S_Zero)
  if (imax %% 2 == 1) sumres <- -sumres
  res_num <- .wig_comb(j123 + 1, dj1 + 1) * .wig_comb(dj1, jpm123)
  res_den <- .wig_comb(j156 + 1, dj1 + 1) * .wig_comb(dj1, jpm156) *
    .wig_comb(j426 + 1, dj4 + 1) * .wig_comb(dj4, jpm426) *
    .wig_comb(j453 + 1, dj4 + 1) * .wig_comb(dj4, jpm453)
  res_den <- res_den * gmp::as.bigz(dj4 + 1)^2
  ressqrt <- .sqrt(.wig_rat(res_num, res_den))
  result <- .mul(list(ressqrt, make_rational(gmp::as.bigq(sumres))))
  .wig_prec(result, prec)
}

#' @rdname wigner_3j
#' @export
racah <- function(aa, bb, cc, dd, ee, ff, prec = NULL) {
  res <- wigner_6j(aa, bb, ee, dd, cc, ff)
  if (identical(res$key, S_Zero$key)) return(S_Zero)
  tot <- .add(lapply(list(aa, bb, cc, dd), as_sym))
  k <- trunc(as.numeric(tot))
  res <- .mul(list(.int(.wig_sign(k)), res))
  .wig_prec(res, prec)
}

#' @rdname wigner_3j
#' @export
wigner_9j <- function(j_1, j_2, j_3, j_4, j_5, j_6, j_7, j_8, j_9, prec = NULL) {
  d <- vapply(list(j_1, j_2, j_3, j_4, j_5, j_6, j_7, j_8, j_9), .wig_doubled_int, 0)
  dj1 <- d[1]; dj2 <- d[2]; dj3 <- d[3]; dj4 <- d[4]; dj5 <- d[5]; dj6 <- d[6]
  dj7 <- d[7]; dj8 <- d[8]; dj9 <- d[9]
  if (!.wig_check_couple(dj1, dj2, dj3)) return(S_Zero)
  if (!.wig_check_couple(dj4, dj5, dj6)) return(S_Zero)
  if (!.wig_check_couple(dj7, dj8, dj9)) return(S_Zero)
  if (!.wig_check_couple(dj1, dj4, dj7)) return(S_Zero)
  if (!.wig_check_couple(dj2, dj5, dj8)) return(S_Zero)
  if (!.wig_check_couple(dj3, dj6, dj9)) return(S_Zero)
  j123 <- (dj1 + dj2 + dj3) %/% 2; j456 <- (dj4 + dj5 + dj6) %/% 2
  j789 <- (dj7 + dj8 + dj9) %/% 2; j147 <- (dj1 + dj4 + dj7) %/% 2
  j258 <- (dj2 + dj5 + dj8) %/% 2; j369 <- (dj3 + dj6 + dj9) %/% 2
  pm123 <- (dj1 + dj2 - dj3) %/% 2; pm132 <- (dj1 + dj3 - dj2) %/% 2; pm231 <- (dj2 + dj3 - dj1) %/% 2
  pm456 <- (dj4 + dj5 - dj6) %/% 2; pm465 <- (dj4 + dj6 - dj5) %/% 2; pm564 <- (dj5 + dj6 - dj4) %/% 2
  pm789 <- (dj7 + dj8 - dj9) %/% 2; pm798 <- (dj7 + dj9 - dj8) %/% 2; pm897 <- (dj8 + dj9 - dj7) %/% 2
  sumres <- gmp::as.bigz(0)
  dtmin <- max(abs(dj2 - dj6), abs(dj4 - dj8), abs(dj1 - dj9))
  dtmax <- min(dj2 + dj6, dj4 + dj8, dj1 + dj9)
  rng <- function(a, b) if (b >= a) a:b else numeric()
  for (dt in if (dtmax >= dtmin) seq(dtmin, dtmax, by = 2) else numeric()) {
    j19t <- (dj1 + dj9 + dt) %/% 2; j26t <- (dj2 + dj6 + dt) %/% 2; j48t <- (dj4 + dj8 + dt) %/% 2
    abc <- gmp::as.bigz(dt + 1)
    xmin <- max(j123, j369, j26t, j19t); xmax <- min(pm123 + j369, pm132 + j26t, pm231 + j19t)
    suma <- gmp::as.bigz(0)
    for (x in rng(xmin, xmax)) {
      ax <- .wig_comb(x + 1, j19t + 1) * .wig_comb((dj1 + dj9 - dt) %/% 2, x - j26t) *
        .wig_comb((dj1 + dt - dj9) %/% 2, x - j369) * .wig_comb((dt + dj9 - dj1) %/% 2, x - j123)
      suma <- ax - suma
    }
    ymin <- max(j456, j258, j26t, j48t); ymax <- min(pm456 + j26t, pm465 + j258, pm564 + j48t)
    sumb <- gmp::as.bigz(0)
    for (y in rng(ymin, ymax)) {
      by <- .wig_comb(y + 1, j26t + 1) * .wig_comb((dj2 + dj6 - dt) %/% 2, y - j48t) *
        .wig_comb((dt + dj6 - dj2) %/% 2, y - j258) * .wig_comb((dj2 + dt - dj6) %/% 2, y - j456)
      sumb <- by - sumb
    }
    zmin <- max(j789, j147, j48t, j19t); zmax <- min(pm789 + j19t, pm798 + j48t, pm897 + j147)
    sumc <- gmp::as.bigz(0)
    for (z in rng(zmin, zmax)) {
      cz <- .wig_comb(z + 1, j48t + 1) * .wig_comb((dj4 + dj8 - dt) %/% 2, z - j19t) *
        .wig_comb((dt + dj8 - dj4) %/% 2, z - j147) * .wig_comb((dj4 + dt - dj8) %/% 2, z - j789)
      sumc <- cz - sumc
    }
    abc <- abc * suma * sumb * sumc
    if ((xmax + ymax + zmax) %% 2 == 1) abc <- -abc
    sumres <- sumres + abc
  }
  if (sumres == 0) return(S_Zero)
  if (dtmax %% 2 == 1) sumres <- -sumres
  res_den <- gmp::as.bigz((dj3 + 1) * (dj6 + 1) * (dj9 + 1) * (dj7 + 1) * (dj8 + 1) * (dj9 + 1))
  res_den <- res_den * .wig_comb(j123 + 1, dj3 + 1) * .wig_comb(dj3, pm231) *
    .wig_comb(j456 + 1, dj6 + 1) * .wig_comb(dj6, pm564) *
    .wig_comb(j789 + 1, dj9 + 1) * .wig_comb(dj9, pm897) *
    .wig_comb(j147 + 1, dj7 + 1) * .wig_comb(dj7, (dj4 + dj7 - dj1) %/% 2) *
    .wig_comb(j258 + 1, dj8 + 1) * .wig_comb(dj8, (dj5 + dj8 - dj2) %/% 2) *
    .wig_comb(j369 + 1, dj9 + 1) * .wig_comb(dj9, (dj6 + dj9 - dj3) %/% 2)
  res <- .mul(list(make_rational(gmp::as.bigq(sumres)), .sqrt(.wig_rat(1, res_den))))
  .wig_prec(res, prec)
}

#' @rdname wigner_3j
#' @export
gaunt <- function(l_1, l_2, l_3, m_1, m_2, m_3, prec = NULL) {
  d <- vapply(list(l_1, l_2, l_3, m_1, m_2, m_3), .wig_as_int, 0)
  l1 <- d[1]; l2 <- d[2]; l3 <- d[3]; m1 <- d[4]; m2 <- d[5]; m3 <- d[6]
  if (!.wig_check_couple(l1, l2, l3)) return(S_Zero)
  if (m1 + m2 + m3 != 0) return(S_Zero)
  if (abs(m1) > l1 || abs(m2) > l2 || abs(m3) > l3) return(S_Zero)
  sumj <- l1 + l2 + l3
  g <- sumj %/% 2
  jm1 <- l2 + l3 - l1; jm2 <- l1 + l3 - l2; jm3 <- l1 + l2 - l3
  j1mm1 <- l1 - m1; j2mm2 <- l2 - m2; j3mm3 <- l3 - m3; j1pm1 <- l1 + m1
  imin <- max(0, j1pm1 - jm2, j2mm2 - jm1)
  imax <- min(jm3, j1pm1, j2mm2)
  sumres <- gmp::as.bigz(0)
  for (ii in seq_len(max(0, imax - imin + 1)) + imin - 1) {
    ti <- .wig_comb(jm3, ii) * .wig_comb(jm2, j1pm1 - ii) * .wig_comb(jm1, j2mm2 - ii)
    sumres <- ti - sumres
  }
  if (sumres == 0) return(S_Zero)
  if ((g + j3mm3 + imax) %% 2 == 1) sumres <- -sumres
  sumres <- sumres * .wig_comb(g, g - l3) * .wig_comb(l3, g - l1)
  res_num <- (2 * l3 + 1) * .wig_comb(2 * l3, jm1)
  res_den <- .wig_comb(2 * l1, j1mm1) * .wig_comb(2 * l2, j2mm2) * .wig_comb(2 * l3, j3mm3) *
    .wig_comb(sumj + 1, 2 * l1 + 1) * .wig_comb(sumj + 1, 2 * l2 + 1)
  ressqrt <- .sqrt(.wig_rat(res_num, res_den))
  res <- .div(.mul(list(ressqrt, make_rational(gmp::as.bigq(sumres)))), .sqrt(.mul(list(.int(4), S_Pi))))
  if (!is.null(prec)) res <- N(res, prec)
  res
}

#' @rdname wigner_3j
#' @export
real_gaunt <- function(l_1, l_2, l_3, mu_1, mu_2, mu_3, prec = NULL) {
  d <- vapply(list(l_1, l_2, l_3, mu_1, mu_2, mu_3), .wig_as_int, 0)
  l_1 <- d[1]; l_2 <- d[2]; l_3 <- d[3]; mu_1 <- d[4]; mu_2 <- d[5]; mu_3 <- d[6]
  if (sum(c(mu_1, mu_2, mu_3) < 0) %% 2) return(S_Zero)
  if ((l_1 + l_2 + l_3) %% 2) return(S_Zero)
  lmax <- l_2 + l_3
  lmin <- max(abs(l_2 - l_3), min(abs(mu_2 + mu_3), abs(mu_2 - mu_3)))
  if ((lmin + lmax) %% 2) lmin <- lmin + 1
  if (lmin > lmax) return(S_Zero)
  kron_del <- function(i, j) if (i == j) 1 else 0
  s <- function(e) if (e %% 2) -1 else 1
  t <- function(x) if (x > 0) 1 else 0
  A <- function(mu, m) t(-mu) * (s(m) * kron_del(m, mu) - kron_del(m, -mu))
  B <- function(mu, m) t(mu) * (kron_del(m, mu) + s(m) * kron_del(m, -mu))
  U <- function(mu, m) .mul(list(.int(kron_del(abs(mu), abs(m))),
                                .add(list(.int(kron_del(mu, 0) * kron_del(m, 0)),
                                          .div(.add(list(.int(B(mu, m)), .mul(list(S_I, .int(A(mu, m)))))), .sqrt(S_Two))))))
  ugnt <- S_Zero
  for (m1 in seq(-l_1, l_1)) {
    U1 <- U(mu_1, m1)
    for (m2 in seq(-l_2, l_2)) {
      U2 <- U(mu_2, m2)
      U3 <- U(mu_3, -m1 - m2)
      ugnt <- .add(list(ugnt, .mul(list(fcall("re", list(.mul(list(U1, U2, U3)))),
                                        gaunt(l_1, l_2, l_3, m1, m2, -m1 - m2, prec = prec)))))
    }
  }
  ugnt
}

# Wigner3j: an unevaluated 3j symbol whose doit() evaluates numeric arguments
def_function("Wigner3j", nargs = 6L)
local({
  old <- .doit_hooks$Function
  register_doit("Function", function(e, deep = TRUE) {
    if (identical(e$name, "Wigner3j")) {
      if (all(vapply(e$args, function(a) .wig_is_number(a), TRUE))) return(do.call(wigner_3j, e$args))
      return(e)
    }
    if (is.null(old)) e else if ("deep" %in% names(formals(old))) old(e, deep = deep) else old(e)
  })
})

.wig_is_number <- function(a) {
  if (is_number_type(a)) return(TRUE)
  !length(free_symbols(a)) && a$type %in% c("Add", "Mul", "Pow", "Constant", "Function")
}

#' Dot product of rotational gradients of spherical harmonics
#'
#' `dot_rot_grad_Ynm(j, p, l, m, theta, phi)` returns the dot product of the
#' rotational gradients of two spherical harmonics as an unevaluated sum
#' over Wigner 3j symbols; [doit()] evaluates it. After
#' `sympy.physics.wigner.dot_rot_grad_Ynm`.
#'
#' @param j,p Degree and order of the first spherical harmonic.
#' @param l,m Degree and order of the second spherical harmonic.
#' @param theta,phi The angles.
#' @return A `sym` expression.
#' @examples
#' doit(dot_rot_grad_Ynm(3, 2, 2, 0, Symbol("theta"), Symbol("phi")))
#' @export
dot_rot_grad_Ynm <- function(j, p, l, m, theta, phi) {
  j <- as_sym(j); p <- as_sym(p); l <- as_sym(l); m <- as_sym(m)
  theta <- as_sym(theta); phi <- as_sym(phi)
  k <- Dummy("k")
  alpha <- function(l, m, j, p, k) {
    .mul(list(.sqrt(.div(.mul(list(.add(list(.mul(list(S_Two, l)), S_One)), .add(list(.mul(list(S_Two, j)), S_One)),
                                   .add(list(.mul(list(S_Two, k)), S_One)))), .mul(list(.int(4), S_Pi)))),
              fcall("Wigner3j", list(j, l, k, S_Zero, S_Zero, S_Zero)),
              fcall("Wigner3j", list(j, l, k, p, m, .sub(.neg(m), p)))))
  }
  summand <- .mul(list(.div(.mul(list(fcall("Ynm", list(k, .add(list(m, p)), theta, phi)), alpha(l, m, j, p, k))), S_Two),
                       .add(list(.pow(k, S_Two), .neg(.pow(j, S_Two)), .neg(.pow(l, S_Two)), k, .neg(j), .neg(l)))))
  .mul(list(.pow(S_NegativeOne, .add(list(m, p))), Sum(summand, list(k, fcall("Abs", list(.sub(l, j))), .add(list(l, j))))))
}

#' Wigner small d and D matrices
#'
#' `wigner_d_small(J, beta)` is the small Wigner d matrix of a rotation by
#' `beta` about the y axis for angular momentum `J`; `wigner_d(J, alpha,
#' beta, gamma)` the Wigner D matrix of the rotation with Euler angles
#' `alpha`, `beta`, `gamma`. After `sympy.physics.wigner`.
#'
#' @param J An integer or half integer.
#' @param alpha,beta,gamma Euler angles.
#' @return A `symmatrix`.
#' @examples
#' wigner_d_small(S_Half, S_Pi / 2)
#' wigner_d(S_Half, Symbol("a"), Symbol("b"), Symbol("c"))
#' @export
wigner_d_small <- function(J, beta) {
  J <- as_sym(J); beta <- as_sym(beta)
  n <- .wig_as_int(.add(list(.mul(list(S_Two, J)), S_One)))
  M <- lapply(seq_len(n) - 1L, function(i) .sub(J, .int(i)))
  d <- zeros(n)
  for (i in seq_len(n)) {
    Mi <- M[[i]]
    for (j in seq_len(n)) {
      Mj <- M[[j]]
      sigmamax <- as.numeric(Min(.sub(J, Mi), .sub(J, Mj)))
      sigmamin <- as.numeric(Max(S_Zero, .sub(.neg(Mi), Mj)))
      dij <- .sqrt(.div(fcall("binomial", list(.mul(list(S_Two, J)), .add(list(J, Mj)))),
                        fcall("binomial", list(.mul(list(S_Two, J)), .add(list(J, Mi))))))
      terms <- lapply(if (sigmamax >= sigmamin) sigmamin:sigmamax else numeric(), function(s) {
        s <- .int(s)
        .mul(list(.pow(S_NegativeOne, .sub(.sub(J, Mi), s)),
                  fcall("binomial", list(.add(list(J, Mj)), .sub(.sub(J, Mi), s))),
                  fcall("binomial", list(.sub(J, Mj), s)),
                  .pow(fcall("cos", list(.div(beta, S_Two))), .add(list(.mul(list(S_Two, s)), Mi, Mj))),
                  .pow(fcall("sin", list(.div(beta, S_Two))), .sub(.sub(.sub(.mul(list(S_Two, J)), .mul(list(S_Two, s))), Mj), Mi))))
      })
      d$data[[(i - 1L) * n + j]] <- .mul(list(dij, .add(terms)))
    }
  }
  d
}

#' @rdname wigner_d_small
#' @export
wigner_d <- function(J, alpha, beta, gamma) {
  J <- as_sym(J); alpha <- as_sym(alpha); gamma <- as_sym(gamma)
  d <- wigner_d_small(J, beta)
  n <- d$nrow
  M <- lapply(seq_len(n) - 1L, function(i) .sub(J, .int(i)))
  D <- zeros(n)
  for (i in seq_len(n)) for (j in seq_len(n))
    D$data[[(i - 1L) * n + j]] <- .mul(list(fcall("exp", list(.mul(list(S_I, M[[i]], alpha)))), .mget(d, i, j),
                                            fcall("exp", list(.mul(list(S_I, M[[j]], gamma))))))
  D
}
