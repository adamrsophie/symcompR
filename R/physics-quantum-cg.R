# Clebsch-Gordan coefficients and Wigner symbols, after
# sympy/physics/quantum/cg.py, with the exact evaluation of
# sympy/physics/wigner.py (wigner_3j, clebsch_gordan, wigner_6j, racah,
# wigner_9j).

.wcomb <- function(n, k) if (k < 0 || k > n) gmp::as.bigz(0) else gmp::chooseZ(n, k)

.wdoubled <- function(value) {
  v <- as_sym(value)
  d <- .mul(list(S_Two, v))
  # as_int(value*2, strict=False) accepts integral floats
  if (d$type == "Float" && is.finite(as.numeric(d)) && as.numeric(d) %% 1 == 0) return(as.numeric(d))
  if (d$type != "Integer") .qstop("ValueError", "expecting integer or half-integer, got ", .str(v))
  as.numeric(d$value)
}

.wcheck_couple <- function(dj1, dj2, dj3) {
  if ((dj1 + dj2 + dj3) %% 2 != 0) return(FALSE)
  if (dj1 > dj2 + dj3 || dj2 > dj1 + dj3 || dj3 > dj1 + dj2) return(FALSE)
  TRUE
}

.wsqrt_q <- function(num, den) .sqrt(as_sym(gmp::as.bigq(num, den)))
.wint <- function(z) as_sym(gmp::as.bigz(z))

#' Wigner symbols and Clebsch-Gordan coefficients
#'
#' Exact values of the Wigner 3j, 6j and 9j symbols, the Clebsch-Gordan
#' coefficients and the Racah W coefficients for integer or half-integer
#' arguments, as SymPy's `sympy.physics.wigner`.
#'
#' @param j_1,j_2,j_3,j_4,j_5,j_6,j_7,j_8,j_9 Angular momenta.
#' @param m_1,m_2,m_3 Projections.
#' @param aa,bb,cc,dd,ee,ff The arguments of the Racah symbol.
#' @return An exact number.
#' @examples
#' wigner_3j(2, 6, 4, 0, 0, 0)
#' clebsch_gordan(S(1)/2, S(1)/2, 1, S(1)/2, S(1)/2, 1)
#' wigner_6j(3, 3, 3, 3, 3, 3)
#' @name wigner-symbols
NULL

#' @rdname wigner-symbols
#' @export
wigner_3j <- function(j_1, j_2, j_3, m_1, m_2, m_3) {
  d <- vapply(list(j_1, j_2, j_3, m_1, m_2, m_3), .wdoubled, 0)
  dj1 <- d[1]; dj2 <- d[2]; dj3 <- d[3]; dm1 <- d[4]; dm2 <- d[5]; dm3 <- d[6]
  if (dm1 + dm2 + dm3 != 0) return(S_Zero)
  if (!.wcheck_couple(dj1, dj2, dj3)) return(S_Zero)
  if (abs(dm1) > dj1 || abs(dm2) > dj2 || abs(dm3) > dj3) return(S_Zero)
  if (!((dj1 - dm1) %% 2 == 0 && (dj2 - dm2) %% 2 == 0 && (dj3 - dm3) %% 2 == 0)) return(S_Zero)
  sumj <- (dj1 + dj2 + dj3) %/% 2
  jm1 <- sumj - dj1; jm2 <- sumj - dj2; jm3 <- sumj - dj3
  j1mm1 <- (dj1 - dm1) %/% 2; j2mm2 <- (dj2 - dm2) %/% 2; j3mm3 <- (dj3 - dm3) %/% 2
  j1pm1 <- (dj1 + dm1) %/% 2
  imin <- max(0, j1pm1 - jm2, j2mm2 - jm1)
  imax <- min(jm3, j1pm1, j2mm2)
  sumres <- gmp::as.bigz(0)
  if (imin <= imax) for (ii in imin:imax)
    sumres <- .wcomb(jm3, ii) * .wcomb(jm2, j1pm1 - ii) * .wcomb(jm1, j2mm2 - ii) - sumres
  if (sumres == 0) return(S_Zero)
  res_num <- .wcomb(dj1, jm2) * .wcomb(dj2, jm1)
  res_den <- .wcomb(sumj, jm3) * .wcomb(dj1, j1mm1) * .wcomb(dj2, j2mm2) * .wcomb(dj3, j3mm3) * (sumj + 1)
  ressqrt <- .wsqrt_q(res_num, res_den)
  phase <- if ((dj1 + (dj3 + dm3) %/% 2 + imax) %% 2 == 0) S_One else S_NegativeOne
  .mul(list(phase, ressqrt, .wint(sumres)))
}

#' @rdname wigner-symbols
#' @export
clebsch_gordan <- function(j_1, j_2, j_3, m_1, m_2, m_3) {
  w <- wigner_3j(j_1, j_2, j_3, m_1, m_2, .neg(as_sym(m_3)))
  .mul(list(.pow(S_NegativeOne, .add(list(as_sym(j_1), .neg(as_sym(j_2)), as_sym(m_3)))),
            .sqrt(.add(list(.mul(list(S_Two, as_sym(j_3))), S_One))), w))
}

#' @rdname wigner-symbols
#' @export
wigner_6j <- function(j_1, j_2, j_3, j_4, j_5, j_6) {
  d <- vapply(list(j_1, j_2, j_3, j_4, j_5, j_6), .wdoubled, 0)
  dj1 <- d[1]; dj2 <- d[2]; dj3 <- d[3]; dj4 <- d[4]; dj5 <- d[5]; dj6 <- d[6]
  if (!.wcheck_couple(dj1, dj2, dj3) || !.wcheck_couple(dj1, dj5, dj6) ||
      !.wcheck_couple(dj4, dj2, dj6) || !.wcheck_couple(dj4, dj5, dj3)) return(S_Zero)
  j123 <- (dj1 + dj2 + dj3) %/% 2; j156 <- (dj1 + dj5 + dj6) %/% 2
  j426 <- (dj4 + dj2 + dj6) %/% 2; j453 <- (dj4 + dj5 + dj3) %/% 2
  jpm123 <- (dj1 + dj2 - dj3) %/% 2; jpm132 <- (dj1 + dj3 - dj2) %/% 2; jpm231 <- (dj2 + dj3 - dj1) %/% 2
  jpm156 <- (dj1 + dj5 - dj6) %/% 2; jpm426 <- (dj4 + dj2 - dj6) %/% 2; jpm453 <- (dj4 + dj5 - dj3) %/% 2
  sumres <- gmp::as.bigz(0)
  imin <- max(j123, j453, j426, j156)
  imax <- max(jpm123 + j453, jpm132 + j426, jpm231 + j156)
  if (imin <= imax) for (ii in imin:imax)
    sumres <- .wcomb(ii + 1, j123 + 1) * .wcomb(jpm123, ii - j453) * .wcomb(jpm132, ii - j426) * .wcomb(jpm231, ii - j156) - sumres
  if (sumres == 0) return(S_Zero)
  if (imax %% 2 == 1) sumres <- -sumres
  res_num <- .wcomb(j123 + 1, dj1 + 1) * .wcomb(dj1, jpm123)
  res_den <- .wcomb(j156 + 1, dj1 + 1) * .wcomb(dj1, jpm156) * .wcomb(j426 + 1, dj4 + 1) * .wcomb(dj4, jpm426) *
    .wcomb(j453 + 1, dj4 + 1) * .wcomb(dj4, jpm453)
  res_den <- res_den * (dj4 + 1)^2
  .mul(list(.wsqrt_q(res_num, res_den), .wint(sumres)))
}

#' @rdname wigner-symbols
#' @export
racah <- function(aa, bb, cc, dd, ee, ff) {
  res <- wigner_6j(aa, bb, ee, dd, cc, ff)
  if (.qeq(res, S_Zero)) return(S_Zero)
  s <- trunc(as.numeric(N(.add(lapply(list(aa, bb, cc, dd), as_sym)))))
  .mul(list(if (s %% 2 == 0) S_One else S_NegativeOne, res))
}

#' @rdname wigner-symbols
#' @export
wigner_9j <- function(j_1, j_2, j_3, j_4, j_5, j_6, j_7, j_8, j_9) {
  d <- vapply(list(j_1, j_2, j_3, j_4, j_5, j_6, j_7, j_8, j_9), .wdoubled, 0)
  dj1 <- d[1]; dj2 <- d[2]; dj3 <- d[3]; dj4 <- d[4]; dj5 <- d[5]; dj6 <- d[6]; dj7 <- d[7]; dj8 <- d[8]; dj9 <- d[9]
  if (!.wcheck_couple(dj1, dj2, dj3) || !.wcheck_couple(dj4, dj5, dj6) || !.wcheck_couple(dj7, dj8, dj9) ||
      !.wcheck_couple(dj1, dj4, dj7) || !.wcheck_couple(dj2, dj5, dj8) || !.wcheck_couple(dj3, dj6, dj9)) return(S_Zero)
  j123 <- (dj1 + dj2 + dj3) %/% 2; j456 <- (dj4 + dj5 + dj6) %/% 2; j789 <- (dj7 + dj8 + dj9) %/% 2
  j147 <- (dj1 + dj4 + dj7) %/% 2; j258 <- (dj2 + dj5 + dj8) %/% 2; j369 <- (dj3 + dj6 + dj9) %/% 2
  pm123 <- (dj1 + dj2 - dj3) %/% 2; pm132 <- (dj1 + dj3 - dj2) %/% 2; pm231 <- (dj2 + dj3 - dj1) %/% 2
  pm456 <- (dj4 + dj5 - dj6) %/% 2; pm465 <- (dj4 + dj6 - dj5) %/% 2; pm564 <- (dj5 + dj6 - dj4) %/% 2
  pm789 <- (dj7 + dj8 - dj9) %/% 2; pm798 <- (dj7 + dj9 - dj8) %/% 2; pm897 <- (dj8 + dj9 - dj7) %/% 2
  sumres <- gmp::as.bigz(0)
  dtmin <- max(abs(dj2 - dj6), abs(dj4 - dj8), abs(dj1 - dj9))
  dtmax <- min(dj2 + dj6, dj4 + dj8, dj1 + dj9)
  if (dtmin <= dtmax) for (dt in seq(dtmin, dtmax, by = 2)) {
    j19t <- (dj1 + dj9 + dt) %/% 2; j26t <- (dj2 + dj6 + dt) %/% 2; j48t <- (dj4 + dj8 + dt) %/% 2
    abc <- gmp::as.bigz(dt + 1)
    xmin <- max(j123, j369, j26t, j19t); xmax <- min(pm123 + j369, pm132 + j26t, pm231 + j19t)
    suma <- gmp::as.bigz(0)
    if (xmin <= xmax) for (x in xmin:xmax)
      suma <- .wcomb(x + 1, j19t + 1) * .wcomb((dj1 + dj9 - dt) %/% 2, x - j26t) * .wcomb((dj1 + dt - dj9) %/% 2, x - j369) *
        .wcomb((dt + dj9 - dj1) %/% 2, x - j123) - suma
    ymin <- max(j456, j258, j26t, j48t); ymax <- min(pm456 + j26t, pm465 + j258, pm564 + j48t)
    sumb <- gmp::as.bigz(0)
    if (ymin <= ymax) for (y in ymin:ymax)
      sumb <- .wcomb(y + 1, j26t + 1) * .wcomb((dj2 + dj6 - dt) %/% 2, y - j48t) * .wcomb((dt + dj6 - dj2) %/% 2, y - j258) *
        .wcomb((dj2 + dt - dj6) %/% 2, y - j456) - sumb
    zmin <- max(j789, j147, j48t, j19t); zmax <- min(pm789 + j19t, pm798 + j48t, pm897 + j147)
    sumc <- gmp::as.bigz(0)
    if (zmin <= zmax) for (z in zmin:zmax)
      sumc <- .wcomb(z + 1, j48t + 1) * .wcomb((dj4 + dj8 - dt) %/% 2, z - j19t) * .wcomb((dt + dj8 - dj4) %/% 2, z - j147) *
        .wcomb((dj4 + dt - dj8) %/% 2, z - j789) - sumc
    abc <- abc * suma * sumb * sumc
    if ((xmax + ymax + zmax) %% 2 == 1) abc <- -abc
    sumres <- sumres + abc
  }
  if (sumres == 0) return(S_Zero)
  if (dtmax %% 2 == 1) sumres <- -sumres
  res_den <- gmp::as.bigz((dj3 + 1) * (dj6 + 1) * (dj9 + 1) * (dj7 + 1) * (dj8 + 1) * (dj9 + 1))
  res_den <- res_den * .wcomb(j123 + 1, dj3 + 1) * .wcomb(dj3, pm231) * .wcomb(j456 + 1, dj6 + 1) * .wcomb(dj6, pm564) *
    .wcomb(j789 + 1, dj9 + 1) * .wcomb(dj9, pm897) * .wcomb(j147 + 1, dj7 + 1) * .wcomb(dj7, (dj4 + dj7 - dj1) %/% 2) *
    .wcomb(j258 + 1, dj8 + 1) * .wcomb(dj8, (dj5 + dj8 - dj2) %/% 2) * .wcomb(j369 + 1, dj9 + 1) * .wcomb(dj9, (dj6 + dj9 - dj3) %/% 2)
  .mul(list(.wint(sumres), .wsqrt_q(1, res_den)))
}

# ---------------------------------------------------------------------------
# CG, Wigner3j, Wigner6j, Wigner9j: commutative Expr subclasses printed as
# Basic (ClassName(args)) by str and srepr

.qbasic_printers <- list(
  `_sympystr` = function(self, p) paste0(.qpyname(self$type), "(", paste(vapply(self$args, .str, ""), collapse = ", "), ")"),
  `_sympyrepr` = function(self, p) paste0(.qpyname(self$type), "(", paste(vapply(self$args, srepr, ""), collapse = ", "), ")"))

.qexpr_new <- .cm(function(cls, args, ...) {
  .qm_state$on <- TRUE
  new_sym(cls, lapply(args, as_sym))
})

.qis_symbolic_args <- function(self) !all(vapply(self$args, function(a) .expr_is_number(a), NA))

.qwigner_attrs <- function(names) {
  out <- lapply(seq_along(names), function(i) .prop(function(self) self$args[[i]]))
  names(out) <- names
  out
}

do.call(.qclass, c(list("Wigner3j", "Expr", `__new__` = .qexpr_new, is_commutative = TRUE,
  is_symbolic = .prop(.qis_symbolic_args),
  `_latex` = function(self, p) {
    l <- vapply(self$args[c(1, 3, 5, 2, 4, 6)], .latex, "")
    sprintf("\\left(\\begin{array}{ccc} %s & %s & %s \\\\ %s & %s & %s \\end{array}\\right)", l[1], l[2], l[3], l[4], l[5], l[6])
  },
  doit = function(self, ...) {
    if (.qa(self, "is_symbolic")) .qstop("ValueError", "Coefficients must be numerical")
    a <- self$args
    wigner_3j(a[[1]], a[[3]], a[[5]], a[[2]], a[[4]], a[[6]])
  }), .qbasic_printers, .qwigner_attrs(c("j1", "m1", "j2", "m2", "j3", "m3"))))

.qclass("CG", "Wigner3j", precedence = 59,
  doit = function(self, ...) {
    if (.qa(self, "is_symbolic")) .qstop("ValueError", "Coefficients must be numerical")
    a <- self$args
    clebsch_gordan(a[[1]], a[[3]], a[[5]], a[[2]], a[[4]], a[[6]])
  },
  `_latex` = function(self, p) {
    l <- vapply(self$args[c(5, 6, 1, 2, 3, 4)], .latex, "")
    sprintf("C^{%s,%s}_{%s,%s,%s,%s}", l[1], l[2], l[3], l[4], l[5], l[6])
  })

do.call(.qclass, c(list("Wigner6j", "Expr", `__new__` = .qexpr_new, is_commutative = TRUE,
  is_symbolic = .prop(.qis_symbolic_args),
  `_latex` = function(self, p) {
    l <- vapply(self$args, .latex, "")
    sprintf("\\left\\{\\begin{array}{ccc} %s & %s & %s \\\\ %s & %s & %s \\end{array}\\right\\}", l[1], l[2], l[3], l[4], l[5], l[6])
  },
  doit = function(self, ...) {
    if (.qa(self, "is_symbolic")) .qstop("ValueError", "Coefficients must be numerical")
    do.call(wigner_6j, self$args)
  }), .qbasic_printers, .qwigner_attrs(c("j1", "j2", "j12", "j3", "j", "j23"))))

do.call(.qclass, c(list("Wigner9j", "Expr", `__new__` = .qexpr_new, is_commutative = TRUE,
  is_symbolic = .prop(.qis_symbolic_args),
  `_latex` = function(self, p) {
    l <- vapply(self$args, .latex, "")
    sprintf("\\left\\{\\begin{array}{ccc} %s & %s & %s \\\\ %s & %s & %s \\\\ %s & %s & %s \\end{array}\\right\\}",
            l[1], l[2], l[3], l[4], l[5], l[6], l[7], l[8], l[9])
  },
  doit = function(self, ...) {
    if (.qa(self, "is_symbolic")) .qstop("ValueError", "Coefficients must be numerical")
    do.call(wigner_9j, self$args)
  }), .qbasic_printers, .qwigner_attrs(c("j1", "j2", "j12", "j3", "j4", "j34", "j13", "j24", "j"))))

for (.t in c("Wigner3j", "CG", "Wigner6j", "Wigner9j"))
  assign(.t, function(e) c(commutative = TRUE), envir = .facts_hooks)
rm(.t)

#' Clebsch-Gordan coefficients and Wigner symbols as expressions
#'
#' `CG(j1, m1, j2, m2, j3, m3)` is the unevaluated Clebsch-Gordan
#' coefficient and `Wigner3j()`, `Wigner6j()` and `Wigner9j()` the Wigner
#' symbols of `sympy.physics.quantum.cg`; [doit()] evaluates them for
#' numeric arguments. `cg_simp()` simplifies sums and products of
#' Clebsch-Gordan coefficients with the Varshalovich identities.
#'
#' @param j1,m1,j2,m2,j3,m3,j12,j23,j,j4,j34,j13,j24 Angular momenta and projections.
#' @param e An expression.
#' @return An expression.
#' @examples
#' cg <- CG(1, 0, 1, 0, 2, 0)
#' doit(cg)
#' doit(Wigner3j(2, 0, 6, 0, 4, 0))
#' @name cg-coefficients
NULL

#' @rdname cg-coefficients
#' @export
Wigner3j <- function(j1, m1, j2, m2, j3, m3) .qfunc("Wigner3j", list(j1, m1, j2, m2, j3, m3))
#' @rdname cg-coefficients
#' @export
CG <- function(j1, m1, j2, m2, j3, m3) .qfunc("CG", list(j1, m1, j2, m2, j3, m3))
#' @rdname cg-coefficients
#' @export
Wigner6j <- function(j1, j2, j12, j3, j, j23) .qfunc("Wigner6j", list(j1, j2, j12, j3, j, j23))
#' @rdname cg-coefficients
#' @export
Wigner9j <- function(j1, j2, j12, j3, j4, j34, j13, j24, j) .qfunc("Wigner9j", list(j1, j2, j12, j3, j4, j34, j13, j24, j))

# ---------------------------------------------------------------------------
# cg_simp

#' @rdname cg-coefficients
#' @export
cg_simp <- function(e) {
  e <- as_sym(e)
  if (e$type == "Add") return(.cg_simp_add(e))
  if (e$type == "Sum") return(.cg_simp_sum(e))
  if (e$type == "Mul") return(.mul(lapply(.sympy_args(e), cg_simp)))
  if (e$type == "Pow") return(.pow(cg_simp(e$args[[1L]]), e$args[[2L]]))
  e
}

.has_cg <- function(e) .has_type(e, "CG")

.cg_simp_add <- function(e) {
  cg_part <- list(); other_part <- list()
  e <- expand(e)
  for (arg in .sympy_args(e)) {
    if (.has_cg(arg)) {
      if (arg$type == "Sum") {
        other_part[[length(other_part) + 1L]] <- .cg_simp_sum(arg)
      } else if (arg$type == "Mul") {
        terms <- S_One
        for (term in .sympy_args(arg)) terms <- .mul(list(terms, if (term$type == "Sum") .cg_simp_sum(term) else term))
        if (.has_cg(terms)) cg_part[[length(cg_part) + 1L]] <- terms else other_part[[length(other_part) + 1L]] <- terms
      } else {
        cg_part[[length(cg_part) + 1L]] <- arg
      }
    } else {
      other_part[[length(other_part) + 1L]] <- arg
    }
  }
  r <- .check_varsh_871_1(cg_part); cg_part <- r[[1L]]; other_part[[length(other_part) + 1L]] <- r[[2L]]
  r <- .check_varsh_871_2(cg_part); cg_part <- r[[1L]]; other_part[[length(other_part) + 1L]] <- r[[2L]]
  r <- .check_varsh_872_9(cg_part); cg_part <- r[[1L]]; other_part[[length(other_part) + 1L]] <- r[[2L]]
  .add(list(.add(cg_part), .add(other_part)))
}

.W <- function(n) Wild(n)

.check_varsh_871_1 <- function(term_list) {
  a <- .W("a"); alpha <- .W("alpha"); b <- .W("b"); lt <- .W("lt")
  expr <- .mul(list(lt, CG(a, alpha, b, 0, a, alpha)))
  simp <- .mul(list(.add(list(.mul(list(S_Two, a)), S_One)), KroneckerDelta(b, 0)))
  sign <- .mul(list(lt, .pow(abs(lt), S_NegativeOne)))
  .check_cg_simp(expr, simp, sign, lt, term_list, list(a, alpha, b, lt), list(a, b),
                 .add(list(.mul(list(S_Two, a)), S_One)), .add(list(a, alpha)))
}

.check_varsh_871_2 <- function(term_list) {
  a <- .W("a"); alpha <- .W("alpha"); c <- .W("c"); lt <- .W("lt")
  expr <- .mul(list(lt, CG(a, alpha, a, .neg(alpha), c, 0)))
  simp <- .mul(list(.sqrt(.add(list(.mul(list(S_Two, a)), S_One))), KroneckerDelta(c, 0)))
  sign <- .mul(list(.pow(S_NegativeOne, .sub(a, alpha)), lt, .pow(abs(lt), S_NegativeOne)))
  .check_cg_simp(expr, simp, sign, lt, term_list, list(a, alpha, c, lt), list(a, c),
                 .add(list(.mul(list(S_Two, a)), S_One)), .add(list(a, alpha)))
}

.check_varsh_872_9 <- function(term_list) {
  a <- .W("a"); alpha <- .W("alpha"); alphap <- .W("alphap"); b <- .W("b"); beta <- .W("beta")
  betap <- .W("betap"); c <- .W("c"); gamma <- .W("gamma"); lt <- .W("lt")
  expr <- .mul(list(lt, .pow(CG(a, alpha, b, beta, c, gamma), S_Two)))
  simp <- S_One
  sign <- .mul(list(lt, .pow(abs(lt), S_NegativeOne)))
  x <- abs(.sub(a, b)); y <- abs(.add(list(alpha, beta)))
  build_expr <- .sub(.add(list(a, b, S_One)), Piecewise(list(x, Gt(x, y)), list(S_Zero, Eq(x, y)), list(y, Gt(y, x))))
  index_expr <- .sub(.add(list(a, b)), c)
  r <- .check_cg_simp(expr, simp, sign, lt, term_list, list(a, alpha, b, beta, c, gamma, lt), list(a, alpha, b, beta),
                      build_expr, index_expr)
  term_list <- r[[1L]]; other1 <- r[[2L]]
  x <- abs(.sub(a, b)); y <- .add(list(a, b))
  build_expr <- .mul(list(.sub(.add(list(y, S_One)), x), .add(list(x, y, S_One))))
  index_expr <- .add(list(.mul(list(.sub(c, x), .add(list(x, c)))), c, gamma))
  r <- .check_cg_simp(expr, simp, sign, lt, term_list, list(a, alpha, b, beta, c, gamma, lt), list(a, alpha, b, beta),
                      build_expr, index_expr)
  term_list <- r[[1L]]; other2 <- r[[2L]]
  expr <- .mul(list(CG(a, alpha, b, beta, c, gamma), CG(a, alphap, b, betap, c, gamma)))
  simp <- .mul(list(KroneckerDelta(alpha, alphap), KroneckerDelta(beta, betap)))
  sign <- S_One
  x <- abs(.sub(a, b)); y <- abs(.add(list(alpha, beta)))
  build_expr <- .sub(.add(list(a, b, S_One)), Piecewise(list(x, Gt(x, y)), list(S_Zero, Eq(x, y)), list(y, Gt(y, x))))
  index_expr <- .sub(.add(list(a, b)), c)
  r <- .check_cg_simp(expr, simp, sign, S_One, term_list, list(a, alpha, alphap, b, beta, betap, c, gamma),
                      list(a, alpha, alphap, b, beta, betap), build_expr, index_expr)
  term_list <- r[[1L]]
  x <- abs(.sub(a, b)); y <- .add(list(a, b))
  build_expr <- .mul(list(.sub(.add(list(y, S_One)), x), .add(list(x, y, S_One))))
  index_expr <- .add(list(.mul(list(.sub(c, x), .add(list(x, c)))), c, gamma))
  r <- .check_cg_simp(expr, simp, sign, S_One, term_list, list(a, alpha, alphap, b, beta, betap, c, gamma),
                      list(a, alpha, alphap, b, beta, betap), build_expr, index_expr)
  term_list <- r[[1L]]; other4 <- r[[2L]]
  list(term_list, .add(list(other1, other2, other4)))
}

.qsubs_d <- function(e, d) {
  if (!length(d)) return(e)
  subs(e, lapply(names(d), function(n) list(Wild(n), d[[n]])), simultaneous = TRUE)
}

.check_cg_simp <- function(expr, simp, sign, lt, term_list, variables, dep_variables, build_index_expr, index_expr) {
  other_part <- S_Zero
  i <- 1L
  while (i <= length(term_list)) {
    sub_1 <- .check_cg(term_list[[i]], expr, length(variables))
    if (is.null(sub_1)) { i <- i + 1L; next }
    bi <- .qsubs_d(build_index_expr, sub_1)
    if (!.expr_is_number(bi)) { i <- i + 1L; next }
    sub_dep <- sub_1[vapply(dep_variables, function(x) x$name, "")]
    cg_index <- vector("list", .qint(bi))
    for (j in seq(i, length(term_list))) {
      sub_2 <- .check_cg(term_list[[j]], .qsubs_d(expr, sub_dep), length(variables) - length(dep_variables),
                         sign = list(.qsubs_d(sign, sub_1), .qsubs_d(sign, sub_dep)))
      if (is.null(sub_2)) next
      ie <- .qsubs_d(.qsubs_d(index_expr, sub_dep), sub_2)
      if (!.expr_is_number(ie)) next
      k <- .qint(ie) + 1
      if (k >= 1 && k <= length(cg_index))
        cg_index[[k]] <- list(j, .qsubs_d(.qsubs_d(subs(expr, lt, S_One), sub_dep), sub_2), .qsubs_d(lt, sub_2),
                              .qsubs_d(.qsubs_d(sign, sub_dep), sub_2))
    }
    if (!any(vapply(cg_index, is.null, NA))) {
      lts <- lapply(cg_index, function(t) abs(t[[3L]]))
      min_lt <- lts[[1L]]
      for (z in lts[-1L]) if (isTRUE(Lt(z, min_lt) == S_true)) min_lt <- z
      indices <- sort(vapply(cg_index, function(t) t[[1L]], 0), decreasing = TRUE)
      for (j in indices) term_list[[j]] <- NULL
      for (term in cg_index)
        if (isTRUE(Gt(abs(term[[3L]]), min_lt) == S_true))
          term_list[[length(term_list) + 1L]] <- .mul(list(.sub(term[[3L]], .mul(list(min_lt, term[[4L]]))), term[[2L]]))
      other_part <- .add(list(other_part, .mul(list(min_lt, .qsubs_d(.mul(list(sign, simp)), sub_1)))))
    } else {
      i <- i + 1L
    }
  }
  list(term_list, other_part)
}

.check_cg <- function(cg_term, expr, length, sign = NULL) {
  matches <- sym_match(cg_term, expr)
  if (is.null(matches)) return(NULL)
  if (!is.null(sign)) {
    if (!.qeq(sign[[1L]], .qsubs_d(sign[[2L]], matches))) return(NULL)
  }
  if (length(matches) == length) return(matches)
  NULL
}

.cg_simp_sum <- function(e) {
  e <- .check_varsh_sum_871_1(e)
  e <- .check_varsh_sum_871_2(e)
  e <- .check_varsh_sum_872_4(e)
  e
}

.check_varsh_sum_871_1 <- function(e) {
  a <- Wild("a"); alpha <- Symbol("alpha"); b <- Wild("b")
  m <- sym_match(e, Sum(CG(a, alpha, b, 0, a, alpha), list(alpha, .neg(a), a)))
  if (!is.null(m) && length(m) == 2L) return(.qsubs_d(.mul(list(.add(list(.mul(list(S_Two, a)), S_One)), KroneckerDelta(b, 0))), m))
  e
}

.check_varsh_sum_871_2 <- function(e) {
  a <- Wild("a"); alpha <- Symbol("alpha"); c <- Wild("c")
  m <- sym_match(e, Sum(.mul(list(.pow(S_NegativeOne, .sub(a, alpha)), CG(a, alpha, a, .neg(alpha), c, 0))), list(alpha, .neg(a), a)))
  if (!is.null(m) && length(m) == 2L) return(.qsubs_d(.mul(list(.sqrt(.add(list(.mul(list(S_Two, a)), S_One))), KroneckerDelta(c, 0))), m))
  e
}

.check_varsh_sum_872_4 <- function(e) {
  alpha <- Symbol("alpha"); beta <- Symbol("beta")
  a <- Wild("a"); b <- Wild("b"); c <- Wild("c"); cp <- Wild("cp"); gamma <- Wild("gamma"); gammap <- Wild("gammap")
  cg1 <- CG(a, alpha, b, beta, c, gamma)
  cg2 <- CG(a, alpha, b, beta, cp, gammap)
  m1 <- sym_match(e, Sum(.mul(list(cg1, cg2)), list(alpha, .neg(a), a), list(beta, .neg(b), b)))
  if (!is.null(m1) && length(m1) == 6L) return(.qsubs_d(.mul(list(KroneckerDelta(c, cp), KroneckerDelta(gamma, gammap))), m1))
  m2 <- sym_match(e, Sum(.pow(cg1, S_Two), list(alpha, .neg(a), a), list(beta, .neg(b), b)))
  if (!is.null(m2) && length(m2) == 4L) return(S_One)
  e
}
