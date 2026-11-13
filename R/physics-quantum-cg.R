# Clebsch-Gordan coefficients and Wigner symbols, after
# sympy/physics/quantum/cg.py, with the exact evaluation of
# sympy/physics/wigner.py (wigner_3j, clebsch_gordan, wigner_6j, racah,
# wigner_9j).


# wigner_3j, clebsch_gordan, wigner_6j, racah and wigner_9j live in
# physics-misc-wigner.R (sympy/physics/wigner.py).

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
