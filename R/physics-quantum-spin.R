# Angular momentum: spin operators, rotations, Wigner D functions and the
# coupled and uncoupled spin states, after sympy/physics/quantum/spin.py.

#' Spin quantum numbers
#'
#' `m_values(j)` gives the dimension `2*j + 1` and the values `j, j - 1,
#' ..., -j` of the projection quantum number, as SymPy's
#' `spin.m_values`.
#'
#' @param j An integer or half-integer.
#' @return A list of the size and the list of m values.
#' @examples
#' m_values(S(3) / 2)
#' @export
m_values <- function(j) {
  j <- as_sym(j)
  size <- .add(list(.mul(list(S_Two, j)), S_One))
  if (!(size$type == "Integer" && num_sign(size) > 0))
    .qstop("ValueError", "Only integer or half-integer values allowed for j, got: : ", .qrepr(j))
  n <- .qint(size)
  list(size, lapply(seq_len(n) - 1L, function(i) .sub(j, as_sym(i))))
}

.qkd <- function(a, b) KroneckerDelta(a, b)

# ---------------------------------------------------------------------------
# Spin operators

.qclass("SpinOpBase", character(), node = FALSE,
  `_eval_hilbert_space` = .cm(function(cls, label) ComplexSpace(S_Infinity)),
  name = .prop(function(self) self$args[[1L]]),
  `_print_contents` = function(self, p) paste0(.str(.qa(self, "name")), .qa(self, "_coord")),
  `_print_contents_latex` = function(self, p) paste0(.str(.qa(self, "name")), "_", .qa(self, "_coord")),
  `_represent_base` = function(self, basis, j = S_Half, ...) {
    mv <- m_values(j)
    size <- .qint(mv[[1L]]); mvals <- mv[[2L]]
    j <- as_sym(j)
    vals <- vector("list", size * size)
    for (p in seq_len(size)) for (q in seq_len(size))
      vals[[(p - 1L) * size + q]] <- .qmeth(self, "matrix_element", j, mvals[[p]], j, mvals[[q]])
    .new_mat(vals, size, size)
  },
  `_apply_op` = function(self, ket, orig_basis, ...) {
    state <- .qrewrite(ket, .qa(self, "basis"))
    if (.qisa(state, "State")) {
      ret <- .mul(list(hbar, .qa(state, "m"), state))
    } else if (state$type == "Sum") {
      ret <- .qmeth(self, "_apply_operator_Sum", state, ...)
    } else {
      ret <- qapply(.mul(list(self, state)))
    }
    if (.qeq(ret, .mul(list(self, state)))) .qni("")
    .qrewrite(ret, orig_basis)
  },
  `_apply_operator_JxKet` = function(self, ket, ...) .qmeth(self, "_apply_op", ket, "Jx", ...),
  `_apply_operator_JxKetCoupled` = function(self, ket, ...) .qmeth(self, "_apply_op", ket, "Jx", ...),
  `_apply_operator_JyKet` = function(self, ket, ...) .qmeth(self, "_apply_op", ket, "Jy", ...),
  `_apply_operator_JyKetCoupled` = function(self, ket, ...) .qmeth(self, "_apply_op", ket, "Jy", ...),
  `_apply_operator_JzKet` = function(self, ket, ...) .qmeth(self, "_apply_op", ket, "Jz", ...),
  `_apply_operator_JzKetCoupled` = function(self, ket, ...) .qmeth(self, "_apply_op", ket, "Jz", ...),
  `_apply_operator_TensorProduct` = function(self, tp, ...) {
    if (!.qisa(self, c("JxOp", "JyOp", "JzOp"))) .qni("")
    result <- list()
    for (n in seq_along(tp$args)) {
      arg <- tp$args
      arg[[n]] <- .qmeth(self, "_apply_operator", tp$args[[n]])
      result[[n]] <- .qtp(arg)
    }
    expand(.add(result))
  },
  `_apply_operator_Sum` = function(self, s, ...) {
    new_func <- qapply(.mul(list(self, s$args[[1L]])))
    if (.qeq(new_func, .mul(list(self, s$args[[1L]])))) .qni("")
    do.call(Sum, c(list(new_func), lapply(s$args[-1L], function(l) l$args)))
  },
  `_eval_trace` = function(self, ...) {
    m <- .qmeth(self, "_represent_default_basis")
    .add(lapply(seq_len(m$nrow), function(i) .mget(m, i, i)))
  })

.qjpm_apply <- function(self, ket, sign, coupled) {
  j <- .qa(ket, "j"); m <- .qa(ket, "m")
  if (is_number_type(m) && is_number_type(j)) {
    if (sign > 0 && num_cmp(m, j) >= 0) return(S_Zero)
    if (sign < 0 && num_cmp(m, .neg(j)) <= 0) return(S_Zero)
  }
  s <- as_sym(sign)
  coeff <- .mul(list(hbar, .sqrt(.sub(.mul(list(j, .add(list(j, S_One)))), .mul(list(m, .add(list(m, s))))))))
  newm <- .add(list(m, s))
  st <- if (coupled) JzKetCoupled(j, newm, .qa(ket, "jn"), .qa(ket, "coupling")) else JzKet(j, newm)
  .mul(list(coeff, st))
}

.qclass("JplusOp", c("SpinOpBase", "Operator"),
  `_coord` = "+", basis = "Jz",
  `_eval_commutator_JminusOp` = function(self, other, ...) .mul(list(S_Two, hbar, JzOp(.qa(self, "name")))),
  `_apply_operator_JzKet` = function(self, ket, ...) .qjpm_apply(self, ket, 1, FALSE),
  `_apply_operator_JzKetCoupled` = function(self, ket, ...) .qjpm_apply(self, ket, 1, TRUE),
  matrix_element = function(self, j, m, jp, mp) {
    result <- .mul(list(hbar, .sqrt(.sub(.mul(list(j, .add(list(j, S_One)))), .mul(list(mp, .add(list(mp, S_One))))))))
    result <- .mul(list(result, .qkd(m, .add(list(mp, S_One)))))
    .mul(list(result, .qkd(j, jp)))
  },
  `_represent_default_basis` = function(self, ...) .qmeth(self, "_represent_JzOp", NULL, ...),
  `_represent_JzOp` = function(self, basis, ...) .qmeth(self, "_represent_base", basis, ...),
  `_eval_rewrite_as_xyz` = function(self, args, ...) .add(list(JxOp(args[[1L]]), .mul(list(sym_I, JyOp(args[[1L]]))))))

.qclass("JminusOp", c("SpinOpBase", "Operator"),
  `_coord` = "-", basis = "Jz",
  `_apply_operator_JzKet` = function(self, ket, ...) .qjpm_apply(self, ket, -1, FALSE),
  `_apply_operator_JzKetCoupled` = function(self, ket, ...) .qjpm_apply(self, ket, -1, TRUE),
  matrix_element = function(self, j, m, jp, mp) {
    result <- .mul(list(hbar, .sqrt(.sub(.mul(list(j, .add(list(j, S_One)))), .mul(list(mp, .sub(mp, S_One)))))))
    result <- .mul(list(result, .qkd(m, .sub(mp, S_One))))
    .mul(list(result, .qkd(j, jp)))
  },
  `_represent_default_basis` = function(self, ...) .qmeth(self, "_represent_JzOp", NULL, ...),
  `_represent_JzOp` = function(self, basis, ...) .qmeth(self, "_represent_base", basis, ...),
  `_eval_rewrite_as_xyz` = function(self, args, ...) .sub(JxOp(args[[1L]]), .mul(list(sym_I, JyOp(args[[1L]])))))

.qjxy_apply <- function(self, ket, nm, y, ...) {
  jp <- .qmeth(JplusOp(.qa(self, "name")), nm, ket, ...)
  jm <- .qmeth(JminusOp(.qa(self, "name")), nm, ket, ...)
  if (y) .mul(list(.sub(jp, jm), .pow(.mul(list(S_Two, sym_I)), S_NegativeOne)))
  else .mul(list(.add(list(jp, jm)), .pow(S_Two, S_NegativeOne)))
}
.qjxy_rep <- function(self, basis, y, ...) {
  jp <- .qmeth(JplusOp(.qa(self, "name")), "_represent_JzOp", basis, ...)
  jm <- .qmeth(JminusOp(.qa(self, "name")), "_represent_JzOp", basis, ...)
  if (y) .mat_map(jp - jm, function(z) .mul(list(z, .pow(.mul(list(S_Two, sym_I)), S_NegativeOne))))
  else .mat_map(jp + jm, function(z) .mul(list(z, .pow(S_Two, S_NegativeOne))))
}

.qclass("JxOp", c("SpinOpBase", "HermitianOperator"),
  `_coord` = "x", basis = "Jx",
  `_eval_commutator_JyOp` = function(self, other, ...) .mul(list(sym_I, hbar, JzOp(.qa(self, "name")))),
  `_eval_commutator_JzOp` = function(self, other, ...) .mul(list(S_NegativeOne, sym_I, hbar, JyOp(.qa(self, "name")))),
  `_apply_operator_JzKet` = function(self, ket, ...) .qjxy_apply(self, ket, "_apply_operator_JzKet", FALSE, ...),
  `_apply_operator_JzKetCoupled` = function(self, ket, ...) .qjxy_apply(self, ket, "_apply_operator_JzKetCoupled", FALSE, ...),
  `_represent_default_basis` = function(self, ...) .qmeth(self, "_represent_JzOp", NULL, ...),
  `_represent_JzOp` = function(self, basis, ...) .qjxy_rep(self, basis, FALSE, ...),
  `_eval_rewrite_as_plusminus` = function(self, args, ...)
    .mul(list(.add(list(JplusOp(args[[1L]]), JminusOp(args[[1L]]))), S_Half)))

.qclass("JyOp", c("SpinOpBase", "HermitianOperator"),
  `_coord` = "y", basis = "Jy",
  `_eval_commutator_JzOp` = function(self, other, ...) .mul(list(sym_I, hbar, JxOp(.qa(self, "name")))),
  `_eval_commutator_JxOp` = function(self, other, ...) .mul(list(S_NegativeOne, sym_I, hbar, J2Op(.qa(self, "name")))),
  `_apply_operator_JzKet` = function(self, ket, ...) .qjxy_apply(self, ket, "_apply_operator_JzKet", TRUE, ...),
  `_apply_operator_JzKetCoupled` = function(self, ket, ...) .qjxy_apply(self, ket, "_apply_operator_JzKetCoupled", TRUE, ...),
  `_represent_default_basis` = function(self, ...) .qmeth(self, "_represent_JzOp", NULL, ...),
  `_represent_JzOp` = function(self, basis, ...) .qjxy_rep(self, basis, TRUE, ...),
  `_eval_rewrite_as_plusminus` = function(self, args, ...)
    .mul(list(.sub(JplusOp(args[[1L]]), JminusOp(args[[1L]])), .pow(.mul(list(S_Two, sym_I)), S_NegativeOne))))

.qclass("JzOp", c("SpinOpBase", "HermitianOperator"),
  `_coord` = "z", basis = "Jz",
  `_eval_commutator_JxOp` = function(self, other, ...) .mul(list(sym_I, hbar, JyOp(.qa(self, "name")))),
  `_eval_commutator_JyOp` = function(self, other, ...) .mul(list(S_NegativeOne, sym_I, hbar, JxOp(.qa(self, "name")))),
  `_eval_commutator_JplusOp` = function(self, other, ...) .mul(list(hbar, JplusOp(.qa(self, "name")))),
  `_eval_commutator_JminusOp` = function(self, other, ...) .mul(list(S_NegativeOne, hbar, JminusOp(.qa(self, "name")))),
  matrix_element = function(self, j, m, jp, mp) .mul(list(hbar, mp, .qkd(m, mp), .qkd(j, jp))),
  `_represent_default_basis` = function(self, ...) .qmeth(self, "_represent_JzOp", NULL, ...),
  `_represent_JzOp` = function(self, basis, ...) .qmeth(self, "_represent_base", basis, ...))

.qj2_apply <- function(self, ket, ...) {
  j <- .qa(ket, "j")
  .mul(list(.pow(hbar, S_Two), j, .add(list(j, S_One)), ket))
}

.qclass("J2Op", c("SpinOpBase", "HermitianOperator"),
  `_coord` = "2",
  `_eval_commutator_JxOp` = function(self, other, ...) S_Zero,
  `_eval_commutator_JyOp` = function(self, other, ...) S_Zero,
  `_eval_commutator_JzOp` = function(self, other, ...) S_Zero,
  `_eval_commutator_JplusOp` = function(self, other, ...) S_Zero,
  `_eval_commutator_JminusOp` = function(self, other, ...) S_Zero,
  `_apply_operator_JxKet` = .qj2_apply, `_apply_operator_JxKetCoupled` = .qj2_apply,
  `_apply_operator_JyKet` = .qj2_apply, `_apply_operator_JyKetCoupled` = .qj2_apply,
  `_apply_operator_JzKet` = .qj2_apply, `_apply_operator_JzKetCoupled` = .qj2_apply,
  matrix_element = function(self, j, m, jp, mp)
    .mul(list(.pow(hbar, S_Two), j, .add(list(j, S_One)), .qkd(m, mp), .qkd(j, jp))),
  `_represent_default_basis` = function(self, ...) .qmeth(self, "_represent_JzOp", NULL, ...),
  `_represent_JzOp` = function(self, basis, ...) .qmeth(self, "_represent_base", basis, ...),
  `_print_contents_latex` = function(self, p) paste0(.str(.qa(self, "name")), "^2"),
  `_eval_rewrite_as_xyz` = function(self, args, ...)
    .add(list(.pow(JxOp(args[[1L]]), S_Two), .pow(JyOp(args[[1L]]), S_Two), .pow(JzOp(args[[1L]]), S_Two))),
  `_eval_rewrite_as_plusminus` = function(self, args, ...) {
    a <- args[[1L]]
    .add(list(.pow(JzOp(a), S_Two), .mul(list(S_Half, .add(list(.mul(list(JplusOp(a), JminusOp(a))), .mul(list(JminusOp(a), JplusOp(a)))))))))
  })

# ---------------------------------------------------------------------------
# Rotation and WignerD

.qclass("Rotation", "UnitaryOperator",
  `_eval_args` = .cm(function(cls, args, ...) {
    args <- .qsuper("QExpr", "_eval_args", cls, args)
    if (length(args) != 3L) .qstop("ValueError", "3 Euler angles required, got: ", .qrepr(args))
    args
  }),
  `_eval_hilbert_space` = .cm(function(cls, label) ComplexSpace(S_Infinity)),
  alpha = .prop(function(self) .qa(self, "label")[[1L]]),
  beta = .prop(function(self) .qa(self, "label")[[2L]]),
  gamma = .prop(function(self) .qa(self, "label")[[3L]]),
  `_print_operator_name` = function(self, p) "R",
  `_print_operator_name_latex` = function(self, p) "\\mathcal{R}",
  `_eval_inverse` = function(self) Rotation(.neg(.qa(self, "gamma")), .neg(.qa(self, "beta")), .neg(.qa(self, "alpha"))),
  D = .cm(function(cls, j, m, mp, alpha, beta, gamma) WignerD(j, m, mp, alpha, beta, gamma)),
  d = .cm(function(cls, j, m, mp, beta) WignerD(j, m, mp, 0, beta, 0)),
  matrix_element = function(self, j, m, jp, mp) {
    result <- WignerD(jp, m, mp, .qa(self, "alpha"), .qa(self, "beta"), .qa(self, "gamma"))
    .mul(list(result, .qkd(j, jp)))
  },
  `_represent_base` = function(self, basis, j = S_Half, doit = NULL, ...) {
    j <- as_sym(j)
    evaluate <- isTRUE(doit) || (is_sym(doit) && .qtruthy(doit))
    mv <- m_values(j)
    size <- .qint(mv[[1L]]); mvals <- mv[[2L]]
    vals <- vector("list", size * size)
    for (p in seq_len(size)) for (q in seq_len(size)) {
      me <- .qmeth(self, "matrix_element", j, mvals[[p]], j, mvals[[q]])
      vals[[(p - 1L) * size + q]] <- if (evaluate) doit(me) else me
    }
    .new_mat(vals, size, size)
  },
  `_represent_default_basis` = function(self, ...) .qmeth(self, "_represent_JzOp", NULL, ...),
  `_represent_JzOp` = function(self, basis, ...) .qmeth(self, "_represent_base", basis, ...),
  `_apply_operator_uncoupled` = function(self, state, ket, dummy = TRUE, ...) {
    a <- .qa(self, "alpha"); b <- .qa(self, "beta"); g <- .qa(self, "gamma")
    j <- .qa(ket, "j"); m <- .qa(ket, "m")
    if (.expr_is_number(j)) {
      s <- lapply(m_values(j)[[2L]], function(mp) .mul(list(doit(WignerD(j, m, mp, a, b, g)), .qfunc(state, list(j, mp)))))
      return(.add(s))
    }
    mp <- if (dummy) Dummy("mp") else Symbol("mp")
    Sum(.mul(list(WignerD(j, m, mp, a, b, g), .qfunc(state, list(j, mp)))), list(mp, .neg(j), j))
  },
  `_apply_operator_JxKet` = function(self, ket, ...) .qmeth(self, "_apply_operator_uncoupled", "JxKet", ket, ...),
  `_apply_operator_JyKet` = function(self, ket, ...) .qmeth(self, "_apply_operator_uncoupled", "JyKet", ket, ...),
  `_apply_operator_JzKet` = function(self, ket, ...) .qmeth(self, "_apply_operator_uncoupled", "JzKet", ket, ...),
  `_apply_operator_coupled` = function(self, state, ket, dummy = TRUE, ...) {
    a <- .qa(self, "alpha"); b <- .qa(self, "beta"); g <- .qa(self, "gamma")
    j <- .qa(ket, "j"); m <- .qa(ket, "m"); jn <- .qa(ket, "jn"); coupling <- .qa(ket, "coupling")
    if (.expr_is_number(j)) {
      s <- lapply(m_values(j)[[2L]], function(mp)
        .mul(list(doit(WignerD(j, m, mp, a, b, g)), .qfunc(state, list(j, mp, jn, coupling)))))
      return(.add(s))
    }
    mp <- if (dummy) Dummy("mp") else Symbol("mp")
    Sum(.mul(list(WignerD(j, m, mp, a, b, g), .qfunc(state, list(j, mp, jn, coupling)))), list(mp, .neg(j), j))
  },
  `_apply_operator_JxKetCoupled` = function(self, ket, ...) .qmeth(self, "_apply_operator_coupled", "JxKetCoupled", ket, ...),
  `_apply_operator_JyKetCoupled` = function(self, ket, ...) .qmeth(self, "_apply_operator_coupled", "JyKetCoupled", ket, ...),
  `_apply_operator_JzKetCoupled` = function(self, ket, ...) .qmeth(self, "_apply_operator_coupled", "JzKetCoupled", ket, ...))

do.call(.qclass, c(list("WignerD", "Expr",
  `__new__` = .cm(function(cls, args, evaluate = FALSE, ...) {
    if (length(args) != 6L) .qstop("ValueError", "6 parameters expected, got ", .qrepr(args))
    .qm_state$on <- TRUE
    obj <- new_sym("WignerD", lapply(args, as_sym))
    if (isTRUE(evaluate)) return(.qmeth(obj, "_eval_wignerd"))
    obj
  }),
  is_commutative = TRUE,
  `_latex` = function(self, p) {
    a <- self$args
    l <- vapply(a, .latex, "")
    if (.qeq(a[[4L]], S_Zero) && .qeq(a[[6L]], S_Zero))
      return(sprintf("d^{%s}_{%s,%s}\\left(%s\\right)", l[1], l[2], l[3], l[5]))
    sprintf("D^{%s}_{%s,%s}\\left(%s,%s,%s\\right)", l[1], l[2], l[3], l[4], l[5], l[6])
  },
  doit = function(self, ...) .qfunc("WignerD", self$args, evaluate = TRUE),
  `_eval_wignerd` = function(self) {
    a <- self$args
    j <- a[[1L]]; m <- a[[2L]]; mp <- a[[3L]]; alpha <- a[[4L]]; beta <- a[[5L]]; gamma <- a[[6L]]
    if (.qeq(alpha, S_Zero) && .qeq(beta, S_Zero) && .qeq(gamma, S_Zero)) return(KroneckerDelta(m, mp))
    if (!.expr_is_number(j)) .qstop("ValueError", "j parameter must be numerical to evaluate, got ", .str(j))
    r <- S_Zero
    if (.qeq(beta, .mul(list(sym_pi, S_Half)))) {
      jpm <- .add(list(j, mp)); jmm <- .sub(j, mp)
      for (k in seq_len(.qint(.add(list(.mul(list(S_Two, j)), S_One)))) - 1L) {
        ks <- as_sym(k)
        if (isTRUE(.qgt(ks, .add(list(j, mp)))) || isTRUE(.qgt(ks, .sub(j, m))) || isTRUE(.qgt(.sub(mp, m), ks))) next
        r <- .add(list(r, .mul(list(.pow(S_NegativeOne, ks), binomial(jpm, ks), binomial(jmm, .sub(.add(list(ks, m)), mp))))))
      }
      r <- .mul(list(r, .pow(S_NegativeOne, .sub(m, mp)), .pow(.pow(S_Two, j), S_NegativeOne),
                     .sqrt(.mul(list(factorial(.add(list(j, m))), factorial(.sub(j, m)),
                                     .pow(.mul(list(factorial(.add(list(j, mp))), factorial(.sub(j, mp)))), S_NegativeOne))))))
    } else {
      for (mpp in m_values(j)[[2L]]) {
        t1 <- doit(WignerD(j, m, mpp, 0, .mul(list(sym_pi, S_Half)), 0))
        t2 <- .add(list(cos(.mul(list(S_NegativeOne, mpp, beta))), .mul(list(sym_I, sin(.mul(list(S_NegativeOne, mpp, beta)))))))
        t3 <- doit(WignerD(j, mpp, .neg(mp), 0, .mul(list(sym_pi, S_Half)), 0))
        r <- .add(list(r, .mul(list(t1, t2, t3))))
      }
      r <- .mul(list(r, .pow(sym_I, .sub(.sub(.mul(list(S_Two, j)), m), mp)), .pow(S_NegativeOne, .mul(list(S_Two, m)))))
      r <- simplify(r)
    }
    .mul(list(r, exp(.mul(list(S_NegativeOne, sym_I, m, alpha))), exp(.mul(list(S_NegativeOne, sym_I, mp, gamma)))))
  }), .qbasic_printers, .qwigner_attrs(c("j", "m", "mp", "alpha", "beta", "gamma"))))
assign("WignerD", function(e) c(commutative = TRUE), envir = .facts_hooks)

.qgt <- function(a, b) {
  d <- .sub(a, b)
  if (is_number_type(d)) return(num_sign(d) > 0)
  NA
}

#' Spin operators, rotations and Wigner D functions
#'
#' The angular momentum operators of `sympy.physics.quantum.spin`:
#' `JxOp()`, `JyOp()`, `JzOp()`, `J2Op()`, the ladder operators `JplusOp()`
#' and `JminusOp()` (each takes a name; the objects `Jx`, `Jy`, `Jz`, `J2`,
#' `Jplus` and `Jminus` are these operators named `J`, as in SymPy), the
#' rotation `Rotation(alpha, beta, gamma)` and the Wigner D function
#' `WignerD(j, m, mp, alpha, beta, gamma)` (evaluated by [doit()]).
#' `Rotation_D()` and `Rotation_d()` are SymPy's `Rotation.D` and
#' `Rotation.d`.
#'
#' @param ... The name of the operator, or the Euler angles.
#' @param j,m,mp Angular momentum quantum numbers.
#' @param alpha,beta,gamma Euler angles.
#' @param evaluate Evaluate the Wigner D function.
#' @return A quantum operator or an expression.
#' @examples
#' qapply(Jz * JzKet(1, 1))
#' represent(JzOp("J"), j = 1)
#' doit(WignerD(1, 1, 0, 0, sym_pi / 2, 0))
#' @name spin-operators
NULL

#' @rdname spin-operators
#' @export
JxOp <- function(...) .qfunc("JxOp", list(...))
#' @rdname spin-operators
#' @export
JyOp <- function(...) .qfunc("JyOp", list(...))
#' @rdname spin-operators
#' @export
JzOp <- function(...) .qfunc("JzOp", list(...))
#' @rdname spin-operators
#' @export
J2Op <- function(...) .qfunc("J2Op", list(...))
#' @rdname spin-operators
#' @export
JplusOp <- function(...) .qfunc("JplusOp", list(...))
#' @rdname spin-operators
#' @export
JminusOp <- function(...) .qfunc("JminusOp", list(...))
#' @rdname spin-operators
#' @export
Jx <- NULL
delayedAssign("Jx", JxOp("J"))
#' @rdname spin-operators
#' @export
Jy <- NULL
delayedAssign("Jy", JyOp("J"))
#' @rdname spin-operators
#' @export
Jz <- NULL
delayedAssign("Jz", JzOp("J"))
#' @rdname spin-operators
#' @export
J2 <- NULL
delayedAssign("J2", J2Op("J"))
#' @rdname spin-operators
#' @export
Jplus <- NULL
delayedAssign("Jplus", JplusOp("J"))
#' @rdname spin-operators
#' @export
Jminus <- NULL
delayedAssign("Jminus", JminusOp("J"))
#' @rdname spin-operators
#' @export
Rotation <- function(...) .qfunc("Rotation", list(...))
#' @rdname spin-operators
#' @export
Rotation_D <- function(j, m, mp, alpha, beta, gamma) WignerD(j, m, mp, alpha, beta, gamma)
#' @rdname spin-operators
#' @export
Rotation_d <- function(j, m, mp, beta) WignerD(j, m, mp, 0, beta, 0)
#' @rdname spin-operators
#' @export
WignerD <- function(j, m, mp, alpha, beta, gamma, evaluate = FALSE)
  .qfunc("WignerD", list(j, m, mp, alpha, beta, gamma), evaluate = evaluate)

# ---------------------------------------------------------------------------
# Spin states

.qhalfint_bad <- function(x) { d <- .mul(list(S_Two, x)); !(d$type == "Integer") }

.qspin_validate <- function(j, m) {
  if (.expr_is_number(j)) {
    if (.qhalfint_bad(j)) .qstop("ValueError", "j must be integer or half-integer, got: ", .str(j))
    if (isTRUE(num_sign(N(j)) < 0)) .qstop("ValueError", "j must be >= 0, got: ", .str(j))
  }
  if (.expr_is_number(m)) {
    if (.qhalfint_bad(m)) .qstop("ValueError", "m must be integer or half-integer, got: ", .str(m))
  }
  if (.expr_is_number(j) && .expr_is_number(m)) {
    if (isTRUE(.qgt(abs(m), j))) .qstop("ValueError", "Allowed values for m are -j <= m <= j, got j, m: ", .str(j), ", ", .str(m))
    if (.sub(j, m)$type != "Integer")
      .qstop("ValueError", "Both j and m must be integer or half-integer, got j, m: ", .str(j), ", ", .str(m))
  }
}

.qrep_angles <- function(self, alpha = S_Zero, beta = S_Zero, gamma = S_Zero, ...) {
  .qmeth(self, "_represent_base", alpha = alpha, beta = beta, gamma = gamma, ...)
}

.qclass("SpinState", "State",
  `_label_separator` = ",",
  `__new__` = .cm(function(cls, args, ...) {
    if (length(args) != 2L) .qstop("TypeError", "SpinState.__new__() takes 3 positional arguments but ", length(args) + 1L, " were given")
    j <- .qS(args[[1L]]); m <- .qS(args[[2L]])
    .qspin_validate(j, m)
    .qsuper("State", "__new__", cls, list(j, m))
  }),
  j = .prop(function(self) .qa(self, "label")[[1L]]),
  m = .prop(function(self) .qa(self, "label")[[2L]]),
  `_eval_hilbert_space` = .cm(function(cls, label) ComplexSpace(.add(list(.mul(list(S_Two, label[[1L]])), S_One)))),
  `_represent_base` = function(self, alpha = S_Zero, beta = S_Zero, gamma = S_Zero, ...) {
    j <- .qa(self, "j"); m <- .qa(self, "m")
    alpha <- as_sym(alpha); beta <- as_sym(beta); gamma <- as_sym(gamma)
    mv <- m_values(j)
    size <- .qint(mv[[1L]])
    vals <- lapply(mv[[2L]], function(mval) {
      w <- WignerD(j, mval, m, alpha, beta, gamma)
      if (.expr_is_number(m)) doit(w) else w
    })
    .new_mat(vals, size, 1L)
  },
  `_eval_rewrite_as_Jx` = function(self, args, ...) {
    if (.qisa(self, "Bra")) return(.qmeth(self, "_rewrite_basis", JxOp("J"), "JxBra", ...))
    .qmeth(self, "_rewrite_basis", JxOp("J"), "JxKet", ...)
  },
  `_eval_rewrite_as_Jy` = function(self, args, ...) {
    if (.qisa(self, "Bra")) return(.qmeth(self, "_rewrite_basis", JyOp("J"), "JyBra", ...))
    .qmeth(self, "_rewrite_basis", JyOp("J"), "JyKet", ...)
  },
  `_eval_rewrite_as_Jz` = function(self, args, ...) {
    if (.qisa(self, "Bra")) return(.qmeth(self, "_rewrite_basis", JzOp("J"), "JzBra", ...))
    .qmeth(self, "_rewrite_basis", JzOp("J"), "JzKet", ...)
  },
  `_rewrite_basis` = function(self, basis, evect, coupled = NULL, ...) {
    j <- .qa(self, "j")
    args <- self$args[-(1:2)]
    if (.expr_is_number(j)) {
      start <- 0
      if (.qisa(self, "CoupledSpinState")) {
        jv <- as.numeric(N(j))
        start <- if (jv == round(jv)) jv^2 else (2 * jv - 1) * (2 * jv + 1) / 4
      }
      vect <- represent(self, basis = basis, ...)
      n <- .qint(.add(list(.mul(list(S_Two, j)), S_One)))
      result <- .add(lapply(seq_len(n) - 1L, function(i)
        .mul(list(vect$data[[start + i + 1]], .qfunc(evect, c(list(j, .sub(j, as_sym(i))), args))))))
      if (.qisa(self, "CoupledSpinState") && identical(coupled, FALSE)) return(uncouple(result))
      return(result)
    }
    mi <- Symbol("mi")
    # make sure not to introduce a symbol already in the state (the loop of
    # SymPy runs once at most)
    if (!.qeq(subs(self, mi, S_Zero), self)) mi <- Symbol("mi1")
    test_args <- if (.qisa(self, "CoupledSpinState")) list(0, mi, list(0, 0)) else list(0, mi)
    r0 <- represent(.qfunc(self$type, test_args), basis = basis)$data[[1L]]
    angles <- if (.qisa(self, "Ket")) r0$args[4:6] else r0$args[[1L]]$args[4:6]
    if (all(vapply(angles, function(z) .qeq(z, S_Zero), NA))) return(self)
    state <- .qfunc(evect, c(list(j, mi), args))
    lt <- WignerD(j, mi, .qa(self, "m"), angles[[1L]], angles[[2L]], angles[[3L]])
    Sum(.mul(list(lt, state)), list(mi, .neg(j), j))
  },
  `_eval_innerproduct_JxBra` = function(self, bra, ...) .qspin_ip(self, bra, "_represent_JxOp"),
  `_eval_innerproduct_JyBra` = function(self, bra, ...) .qspin_ip(self, bra, "_represent_JyOp"),
  `_eval_innerproduct_JzBra` = function(self, bra, ...) .qspin_ip(self, bra, "_represent_JzOp"),
  `_eval_trace` = function(self, bra, ...) doit(.mul(list(bra, self))))

.qspin_ip <- function(self, bra, rep) {
  result <- KroneckerDelta(.qa(self, "j"), .qa(bra, "j"))
  if (!identical(.qmeth(bra$type, "dual_class"), self$type)) {
    v <- .qmeth(self, rep, NULL)
    result <- .mul(list(result, v$data[[.qint(.sub(.qa(bra, "j"), .qa(bra, "m"))) + 1]]))
  } else {
    result <- .mul(list(result, KroneckerDelta(.qa(self, "j"), .qa(bra, "j")), KroneckerDelta(.qa(self, "m"), .qa(bra, "m"))))
  }
  result
}

.q3pi2 <- function() .mul(list(sym_pi, as_sym(gmp::as.bigq(3, 2))))
.qpi2 <- function() .mul(list(sym_pi, S_Half))

.qclass("JxKet", c("SpinState", "Ket"),
  dual_class = .cm(function(cls) "JxBra"), coupled_class = .cm(function(cls) "JxKetCoupled"),
  `_represent_default_basis` = function(self, ...) .qmeth(self, "_represent_JxOp", NULL, ...),
  `_represent_JxOp` = function(self, basis, ...) .qmeth(self, "_represent_base", ...),
  `_represent_JyOp` = function(self, basis, ...) .qmeth(self, "_represent_base", alpha = .q3pi2(), ...),
  `_represent_JzOp` = function(self, basis, ...) .qmeth(self, "_represent_base", beta = .qpi2(), ...))
.qclass("JxBra", c("SpinState", "Bra"),
  dual_class = .cm(function(cls) "JxKet"), coupled_class = .cm(function(cls) "JxBraCoupled"))
.qclass("JyKet", c("SpinState", "Ket"),
  dual_class = .cm(function(cls) "JyBra"), coupled_class = .cm(function(cls) "JyKetCoupled"),
  `_represent_default_basis` = function(self, ...) .qmeth(self, "_represent_JyOp", NULL, ...),
  `_represent_JxOp` = function(self, basis, ...) .qmeth(self, "_represent_base", gamma = .qpi2(), ...),
  `_represent_JyOp` = function(self, basis, ...) .qmeth(self, "_represent_base", ...),
  `_represent_JzOp` = function(self, basis, ...)
    .qmeth(self, "_represent_base", alpha = .q3pi2(), beta = .neg(.qpi2()), gamma = .qpi2(), ...))
.qclass("JyBra", c("SpinState", "Bra"),
  dual_class = .cm(function(cls) "JyKet"), coupled_class = .cm(function(cls) "JyBraCoupled"))
.qclass("JzKet", c("SpinState", "Ket"),
  dual_class = .cm(function(cls) "JzBra"), coupled_class = .cm(function(cls) "JzKetCoupled"),
  `_represent_default_basis` = function(self, ...) .qmeth(self, "_represent_JzOp", NULL, ...),
  `_represent_JxOp` = function(self, basis, ...) .qmeth(self, "_represent_base", beta = .q3pi2(), ...),
  `_represent_JyOp` = function(self, basis, ...)
    .qmeth(self, "_represent_base", alpha = .q3pi2(), beta = .qpi2(), gamma = .qpi2(), ...),
  `_represent_JzOp` = function(self, basis, ...) .qmeth(self, "_represent_base", ...))
.qclass("JzBra", c("SpinState", "Bra"),
  dual_class = .cm(function(cls) "JzKet"), coupled_class = .cm(function(cls) "JzBraCoupled"))

# ---------------------------------------------------------------------------
# Coupled spin states

.qbuild_coupled <- function(jcoupling, length) {
  n_list <- lapply(seq_len(length), function(n) n)
  coupled_jn <- list(); coupled_n <- list()
  for (t in jcoupling) {
    n1 <- .qint(t[[1L]]); n2 <- .qint(t[[2L]])
    coupled_jn[[length(coupled_jn) + 1L]] <- t[[3L]]
    coupled_n[[length(coupled_n) + 1L]] <- list(n_list[[n1]], n_list[[n2]])
    n_sort <- sort(c(n_list[[n1]], n_list[[n2]]))
    n_list[[n_sort[1L]]] <- n_sort
  }
  list(coupled_n, coupled_jn)
}

.qpyclassname <- function(x) {
  if (is_sym(x)) {
    if (x$type == "Integer") return("int")
    return(.qclassname(x))
  }
  if (is.character(x)) return("str")
  if (is.numeric(x)) return(if (x == round(x)) "int" else "float")
  class(x)[1L]
}

.qtuple_args <- function(x) if (is_sym(x) && x$type == "Tuple") x$args else x

.qclass("CoupledSpinState", "SpinState",
  `__new__` = .cm(function(cls, args, ...) {
    if (length(args) < 3L) .qstop("TypeError", "CoupledSpinState.__new__() missing required argument: 'jn'")
    j <- args[[1L]]; m <- args[[2L]]; jn <- args[[3L]]; jcoupling <- args[-(1:3)]
    SpinState(j, m)
    is_seq <- function(x) (is.list(x) && !is_sym(x)) || (is_sym(x) && x$type == "Tuple")
    if (!length(jcoupling)) {
      if (!is_seq(jn)) .qstop("TypeError", "object of type '", .qpyclassname(jn), "' has no len()")
      jnl <- .qtuple_args(jn)
      jcoupling <- list()
      if (length(jnl) >= 3L) for (n in 2:(length(jnl) - 1L))
        jcoupling[[length(jcoupling) + 1L]] <- list(S_One, as_sym(n), .add(lapply(jnl[seq_len(n)], .qS)))
      jcoupling[[length(jcoupling) + 1L]] <- list(S_One, as_sym(length(jnl)), .qS(j))
    } else if (length(jcoupling) == 1L) {
      jcoupling <- jcoupling[[1L]]
    } else {
      .qstop("TypeError", "CoupledSpinState only takes 3 or 4 arguments, got: ", length(jcoupling) + 3L)
    }
    if (!is_seq(jn)) .qstop("TypeError", "jn must be Tuple, list or tuple, got ", .qpyclassname(jn))
    if (!is_seq(jcoupling)) .qstop("TypeError", "jcoupling must be Tuple, list or tuple, got ", .qpyclassname(jcoupling))
    jcoupling <- .qtuple_args(jcoupling)
    if (!all(vapply(jcoupling, is_seq, NA))) .qstop("TypeError", "All elements of jcoupling must be list, tuple or Tuple")
    jnl <- .qtuple_args(jn)
    if (length(jnl) - 1L != length(jcoupling))
      .qstop("ValueError", "jcoupling must have length of ", length(jnl) - 1L, ", got ", length(jcoupling))
    if (!all(vapply(jcoupling, function(x) length(.qtuple_args(x)) == 3L, NA)))
      .qstop("ValueError", "All elements of jcoupling must have length 3")
    j <- .qS(j); m <- .qS(m)
    jnl <- lapply(jnl, .qS)
    jcl <- lapply(jcoupling, function(t) lapply(.qtuple_args(t), .qS))
    if (any(vapply(jnl, function(ji) .expr_is_number(ji) && .qhalfint_bad(ji), NA)))
      .qstop("ValueError", "All elements of jn must be integer or half-integer, got: ", .str(do.call(Tuple, jnl)))
    if (any(vapply(jcl, function(t) t[[1L]]$type != "Integer" || t[[2L]]$type != "Integer", NA)))
      .qstop("ValueError", "Indices in jcoupling must be integers")
    if (any(vapply(jcl, function(t) .qint(t[[1L]]) < 1 || .qint(t[[2L]]) < 1 || .qint(t[[1L]]) > length(jnl) ||
                     .qint(t[[2L]]) > length(jnl), NA)))
      .qstop("ValueError", "Indices must be between 1 and the number of coupled spin spaces")
    if (any(vapply(jcl, function(t) .expr_is_number(t[[3L]]) && .qhalfint_bad(t[[3L]]), NA)))
      .qstop("ValueError", "All coupled j values in coupling scheme must be integer or half-integer")
    bc <- .qbuild_coupled(jcl, length(jnl))
    coupled_n <- bc[[1L]]; coupled_jn <- bc[[2L]]
    jvals <- jnl
    for (n in seq_along(coupled_n)) {
      n1 <- coupled_n[[n]][[1L]]; n2 <- coupled_n[[n]][[2L]]
      j1 <- jvals[[min(n1)]]; j2 <- jvals[[min(n2)]]; j3 <- coupled_jn[[n]]
      if (.expr_is_number(j1) && .expr_is_number(j2) && .expr_is_number(j3)) {
        if (isTRUE(.qgt(j3, .add(list(j1, j2)))))
          .qstop("ValueError", sprintf("All couplings must have j1+j2 >= j3, in coupling number %d got j1,j2,j3: %d,%d,%d",
                                       n, trunc(.qnum(j1)), trunc(.qnum(j2)), trunc(.qnum(j3))))
        if (isTRUE(.qgt(abs(.sub(j1, j2)), j3)))
          .qstop("ValueError", sprintf("All couplings must have |j1+j2| <= j3, in coupling number %d got j1,j2,j3: %d,%d,%d",
                                       n, trunc(.qnum(j1)), trunc(.qnum(j2)), trunc(.qnum(j3))))
      }
      jvals[[min(c(n1, n2))]] <- j3
    }
    if (length(jcl) && !.qeq(jcl[[length(jcl)]][[3L]], j))
      .qstop("ValueError", "Last j value coupled together must be the final j of the state")
    .qsuper("State", "__new__", cls, list(j, m, do.call(Tuple, jnl), do.call(Tuple, lapply(jcl, function(t) do.call(Tuple, t)))))
  }),
  `_print_label` = function(self, p) {
    label <- c(.qpr(p, .qa(self, "j")), .qpr(p, .qa(self, "m")))
    jn <- .qa(self, "jn")$args
    for (i in seq_along(jn)) label <- c(label, sprintf("j%d=%s", i, .qpr(p, jn[[i]])))
    cjn <- .qa(self, "coupled_jn"); cn <- .qa(self, "coupled_n")
    for (k in seq_len(length(cjn) - 1L)) {
      label <- c(label, sprintf("j(%s)=%s", paste(sort(c(cn[[k]][[1L]], cn[[k]][[2L]])), collapse = ","), .qpr(p, cjn[[k]])))
    }
    paste(label, collapse = ",")
  },
  `_print_label_latex` = function(self, p) {
    label <- c(.qpr(p, .qa(self, "j")), .qpr(p, .qa(self, "m")))
    jn <- .qa(self, "jn")$args
    for (i in seq_along(jn)) label <- c(label, sprintf("j_{%d}=%s", i, .qpr(p, jn[[i]])))
    cjn <- .qa(self, "coupled_jn"); cn <- .qa(self, "coupled_n")
    for (k in seq_len(length(cjn) - 1L)) {
      label <- c(label, sprintf("j_{%s}=%s", paste(sort(c(cn[[k]][[1L]], cn[[k]][[2L]])), collapse = ","), .qpr(p, cjn[[k]])))
    }
    paste(label, collapse = .qa(self, "_label_separator"))
  },
  jn = .prop(function(self) .qa(self, "label")[[3L]]),
  coupling = .prop(function(self) .qa(self, "label")[[4L]]),
  coupled_jn = .prop(function(self) .qbuild_coupled(lapply(.qa(self, "label")[[4L]]$args, function(t) t$args), length(.qa(self, "label")[[3L]]$args))[[2L]]),
  coupled_n = .prop(function(self) .qbuild_coupled(lapply(.qa(self, "label")[[4L]]$args, function(t) t$args), length(.qa(self, "label")[[3L]]$args))[[1L]]),
  `_eval_hilbert_space` = .cm(function(cls, label) {
    j <- .add(label[[3L]]$args)
    if (.expr_is_number(j)) {
      top <- .qint(.add(list(.mul(list(S_Two, j)), S_One)))
      return(do.call(DirectSumHilbertSpace, lapply(seq(top, 1, by = -2), function(x) ComplexSpace(x))))
    }
    ComplexSpace(.add(list(.mul(list(S_Two, j)), S_One)))
  }),
  `_represent_coupled_base` = function(self, ...) {
    evect <- .qmeth(self$type, "uncoupled_class")
    j <- .qa(self, "j")
    if (!.expr_is_number(j)) .qstop("ValueError", "State must not have symbolic j value to represent")
    dim <- .qa(.qhs(self), "dimension")
    if (!.expr_is_number(dim)) .qstop("ValueError", "State must not have symbolic j values to represent")
    jv <- as.numeric(N(j))
    start <- if (jv == round(jv)) jv^2 else (2 * jv - 1) * (1 + 2 * jv) / 4
    vals <- rep(list(S_Zero), .qint(dim))
    sub <- .qmeth(.qfunc(evect, list(j, .qa(self, "m"))), "_represent_base", ...)
    for (i in seq_along(sub$data)) vals[[start + i]] <- sub$data[[i]]
    .new_mat(vals, length(vals), 1L)
  },
  `_eval_rewrite_as_Jx` = function(self, args, ...) {
    if (.qisa(self, "Bra")) return(.qmeth(self, "_rewrite_basis", JxOp("J"), "JxBraCoupled", ...))
    .qmeth(self, "_rewrite_basis", JxOp("J"), "JxKetCoupled", ...)
  },
  `_eval_rewrite_as_Jy` = function(self, args, ...) {
    if (.qisa(self, "Bra")) return(.qmeth(self, "_rewrite_basis", JyOp("J"), "JyBraCoupled", ...))
    .qmeth(self, "_rewrite_basis", JyOp("J"), "JyKetCoupled", ...)
  },
  `_eval_rewrite_as_Jz` = function(self, args, ...) {
    if (.qisa(self, "Bra")) return(.qmeth(self, "_rewrite_basis", JzOp("J"), "JzBraCoupled", ...))
    .qmeth(self, "_rewrite_basis", JzOp("J"), "JzKetCoupled", ...)
  })

.qnum <- function(x) as.numeric(N(x))

.qclass("JxKetCoupled", c("CoupledSpinState", "Ket"),
  dual_class = .cm(function(cls) "JxBraCoupled"), uncoupled_class = .cm(function(cls) "JxKet"),
  `_represent_default_basis` = function(self, ...) .qmeth(self, "_represent_JzOp", NULL, ...),
  `_represent_JxOp` = function(self, basis, ...) .qmeth(self, "_represent_coupled_base", ...),
  `_represent_JyOp` = function(self, basis, ...) .qmeth(self, "_represent_coupled_base", alpha = .q3pi2(), ...),
  `_represent_JzOp` = function(self, basis, ...) .qmeth(self, "_represent_coupled_base", beta = .qpi2(), ...))
.qclass("JxBraCoupled", c("CoupledSpinState", "Bra"),
  dual_class = .cm(function(cls) "JxKetCoupled"), uncoupled_class = .cm(function(cls) "JxBra"))
.qclass("JyKetCoupled", c("CoupledSpinState", "Ket"),
  dual_class = .cm(function(cls) "JyBraCoupled"), uncoupled_class = .cm(function(cls) "JyKet"),
  `_represent_default_basis` = function(self, ...) .qmeth(self, "_represent_JzOp", NULL, ...),
  `_represent_JxOp` = function(self, basis, ...) .qmeth(self, "_represent_coupled_base", gamma = .qpi2(), ...),
  `_represent_JyOp` = function(self, basis, ...) .qmeth(self, "_represent_coupled_base", ...),
  `_represent_JzOp` = function(self, basis, ...)
    .qmeth(self, "_represent_coupled_base", alpha = .q3pi2(), beta = .neg(.qpi2()), gamma = .qpi2(), ...))
.qclass("JyBraCoupled", c("CoupledSpinState", "Bra"),
  dual_class = .cm(function(cls) "JyKetCoupled"), uncoupled_class = .cm(function(cls) "JyBra"))
.qclass("JzKetCoupled", c("CoupledSpinState", "Ket"),
  dual_class = .cm(function(cls) "JzBraCoupled"), uncoupled_class = .cm(function(cls) "JzKet"),
  `_represent_default_basis` = function(self, ...) .qmeth(self, "_represent_JzOp", NULL, ...),
  `_represent_JxOp` = function(self, basis, ...) .qmeth(self, "_represent_coupled_base", beta = .q3pi2(), ...),
  `_represent_JyOp` = function(self, basis, ...)
    .qmeth(self, "_represent_coupled_base", alpha = .q3pi2(), beta = .qpi2(), gamma = .qpi2(), ...),
  `_represent_JzOp` = function(self, basis, ...) .qmeth(self, "_represent_coupled_base", ...))
.qclass("JzBraCoupled", c("CoupledSpinState", "Bra"),
  dual_class = .cm(function(cls) "JzKetCoupled"), uncoupled_class = .cm(function(cls) "JzBra"))

#' Spin states
#'
#' The eigenstates of `sympy.physics.quantum.spin`: `JzKet(j, m)`,
#' `JxKet()`, `JyKet()` and their bras, and the coupled states
#' `JzKetCoupled(j, m, jn, jcoupling)` (and the `Jx`, `Jy` and bra
#' versions) of several spin spaces with momenta `jn`. `couple()` writes
#' tensor products of uncoupled states in the coupled basis and
#' `uncouple()` the reverse. [rewrite()] with `"Jx"`, `"Jy"` or `"Jz"`
#' changes the basis.
#'
#' @param j,m The angular momentum and its projection.
#' @param ... `jn` and the optional `jcoupling` of coupled states.
#' @param expr An expression of spin states.
#' @param jcoupling_list The order in which spaces are coupled.
#' @param jn The momenta of the spaces, for uncoupling an uncoupled state.
#' @return A spin state or an expression.
#' @examples
#' JzKet(1, 0)
#' couple(TensorProduct(JzKet(S(1)/2, S(1)/2), JzKet(S(1)/2, S(1)/2)))
#' uncouple(JzKetCoupled(1, 0, list(S(1)/2, S(1)/2)))
#' @name spin-states
NULL

#' @rdname spin-states
#' @export
SpinState <- function(j, m) .qfunc("SpinState", list(j, m))
#' @rdname spin-states
#' @export
JxKet <- function(j, m) .qfunc("JxKet", list(j, m))
#' @rdname spin-states
#' @export
JxBra <- function(j, m) .qfunc("JxBra", list(j, m))
#' @rdname spin-states
#' @export
JyKet <- function(j, m) .qfunc("JyKet", list(j, m))
#' @rdname spin-states
#' @export
JyBra <- function(j, m) .qfunc("JyBra", list(j, m))
#' @rdname spin-states
#' @export
JzKet <- function(j, m) .qfunc("JzKet", list(j, m))
#' @rdname spin-states
#' @export
JzBra <- function(j, m) .qfunc("JzBra", list(j, m))
#' @rdname spin-states
#' @export
CoupledSpinState <- function(j, m, ...) .qfunc("CoupledSpinState", list(j, m, ...))
#' @rdname spin-states
#' @export
JxKetCoupled <- function(j, m, ...) .qfunc("JxKetCoupled", list(j, m, ...))
#' @rdname spin-states
#' @export
JxBraCoupled <- function(j, m, ...) .qfunc("JxBraCoupled", list(j, m, ...))
#' @rdname spin-states
#' @export
JyKetCoupled <- function(j, m, ...) .qfunc("JyKetCoupled", list(j, m, ...))
#' @rdname spin-states
#' @export
JyBraCoupled <- function(j, m, ...) .qfunc("JyBraCoupled", list(j, m, ...))
#' @rdname spin-states
#' @export
JzKetCoupled <- function(j, m, ...) .qfunc("JzKetCoupled", list(j, m, ...))
#' @rdname spin-states
#' @export
JzBraCoupled <- function(j, m, ...) .qfunc("JzBraCoupled", list(j, m, ...))

# all nodes of `e` that are instances of `cls`, in traversal order
.qatoms_of <- function(e, cls) {
  out <- list(); seen <- character()
  walk <- function(z) {
    if (!is_sym(z)) return()
    if (.qisa(z, cls) && !(z$key %in% seen)) { seen <<- c(seen, z$key); out[[length(out) + 1L]] <<- z }
    for (a in z$args) walk(a)
  }
  walk(e)
  out
}

#' @rdname spin-states
#' @export
couple <- function(expr, jcoupling_list = NULL) {
  for (tp in .qatoms_of(expr, "QTensorProduct")) {
    if (!all(vapply(tp$args, function(s) .qisa(s, "SpinState"), NA))) next
    if (!all(vapply(tp$args, function(s) s$type == tp$args[[1L]]$type, NA))) .qstop("TypeError", "All states must be the same basis")
    expr <- subs(expr, tp, .qcouple(tp, jcoupling_list))
  }
  expr
}

.qconfignum_to_difflist <- function(config_num, diff, list_len) {
  diff_list <- numeric()
  bin <- function(n, k) if (k < 0) 0 else if (n < 0) { if (k == 0) 1 else 0 } else choose(n, k)
  for (n in seq_len(list_len) - 1L) {
    prev_diff <- diff
    rem_spots <- list_len - n - 1
    rem_configs <- bin(diff + rem_spots - 1, diff)
    while (config_num >= rem_configs) {
      config_num <- config_num - rem_configs
      diff <- diff - 1
      rem_configs <- bin(diff + rem_spots - 1, diff)
    }
    diff_list <- c(diff_list, prev_diff - diff)
  }
  diff_list
}

.qcouple <- function(tp, jcoupling_list) {
  states <- tp$args
  coupled_evect <- .qmeth(states[[1L]]$type, "coupled_class")
  if (is.null(jcoupling_list)) {
    jcoupling_list <- list()
    for (n in seq_len(length(states) - 1L)) jcoupling_list[[n]] <- list(1, n + 1)
  }
  if (length(jcoupling_list) != length(states) - 1L)
    .qstop("TypeError", "jcoupling_list must be length ", length(states) - 1L, ", got ", length(jcoupling_list))
  if (!all(vapply(jcoupling_list, function(cp) length(cp) == 2L, NA))) .qstop("ValueError", "Each coupling must define 2 spaces")
  jl <- lapply(jcoupling_list, function(cp) vapply(cp, .qint, 0))
  if (any(vapply(jl, function(cp) cp[1L] == cp[2L], NA))) .qstop("ValueError", "Spin spaces cannot couple to themselves")
  j_test <- rep(0, length(states))
  for (cp in jl) {
    if (j_test[cp[1L]] == -1 || j_test[cp[2L]] == -1) .qstop("ValueError", "Spaces coupling j_n's are referenced by smallest n value")
    j_test[max(cp)] <- -1
  }
  jn <- lapply(states, function(s) .qa(s, "j"))
  mn <- lapply(states, function(s) .qa(s, "m"))
  coupling_list <- list()
  n_list <- lapply(seq_along(states), function(i) i)
  for (cp in jl) {
    j1_n <- n_list[[cp[1L]]]; j2_n <- n_list[[cp[2L]]]
    coupling_list[[length(coupling_list) + 1L]] <- list(j1_n, j2_n)
    n_list[[min(cp)]] <- sort(c(j1_n, j2_n))
  }
  if (all(vapply(states, function(s) .expr_is_number(.qa(s, "j")) && .expr_is_number(.qa(s, "m")), NA))) {
    diff_max <- lapply(coupling_list, function(cp) .add(lapply(c(cp[[1L]], cp[[2L]]), function(n) .sub(jn[[n]], mn[[n]]))))
    dmax <- vapply(diff_max, .qnum, 0)
    result <- list()
    for (diff in seq(0, dmax[length(dmax)])) {
      n <- length(coupling_list)
      tot <- choose(diff + n - 1, diff)
      for (config_num in seq_len(tot) - 1) {
        diff_list <- .qconfignum_to_difflist(config_num, diff, n)
        if (any(diff_list > dmax)) next
        cg_terms <- list(); coupled_j <- jn; jcoupling <- list()
        for (k in seq_along(coupling_list)) {
          j1_n <- coupling_list[[k]][[1L]]; j2_n <- coupling_list[[k]][[2L]]
          j1 <- coupled_j[[min(j1_n)]]; j2 <- coupled_j[[min(j2_n)]]
          j3 <- .sub(.add(list(j1, j2)), as_sym(diff_list[k]))
          coupled_j[[min(c(j1_n, j2_n))]] <- j3
          m1 <- .add(mn[j1_n]); m2 <- .add(mn[j2_n]); m3 <- .add(list(m1, m2))
          cg_terms[[k]] <- list(j1, m1, j2, m2, j3, m3)
          jcoupling[[k]] <- list(as_sym(min(j1_n)), as_sym(min(j2_n)), j3)
        }
        if (any(vapply(cg_terms, function(t) isTRUE(.qgt(abs(t[[6L]]), t[[5L]])), NA))) next
        if (any(vapply(cg_terms, function(t) isTRUE(.qgt(t[[5L]], .add(list(t[[1L]], t[[3L]])))), NA))) next
        if (any(vapply(cg_terms, function(t) isTRUE(.qgt(abs(.sub(t[[1L]], t[[3L]])), t[[5L]])), NA))) next
        coeff <- .mul(lapply(cg_terms, function(t) doit(do.call(CG, t))))
        state <- .qfunc(coupled_evect, list(j3, m3, jn, jcoupling))
        result[[length(result) + 1L]] <- .mul(list(coeff, state))
      }
    }
    return(.add(result))
  }
  cg_terms <- list(); jcoupling <- list(); sum_terms <- list(); coupled_j <- jn
  for (cp in coupling_list) {
    j1_n <- cp[[1L]]; j2_n <- cp[[2L]]
    j1 <- coupled_j[[min(j1_n)]]; j2 <- coupled_j[[min(j2_n)]]
    j3 <- if (length(c(j1_n, j2_n)) == length(states)) Symbol("j") else Symbol(paste0("j", paste(c(j1_n, j2_n), collapse = "")))
    coupled_j[[min(c(j1_n, j2_n))]] <- j3
    m1 <- .add(mn[j1_n]); m2 <- .add(mn[j2_n]); m3 <- .add(list(m1, m2))
    cg_terms[[length(cg_terms) + 1L]] <- list(j1, m1, j2, m2, j3, m3)
    jcoupling[[length(jcoupling) + 1L]] <- list(as_sym(min(j1_n)), as_sym(min(j2_n)), j3)
    sum_terms[[length(sum_terms) + 1L]] <- list(j3, m3, .add(list(j1, j2)))
  }
  coeff <- .mul(lapply(cg_terms, function(t) do.call(CG, t)))
  state <- .qfunc(coupled_evect, list(j3, m3, jn, jcoupling))
  do.call(Sum, c(list(.mul(list(coeff, state))), sum_terms))
}

#' @rdname spin-states
#' @export
uncouple <- function(expr, jn = NULL, jcoupling_list = NULL) {
  for (state in .qatoms_of(expr, "SpinState")) expr <- subs(expr, state, .quncouple(state, jn, jcoupling_list))
  expr
}

.quncouple <- function(state, jn, jcoupling_list) {
  if (.qisa(state, "CoupledSpinState")) {
    jn <- .qa(state, "jn")$args
    coupled_n <- .qa(state, "coupled_n"); coupled_jn <- .qa(state, "coupled_jn")
    evect <- .qmeth(state$type, "uncoupled_class")
  } else if (.qisa(state, "SpinState")) {
    if (is.null(jn)) .qstop("ValueError", "Must specify j-values for coupled state")
    if (!(is.list(jn) && !is_sym(jn))) .qstop("TypeError", "jn must be list or tuple")
    jn <- lapply(jn, .qS)
    if (is.null(jcoupling_list)) {
      jcoupling_list <- list()
      for (i in seq_len(length(jn) - 1L)) jcoupling_list[[i]] <- list(S_One, as_sym(1L + i), .add(jn[seq_len(i + 1L)]))
    }
    if (!(is.list(jcoupling_list) && !is_sym(jcoupling_list))) .qstop("TypeError", "jcoupling must be a list or tuple")
    if (length(jcoupling_list) != length(jn) - 1L) .qstop("ValueError", "Must specify 2 fewer coupling terms than the number of j values")
    bc <- .qbuild_coupled(lapply(jcoupling_list, function(t) lapply(t, .qS)), length(jn))
    coupled_n <- bc[[1L]]; coupled_jn <- bc[[2L]]
    evect <- state$type
  } else {
    .qstop("TypeError", "state must be a spin state")
  }
  j <- .qa(state, "j"); m <- .qa(state, "m")
  coupling_list <- list(); j_list <- jn
  for (k in seq_along(coupled_jn)) {
    n1 <- coupled_n[[k]][[1L]]; n2 <- coupled_n[[k]][[2L]]
    j1 <- j_list[[n1[1L]]]; j2 <- j_list[[n2[1L]]]
    coupling_list[[k]] <- list(n1, n2, j1, j2, coupled_jn[[k]])
    j_list[[min(c(n1, n2))]] <- coupled_jn[[k]]
  }
  if (.expr_is_number(j) && .expr_is_number(m)) {
    diff_max <- vapply(jn, function(x) 2 * .qnum(x), 0)
    diff <- .qnum(.sub(.add(jn), m))
    n <- length(jn)
    tot <- choose(diff + n - 1, diff)
    result <- list()
    for (config_num in seq_len(tot) - 1) {
      diff_list <- .qconfignum_to_difflist(config_num, diff, n)
      if (any(diff_list > diff_max)) next
      cg_terms <- lapply(coupling_list, function(cp) {
        m1 <- .add(lapply(cp[[1L]], function(x) .sub(jn[[x]], as_sym(diff_list[x]))))
        m2 <- .add(lapply(cp[[2L]], function(x) .sub(jn[[x]], as_sym(diff_list[x]))))
        list(cp[[3L]], m1, cp[[4L]], m2, cp[[5L]], .add(list(m1, m2)))
      })
      coeff <- .mul(lapply(cg_terms, function(t) doit(do.call(CG, t))))
      st <- .qtp(lapply(seq_along(jn), function(k) .qfunc(evect, list(jn[[k]], .sub(jn[[k]], as_sym(diff_list[k]))))))
      result[[length(result) + 1L]] <- .mul(list(coeff, st))
    }
    return(.add(result))
  }
  mvals <- lapply(seq_along(jn), function(i) Symbol(paste0("m", i)))
  nc <- length(coupling_list)
  cg_terms <- lapply(coupling_list[seq_len(nc - 1L)], function(cp)
    list(cp[[3L]], .add(mvals[cp[[1L]]]), cp[[4L]], .add(mvals[cp[[2L]]]), cp[[5L]], .add(mvals[c(cp[[1L]], cp[[2L]])])))
  cp <- coupling_list[[nc]]
  cg_terms[[nc]] <- list(cp[[3L]], .add(mvals[cp[[1L]]]), cp[[4L]], .add(mvals[cp[[2L]]]), j, m)
  cg_coeff <- .mul(lapply(cg_terms, function(t) do.call(CG, t)))
  sum_terms <- lapply(seq_along(jn), function(k) list(mvals[[k]], .neg(jn[[k]]), jn[[k]]))
  st <- .qtp(lapply(seq_along(jn), function(k) .qfunc(evect, list(jn[[k]], mvals[[k]]))))
  do.call(Sum, c(list(.mul(list(cg_coeff, st))), sum_terms))
}

# ---------------------------------------------------------------------------
# Basic.rewrite for quantum objects: _eval_rewrite_as_<target> methods,
# applied bottom up with the rewritten arguments

.qrewrite <- function(e, target, ...) {
  if (inherits(e, "symmatrix")) return(.mat_map(e, function(z) .qrewrite(z, target, ...)))
  if (!is_sym(e)) return(e)
  args <- lapply(e$args, function(a) .qrewrite(a, target, ...))
  meth <- paste0("_eval_rewrite_as_", target)
  if (!is.null(.qnode_def(e$type)) && .qhas(e, meth)) {
    r <- .qmeth(e, meth, args, ...)
    if (!is.null(r)) return(r)
  } else if (!is.null(.qnode_def(e$type)) && .qhas(e, "_eval_rewrite")) {
    r <- .qmeth(e, "_eval_rewrite", target, args, ...)
    if (!is.null(r)) return(r)
  } else if (e$type == "Function") {
    rules <- .rewrite_rules[[target]]
    if (!is.null(rules[[e$name]])) {
      r <- do.call(rules[[e$name]], args)
      if (!is.null(r)) return(r)
    }
  }
  if (!length(args)) return(e)
  if (identical(vapply(args, .k, ""), vapply(e$args, .k, ""))) return(e)
  sym_rebuild(e, args)
}

# TensorProduct._eval_rewrite: the rewritten factors, expanded
.qdef("QTensorProduct", `_eval_rewrite` = function(self, rule, args, ...)
  expand(.qfunc("QTensorProduct", args), tensorproduct = TRUE))

# does `e` contain a node of a quantum class
.qhas_node <- function(e) {
  if (!is_sym(e)) return(FALSE)
  if (!is.null(.qnode_def(e$type))) return(TRUE)
  for (a in e$args) if (.qhas_node(a)) return(TRUE)
  FALSE
}
