# Pauli operators and the one-dimensional harmonic oscillator, after
# sympy/physics/quantum/pauli.py and sho1d.py.

.qsigma_names <- c(SigmaX = "x", SigmaY = "y", SigmaZ = "z", SigmaMinus = "-", SigmaPlus = "+")

.qsigma_same <- function(self, other) .qeq(.qa(self, "name"), .qa(other, "name"))

.qsigma_mat <- function(rows, format) {
  if (!identical(format, "sympy")) .qni("Representation in format ", format, " not implemented.")
  Matrix(rows)
}

.qclass("SigmaOpBase", "Operator",
  name = .prop(function(self) self$args[[1L]]),
  use_name = .prop(function(self) .qtruthy(self$args[[1L]])),
  default_args = .cm(function(cls) list(FALSE)),
  `_eval_commutator_BosonOp` = function(self, other, ...) S_Zero,
  `_eval_adjoint` = function(self) self,
  `_print_contents` = function(self, p) paste0(.qpyname(self$type), "()"),
  `_print_contents_latex` = function(self, p) {
    s <- .qsigma_names[[self$type]]
    if (.qa(self, "use_name")) sprintf("{\\sigma_%s^{(%s)}}", s, .str(.qa(self, "name"))) else sprintf("{\\sigma_%s}", s)
  },
  `_eval_power` = function(self, e) {
    if (e$type == "Integer" && .qint(e) > 0) {
      if (.qisa(self, c("SigmaMinus", "SigmaPlus"))) return(S_Zero)
      return(.pow(.qfunc(self$type, list(.qa(self, "name"))), as_sym(.qint(e) %% 2)))
    }
    NULL
  })

.qsigma_comm <- function(res) function(self, other, ...) if (!.qsigma_same(self, other)) S_Zero else res(self)

.qclass("SigmaX", "SigmaOpBase",
  `_eval_commutator_SigmaY` = .qsigma_comm(function(s) .mul(list(S_Two, sym_I, SigmaZ(.qa(s, "name"))))),
  `_eval_commutator_SigmaZ` = .qsigma_comm(function(s) .mul(list(as_sym(-2L), sym_I, SigmaY(.qa(s, "name"))))),
  `_eval_anticommutator_SigmaY` = function(self, other, ...) S_Zero,
  `_eval_anticommutator_SigmaZ` = function(self, other, ...) S_Zero,
  `_represent_default_basis` = function(self, format = "sympy", ...) .qsigma_mat(list(list(0, 1), list(1, 0)), format))

.qclass("SigmaY", "SigmaOpBase",
  `_eval_commutator_SigmaZ` = .qsigma_comm(function(s) .mul(list(S_Two, sym_I, SigmaX(.qa(s, "name"))))),
  `_eval_commutator_SigmaX` = .qsigma_comm(function(s) .mul(list(as_sym(-2L), sym_I, SigmaZ(.qa(s, "name"))))),
  `_eval_anticommutator_SigmaX` = function(self, other, ...) S_Zero,
  `_eval_anticommutator_SigmaZ` = function(self, other, ...) S_Zero,
  `_represent_default_basis` = function(self, format = "sympy", ...)
    .qsigma_mat(list(list(0, .neg(sym_I)), list(sym_I, 0)), format))

.qclass("SigmaZ", "SigmaOpBase",
  `_eval_commutator_SigmaX` = .qsigma_comm(function(s) .mul(list(S_Two, sym_I, SigmaY(.qa(s, "name"))))),
  `_eval_commutator_SigmaY` = .qsigma_comm(function(s) .mul(list(as_sym(-2L), sym_I, SigmaX(.qa(s, "name"))))),
  `_eval_anticommutator_SigmaX` = function(self, other, ...) S_Zero,
  `_eval_anticommutator_SigmaY` = function(self, other, ...) S_Zero,
  `_represent_default_basis` = function(self, format = "sympy", ...) .qsigma_mat(list(list(1, 0), list(0, -1)), format))

.qclass("SigmaMinus", "SigmaOpBase",
  is_annihilation = TRUE,
  `_eval_commutator_SigmaX` = .qsigma_comm(function(s) .neg(SigmaZ(.qa(s, "name")))),
  `_eval_commutator_SigmaY` = .qsigma_comm(function(s) .mul(list(sym_I, SigmaZ(.qa(s, "name"))))),
  `_eval_commutator_SigmaZ` = function(self, other, ...) .mul(list(S_Two, self)),
  `_eval_commutator_SigmaMinus` = function(self, other, ...) SigmaZ(.qa(self, "name")),
  `_eval_anticommutator_SigmaZ` = function(self, other, ...) S_Zero,
  `_eval_anticommutator_SigmaX` = function(self, other, ...) S_One,
  `_eval_anticommutator_SigmaY` = function(self, other, ...) .neg(sym_I),
  `_eval_anticommutator_SigmaPlus` = function(self, other, ...) S_One,
  `_eval_adjoint` = function(self) SigmaPlus(.qa(self, "name")),
  `_represent_default_basis` = function(self, format = "sympy", ...) .qsigma_mat(list(list(0, 0), list(1, 0)), format))

.qclass("SigmaPlus", "SigmaOpBase",
  is_annihilation = FALSE,
  `_eval_commutator_SigmaX` = .qsigma_comm(function(s) SigmaZ(.qa(s, "name"))),
  `_eval_commutator_SigmaY` = .qsigma_comm(function(s) .mul(list(sym_I, SigmaZ(.qa(s, "name"))))),
  `_eval_commutator_SigmaZ` = .qsigma_comm(function(s) .mul(list(as_sym(-2L), s))),
  `_eval_commutator_SigmaMinus` = function(self, other, ...) SigmaZ(.qa(self, "name")),
  `_eval_anticommutator_SigmaZ` = function(self, other, ...) S_Zero,
  `_eval_anticommutator_SigmaX` = function(self, other, ...) S_One,
  `_eval_anticommutator_SigmaY` = function(self, other, ...) sym_I,
  `_eval_anticommutator_SigmaMinus` = function(self, other, ...) S_One,
  `_eval_adjoint` = function(self) SigmaMinus(.qa(self, "name")),
  `_eval_mul` = function(self, other) .mul(list(self, other)),
  `_represent_default_basis` = function(self, format = "sympy", ...) .qsigma_mat(list(list(0, 1), list(0, 0)), format))

.qsigmaz_n <- function(args, cls) {
  if (length(args) != 1L) .qstop("TypeError", cls, ".__new__() takes 2 positional arguments but ", length(args) + 1L, " were given")
  n <- args[[1L]]
  ok <- (is.numeric(n) || is.logical(n)) && as.numeric(n) %in% c(0, 1) ||
    is_sym(n) && n$type == "Integer" && .qint(n) %in% c(0, 1)
  if (!ok) .qstop("ValueError", "n must be 0 or 1")
}

.qclass("SigmaZKet", "Ket",
  `__new__` = .cm(function(cls, args, ...) {
    .qsigmaz_n(args, "SigmaZKet")
    .qsuper("Ket", "__new__", cls, args)
  }),
  n = .prop(function(self) .qa(self, "label")[[1L]]),
  dual_class = .cm(function(cls) "SigmaZBra"),
  `_eval_hilbert_space` = .cm(function(cls, label) ComplexSpace(2)),
  `_eval_innerproduct_SigmaZBra` = function(self, bra, ...) KroneckerDelta(.qa(self, "n"), .qa(bra, "n")),
  `_apply_from_right_to_SigmaZ` = function(self, op, ...) if (.qint(.qa(self, "n")) == 0) self else .neg(self),
  `_apply_from_right_to_SigmaX` = function(self, op, ...) if (.qint(.qa(self, "n")) == 0) SigmaZKet(1) else SigmaZKet(0),
  `_apply_from_right_to_SigmaY` = function(self, op, ...)
    if (.qint(.qa(self, "n")) == 0) .mul(list(sym_I, SigmaZKet(1))) else .mul(list(.neg(sym_I), SigmaZKet(0))),
  `_apply_from_right_to_SigmaMinus` = function(self, op, ...) if (.qint(.qa(self, "n")) == 0) SigmaZKet(1) else S_Zero,
  `_apply_from_right_to_SigmaPlus` = function(self, op, ...) if (.qint(.qa(self, "n")) == 0) S_Zero else SigmaZKet(0),
  `_represent_default_basis` = function(self, format = "sympy", ...)
    .qsigma_mat(if (.qint(.qa(self, "n")) == 0) list(list(1), list(0)) else list(list(0), list(1)), format))

.qclass("SigmaZBra", "Bra",
  `__new__` = .cm(function(cls, args, ...) {
    .qsigmaz_n(args, "SigmaZBra")
    .qsuper("Bra", "__new__", cls, args)
  }),
  n = .prop(function(self) .qa(self, "label")[[1L]]),
  dual_class = .cm(function(cls) "SigmaZKet"))

.qsigma_lt <- function(a, b) {
  r <- .qgt(b, a)
  if (is.na(r)) .qstop("TypeError", "cannot determine truth value of Relational")
  r
}

.qsimplify_pauli_product <- function(a, b) {
  if (!(.qisa(a, "SigmaOpBase") && .qisa(b, "SigmaOpBase"))) return(.mul(list(a, b)))
  an <- .qa(a, "name")
  if (!.qeq(an, .qa(b, "name"))) {
    if (.qsigma_lt(an, .qa(b, "name"))) return(.mul(list(a, b))) else return(.mul(list(b, a)))
  }
  half_z <- function(sg) .add(list(S_Half, .mul(list(as_sym(sg), SigmaZ(an), S_Half))))
  i_half_z <- function(c, sg) .mul(list(c, .add(list(S_One, .mul(list(as_sym(sg), SigmaZ(an))))), S_Half))
  bt <- b$type
  switch(a$type,
    SigmaX = switch(bt, SigmaX = S_One, SigmaY = .mul(list(sym_I, SigmaZ(an))), SigmaZ = .mul(list(.neg(sym_I), SigmaY(an))),
                    SigmaMinus = half_z(1), SigmaPlus = half_z(-1)),
    SigmaY = switch(bt, SigmaX = .mul(list(.neg(sym_I), SigmaZ(an))), SigmaY = S_One, SigmaZ = .mul(list(sym_I, SigmaX(an))),
                    SigmaMinus = i_half_z(.neg(sym_I), 1), SigmaPlus = i_half_z(sym_I, -1)),
    SigmaZ = switch(bt, SigmaX = .mul(list(sym_I, SigmaY(an))), SigmaY = .mul(list(.neg(sym_I), SigmaX(an))), SigmaZ = S_One,
                    SigmaMinus = .neg(SigmaMinus(an)), SigmaPlus = SigmaPlus(an)),
    SigmaMinus = switch(bt, SigmaX = i_half_z(S_One, -1), SigmaY = i_half_z(.neg(sym_I), -1), SigmaZ = SigmaMinus(.qa(b, "name")),
                        SigmaMinus = S_Zero, SigmaPlus = half_z(-1)),
    SigmaPlus = switch(bt, SigmaX = i_half_z(S_One, 1), SigmaY = i_half_z(sym_I, 1), SigmaZ = .neg(SigmaPlus(an)),
                       SigmaMinus = i_half_z(S_One, 1), SigmaPlus = S_Zero),
    .mul(list(a, b)))
}

#' Pauli operators
#'
#' The Pauli operators `SigmaX()`, `SigmaY()`, `SigmaZ()` and the ladder
#' operators `SigmaMinus()` and `SigmaPlus()` of
#' `sympy.physics.quantum.pauli`, labelled by an optional name, and the
#' states `SigmaZKet(n)` and `SigmaZBra(n)` with `n` 0 or 1.
#' `qsimplify_pauli()` simplifies products of Pauli operators with the
#' same name.
#'
#' @param ... The optional name of the operator.
#' @param n 0 or 1.
#' @param e An expression.
#' @return A quantum operator, state or expression.
#' @examples
#' sx <- SigmaX(1); sy <- SigmaY(1)
#' doit(Commutator(sx, sy))
#' qsimplify_pauli(sx * sy)
#' qapply(sx * SigmaZKet(0))
#' @name pauli
NULL

#' @rdname pauli
#' @export
SigmaX <- function(...) .qfunc("SigmaX", list(...))
#' @rdname pauli
#' @export
SigmaY <- function(...) .qfunc("SigmaY", list(...))
#' @rdname pauli
#' @export
SigmaZ <- function(...) .qfunc("SigmaZ", list(...))
#' @rdname pauli
#' @export
SigmaMinus <- function(...) .qfunc("SigmaMinus", list(...))
#' @rdname pauli
#' @export
SigmaPlus <- function(...) .qfunc("SigmaPlus", list(...))
#' @rdname pauli
#' @export
SigmaZKet <- function(n) .qfunc("SigmaZKet", list(n))
#' @rdname pauli
#' @export
SigmaZBra <- function(n) .qfunc("SigmaZBra", list(n))

#' @rdname pauli
#' @export
qsimplify_pauli <- function(e) {
  if (!is_sym(e)) return(e)
  if (.qisa(e, "Operator")) return(e)
  if (e$type %in% c("Add", "Pow") || (e$type == "Function" && identical(e$name, "exp")))
    return(sym_rebuild(e, lapply(e$args, qsimplify_pauli)))
  if (e$type == "Mul") {
    cnc <- .args_cnc(e)
    c <- cnc[[1L]]; nc <- cnc[[2L]]
    nc_s <- list()
    while (length(nc)) {
      curr <- nc[[1L]]; nc <- nc[-1L]
      while (length(nc) && .qisa(curr, "SigmaOpBase") && .qisa(nc[[1L]], "SigmaOpBase") &&
             .qeq(.qa(curr, "name"), .qa(nc[[1L]], "name"))) {
        x <- nc[[1L]]; nc <- nc[-1L]
        y <- .qsimplify_pauli_product(curr, x)
        cnc1 <- .args_cnc(y)
        curr <- .mul(cnc1[[2L]])
        c <- c(c, cnc1[[1L]])
      }
      nc_s[[length(nc_s) + 1L]] <- curr
    }
    return(.mul(list(.mul(c), .mul(nc_s))))
  }
  e
}

# ---------------------------------------------------------------------------
# One-dimensional harmonic oscillator (sympy/physics/quantum/sho1d.py)

.qsho <- function(nm) switch(nm,
  ad = RaisingOp("a"), a = LoweringOp("a"), H = Hamiltonian("H"), N = NumberOp("N"),
  omega = Symbol("omega"), m = Symbol("m"), X = XOp("X"), Px = PxOp("Px"))

.qsho_matrix <- function(ndim, format, f, scale = NULL) {
  ndim <- .qint(ndim)
  if (!identical(format, "sympy")) .qstop("NotImplementedError", "The ", format, " format is not available in R.")
  vals <- rep(list(S_Zero), ndim * ndim)
  for (ij in f(ndim)) vals[[(ij[[1L]] - 1L) * ndim + ij[[2L]]]] <- ij[[3L]]
  .new_mat(vals, ndim, ndim)
}

.qclass("SHOOp", "Operator",
  `_eval_args` = .cm(function(cls, args, ...) {
    args <- .qsuper("QExpr", "_eval_args", cls, args)
    if (length(args) == 1L) return(args)
    .qstop("ValueError", "Too many arguments")
  }),
  `_eval_hilbert_space` = .cm(function(cls, label) ComplexSpace(oo)),
  `_represent_default_basis` = function(self, ...) .qmeth(self, "_represent_NumberOp", NULL, ...),
  `_represent_XOp` = function(self, basis, ...) .qni("Position representation is not implemented"))

.qinv_sqrt <- function(x) .pow(x, .neg(S_Half))

.qclass("RaisingOp", "SHOOp",
  `_eval_rewrite_as_xp` = function(self, args, ...) {
    s <- .qsho
    .mul(list(.qinv_sqrt(.mul(list(S_Two, hbar, s("m"), s("omega")))),
              .add(list(.mul(list(S_NegativeOne, sym_I, s("Px"))), .mul(list(s("m"), s("omega"), s("X")))))))
  },
  `_eval_adjoint` = function(self) .qfunc("LoweringOp", self$args),
  `_eval_commutator_LoweringOp` = function(self, other, ...) S_NegativeOne,
  `_eval_commutator_NumberOp` = function(self, other, ...) .mul(list(S_NegativeOne, self)),
  `_apply_operator_SHOKet` = function(self, ket, ...) {
    temp <- .add(list(.qa(ket, "n"), S_One))
    .mul(list(.sqrt(temp), SHOKet(temp)))
  },
  `_represent_NumberOp` = function(self, basis, ndim = 4, format = "sympy", ...)
    .qsho_matrix(ndim, format, function(n) lapply(seq_len(n - 1L), function(i) list(i + 1L, i, .sqrt(as_sym(i))))),
  `_print_contents` = function(self, p) sprintf("%s(%s)", .qpyname(self$type), .qpr(p, self$args[[1L]])),
  `_print_contents_latex` = function(self, p) sprintf("%s^{\\dagger}", .latex(self$args[[1L]])))

.qclass("LoweringOp", "SHOOp",
  `_eval_rewrite_as_xp` = function(self, args, ...) {
    s <- .qsho
    .mul(list(.qinv_sqrt(.mul(list(S_Two, hbar, s("m"), s("omega")))),
              .add(list(.mul(list(sym_I, s("Px"))), .mul(list(s("m"), s("omega"), s("X")))))))
  },
  `_eval_adjoint` = function(self) .qfunc("RaisingOp", self$args),
  `_eval_commutator_RaisingOp` = function(self, other, ...) S_One,
  `_eval_commutator_NumberOp` = function(self, other, ...) self,
  `_apply_operator_SHOKet` = function(self, ket, ...) {
    n <- .qa(ket, "n")
    if (.qeq(n, S_Zero)) return(S_Zero)
    .mul(list(.sqrt(n), SHOKet(.sub(n, S_One))))
  },
  `_represent_NumberOp` = function(self, basis, ndim = 4, format = "sympy", ...)
    .qsho_matrix(ndim, format, function(n) lapply(seq_len(n - 1L), function(i) list(i, i + 1L, .sqrt(as_sym(i))))))

.qclass("NumberOp", "SHOOp",
  `_eval_rewrite_as_a` = function(self, args, ...) .mul(list(.qsho("ad"), .qsho("a"))),
  `_eval_rewrite_as_xp` = function(self, args, ...) {
    s <- .qsho
    .sub(.mul(list(.pow(.mul(list(S_Two, s("m"), hbar, s("omega"))), S_NegativeOne),
                   .add(list(.pow(s("Px"), S_Two), .pow(.mul(list(s("m"), s("omega"), s("X"))), S_Two))))), S_Half)
  },
  `_eval_rewrite_as_H` = function(self, args, ...)
    .sub(.mul(list(.qsho("H"), .pow(.mul(list(hbar, .qsho("omega"))), S_NegativeOne))), S_Half),
  `_apply_operator_SHOKet` = function(self, ket, ...) .mul(list(.qa(ket, "n"), ket)),
  `_eval_commutator_Hamiltonian` = function(self, other, ...) S_Zero,
  `_eval_commutator_RaisingOp` = function(self, other, ...) other,
  `_eval_commutator_LoweringOp` = function(self, other, ...) .mul(list(S_NegativeOne, other)),
  `_represent_NumberOp` = function(self, basis, ndim = 4, format = "sympy", ...)
    .qsho_matrix(ndim, format, function(n) lapply(seq_len(n), function(i) list(i, i, as_sym(i - 1L)))))

.qclass("Hamiltonian", "SHOOp",
  `_eval_rewrite_as_a` = function(self, args, ...)
    .mul(list(hbar, .qsho("omega"), .add(list(.mul(list(.qsho("ad"), .qsho("a"))), S_Half)))),
  `_eval_rewrite_as_xp` = function(self, args, ...) {
    s <- .qsho
    .mul(list(.pow(.mul(list(S_Two, s("m"))), S_NegativeOne),
              .add(list(.pow(s("Px"), S_Two), .pow(.mul(list(s("m"), s("omega"), s("X"))), S_Two)))))
  },
  `_eval_rewrite_as_N` = function(self, args, ...) .mul(list(hbar, .qsho("omega"), .add(list(.qsho("N"), S_Half)))),
  `_apply_operator_SHOKet` = function(self, ket, ...) .mul(list(hbar, .qsho("omega"), .add(list(.qa(ket, "n"), S_Half)), ket)),
  `_eval_commutator_NumberOp` = function(self, other, ...) S_Zero,
  `_represent_NumberOp` = function(self, basis, ndim = 4, format = "sympy", ...) {
    M <- .qsho_matrix(ndim, format, function(n) lapply(seq_len(n), function(i) list(i, i, .add(list(as_sym(i - 1L), S_Half)))))
    .mat_map(M, function(z) .mul(list(hbar, .qsho("omega"), z)))
  })

.qclass("SHOState", "State",
  `_eval_hilbert_space` = .cm(function(cls, label) ComplexSpace(oo)),
  n = .prop(function(self) self$args[[1L]]))

.qsho_vector <- function(self, ndim, format, column) {
  ndim <- .qint(ndim)
  n <- .qa(self, "n")
  if (n$type != "Integer") return(.qerror("ValueError", "Not Numerical State"))
  if (.qint(n) >= ndim) return(.qerror("ValueError", "N-Dimension too small"))
  if (!identical(format, "sympy")) .qstop("NotImplementedError", "The ", format, " format is not available in R.")
  vals <- rep(list(S_Zero), ndim)
  vals[[.qint(n) + 1L]] <- S_One
  if (column) .new_mat(vals, ndim, 1L) else .new_mat(vals, 1L, ndim)
}

.qclass("SHOKet", c("SHOState", "Ket"),
  dual_class = .cm(function(cls) "SHOBra"),
  `_eval_innerproduct_SHOBra` = function(self, bra, ...) KroneckerDelta(.qa(self, "n"), .qa(bra, "n")),
  `_represent_default_basis` = function(self, ...) .qmeth(self, "_represent_NumberOp", NULL, ...),
  `_represent_NumberOp` = function(self, basis, ndim = 4, format = "sympy", ...) .qsho_vector(self, ndim, format, TRUE))

.qclass("SHOBra", c("SHOState", "Bra"),
  dual_class = .cm(function(cls) "SHOKet"),
  `_represent_default_basis` = function(self, ...) .qmeth(self, "_represent_NumberOp", NULL, ...),
  `_represent_NumberOp` = function(self, basis, ndim = 4, format = "sympy", ...) .qsho_vector(self, ndim, format, FALSE))

#' One-dimensional harmonic oscillator
#'
#' The ladder operators `RaisingOp()` and `LoweringOp()`, the number
#' operator `NumberOp()`, the Hamiltonian `Hamiltonian()` and the states
#' `SHOKet(n)` and `SHOBra(n)` of `sympy.physics.quantum.sho1d`. Each
#' operator takes one label. [represent()] gives matrices of size `ndim`
#' (4 by default) in the number basis, and [rewrite()] writes the
#' operators in terms of `"a"`, `"xp"`, `"N"` or `"H"`.
#'
#' @param ... The label.
#' @param n The number of quanta.
#' @return A quantum operator or state.
#' @examples
#' ad <- RaisingOp("a"); a <- LoweringOp("a")
#' qapply(ad * SHOKet(1))
#' represent(NumberOp("N"), ndim = 3)
#' doit(Commutator(a, ad))
#' @name sho1d
NULL

#' @rdname sho1d
#' @export
RaisingOp <- function(...) .qfunc("RaisingOp", list(...))
#' @rdname sho1d
#' @export
LoweringOp <- function(...) .qfunc("LoweringOp", list(...))
#' @rdname sho1d
#' @export
NumberOp <- function(...) .qfunc("NumberOp", list(...))
#' @rdname sho1d
#' @export
Hamiltonian <- function(...) .qfunc("Hamiltonian", list(...))
#' @rdname sho1d
#' @export
SHOKet <- function(n) .qfunc("SHOKet", list(n))
#' @rdname sho1d
#' @export
SHOBra <- function(n) .qfunc("SHOBra", list(n))
