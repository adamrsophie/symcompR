# Quantum gates, after sympy/physics/quantum/gate.py, matrixcache.py and
# matrixutils.py. Only the "sympy" matrix format exists in R; the NumPy and
# SciPy formats of matrixutils.py have no R counterpart.

# ---------------------------------------------------------------------------
# matrixutils.py and matrixcache.py

.qmatrix_cache <- local({
  cache <- NULL
  function(name, format = "sympy") {
    if (is.null(cache)) {
      sqrt2_inv <- .pow(S_Two, S_NegativeOne * S_Half)
      sqrt2_inv <- .pow(.sqrt(S_Two), S_NegativeOne)
      M <- function(...) Matrix(list(...))
      cache <<- list(
        eye2 = M(list(1, 0), list(0, 1)), op11 = M(list(0, 0), list(0, 1)), op00 = M(list(1, 0), list(0, 0)),
        op10 = M(list(0, 0), list(1, 0)), op01 = M(list(0, 1), list(0, 0)), X = M(list(0, 1), list(1, 0)),
        Y = M(list(0, .neg(sym_I)), list(sym_I, 0)), Z = M(list(1, 0), list(0, -1)), S = M(list(1, 0), list(0, sym_I)),
        T = M(list(1, 0), list(0, exp(.mul(list(sym_I, sym_pi, .pow(S(4), S_NegativeOne)))))),
        H = .mat_map(M(list(1, 1), list(1, -1)), function(z) .mul(list(sqrt2_inv, z))),
        Hsqrt2 = M(list(1, 1), list(1, -1)),
        SWAP = M(list(1, 0, 0, 0), list(0, 0, 1, 0), list(0, 1, 0, 0), list(0, 0, 0, 1)),
        ZX = .mat_map(M(list(1, 1), list(1, -1)), function(z) .mul(list(sqrt2_inv, z))),
        ZY = M(list(sym_I, 0), list(0, .neg(sym_I))))
    }
    m <- if (identical(format, "sympy")) cache[[name]] else NULL
    if (is.null(m)) .qni("Matrix with name ", name, " and format ", format, " is not available.")
    m
  }
})

#' Matrix utilities of the quantum module
#'
#' The SymPy matrix helpers of `sympy.physics.quantum.matrixutils`:
#' `matrix_tensor_product()` is the Kronecker product of its arguments,
#' `matrix_eye(n)` and `matrix_zeros(m, n)` the identity and zero matrices,
#' `matrix_dagger(m)` the conjugate transpose, `to_sympy(m)` (one argument)
#' the identity on symbolic matrices and `flatten_scalar(m)` the entry of a 1 by 1 matrix.
#' Only the `"sympy"` format exists here; NumPy and SciPy formats raise an
#' error.
#'
#' @param ... Matrices.
#' @param n,m Dimensions.
#' @param format The matrix format, only `"sympy"`.
#' @param e A matrix or an expression.
#' @return A matrix or an expression.
#' @examples
#' matrix_tensor_product(Matrix(list(list(1, 0), list(0, 1))), Matrix(list(list(0, 1), list(1, 0))))
#' @name quantum-matrixutils
NULL

#' @rdname quantum-matrixutils
#' @export
matrix_tensor_product <- function(...) {
  ps <- list(...)
  if (inherits(ps[[1L]], "symmatrix")) return(.matrix_tensor_product(ps))
  NULL
}

.qformat_check <- function(format) {
  if (!identical(format, "sympy")) .qstop("NotImplementedError", "The ", format, " format is not available in R.")
}

#' @rdname quantum-matrixutils
#' @export
matrix_eye <- function(n, format = "sympy") { .qformat_check(format); eye(.qint(n)) }

#' @rdname quantum-matrixutils
#' @export
matrix_zeros <- function(m, n, format = "sympy") { .qformat_check(format); zeros(.qint(m), .qint(n)) }

#' @rdname quantum-matrixutils
#' @export
matrix_dagger <- function(e) {
  if (inherits(e, "symmatrix")) return(Dagger(e))
  .qstop("TypeError", "Expected sympy/numpy/scipy.sparse matrix, got: ", .qrepr(e))
}

#' @rdname quantum-matrixutils
#' @export
matrix_to_zero <- function(e) {
  if (inherits(e, "symmatrix") &&
      all(vapply(e$data, function(z) .qeq(z, S_Zero) || (is_number_type(z) && num_is_zero(z)), NA))) return(S_Zero)
  e
}

# matrixutils.to_sympy(m); the exported to_sympy() calls it with one argument
.qto_sympy <- function(e) {
  if (inherits(e, "symmatrix") || is_sym(e)) return(e)
  .qstop("TypeError", "Expected sympy/numpy/scipy.sparse matrix, got: ", .qrepr(e))
}

# ---------------------------------------------------------------------------
# helpers

.qsort_by_key <- function(xs) .key_sort(lapply(xs, .sort_key), xs)
.qmax <- function(xs) { s <- .qsort_by_key(lapply(xs, as_sym)); s[[length(s)]] }
.qmin <- function(xs) .qsort_by_key(lapply(xs, as_sym))[[1L]]
.qseq_of <- function(x) if (is_sym(x) && x$type == "Tuple") x$args else if (is.list(x) && !is_sym(x)) x else list(x)
.qis_seq <- function(x) (is_sym(x) && x$type == "Tuple") || (is.list(x) && !is_sym(x) && !inherits(x, "symmatrix"))

.qvalidate_targets_controls <- function(tandc) {
  for (bit in tandc)
    if (!(bit$type == "Integer" || bit$type %in% c("Symbol", "Dummy", "Wild")))
      .qstop("TypeError", "Integer expected, got: ", .qrepr(bit))
  keys <- vapply(tandc, .k, "")
  if (anyDuplicated(keys)) .qstop("QuantumError", "Target/control qubits in a gate cannot be duplicated")
}

.qgate_hs <- function(n) TensorPowerHilbertSpace(ComplexSpace(2), n)
.qnormalized <- new.env()
.qnormalized$on <- TRUE

#' Use normalized Hadamard gates
#'
#' `normalized(FALSE)` makes the Hadamard gate's matrix `[[1, 1], [1, -1]]`
#' instead of the normalized one, as SymPy's `gate.normalized`.
#'
#' @param normalize `TRUE` or `FALSE`.
#' @return `NULL`, invisibly.
#' @examples
#' normalized(TRUE)
#' @export
normalized <- function(normalize) { .qnormalized$on <- isTRUE(normalize); invisible(NULL) }

# ---------------------------------------------------------------------------
# Gates

.qclass("Gate", "UnitaryOperator",
  `_label_separator` = ",",
  gate_name = "G", gate_name_latex = "G",
  `_eval_args` = .cm(function(cls, args, ...) {
    args <- .qsuper("UnitaryOperator", "_eval_args", cls, args)
    .qvalidate_targets_controls(args)
    args
  }),
  `_eval_hilbert_space` = .cm(function(cls, args) .qgate_hs(.add(list(.qmax(args), S_One)))),
  nqubits = .prop(function(self) length(.qa(self, "targets"))),
  min_qubits = .prop(function(self) .add(list(.qmax(.qa(self, "targets")), S_One))),
  targets = .prop(function(self) .qa(self, "label")),
  gate_name_plot = .prop(function(self) paste0("$", .qa(self, "gate_name_latex"), "$")),
  get_target_matrix = function(self, format = "sympy") .qni("get_target_matrix is not implemented in Gate."),
  `_apply_operator_IntQubit` = function(self, qubits, ...) .qmeth(self, "_apply_operator_Qubit", qubits, ...),
  `_apply_operator_Qubit` = function(self, qubits, ...) {
    nq <- .qa(qubits, "nqubits"); mq <- .qa(self, "min_qubits")
    if (nq < .qint(mq))
      .qstop("QuantumError", "Gate needs a minimum of ", .qrepr(mq), " qubits to act on, got: ", nq)
    if (.qisa(self, "CGate") && !.qmeth(self, "eval_controls", qubits)) return(qubits)
    targets <- .qa(self, "targets")
    target_matrix <- .qmeth(self, "get_target_matrix", format = "sympy")
    column_index <- 0; n <- 1
    for (target in targets) {
      column_index <- column_index + n * .qint(.qmeth(qubits, "__getitem__", target))
      n <- n * 2
    }
    result <- S_Zero
    for (index in seq_len(target_matrix$nrow) - 1L) {
      new_qubit <- .qfunc(qubits$type, qubits$args)
      for (bit in seq_along(targets) - 1L) {
        target <- targets[[bit + 1L]]
        if (.qint(.qmeth(new_qubit, "__getitem__", target)) != bitwAnd(bitwShiftR(index, bit), 1L))
          new_qubit <- .qmeth(new_qubit, "flip", target)
      }
      result <- .add(list(result, .mul(list(.mget(target_matrix, index + 1L, column_index + 1), new_qubit))))
    }
    result
  },
  `_represent_default_basis` = function(self, ...) .qmeth(self, "_represent_ZGate", NULL, ...),
  `_represent_ZGate` = function(self, basis, format = "sympy", nqubits = 0, ...) {
    nqubits <- .qint(nqubits)
    if (nqubits == 0) .qstop("QuantumError", "The number of qubits must be given as nqubits.")
    if (nqubits < .qint(.qa(self, "min_qubits")))
      .qstop("QuantumError", "The number of qubits ", nqubits, " is too small for the gate.")
    target_matrix <- .qmeth(self, "get_target_matrix", format)
    targets <- .qa(self, "targets")
    controls <- if (.qisa(self, "CGate")) .qa(self, "controls") else list()
    represent_zbasis(controls, targets, target_matrix, nqubits, format)
  },
  `_sympystr` = function(self, p) paste0(.qa(self, "gate_name"), "(", .qmeth(self, "_print_label", p), ")"),
  `_latex` = function(self, p) paste0(.qa(self, "gate_name_latex"), "_{", .qmeth(self, "_print_label", p), "}"),
  plot_gate = function(self, ...) .qni("plot_gate is not implemented."))

.qclass("CGate", "Gate",
  gate_name = "C", gate_name_latex = "C",
  control_value = S_One, simplify_cgate = FALSE,
  `_eval_args` = .cm(function(cls, args, ...) {
    controls <- args[[1L]]; gate <- args[[2L]]
    controls <- .qseq_of(controls)
    controls <- .qsuper("UnitaryOperator", "_eval_args", cls, controls)
    .qvalidate_targets_controls(c(controls, .qa(gate, "targets")))
    list(do.call(Tuple, controls), gate)
  }),
  `_eval_hilbert_space` = .cm(function(cls, args)
    .qgate_hs(.qmax(list(.add(list(.qmax(args[[1L]]$args), S_One)), .qa(args[[2L]], "min_qubits"))))),
  nqubits = .prop(function(self) length(.qa(self, "targets")) + length(.qa(self, "controls"))),
  min_qubits = .prop(function(self) .add(list(.qmax(list(.qmax(.qa(self, "controls")), .qmax(.qa(self, "targets")))), S_One))),
  targets = .prop(function(self) .qa(.qa(self, "gate"), "targets")),
  controls = .prop(function(self) .qa(self, "label")[[1L]]$args),
  gate = .prop(function(self) .qa(self, "label")[[2L]]),
  get_target_matrix = function(self, format = "sympy") .qmeth(.qa(self, "gate"), "get_target_matrix", format),
  eval_controls = function(self, qubit) {
    cv <- .qa(self, "control_value")
    all(vapply(.qa(self, "controls"), function(bit) .qeq(.qmeth(qubit, "__getitem__", bit), cv), NA))
  },
  decompose = function(self, ...) {
    ctl <- .qa(self, "controls")
    if (length(ctl) == 1L) {
      c <- ctl[[1L]]
      g <- .qa(self, "gate")
      t <- .qa(g, "targets")[[1L]]
      if (.qisa(g, "YGate")) return(.mul(list(PhaseGate(t), CNotGate(c, t), PhaseGate(t), ZGate(t))))
      if (.qisa(g, "ZGate")) return(.mul(list(HadamardGate(t), CNotGate(c, t), HadamardGate(t))))
      return(NULL)
    }
    self
  },
  `_print_label` = function(self, p) {
    controls <- .qmeth(self, "_print_sequence", .qa(self, "controls"), ",", p)
    paste0("(", controls, "),", .qpr(p, .qa(self, "gate")))
  },
  `_latex` = function(self, p) {
    controls <- .qmeth(self, "_print_sequence", .qa(self, "controls"), ",", p)
    paste0(.qa(self, "gate_name_latex"), "_{", controls, "}{\\left(", .qpr(p, .qa(self, "gate")), "\\right)}")
  },
  `_eval_inverse` = function(self) {
    if (.qisa(.qa(self, "gate"), "HermitianOperator")) self else .qsuper("Gate", "_eval_inverse", self)
  },
  `_eval_power` = function(self, ex) {
    if (.qisa(.qa(self, "gate"), "HermitianOperator")) {
      if (.qeq(ex, S_NegativeOne)) return(NULL)
      if (ex$type == "Integer" && .qint(ex) %% 2 == 0) return(.mul(list(self, .qsuper("Gate", "_eval_inverse", self))))
      return(self)
    }
    NULL
  })

.qclass("CGateS", "CGate", simplify_cgate = TRUE)

.qclass("UGate", "Gate",
  gate_name = "U", gate_name_latex = "U",
  `_eval_args` = .cm(function(cls, args, ...) {
    targets <- .qseq_of(args[[1L]])
    targets <- .qsuper("Gate", "_eval_args", cls, targets)
    .qvalidate_targets_controls(targets)
    mat <- args[[2L]]
    if (is_sym(mat) && mat$type == "ImmutableMatrix") mat <- as_explicit(mat)
    if (!inherits(mat, "symmatrix")) .qstop("TypeError", "Matrix expected, got: ", .qrepr(mat))
    dim <- 2^length(targets)
    if (!(mat$nrow == dim && mat$ncol == dim))
      .qstop("IndexError", "Number of targets must match the matrix size: ", .qrepr(do.call(Tuple, targets)), " ", .qmat_str(mat))
    list(do.call(Tuple, targets), .mat_imm(mat))
  }),
  `_eval_hilbert_space` = .cm(function(cls, args) .qgate_hs(.add(list(.qmax(args[[1L]]$args), S_One)))),
  targets = .prop(function(self) .qa(self, "label")[[1L]]$args),
  get_target_matrix = function(self, format = "sympy") as_explicit(.qa(self, "label")[[2L]]),
  `_latex` = function(self, p) paste0(.qa(self, "gate_name_latex"), "_{", .qmeth(self, "_print_sequence", .qa(self, "targets"), ",", p), "}"))

.qclass("OneQubitGate", "Gate",
  nqubits = S_One,
  `_eval_commutator` = function(self, other, ...) {
    if (.qisa(other, "OneQubitGate")) {
      if (!identical(vapply(.qa(self, "targets"), .k, ""), vapply(.qa(other, "targets"), .k, "")) || self$type == other$type)
        return(S_Zero)
    }
    .qsuper("Operator", "_eval_commutator", self, other, ...)
  },
  `_eval_anticommutator` = function(self, other, ...) {
    if (.qisa(other, "OneQubitGate")) {
      if (!identical(vapply(.qa(self, "targets"), .k, ""), vapply(.qa(other, "targets"), .k, "")) || self$type == other$type)
        return(.mul(list(S_Two, self, other)))
    }
    .qsuper("Operator", "_eval_anticommutator", self, other, ...)
  })

.qclass("TwoQubitGate", "Gate", nqubits = S_Two)

.qt0 <- function(self) .qa(self, "targets")[[1L]]

.qclass("IdentityGate", "OneQubitGate",
  is_hermitian = TRUE, gate_name = "1", gate_name_latex = "1",
  `_apply_operator_Qubit` = function(self, qubits, ...) {
    nq <- .qa(qubits, "nqubits"); mq <- .qa(self, "min_qubits")
    if (nq < .qint(mq)) .qstop("QuantumError", "Gate needs a minimum of ", .qrepr(mq), " qubits to act on, got: ", nq)
    qubits
  },
  get_target_matrix = function(self, format = "sympy") .qmatrix_cache("eye2", format),
  `_eval_commutator` = function(self, other, ...) S_Zero,
  `_eval_anticommutator` = function(self, other, ...) .mul(list(S_Two, other)))

.qsqrt2 <- function() .sqrt(S_Two)

.qclass("HadamardGate", c("HermitianOperator", "OneQubitGate"),
  gate_name = "H", gate_name_latex = "H",
  get_target_matrix = function(self, format = "sympy")
    .qmatrix_cache(if (.qnormalized$on) "H" else "Hsqrt2", format),
  `_eval_commutator_XGate` = function(self, other, ...) .mul(list(sym_I, .qsqrt2(), YGate(.qt0(self)))),
  `_eval_commutator_YGate` = function(self, other, ...)
    .mul(list(sym_I, .qsqrt2(), .sub(ZGate(.qt0(self)), XGate(.qt0(self))))),
  `_eval_commutator_ZGate` = function(self, other, ...) .mul(list(S_NegativeOne, sym_I, .qsqrt2(), YGate(.qt0(self)))),
  `_eval_anticommutator_XGate` = function(self, other, ...) .mul(list(.qsqrt2(), IdentityGate(.qt0(self)))),
  `_eval_anticommutator_YGate` = function(self, other, ...) S_Zero,
  `_eval_anticommutator_ZGate` = function(self, other, ...) .mul(list(.qsqrt2(), IdentityGate(.qt0(self)))))

.qclass("XGate", c("HermitianOperator", "OneQubitGate"),
  gate_name = "X", gate_name_latex = "X",
  get_target_matrix = function(self, format = "sympy") .qmatrix_cache("X", format),
  `_eval_commutator_YGate` = function(self, other, ...) .mul(list(S_Two, sym_I, ZGate(.qt0(self)))),
  `_eval_anticommutator_XGate` = function(self, other, ...) .mul(list(S_Two, IdentityGate(.qt0(self)))),
  `_eval_anticommutator_YGate` = function(self, other, ...) S_Zero,
  `_eval_anticommutator_ZGate` = function(self, other, ...) S_Zero)

.qclass("YGate", c("HermitianOperator", "OneQubitGate"),
  gate_name = "Y", gate_name_latex = "Y",
  get_target_matrix = function(self, format = "sympy") .qmatrix_cache("Y", format),
  `_eval_commutator_ZGate` = function(self, other, ...) .mul(list(S_Two, sym_I, XGate(.qt0(self)))),
  `_eval_anticommutator_YGate` = function(self, other, ...) .mul(list(S_Two, IdentityGate(.qt0(self)))),
  `_eval_anticommutator_ZGate` = function(self, other, ...) S_Zero)

.qclass("ZGate", c("HermitianOperator", "OneQubitGate"),
  gate_name = "Z", gate_name_latex = "Z",
  get_target_matrix = function(self, format = "sympy") .qmatrix_cache("Z", format),
  `_eval_commutator_XGate` = function(self, other, ...) .mul(list(S_Two, sym_I, YGate(.qt0(self)))),
  `_eval_anticommutator_YGate` = function(self, other, ...) S_Zero)

.qclass("PhaseGate", "OneQubitGate",
  is_hermitian = FALSE, gate_name = "S", gate_name_latex = "S",
  get_target_matrix = function(self, format = "sympy") .qmatrix_cache("S", format),
  `_eval_commutator_ZGate` = function(self, other, ...) S_Zero,
  `_eval_commutator_TGate` = function(self, other, ...) S_Zero)

.qclass("TGate", "OneQubitGate",
  is_hermitian = FALSE, gate_name = "T", gate_name_latex = "T",
  get_target_matrix = function(self, format = "sympy") .qmatrix_cache("T", format),
  `_eval_commutator_ZGate` = function(self, other, ...) S_Zero,
  `_eval_commutator_PhaseGate` = function(self, other, ...) S_Zero)

.qcnot_comm <- function(self, a, b, other) {
  if (.qeq(a, b)) return(S_Zero)
  .qni("Commutator not implemented: ", .qrepr(other))
}

.qclass("CNotGate", c("HermitianOperator", "CGate", "TwoQubitGate"),
  gate_name = "CNOT", gate_name_latex = "\\text{CNOT}", simplify_cgate = TRUE,
  `_eval_args` = .cm(function(cls, args, ...) .qsuper("Gate", "_eval_args", cls, args)),
  `_eval_hilbert_space` = .cm(function(cls, args) .qgate_hs(.add(list(.qmax(args), S_One)))),
  min_qubits = .prop(function(self) .add(list(.qmax(.qa(self, "label")), S_One))),
  targets = .prop(function(self) list(.qa(self, "label")[[2L]])),
  controls = .prop(function(self) list(.qa(self, "label")[[1L]])),
  gate = .prop(function(self) XGate(.qa(self, "label")[[2L]])),
  `_print_label` = function(self, p) .qsuper("Gate", "_print_label", self, p),
  `_latex` = function(self, p) .qsuper("Gate", "_latex", self, p),
  `_eval_commutator_ZGate` = function(self, other, ...) .qcnot_comm(self, .qa(self, "controls")[[1L]], .qt0(other), other),
  `_eval_commutator_TGate` = function(self, other, ...) .qmeth(self, "_eval_commutator_ZGate", other, ...),
  `_eval_commutator_PhaseGate` = function(self, other, ...) .qmeth(self, "_eval_commutator_ZGate", other, ...),
  `_eval_commutator_XGate` = function(self, other, ...) .qcnot_comm(self, .qt0(self), .qt0(other), other),
  `_eval_commutator_CNotGate` = function(self, other, ...)
    .qcnot_comm(self, .qa(self, "controls")[[1L]], .qa(other, "controls")[[1L]], other))

.qclass("SwapGate", "TwoQubitGate",
  is_hermitian = TRUE, gate_name = "SWAP", gate_name_latex = "\\text{SWAP}",
  get_target_matrix = function(self, format = "sympy") .qmatrix_cache("SWAP", format),
  decompose = function(self, ...) {
    tg <- .qa(self, "targets")
    g1 <- CNotGate(tg[[1L]], tg[[2L]]); g2 <- CNotGate(tg[[2L]], tg[[1L]])
    .mul(list(g1, g2, g1))
  },
  `_represent_ZGate` = function(self, basis, format = "sympy", nqubits = NULL, ...) {
    targets <- vapply(.qa(self, "targets"), .qint, 0)
    min_target <- min(targets); max_target <- max(targets)
    nqubits <- if (is.null(nqubits)) .qint(.qa(self, "min_qubits")) else .qint(nqubits)
    mc <- function(n) .qmatrix_cache(n, format)
    result <- NULL
    for (pr in list(c("op01", "op10"), c("op10", "op01"), c("op00", "op00"), c("op11", "op11"))) {
      product <- rep(list(mc("eye2")), nqubits)
      product[[nqubits - min_target]] <- mc(pr[1L])
      product[[nqubits - max_target]] <- mc(pr[2L])
      new_result <- .matrix_tensor_product(product)
      result <- if (is.null(result)) new_result else result + new_result
    }
    result
  })

#' Quantum gates
#'
#' The gates of `sympy.physics.quantum.gate` acting on qubits: the one qubit
#' gates `IdentityGate()`, `HadamardGate()`, `XGate()`, `YGate()`, `ZGate()`,
#' `PhaseGate()` (SymPy's `S`) and `TGate()` take a target qubit;
#' `CNotGate(control, target)` (also `CNOT()`), `SwapGate(i, j)` (also
#' `SWAP()`), `CPHASE(a, b)`, the controlled gate `CGate(controls, gate)`,
#' `CGateS()` and the general `UGate(targets, matrix)` and `Gate()`.
#' Products of gates applied to [Qubit()] states are evaluated by [qapply()];
#' [represent()] with `nqubits` gives the matrix. `gate_sort()` and
#' `gate_simp()` reorder and simplify circuits and `random_circuit()` builds
#' a random one. `represent_zbasis()` is the matrix of a controlled gate.
#'
#' @param ... Target qubits, or the arguments of the gate.
#' @param controls,targets Control and target qubits.
#' @param gate A gate.
#' @param mat A matrix.
#' @param a,b Qubits.
#' @param circuit A product of gates.
#' @param ngates,nqubits The number of gates and qubits.
#' @param gate_space The gate classes to draw from.
#' @param target_matrix The matrix of the gate on its targets.
#' @param format The matrix format, only `"sympy"`.
#' @return A gate, a circuit or a matrix.
#' @examples
#' qapply(HadamardGate(0) * Qubit("0"))
#' represent(CNOT(1, 0), nqubits = 2)
#' gate_simp(HadamardGate(0) * HadamardGate(0))
#' @name quantum-gates
NULL

#' @rdname quantum-gates
#' @export
Gate <- function(...) .qfunc("Gate", list(...))
#' @rdname quantum-gates
#' @export
CGate <- function(controls, gate) .qfunc("CGate", list(controls, gate))
#' @rdname quantum-gates
#' @export
CGateS <- function(controls, gate) .qfunc("CGateS", list(controls, gate))
#' @rdname quantum-gates
#' @export
UGate <- function(targets, mat) .qfunc("UGate", list(targets, mat))
#' @rdname quantum-gates
#' @export
OneQubitGate <- function(...) .qfunc("OneQubitGate", list(...))
#' @rdname quantum-gates
#' @export
TwoQubitGate <- function(...) .qfunc("TwoQubitGate", list(...))
#' @rdname quantum-gates
#' @export
IdentityGate <- function(...) .qfunc("IdentityGate", list(...))
#' @rdname quantum-gates
#' @export
HadamardGate <- function(...) .qfunc("HadamardGate", list(...))
#' @rdname quantum-gates
#' @export
XGate <- function(...) .qfunc("XGate", list(...))
#' @rdname quantum-gates
#' @export
YGate <- function(...) .qfunc("YGate", list(...))
#' @rdname quantum-gates
#' @export
ZGate <- function(...) .qfunc("ZGate", list(...))
#' @rdname quantum-gates
#' @export
PhaseGate <- function(...) .qfunc("PhaseGate", list(...))
#' @rdname quantum-gates
#' @export
TGate <- function(...) .qfunc("TGate", list(...))
#' @rdname quantum-gates
#' @export
CNotGate <- function(...) .qfunc("CNotGate", list(...))
#' @rdname quantum-gates
#' @export
CNOT <- CNotGate
#' @rdname quantum-gates
#' @export
SwapGate <- function(...) .qfunc("SwapGate", list(...))
#' @rdname quantum-gates
#' @export
SWAP <- SwapGate
#' @rdname quantum-gates
#' @export
CPHASE <- function(a, b) CGateS(list(a), ZGate(b))

#' @rdname quantum-gates
#' @export
represent_zbasis <- function(controls, targets, target_matrix, nqubits, format = "sympy") {
  controls <- vapply(controls, .qint, 0)
  targets <- vapply(targets, .qint, 0)
  nqubits <- .qint(nqubits)
  op11 <- .qmatrix_cache("op11", format)
  eye2 <- .qmatrix_cache("eye2", format)
  if (!length(controls) && length(targets) == 1L) {
    product <- list()
    bit <- targets[1L]
    if (bit != nqubits - 1) product[[length(product) + 1L]] <- matrix_eye(2^(nqubits - bit - 1), format = format)
    product[[length(product) + 1L]] <- target_matrix
    if (bit != 0) product[[length(product) + 1L]] <- matrix_eye(2^bit, format = format)
    return(.matrix_tensor_product(product))
  }
  if (length(targets) == 1L && length(controls) >= 1L) {
    target <- targets[1L]
    product2 <- rep(list(matrix_eye(2, format = format)), nqubits)
    for (control in controls) product2[[nqubits - control]] <- op11
    product2[[nqubits - target]] <- target_matrix - eye2
    return(matrix_eye(2^nqubits, format = format) + .matrix_tensor_product(product2))
  }
  .qni("The representation of multi-target, multi-control gates is not implemented.")
}

.qas_base_exp <- function(e) if (e$type == "Pow") list(e$args[[1L]], e$args[[2L]]) else list(e, S_One)

#' @rdname quantum-gates
#' @export
gate_simp <- function(circuit) {
  circuit <- gate_sort(circuit)
  if (circuit$type == "Add") {
    r <- S_Zero
    for (t in .sympy_args(circuit)) r <- .add(list(r, gate_simp(t)))
    return(r)
  } else if (circuit$type == "Mul") {
    circuit_args <- .qargs(circuit)
  } else if (circuit$type == "Pow") {
    circuit_args <- list(.pow(gate_simp(circuit$args[[1L]]), circuit$args[[2L]]))
  } else {
    return(circuit)
  }
  for (i in seq_along(circuit_args)) {
    a <- circuit_args[[i]]
    if (a$type == "Pow") {
      b <- a$args[[1L]]; ex <- a$args[[2L]]
      pre <- circuit_args[seq_len(i - 1L)]
      post <- circuit_args[-seq_len(i)]
      if (.qisa(b, c("HadamardGate", "XGate", "YGate", "ZGate")) && is_number_type(ex)) {
        circuit <- gate_simp(.mul(c(pre, list(.pow(b, .qmod2(ex))), post)))
        break
      } else if (.qisa(b, "PhaseGate")) {
        mid <- list(.pow(ZGate(b$args[[1L]]), .qint_div2(ex)), .pow(b, .qmod2(ex)))
        circuit <- gate_simp(.mul(c(pre, mid, post)))
        break
      } else if (.qisa(b, "TGate")) {
        mid <- list(.pow(PhaseGate(b$args[[1L]]), .qint_div2(ex)), .pow(b, .qmod2(ex)))
        circuit <- gate_simp(.mul(c(pre, mid, post)))
        break
      }
    }
  }
  circuit
}

# e % 2 for a SymPy number (Python's floor modulo)
.qmod2 <- function(ex) .sub(ex, .mul(list(S_Two, floor(.mul(list(ex, S_Half))))))

# Integer(e/2): truncation toward zero
.qint_div2 <- function(ex) as_sym(trunc(as.numeric(N(ex)) / 2))

#' @rdname quantum-gates
#' @export
gate_sort <- function(circuit) {
  if (circuit$type == "Add") {
    r <- S_Zero
    for (t in .sympy_args(circuit)) r <- .add(list(r, gate_sort(t)))
    return(r)
  }
  if (circuit$type == "Pow") return(.pow(gate_sort(circuit$args[[1L]]), circuit$args[[2L]]))
  if (.qisa(circuit, "Gate")) return(circuit)
  if (circuit$type != "Mul") return(circuit)
  changes <- TRUE
  while (changes) {
    changes <- FALSE
    circ_array <- .qargs(circuit)
    n <- length(circ_array)
    for (i in seq_len(n - 1L)) {
      a <- circ_array[[i]]; b <- circ_array[[i + 1L]]
      if (.qisa(a, c("Gate", "Pow")) && .qisa(b, c("Gate", "Pow"))) {
        fb <- .qas_base_exp(a); sb <- .qas_base_exp(b)
        if (.basic_compare(fb[[1L]], sb[[1L]]) > 0) {
          swapped <- c(circ_array[seq_len(i - 1L)], list(b, a), circ_array[-seq_len(i + 1L)])
          if (.qeq(doit(.qCommutator(fb[[1L]], sb[[1L]])), S_Zero)) {
            circuit <- .mul(swapped)
            changes <- TRUE
            break
          }
          if (.qeq(doit(AntiCommutator(fb[[1L]], sb[[1L]])), S_Zero)) {
            sign <- .pow(S_NegativeOne, .mul(list(fb[[2L]], sb[[2L]])))
            circuit <- .mul(list(sign, .mul(swapped)))
            changes <- TRUE
            break
          }
        }
      }
    }
  }
  circuit
}

#' @rdname quantum-gates
#' @export
random_circuit <- function(ngates, nqubits, gate_space = c("XGate", "YGate", "ZGate", "PhaseGate", "TGate", "HadamardGate", "CNotGate", "SwapGate")) {
  nqubits <- .qint(nqubits)
  result <- list()
  for (i in seq_len(.qint(ngates))) {
    g <- gate_space[sample.int(length(gate_space), 1L)]
    if (g %in% c("CNotGate", "SwapGate")) {
      qb <- sample.int(nqubits, 2L) - 1L
      result[[i]] <- .qfunc(g, lapply(qb, as_sym))
    } else {
      result[[i]] <- .qfunc(g, list(as_sym(sample.int(nqubits, 1L) - 1L)))
    }
  }
  .mul(result)
}
