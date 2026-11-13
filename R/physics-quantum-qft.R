# The quantum Fourier transform and Grover's search, after
# sympy/physics/quantum/qft.py and grover.py.

.qclass("RkGate", "OneQubitGate",
  gate_name = "Rk", gate_name_latex = "R",
  `__new__` = .cm(function(cls, args, ...) {
    if (length(args) != 2L) .qstop("QuantumError", "Rk gates only take two arguments, got: ", .qrepr(args))
    target <- args[[1L]]; k <- as_sym(args[[2L]])
    if (.qeq(k, S_One)) return(ZGate(target))
    if (.qeq(k, S_Two)) return(PhaseGate(target))
    if (.qeq(k, as_sym(3L))) return(TGate(target))
    args <- .qmeth(cls, "_eval_args", args)
    .qraw(cls, .qmeth(cls, "_eval_hilbert_space", args), args)
  }),
  `_eval_args` = .cm(function(cls, args, ...) .qsuper("QExpr", "_eval_args", cls, args)),
  k = .prop(function(self) .qa(self, "label")[[2L]]),
  targets = .prop(function(self) .qa(self, "label")[1L]),
  gate_name_plot = .prop(function(self) sprintf("$%s_%s$", .qa(self, "gate_name_latex"), .str(.qa(self, "k")))),
  get_target_matrix = function(self, format = "sympy") {
    if (!identical(format, "sympy")) .qni("Invalid format for the R_k gate: ", .qrepr(format))
    k <- .qa(self, "k")
    Matrix(list(list(1, 0), list(0, exp(.mul(list(sign(k), S_Two, sym_pi, sym_I, .pow(.pow(S_Two, Abs(k)), S_NegativeOne)))))))
  })

.qclass("Fourier", "Gate",
  `_eval_args` = .cm(function(cls, args, ...) {
    if (length(args) != 2L) .qstop("QuantumError", "QFT/IQFT only takes two arguments, got: ", .qrepr(args))
    a <- lapply(args, as_sym)
    if (!isTRUE(.qgt(a[[2L]], a[[1L]]))) .qstop("QuantumError", "Start must be smaller than finish")
    .qsuper("Gate", "_eval_args", cls, args)
  }),
  `_represent_default_basis` = function(self, ...) .qmeth(self, "_represent_ZGate", NULL, ...),
  `_represent_ZGate` = function(self, basis, nqubits = 0, ...) {
    nqubits <- .qint(nqubits)
    if (nqubits == 0) .qstop("QuantumError", "The number of qubits must be given as nqubits.")
    if (nqubits < .qint(.qa(self, "min_qubits")))
      .qstop("QuantumError", "The number of qubits ", nqubits, " is too small for the gate.")
    size <- .qint(.qa(self, "size"))
    omega <- .qa(self, "omega")
    vals <- vector("list", size * size)
    for (j in seq_len(size) - 1L) for (i in seq_len(size) - 1L)
      vals[[j * size + i + 1L]] <- .mul(list(.pow(omega, as_sym((i * j) %% size)), .pow(as_sym(size), .neg(S_Half))))
    m <- .new_mat(vals, size, size)
    l0 <- .qint(.qa(self, "label")[[1L]])
    if (l0 != 0) m <- matrix_tensor_product(eye(2^l0), m)
    mq <- .qint(.qa(self, "min_qubits"))
    if (mq < nqubits) m <- matrix_tensor_product(m, eye(2^(nqubits - mq)))
    m
  },
  targets = .prop(function(self) {
    l <- .qa(self, "label")
    lapply(seq.int(.qint(l[[1L]]), length.out = .qint(l[[2L]]) - .qint(l[[1L]])), as_sym)
  }),
  min_qubits = .prop(function(self) .qa(self, "label")[[2L]]),
  size = .prop(function(self) .pow(S_Two, .sub(.qa(self, "label")[[2L]], .qa(self, "label")[[1L]]))),
  omega = .prop(function(self) Symbol("omega")))

.qclass("QFT", "Fourier",
  gate_name = "QFT", gate_name_latex = "QFT",
  decompose = function(self) {
    start <- .qint(.qa(self, "label")[[1L]]); finish <- .qint(.qa(self, "label")[[2L]])
    circuit <- S_One
    for (level in rev(seq.int(start, length.out = finish - start))) {
      circuit <- .mul(list(HadamardGate(level), circuit))
      for (i in seq_len(level - start) - 1L)
        circuit <- .mul(list(CGate(level - i - 1L, RkGate(level, i + 2L)), circuit))
    }
    for (i in seq_len((finish - start) %/% 2) - 1L)
      circuit <- .mul(list(SwapGate(i + start, finish - i - 1L), circuit))
    circuit
  },
  `_apply_operator_Qubit` = function(self, qubits, ...) qapply(.mul(list(.qmeth(self, "decompose"), qubits))),
  `_eval_inverse` = function(self) .qfunc("IQFT", self$args),
  omega = .prop(function(self) exp(.mul(list(S_Two, sym_pi, sym_I, .pow(.qa(self, "size"), S_NegativeOne))))))

.qclass("IQFT", "Fourier",
  gate_name = "IQFT", gate_name_latex = "{QFT^{-1}}",
  decompose = function(self) {
    start <- .qint(self$args[[1L]]); finish <- .qint(self$args[[2L]])
    circuit <- S_One
    for (i in seq_len((finish - start) %/% 2) - 1L)
      circuit <- .mul(list(SwapGate(i + start, finish - i - 1L), circuit))
    for (level in seq.int(start, length.out = finish - start)) {
      for (i in rev(seq_len(level - start) - 1L))
        circuit <- .mul(list(CGate(level - i - 1L, RkGate(level, -i - 2L)), circuit))
      circuit <- .mul(list(HadamardGate(level), circuit))
    }
    circuit
  },
  `_eval_inverse` = function(self) .qfunc("QFT", self$args),
  omega = .prop(function(self) exp(.mul(list(as_sym(-2L), sym_pi, sym_I, .pow(.qa(self, "size"), S_NegativeOne))))))

#' Quantum Fourier transform
#'
#' `QFT(start, finish)` and `IQFT(start, finish)` are the quantum Fourier
#' transform and its inverse on the qubits `start` to `finish - 1`, and
#' `RkGate(target, k)` (also `Rk()`) the phase gate they decompose into, as
#' in `sympy.physics.quantum.qft`.
#'
#' @param ... The qubit range, or the target and `k`.
#' @return A gate.
#' @examples
#' qcall(QFT(0, 3), "decompose")
#' represent(QFT(0, 2), nqubits = 2)
#' RkGate(1, 2)
#' @name qft
NULL

#' @rdname qft
#' @export
QFT <- function(...) .qfunc("QFT", list(...))
#' @rdname qft
#' @export
IQFT <- function(...) .qfunc("IQFT", list(...))
#' @rdname qft
#' @export
RkGate <- function(...) .qfunc("RkGate", list(...))
#' @rdname qft
#' @export
Rk <- RkGate

# ---------------------------------------------------------------------------
# Grover

#' Grover's search
#'
#' `superposition_basis(nqubits)` is the equal superposition of all basis
#' states, `OracleGate(nqubits, f)` flips the sign of the states for which
#' the R function `f` is `TRUE`, `WGate(nqubits)` is the Grover diffusion
#' operator, and `grover_iteration()` and `apply_grover()` run the
#' algorithm, as in `sympy.physics.quantum.grover`.
#'
#' @param nqubits The number of qubits.
#' @param ... For `OracleGate()`, the number of qubits and the search
#'   function; for `WGate()`, the number of qubits.
#' @param qstate A state.
#' @param oracle An `OracleGate()` or, for `apply_grover()`, a search function.
#' @param iterations The number of iterations.
#' @return A state or a gate.
#' @examples
#' f <- function(q) isTRUE(q == IntQubit(2, nqubits = 2))
#' apply_grover(f, 2)
#' @name grover
NULL

#' @rdname grover
#' @export
superposition_basis <- function(nqubits) {
  nqubits <- .qint(nqubits)
  amp <- .pow(.pow(S_Two, as_sym(nqubits)), .neg(S_Half))
  .add(lapply(seq_len(2^nqubits) - 1L, function(n) .mul(list(amp, IntQubit(n, nqubits = nqubits)))))
}

.qoracle_function <- function(f) {
  if (is_sym(f) && f$type == "OracleGateFunction") return(f)
  if (!is.function(f)) .qstop("TypeError", "Callable expected, got: ", .qrepr(f))
  id <- paste(c(deparse(f), utils::capture.output(print(environment(f)))), collapse = "\n")
  new_sym("OracleGateFunction", list(), extra = list(fn = f), key = paste0("OracleGateFunction:", id))
}

.qclass("OracleGateFunction", "Basic",
  `__call__` = function(self, ...) self$extra$fn(...),
  `_sympystr` = function(self, p) "OracleGateFunction",
  `_sympyrepr` = function(self, p) "OracleGateFunction()")
# an atom that holds the R function: a rebuild keeps it
register_rebuilder("OracleGateFunction", function(e, args) e)

.qclass("OracleGate", "Gate",
  gate_name = "V", gate_name_latex = "V",
  `_eval_args` = .cm(function(cls, args, ...) {
    if (length(args) != 2L)
      .qstop("QuantumError", "Insufficient/excessive arguments to Oracle.  Please supply the number of qubits and an unknown function.")
    sub <- .qsuper("UnitaryOperator", "_eval_args", cls, args[1L])
    if (sub[[1L]]$type != "Integer") .qstop("TypeError", "Integer expected, got: ", .qrepr(sub[[1L]]))
    list(sub[[1L]], .qoracle_function(args[[2L]]))
  }),
  `_eval_hilbert_space` = .cm(function(cls, args) .qgate_hs(args[[1L]])),
  search_function = .prop(function(self) .qa(self, "label")[[2L]]),
  targets = .prop(function(self) lapply(seq_len(.qint(self$args[[1L]])) - 1L, as_sym)),
  `_apply_operator_Qubit` = function(self, qubits, ...) {
    if (.qint(.qa(qubits, "nqubits")) != .qint(.qa(self, "nqubits")))
      .qstop("QuantumError", "OracleGate operates on ", .qint(.qa(self, "nqubits")), " qubits, got: ", .qint(.qa(qubits, "nqubits")))
    if (.qtruthy(.qmeth(.qa(self, "search_function"), "__call__", qubits))) .neg(qubits) else qubits
  },
  `_represent_ZGate` = function(self, basis, ...) {
    nq <- .qint(.qa(self, "nqubits"))
    nbasis <- 2^nq
    vals <- vector("list", nbasis * nbasis)
    for (i in seq_len(nbasis)) for (j in seq_len(nbasis)) vals[[(i - 1L) * nbasis + j]] <- if (i == j) S_One else S_Zero
    for (i in seq_len(nbasis))
      if (.qtruthy(.qmeth(.qa(self, "search_function"), "__call__", IntQubit(i - 1L, nqubits = nq))))
        vals[[(i - 1L) * nbasis + i]] <- S_NegativeOne
    .new_mat(vals, nbasis, nbasis)
  })

.qclass("WGate", "Gate",
  gate_name = "W", gate_name_latex = "W",
  `_eval_args` = .cm(function(cls, args, ...) {
    if (length(args) != 1L)
      .qstop("QuantumError", "Insufficient/excessive arguments to W gate.  Please supply the number of qubits to operate on.")
    args <- .qsuper("UnitaryOperator", "_eval_args", cls, args)
    if (args[[1L]]$type != "Integer") .qstop("TypeError", "Integer expected, got: ", .qrepr(args[[1L]]))
    args
  }),
  targets = .prop(function(self) lapply(rev(seq_len(.qint(self$args[[1L]])) - 1L), as_sym)),
  `_apply_operator_Qubit` = function(self, qubits, ...) {
    nq <- .qint(.qa(self, "nqubits"))
    if (.qint(.qa(qubits, "nqubits")) != nq)
      .qstop("QuantumError", "WGate operates on ", nq, " qubits, got: ", .qint(.qa(qubits, "nqubits")))
    basis_states <- superposition_basis(nq)
    change <- .mul(list(S_Two, .pow(.pow(S_Two, as_sym(nq)), .neg(S_Half)), basis_states))
    .sub(change, qubits)
  })

#' @rdname grover
#' @export
OracleGate <- function(...) .qfunc("OracleGate", list(...))
#' @rdname grover
#' @export
WGate <- function(...) .qfunc("WGate", list(...))

#' @rdname grover
#' @export
grover_iteration <- function(qstate, oracle) {
  wgate <- WGate(.qa(oracle, "nqubits"))
  .mul(list(wgate, oracle, qstate))
}

#' @rdname grover
#' @export
apply_grover <- function(oracle, nqubits, iterations = NULL) {
  nqubits <- .qint(nqubits)
  if (nqubits <= 0) .qstop("QuantumError", "Grover's algorithm needs nqubits > 0, received ", nqubits, " qubits")
  if (is.null(iterations)) iterations <- fcall("floor", list(.mul(list(.sqrt(.pow(S_Two, as_sym(nqubits))), sym_pi, S_Half, S_Half))))
  v <- OracleGate(nqubits, oracle)
  iterated <- superposition_basis(nqubits)
  for (it in seq_len(.qint(iterations))) iterated <- qapply(grover_iteration(iterated, v))
  iterated
}
