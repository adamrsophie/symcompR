# Qubits and their measurement, after sympy/physics/quantum/qubit.py.

.qbit01 <- function(qb) {
  if (is.character(qb)) {
    if (qb == "0") return(S_Zero)
    if (qb == "1") return(S_One)
    return(Symbol(qb))
  }
  as_sym(qb)
}

.qubit_eval_args <- function(cls, args, ...) {
  if (length(args) == 1L && is_sym(args[[1L]]) && .qisa(args[[1L]], "QubitState")) return(.qa(args[[1L]], "qubit_values"))
  if (length(args) == 1L && is.character(args[[1L]])) {
    args <- lapply(strsplit(args[[1L]], "")[[1L]], function(qb) if (qb == "0") S_Zero else S_One)
  } else {
    args <- lapply(args, .qbit01)
  }
  for (el in args)
    if (!(.qeq(el, S_Zero) || .qeq(el, S_One))) .qstop("ValueError", "Qubit values must be 0 or 1, got: ", .qrepr(el))
  args
}

.qclass("QubitState", "State",
  `_eval_args` = .cm(.qubit_eval_args),
  `_eval_hilbert_space` = .cm(function(cls, args) TensorPowerHilbertSpace(ComplexSpace(2), length(args))),
  dimension = .prop(function(self) length(.qa(self, "qubit_values"))),
  nqubits = .prop(function(self) .qa(self, "dimension")),
  qubit_values = .prop(function(self) .qa(self, "label")),
  `__len__` = function(self) .qa(self, "dimension"),
  `__getitem__` = function(self, bit) {
    i <- .qa(self, "dimension") - .qint(bit) - 1
    v <- .qa(self, "qubit_values")
    if (i < 0 || i >= length(v)) .qstop("IndexError", "tuple index out of range")
    v[[i + 1]]
  },
  flip = function(self, ...) {
    newargs <- .qa(self, "qubit_values")
    for (i in list(...)) {
      bit <- .qa(self, "dimension") - .qint(i)
      newargs[[bit]] <- if (.qeq(newargs[[bit]], S_One)) S_Zero else S_One
    }
    .qfunc(self$type, newargs)
  })

.qclass("Qubit", c("QubitState", "Ket"),
  dual_class = .cm(function(cls) "QubitBra"),
  `_eval_innerproduct_QubitBra` = function(self, bra, ...) {
    if (identical(vapply(.qa(self, "label"), .k, ""), vapply(.qa(bra, "label"), .k, ""))) S_One else S_Zero
  },
  `_represent_default_basis` = function(self, ...) .qmeth(self, "_represent_ZGate", NULL, ...),
  `_represent_ZGate` = function(self, basis, format = "sympy", ...) {
    n <- 1; definite_state <- 0
    for (it in rev(.qa(self, "qubit_values"))) {
      definite_state <- definite_state + n * .qint(it)
      n <- n * 2
    }
    result <- rep(list(S_Zero), 2^.qa(self, "dimension"))
    result[[definite_state + 1]] <- S_One
    if (!identical(format, "sympy")) .qstop("NotImplementedError", "The ", format, " format is not available in R.")
    Matrix(result)
  },
  `_eval_trace` = function(self, bra, indices = list(), ...) {
    sorted_idx <- vapply(indices, .qint, 0)
    if (!length(sorted_idx)) sorted_idx <- seq_len(.qa(self, "nqubits")) - 1
    sorted_idx <- sort(sorted_idx)
    new_mat <- .mul(list(self, bra))
    for (i in rev(seq_along(sorted_idx))) new_mat <- .qmeth(self, "_reduced_density", new_mat, sorted_idx[i])
    if (length(sorted_idx) == .qa(self, "nqubits")) return(new_mat$data[[1L]])
    matrix_to_density(new_mat)
  },
  `_reduced_density` = function(self, matrix, qubit, ...) {
    find_index_that_is_projected <- function(j, k, qubit) {
      bit_mask <- 2^qubit - 1
      bitwShiftL(bitwShiftR(j, qubit), 1 + qubit) + bitwAnd(j, bit_mask) + bitwShiftL(k, qubit)
    }
    old_matrix <- represent(matrix, ...)
    old_size <- old_matrix$ncol
    new_size <- old_size %/% 2
    vals <- vector("list", new_size * new_size)
    for (i in seq_len(new_size) - 1L) for (j in seq_len(new_size) - 1L) {
      acc <- S_Zero
      for (k in 0:1) {
        col <- find_index_that_is_projected(j, k, qubit)
        row <- find_index_that_is_projected(i, k, qubit)
        acc <- .add(list(acc, .mget(old_matrix, row + 1L, col + 1L)))
      }
      vals[[i * new_size + j + 1L]] <- acc
    }
    .new_mat(vals, new_size, new_size)
  })

.qclass("QubitBra", c("QubitState", "Bra"), dual_class = .cm(function(cls) "Qubit"))

# type(a) of Python, as printed in error messages
.qpytype <- function(a) {
  if (is.character(a)) return("<class 'str'>")
  if (inherits(a, "symmatrix")) return("<class 'sympy.matrices.dense.MutableDenseMatrix'>")
  if (is.numeric(a)) return(if (a == round(a)) "<class 'int'>" else "<class 'float'>")
  if (is_sym(a) && a$type == "Integer") return("<class 'int'>")
  if (is_sym(a) && a$type == "Float") return("<class 'float'>")
  if (is_sym(a)) return(paste0("<class '", .qclassname(a), "'>"))
  paste0("<class '", class(a)[1L], "'>")
}
.qis_int <- function(a) (is.numeric(a) && length(a) == 1L && a == round(a)) || (is_sym(a) && a$type == "Integer")

.qclass("IntQubitState", "QubitState",
  `_eval_args` = .cm(function(cls, args, nqubits = NULL, ...) {
    extra <- names(list(...))
    if (length(extra)) .qstop("TypeError", "IntQubitState._eval_args() got an unexpected keyword argument '", extra[1L], "'")
    if (length(args) == 1L && is_sym(args[[1L]]) && .qisa(args[[1L]], "QubitState")) return(.qubit_eval_args(cls, args))
    if (!all(vapply(args, .qis_int, NA)))
      .qstop("ValueError", "values must be integers, got (", .qtuple_str(vapply(args, .qpytype, "")), ")")
    if (!is.null(nqubits)) {
      if (!.qis_int(nqubits)) .qstop("ValueError", "nqubits must be an integer, got (", .qpytype(nqubits), ")")
      if (length(args) != 1L)
        .qstop("ValueError", "too many positional arguments (", .qrepr(args), "). should be (number, nqubits=n)")
      return(.qmeth(cls, "_eval_args_with_nqubits", args[[1L]], nqubits))
    }
    if (length(args) == 1L && .qint(args[[1L]]) > 1) {
      v <- .qint(args[[1L]])
      nb <- floor(log2(v)) + 1
      return(.qubit_eval_args(cls, lapply(rev(seq_len(nb) - 1), function(i) as_sym((v %/% 2^i) %% 2))))
    } else if (length(args) == 2L && .qint(args[[2L]]) > 1) {
      return(.qmeth(cls, "_eval_args_with_nqubits", args[[1L]], args[[2L]]))
    }
    .qubit_eval_args(cls, args)
  }),
  `_eval_args_with_nqubits` = .cm(function(cls, number, nqubits) {
    v <- .qint(number); nq <- .qint(nqubits)
    need <- if (v == 0) 0 else floor(log2(v)) + 1
    if (nq < need) .qstop("ValueError", "cannot represent ", v, " with ", nq, " bits")
    .qubit_eval_args(cls, lapply(rev(seq_len(nq) - 1), function(i) as_sym((v %/% 2^i) %% 2)))
  }),
  as_int = function(self) {
    number <- 0; n <- 1
    for (i in rev(.qa(self, "qubit_values"))) { number <- number + n * .qint(i); n <- n * 2 }
    number
  },
  `_print_label` = function(self, p) format(.qmeth(self, "as_int"), scientific = FALSE),
  `_print_label_repr` = function(self, p) format(.qmeth(self, "as_int"), scientific = FALSE),
  `_print_label_latex` = function(self, p) format(.qmeth(self, "as_int"), scientific = FALSE))

.qtuple_str <- function(x) paste0("(", paste(x, collapse = ", "), if (length(x) == 1L) ",", ")")

.qclass("IntQubit", c("IntQubitState", "Qubit"),
  dual_class = .cm(function(cls) "IntQubitBra"),
  `_eval_innerproduct_IntQubitBra` = function(self, bra, ...) .qsuper("Qubit", "_eval_innerproduct_QubitBra", self, bra))

.qclass("IntQubitBra", c("IntQubitState", "QubitBra"), dual_class = .cm(function(cls) "IntQubit"))

#' Qubits
#'
#' The qubit states of `sympy.physics.quantum.qubit`: `Qubit("0101")` or
#' `Qubit(0, 1, 0, 1)`, its bra `QubitBra()`, and `IntQubit(n)` and
#' `IntQubitBra(n)` labeled by an integer (`nqubits` fixes the width). Gates
#' act on them through [qapply()]. `qubit_to_matrix()` and
#' `matrix_to_qubit()` convert to and from column vectors,
#' `matrix_to_density()` gives a [Density()] and `measure_all()`,
#' `measure_partial()`, `measure_all_oneshot()` and
#' `measure_partial_oneshot()` give measurement outcomes with their
#' probabilities (the one shot versions draw a random outcome).
#'
#' @param ... Qubit values, a string of zeros and ones, or a qubit.
#' @param nqubits The number of qubits of an integer qubit.
#' @param qubit A qubit expression.
#' @param matrix,mat A matrix.
#' @param bits The qubits to measure.
#' @param format The matrix format, only `"sympy"`.
#' @param normalize Normalize the state first.
#' @return A qubit expression, a matrix or a list of outcomes.
#' @examples
#' q <- Qubit("01")
#' qubit_to_matrix(q)
#' matrix_to_qubit(qubit_to_matrix(q))
#' measure_all(HadamardGate(0) * Qubit("0"))
#' @name quantum-qubits
NULL

#' @rdname quantum-qubits
#' @export
QubitState <- function(...) .qfunc("QubitState", list(...))
#' @rdname quantum-qubits
#' @export
Qubit <- function(...) .qfunc("Qubit", list(...))
#' @rdname quantum-qubits
#' @export
QubitBra <- function(...) .qfunc("QubitBra", list(...))
#' @rdname quantum-qubits
#' @export
IntQubit <- function(..., nqubits = NULL) .qint_qubit("IntQubit", list(...), nqubits)

# named arguments other than nqubits are keywords of _eval_args
.qint_qubit <- function(cls, a, nqubits) {
  kw <- if (is.null(names(a))) logical(length(a)) else nzchar(names(a))
  do.call(.qfunc, c(list(cls, unname(a[!kw]), nqubits = nqubits), a[kw]))
}
#' @rdname quantum-qubits
#' @export
IntQubitBra <- function(..., nqubits = NULL) .qint_qubit("IntQubitBra", list(...), nqubits)

#' @rdname quantum-qubits
#' @export
matrix_to_qubit <- function(matrix) {
  if (matrix$nrow == 1L) {
    mlistlen <- matrix$ncol; ket <- FALSE; cls <- "QubitBra"
  } else if (matrix$ncol == 1L) {
    mlistlen <- matrix$nrow; ket <- TRUE; cls <- "Qubit"
  } else {
    .qstop("QuantumError", "Matrix must be a row/column vector, got ", .qmat_str(matrix))
  }
  nqubits <- log2(mlistlen)
  if (nqubits != round(nqubits))
    .qstop("QuantumError", "Matrix must be a row/column vector of size 2**nqubits, got: ", .qmat_str(matrix))
  result <- S_Zero
  for (i in seq_len(mlistlen) - 1L) {
    element <- matrix$data[[i + 1L]]
    if (.qtruthy(element)) {
      qubit_array <- lapply(rev(seq_len(nqubits) - 1L), function(x) as_sym(as.integer(bitwAnd(i, bitwShiftL(1L, x)) != 0)))
      result <- .add(list(result, .mul(list(element, .qfunc(cls, qubit_array)))))
    }
  }
  if (result$type %in% c("Mul", "Add", "Pow")) result <- expand(result)
  result
}

#' @rdname quantum-qubits
#' @export
matrix_to_density <- function(mat) {
  eigen <- eigenvects(mat)
  args <- list()
  for (x in eigen) {
    if (.qeq(x[[1L]], S_Zero)) next
    for (vector in x[[3L]]) args[[length(args) + 1L]] <- list(matrix_to_qubit(.new_mat(vector$data, length(vector$data), 1L)), x[[1L]])
  }
  if (!length(args)) return(S_Zero)
  do.call(Density, args)
}

#' @rdname quantum-qubits
#' @export
qubit_to_matrix <- function(qubit, format = "sympy") represent(qubit, format = format)

.qmat_normalized <- function(m) {
  n <- .sqrt(.add(lapply(m$data, function(z) .pow(abs(z), S_Two))))
  if (.qeq(n, S_Zero)) return(m)
  .mat_map(m, function(z) .mul(list(z, .pow(n, S_NegativeOne))))
}

#' @rdname quantum-qubits
#' @export
measure_all <- function(qubit, format = "sympy", normalize = TRUE) {
  m <- qubit_to_matrix(qubit, format)
  if (!identical(format, "sympy")) .qni("This function cannot handle non-SymPy matrix formats yet")
  results <- list()
  if (normalize) m <- .qmat_normalized(m)
  size <- max(m$nrow, m$ncol)
  nqubits <- as.integer(floor(log(size) / log(2) + 1e-12))
  for (i in seq_len(size) - 1L) {
    mi <- m$data[[i + 1L]]
    if (.qtruthy(mi))
      results[[length(results) + 1L]] <- list(Qubit(IntQubit(i, nqubits = nqubits)), .mul(list(mi, fcall("conjugate", list(mi)))))
  }
  results
}

#' @rdname quantum-qubits
#' @export
measure_partial <- function(qubit, bits, format = "sympy", normalize = TRUE) {
  m <- qubit_to_matrix(qubit, format)
  bits <- if (.qis_seq(bits)) vapply(.qseq_of(bits), .qint, 0) else .qint(bits)
  if (!identical(format, "sympy")) .qni("This function cannot handle non-SymPy matrix formats yet")
  if (normalize) m <- .qmat_normalized(m)
  output <- list()
  for (outcome in .qget_possible_outcomes(m, bits)) {
    prob_of_outcome <- .add(lapply(outcome$data, function(z) .mul(list(fcall("conjugate", list(z)), z))))
    if (!.qeq(prob_of_outcome, S_Zero)) {
      next_matrix <- matrix_to_qubit(if (normalize) .qmat_normalized(outcome) else outcome)
      output[[length(output) + 1L]] <- list(next_matrix, prob_of_outcome)
    }
  }
  output
}

#' @rdname quantum-qubits
#' @export
measure_partial_oneshot <- function(qubit, bits, format = "sympy") {
  m <- qubit_to_matrix(qubit, format)
  bits <- if (.qis_seq(bits)) vapply(.qseq_of(bits), .qint, 0) else .qint(bits)
  if (!identical(format, "sympy")) .qni("This function cannot handle non-SymPy matrix formats yet")
  m <- .qmat_normalized(m)
  random_number <- stats::runif(1L)
  total_prob <- 0
  for (outcome in .qget_possible_outcomes(m, bits)) {
    total_prob <- total_prob + as.numeric(N(.add(lapply(outcome$data, function(z) .mul(list(fcall("conjugate", list(z)), z))))))
    if (total_prob >= random_number) return(matrix_to_qubit(.qmat_normalized(outcome)))
  }
  NULL
}

.qget_possible_outcomes <- function(m, bits) {
  size <- max(m$nrow, m$ncol)
  nqubits <- as.integer(floor(log2(size) + 0.1))
  output <- rep(list(rep(list(S_Zero), 2^nqubits)), 2^length(bits))
  bit_masks <- bitwShiftL(1L, as.integer(bits))
  for (i in seq_len(2^nqubits) - 1L) {
    trueness <- 0
    for (j in seq_along(bit_masks)) if (bitwAnd(i, bit_masks[j])) trueness <- trueness + j
    output[[trueness + 1]][[i + 1L]] <- m$data[[i + 1L]]
  }
  lapply(output, function(v) .new_mat(v, length(v), 1L))
}

#' @rdname quantum-qubits
#' @export
measure_all_oneshot <- function(qubit, format = "sympy") {
  m <- qubit_to_matrix(qubit)
  if (!identical(format, "sympy")) .qni("This function cannot handle non-SymPy matrix formats yet")
  m <- .qmat_normalized(m)
  random_number <- stats::runif(1L)
  total <- 0; result <- 0L
  for (i in m$data) {
    total <- total + as.numeric(N(.mul(list(i, fcall("conjugate", list(i))))))
    if (total > random_number) break
    result <- result + 1L
  }
  Qubit(IntQubit(result, nqubits = as.integer(floor(log2(max(m$nrow, m$ncol)) + 0.1))))
}
