# The trace Tr (sympy/physics/quantum/trace.py) and density operators,
# entropy and fidelity (sympy/physics/quantum/density.py).

.qis_scalar <- function(e) {
  e <- as_sym(e)
  e$type %in% c("Integer", "Rational", "Float") || is_number_type(e) ||
    (e$type %in% c("Symbol", "Dummy") && isTRUE(ask_fact(e, "commutative")))
}

# rotate a list so it starts at its smallest element (by sort key); among
# equal smallest elements the smallest rotation wins
.qcycle_permute <- function(l) {
  n <- length(l)
  if (n == 1L) return(l)
  keys <- vapply(l, .k, "")
  mk <- .k(.qmin(l))
  idx <- which(keys == mk)
  rots <- lapply(idx, function(i) c(l[i:n], l[seq_len(i - 1L)]))
  if (length(rots) == 1L) return(rots[[1L]])
  .key_sort(lapply(rots, function(r) lapply(r, .sort_key)), rots)[[1L]]
}

.qrearrange_args <- function(l) {
  if (length(l) == 1L) return(l)
  .qargs(.mul(c(l[length(l)], l[-length(l)])))
}

.qtr_new <- function(expr, indices) {
  args <- list(expr, indices)
  content <- if (expr$type == "Mul") .qcycle_permute(.qrearrange_args(.qargs(expr))) else list(expr)
  new_sym("Tr", args, key = .compound_key("Tr", c(content, list(indices))))
}

.qtr <- function(args) {
  if (length(args) == 2L) {
    ix <- args[[2L]]
    indices <- if (is_sym(ix) && ix$type == "Tuple") ix
               else if (is.list(ix) && !is_sym(ix)) do.call(Tuple, ix)
               else Tuple(ix)
    expr <- args[[1L]]
  } else if (length(args) == 1L) {
    indices <- Tuple()
    expr <- args[[1L]]
  } else {
    .qstop("ValueError", "Arguments to Tr should be of form (expr[, [indices]])")
  }
  if (inherits(expr, "symmatrix")) return(sym_trace(expr))
  if (is_sym(expr) && expr$type == "ImmutableMatrix") return(sym_trace(as_explicit(expr)))
  # any other object with a trace() method
  if (is_sym(expr) && !is.null(.qnode_def(expr$type)) && .qhas(expr, "trace")) return(.qmeth(expr, "trace"))
  expr <- as_sym(expr)
  if (expr$type == "Add") return(.add(lapply(.qargs(expr), function(a) .qtr(list(a, indices)))))
  if (expr$type == "Mul") {
    cnc <- .args_cnc(expr)
    if (!length(cnc[[2L]])) return(.mul(cnc[[1L]]))
    obj <- .qtr_new(.mul(cnc[[2L]]), indices)
    return(if (length(cnc[[1L]])) .mul(c(cnc[[1L]], list(obj))) else obj)
  }
  if (expr$type == "Pow") {
    if (.qis_scalar(expr$args[[1L]]) && .qis_scalar(expr$args[[2L]])) return(expr)
    return(.qtr_new(expr, indices))
  }
  if (.qis_scalar(expr)) return(expr)
  .qtr_new(expr, indices)
}

do.call(.qclass, c(list("Tr", "Expr",
  `__new__` = .cm(function(cls, args, ...) .qtr(args)),
  is_number = .prop(function(self) TRUE),
  doit = function(self, ...) {
    a <- self$args[[1L]]
    if (!is.null(.qnode_def(a$type)) && .qhas(a, "_eval_trace")) return(.qmeth(a, "_eval_trace", indices = self$args[[2L]]$args))
    self
  },
  permute = function(self, pos) {
    pos <- .qint(pos)
    a <- .qargs(self$args[[1L]])
    n <- length(a)
    start <- (-pos) %% n
    Tr(.mul(c(a[seq.int(start + 1L, length.out = n - start)], a[seq_len(start)])))
  },
  `_sympystr` = function(self, p) paste0("Tr(", .str(self$args[[1L]]), ")"),
  `_latex` = function(self, p) paste0("\\operatorname{tr}\\left(", .latex(self$args[[1L]]), "\\right)")),
  .qbasic_printers[2L]))
# Tr.doit traces its argument as it stands
register_doit("Tr", structure(function(e, deep = TRUE) .qmeth(e, "doit"), shallow = TRUE))

#' Trace of a quantum expression
#'
#' `Tr(expr, indices)` is the trace of `sympy.physics.quantum.trace`:
#' scalars come out of it, a sum gives a sum of traces, a matrix gives its
#' trace, and [doit()] evaluates traces of density operators, outer
#' products and tensor products (a partial trace over `indices`).
#'
#' @param ... The expression and the optional indices to trace over.
#' @return An expression.
#' @examples
#' A <- Symbol("A", commutative = FALSE)
#' a <- Symbol("a")
#' Tr(a * A)
#' Tr(Matrix(list(list(1, 1), list(2, 2))))
#' @export
Tr <- function(...) .qtr(list(...))

# ---------------------------------------------------------------------------
# Density

# the node type is QDensity: sympy.stats has its own Density
.qclass("QDensity", "HermitianOperator", pyname = "Density",
  `_eval_args` = .cm(function(cls, args, ...) {
    args <- .qsuper("HermitianOperator", "_eval_args", cls, args, ...)
    for (a in args)
      if (!(is_sym(a) && a$type == "Tuple" && length(a$args) == 2L))
        .qstop("ValueError", "Each argument should be of form [state,prob] or ( state, prob )")
    args
  }),
  states = function(self) do.call(Tuple, lapply(self$args, function(a) a$args[[1L]])),
  probs = function(self) do.call(Tuple, lapply(self$args, function(a) a$args[[2L]])),
  get_state = function(self, index) self$args[[.qint(index) + 1L]]$args[[1L]],
  get_prob = function(self, index) self$args[[.qint(index) + 1L]]$args[[2L]],
  apply_op = function(self, op) .qfunc("QDensity", lapply(self$args, function(a) list(.mul(list(op, a$args[[1L]])), a$args[[2L]]))),
  doit = function(self, ...) {
    terms <- list()
    for (a in self$args) {
      state <- expand(a$args[[1L]]); prob <- a$args[[2L]]
      if (state$type == "Add") {
        sa <- .qargs(state)
        for (x in sa) for (y in sa) terms[[length(terms) + 1L]] <- .mul(list(prob, .qmeth(self, "_generate_outer_prod", x, y)))
      } else {
        terms[[length(terms) + 1L]] <- .mul(list(prob, .qmeth(self, "_generate_outer_prod", state, state)))
      }
    }
    .add(terms)
  },
  `_generate_outer_prod` = function(self, arg1, arg2) {
    p1 <- .args_cnc(arg1); p2 <- .args_cnc(arg2)
    if (!length(p1[[2L]]) || !length(p2[[2L]]))
      .qstop("ValueError", "Atleast one-pair of Non-commutative instance required for outer product.")
    # Ket*Bra through the operator is an OuterProduct
    op <- .mul(p1[[2L]]) * Dagger(.mul(p2[[2L]]))
    .mul(list(.mul(p1[[1L]]), .mul(p2[[1L]]), op))
  },
  `_represent` = function(self, ...) represent(doit(self), ...),
  `_print_operator_name_latex` = function(self, p) "\\rho",
  `_eval_trace` = function(self, indices = list(), ...) doit(Tr(doit(self), indices)),
  entropy = function(self) .qentropy(self))

#' Density operators, entropy and fidelity
#'
#' `Density(list(state, prob), ...)` is the density operator of
#' `sympy.physics.quantum.density`: [doit()] expands it into outer
#' products, [represent()] gives its matrix, and [Tr()] traces it.
#' `Density()` and `entropy()` are shared with sympy.stats (see [Density()]
#' and [entropy()]); they act as here when the arguments are
#' `list(state, probability)` pairs, a quantum density operator or a matrix.
#' `entropy()` is then the von Neumann entropy. `fidelity()` is the fidelity
#' of two density operators or matrices.
#'
#' @param state1,state2 Density operators or matrices.
#' @return A density operator or an expression.
#' @examples
#' d <- Density(list(JzKet(1, 0), 0.5), list(JzKet(1, 1), 0.5))
#' doit(d)
#' entropy(Density(list(JzKet(S(1)/2, S(1)/2), S(1)/2), list(JzKet(S(1)/2, -S(1)/2), S(1)/2)))
#' @name density-operators
NULL

# Density() and entropy() of sympy.stats (R/stats-aa-rv.R and
# R/stats-symbolic.R) call these for quantum arguments
.qdensity <- function(...) .qfunc("QDensity", list(...))

.qentropy <- function(density) {
  if (is_sym(density) && .qisa(density, "QDensity")) density <- represent(density)
  if (!inherits(density, "symmatrix")) .qstop("ValueError", "numpy.ndarray, scipy.sparse or SymPy matrix expected")
  ev <- .qeigenvals(density)
  expand(.neg(.add(lapply(ev, function(e) .mul(list(e, log(e)))))))
}

# the distinct eigenvalues of a matrix
.qeigenvals <- function(M) {
  ev <- eigenvals(M)
  if (is.list(ev) && !is_sym(ev)) {
    vals <- lapply(ev, function(z) if (is_sym(z)) z else z[[1L]])
    keys <- vapply(vals, .k, "")
    return(vals[!duplicated(keys)])
  }
  list(ev)
}

#' @rdname density-operators
#' @export
fidelity <- function(state1, state2) {
  if (is_sym(state1) && .qisa(state1, "QDensity")) state1 <- represent(state1)
  if (is_sym(state2) && .qisa(state2, "QDensity")) state2 <- represent(state2)
  if (!inherits(state1, "symmatrix") || !inherits(state2, "symmatrix"))
    .qstop("ValueError", sprintf("state1 and state2 must be of type Density or Matrix received type=%s for state1 and type=%s for state2",
                                 .qpytype(state1), .qpytype(state2)))
  sq <- state1$nrow == state1$ncol
  if ((state1$nrow != state2$nrow || state1$ncol != state2$ncol) && sq)
    .qstop("ValueError", "The dimensions of both args should be equal and the matrix obtained should be a square matrix")
  if (!sq) .qstop("NonSquareMatrixError", "")
  sqrt_state1 <- .qmat_pow(state1, S_Half)
  m <- .matmul(.matmul(sqrt_state1, state2), sqrt_state1)
  doit(Tr(.qmat_pow(m, S_Half)))
}

# MatrixBase.pow: a diagonal matrix raises its diagonal entries, any other
# matrix goes through its Jordan form
.qmat_pow <- function(M, ex) {
  n <- M$nrow
  off <- unlist(lapply(seq_len(n), function(i) lapply(setdiff(seq_len(n), i), function(j) .mget(M, i, j))), recursive = FALSE)
  if (!all(vapply(off, function(z) isTRUE(ask_fact(z, "zero")), NA))) {
    r <- .qmat_pow_float(M, ex)
    return(if (is.null(r)) M %^% ex else r)
  }
  out <- M
  for (i in seq_len(n)) out$data[[(i - 1L) * n + i]] <- .pow(.mget(M, i, i), ex)
  out
}

# A matrix of numbers with a Float among them: the Jordan form behind
# Matrix.__pow__ is a numerical eigen decomposition, so the power is taken in
# floating point (NULL for a matrix that is not of that kind)
.qmat_pow_float <- function(M, ex) {
  if (!all(vapply(M$data, .expr_is_number, NA)) || !any(vapply(M$data, function(e) .has_type(e, "Float"), NA)))
    return(NULL)
  z <- tryCatch(vapply(M$data, function(e) as.complex(.neval(e, list())), 0i), error = function(e) NULL)
  p <- tryCatch(as.numeric(.neval(ex, list())), error = function(e) NA_real_)
  if (is.null(z) || anyNA(z) || is.na(p)) return(NULL)
  n <- M$nrow
  A <- matrix(z, n, n, byrow = TRUE)
  herm <- isTRUE(all.equal(A, Conj(t(A))))
  e <- eigen(A, symmetric = herm)
  V <- e$vectors
  R <- V %*% diag(as.complex(e$values)^p, n) %*% (if (herm) Conj(t(V)) else solve(V))
  out <- M
  out$data <- lapply(as.vector(t(R)), function(v)
    if (Im(v) == 0) Float(Re(v)) else .add(list(Float(Re(v)), .mul(list(Float(Im(v)), sym_I)))))
  out
}
