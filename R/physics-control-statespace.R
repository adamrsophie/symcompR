# State space models, after sympy/physics/control/lti.py (StateSpaceBase,
# StateSpace and DiscreteStateSpace).

#' State space models
#'
#' `StateSpace(A, B, C, D)` is the continuous state space model
#' `x' = A x + B u`, `y = C x + D u` and `DiscreteStateSpace(A, B, C, D,
#' sampling_time)` the discrete one. The matrices default to zero matrices
#' of the matching shapes. A model has the properties `$A`, `$B`, `$C`, `$D`
#' (also `$state_matrix`, `$input_matrix`, `$output_matrix`,
#' `$feedforward_matrix`), `$num_states`, `$num_inputs`, `$num_outputs` and
#' `$shape`, and the methods `$append()`, `$observability_matrix()`,
#' `$controllability_matrix()`, `$is_observable()`, `$is_controllable()`,
#' the observable and controllable subspaces and forms,
#' `$apply_similarity_transform()`, `$dsolve()`,
#' `$get_asymptotic_stability_conditions()` and
#' `$rewrite("TransferFunction")`. `+`, `-` and `*` connect models.
#' `create_state_space()` picks the class from the sampling time. After
#' `sympy.physics.control.lti`.
#'
#' @param A,B,C,D The state, input, output and feedforward matrices.
#' @param sampling_time The sampling time.
#' @return A `StateSpace` or `DiscreteStateSpace` object.
#' @examples
#' A <- Matrix(list(list(-1, 0), list(0, -2)))
#' ss <- StateSpace(A, Matrix(list(list(1), list(1))), Matrix(list(list(1, 1))), Matrix(list(list(0))))
#' ss$num_states
#' ss$is_controllable()
#' ss$rewrite("TransferFunction")
#' @usage StateSpace(A = NULL, B = NULL, C = NULL, D = NULL)
#' @export
StateSpace <- .pyclass_fn(function(A = NULL, B = NULL, C = NULL, D = NULL)
  .ssb_new("StateSpace", A, B, C, D, extra = list(is_continuous = TRUE)), "StateSpace")

#' @rdname StateSpace
#' @usage DiscreteStateSpace(A = NULL, B = NULL, C = NULL, D = NULL, sampling_time = 1)
#' @export
DiscreteStateSpace <- .pyclass_fn(function(A = NULL, B = NULL, C = NULL, D = NULL, sampling_time = 1L) {
  if (.lti_st_zero(sampling_time))
    .physics_error("ValueError", .py_filldedent("\n                The sampling time cannot be zero.\n                If you want to create a continuous state space,\n                use the StateSpace class instead."))
  sampling_time <- .lti_sympify(sampling_time)
  .ssb_new("DiscreteStateSpace", A, B, C, D, list(sampling_time), extra = list(is_continuous = FALSE))
}, "DiscreteStateSpace")

#' @rdname create_transfer_function
#' @export
create_state_space <- function(A, B, C, D, sampling_time = 0L) {
  if (.lti_st_zero(sampling_time)) return(StateSpace(A, B, C, D))
  DiscreteStateSpace(A, B, C, D, sampling_time)
}

.ss_as_mat <- function(M) {
  if (inherits(M, "symmatrix")) return(M)
  if (is_sym(M) && M$type == "ImmutableMatrix") return(.lti_mat(M))
  NULL
}

# StateSpaceBase.__new__
.ssb_new <- function(cls, A, B, C, D, extra_args = list(), extra = list()) {
  if (is.null(A)) A <- zeros(1L)
  A <- .ss_as_mat(A) %||% A
  if (is.null(B)) B <- zeros(A$nrow, 1L)
  B <- .ss_as_mat(B) %||% B
  if (is.null(C)) C <- zeros(1L, A$ncol)
  C <- .ss_as_mat(C) %||% C
  if (is.null(D)) D <- zeros(C$nrow, B$ncol)
  D <- .ss_as_mat(D) %||% D
  if (!all(vapply(list(A, B, C, D), function(m) inherits(m, "symmatrix"), TRUE)))
    .physics_error("TypeError", "A, B, C and D inputs must all be sympy Matrices.")
  if (A$nrow != A$ncol) .physics_error("ShapeError", "Matrix A must be a square matrix.")
  if (A$nrow != B$nrow) .physics_error("ShapeError", "Matrices A and B must have the same number of rows.")
  if (C$nrow != D$nrow) .physics_error("ShapeError", "Matrices C and D must have the same number of rows.")
  if (A$ncol != C$ncol) .physics_error("ShapeError", "Matrices A and C must have the same number of columns.")
  if (B$ncol != D$ncol) .physics_error("ShapeError", "Matrices B and D must have the same number of columns.")
  .pyobj_new(cls, c(lapply(list(A, B, C, D), .mat_imm), extra_args),
             extra = c(list(is_StateSpace_object = TRUE, is_SISO = D$nrow == 1L && D$ncol == 1L), extra))
}

# a Python number or Symbol, which the operators treat as a scalar
.ss_scalar <- function(x) (is.numeric(x) && !is.object(x)) || is.complex(x) || (is_sym(x) && x$type == "Symbol")

.ss_add <- function(self, other) {
  if (.ss_scalar(other)) return(create_state_space(self$A, self$B, self$C, self$D + .lti_sympify(other), self$sampling_time))
  if (!inherits(other, "StateSpaceBase")) .physics_error("ValueError", "Addition is only supported for 2 State Space models.")
  if (self$num_inputs != other$num_inputs || self$num_outputs != other$num_outputs)
    .physics_error("ShapeError", "Systems with incompatible inputs and outputs cannot be added.")
  .lti_check_time(list(self, other))
  A <- rbind(cbind(self$A, zeros(self$A$nrow, other$A$ncol)), cbind(zeros(other$A$nrow, self$A$ncol), other$A))
  create_state_space(A, rbind(self$B, other$B), cbind(self$C, other$C), self$D + other$D, self$sampling_time)
}

.ss_mul <- function(self, other) {
  if (.ss_scalar(other)) {
    o <- .lti_sympify(other)
    return(create_state_space(self$A, self$B, self$C * o, self$D * o, self$sampling_time))
  }
  if (!inherits(other, "StateSpaceBase")) .physics_error("ValueError", "Multiplication is only supported for 2 State Space models.")
  if (self$num_inputs != other$num_outputs)
    .physics_error("ShapeError", "Systems with incompatible inputs and outputs cannot be multiplied.")
  .lti_check_time(list(self, other))
  A <- rbind(cbind(other$A, zeros(other$A$nrow, self$A$ncol)), cbind(self$B %*% other$C, self$A))
  B <- rbind(other$B, self$B %*% other$D)
  C <- cbind(self$D %*% other$C, self$C)
  create_state_space(A, B, C, self$D %*% other$D, self$sampling_time)
}

.ss_rmul <- function(self, other) {
  if (.ss_scalar(other)) {
    o <- .lti_sympify(other)
    return(create_state_space(self$A, self$B * o, self$C, self$D * o, self$sampling_time))
  }
  .ss_mul(self, other)
}

.ss_neg <- function(self) create_state_space(self$A, self$B, -self$C, -self$D, self$sampling_time)

.ss_append <- function(self, other) {
  .lti_check_time(list(self, other))
  ns <- self$num_states; ni <- self$num_inputs; no <- self$num_outputs
  n <- ns + other$num_states; m <- ni + other$num_inputs; p <- no + other$num_outputs
  blk <- function(X, Y, nr, nc, r0, c0) {
    Matrix(nr, nc, function(i, j) {
      if (i <= X$nrow && j <= X$ncol && i <= r0 && j <= c0) return(.mget(X, i, j))
      if (i > r0 && j > c0) return(.mget(Y, i - r0, j - c0))
      S_Zero
    })
  }
  create_state_space(blk(self$A, other$A, n, n, ns, ns), blk(self$B, other$B, n, m, ns, ni),
                     blk(self$C, other$C, p, n, no, ns), blk(self$D, other$D, p, m, no, ni), self$sampling_time)
}

.ss_hstack <- function(cols, nrow) if (length(cols)) do.call(cbind, cols) else zeros(nrow, 0L)

.ss_free <- function(M) unique(unlist(lapply(M$data, function(z) vapply(free_symbols(z), .k, ""))))

.ss_observability <- function(self) {
  ob <- self$C
  n <- self$num_states
  if (n > 1L) for (i in seq_len(n - 1L)) ob <- rbind(ob, self$C %*% (self$A %^% i))
  ob
}

.ss_controllability <- function(self) {
  co <- self$B
  n <- self$A$nrow
  if (n > 1L) for (i in seq_len(n - 1L)) co <- cbind(co, (self$A %^% i) %*% self$B)
  co
}

.ss_orth_complement <- function(M, dim) {
  if (M$nrow == 0L) return(columnspace(eye(dim)))
  nullspace(t(M))
}

.ss_unobservable <- function(self) {
  M <- .ss_observability(self)
  if (length(.ss_free(M))) .physics_error("NotImplementedError", "Unbservable subspace cannot be determined for symbolic matrices.")
  nullspace(M)
}

.ss_observable <- function(self) {
  M <- .ss_hstack(.ss_unobservable(self), 0L)
  if (length(.ss_free(M))) .physics_error("NotImplementedError", "Observable subspace cannot be determined for symbolic matrices.")
  .ss_orth_complement(M, self$A$nrow)
}

.ss_controllable <- function(self) columnspace(.ss_controllability(self))

.ss_uncontrollable <- function(self) {
  M <- .ss_hstack(.ss_controllable(self), 0L)
  if (length(.ss_free(M))) .physics_error("NotImplementedError", "Uncontrollable subspace cannot be determined for symbolic matrices.")
  .ss_orth_complement(M, self$A$nrow)
}

.ss_similarity <- function(self, Tm) {
  T_inv <- inverse(Tm)
  StateSpace(T_inv %*% self$A %*% Tm, T_inv %*% self$B, self$C %*% Tm, self$D)
}

.ss_form <- function(self, first, second) {
  a <- .ss_hstack(first(self), self$A$nrow)
  b <- .ss_hstack(second(self), self$A$nrow)
  if (b$ncol == 0L || a$ncol == 0L) return(self)
  .ss_similarity(self, cbind(a, b))
}

# StateSpace._eval_rewrite_as_TransferFunction and the discrete one
.ss_to_tf <- function(self, var, make) {
  n <- self$A$nrow
  G <- self$C %*% solve(var * eye(n) - self$A, self$B) + self$D
  G <- simplify(G)
  lapply(.lti_rows(G), function(r) lapply(r, make))
}

# the characteristic polynomial of A in s over the domain of A
.ss_stability <- function(self, fast, schur) {
  s <- Symbol("s")
  cf <- .charpoly_coeffs(self$A)
  cp <- expand(.add(lapply(seq_along(cf), function(k) .mul(list(cf[[k]], .pow(s, .int(length(cf) - k)))))))
  kind <- if (isTRUE(fast)) "EX" else .rc_domain(self$A$data)
  conds <- if (schur) .lti_schur(cp, s, fast, kind = kind) else .lti_hurwitz(cp, s, fast, kind = kind)
  lapply(conds, function(c) c > 0)
}

# StateSpace.dsolve: the output for the initial state and a constant or
# time dependent input, from the matrix exponential of A t
.ss_dsolve <- function(self, initial_conditions = NULL, input_vector = NULL, var = Symbol("t")) {
  if (!(is_sym(var) && var$type == "Symbol")) .physics_error("ValueError", "Variable for representing time must be a Symbol.")
  n <- self$A$nrow
  if (is.null(initial_conditions) || all(vapply(initial_conditions$data, is_zero_num, TRUE))) {
    initial_conditions <- zeros(n, 1L)
  } else if (initial_conditions$nrow != n || initial_conditions$ncol != 1L) {
    .physics_error("ShapeError", "Initial condition vector should have the same number of rows as the state matrix.")
  }
  m <- self$B$ncol
  if (is.null(input_vector) || all(vapply(input_vector$data, is_zero_num, TRUE))) {
    input_vector <- zeros(m, 1L)
  } else if (input_vector$nrow != m || input_vector$ncol != 1L) {
    .physics_error("ShapeError", "Input vector should have the same number of columns as the input matrix.")
  }
  tau <- Dummy("tau")
  E <- matrix_exp(self$A * var)
  x <- E %*% initial_conditions
  Bu <- self$B %*% input_vector
  if (!all(vapply(Bu$data, is_zero_num, TRUE))) {
    integrand <- subs(E, var, var - tau) %*% subs(Bu, var, tau)
    x <- x + .mat_map(integrand, function(z) integrate(expand(z), list(tau, S_Zero, var)))
  }
  res <- self$C %*% x + self$D %*% input_vector
  .mat_map(res, function(z) expand(simplify(z)))
}

local({
  ssb_props <- list(
    A = function(self) .lti_mat(self$args[[1L]]), B = function(self) .lti_mat(self$args[[2L]]),
    C = function(self) .lti_mat(self$args[[3L]]), D = function(self) .lti_mat(self$args[[4L]]),
    state_matrix = function(self) .lti_mat(self$args[[1L]]), input_matrix = function(self) .lti_mat(self$args[[2L]]),
    output_matrix = function(self) .lti_mat(self$args[[3L]]), feedforward_matrix = function(self) .lti_mat(self$args[[4L]]),
    num_states = function(self) as.integer(.int_val(self$args[[1L]]$args[[1L]])),
    num_inputs = function(self) as.integer(.int_val(self$args[[4L]]$args[[2L]])),
    num_outputs = function(self) as.integer(.int_val(self$args[[4L]]$args[[1L]])),
    shape = function(self) .pytuple(list(self$num_outputs, self$num_inputs)),
    is_SISO = function(self) isTRUE(self$extra$is_SISO),
    is_StateSpace_object = function(self) TRUE)
  .pyclass_def("StateSpaceBase", parent = "LinearTimeInvariant",
    props = ssb_props,
    methods = list(
      evalf = function(self, n = 15L, ...)
        create_state_space(N(self$A, n), N(self$B, n), N(self$C, n), N(self$D, n), self$sampling_time),
      append = .ss_append,
      observability_matrix = .ss_observability,
      unobservable_subspace = .ss_unobservable,
      observable_subspace = .ss_observable,
      is_observable = function(self) rank_matrix(.ss_observability(self)) == self$num_states,
      controllability_matrix = .ss_controllability,
      uncontrollable_subspace = .ss_uncontrollable,
      controllable_subspace = .ss_controllable,
      is_controllable = function(self) rank_matrix(.ss_controllability(self)) == self$num_states,
      apply_similarity_transform = .ss_similarity,
      to_observable_form = function(self) .ss_form(self, .ss_observable, .ss_unobservable),
      to_controllable_form = function(self) .ss_form(self, .ss_controllable, .ss_uncontrollable)),
    ops = list(
      `__add__` = .ss_add,
      `__radd__` = .ss_add,
      `__sub__` = function(self, other) .ss_add(self, -other),
      `__rsub__` = function(self, other) .py_binop("+", other, -self),
      `__neg__` = .ss_neg,
      `__mul__` = .ss_mul,
      `__rmul__` = .ss_rmul))
  .pyclass_def("StateSpace", parent = "StateSpaceBase",
    ctor = function(A, B, C, D) StateSpace(A, B, C, D),
    props = list(sampling_time = function(self) S_Zero, is_continuous = function(self) TRUE),
    methods = list(
      dsolve = .ss_dsolve,
      get_asymptotic_stability_conditions = function(self, fast = FALSE) .ss_stability(self, fast, FALSE),
      `_eval_rewrite_as_TransferFunction` = function(self) {
        s <- Symbol("s")
        .ss_to_tf(self, s, function(e) TransferFunction$from_rational_expression(e, s))
      },
      `_eval_rewrite_as_DiscreteTransferFunction` = function(self)
        .physics_error("TypeError", "\n            The continuous state space model cannot be rewritten as a\n            discrete-time transfer function model.\n            ")))
  .pyclass_def("DiscreteStateSpace", parent = "StateSpaceBase",
    ctor = function(A, B, C, D, sampling_time) DiscreteStateSpace(A, B, C, D, sampling_time),
    props = list(sampling_time = function(self) self$args[[5L]], is_continuous = function(self) FALSE),
    methods = list(
      get_asymptotic_stability_conditions = function(self, fast = FALSE) .ss_stability(self, fast, TRUE),
      `_eval_rewrite_as_DiscreteTransferFunction` = function(self) {
        z <- Symbol("z")
        .ss_to_tf(self, z, function(e) DiscreteTransferFunction$from_rational_expression(e, z, self$sampling_time))
      },
      `_eval_rewrite_as_TransferFunction` = function(self)
        .physics_error("TypeError", "\n            The discrete state space model cannot be rewritten as a\n            continuous-time transfer function model.\n            ")))
})
