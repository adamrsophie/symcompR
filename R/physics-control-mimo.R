# MIMO systems, after sympy/physics/control/lti.py: the transfer function
# matrix and the series, parallel and feedback interconnections of MIMO
# systems.

.mimo_other <- function(f) function(self, other) {
  if (!inherits(other, "MIMOLinearTimeInvariant")) return(NULL)
  f(self, other)
}

# an ImmutableMatrix argument as an explicit matrix
.lti_mat <- function(x) {
  if (inherits(x, "symmatrix")) return(x)
  .new_mat(x$args[[3L]]$args, .int_val(x$args[[1L]]), .int_val(x$args[[2L]]))
}

.lti_rows <- function(M) lapply(seq_len(M$nrow), function(i) lapply(seq_len(M$ncol), function(j) .mget(M, i, j)))

# _to_TFM
.lti_to_tfm <- function(M, var, sampling_time = 0L) {
  to_tf <- if (.lti_st_zero(sampling_time)) function(e) TransferFunction$from_rational_expression(e, var)
           else function(e) DiscreteTransferFunction$from_rational_expression(e, var, sampling_time)
  TransferFunctionMatrix(lapply(.lti_rows(M), function(r) lapply(r, to_tf)))
}

# the transfer function matrix of a MIMO system
.lti_expr_mat <- function(x) x$doit()$`_expr_mat`

#' Transfer function matrices and MIMO interconnections
#'
#' `TransferFunctionMatrix(arg)` is the transfer function matrix of a MIMO
#' system, built from a list of rows of SISO systems.
#' `TransferFunctionMatrix$from_Matrix(matrix, var)` builds one from a
#' matrix of rational expressions. `MIMOSeries()`, `MIMOParallel()` and
#' `MIMOFeedback()` are the interconnections of MIMO systems. A transfer
#' function matrix has the properties `$var`, `$num_inputs`,
#' `$num_outputs` and `$shape` and the methods `$transpose()`,
#' `$elem_poles()`, `$elem_zeros()`, `$eval_frequency()`, `$simplify()` and
#' `$expand()`; `tfm[i, j]` is an element or a submatrix. After
#' `sympy.physics.control.lti`.
#'
#' @param arg A list of rows of SISO systems.
#' @param ... The systems.
#' @param evaluate Whether to evaluate the interconnection.
#' @param sys1,sys2 The systems of the feed forward and feedback paths.
#' @param sign The sign of the feedback.
#' @param x A transfer function matrix.
#' @param i,j The row and column indices.
#' @param drop Whether a single element is returned as a SISO system.
#' @return A `TransferFunctionMatrix`, `MIMOSeries`, `MIMOParallel` or
#'   `MIMOFeedback` object.
#' @examples
#' s <- Symbol("s")
#' G1 <- TransferFunction(1, s + 1, s); G2 <- TransferFunction(s, s + 2, s)
#' M <- TransferFunctionMatrix(list(list(G1, G2), list(G2, G1)))
#' M$shape
#' MIMOSeries(M, M)$doit()
#' @usage TransferFunctionMatrix(arg)
#' @export
TransferFunctionMatrix <- .pyclass_fn(function(arg) {
  if (is_sym(arg) && arg$type == "Tuple") arg <- lapply(arg$args, function(r) r$args)
  bad <- .py_filldedent("\n                `arg` param in TransferFunctionMatrix should\n                strictly be a nested list containing TransferFunctionBase\n                objects.")
  if (!.lti_is_list(arg) || !length(arg) || !.lti_is_list(arg[[1L]])) .physics_error("ValueError", bad)
  if (!length(arg[[1L]])) .physics_error("IndexError", "list index out of range")
  first <- arg[[1L]][[1L]]
  if (!inherits(first, "LinearTimeInvariant"))
    .physics_error("AttributeError", sprintf("'%s' object has no attribute 'var'", .py_typename(first)))
  var <- first$var
  rows <- lapply(arg, function(row) {
    lapply(row, function(el) {
      if (!inherits(el, "SISOLinearTimeInvariant"))
        .physics_error("TypeError", .py_filldedent("\n                        Each element is expected to be of\n                        type `SISOLinearTimeInvariant`."))
      if (!.lti_eq(var, el$var))
        .physics_error("ValueError", .py_filldedent("\n                        Conflicting value(s) found for `var`.\n                        All TransferFunction instances in TransferFunctionMatrix\n                        should use the same complex variable in Laplace domain\n                        or z-domain."))
      el$to_expr()
    })
  })
  .lti_check_time(unlist(arg, recursive = FALSE))
  if (length(unique(lengths(rows))) > 1L) .physics_error("ValueError", "mismatched dimensions")
  M <- Matrix(rows)
  targ <- do.call(Tuple, lapply(arg, function(r) do.call(Tuple, r)))
  .pyobj_new("TransferFunctionMatrix", list(targ),
             extra = list(expr_mat = M, is_StateSpace_object = FALSE, is_continuous = first$is_continuous))
}, "TransferFunctionMatrix")

.lti_is_list <- function(x) (is.list(x) && !is.object(x)) || inherits(x, c("pytuple", "pylist"))

#' @rdname TransferFunctionMatrix
#' @export
`[.TransferFunctionMatrix` <- function(x, i, j, drop = TRUE) {
  M <- x$extra$expr_mat
  if (missing(i)) i <- seq_len(M$nrow)
  if (missing(j)) j <- seq_len(M$ncol)
  if (isTRUE(drop) && length(i) == 1L && length(j) == 1L) {
    e <- .mget(M, i, j)
    if (.lti_st_zero(x$sampling_time)) return(TransferFunction$from_rational_expression(e, x$var))
    return(DiscreteTransferFunction$from_rational_expression(e, x$var, x$sampling_time))
  }
  sub <- Matrix(length(i), length(j), function(a, b) .mget(M, i[[a]], j[[b]]))
  .lti_to_tfm(sub, x$var, x$sampling_time)
}

local({
  .pyclass_def("TransferFunctionMatrix", parent = "MIMOLinearTimeInvariant",
    ctor = function(arg) TransferFunctionMatrix(arg),
    props = list(
      var = function(self) self$args[[1L]]$args[[1L]]$args[[1L]]$var,
      num_inputs = function(self) self$extra$expr_mat$ncol,
      num_outputs = function(self) self$extra$expr_mat$nrow,
      shape = function(self) .pytuple(list(self$extra$expr_mat$nrow, self$extra$expr_mat$ncol)),
      `_expr_mat` = function(self) self$extra$expr_mat,
      sampling_time = function(self) self$args[[1L]]$args[[1L]]$args[[1L]]$sampling_time),
    methods = list(
      transpose = function(self) .lti_to_tfm(t(self$extra$expr_mat), self$var, self$sampling_time),
      elem_poles = function(self) lapply(self$doit()$args[[1L]]$args, function(r) lapply(r$args, function(el) el$poles())),
      elem_zeros = function(self) lapply(self$doit()$args[[1L]]$args, function(r) lapply(r$args, function(el) el$zeros())),
      eval_frequency = function(self, other) expand(subs(self$extra$expr_mat, self$var, other)),
      `_flat` = function(self) unlist(lapply(self$args[[1L]]$args, function(r) r$args), recursive = FALSE),
      evalf = function(self, n = 15L, ...) .lti_to_tfm(N(self$extra$expr_mat, n), self$var, self$sampling_time),
      n = function(self, n = 15L, ...) .lti_to_tfm(N(self$extra$expr_mat, n), self$var, self$sampling_time),
      simplify = function(self, ...) .lti_to_tfm(.mat_map(self$extra$expr_mat, .lti_cancel_noexpand), self$var, self$sampling_time),
      expand = function(self, ...) .lti_to_tfm(expand(self$extra$expr_mat), self$var, self$sampling_time),
      `__getitem__` = function(self, i, j, drop = TRUE) self[i, j, drop = drop]),
    statics = list(from_Matrix = function(matrix, var, sampling_time = 0L) .lti_to_tfm(Matrix(matrix), var, sampling_time)),
    ops = list(
      `__neg__` = function(self) .lti_to_tfm(-self$extra$expr_mat, self$var, self$sampling_time),
      `__add__` = .mimo_other(function(self, other) {
        if (!inherits(other, "MIMOParallel")) return(MIMOParallel(self, other))
        do.call(MIMOParallel, c(list(self), other$args))
      }),
      `__sub__` = .mimo_other(function(self, other) self + (-other)),
      `__mul__` = .mimo_other(function(self, other) {
        if (!inherits(other, "MIMOSeries")) return(MIMOSeries(other, self))
        do.call(MIMOSeries, c(other$args, list(self)))
      })))
})

#' @rdname TransferFunctionMatrix
#' @export
MIMOSeries <- function(..., evaluate = FALSE) {
  args <- .lti_flatten(list(...), "MIMOSeries")
  ss <- FALSE
  if (length(args) && .lti_any_ss(args)) {
    if (!.lti_io_compatible(args))
      .physics_error("ValueError", .py_filldedent("\n                    Systems with incompatible inputs and outputs\n                    cannot be connected in MIMOSeries."))
    ss <- TRUE
  } else {
    .lti_check_args(args, "MIMOLinearTimeInvariant")
    if (!.lti_io_compatible(args))
      .physics_error("ValueError", .py_filldedent("\n                    Number of input signals do not match the number\n                    of output signals of adjacent systems for some args."))
  }
  .lti_check_time(args)
  obj <- .pyobj_new("MIMOSeries", args, extra = list(ss = ss, is_continuous = args[[1L]]$is_continuous))
  if (isTRUE(evaluate)) obj$doit() else obj
}

.mimo_series_doit <- function(self, cancel = FALSE, ...) {
  if (isTRUE(self$extra$ss)) {
    ss_class <- .lti_ss_class(self$is_continuous)
    res <- self$args[[1L]]
    if (!inherits(res, ss_class)) res <- res$doit()$rewrite(ss_class)
    for (arg in self$args[-1L]) {
      if (!inherits(arg, ss_class)) arg <- arg$doit()$rewrite(ss_class)
      res <- arg * res
    }
    return(res)
  }
  mats <- lapply(rev(self$args), .lti_expr_mat)
  if (isTRUE(cancel)) return(.lti_to_tfm(Reduce(`%*%`, mats), self$var))
  dums <- lapply(mats, function(m) Dummy())
  mats <- Map(function(m, d) subs(m, self$var, d), mats, dums)
  res <- .lti_to_tfm(Reduce(`%*%`, mats), self$var, self$sampling_time)
  for (d in dums) res <- subs(res, d, self$var)
  res
}

.mimo_rewrite_tfm <- function(self) {
  if (isTRUE(self$is_StateSpace_object)) return(self$doit()$rewrite(.lti_tf_class(self$is_continuous)))
  self$doit()
}

#' @rdname TransferFunctionMatrix
#' @export
MIMOParallel <- function(..., evaluate = FALSE) {
  args <- .lti_flatten(list(...), "MIMOParallel")
  ss <- FALSE
  if (length(args) && .lti_any_ss(args)) {
    nin <- vapply(args, function(a) as.integer(a$num_inputs), 0L)
    nout <- vapply(args, function(a) as.integer(a$num_outputs), 0L)
    if (!.lti_io_compatible(args))
      .physics_error("ShapeError", "Systems with incompatible inputs and outputscannot be connected in MIMOParallel.")
    ss <- TRUE
  } else {
    .lti_check_args(args, "MIMOLinearTimeInvariant")
    sh <- lapply(args, function(a) unlist(a$shape))
    if (any(vapply(sh, function(x) !identical(x, sh[[1L]]), TRUE)))
      .physics_error("TypeError", "Shape of all the args is not equal.")
  }
  .lti_check_time(args)
  obj <- .pyobj_new("MIMOParallel", args, extra = list(ss = ss, is_continuous = args[[1L]]$is_continuous))
  if (isTRUE(evaluate)) obj$doit() else obj
}

.mimo_parallel_doit <- function(self, ...) {
  if (isTRUE(self$extra$ss)) {
    res <- self$args[[1L]]
    if (!inherits(res, "StateSpace")) res <- res$doit()$rewrite("StateSpace")
    for (arg in self$args[-1L]) {
      arg <- if (!inherits(arg, "StateSpace")) arg$doit()$rewrite("StateSpace") else arg$doit()
      res <- res + arg
    }
    return(res)
  }
  mats <- lapply(self$args, .lti_expr_mat)
  .lti_to_tfm(Reduce(`+`, mats), self$var, self$sampling_time)
}

#' @rdname TransferFunctionMatrix
#' @export
MIMOFeedback <- function(sys1, sys2, sign = -1L) {
  ok <- c("TransferFunctionMatrix", "MIMOSeries", "StateSpaceBase")
  if (!inherits(sys1, ok)) .physics_error("TypeError", "Unsupported type for `sys1` in MIMO Feedback.")
  if (!inherits(sys2, ok)) .physics_error("TypeError", "Unsupported type for `sys2` in MIMO Feedback.")
  if (!identical(as.integer(sys1$num_inputs), as.integer(sys2$num_outputs)) ||
      !identical(as.integer(sys1$num_outputs), as.integer(sys2$num_inputs)))
    .physics_error("ValueError", .py_filldedent("\n                Product of `sys1` and `sys2` must\n                yield a square matrix."))
  sgn <- if (is_sym(sign)) suppressWarnings(as.numeric(N(sign))) else sign
  if (!(is.numeric(sgn) && length(sgn) == 1L && !is.na(sgn) && sgn %in% c(-1, 1)))
    .physics_error("ValueError", .py_filldedent("\n                Unsupported type for feedback. `sign` arg should\n                either be 1 (positive feedback loop) or -1\n                (negative feedback loop)."))
  ss <- isTRUE(sys1$is_StateSpace_object) || isTRUE(sys2$is_StateSpace_object)
  if (!ss) {
    M <- eye(as.integer(sys1$num_outputs)) - (.lti_expr_mat(sys1) %*% .lti_expr_mat(sys2)) * as_sym(sign)
    # SymPy simplifies the 2 by 2 determinant (dotprodsimp) before comparing it with zero
    if (.mimo_det_zero(det(M))) .physics_error("ValueError", "Non-Invertible system inputted.")
    if (!.lti_eq(sys1$var, sys2$var))
      .physics_error("ValueError", .py_filldedent("\n                Both `sys1` and `sys2` should be using the\n                same complex variable."))
  }
  .lti_check_time(list(sys1, sys2))
  .pyobj_new("MIMOFeedback", list(sys1, sys2, as_sym(sign)),
             extra = list(is_StateSpace_object = ss, is_continuous = sys1$is_continuous))
}

.mimo_sensitivity <- function(self)
  inverse(eye(as.integer(self$sys1$num_inputs)) - (.lti_expr_mat(self$sys1) %*% .lti_expr_mat(self$sys2)) * self$sign)

.mimo_feedback_doit <- function(self, cancel = TRUE, expand = FALSE, ...) {
  if (isTRUE(self$is_StateSpace_object)) return(.feedback_ss_doit(self))
  M <- .mimo_sensitivity(self) %*% .lti_expr_mat(self$sys1)
  res <- .lti_to_tfm(M, self$var, self$sampling_time)
  if (isTRUE(cancel)) res <- res$simplify()
  if (isTRUE(expand)) res <- res$expand()
  res
}

local({
  rw <- list(`_eval_rewrite_as_TransferFunctionMatrix` = .mimo_rewrite_tfm)
  shape <- function(self) .pytuple(list(self$num_outputs, self$num_inputs))
  .pyclass_def("MIMOSeries", parent = "MIMOLinearTimeInvariant",
    ctor = function(...) MIMOSeries(...),
    props = list(var = function(self) self$args[[1L]]$var,
                 num_inputs = function(self) self$args[[1L]]$num_inputs,
                 num_outputs = function(self) self$args[[length(self$args)]]$num_outputs,
                 shape = shape,
                 is_StateSpace_object = function(self) isTRUE(self$extra$ss),
                 sampling_time = function(self) self$args[[1L]]$sampling_time),
    methods = c(rw, list(doit = .mimo_series_doit)),
    ops = list(
      `__add__` = .mimo_other(function(self, other) {
        if (inherits(other, "MIMOParallel")) return(do.call(MIMOParallel, c(list(self), other$args)))
        MIMOParallel(self, other)
      }),
      `__radd__` = .mimo_other(function(self, other) {
        if (inherits(other, "MIMOParallel")) return(do.call(MIMOParallel, c(list(self), other$args)))
        MIMOParallel(self, other)
      }),
      `__sub__` = .mimo_other(function(self, other) self + (-other)),
      `__rsub__` = function(self, other) .py_binop("+", -self, other),
      `__mul__` = .mimo_other(function(self, other) {
        if (inherits(other, "MIMOSeries")) return(do.call(MIMOSeries, c(other$args, self$args)))
        do.call(MIMOSeries, c(list(other), self$args))
      }),
      `__neg__` = function(self) {
        a <- self$args
        a[[1L]] <- -a[[1L]]
        do.call(MIMOSeries, a)
      }))
  .pyclass_def("MIMOParallel", parent = "MIMOLinearTimeInvariant",
    ctor = function(...) MIMOParallel(...),
    props = list(var = function(self) self$args[[1L]]$var,
                 num_inputs = function(self) self$args[[1L]]$num_inputs,
                 num_outputs = function(self) self$args[[1L]]$num_outputs,
                 shape = shape,
                 is_StateSpace_object = function(self) isTRUE(self$extra$ss),
                 sampling_time = function(self) self$args[[1L]]$sampling_time),
    methods = c(rw, list(doit = .mimo_parallel_doit)),
    ops = list(
      `__add__` = .mimo_other(function(self, other) do.call(MIMOParallel, c(self$args, list(other)))),
      `__radd__` = .mimo_other(function(self, other) do.call(MIMOParallel, c(self$args, list(other)))),
      `__sub__` = .mimo_other(function(self, other) self + (-other)),
      `__rsub__` = function(self, other) .py_binop("+", -self, other),
      `__mul__` = .mimo_other(function(self, other) {
        if (inherits(other, "MIMOSeries")) return(do.call(MIMOSeries, c(other$args, list(self))))
        MIMOSeries(other, self)
      }),
      `__neg__` = function(self) do.call(MIMOParallel, lapply(self$args, function(a) -a))))
  .pyclass_def("MIMOFeedback", parent = "MIMOLinearTimeInvariant",
    ctor = function(sys1, sys2, sign) MIMOFeedback(sys1, sys2, sign),
    props = list(sys1 = function(self) self$args[[1L]], sys2 = function(self) self$args[[2L]],
                 var = function(self) self$args[[1L]]$var, sign = function(self) self$args[[3L]],
                 sensitivity = .mimo_sensitivity,
                 num_inputs = function(self) self$args[[1L]]$num_inputs,
                 num_outputs = function(self) self$args[[1L]]$num_outputs,
                 sampling_time = function(self) self$args[[1L]]$sampling_time),
    methods = list(doit = .mimo_feedback_doit,
                   `_eval_rewrite_as_TransferFunctionMatrix` = function(self) self$doit()),
    ops = list(`__neg__` = function(self) MIMOFeedback(-self$sys1, -self$sys2, self$sign)))
  # Basic.doit of an expression holding these objects evaluates them
  for (cls in c("Series", "Parallel", "Feedback", "MIMOSeries", "MIMOParallel", "MIMOFeedback"))
    register_doit(cls, function(e) e$doit())
})

# the Python class names, for isinstance() style checks
attr(MIMOSeries, "pyname") <- "MIMOSeries"
attr(MIMOParallel, "pyname") <- "MIMOParallel"
attr(MIMOFeedback, "pyname") <- "MIMOFeedback"

# whether a determinant is zero: a nonzero value at a sample point settles it
# without the slower cancel()
.mimo_det_zero <- function(d) {
  if (is_zero_num(d)) return(TRUE)
  nm <- vapply(free_symbols(d), function(z) z$name, "")
  for (pt in c(7 / 3, 11 / 5)) {
    v <- tryCatch(.neval(d, stats::setNames(as.list(rep(pt, length(nm))), nm)), error = function(e) NA)
    if (length(v) == 1L && is.finite(v) && base::Mod(v) > 1e-9) return(FALSE)
  }
  is_zero_num(cancel(d))
}
