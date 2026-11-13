# Applying operators to states, matrix representations and the maps
# between operators and their eigenstates, after
# sympy/physics/quantum/qapply.py, represent.py and operatorset.py.

# Basic.replace(cls, f) for a class: bottom up, f applied to every instance
.qreplace_type <- function(e, cls, f) {
  if (!is_sym(e)) return(e)
  if (length(e$args)) {
    args <- lapply(e$args, .qreplace_type, cls = cls, f = f)
    if (!identical(vapply(args, .k, ""), vapply(e$args, .k, ""))) e <- sym_rebuild(e, args)
  }
  if (.qisa(e, cls)) e <- f(e)
  e
}

.ip_doit_func <- function(e) .qreplace_type(e, "InnerProduct", function(z) doit(z))
.sum_doit_func <- function(e) .qreplace_type(e, "Sum", function(z) doit(z))

#' Apply operators to states
#'
#' `qapply(e)` applies the operators in `e` to the states they act on, as
#' SymPy's `sympy.physics.quantum.qapply`: products are evaluated from the
#' right, commutators and tensor products are expanded, and inner products
#' are evaluated unless `ip_doit = FALSE`. With `dagger = TRUE` operators are
#' also applied to bras on their left; `sum_doit = TRUE` evaluates sums.
#'
#' @param e An expression of operators and states.
#' @param dagger Apply operators to the left too.
#' @param sum_doit Evaluate `Sum` objects.
#' @param ip_doit Evaluate inner products.
#' @param ... Further options for the operators.
#' @return An expression.
#' @examples
#' qapply(OuterProduct(Ket("a"), Bra("b")) * Ket("c"))
#' @export
qapply <- function(e, dagger = FALSE, sum_doit = FALSE, ip_doit = TRUE, ...) {
  opts <- c(list(dagger = dagger, sum_doit = sum_doit, ip_doit = ip_doit), list(...))
  .qapply(.qS(e), opts)
}

.qapply_opts_call <- function(e, opts) do.call(qapply, c(list(e), opts))

.qapply <- function(e, opts) {
  if (inherits(e, "symmatrix")) return(e)
  if (.qkind(e) == "NumberKind") return(if (isTRUE(opts$ip_doit)) .ip_doit_func(e) else e)
  e <- expand(e, commutator = TRUE, tensorproduct = TRUE)
  if (.qisa(e, "KetBase")) return(e)
  if (e$type == "Add") {
    result <- S_Zero
    for (arg in .sympy_args(e)) result <- .qplus(result, .qapply(arg, opts))
    return(expand(result))
  }
  if (.qisa(e, "QDensity")) {
    new_args <- lapply(e$args, function(sp) Tuple(.qapply(sp$args[[1L]], opts), sp$args[[2L]]))
    return(.qfunc("QDensity", new_args))
  }
  if (e$type == "QTensorProduct") return(.qtp(lapply(e$args, function(t) .qapply(t, opts))))
  if (e$type == "Sum") {
    result <- new_sym("Sum", c(list(.qapply(e$args[[1L]], opts)), e$args[-1L]))
    return(if (isTRUE(opts$sum_doit)) .sum_doit_func(result) else result)
  }
  if (e$type == "Pow") return(.pow(.qapply(e$args[[1L]], opts), e$args[[2L]]))
  if (e$type == "Mul") {
    cn <- .args_cnc(e)
    c_mul <- .mul(cn$c)
    nc_mul <- .mul(cn$nc)
    result <- if (!length(cn$nc)) c_mul
              else if (nc_mul$type == "Mul") .qtimes(c_mul, .qapply_mul(nc_mul, opts))
              else .qtimes(c_mul, .qapply(nc_mul, opts))
    if (.qeq(result, e) && isTRUE(opts$dagger)) result <- Dagger(.qapply_mul(Dagger(e), opts))
    if (isTRUE(opts$ip_doit)) result <- .ip_doit_func(result)
    if (isTRUE(opts$sum_doit)) result <- .sum_doit_func(result)
    return(result)
  }
  e
}

.qop_like <- function(a) .qisa(a, c("Operator", "State", "Mul", "Pow")) || .qeq(a, S_One)

.qapply_mul <- function(e, opts) {
  args <- .qargs(e)
  extra <- S_One
  result <- NULL
  if (length(args) <= 1L || e$type != "Mul") return(e)
  n <- length(args)
  rhs <- args[[n]]; lhs <- args[[n - 1L]]
  args <- args[seq_len(n - 2L)]
  if ((!.qisa(rhs, "Wavefunction") && .is_comm(rhs)) || (!.qisa(lhs, "Wavefunction") && .is_comm(lhs))) return(e)
  if (lhs$type == "Pow" && lhs$args[[2L]]$type == "Integer") {
    args <- c(args, list(.pow(lhs$args[[1L]], .sub(lhs$args[[2L]], S_One))))
    lhs <- lhs$args[[1L]]
  }
  if (.qisa(lhs, "OuterProduct")) {
    args <- c(args, list(lhs$args[[1L]]))
    lhs <- lhs$args[[2L]]
  }
  if (.qisa(rhs, "OuterProduct")) {
    extra <- rhs$args[[2L]]
    rhs <- rhs$args[[1L]]
  }
  if (.qisa(lhs, c("QCommutator", "AntiCommutator"))) {
    comm <- doit(lhs)
    if (comm$type == "Add") {
      ca <- .sympy_args(comm)
      return(.qtimes(.qapply(.add(list(.mul(c(args, list(ca[[1L]], rhs))), .mul(c(args, list(ca[[2L]], rhs))))), opts), extra))
    }
    return(.qtimes(.qapply(.mul(list(.mul(args), comm, rhs)), opts), extra))
  }
  if (lhs$type == "QTensorProduct" && all(vapply(lhs$args, .qop_like, NA)) &&
      rhs$type == "QTensorProduct" && all(vapply(rhs$args, .qop_like, NA)) &&
      length(lhs$args) == length(rhs$args)) {
    result <- expand(.qtp(lapply(seq_along(lhs$args), function(k) .qapply(.mul(list(lhs$args[[k]], rhs$args[[k]])), opts))),
                     tensorproduct = TRUE)
    return(.mul(list(.qapply_mul(.mul(args), opts), result, extra)))
  }
  if (rhs$type == "Sum") {
    if (lhs$type == "Sum") {
      lv <- vapply(lhs$args[-1L], function(l) l$args[[1L]]$key, "")
      rv <- vapply(rhs$args[-1L], function(l) l$args[[1L]]$key, "")
      if (length(intersect(lv, rv))) .qstop("ValueError", "Duplicated dummy indices in separate sums in qapply.")
      result <- new_sym("Sum", c(list(.qapply(.mul(list(lhs$args[[1L]], rhs$args[[1L]])), opts)), lhs$args[-1L], rhs$args[-1L]))
      return(.qapply_mul(.mul(list(.mul(args), result)), opts))
    }
    result <- new_sym("Sum", c(list(.qapply(.mul(list(lhs, rhs$args[[1L]])), opts)), rhs$args[-1L]))
    return(.qapply_mul(.mul(list(.mul(args), result)), opts))
  }
  if (lhs$type == "Sum") {
    result <- new_sym("Sum", c(list(.qapply(.mul(list(lhs$args[[1L]], rhs)), opts)), lhs$args[-1L]))
    return(.qapply_mul(.mul(list(.mul(args), result)), opts))
  }
  ex <- opts[setdiff(names(opts), c("dagger", "sum_doit", "ip_doit"))]
  if (!is.null(.qnode_def(lhs$type)) && .qhas(lhs, "_apply_operator"))
    result <- tryCatch(do.call(.qmeth, c(list(lhs, "_apply_operator", rhs), opts)), NotImplementedError = function(err) NULL)
  if (is.null(result) && !is.null(.qnode_def(rhs$type)) && .qhas(rhs, "_apply_from_right_to"))
    result <- tryCatch(do.call(.qmeth, c(list(rhs, "_apply_from_right_to", lhs), opts)), NotImplementedError = function(err) NULL)
  if (is.null(result) && .qisa(lhs, "BraBase") && .qisa(rhs, "KetBase")) result <- InnerProduct(lhs, rhs)
  if (is.numeric(result)) return(as_sym(result))
  if (is.null(result)) {
    if (!length(args)) return(e)
    return(.qtimes(.qtimes(.qapply_mul(.mul(c(args, list(lhs))), opts), rhs), extra))
  }
  if (.qisa(result, "InnerProduct")) return(.qtimes(.qtimes(result, .qapply_mul(.mul(args), opts)), extra))
  .qtimes(.qapply(.qtimes(.mul(args), result), opts), extra)
}

# ---------------------------------------------------------------------------
# represent

#' Matrix representations
#'
#' `represent(expr, basis = ...)` gives the matrix (or wavefunction) of a
#' quantum expression in a basis, as SymPy's
#' `sympy.physics.quantum.represent`. The basis is an operator, a state or a
#' class name such as `"JzOp"`; without one the default basis of the object
#' is used. `rep_innerproduct()` and `rep_expectation()` build the position
#' representation of states and operators, `get_basis()` finds the basis
#' state, `enumerate_states()` lists states with consecutive labels and
#' `integrate_result()` integrates over the unities inserted by
#' `represent()`.
#'
#' @param expr A quantum expression.
#' @param basis A basis operator, state or class name.
#' @param ... Options such as `nqubits`, `j` or `format`.
#' @param orig_expr,result The represented expression and its matrix.
#' @return A matrix, an expression or a list of states.
#' @examples
#' represent(IdentityOperator(2))
#' @export
represent <- function(expr, ...) {
  opts <- list(...)
  if (is.null(opts$format)) opts$format <- "sympy"
  .represent(.qS(expr), opts)
}

.qcall_opts <- function(obj, nm, opts, ...) do.call(.qmeth, c(list(obj, nm, ...), opts))

.represent <- function(expr, opts) {
  if (inherits(expr, "symmatrix")) return(expr)
  if (.qisa(expr, "QExpr") && !.qisa(expr, "OuterProduct")) {
    opts$replace_none <- FALSE
    temp_basis <- do.call(get_basis, c(list(expr), opts))
    if (!is.null(temp_basis)) opts$basis <- temp_basis
    r <- tryCatch(list(ok = .qcall_opts(expr, "_represent", opts)), NotImplementedError = function(err) list(err = err))
    if (!is.null(r$err)) {
      opts$replace_none <- TRUE
      if (.qisa(expr, c("KetBase", "BraBase"))) {
        return(tryCatch(do.call(rep_innerproduct, c(list(expr), opts)),
                        NotImplementedError = function(e2) stop(r$err)))
      } else if (.qisa(expr, "Operator")) {
        return(tryCatch(do.call(rep_expectation, c(list(expr), opts)),
                        NotImplementedError = function(e2) stop(r$err)))
      }
      stop(r$err)
    }
    return(r$ok)
  }
  if (expr$type == "Add") {
    a <- .sympy_args(expr)
    result <- .represent(a[[1L]], opts)
    for (x in a[-1L]) result <- .qplus(result, .represent(x, opts))
    return(result)
  }
  if (expr$type == "Pow") {
    base <- .represent(expr$args[[1L]], opts)
    ex <- expr$args[[2L]]
    if (inherits(base, "symmatrix")) return(base %^% ex)
    return(.pow(base, ex))
  }
  if (expr$type == "QTensorProduct") return(.qtp(lapply(expr$args, function(a) .represent(a, opts))))
  if (expr$type == "Function" && expr$name == "Dagger") return(Dagger(.represent(expr$args[[1L]], opts)))
  if (expr$type == "QCommutator") {
    A <- expr$args[[1L]]; B <- expr$args[[2L]]
    return(.represent(.sub(.mul(list(A, B)), .mul(list(B, A))), opts))
  }
  if (expr$type == "AntiCommutator") {
    A <- expr$args[[1L]]; B <- expr$args[[2L]]
    return(.represent(.add(list(.mul(list(A, B)), .mul(list(B, A)))), opts))
  }
  if (!(expr$type == "Mul" || .qisa(expr, c("OuterProduct", "InnerProduct")))) return(expr)
  opts$index <- if (!is.null(opts$index)) opts$index + 1L else 1L
  if (is.null(opts$unities)) { opts$unities <- new.env(); opts$unities$v <- integer() }
  args <- if (expr$type == "Mul") .qargs(expr) else expr$args
  n <- length(args)
  result <- .represent(args[[n]], opts)
  last_arg <- args[[n]]
  for (arg in rev(args[-n])) {
    if (.qisa(last_arg, "Operator")) {
      opts$index <- opts$index + 1L
      opts$unities$v <- c(opts$unities$v, opts$index)
    } else if (.qisa(last_arg, "BraBase") && .qisa(arg, "KetBase")) {
      opts$index <- opts$index + 1L
    } else if (.qisa(last_arg, "KetBase") && .qisa(arg, "Operator")) {
      opts$unities$v <- c(opts$unities$v, opts$index)
    } else if (.qisa(last_arg, "KetBase") && .qisa(arg, "BraBase")) {
      opts$unities$v <- c(opts$unities$v, opts$index)
    }
    next_arg <- .represent(arg, opts)
    result <- .qtimes(next_arg, result)
    last_arg <- arg
  }
  result <- flatten_scalar(result)
  do.call(integrate_result, c(list(expr, result), opts))
}

# flatten_scalar: a 1x1 matrix is its entry
flatten_scalar <- function(e) {
  if (inherits(e, "symmatrix") && e$nrow == 1L && e$ncol == 1L) return(e$data[[1L]])
  e
}

#' @rdname represent
#' @export
rep_innerproduct <- function(expr, ...) {
  opts <- list(...)
  if (!.qisa(expr, c("KetBase", "BraBase"))) .qstop("TypeError", "expr passed is not a Bra or Ket")
  basis <- do.call(get_basis, c(list(expr), opts))
  if (!.qisa(basis, "StateBase")) .qni("Can't form this representation!")
  if (is.null(opts$index)) opts$index <- 1L
  basis_kets <- enumerate_states(basis, opts$index, 2L)
  if (.qisa(expr, "BraBase")) {
    bra <- expr
    ket <- if (.qeq(.qa(basis_kets[[1L]], "dual"), expr)) basis_kets[[2L]] else basis_kets[[1L]]
  } else {
    bra <- if (.qeq(basis_kets[[1L]], expr)) .qa(basis_kets[[2L]], "dual") else .qa(basis_kets[[1L]], "dual")
    ket <- expr
  }
  doit(InnerProduct(bra, ket))
}

#' @rdname represent
#' @export
rep_expectation <- function(expr, ...) {
  opts <- list(...)
  if (is.null(opts$index)) opts$index <- 1L
  if (!.qisa(expr, "Operator")) .qstop("TypeError", "The passed expression is not an operator")
  basis_state <- do.call(get_basis, c(list(expr), opts))
  if (is.null(basis_state) || !.qisa(basis_state, "StateBase")) .qni("Could not get basis kets for this operator")
  basis_kets <- enumerate_states(basis_state, opts$index, 2L)
  bra <- .qa(basis_kets[[2L]], "dual")
  ket <- basis_kets[[1L]]
  qapply(.mul(list(bra, expr, ket)))
}

#' @rdname represent
#' @export
integrate_result <- function(orig_expr, result, ...) {
  opts <- list(...)
  if (!is_sym(result)) return(result)
  opts$replace_none <- TRUE
  if (is.null(opts$basis)) {
    arg <- orig_expr$args[[length(orig_expr$args)]]
    if (orig_expr$type == "Mul") { qa <- .qargs(orig_expr); arg <- qa[[length(qa)]] }
    opts$basis <- do.call(get_basis, c(list(arg), opts))
  } else if (!.qisa(opts$basis, "StateBase")) {
    opts$basis <- do.call(get_basis, c(list(orig_expr), opts))
  }
  basis <- opts$basis
  if (is.null(basis)) return(result)
  unities <- if (is.environment(opts$unities)) opts$unities$v else unlist(opts$unities)
  if (!length(unities)) return(result)
  kets <- enumerate_states(basis, as.list(unities))
  coords <- lapply(kets, function(k) .qa(k, "label")[[1L]])
  for (coord in coords) {
    if (!is.null(.free_sym_list(result)[[coord$key]])) {
      basis_op <- state_to_operators(basis)
      iv <- .qa(.qhs(basis_op), "interval")
      result <- integrate(result, list(coord, iv$args[[1L]], iv$args[[2L]]))
    }
  }
  result
}

.qclass_of <- function(x) {
  if (is.character(x)) return(x)
  if (is.function(x) && !is.null(attr(x, "qclass"))) return(attr(x, "qclass"))
  NULL
}

#' @rdname represent
#' @param replace_none Use the default basis when `basis` is `NULL`.
#' @export
get_basis <- function(expr, basis = NULL, replace_none = TRUE, ...) {
  if (is.null(basis) && !replace_none) return(NULL)
  if (is.null(basis)) {
    if (.qisa(expr, "KetBase")) return(.qmake_default(expr$type))
    if (.qisa(expr, "BraBase")) return(.qmake_default(.qmeth(expr, "dual_class")))
    if (.qisa(expr, "Operator")) return(operators_to_state(expr))
    return(NULL)
  }
  bc <- .qclass_of(basis)
  if (.qisa(basis, "Operator") || (!is.null(bc) && .qsubclass(bc, "Operator"))) {
    state <- operators_to_state(if (is.null(bc)) basis else bc)
    if (is.null(state)) return(NULL)
    if (.qisa(state, "StateBase")) return(state)
    return(.qmake_default(state))
  }
  if (.qisa(basis, "StateBase")) return(basis)
  if (!is.null(bc) && .qsubclass(bc, "StateBase")) return(.qmake_default(bc))
  NULL
}

.qmake_default <- function(cls) {
  if (!is.character(cls)) return(cls)
  tryCatch(.qfunc(cls, list()), error = function(e) cls)
}

#' @rdname represent
#' @export
enumerate_states <- function(...) {
  args <- list(...)
  state <- args[[1L]]
  if (!(length(args) %in% c(2L, 3L))) .qni("Wrong number of arguments!")
  if (!.qisa(state, "StateBase")) .qstop("TypeError", "First argument is not a state!")
  if (length(args) == 3L) {
    tryCatch(.qmeth(state, "_enumerate_state", .qint(args[[3L]]), start_index = .qint(args[[2L]])),
             NotImplementedError = function(e) list())
  } else {
    il <- lapply(args[[2L]], .qint)
    tryCatch(.qmeth(state, "_enumerate_state", length(il), index_list = il), NotImplementedError = function(e) list())
  }
}

# ---------------------------------------------------------------------------
# operatorset.py: states and the operators they are eigenstates of. Each
# entry is list(state class, operator classes).

.qstate_mapping <- function() list(
  list("JxKet", c("J2Op", "JxOp")), list("JyKet", c("J2Op", "JyOp")), list("JzKet", c("J2Op", "JzOp")),
  list("Ket", "Operator"), list("PositionKet3D", c("XOp", "YOp", "ZOp")), list("PxKet", "PxOp"),
  list("XKet", "XOp"))

.qset_key <- function(cls) paste(sort(cls, method = "radix"), collapse = "|")

#' Operators and their eigenstates
#'
#' `operators_to_state(ops)` gives the eigenstate of an operator or of a set
#' (an R list) of commuting operators, and `state_to_operators(state)` the
#' operators of a state, as SymPy's `sympy.physics.quantum.operatorset`.
#' Classes are given by name, for example `"XOp"` or `list("J2Op", "JzOp")`.
#'
#' @param operators An operator, a class name or a list of either.
#' @param state A state or a state class name.
#' @param ... Options for the state or operator constructors.
#' @return A state or operator, a list of operators, a class name or `NULL`.
#' @examples
#' operators_to_state(XOp())
#' state_to_operators(XKet())
#' @export
operators_to_state <- function(operators, ...) {
  sm <- .qstate_mapping()
  is_set <- is.list(operators) && !is_sym(operators)
  oc <- if (is.character(operators) && length(operators) > 1L) { is_set <- TRUE; operators <- as.list(operators); NULL } else .qclass_of(operators)
  if (!(is_set || .qisa(operators, "Operator") || (!is.null(oc) && .qsubclass(oc, "Operator"))))
    .qni("Argument is not an Operator or a set!")
  find_state <- function(key) { for (m in sm) if (.qset_key(m[[2L]]) == key) return(m[[1L]]); NULL }
  if (is_set) {
    cl <- lapply(operators, .qclass_of)
    for (k in seq_along(operators))
      if (!(.qisa(operators[[k]], "Operator") || (!is.null(cl[[k]]) && .qsubclass(cl[[k]], "Operator"))))
        .qni("Set is not all Operators!")
    if (all(!vapply(cl, is.null, NA))) {
      st <- find_state(.qset_key(unlist(cl)))
      if (!is.null(st)) {
        return(tryCatch({
          inst <- lapply(cl, function(op) .qfunc(op, list()))
          .qget_state(st, inst, ...)
        }, NotImplementedError = function(e) st))
      }
      return(NULL)
    }
    st <- find_state(.qset_key(vapply(operators, function(o) o$type, "")))
    if (!is.null(st)) return(.qget_state(st, operators, ...))
    return(NULL)
  }
  if (!is.null(oc)) {
    st <- find_state(.qset_key(oc))
    if (!is.null(st)) {
      return(tryCatch(.qget_state(st, .qfunc(oc, list()), ...), NotImplementedError = function(e) st))
    }
  }
  if (is_sym(operators)) {
    st <- find_state(.qset_key(operators$type))
    if (!is.null(st)) return(.qget_state(st, operators, ...))
  }
  NULL
}

#' @rdname operators_to_state
#' @export
state_to_operators <- function(state, ...) {
  sm <- .qstate_mapping()
  sc <- .qclass_of(state)
  if (!(.qisa(state, "StateBase") || (!is.null(sc) && .qsubclass(sc, "StateBase"))))
    .qni("Argument is not a state!")
  find_ops <- function(cls) { for (m in sm) if (m[[1L]] == cls) return(m[[2L]]); NULL }
  ret <- NULL
  if (!is.null(sc) && !is.null(find_ops(sc))) {
    state_inst <- .qmake_default(sc)
    ops <- find_ops(sc)
    ret <- tryCatch(.qget_ops(state_inst, ops, ...), NotImplementedError = function(e) .qset_or_one(ops),
                    TypeError = function(e) .qset_or_one(ops))
  } else if (is_sym(state) && !is.null(find_ops(state$type))) {
    ret <- .qget_ops(state, find_ops(state$type), ...)
  } else if (is_sym(state) && .qisa(state, "BraBase") && !is.null(find_ops(.qmeth(state, "dual_class")))) {
    ret <- .qget_ops(state, find_ops(.qmeth(state, "dual_class")))
  } else if (!is.null(sc) && .qsubclass(sc, "BraBase") && !is.null(find_ops(.qmeth(sc, "dual_class")))) {
    ops <- find_ops(.qmeth(sc, "dual_class"))
    state_inst <- .qmake_default(sc)
    ret <- tryCatch(.qget_ops(state_inst, ops), NotImplementedError = function(e) .qset_or_one(ops),
                    TypeError = function(e) .qset_or_one(ops))
  }
  .qmake_set(ret)
}

.qset_or_one <- function(ops) if (length(ops) == 1L) ops else as.list(ops)

# operatorset._make_set: a tuple or list of operators becomes a set
.qmake_set <- function(ops) if (is.list(ops) && !is_sym(ops)) .qset_of(ops) else ops

.qget_state <- function(state_class, ops, ...) {
  tryCatch(.qmeth(state_class, "_operators_to_state", ops, ...),
           NotImplementedError = function(e) .qmake_default(state_class))
}

.qget_ops <- function(state_inst, op_classes, ...) {
  # an instance method called on the class itself
  if (is.character(state_inst))
    .qstop("TypeError", "_state_to_operators() missing 1 required positional argument: 'op_classes'")
  ret <- tryCatch(.qmeth(state_inst, "_state_to_operators", if (length(op_classes) == 1L) op_classes else as.list(op_classes), ...),
                  NotImplementedError = function(e) {
                    if (length(op_classes) > 1L) lapply(op_classes, .qmake_default) else .qmake_default(op_classes)
                  })
  if (is.list(ret) && !is_sym(ret) && length(ret) == 1L) return(ret[[1L]])
  ret
}
