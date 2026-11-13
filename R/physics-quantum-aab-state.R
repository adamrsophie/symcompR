# States, operators, inner and outer products, tensor products,
# commutators and the constructor transforms of quantum products, after
# sympy/physics/quantum/state.py, operator.py, innerproduct.py,
# tensorproduct.py, commutator.py, anticommutator.py, constants.py and
# transforms.py.

# ---------------------------------------------------------------------------
# States

.qclass("StateBase", "QExpr",
  `_operators_to_state` = .cm(function(cls, ops, ...) .qni("Cannot map operators to states in this class. Method not implemented!")),
  `_state_to_operators` = function(self, op_classes, ...) .qni("Cannot map this state to operators. Method not implemented!"),
  operators = .prop(function(self) state_to_operators(self)),
  `_enumerate_state` = function(self, num_states, ...) .qni("Cannot enumerate this state!"),
  `_represent_default_basis` = function(self, ...) .qmeth(self, "_represent", basis = .qa(self, "operators"), ...),
  `_apply_operator` = function(self, op, ...) NULL,
  dual = .prop(function(self) .qmeth(.qmeth(self$type, "dual_class"), "_new_rawargs", .qhs(self), self$args)),
  dual_class = .cm(function(cls) .qni("dual_class must be implemented in a subclass")),
  `_eval_adjoint` = function(self) .qa(self, "dual"),
  `_sympystr` = function(self, p) {
    contents <- .qmeth(self, "_print_contents", p)
    paste0(.qattr_or(self, "lbracket", ""), contents, .qattr_or(self, "rbracket", ""))
  },
  `_latex` = function(self, p) {
    contents <- .qmeth(self, "_print_contents_latex", p)
    paste0("{", .qattr_or(self, "lbracket_latex", ""), contents, .qattr_or(self, "rbracket_latex", ""), "}")
  })

.qattr_or <- function(obj, nm, default) if (.qhas(obj, nm)) .qa(obj, nm) else default

.qclass("KetBase", "StateBase",
  kind = "KetKind",
  lbracket = "|", rbracket = ">",
  lbracket_latex = "\\left|", rbracket_latex = "\\right\\rangle ",
  default_args = .cm(function(cls) list("psi")),
  dual_class = .cm(function(cls) "BraBase"),
  `_eval_innerproduct` = function(self, bra, ...) .qdispatch(self, "_eval_innerproduct", bra, ...),
  `_apply_from_right_to` = function(self, op, ...) .qdispatch(self, "_apply_from_right_to", op, ...))

.qclass("BraBase", "StateBase",
  kind = "BraKind",
  lbracket = "<", rbracket = "|",
  lbracket_latex = "\\left\\langle ", rbracket_latex = "\\right|",
  `_operators_to_state` = .cm(function(cls, ops, ...) {
    state <- .qmeth(.qmeth(cls, "dual_class"), "_operators_to_state", ops, ...)
    .qa(state, "dual")
  }),
  `_state_to_operators` = function(self, op_classes, ...) .qmeth(.qa(self, "dual"), "_state_to_operators", op_classes, ...),
  `_enumerate_state` = function(self, num_states, ...) {
    lapply(.qmeth(.qa(self, "dual"), "_enumerate_state", num_states, ...), function(x) .qa(x, "dual"))
  },
  default_args = .cm(function(cls) .qmeth(.qmeth(cls, "dual_class"), "default_args")),
  dual_class = .cm(function(cls) "KetBase"),
  `_represent` = function(self, ...) Dagger(.qmeth(.qa(self, "dual"), "_represent", ...)))

.qclass("State", "StateBase")
.qclass("Ket", c("State", "KetBase"), dual_class = .cm(function(cls) "Bra"))
.qclass("Bra", c("State", "BraBase"), dual_class = .cm(function(cls) "Ket"))

.qclass("TimeDepState", "StateBase",
  default_args = .cm(function(cls) list("psi", "t")),
  label = .prop(function(self) self$args[-length(self$args)]),
  time = .prop(function(self) self$args[[length(self$args)]]),
  `_print_time` = function(self, p) .qpr(p, .qa(self, "time")),
  `_print_time_repr` = function(self, p) .qpr(p, .qa(self, "time")),
  `_print_time_latex` = function(self, p) .qpr(p, .qa(self, "time")),
  `_print_contents` = function(self, p) paste0(.qmeth(self, "_print_label", p), ";", .qmeth(self, "_print_time", p)),
  `_print_label_repr` = function(self, p) {
    paste0(.qmeth(self, "_print_sequence", .qa(self, "label"), ",", p), ",", .qmeth(self, "_print_time_repr", p))
  },
  `_print_contents_latex` = function(self, p) {
    paste0(.qmeth(self, "_print_sequence", .qa(self, "label"), .qa(self, "_label_separator"), p), ";",
           .qmeth(self, "_print_time_latex", p))
  })
.qclass("TimeDepKet", c("TimeDepState", "KetBase"), dual_class = .cm(function(cls) "TimeDepBra"))
.qclass("TimeDepBra", c("TimeDepState", "BraBase"), dual_class = .cm(function(cls) "TimeDepKet"))

.qclass("OrthogonalState", "State")
.qclass("OrthogonalKet", c("OrthogonalState", "KetBase"),
  dual_class = .cm(function(cls) "OrthogonalBra"),
  `_eval_innerproduct` = function(self, bra, ...) {
    if (length(self$args) != length(bra$args))
      .qstop("ValueError", "Cannot multiply a ket that has a different number of labels.")
    for (i in seq_along(self$args)) {
      d <- expand(.sub(self$args[[i]], bra$args[[i]]))
      is_zero <- ask_fact(d, "zero")
      if (identical(is_zero, FALSE)) return(S_Zero)
      if (is.na(is_zero)) return(NULL)
    }
    S_One
  })
.qclass("OrthogonalBra", c("OrthogonalState", "BraBase"), dual_class = .cm(function(cls) "OrthogonalKet"))

#' Quantum states
#'
#' Kets and bras of `sympy.physics.quantum.state`. `Ket(...)` and `Bra(...)`
#' take labels (strings become symbols); with no label the label is `psi`.
#' `TimeDepKet()` and `TimeDepBra()` take a time as the last label;
#' `OrthogonalKet()` and `OrthogonalBra()` are orthogonal for different
#' labels. `QExpr()`, `StateBase()`, `KetBase()`, `BraBase()`, `State()` and
#' `TimeDepState()` are the base classes. `Bra * Ket` is an [InnerProduct()],
#' `Ket * Bra` an [OuterProduct()], and `Dagger()` turns one into the other.
#'
#' @param ... Labels.
#' @return A quantum state.
#' @examples
#' k <- Ket("psi"); b <- Bra("phi")
#' b * k
#' k * b
#' Dagger(k)
#' TimeDepKet("psi", "t")
#' @name quantum-states
NULL

#' @rdname quantum-states
#' @export
QExpr <- function(...) .qfunc("QExpr", list(...))
#' @rdname quantum-states
#' @export
StateBase <- function(...) .qfunc("StateBase", list(...))
#' @rdname quantum-states
#' @export
KetBase <- function(...) .qfunc("KetBase", list(...))
#' @rdname quantum-states
#' @export
BraBase <- function(...) .qfunc("BraBase", list(...))
#' @rdname quantum-states
#' @export
State <- function(...) .qfunc("State", list(...))
#' @rdname quantum-states
#' @export
Ket <- function(...) .qfunc("Ket", list(...))
#' @rdname quantum-states
#' @export
Bra <- function(...) .qfunc("Bra", list(...))
#' @rdname quantum-states
#' @export
TimeDepState <- function(...) .qfunc("TimeDepState", list(...))
#' @rdname quantum-states
#' @export
TimeDepKet <- function(...) .qfunc("TimeDepKet", list(...))
#' @rdname quantum-states
#' @export
TimeDepBra <- function(...) .qfunc("TimeDepBra", list(...))
#' @rdname quantum-states
#' @export
OrthogonalState <- function(...) .qfunc("OrthogonalState", list(...))
#' @rdname quantum-states
#' @export
OrthogonalKet <- function(...) .qfunc("OrthogonalKet", list(...))
#' @rdname quantum-states
#' @export
OrthogonalBra <- function(...) .qfunc("OrthogonalBra", list(...))

# ---------------------------------------------------------------------------
# Wavefunction (a Function subclass)

.qclass("Wavefunction", "Function",
  `__new__` = .cm(function(cls, args) {
    args <- lapply(args, function(a) if (is.list(a) && !is_sym(a)) do.call(Tuple, a) else .qS(a))
    .qm_state$on <- TRUE
    .nc_state$any <- TRUE
    new_sym("Wavefunction", args)
  }),
  is_commutative = FALSE,
  `__call__` = function(self, ...) {
    dots <- list(...)
    kw <- if (is.null(names(dots))) logical(length(dots)) else nzchar(names(dots))
    options <- dots[kw]
    args <- unname(dots[!kw])
    var <- .qa(self, "variables")
    if (length(args) != length(var)) .qni("Incorrect number of arguments to function!")
    lims <- .qa(self, "limits")
    ct <- 1L
    for (k in seq_along(var)) {
      lu <- lims[[k]][[2L]]
      lower <- lu[[1L]]; upper <- lu[[2L]]
      a <- args[[ct]]
      # an Expr argument is compared only when a limit is among its free
      # symbols (a plain R number is a Python number and always is)
      if (is_sym(a)) {
        fs <- names(.free_sym_list(a))
        if (!(lower$key %in% fs || upper$key %in% fs)) next
      }
      if (isTRUE(Lt(a, lower) == S_true) || isTRUE(Gt(a, upper) == S_true)) return(S_Zero)
      ct <- ct + 1L
    }
    ex <- .qa(self, "expr")
    # f(2, 4, m=1, n=1): keyword values replace the symbols of the same name
    for (sy in .free_sym_list(ex)) if (sy$name %in% names(options)) ex <- subs(ex, sy, options[[sy$name]])
    subs(ex, lapply(seq_along(var), function(k) list(var[[k]], .qS(args[[k]]))))
  },
  `_eval_derivative` = function(self, symbol) {
    d <- diff(.qa(self, "expr"), symbol)
    .qfunc("Wavefunction", c(list(d), self$args[-1L]))
  },
  `_eval_conjugate` = function(self) .qfunc("Wavefunction", c(list(fcall("conjugate", list(.qa(self, "expr")))), self$args[-1L])),
  `_eval_transpose` = function(self) self,
  variables = .prop(function(self) lapply(self$args[-1L], function(g) if (g$type == "Tuple") g$args[[1L]] else g)),
  # a dict {variable: (lower, upper)}: a list of pairs with attribute qdict
  limits = .prop(function(self) structure(lapply(self$args[-1L], function(g)
    if (g$type == "Tuple") list(g$args[[1L]], list(g$args[[2L]], g$args[[3L]])) else list(g, list(S_NegativeInfinity, S_Infinity))),
    qdict = TRUE)),
  expr = .prop(function(self) self$args[[1L]]),
  is_normalized = .prop(function(self) .qeq(.qa(self, "norm"), S_One)),
  norm = .prop(function(self) {
    ex <- .mul(list(.qa(self, "expr"), fcall("conjugate", list(.qa(self, "expr")))))
    var <- .qa(self, "variables")
    lims <- .qa(self, "limits")
    for (k in seq_along(var)) ex <- integrate(ex, list(var[[k]], lims[[k]][[2L]][[1L]], lims[[k]][[2L]][[2L]]))
    .sqrt(ex)
  }),
  normalize = function(self) {
    const <- .qa(self, "norm")
    if (identical(const$key, S_Infinity$key)) .qni("The function is not normalizable!")
    .qfunc("Wavefunction", c(list(.mul(list(.pow(const, S_NegativeOne), .qa(self, "expr")))), self$args[-1L]))
  },
  prob = function(self) {
    .qfunc("Wavefunction", c(list(.mul(list(.qa(self, "expr"), fcall("conjugate", list(.qa(self, "expr")))))),
                             .qa(self, "variables")))
  },
  `_sympystr` = function(self, p) paste0("Wavefunction(", paste(vapply(self$args, .str, ""), collapse = ", "), ")"),
  `_sympyrepr` = function(self, p) paste0("Wavefunction(", paste(vapply(self$args, srepr, ""), collapse = ", "), ")"),
  `_latex` = function(self, p) paste0("\\operatorname{Wavefunction}{\\left(", paste(vapply(self$args, .latex, ""), collapse = ", "), " \\right)}"))
register_diff_hook("Wavefunction", function(e, x) .qmeth(e, "_eval_derivative", x))

#' Wavefunctions
#'
#' `Wavefunction(expr, ...)` is a function of the coordinates given after the
#' expression, each a symbol or `list(x, lower, upper)`, as SymPy's
#' `sympy.physics.quantum.state.Wavefunction`. `qcall(f, "__call__", ...)`
#' evaluates it (zero outside the limits), `qattr(f, "norm")`,
#' `qcall(f, "normalize")` and `qcall(f, "prob")` normalize it.
#'
#' @param ... The expression, then the coordinates.
#' @return A `Wavefunction`.
#' @examples
#' x <- Symbol("x"); L <- Symbol("L", positive = TRUE)
#' f <- Wavefunction(sqrt(2 / L) * sin(sym_pi * x / L), list(x, 0, L))
#' qattr(f, "norm")
#' @export
Wavefunction <- function(...) .qfunc("Wavefunction", list(...))

# ---------------------------------------------------------------------------
# Operators

.qclass("Operator", "QExpr",
  is_hermitian = NULL, is_unitary = NULL,
  default_args = .cm(function(cls) list("O")),
  kind = "OperatorKind",
  `_label_separator` = ",",
  `_print_operator_name` = function(self, p) .qpyname(self$type),
  `_print_operator_name_latex` = function(self, p) .qpyname(self$type),
  `_print_contents` = function(self, p) {
    if (length(.qa(self, "label")) == 1L) .qmeth(self, "_print_label", p)
    else paste0(.qmeth(self, "_print_operator_name", p), "(", .qmeth(self, "_print_label", p), ")")
  },
  `_print_contents_latex` = function(self, p) {
    if (length(.qa(self, "label")) == 1L) .qmeth(self, "_print_label_latex", p)
    else paste0(.qmeth(self, "_print_operator_name_latex", p), "\\left(", .qmeth(self, "_print_label_latex", p), "\\right)")
  },
  `_eval_commutator` = function(self, other, ...) .qdispatch(self, "_eval_commutator", other, ...),
  `_eval_anticommutator` = function(self, other, ...) .qdispatch(self, "_eval_anticommutator", other, ...),
  `_apply_operator` = function(self, ket, ...) .qdispatch(self, "_apply_operator", ket, ...),
  `_apply_from_right_to` = function(self, bra, ...) NULL,
  matrix_element = function(self, ...) .qni("matrix_elements is not defined"),
  inverse = function(self) .qmeth(self, "_eval_inverse"),
  inv = function(self) .qmeth(self, "_eval_inverse"),
  `_eval_inverse` = function(self) .pow(self, S_NegativeOne))

.qclass("HermitianOperator", "Operator",
  is_hermitian = TRUE,
  `_eval_inverse` = function(self) {
    if (.qisa(self, "UnitaryOperator")) self else .qsuper("Operator", "_eval_inverse", self)
  },
  `_eval_power` = function(self, ex) {
    if (.qisa(self, "UnitaryOperator")) {
      if (isTRUE(ask_fact(ex, "even"))) return(S_One)
      if (isTRUE(ask_fact(ex, "odd"))) return(self)
    }
    NULL
  })

.qclass("UnitaryOperator", "Operator",
  is_unitary = TRUE,
  `_eval_adjoint` = function(self) .qmeth(self, "_eval_inverse"))

.qclass("IdentityOperator", "Operator",
  is_hermitian = TRUE, is_unitary = TRUE,
  `__new__` = .cm(function(cls, args, ...) {
    if (!(length(args) %in% c(0L, 1L))) .qstop("ValueError", "0 or 1 parameters expected, got ", .qrepr(args))
    .qsuper("QExpr", "__new__", cls, args, ...)
  }),
  dimension = .prop(function(self) .qa(self, "N")),
  N = .prop(function(self) {
    if (length(self$args) == 1L && .qtruthy(self$args[[1L]])) self$args[[1L]] else S_Infinity
  }),
  default_args = .cm(function(cls) list(S_Infinity)),
  `_eval_commutator` = function(self, other, ...) S_Zero,
  `_eval_anticommutator` = function(self, other, ...) .mul(list(S_Two, other)),
  `_eval_inverse` = function(self) self,
  `_eval_adjoint` = function(self) self,
  `_apply_operator` = function(self, ket, ...) ket,
  `_apply_from_right_to` = function(self, bra, ...) bra,
  `_eval_power` = function(self, ex) self,
  `_print_contents` = function(self, p) "I",
  `_print_contents_latex` = function(self, p) "{\\mathcal{I}}",
  `_represent_default_basis` = function(self, format = "sympy", ...) {
    N <- .qa(self, "N")
    if (identical(N$key, S_Infinity$key))
      .qni("Cannot represent infinite dimensional identity operator as a matrix")
    if (format != "sympy") .qni("Representation in format ", format, " not implemented.")
    eye(.qint(N))
  })

.qclass("OuterProduct", "Operator",
  is_commutative = FALSE,
  `__new__` = .cm(function(cls, args, ...) {
    if (length(args) != 2L) .qstop("ValueError", "2 parameters expected, got ", length(args))
    ket_expr <- expand(.qS(args[[1L]]))
    bra_expr <- expand(.qS(args[[2L]]))
    if ((.qisa(ket_expr, "KetBase") || ket_expr$type == "Mul") && (.qisa(bra_expr, "BraBase") || bra_expr$type == "Mul")) {
      kc <- .args_cnc(ket_expr); bc <- .args_cnc(bra_expr)
      if (length(kc$nc) != 1L || !.qisa(kc$nc[[1L]], "KetBase"))
        .qstop("TypeError", "KetBase subclass expected, got: ", .qrepr(.mul(kc$nc)))
      if (length(bc$nc) != 1L || !.qisa(bc$nc[[1L]], "BraBase"))
        .qstop("TypeError", "BraBase subclass expected, got: ", .qrepr(.mul(bc$nc)))
      if (!identical(.qmeth(kc$nc[[1L]], "dual_class"), bc$nc[[1L]]$type))
        .qstop("TypeError", "ket and bra are not dual classes: ", .qpyrepr_class(kc$nc[[1L]]$type), ", ",
               .qpyrepr_class(bc$nc[[1L]]$type))
      obj <- .qraw(cls, .qhs(kc$nc[[1L]]), list(kc$nc[[1L]], bc$nc[[1L]]))
      return(.mul(list(.mul(c(kc$c, bc$c)), obj)))
    }
    op_terms <- list()
    if (ket_expr$type == "Add" && bra_expr$type == "Add") {
      for (kt in .sympy_args(ket_expr)) for (bt in .sympy_args(bra_expr))
        op_terms[[length(op_terms) + 1L]] <- OuterProduct(kt, bt)
    } else if (ket_expr$type == "Add") {
      for (kt in .sympy_args(ket_expr)) op_terms[[length(op_terms) + 1L]] <- OuterProduct(kt, bra_expr)
    } else if (bra_expr$type == "Add") {
      for (bt in .sympy_args(bra_expr)) op_terms[[length(op_terms) + 1L]] <- OuterProduct(ket_expr, bt)
    } else {
      .qstop("TypeError", "Expected ket and bra expression, got: ", .qrepr(ket_expr), ", ", .qrepr(bra_expr))
    }
    .add(op_terms)
  }),
  ket = .prop(function(self) self$args[[1L]]),
  bra = .prop(function(self) self$args[[2L]]),
  label = .prop(function(self) self$args),
  `_eval_adjoint` = function(self) OuterProduct(Dagger(.qa(self, "bra")), Dagger(.qa(self, "ket"))),
  `_sympystr` = function(self, p) paste0(.qpr(p, self$args[[1L]]), .qpr(p, self$args[[2L]])),
  `_sympyrepr` = function(self, p) paste0("OuterProduct(", .qpr(p, self$args[[1L]]), ",", .qpr(p, self$args[[2L]]), ")"),
  `_latex` = function(self, p) paste0(.qpr(p, self$args[[1L]]), .qpr(p, self$args[[2L]])),
  `_represent` = function(self, ...) {
    k <- .qmeth(self$args[[1L]], "_represent", ...)
    b <- .qmeth(self$args[[2L]], "_represent", ...)
    .qtimes(k, b)
  },
  `_eval_trace` = function(self, ...) .qmeth(self$args[[1L]], "_eval_trace", self$args[[2L]], ...))

.qpyrepr_class <- function(cls) {
  mod <- .qa(cls, "__module__")
  paste0("<class '", mod, ".", .qpyname(cls), "'>")
}
.qdef("Basic", `__module__` = "sympy.physics.quantum.state")
.qdef("Operator", `__module__` = "sympy.physics.quantum.operator")

.qclass("DifferentialOperator", "Operator",
  variables = .prop(function(self) self$args[[length(self$args)]]$args),
  `function` = .prop(function(self) self$args[[length(self$args)]]),
  expr = .prop(function(self) self$args[[1L]]),
  free_symbols = .prop(function(self) free_symbols(self$args[[1L]])),
  `_apply_operator_Wavefunction` = function(self, func, ...) {
    var <- .qa(self, "variables")
    wf_vars <- func$args[-1L]
    f <- .qa(self, "function")
    new_expr <- subs(.qa(self, "expr"), f, do.call(.qa(func, "__call__"), var))
    new_expr <- doit(new_expr)
    .qfunc("Wavefunction", c(list(new_expr), wf_vars))
  },
  `_eval_derivative` = function(self, symbol) {
    new_expr <- Derivative(.qa(self, "expr"), symbol)
    .qfunc("DifferentialOperator", list(new_expr, self$args[[length(self$args)]]))
  })
register_diff_hook("DifferentialOperator", function(e, x) .qmeth(e, "_eval_derivative", x))
# the free symbols of a DifferentialOperator are those of its expression
register_bound_vars("DifferentialOperator", function(e) {
  inner <- names(.free_sym_list(e$args[[length(e$args)]]))
  setdiff(inner, names(.free_sym_list(e$args[[1L]])))
})

#' Quantum operators
#'
#' The operators of `sympy.physics.quantum.operator`: a general
#' `Operator(label)` (`O` by default), `HermitianOperator()`,
#' `UnitaryOperator()`, the deprecated `IdentityOperator(N)`, the outer
#' product `OuterProduct(ket, bra)` (also `ket * bra`) and
#' `DifferentialOperator(expr, f)`, which applies `expr` with `f` replaced by
#' a [Wavefunction()] under [qapply()]. Operators do not commute; `qcall(A,
#' "inverse")` is `A^-1`.
#'
#' @param ... Labels or arguments.
#' @return An operator.
#' @examples
#' A <- Operator("A"); B <- Operator("B")
#' A * B - B * A
#' Dagger(UnitaryOperator("U")) * UnitaryOperator("U")
#' OuterProduct(Ket("a"), Bra("b"))
#' @name quantum-operators
NULL

#' @rdname quantum-operators
#' @export
Operator <- function(...) .qfunc("Operator", list(...))
#' @rdname quantum-operators
#' @export
HermitianOperator <- function(...) .qfunc("HermitianOperator", list(...))
#' @rdname quantum-operators
#' @export
UnitaryOperator <- function(...) .qfunc("UnitaryOperator", list(...))
#' @rdname quantum-operators
#' @export
IdentityOperator <- function(...) .qfunc("IdentityOperator", list(...))
#' @rdname quantum-operators
#' @export
OuterProduct <- function(...) .qfunc("OuterProduct", list(...))
#' @rdname quantum-operators
#' @export
DifferentialOperator <- function(...) .qfunc("DifferentialOperator", list(...))

# ---------------------------------------------------------------------------
# InnerProduct

.qclass("InnerProduct", "Expr",
  `__new__` = .cm(function(cls, args) {
    bra <- args[[1L]]; ket <- args[[2L]]
    if (!.qisa(ket, "KetBase")) .qstop("TypeError", "KetBase subclass expected, got: ", .qrepr(ket))
    if (!.qisa(bra, "BraBase")) .qstop("TypeError", "BraBase subclass expected, got: ", .qrepr(ket))
    .qm_state$on <- TRUE
    new_sym("InnerProduct", list(bra, ket))
  }),
  kind = "NumberKind",
  is_complex = TRUE,
  is_commutative = TRUE,
  bra = .prop(function(self) self$args[[1L]]),
  ket = .prop(function(self) self$args[[2L]]),
  `_eval_conjugate` = function(self) InnerProduct(Dagger(self$args[[2L]]), Dagger(self$args[[1L]])),
  `_sympyrepr` = function(self, p) paste0("InnerProduct(", .qpr(p, self$args[[1L]]), ",", .qpr(p, self$args[[2L]]), ")"),
  `_sympystr` = function(self, p) {
    sbra <- .qpr(p, self$args[[1L]]); sket <- .qpr(p, self$args[[2L]])
    paste0(substr(sbra, 1L, nchar(sbra) - 1L), "|", substr(sket, 2L, nchar(sket)))
  },
  `_latex` = function(self, p) {
    bra_label <- .qmeth(self$args[[1L]], "_print_contents_latex", p)
    paste0("\\left\\langle ", bra_label, " \\right. ", .qpr(p, self$args[[2L]]))
  },
  doit = function(self, ...) {
    r <- tryCatch(.qmeth(self$args[[2L]], "_eval_innerproduct", self$args[[1L]], ...),
                  NotImplementedError = function(e) {
                    tryCatch(fcall("conjugate", list(.qmeth(.qa(self$args[[1L]], "dual"), "_eval_innerproduct",
                                                         .qa(self$args[[2L]], "dual"), ...))),
                             NotImplementedError = function(e2) NULL)
                  })
    if (!is.null(r)) r else self
  })
# InnerProduct is_complex: complex, so finite as well
.facts_hooks$InnerProduct <- function(e) c(commutative = TRUE, complex = TRUE)

#' Inner products
#'
#' `InnerProduct(bra, ket)` is the unevaluated `<bra|ket>`; [doit()]
#' evaluates it when the states know their inner product.
#'
#' @param bra A bra.
#' @param ket A ket.
#' @return An `InnerProduct`.
#' @examples
#' InnerProduct(Bra("a"), Ket("b"))
#' doit(InnerProduct(OrthogonalBra(1), OrthogonalKet(2)))
#' @export
InnerProduct <- function(bra, ket) .qfunc("InnerProduct", list(bra, ket))

# ---------------------------------------------------------------------------
# TensorProduct (sym type QTensorProduct; diffgeom has its own)

.qtp_flatten <- function(args) {
  c_part <- list(); nc_parts <- list()
  for (arg in args) {
    cn <- .args_cnc(arg)
    c_part <- c(c_part, cn$c)
    nc_parts[[length(nc_parts) + 1L]] <- .qmul_raw(cn$nc)
  }
  list(c = c_part, nc = nc_parts)
}

.qtp <- function(args) {
  if (inherits(args[[1L]], "symmatrix")) return(.matrix_tensor_product(args))
  args <- lapply(args, .qS)
  fl <- .qtp_flatten(args)
  c_part <- .mul(fl$c)
  new_args <- fl$nc
  if (!length(new_args)) return(c_part)
  if (length(new_args) == 1L) return(.mul(list(c_part, new_args[[1L]])))
  .nc_state$any <- TRUE
  .qm_state$on <- TRUE
  tp <- new_sym("QTensorProduct", new_args)
  .mul(list(c_part, tp))
}

.matrix_tensor_product <- function(mats) {
  r <- Matrix(mats[[1L]])
  for (m in mats[-1L]) r <- kronecker_product(r, Matrix(m))
  r
}

.qclass("QTensorProduct", "Expr", pyname = "TensorProduct",
  `__new__` = .cm(function(cls, args) .qtp(args)),
  is_commutative = FALSE,
  kind = .prop(function(self) .kind_dispatch(vapply(self$args, .qkind, ""), .KIND_TP)),
  `_eval_adjoint` = function(self) .qtp(lapply(self$args, Dagger)),
  `_sympystr` = function(self, p) {
    n <- length(self$args)
    s <- ""
    for (i in seq_len(n)) {
      a <- self$args[[i]]
      par <- a$type %in% c("Add", "Pow", "Mul")
      s <- paste0(s, if (par) "(", .qpr(p, a), if (par) ")")
      if (i != n && !.qisa(a, c("KetBase", "BraBase"))) s <- paste0(s, "x")
    }
    s
  },
  `_sympyrepr` = function(self, p) paste0("TensorProduct(", paste(vapply(self$args, srepr, ""), collapse = ", "), ")"),
  `_latex` = function(self, p) {
    if (isTRUE(.qm_state$combined) && (all(vapply(self$args, function(a) .qisa(a, "Ket"), NA)) ||
                                       all(vapply(self$args, function(a) .qisa(a, "Bra"), NA)))) {
      lw <- function(label, n) if (n == 1L) label else paste0("\\left\\{", label, "\\right\\}")
      s <- paste(vapply(self$args, function(a) lw(.qmeth(a, "_print_label_latex", p), length(a$args)), ""), collapse = ", ")
      return(paste0("{", .qa(self$args[[1L]], "lbracket_latex"), s, .qa(self$args[[1L]], "rbracket_latex"), "}"))
    }
    n <- length(self$args)
    s <- ""
    for (i in seq_len(n)) {
      a <- self$args[[i]]
      par <- a$type %in% c("Add", "Mul")
      s <- paste0(s, if (par) "\\left(", "{", .qpr(p, a), "}", if (par) "\\right)")
      if (i != n && !.qisa(a, c("KetBase", "BraBase"))) s <- paste0(s, "\\otimes ")
    }
    s
  },
  doit = function(self, ...) .qtp(lapply(self$args, doit)),
  `_eval_expand_tensorproduct` = function(self, ...) {
    args <- self$args
    add_args <- list()
    for (i in seq_along(args)) {
      if (args[[i]]$type == "Add") {
        for (aa in .sympy_args(args[[i]])) {
          na <- args; na[[i]] <- aa
          tp <- .qtp(na)
          cn <- .args_cnc(tp)
          nc_part <- cn$nc
          if (length(nc_part) == 1L && nc_part[[1L]]$type == "QTensorProduct")
            nc_part <- list(.qmeth(nc_part[[1L]], "_eval_expand_tensorproduct"))
          add_args[[length(add_args) + 1L]] <- .mul(list(.mul(cn$c), .mul(nc_part)))
        }
        break
      }
    }
    if (length(add_args)) .add(add_args) else self
  },
  `_eval_trace` = function(self, indices = NULL, ...) {
    if (is.null(indices) || !length(indices)) return(.mul(lapply(self$args, function(a) doit(Tr(a)))))
    .mul(lapply(seq_along(self$args), function(k) {
      v <- self$args[[k]]
      if ((k - 1L) %in% vapply(indices, .qint, 0)) doit(Tr(v)) else v
    }))
  })
.expand_custom$tensorproduct <- function(e) if (e$type == "QTensorProduct") .qmeth(e, "_eval_expand_tensorproduct") else NULL

#' @noRd
.qTensorProduct <- function(...) .qtp(list(...))

#' Combined printing of tensor products
#'
#' `combined_tensor_printing(TRUE)` prints a tensor product of kets (or of
#' bras) in LaTeX as one ket with all the labels, as SymPy's
#' `combined_tensor_printing`.
#'
#' @param combined `TRUE` or `FALSE`.
#' @return `NULL`, invisibly.
#' @examples
#' combined_tensor_printing(FALSE)
#' @export
combined_tensor_printing <- function(combined) {
  .qm_state$combined <- isTRUE(combined)
  invisible(NULL)
}

#' @rdname combined_tensor_printing
#' @param e An expression.
#' @param ... Ignored hints.
#' @export
tensor_product_simp <- function(e, ...) e

# ---------------------------------------------------------------------------
# Commutator and AntiCommutator (sym types QCommutator and AntiCommutator)

.qcomm_new <- function(cls, a, b, anti) {
  a <- .qS(a); b <- .qS(b)
  if (!(.qtruthy(a) && .qtruthy(b))) return(S_Zero)
  if (identical(a$key, b$key)) return(if (anti) .mul(list(S_Two, .pow(a, S_Two))) else S_Zero)
  if (.is_comm(a) || .is_comm(b)) return(if (anti) .mul(list(S_Two, a, b)) else S_Zero)
  ca <- .args_cnc(a); cb <- .args_cnc(b)
  c_part <- c(ca$c, cb$c)
  if (length(c_part)) return(.mul(list(.mul(c_part), .qcomm_new(cls, .qmul_raw(ca$nc), .qmul_raw(cb$nc), anti))))
  if (.basic_compare(a, b) == 1L) {
    return(if (anti) .qcomm_new(cls, b, a, anti) else .mul(list(S_NegativeOne, .qcomm_new(cls, b, a, anti))))
  }
  .nc_state$any <- TRUE
  .qm_state$on <- TRUE
  new_sym(cls, list(a, b))
}

.qclass("QCommutator", "Expr", pyname = "Commutator",
  `__new__` = .cm(function(cls, args) .qcomm_new("QCommutator", args[[1L]], args[[2L]], FALSE)),
  is_commutative = FALSE,
  kind = .prop(function(self) .kind_dispatch(vapply(self$args, .qkind, ""), .KIND_COMM)),
  `_expand_pow` = function(self, A, B, sign) {
    ex <- A$args[[2L]]
    if (!isTRUE(ask_fact(ex, "integer")) || length(free_symbols(ex)) || abs(.qint(ex)) <= 1) return(self)
    base <- A$args[[1L]]
    if (isTRUE(ask_fact(ex, "negative"))) {
      base <- .pow(A$args[[1L]], S_NegativeOne)
      ex <- .neg(ex)
    }
    comm <- expand(.qcomm_new("QCommutator", base, B, FALSE), commutator = TRUE)
    n <- .qint(ex)
    result <- .mul(list(.pow(base, .int(n - 1)), comm))
    for (i in seq_len(n - 1L))
      result <- .add(list(result, .mul(list(.pow(base, .int(n - 1 - i)), comm, .pow(base, .int(i))))))
    .mul(list(sign, expand(result)))
  },
  `_eval_expand_commutator` = function(self, ...) {
    A <- self$args[[1L]]; B <- self$args[[2L]]
    sub_expand <- function(cm) if (is_sym(cm) && cm$type == "QCommutator") .qmeth(cm, "_eval_expand_commutator") else cm
    if (A$type == "Add") {
      return(.add(lapply(.sympy_args(A), function(term) sub_expand(.qcomm_new("QCommutator", term, B, FALSE)))))
    } else if (B$type == "Add") {
      return(.add(lapply(.sympy_args(B), function(term) sub_expand(.qcomm_new("QCommutator", A, term, FALSE)))))
    } else if (A$type == "Mul") {
      aa <- .qargs(A)
      a <- aa[[1L]]; b <- .mul(aa[-1L]); c <- B
      comm1 <- sub_expand(.qcomm_new("QCommutator", b, c, FALSE))
      comm2 <- sub_expand(.qcomm_new("QCommutator", a, c, FALSE))
      return(.add(list(.mul(list(a, comm1)), .mul(list(comm2, b)))))
    } else if (B$type == "Mul") {
      ba <- .qargs(B)
      a <- A; b <- ba[[1L]]; c <- .mul(ba[-1L])
      comm1 <- sub_expand(.qcomm_new("QCommutator", a, b, FALSE))
      comm2 <- sub_expand(.qcomm_new("QCommutator", a, c, FALSE))
      return(.add(list(.mul(list(comm1, c)), .mul(list(b, comm2)))))
    } else if (A$type == "Pow") {
      return(.qmeth(self, "_expand_pow", A, B, S_One))
    } else if (B$type == "Pow") {
      return(.qmeth(self, "_expand_pow", B, A, S_NegativeOne))
    }
    self
  },
  doit = function(self, ...) {
    A <- self$args[[1L]]; B <- self$args[[2L]]
    if (.qisa(A, "Operator") && .qisa(B, "Operator")) {
      comm <- tryCatch(.qmeth(A, "_eval_commutator", B, ...), NotImplementedError = function(e)
        tryCatch(.mul(list(S_NegativeOne, .qmeth(B, "_eval_commutator", A, ...))), NotImplementedError = function(e2) NULL))
      if (!is.null(comm)) return(.qdoit(comm))
    }
    .qdoit(.sub(.qtimes(A, B), .qtimes(B, A)))
  },
  `_eval_adjoint` = function(self) .qcomm_new("QCommutator", Dagger(self$args[[2L]]), Dagger(self$args[[1L]]), FALSE),
  `_sympyrepr` = function(self, p) paste0("Commutator(", .qpr(p, self$args[[1L]]), ",", .qpr(p, self$args[[2L]]), ")"),
  `_sympystr` = function(self, p) paste0("[", .qpr(p, self$args[[1L]]), ",", .qpr(p, self$args[[2L]]), "]"),
  `_latex` = function(self, p) paste0("\\left[", .qpr(p, self$args[[1L]]), ",", .qpr(p, self$args[[2L]]), "\\right]"))
.expand_custom$commutator <- function(e) if (e$type == "QCommutator") .qmeth(e, "_eval_expand_commutator") else NULL

.qclass("AntiCommutator", "Expr",
  `__new__` = .cm(function(cls, args) .qcomm_new("AntiCommutator", args[[1L]], args[[2L]], TRUE)),
  is_commutative = FALSE,
  kind = .prop(function(self) .kind_dispatch(vapply(self$args, .qkind, ""), .KIND_COMM)),
  doit = function(self, ...) {
    A <- self$args[[1L]]; B <- self$args[[2L]]
    if (.qisa(A, "Operator") && .qisa(B, "Operator")) {
      comm <- tryCatch(.qmeth(A, "_eval_anticommutator", B, ...), NotImplementedError = function(e)
        tryCatch(.qmeth(B, "_eval_anticommutator", A, ...), NotImplementedError = function(e2) NULL))
      if (!is.null(comm)) return(.qdoit(comm))
    }
    .qdoit(.add(list(.qtimes(A, B), .qtimes(B, A))))
  },
  `_eval_adjoint` = function(self) .qcomm_new("AntiCommutator", Dagger(self$args[[1L]]), Dagger(self$args[[2L]]), TRUE),
  `_sympyrepr` = function(self, p) paste0("AntiCommutator(", .qpr(p, self$args[[1L]]), ",", .qpr(p, self$args[[2L]]), ")"),
  `_sympystr` = function(self, p) paste0("{", .qpr(p, self$args[[1L]]), ",", .qpr(p, self$args[[2L]]), "}"),
  `_latex` = function(self, p) paste0("\\left\\{", .qpr(p, self$args[[1L]]), ",", .qpr(p, self$args[[2L]]), "\\right\\}"))

.qdoit <- function(x) if (inherits(x, "symmatrix")) .mat_map(x, doit) else doit(x)

#' @noRd
.qCommutator <- function(A, B) .qcomm_new("QCommutator", A, B, FALSE)

#' Anticommutators
#'
#' `AntiCommutator(A, B)` is the unevaluated `{A, B} = A*B + B*A` of
#' `sympy.physics.quantum.anticommutator`; it is `2*A*B` when one argument
#' commutes and its arguments are put in canonical order. [doit()] evaluates
#' it. The commutator `[A, B]` is [Commutator()].
#'
#' @param A,B Operators or expressions.
#' @return An expression.
#' @examples
#' A <- Operator("A"); B <- Operator("B")
#' AntiCommutator(B, A)
#' doit(AntiCommutator(A, B))
#' @export
AntiCommutator <- function(A, B) .qcomm_new("AntiCommutator", A, B, TRUE)

# ---------------------------------------------------------------------------
# hbar

.qclass("HBar", "Expr",
  `__new__` = .cm(function(cls, args) hbar),
  is_commutative = TRUE, is_real = TRUE, is_positive = TRUE, is_negative = FALSE, is_irrational = TRUE, is_Atom = TRUE,
  `_sympyrepr` = function(self, p) "HBar()",
  `_sympystr` = function(self, p) "hbar",
  `_latex` = function(self, p) "\\hbar")
.facts_hooks$HBar <- function(e) c(commutative = TRUE, real = TRUE, positive = TRUE, negative = FALSE,
                                   irrational = TRUE, finite = TRUE)
register_neval("HBar", function(e, env) 1.05457162e-34)

#' Reduced Planck constant
#'
#' `hbar` is the numerical constant of `sympy.physics.quantum.constants`,
#' a positive real number that prints as `hbar`; `HBar()` returns it.
#'
#' @return The `hbar` constant.
#' @examples
#' 2 * hbar
#' HBar()
#' @export
HBar <- function() hbar
#' @rdname HBar
#' @export
hbar <- new_sym("HBar", list(), key = "HBar()")

# ---------------------------------------------------------------------------
# Constructor transforms (transforms.py)

.qargs <- function(e) if (e$type == "Mul" && !isTRUE(e$extra$raw)) .sympy_args(e) else e$args

.qtransform_pair <- function(a, b) {
  isa <- .qisa
  if (isa(a, "IdentityOperator") && isa(b, "IdentityOperator")) return(list(S_One))
  if (isa(a, "OuterProduct") && isa(b, "OuterProduct"))
    return(list(InnerProduct(a$args[[2L]], b$args[[1L]]), OuterProduct(a$args[[1L]], b$args[[2L]])))
  if (isa(a, "OuterProduct") && isa(b, "KetBase")) return(list(InnerProduct(a$args[[2L]], b), a$args[[1L]]))
  if (isa(a, "BraBase") && isa(b, "OuterProduct")) return(list(InnerProduct(a, b$args[[1L]]), b$args[[2L]]))
  if (isa(a, "QTensorProduct") && isa(b, "QTensorProduct")) {
    if (length(a$args) == length(b$args)) {
      if (.qkind(a) == "BraKind" && .qkind(b) == "KetKind")
        return(lapply(seq_along(a$args), function(k) InnerProduct(a$args[[k]], b$args[[k]])))
      return(list(.qtp(lapply(seq_along(a$args), function(k) .mul(list(a$args[[k]], b$args[[k]]))))))
    }
    return(NULL)
  }
  if (isa(a, "QTensorProduct") && isa(b, "KetBase")) {
    if (.qkind(a) == "KetKind") .qstop("TypeError", "Multiplication of TensorProduct(*kets)*ket is invalid.")
    return(NULL)
  }
  if (isa(a, "KetBase") && isa(b, "QTensorProduct")) {
    if (.qkind(b) == "KetKind") .qstop("TypeError", "Multiplication of ket*TensorProduct(*kets) is invalid.")
    return(NULL)
  }
  if (isa(a, "QTensorProduct") && isa(b, "BraBase")) {
    if (.qkind(a) == "BraKind") .qstop("TypeError", "Multiplication of TensorProduct(*bras)*bra is invalid.")
    return(NULL)
  }
  if (isa(a, "BraBase") && isa(b, "QTensorProduct")) {
    if (.qkind(b) == "BraKind") .qstop("TypeError", "Multiplication of bra*TensorProduct(*bras) is invalid.")
    return(NULL)
  }
  if (isa(a, "BraBase") && isa(b, "KetBase")) return(list(InnerProduct(a, b)))
  if (isa(a, "KetBase") && isa(b, "BraBase")) return(list(OuterProduct(a, b)))
  if (isa(a, "KetBase") && isa(b, "KetBase"))
    .qstop("TypeError", "Multiplication of two kets is not allowed. Use TensorProduct instead.")
  if (isa(a, "BraBase") && isa(b, "BraBase"))
    .qstop("TypeError", "Multiplication of two bras is not allowed. Use TensorProduct instead.")
  if (isa(a, "IdentityOperator") && isa(b, "Expr")) return(list(b))
  if (isa(a, "Expr") && isa(b, "IdentityOperator")) return(list(a))
  NULL
}

.QPOST_TYPES <- c("StateBase", "QTensorProduct", "Operator")

.qmul_post <- function(out) {
  if (!.qm_state$on || !any(vapply(out, function(a) .qisa(a, .QPOST_TYPES), NA))) return(new_sym("Mul", out))
  e <- new_sym("Mul", out)
  args <- .sympy_args(e)
  result <- list()
  changed <- FALSE
  while (length(args) > 1L) {
    first <- args[[1L]]
    args <- args[-1L]
    second <- args[[1L]]
    tr <- .qtransform_pair(first, second)
    if (is.null(tr)) {
      result[[length(result) + 1L]] <- first
    } else {
      changed <- TRUE
      args[[1L]] <- tr[[length(tr)]]
      result <- c(result, tr[-length(tr)])
    }
  }
  if (!changed) return(e)
  result <- c(result, args)
  if (!length(result)) return(S_One)
  if (length(result) == 1L) return(result[[1L]])
  if (all(vapply(result, .is_comm, NA))) return(.mul(result))
  # Mul._from_args: the factors in the order the transforms left them
  flat <- do.call(c, lapply(result, function(a) if (a$type == "Mul") .qargs(a) else list(a)))
  new_sym("Mul", flat, extra = list(raw = TRUE))
}

.qpow_pre <- function(b, e) {
  if (!.qm_state$on || is.null(.qnode_def(b$type))) return(NULL)
  if (is_number_type(e) && (num_is_zero(e) || num_is_one(e))) return(NULL)
  if (is_nan(e)) return(NULL)
  if (.qhas(b, "_eval_power")) {
    obj <- .qmeth(b, "_eval_power", e)
    if (!is.null(obj)) return(obj)
  }
  # the Pow postprocessors of StateBase and TensorProduct
  if (.qisa(b, "StateBase") && .qkind(b) %in% c("KetKind", "BraKind"))
    .qstop("TypeError", "A bra or ket to a power is invalid, use TensorProduct instead.")
  if (.qisa(b, "QTensorProduct") && isTRUE(ask_fact(e, "integer")) && isTRUE(ask_fact(e, "positive")) &&
      .qkind(b) == "OperatorKind")
    return(.qtp(lapply(b$args, function(a) .pow(a, e))))
  NULL
}

.arith_hooks$mul_post <- .qmul_post
.arith_hooks$pow_pre <- .qpow_pre
