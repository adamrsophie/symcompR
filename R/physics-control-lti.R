# Linear time invariant systems, after sympy/physics/control/lti.py:
# transfer functions and their series, parallel and feedback
# interconnections. MIMO systems and state space models are in
# physics-control-mimo.R and physics-control-statespace.R.

.LTI_MOD <- "sympy.physics.control.lti"

# sympy.utilities.misc.filldedent: textwrap.dedent() and textwrap.fill(),
# which keep runs of spaces inside a line and break after hyphens
.py_filldedent <- function(s, w = 70L) {
  lines <- strsplit(s, "\n", fixed = TRUE)[[1L]]
  blank <- !grepl("\\S", lines)
  ind <- if (any(!blank)) min(nchar(sub("^([ \t]*).*$", "\\1", lines[!blank]))) else 0L
  lines <- ifelse(blank, "", substring(lines, ind + 1L))
  txt <- gsub("^\n+|\n+$", "", paste(lines, collapse = "\n"))
  txt <- gsub("[\t\n\r\f\v]", " ", txt)
  chunks <- regmatches(txt, gregexpr("[^ ]*?[A-Za-z]-(?=[A-Za-z])| +|[^ ]+", txt, perl = TRUE))[[1L]]
  ws <- !grepl("\\S", chunks)
  out <- character()
  i <- 1L; n <- length(chunks)
  while (i <= n) {
    if (length(out) && ws[i]) i <- i + 1L
    cur <- character(); len <- 0L
    while (i <= n && len + nchar(chunks[i]) <= w) {
      cur <- c(cur, chunks[i]); len <- len + nchar(chunks[i]); i <- i + 1L
    }
    if (i <= n && nchar(chunks[i]) > w) {
      k <- max(1L, w - len)
      cur <- c(cur, substr(chunks[i], 1L, k)); chunks[i] <- substring(chunks[i], k + 1L)
    }
    if (length(cur) && !grepl("\\S", cur[length(cur)])) cur <- cur[-length(cur)]
    if (length(cur)) out <- c(out, paste(cur, collapse = ""))
  }
  paste0("\n", paste(out, collapse = "\n"))
}

# repr(type(x)) as in SymPy's error messages
.lti_type_repr <- function(x) {
  if (inherits(x, c("pyobj", "pyenv"))) return(sprintf("<class '%s.%s'>", .LTI_MOD, class(x)[1L]))
  if (inherits(x, "symmatrix")) return("<class 'sympy.matrices.dense.MutableDenseMatrix'>")
  if (is_sym(x)) {
    cls <- .sympy_class(x)
    mod <- switch(cls, Integer = , Rational = , Float = , Zero = , One = , NegativeOne = , Half = ,
                  Pi = , Exp1 = , ImaginaryUnit = , Infinity = , NegativeInfinity = , ComplexInfinity = , `NaN` = "sympy.core.numbers",
                  Symbol = , Dummy = "sympy.core.symbol", Add = "sympy.core.add", Mul = "sympy.core.mul",
                  Pow = "sympy.core.power", ImmutableMatrix = "sympy.matrices.immutable", Tuple = "sympy.core.containers",
                  exp = , log = "sympy.functions.elementary.exponential",
                  sin = , cos = , tan = "sympy.functions.elementary.trigonometric",
                  "sympy.core.basic")
    if (cls == "ImmutableMatrix") cls <- "ImmutableDenseMatrix"
    return(sprintf("<class '%s.%s'>", mod, cls))
  }
  sprintf("<class '%s'>", .py_typename(x))
}

# str(x) for a system, as Python prints it
.lti_str <- function(x) {
  if (inherits(x, "PIDController"))
    return(sprintf("PIDController(%s, %s, %s, %s, %s)", .str(x$kp), .str(x$ki), .str(x$kd), .str(x$tf), .str(x$var)))
  if (inherits(x, c("Series", "Parallel")))
    return(paste0(class(x)[1L], "(", paste(vapply(x$args, .lti_str, ""), collapse = ", "), ")"))
  if (inherits(x, "Feedback"))
    return(sprintf("Feedback(%s, %s, %s)", .lti_str(x$sys1), .lti_str(x$sys2), .str(x$sign)))
  if (is_sym(x)) return(.str(x))
  if (inherits(x, "symmatrix")) return(sstr(x))
  format(x)
}

.lti_sympify <- function(x) {
  if (is_sym(x) || inherits(x, "symmatrix")) return(x)
  if (is.complex(x) && length(x) == 1L) return(Float(Re(x)) + Float(Im(x)) * S_I)
  if (is.null(x)) .physics_error("SympifyError", "SympifyError: None")
  as_sym(x)
}

.lti_is_zero <- function(x) is_sym(x) && x$type %in% c("Integer", "Rational", "Float") && is_zero_num(x)

.lti_eq <- function(a, b) identical(.pyobj_k(a), .pyobj_k(b))

.lti_unsupported <- "Unsupported type for numerator or denominator of\n                            TransferFunction."

.lti_same_var_msg <- function(transform = "Laplace domain or z-domain")
  .py_filldedent(paste0("\n                    All the transfer functions should use the same complex\n                    variable of the ", transform, "."))

# _check_time_compatibility
.lti_check_time <- function(systems) {
  if (!length(systems)) .physics_error("IndexError", "tuple index out of range")
  if (!all(vapply(systems, function(s) inherits(s, "LinearTimeInvariant"), TRUE))) return(invisible(NULL))
  cont <- systems[[1L]]$is_continuous
  st <- systems[[1L]]$sampling_time
  for (sys in systems[-1L]) {
    if (!identical(sys$is_continuous, cont))
      .physics_error("TypeError", paste0("\n                Incompatible systems detected.\n                All systems should be either continuous-time or discrete-time.\n                Found: ",
                                         .lti_str(systems[[1L]]), " and ", .lti_str(sys)))
    if (!.lti_eq(sys$sampling_time, st))
      .physics_error("TypeError", paste0("\n                Incompatible sampling times detected.\n                All systems should have the same sampling time.\n                Found sampling time: ",
                                         sstr(st), " and ", sstr(sys$sampling_time)))
  }
  invisible(NULL)
}

# LinearTimeInvariant._check_args
.lti_check_args <- function(args, clstype) {
  if (!length(args)) .physics_error("ValueError", "At least 1 argument must be passed.")
  if (!all(vapply(args, function(a) inherits(a, clstype), TRUE)))
    .physics_error("TypeError", sprintf("All arguments must be of type <class '%s.%s'>.", .LTI_MOD, clstype))
  vars <- unique(vapply(args, function(a) a$var$key, ""))
  if (length(vars) != 1L)
    .physics_error("ValueError", .py_filldedent(sprintf("\n                All transfer functions should use the same complex variable\n                of the Laplace transform or z-transform. %d different\n                values found.", length(vars))))
  invisible(NULL)
}

.lti_is_ss_obj <- function(x) inherits(x, "LinearTimeInvariant") && isTRUE(x$is_StateSpace_object)

#' Create a transfer function
#'
#' `create_transfer_function(num, den, var, sampling_time)` is a continuous
#' [TransferFunction()] when `sampling_time` is 0 and a
#' [DiscreteTransferFunction()] otherwise; `create_state_space()` does the
#' same for state space models. After `sympy.physics.control.lti`.
#'
#' @param num,den The numerator and denominator.
#' @param var The complex variable of the Laplace transform or z-transform.
#' @param sampling_time The sampling time; 0 for a continuous system.
#' @param A,B,C,D The state, input, output and feedforward matrices.
#' @return A transfer function or a state space model.
#' @examples
#' s <- Symbol("s"); z <- Symbol("z")
#' create_transfer_function(s + 1, s^2 + 2, s)
#' create_transfer_function(z, z + 1, z, 0.1)
#' @export
create_transfer_function <- function(num, den, var, sampling_time = 0L) {
  if (.lti_st_zero(sampling_time)) return(TransferFunction(num, den, var))
  DiscreteTransferFunction(num, den, var, sampling_time)
}

.lti_st_zero <- function(st) (is.numeric(st) && !is.object(st) && st == 0) || .lti_is_zero(st)

# TransferFunctionBase.__new__
.tfb_new <- function(cls, num, den, var, extra_args = list(), extra = list()) {
  num <- .lti_sympify(num); den <- .lti_sympify(den)
  if (!(is_sym(var) && var$type %in% c("Symbol", "Dummy"))) .physics_error("TypeError", "Variable input must be a Symbol.")
  if (.lti_is_zero(den)) .physics_error("ValueError", "TransferFunction cannot have a zero denominator.")
  accepted <- function(e) {
    if (!is_sym(e)) return(FALSE)
    if (inherits(e, "pyobj")) return(inherits(e, c("TransferFunctionBase", "Series", "Parallel")) && .has_type(e, c("Symbol", "Dummy")))
    .has_type(e, c("Symbol", "Dummy")) || .expr_is_number(e)
  }
  if (!accepted(num) || !accepted(den)) .physics_error("TypeError", .lti_unsupported)
  .pyobj_new(cls, c(list(num, den, var), extra_args), extra = c(list(is_StateSpace_object = FALSE), extra))
}

.tf_from_rational_expression <- function(make, cls) function(expr, var = NULL, ...) {
  if (cls == "TransferFunctionBase") .tfb_abstract()
  expr <- .lti_sympify(expr)
  if (is.null(var)) {
    fs <- free_symbols(expr)
    if (length(fs) == 1L) {
      var <- fs[[1L]]
    } else if (length(fs) == 0L) {
      .physics_error("ValueError", .py_filldedent("\n                    Positional argument `var` not found in the\n                    TransferFunction defined. Specify it manually."))
    } else {
      .physics_error("ValueError", .py_filldedent(sprintf("\n                    Conflicting values found for positional argument `var` ({%s}).\n                    Specify it manually.",
                                                      paste(vapply(fs, .str, ""), collapse = ", "))))
    }
  }
  nd <- .as_numer_denom(expr)
  if (.lti_is_zero(nd[[2L]]) || .has1(nd[[1L]], zoo))
    .physics_error("ZeroDivisionError", "TransferFunction cannot have a zero denominator.")
  make(nd[[1L]], nd[[2L]], var, ...)
}

.tf_from_coeff_lists <- function(make, cls) function(num_list, den_list, var, ...) {
  if (cls == "TransferFunctionBase") .tfb_abstract()
  num_list <- rev(as.list(unclass(num_list))); den_list <- rev(as.list(unclass(den_list)))
  sumpow <- function(l) Reduce(function(acc, i) acc + .lti_sympify(l[[i]]) * var^(i - 1L), seq_along(l), S_Zero)
  .num <- sumpow(num_list); .den <- sumpow(den_list)
  if (.lti_is_zero(.den)) .physics_error("ZeroDivisionError", "TransferFunction cannot have a zero denominator.")
  make(.num, .den, var, ...)
}

.tf_from_zpk <- function(make, cls) function(zeros, poles, gain, var, ...) {
  if (cls == "TransferFunctionBase") .tfb_abstract()
  num_poly <- S_One; den_poly <- S_One
  for (z in as.list(unclass(zeros))) num_poly <- num_poly * (var - .lti_sympify(z))
  for (p in as.list(unclass(poles))) den_poly <- den_poly * (var - .lti_sympify(p))
  make(.lti_sympify(gain) * num_poly, den_poly, var, ...)
}

.tfb_abstract <- function() .physics_error("NotImplementedError", "\n                The TransferFunctionBase class is not meant to be used directly.\n                ")

# _roots
.lti_roots <- function(poly_expr, var) {
  r <- roots(poly_expr, var, multiple = TRUE)
  if (is.null(r)) r <- list()
  n <- as.integer(as.character(degree(poly_expr, var)))
  if (length(r) != n) {
    # rootof() works over ZZ and QQ only; its roots are those of all_roots()
    coeffs <- all_coeffs(Poly(poly_expr, var))
    kind <- .rc_domain(coeffs)
    if (!kind %in% c("ZZ", "QQ", "RR")) {
      gens <- sort(unique(unlist(lapply(coeffs, function(c) vapply(free_symbols(c), function(z) z$name, "")))))
      ring <- if (all(vapply(coeffs, function(c) all(vapply(.coeff_terms(c), function(t) t$type == "Integer", TRUE)), TRUE))) "ZZ" else "QQ"
      dom <- switch(kind, poly = sprintf("%s[%s]", ring, paste(gens, collapse = ",")),
                    frac = sprintf("%s(%s)", ring, paste(gens, collapse = ",")), "EX")
      .physics_error("NotImplementedError", sprintf("CRootOf is not supported over %s", dom))
    }
    r <- all_roots(poly_expr, var)
  }
  r
}

# Mul.flatten() evaluates the unevaluated products and powers of to_expr()
.lti_eval <- function(e) {
  if (!is_sym(e) || !e$type %in% c("Mul", "Pow")) return(e)
  sym_rebuild(e, lapply(e$args, .lti_eval))
}

# the numeric coefficients of the terms of an expression
.coeff_terms <- function(e) {
  terms <- if (e$type == "Add") e$args else list(e)
  lapply(terms, function(t) if (t$type == "Mul" && t$args[[1L]]$type %in% c("Integer", "Rational", "Float")) t$args[[1L]]
         else if (t$type %in% c("Integer", "Rational", "Float")) t else S_One)
}

.lti_degree <- function(e, var) {
  d <- degree(e, var)
  if (is_sym(d)) d else as_sym(d)
}

.lti_to_expr_tf <- function(self) {
  if (!identical(self$num$key, S_One$key))
    return(new_sym("Mul", list(self$num, new_sym("Pow", list(self$den, S_NegativeOne)))))
  new_sym("Pow", list(self$den, S_NegativeOne))
}

# cancel(f, expand = False): the numerator and the denominator are not
# expanded, so their terms' factors become the polynomial generators and the
# factored forms survive; Poly.cancel(include = False) then gives c*(P/Q).
# Floats, I, Piecewise and noncommutative parts go through cancel().
.lti_cancel_noexpand <- function(f) {
  f <- signsimp(f)
  if (is_number_type(f) || is_relational(f) || !.is_expr_node(f)) return(f)
  f <- factor_terms(f, radical = TRUE)
  nd <- .as_numer_denom(f)
  if (.has_type(f, "Piecewise") || !.is_comm(f) || .has_type(f, "Float")) return(cancel(f))
  gens <- .find_gens(nd)
  if (!length(gens)) return(expand(f))
  if (any(vapply(gens, function(g) identical(g$key, S_I$key), TRUE))) return(cancel(f))
  # composite generators stand in as dummies while the terms are collected
  dums <- gens; rule <- list()
  for (j in which(!vapply(gens, is_symbol_like, TRUE))) {
    dums[[j]] <- Dummy()
    rule[[length(rule) + 1L]] <- list(gens[[j]], dums[[j]])
  }
  polys <- lapply(nd, function(z) expr_to_qpoly(if (length(rule)) xreplace(z, rule) else z, dums, strict = TRUE))
  if (any(vapply(polys, is.null, TRUE))) return(cancel(f))
  p <- polys[[1L]]; q <- polys[[2L]]
  if (qp_is_zero(p)) return(S_Zero)
  # dmp_cancel: clear denominators, divide by the gcd over ZZ, then move the
  # signs of the leading coefficients into c
  cq <- .bq_lcm_den(p); cp <- .bq_lcm_den(q)
  p <- qp_scale(p, gmp::as.bigq(cq)); q <- qp_scale(q, gmp::as.bigq(cp))
  h <- qp_gcd(p, q)
  if (!qp_is_const(h)) { h <- qp_primitive(h)$poly; p <- qp_exquo(p, h); q <- qp_exquo(q, h) }
  k <- gmp::gcd.bigz(.bq_int_content(p), .bq_int_content(q))
  if (k != 1) { p <- qp_scale(p, gmp::as.bigq(1, k)); q <- qp_scale(q, gmp::as.bigq(1, k)) }
  G <- gmp::gcd.bigz(cp, cq)
  c <- gmp::as.bigq(cp / G) / gmp::as.bigq(cq / G)
  pn <- qp_lc(p) < 0; qn <- qp_lc(q) < 0
  if (pn && qn) { p <- qp_neg(p); q <- qp_neg(q) }
  else if (pn) { c <- -c; p <- qp_neg(p) }
  else if (qn) { c <- -c; q <- qp_neg(q) }
  .mul(list(make_rational(c), .div(qpoly_to_expr(p, gens), qpoly_to_expr(q, gens))))
}

.tfb_methods <- list(
  to_expr = .lti_to_expr_tf,
  expand = function(self, ...) create_transfer_function(expand(self$num), expand(self$den), self$var, self$sampling_time),
  simplify = function(self, ...) {
    nd <- .as_numer_denom(.lti_cancel_noexpand(new_sym("Mul", list(self$num, 1 / self$den))))
    create_transfer_function(nd[[1L]], nd[[2L]], self$var, self$sampling_time)
  },
  evalf = function(self, n = 15L, ...) create_transfer_function(N(self$num, n), N(self$den, n), self$var, self$sampling_time),
  n = function(self, n = 15L, ...) create_transfer_function(N(self$num, n), N(self$den, n), self$var, self$sampling_time),
  poles = function(self) .lti_roots(self$den, self$var),
  zeros = function(self) .lti_roots(self$num, self$var),
  eval_frequency = function(self, other) {
    arg_num <- subs(self$num, self$var, other)
    arg_den <- subs(self$den, self$var, other)
    expand(arg_num * (S_One / arg_den))
  },
  is_stable = function(self, cancel_poles_zeros = FALSE) {
    tf <- self$to_standard_form(cancel_poles_zeros)
    conditions <- tf$get_asymptotic_stability_conditions(cancel_poles_zeros = FALSE)
    # reduce_inequalities() raises NotImplementedError for an inequality in
    # more than one symbol
    if (any(vapply(conditions, function(c) length(free_symbols(as_sym(c))) > 1L, TRUE))) return(NULL)
    out <- tryCatch(reduce_inequalities(conditions), NotImplementedError = function(e) NULL,
                    error = function(e) if (grepl("NotImplemented", conditionMessage(e))) NULL else stop(e))
    if (is.null(out)) return(NULL)
    if (is.logical(out)) return(out)
    if (is_sym(out) && out$type == "BooleanAtom") return(out$name == "true")
    NULL
  },
  to_standard_form = function(self, cancel_poles_zeros = FALSE) {
    tf <- self$expand()
    num <- tf$num; den <- tf$den
    if (isTRUE(cancel_poles_zeros)) {
      nd <- .as_numer_denom(cancel(tf$num / tf$den))
      num <- nd[[1L]]; den <- nd[[2L]]
    }
    create_transfer_function(collect(num, self$var), collect(den, self$var), self$var, self$sampling_time)
  },
  `_StateSpace_matrices_equivalent` = function(self) .tf_ss_matrices(self),
  doit = function(self, ...) .pyobj_rebuild_doit(self))

local({
  arg <- function(i) function(self) self$args[[i]]
  deg <- function(op) function(self) op(.lti_degree(self$num, self$var), .lti_degree(self$den, self$var))
  .pyclass_def("LinearTimeInvariant",
    props = list(is_continuous = function(self) self$extra$is_continuous,
                 is_StateSpace_object = function(self) isTRUE(self$extra$is_StateSpace_object)))
  .pyclass_def("SISOLinearTimeInvariant", parent = "LinearTimeInvariant",
    props = list(num_inputs = function(self) 1L, num_outputs = function(self) 1L, is_SISO = function(self) TRUE))
  .pyclass_def("MIMOLinearTimeInvariant", parent = "LinearTimeInvariant",
    props = list(is_SISO = function(self) FALSE))
  .pyclass_def("TransferFunctionBase", parent = "SISOLinearTimeInvariant",
    props = list(num = arg(1L), den = arg(2L), var = arg(3L),
                 is_proper = function(self) isTRUE(.lti_degree(self$num, self$var) <= .lti_degree(self$den, self$var)),
                 is_strictly_proper = function(self) isTRUE(.lti_degree(self$num, self$var) < .lti_degree(self$den, self$var)),
                 is_biproper = function(self) .lti_eq(.lti_degree(self$num, self$var), .lti_degree(self$den, self$var))),
    methods = .tfb_methods,
    ops = list(`__add__` = function(self, other) .tfb_add(self, other),
               `__radd__` = function(self, other) { .lti_check_time(list(self, other)); .py_binop("+", self, other) },
               `__sub__` = function(self, other) .tfb_sub(self, other),
               `__rsub__` = function(self, other) { .lti_check_time(list(self, other)); .py_binop("+", .tfb_neg(self), other) },
               `__mul__` = function(self, other) .tfb_mul(self, other),
               `__rmul__` = function(self, other) .tfb_mul(self, other),
               `__truediv__` = function(self, other) .tfb_div(self, other),
               `__rtruediv__` = function(self, other) .tfb_div(self, other),
               `__pow__` = function(self, p) .tfb_pow(self, p),
               `__neg__` = function(self) .tfb_neg(self)),
    statics = list(
      from_rational_expression = .tf_from_rational_expression(NULL, "TransferFunctionBase"),
      from_coeff_lists = .tf_from_coeff_lists(NULL, "TransferFunctionBase"),
      from_zpk = .tf_from_zpk(NULL, "TransferFunctionBase")))
  .pyclass_def("TransferFunction", parent = "TransferFunctionBase",
    ctor = function(num, den, var) TransferFunction(num, den, var),
    props = list(sampling_time = function(self) S_Zero, is_continuous = function(self) TRUE),
    methods = list(
      dc_gain = function(self) limit(new_sym("Mul", list(self$num, new_sym("Pow", list(self$den, S_NegativeOne)))), self$var, S_Zero),
      get_asymptotic_stability_conditions = function(self, cancel_poles_zeros = FALSE, fast = FALSE) {
        sf <- self$to_standard_form(cancel_poles_zeros)
        lapply(.lti_hurwitz(sf$den, self$var, fast), function(c) c > 0)
      },
      `_eval_rewrite_as_StateSpace` = function(self) {
        m <- self$`_StateSpace_matrices_equivalent`()
        StateSpace(m[[1L]], m[[2L]], m[[3L]], m[[4L]])
      },
      `_eval_rewrite_as_DiscreteStateSpace` = function(self)
        .physics_error("TypeError", "\n            The continuous transfer function model cannot be rewritten as a\n            discrete-time state space model.\n            ")),
    statics = list(
      from_rational_expression = .tf_from_rational_expression(function(n, d, v) TransferFunction(n, d, v), "TransferFunction"),
      from_coeff_lists = .tf_from_coeff_lists(function(n, d, v) TransferFunction(n, d, v), "TransferFunction"),
      from_zpk = .tf_from_zpk(function(n, d, v) TransferFunction(n, d, v), "TransferFunction")))
  .pyclass_def("DiscreteTransferFunction", parent = "TransferFunctionBase",
    ctor = function(num, den, var, sampling_time) DiscreteTransferFunction(num, den, var, sampling_time),
    props = list(sampling_time = arg(4L), is_continuous = function(self) FALSE),
    methods = list(
      dc_gain = function(self) limit(new_sym("Mul", list(self$num, new_sym("Pow", list(self$den, S_NegativeOne)))), self$var, S_One),
      get_asymptotic_stability_conditions = function(self, cancel_poles_zeros = FALSE, fast = FALSE) {
        sf <- self$to_standard_form(cancel_poles_zeros)
        lapply(.lti_schur(sf$den, self$var, fast), function(c) c > 0)
      },
      `_eval_rewrite_as_DiscreteStateSpace` = function(self) {
        m <- self$`_StateSpace_matrices_equivalent`()
        DiscreteStateSpace(m[[1L]], m[[2L]], m[[3L]], m[[4L]], self$sampling_time)
      },
      `_eval_rewrite_as_StateSpace` = function(self)
        .physics_error("TypeError", "\n            The discrete transfer function model cannot be rewritten as a\n            continuous-time state space model.\n            ")),
    statics = list(
      from_rational_expression = function(expr, var = NULL, sampling_time = 1L)
        .tf_from_rational_expression(function(n, d, v) DiscreteTransferFunction(n, d, v, sampling_time), "DiscreteTransferFunction")(expr, var),
      from_coeff_lists = function(num_list, den_list, var, sampling_time = 1L)
        .tf_from_coeff_lists(function(n, d, v) DiscreteTransferFunction(n, d, v, sampling_time), "DiscreteTransferFunction")(num_list, den_list, var),
      from_zpk = function(zeros, poles, gain, var, sampling_time = 1L)
        .tf_from_zpk(function(n, d, v) DiscreteTransferFunction(n, d, v, sampling_time), "DiscreteTransferFunction")(zeros, poles, gain, var),
      from_gbt = function(cont_tf, sampling_time, alpha, var) {
        r <- gbt(cont_tf, sampling_time, alpha)
        DiscreteTransferFunction$from_coeff_lists(r[[1L]], r[[2L]], var, sampling_time)
      },
      from_bilinear = function(cont_tf, sampling_time, var) {
        r <- bilinear(cont_tf, sampling_time)
        DiscreteTransferFunction$from_coeff_lists(r[[1L]], r[[2L]], var, sampling_time)
      },
      from_forward_diff = function(cont_tf, sampling_time, var) {
        r <- forward_diff(cont_tf, sampling_time)
        DiscreteTransferFunction$from_coeff_lists(r[[1L]], r[[2L]], var, sampling_time)
      },
      from_backward_diff = function(cont_tf, sampling_time, var) {
        r <- backward_diff(cont_tf, sampling_time)
        DiscreteTransferFunction$from_coeff_lists(r[[1L]], r[[2L]], var, sampling_time)
      }))
  .pyclass_def("PIDController", parent = "TransferFunction",
    ctor = function(num, den, var) TransferFunction(num, den, var),
    props = list(kp = function(self) self$extra$kp, ki = function(self) self$extra$ki,
                 kd = function(self) self$extra$kd, tf = function(self) self$extra$tf),
    methods = list(doit = function(self, ...) TransferFunction(self$num, self$den, self$var)))
})

# the subs hook: a transfer function keeps its variable
local({
  h <- function(e, old, new) {
    if (identical(old$key, e$var$key)) return(e)
    create_transfer_function(subs(e$num, old, new), subs(e$den, old, new), e$var, e$sampling_time)
  }
  for (t in c("TransferFunction", "DiscreteTransferFunction", "PIDController")) register_subs_hook(t, h)
})

#' Transfer functions
#'
#' `TransferFunction(num, den, var)` is the continuous transfer function
#' `num/den` in the complex variable `var` of the Laplace transform, and
#' `DiscreteTransferFunction(num, den, var, sampling_time)` the discrete one
#' in the variable of the z-transform. `PIDController(kp, ki, kd, tf, var)`
#' is the transfer function of a PID controller with a derivative filter
#' time constant `tf`. The class methods
#' `TransferFunction$from_rational_expression(expr, var)`,
#' `$from_coeff_lists(num_list, den_list, var)` and
#' `$from_zpk(zeros, poles, gain, var)` build transfer functions; the
#' discrete class also has `$from_gbt()`, `$from_bilinear()`,
#' `$from_forward_diff()` and `$from_backward_diff()`. A transfer function has
#' the properties `$num`, `$den`, `$var`, `$sampling_time`, `$is_proper`,
#' `$is_strictly_proper` and `$is_biproper` and the methods `$poles()`,
#' `$zeros()`, `$dc_gain()`, `$eval_frequency(other)`, `$is_stable()`,
#' `$to_expr()`, `$to_standard_form()`, `$simplify()`, `$expand()` and
#' `$get_asymptotic_stability_conditions()`. The operators `+`, `-`, `*`,
#' `/` and `^` build `Series()`, `Parallel()` and `Feedback()`
#' interconnections. After `sympy.physics.control.lti`.
#'
#' @param num,den The numerator and denominator.
#' @param var The complex variable.
#' @param sampling_time The sampling time.
#' @param kp,ki,kd The proportional, integral and derivative gains.
#' @param tf The derivative filter time constant.
#' @return A `TransferFunction`, `DiscreteTransferFunction` or
#'   `PIDController` object.
#' @examples
#' s <- Symbol("s")
#' G <- TransferFunction(s + 1, s^2 + 2, s)
#' G$poles()
#' G * TransferFunction(1, s, s)
#' TransferFunction$from_rational_expression((s + 1) / (s^2 - 1), s)
#' @usage TransferFunction(num, den, var)
#' @export
TransferFunction <- .pyclass_fn(function(num, den, var) .tfb_new("TransferFunction", num, den, var, extra = list(is_continuous = TRUE)),
                                "TransferFunction")

#' @rdname TransferFunction
#' @usage DiscreteTransferFunction(num, den, var, sampling_time = 1)
#' @export
DiscreteTransferFunction <- .pyclass_fn(function(num, den, var, sampling_time = 1L) {
  if (.lti_st_zero(sampling_time))
    .physics_error("ValueError", .py_filldedent("\n                The sampling time cannot be zero.\n                If you want to create a continuous transfer function,\n                use the TransferFunction class instead."))
  sampling_time <- .lti_sympify(sampling_time)
  .tfb_new("DiscreteTransferFunction", num, den, var, list(sampling_time), extra = list(is_continuous = FALSE))
}, "DiscreteTransferFunction")

#' @rdname TransferFunction
#' @usage TransferFunctionBase
#' @format NULL
#' @export
TransferFunctionBase <- .pyclass_fn(function(...) .tfb_abstract(), "TransferFunctionBase")

#' @rdname TransferFunction
#' @export
PIDController <- .pyclass_fn(function(kp = Symbol("kp"), ki = Symbol("ki"), kd = Symbol("kd"), tf = 0L, var = Symbol("s")) {
  kp0 <- kp; ki0 <- ki; kd0 <- kd; tf0 <- tf
  kp <- .lti_sympify(kp); ki <- .lti_sympify(ki); kd <- .lti_sympify(kd); tf <- .lti_sympify(tf)
  num <- kp * tf * var^2 + kp * var + ki * tf * var + ki + kd * var^2
  den <- tf * var^2 + var
  .tfb_new("PIDController", num, den, var,
           extra = list(is_continuous = TRUE, kp = .lti_sympify(kp0), ki = .lti_sympify(ki0), kd = .lti_sympify(kd0), tf = .lti_sympify(tf0)))
}, "PIDController")

#' @export
print.PIDController <- function(x, ...) {
  cat(.lti_str(x), "\n", sep = "")
  invisible(x)
}

# Basic.doit: the class rebuilt from the evaluated arguments
.pyobj_rebuild_doit <- function(self) {
  args <- lapply(self$args, function(a) if (inherits(a, "pyobj")) a$doit() else if (is_sym(a)) doit(a) else a)
  do.call(.PYCLASSES[[class(self)[1L]]]$ctor, args)
}

# ---------------------------------------------------------------------------
# TransferFunctionBase operators

.tfb_same_class <- function(self, other) identical(class(other)[1L], class(self)[1L])

.tfb_add <- function(self, other) {
  .lti_check_time(list(self, other))
  if (.lti_is_ss_obj(other)) return(Parallel(self, other))
  if (.tfb_same_class(self, other) || inherits(other, c("Series", "Feedback"))) {
    if (!.lti_eq(self$var, other$var)) .physics_error("ValueError", .lti_same_var_msg())
    return(Parallel(self, other))
  }
  if (inherits(other, "Parallel")) {
    if (!.lti_eq(self$var, other$var)) .physics_error("ValueError", .lti_same_var_msg())
    return(do.call(Parallel, c(list(self), other$args)))
  }
  .physics_error("ValueError", sprintf("%s cannot be added with %s.", .lti_type_repr(self), .lti_type_repr(other)))
}

.tfb_sub <- function(self, other) {
  .lti_check_time(list(self, other))
  if (.lti_is_ss_obj(other)) return(Parallel(self, -other))
  if (.tfb_same_class(self, other) || inherits(other, "Series")) {
    if (!.lti_eq(self$var, other$var)) .physics_error("ValueError", .lti_same_var_msg())
    return(Parallel(self, -other))
  }
  if (inherits(other, "Parallel")) {
    if (!.lti_eq(self$var, other$var)) .physics_error("ValueError", .lti_same_var_msg("Laplace transform"))
    return(do.call(Parallel, c(list(self), lapply(other$args, function(a) -a))))
  }
  .physics_error("ValueError", sprintf("%s cannot be subtracted from a %s.", .lti_type_repr(self), .lti_type_repr(other)))
}

.tfb_mul <- function(self, other) {
  .lti_check_time(list(self, other))
  if (.lti_is_ss_obj(other)) return(Series(self, other))
  if (.tfb_same_class(self, other) || inherits(other, c("Parallel", "Feedback"))) {
    if (!.lti_eq(self$var, other$var))
      .physics_error("ValueError", .py_filldedent("\n                    All the transfer functions should use the same complex variable\n                    of the Laplace domain or z-domain."))
    return(Series(self, other))
  }
  if (inherits(other, "Series")) {
    if (!.lti_eq(self$var, other$var))
      .physics_error("ValueError", .py_filldedent("\n                    All the transfer functions should use the same complex variable\n                    of the Laplace domain or z-domain."))
    return(do.call(Series, c(list(self), other$args)))
  }
  .physics_error("ValueError", sprintf("%s cannot be multiplied with %s.", .lti_type_repr(self), .lti_type_repr(other)))
}

.tfb_div <- function(self, other) {
  .lti_check_time(list(self, other))
  if (.tfb_same_class(self, other)) {
    if (!.lti_eq(self$var, other$var)) .physics_error("ValueError", .lti_same_var_msg())
    return(Series(self, create_transfer_function(other$den, other$num, self$var, self$sampling_time)))
  }
  if (inherits(other, "Parallel") && length(other$args) == 2L && .tfb_same_class(self, other$args[[1L]]) &&
      (inherits(other$args[[2L]], "Series") || .tfb_same_class(self, other$args[[2L]]))) {
    if (!.lti_eq(self$var, other$var))
      .physics_error("ValueError", .py_filldedent("\n                    Both TransferFunction and Parallel should use the\n                    same complex variable of the Laplace domain or z-domain."))
    a2 <- other$args[[2L]]
    if (.lti_eq(a2, self)) return(Feedback(self, other$args[[1L]]))
    if (!inherits(a2, "Series")) return(Feedback(self, a2))
    lst <- a2$args
    keys <- vapply(lst, .pyobj_k, "")
    i <- match(.pyobj_k(self), keys)
    if (is.na(i)) return(Feedback(self, do.call(Series, lst)))
    lst <- lst[-i]
    if (length(lst) == 1L) return(Feedback(self, lst[[1L]]))
    return(Feedback(self, do.call(Series, lst)))
  }
  .physics_error("ValueError", sprintf("%s cannot be divided by %s.", .lti_type_repr(self), .lti_type_repr(other)))
}

.tfb_pow <- function(self, p) {
  p <- .lti_sympify(p)
  if (!(is_sym(p) && p$type == "Integer")) .physics_error("ValueError", "Exponent must be an integer.")
  if (identical(p$key, S_Zero$key)) return(create_transfer_function(1L, 1L, self$var, self$sampling_time))
  if (isTRUE(p > 0)) {
    num_ <- self$num^p; den_ <- self$den^p
  } else {
    p <- abs(p)
    num_ <- self$den^p; den_ <- self$num^p
  }
  create_transfer_function(num_, den_, self$var, self$sampling_time)
}

.tfb_neg <- function(self) create_transfer_function(-self$num, self$den, self$var, self$sampling_time)

# TransferFunctionBase._StateSpace_matrices_equivalent
.tf_ss_matrices <- function(self) {
  if (!self$is_proper) .physics_error("ValueError", "Transfer Function must be proper.")
  num_coeffs <- all_coeffs(Poly(self$num, self$var))
  den_coeffs <- all_coeffs(Poly(self$den, self$var))
  n <- length(den_coeffs) - 1L
  M1 <- function(v) Matrix(list(list(v)))
  if (n == 0L) return(list(M1(S_Zero), M1(S_Zero), M1(S_Zero), M1(num_coeffs[[1L]] / den_coeffs[[1L]])))
  if (.lti_eq(self$num, self$den)) return(list(M1(S_Zero), M1(S_Zero), M1(S_Zero), M1(S_One)))
  diff <- n - (length(num_coeffs) - 1L)
  num_coeffs <- c(rep(list(S_Zero), diff), num_coeffs)
  a <- den_coeffs[-1L]
  a_mat <- Matrix(list(lapply(rev(a), function(cf) (-1L) * cf / den_coeffs[[1L]])))
  A <- rbind(cbind(zeros(n - 1L, 1L), eye(n - 1L)), a_mat)
  B <- zeros(n, 1L)
  B$data[[n]] <- S_One
  C <- list()
  i <- n
  while (i > 0L) {
    C[[length(C) + 1L]] <- num_coeffs[[i + 1L]] / den_coeffs[[1L]] - den_coeffs[[i + 1L]] * num_coeffs[[1L]] / (den_coeffs[[1L]]^2)
    i <- i - 1L
  }
  list(A, B, Matrix(list(C)), M1(num_coeffs[[1L]] / den_coeffs[[1L]]))
}

# ---------------------------------------------------------------------------
# Series, Parallel and Feedback

.lti_flatten <- function(args, cls) {
  out <- list()
  for (a in args) if (inherits(a, cls)) out <- c(out, a$args) else out[[length(out) + 1L]] <- a
  out
}

.lti_any_ss <- function(args) any(vapply(args, function(a) inherits(a, "StateSpaceBase") || .lti_is_ss_obj(a), TRUE))

.lti_io_compatible <- function(args) {
  if (length(args) > 1L) for (i in 2:length(args)) if (!identical(args[[i]]$num_inputs, args[[i - 1L]]$num_outputs)) return(FALSE)
  TRUE
}

.lti_ss_class <- function(continuous) if (continuous) "StateSpace" else "DiscreteStateSpace"
.lti_tf_class <- function(continuous) if (continuous) "TransferFunction" else "DiscreteTransferFunction"


#' Series, parallel and feedback interconnections
#'
#' `Series(...)` is the series and `Parallel(...)` the parallel connection
#' of SISO systems; `Feedback(sys1, sys2, sign)` is the closed loop with
#' `sys1` in the feed forward path and `sys2` in the feedback path (`sign`
#' -1 for negative feedback). `$doit()` evaluates an interconnection to a
#' transfer function (or a state space model), `$to_expr()` gives it as an
#' expression and `$rewrite("TransferFunction")` rewrites it. After
#' `sympy.physics.control.lti`.
#'
#' @param ... The systems.
#' @param evaluate Whether to evaluate the interconnection.
#' @param sys1,sys2 The systems of the feed forward and feedback paths.
#' @param sign The sign of the feedback.
#' @return A `Series`, `Parallel` or `Feedback` object.
#' @examples
#' s <- Symbol("s")
#' G1 <- TransferFunction(1, s + 1, s); G2 <- TransferFunction(s, s + 2, s)
#' Series(G1, G2)$doit()
#' Parallel(G1, G2)$doit()
#' Feedback(G1, G2)$doit()
#' @rdname control-interconnections
#' @export
Series <- function(..., evaluate = FALSE) {
  args <- .lti_flatten(list(...), "Series")
  .lti_check_time(args)
  ss <- FALSE
  if (length(args) && .lti_any_ss(args)) {
    if (any(vapply(args, function(a) !isTRUE(a$is_SISO), TRUE)))
      .physics_error("ValueError", .py_filldedent("To use Series connection for MIMO systems use\n                               MIMOSeries instead."))
    if (!.lti_io_compatible(args))
      .physics_error("ValueError", .py_filldedent("\n                        Systems with incompatible inputs and outputs\n                        cannot be connected in Series.\n                        "))
    ss <- TRUE
  } else {
    .lti_check_args(args, "SISOLinearTimeInvariant")
  }
  obj <- .pyobj_new("Series", args, extra = list(ss = ss, is_StateSpace_object = ss, is_continuous = args[[1L]]$is_continuous))
  if (isTRUE(evaluate)) obj$doit() else obj
}

#' @rdname control-interconnections
#' @export
Parallel <- function(..., evaluate = FALSE) {
  args <- .lti_flatten(list(...), "Parallel")
  ss <- FALSE
  if (length(args) && .lti_any_ss(args)) {
    if (any(vapply(args, function(a) !isTRUE(a$is_SISO), TRUE)))
      .physics_error("ValueError", .py_filldedent("\n                    To use Parallel connection for MIMO systems use\n                    MIMOParallel instead."))
    ss <- TRUE
  } else {
    .lti_check_args(args, "SISOLinearTimeInvariant")
  }
  .lti_check_time(args)
  obj <- .pyobj_new("Parallel", args, extra = list(ss = ss, is_StateSpace_object = ss, is_continuous = args[[1L]]$is_continuous))
  if (isTRUE(evaluate)) obj$doit() else obj
}

.series_doit <- function(self, ...) {
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
  ds <- lapply(self$args, function(a) a$doit())
  res_num <- .mul(lapply(ds, function(d) d$num))
  res_den <- .mul(lapply(ds, function(d) d$den))
  create_transfer_function(res_num, res_den, self$var, self$args[[1L]]$sampling_time)
}

.parallel_doit <- function(self, ...) {
  if (isTRUE(self$extra$ss)) {
    ss_class <- .lti_ss_class(self$is_continuous)
    res <- self$args[[1L]]$doit()
    if (!inherits(res, ss_class)) res <- res$rewrite(ss_class)
    for (arg in self$args[-1L]) {
      if (!inherits(arg, ss_class)) arg <- arg$doit()$rewrite(ss_class)
      res <- res + arg
    }
    return(res)
  }
  terms <- lapply(self$args, function(a) a$doit()$to_expr())
  nd <- .as_numer_denom(.add(terms))
  create_transfer_function(nd[[1L]], nd[[2L]], self$var, self$args[[1L]]$sampling_time)
}

.lti_rewrite_tf <- function(self, continuous_target) {
  name <- if (continuous_target) "TransferFunction" else "DiscreteTransferFunction"
  if (!identical(self$is_continuous, continuous_target)) {
    if (continuous_target)
      .physics_error("TypeError", sprintf("\n                    Cannot rewrite a discrete-time %s object as a\n                    TransferFunction.", class(self)[1L]))
    .physics_error("TypeError", sprintf("\n                    Cannot rewrite a continuous-time %s object as a\n                    DiscreteTransferFunction.", class(self)[1L]))
  }
  if (isTRUE(self$is_StateSpace_object)) return(self$doit()$rewrite(name)[[1L]][[1L]])
  self$doit()
}

local({
  common_props <- list(
    var = function(self) self$args[[1L]]$var,
    sampling_time = function(self) self$args[[1L]]$sampling_time,
    is_proper = function(self) self$doit()$is_proper,
    is_strictly_proper = function(self) self$doit()$is_strictly_proper,
    is_biproper = function(self) self$doit()$is_biproper,
    is_StateSpace_object = function(self) isTRUE(self$extra$ss))
  rw <- list(
    `_eval_rewrite_as_TransferFunction` = function(self) .lti_rewrite_tf(self, TRUE),
    `_eval_rewrite_as_DiscreteTransferFunction` = function(self) .lti_rewrite_tf(self, FALSE))
  siso_other <- function(f) function(self, other) {
    .lti_check_time(list(self, other))
    if (!inherits(other, "SISOLinearTimeInvariant")) return(NULL)
    f(self, other)
  }
  .pyclass_def("Series", parent = "SISOLinearTimeInvariant",
    ctor = function(...) Series(...),
    props = common_props,
    methods = c(rw, list(
      doit = .series_doit,
      to_expr = function(self) new_sym("Mul", lapply(self$args, function(a) a$to_expr())))),
    ops = list(
      `__add__` = siso_other(function(self, other) {
        if (inherits(other, "Parallel")) return(do.call(Parallel, c(list(self), other$args)))
        Parallel(self, other)
      }),
      `__radd__` = siso_other(function(self, other) {
        if (inherits(other, "Parallel")) return(do.call(Parallel, c(list(self), other$args)))
        Parallel(self, other)
      }),
      `__sub__` = siso_other(function(self, other) self + (-other)),
      `__rsub__` = function(self, other) { .lti_check_time(list(self, other)); .py_binop("+", -self, other) },
      `__mul__` = siso_other(function(self, other) do.call(Series, c(self$args, list(other)))),
      `__truediv__` = function(self, other) .series_div(self, other),
      `__neg__` = function(self) Series(create_transfer_function(-1L, 1L, self$var, self$sampling_time), self)))
  .pyclass_def("Parallel", parent = "SISOLinearTimeInvariant",
    ctor = function(...) Parallel(...),
    props = common_props,
    methods = c(rw, list(
      doit = .parallel_doit,
      to_expr = function(self) new_sym("Add", lapply(self$args, function(a) a$to_expr())))),
    ops = list(
      `__add__` = siso_other(function(self, other) do.call(Parallel, c(self$args, list(other)))),
      `__radd__` = siso_other(function(self, other) do.call(Parallel, c(self$args, list(other)))),
      `__sub__` = siso_other(function(self, other) self + (-other)),
      `__rsub__` = function(self, other) { .lti_check_time(list(self, other)); .py_binop("+", -self, other) },
      `__mul__` = siso_other(function(self, other) {
        if (inherits(other, "Series")) return(do.call(Series, c(list(self), other$args)))
        Series(self, other)
      }),
      `__neg__` = function(self) Series(create_transfer_function(-1L, 1L, self$var, self$sampling_time), self)))
})

.series_div <- function(self, other) {
  .lti_check_time(list(self, other))
  tf_class <- .lti_tf_class(self$is_continuous)
  if (inherits(other, "TransferFunctionBase"))
    return(do.call(Series, c(self$args, list(create_transfer_function(other$den, other$num, other$var, other$sampling_time)))))
  if (inherits(other, "Series")) {
    tf_self <- self$rewrite(tf_class)
    tf_other <- other$rewrite(tf_class)
    return(tf_self / tf_other)
  }
  if (inherits(other, "Parallel") && length(other$args) == 2L && inherits(other$args[[1L]], tf_class) &&
      inherits(other$args[[2L]], "Series")) {
    if (!.lti_eq(self$var, other$var))
      .physics_error("ValueError", .py_filldedent("\n                    All the transfer functions should use the same complex\n                    variable of the Laplace transform or z-transform."))
    a <- unique(self$args); b <- unique(other$args[[2L]]$args)
    ka <- vapply(a, .pyobj_k, ""); kb <- vapply(b, .pyobj_k, "")
    res <- c(a[!ka %in% kb], b[!kb %in% ka])
    if (length(res) == 0L) return(Feedback(self, other$args[[1L]]))
    if (length(res) == 1L) return(Feedback(self, res[[1L]]))
    return(Feedback(self, do.call(Series, res)))
  }
  .physics_error("ValueError", "This transfer function expression is invalid.")
}

#' @rdname control-interconnections
#' @export
Feedback <- function(sys1, sys2 = NULL, sign = -1L) {
  if (is.null(sys2)) sys2 <- create_transfer_function(1L, 1L, sys1$var, sys1$sampling_time)
  ok <- c("TransferFunctionBase", "Series", "StateSpaceBase", "Feedback")
  if (!inherits(sys1, ok)) .physics_error("TypeError", "Unsupported type for `sys1` in Feedback.")
  if (!inherits(sys2, ok)) .physics_error("TypeError", "Unsupported type for `sys2` in Feedback.")
  if (!(identical(sys1$num_inputs, 1L) && identical(sys1$num_outputs, 1L) && identical(sys2$num_inputs, 1L) && identical(sys2$num_outputs, 1L)))
    .physics_error("ValueError", .py_filldedent("To use Feedback connection for MIMO systems\n                            use MIMOFeedback instead."))
  sgn <- if (is_sym(sign)) as.numeric(N(sign)) else sign
  if (!(is.numeric(sgn) && length(sgn) == 1L && sgn %in% c(-1, 1)))
    .physics_error("ValueError", .py_filldedent("\n                Unsupported type for feedback. `sign` arg should\n                either be 1 (positive feedback loop) or -1\n                (negative feedback loop)."))
  ss <- isTRUE(sys1$is_StateSpace_object) || isTRUE(sys2$is_StateSpace_object)
  if (!ss) {
    prod <- simplify(.mul(list(sys1$to_expr(), sys2$to_expr())))
    if (.lti_eq(prod, as_sym(sign))) .physics_error("ValueError", "\nThe equivalent system will have zero\ndenominator.")
    if (!.lti_eq(sys1$var, sys2$var))
      .physics_error("ValueError", "\nBoth `sys1` and `sys2` should be                     using the same\ncomplex variable.")
  }
  .lti_check_time(list(sys1, sys2))
  .pyobj_new("Feedback", list(sys1, sys2, as_sym(sign)),
             extra = list(is_StateSpace_object = ss, is_continuous = sys1$is_continuous))
}

.feedback_ss_doit <- function(self, create) {
  ss_class <- .lti_ss_class(self$is_continuous)
  s1 <- self$sys1$doit()$rewrite(ss_class)
  s2 <- self$sys2$doit()$rewrite(ss_class)
  A1 <- s1$A; B1 <- s1$B; C1 <- s1$C; D1 <- s1$D
  A2 <- s2$A; B2 <- s2$B; C2 <- s2$C; D2 <- s2$D
  sg <- self$sign
  I_in <- eye(self$num_inputs); I_out <- eye(self$num_outputs)
  F <- I_in - (D2 %*% D1) * sg
  E <- inverse(F)
  E_D2 <- E %*% D2; E_C2 <- E %*% C2
  T1 <- I_out + (D1 %*% E_D2) * sg
  T2 <- I_in + (E_D2 %*% D1) * sg
  A <- rbind(cbind(A1 + (B1 %*% E_D2 %*% C1) * sg, (B1 %*% E_C2) * sg),
             cbind(B2 %*% T1 %*% C1, A2 + (B2 %*% D1 %*% E_C2) * sg))
  B <- rbind(B1 %*% T2, B2 %*% D1 %*% T2)
  C <- cbind(T1 %*% C1, (D1 %*% E_C2) * sg)
  D <- D1 %*% T2
  create_state_space(A, B, C, D, self$sampling_time)
}

local({
  .pyclass_def("Feedback", parent = "SISOLinearTimeInvariant",
    ctor = function(sys1, sys2, sign) Feedback(sys1, sys2, sign),
    props = list(
      sys1 = function(self) self$args[[1L]], sys2 = function(self) self$args[[2L]],
      var = function(self) self$args[[1L]]$var, sign = function(self) self$args[[3L]],
      num = function(self) self$args[[1L]],
      den = function(self) {
        unit <- create_transfer_function(1L, 1L, self$var, self$args[[1L]]$sampling_time)
        arg_list <- if (inherits(self$sys1, "Series")) self$sys1$args else list(self$sys1)
        if (identical(self$sign$key, S_One$key)) return(Parallel(unit, -do.call(Series, c(list(self$sys2), arg_list))))
        Parallel(unit, do.call(Series, c(list(self$sys2), arg_list)))
      },
      sensitivity = function(self) 1 / (1 - self$sign * .lti_eval(self$sys1$to_expr()) * .lti_eval(self$sys2$to_expr())),
      sampling_time = function(self) self$sys1$sampling_time),
    methods = list(
      doit = function(self, cancel = FALSE, expand = FALSE, ...) {
        if (isTRUE(self$is_StateSpace_object)) return(.feedback_ss_doit(self))
        arg_list <- if (inherits(self$sys1, "Series")) self$sys1$args else list(self$sys1)
        F_n <- self$sys1$doit()
        unit <- create_transfer_function(1L, 1L, self$sys1$var, self$sys1$sampling_time)
        F_d <- if (identical(self$sign$key, S_NegativeOne$key)) Parallel(unit, do.call(Series, c(list(self$sys2), arg_list)))$doit()
               else Parallel(unit, -do.call(Series, c(list(self$sys2), arg_list)))$doit()
        res <- create_transfer_function(F_n$num * F_d$den, F_n$den * F_d$num, F_n$var, self$sys1$sampling_time)
        if (isTRUE(cancel)) res <- res$simplify()
        if (isTRUE(expand)) res <- res$expand()
        res
      },
      `_eval_rewrite_as_TransferFunction` = function(self) .lti_rewrite_tf(self, TRUE),
      `_eval_rewrite_as_DiscreteTransferFunction` = function(self) .lti_rewrite_tf(self, FALSE),
      to_expr = function(self) self$doit()$to_expr()),
    ops = list(`__neg__` = function(self) Feedback(-self$sys1, -self$sys2, self$sign)))
})

# ---------------------------------------------------------------------------
# discretization and margins

#' Discretization and stability margins
#'
#' `gbt(tf, sample_per, alpha)` gives the numerator and denominator
#' coefficients of the discretization of a continuous transfer function by
#' the generalized bilinear transformation; `bilinear()`, `forward_diff()`
#' and `backward_diff()` are the cases `alpha = 1/2`, 0 and 1.
#' `phase_margin()` and `gain_margin()` are the phase margin (in degrees)
#' and the gain margin (in decibels) of a SISO system. After
#' `sympy.physics.control.lti`.
#'
#' @param tf A continuous transfer function.
#' @param sample_per The sampling period.
#' @param alpha The parameter of the transformation.
#' @param system A SISO system.
#' @return A list of two coefficient lists, or a `sym` expression.
#' @examples
#' s <- Symbol("s"); T <- Symbol("T")
#' bilinear(TransferFunction(1, s + 1, s), T)
#' gain_margin(TransferFunction(1, s^3 + 2*s^2 + 2*s + 1, s))
#' @export
gbt <- function(tf, sample_per, alpha) {
  if (!isTRUE(tf$is_SISO)) .physics_error("NotImplementedError", "Not implemented for MIMO systems.")
  T <- .lti_sympify(sample_per)
  s <- tf$var
  z <- s
  np <- all_coeffs(Poly(tf$num, s))
  dp <- all_coeffs(Poly(tf$den, s))
  alpha <- .lti_limit_denominator(alpha, 1000L)
  N <- max(length(np), length(dp)) - 1L
  build <- function(cs) {
    cs <- rev(cs)
    .add(lapply(seq_along(cs), function(k) {
      i <- k - 1L
      T^(N - i) * cs[[k]] * (z - 1L)^i * (alpha * z + 1L - alpha)^(N - i)
    }))
  }
  num <- build(np); den <- build(dp)
  num_coefs <- all_coeffs(Poly(num, z))
  den_coefs <- all_coeffs(Poly(den, z))
  para <- den_coefs[[1L]]
  .pytuple(list(lapply(num_coefs, function(cf) cf / para), lapply(den_coefs, function(cf) cf / para)))
}

# Rational(alpha).limit_denominator(n)
.lti_limit_denominator <- function(alpha, n) {
  a <- .lti_sympify(alpha)
  if (a$type %in% c("Integer", "Rational")) {
    q <- a$value
  } else {
    q <- gmp::as.bigq(.py_float_to_ratio(as.numeric(N(a))))
  }
  .limit_den_bigq(gmp::as.bigq(q), n)
}

.py_float_to_ratio <- function(x) {
  # the exact binary value of a double
  if (x == 0) return(gmp::as.bigq(0))
  e <- 0
  while (x != floor(x)) { x <- x * 2; e <- e + 1 }
  gmp::as.bigq(gmp::as.bigz(sprintf("%.0f", x)), gmp::pow.bigz(2, e))
}

# fractions.Fraction.limit_denominator
.limit_den_bigq <- function(q, max_den) {
  if (gmp::denominator(q) <= max_den) return(as_sym(q))
  p0 <- gmp::as.bigz(0); q0 <- gmp::as.bigz(1); p1 <- gmp::as.bigz(1); q1 <- gmp::as.bigz(0)
  n <- gmp::numerator(q); d <- gmp::denominator(q)
  repeat {
    a <- n %/% d
    q2 <- q0 + a * q1
    if (q2 > max_den) break
    t <- p0 + a * p1; p0 <- p1; q0 <- q1; p1 <- t; q1 <- q2
    t <- n - a * d; n <- d; d <- t
  }
  k <- (max_den - q0) %/% q1
  b1 <- gmp::as.bigq(p0 + k * p1, q0 + k * q1)
  b2 <- gmp::as.bigq(p1, q1)
  if (abs(b2 - q) <= abs(b1 - q)) as_sym(b2) else as_sym(b1)
}

#' @rdname gbt
#' @export
bilinear <- function(tf, sample_per) gbt(tf, sample_per, S_Half)

#' @rdname gbt
#' @export
forward_diff <- function(tf, sample_per) gbt(tf, sample_per, S_Zero)

#' @rdname gbt
#' @export
backward_diff <- function(tf, sample_per) gbt(tf, sample_per, S_One)

.lti_margin_setup <- function(system) {
  if (!inherits(system, "SISOLinearTimeInvariant")) .physics_error("ValueError", "Margins are only applicable for SISO LTI systems.")
  w <- Dummy("w", real = TRUE)
  expr <- system$to_expr()
  if (.lti_has_fn(expr, "exp")) .physics_error("NotImplementedError", "Margins for systems with Time delay terms are not supported.")
  if (length(free_symbols(expr)) > 1L)
    .physics_error("ValueError", "Extra degree of freedom found. Make sure that there are no free symbols in the dynamical system other than the variable of Laplace transform.")
  w_expr <- subs(expr, system$var, S_I * w)
  list(w = w, w_expr = w_expr, mag = 20L * fcall("log", list(Abs(w_expr), .int(10))))
}

.lti_has_fn <- function(e, name) {
  if (!is_sym(e)) return(FALSE)
  if (e$type == "Function" && identical(e$name, name)) return(TRUE)
  any(vapply(e$args, .lti_has_fn, TRUE, name = name))
}

.lti_solveset_pos <- function(e, w) {
  r <- solveset(e, w, Interval(S_Zero, oo, left_open = TRUE))
  if (is_sym(r) && r$type == "FiniteSet") return(r$args)
  if (is_sym(r) && r$type == "EmptySet") return(list())
  as.list(r$args)
}

#' @rdname gbt
#' @export
phase_margin <- function(system) {
  m <- .lti_margin_setup(system)
  mag_sol <- .lti_solveset_pos(m$mag, m$w)
  if (length(mag_sol) == 0L) {
    pm <- as_sym(-180L)
  } else {
    wcp <- mag_sol[[1L]]
    pm <- Mod(subs(arg(m$w_expr) * 180L / S_Pi, m$w, wcp) + 180L, .int(360))
  }
  if (isTRUE(pm >= 180L)) pm <- pm - 360L
  pm
}

#' @rdname gbt
#' @export
gain_margin <- function(system) {
  m <- .lti_margin_setup(system)
  ri <- .as_real_imag(m$w_expr)
  phase_sol <- .lti_solveset_pos(numer(cancel(ri[[2L]])), m$w)
  if (length(phase_sol) == 0L) return(oo)
  wcg <- phase_sol[[1L]]
  -subs(m$mag, m$w, wcg)
}

# the Python class names, for isinstance() style checks
attr(Series, "pyname") <- "Series"
attr(Parallel, "pyname") <- "Parallel"
attr(Feedback, "pyname") <- "Feedback"
