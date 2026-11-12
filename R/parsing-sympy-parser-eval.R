# Evaluation of the parsed Python expression for parse_expr (the eval_expr
# step of sympy/parsing/sympy_parser.py), the EvaluateFalseTransformer that
# parse_expr(evaluate = FALSE) applies to it, and parse_expr itself. Only
# the operators and the names of the namespaces are interpreted; the input
# is never evaluated as R code.

#' Parse a string into an expression with SymPy's parser
#'
#' `parse_expr()` is SymPy's `parse_expr`: the string is tokenized as
#' Python, rewritten by the token `transformations` (see
#' [standard_transformations]) and evaluated with SymPy's names. Unknown
#' names become symbols or undefined functions, numbers become exact
#' integers and floats, `x!` is a factorial and `0.[3]` a repeating
#' decimal. With `evaluate = FALSE` the operators build unevaluated
#' expressions. [sympify()] is unchanged and does not use this parser.
#'
#' @param s A character string.
#' @param local_dict A named list of names to use first.
#' @param transformations A list of token transformations, or `"all"` or
#'   `"implicit"` for SymPy's groups of them.
#' @param global_dict A named list to use instead of SymPy's names, or
#'   `NULL`.
#' @param evaluate If `FALSE`, arithmetic and the elementary functions are
#'   left unevaluated.
#' @return The value of the expression: usually a `sym`, but Python
#'   numbers, strings, tuples and lists are returned as `pe_int`,
#'   `pe_float`, `pe_str`, `pe_tuple` and `pe_list` objects.
#' @examples
#' parse_expr("x**2 + 3*x - 1/2")
#' parse_expr("2x sin x", transformations = "implicit")
#' parse_expr("0.[3] + x!")
#' parse_expr("1 + 2", evaluate = FALSE)
#' @export
parse_expr <- function(s, local_dict = NULL, transformations = standard_transformations,
                       global_dict = NULL, evaluate = TRUE) {
  if (is.null(local_dict)) local_dict <- list()
  else if (!(is.environment(local_dict) || (is.list(local_dict) && !is.object(local_dict) &&
                                             (!length(local_dict) || !is.null(names(local_dict))))))
    stop(.pe_type_error("expecting local_dict to be a dict"))
  if (!is.environment(local_dict) && "" %in% names(local_dict))
    stop(.pe_value_error("cannot use \"\" in local_dict"))
  if (!is.null(global_dict) && !(is.list(global_dict) && !is.object(global_dict) &&
                                 (!length(global_dict) || !is.null(names(global_dict)))))
    stop(.pe_type_error("expecting global_dict to be a dict"))
  if (is.null(transformations) || (is.logical(transformations) && !any(transformations))) transformations <- list()
  if (is.character(transformations)) {
    transformations <- switch(transformations, all = .pe_transformation_list(),
                              implicit = .pe_transformation_list()[1:6],
                              stop(.pe_value_error("unknown transformation group name")))
  }
  if (is.function(transformations)) stop(.pe_type_error("'function' object is not iterable"))
  loc <- if (is.environment(local_dict)) local_dict else list2env(local_dict, envir = new.env(parent = emptyenv()))
  code <- .pe_stringify(s, loc, global_dict, transformations)
  tree <- .pe_parse(code)
  if (!isTRUE(evaluate)) tree <- .pe_ef(tree)
  .pe_eval(tree, .pe_ns(loc, global_dict))
}

.pe_stringify <- function(s, local_dict, global_dict, transformations) {
  tokens <- .pe_tokenize(trimws(s))
  for (tr in transformations) {
    if (!is.function(tr))
      stop(.pe_type_error(sprintf("'%s' object is not callable", if (is.list(tr)) "tuple" else class(tr)[1L])))
    tokens <- tr(tokens, local_dict, global_dict)
  }
  .pe_untokenize(tokens)
}

# sympify() of a string as SymPy does it: parse_expr with the standard
# transformations and convert_xor, a parse failure becoming a SympifyError
.pe_sympify_string <- function(a) {
  tryCatch(parse_expr(gsub("\n", "", a, fixed = TRUE), transformations = c(standard_transformations, convert_xor)),
           TokenError = function(e) .pe_sympify_fail(a, "TokenError", e),
           SyntaxError = function(e) .pe_sympify_fail(a, "SyntaxError", e))
}
.pe_sympify_fail <- function(a, cls, e) {
  stop(.pe_err("SympifyError", sprintf("Sympify of expression 'could not parse %s' failed, because of exception being raised:\n%s: %s",
                                       .pe_repr_str(a), cls, conditionMessage(e)), "ValueError"))
}

# ---------------------------------------------------------------------------
# Values

.pe_class_name <- function(v) {
  if (is_sym(v)) {
    if (v$type == "Integer") {
      return(switch(as.character(v$value), "0" = "Zero", "1" = "One", "-1" = "NegativeOne", "Integer"))
    }
    if (v$type == "Rational" && identical(v$key, S_Half$key)) return("Half")
    if (v$type %in% c("Function", "AppliedUndef") && !is.null(v$name)) return(v$name)
    if (v$type == "Constant") return(switch(EXPR = v$name, pi = "Pi", E = "Exp1", I = "ImaginaryUnit", oo = "Infinity",
                                            "-oo" = "NegativeInfinity", zoo = "ComplexInfinity", nan = "NaN", v$name))
    return(v$type)
  }
  if (inherits(v, "symmatrix")) return("MutableDenseMatrix")
  if (inherits(v, "pe_int")) return("int")
  if (inherits(v, "pe_float")) return("float")
  if (inherits(v, "pe_str") || is.character(v)) return("str")
  if (is.logical(v)) return("bool")
  if (inherits(v, "pe_tuple")) return("tuple")
  if (inherits(v, "pe_list")) return("list")
  if (inherits(v, "pe_none") || is.null(v)) return("NoneType")
  if (inherits(v, "sym_function")) return("UndefinedFunction")
  if (inherits(v, "pe_callable")) return(if (attr(v, "pe_name") %in% c("abs", "pow", "int", "float", "str", "len", "list", "tuple"))
                                           "builtin_function_or_method" else "FunctionClass")
  if (is.function(v)) return("function")
  class(v)[1L]
}

.pe_to_sym <- function(v) {
  if (is_sym(v) || inherits(v, "symmatrix")) return(v)
  if (inherits(v, "pe_int")) return(make_rational(gmp::as.bigq(v$v)))
  if (inherits(v, "pe_float")) return(Float(v$v))
  if (inherits(v, "pe_complex")) return(.add(list(Float(v$re), .mul(list(Float(v$im), S_I)))))
  if (is.logical(v)) return(if (v) S_true else S_false)
  if (inherits(v, "pe_tuple") || inherits(v, "pe_list")) return(do.call(Tuple, lapply(unclass(v), .pe_to_sym)))
  if (inherits(v, "pe_str")) return(.pe_sympify_string(v$v))
  if (is.numeric(v) || is.character(v)) return(as_sym(v))
  stop(.pe_err("SympifyError", sprintf("SympifyError: %s", .pe_repr(v)), "ValueError"))
}

# an argument for a plain R function supplied in local_dict or global_dict
.pe_to_r <- function(v) {
  if (inherits(v, "pe_str")) return(v$v)
  if (inherits(v, c("pe_int", "pe_float"))) return(.pe_to_sym(v))
  if (inherits(v, c("pe_tuple", "pe_list"))) return(lapply(unclass(v), .pe_to_r))
  if (inherits(v, "pe_none")) return(NULL)
  v
}

.pe_truthy <- function(v) {
  if (is.logical(v)) return(isTRUE(v))
  if (inherits(v, "pe_int")) return(v$v != 0)
  if (inherits(v, "pe_float")) return(v$v != 0)
  if (inherits(v, "pe_str")) return(nchar(v$v) > 0L)
  if (inherits(v, c("pe_tuple", "pe_list", "pe_set"))) return(length(v) > 0L)
  if (inherits(v, "pe_none") || is.null(v)) return(FALSE)
  if (is_sym(v)) {
    if (identical(v$key, S_true$key)) return(TRUE)
    if (identical(v$key, S_false$key)) return(FALSE)
    if (v$type %in% .REL_TYPES) stop(.pe_type_error("cannot determine truth value of Relational"))
    if (v$type %in% c("Integer", "Rational", "Float")) return(!is_zero_num(v))
    return(TRUE)
  }
  TRUE
}

.pe_is_num <- function(v) inherits(v, c("pe_int", "pe_float")) || is.logical(v)
.pe_numval <- function(v) if (is.logical(v)) .pe_int(as.integer(v)) else v

.pe_op_error <- function(op, a, b) {
  stop(.pe_type_error(sprintf("unsupported operand type(s) for %s: '%s' and '%s'",
                              if (op == "**") "** or pow()" else op, .pe_class_name(a), .pe_class_name(b))))
}

.pe_binop <- function(op, a, b) {
  if (is.function(a) || is.function(b)) .pe_op_error(op, a, b)
  if (is_sym(a) || is_sym(b) || inherits(a, "symmatrix") || inherits(b, "symmatrix")) {
    if (inherits(a, c("pe_tuple", "pe_list", "pe_str", "pe_none")) || inherits(b, c("pe_tuple", "pe_list", "pe_str", "pe_none")))
      stop(.pe_type_error(sprintf("unsupported operand type(s) for %s: '%s' and '%s'", op, .pe_class_name(a), .pe_class_name(b))))
    a <- .pe_to_sym(a); b <- .pe_to_sym(b)
    if (inherits(a, "symmatrix") || inherits(b, "symmatrix")) {
      return(switch(op, "+" = a + b, "-" = a - b, "*" = , "@" = a %*% b, "/" = a / b, "**" = a ^ b,
                    stop(.pe_type_error(sprintf("unsupported operand type(s) for %s", op)))))
    }
    return(switch(op,
      "+" = .add(list(a, b)), "-" = .sub(a, b), "*" = .mul(list(a, b)), "/" = .div(a, b),
      "**" = .pow(a, b), "%" = fcall("Mod", list(a, b)), "//" = fcall("floor", list(.div(a, b))),
      "&" = And(a, b), "|" = Or(a, b), "^" = Xor(a, b),
      stop(.pe_type_error(sprintf("unsupported operand type(s) for %s: '%s' and '%s'", op, .pe_class_name(a), .pe_class_name(b))))))
  }
  if (.pe_is_num(a) && .pe_is_num(b)) {
    a <- .pe_numval(a); b <- .pe_numval(b)
    if (inherits(a, "pe_int") && inherits(b, "pe_int")) {
      x <- a$v; y <- b$v
      return(switch(op,
        "+" = .pe_int(x + y), "-" = .pe_int(x - y), "*" = .pe_int(x * y),
        "/" = { if (y == 0) stop(.pe_err("ZeroDivisionError", "division by zero")); .pe_float(as.double(x) / as.double(y)) },
        "//" = { if (y == 0) stop(.pe_err("ZeroDivisionError", "integer division or modulo by zero")); .pe_int(gmp::as.bigz(floor(gmp::as.bigq(x) / gmp::as.bigq(y)))) },
        "%" = { if (y == 0) stop(.pe_err("ZeroDivisionError", "integer modulo by zero")); .pe_int(x - y * gmp::as.bigz(floor(gmp::as.bigq(x) / gmp::as.bigq(y)))) },
        "**" = if (y >= 0) .pe_int(x^as.integer(y)) else .pe_float(as.double(x)^as.double(y)),
        "&" = .pe_int(bitwAnd(as.integer(x), as.integer(y))), "|" = .pe_int(bitwOr(as.integer(x), as.integer(y))),
        "^" = .pe_int(bitwXor(as.integer(x), as.integer(y))),
        stop(.pe_type_error(sprintf("unsupported operand type(s) for %s: 'int' and 'int'", op)))))
    }
    x <- if (inherits(a, "pe_int")) as.double(a$v) else a$v
    y <- if (inherits(b, "pe_int")) as.double(b$v) else b$v
    return(switch(op, "+" = .pe_float(x + y), "-" = .pe_float(x - y), "*" = .pe_float(x * y),
                  "/" = { if (y == 0) stop(.pe_err("ZeroDivisionError", "float division by zero")); .pe_float(x / y) },
                  "**" = .pe_float(x^y), "//" = .pe_float(floor(x / y)), "%" = .pe_float(x - y * floor(x / y)),
                  stop(.pe_type_error(sprintf("unsupported operand type(s) for %s: '%s' and '%s'", op, .pe_class_name(a), .pe_class_name(b))))))
  }
  if (inherits(a, "pe_str") && inherits(b, "pe_str") && op == "+") return(.pe_str_val(paste0(a$v, b$v)))
  if (inherits(a, c("pe_tuple", "pe_list")) && identical(class(a), class(b)) && op == "+")
    return(structure(c(unclass(a), unclass(b)), class = class(a)))
  stop(.pe_type_error(sprintf("unsupported operand type(s) for %s: '%s' and '%s'", op, .pe_class_name(a), .pe_class_name(b))))
}

.pe_unary <- function(op, x) {
  if (op == "Not") return(!.pe_truthy(x))
  if (is_sym(x) || inherits(x, "symmatrix")) {
    if (op == "USub") return(if (inherits(x, "symmatrix")) -x else .neg(x))
    if (op == "UAdd") return(x)
    if (x$type == "Integer") return(make_rational(gmp::as.bigq(-x$value - 1)))
    return(Not(x))
  }
  if (.pe_is_num(x)) {
    x <- .pe_numval(x)
    if (op == "UAdd") return(x)
    if (op == "USub") return(if (inherits(x, "pe_int")) .pe_int(-x$v) else .pe_float(-x$v))
    if (inherits(x, "pe_int")) return(.pe_int(-x$v - 1L))
  }
  stop(.pe_type_error(sprintf("bad operand type for unary %s: '%s'", switch(op, USub = "-", UAdd = "+", Invert = "~"), .pe_class_name(x))))
}

.pe_eq <- function(a, b) {
  if (is_sym(a) || is_sym(b)) {
    if (.pe_is_num(a)) a <- .pe_to_sym(a)
    if (.pe_is_num(b)) b <- .pe_to_sym(b)
    return(is_sym(a) && is_sym(b) && identical(a$key, b$key))
  }
  if (.pe_is_num(a) && .pe_is_num(b)) {
    a <- .pe_numval(a); b <- .pe_numval(b)
    return(as.double(a$v) == as.double(b$v))
  }
  identical(a, b)
}

.pe_compare <- function(op, a, b) {
  if (op == "==") return(.pe_eq(a, b))
  if (op == "!=") return(!.pe_eq(a, b))
  if (op == "is") return(identical(a, b) || (inherits(a, "pe_none") && inherits(b, "pe_none")))
  if (op == "is not") return(!(identical(a, b) || (inherits(a, "pe_none") && inherits(b, "pe_none"))))
  if (op %in% c("in", "not in")) {
    if (!inherits(b, c("pe_tuple", "pe_list", "pe_set")))
      stop(.pe_type_error(sprintf("argument of type '%s' is not iterable", .pe_class_name(b))))
    r <- any(vapply(unclass(b), .pe_eq, TRUE, b = a))
    return(if (op == "in") r else !r)
  }
  if (is_sym(a) || is_sym(b)) {
    a <- .pe_to_sym(a); b <- .pe_to_sym(b)
    return(switch(op, "<" = Lt(a, b), "<=" = Le(a, b), ">" = Gt(a, b), ">=" = Ge(a, b)))
  }
  if (.pe_is_num(a) && .pe_is_num(b)) {
    x <- as.double(.pe_numval(a)$v); y <- as.double(.pe_numval(b)$v)
    return(switch(op, "<" = x < y, "<=" = x <= y, ">" = x > y, ">=" = x >= y))
  }
  stop(.pe_type_error(sprintf("'%s' not supported between instances of '%s' and '%s'", op, .pe_class_name(a), .pe_class_name(b))))
}

# ---------------------------------------------------------------------------
# Names

.pe_lookup <- function(name, ns) {
  if (exists(name, envir = ns$local, inherits = FALSE)) return(get(name, envir = ns$local))
  if (is.null(ns$global)) {
    if (name %in% .PY_GLOBAL_OBJECTS || name %in% .PY_GLOBAL_OTHER) return(.pe_sympy_name(name))
  } else if (name %in% names(ns$global)) return(ns$global[[name]])
  b <- .pe_builtin(name)
  if (!is.null(b)) return(b)
  stop(.pe_err("NameError", sprintf("name '%s' is not defined", name)))
}

.pe_kw_evaluate <- function(kw) {
  v <- kw[["evaluate"]]
  is.null(v) || .pe_truthy(v)
}

.pe_builtin <- function(name) {
  switch(EXPR = name,
    abs = .pe_callable(function(args, kw) {
      a <- args[[1L]]
      if (.pe_is_num(a)) { a <- .pe_numval(a); return(if (inherits(a, "pe_int")) .pe_int(abs(a$v)) else .pe_float(abs(a$v))) }
      fcall("Abs", list(.pe_to_sym(a)))
    }, "abs"),
    pow = .pe_callable(function(args, kw) {
      if (length(args) == 3L) {
        m <- .pe_to_sym(args[[3L]])
        b <- .pe_to_sym(args[[1L]]); e <- .pe_to_sym(args[[2L]])
        if (all(vapply(list(b, e, m), function(v) v$type == "Integer", TRUE))) {
          r <- if (e$value < 0) gmp::powm(gmp::inv.bigz(gmp::as.bigz(b$value), gmp::as.bigz(m$value)), -gmp::as.bigz(e$value), gmp::as.bigz(m$value))
               else gmp::powm(gmp::as.bigz(b$value), gmp::as.bigz(e$value), gmp::as.bigz(m$value))
          return(make_rational(gmp::as.bigq(r)))
        }
        return(fcall("Mod", list(.pow(b, e), m)))
      }
      .pe_binop("**", args[[1L]], args[[2L]])
    }, "pow"),
    int = .pe_callable(function(args, kw) {
      a <- args[[1L]]
      if (inherits(a, "pe_str")) return(.pe_int(.str_bigz(a$v)))
      if (.pe_is_num(a)) { a <- .pe_numval(a); return(.pe_int(if (inherits(a, "pe_int")) a$v else trunc(a$v))) }
      .pe_int(gmp::as.bigz(trunc(num_value_dbl(.pe_to_sym(a)))))
    }, "int"),
    float = .pe_callable(function(args, kw) {
      a <- args[[1L]]
      if (inherits(a, "pe_str")) return(.pe_float(as.double(a$v)))
      if (.pe_is_num(a)) return(.pe_float(as.double(.pe_numval(a)$v)))
      .pe_float(num_value_dbl(.pe_to_sym(a)))
    }, "float"),
    str = .pe_callable(function(args, kw) .pe_str_val(.pe_str(args[[1L]])), "str"),
    len = .pe_callable(function(args, kw) .pe_int(length(args[[1L]])), "len"),
    list = .pe_callable(function(args, kw) .pe_list(if (length(args)) unclass(args[[1L]]) else list()), "list"),
    tuple = .pe_callable(function(args, kw) .pe_tuple(if (length(args)) unclass(args[[1L]]) else list()), "tuple"),
    max = .pe_callable(function(args, kw) do.call(Max, lapply(args, .pe_to_sym)), "max"),
    min = .pe_callable(function(args, kw) do.call(Min, lapply(args, .pe_to_sym)), "min"),
    NULL)
}

.PY_CONSTANTS <- c("pi", "E", "I", "oo", "zoo", "nan", "EulerGamma", "Catalan", "GoldenRatio",
                   "TribonacciConstant", "true", "false")

.pe_sympy_name <- function(name) {
  if (name %in% .PY_CONSTANTS)
    return(switch(EXPR = name, pi = S_Pi, E = S_E, I = S_I, oo = S_Infinity, zoo = S_ComplexInfinity, nan = S_NaN,
                  EulerGamma = S_EulerGamma, Catalan = S_Catalan, GoldenRatio = S_GoldenRatio,
                  TribonacciConstant = S_TribonacciConstant, true = S_true, false = S_false))
  sp <- .pe_special(name)
  if (!is.null(sp)) return(sp)
  b <- .pe_builtin(name)
  if (!is.null(b)) return(b)
  .pe_generic(name)
}

.pe_generic <- function(name) {
  force(name)
  .pe_callable(function(args, kw) {
    args <- lapply(args, .pe_to_sym)
    if (!.pe_kw_evaluate(kw)) return(fnode(name, args))
    pf <- get0(name, envir = .parse_functions, inherits = FALSE)
    if (!is.null(pf)) return(do.call(pf, args))
    if (!is.null(get0(name, envir = .FUNS, inherits = FALSE))) return(fcall(name, args))
    rf <- get0(name, envir = asNamespace("symcompR"), mode = "function", inherits = FALSE)
    if (!is.null(rf)) {
      kw[["evaluate"]] <- NULL
      return(do.call(rf, c(args, lapply(kw, .pe_to_r))))
    }
    stop(.pe_err("NotImplementedError", sprintf("%s is not available in symcompR", name)))
  }, name)
}

.pe_rational_arg <- function(v) {
  if (inherits(v, "pe_str")) {
    s <- trimws(v$v)
    if (grepl("/", s, fixed = TRUE)) {
      p <- strsplit(s, "/", fixed = TRUE)[[1L]]
      return(.div(.pe_rational_arg(.pe_str_val(p[1L])), .pe_rational_arg(.pe_str_val(p[2L]))))
    }
    return(make_rational(.decimal_to_bigq(s)))
  }
  if (inherits(v, "pe_float")) return(make_rational(.decimal_to_bigq(.pe_float_repr(v$v))))
  .pe_to_sym(v)
}

.pe_unevaluated <- function(type, args) {
  if (type %in% c("Add", "Mul") && !length(args)) return(if (type == "Add") S_Zero else S_One)
  # And, Or and Xor with evaluate=False print their arguments in this order
  if (type %in% c("And", "Or", "Xor")) return(new_sym(type, args, extra = list(uneval = TRUE)))
  new_sym(type, args)
}

.pe_special <- function(name) {
  f <- switch(EXPR = name,
    Symbol = function(args, kw) {
      facts <- lapply(kw, function(v) if (is.logical(v)) v else .pe_truthy(v))
      do.call(Symbol, c(list(.pe_to_r(args[[1L]])), facts))
    },
    Function = function(args, kw) Function(.pe_to_r(args[[1L]])),
    Integer = function(args, kw) {
      a <- args[[1L]]
      if (inherits(a, "pe_str")) return(make_rational(gmp::as.bigq(.str_bigz(a$v))))
      if (inherits(a, "pe_float")) return(make_rational(gmp::as.bigq(gmp::as.bigz(trunc(a$v)))))
      a <- .pe_to_sym(a)
      if (a$type == "Integer") return(a)
      make_rational(gmp::as.bigq(gmp::as.bigz(trunc(num_value_dbl(a)))))
    },
    Float = function(args, kw) {
      a <- args[[1L]]
      dps <- if (length(args) >= 2L) .pe_to_r(args[[2L]]) else kw[["dps"]]
      if (!is.null(dps) && is_sym(dps)) dps <- as.integer(dps$value)
      x <- if (inherits(a, "pe_str")) a$v else if (inherits(a, "pe_int")) as.character(a$v)
           else if (inherits(a, "pe_float")) a$v else .pe_to_sym(a)
      if (is.null(dps)) Float(x) else Float(x, dps = dps)
    },
    Rational = function(args, kw) {
      p <- .pe_rational_arg(args[[1L]])
      if (length(args) < 2L) return(p)
      .div(p, .pe_rational_arg(args[[2L]]))
    },
    Number = function(args, kw) {
      a <- args[[1L]]
      if (inherits(a, "pe_str")) return(if (grepl("[.eE]", a$v)) Float(a$v) else make_rational(gmp::as.bigq(.str_bigz(a$v))))
      .pe_to_sym(a)
    },
    Add = , Mul = function(args, kw) {
      args <- lapply(args, .pe_to_sym)
      if (!.pe_kw_evaluate(kw)) return(.pe_unevaluated(name, args))
      if (name == "Add") .add(args) else .mul(args)
    },
    Pow = function(args, kw) {
      args <- lapply(args, .pe_to_sym)
      if (!.pe_kw_evaluate(kw)) return(new_sym("Pow", args))
      .pow(args[[1L]], args[[2L]])
    },
    Eq = , Ne = , Lt = , Le = , Gt = , Ge = function(args, kw) {
      args <- lapply(args, .pe_to_sym)
      if (length(args) < 2L) args[[2L]] <- S_Zero
      do.call(name, list(args[[1L]], args[[2L]], evaluate = .pe_kw_evaluate(kw)))
    },
    And = , Or = , Xor = function(args, kw) {
      args <- lapply(args, .pe_to_sym)
      if (!.pe_kw_evaluate(kw)) return(.pe_unevaluated(name, args))
      do.call(name, args)
    },
    Not = function(args, kw) Not(.pe_to_sym(args[[1L]])),
    sqrt = function(args, kw) {
      a <- .pe_to_sym(args[[1L]])
      if (!.pe_kw_evaluate(kw)) return(new_sym("Pow", list(a, S_Half)))
      .pow(a, S_Half)
    },
    cbrt = function(args, kw) {
      a <- .pe_to_sym(args[[1L]])
      if (!.pe_kw_evaluate(kw)) return(new_sym("Pow", list(a, .rat(1, 3))))
      .pow(a, .rat(1, 3))
    },
    root = function(args, kw) {
      args <- lapply(args, .pe_to_sym)
      do.call(root, args)
    },
    ln = , log = function(args, kw) {
      args <- lapply(args, .pe_to_sym)
      if (!.pe_kw_evaluate(kw)) return(fnode("log", args))
      fcall("log", args)
    },
    Lambda = function(args, kw) {
      v <- args[[1L]]
      vars <- if (inherits(v, c("pe_tuple", "pe_list"))) lapply(unclass(v), .pe_to_sym) else list(.pe_to_sym(v))
      Lambda(vars, .pe_to_sym(args[[2L]]))
    },
    Max = , Min = function(args, kw) do.call(name, lapply(args, .pe_to_sym)),
    S = , sympify = function(args, kw) {
      a <- args[[1L]]
      if (inherits(a, "pe_str")) return(.pe_sympify_string(a$v))
      .pe_to_sym(a)
    },
    Limit = function(args, kw) {
      a <- lapply(args, .pe_to_sym)
      dir <- if (length(args) >= 4L) .pe_to_r(args[[4L]]) else if (!is.null(kw[["dir"]])) .pe_to_r(kw[["dir"]]) else "+"
      Limit(a[[1L]], a[[2L]], a[[3L]], dir)
    },
    Integral = , Sum = , Product = , Derivative = function(args, kw) {
      lims <- lapply(args[-1L], function(l) if (inherits(l, c("pe_tuple", "pe_list"))) lapply(unclass(l), .pe_to_sym) else .pe_to_sym(l))
      do.call(name, c(list(.pe_to_sym(args[[1L]])), lims))
    },
    Matrix = function(args, kw) Matrix(lapply(unclass(args[[1L]]), function(r) if (inherits(r, c("pe_tuple", "pe_list"))) lapply(unclass(r), .pe_to_sym) else list(.pe_to_sym(r)))),
    Tuple = function(args, kw) do.call(Tuple, lapply(args, .pe_to_sym)),
    NULL)
  if (is.null(f)) return(NULL)
  .pe_callable(f, name)
}

# ---------------------------------------------------------------------------
# The evaluator

.pe_eval <- function(node, ns) {
  switch(node$k,
    Const = node$v,
    Name = .pe_lookup(node$id, ns),
    Tuple = .pe_tuple(lapply(node$elts, .pe_eval, ns = ns)),
    List = .pe_list(lapply(node$elts, .pe_eval, ns = ns)),
    Set = structure(lapply(node$elts, .pe_eval, ns = ns), class = "pe_set"),
    BinOp = {
      a <- .pe_eval(node$l, ns)
      b <- .pe_eval(node$r, ns)
      .pe_binop(node$op, a, b)
    },
    UnaryOp = .pe_unary(node$op, .pe_eval(node$x, ns)),
    BoolOp = {
      v <- NULL
      for (x in node$values) {
        v <- .pe_eval(x, ns)
        t <- .pe_truthy(v)
        if ((node$op == "And" && !t) || (node$op == "Or" && t)) break
      }
      v
    },
    Compare = {
      left <- .pe_eval(node$left, ns)
      r <- NULL
      for (j in seq_along(node$ops)) {
        right <- .pe_eval(node$comps[[j]], ns)
        r <- .pe_compare(node$ops[j], left, right)
        if (j < length(node$ops) && !.pe_truthy(r)) return(r)
        left <- right
      }
      r
    },
    IfExp = if (.pe_truthy(.pe_eval(node$test, ns))) .pe_eval(node$body, ns) else .pe_eval(node$orelse, ns),
    Lambda = {
      params <- node$params; body <- node$body
      .pe_callable(function(args, kw) {
        loc <- new.env(parent = emptyenv())
        for (nm in ls(ns$local, all.names = TRUE)) assign(nm, get(nm, envir = ns$local), envir = loc)
        for (j in seq_along(params)) assign(params[j], args[[j]], envir = loc)
        .pe_eval(body, .pe_ns(loc, ns$global))
      }, "<lambda>")
    },
    Call = {
      f <- if (node$f$k == "Builtin") node$f else .pe_eval(node$f, ns)
      args <- list()
      for (a in node$args) {
        if (a$k == "Starred") {
          v <- .pe_eval(a$v, ns)
          args <- c(args, if (is_sym(v)) v$args else unclass(v))
        } else args[[length(args) + 1L]] <- .pe_eval(a, ns)
      }
      kw <- list()
      for (k in node$kw) kw[[k$name]] <- .pe_eval(k$v, ns)
      .pe_apply(f, args, kw)
    },
    Attr = .pe_getattr(.pe_eval(node$v, ns), node$attr),
    Subscript = .pe_getitem(.pe_eval(node$v, ns), .pe_eval(node$idx, ns)),
    stop(.pe_syntax_error()))
}

.pe_apply <- function(f, args, kw) {
  if (is.list(f) && identical(f$k, "Builtin")) {
    # the two lambdas EvaluateFalseTransformer writes for a division
    if (f$name == "ratio")
      return(.pe_unevaluated("Mul", list(.pe_to_sym(args[[1L]]), new_sym("Pow", list(.pe_to_sym(args[[2L]]), S_NegativeOne)))))
    n <- .pe_to_sym(args[[1L]]); d <- .pe_to_sym(args[[2L]])
    if (identical(n$key, S_One$key)) return(d)
    if (n$type == "Mul") return(.pe_unevaluated("Mul", c(n$args, list(d))))
    return(.pe_unevaluated("Mul", list(n, d)))
  }
  if (inherits(f, "pe_callable")) return(f(args, kw))
  if (inherits(f, "sym_function")) return(do.call(f, lapply(args, .pe_to_sym)))
  if (is.function(f)) return(do.call(f, c(lapply(args, .pe_to_r), lapply(kw, .pe_to_r))))
  stop(.pe_type_error(sprintf("'%s' object is not callable", .pe_class_name(f))))
}

.pe_getattr <- function(obj, attr) {
  if (is_sym(obj)) {
    v <- switch(attr,
      args = .pe_tuple(obj$args),
      free_symbols = structure(free_symbols(obj), class = "pe_set"),
      is_Integer = obj$type == "Integer", is_Symbol = obj$type == "Symbol",
      doit = .pe_callable(function(args, kw) doit(obj), "doit"),
      subs = .pe_callable(function(args, kw) subs(obj, .pe_to_sym(args[[1L]]), .pe_to_sym(args[[2L]])), "subs"),
      expand = .pe_callable(function(args, kw) expand(obj), "expand"),
      n = , evalf = .pe_callable(function(args, kw) {
        d <- if (length(args)) args[[1L]] else kw[["n"]]
        N(obj, if (is.null(d)) 15 else as.integer(.pe_numval(d)$v))
      }, attr),
      is_polynomial = .pe_callable(function(args, kw) is_polynomial(obj, lapply(args, .pe_to_sym)), attr),
      NULL)
    if (!is.null(v)) return(v)
  }
  stop(.pe_err("AttributeError", sprintf("'%s' object has no attribute '%s'", .pe_class_name(obj), attr)))
}

.pe_getitem <- function(obj, idx) {
  if (inherits(obj, c("pe_tuple", "pe_list")) && inherits(idx, "pe_int")) {
    i <- as.integer(idx$v)
    n <- length(obj)
    if (i < -n || i >= n) stop(.pe_err("IndexError", sprintf("%s index out of range", .pe_class_name(obj))))
    return(unclass(obj)[[if (i < 0L) n + i + 1L else i + 1L]])
  }
  stop(.pe_type_error(sprintf("'%s' object is not subscriptable", .pe_class_name(obj))))
}

# ---------------------------------------------------------------------------
# EvaluateFalseTransformer

.pe_ef_call <- function(fname, args, evaluate_false = TRUE) {
  kw <- if (evaluate_false) list(list(name = "evaluate", v = list(k = "Const", v = FALSE))) else list()
  list(k = "Call", f = list(k = "Name", id = fname), args = args, kw = kw)
}
.pe_ef_minus_one <- function() list(k = "UnaryOp", op = "USub", x = list(k = "Const", v = .pe_int(1L)))
.pe_ef_is_call_to <- function(node, fname) identical(node$k, "Call") && identical(node$f$k, "Name") && identical(node$f$id, fname)

.pe_ef_integer_literal <- function(node) {
  sign <- 1L
  repeat {
    if (node$k == "UnaryOp" && node$op %in% c("UAdd", "USub")) {
      if (node$op == "USub") sign <- -sign
      node <- node$x
    } else if (.pe_ef_is_call_to(node, "Integer") && length(node$args) == 1L && !length(node$kw)) {
      node <- node$args[[1L]]
    } else break
  }
  if (node$k == "Const" && inherits(node$v, "pe_int")) return(sign * node$v$v)
  NULL
}

.pe_ef_literal_fraction <- function(node) {
  q <- .pe_ef_integer_literal(node$r)
  if (is.null(q) || q < 0) return(NULL)
  p <- .pe_ef_integer_literal(node$l)
  if (!is.null(p) && p == 1 && q <= 1) return(NULL)
  if (is.null(p)) {
    if (q > 1 && node$l$k == "UnaryOp" && node$l$op == "USub" && node$l$x$k != "UnaryOp") {
      coeff <- list(k = "Call", f = list(k = "Name", id = "Rational"),
                    args = list(list(k = "Const", v = .pe_int(-1L)), list(k = "Const", v = .pe_int(q))), kw = list())
      return(.pe_ef_call("Mul", list(coeff, .pe_ef(node$l$x))))
    }
    return(NULL)
  }
  if (p != 0 && q > 1 && gmp::gcd.bigz(p, q) == 1)
    return(list(k = "Call", f = list(k = "Name", id = "Rational"),
                args = list(list(k = "Const", v = .pe_int(p)), list(k = "Const", v = .pe_int(q))), kw = list()))
  list(k = "Call", f = list(k = "Builtin", name = "ratio"), args = list(.pe_ef(node$l), .pe_ef(node$r)), kw = list())
}

.PY_EF_OPERATORS <- c("+" = "Add", "*" = "Mul", "**" = "Pow", "-" = "Add", "/" = "Mul", "|" = "Or", "&" = "And", "^" = "Xor")
.PY_EF_FUNCTIONS <- c("Abs", "im", "re", "sign", "arg", "conjugate", "acos", "acot", "acsc", "asec", "asin", "atan",
                      "acosh", "acoth", "acsch", "asech", "asinh", "atanh", "cos", "cot", "csc", "sec", "sin", "tan",
                      "cosh", "coth", "csch", "sech", "sinh", "tanh", "exp", "ln", "log", "sqrt", "cbrt", "factorial")
.PY_EF_RELATIONS <- c("!=" = "Ne", "<" = "Lt", "<=" = "Le", ">" = "Gt", ">=" = "Ge", "==" = "Eq")

.pe_ef_flatten <- function(args, fname) {
  out <- list()
  for (a in args) {
    if (identical(a$k, "Call")) {
      af <- a$f
      if (identical(af$k, "Call")) af <- af$f
      if (identical(af$k, "Name") && identical(af$id, fname)) { out <- c(out, .pe_ef_flatten(a$args, fname)); next }
    }
    out[[length(out) + 1L]] <- a
  }
  out
}

.pe_ef <- function(node) {
  if (!is.list(node) || is.null(node$k)) return(node)
  switch(node$k,
    Compare = {
      args <- list()
      left <- node$left
      for (j in seq_along(node$ops)) {
        op <- node$ops[j]
        if (!(op %in% names(.PY_EF_RELATIONS))) stop(.pe_value_error("Only equation or inequality operators are supported"))
        right <- node$comps[[j]]
        args[[length(args) + 1L]] <- .pe_ef_call(.PY_EF_RELATIONS[[op]], list(.pe_ef(left), .pe_ef(right)))
        left <- right
      }
      if (length(args) == 1L) args[[1L]] else .pe_ef_call("And", args)
    },
    UnaryOp = {
      if (node$op == "USub" && node$x$k == "BinOp" && node$x$op %in% c("+", "*", "/", "-")) {
        operand <- .pe_ef(node$x)
        if (.pe_ef_is_call_to(operand, "Rational")) return(list(k = "UnaryOp", op = "USub", x = operand))
        return(.pe_ef_call("Mul", list(.pe_ef_minus_one(), operand)))
      }
      node$x <- .pe_ef(node$x)
      node
    },
    BinOp = {
      if (node$op == "/") {
        fr <- .pe_ef_literal_fraction(node)
        if (!is.null(fr)) return(fr)
      }
      if (!(node$op %in% names(.PY_EF_OPERATORS))) return(node)
      cls <- .PY_EF_OPERATORS[[node$op]]
      right <- .pe_ef(node$r)
      left <- .pe_ef(node$l)
      if (node$op == "-") {
        if (node$r$k == "BinOp" && node$r$op == "/" && .pe_ef_is_call_to(right, "Rational"))
          right <- list(k = "UnaryOp", op = "USub", x = right)
        else right <- .pe_ef_call("Mul", list(.pe_ef_minus_one(), right))
      } else if (node$op == "/") {
        rev <- FALSE
        if (node$l$k == "UnaryOp") {
          tmp <- left; left <- right; right <- tmp
          rev <- TRUE
          left <- .pe_ef_call("Pow", list(left, .pe_ef_minus_one()))
        } else right <- .pe_ef_call("Pow", list(right, .pe_ef_minus_one()))
        if (rev) { tmp <- left; left <- right; right <- tmp }
        return(list(k = "Call", f = list(k = "Builtin", name = "divlam"), args = list(left, right), kw = list()))
      }
      new <- .pe_ef_call(cls, list(left, right))
      if (cls %in% c("Add", "Mul")) new$args <- .pe_ef_flatten(new$args, cls)
      new
    },
    Call = {
      if (node$f$k == "Name" && node$f$id %in% .PY_EF_FUNCTIONS) {
        args <- lapply(node$args, .pe_ef)
        kw <- lapply(node$kw, function(k) list(name = k$name, v = .pe_ef(k$v)))
        kw[[length(kw) + 1L]] <- list(name = "evaluate", v = list(k = "Const", v = FALSE))
        return(list(k = "Call", f = node$f, args = args, kw = kw))
      }
      if (node$f$k != "Builtin") node$f <- .pe_ef(node$f)
      node$args <- lapply(node$args, .pe_ef)
      node$kw <- lapply(node$kw, function(k) list(name = k$name, v = .pe_ef(k$v)))
      node
    },
    Tuple = , List = , Set = { node$elts <- lapply(node$elts, .pe_ef); node },
    BoolOp = { node$values <- lapply(node$values, .pe_ef); node },
    Attr = , Starred = { node$v <- .pe_ef(node$v); node },
    Subscript = { node$v <- .pe_ef(node$v); node$idx <- .pe_ef(node$idx); node },
    IfExp = { node$test <- .pe_ef(node$test); node$body <- .pe_ef(node$body); node$orelse <- .pe_ef(node$orelse); node },
    Lambda = { node$body <- .pe_ef(node$body); node },
    node)
}
