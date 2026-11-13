# Python code printers, after sympy/printing/pycode.py, numpy.py and
# lambdarepr.py.

.cg_py_kw <- c("and", "as", "assert", "break", "class", "continue", "def", "del", "elif", "else", "except",
               "finally", "for", "from", "global", "if", "import", "in", "is", "lambda", "not", "or", "pass",
               "raise", "return", "try", "while", "with", "yield", "None", "False", "nonlocal", "True")
.cg_kf_py <- list(Abs = "abs", Min = "min", Max = "max")
.cg_kf_math <- c(acos = "acos", acosh = "acosh", asin = "asin", asinh = "asinh", atan = "atan", atan2 = "atan2",
                 atanh = "atanh", ceiling = "ceil", cos = "cos", cosh = "cosh", erf = "erf", erfc = "erfc",
                 exp = "exp", expm1 = "expm1", factorial = "factorial", floor = "floor", gamma = "gamma",
                 hypot = "hypot", isinf = "isinf", isnan = "isnan", loggamma = "lgamma", log = "log", ln = "log",
                 log10 = "log10", log1p = "log1p", log2 = "log2", sin = "sin", sinh = "sinh", Sqrt = "sqrt",
                 tan = "tan", tanh = "tanh")
.cg_kc_math <- c(Exp1 = "e", Pi = "pi", E = "e", Infinity = "inf", `NaN` = "nan", ComplexInfinity = "nan")
.cg_prefixed <- function(x, pre) as.list(setNames(paste0(pre, x), names(x)))

.cg_kf_numpy <- local({
  m <- .cg_kf_math[!(names(.cg_kf_math) %in% c("erf", "erfc", "factorial", "gamma", "loggamma"))]
  m[c("acos", "acosh", "asin", "asinh", "atan", "atan2", "atanh")] <-
    c("arccos", "arccosh", "arcsin", "arcsinh", "arctan", "arctan2", "arctanh")
  c(m, exp2 = "exp2", sign = "sign", logaddexp = "logaddexp", logaddexp2 = "logaddexp2")
})
.cg_kc_numpy <- c(Exp1 = "e", Pi = "pi", EulerGamma = "euler_gamma", `NaN` = "nan", Infinity = "inf")
.cg_kf_scipy <- c(Ei = "expi", erf = "erf", erfc = "erfc", besselj = "jv", bessely = "yv", besseli = "iv",
                  besselk = "kv", cosm1 = "cosm1", powm1 = "powm1", factorial = "factorial", gamma = "gamma",
                  loggamma = "gammaln", digamma = "psi", polygamma = "polygamma", RisingFactorial = "poch",
                  jacobi = "eval_jacobi", gegenbauer = "eval_gegenbauer", chebyshevt = "eval_chebyt",
                  chebyshevu = "eval_chebyu", legendre = "eval_legendre", hermite = "eval_hermite",
                  laguerre = "eval_laguerre", assoc_laguerre = "eval_genlaguerre", beta = "beta", LambertW = "lambertw")
.cg_kc_scipy <- c(GoldenRatio = "golden_ratio", Pi = "pi")
.cg_kf_mpmath <- c(.cg_kf_math[!(names(.cg_kf_math) %in% c("log1p", "log2"))], beta = "beta", frac = "frac",
                   fresnelc = "fresnelc", fresnels = "fresnels", sign = "sign", loggamma = "loggamma",
                   hyper = "hyper", meijerg = "meijerg", besselj = "besselj", bessely = "bessely",
                   besseli = "besseli", besselk = "besselk")
.cg_kc_mpmath <- c(Exp1 = "e", Pi = "pi", GoldenRatio = "phi", EulerGamma = "euler", Catalan = "catalan",
                   `NaN` = "nan", Infinity = "inf", NegativeInfinity = "ninf")

# module_imports: a defaultdict(set) of module name to imported names
.cg_modimports <- function() structure(list(), class = "cg_modimports")
#' @export
format.cg_modimports <- function(x, ...) {
  items <- vapply(names(x), function(k) sprintf("'%s': %s", k, format(.cg_pyset(x[[k]]))), "")
  sprintf("defaultdict(<class 'set'>, {%s})", paste(items, collapse = ", "))
}
#' @export
print.cg_modimports <- function(x, ...) {
  cat(format(x), "\n")
  invisible(x)
}

.cgp_module_format <- function(p, fqn, register = TRUE) {
  parts <- strsplit(fqn, ".", fixed = TRUE)[[1L]]
  if (register && length(parts) > 1L) {
    mod <- paste(parts[-length(parts)], collapse = ".")
    mi <- unclass(p$module_imports)
    mi[[mod]] <- unique(c(mi[[mod]], parts[length(parts)]))
    p$module_imports <- structure(mi, class = "cg_modimports")
  }
  if (isTRUE(.cgp_s(p, "fully_qualified_modules"))) return(fqn)
  s <- strsplit(strsplit(fqn, "(", fixed = TRUE)[[1L]][1L], "[", fixed = TRUE)[[1L]][1L]
  tail(strsplit(s, ".", fixed = TRUE)[[1L]], 1L)
}
.cgp_print_known_func <- function(p, e) {
  known <- p$known_functions[[.cgp_clsname(e)]]
  sprintf("%s(%s)", .cgp_module_format(p, known), paste(vapply(e$args, function(a) .cgpp(p, a), ""), collapse = ", "))
}
.cgp_print_known_const <- function(p, e) .cgp_module_format(p, p$known_constants[[.cgp_clsname(e)]])

.cgp_hprint_pow <- function(p, e, rational = FALSE, sqrt = "math.sqrt") {
  PREC <- .cgp_prec(e)
  b <- e$args[[1L]]; x <- e$args[[2L]]
  if (identical(x$key, S_Half$key) && !rational)
    return(sprintf("%s(%s)", .cgp_module_format(p, sqrt), .cgpp(p, b)))
  if (!rational && !isFALSE(ask_fact(e, "commutative"))) {
    if (identical(x$key, Rational(-1, 2)$key))
      return(sprintf("%s/%s(%s)", .cgpp(p, S_One), .cgp_module_format(p, sqrt), .cgpp(p, b)))
    if (identical(x$key, S_NegativeOne$key)) return(sprintf("%s/%s", .cgpp(p, S_One), .cgp_paren(p, b, PREC)))
  }
  sprintf("%s**%s", .cgp_paren(p, b, PREC), .cgp_paren(p, x, PREC))
}
.cgp_indent_codestring <- function(p, s) paste(paste0("    ", .py_splitlines(s)), collapse = "\n")

.cgp_py_relational <- function(p, e) {
  op <- .REL_OPS[[e$type]]
  if (op %in% c("==", "!=", "<", "<=", ">", ">="))
    return(sprintf("(%s %s %s)", .cgpp(p, e$args[[1L]]), op, .cgpp(p, e$args[[2L]])))
  .cgp_super(p, "AbstractPythonCodePrinter", "_print_Relational", e)
}

.cgp_def("AbstractPythonCodePrinter", "CodePrinter",
         attrs = list(printmethod = "_pythoncode", language = "Python", reserved_words = .cg_py_kw,
                      module = "sympy.printing.pycode", tab = "    ",
                      `_kf` = c(.cg_kf_py, .cg_prefixed(.cg_kf_math, "math.")), `_kc` = .cg_prefixed(.cg_kc_math, "math."),
                      `_operators` = list(and = "and", or = "or", not = "not"),
                      `_default_settings` = list(order = NULL, full_prec = "auto", error_on_reserved = FALSE,
                                                 reserved_word_suffix = "_", human = TRUE, inline = TRUE,
                                                 allow_unknown_functions = FALSE, strict = NULL, user_functions = list(),
                                                 precision = 17, fully_qualified_modules = TRUE, contract = FALSE,
                                                 standard = "python3")),
         init = function(p, settings, nxt) {
           nxt(settings)
           std <- .cgp_s(p, "standard")
           if (is.null(std)) std <- "python3"
           if (std != "python3") .cg_err("ValueError", "Only Python 3 is supported.")
           p$standard <- std
           p$module_imports <- .cg_modimports()
           kf <- .cgp_a(p, "_kf")
           for (k in names(settings$user_functions)) kf[k] <- list(settings$user_functions[[k]])
           p$known_functions <- kf
           kc <- .cgp_a(p, "_kc")
           for (k in names(settings$user_constants)) kc[k] <- list(settings$user_constants[[k]])
           p$known_constants <- kc
           invisible(NULL)
         },
         methods = list(
           `_declare_number_const` = function(p, name, value) sprintf("%s = %s", .py_str_of(name), .py_str_of(value)),
           `_format_code` = function(p, lines) lines,
           `_get_statement` = function(p, s) s,
           `_get_comment` = function(p, text) sprintf("  # %s", text),
           `_print_NaN` = function(p, e) "float('nan')",
           `_print_Infinity` = function(p, e) "float('inf')",
           `_print_NegativeInfinity` = function(p, e) "float('-inf')",
           `_print_ComplexInfinity` = function(p, e) "float('nan')",
           `_print_Assignment` = function(p, e) {
             if (e$args[[2L]]$type == "Piecewise") return(sprintf("%s = %s", .cgpp(p, e$args[[1L]]), .cgpp(p, e$args[[2L]])))
             .cgp_super(p, "AbstractPythonCodePrinter", "_print_Assignment", e)
           },
           `_print_Mod` = function(p, e) {
             PREC <- .cgp_prec(e)
             num <- e$args[[1L]]; den <- e$args[[2L]]
             den_prec <- if (den$type == "Pow" && isTRUE(ask_fact(den$args[[2L]], "negative"))) .CG_PRECEDENCE[["Pow"]] else PREC
             sprintf("%s %% %s", .cgp_paren(p, num, PREC), .cgp_paren(p, den, den_prec))
           },
           `_print_Piecewise` = function(p, e) {
             res <- character()
             for (i in seq_along(e$args)) {
               if (i == 1L) res <- c(res, "(")
               res <- c(res, "(", .cgpp(p, e$args[[i]]$args[[1L]]), ")", " if ", .cgpp(p, e$args[[i]]$args[[2L]]), " else ")
             }
             res <- res[-length(res)]
             if (res[length(res)] == "True") res <- c(res[seq_len(length(res) - 2L)], ")")
             else res <- c(res, " else None)")
             paste(res, collapse = "")
           },
           `_print_Relational` = .cgp_py_relational,
           `_print_ImaginaryUnit` = function(p, e) "1j",
           `_print_conjugate` = function(p, e) sprintf("(%s).conjugate()", .cgpp(p, e$args[[1L]])),
           `_print_KroneckerDelta` = function(p, e) sprintf("(1 if %s == %s else 0)", .cgpp(p, e$args[[1L]]), .cgpp(p, e$args[[2L]])),
           `_print_MatrixBase` = function(p, e) {
             name <- "ImmutableDenseMatrix"
             func <- if (!is.null(p$known_functions[[name]])) p$known_functions[[name]] else name
             rows <- lapply(seq_len(e$nrow), function(i) lapply(seq_len(e$ncol), function(j) e[i, j]))
             sprintf("%s(%s)", func, .cgpp(p, rows))
           },
           `_print_FunctionDefinition` = function(p, fd) {
             v <- .cg_vals(fd)
             body <- paste(vapply(v$body$args, function(a) .cgpp(p, a), ""), collapse = "\n")
             sprintf("def %s(%s):\n%s", .cgpp(p, v$name),
                     paste(vapply(v$parameters$args, function(a) .cgpp(p, .cg_vals(a)$symbol), ""), collapse = ", "),
                     .cgp_indent_codestring(p, body))
           },
           `_print_While` = function(p, w) {
             v <- .cg_vals(w)
             body <- paste(vapply(v$body$args, function(a) .cgpp(p, a), ""), collapse = "\n")
             sprintf("while %s:\n%s", .cgpp(p, v$condition), .cgp_indent_codestring(p, body))
           },
           `_print_Declaration` = function(p, d) {
             v <- .cg_vals(d$args[[1L]])
             sprintf("%s = %s", .cgpp(p, v$symbol), .cgpp(p, v$value))
           },
           `_print_BreakToken` = function(p, e) "break",
           `_print_Return` = function(p, e) paste0("return ", .cgpp(p, e$args[[1L]])),
           `_print_Raise` = function(p, e) paste0("raise ", .cgpp(p, e$args[[1L]])),
           `_print_RuntimeError_` = function(p, e) sprintf("RuntimeError(%s)", .cgpp(p, e$args[[1L]])),
           `_print_Print` = function(p, e) {
             v <- .cg_vals(e)
             pa <- paste(vapply(v$print_args$args, function(a) .cgpp(p, a), ""), collapse = ", ")
             if (!.is_none(v$format_string)) pa <- sprintf("%s %% (%s), end=\"\"", .cgpp(p, v$format_string), pa)
             if (!is.null(v$file) && !.is_none(v$file)) pa <- paste0(pa, ", file=", .cgpp(p, v$file))
             sprintf("print(%s)", pa)
           },
           `_print_Stream` = function(p, e) {
             nm <- .cg_text(.cg_vals(e)$name)
             if (nm == "stdout") return(.cgp_module_format(p, "sys.stdout"))
             if (nm == "stderr") return(.cgp_module_format(p, "sys.stderr"))
             .cgpp(p, .cg_vals(e)$name)
           },
           `_print_NoneToken` = function(p, e) "None"
         ))

.cgp_def("PythonCodePrinter", "AbstractPythonCodePrinter",
         methods = list(
           `_print_sign` = function(p, e)
             sprintf("(0.0 if %s == 0 else %s(1, %s))", .cgpp(p, e$args[[1L]]), .cgp_module_format(p, "math.copysign"), .cgpp(p, e$args[[1L]])),
           `_print_Not` = function(p, e) paste0(.cgp_a(p, "_operators")$not, " ", .cgp_paren(p, e$args[[1L]], .cgp_prec(e))),
           `_print_IndexedBase` = function(p, e) .cgpp(p, e$args[[1L]]),
           `_print_Indexed` = function(p, e)
             sprintf("%s[%s]", .cgpp(p, e$args[[1L]]), paste(vapply(e$args[-1L], function(i) .cgpp(p, i), ""), collapse = ", ")),
           `_print_Pow` = function(p, e, rational = FALSE) .cgp_hprint_pow(p, e, rational),
           `_print_Rational` = function(p, e) sprintf("%s/%s", as.character(gmp::numerator(e$value)), as.character(gmp::denominator(e$value))),
           `_print_Half` = function(p, e) "1/2",
           `_print_frac` = function(p, e) sprintf("(%s)", .cgp_call(p, "_print_Mod", sym_Mod(e$args[[1L]], S_One))),
           `_print_Symbol` = function(p, e) {
             name <- .cgp_super(p, "PythonCodePrinter", "_print_Symbol", e)
             if (name %in% .cgp_a(p, "reserved_words")) {
               if (isTRUE(.cgp_s(p, "error_on_reserved")))
                 .cg_err("ValueError", sprintf("This expression includes the symbol \"%s\" which is a reserved keyword in this language.", name))
               return(paste0(name, .cgp_s(p, "reserved_word_suffix")))
             }
             gsub("[{}]", "", name)
           },
           `_print_lowergamma` = function(p, e) .cgp_not_supported(p, e),
           `_print_uppergamma` = function(p, e) .cgp_not_supported(p, e),
           `_print_fresnelc` = function(p, e) .cgp_not_supported(p, e),
           `_print_fresnels` = function(p, e) .cgp_not_supported(p, e)
         ))
.cgp_add_known <- function(cls, funs, consts) {
  env <- .cgp_classes[[cls]]$methods
  for (k in funs) assign(paste0("_print_", k), .cgp_print_known_func, envir = env)
  for (k in consts) assign(paste0("_print_", k), .cgp_print_known_const, envir = env)
}
.cgp_add_known("PythonCodePrinter", c(names(.cg_kf_py), names(.cg_kf_math)), names(.cg_kc_math))

.cgp_mpf_tuple <- function(v) {
  if (v == 0) return("(0, 0, 0, 0)")
  q <- gmp::as.bigq(abs(v)); num <- gmp::numerator(q); ex <- -as.integer(log2(as.numeric(gmp::denominator(q))))
  while (gmp::mod.bigz(num, 2) == 0) { num <- num / 2; ex <- ex + 1L }
  sprintf("(%d, %s, %d, %d)", as.integer(v < 0), as.character(num), ex, nchar(as.character(num, b = 2)))
}
.cgp_def("MpmathPrinter", "PythonCodePrinter",
         attrs = list(printmethod = "_mpmathcode", language = "Python with mpmath",
                      `_kf` = c(.cg_kf_py, .cg_prefixed(.cg_kf_mpmath, "mpmath.")), `_kc` = .cg_prefixed(.cg_kc_mpmath, "mpmath.")),
         methods = list(
           `_print_Float` = function(p, e) sprintf("%s(%s)", .cgp_module_format(p, "mpmath.mpf"), .cgp_mpf_tuple(as.numeric(e$value))),
           `_print_Rational` = function(p, e) {
             f <- .cgp_module_format(p, "mpmath.mpf")
             sprintf("%s(%s)/%s(%s)", f, as.character(gmp::numerator(e$value)), f, as.character(gmp::denominator(e$value)))
           },
           `_print_Half` = function(p, e) .cgp_call(p, "_print_Rational", e),
           `_print_log2` = function(p, e) sprintf("%s(%s)/%s(2)", .cgp_module_format(p, "mpmath.log"), .cgpp(p, e$args[[1L]]), .cgp_module_format(p, "mpmath.log")),
           `_print_log1p` = function(p, e) sprintf("%s(%s)", .cgp_module_format(p, "mpmath.log1p"), .cgpp(p, e$args[[1L]])),
           `_print_Pow` = function(p, e, rational = FALSE) .cgp_hprint_pow(p, e, rational, "mpmath.sqrt")
         ))
.cgp_add_known("MpmathPrinter", c(names(.cg_kf_py), names(.cg_kf_mpmath)), names(.cg_kc_mpmath))

.cgp_lambda_join <- function(p, e, op) {
  args <- .cgp_sorted_args(e$args)
  paste0("(", paste(vapply(args, function(a) paste0("(", .cgpp(p, a), ")"), ""), collapse = paste0(" ", op, " ")), ")")
}
.cgp_def("LambdaPrinter", "PythonCodePrinter",
         attrs = list(printmethod = "_lambdacode", module = "sympy.printing.lambdarepr"),
         methods = list(
           `_print_And` = function(p, e) .cgp_lambda_join(p, e, "and"),
           `_print_Or` = function(p, e) .cgp_lambda_join(p, e, "or"),
           `_print_Not` = function(p, e) paste0("(not (", .cgpp(p, e$args[[1L]]), "))"),
           `_print_BooleanTrue` = function(p, e) "True",
           `_print_BooleanFalse` = function(p, e) "False",
           `_print_ITE` = function(p, e) sprintf("((%s) if (%s) else (%s))", .cgpp(p, e$args[[2L]]), .cgpp(p, e$args[[1L]]), .cgpp(p, e$args[[3L]])),
           `_print_NumberSymbol` = function(p, e) sstr(e),
           `_print_Pow` = function(p, e, ...) .cgp_super(p, "PythonCodePrinter", "_print_Pow", e, ...)
         ))

# NumPy and SciPy
.cgp_np <- function(p, name) .cgp_module_format(p, paste0(.cgp_a(p, "_module"), ".", name))
.cgp_np_minmax <- function(p, op, args) {
  if (!length(args)) .cg_err("NotImplementedError", sprintf("Need at least one argument for %s", op))
  if (length(args) == 1L) return(.cgpp(p, args[[1L]]))
  sprintf("%s(%s, [%s])", .cgp_module_format(p, "functools.reduce"), op, paste(vapply(args, function(a) .cgpp(p, a), ""), collapse = ", "))
}
.cgp_np_nest <- function(p, name, args) {
  f <- .cgp_np(p, name)
  n <- length(args)
  paste0(paste(vapply(args[-n], function(a) sprintf("%s(%s, ", f, .cgpp(p, a)), ""), collapse = ""), .cgpp(p, args[[n]]), strrep(")", n - 1L))
}
.cgp_def("NumPyPrinter", "PythonCodePrinter",
         attrs = list(`_module` = "numpy", `_kf` = .cg_prefixed(.cg_kf_numpy, "numpy."), `_kc` = .cg_prefixed(.cg_kc_numpy, "numpy."),
                      module = "sympy.printing.numpy"),
         init = function(p, settings, nxt) {
           p$language <- sprintf("Python with %s", .cgp_a(p, "_module"))
           p$printmethod <- sprintf("_%scode", .cgp_a(p, "_module"))
           p[["_kf"]] <- c(.cgp_classes$AbstractPythonCodePrinter$attrs[["_kf"]], .cgp_a(p, "_kf"))
           p[["_kf"]] <- p[["_kf"]][!duplicated(names(p[["_kf"]]), fromLast = TRUE)]
           nxt(settings)
         },
         methods = list(
           `_print_seq` = function(p, seq) sprintf("(%s,)", paste(vapply(seq, function(a) .cgpp(p, a), ""), collapse = ", ")),
           `_print_NegativeInfinity` = function(p, e) paste0("-", .cgpp(p, oo)),
           `_print_MatPow` = function(p, e) sprintf("%s(%s, %s)", .cgp_np(p, "linalg.matrix_power"), .cgpp(p, e$args[[1L]]), .cgpp(p, e$args[[2L]])),
           `_print_Inverse` = function(p, e) sprintf("%s(%s)", .cgp_np(p, "linalg.inv"), .cgpp(p, e$args[[1L]])),
           `_print_MatrixSolve` = function(p, e) sprintf("%s(%s, %s)", .cgp_np(p, "linalg.solve"), .cgpp(p, e$args[[1L]]), .cgpp(p, e$args[[2L]])),
           `_print_MatMul` = function(p, e) {
             c_m <- .cgp_coeff_matrices(e)
             l <- if (!identical(c_m$coeff$key, S_One$key)) c(c_m$matrices, list(c_m$coeff)) else e$args
             sprintf("(%s)", paste(vapply(l, function(i) .cgpp(p, i), ""), collapse = ").dot("))
           },
           `_print_Piecewise` = function(p, e) {
             exprs <- paste(vapply(e$args, function(a) .cgpp(p, a$args[[1L]]), ""), collapse = ",")
             conds <- paste(vapply(e$args, function(a) .cgpp(p, a$args[[2L]]), ""), collapse = ",")
             sprintf("%s([%s], [%s], default=%s)", .cgp_np(p, "select"), conds, exprs, .cgpp(p, nan))
           },
           `_print_Relational` = function(p, e) {
             ops <- c(`==` = "equal", `!=` = "not_equal", `<` = "less", `<=` = "less_equal", `>` = "greater", `>=` = "greater_equal")
             op <- .REL_OPS[[e$type]]
             if (op %in% names(ops))
               return(sprintf("%s(%s, %s)", .cgp_np(p, ops[[op]]), .cgpp(p, e$args[[1L]]), .cgpp(p, e$args[[2L]])))
             .cgp_super(p, "NumPyPrinter", "_print_Relational", e)
           },
           `_print_And` = function(p, e) sprintf("%s.reduce((%s))", .cgp_np(p, "logical_and"), paste(vapply(e$args, function(a) .cgpp(p, a), ""), collapse = ",")),
           `_print_Or` = function(p, e) sprintf("%s.reduce((%s))", .cgp_np(p, "logical_or"), paste(vapply(e$args, function(a) .cgpp(p, a), ""), collapse = ",")),
           `_print_Not` = function(p, e) sprintf("%s(%s)", .cgp_np(p, "logical_not"), paste(vapply(e$args, function(a) .cgpp(p, a), ""), collapse = ",")),
           `_print_Pow` = function(p, e, rational = FALSE) {
             x <- e$args[[2L]]
             if (isTRUE(ask_fact(x, "integer")) && isTRUE(ask_fact(x, "negative")))
               e <- new_sym("Pow", list(e$args[[1L]], N(x)))
             .cgp_hprint_pow(p, e, rational, paste0(.cgp_a(p, "_module"), ".sqrt"))
           },
           `_print_Min` = function(p, e) .cgp_call(p, "_print_minimum", e),
           `_print_amin` = function(p, e) { v <- .cg_vals(e); sprintf("%s(%s, axis=%s)", .cgp_np(p, "amin"), .cgpp(p, v$array), .cgpp(p, v$axis)) },
           `_print_minimum` = function(p, e) .cgp_np_minmax(p, .cgp_np(p, "minimum"), e$args),
           `_print_Max` = function(p, e) .cgp_call(p, "_print_maximum", e),
           `_print_amax` = function(p, e) { v <- .cg_vals(e); sprintf("%s(%s, axis=%s)", .cgp_np(p, "amax"), .cgpp(p, v$array), .cgpp(p, v$axis)) },
           `_print_maximum` = function(p, e) .cgp_np_minmax(p, .cgp_np(p, "maximum"), e$args),
           `_print_arg` = function(p, e) sprintf("%s(%s)", .cgp_np(p, "angle"), .cgpp(p, e$args[[1L]])),
           `_print_im` = function(p, e) sprintf("%s(%s)", .cgp_np(p, "imag"), .cgpp(p, e$args[[1L]])),
           `_print_Mod` = function(p, e) sprintf("%s(%s)", .cgp_np(p, "mod"), paste(vapply(e$args, function(a) .cgpp(p, a), ""), collapse = ", ")),
           `_print_frac` = function(p, e) sprintf("%s(%s, 1)", .cgp_np(p, "mod"), .cgpp(p, e$args[[1L]])),
           `_print_re` = function(p, e) sprintf("%s(%s)", .cgp_np(p, "real"), .cgpp(p, e$args[[1L]])),
           `_print_conjugate` = function(p, e) sprintf("%s(%s)", .cgp_np(p, "conjugate"), .cgpp(p, e$args[[1L]])),
           `_print_sinc` = function(p, e) sprintf("%s(%s)", .cgp_np(p, "sinc"), .cgpp(p, .div(e$args[[1L]], S_Pi))),
           `_print_MatrixBase` = function(p, e) {
             if (e$nrow == 0L || e$ncol == 0L)
               return(sprintf("%s(%s)", .cgp_np(p, .cgp_a(p, "_zeros")), .cgpp(p, structure(list(e$nrow, e$ncol), class = "pytuple"))))
             rows <- lapply(seq_len(e$nrow), function(i) lapply(seq_len(e$ncol), function(j) e[i, j]))
             sprintf("%s(%s)", .cgp_np(p, "array"), .cgpp(p, rows))
           },
           `_print_ZeroMatrix` = function(p, e) sprintf("%s(%s)", .cgp_np(p, "zeros"), .cgpp(p, Tuple(e$args[[1L]], e$args[[2L]]))),
           `_print_OneMatrix` = function(p, e) sprintf("%s(%s)", .cgp_np(p, "ones"), .cgpp(p, Tuple(e$args[[1L]], e$args[[2L]]))),
           `_print_FunctionMatrix` = function(p, e) {
             L <- e$args[[3L]]
             sprintf("%s(lambda %s: %s, %s)", .cgp_np(p, "fromfunction"),
                     paste(vapply(L$args[[1L]]$args, function(a) .cgpp(p, a), ""), collapse = ", "),
                     .cgpp(p, L$args[[2L]]), .cgpp(p, Tuple(e$args[[1L]], e$args[[2L]])))
           },
           `_print_HadamardProduct` = function(p, e) .cgp_np_nest(p, "multiply", e$args),
           `_print_KroneckerProduct` = function(p, e) .cgp_np_nest(p, "kron", e$args),
           `_print_Adjoint` = function(p, e) sprintf("%s(%s(%s))", .cgp_np(p, "conjugate"), .cgp_np(p, "transpose"), .cgpp(p, e$args[[1L]])),
           `_print_DiagonalOf` = function(p, e) sprintf("%s(%s(%s), (-1, 1))", .cgp_np(p, "reshape"), .cgp_np(p, "diag"), .cgpp(p, e$args[[1L]])),
           `_print_DiagMatrix` = function(p, e) sprintf("%s(%s)", .cgp_np(p, "diagflat"), .cgpp(p, e$args[[1L]])),
           `_print_DiagonalMatrix` = function(p, e) {
             sh <- .mshape(e)
             sprintf("%s(%s, %s(%s, %s))", .cgp_np(p, "multiply"), .cgpp(p, e$args[[1L]]), .cgp_np(p, "eye"), .cgpp(p, sh[[1L]]), .cgpp(p, sh[[2L]]))
           },
           `_print_Identity` = function(p, e) sprintf("%s(%s)", .cgp_np(p, "eye"), .cgpp(p, .mshape(e)[[1L]])),
           `_print_lowergamma` = function(p, e) .cgp_not_supported(p, e),
           `_print_uppergamma` = function(p, e) .cgp_not_supported(p, e)
         ))
.cgp_add_known("NumPyPrinter", names(.cg_kf_numpy), names(.cg_kc_numpy))
.cgp_def("SciPyPrinter", "NumPyPrinter",
         attrs = list(`_kf` = c(.cg_prefixed(.cg_kf_numpy, "numpy."), .cg_prefixed(.cg_kf_scipy, "scipy.special.")),
                      `_kc` = c(.cg_prefixed(.cg_kc_numpy, "numpy."), .cg_prefixed(.cg_kc_scipy, "scipy.constants."))),
         init = function(p, settings, nxt) {
           nxt(settings)
           p$language <- "Python with SciPy and NumPy"
         })
.cgp_add_known("SciPyPrinter", names(.cg_kf_scipy), names(.cg_kc_scipy))

#' @rdname code-printer-objects
#' @export
PythonCodePrinter <- function(settings = NULL) .cgp_new("PythonCodePrinter", settings)
#' @rdname code-printer-objects
#' @export
MpmathPrinter <- function(settings = NULL) .cgp_new("MpmathPrinter", settings)
#' @rdname code-printer-objects
#' @export
NumPyPrinter <- function(settings = NULL) .cgp_new("NumPyPrinter", settings)
#' @rdname code-printer-objects
#' @export
SciPyPrinter <- function(settings = NULL) .cgp_new("SciPyPrinter", settings)
#' @rdname code-printer-objects
#' @export
LambdaPrinter <- function(settings = NULL) .cgp_new("LambdaPrinter", settings)

.cg_pycode <- function(e, assign_to = NULL, ...) .cgp_doprint(.cgp_new("PythonCodePrinter", list(...)), e, assign_to)
