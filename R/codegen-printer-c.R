# C and C++ code printers, after sympy/printing/c.py and cxx.py.

.cg_kf_C89 <- list(Abs = list(list(function(x) !isTRUE(ask_fact(x, "integer")), "fabs"),
                              list(function(x) isTRUE(ask_fact(x, "integer")), "abs")),
                   sin = "sin", cos = "cos", tan = "tan", asin = "asin", acos = "acos", atan = "atan",
                   atan2 = "atan2", exp = "exp", log = "log", log10 = "log10", sinh = "sinh", cosh = "cosh",
                   tanh = "tanh", floor = "floor", ceiling = "ceil", sqrt = "sqrt")
.cg_kf_C99 <- c(.cg_kf_C89, list(exp2 = "exp2", expm1 = "expm1", log2 = "log2", log1p = "log1p", Cbrt = "cbrt",
                                 hypot = "hypot", fma = "fma", loggamma = "lgamma", erfc = "erfc", Max = "fmax",
                                 Min = "fmin", asinh = "asinh", acosh = "acosh", atanh = "atanh", erf = "erf",
                                 gamma = "tgamma"))
.cg_c_reserved <- c("auto", "break", "case", "char", "const", "continue", "default", "do", "double", "else",
                    "enum", "extern", "float", "for", "goto", "if", "int", "long", "register", "return", "short",
                    "signed", "sizeof", "static", "struct", "entry", "switch", "typedef", "union", "unsigned",
                    "void", "volatile", "while")
.cg_c99_reserved <- c("inline", "restrict")

#' C math macros
#'
#' `get_math_macros()` returns the constants that the C printers print as
#' the macros of `math.h`, as SymPy's `sympy.printing.c.get_math_macros`:
#' a dictionary from expressions to macro names, which `py_getitem()`-style
#' lookups read with `cgdict_get()`.
#'
#' @param d A dictionary from `get_math_macros()`.
#' @param key An expression.
#' @return `get_math_macros()` returns a dictionary; `cgdict_get()` the
#'   value for `key`, or `NULL`.
#' @examples
#' m <- get_math_macros()
#' cgdict_get(m, exp(Integer(1)))
#' @export
get_math_macros <- function() {
  .cgd(list(list(S_E, "M_E"), list(sym_log2(S_E), "M_LOG2E"), list(.div(S_One, fcall("log", list(S_Two))), "M_LOG2E"),
            list(fcall("log", list(S_Two)), "M_LN2"), list(fcall("log", list(Integer(10))), "M_LN10"),
            list(S_Pi, "M_PI"), list(.div(S_Pi, S_Two), "M_PI_2"), list(.div(S_Pi, Integer(4)), "M_PI_4"),
            list(.div(S_One, S_Pi), "M_1_PI"), list(.div(S_Two, S_Pi), "M_2_PI"),
            list(.div(S_Two, sqrt(S_Pi)), "M_2_SQRTPI"), list(.div(S_Two, Sqrt(S_Pi)), "M_2_SQRTPI"),
            list(sqrt(S_Two), "M_SQRT2"), list(Sqrt(S_Two), "M_SQRT2"),
            list(.div(S_One, sqrt(S_Two)), "M_SQRT1_2"), list(.div(S_One, Sqrt(S_Two)), "M_SQRT1_2")))
}
#' @rdname get_math_macros
#' @export
cgdict_get <- function(d, key) .cgd_get(d, key)

# the _as_macro_if_defined decorator
.cgp_macro_or <- function(p, e, f) {
  mm <- p$math_macros
  if (is_sym(e) && .cgd_has(mm, e)) return(paste0(.cgd_get(mm, e), .cgp_c_macro_suffix(p, real)))
  f()
}
# the requires decorator
.cgp_requires <- function(p, headers = character(), libraries = character()) {
  if (length(headers)) .cgp_update(p, "headers", headers)
  if (length(libraries)) .cgp_update(p, "libraries", libraries)
}

.cgp_alias <- function(p, t) .cgd_get(p$type_aliases, t, t)
.cgp_c_func_suffix <- function(p, t) {
  r <- .cgd_get(p$type_func_suffixes, .cgp_alias(p, t))
  if (is.null(r)) .cg_err("KeyError", .py_repr_of(.cgp_alias(p, t)))
  r
}
.cgp_c_literal_suffix <- function(p, t) {
  r <- .cgd_get(p$type_literal_suffixes, .cgp_alias(p, t))
  if (is.null(r)) .cg_err("KeyError", .py_repr_of(.cgp_alias(p, t)))
  r
}
.cgp_c_macro_suffix <- function(p, t) {
  alias <- .cgp_alias(p, t)
  dflt <- .cgd_get(p$type_math_macro_suffixes, alias, "")
  .cgd_get(p$type_math_macro_suffixes, t, dflt)
}

.cgp_equal_valued <- function(x, v) {
  if (!is_sym(x) || !is_number_type(x)) return(FALSE)
  if (x$type == "Float") return(isTRUE(as.numeric(x$value) == v))
  isTRUE(x$value == gmp::as.bigq(v))
}

# indent_code for brace languages
.cgp_indent_braces <- function(code, tab = "   ", inc_token = c("{", "(", "{\n", "(\n"), dec_token = c("}", ")")) {
  if (length(code) == 1L && is.character(code) && isTRUE(attr(code, "str"))) {
    lines <- .cgp_indent_braces(strsplit(code, "(?<=\n)", perl = TRUE)[[1L]], tab, inc_token, dec_token)
    return(paste(lines, collapse = ""))
  }
  code <- sub("^[ \t]+", "", code)
  increase <- vapply(code, function(l) any(endsWith(l, inc_token)), TRUE)
  decrease <- vapply(code, function(l) any(startsWith(l, dec_token)), TRUE)
  pretty <- character()
  level <- 0L
  for (n in seq_along(code)) {
    line <- code[n]
    if (line %in% c("", "\n")) { pretty <- c(pretty, line); next }
    level <- level - decrease[n]
    pretty <- c(pretty, paste0(strrep(tab, level), line))
    level <- level + increase[n]
  }
  unname(pretty)
}

# Indexed: the flat index of a C array
.cgp_indexed_shape <- function(e) {
  base <- e$args[[1L]]
  if (!is.null(base$extra$shape)) return(base$extra$shape)
  lapply(e$args[-1L], function(i) {
    r <- if (is_sym(i) && .is_idx(i)) i$extra$range else NULL
    if (is.null(r)) .cg_err("IndexException", sprintf("Range is not defined for all indices in: %s", sstr(e)))
    .add(list(.sub(r[[2L]], r[[1L]]), S_One))
  })
}

.cgp_c_print_indexed <- function(p, e) {
  base <- e$args[[1L]]
  offset <- if (!is.null(base$extra$offset)) base$extra$offset else S_Zero
  strides <- base$extra$strides
  indices <- e$args[-1L]
  if (is.null(strides) || is.character(strides)) {
    dims <- .cgp_indexed_shape(e)
    rank <- length(indices)
    shift <- S_One
    temp <- list()
    if (is.null(strides) || identical(strides, "C")) {
      traversal <- rev(seq_len(rank))
      indices <- rev(indices)
    } else traversal <- seq_len(rank)
    for (i in traversal) {
      temp[[length(temp) + 1L]] <- shift
      shift <- .mul(list(shift, dims[[i]]))
    }
    strides <- temp
  }
  flat <- .add(c(Map(function(a, b) .mul(list(a, b)), indices, strides), list(offset)))
  sprintf("%s[%s]", .cgpp(p, base$args[[1L]]), .cgpp(p, flat))
}

.cgp_c_print_piecewise <- function(p, e) {
  last <- e$args[[length(e$args)]]
  if (!.is_true(last$args[[2L]]))
    .cg_err("ValueError", paste0("All Piecewise expressions must contain an (expr, True) statement to be used as a default condition. ",
                                 "Without one, the generated expression may not evaluate to anything under some condition."))
  lines <- character()
  n <- length(e$args)
  if (.cg_has_type(e, "Assignment")) {
    for (i in seq_len(n)) {
      ec <- e$args[[i]]
      c1 <- ec$args[[2L]]
      if (i == 1L) lines <- c(lines, sprintf("if (%s) {", .cgpp(p, c1)))
      else if (i == n && .is_true(c1)) lines <- c(lines, "else {")
      else lines <- c(lines, sprintf("else if (%s) {", .cgpp(p, c1)))
      lines <- c(lines, .cgpp(p, ec$args[[1L]]), "}")
    }
    return(paste(lines, collapse = "\n"))
  }
  ecpairs <- vapply(e$args[-n], function(ec) sprintf("((%s) ? (\n%s\n)\n", .cgpp(p, ec$args[[2L]]), .cgpp(p, ec$args[[1L]])), "")
  last_line <- sprintf(": (\n%s\n)", .cgpp(p, last$args[[1L]]))
  paste0(paste(ecpairs, collapse = ": "), last_line, paste(strrep(")", length(ecpairs)), collapse = " "))
}

.cgp_c_print_float <- function(p, e) {
  type_ <- .cgd_get(p$type_aliases, real, real)
  .cgp_update(p, "macros", .cgd_get(p$type_macros, type_, character()))
  suffix <- .cgp_c_literal_suffix(p, type_)
  num <- .str_top(N(e, .int_val(.ft_attr(type_, "decimal_dig"))))
  if (!grepl("e", num, fixed = TRUE) && !grepl(".", num, fixed = TRUE)) num <- paste0(num, ".0")
  parts <- strsplit(num, "e", fixed = TRUE)[[1L]]
  parts[1L] <- sub("0+$", "", parts[1L])
  if (endsWith(parts[1L], ".")) parts[1L] <- paste0(parts[1L], "0")
  paste0(paste(parts, collapse = "e"), suffix)
}

.cgp_c_print_type <- function(p, t) {
  .cgp_update(p, "headers", .cgd_get(p$type_headers, t, character()))
  .cgp_update(p, "macros", .cgd_get(p$type_macros, t, character()))
  m <- .cgd_get(p$type_mappings, t)
  .cgpp(p, if (is.null(m)) t$args[[1L]] else m)
}

.cg_in_attrs <- function(v, a) any(vapply(.cg_vals(v)$attrs$args, function(z) identical(z$key, a$key), TRUE))

.cgp_c_print_declaration <- function(p, decl) {
  var <- decl$args[[1L]]
  vv <- .cg_vals(var)
  val <- vv$value
  if (identical(vv$type$key, untyped$key)) .cg_err("ValueError", "C does not support untyped variables")
  if (.cg_isa(var, "Pointer")) {
    result <- sprintf("%s%s *%s %s%s", if (.cg_in_attrs(var, value_const)) "const " else "", .cgpp(p, vv$type),
                      if (.cg_in_attrs(var, pointer_const)) " const" else "",
                      if (.cg_in_attrs(var, restrict)) "restrict " else "", .cgpp(p, vv$symbol))
  } else if (.cg_isa(var, "Variable")) {
    result <- sprintf("%s%s %s", if (.cg_in_attrs(var, value_const)) "const " else "", .cgpp(p, vv$type), .cgpp(p, vv$symbol))
  } else .cg_err("NotImplementedError", sprintf("Unknown type of var: %s", .cgp_type_str(var)))
  if (!.is_none(val)) result <- paste0(result, " = ", .cgpp(p, val))
  result
}

.cgp_c_print_element <- function(p, elem) {
  v <- .cg_vals(elem)
  if (.is_none(v$strides)) {
    if (!.is_none(v$offset)) .cg_err("ValueError", "Expected strides when offset is given")
    idxs <- paste(vapply(v$indices$args, function(a) .cgpp(p, a), ""), collapse = "][")
  } else {
    gi <- .add(Map(function(i, s) .mul(list(i, s)), v$indices$args, v$strides$args))
    if (!.is_none(v$offset)) gi <- .add(list(gi, v$offset))
    idxs <- .cgpp(p, gi)
  }
  sprintf("%s[%s]", .cgpp(p, v$symbol), idxs)
}

.cgp_c_codeblock <- function(p, e) paste(vapply(e$args, function(i) .cgp_call(p, "_get_statement", .cgpp(p, i)), ""), collapse = "\n")

.cgp_c_print_math_func <- function(p, e, nest = FALSE, known = NULL) {
  .cgp_requires(p, headers = "math.h", libraries = "m")
  .cgp_macro_or(p, e, function() {
    if (is.null(known)) known <- p$known_functions[[.cgp_clsname(e)]]
    if (!is.character(known) && !is.function(known)) {
      hit <- NULL
      for (cb in known) if (isTRUE(.cgp_truth(do.call(cb[[1L]], e$args)))) { hit <- cb[[2L]]; break }
      if (is.null(hit)) .cg_err("ValueError", "No matching printer")
      known <- hit
    }
    if (is.function(known)) return(do.call(known, c(list(p), e$args)))
    ns <- .cgp_a(p, "_ns")
    suffix <- if (paste0(ns, known) %in% .cgp_a(p, "_prec_funcs")) .cgp_c_func_suffix(p, real) else ""
    if (nest) {
      args <- .cgpp(p, e$args[[1L]])
      paren_pile <- ""
      if (length(e$args) > 1L) {
        for (curr in e$args[-c(1L, length(e$args))]) {
          paren_pile <- paste0(paren_pile, ")")
          args <- paste0(args, sprintf(", %s%s%s(%s", ns, known, suffix, .cgpp(p, curr)))
        }
        last <- e$args[[length(e$args)]]
        args <- paste0(args, sprintf(", %s%s", .cgpp(p, fcall(e$name, list(last))), paren_pile))
      }
    } else args <- paste(vapply(e$args, function(a) .cgpp(p, a), ""), collapse = ", ")
    sprintf("%s%s%s(%s)", ns, known, suffix, args)
  })
}

.cgp_c_minmax <- function(p, e, op) {
  if (.cgp_fname(e) %in% names(p$known_functions)) return(.cgp_print_function(p, e))
  inner <- function(args) {
    if (length(args) == 1L) return(.cgpp(p, args[[1L]]))
    half <- length(args) %/% 2L
    a <- inner(args[seq_len(half)]); b <- inner(args[-seq_len(half)])
    sprintf("((%s %s %s) ? %s : %s)", a, op, b, a, b)
  }
  inner(e$args)
}

.cgp_def("C89CodePrinter", "CodePrinter",
         attrs = list(printmethod = "_ccode", language = "C", standard = "C89", module = "sympy.printing.c",
                      reserved_words = .cg_c_reserved,
                      `_default_settings` = list(order = NULL, full_prec = "auto", error_on_reserved = FALSE,
                                                 reserved_word_suffix = "_", human = TRUE, inline = FALSE,
                                                 allow_unknown_functions = FALSE, strict = NULL, precision = 17,
                                                 user_functions = list(), contract = TRUE, dereference = list()),
                      `_ns` = "", `_kf` = .cg_kf_C89, math_macros = NULL),
         init = function(p, settings, nxt) {
           base_attr <- function(nm) .cgp_a(p, paste0(".", nm))
           if (is.null(.cgp_a(p, "math_macros"))) {
             p$math_macros <- if ("math_macros" %in% names(settings)) .cgd_as(settings$math_macros) else get_math_macros()
             settings$math_macros <- NULL
           }
           for (nm in c("type_aliases", "type_mappings", "type_headers", "type_macros", "type_func_suffixes",
                        "type_literal_suffixes", "type_math_macro_suffixes")) {
             assign(nm, .cgd_merge(base_attr(nm), settings[[nm]]), envir = p)
             settings[[nm]] <- NULL
           }
           nxt(settings)
           uf <- settings$user_functions
           kf <- .cgp_a(p, "_kf")
           for (k in names(uf)) kf[k] <- list(uf[[k]])
           p$known_functions <- kf
           p$dereference_set <- settings$dereference
           p$headers <- .cg_pyset(); p$libraries <- .cg_pyset(); p$macros <- .cg_pyset()
           invisible(NULL)
         },
         methods = list(
           `_rate_index_position` = function(p, pos) pos * 5L,
           `_get_statement` = function(p, s) if (endsWith(s, ";")) s else paste0(s, ";"),
           `_get_comment` = function(p, text) sprintf("/* %s */", text),
           `_declare_number_const` = function(p, name, value) {
             type_ <- .cgd_get(p$type_aliases, real)
             var <- Variable(name, type = type_, value = N(value, .int_val(.ft_attr(type_, "decimal_dig"))), attrs = list(value_const))
             .cgp_call(p, "_get_statement", .cgpp(p, Declaration(var)))
           },
           `_format_code` = function(p, lines) .cgp_call(p, "indent_code", lines),
           `_traverse_matrix_indices` = function(p, mat) {
             sh <- .cgp_mshape(mat)
             out <- list()
             for (i in seq_len(sh[1L]) - 1L) for (j in seq_len(sh[2L]) - 1L) out[[length(out) + 1L]] <- c(i, j)
             out
           },
           `_print_Mul` = function(p, e, ...) .cgp_macro_or(p, e, function() .cgp_super(p, "C89CodePrinter", "_print_Mul", e)),
           `_print_Pow` = function(p, e) .cgp_macro_or(p, e, function() {
             if ("Pow" %in% names(p$known_functions)) return(.cgp_print_function(p, e))
             PREC <- .cgp_prec(e)
             suffix <- .cgp_c_func_suffix(p, real)
             ns <- .cgp_a(p, "_ns")
             b <- e$args[[1L]]; x <- e$args[[2L]]
             if (.cgp_equal_valued(x, -1)) return(sprintf("%s/%s", .cgp_call(p, "_print_Float", Float(1.0)), .cgp_paren(p, b, PREC)))
             if (.cgp_equal_valued(x, 0.5)) return(sprintf("%ssqrt%s(%s)", ns, suffix, .cgpp(p, b)))
             if (x$type == "Rational" && x$value == gmp::as.bigq(1, 3) && .cgp_a(p, "standard") != "C89")
               return(sprintf("%scbrt%s(%s)", ns, suffix, .cgpp(p, b)))
             sprintf("%spow%s(%s, %s)", ns, suffix, .cgpp(p, b), .cgpp(p, x))
           }),
           `_print_Mod` = function(p, e) {
             num <- e$args[[1L]]; den <- e$args[[2L]]
             if (isTRUE(ask_fact(num, "integer")) && isTRUE(ask_fact(den, "integer"))) {
               PREC <- .cgp_prec(e)
               snum <- .cgp_paren(p, num, PREC); sden <- .cgp_paren(p, den, PREC)
               if ((isTRUE(ask_fact(num, "nonnegative")) && isTRUE(ask_fact(den, "nonnegative"))) ||
                   (isTRUE(ask_fact(num, "nonpositive")) && isTRUE(ask_fact(den, "nonpositive"))))
                 return(sprintf("%s %% %s", snum, sden))
               return(sprintf("((%s %% %s) + %s) %% %s", snum, sden, sden, sden))
             }
             .cgp_call(p, "_print_math_func", e, known = "fmod")
           },
           `_print_Rational` = function(p, e) {
             suffix <- .cgp_c_literal_suffix(p, real)
             sprintf("%s.0%s/%s.0%s", as.character(gmp::numerator(e$value)), suffix, as.character(gmp::denominator(e$value)), suffix)
           },
           `_print_Indexed` = function(p, e) .cgp_c_print_indexed(p, e),
           `_print_NumberSymbol` = function(p, e) .cgp_macro_or(p, e, function() .cgp_super(p, "C89CodePrinter", "_print_NumberSymbol", e)),
           `_print_Infinity` = function(p, e) "HUGE_VAL",
           `_print_NegativeInfinity` = function(p, e) "-HUGE_VAL",
           `_print_Piecewise` = function(p, e) .cgp_c_print_piecewise(p, e),
           `_print_ITE` = function(p, e) .cgpp(p, .cg_ite_piecewise(e)),
           `_print_MatrixElement` = function(p, e) {
             nc <- .cgp_mshape(e$args[[1L]])[2L]
             sprintf("%s[%s]", .cgp_paren(p, e$args[[1L]], .CG_PRECEDENCE[["Atom"]], strict = TRUE),
                     sstr(.add(list(e$args[[3L]], .mul(list(e$args[[2L]], Integer(nc)))))))
           },
           `_print_Symbol` = function(p, e) {
             name <- .cgp_super(p, "C89CodePrinter", "_print_Symbol", e)
             if (any(vapply(.cgp_s(p, "dereference"), function(d) is_sym(d) && identical(d$key, e$key), TRUE))) sprintf("(*%s)", name)
             else name
           },
           `_print_Relational` = function(p, e) sprintf("%s %s %s", .cgpp(p, e$args[[1L]]), .REL_OPS[[e$type]], .cgpp(p, e$args[[2L]])),
           `_print_For` = function(p, e) {
             v <- .cg_vals(e)
             target <- .cgpp(p, v$target)
             if (!(is_sym(v$iterable) && v$iterable$type == "Range"))
               .cg_err("NotImplementedError", "Only iterable currently supported is Range")
             r <- lapply(v$iterable$args, sstr)
             body <- .cgpp(p, v$body)
             sprintf("for (%s = %s; %s < %s; %s += %s) {\n%s\n}", target, r[[1L]], target, r[[2L]], target, r[[3L]], body)
           },
           `_print_sign` = function(p, e) sprintf("(((%s) > 0) - ((%s) < 0))", .cgpp(p, e$args[[1L]]), .cgpp(p, e$args[[1L]])),
           `_print_Max` = function(p, e) .cgp_c_minmax(p, e, ">"),
           `_print_Min` = function(p, e) .cgp_c_minmax(p, e, "<"),
           indent_code = function(p, code) .cgp_indent_braces(code),
           `_print_Tuple` = function(p, e) paste0("{", paste(vapply(e$args, function(a) .cgpp(p, a), ""), collapse = ", "), "}"),
           `_print_List` = function(p, e) paste0("{", paste(vapply(e$args, function(a) .cgpp(p, a), ""), collapse = ", "), "}"),
           `_print_Type` = function(p, e) .cgp_c_print_type(p, e),
           `_print_Declaration` = function(p, e) .cgp_c_print_declaration(p, e),
           `_print_Float` = function(p, e) .cgp_c_print_float(p, e),
           `_print_BooleanTrue` = function(p, e) { .cgp_requires(p, headers = "stdbool.h"); "true" },
           `_print_BooleanFalse` = function(p, e) { .cgp_requires(p, headers = "stdbool.h"); "false" },
           `_print_Element` = function(p, e) .cgp_c_print_element(p, e),
           `_print_CodeBlock` = function(p, e) .cgp_c_codeblock(p, e),
           `_print_While` = function(p, e) {
             v <- .cg_vals(e)
             sprintf("while (%s) {\n%s\n}", .cgpp(p, v$condition), .cgpp(p, v$body))
           },
           `_print_Scope` = function(p, e) sprintf("{\n%s\n}", .cgp_c_codeblock(p, .cg_vals(e)$body)),
           `_print_Print` = function(p, e) {
             .cgp_requires(p, headers = "stdio.h")
             v <- .cg_vals(e)
             template <- if (.is_none(v$file)) "printf(%s, %s)" else paste0("fprintf(", .cgpp(p, v$file), ", %s, %s)")
             fmt <- if (.is_none(v$format_string)) "%s\n" else .cgpp(p, v$format_string)
             sprintf(template, fmt, paste(vapply(v$print_args$args, function(a) .cgpp(p, a), ""), collapse = ", "))
           },
           `_print_Stream` = function(p, e) .cg_text(.cg_vals(e)$name),
           `_print_FunctionPrototype` = function(p, e) {
             v <- .cg_vals(e)
             pars <- paste(vapply(v$parameters$args, function(a) .cgpp(p, Declaration(a)), ""), collapse = ", ")
             sprintf("%s %s(%s)", .cgpp(p, v$return_type), .cgpp(p, v$name), pars)
           },
           `_print_FunctionDefinition` = function(p, e)
             paste0(.cgp_call(p, "_print_FunctionPrototype", e), .cgp_call(p, "_print_Scope", e)),
           `_print_Return` = function(p, e) paste0("return ", .cgpp(p, e$args[[1L]])),
           `_print_CommaOperator` = function(p, e) sprintf("(%s)", paste(vapply(e$args, function(a) .cgpp(p, a), ""), collapse = ", ")),
           `_print_Label` = function(p, e) {
             v <- .cg_vals(e)
             nm <- .cg_text(v$name)
             if (.is_none(v$body)) return(paste0(nm, ":"))
             if (length(v$body$args) == 1L) return(sprintf("%s:\n%s", nm, .cgp_c_codeblock(p, v$body)))
             sprintf("%s:\n{\n%s\n}", nm, .cgp_c_codeblock(p, v$body))
           },
           `_print_goto` = function(p, e) paste0("goto ", .cg_text(.cg_vals(.cg_vals(e)$label)$name)),
           `_print_PreIncrement` = function(p, e) sprintf("++(%s)", .cgpp(p, e$args[[1L]])),
           `_print_PostIncrement` = function(p, e) sprintf("(%s)++", .cgpp(p, e$args[[1L]])),
           `_print_PreDecrement` = function(p, e) sprintf("--(%s)", .cgpp(p, e$args[[1L]])),
           `_print_PostDecrement` = function(p, e) sprintf("(%s)--", .cgpp(p, e$args[[1L]])),
           `_print_struct` = function(p, e) .cgp_c_struct(p, e),
           `_print_union` = function(p, e) .cgp_c_struct(p, e),
           `_print_BreakToken` = function(p, e) "break",
           `_print_ContinueToken` = function(p, e) "continue"
         ))

.cgp_c_struct <- function(p, e) {
  v <- .cg_vals(e)
  lines <- c(vapply(v$declarations$args, function(d) .cgpp(p, d), ""), "")
  sprintf("%s %s {\n%s}", e$type, sstr(v$name), paste(lines, collapse = ";\n"))
}

# class attributes of C89CodePrinter (dictionaries keyed by types)
.cgp_classes$C89CodePrinter$attrs$.type_aliases <- .cgd(list(list(real, float64), list(complex_, complex128), list(sym_integer, intc)))
.cgp_classes$C89CodePrinter$attrs$.type_mappings <- .cgd(list(
  list(real, "double"), list(intc, "int"), list(float32, "float"), list(float64, "double"), list(sym_integer, "int"),
  list(bool_, "bool"), list(int8, "int8_t"), list(int16, "int16_t"), list(int32, "int32_t"), list(int64, "int64_t"),
  list(uint8, "int8_t"), list(uint16, "int16_t"), list(uint32, "int32_t"), list(uint64, "int64_t")))
.cgp_classes$C89CodePrinter$attrs$.type_headers <- .cgd(c(list(list(bool_, "stdbool.h")),
  lapply(list(int8, int16, int32, int64, uint8, uint16, uint32, uint64), function(t) list(t, "stdint.h"))))
.cgp_classes$C89CodePrinter$attrs$.type_macros <- .cgd()
.cgp_classes$C89CodePrinter$attrs$.type_func_suffixes <- .cgd(list(list(float32, "f"), list(float64, ""), list(float80, "l")))
.cgp_classes$C89CodePrinter$attrs$.type_literal_suffixes <- .cgd(list(list(float32, "F"), list(float64, ""), list(float80, "L")))
.cgp_classes$C89CodePrinter$attrs$.type_math_macro_suffixes <- .cgd(list(list(float80, "l")))

.cg_c99_prec_funcs <- strsplit(paste("fabs fmod remainder remquo fma fmax fmin fdim nan exp exp2 expm1 log log10 log2",
                                     "log1p pow sqrt cbrt hypot sin cos tan asin acos atan atan2 sinh cosh tanh asinh",
                                     "acosh atanh erf erfc tgamma lgamma ceil floor trunc round nearbyint rint frexp",
                                     "ldexp modf scalbn ilogb logb nextafter copysign"), " ")[[1L]]

.cgp_def("C99CodePrinter", "C89CodePrinter",
         attrs = list(standard = "C99", reserved_words = c(.cg_c_reserved, .cg_c99_reserved), `_kf` = .cg_kf_C99,
                      `_prec_funcs` = .cg_c99_prec_funcs,
                      .type_mappings = .cgd_merge(.cgp_classes$C89CodePrinter$attrs$.type_mappings,
                                                  .cgd(list(list(complex64, "float complex"), list(complex128, "double complex")))),
                      .type_headers = .cgd_merge(.cgp_classes$C89CodePrinter$attrs$.type_headers,
                                                 .cgd(list(list(complex64, "complex.h"), list(complex128, "complex.h"))))),
         methods = list(
           `_print_Infinity` = function(p, e) "INFINITY",
           `_print_NegativeInfinity` = function(p, e) "-INFINITY",
           `_print_NaN` = function(p, e) "NAN",
           `_print_math_func` = function(p, e, nest = FALSE, known = NULL) .cgp_c_print_math_func(p, e, nest, known),
           `_print_Max` = function(p, e) .cgp_call(p, "_print_math_func", e, nest = TRUE),
           `_print_Min` = function(p, e) .cgp_call(p, "_print_math_func", e, nest = TRUE),
           `_get_loop_opening_ending` = function(p, indices) {
             open <- vapply(indices, function(i) sprintf("for (int %s=%s; %s<%s; %s++){", .cgpp(p, .cg_idx_label(i)),
                                                         .cgpp(p, .cgp_idx_lower(i)), .cgpp(p, .cg_idx_label(i)),
                                                         .cgpp(p, .add(list(.cgp_idx_upper(i), S_One))), .cgpp(p, .cg_idx_label(i))), "")
             list(open, rep("}", length(indices)))
           }
         ))
for (.cg_k in strsplit("Abs Sqrt exp exp2 expm1 log log10 log2 log1p Cbrt hypot fma loggamma sin cos tan asin acos atan atan2 sinh cosh tanh asinh acosh atanh erf erfc loggamma gamma ceiling floor", " ")[[1L]])
  assign(paste0("_print_", .cg_k), function(p, e) .cgp_call(p, "_print_math_func", e), envir = .cgp_classes$C99CodePrinter$methods)
rm(.cg_k)

.cgp_def("C11CodePrinter", "C99CodePrinter",
         methods = list(`_print_alignof` = function(p, e) {
           .cgp_requires(p, headers = "stdalign.h")
           sprintf("alignof(%s)", .cgpp(p, e$args[[1L]]))
         }))

# ---------------------------------------------------------------------------
# C++

.cg_cxx_reserved98 <- c("and", "and_eq", "asm", "auto", "bitand", "bitor", "bool", "break", "case", "catch,", "char",
                        "class", "compl", "const", "const_cast", "continue", "default", "delete", "do", "double",
                        "dynamic_cast", "else", "enum", "explicit", "export", "extern", "false", "float", "for",
                        "friend", "goto", "if", "inline", "int", "long", "mutable", "namespace", "new", "not",
                        "not_eq", "operator", "or", "or_eq", "private", "protected", "public", "register",
                        "reinterpret_cast", "return", "short", "signed", "sizeof", "static", "static_cast",
                        "struct", "switch", "template", "this", "throw", "true", "try", "typedef", "typeid",
                        "typename", "union", "unsigned", "using", "virtual", "void", "volatile", "wchar_t", "while",
                        "xor", "xor_eq")
.cg_cxx_reserved11 <- c(.cg_cxx_reserved98, "alignas", "alignof", "char16_t", "char32_t", "constexpr", "decltype",
                        "noexcept", "nullptr", "static_assert", "thread_local")
.cg_cxx_reserved17 <- setdiff(.cg_cxx_reserved11, "register")

.cg_cxx_math <- list(`C++98` = c(Mod = "fmod", ceiling = "ceil"), `C++11` = c(gamma = "tgamma"),
                     `C++17` = c(beta = "beta", Ei = "expint", zeta = "riemann_zeta"))
for (.cg_k in c("Abs", "exp", "log", "log10", "sqrt", "sin", "cos", "tan", "asin", "acos", "atan", "atan2", "sinh", "cosh", "tanh", "floor"))
  .cg_cxx_math$`C++98`[.cg_k] <- tolower(.cg_k)
for (.cg_k in c("asinh", "acosh", "atanh", "erf", "erfc")) .cg_cxx_math$`C++11`[.cg_k] <- tolower(.cg_k)
rm(.cg_k)

.cgp_cxx_minmax <- function(p, e, fn) {
  .cgp_requires(p, headers = "algorithm")
  if (length(e$args) == 1L) return(.cgpp(p, e$args[[1L]]))
  sprintf("%s%s(%s, %s)", .cgp_a(p, "_ns"), tolower(fn), .cgpp(p, e$args[[1L]]), .cgpp(p, do.call(if (fn == "Max") Max else Min, e$args[-1L])))
}

.cgp_def("_CXXCodePrinterBase",
         attrs = list(printmethod = "_cxxcode", language = "C++", `_ns` = "std::", module = "sympy.printing.cxx"),
         methods = list(
           `_print_Max` = function(p, e) .cgp_cxx_minmax(p, e, "Max"),
           `_print_Min` = function(p, e) .cgp_cxx_minmax(p, e, "Min"),
           `_print_using` = function(p, e) {
             v <- .cg_vals(e)
             if (.is_none(v$alias)) return(paste0("using ", sstr(v$type)))
             .cg_err("ValueError", "C++98 does not support type aliases")
           },
           `_print_Raise` = function(p, e) paste0("throw ", .cgpp(p, e$args[[1L]])),
           `_print_RuntimeError_` = function(p, e) {
             .cgp_requires(p, headers = "stdexcept")
             sprintf("%sruntime_error(%s)", .cgp_a(p, "_ns"), .cgpp(p, e$args[[1L]]))
           }))

.cgp_attach_cxx <- function(cls, std) {
  env <- .cgp_classes[[cls]]$methods
  for (nm in names(.cg_cxx_math[[std]])) {
    local({
      fname <- .cg_cxx_math[[std]][[nm]]
      assign(paste0("_print_", nm), function(p, e)
        sprintf("%s%s(%s)", .cgp_a(p, "_ns"), fname, paste(vapply(e$args, function(a) .cgpp(p, a), ""), collapse = ", ")),
        envir = env)
    })
  }
}

.cgp_def("CXX98CodePrinter", c("_CXXCodePrinterBase", "C89CodePrinter"),
         attrs = list(standard = "C++98", reserved_words = .cg_cxx_reserved98))
.cgp_attach_cxx("CXX98CodePrinter", "C++98")

.cgp_def("CXX11CodePrinter", c("_CXXCodePrinterBase", "C99CodePrinter"),
         attrs = list(standard = "C++11", reserved_words = .cg_cxx_reserved11,
                      .type_mappings = .cgd_merge(.cgp_classes$C89CodePrinter$attrs$.type_mappings, .cgd(list(
                        list(Type("int8"), list("int8_t", "cstdint")), list(Type("int16"), list("int16_t", "cstdint")),
                        list(Type("int32"), list("int32_t", "cstdint")), list(Type("int64"), list("int64_t", "cstdint")),
                        list(Type("uint8"), list("uint8_t", "cstdint")), list(Type("uint16"), list("uint16_t", "cstdint")),
                        list(Type("uint32"), list("uint32_t", "cstdint")), list(Type("uint64"), list("uint64_t", "cstdint")),
                        list(Type("complex64"), list("std::complex<float>", "complex")),
                        list(Type("complex128"), list("std::complex<double>", "complex")),
                        list(Type("bool"), list("bool", NULL)))))),
         methods = list(`_print_using` = function(p, e) {
           v <- .cg_vals(e)
           if (.is_none(v$alias)) return(.cgp_super(p, "CXX11CodePrinter", "_print_using", e))
           sprintf("using %s = %s", .cgpp(p, v$alias), .cgpp(p, v$type))
         }))
.cgp_attach_cxx("CXX11CodePrinter", "C++11")

.cgp_def("CXX17CodePrinter", c("_CXXCodePrinterBase", "C99CodePrinter"),
         attrs = list(standard = "C++17", reserved_words = .cg_cxx_reserved17,
                      `_kf` = c(.cg_kf_C99, as.list(.cg_cxx_math$`C++17`))),
         methods = list(`_print_beta` = function(p, e) .cgp_call(p, "_print_math_func", e),
                        `_print_Ei` = function(p, e) .cgp_call(p, "_print_math_func", e),
                        `_print_zeta` = function(p, e) .cgp_call(p, "_print_math_func", e)))
.cgp_attach_cxx("CXX17CodePrinter", "C++17")

#' @rdname code-printer-objects
#' @export
C89CodePrinter <- function(settings = NULL) .cgp_new("C89CodePrinter", settings)
#' @rdname code-printer-objects
#' @export
C99CodePrinter <- function(settings = NULL) .cgp_new("C99CodePrinter", settings)
#' @rdname code-printer-objects
#' @export
C11CodePrinter <- function(settings = NULL) .cgp_new("C11CodePrinter", settings)
#' @rdname code-printer-objects
#' @export
CXX98CodePrinter <- function(settings = NULL) .cgp_new("CXX98CodePrinter", settings)
#' @rdname code-printer-objects
#' @export
CXX11CodePrinter <- function(settings = NULL) .cgp_new("CXX11CodePrinter", settings)
#' @rdname code-printer-objects
#' @export
CXX17CodePrinter <- function(settings = NULL) .cgp_new("CXX17CodePrinter", settings)

.cg_c_printers <- c(c89 = "C89CodePrinter", c99 = "C99CodePrinter", c11 = "C11CodePrinter")
.cg_cxx_printers <- c(`c++98` = "CXX98CodePrinter", `c++11` = "CXX11CodePrinter", `c++17` = "CXX17CodePrinter")

.cg_ccode <- function(e, assign_to = NULL, standard = "c99", ...) {
  cls <- .cg_c_printers[[tolower(standard)]]
  .cgp_doprint(.cgp_new(cls, list(...)), e, assign_to)
}
.cg_cxxcode <- function(e, assign_to = NULL, standard = "c++11", ...) {
  cls <- .cg_cxx_printers[[tolower(standard)]]
  .cgp_doprint(.cgp_new(cls, list(...)), e, assign_to)
}
