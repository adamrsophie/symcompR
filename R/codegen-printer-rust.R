# Rust code printer, after sympy/printing/rust.py.

def_function("float_floor", facts = function(args, afs) c(integer = FALSE))
def_function("float_ceiling", facts = function(args, afs) c(integer = FALSE))

.cg_rust_kf <- list(
  float_floor = "floor", float_ceiling = "ceil", Abs = "abs",
  Pow = list(list(function(base, exp) .cgp_equal_valued(exp, -1), "recip", 2L),
             list(function(base, exp) .cgp_equal_valued(exp, 0.5), "sqrt", 2L),
             list(function(base, exp) .cgp_equal_valued(exp, -0.5), "sqrt().recip", 2L),
             list(function(base, exp) identical(exp$key, Rational(1L, 3L)$key), "cbrt", 2L),
             list(function(base, exp) .cgp_equal_valued(base, 2), "exp2", 3L),
             list(function(base, exp) ask_fact(exp, "integer"), "powi", 1L),
             list(function(base, exp) .py_not(ask_fact(exp, "integer")), "powf", 1L)),
  exp = list(list(function(exp) TRUE, "exp", 2L)),
  log = "ln", Max = "max", Min = "min", sin = "sin", cos = "cos", tan = "tan", asin = "asin", acos = "acos",
  atan = "atan", atan2 = "atan2", sinh = "sinh", cosh = "cosh", tanh = "tanh", asinh = "asinh", acosh = "acosh",
  atanh = "atanh", sqrt = "sqrt")
.cg_rust_overrides <- list(floor = "float_floor", ceiling = "float_ceiling")
.cg_rust_reserved <- c("abstract", "as", "async", "await", "become", "box", "break", "const", "continue", "crate", "do",
                       "dyn", "else", "enum", "extern", "false", "final", "fn", "for", "gen", "if", "impl", "in", "let",
                       "loop", "macro", "match", "mod", "move", "mut", "override", "priv", "pub", "ref", "return", "Self",
                       "self", "static", "struct", "super", "trait", "true", "try", "type", "typeof", "unsafe", "unsized",
                       "use", "virtual", "where", "while", "yield")

# Python's `not v` for a three valued fact
.py_not <- function(v) !isTRUE(v)

# TypeCast(expr, type_): the expression printed with an `as` cast when an
# integer expression meets a float
.cg_typecast <- function(e, type_) {
  explicit <- isTRUE(ask_fact(e, "integer")) && !identical(type_$key, sym_integer$key)
  new_sym("TypeCast", list(e, type_), extra = list(explicit = explicit))
}
.cgp_rust_cast <- function(e) {
  if (!.expr_is_number(e)) return(.cg_typecast(e, real))
  if (isTRUE(ask_fact(e, "integer"))) return(Float(N(e)))
  e
}
.cgp_rust_floats <- function(args) any(vapply(args, function(a) isTRUE(ask_fact(a, "real")) && .py_not(ask_fact(a, "integer")), TRUE))
.cgp_rust_type <- function(p, t) .cgd_get(p$type_mappings, .cgd_get(p$type_aliases, t))

.cgp_rust_caller_var <- function(p, e) {
  if (length(e$args) > 1L) return(paste0("(", .cgpp(p, e), ")"))
  if (.expr_is_number(e)) {
    if (e$type %in% c("Integer", "Float", "Exp1", "Pi", "Infinity", "NegativeInfinity", "NaN")) return(.cgpp(p, e, `_type` = TRUE))
    return(.cgpp(p, e))
  }
  .cgpp(p, e)
}

.cgp_rust_function <- function(p, e) {
  nm <- .cgp_fname(e)
  kf <- p$known_functions
  imp <- if (e$type == "Function") e$extra$imp else NULL
  if (nm %in% names(kf)) {
    cond_func <- kf[[nm]]
    func <- NULL; style <- 1L
    if (is.character(cond_func)) func <- cond_func
    else for (cf in cond_func) {
      func <- cf[[2L]]; style <- if (length(cf) > 2L) as.integer(cf[[3L]]) else 1L
      if (isTRUE(.cgp_truth(do.call(cf[[1L]], e$args)))) break
    }
    if (!is.null(func)) {
      a <- e$args
      return(switch(as.character(style),
                    `1` = sprintf("%s.%s(%s)", .cgp_rust_caller_var(p, a[[1L]]), func,
                                  if (length(a) > 1L) .cgp_stringify(p, a[-1L], ", ") else ""),
                    `2` = sprintf("%s.%s()", .cgp_rust_caller_var(p, a[[1L]]), func),
                    `3` = sprintf("%s.%s()", .cgp_rust_caller_var(p, a[[2L]]), func),
                    sprintf("%s(%s)", func, .cgp_stringify(p, a, ", "))))
    }
    return(NULL)
  }
  if (!is.null(imp) && is_sym(imp) && imp$type == "Lambda") return(.cgpp(p, .cg_apply_lambda(imp, e$args)))
  .cgp_not_supported(p, e)
}

.cgp_rust_can_print <- function(p, name)
  name %in% names(p$known_functions) || name %in% names(.cg_rust_overrides) || is.function(.cgp_get(p, paste0("_print_", name)))

.cgp_rust_rewrite_known <- function(p, e) {
  if (!is_sym(e)) return(e)
  fns <- character()
  walk <- function(z) {
    if (!is_sym(z)) return()
    if (.cgp_is_function(z)) fns <<- union(fns, .cgp_fname(z))
    for (a in z$args) walk(a)
  }
  walk(e)
  for (f in fns) {
    rw <- .CG_REWRITEABLE[[f]]
    if (is.null(rw)) next
    if (.cgp_rust_can_print(p, rw[[1L]]) && all(vapply(rw[[2L]], function(g) .cgp_rust_can_print(p, g), TRUE)))
      e <- rewrite(e, rw[[1L]])
  }
  e
}

.cgp_rust_override <- function(e) {
  if (!is_sym(e) || !length(e$args)) return(e)
  args <- lapply(e$args, .cgp_rust_override)
  if (e$type == "Function" && e$name %in% names(.cg_rust_overrides) && is.null(e$extra$imp))
    return(fcall(.cg_rust_overrides[[e$name]], args))
  sym_rebuild(e, args)
}

.cgp_def("RustCodePrinter", "CodePrinter",
         attrs = list(printmethod = "_rust_code", language = "Rust", module = "sympy.printing.rust",
                      `_default_settings` = list(order = NULL, full_prec = "auto", error_on_reserved = FALSE,
                                                 reserved_word_suffix = "_", human = TRUE, inline = FALSE,
                                                 allow_unknown_functions = FALSE, strict = NULL, precision = 17,
                                                 user_functions = list(), contract = TRUE, dereference = list())),
         init = function(p, settings, nxt) {
           nxt(settings)
           kf <- .cg_rust_kf
           for (k in names(settings$user_functions)) kf[k] <- list(settings$user_functions[[k]])
           p$known_functions <- kf
           p$dereference_set <- settings$dereference
           p$reserved_words <- .cg_rust_reserved
           p$type_aliases <- .cgd(list(list(sym_integer, int32), list(real, float64)))
           p$type_mappings <- .cgd(list(list(int32, "i32"), list(float32, "f32"), list(float64, "f64"), list(bool_, "bool")))
           invisible(NULL)
         },
         methods = list(
           `_rate_index_position` = function(p, pos) pos * 5L,
           `_get_statement` = function(p, s) paste0(s, ";"),
           `_get_comment` = function(p, text) paste0("// ", text),
           `_declare_number_const` = function(p, name, value)
             sprintf("const %s: %s = %s;", .py_str_of(name), .cgp_rust_type(p, real), .py_str_of(value)),
           `_format_code` = function(p, lines) .cgp_call(p, "indent_code", lines),
           `_traverse_matrix_indices` = function(p, mat) {
             sh <- .cgp_mshape(mat)
             out <- list()
             for (i in seq_len(sh[1L]) - 1L) for (j in seq_len(sh[2L]) - 1L) out[[length(out) + 1L]] <- c(i, j)
             out
           },
           `_get_loop_opening_ending` = function(p, indices) {
             open <- vapply(indices, function(i) sprintf("for %s in %s..%s {", .cgpp(p, i), .cgpp(p, .cgp_idx_lower(i)),
                                                         .cgpp(p, .add(list(.cgp_idx_upper(i), S_One)))), "")
             list(open, rep("}", length(indices)))
           },
           `_print_Function` = function(p, e) {
             r <- .cgp_rust_function(p, e)
             if (is.null(r)) .cgp_not_supported(p, e) else r
           },
           `_print_Mul` = function(p, e) {
             if (.cgp_rust_floats(e$args)) {
               ord <- .cgp_s(p, "order")
               args <- if (identical(ord, "old") || identical(ord, "none")) e$args else .mul_print_order(e$args)
               args <- lapply(args, function(a) if (a$type == "Integer" && a$value == -1) a else .cgp_rust_cast(a))
               old <- p$settings$order
               p$settings$order <- "none"
               on.exit(p$settings$order <- old)
               return(.cgp_print_mul(p, new_sym("Mul", args)))
             }
             .cgp_print_mul(p, e)
           },
           `_print_Add` = function(p, e, order = NULL) {
             if (.cgp_rust_floats(e$args)) {
               terms <- lapply(.cgp_ordered_terms(p, e, order), .cgp_rust_cast)
               return(.cgp_super(p, "RustCodePrinter", "_print_Add", new_sym("Add", terms), order = "none"))
             }
             .cgp_super(p, "RustCodePrinter", "_print_Add", e, order = order)
           },
           `_print_Pow` = function(p, e) {
             b <- e$args[[1L]]; x <- e$args[[2L]]
             if (isTRUE(ask_fact(b, "integer")) && .py_not(ask_fact(x, "integer"))) return(.cgpp(p, .pow(Float(N(b)), x)))
             .cgp_rust_function(p, e) %||% .cgp_not_supported(p, e)
           },
           `_print_TypeCast` = function(p, e) {
             if (!isTRUE(e$extra$explicit)) return(.cgpp(p, e$args[[1L]]))
             paste0(.cgpp(p, e$args[[1L]]), " as ", .cgp_rust_type(p, e$args[[2L]]))
           },
           `_print_Float` = function(p, e, `_type` = FALSE) {
             r <- .cgp_super(p, "RustCodePrinter", "_print_Float", e)
             if (isTRUE(`_type`)) paste0(r, "_", .cgp_rust_type(p, real)) else r
           },
           `_print_Integer` = function(p, e, `_type` = FALSE) {
             r <- .cgp_super(p, "RustCodePrinter", "_print_Integer", e)
             if (isTRUE(`_type`)) paste0(r, "_", .cgp_rust_type(p, sym_integer)) else r
           },
           `_print_Rational` = function(p, e)
             sprintf("%s_%s/%s.0", as.character(gmp::numerator(e$value)), .cgp_rust_type(p, real), as.character(gmp::denominator(e$value))),
           `_print_Relational` = function(p, e) {
             l <- e$args[[1L]]; r <- e$args[[2L]]
             li <- isTRUE(ask_fact(l, "integer")); ri <- isTRUE(ask_fact(r, "integer"))
             if ((li && .py_not(ask_fact(r, "integer"))) || (ri && .py_not(ask_fact(l, "integer")))) {
               l <- .cgp_rust_cast(l); r <- .cgp_rust_cast(r)
             }
             sprintf("%s %s %s", .cgpp(p, l), .REL_OPS[[e$type]], .cgpp(p, r))
           },
           `_print_Indexed` = function(p, e) {
             dims <- .cgp_indexed_shape(e)
             idx <- e$args[-1L]
             elem <- S_Zero; offset <- S_One
             for (i in rev(seq_along(idx))) {
               elem <- .add(list(elem, .mul(list(idx[[i]], offset))))
               offset <- .mul(list(offset, dims[[i]]))
             }
             sprintf("%s[%s]", .cgpp(p, e$args[[1L]]$args[[1L]]), .cgpp(p, elem))
           },
           `_print_Idx` = function(p, e) e$name,
           `_print_Dummy` = function(p, e) e$name,
           `_print_Exp1` = function(p, e, `_type` = FALSE) "E",
           `_print_Pi` = function(p, e, `_type` = FALSE) "PI",
           `_print_Infinity` = function(p, e, `_type` = FALSE) "INFINITY",
           `_print_NegativeInfinity` = function(p, e, `_type` = FALSE) "NEG_INFINITY",
           `_print_BooleanTrue` = function(p, e, `_type` = FALSE) "true",
           `_print_BooleanFalse` = function(p, e, `_type` = FALSE) "false",
           `_print_bool` = function(p, e, `_type` = FALSE) tolower(as.character(e)),
           `_print_NaN` = function(p, e, `_type` = FALSE) "NAN",
           `_print_Piecewise` = function(p, e) {
             n <- length(e$args)
             if (!.is_true(e$args[[n]]$args[[2L]]))
               .cg_err("ValueError", paste0("All Piecewise expressions must contain an (expr, True) statement to be used as a default condition. ",
                                            "Without one, the generated expression may not evaluate to anything under some condition."))
             lines <- character()
             for (i in seq_len(n)) {
               c1 <- e$args[[i]]$args[[2L]]
               if (i == 1L) lines <- c(lines, sprintf("if (%s) {", .cgpp(p, c1)))
               else if (i == n && .is_true(c1)) lines[length(lines)] <- paste0(lines[length(lines)], " else {")
               else lines[length(lines)] <- paste0(lines[length(lines)], sprintf(" else if (%s) {", .cgpp(p, c1)))
               lines <- c(lines, .cgpp(p, e$args[[i]]$args[[1L]]), "}")
             }
             paste(lines, collapse = if (isTRUE(.cgp_s(p, "inline"))) " " else "\n")
           },
           `_print_ITE` = function(p, e) .cgpp(p, .cg_ite_piecewise(e)),
           `_print_MatrixBase` = function(p, A) {
             A <- .unwrap_explicit(A)
             if (A$ncol == 1L) return(sprintf("[%s]", paste(vapply(A$data, function(a) .cgpp(p, a), ""), collapse = ", ")))
             .cg_err("ValueError", "Full Matrix Support in Rust need Crates (https://crates.io/keywords/matrix).")
           },
           `_print_MatrixElement` = function(p, e) {
             nc <- .cgp_mshape(e$args[[1L]])[2L]
             sprintf("%s[%s]", sstr(e$args[[1L]]), sstr(.add(list(e$args[[3L]], .mul(list(e$args[[2L]], Integer(nc)))))))
           },
           `_print_Symbol` = function(p, e) {
             name <- .cgp_super(p, "RustCodePrinter", "_print_Symbol", e)
             if (any(vapply(p$dereference_set, function(d) is_sym(d) && identical(d$key, e$key), TRUE))) sprintf("(*%s)", name)
             else name
           },
           `_print_Assignment` = function(p, e) {
             lhs <- e$args[[1L]]; rhs <- e$args[[2L]]
             if (isTRUE(.cgp_s(p, "contract")) && (.cg_has_type(lhs, "IndexedBase") || .cg_has_type(rhs, "IndexedBase")))
               return(.cgp_doprint_loops(p, rhs, lhs))
             .cgp_call(p, "_get_statement", sprintf("%s = %s", .cgpp(p, lhs), .cgpp(p, rhs)))
           },
           `_print_sign` = function(p, e) {
             a <- .cgpp(p, e$args[[1L]])
             sprintf("(if (%s == 0.0) { 0.0 } else { (%s).signum() })", a, a)
           },
           indent_code = function(p, code) .cgp_indent_braces(code, "    ")
         ))

#' @rdname code-printer-objects
#' @export
RustCodePrinter <- function(settings = NULL) .cgp_new("RustCodePrinter", settings)

.cg_rust_code <- function(e, assign_to = NULL, ...) {
  p <- .cgp_new("RustCodePrinter", list(...))
  e <- .cgp_rust_rewrite_known(p, e)
  if (is_sym(e)) e <- .cgp_rust_override(e)
  .cgp_doprint(p, e, assign_to)
}
