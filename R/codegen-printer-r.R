# R code printer, after sympy/printing/rcode.py.

.cg_kf_r <- list(Abs = "abs", sin = "sin", cos = "cos", tan = "tan", asin = "asin", acos = "acos", atan = "atan",
                 atan2 = "atan2", exp = "exp", log = "log", erf = "erf", sinh = "sinh", cosh = "cosh", tanh = "tanh",
                 asinh = "asinh", acosh = "acosh", atanh = "atanh", floor = "floor", ceiling = "ceiling",
                 sign = "sign", Max = "max", Min = "min", factorial = "factorial", gamma = "gamma",
                 digamma = "digamma", trigamma = "trigamma", beta = "beta", sqrt = "sqrt")

.cg_r_reserved <- c("if", "else", "repeat", "while", "function", "for", "in", "next", "break", "TRUE", "FALSE",
                    "NULL", "Inf", "NaN", "NA", "NA_integer_", "NA_real_", "NA_complex_", "NA_character_",
                    "volatile")

.cgp_def("RCodePrinter", "CodePrinter",
         attrs = list(printmethod = "_rcode", language = "R", module = "sympy.printing.rcode",
                      reserved_words = .cg_r_reserved,
                      `_operators` = list(and = "&", or = "|", not = "!"),
                      `_relationals` = list(),
                      `_default_settings` = list(order = NULL, full_prec = "auto", error_on_reserved = FALSE,
                                                 reserved_word_suffix = "_", human = TRUE, inline = FALSE,
                                                 allow_unknown_functions = FALSE, strict = NULL, precision = 15,
                                                 user_functions = list(), contract = TRUE, dereference = list())),
         init = function(p, settings, nxt) { nxt(settings); .cgp_jo_known(p, .cg_kf_r, settings) },
         methods = list(
           `_rate_index_position` = function(p, pos) pos * 5L,
           `_get_statement` = function(p, s) paste0(s, ";"),
           `_get_comment` = function(p, text) paste0("// ", text),
           `_declare_number_const` = function(p, name, value) sprintf("%s = %s;", .py_str_of(name), .py_str_of(value)),
           `_format_code` = function(p, lines) .cgp_call(p, "indent_code", lines),
           `_traverse_matrix_indices` = function(p, mat) {
             sh <- .cgp_mshape(mat)
             out <- list()
             for (i in seq_len(sh[1L]) - 1L) for (j in seq_len(sh[2L]) - 1L) out[[length(out) + 1L]] <- c(i, j)
             out
           },
           `_get_loop_opening_ending` = function(p, indices) {
             open <- vapply(indices, function(i)
               sprintf("for (%s in %s:%s){", .cgpp(p, .cg_idx_label(i)), .cgpp(p, .add(list(.cgp_idx_lower(i), S_One))),
                       .cgpp(p, .add(list(.cgp_idx_upper(i), S_One)))), "")
             list(open, rep("}", length(indices)))
           },
           `_print_Pow` = function(p, e) {
             if ("Pow" %in% names(p$known_functions)) return(.cgp_print_function(p, e))
             PREC <- .cgp_prec(e)
             b <- e$args[[1L]]; x <- e$args[[2L]]
             if (.cgp_equal_valued(x, -1)) return(paste0("1.0/", .cgp_paren(p, b, PREC)))
             if (.cgp_equal_valued(x, 0.5)) return(sprintf("sqrt(%s)", .cgpp(p, b)))
             sprintf("%s^%s", .cgp_paren(p, b, PREC), .cgp_paren(p, x, PREC))
           },
           `_print_Rational` = function(p, e)
             sprintf("%s.0/%s.0", as.character(gmp::numerator(e$value)), as.character(gmp::denominator(e$value))),
           `_print_Indexed` = function(p, e)
             sprintf("%s[%s]", .cgpp(p, e$args[[1L]]$args[[1L]]), paste(vapply(e$args[-1L], function(i) .cgpp(p, i), ""), collapse = ", ")),
           `_print_Exp1` = function(p, e) "exp(1)",
           `_print_Pi` = function(p, e) "pi",
           `_print_Infinity` = function(p, e) "Inf",
           `_print_NegativeInfinity` = function(p, e) "-Inf",
           `_print_Assignment` = function(p, e) {
             lhs <- e$args[[1L]]; rhs <- e$args[[2L]]
             if (lhs$type == "MatrixSymbol") {
               lines <- character()
               for (ij in .cgp_call(p, "_traverse_matrix_indices", lhs)) {
                 temp <- Assignment(MatrixElement(lhs, ij[1L], ij[2L]), .cgp_mget(rhs, ij[1L], ij[2L]))
                 lines <- c(lines, .cgpp(p, temp))
               }
               return(paste(lines, collapse = "\n"))
             }
             if (isTRUE(.cgp_s(p, "contract")) && (.cg_has_type(lhs, "IndexedBase") || .cg_has_type(rhs, "IndexedBase")))
               return(.cgp_doprint_loops(p, rhs, lhs))
             .cgp_call(p, "_get_statement", sprintf("%s = %s", .cgpp(p, lhs), .cgpp(p, rhs)))
           },
           `_print_Piecewise` = function(p, e) {
             n <- length(e$args)
             last <- e$args[[n]]
             code <- if (.is_true(last$args[[2L]])) .cgpp(p, last$args[[1L]])
                     else sprintf("ifelse(%s,%s,NA)", .cgpp(p, last$args[[2L]]), .cgpp(p, last$args[[1L]]))
             for (pr in rev(e$args[-n]))
               code <- paste0(sprintf("ifelse(%s,%s,", .cgpp(p, pr$args[[2L]]), .cgpp(p, pr$args[[1L]])), code, ")")
             code
           },
           `_print_ITE` = function(p, e) .cgpp(p, .cg_ite_piecewise(e)),
           `_print_MatrixElement` = function(p, e) {
             nc <- .cgp_mshape(e$args[[1L]])[2L]
             sprintf("%s[%s]", .cgp_paren(p, e$args[[1L]], .CG_PRECEDENCE[["Atom"]], strict = TRUE),
                     sstr(.add(list(e$args[[3L]], .mul(list(e$args[[2L]], Integer(nc)))))))
           },
           `_print_Symbol` = function(p, e) {
             name <- .cgp_super(p, "RCodePrinter", "_print_Symbol", e)
             if (any(vapply(.cgp_s(p, "dereference"), function(d) is_sym(d) && identical(d$key, e$key), TRUE))) sprintf("(*%s)", name)
             else name
           },
           `_print_Relational` = function(p, e) sprintf("%s %s %s", .cgpp(p, e$args[[1L]]), .REL_OPS[[e$type]], .cgpp(p, e$args[[2L]])),
           `_print_AugmentedAssignment` = function(p, e)
             sprintf("%s %s %s;", .cgpp(p, e$args[[1L]]), .cg_aug_op(e), .cgpp(p, e$args[[2L]])),
           `_print_For` = function(p, e) {
             v <- .cg_vals(e)
             target <- .cgpp(p, v$target)
             if (!(is_sym(v$iterable) && v$iterable$type == "Range"))
               .cg_err("NotImplementedError", "Only iterable currently supported is Range")
             r <- v$iterable$args
             body <- .cgpp(p, v$body)
             sprintf("for(%s in seq(from=%s, to=%s, by=%s)){\n%s\n}", target, sstr(r[[1L]]), sstr(.sub(r[[2L]], S_One)),
                     sstr(r[[3L]]), body)
           },
           indent_code = function(p, code) .cgp_indent_braces(code, "   ")
         ))

#' @rdname code-printer-objects
#' @export
RCodePrinter <- function(settings = NULL) .cgp_new("RCodePrinter", settings)

.cg_rcode <- function(e, assign_to = NULL, ...) .cgp_doprint(.cgp_new("RCodePrinter", list(...)), e, assign_to)
