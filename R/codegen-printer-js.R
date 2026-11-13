# JavaScript code printer, after sympy/printing/jscode.py.

.cg_kf_js <- list(Abs = "Math.abs", acos = "Math.acos", acosh = "Math.acosh", asin = "Math.asin", asinh = "Math.asinh",
                  atan = "Math.atan", atan2 = "Math.atan2", atanh = "Math.atanh", ceiling = "Math.ceil", cos = "Math.cos",
                  cosh = "Math.cosh", exp = "Math.exp", floor = "Math.floor", log = "Math.log", Max = "Math.max",
                  Min = "Math.min", sign = "Math.sign", sin = "Math.sin", sinh = "Math.sinh", tan = "Math.tan",
                  tanh = "Math.tanh")

.cgp_def("JavascriptCodePrinter", "CodePrinter",
         attrs = list(printmethod = "_javascript", language = "JavaScript", module = "sympy.printing.jscode",
                      `_default_settings` = list(order = NULL, full_prec = "auto", error_on_reserved = FALSE,
                                                 reserved_word_suffix = "_", human = TRUE, inline = FALSE,
                                                 allow_unknown_functions = FALSE, strict = NULL, precision = 17,
                                                 user_functions = list(), contract = TRUE)),
         init = function(p, settings, nxt) { nxt(settings); .cgp_jo_known(p, .cg_kf_js, settings) },
         methods = list(
           `_rate_index_position` = function(p, pos) pos * 5L,
           `_get_statement` = function(p, s) paste0(s, ";"),
           `_get_comment` = function(p, text) paste0("// ", text),
           `_declare_number_const` = function(p, name, value)
             sprintf("var %s = %s;", .py_str_of(name), sstr(N(value, .cgp_s(p, "precision")))),
           `_format_code` = function(p, lines) .cgp_call(p, "indent_code", lines),
           `_traverse_matrix_indices` = function(p, mat) {
             sh <- .cgp_mshape(mat)
             out <- list()
             for (i in seq_len(sh[1L]) - 1L) for (j in seq_len(sh[2L]) - 1L) out[[length(out) + 1L]] <- c(i, j)
             out
           },
           `_get_loop_opening_ending` = function(p, indices) {
             open <- vapply(indices, function(i) {
               v <- .cgpp(p, .cg_idx_label(i))
               sprintf("for (var %s=%s; %s<%s; %s++){", v, .cgpp(p, .cgp_idx_lower(i)), v,
                       .cgpp(p, .add(list(.cgp_idx_upper(i), S_One))), v)
             }, "")
             list(open, rep("}", length(indices)))
           },
           `_print_Pow` = function(p, e) {
             PREC <- .cgp_prec(e)
             b <- e$args[[1L]]; x <- e$args[[2L]]
             if (.cgp_equal_valued(x, -1)) return(paste0("1/", .cgp_paren(p, b, PREC)))
             if (.cgp_equal_valued(x, 0.5)) return(sprintf("Math.sqrt(%s)", .cgpp(p, b)))
             if (identical(x$key, Rational(1L, 3L)$key)) return(sprintf("Math.cbrt(%s)", .cgpp(p, b)))
             sprintf("Math.pow(%s, %s)", .cgpp(p, b), .cgpp(p, x))
           },
           `_print_Rational` = function(p, e) paste0(as.character(gmp::numerator(e$value)), "/", as.character(gmp::denominator(e$value))),
           `_print_Mod` = function(p, e) {
             num <- e$args[[1L]]; den <- e$args[[2L]]
             PREC <- .cgp_prec(e)
             den_prec <- if (den$type == "Pow" && isTRUE(ask_fact(den$args[[2L]], "negative"))) .CG_PRECEDENCE[["Pow"]] else PREC
             snum <- .cgp_paren(p, num, PREC); sden <- .cgp_paren(p, den, den_prec)
             f <- function(a, k) isTRUE(ask_fact(a, k))
             if ((f(num, "nonnegative") && f(den, "nonnegative")) || (f(num, "nonpositive") && f(den, "nonpositive")))
               return(sprintf("%s %% %s", snum, sden))
             sprintf("((%s %% %s) + %s) %% %s", snum, sden, sden, sden)
           },
           `_print_Relational` = function(p, e) sprintf("%s %s %s", .cgpp(p, e$args[[1L]]), .REL_OPS[[e$type]], .cgpp(p, e$args[[2L]])),
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
           `_print_Exp1` = function(p, e) "Math.E",
           `_print_Pi` = function(p, e) "Math.PI",
           `_print_Infinity` = function(p, e) "Number.POSITIVE_INFINITY",
           `_print_NegativeInfinity` = function(p, e) "Number.NEGATIVE_INFINITY",
           `_print_Piecewise` = function(p, e) {
             n <- length(e$args)
             if (!.is_true(e$args[[n]]$args[[2L]]))
               .cg_err("ValueError", paste0("All Piecewise expressions must contain an (expr, True) statement to be used as a default condition. ",
                                            "Without one, the generated expression may not evaluate to anything under some condition."))
             if (.cg_has_type(e, "Assignment")) {
               lines <- character()
               for (i in seq_len(n)) {
                 c1 <- e$args[[i]]$args[[2L]]
                 lines <- c(lines, if (i == 1L) sprintf("if (%s) {", .cgpp(p, c1)) else if (i == n && .is_true(c1)) "else {"
                                   else sprintf("else if (%s) {", .cgpp(p, c1)), .cgpp(p, e$args[[i]]$args[[1L]]), "}")
               }
               return(paste(lines, collapse = "\n"))
             }
             pairs <- vapply(e$args[-n], function(pr) sprintf("((%s) ? (\n%s\n)\n", .cgpp(p, pr$args[[2L]]), .cgpp(p, pr$args[[1L]])), "")
             paste0(paste(pairs, collapse = ": "), sprintf(": (\n%s\n)", .cgpp(p, e$args[[n]]$args[[1L]])), strrep(")", length(pairs)))
           },
           `_print_MatrixElement` = function(p, e) {
             nc <- .cgp_mshape(e$args[[1L]])[2L]
             sprintf("%s[%s]", .cgp_paren(p, e$args[[1L]], .CG_PRECEDENCE[["Atom"]], strict = TRUE),
                     sstr(.add(list(e$args[[3L]], .mul(list(e$args[[2L]], Integer(nc)))))))
           },
           indent_code = function(p, code) .cgp_indent_braces(code, "   ")
         ))

#' @rdname code-printer-objects
#' @export
JavascriptCodePrinter <- function(settings = NULL) .cgp_new("JavascriptCodePrinter", settings)

.cg_jscode <- function(e, assign_to = NULL, ...) .cgp_doprint(.cgp_new("JavascriptCodePrinter", list(...)), e, assign_to)
