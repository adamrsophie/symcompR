# Fortran code printer, after sympy/printing/fortran.py.

.cg_kf_fortran <- list(sin = "sin", cos = "cos", tan = "tan", asin = "asin", acos = "acos", atan = "atan",
                       atan2 = "atan2", sinh = "sinh", cosh = "cosh", tanh = "tanh", log = "log", exp = "exp",
                       erf = "erf", Abs = "abs", conjugate = "conjg", Max = "max", Min = "min")

.cgp_f_lead <- function(p) {
  sf <- .cgp_s(p, "source_format")
  if (identical(sf, "fixed")) return(list(code = "      ", cont = "     @ ", comment = "C     "))
  if (identical(sf, "free")) return(list(code = "", cont = "      ", comment = "! "))
  .cg_err("ValueError", sprintf("Unknown source format: %s", .py_str_of(sf)))
}

# FCodePrinter._eval_function: the function of the numerically evaluated
# arguments
.cgp_f_eval_function <- function(p, e) {
  prec <- .cgp_s(p, "precision")
  args <- lapply(e$args, function(a) N(a, prec))
  if (e$type == "Piecewise") return(e)
  sym_rebuild(e, args)
}
.cgp_is_Function <- function(e) is_sym(e) && "Function" %in% .cgp_mro(e, full = TRUE)

.cgp_f_pad <- function(p, lines) {
  lead <- .cgp_f_lead(p)
  vapply(lines, function(line) {
    if (startsWith(line, "!")) paste0(lead$comment, sub("^\\s+", "", substring(line, 2L)))
    else paste0(lead$code, line)
  }, "", USE.NAMES = FALSE)
}

.cg_f_alnum <- strsplit(paste0("_+-.", "0123456789", "abcdefghijklmnopqrstuvwxyz", "ABCDEFGHIJKLMNOPQRSTUVWXYZ"), "")[[1L]]
.cg_f_white <- c(" ", "\t", "(", ")")

.cgp_f_wrap <- function(p, lines) {
  lead <- .cgp_f_lead(p)
  ch <- function(line, pos) substr(line, pos + 1L, pos + 1L)
  split_pos_code <- function(line, endpos) {
    if (nchar(line) <= endpos) return(nchar(line))
    pos <- endpos
    split <- function(pos) {
      a <- ch(line, pos); b <- ch(line, pos - 1L)
      (a %in% .cg_f_alnum && !(b %in% .cg_f_alnum)) || (!(a %in% .cg_f_alnum) && b %in% .cg_f_alnum) ||
        (a %in% .cg_f_white && !(b %in% .cg_f_white)) || (!(a %in% .cg_f_white) && b %in% .cg_f_white)
    }
    while (!split(pos)) {
      pos <- pos - 1L
      if (pos == 0L) return(endpos)
    }
    pos
  }
  # str.rfind(sub, start, end): last index in [start, end), or -1
  rfind <- function(s, start, end) {
    sub_s <- substr(s, start + 1L, end)
    hits <- gregexpr(" ", sub_s, fixed = TRUE)[[1L]]
    if (hits[1L] == -1L) -1L else start + max(hits) - 1L
  }
  trailing <- if (identical(.cgp_s(p, "source_format"), "free")) " &" else ""
  result <- character()
  for (line in lines) {
    if (startsWith(line, lead$comment)) {
      if (nchar(line) > 72L) {
        pos <- rfind(line, 6L, 72L)
        if (pos == -1L) pos <- 72L
        hunk <- substr(line, 1L, pos)
        line <- sub("^\\s+", "", substring(line, pos + 1L))
        result <- c(result, hunk)
        while (nzchar(line)) {
          pos <- rfind(line, 0L, 66L)
          if (pos == -1L || nchar(line) < 66L) pos <- 66L
          hunk <- substr(line, 1L, pos)
          line <- sub("^\\s+", "", substring(line, pos + 1L))
          result <- c(result, paste0(lead$comment, hunk))
        }
      } else result <- c(result, line)
    } else if (startsWith(line, lead$code)) {
      pos <- split_pos_code(line, 72L)
      hunk <- sub("\\s+$", "", substr(line, 1L, pos))
      line <- sub("^\\s+", "", substring(line, pos + 1L))
      if (nzchar(line)) hunk <- paste0(hunk, trailing)
      result <- c(result, hunk)
      while (nzchar(line)) {
        pos <- split_pos_code(line, 65L)
        hunk <- sub("\\s+$", "", substr(line, 1L, pos))
        line <- sub("^\\s+", "", substring(line, pos + 1L))
        if (nzchar(line)) hunk <- paste0(hunk, trailing)
        result <- c(result, paste0(lead$cont, hunk))
      }
    } else result <- c(result, line)
  }
  result
}

.cgp_f_indent <- function(p, code) {
  free <- identical(.cgp_s(p, "source_format"), "free")
  code <- sub("^[ \t]+", "", code)
  inc <- c("do ", "if(", "if ", "do\n", "else", "program", "interface")
  dec <- c("end do", "enddo", "end if", "endif", "else", "end program", "end interface")
  increase <- vapply(code, function(l) any(startsWith(l, inc)), TRUE)
  decrease <- vapply(code, function(l) any(startsWith(l, dec)), TRUE)
  continuation <- vapply(code, function(l) any(endsWith(l, c("&", "&\n"))), TRUE)
  level <- 0L; cont_padding <- 0L; tabwidth <- 3L
  new_code <- character()
  for (i in seq_along(code)) {
    line <- code[i]
    if (line %in% c("", "\n")) { new_code <- c(new_code, line); next }
    level <- level - decrease[i]
    padding <- if (free) strrep(" ", level * tabwidth + cont_padding) else strrep(" ", level * tabwidth)
    line <- paste0(padding, line)
    if (!free) line <- .cgp_f_pad(p, line)
    new_code <- c(new_code, line)
    cont_padding <- if (continuation[i]) 2L * tabwidth else 0L
    level <- level + increase[i]
  }
  new_code <- unname(new_code)
  if (!free) return(.cgp_f_wrap(p, new_code))
  new_code
}

.cgp_f_head <- function(p, entity, fp) {
  v <- .cg_vals(fp)
  bp <- attr_params(fp, "bind_C")
  bind <- if (is.null(bp)) "" else if (length(bp$args)) sprintf(" bind(C, name=\"%s\")", .py_str_of(bp$args[[1L]])) else " bind(C)"
  rn <- .cgp_s(p, "result_name")
  sprintf("%s%s(%s)%s%s\n%s", entity, .cgpp(p, v$name),
          paste(vapply(v$parameters$args, function(a) .cgpp(p, .cg_vals(a)$symbol), ""), collapse = ", "),
          if (!is.null(rn)) sprintf(" result(%s)", .py_str_of(rn)) else "", bind,
          paste(vapply(v$parameters$args, function(a) .cgpp(p, Declaration(a)), ""), collapse = "\n"))
}

.cgp_f_kw <- function(p, e, exclude = character()) lapply(ast_kwargs(e, exclude = exclude), function(a) .cgpp(p, a))

.cgp_def("FCodePrinter", "CodePrinter",
         attrs = list(printmethod = "_fcode", language = "Fortran", module = "sympy.printing.fortran",
                      `_default_settings` = list(order = NULL, full_prec = "auto", error_on_reserved = FALSE,
                                                 reserved_word_suffix = "_", human = TRUE, inline = FALSE,
                                                 allow_unknown_functions = FALSE, strict = NULL, precision = 17,
                                                 user_functions = list(), source_format = "fixed", contract = TRUE,
                                                 standard = 77, name_mangling = TRUE),
                      `_operators` = list(and = ".and.", or = ".or.", xor = ".neqv.", equivalent = ".eqv.", not = ".not. "),
                      `_relationals` = list(`!=` = "/=")),
         init = function(p, settings, nxt) {
           p$mangled_symbols <- list()
           p$used_name <- character()
           p$type_aliases <- .cgd_merge(.cgd(list(list(sym_integer, int32), list(.cg_real, float64), list(complex_, complex128))),
                                        settings$type_aliases)
           p$type_mappings <- .cgd_merge(.cgd(list(list(intc, "integer(c_int)"), list(float32, "real*4"), list(float64, "real*8"),
                                                   list(float80, "real*10"), list(complex64, "complex*8"),
                                                   list(complex128, "complex*16"), list(int8, "integer*1"),
                                                   list(int16, "integer*2"), list(int32, "integer*4"),
                                                   list(int64, "integer*8"), list(bool_, "logical"))),
                                         settings$type_mappings)
           p$type_modules <- .cgd(list(list(intc, list(iso_c_binding = "c_int"))))
           settings$type_aliases <- NULL; settings$type_mappings <- NULL
           nxt(settings)
           kf <- .cg_kf_fortran
           uf <- settings$user_functions
           for (k in names(uf)) kf[k] <- list(uf[[k]])
           p$known_functions <- kf
           std <- .cgp_s(p, "standard")
           if (!(std %in% c(66, 77, 90, 95, 2003, 2008)))
             .cg_err("ValueError", sprintf("Unknown Fortran standard: %s", .py_str_of(std)))
           p$module_uses <- list()
           p$context <- list()
           invisible(NULL)
         },
         methods = list(
           parenthesize = function(p, item, level, strict = FALSE) {
             if (.cgp_is_Function(item)) item <- .cgp_f_eval_function(p, item)
             pr <- .cgp_prec(item)
             if (pr < level || (!strict && pr <= level)) paste0("(", .cgpp(p, item), ")") else .cgpp(p, item)
           },
           `_print_Symbol` = function(p, e) {
             if (isTRUE(.cgp_s(p, "name_mangling"))) {
               if (!(e$key %in% names(p$mangled_symbols))) {
                 name <- e$name
                 while (tolower(name) %in% p$used_name) name <- paste0(name, "_")
                 p$used_name <- c(p$used_name, tolower(name))
                 ms <- p$mangled_symbols
                 ms[[e$key]] <- if (name == e$name) e else Symbol(name)
                 p$mangled_symbols <- ms
               }
               e <- p$mangled_symbols[[e$key]]
             }
             .cgp_super(p, "FCodePrinter", "_print_Symbol", e)
           },
           `_rate_index_position` = function(p, pos) -pos * 5L,
           `_get_statement` = function(p, s) s,
           `_get_comment` = function(p, text) paste0("! ", text),
           `_declare_number_const` = function(p, name, value) sprintf("parameter (%s = %s)", .py_str_of(name), .cgpp(p, value)),
           `_print_NumberSymbol` = function(p, e) {
             .cgp_set_add(p, "number_symbols", list(e, N(e, .cgp_s(p, "precision"))))
             sstr(e)
           },
           `_format_code` = function(p, lines) .cgp_f_wrap(p, .cgp_call(p, "indent_code", lines)),
           `_traverse_matrix_indices` = function(p, mat) {
             sh <- .cgp_mshape(mat)
             out <- list()
             for (j in seq_len(sh[2L]) - 1L) for (i in seq_len(sh[1L]) - 1L) out[[length(out) + 1L]] <- c(i, j)
             out
           },
           `_get_loop_opening_ending` = function(p, indices) {
             open <- vapply(indices, function(i) sprintf("do %s = %s, %s", .cgpp(p, .cg_idx_label(i)),
                                                         .cgpp(p, .add(list(.cgp_idx_lower(i), S_One))),
                                                         .cgpp(p, .add(list(.cgp_idx_upper(i), S_One)))), "")
             list(open, rep("end do", length(indices)))
           },
           `_print_sign` = function(p, e) {
             arg <- e$args[[1L]]
             if (isTRUE(ask_fact(arg, "integer"))) new <- fortran_merge(0L, isign(1L, arg), Eq(arg, 0L))
             else if (isTRUE(ask_fact(arg, "complex")) || isTRUE(ask_fact(arg, "infinite")))
               new <- fortran_merge(cmplx(literal_dp(0), literal_dp(0)), .div(arg, Abs(arg)), Eq(Abs(arg), literal_dp(0)))
             else new <- fortran_merge(literal_dp(0), dsign(literal_dp(1), arg), Eq(arg, literal_dp(0)))
             .cgpp(p, new)
           },
           `_print_Piecewise` = function(p, e) {
             n <- length(e$args)
             if (!.is_true(e$args[[n]]$args[[2L]]))
               .cg_err("ValueError", paste0("All Piecewise expressions must contain an (expr, True) statement to be used as a default condition. ",
                                            "Without one, the generated expression may not evaluate to anything under some condition."))
             if (.cg_has_type(e, "Assignment")) {
               lines <- character()
               for (i in seq_len(n)) {
                 c1 <- e$args[[i]]$args[[2L]]
                 if (i == 1L) lines <- c(lines, sprintf("if (%s) then", .cgpp(p, c1)))
                 else if (i == n && .is_true(c1)) lines <- c(lines, "else")
                 else lines <- c(lines, sprintf("else if (%s) then", .cgpp(p, c1)))
                 lines <- c(lines, .cgpp(p, e$args[[i]]$args[[1L]]))
               }
               return(paste(c(lines, "end if"), collapse = "\n"))
             }
             if (.cgp_s(p, "standard") >= 95) {
               code <- .cgpp(p, e$args[[n]]$args[[1L]])
               if (n > 1L) for (i in rev(seq_len(n - 1L)))
                 code <- sprintf("merge(%s, %s, %s)", .cgpp(p, e$args[[i]]$args[[1L]]), code, .cgpp(p, e$args[[i]]$args[[2L]]))
               return(code)
             }
             .cg_err("NotImplementedError", "Using Piecewise as an expression using inline operators is not supported in standards earlier than Fortran95.")
           },
           `_print_MatrixElement` = function(p, e)
             sprintf("%s(%s, %s)", .cgp_paren(p, e$args[[1L]], .CG_PRECEDENCE[["Atom"]], strict = TRUE),
                     sstr(.add(list(e$args[[2L]], S_One))), sstr(.add(list(e$args[[3L]], S_One)))),
           `_print_Add` = function(p, e) {
             pr <- list(); pim <- list(); mixed <- list()
             for (a in e$args) {
               if (.expr_is_number(a) && isTRUE(ask_fact(a, "real"))) pr[[length(pr) + 1L]] <- a
               else if (.expr_is_number(a) && isTRUE(ask_fact(a, "imaginary"))) pim[[length(pim) + 1L]] <- a
               else mixed[[length(mixed) + 1L]] <- a
             }
             if (length(pim)) {
               im <- .cgpp(p, .mul(list(.neg(S_I), .add(pim))))
               re_s <- .cgpp(p, .add(pr))
               if (length(mixed)) {
                 PREC <- .cgp_prec(e)
                 term <- .add(mixed)
                 t <- .cgpp(p, term)
                 if (startsWith(t, "-")) { sign <- "-"; t <- substring(t, 2L) } else sign <- "+"
                 if (.cgp_prec(term) < PREC) t <- sprintf("(%s)", t)
                 return(sprintf("cmplx(%s,%s) %s %s", re_s, im, sign, t))
               }
               return(sprintf("cmplx(%s,%s)", re_s, im))
             }
             .cgp_super(p, "FCodePrinter", "_print_Add", e)
           },
           `_print_Function` = function(p, e) {
             ev <- .cgp_f_eval_function(p, e)
             if (!.cgp_is_Function(ev)) return(.cgpp(p, ev))
             .cgp_print_function(p, ev)
           },
           `_print_Mod` = function(p, e) {
             if (.cgp_s(p, "standard") %in% c(66, 77))
               .cg_err("NotImplementedError", "Python % operator and SymPy's Mod() function are not supported by Fortran 66 or 77 standards.")
             sprintf("      modulo(%s, %s)", .cgpp(p, e$args[[1L]]), .cgpp(p, e$args[[2L]]))
           },
           `_print_ImaginaryUnit` = function(p, e) "cmplx(0,1)",
           `_print_int` = function(p, e) as.character(e),
           `_print_Mul` = function(p, e) {
             if (.expr_is_number(e) && isTRUE(ask_fact(e, "imaginary")))
               return(sprintf("cmplx(0,%s)", .cgpp(p, .mul(list(.neg(S_I), e)))))
             .cgp_print_mul(p, e)
           },
           `_print_Pow` = function(p, e) {
             PREC <- .cgp_prec(e)
             b <- e$args[[1L]]; x <- e$args[[2L]]
             if (.cgp_equal_valued(x, -1)) return(sprintf("%s/%s", .cgpp(p, literal_dp(1)), .cgp_paren(p, b, PREC)))
             if (.cgp_equal_valued(x, 0.5)) {
               if (isTRUE(ask_fact(b, "integer"))) {
                 if (is_number_type(b)) return(sprintf("sqrt(%s.0d0)", .cgpp(p, b)))
                 return(sprintf("sqrt(dble(%s))", .cgpp(p, b)))
               }
               return(sprintf("sqrt(%s)", .cgpp(p, b)))
             }
             .cgp_super(p, "FCodePrinter", "_print_Pow", e)
           },
           `_print_Rational` = function(p, e)
             sprintf("%s.0d0/%s.0d0", as.character(gmp::numerator(e$value)), as.character(gmp::denominator(e$value))),
           `_print_Float` = function(p, e) {
             printed <- .cgp_super(p, "FCodePrinter", "_print_Float", e)
             k <- regexpr("e", printed, fixed = TRUE)
             if (k > 0L) return(sprintf("%sd%s", substr(printed, 1L, k - 1L), substring(printed, k + 1L)))
             paste0(printed, "d0")
           },
           `_print_Relational` = function(p, e) {
             op <- .REL_OPS[[e$type]]
             rel <- .cgp_a(p, "_relationals")
             if (op %in% names(rel)) op <- rel[[op]]
             sprintf("%s %s %s", .cgpp(p, e$args[[1L]]), op, .cgpp(p, e$args[[2L]]))
           },
           `_print_Indexed` = function(p, e)
             sprintf("%s(%s)", .cgpp(p, e$args[[1L]]$args[[1L]]), paste(vapply(e$args[-1L], function(i) .cgpp(p, i), ""), collapse = ", ")),
           `_print_AugmentedAssignment` = function(p, e) {
             lhs <- .cgpp(p, e$args[[1L]]); rhs <- .cgpp(p, e$args[[2L]])
             .cgp_call(p, "_get_statement", sprintf("%s = %s %s %s", .cgpp(p, lhs), .cgpp(p, lhs), .cgpp(p, .AUG_OPS[[e$type]]), .cgpp(p, rhs)))
           },
           `_print_sum_` = function(p, e) {
             v <- .cg_vals(e)
             params <- .cgpp(p, v$array)
             if (!.is_none(v$dim)) params <- paste0(params, ", ", .cgpp(p, v$dim))
             if (!.is_none(v$mask)) params <- paste0(params, ", mask=", .cgpp(p, v$mask))
             sprintf("%s(%s)", sub("_+$", "", e$type), params)
           },
           `_print_product_` = function(p, e) .cgp_call(p, "_print_sum_", e),
           `_print_Do` = function(p, e) {
             v <- .cg_vals(e)
             excl <- "concurrent"
             if (identical(v$step$key, S_One$key)) excl <- c(excl, "step")
             kw <- .cgp_f_kw(p, e, excl)
             sprintf("do %s%s = %s, %s%s\n%s\nend do\n", if (.is_true(v$concurrent)) "concurrent " else "", kw$counter,
                     kw$first, kw$last, if ("step" %in% excl) "" else paste0(", ", kw$step), kw$body)
           },
           `_print_ImpliedDoLoop` = function(p, e) {
             v <- .cg_vals(e)
             kw <- .cgp_f_kw(p, e)
             sprintf("(%s, %s = %s, %s%s)", kw$expr, kw$counter, kw$first, kw$last,
                     if (identical(v$step$key, S_One$key)) "" else paste0(", ", kw$step))
           },
           `_print_For` = function(p, e) {
             v <- .cg_vals(e)
             target <- .cgpp(p, v$target)
             if (!(is_sym(v$iterable) && v$iterable$type == "Range"))
               .cg_err("NotImplementedError", "Only iterable currently supported is Range")
             r <- v$iterable$args
             body <- .cgpp(p, v$body)
             sprintf("do %s = %s, %s, %s\n%s\nend do", target, sstr(r[[1L]]), sstr(.sub(r[[2L]], S_One)), sstr(r[[3L]]), body)
           },
           `_print_Type` = function(p, t) {
             t <- .cgd_get(p$type_aliases, t, t)
             ts <- .cgd_get(p$type_mappings, t)
             if (is.null(ts)) ts <- .cg_text(t$args[[1L]])
             mu <- .cgd_get(p$type_modules, t)
             if (length(mu)) .cg_err("ValueError", "too many values to unpack (expected 2)")
             ts
           },
           `_print_Element` = function(p, e) {
             v <- .cg_vals(e)
             sprintf("%s(%s)", .cgpp(p, v$symbol), paste(vapply(v$indices$args, function(a) .cgpp(p, a), ""), collapse = ", "))
           },
           `_print_Extent` = function(p, e) sstr(e),
           `_print_Declaration` = function(p, e) {
             var <- e$args[[1L]]
             vv <- .cg_vals(var)
             val <- vv$value
             dim <- attr_params(var, "dimension")
             intents <- vapply(list(intent_in, intent_out, intent_inout), function(a) .cg_in_attrs(var, a), TRUE)
             if (sum(intents) == 0L) intent <- ""
             else if (sum(intents) == 1L) intent <- sprintf(", intent(%s)", c("in", "out", "inout")[which(intents)])
             else .cg_err("ValueError", sprintf("Multiple intents specified for %s", "<printer>"))
             if (.cg_isa(var, "Pointer")) .cg_err("NotImplementedError", "Pointers are not available by default in Fortran.")
             if (.cgp_s(p, "standard") >= 90) {
               result <- sprintf("%s%s%s%s%s :: %s", .cgpp(p, vv$type), if (.cg_in_attrs(var, value_const)) ", parameter" else "",
                                 if (!is.null(dim) && length(dim$args)) sprintf(", dimension(%s)", paste(vapply(dim$args, function(a) .cgpp(p, a), ""), collapse = ", ")) else "",
                                 intent, if (.cg_in_attrs(var, allocatable)) ", allocatable" else "", .cgpp(p, vv$symbol))
               if (!.is_none(val)) result <- paste0(result, " = ", .cgpp(p, val))
             } else {
               if (.cg_in_attrs(var, value_const) || .cgp_truth(val))
                 .cg_err("NotImplementedError", "F77 init./parameter statem. req. multiple lines.")
               result <- paste(.cgpp(p, vv$type), .cgpp(p, vv$symbol))
             }
             result
           },
           `_print_Infinity` = function(p, e) sprintf("(huge(%s) + 1)", .cgpp(p, literal_dp(0))),
           `_print_While` = function(p, e) {
             kw <- .cgp_f_kw(p, e)
             sprintf("do while (%s)\n%s\nend do", kw$condition, kw$body)
           },
           `_print_BooleanTrue` = function(p, e) ".true.",
           `_print_BooleanFalse` = function(p, e) ".false.",
           indent_code = function(p, code) .cgp_f_indent(p, code),
           `_print_GoTo` = function(p, e) {
             v <- .cg_vals(e)
             if (!.is_none(v$expr) && .cgp_truth(v$expr))
               return(sprintf("go to (%s), %s", paste(vapply(v$labels$args, function(a) .cgpp(p, a), ""), collapse = ", "), .cgpp(p, v$expr)))
             sprintf("go to %s", .cgpp(p, v$labels$args[[1L]]))
           },
           `_print_Program` = function(p, e) {
             kw <- .cgp_f_kw(p, e)
             sprintf("program %s\n%s\nend program\n", kw$name, kw$body)
           },
           `_print_Module` = function(p, e) {
             kw <- .cgp_f_kw(p, e)
             sprintf("module %s\n%s\n\ncontains\n\n%s\nend module\n", kw$name, kw$declarations, kw$definitions)
           },
           `_print_Stream` = function(p, e) {
             nm <- .cg_text(.cg_vals(e)$name)
             std <- .cgp_s(p, "standard")
             if (nm %in% c("stdout", "stderr") && std >= 2003) {
               u <- if (nm == "stdout") "input_unit" else "error_unit"
               mu <- p$module_uses
               mu[["iso_c_binding"]] <- unique(c(mu[["iso_c_binding"]], paste0("stdint=>", u)))
               p$module_uses <- mu
               return(u)
             }
             if (nm == "stdout") "*" else nm
           },
           `_print_Print` = function(p, e) {
             v <- .cg_vals(e)
             iolist <- paste(vapply(v$print_args$args, function(a) .cgpp(p, a), ""), collapse = ", ")
             if (.is_none(v$format_string)) return(sprintf("print %s, %s", "*", iolist))
             out <- if (identical(v$file$key, stderr_$key)) "0" else if (identical(v$file$key, stdout_$key)) "6" else "*"
             sprintf("write(%s, fmt=\"%s\", advance=\"no\"), %s", out, .cgpp(p, v$format_string), iolist)
           },
           `_print_Return` = function(p, e) {
             rn <- p$context$result_name
             sprintf("%s = %s", if (is.null(rn)) "sympy_result" else .py_str_of(rn), .cgpp(p, e$args[[1L]]))
           },
           `_print_FortranReturn` = function(p, e) {
             arg <- e$args[[1L]]
             if (!.is_none(arg) && .cgp_truth(arg)) paste0("return ", .cgpp(p, arg)) else "return"
           },
           `_print_FunctionPrototype` = function(p, e) {
             entity <- sprintf("%s function ", .cgpp(p, .cg_vals(e)$return_type))
             sprintf("interface\n%s\nend function\nend interface", .cgp_f_head(p, entity, e))
           },
           `_print_FunctionDefinition` = function(p, e) {
             v <- .cg_vals(e)
             prefix <- if (.cg_in_attrs(e, elemental)) "elemental " else if (.cg_in_attrs(e, pure)) "pure " else ""
             entity <- sprintf("%s function ", .cgpp(p, v$return_type))
             old <- p$context
             p$context <- c(old, list(result_name = v$name))
             on.exit(p$context <- old)
             sprintf("%s%s\n%s\nend function\n", prefix, .cgp_f_head(p, entity, e), .cgpp(p, v$body))
           },
           `_print_Subroutine` = function(p, e)
             sprintf("%s\n%s\nend subroutine\n", .cgp_f_head(p, "subroutine ", e), .cgpp(p, .cg_vals(e)$body)),
           `_print_SubroutineCall` = function(p, e) {
             v <- .cg_vals(e)
             sprintf("call %s(%s)", .cgpp(p, v$name), paste(vapply(v$subroutine_args$args, function(a) .cgpp(p, a), ""), collapse = ", "))
           },
           `_print_use_rename` = function(p, e) sprintf("%s => %s", .cgpp(p, e$args[[1L]]), .cgpp(p, e$args[[2L]])),
           `_print_use` = function(p, e) {
             v <- .cg_vals(e)
             result <- paste0("use ", .cgpp(p, v$namespace))
             if (!.is_none(v$rename)) result <- paste0(result, ", ", paste(vapply(v$rename$args, function(a) .cgpp(p, a), ""), collapse = ", "))
             if (!.is_none(v$only)) result <- paste0(result, ", only: ", paste(vapply(v$only$args, function(a) .cgpp(p, a), ""), collapse = ", "))
             result
           },
           `_print_BreakToken` = function(p, e) "exit",
           `_print_ContinueToken` = function(p, e) "cycle",
           `_print_ArrayConstructor` = function(p, e) {
             fmt <- if (.cgp_s(p, "standard") >= 2003) "[%s]" else "(/%s/)"
             sprintf(fmt, paste(vapply(e$args[[1L]]$args, function(a) .cgpp(p, a), ""), collapse = ", "))
           },
           `_print_KeywordFunctionCall` = function(p, e) {
             v <- .cg_vals(e)
             args <- vapply(v$function_args$args, function(a) .cgpp(p, a), "")
             kw <- v$keyword_args
             for (pr in kw$args) args <- c(args, sprintf("%s=%s", .cg_text(pr$args[[1L]]), .cgpp(p, pr$args[[2L]])))
             sprintf("%s(%s)", .cgpp(p, v$name), paste(args, collapse = ", "))
           }
         ))

# FFunction._fcode and _literal._fcode
.cgp_nodemethods[["_fcode"]] <- list(
  FFunction = function(p, e) {
    name <- e$name
    req <- .cg_ffunctions[[name]]
    if (.cgp_s(p, "standard") < req) .cg_err("NotImplementedError", sprintf("%s requires Fortran %d or newer", name, as.integer(req)))
    sprintf("%s(%s)", name, paste(vapply(e$args, function(a) .cgpp(p, a), ""), collapse = ", "))
  },
  `_literal` = function(p, e) {
    dec <- if (e$extra$literal == "literal_dp") 17L else 9L
    tok <- if (e$extra$literal == "literal_dp") "d" else "e"
    s <- strsplit(sprintf("%.*e", dec, as.numeric(e$value)), "e", fixed = TRUE)[[1L]]
    mant <- sub("\\.$", "", gsub("^0+|0+$", "", s[1L]))
    ex_sgn <- substr(s[2L], 1L, 1L); ex_num <- sub("^0+", "", substring(s[2L], 2L))
    if (ex_sgn == "+") ex_sgn <- ""
    paste0(if (nzchar(mant)) mant else "0", tok, ex_sgn, if (nzchar(ex_num)) ex_num else "0")
  })

#' @rdname code-printer-objects
#' @export
FCodePrinter <- function(settings = NULL) .cgp_new("FCodePrinter", settings)

.cg_fcode <- function(e, assign_to = NULL, ...) .cgp_doprint(.cgp_new("FCodePrinter", list(...)), e, assign_to)
