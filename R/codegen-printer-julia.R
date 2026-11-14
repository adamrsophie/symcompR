# Julia and Octave code printers, after sympy/printing/julia.py and
# octave.py.

.cg_kf_julia <- local({
  s1 <- c("sin", "cos", "tan", "cot", "sec", "csc", "asin", "acos", "atan", "acot", "asec", "acsc", "sinh", "cosh",
          "tanh", "coth", "sech", "csch", "asinh", "acosh", "atanh", "acoth", "asech", "acsch", "atan2", "sign",
          "floor", "log", "exp", "cbrt", "sqrt", "erf", "erfc", "erfi", "factorial", "gamma", "digamma", "trigamma",
          "polygamma", "beta", "airyai", "airyaiprime", "airybi", "airybiprime", "besselj", "bessely", "besseli",
          "besselk", "erfinv", "erfcinv")
  c(as.list(stats::setNames(s1, s1)), list(Abs = "abs", ceiling = "ceil", conjugate = "conj", hankel1 = "hankelh1",
                                    hankel2 = "hankelh2", im = "imag", re = "real"))
})
.cg_kf_octave <- local({
  s1 <- c("sin", "cos", "tan", "cot", "sec", "csc", "asin", "acos", "acot", "atan", "atan2", "asec", "acsc", "sinh",
          "cosh", "tanh", "coth", "csch", "sech", "asinh", "acosh", "atanh", "acoth", "asech", "acsch", "erfc", "erfi",
          "erf", "erfinv", "erfcinv", "besseli", "besselj", "besselk", "bessely", "bernoulli", "beta", "euler", "exp",
          "factorial", "floor", "fresnelc", "fresnels", "gamma", "harmonic", "log", "polylog", "sign", "zeta", "legendre")
  c(as.list(stats::setNames(s1, s1)), list(Abs = "abs", arg = "angle", binomial = "bincoeff", ceiling = "ceil",
                                    chebyshevu = "chebyshevU", chebyshevt = "chebyshevT", Chi = "coshint", Ci = "cosint",
                                    conjugate = "conj", DiracDelta = "dirac", Heaviside = "heaviside", im = "imag",
                                    laguerre = "laguerreL", LambertW = "lambertw", li = "logint", loggamma = "gammaln",
                                    Max = "max", Min = "min", Mod = "mod", polygamma = "psi", re = "real",
                                    RisingFactorial = "pochhammer", Shi = "sinhint", Si = "sinint"))
})

.cgp_is_number <- function(e) is_sym(e) && .expr_is_number(e)

# _print_Mul of the Julia and Octave printers; `octave` selects the Octave
# spacing and rational handling
.cgp_jo_mul <- function(p, e, octave) {
  if (.cgp_is_number(e) && isTRUE(ask_fact(e, "imaginary"))) {
    q <- .mul(list(.neg(S_I), e))
    ok <- if (octave) .mul(list(S_I, e))$type == "Integer" else isTRUE(ask_fact(.cgp_coeff_mul(e)[[1L]], "integer"))
    if (ok) return(paste0(.cgpp(p, q), if (octave) "i" else "im"))
  }
  prec <- .cgp_prec(e)
  cm <- .cgp_coeff_mul(e)
  if (.cgp_num_neg(cm[[1L]])) { e <- .cgp_keep_coeff(.neg(cm[[1L]]), cm[[2L]]); sign <- "-" } else sign <- ""
  ord <- .cgp_s(p, "order")
  args <- if (!(identical(ord, "old") || identical(ord, "none"))) { if (e$type == "Mul") .mul_print_order(e$args) else list(e) }
          else if (e$type == "Mul") e$args else list(e)
  a <- list(); b <- list(); pow_paren <- list()
  for (item in args) {
    ex <- if (item$type == "Pow") item$args[[2L]] else NULL
    if (!is.null(ex) && .is_comm(item) && ex$type %in% c("Integer", "Rational") && ex$value < 0) {
      if (!(ex$type == "Integer" && ex$value == -1)) b[[length(b) + 1L]] <- new_sym("Pow", list(item$args[[1L]], .neg(ex)))
      else {
        if (length(item$args[[1L]]$args) != 1L && item$args[[1L]]$type == "Mul") pow_paren[[length(pow_paren) + 1L]] <- item
        b[[length(b) + 1L]] <- .pow(item$args[[1L]], .neg(ex))
      }
    } else if (item$type %in% c("Integer", "Rational") && (octave || gmp::numerator(item$value) == 1)) {
      if (octave) {
        if (gmp::numerator(item$value) != 1) a[[length(a) + 1L]] <- Integer(gmp::numerator(item$value))
        if (gmp::denominator(item$value) != 1) b[[length(b) + 1L]] <- Integer(gmp::denominator(item$value))
      } else b[[length(b) + 1L]] <- Integer(gmp::denominator(item$value))
    } else a[[length(a) + 1L]] <- item
  }
  if (!length(a)) a <- list(S_One)
  a_str <- vapply(a, function(z) .cgp_paren(p, z, prec), "")
  b_str <- vapply(b, function(z) .cgp_paren(p, z, prec), "")
  bk <- vapply(b, .k, "")
  for (item in pow_paren) {
    hit <- match(item$args[[1L]]$key, bk)
    if (!is.na(hit)) b_str[hit] <- sprintf("(%s)", b_str[hit])
  }
  multjoin <- function(a, s) {
    r <- s[1L]
    if (length(a) > 1L) for (i in 2:length(a)) {
      ms <- if (.cgp_is_number(a[[i - 1L]])) "*" else ".*"
      r <- if (octave) paste0(r, ms, s[i]) else sprintf("%s %s %s", r, ms, s[i])
    }
    r
  }
  sp <- function(x) if (octave) x else paste0(" ", x, " ")
  if (!length(b)) return(paste0(sign, multjoin(a, a_str)))
  if (length(b) == 1L) {
    ds <- if (.cgp_is_number(b[[1L]])) "/" else "./"
    return(paste0(sign, multjoin(a, a_str), sp(ds), b_str[1L]))
  }
  ds <- if (all(vapply(b, .cgp_is_number, TRUE))) "/" else "./"
  paste0(sign, multjoin(a, a_str), sp(ds), "(", multjoin(b, b_str), ")")
}

.cgp_jo_pow <- function(p, e, octave) {
  b <- e$args[[1L]]; x <- e$args[[2L]]
  powsymbol <- if (all(vapply(e$args, .cgp_is_number, TRUE))) "^" else ".^"
  PREC <- .cgp_prec(e)
  if (.cgp_equal_valued(x, 0.5)) return(sprintf("sqrt(%s)", .cgpp(p, b)))
  if (.is_comm(e)) {
    sym <- if (.cgp_is_number(b)) "/" else "./"
    sp <- function(s) if (octave) s else paste0(" ", s, " ")
    if (.cgp_equal_valued(x, -0.5)) return(sprintf("1%ssqrt(%s)", sp(sym), .cgpp(p, b)))
    if (.cgp_equal_valued(x, -1)) return(sprintf("1%s%s", sp(sym), .cgp_paren(p, b, PREC)))
  }
  sprintf(if (octave) "%s%s%s" else "%s %s %s", .cgp_paren(p, b, PREC), powsymbol, .cgp_paren(p, x, PREC))
}

.cgp_jo_assignment <- function(p, e) {
  lhs <- e$args[[1L]]; rhs <- e$args[[2L]]
  if (!isTRUE(.cgp_s(p, "inline")) && rhs$type == "Piecewise") {
    prs <- lapply(rhs$args, function(pr) list(Assignment(lhs, pr$args[[1L]]), pr$args[[2L]]))
    return(.cgpp(p, new_sym("Piecewise", lapply(prs, function(q) new_sym("ExprCondPair", q)))))
  }
  if (isTRUE(.cgp_s(p, "contract")) && (.cg_has_type(lhs, "IndexedBase") || .cg_has_type(rhs, "IndexedBase")))
    return(.cgp_doprint_loops(p, rhs, lhs))
  .cgp_call(p, "_get_statement", sprintf("%s = %s", .cgpp(p, lhs), .cgpp(p, rhs)))
}

.cgp_jo_piecewise <- function(p, e, octave) {
  n <- length(e$args)
  if (!.is_true(e$args[[n]]$args[[2L]]))
    .cg_err("ValueError", paste0("All Piecewise expressions must contain an (expr, True) statement to be used as a default condition. ",
                                 "Without one, the generated expression may not evaluate to anything under some condition."))
  if (isTRUE(.cgp_s(p, "inline"))) {
    init <- e$args[-n]
    last <- .cgpp(p, e$args[[n]]$args[[1L]])
    if (octave) {
      pairs <- vapply(init, function(pr) sprintf("(%s).*(%s) + (~(%s)).*(", .cgpp(p, pr$args[[2L]]), .cgpp(p, pr$args[[1L]]), .cgpp(p, pr$args[[2L]])), "")
      return(paste0("(", paste(pairs, collapse = " ...\n"), last, strrep(")", length(pairs)), ")"))
    }
    pairs <- vapply(init, function(pr) sprintf("(%s) ? (%s) :", .cgpp(p, pr$args[[2L]]), .cgpp(p, pr$args[[1L]])), "")
    return(paste0("(", paste(pairs, collapse = "\n"), " (", last, "))"))
  }
  lines <- character()
  for (i in seq_len(n)) {
    c1 <- e$args[[i]]$args[[2L]]
    lines <- c(lines, if (i == 1L) sprintf("if (%s)", .cgpp(p, c1)) else if (i == n && .is_true(c1)) "else"
                      else sprintf("elseif (%s)", .cgpp(p, c1)))
    lines <- c(lines, .cgpp(p, e$args[[i]]$args[[1L]]))
    if (i == n) lines <- c(lines, "end")
  }
  paste(lines, collapse = "\n")
}

.cgp_jo_indent <- function(code, tab) {
  if (length(code) == 1L && isTRUE(attr(code, "str"))) {
    lines <- regmatches(code, gregexpr("[^\n]*\n|[^\n]+$", code))[[1L]]
    return(paste(.cgp_jo_indent(lines, tab), collapse = ""))
  }
  code <- sub("^[ \t]+", "", code)
  inc <- c("^function ", "^if ", "^elseif ", "^else$", "^for ")
  dec <- c("^end$", "^elseif ", "^else$")
  hit <- function(l, res) any(vapply(res, function(r) grepl(r, l, perl = TRUE), TRUE))
  level <- 0L; out <- character()
  for (line in code) {
    if (line %in% c("", "\n")) { out <- c(out, line); next }
    level <- level - hit(line, dec)
    out <- c(out, paste0(strrep(tab, level), line))
    level <- level + hit(line, inc)
  }
  out
}

.cgp_jo_loops <- function(p, indices) {
  open <- vapply(indices, function(i) sprintf("for %s = %s:%s", .cgpp(p, .cg_idx_label(i)), .cgpp(p, .add(list(.cgp_idx_lower(i), S_One))),
                                              .cgpp(p, .add(list(.cgp_idx_upper(i), S_One)))), "")
  list(open, rep("end", length(indices)))
}
.cgp_colmajor <- function(p, mat) {
  sh <- .cgp_mshape(mat)
  out <- list()
  for (j in seq_len(sh[2L]) - 1L) for (i in seq_len(sh[1L]) - 1L) out[[length(out) + 1L]] <- c(i, j)
  out
}
.cgp_jo_known <- function(p, kf, settings) {
  for (k in names(settings$user_functions)) kf[k] <- list(settings$user_functions[[k]])
  p$known_functions <- kf
}
.cgp_jo_default <- list(order = NULL, full_prec = "auto", error_on_reserved = FALSE, reserved_word_suffix = "_",
                        human = TRUE, inline = TRUE, allow_unknown_functions = FALSE, strict = NULL, precision = 17,
                        user_functions = list(), contract = TRUE)
.cgp_numsym_or <- function(p, e, s) if (isTRUE(.cgp_s(p, "inline"))) s else .cgp_print_numbersymbol(p, e)

.cgp_jo_slice <- function(p, e, octave) {
  sh <- .cgp_mshape(e$args[[1L]])
  strslice <- function(x, lim) {
    l <- .add(list(x$args[[1L]], S_One)); h <- x$args[[2L]]; step <- x$args[[3L]]
    lstr <- .cgpp(p, l)
    hstr <- if (.cgp_equal_valued(h, lim)) "end" else .cgpp(p, h)
    if (.cgp_equal_valued(step, 1)) {
      if (.cgp_equal_valued(l, 1) && .cgp_equal_valued(h, lim)) return(":")
      if (identical(l$key, h$key)) return(lstr)
      return(paste0(lstr, ":", hstr))
    }
    paste(lstr, .cgpp(p, step), hstr, sep = ":")
  }
  r <- strslice(e$args[[2L]], sh[1L]); cc <- strslice(e$args[[3L]], sh[2L])
  if (octave) sprintf("%s(%s, %s)", .cgpp(p, e$args[[1L]]), r, cc) else sprintf("%s[%s,%s]", .cgpp(p, e$args[[1L]]), r, cc)
}
# jn and yn through besselj and bessely
.cgp_jo_sph <- function(p, e, f) {
  nu <- e$args[[1L]]; x <- e$args[[2L]]
  .cgpp(p, .mul(list(sqrt(.div(S_Pi, .mul(list(S_Two, x)))), fcall(f, list(.add(list(nu, S_Half)), x)))))
}
.cgp_jo_hadamard_power <- function(p, e) {
  PREC <- .cgp_prec(e)
  paste0(.cgp_paren(p, e$args[[1L]], PREC), ".**", .cgp_paren(p, e$args[[2L]], PREC))
}

.cgp_def("JuliaCodePrinter", "CodePrinter",
         attrs = list(printmethod = "_julia", language = "Julia", module = "sympy.printing.julia",
                      `_operators` = list(and = "&&", or = "||", not = "!"), `_default_settings` = .cgp_jo_default),
         init = function(p, settings, nxt) { nxt(settings); .cgp_jo_known(p, .cg_kf_julia, settings) },
         methods = list(
           `_rate_index_position` = function(p, pos) pos * 5L,
           `_get_statement` = function(p, s) s,
           `_get_comment` = function(p, text) paste0("# ", text),
           `_declare_number_const` = function(p, name, value) sprintf("const %s = %s", .py_str_of(name), .py_str_of(value)),
           `_format_code` = function(p, lines) .cgp_call(p, "indent_code", lines),
           `_traverse_matrix_indices` = .cgp_colmajor,
           `_get_loop_opening_ending` = .cgp_jo_loops,
           `_print_Mul` = function(p, e) .cgp_jo_mul(p, e, FALSE),
           `_print_Relational` = function(p, e) sprintf("%s %s %s", .cgpp(p, e$args[[1L]]), .REL_OPS[[e$type]], .cgpp(p, e$args[[2L]])),
           `_print_Pow` = function(p, e) .cgp_jo_pow(p, e, FALSE),
           `_print_MatPow` = function(p, e) { PREC <- .cgp_prec(e); sprintf("%s ^ %s", .cgp_paren(p, e$args[[1L]], PREC), .cgp_paren(p, e$args[[2L]], PREC)) },
           `_print_Pi` = function(p, e) .cgp_numsym_or(p, e, "pi"),
           `_print_ImaginaryUnit` = function(p, e) "im",
           `_print_Exp1` = function(p, e) .cgp_numsym_or(p, e, "e"),
           `_print_EulerGamma` = function(p, e) .cgp_numsym_or(p, e, "eulergamma"),
           `_print_Catalan` = function(p, e) .cgp_numsym_or(p, e, "catalan"),
           `_print_GoldenRatio` = function(p, e) .cgp_numsym_or(p, e, "golden"),
           `_print_Assignment` = .cgp_jo_assignment,
           `_print_Infinity` = function(p, e) "Inf",
           `_print_NegativeInfinity` = function(p, e) "-Inf",
           `_print_NaN` = function(p, e) "NaN",
           `_print_list` = function(p, e) paste0("Any[", paste(vapply(e, function(a) .cgpp(p, a), ""), collapse = ", "), "]"),
           `_print_tuple` = function(p, e) {
             a <- if (is_sym(e)) e$args else e
             if (length(a) == 1L) sprintf("(%s,)", .cgpp(p, a[[1L]])) else sprintf("(%s)", .cgp_stringify(p, a, ", "))
           },
           `_print_Tuple` = function(p, e) .cgp_call(p, "_print_tuple", e),
           `_print_BooleanTrue` = function(p, e) "true",
           `_print_BooleanFalse` = function(p, e) "false",
           `_print_bool` = function(p, e) tolower(as.character(e)),
           `_print_MatrixBase` = function(p, A) {
             A <- .unwrap_explicit(A)
             if (A$nrow == 0L || A$ncol == 0L) return(sprintf("zeros(%d, %d)", A$nrow, A$ncol))
             if (A$nrow == 1L && A$ncol == 1L) return(sprintf("[%s]", sstr(A[1L, 1L])))
             row <- function(i) paste(vapply(seq_len(A$ncol), function(j) .cgpp(p, A[i, j]), ""), collapse = " ")
             if (A$nrow == 1L) return(sprintf("[%s]", row(1L)))
             if (A$ncol == 1L) return(sprintf("[%s]", paste(vapply(A$data, function(a) .cgpp(p, a), ""), collapse = ", ")))
             sprintf("[%s]", paste(.cgp_table(p, A), collapse = ";\n"))
           },
           `_print_MatrixElement` = function(p, e)
             paste0(.cgp_paren(p, e$args[[1L]], .CG_PRECEDENCE[["Atom"]], strict = TRUE),
                    sprintf("[%s,%s]", sstr(.add(list(e$args[[2L]], S_One))), sstr(.add(list(e$args[[3L]], S_One))))),
           `_print_Indexed` = function(p, e)
             sprintf("%s[%s]", .cgpp(p, e$args[[1L]]$args[[1L]]), paste(vapply(e$args[-1L], function(i) .cgpp(p, i), ""), collapse = ",")),
           `_print_Identity` = function(p, e) sprintf("eye(%s)", .cgpp(p, .mshape(e)[[1L]])),
           `_print_HadamardProduct` = function(p, e) paste(vapply(e$args, function(a) .cgp_paren(p, a, .cgp_prec(e)), ""), collapse = " .* "),
           `_print_Rational` = function(p, e) {
             if (gmp::denominator(e$value) == 1) return(as.character(gmp::numerator(e$value)))
             sprintf("%s // %s", as.character(gmp::numerator(e$value)), as.character(gmp::denominator(e$value)))
           },
           `_print_sinc` = function(p, e) sprintf("sinc(%s)", .cgpp(p, .div(e$args[[1L]], S_Pi))),
           `_print_Piecewise` = function(p, e) .cgp_jo_piecewise(p, e, FALSE),
           `_print_MatrixSlice` = function(p, e) .cgp_jo_slice(p, e, FALSE),
           `_print_HadamardPower` = function(p, e) .cgp_jo_hadamard_power(p, e),
           `_print_jn` = function(p, e) .cgp_jo_sph(p, e, "besselj"),
           `_print_yn` = function(p, e) .cgp_jo_sph(p, e, "bessely"),
           `_print_MatMul` = function(p, e) {
             args <- e$args
             sign <- ""
             if (is_number_type(args[[1L]]) && .cgp_num_neg(args[[1L]])) {
               c1 <- .neg(args[[1L]])
               args <- if (.cgp_equal_valued(c1, 1)) args[-1L] else c(list(c1), args[-1L])
               sign <- "-"
             }
             PREC <- .cgp_prec(e)
             paste0(sign, paste(vapply(args, function(a) .cgp_paren(p, a, PREC), ""), collapse = " * "))
           },
           indent_code = function(p, code) .cgp_jo_indent(code, "    ")
         ))

# A.table(self, rowstart='', rowend='', rowsep=..., colsep=' '): the
# columns are right aligned to their widest entry
.cgp_table <- function(p, A, colsep = " ") {
  s <- matrix(vapply(seq_len(A$nrow * A$ncol), function(k) {
    i <- (k - 1L) %/% A$ncol + 1L; j <- (k - 1L) %% A$ncol + 1L
    .cgpp(p, A[i, j])
  }, ""), nrow = A$nrow, byrow = TRUE)
  w <- apply(nchar(s), 2L, max)
  vapply(seq_len(A$nrow), function(i) paste(mapply(function(v, n) formatC(v, width = n), s[i, ], w), collapse = colsep), "")
}

.cgp_def("OctaveCodePrinter", "CodePrinter",
         attrs = list(printmethod = "_octave", language = "Octave", module = "sympy.printing.octave",
                      `_operators` = list(and = "&", or = "|", not = "~"), `_default_settings` = .cgp_jo_default),
         init = function(p, settings, nxt) { nxt(settings); .cgp_jo_known(p, .cg_kf_octave, settings) },
         methods = list(
           `_rate_index_position` = function(p, pos) pos * 5L,
           `_get_statement` = function(p, s) paste0(s, ";"),
           `_get_comment` = function(p, text) paste0("% ", text),
           `_declare_number_const` = function(p, name, value) sprintf("%s = %s;", .py_str_of(name), .py_str_of(value)),
           `_format_code` = function(p, lines) .cgp_call(p, "indent_code", lines),
           `_traverse_matrix_indices` = .cgp_colmajor,
           `_get_loop_opening_ending` = .cgp_jo_loops,
           `_print_Mul` = function(p, e) .cgp_jo_mul(p, e, TRUE),
           `_print_Relational` = function(p, e) sprintf("%s %s %s", .cgpp(p, e$args[[1L]]), .REL_OPS[[e$type]], .cgpp(p, e$args[[2L]])),
           `_print_Pow` = function(p, e) .cgp_jo_pow(p, e, TRUE),
           `_print_MatPow` = function(p, e) { PREC <- .cgp_prec(e); sprintf("%s^%s", .cgp_paren(p, e$args[[1L]], PREC), .cgp_paren(p, e$args[[2L]], PREC)) },
           `_print_MatrixSolve` = function(p, e) {
             PREC <- .cgp_prec(e)
             sprintf("%s \\ %s", .cgp_paren(p, e$args[[1L]], PREC), .cgp_paren(p, e$args[[2L]], PREC))
           },
           `_print_Pi` = function(p, e) "pi",
           `_print_ImaginaryUnit` = function(p, e) "1i",
           `_print_Exp1` = function(p, e) "exp(1)",
           `_print_GoldenRatio` = function(p, e) "(1+sqrt(5))/2",
           `_print_Assignment` = .cgp_jo_assignment,
           `_print_Infinity` = function(p, e) "inf",
           `_print_NegativeInfinity` = function(p, e) "-inf",
           `_print_NaN` = function(p, e) "NaN",
           `_print_list` = function(p, e) paste0("{", paste(vapply(e, function(a) .cgpp(p, a), ""), collapse = ", "), "}"),
           `_print_tuple` = function(p, e) paste0("{", paste(vapply(if (is_sym(e)) e$args else e, function(a) .cgpp(p, a), ""), collapse = ", "), "}"),
           `_print_Tuple` = function(p, e) .cgp_call(p, "_print_tuple", e),
           `_print_BooleanTrue` = function(p, e) "true",
           `_print_BooleanFalse` = function(p, e) "false",
           `_print_bool` = function(p, e) tolower(as.character(e)),
           `_print_MatrixBase` = function(p, A) {
             A <- .unwrap_explicit(A)
             if (A$nrow == 0L && A$ncol == 0L) return("[]")
             if (A$nrow == 0L || A$ncol == 0L) return(sprintf("zeros(%d, %d)", A$nrow, A$ncol))
             if (A$nrow == 1L && A$ncol == 1L) return(.cgpp(p, A[1L, 1L]))
             sprintf("[%s]", paste(vapply(seq_len(A$nrow), function(i)
               paste(vapply(seq_len(A$ncol), function(j) .cgpp(p, A[i, j]), ""), collapse = " "), ""), collapse = "; "))
           },
           `_print_MatrixElement` = function(p, e)
             paste0(.cgp_paren(p, e$args[[1L]], .CG_PRECEDENCE[["Atom"]], strict = TRUE),
                    sprintf("(%s, %s)", sstr(.add(list(e$args[[2L]], S_One))), sstr(.add(list(e$args[[3L]], S_One))))),
           `_print_Indexed` = function(p, e)
             sprintf("%s(%s)", .cgpp(p, e$args[[1L]]$args[[1L]]), paste(vapply(e$args[-1L], function(i) .cgpp(p, i), ""), collapse = ", ")),
           `_print_KroneckerDelta` = function(p, e) {
             prec <- .CG_PRECEDENCE[["Pow"]]
             sprintf("double(%s == %s)", .cgp_paren(p, e$args[[1L]], prec), .cgp_paren(p, e$args[[2L]], prec))
           },
           `_print_HadamardProduct` = function(p, e) paste(vapply(e$args, function(a) .cgp_paren(p, a, .cgp_prec(e)), ""), collapse = ".*"),
           `_print_Identity` = function(p, e) {
             sh <- .mshape(e)
             if (identical(sh[[1L]]$key, sh[[2L]]$key)) sh <- sh[1L]
             sprintf("eye(%s)", paste(vapply(sh, function(n) .cgpp(p, n), ""), collapse = ", "))
           },
           `_print_lowergamma` = function(p, e) sprintf("(gammainc(%s, %s).*gamma(%s))", .cgpp(p, e$args[[2L]]), .cgpp(p, e$args[[1L]]), .cgpp(p, e$args[[1L]])),
           `_print_uppergamma` = function(p, e) sprintf("(gammainc(%s, %s, 'upper').*gamma(%s))", .cgpp(p, e$args[[2L]]), .cgpp(p, e$args[[1L]]), .cgpp(p, e$args[[1L]])),
           `_print_sinc` = function(p, e) sprintf("sinc(%s)", .cgpp(p, .div(e$args[[1L]], S_Pi))),
           `_print_DiracDelta` = function(p, e) sprintf("%s(%s)", p$known_functions[[.cgp_clsname(e)]], paste(vapply(rev(e$args), function(a) .cgpp(p, a), ""), collapse = ", ")),
           `_print_LambertW` = function(p, e) sprintf("%s(%s)", p$known_functions[[.cgp_clsname(e)]], paste(vapply(rev(e$args), function(a) .cgpp(p, a), ""), collapse = ", ")),
           `_print_Max` = function(p, e) sprintf("%s(%s, %s)", p$known_functions[["Max"]], .cgpp(p, e$args[[1L]]), .cgpp(p, if (length(e$args) == 2L) e$args[[2L]] else do.call(Max, e$args[-1L]))),
           `_print_Min` = function(p, e) sprintf("%s(%s, %s)", p$known_functions[["Min"]], .cgpp(p, e$args[[1L]]), .cgpp(p, if (length(e$args) == 2L) e$args[[2L]] else do.call(Min, e$args[-1L]))),
           `_print_Piecewise` = function(p, e) .cgp_jo_piecewise(p, e, TRUE),
           `_print_MatrixSlice` = function(p, e) .cgp_jo_slice(p, e, TRUE),
           `_print_HadamardPower` = function(p, e) .cgp_jo_hadamard_power(p, e),
           `_print_jn` = function(p, e) .cgp_jo_sph(p, e, "besselj"),
           `_print_yn` = function(p, e) .cgp_jo_sph(p, e, "bessely"),
           `_print_airyai` = function(p, e) sprintf("airy(0, %s)", .cgpp(p, e$args[[1L]])),
           `_print_airyaiprime` = function(p, e) sprintf("airy(1, %s)", .cgpp(p, e$args[[1L]])),
           `_print_airybi` = function(p, e) sprintf("airy(2, %s)", .cgpp(p, e$args[[1L]])),
           `_print_airybiprime` = function(p, e) sprintf("airy(3, %s)", .cgpp(p, e$args[[1L]])),
           `_print_expint` = function(p, e) if (!.cgp_equal_valued(e$args[[1L]], 1)) .cgp_not_supported(p, e) else sprintf("expint(%s)", .cgpp(p, e$args[[2L]])),
           `_print_hankel1` = function(p, e) sprintf("besselh(%s, 1, %s)", .cgpp(p, e$args[[1L]]), .cgpp(p, e$args[[2L]])),
           `_print_hankel2` = function(p, e) sprintf("besselh(%s, 2, %s)", .cgpp(p, e$args[[1L]]), .cgpp(p, e$args[[2L]])),
           `_print_List` = function(p, e) .cgp_call(p, "_print_list", e$args),
           `_print_zeta` = function(p, e) if (length(e$args) == 1L) sprintf("zeta(%s)", .cgpp(p, e$args[[1L]])) else .cgp_not_supported(p, e),
           indent_code = function(p, code) .cgp_jo_indent(code, "  ")
         ))

#' @rdname code-printer-objects
#' @export
JuliaCodePrinter <- function(settings = NULL) .cgp_new("JuliaCodePrinter", settings)
#' @rdname code-printer-objects
#' @export
OctaveCodePrinter <- function(settings = NULL) .cgp_new("OctaveCodePrinter", settings)

.cg_julia_code <- function(e, assign_to = NULL, ...) .cgp_doprint(.cgp_new("JuliaCodePrinter", list(...)), e, assign_to)
.cg_octave_code <- function(e, assign_to = NULL, ...) .cgp_doprint(.cgp_new("OctaveCodePrinter", list(...)), e, assign_to)
