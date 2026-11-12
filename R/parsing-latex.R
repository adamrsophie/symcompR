# LaTeX parsing, after sympy/parsing/latex/__init__.py and the lark backend in
# sympy/parsing/latex/lark (latex_parser.py and transformer.py). The parse
# tree comes from the Earley parser in parsing-latex-earley.R; the transformer
# below turns it into an expression rule by rule, as TransformToSymPyExpr does.

.latex_error <- function(message) {
  structure(class = c("LaTeXParsingError", "sympy_error", "error", "condition"),
            list(message = message, call = NULL))
}

.LATEX_IGNORE_L <- "\\s*[{]*\\s*"
.LATEX_IGNORE_R <- "\\s*[}]*\\s*"
.LATEX_NO_LEFT <- "(?<!\\\\left)"

# check_matrix_delimiters: MATRIX_DELIMS as (begin regex, end regex, begin
# representative, end representative), in SymPy's dict order.
.latex_matrix_delims <- function() {
  L <- .LATEX_IGNORE_L; R <- .LATEX_IGNORE_R; N <- .LATEX_NO_LEFT
  bm <- "\\\\begin{matrix}"; em <- "\\\\end{matrix}"
  ba <- "\\\\begin{array}{.*?}"; ea <- "\\\\end{array}"
  d <- list(
    c(paste0("\\\\left\\(", L, bm), paste0(em, R, "\\\\right\\)"), "\\left(\\begin{matrix}", "\\end{matrix}\\right)"),
    c(paste0(N, "\\(", L, bm), paste0(em, R, "\\)"), "(\\begin{matrix}", "\\end{matrix})"),
    c(paste0("\\\\left\\[", L, bm), paste0(em, R, "\\\\right\\]"), "\\left[\\begin{matrix}", "\\end{matrix}\\right]"),
    c(paste0(N, "\\[", L, bm), paste0(em, R, "\\]"), "[\\begin{matrix}", "\\end{matrix}]"),
    c(paste0("\\\\left\\|", L, bm), paste0(em, R, "\\\\right\\|"), "\\left|\\begin{matrix}", "\\end{matrix}\\right|"),
    c(paste0(N, "\\|", L, bm), paste0(em, R, "\\|"), "|\\begin{matrix}", "\\end{matrix}|"),
    c("\\\\begin{pmatrix}", "\\\\end{pmatrix}", "\\begin{pmatrix}", "\\end{pmatrix}"),
    c("\\\\begin{bmatrix}", "\\\\end{bmatrix}", "\\begin{bmatrix}", "\\end{bmatrix}"),
    c("\\\\begin{vmatrix}", "\\\\end{vmatrix}", "\\begin{vmatrix}", "\\end{vmatrix}"),
    c(paste0("\\\\left\\(", L, ba), paste0(ea, R, "\\\\right\\)"), "\\left(\\begin{array}{COLUMN_SPECIFIERS}", "\\end{array}\\right)"),
    c(paste0(N, "\\(", L, ba), paste0(ea, R, "\\)"), "(\\begin{array}{COLUMN_SPECIFIERS}", "\\end{array})"),
    c(paste0("\\\\left\\[", L, ba), paste0(ea, R, "\\\\right\\]"), "\\left[\\begin{array}{COLUMN_SPECIFIERS}", "\\end{array}\\right]"),
    c(paste0(N, "\\[", L, ba), paste0(ea, R, "\\]"), "[\\begin{array}{COLUMN_SPECIFIERS}", "\\end{array}]"),
    c(paste0("\\\\left\\|", L, ba), paste0(ea, R, "\\\\right\\|"), "\\left|\\begin{array}{COLUMN_SPECIFIERS}", "\\end{array}\\right|"),
    c(paste0(N, "\\|", L, ba), paste0(ea, R, "\\|"), "|\\begin{array}{COLUMN_SPECIFIERS}", "\\end{array}|"))
  d
}

.latex_check_matrix_delimiters <- function(s) {
  d <- .latex_matrix_delims()
  n <- nchar(s)
  spans <- list()
  for (k in seq_along(d)) {
    for (w in 1:2) {
      m <- gregexpr(d[[k]][w], s, perl = TRUE)[[1L]]
      if (m[1L] == -1L) next
      for (j in seq_along(m))
        spans[[length(spans) + 1L]] <- list(start = m[j] - 1L, end = m[j] - 1L + attr(m, "match.length")[j],
                                           text = substr(s, m[j], m[j] + attr(m, "match.length")[j] - 1L),
                                           delim = k, is_end = w == 2L)
    }
  }
  if (!length(spans)) return(invisible())
  spans <- spans[order(vapply(spans, `[[`, 0, "start"), method = "radix")]
  if (length(spans) %% 2L == 1L) spans[[length(spans) + 1L]] <- NULL
  sub_ <- function(a, b) substr(s, a + 1L, b)
  for (p in seq(1L, length(spans), by = 2L)) {
    x <- spans[[p]]; y <- if (p + 1L <= length(spans)) spans[[p + 1L]] else NULL
    sl <- if (x$start - 10L < 0L) "" else "..."
    st <- max(x$start - 10L, 0L)
    el <- if (x$end + 10L > n) "" else "..."
    en <- min(x$end + 10L, n)
    if (x$is_end)
      stop(.latex_error(sprintf("Extra '%s' at index %d or missing corresponding '%s' in LaTeX string: %s%s%s",
                                x$text, x$start, d[[x$delim]][3L], sl, sub_(st, en), el)))
    if (is.null(y))
      stop(.latex_error(sprintf("Extra '%s' at index %d or missing corresponding '%s' in LaTeX string: %s%s%s",
                                x$text, x$start, d[[x$delim]][4L], sl, sub_(st, en), el)))
    if (!(y$is_end && y$delim == x$delim)) {
      sl <- if (x$start > 0L) "..." else ""
      el <- if (y$end < n) "..." else ""
      stop(.latex_error(sprintf("Expected '%s' to close the '%s' at index %d but found '%s' at index %d of LaTeX string instead: %s%s%s",
                                d[[x$delim]][4L], x$text, x$start, y$text, y$start, sl, sub_(x$start, y$end), el)))
    }
  }
  invisible()
}

#' Parse a LaTeX string into an expression
#'
#' `parse_latex()` converts a LaTeX math string to an expression with the
#' grammar and transformer of SymPy's lark LaTeX backend: an Earley parser
#' reads the whole grammar, including implicit products, fractions, roots,
#' integrals, sums, limits, derivatives in Leibniz notation, matrices and
#' relations. An input with several readings that give different
#' expressions returns all of them in a `lark_tree` with data `"_ambig"`,
#' printed as SymPy prints it.
#'
#' @param s A character string holding LaTeX.
#' @param strict Ignored; kept for the signature of SymPy's `parse_latex()`.
#' @param backend Only `"lark"` is available.
#' @return A `sym` expression, a matrix, or a `lark_tree` of the readings.
#' @examples
#' parse_latex("\\frac{1 + x}{\\sqrt{y}}")
#' parse_latex("\\int_0^1 x^2 dx")
#' parse_latex("x^2 + 2x + 1 = 0")
#' @export
parse_latex <- function(s, strict = FALSE, backend = "lark") {
  if (grepl("\\begin{cases}", s, fixed = TRUE))
    stop(.latex_error(paste0("The 'cases' environment is not currently supported by parse_latex. \n",
                             "If you're trying to parse a piecewise function or a block of expressions, \n",
                             "consider breaking them into separate parse_latex calls.")))
  .latex_check_matrix_delimiters(s)
  if (backend == "lark") return(parse_latex_lark(s))
  if (backend == "antlr") stop("the antlr backend of parse_latex() is not available; use backend = \"lark\"", call. = FALSE)
  stop(structure(class = c("NotImplementedError", "error", "condition"),
                 list(message = sprintf("Using the '%s' backend in the LaTeX parser is not supported, backend must be one of ('antlr', 'lark')", backend),
                      call = NULL)))
}

#' Parse LaTeX with the lark grammar
#'
#' The lark backend of [parse_latex()], with the readings of constants and
#' functions given in `locals`.
#'
#' @param s A character string holding LaTeX.
#' @param locals A named list: names are commands as written in `s` (the
#'   constants `\\pi`, `\\mathrm{e}`, `\\infty`, `\\imaginaryunit`, the
#'   trigonometric and hyperbolic functions, `\\exp`, `\\ln`, `\\log`, `\\lg`,
#'   or `\\operatorname{name}`); values are expressions for constants and
#'   functions (such as [Function()] results) for functions.
#' @return A `sym` expression, a matrix, or a `lark_tree` of the readings.
#' @examples
#' parse_latex_lark("2\\pi r")
#' parse_latex_lark("2\\pi r", locals = list("\\pi" = Symbol("pi")))
#' @export
parse_latex_lark <- function(s, locals = NULL) {
  ctx <- .latex_ctx(locals)
  .latex_check_commands(s)
  tree <- .lark_parse(s)
  .latex_transform_all(tree, ctx)
}

.latex_check_commands <- function(s) {
  m <- gregexpr("\\\\(?:\\\\|([a-zA-Z]+)|.)", s, perl = TRUE)[[1L]]
  if (m[1L] == -1L) return(invisible())
  cs <- attr(m, "capture.start")[, 1L]; cl <- attr(m, "capture.length")[, 1L]
  for (j in seq_along(m)) {
    if (cl[j] > 0L) {
      w <- substr(s, cs[j], cs[j] + cl[j] - 1L)
      if (!(w %in% .LARK_CONTROL_WORDS))
        stop(.latex_error(paste0("Unknown command ", substr(s, m[j], m[j] + attr(m, "match.length")[j] - 1L))))
    }
  }
  invisible()
}

#' @export
format.lark_tree <- function(x, ...) {
  inner <- vapply(x$children, function(c) if (is_sym(c) || inherits(c, "symmatrix")) sstr(c) else format(c), "")
  paste0("Tree('", x$data, "', [", paste(inner, collapse = ", "), "])")
}

#' @export
print.lark_tree <- function(x, ...) {
  cat(format(x), "\n", sep = "")
  invisible(x)
}

#' @export
format.lark_token <- function(x, ...) paste0("Token('", x$type, "', '", x$value, "')")

# ---------------------------------------------------------------------------
# TransformToSymPyExpr

.LATEX_CONSTANTS <- c("\\pi", "\\mathrm{e}", "\\infty", "\\imaginaryunit")
.LATEX_TRIG <- c("\\sin" = "sin", "\\cos" = "cos", "\\tan" = "tan", "\\csc" = "csc", "\\sec" = "sec",
                 "\\cot" = "cot", "\\arcsin" = "asin", "\\arccos" = "acos", "\\arctan" = "atan",
                 "\\arccsc" = "acsc", "\\arcsec" = "asec", "\\arccot" = "acot", "\\sinh" = "sinh",
                 "\\cosh" = "cosh", "\\tanh" = "tanh", "\\coth" = "coth", "\\sech" = "sech",
                 "\\csch" = "csch", "\\arsinh" = "asinh", "\\arcosh" = "acosh", "\\artanh" = "atanh",
                 "\\arcsinh" = "asinh", "\\arccosh" = "acosh", "\\arctanh" = "atanh")
.LATEX_INVERSE_TRIG <- c(sin = "asin", cos = "acos", tan = "atan", csc = "acsc", sec = "asec", cot = "acot",
                         sinh = "asinh", cosh = "acosh", tanh = "atanh", coth = "acoth", sech = "asech",
                         csch = "acsch")
.LATEX_EXPONENTIAL <- c("\\exp" = "exp", "\\ln" = "log", "\\log" = "log", "\\lg" = "lg")
.LATEX_NAMED <- c(exp = "exp", ln = "log", log = "log", min = "Min", max = "Max", sqrt = "sqrt", abs = "Abs",
                  sign = "sign", sgn = "sign", Re = "re", Im = "im", arg = "arg", erf = "erf", erfc = "erfc",
                  sinc = "sinc")

.latex_ctx <- function(locals) {
  locals <- if (is.null(locals)) list() else locals
  for (cmd in names(locals)) {
    if (!(cmd %in% .LATEX_CONSTANTS || cmd %in% names(.LATEX_TRIG) || cmd %in% names(.LATEX_EXPONENTIAL) ||
          grepl("^\\\\operatorname\\{[a-zA-Z]+\\}$", cmd)))
      stop(structure(class = c("ValueError", "error", "condition"),
                     list(message = sprintf("The reading of %s cannot be overridden.", cmd), call = NULL)))
  }
  list(locals = locals)
}

# A "function" of the transformer is list(kind, name, f): kind "sympy" for a
# SymPy function applied with fcall semantics, "user" for a locals override.
.ltx_fn <- function(name) list(kind = "sympy", name = name)
.ltx_user <- function(f) list(kind = "user", f = f)
.ltx_call <- function(fn, args) {
  if (fn$kind == "user") {
    f <- fn$f
    if (is_sym(f)) stop(structure(class = c("TypeError", "error", "condition"),
                                  list(message = sprintf("'%s' object is not callable", f$type), call = NULL)))
    return(do.call(f, args))
  }
  args <- lapply(args, .ltx_sym)
  switch(fn$name,
    lg = fcall("log", list(args[[1L]], .int(10L))),
    sqrt = .pow(args[[1L]], S_Half),
    Min = do.call(Min, args), Max = do.call(Max, args),
    det = .ltx_determinant(args[[1L]]), trace = .ltx_trace(args[[1L]]), adjugate = .ltx_adjugate(args[[1L]]),
    fcall(fn$name, args))
}
.ltx_same_fn <- function(a, b) a$kind == b$kind && identical(a$name, b$name) && identical(a$f, b$f)

.ltx_is_matrix <- function(x) inherits(x, "symmatrix") || .is_matexpr(x)
.ltx_is_derivop <- function(x) inherits(x, "latex_derivative_operator")
.ltx_is_faa <- function(x) inherits(x, "latex_factors_after_argument")
.ltx_sym <- function(x) {
  if (.ltx_is_derivop(x)) stop(.latex_error("A derivative must be applied to an expression."))
  if (inherits(x, "symmatrix")) return(x)
  as_sym(x)
}
.ltx_tok <- function(x, value = NULL) inherits(x, "lark_token") && (is.null(value) || x$value %in% value)
.ltx_index <- function(tokens, value, from = 1L) {
  for (j in seq.int(from, length(tokens))) if (.ltx_tok(tokens[[j]], value)) return(j)
  stop(structure(class = c("ValueError", "error", "condition"),
                 list(message = sprintf("'%s' is not in list", value), call = NULL)))
}
.ltx_in <- function(tokens, value) any(vapply(tokens, .ltx_tok, TRUE, value = value))
.ltx_add <- function(a, b) .add(list(.ltx_sym(a), .ltx_sym(b)))
.ltx_neg <- function(a) .neg(.ltx_sym(a))
.ltx_multiply <- function(lh, rh) {
  if (.ltx_is_matrix(lh) || .ltx_is_matrix(rh)) return(MatMul(lh, rh))
  .mul(list(.ltx_sym(lh), .ltx_sym(rh)))
}
.ltx_division <- function(num, den) {
  if (.ltx_is_matrix(den))
    stop(.latex_error(paste0("Cannot divide by matrices like this since it is not clear if left or right ",
                             "multiplication by the inverse is intended. Try explicitly multiplying by the ",
                             "inverse instead.")))
  if (.ltx_is_matrix(num)) return(MatMul(num, .pow(.ltx_sym(den), S_NegativeOne)))
  .mul(list(.ltx_sym(num), .pow(.ltx_sym(den), S_NegativeOne)))
}
.ltx_apply <- function(fn, argument) {
  if (.ltx_is_faa(argument)) return(.ltx_multiply(fn(argument$argument), argument$factors))
  fn(argument)
}
.ltx_greek <- function(x) gsub("var(?!pi)", "", x, perl = TRUE)
.ltx_split_sub <- function(v) {
  k <- regexpr("_", v, fixed = TRUE)
  c(substr(v, 1L, k - 1L), substr(v, k + 1L, nchar(v)))
}
.ltx_constant <- function(ctx, cmd) {
  if (!is.null(ctx$locals[[cmd]])) return(ctx$locals[[cmd]])
  switch(cmd, "\\pi" = S_Pi, "\\mathrm{e}" = S_E, "\\infty" = S_Infinity, "\\imaginaryunit" = S_I)
}
.ltx_function <- function(ctx, cmd) {
  if (!is.null(ctx$locals[[cmd]])) return(.ltx_user(ctx$locals[[cmd]]))
  if (cmd %in% names(.LATEX_TRIG)) return(.ltx_fn(.LATEX_TRIG[[cmd]]))
  .ltx_fn(.LATEX_EXPONENTIAL[[cmd]])
}
.ltx_named_function <- function(ctx, name) {
  key <- paste0("\\operatorname{", name, "}")
  if (!is.null(ctx$locals[[key]])) return(.ltx_user(ctx$locals[[key]]))
  if (name %in% names(.LATEX_NAMED)) return(.ltx_fn(.LATEX_NAMED[[name]]))
  if (paste0("\\", name) %in% names(.LATEX_TRIG)) return(.ltx_fn(.LATEX_TRIG[[paste0("\\", name)]]))
  if (name == "lg") return(.ltx_fn("lg"))
  if (name == "det") return(.ltx_fn("det"))
  if (name %in% c("tr", "trace")) return(.ltx_fn("trace"))
  if (name %in% c("adj", "adjugate")) return(.ltx_fn("adjugate"))
  list(kind = "undefined", name = name)
}
.ltx_fcall <- function(fn) {
  force(fn)
  function(a) {
    if (fn$kind == "undefined") return(fcall(fn$name, list(.ltx_sym(a))))
    .ltx_call(fn, list(a))
  }
}
.ltx_determinant <- function(m) {
  if (!.ltx_is_matrix(m)) stop(.latex_error("Cannot take determinant of non-matrix."))
  det(m)
}
.ltx_trace <- function(m) {
  if (!.ltx_is_matrix(m)) stop(.latex_error("Cannot take trace of non-matrix."))
  Trace(m)
}
.ltx_adjugate <- function(m) {
  if (!.ltx_is_matrix(m)) stop(.latex_error("Cannot take adjugate of non-matrix."))
  adjugate(doit(m))
}
.ltx_matrix_from <- function(body) {
  rows <- Filter(function(r) .is_tree(r, "matrix_row"), body$children)
  Matrix(lapply(rows, function(r) Filter(function(y) !.ltx_tok(y) || y$type != "MATRIX_COL_DELIM", r$children)))
}

.ltx_integral <- function(tokens, special) {
  ui <- if (.ltx_in(tokens, "_")) .ltx_index(tokens, "_") else NULL
  ci <- if (.ltx_in(tokens, "^")) .ltx_index(tokens, "^") else NULL
  lower <- if (!is.null(ui)) tokens[[ui + 1L]] else NULL
  upper <- if (!is.null(ci)) tokens[[ci + 1L]] else NULL
  if (special) {
    tup <- tokens[[length(tokens)]]
    integrand <- tup[[1L]]; var <- tup[[2L]]
  } else {
    ds <- NULL
    for (cand in c("d", "\\text{d}", "\\mathrm{d}")) if (.ltx_in(tokens, cand)) { ds <- cand; break }
    if (is.null(ds))
      stop(.latex_error(paste0("Differential symbol was not found in the expression.",
                               "Valid differential symbols are \"d\", \"\\text{d}, and \"\\mathrm{d}\".")))
    dvi <- .ltx_index(tokens, ds) + 1L
    var <- tokens[[dvi]]
  }
  if (!is.null(lower) && is.null(upper))
    stop(.latex_error("Lower bound for the integral was found, but upper bound was not found."))
  if (!is.null(upper) && is.null(lower))
    stop(.latex_error("Upper bound for the integral was found, but lower bound was not found."))
  if (!special) {
    integrand <- if (!is.null(ui) && ui == dvi - 3L) S_One
                 else if (!is.null(ci) && ci == dvi - 3L) S_One
                 else if (dvi == 3L) S_One
                 else tokens[[dvi - 2L]]
  }
  if (!is.null(lower)) return(Integral(.ltx_sym(integrand), list(var, .ltx_sym(lower), .ltx_sym(upper))))
  Integral(.ltx_sym(integrand), var)
}

.ltx_superscript <- function(tokens) {
  base <- tokens[[1L]]
  sup <- if (length(tokens) == 3L) tokens[[3L]] else tokens[[4L]]
  tt <- if (.ltx_tok(sup)) sup$type else ""
  nreps <- function(v, unit) nchar(v) / nchar(unit)
  if (.ltx_is_matrix(base)) {
    if (is_sym(sup) && identical(sup$key, Symbol("T")$key)) return(Transpose(base))
    if (is_sym(sup) && identical(sup$key, Symbol("H")$key)) return(adjoint(base))
    if (tt == "PRIMES") return(if (nchar(sup$value) %% 2 == 0) base else Transpose(base))
    if (tt %in% c("PRIMES_VIA_CMD", "CMD_PRIME"))
      return(if (nreps(sup$value, "\\prime") %% 2 == 0) base else Transpose(base))
    if (tt == "STARS") return(if (nchar(sup$value) %% 2 == 0) doit(base) else adjoint(base))
    if (tt %in% c("STARS_VIA_CMD", "CMD_ASTERISK"))
      return(if (nreps(sup$value, "\\ast") %% 2 == 0) doit(base) else adjoint(base))
  }
  if (tt %in% c("PRIMES", "PRIMES_VIA_CMD", "CMD_PRIME", "STARS", "STARS_VIA_CMD", "CMD_ASTERISK"))
    stop(.latex_error(sprintf("%s with superscript %s is not understood.", .ltx_str(base), sup$value)))
  if (.ltx_is_matrix(base)) return(MatPow(base, .ltx_sym(sup)))
  .pow(.ltx_sym(base), .ltx_sym(sup))
}
.ltx_str <- function(x) if (is_sym(x) || inherits(x, "symmatrix")) sstr(x) else format(x)

.ltx_leibniz <- function(tokens) {
  d <- Symbol("d")
  n <- length(tokens)
  pd <- tokens[[1L]]; fun <- if (n > 2L) tokens[2:(n - 1L)] else list(); pv <- tokens[[n]]
  if (is_sym(pd) && is_sym(pv)) {
    be <- .ltx_base_exp(pd); bv <- .ltx_base_exp(pv)
    order_ <- be[[2L]]
    if (identical(be[[1L]]$key, d$key) && bv[[1L]]$type == "Symbol" && identical(order_$key, bv[[2L]]$key) &&
        ((order_$type == "Integer" && isTRUE(is_positive(order_))) || order_$type == "Symbol")) {
      op <- structure(list(variable = bv[[1L]], order = order_), class = "latex_derivative_operator")
      if (length(fun)) return(.ltx_derive(op, fun[[1L]]))
      return(op)
    }
  }
  numer <- if (length(fun)) .ltx_multiply(pd, fun[[1L]]) else pd
  .ltx_division(numer, .ltx_multiply(d, pv))
}
.ltx_base_exp <- function(e) {
  if (e$type == "Pow") return(e$args)
  if (e$type == "Function" && identical(e$name, "exp")) return(list(S_E, e$args[[1L]]))
  list(e, S_One)
}
.ltx_derive <- function(op, e) Derivative(.ltx_sym(e), list(op$variable, op$order))

.latex_rules <- list(
  "_ambig" = function(t, ctx) {
    readings <- list()
    for (v in t) if (!any(vapply(readings, .latex_same, TRUE, b = v))) readings[[length(readings) + 1L]] <- v
    if (length(readings) == 1L) readings[[1L]] else .lark_tree("_ambig", readings)
  },
  latex_string = function(t, ctx) {
    if (.ltx_is_derivop(t[[1L]])) stop(.latex_error("A derivative must be applied to an expression."))
    t[[1L]]
  },
  group_round_parentheses = function(t, ctx) t[[2L]],
  group_square_brackets = function(t, ctx) t[[2L]],
  group_curly_parentheses = function(t, ctx) t[[2L]],
  group_curly_literal_parentheses = function(t, ctx) t[[2L]],
  eq = function(t, ctx) Eq(.ltx_sym(t[[1L]]), .ltx_sym(t[[3L]])),
  ne = function(t, ctx) Ne(.ltx_sym(t[[1L]]), .ltx_sym(t[[3L]])),
  lt = function(t, ctx) Lt(.ltx_sym(t[[1L]]), .ltx_sym(t[[3L]])),
  lte = function(t, ctx) Le(.ltx_sym(t[[1L]]), .ltx_sym(t[[3L]])),
  gt = function(t, ctx) Gt(.ltx_sym(t[[1L]]), .ltx_sym(t[[3L]])),
  gte = function(t, ctx) Ge(.ltx_sym(t[[1L]]), .ltx_sym(t[[3L]])),
  relation_chain = function(t, ctx) {
    rels <- lapply(seq(2L, length(t), by = 2L), function(i) {
      a <- .ltx_sym(t[[i - 1L]]); b <- .ltx_sym(t[[i + 1L]])
      switch(t[[i]]$type, EQUAL = Eq(a, b), NOT_EQUAL = Ne(a, b), LT = Lt(a, b), LTE = Le(a, b),
             GT = Gt(a, b), GTE = Ge(a, b))
    })
    do.call(And, rels)
  },
  number_set = function(t, ctx) {
    name <- if (is_sym(t[[3L]])) sstr(t[[3L]]) else t[[3L]]$value
    sets <- list(N = "Naturals", Z = "Integers", Q = "Rationals", R = "Reals", C = "Complexes")
    if (is.null(sets[[name]])) stop(.latex_error(sprintf("Unknown number set \\mathbb{%s}.", name)))
    get(sets[[name]], envir = asNamespace("symcompR"))
  },
  membership = function(t, ctx) {
    cont <- Contains(t[[1L]], t[[3L]])
    if (t[[2L]]$type == "CMD_NOTIN") Not(cont) else cont
  },
  add = function(t, ctx) {
    if (length(t) == 2L) return(t[[2L]])
    lh <- t[[1L]]; rh <- t[[3L]]
    if (.ltx_is_matrix(lh) || .ltx_is_matrix(rh)) return(MatAdd(lh, rh))
    .ltx_add(lh, rh)
  },
  sub = function(t, ctx) {
    if (length(t) == 2L) {
      x <- t[[2L]]
      if (.ltx_is_faa(x)) return(structure(list(argument = .ltx_neg(x$argument), factors = x$factors),
                                          class = "latex_factors_after_argument"))
      if (.ltx_is_matrix(x)) return(MatMul(S_NegativeOne, x))
      return(.ltx_neg(x))
    }
    lh <- t[[1L]]; rh <- t[[3L]]
    if (.ltx_is_matrix(lh) || .ltx_is_matrix(rh)) return(MatAdd(lh, MatMul(S_NegativeOne, rh)))
    .ltx_add(lh, .ltx_neg(rh))
  },
  unary_sign = function(t, ctx) if (t[[1L]]$type == "SUB") .latex_rules$sub(t, ctx) else t[[2L]],
  mul = function(t, ctx) .ltx_multiply(t[[1L]], t[[3L]]),
  div = function(t, ctx) .ltx_division(t[[1L]], t[[3L]]),
  adjacent_expressions = function(t, ctx) {
    if (.ltx_is_derivop(t[[1L]])) return(.ltx_derive(t[[1L]], t[[2L]]))
    .ltx_multiply(t[[1L]], t[[2L]])
  },
  factors_after_argument = function(t, ctx)
    structure(list(argument = t[[1L]], factors = t[[2L]]), class = "latex_factors_after_argument"),
  superscript = function(t, ctx) .ltx_superscript(t),
  degree = function(t, ctx) .mul(list(.ltx_sym(t[[1L]]), .div(S_Pi, .int(180L)))),
  matrix_prime = function(t, ctx) {
    base <- t[[1L]]; primes <- t[[2L]]$value
    if (!.ltx_is_matrix(base)) stop(.latex_error(sprintf("(%s)%s is not understood.", .ltx_str(base), primes)))
    if (nchar(primes) %% 2L == 0L) base else Transpose(base)
  },
  power_with_subscript = function(t, ctx) {
    sub <- gsub(" ", "", gsub("^[{}]+|[{}]+$", "", substring(t[[4L]]$value, 2L)), fixed = TRUE)
    .pow(Symbol(sprintf("%s_{%s}", t[[1L]]$name, sub)), .ltx_sym(t[[3L]]))
  },
  symbol_prime = function(t, ctx) Symbol(paste0(t[[1L]]$name, t[[2L]]$value)),
  fraction = function(t, ctx) .ltx_division(t[[2L]], t[[3L]]),
  leibniz_fraction = function(t, ctx) .ltx_leibniz(t),
  binomial = function(t, ctx) fcall("binomial", list(.ltx_sym(t[[2L]]), .ltx_sym(t[[3L]]))),
  normal_integral = function(t, ctx) .ltx_integral(t, FALSE),
  group_curly_parentheses_int = function(t, ctx) {
    if (length(t) == 3L) list(S_One, t[[2L]]) else list(t[[2L]], t[[3L]])
  },
  special_fraction = function(t, ctx) {
    nv <- t[[2L]]
    list(.mul(list(.ltx_sym(nv[[1L]]), .pow(.ltx_sym(t[[3L]]), S_NegativeOne))), nv[[2L]])
  },
  integral_with_special_fraction = function(t, ctx) .ltx_integral(t, TRUE),
  group_curly_parentheses_special = function(t, ctx) {
    ui <- .ltx_index(t, "_"); ci <- .ltx_index(t, "^")
    lb <- .ltx_index(t, "{", ui); rb <- .ltx_index(t, "}", ui)
    bottom <- t[seq_len(rb - lb - 1L) + lb]
    top <- t[seq.int(ci + 1L, length(t))]
    list(bottom[[1L]], bottom[[length(bottom)]], top[[1L]])
  },
  summation = function(t, ctx) Sum(.ltx_sym(t[[3L]]), lapply(t[[2L]], .ltx_sym)),
  product = function(t, ctx) Product(.ltx_sym(t[[3L]]), lapply(t[[2L]], .ltx_sym)),
  limit_dir_expr = function(t, ctx) {
    ci <- .ltx_index(t, "^")
    dirn <- if (.ltx_in(t, "{")) t[[.ltx_index(t, "{", ci) + 1L]] else t[[ci + 1L]]
    v <- if (.ltx_tok(dirn, "+")) "+" else if (.ltx_tok(dirn, "-")) "-" else "+-"
    structure(list(t[[1L]], v), class = "latex_tuple")
  },
  group_curly_parentheses_lim = function(t, ctx) {
    if (inherits(t[[4L]], "latex_tuple")) list(t[[2L]], t[[4L]][[1L]], t[[4L]][[2L]])
    else list(t[[2L]], t[[4L]], "+-")
  },
  limit = function(t, ctx) {
    v <- t[[3L]]
    Limit(.ltx_sym(t[[length(t)]]), v[[1L]], .ltx_sym(v[[2L]]), v[[3L]])
  },
  differential = function(t, ctx) t[[2L]],
  list_of_expressions = function(t, ctx) {
    if (length(t) == 1L) return(t)
    Filter(function(a) {
      if (.ltx_tok(a)) {
        if (a$type != "COMMA") stop(.latex_error("A comma token was expected, but some other token was encountered."))
        return(FALSE)
      }
      TRUE
    }, t)
  },
  function_applied = function(t, ctx) fcall(sstr(t[[1L]]), lapply(t[[3L]], .ltx_sym)),
  inverse_function_applied = function(t, ctx)
    fcall(sprintf("%s^{-1}", sstr(t[[1L]])), lapply(t[[length(t) - 1L]], .ltx_sym)),
  min = function(t, ctx) do.call(Min, lapply(t[[3L]], .ltx_sym)),
  max = function(t, ctx) do.call(Max, lapply(t[[3L]], .ltx_sym)),
  bra = function(t, ctx) .latex_quantum("Bra"),
  ket = function(t, ctx) .latex_quantum("Ket"),
  inner_product = function(t, ctx) .latex_quantum("InnerProduct"),
  trigonometric_function = function(t, ctx) .ltx_apply(.ltx_fcall(.ltx_function(ctx, t[[1L]]$value)), t[[2L]]),
  trigonometric_function_power = function(t, ctx) {
    fn <- .ltx_function(ctx, t[[1L]]$value)
    ex <- t[[3L]]
    if (is_sym(ex) && identical(ex$key, S_NegativeOne$key) && fn$kind == "sympy" && fn$name %in% names(.LATEX_INVERSE_TRIG))
      return(.ltx_apply(.ltx_fcall(.ltx_fn(.LATEX_INVERSE_TRIG[[fn$name]])), t[[length(t)]]))
    f1 <- .ltx_fcall(fn)
    .ltx_apply(function(a) .pow(.ltx_sym(f1(a)), .ltx_sym(ex)), t[[length(t)]])
  },
  named_function = function(t, ctx) {
    name <- sub("^\\\\operatorname\\s*\\{\\s*([a-zA-Z]+)\\s*\\}.*$", "\\1", t[[1L]]$value, perl = TRUE)
    fn <- .ltx_named_function(ctx, name)
    f1 <- .ltx_fcall(fn)
    if (length(t) == 2L) return(.ltx_apply(f1, t[[2L]]))
    if (.ltx_tok(t[[2L]]) && t[[2L]]$type == "CARET") {
      ex <- t[[3L]]
      if (is_sym(ex) && identical(ex$key, S_NegativeOne$key) && fn$kind == "sympy" && fn$name %in% names(.LATEX_INVERSE_TRIG))
        return(.ltx_apply(.ltx_fcall(.ltx_fn(.LATEX_INVERSE_TRIG[[fn$name]])), t[[4L]]))
      return(.ltx_apply(function(a) .pow(.ltx_sym(f1(a)), .ltx_sym(ex)), t[[4L]]))
    }
    if (!((fn$kind == "sympy" && fn$name %in% c("Min", "Max")) || fn$kind %in% c("user", "undefined")))
      stop(.latex_error(sprintf("\\operatorname{%s} takes one argument.", name)))
    args <- t[seq(3L, length(t) - 1L, by = 2L)]
    if (fn$kind == "undefined") return(fcall(fn$name, lapply(args, .ltx_sym)))
    .ltx_call(fn, args)
  },
  abs = function(t, ctx) fcall("Abs", list(.ltx_sym(t[[2L]]))),
  floor = function(t, ctx) fcall("floor", list(.ltx_sym(t[[2L]]))),
  ceil = function(t, ctx) fcall("ceiling", list(.ltx_sym(t[[2L]]))),
  factorial = function(t, ctx) fcall("factorial", list(.ltx_sym(t[[1L]]))),
  conjugate = function(t, ctx) {
    a <- t[[2L]]
    if (inherits(a, "symmatrix")) return(.mat_map(a, function(z) fcall("conjugate", list(z))))
    fcall("conjugate", list(.ltx_sym(a)))
  },
  square_root = function(t, ctx) {
    if (length(t) == 2L) return(.pow(.ltx_sym(t[[2L]]), S_Half))
    root(.ltx_sym(t[[3L]]), .ltx_sym(t[[2L]]))
  },
  exponential = function(t, ctx) .ltx_apply(.ltx_fcall(.ltx_function(ctx, t[[1L]]$value)), t[[2L]]),
  function_power = function(t, ctx) {
    f1 <- .ltx_fcall(.ltx_function(ctx, t[[1L]]$value))
    .ltx_apply(function(a) .pow(.ltx_sym(f1(a)), .ltx_sym(t[[3L]])), t[[4L]])
  },
  log = function(t, ctx) {
    fn <- .ltx_function(ctx, t[[1L]]$value)
    if (.ltx_in(t, "_")) {
      base <- t[[3L]]
      return(.ltx_apply(function(a) .ltx_call(fn, list(a, base)), t[[4L]]))
    }
    .ltx_apply(.ltx_fcall(fn), t[[2L]])
  },
  matrix = function(t, ctx) .ltx_matrix_from(t[[2L]]),
  determinant = function(t, ctx) {
    if (length(t) == 2L) return(.ltx_apply(.ltx_determinant, t[[2L]]))
    det(.ltx_matrix_from(t[[2L]]))
  },
  trace = function(t, ctx) .ltx_apply(.ltx_trace, t[[2L]]),
  adjugate = function(t, ctx) .ltx_apply(.ltx_adjugate, t[[2L]]),
  multi_letter_symbol = function(t, ctx) {
    if (length(t) == 4L) return(Symbol(t[[3L]]$value))
    Symbol(paste0(t[[3L]]$value, t[[5L]]$value))
  },
  number = function(t, ctx) {
    tk <- t[[1L]]
    if (tk$type == "CMD_IMAGINARY_UNIT") return(.ltx_constant(ctx, tk$value))
    if (grepl(".", tk$value, fixed = TRUE)) Float(tk$value) else make_rational(gmp::as.bigq(.str_bigz(tk$value)))
  }
)

.latex_quantum <- function(cls) {
  stop(structure(class = c("NotImplementedError", "error", "condition"),
                 list(message = sprintf("sympy.physics.quantum.%s is not available", cls), call = NULL)))
}

.latex_tokens <- list(
  SYMBOL = function(v, ctx) Symbol(v),
  DIGIT = function(v, ctx) make_rational(gmp::as.bigq(.str_bigz(v))),
  CMD_PI = function(v, ctx) .ltx_constant(ctx, v),
  EULER_NUMBER = function(v, ctx) .ltx_constant(ctx, v),
  CMD_INFTY = function(v, ctx) .ltx_constant(ctx, v),
  GREEK_SYMBOL_WITH_PRIMES = function(v, ctx) Symbol(.ltx_greek(substring(v, 2L))),
  SUBSCRIPTED_GREEK_SYMBOL = function(v, ctx) Symbol(.ltx_greek(substring(v, 2L))),
  LATIN_SYMBOL_WITH_LATIN_SUBSCRIPT = function(v, ctx) {
    p <- .ltx_split_sub(v)
    sub <- if (startsWith(p[2L], "{")) gsub(" ", "", substr(p[2L], 2L, nchar(p[2L]) - 1L), fixed = TRUE) else p[2L]
    Symbol(sprintf("%s_{%s}", p[1L], sub))
  },
  GREEK_SYMBOL_WITH_LATIN_SUBSCRIPT = function(v, ctx) {
    p <- .ltx_split_sub(v)
    g <- .ltx_greek(substring(p[1L], 2L))
    sub <- if (startsWith(p[2L], "{")) gsub(" ", "", substr(p[2L], 2L, nchar(p[2L]) - 1L), fixed = TRUE) else p[2L]
    Symbol(sprintf("%s_{%s}", g, sub))
  },
  LATIN_SYMBOL_WITH_GREEK_SUBSCRIPT = function(v, ctx) {
    p <- .ltx_split_sub(v)
    g <- if (startsWith(p[2L], "{")) substr(p[2L], 3L, nchar(p[2L]) - 1L) else substring(p[2L], 2L)
    Symbol(sprintf("%s_{%s}", p[1L], .ltx_greek(g)))
  },
  GREEK_SYMBOL_WITH_GREEK_SUBSCRIPT = function(v, ctx) {
    p <- .ltx_split_sub(v)
    gb <- .ltx_greek(substring(p[1L], 2L))
    gs <- if (startsWith(p[2L], "{")) substr(p[2L], 3L, nchar(p[2L]) - 1L) else substring(p[2L], 2L)
    Symbol(sprintf("%s_{%s}", gb, .ltx_greek(gs)))
  }
)

# TransformToSymPyExpr.transform: mark the Leibniz fractions, then transform
# bottom up; a reading that fails fails the whole parse, as in lark.
.latex_mark_leibniz <- function(node) {
  if (!.is_tree(node)) return(node)
  node$children <- lapply(node$children, .latex_mark_leibniz)
  if (node$data != "fraction") return(node)
  content <- function(a) if (.is_tree(a, "group_curly_parentheses")) a$children[[2L]] else a
  is_d <- function(n) .ltx_tok(n, "d") && n$type == "SYMBOL"
  is_power_of_d <- function(n) is_d(n) || (.is_tree(n, "superscript") && is_d(n$children[[1L]]))
  is_product <- function(n) .is_tree(n, "adjacent_expressions")
  num <- content(node$children[[2L]]); den <- content(node$children[[3L]])
  if (!(is_product(den) && is_d(den$children[[1L]]))) return(node)
  if (is_power_of_d(num)) factors <- list(num)
  else if (is_product(num) && is_power_of_d(num$children[[1L]])) factors <- num$children
  else return(node)
  .lark_tree("leibniz_fraction", c(factors, list(den$children[[2L]])))
}

.latex_transform_all <- function(tree, ctx) .latex_transform(.latex_mark_leibniz(tree), ctx)

.latex_same <- function(a, b) {
  if (is_sym(a) && is_sym(b)) return(identical(a$key, b$key))
  identical(a, b)
}

.latex_visit_error <- function(rule, e) {
  if (inherits(e, "LaTeXParsingError")) stop(e)
  cls <- setdiff(class(e), c("simpleError", "error", "condition", "sympy_error"))
  msg <- sprintf("Error trying to process rule \"%s\":\n\n%s", rule, conditionMessage(e))
  stop(structure(class = c("VisitError", "error", "condition"), list(message = msg, call = NULL, orig = cls)))
}

.latex_transform <- function(node, ctx) {
  if (inherits(node, "lark_token")) {
    h <- .latex_tokens[[node$type]]
    if (is.null(h)) return(node)
    return(tryCatch(h(node$value, ctx), error = function(e) .latex_visit_error(node$type, e)))
  }
  kids <- lapply(node$children, .latex_transform, ctx = ctx)
  h <- .latex_rules[[node$data]]
  if (is.null(h)) return(.lark_tree(node$data, kids))
  tryCatch(h(kids, ctx), error = function(e) .latex_visit_error(node$data, e))
}

#' All the readings of an ambiguous LaTeX string
#'
#' `parse_latex_lark_all()` returns the list of the expressions `s` may mean
#' under the lark grammar, with one element when `s` is not ambiguous. When
#' the same ambiguity occurs several times, the readings choose the same way
#' at every place where it occurs, as in SymPy's `parse_latex_lark_all`.
#'
#' @inheritParams parse_latex_lark
#' @return A list of expressions.
#' @examples
#' parse_latex_lark_all("\\sin xy")
#' @export
parse_latex_lark_all <- function(s, locals = NULL) {
  ctx <- .latex_ctx(locals)
  .latex_check_commands(s)
  tree <- .lark_parse(s)
  trees <- if (.is_tree(tree, "_ambig")) tree$children else list(tree)
  trees <- .lark_parallel_readings(trees)
  out <- list()
  for (t in trees) {
    v <- .latex_transform_all(t, ctx)
    if (!any(vapply(out, .latex_same, TRUE, b = v))) out[[length(out) + 1L]] <- v
  }
  out
}

# _node_spans, _signature, _choice and _parallel_readings of
# sympy/parsing/latex/lark/latex_parser.py. A tree is flattened into a table
# of nodes; spans are character offsets.
.lark_flatten <- function(tree) {
  nodes <- list()
  visit <- function(node) {
    if (inherits(node, "lark_token")) {
      id <- length(nodes) + 1L
      nodes[[id]] <<- list(tok = TRUE, type = node$type, value = node$value,
                           span = c(node$start, node$start + nchar(node$value)), kids = integer())
      return(id)
    }
    kids <- vapply(node$children, visit, 0L)
    id <- length(nodes) + 1L
    nodes[[id]] <<- list(tok = FALSE, data = node$data, kids = kids,
                         span = c(nodes[[kids[1L]]]$span[1L], nodes[[kids[length(kids)]]]$span[2L]))
    id
  }
  root <- visit(tree)
  spans <- vapply(nodes, function(n) paste(n$span, collapse = ":"), "")
  # node_at: the outermost node covering each span (the last one visited)
  node_at <- tapply(seq_along(nodes), spans, max)
  list(nodes = nodes, spans = spans, node_at = node_at, root = root)
}

.lark_signature <- function(fl, id, sites) {
  n <- fl$nodes[[id]]
  enc <- function(x) paste0(nchar(x), ":", x)
  if (n$tok) return(paste0("T", enc(n$type), enc(n$value)))
  sp <- fl$spans[id]
  parts <- vapply(n$kids, function(k) {
    if (fl$spans[k] != sp && fl$spans[k] %in% sites) "N" else enc(.lark_signature(fl, k, sites))
  }, "")
  paste0("R", enc(n$data), "(", paste(parts, collapse = ""), ")")
}

.lark_choice <- function(fl, id) {
  sp <- fl$spans[id]
  labels <- character()
  repeat {
    n <- fl$nodes[[id]]
    if (n$tok) return(c(labels, n$type))
    labels <- c(labels, n$data)
    inner <- n$kids[fl$spans[n$kids] == sp]
    if (!length(inner)) return(labels)
    id <- inner[1L]
  }
}

.lark_parallel_readings <- function(trees) {
  if (length(trees) < 2L) return(trees)
  fls <- lapply(trees, .lark_flatten)
  common <- Reduce(intersect, lapply(fls, function(f) names(f$node_at)))
  len <- vapply(common, function(s) { p <- as.integer(strsplit(s, ":", fixed = TRUE)[[1L]]); p[2L] - p[1L] }, 0L)
  common <- common[order(len)]
  site_spans <- character()
  site_kind <- character()
  site_choices <- list()
  for (sp in common) {
    ids <- vapply(fls, function(f) unname(f$node_at[[sp]]), 0L)
    sigs <- vapply(seq_along(fls), function(i) .lark_signature(fls[[i]], ids[i], site_spans), "")
    if (length(unique(sigs)) > 1L) {
      ch <- lapply(seq_along(fls), function(i) .lark_choice(fls[[i]], ids[i]))
      k <- 0L
      repeat {
        if (any(lengths(ch) <= k)) break
        if (length(unique(vapply(ch, `[`, "", k + 1L))) != 1L) break
        k <- k + 1L
      }
      ch <- vapply(ch, function(c1) paste(c1[seq_len(length(c1) - k) + k], collapse = "\r"), "")
      kind <- if (length(unique(ch)) == length(unique(sigs))) paste0("K", paste(sort(unique(ch), method = "radix"), collapse = "\t"))
              else paste0("S", sp)
      site_spans <- c(site_spans, sp)
      site_kind <- c(site_kind, kind)
      site_choices[[length(site_choices) + 1L]] <- ch
    }
  }
  keep <- vapply(seq_along(trees), function(i) {
    all(vapply(unique(site_kind), function(kd) {
      length(unique(vapply(site_choices[site_kind == kd], `[`, "", i))) == 1L
    }, TRUE))
  }, TRUE)
  if (any(keep)) trees[keep] else trees
}
