# LaTeX and code printers, lambdify and common subexpression elimination,
# after sympy/printing and sympy/simplify/cse_main.py.

# split_super_sub() of sympy/printing/conventions.py: the name, then parts
# opened by "^" or "__" (superscripts) or "_" (subscripts); trailing digits
# of the name are a first subscript
.split_super_sub <- function(text) {
  if (!nzchar(text)) return(list(name = text, supers = character(), subs = character()))
  pos <- 1L
  n <- nchar(text)
  name <- NULL
  supers <- subs <- character()
  while (pos <= n) {
    start <- pos + 1L
    if (substr(text, pos, pos + 1L) == "__") start <- start + 1L
    rest <- substr(text, start, n)
    ph <- regexpr("^", rest, fixed = TRUE)
    pu <- regexpr("_", rest, fixed = TRUE)
    nxt <- min(if (ph > 0L) start + ph - 1L else n + 1L, if (pu > 0L) start + pu - 1L else n + 1L)
    part <- substr(text, pos, nxt - 1L)
    pos <- nxt
    if (is.null(name)) name <- part
    else if (startsWith(part, "^")) supers <- c(supers, substring(part, 2L))
    else if (startsWith(part, "__")) supers <- c(supers, substring(part, 3L))
    else subs <- c(subs, substring(part, 2L))
  }
  m <- regmatches(name, regexec("^([[:alpha:]]+)([0-9]+)$", name))[[1L]]
  if (length(m)) {
    name <- m[2L]
    subs <- c(m[3L], subs)
  }
  list(name = name, supers = supers, subs = subs)
}

# ---------------------------------------------------------------------------
# Code printers

.code_styles <- list(
  c = list(pow = function(b, e) paste0("pow(", b, ", ", e, ")"), sqrt = "sqrt", pi = "M_PI", e = "M_E",
           num = function(q) paste0(gmp::numerator(q), ".0/", gmp::denominator(q), ".0"), abs = "fabs",
           ln = "log", ceiling = "ceil", I = "I", inf = "HUGE_VAL"),
  fortran = list(pow = function(b, e) paste0(b, "**", e), sqrt = "sqrt", pi = "3.14159265358979d0", e = "2.71828182845905d0",
                 num = function(q) paste0(gmp::numerator(q), ".0d0/", gmp::denominator(q), ".0d0"), abs = "abs",
                 ln = "log", ceiling = "ceiling", I = "cmplx(0,1)", inf = "huge(0d0)", inv = "1d0/"),
  js = list(pow = function(b, e) paste0("Math.pow(", b, ", ", e, ")"), sqrt = "Math.sqrt", pi = "Math.PI", e = "Math.E",
            num = function(q) paste0(gmp::numerator(q), "/", gmp::denominator(q)), abs = "Math.abs",
            ln = "Math.log", ceiling = "Math.ceil", I = "1i", inf = "Number.POSITIVE_INFINITY", prefix = "Math.", inv = "1/"),
  octave = list(plain_frac = TRUE, pow = function(b, e) paste0(b, "^", e), sqrt = "sqrt", pi = "pi", e = "e",
                num = function(q) paste0(gmp::numerator(q), "/", gmp::denominator(q)), abs = "abs",
                ln = "log", ceiling = "ceil", I = "1i", inf = "inf", inv = "1./"),
  r = list(pow = function(b, e) paste0(b, "^", e), sqrt = "sqrt", pi = "pi", e = "exp(1)",
           num = function(q) paste0(gmp::numerator(q), ".0/", gmp::denominator(q), ".0"), abs = "abs",
           ln = "log", ceiling = "ceiling", I = "1i", inf = "Inf"),
  python = list(pow = function(b, e) paste0(b, "**", e), sqrt = "math.sqrt", pi = "math.pi", e = "math.e",
                num = function(q) paste0(gmp::numerator(q), "/", gmp::denominator(q)), abs = "abs",
                ln = "math.log", ceiling = "math.ceil", I = "1j", inf = "math.inf", prefix = "math.", inv = "1/")
)

# A Float prints as SymPy's printers do: str(evalf(17)) with the zeros
# stripped in C and C++; str() in Python and Octave, with a d exponent in
# Fortran (zeros are stripped below the top level); 17 digits elsewhere
.code_float <- function(e, st) {
  if (st$pi %in% c("math.pi", "3.14159265358979d0") || identical(st$e, "e")) {
    s <- .fmt_float_sym(e, strip = !isTRUE(st$top))
    if (identical(st$pi, "3.14159265358979d0")) s <- if (grepl("e", s, fixed = TRUE)) sub("e", "d", s, fixed = TRUE) else paste0(s, "d0")
    return(s)
  }
  if (!identical(st$pi, "M_PI")) return(if (.float_prec(e) <= 53) sprintf("%.17g", e$value) else .fmt_float_sym(e, strip = TRUE))
  s <- .mp_to_str(.float_mp(e), 17, strip_zeros = FALSE)
  if (!grepl("[e.]", s)) s <- paste0(s, ".0")
  p <- strsplit(s, "e", fixed = TRUE)[[1L]]
  p[1L] <- sub("0+$", "", p[1L])
  if (endsWith(p[1L], ".")) p[1L] <- paste0(p[1L], "0")
  paste(p, collapse = "e")
}

.code <- function(e, st) {
  t <- e$type
  atomize <- function(z) {
    s <- .code(z, st)
    if (.prec(z) < PREC$Atom) paste0("(", s, ")") else s
  }
  mulp <- function(z) {
    s <- .code(z, st)
    if (.prec(z) < PREC$Mul) paste0("(", s, ")") else s
  }
  switch(t,
         Integer = as.character(gmp::numerator(e$value)),
         Rational = st$num(e$value),
         Float = .code_float(e, st),
         Symbol = , Dummy = e$name,
         Constant = switch(EXPR = e$name, pi = st$pi, E = st$e, I = st$I, oo = st$inf, `-oo` = paste0("-", st$inf),
                           GoldenRatio = "((1 + sqrt(5))/2)", EulerGamma = "0.57721566490153286",
                           stop("cannot print constant ", e$name)),
         Add = {
           terms <- .add_print_order(e$args)
           out <- character(length(terms))
           for (i in seq_along(terms)) {
             z <- terms[[i]]
             if (i > 1L && could_extract_minus_sign(z)) out[i] <- paste0(" - ", mulp(.neg(z)))
             else out[i] <- paste0(if (i > 1L) " + ", .code(z, st))
           }
           paste(out, collapse = "")
         },
         Mul = {
           if (could_extract_minus_sign(e)) return(paste0("-", mulp(.neg(e))))
           c1 <- e$args[[1L]]
           if (c1$type == "Rational" && !isTRUE(st$plain_frac)) {
             # (p/q)*rest as SymPy prints, which also avoids integer division in C
             rest <- if (length(e$args) == 2L) e$args[[2L]] else new_sym("Mul", e$args[-1L])
             return(paste0("(", st$num(c1$value), ")*", mulp(rest)))
           }
           nd <- fraction(e, exact = TRUE)
           fmt <- function(z) {
             fs <- if (z$type == "Mul") .mul_print_order(z$args) else list(z)
             paste(vapply(fs, mulp, ""), collapse = "*")
           }
           if (is_one_num(nd[[2L]])) return(fmt(nd[[1L]]))
           den <- nd[[2L]]
           ds <- if (den$type == "Mul") paste0("(", .code(den, st), ")") else mulp(den)
           n <- nd[[1L]]
           ns <- if (is_one_num(n)) "1" else fmt(n)
           paste0(ns, "/", ds)
         },
         Pow = {
           b <- e$args[[1L]]; x <- e$args[[2L]]
           if (identical(x$key, S_Half$key)) return(paste0(st$sqrt, "(", .code(b, st), ")"))
           if (identical(x$key, S_NegativeOne$key)) return(paste0(st$inv %||% "1.0/", atomize(b)))
           st$pow(atomize(b), if (x$type == "Rational") paste0("(", st$num(x$value), ")") else atomize(x))
         },
         UnevaluatedExpr = .code(e$args[[1L]], st),
         Function = {
           d <- .FUNS[[e$name]]
           if (!is.null(d$code)) return(paste0("(", .code(d$code(e), st), ")"))
           args <- paste(vapply(e$args, .code, "", st = st), collapse = ", ")
           nm <- switch(EXPR = e$name, Abs = st$abs, log = st$ln, ceiling = st$ceiling,
                        paste0(if (!is.null(st$prefix)) st$prefix, e$name))
           paste0(nm, "(", args, ")")
         },
         stop("cannot generate code for ", t))
}

#' Code generation
#'
#' Print an expression as C (`ccode`), Fortran (`fcode`), JavaScript
#' (`jscode`), Octave/MATLAB (`octave_code`), R (`rcode`) or Python
#' (`pycode`) source.
#'
#' @param e A `sym` expression.
#' @param assign_to Optional name of a variable to assign.
#' @param ... Printer settings, as in SymPy (for example `standard`,
#'   `precision`, `human`, `contract` or `user_functions`).
#' @return A character scalar.
#' @examples
#' x <- Symbol("x")
#' ccode(x^2 + sin(x) / 3)
#' rcode(sqrt(x + 1), assign_to = "y")
#' @name code-printers
NULL

.mk_printer <- function(style, fmt) {
  force(style); force(fmt)
  function(e, assign_to = NULL, ...) {
    eng <- .cg_engine(style)
    if (!is.null(eng)) return(eng(e, assign_to, ...))
    e <- as_sym(e)
    st <- .code_styles[[style]]
    st$top <- e$type == "Float"
    s <- .code(e, st)
    if (!is.null(assign_to)) s <- sprintf(fmt, assign_to, s)
    s
  }
}

#' @rdname code-printers
#' @export
ccode <- .mk_printer("c", "%s = %s;")
#' @rdname code-printers
#' @export
fcode <- .mk_printer("fortran", "%s = %s")
#' @rdname code-printers
#' @export
jscode <- .mk_printer("js", "%s = %s;")
#' @rdname code-printers
#' @export
octave_code <- .mk_printer("octave", "%s = %s;")
#' @rdname code-printers
#' @export
rcode <- .mk_printer("r", "%s <- %s")
#' @rdname code-printers
#' @export
pycode <- .mk_printer("python", "%s = %s")

# ---------------------------------------------------------------------------
# lambdify and cse

#' Turn an expression into an R function
#'
#' `lambdify(args, e)` returns a fast vectorized R function, like SymPy's
#' `lambdify()`. The function body is built as an R call from the expression
#' tree.
#'
#' @param args A symbol or list of symbols.
#' @param e An expression or a list of expressions.
#' @return An R function.
#' @examples
#' x <- Symbol("x")
#' f <- lambdify(x, x^2 + sin(x))
#' f(c(0, 1, 2))
#' @export
lambdify <- function(args, e) {
  args <- .as_gens(list(args))
  nms <- vapply(args, function(s) s$name, "")
  .cg_imp_namespace(e)
  body <- if (inherits(e, "symmatrix"))
            as.call(c(list(as.name("matrix"), as.call(c(list(as.name("c")), lapply(e$data, .to_call)))),
                      list(nrow = e$nrow, ncol = e$ncol, byrow = TRUE)))
          else if (is.list(e) && !inherits(e, "sym")) as.call(c(list(as.name("list")), lapply(e, function(z) .to_call(as_sym(z)))))
          else .to_call(as_sym(e))
  fn <- function() NULL
  formals(fn) <- stats::setNames(rep(list(quote(expr = )), length(nms)), nms)
  body(fn) <- body
  environment(fn) <- list2env(list(), parent = asNamespace("stats"))
  fn
}

.R_FUN_MAP <- c(Abs = "abs", ceiling = "ceiling", floor = "floor", loggamma = "lgamma")

.to_call <- function(e) {
  t <- e$type
  switch(t,
         Integer = , Rational = as.double(e$value),
         Float = e$value,
         Symbol = , Dummy = as.name(e$name),
         Constant = switch(EXPR = e$name, pi = quote(pi), E = exp(1), I = 1i, oo = Inf, `-oo` = -Inf, nan = NaN,
                           .CONST_VALUES[[e$name]]),
         Add = Reduce(function(a, b) call("+", a, b), lapply(e$args, .to_call)),
         Mul = Reduce(function(a, b) call("*", a, b), lapply(e$args, .to_call)),
         Pow = {
           if (identical(e$args[[2L]]$key, S_Half$key)) return(call("sqrt", .to_call(e$args[[1L]])))
           call("^", .to_call(e$args[[1L]]), .to_call(e$args[[2L]]))
         },
         Function = {
           nm <- e$name
           imp <- e$extra$imp
           if (is.function(imp)) return(as.call(c(list(imp), lapply(e$args, .to_call))))
           if (!is.null(imp)) return(.to_call(.cg_apply_lambda(imp, e$args)))
           mapped <- .R_FUN_MAP[nm]
           if (!is.na(mapped)) return(as.call(c(list(as.name(mapped)), lapply(e$args, .to_call))))
           if (exists(nm, envir = baseenv(), mode = "function") && nm %in%
               c("sin", "cos", "tan", "asin", "acos", "atan", "sinh", "cosh", "tanh", "asinh", "acosh", "atanh",
                 "exp", "log", "gamma", "beta", "sign", "sqrt", "atan2", "factorial"))
             return(as.call(c(list(as.name(nm)), lapply(e$args, .to_call))))
           d <- .FUNS[[nm]]
           if (is.null(d) || is.null(d$numeric)) stop("cannot lambdify ", nm, "()")
           as.call(c(list(function(...) mapply(d$numeric, ...)), lapply(e$args, .to_call)))
         },
         Piecewise = {
           out <- NaN
           for (p in rev(.pw_pairs(e))) out <- call("ifelse", .to_call(p[[2L]]), .to_call(p[[1L]]), out)
           out
         },
         BooleanAtom = e$name == "true",
         Equality = call("==", .to_call(e$args[[1L]]), .to_call(e$args[[2L]])),
         Unequality = call("!=", .to_call(e$args[[1L]]), .to_call(e$args[[2L]])),
         StrictLessThan = call("<", .to_call(e$args[[1L]]), .to_call(e$args[[2L]])),
         LessThan = call("<=", .to_call(e$args[[1L]]), .to_call(e$args[[2L]])),
         StrictGreaterThan = call(">", .to_call(e$args[[1L]]), .to_call(e$args[[2L]])),
         GreaterThan = call(">=", .to_call(e$args[[1L]]), .to_call(e$args[[2L]])),
         And = Reduce(function(a, b) call("&", a, b), lapply(e$args, .to_call)),
         Or = Reduce(function(a, b) call("|", a, b), lapply(e$args, .to_call)),
         Not = call("!", .to_call(e$args[[1L]])),
         stop("cannot lambdify ", t))
}

#' Common subexpression elimination
#'
#' `cse(exprs)` finds common subexpressions as SymPy's `cse()`: repeated
#' subtrees, negated terms, negative powers, and the common arguments of
#' sums and products. It returns a list with `replacements` (a list of
#' `list(symbol, expression)` pairs) and `reduced` (the rewritten
#' expressions).
#'
#' @param exprs An expression, a matrix, or a list of them.
#' @param prefix Prefix of the new symbols when `symbols` is `NULL`.
#' @param symbols A list of symbols to use for the replacements, such as
#'   `numbered_symbols("t")`, or `NULL` for `prefix` followed by a number.
#' @param optimizations A list of `list(pre, post)` pairs of functions
#'   applied to each expression before and after the elimination, or `NULL`.
#' @param postprocess A function of the replacements and the reduced
#'   expressions whose value is returned, or `NULL`.
#' @param order `"canonical"` to visit the arguments of sums and products
#'   in sorted order, or `"none"`.
#' @param ignore Symbols; subexpressions containing any of them are not
#'   replaced.
#' @param list If `FALSE`, a single expression gives a single reduced
#'   expression rather than a list.
#' @return A list.
#' @examples
#' x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); t <- Symbol("t")
#' cse(list((x + y)^2 + sin(x + y), (x + y) * 3))
#' cse(list(x*y*z, x*y*t))
#' @export
cse <- function(exprs, prefix = "x", symbols = NULL, optimizations = NULL, postprocess = NULL,
                order = "canonical", ignore = list(), list = TRUE)
  .fs_cse(exprs, symbols, optimizations, postprocess, order, ignore, list, prefix)
