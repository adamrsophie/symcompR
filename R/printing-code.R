# LaTeX and code printers, lambdify and common subexpression elimination,
# after sympy/printing and sympy/simplify/cse_main.py.

.GREEK <- c("alpha", "beta", "gamma", "delta", "epsilon", "zeta", "eta", "theta", "iota",
            "kappa", "lambda", "mu", "nu", "xi", "pi", "rho", "sigma", "tau", "upsilon",
            "phi", "chi", "psi", "omega", "Gamma", "Delta", "Theta", "Lambda", "Xi", "Pi",
            "Sigma", "Upsilon", "Phi", "Psi", "Omega")

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

.latex_symbol <- function(nm) {
  if (grepl("{", nm, fixed = TRUE)) return(nm)
  tr <- function(s) if (s %in% .GREEK) paste0("\\", s) else s
  p <- .split_super_sub(nm)
  out <- tr(p$name)
  if (length(p$supers)) out <- paste0(out, "^{", paste(vapply(p$supers, tr, ""), collapse = " "), "}")
  if (length(p$subs)) out <- paste0(out, "_{", paste(vapply(p$subs, tr, ""), collapse = " "), "}")
  out
}

.LATEX_FUNS <- c(sin = "\\sin", cos = "\\cos", tan = "\\tan", cot = "\\cot", sec = "\\sec", csc = "\\csc",
                 asin = "\\operatorname{asin}", acos = "\\operatorname{acos}", atan = "\\operatorname{atan}",
                 acot = "\\operatorname{acot}", sinh = "\\sinh", cosh = "\\cosh", tanh = "\\tanh", coth = "\\coth",
                 asinh = "\\operatorname{asinh}", acosh = "\\operatorname{acosh}", atanh = "\\operatorname{atanh}",
                 erf = "\\operatorname{erf}", erfc = "\\operatorname{erfc}", zeta = "\\zeta")

.latex_paren <- function(e, level, strict = FALSE) {
  s <- .latex(e)
  p <- .prec(e)
  if (p < level || (!strict && p <= level)) paste0("\\left(", s, "\\right)") else s
}

.latex_mul <- function(e) {
  # a kept coefficient times a sum, as collect_const() builds: 2*(a + y)
  if (length(e$args) == 2L && is_number_type(e$args[[1L]]) && e$args[[2L]]$type == "Add") {
    c1 <- e$args[[1L]]
    inner <- paste0("\\left(", .latex(e$args[[2L]]), "\\right)")
    if (num_sign(c1) < 0) return(if (is_one_num(.neg(c1))) paste0("- ", inner) else paste0("- ", .latex(.neg(c1)), " ", inner))
    return(paste(.latex(c1), inner))
  }
  if (could_extract_minus_sign(e)) {
    ne <- if (e$type == "Mul") .neg_mul_noflat(e) else .neg(e)
    s <- .latex(ne)
    if (.prec(ne) < PREC$Mul) s <- paste0("\\left(", s, "\\right)")
    return(paste0("- ", s))
  }
  nd <- fraction(e, exact = TRUE)
  fmt <- function(z) {
    fs <- if (z$type == "Mul") .mul_print_order(z$args) else list(z)
    paste(vapply(fs, function(f) .latex_paren(f, PREC$Mul, strict = TRUE), ""), collapse = " ")
  }
  if (is_one_num(nd[[2L]])) return(fmt(nd[[1L]]))
  paste0("\\frac{", .latex(nd[[1L]]), "}{", .latex(nd[[2L]]), "}")
}

.latex_pow <- function(e) {
  b <- e$args[[1L]]
  x <- e$args[[2L]]
  if (x$type == "Rational" && gmp::numerator(x$value) == 1) {
    q <- gmp::denominator(x$value)
    if (q == 2) return(paste0("\\sqrt{", .latex(b), "}"))
    return(paste0("\\sqrt[", as.character(q), "]{", .latex(b), "}"))
  }
  # a noncommutative base keeps its negative exponent, A^{-1}
  if (.is_neg_exp(x) && .is_comm(b)) return(paste0("\\frac{1}{", .latex(.pow(b, .neg(x))), "}"))
  bs <- if (b$type == "Function") paste0("\\left(", .latex(b), "\\right)") else .latex_paren(b, PREC$Pow)
  paste0(bs, "^{", .latex(x), "}")
}

.latex_add <- function(e) {
  terms <- .add_print_order(e$args)
  out <- character(length(terms))
  for (i in seq_along(terms)) {
    t <- terms[[i]]
    if (i == 1L) out[i] <- .latex(t)
    else if (could_extract_minus_sign(t)) out[i] <- paste0(" - ", .latex(if (t$type == "Mul") .neg_mul_noflat(t) else .neg(t)))
    else out[i] <- paste0(" + ", .latex(t))
  }
  paste(out, collapse = "")
}

.latex <- function(e) {
  if (!is.null(.latex_mode$hook) && !is.null(r <- .latex_mode$hook(e))) return(r)
  t <- e$type
  switch(t,
         Integer = as.character(gmp::numerator(e$value)),
         Rational = {
           n <- gmp::numerator(e$value); d <- gmp::denominator(e$value)
           if (n < 0) paste0("- \\frac{", as.character(-n), "}{", as.character(d), "}")
           else paste0("\\frac{", as.character(n), "}{", as.character(d), "}")
         },
         Float = .fmt_float_sym(e, strip = TRUE),
         Wild = if (isTRUE(e$extra$wildfunction)) paste0("\\operatorname{WildFunction}{\\left(", .latex_symbol(e$name), " \\right)}")
                else .latex_symbol(e$name),
         Symbol = , Dummy = .latex_symbol(e$name),
         Constant = switch(EXPR = e$name, pi = "\\pi", E = "e", I = "i", oo = "\\infty", `-oo` = "-\\infty",
                           zoo = "\\tilde{\\infty}", nan = "\\text{NaN}", EulerGamma = "\\gamma",
                           GoldenRatio = "\\phi", Catalan = "G", e$name),
         BooleanAtom = if (e$name == "true") "\\text{True}" else "\\text{False}",
         Add = .latex_add(e),
         Mul = .latex_mul(e),
         Pow = .latex_pow(e),
         Function = {
           d <- .FUNS[[e$name]]
           args <- vapply(e$args, .latex, "")
           if (!is.null(d) && !is.null(d$latex)) d$latex(args)
           else if (!is.na(.LATEX_FUNS[e$name])) paste0(.LATEX_FUNS[[e$name]], "{\\left(", paste(args, collapse = ", "), " \\right)}")
           else if (is.null(d)) paste0(.latex_symbol(e$name), "{\\left(", paste(args, collapse = ", "), " \\right)}")
           else paste0("\\operatorname{", e$name, "}{\\left(", paste(args, collapse = ", "), "\\right)}")
         },
         Equality = paste(.latex(e$args[[1L]]), "=", .latex(e$args[[2L]])),
         Unequality = paste(.latex(e$args[[1L]]), "\\neq", .latex(e$args[[2L]])),
         StrictLessThan = paste(.latex(e$args[[1L]]), "<", .latex(e$args[[2L]])),
         LessThan = paste(.latex(e$args[[1L]]), "\\leq", .latex(e$args[[2L]])),
         StrictGreaterThan = paste(.latex(e$args[[1L]]), ">", .latex(e$args[[2L]])),
         GreaterThan = paste(.latex(e$args[[1L]]), "\\geq", .latex(e$args[[2L]])),
         And = paste(vapply(e$args, .latex, ""), collapse = " \\wedge "),
         Or = paste(vapply(e$args, .latex, ""), collapse = " \\vee "),
         Not = paste0("\\neg ", .latex(e$args[[1L]])),
         Derivative = {
           f <- e$args[[1L]]
           tot <- sum(vapply(e$args[-1L], function(tp) as.numeric(tp$args[[2L]]), 0))
           dn <- paste(vapply(e$args[-1L], function(tp) {
             n <- as.numeric(tp$args[[2L]])
             paste0("\\partial ", .latex(tp$args[[1L]]), if (n > 1) paste0("^{", n, "}"))
           }, ""), collapse = " ")
           d <- if (length(e$args) == 2L && !.is_undefined_multi(f)) "d" else "\\partial"
           dn <- gsub("\\partial", d, dn, fixed = TRUE)
           paste0("\\frac{", d, if (tot > 1) paste0("^{", tot, "}"), "}{", dn, "} ", .latex(f))
         },
         Integral = {
           s <- .latex(e$args[[1L]])
           for (lim in rev(e$args[-1L])) {
             v <- .latex(lim$args[[1L]])
             bounds <- if (length(lim$args) == 3L) paste0("_{", .latex(lim$args[[2L]]), "}^{", .latex(lim$args[[3L]]), "}") else ""
             s <- paste0(if (nzchar(bounds)) "\\int\\limits" else "\\int", bounds, " ", s, "\\, d", v)
           }
           s
         },
         Sum = , Product = {
           op <- if (t == "Sum") "\\sum" else "\\prod"
           lim <- e$args[[2L]]
           paste0(op, "_{", .latex(lim$args[[1L]]), "=", .latex(lim$args[[2L]]), "}^{", .latex(lim$args[[3L]]), "} ", .latex(e$args[[1L]]))
         },
         Limit = {
           # LatexPrinter._print_Limit: the direction unless two sided or at
           # infinity; an AssocOp body in parentheses
           z0 <- e$args[[3L]]
           dir <- e$args[[4L]]$name
           tex <- paste0("\\lim_{", .latex(e$args[[2L]]), " \\to ", .latex(z0),
                         if (dir != "+-" && !is_const(z0, "oo") && !is_const(z0, "-oo")) paste0("^", dir), "}")
           b <- e$args[[1L]]
           if (b$type %in% c("Add", "Mul", "And", "Or") || (b$type == "Function" && b$name %in% c("Max", "Min")))
             paste0(tex, "\\left(", .latex(b), "\\right)") else paste0(tex, " ", .latex(b))
         },
         Order = paste0("O\\left(", .latex(e$args[[1L]]), "\\right)"),
         Piecewise = {
           pr <- .pw_pairs(e)
           tx <- vapply(seq_along(pr), function(i) {
             p <- pr[[i]]
             if (i == length(pr) && .is_true(p[[2L]])) paste0(.latex(p[[1L]]), " & \\text{otherwise}")
             else paste0(.latex(p[[1L]]), " & \\text{for}\\: ", .latex(p[[2L]]))
           }, "")
           paste0("\\begin{cases} ", paste(tx, collapse = " \\\\"), " \\end{cases}")
         },
         {
           f <- .latex_printers[[t]]
           if (!is.null(f)) f(e) else paste0("\\text{", .str(e), "}")
         })
}

.is_undefined_multi <- function(f) length(free_symbols(f)) > 1L

#' LaTeX output
#'
#' @param e A `sym` expression or matrix.
#' @return A character scalar.
#' @examples
#' x <- Symbol("x")
#' latex(x^2 / 2 + sin(x))
#' @export
latex <- function(e) {
  if (inherits(e, "TableForm")) return(.tableform_latex(e))
  if (inherits(e, "symmatrix")) {
    rows <- vapply(.mrows(e), function(r) paste(vapply(r, .latex, ""), collapse = " & "), "")
    return(paste0("\\left[\\begin{matrix}", paste(rows, collapse = "\\\\"), "\\end{matrix}\\right]"))
  }
  .latex(as_sym(e))
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
  body <- if (is.list(e) && !inherits(e, "sym")) as.call(c(list(as.name("list")), lapply(e, function(z) .to_call(as_sym(z)))))
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
#' `cse(exprs)` finds repeated subexpressions and returns a list with
#' `replacements` (a list of `list(symbol, expression)` pairs) and `reduced`
#' (the rewritten expressions), like SymPy's `cse()`.
#'
#' @param exprs An expression or list of expressions.
#' @param prefix Prefix of the new symbols.
#' @return A list.
#' @examples
#' x <- Symbol("x"); y <- Symbol("y")
#' cse(list((x + y)^2 + sin(x + y), (x + y) * 3))
#' @export
cse <- function(exprs, prefix = "x") {
  if (!is.list(exprs) || inherits(exprs, "sym")) exprs <- list(exprs)
  exprs <- lapply(exprs, as_sym)
  counts <- new.env(hash = TRUE)
  nodes <- list()
  walk <- function(z) {
    if (!length(z$args) || nchar(z$key) > 9000L || z$type %in% c("MatrixSymbol", "MatrixElement")) return(invisible())
    cur <- counts[[z$key]]
    assign(z$key, if (is.null(cur)) 1L else cur + 1L, envir = counts)
    if (is.null(cur)) {
      nodes[[length(nodes) + 1L]] <<- z
      for (a in z$args) walk(a)
    }
  }
  for (e in exprs) walk(e)
  taken <- unique(unlist(lapply(exprs, function(e) vapply(atoms(e, c("Symbol", "Dummy")), function(s) s$name, ""))))
  rep_nodes <- Filter(function(z) counts[[z$key]] > 1L, nodes)
  # innermost first, so that later replacements use earlier symbols
  rep_nodes <- rep_nodes[order(vapply(rep_nodes, function(z) nchar(z$key), 0L))]
  replacements <- list()
  k <- 0L
  mapping <- list()
  apply_map <- function(z) {
    hit <- mapping[[z$key]]
    if (!is.null(hit)) return(hit)
    if (!length(z$args)) return(z)
    args <- lapply(z$args, apply_map)
    if (identical(vapply(args, .k, ""), vapply(z$args, .k, ""))) return(z)
    sym_rebuild(z, args)
  }
  for (z in rep_nodes) {
    repeat {
      nm <- paste0(prefix, k)
      k <- k + 1L
      if (!(nm %in% taken)) break
    }
    sy <- Symbol(nm)
    replacements[[length(replacements) + 1L]] <- list(sy, apply_map(z))
    mapping[[z$key]] <- sy
  }
  list(replacements = replacements, reduced = lapply(exprs, apply_map))
}
