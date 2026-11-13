# Further printers, after sympy/printing: rust.py, glsl.py, python.py,
# smtlib.py, latex.multiline_latex, mathml.py and dot.py.

.negz <- function(z) if (z$type == "Mul") .neg_mul_noflat(z) else .neg(z)

.rust_num <- function(q) paste0(gmp::numerator(q), "_f64/", gmp::denominator(q), ".0")
.RUST_FUNS <- c(sin = "sin", cos = "cos", tan = "tan", asin = "asin", acos = "acos", atan = "atan",
                sinh = "sinh", cosh = "cosh", tanh = "tanh", asinh = "asinh", acosh = "acosh", atanh = "atanh",
                exp = "exp", log = "ln", Abs = "abs", floor = "floor", ceiling = "ceil", sign = "signum")

.rust <- function(e) {
  recv <- function(z) {
    if (z$type == "Integer") return(paste0(as.character(gmp::numerator(z$value)), ".0_f64"))
    s <- .rust(z)
    if (.prec(z) < PREC$Atom) paste0("(", s, ")") else s
  }
  mulp <- function(z) { s <- .rust(z); if (.prec(z) < PREC$Mul) paste0("(", s, ")") else s }
  switch(e$type,
    Integer = as.character(gmp::numerator(e$value)),
    Rational = .rust_num(e$value),
    Float = format(e$value, digits = 17),
    Symbol = , Dummy = e$name,
    Constant = switch(EXPR = e$name, pi = "PI", E = "E", oo = "INFINITY", `-oo` = "NEG_INFINITY", nan = "NAN",
                      stop("cannot print constant ", e$name, " in Rust", call. = FALSE)),
    Add = {
      terms <- .add_print_order(e$args)
      out <- vapply(seq_along(terms), function(i) {
        z <- terms[[i]]
        if (i > 1L && could_extract_minus_sign(z)) paste0(" - ", mulp(.negz(z)))
        else paste0(if (i > 1L) " + ", .rust(z))
      }, "")
      paste(out, collapse = "")
    },
    Mul = {
      if (could_extract_minus_sign(e)) return(paste0("-", mulp(.negz(e))))
      c1 <- e$args[[1L]]
      if (c1$type == "Rational") return(paste0("(", .rust_num(c1$value), ")*", mulp(.mul(e$args[-1L]))))
      paste(vapply(.mul_print_order(e$args), mulp, ""), collapse = "*")
    },
    Pow = {
      b <- e$args[[1L]]; x <- e$args[[2L]]
      if (!length(free_symbols(e)) && !is_number_type(e)) return(format(signif(Re(as.complex(.neval(e, list()))), 15), digits = 15))
      if (identical(x$key, S_Half$key)) return(paste0(recv(b), ".sqrt()"))
      if (identical(x$key, S_NegativeOne$key)) return(paste0(recv(b), ".recip()"))
      if (x$type == "Integer") return(paste0(recv(b), ".powi(", .rust(x), ")"))
      if (is_const(b, "E")) return(paste0(recv(x), ".exp()"))
      paste0(recv(b), ".powf(", if (x$type == "Rational") .rust_num(x$value) else .rust(x), ")")
    },
    Function = {
      m <- .RUST_FUNS[e$name]
      if (is.na(m) || length(e$args) != 1L) stop("cannot print ", e$name, " in Rust", call. = FALSE)
      paste0(recv(e$args[[1L]]), ".", m, "()")
    },
    stop("cannot generate Rust code for ", e$type, call. = FALSE))
}

.glsl_num <- function(q) paste0(gmp::numerator(q), ".0/", gmp::denominator(q), ".0")
.glsl <- function(e, consts) {
  rec <- function(z) .glsl(z, consts)
  mulp <- function(z) { s <- rec(z); if (.prec(z) < PREC$Mul) paste0("(", s, ")") else s }
  flt <- function(z) {
    if (z$type == "Integer") return(paste0(as.character(gmp::numerator(z$value)), ".0"))
    if (z$type == "Rational") return(format(as.numeric(z$value), digits = 15))
    rec(z)
  }
  switch(e$type,
    Integer = as.character(gmp::numerator(e$value)),
    Rational = .glsl_num(e$value),
    Float = format(e$value, digits = 15),
    Symbol = , Dummy = e$name,
    Constant = {
      v <- switch(EXPR = e$name, pi = list("pi", "3.14159265"), E = list("E", "2.71828183"), NULL)
      if (is.null(v)) stop("cannot print constant ", e$name, " in GLSL", call. = FALSE)
      consts[[v[[1L]]]] <- v[[2L]]
      v[[1L]]
    },
    Add = {
      terms <- .add_print_order(e$args)
      out <- vapply(seq_along(terms), function(i) {
        z <- terms[[i]]
        if (i > 1L && could_extract_minus_sign(z)) paste0(" - ", mulp(.negz(z)))
        else paste0(if (i > 1L) " + ", rec(z))
      }, "")
      paste(out, collapse = "")
    },
    Mul = {
      if (could_extract_minus_sign(e)) return(paste0("-", mulp(.negz(e))))
      c1 <- e$args[[1L]]
      if (c1$type == "Rational") return(paste0("(", .glsl_num(c1$value), ")*", mulp(.mul(e$args[-1L]))))
      paste(vapply(.mul_print_order(e$args), mulp, ""), collapse = "*")
    },
    Pow = {
      b <- e$args[[1L]]; x <- e$args[[2L]]
      if (identical(x$key, S_Half$key)) return(paste0("sqrt(", rec(b), ")"))
      if (identical(x$key, S_NegativeOne$key)) return(paste0("1.0/", mulp(b)))
      paste0("pow(", rec(b), ", ", flt(x), ")")
    },
    Function = paste0(switch(e$name, Abs = "abs", ceiling = "ceil", log = "log", atan2 = "atan", e$name), "(",
                      paste(vapply(e$args, rec, ""), collapse = ", "), ")"),
    stop("cannot generate GLSL code for ", e$type, call. = FALSE))
}

#' More code printers
#'
#' `rust_code()` prints Rust (method call syntax such as `x.powi(2)`),
#' `glsl_code()` GLSL, `smtlib_code()` SMT-LIB s-expressions (relations
#' become `assert`s, with `declare-const` lines for the symbols when
#' `auto_declare` is `TRUE`) and `python()` a Python snippet that rebuilds
#' the expression with SymPy. `multiline_latex(lhs, rhs, terms_per_line)`
#' breaks a long right side over lines of an `align*` environment.
#' `print_glsl()`, `print_python()` and `print_maple_code()` print the
#' strings.
#'
#' @param e An expression, or for `smtlib_code()` a list of them.
#' @param assign_to Optional variable name.
#' @param ... Printer settings of `rust_code()`, as in SymPy (for example
#'   `user_functions` or `contract`).
#' @param auto_declare Declare the symbols as reals.
#' @param lhs,rhs Sides of the equation.
#' @param terms_per_line Terms per line.
#' @param environment The LaTeX environment.
#' @return A character scalar; the `print_` functions print it.
#' @examples
#' x <- Symbol("x")
#' rust_code(x^2 + sin(x)/3)
#' glsl_code(sqrt(x + 1))
#' smtlib_code(Eq(x + 1, 3))
#' @name more-code-printers
NULL

#' @rdname more-code-printers
#' @export
rust_code <- function(e, assign_to = NULL, ...) .cg_rust_code(e, assign_to, ...)

#' @rdname more-code-printers
#' @export
glsl_code <- function(e, assign_to = NULL) {
  consts <- new.env()
  s <- .glsl(as_sym(e), consts)
  if (!is.null(assign_to)) s <- paste0(assign_to, " = ", s, ";")
  decl <- vapply(sort(ls(consts)), function(k) paste0("float ", k, " = ", consts[[k]], ";\n"), "")
  paste0(paste(decl, collapse = ""), s)
}

#' @rdname more-code-printers
#' @export
print_glsl <- function(e, assign_to = NULL) {
  cat(glsl_code(e, assign_to), "\n", sep = "")
  invisible(NULL)
}

#' @rdname more-code-printers
#' @export
print_maple_code <- function(e, assign_to = NULL) {
  cat(maple_code(e, assign_to), "\n", sep = "")
  invisible(NULL)
}

#' @rdname more-code-printers
#' @export
python <- function(e) {
  e <- as_sym(e)
  old <- .str_mode$python
  .str_mode$python <- TRUE
  on.exit(.str_mode$python <- old)
  # Python keywords get a trailing underscore as identifiers, as SymPy's python()
  kw <- c("False", "None", "True", "and", "as", "assert", "async", "await", "break", "class", "continue", "def", "del",
          "elif", "else", "except", "finally", "for", "from", "global", "if", "import", "in", "is", "lambda", "nonlocal",
          "not", "or", "pass", "raise", "return", "try", "while", "with", "yield")
  orig <- atoms(e, "Symbol")
  s0 <- .str(e)
  fnames <- unique(vapply(Filter(.is_undefined_application, atoms(e, "Function")), function(f) f$name, ""))
  # a keyword gets underscores until the identifier collides with no other name
  used <- c(vapply(orig, function(v) v$name, ""), fnames)
  idmap <- character()
  # braces of subscripts are removed; SymPy substitutes the new name as a
  # string, so it is parsed (a^{2} becomes a**2 in the expression)
  repl <- list()
  for (nm in unique(vapply(orig, function(v) v$name, ""))) if (grepl("{", nm, fixed = TRUE)) {
    idmap[nm] <- gsub("[{}]", "", nm)
    repl[[nm]] <- tryCatch(sympify(idmap[[nm]]), error = function(err) Symbol(idmap[[nm]]))
  }
  for (nm in unique(used)) {
    base <- if (nm %in% names(idmap)) idmap[[nm]] else nm
    if (!(base %in% kw)) next
    new <- paste0(base, "_")
    while (new %in% c(used, kw, idmap)) new <- paste0(new, "_")
    idmap[nm] <- new
    repl[[nm]] <- Symbol(new)
  }
  ren <- Filter(function(v) v$name %in% names(repl), orig)
  if (length(ren)) e <- xreplace(e, lapply(ren, function(v) list(v, repl[[v$name]])))
  if (any(fnames %in% names(idmap)))
    e <- .map_bottom_up(e, function(z) if (.is_undefined_application(z) && z$name %in% names(idmap)) fnode(idmap[[z$name]], z$args) else z)
  shown <- function(nm) if (nm %in% names(idmap)) idmap[[nm]] else nm
  s <- .str(e)
  # every symbol, bound ones included, in order of appearance
  syms <- orig
  # declared in their order in the original printed expression, as SymPy
  pos <- vapply(syms, function(v) { m <- regexpr(paste0("(?<![A-Za-z0-9_])", gsub("([][{}()+*^$|\\\\.?])", "\\\\\\1", v$name), "(?![A-Za-z0-9_])"), s0, perl = TRUE); if (m < 0) .Machine$integer.max else m }, 0L)
  syms <- syms[order(pos)]
  funs <- fnames
  lines <- c(vapply(syms, function(v) paste0(shown(v$name), " = Symbol('", v$name, "')"), ""),
             vapply(funs, function(f) paste0(shown(f), " = Function('", f, "')"), ""), paste0("e = ", s))
  paste(lines, collapse = "\n")
}

#' @rdname more-code-printers
#' @export
print_python <- function(e) {
  cat(python(e), "\n", sep = "")
  invisible(NULL)
}

.smt_order <- function(args) {
  num <- vapply(args, is_number_type, TRUE)
  c(args[num], .key_sort(lapply(args[!num], .sort_key), args[!num]))
}
.smt <- function(e) {
  op <- function(name, args) paste0("(", name, " ", paste(vapply(args, .smt, ""), collapse = " "), ")")
  switch(e$type,
    Integer = as.character(gmp::numerator(e$value)),
    Rational = paste0("(/ ", gmp::numerator(e$value), " ", gmp::denominator(e$value), ")"),
    Float = format(e$value, digits = 17),
    Symbol = , Dummy = e$name,
    BooleanAtom = if (.is_true(e)) "true" else "false",
    Add = op("+", .smt_order(e$args)),
    Mul = op("*", .smt_order(e$args)),
    Pow = op("pow", e$args),
    Equality = op("=", e$args),
    Unequality = paste0("(not ", op("=", e$args), ")"),
    StrictLessThan = op("<", e$args), LessThan = op("<=", e$args),
    StrictGreaterThan = op(">", e$args), GreaterThan = op(">=", e$args),
    And = op("and", e$args), Or = op("or", e$args), Not = op("not", e$args),
    Function = op(e$name, e$args),
    Constant = switch(EXPR = e$name, pi = "pi", E = "(exp 1)", stop("cannot print constant ", e$name, " in SMT-LIB", call. = FALSE)),
    stop("cannot print ", e$type, " in SMT-LIB", call. = FALSE))
}

#' @rdname more-code-printers
#' @export
smtlib_code <- function(e, auto_declare = TRUE) {
  es <- if (is.list(e) && !is_sym(e)) lapply(e, as_sym) else list(as_sym(e))
  # Bool for symbols in boolean position (including an equality with a
  # boolean operand), Int for integer symbols
  boolk <- character()
  isbool <- function(z) z$type %in% c("BooleanAtom", "And", "Or", "Not", "Xor", "Implies", "Equivalent") || z$key %in% boolk
  walk <- function(z, bctx) {
    if (is_symbol_like(z) && bctx) boolk <<- c(boolk, z$key)
    b <- z$type %in% c("And", "Or", "Not", "Xor", "Implies", "Equivalent") ||
      (z$type %in% c("Equality", "Unequality") && any(vapply(z$args, isbool, TRUE)))
    for (a in z$args) walk(a, b)
  }
  # repeat until no new boolean symbol appears, so the order of es is irrelevant
  repeat {
    n0 <- length(unique(boolk))
    for (z in es) walk(z, TRUE)
    boolk <- unique(boolk)
    if (length(boolk) == n0) break
  }
  body <- vapply(es, function(z) {
    s <- .smt(z)
    if (is_relational(z) || z$type %in% c("And", "Or", "Not") || (is_symbol_like(z) && z$key %in% boolk)) paste0("(assert ", s, ")") else s
  }, "")
  if (!auto_declare) return(paste(body, collapse = "\n"))
  syms <- unique_syms(lapply(es, free_symbols))
  sort_of <- function(v) if (v$key %in% boolk) "Bool" else if (isTRUE(.facts(v)[["integer"]])) "Int" else "Real"
  decl <- vapply(syms, function(v) paste0("(declare-const ", v$name, " ", sort_of(v), ")"), "")
  paste(c(decl, body), collapse = "\n")
}

#' @rdname more-code-printers
#' @export
multiline_latex <- function(lhs, rhs, terms_per_line = 1, environment = "align*") {
  lhs <- as_sym(lhs); rhs <- as_sym(rhs)
  terms <- if (rhs$type == "Add") .add_print_order(rhs$args) else list(rhs)
  pieces <- vapply(seq_along(terms), function(i) {
    z <- terms[[i]]
    if (i == 1L) return(.latex(z))
    if (could_extract_minus_sign(z)) paste0("- ", .latex(if (z$type == "Mul") .negz(z) else .neg(z)))
    else paste0("+ ", .latex(z))
  }, "")
  lines <- character()
  for (st in seq(1L, length(pieces), by = terms_per_line)) {
    chunk <- pieces[st:min(length(pieces), st + terms_per_line - 1L)]
    lines <- c(lines, paste(chunk, collapse = " "))
  }
  body <- paste0(.latex(lhs), " = & ", lines[1L])
  if (length(lines) > 1L) body <- paste0(body, " \\\\\n", paste0("& ", lines[-1L], collapse = " \\\\\n"))
  paste0("\\begin{", environment, "}\n", body, " \n\\end{", environment, "}")
}

# ---------------------------------------------------------------------------
# MathML

.mml_content <- function(e) {
  ap <- function(op, args) paste0("<apply><", op, "/>", paste(vapply(args, .mml_content, ""), collapse = ""), "</apply>")
  switch(e$type,
    Integer = paste0("<cn>", gmp::numerator(e$value), "</cn>"),
    Rational = paste0("<apply><divide/><cn>", gmp::numerator(e$value), "</cn><cn>", gmp::denominator(e$value), "</cn></apply>"),
    Float = paste0("<cn>", format(e$value, digits = 15), "</cn>"),
    Symbol = , Dummy = paste0("<ci>", e$name, "</ci>"),
    Constant = switch(EXPR = e$name, pi = "<pi/>", E = "<exponentiale/>", I = "<imaginaryi/>", oo = "<infinity/>",
                      `-oo` = "<apply><minus/><infinity/></apply>", nan = "<notanumber/>", paste0("<ci>", e$name, "</ci>")),
    Add = {
      terms <- .add_print_order(e$args)
      paste0("<apply><plus/>", paste(vapply(terms, function(z) {
        if (could_extract_minus_sign(z)) paste0("<apply><minus/>", .mml_content(.negz(z)), "</apply>") else .mml_content(z)
      }, ""), collapse = ""), "</apply>")
    },
    Mul = {
      if (could_extract_minus_sign(e)) return(paste0("<apply><minus/>", .mml_content(.negz(e)), "</apply>"))
      nd <- fraction(e, exact = TRUE)
      if (!is_one_num(nd[[2L]])) return(ap("divide", list(nd[[1L]], nd[[2L]])))
      ap("times", .mul_print_order(e$args))
    },
    Pow = {
      if (identical(e$args[[2L]]$key, S_Half$key)) return(ap("root", list(e$args[[1L]])))
      if (is_number_type(e$args[[2L]]) && num_sign(e$args[[2L]]) < 0)
        return(ap("divide", list(S_One, .pow(e$args[[1L]], .neg(e$args[[2L]])))))
      ap("power", e$args)
    },
    Function = paste0("<apply><", switch(e$name, Abs = "abs", log = "ln", e$name), "/>",
                      paste(vapply(e$args, .mml_content, ""), collapse = ""), "</apply>"),
    Equality = ap("eq", e$args), Unequality = ap("neq", e$args), StrictLessThan = ap("lt", e$args),
    LessThan = ap("leq", e$args), StrictGreaterThan = ap("gt", e$args), GreaterThan = ap("geq", e$args),
    stop("cannot print ", e$type, " in MathML", call. = FALSE))
}

.mml_pres <- function(e) {
  row <- function(...) paste0("<mrow>", ..., "</mrow>")
  mi <- function(s) paste0("<mi>", s, "</mi>")
  mn <- function(s) paste0("<mn>", s, "</mn>")
  mo <- function(s) paste0("<mo>", s, "</mo>")
  par <- function(z) if (z$type %in% c("Add") || (z$type == "Mul" && could_extract_minus_sign(z)))
    row(mo("("), .mml_pres(z), mo(")")) else .mml_pres(z)
  switch(e$type,
    Integer = mn(gmp::numerator(e$value)),
    Rational = paste0("<mfrac>", mn(gmp::numerator(e$value)), mn(gmp::denominator(e$value)), "</mfrac>"),
    Float = mn(format(e$value, digits = 15)),
    Symbol = , Dummy = mi(e$name),
    Constant = switch(EXPR = e$name, pi = "<mi>&pi;</mi>", E = "<mi>&ExponentialE;</mi>", I = "<mi>&ImaginaryI;</mi>",
                      oo = "<mi>&infin;</mi>", mi(e$name)),
    Add = {
      terms <- .add_print_order(e$args)
      parts <- vapply(seq_along(terms), function(i) {
        z <- terms[[i]]
        if (could_extract_minus_sign(z)) {
          inner <- .mml_pres(.negz(z))
          if (i == 1L) row(mo("-"), inner) else paste0(mo("-"), inner)
        } else paste0(if (i > 1L) mo("+"), .mml_pres(z))
      }, "")
      row(paste(parts, collapse = ""))
    },
    Mul = {
      if (could_extract_minus_sign(e)) return(row(mo("-"), .mml_pres(.negz(e))))
      nd <- fraction(e, exact = TRUE)
      if (!is_one_num(nd[[2L]])) return(row(paste0("<mfrac>", .mml_pres(nd[[1L]]), .mml_pres(nd[[2L]]), "</mfrac>")))
      row(paste(vapply(.mul_print_order(e$args), par, ""), collapse = mo("&InvisibleTimes;")))
    },
    Pow = {
      b <- e$args[[1L]]; x <- e$args[[2L]]
      if (identical(x$key, S_Half$key)) return(paste0("<msqrt>", .mml_pres(b), "</msqrt>"))
      if (is_number_type(x) && num_sign(x) < 0) return(row(paste0("<mfrac>", mn(1), .mml_pres(.pow(b, .neg(x))), "</mfrac>")))
      paste0("<msup>", par(b), .mml_pres(x), "</msup>")
    },
    Function = row(mi(e$name), row(mo("("), paste(vapply(e$args, .mml_pres, ""), collapse = mo(",")), mo(")"))),
    stop("cannot print ", e$type, " in MathML", call. = FALSE))
}

#' MathML and Graphviz output
#'
#' `mathml(e, printer)` writes content MathML (default) or presentation
#' MathML (`printer = "presentation"`); `print_mathml()` prints it.
#' `dotprint(e)` writes the expression tree in the Graphviz dot language,
#' with SymPy's node naming.
#'
#' @param e An expression.
#' @param printer `"content"` or `"presentation"`.
#' @return A character scalar.
#' @examples
#' x <- Symbol("x")
#' mathml(x^2 + 1)
#' mathml(x^2 + 1, printer = "presentation")
#' cat(dotprint(x + 2))
#' @export
mathml <- function(e, printer = "content") {
  e <- as_sym(e)
  if (identical(printer, "presentation")) {
    s <- .mml_pres(e)
    if (!startsWith(s, "<mrow>") && !startsWith(s, "<msqrt>")) s <- paste0("<mrow>", s, "</mrow>")
    s
  } else .mml_content(e)
}

#' @rdname mathml
#' @export
print_mathml <- function(e, printer = "content") {
  cat(mathml(e, printer), "\n", sep = "")
  invisible(NULL)
}

#' @rdname mathml
#' @export
dotprint <- function(e) {
  e <- as_sym(e)
  nodes <- character(); edges <- character()
  label <- function(z) if (length(z$args) && !(z$type %in% c("Integer", "Rational", "Float"))) {
    if (z$type == "Function") z$name else z$type
  } else .str(z)
  walk <- function(z, pos) {
    id <- paste0("\"", srepr(z), "_(", paste0(pos, if (length(pos) == 1L) ",", collapse = ", "), ")\"")
    if (!length(pos)) id <- paste0("\"", srepr(z), "_()\"")
    nodes <<- c(nodes, paste0(id, " [\"color\"=\"black\", \"label\"=\"", label(z), "\", \"shape\"=\"ellipse\"];"))
    if (length(z$args) && !(z$type %in% c("Integer", "Rational", "Float"))) for (i in seq_along(z$args)) {
      cid <- walk(z$args[[i]], c(pos, i - 1L))
      edges <<- c(edges, paste0(id, " -> ", cid, ";"))
    }
    id
  }
  walk(e, integer())
  paste0("digraph{\n\n# Graph style\n\"ordering\"=\"out\"\n\"rankdir\"=\"TD\"\n\n#########\n# Nodes #\n#########\n\n",
         paste(nodes, collapse = "\n"), "\n\n#########\n# Edges #\n#########\n\n", paste(edges, collapse = "\n"), "\n}")
}
