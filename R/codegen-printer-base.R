# Code printer engine, after sympy/printing/printer.py, precedence.py, the
# parts of str.py that the code printers inherit, and codeprinter.py.
#
# A printer is an environment of class "cg_printer". Its class is a SymPy
# printer class name (such as "C99CodePrinter"); the classes form the same
# hierarchy as in SymPy, and methods are found along the printer MRO, while
# the expression classes are tried along the expression MRO, as in
# Printer._print.

.cgp_classes <- new.env()
.cgp_mro_cache <- new.env()

# Define a printer class. Methods are functions (p, e, ...); FALSE stands for
# a method set to None, which stops the lookup for that expression class.
.cgp_def <- function(name, parents = character(), attrs = list(), methods = list(), init = NULL) {
  assign(name, list(name = name, parents = parents, attrs = attrs,
                    methods = list2env(methods), init = init), envir = .cgp_classes)
  rm(list = ls(.cgp_mro_cache), envir = .cgp_mro_cache)
  invisible(NULL)
}
.cgp_set <- function(cls, ...) {
  m <- list(...)
  env <- .cgp_classes[[cls]]$methods
  for (k in names(m)) assign(k, m[[k]], envir = env)
  invisible(NULL)
}

# C3 linearization of the printer classes
.cgp_mro_of <- function(name) {
  hit <- .cgp_mro_cache[[name]]
  if (!is.null(hit)) return(hit)
  cl <- .cgp_classes[[name]]
  out <- name
  if (length(cl$parents)) {
    seqs <- c(lapply(cl$parents, .cgp_mro_of), list(cl$parents))
    repeat {
      seqs <- Filter(length, seqs)
      if (!length(seqs)) break
      cand <- NULL
      for (s in seqs) {
        h <- s[1L]
        if (!any(vapply(seqs, function(q) h %in% q[-1L], TRUE))) { cand <- h; break }
      }
      out <- c(out, cand)
      seqs <- lapply(seqs, function(q) if (q[1L] == cand) q[-1L] else q)
    }
  }
  assign(name, out, envir = .cgp_mro_cache)
  out
}

.cgp_new <- function(cls, settings = NULL) {
  p <- new.env()
  p$cls <- cls
  p$chain <- .cgp_mro_of(cls)
  p$print_level <- 0L
  class(p) <- "cg_printer"
  .cgp_init(p, 1L, if (is.null(settings)) list() else as.list(settings))
  p
}
.cgp_init <- function(p, from, settings) {
  ch <- p$chain
  if (from <= length(ch)) for (k in from:length(ch)) {
    f <- .cgp_classes[[ch[k]]]$init
    if (!is.null(f)) return(invisible(f(p, settings, function(s) .cgp_init(p, k + 1L, s))))
  }
  invisible(NULL)
}

#' @export
print.cg_printer <- function(x, ...) {
  cat("<", x$cls, ">\n", sep = "")
  invisible(x)
}

# Attribute lookup: the instance, then the class attributes along the MRO.
.cgp_a <- function(p, name) {
  if (exists(name, envir = p, inherits = FALSE)) return(get(name, envir = p, inherits = FALSE))
  for (pc in p$chain) {
    a <- .cgp_classes[[pc]]$attrs
    if (name %in% names(a)) return(a[[name]])
  }
  NULL
}
.cgp_has <- function(p, name) {
  if (exists(name, envir = p, inherits = FALSE)) return(TRUE)
  for (pc in p$chain) if (name %in% names(.cgp_classes[[pc]]$attrs)) return(TRUE)
  FALSE
}
# getattr(printer, name), for ast_attr()
.cgp_attr_get <- function(p, name) {
  if (name == "_settings") name <- "settings"
  if (.cgp_has(p, name)) return(.cgp_a(p, name))
  f <- .cgp_get(p, name)
  if (is.function(f)) return(function(...) f(p, ...))
  .cg_err("AttributeError", sprintf("'%s' object has no attribute '%s'", p$cls, name))
}
.cgp_get <- function(p, name) {
  for (pc in p$chain) {
    env <- .cgp_classes[[pc]]$methods
    if (exists(name, envir = env, inherits = FALSE)) return(get(name, envir = env, inherits = FALSE))
  }
  NULL
}
# Call a method by name, as self.<name>(...)
.cgp_call <- function(p, name, ...) {
  f <- .cgp_get(p, name)
  if (!is.function(f)) .cg_err("AttributeError", sprintf("'%s' object has no attribute '%s'", p$cls, name))
  f(p, ...)
}
# super().<name>(...) from the class `from`
.cgp_super <- function(p, from, name, ...) {
  ch <- p$chain
  k <- match(from, ch)
  for (pc in ch[-seq_len(k)]) {
    env <- .cgp_classes[[pc]]$methods
    if (exists(name, envir = env, inherits = FALSE)) {
      f <- get(name, envir = env, inherits = FALSE)
      if (is.function(f)) return(f(p, ...))
    }
  }
  .cg_err("AttributeError", sprintf("'super' object has no attribute '%s'", name))
}

# settings.get(key, default)
.cgp_s <- function(p, key, default = NULL) if (key %in% names(p$settings)) p$settings[[key]] else default

# ---------------------------------------------------------------------------
# Typed dictionaries (keys that are syms or strings), sets of strings

.cgd <- function(pairs = list()) structure(list(pairs = pairs), class = "cg_dict")
.cgd_k <- function(k) if (is_sym(k)) k$key else if (is.character(k)) paste0("\001", k) else paste(deparse(k), collapse = "")
.cgd_as <- function(x) {
  if (inherits(x, "cg_dict")) return(x)
  if (is.null(x) || !length(x)) return(.cgd())
  if (is.null(names(x)) && all(vapply(x, function(pr) is.list(pr) && !is_sym(pr) && length(pr) == 2L, TRUE)))
    return(.cgd(x))
  .cgd(Map(function(k, v) list(k, v), as.list(names(x)), unname(as.list(x))))
}
.cgd_get <- function(d, k, default = NULL) {
  kk <- .cgd_k(k)
  for (pr in d$pairs) if (identical(.cgd_k(pr[[1L]]), kk)) return(pr[[2L]])
  default
}
.cgd_has <- function(d, k) {
  kk <- .cgd_k(k)
  for (pr in d$pairs) if (identical(.cgd_k(pr[[1L]]), kk)) return(TRUE)
  FALSE
}
# dict(chain(a.items(), b.items()))
.cgd_merge <- function(a, b) {
  a <- .cgd_as(a); b <- .cgd_as(b)
  out <- a$pairs
  for (pr in b$pairs) {
    kk <- .cgd_k(pr[[1L]])
    hit <- which(vapply(out, function(q) identical(.cgd_k(q[[1L]]), kk), TRUE))
    if (length(hit)) out[[hit[1L]]] <- pr else out[[length(out) + 1L]] <- pr
  }
  .cgd(out)
}
.cgd_keys <- function(d) lapply(d$pairs, `[[`, 1L)

# A Python set of strings
.cg_pyset <- function(x = character()) structure(unique(as.character(x)), class = "cg_pyset")
#' @export
format.cg_pyset <- function(x, ...) {
  x <- unclass(x)
  if (!length(x)) return("set()")
  paste0("{", paste(vapply(sort(x), .py_repr_str, ""), collapse = ", "), "}")
}
#' @export
print.cg_pyset <- function(x, ...) {
  cat(format(x), "\n")
  invisible(x)
}
.cgp_update <- function(p, field, values) {
  cur <- p[[field]]
  assign(field, .cg_pyset(c(unclass(cur), values)), envir = p)
}

# ---------------------------------------------------------------------------
# The SymPy class MRO of an expression

.CGP_NUM <- c("Number", "AtomicExpr", "Atom", "Expr", "Basic")
.CGP_FUN <- c("Function", "Application", "Expr", "Basic")
.CGP_REL <- c("Relational", "Boolean", "Basic")
.CGP_BOOLF <- c("BooleanFunction", "Application", "Boolean", "Basic")
.CGP_NUMSYM <- c("NumberSymbol", "AtomicExpr", "Atom", "Expr", "Basic")

# Classes between a function class and Function in SymPy (names ending in
# Base survive the filter in Printer._print), and function classes that do
# not derive from Function.
.cgp_fun_bases <- new.env()
.cgp_fun_bases$exp <- "ExpBase"
.cgp_fun_bases$Max <- c("MinMaxBase", "Expr", "LatticeOp", "AssocOp", "Application", "Basic", "NOFUN")
.cgp_fun_bases$Min <- .cgp_fun_bases$Max
for (.cg_n in c("besselj", "bessely", "besseli", "besselk", "hankel1", "hankel2", "jn", "yn"))
  assign(.cg_n, "BesselBase", envir = .cgp_fun_bases)
# the FFunction classes of R/codegen-nodes.R
for (.cg_n in names(.cg_ffunctions))
  assign(.cg_n, c(if (.cg_ffunctions[[.cg_n]] == 95) "F95Function", "FFunction"), envir = .cgp_fun_bases)
rm(.cg_n)
# Function classes with SymPy names different from the internal names
.cgp_fun_names <- c(expm1 = "expm1", log1p = "log1p")

.cgp_is_undef <- function(name) is.null(.FUNS[[name]]) && is.null(.cgp_fun_bases[[name]])

.cgp_mro <- function(e, full = FALSE) {
  if (!is_sym(e)) {
    if (inherits(e, "symmatrix")) return(c("MutableDenseMatrix", "DenseMatrix", "MutableRepMatrix", "RepMatrix", "MatrixBase"))
    if (inherits(e, "cg_dict")) return("dict")
    if (inherits(e, "pytuple")) return("tuple")
    if (is.null(e)) return("NoneType")
    if (is.character(e)) return("str")
    if (is.logical(e)) return(c("bool", "int"))
    if (is.integer(e)) return("int")
    if (is.double(e)) return("float")
    if (is.list(e)) return("list")
    return(class(e)[1L])
  }
  t <- e$type
  cl <- .cg_classes[[t]]
  if (!is.null(cl)) {
    isa <- cl$isa
    if ("String" %in% isa) isa <- append(isa, "Atom", match("String", isa))
    tail <- if (!is.null(cl$tail)) cl$tail else c("CodegenAST", if (isTRUE(cl$expr)) "Expr", "Basic")
    return(c(isa, tail))
  }
  switch(t,
         Integer = {
           v <- e$value
           c(if (v == 0) c("Zero", "IntegerConstant") else if (v == 1) c("One", "IntegerConstant")
             else if (v == -1) c("NegativeOne", "IntegerConstant"), "Integer", "Rational", .CGP_NUM)
         },
         Rational = c(if (e$value == gmp::as.bigq(1, 2)) c("Half", "RationalConstant"), "Rational", .CGP_NUM),
         Float = c(if (!is.null(e$extra$literal)) c(e$extra$literal, "_literal"), "Float", .CGP_NUM),
         Symbol = if (.is_idx(e)) c("Idx", "Expr", "Basic") else c("Symbol", "AtomicExpr", "Atom", "Expr", "Boolean", "Basic"),
         Dummy = c("Dummy", "Symbol", "AtomicExpr", "Atom", "Expr", "Boolean", "Basic"),
         Wild = c("Wild", "Symbol", "AtomicExpr", "Atom", "Expr", "Boolean", "Basic"),
         Constant = switch(e$name,
                           pi = c("Pi", .CGP_NUMSYM), E = c("Exp1", .CGP_NUMSYM),
                           EulerGamma = c("EulerGamma", .CGP_NUMSYM), Catalan = c("Catalan", .CGP_NUMSYM),
                           GoldenRatio = c("GoldenRatio", .CGP_NUMSYM),
                           TribonacciConstant = c("TribonacciConstant", .CGP_NUMSYM),
                           I = c("ImaginaryUnit", "AtomicExpr", "Atom", "Expr", "Basic"),
                           oo = c("Infinity", .CGP_NUM), `-oo` = c("NegativeInfinity", .CGP_NUM),
                           zoo = c("ComplexInfinity", "AtomicExpr", "Atom", "Expr", "Basic"),
                           nan = c("NaN", .CGP_NUM),
                           c(e$name, .CGP_NUMSYM)),
         BooleanAtom = c(if (e$name == "true") "BooleanTrue" else "BooleanFalse", "BooleanAtom", "Boolean", "Basic"),
         Add = c("Add", "Expr", "AssocOp", "Basic"),
         Mul = c("Mul", "Expr", "AssocOp", "Basic"),
         Pow = c("Pow", "Expr", "Basic"),
         Function = .cgp_fun_mro(e, full),
         Piecewise = c("Piecewise", .CGP_FUN),
         Equality = c("Equality", .CGP_REL), Unequality = c("Unequality", .CGP_REL),
         StrictLessThan = c("StrictLessThan", "_Less", "_Inequality", .CGP_REL),
         LessThan = c("LessThan", "_Less", "_Inequality", .CGP_REL),
         StrictGreaterThan = c("StrictGreaterThan", "_Greater", "_Inequality", .CGP_REL),
         GreaterThan = c("GreaterThan", "_Greater", "_Inequality", .CGP_REL),
         And = c("And", "LatticeOp", "AssocOp", .CGP_BOOLF), Or = c("Or", "LatticeOp", "AssocOp", .CGP_BOOLF),
         Not = c("Not", .CGP_BOOLF), Xor = c("Xor", .CGP_BOOLF), Implies = c("Implies", .CGP_BOOLF),
         Equivalent = c("Equivalent", .CGP_BOOLF), ITE = c("ITE", .CGP_BOOLF),
         Tuple = c("Tuple", "Basic"),
         ExprCondPair = c("ExprCondPair", "Tuple", "Basic"),
         Range = c("Range", "Set", "Basic"), Interval = c("Interval", "Set", "Basic"),
         Lambda = c("Lambda", "Expr", "Basic"),
         Indexed = c("Indexed", "Expr", "Basic"), IndexedBase = c("IndexedBase", "Expr", "Basic"),
         MatrixSymbol = c("MatrixSymbol", "MatrixExpr", "Expr", "Basic"),
         MatrixElement = c("MatrixElement", "Expr", "Basic"),
         ImmutableMatrix = c("ImmutableDenseMatrix", "DenseMatrix", "ImmutableRepMatrix", "RepMatrix", "MatrixBase", "MatrixExpr", "Expr", "Basic"),
         MatMul = c("MatMul", "MatrixExpr", "Mul", "Expr", "AssocOp", "Basic"), MatAdd = c("MatAdd", "MatrixExpr", "Add", "Expr", "AssocOp", "Basic"),
         MatPow = c("MatPow", "MatrixExpr", "Expr", "Basic"), Transpose = c("Transpose", "MatrixExpr", "Expr", "Basic"),
         Inverse = c("Inverse", "MatPow", "MatrixExpr", "Expr", "Basic"),
         Derivative = c("Derivative", "Expr", "Basic"),
         Integral = c("Integral", "AddWithLimits", "ExprWithLimits", "Expr", "Basic"),
         Sum = c("Sum", "AddWithLimits", "ExprWithIntLimits", "ExprWithLimits", "Expr", "Basic"),
         UnevaluatedExpr = c("UnevaluatedExpr", "Expr", "Basic"),
         c(t, "Basic"))
}

.cgp_fun_mro <- function(e, full = FALSE) {
  nm <- e$name
  b <- .cgp_fun_bases[[nm]]
  if (!is.null(b) && "NOFUN" %in% b) return(c(nm, setdiff(b, "NOFUN")))
  if (.cgp_is_undef(nm)) return(c(if (full) nm, "AppliedUndef", .CGP_FUN))
  if (full) return(c(nm, b, .CGP_FUN))
  c(nm, b[grepl("Base$", b)], .CGP_FUN)
}

# The class name of an expression: type(expr).__name__
.cgp_clsname <- function(e) if (is_sym(e) && e$type == "Function") e$name else .cgp_mro(e, full = TRUE)[1L]
# expr.func.__name__
.cgp_fname <- .cgp_clsname

# ---------------------------------------------------------------------------
# precedence.py

.CG_PRECEDENCE <- c(Lambda = 1, Xor = 10, Or = 20, And = 30, Relational = 35, Add = 40, Mul = 50,
                    Pow = 60, Func = 70, Not = 100, Atom = 1000, BitwiseOr = 36, BitwiseXor = 37,
                    BitwiseAnd = 38)
.CG_PREC_VALUES <- c(Equivalent = 10, Xor = 10, Implies = 10, Or = 20, And = 30, Add = 40, Pow = 60,
                     Relational = 35, Sub = 40, Not = 100, Function = 70, NegativeInfinity = 40,
                     MatAdd = 40, MatPow = 60, MatrixSolve = 50, Mod = 50, TensAdd = 40, TensMul = 50,
                     HadamardProduct = 50, HadamardPower = 60, KroneckerProduct = 50,
                     Equality = 50, Unequality = 50)

.cgp_num_neg <- function(e) {
  if (e$type %in% c("Integer", "Rational")) return(e$value < 0)
  if (e$type == "Float") return(num_sign(e) < 0)
  FALSE
}

.cgp_prec <- function(e) {
  if (!is_sym(e)) return(.CG_PRECEDENCE[["Atom"]])
  if (e$type == "TypeCast") return(if (isTRUE(e$extra$explicit)) .CG_PRECEDENCE[["Func"]] + 10 else .cgp_prec(e$args[[1L]]))
  cl <- .cg_classes[[e$type]]
  if (!is.null(cl) && !is.null(cl$precedence)) return(cl$precedence)
  for (n in .cgp_mro(e, full = TRUE)) {
    if (n == "Integer") return(if (.cgp_num_neg(e)) 40 else 1000)
    if (n == "Rational") return(if (.cgp_num_neg(e)) 40 else 50)
    if (n == "Float") return(if (.cgp_num_neg(e)) 40 else 1000)
    if (n == "Mul") return(if (could_extract_minus_sign(e)) 40 else 50)
    if (n == "UnevaluatedExpr") return(.cgp_prec(e$args[[1L]]))
    v <- .CG_PREC_VALUES[n]
    if (!is.na(v)) return(unname(v))
  }
  1000
}

# ---------------------------------------------------------------------------
# Printer._print

# Print methods that nodes define for one printer (printmethod such as
# _fcode), keyed by printmethod and class name.
.cgp_nodemethods <- new.env()

.CG_CODE_PM <- c("_ccode", "_cupycode", "_cxxcode", "_fcode", "_jaxcode", "_lambdacode", "_mpmathcode", "_octave",
                 "_pythoncode", "_numpycode", "_scipycode")

.cgpp <- function(p, e, ...) {
  if (is.function(e) && !is.null(attr(e, "base"))) e <- attr(e, "base")
  p$print_level <- p$print_level + 1L
  on.exit(p$print_level <- p$print_level - 1L)
  pm <- .cgp_a(p, "printmethod")
  if (!is.null(pm) && is_sym(e)) {
    tab <- .cgp_nodemethods[[pm]]
    if (!is.null(tab)) for (cn in .cgp_mro(e, full = TRUE)) {
      f <- tab[[cn]]
      if (!is.null(f)) return(f(p, e, ...))
    }
    # a function with a code form (the _print_code printmethods of SymPy's
    # biomechanics curves): that form, parenthesized as an atom
    if (is_sym(e) && e$type == "Function" && pm %in% .CG_CODE_PM) {
      cf <- .FUNS[[e$name]]$code
      if (!is.null(cf)) return(.cgp_paren(p, cf(e), .CG_PRECEDENCE[["Atom"]]))
    }
  }
  for (cn in .cgp_mro(e)) {
    f <- .cgp_get(p, paste0("_print_", cn))
    if (is.function(f)) return(f(p, e, ...))
  }
  .cgp_call(p, "emptyPrinter", e)
}

.cgp_paren <- function(p, item, level, strict = FALSE) {
  pr <- .cgp_prec(item)
  if (pr < level || (!strict && pr <= level)) paste0("(", .cgpp(p, item), ")") else .cgpp(p, item)
}
.cgp_stringify <- function(p, args, sep, level = 0) {
  if (is_sym(args)) args <- args$args
  paste(vapply(args, function(a) .cgp_paren(p, a, level), ""), collapse = sep)
}

# Python's repr() of a float: the shortest digits that read back the value
.cg_float_repr <- function(v) {
  if (is.nan(v)) return("nan")
  if (is.infinite(v)) return(if (v > 0) "inf" else "-inf")
  if (v == 0) return(if (1 / v < 0) "-0.0" else "0.0")
  for (d in 1:17) {
    s <- sprintf("%.*e", d - 1L, v)
    if (as.numeric(s) == v) break
  }
  m <- regmatches(s, regexec("^(-?)([0-9])(?:\\.([0-9]+))?e([+-][0-9]+)$", s))[[1L]]
  sgn <- m[2L]; digits <- paste0(m[3L], m[4L]); ex <- as.integer(m[5L])
  digits <- sub("0+$", "", digits)
  if (!nzchar(digits)) digits <- "0"
  if (ex < -4L || ex >= 16L) {
    mant <- if (nchar(digits) > 1L) paste0(substr(digits, 1L, 1L), ".", substring(digits, 2L)) else digits
    return(sprintf("%s%se%s%02d", sgn, mant, if (ex < 0L) "-" else "+", abs(ex)))
  }
  if (ex >= 0L) {
    if (nchar(digits) <= ex + 1L) {
      ip <- paste0(digits, strrep("0", ex + 1L - nchar(digits))); fp <- "0"
    } else {
      ip <- substr(digits, 1L, ex + 1L); fp <- substring(digits, ex + 2L)
    }
  } else {
    ip <- "0"; fp <- paste0(strrep("0", -ex - 1L), digits)
  }
  paste0(sgn, ip, ".", fp)
}

# str() of a Python value or of a node
.py_str_of <- function(v) {
  if (is.null(v)) return("None")
  if (is_sym(v)) return(sstr(v))
  if (is.character(v)) return(v)
  if (is.logical(v)) return(if (isTRUE(v)) "True" else "False")
  if (is.integer(v)) return(as.character(v))
  if (is.double(v)) return(.cg_float_repr(v))
  if (inherits(v, "symmatrix")) return(sstr(v))
  paste(format(v), collapse = "")
}
# repr() of a node: Token nodes have their own repr, other nodes print as str
.py_repr_of <- function(v) {
  if (is_sym(v)) {
    cl <- .cg_classes[[v$type]]
    if (!is.null(cl) && isTRUE(cl$token)) return(.cg_print(v, TRUE))
    return(sstr(v))
  }
  if (is.character(v)) return(.py_repr_str(v))
  .py_str_of(v)
}

# ---------------------------------------------------------------------------
# Printer, and the StrPrinter methods that the code printers inherit

.cgp_def("Printer", attrs = list(`_default_settings` = list(), printmethod = NULL),
         init = function(p, settings, nxt) {
           defaults <- .cgp_a(p, "_default_settings")
           s <- defaults
           for (k in names(settings)) s[k] <- list(settings[[k]])
           if (length(s) > length(defaults)) for (k in names(s))
             if (!(k %in% names(defaults))) .cg_err("TypeError", sprintf("Unknown setting '%s'.", k))
           p$settings <- s
           p$print_level <- 0L
           invisible(NULL)
         },
         methods = list(
           emptyPrinter = function(p, e) .py_str_of(e),
           doprint = function(p, e) .cgpp(p, e)
         ))

.cgp_ordered_terms <- function(p, e, order = NULL) {
  if (is.null(order)) order <- .cgp_s(p, "order")
  if (identical(order, "none") || !(e$type %in% c("Add", "MatAdd"))) return(if (e$type %in% c("Add", "MatAdd")) e$args else list(e))
  .add_print_order(e$args)
}

# MatMul.as_coeff_matrices: the product of the scalar factors and the matrices
.cgp_coeff_matrices <- function(e) {
  ism <- vapply(e$args, function(a) .is_matexpr(a) || inherits(a, "symmatrix"), NA)
  list(coeff = if (any(!ism)) .mul(e$args[!ism]) else S_One, matrices = e$args[ism])
}

.cgp_is_add <- function(e) is_sym(e) && (e$type == "Add" || (e$type == "UnevaluatedExpr" && e$args[[1L]]$type == "Add"))

.cgp_def("StrPrinter", "Printer",
         attrs = list(printmethod = "_sympystr",
                      `_default_settings` = list(order = NULL, full_prec = "auto", sympy_integers = FALSE,
                                                 abbrev = FALSE, perm_cyclic = TRUE, min = NULL, max = NULL,
                                                 dps = NULL),
                      `_relationals` = list()),
         methods = list(
           emptyPrinter = function(p, e) {
             if (is.character(e)) return(e)
             if (is_sym(e)) return(.py_repr_of(e))
             .py_str_of(e)
           },
           `_print_Add` = function(p, e, order = NULL) {
             terms <- .cgp_ordered_terms(p, e, order)
             prec <- .cgp_prec(e)
             l <- character()
             for (term in terms) {
               t <- .cgpp(p, term)
               if (startsWith(t, "-") && !.cgp_is_add(term)) { sign <- "-"; t <- substring(t, 2L) } else sign <- "+"
               if (.cgp_prec(term) < prec || .cgp_is_add(term)) l <- c(l, sign, paste0("(", t, ")"))
               else l <- c(l, sign, t)
             }
             sign <- l[1L]
             l <- l[-1L]
             if (sign == "+") sign <- ""
             paste0(sign, paste(l, collapse = " "))
           },
           `_print_MatMul` = function(p, e) {
             c_m <- .cgp_coeff_matrices(e)
             args <- e$args
             sign <- ""
             prec <- .cgp_prec(e)
             if (.expr_is_number(c_m$coeff)) {
               v <- tryCatch(.neval(c_m$coeff, list()), error = function(err) NA)
               if (!is.na(v) && ((Re(v) < 0 && Im(v) == 0) || (Re(v) == 0 && Im(v) < 0))) {
                 nc <- .neg(c_m$coeff)
                 args <- if (identical(nc$key, S_One$key)) c_m$matrices else c(list(nc), c_m$matrices)
                 sign <- "-"
                 prec <- .CG_PRECEDENCE[["Mul"]]
               }
             }
             paste0(sign, paste(vapply(args, function(a) .cgp_paren(p, a, prec), ""), collapse = "*"))
           },
           `_print_BooleanTrue` = function(p, e) "True",
           `_print_BooleanFalse` = function(p, e) "False",
           `_print_Not` = function(p, e) paste0("~", .cgp_paren(p, e$args[[1L]], .CG_PRECEDENCE[["Not"]])),
           `_print_And` = function(p, e) .cgp_stringify(p, e$args, " & ", .CG_PRECEDENCE[["BitwiseAnd"]]),
           `_print_Or` = function(p, e) .cgp_stringify(p, e$args, " | ", .CG_PRECEDENCE[["BitwiseOr"]]),
           `_print_Xor` = function(p, e) .cgp_stringify(p, e$args, " ^ ", .CG_PRECEDENCE[["BitwiseXor"]]),
           `_print_Basic` = function(p, e) paste0(.cgp_clsname(e), "(", paste(vapply(e$args, function(o) .cgpp(p, o), ""), collapse = ", "), ")"),
           `_print_Catalan` = function(p, e) "Catalan",
           `_print_ComplexInfinity` = function(p, e) "zoo",
           `_print_Dummy` = function(p, e) paste0("_", e$name),
           `_print_EulerGamma` = function(p, e) "EulerGamma",
           `_print_Exp1` = function(p, e) "E",
           `_print_ExprCondPair` = function(p, e) paste0("(", .cgpp(p, e$args[[1L]]), ", ", .cgpp(p, e$args[[2L]]), ")"),
           `_print_Function` = function(p, e) paste0(.cgp_fname(e), "(", .cgp_stringify(p, e$args, ", "), ")"),
           `_print_GoldenRatio` = function(p, e) "GoldenRatio",
           `_print_TribonacciConstant` = function(p, e) "TribonacciConstant",
           `_print_ImaginaryUnit` = function(p, e) "I",
           `_print_Infinity` = function(p, e) "oo",
           `_print_Lambda` = function(p, e) {
             sig <- e$args[[1L]]
             if (is_sym(sig) && sig$type == "Tuple" && length(sig$args) == 1L) sig <- sig$args[[1L]]
             paste0("Lambda(", .cgpp(p, sig), ", ", .cgpp(p, e$args[[2L]]), ")")
           },
           `_print_list` = function(p, e) paste0("[", .cgp_stringify(p, e, ", "), "]"),
           `_print_MatrixBase` = function(p, e) sstr(e),
           `_print_MatrixElement` = function(p, e)
             paste0(.cgp_paren(p, e$args[[1L]], .CG_PRECEDENCE[["Atom"]], strict = TRUE), "[",
                    .cgpp(p, e$args[[2L]]), ", ", .cgpp(p, e$args[[3L]]), "]"),
           `_print_NaN` = function(p, e) "nan",
           `_print_NegativeInfinity` = function(p, e) "-oo",
           `_print_Pi` = function(p, e) "pi",
           `_print_Pow` = function(p, e, rational = FALSE) {
             PREC <- .cgp_prec(e)
             b <- e$args[[1L]]; x <- e$args[[2L]]
             if (x$type == "Rational" && x$value == gmp::as.bigq(1, 2) && !rational) return(paste0("sqrt(", .cgpp(p, b), ")"))
             if (.is_comm(e)) {
               if (x$type == "Rational" && x$value == gmp::as.bigq(-1, 2) && !rational)
                 return(paste0(.cgpp(p, S_One), "/sqrt(", .cgpp(p, b), ")"))
               if (x$type == "Integer" && x$value == -1)
                 return(paste0(.cgpp(p, S_One), "/", .cgp_paren(p, b, PREC, strict = FALSE)))
             }
             paste0(.cgp_paren(p, b, PREC, strict = FALSE), "**", .cgp_paren(p, x, PREC, strict = FALSE))
           },
           `_print_UnevaluatedExpr` = function(p, e) .cgpp(p, e$args[[1L]]),
           `_print_Integer` = function(p, e) {
             if (isTRUE(.cgp_s(p, "sympy_integers", FALSE))) return(paste0("S(", as.character(e$value), ")"))
             as.character(gmp::numerator(e$value))
           },
           `_print_int` = function(p, e) as.character(e),
           `_print_float` = function(p, e) .cg_float_repr(e),
           `_print_bool` = function(p, e) if (isTRUE(e)) "True" else "False",
           `_print_NoneType` = function(p, e) "None",
           `_print_Rational` = function(p, e) {
             if (gmp::denominator(e$value) == 1) return(as.character(gmp::numerator(e$value)))
             if (isTRUE(.cgp_s(p, "sympy_integers", FALSE)))
               return(sprintf("S(%s)/%s", as.character(gmp::numerator(e$value)), as.character(gmp::denominator(e$value))))
             paste0(as.character(gmp::numerator(e$value)), "/", as.character(gmp::denominator(e$value)))
           },
           `_print_Float` = function(p, e) {
             fp <- .cgp_s(p, "full_prec", "auto")
             strip <- if (isTRUE(fp)) FALSE else if (isFALSE(fp)) TRUE else p$print_level > 1L
             .fmt_float_sym(e, strip = strip)
           },
           `_print_Relational` = function(p, e) {
             charmap <- c(Equality = "Eq", Unequality = "Ne")
             t <- e$type
             if (t %in% names(charmap))
               return(sprintf("%s(%s, %s)", charmap[[t]], .cgpp(p, e$args[[1L]]), .cgpp(p, e$args[[2L]])))
             op <- .REL_OPS[[t]]
             rel <- .cgp_a(p, "_relationals")
             if (!is.null(rel[[op]])) op <- rel[[op]]
             paste(.cgp_paren(p, e$args[[1L]], .cgp_prec(e)), op, .cgp_paren(p, e$args[[2L]], .cgp_prec(e)))
           },
           `_print_Symbol` = function(p, e) e$name,
           `_print_MatrixSymbol` = function(p, e) e$name,
           `_print_str` = function(p, e) e,
           `_print_tuple` = function(p, e) {
             if (length(e) == 1L) return(paste0("(", .cgpp(p, e[[1L]]), ",)"))
             paste0("(", .cgp_stringify(p, unclass(e), ", "), ")")
           },
           `_print_Tuple` = function(p, e) {
             if (length(e$args) == 1L) return(paste0("(", .cgpp(p, e$args[[1L]]), ",)"))
             paste0("(", .cgp_stringify(p, e$args, ", "), ")")
           },
           `_print_Zero` = function(p, e) {
             if (isTRUE(.cgp_s(p, "sympy_integers", FALSE))) return("S(0)")
             "0"
           },
           `_print_dict` = function(p, e) {
             ks <- .cgd_keys(e)
             ks <- ks[order(vapply(ks, function(k) .py_str_of(k), ""))]
             paste0("{", paste(vapply(ks, function(k) paste0(.cgpp(p, k), ": ", .cgpp(p, .cgd_get(e, k))), ""), collapse = ", "), "}")
           }
         ))

# ---------------------------------------------------------------------------
# CodePrinter

.CG_REWRITEABLE <- list(cot = list("tan", character()), csc = list("sin", character()), sec = list("cos", character()),
                        acot = list("atan", character()), acsc = list("asin", character()), asec = list("acos", character()),
                        coth = list("exp", character()), csch = list("exp", character()), sech = list("exp", character()),
                        acoth = list("log", character()), acsch = list("log", character()), asech = list("log", character()),
                        catalan = list("gamma", character()), fibonacci = list("sqrt", character()),
                        lucas = list("sqrt", character()), beta = list("gamma", character()),
                        sinc = list("sin", "Piecewise"), Mod = list("floor", character()),
                        factorial = list("gamma", character()), factorial2 = list("gamma", "Piecewise"),
                        subfactorial = list("uppergamma", character()), RisingFactorial = list("gamma", "Piecewise"),
                        FallingFactorial = list("gamma", "Piecewise"), binomial = list("gamma", character()),
                        frac = list("floor", character()), Max = list("Piecewise", character()),
                        Min = list("Piecewise", character()), Heaviside = list("Piecewise", character()),
                        erf2 = list("erf", character()), erfc = list("erf", character()), Li = list("li", character()),
                        Ei = list("li", character()), dirichlet_eta = list("zeta", character()),
                        riemann_xi = list("zeta", "gamma"), SingularityFunction = list("Piecewise", character()))

# The module path of a class, for messages that print type(x)
.cgp_type_str <- function(x) {
  if (inherits(x, "cg_printer")) return(sprintf("<class '%s.%s'>", .cgp_a(x, "module"), x$cls))
  if (is_sym(x)) {
    cl <- .cg_classes[[x$type]]
    if (!is.null(cl)) return(sprintf("<class '%s.%s'>", if (is.null(cl$module)) "sympy.codegen.ast" else cl$module, x$type))
    m <- .cg_pytype(x)
    if (!grepl(".", m, fixed = TRUE)) {
      cn <- .cgp_clsname(x)
      if (x$type == "Function" && .cgp_is_undef(x$name)) return(sprintf("%s", cn))
      return(sprintf("<class '%s'>", cn))
    }
    return(sprintf("<class '%s'>", m))
  }
  sprintf("<class '%s'>", .cgp_mro(x)[1L])
}

# Python lists in the expression become List nodes
.cgp_convert_lists <- function(arg) {
  if (is_sym(arg) || inherits(arg, "symmatrix") || is.function(arg)) return(arg)
  if (inherits(arg, "pytuple")) return(structure(lapply(unclass(arg), .cgp_convert_lists), class = "pytuple"))
  if (is.list(arg) && !inherits(arg, "cg_dict")) return(do.call(cg_List, lapply(arg, .cgp_convert_lists)))
  arg
}

.cgp_sympify <- function(x) {
  if (is_sym(x) || inherits(x, "symmatrix")) return(x)
  if (is.function(x) && !is.null(attr(x, "base"))) return(attr(x, "base"))
  if (inherits(x, "pytuple")) return(do.call(Tuple, lapply(unclass(x), .cgp_sympify)))
  if (is.character(x)) return(sympify(x))
  as_sym(x)
}

.cgp_is_matrix <- function(e) inherits(e, "symmatrix") || (is_sym(e) && (e$type == "ImmutableMatrix" || .is_matexpr(e)))
.cgp_mshape <- function(e) {
  if (inherits(e, "symmatrix")) return(c(e$nrow, e$ncol))
  if (e$type == "ImmutableMatrix") return(c(.int_val(e$args[[1L]]), .int_val(e$args[[2L]])))
  vapply(.mshape(e), function(d) as.integer(.int_val(d)), 0L)
}
# m[i, j] of a matrix or matrix expression (0-based)
.cgp_mget <- function(m, i, j) {
  if (inherits(m, "symmatrix")) return(m[i + 1L, j + 1L])
  if (m$type == "ImmutableMatrix") {
    nc <- .int_val(m$args[[2L]])
    return(m$args[[3L]]$args[[i * nc + j + 1L]])
  }
  MatrixElement(m, i, j)
}

.cgp_handle_assign_to <- function(p, expr, assign_to) {
  if (is.null(assign_to)) return(.cgp_sympify(expr))
  if (is.list(assign_to) && !is_sym(assign_to)) {
    ex <- if (is_sym(expr)) expr$args else as.list(expr)
    if (length(ex) != length(assign_to))
      .cg_err("ValueError", sprintf("Failed to assign an expression of length %d to %d variables", length(ex), length(assign_to)))
    return(do.call(CodeBlock, Map(function(l, r) .cgp_handle_assign_to(p, l, r), ex, assign_to)))
  }
  if (is.function(assign_to) && !is.null(attr(assign_to, "base"))) assign_to <- attr(assign_to, "base")
  if (is.character(assign_to)) {
    if (.cgp_is_matrix(expr)) {
      sh <- .cgp_mshape(expr)
      assign_to <- MatrixSymbol(assign_to, sh[1L], sh[2L])
    } else assign_to <- Symbol(assign_to)
  } else if (!is_sym(assign_to)) {
    .cg_err("TypeError", sprintf("%s cannot assign to object of type %s", p$cls, .cgp_type_str(assign_to)))
  }
  Assignment(assign_to, expr)
}

.py_splitlines <- function(s) if (!nzchar(s)) character() else strsplit(s, "\n", fixed = TRUE)[[1L]]

.cgp_set_add <- function(p, field, item) {
  cur <- p[[field]]
  k <- if (is_sym(item)) item$key else paste(vapply(item, .cgd_k, ""), collapse = "|")
  if (!(k %in% names(cur))) { cur[[k]] <- item; assign(field, cur, envir = p) }
}

.cgp_doprint <- function(p, expr, assign_to = NULL) {
  expr <- .cgp_convert_lists(expr)
  expr <- .cgp_handle_assign_to(p, expr, assign_to)
  p$not_supported <- list()
  p$number_symbols <- list()
  lines <- .py_splitlines(.cgpp(p, expr))
  if (isTRUE(.cgp_s(p, "human", TRUE))) {
    front <- character()
    ns <- p$not_supported
    if (length(ns)) {
      front <- c(front, .cgp_call(p, "_get_comment", sprintf("Not supported in %s:", .cgp_a(p, "language"))))
      ns <- ns[order(vapply(ns, sstr, ""))]
      for (z in ns) front <- c(front, .cgp_call(p, "_get_comment", .cgp_clsname(z)))
    }
    nsy <- p$number_symbols
    nsy <- nsy[order(vapply(nsy, function(q) paste0("(", sstr(q[[1L]]), ", ", sstr(q[[2L]]), ")"), ""))]
    for (q in nsy) front <- c(front, .cgp_call(p, "_declare_number_const", q[[1L]], q[[2L]]))
    lines <- c(front, lines)
    lines <- .cgp_call(p, "_format_code", lines)
    result <- paste(lines, collapse = "\n")
  } else {
    lines <- .cgp_call(p, "_format_code", lines)
    num_syms <- lapply(p$number_symbols, function(q) structure(list(q[[1L]], .cgpp(p, q[[2L]])), class = "pytuple"))
    result <- structure(list(structure(unname(num_syms), class = "cg_pysetobj"),
                             structure(unname(p$not_supported), class = "cg_pysetobj"),
                             paste(lines, collapse = "\n")), class = "pytuple")
  }
  p$not_supported <- list()
  p$number_symbols <- list()
  result
}

# a Python set of objects
#' @export
format.cg_pysetobj <- function(x, ...) {
  x <- unclass(x)
  if (!length(x)) return("set()")
  els <- vapply(x, function(v) {
    if (inherits(v, "pytuple")) paste0("(", paste(vapply(unclass(v), .py_repr_of, ""), collapse = ", "), ")")
    else .py_repr_of(v)
  }, "")
  paste0("{", paste(sort(els), collapse = ", "), "}")
}

.cgp_not_supported <- function(p, e) {
  if (isTRUE(.cgp_s(p, "strict", FALSE)))
    .cg_err("PrintMethodNotImplementedError",
            paste0("Unsupported by ", .cgp_type_str(p), ": ", .cgp_type_str(e),
                   "\nSet the printer option 'strict' to False in order to generate partially printed code."))
  if (is_sym(e)) .cgp_set_add(p, "not_supported", e)
  .cgp_call(p, "emptyPrinter", e)
}

.cgp_can_print <- function(p, name) name %in% names(p$known_functions) || is.function(.cgp_get(p, paste0("_print_", name)))

# expr.as_coeff_Mul()
.cgp_coeff_mul <- function(e) {
  if (e$type == "Mul" && is_number_type(e$args[[1L]])) {
    rest <- e$args[-1L]
    return(list(e$args[[1L]], if (length(rest) == 1L) rest[[1L]] else new_sym("Mul", rest)))
  }
  if (is_number_type(e)) return(list(e, S_One))
  list(S_One, e)
}
# core.mul._keep_coeff
.cgp_keep_coeff <- function(coeff, factors) {
  if (!is_number_type(coeff)) {
    if (is_number_type(factors)) { tmp <- factors; factors <- coeff; coeff <- tmp }
    else return(.mul(list(coeff, factors)))
  }
  if (is_one_num(factors) && factors$type == "Integer") return(coeff)
  if (is_one_num(coeff) && coeff$type == "Integer") return(factors)
  if (coeff$type == "Integer" && coeff$value == -1) return(.neg(factors))
  if (factors$type == "Add") return(new_sym("Mul", list(coeff, factors)))
  if (factors$type == "Mul") {
    margs <- factors$args
    if (is_number_type(margs[[1L]])) {
      margs[[1L]] <- .mul(list(margs[[1L]], coeff))
      if (is_one_num(margs[[1L]]) && margs[[1L]]$type == "Integer") margs <- margs[-1L]
    } else margs <- c(list(coeff), margs)
    return(if (length(margs) == 1L) margs[[1L]] else new_sym("Mul", margs))
  }
  m <- .mul(list(coeff, factors))
  if (is_number_type(m) && !is_number_type(factors)) m <- new_sym("Mul", list(coeff, factors))
  m
}

.cgp_print_mul <- function(p, e) {
  prec <- .CG_PRECEDENCE[["Mul"]]
  cm <- .cgp_coeff_mul(e)
  c1 <- cm[[1L]]
  if (.cgp_num_neg(c1)) {
    e <- .cgp_keep_coeff(.neg(c1), cm[[2L]])
    sign <- "-"
  } else sign <- ""
  ord <- .cgp_s(p, "order")
  args <- if (!(identical(ord, "old") || identical(ord, "none"))) {
    if (e$type == "Mul") .mul_print_order(e$args) else list(e)
  } else if (e$type == "Mul") e$args else list(e)
  a <- list(); b <- list()
  for (item in args) {
    if (item$type == "Pow" && .is_comm(item) && item$args[[2L]]$type %in% c("Integer", "Rational") &&
        item$args[[2L]]$value < 0) {
      ex <- item$args[[2L]]
      if (!(ex$type == "Integer" && ex$value == -1)) b[[length(b) + 1L]] <- new_sym("Pow", list(item$args[[1L]], .neg(ex)))
      else b[[length(b) + 1L]] <- item$args[[1L]]
    } else a[[length(a) + 1L]] <- item
  }
  if (!length(a)) a <- list(S_One)
  if (length(a) == 1L && sign == "-") a_str <- .cgp_paren(p, a[[1L]], 0.5 * (.CG_PRECEDENCE[["Pow"]] + .CG_PRECEDENCE[["Mul"]]))
  else a_str <- vapply(a, function(z) .cgp_paren(p, z, prec), "")
  b_str <- vapply(b, function(z) .cgp_paren(p, z, prec), "")
  if (!length(b)) return(paste0(sign, paste(a_str, collapse = "*")))
  if (length(b) == 1L) return(paste0(sign, paste(a_str, collapse = "*"), "/", b_str[1L]))
  paste0(sign, paste(a_str, collapse = "*"), "/(", paste(b_str, collapse = "*"), ")")
}

# Look a function up in known_functions: the printed call, or NULL
.cgp_known_call <- function(p, e, cond_func) {
  if (is.character(cond_func)) return(paste0(cond_func, "(", .cgp_stringify(p, e$args, ", "), ")"))
  func <- NULL
  for (cf in cond_func) {
    func <- cf[[2L]]
    if (isTRUE(.cgp_truth(do.call(cf[[1L]], e$args)))) break
  }
  if (is.null(func)) return(NULL)
  if (is.function(func)) {
    r <- do.call(func, lapply(e$args, function(item) .cgp_paren(p, item, 0)))
    return(if (is_sym(r)) sstr(r) else r)
  }
  paste0(func, "(", .cgp_stringify(p, e$args, ", "), ")")
}
# bool() of a Python value
# Python truthiness: Basic objects are true, except false, zero and empty
# containers
.cgp_truth <- function(v) {
  if (is.logical(v)) return(isTRUE(v))
  if (is.null(v)) return(FALSE)
  if (is_sym(v)) {
    if (v$type == "BooleanAtom") return(identical(v$key, S_true$key))
    if (is_number_type(v)) return(!is_zero_num(v))
    if (v$type %in% c("Tuple", "List")) return(length(v$args) > 0L)
    return(TRUE)
  }
  if (is.numeric(v)) return(v != 0)
  length(v) > 0
}

.cgp_print_function <- function(p, e) {
  nm <- .cgp_fname(e)
  kf <- p$known_functions
  imp <- if (is_sym(e) && e$type == "Function") e$extra$imp else NULL
  if (nm %in% names(kf)) {
    r <- .cgp_known_call(p, e, kf[[nm]])
    if (!is.null(r)) return(r)
  } else if (!is.null(imp) && is_sym(imp) && imp$type == "Lambda") {
    return(.cgpp(p, .cg_apply_lambda(imp, e$args)))
  } else if (nm %in% names(.CG_REWRITEABLE)) {
    rw <- .CG_REWRITEABLE[[nm]]
    if (.cgp_can_print(p, rw[[1L]]) && all(vapply(rw[[2L]], function(f) .cgp_can_print(p, f), TRUE)))
      return(paste0("(", .cgpp(p, rewrite(e, rw[[1L]])), ")"))
  }
  if (.cgp_is_function(e) && isTRUE(.cgp_s(p, "allow_unknown_functions", FALSE)))
    return(paste0(nm, "(", paste(vapply(e$args, function(a) .cgpp(p, a), ""), collapse = ", "), ")"))
  .cgp_not_supported(p, e)
}
.cgp_is_function <- function(e) is_sym(e) && ("Application" %in% .cgp_mro(e, full = TRUE))

# Lambda(sig, body)(*args)
.cg_apply_lambda <- function(L, args) {
  sig <- L$args[[1L]]
  vars <- if (sig$type == "Tuple") sig$args else list(sig)
  subs(L$args[[2L]], Map(function(v, a) list(v, a), vars, args), simultaneous = TRUE)
}

.cgp_print_numbersymbol <- function(p, e) {
  if (isTRUE(.cgp_s(p, "inline", FALSE))) return(.cgpp(p, N(e, .cgp_s(p, "precision"))))
  .cgp_set_add(p, "number_symbols", list(e, N(e, .cgp_s(p, "precision"))))
  sstr(e)
}

.cgp_sorted_args <- function(args) .key_sort(lapply(args, .sort_key), args)

.cgp_def("CodePrinter", "StrPrinter",
         attrs = list(`_operators` = list(and = "&&", or = "||", not = "!"),
                      `_default_settings` = list(order = NULL, full_prec = "auto", error_on_reserved = FALSE,
                                                 reserved_word_suffix = "_", human = TRUE, inline = FALSE,
                                                 allow_unknown_functions = FALSE, strict = NULL),
                      reserved_words = character(), known_functions = list(), module = "sympy.printing.codeprinter"),
         init = function(p, settings, nxt) {
           nxt(settings)
           if ("strict" %in% names(p$settings) && is.null(p$settings$strict))
             p$settings["strict"] <- list(isTRUE(.cgp_s(p, "human", TRUE)))
           else if (!("strict" %in% names(p$settings))) p$settings["strict"] <- list(TRUE)
           p$not_supported <- list()
           p$number_symbols <- list()
           invisible(NULL)
         },
         methods = list(
           doprint = function(p, e, assign_to = NULL) .cgp_doprint(p, e, assign_to),
           `_print_Dummy` = function(p, e) {
             if (startsWith(e$name, "Dummy_")) paste0("_", e$name) else sprintf("%s_%d", e$name, e$extra$id)
           },
           `_print_Idx` = function(p, e) .cgpp(p, .cg_idx_label(e)),
           `_print_CodeBlock` = function(p, e) paste(vapply(e$args, function(i) .cgpp(p, i), ""), collapse = "\n"),
           `_print_String` = function(p, e) e$name,
           `_print_QuotedString` = function(p, e) sprintf("\"%s\"", e$name),
           `_print_Comment` = function(p, e) .cgp_call(p, "_get_comment", e$name),
           `_print_Assignment` = function(p, e) {
             lhs <- e$args[[1L]]; rhs <- e$args[[2L]]
             if (rhs$type == "Piecewise") {
               pairs <- lapply(rhs$args, function(pr) list(Assignment(lhs, pr$args[[1L]]), pr$args[[2L]]))
               return(.cgpp(p, do.call(Piecewise, pairs)))
             }
             if (lhs$type == "MatrixSymbol") {
               lines <- character()
               for (ij in .cgp_call(p, "_traverse_matrix_indices", lhs)) {
                 temp <- Assignment(MatrixElement(lhs, ij[1L], ij[2L]), .cgp_mget(rhs, ij[1L], ij[2L]))
                 lines <- c(lines, .cgpp(p, temp))
               }
               return(paste(lines, collapse = "\n"))
             }
             if (isTRUE(.cgp_s(p, "contract", FALSE)) && (.cg_has_type(lhs, "IndexedBase") || .cg_has_type(rhs, "IndexedBase")))
               return(.cgp_doprint_loops(p, rhs, lhs))
             .cgp_call(p, "_get_statement", sprintf("%s = %s", .cgpp(p, lhs), .cgpp(p, rhs)))
           },
           `_print_AugmentedAssignment` = function(p, e) {
             lhs_code <- .cgpp(p, e$args[[1L]]); rhs_code <- .cgpp(p, e$args[[2L]])
             .cgp_call(p, "_get_statement", paste(.cgpp(p, lhs_code), .cgpp(p, .cg_aug_op(e)), .cgpp(p, rhs_code)))
           },
           `_print_FunctionCall` = function(p, e)
             sprintf("%s(%s)", .cg_text(e$args[[1L]]), paste(vapply(e$args[[2L]]$args, function(a) .cgpp(p, a), ""), collapse = ", ")),
           `_print_Variable` = function(p, e) .cgpp(p, e$args[[1L]]),
           `_print_Symbol` = function(p, e) {
             name <- .cgp_super(p, "CodePrinter", "_print_Symbol", e)
             if (name %in% .cgp_a(p, "reserved_words")) {
               if (isTRUE(.cgp_s(p, "error_on_reserved")))
                 .cg_err("ValueError", sprintf("This expression includes the symbol \"%s\" which is a reserved keyword in this language.", name))
               return(paste0(name, .cgp_s(p, "reserved_word_suffix")))
             }
             name
           },
           `_print_Function` = function(p, e) .cgp_print_function(p, e),
           `_print_Expr` = function(p, e) .cgp_print_function(p, e),
           `_print_Derivative` = function(p, e) {
             obj <- e$args[[1L]]
             for (fa in obj$args) if (fa$type != "Symbol")
               .cg_err("ValueError", sprintf("%s._print_Derivative(...) only supports functions with symbols as arguments.", p$cls))
             meth <- paste0("_print_Derivative_", .cgp_fname(obj))
             pm <- .cgp_get(p, meth)
             if (!is.function(pm)) {
               if (isTRUE(.cgp_s(p, "strict", FALSE)))
                 .cg_err("PrintMethodNotImplementedError",
                         paste0("Unsupported by ", .cgp_type_str(p), ": ", .cgp_type_str(e), "\nPrinter has no method: ", meth,
                                "\nSet the printer option 'strict' to False in order to generate partially printed code."))
               return(.cgp_not_supported(p, e))
             }
             orders <- list()
             for (pr in e$args[-1L]) orders[[pr$args[[1L]]$key]] <- pr$args[[2L]]
             pm(p, obj$args, lapply(obj$args, function(a) orders[[a$key]]))
           },
           `_print_Heaviside` = FALSE,
           `_print_NumberSymbol` = function(p, e) .cgp_print_numbersymbol(p, e),
           `_print_Catalan` = function(p, e) .cgp_call(p, "_print_NumberSymbol", e),
           `_print_EulerGamma` = function(p, e) .cgp_call(p, "_print_NumberSymbol", e),
           `_print_GoldenRatio` = function(p, e) .cgp_call(p, "_print_NumberSymbol", e),
           `_print_TribonacciConstant` = function(p, e) .cgp_call(p, "_print_NumberSymbol", e),
           `_print_Exp1` = function(p, e) .cgp_call(p, "_print_NumberSymbol", e),
           `_print_Pi` = function(p, e) .cgp_call(p, "_print_NumberSymbol", e),
           `_print_And` = function(p, e) {
             PREC <- .cgp_prec(e)
             paste(vapply(.cgp_sorted_args(e$args), function(a) .cgp_paren(p, a, PREC), ""),
                   collapse = sprintf(" %s ", .cgp_a(p, "_operators")[["and"]]))
           },
           `_print_Or` = function(p, e) {
             PREC <- .cgp_prec(e)
             paste(vapply(.cgp_sorted_args(e$args), function(a) .cgp_paren(p, a, PREC), ""),
                   collapse = sprintf(" %s ", .cgp_a(p, "_operators")[["or"]]))
           },
           `_print_Xor` = function(p, e) {
             op <- .cgp_a(p, "_operators")[["xor"]]
             if (is.null(op)) return(.cgpp(p, to_nnf(e)))
             PREC <- .cgp_prec(e)
             paste(vapply(e$args, function(a) .cgp_paren(p, a, PREC), ""), collapse = sprintf(" %s ", op))
           },
           `_print_Equivalent` = function(p, e) {
             op <- .cgp_a(p, "_operators")[["equivalent"]]
             if (is.null(op)) return(.cgpp(p, to_nnf(e)))
             PREC <- .cgp_prec(e)
             paste(vapply(e$args, function(a) .cgp_paren(p, a, PREC), ""), collapse = sprintf(" %s ", op))
           },
           `_print_Not` = function(p, e) paste0(.cgp_a(p, "_operators")[["not"]], .cgp_paren(p, e$args[[1L]], .cgp_prec(e))),
           `_print_BooleanFunction` = function(p, e) .cgpp(p, to_nnf(e)),
           `_print_isnan` = function(p, e) sprintf("isnan(%s)", .cgpp(p, e$args[[1L]])),
           `_print_isinf` = function(p, e) sprintf("isinf(%s)", .cgpp(p, e$args[[1L]])),
           `_print_Mul` = function(p, e) .cgp_print_mul(p, e),
           `_print_not_supported` = function(p, e) .cgp_not_supported(p, e),
           `_print_Basic` = function(p, e) .cgp_not_supported(p, e),
           `_print_ComplexInfinity` = function(p, e) .cgp_not_supported(p, e),
           `_print_ExprCondPair` = function(p, e) .cgp_not_supported(p, e),
           `_print_GeometryEntity` = function(p, e) .cgp_not_supported(p, e),
           `_print_Infinity` = function(p, e) .cgp_not_supported(p, e),
           `_print_Integral` = function(p, e) .cgp_not_supported(p, e),
           `_print_Interval` = function(p, e) .cgp_not_supported(p, e),
           `_print_AccumulationBounds` = function(p, e) .cgp_not_supported(p, e),
           `_print_Limit` = function(p, e) .cgp_not_supported(p, e),
           `_print_MatrixBase` = function(p, e) .cgp_not_supported(p, e),
           `_print_DeferredVector` = function(p, e) .cgp_not_supported(p, e),
           `_print_NaN` = function(p, e) .cgp_not_supported(p, e),
           `_print_NegativeInfinity` = function(p, e) .cgp_not_supported(p, e),
           `_print_Order` = function(p, e) .cgp_not_supported(p, e),
           `_print_RootOf` = function(p, e) .cgp_not_supported(p, e),
           `_print_RootsOf` = function(p, e) .cgp_not_supported(p, e),
           `_print_RootSum` = function(p, e) .cgp_not_supported(p, e),
           `_print_Uniform` = function(p, e) .cgp_not_supported(p, e),
           `_print_Unit` = function(p, e) .cgp_not_supported(p, e),
           `_print_Wild` = function(p, e) .cgp_not_supported(p, e),
           `_print_WildFunction` = function(p, e) .cgp_not_supported(p, e),
           `_print_Relational` = function(p, e) .cgp_not_supported(p, e),
           `_rate_index_position` = function(p, pos) .cg_err("NotImplementedError", "This function must be implemented by subclass of CodePrinter."),
           `_get_statement` = function(p, s) .cg_err("NotImplementedError", "This function must be implemented by subclass of CodePrinter."),
           `_get_comment` = function(p, s) .cg_err("NotImplementedError", "This function must be implemented by subclass of CodePrinter."),
           `_declare_number_const` = function(p, name, value) .cg_err("NotImplementedError", "This function must be implemented by subclass of CodePrinter."),
           `_format_code` = function(p, lines) .cg_err("NotImplementedError", "This function must be implemented by subclass of CodePrinter."),
           `_get_loop_opening_ending` = function(p, indices) .cg_err("NotImplementedError", "This function must be implemented by subclass of CodePrinter.")
         ))

# does the tree contain a node of this type
.cg_has_type <- function(e, type) {
  if (is.function(e)) return(type == "IndexedBase" && !is.null(attr(e, "base")))
  if (!is_sym(e)) return(FALSE)
  if (e$type == type) return(TRUE)
  for (a in e$args) if (.cg_has_type(a, type)) return(TRUE)
  FALSE
}

# the text of a String node, or a plain string
.cg_text <- function(s) if (is.character(s)) s else s$name

# ---------------------------------------------------------------------------
# Loops over Indexed objects (CodePrinter._doprint_loops)

.cgp_idx_lower <- function(i) i$extra$range[[1L]]
.cgp_idx_upper <- function(i) i$extra$range[[2L]]

# AugmentedAssignment.op
.cg_aug_op <- function(e) paste0(.AUG_OPS[[e$type]], "=")

.cgp_doprint_loops <- function(p, expr, assign_to = NULL) {
  if (isTRUE(.cgp_s(p, "contract", TRUE))) {
    indices <- .cgp_expression_indices(p, expr, assign_to)
    dummies <- get_contraction_structure(expr)
  } else {
    indices <- list()
    dummies <- list(list(key = NULL, value = list(expr)))
  }
  oc <- .cgp_call(p, "_get_loop_opening_ending", indices)
  none_terms <- NULL
  for (d in dummies) if (is.null(d$key)) none_terms <- d$value
  text <- if (!is.null(none_terms)) .cgp_strprint(p, .add(none_terms)) else .cgp_strprint(p, S_Zero)
  lhs_printed <- .cgpp(p, assign_to)
  lines <- character()
  if (text != lhs_printed) {
    lines <- c(lines, oc[[1L]])
    if (!is.null(assign_to)) text <- .cgp_call(p, "_get_statement", sprintf("%s = %s", lhs_printed, text))
    lines <- c(lines, text, oc[[2L]])
  }
  for (d in dummies) {
    if (is.null(d$key) || is_sym(d$key)) next
    ind <- .cgp_sort_optimized(p, d$key, expr)
    ocd <- .cgp_call(p, "_get_loop_opening_ending", ind)
    for (term in d$value) {
      for (d2 in dummies) if (is_sym(d2$key) && identical(d2$key$key, term$key))
        .cg_err("NotImplementedError", "FIXME: no support for contractions in factor yet")
      if (is.null(assign_to)) .cg_err("AssignmentError", "need assignment variable for loops")
      if (.cg_contains(term, assign_to))
        .cg_err("ValueError", "FIXME: lhs present in rhs,                                this is undefined in CodePrinter")
      lines <- c(lines, oc[[1L]], ocd[[1L]])
      txt <- sprintf("%s = %s", lhs_printed, .cgp_strprint(p, .add(list(assign_to, term))))
      lines <- c(lines, .cgp_call(p, "_get_statement", txt), ocd[[2L]], oc[[2L]])
    }
  }
  paste(lines, collapse = "\n")
}

# StrPrinter.doprint(self, expr): the printer's own _print at a fresh level
.cgp_strprint <- function(p, e) .cgpp(p, e)

.cg_contains <- function(e, sub) {
  if (!is_sym(e)) return(FALSE)
  if (identical(e$key, sub$key)) return(TRUE)
  for (a in e$args) if (.cg_contains(a, sub)) return(TRUE)
  FALSE
}

.cgp_expression_indices <- function(p, expr, assign_to) {
  rinds <- get_indices(expr)$indices
  linds <- get_indices(assign_to)$indices
  if (length(linds) && !length(rinds)) rinds <- linds
  rk <- sort(vapply(rinds, .k, "")); lk <- sort(vapply(linds, .k, ""))
  if (!identical(rk, lk)) .cg_err("ValueError", sprintf("lhs indices must match non-dummy rhs indices in %s", sstr(expr)))
  .cgp_sort_optimized(p, rinds, assign_to)
}

.cgp_sort_optimized <- function(p, indices, expr) {
  if (!length(indices)) return(list())
  score <- setNames(rep(0, length(indices)), vapply(indices, .k, ""))
  arrays <- .cg_collect(expr, "Indexed")
  for (arr in arrays) {
    inds <- arr$args[-1L]
    for (pos in seq_along(inds)) {
      k <- inds[[pos]]$key
      if (k %in% names(score)) score[[k]] <- score[[k]] + .cgp_call(p, "_rate_index_position", pos - 1L)
    }
  }
  indices[order(score[vapply(indices, .k, "")])]
}

# the distinct subtrees of a type, in the order of first appearance
.cg_collect <- function(e, type) {
  out <- list()
  walk <- function(z) {
    if (!is_sym(z)) return(invisible())
    if (z$type == type && !(z$key %in% names(out))) out[[z$key]] <<- z
    for (a in z$args) walk(a)
  }
  walk(e)
  unname(out)
}

# ---------------------------------------------------------------------------
# Public interface

#' Code printer objects
#'
#' Printer objects after SymPy's code printer classes. `C89CodePrinter()`,
#' `C99CodePrinter()` and `C11CodePrinter()` print C; `CXX98CodePrinter()`,
#' `CXX11CodePrinter()` and `CXX17CodePrinter()` print C++;
#' `FCodePrinter()` prints Fortran; `PythonCodePrinter()`,
#' `MpmathPrinter()`, `NumPyPrinter()` and `SciPyPrinter()` print Python;
#' `JuliaCodePrinter()`, `OctaveCodePrinter()`, `RustCodePrinter()` and
#' `JavascriptCodePrinter()` print those languages. `printer_doprint(p, e)`
#' prints an expression with a printer, as the `doprint` method does; the
#' printer keeps what the printed code needs (`headers`, `libraries`,
#' `macros`, `module_imports`), which [ast_attr()] reads.
#'
#' @param settings A named list of printer settings, such as
#'   `list(type_aliases = ...)` or `list(standard = 95)`.
#' @param p A printer.
#' @param e An expression or codegen node.
#' @param assign_to Optional variable to assign to.
#' @return A printer, or for `printer_doprint()` a character scalar.
#' @examples
#' x <- Symbol("x")
#' p <- C99CodePrinter()
#' printer_doprint(p, sym_expm1(x))
#' printer_doprint(FCodePrinter(list(source_format = "free")), x^2)
#' @name code-printer-objects
NULL

#' @rdname code-printer-objects
#' @export
printer_doprint <- function(p, e, assign_to = NULL) {
  f <- .cgp_get(p, "doprint")
  if (length(formals(f)) >= 3L) f(p, e, assign_to = assign_to) else f(p, e)
}

# The printers of R/printing-code.R that run on this engine
.cg_engine <- function(style) switch(style, c = .cg_ccode, fortran = .cg_fcode, cxx = .cg_cxxcode, python = .cg_pycode,
                                         julia = .cg_julia_code, octave = .cg_octave_code, js = .cg_jscode, r = .cg_rcode, NULL)

# Idx.label: the integer symbol an Idx is named after
.cg_idx_label <- function(i) Symbol(i$name, integer = TRUE)

# The inverse reciprocal functions rewritten for printers that lack them
# (asec._eval_rewrite_as_acos and its siblings)
.cg_add_rewrite("acos", "asec", function(a) fcall("acos", list(.inv(a))))
.cg_add_rewrite("asin", "acsc", function(a) fcall("asin", list(.inv(a))))
.cg_add_rewrite("atan", "acot", function(a) fcall("atan", list(.inv(a))))

# ITE._eval_rewrite_as_Piecewise
.cg_ite_piecewise <- function(e) Piecewise(list(e$args[[2L]], e$args[[1L]]), list(e$args[[3L]], TRUE))
# sinc._eval_rewrite_as_sin
.cg_add_rewrite("sin", "sinc", function(a) Piecewise(list(.div(fcall("sin", list(a)), a), Ne(a, S_Zero)), list(S_One, TRUE)))
