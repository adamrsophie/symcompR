# Abstract syntax tree for code generation, after sympy/codegen/ast.py.
#
# Every node is a sym whose type is the SymPy class name. The fields of a
# Token are stored in order as its args (a String keeps its text in the
# name slot, as SymPy keeps it out of args), so subs, xreplace and the
# generic printers work on these nodes. The class table below records the
# fields, defaults and attribute constructors of each class, following the
# class attributes of the SymPy source; subclasses inherit them.

.cg_classes <- new.env(hash = TRUE)
.cg_ctx <- new.env()
.cg_ctx$il <- 0L
.srepr_printers <- new.env()

.cg_err <- function(cls, msg) {
  stop(structure(class = c(cls, "error", "condition"), list(message = msg, call = NULL)))
}

# Python repr of a str: single quotes unless the text holds a single quote
# and no double quote.
.py_repr_str <- function(s) {
  q <- if (grepl("'", s, fixed = TRUE) && !grepl("\"", s, fixed = TRUE)) "\"" else "'"
  s <- gsub("\\", "\\\\", s, fixed = TRUE)
  if (q == "'") s <- gsub("'", "\\'", s, fixed = TRUE)
  s <- gsub("\n", "\\n", s, fixed = TRUE)
  s <- gsub("\t", "\\t", s, fixed = TRUE)
  paste0(q, s, q)
}

.cg_tuple <- function(args = list()) new_sym("Tuple", unname(args))

# Register a class. Unset properties are inherited from `base`; `construct`
# merges with the base constructors (the _construct_<attr> methods).
.cg_def <- function(name, base = NULL, fields = NULL, defaults = NULL, not_in_args = NULL,
                    indented_args = NULL, construct = list(), token = NULL, ...) {
  b <- if (is.null(base)) NULL else .cg_classes[[base]]
  extra <- list(...)
  pick <- function(v, nm, dflt) if (!is.null(v)) v else if (!is.null(b)) b[[nm]] else dflt
  cl <- list(
    name = name,
    isa = c(name, if (!is.null(b)) b$isa),
    fields = pick(fields, "fields", character()),
    defaults = pick(defaults, "defaults", list()),
    not_in_args = pick(not_in_args, "not_in_args", character()),
    indented_args = pick(indented_args, "indented_args", "body"),
    construct = utils::modifyList(if (!is.null(b)) b$construct else list(), construct),
    token = pick(token, "token", TRUE))
  if (!is.null(b)) for (nm in setdiff(names(b), names(cl))) cl[[nm]] <- b[[nm]]
  for (nm in names(extra)) cl[[nm]] <- extra[[nm]]
  assign(name, cl, envir = .cg_classes)
  register_str_printer(name, function(e) .cg_print(e, FALSE))
  assign(name, function(e) .cg_print(e, TRUE), envir = .srepr_printers)
  invisible(cl)
}

.cg_cls <- function(e) if (is_sym(e)) .cg_classes[[e$type]] else NULL
.cg_isa <- function(e, cls) {
  if (!is_sym(e)) return(FALSE)
  cl <- .cg_classes[[e$type]]
  !is.null(cl) && cls %in% cl$isa
}
.is_none <- function(v) is.null(v) || (is_sym(v) && v$type == "NoneToken")

# Field values of a node as a named list.
.cg_vals <- function(e) {
  cl <- .cg_cls(e)
  out <- list()
  k <- 0L
  for (f in cl$fields) {
    if (f %in% cl$not_in_args) out[f] <- list(e$name)
    else { k <- k + 1L; out[f] <- list(e$args[[k]]) }
  }
  out
}

.cg_eq <- function(a, b) {
  if (.is_none(a) || .is_none(b)) return(.is_none(a) && .is_none(b))
  if (is_sym(a) && is_sym(b)) return(identical(a$key, b$key))
  identical(a, b)
}

# Token._construct
.cg_construct <- function(cl, attr, arg) {
  if (.is_none(arg)) {
    d <- cl$defaults[[attr]]
    return(if (is.null(d)) none else d)
  }
  if (is_sym(arg) && arg$type == "Dummy") return(arg)
  f <- cl$construct[[attr]]
  if (is.null(f)) arg else f(arg)
}

.cg_make <- function(cl, vals) {
  args <- list()
  nm <- NULL
  for (f in cl$fields) {
    if (f %in% cl$not_in_args) nm <- vals[[f]]
    else args[[length(args) + 1L]] <- .cg_as_arg(vals[[f]])
  }
  new_sym(cl$name, args, name = nm)
}

# A field value that is not a sym (a raw R value passed to a field without
# a constructor) is sympified, as Basic.__new__ would store it.
.cg_as_arg <- function(v) {
  if (is_sym(v)) return(v)
  if (inherits(v, "symmatrix")) return(.as_block(v))
  if (is.character(v) && length(v) == 1L) return(String(v))
  as_sym(v)
}

# Token.__new__
.cg_new <- function(name, dots) {
  cl <- .cg_classes[[name]]
  nms <- names(dots)
  if (is.null(nms)) nms <- rep("", length(dots))
  pos <- dots[nms == ""]
  kw <- dots[nms != ""]
  if (length(pos) == 1L && !length(kw) && .cg_isa(pos[[1L]], name)) return(pos[[1L]])
  fields <- cl$fields
  if (length(pos) > length(fields))
    .cg_err("ValueError", sprintf("Too many arguments (%d), expected at most %d", length(pos), length(fields)))
  vals <- list()
  for (i in seq_along(pos)) {
    a <- fields[i]
    if (a %in% names(kw)) .cg_err("TypeError", sprintf("Got multiple values for attribute %s", .py_repr_str(a)))
    vals[a] <- list(.cg_construct(cl, a, pos[[i]]))
  }
  for (a in fields[seq_along(fields) > length(pos)]) {
    if (a %in% names(kw)) {
      v <- kw[[a]]
      kw[[a]] <- NULL
    } else if (a %in% names(cl$defaults)) {
      v <- cl$defaults[[a]]
    } else {
      .cg_err("TypeError", sprintf("No value for %s given and attribute has no default", .py_repr_str(a)))
    }
    vals[a] <- list(.cg_construct(cl, a, v))
  }
  if (length(kw)) .cg_err("ValueError", paste0("Unknown keyword arguments: ", paste(names(kw), collapse = " ")))
  if (!is.null(cl$post)) return(cl$post(.cg_make(cl, vals)))
  .cg_make(cl, vals)
}

# ---------------------------------------------------------------------------
# printing: Token._sympystr and Token._sympyrepr

.cg_p <- function(v, repr) {
  if (is.character(v)) return(if (repr) .py_repr_str(v) else v)
  if (inherits(v, "symmatrix")) v <- .as_block(v)
  if (repr) srepr(v) else .str(v)
}

.cg_with_il <- function(il, expr) {
  old <- .cg_ctx$il
  .cg_ctx$il <- il
  on.exit(.cg_ctx$il <- old)
  expr
}

.cg_print <- function(e, repr, joiner = ", ") {
  cl <- .cg_cls(e)
  f <- if (repr) cl$repr else cl$str
  if (!is.null(f)) return(f(e, repr = repr, joiner = joiner))
  if (!cl$token) return(.cg_basic_print(e, repr))
  .cg_token_print(e, repr, joiner)
}

.cg_basic_print <- function(e, repr) {
  paste0(e$type, "(", paste(vapply(e$args, .cg_p, "", repr = repr), collapse = ", "), ")")
}

.cg_joiner <- function(cl, k, il) if (k %in% cl$indented_args) paste0(",\n", strrep(" ", il)) else ", "

.cg_indented <- function(cl, k, v, repr) {
  il <- .cg_ctx$il
  pr <- function(arg) {
    if (is_sym(arg) && isTRUE(.cg_cls(arg)$token)) .cg_print(arg, repr, joiner = .cg_joiner(cl, k, il))
    else .cg_p(arg, repr)
  }
  if (is_sym(v) && v$type == "Tuple") {
    joined <- paste(vapply(v$args, pr, ""), collapse = .cg_joiner(cl, k, il))
    if (k %in% cl$indented_args)
      return(paste0("(\n", strrep(" ", il), joined, ",\n", strrep(" ", max(il - 4L, 0L)), ")"))
    return(if (length(v$args) == 1L) paste0("(", joined, ",)") else paste0("(", joined, ")"))
  }
  pr(v)
}

.cg_token_print <- function(e, repr, joiner = ", ", exclude = character()) {
  cl <- .cg_cls(e)
  vals <- .cg_vals(e)
  il <- .cg_ctx$il
  parts <- character()
  for (i in seq_along(cl$fields)) {
    a <- cl$fields[i]
    if (a %in% exclude) next
    v <- vals[[a]]
    if (a %in% names(cl$defaults) && .cg_eq(v, cl$defaults[[a]])) next
    ilvl <- if (a %in% cl$indented_args) il + 4L else 0L
    s <- .cg_with_il(ilvl, .cg_indented(cl, a, v, repr))
    s <- sub("^\\s+", "", s, perl = TRUE)
    parts <- c(parts, if (i == 1L) s else paste0(a, "=", s))
  }
  paste0(cl$name, "(", paste(parts, collapse = joiner), ")")
}

# CodeBlock._sympystr and _sympyrepr
.cg_codeblock_print <- function(e, repr, ...) {
  il <- .cg_ctx$il
  joined <- paste(vapply(e$args, .cg_p, "", repr = repr), collapse = paste0(",\n", strrep(" ", il)))
  paste0(strrep(" ", max(il - 4L, 0L)), "CodeBlock(\n", strrep(" ", il), joined, "\n",
         strrep(" ", max(il - 4L, 0L)), ")")
}

# ---------------------------------------------------------------------------
# attribute constructors

.cg_mk_tuple <- function(args) {
  if (is_sym(args) && args$type == "Tuple") args <- args$args
  .cg_tuple(lapply(.cg_iter(args), function(a) if (is.character(a)) String(a) else .cg_as_arg(a)))
}

# The elements of an iterable argument: an R list or vector, or a sym Tuple.
.cg_iter <- function(x) {
  if (is_sym(x)) {
    if (x$type %in% c("Tuple", "List", "CodeBlock")) return(x$args)
    .cg_err("TypeError", paste0("'", x$type, "' object is not iterable"))
  }
  if (inherits(x, "symmatrix")) return(x$data)
  if (is.list(x)) return(unname(x))
  if (is.atomic(x) && !is.character(x)) return(as.list(x))
  if (is.character(x)) return(as.list(strsplit(x, "")[[1L]]))
  .cg_err("TypeError", "object is not iterable")
}
.cg_is_iterable <- function(x) {
  if (is_sym(x)) return(x$type %in% c("Tuple", "List", "Range", "FiniteSet", "Interval", "CodeBlock", "Dict") ||
                          .is_matexpr(x) || x$type == "ImmutableMatrix")
  inherits(x, "symmatrix") || is.list(x) || (is.atomic(x) && length(x) != 1L)
}

.cg_sympify <- function(x) {
  if (is_sym(x)) return(x)
  if (inherits(x, "symmatrix")) return(.as_block(x))
  if (is.function(x) && !is.null(attr(x, "base"))) return(attr(x, "base"))
  if (is.list(x)) return(.cg_tuple(lapply(x, .cg_sympify)))
  as_sym(x)
}

# a comparison that is True: Python's bool() of a decided relational
.cg_true <- function(x) isTRUE(x) || .is_true(x)

# _sympify: strings are rejected
.cg_strict_sympify <- function(x) {
  if (is.character(x)) stop(SympifyError(paste0("SympifyError: ", .py_repr_str(x))))
  .cg_sympify(x)
}

.cg_body <- function(itr) if (is_sym(itr) && itr$type == "CodeBlock") itr else do.call(CodeBlock, .cg_iter(itr))

# ---------------------------------------------------------------------------
# Token and its basic subclasses

.cg_def("Token")
.cg_def("BreakToken", "Token")
.cg_def("ContinueToken", "Token")
.cg_def("NoneToken", "Token")

#' Code generation abstract syntax tree
#'
#' Nodes of an abstract syntax tree for code generation, after SymPy's
#' `sympy.codegen.ast`. Each constructor takes its fields positionally or by
#' name, as the SymPy classes do: for example `Variable(x, type = real)`,
#' `Assignment(x, y + 1)`, `For(i, Range(10), list(aug_assign(s, "+", i)))`
#' or `FunctionDefinition(real, "f", list(x), list(Return(x^2)))`. The nodes
#' print through the code printers such as [ccode()] and [fcode()].
#' `aug_assign(lhs, op, rhs)` builds the augmented assignment for the
#' operator `op`. The objects `real`, `sym_integer`, `complex_`, `bool_`,
#' `intc`, `float32`, `float64` and the rest are the predefined types;
#' `value_const` and `pointer_const` are attributes; `none` is the `None`
#' token.
#'
#' @param ... The fields of the node, positionally or by name.
#' @param lhs,rhs Left and right hand sides of an assignment.
#' @param op An operator: `"+"`, `"-"`, `"*"`, `"/"` or `"%"`.
#' @return A `sym` node.
#' @examples
#' x <- Symbol("x"); y <- Symbol("y")
#' Assignment(x, y + 1)
#' Variable(x, type = real)
#' CodeBlock(Assignment(x, 1), aug_assign(x, "+", y))
#' @name codegen-ast
NULL

#' @rdname codegen-ast
#' @export
Token <- function(...) .cg_new("Token", list(...))
#' @rdname codegen-ast
#' @export
BreakToken <- function(...) .cg_new("BreakToken", list(...))
#' @rdname codegen-ast
#' @export
ContinueToken <- function(...) .cg_new("ContinueToken", list(...))
#' @rdname codegen-ast
#' @export
NoneToken <- function(...) .cg_new("NoneToken", list(...))
#' @rdname codegen-ast
#' @export
none <- new_sym("NoneToken")
#' @rdname codegen-ast
#' @export
break_ <- new_sym("BreakToken")
#' @rdname codegen-ast
#' @export
continue_ <- new_sym("ContinueToken")

.cg_def("String", "Token", fields = "text", not_in_args = "text",
        construct = list(text = function(text) {
          if (!is.character(text) || length(text) != 1L) .cg_err("TypeError", "Argument text is not a string type.")
          text
        }),
        str = function(e, ...) e$name)
.cg_def("QuotedString", "String")
.cg_def("Comment", "String")

#' @rdname codegen-ast
#' @export
String <- function(...) .cg_new("String", list(...))
#' @rdname codegen-ast
#' @export
QuotedString <- function(...) .cg_new("QuotedString", list(...))
#' @rdname codegen-ast
#' @export
Comment <- function(...) .cg_new("Comment", list(...))

.cg_def("Node", "Token", fields = "attrs", defaults = list(attrs = new_sym("Tuple")),
        construct = list(attrs = function(a) .cg_mk_tuple(a)))
#' @rdname codegen-ast
#' @export
Node <- function(...) .cg_new("Node", list(...))

# ---------------------------------------------------------------------------
# Types

.cg_def("Type", "Token", fields = "name", construct = list(name = function(a) String(a)),
        str = function(e, ...) .str(e$args[[1L]]), repr = NULL)
.cg_def("IntBaseType", "Type")
.cg_def("_SizedIntType", "IntBaseType", fields = c("name", "nbits"),
        construct = list(nbits = function(a) Integer(a)))
.cg_def("SignedIntType", "_SizedIntType")
.cg_def("UnsignedIntType", "_SizedIntType")
.cg_def("FloatBaseType", "Type")
.cg_def("FloatType", "FloatBaseType", fields = c("name", "nbits", "nmant", "nexp"),
        construct = list(nbits = function(a) Integer(a), nmant = function(a) Integer(a), nexp = function(a) Integer(a)))
.cg_def("ComplexBaseType", "FloatBaseType")
.cg_def("ComplexType", "FloatType")
assign("ComplexType", within(.cg_classes$ComplexType, isa <- c("ComplexType", "ComplexBaseType", "FloatType", "FloatBaseType", "Type", "Token")),
       envir = .cg_classes)

#' @rdname codegen-ast
#' @export
Type <- function(...) .cg_new("Type", list(...))
#' @rdname codegen-ast
#' @export
IntBaseType <- function(...) .cg_new("IntBaseType", list(...))
#' @rdname codegen-ast
#' @export
SignedIntType <- function(...) .cg_new("SignedIntType", list(...))
#' @rdname codegen-ast
#' @export
UnsignedIntType <- function(...) .cg_new("UnsignedIntType", list(...))
#' @rdname codegen-ast
#' @export
FloatBaseType <- function(...) .cg_new("FloatBaseType", list(...))
#' @rdname codegen-ast
#' @export
FloatType <- function(...) .cg_new("FloatType", list(...))
#' @rdname codegen-ast
#' @export
ComplexBaseType <- function(...) .cg_new("ComplexBaseType", list(...))
#' @rdname codegen-ast
#' @export
ComplexType <- function(...) .cg_new("ComplexType", list(...))

#' @rdname codegen-ast
#' @export
intc <- IntBaseType("intc")
#' @rdname codegen-ast
#' @export
intp <- IntBaseType("intp")
#' @rdname codegen-ast
#' @export
int8 <- SignedIntType("int8", 8)
#' @rdname codegen-ast
#' @export
int16 <- SignedIntType("int16", 16)
#' @rdname codegen-ast
#' @export
int32 <- SignedIntType("int32", 32)
#' @rdname codegen-ast
#' @export
int64 <- SignedIntType("int64", 64)
#' @rdname codegen-ast
#' @export
uint8 <- UnsignedIntType("uint8", 8)
#' @rdname codegen-ast
#' @export
uint16 <- UnsignedIntType("uint16", 16)
#' @rdname codegen-ast
#' @export
uint32 <- UnsignedIntType("uint32", 32)
#' @rdname codegen-ast
#' @export
uint64 <- UnsignedIntType("uint64", 64)
#' @rdname codegen-ast
#' @export
float16 <- FloatType("float16", 16, nexp = 5, nmant = 10)
#' @rdname codegen-ast
#' @export
float32 <- FloatType("float32", 32, nexp = 8, nmant = 23)
#' @rdname codegen-ast
#' @export
float64 <- FloatType("float64", 64, nexp = 11, nmant = 52)
#' @rdname codegen-ast
#' @export
float80 <- FloatType("float80", 80, nexp = 15, nmant = 63)
#' @rdname codegen-ast
#' @export
float128 <- FloatType("float128", 128, nexp = 15, nmant = 112)
#' @rdname codegen-ast
#' @export
float256 <- FloatType("float256", 256, nexp = 19, nmant = 236)
#' @rdname codegen-ast
#' @export
complex64 <- ComplexType("complex64", nbits = 64, nmant = 23, nexp = 8)
#' @rdname codegen-ast
#' @export
complex128 <- ComplexType("complex128", nbits = 128, nmant = 52, nexp = 11)
#' @rdname codegen-ast
#' @export
untyped <- Type("untyped")
#' @rdname codegen-ast
#' @export
real <- FloatBaseType("real")
# package code refers to the type by this name; base R has a defunct real()
.cg_real <- real
#' @rdname codegen-ast
#' @export
sym_integer <- IntBaseType("integer")
#' @rdname codegen-ast
#' @export
complex_ <- ComplexBaseType("complex")
#' @rdname codegen-ast
#' @export
bool_ <- Type("bool")

# FloatType derived attributes, from the bit counts
.ft_num <- function(e, nm) as.numeric(gmp::numerator(e$args[[which(.cg_classes$FloatType$fields == nm)]]$value))
.ft_attr <- function(e, a) {
  nmant <- .ft_num(e, "nmant"); nexp <- .ft_num(e, "nexp")
  two <- function(p) make_rational(if (p >= 0) gmp::as.bigq(gmp::pow.bigz(2, p)) else gmp::as.bigq(1, gmp::pow.bigz(2, -p)))
  maxe <- 2^(nexp - 1)
  switch(a,
         max_exponent = two(nexp - 1),
         min_exponent = Integer(3 - maxe),
         max = make_rational((1 - gmp::as.bigq(1, gmp::pow.bigz(2, nmant + 1))) * gmp::as.bigq(gmp::pow.bigz(2, maxe))),
         tiny = two(3 - maxe - 1),
         eps = two(-nmant),
         dig = Integer(floor(nmant * log10(2))),
         decimal_dig = Integer(ceiling((nmant + 1) * log10(2) + 1)),
         NULL)
}

#' Codegen type helpers
#'
#' `type_from_expr(e)` deduces a codegen type from an expression, as SymPy's
#' `Type.from_expr`. `type_cast_check(t, value)` casts `value` to type `t`
#' and checks that the value did not change by more than the tolerance and
#' that it lies in the type's range; `type_cast_nocheck(t, value)` casts
#' without the checks. `variable_deduced(symbol, value)` builds a
#' [Variable()] whose type is deduced from the symbol or the value, and
#' `as_Declaration(v, ...)` turns a variable into a [Declaration()] with
#' some fields replaced.
#'
#' @param e,value An expression or number.
#' @param t A codegen type such as `real` or `float32`.
#' @param rtol,atol Relative and absolute tolerances.
#' @param symbol A symbol.
#' @param attrs Attributes of the variable.
#' @param cast_check Check the cast of `value`.
#' @param v A [Variable()].
#' @param ... Fields to replace.
#' @return A type, a number, or a node.
#' @examples
#' type_from_expr(Symbol("n", integer = TRUE))
#' type_cast_check(int8, 127)
#' variable_deduced(Symbol("x", real = TRUE))
#' @name codegen-types
NULL

#' @rdname codegen-types
#' @export
type_from_expr <- function(e) {
  if (is.double(e) && length(e) == 1L && !is.integer(e)) return(.cg_real)
  if (is.integer(e) && length(e) == 1L) return(sym_integer)
  if (is.complex(e) && length(e) == 1L) return(complex_)
  if (is.logical(e) && length(e) == 1L) return(bool_)
  if (is_sym(e)) {
    if (e$type == "Float") return(.cg_real)
    if (e$type == "Integer") return(sym_integer)
    f <- tryCatch(.facts(e), error = function(err) NULL)
    if (!is.null(f)) {
      if (isTRUE(f[["integer"]])) return(sym_integer)
      if (isTRUE(f[["real"]])) return(.cg_real)
      if (isTRUE(f[["complex"]])) return(complex_)
    }
    if (e$type %in% .REL_TYPES) return(bool_)
  }
  .cg_err("ValueError", "Could not deduce type from expr.")
}

.cg_type_kind <- function(t) {
  isa <- .cg_classes[[t$type]]$isa
  if ("ComplexType" %in% isa) "ComplexType"
  else if ("ComplexBaseType" %in% isa) "ComplexBaseType"
  else if ("FloatType" %in% isa) "FloatType"
  else if ("FloatBaseType" %in% isa) "FloatBaseType"
  else if ("_SizedIntType" %in% isa) "SizedInt"
  else if ("IntBaseType" %in% isa) "IntBaseType"
  else "Type"
}

.cg_float_nocheck <- function(t, value) {
  v <- as_sym(value)
  if (is_const(v, "oo")) return(Inf)
  if (is_const(v, "-oo")) return(-Inf)
  if (is_nan(v)) return(v)
  dd <- .int_val(.ft_attr(t, "decimal_dig"))
  Float(.str_top(N(v, dd)), dd)
}

#' @rdname codegen-types
#' @export
type_cast_nocheck <- function(t, value) {
  switch(.cg_type_kind(t),
         IntBaseType = , SizedInt = {
           v <- if (is.character(value)) as_sym(value) else as_sym(value)
           if (v$type == "Float") Integer(trunc(v$value)) else if (v$type == "Rational") Integer(gmp::numerator(v$value) %/% gmp::denominator(v$value)) else Integer(v)
         },
         FloatBaseType = Float(value),
         FloatType = .cg_float_nocheck(t, value),
         ComplexBaseType = .add(list(Float(re(as_sym(value))), .mul(list(Float(im(as_sym(value))), S_I)))),
         ComplexType = {
           r <- .cg_float_nocheck(t, re(as_sym(value)))
           i <- .cg_float_nocheck(t, im(as_sym(value)))
           .add(list(as_sym(r), .mul(list(as_sym(i), S_I))))
         },
         .cg_err("AttributeError", paste0("'", t$type, "' object has no attribute 'cast_nocheck'")))
}

.cg_type_check <- function(t, value) {
  kind <- .cg_type_kind(t)
  fmt_d <- function(z) { z <- as_sym(z); as.character(if (z$type == "Float") trunc(z$value) else gmp::numerator(z$value) %/% gmp::denominator(z$value)) }
  float_check <- function(v) {
    mx <- .ft_attr(t, "max")
    if (.cg_true(as_sym(v) < .neg(mx))) .cg_err("ValueError", sprintf("Value is too small: %s < %s", fmt_d(v), fmt_d(.neg(mx))))
    if (.cg_true(as_sym(v) > mx)) .cg_err("ValueError", sprintf("Value is too big: %s > %s", fmt_d(v), fmt_d(mx)))
    if (.cg_true(Abs(as_sym(v)) < .ft_attr(t, "tiny")))
      .cg_err("ValueError", "Smallest (absolute) value for data type bigger than new value.")
  }
  if (kind == "SizedInt") {
    nb <- .int_val(t$args[[2L]])
    signed <- "SignedIntType" %in% .cg_classes[[t$type]]$isa
    mn <- if (signed) -gmp::pow.bigz(2, nb - 1) else gmp::as.bigz(0)
    mx <- if (signed) gmp::pow.bigz(2, nb - 1) - 1 else gmp::pow.bigz(2, nb) - 1
    vz <- gmp::numerator(value$value)
    if (vz < mn) .cg_err("ValueError", sprintf("Value is too small: %s < %s", as.character(vz), as.character(mn)))
    if (vz > mx) .cg_err("ValueError", sprintf("Value is too big: %s > %s", as.character(vz), as.character(mx)))
  } else if (kind == "FloatType") {
    float_check(value)
  } else if (kind == "ComplexType") {
    float_check(re(value)); float_check(im(value))
  }
  invisible(NULL)
}

#' @rdname codegen-types
#' @export
type_cast_check <- function(t, value, rtol = NULL, atol = 0) {
  val <- if (is.character(value)) sympify(value) else as_sym(value)
  dd <- if (.cg_type_kind(t) %in% c("FloatType", "ComplexType")) .ft_attr(t, "decimal_dig") else NULL
  if (is.null(rtol)) rtol <- if (is.null(dd)) Float(1e-15) else .mul(list(Float(2), .pow(Integer(10), .neg(dd))))
  tol <- function(num) .add(list(as_sym(atol), .mul(list(as_sym(rtol), Abs(num)))))
  new_val <- type_cast_nocheck(t, value)
  .cg_type_check(t, as_sym(new_val))
  delta <- .sub(as_sym(new_val), val)
  if (.cg_true(Abs(delta) > tol(val))) .cg_err("ValueError", "Casting gives a significantly different value.")
  new_val
}

# ---------------------------------------------------------------------------
# Attribute, Variable, Pointer, Element, Declaration

.cg_def("Attribute", "Token", fields = c("name", "parameters"), defaults = list(parameters = new_sym("Tuple")),
        construct = list(name = function(a) String(a), parameters = function(a) .cg_mk_tuple(a)),
        str = function(e, repr, ...) {
          r <- .str(e$args[[1L]])
          ps <- e$args[[2L]]$args
          if (length(ps)) r <- paste0(r, "(", paste(vapply(ps, .cg_p, "", repr = repr), collapse = ", "), ")")
          r
        })
#' @rdname codegen-ast
#' @export
Attribute <- function(...) .cg_new("Attribute", list(...))
#' @rdname codegen-ast
#' @export
value_const <- Attribute("value_const")
#' @rdname codegen-ast
#' @export
pointer_const <- Attribute("pointer_const")

.cg_def("Variable", "Node", fields = c("symbol", "type", "value", "attrs"),
        defaults = list(attrs = new_sym("Tuple"), type = untyped, value = none),
        construct = list(symbol = function(a) .cg_sympify(a), value = function(a) .cg_sympify(a)))
.cg_def("Pointer", "Variable")
#' @rdname codegen-ast
#' @export
Variable <- function(...) .cg_new("Variable", list(...))
#' @rdname codegen-ast
#' @export
Pointer <- function(...) .cg_new("Pointer", list(...))

.cg_deduced <- function(cls, symbol, value = NULL, attrs = new_sym("Tuple"), cast_check = TRUE) {
  if (.cg_isa(symbol, "Variable")) return(symbol)
  type_ <- tryCatch(type_from_expr(if (is.character(symbol)) sympify(symbol) else symbol),
                    ValueError = function(err) type_from_expr(value))
  if (!is.null(value) && cast_check) value <- type_cast_check(type_, value)
  .cg_new(cls, list(symbol, type = type_, value = value, attrs = attrs))
}

#' @rdname codegen-types
#' @export
variable_deduced <- function(symbol, value = NULL, attrs = list(), cast_check = TRUE)
  .cg_deduced("Variable", symbol, value, attrs, cast_check)

#' @rdname codegen-types
#' @export
pointer_deduced <- function(symbol, value = NULL, attrs = list(), cast_check = TRUE)
  .cg_deduced("Pointer", symbol, value, attrs, cast_check)

#' @rdname codegen-types
#' @export
as_Declaration <- function(v, ...) {
  kw <- ast_kwargs(v)
  upd <- list(...)
  for (nm in names(upd)) kw[nm] <- list(upd[[nm]])
  Declaration(do.call(ast_func(v), kw))
}

.cg_def("Element", "Token", fields = c("symbol", "indices", "strides", "offset"),
        defaults = list(strides = none, offset = none),
        construct = list(symbol = function(a) .cg_sympify(a),
                         indices = function(a) .cg_tuple(lapply(.cg_iter(a), .cg_sympify)),
                         strides = function(a) .cg_tuple(lapply(.cg_iter(a), .cg_sympify)),
                         offset = function(a) .cg_sympify(a)))
#' @rdname codegen-ast
#' @export
Element <- function(...) .cg_new("Element", list(...))

.cg_def("Declaration", "Token", fields = "variable", construct = list(variable = function(a) Variable(a)))
#' @rdname codegen-ast
#' @export
Declaration <- function(...) .cg_new("Declaration", list(...))

# ---------------------------------------------------------------------------
# Assignments

.cg_assignable <- function(lhs) {
  if (!is_sym(lhs)) return(FALSE)
  if (lhs$type %in% c("Symbol", "Dummy")) return(!.is_idx(lhs))
  lhs$type %in% c("MatrixSymbol", "MatrixElement", "Indexed", "ArrayElement") ||
    .cg_isa(lhs, "Element") || .cg_isa(lhs, "Variable")
}
.cg_shape <- function(e) {
  if (inherits(e, "symmatrix")) return(list(.int(e$nrow), .int(e$ncol)))
  if (!is_sym(e)) return(NULL)
  if (e$type == "Indexed") return(NULL)
  if (.is_matexpr(e)) return(.mshape(e))
  if (e$type == "MatrixSolve") return(.cg_shape(e$args[[2L]]))
  NULL
}

.cg_check_assign <- function(lhs, rhs) {
  if (!.cg_assignable(lhs)) {
    tp <- if (inherits(lhs, "symmatrix")) "ImmutableDenseMatrix" else .cg_pytype(lhs)
    .cg_err("TypeError", sprintf("Cannot assign to lhs of type <class '%s'>.", tp))
  }
  ls <- .cg_shape(lhs); rs <- .cg_shape(rhs)
  if (!is.null(ls)) {
    if (is.null(rs)) .cg_err("ValueError", "Cannot assign a scalar to a matrix.")
    if (!(identical(ls[[1L]]$key, rs[[1L]]$key) && identical(ls[[2L]]$key, rs[[2L]]$key)))
      .cg_err("ValueError", "Dimensions of lhs and rhs do not align.")
  } else if (!is.null(rs)) .cg_err("ValueError", "Cannot assign a matrix to a scalar.")
}

# The SymPy module path of a node's class, for messages that print type(x).
.cg_pytype <- function(e) {
  t <- e$type
  mods <- c(Integer = "sympy.core.numbers.Integer", Mul = "sympy.core.mul.Mul", Add = "sympy.core.add.Add",
            Pow = "sympy.core.power.Pow", MatAdd = "sympy.matrices.expressions.matadd.MatAdd",
            IndexedBase = "sympy.tensor.indexed.IndexedBase", Symbol = "sympy.core.symbol.Symbol",
            Rational = "sympy.core.numbers.Rational", Float = "sympy.core.numbers.Float",
            Idx = "sympy.tensor.indexed.Idx", One = "sympy.core.numbers.One", Zero = "sympy.core.numbers.Zero",
            MatMul = "sympy.matrices.expressions.matmul.MatMul", ImmutableMatrix = "sympy.matrices.immutable.ImmutableDenseMatrix",
            Indexed = "sympy.tensor.indexed.Indexed", MatrixSymbol = "sympy.matrices.expressions.matexpr.MatrixSymbol",
            Derivative = "sympy.core.function.Derivative")
  if (t == "Integer") {
    if (e$value == 0) t <- "Zero" else if (e$value == 1) t <- "One" else if (e$value == -1) t <- "NegativeOne"
  }
  if (t == "Symbol" && .is_idx(e)) t <- "Idx"
  m <- mods[t]
  if (is.na(m)) {
    if (t %in% c("Zero", "One", "NegativeOne")) return(paste0("sympy.core.numbers.", t))
    return(t)
  }
  unname(m)
}

.cg_assign_new <- function(cls, lhs, rhs) {
  lhs <- .cg_strict_sympify(lhs)
  rhs <- .cg_strict_sympify(rhs)
  .cg_check_assign(lhs, rhs)
  new_sym(cls, list(lhs, rhs))
}

.AUG_OPS <- c(AddAugmentedAssignment = "+", SubAugmentedAssignment = "-", MulAugmentedAssignment = "*",
              DivAugmentedAssignment = "/", ModAugmentedAssignment = "%")

.cg_def("AssignmentBase", token = FALSE)
.cg_def("Assignment", "AssignmentBase", op = ":=")
.cg_def("AugmentedAssignment", "AssignmentBase")
for (.cg_n in names(.AUG_OPS)) .cg_def(.cg_n, "AugmentedAssignment", binop = .AUG_OPS[[.cg_n]])
rm(.cg_n)

#' @rdname codegen-ast
#' @export
Assignment <- function(lhs, rhs) .cg_assign_new("Assignment", lhs, rhs)
#' @rdname codegen-ast
#' @export
AddAugmentedAssignment <- function(lhs, rhs) .cg_assign_new("AddAugmentedAssignment", lhs, rhs)
#' @rdname codegen-ast
#' @export
SubAugmentedAssignment <- function(lhs, rhs) .cg_assign_new("SubAugmentedAssignment", lhs, rhs)
#' @rdname codegen-ast
#' @export
MulAugmentedAssignment <- function(lhs, rhs) .cg_assign_new("MulAugmentedAssignment", lhs, rhs)
#' @rdname codegen-ast
#' @export
DivAugmentedAssignment <- function(lhs, rhs) .cg_assign_new("DivAugmentedAssignment", lhs, rhs)
#' @rdname codegen-ast
#' @export
ModAugmentedAssignment <- function(lhs, rhs) .cg_assign_new("ModAugmentedAssignment", lhs, rhs)

#' @rdname codegen-ast
#' @export
aug_assign <- function(lhs, op, rhs) {
  cls <- names(.AUG_OPS)[match(op, .AUG_OPS)]
  if (is.na(cls)) .cg_err("ValueError", paste0("Unrecognized operator ", op))
  .cg_assign_new(cls, lhs, rhs)
}

# ---------------------------------------------------------------------------
# CodeBlock

.cg_def("CodeBlock", token = FALSE, str = function(e, repr, ...) .cg_codeblock_print(e, repr),
        repr = function(e, repr, ...) .cg_codeblock_print(e, repr))

#' @rdname codegen-ast
#' @export
CodeBlock <- function(...) new_sym("CodeBlock", lapply(list(...), .cg_as_arg))

.cg_lhs <- function(e) Filter(Negate(is.null), lapply(e$args, function(a) if (a$type == "Assignment") a$args[[1L]]))
.cg_rhs <- function(e) Filter(Negate(is.null), lapply(e$args, function(a) if (a$type == "Assignment") a$args[[2L]]))

register_bound_vars("CodeBlock", function(e) vapply(.cg_lhs(e), .k, ""))

#' CodeBlock methods
#'
#' `codeblock_topological_sort(assignments)` orders a list of
#' [Assignment()] nodes so that every variable is assigned before it is
#' used, keeping the given order where it is free, as SymPy's
#' `CodeBlock.topological_sort`. `codeblock_cse(block)` runs common
#' subexpression elimination over the right hand sides of a block of
#' assignments.
#'
#' @param assignments A list of [Assignment()] nodes.
#' @param block A [CodeBlock()] of assignments.
#' @param symbols Optional list of symbols for the new variables.
#' @param optimizations,postprocess Passed on to [cse()].
#' @param order Ordering of the subexpressions.
#' @return A [CodeBlock()].
#' @examples
#' x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
#' codeblock_topological_sort(list(Assignment(x, y + z), Assignment(z, 1), Assignment(y, 2)))
#' @name codeblock-methods
NULL

#' @rdname codeblock-methods
#' @export
codeblock_topological_sort <- function(assignments) {
  assignments <- .cg_iter(assignments)
  if (!all(vapply(assignments, function(a) is_sym(a) && a$type == "Assignment", TRUE)))
    .cg_err("NotImplementedError", "CodeBlock.topological_sort only supports Assignments")
  n <- length(assignments)
  lhs_keys <- vapply(assignments, function(a) a$args[[1L]]$key, "")
  E <- list()
  for (j in seq_len(n)) {
    for (s in .cg_free_symbol_list(assignments[[j]]$args[[2L]])) {
      for (i in which(lhs_keys == s$key)) E[[length(E) + 1L]] <- c(i, j)
    }
  }
  .cg_topo(seq_len(n), E, assignments)
}

.cg_free_symbol_list <- function(e) free_symbols(e)

# sympy.utilities.iterables.topological_sort on integer vertices, key the
# vertex itself
.cg_topo <- function(V, E, items) {
  L <- integer()
  targets <- vapply(E, `[`, 0, 2L)
  S <- sort(setdiff(V, targets), decreasing = TRUE)
  while (length(S)) {
    node <- S[length(S)]
    S <- S[-length(S)]
    L <- c(L, node)
    k <- 1L
    for (edge in E) {
      if (edge[1L] == node) {
        idx <- which(vapply(E, function(z) identical(z, edge), TRUE))[1L]
        E <- E[-idx]
        v <- edge[2L]
        if (!any(vapply(E, function(z) z[2L] == v, TRUE))) {
          pos <- which(v > S)
          if (length(pos)) S <- append(S, v, after = pos[1L] - 1L) else S <- c(S, v)
        }
      }
    }
  }
  if (length(E)) .cg_err("ValueError", "cycle detected")
  do.call(CodeBlock, items[L])
}

#' @rdname codeblock-methods
#' @export
codeblock_cse <- function(block, symbols = NULL, optimizations = NULL, postprocess = NULL, order = "canonical") {
  if (!all(vapply(block$args, function(a) a$type == "Assignment", TRUE)))
    .cg_err("NotImplementedError", "CodeBlock.cse only supports Assignments")
  lhs <- .cg_lhs(block)
  for (i in seq_along(lhs)) {
    if (i > 1L && lhs[[i]]$key %in% vapply(lhs[seq_len(i - 1L)], .k, ""))
      .cg_err("NotImplementedError", sprintf("Duplicate assignments to the same variable are not yet supported (%s)", .str(lhs[[i]])))
  }
  existing <- vapply(atoms(block, "Symbol"), function(s) s$name, "")
  gen <- function(k) {
    out <- list(); i <- 0L
    while (length(out) < k) {
      nm <- paste0("x", i); i <- i + 1L
      if (!(nm %in% existing)) out[[length(out) + 1L]] <- Symbol(nm)
    }
    out
  }
  rhs <- .cg_rhs(block)
  res <- .cg_run_cse(rhs, if (is.null(symbols)) gen else function(k) symbols[seq_len(k)], optimizations, postprocess, order)
  new_block <- Map(function(v, ex) Assignment(v, ex), lhs, res$reduced)
  new_assign <- lapply(res$replacements, function(p) Assignment(p[[1L]], p[[2L]]))
  codeblock_topological_sort(c(new_assign, new_block))
}

.cg_run_cse <- function(exprs, symgen, optimizations, postprocess, order) {
  exprs <- lapply(exprs, function(z) if (is_sym(z) && z$type == "ImmutableMatrix") .unwrap_explicit(z) else z)
  # a matrix enters as its entries (SymPy's cse takes Tuple(*m.flat()))
  lens <- vapply(exprs, function(z) if (inherits(z, "symmatrix")) length(z$data) else 1L, 1L)
  flat <- unlist(lapply(exprs, function(z) if (inherits(z, "symmatrix")) z$data else list(z)), recursive = FALSE)
  # ponytail: the package cse() (no optimizations, order or ignore), then
  # its x<k> symbols renamed to the requested ones
  r <- cse(flat, prefix = "\001cse")
  ends <- cumsum(lens)
  red <- r$reduced
  r$reduced <- Map(function(z, e, n) {
    if (!inherits(z, "symmatrix")) return(red[[e]])
    z$data <- red[seq.int(e - n + 1L, e)]
    z
  }, exprs, ends, lens)
  syms <- symgen(length(r$replacements))
  ren <- Map(function(pr, s) list(pr[[1L]], s), r$replacements, syms)
  r$replacements <- Map(function(pr, s) list(s, xreplace(pr[[2L]], ren)), r$replacements, syms)
  r$reduced <- lapply(r$reduced, function(z) {
    if (!length(ren)) return(z)
    if (inherits(z, "symmatrix")) { z$data <- lapply(z$data, xreplace, rule = ren); return(z) }
    xreplace(z, ren)
  })
  r$reduced <- lapply(r$reduced, function(z) if (inherits(z, "symmatrix")) .as_block(z) else z)
  r
}

# ---------------------------------------------------------------------------
# Loops, scopes and functions

.cg_def("For", "Token", fields = c("target", "iterable", "body"),
        construct = list(target = function(a) .cg_strict_sympify(a), body = .cg_body,
                         iterable = function(itr) {
                           if (!.cg_is_iterable(itr)) .cg_err("TypeError", "iterable must be an iterable")
                           .cg_sympify(itr)
                         }))
#' @rdname codegen-ast
#' @export
For <- function(...) .cg_new("For", list(...))

.cg_def("While", "Token", fields = c("condition", "body"),
        construct = list(condition = function(a) .cg_strict_sympify(a), body = .cg_body))
#' @rdname codegen-ast
#' @export
While <- function(...) .cg_new("While", list(...))

.cg_def("Scope", "Token", fields = "body", construct = list(body = .cg_body))
#' @rdname codegen-ast
#' @export
Scope <- function(...) .cg_new("Scope", list(...))

.cg_def("Stream", "Token", fields = "name", construct = list(name = function(a) String(a)))
#' @rdname codegen-ast
#' @export
Stream <- function(...) .cg_new("Stream", list(...))
#' @rdname codegen-ast
#' @export
stdout_ <- Stream("stdout")
#' @rdname codegen-ast
#' @export
stderr_ <- Stream("stderr")

.cg_def("Print", "Token", fields = c("print_args", "format_string", "file"),
        defaults = list(format_string = none, file = none),
        construct = list(print_args = function(a) .cg_mk_tuple(a), format_string = function(a) QuotedString(a),
                         file = function(a) Stream(a)))
#' @rdname codegen-ast
#' @export
Print <- function(...) .cg_new("Print", list(...))

.cg_fp_params <- function(args) {
  .cg_tuple(lapply(.cg_iter(args), function(a) {
    if (.cg_isa(a, "Declaration")) a$args[[1L]]
    else if (.cg_isa(a, "Variable")) a
    else variable_deduced(a)
  }))
}

.cg_def("FunctionPrototype", "Node", fields = c("return_type", "name", "parameters", "attrs"),
        construct = list(return_type = function(a) Type(a), name = function(a) String(a), parameters = .cg_fp_params))
.cg_def("FunctionDefinition", "FunctionPrototype", fields = c("return_type", "name", "parameters", "body", "attrs"),
        construct = list(body = .cg_body))
#' @rdname codegen-ast
#' @export
FunctionPrototype <- function(...) .cg_new("FunctionPrototype", list(...))
#' @rdname codegen-ast
#' @export
FunctionDefinition <- function(...) .cg_new("FunctionDefinition", list(...))

#' @rdname codegen-types
#' @param func_def A [FunctionDefinition()].
#' @export
from_FunctionDefinition <- function(func_def) {
  if (!.cg_isa(func_def, "FunctionDefinition"))
    .cg_err("TypeError", "func_def is not an instance of FunctionDefinition")
  do.call(FunctionPrototype, ast_kwargs(func_def, exclude = "body"))
}

#' @rdname codegen-types
#' @param func_proto A [FunctionPrototype()].
#' @param body The body statements.
#' @export
from_FunctionPrototype <- function(func_proto, body) {
  if (!.cg_isa(func_proto, "FunctionPrototype"))
    .cg_err("TypeError", "func_proto is not an instance of FunctionPrototype")
  do.call(FunctionDefinition, c(list(body = body), ast_kwargs(func_proto)))
}

.cg_def("Return", "Token", fields = "return", construct = list(`return` = function(a) .cg_strict_sympify(a)))
#' @rdname codegen-ast
#' @export
Return <- function(...) .cg_new("Return", list(...))

.cg_def("FunctionCall", "Token", fields = c("name", "function_args"),
        construct = list(name = function(a) String(a),
                         function_args = function(a) .cg_tuple(lapply(.cg_iter(a), .cg_sympify))))
#' @rdname codegen-ast
#' @export
FunctionCall <- function(...) .cg_new("FunctionCall", list(...))

.cg_def("KeywordFunctionCall", "FunctionCall", fields = c("name", "function_args", "keyword_args"),
        defaults = list(keyword_args = new_sym("Dict")),
        construct = list(keyword_args = function(kw) {
          if (is_sym(kw)) return(kw)
          if (!length(kw)) return(new_sym("Dict"))
          Dict(lapply(names(kw), String), unname(kw))
        }))
#' @rdname codegen-ast
#' @export
KeywordFunctionCall <- function(...) .cg_new("KeywordFunctionCall", list(...))

.cg_def("Raise", "Token", fields = "exception")
#' @rdname codegen-ast
#' @export
Raise <- function(...) .cg_new("Raise", list(...))
.cg_def("RuntimeError_", "Token", fields = "message", construct = list(message = function(a) String(a)))
#' @rdname codegen-ast
#' @export
RuntimeError_ <- function(...) .cg_new("RuntimeError_", list(...))

# ---------------------------------------------------------------------------
# Generic access to the fields and methods of nodes

#' Fields of code generation nodes
#'
#' `ast_attr(node, name)` returns a field or derived attribute of a codegen
#' node, as the Python attribute `node.name`: for example
#' `ast_attr(v, "symbol")`, `ast_attr(block, "left_hand_sides")` or
#' `ast_attr(float64, "eps")`. `ast_kwargs(node)` returns the fields as a
#' named list (`Token.kwargs()`), and `ast_func(node)` the constructor, so
#' `do.call(ast_func(node), sym_args(node))` rebuilds the node.
#' `ast_getitem(p, i)` indexes a [Pointer()], giving an [Element()].
#'
#' @param node A codegen node.
#' @param name The attribute name.
#' @param exclude Field names to leave out.
#' @param apply Optional function applied to every value.
#' @param p A [Pointer()].
#' @param key An index or a list of indices.
#' @return The attribute value, a list, or a function.
#' @examples
#' x <- Symbol("x")
#' v <- Variable(x, type = real)
#' ast_attr(v, "type")
#' ast_kwargs(v)
#' @name ast-attributes
NULL

#' @rdname ast-attributes
#' @export
ast_attr <- function(node, name) {
  if (inherits(node, "cg_printer")) return(.cgp_attr_get(node, name))
  if (inherits(node, c("cg_cgobj", "cg_codegen"))) return(.cgu_attr(node, name))
  cl <- .cg_cls(node)
  if (is.null(cl) || (startsWith(name, "is_") && name != "is_Atom")) {
    r <- .py_attr_generic(node, name)
    if (!is.null(r)) return(r)
    .cg_err("AttributeError", sprintf("'%s' object has no attribute '%s'", if (is_sym(node)) node$type else class(node)[1L], name))
  }
  if (name %in% cl$fields) return(.cg_vals(node)[[name]])
  if (!is.null(cl$attrs[[name]])) return(cl$attrs[[name]](node))
  r <- switch(name,
              is_Atom = if (node$type %in% c("String", "QuotedString", "Comment")) TRUE else length(cl$fields) == 0L,
              lhs = if ("AssignmentBase" %in% cl$isa) node$args[[1L]],
              rhs = if ("AssignmentBase" %in% cl$isa) node$args[[2L]],
              op = if (node$type == "Assignment") ":=" else if ("AugmentedAssignment" %in% cl$isa) paste0(cl$binop, "="),
              binop = if ("AugmentedAssignment" %in% cl$isa) cl$binop,
              left_hand_sides = if (node$type == "CodeBlock") .cg_tuple(.cg_lhs(node)),
              right_hand_sides = if (node$type == "CodeBlock") .cg_tuple(.cg_rhs(node)),
              min = , max = if ("_SizedIntType" %in% cl$isa) {
                nb <- .int_val(node$args[[2L]])
                signed <- "SignedIntType" %in% cl$isa
                v <- if (name == "min") (if (signed) -gmp::pow.bigz(2, nb - 1) else gmp::as.bigz(0))
                     else (if (signed) gmp::pow.bigz(2, nb - 1) - 1 else gmp::pow.bigz(2, nb) - 1)
                make_rational(gmp::as.bigq(v))
              } else if ("FloatType" %in% cl$isa && name == "max") .ft_attr(node, "max"),
              max_exponent = , min_exponent = , tiny = , eps = , dig = , decimal_dig =
                if ("FloatType" %in% cl$isa) .ft_attr(node, name),
              NULL)
  if (is.null(r) && !(name %in% c("lhs", "rhs")))
    .cg_err("AttributeError", sprintf("'%s' object has no attribute '%s'", node$type, name))
  r
}

# Attributes of other objects: assumptions, Pow parts and a few more
.py_attr_generic <- function(node, name) {
  if (inherits(node, "symmatrix")) {
    return(switch(name, shape = list(node$nrow, node$ncol), rows = node$nrow, cols = node$ncol, NULL))
  }
  if (!is_sym(node)) return(NULL)
  if (startsWith(name, "is_")) {
    nm <- substring(name, 4L)
    cls <- c(Pow = "Pow", Add = "Add", Mul = "Mul", Symbol = "Symbol", Function = "Function",
             Float = "Float", Integer = "Integer", Rational = "Rational", Relational = "Relational")
    if (nm %in% names(cls)) {
      return(switch(nm, Symbol = node$type %in% c("Symbol", "Dummy"),
                    Function = .cgp_is_function(node),
                    Integer = node$type == "Integer",
                    Rational = node$type %in% c("Integer", "Rational"),
                    Relational = node$type %in% names(.REL_OPS),
                    node$type == cls[[nm]]))
    }
    if (nm == "Number") return(is_number_type(node))
    if (nm == "number") return(.expr_is_number(node))
    if (nm == "Atom") return(!length(node$args))
    v <- tryCatch(ask_fact(node, nm), error = function(err) NA)
    return(v)
  }
  switch(name,
         args = node$args,
         func = ast_func(node),
         base = if (node$type == "Pow") node$args[[1L]] else .base_exp(node)[[1L]],
         exp = if (node$type == "Pow") node$args[[2L]] else .base_exp(node)[[2L]],
         name = node$name,
         lhs = if (node$type %in% names(.REL_OPS)) node$args[[1L]],
         rhs = if (node$type %in% names(.REL_OPS)) node$args[[2L]],
         label = if (node$type == "IndexedBase") node$args[[1L]] else if (.is_idx(node)) Symbol(node$name),
         free_symbols = free_symbols(node),
         shape = if (node$type == "MatrixSymbol" || .is_matexpr(node)) .mshape(node),
         NULL)
}

#' @rdname ast-attributes
#' @export
ast_kwargs <- function(node, exclude = character(), apply = NULL) {
  cl <- .cg_cls(node)
  if (node$type %in% c("String", "QuotedString", "Comment")) return(list())
  vals <- .cg_vals(node)
  vals <- vals[setdiff(cl$fields, exclude)]
  if (!is.null(apply)) vals <- lapply(vals, apply)
  vals
}

#' @rdname ast-attributes
#' @export
ast_func <- function(node) {
  if (!is_sym(node)) stop("not a node")
  if (node$type %in% c("String", "QuotedString", "Comment")) return(function(...) node)
  f <- .cg_ctor(node$type)
  if (!is.null(f)) return(f)
  if (node$type == "Function" && is.null(node$extra$imp))
    for (nm in c(node$name, paste0("sym_", node$name))) {
      f <- get0(nm, envir = environment(ast_func), inherits = FALSE, mode = "function")
      if (!is.null(f)) return(f)
    }
  function(...) sym_rebuild(node, lapply(list(...), as_sym))
}

.cg_ctor <- function(type) {
  ns <- environment(.cg_ctor)
  cl <- .cg_classes[[type]]
  nm <- if (!is.null(cl$ctor)) cl$ctor else type
  if (!is.null(cl) && exists(nm, envir = ns, inherits = FALSE)) {
    f <- get(nm, envir = ns)
    if (is.function(f)) return(f)
  }
  NULL
}

#' @rdname ast-attributes
#' @export
ast_getitem <- function(p, key) {
  sym <- .cg_vals(p)$symbol
  if (.cg_is_iterable(key) && !is.character(key)) Element(sym, key) else Element(sym, list(key))
}

#' @rdname ast-attributes
#' @param looking_for The attribute name.
#' @export
attr_params <- function(node, looking_for) {
  for (a in .cg_vals(node)$attrs$args) {
    if (identical(.str(a$args[[1L]]), as.character(if (is_sym(looking_for)) .str(looking_for) else looking_for)))
      return(a$args[[2L]])
  }
  NULL
}
