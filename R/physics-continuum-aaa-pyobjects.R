# Objects of the physics modules. The file name makes it collate before the
# control, continuum mechanics, optics and misc modules that build on it.
#
# SymPy's physics modules define classes with properties and methods. An
# immutable one (a Basic subclass such as TransferFunction or BeamParameter)
# is a sym node of class c(<class>, <parent classes>, "pyobj", "sym"); its
# arguments are sym nodes, so subs(), free_symbols() and printing work as for
# any expression. A mutable one (Beam, Truss, Cable) is an environment of
# class c(<class>, <parent classes>, "pyenv"). For both, `obj$name` reads a
# property or returns the bound method, as `obj.name` does in Python.
# A class with class methods (TransferFunction.from_rational_expression) is
# exported as a constructor function of class "pyclass", so
# `TransferFunction$from_rational_expression(...)` works.

.PYCLASSES <- new.env(hash = TRUE)

# Define a class. `props` and `methods` are named lists of functions whose
# first argument is the object; `ops` maps Python operator names ("__add__",
# "__radd__", "__neg__", ...) to functions of (self, other); `str` gives the
# printed form; `ctor` rebuilds the object from new arguments; `expr` marks a
# subclass of Expr, which falls back on expression arithmetic.
.pyclass_def <- function(name, parent = NULL, props = list(), methods = list(), ops = list(),
                         str = NULL, ctor = NULL, statics = list(), setters = list(), expr = FALSE) {
  assign(name, list(name = name, parent = parent, props = props, methods = methods, ops = ops,
                    str = str, ctor = ctor, statics = statics, setters = setters, expr = expr),
         envir = .PYCLASSES)
  if (!is.null(str)) register_str_printer(name, str)
  if (!is.null(ctor)) register_rebuilder(name, function(e, args) do.call(ctor, args))
  invisible(name)
}

.pyclass_chain <- function(name) {
  out <- character()
  while (!is.null(name)) {
    out <- c(out, name)
    d <- .PYCLASSES[[name]]
    name <- if (is.null(d)) NULL else d$parent
  }
  out
}

.pyclass_lookup <- function(cls, field, name) {
  for (c in .pyclass_chain(cls)) {
    d <- .PYCLASSES[[c]]
    if (!is.null(d) && !is.null(d[[field]][[name]])) return(d[[field]][[name]])
  }
  NULL
}

# the structural key of an argument that may not be a sym node
.pyobj_k <- function(a) {
  if (is_sym(a)) return(a$key)
  if (inherits(a, "symmatrix")) return(paste0("M", a$nrow, "x", a$ncol, "[", paste(vapply(a$data, .pyobj_k, ""), collapse = ","), "]"))
  if (is.null(a)) return("None")
  if (is.environment(a)) return(paste0("env:", format(a)))
  if (is.list(a)) return(paste0("(", paste(vapply(a, .pyobj_k, ""), collapse = ","), ")"))
  paste0(class(a)[1L], ":", paste(format(a, digits = 17), collapse = ","))
}

# A new immutable object.
.pyobj_new <- function(type, args, name = NULL, extra = NULL) {
  key <- paste0(type, if (!is.null(name)) paste0("[", name, "]"), "(", paste(vapply(args, .pyobj_k, ""), collapse = ","), ")")
  o <- new_sym(type, args, name = name, extra = extra, key = key)
  class(o) <- c(.pyclass_chain(type), "pyobj", "sym")
  o
}

# A new mutable object.
.pyenv_new <- function(type) {
  e <- new.env(parent = emptyenv())
  class(e) <- c(.pyclass_chain(type), "pyenv")
  e
}

.SYM_FIELDS <- c("type", "args", "name", "value", "extra", "key")

#' @export
`$.pyobj` <- function(x, name) {
  if (name %in% .SYM_FIELDS) return(.subset2(x, name))
  .pyobj_get(x, name)
}

#' @export
`$.pyenv` <- function(x, name) .pyobj_get(x, name)

#' @export
`$<-.pyenv` <- function(x, name, value) {
  .pyobj_set(x, name, value)
  x
}

#' @export
`$.pyclass` <- function(x, name) .pyobj_get(x, name)

# Python's type(x).__name__
.py_typename <- function(x) {
  if (inherits(x, c("pyobj", "pyenv"))) return(class(x)[1L])
  if (inherits(x, "pyclass")) return("type")
  if (inherits(x, "symmatrix")) return(if (!is.null(attr(x, "pyclass"))) attr(x, "pyclass") else "MutableDenseMatrix")
  if (is_sym(x)) return(.sympy_class(x))
  if (is.null(x)) return("NoneType")
  if (is.logical(x)) return("bool")
  if (is.integer(x)) return("int")
  if (is.double(x)) return(if (length(x) == 1L && x == round(x) && abs(x) < 2^53) "int" else "float")
  if (is.complex(x)) return("complex")
  if (is.character(x)) return("str")
  if (is.function(x)) return("function")
  if (is.list(x)) return("list")
  class(x)[1L]
}

.pyobj_get <- function(x, name) {
  if (inherits(x, "pyclass")) {
    f <- .pyclass_lookup(attr(x, "pyname"), "statics", name)
    if (!is.null(f)) return(f)
    .physics_error("AttributeError", sprintf("type object '%s' has no attribute '%s'", attr(x, "pyname"), name))
  }
  cls <- class(x)[1L]
  f <- .pyclass_lookup(cls, "props", name)
  if (!is.null(f)) return(f(x))
  m <- .pyclass_lookup(cls, "methods", name)
  if (!is.null(m)) return(function(...) m(x, ...))
  if (is.environment(x)) {
    if (exists(name, envir = x, inherits = FALSE)) return(get(name, envir = x, inherits = FALSE))
  } else {
    m <- .pyobj_basic_method(x, name)
    if (!is.null(m)) return(m)
  }
  .physics_error("AttributeError", sprintf("'%s' object has no attribute '%s'", cls, name))
}

.pyobj_set <- function(x, name, value) {
  if (!is.environment(x)) .physics_error("AttributeError", sprintf("'%s' object attribute '%s' is read-only", class(x)[1L], name))
  f <- .pyclass_lookup(class(x)[1L], "setters", name)
  if (!is.null(f)) f(x, value) else assign(name, value, envir = x)
  invisible(x)
}

# The methods every Basic object has.
.pyobj_basic_method <- function(x, name) {
  switch(name,
    func = .pyobj_ctor(class(x)[1L]),
    args = .pytuple(.subset2(x, "args")),
    free_symbols = free_symbols(x),
    subs = function(...) {
      for (p in .pyobj_pairs(...)) x <- subs(x, p[[1L]], p[[2L]])
      x
    },
    xreplace = function(...) .pyobj_xreplace(x, .pyobj_pairs(...)),
    doit = function(...) doit(x),
    has = function(...) has(x, ...),
    atoms = function(...) atoms(x, ...),
    evalf = , n = function(...) .pyobj_evalf(x, ...),
    rewrite = function(rule, ...) .pyobj_rewrite(x, rule),
    NULL)
}

# the old and new values of subs(old, new), subs(dict) or subs(pairs)
.pyobj_pairs <- function(...) {
  a <- list(...)
  if (length(a) == 1L && inherits(a[[1L]], "pydict")) return(Map(list, a[[1L]]$keys, a[[1L]]$values))
  if (length(a) == 1L && is.list(a[[1L]]) && !is_sym(a[[1L]])) return(lapply(a[[1L]], function(p) list(p[[1L]], p[[2L]])))
  list(list(a[[1L]], a[[2L]]))
}

.pyobj_ctor <- function(cls) .pyclass_lookup(cls, "statics", "__call__") %||% .PYCLASSES[[cls]]$ctor

# Basic.xreplace: rebuild through the constructor when an argument changed
.pyobj_xreplace <- function(x, pairs) {
  keys <- vapply(pairs, function(p) .pyobj_k(as_sym(p[[1L]])), "")
  walk <- function(e) {
    if (!is_sym(e)) return(e)
    i <- match(.pyobj_k(e), keys)
    if (!is.na(i)) return(as_sym(pairs[[i]][[2L]]))
    if (!inherits(e, "pyobj")) return(xreplace(e, pairs))
    old <- .subset2(e, "args")
    args <- lapply(old, walk)
    if (identical(vapply(args, .pyobj_k, ""), vapply(old, .pyobj_k, ""))) return(e)
    do.call(.pyobj_ctor(class(e)[1L]), args)
  }
  walk(x)
}

# factor() of an object that is not an expression factors its arguments
.pyobj_factor <- function(x, ...) {
  args <- lapply(.subset2(x, "args"), function(a) {
    if (inherits(a, "pyobj")) .pyobj_factor(a, ...) else if (is_sym(a)) factor_sym(a, ...) else a
  })
  do.call(.pyobj_ctor(class(x)[1L]), args)
}

# Python dictionaries and lists with keys and items compared as Python does
.pd_empty <- function() .pydict(list(), list())
.pd_idx <- function(d, k) match(.pyobj_k(k), vapply(d$keys, .pyobj_k, ""))
.pd_has <- function(d, k) !is.na(.pd_idx(d, k))
.pd_get <- function(d, k) {
  i <- .pd_idx(d, k)
  if (is.na(i)) .physics_error("KeyError", .pyfmt(k))
  d$values[[i]]
}
.pd_set <- function(d, k, v) {
  i <- .pd_idx(d, k)
  if (is.na(i)) {
    d$keys <- c(d$keys, list(k))
    d$values <- c(d$values, list(v))
  } else d$values[i] <- list(v)
  d
}
.pd_pop <- function(d, k) {
  i <- .pd_idx(d, k)
  if (is.na(i)) .physics_error("KeyError", .pyfmt(k))
  d$keys <- d$keys[-i]
  d$values <- d$values[-i]
  d
}
.pl_index <- function(l, v) match(.pyobj_k(v), vapply(l, .pyobj_k, ""))
.pl_remove <- function(l, v) {
  i <- .pl_index(l, v)
  if (is.na(i)) .physics_error("ValueError", "list.remove(x): x not in list")
  l[-i]
}

# Basic.rewrite: the _eval_rewrite_as_<rule> method of the class, if any
# (deep: the arguments are rewritten first, as Basic._rewrite does)
.pyobj_rewrite <- function(x, rule) {
  if (is.function(rule)) rule <- attr(rule, "pyname")
  if (!is.character(rule)) return(x)
  deep <- function(a) {
    if (inherits(a, "pyobj")) return(.pyobj_rewrite(a, rule))
    if (is_sym(a) && a$type == "Tuple") {
      new <- lapply(a$args, deep)
      if (!identical(vapply(new, .pyobj_k, ""), vapply(a$args, .pyobj_k, ""))) return(do.call(Tuple, new))
    }
    a
  }
  old <- .subset2(x, "args")
  args <- lapply(old, deep)
  m <- .pyclass_lookup(class(x)[1L], "methods", paste0("_eval_rewrite_as_", rule))
  if (!is.null(m)) return(m(x))
  if (identical(vapply(args, .pyobj_k, ""), vapply(old, .pyobj_k, ""))) return(x)
  do.call(.pyobj_ctor(class(x)[1L]), args)
}

# simplify(), expand() and N() of an object: its own method, if any
.pyobj_simplify <- function(x) {
  m <- .pyclass_lookup(class(x)[1L], "methods", "simplify")
  if (is.null(m)) x$doit() else m(x)
}

.pyobj_expand <- function(x) {
  m <- .pyclass_lookup(class(x)[1L], "methods", "expand")
  if (is.null(m) && isTRUE(.PYCLASSES[[class(x)[1L]]]$expr))
    return(do.call(.PYCLASSES[[class(x)[1L]]]$ctor, lapply(x$args, function(a) if (is_sym(a)) expand(a) else a)))
  if (is.null(m)) .physics_error("AttributeError", sprintf("'%s' object has no attribute 'expand'", class(x)[1L]))
  m(x)
}

# Basic.evalf for an object that is not an expression: evaluate the arguments
.pyobj_evalf <- function(x, ...) {
  m <- .pyclass_lookup(class(x)[1L], "methods", "evalf")
  if (!is.null(m)) return(m(x, ...))
  f <- .PYCLASSES[[class(x)[1L]]]$ctor
  args <- lapply(x$args, function(a) if (is_sym(a)) N(a, ...) else a)
  do.call(f, args)
}

.pyobj_isinstance <- function(x, cls) inherits(x, cls)

.pyobj_str <- function(x) {
  if (is.environment(x)) {
    f <- .pyclass_lookup(class(x)[1L], "methods", "__str__")
    return(if (is.null(f)) paste0("<", class(x)[1L], " object>") else f(x))
  }
  sstr(x)
}

#' @export
print.pyenv <- function(x, ...) {
  cat(.pyobj_str(x), "\n", sep = "")
  invisible(x)
}

#' @export
format.pyenv <- function(x, ...) .pyobj_str(x)

#' @export
print.pyclass <- function(x, ...) {
  cat("<class '", attr(x, "pyname"), "'>\n", sep = "")
  invisible(x)
}

# A constructor that also carries the class methods.
.pyclass_fn <- function(fn, name) structure(fn, class = c("pyclass", "function"), pyname = name)

#' @export
chooseOpsMethod.pyobj <- function(x, y, mx, my, cl, reverse) TRUE

.OPS_DUNDER <- c("+" = "add", "-" = "sub", "*" = "mul", "/" = "truediv", "^" = "pow", "%%" = "mod")

.py_binop_error <- function(op, a, b) {
  pyop <- c("+" = "+", "-" = "-", "*" = "*", "/" = "/", "^" = "** or pow()", "%%" = "%")[[op]]
  .physics_error("TypeError", sprintf("unsupported operand type(s) for %s: '%s' and '%s'", pyop, .py_typename(a), .py_typename(b)))
}

# Python's binary operator protocol: a.__op__(b), then b.__rop__(a).
.py_binop <- function(op, a, b) {
  d <- .OPS_DUNDER[[op]]
  if (inherits(a, c("pyobj", "pyenv"))) {
    f <- .pyclass_lookup(class(a)[1L], "ops", paste0("__", d, "__"))
    if (!is.null(f)) { r <- f(a, b); if (!is.null(r)) return(r) }
  }
  if (inherits(b, c("pyobj", "pyenv"))) {
    f <- .pyclass_lookup(class(b)[1L], "ops", paste0("__r", d, "__"))
    if (!is.null(f)) { r <- f(b, a); if (!is.null(r)) return(r) }
  }
  NULL
}

#' @export
Ops.pyobj <- function(e1, e2) {
  if (missing(e2)) {
    if (.Generic == "+") return(e1)
    f <- .pyclass_lookup(class(e1)[1L], "ops", "__neg__")
    if (!is.null(f)) return(f(e1))
    if (isTRUE(.PYCLASSES[[class(e1)[1L]]]$expr)) return(.neg(e1))
    .physics_error("TypeError", sprintf("bad operand type for unary -: '%s'", class(e1)[1L]))
  }
  if (.Generic %in% c("==", "!=")) {
    eq <- is_sym(e1) && is_sym(e2) && identical(.subset2(e1, "key"), .subset2(as_sym(e2), "key"))
    return(if (.Generic == "==") eq else !eq)
  }
  if (.Generic %in% names(.OPS_DUNDER)) {
    r <- .py_binop(.Generic, e1, e2)
    if (!is.null(r)) return(r)
    expr1 <- !inherits(e1, "pyobj") || isTRUE(.PYCLASSES[[class(e1)[1L]]]$expr)
    expr2 <- !inherits(e2, "pyobj") || isTRUE(.PYCLASSES[[class(e2)[1L]]]$expr)
    if (expr1 && expr2 && !inherits(e1, "symmatrix") && !inherits(e2, "symmatrix")) {
      a <- as_sym(e1); b <- as_sym(e2)
      return(switch(.Generic, "+" = .add(list(a, b)), "-" = .sub(a, b), "*" = .mul(list(a, b)),
                    "/" = .div(a, b), "^" = .pow(a, b), "%%" = fcall("Mod", list(a, b))))
    }
    .py_binop_error(.Generic, e1, e2)
  }
  .physics_error("TypeError", sprintf("'%s' not supported between instances of '%s' and '%s'",
                                      .Generic, .py_typename(e1), .py_typename(e2)))
}
