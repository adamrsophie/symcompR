# Vectors and dyadics, after sympy/vector/basisdependent.py, vector.py and
# dyadic.py. Vectors and dyadics are sym nodes of the types BaseVector,
# VectorAdd, VectorMul and VectorZero (BaseDyadic, DyadicAdd, DyadicMul and
# DyadicZero), with the extra S3 class "sym_vec" so that R operators reach
# the basis dependent arithmetic. A VectorMul keeps its measure number and
# base vector in `extra`; a VectorAdd keeps its components in insertion order,
# as SymPy's `_components` dictionary does.

.VEC_TYPES <- c("BaseVector", "VectorAdd", "VectorMul", "VectorZero")
.DYA_TYPES <- c("BaseDyadic", "DyadicAdd", "DyadicMul", "DyadicZero")
.BD_TYPES <- c(.VEC_TYPES, .DYA_TYPES)
.VFAM <- list(V = c(base = "BaseVector", add = "VectorAdd", mul = "VectorMul", zero = "VectorZero"),
              D = c(base = "BaseDyadic", add = "DyadicAdd", mul = "DyadicMul", zero = "DyadicZero"))

.vec_stop <- function(cls, msg) stop(.sympy_condition(cls, msg))

.is_vec <- function(e) is_sym(e) && e$type %in% .VEC_TYPES
.is_dya <- function(e) is_sym(e) && e$type %in% .DYA_TYPES
.is_bd <- function(e) is_sym(e) && e$type %in% .BD_TYPES
.bd_fam <- function(e) if (e$type %in% .VEC_TYPES) "V" else "D"
.is_bd_zero <- function(e) is_sym(e) && e$type %in% c("VectorZero", "DyadicZero")

.as_symvec <- function(o) {
  class(o) <- c("sym_vec", "sym")
  o
}

.VZERO <- .as_symvec(new_sym("VectorZero", key = "VectorZero()"))
.DZERO <- .as_symvec(new_sym("DyadicZero", key = "DyadicZero()"))
.bd_zero <- function(fam) if (fam == "V") .VZERO else .DZERO

# Python's `expr == 0` for a scalar: 0 and 0.0
.is_scalar_zero <- function(e) is_sym(e) && (is_zero_num(e) || (e$type == "Float" && num_is_zero(e)))

# _base_instance and _measure_number
.bd_base <- function(e) if (e$type %in% c("VectorMul", "DyadicMul")) e$extra$base else e
.bd_measure <- function(e) if (e$type %in% c("VectorMul", "DyadicMul")) e$extra$measure else S_One

# the components dictionary as a list of list(base, measure), in order
.bd_components <- function(e) {
  switch(e$type,
         BaseVector = , BaseDyadic = list(list(e, S_One)),
         VectorMul = , DyadicMul = list(list(e$extra$base, e$extra$measure)),
         VectorAdd = , DyadicAdd = e$extra$components,
         list())
}

# ---------------------------------------------------------------------------
# Constructors (BasisDependentMul._new and BasisDependentAdd.__new__)

.bd_mul_node <- function(fam, measure, base, extras = list()) {
  if (is_one_num(measure) && !length(extras)) return(base)
  facs <- if (measure$type == "Mul") measure$args else if (is_one_num(measure)) list() else list(measure)
  .as_symvec(new_sym(.VFAM[[fam]][["mul"]], c(facs, list(base), extras),
                     extra = list(measure = measure, base = base, extras = extras)))
}

.bd_mul <- function(args, fam) {
  F <- .VFAM[[fam]]
  count <- 0L
  measure <- S_One
  zeroflag <- FALSE
  extra_args <- list()
  expr <- NULL
  for (arg in args) {
    arg <- as_sym(arg)
    if (arg$type == F[["zero"]]) {
      count <- count + 1L
      zeroflag <- TRUE
    } else if (.is_scalar_zero(arg)) {
      zeroflag <- TRUE
    } else if (arg$type %in% c(F[["base"]], F[["mul"]])) {
      count <- count + 1L
      expr <- .bd_base(arg)
      measure <- .mul(list(measure, .bd_measure(arg)))
    } else if (arg$type == F[["add"]]) {
      count <- count + 1L
      expr <- arg
    } else if (arg$type %in% c("Cross", "Dot", "Curl", "Gradient")) {
      extra_args[[length(extra_args) + 1L]] <- arg
    } else {
      measure <- .mul(list(measure, arg))
    }
  }
  if (count > 1L) .vec_stop("ValueError", "Invalid multiplication")
  if (count == 0L) return(.mul(args))
  if (zeroflag) return(.bd_zero(fam))
  if (expr$type == F[["add"]])
    return(.bd_add(lapply(expr$args, function(x) .bd_mul(list(measure, x), fam)), fam))
  .bd_mul_node(fam, measure, .bd_base(expr), extra_args)
}

.bd_add <- function(args, fam) {
  F <- .VFAM[[fam]]
  comps <- list()
  for (arg in args) {
    arg <- as_sym(arg)
    if (!(arg$type %in% .VEC_TYPES && fam == "V") && !(arg$type %in% .DYA_TYPES && fam == "D")) {
      if (arg$type == "Mul") arg <- .bd_mul(arg$args, fam)
      else if (arg$type == "Add") arg <- .bd_add(arg$args, fam)
      else .vec_stop("TypeError", paste0(.str_top(arg), " cannot be interpreted correctly"))
      # a product that is not basis dependent has no components
      if (!.is_bd(arg)) .vec_stop("AttributeError", paste0("'", if (arg$type == "Function") arg$name else arg$type,
                                                         "' object has no attribute 'components'"))
    }
    if (arg$type == F[["zero"]]) next
    for (cm in .bd_components(arg)) {
      k <- cm[[1L]]$key
      comps[[k]] <- if (is.null(comps[[k]])) cm else list(cm[[1L]], .add(list(comps[[k]][[2L]], cm[[2L]])))
    }
  }
  comps <- Filter(function(cm) !.is_scalar_zero(cm[[2L]]), comps)
  if (!length(comps)) return(.bd_zero(fam))
  comps <- unname(comps)
  terms <- lapply(comps, function(cm) .bd_mul(list(cm[[1L]], cm[[2L]]), fam))
  if (length(terms) == 1L) return(terms[[1L]])
  key <- paste0(F[["add"]], "{", paste(sort(vapply(terms, .k, ""), method = "radix"), collapse = ","), "}")
  .as_symvec(new_sym(F[["add"]], terms, extra = list(components = comps), key = key))
}

register_rebuilder("VectorMul", function(e, args) .bd_mul(args, "V"))
register_rebuilder("DyadicMul", function(e, args) .bd_mul(args, "D"))
register_rebuilder("VectorAdd", function(e, args) .bd_add(args, "V"))
register_rebuilder("DyadicAdd", function(e, args) .bd_add(args, "D"))
register_rebuilder("BaseDyadic", function(e, args) BaseDyadic(args[[1L]], args[[2L]]))

#' Vectors and dyadics of a coordinate system
#'
#' `VectorAdd()`, `VectorMul()`, `DyadicAdd()` and `DyadicMul()` build sums
#' and scalar multiples of vectors and dyadics, as the arithmetic operators
#' do. `VectorZero()` and `DyadicZero()` are the zero vector and the zero
#' dyadic (`Vector.zero` and `Dyadic.zero` in SymPy). `BaseVector(index,
#' system)` is a basis vector of a [CoordSys3D()] and `BaseDyadic(v1, v2)` the
#' outer product of two basis vectors. `measure_number()` and `base_vector()`
#' give the parts of a `VectorMul`; `components()` the components of a vector
#' or dyadic as a named list.
#'
#' Vectors support `+`, `-`, multiplication and division by scalars, `&` for
#' the dot product and `|` for the outer product; `^` raises an error as
#' SymPy's `**` does. Use [cross()] for the cross product.
#'
#' @param ... Vectors or dyadics, and scalars.
#' @param index An integer 0, 1 or 2.
#' @param system A [CoordSys3D()].
#' @param vector1,vector2 Base vectors.
#' @param v A vector or dyadic.
#' @return A vector or dyadic.
#' @examples
#' C <- CoordSys3D("C")
#' v <- 3 * C$i + C$x * C$j
#' v
#' VectorAdd(C$i, C$j, C$i)
#' VectorMul(2, C$k)
#' C$i | C$j
#' @export
VectorAdd <- function(...) .bd_add(list(...), "V")

#' @rdname VectorAdd
#' @export
VectorMul <- function(...) .bd_mul(list(...), "V")

#' @rdname VectorAdd
#' @export
VectorZero <- function() .VZERO

#' @rdname VectorAdd
#' @export
DyadicAdd <- function(...) .bd_add(list(...), "D")

#' @rdname VectorAdd
#' @export
DyadicMul <- function(...) .bd_mul(list(...), "D")

#' @rdname VectorAdd
#' @export
DyadicZero <- function() .DZERO

#' @rdname VectorAdd
#' @export
BaseVector <- function(index, system) {
  if (!(as.numeric(index) %in% 0:2)) .vec_stop("ValueError", "index must be 0, 1 or 2")
  if (!.is_coordsys(system)) .vec_stop("TypeError", "system should be a CoordSys3D")
  .cs(system)$base_vectors[[as.integer(index) + 1L]]
}

.new_base_vector <- function(index, system, vname) {
  .as_symvec(new_sym("BaseVector", name = paste0(system$name, ".", vname),
                     extra = list(system = system, index = index),
                     key = paste0("BaseVector[", index, "|", system$key, "]")))
}

#' @rdname VectorAdd
#' @export
BaseDyadic <- function(vector1, vector2) {
  ok <- function(v) is_sym(v) && v$type %in% c("BaseVector", "VectorZero")
  if (!ok(vector1) || !ok(vector2))
    .vec_stop("TypeError", "BaseDyadic cannot be composed of non-base vectors")
  if (vector1$type == "VectorZero" || vector2$type == "VectorZero") return(.DZERO)
  .as_symvec(new_sym("BaseDyadic", list(vector1, vector2)))
}

#' @rdname VectorAdd
#' @export
measure_number <- function(v) .bd_measure(v)

.vec_components <- function(v) .pydict(lapply(.bd_components(v), `[[`, 1L), lapply(.bd_components(v), `[[`, 2L))

# A Python dictionary for printing: keys and values in insertion order.
.pydict <- function(keys, values) structure(list(keys = keys, values = values), class = "pydict")

#' @export
format.pydict <- function(x, ...) {
  if (!length(x$keys)) return("{}")
  paste0("{", paste(vapply(seq_along(x$keys), function(i)
    paste0(.pyfmt(x$keys[[i]]), ": ", .pyfmt(x$values[[i]])), ""), collapse = ", "), "}")
}

#' @export
print.pydict <- function(x, ...) {
  cat(format(x), "\n")
  invisible(x)
}

# Python's str() of a value: sym, matrix, list, tuple, dict, bool or None.
.pyfmt <- function(v) {
  if (is.null(v)) return("None")
  if (inherits(v, "pydict")) return(format(v))
  if (inherits(v, "pytuple")) {
    s <- vapply(v, .pyfmt, "")
    return(paste0("(", paste(s, collapse = ", "), if (length(s) == 1L) ",", ")"))
  }
  if (is_sym(v) || inherits(v, "symmatrix")) return(sstr(v))
  if (is.logical(v) && length(v) == 1L) return(if (is.na(v)) "None" else if (v) "True" else "False")
  if (is.character(v) && length(v) == 1L) return(v)
  if (is.numeric(v) && length(v) == 1L) return(as.character(v))
  if (is.list(v) && !is.object(v)) return(paste0("[", paste(vapply(v, .pyfmt, ""), collapse = ", "), "]"))
  paste(format(v), collapse = " ")
}

.pytuple <- function(x) structure(as.list(x), class = "pytuple")

#' @export
format.pytuple <- function(x, ...) .pyfmt(x)

#' @export
print.pytuple <- function(x, ...) {
  cat(format(x), "\n")
  invisible(x)
}

# ---------------------------------------------------------------------------
# Printing

.vmul_str <- function(e) {
  ms <- .str(e$extra$measure)
  if (grepl("[(+-]", ms)) ms <- paste0("(", ms, ")")
  paste0(ms, "*", .str(e$extra$base))
}

register_str_printer("BaseVector", function(e) e$name)
register_str_printer("VectorZero", function(e) "0")
register_str_printer("DyadicZero", function(e) "0")
register_str_printer("VectorMul", .vmul_str)
register_str_printer("DyadicMul", .vmul_str)
register_str_printer("BaseDyadic", function(e) paste0("(", .str(e$args[[1L]]), "|", .str(e$args[[2L]]), ")"))

# VectorAdd._sympystr: systems sorted by name, then the base vectors in order
register_str_printer("VectorAdd", function(e) {
  comps <- .bd_components(e)
  ck <- vapply(comps, function(cm) cm[[1L]]$key, "")
  parts <- .vec_separate(e)
  out <- character()
  for (i in .codepoint_order(vapply(parts$keys, .str, ""))) {
    vk <- vapply(.bd_components(parts$values[[i]]), function(cm) cm[[1L]]$key, "")
    for (bv in .cs(parts$keys[[i]])$base_vectors) if (bv$key %in% vk)
      out <- c(out, .str(.bd_mul(list(comps[[match(bv$key, ck)]][[2L]], bv), "V")))
  }
  paste(out, collapse = " + ")
})

# DyadicAdd._sympystr: components sorted by the string of the base dyadic
register_str_printer("DyadicAdd", function(e) {
  comps <- .bd_components(e)
  o <- .codepoint_order(vapply(comps, function(cm) .str(cm[[1L]]), ""))
  paste(vapply(comps[o], function(cm) .str(.bd_mul(list(cm[[1L]], cm[[2L]]), "D")), ""), collapse = " + ")
})

# stable sort of strings by code point, as Python's list.sort
.codepoint_order <- function(s) order(s, seq_along(s), method = "radix")

# ---------------------------------------------------------------------------
# SymPy's default_sort_key for expressions that hold vectors



.vsort_key <- function(e) {
  e <- as_sym(e)
  t <- e$type
  if (t %in% c("BaseVector", "BaseDyadic"))
    return(list(list(2, 0, t), list(1, list(.str(e))), .ONE_SORT_KEY, 1))
  if (t %in% c("VectorZero", "DyadicZero"))
    return(list(list(5, 0, t), list(0, list()), .ONE_SORT_KEY, 1))
  if (t %in% c("VectorMul", "DyadicMul")) {
    ks <- .key_sort(lapply(e$args, .vsort_key))
    return(list(list(3, 0, t), list(length(ks), ks), .ONE_SORT_KEY, 1))
  }
  if (t %in% c("VectorAdd", "DyadicAdd")) {
    ks <- lapply(.bd_ordered_terms(e), .vsort_key)
    return(list(list(3, 1, t), list(length(ks), ks), .ONE_SORT_KEY, 1))
  }
  if (!.has_type(e, c(.BD_TYPES, "Cross", "Dot", "Gradient", "Divergence", "Curl", "Laplacian"))) return(.sort_key(e))
  ks <- lapply(e$args, .vsort_key)
  if (t == "Mul") ks <- .key_sort(ks)
  list(.class_key(e), list(length(ks), ks), .ONE_SORT_KEY, 1)
}

# Expr.as_ordered_terms of a vector sum: each base vector is a generator of
# its term, ordered with the scalars by default_sort_key.
.bd_ordered_terms <- function(e) {
  terms <- e$args
  proxies <- lapply(terms, function(tm) {
    b <- .bd_base(tm)
    nm <- .str(b)
    px <- new_sym("Symbol", name = nm, extra = list(facts = .symbol_facts(logical()), given = logical(),
                                                     sympy_class = b$type),
                  key = paste0("S", nchar(nm), ":", nm, "{proxy:", b$key, "}"))
    .mul(list(.bd_measure(tm), px))
  })
  ord <- .add_print_order(proxies)
  terms[match(vapply(ord, .k, ""), vapply(proxies, .k, ""))]
}

# ---------------------------------------------------------------------------
# Arithmetic: BasisDependent.__add__, __sub__, __mul__, __truediv__, __neg__

.bd_neg <- function(e) {
  if (.is_bd_zero(e)) return(e)
  .bd_mul(list(S_NegativeOne, e), .bd_fam(e))
}

.bd_div <- function(one, other) {
  if (.is_dya(one) || .is_dya(other)) {
    if (.is_dya(one) && .is_dya(other)) .vec_stop("TypeError", "Cannot divide two dyadics")
    if (.is_dya(one)) return(.bd_mul(list(one, .pow(as_sym(other), S_NegativeOne)), "D"))
    .vec_stop("TypeError", "Invalid divisor for division")
  }
  if (.is_vec(one) && .is_vec(other)) .vec_stop("TypeError", "Cannot divide two vectors")
  if (.is_vec(one)) {
    other <- as_sym(other)
    if (is_zero_num(other)) .vec_stop("ValueError", "Cannot divide a vector by zero")
    return(.bd_mul(list(one, .pow(other, S_NegativeOne)), "V"))
  }
  .vec_stop("TypeError", "Invalid divisor for division")
}

# the operand whose method Python calls: the higher _op_priority wins
# (Vector 12, VectorZero 12.1, Dyadic 13, DyadicZero 13.1)
.bd_prio <- function(e) {
  if (!.is_bd(e)) return(0)
  switch(e$type, VectorZero = 12.1, DyadicZero = 13.1, if (.is_vec(e)) 12 else 13)
}

.bd_add2 <- function(a, b, sub = FALSE) {
  # BasisDependentZero.__add__ and __radd__ (and __sub__, __rsub__)
  self_left <- .bd_prio(a) >= .bd_prio(b)
  self <- if (self_left) a else b
  other <- if (self_left) b else a
  if (.is_bd_zero(self)) {
    if (!.is_bd(other) || .bd_fam(other) != .bd_fam(self))
      .vec_stop("TypeError", if (sub) "Invalid argument types for subtraction" else "Invalid argument types for addition")
    if (sub && self_left) return(.bd_neg(other))
    return(other)
  }
  fam <- .bd_fam(self)
  if (sub) b <- if (.is_bd(b)) .bd_neg(b) else .neg(as_sym(b))
  .bd_add(list(a, b), fam)
}

.vec_ops <- function(gen, e1, e2) {
  if (inherits(e1, "Del") || inherits(e2, "Del")) return(.del_ops(gen, e1, e2))
  b1 <- is_sym(e1) && .is_bd(e1)
  b2 <- is_sym(e2) && .is_bd(e2)
  if (!b1 && !b2) return(NULL)
  if (gen %in% c("==", "!=")) return(NULL)
  a <- as_sym(e1)
  b <- as_sym(e2)
  switch(gen,
         "+" = .bd_add2(a, b),
         "-" = .bd_add2(a, b, sub = TRUE),
         "*" = {
           self <- if (.bd_prio(a) >= .bd_prio(b)) a else b
           .bd_mul(list(a, b), .bd_fam(self))
         },
         "/" = .bd_div(a, b),
         "^" = .vec_stop("TypeError", "Power operation is not supported for vectors"),
         "&" = .bd_dot_method(a, b),
         "|" = {
           if (!.is_vec(a)) .vec_stop("TypeError", paste0("unsupported operand type(s) for |"))
           .vec_outer(a, b)
         },
         .vec_stop("TypeError", paste0("unsupported operand for vectors: ", gen)))
}

local({
  old <- .ops_hooks$dispatch
  .ops_hooks$dispatch <- function(gen, e1, e2) {
    r <- .vec_ops(gen, e1, e2)
    if (!is.null(r)) return(r)
    if (is.null(old)) NULL else old(gen, e1, e2)
  }
})

#' @export
Ops.sym_vec <- function(e1, e2) {
  if (missing(e2)) {
    if (.Generic == "-") return(.bd_neg(e1))
    if (.Generic == "+") return(e1)
    .vec_stop("TypeError", paste0("bad operand type for unary ", .Generic))
  }
  r <- .vec_ops(.Generic, e1, e2)
  if (!is.null(r)) return(r)
  NextMethod()
}

#' @export
diff.sym_vec <- function(x, ..., evaluate = TRUE) {
  vars <- list(...)
  for (v in vars) if (.is_bd(v)) .vec_stop("TypeError", "Invalid arg for differentiation")
  .bd_map(x, function(m) do.call(diff, c(list(m), vars, list(evaluate = evaluate))))
}

# BasisDependent.simplify, trigsimp, factor, evalf and diff: act on each
# component and add the parts back up.
.bd_map <- function(e, f) {
  fam <- .bd_fam(e)
  .bd_add(lapply(.bd_components(e), function(cm) .bd_mul(list(f(cm[[2L]]), cm[[1L]]), fam)), fam)
}

# Basic.conjugate of a vector conjugates each measure number
local({
  d <- .FUNS[["conjugate"]]
  old <- d$eval
  d$eval <- function(args) {
    a <- args[[1L]]
    if (.is_bd(a)) {
      if (.is_bd_zero(a)) return(a)
      return(.bd_add(lapply(.bd_components(a), function(cm) .mul_any(fcall("conjugate", list(cm[[2L]])), cm[[1L]])), .bd_fam(a)))
    }
    old(args)
  }
  assign("conjugate", d, envir = .FUNS)
})

register_diff_hook("VectorAdd", function(e, x) .bd_map(e, function(m) .diff1(m, x)))
register_diff_hook("VectorMul", function(e, x) .bd_map(e, function(m) .diff1(m, x)))
register_diff_hook("DyadicAdd", function(e, x) .bd_map(e, function(m) .diff1(m, x)))
register_diff_hook("DyadicMul", function(e, x) .bd_map(e, function(m) .diff1(m, x)))
register_diff_hook("BaseVector", function(e, x) .VZERO)
register_diff_hook("BaseDyadic", function(e, x) .DZERO)
register_diff_hook("VectorZero", function(e, x) .VZERO)
register_diff_hook("DyadicZero", function(e, x) .DZERO)

# ---------------------------------------------------------------------------
# Products: Vector.dot, Vector.cross, Vector.outer and the module functions

# the method of the left operand: Vector.dot, Dyadic.dot, or Del
.bd_dot_method <- function(a, b) {
  if (.is_dya(a)) return(.dya_dot(a, b))
  if (.is_vec(a)) return(.vec_dot_method(a, b))
  .vec_stop("TypeError", paste0("unsupported operand type(s) for &"))
}

.vec_dot_method <- function(self, other) {
  if (.is_dya(other)) {
    if (self$type == "VectorZero") return(.VZERO)
    out <- .VZERO
    for (cm in .bd_components(other)) {
      vd <- .vec_dot_method(cm[[1L]]$args[[1L]], self)
      out <- .bd_add2(out, .bd_mul(list(vd, cm[[2L]], cm[[1L]]$args[[2L]]), "V"))
    }
    return(out)
  }
  if (inherits(other, "Del")) {
    v <- self
    return(function(field) directional_derivative(field, v))
  }
  if (!.is_vec(other)) .vec_stop("TypeError", paste0(.pystr(other), " is not a vector, dyadic or del operator"))
  .vec_dot(self, other)
}

.pystr <- function(v) if (is_sym(v)) .str(v) else .pyfmt(v)

.dya_dot <- function(self, other) {
  if (.is_bd_zero(other)) return(.VZERO)
  if (.is_vec(other)) {
    out <- .VZERO
    for (cm in .bd_components(self)) {
      vd <- .vec_dot_method(cm[[1L]]$args[[2L]], other)
      out <- .bd_add2(out, .bd_mul(list(vd, cm[[2L]], cm[[1L]]$args[[1L]]), "V"))
    }
    return(out)
  }
  if (.is_dya(other)) {
    out <- .DZERO
    for (c1 in .bd_components(self)) for (c2 in .bd_components(other)) {
      vd <- .vec_dot_method(c1[[1L]]$args[[2L]], c2[[1L]]$args[[1L]])
      op <- .vec_outer(c1[[1L]]$args[[1L]], c2[[1L]]$args[[2L]])
      out <- .bd_add2(out, .bd_mul(list(.mul(list(vd, c1[[2L]], c2[[2L]])), op), "D"))
    }
    return(out)
  }
  .vec_stop("TypeError", paste0("Inner product is not defined for ", .vec_py_class(other), " and Dyadics."))
}

.vec_py_class <- function(v) {
  if (is.numeric(v)) return(if (v == round(v)) "<class 'int'>" else "<class 'float'>")
  if (is_sym(v)) return(paste0("<class 'sympy.core.numbers.", v$type, "'>"))
  paste0("<class '", class(v)[1L], "'>")
}

.vec_cross_method <- function(self, other) {
  if (.is_dya(self)) {
    if (.is_vec(other) && other$type == "VectorZero") return(.DZERO)
    if (!.is_vec(other)) .vec_stop("TypeError", paste0(.vec_py_class(other), " not supported for cross with dyadics"))
    out <- .DZERO
    for (cm in .bd_components(self)) {
      cp <- .vec_cross_method(cm[[1L]]$args[[2L]], other)
      out <- .bd_add2(out, .bd_mul(list(cm[[2L]], .vec_outer(cm[[1L]]$args[[1L]], cp)), "D"))
    }
    return(out)
  }
  if (.is_dya(other)) {
    if (self$type == "VectorZero") return(.DZERO)
    out <- .DZERO
    for (cm in .bd_components(other)) {
      cp <- .vec_cross_method(self, cm[[1L]]$args[[1L]])
      out <- .bd_add2(out, .bd_mul(list(cm[[2L]], .vec_outer(cp, cm[[1L]]$args[[2L]])), "D"))
    }
    return(out)
  }
  .vec_cross(self, other)
}

.vec_outer <- function(self, other) {
  if (!.is_vec(other)) .vec_stop("TypeError", "Invalid operand for outer product")
  if (self$type == "VectorZero" || other$type == "VectorZero") return(.DZERO)
  args <- list()
  for (c1 in .bd_components(self)) for (c2 in .bd_components(other))
    args[[length(args) + 1L]] <- .bd_mul(list(.mul(list(c1[[2L]], c2[[2L]])), BaseDyadic(c1[[1L]], c2[[1L]])), "D")
  .bd_add(args, "D")
}

.is_add_like <- function(e) is_sym(e) && e$type %in% c("Add", "VectorAdd", "DyadicAdd")

# sympy.vector.vector.dot
.vec_dot <- function(v1, v2) {
  v1 <- as_sym(v1); v2 <- as_sym(v2)
  if (.is_add_like(v1)) return(.add(lapply(v1$args, function(i) .vec_dot(i, v2))))
  if (.is_add_like(v2)) return(.add(lapply(v2$args, function(i) .vec_dot(v1, i))))
  if (v1$type == "BaseVector" && v2$type == "BaseVector") {
    if (identical(v1$extra$system$key, v2$extra$system$key)) return(if (identical(v1$key, v2$key)) S_One else S_Zero)
    v <- tryCatch(express(v2, v1$extra$system), ValueError = function(err) NULL)
    if (is.null(v)) return(Dot(v1, v2))
    return(.vec_dot(v1, v))
  }
  if (v1$type == "VectorZero" || v2$type == "VectorZero") return(S_Zero)
  if (v1$type == "VectorMul") return(.mul_any(v1$extra$measure, .vec_dot(v1$extra$base, v2)))
  if (v2$type == "VectorMul") return(.mul_any(v2$extra$measure, .vec_dot(v1, v2$extra$base)))
  Dot(v1, v2)
}

# sympy.vector.vector.cross
.vec_cross <- function(v1, v2) {
  v1 <- as_sym(v1); v2 <- as_sym(v2)
  if (.is_add_like(v1)) return(.bd_add(lapply(v1$args, function(i) .vec_cross(i, v2)), "V"))
  if (.is_add_like(v2)) return(.bd_add(lapply(v2$args, function(i) .vec_cross(v1, i)), "V"))
  if (v1$type == "BaseVector" && v2$type == "BaseVector") {
    if (identical(v1$extra$system$key, v2$extra$system$key)) {
      n1 <- v1$extra$index; n2 <- v2$extra$index
      if (n1 == n2) return(.VZERO)
      n3 <- setdiff(0:2, c(n1, n2))
      sgn <- if ((n1 + 1L) %% 3L == n2) S_One else S_NegativeOne
      return(.bd_mul(list(sgn, .cs(v1$extra$system)$base_vectors[[n3 + 1L]]), "V"))
    }
    v <- tryCatch(express(v1, v2$extra$system), ValueError = function(err) NULL)
    if (is.null(v)) return(Cross(v1, v2))
    return(.vec_cross(v, v2))
  }
  if (v1$type == "VectorZero" || v2$type == "VectorZero") return(.VZERO)
  if (v1$type == "VectorMul") return(.mul_any(v1$extra$measure, .vec_cross(v1$extra$base, v2)))
  if (v2$type == "VectorMul") return(.mul_any(v2$extra$measure, .vec_cross(v1, v2$extra$base)))
  Cross(v1, v2)
}

# Python's a*b for a scalar a and a scalar or basis dependent b
.mul_any <- function(a, b) {
  a <- as_sym(a); b <- as_sym(b)
  if (.is_bd(a) || .is_bd(b)) {
    self <- if (.bd_prio(a) >= .bd_prio(b)) a else b
    return(.bd_mul(list(a, b), .bd_fam(self)))
  }
  .mul(list(a, b))
}

.add_any <- function(a, b) {
  a <- as_sym(a); b <- as_sym(b)
  if (.is_bd(a) || .is_bd(b)) return(.bd_add2(a, b))
  .add(list(a, b))
}

# ---------------------------------------------------------------------------
# Cross and Dot (unevaluated products)

#' Unevaluated dot and cross products
#'
#' `Dot(v1, v2)` and `Cross(v1, v2)` represent products without computing
#' them; [doit()] evaluates. As in SymPy, `Dot` orders its arguments and
#' `Cross(v2, v1)` becomes `-Cross(v1, v2)` when `v1` sorts first.
#'
#' @param expr1,expr2 Vectors.
#' @return A `sym` expression.
#' @examples
#' C <- CoordSys3D("C")
#' Dot(C$i + C$j, C$k)
#' doit(Cross(C$i, C$j))
#' @export
Cross <- function(expr1, expr2) {
  expr1 <- as_sym(expr1); expr2 <- as_sym(expr2)
  if (.key_cmp(.vsort_key(expr1), .vsort_key(expr2)) > 0) return(.neg(Cross(expr2, expr1)))
  new_sym("Cross", list(expr1, expr2))
}

#' @rdname Cross
#' @export
Dot <- function(expr1, expr2) {
  expr1 <- as_sym(expr1); expr2 <- as_sym(expr2)
  if (.key_cmp(.vsort_key(expr1), .vsort_key(expr2)) > 0) { tmp <- expr1; expr1 <- expr2; expr2 <- tmp }
  new_sym("Dot", list(expr1, expr2))
}

register_rebuilder("Cross", function(e, args) Cross(args[[1L]], args[[2L]]))
register_rebuilder("Dot", function(e, args) Dot(args[[1L]], args[[2L]]))
register_doit("Cross", function(e) .vec_cross(e$args[[1L]], e$args[[2L]]))
register_doit("Dot", function(e) .vec_dot(e$args[[1L]], e$args[[2L]]))

# ---------------------------------------------------------------------------
# Vector methods

#' Vector operations
#'
#' `magnitude()`, `normalize()`, `projection()`, `to_matrix()`,
#' `separate()` and `vector_equals()` port the methods of SymPy's `Vector`.
#' `outer_sym()` is the outer product (also `v1 | v2`). [dot()] and [cross()]
#' accept vectors and dyadics as well as matrices.
#'
#' @param v,other Vectors.
#' @param system,second_system Coordinate systems ([CoordSys3D()]).
#' @return An expression, vector, dyadic, matrix or list.
#' @examples
#' C <- CoordSys3D("C")
#' v <- 3 * C$i + 4 * C$j
#' magnitude(v)
#' normalize(v)
#' to_matrix(v, C)
#' outer_sym(C$i, C$j)
#' @export
magnitude <- function(v) .sqrt(.vec_dot_method(v, v))

#' @rdname magnitude
#' @export
normalize <- function(v) {
  if (v$type == "VectorZero") return(v)
  .bd_div(v, magnitude(v))
}

#' @rdname magnitude
#' @export
outer_sym <- function(v, other) if (.pv_dispatch(v, other)) .pv_outer_dispatch(v, other) else .vec_outer(v, other)

#' @rdname magnitude
#' @export
vector_equals <- function(v, other) {
  d <- .bd_add2(v, other, sub = TRUE)
  .expr_equals(.vec_dot_method(d, d), S_Zero)
}

#' @rdname magnitude
#' @export
to_matrix <- function(v, system, second_system = NULL) {
  if (.is_dya(v)) {
    if (is.null(second_system)) second_system <- system
    d <- list()
    for (i in .cs(system)$base_vectors) for (j in .cs(second_system)$base_vectors)
      d[[length(d) + 1L]] <- .vec_dot_method(.vec_dot_method(i, v), j)
    return(.new_mat(d, 3L, 3L))
  }
  Matrix(lapply(.cs(system)$base_vectors, function(u) .vec_dot_method(v, u)))
}

#' @rdname magnitude
#' @export
separate <- function(v) .vec_separate(v)

.vec_separate <- function(v) {
  keys <- list(); vals <- list()
  for (cm in .bd_components(v)) {
    sys <- cm[[1L]]$extra$system
    i <- match(sys$key, vapply(keys, .k, ""))
    term <- .bd_mul(list(cm[[1L]], cm[[2L]]), "V")
    if (is.na(i)) {
      keys[[length(keys) + 1L]] <- sys
      vals[[length(vals) + 1L]] <- .bd_add2(.VZERO, term)
    } else {
      vals[[i]] <- .bd_add2(vals[[i]], term)
    }
  }
  .pydict(keys, vals)
}

.vec_projection <- function(v, other, scalar = FALSE) {
  if (isTRUE(vector_equals(v, .VZERO))) return(if (scalar) S_Zero else .VZERO)
  r <- .div(.vec_dot_method(v, other), .vec_dot_method(v, v))
  if (scalar) r else .bd_mul(list(r, v), "V")
}

# Vector._projections
.vec_projections <- function(v) {
  if (v$type == "VectorZero") return(list(S_Zero, S_Zero, S_Zero))
  sys <- .get_coord_systems(v)[[1L]]
  lapply(.cs(sys)$base_vectors, function(i) .vec_dot_method(v, i))
}

# Expr.equals(other): TRUE, FALSE or NULL when undecided.
.expr_equals <- function(a, b) {
  a <- as_sym(a); b <- as_sym(b)
  if (identical(a$key, b$key)) return(TRUE)
  d <- factor_terms(simplify(.sub(a, b)))
  if (is_zero_num(d)) return(TRUE)
  fs <- free_symbols(d)
  pts <- if (length(fs)) {
    lapply(1:2, function(k) stats::setNames(as.list(c(0.3137, 0.7219)[k] + seq_along(fs) / 7), vapply(fs, .k, "")))
  } else list(list())
  vals <- vapply(pts, function(p) {
    z <- if (length(fs)) xreplace(d, lapply(seq_along(fs), function(i) list(fs[[i]], Float(p[[i]])))) else d
    v <- tryCatch(as.complex(.neval(z, list())), error = function(err) NA_complex_)
    if (length(v) != 1L) NA_complex_ else v
  }, 0i)
  if (any(is.na(vals))) return(NULL)
  if (any(Mod(vals) > 1e-9)) return(FALSE)
  if (!length(fs)) return(TRUE)
  NULL
}

# ---------------------------------------------------------------------------
# Componentwise simplification (BasisDependent.simplify and friends)

.bd_simplify <- function(e, ...) {
  if (.is_bd_zero(e)) return(e)
  .bd_add(lapply(.bd_components(e), function(cm) .mul_any(simplify(cm[[2L]], ...), cm[[1L]])), .bd_fam(e))
}

.bd_trigsimp <- function(e, ...) {
  if (.is_bd_zero(e)) return(e)
  .bd_add(lapply(.bd_components(e), function(cm) .mul_any(trigsimp(cm[[2L]], ...), cm[[1L]])), .bd_fam(e))
}

.bd_factor <- function(e, ...) {
  if (.is_bd_zero(e)) return(e)
  .bd_add(lapply(.bd_components(e), function(cm) .mul_any(factor(cm[[2L]], ...), cm[[1L]])), .bd_fam(e))
}

.bd_evalf <- function(e, ...) {
  out <- .bd_zero(.bd_fam(e))
  for (cm in .bd_components(e)) out <- .bd_add2(out, .mul_any(N(cm[[2L]], ...), cm[[1L]]))
  out
}

.bd_integral <- function(e, ...) {
  .bd_add(lapply(.bd_components(e), function(cm) .mul_any(Integral(cm[[2L]], ...), cm[[1L]])), .bd_fam(e))
}

# Expr.expand on a vector sum rebuilds it from its expanded terms; a
# VectorMul expands its measure number.
.bd_expand <- function(e, ...) {
  if (.is_bd_zero(e) || e$type %in% c("BaseVector", "BaseDyadic")) return(e)
  .bd_add(lapply(.bd_components(e), function(cm) .mul_any(expand(cm[[2L]], ...), cm[[1L]])), .bd_fam(e))
}
