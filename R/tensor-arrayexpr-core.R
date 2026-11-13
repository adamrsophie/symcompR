# Array expressions, after sympy/tensor/array/expressions/array_expressions.py:
# symbolic arrays (ArraySymbol), tensor products, sums, axis permutations,
# diagonals and contractions, with the canonicalization rules SymPy applies
# when an operation is built with canonicalize = TRUE. Axis indices follow
# SymPy and start at 0.

.AE_TYPES <- c("ArraySymbol", "ZeroArray", "OneArray", "ArrayTensorProduct", "ArrayAdd", "PermuteDims",
               "ArrayDiagonal", "ArrayElementwiseApplyFunc", "ArrayContraction", "Reshape", "ArraySum")

# isinstance(e, _ArrayExpr)
.is_arrexpr <- function(e) is_sym(e) && e$type %in% .AE_TYPES

.ae_new <- function(type, args, shape, sub = NULL, name = NULL)
  new_sym(type, args, name = name, extra = list(shape = shape, sub = sub))

# integer vector <-> Tuple of Integers
.ae_t <- function(v) new_sym("Tuple", lapply(as.integer(v), .int))
.ae_v <- function(t) {
  if (is.numeric(t)) return(as.integer(t))
  if (inherits(t, "Permutation")) return(t$array)
  if (is_sym(t) && t$type == "Tuple") return(vapply(t$args, function(a) as.integer(.int_val(a)), 1L))
  if (is.list(t)) return(vapply(t, function(a) if (is_sym(a)) as.integer(.int_val(a)) else as.integer(a), 1L))
  stop("expected a tuple of integers", call. = FALSE)
}

.ae_skey <- function(s) paste(vapply(s, .k, ""), collapse = ",")
.ae_shape_eq <- function(a, b) identical(.ae_skey(a), .ae_skey(b))
.ae_eq <- function(a, b) {
  if (is_sym(a) && is_sym(b)) return(identical(a$key, b$key))
  identical(a, b)
}

.ae_sympify <- function(a) {
  if (inherits(a, "symmatrix")) return(.as_block(a))
  if (inherits(a, "SymArray")) return(.ae_ndarray_sym(a))
  if (is_sym(a)) return(a)
  as_sym(a)
}

# An explicit array held inside an expression tree.
.ae_ndarray_sym <- function(a) new_sym("ImmutableDenseNDimArray", a$data, extra = list(shape = a$shape))
.ae_ndarray_r <- function(e) Array(e$args, e$extra$shape)

# commutative scalar
.ae_comm <- function(e) !identical(.facts(e)[["commutative"]], FALSE)

# ---------------------------------------------------------------------------
# Permutations stored inside expressions, printed in cyclic notation

.ae_perm_sym <- function(af) new_sym("Permutation", lapply(as.integer(af), .int))
.ae_paf <- function(p) if (is_sym(p)) vapply(p$args, function(a) as.integer(.int_val(a)), 1L) else as.integer(p)
register_str_printer("Permutation", function(e) .perm_str(.perm_af(.ae_paf(e))))
register_rebuilder("Permutation", function(e, args) e)

# SymPy's Permutation(p)*Permutation(q): apply p, then q
.ae_pmul <- function(a, b) {
  n <- max(length(a), length(b))
  if (length(a) < n) a <- c(a, seq.int(length(a), length.out = n - length(a)))
  if (length(b) < n) b <- c(b, seq.int(length(b), length.out = n - length(b)))
  b[a + 1L]
}
.ae_inv <- function(a) .af_invert(as.integer(a))
.ae_id <- function(n) seq_len(n) - 1L
# Permutation from disjoint cycles
.ae_cycles_af <- function(cycles, size = 0L) {
  m <- max(c(size, vapply(cycles, function(c) if (length(c)) max(c) + 1L else 0L, 0L)))
  af <- .ae_id(m)
  for (c in cycles) if (length(c) > 1L) af[c + 1L] <- c[c(seq_along(c)[-1L], 1L)]
  af
}
.ae_cyclic <- function(af) cyclic_form(.perm_af(af))

# Permutation(...) argument: an array form, a Permutation or a list of cycles
.ae_perm_input <- function(p) {
  if (inherits(p, "Permutation")) return(p$array)
  if (is_sym(p) && p$type == "Permutation") return(.ae_paf(p))
  if (is_sym(p) && p$type == "Tuple") return(.ae_v(p))
  if (is.list(p) && length(p) && all(vapply(p, function(z) length(z) != 1L || is.list(z), TRUE)))
    return(.ae_cycles_af(lapply(p, .ae_v)))
  if (is.list(p)) return(.ae_v(p))
  as.integer(p)
}

# ---------------------------------------------------------------------------
# shapes and ranks

#' Shape and rank of array and matrix expressions
#'
#' `get_shape()` returns the shape of an array expression, a matrix
#' expression or an explicit array, as a list of dimensions (an empty list
#' for a scalar). `get_ndim()` returns the number of axes, and `get_rank()`
#' is its older name.
#'
#' @param expr An expression.
#' @return `get_shape()` returns a list of `sym` dimensions; `get_ndim()`
#'   and `get_rank()` return an integer.
#' @examples
#' A <- ArraySymbol("A", c(3, 3, 4))
#' get_shape(A)
#' get_ndim(ArrayContraction(A, c(0, 1)))
#' @export
get_shape <- function(expr) .ae_shape(expr)

#' @rdname get_shape
#' @export
get_ndim <- function(expr) {
  if (inherits(expr, "symmatrix") || .is_matexpr(expr)) return(2L)
  if (inherits(expr, "SymArray")) return(length(expr$shape))
  if (!is_sym(expr)) return(0L)
  if (expr$type == "MatrixElement") return(0L)
  if (expr$type == "Indexed") return(length(expr$args) - 1L)
  length(.ae_shape(expr))
}

#' @rdname get_shape
#' @export
get_rank <- function(expr) get_ndim(expr)

.ae_shape <- function(e) {
  if (inherits(e, "symmatrix")) return(list(.int(e$nrow), .int(e$ncol)))
  if (inherits(e, "SymArray")) return(lapply(e$shape, .int))
  if (!is_sym(e)) return(list())
  if (e$type %in% c(.AE_TYPES, "ArrayComprehension", "ArrayComprehensionMap")) return(e$extra$shape)
  if (e$type == "ImmutableDenseNDimArray") return(lapply(e$extra$shape, .int))
  if (e$type == "MatrixSet") return(e$args[1:2])
  if (.is_matexpr(e)) return(.mshape(e))
  if (e$type == "Derivative" && .ae_derivative_is_array(e)) return(.ae_derivative_shape(e))
  list()
}

.ae_sub_ndim_list <- function(e) {
  if (is_sym(e) && !is.null(e$extra$sub)) return(e$extra$sub)
  get_ndim(e)
}
.ae_sub_ndim <- function(e) sum(.ae_sub_ndim_list(e))

.ae_is_plain_scalar <- function(a) !length(.ae_shape(a)) && !.is_arrexpr(a)

# ---------------------------------------------------------------------------
# ArraySymbol, ArrayElement, ZeroArray, OneArray

#' Array expressions
#'
#' Symbolic arrays and the operations on them, following
#' `sympy.tensor.array.expressions`. `ArraySymbol()` is an array of a given
#' shape, `ArrayElement()` one of its entries, `ZeroArray()` and `OneArray()`
#' arrays of zeros and ones. `ArrayTensorProduct()`, `ArrayAdd()`,
#' `PermuteDims()`, `ArrayDiagonal()`, `ArrayContraction()`, `Reshape()`,
#' `ArrayElementwiseApplyFunc()` and `ArraySum()` build expressions; they keep
#' their arguments as given unless `canonicalize = TRUE`, as in SymPy. Axes
#' are numbered from 0, as in SymPy.
#'
#' @param name A name (character) or an expression.
#' @param shape A vector or list of dimensions.
#' @param indices A vector or list of indices.
#' @param ... Dimensions, arguments, or axis groups (integer vectors).
#' @param expr An array or matrix expression.
#' @param permutation A permutation: a `Permutation`, an array form or a list
#'   of cycles.
#' @param index_order_old,index_order_new Index orders defining the
#'   permutation, as an alternative to `permutation`.
#' @param canonicalize Put the result in canonical form.
#' @param allow_trivial_diags Allow diagonal groups of one axis.
#' @param fun A function or `Lambda` applied to every element.
#' @param element An array expression.
#' @param f An array expression to sum.
#' @return A `sym` array expression.
#' @examples
#' A <- ArraySymbol("A", c(3, 3))
#' B <- ArraySymbol("B", c(3, 3))
#' ArrayContraction(ArrayTensorProduct(A, B), c(1, 2))
#' PermuteDims(A, c(1, 0))
#' @name array-expressions
NULL

#' @rdname array-expressions
#' @export
ArraySymbol <- function(name, shape) {
  sym <- if (is.character(name)) Symbol(name) else as_sym(name)
  shape <- .ae_seq(shape)
  .ae_new("ArraySymbol", list(sym, new_sym("Tuple", shape)), shape)
}
register_str_printer("ArraySymbol", function(e) .str(e$args[[1L]]))
register_latex_printer("ArraySymbol", function(e) .latex(e$args[[1L]]))
register_rebuilder("ArraySymbol", function(e, args) ArraySymbol(args[[1L]], args[[2L]]$args))

.ae_rel_true <- function(d, fact) isTRUE(.facts(d)[[fact]]) || (is_exact_number(d) && switch(fact,
  nonnegative = d$value >= 0, negative = d$value < 0, positive = d$value > 0, FALSE))

.ae_check_shape <- function(name, indices) {
  sh <- .ae_shape(name)
  if (is_sym(name) && (length(sh) || .is_arrexpr(name) || .is_matexpr(name))) {
    if (length(indices) != length(sh))
      .perr("IndexError", "number of indices does not match shape of the array")
    # (i >= s) == True: is_ge decides only for extended real sides
    ge <- function(a, b) isTRUE(ask_fact(a, "extended_real")) && isTRUE(ask_fact(b, "extended_real")) &&
      .ae_rel_true(.sub(a, b), "nonnegative")
    if (any(vapply(seq_along(indices), function(i) ge(indices[[i]], sh[[i]]), NA)))
      .perr("ValueError", "shape is out of bounds")
  }
  if (any(vapply(indices, function(i) .ae_rel_true(i, "negative"), NA)))
    .perr("ValueError", "shape contains negative values")
}

#' @rdname array-expressions
#' @export
ArrayElement <- function(name, indices) {
  name <- if (is.character(name)) Symbol(name) else .ae_sympify(name)
  indices <- .ae_seq(indices)
  .ae_check_shape(name, indices)
  new_sym("ArrayElement", list(name, new_sym("Tuple", indices)))
}
# _ArrayExpr.__getitem__: ZeroArray and OneArray give their entries
.ae_getitem <- function(x, indices) {
  indices <- .ae_seq(indices)
  .ae_check_shape(x, indices)
  if (x$type == "ZeroArray") return(S_Zero)
  if (x$type == "OneArray") return(S_One)
  ArrayElement(x, indices)
}
register_str_printer("ArrayElement", function(e)
  paste0(.paren(e$args[[1L]], PREC$Func, TRUE), "[", paste(vapply(e$args[[2L]]$args, .str, ""), collapse = ", "), "]"))
register_rebuilder("ArrayElement", function(e, args) ArrayElement(args[[1L]], args[[2L]]$args))

# a Tuple, vector or list as a list of sym
.ae_seq <- function(x) {
  if (is_sym(x)) return(if (x$type == "Tuple") x$args else list(x))
  lapply(as.list(x), as_sym)
}

# shape arguments: dimensions, or vectors or lists of them
.ae_dims <- function(...) {
  out <- list()
  for (a in list(...)) {
    if (is_sym(a) && a$type != "Tuple") out[[length(out) + 1L]] <- a
    else out <- c(out, .ae_seq(a))
  }
  out
}

#' @rdname array-expressions
#' @export
ZeroArray <- function(...) {
  shape <- .ae_dims(...)
  if (!length(shape)) return(S_Zero)
  .ae_new("ZeroArray", shape, shape)
}
#' @rdname array-expressions
#' @export
OneArray <- function(...) {
  shape <- .ae_dims(...)
  if (!length(shape)) return(S_One)
  .ae_new("OneArray", shape, shape)
}
register_rebuilder("ZeroArray", function(e, args) do.call(ZeroArray, args))
register_rebuilder("OneArray", function(e, args) do.call(OneArray, args))

# ---------------------------------------------------------------------------
# scalar coefficients

.ae_as_coeff_matrices <- function(m) {
  sc <- Filter(function(a) !.is_matexpr(a), m$args)
  list(.mul(sc), Filter(.is_matexpr, m$args))
}

# _split_scalar_coefficient: list(coeff, array or NULL)
.ae_split_coeff <- function(arg) {
  if (.is_arrexpr(arg)) return(list(S_One, arg))
  if (is_sym(arg) && arg$type == "MatMul") {
    cm <- .ae_as_coeff_matrices(arg)
    if (is_one_num(cm[[1L]])) return(list(S_One, arg))
    rest <- if (length(cm[[2L]]) > 1L) do.call(MatMul, cm[[2L]]) else cm[[2L]][[1L]]
    return(list(cm[[1L]], rest))
  }
  if (is_sym(arg) && arg$type == "Mul") {
    sh <- vapply(arg$args, function(f) length(.ae_shape(f)) == 0L, NA)
    if (sum(!sh) == 1L) return(list(.mul(arg$args[sh]), arg$args[!sh][[1L]]))
    if (sum(!sh) > 1L) return(list(S_One, arg))
    return(list(arg, NULL))
  }
  if (!length(.ae_shape(arg))) return(list(arg, NULL))
  list(S_One, arg)
}

# _array_term_as_coeff_arrays: list(coeff, list of arrays)
.ae_term_coeff_arrays <- function(term) {
  if (is_sym(term) && term$type == "ArrayTensorProduct") {
    coeff <- S_One
    arrays <- list()
    for (a in term$args) {
      sp <- .ae_split_coeff(a)
      coeff <- .mul(list(coeff, sp[[1L]]))
      if (!is.null(sp[[2L]])) arrays[[length(arrays) + 1L]] <- sp[[2L]]
    }
    return(list(coeff, arrays))
  }
  sp <- .ae_split_coeff(term)
  if (is.null(sp[[2L]])) return(list(sp[[1L]], list()))
  list(sp[[1L]], list(sp[[2L]]))
}

.ae_term_from_coeff_arrays <- function(coeff, arrays) {
  if (!length(arrays)) return(coeff)
  if (is_one_num(coeff)) {
    if (length(arrays) == 1L) return(arrays[[1L]])
    return(do.call(.ae_tensor_product, arrays))
  }
  if (length(arrays) == 1L && .is_matexpr(arrays[[1L]]) && .ae_comm(coeff)) return(.ae_mmul(coeff, arrays[[1L]]))
  do.call(.ae_tensor_product, c(list(coeff), arrays))
}

# ---------------------------------------------------------------------------
# ArrayTensorProduct

#' @rdname array-expressions
#' @export
ArrayTensorProduct <- function(..., canonicalize = FALSE) {
  args <- lapply(list(...), .ae_sympify)
  norm <- list()
  for (a in args) {
    if (a$type == "Mul") {
      sp <- .ae_split_coeff(a)
      if (!is.null(sp[[2L]])) {
        if (!is_one_num(sp[[1L]])) norm[[length(norm) + 1L]] <- sp[[1L]]
        norm[[length(norm) + 1L]] <- sp[[2L]]
        next
      }
    }
    norm[[length(norm) + 1L]] <- a
  }
  args <- norm
  ndims <- vapply(args, get_ndim, 0L)
  shape <- do.call(c, c(list(list()), lapply(args, .ae_shape)))
  obj <- .ae_new("ArrayTensorProduct", args, shape, ndims)
  if (canonicalize) .ae_tp_canon(obj) else obj
}

.ae_tp_canon <- function(self) {
  args <- list()
  for (a in self$args) if (a$type == "ArrayTensorProduct") args <- c(args, a$args) else args[[length(args) + 1L]] <- a
  split <- list()
  for (a in args) {
    if (.is_matexpr(a)) {
      sp <- .ae_split_coeff(a)
      if (!is_one_num(sp[[1L]])) split[[length(split) + 1L]] <- sp[[1L]]
      split[[length(split) + 1L]] <- sp[[2L]]
    } else split[[length(split) + 1L]] <- a
  }
  args <- split
  isc <- vapply(args, .ae_is_plain_scalar, NA)
  if (any(isc)) {
    comm <- list(); nc <- list()
    for (a in args[isc]) {
      if (a$type == "Mul") {
        for (f in a$args) if (.ae_comm(f)) comm[[length(comm) + 1L]] <- f else nc[[length(nc) + 1L]] <- f
      } else if (!.ae_comm(a)) nc[[length(nc) + 1L]] <- a else comm[[length(comm) + 1L]] <- a
    }
    coeff <- .mul(comm)
    arrs <- args[!isc]
    if (is_zero_num(coeff)) return(do.call(ZeroArray, do.call(c, c(list(list()), lapply(arrs, .ae_shape)))))
    if (!length(arrs) && !length(nc)) return(coeff)
    args <- if (is_one_num(coeff)) c(nc, arrs) else c(list(coeff), nc, arrs)
  }
  ndim_list <- vapply(args, get_ndim, 0L)
  sums <- which(vapply(args, function(a) a$type == "ArraySum", NA))
  if (length(sums)) {
    new_limits <- list(); new_args <- list(); last <- 0L
    for (i in sums) {
      s <- args[[i]]
      new_args <- c(new_args, args[seq_len(i - 1L)][seq_len(i - 1L) > last])
      lims <- s$args[-1L]
      repl <- lapply(lims, function(l) Dummy(.str(l$args[[1L]])))
      new_limits <- c(new_limits, lapply(seq_along(lims), function(k) Tuple(repl[[k]], lims[[k]]$args[[2L]], lims[[k]]$args[[3L]])))
      new_args[[length(new_args) + 1L]] <- subs(s$args[[1L]], lapply(seq_along(lims), function(k) list(lims[[k]]$args[[1L]], repl[[k]])))
      last <- i
    }
    if (last < length(args)) new_args <- c(new_args, args[(last + 1L):length(args)])
    return(do.call(ArraySum, c(list(do.call(.ae_tensor_product, new_args)), new_limits)))
  }
  cycles <- list()
  for (i in seq_along(args)) {
    a <- args[[i]]
    if (a$type != "PermuteDims") next
    off <- sum(ndim_list[seq_len(i - 1L)])
    for (cy in .ae_cyclic(.ae_paf(a$args[[2L]]))) cycles[[length(cycles) + 1L]] <- cy + off
    args[[i]] <- a$args[[1L]]
  }
  if (length(cycles))
    return(.ae_permute_dims(do.call(.ae_tensor_product, args), .ae_pmul(.ae_id(sum(ndim_list)), .ae_cycles_af(cycles))))
  if (length(args) == 1L) return(args[[1L]])
  if (any(vapply(args, function(a) a$type %in% c("ZeroArray", "ZeroMatrix"), NA)))
    return(do.call(ZeroArray, do.call(c, c(list(list()), lapply(args, .ae_shape)))))
  isc <- vapply(args, function(a) a$type == "ArrayContraction", NA)
  if (any(isc)) {
    nd <- vapply(seq_along(args), function(i) if (isc[i]) .ae_sub_ndim(args[[i]]) else get_ndim(args[[i]]), 0L)
    cum <- cumsum(c(0L, nd))[seq_along(nd)]
    tp <- do.call(.ae_tensor_product, lapply(seq_along(args), function(i) if (isc[i]) args[[i]]$args[[1L]] else args[[i]]))
    ci <- list()
    for (i in which(isc)) for (g in .ae_contr_ind(args[[i]])) ci[[length(ci) + 1L]] <- cum[i] + g
    return(do.call(.ae_contraction, c(list(tp), ci)))
  }
  isd <- vapply(args, function(a) a$type == "ArrayDiagonal", NA)
  if (any(isd)) {
    inv <- integer(); last_perm <- integer()
    nd <- vapply(args, get_ndim, 0L)
    cum <- cumsum(c(0L, nd))[seq_along(nd)]
    for (i in seq_along(args)) {
      if (isd[i]) {
        i2 <- length(args[[i]]$args) - 1L
        i1 <- get_ndim(args[[i]]) - i2
        inv <- c(inv, cum[i] + seq_len(i1) - 1L)
        last_perm <- c(last_perm, cum[i] + i1 + seq_len(i2) - 1L)
      } else inv <- c(inv, cum[i] + seq_len(nd[i]) - 1L)
    }
    inv <- c(inv, last_perm)
    tp <- do.call(.ae_tensor_product, lapply(seq_along(args), function(i) if (isd[i]) args[[i]]$args[[1L]] else args[[i]]))
    nd2 <- vapply(seq_along(args), function(i) if (isd[i]) .ae_sub_ndim(args[[i]]) else get_ndim(args[[i]]), 0L)
    cum2 <- cumsum(c(0L, nd2))[seq_along(nd2)]
    di <- list()
    for (i in which(isd)) for (g in .ae_diag_ind(args[[i]])) di[[length(di) + 1L]] <- cum2[i] + g
    return(.ae_permute_dims(do.call(.ae_diagonal, c(list(tp), di)), .ae_inv(inv)))
  }
  do.call(ArrayTensorProduct, args)
}

# ---------------------------------------------------------------------------
# ArrayAdd

#' @rdname array-expressions
#' @export
ArrayAdd <- function(..., canonicalize = FALSE) {
  args <- lapply(list(...), .ae_sympify)
  args <- lapply(args, function(a) {
    if (a$type == "Mul") {
      sp <- .ae_split_coeff(a)
      if (!is.null(sp[[2L]])) return(ArrayTensorProduct(sp[[1L]], sp[[2L]]))
    }
    a
  })
  if (any(vapply(args, .ae_is_plain_scalar, NA))) .perr("TypeError", "scalar arguments are not allowed in ArrayAdd")
  if (length(unique(vapply(args, get_ndim, 0L))) != 1L) .perr("ValueError", "summing arrays of different number of dims")
  shapes <- lapply(args, .ae_shape)
  if (length(unique(vapply(shapes, .ae_skey, ""))) > 1L) .perr("ValueError", "mismatching shapes in addition")
  obj <- .ae_new("ArrayAdd", args, shapes[[1L]], get_ndim(args[[1L]]))
  if (canonicalize) .ae_add_canon(obj) else obj
}

.ae_add_canon <- function(self) {
  args <- list()
  for (a in self$args) if (a$type == "ArrayAdd") args <- c(args, a$args) else args[[length(args) + 1L]] <- a
  shapes <- lapply(args, .ae_shape)
  args <- Filter(function(a) !(a$type %in% c("ZeroArray", "ZeroMatrix")), args)
  if (length(args) > 1L) args <- .ae_collect_scalar_coefficients(args)
  if (!length(args)) return(do.call(ZeroArray, shapes[[1L]]))
  if (length(args) == 1L) return(args[[1L]])
  if (any(vapply(args, .ae_is_plain_scalar, NA))) return(.add(args))
  do.call(ArrayAdd, args)
}

.ae_collect_scalar_coefficients <- function(args) {
  keys <- character(); coeffs <- list(); arrs <- list()
  for (a in args) {
    ca <- .ae_term_coeff_arrays(a)
    k <- paste(vapply(ca[[2L]], .k, ""), collapse = "|")
    w <- match(k, keys)
    if (is.na(w)) {
      keys <- c(keys, k); coeffs[[length(coeffs) + 1L]] <- ca[[1L]]; arrs[[length(arrs) + 1L]] <- ca[[2L]]
    } else coeffs[[w]] <- .add(list(coeffs[[w]], ca[[1L]]))
  }
  if (all(vapply(coeffs, is_one_num, NA)) && length(keys) == length(args)) return(args)
  out <- list()
  for (i in seq_along(keys)) {
    if (is_zero_num(coeffs[[i]])) next
    out[[length(out) + 1L]] <- .ae_term_from_coeff_arrays(coeffs[[i]], arrs[[i]])
  }
  out
}

# ---------------------------------------------------------------------------
# PermuteDims

#' @rdname array-expressions
#' @export
PermuteDims <- function(expr, permutation = NULL, index_order_old = NULL, index_order_new = NULL, canonicalize = FALSE) {
  expr <- .ae_sympify(expr)
  nd <- get_ndim(expr)
  if (is.null(permutation)) {
    if (is.null(index_order_new) || is.null(index_order_old)) .perr("ValueError", "Permutation not defined")
    old <- vapply(as.list(index_order_old), function(z) as_sym(z)$key, "")
    new <- vapply(as.list(index_order_new), function(z) as_sym(z)$key, "")
    if (length(unique(new)) != nd) .perr("ValueError", "wrong number of indices in index_order_new")
    if (length(unique(old)) != nd) .perr("ValueError", "wrong number of indices in index_order_old")
    if (length(setdiff(union(new, old), intersect(new, old)))) .perr("ValueError", "index_order_new and index_order_old must have the same indices")
    af <- match(new, old) - 1L
  } else {
    if (!is.null(index_order_new)) .perr("ValueError", "index_order_new cannot be defined with permutation")
    if (!is.null(index_order_old)) .perr("ValueError", "index_order_old cannot be defined with permutation")
    af <- .ae_perm_input(permutation)
  }
  if (length(af) != nd) .perr("ValueError", "Permutation size must be the length of the shape of expr")
  shape <- .ae_shape(expr)[af + 1L]
  obj <- .ae_new("PermuteDims", list(expr, .ae_perm_sym(af)), shape, nd)
  if (canonicalize) .ae_pd_canon(obj) else obj
}

.ae_pd_canon <- function(self) {
  expr <- self$args[[1L]]
  perm <- .ae_paf(self$args[[2L]])
  if (expr$type == "PermuteDims") {
    perm <- .ae_pmul(perm, .ae_paf(expr$args[[2L]]))
    expr <- expr$args[[1L]]
  }
  if (expr$type == "ArrayContraction") {
    r <- .ae_pd_denest_contraction(expr, perm); expr <- r[[1L]]; perm <- r[[2L]]
  }
  if (expr$type == "ArrayTensorProduct") {
    r <- .ae_pd_denest_tp(expr, perm); expr <- r[[1L]]; perm <- r[[2L]]
  }
  if (expr$type %in% c("ZeroArray", "ZeroMatrix")) return(do.call(ZeroArray, .ae_shape(expr)[perm + 1L]))
  if (identical(perm, sort(perm))) return(expr)
  PermuteDims(expr, perm)
}

.ae_lexsort <- function(keys) {
  # order of a list of integer vectors, lexicographically (Python list sort)
  n <- length(keys)
  idx <- seq_len(n)
  cmp <- function(a, b) {
    for (i in seq_len(min(length(a), length(b)))) if (a[i] != b[i]) return(if (a[i] < b[i]) -1L else 1L)
    sign(length(a) - length(b))
  }
  if (n < 2L) return(idx)
  for (i in 2:n) {
    j <- i
    while (j > 1L && cmp(keys[[idx[j - 1L]]], keys[[idx[j]]]) > 0L) {
      tmp <- idx[j]; idx[j] <- idx[j - 1L]; idx[j - 1L] <- tmp; j <- j - 1L
    }
  }
  idx
}

.ae_pd_denest_tp <- function(expr, perm) {
  img <- .ae_inv(perm)
  args <- expr$args
  cum <- cumsum(c(0L, .ae_sub_ndim_list(expr)))
  comps <- lapply(seq_along(args), function(i) img[seq_len(cum[i + 1L] - cum[i]) + cum[i]])
  ord <- .ae_lexsort(lapply(comps, sort))
  new_perm <- .ae_inv(unlist(comps[ord]))
  list(do.call(.ae_tensor_product, args[ord]), as.integer(new_perm))
}

.ae_pd_denest_contraction <- function(expr, perm) {
  if (expr$type != "ArrayContraction") return(list(expr, perm))
  if (expr$args[[1L]]$type != "ArrayTensorProduct") return(list(expr, perm))
  inner <- expr$args[[1L]]
  args <- inner$args
  sub <- vapply(args, get_ndim, 0L)
  ci <- .ae_contr_ind(expr)
  flat <- unlist(ci)
  cum <- cumsum(c(0L, sub))
  img <- .ae_inv(perm)
  blocks_up <- list(); counter <- 0L
  for (i in seq_along(sub)) {
    cur <- integer()
    for (j in seq_len(sub[i]) + cum[i] - 1L) {
      if (j %in% flat) next
      cur <- c(cur, img[counter + 1L]); counter <- counter + 1L
    }
    blocks_up[[i]] <- cur
  }
  sl <- .ae_sub_ndim_list(expr)
  cs <- cumsum(c(0L, sl))
  index_blocks <- lapply(seq_along(sl), function(i) seq_len(sl[i]) + cs[i] - 1L)
  up <- .ae_contr_push_up(ci, lapply(index_blocks, as.list))
  invp <- .ae_inv(perm)
  permuted <- lapply(up, function(b) { b <- unlist(b); b <- b[!is.na(b)]; invp[b + 1L] })
  ord <- .ae_lexsort(permuted)
  new_index_blocks <- index_blocks[ord]
  nip <- .ae_inv(unlist(new_index_blocks))
  new_ci <- lapply(ci, function(g) nip[g + 1L])
  new_expr <- do.call(.ae_contraction, c(list(do.call(.ae_tensor_product, args[ord])), new_ci))
  new_perm <- .ae_inv(unlist(blocks_up[ord]))
  list(new_expr, as.integer(new_perm))
}

.ae_closed_cycles <- function(expr, perm) {
  args <- expr$args
  sub <- .ae_sub_ndim_list(expr)
  cy <- .ae_cyclic(perm)
  cum <- c(0L, cumsum(sub))
  keep <- list()
  for (c in cy) {
    flag <- TRUE
    for (j in seq_len(length(cum) - 1L)) {
      if (min(c) >= cum[j] && max(c) < cum[j + 1L]) {
        args[[j]] <- .ae_permute_dims(args[[j]], .ae_cycles_af(list(c - cum[j]), sub[j]))
        flag <- FALSE
        break
      }
    }
    if (flag) keep[[length(keep) + 1L]] <- c
  }
  list(do.call(.ae_tensor_product, args), .ae_cycles_af(keep, length(perm)))
}

#' Nest the permutation of a PermuteDims
#'
#' Moves the axis permutation of a `PermuteDims` expression inside its
#' argument where possible, as SymPy's `nest_permutation()`.
#'
#' @param expr An array expression.
#' @return An array expression.
#' @examples
#' M <- ArraySymbol("M", c(3, 3)); N <- ArraySymbol("N", c(3, 3))
#' nest_permutation(PermuteDims(ArrayTensorProduct(M, N), c(1, 0, 3, 2)))
#' @export
nest_permutation <- function(expr) {
  if (!is_sym(expr) || expr$type != "PermuteDims") return(expr)
  r <- .ae_nest_permutation(expr$args[[1L]], .ae_paf(expr$args[[2L]]))
  if (is.null(r)) expr else r
}

.ae_nest_permutation <- function(expr, perm) {
  if (expr$type == "ArrayTensorProduct") {
    r <- .ae_closed_cycles(expr, perm)
    return(.ae_permute_dims(r[[1L]], r[[2L]]))
  }
  if (expr$type == "ArrayContraction") {
    cycles <- .ae_cyclic(perm)
    newcycles <- .ae_outer_to_inner(expr, cycles)
    np <- .ae_cycles_af(newcycles)
    new_ci <- lapply(.ae_contr_ind(expr), function(g) np[g + 1L])
    return(do.call(.ae_contraction, c(list(PermuteDims(expr$args[[1L]], np)), new_ci)))
  }
  if (expr$type == "ArrayAdd") return(do.call(.ae_array_add, lapply(expr$args, function(a) PermuteDims(a, perm))))
  NULL
}

# ---------------------------------------------------------------------------
# nested index lists: an int is an integer scalar, a tuple an R list

.ae_rec <- function(f, arr) if (is.list(arr)) lapply(arr, function(z) .ae_rec(f, z)) else f(arr)

.ae_positions_shape <- function(shape, diag) {
  inany <- function(i) any(vapply(diag, function(g) i %in% g, NA))
  pos1 <- list(); shp1 <- list()
  for (i in seq_along(shape) - 1L) if (!inany(i)) { pos1[[length(pos1) + 1L]] <- i; shp1[[length(shp1) + 1L]] <- shape[[i + 1L]] }
  pos2 <- lapply(diag, function(g) as.list(as.integer(g)))
  shp2 <- lapply(diag, function(g) shape[[g[1L] + 1L]])
  list(c(pos1, pos2), c(shp1, shp2))
}

.ae_diag_push_down <- function(diag, indices, ndim) {
  positions <- .ae_positions_shape(as.list(seq_len(ndim) - 1L), diag)[[1L]]
  .ae_rec(function(x) if (!is.na(x) && x < length(positions)) positions[[x + 1L]] else NA_integer_, indices)
}

.ae_diag_push_up <- function(diag, indices, ndim) {
  positions <- .ae_positions_shape(as.list(seq_len(ndim) - 1L), diag)[[1L]]
  .ae_rec(function(x) {
    for (i in seq_along(positions)) {
      e <- positions[[i]]
      if (!is.list(e) && identical(as.integer(x), e)) return(i - 1L)
      if (is.list(e) && x %in% unlist(e)) return(i - 1L)
    }
    NA_integer_
  }, indices)
}

# flattened contraction indices push transformations, from utils.py
.ae_push_up_fn <- function(flat) {
  shifts_k <- 0L; shifts_v <- 0L
  i <- 1L; cumulative <- 0L
  while (i <= length(flat)) {
    j <- 1L
    while (i + j <= length(flat)) {
      if (flat[i] + j != flat[i + j]) break
      j <- j + 1L
    }
    cumulative <- cumulative + j
    w <- match(flat[i], shifts_k)
    if (is.na(w)) { shifts_k <- c(shifts_k, flat[i]); shifts_v <- c(shifts_v, cumulative) } else shifts_v[w] <- cumulative
    i <- i + j
  }
  o <- order(shifts_k)
  sk <- shifts_k[o]; sv <- shifts_v[o]
  function(j) {
    if (is.na(j) || j %in% flat) return(NA_integer_)
    as.integer(j - sv[findInterval(j, sk)])
  }
}

.ae_push_down_fn <- function(flat) {
  N <- flat[length(flat)] + 2L
  shifts <- setdiff(seq_len(N) - 1L, flat)
  function(j) {
    if (is.na(j)) return(NA_integer_)
    if (j < length(shifts)) shifts[j + 1L] else as.integer(j + shifts[length(shifts)] - length(shifts) + 1L)
  }
}

.ae_contr_push_up <- function(ci, indices) .ae_rec(.ae_push_up_fn(sort(unlist(ci))), indices)
.ae_contr_push_down <- function(ci, indices) .ae_rec(.ae_push_down_fn(sort(unlist(ci))), indices)

# ---------------------------------------------------------------------------
# ArrayDiagonal

.ae_diag_ind <- function(e) lapply(e$args[-1L], .ae_v)
.ae_contr_ind <- function(e) lapply(e$args[-1L], .ae_v)

#' @rdname array-expressions
#' @export
ArrayDiagonal <- function(expr, ..., canonicalize = FALSE, allow_trivial_diags = FALSE) {
  expr <- .ae_sympify(expr)
  diag <- lapply(list(...), function(g) sort(.ae_v(g)))
  shape <- .ae_shape(expr)
  for (g in diag) {
    if (any(g >= length(shape))) .perr("ValueError", "index is larger than expression shape")
    if (length(unique(vapply(shape[g + 1L], .k, ""))) != 1L) .perr("ValueError", "diagonalizing indices of different dimensions")
    if (!allow_trivial_diags && length(g) <= 1L) .perr("ValueError", "need at least two axes to diagonalize")
    if (anyDuplicated(g)) .perr("ValueError", "axis index cannot be repeated")
  }
  if (!length(diag)) return(expr)
  ps <- .ae_positions_shape(shape, diag)
  obj <- .ae_new("ArrayDiagonal", c(list(expr), lapply(diag, .ae_t)), ps[[2L]], .ae_sub_ndim_list(expr))
  if (canonicalize) .ae_diag_canon(obj) else obj
}

.ae_diag_canon <- function(self) {
  expr <- self$args[[1L]]
  diag <- .ae_diag_ind(self)
  if (any(lengths(diag) == 1L)) {
    triv <- which(lengths(diag) == 1L)
    triv_first <- vapply(diag[triv], `[`, 0L, 1L)
    short <- diag[lengths(diag) > 1L]
    short_keys <- vapply(short, paste, "", collapse = ",")
    full_keys <- vapply(diag, paste, "", collapse = ",")
    ndim1 <- get_ndim(self); ndim2 <- length(diag); ndim3 <- ndim1 - ndim2
    inv <- integer(); counter1 <- 0L
    down <- .ae_diag_push_down(short, as.list(seq_len(ndim1) - 1L), get_ndim(expr))
    for (i in down) {
      if (!is.list(i) && !is.na(i) && i %in% triv_first) inv <- c(inv, ndim3 + triv[match(i, triv_first)] - 1L)
      else if (!is.list(i)) { inv <- c(inv, counter1); counter1 <- counter1 + 1L }
      else inv <- c(inv, ndim3 + match(paste(unlist(i), collapse = ","), full_keys) - 1L)
    }
    p <- .ae_inv(inv)
    if (length(short)) return(.ae_permute_dims(do.call(.ae_diagonal, c(list(expr), short)), p))
    return(.ae_permute_dims(expr, p))
  }
  if (expr$type == "ArrayAdd") return(do.call(.ae_array_add, lapply(expr$args, function(a) do.call(.ae_diagonal, c(list(a), diag)))))
  if (expr$type == "ArrayDiagonal") return(.ae_diag_flatten(expr, get_ndim(self), diag))
  if (expr$type == "PermuteDims") return(.ae_diag_denest_pd(expr, diag))
  if (expr$type %in% c("ZeroArray", "ZeroMatrix")) return(do.call(ZeroArray, .ae_positions_shape(.ae_shape(expr), diag)[[2L]]))
  do.call(ArrayDiagonal, c(list(expr), diag))
}

.ae_diag_flatten <- function(expr, outer_dims, outer_diag) {
  inddown <- .ae_diag_push_down(outer_diag, as.list(seq_len(outer_dims) - 1L), get_ndim(expr))
  inddown2 <- .ae_diag_push_down(.ae_diag_ind(expr), inddown, get_ndim(expr$args[[1L]]))
  new_diag <- list()
  for (i in inddown2) {
    if (!is.list(i)) next
    new_diag[[length(new_diag) + 1L]] <- as.integer(unlist(i))
  }
  do.call(.ae_diagonal, c(list(expr$args[[1L]]), new_diag))
}

.ae_diag_denest_pd <- function(expr, diag) {
  p <- .ae_paf(expr$args[[2L]])
  back_diag <- lapply(diag, function(g) p[g + 1L])
  nondiag <- setdiff(seq_len(get_ndim(expr)) - 1L, unlist(diag))
  back_nondiag <- p[nondiag + 1L]
  remap <- match(back_nondiag, sort(back_nondiag)) - 1L
  np <- c(remap, length(remap) + seq_along(back_diag) - 1L)
  .ae_permute_dims(do.call(.ae_diagonal, c(list(expr$args[[1L]]), back_diag)), np)
}

# ---------------------------------------------------------------------------
# ArrayElementwiseApplyFunc

.ae_as_lambda <- function(f) {
  if (is_sym(f) && f$type == "Lambda") return(f)
  d <- Dummy("d")
  if (is.character(f)) return(Lambda(d, fcall(f, list(d))))
  Lambda(d, f(d))
}
.ae_lambda_call <- function(lam, x) apply_lambda(lam, x)

#' @rdname array-expressions
#' @export
ArrayElementwiseApplyFunc <- function(fun, element) {
  fun <- .ae_as_lambda(fun)
  element <- .ae_sympify(element)
  .ae_new("ArrayElementwiseApplyFunc", list(fun, element), .ae_shape(element), .ae_sub_ndim_list(element))
}
register_rebuilder("ArrayElementwiseApplyFunc", function(e, args) ArrayElementwiseApplyFunc(args[[1L]], args[[2L]]))

# _get_function_fdiff: the derivative of the applied function, as a Lambda
.ae_function_fdiff <- function(lam) {
  d <- Dummy("d")
  fd <- .diff1(.ae_lambda_call(lam, d), d)
  if (fd$type == "Function" && length(fd$args) == 1L && identical(fd$args[[1L]]$key, d$key)) {
    d2 <- Dummy("d")
    return(Lambda(d2, fcall(fd$name, list(d2))))
  }
  Lambda(d, fd)
}

# ---------------------------------------------------------------------------
# ArrayContraction

.ae_sort_contr <- function(ci) {
  ci <- lapply(ci, function(g) sort(as.integer(g)))
  if (length(ci) > 1L) ci <- ci[order(vapply(ci, min, 0L), method = "radix")]
  ci
}

#' @rdname array-expressions
#' @export
ArrayContraction <- function(expr, ..., canonicalize = FALSE) {
  ci <- .ae_sort_contr(lapply(list(...), .ae_v))
  expr <- .ae_sympify(expr)
  shape <- .ae_shape(expr)
  for (g in ci) if (length(unique(vapply(shape[g + 1L], .k, ""))) != 1L) .perr("ValueError", "contracting indices of different dimensions")
  if (length(shape)) shape <- shape[!(seq_along(shape) - 1L) %in% unlist(ci)]
  obj <- .ae_new("ArrayContraction", c(list(expr), lapply(ci, .ae_t)), shape, .ae_sub_ndim_list(expr))
  if (canonicalize) .ae_contr_canon(obj) else obj
}

.ae_contr_canon <- function(self) {
  expr <- self$args[[1L]]
  ci <- .ae_contr_ind(self)
  if (!length(ci)) return(expr)
  if (expr$type == "ArraySum")
    return(do.call(ArraySum, c(list(do.call(.ae_contraction, c(list(expr$args[[1L]]), ci))), expr$args[-1L])))
  if (expr$type == "ArrayContraction") return(.ae_contr_flatten(expr, ci))
  if (expr$type %in% c("ZeroArray", "ZeroMatrix")) {
    sh <- .ae_shape(expr)
    return(do.call(ZeroArray, sh[!(seq_along(sh) - 1L) %in% unlist(ci)]))
  }
  if (expr$type == "PermuteDims") return(.ae_contr_denest_pd(expr, ci))
  if (expr$type == "ArrayTensorProduct") {
    r <- .ae_sort_fully_contracted_args(expr, ci); expr <- r[[1L]]; ci <- r[[2L]]
    r <- .ae_lower_contraction_to_addends(expr, ci); expr <- r[[1L]]; ci <- r[[2L]]
    if (!length(ci)) return(expr)
  }
  if (expr$type == "ArrayDiagonal") return(.ae_contr_denest_diag(expr, ci))
  if (expr$type == "ArrayAdd") return(do.call(.ae_array_add, lapply(expr$args, function(a) do.call(.ae_contraction, c(list(a), ci)))))
  shape <- .ae_shape(expr)
  triv <- Filter(function(g) length(g) == 1L && is_one_num(shape[[g + 1L]]), ci)
  if (length(triv)) {
    new_shape <- shape[!(seq_along(shape) - 1L) %in% unlist(ci)]
    if (length(new_shape)) {
      tk <- vapply(triv, paste, "", collapse = ",")
      remaining <- Filter(function(g) !(paste(g, collapse = ",") %in% tk), ci)
      newexpr <- do.call(.ae_contraction, c(list(expr), remaining))
      if (.ae_shape_eq(.ae_shape(newexpr), new_shape)) return(newexpr)
      return(Reshape(newexpr, new_shape))
    }
  }
  do.call(ArrayContraction, c(list(expr), ci))
}

.ae_lower_contraction_to_addends <- function(expr, ci) {
  if (expr$type != "ArrayTensorProduct") return(list(expr, ci))
  cum <- cumsum(c(0L, .ae_sub_ndim_list(expr)))
  remaining <- list()
  per_arg <- rep(list(list()), length(expr$args))
  backshift <- integer()
  for (g in ci) {
    placed <- FALSE
    for (j in seq_along(expr$args)) {
      if (expr$args[[j]]$type != "ArrayAdd") next
      if (all(cum[j] <= g & g < cum[j + 1L])) {
        per_arg[[j]][[length(per_arg[[j]]) + 1L]] <- g - cum[j]
        backshift <- union(backshift, g)
        placed <- TRUE
        break
      }
    }
    if (!placed) remaining[[length(remaining) + 1L]] <- g
  }
  if (length(remaining) == length(ci)) return(list(expr, ci))
  total <- get_ndim(expr)
  shifts <- cumsum(as.integer((seq_len(total) - 1L) %in% backshift))
  remaining <- lapply(remaining, function(g) g - shifts[g + 1L])
  ret <- do.call(.ae_tensor_product, lapply(seq_along(expr$args), function(j) do.call(.ae_contraction, c(list(expr$args[[j]]), per_arg[[j]]))))
  list(ret, remaining)
}

.ae_index_shifts <- function(expr) {
  all_inner <- sort(unlist(.ae_contr_ind(expr)))
  total <- .ae_sub_ndim(expr)
  inner <- length(all_inner)
  outer <- total - inner
  shifts <- integer(outer)
  counter <- 0L; pointer <- 0L
  for (i in seq_len(outer)) {
    while (pointer < inner && counter >= all_inner[pointer + 1L]) { counter <- counter + 1L; pointer <- pointer + 1L }
    shifts[i] <- shifts[i] + pointer
    counter <- counter + 1L
  }
  shifts
}

.ae_outer_to_inner <- function(expr, outer) {
  shifts <- .ae_index_shifts(expr)
  lapply(outer, function(g) as.integer(shifts[g + 1L] + g))
}

.ae_contr_flatten <- function(expr, outer) {
  inner <- .ae_contr_ind(expr)
  do.call(.ae_contraction, c(list(expr$args[[1L]]), inner, .ae_outer_to_inner(expr, outer)))
}

.ae_contr_denest_pd <- function(expr, ci) {
  p <- .ae_paf(expr$args[[2L]])
  new_ci <- lapply(ci, function(g) p[g + 1L])
  new_plist <- p[!p %in% unlist(new_ci)]
  new_plist <- unlist(.ae_contr_push_up(new_ci, as.list(new_plist)))
  .ae_permute_dims(do.call(.ae_contraction, c(list(expr$args[[1L]]), new_ci)), as.integer(new_plist))
}

.ae_contr_denest_diag <- function(expr, ci) {
  diag <- .ae_diag_ind(expr)
  down <- .ae_diag_push_down(diag, lapply(ci, as.list), get_ndim(expr$args[[1L]]))
  down <- lapply(down, function(g) as.integer(unlist(g)))
  new_ci <- list()
  for (g in down) {
    ind <- g
    for (j in seq_along(diag)) {
      if (is.null(diag[[j]])) next
      if (any(g %in% diag[[j]])) { ind <- c(ind, diag[[j]]); diag[j] <- list(NULL) }
    }
    new_ci[[length(new_ci) + 1L]] <- sort(unique(ind))
  }
  down_diag <- Filter(Negate(is.null), diag)
  new_diag <- lapply(.ae_contr_push_up(new_ci, lapply(down_diag, as.list)), function(g) as.integer(unlist(g)))
  do.call(.ae_diagonal, c(list(do.call(.ae_contraction, c(list(expr$args[[1L]]), new_ci))), new_diag))
}

.ae_default_sort_order <- function(items) {
  keys <- lapply(items, .sort_key)
  idx <- seq_along(items)
  if (length(idx) < 2L) return(idx)
  for (i in 2:length(idx)) {
    j <- i
    while (j > 1L && .key_cmp(keys[[idx[j - 1L]]], keys[[idx[j]]]) > 0) {
      tmp <- idx[j]; idx[j] <- idx[j - 1L]; idx[j - 1L] <- tmp; j <- j - 1L
    }
  }
  idx
}

.ae_sort_fully_contracted_args <- function(expr, ci) {
  cum <- cumsum(c(0L, .ae_sub_ndim_list(expr)))
  n <- length(expr$args)
  blocks <- lapply(seq_len(n), function(i) seq_len(cum[i + 1L] - cum[i]) + cum[i] - 1L)
  flat <- unlist(ci)
  full <- vapply(blocks, function(b) all(b %in% flat), NA)
  fc <- which(full)
  fc <- fc[.ae_default_sort_order(expr$args[fc])]
  new_pos <- c(fc, which(!full))
  ipa <- .ae_inv(unlist(blocks[new_pos]))
  new_ci <- .ae_sort_contr(lapply(ci, function(g) ipa[g + 1L]))
  list(do.call(.ae_tensor_product, expr$args[new_pos]), new_ci)
}

# ---------------------------------------------------------------------------
# Reshape

#' @rdname array-expressions
#' @export
Reshape <- function(expr, shape) {
  expr <- .ae_sympify(expr)
  shape <- if (is_sym(shape) && shape$type == "Tuple") shape$args else lapply(as.list(shape), as_sym)
  d <- .sub(.mul(.ae_shape(expr)), .mul(shape))
  if (isTRUE(.facts(d)[["nonzero"]]) || (is_exact_number(d) && !is_zero_num(d))) .perr("ValueError", "shape mismatch")
  .ae_new("Reshape", list(expr, new_sym("Tuple", shape)), shape)
}
register_rebuilder("Reshape", function(e, args) Reshape(args[[1L]], args[[2L]]))

# ---------------------------------------------------------------------------
# ArraySum: a Sum whose summand is an array expression

#' @rdname array-expressions
#' @export
ArraySum <- function(f, ...) {
  f <- .ae_sympify(f)
  lims <- lapply(list(...), function(l) if (is_sym(l) && l$type == "Tuple") l else do.call(Tuple, as.list(l)))
  .ae_new("ArraySum", c(list(f), lims), .ae_shape(f), name = "Sum")
}
register_str_printer("ArraySum", function(e) .str(new_sym("Sum", e$args)))
register_bound_vars("ArraySum", function(e) vapply(e$args[-1L], function(l) l$args[[1L]]$key, ""))
register_rebuilder("ArraySum", function(e, args) do.call(ArraySum, args))

# ---------------------------------------------------------------------------
# canonicalizing constructors (_array_tensor_product and friends)

.ae_tensor_product <- function(...) {
  args <- lapply(list(...), .ae_sympify)
  if (all(vapply(args, .ae_is_plain_scalar, NA))) return(.mul(args))
  do.call(ArrayTensorProduct, c(args, list(canonicalize = TRUE)))
}
.ae_contraction <- function(expr, ...) ArrayContraction(expr, ..., canonicalize = TRUE)
.ae_diagonal <- function(expr, ..., allow_trivial_diags = FALSE)
  ArrayDiagonal(expr, ..., canonicalize = TRUE, allow_trivial_diags = allow_trivial_diags)
.ae_permute_dims <- function(expr, permutation) PermuteDims(expr, permutation, canonicalize = TRUE)
.ae_array_add <- function(...) {
  args <- lapply(list(...), .ae_sympify)
  if (any(vapply(args, .ae_is_plain_scalar, NA)) && all(vapply(args, get_ndim, 0L) == 0L)) return(.add(args))
  do.call(ArrayAdd, c(args, list(canonicalize = TRUE)))
}

# ---------------------------------------------------------------------------
# rebuilders, doit and printing

register_rebuilder("ArrayTensorProduct", function(e, args) do.call(ArrayTensorProduct, args))
register_rebuilder("ArrayAdd", function(e, args) do.call(ArrayAdd, args))
register_rebuilder("PermuteDims", function(e, args) PermuteDims(args[[1L]], args[[2L]]))
register_rebuilder("ArrayDiagonal", function(e, args) do.call(ArrayDiagonal, c(list(args[[1L]]), args[-1L], list(allow_trivial_diags = TRUE))))
register_rebuilder("ArrayContraction", function(e, args) do.call(ArrayContraction, c(list(args[[1L]]), args[-1L])))

.ae_canonicalize <- function(e) switch(e$type,
  ArrayTensorProduct = .ae_tp_canon(e), ArrayAdd = .ae_add_canon(e), PermuteDims = .ae_pd_canon(e),
  ArrayDiagonal = .ae_diag_canon(e), ArrayContraction = .ae_contr_canon(e), e)

.ae_doit_deep <- function(e, deep) if (deep) lapply(e$args, function(a) doit(a)) else e$args
for (.t in c("ArrayTensorProduct", "PermuteDims", "ArrayContraction"))
  register_doit(.t, function(e, deep = TRUE) .ae_canonicalize(sym_rebuild(e, .ae_doit_deep(e, deep))))
register_doit("ArrayAdd", function(e, deep = TRUE) {
  if (deep) return(do.call(.ae_array_add, lapply(e$args, doit)))
  .ae_add_canon(e)
})
register_doit("ArrayDiagonal", function(e, deep = TRUE) {
  a <- .ae_doit_deep(e, deep)
  .ae_diag_canon(do.call(ArrayDiagonal, c(list(a[[1L]]), a[-1L], list(allow_trivial_diags = TRUE))))
})
register_doit("Reshape", function(e, deep = TRUE) {
  x <- if (deep) doit(e$args[[1L]]) else e$args[[1L]]
  if (is_sym(x) && x$type == "ImmutableDenseNDimArray") x <- .ae_ndarray_r(x)
  if (inherits(x, "symmatrix") || (is_sym(x) && x$type == "ImmutableMatrix")) x <- Array(as_explicit(x)$data, c(as_explicit(x)$nrow, as_explicit(x)$ncol))
  if (inherits(x, "SymArray")) return(Array(x$data, vapply(e$args[[2L]]$args, .int_val, 0)))
  Reshape(x, e$args[[2L]])
})
register_doit("ArraySum", function(e, deep = TRUE) .ae_arraysum_doit(e, deep))
rm(.t)

.ae_arraysum_doit <- function(e, deep) {
  f <- if (deep) doit(e$args[[1L]]) else e$args[[1L]]
  lims <- e$args[-1L]
  if (all(vapply(lims, function(l) length(l$args) == 3L && l$args[[2L]]$type == "Integer" && l$args[[3L]]$type == "Integer", NA))) {
    lo <- vapply(lims, function(l) .int_val(l$args[[2L]]), 0)
    hi <- vapply(lims, function(l) .int_val(l$args[[3L]]), 0)
    if (any(hi < lo)) return(do.call(ZeroArray, .ae_shape(e)))
    grid <- as.matrix(expand.grid(rev(lapply(seq_along(lo), function(i) lo[i]:hi[i]))))[, rev(seq_along(lo)), drop = FALSE]
    terms <- lapply(seq_len(nrow(grid)), function(r)
      subs(f, lapply(seq_along(lims), function(i) list(lims[[i]]$args[[1L]], .int(grid[r, i])))))
    if (all(vapply(terms, .is_matexpr, NA))) return(Reduce(`+`, terms))
    return(do.call(.ae_array_add, terms))
  }
  done <- doit(new_sym("Sum", c(list(f), lims)))
  if (is_sym(done) && is_zero_num(done)) return(do.call(ZeroArray, .ae_shape(e)))
  if (is_sym(done) && done$type == "Sum") return(do.call(ArraySum, done$args))
  done
}

.ae_print_basic <- function(name) { force(name); function(e) paste0(name, "(", paste(vapply(e$args, .str, ""), collapse = ", "), ")") }
for (.t in c("ZeroArray", "OneArray", "ArrayTensorProduct", "ArrayAdd", "PermuteDims", "ArrayDiagonal",
             "ArrayContraction", "ArrayElementwiseApplyFunc"))
  register_str_printer(.t, .ae_print_basic(.t))
register_str_printer("Reshape", function(e) paste0("Reshape(", .str(e$args[[1L]]), ", ", .str(e$args[[2L]]), ")"))
rm(.t)

# a*b for matrix expressions: MatMul(a, b).doit()
.ae_mmul <- function(...) doit(do.call(MatMul, list(...)))

#' Collect sums of tensor products
#'
#' Rewrites every [ArrayAdd()] in `expr` whose terms are tensor products that
#' differ in at most one factor, so that
#' `c1*(A x B) + c2*(C x B)` becomes `(c1*A + c2*C) x B`.
#'
#' @param expr An array expression.
#' @return The rewritten expression.
#' @examples
#' k <- Symbol("k")
#' A <- MatrixSymbol("A", k, k); B <- MatrixSymbol("B", k, k); C <- MatrixSymbol("C", k, k)
#' collect_tensor_products(ArrayAdd(ArrayTensorProduct(A, B), ArrayTensorProduct(A, C)))
#' @export
collect_tensor_products <- function(expr) {
  if (!is_sym(expr)) return(expr)
  if (length(expr$args)) {
    newargs <- lapply(expr$args, collect_tensor_products)
    if (!all(mapply(function(a, b) identical(a$key, b$key), newargs, expr$args))) expr <- sym_rebuild(expr, newargs)
  }
  if (is_sym(expr) && expr$type == "ArrayAdd") return(.ae_collect_array_add(expr))
  expr
}

.ae_scale_slot <- function(coeff, arg) if (is_one_num(coeff)) arg else .ae_term_from_coeff_arrays(coeff, list(arg))

.ae_merge_single_slot <- function(t1, t2) {
  a1 <- t1[[2L]]; a2 <- t2[[2L]]
  if (length(a1) != length(a2) || !length(a1)) return(NULL)
  if (any(mapply(function(a, b) .ae_skey(get_shape(a)) != .ae_skey(get_shape(b)), a1, a2))) return(NULL)
  d <- which(!mapply(.ae_eq, a1, a2))
  if (!length(d)) return(list(.add(list(t1[[1L]], t2[[1L]])), a1))
  if (length(d) != 1L) return(NULL)
  sa <- .ae_scale_slot(t1[[1L]], a1[[d]]); sb <- .ae_scale_slot(t2[[1L]], a2[[d]])
  a1[[d]] <- if (.is_matexpr(sa) && .is_matexpr(sb)) sa + sb else .ae_array_add(sa, sb)
  list(S_One, a1)
}

.ae_collect_array_add <- function(add) {
  terms <- lapply(add$args, .ae_term_coeff_arrays)
  repeat {
    merged <- list()
    for (term in terms) {
      done <- FALSE
      for (i in seq_along(merged)) {
        comb <- .ae_merge_single_slot(merged[[i]], term)
        if (!is.null(comb)) { merged[[i]] <- comb; done <- TRUE; break }
      }
      if (!done) merged[[length(merged) + 1L]] <- term
    }
    if (length(merged) == length(terms)) break
    terms <- merged
  }
  if (length(terms) == length(add$args)) return(add)
  do.call(.ae_array_add, lapply(terms, function(t) .ae_term_from_coeff_arrays(t[[1L]], t[[2L]])))
}

# NDimArray.__str__: nested lists of the entries, or the entry for rank 0
.ae_str_ndarray <- function(entries, shape) {
  s <- vapply(entries, .str, "")
  if (!length(shape)) return(s[[1L]])
  rec <- function(d, off) {
    step <- prod(shape[-seq_len(d)])
    parts <- if (d == length(shape)) s[off + seq_len(shape[d])]
             else vapply(seq_len(shape[d]) - 1L, function(k) rec(d + 1L, off + k * step), "")
    paste0("[", paste(parts, collapse = ", "), "]")
  }
  rec(1L, 0L)
}
register_str_printer("ImmutableDenseNDimArray", function(e) .ae_str_ndarray(e$args, e$extra$shape))
