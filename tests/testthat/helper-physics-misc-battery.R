# Helpers for the test-physics-*-sympy.R files: the Python semantics the
# battery cases are written in, loaded into an environment whose parent is the
# package namespace.

.physmisc_battery_env <- local({
  env <- NULL
  function() {
    if (is.null(env)) {
      env <<- new.env(parent = asNamespace("symcompR"))
      eval(parse(text = r"---[# Python semantics for the transpiled battery (R side). Sourced into the
# global environment after the package is attached.

.ns <- asNamespace("symcompR")
.I <- function(name) get(name, envir = .ns)

# ---------------------------------------------------------------- values
# a Python float literal: a whole one becomes a Float (R would make it an
# Integer in sym arithmetic), any other stays an R double
.py_flt <- function(x) if (is.finite(x) && x == round(x)) Float(x) else x
.py_abs <- function(x) if (is_sym(x) || inherits(x, "symmatrix")) Abs(x) else abs(x)
.py_np_float64 <- function(x) .py_float(x)
.py_np_float32 <- function(x) .py_flt(as.numeric(x))
.py_np_int64 <- function(x) as.integer(x)
.py_np_zeros <- function(n) rep(0, n)
.py_np_array <- function(x) x
.py_np_sum <- function(l) Reduce(.py_add, l)
.py_np_isnan <- function(v) is.nan(unlist(v))
.py_bigint <- function(s) Integer(s)
.py_S <- function(x, ...) {
  if (is.character(x)) return(sympify(x))
  if (inherits(x, "symmatrix")) return(x)
  if (is.list(x) && !is_sym(x)) return(lapply(x, .py_S))
  as_sym(x)
}
.py_symbols <- function(names, ...) {
  r <- .I("symbols")(names, ...)
  if (is_sym(r)) return(r)
  .I(".pytuple")(unname(r))
}
.py_unit <- function(name) {
  u <- physics_units()
  v <- u[[name]]
  if (is.null(v)) stop("no unit ", name)
  v
}
.py_hbar <- tryCatch(.I(".physics_hbar")(), error = function(e) Symbol("hbar"))
.S1 <- function(f) function(x, ...) f(.py_S(x), ...)
.py_sqrt <- .S1(sqrt)
.py_exp <- .S1(exp)
.py_log <- function(x, b = NULL) if (is.null(b)) log(.py_S(x)) else log(.py_S(x), .py_S(b))
.py_sin <- .S1(sin); .py_cos <- .S1(cos); .py_tan <- .S1(tan)
.py_asin <- .S1(asin); .py_acos <- .S1(acos); .py_atan <- .S1(atan)
.py_sinh <- .S1(sinh); .py_cosh <- .S1(cosh); .py_tanh <- .S1(tanh)
.py_gamma <- .S1(gamma)
.py_floor <- .S1(floor); .py_ceiling <- .S1(ceiling)
.py_sign <- .S1(sign)
.py_arg <- function(x) arg(.py_S(x))
.py_deg <- function() NULL
.py_numer <- function(x) numer(.py_S(x))
.py_denom <- function(x) denom(.py_S(x))
.py_N <- function(x, n = 15, ...) {
  if (.is_obj(x)) return(.py_method(x, "evalf", n, ...))
  if (inherits(x, "symmatrix")) return(.I(".mat_map")(x, function(z) N(z, n)))
  N(.py_S(x), n)
}
.py_diff <- function(f, ...) {
  a <- list(...)
  do.call(diff, c(list(.py_S(f)), a))
}
.py_Matrix <- function(...) {
  a <- list(...)
  if (length(a) == 0L) return(Matrix(list()))
  if (length(a) == 1L) {
    d <- a[[1L]]
    if (inherits(d, "symmatrix")) return(.I(".new_mat")(d$data, d$nrow, d$ncol))
    if (is_sym(d) && d$type %in% c("ImmutableMatrix", "MatrixBlock")) return(Matrix(d))
    d <- .py_iter(d)
    if (!length(d)) return(Matrix(list()))
    if (all(vapply(d, function(r) inherits(r, "symmatrix"), TRUE))) {
      out <- d[[1L]]
      if (length(d) > 1L) for (k in 2:length(d)) out <- rbind(out, d[[k]])
      return(out)
    }
    if (all(vapply(d, function(r) .is_seq(r), TRUE))) {
      rows <- lapply(d, function(r) lapply(.py_iter(r), as_sym))
      return(Matrix(rows))
    }
    return(Matrix(lapply(d, as_sym)))
  }
  nr <- as.integer(a[[1L]]); nc <- as.integer(a[[2L]])
  f <- a[[3L]]
  if (is.function(f)) return(Matrix(nr, nc, function(i, j) f(i - 1L, j - 1L)))
  Matrix(nr, nc, f = lapply(.py_iter(f), as_sym))
}
.py_Mul <- function(..., evaluate = TRUE) {
  a <- lapply(list(...), .py_S)
  if (isTRUE(evaluate)) return(do.call(Mul, a))
  if (!length(a)) return(S_One)
  if (length(a) == 1L) return(a[[1L]])
  .I("new_sym")("Mul", a)
}
.py_Add <- function(..., evaluate = TRUE) {
  a <- lapply(list(...), .py_S)
  if (isTRUE(evaluate)) return(do.call(Add, a))
  if (!length(a)) return(S_Zero)
  if (length(a) == 1L) return(a[[1L]])
  .I("new_sym")("Add", a)
}
.py_Pow <- function(a, b, evaluate = TRUE) {
  if (isTRUE(evaluate)) return(Pow(a, b))
  .I("new_sym")("Pow", list(.py_S(a), .py_S(b)))
}
.py_Tuple <- function(...) do.call(Tuple, lapply(list(...), .py_S))
.py_all_close <- function(a, b, rtol = 1e-5, atol = 1e-8) {
  fa <- .flatnum(a); fb <- .flatnum(b)
  length(fa) == length(fb) && all(base::Mod(fa - fb) <= atol + rtol * base::Mod(fb))
}
.flatnum <- function(x) {
  if (is.list(x) && !is_sym(x)) return(as.complex(unlist(lapply(x, .flatnum))))
  if (is_sym(x)) return(as.complex(x))
  as.complex(x)
}
.py_isclose <- function(a, b, rel_tol = 1e-9, abs_tol = 0) {
  a <- .py_float(a); b <- .py_float(b)
  abs(a - b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)
}

# ---------------------------------------------------------------- objects
.is_obj <- function(x) inherits(x, c("pyobj", "pyenv"))
.is_seq <- function(x) (is.list(x) && !is.object(x)) || inherits(x, c("pytuple", "pylist"))
.tup <- function(x) .I(".pytuple")(x)
.py_tup <- function(...) .tup(list(...))
.py_setlit <- function(x) {
  x <- .py_iter(x)
  keys <- vapply(x, .hkey, "")
  structure(x[!duplicated(keys)], class = "pyset")
}
.py_dict <- function(keys, values) structure(list(keys = keys, values = values), class = "pydict")
.py_dict_from_pairs <- function(pairs) .py_dict(lapply(pairs, `[[`, 1L), lapply(pairs, `[[`, 2L))
.py_dict_fn <- function(x = NULL, ...) {
  if (is.null(x)) return(.py_dict(list(), list()))
  if (inherits(x, "pydict")) return(x)
  .py_dict_from_pairs(lapply(.py_iter(x), .py_iter))
}
.hkey <- function(x) {
  if (is.null(x)) return("None")
  if (is_sym(x)) return(x$key)
  if (is.numeric(x) && length(x) == 1L && x == round(x)) return(paste0("Z", format(x, scientific = FALSE)))
  if (.is_seq(x)) return(paste0("(", paste(vapply(x, .hkey, ""), collapse = ","), ")"))
  if (inherits(x, "symmatrix")) return(paste0("M", paste(vapply(x$data, .hkey, ""), collapse = ",")))
  if (is.character(x)) return(paste0("s", x))
  paste(format(x), collapse = " ")
}
.py_comp <- function(it, f) unlist(lapply(.py_iter(it), f), recursive = FALSE)
.py_splat <- function(x) structure(list(x), class = "pysplat")
.py_call <- function(fn, ...) {
  a <- list(...)
  out <- list()
  for (i in seq_along(a)) {
    nm <- names(a)[i]
    v <- a[[i]]
    if (inherits(v, "pysplat")) {
      v <- v[[1L]]
      if (inherits(v, "pydict")) {
        vals <- v$values
        names(vals) <- vapply(v$keys, function(k) if (is.character(k)) k else sstr(k), "")
        out <- c(out, vals)
      } else out <- c(out, .py_iter(v))
    } else {
      l <- list(v)
      if (!is.null(nm) && nzchar(nm)) names(l) <- nm
      out <- c(out, l)
    }
  }
  do.call(fn, out)
}
.py_iter <- function(x) {
  if (is.null(x)) return(list())
  if (inherits(x, "pydict")) return(x$keys)
  if (inherits(x, c("pytuple", "pyset", "pylist"))) return(unclass(x))
  if (inherits(x, "symmatrix")) return(x$data)
  if (is_sym(x)) {
    if (x$type %in% c("Tuple", "FiniteSet")) return(x$args)
    if (.is_obj(x)) return(as.list(.py_method(x, "__iter__")))
    stop("TypeError: object is not iterable")
  }
  if (is.environment(x) && .is_obj(x)) return(as.list(.py_method(x, "__iter__")))
  if (is.list(x)) return(unclass(x))
  if (is.character(x) && length(x) == 1L) return(as.list(strsplit(x, "")[[1L]]))
  as.list(x)
}
.py_cap <- function(x, n) { x <- .py_iter(x); x[seq_len(min(n, length(x)))] }
.py_range <- function(a, b = NULL, step = 1L) {
  if (is.null(b)) { b <- a; a <- 0L }
  a <- as.integer(.py_int(a)); b <- as.integer(.py_int(b)); step <- as.integer(step)
  if ((step > 0 && a >= b) || (step < 0 && a <= b)) return(integer())
  v <- seq(a, b, by = step)
  v[if (step > 0) v < b else v > b]
}
.py_list <- function(x = NULL) if (is.null(x)) list() else .py_iter(x)
.py_tuple <- function(x = NULL) .tup(.py_list(x))
.py_set <- function(x = NULL) .py_setlit(.py_list(x))
.py_len <- function(x) {
  if (inherits(x, "symmatrix")) return(x$nrow * x$ncol)
  if (inherits(x, "pydict")) return(length(x$keys))
  if (.is_obj(x)) return(.py_method(x, "__len__"))
  if (is_sym(x)) return(length(x$args))
  if (is.character(x) && length(x) == 1L) return(nchar(x))
  length(.py_iter(x))
}
.py_sum <- function(x, start = 0L) Reduce(.py_add, .py_iter(x), start)
.py_all <- function(x) all(vapply(.py_iter(x), .py_truthy, TRUE))
.py_any <- function(x) any(vapply(.py_iter(x), .py_truthy, TRUE))
.py_float <- function(x) {
  if (is_sym(x)) {
    if (x$type %in% c("Float", "Integer", "Rational")) return(as.numeric(x))
    v <- as.complex(x)
    if (Im(v) != 0) stop("TypeError: Cannot convert complex to float")
    return(Re(v))
  }
  as.numeric(x)
}
.py_int <- function(x) { v <- .py_float(x); trunc(v) }
.py_complex <- function(x, y = 0) if (is_sym(x)) as.complex(x) else complex(real = x, imaginary = y)
.py_str <- function(x) {
  if (is.character(x)) return(x)
  .fmt(x, py_str = TRUE)
}
.py_bool <- function(x) .py_truthy(x)
.py_round <- function(x, n = 0) round(.py_float(x), n)
.py_min <- function(...) { a <- list(...); if (length(a) == 1L) a <- .py_iter(a[[1L]]); Reduce(function(p, q) if (.py_truthy(.py_lt(q, p))) q else p, a) }
.py_max <- function(...) { a <- list(...); if (length(a) == 1L) a <- .py_iter(a[[1L]]); Reduce(function(p, q) if (.py_truthy(.py_gt(q, p))) q else p, a) }
.py_zip <- function(...) { a <- lapply(list(...), .py_iter); n <- min(vapply(a, length, 0L)); lapply(seq_len(n), function(i) .tup(lapply(a, `[[`, i))) }
.py_enumerate <- function(x, start = 0L) { x <- .py_iter(x); lapply(seq_along(x), function(i) .tup(list(i - 1L + start, x[[i]]))) }
.py_reversed <- function(x) rev(.py_iter(x))
.py_sorted <- function(x, key = NULL, reverse = FALSE) {
  x <- .py_iter(x)
  k <- if (is.null(key)) x else lapply(x, key)
  o <- seq_along(x)
  if (length(x) > 1L) for (i in 2:length(x)) { j <- i; while (j > 1L && .py_truthy(.py_lt(k[[o[j]]], k[[o[j - 1L]]]))) { t <- o[j]; o[j] <- o[j - 1L]; o[j - 1L] <- t; j <- j - 1L } }
  r <- x[o]
  if (reverse) rev(r) else r
}
.py_map <- function(f, ...) { a <- lapply(list(...), .py_iter); lapply(seq_along(a[[1L]]), function(i) do.call(f, lapply(a, `[[`, i))) }
.py_type <- function(x) if (is_sym(x)) x$type else class(x)[1L]
.py_hasattr <- function(x, nm) !inherits(tryCatch(.py_attr(x, nm), error = function(e) e), "error")
.py_isinstance <- function(x, cls) {
  if (is.function(cls)) {
    nm <- attr(cls, "pyclass") %||% attr(cls, "pyname")
    if (is.null(nm)) for (n in c("Rational", "Integer", "Float", "Symbol", "Matrix", "ImmutableMatrix", "Mul", "Add", "Pow"))
      if (identical(cls, get0(n, envir = .ns))) { nm <- n; break }
    cls <- nm
  }
  if (inherits(cls, "pytuple") || is.list(cls)) return(any(vapply(cls, function(c) .py_isinstance(x, c), TRUE)))
  if (is.null(cls)) return(FALSE)
  switch(cls,
    int = is.numeric(x) && length(x) == 1L && x == round(x),
    float = is.double(x),
    tuple = inherits(x, "pytuple"), list = is.list(x) && !is.object(x),
    Float = is_sym(x) && x$type == "Float",
    Rational = is_sym(x) && x$type %in% c("Rational", "Integer"),
    Integer = is_sym(x) && x$type == "Integer",
    Number = is_sym(x) && x$type %in% c("Integer", "Rational", "Float"),
    Symbol = is_sym(x) && x$type %in% c("Symbol", "Dummy"),
    Matrix = , MatrixBase = , MutableDenseMatrix = inherits(x, "symmatrix"),
    ImmutableMatrix = , ImmutableDenseMatrix = inherits(x, "symmatrix") || (is_sym(x) && x$type == "ImmutableMatrix"),
    Mul = is_sym(x) && x$type == "Mul", Add = is_sym(x) && x$type == "Add", Pow = is_sym(x) && x$type == "Pow",
    Basic = , Expr = is_sym(x),
    LinearTimeInvariant = , SISOLinearTimeInvariant = , MIMOLinearTimeInvariant = , StateSpaceBase = ,
    TransferFunctionBase = inherits(x, cls),
    if (.is_obj(x)) .obj_isinstance(x, cls) else stop("isinstance: unknown class ", cls))
}
.obj_isinstance <- function(x, cls) {
  f <- tryCatch(.I(".pyobj_isinstance"), error = function(e) NULL)
  if (!is.null(f)) return(f(x, cls))
  inherits(x, cls)
}

# ---------------------------------------------------------------- truthiness, comparisons
.py_truthy <- function(x) {
  if (is.null(x)) return(FALSE)
  if (is.logical(x)) return(length(x) == 1L && !is.na(x) && x)
  if (is.numeric(x) || is.complex(x)) return(length(x) == 1L && x != 0)
  if (is.character(x)) return(nzchar(x))
  if (is_sym(x)) {
    if (x$type == "BooleanAtom") return(x$name == "true")
    if (is_number_type_(x)) return(!(x$type %in% c("Integer", "Rational") && x$key %in% c("Z0")) && !(x$type == "Float" && as.numeric(x) == 0))
    if (x$type %in% .I(".REL_TYPES") || x$type %in% c("And", "Or", "Not")) stop("TypeError: cannot determine truth value of Relational")
    return(TRUE)
  }
  if (inherits(x, "pydict")) return(length(x$keys) > 0L)
  if (inherits(x, "symmatrix")) return(TRUE)
  if (is.list(x)) return(length(x) > 0L)
  TRUE
}
is_number_type_ <- function(x) x$type %in% c("Integer", "Rational", "Float")
.py_eq <- function(a, b) {
  if (is.null(a) || is.null(b)) return(is.null(a) && is.null(b))
  if (inherits(a, "symmatrix") || inherits(b, "symmatrix")) {
    if (!(inherits(a, "symmatrix") && inherits(b, "symmatrix"))) return(FALSE)
    return(a$nrow == b$nrow && a$ncol == b$ncol && all(vapply(seq_along(a$data), function(k) identical(a$data[[k]]$key, b$data[[k]]$key), TRUE)))
  }
  if (.is_seq(a) || .is_seq(b)) {
    if (!(.is_seq(a) && .is_seq(b)) || length(a) != length(b)) return(FALSE)
    return(all(vapply(seq_along(a), function(i) .py_eq(a[[i]], b[[i]]), TRUE)))
  }
  if (is.character(a) || is.character(b)) return(identical(a, b))
  if (is.logical(a) && is.logical(b)) return(identical(a, b))
  if (is_sym(a) || is_sym(b)) {
    if (.is_obj(a) || .is_obj(b)) return(is_sym(a) && is_sym(b) && identical(a$key, b$key))
    if (is.logical(a)) a <- as_sym(a)
    if (is.logical(b)) b <- as_sym(b)
    return(identical(as_sym(a)$key, as_sym(b)$key))
  }
  if (is.environment(a) || is.environment(b)) return(identical(a, b))
  isTRUE(all.equal(a, b, tolerance = 0))
}
.py_ne <- function(a, b) !.py_eq(a, b)
.cmp <- function(op) function(a, b) {
  if (is_sym(a) || is_sym(b)) {
    r <- get(op)(as_sym(a), as_sym(b))
    return(r)
  }
  get(op)(a, b)
}
.py_lt <- .cmp("<"); .py_le <- .cmp("<="); .py_gt <- .cmp(">"); .py_ge <- .cmp(">=")
.py_is <- function(a, b) {
  if (is.null(b)) return(is.null(a))
  if (is.logical(b)) return(is.logical(a) && identical(a, b))
  if (is_sym(b) && is_sym(a)) return(identical(a$key, b$key) && identical(a$type, b$type))
  .py_eq(a, b)
}
.py_isnot <- function(a, b) !.py_is(a, b)
.py_in <- function(a, b) {
  if (inherits(b, "pydict")) b <- b$keys
  if (is.character(b) && length(b) == 1L) return(grepl(a, b, fixed = TRUE))
  any(vapply(.py_iter(b), function(z) .py_eq(a, z), TRUE))
}
.py_notin <- function(a, b) !.py_in(a, b)

# ---------------------------------------------------------------- arithmetic
.plainmat <- function(x) inherits(x, "symmatrix") && identical(class(x), "symmatrix")
.num_or_sym <- function(x) if (is.numeric(x) || is.logical(x)) x else x
.py_add <- function(a, b) {
  if (.is_seq(a) && .is_seq(b)) { r <- c(unclass(a), unclass(b)); if (inherits(a, "pytuple")) r <- .tup(r); return(r) }
  if (is.character(a) && is.character(b)) return(paste0(a, b))
  a + b
}
.py_sub <- function(a, b) a - b
.py_mul <- function(a, b) {
  if (.plainmat(a) && .plainmat(b)) return(a %*% b)
  if (inherits(a, "RayTransferMatrix")) return(a %*% b)
  if (.is_seq(a) && is.numeric(b)) { r <- rep(unclass(a), b); if (inherits(a, "pytuple")) r <- .tup(r); return(r) }
  if (is.numeric(a) && is.numeric(b) && !is_sym(a) && !is_sym(b)) return(a * b)
  a * b
}
.py_div <- function(a, b) {
  if (is.numeric(a) && is.numeric(b) && !is_sym(a) && !is_sym(b)) {
    if (b == 0) stop(structure(class = c("ZeroDivisionError", "error", "condition"), list(message = "division by zero", call = NULL)))
    return(a / b)
  }
  a / b
}
.py_floordiv <- function(a, b) if (is.numeric(a) && is.numeric(b)) a %/% b else floor(a / b)
.py_mod <- function(a, b) a %% b
.py_pow <- function(a, b) {
  if (.plainmat(a)) return(a %^% b)
  if (is.numeric(a) && is.numeric(b) && !is_sym(a) && !is_sym(b)) {
    if (b >= 0 && a == round(a) && b == round(b)) return(a^b)
    return(a^b)
  }
  a^b
}
.py_matmul <- function(a, b) a %*% b
.py_neg <- function(a) -a
.py_xor <- function(a, b) xor(a, b)
.py_and <- function(a, b) a & b
.py_or <- function(a, b) a | b

# ---------------------------------------------------------------- indexing
.py_slice <- function(a, b, c) structure(list(a, b, c), class = "pyslice")
.slice_idx <- function(s, n) {
  st <- if (is.null(s[[3L]])) 1L else as.integer(s[[3L]])
  norm <- function(v, d) { if (is.null(v)) return(d); v <- as.integer(v); if (v < 0) v <- v + n; v }
  if (st > 0) { a <- max(0L, min(n, norm(s[[1L]], 0L))); b <- max(0L, min(n, norm(s[[2L]], n))) }
  else { a <- min(n - 1L, norm(s[[1L]], n - 1L)); b <- if (is.null(s[[2L]])) -1L else max(-1L, norm(s[[2L]], -1L)) }
  if ((st > 0 && a >= b) || (st < 0 && a <= b)) return(integer())
  v <- seq(a, b, by = st)
  v <- v[if (st > 0) v < b else v > b]
  v + 1L
}
.idx1 <- function(i, n) { i <- as.integer(.py_int(i)); if (i < 0) i <- i + n; if (i < 0 || i >= n) stop(structure(class = c("IndexError", "error", "condition"), list(message = "index out of range", call = NULL))); i + 1L }
.py_getitem <- function(x, k) {
  if (inherits(x, "TransferFunctionMatrix")) {
    nr <- x$num_outputs; nc <- x$num_inputs
    if (inherits(k, "pytuple")) {
      s1 <- inherits(k[[1L]], "pyslice"); s2 <- inherits(k[[2L]], "pyslice")
      ri <- if (s1) .slice_idx(k[[1L]], nr) else .idx1(k[[1L]], nr)
      ci <- if (s2) .slice_idx(k[[2L]], nc) else .idx1(k[[2L]], nc)
      return(x[ri, ci, drop = !(s1 || s2)])
    }
    i <- .idx1(k, nr * nc) - 1L
    return(x[i %/% nc + 1L, i %% nc + 1L])
  }
  if (.is_obj(x)) return(.py_method(x, "__getitem__", k))
  if (inherits(x, "pydict")) {
    i <- match(.hkey(k), vapply(x$keys, .hkey, ""))
    if (is.na(i)) stop(structure(class = c("KeyError", "error", "condition"), list(message = .py_str(k), call = NULL)))
    return(x$values[[i]])
  }
  if (inherits(x, "symmatrix")) {
    if (inherits(k, "pytuple")) {
      ri <- if (inherits(k[[1L]], "pyslice")) .slice_idx(k[[1L]], x$nrow) else .idx1(k[[1L]], x$nrow)
      ci <- if (inherits(k[[2L]], "pyslice")) .slice_idx(k[[2L]], x$ncol) else .idx1(k[[2L]], x$ncol)
      if (!inherits(k[[1L]], "pyslice") && !inherits(k[[2L]], "pyslice")) return(.I(".mget")(x, ri, ci))
      d <- list()
      for (i in ri) for (j in ci) d[[length(d) + 1L]] <- .I(".mget")(x, i, j)
      return(.I(".new_mat")(d, length(ri), length(ci)))
    }
    if (inherits(k, "pyslice")) return(x$data[.slice_idx(k, length(x$data))])
    return(x$data[[.idx1(k, length(x$data))]])
  }
  if (is_sym(x)) {
    if (x$type %in% c("Tuple")) return(.py_getitem(x$args, k))
    if (x$type %in% c("ImmutableMatrix")) return(.py_getitem(Matrix(x), k))
    stop("TypeError: object is not subscriptable")
  }
  if (is.character(x) && length(x) == 1L) {
    ch <- strsplit(x, "")[[1L]]
    if (inherits(k, "pyslice")) return(paste(ch[.slice_idx(k, length(ch))], collapse = ""))
    return(ch[[.idx1(k, length(ch))]])
  }
  if (inherits(k, "pyslice")) {
    r <- unclass(x)[.slice_idx(k, length(x))]
    if (inherits(x, "pytuple")) r <- .tup(r)
    return(r)
  }
  if (is.list(x)) return(unclass(x)[[.idx1(k, length(x))]])
  if (is.atomic(x)) return(x[[.idx1(k, length(x))]])
  stop("cannot index")
}
.py_setitem <- function(x, k, v) {
  if (.is_obj(x)) { .py_method(x, "__setitem__", k, v); return(x) }
  if (inherits(x, "pydict")) {
    i <- match(.hkey(k), vapply(x$keys, .hkey, ""))
    if (is.na(i)) { x$keys[[length(x$keys) + 1L]] <- k; x$values[[length(x$values) + 1L]] <- v }
    else x$values[[i]] <- v
    return(x)
  }
  if (inherits(x, "symmatrix")) {
    if (inherits(k, "pytuple")) {
      ri <- if (inherits(k[[1L]], "pyslice")) .slice_idx(k[[1L]], x$nrow) else .idx1(k[[1L]], x$nrow)
      ci <- if (inherits(k[[2L]], "pyslice")) .slice_idx(k[[2L]], x$ncol) else .idx1(k[[2L]], x$ncol)
      if (inherits(v, "symmatrix")) {
        for (a in seq_along(ri)) for (b in seq_along(ci)) x$data[[(ri[a] - 1L) * x$ncol + ci[b]]] <- .I(".mget")(v, a, b)
      } else for (i in ri) for (j in ci) x$data[[(i - 1L) * x$ncol + j]] <- as_sym(v)
      return(x)
    }
    x$data[[.idx1(k, length(x$data))]] <- as_sym(v)
    return(x)
  }
  cl <- class(x)
  x <- unclass(x)
  if (is.null(v)) x[.idx1(k, length(x))] <- list(NULL) else x[[.idx1(k, length(x))]] <- v
  if (!identical(cl, "list")) class(x) <- cl
  x
}

# ---------------------------------------------------------------- attributes and methods
.py_attr <- function(x, name) {
  if (.is_obj(x)) return(.I(".pyobj_get")(x, name))
  if (inherits(x, "pyclass")) return(.I(".pyobj_get")(x, name))
  if (is.complex(x)) return(switch(name, real = Re(x), imag = Im(x)))
  if (inherits(x, "symmatrix")) {
    return(switch(name,
      shape = .py_tup(x$nrow, x$ncol), rows = x$nrow, cols = x$ncol, T = t(x),
      free_symbols = .py_setlit(free_symbols(x)),
      is_square = x$nrow == x$ncol,
      if (length(class(x)) > 1L) do.call(`$`, list(x, name)) else stop("AttributeError: Matrix has no attribute ", name)))
  }
  if (is_sym(x)) {
    if (name == "args") return(.tup(.I(".sympy_args")(x)))
    if (name == "free_symbols") return(.py_setlit(free_symbols(x)))
    if (name == "func") return(x$type)
    if (name == "p") return(as_sym(gmp::numerator(x$value)))
    if (name == "q") return(as_sym(gmp::denominator(x$value)))
    if (name == "name") return(x$name)
    if (name == "lhs") return(x$args[[1L]])
    if (name == "rhs") return(x$args[[2L]])
    if (name == "exp") return(x$args[[2L]])
    if (name == "base") return(x$args[[1L]])
    if (name %in% c("is_Float", "is_Rational", "is_Integer", "is_Symbol", "is_Add", "is_Mul", "is_Pow", "is_Number", "is_Function", "is_Matrix")) {
      t <- sub("^is_", "", name)
      return(switch(t, Number = x$type %in% c("Integer", "Rational", "Float"), Rational = x$type %in% c("Integer", "Rational"),
                    Function = x$type == "Function", Matrix = FALSE, x$type == t))
    }
    if (name == "is_number") return(length(free_symbols(x)) == 0L && !x$type %in% c("Symbol"))
    if (startsWith(name, "is_")) {
      r <- .I("ask_fact")(x, sub("^is_", "", name))
      return(if (is.na(r)) NULL else r)
    }
  }
  stop("AttributeError: no attribute ", name)
}
.py_setattr <- function(x, name, v) {
  if (.is_obj(x)) return(invisible(.I(".pyobj_set")(x, name, v)))
  stop("cannot set attribute")
}
.py_method <- function(x, name, ...) {
  if (is.atomic(x) && !is.null(x) && !is.character(x)) {
    if (name == "any") return(any(x))
    if (name == "astype") return(list(...)[[1L]](x))
  }
  if (.is_obj(x) || inherits(x, "pyclass")) {
    f <- .I(".pyobj_get")(x, name)
    if (name == "rewrite") return(.py_rewrite_obj(f, ...))
    return(f(...))
  }
  a <- list(...)
  if (inherits(x, "symmatrix")) {
    return(switch(name,
      subs = .py_subs(x, ...), doit = doit(x), simplify = simplify(x), expand = expand(x, ...),
      inv = inverse(x), det = sym_det(x), rank = rank_matrix(x), transpose = t(x),
      tolist = lapply(.I(".mrows")(x), function(r) r), applyfunc = .I(".mat_map")(x, function(z) as_sym(a[[1L]](z))),
      evalf = , n = .py_N(x, ...), nullspace = nullspace(x), columnspace = columnspace(x),
      row_join = cbind(x, a[[1L]]), col_join = rbind(x, a[[1L]]), copy = x,
      eigenvals = { ev <- eigenvals(x); .I(".pydict")(unclass(unname(ev))[seq_along(ev)], as.list(attr(ev, "multiplicity"))) },
      equals = .py_mat_equals(x, a[[1L]]),
      as_real_imag = .tup(list(.I(".mat_map")(x, re), .I(".mat_map")(x, im))), replace = .I(".mat_map")(x, function(z) .py_method(z, "replace", ...)),
      stop("AttributeError: Matrix method ", name)))
  }
  if (is.list(x) && !is.object(x)) {
    return(switch(name, index = match(.hkey(a[[1L]]), vapply(x, .hkey, "")) - 1L,
                  count = sum(vapply(x, function(z) .py_eq(z, a[[1L]]), TRUE)),
                  stop("AttributeError: list method ", name)))
  }
  if (inherits(x, "pydict")) {
    return(switch(name, keys = x$keys, values = x$values, items = lapply(seq_along(x$keys), function(i) .tup(list(x$keys[[i]], x$values[[i]]))),
                  get = { i <- match(.hkey(a[[1L]]), vapply(x$keys, .hkey, "")); if (is.na(i)) (if (length(a) > 1L) a[[2L]] else NULL) else x$values[[i]] },
                  stop("AttributeError: dict method ", name)))
  }
  if (is.character(x)) return(switch(name, format = .py_strformat(x, ...), stop("str method")))
  x <- as_sym(x)
  switch(name,
    subs = .py_subs(x, ...), xreplace = .py_subs(x, ...),
    doit = doit(x), simplify = simplify(x, ...), expand = expand(x, ...), factor = sym_factor(x, ...),
    cancel = cancel(x, ...), together = together(x), convert_to = convert_to(x, ...), evalf = , n = .py_N(x, ...),
    rewrite = .py_rewrite(x, ...), as_real_imag = .tup((if (exists("as_real_imag", envir = .ns)) .I("as_real_imag") else .I(".as_real_imag"))(x)),
    as_numer_denom = .tup(fraction(x)), has = .I(".has1")(x, as_sym(a[[1L]])),
    equals = .py_equals(x, a[[1L]]), diff = .py_diff(x, ...), integrate = sym_integrate(x, ...),
    coeff = coeff(x, ...), collect = collect(x, ...), atoms = .py_setlit(atoms(x, ...)),
    replace = sym_replace(x, ...), trigsimp = trigsimp(x), radsimp = radsimp(x), nsimplify = nsimplify(x, ...),
    as_poly = Poly(x, ...), all_coeffs = .I("all_coeffs")(x), is_constant = .I("is_constant")(x, ...),
    could_extract_minus_sign = .I(".could_extract_minus_sign")(x), count_ops = count_ops(x),
    round = .py_round_sym(x, ...), limit = limit(x, ...), series = series(x, ...),
    is_same = {
      o <- as_sym(a[[1L]])
      if (length(a) > 1L && x$type == "Float" && o$type == "Float") a[[2L]](as.numeric(N(x)), as.numeric(N(o)))
      else identical(x$key, o$key)
    },
    stop("AttributeError: sym method ", name))
}
.py_round_sym <- function(x, n = 0) round(x, n)
# MatrixBase.equals: False on a shape mismatch or a False entry, else True
# when every entry is equal
.py_mat_equals <- function(a, b) {
  if (!inherits(b, "symmatrix") || a$nrow != b$nrow || a$ncol != b$ncol) return(FALSE)
  rv <- TRUE
  for (i in seq_len(a$nrow)) for (j in seq_len(a$ncol)) {
    ans <- .py_equals(.I(".mget")(a, i, j), .I(".mget")(b, i, j))
    if (isFALSE(ans)) return(FALSE)
  }
  rv
}
# Expr.equals: True when the simplified difference is a zero number
.py_equals <- function(a, b) { d <- simplify(.py_sub(a, b)); identical(d$key, S_Zero$key) || (is_sym(d) && d$type == "Float" && d$value == 0) }
.py_subs <- function(x, ...) {
  a <- list(...)
  if (length(a) == 1L && inherits(a[[1L]], "pydict")) {
    d <- a[[1L]]
    pairs <- lapply(seq_along(d$keys), function(i) list(d$keys[[i]], d$values[[i]]))
    return(subs(x, pairs))
  }
  if (length(a) == 1L && .is_seq(a[[1L]])) {
    pairs <- lapply(.py_iter(a[[1L]]), function(p) as.list(.py_iter(p)))
    return(subs(x, pairs))
  }
  subs(x, a[[1L]], a[[2L]])
}
.py_rewrite_obj <- function(f, target, ...) {
  if (is.function(target)) target <- attr(target, "pyname") %||% .fname(target)
  f(target, ...)
}
.py_rewrite <- function(x, target, ...) {
  if (is.function(target)) target <- attr(target, "pyname") %||% .fname(target)
  rewrite(x, target)
}
`%||%` <- function(a, b) if (is.null(a)) b else a
.fname <- function(f) {
  for (nm in c("sin", "cos", "exp", "tan", "sqrt", "log")) if (identical(f, get(paste0(".py_", nm)))) return(nm)
  stop("unknown rewrite target")
}
.py_unsupported <- function(msg) stop("UNSUPPORTED ", msg)
.py_invert <- function(x) Not(x)
.py_strformat <- function(x, ...) x

# ---------------------------------------------------------------- formatting
.pyrepr_float <- function(x) sprintf("%.12g", x)
.fmt <- function(v, py_str = FALSE) {
  if (is.null(v)) return("None")
  if (inherits(v, "pydict")) {
    if (!length(v$keys)) return("{}")
    return(paste0("{", paste(vapply(seq_along(v$keys), function(i) paste0(.fmt(v$keys[[i]]), ": ", .fmt(v$values[[i]])), ""), collapse = ", "), "}"))
  }
  if (inherits(v, "pyset")) {
    if (!length(v)) return("set()")
    return(paste0("{", paste(sort(vapply(unclass(v), .fmt, ""), method = "radix"), collapse = ", "), "}"))
  }
  if (inherits(v, "pytuple")) {
    s <- vapply(unclass(v), .fmt, "")
    return(paste0("(", paste(s, collapse = ", "), if (length(s) == 1L) ",", ")"))
  }
  # Python fmt(): a matrix in a tuple or a list prints on one line too
  if (is_sym(v) && identical(v$type, "ImmutableMatrix")) return(.fmt_top(v))
  if (inherits(v, "symmatrix") && v$nrow > 0L && v$ncol > 0L) return(.fmt_top(v))
  if (is_sym(v) || inherits(v, "symmatrix")) return(sstr(v))
  if (inherits(v, "Poly")) return(paste(utils::capture.output(print(v)), collapse = "\n"))
  if (.is_obj(v)) return(.I(".pyobj_str")(v))
  if (is.character(v) && length(v) == 1L) return(if (py_str) v else paste0("'", v, "'"))
  if (is.logical(v) && length(v) == 1L) return(if (is.na(v)) "None" else if (v) "True" else "False")
  if (is.integer(v) && length(v) == 1L) return(as.character(v))
  if (is.double(v) && length(v) == 1L) return(.pyrepr_float(v))
  if (is.complex(v) && length(v) == 1L) return(sprintf("(%.12g%+.12gj)", Re(v), Im(v)))
  if (is.list(v) && !is.object(v)) return(paste0("[", paste(vapply(v, .fmt, ""), collapse = ", "), "]"))
  if (is.atomic(v) && length(v) != 1L) return(paste0("[", paste(vapply(as.list(v), .fmt, ""), collapse = ", "), "]"))
  paste(format(v), collapse = " ")
}

# Python's str() of a top level Matrix: str(M.tolist()), so each element is
# printed on its own (a bare Float at full precision)
.fmt_top <- function(v) {
  if (is_sym(v) && identical(v$type, "ImmutableMatrix")) v <- .I(".lti_mat")(v)
  if (inherits(v, "symmatrix") && v$nrow > 0L && v$ncol > 0L) {
    rows <- vapply(seq_len(v$nrow), function(i)
      paste0("[", paste(vapply(seq_len(v$ncol), function(j) sstr(.I(".mget")(v, i, j)), ""), collapse = ", "), "]"), "")
    return(paste0("Matrix([", paste(rows, collapse = ", "), "])"))
  }
  .fmt(v)
}

.err_class <- function(e) {
  cl <- setdiff(class(e), c("simpleError", "error", "condition", "sympy_error", "rlang_error"))
  if (length(cl)) cl[1L] else "Error"
}

# ---------------------------------------------------------------- recorder
.rec_state <- new.env()
.py_rec <- function(cid, thunk) {
  n <- .rec_state$counts[[cid]]
  n <- if (is.null(n)) 1L else n + 1L
  .rec_state$counts[[cid]] <- n
  key <- if (n == 1L) cid else paste0(cid, "#", n)
  res <- if (cid %in% .rec_state$fork_ids) .py_rec_fork(thunk) else tryCatch({
    setTimeLimit(cpu = .rec_state$cpu, transient = TRUE)
    v <- suppressWarnings(thunk())
    setTimeLimit(cpu = Inf)
    .fmt_top(v)
  }, error = function(e) {
    setTimeLimit(cpu = Inf)
    msg <- conditionMessage(e)
    if (grepl("reached (CPU|elapsed) time limit", msg)) "TIMEOUT" else paste0("ERR ", .err_class(e), ": ", msg)
  })
  res <- gsub("**", "^", res, fixed = TRUE)
  res <- gsub("\n", "\\n", res, fixed = TRUE)
  res <- gsub("\t", " ", res, fixed = TRUE)
  writeLines(paste0(key, "\t", res), .rec_state$out)
  flush(.rec_state$out)
  invisible(NULL)
}
# cases known to swallow the time limit error inside package code run in a
# forked child that is killed when the limit passes
.rec_state$fork_ids <- character()
.py_rec_fork <- function(thunk) {
  job <- parallel::mcparallel(tryCatch(.fmt_top(suppressWarnings(thunk())),
    error = function(e) paste0("ERR ", .err_class(e), ": ", conditionMessage(e))))
  r <- parallel::mccollect(job, wait = FALSE, timeout = .rec_state$cpu * 4)
  if (is.null(r)) {
    tools::pskill(job$pid, tools::SIGKILL)
    parallel::mccollect(job, wait = TRUE)
    return("TIMEOUT")
  }
  r[[1L]]
}
.py_try <- function(expr) {
  tryCatch({
    setTimeLimit(cpu = .rec_state$cpu * 2, transient = TRUE)
    force(expr)
    setTimeLimit(cpu = Inf)
    invisible(NULL)
  }, error = function(e) {
    setTimeLimit(cpu = Inf)
    .rec_state$fails[[length(.rec_state$fails) + 1L]] <- conditionMessage(e)
    invisible(NULL)
  })
}

# class names of the control module that only appear in isinstance()
LinearTimeInvariant <- "LinearTimeInvariant"
SISOLinearTimeInvariant <- "SISOLinearTimeInvariant"
MIMOLinearTimeInvariant <- "MIMOLinearTimeInvariant"
StateSpaceBase <- "StateSpaceBase"

# list.append(), extend(), insert(), remove() and clear() give the new list
.py_list_mut <- function(x, op, ...) {
  a <- list(...)
  if (.is_obj(x) || inherits(x, "pydict")) { .py_method(x, op, ...); return(x) }
  l <- if (is.null(x)) list() else unclass(x)
  if (!is.list(l)) l <- as.list(l)
  l <- switch(op,
    append = c(l, list(a[[1L]])),
    extend = c(l, .py_iter(a[[1L]])),
    insert = { i <- as.integer(a[[1L]]); if (i < 0) i <- max(0L, length(l) + i); append(l, list(a[[2L]]), after = min(i, length(l))) },
    remove = {
      i <- which(vapply(l, function(e) isTRUE(.py_truthy(.py_eq(e, a[[1L]]))), TRUE))
      if (!length(i)) stop(structure(class = c(if (inherits(x, "pyset")) "KeyError" else "ValueError", "error", "condition"),
                                     list(message = "list.remove(x): x not in list", call = NULL)))
      l[-i[[1L]]]
    },
    clear = list())
  if (is.object(x)) class(l) <- class(x)
  l
}

# the builtin pow()
pow <- function(a, b, m = NULL) if (is.null(m)) .py_pow(a, b) else .py_mod(.py_pow(a, b), m)
]---"), envir = env)
    }
    env
  }
})

# run the statements of one case block and return its recorded rows, named
# by case id
.physmisc_battery_run <- function(code) {
  ph <- .physmisc_battery_env()
  lines <- character()
  con <- textConnection("lines", "w", local = TRUE)
  on.exit(close(con))
  st <- ph$.rec_state
  st$out <- con
  st$cpu <- 60
  st$counts <- list()
  st$fails <- list()
  st$fork_ids <- character()
  env <- new.env(parent = ph)
  for (e in parse(text = code, keep.source = FALSE)) try(eval(e, env), silent = TRUE)
  kv <- strsplit(textConnectionValue(con), "\t", fixed = TRUE)
  stats::setNames(vapply(kv, function(p) paste(p[-1L], collapse = "\t"), ""), vapply(kv, `[`, "", 1L))
}
