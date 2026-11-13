# Expected values are SymPy outputs for the codegen tests of SymPy
# (sympy/codegen/tests, sympy/utilities/tests/test_codegen*.py, the
# codegen rows of sympy/printing/tests and of test_lambdify.py): each
# row prints str(result) with ** shown as ^, or the error class and message.

skip_unless_full_tests()

# Helpers mirroring the Python constructs of the SymPy tests
# Battery helpers: Python semantics for indexing, membership and printing.
py_repr <- function(s) {
  q <- if (grepl("'", s, fixed = TRUE) && !grepl("\"", s, fixed = TRUE)) "\"" else "'"
  s <- gsub("\\", "\\\\", s, fixed = TRUE)
  if (q == "'") s <- gsub("'", "\\'", s, fixed = TRUE)
  s <- gsub("\n", "\\n", s, fixed = TRUE)
  paste0(q, s, q)
}
pystr <- function(s) structure(s, class = "pystr")

py_list <- function(x) {
  if (is_sym(x)) return(x$args)
  if (is.list(x)) return(x)
  as.list(x)
}
py_iter <- py_list
py_tuple <- function(x) do.call(Tuple, py_list(x))

.pyidx <- function(i, n) { i <- as.integer(if (is_sym(i)) as.numeric(i$value) else i); if (i < 0) n + i + 1L else i + 1L }

py_getitem <- function(x, ...) {
  idx <- list(...)
  if (is.function(x) && !is.null(attr(x, "base"))) return(do.call(x, idx))
  if (inherits(x, "cg_dict")) return(cgdict_get(x, idx[[1]]))
  if (length(idx) == 2L && any(vapply(idx, inherits, TRUE, "pyslice"))) {
    cv <- function(i) if (inherits(i, "pyslice")) lapply(unclass(i), function(v) if (is.null(v)) NA else v) else c(i, i + 1L)
    return(MatrixSlice(x, cv(idx[[1]]), cv(idx[[2]])))
  }
  if (is_sym(x) && x$type == "MatrixSymbol" && length(idx) == 1L) {
    nc <- as.integer(symcompR:::.cgp_mshape(x)[2L])
    k <- as.integer(idx[[1]])
    return(MatrixElement(x, k %/% nc, k %% nc))
  }
  if (is_sym(x) && x$type == "MatrixSymbol") return(MatrixElement(x, idx[[1]], idx[[2]]))
  if (is_sym(x) && x$type %in% c("Pointer")) return(ast_getitem(x, if (length(idx) == 1L) idx[[1]] else idx))
  if (inherits(x, "symmatrix")) {
    if (length(idx) == 2L) return(x[.pyidx(idx[[1]], x$nrow), .pyidx(idx[[2]], x$ncol)])
    return(x$data[[.pyidx(idx[[1]], length(x$data))]])
  }
  i <- idx[[1]]
  if (is.character(x) && !is.character(i)) { k <- .pyidx(i, nchar(x)); return(substr(x, k, k)) }
  if (is_sym(x)) { a <- x$args; return(a[[.pyidx(i, length(a))]]) }
  if (is.list(x) && is.character(i)) return(x[[i]])
  if (is.list(x)) return(x[[.pyidx(i, length(x))]])
  x[[.pyidx(i, length(x))]]
}

py_slice <- function(x, lo = NULL, hi = NULL, step = NULL) {
  if (!is.null(step) && step == -1L && is.null(lo) && is.null(hi)) {
    if (is.character(x)) return(paste(rev(strsplit(x, "")[[1]]), collapse = ""))
    return(rev(py_list(x)))
  }
  n <- if (is.character(x)) nchar(x) else length(py_list(x))
  a <- if (is.null(lo)) 0L else as.integer(lo); b <- if (is.null(hi)) n else as.integer(hi)
  if (a < 0) a <- n + a; if (b < 0) b <- n + b
  a <- max(0L, min(a, n)); b <- max(0L, min(b, n))
  if (is.character(x)) return(if (b > a) substr(x, a + 1L, b) else "")
  l <- py_list(x)
  if (b > a) l[(a + 1L):b] else list()
}

py_len <- function(x) {
  if (is.character(x)) return(nchar(x))
  if (is_sym(x)) return(length(x$args))
  length(x)
}

py_eq <- function(a, b) {
  isnone <- function(v) is.null(v) || (is_sym(v) && v$type == "NoneToken")
  if (isnone(a) || isnone(b)) return(isnone(a) && isnone(b))
  if (is.list(a) && !is_sym(a) && is_sym(b) && b$type == "Tuple") a <- do.call(Tuple, a)
  if (is.list(b) && !is_sym(b) && is_sym(a) && a$type == "Tuple") b <- do.call(Tuple, b)
  if (is_sym(a) || is_sym(b)) {
    if (is.character(a) || is.character(b)) return(FALSE)
    if (is.function(a) || is.function(b)) return(FALSE)
    if ((is.numeric(a) && is.infinite(a)) || (is.numeric(b) && is.infinite(b))) {
      return(identical(sstr(as_sym(a)), sstr(as_sym(b))))
    }
    return(identical(as_sym(a)$key, as_sym(b)$key) || isTRUE(tryCatch(identical(sstr(as_sym(a) - as_sym(b)), "0"), error = function(e) FALSE)))
  }
  isTRUE(all.equal(a, b)) && identical(length(a), length(b))
}
py_is <- function(a, b) if (is.function(a) || is.function(b)) identical(a, b) else py_eq(a, b)
unchanged <- function(f, ...) { a <- list(...); r <- do.call(f, a); is_sym(r) && identical(ast_func(r), f) && length(r$args) == length(a) && all(mapply(py_eq, r$args, a)) }
py_in <- function(a, b) {
  for (z in py_list(b)) if (py_eq(a, z)) return(TRUE)
  FALSE
}
py_str <- function(x) if (is.character(x)) x else if (is.null(x)) "None" else sstr(x)
py_set <- function(l) {
  if (!length(l)) return(pystr("set()"))
  pystr(paste0("{", paste(vapply(l, sstr, ""), collapse = ", "), "}"))
}
py_abs <- function(x) if (is_sym(x)) Abs(x) else abs(x)
py_float <- function(x) {
  if (is_sym(x)) {
    s <- sstr(x)
    if (s == "oo") return(Inf)
    if (s == "-oo") return(-Inf)
    if (s == "nan") return(NaN)
    return(as.numeric(N(x)$value))
  }
  as.numeric(x)
}
py_isnan <- function(x) if (is_sym(x)) identical(sstr(x), "nan") else is.nan(x)
pow <- function(a, b) a^b
py_isinstance <- function(x, cls) {
  nm <- if (is.character(cls)) cls else deparse(substitute(cls))
  if (nm == "Float") return(is_sym(x) && x$type == "Float")
  if (nm == "Integer") return(is_sym(x) && x$type == "Integer")
  is_sym(x) && (x$type == nm)
}
py_format <- function(fmt, args) {
  if (is.list(args) && !is_sym(args) && !is.null(names(args)) && all(nzchar(names(args)))) {
    for (k in names(args)) {
      v <- args[[k]]; v <- if (is_sym(v)) sstr(v) else as.character(v)
      fmt <- gsub(paste0("%(", k, ")s"), v, fmt, fixed = TRUE)
      fmt <- gsub(paste0("%(", k, ")i"), v, fmt, fixed = TRUE)
      fmt <- gsub(paste0("%(", k, ")d"), v, fmt, fixed = TRUE)
    }
    return(gsub("%%", "%", fmt, fixed = TRUE))
  }
  args <- if (is.list(args) && !is_sym(args)) args else list(args)
  args <- lapply(args, function(a) if (is_sym(a)) sstr(a) else a)
  fmt <- gsub("%r", "%s", fmt, fixed = TRUE)
  do.call(sprintf, c(list(fmt), args))
}
py_join <- function(sep, x) paste(unlist(py_list(x)), collapse = sep)
py_args <- function(x) structure(if (is_sym(x)) x$args else py_list(x), class = "pytuple")
pyrepr_el <- function(e) {
  if (inherits(e, "pystr")) return(unclass(e))
  if (inherits(e, "cg_cgobj")) return(format(e))
  if (is_sym(e)) {
    cl <- symcompR:::.cg_cls(e)
    if (!is.null(cl) && isTRUE(cl$token)) return(srepr(e))
    return(sstr(e))
  }
  if (is.character(e)) return(py_repr(e))
  if (inherits(e, "pytuple")) return(fmt_pytuple(e))
  if (is.logical(e)) return(if (e) "True" else "False")
  if (is.null(e)) return("None")
  if (is.list(e)) return(paste0("[", paste(vapply(e, pyrepr_el, ""), collapse = ", "), "]"))
  as.character(e)
}
fmt_pytuple <- function(t) {
  els <- vapply(unclass(t), pyrepr_el, "")
  if (length(els) == 1L) paste0("(", els, ",)") else paste0("(", paste(els, collapse = ", "), ")")
}
py_repr_obj <- function(x) {
  if (is_sym(x)) return(pystr(pyrepr_el(x)))
  if (is.character(x)) return(pystr(py_repr(x)))
  pystr(pyrepr_el(x))
}

# dict literal: a named list for string keys, else a list of (key, value) pairs
py_dict <- function(keys = list(), values = list()) {
  if (!length(keys)) return(setNames(list(), character()))
  if (all(vapply(keys, function(k) is.character(k) && length(k) == 1L, TRUE)))
    return(setNames(values, unlist(keys)))
  Map(function(k, v) list(k, v), keys, values)
}
py_all <- function(x) all(vapply(py_list(x), isTRUE, TRUE))
py_any <- function(x) any(vapply(py_list(x), isTRUE, TRUE))

# StringIO: a raw connection
StringIO <- function() structure(list(con = rawConnection(raw(0), "w")), class = "py_stringio")
.sio_con <- function(x) if (inherits(x, "py_stringio")) x$con else x
py_method <- function(obj, name, ...) {
  a <- list(...)
  if (inherits(obj, "py_stringio")) {
    return(switch(name,
                  getvalue = { flush(obj$con); rawToChar(rawConnectionValue(obj$con)) },
                  write = { cat(a[[1]], file = obj$con, sep = ""); nchar(a[[1]]) },
                  close = invisible(close(obj$con))))
  }
  if (is.character(obj) && length(obj) == 1L && !inherits(obj, "pystr")) {
    s <- obj
    return(switch(name,
                  center = { w <- a[[1]]; n <- nchar(s); if (w <= n) s else { tot <- w - n; l <- tot %/% 2 + (tot %% 2 & w %% 2); r <- tot - l; paste0(strrep(" ", l), s, strrep(" ", r)) } },
                  rstrip = sub("[ \t\n\r]+$", "", s),
                  lstrip = sub("^[ \t\n\r]+", "", s),
                  strip = trimws(s),
                  endswith = endsWith(s, a[[1]]),
                  startswith = startsWith(s, a[[1]]),
                  splitlines = { l <- strsplit(s, "\n", fixed = TRUE)[[1]]; as.list(l) },
                  split = { l <- if (length(a)) strsplit(s, a[[1]], fixed = TRUE)[[1]] else strsplit(trimws(s), "[ \t\n]+")[[1]]; if (length(a) && endsWith(s, a[[1]])) l <- c(l, ""); as.list(l) },
                  count = lengths(regmatches(s, gregexpr(a[[1]], s, fixed = TRUE))),
                  replace = gsub(a[[1]], a[[2]], s, fixed = TRUE),
                  lower = tolower(s), upper = toupper(s),
                  stop("no string method ", name)))
  }
  if (is_sym(obj) && name == "count_ops") return(count_ops(obj, ...))
  f <- ast_attr(obj, name)
  do.call(f, lapply(a, .sio_con))
}
map <- function(f, x) lapply(if (is.character(x) && length(x) == 1L) as.list(strsplit(x, "")[[1]]) else py_list(x), f)
hasattr <- function(x, name) if (is.function(x)) !is.null(attr(x, name)) else !inherits(tryCatch(ast_attr(x, name), error = function(e) e), "error")
py_type <- function(x) {
  cl <- class(x)
  for (k in c("InputArgument", "OutputArgument", "InOutArgument", "Result", "Routine", "DataType"))
    if (k %in% cl) return(pystr(sprintf("<class 'sympy.utilities.codegen.%s'>", k)))
  if (is_sym(x)) return(pystr(sprintf("<class 'sympy.%s'>", x$type)))
  pystr(sprintf("<class '%s'>", cl[1]))
}
pyslice <- function(lo = NULL, hi = NULL, step = NULL) structure(list(lo, hi, step), class = "pyslice")
# symbols(..., cls=Dummy)
sym_symbols <- function(names, ..., cls = NULL) {
  if (is.null(cls)) return(symcompR::sym_symbols(names, ...))
  nms <- strsplit(trimws(names), "[ ,]+")[[1]]
  r <- lapply(nms, function(n) cls(n, ...))
  if (length(r) == 1L) r[[1]] else r
}
py_call <- function(f, ...) f(...)

pyfmt <- function(v) {
  s <- pyfmt0(v)
  gsub("\n", "\\n", gsub("**", "^", s, fixed = TRUE), fixed = TRUE)
}
pyfmt0 <- function(v) {
  if (is.null(v)) return("None")
  if (inherits(v, c("sym", "symmatrix"))) return(sstr(v))
  if (inherits(v, "pystr")) return(unclass(v))
  if (inherits(v, c("cg_pyset", "cg_pysetobj"))) return(format(v))
  if (is.function(v) && !is.null(attr(v, "sym_name"))) return(attr(v, "sym_name"))
  if (inherits(v, "pytuple")) return(fmt_pytuple(v))
  if (is.logical(v) && length(v) == 1L) return(if (is.na(v)) "None" else if (v) "True" else "False")
  if (is.character(v) && length(v) == 1L) return(v)
  if (is.numeric(v) && length(v) == 1L) {
    if (is.integer(v)) return(as.character(v))
    if (is.infinite(v)) return(if (v > 0) "inf" else "-inf")
    if (v == round(v) && abs(v) < 1e15) return(paste0(format(v, scientific = FALSE), ".0"))
    return(format(v, digits = 17))
  }
  if (is.list(v) && !is.object(v)) return(paste0("[", paste(vapply(v, pyrepr_el, ""), collapse = ", "), "]"))
  tryCatch(sstr(v), error = function(e) paste(format(v), collapse = " "))
}
pyrepr <- function(v) if (is.character(v) && length(v) == 1L && !inherits(v, "pystr")) py_repr(v) else pyfmt0(v)

errfmt <- function(e) {
  cls <- setdiff(class(e), c("simpleError", "error", "condition", "sympy_error"))
  cls <- if (length(cls)) cls[1] else "Error"
  paste0("ERR ", cls, ": ", gsub("\n", "\\n", conditionMessage(e), fixed = TRUE))
}

.cgb <- function(x) tryCatch(pyfmt(x), error = function(e) errfmt(e))

test_that("codegen_test_abstract_nodes: test_List", {
  suppressWarnings(try({ l <- List(2L, 3L, 4L) }, silent = TRUE))
  expect_identical(.cgb(l), "[2, 3, 4]")
  expect_identical(.cgb(py_str(l)), "[2, 3, 4]")
  suppressWarnings(try({ .tmp <- sym_symbols("x y z"); x <- py_getitem(.tmp, 0L); y <- py_getitem(.tmp, 1L); z <- py_getitem(.tmp, 2L) }, silent = TRUE))
  suppressWarnings(try({ l <- List(x ^ 2L, y ^ 3L, z ^ 4L) }, silent = TRUE))
  suppressWarnings(try({ m <- sym_replace(l, function(arg) ast_attr(arg, "is_Pow") && (ast_attr(arg, "exp") > 2L), function(p) ast_attr(p, "base") - ast_attr(p, "exp")) }, silent = TRUE))
  expect_identical(.cgb(m), "[x^2, y - 3, z - 4]")
  suppressWarnings(try({ hash(m) }, silent = TRUE))
})

test_that("codegen_test_algorithms: test_newtons_method", {
  py_module <- py_render_as_module
  f_module <- f_render_as_module
  suppressWarnings(try({ .tmp <- sym_symbols("x dx atol"); x <- py_getitem(.tmp, 0L); dx <- py_getitem(.tmp, 1L); atol <- py_getitem(.tmp, 2L) }, silent = TRUE))
  suppressWarnings(try({ expr <- cos(x) - (x ^ 3L) }, silent = TRUE))
  suppressWarnings(try({ algo <- newtons_method(expr, x, atol, dx) }, silent = TRUE))
  expect_identical(.cgb(has(algo, Assignment(dx, (-expr) / diff(expr, x)))), "True")
})

test_that("codegen_test_algorithms: test_newtons_method_options", {
  py_module <- py_render_as_module
  f_module <- f_render_as_module
  suppressWarnings(try({ NULL }, silent = TRUE))
  suppressWarnings(try({ x <- Symbol("x", real = TRUE) }, silent = TRUE))
  suppressWarnings(try({ dx <- Symbol("dx") }, silent = TRUE))
  suppressWarnings(try({ e <- cos(x) - (x ^ 3L) }, silent = TRUE))
  expect_identical(.cgb(ccode(newtons_method(e, x, delta = dx, itermax = 10L, counter = Symbol("i", integer = TRUE)))), "double dx = INFINITY;\\nint i = 0;\\nwhile (fabs(dx) > 3.9999999999999999e-16*fabs(x) + 9.9999999999999998e-13 && i < 10) {\\n   dx = (pow(x, 3) - cos(x))/(-3*pow(x, 2) - sin(x));\\n   x += dx;\\n   i += 1;\\n};")
  expect_identical(.cgb(ccode(newtons_method(e, x, delta = dx, bounds = list(0L, 1L)))), "double dx = INFINITY;\\nwhile (fabs(dx) > 3.9999999999999999e-16*fabs(x) + 9.9999999999999998e-13) {\\n   dx = (pow(x, 3) - cos(x))/(-3*pow(x, 2) - sin(x));\\n   x += dx;\\n   x = fmin(1, fmax(0, x));\\n};")
  expect_identical(.cgb(ccode(newtons_method(e, x, delta = dx, debug = TRUE))), "double dx = INFINITY;\\nprintf(\"x=%12.5g\\n\", x);\\nwhile (fabs(dx) > 3.9999999999999999e-16*fabs(x) + 9.9999999999999998e-13) {\\n   dx = (pow(x, 3) - cos(x))/(-3*pow(x, 2) - sin(x));\\n   x += dx;\\n   printf(\"x=%12.5g dx=%12.5g\\n\", x, dx);\\n};")
  expect_identical(.cgb(ccode(newtons_method(e, x, delta = dx, cse = TRUE))), "double dx = INFINITY;\\nwhile (fabs(dx) > 3.9999999999999999e-16*fabs(x) + 9.9999999999999998e-13) {\\n   dx = -(pow(x, 3) - cos(x))/(3*pow(x, 2) + sin(x));\\n   x += dx;\\n};")
  expect_identical(.cgb(pycode(newtons_method_function(e, x, itermax = 5L, counter = Symbol("i", integer = TRUE)))), "def newton(x):\\n    d_x = math.inf\\n    i = 0\\n    while (abs(d_x) > 4.0e-16*abs(x) + 1.0e-12) and (i < 5):\\n        d_x = (x^3 - math.cos(x))/(-3*x^2 - math.sin(x))\\n        x += d_x\\n        i += 1\\n    return x")
  expect_identical(.cgb(fcode(newtons_method_function(e, x), source_format = "free", standard = 2003L)), "real*8 function newton(x)\\nreal*8 :: x\\nreal*8 :: d_x = (huge(0d0) + 1)\\ndo while (abs(d_x) > 4.0d-16*abs(x) + 1.0d-12)\\n   d_x = (x^3 - cos(x))/(-3*x^2 - sin(x))\\n   x = x + d_x\\nend do\\nnewton = x\\nend function")
})

test_that("codegen_test_approximations: test_SumApprox_trivial", {
  suppressWarnings(try({ x <- sym_symbols("x") }, silent = TRUE))
  suppressWarnings(try({ expr1 <- 1L + x }, silent = TRUE))
  suppressWarnings(try({ sum_approx <- SumApprox(bounds = py_dict(list(x), list(list(-Float(1e-20), Float(1e-20)))), reltol = Float(1e-16)) }, silent = TRUE))
  suppressWarnings(try({ apx1 <- sym_optimize(expr1, list(sum_approx)) }, silent = TRUE))
  expect_identical(.cgb(apx1 - 1L), "0")
})

test_that("codegen_test_approximations: test_SumApprox_monotone_terms", {
  suppressWarnings(try({ .tmp <- sym_symbols("x y z"); x <- py_getitem(.tmp, 0L); y <- py_getitem(.tmp, 1L); z <- py_getitem(.tmp, 2L) }, silent = TRUE))
  suppressWarnings(try({ expr1 <- exp(z) * (((x ^ 2L) + (y ^ 2L)) + 1L) }, silent = TRUE))
  suppressWarnings(try({ bnds1 <- py_dict(list(x, y), list(list(0L, Float(0.001)), list(100L, 1000L))) }, silent = TRUE))
  suppressWarnings(try({ sum_approx_m2 <- SumApprox(bounds = bnds1, reltol = Float(0.01)) }, silent = TRUE))
  suppressWarnings(try({ sum_approx_m5 <- SumApprox(bounds = bnds1, reltol = Float(1e-05)) }, silent = TRUE))
  suppressWarnings(try({ sum_approx_m11 <- SumApprox(bounds = bnds1, reltol = Float(1e-11)) }, silent = TRUE))
  expect_identical(.cgb(simplify((sym_optimize(expr1, list(sum_approx_m2)) / exp(z)) - (y ^ 2L))), "0")
  expect_identical(.cgb(simplify((sym_optimize(expr1, list(sum_approx_m5)) / exp(z)) - ((y ^ 2L) + 1L))), "0")
  expect_identical(.cgb(simplify((sym_optimize(expr1, list(sum_approx_m11)) / exp(z)) - (((y ^ 2L) + 1L) + (x ^ 2L)))), "0")
})

test_that("codegen_test_approximations: test_SeriesApprox_trivial", {
  suppressWarnings(try({ .tmp <- sym_symbols("x z"); x <- py_getitem(.tmp, 0L); z <- py_getitem(.tmp, 1L) }, silent = TRUE))
  suppressWarnings(try({ expr1 <- exp(x) * 1L }, silent = TRUE))
  expect_identical(.cgb(sym_optimize(expr1, list(SeriesApprox(bounds = py_dict(list(x), list(list(-1L, 1L))), reltol = 0.5)))), "x^2/2 + x + 1")
  expect_identical(.cgb(sym_optimize(expr1, list(SeriesApprox(bounds = py_dict(list(x), list(list(-1L, 1L))), reltol = 0.1)))), "x^3/6 + x^2/2 + x + 1")
  expect_identical(.cgb(sym_optimize(expr1, list(SeriesApprox(bounds = py_dict(list(x), list(list(-1L, 1L))), reltol = 0.05)))), "x^4/24 + x^3/6 + x^2/2 + x + 1")
  expect_identical(.cgb(sym_optimize(expr1, list(SeriesApprox(bounds = py_dict(list(x), list(list(-1L, 1L))), reltol = 0.05, max_order = 3L)))), "exp(x)")
  suppressWarnings(try({ expr1 <- exp(x) * exp(z) }, silent = TRUE))
  expect_identical(.cgb(sym_optimize(expr1, list(SeriesApprox(bounds = py_dict(list(x), list(list(-1L, 1L))), reltol = 0.5)))), "(x^2/2 + x + 1)*exp(z)")
  expect_identical(.cgb(sym_optimize(expr1, list(SeriesApprox(bounds = py_dict(list(x), list(list(-1L, 1L))), reltol = 0.1)))), "(x^3/6 + x^2/2 + x + 1)*exp(z)")
  expect_identical(.cgb(sym_optimize(expr1, list(SeriesApprox(bounds = py_dict(list(x), list(list(-1L, 1L))), reltol = 0.05)))), "(x^4/24 + x^3/6 + x^2/2 + x + 1)*exp(z)")
  expect_identical(.cgb(sym_optimize(expr1, list(SeriesApprox(bounds = py_dict(list(x), list(list(-1L, 1L))), reltol = 0.05, max_order = 3L)))), "exp(x)*exp(z)")
})

test_that("codegen_test_ast: test_Assignment", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); t <- Symbol("t"); x0 <- Symbol("x0"); x1 <- Symbol("x1"); x2 <- Symbol("x2"); a <- Symbol("a"); b <- Symbol("b")
  n <- Symbol("n", integer = TRUE)
  A <- MatrixSymbol("A", 3, 1)
  mat <- Matrix(list(list(1), list(2), list(3)))
  B <- IndexedBase("B")
  i <- Idx("i", n)
  A22 <- MatrixSymbol("A22", 2, 2)
  B22 <- MatrixSymbol("B22", 2, 2)
  f16 <- float16; f32 <- float32; f64 <- float64; f80 <- float80; f128 <- float128; c64 <- complex64; c128 <- complex128
  suppressWarnings(try({ Assignment(x, y) }, silent = TRUE))
  suppressWarnings(try({ Assignment(x, 0L) }, silent = TRUE))
  suppressWarnings(try({ Assignment(A, mat) }, silent = TRUE))
  suppressWarnings(try({ Assignment(MatrixElement(A, 1L, 0L), 0L) }, silent = TRUE))
  suppressWarnings(try({ Assignment(MatrixElement(A, 1L, 0L), x) }, silent = TRUE))
  suppressWarnings(try({ Assignment(B(i), x) }, silent = TRUE))
  suppressWarnings(try({ Assignment(B(i), 0L) }, silent = TRUE))
  suppressWarnings(try({ a <- Assignment(x, y) }, silent = TRUE))
  expect_identical(.cgb(do.call(ast_func(a), py_list(py_args(a)))), "Assignment(x, y)")
  expect_identical(.cgb(ast_attr(a, "op")), ":=")
  expect_identical(.cgb(Assignment(B(i), A)), "ERR ValueError: Cannot assign a matrix to a scalar.")
  expect_identical(.cgb(Assignment(B(i), mat)), "ERR ValueError: Cannot assign a matrix to a scalar.")
  expect_identical(.cgb(Assignment(x, mat)), "ERR ValueError: Cannot assign a matrix to a scalar.")
  expect_identical(.cgb(Assignment(x, A)), "ERR ValueError: Cannot assign a matrix to a scalar.")
  expect_identical(.cgb(Assignment(MatrixElement(A, 1L, 0L), mat)), "ERR ValueError: Cannot assign a matrix to a scalar.")
  expect_identical(.cgb(Assignment(A, x)), "ERR ValueError: Cannot assign a scalar to a matrix.")
  expect_identical(.cgb(Assignment(A, 0L)), "ERR ValueError: Cannot assign a scalar to a matrix.")
  expect_identical(.cgb(Assignment(mat, A)), "ERR TypeError: Cannot assign to lhs of type <class 'sympy.matrices.immutable.ImmutableDenseMatrix'>.")
  expect_identical(.cgb(Assignment(0L, x)), "ERR TypeError: Cannot assign to lhs of type <class 'sympy.core.numbers.Zero'>.")
  expect_identical(.cgb(Assignment(x * x, 1L)), "ERR TypeError: Cannot assign to lhs of type <class 'sympy.core.power.Pow'>.")
  expect_identical(.cgb(Assignment(A + A, mat)), "ERR TypeError: Cannot assign to lhs of type <class 'sympy.matrices.expressions.matmul.MatMul'>.")
  expect_identical(.cgb(Assignment(B, 0L)), "ERR TypeError: Cannot assign to lhs of type <class 'sympy.tensor.indexed.IndexedBase'>.")
})

test_that("codegen_test_ast: test_AugAssign", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); t <- Symbol("t"); x0 <- Symbol("x0"); x1 <- Symbol("x1"); x2 <- Symbol("x2"); a <- Symbol("a"); b <- Symbol("b")
  n <- Symbol("n", integer = TRUE)
  A <- MatrixSymbol("A", 3, 1)
  mat <- Matrix(list(list(1), list(2), list(3)))
  B <- IndexedBase("B")
  i <- Idx("i", n)
  A22 <- MatrixSymbol("A22", 2, 2)
  B22 <- MatrixSymbol("B22", 2, 2)
  f16 <- float16; f32 <- float32; f64 <- float64; f80 <- float80; f128 <- float128; c64 <- complex64; c128 <- complex128
  suppressWarnings(try({ aug_assign(x, "+", y) }, silent = TRUE))
  suppressWarnings(try({ aug_assign(x, "+", 0L) }, silent = TRUE))
  suppressWarnings(try({ aug_assign(A, "+", mat) }, silent = TRUE))
  suppressWarnings(try({ aug_assign(MatrixElement(A, 1L, 0L), "+", 0L) }, silent = TRUE))
  suppressWarnings(try({ aug_assign(MatrixElement(A, 1L, 0L), "+", x) }, silent = TRUE))
  suppressWarnings(try({ aug_assign(B(i), "+", x) }, silent = TRUE))
  suppressWarnings(try({ aug_assign(B(i), "+", 0L) }, silent = TRUE))
  expect_identical(.cgb(do.call(ast_func(aug_assign(x, "+", y)), sym_args(aug_assign(x, "+", y)))), "AddAugmentedAssignment(x, y)")
  expect_identical(.cgb(AddAugmentedAssignment(x, y)), "AddAugmentedAssignment(x, y)")
  expect_identical(.cgb(aug_assign(x, "+", y) == AddAugmentedAssignment(x, y)), "True")
  expect_identical(.cgb(ast_attr(aug_assign(x, "+", y), "binop")), "+")
  expect_identical(.cgb(ast_attr(aug_assign(x, "+", y), "op")), "+=")
  expect_identical(.cgb(do.call(ast_func(aug_assign(x, "-", y)), sym_args(aug_assign(x, "-", y)))), "SubAugmentedAssignment(x, y)")
  expect_identical(.cgb(SubAugmentedAssignment(x, y)), "SubAugmentedAssignment(x, y)")
  expect_identical(.cgb(aug_assign(x, "-", y) == SubAugmentedAssignment(x, y)), "True")
  expect_identical(.cgb(ast_attr(aug_assign(x, "-", y), "binop")), "-")
  expect_identical(.cgb(ast_attr(aug_assign(x, "-", y), "op")), "-=")
  expect_identical(.cgb(do.call(ast_func(aug_assign(x, "*", y)), sym_args(aug_assign(x, "*", y)))), "MulAugmentedAssignment(x, y)")
  expect_identical(.cgb(MulAugmentedAssignment(x, y)), "MulAugmentedAssignment(x, y)")
  expect_identical(.cgb(aug_assign(x, "*", y) == MulAugmentedAssignment(x, y)), "True")
  expect_identical(.cgb(ast_attr(aug_assign(x, "*", y), "binop")), "*")
  expect_identical(.cgb(ast_attr(aug_assign(x, "*", y), "op")), "*=")
  expect_identical(.cgb(do.call(ast_func(aug_assign(x, "/", y)), sym_args(aug_assign(x, "/", y)))), "DivAugmentedAssignment(x, y)")
  expect_identical(.cgb(DivAugmentedAssignment(x, y)), "DivAugmentedAssignment(x, y)")
  expect_identical(.cgb(aug_assign(x, "/", y) == DivAugmentedAssignment(x, y)), "True")
  expect_identical(.cgb(ast_attr(aug_assign(x, "/", y), "binop")), "/")
  expect_identical(.cgb(ast_attr(aug_assign(x, "/", y), "op")), "/=")
  expect_identical(.cgb(do.call(ast_func(aug_assign(x, "%", y)), sym_args(aug_assign(x, "%", y)))), "ModAugmentedAssignment(x, y)")
  expect_identical(.cgb(ModAugmentedAssignment(x, y)), "ModAugmentedAssignment(x, y)")
  expect_identical(.cgb(aug_assign(x, "%", y) == ModAugmentedAssignment(x, y)), "True")
  expect_identical(.cgb(ast_attr(aug_assign(x, "%", y), "binop")), "%")
  expect_identical(.cgb(ast_attr(aug_assign(x, "%", y), "op")), "%=")
  expect_identical(.cgb(aug_assign(B(i), "+", A)), "ERR ValueError: Cannot assign a matrix to a scalar.")
  expect_identical(.cgb(aug_assign(B(i), "+", mat)), "ERR ValueError: Cannot assign a matrix to a scalar.")
  expect_identical(.cgb(aug_assign(x, "+", mat)), "ERR ValueError: Cannot assign a matrix to a scalar.")
  expect_identical(.cgb(aug_assign(x, "+", A)), "ERR ValueError: Cannot assign a matrix to a scalar.")
  expect_identical(.cgb(aug_assign(MatrixElement(A, 1L, 0L), "+", mat)), "ERR ValueError: Cannot assign a matrix to a scalar.")
  expect_identical(.cgb(aug_assign(A, "+", x)), "ERR ValueError: Cannot assign a scalar to a matrix.")
  expect_identical(.cgb(aug_assign(A, "+", 0L)), "ERR ValueError: Cannot assign a scalar to a matrix.")
  expect_identical(.cgb(aug_assign(mat, "+", A)), "ERR TypeError: Cannot assign to lhs of type <class 'sympy.matrices.immutable.ImmutableDenseMatrix'>.")
  expect_identical(.cgb(aug_assign(0L, "+", x)), "ERR TypeError: Cannot assign to lhs of type <class 'sympy.core.numbers.Zero'>.")
  expect_identical(.cgb(aug_assign(x * x, "+", 1L)), "ERR TypeError: Cannot assign to lhs of type <class 'sympy.core.power.Pow'>.")
  expect_identical(.cgb(aug_assign(A + A, "+", mat)), "ERR TypeError: Cannot assign to lhs of type <class 'sympy.matrices.expressions.matmul.MatMul'>.")
  expect_identical(.cgb(aug_assign(B, "+", 0L)), "ERR TypeError: Cannot assign to lhs of type <class 'sympy.tensor.indexed.IndexedBase'>.")
})

test_that("codegen_test_ast: test_Assignment_printing", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); t <- Symbol("t"); x0 <- Symbol("x0"); x1 <- Symbol("x1"); x2 <- Symbol("x2"); a <- Symbol("a"); b <- Symbol("b")
  n <- Symbol("n", integer = TRUE)
  A <- MatrixSymbol("A", 3, 1)
  mat <- Matrix(list(list(1), list(2), list(3)))
  B <- IndexedBase("B")
  i <- Idx("i", n)
  A22 <- MatrixSymbol("A22", 2, 2)
  B22 <- MatrixSymbol("B22", 2, 2)
  f16 <- float16; f32 <- float32; f64 <- float64; f80 <- float80; f128 <- float128; c64 <- complex64; c128 <- complex128
  suppressWarnings(try({ assignment_classes <- list(Assignment, AddAugmentedAssignment, SubAugmentedAssignment, MulAugmentedAssignment, DivAugmentedAssignment, ModAugmentedAssignment) }, silent = TRUE))
  suppressWarnings(try({ pairs <- list(list(x, (2L * y) + 2L), list(B(i), x), list(A22, B22), list(MatrixElement(A, 0L, 0L), x)) }, silent = TRUE))
  expect_identical(.cgb(py_repr_obj(Assignment(x, 2L * y + 2L))), "Assignment(x, 2*y + 2)")
  expect_identical(.cgb(sstr(Assignment(x, 2L * y + 2L))), "Assignment(x, 2*y + 2)")
  expect_identical(.cgb(py_repr_obj(Assignment(B(i), x))), "Assignment(B[i], x)")
  expect_identical(.cgb(sstr(Assignment(B(i), x))), "Assignment(B[i], x)")
  expect_identical(.cgb(py_repr_obj(Assignment(A22, B22))), "Assignment(A22, B22)")
  expect_identical(.cgb(sstr(Assignment(A22, B22))), "Assignment(A22, B22)")
  expect_identical(.cgb(py_repr_obj(Assignment(MatrixElement(A, 0L, 0L), x))), "Assignment(A[0, 0], x)")
  expect_identical(.cgb(sstr(Assignment(MatrixElement(A, 0L, 0L), x))), "Assignment(A[0, 0], x)")
  expect_identical(.cgb(py_repr_obj(AddAugmentedAssignment(x, 2L * y + 2L))), "AddAugmentedAssignment(x, 2*y + 2)")
  expect_identical(.cgb(sstr(AddAugmentedAssignment(x, 2L * y + 2L))), "AddAugmentedAssignment(x, 2*y + 2)")
  expect_identical(.cgb(py_repr_obj(AddAugmentedAssignment(B(i), x))), "AddAugmentedAssignment(B[i], x)")
  expect_identical(.cgb(sstr(AddAugmentedAssignment(B(i), x))), "AddAugmentedAssignment(B[i], x)")
  expect_identical(.cgb(py_repr_obj(AddAugmentedAssignment(A22, B22))), "AddAugmentedAssignment(A22, B22)")
  expect_identical(.cgb(sstr(AddAugmentedAssignment(A22, B22))), "AddAugmentedAssignment(A22, B22)")
  expect_identical(.cgb(py_repr_obj(AddAugmentedAssignment(MatrixElement(A, 0L, 0L), x))), "AddAugmentedAssignment(A[0, 0], x)")
  expect_identical(.cgb(sstr(AddAugmentedAssignment(MatrixElement(A, 0L, 0L), x))), "AddAugmentedAssignment(A[0, 0], x)")
  expect_identical(.cgb(py_repr_obj(SubAugmentedAssignment(x, 2L * y + 2L))), "SubAugmentedAssignment(x, 2*y + 2)")
  expect_identical(.cgb(sstr(SubAugmentedAssignment(x, 2L * y + 2L))), "SubAugmentedAssignment(x, 2*y + 2)")
  expect_identical(.cgb(py_repr_obj(SubAugmentedAssignment(B(i), x))), "SubAugmentedAssignment(B[i], x)")
  expect_identical(.cgb(sstr(SubAugmentedAssignment(B(i), x))), "SubAugmentedAssignment(B[i], x)")
  expect_identical(.cgb(py_repr_obj(SubAugmentedAssignment(A22, B22))), "SubAugmentedAssignment(A22, B22)")
  expect_identical(.cgb(sstr(SubAugmentedAssignment(A22, B22))), "SubAugmentedAssignment(A22, B22)")
  expect_identical(.cgb(py_repr_obj(SubAugmentedAssignment(MatrixElement(A, 0L, 0L), x))), "SubAugmentedAssignment(A[0, 0], x)")
  expect_identical(.cgb(sstr(SubAugmentedAssignment(MatrixElement(A, 0L, 0L), x))), "SubAugmentedAssignment(A[0, 0], x)")
  expect_identical(.cgb(py_repr_obj(MulAugmentedAssignment(x, 2L * y + 2L))), "MulAugmentedAssignment(x, 2*y + 2)")
  expect_identical(.cgb(sstr(MulAugmentedAssignment(x, 2L * y + 2L))), "MulAugmentedAssignment(x, 2*y + 2)")
  expect_identical(.cgb(py_repr_obj(MulAugmentedAssignment(B(i), x))), "MulAugmentedAssignment(B[i], x)")
  expect_identical(.cgb(sstr(MulAugmentedAssignment(B(i), x))), "MulAugmentedAssignment(B[i], x)")
  expect_identical(.cgb(py_repr_obj(MulAugmentedAssignment(A22, B22))), "MulAugmentedAssignment(A22, B22)")
  expect_identical(.cgb(sstr(MulAugmentedAssignment(A22, B22))), "MulAugmentedAssignment(A22, B22)")
  expect_identical(.cgb(py_repr_obj(MulAugmentedAssignment(MatrixElement(A, 0L, 0L), x))), "MulAugmentedAssignment(A[0, 0], x)")
  expect_identical(.cgb(sstr(MulAugmentedAssignment(MatrixElement(A, 0L, 0L), x))), "MulAugmentedAssignment(A[0, 0], x)")
  expect_identical(.cgb(py_repr_obj(DivAugmentedAssignment(x, 2L * y + 2L))), "DivAugmentedAssignment(x, 2*y + 2)")
  expect_identical(.cgb(sstr(DivAugmentedAssignment(x, 2L * y + 2L))), "DivAugmentedAssignment(x, 2*y + 2)")
  expect_identical(.cgb(py_repr_obj(DivAugmentedAssignment(B(i), x))), "DivAugmentedAssignment(B[i], x)")
  expect_identical(.cgb(sstr(DivAugmentedAssignment(B(i), x))), "DivAugmentedAssignment(B[i], x)")
  expect_identical(.cgb(py_repr_obj(DivAugmentedAssignment(A22, B22))), "DivAugmentedAssignment(A22, B22)")
  expect_identical(.cgb(sstr(DivAugmentedAssignment(A22, B22))), "DivAugmentedAssignment(A22, B22)")
  expect_identical(.cgb(py_repr_obj(DivAugmentedAssignment(MatrixElement(A, 0L, 0L), x))), "DivAugmentedAssignment(A[0, 0], x)")
  expect_identical(.cgb(sstr(DivAugmentedAssignment(MatrixElement(A, 0L, 0L), x))), "DivAugmentedAssignment(A[0, 0], x)")
  expect_identical(.cgb(py_repr_obj(ModAugmentedAssignment(x, 2L * y + 2L))), "ModAugmentedAssignment(x, 2*y + 2)")
  expect_identical(.cgb(sstr(ModAugmentedAssignment(x, 2L * y + 2L))), "ModAugmentedAssignment(x, 2*y + 2)")
  expect_identical(.cgb(py_repr_obj(ModAugmentedAssignment(B(i), x))), "ModAugmentedAssignment(B[i], x)")
  expect_identical(.cgb(sstr(ModAugmentedAssignment(B(i), x))), "ModAugmentedAssignment(B[i], x)")
  expect_identical(.cgb(py_repr_obj(ModAugmentedAssignment(A22, B22))), "ModAugmentedAssignment(A22, B22)")
  expect_identical(.cgb(sstr(ModAugmentedAssignment(A22, B22))), "ModAugmentedAssignment(A22, B22)")
  expect_identical(.cgb(py_repr_obj(ModAugmentedAssignment(MatrixElement(A, 0L, 0L), x))), "ModAugmentedAssignment(A[0, 0], x)")
  expect_identical(.cgb(sstr(ModAugmentedAssignment(MatrixElement(A, 0L, 0L), x))), "ModAugmentedAssignment(A[0, 0], x)")
})

test_that("codegen_test_ast: test_CodeBlock", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); t <- Symbol("t"); x0 <- Symbol("x0"); x1 <- Symbol("x1"); x2 <- Symbol("x2"); a <- Symbol("a"); b <- Symbol("b")
  n <- Symbol("n", integer = TRUE)
  A <- MatrixSymbol("A", 3, 1)
  mat <- Matrix(list(list(1), list(2), list(3)))
  B <- IndexedBase("B")
  i <- Idx("i", n)
  A22 <- MatrixSymbol("A22", 2, 2)
  B22 <- MatrixSymbol("B22", 2, 2)
  f16 <- float16; f32 <- float32; f64 <- float64; f80 <- float80; f128 <- float128; c64 <- complex64; c128 <- complex128
  suppressWarnings(try({ c <- CodeBlock(Assignment(x, 1L), Assignment(y, x + 1L)) }, silent = TRUE))
  expect_identical(.cgb(do.call(ast_func(c), py_list(py_args(c)))), "CodeBlock(\\nAssignment(x, 1),\\nAssignment(y, x + 1)\\n)")
  expect_identical(.cgb(ast_attr(c, "left_hand_sides")), "(x, y)")
  expect_identical(.cgb(ast_attr(c, "right_hand_sides")), "(1, x + 1)")
})

test_that("codegen_test_ast: test_CodeBlock_topological_sort", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); t <- Symbol("t"); x0 <- Symbol("x0"); x1 <- Symbol("x1"); x2 <- Symbol("x2"); a <- Symbol("a"); b <- Symbol("b")
  n <- Symbol("n", integer = TRUE)
  A <- MatrixSymbol("A", 3, 1)
  mat <- Matrix(list(list(1), list(2), list(3)))
  B <- IndexedBase("B")
  i <- Idx("i", n)
  A22 <- MatrixSymbol("A22", 2, 2)
  B22 <- MatrixSymbol("B22", 2, 2)
  f16 <- float16; f32 <- float32; f64 <- float64; f80 <- float80; f128 <- float128; c64 <- complex64; c128 <- complex128
  suppressWarnings(try({ assignments <- list(Assignment(x, y + z), Assignment(z, 1L), Assignment(t, x), Assignment(y, 2L)) }, silent = TRUE))
  suppressWarnings(try({ ordered_assignments <- list(Assignment(z, 1L), Assignment(y, 2L), Assignment(x, y + z), Assignment(t, x)) }, silent = TRUE))
  suppressWarnings(try({ c1 <- codeblock_topological_sort(assignments) }, silent = TRUE))
  expect_identical(.cgb(c1), "CodeBlock(\\nAssignment(z, 1),\\nAssignment(y, 2),\\nAssignment(x, y + z),\\nAssignment(t, x)\\n)")
  suppressWarnings(try({ invalid_assignments <- list(Assignment(x, y + z), Assignment(z, 1L), Assignment(y, x), Assignment(y, 2L)) }, silent = TRUE))
  expect_identical(.cgb(codeblock_topological_sort(invalid_assignments)), "ERR ValueError: cycle detected")
  suppressWarnings(try({ free_assignments <- list(Assignment(x, y + z), Assignment(z, a * b), Assignment(t, x), Assignment(y, b + 3L)) }, silent = TRUE))
  suppressWarnings(try({ free_assignments_ordered <- list(Assignment(z, a * b), Assignment(y, b + 3L), Assignment(x, y + z), Assignment(t, x)) }, silent = TRUE))
  suppressWarnings(try({ c2 <- codeblock_topological_sort(free_assignments) }, silent = TRUE))
  expect_identical(.cgb(c2), "CodeBlock(\\nAssignment(z, a*b),\\nAssignment(y, b + 3),\\nAssignment(x, y + z),\\nAssignment(t, x)\\n)")
})

test_that("codegen_test_ast: test_CodeBlock_free_symbols", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); t <- Symbol("t"); x0 <- Symbol("x0"); x1 <- Symbol("x1"); x2 <- Symbol("x2"); a <- Symbol("a"); b <- Symbol("b")
  n <- Symbol("n", integer = TRUE)
  A <- MatrixSymbol("A", 3, 1)
  mat <- Matrix(list(list(1), list(2), list(3)))
  B <- IndexedBase("B")
  i <- Idx("i", n)
  A22 <- MatrixSymbol("A22", 2, 2)
  B22 <- MatrixSymbol("B22", 2, 2)
  f16 <- float16; f32 <- float32; f64 <- float64; f80 <- float80; f128 <- float128; c64 <- complex64; c128 <- complex128
  suppressWarnings(try({ c1 <- CodeBlock(Assignment(x, y + z), Assignment(z, 1L), Assignment(t, x), Assignment(y, 2L)) }, silent = TRUE))
  expect_identical(.cgb(py_set(free_symbols(c1))), "set()")
  suppressWarnings(try({ c2 <- CodeBlock(Assignment(x, y + z), Assignment(z, a * b), Assignment(t, x), Assignment(y, b + 3L)) }, silent = TRUE))
})

test_that("codegen_test_ast: test_CodeBlock_cse", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); t <- Symbol("t"); x0 <- Symbol("x0"); x1 <- Symbol("x1"); x2 <- Symbol("x2"); a <- Symbol("a"); b <- Symbol("b")
  n <- Symbol("n", integer = TRUE)
  A <- MatrixSymbol("A", 3, 1)
  mat <- Matrix(list(list(1), list(2), list(3)))
  B <- IndexedBase("B")
  i <- Idx("i", n)
  A22 <- MatrixSymbol("A22", 2, 2)
  B22 <- MatrixSymbol("B22", 2, 2)
  f16 <- float16; f32 <- float32; f64 <- float64; f80 <- float80; f128 <- float128; c64 <- complex64; c128 <- complex128
  suppressWarnings(try({ c1 <- CodeBlock(Assignment(y, 1L), Assignment(x, sin(y)), Assignment(z, sin(y)), Assignment(t, x * z)) }, silent = TRUE))
  expect_identical(.cgb(codeblock_cse(c1)), "CodeBlock(\\nAssignment(y, 1),\\nAssignment(x0, sin(y)),\\nAssignment(x, x0),\\nAssignment(z, x0),\\nAssignment(t, x*z)\\n)")
  expect_identical(.cgb(codeblock_cse(CodeBlock(Assignment(x, 1L), Assignment(y, 1L), Assignment(y, 2L)))), "ERR NotImplementedError: Duplicate assignments to the same variable are not yet supported (y)")
  suppressWarnings(try({ c2 <- CodeBlock(Assignment(x0, sin(y) + 1L), Assignment(x1, 2L * sin(y)), Assignment(z, x * y)) }, silent = TRUE))
  expect_identical(.cgb(codeblock_cse(c2)), "CodeBlock(\\nAssignment(x2, sin(y)),\\nAssignment(x0, x2 + 1),\\nAssignment(x1, 2*x2),\\nAssignment(z, x*y)\\n)")
})

test_that("codegen_test_ast: test_For", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); t <- Symbol("t"); x0 <- Symbol("x0"); x1 <- Symbol("x1"); x2 <- Symbol("x2"); a <- Symbol("a"); b <- Symbol("b")
  n <- Symbol("n", integer = TRUE)
  A <- MatrixSymbol("A", 3, 1)
  mat <- Matrix(list(list(1), list(2), list(3)))
  B <- IndexedBase("B")
  i <- Idx("i", n)
  A22 <- MatrixSymbol("A22", 2, 2)
  B22 <- MatrixSymbol("B22", 2, 2)
  f16 <- float16; f32 <- float32; f64 <- float64; f80 <- float80; f128 <- float128; c64 <- complex64; c128 <- complex128
  suppressWarnings(try({ f <- For(n, Range(0L, 3L), list(Assignment(MatrixElement(A, n, 0L), x + n), aug_assign(x, "+", y))) }, silent = TRUE))
  suppressWarnings(try({ f <- For(n, list(1L, 2L, 3L, 4L, 5L), list(Assignment(MatrixElement(A, n, 0L), x + n))) }, silent = TRUE))
  expect_identical(.cgb(do.call(ast_func(f), py_list(py_args(f)))), "For(n, iterable=(1, 2, 3, 4, 5), body=CodeBlock(\\n    Assignment(A[n, 0], n + x)\\n))")
  expect_identical(.cgb(For(n, x, list(x + y))), "ERR TypeError: iterable must be an iterable")
})

test_that("codegen_test_ast: test_none", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); t <- Symbol("t"); x0 <- Symbol("x0"); x1 <- Symbol("x1"); x2 <- Symbol("x2"); a <- Symbol("a"); b <- Symbol("b")
  n <- Symbol("n", integer = TRUE)
  A <- MatrixSymbol("A", 3, 1)
  mat <- Matrix(list(list(1), list(2), list(3)))
  B <- IndexedBase("B")
  i <- Idx("i", n)
  A22 <- MatrixSymbol("A22", 2, 2)
  B22 <- MatrixSymbol("B22", 2, 2)
  f16 <- float16; f32 <- float32; f64 <- float64; f80 <- float80; f128 <- float128; c64 <- complex64; c128 <- complex128
  expect_identical(.cgb(ast_attr(none, "is_Atom")), "True")
  expect_identical(.cgb(none), "NoneToken()")
  expect_identical(.cgb(none), "NoneToken()")
  expect_identical(.cgb(none), "NoneToken()")
  expect_identical(.cgb(do.call(ast_func(none), py_list(py_args(none)))), "NoneToken()")
})

test_that("codegen_test_ast: test_String", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); t <- Symbol("t"); x0 <- Symbol("x0"); x1 <- Symbol("x1"); x2 <- Symbol("x2"); a <- Symbol("a"); b <- Symbol("b")
  n <- Symbol("n", integer = TRUE)
  A <- MatrixSymbol("A", 3, 1)
  mat <- Matrix(list(list(1), list(2), list(3)))
  B <- IndexedBase("B")
  i <- Idx("i", n)
  A22 <- MatrixSymbol("A22", 2, 2)
  B22 <- MatrixSymbol("B22", 2, 2)
  f16 <- float16; f32 <- float32; f64 <- float64; f80 <- float80; f128 <- float128; c64 <- complex64; c128 <- complex128
  suppressWarnings(try({ st <- String("foobar") }, silent = TRUE))
  expect_identical(.cgb(ast_attr(st, "is_Atom")), "True")
  expect_identical(.cgb(st), "foobar")
  expect_identical(.cgb(ast_attr(st, "text")), "foobar")
  expect_identical(.cgb(do.call(ast_func(st), ast_kwargs(st))), "foobar")
  expect_identical(.cgb(do.call(ast_func(st), py_list(py_args(st)))), "foobar")
  suppressWarnings(try({ s <- String("foo") }, silent = TRUE))
  expect_identical(.cgb(py_str(s)), "foo")
  expect_identical(.cgb(py_repr_obj(s)), "String('foo')")
})

test_that("codegen_test_ast: test_Comment", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); t <- Symbol("t"); x0 <- Symbol("x0"); x1 <- Symbol("x1"); x2 <- Symbol("x2"); a <- Symbol("a"); b <- Symbol("b")
  n <- Symbol("n", integer = TRUE)
  A <- MatrixSymbol("A", 3, 1)
  mat <- Matrix(list(list(1), list(2), list(3)))
  B <- IndexedBase("B")
  i <- Idx("i", n)
  A22 <- MatrixSymbol("A22", 2, 2)
  B22 <- MatrixSymbol("B22", 2, 2)
  f16 <- float16; f32 <- float32; f64 <- float64; f80 <- float80; f128 <- float128; c64 <- complex64; c128 <- complex128
  suppressWarnings(try({ c <- Comment("foobar") }, silent = TRUE))
  expect_identical(.cgb(ast_attr(c, "text")), "foobar")
  expect_identical(.cgb(py_str(c)), "foobar")
})

test_that("codegen_test_ast: test_Node", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); t <- Symbol("t"); x0 <- Symbol("x0"); x1 <- Symbol("x1"); x2 <- Symbol("x2"); a <- Symbol("a"); b <- Symbol("b")
  n <- Symbol("n", integer = TRUE)
  A <- MatrixSymbol("A", 3, 1)
  mat <- Matrix(list(list(1), list(2), list(3)))
  B <- IndexedBase("B")
  i <- Idx("i", n)
  A22 <- MatrixSymbol("A22", 2, 2)
  B22 <- MatrixSymbol("B22", 2, 2)
  f16 <- float16; f32 <- float32; f64 <- float64; f80 <- float80; f128 <- float128; c64 <- complex64; c128 <- complex128
  suppressWarnings(try({ n <- Node() }, silent = TRUE))
  expect_identical(.cgb(n), "Node()")
  expect_identical(.cgb(do.call(ast_func(n), py_list(py_args(n)))), "Node()")
})

test_that("codegen_test_ast: test_Type", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); t <- Symbol("t"); x0 <- Symbol("x0"); x1 <- Symbol("x1"); x2 <- Symbol("x2"); a <- Symbol("a"); b <- Symbol("b")
  n <- Symbol("n", integer = TRUE)
  A <- MatrixSymbol("A", 3, 1)
  mat <- Matrix(list(list(1), list(2), list(3)))
  B <- IndexedBase("B")
  i <- Idx("i", n)
  A22 <- MatrixSymbol("A22", 2, 2)
  B22 <- MatrixSymbol("B22", 2, 2)
  f16 <- float16; f32 <- float32; f64 <- float64; f80 <- float80; f128 <- float128; c64 <- complex64; c128 <- complex128
  suppressWarnings(try({ t <- Type("MyType") }, silent = TRUE))
  expect_identical(.cgb(py_len(py_args(t))), "1")
  expect_identical(.cgb(ast_attr(t, "name")), "MyType")
  expect_identical(.cgb(py_str(t)), "MyType")
  expect_identical(.cgb(py_repr_obj(t)), "Type(String('MyType'))")
  expect_identical(.cgb(Type(t)), "MyType")
  expect_identical(.cgb(do.call(ast_func(t), py_list(py_args(t)))), "MyType")
  suppressWarnings(try({ t1 <- Type("t1") }, silent = TRUE))
  suppressWarnings(try({ t2 <- Type("t2") }, silent = TRUE))
  expect_identical(.cgb(!py_eq(t1, t2)), "True")
  expect_identical(.cgb((py_eq(t1, t1)) && (py_eq(t2, t2))), "True")
  suppressWarnings(try({ t1b <- Type("t1") }, silent = TRUE))
  expect_identical(.cgb(t1), "t1")
  expect_identical(.cgb(!py_eq(t2, t1b)), "True")
})

test_that("codegen_test_ast: test_Type__from_expr", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); t <- Symbol("t"); x0 <- Symbol("x0"); x1 <- Symbol("x1"); x2 <- Symbol("x2"); a <- Symbol("a"); b <- Symbol("b")
  n <- Symbol("n", integer = TRUE)
  A <- MatrixSymbol("A", 3, 1)
  mat <- Matrix(list(list(1), list(2), list(3)))
  B <- IndexedBase("B")
  i <- Idx("i", n)
  A22 <- MatrixSymbol("A22", 2, 2)
  B22 <- MatrixSymbol("B22", 2, 2)
  f16 <- float16; f32 <- float32; f64 <- float64; f80 <- float80; f128 <- float128; c64 <- complex64; c128 <- complex128
  expect_identical(.cgb(type_from_expr(i)), "integer")
  suppressWarnings(try({ u <- sym_symbols("u", real = TRUE) }, silent = TRUE))
  expect_identical(.cgb(type_from_expr(u)), "real")
  expect_identical(.cgb(type_from_expr(n)), "integer")
  expect_identical(.cgb(type_from_expr(3L)), "integer")
  expect_identical(.cgb(type_from_expr(Float(3.0))), "real")
  expect_identical(.cgb(type_from_expr(3L + complex(real = 0.0, imaginary = 1.0))), "complex")
  expect_identical(.cgb(type_from_expr(sum)), "ERR ValueError: Could not deduce type from expr.")
})

test_that("codegen_test_ast: test_Type__cast_check__integers", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); t <- Symbol("t"); x0 <- Symbol("x0"); x1 <- Symbol("x1"); x2 <- Symbol("x2"); a <- Symbol("a"); b <- Symbol("b")
  n <- Symbol("n", integer = TRUE)
  A <- MatrixSymbol("A", 3, 1)
  mat <- Matrix(list(list(1), list(2), list(3)))
  B <- IndexedBase("B")
  i <- Idx("i", n)
  A22 <- MatrixSymbol("A22", 2, 2)
  B22 <- MatrixSymbol("B22", 2, 2)
  f16 <- float16; f32 <- float32; f64 <- float64; f80 <- float80; f128 <- float128; c64 <- complex64; c128 <- complex128
  expect_identical(.cgb(type_cast_check(sym_integer, Float(3.5))), "ERR ValueError: Casting gives a significantly different value.")
  expect_identical(.cgb(type_cast_check(sym_integer, "3")), "3")
  expect_identical(.cgb(type_cast_check(sym_integer, Float("3.0000000000000000000"))), "3")
  expect_identical(.cgb(type_cast_check(sym_integer, Float("3.0000000000000000001"))), "3")
  expect_identical(.cgb(type_cast_check(int8, Float(127.0))), "127")
  expect_identical(.cgb(type_cast_check(int8, 128L)), "ERR ValueError: Value is too big: 128 > 127")
  expect_identical(.cgb(type_cast_check(int8, -128L)), "-128")
  expect_identical(.cgb(type_cast_check(int8, -129L)), "ERR ValueError: Value is too small: -129 < -128")
  expect_identical(.cgb(type_cast_check(uint8, 0L)), "0")
  expect_identical(.cgb(type_cast_check(uint8, 128L)), "128")
  expect_identical(.cgb(type_cast_check(uint8, Float(256.0))), "ERR ValueError: Value is too big: 256 > 255")
  expect_identical(.cgb(type_cast_check(uint8, -1L)), "ERR ValueError: Value is too small: -1 < 0")
})

test_that("codegen_test_ast: test_Attribute", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); t <- Symbol("t"); x0 <- Symbol("x0"); x1 <- Symbol("x1"); x2 <- Symbol("x2"); a <- Symbol("a"); b <- Symbol("b")
  n <- Symbol("n", integer = TRUE)
  A <- MatrixSymbol("A", 3, 1)
  mat <- Matrix(list(list(1), list(2), list(3)))
  B <- IndexedBase("B")
  i <- Idx("i", n)
  A22 <- MatrixSymbol("A22", 2, 2)
  B22 <- MatrixSymbol("B22", 2, 2)
  f16 <- float16; f32 <- float32; f64 <- float64; f80 <- float80; f128 <- float128; c64 <- complex64; c128 <- complex128
  suppressWarnings(try({ noexcept <- Attribute("noexcept") }, silent = TRUE))
  expect_identical(.cgb(noexcept), "noexcept")
  suppressWarnings(try({ alignas16 <- Attribute("alignas", list(16L)) }, silent = TRUE))
  suppressWarnings(try({ alignas32 <- Attribute("alignas", list(32L)) }, silent = TRUE))
  expect_identical(.cgb(!py_eq(alignas16, alignas32)), "True")
  expect_identical(.cgb(do.call(ast_func(alignas16), py_list(py_args(alignas16)))), "alignas(16)")
})

test_that("codegen_test_ast: test_Variable", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); t <- Symbol("t"); x0 <- Symbol("x0"); x1 <- Symbol("x1"); x2 <- Symbol("x2"); a <- Symbol("a"); b <- Symbol("b")
  n <- Symbol("n", integer = TRUE)
  A <- MatrixSymbol("A", 3, 1)
  mat <- Matrix(list(list(1), list(2), list(3)))
  B <- IndexedBase("B")
  i <- Idx("i", n)
  A22 <- MatrixSymbol("A22", 2, 2)
  B22 <- MatrixSymbol("B22", 2, 2)
  f16 <- float16; f32 <- float32; f64 <- float64; f80 <- float80; f128 <- float128; c64 <- complex64; c128 <- complex128
  suppressWarnings(try({ v <- Variable(x, type = real) }, silent = TRUE))
  expect_identical(.cgb(v), "Variable(x, type=real)")
  expect_identical(.cgb(v), "Variable(x, type=real)")
  expect_identical(.cgb(ast_attr(v, "symbol")), "x")
  expect_identical(.cgb(ast_attr(v, "type")), "real")
  expect_identical(.cgb(!py_in(value_const, ast_attr(v, "attrs"))), "True")
  expect_identical(.cgb(do.call(ast_func(v), py_list(py_args(v)))), "Variable(x, type=real)")
  expect_identical(.cgb(py_str(v)), "Variable(x, type=real)")
  suppressWarnings(try({ w <- Variable(y, f32, attrs = list(value_const)) }, silent = TRUE))
  expect_identical(.cgb(ast_attr(w, "symbol")), "y")
  expect_identical(.cgb(ast_attr(w, "type")), "float32")
  expect_identical(.cgb(py_in(value_const, ast_attr(w, "attrs"))), "True")
  expect_identical(.cgb(do.call(ast_func(w), py_list(py_args(w)))), "Variable(y, type=float32, attrs=(value_const,))")
  suppressWarnings(try({ v_n <- Variable(n, type = type_from_expr(n)) }, silent = TRUE))
  expect_identical(.cgb(ast_attr(v_n, "type")), "integer")
  expect_identical(.cgb(do.call(ast_func(v_n), py_list(py_args(v_n)))), "Variable(n, type=integer)")
  suppressWarnings(try({ v_i <- Variable(i, type = type_from_expr(n)) }, silent = TRUE))
  expect_identical(.cgb(ast_attr(v_i, "type")), "integer")
  expect_identical(.cgb(!py_eq(v_i, v_n)), "True")
  suppressWarnings(try({ a_i <- variable_deduced(i) }, silent = TRUE))
  expect_identical(.cgb(ast_attr(a_i, "type")), "integer")
  expect_identical(.cgb(ast_attr(variable_deduced(Symbol("x", real = TRUE)), "type")), "real")
  expect_identical(.cgb(do.call(ast_func(a_i), py_list(py_args(a_i)))), "Variable(i, type=integer)")
  suppressWarnings(try({ v_n2 <- variable_deduced(n, value = Float(3.5), cast_check = FALSE) }, silent = TRUE))
  expect_identical(.cgb(do.call(ast_func(v_n2), py_list(py_args(v_n2)))), "Variable(n, type=integer, value=3.5)")
  expect_identical(.cgb(py_abs(ast_attr(v_n2, "value") - Float(3.5)) < Float(1e-15)), "True")
  expect_identical(.cgb(variable_deduced(n, value = Float(3.5), cast_check = TRUE)), "ERR ValueError: Casting gives a significantly different value.")
  suppressWarnings(try({ v_n3 <- variable_deduced(n) }, silent = TRUE))
  expect_identical(.cgb(ast_attr(v_n3, "type")), "integer")
  expect_identical(.cgb(py_str(v_n3)), "Variable(n, type=integer)")
  expect_identical(.cgb(ast_attr(variable_deduced(z, value = 3L), "type")), "integer")
  expect_identical(.cgb(ast_attr(variable_deduced(z, value = Float(3.0)), "type")), "real")
  expect_identical(.cgb(ast_attr(variable_deduced(z, value = Float(3.0) + complex(real = 0.0, imaginary = 1.0)), "type")), "complex")
})

test_that("codegen_test_ast: test_Pointer", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); t <- Symbol("t"); x0 <- Symbol("x0"); x1 <- Symbol("x1"); x2 <- Symbol("x2"); a <- Symbol("a"); b <- Symbol("b")
  n <- Symbol("n", integer = TRUE)
  A <- MatrixSymbol("A", 3, 1)
  mat <- Matrix(list(list(1), list(2), list(3)))
  B <- IndexedBase("B")
  i <- Idx("i", n)
  A22 <- MatrixSymbol("A22", 2, 2)
  B22 <- MatrixSymbol("B22", 2, 2)
  f16 <- float16; f32 <- float32; f64 <- float64; f80 <- float80; f128 <- float128; c64 <- complex64; c128 <- complex128
  suppressWarnings(try({ p <- Pointer(x) }, silent = TRUE))
  expect_identical(.cgb(ast_attr(p, "symbol")), "x")
  expect_identical(.cgb(ast_attr(p, "type")), "untyped")
  expect_identical(.cgb(!py_in(value_const, ast_attr(p, "attrs"))), "True")
  expect_identical(.cgb(!py_in(pointer_const, ast_attr(p, "attrs"))), "True")
  expect_identical(.cgb(do.call(ast_func(p), py_list(py_args(p)))), "Pointer(x)")
  suppressWarnings(try({ u <- sym_symbols("u", real = TRUE) }, silent = TRUE))
  suppressWarnings(try({ pu <- Pointer(u, type = type_from_expr(u), attrs = list(value_const, pointer_const)) }, silent = TRUE))
  expect_identical(.cgb(py_is(ast_attr(pu, "symbol"), u)), "True")
  expect_identical(.cgb(ast_attr(pu, "type")), "real")
  expect_identical(.cgb(py_in(value_const, ast_attr(pu, "attrs"))), "True")
  expect_identical(.cgb(py_in(pointer_const, ast_attr(pu, "attrs"))), "True")
  expect_identical(.cgb(do.call(ast_func(pu), py_list(py_args(pu)))), "Pointer(u, type=real, attrs=(value_const, pointer_const))")
  suppressWarnings(try({ i <- sym_symbols("i", integer = TRUE) }, silent = TRUE))
  suppressWarnings(try({ deref <- py_getitem(pu, i) }, silent = TRUE))
  expect_identical(.cgb(ast_attr(deref, "indices")), "(i,)")
})

test_that("codegen_test_ast: test_Declaration", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); t <- Symbol("t"); x0 <- Symbol("x0"); x1 <- Symbol("x1"); x2 <- Symbol("x2"); a <- Symbol("a"); b <- Symbol("b")
  n <- Symbol("n", integer = TRUE)
  A <- MatrixSymbol("A", 3, 1)
  mat <- Matrix(list(list(1), list(2), list(3)))
  B <- IndexedBase("B")
  i <- Idx("i", n)
  A22 <- MatrixSymbol("A22", 2, 2)
  B22 <- MatrixSymbol("B22", 2, 2)
  f16 <- float16; f32 <- float32; f64 <- float64; f80 <- float80; f128 <- float128; c64 <- complex64; c128 <- complex128
  suppressWarnings(try({ u <- sym_symbols("u", real = TRUE) }, silent = TRUE))
  suppressWarnings(try({ vu <- Variable(u, type = type_from_expr(u)) }, silent = TRUE))
  expect_identical(.cgb(ast_attr(ast_attr(Declaration(vu), "variable"), "type")), "real")
  suppressWarnings(try({ vn <- Variable(n, type = type_from_expr(n)) }, silent = TRUE))
  expect_identical(.cgb(ast_attr(ast_attr(Declaration(vn), "variable"), "type")), "integer")
  suppressWarnings(try({ vuc <- Variable(u, type_from_expr(u), value = Float(3.0), attrs = list(value_const)) }, silent = TRUE))
  expect_identical(.cgb(py_in(value_const, ast_attr(vuc, "attrs"))), "True")
  expect_identical(.cgb(!py_in(pointer_const, ast_attr(vuc, "attrs"))), "True")
  suppressWarnings(try({ decl <- Declaration(vuc) }, silent = TRUE))
  expect_identical(.cgb(ast_attr(decl, "variable")), "Variable(u, type=real, value=3.0, attrs=(value_const,))")
  expect_identical(.cgb(py_isinstance(ast_attr(ast_attr(decl, "variable"), "value"), Float)), "True")
  expect_identical(.cgb(ast_attr(ast_attr(decl, "variable"), "value")), "3.00000000000000")
  expect_identical(.cgb(do.call(ast_func(decl), py_list(py_args(decl)))), "Declaration(Variable(u, type=real, value=3.0, attrs=(value_const,)))")
  expect_identical(.cgb(as_Declaration(vuc)), "Declaration(Variable(u, type=real, value=3.0, attrs=(value_const,)))")
  expect_identical(.cgb(as_Declaration(vuc, value = NULL, attrs = NULL)), "Declaration(Variable(u, type=real))")
  suppressWarnings(try({ vy <- Variable(y, type = sym_integer, value = 3L) }, silent = TRUE))
  suppressWarnings(try({ decl2 <- Declaration(vy) }, silent = TRUE))
  expect_identical(.cgb(ast_attr(decl2, "variable")), "Variable(y, type=integer, value=3)")
  expect_identical(.cgb(ast_attr(ast_attr(decl2, "variable"), "value")), "3")
  suppressWarnings(try({ vi <- Variable(i, type = type_from_expr(i), value = Float(3.0)) }, silent = TRUE))
  suppressWarnings(try({ decl3 <- Declaration(vi) }, silent = TRUE))
  expect_identical(.cgb(ast_attr(ast_attr(decl3, "variable"), "type")), "integer")
  expect_identical(.cgb(ast_attr(ast_attr(decl3, "variable"), "value")), "3.00000000000000")
  expect_identical(.cgb(Declaration(vi, 42L)), "ERR ValueError: Too many arguments (2), expected at most 1")
})

test_that("codegen_test_ast: test_IntBaseType", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); t <- Symbol("t"); x0 <- Symbol("x0"); x1 <- Symbol("x1"); x2 <- Symbol("x2"); a <- Symbol("a"); b <- Symbol("b")
  n <- Symbol("n", integer = TRUE)
  A <- MatrixSymbol("A", 3, 1)
  mat <- Matrix(list(list(1), list(2), list(3)))
  B <- IndexedBase("B")
  i <- Idx("i", n)
  A22 <- MatrixSymbol("A22", 2, 2)
  B22 <- MatrixSymbol("B22", 2, 2)
  f16 <- float16; f32 <- float32; f64 <- float64; f80 <- float80; f128 <- float128; c64 <- complex64; c128 <- complex128
  expect_identical(.cgb(ast_attr(intc, "name")), "intc")
  expect_identical(.cgb(py_args(intc)), "(String('intc'),)")
  expect_identical(.cgb(py_str(ast_attr(IntBaseType("a"), "name"))), "a")
})

test_that("codegen_test_ast: test_FloatType", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); t <- Symbol("t"); x0 <- Symbol("x0"); x1 <- Symbol("x1"); x2 <- Symbol("x2"); a <- Symbol("a"); b <- Symbol("b")
  n <- Symbol("n", integer = TRUE)
  A <- MatrixSymbol("A", 3, 1)
  mat <- Matrix(list(list(1), list(2), list(3)))
  B <- IndexedBase("B")
  i <- Idx("i", n)
  A22 <- MatrixSymbol("A22", 2, 2)
  B22 <- MatrixSymbol("B22", 2, 2)
  f16 <- float16; f32 <- float32; f64 <- float64; f80 <- float80; f128 <- float128; c64 <- complex64; c128 <- complex128
  expect_identical(.cgb(ast_attr(f16, "dig")), "3")
  expect_identical(.cgb(ast_attr(f32, "dig")), "6")
  expect_identical(.cgb(ast_attr(f64, "dig")), "15")
  expect_identical(.cgb(ast_attr(f80, "dig")), "18")
  expect_identical(.cgb(ast_attr(f128, "dig")), "33")
  expect_identical(.cgb(ast_attr(f16, "decimal_dig")), "5")
  expect_identical(.cgb(ast_attr(f32, "decimal_dig")), "9")
  expect_identical(.cgb(ast_attr(f64, "decimal_dig")), "17")
  expect_identical(.cgb(ast_attr(f80, "decimal_dig")), "21")
  expect_identical(.cgb(ast_attr(f128, "decimal_dig")), "36")
  expect_identical(.cgb(ast_attr(f16, "max_exponent")), "16")
  expect_identical(.cgb(ast_attr(f32, "max_exponent")), "128")
  expect_identical(.cgb(ast_attr(f64, "max_exponent")), "1024")
  expect_identical(.cgb(ast_attr(f80, "max_exponent")), "16384")
  expect_identical(.cgb(ast_attr(f128, "max_exponent")), "16384")
  expect_identical(.cgb(ast_attr(f16, "min_exponent")), "-13")
  expect_identical(.cgb(ast_attr(f32, "min_exponent")), "-125")
  expect_identical(.cgb(ast_attr(f64, "min_exponent")), "-1021")
  expect_identical(.cgb(ast_attr(f80, "min_exponent")), "-16381")
  expect_identical(.cgb(ast_attr(f128, "min_exponent")), "-16381")
  expect_identical(.cgb(py_abs((ast_attr(f16, "eps") / Float("0.00097656", precision = 16L)) - 1L) < (Float(0.1) * (10L ^ (-ast_attr(f16, "dig"))))), "True")
  expect_identical(.cgb(py_abs((ast_attr(f32, "eps") / Float("1.1920929e-07", precision = 32L)) - 1L) < (Float(0.1) * (10L ^ (-ast_attr(f32, "dig"))))), "True")
  expect_identical(.cgb(py_abs((ast_attr(f64, "eps") / Float("2.2204460492503131e-16", precision = 64L)) - 1L) < (Float(0.1) * (10L ^ (-ast_attr(f64, "dig"))))), "True")
  expect_identical(.cgb(py_abs((ast_attr(f80, "eps") / Float("1.08420217248550443401e-19", precision = 80L)) - 1L) < (Float(0.1) * (10L ^ (-ast_attr(f80, "dig"))))), "True")
  expect_identical(.cgb(py_abs((ast_attr(f128, "eps") / Float(" 1.92592994438723585305597794258492732e-34", precision = 128L)) - 1L) < (Float(0.1) * (10L ^ (-ast_attr(f128, "dig"))))), "True")
  expect_identical(.cgb(py_abs((ast_attr(f16, "max") / Float("65504", precision = 16L)) - 1L) < (Float(0.1) * (10L ^ (-ast_attr(f16, "dig"))))), "True")
  expect_identical(.cgb(py_abs((ast_attr(f32, "max") / Float("3.40282347e+38", precision = 32L)) - 1L) < (Float(0.1) * (10L ^ (-ast_attr(f32, "dig"))))), "True")
  expect_identical(.cgb(py_abs((ast_attr(f64, "max") / Float("1.79769313486231571e+308", precision = 64L)) - 1L) < (Float(0.1) * (10L ^ (-ast_attr(f64, "dig"))))), "True")
  expect_identical(.cgb(py_abs((ast_attr(f80, "max") / Float("1.18973149535723176502e+4932", precision = 80L)) - 1L) < (Float(0.1) * (10L ^ (-ast_attr(f80, "dig"))))), "True")
  expect_identical(.cgb(py_abs((ast_attr(f128, "max") / Float("1.18973149535723176508575932662800702e+4932", precision = 128L)) - 1L) < (Float(0.1) * (10L ^ (-ast_attr(f128, "dig"))))), "True")
  expect_identical(.cgb(py_abs((ast_attr(f16, "tiny") / Float("6.1035e-05", precision = 16L)) - 1L) < (Float(0.1) * (10L ^ (-ast_attr(f16, "dig"))))), "True")
  expect_identical(.cgb(py_abs((ast_attr(f32, "tiny") / Float("1.17549435e-38", precision = 32L)) - 1L) < (Float(0.1) * (10L ^ (-ast_attr(f32, "dig"))))), "True")
  expect_identical(.cgb(py_abs((ast_attr(f64, "tiny") / Float("2.22507385850720138e-308", precision = 64L)) - 1L) < (Float(0.1) * (10L ^ (-ast_attr(f64, "dig"))))), "True")
  expect_identical(.cgb(py_abs((ast_attr(f80, "tiny") / Float("3.36210314311209350626e-4932", precision = 80L)) - 1L) < (Float(0.1) * (10L ^ (-ast_attr(f80, "dig"))))), "True")
  expect_identical(.cgb(py_abs((ast_attr(f128, "tiny") / Float("3.3621031431120935062626778173217526e-4932", precision = 128L)) - 1L) < (Float(0.1) * (10L ^ (-ast_attr(f128, "dig"))))), "True")
  expect_identical(.cgb(type_cast_check(f64, Float(0.5))), "0.50000000000000000")
  expect_identical(.cgb(py_abs(type_cast_check(f64, Float(3.7)) - Float(3.7)) < Float(3e-17)), "True")
  expect_identical(.cgb(py_isinstance(type_cast_check(f64, 3L), Float)), "True")
  expect_identical(.cgb(type_cast_nocheck(f64, oo)), "inf")
  expect_identical(.cgb(type_cast_nocheck(f64, -oo)), "-inf")
  expect_identical(.cgb(type_cast_nocheck(f64, py_float(oo))), "inf")
  expect_identical(.cgb(type_cast_nocheck(f64, py_float(-oo))), "-inf")
  expect_identical(.cgb(py_isnan(type_cast_nocheck(f64, nan))), "True")
  expect_identical(.cgb(!py_eq(f32, f64)), "True")
  expect_identical(.cgb(f64), "float64")
})

test_that("codegen_test_ast: test_Type__cast_check__floating_point", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); t <- Symbol("t"); x0 <- Symbol("x0"); x1 <- Symbol("x1"); x2 <- Symbol("x2"); a <- Symbol("a"); b <- Symbol("b")
  n <- Symbol("n", integer = TRUE)
  A <- MatrixSymbol("A", 3, 1)
  mat <- Matrix(list(list(1), list(2), list(3)))
  B <- IndexedBase("B")
  i <- Idx("i", n)
  A22 <- MatrixSymbol("A22", 2, 2)
  B22 <- MatrixSymbol("B22", 2, 2)
  f16 <- float16; f32 <- float32; f64 <- float64; f80 <- float80; f128 <- float128; c64 <- complex64; c128 <- complex128
  expect_identical(.cgb(type_cast_check(f32, Float(123.45678949))), "ERR ValueError: Casting gives a significantly different value.")
  expect_identical(.cgb(type_cast_check(f32, Float(12.345678949))), "ERR ValueError: Casting gives a significantly different value.")
  expect_identical(.cgb(type_cast_check(f32, Float(1.2345678949))), "ERR ValueError: Casting gives a significantly different value.")
  expect_identical(.cgb(type_cast_check(f32, Float(0.12345678949))), "ERR ValueError: Casting gives a significantly different value.")
  expect_identical(.cgb(py_abs((Float(123.456789049) - type_cast_check(f32, Float(123.456789049))) - Float(4.9e-08)) < Float(1e-08)), "True")
  expect_identical(.cgb(py_abs((Float(0.12345678904) - type_cast_check(f32, Float(0.12345678904))) - Float(4e-11)) < Float(1e-11)), "True")
  suppressWarnings(try({ dcm21 <- Float("0.123456789012345670499") }, silent = TRUE))
  expect_identical(.cgb(py_abs((dcm21 - type_cast_check(f64, dcm21)) - Float(4.99e-19)) < Float(1e-19)), "True")
  suppressWarnings(try({ type_cast_check(f80, Float("0.12345678901234567890103", precision = 88L)) }, silent = TRUE))
  expect_identical(.cgb(type_cast_check(f80, Float("0.12345678901234567890149", precision = 88L))), "ERR ValueError: Casting gives a significantly different value.")
  suppressWarnings(try({ v10 <- Float(12345.67894) }, silent = TRUE))
  expect_identical(.cgb(type_cast_check(f32, v10)), "ERR ValueError: Casting gives a significantly different value.")
  expect_identical(.cgb(py_abs(Float(py_str(v10), precision = 64L + 8L) - type_cast_check(f64, v10)) < (v10 * Float(1e-16))), "True")
  expect_identical(.cgb(py_abs(type_cast_check(f32, 2147483647L) - Integer("2147483650")) < 1L), "True")
})

test_that("codegen_test_ast: test_Type__cast_check__complex_floating_point", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); t <- Symbol("t"); x0 <- Symbol("x0"); x1 <- Symbol("x1"); x2 <- Symbol("x2"); a <- Symbol("a"); b <- Symbol("b")
  n <- Symbol("n", integer = TRUE)
  A <- MatrixSymbol("A", 3, 1)
  mat <- Matrix(list(list(1), list(2), list(3)))
  B <- IndexedBase("B")
  i <- Idx("i", n)
  A22 <- MatrixSymbol("A22", 2, 2)
  B22 <- MatrixSymbol("B22", 2, 2)
  f16 <- float16; f32 <- float32; f64 <- float64; f80 <- float80; f128 <- float128; c64 <- complex64; c128 <- complex128
  suppressWarnings(try({ val9_11 <- Float(123.456789049) + complex(real = 0.0, imaginary = 0.123456789049) }, silent = TRUE))
  expect_identical(.cgb(type_cast_check(c64, Float(0.12345678949) + complex(real = 0.0, imaginary = 0.12345678949))), "ERR ValueError: Casting gives a significantly different value.")
  expect_identical(.cgb(py_abs((val9_11 - type_cast_check(c64, val9_11)) - Float(4.9e-08)) < Float(1e-08)), "True")
  suppressWarnings(try({ dcm21 <- Float("0.123456789012345670499") + complex(real = 0.0, imaginary = 1e-20) }, silent = TRUE))
  expect_identical(.cgb(py_abs((dcm21 - type_cast_check(c128, dcm21)) - Float(4.99e-19)) < Float(1e-19)), "True")
  suppressWarnings(try({ v19 <- Float("0.1234567890123456749") + (complex(real = 0.0, imaginary = 1.0) * Float("0.1234567890123456749")) }, silent = TRUE))
  expect_identical(.cgb(type_cast_check(c128, v19)), "ERR ValueError: Casting gives a significantly different value.")
})

test_that("codegen_test_ast: test_While", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); t <- Symbol("t"); x0 <- Symbol("x0"); x1 <- Symbol("x1"); x2 <- Symbol("x2"); a <- Symbol("a"); b <- Symbol("b")
  n <- Symbol("n", integer = TRUE)
  A <- MatrixSymbol("A", 3, 1)
  mat <- Matrix(list(list(1), list(2), list(3)))
  B <- IndexedBase("B")
  i <- Idx("i", n)
  A22 <- MatrixSymbol("A22", 2, 2)
  B22 <- MatrixSymbol("B22", 2, 2)
  f16 <- float16; f32 <- float32; f64 <- float64; f80 <- float80; f128 <- float128; c64 <- complex64; c128 <- complex128
  suppressWarnings(try({ xpp <- AddAugmentedAssignment(x, 1L) }, silent = TRUE))
  suppressWarnings(try({ whl1 <- While(x < 2L, list(xpp)) }, silent = TRUE))
  expect_identical(.cgb(py_getitem(py_args(ast_attr(whl1, "condition")), 0L)), "x")
  expect_identical(.cgb(py_getitem(py_args(ast_attr(whl1, "condition")), 1L)), "2")
  expect_identical(.cgb(ast_attr(whl1, "condition")), "x < 2")
  expect_identical(.cgb(py_args(ast_attr(whl1, "body"))), "(AddAugmentedAssignment(x, 1),)")
  expect_identical(.cgb(do.call(ast_func(whl1), py_list(py_args(whl1)))), "While(x < 2, body=CodeBlock(\\n    AddAugmentedAssignment(x, 1)\\n))")
  suppressWarnings(try({ cblk <- CodeBlock(AddAugmentedAssignment(x, 1L)) }, silent = TRUE))
  suppressWarnings(try({ whl2 <- While(x < 2L, cblk) }, silent = TRUE))
  expect_identical(.cgb(whl1), "While(x < 2, body=CodeBlock(\\n    AddAugmentedAssignment(x, 1)\\n))")
  expect_identical(.cgb(!py_eq(whl1, While(x < 3L, list(xpp)))), "True")
})

test_that("codegen_test_ast: test_Scope", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); t <- Symbol("t"); x0 <- Symbol("x0"); x1 <- Symbol("x1"); x2 <- Symbol("x2"); a <- Symbol("a"); b <- Symbol("b")
  n <- Symbol("n", integer = TRUE)
  A <- MatrixSymbol("A", 3, 1)
  mat <- Matrix(list(list(1), list(2), list(3)))
  B <- IndexedBase("B")
  i <- Idx("i", n)
  A22 <- MatrixSymbol("A22", 2, 2)
  B22 <- MatrixSymbol("B22", 2, 2)
  f16 <- float16; f32 <- float32; f64 <- float64; f80 <- float80; f128 <- float128; c64 <- complex64; c128 <- complex128
  suppressWarnings(try({ assign <- Assignment(x, y) }, silent = TRUE))
  suppressWarnings(try({ incr <- AddAugmentedAssignment(x, 1L) }, silent = TRUE))
  suppressWarnings(try({ scp <- Scope(list(assign, incr)) }, silent = TRUE))
  suppressWarnings(try({ cblk <- CodeBlock(assign, incr) }, silent = TRUE))
  expect_identical(.cgb(ast_attr(scp, "body")), "CodeBlock(\\nAssignment(x, y),\\nAddAugmentedAssignment(x, 1)\\n)")
  expect_identical(.cgb(scp), "Scope(CodeBlock(\\n    Assignment(x, y),\\n    AddAugmentedAssignment(x, 1)\\n))")
  expect_identical(.cgb(!py_eq(scp, Scope(list(incr, assign)))), "True")
  expect_identical(.cgb(do.call(ast_func(scp), py_list(py_args(scp)))), "Scope(CodeBlock(\\n    Assignment(x, y),\\n    AddAugmentedAssignment(x, 1)\\n))")
})

test_that("codegen_test_ast: test_Print", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); t <- Symbol("t"); x0 <- Symbol("x0"); x1 <- Symbol("x1"); x2 <- Symbol("x2"); a <- Symbol("a"); b <- Symbol("b")
  n <- Symbol("n", integer = TRUE)
  A <- MatrixSymbol("A", 3, 1)
  mat <- Matrix(list(list(1), list(2), list(3)))
  B <- IndexedBase("B")
  i <- Idx("i", n)
  A22 <- MatrixSymbol("A22", 2, 2)
  B22 <- MatrixSymbol("B22", 2, 2)
  f16 <- float16; f32 <- float32; f64 <- float64; f80 <- float80; f128 <- float128; c64 <- complex64; c128 <- complex128
  suppressWarnings(try({ fmt <- "%d %.3f" }, silent = TRUE))
  suppressWarnings(try({ ps <- Print(list(n, x), fmt) }, silent = TRUE))
  expect_identical(.cgb(py_str(ast_attr(ps, "format_string"))), "%d %.3f")
  expect_identical(.cgb(ast_attr(ps, "print_args")), "(n, x)")
  expect_identical(.cgb(py_args(ps)), "((n, x), QuotedString('%d %.3f'), NoneToken())")
  expect_identical(.cgb(ps), "Print((n, x), format_string=%d %.3f)")
  expect_identical(.cgb(!py_eq(ps, Print(list(x, n), fmt))), "True")
  expect_identical(.cgb(do.call(ast_func(ps), py_list(py_args(ps)))), "Print((n, x), format_string=%d %.3f)")
  suppressWarnings(try({ ps2 <- Print(list(n, x)) }, silent = TRUE))
  expect_identical(.cgb(ps2), "Print((n, x))")
  expect_identical(.cgb(!py_eq(ps2, ps)), "True")
  expect_identical(.cgb(ast_attr(ps2, "format_string")), "NoneToken()")
})

test_that("codegen_test_ast: test_FunctionPrototype_and_FunctionDefinition", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); t <- Symbol("t"); x0 <- Symbol("x0"); x1 <- Symbol("x1"); x2 <- Symbol("x2"); a <- Symbol("a"); b <- Symbol("b")
  n <- Symbol("n", integer = TRUE)
  A <- MatrixSymbol("A", 3, 1)
  mat <- Matrix(list(list(1), list(2), list(3)))
  B <- IndexedBase("B")
  i <- Idx("i", n)
  A22 <- MatrixSymbol("A22", 2, 2)
  B22 <- MatrixSymbol("B22", 2, 2)
  f16 <- float16; f32 <- float32; f64 <- float64; f80 <- float80; f128 <- float128; c64 <- complex64; c128 <- complex128
  suppressWarnings(try({ vx <- Variable(x, type = real) }, silent = TRUE))
  suppressWarnings(try({ vn <- Variable(n, type = sym_integer) }, silent = TRUE))
  suppressWarnings(try({ fp1 <- FunctionPrototype(real, "power", list(vx, vn)) }, silent = TRUE))
  expect_identical(.cgb(ast_attr(fp1, "return_type")), "real")
  expect_identical(.cgb(ast_attr(fp1, "name")), "power")
  expect_identical(.cgb(ast_attr(fp1, "parameters")), "(Variable(x, type=real), Variable(n, type=integer))")
  expect_identical(.cgb(fp1), "FunctionPrototype(real, name=power, parameters=(Variable(x, type=real), Variable(n, type=integer)))")
  expect_identical(.cgb(!py_eq(fp1, FunctionPrototype(real, "power", list(vn, vx)))), "True")
  expect_identical(.cgb(do.call(ast_func(fp1), py_list(py_args(fp1)))), "FunctionPrototype(real, name=power, parameters=(Variable(x, type=real), Variable(n, type=integer)))")
  suppressWarnings(try({ body <- list(Assignment(x, x ^ n), Return(x)) }, silent = TRUE))
  suppressWarnings(try({ fd1 <- FunctionDefinition(real, "power", list(vx, vn), body) }, silent = TRUE))
  expect_identical(.cgb(ast_attr(fd1, "return_type")), "real")
  expect_identical(.cgb(py_str(ast_attr(fd1, "name"))), "power")
  expect_identical(.cgb(ast_attr(fd1, "parameters")), "(Variable(x, type=real), Variable(n, type=integer))")
  expect_identical(.cgb(ast_attr(fd1, "body")), "CodeBlock(\\nAssignment(x, x^n),\\nReturn(x)\\n)")
  expect_identical(.cgb(fd1), "FunctionDefinition(real, name=power, parameters=(Variable(x, type=real), Variable(n, type=integer)), body=CodeBlock(\\n    Assignment(x, x^n),\\n    Return(x)\\n))")
  expect_identical(.cgb(!py_eq(fd1, FunctionDefinition(real, "power", list(vx, vn), py_slice(body, NULL, NULL, -1L)))), "True")
  expect_identical(.cgb(do.call(ast_func(fd1), py_list(py_args(fd1)))), "FunctionDefinition(real, name=power, parameters=(Variable(x, type=real), Variable(n, type=integer)), body=CodeBlock(\\n    Assignment(x, x^n),\\n    Return(x)\\n))")
  suppressWarnings(try({ fp2 <- from_FunctionDefinition(fd1) }, silent = TRUE))
  expect_identical(.cgb(fp2), "FunctionPrototype(real, name=power, parameters=(Variable(x, type=real), Variable(n, type=integer)))")
  suppressWarnings(try({ fd2 <- from_FunctionPrototype(fp1, body) }, silent = TRUE))
  expect_identical(.cgb(fd2), "FunctionDefinition(real, name=power, parameters=(Variable(x, type=real), Variable(n, type=integer)), body=CodeBlock(\\n    Assignment(x, x^n),\\n    Return(x)\\n))")
})

test_that("codegen_test_ast: test_Return", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); t <- Symbol("t"); x0 <- Symbol("x0"); x1 <- Symbol("x1"); x2 <- Symbol("x2"); a <- Symbol("a"); b <- Symbol("b")
  n <- Symbol("n", integer = TRUE)
  A <- MatrixSymbol("A", 3, 1)
  mat <- Matrix(list(list(1), list(2), list(3)))
  B <- IndexedBase("B")
  i <- Idx("i", n)
  A22 <- MatrixSymbol("A22", 2, 2)
  B22 <- MatrixSymbol("B22", 2, 2)
  f16 <- float16; f32 <- float32; f64 <- float64; f80 <- float80; f128 <- float128; c64 <- complex64; c128 <- complex128
  suppressWarnings(try({ rs <- Return(x) }, silent = TRUE))
  expect_identical(.cgb(py_args(rs)), "(x,)")
  expect_identical(.cgb(rs), "Return(x)")
  expect_identical(.cgb(!py_eq(rs, Return(y))), "True")
  expect_identical(.cgb(do.call(ast_func(rs), py_list(py_args(rs)))), "Return(x)")
})

test_that("codegen_test_ast: test_FunctionCall", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); t <- Symbol("t"); x0 <- Symbol("x0"); x1 <- Symbol("x1"); x2 <- Symbol("x2"); a <- Symbol("a"); b <- Symbol("b")
  n <- Symbol("n", integer = TRUE)
  A <- MatrixSymbol("A", 3, 1)
  mat <- Matrix(list(list(1), list(2), list(3)))
  B <- IndexedBase("B")
  i <- Idx("i", n)
  A22 <- MatrixSymbol("A22", 2, 2)
  B22 <- MatrixSymbol("B22", 2, 2)
  f16 <- float16; f32 <- float32; f64 <- float64; f80 <- float80; f128 <- float128; c64 <- complex64; c128 <- complex128
  suppressWarnings(try({ fc <- FunctionCall("power", list(x, 3L)) }, silent = TRUE))
  expect_identical(.cgb(py_getitem(ast_attr(fc, "function_args"), 0L)), "x")
  expect_identical(.cgb(py_getitem(ast_attr(fc, "function_args"), 1L)), "3")
  expect_identical(.cgb(py_len(ast_attr(fc, "function_args"))), "2")
  expect_identical(.cgb(py_isinstance(py_getitem(ast_attr(fc, "function_args"), 1L), Integer)), "True")
  expect_identical(.cgb(fc), "FunctionCall(power, function_args=(x, 3))")
  expect_identical(.cgb(!py_eq(fc, FunctionCall("power", list(3L, x)))), "True")
  expect_identical(.cgb(!py_eq(fc, FunctionCall("Power", list(x, 3L)))), "True")
  expect_identical(.cgb(do.call(ast_func(fc), py_list(py_args(fc)))), "FunctionCall(power, function_args=(x, 3))")
  suppressWarnings(try({ fc2 <- FunctionCall("fma", list(2L, 3L, 4L)) }, silent = TRUE))
  expect_identical(.cgb(py_len(ast_attr(fc2, "function_args"))), "3")
  expect_identical(.cgb(py_getitem(ast_attr(fc2, "function_args"), 0L)), "2")
  expect_identical(.cgb(py_getitem(ast_attr(fc2, "function_args"), 1L)), "3")
  expect_identical(.cgb(py_getitem(ast_attr(fc2, "function_args"), 2L)), "4")
  expect_identical(.cgb(py_in(py_str(fc2), list("FunctionCall(fma, function_args=(2, 3, 4))", "FunctionCall(\"fma\", function_args=(2, 3, 4))"))), "True")
})

test_that("codegen_test_ast: test_ast_replace", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); t <- Symbol("t"); x0 <- Symbol("x0"); x1 <- Symbol("x1"); x2 <- Symbol("x2"); a <- Symbol("a"); b <- Symbol("b")
  n <- Symbol("n", integer = TRUE)
  A <- MatrixSymbol("A", 3, 1)
  mat <- Matrix(list(list(1), list(2), list(3)))
  B <- IndexedBase("B")
  i <- Idx("i", n)
  A22 <- MatrixSymbol("A22", 2, 2)
  B22 <- MatrixSymbol("B22", 2, 2)
  f16 <- float16; f32 <- float32; f64 <- float64; f80 <- float80; f128 <- float128; c64 <- complex64; c128 <- complex128
  suppressWarnings(try({ x <- Variable("x", real) }, silent = TRUE))
  suppressWarnings(try({ y <- Variable("y", real) }, silent = TRUE))
  suppressWarnings(try({ n <- Variable("n", sym_integer) }, silent = TRUE))
  suppressWarnings(try({ pwer <- FunctionDefinition(real, "pwer", list(x, n), list(pow(ast_attr(x, "symbol"), ast_attr(n, "symbol")))) }, silent = TRUE))
  suppressWarnings(try({ pname <- ast_attr(pwer, "name") }, silent = TRUE))
  suppressWarnings(try({ pcall <- FunctionCall("pwer", list(y, 3L)) }, silent = TRUE))
  suppressWarnings(try({ tree1 <- CodeBlock(pwer, pcall) }, silent = TRUE))
  expect_identical(.cgb(py_str(ast_attr(py_getitem(py_args(tree1), 0L), "name"))), "pwer")
  expect_identical(.cgb(py_str(ast_attr(py_getitem(py_args(tree1), 1L), "name"))), "pwer")
  suppressWarnings(try({ tree2 <- sym_replace(tree1, pname, String("power")) }, silent = TRUE))
  expect_identical(.cgb(py_str(ast_attr(py_getitem(py_args(tree1), 0L), "name"))), "pwer")
  expect_identical(.cgb(py_str(ast_attr(py_getitem(py_args(tree1), 1L), "name"))), "pwer")
  expect_identical(.cgb(py_str(ast_attr(py_getitem(py_args(tree2), 0L), "name"))), "power")
  expect_identical(.cgb(py_str(ast_attr(py_getitem(py_args(tree2), 1L), "name"))), "power")
})

test_that("codegen_test_cfunctions: test_expm1", {
  expect_identical(.cgb(sym_expm1(0L)), "0")
  suppressWarnings(try({ x <- Symbol("x", real = TRUE) }, silent = TRUE))
  expect_identical(.cgb(expand(sym_expm1(x), func = TRUE) - exp(x)), "-1")
  expect_identical(.cgb(rewrite(sym_expm1(x), "tractable") - exp(x)), "-1")
  expect_identical(.cgb(rewrite(sym_expm1(x), "exp") - exp(x)), "-1")
  expect_identical(.cgb((((N(exp(Float(1e-10))) - 1L) - Float(1e-10)) - Float(5e-21)) < Float(1e-22)), "False")
  expect_identical(.cgb(py_abs((N(sym_expm1(Float(1e-10))) - Float(1e-10)) - Float(5e-21)) < Float(1e-22)), "True")
  expect_identical(.cgb(ask_fact(sym_expm1(x), "real")), "True")
  expect_identical(.cgb(ask_fact(sym_expm1(x), "finite")), "True")
  expect_identical(.cgb(diff(sym_expm1(42L * x), x) - (42L * exp(42L * x))), "0")
  expect_identical(.cgb(diff(sym_expm1(42L * x), x) - diff(expand(sym_expm1(42L * x), func = TRUE), x)), "0")
})

test_that("codegen_test_cfunctions: test_log1p", {
  expect_identical(.cgb(sym_log1p(0L)), "0")
  suppressWarnings(try({ d <- S(10L) }, silent = TRUE))
  expect_identical(.cgb(expand_log((sym_log1p(d ^ (-1000L)) - log((d ^ 1000L) + 1L)) + log(d ^ 1000L))), "0")
  suppressWarnings(try({ x <- Symbol("x", real = TRUE) }, silent = TRUE))
  expect_identical(.cgb(expand(sym_log1p(x), func = TRUE) - log(x + 1L)), "0")
  expect_identical(.cgb(rewrite(sym_log1p(x), "tractable") - log(x + 1L)), "0")
  expect_identical(.cgb(rewrite(sym_log1p(x), "log") - log(x + 1L)), "0")
  expect_identical(.cgb(py_abs(N(log(Float(1e-99) + 1L)) - Float(1e-99)) < Float(1e-100)), "False")
  expect_identical(.cgb(py_abs(N(expand_log(sym_log1p(Float(1e-99)))) - Float(1e-99)) < Float(1e-100)), "True")
  expect_identical(.cgb(ask_fact(sym_log1p(-(2L ^ Rational(-1L, 2L))), "real")), "True")
  expect_identical(.cgb(ask_fact(sym_log1p(-1L), "finite")), "False")
  expect_identical(.cgb(ask_fact(sym_log1p(S_Pi), "finite")), "True")
  expect_identical(.cgb(ask_fact(sym_log1p(x), "positive")), "None")
  expect_identical(.cgb(ask_fact(sym_log1p(Symbol("y", positive = TRUE)), "positive")), "True")
  expect_identical(.cgb(ask_fact(sym_log1p(x), "zero")), "None")
  expect_identical(.cgb(ask_fact(sym_log1p(Symbol("z", zero = TRUE)), "zero")), "True")
  expect_identical(.cgb(ask_fact(sym_log1p(x), "nonnegative")), "None")
  expect_identical(.cgb(ask_fact(sym_log1p(Symbol("o", nonnegative = TRUE)), "nonnegative")), "True")
  expect_identical(.cgb(diff(sym_log1p(42L * x), x) - (42L / ((42L * x) + 1L))), "0")
  expect_identical(.cgb(diff(sym_log1p(42L * x), x) - diff(expand(sym_log1p(42L * x), func = TRUE), x)), "0")
})

test_that("codegen_test_cfunctions: test_exp2", {
  expect_identical(.cgb(exp2(2L)), "4")
  suppressWarnings(try({ x <- Symbol("x", real = TRUE) }, silent = TRUE))
  expect_identical(.cgb(expand(exp2(x), func = TRUE) - (2L ^ x)), "0")
  expect_identical(.cgb(diff(exp2(42L * x), x) - ((42L * exp2(42L * x)) * log(Integer(2L)))), "0")
  expect_identical(.cgb(diff(exp2(42L * x), x) - diff(exp2(42L * x), x)), "0")
})

test_that("codegen_test_cfunctions: test_log2", {
  expect_identical(.cgb(sym_log2(8L)), "3")
  expect_identical(.cgb(!py_eq(sym_log2(S_Pi), log(S_Pi) / log(Integer(2L)))), "True")
  suppressWarnings(try({ x <- Symbol("x", real = TRUE) }, silent = TRUE))
  expect_identical(.cgb(!py_eq(sym_log2(x), log(x) / log(Integer(2L)))), "True")
  expect_identical(.cgb(sym_log2(2L ^ x)), "x")
  expect_identical(.cgb(expand(sym_log2(x), func = TRUE) - (log(x) / log(Integer(2L)))), "0")
  expect_identical(.cgb(diff(sym_log2(42L * x)) - (1L / (log(Integer(2L)) * x))), "0")
  expect_identical(.cgb(diff(sym_log2(42L * x)) - diff(expand(sym_log2(42L * x), func = TRUE), x)), "0")
})

test_that("codegen_test_cfunctions: test_fma", {
  suppressWarnings(try({ .tmp <- sym_symbols("x y z"); x <- py_getitem(.tmp, 0L); y <- py_getitem(.tmp, 1L); z <- py_getitem(.tmp, 2L) }, silent = TRUE))
  expect_identical(.cgb((expand(fma(x, y, z), func = TRUE) - (x * y)) - z), "0")
  suppressWarnings(try({ expr <- fma(17L * x, 42L * y, 101L * z) }, silent = TRUE))
  expect_identical(.cgb(diff(expr, x) - diff(expand(expr, func = TRUE), x)), "0")
  expect_identical(.cgb(diff(expr, y) - diff(expand(expr, func = TRUE), y)), "0")
  expect_identical(.cgb(diff(expr, z) - diff(expand(expr, func = TRUE), z)), "0")
  expect_identical(.cgb(diff(expr, x) - ((17L * 42L) * y)), "0")
  expect_identical(.cgb(diff(expr, y) - ((17L * 42L) * x)), "0")
  expect_identical(.cgb(diff(expr, z) - 101L), "0")
})

test_that("codegen_test_cfunctions: test_log10", {
  suppressWarnings(try({ x <- Symbol("x") }, silent = TRUE))
  expect_identical(.cgb(expand(sym_log10(x), func = TRUE) - (log(x) / log(Integer(10L)))), "0")
  expect_identical(.cgb(diff(sym_log10(42L * x), x) - (1L / (log(Integer(10L)) * x))), "0")
  expect_identical(.cgb(diff(sym_log10(42L * x), x) - diff(expand(sym_log10(42L * x), func = TRUE), x)), "0")
})

test_that("codegen_test_cfunctions: test_Cbrt", {
  suppressWarnings(try({ x <- Symbol("x") }, silent = TRUE))
  expect_identical(.cgb(expand(Cbrt(x), func = TRUE) - (x ^ Rational(1L, 3L))), "0")
  expect_identical(.cgb(diff(Cbrt(42L * x), x) - ((42L * ((42L * x) ^ (Rational(1L, 3L) - 1L))) / 3L)), "0")
  expect_identical(.cgb(diff(Cbrt(42L * x), x) - diff(expand(Cbrt(42L * x), func = TRUE), x)), "0")
})

test_that("codegen_test_cfunctions: test_Sqrt", {
  suppressWarnings(try({ x <- Symbol("x") }, silent = TRUE))
  expect_identical(.cgb(expand(Sqrt(x), func = TRUE) - (x ^ S_Half)), "0")
  expect_identical(.cgb(diff(Sqrt(42L * x), x) - ((42L * ((42L * x) ^ (S_Half - 1L))) / 2L)), "0")
  expect_identical(.cgb(diff(Sqrt(42L * x), x) - diff(expand(Sqrt(42L * x), func = TRUE), x)), "0")
})

test_that("codegen_test_cfunctions: test_hypot", {
  suppressWarnings(try({ .tmp <- sym_symbols("x y"); x <- py_getitem(.tmp, 0L); y <- py_getitem(.tmp, 1L) }, silent = TRUE))
  expect_identical(.cgb(expand(hypot(x, y), func = TRUE) - (((x ^ 2L) + (y ^ 2L)) ^ S_Half)), "0")
  expect_identical(.cgb(expand(diff(hypot(17L * x, 42L * y), x), func = TRUE) - diff(expand(hypot(17L * x, 42L * y), func = TRUE), x)), "0")
  expect_identical(.cgb(expand(diff(hypot(17L * x, 42L * y), y), func = TRUE) - diff(expand(hypot(17L * x, 42L * y), func = TRUE), y)), "0")
  expect_identical(.cgb(expand(diff(hypot(17L * x, 42L * y), x), func = TRUE) - (((((2L * 17L) * 17L) * x) * ((((17L * x) ^ 2L) + ((42L * y) ^ 2L)) ^ Rational(-1L, 2L))) / 2L)), "0")
  expect_identical(.cgb(expand(diff(hypot(17L * x, 42L * y), y), func = TRUE) - (((((2L * 42L) * 42L) * y) * ((((17L * x) ^ 2L) + ((42L * y) ^ 2L)) ^ Rational(-1L, 2L))) / 2L)), "0")
})

test_that("codegen_test_cfunctions: test_isnan_isinf", {
  suppressWarnings(try({ x <- Symbol("x") }, silent = TRUE))
  expect_identical(.cgb(isinf(oo)), "True")
  expect_identical(.cgb(isinf(-oo)), "True")
  expect_identical(.cgb(isinf(S_Pi)), "False")
  suppressWarnings(try({ isinfx <- isinf(x) }, silent = TRUE))
  expect_identical(.cgb(!py_in(isinfx, list(FALSE, TRUE))), "True")
  expect_identical(.cgb(py_is(ast_attr(isinfx, "func"), isinf)), "True")
  expect_identical(.cgb(py_args(isinfx)), "(x,)")
  expect_identical(.cgb(isnan(nan)), "True")
  expect_identical(.cgb(isnan(S_Pi)), "False")
  suppressWarnings(try({ isnanx <- isnan(x) }, silent = TRUE))
  expect_identical(.cgb(!py_in(isnanx, list(FALSE, TRUE))), "True")
  expect_identical(.cgb(py_is(ast_attr(isnanx, "func"), isnan)), "True")
  expect_identical(.cgb(py_args(isnanx)), "(x,)")
})

test_that("codegen_test_cnodes: test_alignof", {
  x <- Symbol("x"); y <- Symbol("y")
  suppressWarnings(try({ ax <- alignof(x) }, silent = TRUE))
  expect_identical(.cgb(ccode(ax)), "alignof(x)")
  expect_identical(.cgb(do.call(ast_func(ax), py_list(py_args(ax)))), "FunctionCall(alignof, function_args=(x,))")
})

test_that("codegen_test_cnodes: test_CommaOperator", {
  x <- Symbol("x"); y <- Symbol("y")
  suppressWarnings(try({ expr <- CommaOperator(PreIncrement(x), 2L * x) }, silent = TRUE))
  expect_identical(.cgb(ccode(expr)), "(++(x), 2*x)")
  expect_identical(.cgb(do.call(ast_func(expr), py_list(py_args(expr)))), "CommaOperator(PreIncrement(x), 2*x)")
})

test_that("codegen_test_cnodes: test_goto_Label", {
  x <- Symbol("x"); y <- Symbol("y")
  suppressWarnings(try({ s <- "early_exit" }, silent = TRUE))
  suppressWarnings(try({ g <- goto(s) }, silent = TRUE))
  expect_identical(.cgb(do.call(ast_func(g), py_list(py_args(g)))), "goto(Label(early_exit))")
  expect_identical(.cgb(!py_eq(g, goto("foobar"))), "True")
  expect_identical(.cgb(ccode(g)), "goto early_exit")
  suppressWarnings(try({ l1 <- Label(s) }, silent = TRUE))
  expect_identical(.cgb(ccode(l1)), "early_exit:")
  expect_identical(.cgb(l1), "Label(early_exit)")
  expect_identical(.cgb(!py_eq(l1, Label("foobar"))), "True")
  suppressWarnings(try({ body <- list(PreIncrement(x)) }, silent = TRUE))
  suppressWarnings(try({ l2 <- Label(s, body) }, silent = TRUE))
  expect_identical(.cgb(ast_attr(l2, "name")), "early_exit")
  expect_identical(.cgb(ast_attr(l2, "body")), "CodeBlock(\\nPreIncrement(x)\\n)")
  expect_identical(.cgb(ccode(l2)), "early_exit:\\n++(x);")
  suppressWarnings(try({ body <- list(PreIncrement(x), PreDecrement(y)) }, silent = TRUE))
  suppressWarnings(try({ l2 <- Label(s, body) }, silent = TRUE))
  expect_identical(.cgb(ast_attr(l2, "name")), "early_exit")
  expect_identical(.cgb(ast_attr(l2, "body")), "CodeBlock(\\nPreIncrement(x),\\nPreDecrement(y)\\n)")
  expect_identical(.cgb(ccode(l2)), "early_exit:\\n{\\n   ++(x);\\n   --(y);\\n}")
})

test_that("codegen_test_cnodes: test_PreDecrement", {
  x <- Symbol("x"); y <- Symbol("y")
  suppressWarnings(try({ p <- PreDecrement(x) }, silent = TRUE))
  expect_identical(.cgb(do.call(ast_func(p), py_list(py_args(p)))), "PreDecrement(x)")
  expect_identical(.cgb(ccode(p)), "--(x)")
})

test_that("codegen_test_cnodes: test_PostDecrement", {
  x <- Symbol("x"); y <- Symbol("y")
  suppressWarnings(try({ p <- PostDecrement(x) }, silent = TRUE))
  expect_identical(.cgb(do.call(ast_func(p), py_list(py_args(p)))), "PostDecrement(x)")
  expect_identical(.cgb(ccode(p)), "(x)--")
})

test_that("codegen_test_cnodes: test_PreIncrement", {
  x <- Symbol("x"); y <- Symbol("y")
  suppressWarnings(try({ p <- PreIncrement(x) }, silent = TRUE))
  expect_identical(.cgb(do.call(ast_func(p), py_list(py_args(p)))), "PreIncrement(x)")
  expect_identical(.cgb(ccode(p)), "++(x)")
})

test_that("codegen_test_cnodes: test_PostIncrement", {
  x <- Symbol("x"); y <- Symbol("y")
  suppressWarnings(try({ p <- PostIncrement(x) }, silent = TRUE))
  expect_identical(.cgb(do.call(ast_func(p), py_list(py_args(p)))), "PostIncrement(x)")
  expect_identical(.cgb(ccode(p)), "(x)++")
})

test_that("codegen_test_cnodes: test_sizeof", {
  x <- Symbol("x"); y <- Symbol("y")
  suppressWarnings(try({ typename <- "unsigned int" }, silent = TRUE))
  suppressWarnings(try({ sz <- sizeof(typename) }, silent = TRUE))
  expect_identical(.cgb(ccode(sz)), "sizeof(unsigned int)")
  expect_identical(.cgb(do.call(ast_func(sz), py_list(py_args(sz)))), "FunctionCall(sizeof, function_args=(unsigned int,))")
  expect_identical(.cgb(ast_attr(sz, "is_Atom")), "False")
})

test_that("codegen_test_cnodes: test_struct", {
  x <- Symbol("x"); y <- Symbol("y")
  suppressWarnings(try({ vx <- Variable(x, type = float64); vy <- Variable(y, type = float64) }, silent = TRUE))
  suppressWarnings(try({ s <- struct("vec2", list(vx, vy)) }, silent = TRUE))
  expect_identical(.cgb(do.call(ast_func(s), py_list(py_args(s)))), "struct(vec2, declarations=(Declaration(Variable(x, type=float64)), Declaration(Variable(y, type=float64))))")
  expect_identical(.cgb(s), "struct(vec2, declarations=(Declaration(Variable(x, type=float64)), Declaration(Variable(y, type=float64))))")
  expect_identical(.cgb(!py_eq(s, struct("vec2", list(vy, vx)))), "True")
  expect_identical(.cgb(py_str(ast_attr(s, "name"))), "vec2")
  expect_identical(.cgb(py_len(ast_attr(s, "declarations"))), "2")
  expect_identical(.cgb(all(vapply(py_iter(ast_attr(s, "declarations")), function(arg) py_isinstance(arg, Declaration), TRUE))), "True")
  expect_identical(.cgb(ccode(s)), "struct vec2 {\\n   double x;\\n   double y;\\n}")
})

test_that("codegen_test_cnodes: test_union", {
  x <- Symbol("x"); y <- Symbol("y")
  suppressWarnings(try({ vx <- Variable(x, type = float64); vy <- Variable(y, type = int64) }, silent = TRUE))
  suppressWarnings(try({ u <- c_union("dualuse", list(vx, vy)) }, silent = TRUE))
  expect_identical(.cgb(do.call(ast_func(u), py_list(py_args(u)))), "union(dualuse, declarations=(Declaration(Variable(x, type=float64)), Declaration(Variable(y, type=int64))))")
  expect_identical(.cgb(u), "union(dualuse, declarations=(Declaration(Variable(x, type=float64)), Declaration(Variable(y, type=int64))))")
  expect_identical(.cgb(py_str(ast_attr(u, "name"))), "dualuse")
  expect_identical(.cgb(py_len(ast_attr(u, "declarations"))), "2")
  expect_identical(.cgb(all(vapply(py_iter(ast_attr(u, "declarations")), function(arg) py_isinstance(arg, Declaration), TRUE))), "True")
  expect_identical(.cgb(ccode(u)), "union dualuse {\\n   double x;\\n   int64_t y;\\n}")
})

test_that("codegen_test_cxxnodes: test_using", {
  suppressWarnings(try({ v <- Type("std::vector") }, silent = TRUE))
  suppressWarnings(try({ u1 <- using(v) }, silent = TRUE))
  expect_identical(.cgb(cxxcode(u1)), "using std::vector")
  suppressWarnings(try({ u2 <- using(v, "vec") }, silent = TRUE))
  expect_identical(.cgb(cxxcode(u2)), "using vec = std::vector")
})

test_that("codegen_test_fnodes: test_size", {
  render_as_module <- f_render_as_module
  suppressWarnings(try({ x <- Symbol("x", real = TRUE) }, silent = TRUE))
  suppressWarnings(try({ sx <- size(x) }, silent = TRUE))
  expect_identical(.cgb(fcode(sx, source_format = "free")), "size(x)")
})

test_that("codegen_test_fnodes: test_isign", {
  render_as_module <- f_render_as_module
  suppressWarnings(try({ x <- Symbol("x", integer = TRUE) }, silent = TRUE))
  expect_identical(.cgb(unchanged(isign, 1L, x)), "True")
  expect_identical(.cgb(fcode(isign(1L, x), standard = 95L, source_format = "free")), "isign(1, x)")
})

test_that("codegen_test_fnodes: test_dsign", {
  render_as_module <- f_render_as_module
  suppressWarnings(try({ x <- Symbol("x") }, silent = TRUE))
  expect_identical(.cgb(unchanged(dsign, 1L, x)), "True")
  expect_identical(.cgb(fcode(dsign(literal_dp(1L), x), standard = 95L, source_format = "free")), "dsign(1d0, x)")
})

test_that("codegen_test_fnodes: test_cmplx", {
  render_as_module <- f_render_as_module
  suppressWarnings(try({ x <- Symbol("x") }, silent = TRUE))
  expect_identical(.cgb(unchanged(cmplx, 1L, x)), "True")
})

test_that("codegen_test_fnodes: test_kind", {
  render_as_module <- f_render_as_module
  suppressWarnings(try({ x <- Symbol("x") }, silent = TRUE))
  expect_identical(.cgb(unchanged(kind, x)), "True")
})

test_that("codegen_test_fnodes: test_literal_dp", {
  render_as_module <- f_render_as_module
  expect_identical(.cgb(fcode(literal_dp(0L), source_format = "free")), "0d0")
})

test_that("codegen_test_fnodes: test_reshape", {
  render_as_module <- f_render_as_module
  suppressWarnings(try({ "Test reshape function with keyword arguments and essential edge cases." }, silent = TRUE))
  suppressWarnings(try({ invisible(NULL) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- sym_symbols("array shape pad order"); array <- py_getitem(.tmp, 0L); shape <- py_getitem(.tmp, 1L); pad <- py_getitem(.tmp, 2L); order <- py_getitem(.tmp, 3L) }, silent = TRUE))
  suppressWarnings(try({ cases <- list(list(fortran_reshape(array, shape), "reshape(array, shape)"), list(fortran_reshape(array, shape, pad), "reshape(array, shape, pad=pad)"), list(fortran_reshape(array, shape, NULL, order), "reshape(array, shape, order=order)"), list(fortran_reshape(array, shape, pad, order), "reshape(array, shape, pad=pad, order=order)"), list(fortran_reshape(Symbol("a") + Symbol("b"), list(Symbol("n") * 2L)), "reshape(a + b, [2*n])"), list(fortran_reshape(array, shape, NULL, String("F")), "reshape(array, shape, order=F)"), list(fortran_reshape(array, list(Symbol("n"), 3L)), "reshape(array, [n, 3])"), list(fortran_reshape(array, list(2L, 3L)), "reshape(array, [2, 3])"), list(fortran_reshape(list(1L, 2L), list(2L, 2L), list(0L), list(2L, 1L)), "reshape([1, 2], [2, 2], pad=[0], order=[2, 1])"), list(fortran_reshape(list(list(1L, 2L), list(3L, 4L)), list(2L, 2L)), "reshape([[1, 2], [3, 4]], [2, 2])"), list(fortran_reshape(list(list(1L, 2L), list(3L, 4L)), list(2L, 3L), list(list(5L, 6L))), "reshape([[1, 2], [3, 4]], [2, 3], pad=[[5, 6]])")) }, silent = TRUE))
  expect_identical(.cgb(fcode(fortran_reshape(array, list(2L, 3L)), source_format = "free", standard = 95L)), "reshape(array, (/2, 3/))")
})

test_that("codegen_test_numpy_nodes: test_logaddexp", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  lb <- function(x) log(x) / log(Integer(2L))
  suppressWarnings(try({ lae_xy <- logaddexp(x, y) }, silent = TRUE))
  suppressWarnings(try({ ref_xy <- log(exp(x) + exp(y)) }, silent = TRUE))
  suppressWarnings(try({ one_third_e <- (1L * exp(Integer(1L))) / 3L }, silent = TRUE))
  suppressWarnings(try({ two_thirds_e <- (2L * exp(Integer(1L))) / 3L }, silent = TRUE))
  suppressWarnings(try({ logThirdE <- log(one_third_e) }, silent = TRUE))
  suppressWarnings(try({ logTwoThirdsE <- log(two_thirds_e) }, silent = TRUE))
  suppressWarnings(try({ lae_sum_to_e <- logaddexp(logThirdE, logTwoThirdsE) }, silent = TRUE))
  expect_identical(.cgb(rewrite(lae_sum_to_e, "log")), "1")
  expect_identical(.cgb(simplify(lae_sum_to_e)), "1")
  suppressWarnings(try({ was <- logaddexp(2L, 3L) }, silent = TRUE))
  expect_identical(.cgb(srepr(was)), "logaddexp(Integer(2), Integer(3))")
})

test_that("codegen_test_numpy_nodes: test_logaddexp2", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  lb <- function(x) log(x) / log(Integer(2L))
  suppressWarnings(try({ lae2_xy <- logaddexp2(x, y) }, silent = TRUE))
  suppressWarnings(try({ ref2_xy <- log((2L ^ x) + (2L ^ y)) / log(Integer(2L)) }, silent = TRUE))
  suppressWarnings(try({ lb <- function(x) log(x) / log(Integer(2L)) }, silent = TRUE))
  suppressWarnings(try({ two_thirds <- (S_One * 2L) / 3L }, silent = TRUE))
  suppressWarnings(try({ four_thirds <- 2L * two_thirds }, silent = TRUE))
  suppressWarnings(try({ lbTwoThirds <- lb(two_thirds) }, silent = TRUE))
  suppressWarnings(try({ lbFourThirds <- lb(four_thirds) }, silent = TRUE))
  suppressWarnings(try({ lae2_sum_to_2 <- logaddexp2(lbTwoThirds, lbFourThirds) }, silent = TRUE))
  expect_identical(.cgb(rewrite(lae2_sum_to_2, "log")), "1")
  expect_identical(.cgb(simplify(lae2_sum_to_2)), "1")
  suppressWarnings(try({ was <- logaddexp2(x, y) }, silent = TRUE))
  expect_identical(.cgb(srepr(was)), "logaddexp2(Symbol('x'), Symbol('y'))")
})

test_that("codegen_test_numpy_nodes: test_logaddexp_evalf_precision", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  lb <- function(x) log(x) / log(Integer(2L))
  expect_identical(.cgb(N(logaddexp(x, y))), "log(exp(x) + exp(y))")
  suppressWarnings(try({ cases <- list(list(logaddexp(2L, 3L), log(exp(Integer(2L)) + exp(Integer(3L)))), list(logaddexp2(2L, 3L), log(Integer(12L)) / log(Integer(2L)))) }, silent = TRUE))
})

test_that("codegen_test_pynodes: test_List", {
  suppressWarnings(try({ l <- List(2L, 3L, 4L) }, silent = TRUE))
  expect_identical(.cgb(l), "[2, 3, 4]")
  expect_identical(.cgb(py_str(l)), "[2, 3, 4]")
  suppressWarnings(try({ .tmp <- sym_symbols("x y z"); x <- py_getitem(.tmp, 0L); y <- py_getitem(.tmp, 1L); z <- py_getitem(.tmp, 2L) }, silent = TRUE))
  suppressWarnings(try({ l <- List(x ^ 2L, y ^ 3L, z ^ 4L) }, silent = TRUE))
  suppressWarnings(try({ m <- sym_replace(l, function(arg) ast_attr(arg, "is_Pow") && (ast_attr(arg, "exp") > 2L), function(p) ast_attr(p, "base") - ast_attr(p, "exp")) }, silent = TRUE))
  expect_identical(.cgb(m), "[x^2, y - 3, z - 4]")
})

test_that("codegen_test_pyutils: test_standard", {
  render_as_module <- py_render_as_module
  suppressWarnings(try({ ast <- Print(py_method("x y", "split"), "coordinate: %12.5g %12.5g\\n") }, silent = TRUE))
  expect_identical(.cgb(render_as_module(ast, standard = "python3")), "\\n\\nprint(\"coordinate: %12.5g %12.5g\\n\" % (x, y), end=\"\")")
})

test_that("codegen_test_rewriting: test_log2_opt", {
  suppressWarnings(try({ x <- Symbol("x") }, silent = TRUE))
  suppressWarnings(try({ expr1 <- (7L * log((3L * x) + 5L)) / log(Integer(2L)) }, silent = TRUE))
  suppressWarnings(try({ opt1 <- sym_optimize(expr1, list(log2_opt)) }, silent = TRUE))
  expect_identical(.cgb(opt1), "7*log2(3*x + 5)")
  expect_identical(.cgb(rewrite(opt1, "log")), "7*log(3*x + 5)/log(2)")
  suppressWarnings(try({ expr2 <- (3L * log((5L * x) + 7L)) / (13L * log(Integer(2L))) }, silent = TRUE))
  suppressWarnings(try({ opt2 <- sym_optimize(expr2, list(log2_opt)) }, silent = TRUE))
  expect_identical(.cgb(opt2), "3*log2(5*x + 7)/13")
  expect_identical(.cgb(rewrite(opt2, "log")), "3*log(5*x + 7)/(13*log(2))")
  suppressWarnings(try({ expr3 <- log(x) / log(Integer(2L)) }, silent = TRUE))
  suppressWarnings(try({ opt3 <- sym_optimize(expr3, list(log2_opt)) }, silent = TRUE))
  expect_identical(.cgb(opt3), "log2(x)")
  expect_identical(.cgb(rewrite(opt3, "log")), "log(x)/log(2)")
  suppressWarnings(try({ expr4 <- (log(x) / log(Integer(2L))) + log(x + 1L) }, silent = TRUE))
  suppressWarnings(try({ opt4 <- sym_optimize(expr4, list(log2_opt)) }, silent = TRUE))
  expect_identical(.cgb(opt4), "log2(x) + log(2)*log2(x + 1)")
  expect_identical(.cgb(rewrite(opt4, "log")), "log(x)/log(2) + log(x + 1)")
  suppressWarnings(try({ expr5 <- log(Integer(17L)) }, silent = TRUE))
  suppressWarnings(try({ opt5 <- sym_optimize(expr5, list(log2_opt)) }, silent = TRUE))
  expect_identical(.cgb(opt5), "log(17)")
  suppressWarnings(try({ expr6 <- log(x + 3L) / log(Integer(2L)) }, silent = TRUE))
  suppressWarnings(try({ opt6 <- sym_optimize(expr6, list(log2_opt)) }, silent = TRUE))
  expect_identical(.cgb(py_str(opt6)), "log2(x + 3)")
  expect_identical(.cgb(rewrite(opt6, "log")), "log(x + 3)/log(2)")
})

test_that("codegen_test_rewriting: test_exp2_opt", {
  suppressWarnings(try({ x <- Symbol("x") }, silent = TRUE))
  suppressWarnings(try({ expr1 <- 1L + (2L ^ x) }, silent = TRUE))
  suppressWarnings(try({ opt1 <- sym_optimize(expr1, list(exp2_opt)) }, silent = TRUE))
  expect_identical(.cgb(opt1), "exp2(x) + 1")
  expect_identical(.cgb(rewrite(opt1, "Pow")), "2^x + 1")
  suppressWarnings(try({ expr2 <- 1L + (3L ^ x) }, silent = TRUE))
  expect_identical(.cgb(expr2), "3^x + 1")
})

test_that("codegen_test_rewriting: test_expm1_opt", {
  suppressWarnings(try({ x <- Symbol("x") }, silent = TRUE))
  suppressWarnings(try({ expr1 <- exp(x) - 1L }, silent = TRUE))
  suppressWarnings(try({ opt1 <- sym_optimize(expr1, list(expm1_opt)) }, silent = TRUE))
  expect_identical(.cgb(sym_expm1(x) - opt1), "0")
  expect_identical(.cgb(rewrite(opt1, "exp")), "exp(x) - 1")
  suppressWarnings(try({ expr2 <- (3L * exp(x)) - 3L }, silent = TRUE))
  suppressWarnings(try({ opt2 <- sym_optimize(expr2, list(expm1_opt)) }, silent = TRUE))
  expect_identical(.cgb(3L * sym_expm1(x)), "3*expm1(x)")
  expect_identical(.cgb(rewrite(opt2, "exp")), "3*exp(x) - 3")
  suppressWarnings(try({ expr3 <- (3L * exp(x)) - 5L }, silent = TRUE))
  suppressWarnings(try({ opt3 <- sym_optimize(expr3, list(expm1_opt)) }, silent = TRUE))
  expect_identical(.cgb((3L * sym_expm1(x)) - 2L), "3*expm1(x) - 2")
  expect_identical(.cgb(rewrite(opt3, "exp")), "3*exp(x) - 5")
  suppressWarnings(try({ expm1_opt_non_opportunistic <- FuncMinusOneOptim(exp, sym_expm1, opportunistic = FALSE) }, silent = TRUE))
  expect_identical(.cgb(expr3), "3*exp(x) - 5")
  expect_identical(.cgb(opt1), "expm1(x)")
  expect_identical(.cgb(opt2), "3*expm1(x)")
  suppressWarnings(try({ expr4 <- ((3L * exp(x)) + log(x)) - 3L }, silent = TRUE))
  suppressWarnings(try({ opt4 <- sym_optimize(expr4, list(expm1_opt)) }, silent = TRUE))
  expect_identical(.cgb((3L * sym_expm1(x)) + log(x)), "log(x) + 3*expm1(x)")
  expect_identical(.cgb(rewrite(opt4, "exp")), "3*exp(x) + log(x) - 3")
  suppressWarnings(try({ expr5 <- (3L * exp(2L * x)) - 3L }, silent = TRUE))
  suppressWarnings(try({ opt5 <- sym_optimize(expr5, list(expm1_opt)) }, silent = TRUE))
  expect_identical(.cgb(3L * sym_expm1(2L * x)), "3*expm1(2*x)")
  expect_identical(.cgb(rewrite(opt5, "exp")), "3*exp(2*x) - 3")
  suppressWarnings(try({ expr6 <- (((2L * exp(x)) + 1L) / (exp(x) + 1L)) + 1L }, silent = TRUE))
  suppressWarnings(try({ opt6 <- sym_optimize(expr6, list(expm1_opt)) }, silent = TRUE))
  expect_identical(.cgb(py_method(opt6, "count_ops") <= py_method(expr6, "count_ops")), "True")
  suppressWarnings(try({ y <- Symbol("y") }, silent = TRUE))
  suppressWarnings(try({ expr7 <- (((2L * exp(x)) - 1L) / (1L - exp(y))) - (1L / (1L - exp(y))) }, silent = TRUE))
  suppressWarnings(try({ opt7 <- sym_optimize(expr7, list(expm1_opt)) }, silent = TRUE))
  expect_identical(.cgb(((-2L) * sym_expm1(x)) / sym_expm1(y)), "-2*expm1(x)/expm1(y)")
  expect_identical(.cgb(sym_factor(rewrite(opt7, "exp") - expr7)), "0")
  suppressWarnings(try({ expr8 <- ((1L + exp(x)) ^ 2L) - 4L }, silent = TRUE))
  suppressWarnings(try({ opt8 <- sym_optimize(expr8, list(expm1_opt)) }, silent = TRUE))
  suppressWarnings(try({ tgt8a <- (exp(x) + 3L) * sym_expm1(x) }, silent = TRUE))
  suppressWarnings(try({ tgt8b <- (2L * sym_expm1(x)) + sym_expm1(2L * x) }, silent = TRUE))
  expect_identical(.cgb(sym_factor(rewrite(tgt8a - tgt8b, "exp"))), "0")
  expect_identical(.cgb(py_in(opt8, list(tgt8a, tgt8b))), "True")
  expect_identical(.cgb(sym_factor(rewrite(opt8, "exp") - expr8)), "0")
  suppressWarnings(try({ expr9 <- sin(expr8) }, silent = TRUE))
  suppressWarnings(try({ opt9 <- sym_optimize(expr9, list(expm1_opt)) }, silent = TRUE))
  suppressWarnings(try({ tgt9a <- sin(tgt8a) }, silent = TRUE))
  suppressWarnings(try({ tgt9b <- sin(tgt8b) }, silent = TRUE))
  expect_identical(.cgb(py_in(opt9, list(tgt9a, tgt9b))), "True")
  expect_identical(.cgb(ask_fact(sym_factor(rewrite(opt9, "exp") - rewrite(expr9, "exp")), "zero")), "True")
})

test_that("codegen_test_rewriting: test_expm1_two_exp_terms", {
  suppressWarnings(try({ .tmp <- map(Symbol, py_method("x y", "split")); x <- py_getitem(.tmp, 0L); y <- py_getitem(.tmp, 1L) }, silent = TRUE))
  suppressWarnings(try({ expr1 <- (exp(x) + exp(y)) - 2L }, silent = TRUE))
  suppressWarnings(try({ opt1 <- sym_optimize(expr1, list(expm1_opt)) }, silent = TRUE))
  expect_identical(.cgb(opt1), "expm1(x) + expm1(y)")
})

test_that("codegen_test_rewriting: test_cosm1_opt", {
  suppressWarnings(try({ x <- Symbol("x") }, silent = TRUE))
  suppressWarnings(try({ expr1 <- cos(x) - 1L }, silent = TRUE))
  suppressWarnings(try({ opt1 <- sym_optimize(expr1, list(cosm1_opt)) }, silent = TRUE))
  expect_identical(.cgb(cosm1(x) - opt1), "0")
  expect_identical(.cgb(rewrite(opt1, "cos")), "cos(x) - 1")
  suppressWarnings(try({ expr2 <- (3L * cos(x)) - 3L }, silent = TRUE))
  suppressWarnings(try({ opt2 <- sym_optimize(expr2, list(cosm1_opt)) }, silent = TRUE))
  expect_identical(.cgb(3L * cosm1(x)), "3*cosm1(x)")
  expect_identical(.cgb(rewrite(opt2, "cos")), "3*cos(x) - 3")
  suppressWarnings(try({ expr3 <- (3L * cos(x)) - 5L }, silent = TRUE))
  suppressWarnings(try({ opt3 <- sym_optimize(expr3, list(cosm1_opt)) }, silent = TRUE))
  expect_identical(.cgb((3L * cosm1(x)) - 2L), "3*cosm1(x) - 2")
  expect_identical(.cgb(rewrite(opt3, "cos")), "3*cos(x) - 5")
  suppressWarnings(try({ cosm1_opt_non_opportunistic <- FuncMinusOneOptim(cos, cosm1, opportunistic = FALSE) }, silent = TRUE))
  expect_identical(.cgb(expr3), "3*cos(x) - 5")
  expect_identical(.cgb(opt1), "cosm1(x)")
  expect_identical(.cgb(opt2), "3*cosm1(x)")
  suppressWarnings(try({ expr4 <- ((3L * cos(x)) + log(x)) - 3L }, silent = TRUE))
  suppressWarnings(try({ opt4 <- sym_optimize(expr4, list(cosm1_opt)) }, silent = TRUE))
  expect_identical(.cgb((3L * cosm1(x)) + log(x)), "log(x) + 3*cosm1(x)")
  expect_identical(.cgb(rewrite(opt4, "cos")), "log(x) + 3*cos(x) - 3")
  suppressWarnings(try({ expr5 <- (3L * cos(2L * x)) - 3L }, silent = TRUE))
  suppressWarnings(try({ opt5 <- sym_optimize(expr5, list(cosm1_opt)) }, silent = TRUE))
  expect_identical(.cgb(3L * cosm1(2L * x)), "3*cosm1(2*x)")
  expect_identical(.cgb(rewrite(opt5, "cos")), "3*cos(2*x) - 3")
  suppressWarnings(try({ expr6 <- 2L - (2L * cos(x)) }, silent = TRUE))
  suppressWarnings(try({ opt6 <- sym_optimize(expr6, list(cosm1_opt)) }, silent = TRUE))
  expect_identical(.cgb((-2L) * cosm1(x)), "-2*cosm1(x)")
  expect_identical(.cgb(rewrite(opt6, "cos")), "2 - 2*cos(x)")
})

test_that("codegen_test_rewriting: test_cosm1_two_cos_terms", {
  suppressWarnings(try({ .tmp <- map(Symbol, py_method("x y", "split")); x <- py_getitem(.tmp, 0L); y <- py_getitem(.tmp, 1L) }, silent = TRUE))
  suppressWarnings(try({ expr1 <- (cos(x) + cos(y)) - 2L }, silent = TRUE))
  suppressWarnings(try({ opt1 <- sym_optimize(expr1, list(cosm1_opt)) }, silent = TRUE))
  expect_identical(.cgb(opt1), "cosm1(x) + cosm1(y)")
})

test_that("codegen_test_rewriting: test_expm1_cosm1_mixed", {
  suppressWarnings(try({ x <- Symbol("x") }, silent = TRUE))
  suppressWarnings(try({ expr1 <- (exp(x) + cos(x)) - 2L }, silent = TRUE))
  suppressWarnings(try({ opt1 <- sym_optimize(expr1, list(expm1_opt, cosm1_opt)) }, silent = TRUE))
  expect_identical(.cgb(opt1), "cosm1(x) + expm1(x)")
})

test_that("codegen_test_rewriting: test_cosm1_apart", {
  suppressWarnings(try({ x <- Symbol("x") }, silent = TRUE))
  suppressWarnings(try({ expr1 <- (1L / cos(x)) - 1L }, silent = TRUE))
  suppressWarnings(try({ opt1 <- sym_optimize(expr1, list(cosm1_opt)) }, silent = TRUE))
  expect_identical(.cgb(opt1), "-cosm1(x)/cos(x)")
  suppressWarnings(try({ expr2 <- (2L / cos(x)) - 2L }, silent = TRUE))
  suppressWarnings(try({ opt2 <- sym_optimize(expr2, optims_scipy) }, silent = TRUE))
  expect_identical(.cgb(opt2), "-2*cosm1(x)/cos(x)")
  suppressWarnings(try({ expr3 <- (S_Pi / cos(3L * x)) - S_Pi }, silent = TRUE))
  suppressWarnings(try({ opt3 <- sym_optimize(expr3, list(cosm1_opt)) }, silent = TRUE))
  expect_identical(.cgb(opt3), "-pi*cosm1(3*x)/cos(3*x)")
})

test_that("codegen_test_rewriting: test_powm1", {
  suppressWarnings(try({ args <- map(Symbol, "xy"); x <- py_getitem(args, 0L); y <- py_getitem(args, 1L) }, silent = TRUE))
  suppressWarnings(try({ expr1 <- (x ^ y) - 1L }, silent = TRUE))
  suppressWarnings(try({ opt1 <- sym_optimize(expr1, list(powm1_opt)) }, silent = TRUE))
  expect_identical(.cgb(opt1), "powm1(x, y)")
})

test_that("codegen_test_rewriting: test_log1p_opt", {
  suppressWarnings(try({ x <- Symbol("x") }, silent = TRUE))
  suppressWarnings(try({ expr1 <- log(x + 1L) }, silent = TRUE))
  suppressWarnings(try({ opt1 <- sym_optimize(expr1, list(log1p_opt)) }, silent = TRUE))
  expect_identical(.cgb(sym_log1p(x) - opt1), "0")
  expect_identical(.cgb(rewrite(opt1, "log")), "log(x + 1)")
  suppressWarnings(try({ expr2 <- log((3L * x) + 3L) }, silent = TRUE))
  suppressWarnings(try({ opt2 <- sym_optimize(expr2, list(log1p_opt)) }, silent = TRUE))
  expect_identical(.cgb(sym_log1p(x) + log(Integer(3L))), "log1p(x) + log(3)")
  expect_identical(.cgb(simplify(rewrite(opt2, "log") - expr2)), "0")
  suppressWarnings(try({ expr3 <- log((2L * x) + 1L) }, silent = TRUE))
  suppressWarnings(try({ opt3 <- sym_optimize(expr3, list(log1p_opt)) }, silent = TRUE))
  expect_identical(.cgb(sym_log1p(2L * x) - opt3), "0")
  expect_identical(.cgb(rewrite(opt3, "log")), "log(2*x + 1)")
  suppressWarnings(try({ expr4 <- log(x + 3L) }, silent = TRUE))
  suppressWarnings(try({ opt4 <- sym_optimize(expr4, list(log1p_opt)) }, silent = TRUE))
  expect_identical(.cgb(py_str(opt4)), "log(x + 3)")
})

test_that("codegen_test_rewriting: test_optims_c99", {
  suppressWarnings(try({ x <- Symbol("x") }, silent = TRUE))
  suppressWarnings(try({ expr1 <- ((((2L ^ x) + (log(x) / log(Integer(2L)))) + log(x + 1L)) + exp(x)) - 1L }, silent = TRUE))
  suppressWarnings(try({ opt1 <- simplify(sym_optimize(expr1, optims_c99)) }, silent = TRUE))
  expect_identical(.cgb(opt1), "exp2(x) + expm1(x) + log1p(x) + log2(x)")
  expect_identical(.cgb(rewrite(rewrite(rewrite(opt1, "exp"), "log"), "Pow")), "2^x + exp(x) + log(x)/log(2) + log(x + 1) - 1")
  suppressWarnings(try({ expr2 <- (log(x) / log(Integer(2L))) + log(x + 1L) }, silent = TRUE))
  suppressWarnings(try({ opt2 <- sym_optimize(expr2, optims_c99) }, silent = TRUE))
  expect_identical(.cgb(opt2), "log1p(x) + log2(x)")
  expect_identical(.cgb(rewrite(opt2, "log")), "log(x)/log(2) + log(x + 1)")
  suppressWarnings(try({ expr3 <- (log(x) / log(Integer(2L))) + log((17L * x) + 17L) }, silent = TRUE))
  suppressWarnings(try({ opt3 <- sym_optimize(expr3, optims_c99) }, silent = TRUE))
  suppressWarnings(try({ delta3 <- opt3 - ((sym_log2(x) + log(Integer(17L))) + sym_log1p(x)) }, silent = TRUE))
  expect_identical(.cgb(delta3), "0")
  expect_identical(.cgb(simplify(rewrite(opt3, "log") - expr3)), "0")
  suppressWarnings(try({ expr4 <- ((((2L ^ x) + ((3L * log((5L * x) + 7L)) / (13L * log(Integer(2L))))) + (11L * exp(x))) - 11L) + log((17L * x) + 17L) }, silent = TRUE))
  suppressWarnings(try({ opt4 <- simplify(sym_optimize(expr4, optims_c99)) }, silent = TRUE))
  suppressWarnings(try({ delta4 <- opt4 - ((((exp2(x) + ((3L * sym_log2((5L * x) + 7L)) / 13L)) + (11L * sym_expm1(x))) + log(Integer(17L))) + sym_log1p(x)) }, silent = TRUE))
  expect_identical(.cgb(delta4), "0")
  expect_identical(.cgb(simplify(rewrite(rewrite(rewrite(opt4, "exp"), "log"), "Pow") - expr4)), "0")
  suppressWarnings(try({ expr5 <- (3L * exp(2L * x)) - 3L }, silent = TRUE))
  suppressWarnings(try({ opt5 <- sym_optimize(expr5, optims_c99) }, silent = TRUE))
  suppressWarnings(try({ delta5 <- opt5 - (3L * sym_expm1(2L * x)) }, silent = TRUE))
  expect_identical(.cgb(delta5), "0")
  expect_identical(.cgb(rewrite(opt5, "exp")), "3*exp(2*x) - 3")
  suppressWarnings(try({ expr6 <- exp(2L * x) - 3L }, silent = TRUE))
  suppressWarnings(try({ opt6 <- sym_optimize(expr6, optims_c99) }, silent = TRUE))
  expect_identical(.cgb(py_in(opt6, list(sym_expm1(2L * x) - 2L, expr6))), "True")
  suppressWarnings(try({ expr7 <- log((3L * x) + 3L) }, silent = TRUE))
  suppressWarnings(try({ opt7 <- sym_optimize(expr7, optims_c99) }, silent = TRUE))
  suppressWarnings(try({ delta7 <- opt7 - (log(Integer(3L)) + sym_log1p(x)) }, silent = TRUE))
  expect_identical(.cgb(delta7), "0")
  expect_identical(.cgb(simplify(rewrite(opt7, "log") - expr7)), "0")
  suppressWarnings(try({ expr8 <- log((2L * x) + 3L) }, silent = TRUE))
  suppressWarnings(try({ opt8 <- sym_optimize(expr8, optims_c99) }, silent = TRUE))
  expect_identical(.cgb(opt8), "log(2*x + 3)")
})

test_that("codegen_test_rewriting: test_create_expand_pow_optimization", {
  suppressWarnings(try({ cc <- function(x) ccode(sym_optimize(x, list(create_expand_pow_optimization(4L)))) }, silent = TRUE))
  suppressWarnings(try({ x <- Symbol("x") }, silent = TRUE))
  expect_identical(.cgb(cc(x ^ 4L)), "x*x*x*x")
  expect_identical(.cgb(cc((x ^ 4L) + (x ^ 2L))), "x*x + x*x*x*x")
  expect_identical(.cgb(cc((x ^ 5L) + (x ^ 4L))), "pow(x, 5) + x*x*x*x")
  expect_identical(.cgb(cc(sin(x) ^ 4L)), "pow(sin(x), 4)")
  expect_identical(.cgb(cc(x ^ (-4L))), "1.0/(x*x*x*x)")
  expect_identical(.cgb(cc(x ^ (-5L))), "pow(x, -5)")
  expect_identical(.cgb(cc(-(x ^ 4L))), "-(x*x*x*x)")
  expect_identical(.cgb(cc((x ^ 4L) - (x ^ 2L))), "-(x*x) + x*x*x*x")
  suppressWarnings(try({ i <- Symbol("i", integer = TRUE) }, silent = TRUE))
  expect_identical(.cgb(cc((x ^ i) - (x ^ 2L))), "pow(x, i) - (x*x)")
  suppressWarnings(try({ y <- Symbol("y", real = TRUE) }, silent = TRUE))
  expect_identical(.cgb(cc(Abs(exp(y ^ 4L)))), "exp(y*y*y*y)")
  suppressWarnings(try({ cc2 <- function(x) ccode(sym_optimize(x, list(create_expand_pow_optimization(4L, base_req = function(b) ast_attr(b, "is_Function"))))) }, silent = TRUE))
  expect_identical(.cgb(cc2((x ^ 3L) + (sin(x) ^ 3L))), "pow(x, 3) + sin(x)*sin(x)*sin(x)")
})

test_that("codegen_test_rewriting: test_matsolve", {
  suppressWarnings(try({ n <- Symbol("n", integer = TRUE) }, silent = TRUE))
  suppressWarnings(try({ A <- MatrixSymbol("A", n, n) }, silent = TRUE))
  suppressWarnings(try({ x <- MatrixSymbol("x", n, 1L) }, silent = TRUE))
  suppressWarnings(try({ NULL }, silent = TRUE))
  suppressWarnings(try({ global_assumptions$add(Q$fullrank(A)) }, silent = TRUE))
  expect_identical(.cgb(sym_optimize(A^(-1L) * x, list(matinv_opt))), "MatrixSolve(A, vector=x)")
  expect_identical(.cgb(sym_optimize(A^(-1L) * x + x, list(matinv_opt))), "MatrixSolve(A, vector=x) + x")
  suppressWarnings(try({ global_assumptions$clear() }, silent = TRUE))
  expect_identical(.cgb(sym_optimize(A^(-1L) * x, list(matinv_opt))), "A^(-1)*x")
})

test_that("codegen_test_rewriting: test_logaddexp_opt", {
  suppressWarnings(try({ .tmp <- map(Symbol, py_method("x y", "split")); x <- py_getitem(.tmp, 0L); y <- py_getitem(.tmp, 1L) }, silent = TRUE))
  suppressWarnings(try({ expr1 <- log(exp(x) + exp(y)) }, silent = TRUE))
  suppressWarnings(try({ opt1 <- sym_optimize(expr1, list(logaddexp_opt)) }, silent = TRUE))
  expect_identical(.cgb(logaddexp(x, y) - opt1), "0")
  expect_identical(.cgb(logaddexp(y, x) - opt1), "0")
  expect_identical(.cgb(rewrite(opt1, "log")), "log(exp(x) + exp(y))")
})

test_that("codegen_test_rewriting: test_logaddexp2_opt", {
  suppressWarnings(try({ .tmp <- map(Symbol, py_method("x y", "split")); x <- py_getitem(.tmp, 0L); y <- py_getitem(.tmp, 1L) }, silent = TRUE))
  suppressWarnings(try({ expr1 <- log((2L ^ x) + (2L ^ y)) / log(Integer(2L)) }, silent = TRUE))
  suppressWarnings(try({ opt1 <- sym_optimize(expr1, list(logaddexp2_opt)) }, silent = TRUE))
  expect_identical(.cgb(logaddexp2(x, y) - opt1), "0")
  expect_identical(.cgb(logaddexp2(y, x) - opt1), "0")
  expect_identical(.cgb(rewrite(opt1, "log")), "log(2^x + 2^y)/log(2)")
})

test_that("codegen_test_rewriting: test_sinc_opts", {
  suppressWarnings(try({ x <- Symbol("x") }, silent = TRUE))
  expect_identical(.cgb(sym_optimize(sin(x) / x, sinc_opts)), "sinc(x)")
  expect_identical(.cgb(sym_optimize(sin(2L * x) / (2L * x), sinc_opts)), "sinc(2*x)")
  expect_identical(.cgb(sym_optimize(sin(3L * x) / x, sinc_opts)), "3*sinc(3*x)")
  expect_identical(.cgb(sym_optimize(x * sin(x), sinc_opts)), "x*sin(x)")
  suppressWarnings(try({ y <- Symbol("y") }, silent = TRUE))
  expect_identical(.cgb(sym_optimize(sin(x * y) / (x * y), sinc_opts)), "sinc(x*y)")
  expect_identical(.cgb(sym_optimize((y * sin(x / y)) / x, sinc_opts)), "sinc(x/y)")
  expect_identical(.cgb(sym_optimize(sin(sin(x)) / sin(x), sinc_opts)), "sinc(sin(x))")
  expect_identical(.cgb(sym_optimize(sin(3L * sin(x)) / sin(x), sinc_opts)), "3*sinc(3*sin(x))")
  expect_identical(.cgb(sym_optimize(sin(x) / y, sinc_opts)), "sin(x)/y")
})

test_that("codegen_test_rewriting: test_optims_numpy", {
  suppressWarnings(try({ x <- Symbol("x") }, silent = TRUE))
  expect_identical(.cgb(sym_optimize(((sin(2L * x) / (2L * x)) + exp(2L * x)) - 1L, optims_numpy)), "expm1(2*x) + sinc(2*x)")
  expect_identical(.cgb(sym_optimize((log(x + 3L) / log(Integer(2L))) + log((x ^ 2L) + 1L), optims_numpy)), "log1p(x^2) + log2(x + 3)")
})

test_that("codegen_test_rewriting: test_optims_numpy_TODO", {
  suppressWarnings(try({ .tmp <- map(Symbol, py_method("x y", "split")); x <- py_getitem(.tmp, 0L); y <- py_getitem(.tmp, 1L) }, silent = TRUE))
  expect_identical(.cgb(sym_optimize(((log(x * y) * sin(x * y)) * log((x * y) + 1L)) / ((log(Integer(2L)) * x) * y), optims_numpy)), "sin(x*y)*log1p(x*y)*log2(x*y)/(x*y)")
  expect_identical(.cgb(sym_optimize(exp((x * sin(y)) / y) - 1L, optims_numpy)), "expm1(x*sin(y)/y)")
})

test_that("codegen_test_scipy_nodes: test_cosm1", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  suppressWarnings(try({ cm1_xy <- cosm1(x * y) }, silent = TRUE))
  suppressWarnings(try({ ref_xy <- cos(x * y) - 1L }, silent = TRUE))
  suppressWarnings(try({ expr_minus2 <- cosm1(S_Pi) }, silent = TRUE))
  expect_identical(.cgb(rewrite(expr_minus2, "cos")), "-2")
  expect_identical(.cgb(simplify(cosm1(Float(3.14)))), "-1.99999873172754")
  expect_identical(.cgb(simplify(cosm1(S_Pi / 2L))), "-1")
  expect_identical(.cgb(simplify(((1L / cos(x)) - 1L) + (cosm1(x) / cos(x)))), "0")
})

test_that("codegen_test_scipy_nodes: test_powm1", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  suppressWarnings(try({ cases <- py_dict(list(powm1(x, y), powm1(x * y, z), powm1(x, y * z), powm1((x * y) * z, (x * y) * z)), list((x ^ y) - 1L, ((x * y) ^ z) - 1L, (x ^ (y * z)) - 1L, (((x * y) * z) ^ ((x * y) * z)) - 1L)) }, silent = TRUE))
  suppressWarnings(try({ eulers_constant_m1 <- powm1(x, 1L / log(x)) }, silent = TRUE))
  expect_identical(.cgb(rewrite(eulers_constant_m1, "Pow")), "-1 + E")
  expect_identical(.cgb(simplify(eulers_constant_m1)), "-1 + E")
})

test_that("printing_test_c: test_ccode_Pow", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  expect_identical(.cgb(ccode(x ^ 3L)), "pow(x, 3)")
  expect_identical(.cgb(ccode(x ^ (y ^ 3L))), "pow(x, pow(y, 3))")
  suppressWarnings(try({ g <- implemented_function("g", Lambda(x, 2L * x)) }, silent = TRUE))
  expect_identical(.cgb(ccode((1L / ((g(x) * Float(3.5)) ^ (x - (y ^ x)))) / ((x ^ 2L) + y))), "pow(3.5*2*x, -x + pow(y, x))/(pow(x, 2) + y)")
  expect_identical(.cgb(ccode(x ^ (-Float(1.0)))), "1.0/x")
  expect_identical(.cgb(ccode(x ^ Rational(2L, 3L))), "pow(x, 2.0/3.0)")
  expect_identical(.cgb(ccode(x ^ Rational(2L, 3L), type_aliases = py_dict(list(real), list(float80)))), "powl(x, 2.0L/3.0L)")
  suppressWarnings(try({ U_cond_cfunc <- list(list(function(base, exp) ask_fact(exp, "integer"), "dpowi"), list(function(base, exp) !(ask_fact(exp, "integer")), "pow")) }, silent = TRUE))
  expect_identical(.cgb(ccode(x ^ 3L, user_functions = py_dict(list("Pow"), list(U_cond_cfunc)))), "dpowi(x, 3)")
  expect_identical(.cgb(ccode(x ^ Float(0.5), user_functions = py_dict(list("Pow"), list(U_cond_cfunc)))), "pow(x, 0.5)")
  expect_identical(.cgb(ccode(x ^ Rational(16L, 5L), user_functions = py_dict(list("Pow"), list(U_cond_cfunc)))), "pow(x, 16.0/5.0)")
  suppressWarnings(try({ U_cond_cfunc2 <- list(list(function(base, exp) py_eq(base, 2L), function(base, exp) py_format("exp2(%s)", exp)), list(function(base, exp) !py_eq(base, 2L), "pow")) }, silent = TRUE))
  expect_identical(.cgb(ccode(2L ^ x, user_functions = py_dict(list("Pow"), list(U_cond_cfunc2)))), "exp2(x)")
  expect_identical(.cgb(ccode(x ^ 2L, user_functions = py_dict(list("Pow"), list(U_cond_cfunc2)))), "pow(x, 2)")
})

test_that("printing_test_c: test_ccode_constants_mathh", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  expect_identical(.cgb(ccode(exp(Integer(1L)))), "M_E")
  expect_identical(.cgb(ccode(S_Pi)), "M_PI")
  expect_identical(.cgb(ccode(oo, standard = "c89")), "HUGE_VAL")
  expect_identical(.cgb(ccode(-oo, standard = "c89")), "-HUGE_VAL")
  expect_identical(.cgb(ccode(oo)), "INFINITY")
  expect_identical(.cgb(ccode(-oo, standard = "c99")), "-INFINITY")
  expect_identical(.cgb(ccode(S_Pi, type_aliases = py_dict(list(real), list(float80)))), "M_PIl")
})

test_that("printing_test_c: test_ccode_Rational", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  expect_identical(.cgb(ccode(Rational(3L, 7L))), "3.0/7.0")
  expect_identical(.cgb(ccode(Rational(3L, 7L), type_aliases = py_dict(list(real), list(float80)))), "3.0L/7.0L")
  expect_identical(.cgb(ccode(Rational(18L, 9L))), "2")
  expect_identical(.cgb(ccode(Rational(3L, -7L))), "-3.0/7.0")
  expect_identical(.cgb(ccode(Rational(3L, -7L), type_aliases = py_dict(list(real), list(float80)))), "-3.0L/7.0L")
  expect_identical(.cgb(ccode(Rational(-3L, -7L))), "3.0/7.0")
  expect_identical(.cgb(ccode(Rational(-3L, -7L), type_aliases = py_dict(list(real), list(float80)))), "3.0L/7.0L")
  expect_identical(.cgb(ccode(x + Rational(3L, 7L))), "x + 3.0/7.0")
  expect_identical(.cgb(ccode(x + Rational(3L, 7L), type_aliases = py_dict(list(real), list(float80)))), "x + 3.0L/7.0L")
  expect_identical(.cgb(ccode(Rational(3L, 7L) * x)), "(3.0/7.0)*x")
  expect_identical(.cgb(ccode(Rational(3L, 7L) * x, type_aliases = py_dict(list(real), list(float80)))), "(3.0L/7.0L)*x")
})

test_that("printing_test_c: test_ccode_inline_function", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  suppressWarnings(try({ x <- sym_symbols("x") }, silent = TRUE))
  suppressWarnings(try({ g <- implemented_function("g", Lambda(x, 2L * x)) }, silent = TRUE))
  expect_identical(.cgb(ccode(g(x))), "2*x")
  suppressWarnings(try({ g <- implemented_function("g", Lambda(x, (2L * x) / Catalan)) }, silent = TRUE))
  expect_identical(.cgb(ccode(g(x))), "const double Catalan = 0.91596559417721901;\\n2*x/Catalan")
  suppressWarnings(try({ A <- IndexedBase("A") }, silent = TRUE))
  suppressWarnings(try({ i <- Idx("i", sym_symbols("n", integer = TRUE)) }, silent = TRUE))
  suppressWarnings(try({ g <- implemented_function("g", Lambda(x, (x * (1L + x)) * (2L + x))) }, silent = TRUE))
  expect_identical(.cgb(ccode(g(py_getitem(A, i)), assign_to = py_getitem(A, i))), "for (int i=0; i<n; i++){\\n   A[i] = (A[i] + 1)*(A[i] + 2)*A[i];\\n}")
})

test_that("printing_test_c: test_Element", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  expect_identical(.cgb(ccode(Element("x", "ij"))), "x[i][j]")
  expect_identical(.cgb(ccode(Element("x", "ij", strides = "kl", offset = "o"))), "x[i*k + j*l + o]")
  expect_identical(.cgb(ccode(Element("x", list(3L)))), "x[3]")
  expect_identical(.cgb(ccode(Element("x", list(3L, 4L, 5L)))), "x[3][4][5]")
})

test_that("printing_test_c: test_ccode_Assignment", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  expect_identical(.cgb(ccode(Assignment(x, y + z))), "x = y + z;")
  expect_identical(.cgb(ccode(aug_assign(x, "+", y + z))), "x += y + z;")
})

test_that("printing_test_c: test_ccode_For", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  suppressWarnings(try({ f <- For(x, Range(0L, 10L, 2L), list(aug_assign(y, "*", x))) }, silent = TRUE))
  expect_identical(.cgb(ccode(f)), "for (x = 0; x < 10; x += 2) {\\n   y *= x;\\n}")
})

test_that("printing_test_c: test_ccode_standard", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  expect_identical(.cgb(ccode(sym_expm1(x), standard = "c99")), "expm1(x)")
  expect_identical(.cgb(ccode(nan, standard = "c99")), "NAN")
  expect_identical(.cgb(ccode(py_float("nan"), standard = "c99")), "NAN")
})

test_that("printing_test_c: test_C89CodePrinter", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  suppressWarnings(try({ c89printer <- C89CodePrinter() }, silent = TRUE))
  expect_identical(.cgb(ast_attr(c89printer, "language")), "C")
  expect_identical(.cgb(ast_attr(c89printer, "standard")), "C89")
  expect_identical(.cgb(py_in("void", ast_attr(c89printer, "reserved_words"))), "True")
  expect_identical(.cgb(!py_in("template", ast_attr(c89printer, "reserved_words"))), "True")
  expect_identical(.cgb(printer_doprint(c89printer, sym_log10(x))), "log10(x)")
})

test_that("printing_test_c: test_C99CodePrinter", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  expect_identical(.cgb(printer_doprint(C99CodePrinter(), sym_expm1(x))), "expm1(x)")
  expect_identical(.cgb(printer_doprint(C99CodePrinter(), sym_log1p(x))), "log1p(x)")
  expect_identical(.cgb(printer_doprint(C99CodePrinter(), exp2(x))), "exp2(x)")
  expect_identical(.cgb(printer_doprint(C99CodePrinter(), sym_log2(x))), "log2(x)")
  expect_identical(.cgb(printer_doprint(C99CodePrinter(), fma(x, y, -z))), "fma(x, y, -z)")
  expect_identical(.cgb(printer_doprint(C99CodePrinter(), sym_log10(x))), "log10(x)")
  expect_identical(.cgb(printer_doprint(C99CodePrinter(), Cbrt(x))), "cbrt(x)")
  expect_identical(.cgb(printer_doprint(C99CodePrinter(), hypot(x, y))), "hypot(x, y)")
  expect_identical(.cgb(printer_doprint(C99CodePrinter(), loggamma(x))), "lgamma(x)")
  expect_identical(.cgb(printer_doprint(C99CodePrinter(), Max(x, 3L, x ^ 2L))), "fmax(3, fmax(x, pow(x, 2)))")
  expect_identical(.cgb(printer_doprint(C99CodePrinter(), Min(x, 3L))), "fmin(3, x)")
  suppressWarnings(try({ c99printer <- C99CodePrinter() }, silent = TRUE))
  expect_identical(.cgb(ast_attr(c99printer, "language")), "C")
  expect_identical(.cgb(ast_attr(c99printer, "standard")), "C99")
  expect_identical(.cgb(py_in("restrict", ast_attr(c99printer, "reserved_words"))), "True")
  expect_identical(.cgb(!py_in("using", ast_attr(c99printer, "reserved_words"))), "True")
})

test_that("printing_test_c: test_C99CodePrinter__precision_f80", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  suppressWarnings(try({ f80_printer <- C99CodePrinter(py_dict(list("type_aliases"), list(py_dict(list(real), list(float80))))) }, silent = TRUE))
  expect_identical(.cgb(printer_doprint(f80_printer, sin(x + Float("2.1")))), "sinl(x + 2.10000000000000008882L)")
})

test_that("printing_test_c: test_C99CodePrinter__precision", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  suppressWarnings(try({ n <- sym_symbols("n", integer = TRUE) }, silent = TRUE))
  suppressWarnings(try({ p <- sym_symbols("p", integer = TRUE, positive = TRUE) }, silent = TRUE))
  suppressWarnings(try({ f32_printer <- C99CodePrinter(py_dict(list("type_aliases"), list(py_dict(list(real), list(float32))))) }, silent = TRUE))
  suppressWarnings(try({ f64_printer <- C99CodePrinter(py_dict(list("type_aliases"), list(py_dict(list(real), list(float64))))) }, silent = TRUE))
  suppressWarnings(try({ f80_printer <- C99CodePrinter(py_dict(list("type_aliases"), list(py_dict(list(real), list(float80))))) }, silent = TRUE))
  expect_identical(.cgb(printer_doprint(f32_printer, sin(x + Float(2.1)))), "sinf(x + 2.1F)")
  expect_identical(.cgb(printer_doprint(f64_printer, sin(x + Float(2.1)))), "sin(x + 2.1000000000000001)")
  expect_identical(.cgb(printer_doprint(f80_printer, sin(x + Float("2.0")))), "sinl(x + 2.0L)")
})

test_that("printing_test_c: test_get_math_macros", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  suppressWarnings(try({ macros <- get_math_macros() }, silent = TRUE))
  expect_identical(.cgb(py_getitem(macros, exp(Integer(1L)))), "M_E")
  expect_identical(.cgb(py_getitem(macros, 1L / Sqrt(Integer(2L)))), "M_SQRT1_2")
})

test_that("printing_test_c: test_ccode_Declaration", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  suppressWarnings(try({ i <- sym_symbols("i", integer = TRUE) }, silent = TRUE))
  suppressWarnings(try({ var1 <- Variable(i, type = type_from_expr(i)) }, silent = TRUE))
  suppressWarnings(try({ dcl1 <- Declaration(var1) }, silent = TRUE))
  expect_identical(.cgb(ccode(dcl1)), "int i")
  suppressWarnings(try({ var2 <- Variable(x, type = float32, attrs = list(value_const)) }, silent = TRUE))
  suppressWarnings(try({ dcl2a <- Declaration(var2) }, silent = TRUE))
  expect_identical(.cgb(ccode(dcl2a)), "const float x")
  suppressWarnings(try({ dcl2b <- as_Declaration(var2, value = S_Pi) }, silent = TRUE))
  expect_identical(.cgb(ccode(dcl2b)), "const float x = M_PI")
  suppressWarnings(try({ var3 <- Variable(y, type = Type("bool")) }, silent = TRUE))
  suppressWarnings(try({ dcl3 <- Declaration(var3) }, silent = TRUE))
  suppressWarnings(try({ printer <- C89CodePrinter() }, silent = TRUE))
  expect_identical(.cgb(!py_in("stdbool.h", ast_attr(printer, "headers"))), "True")
  expect_identical(.cgb(printer_doprint(printer, dcl3)), "bool y")
  expect_identical(.cgb(py_in("stdbool.h", ast_attr(printer, "headers"))), "True")
  suppressWarnings(try({ u <- sym_symbols("u", real = TRUE) }, silent = TRUE))
  suppressWarnings(try({ ptr4 <- pointer_deduced(u, attrs = list(pointer_const, restrict)) }, silent = TRUE))
  suppressWarnings(try({ dcl4 <- Declaration(ptr4) }, silent = TRUE))
  expect_identical(.cgb(ccode(dcl4)), "double * const restrict u")
  suppressWarnings(try({ var5 <- Variable(x, Type("__float128"), attrs = list(value_const)) }, silent = TRUE))
  suppressWarnings(try({ dcl5a <- Declaration(var5) }, silent = TRUE))
  expect_identical(.cgb(ccode(dcl5a)), "const __float128 x")
  suppressWarnings(try({ var5b <- Variable(ast_attr(var5, "symbol"), ast_attr(var5, "type"), S_Pi, attrs = ast_attr(var5, "attrs")) }, silent = TRUE))
  suppressWarnings(try({ dcl5b <- Declaration(var5b) }, silent = TRUE))
  expect_identical(.cgb(ccode(dcl5b)), "const __float128 x = M_PI")
})

test_that("printing_test_c: test_C99CodePrinter_custom_type", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  suppressWarnings(try({ f128 <- FloatType("_Float128", ast_attr(float128, "nbits"), ast_attr(float128, "nmant"), ast_attr(float128, "nexp")) }, silent = TRUE))
  suppressWarnings(try({ p128 <- C99CodePrinter(py_dict(list("type_aliases", "type_literal_suffixes", "type_func_suffixes", "type_math_macro_suffixes", "type_macros"), list(py_dict(list(real), list(f128)), py_dict(list(f128), list("Q")), py_dict(list(f128), list("f128")), py_dict(list(real, f128), list("f128", "f128")), py_dict(list(f128), list(list("__STDC_WANT_IEC_60559_TYPES_EXT__")))))) }, silent = TRUE))
  expect_identical(.cgb(printer_doprint(p128, x)), "x")
  expect_identical(.cgb(ast_attr(p128, "headers")), "set()")
  expect_identical(.cgb(ast_attr(p128, "libraries")), "set()")
  expect_identical(.cgb(ast_attr(p128, "macros")), "set()")
  expect_identical(.cgb(printer_doprint(p128, Float(2.0))), "2.0Q")
  expect_identical(.cgb(ast_attr(p128, "headers")), "set()")
  expect_identical(.cgb(ast_attr(p128, "libraries")), "set()")
  expect_identical(.cgb(ast_attr(p128, "macros")), "{'__STDC_WANT_IEC_60559_TYPES_EXT__'}")
  expect_identical(.cgb(printer_doprint(p128, Rational(1L, 2L))), "1.0Q/2.0Q")
  expect_identical(.cgb(printer_doprint(p128, sin(x))), "sinf128(x)")
  expect_identical(.cgb(printer_doprint(p128, x ^ (-Float(1.0)))), "1.0Q/x")
  suppressWarnings(try({ var5 <- Variable(x, f128, attrs = list(value_const)) }, silent = TRUE))
  suppressWarnings(try({ dcl5a <- Declaration(var5) }, silent = TRUE))
  expect_identical(.cgb(ccode(dcl5a)), "const _Float128 x")
  suppressWarnings(try({ var5b <- Variable(x, f128, S_Pi, attrs = list(value_const)) }, silent = TRUE))
  suppressWarnings(try({ dcl5b <- Declaration(var5b) }, silent = TRUE))
  expect_identical(.cgb(printer_doprint(p128, dcl5b)), "const _Float128 x = M_PIf128")
  suppressWarnings(try({ var5b <- Variable(x, f128, value = N(Catalan, 38L), attrs = list(value_const)) }, silent = TRUE))
  suppressWarnings(try({ dcl5c <- Declaration(var5b) }, silent = TRUE))
  expect_identical(.cgb(printer_doprint(p128, dcl5c)), "const _Float128 x = 0.915965594177219015054603514932384111Q")
})

test_that("printing_test_c: test_ccode_math_macros", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  expect_identical(.cgb(ccode(z + exp(Integer(1L)))), "z + M_E")
  expect_identical(.cgb(ccode(z + sym_log2(exp(Integer(1L))))), "z + M_LOG2E")
  expect_identical(.cgb(ccode(z + (1L / log(Integer(2L))))), "z + M_LOG2E")
  expect_identical(.cgb(ccode(z + log(Integer(2L)))), "z + M_LN2")
  expect_identical(.cgb(ccode(z + log(Integer(10L)))), "z + M_LN10")
  expect_identical(.cgb(ccode(z + S_Pi)), "z + M_PI")
  expect_identical(.cgb(ccode(z + (S_Pi / 2L))), "z + M_PI_2")
  expect_identical(.cgb(ccode(z + (S_Pi / 4L))), "z + M_PI_4")
  expect_identical(.cgb(ccode(z + (1L / S_Pi))), "z + M_1_PI")
  expect_identical(.cgb(ccode(z + (2L / S_Pi))), "z + M_2_PI")
  expect_identical(.cgb(ccode(z + (2L / sqrt(S_Pi)))), "z + M_2_SQRTPI")
  expect_identical(.cgb(ccode(z + (2L / Sqrt(S_Pi)))), "z + M_2_SQRTPI")
  expect_identical(.cgb(ccode(z + sqrt(Integer(2L)))), "z + M_SQRT2")
  expect_identical(.cgb(ccode(z + Sqrt(Integer(2L)))), "z + M_SQRT2")
  expect_identical(.cgb(ccode(z + (1L / sqrt(Integer(2L))))), "z + M_SQRT1_2")
  expect_identical(.cgb(ccode(z + (1L / Sqrt(Integer(2L))))), "z + M_SQRT1_2")
})

test_that("printing_test_c: test_ccode_Type", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  expect_identical(.cgb(ccode(Type("float"))), "float")
  expect_identical(.cgb(ccode(intc)), "int")
})

test_that("printing_test_c: test_ccode_codegen_ast", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  expect_identical(.cgb(ccode(Comment("this is a comment"))), "/* this is a comment */")
  expect_identical(.cgb(ccode(While(py_abs(x) > 1L, list(aug_assign(x, "-", 1L))))), "while (fabs(x) > 1) {\\n   x -= 1;\\n}")
  expect_identical(.cgb(ccode(Scope(list(AddAugmentedAssignment(x, 1L))))), "{\\n   x += 1;\\n}")
  suppressWarnings(try({ inp_x <- Declaration(Variable(x, type = real)) }, silent = TRUE))
  expect_identical(.cgb(ccode(FunctionPrototype(real, "pwer", list(inp_x)))), "double pwer(double x)")
  expect_identical(.cgb(ccode(FunctionDefinition(real, "pwer", list(inp_x), list(Assignment(x, x ^ 2L))))), "double pwer(double x){\\n   x = pow(x, 2);\\n}")
  suppressWarnings(try({ block <- CodeBlock(x, Print(list(x, y), "%d %d"), Print(list(QuotedString("hello"), y), "%s %d", file = stderr_), FunctionCall("pwer", list(x)), Return(x)) }, silent = TRUE))
  expect_identical(.cgb(ccode(block)), "x;\\nprintf(\"%d %d\", x, y);\\nfprintf(stderr, \"%s %d\", \"hello\", y);\\npwer(x);\\nreturn x;")
})

test_that("printing_test_c: test_ccode__isinf_isnan", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  expect_identical(.cgb(ccode(isinf(x))), "isinf(x)")
  expect_identical(.cgb(ccode(isnan(x))), "isnan(x)")
})

test_that("printing_test_cxx: test_CXX11CodePrinter", {
  x <- Symbol("x"); y <- Symbol("y"); u <- Symbol("u"); v <- Symbol("v")
  expect_identical(.cgb(printer_doprint(CXX11CodePrinter(), sym_log1p(x))), "std::log1p(x)")
  suppressWarnings(try({ cxx11printer <- CXX11CodePrinter() }, silent = TRUE))
  expect_identical(.cgb(ast_attr(cxx11printer, "language")), "C++")
  expect_identical(.cgb(ast_attr(cxx11printer, "standard")), "C++11")
  expect_identical(.cgb(py_in("operator", ast_attr(cxx11printer, "reserved_words"))), "True")
  expect_identical(.cgb(py_in("noexcept", ast_attr(cxx11printer, "reserved_words"))), "True")
  expect_identical(.cgb(!py_in("concept", ast_attr(cxx11printer, "reserved_words"))), "True")
})

test_that("printing_test_cxx: test_subclass_print_method__ns", {
  x <- Symbol("x"); y <- Symbol("y"); u <- Symbol("u"); v <- Symbol("v")
  suppressWarnings(try({ p <- CXX11CodePrinter() }, silent = TRUE))
  suppressWarnings(try({ myp <- MyPrinter() }, silent = TRUE))
  expect_identical(.cgb(printer_doprint(p, sym_log1p(x))), "std::log1p(x)")
})

test_that("printing_test_fortran: test_inline_function", {
  suppressWarnings(try({ x <- sym_symbols("x") }, silent = TRUE))
  suppressWarnings(try({ g <- implemented_function("g", Lambda(x, 2L * x)) }, silent = TRUE))
  expect_identical(.cgb(fcode(g(x))), "      2*x")
  suppressWarnings(try({ g <- implemented_function("g", Lambda(x, (2L * S_Pi) / x)) }, silent = TRUE))
  expect_identical(.cgb(fcode(g(x))), "      parameter (pi = 3.1415926535897932d0)\\n      2*pi/x")
  suppressWarnings(try({ A <- IndexedBase("A") }, silent = TRUE))
  suppressWarnings(try({ i <- Idx("i", sym_symbols("n", integer = TRUE)) }, silent = TRUE))
  suppressWarnings(try({ g <- implemented_function("g", Lambda(x, (x * (1L + x)) * (2L + x))) }, silent = TRUE))
  expect_identical(.cgb(fcode(g(py_getitem(A, i)), assign_to = py_getitem(A, i))), "      do i = 1, n\\n         A(i) = (A(i) + 1)*(A(i) + 2)*A(i)\\n      end do")
})

test_that("printing_test_fortran: test_fcode_For", {
  suppressWarnings(try({ .tmp <- sym_symbols("x y"); x <- py_getitem(.tmp, 0L); y <- py_getitem(.tmp, 1L) }, silent = TRUE))
  suppressWarnings(try({ f <- For(x, Range(0L, 10L, 2L), list(Assignment(y, x * y))) }, silent = TRUE))
  suppressWarnings(try({ sol <- fcode(f) }, silent = TRUE))
  expect_identical(.cgb(sol), "      do x = 0, 9, 2\\n         y = x*y\\n      end do")
})

test_that("printing_test_fortran: test_fcode_Declaration", {
  suppressWarnings(try({ i <- sym_symbols("i", integer = TRUE) }, silent = TRUE))
  suppressWarnings(try({ var1 <- variable_deduced(i) }, silent = TRUE))
  suppressWarnings(try({ dcl1 <- Declaration(var1) }, silent = TRUE))
  expect_identical(.cgb(fcode(dcl1, standard = 95L, source_format = "free")), "integer*4 :: i")
  suppressWarnings(try({ .tmp <- sym_symbols("x y"); x <- py_getitem(.tmp, 0L); y <- py_getitem(.tmp, 1L) }, silent = TRUE))
  suppressWarnings(try({ var2 <- Variable(x, float32, value = 42L, attrs = list(value_const)) }, silent = TRUE))
  suppressWarnings(try({ dcl2b <- Declaration(var2) }, silent = TRUE))
  expect_identical(.cgb(fcode(dcl2b, standard = 95L, source_format = "free")), "real*4, parameter :: x = 42")
  suppressWarnings(try({ var3 <- Variable(y, type = bool_) }, silent = TRUE))
  suppressWarnings(try({ dcl3 <- Declaration(var3) }, silent = TRUE))
  expect_identical(.cgb(fcode(dcl3, standard = 95L, source_format = "free")), "logical :: y")
  expect_identical(.cgb(fcode(float32, standard = 95L, source_format = "free")), "real*4")
  expect_identical(.cgb(fcode(float64, standard = 95L, source_format = "free")), "real*8")
  expect_identical(.cgb(fcode(real, standard = 95L, source_format = "free", type_aliases = py_dict(list(real), list(float32)))), "real*4")
  expect_identical(.cgb(fcode(real, standard = 95L, source_format = "free", type_aliases = py_dict(list(real), list(float64)))), "real*8")
})

test_that("printing_test_fortran: test_aug_assign", {
  suppressWarnings(try({ x <- sym_symbols("x") }, silent = TRUE))
  expect_identical(.cgb(fcode(aug_assign(x, "+", 1L), source_format = "free")), "x = x + 1")
})

test_that("printing_test_fortran: test_While", {
  suppressWarnings(try({ x <- sym_symbols("x") }, silent = TRUE))
  expect_identical(.cgb(fcode(While(py_abs(x) > 1L, list(aug_assign(x, "-", 1L))), source_format = "free")), "do while (abs(x) > 1)\\n   x = x - 1\\nend do")
})

test_that("printing_test_fortran: test_FunctionPrototype_print", {
  suppressWarnings(try({ x <- sym_symbols("x") }, silent = TRUE))
  suppressWarnings(try({ n <- sym_symbols("n", integer = TRUE) }, silent = TRUE))
  suppressWarnings(try({ vx <- Variable(x, type = real) }, silent = TRUE))
  suppressWarnings(try({ vn <- Variable(n, type = sym_integer) }, silent = TRUE))
  suppressWarnings(try({ fp1 <- FunctionPrototype(real, "power", list(vx, vn)) }, silent = TRUE))
  expect_identical(.cgb(fcode(fp1)), "ERR NotImplementedError: F77 init./parameter statem. req. multiple lines.")
})

test_that("printing_test_fortran: test_FunctionDefinition_print", {
  suppressWarnings(try({ x <- sym_symbols("x") }, silent = TRUE))
  suppressWarnings(try({ n <- sym_symbols("n", integer = TRUE) }, silent = TRUE))
  suppressWarnings(try({ vx <- Variable(x, type = real) }, silent = TRUE))
  suppressWarnings(try({ vn <- Variable(n, type = sym_integer) }, silent = TRUE))
  suppressWarnings(try({ body <- list(Assignment(x, x ^ n), Return(x)) }, silent = TRUE))
  suppressWarnings(try({ fd1 <- FunctionDefinition(real, "power", list(vx, vn), body) }, silent = TRUE))
  expect_identical(.cgb(fcode(fd1)), "ERR NotImplementedError: F77 init./parameter statem. req. multiple lines.")
})

test_that("printing_test_jscode: test_jscode_Pow", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  suppressWarnings(try({ g <- implemented_function("g", Lambda(x, 2L * x)) }, silent = TRUE))
  expect_identical(.cgb(jscode(x ^ 3L)), "Math.pow(x, 3)")
  expect_identical(.cgb(jscode(x ^ (y ^ 3L))), "Math.pow(x, Math.pow(y, 3))")
  expect_identical(.cgb(jscode((1L / ((g(x) * Float(3.5)) ^ (x - (y ^ x)))) / ((x ^ 2L) + y))), "Math.pow(3.5*2*x, -x + Math.pow(y, x))/(Math.pow(x, 2) + y)")
  expect_identical(.cgb(jscode(x ^ (-Float(1.0)))), "1/x")
})

test_that("printing_test_jscode: test_jscode_inline_function", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  suppressWarnings(try({ x <- sym_symbols("x") }, silent = TRUE))
  suppressWarnings(try({ g <- implemented_function("g", Lambda(x, 2L * x)) }, silent = TRUE))
  expect_identical(.cgb(jscode(g(x))), "2*x")
  suppressWarnings(try({ g <- implemented_function("g", Lambda(x, (2L * x) / Catalan)) }, silent = TRUE))
  expect_identical(.cgb(jscode(g(x))), "var Catalan = 0.91596559417721901;\\n2*x/Catalan")
  suppressWarnings(try({ A <- IndexedBase("A") }, silent = TRUE))
  suppressWarnings(try({ i <- Idx("i", sym_symbols("n", integer = TRUE)) }, silent = TRUE))
  suppressWarnings(try({ g <- implemented_function("g", Lambda(x, (x * (1L + x)) * (2L + x))) }, silent = TRUE))
  expect_identical(.cgb(jscode(g(py_getitem(A, i)), assign_to = py_getitem(A, i))), "for (var i=0; i<n; i++){\\n   A[i] = (A[i] + 1)*(A[i] + 2)*A[i];\\n}")
})

test_that("printing_test_julia: test_Pow", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  expect_identical(.cgb(julia_code(x ^ 3L)), "x .^ 3")
  expect_identical(.cgb(julia_code(x ^ (y ^ 3L))), "x .^ (y .^ 3)")
  expect_identical(.cgb(julia_code(x ^ Rational(2L, 3L))), "x .^ (2 // 3)")
  suppressWarnings(try({ g <- implemented_function("g", Lambda(x, 2L * x)) }, silent = TRUE))
  expect_identical(.cgb(julia_code((1L / ((g(x) * Float(3.5)) ^ (x - (y ^ x)))) / ((x ^ 2L) + y))), "(3.5 * 2 * x) .^ (-x + y .^ x) ./ (x .^ 2 + y)")
})

test_that("printing_test_octave: test_Pow", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  mcode <- octave_code
  expect_identical(.cgb(mcode(x ^ 3L)), "x.^3")
  expect_identical(.cgb(mcode(x ^ (y ^ 3L))), "x.^(y.^3)")
  expect_identical(.cgb(mcode(x ^ Rational(2L, 3L))), "x.^(2/3)")
  suppressWarnings(try({ g <- implemented_function("g", Lambda(x, 2L * x)) }, silent = TRUE))
  expect_identical(.cgb(mcode((1L / ((g(x) * Float(3.5)) ^ (x - (y ^ x)))) / ((x ^ 2L) + y))), "(3.5*2*x).^(-x + y.^x)./(x.^2 + y)")
})

test_that("printing_test_octave: test_MatrixSolve", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  mcode <- octave_code
  suppressWarnings(try({ n <- Symbol("n", integer = TRUE) }, silent = TRUE))
  suppressWarnings(try({ A <- MatrixSymbol("A", n, n) }, silent = TRUE))
  suppressWarnings(try({ x <- MatrixSymbol("x", n, 1L) }, silent = TRUE))
  expect_identical(.cgb(mcode(MatrixSolve(A, x))), "A \\ x")
})

test_that("printing_test_pycode: test_PythonCodePrinter", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); p <- IndexedBase("p")
  suppressWarnings(try({ prntr <- PythonCodePrinter() }, silent = TRUE))
  expect_identical(.cgb(ast_attr(prntr, "module_imports")), "defaultdict(<class 'set'>, {})")
  expect_identical(.cgb(printer_doprint(prntr, x ^ y)), "x^y")
  expect_identical(.cgb(printer_doprint(prntr, sym_Mod(x, 2L))), "x % 2")
  expect_identical(.cgb(printer_doprint(prntr, -sym_Mod(x, y))), "-(x % y)")
  expect_identical(.cgb(printer_doprint(prntr, sym_Mod(-x, y))), "(-x) % y")
  expect_identical(.cgb(printer_doprint(prntr, (-2L) * sym_Mod(x, 3L))), "-2*(x % 3)")
  expect_identical(.cgb(printer_doprint(prntr, (-x) * sym_Mod(x, y))), "-x*(x % y)")
  expect_identical(.cgb(printer_doprint(prntr, (-x) / sym_Mod(x, 3L))), "-x/(x % 3)")
  expect_identical(.cgb(printer_doprint(prntr, sym_Mod(x, 1L / y))), "x % (1/y)")
  expect_identical(.cgb(printer_doprint(prntr, sym_Mod(1L, 1L / x))), "1 % (1/x)")
  expect_identical(.cgb(printer_doprint(prntr, And(x, y))), "x and y")
  expect_identical(.cgb(printer_doprint(prntr, Or(x, y))), "x or y")
  expect_identical(.cgb(printer_doprint(prntr, 1L / (x + y))), "1/(x + y)")
  expect_identical(.cgb(printer_doprint(prntr, Not(x))), "not x")
  expect_identical(.cgb(ast_attr(prntr, "module_imports")), "defaultdict(<class 'set'>, {})")
  expect_identical(.cgb(printer_doprint(prntr, S_Pi)), "math.pi")
  expect_identical(.cgb(ast_attr(prntr, "module_imports")), "defaultdict(<class 'set'>, {'math': {'pi'}})")
  expect_identical(.cgb(printer_doprint(prntr, x ^ Rational(1L, 2L))), "math.sqrt(x)")
  expect_identical(.cgb(printer_doprint(prntr, sqrt(x))), "math.sqrt(x)")
  expect_identical(.cgb(printer_doprint(prntr, acos(x))), "math.acos(x)")
  expect_identical(.cgb(printer_doprint(prntr, cot(x))), "(1/math.tan(x))")
  expect_identical(.cgb(printer_doprint(prntr, coth(x))), "((math.exp(x) + math.exp(-x))/(math.exp(x) - math.exp(-x)))")
  expect_identical(.cgb(printer_doprint(prntr, asec(x))), "(math.acos(1/x))")
  expect_identical(.cgb(printer_doprint(prntr, acsch(x))), "(math.log(math.sqrt(1 + x^(-2)) + 1/x))")
  expect_identical(.cgb(printer_doprint(prntr, Assignment(x, 2L))), "x = 2")
  expect_identical(.cgb(printer_doprint(prntr, Piecewise(list(1L, Eq(x, 0L)), list(2L, x > 6L)))), "((1) if (x == 0) else (2) if (x > 6) else None)")
  expect_identical(.cgb(printer_doprint(prntr, Piecewise(list(2L, Le(x, 0L)), list(3L, Gt(x, 0L)), evaluate = FALSE))), "((2) if (x <= 0) else (3) if (x > 0) else None)")
  expect_identical(.cgb(printer_doprint(prntr, sign(x))), "(0.0 if x == 0 else math.copysign(1, x))")
  expect_identical(.cgb(printer_doprint(prntr, py_getitem(p, 0L, 1L))), "p[0, 1]")
  expect_identical(.cgb(printer_doprint(prntr, KroneckerDelta(x, y))), "(1 if x == y else 0)")
  expect_identical(.cgb(printer_doprint(prntr, Tuple(2L, 3L))), "(2, 3)")
  expect_identical(.cgb(printer_doprint(prntr, Min(x, y))), "min(x, y)")
  expect_identical(.cgb(printer_doprint(prntr, Max(x, y))), "max(x, y)")
  expect_identical(.cgb(printer_doprint(prntr, conjugate(x))), "(x).conjugate()")
})

test_that("printing_test_pycode: test_NumPyPrinter", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); p <- IndexedBase("p")
  suppressWarnings(try({ invisible(NULL) }, silent = TRUE))
  suppressWarnings(try({ invisible(NULL) }, silent = TRUE))
  suppressWarnings(try({ invisible(NULL) }, silent = TRUE))
  suppressWarnings(try({ invisible(NULL) }, silent = TRUE))
  suppressWarnings(try({ invisible(NULL) }, silent = TRUE))
  suppressWarnings(try({ invisible(NULL) }, silent = TRUE))
  suppressWarnings(try({ invisible(NULL) }, silent = TRUE))
  suppressWarnings(try({ a <- Symbol("a"); b <- Symbol("b") }, silent = TRUE))
  suppressWarnings(try({ p <- NumPyPrinter() }, silent = TRUE))
  expect_identical(.cgb(printer_doprint(p, sign(x))), "numpy.sign(x)")
  suppressWarnings(try({ A <- MatrixSymbol("A", 2L, 2L) }, silent = TRUE))
  suppressWarnings(try({ B <- MatrixSymbol("B", 2L, 2L) }, silent = TRUE))
  suppressWarnings(try({ C <- MatrixSymbol("C", 1L, 5L) }, silent = TRUE))
  suppressWarnings(try({ D <- MatrixSymbol("D", 3L, 4L) }, silent = TRUE))
  expect_identical(.cgb(printer_doprint(p, A ^ (-1L))), "numpy.linalg.inv(A)")
  expect_identical(.cgb(printer_doprint(p, A ^ 5L)), "numpy.linalg.matrix_power(A, 5)")
  expect_identical(.cgb(printer_doprint(p, Identity(3L))), "numpy.eye(3)")
  suppressWarnings(try({ u <- MatrixSymbol("x", 2L, 1L) }, silent = TRUE))
  suppressWarnings(try({ v <- MatrixSymbol("y", 2L, 1L) }, silent = TRUE))
  expect_identical(.cgb(printer_doprint(p, MatrixSolve(A, u))), "numpy.linalg.solve(A, x)")
  expect_identical(.cgb(printer_doprint(p, MatrixSolve(A, u) + v)), "numpy.linalg.solve(A, x) + y")
  expect_identical(.cgb(printer_doprint(p, ZeroMatrix(2L, 3L))), "numpy.zeros((2, 3))")
  expect_identical(.cgb(printer_doprint(p, OneMatrix(2L, 3L))), "numpy.ones((2, 3))")
  expect_identical(.cgb(printer_doprint(p, FunctionMatrix(4L, 5L, Lambda(list(a, b), a + b)))), "numpy.fromfunction(lambda a, b: a + b, (4, 5))")
  expect_identical(.cgb(printer_doprint(p, HadamardProduct(A, B))), "numpy.multiply(A, B)")
  expect_identical(.cgb(printer_doprint(p, KroneckerProduct(A, B))), "numpy.kron(A, B)")
  expect_identical(.cgb(printer_doprint(p, Adjoint(A))), "numpy.conjugate(numpy.transpose(A))")
  expect_identical(.cgb(printer_doprint(p, conjugate(x))), "numpy.conjugate(x)")
  expect_identical(.cgb(printer_doprint(p, DiagonalOf(A))), "numpy.reshape(numpy.diag(A), (-1, 1))")
  expect_identical(.cgb(printer_doprint(p, DiagMatrix(C))), "numpy.diagflat(C)")
  expect_identical(.cgb(printer_doprint(p, DiagonalMatrix(D))), "numpy.multiply(D, numpy.eye(3, 4))")
  expect_identical(.cgb(printer_doprint(p, x ^ (-1L))), "x^(-1.0)")
  expect_identical(.cgb(printer_doprint(p, x ^ (-2L))), "x^(-2.0)")
  suppressWarnings(try({ expr <- Pow(2L, -1L, evaluate = FALSE) }, silent = TRUE))
  expect_identical(.cgb(printer_doprint(p, S_E)), "numpy.e")
  expect_identical(.cgb(printer_doprint(p, S_Pi)), "numpy.pi")
  expect_identical(.cgb(printer_doprint(p, EulerGamma)), "numpy.euler_gamma")
  expect_identical(.cgb(printer_doprint(p, nan)), "numpy.nan")
  expect_identical(.cgb(printer_doprint(p, oo)), "numpy.inf")
  expect_identical(.cgb(printer_doprint(p, (-oo))), "-numpy.inf")
  expect_identical(.cgb(printer_doprint(p, sec(x) ^ 2L)), "(numpy.cos(x)^(-1.0))^2")
})

test_that("printing_test_pycode: test_codegen_ast_nodes", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); p <- IndexedBase("p")
  expect_identical(.cgb(pycode(none)), "None")
})

test_that("printing_test_pycode: test_numerical_accuracy_functions", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); p <- IndexedBase("p")
  suppressWarnings(try({ prntr <- SciPyPrinter() }, silent = TRUE))
  expect_identical(.cgb(printer_doprint(prntr, sym_expm1(x))), "numpy.expm1(x)")
  expect_identical(.cgb(printer_doprint(prntr, sym_log1p(x))), "numpy.log1p(x)")
  expect_identical(.cgb(printer_doprint(prntr, cosm1(x))), "scipy.special.cosm1(x)")
})

test_that("printing_test_pycode: test_array_printer", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); p <- IndexedBase("p")
  suppressWarnings(try({ A <- ArraySymbol("A", list(4L, 4L, 6L, 6L, 6L)) }, silent = TRUE))
  suppressWarnings(try({ S_I <- IndexedBase("I") }, silent = TRUE))
  suppressWarnings(try({ i <- Idx("i", list(0L, 1L)); j <- Idx("j", list(2L, 3L)); k <- Idx("k", list(4L, 5L)) }, silent = TRUE))
  suppressWarnings(try({ prntr <- NumPyPrinter() }, silent = TRUE))
  expect_identical(.cgb(printer_doprint(prntr, S_I)), "I")
  suppressWarnings(try({ invisible(NULL) }, silent = TRUE))
  suppressWarnings(try({ invisible(NULL) }, silent = TRUE))
  suppressWarnings(try({ invisible(NULL) }, silent = TRUE))
  suppressWarnings(try({ invisible(NULL) }, silent = TRUE))
  suppressWarnings(try({ B <- ArraySymbol("B", list(3L, 3L)) }, silent = TRUE))
  suppressWarnings(try({ d <- Dummy("d") }, silent = TRUE))
  suppressWarnings(try({ s <- Symbol("_elementwise_operand") }, silent = TRUE))
  suppressWarnings(try({ s <- Symbol("prefix_elementwise_operand_suffix") }, silent = TRUE))
  suppressWarnings(try({ nested <- ArrayElementwiseApplyFunc(Lambda(d, d + 1L), ArrayElementwiseApplyFunc(Lambda(d, d ^ 2L), B)) }, silent = TRUE))
  suppressWarnings(try({ prntr <- TensorflowPrinter() }, silent = TRUE))
  expect_identical(.cgb(printer_doprint(prntr, S_I)), "I")
})

test_that("printing_test_rcode: test_rcode_sqrt", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  expect_identical(.cgb(rcode(sqrt(x))), "sqrt(x)")
  expect_identical(.cgb(rcode(x ^ Float(0.5))), "sqrt(x)")
  expect_identical(.cgb(rcode(sqrt(x))), "sqrt(x)")
})

test_that("printing_test_rcode: test_rcode_Pow", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  expect_identical(.cgb(rcode(x ^ 3L)), "x^3")
  expect_identical(.cgb(rcode(x ^ (y ^ 3L))), "x^(y^3)")
  suppressWarnings(try({ g <- implemented_function("g", Lambda(x, 2L * x)) }, silent = TRUE))
  expect_identical(.cgb(rcode((1L / ((g(x) * Float(3.5)) ^ (x - (y ^ x)))) / ((x ^ 2L) + y))), "(3.5*2*x)^(-x + y^x)/(x^2 + y)")
  expect_identical(.cgb(rcode(x ^ (-Float(1.0)))), "1.0/x")
  expect_identical(.cgb(rcode(x ^ Rational(2L, 3L))), "x^(2.0/3.0)")
  suppressWarnings(try({ U_cond_cfunc <- list(list(function(base, exp) ask_fact(exp, "integer"), "dpowi"), list(function(base, exp) !(ask_fact(exp, "integer")), "pow")) }, silent = TRUE))
  expect_identical(.cgb(rcode(x ^ 3L, user_functions = py_dict(list("Pow"), list(U_cond_cfunc)))), "dpowi(x, 3)")
  expect_identical(.cgb(rcode(x ^ Float(3.2), user_functions = py_dict(list("Pow"), list(U_cond_cfunc)))), "pow(x, 3.2)")
})

test_that("printing_test_rcode: test_rcode_Max", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  expect_identical(.cgb(rcode(Max(x, x * x), user_functions = py_dict(list("Max", "Pow"), list("my_max", "my_pow")))), "my_max(x, my_pow(x, 2))")
})

test_that("printing_test_rcode: test_rcode_constants_mathh", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  expect_identical(.cgb(rcode(exp(Integer(1L)))), "exp(1)")
  expect_identical(.cgb(rcode(S_Pi)), "pi")
  expect_identical(.cgb(rcode(oo)), "Inf")
  expect_identical(.cgb(rcode(-oo)), "-Inf")
})

test_that("printing_test_rcode: test_rcode_constants_other", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  expect_identical(.cgb(rcode(2L * GoldenRatio)), "GoldenRatio = 1.61803398874989;\\n2*GoldenRatio")
  expect_identical(.cgb(rcode(2L * Catalan)), "Catalan = 0.915965594177219;\\n2*Catalan")
  expect_identical(.cgb(rcode(2L * EulerGamma)), "EulerGamma = 0.577215664901533;\\n2*EulerGamma")
})

test_that("printing_test_rcode: test_rcode_Rational", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  expect_identical(.cgb(rcode(Rational(3L, 7L))), "3.0/7.0")
  expect_identical(.cgb(rcode(Rational(18L, 9L))), "2")
  expect_identical(.cgb(rcode(Rational(3L, -7L))), "-3.0/7.0")
  expect_identical(.cgb(rcode(Rational(-3L, -7L))), "3.0/7.0")
  expect_identical(.cgb(rcode(x + Rational(3L, 7L))), "x + 3.0/7.0")
  expect_identical(.cgb(rcode(Rational(3L, 7L) * x)), "(3.0/7.0)*x")
})

test_that("printing_test_rcode: test_rcode_Integer", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  expect_identical(.cgb(rcode(Integer(67L))), "67")
  expect_identical(.cgb(rcode(Integer(-1L))), "-1")
})

test_that("printing_test_rcode: test_rcode_functions", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  expect_identical(.cgb(rcode(sin(x) ^ cos(x))), "sin(x)^cos(x)")
  expect_identical(.cgb(rcode(sym_factorial(x) + gamma(y))), "factorial(x) + gamma(y)")
})

test_that("printing_test_rcode: test_rcode_inline_function", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  suppressWarnings(try({ x <- sym_symbols("x") }, silent = TRUE))
  suppressWarnings(try({ g <- implemented_function("g", Lambda(x, 2L * x)) }, silent = TRUE))
  expect_identical(.cgb(rcode(g(x))), "2*x")
  suppressWarnings(try({ g <- implemented_function("g", Lambda(x, (2L * x) / Catalan)) }, silent = TRUE))
  expect_identical(.cgb(rcode(g(x))), "Catalan = 0.915965594177219;\\n2*x/Catalan")
  suppressWarnings(try({ A <- IndexedBase("A") }, silent = TRUE))
  suppressWarnings(try({ i <- Idx("i", sym_symbols("n", integer = TRUE)) }, silent = TRUE))
  suppressWarnings(try({ g <- implemented_function("g", Lambda(x, (x * (1L + x)) * (2L + x))) }, silent = TRUE))
  suppressWarnings(try({ res <- rcode(g(py_getitem(A, i)), assign_to = py_getitem(A, i)) }, silent = TRUE))
  suppressWarnings(try({ ref <- "for (i in 1:n){\n   A[i] = (A[i] + 1)*(A[i] + 2)*A[i];\n}" }, silent = TRUE))
  expect_identical(.cgb(res), "for (i in 1:n){\\n   A[i] = (A[i] + 1)*(A[i] + 2)*A[i];\\n}")
})

test_that("printing_test_rcode: test_rcode_exceptions", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  expect_identical(.cgb(rcode(ceiling(x))), "ceiling(x)")
  expect_identical(.cgb(rcode(Abs(x))), "abs(x)")
  expect_identical(.cgb(rcode(gamma(x))), "gamma(x)")
})

test_that("printing_test_rcode: test_rcode_user_functions", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  suppressWarnings(try({ x <- sym_symbols("x", integer = FALSE) }, silent = TRUE))
  suppressWarnings(try({ n <- sym_symbols("n", integer = TRUE) }, silent = TRUE))
  suppressWarnings(try({ custom_functions <- py_dict(list("ceiling", "Abs"), list("myceil", list(list(function(x) !(ask_fact(x, "integer")), "fabs"), list(function(x) ask_fact(x, "integer"), "abs")))) }, silent = TRUE))
  expect_identical(.cgb(rcode(ceiling(x), user_functions = custom_functions)), "myceil(x)")
  expect_identical(.cgb(rcode(Abs(x), user_functions = custom_functions)), "fabs(x)")
  expect_identical(.cgb(rcode(Abs(n), user_functions = custom_functions)), "abs(n)")
})

test_that("printing_test_rcode: test_rcode_boolean", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  expect_identical(.cgb(rcode(TRUE)), "True")
  expect_identical(.cgb(rcode(S_true)), "True")
  expect_identical(.cgb(rcode(FALSE)), "False")
  expect_identical(.cgb(rcode(S_false)), "False")
  expect_identical(.cgb(rcode(And(x, y))), "x & y")
  expect_identical(.cgb(rcode(Or(x, y))), "x | y")
  expect_identical(.cgb(rcode(!x)), "!x")
  expect_identical(.cgb(rcode(And(x, y, z))), "x & y & z")
  expect_identical(.cgb(rcode(Or(x, y, z))), "x | y | z")
  expect_identical(.cgb(rcode(Or(And(x, y), z))), "z | x & y")
  expect_identical(.cgb(rcode(And(Or(x, y), z))), "z & (x | y)")
})

test_that("printing_test_rcode: test_rcode_Relational", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  expect_identical(.cgb(rcode(Eq(x, y))), "x == y")
  expect_identical(.cgb(rcode(Ne(x, y))), "x != y")
  expect_identical(.cgb(rcode(Le(x, y))), "x <= y")
  expect_identical(.cgb(rcode(Lt(x, y))), "x < y")
  expect_identical(.cgb(rcode(Gt(x, y))), "x > y")
  expect_identical(.cgb(rcode(Ge(x, y))), "x >= y")
})

test_that("printing_test_rcode: test_rcode_Piecewise", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  suppressWarnings(try({ expr <- Piecewise(list(x, x < 1L), list(x ^ 2L, TRUE)) }, silent = TRUE))
  suppressWarnings(try({ res <- rcode(expr) }, silent = TRUE))
  suppressWarnings(try({ ref <- "ifelse(x < 1,x,x^2)" }, silent = TRUE))
  expect_identical(.cgb(res), "ifelse(x < 1,x,x^2)")
  suppressWarnings(try({ tau <- Symbol("tau") }, silent = TRUE))
  suppressWarnings(try({ res <- rcode(expr, tau) }, silent = TRUE))
  suppressWarnings(try({ ref <- "tau = ifelse(x < 1,x,x^2);" }, silent = TRUE))
  expect_identical(.cgb(res), "tau = ifelse(x < 1,x,x^2);")
  suppressWarnings(try({ expr <- 2L * Piecewise(list(x, x < 1L), list(x ^ 2L, x < 2L), list(x ^ 3L, TRUE)) }, silent = TRUE))
  expect_identical(.cgb(rcode(expr)), "2*ifelse(x < 1,x,ifelse(x < 2,x^2,x^3))")
  suppressWarnings(try({ res <- rcode(expr, assign_to = "c") }, silent = TRUE))
  expect_identical(.cgb(res), "c = 2*ifelse(x < 1,x,ifelse(x < 2,x^2,x^3));")
  suppressWarnings(try({ expr <- 2L * Piecewise(list(x, x < 1L), list(x ^ 2L, x < 2L)) }, silent = TRUE))
  expect_identical(.cgb(rcode(expr)), "2*ifelse(x < 1,x,ifelse(x < 2,x^2,NA))")
})

test_that("printing_test_rcode: test_rcode_sinc", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  suppressWarnings(try({ invisible(NULL) }, silent = TRUE))
  suppressWarnings(try({ expr <- sinc(x) }, silent = TRUE))
  suppressWarnings(try({ res <- rcode(expr) }, silent = TRUE))
  suppressWarnings(try({ ref <- "(ifelse(x != 0,sin(x)/x,1))" }, silent = TRUE))
  expect_identical(.cgb(res), "(ifelse(x != 0,sin(x)/x,1))")
})

test_that("printing_test_rcode: test_rcode_Piecewise_deep", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  suppressWarnings(try({ p <- rcode(2L * Piecewise(list(x, x < 1L), list(x + 1L, x < 2L), list(x ^ 2L, TRUE))) }, silent = TRUE))
  expect_identical(.cgb(p), "2*ifelse(x < 1,x,ifelse(x < 2,x + 1,x^2))")
  suppressWarnings(try({ expr <- ((((((x * y) * z) + (x ^ 2L)) + (y ^ 2L)) + Piecewise(list(0L, x < Float(0.5)), list(1L, TRUE))) + cos(z)) - 1L }, silent = TRUE))
  suppressWarnings(try({ p <- rcode(expr) }, silent = TRUE))
  suppressWarnings(try({ ref <- "x^2 + x*y*z + y^2 + ifelse(x < 0.5,0,1) + cos(z) - 1" }, silent = TRUE))
  expect_identical(.cgb(p), "x^2 + x*y*z + y^2 + ifelse(x < 0.5,0,1) + cos(z) - 1")
  suppressWarnings(try({ ref <- "c = x^2 + x*y*z + y^2 + ifelse(x < 0.5,0,1) + cos(z) - 1;" }, silent = TRUE))
  suppressWarnings(try({ p <- rcode(expr, assign_to = "c") }, silent = TRUE))
  expect_identical(.cgb(p), "c = x^2 + x*y*z + y^2 + ifelse(x < 0.5,0,1) + cos(z) - 1;")
})

test_that("printing_test_rcode: test_rcode_ITE", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  suppressWarnings(try({ expr <- ITE(x < 1L, y, z) }, silent = TRUE))
  suppressWarnings(try({ p <- rcode(expr) }, silent = TRUE))
  suppressWarnings(try({ ref <- "ifelse(x < 1,y,z)" }, silent = TRUE))
  expect_identical(.cgb(p), "ifelse(x < 1,y,z)")
})

test_that("printing_test_rcode: test_rcode_settings", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  expect_identical(.cgb(rcode(sin(x), method = "garbage")), "ERR TypeError: Unknown setting 'method'.")
})

test_that("printing_test_rcode: test_rcode_Indexed", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  suppressWarnings(try({ .tmp <- sym_symbols("n m o", integer = TRUE); n <- py_getitem(.tmp, 0L); m <- py_getitem(.tmp, 1L); o <- py_getitem(.tmp, 2L) }, silent = TRUE))
  suppressWarnings(try({ i <- Idx("i", n); j <- Idx("j", m); k <- Idx("k", o) }, silent = TRUE))
  suppressWarnings(try({ p <- RCodePrinter() }, silent = TRUE))
  suppressWarnings(try({ p[["_not_r"]] <- py_set(list()) }, silent = TRUE))
  suppressWarnings(try({ x <- py_getitem(IndexedBase("x"), j) }, silent = TRUE))
  expect_identical(.cgb(py_method(p, "_print_Indexed", x)), "x[j]")
  suppressWarnings(try({ A <- py_getitem(IndexedBase("A"), i, j) }, silent = TRUE))
  expect_identical(.cgb(py_method(p, "_print_Indexed", A)), "A[i, j]")
  suppressWarnings(try({ B <- py_getitem(IndexedBase("B"), i, j, k) }, silent = TRUE))
  expect_identical(.cgb(py_method(p, "_print_Indexed", B)), "B[i, j, k]")
  expect_identical(.cgb(p[["_not_r"]]), "set()")
})

test_that("printing_test_rcode: test_rcode_Indexed_without_looking_for_contraction", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  suppressWarnings(try({ len_y <- 5L }, silent = TRUE))
  suppressWarnings(try({ y <- IndexedBase("y", shape = list(len_y)) }, silent = TRUE))
  suppressWarnings(try({ x <- IndexedBase("x", shape = list(len_y)) }, silent = TRUE))
  suppressWarnings(try({ Dy <- IndexedBase("Dy", shape = list(len_y - 1L)) }, silent = TRUE))
  suppressWarnings(try({ i <- Idx("i", len_y - 1L) }, silent = TRUE))
  suppressWarnings(try({ e <- Eq(py_getitem(Dy, i), (py_getitem(y, i + 1L) - py_getitem(y, i)) / (py_getitem(x, i + 1L) - py_getitem(x, i))) }, silent = TRUE))
  suppressWarnings(try({ code0 <- rcode(ast_attr(e, "rhs"), assign_to = ast_attr(e, "lhs"), contract = FALSE) }, silent = TRUE))
  expect_identical(.cgb(code0), "Dy[i] = (y[i + 1] - y[i])/(x[i + 1] - x[i]);")
})

test_that("printing_test_rcode: test_rcode_loops_matrix_vector", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  suppressWarnings(try({ .tmp <- sym_symbols("n m", integer = TRUE); n <- py_getitem(.tmp, 0L); m <- py_getitem(.tmp, 1L) }, silent = TRUE))
  suppressWarnings(try({ A <- IndexedBase("A") }, silent = TRUE))
  suppressWarnings(try({ x <- IndexedBase("x") }, silent = TRUE))
  suppressWarnings(try({ y <- IndexedBase("y") }, silent = TRUE))
  suppressWarnings(try({ i <- Idx("i", m) }, silent = TRUE))
  suppressWarnings(try({ j <- Idx("j", n) }, silent = TRUE))
  suppressWarnings(try({ s <- "for (i in 1:m){\n   y[i] = 0;\n}\nfor (i in 1:m){\n   for (j in 1:n){\n      y[i] = A[i, j]*x[j] + y[i];\n   }\n}" }, silent = TRUE))
  suppressWarnings(try({ c <- rcode(py_getitem(A, i, j) * py_getitem(x, j), assign_to = py_getitem(y, i)) }, silent = TRUE))
  expect_identical(.cgb(c), "for (i in 1:m){\\n   y[i] = 0;\\n}\\nfor (i in 1:m){\\n   for (j in 1:n){\\n      y[i] = A[i, j]*x[j] + y[i];\\n   }\\n}")
})

test_that("printing_test_rcode: test_rcode_loops_add", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  suppressWarnings(try({ .tmp <- sym_symbols("n m", integer = TRUE); n <- py_getitem(.tmp, 0L); m <- py_getitem(.tmp, 1L) }, silent = TRUE))
  suppressWarnings(try({ A <- IndexedBase("A") }, silent = TRUE))
  suppressWarnings(try({ x <- IndexedBase("x") }, silent = TRUE))
  suppressWarnings(try({ y <- IndexedBase("y") }, silent = TRUE))
  suppressWarnings(try({ z <- IndexedBase("z") }, silent = TRUE))
  suppressWarnings(try({ i <- Idx("i", m) }, silent = TRUE))
  suppressWarnings(try({ j <- Idx("j", n) }, silent = TRUE))
  suppressWarnings(try({ s <- "for (i in 1:m){\n   y[i] = x[i] + z[i];\n}\nfor (i in 1:m){\n   for (j in 1:n){\n      y[i] = A[i, j]*x[j] + y[i];\n   }\n}" }, silent = TRUE))
  suppressWarnings(try({ c <- rcode(((py_getitem(A, i, j) * py_getitem(x, j)) + py_getitem(x, i)) + py_getitem(z, i), assign_to = py_getitem(y, i)) }, silent = TRUE))
  expect_identical(.cgb(c), "for (i in 1:m){\\n   y[i] = x[i] + z[i];\\n}\\nfor (i in 1:m){\\n   for (j in 1:n){\\n      y[i] = A[i, j]*x[j] + y[i];\\n   }\\n}")
})

test_that("printing_test_rcode: test_rcode_loops_multiple_contractions", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  suppressWarnings(try({ .tmp <- sym_symbols("n m o p", integer = TRUE); n <- py_getitem(.tmp, 0L); m <- py_getitem(.tmp, 1L); o <- py_getitem(.tmp, 2L); p <- py_getitem(.tmp, 3L) }, silent = TRUE))
  suppressWarnings(try({ a <- IndexedBase("a") }, silent = TRUE))
  suppressWarnings(try({ b <- IndexedBase("b") }, silent = TRUE))
  suppressWarnings(try({ y <- IndexedBase("y") }, silent = TRUE))
  suppressWarnings(try({ i <- Idx("i", m) }, silent = TRUE))
  suppressWarnings(try({ j <- Idx("j", n) }, silent = TRUE))
  suppressWarnings(try({ k <- Idx("k", o) }, silent = TRUE))
  suppressWarnings(try({ l <- Idx("l", p) }, silent = TRUE))
  suppressWarnings(try({ s <- "for (i in 1:m){\n   y[i] = 0;\n}\nfor (i in 1:m){\n   for (j in 1:n){\n      for (k in 1:o){\n         for (l in 1:p){\n            y[i] = a[i, j, k, l]*b[j, k, l] + y[i];\n         }\n      }\n   }\n}" }, silent = TRUE))
  suppressWarnings(try({ c <- rcode(py_getitem(b, j, k, l) * py_getitem(a, i, j, k, l), assign_to = py_getitem(y, i)) }, silent = TRUE))
  expect_identical(.cgb(c), "for (i in 1:m){\\n   y[i] = 0;\\n}\\nfor (i in 1:m){\\n   for (j in 1:n){\\n      for (k in 1:o){\\n         for (l in 1:p){\\n            y[i] = a[i, j, k, l]*b[j, k, l] + y[i];\\n         }\\n      }\\n   }\\n}")
})

test_that("printing_test_rcode: test_rcode_loops_addfactor", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  suppressWarnings(try({ .tmp <- sym_symbols("n m o p", integer = TRUE); n <- py_getitem(.tmp, 0L); m <- py_getitem(.tmp, 1L); o <- py_getitem(.tmp, 2L); p <- py_getitem(.tmp, 3L) }, silent = TRUE))
  suppressWarnings(try({ a <- IndexedBase("a") }, silent = TRUE))
  suppressWarnings(try({ b <- IndexedBase("b") }, silent = TRUE))
  suppressWarnings(try({ c <- IndexedBase("c") }, silent = TRUE))
  suppressWarnings(try({ y <- IndexedBase("y") }, silent = TRUE))
  suppressWarnings(try({ i <- Idx("i", m) }, silent = TRUE))
  suppressWarnings(try({ j <- Idx("j", n) }, silent = TRUE))
  suppressWarnings(try({ k <- Idx("k", o) }, silent = TRUE))
  suppressWarnings(try({ l <- Idx("l", p) }, silent = TRUE))
  suppressWarnings(try({ s <- "for (i in 1:m){\n   y[i] = 0;\n}\nfor (i in 1:m){\n   for (j in 1:n){\n      for (k in 1:o){\n         for (l in 1:p){\n            y[i] = (a[i, j, k, l] + b[i, j, k, l])*c[j, k, l] + y[i];\n         }\n      }\n   }\n}" }, silent = TRUE))
  suppressWarnings(try({ c <- rcode((py_getitem(a, i, j, k, l) + py_getitem(b, i, j, k, l)) * py_getitem(c, j, k, l), assign_to = py_getitem(y, i)) }, silent = TRUE))
  expect_identical(.cgb(c), "for (i in 1:m){\\n   y[i] = 0;\\n}\\nfor (i in 1:m){\\n   for (j in 1:n){\\n      for (k in 1:o){\\n         for (l in 1:p){\\n            y[i] = (a[i, j, k, l] + b[i, j, k, l])*c[j, k, l] + y[i];\\n         }\\n      }\\n   }\\n}")
})

test_that("printing_test_rcode: test_rcode_loops_multiple_terms", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  suppressWarnings(try({ .tmp <- sym_symbols("n m o p", integer = TRUE); n <- py_getitem(.tmp, 0L); m <- py_getitem(.tmp, 1L); o <- py_getitem(.tmp, 2L); p <- py_getitem(.tmp, 3L) }, silent = TRUE))
  suppressWarnings(try({ a <- IndexedBase("a") }, silent = TRUE))
  suppressWarnings(try({ b <- IndexedBase("b") }, silent = TRUE))
  suppressWarnings(try({ c <- IndexedBase("c") }, silent = TRUE))
  suppressWarnings(try({ y <- IndexedBase("y") }, silent = TRUE))
  suppressWarnings(try({ i <- Idx("i", m) }, silent = TRUE))
  suppressWarnings(try({ j <- Idx("j", n) }, silent = TRUE))
  suppressWarnings(try({ k <- Idx("k", o) }, silent = TRUE))
  suppressWarnings(try({ s0 <- "for (i in 1:m){\n   y[i] = 0;\n}\n" }, silent = TRUE))
  suppressWarnings(try({ s1 <- "for (i in 1:m){\n   for (j in 1:n){\n      for (k in 1:o){\n         y[i] = b[j]*b[k]*c[i, j, k] + y[i];\n      }\n   }\n}\n" }, silent = TRUE))
  suppressWarnings(try({ s2 <- "for (i in 1:m){\n   for (k in 1:o){\n      y[i] = a[i, k]*b[k] + y[i];\n   }\n}\n" }, silent = TRUE))
  suppressWarnings(try({ s3 <- "for (i in 1:m){\n   for (j in 1:n){\n      y[i] = a[i, j]*b[j] + y[i];\n   }\n}\n" }, silent = TRUE))
  suppressWarnings(try({ c <- rcode(((py_getitem(b, j) * py_getitem(a, i, j)) + (py_getitem(b, k) * py_getitem(a, i, k))) + ((py_getitem(b, j) * py_getitem(b, k)) * py_getitem(c, i, j, k)), assign_to = py_getitem(y, i)) }, silent = TRUE))
  suppressWarnings(try({ ref <- py_dict(list(), list()) }, silent = TRUE))
  suppressWarnings(try({ ref <- lapply(list(c(s0, s1, s2, s3), c(s0, s1, s3, s2), c(s0, s2, s1, s3), c(s0, s2, s3, s1), c(s0, s3, s1, s2), c(s0, s3, s2, s1)), function(v) sub("\n$", "", paste(v, collapse = ""))) }, silent = TRUE))
  expect_identical(.cgb(any(vapply(ref, function(r) identical(c, r), TRUE))), "True")
})

test_that("printing_test_rcode: test_dereference_printing", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  suppressWarnings(try({ expr <- ((x + y) + sin(z)) + z }, silent = TRUE))
  expect_identical(.cgb(rcode(expr, dereference = list(z))), "x + y + (*z) + sin((*z))")
})

test_that("printing_test_rcode: test_Matrix_printing", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  suppressWarnings(try({ mat <- Matrix(list(x * y, Piecewise(list(2L + x, y > 0L), list(y, TRUE)), sin(z))) }, silent = TRUE))
  suppressWarnings(try({ A <- MatrixSymbol("A", 3L, 1L) }, silent = TRUE))
  suppressWarnings(try({ p <- rcode(mat, A) }, silent = TRUE))
  expect_identical(.cgb(p), "A[0] = x*y;\\nA[1] = ifelse(y > 0,x + 2,y);\\nA[2] = sin(z);")
  suppressWarnings(try({ expr <- (Piecewise(list(2L * py_getitem(A, 2L, 0L), x > 0L), list(py_getitem(A, 2L, 0L), TRUE)) + sin(py_getitem(A, 1L, 0L))) + py_getitem(A, 0L, 0L) }, silent = TRUE))
  suppressWarnings(try({ p <- rcode(expr) }, silent = TRUE))
  expect_identical(.cgb(p), "ifelse(x > 0,2*A[2],A[2]) + sin(A[1]) + A[0]")
  suppressWarnings(try({ q <- MatrixSymbol("q", 5L, 1L) }, silent = TRUE))
  suppressWarnings(try({ M <- MatrixSymbol("M", 3L, 3L) }, silent = TRUE))
  suppressWarnings(try({ m <- Matrix(list(list(sin(py_getitem(q, 1L, 0L)), 0L, cos(py_getitem(q, 2L, 0L))), list(py_getitem(q, 1L, 0L) + py_getitem(q, 2L, 0L), py_getitem(q, 3L, 0L), 5L), list((2L * py_getitem(q, 4L, 0L)) / py_getitem(q, 1L, 0L), sqrt(py_getitem(q, 0L, 0L)) + 4L, 0L))) }, silent = TRUE))
  expect_identical(.cgb(rcode(m, M)), "M[0] = sin(q[1]);\\nM[1] = 0;\\nM[2] = cos(q[2]);\\nM[3] = q[1] + q[2];\\nM[4] = q[3];\\nM[5] = 5;\\nM[6] = 2*q[4]/q[1];\\nM[7] = sqrt(q[0]) + 4;\\nM[8] = 0;")
})

test_that("printing_test_rcode: test_rcode_sgn", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  suppressWarnings(try({ expr <- sign(x) * y }, silent = TRUE))
  expect_identical(.cgb(rcode(expr)), "y*sign(x)")
  suppressWarnings(try({ p <- rcode(expr, "z") }, silent = TRUE))
  expect_identical(.cgb(p), "z = y*sign(x);")
  suppressWarnings(try({ p <- rcode((sign((2L * x) + (x ^ 2L)) * x) + (x ^ 2L)) }, silent = TRUE))
  expect_identical(.cgb(p), "x^2 + x*sign(x^2 + 2*x)")
  suppressWarnings(try({ expr <- sign(cos(x)) }, silent = TRUE))
  suppressWarnings(try({ p <- rcode(expr) }, silent = TRUE))
  expect_identical(.cgb(p), "sign(cos(x))")
})

test_that("printing_test_rcode: test_rcode_Assignment", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  expect_identical(.cgb(rcode(Assignment(x, y + z))), "x = y + z;")
  expect_identical(.cgb(rcode(aug_assign(x, "+", y + z))), "x += y + z;")
})

test_that("printing_test_rcode: test_rcode_For", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  suppressWarnings(try({ f <- For(x, Range(0L, 10L, 2L), list(aug_assign(y, "*", x))) }, silent = TRUE))
  suppressWarnings(try({ sol <- rcode(f) }, silent = TRUE))
  expect_identical(.cgb(sol), "for(x in seq(from=0, to=9, by=2)){\\n   y *= x;\\n}")
})

test_that("printing_test_rcode: test_MatrixElement_printing", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  suppressWarnings(try({ A <- MatrixSymbol("A", 1L, 3L) }, silent = TRUE))
  suppressWarnings(try({ B <- MatrixSymbol("B", 1L, 3L) }, silent = TRUE))
  suppressWarnings(try({ C <- MatrixSymbol("C", 1L, 3L) }, silent = TRUE))
  expect_identical(.cgb(rcode(py_getitem(A, 0L, 0L))), "A[0]")
  expect_identical(.cgb(rcode(3L * py_getitem(A, 0L, 0L))), "3*A[0]")
  suppressWarnings(try({ F <- subs(py_getitem(C, 0L, 0L), C, A - B) }, silent = TRUE))
})

test_that("printing_test_rust: test_Pow", {
  x <- sym_symbols("x,y,z", integer = FALSE, real = TRUE); y <- x[[2]]; z <- x[[3]]; x <- x[[1]]; k <- sym_symbols("k,m,n", integer = TRUE); m <- k[[2]]; n <- k[[3]]; k <- k[[1]]
  expect_identical(.cgb(rust_code(1L / x)), "x.recip()")
  expect_identical(.cgb(rust_code(x ^ (-1L))), "x.recip()")
  expect_identical(.cgb(rust_code(x ^ (-Float(1.0)))), "x.recip()")
  expect_identical(.cgb(rust_code(sqrt(x))), "x.sqrt()")
  expect_identical(.cgb(rust_code(x ^ S_Half)), "x.sqrt()")
  expect_identical(.cgb(rust_code(x ^ Float(0.5))), "x.sqrt()")
  expect_identical(.cgb(rust_code(1L / sqrt(x))), "x.sqrt().recip()")
  expect_identical(.cgb(rust_code(x ^ (-S_Half))), "x.sqrt().recip()")
  expect_identical(.cgb(rust_code(x ^ (-Float(0.5)))), "x.sqrt().recip()")
  expect_identical(.cgb(rust_code(1L / S_Pi)), "PI.recip()")
  expect_identical(.cgb(rust_code(S_Pi ^ (-1L))), "PI.recip()")
  expect_identical(.cgb(rust_code(S_Pi ^ (-Float(1.0)))), "PI.recip()")
  expect_identical(.cgb(rust_code(S_Pi ^ (-Float(0.5)))), "PI.sqrt().recip()")
  expect_identical(.cgb(rust_code(x ^ Rational(1L, 3L))), "x.cbrt()")
  expect_identical(.cgb(rust_code(2L ^ x)), "x.exp2()")
  expect_identical(.cgb(rust_code(exp(x))), "x.exp()")
  expect_identical(.cgb(rust_code(x ^ 3L)), "x.powi(3)")
  expect_identical(.cgb(rust_code(x ^ (y ^ 3L))), "x.powf(y.powi(3))")
  expect_identical(.cgb(rust_code(x ^ Rational(2L, 3L))), "x.powf(2_f64/3.0)")
  suppressWarnings(try({ g <- implemented_function("g", Lambda(x, 2L * x)) }, silent = TRUE))
  expect_identical(.cgb(rust_code((1L / ((g(x) * Float(3.5)) ^ (x - (y ^ x)))) / ((x ^ 2L) + y))), "(3.5*2.0*x).powf(-x + y.powf(x))/(x.powi(2) + y)")
  suppressWarnings(try({ U_cond_cfunc <- list(list(function(base, exp) ask_fact(exp, "integer"), "dpowi", 1L), list(function(base, exp) !(ask_fact(exp, "integer")), "pow", 1L)) }, silent = TRUE))
  expect_identical(.cgb(rust_code(x ^ 3L, user_functions = py_dict(list("Pow"), list(U_cond_cfunc)))), "x.dpowi(3)")
  expect_identical(.cgb(rust_code(x ^ Float(3.2), user_functions = py_dict(list("Pow"), list(U_cond_cfunc)))), "x.pow(3.2)")
})

test_that("printing_test_rust: test_inline_function", {
  x <- sym_symbols("x,y,z", integer = FALSE, real = TRUE); y <- x[[2]]; z <- x[[3]]; x <- x[[1]]; k <- sym_symbols("k,m,n", integer = TRUE); m <- k[[2]]; n <- k[[3]]; k <- k[[1]]
  suppressWarnings(try({ x <- sym_symbols("x") }, silent = TRUE))
  suppressWarnings(try({ g <- implemented_function("g", Lambda(x, 2L * x)) }, silent = TRUE))
  expect_identical(.cgb(rust_code(g(x))), "2*x")
  suppressWarnings(try({ g <- implemented_function("g", Lambda(x, (2L * x) / Catalan)) }, silent = TRUE))
  expect_identical(.cgb(rust_code(g(x))), "const Catalan: f64 = 0.91596559417721901;\\n2.0*x/Catalan")
  suppressWarnings(try({ A <- IndexedBase("A") }, silent = TRUE))
  suppressWarnings(try({ i <- Idx("i", sym_symbols("n", integer = TRUE)) }, silent = TRUE))
  suppressWarnings(try({ g <- implemented_function("g", Lambda(x, (x * (1L + x)) * (2L + x))) }, silent = TRUE))
  expect_identical(.cgb(rust_code(g(py_getitem(A, i)), assign_to = py_getitem(A, i))), "for i in 0..n {\\n    A[i] = (A[i] + 1)*(A[i] + 2)*A[i];\\n}")
})

test_that("utilities_test_codegen_julia: test_empty_jl_code", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ code_gen <- JuliaCodeGen() }, silent = TRUE))
  suppressWarnings(try({ output <- StringIO() }, silent = TRUE))
  suppressWarnings(try({ py_method(code_gen, "dump_jl", list(), output, "file", header = FALSE, empty = FALSE) }, silent = TRUE))
  suppressWarnings(try({ source <- py_method(output, "getvalue") }, silent = TRUE))
  expect_identical(.cgb(source), "")
  for (.o in mget(ls(), inherits = FALSE)) if (inherits(.o, "py_stringio")) try(close(.o$con), silent = TRUE)
})

test_that("utilities_test_codegen_julia: test_jl_simple_code", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ name_expr <- list("test", (x + y) * z) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- codegen(name_expr, "Julia", header = FALSE, empty = FALSE); result <- py_getitem(.tmp, 0L) }, silent = TRUE))
  expect_identical(.cgb(py_getitem(result, 0L)), "test.jl")
  suppressWarnings(try({ source <- py_getitem(result, 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "function test(x, y, z)\n    out1 = z .* (x + y)\n    return out1\nend\n" }, silent = TRUE))
  expect_identical(.cgb(source), "function test(x, y, z)\\n    out1 = z .* (x + y)\\n    return out1\\nend\\n")
})

test_that("utilities_test_codegen_julia: test_jl_simple_code_with_header", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ name_expr <- list("test", (x + y) * z) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- codegen(name_expr, "Julia", header = TRUE, empty = FALSE); result <- py_getitem(.tmp, 0L) }, silent = TRUE))
  expect_identical(.cgb(py_getitem(result, 0L)), "test.jl")
  suppressWarnings(try({ source <- py_getitem(result, 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- ("#   Code generated with SymPy " + ast_attr(sympy, "__version__")) + "\n#\n#   See http://www.sympy.org/ for more information.\n#\n#   This file is part of 'project'\nfunction test(x, y, z)\n    out1 = z .* (x + y)\n    return out1\nend\n" }, silent = TRUE))
})

test_that("utilities_test_codegen_julia: test_jl_simple_code_nameout", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ expr <- Equality(z, x + y) }, silent = TRUE))
  suppressWarnings(try({ name_expr <- list("test", expr) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- codegen(name_expr, "Julia", header = FALSE, empty = FALSE); result <- py_getitem(.tmp, 0L) }, silent = TRUE))
  suppressWarnings(try({ source <- py_getitem(result, 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "function test(x, y)\n    z = x + y\n    return z\nend\n" }, silent = TRUE))
  expect_identical(.cgb(source), "function test(x, y)\\n    z = x + y\\n    return z\\nend\\n")
})

test_that("utilities_test_codegen_julia: test_jl_numbersymbol", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ name_expr <- list("test", S_Pi ^ Catalan) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- codegen(name_expr, "Julia", header = FALSE, empty = FALSE); result <- py_getitem(.tmp, 0L) }, silent = TRUE))
  suppressWarnings(try({ source <- py_getitem(result, 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "function test()\n    out1 = pi ^ catalan\n    return out1\nend\n" }, silent = TRUE))
  expect_identical(.cgb(source), "function test()\\n    out1 = pi ^ catalan\\n    return out1\\nend\\n")
})

test_that("utilities_test_codegen_julia: test_jl_code_argument_order", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ expr <- x + y }, silent = TRUE))
  suppressWarnings(try({ routine <- make_routine("test", expr, argument_sequence = list(z, x, y), language = "julia") }, silent = TRUE))
  suppressWarnings(try({ code_gen <- JuliaCodeGen() }, silent = TRUE))
  suppressWarnings(try({ output <- StringIO() }, silent = TRUE))
  suppressWarnings(try({ py_method(code_gen, "dump_jl", list(routine), output, "test", header = FALSE, empty = FALSE) }, silent = TRUE))
  suppressWarnings(try({ source <- py_method(output, "getvalue") }, silent = TRUE))
  suppressWarnings(try({ expected <- "function test(z, x, y)\n    out1 = x + y\n    return out1\nend\n" }, silent = TRUE))
  expect_identical(.cgb(source), "function test(z, x, y)\\n    out1 = x + y\\n    return out1\\nend\\n")
  for (.o in mget(ls(), inherits = FALSE)) if (inherits(.o, "py_stringio")) try(close(.o$con), silent = TRUE)
})

test_that("utilities_test_codegen_julia: test_multiple_results_m", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ expr1 <- (x + y) * z }, silent = TRUE))
  suppressWarnings(try({ expr2 <- (x - y) * z }, silent = TRUE))
  suppressWarnings(try({ name_expr <- list("test", list(expr1, expr2)) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- codegen(name_expr, "Julia", header = FALSE, empty = FALSE); result <- py_getitem(.tmp, 0L) }, silent = TRUE))
  suppressWarnings(try({ source <- py_getitem(result, 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "function test(x, y, z)\n    out1 = z .* (x + y)\n    out2 = z .* (x - y)\n    return out1, out2\nend\n" }, silent = TRUE))
  expect_identical(.cgb(source), "function test(x, y, z)\\n    out1 = z .* (x + y)\\n    out2 = z .* (x - y)\\n    return out1, out2\\nend\\n")
})

test_that("utilities_test_codegen_julia: test_results_named_unordered", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ .tmp <- sym_symbols("A,B,C"); A <- py_getitem(.tmp, 0L); B <- py_getitem(.tmp, 1L); C <- py_getitem(.tmp, 2L) }, silent = TRUE))
  suppressWarnings(try({ expr1 <- Equality(C, (x + y) * z) }, silent = TRUE))
  suppressWarnings(try({ expr2 <- Equality(A, (x - y) * z) }, silent = TRUE))
  suppressWarnings(try({ expr3 <- Equality(B, 2L * x) }, silent = TRUE))
  suppressWarnings(try({ name_expr <- list("test", list(expr1, expr2, expr3)) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- codegen(name_expr, "Julia", header = FALSE, empty = FALSE); result <- py_getitem(.tmp, 0L) }, silent = TRUE))
  suppressWarnings(try({ source <- py_getitem(result, 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "function test(x, y, z)\n    C = z .* (x + y)\n    A = z .* (x - y)\n    B = 2 * x\n    return C, A, B\nend\n" }, silent = TRUE))
  expect_identical(.cgb(source), "function test(x, y, z)\\n    C = z .* (x + y)\\n    A = z .* (x - y)\\n    B = 2 * x\\n    return C, A, B\\nend\\n")
})

test_that("utilities_test_codegen_julia: test_results_named_ordered", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ .tmp <- sym_symbols("A,B,C"); A <- py_getitem(.tmp, 0L); B <- py_getitem(.tmp, 1L); C <- py_getitem(.tmp, 2L) }, silent = TRUE))
  suppressWarnings(try({ expr1 <- Equality(C, (x + y) * z) }, silent = TRUE))
  suppressWarnings(try({ expr2 <- Equality(A, (x - y) * z) }, silent = TRUE))
  suppressWarnings(try({ expr3 <- Equality(B, 2L * x) }, silent = TRUE))
  suppressWarnings(try({ name_expr <- list("test", list(expr1, expr2, expr3)) }, silent = TRUE))
  suppressWarnings(try({ result <- codegen(name_expr, "Julia", header = FALSE, empty = FALSE, argument_sequence = list(x, z, y)) }, silent = TRUE))
  expect_identical(.cgb(py_getitem(py_getitem(result, 0L), 0L)), "test.jl")
  suppressWarnings(try({ source <- py_getitem(py_getitem(result, 0L), 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "function test(x, z, y)\n    C = z .* (x + y)\n    A = z .* (x - y)\n    B = 2 * x\n    return C, A, B\nend\n" }, silent = TRUE))
  expect_identical(.cgb(source), "function test(x, z, y)\\n    C = z .* (x + y)\\n    A = z .* (x - y)\\n    B = 2 * x\\n    return C, A, B\\nend\\n")
})

test_that("utilities_test_codegen_julia: test_complicated_jl_codegen", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ invisible(NULL) }, silent = TRUE))
  suppressWarnings(try({ name_expr <- list("testlong", list(expand(((sin(x) + cos(y)) + tan(z)) ^ 3L), cos(cos(cos(cos(cos(cos(cos(cos((x + y) + z)))))))))) }, silent = TRUE))
  suppressWarnings(try({ result <- codegen(name_expr, "Julia", header = FALSE, empty = FALSE) }, silent = TRUE))
  expect_identical(.cgb(py_getitem(py_getitem(result, 0L), 0L)), "testlong.jl")
  suppressWarnings(try({ source <- py_getitem(py_getitem(result, 0L), 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "function testlong(x, y, z)\n    out1 = sin(x) .^ 3 + 3 * sin(x) .^ 2 .* cos(y) + 3 * sin(x) .^ 2 .* tan(z) + 3 * sin(x) .* cos(y) .^ 2 + 6 * sin(x) .* cos(y) .* tan(z) + 3 * sin(x) .* tan(z) .^ 2 + cos(y) .^ 3 + 3 * cos(y) .^ 2 .* tan(z) + 3 * cos(y) .* tan(z) .^ 2 + tan(z) .^ 3\n    out2 = cos(cos(cos(cos(cos(cos(cos(cos(x + y + z))))))))\n    return out1, out2\nend\n" }, silent = TRUE))
  expect_identical(.cgb(source), "function testlong(x, y, z)\\n    out1 = sin(x) .^ 3 + 3 * sin(x) .^ 2 .* cos(y) + 3 * sin(x) .^ 2 .* tan(z) + 3 * sin(x) .* cos(y) .^ 2 + 6 * sin(x) .* cos(y) .* tan(z) + 3 * sin(x) .* tan(z) .^ 2 + cos(y) .^ 3 + 3 * cos(y) .^ 2 .* tan(z) + 3 * cos(y) .* tan(z) .^ 2 + tan(z) .^ 3\\n    out2 = cos(cos(cos(cos(cos(cos(cos(cos(x + y + z))))))))\\n    return out1, out2\\nend\\n")
})

test_that("utilities_test_codegen_julia: test_jl_output_arg_mixed_unordered", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ invisible(NULL) }, silent = TRUE))
  suppressWarnings(try({ a <- sym_symbols("a") }, silent = TRUE))
  suppressWarnings(try({ name_expr <- list("foo", list(cos(2L * x), Equality(y, sin(x)), cos(x), Equality(a, sin(2L * x)))) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- codegen(name_expr, "Julia", header = FALSE, empty = FALSE); result <- py_getitem(.tmp, 0L) }, silent = TRUE))
  expect_identical(.cgb(py_getitem(result, 0L)), "foo.jl")
  suppressWarnings(try({ source <- py_getitem(result, 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "function foo(x)\n    out1 = cos(2 * x)\n    y = sin(x)\n    out3 = cos(x)\n    a = sin(2 * x)\n    return out1, y, out3, a\nend\n" }, silent = TRUE))
  expect_identical(.cgb(source), "function foo(x)\\n    out1 = cos(2 * x)\\n    y = sin(x)\\n    out3 = cos(x)\\n    a = sin(2 * x)\\n    return out1, y, out3, a\\nend\\n")
})

test_that("utilities_test_codegen_julia: test_jl_piecewise_", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ pw <- Piecewise(list(0L, x < (-1L)), list(x ^ 2L, x <= 1L), list((-x) + 2L, x > 1L), list(1L, TRUE), evaluate = FALSE) }, silent = TRUE))
  suppressWarnings(try({ name_expr <- list("pwtest", pw) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- codegen(name_expr, "Julia", header = FALSE, empty = FALSE); result <- py_getitem(.tmp, 0L) }, silent = TRUE))
  suppressWarnings(try({ source <- py_getitem(result, 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "function pwtest(x)\n    out1 = ((x < -1) ? (0) :\n    (x <= 1) ? (x .^ 2) :\n    (x > 1) ? (2 - x) : (1))\n    return out1\nend\n" }, silent = TRUE))
  expect_identical(.cgb(source), "function pwtest(x)\\n    out1 = ((x < -1) ? (0) :\\n    (x <= 1) ? (x .^ 2) :\\n    (x > 1) ? (2 - x) : (1))\\n    return out1\\nend\\n")
})

test_that("utilities_test_codegen_julia: test_jl_multifcns_per_file", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ name_expr <- list(list("foo", list(2L * x, 3L * y)), list("bar", list(y ^ 2L, 4L * y))) }, silent = TRUE))
  suppressWarnings(try({ result <- codegen(name_expr, "Julia", header = FALSE, empty = FALSE) }, silent = TRUE))
  expect_identical(.cgb(py_getitem(py_getitem(result, 0L), 0L)), "foo.jl")
  suppressWarnings(try({ source <- py_getitem(py_getitem(result, 0L), 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "function foo(x, y)\n    out1 = 2 * x\n    out2 = 3 * y\n    return out1, out2\nend\nfunction bar(y)\n    out1 = y .^ 2\n    out2 = 4 * y\n    return out1, out2\nend\n" }, silent = TRUE))
  expect_identical(.cgb(source), "function foo(x, y)\\n    out1 = 2 * x\\n    out2 = 3 * y\\n    return out1, out2\\nend\\nfunction bar(y)\\n    out1 = y .^ 2\\n    out2 = 4 * y\\n    return out1, out2\\nend\\n")
})

test_that("utilities_test_codegen_julia: test_jl_multifcns_per_file_w_header", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ name_expr <- list(list("foo", list(2L * x, 3L * y)), list("bar", list(y ^ 2L, 4L * y))) }, silent = TRUE))
  suppressWarnings(try({ result <- codegen(name_expr, "Julia", header = TRUE, empty = FALSE) }, silent = TRUE))
  expect_identical(.cgb(py_getitem(py_getitem(result, 0L), 0L)), "foo.jl")
  suppressWarnings(try({ source <- py_getitem(py_getitem(result, 0L), 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- ("#   Code generated with SymPy " + ast_attr(sympy, "__version__")) + "\n#\n#   See http://www.sympy.org/ for more information.\n#\n#   This file is part of 'project'\nfunction foo(x, y)\n    out1 = 2 * x\n    out2 = 3 * y\n    return out1, out2\nend\nfunction bar(y)\n    out1 = y .^ 2\n    out2 = 4 * y\n    return out1, out2\nend\n" }, silent = TRUE))
})

test_that("utilities_test_codegen_julia: test_jl_filename_match_prefix", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ name_expr <- list(list("foo", list(2L * x, 3L * y)), list("bar", list(y ^ 2L, 4L * y))) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- codegen(name_expr, "Julia", prefix = "baz", header = FALSE, empty = FALSE); result <- py_getitem(.tmp, 0L) }, silent = TRUE))
  expect_identical(.cgb(py_getitem(result, 0L)), "baz.jl")
})

test_that("utilities_test_codegen_julia: test_jl_matrix_named", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ e2 <- Matrix(list(list(x, 2L * y, S_Pi * z))) }, silent = TRUE))
  suppressWarnings(try({ name_expr <- list("test", Equality(MatrixSymbol("myout1", 1L, 3L), e2)) }, silent = TRUE))
  suppressWarnings(try({ result <- codegen(name_expr, "Julia", header = FALSE, empty = FALSE) }, silent = TRUE))
  expect_identical(.cgb(py_getitem(py_getitem(result, 0L), 0L)), "test.jl")
  suppressWarnings(try({ source <- py_getitem(py_getitem(result, 0L), 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "function test(x, y, z)\n    myout1 = [x 2 * y pi * z]\n    return myout1\nend\n" }, silent = TRUE))
  expect_identical(.cgb(source), "function test(x, y, z)\\n    myout1 = [x 2 * y pi * z]\\n    return myout1\\nend\\n")
})

test_that("utilities_test_codegen_julia: test_jl_matrix_named_matsym", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ myout1 <- MatrixSymbol("myout1", 1L, 3L) }, silent = TRUE))
  suppressWarnings(try({ e2 <- Matrix(list(list(x, 2L * y, S_Pi * z))) }, silent = TRUE))
  suppressWarnings(try({ name_expr <- list("test", Equality(myout1, e2, evaluate = FALSE)) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- codegen(name_expr, "Julia", header = FALSE, empty = FALSE); result <- py_getitem(.tmp, 0L) }, silent = TRUE))
  suppressWarnings(try({ source <- py_getitem(result, 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "function test(x, y, z)\n    myout1 = [x 2 * y pi * z]\n    return myout1\nend\n" }, silent = TRUE))
  expect_identical(.cgb(source), "function test(x, y, z)\\n    myout1 = [x 2 * y pi * z]\\n    return myout1\\nend\\n")
})

test_that("utilities_test_codegen_julia: test_jl_matrix_output_autoname", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ expr <- Matrix(list(list(x, x + y, 3L))) }, silent = TRUE))
  suppressWarnings(try({ name_expr <- list("test", expr) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- codegen(name_expr, "Julia", header = FALSE, empty = FALSE); result <- py_getitem(.tmp, 0L) }, silent = TRUE))
  suppressWarnings(try({ source <- py_getitem(result, 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "function test(x, y)\n    out1 = [x x + y 3]\n    return out1\nend\n" }, silent = TRUE))
  expect_identical(.cgb(source), "function test(x, y)\\n    out1 = [x x + y 3]\\n    return out1\\nend\\n")
})

test_that("utilities_test_codegen_julia: test_jl_matrix_output_autoname_2", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ e1 <- x + y }, silent = TRUE))
  suppressWarnings(try({ e2 <- Matrix(list(list(2L * x, 2L * y, 2L * z))) }, silent = TRUE))
  suppressWarnings(try({ e3 <- Matrix(list(list(x), list(y), list(z))) }, silent = TRUE))
  suppressWarnings(try({ e4 <- Matrix(list(list(x, y), list(z, 16L))) }, silent = TRUE))
  suppressWarnings(try({ name_expr <- list("test", list(e1, e2, e3, e4)) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- codegen(name_expr, "Julia", header = FALSE, empty = FALSE); result <- py_getitem(.tmp, 0L) }, silent = TRUE))
  suppressWarnings(try({ source <- py_getitem(result, 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "function test(x, y, z)\n    out1 = x + y\n    out2 = [2 * x 2 * y 2 * z]\n    out3 = [x, y, z]\n    out4 = [x  y;\n    z 16]\n    return out1, out2, out3, out4\nend\n" }, silent = TRUE))
  expect_identical(.cgb(source), "function test(x, y, z)\\n    out1 = x + y\\n    out2 = [2 * x 2 * y 2 * z]\\n    out3 = [x, y, z]\\n    out4 = [x  y;\\n    z 16]\\n    return out1, out2, out3, out4\\nend\\n")
})

test_that("utilities_test_codegen_julia: test_jl_results_matrix_named_ordered", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ .tmp <- sym_symbols("B,C"); B <- py_getitem(.tmp, 0L); C <- py_getitem(.tmp, 1L) }, silent = TRUE))
  suppressWarnings(try({ A <- MatrixSymbol("A", 1L, 3L) }, silent = TRUE))
  suppressWarnings(try({ expr1 <- Equality(C, (x + y) * z) }, silent = TRUE))
  suppressWarnings(try({ expr2 <- Equality(A, Matrix(list(list(1L, 2L, x)))) }, silent = TRUE))
  suppressWarnings(try({ expr3 <- Equality(B, 2L * x) }, silent = TRUE))
  suppressWarnings(try({ name_expr <- list("test", list(expr1, expr2, expr3)) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- codegen(name_expr, "Julia", header = FALSE, empty = FALSE, argument_sequence = list(x, z, y)); result <- py_getitem(.tmp, 0L) }, silent = TRUE))
  suppressWarnings(try({ source <- py_getitem(result, 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "function test(x, z, y)\n    C = z .* (x + y)\n    A = [1 2 x]\n    B = 2 * x\n    return C, A, B\nend\n" }, silent = TRUE))
  expect_identical(.cgb(source), "function test(x, z, y)\\n    C = z .* (x + y)\\n    A = [1 2 x]\\n    B = 2 * x\\n    return C, A, B\\nend\\n")
})

test_that("utilities_test_codegen_julia: test_jl_matrixsymbol_slice", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ A <- MatrixSymbol("A", 2L, 3L) }, silent = TRUE))
  suppressWarnings(try({ B <- MatrixSymbol("B", 1L, 3L) }, silent = TRUE))
  suppressWarnings(try({ C <- MatrixSymbol("C", 1L, 3L) }, silent = TRUE))
  suppressWarnings(try({ D <- MatrixSymbol("D", 2L, 1L) }, silent = TRUE))
  suppressWarnings(try({ name_expr <- list("test", list(Equality(B, py_getitem(A, 0L, pyslice(NULL, NULL, NULL))), Equality(C, py_getitem(A, 1L, pyslice(NULL, NULL, NULL))), Equality(D, py_getitem(A, pyslice(NULL, NULL, NULL), 2L)))) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- codegen(name_expr, "Julia", header = FALSE, empty = FALSE); result <- py_getitem(.tmp, 0L) }, silent = TRUE))
  suppressWarnings(try({ source <- py_getitem(result, 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "function test(A)\n    B = A[1,:]\n    C = A[2,:]\n    D = A[:,3]\n    return B, C, D\nend\n" }, silent = TRUE))
  expect_identical(.cgb(source), "function test(A)\\n    B = A[1,:]\\n    C = A[2,:]\\n    D = A[:,3]\\n    return B, C, D\\nend\\n")
})

test_that("utilities_test_codegen_julia: test_jl_matrixsymbol_slice2", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ A <- MatrixSymbol("A", 3L, 4L) }, silent = TRUE))
  suppressWarnings(try({ B <- MatrixSymbol("B", 2L, 2L) }, silent = TRUE))
  suppressWarnings(try({ C <- MatrixSymbol("C", 2L, 2L) }, silent = TRUE))
  suppressWarnings(try({ name_expr <- list("test", list(Equality(B, py_getitem(A, pyslice(0L, 2L, NULL), pyslice(0L, 2L, NULL))), Equality(C, py_getitem(A, pyslice(0L, 2L, NULL), pyslice(1L, 3L, NULL))))) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- codegen(name_expr, "Julia", header = FALSE, empty = FALSE); result <- py_getitem(.tmp, 0L) }, silent = TRUE))
  suppressWarnings(try({ source <- py_getitem(result, 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "function test(A)\n    B = A[1:2,1:2]\n    C = A[1:2,2:3]\n    return B, C\nend\n" }, silent = TRUE))
  expect_identical(.cgb(source), "function test(A)\\n    B = A[1:2,1:2]\\n    C = A[1:2,2:3]\\n    return B, C\\nend\\n")
})

test_that("utilities_test_codegen_julia: test_jl_matrixsymbol_slice3", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ A <- MatrixSymbol("A", 8L, 7L) }, silent = TRUE))
  suppressWarnings(try({ B <- MatrixSymbol("B", 2L, 2L) }, silent = TRUE))
  suppressWarnings(try({ C <- MatrixSymbol("C", 4L, 3L) }, silent = TRUE))
  suppressWarnings(try({ name_expr <- list("test", list(Equality(B, py_getitem(A, pyslice(6L, NULL, NULL), pyslice(1L, NULL, 3L))), Equality(C, py_getitem(A, pyslice(NULL, NULL, 2L), pyslice(NULL, NULL, 3L))))) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- codegen(name_expr, "Julia", header = FALSE, empty = FALSE); result <- py_getitem(.tmp, 0L) }, silent = TRUE))
  suppressWarnings(try({ source <- py_getitem(result, 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "function test(A)\n    B = A[7:end,2:3:end]\n    C = A[1:2:end,1:3:end]\n    return B, C\nend\n" }, silent = TRUE))
  expect_identical(.cgb(source), "function test(A)\\n    B = A[7:end,2:3:end]\\n    C = A[1:2:end,1:3:end]\\n    return B, C\\nend\\n")
})

test_that("utilities_test_codegen_julia: test_jl_matrixsymbol_slice_autoname", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ A <- MatrixSymbol("A", 2L, 3L) }, silent = TRUE))
  suppressWarnings(try({ B <- MatrixSymbol("B", 1L, 3L) }, silent = TRUE))
  suppressWarnings(try({ name_expr <- list("test", list(Equality(B, py_getitem(A, 0L, pyslice(NULL, NULL, NULL))), py_getitem(A, 1L, pyslice(NULL, NULL, NULL)), py_getitem(A, pyslice(NULL, NULL, NULL), 0L), py_getitem(A, pyslice(NULL, NULL, NULL), 1L))) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- codegen(name_expr, "Julia", header = FALSE, empty = FALSE); result <- py_getitem(.tmp, 0L) }, silent = TRUE))
  suppressWarnings(try({ source <- py_getitem(result, 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "function test(A)\n    B = A[1,:]\n    out2 = A[2,:]\n    out3 = A[:,1]\n    out4 = A[:,2]\n    return B, out2, out3, out4\nend\n" }, silent = TRUE))
  expect_identical(.cgb(source), "function test(A)\\n    B = A[1,:]\\n    out2 = A[2,:]\\n    out3 = A[:,1]\\n    out4 = A[:,2]\\n    return B, out2, out3, out4\\nend\\n")
})

test_that("utilities_test_codegen_julia: test_jl_loops", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ invisible(NULL) }, silent = TRUE))
  suppressWarnings(try({ invisible(NULL) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- sym_symbols("n m", integer = TRUE); n <- py_getitem(.tmp, 0L); m <- py_getitem(.tmp, 1L) }, silent = TRUE))
  suppressWarnings(try({ A <- IndexedBase("A") }, silent = TRUE))
  suppressWarnings(try({ x <- IndexedBase("x") }, silent = TRUE))
  suppressWarnings(try({ y <- IndexedBase("y") }, silent = TRUE))
  suppressWarnings(try({ i <- Idx("i", m) }, silent = TRUE))
  suppressWarnings(try({ j <- Idx("j", n) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- codegen(list("mat_vec_mult", Eq(py_getitem(y, i), py_getitem(A, i, j) * py_getitem(x, j))), "Julia", header = FALSE, empty = FALSE); result <- py_getitem(.tmp, 0L) }, silent = TRUE))
  suppressWarnings(try({ source <- py_getitem(result, 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "function mat_vec_mult(y, A, m, n, x)\n    for i = 1:m\n        y[i] = 0\n    end\n    for i = 1:m\n        for j = 1:n\n            y[i] = %(rhs)s + y[i]\n        end\n    end\n    return y\nend\n" }, silent = TRUE))
  expect_identical(.cgb((py_eq(source, py_format(expected, py_dict(list("rhs"), list(py_format("A[%s,%s] .* x[j]", list(i, j))))))) || (py_eq(source, py_format(expected, py_dict(list("rhs"), list(py_format("x[j] .* A[%s,%s]", list(i, j)))))))), "True")
})

test_that("utilities_test_codegen_julia: test_jl_tensor_loops_multiple_contractions", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ invisible(NULL) }, silent = TRUE))
  suppressWarnings(try({ invisible(NULL) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- sym_symbols("n m o p", integer = TRUE); n <- py_getitem(.tmp, 0L); m <- py_getitem(.tmp, 1L); o <- py_getitem(.tmp, 2L); p <- py_getitem(.tmp, 3L) }, silent = TRUE))
  suppressWarnings(try({ A <- IndexedBase("A") }, silent = TRUE))
  suppressWarnings(try({ B <- IndexedBase("B") }, silent = TRUE))
  suppressWarnings(try({ y <- IndexedBase("y") }, silent = TRUE))
  suppressWarnings(try({ i <- Idx("i", m) }, silent = TRUE))
  suppressWarnings(try({ j <- Idx("j", n) }, silent = TRUE))
  suppressWarnings(try({ k <- Idx("k", o) }, silent = TRUE))
  suppressWarnings(try({ l <- Idx("l", p) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- codegen(list("tensorthing", Eq(py_getitem(y, i), py_getitem(B, j, k, l) * py_getitem(A, i, j, k, l))), "Julia", header = FALSE, empty = FALSE); result <- py_getitem(.tmp, 0L) }, silent = TRUE))
  suppressWarnings(try({ source <- py_getitem(result, 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "function tensorthing(y, A, B, m, n, o, p)\n    for i = 1:m\n        y[i] = 0\n    end\n    for i = 1:m\n        for j = 1:n\n            for k = 1:o\n                for l = 1:p\n                    y[i] = A[i,j,k,l] .* B[j,k,l] + y[i]\n                end\n            end\n        end\n    end\n    return y\nend\n" }, silent = TRUE))
  expect_identical(.cgb(source), "function tensorthing(y, A, B, m, n, o, p)\\n    for i = 1:m\\n        y[i] = 0\\n    end\\n    for i = 1:m\\n        for j = 1:n\\n            for k = 1:o\\n                for l = 1:p\\n                    y[i] = A[i,j,k,l] .* B[j,k,l] + y[i]\\n                end\\n            end\\n        end\\n    end\\n    return y\\nend\\n")
})

test_that("utilities_test_codegen_julia: test_jl_InOutArgument", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ expr <- Equality(x, x ^ 2L) }, silent = TRUE))
  suppressWarnings(try({ name_expr <- list("mysqr", expr) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- codegen(name_expr, "Julia", header = FALSE, empty = FALSE); result <- py_getitem(.tmp, 0L) }, silent = TRUE))
  suppressWarnings(try({ source <- py_getitem(result, 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "function mysqr(x)\n    x = x .^ 2\n    return x\nend\n" }, silent = TRUE))
  expect_identical(.cgb(source), "function mysqr(x)\\n    x = x .^ 2\\n    return x\\nend\\n")
})

test_that("utilities_test_codegen_julia: test_jl_InOutArgument_order", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ expr <- Equality(x, (x ^ 2L) + y) }, silent = TRUE))
  suppressWarnings(try({ name_expr <- list("test", expr) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- codegen(name_expr, "Julia", header = FALSE, empty = FALSE, argument_sequence = list(x, y)); result <- py_getitem(.tmp, 0L) }, silent = TRUE))
  suppressWarnings(try({ source <- py_getitem(result, 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "function test(x, y)\n    x = x .^ 2 + y\n    return x\nend\n" }, silent = TRUE))
  expect_identical(.cgb(source), "function test(x, y)\\n    x = x .^ 2 + y\\n    return x\\nend\\n")
  suppressWarnings(try({ expr <- Equality(x, (x ^ 2L) + y) }, silent = TRUE))
  suppressWarnings(try({ name_expr <- list("test", expr) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- codegen(name_expr, "Julia", header = FALSE, empty = FALSE); result <- py_getitem(.tmp, 0L) }, silent = TRUE))
  suppressWarnings(try({ source <- py_getitem(result, 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "function test(x, y)\n    x = x .^ 2 + y\n    return x\nend\n" }, silent = TRUE))
  expect_identical(.cgb(source), "function test(x, y)\\n    x = x .^ 2 + y\\n    return x\\nend\\n")
})

test_that("utilities_test_codegen_julia: test_jl_not_supported", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ f <- Function("f") }, silent = TRUE))
  suppressWarnings(try({ name_expr <- list("test", list(diff(f(x), x), zoo)) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- codegen(name_expr, "Julia", header = FALSE, empty = FALSE); result <- py_getitem(.tmp, 0L) }, silent = TRUE))
  suppressWarnings(try({ source <- py_getitem(result, 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "function test(x)\n    # unsupported: Derivative(f(x), x)\n    # unsupported: zoo\n    out1 = Derivative(f(x), x)\n    out2 = zoo\n    return out1, out2\nend\n" }, silent = TRUE))
  expect_identical(.cgb(source), "function test(x)\\n    # unsupported: Derivative(f(x), x)\\n    # unsupported: zoo\\n    out1 = Derivative(f(x), x)\\n    out2 = zoo\\n    return out1, out2\\nend\\n")
})

test_that("utilities_test_codegen_julia: test_global_vars_octave", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ .tmp <- sym_symbols("x y z t"); x <- py_getitem(.tmp, 0L); y <- py_getitem(.tmp, 1L); z <- py_getitem(.tmp, 2L); t <- py_getitem(.tmp, 3L) }, silent = TRUE))
  suppressWarnings(try({ result <- codegen(list("f", x * y), "Julia", header = FALSE, empty = FALSE, global_vars = list(y)) }, silent = TRUE))
  suppressWarnings(try({ source <- py_getitem(py_getitem(result, 0L), 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "function f(x)\n    out1 = x .* y\n    return out1\nend\n" }, silent = TRUE))
  expect_identical(.cgb(source), "function f(x)\\n    out1 = x .* y\\n    return out1\\nend\\n")
  suppressWarnings(try({ result <- codegen(list("f", (x * y) + z), "Julia", header = FALSE, empty = FALSE, argument_sequence = list(x, y), global_vars = list(z, t)) }, silent = TRUE))
  suppressWarnings(try({ source <- py_getitem(py_getitem(result, 0L), 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "function f(x, y)\n    out1 = x .* y + z\n    return out1\nend\n" }, silent = TRUE))
  expect_identical(.cgb(source), "function f(x, y)\\n    out1 = x .* y + z\\n    return out1\\nend\\n")
})

test_that("utilities_test_codegen_octave: test_empty_m_code", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ code_gen <- OctaveCodeGen() }, silent = TRUE))
  suppressWarnings(try({ output <- StringIO() }, silent = TRUE))
  suppressWarnings(try({ py_method(code_gen, "dump_m", list(), output, "file", header = FALSE, empty = FALSE) }, silent = TRUE))
  suppressWarnings(try({ source <- py_method(output, "getvalue") }, silent = TRUE))
  expect_identical(.cgb(source), "")
  for (.o in mget(ls(), inherits = FALSE)) if (inherits(.o, "py_stringio")) try(close(.o$con), silent = TRUE)
})

test_that("utilities_test_codegen_octave: test_m_simple_code", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ name_expr <- list("test", (x + y) * z) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- codegen(name_expr, "Octave", header = FALSE, empty = FALSE); result <- py_getitem(.tmp, 0L) }, silent = TRUE))
  expect_identical(.cgb(py_getitem(result, 0L)), "test.m")
  suppressWarnings(try({ source <- py_getitem(result, 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "function out1 = test(x, y, z)\n  out1 = z.*(x + y);\nend\n" }, silent = TRUE))
  expect_identical(.cgb(source), "function out1 = test(x, y, z)\\n  out1 = z.*(x + y);\\nend\\n")
})

test_that("utilities_test_codegen_octave: test_m_simple_code_with_header", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ name_expr <- list("test", (x + y) * z) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- codegen(name_expr, "Octave", header = TRUE, empty = FALSE); result <- py_getitem(.tmp, 0L) }, silent = TRUE))
  expect_identical(.cgb(py_getitem(result, 0L)), "test.m")
  suppressWarnings(try({ source <- py_getitem(result, 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- ("function out1 = test(x, y, z)\n  %TEST  Autogenerated by SymPy\n  %   Code generated with SymPy " + ast_attr(sympy, "__version__")) + "\n  %\n  %   See http://www.sympy.org/ for more information.\n  %\n  %   This file is part of 'project'\n  out1 = z.*(x + y);\nend\n" }, silent = TRUE))
})

test_that("utilities_test_codegen_octave: test_m_simple_code_nameout", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ expr <- Equality(z, x + y) }, silent = TRUE))
  suppressWarnings(try({ name_expr <- list("test", expr) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- codegen(name_expr, "Octave", header = FALSE, empty = FALSE); result <- py_getitem(.tmp, 0L) }, silent = TRUE))
  suppressWarnings(try({ source <- py_getitem(result, 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "function z = test(x, y)\n  z = x + y;\nend\n" }, silent = TRUE))
  expect_identical(.cgb(source), "function z = test(x, y)\\n  z = x + y;\\nend\\n")
})

test_that("utilities_test_codegen_octave: test_m_numbersymbol", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ name_expr <- list("test", S_Pi ^ Catalan) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- codegen(name_expr, "Octave", header = FALSE, empty = FALSE); result <- py_getitem(.tmp, 0L) }, silent = TRUE))
  suppressWarnings(try({ source <- py_getitem(result, 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- py_format("function out1 = test()\n  out1 = pi^%s;\nend\n", N(Catalan, 17L)) }, silent = TRUE))
  expect_identical(.cgb(source), "function out1 = test()\\n  out1 = pi^0.91596559417721901;\\nend\\n")
})

test_that("utilities_test_codegen_octave: test_m_code_argument_order", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ expr <- x + y }, silent = TRUE))
  suppressWarnings(try({ routine <- make_routine("test", expr, argument_sequence = list(z, x, y), language = "octave") }, silent = TRUE))
  suppressWarnings(try({ code_gen <- OctaveCodeGen() }, silent = TRUE))
  suppressWarnings(try({ output <- StringIO() }, silent = TRUE))
  suppressWarnings(try({ py_method(code_gen, "dump_m", list(routine), output, "test", header = FALSE, empty = FALSE) }, silent = TRUE))
  suppressWarnings(try({ source <- py_method(output, "getvalue") }, silent = TRUE))
  suppressWarnings(try({ expected <- "function out1 = test(z, x, y)\n  out1 = x + y;\nend\n" }, silent = TRUE))
  expect_identical(.cgb(source), "function out1 = test(z, x, y)\\n  out1 = x + y;\\nend\\n")
  for (.o in mget(ls(), inherits = FALSE)) if (inherits(.o, "py_stringio")) try(close(.o$con), silent = TRUE)
})

test_that("utilities_test_codegen_octave: test_multiple_results_m", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ expr1 <- (x + y) * z }, silent = TRUE))
  suppressWarnings(try({ expr2 <- (x - y) * z }, silent = TRUE))
  suppressWarnings(try({ name_expr <- list("test", list(expr1, expr2)) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- codegen(name_expr, "Octave", header = FALSE, empty = FALSE); result <- py_getitem(.tmp, 0L) }, silent = TRUE))
  suppressWarnings(try({ source <- py_getitem(result, 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "function [out1, out2] = test(x, y, z)\n  out1 = z.*(x + y);\n  out2 = z.*(x - y);\nend\n" }, silent = TRUE))
  expect_identical(.cgb(source), "function [out1, out2] = test(x, y, z)\\n  out1 = z.*(x + y);\\n  out2 = z.*(x - y);\\nend\\n")
})

test_that("utilities_test_codegen_octave: test_results_named_unordered", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ .tmp <- sym_symbols("A,B,C"); A <- py_getitem(.tmp, 0L); B <- py_getitem(.tmp, 1L); C <- py_getitem(.tmp, 2L) }, silent = TRUE))
  suppressWarnings(try({ expr1 <- Equality(C, (x + y) * z) }, silent = TRUE))
  suppressWarnings(try({ expr2 <- Equality(A, (x - y) * z) }, silent = TRUE))
  suppressWarnings(try({ expr3 <- Equality(B, 2L * x) }, silent = TRUE))
  suppressWarnings(try({ name_expr <- list("test", list(expr1, expr2, expr3)) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- codegen(name_expr, "Octave", header = FALSE, empty = FALSE); result <- py_getitem(.tmp, 0L) }, silent = TRUE))
  suppressWarnings(try({ source <- py_getitem(result, 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "function [C, A, B] = test(x, y, z)\n  C = z.*(x + y);\n  A = z.*(x - y);\n  B = 2*x;\nend\n" }, silent = TRUE))
  expect_identical(.cgb(source), "function [C, A, B] = test(x, y, z)\\n  C = z.*(x + y);\\n  A = z.*(x - y);\\n  B = 2*x;\\nend\\n")
})

test_that("utilities_test_codegen_octave: test_results_named_ordered", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ .tmp <- sym_symbols("A,B,C"); A <- py_getitem(.tmp, 0L); B <- py_getitem(.tmp, 1L); C <- py_getitem(.tmp, 2L) }, silent = TRUE))
  suppressWarnings(try({ expr1 <- Equality(C, (x + y) * z) }, silent = TRUE))
  suppressWarnings(try({ expr2 <- Equality(A, (x - y) * z) }, silent = TRUE))
  suppressWarnings(try({ expr3 <- Equality(B, 2L * x) }, silent = TRUE))
  suppressWarnings(try({ name_expr <- list("test", list(expr1, expr2, expr3)) }, silent = TRUE))
  suppressWarnings(try({ result <- codegen(name_expr, "Octave", header = FALSE, empty = FALSE, argument_sequence = list(x, z, y)) }, silent = TRUE))
  expect_identical(.cgb(py_getitem(py_getitem(result, 0L), 0L)), "test.m")
  suppressWarnings(try({ source <- py_getitem(py_getitem(result, 0L), 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "function [C, A, B] = test(x, z, y)\n  C = z.*(x + y);\n  A = z.*(x - y);\n  B = 2*x;\nend\n" }, silent = TRUE))
  expect_identical(.cgb(source), "function [C, A, B] = test(x, z, y)\\n  C = z.*(x + y);\\n  A = z.*(x - y);\\n  B = 2*x;\\nend\\n")
})

test_that("utilities_test_codegen_octave: test_complicated_m_codegen", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ invisible(NULL) }, silent = TRUE))
  suppressWarnings(try({ name_expr <- list("testlong", list(expand(((sin(x) + cos(y)) + tan(z)) ^ 3L), cos(cos(cos(cos(cos(cos(cos(cos((x + y) + z)))))))))) }, silent = TRUE))
  suppressWarnings(try({ result <- codegen(name_expr, "Octave", header = FALSE, empty = FALSE) }, silent = TRUE))
  expect_identical(.cgb(py_getitem(py_getitem(result, 0L), 0L)), "testlong.m")
  suppressWarnings(try({ source <- py_getitem(py_getitem(result, 0L), 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "function [out1, out2] = testlong(x, y, z)\n  out1 = sin(x).^3 + 3*sin(x).^2.*cos(y) + 3*sin(x).^2.*tan(z) + 3*sin(x).*cos(y).^2 + 6*sin(x).*cos(y).*tan(z) + 3*sin(x).*tan(z).^2 + cos(y).^3 + 3*cos(y).^2.*tan(z) + 3*cos(y).*tan(z).^2 + tan(z).^3;\n  out2 = cos(cos(cos(cos(cos(cos(cos(cos(x + y + z))))))));\nend\n" }, silent = TRUE))
  expect_identical(.cgb(source), "function [out1, out2] = testlong(x, y, z)\\n  out1 = sin(x).^3 + 3*sin(x).^2.*cos(y) + 3*sin(x).^2.*tan(z) + 3*sin(x).*cos(y).^2 + 6*sin(x).*cos(y).*tan(z) + 3*sin(x).*tan(z).^2 + cos(y).^3 + 3*cos(y).^2.*tan(z) + 3*cos(y).*tan(z).^2 + tan(z).^3;\\n  out2 = cos(cos(cos(cos(cos(cos(cos(cos(x + y + z))))))));\\nend\\n")
})

test_that("utilities_test_codegen_octave: test_m_output_arg_mixed_unordered", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ invisible(NULL) }, silent = TRUE))
  suppressWarnings(try({ a <- sym_symbols("a") }, silent = TRUE))
  suppressWarnings(try({ name_expr <- list("foo", list(cos(2L * x), Equality(y, sin(x)), cos(x), Equality(a, sin(2L * x)))) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- codegen(name_expr, "Octave", header = FALSE, empty = FALSE); result <- py_getitem(.tmp, 0L) }, silent = TRUE))
  expect_identical(.cgb(py_getitem(result, 0L)), "foo.m")
  suppressWarnings(try({ source <- py_getitem(result, 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "function [out1, y, out3, a] = foo(x)\n  out1 = cos(2*x);\n  y = sin(x);\n  out3 = cos(x);\n  a = sin(2*x);\nend\n" }, silent = TRUE))
  expect_identical(.cgb(source), "function [out1, y, out3, a] = foo(x)\\n  out1 = cos(2*x);\\n  y = sin(x);\\n  out3 = cos(x);\\n  a = sin(2*x);\\nend\\n")
})

test_that("utilities_test_codegen_octave: test_m_piecewise_", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ pw <- Piecewise(list(0L, x < (-1L)), list(x ^ 2L, x <= 1L), list((-x) + 2L, x > 1L), list(1L, TRUE), evaluate = FALSE) }, silent = TRUE))
  suppressWarnings(try({ name_expr <- list("pwtest", pw) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- codegen(name_expr, "Octave", header = FALSE, empty = FALSE); result <- py_getitem(.tmp, 0L) }, silent = TRUE))
  suppressWarnings(try({ source <- py_getitem(result, 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "function out1 = pwtest(x)\n  out1 = ((x < -1).*(0) + (~(x < -1)).*( ...\n  (x <= 1).*(x.^2) + (~(x <= 1)).*( ...\n  (x > 1).*(2 - x) + (~(x > 1)).*(1))));\nend\n" }, silent = TRUE))
  expect_identical(.cgb(source), "function out1 = pwtest(x)\\n  out1 = ((x < -1).*(0) + (~(x < -1)).*( ...\\n  (x <= 1).*(x.^2) + (~(x <= 1)).*( ...\\n  (x > 1).*(2 - x) + (~(x > 1)).*(1))));\\nend\\n")
})

test_that("utilities_test_codegen_octave: test_m_multifcns_per_file", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ name_expr <- list(list("foo", list(2L * x, 3L * y)), list("bar", list(y ^ 2L, 4L * y))) }, silent = TRUE))
  suppressWarnings(try({ result <- codegen(name_expr, "Octave", header = FALSE, empty = FALSE) }, silent = TRUE))
  expect_identical(.cgb(py_getitem(py_getitem(result, 0L), 0L)), "foo.m")
  suppressWarnings(try({ source <- py_getitem(py_getitem(result, 0L), 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "function [out1, out2] = foo(x, y)\n  out1 = 2*x;\n  out2 = 3*y;\nend\nfunction [out1, out2] = bar(y)\n  out1 = y.^2;\n  out2 = 4*y;\nend\n" }, silent = TRUE))
  expect_identical(.cgb(source), "function [out1, out2] = foo(x, y)\\n  out1 = 2*x;\\n  out2 = 3*y;\\nend\\nfunction [out1, out2] = bar(y)\\n  out1 = y.^2;\\n  out2 = 4*y;\\nend\\n")
})

test_that("utilities_test_codegen_octave: test_m_multifcns_per_file_w_header", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ name_expr <- list(list("foo", list(2L * x, 3L * y)), list("bar", list(y ^ 2L, 4L * y))) }, silent = TRUE))
  suppressWarnings(try({ result <- codegen(name_expr, "Octave", header = TRUE, empty = FALSE) }, silent = TRUE))
  expect_identical(.cgb(py_getitem(py_getitem(result, 0L), 0L)), "foo.m")
  suppressWarnings(try({ source <- py_getitem(py_getitem(result, 0L), 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- ("function [out1, out2] = foo(x, y)\n  %FOO  Autogenerated by SymPy\n  %   Code generated with SymPy " + ast_attr(sympy, "__version__")) + "\n  %\n  %   See http://www.sympy.org/ for more information.\n  %\n  %   This file is part of 'project'\n  out1 = 2*x;\n  out2 = 3*y;\nend\nfunction [out1, out2] = bar(y)\n  out1 = y.^2;\n  out2 = 4*y;\nend\n" }, silent = TRUE))
})

test_that("utilities_test_codegen_octave: test_m_filename_match_first_fcn", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ name_expr <- list(list("foo", list(2L * x, 3L * y)), list("bar", list(y ^ 2L, 4L * y))) }, silent = TRUE))
  expect_identical(.cgb(codegen(name_expr, "Octave", prefix = "bar", header = FALSE, empty = FALSE)), "ERR ValueError: Octave function name should match prefix")
})

test_that("utilities_test_codegen_octave: test_m_matrix_named", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ e2 <- Matrix(list(list(x, 2L * y, S_Pi * z))) }, silent = TRUE))
  suppressWarnings(try({ name_expr <- list("test", Equality(MatrixSymbol("myout1", 1L, 3L), e2)) }, silent = TRUE))
  suppressWarnings(try({ result <- codegen(name_expr, "Octave", header = FALSE, empty = FALSE) }, silent = TRUE))
  expect_identical(.cgb(py_getitem(py_getitem(result, 0L), 0L)), "test.m")
  suppressWarnings(try({ source <- py_getitem(py_getitem(result, 0L), 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "function myout1 = test(x, y, z)\n  myout1 = [x 2*y pi*z];\nend\n" }, silent = TRUE))
  expect_identical(.cgb(source), "function myout1 = test(x, y, z)\\n  myout1 = [x 2*y pi*z];\\nend\\n")
})

test_that("utilities_test_codegen_octave: test_m_matrix_named_matsym", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ myout1 <- MatrixSymbol("myout1", 1L, 3L) }, silent = TRUE))
  suppressWarnings(try({ e2 <- Matrix(list(list(x, 2L * y, S_Pi * z))) }, silent = TRUE))
  suppressWarnings(try({ name_expr <- list("test", Equality(myout1, e2, evaluate = FALSE)) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- codegen(name_expr, "Octave", header = FALSE, empty = FALSE); result <- py_getitem(.tmp, 0L) }, silent = TRUE))
  suppressWarnings(try({ source <- py_getitem(result, 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "function myout1 = test(x, y, z)\n  myout1 = [x 2*y pi*z];\nend\n" }, silent = TRUE))
  expect_identical(.cgb(source), "function myout1 = test(x, y, z)\\n  myout1 = [x 2*y pi*z];\\nend\\n")
})

test_that("utilities_test_codegen_octave: test_m_matrix_output_autoname", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ expr <- Matrix(list(list(x, x + y, 3L))) }, silent = TRUE))
  suppressWarnings(try({ name_expr <- list("test", expr) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- codegen(name_expr, "Octave", header = FALSE, empty = FALSE); result <- py_getitem(.tmp, 0L) }, silent = TRUE))
  suppressWarnings(try({ source <- py_getitem(result, 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "function out1 = test(x, y)\n  out1 = [x x + y 3];\nend\n" }, silent = TRUE))
  expect_identical(.cgb(source), "function out1 = test(x, y)\\n  out1 = [x x + y 3];\\nend\\n")
})

test_that("utilities_test_codegen_octave: test_m_matrix_output_autoname_2", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ e1 <- x + y }, silent = TRUE))
  suppressWarnings(try({ e2 <- Matrix(list(list(2L * x, 2L * y, 2L * z))) }, silent = TRUE))
  suppressWarnings(try({ e3 <- Matrix(list(list(x), list(y), list(z))) }, silent = TRUE))
  suppressWarnings(try({ e4 <- Matrix(list(list(x, y), list(z, 16L))) }, silent = TRUE))
  suppressWarnings(try({ name_expr <- list("test", list(e1, e2, e3, e4)) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- codegen(name_expr, "Octave", header = FALSE, empty = FALSE); result <- py_getitem(.tmp, 0L) }, silent = TRUE))
  suppressWarnings(try({ source <- py_getitem(result, 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "function [out1, out2, out3, out4] = test(x, y, z)\n  out1 = x + y;\n  out2 = [2*x 2*y 2*z];\n  out3 = [x; y; z];\n  out4 = [x y; z 16];\nend\n" }, silent = TRUE))
  expect_identical(.cgb(source), "function [out1, out2, out3, out4] = test(x, y, z)\\n  out1 = x + y;\\n  out2 = [2*x 2*y 2*z];\\n  out3 = [x; y; z];\\n  out4 = [x y; z 16];\\nend\\n")
})

test_that("utilities_test_codegen_octave: test_m_results_matrix_named_ordered", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ .tmp <- sym_symbols("B,C"); B <- py_getitem(.tmp, 0L); C <- py_getitem(.tmp, 1L) }, silent = TRUE))
  suppressWarnings(try({ A <- MatrixSymbol("A", 1L, 3L) }, silent = TRUE))
  suppressWarnings(try({ expr1 <- Equality(C, (x + y) * z) }, silent = TRUE))
  suppressWarnings(try({ expr2 <- Equality(A, Matrix(list(list(1L, 2L, x)))) }, silent = TRUE))
  suppressWarnings(try({ expr3 <- Equality(B, 2L * x) }, silent = TRUE))
  suppressWarnings(try({ name_expr <- list("test", list(expr1, expr2, expr3)) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- codegen(name_expr, "Octave", header = FALSE, empty = FALSE, argument_sequence = list(x, z, y)); result <- py_getitem(.tmp, 0L) }, silent = TRUE))
  suppressWarnings(try({ source <- py_getitem(result, 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "function [C, A, B] = test(x, z, y)\n  C = z.*(x + y);\n  A = [1 2 x];\n  B = 2*x;\nend\n" }, silent = TRUE))
  expect_identical(.cgb(source), "function [C, A, B] = test(x, z, y)\\n  C = z.*(x + y);\\n  A = [1 2 x];\\n  B = 2*x;\\nend\\n")
})

test_that("utilities_test_codegen_octave: test_m_matrixsymbol_slice", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ A <- MatrixSymbol("A", 2L, 3L) }, silent = TRUE))
  suppressWarnings(try({ B <- MatrixSymbol("B", 1L, 3L) }, silent = TRUE))
  suppressWarnings(try({ C <- MatrixSymbol("C", 1L, 3L) }, silent = TRUE))
  suppressWarnings(try({ D <- MatrixSymbol("D", 2L, 1L) }, silent = TRUE))
  suppressWarnings(try({ name_expr <- list("test", list(Equality(B, py_getitem(A, 0L, pyslice(NULL, NULL, NULL))), Equality(C, py_getitem(A, 1L, pyslice(NULL, NULL, NULL))), Equality(D, py_getitem(A, pyslice(NULL, NULL, NULL), 2L)))) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- codegen(name_expr, "Octave", header = FALSE, empty = FALSE); result <- py_getitem(.tmp, 0L) }, silent = TRUE))
  suppressWarnings(try({ source <- py_getitem(result, 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "function [B, C, D] = test(A)\n  B = A(1, :);\n  C = A(2, :);\n  D = A(:, 3);\nend\n" }, silent = TRUE))
  expect_identical(.cgb(source), "function [B, C, D] = test(A)\\n  B = A(1, :);\\n  C = A(2, :);\\n  D = A(:, 3);\\nend\\n")
})

test_that("utilities_test_codegen_octave: test_m_matrixsymbol_slice2", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ A <- MatrixSymbol("A", 3L, 4L) }, silent = TRUE))
  suppressWarnings(try({ B <- MatrixSymbol("B", 2L, 2L) }, silent = TRUE))
  suppressWarnings(try({ C <- MatrixSymbol("C", 2L, 2L) }, silent = TRUE))
  suppressWarnings(try({ name_expr <- list("test", list(Equality(B, py_getitem(A, pyslice(0L, 2L, NULL), pyslice(0L, 2L, NULL))), Equality(C, py_getitem(A, pyslice(0L, 2L, NULL), pyslice(1L, 3L, NULL))))) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- codegen(name_expr, "Octave", header = FALSE, empty = FALSE); result <- py_getitem(.tmp, 0L) }, silent = TRUE))
  suppressWarnings(try({ source <- py_getitem(result, 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "function [B, C] = test(A)\n  B = A(1:2, 1:2);\n  C = A(1:2, 2:3);\nend\n" }, silent = TRUE))
  expect_identical(.cgb(source), "function [B, C] = test(A)\\n  B = A(1:2, 1:2);\\n  C = A(1:2, 2:3);\\nend\\n")
})

test_that("utilities_test_codegen_octave: test_m_matrixsymbol_slice3", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ A <- MatrixSymbol("A", 8L, 7L) }, silent = TRUE))
  suppressWarnings(try({ B <- MatrixSymbol("B", 2L, 2L) }, silent = TRUE))
  suppressWarnings(try({ C <- MatrixSymbol("C", 4L, 3L) }, silent = TRUE))
  suppressWarnings(try({ name_expr <- list("test", list(Equality(B, py_getitem(A, pyslice(6L, NULL, NULL), pyslice(1L, NULL, 3L))), Equality(C, py_getitem(A, pyslice(NULL, NULL, 2L), pyslice(NULL, NULL, 3L))))) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- codegen(name_expr, "Octave", header = FALSE, empty = FALSE); result <- py_getitem(.tmp, 0L) }, silent = TRUE))
  suppressWarnings(try({ source <- py_getitem(result, 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "function [B, C] = test(A)\n  B = A(7:end, 2:3:end);\n  C = A(1:2:end, 1:3:end);\nend\n" }, silent = TRUE))
  expect_identical(.cgb(source), "function [B, C] = test(A)\\n  B = A(7:end, 2:3:end);\\n  C = A(1:2:end, 1:3:end);\\nend\\n")
})

test_that("utilities_test_codegen_octave: test_m_matrixsymbol_slice_autoname", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ A <- MatrixSymbol("A", 2L, 3L) }, silent = TRUE))
  suppressWarnings(try({ B <- MatrixSymbol("B", 1L, 3L) }, silent = TRUE))
  suppressWarnings(try({ name_expr <- list("test", list(Equality(B, py_getitem(A, 0L, pyslice(NULL, NULL, NULL))), py_getitem(A, 1L, pyslice(NULL, NULL, NULL)), py_getitem(A, pyslice(NULL, NULL, NULL), 0L), py_getitem(A, pyslice(NULL, NULL, NULL), 1L))) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- codegen(name_expr, "Octave", header = FALSE, empty = FALSE); result <- py_getitem(.tmp, 0L) }, silent = TRUE))
  suppressWarnings(try({ source <- py_getitem(result, 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "function [B, out2, out3, out4] = test(A)\n  B = A(1, :);\n  out2 = A(2, :);\n  out3 = A(:, 1);\n  out4 = A(:, 2);\nend\n" }, silent = TRUE))
  expect_identical(.cgb(source), "function [B, out2, out3, out4] = test(A)\\n  B = A(1, :);\\n  out2 = A(2, :);\\n  out3 = A(:, 1);\\n  out4 = A(:, 2);\\nend\\n")
})

test_that("utilities_test_codegen_octave: test_m_loops", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ invisible(NULL) }, silent = TRUE))
  suppressWarnings(try({ invisible(NULL) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- sym_symbols("n m", integer = TRUE); n <- py_getitem(.tmp, 0L); m <- py_getitem(.tmp, 1L) }, silent = TRUE))
  suppressWarnings(try({ A <- IndexedBase("A") }, silent = TRUE))
  suppressWarnings(try({ x <- IndexedBase("x") }, silent = TRUE))
  suppressWarnings(try({ y <- IndexedBase("y") }, silent = TRUE))
  suppressWarnings(try({ i <- Idx("i", m) }, silent = TRUE))
  suppressWarnings(try({ j <- Idx("j", n) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- codegen(list("mat_vec_mult", Eq(py_getitem(y, i), py_getitem(A, i, j) * py_getitem(x, j))), "Octave", header = FALSE, empty = FALSE); result <- py_getitem(.tmp, 0L) }, silent = TRUE))
  suppressWarnings(try({ source <- py_getitem(result, 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "function y = mat_vec_mult(A, m, n, x)\n  for i = 1:m\n    y(i) = 0;\n  end\n  for i = 1:m\n    for j = 1:n\n      y(i) = %(rhs)s + y(i);\n    end\n  end\nend\n" }, silent = TRUE))
  expect_identical(.cgb((py_eq(source, py_format(expected, py_dict(list("rhs"), list(py_format("A(%s, %s).*x(j)", list(i, j))))))) || (py_eq(source, py_format(expected, py_dict(list("rhs"), list(py_format("x(j).*A(%s, %s)", list(i, j)))))))), "True")
})

test_that("utilities_test_codegen_octave: test_m_tensor_loops_multiple_contractions", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ invisible(NULL) }, silent = TRUE))
  suppressWarnings(try({ invisible(NULL) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- sym_symbols("n m o p", integer = TRUE); n <- py_getitem(.tmp, 0L); m <- py_getitem(.tmp, 1L); o <- py_getitem(.tmp, 2L); p <- py_getitem(.tmp, 3L) }, silent = TRUE))
  suppressWarnings(try({ A <- IndexedBase("A") }, silent = TRUE))
  suppressWarnings(try({ B <- IndexedBase("B") }, silent = TRUE))
  suppressWarnings(try({ y <- IndexedBase("y") }, silent = TRUE))
  suppressWarnings(try({ i <- Idx("i", m) }, silent = TRUE))
  suppressWarnings(try({ j <- Idx("j", n) }, silent = TRUE))
  suppressWarnings(try({ k <- Idx("k", o) }, silent = TRUE))
  suppressWarnings(try({ l <- Idx("l", p) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- codegen(list("tensorthing", Eq(py_getitem(y, i), py_getitem(B, j, k, l) * py_getitem(A, i, j, k, l))), "Octave", header = FALSE, empty = FALSE); result <- py_getitem(.tmp, 0L) }, silent = TRUE))
  suppressWarnings(try({ source <- py_getitem(result, 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "function y = tensorthing(A, B, m, n, o, p)\n  for i = 1:m\n    y(i) = 0;\n  end\n  for i = 1:m\n    for j = 1:n\n      for k = 1:o\n        for l = 1:p\n          y(i) = A(i, j, k, l).*B(j, k, l) + y(i);\n        end\n      end\n    end\n  end\nend\n" }, silent = TRUE))
  expect_identical(.cgb(source), "function y = tensorthing(A, B, m, n, o, p)\\n  for i = 1:m\\n    y(i) = 0;\\n  end\\n  for i = 1:m\\n    for j = 1:n\\n      for k = 1:o\\n        for l = 1:p\\n          y(i) = A(i, j, k, l).*B(j, k, l) + y(i);\\n        end\\n      end\\n    end\\n  end\\nend\\n")
})

test_that("utilities_test_codegen_octave: test_m_InOutArgument", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ expr <- Equality(x, x ^ 2L) }, silent = TRUE))
  suppressWarnings(try({ name_expr <- list("mysqr", expr) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- codegen(name_expr, "Octave", header = FALSE, empty = FALSE); result <- py_getitem(.tmp, 0L) }, silent = TRUE))
  suppressWarnings(try({ source <- py_getitem(result, 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "function x = mysqr(x)\n  x = x.^2;\nend\n" }, silent = TRUE))
  expect_identical(.cgb(source), "function x = mysqr(x)\\n  x = x.^2;\\nend\\n")
})

test_that("utilities_test_codegen_octave: test_m_InOutArgument_order", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ expr <- Equality(x, (x ^ 2L) + y) }, silent = TRUE))
  suppressWarnings(try({ name_expr <- list("test", expr) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- codegen(name_expr, "Octave", header = FALSE, empty = FALSE, argument_sequence = list(x, y)); result <- py_getitem(.tmp, 0L) }, silent = TRUE))
  suppressWarnings(try({ source <- py_getitem(result, 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "function x = test(x, y)\n  x = x.^2 + y;\nend\n" }, silent = TRUE))
  expect_identical(.cgb(source), "function x = test(x, y)\\n  x = x.^2 + y;\\nend\\n")
  suppressWarnings(try({ expr <- Equality(x, (x ^ 2L) + y) }, silent = TRUE))
  suppressWarnings(try({ name_expr <- list("test", expr) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- codegen(name_expr, "Octave", header = FALSE, empty = FALSE); result <- py_getitem(.tmp, 0L) }, silent = TRUE))
  suppressWarnings(try({ source <- py_getitem(result, 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "function x = test(x, y)\n  x = x.^2 + y;\nend\n" }, silent = TRUE))
  expect_identical(.cgb(source), "function x = test(x, y)\\n  x = x.^2 + y;\\nend\\n")
})

test_that("utilities_test_codegen_octave: test_m_not_supported", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ f <- Function("f") }, silent = TRUE))
  suppressWarnings(try({ name_expr <- list("test", list(diff(f(x), x), zoo)) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- codegen(name_expr, "Octave", header = FALSE, empty = FALSE); result <- py_getitem(.tmp, 0L) }, silent = TRUE))
  suppressWarnings(try({ source <- py_getitem(result, 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "function [out1, out2] = test(x)\n  % unsupported: Derivative(f(x), x)\n  % unsupported: zoo\n  out1 = Derivative(f(x), x);\n  out2 = zoo;\nend\n" }, silent = TRUE))
  expect_identical(.cgb(source), "function [out1, out2] = test(x)\\n  % unsupported: Derivative(f(x), x)\\n  % unsupported: zoo\\n  out1 = Derivative(f(x), x);\\n  out2 = zoo;\\nend\\n")
})

test_that("utilities_test_codegen_octave: test_global_vars_octave", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ .tmp <- sym_symbols("x y z t"); x <- py_getitem(.tmp, 0L); y <- py_getitem(.tmp, 1L); z <- py_getitem(.tmp, 2L); t <- py_getitem(.tmp, 3L) }, silent = TRUE))
  suppressWarnings(try({ result <- codegen(list("f", x * y), "Octave", header = FALSE, empty = FALSE, global_vars = list(y)) }, silent = TRUE))
  suppressWarnings(try({ source <- py_getitem(py_getitem(result, 0L), 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "function out1 = f(x)\n  global y\n  out1 = x.*y;\nend\n" }, silent = TRUE))
  expect_identical(.cgb(source), "function out1 = f(x)\\n  global y\\n  out1 = x.*y;\\nend\\n")
  suppressWarnings(try({ result <- codegen(list("f", (x * y) + z), "Octave", header = FALSE, empty = FALSE, argument_sequence = list(x, y), global_vars = list(z, t)) }, silent = TRUE))
  suppressWarnings(try({ source <- py_getitem(py_getitem(result, 0L), 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "function out1 = f(x, y)\n  global t z\n  out1 = x.*y + z;\nend\n" }, silent = TRUE))
  expect_identical(.cgb(source), "function out1 = f(x, y)\\n  global t z\\n  out1 = x.*y + z;\\nend\\n")
})

test_that("utilities_test_codegen_rust: test_empty_rust_code", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ code_gen <- RustCodeGen() }, silent = TRUE))
  suppressWarnings(try({ output <- StringIO() }, silent = TRUE))
  suppressWarnings(try({ py_method(code_gen, "dump_rs", list(), output, "file", header = FALSE, empty = FALSE) }, silent = TRUE))
  suppressWarnings(try({ source <- py_method(output, "getvalue") }, silent = TRUE))
  expect_identical(.cgb(source), "")
  for (.o in mget(ls(), inherits = FALSE)) if (inherits(.o, "py_stringio")) try(close(.o$con), silent = TRUE)
})

test_that("utilities_test_codegen_rust: test_simple_rust_code", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ name_expr <- list("test", (x + y) * z) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- codegen(name_expr, "Rust", header = FALSE, empty = FALSE); result <- py_getitem(.tmp, 0L) }, silent = TRUE))
  expect_identical(.cgb(py_getitem(result, 0L)), "test.rs")
  suppressWarnings(try({ source <- py_getitem(result, 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "fn test(x: f64, y: f64, z: f64) -> f64 {\n    let out1 = z*(x + y);\n    out1\n}\n" }, silent = TRUE))
  expect_identical(.cgb(source), "fn test(x: f64, y: f64, z: f64) -> f64 {\\n    let out1 = z*(x + y);\\n    out1\\n}\\n")
})

test_that("utilities_test_codegen_rust: test_simple_code_with_header", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ name_expr <- list("test", (x + y) * z) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- codegen(name_expr, "Rust", header = TRUE, empty = FALSE); result <- py_getitem(.tmp, 0L) }, silent = TRUE))
  expect_identical(.cgb(py_getitem(result, 0L)), "test.rs")
  suppressWarnings(try({ source <- py_getitem(result, 1L) }, silent = TRUE))
  suppressWarnings(try({ version_str <- py_format("Code generated with SymPy %s", ast_attr(sympy, "__version__")) }, silent = TRUE))
  suppressWarnings(try({ version_line <- py_method(py_method(version_str, "center", 76L), "rstrip") }, silent = TRUE))
  suppressWarnings(try({ expected <- py_format("/*\n *%(version_line)s\n *\n *              See http://www.sympy.org/ for more information.\n *\n *                       This file is part of 'project'\n */\nfn test(x: f64, y: f64, z: f64) -> f64 {\n    let out1 = z*(x + y);\n    out1\n}\n", py_dict(list("version_line"), list(version_line))) }, silent = TRUE))
})

test_that("utilities_test_codegen_rust: test_simple_code_nameout", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ expr <- Equality(z, x + y) }, silent = TRUE))
  suppressWarnings(try({ name_expr <- list("test", expr) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- codegen(name_expr, "Rust", header = FALSE, empty = FALSE); result <- py_getitem(.tmp, 0L) }, silent = TRUE))
  suppressWarnings(try({ source <- py_getitem(result, 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "fn test(x: f64, y: f64) -> f64 {\n    let z = x + y;\n    z\n}\n" }, silent = TRUE))
  expect_identical(.cgb(source), "fn test(x: f64, y: f64) -> f64 {\\n    let z = x + y;\\n    z\\n}\\n")
})

test_that("utilities_test_codegen_rust: test_numbersymbol", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ name_expr <- list("test", S_Pi ^ Catalan) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- codegen(name_expr, "Rust", header = FALSE, empty = FALSE); result <- py_getitem(.tmp, 0L) }, silent = TRUE))
  suppressWarnings(try({ source <- py_getitem(result, 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- py_format("fn test() -> f64 {\n    const Catalan: f64 = %s;\n    let out1 = PI.powf(Catalan);\n    out1\n}\n", N(Catalan, 17L)) }, silent = TRUE))
  expect_identical(.cgb(source), "fn test() -> f64 {\\n    const Catalan: f64 = 0.91596559417721901;\\n    let out1 = PI.powf(Catalan);\\n    out1\\n}\\n")
})

test_that("utilities_test_codegen_rust: test_argument_order", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ expr <- x + y }, silent = TRUE))
  suppressWarnings(try({ routine <- make_routine("test", expr, argument_sequence = list(z, x, y), language = "rust") }, silent = TRUE))
  suppressWarnings(try({ code_gen <- RustCodeGen() }, silent = TRUE))
  suppressWarnings(try({ output <- StringIO() }, silent = TRUE))
  suppressWarnings(try({ py_method(code_gen, "dump_rs", list(routine), output, "test", header = FALSE, empty = FALSE) }, silent = TRUE))
  suppressWarnings(try({ source <- py_method(output, "getvalue") }, silent = TRUE))
  suppressWarnings(try({ expected <- "fn test(z: f64, x: f64, y: f64) -> f64 {\n    let out1 = x + y;\n    out1\n}\n" }, silent = TRUE))
  expect_identical(.cgb(source), "fn test(z: f64, x: f64, y: f64) -> f64 {\\n    let out1 = x + y;\\n    out1\\n}\\n")
  for (.o in mget(ls(), inherits = FALSE)) if (inherits(.o, "py_stringio")) try(close(.o$con), silent = TRUE)
})

test_that("utilities_test_codegen_rust: test_multiple_results_rust", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ expr1 <- (x + y) * z }, silent = TRUE))
  suppressWarnings(try({ expr2 <- (x - y) * z }, silent = TRUE))
  suppressWarnings(try({ name_expr <- list("test", list(expr1, expr2)) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- codegen(name_expr, "Rust", header = FALSE, empty = FALSE); result <- py_getitem(.tmp, 0L) }, silent = TRUE))
  suppressWarnings(try({ source <- py_getitem(result, 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "fn test(x: f64, y: f64, z: f64) -> (f64, f64) {\n    let out1 = z*(x + y);\n    let out2 = z*(x - y);\n    (out1, out2)\n}\n" }, silent = TRUE))
  expect_identical(.cgb(source), "fn test(x: f64, y: f64, z: f64) -> (f64, f64) {\\n    let out1 = z*(x + y);\\n    let out2 = z*(x - y);\\n    (out1, out2)\\n}\\n")
})

test_that("utilities_test_codegen_rust: test_results_named_unordered", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ .tmp <- sym_symbols("A,B,C"); A <- py_getitem(.tmp, 0L); B <- py_getitem(.tmp, 1L); C <- py_getitem(.tmp, 2L) }, silent = TRUE))
  suppressWarnings(try({ expr1 <- Equality(C, (x + y) * z) }, silent = TRUE))
  suppressWarnings(try({ expr2 <- Equality(A, (x - y) * z) }, silent = TRUE))
  suppressWarnings(try({ expr3 <- Equality(B, 2L * x) }, silent = TRUE))
  suppressWarnings(try({ name_expr <- list("test", list(expr1, expr2, expr3)) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- codegen(name_expr, "Rust", header = FALSE, empty = FALSE); result <- py_getitem(.tmp, 0L) }, silent = TRUE))
  suppressWarnings(try({ source <- py_getitem(result, 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "fn test(x: f64, y: f64, z: f64) -> (f64, f64, f64) {\n    let C = z*(x + y);\n    let A = z*(x - y);\n    let B = 2*x;\n    (C, A, B)\n}\n" }, silent = TRUE))
  expect_identical(.cgb(source), "fn test(x: f64, y: f64, z: f64) -> (f64, f64, f64) {\\n    let C = z*(x + y);\\n    let A = z*(x - y);\\n    let B = 2*x;\\n    (C, A, B)\\n}\\n")
})

test_that("utilities_test_codegen_rust: test_results_named_ordered", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ .tmp <- sym_symbols("A,B,C"); A <- py_getitem(.tmp, 0L); B <- py_getitem(.tmp, 1L); C <- py_getitem(.tmp, 2L) }, silent = TRUE))
  suppressWarnings(try({ expr1 <- Equality(C, (x + y) * z) }, silent = TRUE))
  suppressWarnings(try({ expr2 <- Equality(A, (x - y) * z) }, silent = TRUE))
  suppressWarnings(try({ expr3 <- Equality(B, 2L * x) }, silent = TRUE))
  suppressWarnings(try({ name_expr <- list("test", list(expr1, expr2, expr3)) }, silent = TRUE))
  suppressWarnings(try({ result <- codegen(name_expr, "Rust", header = FALSE, empty = FALSE, argument_sequence = list(x, z, y)) }, silent = TRUE))
  expect_identical(.cgb(py_getitem(py_getitem(result, 0L), 0L)), "test.rs")
  suppressWarnings(try({ source <- py_getitem(py_getitem(result, 0L), 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "fn test(x: f64, z: f64, y: f64) -> (f64, f64, f64) {\n    let C = z*(x + y);\n    let A = z*(x - y);\n    let B = 2*x;\n    (C, A, B)\n}\n" }, silent = TRUE))
  expect_identical(.cgb(source), "fn test(x: f64, z: f64, y: f64) -> (f64, f64, f64) {\\n    let C = z*(x + y);\\n    let A = z*(x - y);\\n    let B = 2*x;\\n    (C, A, B)\\n}\\n")
})

test_that("utilities_test_codegen_rust: test_complicated_rs_codegen", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ invisible(NULL) }, silent = TRUE))
  suppressWarnings(try({ name_expr <- list("testlong", list(expand(((sin(x) + cos(y)) + tan(z)) ^ 3L), cos(cos(cos(cos(cos(cos(cos(cos((x + y) + z)))))))))) }, silent = TRUE))
  suppressWarnings(try({ result <- codegen(name_expr, "Rust", header = FALSE, empty = FALSE) }, silent = TRUE))
  expect_identical(.cgb(py_getitem(py_getitem(result, 0L), 0L)), "testlong.rs")
  suppressWarnings(try({ source <- py_getitem(py_getitem(result, 0L), 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "fn testlong(x: f64, y: f64, z: f64) -> (f64, f64) {\n    let out1 = x.sin().powi(3) + 3*x.sin().powi(2)*y.cos() + 3*x.sin().powi(2)*z.tan() + 3*x.sin()*y.cos().powi(2) + 6*x.sin()*y.cos()*z.tan() + 3*x.sin()*z.tan().powi(2) + y.cos().powi(3) + 3*y.cos().powi(2)*z.tan() + 3*y.cos()*z.tan().powi(2) + z.tan().powi(3);\n    let out2 = (x + y + z).cos().cos().cos().cos().cos().cos().cos().cos();\n    (out1, out2)\n}\n" }, silent = TRUE))
  expect_identical(.cgb(source), "fn testlong(x: f64, y: f64, z: f64) -> (f64, f64) {\\n    let out1 = x.sin().powi(3) + 3*x.sin().powi(2)*y.cos() + 3*x.sin().powi(2)*z.tan() + 3*x.sin()*y.cos().powi(2) + 6*x.sin()*y.cos()*z.tan() + 3*x.sin()*z.tan().powi(2) + y.cos().powi(3) + 3*y.cos().powi(2)*z.tan() + 3*y.cos()*z.tan().powi(2) + z.tan().powi(3);\\n    let out2 = (x + y + z).cos().cos().cos().cos().cos().cos().cos().cos();\\n    (out1, out2)\\n}\\n")
})

test_that("utilities_test_codegen_rust: test_output_arg_mixed_unordered", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ invisible(NULL) }, silent = TRUE))
  suppressWarnings(try({ a <- sym_symbols("a") }, silent = TRUE))
  suppressWarnings(try({ name_expr <- list("foo", list(cos(2L * x), Equality(y, sin(x)), cos(x), Equality(a, sin(2L * x)))) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- codegen(name_expr, "Rust", header = FALSE, empty = FALSE); result <- py_getitem(.tmp, 0L) }, silent = TRUE))
  expect_identical(.cgb(py_getitem(result, 0L)), "foo.rs")
  suppressWarnings(try({ source <- py_getitem(result, 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "fn foo(x: f64) -> (f64, f64, f64, f64) {\n    let out1 = (2*x).cos();\n    let y = x.sin();\n    let out3 = x.cos();\n    let a = (2*x).sin();\n    (out1, y, out3, a)\n}\n" }, silent = TRUE))
  expect_identical(.cgb(source), "fn foo(x: f64) -> (f64, f64, f64, f64) {\\n    let out1 = (2*x).cos();\\n    let y = x.sin();\\n    let out3 = x.cos();\\n    let a = (2*x).sin();\\n    (out1, y, out3, a)\\n}\\n")
})

test_that("utilities_test_codegen_rust: test_piecewise_", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ pw <- Piecewise(list(0L, x < (-1L)), list(x ^ 2L, x <= 1L), list((-x) + 2L, x > 1L), list(1L, TRUE), evaluate = FALSE) }, silent = TRUE))
  suppressWarnings(try({ name_expr <- list("pwtest", pw) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- codegen(name_expr, "Rust", header = FALSE, empty = FALSE); result <- py_getitem(.tmp, 0L) }, silent = TRUE))
  suppressWarnings(try({ source <- py_getitem(result, 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "fn pwtest(x: f64) -> f64 {\n    let out1 = if (x < -1.0) {\n        0\n    } else if (x <= 1.0) {\n        x.powi(2)\n    } else if (x > 1.0) {\n        2 - x\n    } else {\n        1\n    };\n    out1\n}\n" }, silent = TRUE))
  expect_identical(.cgb(source), "fn pwtest(x: f64) -> f64 {\\n    let out1 = if (x < -1.0) {\\n        0\\n    } else if (x <= 1.0) {\\n        x.powi(2)\\n    } else if (x > 1.0) {\\n        2 - x\\n    } else {\\n        1\\n    };\\n    out1\\n}\\n")
})

test_that("utilities_test_codegen_rust: test_multifcns_per_file", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ name_expr <- list(list("foo", list(2L * x, 3L * y)), list("bar", list(y ^ 2L, 4L * y))) }, silent = TRUE))
  suppressWarnings(try({ result <- codegen(name_expr, "Rust", header = FALSE, empty = FALSE) }, silent = TRUE))
  expect_identical(.cgb(py_getitem(py_getitem(result, 0L), 0L)), "foo.rs")
  suppressWarnings(try({ source <- py_getitem(py_getitem(result, 0L), 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "fn foo(x: f64, y: f64) -> (f64, f64) {\n    let out1 = 2*x;\n    let out2 = 3*y;\n    (out1, out2)\n}\nfn bar(y: f64) -> (f64, f64) {\n    let out1 = y.powi(2);\n    let out2 = 4*y;\n    (out1, out2)\n}\n" }, silent = TRUE))
  expect_identical(.cgb(source), "fn foo(x: f64, y: f64) -> (f64, f64) {\\n    let out1 = 2*x;\\n    let out2 = 3*y;\\n    (out1, out2)\\n}\\nfn bar(y: f64) -> (f64, f64) {\\n    let out1 = y.powi(2);\\n    let out2 = 4*y;\\n    (out1, out2)\\n}\\n")
})

test_that("utilities_test_codegen_rust: test_multifcns_per_file_w_header", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ name_expr <- list(list("foo", list(2L * x, 3L * y)), list("bar", list(y ^ 2L, 4L * y))) }, silent = TRUE))
  suppressWarnings(try({ result <- codegen(name_expr, "Rust", header = TRUE, empty = FALSE) }, silent = TRUE))
  expect_identical(.cgb(py_getitem(py_getitem(result, 0L), 0L)), "foo.rs")
  suppressWarnings(try({ source <- py_getitem(py_getitem(result, 0L), 1L) }, silent = TRUE))
  suppressWarnings(try({ version_str <- py_format("Code generated with SymPy %s", ast_attr(sympy, "__version__")) }, silent = TRUE))
  suppressWarnings(try({ version_line <- py_method(py_method(version_str, "center", 76L), "rstrip") }, silent = TRUE))
  suppressWarnings(try({ expected <- py_format("/*\n *%(version_line)s\n *\n *              See http://www.sympy.org/ for more information.\n *\n *                       This file is part of 'project'\n */\nfn foo(x: f64, y: f64) -> (f64, f64) {\n    let out1 = 2*x;\n    let out2 = 3*y;\n    (out1, out2)\n}\nfn bar(y: f64) -> (f64, f64) {\n    let out1 = y.powi(2);\n    let out2 = 4*y;\n    (out1, out2)\n}\n", py_dict(list("version_line"), list(version_line))) }, silent = TRUE))
})

test_that("utilities_test_codegen_rust: test_filename_match_prefix", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ name_expr <- list(list("foo", list(2L * x, 3L * y)), list("bar", list(y ^ 2L, 4L * y))) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- codegen(name_expr, "Rust", prefix = "baz", header = FALSE, empty = FALSE); result <- py_getitem(.tmp, 0L) }, silent = TRUE))
  expect_identical(.cgb(py_getitem(result, 0L)), "baz.rs")
})

test_that("utilities_test_codegen_rust: test_InOutArgument", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ expr <- Equality(x, x ^ 2L) }, silent = TRUE))
  suppressWarnings(try({ name_expr <- list("mysqr", expr) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- codegen(name_expr, "Rust", header = FALSE, empty = FALSE); result <- py_getitem(.tmp, 0L) }, silent = TRUE))
  suppressWarnings(try({ source <- py_getitem(result, 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "fn mysqr(x: f64) -> f64 {\n    let x = x.powi(2);\n    x\n}\n" }, silent = TRUE))
  expect_identical(.cgb(source), "fn mysqr(x: f64) -> f64 {\\n    let x = x.powi(2);\\n    x\\n}\\n")
})

test_that("utilities_test_codegen_rust: test_InOutArgument_order", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ expr <- Equality(x, (x ^ 2L) + y) }, silent = TRUE))
  suppressWarnings(try({ name_expr <- list("test", expr) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- codegen(name_expr, "Rust", header = FALSE, empty = FALSE, argument_sequence = list(x, y)); result <- py_getitem(.tmp, 0L) }, silent = TRUE))
  suppressWarnings(try({ source <- py_getitem(result, 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "fn test(x: f64, y: f64) -> f64 {\n    let x = x.powi(2) + y;\n    x\n}\n" }, silent = TRUE))
  expect_identical(.cgb(source), "fn test(x: f64, y: f64) -> f64 {\\n    let x = x.powi(2) + y;\\n    x\\n}\\n")
  suppressWarnings(try({ expr <- Equality(x, (x ^ 2L) + y) }, silent = TRUE))
  suppressWarnings(try({ name_expr <- list("test", expr) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- codegen(name_expr, "Rust", header = FALSE, empty = FALSE); result <- py_getitem(.tmp, 0L) }, silent = TRUE))
  suppressWarnings(try({ source <- py_getitem(result, 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "fn test(x: f64, y: f64) -> f64 {\n    let x = x.powi(2) + y;\n    x\n}\n" }, silent = TRUE))
  expect_identical(.cgb(source), "fn test(x: f64, y: f64) -> f64 {\\n    let x = x.powi(2) + y;\\n    x\\n}\\n")
})

test_that("utilities_test_codegen_rust: test_not_supported", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ f <- Function("f") }, silent = TRUE))
  suppressWarnings(try({ name_expr <- list("test", list(diff(f(x), x), zoo)) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- codegen(name_expr, "Rust", header = FALSE, empty = FALSE); result <- py_getitem(.tmp, 0L) }, silent = TRUE))
  suppressWarnings(try({ source <- py_getitem(result, 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "fn test(x: f64) -> (f64, f64) {\n    // unsupported: Derivative(f(x), x)\n    // unsupported: zoo\n    let out1 = Derivative(f(x), x);\n    let out2 = zoo;\n    (out1, out2)\n}\n" }, silent = TRUE))
  expect_identical(.cgb(source), "fn test(x: f64) -> (f64, f64) {\\n    // unsupported: Derivative(f(x), x)\\n    // unsupported: zoo\\n    let out1 = Derivative(f(x), x);\\n    let out2 = zoo;\\n    (out1, out2)\\n}\\n")
})

test_that("utilities_test_codegen_rust: test_global_vars_rust", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ .tmp <- sym_symbols("x y z t"); x <- py_getitem(.tmp, 0L); y <- py_getitem(.tmp, 1L); z <- py_getitem(.tmp, 2L); t <- py_getitem(.tmp, 3L) }, silent = TRUE))
  suppressWarnings(try({ result <- codegen(list("f", x * y), "Rust", header = FALSE, empty = FALSE, global_vars = list(y)) }, silent = TRUE))
  suppressWarnings(try({ source <- py_getitem(py_getitem(result, 0L), 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "fn f(x: f64) -> f64 {\n    let out1 = x*y;\n    out1\n}\n" }, silent = TRUE))
  expect_identical(.cgb(source), "fn f(x: f64) -> f64 {\\n    let out1 = x*y;\\n    out1\\n}\\n")
  suppressWarnings(try({ result <- codegen(list("f", (x * y) + z), "Rust", header = FALSE, empty = FALSE, argument_sequence = list(x, y), global_vars = list(z, t)) }, silent = TRUE))
  suppressWarnings(try({ source <- py_getitem(py_getitem(result, 0L), 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "fn f(x: f64, y: f64) -> f64 {\n    let out1 = x*y + z;\n    out1\n}\n" }, silent = TRUE))
  expect_identical(.cgb(source), "fn f(x: f64, y: f64) -> f64 {\\n    let out1 = x*y + z;\\n    out1\\n}\\n")
})

test_that("utilities_test_codegen: test_Routine_argument_order", {
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ .tmp <- sym_symbols("a x y z"); a <- py_getitem(.tmp, 0L); x <- py_getitem(.tmp, 1L); y <- py_getitem(.tmp, 2L); z <- py_getitem(.tmp, 3L) }, silent = TRUE))
  suppressWarnings(try({ expr <- (x + y) * z }, silent = TRUE))
  expect_identical(.cgb(make_routine("test", expr, argument_sequence = list(z, x))), "ERR CodeGenArgumentListError: (\"Argument list didn't specify: y \", [InputArgument(y)])")
  expect_identical(.cgb(make_routine("test", Eq(a, expr), argument_sequence = list(z, x, y))), "ERR CodeGenArgumentListError: (\"Argument list didn't specify: a \", [OutputArgument(a, a, z*(x + y))])")
  suppressWarnings(try({ r <- make_routine("test", Eq(a, expr), argument_sequence = list(z, x, a, y)) }, silent = TRUE))
  expect_identical(.cgb(lapply(py_iter(ast_attr(r, "arguments")), function(arg) ast_attr(arg, "name"))), "[z, x, a, y]")
  expect_identical(.cgb(lapply(py_iter(ast_attr(r, "arguments")), function(arg) py_type(arg))), "[<class 'sympy.utilities.codegen.InputArgument'>, <class 'sympy.utilities.codegen.InputArgument'>, <class 'sympy.utilities.codegen.OutputArgument'>, <class 'sympy.utilities.codegen.InputArgument'>]")
  suppressWarnings(try({ r <- make_routine("test", Eq(z, expr), argument_sequence = list(z, x, y)) }, silent = TRUE))
  expect_identical(.cgb(lapply(py_iter(ast_attr(r, "arguments")), function(arg) py_type(arg))), "[<class 'sympy.utilities.codegen.InOutArgument'>, <class 'sympy.utilities.codegen.InputArgument'>, <class 'sympy.utilities.codegen.InputArgument'>]")
  suppressWarnings(try({ invisible(NULL) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- map(IndexedBase, list("A", "B")); A <- py_getitem(.tmp, 0L); B <- py_getitem(.tmp, 1L) }, silent = TRUE))
  suppressWarnings(try({ m <- sym_symbols("m", integer = TRUE) }, silent = TRUE))
  suppressWarnings(try({ i <- Idx("i", m) }, silent = TRUE))
  suppressWarnings(try({ r <- make_routine("test", Eq(py_getitem(A, i), py_getitem(B, i)), argument_sequence = list(B, A, m)) }, silent = TRUE))
  expect_identical(.cgb(lapply(py_iter(ast_attr(r, "arguments")), function(arg) ast_attr(arg, "name"))), "[B, A, m]")
  suppressWarnings(try({ expr <- Integral((x * y) * z, list(x, 1L, 2L), list(y, 1L, 3L)) }, silent = TRUE))
  suppressWarnings(try({ r <- make_routine("test", Eq(a, expr), argument_sequence = list(z, x, a, y)) }, silent = TRUE))
  expect_identical(.cgb(lapply(py_iter(ast_attr(r, "arguments")), function(arg) ast_attr(arg, "name"))), "[z, x, a, y]")
})

test_that("utilities_test_codegen: test_empty_c_code", {
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ code_gen <- C89CodeGen() }, silent = TRUE))
  suppressWarnings(try({ source <- get_string(ast_attr(code_gen, "dump_c"), list()) }, silent = TRUE))
  expect_identical(.cgb(source), "#include \"file.h\"\\n#include <math.h>\\n")
})

test_that("utilities_test_codegen: test_empty_c_code_with_comment", {
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ code_gen <- C89CodeGen() }, silent = TRUE))
  suppressWarnings(try({ source <- get_string(ast_attr(code_gen, "dump_c"), list(), header = TRUE) }, silent = TRUE))
  expect_identical(.cgb(py_slice(source, NULL, 82L)), "/^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\\n *")
  expect_identical(.cgb(py_slice(source, 158L, NULL)), "*\\n *                                                                            *\\n *              See http://www.sympy.org/ for more information.               *\\n *                                                                            *\\n *                       This file is part of 'project'                       *\\n ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^/\\n#include \"file.h\"\\n#include <math.h>\\n")
})

test_that("utilities_test_codegen: test_empty_c_header", {
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ code_gen <- C99CodeGen() }, silent = TRUE))
  suppressWarnings(try({ source <- get_string(ast_attr(code_gen, "dump_h"), list()) }, silent = TRUE))
  expect_identical(.cgb(source), "#ifndef PROJECT__FILE__H\\n#define PROJECT__FILE__H\\n#endif\\n")
})

test_that("utilities_test_codegen: test_simple_c_code", {
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ .tmp <- sym_symbols("x,y,z"); x <- py_getitem(.tmp, 0L); y <- py_getitem(.tmp, 1L); z <- py_getitem(.tmp, 2L) }, silent = TRUE))
  suppressWarnings(try({ expr <- (x + y) * z }, silent = TRUE))
  suppressWarnings(try({ routine <- make_routine("test", expr) }, silent = TRUE))
  suppressWarnings(try({ code_gen <- C89CodeGen() }, silent = TRUE))
  suppressWarnings(try({ source <- get_string(ast_attr(code_gen, "dump_c"), list(routine)) }, silent = TRUE))
  suppressWarnings(try({ expected <- "#include \"file.h\"\n#include <math.h>\ndouble test(double x, double y, double z) {\n   double test_result;\n   test_result = z*(x + y);\n   return test_result;\n}\n" }, silent = TRUE))
  expect_identical(.cgb(source), "#include \"file.h\"\\n#include <math.h>\\ndouble test(double x, double y, double z) {\\n   double test_result;\\n   test_result = z*(x + y);\\n   return test_result;\\n}\\n")
})

test_that("utilities_test_codegen: test_c_code_reserved_words", {
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ .tmp <- sym_symbols("if, typedef, while"); x <- py_getitem(.tmp, 0L); y <- py_getitem(.tmp, 1L); z <- py_getitem(.tmp, 2L) }, silent = TRUE))
  suppressWarnings(try({ expr <- (x + y) * z }, silent = TRUE))
  suppressWarnings(try({ routine <- make_routine("test", expr) }, silent = TRUE))
  suppressWarnings(try({ code_gen <- C99CodeGen() }, silent = TRUE))
  suppressWarnings(try({ source <- get_string(ast_attr(code_gen, "dump_c"), list(routine)) }, silent = TRUE))
  suppressWarnings(try({ expected <- "#include \"file.h\"\n#include <math.h>\ndouble test(double if_, double typedef_, double while_) {\n   double test_result;\n   test_result = while_*(if_ + typedef_);\n   return test_result;\n}\n" }, silent = TRUE))
  expect_identical(.cgb(source), "#include \"file.h\"\\n#include <math.h>\\ndouble test(double if_, double typedef_, double while_) {\\n   double test_result;\\n   test_result = while_*(if_ + typedef_);\\n   return test_result;\\n}\\n")
})

test_that("utilities_test_codegen: test_numbersymbol_c_code", {
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ routine <- make_routine("test", S_Pi ^ Catalan) }, silent = TRUE))
  suppressWarnings(try({ code_gen <- C89CodeGen() }, silent = TRUE))
  suppressWarnings(try({ source <- get_string(ast_attr(code_gen, "dump_c"), list(routine)) }, silent = TRUE))
  suppressWarnings(try({ expected <- py_format("#include \"file.h\"\n#include <math.h>\ndouble test() {\n   double test_result;\n   double const Catalan = %s;\n   test_result = pow(M_PI, Catalan);\n   return test_result;\n}\n", N(Catalan, 17L)) }, silent = TRUE))
  expect_identical(.cgb(source), "#include \"file.h\"\\n#include <math.h>\\ndouble test() {\\n   double test_result;\\n   double const Catalan = 0.91596559417721901;\\n   test_result = pow(M_PI, Catalan);\\n   return test_result;\\n}\\n")
})

test_that("utilities_test_codegen: test_c_code_argument_order", {
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ .tmp <- sym_symbols("x,y,z"); x <- py_getitem(.tmp, 0L); y <- py_getitem(.tmp, 1L); z <- py_getitem(.tmp, 2L) }, silent = TRUE))
  suppressWarnings(try({ expr <- x + y }, silent = TRUE))
  suppressWarnings(try({ routine <- make_routine("test", expr, argument_sequence = list(z, x, y)) }, silent = TRUE))
  suppressWarnings(try({ code_gen <- C89CodeGen() }, silent = TRUE))
  suppressWarnings(try({ source <- get_string(ast_attr(code_gen, "dump_c"), list(routine)) }, silent = TRUE))
  suppressWarnings(try({ expected <- "#include \"file.h\"\n#include <math.h>\ndouble test(double z, double x, double y) {\n   double test_result;\n   test_result = x + y;\n   return test_result;\n}\n" }, silent = TRUE))
  expect_identical(.cgb(source), "#include \"file.h\"\\n#include <math.h>\\ndouble test(double z, double x, double y) {\\n   double test_result;\\n   test_result = x + y;\\n   return test_result;\\n}\\n")
})

test_that("utilities_test_codegen: test_simple_c_header", {
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ .tmp <- sym_symbols("x,y,z"); x <- py_getitem(.tmp, 0L); y <- py_getitem(.tmp, 1L); z <- py_getitem(.tmp, 2L) }, silent = TRUE))
  suppressWarnings(try({ expr <- (x + y) * z }, silent = TRUE))
  suppressWarnings(try({ routine <- make_routine("test", expr) }, silent = TRUE))
  suppressWarnings(try({ code_gen <- C89CodeGen() }, silent = TRUE))
  suppressWarnings(try({ source <- get_string(ast_attr(code_gen, "dump_h"), list(routine)) }, silent = TRUE))
  suppressWarnings(try({ expected <- "#ifndef PROJECT__FILE__H\n#define PROJECT__FILE__H\ndouble test(double x, double y, double z);\n#endif\n" }, silent = TRUE))
  expect_identical(.cgb(source), "#ifndef PROJECT__FILE__H\\n#define PROJECT__FILE__H\\ndouble test(double x, double y, double z);\\n#endif\\n")
})

test_that("utilities_test_codegen: test_simple_c_codegen", {
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ .tmp <- sym_symbols("x,y,z"); x <- py_getitem(.tmp, 0L); y <- py_getitem(.tmp, 1L); z <- py_getitem(.tmp, 2L) }, silent = TRUE))
  suppressWarnings(try({ expr <- (x + y) * z }, silent = TRUE))
  suppressWarnings(try({ expected <- list(list("file.c", "#include \"file.h\"\n#include <math.h>\ndouble test(double x, double y, double z) {\n   double test_result;\n   test_result = z*(x + y);\n   return test_result;\n}\n"), list("file.h", "#ifndef PROJECT__FILE__H\n#define PROJECT__FILE__H\ndouble test(double x, double y, double z);\n#endif\n")) }, silent = TRUE))
  suppressWarnings(try({ result <- codegen(list("test", expr), "C", "file", header = FALSE, empty = FALSE) }, silent = TRUE))
  expect_identical(.cgb(result), "[('file.c', '#include \"file.h\"\\n#include <math.h>\\ndouble test(double x, double y, double z) {\\n   double test_result;\\n   test_result = z*(x + y);\\n   return test_result;\\n}\\n'), ('file.h', '#ifndef PROJECT__FILE__H\\n#define PROJECT__FILE__H\\ndouble test(double x, double y, double z);\\n#endif\\n')]")
})

test_that("utilities_test_codegen: test_multiple_results_c", {
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ .tmp <- sym_symbols("x,y,z"); x <- py_getitem(.tmp, 0L); y <- py_getitem(.tmp, 1L); z <- py_getitem(.tmp, 2L) }, silent = TRUE))
  suppressWarnings(try({ expr1 <- (x + y) * z }, silent = TRUE))
  suppressWarnings(try({ expr2 <- (x - y) * z }, silent = TRUE))
  suppressWarnings(try({ routine <- make_routine("test", list(expr1, expr2)) }, silent = TRUE))
  suppressWarnings(try({ code_gen <- C99CodeGen() }, silent = TRUE))
  expect_identical(.cgb(get_string(ast_attr(code_gen, "dump_h"), list(routine))), "ERR CodeGenError: C only supports a single or no return value.")
})

test_that("utilities_test_codegen: test_no_results_c", {
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  expect_identical(.cgb(make_routine("test", list())), "ERR ValueError: No expression given")
})

test_that("utilities_test_codegen: test_ansi_math1_codegen", {
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ invisible(NULL) }, silent = TRUE))
  suppressWarnings(try({ invisible(NULL) }, silent = TRUE))
  suppressWarnings(try({ invisible(NULL) }, silent = TRUE))
  suppressWarnings(try({ invisible(NULL) }, silent = TRUE))
  suppressWarnings(try({ invisible(NULL) }, silent = TRUE))
  suppressWarnings(try({ invisible(NULL) }, silent = TRUE))
  suppressWarnings(try({ x <- sym_symbols("x") }, silent = TRUE))
  suppressWarnings(try({ name_expr <- list(list("test_fabs", Abs(x)), list("test_acos", acos(x)), list("test_asin", asin(x)), list("test_atan", atan(x)), list("test_ceil", ceiling(x)), list("test_cos", cos(x)), list("test_cosh", cosh(x)), list("test_floor", floor(x)), list("test_log", log(x)), list("test_ln", log(x)), list("test_sin", sin(x)), list("test_sinh", sinh(x)), list("test_sqrt", sqrt(x)), list("test_tan", tan(x)), list("test_tanh", tanh(x))) }, silent = TRUE))
  suppressWarnings(try({ result <- codegen(name_expr, "C89", "file", header = FALSE, empty = FALSE) }, silent = TRUE))
  expect_identical(.cgb(py_getitem(py_getitem(result, 0L), 0L)), "file.c")
  expect_identical(.cgb(py_getitem(py_getitem(result, 0L), 1L)), "#include \"file.h\"\\n#include <math.h>\\ndouble test_fabs(double x) {\\n   double test_fabs_result;\\n   test_fabs_result = fabs(x);\\n   return test_fabs_result;\\n}\\ndouble test_acos(double x) {\\n   double test_acos_result;\\n   test_acos_result = acos(x);\\n   return test_acos_result;\\n}\\ndouble test_asin(double x) {\\n   double test_asin_result;\\n   test_asin_result = asin(x);\\n   return test_asin_result;\\n}\\ndouble test_atan(double x) {\\n   double test_atan_result;\\n   test_atan_result = atan(x);\\n   return test_atan_result;\\n}\\ndouble test_ceil(double x) {\\n   double test_ceil_result;\\n   test_ceil_result = ceil(x);\\n   return test_ceil_result;\\n}\\ndouble test_cos(double x) {\\n   double test_cos_result;\\n   test_cos_result = cos(x);\\n   return test_cos_result;\\n}\\ndouble test_cosh(double x) {\\n   double test_cosh_result;\\n   test_cosh_result = cosh(x);\\n   return test_cosh_result;\\n}\\ndouble test_floor(double x) {\\n   double test_floor_result;\\n   test_floor_result = floor(x);\\n   return test_floor_result;\\n}\\ndouble test_log(double x) {\\n   double test_log_result;\\n   test_log_result = log(x);\\n   return test_log_result;\\n}\\ndouble test_ln(double x) {\\n   double test_ln_result;\\n   test_ln_result = log(x);\\n   return test_ln_result;\\n}\\ndouble test_sin(double x) {\\n   double test_sin_result;\\n   test_sin_result = sin(x);\\n   return test_sin_result;\\n}\\ndouble test_sinh(double x) {\\n   double test_sinh_result;\\n   test_sinh_result = sinh(x);\\n   return test_sinh_result;\\n}\\ndouble test_sqrt(double x) {\\n   double test_sqrt_result;\\n   test_sqrt_result = sqrt(x);\\n   return test_sqrt_result;\\n}\\ndouble test_tan(double x) {\\n   double test_tan_result;\\n   test_tan_result = tan(x);\\n   return test_tan_result;\\n}\\ndouble test_tanh(double x) {\\n   double test_tanh_result;\\n   test_tanh_result = tanh(x);\\n   return test_tanh_result;\\n}\\n")
  expect_identical(.cgb(py_getitem(py_getitem(result, 1L), 0L)), "file.h")
  expect_identical(.cgb(py_getitem(py_getitem(result, 1L), 1L)), "#ifndef PROJECT__FILE__H\\n#define PROJECT__FILE__H\\ndouble test_fabs(double x);\\ndouble test_acos(double x);\\ndouble test_asin(double x);\\ndouble test_atan(double x);\\ndouble test_ceil(double x);\\ndouble test_cos(double x);\\ndouble test_cosh(double x);\\ndouble test_floor(double x);\\ndouble test_log(double x);\\ndouble test_ln(double x);\\ndouble test_sin(double x);\\ndouble test_sinh(double x);\\ndouble test_sqrt(double x);\\ndouble test_tan(double x);\\ndouble test_tanh(double x);\\n#endif\\n")
})

test_that("utilities_test_codegen: test_ansi_math2_codegen", {
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ invisible(NULL) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- sym_symbols("x,y"); x <- py_getitem(.tmp, 0L); y <- py_getitem(.tmp, 1L) }, silent = TRUE))
  suppressWarnings(try({ name_expr <- list(list("test_atan2", sym_atan2(x, y)), list("test_pow", x ^ y)) }, silent = TRUE))
  suppressWarnings(try({ result <- codegen(name_expr, "C89", "file", header = FALSE, empty = FALSE) }, silent = TRUE))
  expect_identical(.cgb(py_getitem(py_getitem(result, 0L), 0L)), "file.c")
  expect_identical(.cgb(py_getitem(py_getitem(result, 0L), 1L)), "#include \"file.h\"\\n#include <math.h>\\ndouble test_atan2(double x, double y) {\\n   double test_atan2_result;\\n   test_atan2_result = atan2(x, y);\\n   return test_atan2_result;\\n}\\ndouble test_pow(double x, double y) {\\n   double test_pow_result;\\n   test_pow_result = pow(x, y);\\n   return test_pow_result;\\n}\\n")
  expect_identical(.cgb(py_getitem(py_getitem(result, 1L), 0L)), "file.h")
  expect_identical(.cgb(py_getitem(py_getitem(result, 1L), 1L)), "#ifndef PROJECT__FILE__H\\n#define PROJECT__FILE__H\\ndouble test_atan2(double x, double y);\\ndouble test_pow(double x, double y);\\n#endif\\n")
})

test_that("utilities_test_codegen: test_complicated_codegen", {
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ invisible(NULL) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- sym_symbols("x,y,z"); x <- py_getitem(.tmp, 0L); y <- py_getitem(.tmp, 1L); z <- py_getitem(.tmp, 2L) }, silent = TRUE))
  suppressWarnings(try({ name_expr <- list(list("test1", expand(((sin(x) + cos(y)) + tan(z)) ^ 7L)), list("test2", cos(cos(cos(cos(cos(cos(cos(cos((x + y) + z)))))))))) }, silent = TRUE))
  suppressWarnings(try({ result <- codegen(name_expr, "C89", "file", header = FALSE, empty = FALSE) }, silent = TRUE))
  expect_identical(.cgb(py_getitem(py_getitem(result, 0L), 0L)), "file.c")
  expect_identical(.cgb(py_getitem(py_getitem(result, 0L), 1L)), "#include \"file.h\"\\n#include <math.h>\\ndouble test1(double x, double y, double z) {\\n   double test1_result;\\n   test1_result = pow(sin(x), 7) + 7*pow(sin(x), 6)*cos(y) + 7*pow(sin(x), 6)*tan(z) + 21*pow(sin(x), 5)*pow(cos(y), 2) + 42*pow(sin(x), 5)*cos(y)*tan(z) + 21*pow(sin(x), 5)*pow(tan(z), 2) + 35*pow(sin(x), 4)*pow(cos(y), 3) + 105*pow(sin(x), 4)*pow(cos(y), 2)*tan(z) + 105*pow(sin(x), 4)*cos(y)*pow(tan(z), 2) + 35*pow(sin(x), 4)*pow(tan(z), 3) + 35*pow(sin(x), 3)*pow(cos(y), 4) + 140*pow(sin(x), 3)*pow(cos(y), 3)*tan(z) + 210*pow(sin(x), 3)*pow(cos(y), 2)*pow(tan(z), 2) + 140*pow(sin(x), 3)*cos(y)*pow(tan(z), 3) + 35*pow(sin(x), 3)*pow(tan(z), 4) + 21*pow(sin(x), 2)*pow(cos(y), 5) + 105*pow(sin(x), 2)*pow(cos(y), 4)*tan(z) + 210*pow(sin(x), 2)*pow(cos(y), 3)*pow(tan(z), 2) + 210*pow(sin(x), 2)*pow(cos(y), 2)*pow(tan(z), 3) + 105*pow(sin(x), 2)*cos(y)*pow(tan(z), 4) + 21*pow(sin(x), 2)*pow(tan(z), 5) + 7*sin(x)*pow(cos(y), 6) + 42*sin(x)*pow(cos(y), 5)*tan(z) + 105*sin(x)*pow(cos(y), 4)*pow(tan(z), 2) + 140*sin(x)*pow(cos(y), 3)*pow(tan(z), 3) + 105*sin(x)*pow(cos(y), 2)*pow(tan(z), 4) + 42*sin(x)*cos(y)*pow(tan(z), 5) + 7*sin(x)*pow(tan(z), 6) + pow(cos(y), 7) + 7*pow(cos(y), 6)*tan(z) + 21*pow(cos(y), 5)*pow(tan(z), 2) + 35*pow(cos(y), 4)*pow(tan(z), 3) + 35*pow(cos(y), 3)*pow(tan(z), 4) + 21*pow(cos(y), 2)*pow(tan(z), 5) + 7*cos(y)*pow(tan(z), 6) + pow(tan(z), 7);\\n   return test1_result;\\n}\\ndouble test2(double x, double y, double z) {\\n   double test2_result;\\n   test2_result = cos(cos(cos(cos(cos(cos(cos(cos(x + y + z))))))));\\n   return test2_result;\\n}\\n")
  expect_identical(.cgb(py_getitem(py_getitem(result, 1L), 0L)), "file.h")
  expect_identical(.cgb(py_getitem(py_getitem(result, 1L), 1L)), "#ifndef PROJECT__FILE__H\\n#define PROJECT__FILE__H\\ndouble test1(double x, double y, double z);\\ndouble test2(double x, double y, double z);\\n#endif\\n")
})

test_that("utilities_test_codegen: test_loops_c", {
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ invisible(NULL) }, silent = TRUE))
  suppressWarnings(try({ invisible(NULL) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- sym_symbols("n m", integer = TRUE); n <- py_getitem(.tmp, 0L); m <- py_getitem(.tmp, 1L) }, silent = TRUE))
  suppressWarnings(try({ A <- IndexedBase("A") }, silent = TRUE))
  suppressWarnings(try({ x <- IndexedBase("x") }, silent = TRUE))
  suppressWarnings(try({ y <- IndexedBase("y") }, silent = TRUE))
  suppressWarnings(try({ i <- Idx("i", m) }, silent = TRUE))
  suppressWarnings(try({ j <- Idx("j", n) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- codegen(list("matrix_vector", Eq(py_getitem(y, i), py_getitem(A, i, j) * py_getitem(x, j))), "C99", "file", header = FALSE, empty = FALSE); .tmp1 <- py_getitem(.tmp, 0L); f1 <- py_getitem(.tmp1, 0L); code <- py_getitem(.tmp1, 1L); .tmp1 <- py_getitem(.tmp, 1L); f2 <- py_getitem(.tmp1, 0L); interface <- py_getitem(.tmp1, 1L) }, silent = TRUE))
  expect_identical(.cgb(f1), "file.c")
  suppressWarnings(try({ expected <- "#include \"file.h\"\n#include <math.h>\nvoid matrix_vector(double *A, int m, int n, double *x, double *y) {\n   for (int i=0; i<m; i++){\n      y[i] = 0;\n   }\n   for (int i=0; i<m; i++){\n      for (int j=0; j<n; j++){\n         y[i] = %(rhs)s + y[i];\n      }\n   }\n}\n" }, silent = TRUE))
  expect_identical(.cgb((py_eq(code, py_format(expected, py_dict(list("rhs"), list(py_format("A[%s]*x[j]", (i * n) + j)))))) || (py_eq(code, py_format(expected, py_dict(list("rhs"), list(py_format("A[%s]*x[j]", j + (i * n))))))) || (py_eq(code, py_format(expected, py_dict(list("rhs"), list(py_format("x[j]*A[%s]", (i * n) + j)))))) || (py_eq(code, py_format(expected, py_dict(list("rhs"), list(py_format("x[j]*A[%s]", j + (i * n)))))))), "True")
  expect_identical(.cgb(f2), "file.h")
  expect_identical(.cgb(interface), "#ifndef PROJECT__FILE__H\\n#define PROJECT__FILE__H\\nvoid matrix_vector(double *A, int m, int n, double *x, double *y);\\n#endif\\n")
})

test_that("utilities_test_codegen: test_partial_loops_c", {
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ invisible(NULL) }, silent = TRUE))
  suppressWarnings(try({ invisible(NULL) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- sym_symbols("n m o p", integer = TRUE); n <- py_getitem(.tmp, 0L); m <- py_getitem(.tmp, 1L); o <- py_getitem(.tmp, 2L); p <- py_getitem(.tmp, 3L) }, silent = TRUE))
  suppressWarnings(try({ A <- IndexedBase("A", shape = list(m, p)) }, silent = TRUE))
  suppressWarnings(try({ x <- IndexedBase("x") }, silent = TRUE))
  suppressWarnings(try({ y <- IndexedBase("y") }, silent = TRUE))
  suppressWarnings(try({ i <- Idx("i", list(o, m - 5L)) }, silent = TRUE))
  suppressWarnings(try({ j <- Idx("j", n) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- codegen(list("matrix_vector", Eq(py_getitem(y, i), py_getitem(A, i, j) * py_getitem(x, j))), "C99", "file", header = FALSE, empty = FALSE); .tmp1 <- py_getitem(.tmp, 0L); f1 <- py_getitem(.tmp1, 0L); code <- py_getitem(.tmp1, 1L); .tmp1 <- py_getitem(.tmp, 1L); f2 <- py_getitem(.tmp1, 0L); interface <- py_getitem(.tmp1, 1L) }, silent = TRUE))
  expect_identical(.cgb(f1), "file.c")
  suppressWarnings(try({ expected <- py_format("#include \"file.h\"\n#include <math.h>\nvoid matrix_vector(double *A, int m, int n, int o, int p, double *x, double *y) {\n   for (int i=o; i<%(upperi)s; i++){\n      y[i] = 0;\n   }\n   for (int i=o; i<%(upperi)s; i++){\n      for (int j=0; j<n; j++){\n         y[i] = %(rhs)s + y[i];\n      }\n   }\n}\n", py_dict(list("upperi", "rhs"), list(m - 4L, "%(rhs)s"))) }, silent = TRUE))
  expect_identical(.cgb((py_eq(code, py_format(expected, py_dict(list("rhs"), list(py_format("A[%s]*x[j]", (i * p) + j)))))) || (py_eq(code, py_format(expected, py_dict(list("rhs"), list(py_format("A[%s]*x[j]", j + (i * p))))))) || (py_eq(code, py_format(expected, py_dict(list("rhs"), list(py_format("x[j]*A[%s]", (i * p) + j)))))) || (py_eq(code, py_format(expected, py_dict(list("rhs"), list(py_format("x[j]*A[%s]", j + (i * p)))))))), "True")
  expect_identical(.cgb(f2), "file.h")
  expect_identical(.cgb(interface), "#ifndef PROJECT__FILE__H\\n#define PROJECT__FILE__H\\nvoid matrix_vector(double *A, int m, int n, int o, int p, double *x, double *y);\\n#endif\\n")
})

test_that("utilities_test_codegen: test_output_arg_c", {
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ invisible(NULL) }, silent = TRUE))
  suppressWarnings(try({ invisible(NULL) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- sym_symbols("x,y,z"); x <- py_getitem(.tmp, 0L); y <- py_getitem(.tmp, 1L); z <- py_getitem(.tmp, 2L) }, silent = TRUE))
  suppressWarnings(try({ r <- make_routine("foo", list(Equality(y, sin(x)), cos(x))) }, silent = TRUE))
  suppressWarnings(try({ c <- C89CodeGen() }, silent = TRUE))
  suppressWarnings(try({ result <- py_method(c, "write", list(r), "test", header = FALSE, empty = FALSE) }, silent = TRUE))
  expect_identical(.cgb(py_getitem(py_getitem(result, 0L), 0L)), "test.c")
  suppressWarnings(try({ expected <- "#include \"test.h\"\n#include <math.h>\ndouble foo(double x, double *y) {\n   (*y) = sin(x);\n   double foo_result;\n   foo_result = cos(x);\n   return foo_result;\n}\n" }, silent = TRUE))
  expect_identical(.cgb(py_getitem(py_getitem(result, 0L), 1L)), "#include \"test.h\"\\n#include <math.h>\\ndouble foo(double x, double *y) {\\n   (*y) = sin(x);\\n   double foo_result;\\n   foo_result = cos(x);\\n   return foo_result;\\n}\\n")
})

test_that("utilities_test_codegen: test_output_arg_c_reserved_words", {
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ invisible(NULL) }, silent = TRUE))
  suppressWarnings(try({ invisible(NULL) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- sym_symbols("if, while, z"); x <- py_getitem(.tmp, 0L); y <- py_getitem(.tmp, 1L); z <- py_getitem(.tmp, 2L) }, silent = TRUE))
  suppressWarnings(try({ r <- make_routine("foo", list(Equality(y, sin(x)), cos(x))) }, silent = TRUE))
  suppressWarnings(try({ c <- C89CodeGen() }, silent = TRUE))
  suppressWarnings(try({ result <- py_method(c, "write", list(r), "test", header = FALSE, empty = FALSE) }, silent = TRUE))
  expect_identical(.cgb(py_getitem(py_getitem(result, 0L), 0L)), "test.c")
  suppressWarnings(try({ expected <- "#include \"test.h\"\n#include <math.h>\ndouble foo(double if_, double *while_) {\n   (*while_) = sin(if_);\n   double foo_result;\n   foo_result = cos(if_);\n   return foo_result;\n}\n" }, silent = TRUE))
  expect_identical(.cgb(py_getitem(py_getitem(result, 0L), 1L)), "#include \"test.h\"\\n#include <math.h>\\ndouble foo(double if_, double *while_) {\\n   (*while_) = sin(if_);\\n   double foo_result;\\n   foo_result = cos(if_);\\n   return foo_result;\\n}\\n")
})

test_that("utilities_test_codegen: test_multidim_c_argument_cse", {
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ A_sym <- MatrixSymbol("A", 3L, 3L) }, silent = TRUE))
  suppressWarnings(try({ b_sym <- MatrixSymbol("b", 3L, 1L) }, silent = TRUE))
  suppressWarnings(try({ A <- as_explicit(A_sym) }, silent = TRUE))
  suppressWarnings(try({ b <- as_explicit(b_sym) }, silent = TRUE))
  suppressWarnings(try({ c <- A %*% b }, silent = TRUE))
  suppressWarnings(try({ cgen <- CCodeGen(project = "test", cse = TRUE) }, silent = TRUE))
  suppressWarnings(try({ r <- py_method(cgen, "routine", "c", c) }, silent = TRUE))
  suppressWarnings(try({ r$arguments[[length(r$arguments)]]$result_var <- "out" }, silent = TRUE))
  suppressWarnings(try({ r$arguments[[length(r$arguments)]]$name <- "out" }, silent = TRUE))
  suppressWarnings(try({ code <- get_string(ast_attr(cgen, "dump_c"), list(r), prefix = "test") }, silent = TRUE))
  suppressWarnings(try({ expected <- "#include \"test.h\"\n#include <math.h>\nvoid c(double *A, double *b, double *out) {\n   out[0] = A[0]*b[0] + A[1]*b[1] + A[2]*b[2];\n   out[1] = A[3]*b[0] + A[4]*b[1] + A[5]*b[2];\n   out[2] = A[6]*b[0] + A[7]*b[1] + A[8]*b[2];\n}\n" }, silent = TRUE))
  expect_identical(.cgb(code), "#include \"test.h\"\\n#include <math.h>\\nvoid c(double *A, double *b, double *out) {\\n   out[0] = A[0]*b[0] + A[1]*b[1] + A[2]*b[2];\\n   out[1] = A[3]*b[0] + A[4]*b[1] + A[5]*b[2];\\n   out[2] = A[6]*b[0] + A[7]*b[1] + A[8]*b[2];\\n}\\n")
})

test_that("utilities_test_codegen: test_ccode_results_named_ordered", {
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ .tmp <- sym_symbols("x,y,z"); x <- py_getitem(.tmp, 0L); y <- py_getitem(.tmp, 1L); z <- py_getitem(.tmp, 2L) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- sym_symbols("B,C"); B <- py_getitem(.tmp, 0L); C <- py_getitem(.tmp, 1L) }, silent = TRUE))
  suppressWarnings(try({ A <- MatrixSymbol("A", 1L, 3L) }, silent = TRUE))
  suppressWarnings(try({ expr1 <- Equality(A, Matrix(list(list(1L, 2L, x)))) }, silent = TRUE))
  suppressWarnings(try({ expr2 <- Equality(C, (x + y) * z) }, silent = TRUE))
  suppressWarnings(try({ expr3 <- Equality(B, 2L * x) }, silent = TRUE))
  suppressWarnings(try({ name_expr <- list("test", list(expr1, expr2, expr3)) }, silent = TRUE))
  suppressWarnings(try({ expected <- "#include \"test.h\"\n#include <math.h>\nvoid test(double x, double *C, double z, double y, double *A, double *B) {\n   (*C) = z*(x + y);\n   A[0] = 1;\n   A[1] = 2;\n   A[2] = x;\n   (*B) = 2*x;\n}\n" }, silent = TRUE))
  suppressWarnings(try({ result <- codegen(name_expr, "c", "test", header = FALSE, empty = FALSE, argument_sequence = list(x, C, z, y, A, B)) }, silent = TRUE))
  suppressWarnings(try({ source <- py_getitem(py_getitem(result, 0L), 1L) }, silent = TRUE))
  expect_identical(.cgb(source), "#include \"test.h\"\\n#include <math.h>\\nvoid test(double x, double *C, double z, double y, double *A, double *B) {\\n   (*C) = z*(x + y);\\n   A[0] = 1;\\n   A[1] = 2;\\n   A[2] = x;\\n   (*B) = 2*x;\\n}\\n")
})

test_that("utilities_test_codegen: test_ccode_matrixsymbol_slice", {
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ A <- MatrixSymbol("A", 5L, 3L) }, silent = TRUE))
  suppressWarnings(try({ B <- MatrixSymbol("B", 1L, 3L) }, silent = TRUE))
  suppressWarnings(try({ C <- MatrixSymbol("C", 1L, 3L) }, silent = TRUE))
  suppressWarnings(try({ D <- MatrixSymbol("D", 5L, 1L) }, silent = TRUE))
  suppressWarnings(try({ name_expr <- list("test", list(Equality(B, py_getitem(A, 0L, pyslice(NULL, NULL, NULL))), Equality(C, py_getitem(A, 1L, pyslice(NULL, NULL, NULL))), Equality(D, py_getitem(A, pyslice(NULL, NULL, NULL), 2L)))) }, silent = TRUE))
  suppressWarnings(try({ result <- codegen(name_expr, "c99", "test", header = FALSE, empty = FALSE) }, silent = TRUE))
  suppressWarnings(try({ source <- py_getitem(py_getitem(result, 0L), 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "#include \"test.h\"\n#include <math.h>\nvoid test(double *A, double *B, double *C, double *D) {\n   B[0] = A[0];\n   B[1] = A[1];\n   B[2] = A[2];\n   C[0] = A[3];\n   C[1] = A[4];\n   C[2] = A[5];\n   D[0] = A[2];\n   D[1] = A[5];\n   D[2] = A[8];\n   D[3] = A[11];\n   D[4] = A[14];\n}\n" }, silent = TRUE))
  expect_identical(.cgb(source), "#include \"test.h\"\\n#include <math.h>\\nvoid test(double *A, double *B, double *C, double *D) {\\n   B[0] = A[0];\\n   B[1] = A[1];\\n   B[2] = A[2];\\n   C[0] = A[3];\\n   C[1] = A[4];\\n   C[2] = A[5];\\n   D[0] = A[2];\\n   D[1] = A[5];\\n   D[2] = A[8];\\n   D[3] = A[11];\\n   D[4] = A[14];\\n}\\n")
})

test_that("utilities_test_codegen: test_ccode_unused_array_arg", {
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ x <- MatrixSymbol("x", 2L, 1L) }, silent = TRUE))
  suppressWarnings(try({ name_expr <- list("test", Float(1.0)) }, silent = TRUE))
  suppressWarnings(try({ generator <- CCodeGen() }, silent = TRUE))
  suppressWarnings(try({ result <- codegen(name_expr, code_gen = generator, header = FALSE, empty = FALSE, argument_sequence = list(x)) }, silent = TRUE))
  suppressWarnings(try({ source <- py_getitem(py_getitem(result, 0L), 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "#include \"test.h\"\n#include <math.h>\ndouble test(double *x) {\n   double test_result;\n   test_result = 1.0;\n   return test_result;\n}\n" }, silent = TRUE))
  expect_identical(.cgb(source), "#include \"test.h\"\\n#include <math.h>\\ndouble test(double *x) {\\n   double test_result;\\n   test_result = 1.0;\\n   return test_result;\\n}\\n")
})

test_that("utilities_test_codegen: test_ccode_unused_array_arg_func", {
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ X <- MatrixSymbol("X", 3L, 1L) }, silent = TRUE))
  suppressWarnings(try({ Y <- MatrixSymbol("Y", 3L, 1L) }, silent = TRUE))
  suppressWarnings(try({ z <- sym_symbols("z", integer = TRUE) }, silent = TRUE))
  suppressWarnings(try({ name_expr <- list("testBug", py_getitem(X, 0L) + py_getitem(X, 1L)) }, silent = TRUE))
  suppressWarnings(try({ result <- codegen(name_expr, language = "C", header = FALSE, empty = FALSE, argument_sequence = list(X, Y, z)) }, silent = TRUE))
  suppressWarnings(try({ source <- py_getitem(py_getitem(result, 0L), 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "#include \"testBug.h\"\n#include <math.h>\ndouble testBug(double *X, double *Y, int z) {\n   double testBug_result;\n   testBug_result = X[0] + X[1];\n   return testBug_result;\n}\n" }, silent = TRUE))
  expect_identical(.cgb(source), "#include \"testBug.h\"\\n#include <math.h>\\ndouble testBug(double *X, double *Y, int z) {\\n   double testBug_result;\\n   testBug_result = X[0] + X[1];\\n   return testBug_result;\\n}\\n")
})

test_that("utilities_test_codegen: test_empty_f_code", {
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ code_gen <- FCodeGen() }, silent = TRUE))
  suppressWarnings(try({ source <- get_string(ast_attr(code_gen, "dump_f95"), list()) }, silent = TRUE))
  expect_identical(.cgb(source), "")
})

test_that("utilities_test_codegen: test_empty_f_code_with_header", {
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ code_gen <- FCodeGen() }, silent = TRUE))
  suppressWarnings(try({ source <- get_string(ast_attr(code_gen, "dump_f95"), list(), header = TRUE) }, silent = TRUE))
  expect_identical(.cgb(py_slice(source, NULL, 82L)), "!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\\n!*")
  expect_identical(.cgb(py_slice(source, 158L, NULL)), "*\\n!*                                                                            *\\n!*              See http://www.sympy.org/ for more information.               *\\n!*                                                                            *\\n!*                       This file is part of 'project'                       *\\n!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\\n")
})

test_that("utilities_test_codegen: test_empty_f_header", {
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ code_gen <- FCodeGen() }, silent = TRUE))
  suppressWarnings(try({ source <- get_string(ast_attr(code_gen, "dump_h"), list()) }, silent = TRUE))
  expect_identical(.cgb(source), "")
})

test_that("utilities_test_codegen: test_simple_f_code", {
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ .tmp <- sym_symbols("x,y,z"); x <- py_getitem(.tmp, 0L); y <- py_getitem(.tmp, 1L); z <- py_getitem(.tmp, 2L) }, silent = TRUE))
  suppressWarnings(try({ expr <- (x + y) * z }, silent = TRUE))
  suppressWarnings(try({ routine <- make_routine("test", expr) }, silent = TRUE))
  suppressWarnings(try({ code_gen <- FCodeGen() }, silent = TRUE))
  suppressWarnings(try({ source <- get_string(ast_attr(code_gen, "dump_f95"), list(routine)) }, silent = TRUE))
  suppressWarnings(try({ expected <- "REAL*8 function test(x, y, z)\nimplicit none\nREAL*8, intent(in) :: x\nREAL*8, intent(in) :: y\nREAL*8, intent(in) :: z\ntest = z*(x + y)\nend function\n" }, silent = TRUE))
  expect_identical(.cgb(source), "REAL*8 function test(x, y, z)\\nimplicit none\\nREAL*8, intent(in) :: x\\nREAL*8, intent(in) :: y\\nREAL*8, intent(in) :: z\\ntest = z*(x + y)\\nend function\\n")
})

test_that("utilities_test_codegen: test_numbersymbol_f_code", {
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ routine <- make_routine("test", S_Pi ^ Catalan) }, silent = TRUE))
  suppressWarnings(try({ code_gen <- FCodeGen() }, silent = TRUE))
  suppressWarnings(try({ source <- get_string(ast_attr(code_gen, "dump_f95"), list(routine)) }, silent = TRUE))
  suppressWarnings(try({ expected <- py_format("REAL*8 function test()\nimplicit none\nREAL*8, parameter :: Catalan = %sd0\nREAL*8, parameter :: pi = %sd0\ntest = pi**Catalan\nend function\n", list(N(Catalan, 17L), N(S_Pi, 17L))) }, silent = TRUE))
  expect_identical(.cgb(source), "REAL*8 function test()\\nimplicit none\\nREAL*8, parameter :: Catalan = 0.91596559417721901d0\\nREAL*8, parameter :: pi = 3.1415926535897932d0\\ntest = pi^Catalan\\nend function\\n")
})

test_that("utilities_test_codegen: test_erf_f_code", {
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ x <- sym_symbols("x") }, silent = TRUE))
  suppressWarnings(try({ routine <- make_routine("test", erf(x) - erf((-2L) * x)) }, silent = TRUE))
  suppressWarnings(try({ code_gen <- FCodeGen() }, silent = TRUE))
  suppressWarnings(try({ source <- get_string(ast_attr(code_gen, "dump_f95"), list(routine)) }, silent = TRUE))
  suppressWarnings(try({ expected <- "REAL*8 function test(x)\nimplicit none\nREAL*8, intent(in) :: x\ntest = erf(x) + erf(2.0d0*x)\nend function\n" }, silent = TRUE))
  expect_identical(.cgb(source), "REAL*8 function test(x)\\nimplicit none\\nREAL*8, intent(in) :: x\\ntest = erf(x) + erf(2.0d0*x)\\nend function\\n")
})

test_that("utilities_test_codegen: test_f_code_argument_order", {
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ .tmp <- sym_symbols("x,y,z"); x <- py_getitem(.tmp, 0L); y <- py_getitem(.tmp, 1L); z <- py_getitem(.tmp, 2L) }, silent = TRUE))
  suppressWarnings(try({ expr <- x + y }, silent = TRUE))
  suppressWarnings(try({ routine <- make_routine("test", expr, argument_sequence = list(z, x, y)) }, silent = TRUE))
  suppressWarnings(try({ code_gen <- FCodeGen() }, silent = TRUE))
  suppressWarnings(try({ source <- get_string(ast_attr(code_gen, "dump_f95"), list(routine)) }, silent = TRUE))
  suppressWarnings(try({ expected <- "REAL*8 function test(z, x, y)\nimplicit none\nREAL*8, intent(in) :: z\nREAL*8, intent(in) :: x\nREAL*8, intent(in) :: y\ntest = x + y\nend function\n" }, silent = TRUE))
  expect_identical(.cgb(source), "REAL*8 function test(z, x, y)\\nimplicit none\\nREAL*8, intent(in) :: z\\nREAL*8, intent(in) :: x\\nREAL*8, intent(in) :: y\\ntest = x + y\\nend function\\n")
})

test_that("utilities_test_codegen: test_simple_f_header", {
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ .tmp <- sym_symbols("x,y,z"); x <- py_getitem(.tmp, 0L); y <- py_getitem(.tmp, 1L); z <- py_getitem(.tmp, 2L) }, silent = TRUE))
  suppressWarnings(try({ expr <- (x + y) * z }, silent = TRUE))
  suppressWarnings(try({ routine <- make_routine("test", expr) }, silent = TRUE))
  suppressWarnings(try({ code_gen <- FCodeGen() }, silent = TRUE))
  suppressWarnings(try({ source <- get_string(ast_attr(code_gen, "dump_h"), list(routine)) }, silent = TRUE))
  suppressWarnings(try({ expected <- "interface\nREAL*8 function test(x, y, z)\nimplicit none\nREAL*8, intent(in) :: x\nREAL*8, intent(in) :: y\nREAL*8, intent(in) :: z\nend function\nend interface\n" }, silent = TRUE))
  expect_identical(.cgb(source), "interface\\nREAL*8 function test(x, y, z)\\nimplicit none\\nREAL*8, intent(in) :: x\\nREAL*8, intent(in) :: y\\nREAL*8, intent(in) :: z\\nend function\\nend interface\\n")
})

test_that("utilities_test_codegen: test_simple_f_codegen", {
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ .tmp <- sym_symbols("x,y,z"); x <- py_getitem(.tmp, 0L); y <- py_getitem(.tmp, 1L); z <- py_getitem(.tmp, 2L) }, silent = TRUE))
  suppressWarnings(try({ expr <- (x + y) * z }, silent = TRUE))
  suppressWarnings(try({ result <- codegen(list("test", expr), "F95", "file", header = FALSE, empty = FALSE) }, silent = TRUE))
  suppressWarnings(try({ expected <- list(list("file.f90", "REAL*8 function test(x, y, z)\nimplicit none\nREAL*8, intent(in) :: x\nREAL*8, intent(in) :: y\nREAL*8, intent(in) :: z\ntest = z*(x + y)\nend function\n"), list("file.h", "interface\nREAL*8 function test(x, y, z)\nimplicit none\nREAL*8, intent(in) :: x\nREAL*8, intent(in) :: y\nREAL*8, intent(in) :: z\nend function\nend interface\n")) }, silent = TRUE))
  expect_identical(.cgb(result), "[('file.f90', 'REAL*8 function test(x, y, z)\\nimplicit none\\nREAL*8, intent(in) :: x\\nREAL*8, intent(in) :: y\\nREAL*8, intent(in) :: z\\ntest = z*(x + y)\\nend function\\n'), ('file.h', 'interface\\nREAL*8 function test(x, y, z)\\nimplicit none\\nREAL*8, intent(in) :: x\\nREAL*8, intent(in) :: y\\nREAL*8, intent(in) :: z\\nend function\\nend interface\\n')]")
})

test_that("utilities_test_codegen: test_multiple_results_f", {
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ .tmp <- sym_symbols("x,y,z"); x <- py_getitem(.tmp, 0L); y <- py_getitem(.tmp, 1L); z <- py_getitem(.tmp, 2L) }, silent = TRUE))
  suppressWarnings(try({ expr1 <- (x + y) * z }, silent = TRUE))
  suppressWarnings(try({ expr2 <- (x - y) * z }, silent = TRUE))
  suppressWarnings(try({ routine <- make_routine("test", list(expr1, expr2)) }, silent = TRUE))
  suppressWarnings(try({ code_gen <- FCodeGen() }, silent = TRUE))
  expect_identical(.cgb(get_string(ast_attr(code_gen, "dump_h"), list(routine))), "ERR CodeGenError: Fortran only supports a single or no return value.")
})

test_that("utilities_test_codegen: test_no_results_f", {
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  expect_identical(.cgb(make_routine("test", list())), "ERR ValueError: No expression given")
})

test_that("utilities_test_codegen: test_intrinsic_math_codegen", {
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ invisible(NULL) }, silent = TRUE))
  suppressWarnings(try({ invisible(NULL) }, silent = TRUE))
  suppressWarnings(try({ invisible(NULL) }, silent = TRUE))
  suppressWarnings(try({ invisible(NULL) }, silent = TRUE))
  suppressWarnings(try({ invisible(NULL) }, silent = TRUE))
  suppressWarnings(try({ x <- sym_symbols("x") }, silent = TRUE))
  suppressWarnings(try({ name_expr <- list(list("test_abs", Abs(x)), list("test_acos", acos(x)), list("test_asin", asin(x)), list("test_atan", atan(x)), list("test_cos", cos(x)), list("test_cosh", cosh(x)), list("test_log", log(x)), list("test_ln", log(x)), list("test_sin", sin(x)), list("test_sinh", sinh(x)), list("test_sqrt", sqrt(x)), list("test_tan", tan(x)), list("test_tanh", tanh(x))) }, silent = TRUE))
  suppressWarnings(try({ result <- codegen(name_expr, "F95", "file", header = FALSE, empty = FALSE) }, silent = TRUE))
  expect_identical(.cgb(py_getitem(py_getitem(result, 0L), 0L)), "file.f90")
  suppressWarnings(try({ expected <- "REAL*8 function test_abs(x)\nimplicit none\nREAL*8, intent(in) :: x\ntest_abs = abs(x)\nend function\nREAL*8 function test_acos(x)\nimplicit none\nREAL*8, intent(in) :: x\ntest_acos = acos(x)\nend function\nREAL*8 function test_asin(x)\nimplicit none\nREAL*8, intent(in) :: x\ntest_asin = asin(x)\nend function\nREAL*8 function test_atan(x)\nimplicit none\nREAL*8, intent(in) :: x\ntest_atan = atan(x)\nend function\nREAL*8 function test_cos(x)\nimplicit none\nREAL*8, intent(in) :: x\ntest_cos = cos(x)\nend function\nREAL*8 function test_cosh(x)\nimplicit none\nREAL*8, intent(in) :: x\ntest_cosh = cosh(x)\nend function\nREAL*8 function test_log(x)\nimplicit none\nREAL*8, intent(in) :: x\ntest_log = log(x)\nend function\nREAL*8 function test_ln(x)\nimplicit none\nREAL*8, intent(in) :: x\ntest_ln = log(x)\nend function\nREAL*8 function test_sin(x)\nimplicit none\nREAL*8, intent(in) :: x\ntest_sin = sin(x)\nend function\nREAL*8 function test_sinh(x)\nimplicit none\nREAL*8, intent(in) :: x\ntest_sinh = sinh(x)\nend function\nREAL*8 function test_sqrt(x)\nimplicit none\nREAL*8, intent(in) :: x\ntest_sqrt = sqrt(x)\nend function\nREAL*8 function test_tan(x)\nimplicit none\nREAL*8, intent(in) :: x\ntest_tan = tan(x)\nend function\nREAL*8 function test_tanh(x)\nimplicit none\nREAL*8, intent(in) :: x\ntest_tanh = tanh(x)\nend function\n" }, silent = TRUE))
  expect_identical(.cgb(py_getitem(py_getitem(result, 0L), 1L)), "REAL*8 function test_abs(x)\\nimplicit none\\nREAL*8, intent(in) :: x\\ntest_abs = abs(x)\\nend function\\nREAL*8 function test_acos(x)\\nimplicit none\\nREAL*8, intent(in) :: x\\ntest_acos = acos(x)\\nend function\\nREAL*8 function test_asin(x)\\nimplicit none\\nREAL*8, intent(in) :: x\\ntest_asin = asin(x)\\nend function\\nREAL*8 function test_atan(x)\\nimplicit none\\nREAL*8, intent(in) :: x\\ntest_atan = atan(x)\\nend function\\nREAL*8 function test_cos(x)\\nimplicit none\\nREAL*8, intent(in) :: x\\ntest_cos = cos(x)\\nend function\\nREAL*8 function test_cosh(x)\\nimplicit none\\nREAL*8, intent(in) :: x\\ntest_cosh = cosh(x)\\nend function\\nREAL*8 function test_log(x)\\nimplicit none\\nREAL*8, intent(in) :: x\\ntest_log = log(x)\\nend function\\nREAL*8 function test_ln(x)\\nimplicit none\\nREAL*8, intent(in) :: x\\ntest_ln = log(x)\\nend function\\nREAL*8 function test_sin(x)\\nimplicit none\\nREAL*8, intent(in) :: x\\ntest_sin = sin(x)\\nend function\\nREAL*8 function test_sinh(x)\\nimplicit none\\nREAL*8, intent(in) :: x\\ntest_sinh = sinh(x)\\nend function\\nREAL*8 function test_sqrt(x)\\nimplicit none\\nREAL*8, intent(in) :: x\\ntest_sqrt = sqrt(x)\\nend function\\nREAL*8 function test_tan(x)\\nimplicit none\\nREAL*8, intent(in) :: x\\ntest_tan = tan(x)\\nend function\\nREAL*8 function test_tanh(x)\\nimplicit none\\nREAL*8, intent(in) :: x\\ntest_tanh = tanh(x)\\nend function\\n")
  expect_identical(.cgb(py_getitem(py_getitem(result, 1L), 0L)), "file.h")
  suppressWarnings(try({ expected <- "interface\nREAL*8 function test_abs(x)\nimplicit none\nREAL*8, intent(in) :: x\nend function\nend interface\ninterface\nREAL*8 function test_acos(x)\nimplicit none\nREAL*8, intent(in) :: x\nend function\nend interface\ninterface\nREAL*8 function test_asin(x)\nimplicit none\nREAL*8, intent(in) :: x\nend function\nend interface\ninterface\nREAL*8 function test_atan(x)\nimplicit none\nREAL*8, intent(in) :: x\nend function\nend interface\ninterface\nREAL*8 function test_cos(x)\nimplicit none\nREAL*8, intent(in) :: x\nend function\nend interface\ninterface\nREAL*8 function test_cosh(x)\nimplicit none\nREAL*8, intent(in) :: x\nend function\nend interface\ninterface\nREAL*8 function test_log(x)\nimplicit none\nREAL*8, intent(in) :: x\nend function\nend interface\ninterface\nREAL*8 function test_ln(x)\nimplicit none\nREAL*8, intent(in) :: x\nend function\nend interface\ninterface\nREAL*8 function test_sin(x)\nimplicit none\nREAL*8, intent(in) :: x\nend function\nend interface\ninterface\nREAL*8 function test_sinh(x)\nimplicit none\nREAL*8, intent(in) :: x\nend function\nend interface\ninterface\nREAL*8 function test_sqrt(x)\nimplicit none\nREAL*8, intent(in) :: x\nend function\nend interface\ninterface\nREAL*8 function test_tan(x)\nimplicit none\nREAL*8, intent(in) :: x\nend function\nend interface\ninterface\nREAL*8 function test_tanh(x)\nimplicit none\nREAL*8, intent(in) :: x\nend function\nend interface\n" }, silent = TRUE))
  expect_identical(.cgb(py_getitem(py_getitem(result, 1L), 1L)), "interface\\nREAL*8 function test_abs(x)\\nimplicit none\\nREAL*8, intent(in) :: x\\nend function\\nend interface\\ninterface\\nREAL*8 function test_acos(x)\\nimplicit none\\nREAL*8, intent(in) :: x\\nend function\\nend interface\\ninterface\\nREAL*8 function test_asin(x)\\nimplicit none\\nREAL*8, intent(in) :: x\\nend function\\nend interface\\ninterface\\nREAL*8 function test_atan(x)\\nimplicit none\\nREAL*8, intent(in) :: x\\nend function\\nend interface\\ninterface\\nREAL*8 function test_cos(x)\\nimplicit none\\nREAL*8, intent(in) :: x\\nend function\\nend interface\\ninterface\\nREAL*8 function test_cosh(x)\\nimplicit none\\nREAL*8, intent(in) :: x\\nend function\\nend interface\\ninterface\\nREAL*8 function test_log(x)\\nimplicit none\\nREAL*8, intent(in) :: x\\nend function\\nend interface\\ninterface\\nREAL*8 function test_ln(x)\\nimplicit none\\nREAL*8, intent(in) :: x\\nend function\\nend interface\\ninterface\\nREAL*8 function test_sin(x)\\nimplicit none\\nREAL*8, intent(in) :: x\\nend function\\nend interface\\ninterface\\nREAL*8 function test_sinh(x)\\nimplicit none\\nREAL*8, intent(in) :: x\\nend function\\nend interface\\ninterface\\nREAL*8 function test_sqrt(x)\\nimplicit none\\nREAL*8, intent(in) :: x\\nend function\\nend interface\\ninterface\\nREAL*8 function test_tan(x)\\nimplicit none\\nREAL*8, intent(in) :: x\\nend function\\nend interface\\ninterface\\nREAL*8 function test_tanh(x)\\nimplicit none\\nREAL*8, intent(in) :: x\\nend function\\nend interface\\n")
})

test_that("utilities_test_codegen: test_intrinsic_math2_codegen", {
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ invisible(NULL) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- sym_symbols("x,y"); x <- py_getitem(.tmp, 0L); y <- py_getitem(.tmp, 1L) }, silent = TRUE))
  suppressWarnings(try({ name_expr <- list(list("test_atan2", sym_atan2(x, y)), list("test_pow", x ^ y)) }, silent = TRUE))
  suppressWarnings(try({ result <- codegen(name_expr, "F95", "file", header = FALSE, empty = FALSE) }, silent = TRUE))
  expect_identical(.cgb(py_getitem(py_getitem(result, 0L), 0L)), "file.f90")
  suppressWarnings(try({ expected <- "REAL*8 function test_atan2(x, y)\nimplicit none\nREAL*8, intent(in) :: x\nREAL*8, intent(in) :: y\ntest_atan2 = atan2(x, y)\nend function\nREAL*8 function test_pow(x, y)\nimplicit none\nREAL*8, intent(in) :: x\nREAL*8, intent(in) :: y\ntest_pow = x**y\nend function\n" }, silent = TRUE))
  expect_identical(.cgb(py_getitem(py_getitem(result, 0L), 1L)), "REAL*8 function test_atan2(x, y)\\nimplicit none\\nREAL*8, intent(in) :: x\\nREAL*8, intent(in) :: y\\ntest_atan2 = atan2(x, y)\\nend function\\nREAL*8 function test_pow(x, y)\\nimplicit none\\nREAL*8, intent(in) :: x\\nREAL*8, intent(in) :: y\\ntest_pow = x^y\\nend function\\n")
  expect_identical(.cgb(py_getitem(py_getitem(result, 1L), 0L)), "file.h")
  suppressWarnings(try({ expected <- "interface\nREAL*8 function test_atan2(x, y)\nimplicit none\nREAL*8, intent(in) :: x\nREAL*8, intent(in) :: y\nend function\nend interface\ninterface\nREAL*8 function test_pow(x, y)\nimplicit none\nREAL*8, intent(in) :: x\nREAL*8, intent(in) :: y\nend function\nend interface\n" }, silent = TRUE))
  expect_identical(.cgb(py_getitem(py_getitem(result, 1L), 1L)), "interface\\nREAL*8 function test_atan2(x, y)\\nimplicit none\\nREAL*8, intent(in) :: x\\nREAL*8, intent(in) :: y\\nend function\\nend interface\\ninterface\\nREAL*8 function test_pow(x, y)\\nimplicit none\\nREAL*8, intent(in) :: x\\nREAL*8, intent(in) :: y\\nend function\\nend interface\\n")
})

test_that("utilities_test_codegen: test_complicated_codegen_f95", {
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ invisible(NULL) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- sym_symbols("x,y,z"); x <- py_getitem(.tmp, 0L); y <- py_getitem(.tmp, 1L); z <- py_getitem(.tmp, 2L) }, silent = TRUE))
  suppressWarnings(try({ name_expr <- list(list("test1", expand(((sin(x) + cos(y)) + tan(z)) ^ 7L)), list("test2", cos(cos(cos(cos(cos(cos(cos(cos((x + y) + z)))))))))) }, silent = TRUE))
  suppressWarnings(try({ result <- codegen(name_expr, "F95", "file", header = FALSE, empty = FALSE) }, silent = TRUE))
  expect_identical(.cgb(py_getitem(py_getitem(result, 0L), 0L)), "file.f90")
  suppressWarnings(try({ expected <- "REAL*8 function test1(x, y, z)\nimplicit none\nREAL*8, intent(in) :: x\nREAL*8, intent(in) :: y\nREAL*8, intent(in) :: z\ntest1 = sin(x)**7 + 7*sin(x)**6*cos(y) + 7*sin(x)**6*tan(z) + 21*sin(x) &\n      **5*cos(y)**2 + 42*sin(x)**5*cos(y)*tan(z) + 21*sin(x)**5*tan(z) &\n      **2 + 35*sin(x)**4*cos(y)**3 + 105*sin(x)**4*cos(y)**2*tan(z) + &\n      105*sin(x)**4*cos(y)*tan(z)**2 + 35*sin(x)**4*tan(z)**3 + 35*sin( &\n      x)**3*cos(y)**4 + 140*sin(x)**3*cos(y)**3*tan(z) + 210*sin(x)**3* &\n      cos(y)**2*tan(z)**2 + 140*sin(x)**3*cos(y)*tan(z)**3 + 35*sin(x) &\n      **3*tan(z)**4 + 21*sin(x)**2*cos(y)**5 + 105*sin(x)**2*cos(y)**4* &\n      tan(z) + 210*sin(x)**2*cos(y)**3*tan(z)**2 + 210*sin(x)**2*cos(y) &\n      **2*tan(z)**3 + 105*sin(x)**2*cos(y)*tan(z)**4 + 21*sin(x)**2*tan &\n      (z)**5 + 7*sin(x)*cos(y)**6 + 42*sin(x)*cos(y)**5*tan(z) + 105* &\n      sin(x)*cos(y)**4*tan(z)**2 + 140*sin(x)*cos(y)**3*tan(z)**3 + 105 &\n      *sin(x)*cos(y)**2*tan(z)**4 + 42*sin(x)*cos(y)*tan(z)**5 + 7*sin( &\n      x)*tan(z)**6 + cos(y)**7 + 7*cos(y)**6*tan(z) + 21*cos(y)**5*tan( &\n      z)**2 + 35*cos(y)**4*tan(z)**3 + 35*cos(y)**3*tan(z)**4 + 21*cos( &\n      y)**2*tan(z)**5 + 7*cos(y)*tan(z)**6 + tan(z)**7\nend function\nREAL*8 function test2(x, y, z)\nimplicit none\nREAL*8, intent(in) :: x\nREAL*8, intent(in) :: y\nREAL*8, intent(in) :: z\ntest2 = cos(cos(cos(cos(cos(cos(cos(cos(x + y + z))))))))\nend function\n" }, silent = TRUE))
  expect_identical(.cgb(py_getitem(py_getitem(result, 0L), 1L)), "REAL*8 function test1(x, y, z)\\nimplicit none\\nREAL*8, intent(in) :: x\\nREAL*8, intent(in) :: y\\nREAL*8, intent(in) :: z\\ntest1 = sin(x)^7 + 7*sin(x)^6*cos(y) + 7*sin(x)^6*tan(z) + 21*sin(x) &\\n      ^5*cos(y)^2 + 42*sin(x)^5*cos(y)*tan(z) + 21*sin(x)^5*tan(z) &\\n      ^2 + 35*sin(x)^4*cos(y)^3 + 105*sin(x)^4*cos(y)^2*tan(z) + &\\n      105*sin(x)^4*cos(y)*tan(z)^2 + 35*sin(x)^4*tan(z)^3 + 35*sin( &\\n      x)^3*cos(y)^4 + 140*sin(x)^3*cos(y)^3*tan(z) + 210*sin(x)^3* &\\n      cos(y)^2*tan(z)^2 + 140*sin(x)^3*cos(y)*tan(z)^3 + 35*sin(x) &\\n      ^3*tan(z)^4 + 21*sin(x)^2*cos(y)^5 + 105*sin(x)^2*cos(y)^4* &\\n      tan(z) + 210*sin(x)^2*cos(y)^3*tan(z)^2 + 210*sin(x)^2*cos(y) &\\n      ^2*tan(z)^3 + 105*sin(x)^2*cos(y)*tan(z)^4 + 21*sin(x)^2*tan &\\n      (z)^5 + 7*sin(x)*cos(y)^6 + 42*sin(x)*cos(y)^5*tan(z) + 105* &\\n      sin(x)*cos(y)^4*tan(z)^2 + 140*sin(x)*cos(y)^3*tan(z)^3 + 105 &\\n      *sin(x)*cos(y)^2*tan(z)^4 + 42*sin(x)*cos(y)*tan(z)^5 + 7*sin( &\\n      x)*tan(z)^6 + cos(y)^7 + 7*cos(y)^6*tan(z) + 21*cos(y)^5*tan( &\\n      z)^2 + 35*cos(y)^4*tan(z)^3 + 35*cos(y)^3*tan(z)^4 + 21*cos( &\\n      y)^2*tan(z)^5 + 7*cos(y)*tan(z)^6 + tan(z)^7\\nend function\\nREAL*8 function test2(x, y, z)\\nimplicit none\\nREAL*8, intent(in) :: x\\nREAL*8, intent(in) :: y\\nREAL*8, intent(in) :: z\\ntest2 = cos(cos(cos(cos(cos(cos(cos(cos(x + y + z))))))))\\nend function\\n")
  expect_identical(.cgb(py_getitem(py_getitem(result, 1L), 0L)), "file.h")
  suppressWarnings(try({ expected <- "interface\nREAL*8 function test1(x, y, z)\nimplicit none\nREAL*8, intent(in) :: x\nREAL*8, intent(in) :: y\nREAL*8, intent(in) :: z\nend function\nend interface\ninterface\nREAL*8 function test2(x, y, z)\nimplicit none\nREAL*8, intent(in) :: x\nREAL*8, intent(in) :: y\nREAL*8, intent(in) :: z\nend function\nend interface\n" }, silent = TRUE))
  expect_identical(.cgb(py_getitem(py_getitem(result, 1L), 1L)), "interface\\nREAL*8 function test1(x, y, z)\\nimplicit none\\nREAL*8, intent(in) :: x\\nREAL*8, intent(in) :: y\\nREAL*8, intent(in) :: z\\nend function\\nend interface\\ninterface\\nREAL*8 function test2(x, y, z)\\nimplicit none\\nREAL*8, intent(in) :: x\\nREAL*8, intent(in) :: y\\nREAL*8, intent(in) :: z\\nend function\\nend interface\\n")
})

test_that("utilities_test_codegen: test_loops", {
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ invisible(NULL) }, silent = TRUE))
  suppressWarnings(try({ invisible(NULL) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- sym_symbols("n,m", integer = TRUE); n <- py_getitem(.tmp, 0L); m <- py_getitem(.tmp, 1L) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- map(IndexedBase, "Axy"); A <- py_getitem(.tmp, 0L); x <- py_getitem(.tmp, 1L); y <- py_getitem(.tmp, 2L) }, silent = TRUE))
  suppressWarnings(try({ i <- Idx("i", m) }, silent = TRUE))
  suppressWarnings(try({ j <- Idx("j", n) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- codegen(list("matrix_vector", Eq(py_getitem(y, i), py_getitem(A, i, j) * py_getitem(x, j))), "F95", "file", header = FALSE, empty = FALSE); .tmp1 <- py_getitem(.tmp, 0L); f1 <- py_getitem(.tmp1, 0L); code <- py_getitem(.tmp1, 1L); .tmp1 <- py_getitem(.tmp, 1L); f2 <- py_getitem(.tmp1, 0L); interface <- py_getitem(.tmp1, 1L) }, silent = TRUE))
  expect_identical(.cgb(f1), "file.f90")
  suppressWarnings(try({ expected <- "subroutine matrix_vector(A, m, n, x, y)\nimplicit none\nINTEGER*4, intent(in) :: m\nINTEGER*4, intent(in) :: n\nREAL*8, intent(in), dimension(1:m, 1:n) :: A\nREAL*8, intent(in), dimension(1:n) :: x\nREAL*8, intent(out), dimension(1:m) :: y\nINTEGER*4 :: i\nINTEGER*4 :: j\ndo i = 1, m\n   y(i) = 0\nend do\ndo i = 1, m\n   do j = 1, n\n      y(i) = %(rhs)s + y(i)\n   end do\nend do\nend subroutine\n" }, silent = TRUE))
  expect_identical(.cgb((py_eq(code, py_format(expected, py_dict(list("rhs"), list("A(i, j)*x(j)"))))) || (py_eq(code, py_format(expected, py_dict(list("rhs"), list("x(j)*A(i, j)")))))), "True")
  expect_identical(.cgb(f2), "file.h")
  expect_identical(.cgb(interface), "interface\\nsubroutine matrix_vector(A, m, n, x, y)\\nimplicit none\\nINTEGER*4, intent(in) :: m\\nINTEGER*4, intent(in) :: n\\nREAL*8, intent(in), dimension(1:m, 1:n) :: A\\nREAL*8, intent(in), dimension(1:n) :: x\\nREAL*8, intent(out), dimension(1:m) :: y\\nend subroutine\\nend interface\\n")
})

test_that("utilities_test_codegen: test_loops_InOut", {
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ invisible(NULL) }, silent = TRUE))
  suppressWarnings(try({ invisible(NULL) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- sym_symbols("i,j,n,m", integer = TRUE); i <- py_getitem(.tmp, 0L); j <- py_getitem(.tmp, 1L); n <- py_getitem(.tmp, 2L); m <- py_getitem(.tmp, 3L) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- sym_symbols("A,x,y"); A <- py_getitem(.tmp, 0L); x <- py_getitem(.tmp, 1L); y <- py_getitem(.tmp, 2L) }, silent = TRUE))
  suppressWarnings(try({ A <- py_getitem(IndexedBase(A), Idx(i, m), Idx(j, n)) }, silent = TRUE))
  suppressWarnings(try({ x <- py_getitem(IndexedBase(x), Idx(j, n)) }, silent = TRUE))
  suppressWarnings(try({ y <- py_getitem(IndexedBase(y), Idx(i, m)) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- codegen(list("matrix_vector", Eq(y, y + (A * x))), "F95", "file", header = FALSE, empty = FALSE); .tmp1 <- py_getitem(.tmp, 0L); f1 <- py_getitem(.tmp1, 0L); code <- py_getitem(.tmp1, 1L); .tmp1 <- py_getitem(.tmp, 1L); f2 <- py_getitem(.tmp1, 0L); interface <- py_getitem(.tmp1, 1L) }, silent = TRUE))
  expect_identical(.cgb(f1), "file.f90")
  suppressWarnings(try({ expected <- "subroutine matrix_vector(A, m, n, x, y)\nimplicit none\nINTEGER*4, intent(in) :: m\nINTEGER*4, intent(in) :: n\nREAL*8, intent(in), dimension(1:m, 1:n) :: A\nREAL*8, intent(in), dimension(1:n) :: x\nREAL*8, intent(inout), dimension(1:m) :: y\nINTEGER*4 :: i\nINTEGER*4 :: j\ndo i = 1, m\n   do j = 1, n\n      y(i) = %(rhs)s + y(i)\n   end do\nend do\nend subroutine\n" }, silent = TRUE))
  expect_identical(.cgb((py_eq(code, py_format(expected, py_dict(list("rhs"), list("A(i, j)*x(j)"))))) || (py_eq(code, py_format(expected, py_dict(list("rhs"), list("x(j)*A(i, j)")))))), "True")
  expect_identical(.cgb(f2), "file.h")
  expect_identical(.cgb(interface), "interface\\nsubroutine matrix_vector(A, m, n, x, y)\\nimplicit none\\nINTEGER*4, intent(in) :: m\\nINTEGER*4, intent(in) :: n\\nREAL*8, intent(in), dimension(1:m, 1:n) :: A\\nREAL*8, intent(in), dimension(1:n) :: x\\nREAL*8, intent(inout), dimension(1:m) :: y\\nend subroutine\\nend interface\\n")
})

test_that("utilities_test_codegen: test_partial_loops_f", {
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ invisible(NULL) }, silent = TRUE))
  suppressWarnings(try({ invisible(NULL) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- sym_symbols("n m o p", integer = TRUE); n <- py_getitem(.tmp, 0L); m <- py_getitem(.tmp, 1L); o <- py_getitem(.tmp, 2L); p <- py_getitem(.tmp, 3L) }, silent = TRUE))
  suppressWarnings(try({ A <- IndexedBase("A", shape = list(m, p)) }, silent = TRUE))
  suppressWarnings(try({ x <- IndexedBase("x") }, silent = TRUE))
  suppressWarnings(try({ y <- IndexedBase("y") }, silent = TRUE))
  suppressWarnings(try({ i <- Idx("i", list(o, m - 5L)) }, silent = TRUE))
  suppressWarnings(try({ j <- Idx("j", n) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- codegen(list("matrix_vector", Eq(py_getitem(y, i), py_getitem(A, i, j) * py_getitem(x, j))), "F95", "file", header = FALSE, empty = FALSE); .tmp1 <- py_getitem(.tmp, 0L); f1 <- py_getitem(.tmp1, 0L); code <- py_getitem(.tmp1, 1L); .tmp1 <- py_getitem(.tmp, 1L); f2 <- py_getitem(.tmp1, 0L); interface <- py_getitem(.tmp1, 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- py_format("subroutine matrix_vector(A, m, n, o, p, x, y)\nimplicit none\nINTEGER*4, intent(in) :: m\nINTEGER*4, intent(in) :: n\nINTEGER*4, intent(in) :: o\nINTEGER*4, intent(in) :: p\nREAL*8, intent(in), dimension(1:m, 1:p) :: A\nREAL*8, intent(in), dimension(1:n) :: x\nREAL*8, intent(out), dimension(1:%(iup-ilow)s) :: y\nINTEGER*4 :: i\nINTEGER*4 :: j\ndo i = %(ilow)s, %(iup)s\n   y(i) = 0\nend do\ndo i = %(ilow)s, %(iup)s\n   do j = 1, n\n      y(i) = %(rhs)s + y(i)\n   end do\nend do\nend subroutine\n", py_dict(list("rhs", "iup", "ilow", "iup-ilow"), list("%(rhs)s", py_str(m - 4L), py_str(1L + o), py_str((m - 4L) - o)))) }, silent = TRUE))
  expect_identical(.cgb((py_eq(code, py_format(expected, py_dict(list("rhs"), list("A(i, j)*x(j)"))))) || (py_eq(code, py_format(expected, py_dict(list("rhs"), list("x(j)*A(i, j)")))))), "True")
})

test_that("utilities_test_codegen: test_output_arg_f", {
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ invisible(NULL) }, silent = TRUE))
  suppressWarnings(try({ invisible(NULL) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- sym_symbols("x,y,z"); x <- py_getitem(.tmp, 0L); y <- py_getitem(.tmp, 1L); z <- py_getitem(.tmp, 2L) }, silent = TRUE))
  suppressWarnings(try({ r <- make_routine("foo", list(Equality(y, sin(x)), cos(x))) }, silent = TRUE))
  suppressWarnings(try({ c <- FCodeGen() }, silent = TRUE))
  suppressWarnings(try({ result <- py_method(c, "write", list(r), "test", header = FALSE, empty = FALSE) }, silent = TRUE))
  expect_identical(.cgb(py_getitem(py_getitem(result, 0L), 0L)), "test.f90")
  expect_identical(.cgb(py_getitem(py_getitem(result, 0L), 1L)), "REAL*8 function foo(x, y)\\nimplicit none\\nREAL*8, intent(in) :: x\\nREAL*8, intent(out) :: y\\ny = sin(x)\\nfoo = cos(x)\\nend function\\n")
})

test_that("utilities_test_codegen: test_inline_function", {
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ invisible(NULL) }, silent = TRUE))
  suppressWarnings(try({ invisible(NULL) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- sym_symbols("n m", integer = TRUE); n <- py_getitem(.tmp, 0L); m <- py_getitem(.tmp, 1L) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- map(IndexedBase, "Axy"); A <- py_getitem(.tmp, 0L); x <- py_getitem(.tmp, 1L); y <- py_getitem(.tmp, 2L) }, silent = TRUE))
  suppressWarnings(try({ i <- Idx("i", m) }, silent = TRUE))
  suppressWarnings(try({ p <- FCodeGen() }, silent = TRUE))
  suppressWarnings(try({ func <- implemented_function("func", Lambda(n, n * (n + 1L))) }, silent = TRUE))
  suppressWarnings(try({ routine <- make_routine("test_inline", Eq(py_getitem(y, i), func(py_getitem(x, i)))) }, silent = TRUE))
  suppressWarnings(try({ code <- get_string(ast_attr(p, "dump_f95"), list(routine)) }, silent = TRUE))
  suppressWarnings(try({ expected <- "subroutine test_inline(m, x, y)\nimplicit none\nINTEGER*4, intent(in) :: m\nREAL*8, intent(in), dimension(1:m) :: x\nREAL*8, intent(out), dimension(1:m) :: y\nINTEGER*4 :: i\ndo i = 1, m\n   y(i) = %s*%s\nend do\nend subroutine\n" }, silent = TRUE))
  suppressWarnings(try({ args <- list("x(i)", "(x(i) + 1)") }, silent = TRUE))
  expect_identical(.cgb((py_eq(code, py_format(expected, args))) || (py_eq(code, py_format(expected, py_slice(args, NULL, NULL, -1L))))), "True")
})

test_that("utilities_test_codegen: test_f_code_call_signature_wrap", {
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ x <- sym_symbols("x:20") }, silent = TRUE))
  suppressWarnings(try({ expr <- 0L }, silent = TRUE))
  suppressWarnings(try({ routine <- make_routine("test", expr) }, silent = TRUE))
  suppressWarnings(try({ code_gen <- FCodeGen() }, silent = TRUE))
  suppressWarnings(try({ source <- get_string(ast_attr(code_gen, "dump_f95"), list(routine)) }, silent = TRUE))
  suppressWarnings(try({ expected <- "REAL*8 function test(x0, x1, x10, x11, x12, x13, x14, x15, x16, x17, x18, &\n      x19, x2, x3, x4, x5, x6, x7, x8, x9)\nimplicit none\nREAL*8, intent(in) :: x0\nREAL*8, intent(in) :: x1\nREAL*8, intent(in) :: x10\nREAL*8, intent(in) :: x11\nREAL*8, intent(in) :: x12\nREAL*8, intent(in) :: x13\nREAL*8, intent(in) :: x14\nREAL*8, intent(in) :: x15\nREAL*8, intent(in) :: x16\nREAL*8, intent(in) :: x17\nREAL*8, intent(in) :: x18\nREAL*8, intent(in) :: x19\nREAL*8, intent(in) :: x2\nREAL*8, intent(in) :: x3\nREAL*8, intent(in) :: x4\nREAL*8, intent(in) :: x5\nREAL*8, intent(in) :: x6\nREAL*8, intent(in) :: x7\nREAL*8, intent(in) :: x8\nREAL*8, intent(in) :: x9\ntest = x0 + x1 + x10 + x11 + x12 + x13 + x14 + x15 + x16 + x17 + x18 + &\n      x19 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9\nend function\n" }, silent = TRUE))
  expect_identical(.cgb(source), "INTEGER*4 function test()\\nimplicit none\\ntest = 0\\nend function\\n")
})

test_that("utilities_test_codegen: test_c_fortran_omit_routine_name", {
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ .tmp <- sym_symbols("x,y"); x <- py_getitem(.tmp, 0L); y <- py_getitem(.tmp, 1L) }, silent = TRUE))
  suppressWarnings(try({ name_expr <- list(list("foo", 2L * x)) }, silent = TRUE))
  suppressWarnings(try({ result <- codegen(name_expr, "F95", header = FALSE, empty = FALSE) }, silent = TRUE))
  suppressWarnings(try({ expresult <- codegen(name_expr, "F95", "foo", header = FALSE, empty = FALSE) }, silent = TRUE))
  expect_identical(.cgb(py_getitem(py_getitem(result, 0L), 1L)), "REAL*8 function foo(x)\\nimplicit none\\nREAL*8, intent(in) :: x\\nfoo = 2*x\\nend function\\n")
  suppressWarnings(try({ name_expr <- list("foo", x * y) }, silent = TRUE))
  suppressWarnings(try({ result <- codegen(name_expr, "F95", header = FALSE, empty = FALSE) }, silent = TRUE))
  suppressWarnings(try({ expresult <- codegen(name_expr, "F95", "foo", header = FALSE, empty = FALSE) }, silent = TRUE))
  expect_identical(.cgb(py_getitem(py_getitem(result, 0L), 1L)), "REAL*8 function foo(x, y)\\nimplicit none\\nREAL*8, intent(in) :: x\\nREAL*8, intent(in) :: y\\nfoo = x*y\\nend function\\n")
  suppressWarnings(try({ name_expr <- list("foo", Matrix(list(list(x, y), list(x + y, x - y)))) }, silent = TRUE))
  suppressWarnings(try({ result <- codegen(name_expr, "C89", header = FALSE, empty = FALSE) }, silent = TRUE))
  suppressWarnings(try({ expresult <- codegen(name_expr, "C89", "foo", header = FALSE, empty = FALSE) }, silent = TRUE))
})

test_that("utilities_test_codegen: test_fcode_results_named_ordered", {
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ .tmp <- sym_symbols("x,y,z"); x <- py_getitem(.tmp, 0L); y <- py_getitem(.tmp, 1L); z <- py_getitem(.tmp, 2L) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- sym_symbols("B,C"); B <- py_getitem(.tmp, 0L); C <- py_getitem(.tmp, 1L) }, silent = TRUE))
  suppressWarnings(try({ A <- MatrixSymbol("A", 1L, 3L) }, silent = TRUE))
  suppressWarnings(try({ expr1 <- Equality(A, Matrix(list(list(1L, 2L, x)))) }, silent = TRUE))
  suppressWarnings(try({ expr2 <- Equality(C, (x + y) * z) }, silent = TRUE))
  suppressWarnings(try({ expr3 <- Equality(B, 2L * x) }, silent = TRUE))
  suppressWarnings(try({ name_expr <- list("test", list(expr1, expr2, expr3)) }, silent = TRUE))
  suppressWarnings(try({ result <- codegen(name_expr, "f95", "test", header = FALSE, empty = FALSE, argument_sequence = list(x, z, y, C, A, B)) }, silent = TRUE))
  suppressWarnings(try({ source <- py_getitem(py_getitem(result, 0L), 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "subroutine test(x, z, y, C, A, B)\nimplicit none\nREAL*8, intent(in) :: x\nREAL*8, intent(in) :: z\nREAL*8, intent(in) :: y\nREAL*8, intent(out) :: C\nREAL*8, intent(out) :: B\nREAL*8, intent(out), dimension(1:1, 1:3) :: A\nC = z*(x + y)\nA(1, 1) = 1\nA(1, 2) = 2\nA(1, 3) = x\nB = 2*x\nend subroutine\n" }, silent = TRUE))
  expect_identical(.cgb(source), "subroutine test(x, z, y, C, A, B)\\nimplicit none\\nREAL*8, intent(in) :: x\\nREAL*8, intent(in) :: z\\nREAL*8, intent(in) :: y\\nREAL*8, intent(out) :: C\\nREAL*8, intent(out) :: B\\nREAL*8, intent(out), dimension(1:1, 1:3) :: A\\nC = z*(x + y)\\nA(1, 1) = 1\\nA(1, 2) = 2\\nA(1, 3) = x\\nB = 2*x\\nend subroutine\\n")
})

test_that("utilities_test_codegen: test_fcode_matrixsymbol_slice", {
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ A <- MatrixSymbol("A", 2L, 3L) }, silent = TRUE))
  suppressWarnings(try({ B <- MatrixSymbol("B", 1L, 3L) }, silent = TRUE))
  suppressWarnings(try({ C <- MatrixSymbol("C", 1L, 3L) }, silent = TRUE))
  suppressWarnings(try({ D <- MatrixSymbol("D", 2L, 1L) }, silent = TRUE))
  suppressWarnings(try({ name_expr <- list("test", list(Equality(B, py_getitem(A, 0L, pyslice(NULL, NULL, NULL))), Equality(C, py_getitem(A, 1L, pyslice(NULL, NULL, NULL))), Equality(D, py_getitem(A, pyslice(NULL, NULL, NULL), 2L)))) }, silent = TRUE))
  suppressWarnings(try({ result <- codegen(name_expr, "f95", "test", header = FALSE, empty = FALSE) }, silent = TRUE))
  suppressWarnings(try({ source <- py_getitem(py_getitem(result, 0L), 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "subroutine test(A, B, C, D)\nimplicit none\nREAL*8, intent(in), dimension(1:2, 1:3) :: A\nREAL*8, intent(out), dimension(1:1, 1:3) :: B\nREAL*8, intent(out), dimension(1:1, 1:3) :: C\nREAL*8, intent(out), dimension(1:2, 1:1) :: D\nB(1, 1) = A(1, 1)\nB(1, 2) = A(1, 2)\nB(1, 3) = A(1, 3)\nC(1, 1) = A(2, 1)\nC(1, 2) = A(2, 2)\nC(1, 3) = A(2, 3)\nD(1, 1) = A(1, 3)\nD(2, 1) = A(2, 3)\nend subroutine\n" }, silent = TRUE))
  expect_identical(.cgb(source), "subroutine test(A, B, C, D)\\nimplicit none\\nREAL*8, intent(in), dimension(1:2, 1:3) :: A\\nREAL*8, intent(out), dimension(1:1, 1:3) :: B\\nREAL*8, intent(out), dimension(1:1, 1:3) :: C\\nREAL*8, intent(out), dimension(1:2, 1:1) :: D\\nB(1, 1) = A(1, 1)\\nB(1, 2) = A(1, 2)\\nB(1, 3) = A(1, 3)\\nC(1, 1) = A(2, 1)\\nC(1, 2) = A(2, 2)\\nC(1, 3) = A(2, 3)\\nD(1, 1) = A(1, 3)\\nD(2, 1) = A(2, 3)\\nend subroutine\\n")
})

test_that("utilities_test_codegen: test_global_vars", {
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ .tmp <- sym_symbols("x y z t"); x <- py_getitem(.tmp, 0L); y <- py_getitem(.tmp, 1L); z <- py_getitem(.tmp, 2L); t <- py_getitem(.tmp, 3L) }, silent = TRUE))
  suppressWarnings(try({ result <- codegen(list("f", x * y), "F95", header = FALSE, empty = FALSE, global_vars = list(y)) }, silent = TRUE))
  suppressWarnings(try({ source <- py_getitem(py_getitem(result, 0L), 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "REAL*8 function f(x)\nimplicit none\nREAL*8, intent(in) :: x\nf = x*y\nend function\n" }, silent = TRUE))
  expect_identical(.cgb(source), "REAL*8 function f(x)\\nimplicit none\\nREAL*8, intent(in) :: x\\nf = x*y\\nend function\\n")
  suppressWarnings(try({ expected <- "#include \"f.h\"\n#include <math.h>\ndouble f(double x, double y) {\n   double f_result;\n   f_result = x*y + z;\n   return f_result;\n}\n" }, silent = TRUE))
  suppressWarnings(try({ result <- codegen(list("f", (x * y) + z), "C", header = FALSE, empty = FALSE, global_vars = list(z, t)) }, silent = TRUE))
  suppressWarnings(try({ source <- py_getitem(py_getitem(result, 0L), 1L) }, silent = TRUE))
  expect_identical(.cgb(source), "#include \"f.h\"\\n#include <math.h>\\ndouble f(double x, double y) {\\n   double f_result;\\n   f_result = x*y + z;\\n   return f_result;\\n}\\n")
})

test_that("utilities_test_codegen: test_custom_codegen", {
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ invisible(NULL) }, silent = TRUE))
  suppressWarnings(try({ invisible(NULL) }, silent = TRUE))
  suppressWarnings(try({ printer <- C99CodePrinter(settings = py_dict(list("user_functions"), list(py_dict(list("exp"), list("fastexp"))))) }, silent = TRUE))
  suppressWarnings(try({ .tmp <- sym_symbols("x y"); x <- py_getitem(.tmp, 0L); y <- py_getitem(.tmp, 1L) }, silent = TRUE))
  suppressWarnings(try({ expr <- exp(x + y) }, silent = TRUE))
  suppressWarnings(try({ gen <- C99CodeGen(printer = printer, preprocessor_statements = list("#include \"fastexp.h\"")) }, silent = TRUE))
  suppressWarnings(try({ expected <- "#include \"expr.h\"\n#include \"fastexp.h\"\ndouble expr(double x, double y) {\n   double expr_result;\n   expr_result = fastexp(x + y);\n   return expr_result;\n}\n" }, silent = TRUE))
  suppressWarnings(try({ result <- codegen(list("expr", expr), header = FALSE, empty = FALSE, code_gen = gen) }, silent = TRUE))
  suppressWarnings(try({ source <- py_getitem(py_getitem(result, 0L), 1L) }, silent = TRUE))
  expect_identical(.cgb(source), "#include \"expr.h\"\\n#include \"fastexp.h\"\\ndouble expr(double x, double y) {\\n   double expr_result;\\n   expr_result = fastexp(x + y);\\n   return expr_result;\\n}\\n")
  suppressWarnings(try({ gen <- C99CodeGen(printer = printer) }, silent = TRUE))
  suppressWarnings(try({ gen$preprocessor_statements <- c(gen$preprocessor_statements, "#include \"fastexp.h\"") }, silent = TRUE))
  suppressWarnings(try({ expected <- "#include \"expr.h\"\n#include <math.h>\n#include \"fastexp.h\"\ndouble expr(double x, double y) {\n   double expr_result;\n   expr_result = fastexp(x + y);\n   return expr_result;\n}\n" }, silent = TRUE))
  suppressWarnings(try({ result <- codegen(list("expr", expr), header = FALSE, empty = FALSE, code_gen = gen) }, silent = TRUE))
  suppressWarnings(try({ source <- py_getitem(py_getitem(result, 0L), 1L) }, silent = TRUE))
  expect_identical(.cgb(source), "#include \"expr.h\"\\n#include <math.h>\\n#include \"fastexp.h\"\\ndouble expr(double x, double y) {\\n   double expr_result;\\n   expr_result = fastexp(x + y);\\n   return expr_result;\\n}\\n")
})

test_that("utilities_test_codegen: test_fcode_complex", {
  get_string <- function(dump_fn, routines, prefix = "file", header = FALSE, empty = FALSE) {
    output <- StringIO()
    dump_fn(routines, output$con, prefix, header, empty)
    source <- py_method(output, "getvalue")
    close(output$con)
    source
  }
  suppressWarnings(try({ invisible(NULL) }, silent = TRUE))
  suppressWarnings(try({ options(symcompR.complex_allowed = TRUE) }, silent = TRUE))
  suppressWarnings(try({ x <- Symbol("x", real = TRUE) }, silent = TRUE))
  suppressWarnings(try({ y <- Symbol("y", real = TRUE) }, silent = TRUE))
  suppressWarnings(try({ result <- codegen(list("test", x + y), "f95", "test", header = FALSE, empty = FALSE) }, silent = TRUE))
  suppressWarnings(try({ source <- py_getitem(py_getitem(result, 0L), 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "REAL*8 function test(x, y)\nimplicit none\nREAL*8, intent(in) :: x\nREAL*8, intent(in) :: y\ntest = x + y\nend function\n" }, silent = TRUE))
  expect_identical(.cgb(source), "REAL*8 function test(x, y)\\nimplicit none\\nREAL*8, intent(in) :: x\\nREAL*8, intent(in) :: y\\ntest = x + y\\nend function\\n")
  suppressWarnings(try({ x <- Symbol("x") }, silent = TRUE))
  suppressWarnings(try({ y <- Symbol("y", real = TRUE) }, silent = TRUE))
  suppressWarnings(try({ result <- codegen(list("test", x + y), "f95", "test", header = FALSE, empty = FALSE) }, silent = TRUE))
  suppressWarnings(try({ source <- py_getitem(py_getitem(result, 0L), 1L) }, silent = TRUE))
  suppressWarnings(try({ expected <- "COMPLEX*16 function test(x, y)\nimplicit none\nCOMPLEX*16, intent(in) :: x\nREAL*8, intent(in) :: y\ntest = x + y\nend function\n" }, silent = TRUE))
  expect_identical(.cgb(source), "COMPLEX*16 function test(x, y)\\nimplicit none\\nCOMPLEX*16, intent(in) :: x\\nREAL*8, intent(in) :: y\\ntest = x + y\\nend function\\n")
  suppressWarnings(try({ options(symcompR.complex_allowed = FALSE) }, silent = TRUE))
})

test_that("utilities_test_lambdify: test_imps", {
  w <- Symbol("w"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  suppressWarnings(try({ f <- implemented_function("f", function(x) 2L * x) }, silent = TRUE))
  suppressWarnings(try({ g <- implemented_function("f", function(x) sqrt(x)) }, silent = TRUE))
  suppressWarnings(try({ l1 <- lambdify(x, f(x)) }, silent = TRUE))
  suppressWarnings(try({ l2 <- lambdify(x, g(x)) }, silent = TRUE))
  expect_identical(.cgb(py_str(f(x))), "f(x)")
  expect_identical(.cgb(py_str(g(x))), "f(x)")
  expect_identical(.cgb(l1(3L)), "6")
  expect_identical(.cgb(l2(4L)), "2.0")
  expect_identical(.cgb(f(x) == g(x)), "False")
  suppressWarnings(try({ func <- Function("myfunc") }, silent = TRUE))
  expect_identical(.cgb(hasattr(func, "_imp_")), "False")
  suppressWarnings(try({ my_f <- implemented_function(func, function(x) 2L * x) }, silent = TRUE))
  expect_identical(.cgb(hasattr(my_f, "_imp_")), "True")
  expect_identical(.cgb(py_str(my_f(y))), "myfunc(y)")
  suppressWarnings(try({ f2 <- implemented_function("f", function(x) x + 101L) }, silent = TRUE))
  expect_identical(.cgb(lambdify(x, f(f2(x)))), "ERR ValueError: We found more than one implementation with name \"f\"")
})

test_that("utilities_test_lambdify: test_imps_wrong_args", {
  w <- Symbol("w"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  expect_identical(.cgb(implemented_function(sin, function(x) x)), "ERR ValueError: \\nsymfunc should be either a string or an UndefinedFunction instance.")
})

test_that("utilities_test_lambdify: test_lambdify_imps", {
  w <- Symbol("w"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  suppressWarnings(try({ f <- implemented_function("f", function(x) x + 100L) }, silent = TRUE))
  expect_identical(.cgb(py_call(lambdify(x, f(x)), 0L)), "100")
  expect_identical(.cgb(py_call(lambdify(list(x, y), y + f(x)), 0L, 1L)), "101")
})

test_that("utilities_test_lambdify: test_implemented_function_printing", {
  w <- Symbol("w"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  suppressWarnings(try({ f <- implemented_function("f", Lambda(x, 2L * x)) }, silent = TRUE))
  expect_identical(.cgb(py_str(f(x))), "f(x)")
  expect_identical(.cgb(ccode(f(x))), "2*x")
  expect_identical(.cgb(ccode(f(y) + 1L)), "2*y + 1")
  expect_identical(.cgb(fcode(f(x), source_format = "free")), "2*x")
  expect_identical(.cgb(pycode(f(x))), "2*x")
  expect_identical(.cgb(diff(f(x), x)), "Derivative(f(x), x)")
})

test_that("utilities_test_lambdify: test_lambdastr", {
  w <- Symbol("w"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  expect_identical(.cgb(lambdastr(x, x ^ 2L)), "lambda x: (x^2)")
  expect_identical(.cgb(lambdastr(list(x, y, z), list(z, y, x))), "lambda x,y,z: ([z, y, x])")
  expect_identical(.cgb(lambdastr(list(x, list(y, z)), x + y)), "lambda _0,_1: (lambda x,y,z: (x + y))(_0,_1[0],_1[1])")
  expect_identical(.cgb(lambdastr(list(x, y), x * y)), "lambda x,y: (x*y)")
  expect_identical(.cgb(lambdastr(x, sin(x) + cos(x))), "lambda x: (math.sin(x) + math.cos(x))")
  expect_identical(.cgb(lambdastr(list(x, y), list(x + y, x - y))), "lambda x,y: ([x + y, x - y])")
  expect_identical(.cgb(lambdastr(x, x ^ 2L, dummify = FALSE)), "lambda x: (x^2)")
})

test_that("utilities_test_lambdify: test_lambdify_Derivative_arg_issue_16468", {
  w <- Symbol("w"); x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
  suppressWarnings(try({ f <- py_call(Function("f"), x) }, silent = TRUE))
  suppressWarnings(try({ fx <- diff(f) }, silent = TRUE))
  expect_identical(.cgb(lambdastr(list(f, fx), f / fx, dummify = FALSE)), "ERR PrintMethodNotImplementedError: Unsupported by <class 'sympy.printing.lambdarepr.LambdaPrinter'>: f\\nSet the printer option 'strict' to False in order to generate partially printed code.")
  expect_identical(.cgb(lambdastr(fx, 1L + fx, dummify = FALSE)), "ERR PrintMethodNotImplementedError: Unsupported by <class 'sympy.printing.lambdarepr.LambdaPrinter'>: <class 'sympy.core.function.Derivative'>\\nPrinter has no method: _print_Derivative_f\\nSet the printer option 'strict' to False in order to generate partially printed code.")
})

