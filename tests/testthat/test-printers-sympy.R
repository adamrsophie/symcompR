skip_unless_full_tests()

# Expected values are SymPy outputs for the tests of SymPy in
# sympy/printing/tests (test_numpy.py, test_tensorflow.py, test_torch.py,
# test_lambdarepr.py, the array printer rows of test_pycode.py) and
# sympy/strategies/branch/tests: each row prints str(result) with ** shown
# as ^, or the error class and message.

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

test_that("printing_test_numpy: test_numpy_piecewise_regression", {
  x <- Symbol("x"); i <- Symbol("i"); j <- Symbol("j"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d")
  U_numpy_known_constants <- .cg_prefixed(.cg_kc_numpy, "numpy.")
  U_numpy_known_functions <- .cg_prefixed(.cg_kf_numpy, "numpy.")
  U_scipy_known_constants <- .cg_prefixed(.cg_kc_scipy, "scipy.constants.")
  U_scipy_known_functions <- .cg_prefixed(.cg_kf_scipy, "scipy.special.")
  mi_sorted <- function(p) {
    mi <- unclass(ast_attr(p, "module_imports"))
    pystr(paste0("{", paste(vapply(names(mi), function(k) sprintf("'%s': [%s]", k, paste0("'", sort(mi[[k]], method = "radix"), "'", collapse = ", ")), ""), collapse = ", "), "}"))
  }
  suppressWarnings(try({ printer <- NumPyPrinter() }, silent = TRUE))
  suppressWarnings(try({ p <- Piecewise(list(1L, x < 0L), list(0L, TRUE)) }, silent = TRUE))
  expect_identical(.cgb(printer_doprint(printer, p)), "numpy.select([numpy.less(x, 0),True], [1,0], default=numpy.nan)")
  expect_identical(.cgb(mi_sorted(printer)), "{'numpy': ['less', 'nan', 'select']}")
})

test_that("printing_test_numpy: test_numpy_logaddexp", {
  x <- Symbol("x"); i <- Symbol("i"); j <- Symbol("j"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d")
  U_numpy_known_constants <- .cg_prefixed(.cg_kc_numpy, "numpy.")
  U_numpy_known_functions <- .cg_prefixed(.cg_kf_numpy, "numpy.")
  U_scipy_known_constants <- .cg_prefixed(.cg_kc_scipy, "scipy.constants.")
  U_scipy_known_functions <- .cg_prefixed(.cg_kf_scipy, "scipy.special.")
  mi_sorted <- function(p) {
    mi <- unclass(ast_attr(p, "module_imports"))
    pystr(paste0("{", paste(vapply(names(mi), function(k) sprintf("'%s': [%s]", k, paste0("'", sort(mi[[k]], method = "radix"), "'", collapse = ", ")), ""), collapse = ", "), "}"))
  }
  suppressWarnings(try({ lae <- logaddexp(a, b) }, silent = TRUE))
  expect_identical(.cgb(printer_doprint(NumPyPrinter(), lae)), "numpy.logaddexp(a, b)")
  suppressWarnings(try({ lae2 <- logaddexp2(a, b) }, silent = TRUE))
  expect_identical(.cgb(printer_doprint(NumPyPrinter(), lae2)), "numpy.logaddexp2(a, b)")
})

test_that("printing_test_numpy: test_numpy_conjugate", {
  x <- Symbol("x"); i <- Symbol("i"); j <- Symbol("j"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d")
  U_numpy_known_constants <- .cg_prefixed(.cg_kc_numpy, "numpy.")
  U_numpy_known_functions <- .cg_prefixed(.cg_kf_numpy, "numpy.")
  U_scipy_known_constants <- .cg_prefixed(.cg_kc_scipy, "scipy.constants.")
  U_scipy_known_functions <- .cg_prefixed(.cg_kf_scipy, "scipy.special.")
  mi_sorted <- function(p) {
    mi <- unclass(ast_attr(p, "module_imports"))
    pystr(paste0("{", paste(vapply(names(mi), function(k) sprintf("'%s': [%s]", k, paste0("'", sort(mi[[k]], method = "radix"), "'", collapse = ", ")), ""), collapse = ", "), "}"))
  }
  expect_identical(.cgb(printer_doprint(NumPyPrinter(), conjugate(x))), "numpy.conjugate(x)")
})

test_that("printing_test_numpy: test_sum", {
  x <- Symbol("x"); i <- Symbol("i"); j <- Symbol("j"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d")
  U_numpy_known_constants <- .cg_prefixed(.cg_kc_numpy, "numpy.")
  U_numpy_known_functions <- .cg_prefixed(.cg_kf_numpy, "numpy.")
  U_scipy_known_constants <- .cg_prefixed(.cg_kc_scipy, "scipy.constants.")
  U_scipy_known_functions <- .cg_prefixed(.cg_kf_scipy, "scipy.special.")
  mi_sorted <- function(p) {
    mi <- unclass(ast_attr(p, "module_imports"))
    pystr(paste0("{", paste(vapply(names(mi), function(k) sprintf("'%s': [%s]", k, paste0("'", sort(mi[[k]], method = "radix"), "'", collapse = ", ")), ""), collapse = ", "), "}"))
  }
  suppressWarnings(try({ s <- Sum(x ^ i, list(i, a, b)) }, silent = TRUE))
  expect_identical(.cgb(printer_doprint(NumPyPrinter(), s)), "(builtins.sum(x^i for i in range(a, b+1)))")
  suppressWarnings(try({ s <- Sum(i * x, list(i, a, b)) }, silent = TRUE))
  expect_identical(.cgb(printer_doprint(NumPyPrinter(), s)), "(builtins.sum(i*x for i in range(a, b+1)))")
})

test_that("printing_test_numpy: test_multiple_sums", {
  x <- Symbol("x"); i <- Symbol("i"); j <- Symbol("j"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d")
  U_numpy_known_constants <- .cg_prefixed(.cg_kc_numpy, "numpy.")
  U_numpy_known_functions <- .cg_prefixed(.cg_kf_numpy, "numpy.")
  U_scipy_known_constants <- .cg_prefixed(.cg_kc_scipy, "scipy.constants.")
  U_scipy_known_functions <- .cg_prefixed(.cg_kf_scipy, "scipy.special.")
  mi_sorted <- function(p) {
    mi <- unclass(ast_attr(p, "module_imports"))
    pystr(paste0("{", paste(vapply(names(mi), function(k) sprintf("'%s': [%s]", k, paste0("'", sort(mi[[k]], method = "radix"), "'", collapse = ", ")), ""), collapse = ", "), "}"))
  }
  suppressWarnings(try({ s <- Sum((x + j) * i, list(i, a, b), list(j, c, d)) }, silent = TRUE))
  expect_identical(.cgb(printer_doprint(NumPyPrinter(), s)), "(builtins.sum(i*(j + x) for j in range(c, d+1) for i in range(a, b+1)))")
})

test_that("printing_test_numpy: test_codegen_einsum", {
  x <- Symbol("x"); i <- Symbol("i"); j <- Symbol("j"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d")
  U_numpy_known_constants <- .cg_prefixed(.cg_kc_numpy, "numpy.")
  U_numpy_known_functions <- .cg_prefixed(.cg_kf_numpy, "numpy.")
  U_scipy_known_constants <- .cg_prefixed(.cg_kc_scipy, "scipy.constants.")
  U_scipy_known_functions <- .cg_prefixed(.cg_kf_scipy, "scipy.special.")
  mi_sorted <- function(p) {
    mi <- unclass(ast_attr(p, "module_imports"))
    pystr(paste0("{", paste(vapply(names(mi), function(k) sprintf("'%s': [%s]", k, paste0("'", sort(mi[[k]], method = "radix"), "'", collapse = ", ")), ""), collapse = ", "), "}"))
  }
  suppressWarnings(try({ M <- MatrixSymbol("M", 2L, 2L) }, silent = TRUE))
  suppressWarnings(try({ N <- MatrixSymbol("N", 2L, 2L) }, silent = TRUE))
  suppressWarnings(try({ cg <- convert_matrix_to_array(M * N) }, silent = TRUE))
  expect_identical(.cgb(cg), "ArrayContraction(ArrayTensorProduct(M, N), (1, 2))")
  expect_identical(.cgb(printer_doprint(NumPyPrinter(), cg)), "numpy.einsum(\"ab,bc->ac\", M,N)")
})

test_that("printing_test_numpy: test_codegen_extra", {
  x <- Symbol("x"); i <- Symbol("i"); j <- Symbol("j"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d")
  U_numpy_known_constants <- .cg_prefixed(.cg_kc_numpy, "numpy.")
  U_numpy_known_functions <- .cg_prefixed(.cg_kf_numpy, "numpy.")
  U_scipy_known_constants <- .cg_prefixed(.cg_kc_scipy, "scipy.constants.")
  U_scipy_known_functions <- .cg_prefixed(.cg_kf_scipy, "scipy.special.")
  mi_sorted <- function(p) {
    mi <- unclass(ast_attr(p, "module_imports"))
    pystr(paste0("{", paste(vapply(names(mi), function(k) sprintf("'%s': [%s]", k, paste0("'", sort(mi[[k]], method = "radix"), "'", collapse = ", ")), ""), collapse = ", "), "}"))
  }
  suppressWarnings(try({ M <- MatrixSymbol("M", 2L, 2L) }, silent = TRUE))
  suppressWarnings(try({ N <- MatrixSymbol("N", 2L, 2L) }, silent = TRUE))
  suppressWarnings(try({ P <- MatrixSymbol("P", 2L, 2L) }, silent = TRUE))
  suppressWarnings(try({ Q <- MatrixSymbol("Q", 2L, 2L) }, silent = TRUE))
  expect_identical(.cgb(printer_doprint(NumPyPrinter(), ArrayTensorProduct(M, N))), "numpy.einsum(\"ab,cd\", M, N)")
  expect_identical(.cgb(printer_doprint(NumPyPrinter(), ArrayAdd(M, N))), "numpy.add(M, N)")
  expect_identical(.cgb(printer_doprint(NumPyPrinter(), ArrayAdd(M, N, P))), "numpy.add(numpy.add(M, N), P)")
  expect_identical(.cgb(printer_doprint(NumPyPrinter(), ArrayAdd(M, N, P, Q))), "numpy.add(numpy.add(numpy.add(M, N), P), Q)")
  expect_identical(.cgb(printer_doprint(NumPyPrinter(), PermuteDims(M, list(1L, 0L)))), "numpy.transpose(M, [1, 0])")
  expect_identical(.cgb(printer_doprint(NumPyPrinter(), PermuteDims(ArrayTensorProduct(M, N), list(1L, 2L, 3L, 0L)))), "numpy.transpose(numpy.einsum(\"ab,cd\", M, N), [1, 2, 3, 0])")
  expect_identical(.cgb(printer_doprint(NumPyPrinter(), ArrayDiagonal(ArrayTensorProduct(M, N), list(1L, 2L)))), "numpy.einsum(\"ab,bc->acb\", M, N)")
})

test_that("printing_test_numpy: test_relational", {
  x <- Symbol("x"); i <- Symbol("i"); j <- Symbol("j"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d")
  U_numpy_known_constants <- .cg_prefixed(.cg_kc_numpy, "numpy.")
  U_numpy_known_functions <- .cg_prefixed(.cg_kf_numpy, "numpy.")
  U_scipy_known_constants <- .cg_prefixed(.cg_kc_scipy, "scipy.constants.")
  U_scipy_known_functions <- .cg_prefixed(.cg_kf_scipy, "scipy.special.")
  mi_sorted <- function(p) {
    mi <- unclass(ast_attr(p, "module_imports"))
    pystr(paste0("{", paste(vapply(names(mi), function(k) sprintf("'%s': [%s]", k, paste0("'", sort(mi[[k]], method = "radix"), "'", collapse = ", ")), ""), collapse = ", "), "}"))
  }
  expect_identical(.cgb(printer_doprint(NumPyPrinter(), Equality(x, 1L))), "numpy.equal(x, 1)")
  expect_identical(.cgb(printer_doprint(NumPyPrinter(), Unequality(x, 1L))), "numpy.not_equal(x, 1)")
  expect_identical(.cgb(printer_doprint(NumPyPrinter(), x < 1L)), "numpy.less(x, 1)")
  expect_identical(.cgb(printer_doprint(NumPyPrinter(), x <= 1L)), "numpy.less_equal(x, 1)")
  expect_identical(.cgb(printer_doprint(NumPyPrinter(), x > 1L)), "numpy.greater(x, 1)")
  expect_identical(.cgb(printer_doprint(NumPyPrinter(), x >= 1L)), "numpy.greater_equal(x, 1)")
})

test_that("printing_test_numpy: test_mod", {
  x <- Symbol("x"); i <- Symbol("i"); j <- Symbol("j"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d")
  U_numpy_known_constants <- .cg_prefixed(.cg_kc_numpy, "numpy.")
  U_numpy_known_functions <- .cg_prefixed(.cg_kf_numpy, "numpy.")
  U_scipy_known_constants <- .cg_prefixed(.cg_kc_scipy, "scipy.constants.")
  U_scipy_known_functions <- .cg_prefixed(.cg_kf_scipy, "scipy.special.")
  mi_sorted <- function(p) {
    mi <- unclass(ast_attr(p, "module_imports"))
    pystr(paste0("{", paste(vapply(names(mi), function(k) sprintf("'%s': [%s]", k, paste0("'", sort(mi[[k]], method = "radix"), "'", collapse = ", ")), ""), collapse = ", "), "}"))
  }
  expect_identical(.cgb(printer_doprint(NumPyPrinter(), sym_Mod(a, b))), "numpy.mod(a, b)")
})

test_that("printing_test_numpy: test_pow", {
  x <- Symbol("x"); i <- Symbol("i"); j <- Symbol("j"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d")
  U_numpy_known_constants <- .cg_prefixed(.cg_kc_numpy, "numpy.")
  U_numpy_known_functions <- .cg_prefixed(.cg_kf_numpy, "numpy.")
  U_scipy_known_constants <- .cg_prefixed(.cg_kc_scipy, "scipy.constants.")
  U_scipy_known_functions <- .cg_prefixed(.cg_kf_scipy, "scipy.special.")
  mi_sorted <- function(p) {
    mi <- unclass(ast_attr(p, "module_imports"))
    pystr(paste0("{", paste(vapply(names(mi), function(k) sprintf("'%s': [%s]", k, paste0("'", sort(mi[[k]], method = "radix"), "'", collapse = ", ")), ""), collapse = ", "), "}"))
  }
  expect_identical(.cgb(printer_doprint(NumPyPrinter(), new_sym("Pow", list(S(2L), S(-1L))))), "2^(-1.0)")
})

test_that("printing_test_numpy: test_expm1", {
  x <- Symbol("x"); i <- Symbol("i"); j <- Symbol("j"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d")
  U_numpy_known_constants <- .cg_prefixed(.cg_kc_numpy, "numpy.")
  U_numpy_known_functions <- .cg_prefixed(.cg_kf_numpy, "numpy.")
  U_scipy_known_constants <- .cg_prefixed(.cg_kc_scipy, "scipy.constants.")
  U_scipy_known_functions <- .cg_prefixed(.cg_kf_scipy, "scipy.special.")
  mi_sorted <- function(p) {
    mi <- unclass(ast_attr(p, "module_imports"))
    pystr(paste0("{", paste(vapply(names(mi), function(k) sprintf("'%s': [%s]", k, paste0("'", sort(mi[[k]], method = "radix"), "'", collapse = ", ")), ""), collapse = ", "), "}"))
  }
  expect_identical(.cgb(printer_doprint(NumPyPrinter(), sym_expm1(a))), "numpy.expm1(a)")
})

test_that("printing_test_numpy: test_log1p", {
  x <- Symbol("x"); i <- Symbol("i"); j <- Symbol("j"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d")
  U_numpy_known_constants <- .cg_prefixed(.cg_kc_numpy, "numpy.")
  U_numpy_known_functions <- .cg_prefixed(.cg_kf_numpy, "numpy.")
  U_scipy_known_constants <- .cg_prefixed(.cg_kc_scipy, "scipy.constants.")
  U_scipy_known_functions <- .cg_prefixed(.cg_kf_scipy, "scipy.special.")
  mi_sorted <- function(p) {
    mi <- unclass(ast_attr(p, "module_imports"))
    pystr(paste0("{", paste(vapply(names(mi), function(k) sprintf("'%s': [%s]", k, paste0("'", sort(mi[[k]], method = "radix"), "'", collapse = ", ")), ""), collapse = ", "), "}"))
  }
  expect_identical(.cgb(printer_doprint(NumPyPrinter(), sym_log1p(a))), "numpy.log1p(a)")
})

test_that("printing_test_numpy: test_hypot", {
  x <- Symbol("x"); i <- Symbol("i"); j <- Symbol("j"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d")
  U_numpy_known_constants <- .cg_prefixed(.cg_kc_numpy, "numpy.")
  U_numpy_known_functions <- .cg_prefixed(.cg_kf_numpy, "numpy.")
  U_scipy_known_constants <- .cg_prefixed(.cg_kc_scipy, "scipy.constants.")
  U_scipy_known_functions <- .cg_prefixed(.cg_kf_scipy, "scipy.special.")
  mi_sorted <- function(p) {
    mi <- unclass(ast_attr(p, "module_imports"))
    pystr(paste0("{", paste(vapply(names(mi), function(k) sprintf("'%s': [%s]", k, paste0("'", sort(mi[[k]], method = "radix"), "'", collapse = ", ")), ""), collapse = ", "), "}"))
  }
  expect_identical(.cgb(printer_doprint(NumPyPrinter(), hypot(a, b))), "numpy.hypot(a, b)")
})

test_that("printing_test_numpy: test_log10", {
  x <- Symbol("x"); i <- Symbol("i"); j <- Symbol("j"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d")
  U_numpy_known_constants <- .cg_prefixed(.cg_kc_numpy, "numpy.")
  U_numpy_known_functions <- .cg_prefixed(.cg_kf_numpy, "numpy.")
  U_scipy_known_constants <- .cg_prefixed(.cg_kc_scipy, "scipy.constants.")
  U_scipy_known_functions <- .cg_prefixed(.cg_kf_scipy, "scipy.special.")
  mi_sorted <- function(p) {
    mi <- unclass(ast_attr(p, "module_imports"))
    pystr(paste0("{", paste(vapply(names(mi), function(k) sprintf("'%s': [%s]", k, paste0("'", sort(mi[[k]], method = "radix"), "'", collapse = ", ")), ""), collapse = ", "), "}"))
  }
  expect_identical(.cgb(printer_doprint(NumPyPrinter(), sym_log10(a))), "numpy.log10(a)")
})

test_that("printing_test_numpy: test_exp2", {
  x <- Symbol("x"); i <- Symbol("i"); j <- Symbol("j"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d")
  U_numpy_known_constants <- .cg_prefixed(.cg_kc_numpy, "numpy.")
  U_numpy_known_functions <- .cg_prefixed(.cg_kf_numpy, "numpy.")
  U_scipy_known_constants <- .cg_prefixed(.cg_kc_scipy, "scipy.constants.")
  U_scipy_known_functions <- .cg_prefixed(.cg_kf_scipy, "scipy.special.")
  mi_sorted <- function(p) {
    mi <- unclass(ast_attr(p, "module_imports"))
    pystr(paste0("{", paste(vapply(names(mi), function(k) sprintf("'%s': [%s]", k, paste0("'", sort(mi[[k]], method = "radix"), "'", collapse = ", ")), ""), collapse = ", "), "}"))
  }
  expect_identical(.cgb(printer_doprint(NumPyPrinter(), exp2(a))), "numpy.exp2(a)")
})

test_that("printing_test_numpy: test_log2", {
  x <- Symbol("x"); i <- Symbol("i"); j <- Symbol("j"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d")
  U_numpy_known_constants <- .cg_prefixed(.cg_kc_numpy, "numpy.")
  U_numpy_known_functions <- .cg_prefixed(.cg_kf_numpy, "numpy.")
  U_scipy_known_constants <- .cg_prefixed(.cg_kc_scipy, "scipy.constants.")
  U_scipy_known_functions <- .cg_prefixed(.cg_kf_scipy, "scipy.special.")
  mi_sorted <- function(p) {
    mi <- unclass(ast_attr(p, "module_imports"))
    pystr(paste0("{", paste(vapply(names(mi), function(k) sprintf("'%s': [%s]", k, paste0("'", sort(mi[[k]], method = "radix"), "'", collapse = ", ")), ""), collapse = ", "), "}"))
  }
  expect_identical(.cgb(printer_doprint(NumPyPrinter(), sym_log2(a))), "numpy.log2(a)")
})

test_that("printing_test_numpy: test_Sqrt", {
  x <- Symbol("x"); i <- Symbol("i"); j <- Symbol("j"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d")
  U_numpy_known_constants <- .cg_prefixed(.cg_kc_numpy, "numpy.")
  U_numpy_known_functions <- .cg_prefixed(.cg_kf_numpy, "numpy.")
  U_scipy_known_constants <- .cg_prefixed(.cg_kc_scipy, "scipy.constants.")
  U_scipy_known_functions <- .cg_prefixed(.cg_kf_scipy, "scipy.special.")
  mi_sorted <- function(p) {
    mi <- unclass(ast_attr(p, "module_imports"))
    pystr(paste0("{", paste(vapply(names(mi), function(k) sprintf("'%s': [%s]", k, paste0("'", sort(mi[[k]], method = "radix"), "'", collapse = ", ")), ""), collapse = ", "), "}"))
  }
  expect_identical(.cgb(printer_doprint(NumPyPrinter(), Sqrt(a))), "numpy.sqrt(a)")
})

test_that("printing_test_numpy: test_sqrt", {
  x <- Symbol("x"); i <- Symbol("i"); j <- Symbol("j"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d")
  U_numpy_known_constants <- .cg_prefixed(.cg_kc_numpy, "numpy.")
  U_numpy_known_functions <- .cg_prefixed(.cg_kf_numpy, "numpy.")
  U_scipy_known_constants <- .cg_prefixed(.cg_kc_scipy, "scipy.constants.")
  U_scipy_known_functions <- .cg_prefixed(.cg_kf_scipy, "scipy.special.")
  mi_sorted <- function(p) {
    mi <- unclass(ast_attr(p, "module_imports"))
    pystr(paste0("{", paste(vapply(names(mi), function(k) sprintf("'%s': [%s]", k, paste0("'", sort(mi[[k]], method = "radix"), "'", collapse = ", ")), ""), collapse = ", "), "}"))
  }
  expect_identical(.cgb(printer_doprint(NumPyPrinter(), sqrt(a))), "numpy.sqrt(a)")
})

test_that("printing_test_numpy: test_matsolve", {
  x <- Symbol("x"); i <- Symbol("i"); j <- Symbol("j"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d")
  U_numpy_known_constants <- .cg_prefixed(.cg_kc_numpy, "numpy.")
  U_numpy_known_functions <- .cg_prefixed(.cg_kf_numpy, "numpy.")
  U_scipy_known_constants <- .cg_prefixed(.cg_kc_scipy, "scipy.constants.")
  U_scipy_known_functions <- .cg_prefixed(.cg_kf_scipy, "scipy.special.")
  mi_sorted <- function(p) {
    mi <- unclass(ast_attr(p, "module_imports"))
    pystr(paste0("{", paste(vapply(names(mi), function(k) sprintf("'%s': [%s]", k, paste0("'", sort(mi[[k]], method = "radix"), "'", collapse = ", ")), ""), collapse = ", "), "}"))
  }
  suppressWarnings(try({ M <- MatrixSymbol("M", 3L, 3L) }, silent = TRUE))
  suppressWarnings(try({ x <- MatrixSymbol("x", 3L, 1L) }, silent = TRUE))
  expect_identical(.cgb(printer_doprint(NumPyPrinter(), ((M ^ (-1L)) * x) + x)), "(numpy.linalg.inv(M)).dot(x) + x")
  expect_identical(.cgb(printer_doprint(NumPyPrinter(), MatrixSolve(M, x) + x)), "numpy.linalg.solve(M, x) + x")
})

test_that("printing_test_numpy: test_16857", {
  x <- Symbol("x"); i <- Symbol("i"); j <- Symbol("j"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d")
  U_numpy_known_constants <- .cg_prefixed(.cg_kc_numpy, "numpy.")
  U_numpy_known_functions <- .cg_prefixed(.cg_kf_numpy, "numpy.")
  U_scipy_known_constants <- .cg_prefixed(.cg_kc_scipy, "scipy.constants.")
  U_scipy_known_functions <- .cg_prefixed(.cg_kf_scipy, "scipy.special.")
  mi_sorted <- function(p) {
    mi <- unclass(ast_attr(p, "module_imports"))
    pystr(paste0("{", paste(vapply(names(mi), function(k) sprintf("'%s': [%s]", k, paste0("'", sort(mi[[k]], method = "radix"), "'", collapse = ", ")), ""), collapse = ", "), "}"))
  }
  suppressWarnings(try({ a_1 <- MatrixSymbol("a_1", 10L, 3L) }, silent = TRUE))
  suppressWarnings(try({ a_2 <- MatrixSymbol("a_2", 10L, 3L) }, silent = TRUE))
  suppressWarnings(try({ a_3 <- MatrixSymbol("a_3", 10L, 3L) }, silent = TRUE))
  suppressWarnings(try({ a_4 <- MatrixSymbol("a_4", 10L, 3L) }, silent = TRUE))
  suppressWarnings(try({ A <- BlockMatrix(list(list(a_1, a_2), list(a_3, a_4))) }, silent = TRUE))
  expect_identical(.cgb(structure(.mshape(A), class = "pytuple")), "(20, 6)")
  suppressWarnings(try({ printer <- NumPyPrinter() }, silent = TRUE))
  expect_identical(.cgb(printer_doprint(printer, A)), "numpy.block([[a_1, a_2], [a_3, a_4]])")
})

test_that("printing_test_numpy: test_issue_17006", {
  x <- Symbol("x"); i <- Symbol("i"); j <- Symbol("j"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d")
  U_numpy_known_constants <- .cg_prefixed(.cg_kc_numpy, "numpy.")
  U_numpy_known_functions <- .cg_prefixed(.cg_kf_numpy, "numpy.")
  U_scipy_known_constants <- .cg_prefixed(.cg_kc_scipy, "scipy.constants.")
  U_scipy_known_functions <- .cg_prefixed(.cg_kf_scipy, "scipy.special.")
  mi_sorted <- function(p) {
    mi <- unclass(ast_attr(p, "module_imports"))
    pystr(paste0("{", paste(vapply(names(mi), function(k) sprintf("'%s': [%s]", k, paste0("'", sort(mi[[k]], method = "radix"), "'", collapse = ", ")), ""), collapse = ", "), "}"))
  }
  suppressWarnings(try({ M <- MatrixSymbol("M", 2L, 2L) }, silent = TRUE))
  suppressWarnings(try({ n <- sym_symbols("n", integer = TRUE) }, silent = TRUE))
  suppressWarnings(try({ N <- MatrixSymbol("M", n, n) }, silent = TRUE))
  expect_identical(.cgb(printer_doprint(NumPyPrinter(), N + Identity(n))), "ERR NotImplementedError: Symbolic matrix dimensions are not yet supported for identity matrices")
})

test_that("printing_test_numpy: test_jax_tuple_compatibility", {
  x <- Symbol("x"); i <- Symbol("i"); j <- Symbol("j"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d")
  U_numpy_known_constants <- .cg_prefixed(.cg_kc_numpy, "numpy.")
  U_numpy_known_functions <- .cg_prefixed(.cg_kf_numpy, "numpy.")
  U_scipy_known_constants <- .cg_prefixed(.cg_kc_scipy, "scipy.constants.")
  U_scipy_known_functions <- .cg_prefixed(.cg_kf_scipy, "scipy.special.")
  mi_sorted <- function(p) {
    mi <- unclass(ast_attr(p, "module_imports"))
    pystr(paste0("{", paste(vapply(names(mi), function(k) sprintf("'%s': [%s]", k, paste0("'", sort(mi[[k]], method = "radix"), "'", collapse = ", ")), ""), collapse = ", "), "}"))
  }
  suppressWarnings(try({ .tmp <- sym_symbols("x y z"); x <- py_getitem(.tmp, 0L); y <- py_getitem(.tmp, 1L); z <- py_getitem(.tmp, 2L) }, silent = TRUE))
  suppressWarnings(try({ expr <- Max(x, y, z) + Min(x, y, z) }, silent = TRUE))
  expect_identical(.cgb(printer_doprint(JaxPrinter(), expr)), "functools.reduce(jax.numpy.maximum, [x, y, z]) + functools.reduce(jax.numpy.minimum, [x, y, z])")
  expect_identical(.cgb(printer_doprint(NumPyPrinter(), expr)), "functools.reduce(numpy.maximum, [x, y, z]) + functools.reduce(numpy.minimum, [x, y, z])")
})

test_that("printing_test_numpy: test_numpy_array", {
  x <- Symbol("x"); i <- Symbol("i"); j <- Symbol("j"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d")
  U_numpy_known_constants <- .cg_prefixed(.cg_kc_numpy, "numpy.")
  U_numpy_known_functions <- .cg_prefixed(.cg_kf_numpy, "numpy.")
  U_scipy_known_constants <- .cg_prefixed(.cg_kc_scipy, "scipy.constants.")
  U_scipy_known_functions <- .cg_prefixed(.cg_kf_scipy, "scipy.special.")
  mi_sorted <- function(p) {
    mi <- unclass(ast_attr(p, "module_imports"))
    pystr(paste0("{", paste(vapply(names(mi), function(k) sprintf("'%s': [%s]", k, paste0("'", sort(mi[[k]], method = "radix"), "'", collapse = ", ")), ""), collapse = ", "), "}"))
  }
  suppressWarnings(try({ p <- NumPyPrinter() }, silent = TRUE))
  expect_identical(.cgb(printer_doprint(p, Array(list(list(1L, 2L), list(3L, 5L))))), "numpy.array([[1, 2], [3, 5]])")
  expect_identical(.cgb(printer_doprint(p, Array(list(1L, 2L)))), "numpy.array([1, 2])")
  expect_identical(.cgb(printer_doprint(p, Array(list(list(list(1L, 2L, 3L)))))), "numpy.array([[[1, 2, 3]]])")
  expect_identical(.cgb(printer_doprint(p, Array(list(), list(0L)))), "numpy.zeros((0,))")
  expect_identical(.cgb(printer_doprint(p, Array(list(), list(0L, 0L)))), "numpy.zeros((0, 0))")
  expect_identical(.cgb(printer_doprint(p, Array(list(), list(0L, 1L)))), "numpy.zeros((0, 1))")
  expect_identical(.cgb(printer_doprint(p, Array(list(), list(1L, 0L)))), "numpy.zeros((1, 0))")
  expect_identical(.cgb(printer_doprint(p, Array(list(1L), list()))), "numpy.array(1)")
})

test_that("printing_test_numpy: test_numpy_matrix", {
  x <- Symbol("x"); i <- Symbol("i"); j <- Symbol("j"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d")
  U_numpy_known_constants <- .cg_prefixed(.cg_kc_numpy, "numpy.")
  U_numpy_known_functions <- .cg_prefixed(.cg_kf_numpy, "numpy.")
  U_scipy_known_constants <- .cg_prefixed(.cg_kc_scipy, "scipy.constants.")
  U_scipy_known_functions <- .cg_prefixed(.cg_kf_scipy, "scipy.special.")
  mi_sorted <- function(p) {
    mi <- unclass(ast_attr(p, "module_imports"))
    pystr(paste0("{", paste(vapply(names(mi), function(k) sprintf("'%s': [%s]", k, paste0("'", sort(mi[[k]], method = "radix"), "'", collapse = ", ")), ""), collapse = ", "), "}"))
  }
  suppressWarnings(try({ p <- NumPyPrinter() }, silent = TRUE))
  expect_identical(.cgb(printer_doprint(p, Matrix(list(list(1L, 2L), list(3L, 5L))))), "numpy.array([[1, 2], [3, 5]])")
  expect_identical(.cgb(printer_doprint(p, Matrix(list(1L, 2L)))), "numpy.array([[1], [2]])")
  expect_identical(.cgb(printer_doprint(p, Matrix(0L, 0L, list()))), "numpy.zeros((0, 0))")
  expect_identical(.cgb(printer_doprint(p, Matrix(0L, 1L, list()))), "numpy.zeros((0, 1))")
  expect_identical(.cgb(printer_doprint(p, Matrix(1L, 0L, list()))), "numpy.zeros((1, 0))")
})

test_that("printing_test_numpy: test_numpy_known_funcs_consts", {
  x <- Symbol("x"); i <- Symbol("i"); j <- Symbol("j"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d")
  U_numpy_known_constants <- .cg_prefixed(.cg_kc_numpy, "numpy.")
  U_numpy_known_functions <- .cg_prefixed(.cg_kf_numpy, "numpy.")
  U_scipy_known_constants <- .cg_prefixed(.cg_kc_scipy, "scipy.constants.")
  U_scipy_known_functions <- .cg_prefixed(.cg_kf_scipy, "scipy.special.")
  mi_sorted <- function(p) {
    mi <- unclass(ast_attr(p, "module_imports"))
    pystr(paste0("{", paste(vapply(names(mi), function(k) sprintf("'%s': [%s]", k, paste0("'", sort(mi[[k]], method = "radix"), "'", collapse = ", ")), ""), collapse = ", "), "}"))
  }
  expect_identical(.cgb(py_getitem(U_numpy_known_constants, "NaN")), "numpy.nan")
  expect_identical(.cgb(py_getitem(U_numpy_known_constants, "EulerGamma")), "numpy.euler_gamma")
  expect_identical(.cgb(py_getitem(U_numpy_known_functions, "acos")), "numpy.arccos")
  expect_identical(.cgb(py_getitem(U_numpy_known_functions, "log")), "numpy.log")
})

test_that("printing_test_numpy: test_scipy_known_funcs_consts", {
  x <- Symbol("x"); i <- Symbol("i"); j <- Symbol("j"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d")
  U_numpy_known_constants <- .cg_prefixed(.cg_kc_numpy, "numpy.")
  U_numpy_known_functions <- .cg_prefixed(.cg_kf_numpy, "numpy.")
  U_scipy_known_constants <- .cg_prefixed(.cg_kc_scipy, "scipy.constants.")
  U_scipy_known_functions <- .cg_prefixed(.cg_kf_scipy, "scipy.special.")
  mi_sorted <- function(p) {
    mi <- unclass(ast_attr(p, "module_imports"))
    pystr(paste0("{", paste(vapply(names(mi), function(k) sprintf("'%s': [%s]", k, paste0("'", sort(mi[[k]], method = "radix"), "'", collapse = ", ")), ""), collapse = ", "), "}"))
  }
  expect_identical(.cgb(py_getitem(U_scipy_known_constants, "GoldenRatio")), "scipy.constants.golden_ratio")
  expect_identical(.cgb(py_getitem(U_scipy_known_constants, "Pi")), "scipy.constants.pi")
  expect_identical(.cgb(py_getitem(U_scipy_known_functions, "erf")), "scipy.special.erf")
  expect_identical(.cgb(py_getitem(U_scipy_known_functions, "factorial")), "scipy.special.factorial")
})

test_that("printing_test_numpy: test_numpy_print_methods", {
  x <- Symbol("x"); i <- Symbol("i"); j <- Symbol("j"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d")
  U_numpy_known_constants <- .cg_prefixed(.cg_kc_numpy, "numpy.")
  U_numpy_known_functions <- .cg_prefixed(.cg_kf_numpy, "numpy.")
  U_scipy_known_constants <- .cg_prefixed(.cg_kc_scipy, "scipy.constants.")
  U_scipy_known_functions <- .cg_prefixed(.cg_kf_scipy, "scipy.special.")
  mi_sorted <- function(p) {
    mi <- unclass(ast_attr(p, "module_imports"))
    pystr(paste0("{", paste(vapply(names(mi), function(k) sprintf("'%s': [%s]", k, paste0("'", sort(mi[[k]], method = "radix"), "'", collapse = ", ")), ""), collapse = ", "), "}"))
  }
  suppressWarnings(try({ prntr <- NumPyPrinter() }, silent = TRUE))
  expect_identical(.cgb(hasattr(prntr, "_print_acos")), "True")
  expect_identical(.cgb(hasattr(prntr, "_print_log")), "True")
})

test_that("printing_test_numpy: test_scipy_print_methods", {
  x <- Symbol("x"); i <- Symbol("i"); j <- Symbol("j"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d")
  U_numpy_known_constants <- .cg_prefixed(.cg_kc_numpy, "numpy.")
  U_numpy_known_functions <- .cg_prefixed(.cg_kf_numpy, "numpy.")
  U_scipy_known_constants <- .cg_prefixed(.cg_kc_scipy, "scipy.constants.")
  U_scipy_known_functions <- .cg_prefixed(.cg_kf_scipy, "scipy.special.")
  mi_sorted <- function(p) {
    mi <- unclass(ast_attr(p, "module_imports"))
    pystr(paste0("{", paste(vapply(names(mi), function(k) sprintf("'%s': [%s]", k, paste0("'", sort(mi[[k]], method = "radix"), "'", collapse = ", ")), ""), collapse = ", "), "}"))
  }
  suppressWarnings(try({ prntr <- SciPyPrinter() }, silent = TRUE))
  expect_identical(.cgb(hasattr(prntr, "_print_acos")), "True")
  expect_identical(.cgb(hasattr(prntr, "_print_log")), "True")
  expect_identical(.cgb(hasattr(prntr, "_print_erf")), "True")
  expect_identical(.cgb(hasattr(prntr, "_print_factorial")), "True")
  expect_identical(.cgb(hasattr(prntr, "_print_chebyshevt")), "True")
  suppressWarnings(try({ k <- Symbol("k", integer = TRUE, nonnegative = TRUE) }, silent = TRUE))
  suppressWarnings(try({ x <- Symbol("x", real = TRUE) }, silent = TRUE))
  expect_identical(.cgb(printer_doprint(prntr, polygamma(k, x))), "scipy.special.polygamma(k, x)")
  expect_identical(.cgb(printer_doprint(prntr, Si(x))), "scipy.special.sici(x)[0]")
  expect_identical(.cgb(printer_doprint(prntr, Ci(x))), "scipy.special.sici(x)[1]")
  expect_identical(.cgb(printer_doprint(prntr, owens_t(x, a))), "scipy.special.owens_t(x, a)")
})

test_that("printing_test_numpy: test_numpy_contains_integers", {
  x <- Symbol("x"); i <- Symbol("i"); j <- Symbol("j"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d")
  U_numpy_known_constants <- .cg_prefixed(.cg_kc_numpy, "numpy.")
  U_numpy_known_functions <- .cg_prefixed(.cg_kf_numpy, "numpy.")
  U_scipy_known_constants <- .cg_prefixed(.cg_kc_scipy, "scipy.constants.")
  U_scipy_known_functions <- .cg_prefixed(.cg_kf_scipy, "scipy.special.")
  mi_sorted <- function(p) {
    mi <- unclass(ast_attr(p, "module_imports"))
    pystr(paste0("{", paste(vapply(names(mi), function(k) sprintf("'%s': [%s]", k, paste0("'", sort(mi[[k]], method = "radix"), "'", collapse = ", ")), ""), collapse = ", "), "}"))
  }
  suppressWarnings(try({ x <- sym_symbols("x") }, silent = TRUE))
  expect_identical(.cgb(printer_doprint(NumPyPrinter(), Contains(x, Integers))), "equal(mod(x, 1), 0)")
})

test_that("printing_test_numpy: test_issue_30395_scalar_assignment", {
  x <- Symbol("x"); i <- Symbol("i"); j <- Symbol("j"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d")
  U_numpy_known_constants <- .cg_prefixed(.cg_kc_numpy, "numpy.")
  U_numpy_known_functions <- .cg_prefixed(.cg_kf_numpy, "numpy.")
  U_scipy_known_constants <- .cg_prefixed(.cg_kc_scipy, "scipy.constants.")
  U_scipy_known_functions <- .cg_prefixed(.cg_kf_scipy, "scipy.special.")
  mi_sorted <- function(p) {
    mi <- unclass(ast_attr(p, "module_imports"))
    pystr(paste0("{", paste(vapply(names(mi), function(k) sprintf("'%s': [%s]", k, paste0("'", sort(mi[[k]], method = "radix"), "'", collapse = ", ")), ""), collapse = ", "), "}"))
  }
  suppressWarnings(try({ .tmp <- sym_symbols("x y z"); x <- py_getitem(.tmp, 0L); y <- py_getitem(.tmp, 1L); z <- py_getitem(.tmp, 2L) }, silent = TRUE))
  suppressWarnings(try({ np_p <- NumPyPrinter() }, silent = TRUE))
  suppressWarnings(try({ jp <- JaxPrinter() }, silent = TRUE))
  expect_identical(.cgb(printer_doprint(np_p, (x + y) + z, "out")), "out = x + y + z")
  expect_identical(.cgb(printer_doprint(jp, (x + y) + z, "out")), "out = x + y + z")
  suppressWarnings(try({ expr <- ((x ^ 2L) + sqrt(y)) - sin(x * z) }, silent = TRUE))
  expect_identical(.cgb(printer_doprint(np_p, expr, "out")), "out = x^2 + numpy.sqrt(y) - numpy.sin(x*z)")
  expect_identical(.cgb(printer_doprint(jp, expr, "out")), "out = x^2 + jax.numpy.sqrt(y) - jax.numpy.sin(x*z)")
  expect_identical(.cgb(printer_doprint(np_p, x * z, "out")), "out = x*z")
  suppressWarnings(try({ A <- ArraySymbol("A", list(3L, 3L)) }, silent = TRUE))
  expect_identical(.cgb(printer_doprint(np_p, A + A, "out")), "out = numpy.einsum(\",ab\", 2, A)")
  expect_identical(.cgb(printer_doprint(jp, ArrayElementwiseApplyFunc(sin, A))), "jax.numpy.sin(A)")
  suppressWarnings(try({ d <- Dummy("d") }, silent = TRUE))
  suppressWarnings(try({ M <- MatrixSymbol("M", 3L, 3L) }, silent = TRUE))
  expect_identical(.cgb(printer_doprint(np_p, ArrayElementwiseApplyFunc(Lambda(d, d ^ 2L), M))), "M^2")
  suppressWarnings(try({ N <- MatrixSymbol("N", 3L, 3L) }, silent = TRUE))
  expect_identical(.cgb(printer_doprint(np_p, ArrayElementwiseApplyFunc(Lambda(d, (2L * d) + 1L), M + N))), "2*(M + N) + 1")
})

test_that("printing_test_numpy: cupy_jax_printers", {
  x <- Symbol("x"); i <- Symbol("i"); j <- Symbol("j"); a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d")
  U_numpy_known_constants <- .cg_prefixed(.cg_kc_numpy, "numpy.")
  U_numpy_known_functions <- .cg_prefixed(.cg_kf_numpy, "numpy.")
  U_scipy_known_constants <- .cg_prefixed(.cg_kc_scipy, "scipy.constants.")
  U_scipy_known_functions <- .cg_prefixed(.cg_kf_scipy, "scipy.special.")
  mi_sorted <- function(p) {
    mi <- unclass(ast_attr(p, "module_imports"))
    pystr(paste0("{", paste(vapply(names(mi), function(k) sprintf("'%s': [%s]", k, paste0("'", sort(mi[[k]], method = "radix"), "'", collapse = ", ")), ""), collapse = ", "), "}"))
  }
  suppressWarnings(try({ .tmp <- sym_symbols("x y z"); x <- py_getitem(.tmp, 0L); y <- py_getitem(.tmp, 1L); z <- py_getitem(.tmp, 2L) }, silent = TRUE))
  suppressWarnings(try({ cp <- CuPyPrinter() }, silent = TRUE))
  suppressWarnings(try({ jp <- JaxPrinter() }, silent = TRUE))
  expect_identical(.cgb(printer_doprint(cp, Piecewise(list(1L, x < 0L), list(0L, TRUE)))), "cupy.select([cupy.less(x, 0),True], [1,0], default=cupy.nan)")
  expect_identical(.cgb(mi_sorted(cp)), "{'cupy': ['less', 'nan', 'select']}")
  expect_identical(.cgb(printer_doprint(cp, logaddexp(a, b))), "cupy.logaddexp(a, b)")
  expect_identical(.cgb(printer_doprint(cp, conjugate(x))), "cupy.conjugate(x)")
  expect_identical(.cgb(printer_doprint(cp, sqrt(a) + sin(a))), "cupy.sqrt(a) + cupy.sin(a)")
  expect_identical(.cgb(printer_doprint(cp, Max(x, y, z))), "functools.reduce(cupy.maximum, [x, y, z])")
  expect_identical(.cgb(printer_doprint(cp, And(x < 1L, y > 2L))), "cupy.logical_and.reduce((cupy.greater(y, 2),cupy.less(x, 1)))")
  expect_identical(.cgb(printer_doprint(cp, Or(x < 1L, y > 2L))), "cupy.logical_or.reduce((cupy.greater(y, 2),cupy.less(x, 1)))")
  expect_identical(.cgb(printer_doprint(jp, And(x < 1L, y > 2L))), "jax.numpy.all(jax.numpy.asarray([jax.numpy.greater(y, 2),jax.numpy.less(x, 1)]), axis=0)")
  expect_identical(.cgb(printer_doprint(jp, Or(x < 1L, y > 2L))), "jax.numpy.any(jax.numpy.asarray([jax.numpy.greater(y, 2),jax.numpy.less(x, 1)]), axis=0)")
  expect_identical(.cgb(printer_doprint(jp, Not(x < 1L))), "jax.numpy.greater_equal(x, 1)")
  expect_identical(.cgb(printer_doprint(jp, Piecewise(list(1L, x < 0L), list(0L, TRUE)))), "jax.numpy.select([jax.numpy.less(x, 0),True], [1,0], default=jax.numpy.nan)")
  expect_identical(.cgb(mi_sorted(jp)), "{'jax.numpy': ['all', 'any', 'greater', 'greater_equal', 'less', 'nan', 'select'], 'jax': ['numpy']}")
  expect_identical(.cgb(printer_doprint(cp, Matrix(list(list(1L, 2L), list(3L, 5L))))), "cupy.array([[1, 2], [3, 5]])")
  expect_identical(.cgb(printer_doprint(jp, Matrix(list(list(1L, 2L), list(3L, 5L))))), "jax.numpy.array([[1, 2], [3, 5]])")
  expect_identical(.cgb(printer_doprint(jp, (acos(x) + exp(x)) + S_Pi)), "jax.numpy.exp(x) + jax.numpy.arccos(x) + jax.numpy.pi")
  expect_identical(.cgb(printer_doprint(cp, ((acos(x) + exp(x)) + S_Pi) + EulerGamma)), "cupy.exp(x) + cupy.arccos(x) + cupy.euler_gamma + cupy.pi")
  expect_identical(.cgb(ast_attr(cp, "language")), "Python with cupy")
  expect_identical(.cgb(ast_attr(jp, "language")), "Python with jax.numpy")
  expect_identical(.cgb(ast_attr(cp, "printmethod")), "_cupycode")
  expect_identical(.cgb(ast_attr(jp, "printmethod")), "_jaxcode")
  expect_identical(.cgb(printer_doprint(cp, MatrixSymbol("M", 2L, 2L) ^ 3L)), "cupy.linalg.matrix_power(M, 3)")
  expect_identical(.cgb(printer_doprint(jp, Inverse(MatrixSymbol("M", 2L, 2L)))), "jax.numpy.linalg.inv(M)")
  expect_identical(.cgb(printer_doprint(cp, ArrayTensorProduct(MatrixSymbol("M", 2L, 2L), MatrixSymbol("N", 2L, 2L)))), "cupy.einsum(\"ab,cd\", M, N)")
  expect_identical(.cgb(printer_doprint(jp, ArrayAdd(MatrixSymbol("M", 2L, 2L), MatrixSymbol("N", 2L, 2L), MatrixSymbol("P", 2L, 2L)))), "jax.numpy.add(jax.numpy.add(M, N), P)")
  expect_identical(.cgb(printer_doprint(jp, PermuteDims(MatrixSymbol("M", 2L, 2L), list(1L, 0L)))), "jax.numpy.transpose(M, [1, 0])")
  expect_identical(.cgb(printer_doprint(cp, ArrayDiagonal(ArrayTensorProduct(MatrixSymbol("M", 2L, 2L), MatrixSymbol("N", 2L, 2L)), list(1L, 2L)))), "cupy.einsum(\"ab,bc->acb\", M, N)")
  expect_identical(.cgb(printer_doprint(jp, Array(list(), list(0L, 2L)))), "jax.numpy.zeros((0, 2))")
  expect_identical(.cgb(printer_doprint(cp, x ^ (-2L))), "x^(-2.0)")
  expect_identical(.cgb(printer_doprint(jp, sym_Mod(x, y))), "jax.numpy.mod(x, y)")
  expect_identical(.cgb(printer_doprint(jp, re(x) + im(x))), "jax.numpy.real(x) + jax.numpy.imag(x)")
})

test_that("printing_test_tensorflow: test_tensorflow_printing", {
  M <- MatrixSymbol("M", 3L, 3L); N <- MatrixSymbol("N", 3L, 3L); P <- MatrixSymbol("P", 3L, 3L); Q <- MatrixSymbol("Q", 3L, 3L)
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); t <- Symbol("t")
  expect_identical(.cgb(tensorflow_code(eye(3L))), "tensorflow.constant([[1, 0, 0], [0, 1, 0], [0, 0, 1]])")
  suppressWarnings(try({ expr <- Matrix(list(list(x, sin(y)), list(exp(z), -t))) }, silent = TRUE))
  expect_identical(.cgb(tensorflow_code(expr)), "tensorflow.Variable([[x, tensorflow.math.sin(y)], [tensorflow.math.exp(z), -t]])")
})

test_that("printing_test_tensorflow: test_tensorflow_math", {
  M <- MatrixSymbol("M", 3L, 3L); N <- MatrixSymbol("N", 3L, 3L); P <- MatrixSymbol("P", 3L, 3L); Q <- MatrixSymbol("Q", 3L, 3L)
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); t <- Symbol("t")
  expect_identical(.cgb(tensorflow_code(Abs(x))), "tensorflow.math.abs(x)")
  expect_identical(.cgb(tensorflow_code(sign(x))), "tensorflow.math.sign(x)")
  expect_identical(.cgb(tensorflow_code(ceiling(x))), "tensorflow.math.ceil(x)")
  expect_identical(.cgb(tensorflow_code(floor(x))), "tensorflow.math.floor(x)")
  expect_identical(.cgb(tensorflow_code(exp(x))), "tensorflow.math.exp(x)")
  expect_identical(.cgb(tensorflow_code(sqrt(x))), "tensorflow.math.sqrt(x)")
  expect_identical(.cgb(tensorflow_code(x ^ 4L)), "tensorflow.math.pow(x, 4)")
  expect_identical(.cgb(tensorflow_code(cos(x))), "tensorflow.math.cos(x)")
  expect_identical(.cgb(tensorflow_code(acos(x))), "tensorflow.math.acos(x)")
  expect_identical(.cgb(tensorflow_code(sin(x))), "tensorflow.math.sin(x)")
  expect_identical(.cgb(tensorflow_code(asin(x))), "tensorflow.math.asin(x)")
  expect_identical(.cgb(tensorflow_code(tan(x))), "tensorflow.math.tan(x)")
  expect_identical(.cgb(tensorflow_code(atan(x))), "tensorflow.math.atan(x)")
  expect_identical(.cgb(tensorflow_code(sym_atan2(y, x))), "tensorflow.math.atan2(y, x)")
  expect_identical(.cgb(tensorflow_code(cosh(x))), "tensorflow.math.cosh(x)")
  expect_identical(.cgb(tensorflow_code(acosh(x))), "tensorflow.math.acosh(x)")
  expect_identical(.cgb(tensorflow_code(sinh(x))), "tensorflow.math.sinh(x)")
  expect_identical(.cgb(tensorflow_code(asinh(x))), "tensorflow.math.asinh(x)")
  expect_identical(.cgb(tensorflow_code(tanh(x))), "tensorflow.math.tanh(x)")
  expect_identical(.cgb(tensorflow_code(atanh(x))), "tensorflow.math.atanh(x)")
  expect_identical(.cgb(tensorflow_code(erf(x))), "tensorflow.math.erf(x)")
  expect_identical(.cgb(tensorflow_code(loggamma(x))), "tensorflow.math.lgamma(x)")
})

test_that("printing_test_tensorflow: test_tensorflow_complexes", {
  M <- MatrixSymbol("M", 3L, 3L); N <- MatrixSymbol("N", 3L, 3L); P <- MatrixSymbol("P", 3L, 3L); Q <- MatrixSymbol("Q", 3L, 3L)
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); t <- Symbol("t")
  expect_identical(.cgb(tensorflow_code(re(x))), "tensorflow.math.real(x)")
  expect_identical(.cgb(tensorflow_code(im(x))), "tensorflow.math.imag(x)")
  expect_identical(.cgb(tensorflow_code(arg(x))), "tensorflow.math.angle(x)")
})

test_that("printing_test_tensorflow: test_tensorflow_relational", {
  M <- MatrixSymbol("M", 3L, 3L); N <- MatrixSymbol("N", 3L, 3L); P <- MatrixSymbol("P", 3L, 3L); Q <- MatrixSymbol("Q", 3L, 3L)
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); t <- Symbol("t")
  expect_identical(.cgb(tensorflow_code(Eq(x, y))), "tensorflow.math.equal(x, y)")
  expect_identical(.cgb(tensorflow_code(Ne(x, y))), "tensorflow.math.not_equal(x, y)")
  expect_identical(.cgb(tensorflow_code(Ge(x, y))), "tensorflow.math.greater_equal(x, y)")
  expect_identical(.cgb(tensorflow_code(Gt(x, y))), "tensorflow.math.greater(x, y)")
  expect_identical(.cgb(tensorflow_code(Le(x, y))), "tensorflow.math.less_equal(x, y)")
  expect_identical(.cgb(tensorflow_code(Lt(x, y))), "tensorflow.math.less(x, y)")
})

test_that("printing_test_tensorflow: test_tensorflow_matrices", {
  M <- MatrixSymbol("M", 3L, 3L); N <- MatrixSymbol("N", 3L, 3L); P <- MatrixSymbol("P", 3L, 3L); Q <- MatrixSymbol("Q", 3L, 3L)
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); t <- Symbol("t")
  expect_identical(.cgb(tensorflow_code(M)), "M")
  expect_identical(.cgb(tensorflow_code(M + N)), "tensorflow.math.add(M, N)")
  expect_identical(.cgb(tensorflow_code(M * N)), "tensorflow.linalg.matmul(M, N)")
  expect_identical(.cgb(tensorflow_code(HadamardProduct(M, N))), "tensorflow.math.multiply(M, N)")
  expect_identical(.cgb(tensorflow_code(((M * N) * P) * Q)), "tensorflow.linalg.matmul(tensorflow.linalg.matmul(tensorflow.linalg.matmul(M, N), P), Q)")
  expect_identical(.cgb(tensorflow_code(M ^ 3L)), "tensorflow.linalg.matmul(tensorflow.linalg.matmul(M, M), M)")
  expect_identical(.cgb(tensorflow_code(Trace(M))), "tensorflow.linalg.trace(M)")
  expect_identical(.cgb(tensorflow_code(Determinant(M))), "tensorflow.linalg.det(M)")
  expect_identical(.cgb(tensorflow_code(Inverse(M))), "tensorflow.linalg.inv(M)")
  expect_identical(.cgb(tensorflow_code(Transpose(M), tensorflow_version = "1.14")), "tensorflow.linalg.matrix_transpose(M)")
  expect_identical(.cgb(tensorflow_code(Transpose(M), tensorflow_version = "1.13")), "tensorflow.matrix_transpose(M)")
})

test_that("printing_test_tensorflow: test_codegen_einsum", {
  M <- MatrixSymbol("M", 3L, 3L); N <- MatrixSymbol("N", 3L, 3L); P <- MatrixSymbol("P", 3L, 3L); Q <- MatrixSymbol("Q", 3L, 3L)
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); t <- Symbol("t")
  suppressWarnings(try({ M <- MatrixSymbol("M", 2L, 2L) }, silent = TRUE))
  suppressWarnings(try({ N <- MatrixSymbol("N", 2L, 2L) }, silent = TRUE))
  suppressWarnings(try({ cg <- convert_matrix_to_array(M * N) }, silent = TRUE))
  expect_identical(.cgb(tensorflow_code(cg)), "tensorflow.linalg.einsum(\"ab,bc->ac\", M,N)")
})

test_that("printing_test_tensorflow: test_codegen_extra", {
  M <- MatrixSymbol("M", 3L, 3L); N <- MatrixSymbol("N", 3L, 3L); P <- MatrixSymbol("P", 3L, 3L); Q <- MatrixSymbol("Q", 3L, 3L)
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); t <- Symbol("t")
  suppressWarnings(try({ M <- MatrixSymbol("M", 2L, 2L) }, silent = TRUE))
  suppressWarnings(try({ N <- MatrixSymbol("N", 2L, 2L) }, silent = TRUE))
  suppressWarnings(try({ P <- MatrixSymbol("P", 2L, 2L) }, silent = TRUE))
  suppressWarnings(try({ Q <- MatrixSymbol("Q", 2L, 2L) }, silent = TRUE))
  expect_identical(.cgb(tensorflow_code(ArrayTensorProduct(M, N))), "tensorflow.linalg.einsum(\"ab,cd\", M, N)")
  expect_identical(.cgb(tensorflow_code(ArrayAdd(M, N))), "tensorflow.math.add(M, N)")
  expect_identical(.cgb(tensorflow_code(ArrayAdd(M, N, P))), "tensorflow.math.add(tensorflow.math.add(M, N), P)")
  expect_identical(.cgb(tensorflow_code(ArrayAdd(M, N, P, Q))), "tensorflow.math.add(tensorflow.math.add(tensorflow.math.add(M, N), P), Q)")
  expect_identical(.cgb(tensorflow_code(PermuteDims(M, list(1L, 0L)))), "tensorflow.transpose(M, [1, 0])")
  expect_identical(.cgb(tensorflow_code(PermuteDims(ArrayTensorProduct(M, N), list(1L, 2L, 3L, 0L)))), "tensorflow.transpose(tensorflow.linalg.einsum(\"ab,cd\", M, N), [1, 2, 3, 0])")
  expect_identical(.cgb(tensorflow_code(ArrayDiagonal(ArrayTensorProduct(M, N), list(1L, 2L)))), "tensorflow.linalg.einsum(\"ab,bc->acb\", M, N)")
})

test_that("printing_test_tensorflow: test_MatrixElement_printing", {
  M <- MatrixSymbol("M", 3L, 3L); N <- MatrixSymbol("N", 3L, 3L); P <- MatrixSymbol("P", 3L, 3L); Q <- MatrixSymbol("Q", 3L, 3L)
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); t <- Symbol("t")
  suppressWarnings(try({ A <- MatrixSymbol("A", 1L, 3L) }, silent = TRUE))
  suppressWarnings(try({ B <- MatrixSymbol("B", 1L, 3L) }, silent = TRUE))
  suppressWarnings(try({ C <- MatrixSymbol("C", 1L, 3L) }, silent = TRUE))
  expect_identical(.cgb(tensorflow_code(py_getitem(A, 0L, 0L))), "A[0, 0]")
  expect_identical(.cgb(tensorflow_code(3L * py_getitem(A, 0L, 0L))), "3*A[0, 0]")
  suppressWarnings(try({ F <- new_sym("MatrixElement", list(A - B, S(0L), S(0L))) }, silent = TRUE))
})

test_that("printing_test_tensorflow: test_tensorflow_Derivative", {
  M <- MatrixSymbol("M", 3L, 3L); N <- MatrixSymbol("N", 3L, 3L); P <- MatrixSymbol("P", 3L, 3L); Q <- MatrixSymbol("Q", 3L, 3L)
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); t <- Symbol("t")
  suppressWarnings(try({ expr <- Derivative(sin(x), x) }, silent = TRUE))
  expect_identical(.cgb(tensorflow_code(expr)), "tensorflow.gradients(tensorflow.math.sin(x), x)[0]")
})

test_that("printing_test_tensorflow: test_tensorflow_isnan_isinf", {
  M <- MatrixSymbol("M", 3L, 3L); N <- MatrixSymbol("N", 3L, 3L); P <- MatrixSymbol("P", 3L, 3L); Q <- MatrixSymbol("Q", 3L, 3L)
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); t <- Symbol("t")
  suppressWarnings(try({ x <- sym_symbols("x") }, silent = TRUE))
  expect_identical(.cgb(tensorflow_code(Piecewise(list(Float(0.0), isnan(x)), list(Float(1.0), TRUE)))), "tensorflow.where(tensorflow.math.is_nan(x), 0.0, 1.0)")
  expect_identical(.cgb(tensorflow_code(Piecewise(list(Float(0.0), isinf(x)), list(Float(1.0), TRUE)))), "tensorflow.where(tensorflow.math.is_inf(x), 0.0, 1.0)")
})

test_that("printing_test_tensorflow: test_tensorflow_elementwise_apply_func", {
  M <- MatrixSymbol("M", 3L, 3L); N <- MatrixSymbol("N", 3L, 3L); P <- MatrixSymbol("P", 3L, 3L); Q <- MatrixSymbol("Q", 3L, 3L)
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); t <- Symbol("t")
  suppressWarnings(try({ d <- Dummy("d") }, silent = TRUE))
  expect_identical(.cgb(tensorflow_code(ArrayElementwiseApplyFunc(sin, M))), "tensorflow.math.sin(M)")
  expect_identical(.cgb(tensorflow_code(ArrayElementwiseApplyFunc(Lambda(d, d ^ 2L), M))), "tensorflow.math.pow(M, 2)")
})

test_that("printing_test_tensorflow: tensorflow_extra", {
  M <- MatrixSymbol("M", 3L, 3L); N <- MatrixSymbol("N", 3L, 3L); P <- MatrixSymbol("P", 3L, 3L); Q <- MatrixSymbol("Q", 3L, 3L)
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); t <- Symbol("t")
  expect_identical(.cgb(tensorflow_code(Piecewise(list(x, x < 0L), list(y, x < 1L), list(z, TRUE)))), "tensorflow.where(tensorflow.math.less(x, 0), x, tensorflow.where(tensorflow.math.less(x, 1), y, z))")
  expect_identical(.cgb(tensorflow_code(Piecewise(list(x, x < 0L)))), "tensorflow.where(tensorflow.math.less(x, 0), x, 0)")
  expect_identical(.cgb(tensorflow_code(Piecewise(list(x, x < 0L), list(y, TRUE)), tensorflow_version = "0.9")), "tensorflow.select(tensorflow.math.less(x, 0), x, y)")
  expect_identical(.cgb(tensorflow_code(And(x < 1L, y > 2L))), "tensorflow.math.logical_and(tensorflow.math.greater(y, 2), tensorflow.math.less(x, 1))")
  expect_identical(.cgb(tensorflow_code(Or(x < 1L, y > 2L, z > 3L))), "tensorflow.math.logical_or(tensorflow.math.logical_or(tensorflow.math.greater(y, 2), tensorflow.math.greater(z, 3)), tensorflow.math.less(x, 1))")
  expect_identical(.cgb(tensorflow_code(Not(x < 1L))), "tensorflow.math.greater_equal(x, 1)")
  expect_identical(.cgb(tensorflow_code(Max(x, y, z))), "tensorflow.math.maximum(tensorflow.math.maximum(x, y), z)")
  expect_identical(.cgb(tensorflow_code(Min(x, y))), "tensorflow.math.minimum(x, y)")
  expect_identical(.cgb(tensorflow_code(log(x) + exp(y))), "tensorflow.math.exp(y) + tensorflow.math.log(x)")
  expect_identical(.cgb(tensorflow_code(x ^ y)), "tensorflow.math.pow(x, y)")
  expect_identical(.cgb(tensorflow_code((2L * M) * N)), "2*tensorflow.linalg.matmul(M, N)")
  expect_identical(.cgb(tensorflow_code(Matrix(list(list(1L, 2L), list(3L, 4L))))), "tensorflow.constant([[1, 2], [3, 4]])")
  expect_identical(.cgb(tensorflow_code(Derivative(sin(x), x, y))), "tensorflow.gradients(tensorflow.gradients(tensorflow.math.sin(x), x)[0], y)[0]")
  expect_identical(.cgb(tensorflow_code(Assignment(y, sin(x)))), "y = tensorflow.math.sin(x)")
  expect_identical(.cgb(tensorflow_code((x * y) + z)), "x*y + z")
  expect_identical(.cgb(tensorflow_code(1L / x)), "tensorflow.math.pow(x, -1)")
  expect_identical(.cgb(tensorflow_code(gamma(x))), "ERR PrintMethodNotImplementedError: Unsupported by <class 'sympy.printing.tensorflow.TensorflowPrinter'>: gamma\\nSet the printer option 'strict' to False in order to generate partially printed code.")
  expect_identical(.cgb(tensorflow_code(ZeroArray(3L, 2L))), "tensorflow.zeros((3,2,))")
  expect_identical(.cgb(tensorflow_code(OneArray(4L))), "tensorflow.ones((4,))")
  expect_identical(.cgb(printer_doprint(TensorflowPrinter(), ArrayContraction(ArrayTensorProduct(M, N), list(1L, 2L)))), "tensorflow.linalg.einsum(\"ab,bc->ac\", M,N)")
  expect_identical(.cgb(ast_attr(TensorflowPrinter(), "printmethod")), "_tensorflowcode")
  expect_identical(.cgb(ast_attr(TensorflowPrinter(py_dict(list("tensorflow_version"), list("2.0"))), "tensorflow_version")), "2.0")
})

test_that("printing_test_torch: test_torch_math", {
  M <- MatrixSymbol("M", 3L, 3L); N <- MatrixSymbol("N", 3L, 3L); P <- MatrixSymbol("P", 3L, 3L); Q <- MatrixSymbol("Q", 3L, 3L)
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); t <- Symbol("t")
  expect_identical(.cgb(torch_code(Abs(x))), "torch.abs(x)")
  expect_identical(.cgb(torch_code(sign(x))), "torch.sign(x)")
  expect_identical(.cgb(torch_code(ceiling(x))), "torch.ceil(x)")
  expect_identical(.cgb(torch_code(floor(x))), "torch.floor(x)")
  expect_identical(.cgb(torch_code(exp(x))), "torch.exp(x)")
  expect_identical(.cgb(torch_code(sqrt(x))), "torch.sqrt(x)")
  expect_identical(.cgb(torch_code(x ^ 4L)), "torch.pow(x, 4)")
  expect_identical(.cgb(torch_code(cos(x))), "torch.cos(x)")
  expect_identical(.cgb(torch_code(acos(x))), "torch.acos(x)")
  expect_identical(.cgb(torch_code(sin(x))), "torch.sin(x)")
  expect_identical(.cgb(torch_code(asin(x))), "torch.asin(x)")
  expect_identical(.cgb(torch_code(tan(x))), "torch.tan(x)")
  expect_identical(.cgb(torch_code(atan(x))), "torch.atan(x)")
  expect_identical(.cgb(torch_code(sym_atan2(y, x))), "torch.atan2(y, x)")
  expect_identical(.cgb(torch_code(cosh(x))), "torch.cosh(x)")
  expect_identical(.cgb(torch_code(acosh(x))), "torch.acosh(x)")
  expect_identical(.cgb(torch_code(sinh(x))), "torch.sinh(x)")
  expect_identical(.cgb(torch_code(asinh(x))), "torch.asinh(x)")
  expect_identical(.cgb(torch_code(tanh(x))), "torch.tanh(x)")
  expect_identical(.cgb(torch_code(atanh(x))), "torch.atanh(x)")
  expect_identical(.cgb(torch_code(erf(x))), "torch.erf(x)")
  expect_identical(.cgb(torch_code(loggamma(x))), "torch.lgamma(x)")
})

test_that("printing_test_torch: test_torch_complexes", {
  M <- MatrixSymbol("M", 3L, 3L); N <- MatrixSymbol("N", 3L, 3L); P <- MatrixSymbol("P", 3L, 3L); Q <- MatrixSymbol("Q", 3L, 3L)
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); t <- Symbol("t")
  expect_identical(.cgb(torch_code(re(x))), "torch.real(x)")
  expect_identical(.cgb(torch_code(im(x))), "torch.imag(x)")
  expect_identical(.cgb(torch_code(arg(x))), "torch.angle(x)")
})

test_that("printing_test_torch: test_torch_relational", {
  M <- MatrixSymbol("M", 3L, 3L); N <- MatrixSymbol("N", 3L, 3L); P <- MatrixSymbol("P", 3L, 3L); Q <- MatrixSymbol("Q", 3L, 3L)
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); t <- Symbol("t")
  expect_identical(.cgb(torch_code(Eq(x, y))), "torch.eq(x, y)")
  expect_identical(.cgb(torch_code(Ne(x, y))), "torch.ne(x, y)")
  expect_identical(.cgb(torch_code(Ge(x, y))), "torch.ge(x, y)")
  expect_identical(.cgb(torch_code(Gt(x, y))), "torch.gt(x, y)")
  expect_identical(.cgb(torch_code(Le(x, y))), "torch.le(x, y)")
  expect_identical(.cgb(torch_code(Lt(x, y))), "torch.lt(x, y)")
})

test_that("printing_test_torch: test_torch_max_min", {
  M <- MatrixSymbol("M", 3L, 3L); N <- MatrixSymbol("N", 3L, 3L); P <- MatrixSymbol("P", 3L, 3L); Q <- MatrixSymbol("Q", 3L, 3L)
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); t <- Symbol("t")
  suppressWarnings(try({ .tmp <- sym_symbols("x y"); x <- py_getitem(.tmp, 0L); y <- py_getitem(.tmp, 1L) }, silent = TRUE))
  expect_identical(.cgb(torch_code(Max(x, y) + Min(x, y))), "torch.maximum(x, y) + torch.minimum(x, y)")
})

test_that("printing_test_torch: test_torch_matrix", {
  M <- MatrixSymbol("M", 3L, 3L); N <- MatrixSymbol("N", 3L, 3L); P <- MatrixSymbol("P", 3L, 3L); Q <- MatrixSymbol("Q", 3L, 3L)
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); t <- Symbol("t")
  expect_identical(.cgb(torch_code(M)), "M")
  expect_identical(.cgb(torch_code(M * N)), "torch.matmul(M, N)")
  expect_identical(.cgb(torch_code(M ^ 3L)), "torch.mm(torch.mm(M, M), M)")
  expect_identical(.cgb(torch_code(((M * N) * P) * Q)), "torch.matmul(torch.matmul(torch.matmul(M, N), P), Q)")
  expect_identical(.cgb(torch_code(Trace(M))), "torch.trace(M)")
  expect_identical(.cgb(torch_code(Determinant(M))), "torch.det(M)")
  expect_identical(.cgb(torch_code(HadamardProduct(M, N))), "torch.mul(M, N)")
  expect_identical(.cgb(torch_code(Inverse(M))), "torch.linalg.inv(M)")
})

test_that("printing_test_torch: test_torch_array_operations", {
  M <- MatrixSymbol("M", 3L, 3L); N <- MatrixSymbol("N", 3L, 3L); P <- MatrixSymbol("P", 3L, 3L); Q <- MatrixSymbol("Q", 3L, 3L)
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); t <- Symbol("t")
  suppressWarnings(try({ M <- MatrixSymbol("M", 2L, 2L) }, silent = TRUE))
  suppressWarnings(try({ N <- MatrixSymbol("N", 2L, 2L) }, silent = TRUE))
  suppressWarnings(try({ P <- MatrixSymbol("P", 2L, 2L) }, silent = TRUE))
  suppressWarnings(try({ Q <- MatrixSymbol("Q", 2L, 2L) }, silent = TRUE))
  expect_identical(.cgb(torch_code(ArrayTensorProduct(M, N))), "torch.einsum(\"ab,cd\", M, N)")
  expect_identical(.cgb(torch_code(ArrayAdd(M, N))), "torch.add(M, N)")
  expect_identical(.cgb(torch_code(ArrayAdd(M, N, P))), "torch.add(torch.add(M, N), P)")
  expect_identical(.cgb(torch_code(ArrayAdd(M, N, P, Q))), "torch.add(torch.add(torch.add(M, N), P), Q)")
  expect_identical(.cgb(torch_code(PermuteDims(M, list(1L, 0L)))), "M.permute(1, 0)")
  expect_identical(.cgb(torch_code(PermuteDims(ArrayTensorProduct(M, N), list(1L, 2L, 3L, 0L)))), "torch.einsum(\"ab,cd\", M, N).permute(1, 2, 3, 0)")
  expect_identical(.cgb(torch_code(ArrayDiagonal(ArrayTensorProduct(M, N), list(1L, 2L)))), "torch.einsum(\"ab,bc->acb\", M, N)")
})

test_that("printing_test_torch: test_torch_derivative", {
  M <- MatrixSymbol("M", 3L, 3L); N <- MatrixSymbol("N", 3L, 3L); P <- MatrixSymbol("P", 3L, 3L); Q <- MatrixSymbol("Q", 3L, 3L)
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); t <- Symbol("t")
  expect_identical(.cgb(torch_code(Derivative(sin(x), x))), "torch.autograd.grad(torch.sin(x), x)[0]")
})

test_that("printing_test_torch: test_torch_printing_dtype", {
  M <- MatrixSymbol("M", 3L, 3L); N <- MatrixSymbol("N", 3L, 3L); P <- MatrixSymbol("P", 3L, 3L); Q <- MatrixSymbol("Q", 3L, 3L)
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); t <- Symbol("t")
  suppressWarnings(try({ expr <- Matrix(list(list(x, sin(y)), list(exp(z), -t))) }, silent = TRUE))
  expect_identical(.cgb(torch_code(expr)), "torch.tensor([[x, torch.sin(y)], [torch.exp(z), -t]], dtype=torch.float64)")
  expect_identical(.cgb(torch_code(expr, dtype = "torch.float32")), "torch.tensor([[x, torch.sin(y)], [torch.exp(z), -t]], dtype=torch.float32)")
  expect_identical(.cgb(torch_code(expr, requires_grad = TRUE)), "torch.tensor([[x, torch.sin(y)], [torch.exp(z), -t]], dtype=torch.float64, requires_grad=True)")
  expect_identical(.cgb(torch_code(expr, requires_grad = TRUE, dtype = "torch.float32")), "torch.tensor([[x, torch.sin(y)], [torch.exp(z), -t]], dtype=torch.float32, requires_grad=True)")
})

test_that("printing_test_torch: test_requires_grad", {
  M <- MatrixSymbol("M", 3L, 3L); N <- MatrixSymbol("N", 3L, 3L); P <- MatrixSymbol("P", 3L, 3L); Q <- MatrixSymbol("Q", 3L, 3L)
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); t <- Symbol("t")
  expect_identical(.cgb(torch_code(sin(x) + cos(y))), "torch.sin(x) + torch.cos(y)")
})

test_that("printing_test_torch: test_torch_multi_variable_derivatives", {
  M <- MatrixSymbol("M", 3L, 3L); N <- MatrixSymbol("N", 3L, 3L); P <- MatrixSymbol("P", 3L, 3L); Q <- MatrixSymbol("Q", 3L, 3L)
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); t <- Symbol("t")
  suppressWarnings(try({ .tmp <- sym_symbols("x y z"); x <- py_getitem(.tmp, 0L); y <- py_getitem(.tmp, 1L); z <- py_getitem(.tmp, 2L) }, silent = TRUE))
  expect_identical(.cgb(torch_code(Derivative(sin(x), x))), "torch.autograd.grad(torch.sin(x), x)[0]")
  expect_identical(.cgb(torch_code(Derivative(sin(x), list(x, 2L)))), "torch.autograd.grad(torch.autograd.grad(torch.sin(x), x, create_graph=True)[0], x, create_graph=True)[0]")
  expect_identical(.cgb(torch_code(Derivative(sin(x * y), x, y))), "torch.autograd.grad(torch.autograd.grad(torch.sin(x*y), x, create_graph=True)[0], y, create_graph=True)[0]")
  expect_identical(.cgb(torch_code(Derivative(sin(x), x, x))), "torch.autograd.grad(torch.autograd.grad(torch.sin(x), x, create_graph=True)[0], x, create_graph=True)[0]")
  expect_identical(.cgb(torch_code(Derivative(sin((x * y) * z), x, list(y, 2L), z))), "torch.autograd.grad(torch.autograd.grad(torch.autograd.grad(torch.autograd.grad(torch.sin(x*y*z), x, create_graph=True)[0], y, create_graph=True)[0], y, create_graph=True)[0], z, create_graph=True)[0]")
})

test_that("printing_test_torch: test_torch_derivative_lambdify", {
  M <- MatrixSymbol("M", 3L, 3L); N <- MatrixSymbol("N", 3L, 3L); P <- MatrixSymbol("P", 3L, 3L); Q <- MatrixSymbol("Q", 3L, 3L)
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); t <- Symbol("t")
  suppressWarnings(try({ x <- sym_symbols("x") }, silent = TRUE))
  suppressWarnings(try({ y <- sym_symbols("y") }, silent = TRUE))
  expect_identical(.cgb(torch_code(Derivative(x ^ 2L, x))), "torch.autograd.grad(torch.pow(x, 2), x)[0]")
  expect_identical(.cgb(torch_code(Derivative(sin(x), list(x, 2L)))), "torch.autograd.grad(torch.autograd.grad(torch.sin(x), x, create_graph=True)[0], x, create_graph=True)[0]")
  expect_identical(.cgb(torch_code(Derivative(x * (y ^ 2L), x, y))), "torch.autograd.grad(torch.autograd.grad(x*torch.pow(y, 2), x, create_graph=True)[0], y, create_graph=True)[0]")
})

test_that("printing_test_torch: test_torch_special_matrices", {
  M <- MatrixSymbol("M", 3L, 3L); N <- MatrixSymbol("N", 3L, 3L); P <- MatrixSymbol("P", 3L, 3L); Q <- MatrixSymbol("Q", 3L, 3L)
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); t <- Symbol("t")
  expect_identical(.cgb(torch_code(Identity(3L))), "torch.eye(3)")
  suppressWarnings(try({ n <- sym_symbols("n") }, silent = TRUE))
  expect_identical(.cgb(torch_code(Identity(n))), "torch.eye(n, n)")
  expect_identical(.cgb(torch_code(ZeroMatrix(2L, 3L))), "torch.zeros((2, 3))")
  suppressWarnings(try({ .tmp <- sym_symbols("m n"); m <- py_getitem(.tmp, 0L); n <- py_getitem(.tmp, 1L) }, silent = TRUE))
  expect_identical(.cgb(torch_code(ZeroMatrix(m, n))), "torch.zeros((m, n))")
  expect_identical(.cgb(torch_code(OneMatrix(2L, 3L))), "torch.ones((2, 3))")
  expect_identical(.cgb(torch_code(OneMatrix(m, n))), "torch.ones((m, n))")
})

test_that("printing_test_torch: test_torch_complex_operations", {
  M <- MatrixSymbol("M", 3L, 3L); N <- MatrixSymbol("N", 3L, 3L); P <- MatrixSymbol("P", 3L, 3L); Q <- MatrixSymbol("Q", 3L, 3L)
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); t <- Symbol("t")
  expect_identical(.cgb(torch_code(conjugate(x))), "torch.conj(x)")
  expect_identical(.cgb(torch_code(conjugate(sin(x) + (S_I * cos(y))))), "torch.sin(torch.conj(x)) - 1j*torch.cos(torch.conj(y))")
  expect_identical(.cgb(torch_code(S_I)), "1j")
  expect_identical(.cgb(torch_code((2L * S_I) + x)), "x + 2*1j")
  expect_identical(.cgb(torch_code(exp(S_I * x))), "torch.exp(1j*x)")
})

test_that("printing_test_torch: test_torch_special_functions", {
  M <- MatrixSymbol("M", 3L, 3L); N <- MatrixSymbol("N", 3L, 3L); P <- MatrixSymbol("P", 3L, 3L); Q <- MatrixSymbol("Q", 3L, 3L)
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); t <- Symbol("t")
  expect_identical(.cgb(torch_code(Heaviside(x))), "torch.heaviside(x, 1/2)")
  expect_identical(.cgb(torch_code(Heaviside(x, 0L))), "torch.heaviside(x, 0)")
  expect_identical(.cgb(torch_code(gamma(x))), "torch.special.gamma(x)")
  expect_identical(.cgb(torch_code(polygamma(0L, x))), "torch.special.digamma(x)")
  expect_identical(.cgb(torch_code(gamma(sin(x)))), "torch.special.gamma(torch.sin(x))")
  expect_identical(.cgb(torch_code(polygamma(1L, x))), "ERR NotImplementedError: PyTorch only supports digamma (0th order polygamma)")
})

test_that("printing_test_torch: torch_extra", {
  M <- MatrixSymbol("M", 3L, 3L); N <- MatrixSymbol("N", 3L, 3L); P <- MatrixSymbol("P", 3L, 3L); Q <- MatrixSymbol("Q", 3L, 3L)
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); t <- Symbol("t")
  expect_identical(.cgb(torch_code(Piecewise(list(x, x < 0L), list(y, x < 1L), list(z, TRUE)))), "torch.where(torch.lt(x, 0), x, torch.where(torch.lt(x, 1), y, z))")
  expect_identical(.cgb(torch_code(Piecewise(list(x, x < 0L)))), "torch.where(torch.lt(x, 0), x, 0)")
  expect_identical(.cgb(torch_code(And(x < 1L, y > 2L))), "torch.logical_and(torch.gt(y, 2), torch.lt(x, 1))")
  expect_identical(.cgb(torch_code(Or(x < 1L, y > 2L, z > 3L))), "torch.logical_or(torch.logical_or(torch.gt(y, 2), torch.gt(z, 3)), torch.lt(x, 1))")
  expect_identical(.cgb(torch_code(Not(x < 1L))), "torch.ge(x, 1)")
  expect_identical(.cgb(torch_code(Max(x, y, z))), "torch.maximum(torch.maximum(x, y), z)")
  expect_identical(.cgb(torch_code(x ^ y)), "torch.pow(x, y)")
  expect_identical(.cgb(torch_code((2L * M) * N)), "2*torch.matmul(M, N)")
  expect_identical(.cgb(torch_code(Transpose(M))), "torch.transpose(M).t()")
  expect_identical(.cgb(torch_code(Transpose(MatrixSymbol("R", 2L, 3L)))), "R.permute(1 ,0)")
  expect_identical(.cgb(torch_code(Piecewise(list(Float(0.0), isnan(x)), list(Float(1.0), TRUE)))), "torch.where(torch.isnan(x), 0.0, 1.0)")
  expect_identical(.cgb(torch_code(Piecewise(list(Float(0.0), isinf(x)), list(Float(1.0), TRUE)))), "torch.where(torch.isinf(x), 0.0, 1.0)")
  expect_identical(.cgb(torch_code(Assignment(y, sin(x)))), "y = torch.sin(x)")
  expect_identical(.cgb(torch_code(1L / x)), "torch.pow(x, -1)")
  expect_identical(.cgb(torch_code(log(x) + exp(y))), "torch.exp(y) + torch.log(x)")
  expect_identical(.cgb(torch_code(ZeroArray(3L, 2L))), "torch.zeros((3,2,))")
  expect_identical(.cgb(torch_code(OneArray(4L))), "torch.ones((4,))")
  expect_identical(.cgb(torch_code(ArrayContraction(ArrayTensorProduct(M, N), list(1L, 2L)))), "torch.einsum(\"ab,bc->ac\", M,N)")
  expect_identical(.cgb(torch_code(ArrayElementwiseApplyFunc(sin, M))), "torch.sin(M)")
  expect_identical(.cgb(ast_attr(TorchPrinter(), "printmethod")), "_torchcode")
  expect_identical(.cgb(ast_attr(TorchPrinter(py_dict(list("requires_grad"), list(TRUE))), "requires_grad")), "True")
  expect_identical(.cgb(ast_attr(TorchPrinter(), "dtype")), "torch.float64")
  expect_identical(.cgb(torch_code(Identity(3L) + M)), "torch.add(torch.eye(3), M)")
})

test_that("printing_test_lambdarepr: test_basic", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); i <- Symbol("i"); a <- Symbol("a"); b <- Symbol("b"); j <- Symbol("j"); c <- Symbol("c"); d <- Symbol("d")
  expect_identical(.cgb(lambdarepr(x * y)), "x*y")
  expect_identical(.cgb(lambdarepr(x + y)), "x + y")
  expect_identical(.cgb(lambdarepr(x ^ y)), "x^y")
})

test_that("printing_test_lambdarepr: test_matrix", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); i <- Symbol("i"); a <- Symbol("a"); b <- Symbol("b"); j <- Symbol("j"); c <- Symbol("c"); d <- Symbol("d")
  suppressWarnings(try({ e <- x %% 2L }, silent = TRUE))
  expect_identical(.cgb(lambdarepr(e)), "x % 2")
  expect_identical(.cgb(py_str(e)), "Mod(x, 2)")
  expect_identical(.cgb(lambdarepr(Matrix(list(e)))), "ImmutableDenseMatrix([[x % 2]])")
})

test_that("printing_test_lambdarepr: test_piecewise", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); i <- Symbol("i"); a <- Symbol("a"); b <- Symbol("b"); j <- Symbol("j"); c <- Symbol("c"); d <- Symbol("d")
  expect_identical(.cgb(lambdarepr(Piecewise(list(x, x < 0L)))), "((x) if (x < 0) else None)")
  expect_identical(.cgb(lambdarepr(Piecewise(list(1L, x < 1L), list(2L, x < 2L), list(0L, TRUE)))), "((1) if (x < 1) else (2) if (x < 2) else (0))")
  expect_identical(.cgb(lambdarepr(Piecewise(list(1L, x < 1L), list(2L, x < 2L)))), "((1) if (x < 1) else (2) if (x < 2) else None)")
  expect_identical(.cgb(lambdarepr(Piecewise(list(x, x < 1L), list(x ^ 2L, contains(Interval(3L, 4L, TRUE, FALSE), x)), list(0L, TRUE)))), "((x) if (x < 1) else (x^2) if (((x <= 4)) and ((x > 3))) else (0))")
  expect_identical(.cgb(lambdarepr(Piecewise(list(x ^ 2L, x < 0L), list(x, x < 1L), list(2L - x, x >= 1L), list(0L, TRUE), evaluate = FALSE))), "((x^2) if (x < 0) else (x) if (x < 1) else (2 - x) if (x >= 1) else (0))")
  expect_identical(.cgb(lambdarepr(Piecewise(list(x ^ 2L, x < 0L), list(x, x < 1L), list(2L - x, x >= 1L), evaluate = FALSE))), "((x^2) if (x < 0) else (x) if (x < 1) else (2 - x) if (x >= 1) else None)")
  expect_identical(.cgb(lambdarepr(Piecewise(list(1L, x >= 1L), list(2L, x >= 2L), list(3L, x >= 3L), list(4L, x >= 4L), list(5L, x >= 5L), list(6L, TRUE)))), "((1) if (x >= 1) else (2) if (x >= 2) else (3) if (x >= 3) else (4) if (x >= 4) else (5) if (x >= 5) else (6))")
  expect_identical(.cgb(lambdarepr(Piecewise(list(1L, x <= 1L), list(2L, x <= 2L), list(3L, x <= 3L), list(4L, x <= 4L), list(5L, x <= 5L), list(6L, TRUE)))), "((1) if (x <= 1) else (2) if (x <= 2) else (3) if (x <= 3) else (4) if (x <= 4) else (5) if (x <= 5) else (6))")
  expect_identical(.cgb(lambdarepr(Piecewise(list(1L, x > 1L), list(2L, x > 2L), list(3L, x > 3L), list(4L, x > 4L), list(5L, x > 5L), list(6L, TRUE)))), "((1) if (x > 1) else (2) if (x > 2) else (3) if (x > 3) else (4) if (x > 4) else (5) if (x > 5) else (6))")
  expect_identical(.cgb(lambdarepr(Piecewise(list(1L, x < 1L), list(2L, x < 2L), list(3L, x < 3L), list(4L, x < 4L), list(5L, x < 5L), list(6L, TRUE)))), "((1) if (x < 1) else (2) if (x < 2) else (3) if (x < 3) else (4) if (x < 4) else (5) if (x < 5) else (6))")
  expect_identical(.cgb(lambdarepr(Piecewise(list(Piecewise(list(1L, x > 0L), list(2L, TRUE)), y > 0L), list(3L, TRUE)))), "((((1) if (x > 0) else (2))) if (y > 0) else (3))")
})

test_that("printing_test_lambdarepr: test_sum__1", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); i <- Symbol("i"); a <- Symbol("a"); b <- Symbol("b"); j <- Symbol("j"); c <- Symbol("c"); d <- Symbol("d")
  expect_identical(.cgb(lambdarepr(Sum(x ^ i, list(i, a, b)))), "(builtins.sum(x^i for i in range(a, b+1)))")
})

test_that("printing_test_lambdarepr: test_sum__2", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); i <- Symbol("i"); a <- Symbol("a"); b <- Symbol("b"); j <- Symbol("j"); c <- Symbol("c"); d <- Symbol("d")
  expect_identical(.cgb(lambdarepr(Sum(i * x, list(i, a, b)))), "(builtins.sum(i*x for i in range(a, b+1)))")
})

test_that("printing_test_lambdarepr: test_multiple_sums", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); i <- Symbol("i"); a <- Symbol("a"); b <- Symbol("b"); j <- Symbol("j"); c <- Symbol("c"); d <- Symbol("d")
  expect_identical(.cgb(lambdarepr(Sum((i * x) + j, list(i, a, b), list(j, c, d)))), "(builtins.sum(i*x + j for j in range(c, d+1) for i in range(a, b+1)))")
})

test_that("printing_test_lambdarepr: test_sqrt", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); i <- Symbol("i"); a <- Symbol("a"); b <- Symbol("b"); j <- Symbol("j"); c <- Symbol("c"); d <- Symbol("d")
  suppressWarnings(try({ prntr <- LambdaPrinter(py_dict(list("standard"), list("python3"))) }, silent = TRUE))
  expect_identical(.cgb(py_method(prntr, "_print_Pow", sqrt(x), rational = FALSE)), "sqrt(x)")
  expect_identical(.cgb(py_method(prntr, "_print_Pow", sqrt(x), rational = TRUE)), "x^(1/2)")
})

test_that("printing_test_lambdarepr: test_settings", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); i <- Symbol("i"); a <- Symbol("a"); b <- Symbol("b"); j <- Symbol("j"); c <- Symbol("c"); d <- Symbol("d")
  expect_identical(.cgb(lambdarepr(sin(x), method = "garbage")), "ERR TypeError: Unknown setting 'method'.")
})

test_that("printing_test_lambdarepr: test_numexpr", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); i <- Symbol("i"); a <- Symbol("a"); b <- Symbol("b"); j <- Symbol("j"); c <- Symbol("c"); d <- Symbol("d")
  suppressWarnings(try({ expr <- new_sym("ITE", list(x > 0L, .as_bool(TRUE), .as_bool(FALSE))) }, silent = TRUE))
  expect_identical(.cgb(printer_doprint(NumExprPrinter(), expr)), "numexpr.evaluate('where((x > 0), True, False)', truediv=True)")
  suppressWarnings(try({ func_def <- FunctionDefinition(NULL, "foo", list(Variable(x)), list(Assignment(y, x), Return(y ^ 2L))) }, silent = TRUE))
  expect_identical(.cgb(printer_doprint(NumExprPrinter(), func_def)), "def foo(x):\\n    y = numexpr.evaluate('x', truediv=True)\\n    return numexpr.evaluate('y^2', truediv=True)")
})

test_that("printing_test_lambdarepr: test_printmethod", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); i <- Symbol("i"); a <- Symbol("a"); b <- Symbol("b"); j <- Symbol("j"); c <- Symbol("c"); d <- Symbol("d")
  expect_identical(.cgb(printer_doprint(NumExprPrinter(), Piecewise(list(y, x >= 0L), list(z, x < 0L)))), "numexpr.evaluate('where((x >= 0), y, z)', truediv=True)")
})

test_that("printing_test_lambdarepr: lambdarepr_extra", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); i <- Symbol("i"); a <- Symbol("a"); b <- Symbol("b"); j <- Symbol("j"); c <- Symbol("c"); d <- Symbol("d")
  expect_identical(.cgb(lambdarepr(And(x, y))), "((x) and (y))")
  expect_identical(.cgb(lambdarepr(Or(x, y, z))), "((x) or (y) or (z))")
  expect_identical(.cgb(lambdarepr(Not(x))), "(not (x))")
  expect_identical(.cgb(lambdarepr(ITE(x > 0L, y, z))), "((y) if ((x > 0)) else (z))")
  expect_identical(.cgb(lambdarepr(S_Pi + S_E)), "math.e + math.pi")
  expect_identical(.cgb(lambdarepr(S_true)), "True")
  expect_identical(.cgb(lambdarepr(sqrt(x))), "sqrt(x)")
  expect_identical(.cgb(lambdarepr(x ^ Rational(1L, 3L))), "x^(1/3)")
  expect_identical(.cgb(lambdarepr(sin(x) + exp(y))), "math.exp(y) + math.sin(x)")
  expect_identical(.cgb(printer_doprint(NumExprPrinter(), sin(x) + cos(y))), "numexpr.evaluate('sin(x) + cos(y)', truediv=True)")
  expect_identical(.cgb(printer_doprint(NumExprPrinter(), asin(x) * sym_atan2(x, y))), "numexpr.evaluate('arcsin(x)*arctan2(x, y)', truediv=True)")
  expect_identical(.cgb(printer_doprint(NumExprPrinter(), S_I * x)), "numexpr.evaluate('1j*x', truediv=True)")
  expect_identical(.cgb(printer_doprint(NumExprPrinter(), ((Abs(x) + conjugate(y)) + re(x)) + im(y))), "numexpr.evaluate('conj(y) + real(x) + imag(y) + abs(x)', truediv=True)")
  expect_identical(.cgb(printer_doprint(NumExprPrinter(), (log(x) + exp(y)) + sqrt(z))), "numexpr.evaluate('sqrt(z) + exp(y) + log(x)', truediv=True)")
  expect_identical(.cgb(printer_doprint(NumExprPrinter(), gamma(x))), "numexpr.evaluate('math.gamma(x)', truediv=True)")
  expect_identical(.cgb(printer_doprint(NumExprPrinter(), Piecewise(list(x, x < 0L), list(y, x < 1L)))), "numexpr.evaluate('where((x < 0), x, where((x < 1), y, log(-1)))', truediv=True)")
  expect_identical(.cgb(printer_doprint(NumExprPrinter(), Piecewise(list(x, x < 0L), list(y, x < 1L), list(z, TRUE)))), "numexpr.evaluate('where((x < 0), x, where((x < 1), y, z))', truediv=True)")
  expect_identical(.cgb(printer_doprint(NumExprPrinter(), Assignment(y, x + 1L))), "y = numexpr.evaluate('x + 1', truediv=True)")
  expect_identical(.cgb(printer_doprint(NumExprPrinter(), (x ^ 2L) / y)), "numexpr.evaluate('x^2/y', truediv=True)")
  expect_identical(.cgb(ast_attr(NumExprPrinter(), "printmethod")), "_numexprcode")
  expect_identical(.cgb(printer_doprint(IntervalPrinter(), x + 1L)), "x + mpi('1')")
  expect_identical(.cgb(printer_doprint(IntervalPrinter(), (x / 2L) + Rational(3L, 4L))), "(mpi('1/2'))*x + mpi('3/4')")
  expect_identical(.cgb(printer_doprint(IntervalPrinter(), x ^ 2L)), "x^mpi('2')")
  expect_identical(.cgb(printer_doprint(IntervalPrinter(), sqrt(x))), "x^(mpi('1/2'))")
  expect_identical(.cgb(printer_doprint(IntervalPrinter(), (x ^ Rational(1L, 2L)) + S_Half)), "x^(mpi('1/2')) + mpi('1/2')")
  expect_identical(.cgb(printer_doprint(IntervalPrinter(), sin(x) * 3L)), "mpi('3')*mpmath.sin(x)")
  expect_identical(.cgb(printer_doprint(IntervalPrinter(), S_Pi * x)), "mpmath.pi*x")
  expect_identical(.cgb(printer_doprint(IntervalPrinter(), 1L / x)), "mpi('1')/x")
  expect_identical(.cgb(printer_doprint(IntervalPrinter(), And(x > 1L, y < 2L))), "(((x > mpi('1'))) and ((y < mpi('2'))))")
  expect_identical(.cgb(printer_doprint(IntervalPrinter(), exp(x) + Float(2.5))), "mpmath.exp(x) + mpmath.mpf((0, 5, -1, 3))")
  expect_identical(.cgb(ast_attr(IntervalPrinter(), "printmethod")), "_mpmathcode")
})

test_that("printing_test_pycode: test_array_printer", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); p <- IndexedBase("p")
  suppressWarnings(try({ A <- ArraySymbol("A", list(4L, 4L, 6L, 6L, 6L)) }, silent = TRUE))
  suppressWarnings(try({ S_I <- IndexedBase("I") }, silent = TRUE))
  suppressWarnings(try({ i <- Idx("i", list(0L, 1L)); j <- Idx("j", list(2L, 3L)); k <- Idx("k", list(4L, 5L)) }, silent = TRUE))
  suppressWarnings(try({ prntr <- NumPyPrinter() }, silent = TRUE))
  expect_identical(.cgb(printer_doprint(prntr, ZeroArray(5L))), "numpy.zeros((5,))")
  expect_identical(.cgb(printer_doprint(prntr, OneArray(5L))), "numpy.ones((5,))")
  expect_identical(.cgb(printer_doprint(prntr, ArrayContraction(A, list(2L, 3L)))), "numpy.einsum(\"abccd->abd\", A)")
  expect_identical(.cgb(printer_doprint(prntr, S_I)), "I")
  expect_identical(.cgb(printer_doprint(prntr, ArrayDiagonal(A, list(2L, 3L, 4L)))), "numpy.einsum(\"abccc->abc\", A)")
  expect_identical(.cgb(printer_doprint(prntr, ArrayDiagonal(A, list(0L, 1L), list(2L, 3L)))), "numpy.einsum(\"aabbc->cab\", A)")
  expect_identical(.cgb(printer_doprint(prntr, ArrayContraction(A, list(2L), list(3L)))), "numpy.einsum(\"abcde->abe\", A)")
  expect_identical(.cgb(printer_doprint(prntr, Assignment(py_getitem(S_I, i, j, k), py_getitem(S_I, i, j, k)))), "I = I")
  suppressWarnings(try({ B <- ArraySymbol("B", list(3L, 3L)) }, silent = TRUE))
  suppressWarnings(try({ d <- Dummy("d") }, silent = TRUE))
  expect_identical(.cgb(printer_doprint(prntr, ArrayElementwiseApplyFunc(sin, B))), "numpy.sin(B)")
  expect_identical(.cgb(printer_doprint(prntr, ArrayElementwiseApplyFunc(Lambda(d, d ^ 2L), B))), "B^2")
  expect_identical(.cgb(printer_doprint(prntr, ArrayElementwiseApplyFunc(Lambda(d, 1L / d), ArrayContraction(ArrayTensorProduct(B, B), list(1L, 2L))))), "numpy.einsum(\"ab,bc->ac\", B,B)^(-1.0)")
  suppressWarnings(try({ s <- Symbol("_elementwise_operand") }, silent = TRUE))
  expect_identical(.cgb(printer_doprint(prntr, ArrayElementwiseApplyFunc(Lambda(d, d + s), B))), "B + _elementwise_operand")
  suppressWarnings(try({ s <- Symbol("prefix_elementwise_operand_suffix") }, silent = TRUE))
  expect_identical(.cgb(printer_doprint(prntr, ArrayElementwiseApplyFunc(Lambda(d, d + s), B))), "B + prefix_elementwise_operand_suffix")
  suppressWarnings(try({ nested <- ArrayElementwiseApplyFunc(Lambda(d, d + 1L), ArrayElementwiseApplyFunc(Lambda(d, d ^ 2L), B)) }, silent = TRUE))
  expect_identical(.cgb(printer_doprint(prntr, nested)), "B^2 + 1")
  expect_identical(.cgb(printer_doprint(prntr, ArrayElementwiseApplyFunc(Lambda(d, d ^ 3L), ArrayElementwiseApplyFunc(Lambda(d, d ^ 2L), B)))), "(B^2)^3")
  expect_identical(.cgb(printer_doprint(prntr, ArrayElementwiseApplyFunc(sin, ArrayElementwiseApplyFunc(Lambda(d, d + 1L), B)))), "numpy.sin(B + 1)")
  suppressWarnings(try({ prntr <- TensorflowPrinter() }, silent = TRUE))
  expect_identical(.cgb(printer_doprint(prntr, ZeroArray(5L))), "tensorflow.zeros((5,))")
  expect_identical(.cgb(printer_doprint(prntr, OneArray(5L))), "tensorflow.ones((5,))")
  expect_identical(.cgb(printer_doprint(prntr, ArrayContraction(A, list(2L, 3L)))), "tensorflow.linalg.einsum(\"abccd->abd\", A)")
  expect_identical(.cgb(printer_doprint(prntr, S_I)), "I")
  expect_identical(.cgb(printer_doprint(prntr, ArrayDiagonal(A, list(2L, 3L, 4L)))), "tensorflow.linalg.einsum(\"abccc->abc\", A)")
  expect_identical(.cgb(printer_doprint(prntr, ArrayDiagonal(A, list(0L, 1L), list(2L, 3L)))), "tensorflow.linalg.einsum(\"aabbc->cab\", A)")
  expect_identical(.cgb(printer_doprint(prntr, ArrayContraction(A, list(2L), list(3L)))), "tensorflow.linalg.einsum(\"abcde->abe\", A)")
  expect_identical(.cgb(printer_doprint(prntr, Assignment(py_getitem(S_I, i, j, k), py_getitem(S_I, i, j, k)))), "I = I")
})

test_that("printing_test_pycode: array_printer_extra", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); p <- IndexedBase("p")
  suppressWarnings(try({ A <- ArraySymbol("A", list(4L, 4L, 6L, 6L, 6L)) }, silent = TRUE))
  suppressWarnings(try({ B <- ArraySymbol("B", list(3L, 3L)) }, silent = TRUE))
  suppressWarnings(try({ C <- ArraySymbol("C", list(3L, 3L, 3L)) }, silent = TRUE))
  suppressWarnings(try({ S_I <- IndexedBase("I") }, silent = TRUE))
  suppressWarnings(try({ J <- IndexedBase("J") }, silent = TRUE))
  suppressWarnings(try({ i <- Idx("i", list(0L, 2L)); j <- Idx("j", list(0L, 2L)); k <- Idx("k", list(0L, 2L)) }, silent = TRUE))
  suppressWarnings(try({ M <- MatrixSymbol("M", 3L, 3L) }, silent = TRUE))
  suppressWarnings(try({ N <- MatrixSymbol("N", 3L, 3L) }, silent = TRUE))
  suppressWarnings(try({ np_p <- NumPyPrinter() }, silent = TRUE))
  suppressWarnings(try({ cp <- CuPyPrinter() }, silent = TRUE))
  suppressWarnings(try({ jp <- JaxPrinter() }, silent = TRUE))
  suppressWarnings(try({ tp <- TorchPrinter() }, silent = TRUE))
  expect_identical(.cgb(printer_doprint(cp, ZeroArray(5L))), "cupy.zeros((5,))")
  expect_identical(.cgb(printer_doprint(jp, OneArray(5L, 3L))), "jax.numpy.ones((5,3,))")
  expect_identical(.cgb(printer_doprint(tp, ZeroArray(2L, 3L))), "torch.zeros((2,3,))")
  expect_identical(.cgb(printer_doprint(tp, OneArray(5L))), "torch.ones((5,))")
  expect_identical(.cgb(printer_doprint(cp, ArrayContraction(A, list(2L, 3L)))), "cupy.einsum(\"abccd->abd\", A)")
  expect_identical(.cgb(printer_doprint(jp, ArrayDiagonal(A, list(0L, 1L), list(2L, 3L)))), "jax.numpy.einsum(\"aabbc->cab\", A)")
  expect_identical(.cgb(printer_doprint(tp, ArrayContraction(A, list(2L, 3L)))), "torch.einsum(\"abccd->abd\", A)")
  expect_identical(.cgb(printer_doprint(tp, ArrayDiagonal(A, list(2L, 3L, 4L)))), "torch.einsum(\"abccc->abc\", A)")
  expect_identical(.cgb(printer_doprint(np_p, ArrayTensorProduct(B, C))), "numpy.einsum(\"ab,cde\", B, C)")
  expect_identical(.cgb(printer_doprint(np_p, ArrayTensorProduct(B, C, B))), "numpy.einsum(\"ab,cde,fg\", B, C, B)")
  expect_identical(.cgb(printer_doprint(np_p, ArrayContraction(ArrayTensorProduct(B, C), list(0L, 2L), list(1L, 4L)))), "numpy.einsum(\"ab,acb->c\", B,C)")
  expect_identical(.cgb(printer_doprint(np_p, ArrayContraction(ArrayTensorProduct(B, B, B), list(1L, 2L), list(3L, 4L)))), "numpy.einsum(\"ab,bc,cd->ad\", B,B,B)")
  expect_identical(.cgb(printer_doprint(np_p, ArrayContraction(B, list(0L, 1L)))), "numpy.einsum(\"aa->\", B)")
  expect_identical(.cgb(printer_doprint(np_p, ArrayDiagonal(ArrayTensorProduct(B, C), list(0L, 2L)))), "numpy.einsum(\"ab,acd->bcda\", B, C)")
  expect_identical(.cgb(printer_doprint(np_p, ArrayDiagonal(C, list(0L, 1L, 2L)))), "numpy.einsum(\"aaa->a\", C)")
  expect_identical(.cgb(printer_doprint(np_p, PermuteDims(C, list(2L, 0L, 1L)))), "numpy.transpose(C, [2, 0, 1])")
  expect_identical(.cgb(printer_doprint(np_p, PermuteDims(ArrayTensorProduct(B, C), list(4L, 3L, 2L, 1L, 0L)))), "numpy.transpose(numpy.einsum(\"ab,cde\", B, C), [4, 3, 2, 1, 0])")
  expect_identical(.cgb(printer_doprint(np_p, ArrayAdd(B, PermuteDims(B, list(1L, 0L))))), "numpy.add(B, numpy.transpose(B, [1, 0]))")
  expect_identical(.cgb(printer_doprint(cp, ArrayAdd(B, B, PermuteDims(B, list(1L, 0L))))), "cupy.add(cupy.add(B, B), cupy.transpose(B, [1, 0]))")
  expect_identical(.cgb(printer_doprint(tp, PermuteDims(C, list(2L, 0L, 1L)))), "C.permute(2, 0, 1)")
  expect_identical(.cgb(printer_doprint(tp, ArrayAdd(M, N))), "torch.add(M, N)")
  expect_identical(.cgb(printer_doprint(np_p, convert_matrix_to_array((M * N) * M))), "numpy.einsum(\"ab,ca,bd->cd\", N,M,M)")
  expect_identical(.cgb(printer_doprint(np_p, convert_matrix_to_array(Transpose(M) * N))), "numpy.einsum(\"ab,ac->bc\", M,N)")
  expect_identical(.cgb(printer_doprint(np_p, convert_matrix_to_array(M + Transpose(N)))), "numpy.add(M, numpy.transpose(N, [1, 0]))")
  expect_identical(.cgb(printer_doprint(np_p, convert_matrix_to_array(Trace(M * N)))), "numpy.einsum(\"ab,ba->\", M,N)")
  expect_identical(.cgb(printer_doprint(np_p, Assignment(py_getitem(S_I, i, j), py_getitem(J, j, i)))), "I = J")
  expect_identical(.cgb(printer_doprint(np_p, S_I)), "I")
  expect_identical(.cgb(printer_doprint(np_p, py_getitem(S_I, i, j))), "I[i, j]")
  expect_identical(.cgb(printer_doprint(np_p, Assignment(x, y + z))), "x = y + z")
  expect_identical(.cgb(printer_doprint(np_p, ArrayElementwiseApplyFunc(cos, ArrayTensorProduct(B, B)))), "numpy.cos(numpy.einsum(\"ab,cd\", B, B))")
  expect_identical(.cgb(printer_doprint(cp, ArrayElementwiseApplyFunc(exp, B))), "cupy.exp(B)")
  expect_identical(.cgb(printer_doprint(np_p, ArrayTensorProduct(2L, B))), "numpy.einsum(\",ab\", 2, B)")
  expect_identical(.cgb(printer_doprint(np_p, ArrayContraction(ArrayTensorProduct(A, B), list(5L, 6L)))), "numpy.einsum(\"abcde,ff->abcde\", A,B)")
})

test_that("printing_test_pycode: list_minmax_extra", {
  x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); p <- IndexedBase("p")
  suppressWarnings(try({ x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z") }, silent = TRUE))
  expect_identical(.cgb(printer_doprint(NumPyPrinter(), np_maximum(x, y))), "functools.reduce(numpy.maximum, [x, y])")
  expect_identical(.cgb(printer_doprint(NumPyPrinter(), np_minimum(x, y))), "functools.reduce(numpy.minimum, [x, y])")
  expect_identical(.cgb(printer_doprint(NumPyPrinter(), np_maximum(x, y, z))), "functools.reduce(numpy.maximum, [x, y, z])")
  expect_identical(.cgb(pycode(List(1L, 2L, 3L))), "[1, 2, 3]")
  expect_identical(.cgb(pycode(List())), "[]")
  expect_identical(.cgb(printer_doprint(NumPyPrinter(), List(x, y))), "[x, y]")
  expect_identical(.cgb(pycode(List(x, List(1L, 2L)))), "[x, [1, 2]]")
})

test_that("branch_test_core: test_exhaust", {
  posdec <- function(x) if (x > 0) list(x - 1L) else list(x)
  branch5 <- function(x) {
    if (0 < x && x < 5) list(x - 1L)
    else if (5 < x && x < 10) list(x + 1L)
    else if (x == 5) list(x + 1L, x - 1L)
    else list(x)
  }
  even <- function(x) x %% 2L == 0L
  inc <- function(x) list(x + 1L)
  one_to_n <- function(n) as.list(seq_len(n) - 1L)
  ident_if_even <- function(x) if (even(x)) list(x) else list()
  record <- function(L0) function(fn, input, output) L <<- c(L, list(structure(list(input, output), class = "pytuple")))
  bad <- function(expr) stop("ValueError")
  suppressWarnings(try({ brl <- branch_exhaust(branch5) }, silent = TRUE))
  expect_identical(.cgb(py_list(brl(3L))), "[0]")
  expect_identical(.cgb(py_list(brl(7L))), "[10]")
  expect_identical(.cgb(py_list(brl(5L))), "[10, 0]")
})

test_that("branch_test_core: test_debug", {
  posdec <- function(x) if (x > 0) list(x - 1L) else list(x)
  branch5 <- function(x) {
    if (0 < x && x < 5) list(x - 1L)
    else if (5 < x && x < 10) list(x + 1L)
    else if (x == 5) list(x + 1L, x - 1L)
    else list(x)
  }
  even <- function(x) x %% 2L == 0L
  inc <- function(x) list(x + 1L)
  one_to_n <- function(n) as.list(seq_len(n) - 1L)
  ident_if_even <- function(x) if (even(x)) list(x) else list()
  record <- function(L0) function(fn, input, output) L <<- c(L, list(structure(list(input, output), class = "pytuple")))
  bad <- function(expr) stop("ValueError")
  suppressWarnings(try({ file <- StringIO() }, silent = TRUE))
  suppressWarnings(try({ rl <- branch_debug(posdec, file$con) }, silent = TRUE))
  expect_identical(.cgb(py_list(rl(5L))), "[4]")
  expect_identical(.cgb(py_method(file, "getvalue")), "Rule: posdec\\nIn: 5\\nOut: 4\\n\\n")
  for (.o in mget(ls(), inherits = FALSE)) if (inherits(.o, "py_stringio")) try(close(.o$con), silent = TRUE)
})

test_that("branch_test_core: test_multiplex", {
  posdec <- function(x) if (x > 0) list(x - 1L) else list(x)
  branch5 <- function(x) {
    if (0 < x && x < 5) list(x - 1L)
    else if (5 < x && x < 10) list(x + 1L)
    else if (x == 5) list(x + 1L, x - 1L)
    else list(x)
  }
  even <- function(x) x %% 2L == 0L
  inc <- function(x) list(x + 1L)
  one_to_n <- function(n) as.list(seq_len(n) - 1L)
  ident_if_even <- function(x) if (even(x)) list(x) else list()
  record <- function(L0) function(fn, input, output) L <<- c(L, list(structure(list(input, output), class = "pytuple")))
  bad <- function(expr) stop("ValueError")
  suppressWarnings(try({ brl <- multiplex(posdec, branch5) }, silent = TRUE))
  expect_identical(.cgb(py_list(brl(3L))), "[2]")
  expect_identical(.cgb(py_list(brl(7L))), "[6, 8]")
  expect_identical(.cgb(py_list(brl(5L))), "[4, 6]")
})

test_that("branch_test_core: test_condition", {
  posdec <- function(x) if (x > 0) list(x - 1L) else list(x)
  branch5 <- function(x) {
    if (0 < x && x < 5) list(x - 1L)
    else if (5 < x && x < 10) list(x + 1L)
    else if (x == 5) list(x + 1L, x - 1L)
    else list(x)
  }
  even <- function(x) x %% 2L == 0L
  inc <- function(x) list(x + 1L)
  one_to_n <- function(n) as.list(seq_len(n) - 1L)
  ident_if_even <- function(x) if (even(x)) list(x) else list()
  record <- function(L0) function(fn, input, output) L <<- c(L, list(structure(list(input, output), class = "pytuple")))
  bad <- function(expr) stop("ValueError")
  suppressWarnings(try({ brl <- branch_condition(even, branch5) }, silent = TRUE))
  expect_identical(.cgb(py_list(brl(4L))), "[3]")
  expect_identical(.cgb(py_list(branch5(4L))), "[3]")
  expect_identical(.cgb(py_list(brl(5L))), "[]")
})

test_that("branch_test_core: test_sfilter", {
  posdec <- function(x) if (x > 0) list(x - 1L) else list(x)
  branch5 <- function(x) {
    if (0 < x && x < 5) list(x - 1L)
    else if (5 < x && x < 10) list(x + 1L)
    else if (x == 5) list(x + 1L, x - 1L)
    else list(x)
  }
  even <- function(x) x %% 2L == 0L
  inc <- function(x) list(x + 1L)
  one_to_n <- function(n) as.list(seq_len(n) - 1L)
  ident_if_even <- function(x) if (even(x)) list(x) else list()
  record <- function(L0) function(fn, input, output) L <<- c(L, list(structure(list(input, output), class = "pytuple")))
  bad <- function(expr) stop("ValueError")
  suppressWarnings(try({ brl <- sfilter(even, one_to_n) }, silent = TRUE))
  expect_identical(.cgb(py_list(brl(10L))), "[0, 2, 4, 6, 8]")
})

test_that("branch_test_core: test_notempty", {
  posdec <- function(x) if (x > 0) list(x - 1L) else list(x)
  branch5 <- function(x) {
    if (0 < x && x < 5) list(x - 1L)
    else if (5 < x && x < 10) list(x + 1L)
    else if (x == 5) list(x + 1L, x - 1L)
    else list(x)
  }
  even <- function(x) x %% 2L == 0L
  inc <- function(x) list(x + 1L)
  one_to_n <- function(n) as.list(seq_len(n) - 1L)
  ident_if_even <- function(x) if (even(x)) list(x) else list()
  record <- function(L0) function(fn, input, output) L <<- c(L, list(structure(list(input, output), class = "pytuple")))
  bad <- function(expr) stop("ValueError")
  suppressWarnings(try({ brl <- notempty(ident_if_even) }, silent = TRUE))
  expect_identical(.cgb(py_list(brl(4L))), "[4]")
  expect_identical(.cgb(py_list(brl(5L))), "[5]")
})

test_that("branch_test_core: test_chain", {
  posdec <- function(x) if (x > 0) list(x - 1L) else list(x)
  branch5 <- function(x) {
    if (0 < x && x < 5) list(x - 1L)
    else if (5 < x && x < 10) list(x + 1L)
    else if (x == 5) list(x + 1L, x - 1L)
    else list(x)
  }
  even <- function(x) x %% 2L == 0L
  inc <- function(x) list(x + 1L)
  one_to_n <- function(n) as.list(seq_len(n) - 1L)
  ident_if_even <- function(x) if (even(x)) list(x) else list()
  record <- function(L0) function(fn, input, output) L <<- c(L, list(structure(list(input, output), class = "pytuple")))
  bad <- function(expr) stop("ValueError")
  expect_identical(.cgb(py_list(py_call(branch_chain(), 2L))), "[2]")
  expect_identical(.cgb(py_list(py_call(branch_chain(inc, inc), 2L))), "[4]")
  expect_identical(.cgb(py_list(py_call(branch_chain(branch5, inc), 4L))), "[4]")
  expect_identical(.cgb(py_list(py_call(branch_chain(branch5, inc), 5L))), "[7, 5]")
  expect_identical(.cgb(py_list(py_call(branch_chain(inc, branch5), 5L))), "[7]")
})

test_that("branch_test_core: test_onaction", {
  posdec <- function(x) if (x > 0) list(x - 1L) else list(x)
  branch5 <- function(x) {
    if (0 < x && x < 5) list(x - 1L)
    else if (5 < x && x < 10) list(x + 1L)
    else if (x == 5) list(x + 1L, x - 1L)
    else list(x)
  }
  even <- function(x) x %% 2L == 0L
  inc <- function(x) list(x + 1L)
  one_to_n <- function(n) as.list(seq_len(n) - 1L)
  ident_if_even <- function(x) if (even(x)) list(x) else list()
  record <- function(L0) function(fn, input, output) L <<- c(L, list(structure(list(input, output), class = "pytuple")))
  bad <- function(expr) stop("ValueError")
  suppressWarnings(try({ L <- list() }, silent = TRUE))
  expect_identical(.cgb(py_list(py_call(onaction(inc, record(L)), 2L))), "[3]")
  expect_identical(.cgb(L), "[(2, 3)]")
  expect_identical(.cgb(py_list(py_call(onaction(branch_identity, record(L)), 2L))), "[2]")
  expect_identical(.cgb(L), "[(2, 3)]")
})

test_that("branch_test_core: test_yieldify", {
  posdec <- function(x) if (x > 0) list(x - 1L) else list(x)
  branch5 <- function(x) {
    if (0 < x && x < 5) list(x - 1L)
    else if (5 < x && x < 10) list(x + 1L)
    else if (x == 5) list(x + 1L, x - 1L)
    else list(x)
  }
  even <- function(x) x %% 2L == 0L
  inc <- function(x) list(x + 1L)
  one_to_n <- function(n) as.list(seq_len(n) - 1L)
  ident_if_even <- function(x) if (even(x)) list(x) else list()
  record <- function(L0) function(fn, input, output) L <<- c(L, list(structure(list(input, output), class = "pytuple")))
  bad <- function(expr) stop("ValueError")
  suppressWarnings(try({ yinc <- yieldify(function(x) x + 1L) }, silent = TRUE))
  expect_identical(.cgb(py_list(yinc(3L))), "[4]")
})

test_that("branch_test_core: test_do_one", {
  posdec <- function(x) if (x > 0) list(x - 1L) else list(x)
  branch5 <- function(x) {
    if (0 < x && x < 5) list(x - 1L)
    else if (5 < x && x < 10) list(x + 1L)
    else if (x == 5) list(x + 1L, x - 1L)
    else list(x)
  }
  even <- function(x) x %% 2L == 0L
  inc <- function(x) list(x + 1L)
  one_to_n <- function(n) as.list(seq_len(n) - 1L)
  ident_if_even <- function(x) if (even(x)) list(x) else list()
  record <- function(L0) function(fn, input, output) L <<- c(L, list(structure(list(input, output), class = "pytuple")))
  bad <- function(expr) stop("ValueError")
  expect_identical(.cgb(py_list(py_call(branch_do_one(inc), 3L))), "[4]")
  expect_identical(.cgb(py_list(py_call(branch_do_one(inc, bad), 3L))), "[4]")
  expect_identical(.cgb(py_list(py_call(branch_do_one(inc, posdec), 3L))), "[4]")
})

test_that("branch_test_core: branch_core_extra", {
  posdec <- function(x) if (x > 0) list(x - 1L) else list(x)
  branch5 <- function(x) {
    if (0 < x && x < 5) list(x - 1L)
    else if (5 < x && x < 10) list(x + 1L)
    else if (x == 5) list(x + 1L, x - 1L)
    else list(x)
  }
  even <- function(x) x %% 2L == 0L
  inc <- function(x) list(x + 1L)
  one_to_n <- function(n) as.list(seq_len(n) - 1L)
  ident_if_even <- function(x) if (even(x)) list(x) else list()
  record <- function(L0) function(fn, input, output) L <<- c(L, list(structure(list(input, output), class = "pytuple")))
  bad <- function(expr) stop("ValueError")
  expect_identical(.cgb(py_list(branch_identity(7L))), "[7]")
  expect_identical(.cgb(py_list(py_call(branch_exhaust(posdec), 4L))), "[0]")
  expect_identical(.cgb(py_list(py_call(branch_exhaust(branch_identity), 4L))), "[4]")
  expect_identical(.cgb(py_list(py_call(multiplex(branch5, posdec), 5L))), "[6, 4]")
  expect_identical(.cgb(py_list(py_call(multiplex(), 5L))), "[]")
  expect_identical(.cgb(py_list(py_call(branch_do_one(ident_if_even, posdec), 5L))), "[4]")
  expect_identical(.cgb(py_list(py_call(branch_do_one(ident_if_even, posdec), 4L))), "[4]")
  expect_identical(.cgb(py_list(py_call(branch_do_one(), 4L))), "[]")
  expect_identical(.cgb(py_list(py_call(branch_chain(posdec, posdec, posdec), 5L))), "[2]")
  expect_identical(.cgb(py_list(py_call(branch_chain(branch5, branch5), 5L))), "[7, 3]")
  expect_identical(.cgb(py_list(py_call(sfilter(even, branch5), 5L))), "[6, 4]")
  expect_identical(.cgb(py_list(py_call(notempty(sfilter(even, branch5)), 3L))), "[2]")
  expect_identical(.cgb(py_list(py_call(branch_condition(even, inc), 2L))), "[3]")
  expect_identical(.cgb(py_list(py_call(branch_condition(even, inc), 3L))), "[]")
})

test_that("branch_test_tools: test_zero_ints", {
  Basic <- function(...) new_sym("Basic", list(...))
  .isint <- function(x) is_sym(x) && x$type == "Integer"
  .iv <- function(x) as.numeric(x$value)
  posdec <- function(x) if (.isint(x) && .iv(x) > 0) list(x - 1L) else list(x)
  branch5 <- function(x) {
    if (!.isint(x)) return(list())
    v <- .iv(x)
    if (0 < v && v < 5) list(x - 1L)
    else if (5 < v && v < 10) list(x + 1L)
    else if (v == 5) list(x + 1L, x - 1L)
    else list(x)
  }
  suppressWarnings(try({ expr <- Basic(S(2L), Basic(S(5L), S(3L)), S(8L)) }, silent = TRUE))
  suppressWarnings(try({ brl <- branch_canon(posdec) }, silent = TRUE))
  expect_identical(.cgb(py_list(brl(expr))), "[Basic(0, Basic(0, 0), 0)]")
})

test_that("branch_test_tools: test_split5", {
  Basic <- function(...) new_sym("Basic", list(...))
  .isint <- function(x) is_sym(x) && x$type == "Integer"
  .iv <- function(x) as.numeric(x$value)
  posdec <- function(x) if (.isint(x) && .iv(x) > 0) list(x - 1L) else list(x)
  branch5 <- function(x) {
    if (!.isint(x)) return(list())
    v <- .iv(x)
    if (0 < v && v < 5) list(x - 1L)
    else if (5 < v && v < 10) list(x + 1L)
    else if (v == 5) list(x + 1L, x - 1L)
    else list(x)
  }
  suppressWarnings(try({ expr <- Basic(S(2L), Basic(S(5L), S(3L)), S(8L)) }, silent = TRUE))
  suppressWarnings(try({ brl <- branch_canon(branch5) }, silent = TRUE))
  expect_identical(.cgb(py_list(brl(expr))), "[Basic(0, Basic(10, 0), 10), Basic(0, Basic(0, 0), 10)]")
})

test_that("branch_test_traverse: test_top_down_easy", {
  Basic <- function(...) new_sym("Basic", list(...))
  inc <- function(x) if (is_sym(x) && x$type == "Integer") list(x + 1L) else list()
  split5 <- function(x) if (is_sym(x) && x$type == "Integer" && as.numeric(x$value) == 5) list(x - 1L, x + 1L) else list()
  suppressWarnings(try({ expr <- Basic(S(1L), S(2L)) }, silent = TRUE))
  suppressWarnings(try({ brl <- branch_top_down(inc) }, silent = TRUE))
  expect_identical(.cgb(py_list(brl(expr))), "[Basic(2, 3)]")
})

test_that("branch_test_traverse: test_top_down_big_tree", {
  Basic <- function(...) new_sym("Basic", list(...))
  inc <- function(x) if (is_sym(x) && x$type == "Integer") list(x + 1L) else list()
  split5 <- function(x) if (is_sym(x) && x$type == "Integer" && as.numeric(x$value) == 5) list(x - 1L, x + 1L) else list()
  suppressWarnings(try({ expr <- Basic(S(1L), Basic(S(2L)), Basic(S(3L), Basic(S(4L)), S(5L))) }, silent = TRUE))
  suppressWarnings(try({ brl <- branch_top_down(inc) }, silent = TRUE))
  expect_identical(.cgb(py_list(brl(expr))), "[Basic(2, Basic(3), Basic(4, Basic(5), 6))]")
})

test_that("branch_test_traverse: test_top_down_harder_function", {
  Basic <- function(...) new_sym("Basic", list(...))
  inc <- function(x) if (is_sym(x) && x$type == "Integer") list(x + 1L) else list()
  split5 <- function(x) if (is_sym(x) && x$type == "Integer" && as.numeric(x$value) == 5) list(x - 1L, x + 1L) else list()
  suppressWarnings(try({ expr <- Basic(Basic(S(5L), S(6L)), S(1L)) }, silent = TRUE))
  suppressWarnings(try({ brl <- branch_top_down(split5) }, silent = TRUE))
  expect_identical(.cgb(py_list(brl(expr))), "[Basic(Basic(4, 6), 1), Basic(Basic(6, 6), 1)]")
})

test_that("branch_test_traverse: test_sall", {
  Basic <- function(...) new_sym("Basic", list(...))
  inc <- function(x) if (is_sym(x) && x$type == "Integer") list(x + 1L) else list()
  split5 <- function(x) if (is_sym(x) && x$type == "Integer" && as.numeric(x$value) == 5) list(x - 1L, x + 1L) else list()
  suppressWarnings(try({ expr <- Basic(S(1L), S(2L)) }, silent = TRUE))
  suppressWarnings(try({ brl <- branch_sall(inc) }, silent = TRUE))
  expect_identical(.cgb(py_list(brl(expr))), "[Basic(2, 3)]")
  suppressWarnings(try({ expr <- Basic(S(1L), S(2L), Basic(S(3L), S(4L))) }, silent = TRUE))
  suppressWarnings(try({ brl <- branch_sall(branch_do_one(inc, branch_identity)) }, silent = TRUE))
  expect_identical(.cgb(py_list(brl(expr))), "[Basic(2, 3, Basic(3, 4))]")
})

test_that("branch_test_traverse: branch_traverse_extra", {
  Basic <- function(...) new_sym("Basic", list(...))
  inc <- function(x) if (is_sym(x) && x$type == "Integer") list(x + 1L) else list()
  split5 <- function(x) if (is_sym(x) && x$type == "Integer" && as.numeric(x$value) == 5) list(x - 1L, x + 1L) else list()
  expect_identical(.cgb(py_list(py_call(branch_top_down(split5), Basic(S(5L), S(5L))))), "[Basic(4, 4), Basic(4, 6), Basic(6, 4), Basic(6, 6)]")
  expect_identical(.cgb(py_list(py_call(branch_sall(split5), Basic(S(5L), S(5L))))), "[Basic(4, 4), Basic(4, 6), Basic(6, 4), Basic(6, 6)]")
  expect_identical(.cgb(py_list(py_call(branch_sall(inc), S(3L)))), "[3]")
  expect_identical(.cgb(py_list(py_call(branch_top_down(branch_identity), Basic(S(1L), Basic(S(2L)))))), "[Basic(1, Basic(2))]")
})

