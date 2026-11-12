skip_unless_full_tests()
# Expected values are the outputs of SymPy for the same calls (sympy/physics/vector/tests).
# The calls are SymPy test code translated to R; the helpers below map
# Python syntax onto the package API.

ns <- asNamespace("symcompR")

# Shims mapping Python syntax used by the tests onto the package API. They
# only route calls; all mathematics is done by package functions.
H <- new.env(parent = ns)
local(envir = H, {
  `%%` <- function(e1, e2) {
    if (is.character(e1)) {
      args <- if (inherits(e2, "pytuple") || (is.list(e2) && !is_sym(e2) && !is.object(e2))) unclass(e2) else list(e2)
      args <- lapply(args, function(z) if (is.character(z)) z else .pv_any_str(z))
      return(do.call(sprintf, c(list(gsub("%r", "%s", e1, fixed = TRUE)), args)))
    }
    base::`%%`(e1, e2)
  }
  `+` <- function(e1, e2) {
    if (missing(e2)) return(base::`+`(e1))
    if (is.character(e1) && is.character(e2)) return(paste0(e1, e2))
    base::`+`(e1, e2)
  }
  is_phys <- function(x) inherits(x, "physics_obj")
  tup <- function(...) .pytuple(list(...))
  pydict <- function(keys, values) .pydict(keys, values)
  pyset <- function(x) structure(pylist(x), pyset = TRUE)
  pylist <- function(x) {
    if (is.null(x)) return(list())
    if (inherits(x, "pydict")) return(x$keys)
    if (inherits(x, "symmatrix")) return(x$data)
    if (is_sym(x) && x$type == "Tuple") return(x$args)
    if (is.list(x) && !is_sym(x) && !is_phys(x)) { a <- attributes(x); y <- x; attributes(y) <- NULL; return(y) }
    if (is.atomic(x) && !is.character(x)) return(as.list(x))
    if (is.character(x) && length(x) == 1L) return(as.list(strsplit(x, "")[[1L]]))
    if (is.character(x)) return(as.list(x))
    if (inherits(x, "pynamedtuple")) { y <- x; attributes(y) <- NULL; return(y) }
    if (inherits(x, "ReferenceFrame")) return(.pv_frame_iter(x))
    stop("cannot iterate")
  }
  pytup <- function(x = list()) .pytuple(pylist(x))
  pyiter <- function(x) pylist(x)
  pyrange <- function(a, b = NULL, s = 1) { if (is.null(b)) { b <- a; a <- 0 }; if (b <= a) list() else as.list(seq(a, b - 1, by = s)) }
  pyzip <- function(...) { xs <- lapply(list(...), pylist); n <- min(vapply(xs, length, 0L)); lapply(seq_len(n), function(i) .pytuple(lapply(xs, `[[`, i))) }
  pyenumerate <- function(x) { x <- pylist(x); lapply(seq_along(x), function(i) .pytuple(list(i - 1, x[[i]]))) }
  pysum <- function(x, start = 0) Reduce(function(a, b) a + b, pylist(x), start)
  pyall <- function(x) all(vapply(pylist(x), pybool, TRUE))
  pyany <- function(x) any(vapply(pylist(x), pybool, TRUE))
  pybool <- function(x) {
    if (is.logical(x)) return(isTRUE(x))
    if (is.numeric(x)) return(x != 0)
    if (is_sym(x)) { if (x$type == "BooleanAtom") return(x$name == "true"); return(!is_zero_num(x)) }
    if (inherits(x, c("PhysVector", "PhysDyadic"))) return(length(.subset2(x, "args")) > 0)
    if (is.list(x)) return(length(x) > 0)
    TRUE
  }
  pyabs <- function(x) if (is.numeric(x)) abs(x) else fcall("Abs", list(as_sym(x)))
  pymax <- function(...) { a <- list(...); if (length(a) == 1L) a <- pylist(a[[1L]]); Reduce(function(x, y) if (pycmp(">", y, x)) y else x, a) }
  pymin <- function(...) { a <- list(...); if (length(a) == 1L) a <- pylist(a[[1L]]); Reduce(function(x, y) if (pycmp("<", y, x)) y else x, a) }
  pyvsimp <- function() .pv_state$simp
  pyeq <- function(a, b) {
    if (is_phys(a) && !is.null(.pv_lookup(a, ".pvM_", "__eq__"))) return(a$`__eq__`(b))
    if (inherits(a, c("PhysVector", "PhysDyadic")) || inherits(b, c("PhysVector", "PhysDyadic"))) return(a == b)
    if (is.environment(a) || is.environment(b)) return(identical(a, b))
    if (inherits(a, "symmatrix") || inherits(b, "symmatrix"))
      return(inherits(a, "symmatrix") && inherits(b, "symmatrix") && .pv_meq(a, b))
    if (is.null(a) || is.null(b)) return(is.null(a) && is.null(b))
    if (is.character(a) || is.character(b)) return(identical(a, b))
    if (is.logical(a) && is.logical(b)) return(identical(a, b))
    if ((is_sym(a) || is.numeric(a) || is.logical(a)) && (is_sym(b) || is.numeric(b) || is.logical(b)))
      return(identical(as_sym(a)$key, as_sym(b)$key))
    if (inherits(a, "pydict") && inherits(b, "pydict"))
      return(length(a$keys) == length(b$keys) && all(vapply(seq_along(a$keys), function(i) {
        j <- which(vapply(b$keys, function(k) pyeq(k, a$keys[[i]]), TRUE))
        length(j) == 1L && pyeq(a$values[[i]], b$values[[j]])
      }, TRUE)))
    if (isTRUE(attr(a, "pyset")) || isTRUE(attr(b, "pyset"))) {
      a <- pylist(a); b <- pylist(b)
      return(all(vapply(a, pyin, TRUE, b = b)) && all(vapply(b, pyin, TRUE, b = a)))
    }
    if (is.list(a) && is.list(b)) return(length(a) == length(b) && all(vapply(seq_along(a), function(i) pyeq(a[[i]], b[[i]]), TRUE)))
    identical(a, b)
  }
  pyin <- function(a, b) {
    if (inherits(b, "pydict")) return(any(vapply(b$keys, function(k) pyeq(k, a), TRUE)))
    if (is.character(b) && length(b) == 1L) return(grepl(a, b, fixed = TRUE))
    any(vapply(pylist(b), function(k) pyeq(k, a), TRUE))
  }
  pyis <- function(a, b) {
    if (is.null(a) || is.null(b)) return(is.null(a) && is.null(b))
    if (is.logical(a) || is.logical(b)) return(identical(a, b))
    if (is_sym(a) && is_sym(b)) return(identical(a$key, b$key))
    identical(a, b)
  }
  pycmp <- function(op, a, b) {
    if (is.numeric(a) && is.numeric(b)) return(do.call(op, list(a, b)))
    r <- do.call(op, list(as_sym(a), as_sym(b)))
    if (is.logical(r)) r else stop("TypeError: cannot determine truth value of Relational")
  }
  pypow <- function(a, b) {
    if (is_phys(a) || is_phys(b)) stop(.sympy_condition("TypeError", "unsupported operand type(s) for ** or pow()"))
    a ^ b
  }
  pyxor <- function(a, b) a ^ b
  pystr <- function(x) .pv_any_str(x)
  pyrepr <- function(x) if (is_phys(x)) .pv_objrepr(x) else .pv_any_str(x, TRUE)
  pytype <- function(x) class(x)[1L]
  pyinplace <- function(x, m) { r <- pycall(x, m); if (inherits(x, "symmatrix")) r else x }
  pyhasattr <- function(x, n) !inherits(tryCatch(pyattr(x, n), error = function(e) e), "error")
  pyisinstance <- function(x, cls) {
    e <- substitute(cls)
    if (is.call(e) && identical(e[[1L]], as.name("tup")))
      return(any(vapply(as.list(e)[-1L], function(a) .isinst(x, deparse(a)), TRUE)))
    v <- tryCatch(cls, error = function(err) NULL)
    .isinst(x, if (is.character(v)) v else deparse(e))
  }
  .isinst <- function(x, nm) {
    switch(nm,
           pytup = inherits(x, c("pytuple", "pynamedtuple")), pylist = is.list(x) && !is.object(x),
           BodyBase = inherits(x, "BodyBase"),
           PhysVector = inherits(x, "PhysVector"), PhysDyadic = inherits(x, "PhysDyadic"),
           ReferenceFrame = inherits(x, "ReferenceFrame"), PhysPoint = inherits(x, "PhysPoint"),
           CoordinateSym = .pv_is_coordsym(x), Symbol = is_sym(x) && x$type == "Symbol",
           Matrix = , ImmutableMatrix = , MatrixBase = inherits(x, "symmatrix"),
           Expr = , ExprType = is_sym(x),
           if (is_sym(x) && x$type == "Function") identical(x$name, nm) else inherits(x, nm))
  }
  pysymbols <- function(names, cls = NULL, ...) {
    nms <- .expand_symbol_names(names)
    mk <- if (is.null(cls)) function(n) Symbol(n, ...) else if (identical(cls, Function)) function(n) Function(n) else function(n) cls(n)
    res <- lapply(nms, mk)
    if (length(res) == 1L && !endsWith(trimws(names), ",")) return(res[[1L]])
    res
  }
  pyzeros <- function(n, m = n) zeros(n, m)
  pyeye <- function(n) eye(n)
  pydiff <- function(e, ...) {
    if (is_phys(e)) return(pycall(e, "diff", ...))
    a <- list(...)
    if (inherits(e, "symmatrix")) return(.pv_mmap(e, function(z) do.call(pydiff, c(list(z), a))))
    do.call(diff, c(list(as_sym(e)), a))
  }
  pysimplify <- function(e, ...) {
    if (is_phys(e)) return(pycall(e, "simplify"))
    if (inherits(e, "symmatrix")) return(.pv_mmap(e, simplify))
    simplify(as_sym(e), ...)
  }
  pyexpand <- function(e, ...) {
    if (inherits(e, "symmatrix")) return(.pv_mmap(e, function(z) expand(z, ...)))
    expand(as_sym(e), ...)
  }
  pytrigsimp <- function(e, ...) {
    if (inherits(e, "symmatrix")) return(.pv_mmap(e, function(z) trigsimp(z, ...)))
    trigsimp(as_sym(e), ...)
  }
  pycancel <- function(e, ...) cancel(as_sym(e), ...)
  pyintegrate <- function(f, ...) integrate(f, ...)
  pysolve <- function(...) solve(...)
  pyMul <- function(...) .mul(lapply(list(...), as_sym))
  pyAdd <- function(...) .add(lapply(list(...), as_sym))
  pyPow <- function(a, b) .pow(as_sym(a), as_sym(b))
  pylinear_eq_to_matrix <- function(eqs, syms) .pytuple(.pv_linear_eq_to_matrix(pylist(eqs), pylist(syms)))
  pysorted <- function(x, key = NULL, ...) {
    x <- pylist(x)
    if (!is.null(key)) { ks <- lapply(x, key); return(x[order(vapply(ks, function(k) if (is.numeric(k)) k else NA_real_, 0))]) }
    x[order(vapply(x, .pv_any_str, ""))]
  }
  pyordered <- function(x, ...) { x <- pylist(x); if (all(vapply(x, is_sym, TRUE))) .bool_print_order(x) else x }
  pydictf <- function(x = NULL) if (is.null(x)) .pydict(list(), list()) else x
  pydocall <- function(f, args, star) do.call(f, c(args, pylist(star)))
  pyunbound <- function(name) function(self, ...) pycall(self, name, ...)
  pyclassmethod <- function(cls, meth) {
    f <- get0(paste0(cls, "_", meth), envir = H)
    if (is.function(f)) return(f)
    # members of an enumeration, such as MusculotendonFormulation.RIGID_TENDON
    g <- get0(cls, envir = H)
    if (is.function(g) && grepl("^[A-Z][A-Z_]*$", meth)) return(g(meth))
    function(self, ...) pycall(self, meth, ...)
  }
  pysetattr <- function(obj, name, val) {
    if (inherits(obj, "pyself")) { assign(name, val, envir = obj); return(invisible(NULL)) }
    if (is_phys(obj)) { eval(call("$<-", obj, name, val)); return(invisible(NULL)) }
    if (identical(name, "simp")) { .pv_state$simp <- isTRUE(val); return(invisible(NULL)) }
    if (identical(name, "_t")) { .pv_state$t <- val; return(invisible(NULL)) }
    stop("cannot set attribute ", name)
  }
  pyget <- function(x, i, j = NULL) {
    if (!is.null(j) || inherits(i, "pysl")) {
      if (inherits(x, "symmatrix")) {
        ri <- if (inherits(i, "pysl")) .sl_range(i, x$nrow) else i + 1L
        ci <- if (is.null(j)) seq_len(x$ncol) else if (inherits(j, "pysl")) .sl_range(j, x$ncol) else j + 1L
        if (length(ri) == 1L && length(ci) == 1L && !inherits(i, "pysl") && !inherits(j, "pysl")) return(.mget(x, ri, ci))
        return(.new_mat(unlist(lapply(ri, function(r) lapply(ci, function(c) .mget(x, r, c))), recursive = FALSE), length(ri), length(ci)))
      }
      stop("2-d index")
    }
    if (inherits(x, "ReferenceFrame")) return(if (is.character(i)) x[i] else x[i + 1L])
    if (inherits(x, "pydict")) {
      k <- which(vapply(x$keys, function(z) pyeq(z, i), TRUE))
      if (!length(k)) stop(.sympy_condition("KeyError", .pv_any_str(i, TRUE)))
      return(x$values[[k[1L]]])
    }
    if (inherits(x, "symmatrix")) {
      if (i < 0) i <- length(x$data) + i
      if (x$ncol == 1L || x$nrow == 1L || TRUE) return(x$data[[i + 1L]])
    }
    if (is_sym(x)) { a <- .sympy_args(x); if (i < 0) i <- length(a) + i; return(a[[i + 1L]]) }
    if (is.character(x) && length(x) == 1L) return(substr(x, i + 1L, i + 1L))
    if (is.character(x)) return(x[[i + 1L]])
    if (is.list(x)) { if (i < 0) i <- length(x) + i; if (i >= length(x)) stop(.sympy_condition("IndexError", "list index out of range")); return(x[[i + 1L]]) }
    stop("cannot index")
  }
  pysl <- function(lo, hi) structure(list(lo = lo, hi = hi), class = "pysl")
  .sl_range <- function(s, n) {
    lo <- if (is.null(s$lo)) 0 else if (s$lo < 0) n + s$lo else s$lo
    hi <- if (is.null(s$hi)) n else if (s$hi < 0) n + s$hi else min(s$hi, n)
    if (hi <= lo) integer() else (lo + 1L):hi
  }
  pyslice <- function(x, lo, hi, st) {
    n <- if (inherits(x, "symmatrix")) length(x$data) else length(x)
    r <- .sl_range(list(lo = lo, hi = hi), n)
    if (!is.null(st)) r <- r[seq(1L, length(r), by = st)]
    if (inherits(x, "symmatrix")) return(x$data[r])
    if (is.character(x) && length(x) == 1L) return(paste(strsplit(x, "")[[1L]][r], collapse = ""))
    y <- x[r]
    if (inherits(x, "pytuple")) class(y) <- "pytuple"
    y
  }
  # Attribute access: physics objects use `$`; others map to package functions
  pyself <- function() structure(new.env(), class = "pyself")
  pyattr <- function(x, name) {
    if (inherits(x, "pyself")) {
      if (!exists(name, envir = x, inherits = FALSE))
        stop(.sympy_condition("AttributeError", paste0("'S' object has no attribute '", name, "'")))
      return(get(name, envir = x))
    }
    if (is_phys(x)) return(eval(call("$", x, name)))
    if (inherits(x, "symmatrix")) {
      return(switch(name, T = transpose(x), shape = .pytuple(list(x$nrow, x$ncol)), rows = x$nrow, cols = x$ncol,
                    free_symbols = structure(.pv_unique_syms(x), pyset = TRUE),
                    is_zero_matrix = all(vapply(x$data, is_zero_num, TRUE)),
                    stop(.sympy_condition("AttributeError", paste0("Matrix has no attribute ", name)))))
    }
    if (is_sym(x)) {
      return(switch(name,
        args = .pytuple(.sympy_args(x)), func = x$type, free_symbols = structure(free_symbols(x), pyset = TRUE),
        name = x$name, expr = x$args[[1L]], is_zero = { r <- is_zero(x); if (is.na(r)) NULL else r },
        is_Add = x$type == "Add", is_Mul = x$type == "Mul", is_Pow = x$type == "Pow",
        is_Symbol = x$type == "Symbol", is_number = .expr_is_number(x), lhs = x$args[[1L]], rhs = x$args[[2L]],
        variables = .pytuple(.pv_deriv_vars(x)), frame = x$extra$frame,
        if (startsWith(name, "is_")) { r <- ask_fact(x, substring(name, 4L)); if (is.na(r)) NULL else r } else
        stop(.sympy_condition("AttributeError", paste0("'", x$type, "' object has no attribute '", name, "'")))))
    }
    if (inherits(x, "pydict")) stop("dict attribute ", name)
    stop(.sympy_condition("AttributeError", paste0("object has no attribute '", name, "'")))
  }
  pycall <- function(x, name, ...) {
    a <- list(...)
    if (inherits(x, "pyprinter")) return(do.call(x[[name]], a))
    if (is.function(x) && name == "with_defaults") return(ns$.bm_with_defaults(x, a[[1L]]))
    if (is_phys(x)) return(do.call(eval(call("$", x, name)), a))
    if (is.character(x) && length(x) == 1L) {
      return(switch(name, lower = tolower(x), upper = toupper(x), startswith = startsWith(x, a[[1L]]),
                    endswith = endsWith(x, a[[1L]]), strip = trimws(x), split = as.list(strsplit(x, if (length(a)) a[[1L]] else " +")[[1L]]),
                    replace = gsub(a[[1L]], a[[2L]], x, fixed = TRUE), join = paste(unlist(pylist(a[[1L]])), collapse = x),
                    stop(.sympy_condition("AttributeError", paste0("'str' object has no attribute '", name, "'")))))
    }
    if (inherits(x, "pydict")) {
      return(switch(name, keys = x$keys, copy = x, values = x$values, items = Map(function(k, v) .pytuple(list(k, v)), x$keys, x$values),
                    get = { r <- tryCatch(pyget(x, a[[1L]]), error = function(e) if (length(a) > 1L) a[[2L]] else NULL); r },
                    stop("dict method ", name)))
    }
    if (inherits(x, "symmatrix")) {
      return(switch(name,
        subs = .pv_mmap(x, function(z) do.call(.pv_subs, c(list(z), a))),
        xreplace = .pv_mmap(x, function(z) .pv_xreplace(z, a[[1L]])),
        diff = do.call(pydiff, c(list(x), a)), simplify = .pv_mmap(x, simplify), expand = .pv_mmap(x, expand),
        doit = .pv_mmap(x, doit), applyfunc = .pv_mmap(x, function(z) as_sym(a[[1L]](z))),
        evalf = , n = .pv_mmap(x, function(z) N(z, if (length(a)) a[[1L]] else 15)),
        trigsimp = .pv_mmap(x, function(z) do.call(trigsimp, c(list(z), a))),
        jacobian = jacobian(x, if (inherits(a[[1L]], "symmatrix")) a[[1L]] else Matrix(pylist(a[[1L]]))),
        reshape = .new_mat(x$data, a[[1L]], a[[2L]]), transpose = transpose(x), tolist = .mrows(x),
        dot = .add(Map(function(p, q) .mul(list(p, q)), x$data, Matrix(a[[1L]])$data)),
        inv = inverse(x), det = det(x), LUsolve = .pv_mat_lusolve(x, a[[1L]]), solve = solve(x, a[[1L]]),
        col = .new_mat(lapply(seq_len(x$nrow), function(i) .mget(x, i, a[[1L]] + 1L)), x$nrow, 1L),
        row = .new_mat(lapply(seq_len(x$ncol), function(j) .mget(x, a[[1L]] + 1L, j)), 1L, x$ncol),
        has = any(vapply(x$data, function(z) .has1(z, as_sym(a[[1L]])), TRUE)),
        atoms = structure(Reduce(function(p, q) c(p, q), lapply(x$data, function(z) do.call(atoms, c(list(z), a))), list()), pyset = TRUE),
        stop(.sympy_condition("AttributeError", paste0("Matrix has no attribute '", name, "'")))))
    }
    if (is_sym(x) || is.numeric(x)) {
      x <- as_sym(x)
      return(switch(name,
        subs = do.call(.pv_subs, c(list(x), a)), xreplace = .pv_xreplace(x, a[[1L]]),
        diff = do.call(diff, c(list(x), a)), simplify = simplify(x), expand = do.call(expand, c(list(x), a)),
        doit = do.call(doit, c(list(x), a)), inverse = inverse(x),
        evalf = , n = N(x, if (length(a)) a[[1L]] else 15), trigsimp = do.call(trigsimp, c(list(x), a)),
        factor = factor_sym(x), together = together(x), cancel = cancel(x), rewrite = do.call(rewrite, c(list(x), a)),
        has = any(vapply(a, function(z) .has1(x, as_sym(z)), TRUE)),
        atoms = structure(do.call(atoms, c(list(x), a)), pyset = TRUE),
        as_ordered_terms = .add_print_order(if (x$type == "Add") x$args else list(x)),
        coeff = coeff(x, a[[1L]]), is_constant = is_constant(x),
        stop(.sympy_condition("AttributeError", paste0("'", x$type, "' object has no attribute '", name, "'")))))
    }
    if (is.list(x) && name == "index") return(which(vapply(x, function(z) pyeq(z, a[[1L]]), TRUE))[1L] - 1)
    if (is.list(x) && name == "issubset") return(all(vapply(x, function(z) any(vapply(a[[1L]], function(w) pyeq(z, w), TRUE)), TRUE)))
    if (is.list(x) && name == "count") return(sum(vapply(x, function(z) pyeq(z, a[[1L]]), TRUE)))
    stop(.sympy_condition("AttributeError", paste0("object has no attribute '", name, "'")))
  }
})
# printer classes mapped onto the package printers
local(envir = H, {
  .printer <- function(f) { force(f); function(...) structure(list(doprint = f), class = "pyprinter") }
  C89CodePrinter <- C99CodePrinter <- C11CodePrinter <- .printer(function(e) ccode(e))
  CXX98CodePrinter <- CXX11CodePrinter <- CXX17CodePrinter <- .printer(function(e) cxxcode(e))
  FCodePrinter <- .printer(function(e) fcode(e))
  OctaveCodePrinter <- .printer(function(e) octave_code(e))
  PythonCodePrinter <- .printer(function(e) pycode(e))
  LatexPrinter <- .printer(function(e) latex(e))
})
# private SymPy helpers some tests call directly
assign("_check_dyadic", function(x) ns$.pv_check_dyadic(x), envir = H)
assign("_check_vector", function(x) ns$.pv_check_vector(x), envir = H)
assign("_check_frame", function(x) ns$.pv_check_frame(x), envir = H)
assign("_validate_coordinates", function(coordinates = NULL, speeds = NULL, check_duplicates = TRUE,
                                         is_dynamicsymbols = TRUE, u_auxiliary = NULL)
  ns$.pm_validate_coordinates(coordinates, speeds, check_duplicates, is_dynamicsymbols, u_auxiliary), envir = H)
assign("_parse_linear_solver", function(x) ns$.pm_parse_linear_solver(x), envir = H)
# a lambdified function takes Python floats (Float here) and returns a float
assign("lambdify", function(...) {
  f <- ns$lambdify(...)
  function(...) {
    a <- lapply(list(...), function(z) if (ns$is_sym(z)) ns$num_value_dbl(z) else z)
    r <- do.call(f, a)
    if (is.numeric(r) && length(r) == 1L) ns$Float(r) else r
  }
}, envir = H)
assign("_parse_load", function(x) ns$.pm_parse_load(x), envir = H)
assign("Joint__create_aligned_interframe", function(...) ns$.pm_create_aligned_interframe(...), envir = H)
assign("System_from_newtonian", function(x) ns$System_from_newtonian(x), envir = H)
assign("TorqueActuator_at_pin_joint", function(...) ns$TorqueActuator_at_pin_joint(...), envir = H)
assign("Inertia_from_inertia_scalars", function(...) ns$Inertia_from_inertia_scalars(...), envir = H)
# a matrix LU solve for the shims (only used when tests call M.LUsolve)
assign(".pv_mat_lusolve", function(A, b) solve(A, b), envir = H)

E0 <- new.env(parent = H)
# Python's warnings module: simplefilter("error") turns warnings into errors
.wstate <- new.env()
.wstate$mode <- "default"
assign("warnings", structure(list(simplefilter = function(action, ...) { .wstate$mode <- action; invisible(NULL) }),
                             class = "pyprinter"), envir = H)
.pywarn <- function(expr) withCallingHandlers(expr, warning = function(w) {
  if (identical(.wstate$mode, "error"))
    stop(structure(list(message = conditionMessage(w), call = NULL), class = c("UserWarning", "error", "condition")))
  invokeRestart("muffleWarning")
})

.fmt <- function(v, top = TRUE) {
  if (is.null(v)) return("None")
  if (isTRUE(attr(v, "pyset"))) {
    vals <- attr(v, "pyset"); attr(v, "pyset") <- NULL
    if (!length(v)) return("set()")
    return(paste0("{", paste(sort(vapply(v, .fmt, "", top = FALSE), method = "radix"), collapse = ", "), "}"))
  }
  if (inherits(v, "physics_obj") && !inherits(v, "pynamedtuple")) return(format(v))
  if (inherits(v, c("pytuple", "pynamedtuple"))) {
    s <- vapply(v, .fmt, "", top = FALSE)
    return(paste0("(", paste(s, collapse = ", "), if (length(s) == 1L) ",", ")"))
  }
  if (inherits(v, "pydict")) {
    if (!length(v$keys)) return("{}")
    return(paste0("{", paste(vapply(seq_along(v$keys), function(i)
      paste0(.fmt(v$keys[[i]], FALSE), ": ", .fmt(v$values[[i]], FALSE)), ""), collapse = ", "), "}"))
  }
  if (inherits(v, "sym") || inherits(v, "symmatrix")) return(sstr(v))
  if (is.character(v) && length(v) == 1L) return(if (top) v else paste0("'", gsub("\\", "\\\\", v, fixed = TRUE), "'"))
  if (is.logical(v) && length(v) == 1L) return(if (is.na(v)) "None" else if (v) "True" else "False")
  if (is.numeric(v) && length(v) == 1L) return(if (v == round(v) && abs(v) < 1e15) format(v, scientific = FALSE) else format(v, digits = 17))
  if (is.function(v)) return("<function>")
  if (is.list(v) && !is.object(v)) return(paste0("[", paste(vapply(v, .fmt, "", top = FALSE), collapse = ", "), "]"))
  paste(format(v), collapse = " ")
}

.err_class <- function(e) {
  cl <- setdiff(class(e), c("simpleError", "error", "condition", "sympy_error"))
  if (length(cl)) cl[1L] else "Error"
}

.dec <- function(s) {
  s <- gsub("\\\\n", "\n", s)
  gsub("\\\\\\\\", "\\\\", s)
}
.case <- function(code, env) {
  res <- tryCatch({
    v <- .pywarn(eval(parse(text = .dec(code)), envir = env))
    gsub("**", "^", .fmt(v), fixed = TRUE)
  }, error = function(e) paste0("ERR ", .err_class(e), ": ", conditionMessage(e)))
  gsub("\n", "\\n", res, fixed = TRUE)
}
.run <- function(code, env) .pywarn(eval(parse(text = .dec(code)), envir = env))

.menv <- new.env(parent = H)
local({ s <- get(".pv_state", envir = ns); s$simp <- FALSE; s$t <- NULL })
try(.run("A <- ReferenceFrame(\"A\")", .menv), silent = TRUE)
test_that("vector/tests/test_dyadic.py::test_dyadic matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run("d1 <- (pyattr(A, \"x\") | pyattr(A, \"x\"))", e)
  .run("d2 <- (pyattr(A, \"y\") | pyattr(A, \"y\"))", e)
  .run("d3 <- (pyattr(A, \"x\") | pyattr(A, \"y\"))", e)
  expect_equal(.case("(d1 * 0)", e), "0")
  expect_equal(.case("d1", e), "(A.x|A.x)")
  expect_equal(.case("(d1 * 2)", e), "2*(A.x|A.x)")
  expect_equal(.case("(d1 / Float(2.0))", e), "0.5*(A.x|A.x)")
  expect_equal(.case("(d1 / 3)", e), "1/3*(A.x|A.x)")
  expect_equal(.case("(d1 & (0 * d1))", e), "0")
  expect_equal(.case("(d1 & d2)", e), "0")
  expect_equal(.case("(d1 & pyattr(A, \"x\"))", e), "A.x")
  expect_equal(.case("pyxor(d1, pyattr(A, \"x\"))", e), "0")
  expect_equal(.case("pyxor(d1, pyattr(A, \"y\"))", e), "(A.x|A.z)")
  expect_equal(.case("pyxor(d1, pyattr(A, \"z\"))", e), "- (A.x|A.y)")
  expect_equal(.case("pyxor(d2, pyattr(A, \"x\"))", e), "- (A.y|A.z)")
  expect_equal(.case("pyxor(pyattr(A, \"x\"), d1)", e), "0")
  expect_equal(.case("pyxor(pyattr(A, \"y\"), d1)", e), "- (A.z|A.x)")
  expect_equal(.case("pyxor(pyattr(A, \"z\"), d1)", e), "(A.y|A.x)")
  expect_equal(.case("(pyattr(A, \"x\") & d1)", e), "A.x")
  expect_equal(.case("(pyattr(A, \"y\") & d1)", e), "0")
  expect_equal(.case("(pyattr(A, \"y\") & d2)", e), "A.y")
  expect_equal(.case("(d1 & d3)", e), "(A.x|A.y)")
  expect_equal(.case("(d3 & d1)", e), "0")
  expect_equal(.case("pycall(d1, \"dt\", A)", e), "0")
  .run("q <- dynamicsymbols(\"q\")", e)
  .run("qd <- dynamicsymbols(\"q\", 1)", e)
  .run("B <- pycall(A, \"orientnew\", \"B\", \"Axis\", list(q, pyattr(A, \"z\")))", e)
  expect_equal(.case("pycall(d1, \"express\", B)", e), "cos(q(t))^2*(B.x|B.x) - sin(q(t))*cos(q(t))*(B.x|B.y) - sin(q(t))*cos(q(t))*(B.y|B.x) + sin(q(t))^2*(B.y|B.y)")
  expect_equal(.case("pycall(d1, \"express\", B)", e), "cos(q(t))^2*(B.x|B.x) - sin(q(t))*cos(q(t))*(B.x|B.y) - sin(q(t))*cos(q(t))*(B.y|B.x) + sin(q(t))^2*(B.y|B.y)")
  expect_equal(.case("pycall(d1, \"express\", B, A)", e), "cos(q(t))*(B.x|A.x) - sin(q(t))*(B.y|A.x)")
  expect_equal(.case("pycall(d1, \"express\", A, B)", e), "cos(q(t))*(A.x|B.x) - sin(q(t))*(A.x|B.y)")
  expect_equal(.case("pycall(d1, \"dt\", B)", e), "- Derivative(q(t), t)*(A.y|A.x) - Derivative(q(t), t)*(A.x|A.y)")
  expect_equal(.case("pycall(d1, \"to_matrix\", A)", e), "Matrix([[1, 0, 0], [0, 0, 0], [0, 0, 0]])")
  expect_equal(.case("pycall(d1, \"to_matrix\", A, B)", e), "Matrix([[cos(q(t)), -sin(q(t)), 0], [0, 0, 0], [0, 0, 0]])")
  expect_equal(.case("pycall(d3, \"to_matrix\", A)", e), "Matrix([[0, 1, 0], [0, 0, 0], [0, 0, 0]])")
  .run(".tmp <- pylist(pysymbols(\"a, b, c, d, e, f\")); a <- .tmp[[1]]; b <- .tmp[[2]]; c <- .tmp[[3]]; d <- .tmp[[4]]; e <- .tmp[[5]]; f <- .tmp[[6]]", e)
  .run("v1 <- (((a * pyattr(A, \"x\")) + (b * pyattr(A, \"y\"))) + (c * pyattr(A, \"z\")))", e)
  .run("v2 <- (((d * pyattr(A, \"x\")) + (e * pyattr(A, \"y\"))) + (f * pyattr(A, \"z\")))", e)
  .run("d4 <- pycall(v1, \"outer\", v2)", e)
  expect_equal(.case("pycall(d4, \"to_matrix\", A)", e), "Matrix([[a*d, a*e, a*f], [b*d, b*e, b*f], [c*d, c*e, c*f]])")
  .run("d5 <- pycall(v1, \"outer\", v1)", e)
  .run("C <- pycall(A, \"orientnew\", \"C\", \"Axis\", list(q, pyattr(A, \"x\")))", e)
  .run("for (.it in pyiter(pyzip(((pycall(C, \"dcm\", A) * pycall(d5, \"to_matrix\", A)) * pyattr(pycall(C, \"dcm\", A), \"T\")), pycall(d5, \"to_matrix\", C)))) { expected <- pyget(.it, 0); actual <- pyget(.it, 1); NULL }", e)
  expect_equal(.case("pycall(d1, \"applyfunc\", 0)", e), "ERR TypeError: `f` must be callable.")
})

test_that("vector/tests/test_dyadic.py::test_dyadic_simplify matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run(".tmp <- pylist(pysymbols(\"x, y, z, k, n, m, w, f, s, A\")); x <- .tmp[[1]]; y <- .tmp[[2]]; z <- .tmp[[3]]; k <- .tmp[[4]]; n <- .tmp[[5]]; m <- .tmp[[6]]; w <- .tmp[[7]]; f <- .tmp[[8]]; s <- .tmp[[9]]; A <- .tmp[[10]]", e)
  .run("N <- ReferenceFrame(\"N\")", e)
  .run("dy <- (pyattr(N, \"x\") | pyattr(N, \"x\"))", e)
  .run("test1 <- (((1 / x) + (1 / y)) * dy)", e)
  expect_equal(.case("((pyattr(N, \"x\") & test1) & pyattr(N, \"x\"))", e), "1/y + 1/x")
  .run("test1 <- pycall(test1, \"simplify\")", e)
  expect_equal(.case("((pyattr(N, \"x\") & test1) & pyattr(N, \"x\"))", e), "(x + y)/(x*y)")
  .run("test2 <- (((pypow(A, 2) * pypow(s, 4)) / (((4 * S_Pi) * k) * pypow(m, 3))) * dy)", e)
  .run("test2 <- pycall(test2, \"simplify\")", e)
  expect_equal(.case("((pyattr(N, \"x\") & test2) & pyattr(N, \"x\"))", e), "A^2*s^4/(4*pi*k*m^3)")
  .run("test3 <- ((((4 + (4 * x)) - (2 * (2 + (2 * x)))) / (2 + (2 * x))) * dy)", e)
  .run("test3 <- pycall(test3, \"simplify\")", e)
  expect_equal(.case("((pyattr(N, \"x\") & test3) & pyattr(N, \"x\"))", e), "0")
  .run("test4 <- (((((((-4) * x) * pypow(y, 2)) - (2 * pypow(y, 3))) - ((2 * pypow(x, 2)) * y)) / pypow((x + y), 2)) * dy)", e)
  .run("test4 <- pycall(test4, \"simplify\")", e)
  expect_equal(.case("((pyattr(N, \"x\") & test4) & pyattr(N, \"x\"))", e), "-2*y")
})

test_that("vector/tests/test_dyadic.py::test_dyadic_subs matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run("N <- ReferenceFrame(\"N\")", e)
  .run("s <- pysymbols(\"s\")", e)
  .run("a <- (s * (pyattr(N, \"x\") | pyattr(N, \"x\")))", e)
  expect_equal(.case("pycall(a, \"subs\", pydict(list(s), list(2)))", e), "2*(N.x|N.x)")
})

test_that("vector/tests/test_dyadic.py::test_check_dyadic matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  expect_equal(.case("`_check_dyadic`(0)", e), "ERR TypeError: A Dyadic must be supplied")
})

test_that("vector/tests/test_dyadic.py::test_dyadic_evalf matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run("N <- ReferenceFrame(\"N\")", e)
  .run("a <- (S_Pi * (pyattr(N, \"x\") | pyattr(N, \"x\")))", e)
  expect_equal(.case("pycall(a, \"evalf\", 3)", e), "3.14*(N.x|N.x)")
  .run("s <- pysymbols(\"s\")", e)
  .run("a <- (((5 * s) * S_Pi) * (pyattr(N, \"x\") | pyattr(N, \"x\")))", e)
  expect_equal(.case("pycall(a, \"evalf\", 2)", e), "16.0*s*(N.x|N.x)")
  expect_equal(.case("pycall(a, \"evalf\", 9, subs = pydict(list(s), list(Float(5.124))))", e), "80.4876038*(N.x|N.x)")
})

test_that("vector/tests/test_dyadic.py::test_dyadic_xreplace matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run(".tmp <- pylist(pysymbols(\"x y z\")); x <- .tmp[[1]]; y <- .tmp[[2]]; z <- .tmp[[3]]", e)
  .run("N <- ReferenceFrame(\"N\")", e)
  .run("D <- outer_sym(pyattr(N, \"x\"), pyattr(N, \"x\"))", e)
  .run("v <- ((x * y) * D)", e)
  expect_equal(.case("pycall(v, \"xreplace\", pydict(list(x), list(cos(x))))", e), "y*cos(x)*(N.x|N.x)")
  expect_equal(.case("pycall(v, \"xreplace\", pydict(list((x * y)), list(S_Pi)))", e), "pi*(N.x|N.x)")
  .run("v <- (pypow((x * y), z) * D)", e)
  expect_equal(.case("pycall(v, \"xreplace\", pydict(list(pypow((x * y), z)), list(1)))", e), "(N.x|N.x)")
  expect_equal(.case("pycall(v, \"xreplace\", pydict(list(x, z), list(1, 0)))", e), "(N.x|N.x)")
  expect_equal(.case("pycall(v, \"xreplace\")", e), "ERR TypeError: Dyadic.xreplace() missing 1 required positional argument: 'rule'")
  expect_equal(.case("pycall(v, \"xreplace\", list(x, y))", e), "ERR TypeError: list indices must be integers or slices, not Symbol")
})

.menv <- new.env(parent = H)
local({ s <- get(".pv_state", envir = ns); s$simp <- FALSE; s$t <- NULL })
try(.run("R <- ReferenceFrame(\"R\"); q <- dynamicsymbols(\"q\"); P <- pycall(R, \"orientnew\", \"P\", \"Axis\", list(q, pyattr(R, \"z\"))); scalar_field <- (((2 * pypow(pyget(R, 0), 2)) * pyget(R, 1)) * pyget(R, 2)); grad_field <- gradient(scalar_field, R); vector_field <- (((pypow(pyget(R, 1), 2) * pyattr(R, \"x\")) + ((3 * pyget(R, 0)) * pyattr(R, \"y\"))) + (((5 * pyget(R, 1)) * pyget(R, 2)) * pyattr(R, \"z\"))); curl_field <- curl(vector_field, R)", .menv), silent = TRUE)
test_that("vector/tests/test_fieldfunctions.py::test_curl matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  expect_equal(.case("curl(PhysVector(0), R)", e), "0")
  expect_equal(.case("curl(pyattr(R, \"x\"), R)", e), "0")
  expect_equal(.case("curl(((2 * pypow(pyget(R, 1), 2)) * pyattr(R, \"y\")), R)", e), "0")
  expect_equal(.case("curl(((pyget(R, 0) * pyget(R, 1)) * pyattr(R, \"z\")), R)", e), "R_x*R.x - R_y*R.y")
  expect_equal(.case("curl((((pyget(R, 0) * pyget(R, 1)) * pyget(R, 2)) * ((pyattr(R, \"x\") + pyattr(R, \"y\")) + pyattr(R, \"z\"))), R)", e), "(-R_x*R_y + R_x*R_z)*R.x + (R_x*R_y - R_y*R_z)*R.y + (-R_x*R_z + R_y*R_z)*R.z")
  expect_equal(.case("curl(((2 * pypow(pyget(R, 0), 2)) * pyattr(R, \"y\")), R)", e), "4*R_x*R.z")
  expect_equal(.case("curl(((pypow(pyget(P, 0), 2) * pyattr(R, \"x\")) + pyattr(P, \"y\")), R)", e), "- 2*(R_x*cos(q(t)) + R_y*sin(q(t)))*sin(q(t))*R.z")
  expect_equal(.case("curl((pyget(P, 0) * pyattr(R, \"y\")), P)", e), "cos(q(t))*P.z")
})

test_that("vector/tests/test_fieldfunctions.py::test_divergence matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  expect_equal(.case("pyis(divergence(PhysVector(0), R), S_Zero)", e), "True")
  expect_equal(.case("pyis(divergence(pyattr(R, \"x\"), R), S_Zero)", e), "True")
  expect_equal(.case("divergence((pypow(pyget(R, 0), 2) * pyattr(R, \"x\")), R)", e), "2*R_x")
  expect_equal(.case("divergence((((pyget(R, 0) * pyget(R, 1)) * pyget(R, 2)) * ((pyattr(R, \"x\") + pyattr(R, \"y\")) + pyattr(R, \"z\"))), R)", e), "R_x*R_y + R_x*R_z + R_y*R_z")
  expect_equal(.case("divergence(((1 / ((pyget(R, 0) * pyget(R, 1)) * pyget(R, 2))) * ((pyattr(R, \"x\") + pyattr(R, \"y\")) + pyattr(R, \"z\"))), R)", e), "-1/(R_x*R_y*R_z^2) - 1/(R_x*R_y^2*R_z) - 1/(R_x^2*R_y*R_z)")
  .run("v <- (((pyget(P, 0) * pyattr(P, \"x\")) + (pyget(P, 1) * pyattr(P, \"y\"))) + (pyget(P, 2) * pyattr(P, \"z\")))", e)
  expect_equal(.case("divergence(v, P)", e), "3")
  expect_equal(.case("pycall(divergence(v, R), \"simplify\")", e), "3")
  expect_equal(.case("divergence(((pyget(P, 0) * pyattr(R, \"x\")) + (pyget(R, 0) * pyattr(P, \"x\"))), R)", e), "2*cos(q(t))")
})

test_that("vector/tests/test_fieldfunctions.py::test_gradient matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run("a <- Symbol(\"a\")", e)
  expect_equal(.case("gradient(0, R)", e), "0")
  expect_equal(.case("gradient(pyget(R, 0), R)", e), "R.x")
  expect_equal(.case("gradient(((pyget(R, 0) * pyget(R, 1)) * pyget(R, 2)), R)", e), "R_y*R_z*R.x + R_x*R_z*R.y + R_x*R_y*R.z")
  expect_equal(.case("gradient((2 * pypow(pyget(R, 0), 2)), R)", e), "4*R_x*R.x")
  expect_equal(.case("gradient(((a * sin(pyget(R, 1))) / pyget(R, 0)), R)", e), "- a*sin(R_y)/R_x^2*R.x + a*cos(R_y)/R_x*R.y")
  expect_equal(.case("gradient((pyget(P, 0) * pyget(P, 1)), R)", e), "((-R_x*sin(q(t)) + R_y*cos(q(t)))*cos(q(t)) - (R_x*cos(q(t)) + R_y*sin(q(t)))*sin(q(t)))*R.x + ((-R_x*sin(q(t)) + R_y*cos(q(t)))*sin(q(t)) + (R_x*cos(q(t)) + R_y*sin(q(t)))*cos(q(t)))*R.y")
  expect_equal(.case("gradient((pyget(P, 0) * pyget(R, 2)), P)", e), "P_z*P.x + P_x*P.z")
})

test_that("vector/tests/test_fieldfunctions.py::test_conservative matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  expect_equal(.case("pyis(is_conservative(0), TRUE)", e), "True")
  expect_equal(.case("pyis(is_conservative(pyattr(R, \"x\")), TRUE)", e), "True")
  expect_equal(.case("pyis(is_conservative((((2 * pyattr(R, \"x\")) + (3 * pyattr(R, \"y\"))) + (4 * pyattr(R, \"z\")))), TRUE)", e), "True")
  expect_equal(.case("pyis(is_conservative(((((pyget(R, 1) * pyget(R, 2)) * pyattr(R, \"x\")) + ((pyget(R, 0) * pyget(R, 2)) * pyattr(R, \"y\"))) + ((pyget(R, 0) * pyget(R, 1)) * pyattr(R, \"z\")))), TRUE)", e), "True")
  expect_equal(.case("pyis(is_conservative((pyget(R, 0) * pyattr(R, \"y\"))), FALSE)", e), "True")
  expect_equal(.case("pyis(is_conservative(grad_field), TRUE)", e), "True")
  expect_equal(.case("pyis(is_conservative(curl_field), FALSE)", e), "True")
  expect_equal(.case("pyis(is_conservative((((((4 * pyget(R, 0)) * pyget(R, 1)) * pyget(R, 2)) * pyattr(R, \"x\")) + (((2 * pypow(pyget(R, 0), 2)) * pyget(R, 2)) * pyattr(R, \"y\")))), FALSE)", e), "True")
  expect_equal(.case("pyis(is_conservative(((pyget(R, 2) * pyattr(P, \"x\")) + (pyget(P, 0) * pyattr(R, \"z\")))), TRUE)", e), "True")
})

test_that("vector/tests/test_fieldfunctions.py::test_solenoidal matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  expect_equal(.case("pyis(is_solenoidal(0), TRUE)", e), "True")
  expect_equal(.case("pyis(is_solenoidal(pyattr(R, \"x\")), TRUE)", e), "True")
  expect_equal(.case("pyis(is_solenoidal((((2 * pyattr(R, \"x\")) + (3 * pyattr(R, \"y\"))) + (4 * pyattr(R, \"z\")))), TRUE)", e), "True")
  expect_equal(.case("pyis(is_solenoidal(((((pyget(R, 1) * pyget(R, 2)) * pyattr(R, \"x\")) + ((pyget(R, 0) * pyget(R, 2)) * pyattr(R, \"y\"))) + ((pyget(R, 0) * pyget(R, 1)) * pyattr(R, \"z\")))), TRUE)", e), "True")
  expect_equal(.case("pyis(is_solenoidal((pyget(R, 1) * pyattr(R, \"y\"))), FALSE)", e), "True")
  expect_equal(.case("pyis(is_solenoidal(grad_field), FALSE)", e), "True")
  expect_equal(.case("pyis(is_solenoidal(curl_field), TRUE)", e), "True")
  expect_equal(.case("pyis(is_solenoidal(((((-2) * pyget(R, 1)) + 3) * pyattr(R, \"z\"))), TRUE)", e), "True")
  expect_equal(.case("pyis(is_solenoidal((((cos(q) * pyattr(R, \"x\")) + (sin(q) * pyattr(R, \"y\"))) + (cos(q) * pyattr(P, \"z\")))), TRUE)", e), "True")
  expect_equal(.case("pyis(is_solenoidal(((pyget(R, 2) * pyattr(P, \"x\")) + (pyget(P, 0) * pyattr(R, \"z\")))), TRUE)", e), "True")
})

test_that("vector/tests/test_fieldfunctions.py::test_scalar_potential matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  expect_equal(.case("scalar_potential(0, R)", e), "0")
  expect_equal(.case("scalar_potential(pyattr(R, \"x\"), R)", e), "R_x")
  expect_equal(.case("scalar_potential(pyattr(R, \"y\"), R)", e), "R_y")
  expect_equal(.case("scalar_potential(pyattr(R, \"z\"), R)", e), "R_z")
  expect_equal(.case("scalar_potential(((((pyget(R, 1) * pyget(R, 2)) * pyattr(R, \"x\")) + ((pyget(R, 0) * pyget(R, 2)) * pyattr(R, \"y\"))) + ((pyget(R, 0) * pyget(R, 1)) * pyattr(R, \"z\"))), R)", e), "R_x*R_y*R_z")
  expect_equal(.case("scalar_potential(grad_field, R)", e), "2*R_x^2*R_y*R_z")
  expect_equal(.case("scalar_potential(((pyget(R, 2) * pyattr(P, \"x\")) + (pyget(P, 0) * pyattr(R, \"z\"))), R)", e), "R_x*R_z*cos(q(t)) + R_y*R_z*sin(q(t))")
  expect_equal(.case("scalar_potential(((pyget(R, 2) * pyattr(P, \"x\")) + (pyget(P, 0) * pyattr(R, \"z\"))), P)", e), "P_x*P_z")
  expect_equal(.case("scalar_potential((pyget(R, 0) * pyattr(R, \"y\")), R)", e), "ERR ValueError: Field is not conservative")
})

test_that("vector/tests/test_fieldfunctions.py::test_scalar_potential_difference matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run("origin <- PhysPoint(\"O\")", e)
  .run("point1 <- pycall(origin, \"locatenew\", \"P1\", (((1 * pyattr(R, \"x\")) + (2 * pyattr(R, \"y\"))) + (3 * pyattr(R, \"z\"))))", e)
  .run("point2 <- pycall(origin, \"locatenew\", \"P2\", (((4 * pyattr(R, \"x\")) + (5 * pyattr(R, \"y\"))) + (6 * pyattr(R, \"z\"))))", e)
  .run("genericpointR <- pycall(origin, \"locatenew\", \"RP\", (((pyget(R, 0) * pyattr(R, \"x\")) + (pyget(R, 1) * pyattr(R, \"y\"))) + (pyget(R, 2) * pyattr(R, \"z\"))))", e)
  .run("genericpointP <- pycall(origin, \"locatenew\", \"PP\", (((pyget(P, 0) * pyattr(P, \"x\")) + (pyget(P, 1) * pyattr(P, \"y\"))) + (pyget(P, 2) * pyattr(P, \"z\"))))", e)
  expect_equal(.case("scalar_potential_difference(S_Zero, R, point1, point2, origin)", e), "0")
  expect_equal(.case("scalar_potential_difference(scalar_field, R, origin, genericpointR, origin)", e), "2*R_x^2*R_y*R_z")
  expect_equal(.case("scalar_potential_difference(grad_field, R, origin, genericpointR, origin)", e), "2*R_x^2*R_y*R_z")
  expect_equal(.case("scalar_potential_difference(grad_field, R, point1, point2, origin)", e), "948")
  expect_equal(.case("scalar_potential_difference(((((pyget(R, 1) * pyget(R, 2)) * pyattr(R, \"x\")) + ((pyget(R, 0) * pyget(R, 2)) * pyattr(R, \"y\"))) + ((pyget(R, 0) * pyget(R, 1)) * pyattr(R, \"z\"))), R, point1, genericpointR, origin)", e), "R_x*R_y*R_z - 6")
  .run("potential_diff_P <- (((2 * pyget(P, 2)) * ((pyget(P, 0) * sin(q)) + (pyget(P, 1) * cos(q)))) * pypow(((pyget(P, 0) * cos(q)) - (pyget(P, 1) * sin(q))), 2))", e)
  expect_equal(.case("pycall(scalar_potential_difference(grad_field, P, origin, genericpointP, origin), \"simplify\")", e), "2*P_z*(P_x*sin(q(t)) + P_y*cos(q(t)))*(P_x*cos(q(t)) - P_y*sin(q(t)))^2")
})

.menv <- new.env(parent = H)
local({ s <- get(".pv_state", envir = ns); s$simp <- FALSE; s$t <- NULL })
try(.run("NULL", .menv), silent = TRUE)
test_that("vector/tests/test_frame.py::test_dict_list matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run("A <- ReferenceFrame(\"A\")", e)
  .run("B <- ReferenceFrame(\"B\")", e)
  .run("C <- ReferenceFrame(\"C\")", e)
  .run("D <- ReferenceFrame(\"D\")", e)
  .run("E <- ReferenceFrame(\"E\")", e)
  .run("`F` <- ReferenceFrame(\"F\")", e)
  .run("pycall(B, \"orient_axis\", A, pyattr(A, \"x\"), Float(1.0))", e)
  .run("pycall(C, \"orient_axis\", B, pyattr(B, \"x\"), Float(1.0))", e)
  .run("pycall(D, \"orient_axis\", C, pyattr(C, \"x\"), Float(1.0))", e)
  expect_equal(.case("pycall(D, \"_dict_list\", A, 0)", e), "[D, C, B, A]")
  .run("pycall(E, \"orient_axis\", D, pyattr(D, \"x\"), Float(1.0))", e)
  expect_equal(.case("pycall(C, \"_dict_list\", A, 0)", e), "[C, B, A]")
  expect_equal(.case("pycall(C, \"_dict_list\", E, 0)", e), "[C, D, E]")
  expect_equal(.case("pycall(C, \"_dict_list\", E, 5)", e), "ERR ValueError: Valid values for num are 0, 1, or 2.")
  expect_equal(.case("pycall(`F`, \"_dict_list\", A, 0)", e), "ERR ValueError: No connecting orientation path found between F and A.")
})

test_that("vector/tests/test_frame.py::test_coordinate_vars matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run("\"Tests the coordinate variables functionality\"", e)
  .run("A <- ReferenceFrame(\"A\")", e)
  expect_equal(.case("CoordinateSym(\"Ax\", A, 0)", e), "Ax")
  expect_equal(.case("CoordinateSym(\"Ax\", A, 1)", e), "Ax")
  expect_equal(.case("CoordinateSym(\"Ax\", A, 2)", e), "Ax")
  expect_equal(.case("CoordinateSym(\"Ax\", A, 3)", e), "ERR ValueError: Invalid index specified")
  .run("q <- dynamicsymbols(\"q\")", e)
  .run("qd <- dynamicsymbols(\"q\", 1)", e)
  expect_equal(.case("(pybool(pyisinstance(pyget(A, 0), CoordinateSym)) && pybool(pyisinstance(pyget(A, 0), CoordinateSym)) && pybool(pyisinstance(pyget(A, 0), CoordinateSym)))", e), "True")
  expect_equal(.case("pycall(A, \"variable_map\", A)", e), "{A_x: A_x, A_y: A_y, A_z: A_z}")
  expect_equal(.case("pyattr(pyget(A, 0), \"frame\")", e), "A")
  .run("B <- pycall(A, \"orientnew\", \"B\", \"Axis\", list(q, pyattr(A, \"z\")))", e)
  expect_equal(.case("pycall(B, \"variable_map\", A)", e), "{B_x: A_x*cos(q(t)) + A_y*sin(q(t)), B_y: -A_x*sin(q(t)) + A_y*cos(q(t)), B_z: A_z}")
  expect_equal(.case("pycall(A, \"variable_map\", B)", e), "{A_x: B_x*cos(q(t)) - B_y*sin(q(t)), A_y: B_x*sin(q(t)) + B_y*cos(q(t)), A_z: B_z}")
  expect_equal(.case("time_derivative(pyget(B, 0), A)", e), "-A_x*sin(q(t))*Derivative(q(t), t) + A_y*cos(q(t))*Derivative(q(t), t)")
  expect_equal(.case("time_derivative(pyget(B, 1), A)", e), "-A_x*cos(q(t))*Derivative(q(t), t) - A_y*sin(q(t))*Derivative(q(t), t)")
  expect_equal(.case("time_derivative(pyget(B, 2), A)", e), "0")
  expect_equal(.case("express(pyget(B, 0), A, variables = TRUE)", e), "A_x*cos(q(t)) + A_y*sin(q(t))")
  expect_equal(.case("express(pyget(B, 1), A, variables = TRUE)", e), "-A_x*sin(q(t)) + A_y*cos(q(t))")
  expect_equal(.case("express(pyget(B, 2), A, variables = TRUE)", e), "A_z")
  expect_equal(.case("time_derivative((((pyget(A, 0) * pyattr(A, \"x\")) + (pyget(A, 1) * pyattr(A, \"y\"))) + (pyget(A, 2) * pyattr(A, \"z\"))), B)", e), "A_y*Derivative(q(t), t)*A.x - A_x*Derivative(q(t), t)*A.y")
  expect_equal(.case("time_derivative((((pyget(B, 0) * pyattr(B, \"x\")) + (pyget(B, 1) * pyattr(B, \"y\"))) + (pyget(B, 2) * pyattr(B, \"z\"))), A)", e), "- B_y*Derivative(q(t), t)*B.x + B_x*Derivative(q(t), t)*B.y")
  expect_equal(.case("express(((pyget(B, 0) * pyget(B, 1)) * pyget(B, 2)), A, variables = TRUE)", e), "A_z*(-A_x*sin(q(t)) + A_y*cos(q(t)))*(A_x*cos(q(t)) + A_y*sin(q(t)))")
  expect_equal(.case("pycall((time_derivative(((pyget(B, 0) * pyget(B, 1)) * pyget(B, 2)), A) - ((pyget(A, 2) * ((((-pypow(pyget(A, 0), 2)) * cos((2 * q))) - (((2 * pyget(A, 0)) * pyget(A, 1)) * sin((2 * q)))) + (pypow(pyget(A, 1), 2) * cos((2 * q))))) * qd)), \"trigsimp\")", e), "0")
  expect_equal(.case("express((((pyget(B, 0) * pyattr(B, \"x\")) + (pyget(B, 1) * pyattr(B, \"y\"))) + (pyget(B, 2) * pyattr(B, \"z\"))), A)", e), "(B_x*cos(q(t)) - B_y*sin(q(t)))*A.x + (B_x*sin(q(t)) + B_y*cos(q(t)))*A.y + B_z*A.z")
  expect_equal(.case("pycall(express((((pyget(B, 0) * pyattr(B, \"x\")) + (pyget(B, 1) * pyattr(B, \"y\"))) + (pyget(B, 2) * pyattr(B, \"z\"))), A, variables = TRUE), \"simplify\")", e), "A_x*A.x + A_y*A.y + A_z*A.z")
  expect_equal(.case("express((((pyget(A, 0) * pyattr(A, \"x\")) + (pyget(A, 1) * pyattr(A, \"y\"))) + (pyget(A, 2) * pyattr(A, \"z\"))), B)", e), "(A_x*cos(q(t)) + A_y*sin(q(t)))*B.x + (-A_x*sin(q(t)) + A_y*cos(q(t)))*B.y + A_z*B.z")
  expect_equal(.case("pycall(express((((pyget(A, 0) * pyattr(A, \"x\")) + (pyget(A, 1) * pyattr(A, \"y\"))) + (pyget(A, 2) * pyattr(A, \"z\"))), B, variables = TRUE), \"simplify\")", e), "B_x*B.x + B_y*B.y + B_z*B.z")
  .run("N <- pycall(B, \"orientnew\", \"N\", \"Axis\", list((-q), pyattr(B, \"z\")))", e)
  expect_equal(.case("pydict(lapply(pyiter(pycall(pycall(N, \"variable_map\", A), \"items\")), function(.it) { k <- pyget(.it, 0); v <- pyget(.it, 1); k }), lapply(pyiter(pycall(pycall(N, \"variable_map\", A), \"items\")), function(.it) { k <- pyget(.it, 0); v <- pyget(.it, 1); pycall(v, \"simplify\") }))", e), "{N_x: A_x, N_y: A_y, N_z: A_z}")
  .run("C <- pycall(A, \"orientnew\", \"C\", \"Axis\", list(q, ((pyattr(A, \"x\") + pyattr(A, \"y\")) + pyattr(A, \"z\"))))", e)
  .run("mapping <- pycall(A, \"variable_map\", C)", e)
  expect_equal(.case("pytrigsimp(pyget(mapping, pyget(A, 0)))", e), "2*C_x*cos(q(t))/3 + C_x/3 - 2*C_y*sin(q(t) + pi/6)/3 + C_y/3 - 2*C_z*cos(q(t) + pi/3)/3 + C_z/3")
  expect_equal(.case("pytrigsimp(pyget(mapping, pyget(A, 1)))", e), "-2*C_x*cos(q(t) + pi/3)/3 + C_x/3 + 2*C_y*cos(q(t))/3 + C_y/3 - 2*C_z*sin(q(t) + pi/6)/3 + C_z/3")
  expect_equal(.case("pytrigsimp(pyget(mapping, pyget(A, 2)))", e), "-2*C_x*sin(q(t) + pi/6)/3 + C_x/3 - 2*C_y*cos(q(t) + pi/3)/3 + C_y/3 + 2*C_z*cos(q(t))/3 + C_z/3")
})

test_that("vector/tests/test_frame.py::test_ang_vel matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run(".tmp <- pylist(dynamicsymbols(\"q1 q2 q3 q4\")); q1 <- .tmp[[1]]; q2 <- .tmp[[2]]; q3 <- .tmp[[3]]; q4 <- .tmp[[4]]", e)
  .run(".tmp <- pylist(dynamicsymbols(\"q1 q2 q3 q4\", 1)); q1d <- .tmp[[1]]; q2d <- .tmp[[2]]; q3d <- .tmp[[3]]; q4d <- .tmp[[4]]", e)
  .run("N <- ReferenceFrame(\"N\")", e)
  .run("A <- pycall(N, \"orientnew\", \"A\", \"Axis\", list(q1, pyattr(N, \"z\")))", e)
  .run("B <- pycall(A, \"orientnew\", \"B\", \"Axis\", list(q2, pyattr(A, \"x\")))", e)
  .run("C <- pycall(B, \"orientnew\", \"C\", \"Axis\", list(q3, pyattr(B, \"y\")))", e)
  .run("D <- pycall(N, \"orientnew\", \"D\", \"Axis\", list(q4, pyattr(N, \"y\")))", e)
  .run(".tmp <- pylist(dynamicsymbols(\"u1 u2 u3\")); u1 <- .tmp[[1]]; u2 <- .tmp[[2]]; u3 <- .tmp[[3]]", e)
  expect_equal(.case("pycall(A, \"ang_vel_in\", N)", e), "Derivative(q1(t), t)*N.z")
  expect_equal(.case("pycall(B, \"ang_vel_in\", N)", e), "Derivative(q1(t), t)*N.z + Derivative(q2(t), t)*A.x")
  expect_equal(.case("pycall(C, \"ang_vel_in\", N)", e), "Derivative(q1(t), t)*N.z + Derivative(q2(t), t)*A.x + Derivative(q3(t), t)*B.y")
  .run("A2 <- pycall(N, \"orientnew\", \"A2\", \"Axis\", list(q4, pyattr(N, \"y\")))", e)
  expect_equal(.case("pycall(N, \"ang_vel_in\", N)", e), "0")
  expect_equal(.case("pycall(N, \"ang_vel_in\", A)", e), "- Derivative(q1(t), t)*N.z")
  expect_equal(.case("pycall(N, \"ang_vel_in\", B)", e), "- Derivative(q1(t), t)*N.z - Derivative(q2(t), t)*A.x")
  expect_equal(.case("pycall(N, \"ang_vel_in\", C)", e), "- Derivative(q1(t), t)*N.z - Derivative(q2(t), t)*A.x - Derivative(q3(t), t)*B.y")
  expect_equal(.case("pycall(N, \"ang_vel_in\", A2)", e), "- Derivative(q4(t), t)*N.y")
  expect_equal(.case("pycall(A, \"ang_vel_in\", N)", e), "Derivative(q1(t), t)*N.z")
  expect_equal(.case("pycall(A, \"ang_vel_in\", A)", e), "0")
  expect_equal(.case("pycall(A, \"ang_vel_in\", B)", e), "- Derivative(q2(t), t)*A.x")
  expect_equal(.case("pycall(A, \"ang_vel_in\", C)", e), "- Derivative(q2(t), t)*A.x - Derivative(q3(t), t)*B.y")
  expect_equal(.case("pycall(A, \"ang_vel_in\", A2)", e), "- Derivative(q4(t), t)*N.y + Derivative(q1(t), t)*N.z")
  expect_equal(.case("pycall(B, \"ang_vel_in\", N)", e), "Derivative(q1(t), t)*N.z + Derivative(q2(t), t)*A.x")
  expect_equal(.case("pycall(B, \"ang_vel_in\", A)", e), "Derivative(q2(t), t)*A.x")
  expect_equal(.case("pycall(B, \"ang_vel_in\", B)", e), "0")
  expect_equal(.case("pycall(B, \"ang_vel_in\", C)", e), "- Derivative(q3(t), t)*B.y")
  expect_equal(.case("pycall(B, \"ang_vel_in\", A2)", e), "- Derivative(q4(t), t)*N.y + Derivative(q1(t), t)*N.z + Derivative(q2(t), t)*A.x")
  expect_equal(.case("pycall(C, \"ang_vel_in\", N)", e), "Derivative(q1(t), t)*N.z + Derivative(q2(t), t)*A.x + Derivative(q3(t), t)*B.y")
  expect_equal(.case("pycall(C, \"ang_vel_in\", A)", e), "Derivative(q2(t), t)*A.x + Derivative(q3(t), t)*B.y")
  expect_equal(.case("pycall(C, \"ang_vel_in\", B)", e), "Derivative(q3(t), t)*B.y")
  expect_equal(.case("pycall(C, \"ang_vel_in\", C)", e), "0")
  expect_equal(.case("pycall(C, \"ang_vel_in\", A2)", e), "- Derivative(q4(t), t)*N.y + Derivative(q1(t), t)*N.z + Derivative(q2(t), t)*A.x + Derivative(q3(t), t)*B.y")
  expect_equal(.case("pycall(A2, \"ang_vel_in\", N)", e), "Derivative(q4(t), t)*N.y")
  expect_equal(.case("pycall(A2, \"ang_vel_in\", A)", e), "Derivative(q4(t), t)*N.y - Derivative(q1(t), t)*N.z")
  expect_equal(.case("pycall(A2, \"ang_vel_in\", B)", e), "Derivative(q4(t), t)*N.y - Derivative(q1(t), t)*N.z - Derivative(q2(t), t)*A.x")
  expect_equal(.case("pycall(A2, \"ang_vel_in\", C)", e), "Derivative(q4(t), t)*N.y - Derivative(q1(t), t)*N.z - Derivative(q2(t), t)*A.x - Derivative(q3(t), t)*B.y")
  expect_equal(.case("pycall(A2, \"ang_vel_in\", A2)", e), "0")
  .run("pycall(C, \"set_ang_vel\", N, (((u1 * pyattr(C, \"x\")) + (u2 * pyattr(C, \"y\"))) + (u3 * pyattr(C, \"z\"))))", e)
  expect_equal(.case("pycall(C, \"ang_vel_in\", N)", e), "u1(t)*C.x + u2(t)*C.y + u3(t)*C.z")
  expect_equal(.case("pycall(N, \"ang_vel_in\", C)", e), "- u1(t)*C.x - u2(t)*C.y - u3(t)*C.z")
  expect_equal(.case("pycall(C, \"ang_vel_in\", D)", e), "- Derivative(q4(t), t)*N.y + u1(t)*C.x + u2(t)*C.y + u3(t)*C.z")
  expect_equal(.case("pycall(D, \"ang_vel_in\", C)", e), "Derivative(q4(t), t)*N.y - u1(t)*C.x - u2(t)*C.y - u3(t)*C.z")
  .run("q0 <- dynamicsymbols(\"q0\")", e)
  .run("q0d <- dynamicsymbols(\"q0\", 1)", e)
  .run("E <- pycall(N, \"orientnew\", \"E\", \"Quaternion\", tup(q0, q1, q2, q3))", e)
  expect_equal(.case("pycall(E, \"ang_vel_in\", N)", e), "(2*q0(t)*Derivative(q1(t), t) - 2*q1(t)*Derivative(q0(t), t) - 2*q2(t)*Derivative(q3(t), t) + 2*q3(t)*Derivative(q2(t), t))*E.x + (2*q0(t)*Derivative(q2(t), t) + 2*q1(t)*Derivative(q3(t), t) - 2*q2(t)*Derivative(q0(t), t) - 2*q3(t)*Derivative(q1(t), t))*E.y + (2*q0(t)*Derivative(q3(t), t) - 2*q1(t)*Derivative(q2(t), t) + 2*q2(t)*Derivative(q1(t), t) - 2*q3(t)*Derivative(q0(t), t))*E.z")
  .run("`F` <- pycall(N, \"orientnew\", \"F\", \"Body\", tup(q1, q2, q3), 313)", e)
  expect_equal(.case("pycall(`F`, \"ang_vel_in\", N)", e), "(sin(q2(t))*sin(q3(t))*Derivative(q1(t), t) + cos(q3(t))*Derivative(q2(t), t))*F.x + (sin(q2(t))*cos(q3(t))*Derivative(q1(t), t) - sin(q3(t))*Derivative(q2(t), t))*F.y + (cos(q2(t))*Derivative(q1(t), t) + Derivative(q3(t), t))*F.z")
  .run("G <- pycall(N, \"orientnew\", \"G\", \"Axis\", tup(q1, (pyattr(N, \"x\") + pyattr(N, \"y\"))))", e)
  expect_equal(.case("pycall(G, \"ang_vel_in\", N)", e), "sqrt(2)*Derivative(q1(t), t)/2*N.x + sqrt(2)*Derivative(q1(t), t)/2*N.y")
  expect_equal(.case("pycall(N, \"ang_vel_in\", G)", e), "- sqrt(2)*Derivative(q1(t), t)/2*N.x - sqrt(2)*Derivative(q1(t), t)/2*N.y")
})

test_that("vector/tests/test_frame.py::test_dcm matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run(".tmp <- pylist(dynamicsymbols(\"q1 q2 q3 q4\")); q1 <- .tmp[[1]]; q2 <- .tmp[[2]]; q3 <- .tmp[[3]]; q4 <- .tmp[[4]]", e)
  .run("N <- ReferenceFrame(\"N\")", e)
  .run("A <- pycall(N, \"orientnew\", \"A\", \"Axis\", list(q1, pyattr(N, \"z\")))", e)
  .run("B <- pycall(A, \"orientnew\", \"B\", \"Axis\", list(q2, pyattr(A, \"x\")))", e)
  .run("C <- pycall(B, \"orientnew\", \"C\", \"Axis\", list(q3, pyattr(B, \"y\")))", e)
  .run("D <- pycall(N, \"orientnew\", \"D\", \"Axis\", list(q4, pyattr(N, \"y\")))", e)
  .run("E <- pycall(N, \"orientnew\", \"E\", \"Space\", list(q1, q2, q3), \"123\")", e)
  expect_equal(.case("pycall(N, \"dcm\", C)", e), "Matrix([[-sin(q1(t))*sin(q2(t))*sin(q3(t)) + cos(q1(t))*cos(q3(t)), -sin(q1(t))*cos(q2(t)), sin(q1(t))*sin(q2(t))*cos(q3(t)) + sin(q3(t))*cos(q1(t))], [sin(q1(t))*cos(q3(t)) + sin(q2(t))*sin(q3(t))*cos(q1(t)), cos(q1(t))*cos(q2(t)), sin(q1(t))*sin(q3(t)) - sin(q2(t))*cos(q1(t))*cos(q3(t))], [-sin(q3(t))*cos(q2(t)), sin(q2(t)), cos(q2(t))*cos(q3(t))]])")
  .run("test_mat <- (pycall(D, \"dcm\", C) - Matrix(list(list((((cos(q1) * cos(q3)) * cos(q4)) - (sin(q3) * (((-sin(q4)) * cos(q2)) + ((sin(q1) * sin(q2)) * cos(q4))))), (((-sin(q2)) * sin(q4)) - ((sin(q1) * cos(q2)) * cos(q4))), (((sin(q3) * cos(q1)) * cos(q4)) + (cos(q3) * (((-sin(q4)) * cos(q2)) + ((sin(q1) * sin(q2)) * cos(q4)))))), list(((sin(q1) * cos(q3)) + ((sin(q2) * sin(q3)) * cos(q1))), (cos(q1) * cos(q2)), ((sin(q1) * sin(q3)) - ((sin(q2) * cos(q1)) * cos(q3)))), list((((sin(q4) * cos(q1)) * cos(q3)) - (sin(q3) * ((cos(q2) * cos(q4)) + ((sin(q1) * sin(q2)) * sin(q4))))), ((sin(q2) * cos(q4)) - ((sin(q1) * sin(q4)) * cos(q2))), (((sin(q3) * sin(q4)) * cos(q1)) + (cos(q3) * ((cos(q2) * cos(q4)) + ((sin(q1) * sin(q2)) * sin(q4)))))))))", e)
  expect_equal(.case("pycall(test_mat, \"expand\")", e), "Matrix([[0, 0, 0], [0, 0, 0], [0, 0, 0]])")
  expect_equal(.case("pycall(E, \"dcm\", N)", e), "Matrix([[cos(q2(t))*cos(q3(t)), sin(q3(t))*cos(q2(t)), -sin(q2(t))], [sin(q1(t))*sin(q2(t))*cos(q3(t)) - sin(q3(t))*cos(q1(t)), sin(q1(t))*sin(q2(t))*sin(q3(t)) + cos(q1(t))*cos(q3(t)), sin(q1(t))*cos(q2(t))], [sin(q1(t))*sin(q3(t)) + sin(q2(t))*cos(q1(t))*cos(q3(t)), -sin(q1(t))*cos(q3(t)) + sin(q2(t))*sin(q3(t))*cos(q1(t)), cos(q1(t))*cos(q2(t))]])")
})

test_that("vector/tests/test_frame.py::test_w_diff_dcm1 matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run("A <- ReferenceFrame(\"A\")", e)
  .run("B <- ReferenceFrame(\"B\")", e)
  .run(".tmp <- pylist(dynamicsymbols(\"C11 C12 C13\")); c11 <- .tmp[[1]]; c12 <- .tmp[[2]]; c13 <- .tmp[[3]]", e)
  .run(".tmp <- pylist(dynamicsymbols(\"C21 C22 C23\")); c21 <- .tmp[[1]]; c22 <- .tmp[[2]]; c23 <- .tmp[[3]]", e)
  .run(".tmp <- pylist(dynamicsymbols(\"C31 C32 C33\")); c31 <- .tmp[[1]]; c32 <- .tmp[[2]]; c33 <- .tmp[[3]]", e)
  .run(".tmp <- pylist(dynamicsymbols(\"C11 C12 C13\", level = 1)); c11d <- .tmp[[1]]; c12d <- .tmp[[2]]; c13d <- .tmp[[3]]", e)
  .run(".tmp <- pylist(dynamicsymbols(\"C21 C22 C23\", level = 1)); c21d <- .tmp[[1]]; c22d <- .tmp[[2]]; c23d <- .tmp[[3]]", e)
  .run(".tmp <- pylist(dynamicsymbols(\"C31 C32 C33\", level = 1)); c31d <- .tmp[[1]]; c32d <- .tmp[[2]]; c33d <- .tmp[[3]]", e)
  .run("DCM <- Matrix(list(list(c11, c12, c13), list(c21, c22, c23), list(c31, c32, c33)))", e)
  .run("pycall(B, \"orient\", A, \"DCM\", DCM)", e)
  .run("b1a <- pycall(pyattr(B, \"x\"), \"express\", A)", e)
  .run("b2a <- pycall(pyattr(B, \"y\"), \"express\", A)", e)
  .run("b3a <- pycall(pyattr(B, \"z\"), \"express\", A)", e)
  .run("pycall(B, \"set_ang_vel\", A, (((pyattr(B, \"x\") * dot(pycall(b3a, \"dt\", A), pyattr(B, \"y\"))) + (pyattr(B, \"y\") * dot(pycall(b1a, \"dt\", A), pyattr(B, \"z\")))) + (pyattr(B, \"z\") * dot(pycall(b2a, \"dt\", A), pyattr(B, \"x\")))))", e)
  .run("expr <- ((((((c12 * c13d) + (c22 * c23d)) + (c32 * c33d)) * pyattr(B, \"x\")) + ((((c13 * c11d) + (c23 * c21d)) + (c33 * c31d)) * pyattr(B, \"y\"))) + ((((c11 * c12d) + (c21 * c22d)) + (c31 * c32d)) * pyattr(B, \"z\")))", e)
  expect_equal(.case("(pycall(B, \"ang_vel_in\", A) - expr)", e), "0")
})

test_that("vector/tests/test_frame.py::test_w_diff_dcm2 matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run(".tmp <- pylist(dynamicsymbols(\"q1:4\")); q1 <- .tmp[[1]]; q2 <- .tmp[[2]]; q3 <- .tmp[[3]]", e)
  .run("N <- ReferenceFrame(\"N\")", e)
  .run("A <- pycall(N, \"orientnew\", \"A\", \"axis\", list(q1, pyattr(N, \"x\")))", e)
  .run("B <- pycall(A, \"orientnew\", \"B\", \"axis\", list(q2, pyattr(A, \"y\")))", e)
  .run("C <- pycall(B, \"orientnew\", \"C\", \"axis\", list(q3, pyattr(B, \"z\")))", e)
  .run("DCM <- pyattr(pycall(C, \"dcm\", N), \"T\")", e)
  .run("D <- pycall(N, \"orientnew\", \"D\", \"DCM\", DCM)", e)
  expect_equal(.case("pycall(D, \"dcm\", N)", e), "Matrix([[cos(q2(t))*cos(q3(t)), sin(q1(t))*sin(q2(t))*cos(q3(t)) + sin(q3(t))*cos(q1(t)), sin(q1(t))*sin(q3(t)) - sin(q2(t))*cos(q1(t))*cos(q3(t))], [-sin(q3(t))*cos(q2(t)), -sin(q1(t))*sin(q2(t))*sin(q3(t)) + cos(q1(t))*cos(q3(t)), sin(q1(t))*cos(q3(t)) + sin(q2(t))*sin(q3(t))*cos(q1(t))], [sin(q2(t)), -sin(q1(t))*cos(q2(t)), cos(q1(t))*cos(q2(t))]])")
  expect_equal(.case("pycall(pycall((pycall(D, \"ang_vel_in\", N) - pycall(C, \"ang_vel_in\", N)), \"express\", N), \"simplify\")", e), "0")
})

test_that("vector/tests/test_frame.py::test_orientnew_respects_input_indices matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run("N <- ReferenceFrame(\"N\")", e)
  .run("q1 <- dynamicsymbols(\"q1\")", e)
  .run("A <- pycall(N, \"orientnew\", \"a\", \"Axis\", list(q1, pyattr(N, \"z\")))", e)
  .run("minds <- lapply(pyiter(pyattr(N, \"indices\")), function(x) (x + \"1\"))", e)
  .run("B <- pycall(N, \"orientnew\", \"b\", \"Axis\", list(q1, pyattr(N, \"z\")), indices = minds)", e)
  expect_equal(.case("pyattr(N, \"indices\")", e), "['x', 'y', 'z']")
  expect_equal(.case("pyattr(B, \"indices\")", e), "['x1', 'y1', 'z1']")
})

test_that("vector/tests/test_frame.py::test_orientnew_respects_input_latexs matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run("N <- ReferenceFrame(\"N\")", e)
  .run("q1 <- dynamicsymbols(\"q1\")", e)
  .run("A <- pycall(N, \"orientnew\", \"a\", \"Axis\", list(q1, pyattr(N, \"z\")))", e)
  .run("def_latex_vecs <- list((\"\\\\\\\\mathbf{\\\\\\\\hat{%s}_%s}\" %% tup(pycall(pyattr(A, \"name\"), \"lower\"), pyget(pyattr(A, \"indices\"), 0))), (\"\\\\\\\\mathbf{\\\\\\\\hat{%s}_%s}\" %% tup(pycall(pyattr(A, \"name\"), \"lower\"), pyget(pyattr(A, \"indices\"), 1))), (\"\\\\\\\\mathbf{\\\\\\\\hat{%s}_%s}\" %% tup(pycall(pyattr(A, \"name\"), \"lower\"), pyget(pyattr(A, \"indices\"), 2))))", e)
  .run("name <- \"b\"", e)
  .run("indices <- lapply(pyiter(pyattr(N, \"indices\")), function(x) (x + \"1\"))", e)
  .run("new_latex_vecs <- list((\"\\\\\\\\mathbf{\\\\\\\\hat{%s}_{%s}}\" %% tup(pycall(name, \"lower\"), pyget(indices, 0))), (\"\\\\\\\\mathbf{\\\\\\\\hat{%s}_{%s}}\" %% tup(pycall(name, \"lower\"), pyget(indices, 1))), (\"\\\\\\\\mathbf{\\\\\\\\hat{%s}_{%s}}\" %% tup(pycall(name, \"lower\"), pyget(indices, 2))))", e)
  .run("B <- pycall(N, \"orientnew\", name, \"Axis\", list(q1, pyattr(N, \"z\")), latexs = new_latex_vecs)", e)
  expect_equal(.case("pyattr(A, \"latex_vecs\")", e), "['\\\\mathbf{\\\\hat{a}_x}', '\\\\mathbf{\\\\hat{a}_y}', '\\\\mathbf{\\\\hat{a}_z}']")
  expect_equal(.case("pyattr(B, \"latex_vecs\")", e), "['\\\\mathbf{\\\\hat{b}_{x1}}', '\\\\mathbf{\\\\hat{b}_{y1}}', '\\\\mathbf{\\\\hat{b}_{z1}}']")
  expect_equal(.case("pyattr(B, \"indices\")", e), "['x', 'y', 'z']")
})

test_that("vector/tests/test_frame.py::test_partial_velocity matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run("N <- ReferenceFrame(\"N\")", e)
  .run("A <- ReferenceFrame(\"A\")", e)
  .run(".tmp <- pylist(dynamicsymbols(\"u1, u2\")); u1 <- .tmp[[1]]; u2 <- .tmp[[2]]", e)
  .run("pycall(A, \"set_ang_vel\", N, ((u1 * pyattr(A, \"x\")) + (u2 * pyattr(N, \"y\"))))", e)
  expect_equal(.case("pycall(N, \"partial_velocity\", A, u1)", e), "- A.x")
  expect_equal(.case("pycall(N, \"partial_velocity\", A, u1, u2)", e), "(- A.x, - N.y)")
  expect_equal(.case("pycall(A, \"partial_velocity\", N, u1)", e), "A.x")
  expect_equal(.case("pycall(A, \"partial_velocity\", N, u1, u2)", e), "(A.x, N.y)")
  expect_equal(.case("pycall(N, \"partial_velocity\", N, u1)", e), "0")
  expect_equal(.case("pycall(A, \"partial_velocity\", A, u1)", e), "0")
})

test_that("vector/tests/test_frame.py::test_issue_11498 matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run("A <- ReferenceFrame(\"A\")", e)
  .run("B <- ReferenceFrame(\"B\")", e)
  .run("pycall(A, \"orient\", B, \"DCM\", pyeye(3))", e)
  expect_equal(.case("pycall(A, \"dcm\", B)", e), "Matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])")
  expect_equal(.case("pycall(B, \"dcm\", A)", e), "Matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])")
  .run("pycall(A, \"orient\", B, \"DCM\", Matrix(list(list(0, 1, 0), list(0, 0, (-1)), list((-1), 0, 0))))", e)
  expect_equal(.case("pycall(B, \"dcm\", A)", e), "Matrix([[0, 1, 0], [0, 0, -1], [-1, 0, 0]])")
  expect_equal(.case("pycall(A, \"dcm\", B)", e), "Matrix([[0, 0, -1], [1, 0, 0], [0, -1, 0]])")
  expect_equal(.case("pyattr(pycall(B, \"dcm\", A), \"T\")", e), "Matrix([[0, 0, -1], [1, 0, 0], [0, -1, 0]])")
})

test_that("vector/tests/test_frame.py::test_reference_frame matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  expect_equal(.case("ReferenceFrame(0)", e), "ERR TypeError: Need to supply a valid name")
  expect_equal(.case("ReferenceFrame(\"N\", 0)", e), "ERR TypeError: Supply the indices as a list")
  expect_equal(.case("ReferenceFrame(\"N\", list(0, 1))", e), "ERR ValueError: Supply 3 indices")
  expect_equal(.case("ReferenceFrame(\"N\", list(0, 1, 2))", e), "ERR TypeError: Indices must be strings")
  expect_equal(.case("ReferenceFrame(\"N\", list(\"a\", \"b\", \"c\"), 0)", e), "ERR TypeError: Supply the indices as a list")
  expect_equal(.case("ReferenceFrame(\"N\", list(\"a\", \"b\", \"c\"), list(0, 1))", e), "ERR ValueError: Supply 3 indices")
  expect_equal(.case("ReferenceFrame(\"N\", list(\"a\", \"b\", \"c\"), list(0, 1, 2))", e), "ERR TypeError: Latex entries must be strings")
  expect_equal(.case("ReferenceFrame(\"N\", list(\"a\", \"b\", \"c\"), list(\"a\", \"b\", \"c\"), 0)", e), "ERR TypeError: Supply the variable names as a list/tuple")
  expect_equal(.case("ReferenceFrame(\"N\", list(\"a\", \"b\", \"c\"), list(\"a\", \"b\", \"c\"), list(0, 1))", e), "ERR ValueError: Supply 3 variable names")
  expect_equal(.case("ReferenceFrame(\"N\", list(\"a\", \"b\", \"c\"), list(\"a\", \"b\", \"c\"), list(0, 1, 2))", e), "ERR TypeError: Variable names must be strings")
  .run("N <- ReferenceFrame(\"N\")", e)
  expect_equal(.case("pyget(N, 0)", e), "N_x")
  expect_equal(.case("pyget(N, 1)", e), "N_y")
  expect_equal(.case("pyget(N, 2)", e), "N_z")
  expect_equal(.case("pyget(N, 3)", e), "ERR ValueError: Invalid index provided")
  .run("N <- ReferenceFrame(\"N\", list(\"a\", \"b\", \"c\"))", e)
  expect_equal(.case("pyget(N, \"a\")", e), "N['a']")
  expect_equal(.case("pyget(N, \"b\")", e), "N['b']")
  expect_equal(.case("pyget(N, \"c\")", e), "N['c']")
  expect_equal(.case("pyget(N, \"d\")", e), "ERR ValueError: Not a defined index")
  expect_equal(.case("pystr(N)", e), "N")
  .run("A <- ReferenceFrame(\"A\")", e)
  .run("B <- ReferenceFrame(\"B\")", e)
  .run(".tmp <- pylist(pysymbols(\"q0 q1 q2 q3\")); q0 <- .tmp[[1]]; q1 <- .tmp[[2]]; q2 <- .tmp[[3]]; q3 <- .tmp[[4]]", e)
  expect_equal(.case("pycall(A, \"orient\", B, \"DCM\", 0)", e), "ERR TypeError: Amounts must be a SymPy Matrix type object.")
  expect_equal(.case("pycall(B, \"orient\", N, \"Space\", list(q1, q2, q3), \"222\")", e), "ERR TypeError: The supplied order is not an approved type")
  expect_equal(.case("pycall(B, \"orient\", N, \"Axis\", list(q1, (pyattr(N, \"x\") + (2 * pyattr(N, \"y\")))), \"222\")", e), "ERR TypeError: The supplied order is not an approved type")
  expect_equal(.case("pycall(B, \"orient\", N, \"Axis\", q1)", e), "ERR TypeError: 'Symbol' object is not subscriptable")
  expect_equal(.case("pycall(B, \"orient\", N, \"Axis\", list(q1))", e), "ERR IndexError: list index out of range")
  expect_equal(.case("pycall(B, \"orient\", N, \"Quaternion\", list(q0, q1, q2, q3), \"222\")", e), "ERR TypeError: The supplied order is not an approved type")
  expect_equal(.case("pycall(B, \"orient\", N, \"Quaternion\", q0)", e), "ERR TypeError: 'Symbol' object is not iterable")
  expect_equal(.case("pycall(B, \"orient\", N, \"Quaternion\", list(q0, q1, q2))", e), "ERR TypeError: Amounts are a list or tuple of length 4")
  expect_equal(.case("pycall(B, \"orient\", N, \"Foo\", list(q0, q1, q2))", e), "ERR NotImplementedError: That is not an implemented rotation")
  expect_equal(.case("pycall(B, \"orient\", N, \"Body\", list(q1, q2), \"232\")", e), "ERR TypeError: Body orientation takes 3 values & 3 orders")
  expect_equal(.case("pycall(B, \"orient\", N, \"Space\", list(q1, q2), \"232\")", e), "ERR TypeError: Body orientation takes 3 values & 3 orders")
  .run("pycall(N, \"set_ang_acc\", B, 0)", e)
  expect_equal(.case("pycall(N, \"ang_acc_in\", B)", e), "0")
  .run("pycall(N, \"set_ang_vel\", B, 0)", e)
  expect_equal(.case("pycall(N, \"ang_vel_in\", B)", e), "0")
})

test_that("vector/tests/test_frame.py::test_check_frame matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  expect_equal(.case("`_check_frame`(0)", e), "ERR VectorTypeError: \\nExpected an instance of <class\\n'sympy.physics.vector.frame.ReferenceFrame'>, but received object '0'\\nof <class 'int'>.")
})

test_that("vector/tests/test_frame.py::test_dcm_diff_16824 matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run(".tmp <- pylist(dynamicsymbols(\"q1:4\")); q1 <- .tmp[[1]]; q2 <- .tmp[[2]]; q3 <- .tmp[[3]]", e)
  .run("s1 <- sin(q1)", e)
  .run("c1 <- cos(q1)", e)
  .run("s2 <- sin(q2)", e)
  .run("c2 <- cos(q2)", e)
  .run("s3 <- sin(q3)", e)
  .run("c3 <- cos(q3)", e)
  .run("dcm <- Matrix(list(list((c2 * c3), (((s1 * s2) * c3) - (s3 * c1)), (((c1 * s2) * c3) + (s3 * s1))), list((c2 * s3), (((s1 * s2) * s3) + (c3 * c1)), (((c1 * s2) * s3) - (c3 * s1))), list((-s2), (s1 * c2), (c1 * c2))))", e)
  .run("A <- ReferenceFrame(\"A\")", e)
  .run("B <- ReferenceFrame(\"B\")", e)
  .run("pycall(B, \"orient\", A, \"DCM\", dcm)", e)
  .run("AwB <- pycall(B, \"ang_vel_in\", A)", e)
  .run("alpha2 <- (((s3 * c2) * pycall(q1, \"diff\")) + (c3 * pycall(q2, \"diff\")))", e)
  .run("beta2 <- (((s1 * c2) * pycall(q3, \"diff\")) + (c1 * pycall(q2, \"diff\")))", e)
  expect_equal(.case("pysimplify((pycall(AwB, \"dot\", pyattr(A, \"y\")) - alpha2))", e), "0")
  expect_equal(.case("pysimplify((pycall(AwB, \"dot\", pyattr(B, \"y\")) - beta2))", e), "0")
})

test_that("vector/tests/test_frame.py::test_orient_explicit matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run(".tmp <- pylist(dynamicsymbols(\"c_{xx}, c_{yy}, c_{zz}\")); cxx <- .tmp[[1]]; cyy <- .tmp[[2]]; czz <- .tmp[[3]]", e)
  .run(".tmp <- pylist(dynamicsymbols(\"c_{xy}, c_{xz}, c_{yx}\")); cxy <- .tmp[[1]]; cxz <- .tmp[[2]]; cyx <- .tmp[[3]]", e)
  .run(".tmp <- pylist(dynamicsymbols(\"c_{yz}, c_{zx}, c_{zy}\")); cyz <- .tmp[[1]]; czx <- .tmp[[2]]; czy <- .tmp[[3]]", e)
  .run(".tmp <- pylist(dynamicsymbols(\"c_{xx}, c_{yy}, c_{zz}\", 1)); dcxx <- .tmp[[1]]; dcyy <- .tmp[[2]]; dczz <- .tmp[[3]]", e)
  .run(".tmp <- pylist(dynamicsymbols(\"c_{xy}, c_{xz}, c_{yx}\", 1)); dcxy <- .tmp[[1]]; dcxz <- .tmp[[2]]; dcyx <- .tmp[[3]]", e)
  .run(".tmp <- pylist(dynamicsymbols(\"c_{yz}, c_{zx}, c_{zy}\", 1)); dcyz <- .tmp[[1]]; dczx <- .tmp[[2]]; dczy <- .tmp[[3]]", e)
  .run("A <- ReferenceFrame(\"A\")", e)
  .run("B <- ReferenceFrame(\"B\")", e)
  .run("B_C_A <- Matrix(list(list(cxx, cxy, cxz), list(cyx, cyy, cyz), list(czx, czy, czz)))", e)
  .run("B_w_A <- ((((((cyx * dczx) + (cyy * dczy)) + (cyz * dczz)) * pyattr(B, \"x\")) + ((((czx * dcxx) + (czy * dcxy)) + (czz * dcxz)) * pyattr(B, \"y\"))) + ((((cxx * dcyx) + (cxy * dcyy)) + (cxz * dcyz)) * pyattr(B, \"z\")))", e)
  .run("pycall(A, \"orient_explicit\", B, B_C_A)", e)
  expect_equal(.case("pycall(B, \"dcm\", A)", e), "Matrix([[c_{xx}(t), c_{xy}(t), c_{xz}(t)], [c_{yx}(t), c_{yy}(t), c_{yz}(t)], [c_{zx}(t), c_{zy}(t), c_{zz}(t)]])")
  expect_equal(.case("pycall(A, \"ang_vel_in\", B)", e), "(c_{yx}(t)*Derivative(c_{zx}(t), t) + c_{yy}(t)*Derivative(c_{zy}(t), t) + c_{yz}(t)*Derivative(c_{zz}(t), t))*B.x + (c_{zx}(t)*Derivative(c_{xx}(t), t) + c_{zy}(t)*Derivative(c_{xy}(t), t) + c_{zz}(t)*Derivative(c_{xz}(t), t))*B.y + (c_{xx}(t)*Derivative(c_{yx}(t), t) + c_{xy}(t)*Derivative(c_{yy}(t), t) + c_{xz}(t)*Derivative(c_{yz}(t), t))*B.z")
  expect_equal(.case("pycall(B, \"ang_vel_in\", A)", e), "(-c_{yx}(t)*Derivative(c_{zx}(t), t) - c_{yy}(t)*Derivative(c_{zy}(t), t) - c_{yz}(t)*Derivative(c_{zz}(t), t))*B.x + (-c_{zx}(t)*Derivative(c_{xx}(t), t) - c_{zy}(t)*Derivative(c_{xy}(t), t) - c_{zz}(t)*Derivative(c_{xz}(t), t))*B.y + (-c_{xx}(t)*Derivative(c_{yx}(t), t) - c_{xy}(t)*Derivative(c_{yy}(t), t) - c_{xz}(t)*Derivative(c_{yz}(t), t))*B.z")
})

test_that("vector/tests/test_frame.py::test_orient_dcm matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run(".tmp <- pylist(dynamicsymbols(\"c_{xx}, c_{yy}, c_{zz}\")); cxx <- .tmp[[1]]; cyy <- .tmp[[2]]; czz <- .tmp[[3]]", e)
  .run(".tmp <- pylist(dynamicsymbols(\"c_{xy}, c_{xz}, c_{yx}\")); cxy <- .tmp[[1]]; cxz <- .tmp[[2]]; cyx <- .tmp[[3]]", e)
  .run(".tmp <- pylist(dynamicsymbols(\"c_{yz}, c_{zx}, c_{zy}\")); cyz <- .tmp[[1]]; czx <- .tmp[[2]]; czy <- .tmp[[3]]", e)
  .run("B_C_A <- Matrix(list(list(cxx, cxy, cxz), list(cyx, cyy, cyz), list(czx, czy, czz)))", e)
  .run("A <- ReferenceFrame(\"A\")", e)
  .run("B <- ReferenceFrame(\"B\")", e)
  .run("pycall(B, \"orient_dcm\", A, B_C_A)", e)
  expect_equal(.case("pycall(B, \"dcm\", A)", e), "Matrix([[c_{xx}(t), c_{xy}(t), c_{xz}(t)], [c_{yx}(t), c_{yy}(t), c_{yz}(t)], [c_{zx}(t), c_{zy}(t), c_{zz}(t)]])")
})

test_that("vector/tests/test_frame.py::test_orient_axis matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run("A <- ReferenceFrame(\"A\")", e)
  .run("B <- ReferenceFrame(\"B\")", e)
  .run("pycall(A, \"orient_axis\", B, (-pyattr(B, \"x\")), 1)", e)
  .run("A1 <- pycall(A, \"dcm\", B)", e)
  .run("pycall(A, \"orient_axis\", B, pyattr(B, \"x\"), (-1))", e)
  .run("A2 <- pycall(A, \"dcm\", B)", e)
  .run("pycall(A, \"orient_axis\", B, 1, (-pyattr(B, \"x\")))", e)
  .run("A3 <- pycall(A, \"dcm\", B)", e)
  expect_equal(.case("A1", e), "Matrix([[1, 0, 0], [0, cos(1), -sin(1)], [0, sin(1), cos(1)]])")
  expect_equal(.case("A2", e), "Matrix([[1, 0, 0], [0, cos(1), -sin(1)], [0, sin(1), cos(1)]])")
  expect_equal(.case("pycall(A, \"orient_axis\", B, 1, 1)", e), "ERR TypeError: A Vector must be supplied")
})

test_that("vector/tests/test_frame.py::test_orient_body matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run("A <- ReferenceFrame(\"A\")", e)
  .run("B <- ReferenceFrame(\"B\")", e)
  .run("pycall(B, \"orient_body_fixed\", A, tup(1, 1, 0), \"XYX\")", e)
  expect_equal(.case("pycall(B, \"dcm\", A)", e), "Matrix([[cos(1), sin(1)^2, -sin(1)*cos(1)], [0, cos(1), sin(1)], [sin(1), -sin(1)*cos(1), cos(1)^2]])")
})

test_that("vector/tests/test_frame.py::test_orient_body_advanced matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run(".tmp <- pylist(dynamicsymbols(\"q1:4\")); q1 <- .tmp[[1]]; q2 <- .tmp[[2]]; q3 <- .tmp[[3]]", e)
  .run(".tmp <- pylist(pysymbols(\"c1:4\")); c1 <- .tmp[[1]]; c2 <- .tmp[[2]]; c3 <- .tmp[[3]]", e)
  .run(".tmp <- pylist(dynamicsymbols(\"q1:4\", 1)); u1 <- .tmp[[1]]; u2 <- .tmp[[2]]; u3 <- .tmp[[3]]", e)
  .run(".tmp <- pylist(tup(ReferenceFrame(\"A\"), ReferenceFrame(\"B\"))); A <- .tmp[[1]]; B <- .tmp[[2]]", e)
  .run("pycall(B, \"orient_body_fixed\", A, tup(q1, q2, q3), \"zxy\")", e)
  expect_equal(.case("pycall(A, \"dcm\", B)", e), "Matrix([[-sin(q1(t))*sin(q2(t))*sin(q3(t)) + cos(q1(t))*cos(q3(t)), -sin(q1(t))*cos(q2(t)), sin(q1(t))*sin(q2(t))*cos(q3(t)) + sin(q3(t))*cos(q1(t))], [sin(q1(t))*cos(q3(t)) + sin(q2(t))*sin(q3(t))*cos(q1(t)), cos(q1(t))*cos(q2(t)), sin(q1(t))*sin(q3(t)) - sin(q2(t))*cos(q1(t))*cos(q3(t))], [-sin(q3(t))*cos(q2(t)), sin(q2(t)), cos(q2(t))*cos(q3(t))]])")
  expect_equal(.case("pycall(pycall(B, \"ang_vel_in\", A), \"to_matrix\", B)", e), "Matrix([[-sin(q3(t))*cos(q2(t))*Derivative(q1(t), t) + cos(q3(t))*Derivative(q2(t), t)], [sin(q2(t))*Derivative(q1(t), t) + Derivative(q3(t), t)], [sin(q3(t))*Derivative(q2(t), t) + cos(q2(t))*cos(q3(t))*Derivative(q1(t), t)]])")
  .run(".tmp <- pylist(tup(ReferenceFrame(\"A\"), ReferenceFrame(\"B\"))); A <- .tmp[[1]]; B <- .tmp[[2]]", e)
  .run("pycall(B, \"orient_body_fixed\", A, tup(q1, c2, q3), 131)", e)
  expect_equal(.case("pycall(A, \"dcm\", B)", e), "Matrix([[cos(c2), -sin(c2)*cos(q3(t)), sin(c2)*sin(q3(t))], [sin(c2)*cos(q1(t)), -sin(q1(t))*sin(q3(t)) + cos(c2)*cos(q1(t))*cos(q3(t)), -sin(q1(t))*cos(q3(t)) - sin(q3(t))*cos(c2)*cos(q1(t))], [sin(c2)*sin(q1(t)), sin(q1(t))*cos(c2)*cos(q3(t)) + sin(q3(t))*cos(q1(t)), -sin(q1(t))*sin(q3(t))*cos(c2) + cos(q1(t))*cos(q3(t))]])")
  expect_equal(.case("pycall(pycall(B, \"ang_vel_in\", A), \"to_matrix\", B)", e), "Matrix([[cos(c2)*Derivative(q1(t), t) + Derivative(q3(t), t)], [-sin(c2)*cos(q3(t))*Derivative(q1(t), t)], [sin(c2)*sin(q3(t))*Derivative(q1(t), t)]])")
  .run(".tmp <- pylist(tup(ReferenceFrame(\"A\"), ReferenceFrame(\"B\"))); A <- .tmp[[1]]; B <- .tmp[[2]]", e)
  .run("pycall(B, \"orient_body_fixed\", A, tup(c1, c2, c3), 123)", e)
  expect_equal(.case("pycall(B, \"ang_vel_in\", A)", e), "0")
})

test_that("vector/tests/test_frame.py::test_orient_space_advanced matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run(".tmp <- pylist(dynamicsymbols(\"q1:4\")); q1 <- .tmp[[1]]; q2 <- .tmp[[2]]; q3 <- .tmp[[3]]", e)
  .run(".tmp <- pylist(pysymbols(\"c1:4\")); c1 <- .tmp[[1]]; c2 <- .tmp[[2]]; c3 <- .tmp[[3]]", e)
  .run(".tmp <- pylist(dynamicsymbols(\"q1:4\", 1)); u1 <- .tmp[[1]]; u2 <- .tmp[[2]]; u3 <- .tmp[[3]]", e)
  .run(".tmp <- pylist(tup(ReferenceFrame(\"A\"), ReferenceFrame(\"B\"))); A <- .tmp[[1]]; B <- .tmp[[2]]", e)
  .run("pycall(B, \"orient_space_fixed\", A, tup(q3, q2, q1), \"yxz\")", e)
  expect_equal(.case("pycall(A, \"dcm\", B)", e), "Matrix([[-sin(q1(t))*sin(q2(t))*sin(q3(t)) + cos(q1(t))*cos(q3(t)), -sin(q1(t))*cos(q2(t)), sin(q1(t))*sin(q2(t))*cos(q3(t)) + sin(q3(t))*cos(q1(t))], [sin(q1(t))*cos(q3(t)) + sin(q2(t))*sin(q3(t))*cos(q1(t)), cos(q1(t))*cos(q2(t)), sin(q1(t))*sin(q3(t)) - sin(q2(t))*cos(q1(t))*cos(q3(t))], [-sin(q3(t))*cos(q2(t)), sin(q2(t)), cos(q2(t))*cos(q3(t))]])")
  expect_equal(.case("pycall(pycall(B, \"ang_vel_in\", A), \"to_matrix\", B)", e), "Matrix([[-sin(q3(t))*cos(q2(t))*Derivative(q1(t), t) + cos(q3(t))*Derivative(q2(t), t)], [sin(q2(t))*Derivative(q1(t), t) + Derivative(q3(t), t)], [sin(q3(t))*Derivative(q2(t), t) + cos(q2(t))*cos(q3(t))*Derivative(q1(t), t)]])")
  .run(".tmp <- pylist(tup(ReferenceFrame(\"A\"), ReferenceFrame(\"B\"))); A <- .tmp[[1]]; B <- .tmp[[2]]", e)
  .run("pycall(B, \"orient_space_fixed\", A, tup(q3, c2, q1), 131)", e)
  expect_equal(.case("pycall(A, \"dcm\", B)", e), "Matrix([[cos(c2), -sin(c2)*cos(q3(t)), sin(c2)*sin(q3(t))], [sin(c2)*cos(q1(t)), -sin(q1(t))*sin(q3(t)) + cos(c2)*cos(q1(t))*cos(q3(t)), -sin(q1(t))*cos(q3(t)) - sin(q3(t))*cos(c2)*cos(q1(t))], [sin(c2)*sin(q1(t)), sin(q1(t))*cos(c2)*cos(q3(t)) + sin(q3(t))*cos(q1(t)), -sin(q1(t))*sin(q3(t))*cos(c2) + cos(q1(t))*cos(q3(t))]])")
  expect_equal(.case("pycall(pycall(B, \"ang_vel_in\", A), \"to_matrix\", B)", e), "Matrix([[cos(c2)*Derivative(q1(t), t) + Derivative(q3(t), t)], [-sin(c2)*cos(q3(t))*Derivative(q1(t), t)], [sin(c2)*sin(q3(t))*Derivative(q1(t), t)]])")
  .run(".tmp <- pylist(tup(ReferenceFrame(\"A\"), ReferenceFrame(\"B\"))); A <- .tmp[[1]]; B <- .tmp[[2]]", e)
  .run("pycall(B, \"orient_space_fixed\", A, tup(c1, c2, c3), 123)", e)
  expect_equal(.case("pycall(B, \"ang_vel_in\", A)", e), "0")
})

test_that("vector/tests/test_frame.py::test_orient_body_simple_ang_vel matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run("\"This test ensures that the simplest form of that linear system solution\\\\n    is returned, thus the == for the expression comparison.\"", e)
  .run(".tmp <- pylist(dynamicsymbols(\"psi, theta, varphi\")); psi <- .tmp[[1]]; theta <- .tmp[[2]]; phi <- .tmp[[3]]", e)
  .run("t <- .pv_t()", e)
  .run("A <- ReferenceFrame(\"A\")", e)
  .run("B <- ReferenceFrame(\"B\")", e)
  .run("pycall(B, \"orient_body_fixed\", A, tup(psi, theta, phi), \"ZXZ\")", e)
  .run("A_w_B <- pycall(B, \"ang_vel_in\", A)", e)
  expect_equal(.case("pyget(pyget(pyattr(A_w_B, \"args\"), 0), 1)", e), "B")
  expect_equal(.case("pyget(pyget(pyget(pyattr(A_w_B, \"args\"), 0), 0), 0)", e), "sin(theta(t))*sin(varphi(t))*Derivative(psi(t), t) + cos(varphi(t))*Derivative(theta(t), t)")
  expect_equal(.case("pyget(pyget(pyget(pyattr(A_w_B, \"args\"), 0), 0), 1)", e), "sin(theta(t))*cos(varphi(t))*Derivative(psi(t), t) - sin(varphi(t))*Derivative(theta(t), t)")
  expect_equal(.case("pyget(pyget(pyget(pyattr(A_w_B, \"args\"), 0), 0), 2)", e), "cos(theta(t))*Derivative(psi(t), t) + Derivative(varphi(t), t)")
})

test_that("vector/tests/test_frame.py::test_orient_space matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run("A <- ReferenceFrame(\"A\")", e)
  .run("B <- ReferenceFrame(\"B\")", e)
  .run("pycall(B, \"orient_space_fixed\", A, tup(0, 0, 0), \"123\")", e)
  expect_equal(.case("pycall(B, \"dcm\", A)", e), "Matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])")
})

test_that("vector/tests/test_frame.py::test_orient_quaternion matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run("A <- ReferenceFrame(\"A\")", e)
  .run("B <- ReferenceFrame(\"B\")", e)
  .run("pycall(B, \"orient_quaternion\", A, tup(0, 0, 0, 0))", e)
  expect_equal(.case("pycall(B, \"dcm\", A)", e), "Matrix([[0, 0, 0], [0, 0, 0], [0, 0, 0]])")
})

test_that("vector/tests/test_frame.py::test_frame_dict matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run("A <- ReferenceFrame(\"A\")", e)
  .run("B <- ReferenceFrame(\"B\")", e)
  .run("C <- ReferenceFrame(\"C\")", e)
  .run(".tmp <- pylist(pysymbols(\"a b c\")); a <- .tmp[[1]]; b <- .tmp[[2]]; c <- .tmp[[3]]", e)
  .run("pycall(B, \"orient_axis\", A, pyattr(A, \"x\"), a)", e)
  expect_equal(.case("pyattr(A, \"_dcm_dict\")", e), "{B: Matrix([[1, 0, 0], [0, cos(a), -sin(a)], [0, sin(a), cos(a)]])}")
  expect_equal(.case("pyattr(B, \"_dcm_dict\")", e), "{A: Matrix([[1, 0, 0], [0, cos(a), sin(a)], [0, -sin(a), cos(a)]])}")
  expect_equal(.case("pyattr(C, \"_dcm_dict\")", e), "{}")
  .run("pycall(B, \"orient_axis\", C, pyattr(C, \"x\"), b)", e)
  expect_equal(.case("pyattr(A, \"_dcm_dict\")", e), "{B: Matrix([[1, 0, 0], [0, cos(a), -sin(a)], [0, sin(a), cos(a)]])}")
  expect_equal(.case("pyattr(B, \"_dcm_dict\")", e), "{A: Matrix([[1, 0, 0], [0, cos(a), sin(a)], [0, -sin(a), cos(a)]]), C: Matrix([[1, 0, 0], [0, cos(b), sin(b)], [0, -sin(b), cos(b)]])}")
  expect_equal(.case("pyattr(C, \"_dcm_dict\")", e), "{B: Matrix([[1, 0, 0], [0, cos(b), -sin(b)], [0, sin(b), cos(b)]])}")
  .run("pycall(A, \"orient_axis\", B, pyattr(B, \"x\"), c)", e)
  expect_equal(.case("pyattr(B, \"_dcm_dict\")", e), "{C: Matrix([[1, 0, 0], [0, cos(b), sin(b)], [0, -sin(b), cos(b)]]), A: Matrix([[1, 0, 0], [0, cos(c), -sin(c)], [0, sin(c), cos(c)]])}")
  expect_equal(.case("pyattr(A, \"_dcm_dict\")", e), "{B: Matrix([[1, 0, 0], [0, cos(c), sin(c)], [0, -sin(c), cos(c)]])}")
  expect_equal(.case("pyattr(C, \"_dcm_dict\")", e), "{B: Matrix([[1, 0, 0], [0, cos(b), -sin(b)], [0, sin(b), cos(b)]])}")
})

test_that("vector/tests/test_frame.py::test_dcm_cache_dict matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run("A <- ReferenceFrame(\"A\")", e)
  .run("B <- ReferenceFrame(\"B\")", e)
  .run("C <- ReferenceFrame(\"C\")", e)
  .run("D <- ReferenceFrame(\"D\")", e)
  .run(".tmp <- pylist(pysymbols(\"a b c\")); a <- .tmp[[1]]; b <- .tmp[[2]]; c <- .tmp[[3]]", e)
  .run("pycall(B, \"orient_axis\", A, pyattr(A, \"x\"), a)", e)
  .run("pycall(C, \"orient_axis\", B, pyattr(B, \"x\"), b)", e)
  .run("pycall(D, \"orient_axis\", C, pyattr(C, \"x\"), c)", e)
  expect_equal(.case("pyattr(D, \"_dcm_dict\")", e), "{C: Matrix([[1, 0, 0], [0, cos(c), sin(c)], [0, -sin(c), cos(c)]])}")
  expect_equal(.case("pyattr(C, \"_dcm_dict\")", e), "{B: Matrix([[1, 0, 0], [0, cos(b), sin(b)], [0, -sin(b), cos(b)]]), D: Matrix([[1, 0, 0], [0, cos(c), -sin(c)], [0, sin(c), cos(c)]])}")
  expect_equal(.case("pyattr(B, \"_dcm_dict\")", e), "{A: Matrix([[1, 0, 0], [0, cos(a), sin(a)], [0, -sin(a), cos(a)]]), C: Matrix([[1, 0, 0], [0, cos(b), -sin(b)], [0, sin(b), cos(b)]])}")
  expect_equal(.case("pyattr(A, \"_dcm_dict\")", e), "{B: Matrix([[1, 0, 0], [0, cos(a), -sin(a)], [0, sin(a), cos(a)]])}")
  expect_equal(.case("pyattr(D, \"_dcm_dict\")", e), "{C: Matrix([[1, 0, 0], [0, cos(c), sin(c)], [0, -sin(c), cos(c)]])}")
  .run("pycall(D, \"dcm\", A)", e)
  expect_equal(.case("pylist(pycall(pyattr(A, \"_dcm_cache\"), \"keys\"))", e), "[A, B, D]")
  expect_equal(.case("pylist(pycall(pyattr(D, \"_dcm_cache\"), \"keys\"))", e), "[C, A]")
  expect_equal(.case("pylist(pycall(pyattr(A, \"_dcm_dict\"), \"keys\"))", e), "[B]")
  expect_equal(.case("pylist(pycall(pyattr(D, \"_dcm_dict\"), \"keys\"))", e), "[C]")
  expect_equal(.case("pyattr(A, \"_dcm_dict\")", e), "{B: Matrix([[1, 0, 0], [0, cos(a), -sin(a)], [0, sin(a), cos(a)]])}")
  .run("pycall(A, \"orient_axis\", B, pyattr(B, \"x\"), b)", e)
  expect_equal(.case("pyattr(A, \"_dcm_dict\")", e), "{B: Matrix([[1, 0, 0], [0, cos(b), sin(b)], [0, -sin(b), cos(b)]])}")
  expect_equal(.case("pyattr(A, \"_dcm_dict\")", e), "{B: Matrix([[1, 0, 0], [0, cos(b), sin(b)], [0, -sin(b), cos(b)]])}")
  expect_equal(.case("pyattr(B, \"_dcm_dict\")", e), "{C: Matrix([[1, 0, 0], [0, cos(b), -sin(b)], [0, sin(b), cos(b)]]), A: Matrix([[1, 0, 0], [0, cos(b), -sin(b)], [0, sin(b), cos(b)]])}")
})

test_that("vector/tests/test_frame.py::test_xx_dyad matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run("N <- ReferenceFrame(\"N\")", e)
  .run("`F` <- ReferenceFrame(\"F\", indices = list(\"1\", \"2\", \"3\"))", e)
  expect_equal(.case("pyattr(N, \"xx\")", e), "(N.x|N.x)")
  expect_equal(.case("pyattr(`F`, \"xx\")", e), "(F['1']|F['1'])")
})

test_that("vector/tests/test_frame.py::test_xy_dyad matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run("N <- ReferenceFrame(\"N\")", e)
  .run("`F` <- ReferenceFrame(\"F\", indices = list(\"1\", \"2\", \"3\"))", e)
  expect_equal(.case("pyattr(N, \"xy\")", e), "(N.x|N.y)")
  expect_equal(.case("pyattr(`F`, \"xy\")", e), "(F['1']|F['2'])")
})

test_that("vector/tests/test_frame.py::test_xz_dyad matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run("N <- ReferenceFrame(\"N\")", e)
  .run("`F` <- ReferenceFrame(\"F\", indices = list(\"1\", \"2\", \"3\"))", e)
  expect_equal(.case("pyattr(N, \"xz\")", e), "(N.x|N.z)")
  expect_equal(.case("pyattr(`F`, \"xz\")", e), "(F['1']|F['3'])")
})

test_that("vector/tests/test_frame.py::test_yx_dyad matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run("N <- ReferenceFrame(\"N\")", e)
  .run("`F` <- ReferenceFrame(\"F\", indices = list(\"1\", \"2\", \"3\"))", e)
  expect_equal(.case("pyattr(N, \"yx\")", e), "(N.y|N.x)")
  expect_equal(.case("pyattr(`F`, \"yx\")", e), "(F['2']|F['1'])")
})

test_that("vector/tests/test_frame.py::test_yy_dyad matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run("N <- ReferenceFrame(\"N\")", e)
  .run("`F` <- ReferenceFrame(\"F\", indices = list(\"1\", \"2\", \"3\"))", e)
  expect_equal(.case("pyattr(N, \"yy\")", e), "(N.y|N.y)")
  expect_equal(.case("pyattr(`F`, \"yy\")", e), "(F['2']|F['2'])")
})

test_that("vector/tests/test_frame.py::test_yz_dyad matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run("N <- ReferenceFrame(\"N\")", e)
  .run("`F` <- ReferenceFrame(\"F\", indices = list(\"1\", \"2\", \"3\"))", e)
  expect_equal(.case("pyattr(N, \"yz\")", e), "(N.y|N.z)")
  expect_equal(.case("pyattr(`F`, \"yz\")", e), "(F['2']|F['3'])")
})

test_that("vector/tests/test_frame.py::test_zx_dyad matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run("N <- ReferenceFrame(\"N\")", e)
  .run("`F` <- ReferenceFrame(\"F\", indices = list(\"1\", \"2\", \"3\"))", e)
  expect_equal(.case("pyattr(N, \"zx\")", e), "(N.z|N.x)")
  expect_equal(.case("pyattr(`F`, \"zx\")", e), "(F['3']|F['1'])")
})

test_that("vector/tests/test_frame.py::test_zy_dyad matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run("N <- ReferenceFrame(\"N\")", e)
  .run("`F` <- ReferenceFrame(\"F\", indices = list(\"1\", \"2\", \"3\"))", e)
  expect_equal(.case("pyattr(N, \"zy\")", e), "(N.z|N.y)")
  expect_equal(.case("pyattr(`F`, \"zy\")", e), "(F['3']|F['2'])")
})

test_that("vector/tests/test_frame.py::test_zz_dyad matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run("N <- ReferenceFrame(\"N\")", e)
  .run("`F` <- ReferenceFrame(\"F\", indices = list(\"1\", \"2\", \"3\"))", e)
  expect_equal(.case("pyattr(N, \"zz\")", e), "(N.z|N.z)")
  expect_equal(.case("pyattr(`F`, \"zz\")", e), "(F['3']|F['3'])")
})

test_that("vector/tests/test_frame.py::test_unit_dyadic matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run("N <- ReferenceFrame(\"N\")", e)
  .run("`F` <- ReferenceFrame(\"F\", indices = list(\"1\", \"2\", \"3\"))", e)
  expect_equal(.case("pyattr(N, \"u\")", e), "(N.x|N.x) + (N.y|N.y) + (N.z|N.z)")
  expect_equal(.case("pyattr(`F`, \"u\")", e), "(F['1']|F['1']) + (F['2']|F['2']) + (F['3']|F['3'])")
})

.menv <- new.env(parent = H)
local({ s <- get(".pv_state", envir = ns); s$simp <- FALSE; s$t <- NULL })
try(.run(".tmp <- pylist(pysymbols(\"q1 q2 q3 q4 q5\")); q1 <- .tmp[[1]]; q2 <- .tmp[[2]]; q3 <- .tmp[[3]]; q4 <- .tmp[[4]]; q5 <- .tmp[[5]]; N <- ReferenceFrame(\"N\"); A <- pycall(N, \"orientnew\", \"A\", \"Axis\", list(q1, pyattr(N, \"z\"))); B <- pycall(A, \"orientnew\", \"B\", \"Axis\", list(q2, pyattr(A, \"x\"))); C <- pycall(B, \"orientnew\", \"C\", \"Axis\", list(q3, pyattr(B, \"y\")))", .menv), silent = TRUE)
test_that("vector/tests/test_functions.py::test_dot matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  expect_equal(.case("dot(pyattr(A, \"x\"), pyattr(A, \"x\"))", e), "1")
  expect_equal(.case("dot(pyattr(A, \"x\"), pyattr(A, \"y\"))", e), "0")
  expect_equal(.case("dot(pyattr(A, \"x\"), pyattr(A, \"z\"))", e), "0")
  expect_equal(.case("dot(pyattr(A, \"y\"), pyattr(A, \"x\"))", e), "0")
  expect_equal(.case("dot(pyattr(A, \"y\"), pyattr(A, \"y\"))", e), "1")
  expect_equal(.case("dot(pyattr(A, \"y\"), pyattr(A, \"z\"))", e), "0")
  expect_equal(.case("dot(pyattr(A, \"z\"), pyattr(A, \"x\"))", e), "0")
  expect_equal(.case("dot(pyattr(A, \"z\"), pyattr(A, \"y\"))", e), "0")
  expect_equal(.case("dot(pyattr(A, \"z\"), pyattr(A, \"z\"))", e), "1")
})

test_that("vector/tests/test_functions.py::test_dot_different_frames matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  expect_equal(.case("dot(pyattr(N, \"x\"), pyattr(A, \"x\"))", e), "cos(q1)")
  expect_equal(.case("dot(pyattr(N, \"x\"), pyattr(A, \"y\"))", e), "-sin(q1)")
  expect_equal(.case("dot(pyattr(N, \"x\"), pyattr(A, \"z\"))", e), "0")
  expect_equal(.case("dot(pyattr(N, \"y\"), pyattr(A, \"x\"))", e), "sin(q1)")
  expect_equal(.case("dot(pyattr(N, \"y\"), pyattr(A, \"y\"))", e), "cos(q1)")
  expect_equal(.case("dot(pyattr(N, \"y\"), pyattr(A, \"z\"))", e), "0")
  expect_equal(.case("dot(pyattr(N, \"z\"), pyattr(A, \"x\"))", e), "0")
  expect_equal(.case("dot(pyattr(N, \"z\"), pyattr(A, \"y\"))", e), "0")
  expect_equal(.case("dot(pyattr(N, \"z\"), pyattr(A, \"z\"))", e), "1")
  expect_equal(.case("pytrigsimp(dot(pyattr(N, \"x\"), (pyattr(A, \"x\") + pyattr(A, \"y\"))))", e), "sqrt(2)*cos(q1 + pi/4)")
  expect_equal(.case("pytrigsimp(dot(pyattr(N, \"x\"), (pyattr(A, \"x\") + pyattr(A, \"y\"))))", e), "sqrt(2)*cos(q1 + pi/4)")
  expect_equal(.case("dot(pyattr(A, \"x\"), pyattr(C, \"x\"))", e), "cos(q3)")
  expect_equal(.case("dot(pyattr(A, \"x\"), pyattr(C, \"y\"))", e), "0")
  expect_equal(.case("dot(pyattr(A, \"x\"), pyattr(C, \"z\"))", e), "sin(q3)")
  expect_equal(.case("dot(pyattr(A, \"y\"), pyattr(C, \"x\"))", e), "sin(q2)*sin(q3)")
  expect_equal(.case("dot(pyattr(A, \"y\"), pyattr(C, \"y\"))", e), "cos(q2)")
  expect_equal(.case("dot(pyattr(A, \"y\"), pyattr(C, \"z\"))", e), "-sin(q2)*cos(q3)")
  expect_equal(.case("dot(pyattr(A, \"z\"), pyattr(C, \"x\"))", e), "-sin(q3)*cos(q2)")
  expect_equal(.case("dot(pyattr(A, \"z\"), pyattr(C, \"y\"))", e), "sin(q2)")
  expect_equal(.case("dot(pyattr(A, \"z\"), pyattr(C, \"z\"))", e), "cos(q2)*cos(q3)")
})

test_that("vector/tests/test_functions.py::test_cross matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  expect_equal(.case("cross(pyattr(A, \"x\"), pyattr(A, \"x\"))", e), "0")
  expect_equal(.case("cross(pyattr(A, \"x\"), pyattr(A, \"y\"))", e), "A.z")
  expect_equal(.case("cross(pyattr(A, \"x\"), pyattr(A, \"z\"))", e), "- A.y")
  expect_equal(.case("cross(pyattr(A, \"y\"), pyattr(A, \"x\"))", e), "- A.z")
  expect_equal(.case("cross(pyattr(A, \"y\"), pyattr(A, \"y\"))", e), "0")
  expect_equal(.case("cross(pyattr(A, \"y\"), pyattr(A, \"z\"))", e), "A.x")
  expect_equal(.case("cross(pyattr(A, \"z\"), pyattr(A, \"x\"))", e), "A.y")
  expect_equal(.case("cross(pyattr(A, \"z\"), pyattr(A, \"y\"))", e), "- A.x")
  expect_equal(.case("cross(pyattr(A, \"z\"), pyattr(A, \"z\"))", e), "0")
})

test_that("vector/tests/test_functions.py::test_cross_different_frames matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  expect_equal(.case("cross(pyattr(N, \"x\"), pyattr(A, \"x\"))", e), "sin(q1)*A.z")
  expect_equal(.case("cross(pyattr(N, \"x\"), pyattr(A, \"y\"))", e), "cos(q1)*A.z")
  expect_equal(.case("cross(pyattr(N, \"x\"), pyattr(A, \"z\"))", e), "- sin(q1)*A.x - cos(q1)*A.y")
  expect_equal(.case("cross(pyattr(N, \"y\"), pyattr(A, \"x\"))", e), "- cos(q1)*A.z")
  expect_equal(.case("cross(pyattr(N, \"y\"), pyattr(A, \"y\"))", e), "sin(q1)*A.z")
  expect_equal(.case("cross(pyattr(N, \"y\"), pyattr(A, \"z\"))", e), "cos(q1)*A.x - sin(q1)*A.y")
  expect_equal(.case("cross(pyattr(N, \"z\"), pyattr(A, \"x\"))", e), "A.y")
  expect_equal(.case("cross(pyattr(N, \"z\"), pyattr(A, \"y\"))", e), "- A.x")
  expect_equal(.case("cross(pyattr(N, \"z\"), pyattr(A, \"z\"))", e), "0")
  expect_equal(.case("cross(pyattr(N, \"x\"), pyattr(A, \"x\"))", e), "sin(q1)*A.z")
  expect_equal(.case("cross(pyattr(N, \"x\"), pyattr(A, \"y\"))", e), "cos(q1)*A.z")
  expect_equal(.case("cross(pyattr(N, \"x\"), (pyattr(A, \"x\") + pyattr(A, \"y\")))", e), "(sin(q1) + cos(q1))*A.z")
  expect_equal(.case("cross((pyattr(A, \"x\") + pyattr(A, \"y\")), pyattr(N, \"x\"))", e), "(-sin(q1) - cos(q1))*N.z")
  expect_equal(.case("cross(pyattr(A, \"x\"), pyattr(C, \"x\"))", e), "sin(q3)*C.y")
  expect_equal(.case("cross(pyattr(A, \"x\"), pyattr(C, \"y\"))", e), "- sin(q3)*C.x + cos(q3)*C.z")
  expect_equal(.case("cross(pyattr(A, \"x\"), pyattr(C, \"z\"))", e), "- cos(q3)*C.y")
  expect_equal(.case("cross(pyattr(C, \"x\"), pyattr(A, \"x\"))", e), "- sin(q3)*cos(q2)*A.y - sin(q2)*sin(q3)*A.z")
  expect_equal(.case("pycall(pycall(cross(pyattr(C, \"y\"), pyattr(A, \"x\")), \"express\", C), \"simplify\")", e), "sin(q3)*C.x - cos(q3)*C.z")
  expect_equal(.case("cross(pyattr(C, \"z\"), pyattr(A, \"x\"))", e), "cos(q2)*cos(q3)*A.y + sin(q2)*cos(q3)*A.z")
})

test_that("vector/tests/test_functions.py::test_operator_match matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run("\"Test that the output of dot, cross, outer functions match\\\\n    operator behavior.\\\\n    \"", e)
  .run("A <- ReferenceFrame(\"A\")", e)
  .run("v <- (pyattr(A, \"x\") + pyattr(A, \"y\"))", e)
  .run("d <- (v | v)", e)
  .run("zerov <- PhysVector(0)", e)
  .run("zerod <- PhysDyadic(0)", e)
  expect_equal(.case("(d & d)", e), "2*(A.x|A.x) + 2*(A.x|A.y) + 2*(A.y|A.x) + 2*(A.y|A.y)")
  expect_equal(.case("(d & zerod)", e), "0")
  expect_equal(.case("(zerod & d)", e), "0")
  expect_equal(.case("(d & v)", e), "2*A.x + 2*A.y")
  expect_equal(.case("(v & d)", e), "2*A.x + 2*A.y")
  expect_equal(.case("(d & zerov)", e), "0")
  expect_equal(.case("(zerov & d)", e), "0")
  expect_equal(.case("dot(d, S_Zero)", e), "ERR TypeError: A Vector must be supplied")
  expect_equal(.case("dot(S_Zero, d)", e), "ERR TypeError: Dot product is between two vectors")
  expect_equal(.case("dot(d, 0)", e), "ERR TypeError: A Vector must be supplied")
  expect_equal(.case("dot(0, d)", e), "ERR TypeError: Dot product is between two vectors")
  expect_equal(.case("(v & v)", e), "2")
  expect_equal(.case("(v & zerov)", e), "0")
  expect_equal(.case("(zerov & v)", e), "0")
  expect_equal(.case("dot(v, S_Zero)", e), "ERR TypeError: A Vector must be supplied")
  expect_equal(.case("dot(S_Zero, v)", e), "ERR TypeError: Dot product is between two vectors")
  expect_equal(.case("dot(v, 0)", e), "ERR TypeError: A Vector must be supplied")
  expect_equal(.case("dot(0, v)", e), "ERR TypeError: Dot product is between two vectors")
  expect_equal(.case("cross(d, d)", e), "ERR TypeError: A Vector must be supplied")
  expect_equal(.case("cross(d, zerod)", e), "ERR TypeError: A Vector must be supplied")
  expect_equal(.case("cross(zerod, d)", e), "ERR TypeError: A Vector must be supplied")
  expect_equal(.case("pyxor(d, v)", e), "0")
  expect_equal(.case("pyxor(v, d)", e), "0")
  expect_equal(.case("pyxor(d, zerov)", e), "0")
  expect_equal(.case("pyxor(zerov, d)", e), "0")
  expect_equal(.case("pyxor(zerov, d)", e), "0")
  expect_equal(.case("cross(d, S_Zero)", e), "ERR TypeError: A Vector must be supplied")
  expect_equal(.case("cross(S_Zero, d)", e), "ERR TypeError: Cross product is between two vectors")
  expect_equal(.case("cross(d, 0)", e), "ERR TypeError: A Vector must be supplied")
  expect_equal(.case("cross(0, d)", e), "ERR TypeError: Cross product is between two vectors")
  expect_equal(.case("pyxor(v, v)", e), "0")
  expect_equal(.case("pyxor(v, zerov)", e), "0")
  expect_equal(.case("pyxor(zerov, v)", e), "0")
  expect_equal(.case("cross(v, S_Zero)", e), "ERR TypeError: A Vector must be supplied")
  expect_equal(.case("cross(S_Zero, v)", e), "ERR TypeError: Cross product is between two vectors")
  expect_equal(.case("cross(v, 0)", e), "ERR TypeError: A Vector must be supplied")
  expect_equal(.case("cross(0, v)", e), "ERR TypeError: Cross product is between two vectors")
  expect_equal(.case("outer_sym(d, d)", e), "ERR TypeError: Outer product is between two Vectors")
  expect_equal(.case("outer_sym(d, zerod)", e), "ERR TypeError: Outer product is between two Vectors")
  expect_equal(.case("outer_sym(zerod, d)", e), "ERR TypeError: Outer product is between two Vectors")
  expect_equal(.case("outer_sym(d, v)", e), "ERR TypeError: Outer product is between two Vectors")
  expect_equal(.case("outer_sym(v, d)", e), "ERR TypeError: A Vector must be supplied")
  expect_equal(.case("outer_sym(d, zerov)", e), "ERR TypeError: Outer product is between two Vectors")
  expect_equal(.case("outer_sym(zerov, d)", e), "ERR TypeError: A Vector must be supplied")
  expect_equal(.case("outer_sym(zerov, d)", e), "ERR TypeError: A Vector must be supplied")
  expect_equal(.case("outer_sym(d, S_Zero)", e), "ERR TypeError: Outer product is between two Vectors")
  expect_equal(.case("outer_sym(S_Zero, d)", e), "ERR TypeError: Outer product is between two Vectors")
  expect_equal(.case("outer_sym(d, 0)", e), "ERR TypeError: Outer product is between two Vectors")
  expect_equal(.case("outer_sym(0, d)", e), "ERR TypeError: Outer product is between two Vectors")
  expect_equal(.case("(v | v)", e), "(A.x|A.x) + (A.x|A.y) + (A.y|A.x) + (A.y|A.y)")
  expect_equal(.case("(v | zerov)", e), "0")
  expect_equal(.case("(zerov | v)", e), "0")
  expect_equal(.case("outer_sym(v, S_Zero)", e), "ERR TypeError: A Vector must be supplied")
  expect_equal(.case("outer_sym(S_Zero, v)", e), "ERR TypeError: Outer product is between two Vectors")
  expect_equal(.case("outer_sym(v, 0)", e), "ERR TypeError: A Vector must be supplied")
  expect_equal(.case("outer_sym(0, v)", e), "ERR TypeError: Outer product is between two Vectors")
})

test_that("vector/tests/test_functions.py::test_express matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  expect_equal(.case("express(PhysVector(0), N)", e), "0")
  expect_equal(.case("pyis(express(S_Zero, N), S_Zero)", e), "True")
  expect_equal(.case("express(pyattr(A, \"x\"), C)", e), "cos(q3)*C.x + sin(q3)*C.z")
  expect_equal(.case("express(pyattr(A, \"y\"), C)", e), "sin(q2)*sin(q3)*C.x + cos(q2)*C.y - sin(q2)*cos(q3)*C.z")
  expect_equal(.case("express(pyattr(A, \"z\"), C)", e), "- sin(q3)*cos(q2)*C.x + sin(q2)*C.y + cos(q2)*cos(q3)*C.z")
  expect_equal(.case("express(pyattr(A, \"x\"), N)", e), "cos(q1)*N.x + sin(q1)*N.y")
  expect_equal(.case("express(pyattr(A, \"y\"), N)", e), "- sin(q1)*N.x + cos(q1)*N.y")
  expect_equal(.case("express(pyattr(A, \"z\"), N)", e), "N.z")
  expect_equal(.case("express(pyattr(A, \"x\"), A)", e), "A.x")
  expect_equal(.case("express(pyattr(A, \"y\"), A)", e), "A.y")
  expect_equal(.case("express(pyattr(A, \"z\"), A)", e), "A.z")
  expect_equal(.case("express(pyattr(A, \"x\"), B)", e), "B.x")
  expect_equal(.case("express(pyattr(A, \"y\"), B)", e), "cos(q2)*B.y - sin(q2)*B.z")
  expect_equal(.case("express(pyattr(A, \"z\"), B)", e), "sin(q2)*B.y + cos(q2)*B.z")
  expect_equal(.case("express(pyattr(A, \"x\"), C)", e), "cos(q3)*C.x + sin(q3)*C.z")
  expect_equal(.case("express(pyattr(A, \"y\"), C)", e), "sin(q2)*sin(q3)*C.x + cos(q2)*C.y - sin(q2)*cos(q3)*C.z")
  expect_equal(.case("express(pyattr(A, \"z\"), C)", e), "- sin(q3)*cos(q2)*C.x + sin(q2)*C.y + cos(q2)*cos(q3)*C.z")
  expect_equal(.case("express(pyattr(N, \"x\"), N)", e), "N.x")
  expect_equal(.case("express(pyattr(N, \"y\"), N)", e), "N.y")
  expect_equal(.case("express(pyattr(N, \"z\"), N)", e), "N.z")
  expect_equal(.case("express(pyattr(N, \"x\"), A)", e), "cos(q1)*A.x - sin(q1)*A.y")
  expect_equal(.case("express(pyattr(N, \"y\"), A)", e), "sin(q1)*A.x + cos(q1)*A.y")
  expect_equal(.case("express(pyattr(N, \"z\"), A)", e), "A.z")
  expect_equal(.case("express(pyattr(N, \"x\"), B)", e), "cos(q1)*B.x - sin(q1)*cos(q2)*B.y + sin(q1)*sin(q2)*B.z")
  expect_equal(.case("express(pyattr(N, \"y\"), B)", e), "sin(q1)*B.x + cos(q1)*cos(q2)*B.y - sin(q2)*cos(q1)*B.z")
  expect_equal(.case("express(pyattr(N, \"z\"), B)", e), "sin(q2)*B.y + cos(q2)*B.z")
  expect_equal(.case("express(pyattr(N, \"x\"), C)", e), "(-sin(q1)*sin(q2)*sin(q3) + cos(q1)*cos(q3))*C.x - sin(q1)*cos(q2)*C.y + (sin(q1)*sin(q2)*cos(q3) + sin(q3)*cos(q1))*C.z")
  expect_equal(.case("express(pyattr(N, \"y\"), C)", e), "(sin(q1)*cos(q3) + sin(q2)*sin(q3)*cos(q1))*C.x + cos(q1)*cos(q2)*C.y + (sin(q1)*sin(q3) - sin(q2)*cos(q1)*cos(q3))*C.z")
  expect_equal(.case("express(pyattr(N, \"z\"), C)", e), "- sin(q3)*cos(q2)*C.x + sin(q2)*C.y + cos(q2)*cos(q3)*C.z")
  expect_equal(.case("express(pyattr(A, \"x\"), N)", e), "cos(q1)*N.x + sin(q1)*N.y")
  expect_equal(.case("express(pyattr(A, \"y\"), N)", e), "- sin(q1)*N.x + cos(q1)*N.y")
  expect_equal(.case("express(pyattr(A, \"z\"), N)", e), "N.z")
  expect_equal(.case("express(pyattr(A, \"x\"), A)", e), "A.x")
  expect_equal(.case("express(pyattr(A, \"y\"), A)", e), "A.y")
  expect_equal(.case("express(pyattr(A, \"z\"), A)", e), "A.z")
  expect_equal(.case("express(pyattr(A, \"x\"), B)", e), "B.x")
  expect_equal(.case("express(pyattr(A, \"y\"), B)", e), "cos(q2)*B.y - sin(q2)*B.z")
  expect_equal(.case("express(pyattr(A, \"z\"), B)", e), "sin(q2)*B.y + cos(q2)*B.z")
  expect_equal(.case("express(pyattr(A, \"x\"), C)", e), "cos(q3)*C.x + sin(q3)*C.z")
  expect_equal(.case("express(pyattr(A, \"y\"), C)", e), "sin(q2)*sin(q3)*C.x + cos(q2)*C.y - sin(q2)*cos(q3)*C.z")
  expect_equal(.case("express(pyattr(A, \"z\"), C)", e), "- sin(q3)*cos(q2)*C.x + sin(q2)*C.y + cos(q2)*cos(q3)*C.z")
  expect_equal(.case("express(pyattr(B, \"x\"), N)", e), "cos(q1)*N.x + sin(q1)*N.y")
  expect_equal(.case("express(pyattr(B, \"y\"), N)", e), "- sin(q1)*cos(q2)*N.x + cos(q1)*cos(q2)*N.y + sin(q2)*N.z")
  expect_equal(.case("express(pyattr(B, \"z\"), N)", e), "sin(q1)*sin(q2)*N.x - sin(q2)*cos(q1)*N.y + cos(q2)*N.z")
  expect_equal(.case("express(pyattr(B, \"x\"), A)", e), "A.x")
  expect_equal(.case("express(pyattr(B, \"y\"), A)", e), "cos(q2)*A.y + sin(q2)*A.z")
  expect_equal(.case("express(pyattr(B, \"z\"), A)", e), "- sin(q2)*A.y + cos(q2)*A.z")
  expect_equal(.case("express(pyattr(B, \"x\"), B)", e), "B.x")
  expect_equal(.case("express(pyattr(B, \"y\"), B)", e), "B.y")
  expect_equal(.case("express(pyattr(B, \"z\"), B)", e), "B.z")
  expect_equal(.case("express(pyattr(B, \"x\"), C)", e), "cos(q3)*C.x + sin(q3)*C.z")
  expect_equal(.case("express(pyattr(B, \"y\"), C)", e), "C.y")
  expect_equal(.case("express(pyattr(B, \"z\"), C)", e), "- sin(q3)*C.x + cos(q3)*C.z")
  expect_equal(.case("express(pyattr(C, \"x\"), N)", e), "(-sin(q1)*sin(q2)*sin(q3) + cos(q1)*cos(q3))*N.x + (sin(q1)*cos(q3) + sin(q2)*sin(q3)*cos(q1))*N.y - sin(q3)*cos(q2)*N.z")
  expect_equal(.case("express(pyattr(C, \"y\"), N)", e), "- sin(q1)*cos(q2)*N.x + cos(q1)*cos(q2)*N.y + sin(q2)*N.z")
  expect_equal(.case("express(pyattr(C, \"z\"), N)", e), "(sin(q1)*sin(q2)*cos(q3) + sin(q3)*cos(q1))*N.x + (sin(q1)*sin(q3) - sin(q2)*cos(q1)*cos(q3))*N.y + cos(q2)*cos(q3)*N.z")
  expect_equal(.case("express(pyattr(C, \"x\"), A)", e), "cos(q3)*A.x + sin(q2)*sin(q3)*A.y - sin(q3)*cos(q2)*A.z")
  expect_equal(.case("express(pyattr(C, \"y\"), A)", e), "cos(q2)*A.y + sin(q2)*A.z")
  expect_equal(.case("express(pyattr(C, \"z\"), A)", e), "sin(q3)*A.x - sin(q2)*cos(q3)*A.y + cos(q2)*cos(q3)*A.z")
  expect_equal(.case("express(pyattr(C, \"x\"), B)", e), "cos(q3)*B.x - sin(q3)*B.z")
  expect_equal(.case("express(pyattr(C, \"y\"), B)", e), "B.y")
  expect_equal(.case("express(pyattr(C, \"z\"), B)", e), "sin(q3)*B.x + cos(q3)*B.z")
  expect_equal(.case("express(pyattr(C, \"x\"), C)", e), "C.x")
  expect_equal(.case("express(pyattr(C, \"y\"), C)", e), "C.y")
  expect_equal(.case("express(pyattr(C, \"z\"), C)", e), "C.z")
  expect_equal(.case("pyattr(N, \"x\")", e), "N.x")
  expect_equal(.case("pyattr(N, \"y\")", e), "N.y")
  expect_equal(.case("pyattr(N, \"x\")", e), "N.x")
  expect_equal(.case("pyattr(N, \"y\")", e), "N.y")
  expect_equal(.case("pyattr(N, \"z\")", e), "N.z")
  .run("\"\\\\n    These don't really test our code, they instead test the auto simplification\\\\n    (or lack thereof) of SymPy.\\\\n    assert N.x == express((\\\\n            (cos(q1)*cos(q3)-sin(q1)*sin(q2)*sin(q3))*C.x -\\\\n            sin(q1)*cos(q2)*C.y +\\\\n            (sin(q3)*cos(q1)+sin(q1)*sin(q2)*cos(q3))*C.z), N)\\\\n    assert N.y == express((\\\\n            (sin(q1)*cos(q3) + sin(q2)*sin(q3)*cos(q1))*C.x +\\\\n            cos(q1)*cos(q2)*C.y +\\\\n            (sin(q1)*sin(q3) - sin(q2)*cos(q1)*cos(q3))*C.z), N)\\\\n    assert N.z == express((-sin(q3)*cos(q2)*C.x + sin(q2)*C.y +\\\\n            cos(q2)*cos(q3)*C.z), N)\\\\n    \"", e)
  expect_equal(.case("pyattr(A, \"x\")", e), "A.x")
  expect_equal(.case("pyattr(A, \"y\")", e), "A.y")
  expect_equal(.case("pyattr(A, \"y\")", e), "A.y")
  expect_equal(.case("pyattr(A, \"z\")", e), "A.z")
  expect_equal(.case("pyattr(A, \"x\")", e), "A.x")
  expect_equal(.case("pyattr(A, \"y\")", e), "A.y")
  expect_equal(.case("pyattr(A, \"z\")", e), "A.z")
  expect_equal(.case("pyattr(B, \"x\")", e), "B.x")
  expect_equal(.case("pyattr(B, \"y\")", e), "B.y")
  expect_equal(.case("pyattr(B, \"z\")", e), "B.z")
  expect_equal(.case("pyattr(B, \"y\")", e), "B.y")
  expect_equal(.case("pyattr(B, \"z\")", e), "B.z")
  expect_equal(.case("pyattr(B, \"x\")", e), "B.x")
  expect_equal(.case("pyattr(B, \"z\")", e), "B.z")
  .run("\"\\\\n    assert C.x == express((\\\\n            (cos(q1)*cos(q3)-sin(q1)*sin(q2)*sin(q3))*N.x +\\\\n            (sin(q1)*cos(q3)+sin(q2)*sin(q3)*cos(q1))*N.y -\\\\n                sin(q3)*cos(q2)*N.z), C)\\\\n    assert C.y == express((\\\\n            -sin(q1)*cos(q2)*N.x + cos(q1)*cos(q2)*N.y + sin(q2)*N.z), C)\\\\n    assert C.z == express((\\\\n            (sin(q3)*cos(q1)+sin(q1)*sin(q2)*cos(q3))*N.x +\\\\n            (sin(q1)*sin(q3)-sin(q2)*cos(q1)*cos(q3))*N.y +\\\\n            cos(q2)*cos(q3)*N.z), C)\\\\n    \"", e)
  expect_equal(.case("pyattr(C, \"x\")", e), "C.x")
  expect_equal(.case("pyattr(C, \"y\")", e), "C.y")
  expect_equal(.case("pyattr(C, \"z\")", e), "C.z")
  expect_equal(.case("pyattr(C, \"x\")", e), "C.x")
  expect_equal(.case("pyattr(C, \"z\")", e), "C.z")
})

test_that("vector/tests/test_functions.py::test_time_derivative matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run("A <- ReferenceFrame(\"A\")", e)
  .run("q <- dynamicsymbols(\"q\")", e)
  .run("qd <- dynamicsymbols(\"q\", 1)", e)
  .run("B <- pycall(A, \"orientnew\", \"B\", \"Axis\", list(q, pyattr(A, \"z\")))", e)
  .run("d <- (pyattr(A, \"x\") | pyattr(A, \"x\"))", e)
  expect_equal(.case("time_derivative(d, B)", e), "- Derivative(q(t), t)*(A.y|A.x) - Derivative(q(t), t)*(A.x|A.y)")
  .run("d1 <- (pyattr(A, \"x\") | pyattr(B, \"y\"))", e)
  expect_equal(.case("time_derivative(d1, A)", e), "- Derivative(q(t), t)*(A.x|B.x)")
  expect_equal(.case("time_derivative(d1, B)", e), "- Derivative(q(t), t)*(A.y|B.y)")
  .run("d2 <- (pyattr(A, \"x\") | pyattr(B, \"x\"))", e)
  expect_equal(.case("time_derivative(d2, A)", e), "Derivative(q(t), t)*(A.x|B.y)")
  expect_equal(.case("time_derivative(d2, B)", e), "- Derivative(q(t), t)*(A.y|B.x)")
  .run("d3 <- (pyattr(A, \"x\") | pyattr(B, \"z\"))", e)
  expect_equal(.case("time_derivative(d3, A)", e), "0")
  expect_equal(.case("time_derivative(d3, B)", e), "- Derivative(q(t), t)*(A.y|B.z)")
  .run(".tmp <- pylist(dynamicsymbols(\"q1 q2 q3 q4\")); q1 <- .tmp[[1]]; q2 <- .tmp[[2]]; q3 <- .tmp[[3]]; q4 <- .tmp[[4]]", e)
  .run(".tmp <- pylist(dynamicsymbols(\"q1 q2 q3 q4\", 1)); q1d <- .tmp[[1]]; q2d <- .tmp[[2]]; q3d <- .tmp[[3]]; q4d <- .tmp[[4]]", e)
  .run(".tmp <- pylist(dynamicsymbols(\"q1 q2 q3 q4\", 2)); q1dd <- .tmp[[1]]; q2dd <- .tmp[[2]]; q3dd <- .tmp[[3]]; q4dd <- .tmp[[4]]", e)
  .run("C <- pycall(B, \"orientnew\", \"C\", \"Axis\", list(q4, pyattr(B, \"x\")))", e)
  .run("v1 <- (q1 * pyattr(A, \"z\"))", e)
  .run("v2 <- ((q2 * pyattr(A, \"x\")) + (q3 * pyattr(B, \"y\")))", e)
  .run("v3 <- (((q1 * pyattr(A, \"x\")) + (q2 * pyattr(A, \"y\"))) + (q3 * pyattr(A, \"z\")))", e)
  expect_equal(.case("time_derivative(pyattr(B, \"x\"), C)", e), "0")
  expect_equal(.case("time_derivative(pyattr(B, \"y\"), C)", e), "- Derivative(q4(t), t)*B.z")
  expect_equal(.case("time_derivative(pyattr(B, \"z\"), C)", e), "Derivative(q4(t), t)*B.y")
  expect_equal(.case("time_derivative(v1, B)", e), "Derivative(q1(t), t)*A.z")
  expect_equal(.case("time_derivative(v1, C)", e), "- q1(t)*sin(q(t))*Derivative(q4(t), t)*A.x + q1(t)*cos(q(t))*Derivative(q4(t), t)*A.y + Derivative(q1(t), t)*A.z")
  expect_equal(.case("time_derivative(v2, A)", e), "Derivative(q2(t), t)*A.x - q3(t)*Derivative(q(t), t)*B.x + Derivative(q3(t), t)*B.y")
  expect_equal(.case("time_derivative(v2, C)", e), "Derivative(q2(t), t)*A.x - q2(t)*Derivative(q(t), t)*A.y + q2(t)*sin(q(t))*Derivative(q4(t), t)*A.z + Derivative(q3(t), t)*B.y - q3(t)*Derivative(q4(t), t)*B.z")
  expect_equal(.case("time_derivative(v3, B)", e), "(q2(t)*Derivative(q(t), t) + Derivative(q1(t), t))*A.x + (-q1(t)*Derivative(q(t), t) + Derivative(q2(t), t))*A.y + Derivative(q3(t), t)*A.z")
  expect_equal(.case("time_derivative(d, C)", e), "- Derivative(q(t), t)*(A.y|A.x) + sin(q(t))*Derivative(q4(t), t)*(A.z|A.x) - Derivative(q(t), t)*(A.x|A.y) + sin(q(t))*Derivative(q4(t), t)*(A.x|A.z)")
  expect_equal(.case("time_derivative(pyattr(B, \"x\"), C, order = Float(0.5))", e), "ERR ValueError: Unsupported value of order entered")
  expect_equal(.case("time_derivative(pyattr(B, \"x\"), C, order = (-1))", e), "ERR ValueError: Unsupported value of order entered")
})

test_that("vector/tests/test_functions.py::test_get_motion_methods matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run("t <- .pv_t()", e)
  .run(".tmp <- pylist(pysymbols(\"s1 s2 s3\")); s1 <- .tmp[[1]]; s2 <- .tmp[[2]]; s3 <- .tmp[[3]]", e)
  .run(".tmp <- pylist(pysymbols(\"S1 S2 S3\")); S1 <- .tmp[[1]]; S2 <- .tmp[[2]]; S3 <- .tmp[[3]]", e)
  .run(".tmp <- pylist(pysymbols(\"S4 S5 S6\")); S4 <- .tmp[[1]]; S5 <- .tmp[[2]]; S6 <- .tmp[[3]]", e)
  .run(".tmp <- pylist(pysymbols(\"t1 t2\")); t1 <- .tmp[[1]]; t2 <- .tmp[[2]]", e)
  .run(".tmp <- pylist(dynamicsymbols(\"a b c\")); a <- .tmp[[1]]; b <- .tmp[[2]]; c <- .tmp[[3]]", e)
  .run(".tmp <- pylist(dynamicsymbols(\"a b c\", 1)); ad <- .tmp[[1]]; bd <- .tmp[[2]]; cd <- .tmp[[3]]", e)
  .run(".tmp <- pylist(dynamicsymbols(\"a b c\", 2)); a2d <- .tmp[[1]]; b2d <- .tmp[[2]]; c2d <- .tmp[[3]]", e)
  .run("v0 <- (((S1 * pyattr(N, \"x\")) + (S2 * pyattr(N, \"y\"))) + (S3 * pyattr(N, \"z\")))", e)
  .run("v01 <- (((S4 * pyattr(N, \"x\")) + (S5 * pyattr(N, \"y\"))) + (S6 * pyattr(N, \"z\")))", e)
  .run("v1 <- (((s1 * pyattr(N, \"x\")) + (s2 * pyattr(N, \"y\"))) + (s3 * pyattr(N, \"z\")))", e)
  .run("v2 <- (((a * pyattr(N, \"x\")) + (b * pyattr(N, \"y\"))) + (c * pyattr(N, \"z\")))", e)
  .run("v2d <- (((ad * pyattr(N, \"x\")) + (bd * pyattr(N, \"y\"))) + (cd * pyattr(N, \"z\")))", e)
  .run("v2dd <- (((a2d * pyattr(N, \"x\")) + (b2d * pyattr(N, \"y\"))) + (c2d * pyattr(N, \"z\")))", e)
  expect_equal(.case("get_motion_params(frame = N)", e), "(0, 0, 0)")
  expect_equal(.case("get_motion_params(N, position = v1)", e), "(0, 0, s1*N.x + s2*N.y + s3*N.z)")
  expect_equal(.case("get_motion_params(N, position = v2)", e), "(Derivative(a(t), (t, 2))*N.x + Derivative(b(t), (t, 2))*N.y + Derivative(c(t), (t, 2))*N.z, Derivative(a(t), t)*N.x + Derivative(b(t), t)*N.y + Derivative(c(t), t)*N.z, a(t)*N.x + b(t)*N.y + c(t)*N.z)")
  expect_equal(.case("get_motion_params(N, velocity = v1)", e), "(0, s1*N.x + s2*N.y + s3*N.z, s1*t*N.x + s2*t*N.y + s3*t*N.z)")
  expect_equal(.case("get_motion_params(N, velocity = v1, position = v0, timevalue1 = t1)", e), "(0, s1*N.x + s2*N.y + s3*N.z, (S1 + s1*t - s1*t1)*N.x + (S2 + s2*t - s2*t1)*N.y + (S3 + s3*t - s3*t1)*N.z)")
  .run("answer <- get_motion_params(N, velocity = v1, position = v2, timevalue1 = t1)", e)
  .run("answer_expected <- tup(0, v1, (((v1 * t) - (v1 * t1)) + pycall(v2, \"subs\", t, t1)))", e)
  expect_equal(.case("answer", e), "(0, s1*N.x + s2*N.y + s3*N.z, (s1*t - s1*t1 + a(t1))*N.x + (s2*t - s2*t1 + b(t1))*N.y + (s3*t - s3*t1 + c(t1))*N.z)")
  .run("answer <- get_motion_params(N, velocity = v2, position = v0, timevalue1 = t1)", e)
  .run("integral_vector <- (((Integral(a, tup(t, t1, t)) * pyattr(N, \"x\")) + (Integral(b, tup(t, t1, t)) * pyattr(N, \"y\"))) + (Integral(c, tup(t, t1, t)) * pyattr(N, \"z\")))", e)
  .run("answer_expected <- tup(v2d, v2, (v0 + integral_vector))", e)
  expect_equal(.case("answer", e), "(Derivative(a(t), t)*N.x + Derivative(b(t), t)*N.y + Derivative(c(t), t)*N.z, a(t)*N.x + b(t)*N.y + c(t)*N.z, (S1 + Integral(a(t), (t, t1, t)))*N.x + (S2 + Integral(b(t), (t, t1, t)))*N.y + (S3 + Integral(c(t), (t, t1, t)))*N.z)")
  expect_equal(.case("get_motion_params(N, acceleration = v1)", e), "(s1*N.x + s2*N.y + s3*N.z, s1*t*N.x + s2*t*N.y + s3*t*N.z, s1*t^2/2*N.x + s2*t^2/2*N.y + s3*t^2/2*N.z)")
  expect_equal(.case("get_motion_params(N, acceleration = v1, velocity = v0, position = v2, timevalue1 = t1, timevalue2 = t2)", e), "(s1*N.x + s2*N.y + s3*N.z, (S1 + s1*t - s1*t2)*N.x + (S2 + s2*t - s2*t2)*N.y + (S3 + s3*t - s3*t2)*N.z, (s1*t^2/2 - s1*t1^2/2 + t*(S1 - s1*t2) - t1*(S1 - s1*t2) + a(t1))*N.x + (s2*t^2/2 - s2*t1^2/2 + t*(S2 - s2*t2) - t1*(S2 - s2*t2) + b(t1))*N.y + (s3*t^2/2 - s3*t1^2/2 + t*(S3 - s3*t2) - t1*(S3 - s3*t2) + c(t1))*N.z)")
  expect_equal(.case("get_motion_params(N, acceleration = v1, velocity = v0, position = v01, timevalue1 = t1, timevalue2 = t2)", e), "(s1*N.x + s2*N.y + s3*N.z, (S1 + s1*t - s1*t2)*N.x + (S2 + s2*t - s2*t2)*N.y + (S3 + s3*t - s3*t2)*N.z, (S4 + s1*t^2/2 - s1*t1^2/2 + t*(S1 - s1*t2) - t1*(S1 - s1*t2))*N.x + (S5 + s2*t^2/2 - s2*t1^2/2 + t*(S2 - s2*t2) - t1*(S2 - s2*t2))*N.y + (S6 + s3*t^2/2 - s3*t1^2/2 + t*(S3 - s3*t2) - t1*(S3 - s3*t2))*N.z)")
  .run("answer <- get_motion_params(N, acceleration = (a * pyattr(N, \"x\")), velocity = (S1 * pyattr(N, \"x\")), position = (S2 * pyattr(N, \"x\")), timevalue1 = t1, timevalue2 = t2)", e)
  .run("i1 <- Integral(a, tup(t, t2, t))", e)
  .run("answer_expected <- tup((a * pyattr(N, \"x\")), ((S1 + i1) * pyattr(N, \"x\")), ((S2 + Integral((S1 + i1), tup(t, t1, t))) * pyattr(N, \"x\")))", e)
  expect_equal(.case("answer", e), "(a(t)*N.x, (S1 + Integral(a(t), (t, t2, t)))*N.x, (S2 + Integral(S1 + Integral(a(t), (t, t2, t)), (t, t1, t)))*N.x)")
})

test_that("vector/tests/test_functions.py::test_kin_eqs matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run(".tmp <- pylist(dynamicsymbols(\"q0 q1 q2 q3\")); q0 <- .tmp[[1]]; q1 <- .tmp[[2]]; q2 <- .tmp[[3]]; q3 <- .tmp[[4]]", e)
  .run(".tmp <- pylist(dynamicsymbols(\"q0 q1 q2 q3\", 1)); q0d <- .tmp[[1]]; q1d <- .tmp[[2]]; q2d <- .tmp[[3]]; q3d <- .tmp[[4]]", e)
  .run(".tmp <- pylist(dynamicsymbols(\"u1 u2 u3\")); u1 <- .tmp[[1]]; u2 <- .tmp[[2]]; u3 <- .tmp[[3]]", e)
  .run("ke <- kinematic_equations(list(u1, u2, u3), list(q1, q2, q3), \"body\", 313)", e)
  expect_equal(.case("ke", e), "[-(u1(t)*sin(q3(t)) + u2(t)*cos(q3(t)))/sin(q2(t)) + Derivative(q1(t), t), -u1(t)*cos(q3(t)) + u2(t)*sin(q3(t)) + Derivative(q2(t), t), (u1(t)*sin(q3(t)) + u2(t)*cos(q3(t)))*cos(q2(t))/sin(q2(t)) - u3(t) + Derivative(q3(t), t)]")
  .run("kds <- kinematic_equations(list(u1, u2, u3), list(q0, q1, q2, q3), \"quaternion\")", e)
  expect_equal(.case("kds", e), "[-0.5*q0(t)*u1(t) - 0.5*q2(t)*u3(t) + 0.5*q3(t)*u2(t) + Derivative(q1(t), t), -0.5*q0(t)*u2(t) + 0.5*q1(t)*u3(t) - 0.5*q3(t)*u1(t) + Derivative(q2(t), t), -0.5*q0(t)*u3(t) - 0.5*q1(t)*u2(t) + 0.5*q2(t)*u1(t) + Derivative(q3(t), t), 0.5*q1(t)*u1(t) + 0.5*q2(t)*u2(t) + 0.5*q3(t)*u3(t) + Derivative(q0(t), t)]")
  expect_equal(.case("kinematic_equations(list(u1, u2, u3), list(q0, q1, q2), \"quaternion\")", e), "ERR ValueError: Need 4 coordinates for quaternion")
  expect_equal(.case("kinematic_equations(list(u1, u2, u3), list(q0, q1, q2, q3), \"quaternion\", \"123\")", e), "ERR ValueError: Cannot have rotation order for quaternion")
  expect_equal(.case("kinematic_equations(list(u1, u2, u3), list(q0, q1, q2, q3), \"foo\")", e), "ERR ValueError: Not an approved rotation type for this function")
  expect_equal(.case("kinematic_equations(u1, list(q0, q1, q2, q3), \"quaternion\")", e), "ERR TypeError: Need to supply speeds in a list")
  expect_equal(.case("kinematic_equations(list(u1), list(q0, q1, q2, q3), \"quaternion\")", e), "ERR TypeError: Need to supply 3 body-fixed speeds")
  expect_equal(.case("kinematic_equations(list(u1, u2, u3), q0, \"quaternion\")", e), "ERR TypeError: Need to supply coordinates in a list")
  expect_equal(.case("kinematic_equations(list(u1, u2, u3), list(q0, q1, q2, q3), \"body\")", e), "ERR ValueError: Need 3 coordinates for body or space")
  expect_equal(.case("kinematic_equations(list(u1, u2, u3), list(q0, q1, q2, q3), \"space\")", e), "ERR ValueError: Need 3 coordinates for body or space")
  expect_equal(.case("kinematic_equations(list(u1, u2, u3), list(q0, q1, q2), \"body\", \"222\")", e), "ERR ValueError: Not an acceptable rotation order")
  expect_equal(.case("kinematic_equations(list(0, 0, 0), list(q0, q1, q2), \"space\")", e), "[0, 0, 0]")
})

test_that("vector/tests/test_functions.py::test_partial_velocity matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run(".tmp <- pylist(dynamicsymbols(\"q1 q2 q3 u1 u2 u3\")); q1 <- .tmp[[1]]; q2 <- .tmp[[2]]; q3 <- .tmp[[3]]; u1 <- .tmp[[4]]; u2 <- .tmp[[5]]; u3 <- .tmp[[6]]", e)
  .run(".tmp <- pylist(dynamicsymbols(\"u4, u5\")); u4 <- .tmp[[1]]; u5 <- .tmp[[2]]", e)
  .run("r <- pysymbols(\"r\")", e)
  .run("N <- ReferenceFrame(\"N\")", e)
  .run("Y <- pycall(N, \"orientnew\", \"Y\", \"Axis\", list(q1, pyattr(N, \"z\")))", e)
  .run("L <- pycall(Y, \"orientnew\", \"L\", \"Axis\", list(q2, pyattr(Y, \"x\")))", e)
  .run("R <- pycall(L, \"orientnew\", \"R\", \"Axis\", list(q3, pyattr(L, \"y\")))", e)
  .run("pycall(R, \"set_ang_vel\", N, (((u1 * pyattr(L, \"x\")) + (u2 * pyattr(L, \"y\"))) + (u3 * pyattr(L, \"z\"))))", e)
  .run("C <- PhysPoint(\"C\")", e)
  .run("pycall(C, \"set_vel\", N, ((u4 * pyattr(L, \"x\")) + (u5 * pyxor(pyattr(Y, \"z\"), pyattr(L, \"x\")))))", e)
  .run("Dmc <- pycall(C, \"locatenew\", \"Dmc\", (r * pyattr(L, \"z\")))", e)
  .run("pycall(Dmc, \"v2pt_theory\", C, N, R)", e)
  .run("vel_list <- list(pycall(Dmc, \"vel\", N), pycall(C, \"vel\", N), pycall(R, \"ang_vel_in\", N))", e)
  .run("u_list <- list(u1, u2, u3, u4, u5)", e)
  expect_equal(.case("partial_velocity(vel_list, u_list, N)", e), "[[- r*L.y, r*L.x, 0, L.x, cos(q2(t))*L.y - sin(q2(t))*L.z], [0, 0, 0, L.x, cos(q2(t))*L.y - sin(q2(t))*L.z], [L.x, L.y, L.z, 0, 0]]")
  .run("A <- ReferenceFrame(\"A\")", e)
  .run("B <- ReferenceFrame(\"B\")", e)
  .run("v <- ((u4 * pyattr(A, \"x\")) + (u5 * pyattr(B, \"y\")))", e)
  expect_equal(.case("partial_velocity(tup(v), tup(u4, u5), A)", e), "[[A.x, B.y]]")
  expect_equal(.case("partial_velocity(pycall(Dmc, \"vel\", N), u_list, N)", e), "ERR TypeError: Velocity vectors must be contained in an iterable.")
  expect_equal(.case("partial_velocity(vel_list, u1, N)", e), "ERR TypeError: Generalized speeds must be contained in an iterable")
})

test_that("vector/tests/test_functions.py::test_dynamicsymbols matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run("f1 <- dynamicsymbols(\"f1\")", e)
  .run("f2 <- dynamicsymbols(\"f2\", real = TRUE)", e)
  .run("f3 <- dynamicsymbols(\"f3\", positive = TRUE)", e)
  .run(".tmp <- pylist(dynamicsymbols(\"f4,f5\", commutative = FALSE)); f4 <- .tmp[[1]]; f5 <- .tmp[[2]]", e)
  .run("f6 <- dynamicsymbols(\"f6\", integer = TRUE)", e)
  expect_equal(.case("pyis(pyattr(f1, \"is_real\"), NULL)", e), "True")
  expect_equal(.case("pyattr(f2, \"is_real\")", e), "True")
  expect_equal(.case("pyattr(f3, \"is_positive\")", e), "True")
  expect_equal(.case("(f4 * f5)", e), "f4(t)*f5(t)")
  expect_equal(.case("pyattr(f6, \"is_integer\")", e), "True")
})

.menv <- new.env(parent = H)
local({ s <- get(".pv_state", envir = ns); s$simp <- FALSE; s$t <- NULL })
try(.run("A <- ReferenceFrame(\"A\")", .menv), silent = TRUE)
test_that("vector/tests/test_output.py::test_output_type matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run("A <- ReferenceFrame(\"A\")", e)
  .run("v <- (pyattr(A, \"x\") + pyattr(A, \"y\"))", e)
  .run("d <- (v | v)", e)
  .run("zerov <- PhysVector(0)", e)
  .run("zerod <- PhysDyadic(0)", e)
  expect_equal(.case("pyisinstance((d & d), PhysDyadic)", e), "True")
  expect_equal(.case("pyisinstance((d & zerod), PhysDyadic)", e), "True")
  expect_equal(.case("pyisinstance((zerod & d), PhysDyadic)", e), "True")
  expect_equal(.case("pyisinstance((d & v), PhysVector)", e), "True")
  expect_equal(.case("pyisinstance((v & d), PhysVector)", e), "True")
  expect_equal(.case("pyisinstance((d & zerov), PhysVector)", e), "True")
  expect_equal(.case("pyisinstance((zerov & d), PhysVector)", e), "True")
  expect_equal(.case("(d & S_Zero)", e), "ERR TypeError: A Vector must be supplied")
  expect_equal(.case("(S_Zero & d)", e), "ERR TypeError: unsupported operand type(s) for &: 'Zero' and 'Dyadic'")
  expect_equal(.case("(d & 0)", e), "ERR TypeError: A Vector must be supplied")
  expect_equal(.case("(0 & d)", e), "ERR TypeError: unsupported operand type(s) for &: 'int' and 'Dyadic'")
  expect_equal(.case("(!pybool(pyisinstance((v & v), tup(PhysVector, PhysDyadic))))", e), "True")
  expect_equal(.case("(!pybool(pyisinstance((v & zerov), tup(PhysVector, PhysDyadic))))", e), "True")
  expect_equal(.case("(!pybool(pyisinstance((zerov & v), tup(PhysVector, PhysDyadic))))", e), "True")
  expect_equal(.case("(v & S_Zero)", e), "ERR TypeError: A Vector must be supplied")
  expect_equal(.case("(S_Zero & v)", e), "ERR TypeError: A Vector must be supplied")
  expect_equal(.case("(v & 0)", e), "ERR TypeError: A Vector must be supplied")
  expect_equal(.case("(0 & v)", e), "ERR TypeError: A Vector must be supplied")
  expect_equal(.case("pyxor(d, d)", e), "ERR TypeError: A Vector must be supplied")
  expect_equal(.case("pyxor(d, zerod)", e), "ERR TypeError: A Vector must be supplied")
  expect_equal(.case("pyxor(zerod, d)", e), "ERR TypeError: A Vector must be supplied")
  expect_equal(.case("pyisinstance(pyxor(d, v), PhysDyadic)", e), "True")
  expect_equal(.case("pyisinstance(pyxor(v, d), PhysDyadic)", e), "True")
  expect_equal(.case("pyisinstance(pyxor(d, zerov), PhysDyadic)", e), "True")
  expect_equal(.case("pyisinstance(pyxor(zerov, d), PhysDyadic)", e), "True")
  expect_equal(.case("pyisinstance(pyxor(zerov, d), PhysDyadic)", e), "True")
  expect_equal(.case("pyxor(d, S_Zero)", e), "ERR TypeError: A Vector must be supplied")
  expect_equal(.case("pyxor(S_Zero, d)", e), "ERR TypeError: unsupported operand type(s) for ^: 'Zero' and 'Dyadic'")
  expect_equal(.case("pyxor(d, 0)", e), "ERR TypeError: A Vector must be supplied")
  expect_equal(.case("pyxor(0, d)", e), "ERR TypeError: unsupported operand type(s) for ^: 'int' and 'Dyadic'")
  expect_equal(.case("pyisinstance(pyxor(v, v), PhysVector)", e), "True")
  expect_equal(.case("pyisinstance(pyxor(v, zerov), PhysVector)", e), "True")
  expect_equal(.case("pyisinstance(pyxor(zerov, v), PhysVector)", e), "True")
  expect_equal(.case("pyxor(v, S_Zero)", e), "ERR TypeError: A Vector must be supplied")
  expect_equal(.case("pyxor(S_Zero, v)", e), "ERR TypeError: unsupported operand type(s) for ^: 'Zero' and 'Vector'")
  expect_equal(.case("pyxor(v, 0)", e), "ERR TypeError: A Vector must be supplied")
  expect_equal(.case("pyxor(0, v)", e), "ERR TypeError: unsupported operand type(s) for ^: 'int' and 'Vector'")
  expect_equal(.case("(d | d)", e), "ERR TypeError: unsupported operand type(s) for |: 'Dyadic' and 'Dyadic'")
  expect_equal(.case("(d | zerod)", e), "ERR TypeError: unsupported operand type(s) for |: 'Dyadic' and 'Dyadic'")
  expect_equal(.case("(zerod | d)", e), "ERR TypeError: unsupported operand type(s) for |: 'Dyadic' and 'Dyadic'")
  expect_equal(.case("(d | v)", e), "ERR TypeError: unsupported operand type(s) for |: 'Dyadic' and 'Vector'")
  expect_equal(.case("(v | d)", e), "ERR TypeError: A Vector must be supplied")
  expect_equal(.case("(d | zerov)", e), "ERR TypeError: unsupported operand type(s) for |: 'Dyadic' and 'Vector'")
  expect_equal(.case("(zerov | d)", e), "ERR TypeError: A Vector must be supplied")
  expect_equal(.case("(zerov | d)", e), "ERR TypeError: A Vector must be supplied")
  expect_equal(.case("(d | S_Zero)", e), "ERR TypeError: unsupported operand type(s) for |: 'Dyadic' and 'Zero'")
  expect_equal(.case("(S_Zero | d)", e), "ERR TypeError: unsupported operand type(s) for |: 'Zero' and 'Dyadic'")
  expect_equal(.case("(d | 0)", e), "ERR TypeError: unsupported operand type(s) for |: 'Dyadic' and 'int'")
  expect_equal(.case("(0 | d)", e), "ERR TypeError: unsupported operand type(s) for |: 'int' and 'Dyadic'")
  expect_equal(.case("pyisinstance((v | v), PhysDyadic)", e), "True")
  expect_equal(.case("pyisinstance((v | zerov), PhysDyadic)", e), "True")
  expect_equal(.case("pyisinstance((zerov | v), PhysDyadic)", e), "True")
  expect_equal(.case("(v | S_Zero)", e), "ERR TypeError: A Vector must be supplied")
  expect_equal(.case("(S_Zero | v)", e), "ERR TypeError: unsupported operand type(s) for |: 'Zero' and 'Vector'")
  expect_equal(.case("(v | 0)", e), "ERR TypeError: A Vector must be supplied")
  expect_equal(.case("(0 | v)", e), "ERR TypeError: unsupported operand type(s) for |: 'int' and 'Vector'")
})

.menv <- new.env(parent = H)
local({ s <- get(".pv_state", envir = ns); s$simp <- FALSE; s$t <- NULL })
try(.run("NULL", .menv), silent = TRUE)
test_that("vector/tests/test_point.py::test_point_v1pt_theorys matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run(".tmp <- pylist(dynamicsymbols(\"q q2\")); q <- .tmp[[1]]; q2 <- .tmp[[2]]", e)
  .run(".tmp <- pylist(dynamicsymbols(\"q q2\", 1)); qd <- .tmp[[1]]; q2d <- .tmp[[2]]", e)
  .run(".tmp <- pylist(dynamicsymbols(\"q q2\", 2)); qdd <- .tmp[[1]]; q2dd <- .tmp[[2]]", e)
  .run("N <- ReferenceFrame(\"N\")", e)
  .run("B <- ReferenceFrame(\"B\")", e)
  .run("pycall(B, \"set_ang_vel\", N, (qd * pyattr(B, \"z\")))", e)
  .run("O <- PhysPoint(\"O\")", e)
  .run("P <- pycall(O, \"locatenew\", \"P\", pyattr(B, \"x\"))", e)
  .run("pycall(P, \"set_vel\", B, 0)", e)
  .run("pycall(O, \"set_vel\", N, 0)", e)
  expect_equal(.case("pycall(P, \"v1pt_theory\", O, N, B)", e), "Derivative(q(t), t)*B.y")
  .run("pycall(O, \"set_vel\", N, pyattr(N, \"x\"))", e)
  expect_equal(.case("pycall(P, \"v1pt_theory\", O, N, B)", e), "N.x + Derivative(q(t), t)*B.y")
  .run("pycall(P, \"set_vel\", B, pyattr(B, \"z\"))", e)
  expect_equal(.case("pycall(P, \"v1pt_theory\", O, N, B)", e), "N.x + Derivative(q(t), t)*B.y + B.z")
})

test_that("vector/tests/test_point.py::test_point_a1pt_theorys matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run(".tmp <- pylist(dynamicsymbols(\"q q2\")); q <- .tmp[[1]]; q2 <- .tmp[[2]]", e)
  .run(".tmp <- pylist(dynamicsymbols(\"q q2\", 1)); qd <- .tmp[[1]]; q2d <- .tmp[[2]]", e)
  .run(".tmp <- pylist(dynamicsymbols(\"q q2\", 2)); qdd <- .tmp[[1]]; q2dd <- .tmp[[2]]", e)
  .run("N <- ReferenceFrame(\"N\")", e)
  .run("B <- ReferenceFrame(\"B\")", e)
  .run("pycall(B, \"set_ang_vel\", N, (qd * pyattr(B, \"z\")))", e)
  .run("O <- PhysPoint(\"O\")", e)
  .run("P <- pycall(O, \"locatenew\", \"P\", pyattr(B, \"x\"))", e)
  .run("pycall(P, \"set_vel\", B, 0)", e)
  .run("pycall(O, \"set_vel\", N, 0)", e)
  expect_equal(.case("pycall(P, \"a1pt_theory\", O, N, B)", e), "- Derivative(q(t), t)^2*B.x + Derivative(q(t), (t, 2))*B.y")
  .run("pycall(P, \"set_vel\", B, (q2d * pyattr(B, \"z\")))", e)
  expect_equal(.case("pycall(P, \"a1pt_theory\", O, N, B)", e), "- Derivative(q(t), t)^2*B.x + Derivative(q(t), (t, 2))*B.y + Derivative(q2(t), (t, 2))*B.z")
  .run("pycall(O, \"set_vel\", N, (q2d * pyattr(B, \"x\")))", e)
  expect_equal(.case("pycall(P, \"a1pt_theory\", O, N, B)", e), "(-Derivative(q(t), t)^2 + Derivative(q2(t), (t, 2)))*B.x + (Derivative(q(t), t)*Derivative(q2(t), t) + Derivative(q(t), (t, 2)))*B.y + Derivative(q2(t), (t, 2))*B.z")
})

test_that("vector/tests/test_point.py::test_point_v2pt_theorys matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run("q <- dynamicsymbols(\"q\")", e)
  .run("qd <- dynamicsymbols(\"q\", 1)", e)
  .run("N <- ReferenceFrame(\"N\")", e)
  .run("B <- pycall(N, \"orientnew\", \"B\", \"Axis\", list(q, pyattr(N, \"z\")))", e)
  .run("O <- PhysPoint(\"O\")", e)
  .run("P <- pycall(O, \"locatenew\", \"P\", 0)", e)
  .run("pycall(O, \"set_vel\", N, 0)", e)
  expect_equal(.case("pycall(P, \"v2pt_theory\", O, N, B)", e), "0")
  .run("P <- pycall(O, \"locatenew\", \"P\", pyattr(B, \"x\"))", e)
  expect_equal(.case("pycall(P, \"v2pt_theory\", O, N, B)", e), "Derivative(q(t), t)*B.y")
  .run("pycall(O, \"set_vel\", N, pyattr(N, \"x\"))", e)
  expect_equal(.case("pycall(P, \"v2pt_theory\", O, N, B)", e), "N.x + Derivative(q(t), t)*B.y")
})

test_that("vector/tests/test_point.py::test_point_a2pt_theorys matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run("q <- dynamicsymbols(\"q\")", e)
  .run("qd <- dynamicsymbols(\"q\", 1)", e)
  .run("qdd <- dynamicsymbols(\"q\", 2)", e)
  .run("N <- ReferenceFrame(\"N\")", e)
  .run("B <- pycall(N, \"orientnew\", \"B\", \"Axis\", list(q, pyattr(N, \"z\")))", e)
  .run("O <- PhysPoint(\"O\")", e)
  .run("P <- pycall(O, \"locatenew\", \"P\", 0)", e)
  .run("pycall(O, \"set_vel\", N, 0)", e)
  expect_equal(.case("pycall(P, \"a2pt_theory\", O, N, B)", e), "0")
  .run("pycall(P, \"set_pos\", O, pyattr(B, \"x\"))", e)
  expect_equal(.case("pycall(P, \"a2pt_theory\", O, N, B)", e), "- Derivative(q(t), t)^2*B.x + Derivative(q(t), (t, 2))*B.y")
})

test_that("vector/tests/test_point.py::test_point_funcs matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run(".tmp <- pylist(dynamicsymbols(\"q q2\")); q <- .tmp[[1]]; q2 <- .tmp[[2]]", e)
  .run(".tmp <- pylist(dynamicsymbols(\"q q2\", 1)); qd <- .tmp[[1]]; q2d <- .tmp[[2]]", e)
  .run(".tmp <- pylist(dynamicsymbols(\"q q2\", 2)); qdd <- .tmp[[1]]; q2dd <- .tmp[[2]]", e)
  .run("N <- ReferenceFrame(\"N\")", e)
  .run("B <- ReferenceFrame(\"B\")", e)
  .run("pycall(B, \"set_ang_vel\", N, (5 * pyattr(B, \"y\")))", e)
  .run("O <- PhysPoint(\"O\")", e)
  .run("P <- pycall(O, \"locatenew\", \"P\", ((q * pyattr(B, \"x\")) + (q2 * pyattr(B, \"y\"))))", e)
  expect_equal(.case("pycall(P, \"pos_from\", O)", e), "q(t)*B.x + q2(t)*B.y")
  .run("pycall(P, \"set_vel\", B, ((qd * pyattr(B, \"x\")) + (q2d * pyattr(B, \"y\"))))", e)
  expect_equal(.case("pycall(P, \"vel\", B)", e), "Derivative(q(t), t)*B.x + Derivative(q2(t), t)*B.y")
  .run("pycall(O, \"set_vel\", N, 0)", e)
  expect_equal(.case("pycall(O, \"vel\", N)", e), "0")
  expect_equal(.case("pycall(P, \"a1pt_theory\", O, N, B)", e), "(-25*q(t) + Derivative(q(t), (t, 2)))*B.x + Derivative(q2(t), (t, 2))*B.y - 10*Derivative(q(t), t)*B.z")
  .run("B <- pycall(N, \"orientnew\", \"B\", \"Axis\", list(q, pyattr(N, \"z\")))", e)
  .run("O <- PhysPoint(\"O\")", e)
  .run("P <- pycall(O, \"locatenew\", \"P\", (10 * pyattr(B, \"x\")))", e)
  .run("pycall(O, \"set_vel\", N, (5 * pyattr(N, \"x\")))", e)
  expect_equal(.case("pycall(O, \"vel\", N)", e), "5*N.x")
  expect_equal(.case("pycall(P, \"a2pt_theory\", O, N, B)", e), "- 10*Derivative(q(t), t)^2*B.x + 10*Derivative(q(t), (t, 2))*B.y")
  .run("pycall(B, \"set_ang_vel\", N, (5 * pyattr(B, \"y\")))", e)
  .run("O <- PhysPoint(\"O\")", e)
  .run("P <- pycall(O, \"locatenew\", \"P\", ((q * pyattr(B, \"x\")) + (q2 * pyattr(B, \"y\"))))", e)
  .run("pycall(P, \"set_vel\", B, ((qd * pyattr(B, \"x\")) + (q2d * pyattr(B, \"y\"))))", e)
  .run("pycall(O, \"set_vel\", N, 0)", e)
  expect_equal(.case("pycall(P, \"v1pt_theory\", O, N, B)", e), "Derivative(q(t), t)*B.x + Derivative(q2(t), t)*B.y - 5*q(t)*B.z")
})

test_that("vector/tests/test_point.py::test_point_pos matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run("q <- dynamicsymbols(\"q\")", e)
  .run("N <- ReferenceFrame(\"N\")", e)
  .run("B <- pycall(N, \"orientnew\", \"B\", \"Axis\", list(q, pyattr(N, \"z\")))", e)
  .run("O <- PhysPoint(\"O\")", e)
  .run("P <- pycall(O, \"locatenew\", \"P\", ((10 * pyattr(N, \"x\")) + (5 * pyattr(B, \"x\"))))", e)
  expect_equal(.case("pycall(P, \"pos_from\", O)", e), "10*N.x + 5*B.x")
  .run("Q <- pycall(P, \"locatenew\", \"Q\", ((10 * pyattr(N, \"y\")) + (5 * pyattr(B, \"y\"))))", e)
  expect_equal(.case("pycall(Q, \"pos_from\", P)", e), "10*N.y + 5*B.y")
  expect_equal(.case("pycall(Q, \"pos_from\", O)", e), "10*N.x + 10*N.y + 5*B.x + 5*B.y")
  expect_equal(.case("pycall(O, \"pos_from\", Q)", e), "- 10*N.x - 10*N.y - 5*B.x - 5*B.y")
})

test_that("vector/tests/test_point.py::test_point_partial_velocity matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run("N <- ReferenceFrame(\"N\")", e)
  .run("A <- ReferenceFrame(\"A\")", e)
  .run("p <- PhysPoint(\"p\")", e)
  .run(".tmp <- pylist(dynamicsymbols(\"u1, u2\")); u1 <- .tmp[[1]]; u2 <- .tmp[[2]]", e)
  .run("pycall(p, \"set_vel\", N, ((u1 * pyattr(A, \"x\")) + (u2 * pyattr(N, \"y\"))))", e)
  expect_equal(.case("pycall(p, \"partial_velocity\", N, u1)", e), "A.x")
  expect_equal(.case("pycall(p, \"partial_velocity\", N, u1, u2)", e), "(A.x, N.y)")
  expect_equal(.case("pycall(p, \"partial_velocity\", A, u1)", e), "ERR ValueError: \\nVelocity of point p has not been defined in ReferenceFrame A.")
})

test_that("vector/tests/test_point.py::test_point_vel matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run(".tmp <- pylist(dynamicsymbols(\"q1 q2\")); q1 <- .tmp[[1]]; q2 <- .tmp[[2]]", e)
  .run("N <- ReferenceFrame(\"N\")", e)
  .run("B <- ReferenceFrame(\"B\")", e)
  .run("Q <- PhysPoint(\"Q\")", e)
  .run("O <- PhysPoint(\"O\")", e)
  .run("pycall(Q, \"set_pos\", O, (q1 * pyattr(N, \"x\")))", e)
  expect_equal(.case("pycall(Q, \"vel\", N)", e), "ERR ValueError: \\nVelocity of point Q has not been defined in ReferenceFrame N.")
  .run("pycall(O, \"set_vel\", N, (q2 * pyattr(N, \"y\")))", e)
  expect_equal(.case("pycall(O, \"vel\", N)", e), "q2(t)*N.y")
  expect_equal(.case("pycall(O, \"vel\", B)", e), "ERR ValueError: \\nVelocity of point O has not been defined in ReferenceFrame B.")
})

test_that("vector/tests/test_point.py::test_auto_point_vel matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run("t <- .pv_t()", e)
  .run(".tmp <- pylist(dynamicsymbols(\"q1 q2\")); q1 <- .tmp[[1]]; q2 <- .tmp[[2]]", e)
  .run("N <- ReferenceFrame(\"N\")", e)
  .run("B <- ReferenceFrame(\"B\")", e)
  .run("O <- PhysPoint(\"O\")", e)
  .run("Q <- PhysPoint(\"Q\")", e)
  .run("pycall(Q, \"set_pos\", O, (q1 * pyattr(N, \"x\")))", e)
  .run("pycall(O, \"set_vel\", N, (q2 * pyattr(N, \"y\")))", e)
  expect_equal(.case("pycall(Q, \"vel\", N)", e), "Derivative(q1(t), t)*N.x + q2(t)*N.y")
  .run("P1 <- PhysPoint(\"P1\")", e)
  .run("pycall(P1, \"set_pos\", O, (q1 * pyattr(B, \"x\")))", e)
  .run("P2 <- PhysPoint(\"P2\")", e)
  .run("pycall(P2, \"set_pos\", P1, (q2 * pyattr(B, \"z\")))", e)
  expect_equal(.case("pycall(P2, \"vel\", B)", e), "ERR ValueError: \\nVelocity of point P2 has not been defined in ReferenceFrame B.")
  expect_equal(.case("pycall(P2, \"vel\", N)", e), "ERR ValueError: \\nVelocity of point P2 has not been defined in ReferenceFrame N.")
})

test_that("vector/tests/test_point.py::test_auto_point_vel_multiple_point_path matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run("t <- .pv_t()", e)
  .run(".tmp <- pylist(dynamicsymbols(\"q1 q2\")); q1 <- .tmp[[1]]; q2 <- .tmp[[2]]", e)
  .run("B <- ReferenceFrame(\"B\")", e)
  .run("P <- PhysPoint(\"P\")", e)
  .run("pycall(P, \"set_vel\", B, (q1 * pyattr(B, \"x\")))", e)
  .run("P1 <- PhysPoint(\"P1\")", e)
  .run("pycall(P1, \"set_pos\", P, (q2 * pyattr(B, \"y\")))", e)
  .run("pycall(P1, \"set_vel\", B, (q1 * pyattr(B, \"z\")))", e)
  .run("P2 <- PhysPoint(\"P2\")", e)
  .run("pycall(P2, \"set_pos\", P1, (q1 * pyattr(B, \"z\")))", e)
  .run("P3 <- PhysPoint(\"P3\")", e)
  .run("pycall(P3, \"set_pos\", P2, ((10 * q1) * pyattr(B, \"y\")))", e)
  expect_equal(.case("pycall(P3, \"vel\", B)", e), "10*Derivative(q1(t), t)*B.y + (q1(t) + Derivative(q1(t), t))*B.z")
})

test_that("vector/tests/test_point.py::test_auto_vel_dont_overwrite matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run("t <- .pv_t()", e)
  .run(".tmp <- pylist(dynamicsymbols(\"q1, q2, u1\")); q1 <- .tmp[[1]]; q2 <- .tmp[[2]]; u1 <- .tmp[[3]]", e)
  .run("N <- ReferenceFrame(\"N\")", e)
  .run("P <- PhysPoint(\"P1\")", e)
  .run("pycall(P, \"set_vel\", N, (u1 * pyattr(N, \"x\")))", e)
  .run("P1 <- PhysPoint(\"P1\")", e)
  .run("pycall(P1, \"set_pos\", P, (q2 * pyattr(N, \"y\")))", e)
  expect_equal(.case("pycall(P1, \"vel\", N)", e), "u1(t)*N.x + Derivative(q2(t), t)*N.y")
  expect_equal(.case("pycall(P, \"vel\", N)", e), "u1(t)*N.x")
  .run("pycall(P1, \"set_vel\", N, (u1 * pyattr(N, \"z\")))", e)
  expect_equal(.case("pycall(P1, \"vel\", N)", e), "u1(t)*N.z")
})

test_that("vector/tests/test_point.py::test_auto_point_vel_if_tree_has_vel_but_inappropriate_pos_vector matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run(".tmp <- pylist(dynamicsymbols(\"q1 q2\")); q1 <- .tmp[[1]]; q2 <- .tmp[[2]]", e)
  .run("B <- ReferenceFrame(\"B\")", e)
  .run("S <- ReferenceFrame(\"S\")", e)
  .run("P <- PhysPoint(\"P\")", e)
  .run("pycall(P, \"set_vel\", B, (q1 * pyattr(B, \"x\")))", e)
  .run("P1 <- PhysPoint(\"P1\")", e)
  .run("pycall(P1, \"set_pos\", P, pyattr(S, \"y\"))", e)
  expect_equal(.case("pycall(P1, \"vel\", B)", e), "ERR ValueError: \\nVelocity of point P1 has not been defined in ReferenceFrame B.")
  expect_equal(.case("pycall(P1, \"vel\", S)", e), "ERR ValueError: \\nVelocity of point P1 has not been defined in ReferenceFrame S.")
})

test_that("vector/tests/test_point.py::test_auto_point_vel_shortest_path matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run("t <- .pv_t()", e)
  .run(".tmp <- pylist(dynamicsymbols(\"q1 q2 u1 u2\")); q1 <- .tmp[[1]]; q2 <- .tmp[[2]]; u1 <- .tmp[[3]]; u2 <- .tmp[[4]]", e)
  .run("B <- ReferenceFrame(\"B\")", e)
  .run("P <- PhysPoint(\"P\")", e)
  .run("pycall(P, \"set_vel\", B, (u1 * pyattr(B, \"x\")))", e)
  .run("P1 <- PhysPoint(\"P1\")", e)
  .run("pycall(P1, \"set_pos\", P, (q2 * pyattr(B, \"y\")))", e)
  .run("pycall(P1, \"set_vel\", B, (q1 * pyattr(B, \"z\")))", e)
  .run("P2 <- PhysPoint(\"P2\")", e)
  .run("pycall(P2, \"set_pos\", P1, (q1 * pyattr(B, \"z\")))", e)
  .run("P3 <- PhysPoint(\"P3\")", e)
  .run("pycall(P3, \"set_pos\", P2, ((10 * q1) * pyattr(B, \"y\")))", e)
  .run("P4 <- PhysPoint(\"P4\")", e)
  .run("pycall(P4, \"set_pos\", P3, (q1 * pyattr(B, \"x\")))", e)
  .run("O <- PhysPoint(\"O\")", e)
  .run("pycall(O, \"set_vel\", B, (u2 * pyattr(B, \"y\")))", e)
  .run("O1 <- PhysPoint(\"O1\")", e)
  .run("pycall(O1, \"set_pos\", O, (q2 * pyattr(B, \"z\")))", e)
  .run("pycall(P4, \"set_pos\", O1, ((q1 * pyattr(B, \"x\")) + (q2 * pyattr(B, \"z\"))))", e)
  .run("pycall(warnings, \"simplefilter\", \"error\")", e)
  expect_equal(.case("pycall(P4, \"vel\", B)", e), "ERR UserWarning: \\nVelocity of P4 automatically calculated based on point O but it is\\nalso possible from points(s): [P1]. Velocities from these points are\\nnot necessarily the same. This may cause errors in your calculations.")
})

test_that("vector/tests/test_point.py::test_auto_point_vel_connected_frames matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run("t <- .pv_t()", e)
  .run(".tmp <- pylist(dynamicsymbols(\"q q1 q2 u\")); q <- .tmp[[1]]; q1 <- .tmp[[2]]; q2 <- .tmp[[3]]; u <- .tmp[[4]]", e)
  .run("N <- ReferenceFrame(\"N\")", e)
  .run("B <- ReferenceFrame(\"B\")", e)
  .run("O <- PhysPoint(\"O\")", e)
  .run("pycall(O, \"set_vel\", N, (u * pyattr(N, \"x\")))", e)
  .run("P <- PhysPoint(\"P\")", e)
  .run("pycall(P, \"set_pos\", O, ((q1 * pyattr(N, \"x\")) + (q2 * pyattr(B, \"y\"))))", e)
  expect_equal(.case("pycall(P, \"vel\", N)", e), "ERR ValueError: \\nVelocity of point P has not been defined in ReferenceFrame N.")
  .run("pycall(N, \"orient\", B, \"Axis\", tup(q, pyattr(B, \"x\")))", e)
  expect_equal(.case("pycall(P, \"vel\", N)", e), "(u(t) + Derivative(q1(t), t))*N.x + Derivative(q2(t), t)*B.y - q2(t)*Derivative(q(t), t)*B.z")
})

test_that("vector/tests/test_point.py::test_auto_point_vel_multiple_paths_warning_arises matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run(".tmp <- pylist(dynamicsymbols(\"q u\")); q <- .tmp[[1]]; u <- .tmp[[2]]", e)
  .run("N <- ReferenceFrame(\"N\")", e)
  .run("O <- PhysPoint(\"O\")", e)
  .run("P <- PhysPoint(\"P\")", e)
  .run("Q <- PhysPoint(\"Q\")", e)
  .run("R <- PhysPoint(\"R\")", e)
  .run("pycall(P, \"set_vel\", N, (u * pyattr(N, \"x\")))", e)
  .run("pycall(Q, \"set_vel\", N, (u * pyattr(N, \"y\")))", e)
  .run("pycall(R, \"set_vel\", N, (u * pyattr(N, \"z\")))", e)
  .run("pycall(O, \"set_pos\", P, (q * pyattr(N, \"z\")))", e)
  .run("pycall(O, \"set_pos\", Q, (q * pyattr(N, \"y\")))", e)
  .run("pycall(O, \"set_pos\", R, (q * pyattr(N, \"x\")))", e)
  .run("pycall(warnings, \"simplefilter\", \"error\")", e)
  expect_equal(.case("pycall(O, \"vel\", N)", e), "ERR UserWarning: \\nVelocity of O automatically calculated based on point P but it is also\\npossible from points(s): [Q, R]. Velocities from these points are not\\nnecessarily the same. This may cause errors in your calculations.")
})

test_that("vector/tests/test_point.py::test_auto_vel_cyclic_warning_arises matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run("P <- PhysPoint(\"P\")", e)
  .run("P1 <- PhysPoint(\"P1\")", e)
  .run("P2 <- PhysPoint(\"P2\")", e)
  .run("P3 <- PhysPoint(\"P3\")", e)
  .run("N <- ReferenceFrame(\"N\")", e)
  .run("pycall(P, \"set_vel\", N, pyattr(N, \"x\"))", e)
  .run("pycall(P1, \"set_pos\", P, pyattr(N, \"x\"))", e)
  .run("pycall(P2, \"set_pos\", P1, pyattr(N, \"y\"))", e)
  .run("pycall(P3, \"set_pos\", P2, pyattr(N, \"z\"))", e)
  .run("pycall(P1, \"set_pos\", P3, (pyattr(N, \"x\") + pyattr(N, \"y\")))", e)
  .run("pycall(warnings, \"simplefilter\", \"error\")", e)
  expect_equal(.case("pycall(P2, \"vel\", N)", e), "ERR UserWarning: \\nKinematic loops are defined among the positions of points. This is\\nlikely not desired and may cause errors in your calculations.")
})

test_that("vector/tests/test_point.py::test_auto_vel_derivative matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run(".tmp <- pylist(dynamicsymbols(\"q1:3\")); q1 <- .tmp[[1]]; q2 <- .tmp[[2]]", e)
  .run(".tmp <- pylist(dynamicsymbols(\"u1:3\", 1)); u1 <- .tmp[[1]]; u2 <- .tmp[[2]]", e)
  .run("A <- ReferenceFrame(\"A\")", e)
  .run("B <- ReferenceFrame(\"B\")", e)
  .run("C <- ReferenceFrame(\"C\")", e)
  .run("pycall(B, \"orient_axis\", A, pyattr(A, \"z\"), q1)", e)
  .run("pycall(B, \"set_ang_vel\", A, (u1 * pyattr(A, \"z\")))", e)
  .run("pycall(C, \"orient_axis\", B, pyattr(B, \"z\"), q2)", e)
  .run("pycall(C, \"set_ang_vel\", B, (u2 * pyattr(B, \"z\")))", e)
  .run("Am <- PhysPoint(\"Am\")", e)
  .run("pycall(Am, \"set_vel\", A, 0)", e)
  .run("Bm <- PhysPoint(\"Bm\")", e)
  .run("pycall(Bm, \"set_pos\", Am, pyattr(B, \"x\"))", e)
  .run("pycall(Bm, \"set_vel\", B, 0)", e)
  .run("pycall(Bm, \"set_vel\", C, 0)", e)
  .run("Cm <- PhysPoint(\"Cm\")", e)
  .run("pycall(Cm, \"set_pos\", Bm, pyattr(C, \"x\"))", e)
  .run("pycall(Cm, \"set_vel\", C, 0)", e)
  .run("temp <- pycall(pyattr(Cm, \"_vel_dict\"), \"copy\")", e)
  expect_equal(.case("pycall(Cm, \"vel\", A)", e), "Derivative(u1(t), t)*B.y + (Derivative(u1(t), t) + Derivative(u2(t), t))*C.y")
  .run("pysetattr(Cm, \"_vel_dict\", temp)", e)
  .run("pycall(Cm, \"v2pt_theory\", Bm, B, C)", e)
  expect_equal(.case("pycall(Cm, \"vel\", A)", e), "Derivative(u1(t), t)*B.y + (Derivative(u1(t), t) + Derivative(u2(t), t))*C.y")
})

test_that("vector/tests/test_point.py::test_auto_point_acc_zero_vel matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run("N <- ReferenceFrame(\"N\")", e)
  .run("O <- PhysPoint(\"O\")", e)
  .run("pycall(O, \"set_vel\", N, 0)", e)
  expect_equal(.case("pycall(O, \"acc\", N)", e), "0")
})

test_that("vector/tests/test_point.py::test_auto_point_acc_compute_vel matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run("t <- .pv_t()", e)
  .run("q1 <- dynamicsymbols(\"q1\")", e)
  .run("N <- ReferenceFrame(\"N\")", e)
  .run("A <- ReferenceFrame(\"A\")", e)
  .run("pycall(A, \"orient_axis\", N, pyattr(N, \"z\"), q1)", e)
  .run("O <- PhysPoint(\"O\")", e)
  .run("pycall(O, \"set_vel\", N, 0)", e)
  .run("P <- PhysPoint(\"P\")", e)
  .run("pycall(P, \"set_pos\", O, pyattr(A, \"x\"))", e)
  expect_equal(.case("pycall(P, \"acc\", N)", e), "- Derivative(q1(t), t)^2*A.x + Derivative(q1(t), (t, 2))*A.y")
})

test_that("vector/tests/test_point.py::test_auto_acc_derivative matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run(".tmp <- pylist(dynamicsymbols(\"q1:3\")); q1 <- .tmp[[1]]; q2 <- .tmp[[2]]", e)
  .run(".tmp <- pylist(dynamicsymbols(\"q1:3\", 1)); u1 <- .tmp[[1]]; u2 <- .tmp[[2]]", e)
  .run(".tmp <- pylist(dynamicsymbols(\"q1:3\", 2)); v1 <- .tmp[[1]]; v2 <- .tmp[[2]]", e)
  .run("A <- ReferenceFrame(\"A\")", e)
  .run("B <- ReferenceFrame(\"B\")", e)
  .run("C <- ReferenceFrame(\"C\")", e)
  .run("pycall(B, \"orient_axis\", A, pyattr(A, \"z\"), q1)", e)
  .run("pycall(C, \"orient_axis\", B, pyattr(B, \"z\"), q2)", e)
  .run("Am <- PhysPoint(\"Am\")", e)
  .run("pycall(Am, \"set_vel\", A, 0)", e)
  .run("Bm <- PhysPoint(\"Bm\")", e)
  .run("pycall(Bm, \"set_pos\", Am, pyattr(B, \"x\"))", e)
  .run("pycall(Bm, \"set_vel\", B, 0)", e)
  .run("pycall(Bm, \"set_vel\", C, 0)", e)
  .run("Cm <- PhysPoint(\"Cm\")", e)
  .run("pycall(Cm, \"set_pos\", Bm, pyattr(C, \"x\"))", e)
  .run("pycall(Cm, \"set_vel\", C, 0)", e)
  .run(".tmp <- pylist(tup(pycall(pyattr(Bm, \"_vel_dict\"), \"copy\"), pycall(pyattr(Cm, \"_vel_dict\"), \"copy\"))); Bm_vel_dict <- .tmp[[1]]; Cm_vel_dict <- .tmp[[2]]", e)
  .run(".tmp <- pylist(tup(pycall(pyattr(Bm, \"_acc_dict\"), \"copy\"), pycall(pyattr(Cm, \"_acc_dict\"), \"copy\"))); Bm_acc_dict <- .tmp[[1]]; Cm_acc_dict <- .tmp[[2]]", e)
  .run("check <- (((((-pypow(u1, 2)) * pyattr(B, \"x\")) + (v1 * pyattr(B, \"y\"))) - (pypow((u1 + u2), 2) * pyattr(C, \"x\"))) + ((v1 + v2) * pyattr(C, \"y\")))", e)
  expect_equal(.case("pycall(Cm, \"acc\", A)", e), "- Derivative(q1(t), t)^2*B.x + Derivative(q1(t), (t, 2))*B.y - (Derivative(q1(t), t) + Derivative(q2(t), t))^2*C.x + (Derivative(q1(t), (t, 2)) + Derivative(q2(t), (t, 2)))*C.y")
  .run(".tmp <- pylist(tup(Bm_vel_dict, Cm_vel_dict)); pysetattr(Bm, \"_vel_dict\", .tmp[[1]]); pysetattr(Cm, \"_vel_dict\", .tmp[[2]])", e)
  .run(".tmp <- pylist(tup(Bm_acc_dict, Cm_acc_dict)); pysetattr(Bm, \"_acc_dict\", .tmp[[1]]); pysetattr(Cm, \"_acc_dict\", .tmp[[2]])", e)
  .run("pycall(Bm, \"v2pt_theory\", Am, A, B)", e)
  .run("pycall(Cm, \"v2pt_theory\", Bm, A, C)", e)
  .run("pycall(Bm, \"a2pt_theory\", Am, A, B)", e)
  expect_equal(.case("pycall(Cm, \"a2pt_theory\", Bm, A, C)", e), "- Derivative(q1(t), t)^2*B.x + Derivative(q1(t), (t, 2))*B.y - (Derivative(q1(t), t) + Derivative(q2(t), t))^2*C.x + (Derivative(q1(t), (t, 2)) + Derivative(q2(t), (t, 2)))*C.y")
})

.menv <- new.env(parent = H)
local({ s <- get(".pv_state", envir = ns); s$simp <- FALSE; s$t <- NULL })
try(.run(".tmp <- pylist(pysymbols(\"a, b, c\")); a <- .tmp[[1]]; b <- .tmp[[2]]; c <- .tmp[[3]]; .tmp <- pylist(dynamicsymbols(\"alpha, omega, beta\")); alpha <- .tmp[[1]]; omega <- .tmp[[2]]; beta <- .tmp[[3]]; A <- ReferenceFrame(\"A\"); N <- ReferenceFrame(\"N\"); v <- (((pypow(a, 2) * pyattr(N, \"x\")) + (b * pyattr(N, \"y\"))) + ((c * sin(alpha)) * pyattr(N, \"z\"))); w <- (((alpha * pyattr(N, \"x\")) + (sin(omega) * pyattr(N, \"y\"))) + ((alpha * beta) * pyattr(N, \"z\"))); ww <- (((alpha * pyattr(N, \"x\")) + (asin(omega) * pyattr(N, \"y\"))) - ((pycall(alpha, \"diff\") * beta) * pyattr(N, \"z\"))); o <- ((((a / b) * pyattr(N, \"x\")) + (((c + b) / a) * pyattr(N, \"y\"))) + ((pypow(c, 2) / b) * pyattr(N, \"z\"))); y <- (((pypow(a, 2) * (pyattr(N, \"x\") | pyattr(N, \"y\"))) + (b * (pyattr(N, \"y\") | pyattr(N, \"y\")))) + ((c * sin(alpha)) * (pyattr(N, \"z\") | pyattr(N, \"y\")))); x <- (((alpha * (pyattr(N, \"x\") | pyattr(N, \"x\"))) + (sin(omega) * (pyattr(N, \"y\") | pyattr(N, \"z\")))) + ((alpha * beta) * (pyattr(N, \"z\") | pyattr(N, \"x\")))); xx <- (pyattr(N, \"x\") | ((-pyattr(N, \"y\")) - pyattr(N, \"z\"))); xx2 <- (pyattr(N, \"x\") | (pyattr(N, \"y\") + pyattr(N, \"z\"))); ascii_vpretty <- function(expr) { return(vpprint(expr, use_unicode = FALSE, wrap_line = FALSE)) }; unicode_vpretty <- function(expr) { return(vpprint(expr, use_unicode = TRUE, wrap_line = FALSE)) }", .menv), silent = TRUE)
test_that("vector/tests/test_printing.py::test_latex_printer matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run("r <- (Function(\"r\"))(\"t\")", e)
  expect_equal(.case("pycall(VectorLatexPrinter(), \"doprint\", pypow(r, 2))", e), "r^{2}")
  .run("r2 <- (Function(\"r^2\"))(\"t\")", e)
  expect_equal(.case("pycall(VectorLatexPrinter(), \"doprint\", pycall(r2, \"diff\"))", e), "\\dot{r^{2}}")
  .run("ra <- (Function(\"r__a\"))(\"t\")", e)
  expect_equal(.case("pycall(VectorLatexPrinter(), \"doprint\", pycall(pycall(ra, \"diff\"), \"diff\"))", e), "\\ddot{r^{a}}")
})

test_that("vector/tests/test_printing.py::test_vector_pretty_print matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run("expected <- \" 2                               \\\\na  n_x + b n_y + c*sin(alpha) n_z\"", e)
  .run("uexpected <- \" 2                           \\\\na  n_x + b n_y + c\u{22c5}sin(\u{03b1}) n_z\"", e)
  expect_equal(.case("ascii_vpretty(v)", e), " 2                               \\na  n_x + b n_y + c*sin(alpha) n_z")
  .run("expected <- \"alpha n_x + sin(omega) n_y + alpha*beta n_z\"", e)
  .run("uexpected <- \"\u{03b1} n_x + sin(\u{03c9}) n_y + \u{03b1}\u{22c5}\u{03b2} n_z\"", e)
  expect_equal(.case("ascii_vpretty(w)", e), "alpha n_x + sin(omega) n_y + alpha*beta n_z")
  .run("expected <- \"                     2    \\\\na       b + c       c     \\\\n- n_x + ----- n_y + -- n_z\\\\nb         a         b     \"", e)
  .run("uexpected <- \"                     2    \\\\na       b + c       c     \\\\n\u{2500} n_x + \u{2500}\u{2500}\u{2500}\u{2500}\u{2500} n_y + \u{2500}\u{2500} n_z\\\\nb         a         b     \"", e)
  expect_equal(.case("ascii_vpretty(o)", e), "                     2    \\na       b + c       c     \\n- n_x + ----- n_y + -- n_z\\nb         a         b     ")
  expect_equal(.case("unicode_vpretty(o)", e), "                     2    \\na       b + c       c     \\n\u{2500} n_x + \u{2500}\u{2500}\u{2500}\u{2500}\u{2500} n_y + \u{2500}\u{2500} n_z\\nb         a         b     ")
  expect_equal(.case("ascii_vpretty((-pyattr(A, \"x\")))", e), "-a_x")
  expect_equal(.case("unicode_vpretty((-pyattr(A, \"x\")))", e), "-a_x")
  expect_equal(.case("ascii_vpretty((0 * pyattr(A, \"x\")))", e), "0")
  expect_equal(.case("unicode_vpretty((0 * pyattr(A, \"x\")))", e), "0")
})

test_that("vector/tests/test_printing.py::test_vector_latex matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run(".tmp <- pylist(pysymbols(\"a, b, c, d, omega\")); a <- .tmp[[1]]; b <- .tmp[[2]]; c <- .tmp[[3]]; d <- .tmp[[4]]; omega <- .tmp[[5]]", e)
  .run("v <- ((((pypow(a, 2) + (b / c)) * pyattr(A, \"x\")) + (sqrt(d) * pyattr(A, \"y\"))) + (cos(omega) * pyattr(A, \"z\")))", e)
  expect_equal(.case("vlatex(v)", e), "(a^{2} + \\frac{b}{c})\\mathbf{\\hat{a}_x} + \\sqrt{d}\\mathbf{\\hat{a}_y} + \\cos{\\left(\\omega \\right)}\\mathbf{\\hat{a}_z}")
  .run(".tmp <- pylist(dynamicsymbols(\"theta, omega, alpha, q\")); theta <- .tmp[[1]]; omega <- .tmp[[2]]; alpha <- .tmp[[3]]; q <- .tmp[[4]]", e)
  .run("v <- (((theta * pyattr(A, \"x\")) + ((omega * omega) * pyattr(A, \"y\"))) + ((q * alpha) * pyattr(A, \"z\")))", e)
  expect_equal(.case("vlatex(v)", e), "\\theta\\mathbf{\\hat{a}_x} + \\omega^{2}\\mathbf{\\hat{a}_y} + \\alpha q\\mathbf{\\hat{a}_z}")
  .run(".tmp <- pylist(dynamicsymbols(\"phi1, phi2, phi3\")); phi1 <- .tmp[[1]]; phi2 <- .tmp[[2]]; phi3 <- .tmp[[3]]", e)
  .run(".tmp <- pylist(pysymbols(\"theta1, theta2, theta3\")); theta1 <- .tmp[[1]]; theta2 <- .tmp[[2]]; theta3 <- .tmp[[3]]", e)
  .run("v <- (((sin(theta1) * pyattr(A, \"x\")) + ((cos(phi1) * cos(phi2)) * pyattr(A, \"y\"))) + (cos((theta1 + phi3)) * pyattr(A, \"z\")))", e)
  expect_equal(.case("vlatex(v)", e), "\\sin{\\left(\\theta_{1} \\right)}\\mathbf{\\hat{a}_x} + \\cos{\\left(\\phi_{1} \\right)} \\cos{\\left(\\phi_{2} \\right)}\\mathbf{\\hat{a}_y} + \\cos{\\left(\\theta_{1} + \\phi_{3} \\right)}\\mathbf{\\hat{a}_z}")
  .run("N <- ReferenceFrame(\"N\")", e)
  .run(".tmp <- pylist(pysymbols(\"a, b, c, d, omega\")); a <- .tmp[[1]]; b <- .tmp[[2]]; c <- .tmp[[3]]; d <- .tmp[[4]]; omega <- .tmp[[5]]", e)
  .run("v <- ((((pypow(a, 2) + (b / c)) * pyattr(N, \"x\")) + (sqrt(d) * pyattr(N, \"y\"))) + (cos(omega) * pyattr(N, \"z\")))", e)
  .run("expected <- \"(a^{2} + \\\\\\\\frac{b}{c})\\\\\\\\mathbf{\\\\\\\\hat{n}_x} + \\\\\\\\sqrt{d}\\\\\\\\mathbf{\\\\\\\\hat{n}_y} + \\\\\\\\cos{\\\\\\\\left(\\\\\\\\omega \\\\\\\\right)}\\\\\\\\mathbf{\\\\\\\\hat{n}_z}\"", e)
  expect_equal(.case("vlatex(v)", e), "(a^{2} + \\frac{b}{c})\\mathbf{\\hat{n}_x} + \\sqrt{d}\\mathbf{\\hat{n}_y} + \\cos{\\left(\\omega \\right)}\\mathbf{\\hat{n}_z}")
  .run("N <- ReferenceFrame(\"N\", latexs = tup(\"\\\\\\\\hat{i}\", \"\\\\\\\\hat{j}\", \"\\\\\\\\hat{k}\"))", e)
  .run("v <- ((((pypow(a, 2) + (b / c)) * pyattr(N, \"x\")) + (sqrt(d) * pyattr(N, \"y\"))) + (cos(omega) * pyattr(N, \"z\")))", e)
  .run("expected <- \"(a^{2} + \\\\\\\\frac{b}{c})\\\\\\\\hat{i} + \\\\\\\\sqrt{d}\\\\\\\\hat{j} + \\\\\\\\cos{\\\\\\\\left(\\\\\\\\omega \\\\\\\\right)}\\\\\\\\hat{k}\"", e)
  expect_equal(.case("vlatex(v)", e), "(a^{2} + \\frac{b}{c})\\hat{i} + \\sqrt{d}\\hat{j} + \\cos{\\left(\\omega \\right)}\\hat{k}")
  .run("expected <- \"\\\\\\\\alpha\\\\\\\\mathbf{\\\\\\\\hat{n}_x} + \\\\\\\\operatorname{asin}{\\\\\\\\left(\\\\\\\\omega \\\\\\\\right)}\\\\\\\\mathbf{\\\\\\\\hat{n}_y} -  \\\\\\\\beta \\\\\\\\dot{\\\\\\\\alpha}\\\\\\\\mathbf{\\\\\\\\hat{n}_z}\"", e)
  expect_equal(.case("vlatex(ww)", e), "\\alpha\\mathbf{\\hat{n}_x} + \\operatorname{asin}{\\left(\\omega \\right)}\\mathbf{\\hat{n}_y} -  \\beta \\dot{\\alpha}\\mathbf{\\hat{n}_z}")
  .run("expected <- \"- \\\\\\\\mathbf{\\\\\\\\hat{n}_x}\\\\\\\\otimes \\\\\\\\mathbf{\\\\\\\\hat{n}_y} - \\\\\\\\mathbf{\\\\\\\\hat{n}_x}\\\\\\\\otimes \\\\\\\\mathbf{\\\\\\\\hat{n}_z}\"", e)
  expect_equal(.case("vlatex(xx)", e), "- \\mathbf{\\hat{n}_x}\\otimes \\mathbf{\\hat{n}_y} - \\mathbf{\\hat{n}_x}\\otimes \\mathbf{\\hat{n}_z}")
  .run("expected <- \"\\\\\\\\mathbf{\\\\\\\\hat{n}_x}\\\\\\\\otimes \\\\\\\\mathbf{\\\\\\\\hat{n}_y} + \\\\\\\\mathbf{\\\\\\\\hat{n}_x}\\\\\\\\otimes \\\\\\\\mathbf{\\\\\\\\hat{n}_z}\"", e)
  expect_equal(.case("vlatex(xx2)", e), "\\mathbf{\\hat{n}_x}\\otimes \\mathbf{\\hat{n}_y} + \\mathbf{\\hat{n}_x}\\otimes \\mathbf{\\hat{n}_z}")
})

test_that("vector/tests/test_printing.py::test_vector_latex_arguments matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  expect_equal(.case("vlatex((pyattr(N, \"x\") * Float(3.0)), full_prec = FALSE)", e), "3.0\\mathbf{\\hat{n}_x}")
})

test_that("vector/tests/test_printing.py::test_vector_latex_with_functions matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run("N <- ReferenceFrame(\"N\")", e)
  .run(".tmp <- pylist(dynamicsymbols(\"omega, alpha\")); omega <- .tmp[[1]]; alpha <- .tmp[[2]]", e)
  .run("v <- (pycall(omega, \"diff\") * pyattr(N, \"x\"))", e)
  expect_equal(.case("vlatex(v)", e), "\\dot{\\omega}\\mathbf{\\hat{n}_x}")
  .run("v <- (pypow(pycall(omega, \"diff\"), alpha) * pyattr(N, \"x\"))", e)
  expect_equal(.case("vlatex(v)", e), "\\dot{\\omega}^{\\alpha}\\mathbf{\\hat{n}_x}")
})

test_that("vector/tests/test_printing.py::test_dyadic_latex matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run("expected <- \"a^{2}\\\\\\\\mathbf{\\\\\\\\hat{n}_x}\\\\\\\\otimes \\\\\\\\mathbf{\\\\\\\\hat{n}_y} + b\\\\\\\\mathbf{\\\\\\\\hat{n}_y}\\\\\\\\otimes \\\\\\\\mathbf{\\\\\\\\hat{n}_y} + c \\\\\\\\sin{\\\\\\\\left(\\\\\\\\alpha \\\\\\\\right)}\\\\\\\\mathbf{\\\\\\\\hat{n}_z}\\\\\\\\otimes \\\\\\\\mathbf{\\\\\\\\hat{n}_y}\"", e)
  expect_equal(.case("vlatex(y)", e), "a^{2}\\mathbf{\\hat{n}_x}\\otimes \\mathbf{\\hat{n}_y} + b\\mathbf{\\hat{n}_y}\\otimes \\mathbf{\\hat{n}_y} + c \\sin{\\left(\\alpha \\right)}\\mathbf{\\hat{n}_z}\\otimes \\mathbf{\\hat{n}_y}")
  .run("expected <- \"\\\\\\\\alpha\\\\\\\\mathbf{\\\\\\\\hat{n}_x}\\\\\\\\otimes \\\\\\\\mathbf{\\\\\\\\hat{n}_x} + \\\\\\\\sin{\\\\\\\\left(\\\\\\\\omega \\\\\\\\right)}\\\\\\\\mathbf{\\\\\\\\hat{n}_y}\\\\\\\\otimes \\\\\\\\mathbf{\\\\\\\\hat{n}_z} + \\\\\\\\alpha \\\\\\\\beta\\\\\\\\mathbf{\\\\\\\\hat{n}_z}\\\\\\\\otimes \\\\\\\\mathbf{\\\\\\\\hat{n}_x}\"", e)
  expect_equal(.case("vlatex(x)", e), "\\alpha\\mathbf{\\hat{n}_x}\\otimes \\mathbf{\\hat{n}_x} + \\sin{\\left(\\omega \\right)}\\mathbf{\\hat{n}_y}\\otimes \\mathbf{\\hat{n}_z} + \\alpha \\beta\\mathbf{\\hat{n}_z}\\otimes \\mathbf{\\hat{n}_x}")
  expect_equal(.case("vlatex(PhysDyadic(list()))", e), "0")
})

test_that("vector/tests/test_printing.py::test_dyadic_str matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  expect_equal(.case("vsprint(PhysDyadic(list()))", e), "0")
  expect_equal(.case("vsprint(y)", e), "a^2*(N.x|N.y) + b*(N.y|N.y) + c*sin(alpha)*(N.z|N.y)")
  expect_equal(.case("vsprint(x)", e), "alpha*(N.x|N.x) + sin(omega)*(N.y|N.z) + alpha*beta*(N.z|N.x)")
  expect_equal(.case("vsprint(ww)", e), "alpha*N.x + asin(omega)*N.y - beta*alpha'*N.z")
  expect_equal(.case("vsprint(xx)", e), "- (N.x|N.y) - (N.x|N.z)")
  expect_equal(.case("vsprint(xx2)", e), "(N.x|N.y) + (N.x|N.z)")
})

test_that("vector/tests/test_printing.py::test_vlatex matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run("NULL", e)
  .run("x <- pysymbols(\"x\")", e)
  .run("J <- pysymbols(\"J\")", e)
  .run("f <- Function(\"f\")", e)
  .run("g <- Function(\"g\")", e)
  .run("h <- Function(\"h\")", e)
  .run("expected <- \"J \\\\\\\\left(\\\\\\\\frac{d}{d x} g{\\\\\\\\left(x \\\\\\\\right)} - \\\\\\\\frac{d}{d x} h{\\\\\\\\left(x \\\\\\\\right)}\\\\\\\\right)\"", e)
  .run("expr <- (J * pycall(pycall(f(x), \"diff\", x), \"subs\", f(x), (g(x) - h(x))))", e)
  expect_equal(.case("vlatex(expr)", e), "J \\left(\\frac{d}{d x} g{\\left(x \\right)} - \\frac{d}{d x} h{\\left(x \\right)}\\right)")
})

test_that("vector/tests/test_printing.py::test_issue_13354 matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run("\"\\\\n    Test for proper pretty printing of physics vectors with ADD\\\\n    instances in arguments.\\\\n\\\\n    Test is exactly the one suggested in the original bug report by\\\\n    @moorepants.\\\\n    \"", e)
  .run(".tmp <- pylist(pysymbols(\"a, b, c\")); a <- .tmp[[1]]; b <- .tmp[[2]]; c <- .tmp[[3]]", e)
  .run("A <- ReferenceFrame(\"A\")", e)
  .run("v <- (((a * pyattr(A, \"x\")) + (b * pyattr(A, \"y\"))) + (c * pyattr(A, \"z\")))", e)
  .run("w <- (((b * pyattr(A, \"x\")) + (c * pyattr(A, \"y\"))) + (a * pyattr(A, \"z\")))", e)
  .run("z <- (w + v)", e)
  .run("expected <- \"(a + b) a_x + (b + c) a_y + (a + c) a_z\"", e)
  expect_equal(.case("ascii_vpretty(z)", e), "(a + b) a_x + (b + c) a_y + (a + c) a_z")
})

test_that("vector/tests/test_printing.py::test_vector_derivative_printing matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run("v <- (pycall(omega, \"diff\") * pyattr(N, \"x\"))", e)
  expect_equal(.case("ascii_vpretty(v)", e), "omega'(t) n_x")
  .run("v <- (pycall(pycall(omega, \"diff\"), \"diff\") * pyattr(N, \"x\"))", e)
  expect_equal(.case("vlatex(v)", e), "\\ddot{\\omega}\\mathbf{\\hat{n}_x}")
  expect_equal(.case("ascii_vpretty(v)", e), "omega''(t) n_x")
  .run("v <- (pycall(pycall(pycall(omega, \"diff\"), \"diff\"), \"diff\") * pyattr(N, \"x\"))", e)
  expect_equal(.case("vlatex(v)", e), "\\dddot{\\omega}\\mathbf{\\hat{n}_x}")
  expect_equal(.case("ascii_vpretty(v)", e), "omega'''(t) n_x")
  .run("v <- (pycall(pycall(pycall(pycall(omega, \"diff\"), \"diff\"), \"diff\"), \"diff\") * pyattr(N, \"x\"))", e)
  expect_equal(.case("vlatex(v)", e), "\\ddddot{\\omega}\\mathbf{\\hat{n}_x}")
  expect_equal(.case("ascii_vpretty(v)", e), "omega''''(t) n_x")
  .run("v <- (pycall(pycall(pycall(pycall(pycall(omega, \"diff\"), \"diff\"), \"diff\"), \"diff\"), \"diff\") * pyattr(N, \"x\"))", e)
  expect_equal(.case("vlatex(v)", e), "\\frac{d^{5}}{d t^{5}} \\omega\\mathbf{\\hat{n}_x}")
})

test_that("vector/tests/test_printing.py::test_vector_str_printing matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  expect_equal(.case("vsprint(w)", e), "alpha*N.x + sin(omega)*N.y + alpha*beta*N.z")
  expect_equal(.case("vsprint((pycall(omega, \"diff\") * pyattr(N, \"x\")))", e), "omega'*N.x")
  expect_equal(.case("vsstrrepr(w)", e), "alpha*N.x + sin(omega)*N.y + alpha*beta*N.z")
})

test_that("vector/tests/test_printing.py::test_vector_str_arguments matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  expect_equal(.case("vsprint((pyattr(N, \"x\") * Float(3.0)), full_prec = FALSE)", e), "3.0*N.x")
})

.menv <- new.env(parent = H)
local({ s <- get(".pv_state", envir = ns); s$simp <- FALSE; s$t <- NULL })
try(.run("x <- Symbol(\"x\"); y <- Symbol(\"y\"); z <- Symbol(\"z\"); A <- ReferenceFrame(\"A\")", .menv), silent = TRUE)
test_that("vector/tests/test_vector.py::test_free_dynamicsymbols matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run(".tmp <- pylist(pysymbols(\"A, B, C, D\", cls = ReferenceFrame)); A <- .tmp[[1]]; B <- .tmp[[2]]; C <- .tmp[[3]]; D <- .tmp[[4]]", e)
  .run(".tmp <- pylist(dynamicsymbols(\"a, b, c, d, e, f\")); a <- .tmp[[1]]; b <- .tmp[[2]]; c <- .tmp[[3]]; d <- .tmp[[4]]; e <- .tmp[[5]]; f <- .tmp[[6]]", e)
  .run("pycall(B, \"orient_axis\", A, a, pyattr(A, \"x\"))", e)
  .run("pycall(C, \"orient_axis\", B, b, pyattr(B, \"y\"))", e)
  .run("pycall(D, \"orient_axis\", C, c, pyattr(C, \"x\"))", e)
  .run("v <- (((d * pyattr(D, \"x\")) + (e * pyattr(D, \"y\"))) + (f * pyattr(D, \"z\")))", e)
  expect_equal(.case("pyset(pyordered(pycall(v, \"free_dynamicsymbols\", A)))", e), "{a(t), b(t), c(t), d(t), e(t), f(t)}")
  expect_equal(.case("pyset(pyordered(pycall(v, \"free_dynamicsymbols\", B)))", e), "{b(t), c(t), d(t), e(t), f(t)}")
  expect_equal(.case("pyset(pyordered(pycall(v, \"free_dynamicsymbols\", C)))", e), "{c(t), d(t), e(t), f(t)}")
  expect_equal(.case("pyset(pyordered(pycall(v, \"free_dynamicsymbols\", D)))", e), "{d(t), e(t), f(t)}")
})

test_that("vector/tests/test_vector.py::test_Vector matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  expect_equal(.case("pyattr(A, \"x\")", e), "A.x")
  expect_equal(.case("pyattr(A, \"y\")", e), "A.y")
  expect_equal(.case("pyattr(A, \"z\")", e), "A.z")
  expect_equal(.case("(pyattr(A, \"x\") + 0)", e), "A.x")
  .run("v1 <- (((x * pyattr(A, \"x\")) + (y * pyattr(A, \"y\"))) + (z * pyattr(A, \"z\")))", e)
  .run("v2 <- (((pypow(x, 2) * pyattr(A, \"x\")) + (pypow(y, 2) * pyattr(A, \"y\"))) + (pypow(z, 2) * pyattr(A, \"z\")))", e)
  .run("v3 <- (v1 + v2)", e)
  .run("v4 <- (v1 - v2)", e)
  expect_equal(.case("pyisinstance(v1, PhysVector)", e), "True")
  expect_equal(.case("dot(v1, pyattr(A, \"x\"))", e), "x")
  expect_equal(.case("dot(v1, pyattr(A, \"y\"))", e), "y")
  expect_equal(.case("dot(v1, pyattr(A, \"z\"))", e), "z")
  expect_equal(.case("pyisinstance(v2, PhysVector)", e), "True")
  expect_equal(.case("dot(v2, pyattr(A, \"x\"))", e), "x^2")
  expect_equal(.case("dot(v2, pyattr(A, \"y\"))", e), "y^2")
  expect_equal(.case("dot(v2, pyattr(A, \"z\"))", e), "z^2")
  expect_equal(.case("pyisinstance(v3, PhysVector)", e), "True")
  expect_equal(.case("dot(v3, pyattr(A, \"x\"))", e), "x^2 + x")
  expect_equal(.case("dot(v3, pyattr(A, \"y\"))", e), "y^2 + y")
  expect_equal(.case("dot(v3, pyattr(A, \"z\"))", e), "z^2 + z")
  expect_equal(.case("pyisinstance(v4, PhysVector)", e), "True")
  expect_equal(.case("dot(v4, pyattr(A, \"x\"))", e), "-x^2 + x")
  expect_equal(.case("dot(v4, pyattr(A, \"y\"))", e), "-y^2 + y")
  expect_equal(.case("dot(v4, pyattr(A, \"z\"))", e), "-z^2 + z")
  expect_equal(.case("pycall(v1, \"to_matrix\", A)", e), "Matrix([[x], [y], [z]])")
  .run("q <- pysymbols(\"q\")", e)
  .run("B <- pycall(A, \"orientnew\", \"B\", \"Axis\", tup(q, pyattr(A, \"x\")))", e)
  expect_equal(.case("pycall(v1, \"to_matrix\", B)", e), "Matrix([[x], [y*cos(q) + z*sin(q)], [-y*sin(q) + z*cos(q)]])")
  .run("B <- ReferenceFrame(\"B\")", e)
  .run("v5 <- (((x * pyattr(A, \"x\")) + (y * pyattr(A, \"y\"))) + (z * pyattr(B, \"z\")))", e)
  expect_equal(.case("pycall(PhysVector(0), \"separate\")", e), "{}")
  expect_equal(.case("pycall(v1, \"separate\")", e), "{A: x*A.x + y*A.y + z*A.z}")
  expect_equal(.case("pycall(v5, \"separate\")", e), "{A: x*A.x + y*A.y, B: z*B.z}")
  .run("v6 <- (((x * pyattr(A, \"x\")) + (y * pyattr(A, \"y\"))) + (z * pyattr(A, \"z\")))", e)
  expect_equal(.case("pycall(v6, \"free_symbols\", A)", e), "{x, y, z}")
  expect_equal(.case("pycall(v3, \"applyfunc\", v1)", e), "ERR TypeError: `f` must be callable.")
})

test_that("vector/tests/test_vector.py::test_Vector_diffs matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run(".tmp <- pylist(dynamicsymbols(\"q1 q2 q3 q4\")); q1 <- .tmp[[1]]; q2 <- .tmp[[2]]; q3 <- .tmp[[3]]; q4 <- .tmp[[4]]", e)
  .run(".tmp <- pylist(dynamicsymbols(\"q1 q2 q3 q4\", 1)); q1d <- .tmp[[1]]; q2d <- .tmp[[2]]; q3d <- .tmp[[3]]; q4d <- .tmp[[4]]", e)
  .run(".tmp <- pylist(dynamicsymbols(\"q1 q2 q3 q4\", 2)); q1dd <- .tmp[[1]]; q2dd <- .tmp[[2]]; q3dd <- .tmp[[3]]; q4dd <- .tmp[[4]]", e)
  .run("N <- ReferenceFrame(\"N\")", e)
  .run("A <- pycall(N, \"orientnew\", \"A\", \"Axis\", list(q3, pyattr(N, \"z\")))", e)
  .run("B <- pycall(A, \"orientnew\", \"B\", \"Axis\", list(q2, pyattr(A, \"x\")))", e)
  .run("v1 <- ((q2 * pyattr(A, \"x\")) + (q3 * pyattr(N, \"y\")))", e)
  .run("v2 <- ((q3 * pyattr(B, \"x\")) + v1)", e)
  .run("v3 <- pycall(v1, \"dt\", B)", e)
  .run("v4 <- pycall(v2, \"dt\", B)", e)
  .run("v5 <- (((q1 * pyattr(A, \"x\")) + (q2 * pyattr(A, \"y\"))) + (q3 * pyattr(A, \"z\")))", e)
  expect_equal(.case("pycall(v1, \"dt\", N)", e), "Derivative(q3(t), t)*N.y + Derivative(q2(t), t)*A.x + q2(t)*Derivative(q3(t), t)*A.y")
  expect_equal(.case("pycall(v1, \"dt\", A)", e), "q3(t)*Derivative(q3(t), t)*N.x + Derivative(q3(t), t)*N.y + Derivative(q2(t), t)*A.x")
  expect_equal(.case("pycall(v1, \"dt\", B)", e), "q3(t)*Derivative(q3(t), t)*N.x + Derivative(q3(t), t)*N.y - q3(t)*cos(q3(t))*Derivative(q2(t), t)*N.z + Derivative(q2(t), t)*A.x")
  expect_equal(.case("pycall(v2, \"dt\", N)", e), "Derivative(q3(t), t)*N.y + Derivative(q2(t), t)*A.x + q2(t)*Derivative(q3(t), t)*A.y + Derivative(q3(t), t)*B.x + q3(t)*cos(q2(t))*Derivative(q3(t), t)*B.y - q3(t)*sin(q2(t))*Derivative(q3(t), t)*B.z")
  expect_equal(.case("pycall(v2, \"dt\", A)", e), "q3(t)*Derivative(q3(t), t)*N.x + Derivative(q3(t), t)*N.y + Derivative(q2(t), t)*A.x + Derivative(q3(t), t)*B.x")
  expect_equal(.case("pycall(v2, \"dt\", B)", e), "q3(t)*Derivative(q3(t), t)*N.x + Derivative(q3(t), t)*N.y - q3(t)*cos(q3(t))*Derivative(q2(t), t)*N.z + Derivative(q2(t), t)*A.x + Derivative(q3(t), t)*B.x")
  expect_equal(.case("pycall(v3, \"dt\", N)", e), "(q3(t)*Derivative(q3(t), (t, 2)) + Derivative(q3(t), t)^2)*N.x + Derivative(q3(t), (t, 2))*N.y + (q3(t)*sin(q3(t))*Derivative(q2(t), t)*Derivative(q3(t), t) - q3(t)*cos(q3(t))*Derivative(q2(t), (t, 2)) - cos(q3(t))*Derivative(q2(t), t)*Derivative(q3(t), t))*N.z + Derivative(q2(t), (t, 2))*A.x + Derivative(q2(t), t)*Derivative(q3(t), t)*A.y")
  expect_equal(.case("pycall(v3, \"dt\", A)", e), "(q3(t)*Derivative(q3(t), (t, 2)) + 2*Derivative(q3(t), t)^2)*N.x + (-q3(t)*Derivative(q3(t), t)^2 + Derivative(q3(t), (t, 2)))*N.y + (q3(t)*sin(q3(t))*Derivative(q2(t), t)*Derivative(q3(t), t) - q3(t)*cos(q3(t))*Derivative(q2(t), (t, 2)) - cos(q3(t))*Derivative(q2(t), t)*Derivative(q3(t), t))*N.z + Derivative(q2(t), (t, 2))*A.x")
  expect_equal(.case("pycall(pycall((pycall(v3, \"dt\", B) - (((((q2dd * pyattr(A, \"x\")) - (((q3 * cos(q3)) * pypow(q2d, 2)) * pyattr(A, \"y\"))) + (((2 * pypow(q3d, 2)) + (q3 * q3dd)) * pyattr(N, \"x\"))) + ((q3dd - (q3 * pypow(q3d, 2))) * pyattr(N, \"y\"))) + (((((((2 * q3) * sin(q3)) * q2d) * q3d) - (((2 * cos(q3)) * q2d) * q3d)) - ((q3 * cos(q3)) * q2dd)) * pyattr(N, \"z\")))), \"express\", B), \"simplify\")", e), "0")
  expect_equal(.case("pycall(v4, \"dt\", N)", e), "(q3(t)*Derivative(q3(t), (t, 2)) + Derivative(q3(t), t)^2)*N.x + Derivative(q3(t), (t, 2))*N.y + (q3(t)*sin(q3(t))*Derivative(q2(t), t)*Derivative(q3(t), t) - q3(t)*cos(q3(t))*Derivative(q2(t), (t, 2)) - cos(q3(t))*Derivative(q2(t), t)*Derivative(q3(t), t))*N.z + Derivative(q2(t), (t, 2))*A.x + Derivative(q2(t), t)*Derivative(q3(t), t)*A.y + Derivative(q3(t), (t, 2))*B.x + cos(q2(t))*Derivative(q3(t), t)^2*B.y - sin(q2(t))*Derivative(q3(t), t)^2*B.z")
  expect_equal(.case("pycall(v4, \"dt\", A)", e), "(q3(t)*Derivative(q3(t), (t, 2)) + 2*Derivative(q3(t), t)^2)*N.x + (-q3(t)*Derivative(q3(t), t)^2 + Derivative(q3(t), (t, 2)))*N.y + (q3(t)*sin(q3(t))*Derivative(q2(t), t)*Derivative(q3(t), t) - q3(t)*cos(q3(t))*Derivative(q2(t), (t, 2)) - cos(q3(t))*Derivative(q2(t), t)*Derivative(q3(t), t))*N.z + Derivative(q2(t), (t, 2))*A.x + Derivative(q3(t), (t, 2))*B.x")
  expect_equal(.case("pycall(pycall((pycall(v4, \"dt\", B) - ((((((q2dd * pyattr(A, \"x\")) - (((q3 * cos(q3)) * pypow(q2d, 2)) * pyattr(A, \"y\"))) + (q3dd * pyattr(B, \"x\"))) + (((2 * pypow(q3d, 2)) + (q3 * q3dd)) * pyattr(N, \"x\"))) + ((q3dd - (q3 * pypow(q3d, 2))) * pyattr(N, \"y\"))) + (((((((2 * q3) * sin(q3)) * q2d) * q3d) - (((2 * cos(q3)) * q2d) * q3d)) - ((q3 * cos(q3)) * q2dd)) * pyattr(N, \"z\")))), \"express\", B), \"simplify\")", e), "0")
  expect_equal(.case("pycall(v5, \"dt\", B)", e), "Derivative(q1(t), t)*A.x + (q3(t)*Derivative(q2(t), t) + Derivative(q2(t), t))*A.y + (-q2(t)*Derivative(q2(t), t) + Derivative(q3(t), t))*A.z")
  expect_equal(.case("pycall(v5, \"dt\", A)", e), "Derivative(q1(t), t)*A.x + Derivative(q2(t), t)*A.y + Derivative(q3(t), t)*A.z")
  expect_equal(.case("pycall(v5, \"dt\", N)", e), "(-q2(t)*Derivative(q3(t), t) + Derivative(q1(t), t))*A.x + (q1(t)*Derivative(q3(t), t) + Derivative(q2(t), t))*A.y + Derivative(q3(t), t)*A.z")
  expect_equal(.case("pycall(v3, \"diff\", q1d, N)", e), "0")
  expect_equal(.case("pycall(v3, \"diff\", q2d, N)", e), "- q3(t)*cos(q3(t))*N.z + A.x")
  expect_equal(.case("pycall(v3, \"diff\", q3d, N)", e), "q3(t)*N.x + N.y")
  expect_equal(.case("pycall(v3, \"diff\", q1d, A)", e), "0")
  expect_equal(.case("pycall(v3, \"diff\", q2d, A)", e), "- q3(t)*cos(q3(t))*N.z + A.x")
  expect_equal(.case("pycall(v3, \"diff\", q3d, A)", e), "q3(t)*N.x + N.y")
  expect_equal(.case("pycall(v3, \"diff\", q1d, B)", e), "0")
  expect_equal(.case("pycall(v3, \"diff\", q2d, B)", e), "- q3(t)*cos(q3(t))*N.z + A.x")
  expect_equal(.case("pycall(v3, \"diff\", q3d, B)", e), "q3(t)*N.x + N.y")
  expect_equal(.case("pycall(v4, \"diff\", q1d, N)", e), "0")
  expect_equal(.case("pycall(v4, \"diff\", q2d, N)", e), "- q3(t)*cos(q3(t))*N.z + A.x")
  expect_equal(.case("pycall(v4, \"diff\", q3d, N)", e), "q3(t)*N.x + N.y + B.x")
  expect_equal(.case("pycall(v4, \"diff\", q1d, A)", e), "0")
  expect_equal(.case("pycall(v4, \"diff\", q2d, A)", e), "- q3(t)*cos(q3(t))*N.z + A.x")
  expect_equal(.case("pycall(v4, \"diff\", q3d, A)", e), "q3(t)*N.x + N.y + B.x")
  expect_equal(.case("pycall(v4, \"diff\", q1d, B)", e), "0")
  expect_equal(.case("pycall(v4, \"diff\", q2d, B)", e), "- q3(t)*cos(q3(t))*N.z + A.x")
  expect_equal(.case("pycall(v4, \"diff\", q3d, B)", e), "q3(t)*N.x + N.y + B.x")
  .run("v6 <- (((pypow(q2, 2) * pyattr(N, \"y\")) + (pypow(q2, 2) * pyattr(A, \"y\"))) + (pypow(q2, 2) * pyattr(B, \"y\")))", e)
  .run("n_measy <- (2 * q2)", e)
  .run("a_measy <- (2 * q2)", e)
  .run("b_measx <- pycall(pycall((pypow(q2, 2) * pyattr(B, \"y\")), \"dot\", pyattr(N, \"x\")), \"diff\", q2)", e)
  .run("b_measy <- pycall(pycall((pypow(q2, 2) * pyattr(B, \"y\")), \"dot\", pyattr(N, \"y\")), \"diff\", q2)", e)
  .run("b_measz <- pycall(pycall((pypow(q2, 2) * pyattr(B, \"y\")), \"dot\", pyattr(N, \"z\")), \"diff\", q2)", e)
  .run(".tmp <- pylist(pyattr(pycall(v6, \"diff\", q2, N), \"args\")); n_comp <- .tmp[[1]]; a_comp <- .tmp[[2]]", e)
  expect_equal(.case("length(pyattr(pycall(v6, \"diff\", q2, N), \"args\"))", e), "2")
  expect_equal(.case("pyget(n_comp, 1)", e), "N")
  expect_equal(.case("pyget(a_comp, 1)", e), "A")
  expect_equal(.case("pyget(n_comp, 0)", e), "Matrix([[q2(t)^2*sin(q2(t))*sin(q3(t)) - 2*q2(t)*sin(q3(t))*cos(q2(t))], [-q2(t)^2*sin(q2(t))*cos(q3(t)) + 2*q2(t)*cos(q2(t))*cos(q3(t)) + 2*q2(t)], [q2(t)^2*cos(q2(t)) + 2*q2(t)*sin(q2(t))]])")
  expect_equal(.case("pyget(a_comp, 0)", e), "Matrix([[0], [2*q2(t)], [0]])")
})

test_that("vector/tests/test_vector.py::test_vector_var_in_dcm matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run("N <- ReferenceFrame(\"N\")", e)
  .run("A <- ReferenceFrame(\"A\")", e)
  .run("B <- ReferenceFrame(\"B\")", e)
  .run(".tmp <- pylist(dynamicsymbols(\"u1 u2 u3 u4\")); u1 <- .tmp[[1]]; u2 <- .tmp[[2]]; u3 <- .tmp[[3]]; u4 <- .tmp[[4]]", e)
  .run("v <- ((((u1 * u2) * pyattr(A, \"x\")) + (u3 * pyattr(N, \"y\"))) + (pypow(u4, 2) * pyattr(N, \"z\")))", e)
  expect_equal(.case("pycall(v, \"diff\", u1, N, var_in_dcm = FALSE)", e), "u2(t)*A.x")
  expect_equal(.case("pycall(v, \"diff\", u1, A, var_in_dcm = FALSE)", e), "u2(t)*A.x")
  expect_equal(.case("pycall(v, \"diff\", u3, N, var_in_dcm = FALSE)", e), "N.y")
  expect_equal(.case("pycall(v, \"diff\", u3, A, var_in_dcm = FALSE)", e), "N.y")
  expect_equal(.case("pycall(v, \"diff\", u3, B, var_in_dcm = FALSE)", e), "N.y")
  expect_equal(.case("pycall(v, \"diff\", u4, N, var_in_dcm = FALSE)", e), "2*u4(t)*N.z")
  expect_equal(.case("pycall(v, \"diff\", u1, N)", e), "ERR ValueError: No connecting orientation path found between N and A.")
})

test_that("vector/tests/test_vector.py::test_vector_simplify matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run(".tmp <- pylist(pysymbols(\"x, y, z, k, n, m, w, f, s, A\")); x <- .tmp[[1]]; y <- .tmp[[2]]; z <- .tmp[[3]]; k <- .tmp[[4]]; n <- .tmp[[5]]; m <- .tmp[[6]]; w <- .tmp[[7]]; f <- .tmp[[8]]; s <- .tmp[[9]]; A <- .tmp[[10]]", e)
  .run("N <- ReferenceFrame(\"N\")", e)
  .run("test1 <- (((1 / x) + (1 / y)) * pyattr(N, \"x\"))", e)
  expect_equal(.case("(test1 & pyattr(N, \"x\"))", e), "1/y + 1/x")
  .run("test1 <- pycall(test1, \"simplify\")", e)
  expect_equal(.case("(test1 & pyattr(N, \"x\"))", e), "(x + y)/(x*y)")
  .run("test2 <- (((pypow(A, 2) * pypow(s, 4)) / (((4 * S_Pi) * k) * pypow(m, 3))) * pyattr(N, \"x\"))", e)
  .run("test2 <- pycall(test2, \"simplify\")", e)
  expect_equal(.case("(test2 & pyattr(N, \"x\"))", e), "A^2*s^4/(4*pi*k*m^3)")
  .run("test3 <- ((((4 + (4 * x)) - (2 * (2 + (2 * x)))) / (2 + (2 * x))) * pyattr(N, \"x\"))", e)
  .run("test3 <- pycall(test3, \"simplify\")", e)
  expect_equal(.case("(test3 & pyattr(N, \"x\"))", e), "0")
  .run("test4 <- (((((((-4) * x) * pypow(y, 2)) - (2 * pypow(y, 3))) - ((2 * pypow(x, 2)) * y)) / pypow((x + y), 2)) * pyattr(N, \"x\"))", e)
  .run("test4 <- pycall(test4, \"simplify\")", e)
  expect_equal(.case("(test4 & pyattr(N, \"x\"))", e), "-2*y")
})

test_that("vector/tests/test_vector.py::test_vector_evalf matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run(".tmp <- pylist(pysymbols(\"a b\")); a <- .tmp[[1]]; b <- .tmp[[2]]", e)
  .run("v <- (S_Pi * pyattr(A, \"x\"))", e)
  expect_equal(.case("pycall(v, \"evalf\", 2)", e), "3.1*A.x")
  .run("v <- (((S_Pi * pyattr(A, \"x\")) + ((5 * a) * pyattr(A, \"y\"))) - (b * pyattr(A, \"z\")))", e)
  expect_equal(.case("pycall(v, \"evalf\", 5, subs = pydict(list(a, b), list(Float(1.234), Float(5.8973))))", e), "3.1416*A.x + 6.17*A.y - 5.8973*A.z")
})

test_that("vector/tests/test_vector.py::test_vector_angle matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run("A <- ReferenceFrame(\"A\")", e)
  .run("v1 <- (pyattr(A, \"x\") + pyattr(A, \"y\"))", e)
  .run("v2 <- pyattr(A, \"z\")", e)
  expect_equal(.case("pycall(v1, \"angle_between\", v2)", e), "pi/2")
  .run("B <- ReferenceFrame(\"B\")", e)
  .run("pycall(B, \"orient_axis\", A, pyattr(A, \"x\"), S_Pi)", e)
  .run("v3 <- pyattr(A, \"x\")", e)
  .run("v4 <- pyattr(B, \"x\")", e)
  expect_equal(.case("pycall(v3, \"angle_between\", v4)", e), "0")
})

test_that("vector/tests/test_vector.py::test_vector_xreplace matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run(".tmp <- pylist(pysymbols(\"x y z\")); x <- .tmp[[1]]; y <- .tmp[[2]]; z <- .tmp[[3]]", e)
  .run("v <- (((pypow(x, 2) * pyattr(A, \"x\")) + ((x * y) * pyattr(A, \"y\"))) + (((x * y) * z) * pyattr(A, \"z\")))", e)
  expect_equal(.case("pycall(v, \"xreplace\", pydict(list(x), list(cos(x))))", e), "cos(x)^2*A.x + y*cos(x)*A.y + y*z*cos(x)*A.z")
  expect_equal(.case("pycall(v, \"xreplace\", pydict(list((x * y)), list(S_Pi)))", e), "x^2*A.x + pi*A.y + x*y*z*A.z")
  expect_equal(.case("pycall(v, \"xreplace\", pydict(list(((x * y) * z)), list(1)))", e), "x^2*A.x + x*y*A.y + A.z")
  expect_equal(.case("pycall(v, \"xreplace\", pydict(list(x, z), list(1, 0)))", e), "A.x + y*A.y")
  expect_equal(.case("pycall(v, \"xreplace\")", e), "ERR TypeError: Vector.xreplace() missing 1 required positional argument: 'rule'")
  expect_equal(.case("pycall(v, \"xreplace\", list(x, y))", e), "ERR TypeError: list indices must be integers or slices, not Symbol")
})

test_that("vector/tests/test_vector.py::test_issue_23366 matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run("u1 <- dynamicsymbols(\"u1\")", e)
  .run("N <- ReferenceFrame(\"N\")", e)
  .run("N_v_A <- (u1 * pyattr(N, \"x\"))", e)
  expect_equal(.case("pycall(N_v_A, \"diff\", N, u1)", e), "ERR VectorTypeError: \\nExpected an instance of <class\\n'sympy.physics.vector.frame.ReferenceFrame'>, but received object\\n'u1(t)' of u1.")
})

test_that("vector/tests/test_vector.py::test_vector_outer matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run(".tmp <- pylist(pysymbols(\"a, b, c, d, e, f\")); a <- .tmp[[1]]; b <- .tmp[[2]]; c <- .tmp[[3]]; d <- .tmp[[4]]; e <- .tmp[[5]]; f <- .tmp[[6]]", e)
  .run("N <- ReferenceFrame(\"N\")", e)
  .run("v1 <- (((a * pyattr(N, \"x\")) + (b * pyattr(N, \"y\"))) + (c * pyattr(N, \"z\")))", e)
  .run("v2 <- (((d * pyattr(N, \"x\")) + (e * pyattr(N, \"y\"))) + (f * pyattr(N, \"z\")))", e)
  .run("v1v2 <- Matrix(list(list((a * d), (a * e), (a * f)), list((b * d), (b * e), (b * f)), list((c * d), (c * e), (c * f))))", e)
  expect_equal(.case("pycall(pycall(v1, \"outer\", v2), \"to_matrix\", N)", e), "Matrix([[a*d, a*e, a*f], [b*d, b*e, b*f], [c*d, c*e, c*f]])")
  expect_equal(.case("pycall((v1 | v2), \"to_matrix\", N)", e), "Matrix([[a*d, a*e, a*f], [b*d, b*e, b*f], [c*d, c*e, c*f]])")
  .run("v2v1 <- Matrix(list(list((d * a), (d * b), (d * c)), list((e * a), (e * b), (e * c)), list((f * a), (f * b), (f * c))))", e)
  expect_equal(.case("pycall(pycall(v2, \"outer\", v1), \"to_matrix\", N)", e), "Matrix([[a*d, b*d, c*d], [a*e, b*e, c*e], [a*f, b*f, c*f]])")
  expect_equal(.case("pycall((v2 | v1), \"to_matrix\", N)", e), "Matrix([[a*d, b*d, c*d], [a*e, b*e, c*e], [a*f, b*f, c*f]])")
})

test_that("vector/tests/test_vector.py::test_overloaded_operators matches SymPy", {
  e <- new.env(parent = .menv)
  .wstate$mode <- "default"
  .run(".tmp <- pylist(pysymbols(\"a, b, c, d, e, f\")); a <- .tmp[[1]]; b <- .tmp[[2]]; c <- .tmp[[3]]; d <- .tmp[[4]]; e <- .tmp[[5]]; f <- .tmp[[6]]", e)
  .run("N <- ReferenceFrame(\"N\")", e)
  .run("v1 <- (((a * pyattr(N, \"x\")) + (b * pyattr(N, \"y\"))) + (c * pyattr(N, \"z\")))", e)
  .run("v2 <- (((d * pyattr(N, \"x\")) + (e * pyattr(N, \"y\"))) + (f * pyattr(N, \"z\")))", e)
  expect_equal(.case("(v1 + v2)", e), "(a + d)*N.x + (b + e)*N.y + (c + f)*N.z")
  expect_equal(.case("(v1 - v2)", e), "(a - d)*N.x + (b - e)*N.y + (c - f)*N.z")
  expect_equal(.case("(v1 & v2)", e), "a*d + b*e + c*f")
  expect_equal(.case("pyxor(v1, v2)", e), "(b*f - c*e)*N.x + (-a*f + c*d)*N.y + (a*e - b*d)*N.z")
  expect_equal(.case("pyxor(v2, v1)", e), "(-b*f + c*e)*N.x + (a*f - c*d)*N.y + (-a*e + b*d)*N.z")
})

