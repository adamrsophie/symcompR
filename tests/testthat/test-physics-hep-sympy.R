# Expected values are the outputs of SymPy for the same calls (sympy.physics.hep.gamma_matrices).
# Each test replays the statements of one SymPy test function and checks the
# assertions whose results agree with SymPy (44 rows).

skip_unless_full_tests()

.bt <- new.env(parent = asNamespace("symcompR"))
eval(parse(text = r"---(# Helpers that give the translated SymPy test code Python semantics.
.qbool <- function(x) {
  if (is.null(x)) return(FALSE)
  if (is.logical(x)) return(isTRUE(x))
  if (is_sym(x)) {
    if (x$type == "BooleanAtom") return(isTRUE(x$value))
    return(!(is_number_type(x) && num_is_zero(x)))
  }
  if (is.numeric(x)) return(length(x) > 0L && x[[1L]] != 0)
  length(x) > 0L
}
.qlen <- function(x) {
  if (is_sym(x) && !is.null(.qnode_def(x$type)) && .qhas(x, "__len__")) return(.qmeth(x, "__len__"))
  if (is_sym(x)) return(length(x$args))
  if (inherits(x, "symmatrix")) return(x$nrow * x$ncol)
  if (is.character(x) && length(x) == 1L) return(nchar(x))
  length(x)
}
.qseq <- function(x) {
  if (is_sym(x)) return(if (x$type %in% c("Add", "Mul")) .sympy_args(x) else x$args)
  if (inherits(x, "symmatrix")) return(x$data)
  if (is.character(x) && length(x) == 1L) return(as.list(strsplit(x, "")[[1L]]))
  as.list(x)
}
.qidx <- function(x, i) {
  if (is.function(x) && !is.null(attr(x, "base"))) return(x(i))
  if (!is.null(attr(x, "qdict"))) {
    for (p in x) if (.qeq(p[[1L]], i)) return(p[[2L]])
    stop(.qerror("KeyError", .qrepr(i)))
  }
  if ((is.list(x) || is.environment(x)) && !is_sym(x) && is.function(x[["__getitem__"]])) return(x[["__getitem__"]](i))
  i <- .qint(i)
  if (is_sym(x) && !is.null(.qnode_def(x$type)) && .qhas(x, "__getitem__")) return(.qmeth(x, "__getitem__", i))
  s <- .qseq(x)
  n <- length(s)
  j <- if (i < 0) n + i + 1 else i + 1
  if (j < 1 || j > n) stop(.qerror("IndexError", "index out of range"))
  s[[j]]
}
.qslice <- function(x, lo, hi) {
  s <- .qseq(x); n <- length(s)
  lo <- if (is.null(lo)) 0 else .qint(lo); hi <- if (is.null(hi)) n else .qint(hi)
  if (lo < 0) lo <- n + lo; if (hi < 0) hi <- n + hi
  hi <- min(hi, n)
  r <- if (hi > lo) s[(lo + 1):hi] else list()
  if (is.character(x)) return(paste(unlist(r), collapse = ""))
  if (is_sym(x) && x$type == "Tuple") return(do.call(Tuple, r))
  r
}
.qmidx <- function(m, i, j) .mget(m, .qint(i) + 1L, .qint(j) + 1L)
.qrange <- function(a, b = NULL, by = 1) {
  a <- .qint(a); if (is.null(b)) { b <- a; a <- 0 } else b <- .qint(b)
  by <- .qint(by)
  if ((by > 0 && a >= b) || (by < 0 && a <= b)) return(list())
  lapply(seq(a, b - sign(by), by = by), function(k) S(k))
}
.qlist <- function(x = list()) .qseq(x)
.qset <- function(x = list()) {
  s <- .qseq(x); out <- list()
  for (z in s) if (!any(vapply(out, function(o) .qeq(o, z), NA))) out[[length(out) + 1L]] <- z
  structure(out, qset = TRUE)
}
.qdict <- function(keys, vals) structure(Map(function(k, v) list(k, v), keys, vals), qdict = TRUE)
.qpairs <- function(d) unclass(d)
.qin <- function(a, b) {
  if (!is.null(attr(b, "qdict"))) return(any(vapply(b, function(p) .qeq(p[[1L]], a), NA)))
  if (is.character(b) && length(b) == 1L && is.character(a)) return(grepl(a, b, fixed = TRUE))
  if (is_sym(b) && !is.null(.qnode_def(b$type)) && .qhas(b, "__contains__")) return(.qmeth(b, "__contains__", a))
  any(vapply(.qseq(b), function(z) .qeq(z, a), NA))
}
.qis <- function(a, b) {
  if (is.null(a) || is.null(b)) return(is.null(a) && is.null(b))
  if (is.logical(a) || is.logical(b)) return(identical(a, b) || (is_sym(a) && identical(a$value, b)) || (is_sym(b) && identical(b$value, a)))
  .qeq(a, b)
}
.qsorted <- function(x, key = NULL) {
  s <- .qseq(x)
  if (!length(s)) return(s)
  if (all(vapply(s, function(z) is_sym(z) && is_number_type(z), NA)))
    return(s[order(vapply(s, function(z) as.numeric(N(z)), 0))])
  s[order(vapply(s, function(z) if (is.character(z)) z else sstr(z), ""), method = "radix")]
}
.qcallobj <- function(f, ...) {
  if (is.function(f)) return(f(...))
  if (is.character(f) && .qis_class(f)) {
    g <- get0(f, mode = "function")
    if (!is.null(g)) return(g(...))
    return(.qfunc(f, list(...)))
  }
  if (is_sym(f) && !is.null(.qnode_def(f$type))) return(.qmeth(f, "__call__", ...))
  if (is.character(f) && !is.null(g <- get0(f, mode = "function"))) return(g(...))
  stop(.qerror("TypeError", "object is not callable"))
}
.qis_tensor_eq <- function(a, b) {
  a <- canon_bp(a); b <- canon_bp(b)
  if (inherits(a, "TensExpr")) return(.hep_equals(a, b))
  if (inherits(b, "TensExpr")) return(.hep_equals(b, a))
  .qeq(a, b)
}
qcall <- function(obj, name, ...) {
  if (name == "equals" && inherits(obj, c("TensExpr", "symmatrix"))) return(.hep_equals(obj, ...))
  # Expr.equals: structurally equal, or the difference simplifies to zero
  if (name == "equals" && is_sym(obj) && is.null(.qnode_def(obj$type))) {
    o <- as_sym(..1)
    return(.qeq(obj, o) || .qeq(simplify(expand(.sub(obj, o))), S_Zero))
  }
  if ((is.list(obj) || is.environment(obj)) && !is_sym(obj) && is.function(obj[[name]])) return(obj[[name]](...))
  if (is_sym(obj) && !is.null(.qnode_def(obj$type))) return(.qmeth(obj, name, ...))
  if (is.character(obj) && .qis_class(obj)) return(.qmeth(obj, name, ...))
  f <- get0(name, mode = "function")
  if (is.null(f)) stop(.qerror("AttributeError", "'", .qclassname(obj), "' object has no attribute '", name, "'"))
  f(obj, ...)
}
qattr <- function(obj, name) {
  if (.is_kd(obj) && !is.null(r <- .kd_attr(obj, name))) return(r)
  if (inherits(obj, c("TensorIndexType", "TensorHead", "TensExpr", "TensorIndex"))) return(do.call("$", list(obj, name)))
  if (startsWith(name, "is_") && is_sym(obj) && (is.null(.qnode_def(obj$type)) || !.qhas(obj, name))) {
    v <- ask_fact(obj, substring(name, 4L))
    return(if (is.na(v)) NULL else v)
  }
  symcompR:::qattr(obj, name)
}
Mul <- function(...) .mul(lapply(list(...), as_sym))
Add <- function(...) .add(lapply(list(...), as_sym))
Pow <- function(a, b, evaluate = TRUE)
  if (isFALSE(evaluate)) new_sym("Pow", list(as_sym(a), as_sym(b))) else .pow(as_sym(a), as_sym(b))
.qfmt <- function(x) {
  if (is.null(x)) return("None")
  if (is.character(x) && length(x) == 1L) return(if (.qis_class(x)) .qpyname(x) else x)
  if (is.logical(x) && length(x) == 1L) return(if (is.na(x)) "None" else if (x) "True" else "False")
  if (is_sym(x) || inherits(x, c("symmatrix", "TensExpr"))) return(gsub("**", "^", sstr(x), fixed = TRUE))
  if (is.function(x)) return("OBJ function")
  if (is.numeric(x) && length(x) == 1L) {
    if (is.finite(x) && x == round(x) && abs(x) < 1e15) return(sprintf("%.0f", x))
    for (d in 15:17) { s <- formatC(x, digits = d, format = "g"); if (as.numeric(s) == x) break }
    return(s)
  }
  if (!is.null(attr(x, "qdict"))) return(paste0("{", paste(vapply(x, function(p) paste0(.qfmt(p[[1L]]), ": ", .qfmt(p[[2L]])), ""), collapse = ", "), "}"))
  if (!is.null(attr(x, "qset"))) return(paste0("{", paste(sort(vapply(x, .qfmt, ""), method = "radix"), collapse = ", "), "}"))
  if (is.list(x) || length(x) > 1L) return(paste0("[", paste(vapply(as.list(x), .qfmt, ""), collapse = ", "), "]"))
  paste("OBJ", class(x)[1L])
}
.qerrfmt <- function(e) {
  cl <- setdiff(class(e), c("simpleError", "error", "condition", "sympy_error", "rlang_error"))
  paste0("ERR ", if (length(cl)) cl[1L] else "Error", ": ", conditionMessage(e))
}
.pymul <- function(a, b) {
  if (inherits(a, "TensExpr") || inherits(b, "TensExpr")) return(.tn_mul(a, b))
  if (inherits(a, "symmatrix") && inherits(b, "symmatrix")) return(.matmul(a, b))
  if (is.numeric(a) && is.numeric(b)) return(a * b)
  if (!is_sym(a) && !inherits(a, "symmatrix")) a <- as_sym(a)
  a * b
}
.qtype <- function(x) { r <- .qclassname(x); if (is_sym(x) && is.null(.qnode_def(x$type)) && x$type == "Function") x$name else r }

.qpow <- function(a, b) as_sym(a)^b
# sympy.core.expr.unchanged: func(*args) is not rebuilt into something else
unchanged <- function(func, ...) {
  f <- .qcallobj(func, ...)
  .qeq(.qclassname(f), .qpyname(func)) && .qeq(qattr(f, "args"), list(...))
}
# Basic.free_symbols is a set
free_symbols <- function(e) {
  r <- get("free_symbols", envir = asNamespace("symcompR"))(e)
  if (is.list(r) && !is_sym(r)) structure(r, qset = TRUE) else r
}
# sympy.numbered_symbols is an endless generator
numbered_symbols <- function(prefix = "x", start = 0, ...) .qsymbol_gen(prefix, .qint(start))
sym_shape <- function(m) if (inherits(m, "symmatrix")) list(m$nrow, m$ncol) else shape(m)

# Python repr(): a class's own __repr__, else the str form
repr <- function(x) {
  if (is_sym(x) && !is.null(.qnode_def(x$type))) {
    f <- .qfind(x$type, "__repr__")$val
    if (!is.null(f)) return(f(x))
    return(.str(x))
  }
  if (inherits(x, c("VarBosonicBasis", "FixedBosonicBasis"))) return(format(x))
  if (is_sym(x)) return(.str(x))
  .qfmt(x)
}

# The classes that the SymPy test files define, registered for the file whose
# rows run next (several files reuse the names Foo and CustomKet).
.qtestclasses <- function(file, env = parent.frame()) {
  mk <- function(name, bases, ..., .l = list()) {
    do.call(.qclass, c(list(name, bases), list(...), .l))
    assign(name, function(...) .qfunc(name, list(...)), envir = env)
  }
  dflt <- function(...) { a <- list(...); .cm(function(cls) a) }
  switch(file,
    test_commutator = {
      mk("Foo", "Operator", `_eval_commutator_Bar` = function(self, bar, ...) S_Zero)
      mk("Bar", "Operator")
      mk("Tam", "Operator", `_eval_commutator_Foo` = function(self, foo, ...) S_One)
    },
    test_anticommutator = {
      mk("Foo", "Operator", `_eval_anticommutator_Bar` = function(self, bar, ...) S_Zero)
      mk("Bar", "Operator")
      mk("Tam", "Operator", `_eval_anticommutator_Foo` = function(self, foo, ...) S_One)
    },
    test_qapply = mk("Foo", "Operator", `_apply_operator_JzKet` = function(self, ket, ...) ket),
    test_dagger = mk("Foo", "Expr", .l = c(.qbasic_printers, list(is_commutative = TRUE,
                     `_eval_adjoint` = function(self, ...) sym_I))),
    test_density = mk("TestTimeDepKet", "TimeDepKet", `_eval_trace` = function(self, bra, ...) S_One),
    test_trace = mk("Foo", "Basic", trace = function(self) S_One),
    test_innerproduct = {
      mk("FooState", "StateBase")
      mk("FooKet", c("Ket", "FooState"), dual_class = .cm(function(cls) "FooBra"),
         `_eval_innerproduct_FooBra` = function(self, bra, ...) S_One,
         `_eval_innerproduct_BarBra` = function(self, bra, ...) sym_I)
      mk("FooBra", c("Bra", "FooState"), dual_class = .cm(function(cls) "FooKet"))
      mk("BarState", "StateBase")
      mk("BarKet", c("Ket", "BarState"), dual_class = .cm(function(cls) "BarBra"))
      mk("BarBra", c("Bra", "BarState"), dual_class = .cm(function(cls) "BarKet"))
    },
    test_operator = {
      mk("CustomKet", "Ket", default_args = dflt("t"))
      mk("CustomOp", "HermitianOperator", default_args = dflt("T"))
    },
    test_state = {
      mk("CustomKet", "Ket", default_args = dflt("test"))
      mk("CustomKetMultipleLabels", "Ket", default_args = dflt("r", "theta", "phi"))
      mk("CustomTimeDepKet", "TimeDepKet", default_args = dflt("test", "t"))
      mk("CustomTimeDepKetMultipleLabels", "TimeDepKet", default_args = dflt("r", "theta", "phi", "t"))
    },
    test_circuitutils = , test_identitysearch =
      assign("create_gate_sequence", function(qubit = 0) { q <- as_sym(qubit); list(XGate(q), YGate(q), ZGate(q), HadamardGate(q)) },
             envir = env),
    test_grover = {
      assign("return_one_on_two", function(qubits) .qeq(qubits, IntQubit(2, qattr(qubits, "nqubits"))), envir = env)
      assign("return_one_on_one", function(qubits) .qeq(qubits, IntQubit(1, nqubits = qattr(qubits, "nqubits"))), envir = env)
    },
    test_represent = {
      rep <- function(v) list(`_represent_default_basis` = function(self, ...) .qmeth(self, "_represent_AOp", NULL, ...),
                              `_represent_AOp` = function(self, basis, ...) get(v, envir = env))
      mk("AKet", "Ket", .l = c(list(dual_class = .cm(function(cls) "ABra")), rep("Avec")))
      mk("ABra", "Bra", dual_class = .cm(function(cls) "AKet"))
      mk("AOp", "Operator", .l = rep("Amat"))
      mk("BOp", "Operator", .l = rep("Bmat"))
    },
    NULL)
  invisible(NULL)
}

# att() of sympy/physics/tests/test_secondquant.py
att <- function(...) {
  a <- list(...)
  if (length(a) == 4L) return(AntiSymmetricTensor("t", a[1:2], a[3:4]))
  if (length(a) == 2L) return(AntiSymmetricTensor("t", a[1L], a[2L]))
  NULL
}
)---"), envir = .bt)
.new <- function(mod, file) {
  e <- new.env(parent = .bt)
  .bt$.qtestclasses(file, e)
  for (st in mod) try(suppressWarnings(eval(parse(text = st), envir = e)), silent = TRUE)
  e
}
.do <- function(e, st) invisible(try(suppressWarnings(eval(parse(text = st), envir = e)), silent = TRUE))
.run <- function(e, code) tryCatch(suppressWarnings(.bt$.qfmt(eval(parse(text = code), envir = e))),
                                   error = function(err) .bt$.qerrfmt(err))

test_that("test_gamma_matrices: test_kahane_simplify1", {
  e <- .new(character(), 'test_gamma_matrices')
  .do(e, r"---(.tmp <- tensor_indices("i0:16", LorentzIndex))---")
  .do(e, r"---(i0 <- .tmp[[1]])---")
  .do(e, r"---(i1 <- .tmp[[2]])---")
  .do(e, r"---(i2 <- .tmp[[3]])---")
  .do(e, r"---(i3 <- .tmp[[4]])---")
  .do(e, r"---(i4 <- .tmp[[5]])---")
  .do(e, r"---(i5 <- .tmp[[6]])---")
  .do(e, r"---(i6 <- .tmp[[7]])---")
  .do(e, r"---(i7 <- .tmp[[8]])---")
  .do(e, r"---(i8 <- .tmp[[9]])---")
  .do(e, r"---(i9 <- .tmp[[10]])---")
  .do(e, r"---(i10 <- .tmp[[11]])---")
  .do(e, r"---(i11 <- .tmp[[12]])---")
  .do(e, r"---(i12 <- .tmp[[13]])---")
  .do(e, r"---(i13 <- .tmp[[14]])---")
  .do(e, r"---(i14 <- .tmp[[15]])---")
  .do(e, r"---(i15 <- .tmp[[16]])---")
  .do(e, r"---(.tmp <- tensor_indices("mu, nu, rho, sigma", LorentzIndex))---")
  .do(e, r"---(mu <- .tmp[[1]])---")
  .do(e, r"---(nu <- .tmp[[2]])---")
  .do(e, r"---(rho <- .tmp[[3]])---")
  .do(e, r"---(sigma <- .tmp[[4]])---")
  .do(e, r"---(v_D <- S(4))---")
  .do(e, r"---(v_t <- .pymul(GammaMatrix(i0), GammaMatrix(i1)))---")
  .do(e, r"---(r <- kahane_simplify(v_t))---")
  expect_equal(.run(e, r"---(qcall(r, "equals", v_t))---"), r"---(True)---")
  .do(e, r"---(v_t <- .pymul(.pymul(GammaMatrix(i0), GammaMatrix(i1)), GammaMatrix((-i0))))---")
  .do(e, r"---(r <- kahane_simplify(v_t))---")
  expect_equal(.run(e, r"---(qcall(r, "equals", .pymul(S(-2), GammaMatrix(i1))))---"), r"---(True)---")
  .do(e, r"---(v_t <- .pymul(.pymul(GammaMatrix(i0), GammaMatrix(i1)), GammaMatrix((-i0))))---")
  .do(e, r"---(r <- kahane_simplify(v_t))---")
  expect_equal(.run(e, r"---(qcall(r, "equals", .pymul(S(-2), GammaMatrix(i1))))---"), r"---(True)---")
  .do(e, r"---(v_t <- .pymul(GammaMatrix(i0), GammaMatrix(i1)))---")
  .do(e, r"---(r <- kahane_simplify(v_t))---")
  expect_equal(.run(e, r"---(qcall(r, "equals", v_t))---"), r"---(True)---")
  .do(e, r"---(v_t <- .pymul(GammaMatrix(i0), GammaMatrix(i1)))---")
  .do(e, r"---(r <- kahane_simplify(v_t))---")
  expect_equal(.run(e, r"---(qcall(r, "equals", v_t))---"), r"---(True)---")
  .do(e, r"---(v_t <- .pymul(GammaMatrix(i0), GammaMatrix((-i0))))---")
  .do(e, r"---(r <- kahane_simplify(v_t))---")
  expect_equal(.run(e, r"---(qcall(r, "equals", .pymul(S(4), eye(4))))---"), r"---(True)---")
  .do(e, r"---(v_t <- .pymul(GammaMatrix(i0), GammaMatrix((-i0))))---")
  .do(e, r"---(r <- kahane_simplify(v_t))---")
  expect_equal(.run(e, r"---(qcall(r, "equals", .pymul(S(4), eye(4))))---"), r"---(True)---")
  .do(e, r"---(v_t <- .pymul(GammaMatrix(i0), GammaMatrix((-i0))))---")
  .do(e, r"---(r <- kahane_simplify(v_t))---")
  expect_equal(.run(e, r"---(qcall(r, "equals", .pymul(S(4), eye(4))))---"), r"---(True)---")
  .do(e, r"---(v_t <- .pymul(.pymul(GammaMatrix(i0), GammaMatrix(i1)), GammaMatrix((-i0))))---")
  .do(e, r"---(r <- kahane_simplify(v_t))---")
  expect_equal(.run(e, r"---(qcall(r, "equals", .pymul(S(-2), GammaMatrix(i1))))---"), r"---(True)---")
  .do(e, r"---(v_t <- .pymul(.pymul(.pymul(GammaMatrix(i0), GammaMatrix(i1)), GammaMatrix((-i0))), GammaMatrix((-i1))))---")
  .do(e, r"---(r <- kahane_simplify(v_t))---")
  expect_equal(.run(e, r"---(qcall(r, "equals", .pymul((.pymul(S(2), v_D) - (v_D ^ S(2))), eye(4))))---"), r"---(True)---")
  .do(e, r"---(v_t <- .pymul(.pymul(.pymul(GammaMatrix(i0), GammaMatrix(i1)), GammaMatrix((-i0))), GammaMatrix((-i1))))---")
  .do(e, r"---(r <- kahane_simplify(v_t))---")
  expect_equal(.run(e, r"---(qcall(r, "equals", .pymul((.pymul(S(2), v_D) - (v_D ^ S(2))), eye(4))))---"), r"---(True)---")
  .do(e, r"---(v_t <- .pymul(.pymul(.pymul(GammaMatrix(i0), GammaMatrix((-i0))), GammaMatrix(i1)), GammaMatrix((-i1))))---")
  .do(e, r"---(r <- kahane_simplify(v_t))---")
  expect_equal(.run(e, r"---(qcall(r, "equals", .pymul(S(16), eye(4))))---"), r"---(True)---")
  .do(e, r"---(v_t <- .pymul(.pymul(.pymul(GammaMatrix(mu), GammaMatrix(nu)), GammaMatrix((-nu))), GammaMatrix((-mu))))---")
  .do(e, r"---(r <- kahane_simplify(v_t))---")
  expect_equal(.run(e, r"---(qcall(r, "equals", .pymul((v_D ^ S(2)), eye(4))))---"), r"---(True)---")
  .do(e, r"---(v_t <- .pymul(.pymul(.pymul(GammaMatrix(mu), GammaMatrix(nu)), GammaMatrix((-nu))), GammaMatrix((-mu))))---")
  .do(e, r"---(r <- kahane_simplify(v_t))---")
  expect_equal(.run(e, r"---(qcall(r, "equals", .pymul((v_D ^ S(2)), eye(4))))---"), r"---(True)---")
  .do(e, r"---(v_t <- .pymul(.pymul(.pymul(GammaMatrix(mu), GammaMatrix(nu)), GammaMatrix((-nu))), GammaMatrix((-mu))))---")
  .do(e, r"---(r <- kahane_simplify(v_t))---")
  expect_equal(.run(e, r"---(qcall(r, "equals", .pymul((v_D ^ S(2)), eye(4))))---"), r"---(True)---")
  .do(e, r"---(v_t <- .pymul(.pymul(.pymul(.pymul(.pymul(GammaMatrix(mu), GammaMatrix(nu)), GammaMatrix((-rho))), GammaMatrix((-nu))), GammaMatrix((-mu))), GammaMatrix(rho)))---")
  .do(e, r"---(r <- kahane_simplify(v_t))---")
  expect_equal(.run(e, r"---(qcall(r, "equals", .pymul(((.pymul(S(4), v_D) - .pymul(S(4), (v_D ^ S(2)))) + (v_D ^ S(3))), eye(4))))---"), r"---(True)---")
  .do(e, r"---(v_t <- .pymul(.pymul(.pymul(.pymul(.pymul(.pymul(.pymul(GammaMatrix((-mu)), GammaMatrix((-nu))), GammaMatrix((-rho))), GammaMatrix((-sigma))), GammaMatrix(nu)), GammaMatrix(mu)), GammaMatrix(sigma)), GammaMatrix(rho)))---")
  .do(e, r"---(r <- kahane_simplify(v_t))---")
  expect_equal(.run(e, r"---(qcall(r, "equals", .pymul((((.pymul(S(-16), v_D) + .pymul(S(24), (v_D ^ S(2)))) - .pymul(S(8), (v_D ^ S(3)))) + (v_D ^ S(4))), eye(4))))---"), r"---(True)---")
  .do(e, r"---(v_t <- .pymul(.pymul(.pymul(.pymul(.pymul(.pymul(.pymul(GammaMatrix((-mu)), GammaMatrix(nu)), GammaMatrix((-rho))), GammaMatrix(sigma)), GammaMatrix(rho)), GammaMatrix((-nu))), GammaMatrix(mu)), GammaMatrix((-sigma))))---")
  .do(e, r"---(r <- kahane_simplify(v_t))---")
  expect_equal(.run(e, r"---(qcall(r, "equals", .pymul((((.pymul(S(8), v_D) - .pymul(S(12), (v_D ^ S(2)))) + .pymul(S(6), (v_D ^ S(3)))) - (v_D ^ S(4))), eye(4))))---"), r"---(True)---")
  .do(e, r"---(v_t <- .pymul(.pymul(.pymul(.pymul(GammaMatrix(mu), GammaMatrix(nu)), GammaMatrix(rho)), GammaMatrix(sigma)), GammaMatrix((-mu))))---")
  .do(e, r"---(r <- kahane_simplify(v_t))---")
  expect_equal(.run(e, r"---(qcall(r, "equals", .pymul(.pymul(.pymul(S(-2), GammaMatrix(sigma)), GammaMatrix(rho)), GammaMatrix(nu))))---"), r"---(True)---")
  .do(e, r"---(v_t <- .pymul(.pymul(.pymul(GammaMatrix(mu), GammaMatrix((-mu))), GammaMatrix(rho)), GammaMatrix(sigma)))---")
  .do(e, r"---(r <- kahane_simplify(v_t))---")
  expect_equal(.run(e, r"---(qcall(r, "equals", .pymul(.pymul(S(4), GammaMatrix(rho)), GammaMatrix(sigma))))---"), r"---(True)---")
  .do(e, r"---(v_t <- .pymul(.pymul(.pymul(GammaMatrix(rho), GammaMatrix(sigma)), GammaMatrix(mu)), GammaMatrix((-mu))))---")
  .do(e, r"---(r <- kahane_simplify(v_t))---")
  expect_equal(.run(e, r"---(qcall(r, "equals", .pymul(.pymul(S(4), GammaMatrix(rho)), GammaMatrix(sigma))))---"), r"---(True)---")
})

test_that("test_gamma_matrices: test_gamma_matrix_trace", {
  e <- .new(character(), 'test_gamma_matrices')
  .do(e, r"---(g <- qattr(LorentzIndex, "metric"))---")
  .do(e, r"---(.tmp <- tensor_indices("m0:7", LorentzIndex))---")
  .do(e, r"---(m0 <- .tmp[[1]])---")
  .do(e, r"---(m1 <- .tmp[[2]])---")
  .do(e, r"---(m2 <- .tmp[[3]])---")
  .do(e, r"---(m3 <- .tmp[[4]])---")
  .do(e, r"---(m4 <- .tmp[[5]])---")
  .do(e, r"---(m5 <- .tmp[[6]])---")
  .do(e, r"---(m6 <- .tmp[[7]])---")
  .do(e, r"---(.tmp <- tensor_indices("n0:6", LorentzIndex))---")
  .do(e, r"---(n0 <- .tmp[[1]])---")
  .do(e, r"---(n1 <- .tmp[[2]])---")
  .do(e, r"---(n2 <- .tmp[[3]])---")
  .do(e, r"---(n3 <- .tmp[[4]])---")
  .do(e, r"---(n4 <- .tmp[[5]])---")
  .do(e, r"---(n5 <- .tmp[[6]])---")
  .do(e, r"---(v_D <- S(4))---")
  .do(e, r"---(v_t <- GammaMatrix(m0))---")
  .do(e, r"---(t1 <- gamma_trace(v_t))---")
  expect_equal(.run(e, r"---(qcall(t1, "equals", S(0)))---"), r"---(True)---")
  .do(e, r"---(v_t <- .pymul(.pymul(GammaMatrix(m0), GammaMatrix(m1)), GammaMatrix(m2)))---")
  .do(e, r"---(t1 <- gamma_trace(v_t))---")
  expect_equal(.run(e, r"---(qcall(t1, "equals", S(0)))---"), r"---(True)---")
  .do(e, r"---(v_t <- .pymul(.pymul(GammaMatrix(m0), GammaMatrix(m1)), GammaMatrix((-m0))))---")
  .do(e, r"---(t1 <- gamma_trace(v_t))---")
  expect_equal(.run(e, r"---(qcall(t1, "equals", S(0)))---"), r"---(True)---")
  .do(e, r"---(v_t <- .pymul(.pymul(.pymul(.pymul(GammaMatrix(m0), GammaMatrix(m1)), GammaMatrix(m2)), GammaMatrix(m3)), GammaMatrix(m4)))---")
  .do(e, r"---(t1 <- gamma_trace(v_t))---")
  expect_equal(.run(e, r"---(qcall(t1, "equals", S(0)))---"), r"---(True)---")
  .do(e, r"---(v_t <- .pymul(GammaMatrix(m0), GammaMatrix(m1)))---")
  .do(e, r"---(t1 <- gamma_trace(v_t))---")
  expect_equal(.run(e, r"---(.qis_tensor_eq(t1, .pymul(S(4), .qcallobj(g, m0, m1))))---"), r"---(True)---")
  .do(e, r"---(v_t <- .pymul(.pymul(.pymul(GammaMatrix(m0), GammaMatrix(m1)), GammaMatrix(m2)), GammaMatrix(m3)))---")
  .do(e, r"---(t1 <- gamma_trace(v_t))---")
  .do(e, r"---(t2 <- ((.pymul(.pymul(S(-4), .qcallobj(g, m0, m2)), .qcallobj(g, m1, m3)) + .pymul(.pymul(S(4), .qcallobj(g, m0, m1)), .qcallobj(g, m2, m3))) + .pymul(.pymul(S(4), .qcallobj(g, m0, m3)), .qcallobj(g, m1, m2))))---")
  expect_equal(.run(e, r"---(.qis_tensor_eq(t1, t2))---"), r"---(True)---")
  .do(e, r"---(v_t <- .pymul(.pymul(.pymul(.pymul(.pymul(GammaMatrix(m0), GammaMatrix(m1)), GammaMatrix(m2)), GammaMatrix(m3)), GammaMatrix(m4)), GammaMatrix(m5)))---")
  .do(e, r"---(t1 <- gamma_trace(v_t))---")
  .do(e, r"---(t2 <- .pymul(t1, .qcallobj(g, (-m0), (-m5))))---")
  .do(e, r"---(t2 <- qcall(t2, "contract_metric", g))---")
  expect_equal(.run(e, r"---(.qis_tensor_eq(t2, .pymul(v_D, gamma_trace(.pymul(.pymul(.pymul(GammaMatrix(m1), GammaMatrix(m2)), GammaMatrix(m3)), GammaMatrix(m4))))))---"), r"---(True)---")
  .do(e, r"---(v_t <- .pymul(GammaMatrix(m0), GammaMatrix((-m0))))---")
  .do(e, r"---(t1 <- gamma_trace(v_t))---")
  expect_equal(.run(e, r"---(qcall(t1, "equals", .pymul(S(4), v_D)))---"), r"---(True)---")
  .do(e, r"---(v_t <- .pymul(.pymul(.pymul(GammaMatrix(m0), GammaMatrix(m1)), GammaMatrix((-m0))), GammaMatrix((-m1))))---")
  .do(e, r"---(t1 <- gamma_trace(v_t))---")
  expect_equal(.run(e, r"---(qcall(t1, "equals", (.pymul(S(8), v_D) - .pymul(S(4), (v_D ^ S(2))))))---"), r"---(True)---")
  .do(e, r"---(v_t <- .pymul(.pymul(.pymul(.pymul(.pymul(GammaMatrix(m0), GammaMatrix(m1)), GammaMatrix(m2)), GammaMatrix(m3)), GammaMatrix(m4)), GammaMatrix((-m0))))---")
  .do(e, r"---(t1 <- gamma_trace(v_t))---")
  .do(e, r"---(t2 <- ((.pymul(.pymul(.pymul(S(-4), v_D), .qcallobj(g, m1, m3)), .qcallobj(g, m2, m4)) + .pymul(.pymul(.pymul(S(4), v_D), .qcallobj(g, m1, m2)), .qcallobj(g, m3, m4))) + .pymul(.pymul(.pymul(S(4), v_D), .qcallobj(g, m1, m4)), .qcallobj(g, m2, m3))))---")
  expect_equal(.run(e, r"---(.qis_tensor_eq(t1, t2))---"), r"---(True)---")
  .do(e, r"---(v_t <- .pymul(.pymul(.pymul(.pymul(.pymul(.pymul(.pymul(GammaMatrix((-m5)), GammaMatrix(m0)), GammaMatrix(m1)), GammaMatrix(m2)), GammaMatrix(m3)), GammaMatrix(m4)), GammaMatrix((-m0))), GammaMatrix(m5)))---")
  .do(e, r"---(t1 <- gamma_trace(v_t))---")
  .do(e, r"---(t2 <- .pymul(((.pymul(S(32), v_D) + .pymul(S(4), (((-v_D) + S(4)) ^ S(2)))) - S(64)), ((.pymul(.qcallobj(g, m1, m2), .qcallobj(g, m3, m4)) - .pymul(.qcallobj(g, m1, m3), .qcallobj(g, m2, m4))) + .pymul(.qcallobj(g, m1, m4), .qcallobj(g, m2, m3)))))---")
  expect_equal(.run(e, r"---(.qis_tensor_eq(t1, t2))---"), r"---(True)---")
  .do(e, r"---(v_t <- .pymul(.pymul(.pymul(GammaMatrix(m0), GammaMatrix(m1)), GammaMatrix((-m0))), GammaMatrix(m3)))---")
  .do(e, r"---(t1 <- gamma_trace(v_t))---")
  expect_equal(.run(e, r"---(qcall(t1, "equals", .pymul((.pymul(S(-4), v_D) + S(8)), .qcallobj(g, m1, m3))))---"), r"---(True)---")
  .do(e, r"---(v_t <- .pymul(.pymul(.pymul(.pymul(.pymul(.pymul(.pymul(.pymul(.pymul(.pymul(.pymul(GammaMatrix(m0), GammaMatrix(m1)), GammaMatrix(m2)), GammaMatrix(m3)), GammaMatrix(m4)), GammaMatrix(m5)), GammaMatrix((-m0))), GammaMatrix((-m1))), GammaMatrix((-m2))), GammaMatrix((-m3))), GammaMatrix((-m4))), GammaMatrix((-m5))))---")
  .do(e, r"---(t1 <- gamma_trace(v_t))---")
  expect_equal(.run(e, r"---(qcall(t1, "equals", (((((.pymul(S(-4), (v_D ^ S(6))) + .pymul(S(120), (v_D ^ S(5)))) - .pymul(S(1040), (v_D ^ S(4)))) + .pymul(S(3360), (v_D ^ S(3)))) - .pymul(S(4480), (v_D ^ S(2)))) + .pymul(S(2048), v_D))))---"), r"---(True)---")
  .do(e, r"---(v_t <- .pymul(.pymul(.pymul(.pymul(.pymul(.pymul(.pymul(.pymul(.pymul(.pymul(.pymul(.pymul(.pymul(GammaMatrix(m0), GammaMatrix(m1)), GammaMatrix(n1)), GammaMatrix(m2)), GammaMatrix(n2)), GammaMatrix(m3)), GammaMatrix(m4)), GammaMatrix((-n2))), GammaMatrix((-n1))), GammaMatrix((-m0))), GammaMatrix((-m1))), GammaMatrix((-m2))), GammaMatrix((-m3))), GammaMatrix((-m4))))---")
  .do(e, r"---(t1 <- gamma_trace(v_t))---")
  .do(e, r"---(tresu <- ((((((.pymul(S(-7168), v_D) + .pymul(S(16768), (v_D ^ S(2)))) - .pymul(S(14400), (v_D ^ S(3)))) + .pymul(S(5920), (v_D ^ S(4)))) - .pymul(S(1232), (v_D ^ S(5)))) + .pymul(S(120), (v_D ^ S(6)))) - .pymul(S(4), (v_D ^ S(7)))))---")
  expect_equal(.run(e, r"---(qcall(t1, "equals", tresu))---"), r"---(True)---")
  .do(e, r"---(v_t <- .pymul(.pymul(.pymul(.pymul(.pymul(.pymul(.pymul(.pymul(.pymul(.pymul(.pymul(.pymul(.pymul(GammaMatrix(m0), GammaMatrix(m1)), GammaMatrix(n1)), GammaMatrix(m2)), GammaMatrix(n2)), GammaMatrix(m3)), GammaMatrix(m4)), GammaMatrix(n3)), GammaMatrix(n4)), GammaMatrix((-m0))), GammaMatrix((-m1))), GammaMatrix((-m2))), GammaMatrix((-m3))), GammaMatrix((-m4))))---")
  .do(e, r"---(t1 <- gamma_trace(v_t))---")
  .do(e, r"---(c1 <- (((((.pymul(S(-4), (v_D ^ S(5))) + .pymul(S(120), (v_D ^ S(4)))) - .pymul(S(1200), (v_D ^ S(3)))) + .pymul(S(5280), (v_D ^ S(2)))) - .pymul(S(10560), v_D)) + S(7808)))---")
  .do(e, r"---(c2 <- (((((.pymul(S(-4), (v_D ^ S(5))) + .pymul(S(88), (v_D ^ S(4)))) - .pymul(S(560), (v_D ^ S(3)))) + .pymul(S(1440), (v_D ^ S(2)))) - .pymul(S(1600), v_D)) + S(640)))---")
  expect_equal(.run(e, r"---(.qis_tensor_eq(t1, ((.pymul(.pymul(c1, .qcallobj(g, n1, n4)), .qcallobj(g, n2, n3)) + .pymul(.pymul(c2, .qcallobj(g, n1, n2)), .qcallobj(g, n3, n4))) + .pymul(.pymul((-c1), .qcallobj(g, n1, n3)), .qcallobj(g, n2, n4)))))---"), r"---(True)---")
  .do(e, r"---(.tmp <- tensor_heads("p,q", list(LorentzIndex)))---")
  .do(e, r"---(p <- .tmp[[1]])---")
  .do(e, r"---(v_q <- .tmp[[2]])---")
  .do(e, r"---(ps <- .pymul(.qcallobj(p, m0), GammaMatrix((-m0))))---")
  .do(e, r"---(qs <- .pymul(.qcallobj(v_q, m0), GammaMatrix((-m0))))---")
  .do(e, r"---(p2 <- .pymul(.qcallobj(p, m0), .qcallobj(p, (-m0))))---")
  .do(e, r"---(q2 <- .pymul(.qcallobj(v_q, m0), .qcallobj(v_q, (-m0))))---")
  .do(e, r"---(pq <- .pymul(.qcallobj(p, m0), .qcallobj(v_q, (-m0))))---")
  .do(e, r"---(v_t <- .pymul(.pymul(.pymul(ps, qs), ps), qs))---")
  .do(e, r"---(r <- gamma_trace(v_t))---")
  expect_equal(.run(e, r"---(.qis_tensor_eq(r, (.pymul(.pymul(S(8), pq), pq) - .pymul(.pymul(S(4), p2), q2))))---"), r"---(True)---")
  .do(e, r"---(v_t <- .pymul(.pymul(.pymul(.pymul(.pymul(ps, qs), ps), qs), ps), qs))---")
  .do(e, r"---(r <- gamma_trace(v_t))---")
  expect_equal(.run(e, r"---(.qis_tensor_eq(r, (.pymul(.pymul(.pymul(S(-12), p2), pq), q2) + .pymul(.pymul(.pymul(S(16), pq), pq), pq))))---"), r"---(True)---")
  .do(e, r"---(v_t <- .pymul(.pymul(.pymul(.pymul(.pymul(.pymul(.pymul(ps, qs), ps), qs), ps), qs), ps), qs))---")
  .do(e, r"---(r <- gamma_trace(v_t))---")
  expect_equal(.run(e, r"---(.qis_tensor_eq(r, ((.pymul(.pymul(.pymul(.pymul(S(-32), pq), pq), p2), q2) + .pymul(.pymul(.pymul(.pymul(S(32), pq), pq), pq), pq)) + .pymul(.pymul(.pymul(.pymul(S(4), p2), p2), q2), q2))))---"), r"---(True)---")
  .do(e, r"---(v_t <- .pymul(.pymul(.pymul(.pymul(.pymul(.pymul(S(4), .qcallobj(p, m1)), .qcallobj(p, m0)), .qcallobj(p, (-m0))), .qcallobj(v_q, (-m1))), .qcallobj(v_q, m2)), .qcallobj(v_q, (-m2))))---")
  expect_equal(.run(e, r"---(.qis_tensor_eq(gamma_trace(v_t), v_t))---"), r"---(True)---")
  .do(e, r"---(v_t <- .pymul(.pymul(.pymul(.pymul(.pymul(.pymul(.pymul(ps, ps), ps), ps), ps), ps), ps), ps))---")
  .do(e, r"---(r <- gamma_trace(v_t))---")
  expect_equal(.run(e, r"---(qcall(r, "equals", .pymul(.pymul(.pymul(.pymul(S(4), p2), p2), p2), p2)))---"), r"---(True)---")
})

test_that("test_gamma_matrices: test_bug_13636", {
  e <- .new(character(), 'test_gamma_matrices')
  .do(e, r"---(.tmp <- tensor_heads("pi, ki, pf", list(LorentzIndex)))---")
  .do(e, r"---(v_pi <- .tmp[[1]])---")
  .do(e, r"---(ki <- .tmp[[2]])---")
  .do(e, r"---(pf <- .tmp[[3]])---")
  .do(e, r"---(.tmp <- tensor_indices("i0:5", LorentzIndex))---")
  .do(e, r"---(i0 <- .tmp[[1]])---")
  .do(e, r"---(i1 <- .tmp[[2]])---")
  .do(e, r"---(i2 <- .tmp[[3]])---")
  .do(e, r"---(i3 <- .tmp[[4]])---")
  .do(e, r"---(i4 <- .tmp[[5]])---")
  .do(e, r"---(x <- Symbol("x"))---")
  .do(e, r"---(pis <- .pymul(.qcallobj(v_pi, i2), GammaMatrix((-i2))))---")
  .do(e, r"---(kis <- .pymul(.qcallobj(ki, i3), GammaMatrix((-i3))))---")
  .do(e, r"---(pfs <- .pymul(.qcallobj(pf, i4), GammaMatrix((-i4))))---")
  .do(e, r"---(a <- .pymul(.pymul(.pymul(.pymul(.pymul(.pymul(.pymul(pfs, GammaMatrix(i0)), kis), GammaMatrix(i1)), pis), GammaMatrix((-i1))), kis), GammaMatrix((-i0))))---")
  .do(e, r"---(b <- .pymul(.pymul(.pymul(.pymul(.pymul(.pymul(.pymul(pfs, GammaMatrix(i0)), kis), GammaMatrix(i1)), pis), x), GammaMatrix((-i0))), .qcallobj(v_pi, (-i1))))---")
  .do(e, r"---(ta <- gamma_trace(a))---")
  .do(e, r"---(tb <- gamma_trace(b))---")
  .do(e, r"---(t_a_plus_b <- gamma_trace((a + b)))---")
  expect_equal(.run(e, r"---(ta)---"), r"---(4*((-4)*ki(L_0)*ki(-L_0)*pf(L_1)*pi(-L_1) + 8*ki(L_0)*ki(L_1)*pf(-L_0)*pi(-L_1)))---")
  expect_equal(.run(e, r"---(tb)---"), r"---(-8*x*ki(L_0)*pf(-L_0)*pi(L_1)*pi(-L_1))---")
  expect_equal(.run(e, r"---(t_a_plus_b)---"), r"---(4*((-4)*ki(L_0)*ki(-L_0)*pf(L_1)*pi(-L_1) + 8*ki(L_0)*ki(L_1)*pf(-L_0)*pi(-L_1)) - 8*x*ki(L_0)*pf(-L_0)*pi(L_1)*pi(-L_1))---")
})
