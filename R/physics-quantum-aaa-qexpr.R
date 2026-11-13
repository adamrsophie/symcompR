# Quantum mechanics core, after sympy/physics/quantum/qexpr.py, kind.py,
# dagger.py and hilbert.py.
#
# SymPy's quantum objects form a Python class hierarchy with multiple
# inheritance. Here each class is an entry of .QC holding its bases, its
# method resolution order (C3, as Python) and its attributes: plain values,
# methods f(self, ...), class methods f(cls, ...) marked with .cm() and
# properties f(self) marked with .prop(). An instance is a sym node whose
# type is the class name, whose args are SymPy's args and whose extra$hs is
# the Hilbert space slot of QExpr. A class is passed around by its name.

.QC <- new.env(hash = TRUE)
.qm_state <- new.env()
.qm_state$on <- FALSE

.prop <- function(f) structure(f, qkind = "prop")
.cm <- function(f) structure(f, qkind = "cm")

.qc3 <- function(name, bases) {
  if (!length(bases)) return(name)
  seqs <- c(lapply(bases, function(b) {
    d <- .QC[[b]]
    if (is.null(d)) stop("unknown quantum class ", b)
    d$mro
  }), list(bases))
  out <- name
  repeat {
    seqs <- seqs[lengths(seqs) > 0L]
    if (!length(seqs)) return(out)
    cand <- NULL
    for (s in seqs) {
      h <- s[1L]
      if (!any(vapply(seqs, function(t) h %in% t[-1L], NA))) { cand <- h; break }
    }
    if (is.null(cand)) stop("inconsistent method resolution order for ", name)
    out <- c(out, cand)
    seqs <- lapply(seqs, function(s) if (s[1L] == cand) s[-1L] else s)
  }
}

# Define a class. `node` registers the sym hooks (printing, facts, rebuild,
# doit) for instances; the abstract roots Basic, Expr and Function have none.
.qclass <- function(.name, .bases = "QExpr", ..., pyname = .name, node = TRUE) {
  assign(.name, list(name = .name, pyname = pyname, bases = .bases, mro = .qc3(.name, .bases),
                    attrs = list(...), node = node), envir = .QC)
  if (node) .qhooks(.name)
  invisible(.name)
}

# Add or replace attributes of a class.
.qdef <- function(cls, ...) {
  d <- .QC[[cls]]
  a <- list(...)
  for (nm in names(a)) d$attrs[nm] <- list(a[[nm]])
  assign(cls, d, envir = .QC)
  invisible(cls)
}

# The class of an expression node type, or NULL; the abstract roots (Basic,
# Expr, Function) never own a node, so an applied f(x) is not a quantum class.
.qnode_def <- function(t) { d <- .QC[[t]]; if (is.null(d) || isFALSE(d$node)) NULL else d }

.qis_class <- function(x) is.character(x) && length(x) == 1L && !is.null(.QC[[x]])
.qpyname <- function(cls) { d <- .QC[[cls]]; if (is.null(d)) cls else d$pyname }

# Find attribute `nm` along the MRO of `cls`; list(val) or NULL.
.qfind <- function(cls, nm) {
  d <- .QC[[cls]]
  if (is.null(d)) return(NULL)
  for (c in d$mro) {
    a <- .QC[[c]]$attrs
    if (nm %in% names(a)) return(list(val = a[[nm]], cls = c))
  }
  NULL
}

.qhas <- function(obj, nm) {
  cls <- if (is.character(obj)) obj else obj$type
  !is.null(.qfind(cls, nm))
}

.qerror <- function(name, message, bases = character()) {
  base <- switch(name, QuantumError = character(), HilbertSpaceError = "QuantumError",
                 bases)
  structure(list(message = paste(as.character(message), collapse = ""), call = NULL),
            class = c(name, base, "sympy_error", "error", "condition"))
}
.qstop <- function(name, ...) stop(.qerror(name, paste0(...)))
.qni <- function(...) .qstop("NotImplementedError", ...)
.qtry_ni <- function(expr, default = NULL) tryCatch(expr, NotImplementedError = function(e) default)

# Attribute access: obj.nm. A property is evaluated; a method is returned as a
# function bound to obj.
.qa <- function(obj, nm) {
  cls <- if (is.character(obj)) obj else obj$type
  r <- .qfind(cls, nm)
  if (is.null(r)) .qstop("AttributeError", "'", .qpyname(cls), "' object has no attribute '", nm, "'")
  v <- r$val
  if (is.function(v)) {
    k <- attr(v, "qkind")
    if (identical(k, "prop")) return(v(obj))
    if (identical(k, "cm") && !is.character(obj)) obj <- obj$type
    return(function(...) v(obj, ...))
  }
  v
}

# Method call: obj.nm(...)
.qmeth <- function(.obj, .nm, ...) {
  cls <- if (is.character(.obj)) .obj else .obj$type
  r <- .qfind(cls, .nm)
  if (is.null(r)) .qstop("AttributeError", "'", .qpyname(cls), "' object has no attribute '", .nm, "'")
  f <- r$val
  if (identical(attr(f, "qkind"), "cm") && !is.character(.obj)) .obj <- .obj$type
  f(.obj, ...)
}

# Cls.nm(self, ...): the method as defined for class `cls`
.qsuper <- function(.cls, .nm, .obj, ...) {
  r <- .qfind(.cls, .nm)
  if (is.null(r)) .qstop("AttributeError", "type object '", .qpyname(.cls), "' has no attribute '", .nm, "'")
  f <- r$val
  if (identical(attr(f, "qkind"), "cm") && !is.character(.obj)) .obj <- .obj$type
  f(.obj, ...)
}

.qmro <- function(x) {
  if (is.character(x)) return(.QC[[x]]$mro)
  if (is_sym(x)) { d <- .qnode_def(x$type); if (!is.null(d)) return(d$mro) }
  NULL
}

# isinstance(e, cls); `cls` may be a vector (a tuple of classes)
.qisa <- function(e, cls) {
  if (!is_sym(e)) return(inherits(e, "symmatrix") && any(cls %in% c("Matrix", "MatrixBase", "DenseMatrix", "MutableDenseMatrix")))
  t <- e$type
  d <- .QC[[t]]
  # an applied function is an instance of its own class, e.g. conjugate
  if (!is.null(d) && t != "Function") return(any(cls %in% d$mro))
  if (any(cls == t)) return(TRUE)
  if (t == "Function" && any(cls %in% c("Function", e$name))) return(TRUE)
  if (any(cls == "Basic")) return(TRUE)
  if (any(cls == "Expr")) return(!is_set(e) && !(t %in% c("Tuple", "BooleanAtom", "Str")))
  if (any(cls == "Number")) return(is_number_type(e))
  if (any(cls == "Symbol")) return(t %in% c("Symbol", "Dummy"))
  FALSE
}

# issubclass(cls, base)
.qsubclass <- function(cls, base) .qis_class(cls) && any(base %in% .QC[[cls]]$mro)

.qclassname <- function(x) {
  if (is.character(x)) return("type")
  if (inherits(x, "symmatrix")) return("MutableDenseMatrix")
  if (!is_sym(x)) return(class(x)[1L])
  if (!is.null(.qnode_def(x$type))) return(.qpyname(x$type))
  if (x$type == "Function" && x$name == "Dagger") return("Dagger")
  .sympy_class(x)
}

.qrepr <- function(x) {
  if (is.character(x)) return(paste0("<class 'sympy.physics.quantum.", .qpyname(x), "'>"))
  if (inherits(x, "symmatrix")) return(sstr(x))
  if (is.list(x) && !is_sym(x)) return(paste0("(", paste(vapply(x, .qrepr, ""), collapse = ", "), if (length(x) == 1L) ",", ")"))
  if (!is_sym(x)) return(format(x))
  .str_top(x)
}

# StrPrinter._print_MatrixBase: the rows as a right aligned table
.qmat_str <- function(m) {
  if (m$nrow == 0L || m$ncol == 0L) return(sprintf("Matrix(%d, %d, [])", m$nrow, m$ncol))
  cells <- matrix(vapply(seq_len(m$nrow * m$ncol), function(k) {
    i <- (k - 1L) %/% m$ncol + 1L; j <- (k - 1L) %% m$ncol + 1L
    .str(.mget(m, i, j))
  }, ""), nrow = m$nrow, byrow = TRUE)
  w <- apply(cells, 2L, function(col) max(nchar(col)))
  rows <- vapply(seq_len(m$nrow), function(i)
    paste0("[", paste(mapply(function(z, wd) paste0(strrep(" ", wd - nchar(z)), z), cells[i, ], w), collapse = ", "), "]"), "")
  paste0(if (m$nrow == 1L) "Matrix([" else "Matrix([\n", paste(rows, collapse = ",\n"), "])")
}

# ReprPrinter._print_MatrixBase of an immutable matrix
.qmat_repr <- function(m) {
  rows <- vapply(.mrows(m), function(r) paste0("[", paste(vapply(r, srepr, ""), collapse = ", "), "]"), "")
  paste0("ImmutableDenseMatrix([", paste(rows, collapse = ", "), "])")
}

# ---------------------------------------------------------------------------
# Hooks for instance nodes

.qraw <- function(cls, hs, args) {
  .nc_state$any <- TRUE
  .qm_state$on <- TRUE
  new_sym(cls, args, extra = if (is.null(hs)) NULL else list(hs = hs))
}

.qhooks <- function(name) {
  register_str_printer(name, function(e) .qprint(e, "str"))
  register_latex_printer(name, function(e) .qprint(e, "latex"))
  assign(name, function(e) .qprint(e, "repr"), envir = .srepr_printers)
  assign(name, function(e) .qfacts(e), envir = .facts_hooks)
  register_rebuilder(name, function(e, args) .qfunc(e$type, args))
  register_doit(name, function(e, deep = TRUE, ...) {
    if (!.qhas(e, "doit")) return(e)
    # doit(**hints): a method that does not take a hint ignores it
    h <- list(...)
    fn <- names(formals(.qfind(e$type, "doit")$val))
    if (!("..." %in% fn)) h <- h[names(h) %in% fn]
    do.call(.qmeth, c(list(e, "doit"), h))
  })
}

# self.func(*args)
.qfunc <- function(cls, args, ...) {
  f <- .qfind(cls, "__new__")$val
  if (is.null(f)) return(new_sym(cls, args))
  f(cls, args, ...)
}

.QFACTS <- c("hermitian", "antihermitian", "complex", "real", "extended_real", "positive", "negative",
             "nonnegative", "integer", "zero", "finite", "irrational", "imaginary")

.qfacts <- function(e) {
  # Basic.is_commutative asks the facts again: unknown (None) for a class
  # that does not set it, as Tr
  r <- .qfind(e$type, "is_commutative")
  com <- if (is.null(r) || identical(r$cls, "Basic")) NA else .qa(e, "is_commutative")
  f <- c(commutative = isTRUE(com))
  for (nm in .QFACTS) {
    r <- .qfind(e$type, paste0("is_", nm))
    if (is.null(r) || is.null(r$val)) next
    v <- r$val
    if (is.function(v)) v <- v(e)
    if (is.logical(v) && length(v) == 1L && !is.na(v)) f[nm] <- v
  }
  f
}

# LatexPrinter._print_Tuple
register_latex_printer("Tuple", function(e) {
  a <- vapply(e$args, .latex, "")
  paste0("\\left( ", if (length(a) == 1L) paste0(a, ",") else paste(a, collapse = ", \\  "), "\\right)")
})

# printer._print(x) for printer p in "str", "latex", "repr"
.qpr <- function(p, x) {
  if (is.character(x)) return(x)
  if (is_sym(x) && x$type == "ImmutableMatrix") x <- as_explicit(x)
  if (inherits(x, "symmatrix")) return(switch(p, latex = latex(x), str = .qmat_str(x), .qmat_repr(x)))
  if (is.numeric(x)) x <- as_sym(x)
  switch(p, str = .str(x), latex = .latex(x), repr = srepr(x))
}

.qprint <- function(e, p) {
  nm <- switch(p, str = "_sympystr", latex = "_latex", repr = "_sympyrepr")
  .qmeth(e, nm, p)
}

.qhs <- function(e) {
  if (!is.null(e$extra$hs)) return(e$extra$hs)
  if (.qhas(e, "hilbert_space")) return(.qa(e, "hilbert_space"))
  HilbertSpace()
}

# ---------------------------------------------------------------------------
# Small helpers on SymPy objects

.qS <- function(x) {
  if (is_sym(x) || inherits(x, "symmatrix")) return(x)
  if (is.character(x)) return(Symbol(x))
  as_sym(x)
}

# complex(e) fails: e holds a quantum object without a numeric value (CG,
# WignerD, states), so Add.as_terms keeps it a generator
.qno_value <- function(e) {
  if (!is_sym(e)) return(FALSE)
  if (!is.null(.qnode_def(e$type)) && is.null(.neval_hooks[[e$type]])) return(TRUE)
  for (a in e$args) if (.qno_value(a)) return(TRUE)
  FALSE
}

# expr.args_cnc(): the factors before the first noncommutative one, and the
# rest; a leading negative number other than -1 splits off -1
.args_cnc <- function(e, split_1 = TRUE) {
  args <- if (e$type == "Mul") .sympy_args(e) else list(e)
  k <- which(!vapply(args, .is_comm, NA))
  if (length(k)) {
    c <- args[seq_len(k[1L] - 1L)]
    nc <- args[k[1L]:length(args)]
  } else {
    c <- args
    nc <- list()
  }
  if (length(c) && split_1 && is_number_type(c[[1L]]) && num_sign(c[[1L]]) < 0 && !is_negone_num(c[[1L]]))
    c <- c(list(S_NegativeOne, num_neg(c[[1L]])), c[-1L])
  list(c = c, nc = nc)
}

# Mul._from_args(args)
.qmul_raw <- function(args) {
  if (!length(args)) return(S_One)
  if (length(args) == 1L) return(args[[1L]])
  new_sym("Mul", args)
}

# Python truthiness of a SymPy object (Basic is always true but zero numbers)
.qtruthy <- function(x) {
  if (is.null(x)) return(FALSE)
  if (is_sym(x) && x$type == "BooleanAtom") return(identical(x$name, "true"))
  if (is_sym(x)) return(!(is_number_type(x) && num_is_zero(x)))
  if (is.logical(x)) return(isTRUE(x))
  if (is.numeric(x)) return(length(x) && x != 0)
  length(x) > 0L
}

.qeq <- function(a, b) {
  if (is.null(a) || is.null(b)) return(is.null(a) && is.null(b))
  if (is.character(a) || is.character(b)) return(identical(a, b))
  if (inherits(a, "symmatrix") || inherits(b, "symmatrix"))
    return(inherits(a, "symmatrix") && inherits(b, "symmatrix") && a$nrow == b$nrow && a$ncol == b$ncol &&
             identical(vapply(a$data, .k, ""), vapply(b$data, .k, "")))
  # Python tuple ==: same length and equal items
  la <- is.list(a) && !is_sym(a); lb <- is.list(b) && !is_sym(b)
  if (la || lb)
    return(la && lb && length(a) == length(b) && all(vapply(seq_along(a), function(i) .qeq(a[[i]], b[[i]]), NA)))
  a <- as_sym(a); b <- as_sym(b)
  if (identical(a$key, b$key)) return(TRUE)
  # Python ==: numbers compare by value
  if (is_number_type(a) && is_number_type(b)) return(num_cmp(a, b) == 0)
  FALSE
}

.qint <- function(x) {
  if (is_sym(x)) return(as.numeric(gmp::numerator(x$value)))
  as.numeric(x)
}

# product that is a matrix product for two matrices (SymPy's *)
.qtimes <- function(a, b) {
  am <- inherits(a, "symmatrix"); bm <- inherits(b, "symmatrix")
  if (am && bm) return(.matmul(a, b))
  if (am) return(.mat_map(a, function(z) .mul(list(z, as_sym(b)))))
  if (bm) return(.mat_map(b, function(z) .mul(list(as_sym(a), z))))
  .mul(list(as_sym(a), as_sym(b)))
}
.qplus <- function(a, b) {
  if (inherits(a, "symmatrix") && inherits(b, "symmatrix")) return(a + b)
  if (inherits(a, "symmatrix") || inherits(b, "symmatrix")) stop(.qerror("TypeError", "cannot add matrix and scalar"))
  .add(list(as_sym(a), as_sym(b)))
}

# ---------------------------------------------------------------------------
# Kinds (sympy/physics/quantum/kind.py and the core kind dispatchers)

.KIND_MUL <- c("NumberKind|KetKind" = "KetKind", "NumberKind|BraKind" = "BraKind",
               "NumberKind|OperatorKind" = "OperatorKind", "NumberKind|MatrixKind" = "MatrixKind",
               "MatrixKind|MatrixKind" = "MatrixKind")
.KIND_TP <- c("OperatorKind|OperatorKind" = "OperatorKind", "OperatorKind|NumberKind" = "OperatorKind",
              "KetKind|KetKind" = "KetKind", "NumberKind|KetKind" = "KetKind",
              "BraKind|BraKind" = "BraKind", "NumberKind|BraKind" = "BraKind")
.KIND_COMM <- c("OperatorKind|OperatorKind" = "OperatorKind")

# KindDispatcher(commutative=True): the distinct kinds reduced pairwise
.kind_dispatch <- function(kinds, table) {
  u <- unique(kinds)
  if (length(u) == 1L) return(u)
  res <- u[1L]
  for (k in u[-1L]) {
    r <- table[paste(res, k, sep = "|")]
    if (is.na(r)) r <- table[paste(k, res, sep = "|")]
    res <- if (is.na(r)) "UndefinedKind" else unname(r)
  }
  res
}

.qkind <- function(e) {
  if (inherits(e, "symmatrix")) return("MatrixKind")
  if (is.character(e)) return("UndefinedKind")
  if (!is.null(.qnode_def(e$type))) {
    if (.qhas(e, "kind")) return(.qa(e, "kind"))
    return(if (.is_comm(e)) "NumberKind" else "UndefinedKind")
  }
  switch(e$type,
         Mul = .kind_dispatch(vapply(e$args, .qkind, ""), .KIND_MUL),
         Add = { u <- unique(vapply(e$args, .qkind, "")); if (length(u) == 1L) u else "UndefinedKind" },
         Pow = if (.qkind(e$args[[2L]]) == "NumberKind") .qkind(e$args[[1L]]) else "UndefinedKind",
         Function = if (e$name == "Dagger") .qkind(e$args[[1L]]) else if (.is_comm(e)) "NumberKind" else "UndefinedKind",
         if (.is_matexpr(e)) "MatrixKind" else if (.is_comm(e)) "NumberKind" else "UndefinedKind")
}

# ---------------------------------------------------------------------------
# Root classes

.qclass("Basic", character(), node = FALSE, is_commutative = .prop(function(self) .is_comm(self)))
.qclass("Expr", "Basic", node = FALSE)
.qclass("Function", "Expr", node = FALSE)

# ---------------------------------------------------------------------------
# QExpr

# _qsympify_sequence: nested sequences become Tuples, strings Symbols
.qsympify_seq <- function(seq) lapply(seq, .qsympify_helper)
.qsympify_helper <- function(x) {
  if (inherits(x, "symmatrix")) return(x)
  if (is_sym(x)) return(x)
  if (is.character(x) && length(x) == 1L) return(Symbol(x))
  if (is.list(x) || (is.atomic(x) && length(x) > 1L)) return(do.call(Tuple, lapply(as.list(x), .qsympify_helper)))
  as_sym(x)
}

.qclass("QExpr", "Expr",
  is_commutative = FALSE,
  `_label_separator` = "",
  `__new__` = .cm(function(cls, args, ...) {
    args <- .qmeth(cls, "_eval_args", args, ...)
    if (!length(args)) args <- .qmeth(cls, "_eval_args", as.list(.qmeth(cls, "default_args")), ...)
    .qraw(cls, .qmeth(cls, "_eval_hilbert_space", args), args)
  }),
  # _new_rawargs(hilbert_space, *args); one unnamed list is the args themselves
  `_new_rawargs` = .cm(function(cls, hilbert_space, ...) {
    a <- list(...)
    if (length(a) == 1L && is.list(a[[1L]]) && !is_sym(a[[1L]])) a <- a[[1L]]
    .qraw(cls, hilbert_space, a)
  }),
  label = .prop(function(self) {
    if (!length(self$args)) return(.qmeth(self$type, "_eval_args", as.list(.qmeth(self$type, "default_args"))))
    self$args
  }),
  is_symbolic = TRUE,
  hilbert_space = .prop(function(self) if (!is.null(self$extra$hs)) self$extra$hs else HilbertSpace()),
  default_args = .cm(function(cls) .qni("No default arguments for this class!")),
  `_eval_adjoint` = function(self) {
    obj <- .expr_eval_adjoint(self)
    if (is.null(obj)) obj <- fnode("Dagger", list(self))
    obj
  },
  `_eval_args` = .cm(function(cls, args, ...) .qsympify_seq(args)),
  `_eval_hilbert_space` = .cm(function(cls, args) HilbertSpace()),
  `_print_sequence` = function(self, seq, sep, p) paste(vapply(seq, function(i) .qpr(p, i), ""), collapse = sep),
  `_print_label` = function(self, p) .qmeth(self, "_print_sequence", .qa(self, "label"), .qa(self, "_label_separator"), p),
  `_print_label_repr` = function(self, p) .qmeth(self, "_print_sequence", .qa(self, "label"), ",", p),
  `_print_label_latex` = function(self, p) .qmeth(self, "_print_sequence", .qa(self, "label"), .qa(self, "_label_separator"), p),
  `_print_contents` = function(self, p) .qmeth(self, "_print_label", p),
  `_print_contents_latex` = function(self, p) .qmeth(self, "_print_label_latex", p),
  `_sympystr` = function(self, p) .qmeth(self, "_print_contents", p),
  `_sympyrepr` = function(self, p) paste0(.qpyname(self$type), "(", .qmeth(self, "_print_label_repr", p), ")"),
  `_latex` = function(self, p) .qmeth(self, "_print_contents_latex", p),
  `_represent_default_basis` = function(self, ...) .qni("This object does not have a default basis"),
  `_represent` = function(self, basis = NULL, ...) {
    result <- if (is.null(basis)) .qmeth(self, "_represent_default_basis", ...)
              else .qdispatch(self, "_represent", basis, ...)
    .qformat_represent(result)
  }
)

.qformat_represent <- function(result) {
  if (inherits(result, "symmatrix")) return(result)
  result
}

.qdispatch <- function(self, basename, arg, ...) {
  nm <- paste0(basename, "_", .qclassname(arg))
  if (.qhas(self, nm)) {
    r <- .qmeth(self, nm, arg, ...)
    if (!is.null(r)) return(r)
  }
  .qni(.qpyname(self$type), ".", basename, " cannot handle: ", .qrepr(arg))
}

# ---------------------------------------------------------------------------
# Dagger and adjoint (sympy/functions/elementary/complexes.py adjoint, whose
# printed name is Dagger once sympy.physics.quantum is loaded)

def_function("Dagger",
             eval = function(args) .adjoint_eval(args[[1L]]),
             facts = function(args, afs) {
               f <- afs[[1L]]
               out <- c(commutative = !identical(f[["commutative"]], FALSE))
               for (nm in c("hermitian", "antihermitian", "real", "finite"))
                 if (!is.na(f[[nm]])) out[nm] <- f[[nm]]
               out
             },
             latex = function(a) {
               s <- a[1L]
               if (grepl("^", s, fixed = TRUE)) paste0("\\left(", s, "\\right)^{\\dagger}") else paste0(s, "^{\\dagger}")
             })

.adjoint_eval <- function(arg) {
  obj <- .eval_adjoint(arg)
  if (!is.null(obj)) return(obj)
  obj <- .eval_transpose(arg)
  if (!is.null(obj)) return(fcall("conjugate", list(obj)))
  NULL
}

# adjoint(e)
.adjoint <- function(e) {
  if (inherits(e, "symmatrix")) return(t(.mat_map(e, function(z) fcall("conjugate", list(z)))))
  fcall("Dagger", list(e))
}

.eval_adjoint <- function(e) {
  if (!is.null(.qnode_def(e$type)) && .qhas(e, "_eval_adjoint")) return(.qmeth(e, "_eval_adjoint"))
  switch(e$type,
         Add = .add(lapply(.sympy_args(e), .adjoint)),
         Mul = .mul(lapply(rev(.sympy_args(e)), .adjoint)),
         Pow = {
           b <- e$args[[1L]]; x <- e$args[[2L]]
           i <- ask_fact(x, "integer"); p <- ask_fact(b, "positive")
           if (isTRUE(i)) return(.pow(.adjoint(b), x))
           if (isTRUE(p)) return(.pow(b, .adjoint(x)))
           if (identical(i, FALSE) && identical(p, FALSE)) {
             ex <- expand_complex(e)
             if (!identical(ex$key, e$key)) return(.adjoint(ex))
           }
           NULL
         },
         Function = switch(e$name,
                           Dagger = e$args[[1L]],
                           conjugate = .transpose_fn(e$args[[1L]]),
                           transpose = fcall("conjugate", list(e$args[[1L]])),
                           exp = fcall("exp", list(.adjoint(e$args[[1L]]))),
                           .expr_eval_adjoint(e)),
         .expr_eval_adjoint(e))
}

# Expr._eval_adjoint
.expr_eval_adjoint <- function(e) {
  if (isTRUE(ask_fact(e, "hermitian"))) return(e)
  if (isTRUE(ask_fact(e, "antihermitian"))) return(.neg(e))
  obj <- .eval_conjugate(e)
  if (!is.null(obj)) return(.transpose_fn(obj))
  obj <- .eval_transpose(e)
  if (!is.null(obj)) return(fcall("conjugate", list(obj)))
  NULL
}

.eval_conjugate <- function(e) {
  if (!is.null(.qnode_def(e$type))) {
    if (.qhas(e, "_eval_conjugate")) return(.qmeth(e, "_eval_conjugate"))
    if (isTRUE(ask_fact(e, "extended_real"))) return(e)
    if (isTRUE(ask_fact(e, "imaginary"))) return(.neg(e))
    return(NULL)
  }
  if (!.is_comm(e)) return(NULL)
  r <- fcall("conjugate", list(e))
  if (r$type == "Function" && r$name == "conjugate" && identical(r$args[[1L]]$key, e$key)) NULL else r
}

.eval_transpose <- function(e) {
  if (!is.null(.qnode_def(e$type)) && .qhas(e, "_eval_transpose")) return(.qmeth(e, "_eval_transpose"))
  if (.is_comm(e)) return(e)
  if (isTRUE(ask_fact(e, "hermitian"))) return(fcall("conjugate", list(e)))
  if (isTRUE(ask_fact(e, "antihermitian"))) return(.neg(fcall("conjugate", list(e))))
  NULL
}

.transpose_fn <- function(e) {
  r <- .eval_transpose(e)
  if (!is.null(r)) r else fnode("transpose", list(e))
}

#' Hermitian conjugate
#'
#' `Dagger(arg)` takes the Hermitian conjugate of a quantum object, an
#' operator expression, a matrix or a scalar, as SymPy's
#' `sympy.physics.quantum.Dagger`: kets become bras, products reverse their
#' order and numbers are conjugated. With `evaluate = FALSE` the conjugate
#' stays unevaluated.
#'
#' @param arg An expression or a matrix.
#' @param evaluate Evaluate the conjugate.
#' @return An expression or a matrix.
#' @examples
#' A <- Operator("A"); B <- Operator("B")
#' Dagger(A * B)
#' Dagger(Ket("psi"))
#' Dagger(sym_I * A)
#' @export
Dagger <- function(arg, evaluate = TRUE) {
  if (inherits(arg, "symmatrix")) {
    if (evaluate) return(t(.mat_map(arg, function(z) fcall("conjugate", list(z)))))
  }
  arg <- .qS(arg)
  if (evaluate) return(.adjoint(arg))
  fnode("Dagger", list(arg))
}

# ---------------------------------------------------------------------------
# Hilbert spaces (sympy/physics/quantum/hilbert.py)

.HS_TYPES <- c("HilbertSpace", "ComplexSpace", "L2", "FockSpace", "TensorProductHilbertSpace",
               "DirectSumHilbertSpace", "TensorPowerHilbertSpace")
.is_hs <- function(x) is_sym(x) && x$type %in% .HS_TYPES

.hs_raw <- function(cls, args = list()) new_sym(cls, args)

.qclass("HilbertSpace", "Basic",
  `__new__` = .cm(function(cls, args) .hs_raw("HilbertSpace")),
  is_commutative = TRUE,
  dimension = .prop(function(self) .qni("This Hilbert space has no dimension.")),
  `_sympystr` = function(self, p) "H",
  `_sympyrepr` = function(self, p) "HilbertSpace()",
  `_latex` = function(self, p) "\\mathcal{H}")

.qclass("ComplexSpace", "HilbertSpace",
  `__new__` = .cm(function(cls, args) {
    dimension <- as_sym(args[[1L]])
    ats <- atoms(dimension)
    if (length(ats) == 1L) {
      if (!((dimension$type == "Integer" && num_sign(dimension) > 0) || identical(dimension$key, S_Infinity$key) ||
            dimension$type %in% c("Symbol", "Dummy")))
        .qstop("TypeError", "The dimension of a ComplexSpace can onlybe a positive integer, oo, or a Symbol: ", .qrepr(dimension))
    } else {
      for (dim in ats) if (!(dim$type == "Integer" || identical(dim$key, S_Infinity$key) || dim$type %in% c("Symbol", "Dummy")))
        .qstop("TypeError", "The dimension of a ComplexSpace can only contain integers, oo, or a Symbol: ", .qrepr(dim))
    }
    .hs_raw("ComplexSpace", list(dimension))
  }),
  dimension = .prop(function(self) self$args[[1L]]),
  `_sympyrepr` = function(self, p) paste0("ComplexSpace(", .qpr(p, self$args[[1L]]), ")"),
  `_sympystr` = function(self, p) paste0("C(", .qpr(p, self$args[[1L]]), ")"),
  `_latex` = function(self, p) paste0("\\mathcal{C}^{", .qpr(p, self$args[[1L]]), "}"))

.qclass("L2", "HilbertSpace",
  `__new__` = .cm(function(cls, args) {
    interval <- args[[1L]]
    if (!(is_sym(interval) && interval$type == "Interval"))
      .qstop("TypeError", "L2 interval must be an Interval instance: ", .qrepr(interval))
    .hs_raw("L2", list(interval))
  }),
  dimension = .prop(function(self) S_Infinity),
  interval = .prop(function(self) self$args[[1L]]),
  `_sympyrepr` = function(self, p) paste0("L2(", .qpr(p, self$args[[1L]]), ")"),
  `_sympystr` = function(self, p) paste0("L2(", .qpr(p, self$args[[1L]]), ")"),
  `_latex` = function(self, p) paste0("{\\mathcal{L}^2}\\left( ", .qpr(p, self$args[[1L]]), " \\right)"))

.qclass("FockSpace", "HilbertSpace",
  `__new__` = .cm(function(cls, args) .hs_raw("FockSpace")),
  dimension = .prop(function(self) S_Infinity),
  `_sympyrepr` = function(self, p) "FockSpace()",
  `_sympystr` = function(self, p) "F",
  `_latex` = function(self, p) "\\mathcal{F}")

.hs_spaces_printer <- function(self, p) vapply(self$args, function(a) {
  s <- .qpr(p, a)
  if (a$type == "DirectSumHilbertSpace") paste0("(", s, ")") else s
}, "")

.qclass("TensorProductHilbertSpace", "HilbertSpace",
  `__new__` = .cm(function(cls, args) {
    new_args <- list()
    recall <- FALSE
    for (arg in args) {
      if (is_sym(arg) && arg$type == "TensorProductHilbertSpace") {
        new_args <- c(new_args, arg$args)
        recall <- TRUE
      } else if (.is_hs(arg)) {
        new_args[[length(new_args) + 1L]] <- arg
      } else {
        .qstop("TypeError", "Hilbert spaces can only be multiplied by                 other Hilbert spaces: ", .qrepr(arg))
      }
    }
    is_pow <- function(z) z$type == "TensorPowerHilbertSpace"
    comb_args <- list()
    prev_arg <- NULL
    for (new_arg in new_args) {
      if (!is.null(prev_arg)) {
        if (is_pow(new_arg) && is_pow(prev_arg) && identical(new_arg$args[[1L]]$key, prev_arg$args[[1L]]$key))
          prev_arg <- TensorPowerHilbertSpace(new_arg$args[[1L]], .add(list(new_arg$args[[2L]], prev_arg$args[[2L]])))
        else if (is_pow(new_arg) && identical(new_arg$args[[1L]]$key, prev_arg$key))
          prev_arg <- TensorPowerHilbertSpace(prev_arg, .add(list(new_arg$args[[2L]], S_One)))
        else if (is_pow(prev_arg) && identical(new_arg$key, prev_arg$args[[1L]]$key))
          prev_arg <- TensorPowerHilbertSpace(new_arg, .add(list(prev_arg$args[[2L]], S_One)))
        else if (identical(new_arg$key, prev_arg$key))
          prev_arg <- TensorPowerHilbertSpace(new_arg, S_Two)
        else {
          comb_args[[length(comb_args) + 1L]] <- prev_arg
          prev_arg <- new_arg
        }
      } else {
        prev_arg <- new_arg
      }
    }
    comb_args[[length(comb_args) + 1L]] <- prev_arg
    if (recall) return(do.call(TensorProductHilbertSpace, comb_args))
    if (length(comb_args) == 1L) return(TensorPowerHilbertSpace(comb_args[[1L]]$args[[1L]], comb_args[[1L]]$args[[2L]]))
    .hs_raw("TensorProductHilbertSpace", args)
  }),
  dimension = .prop(function(self) {
    dl <- lapply(self$args, function(a) .qa(a, "dimension"))
    if (any(vapply(dl, function(d) identical(d$key, S_Infinity$key), NA))) return(S_Infinity)
    .mul(dl)
  }),
  spaces = .prop(function(self) self$args),
  `_sympyrepr` = function(self, p) paste0("TensorProductHilbertSpace(", paste(.hs_spaces_printer(self, p), collapse = ","), ")"),
  `_sympystr` = function(self, p) paste(.hs_spaces_printer(self, p), collapse = "*"),
  `_latex` = function(self, p) {
    s <- vapply(self$args, function(a) {
      x <- .qpr(p, a)
      if (a$type %in% c("DirectSumHilbertSpace", "TensorProductHilbertSpace")) paste0("\\left(", x, "\\right)") else x
    }, "")
    paste(s, collapse = "\\otimes ")
  })

.qclass("DirectSumHilbertSpace", "HilbertSpace",
  `__new__` = .cm(function(cls, args) {
    new_args <- list()
    recall <- FALSE
    for (arg in args) {
      if (is_sym(arg) && arg$type == "DirectSumHilbertSpace") {
        new_args <- c(new_args, arg$args)
        recall <- TRUE
      } else if (.is_hs(arg)) {
        new_args[[length(new_args) + 1L]] <- arg
      } else {
        .qstop("TypeError", "Hilbert spaces can only be summed with other                 Hilbert spaces: ", .qrepr(arg))
      }
    }
    if (recall) return(do.call(DirectSumHilbertSpace, new_args))
    .hs_raw("DirectSumHilbertSpace", args)
  }),
  dimension = .prop(function(self) {
    dl <- lapply(self$args, function(a) .qa(a, "dimension"))
    if (any(vapply(dl, function(d) identical(d$key, S_Infinity$key), NA))) return(S_Infinity)
    .add(dl)
  }),
  spaces = .prop(function(self) self$args),
  `_sympyrepr` = function(self, p) paste0("DirectSumHilbertSpace(", paste(vapply(self$args, function(a) .qpr(p, a), ""), collapse = ","), ")"),
  `_sympystr` = function(self, p) paste(vapply(self$args, function(a) .qpr(p, a), ""), collapse = "+"),
  `_latex` = function(self, p) {
    s <- vapply(self$args, function(a) {
      x <- .qpr(p, a)
      if (a$type %in% c("DirectSumHilbertSpace", "TensorProductHilbertSpace")) paste0("\\left(", x, "\\right)") else x
    }, "")
    paste(s, collapse = "\\oplus ")
  })

.qclass("TensorPowerHilbertSpace", "HilbertSpace",
  `__new__` = .cm(function(cls, args) {
    b <- args[[1L]]
    ex <- as_sym(args[[2L]])
    if (is_one_num(ex)) return(b)
    if (is_zero_num(ex)) return(S_One)
    ats <- atoms(ex)
    if (length(ats) == 1L) {
      if (!((ex$type == "Integer" && num_sign(ex) >= 0) || ex$type %in% c("Symbol", "Dummy", "Wild")))
        .qstop("ValueError", "Hilbert spaces can only be raised to                 positive integers or Symbols: ", .qrepr(ex))
    } else {
      for (power in ats) if (!(power$type == "Integer" || power$type %in% c("Symbol", "Dummy", "Wild")))
        .qstop("ValueError", "Tensor powers can only contain integers                     or Symbols: ", .qrepr(power))
    }
    .hs_raw("TensorPowerHilbertSpace", list(b, ex))
  }),
  base = .prop(function(self) self$args[[1L]]),
  exp = .prop(function(self) self$args[[2L]]),
  dimension = .prop(function(self) {
    d <- .qa(self$args[[1L]], "dimension")
    if (identical(d$key, S_Infinity$key)) S_Infinity else .pow(d, self$args[[2L]])
  }),
  `_sympyrepr` = function(self, p) paste0("TensorPowerHilbertSpace(", .qpr(p, self$args[[1L]]), ",", .qpr(p, self$args[[2L]]), ")"),
  `_sympystr` = function(self, p) paste0(.qpr(p, self$args[[1L]]), "^", .qpr(p, self$args[[2L]])),
  `_latex` = function(self, p) paste0("{", .qpr(p, self$args[[1L]]), "}^{\\otimes ", .qpr(p, self$args[[2L]]), "}"))

for (.t in .HS_TYPES) assign(.t, function(e) c(commutative = TRUE), envir = .facts_hooks)
# ExprWithLimits: obj.is_commutative = function.is_commutative, so a sum of
# states stays noncommutative
for (.t in c("Sum", "Product", "Integral"))
  assign(.t, function(e) c(commutative = !identical(ask_fact(e$args[[1L]], "commutative"), FALSE)), envir = .facts_hooks)
rm(.t)

#' Hilbert spaces
#'
#' The Hilbert spaces of `sympy.physics.quantum.hilbert`: the abstract
#' `HilbertSpace()`, the finite or infinite `ComplexSpace(dimension)`, the
#' square integrable functions `L2(interval)` on an [Interval()], the Fock
#' space `FockSpace()`, and the tensor products, direct sums and tensor powers
#' built with `*`, `+` and `^` or with the constructors.
#' `hs_contains(space, obj)` is SymPy's `obj in space`.
#'
#' @param dimension A positive integer, `oo` or a symbol.
#' @param interval An [Interval()].
#' @param ... Hilbert spaces.
#' @param base,exp A Hilbert space and its power.
#' @param space A Hilbert space.
#' @param obj A quantum object.
#' @return A Hilbert space, or `TRUE`/`FALSE` for `hs_contains()`.
#' @examples
#' ComplexSpace(2) * FockSpace()
#' ComplexSpace(2)^3
#' qattr(ComplexSpace(2) + FockSpace(), "dimension")
#' @name hilbert-spaces
NULL

#' @rdname hilbert-spaces
#' @export
HilbertSpace <- function() .qfunc("HilbertSpace", list())
#' @rdname hilbert-spaces
#' @export
ComplexSpace <- function(dimension) .qfunc("ComplexSpace", list(dimension))
#' @rdname hilbert-spaces
#' @export
L2 <- function(interval) .qfunc("L2", list(interval))
#' @rdname hilbert-spaces
#' @export
FockSpace <- function() .qfunc("FockSpace", list())
#' @rdname hilbert-spaces
#' @export
TensorProductHilbertSpace <- function(...) .qfunc("TensorProductHilbertSpace", list(...))
#' @rdname hilbert-spaces
#' @export
DirectSumHilbertSpace <- function(...) .qfunc("DirectSumHilbertSpace", list(...))
#' @rdname hilbert-spaces
#' @export
TensorPowerHilbertSpace <- function(base, exp) .qfunc("TensorPowerHilbertSpace", list(base, exp))
#' @rdname hilbert-spaces
#' @export
hs_contains <- function(space, obj) identical(.qhs(obj)$type, space$type)

.hs_ops <- function(gen, e1, e2) {
  if (!(.is_hs(e1) || .is_hs(e2))) return(NULL)
  switch(gen,
         "*" = TensorProductHilbertSpace(e1, e2),
         "+" = DirectSumHilbertSpace(e1, e2),
         "^" = TensorPowerHilbertSpace(e1, e2),
         "==" = .qeq(e1, e2),
         "!=" = !.qeq(e1, e2),
         NULL)
}

local({
  old <- .ops_hooks$dispatch
  .ops_hooks$dispatch <- function(gen, e1, e2) {
    r <- .hs_ops(gen, e1, e2)
    if (!is.null(r)) return(r)
    if (is.null(old)) NULL else old(gen, e1, e2)
  }
})

# ---------------------------------------------------------------------------
# Generic access to quantum objects from R

#' Attributes and methods of quantum objects
#'
#' Quantum objects keep the Python attributes and methods of their SymPy
#' classes. `qattr(obj, name)` is SymPy's `obj.name` (a property is
#' evaluated), `qcall(obj, name, ...)` calls the method `obj.name(...)`,
#' `qisinstance(obj, cls)` is `isinstance(obj, cls)` with the class given by
#' its name (or a vector of names), and `qclass(obj)` is the class name.
#' For other expressions `qattr()` knows `args`, `is_commutative` and
#' `free_symbols`.
#'
#' @param obj A quantum object, a class name or an expression.
#' @param name An attribute or method name.
#' @param cls A class name or a character vector of names.
#' @param ... Arguments of the method.
#' @return The attribute value, the method result, or a logical.
#' @examples
#' A <- Operator("A")
#' qattr(A, "label")
#' qisinstance(A, "QExpr")
#' qcall(HermitianOperator("H"), "inverse")
#' @name quantum-attributes
NULL

#' @rdname quantum-attributes
#' @export
qattr <- function(obj, name) {
  if (is_sym(obj) && is.null(.qnode_def(obj$type))) {
    if (.is_kd(obj) && !is.null(r <- .kd_attr(obj, name))) return(r)
    return(switch(name,
                  args = if (obj$type %in% c("Add", "Mul")) .sympy_args(obj) else obj$args,
                  func = .qclassname(obj),
                  is_commutative = .is_comm(obj),
                  free_symbols = free_symbols(obj),
                  kind = .qkind(obj),
                  label = .qstop("AttributeError", "'", .qclassname(obj), "' object has no attribute 'label'"),
                  .qstop("AttributeError", "'", .qclassname(obj), "' object has no attribute '", name, "'")))
  }
  if (name == "args" && is_sym(obj)) return(obj$args)
  if (name == "hilbert_space" && is_sym(obj)) return(.qhs(obj))
  .qa(obj, name)
}

#' @rdname quantum-attributes
#' @export
qcall <- function(obj, name, ...) .qmeth(obj, name, ...)

#' @rdname quantum-attributes
#' @export
qisinstance <- function(obj, cls) .qisa(obj, cls)

#' @rdname quantum-attributes
#' @export
qclass <- function(obj) .qclassname(obj)
