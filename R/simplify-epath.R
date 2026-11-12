# Expression paths: a port of sympy/simplify/epathtools.py (EPath, epath).

# The SymPy class ancestry (__mro__) of each class, used by type selectors
# such as "/Symbol" or "/TrigonometricFunction". ~F, ~A, ~E and ~B stand for
# the common tails of function, atomic, expression and plain Basic classes.
.EPATH_MRO_TABLE <- c(
  "Abs,~F", "AccumulationBounds,~E", "acos,InverseTrigonometricFunction,~F",
  "acosh,InverseHyperbolicFunction,~F", "acot,InverseTrigonometricFunction,~F",
  "acoth,InverseHyperbolicFunction,~F", "acsc,InverseTrigonometricFunction,~F",
  "acsch,InverseHyperbolicFunction,~F", "Add,Expr,AssocOp,Basic,Printable,EvalfMixin,object",
  "airyai,AiryBase,~F", "airyaiprime,AiryBase,~F", "airybi,AiryBase,~F",
  "airybiprime,AiryBase,~F", "And,LatticeOp,AssocOp,BooleanFunction,Application,Boolean,~B",
  "andre,~F", "appellf1,~F", "arg,~F", "asec,InverseTrigonometricFunction,~F",
  "asech,InverseHyperbolicFunction,~F", "asin,InverseTrigonometricFunction,~F",
  "asinh,InverseHyperbolicFunction,~F", "assoc_laguerre,OrthogonalPolynomial,~F",
  "assoc_legendre,~F", "atan,InverseTrigonometricFunction,~F",
  "atan2,InverseTrigonometricFunction,~F", "atanh,InverseHyperbolicFunction,~F", "bell,~F",
  "bernoulli,~F", "besseli,BesselBase,~F", "besselj,BesselBase,~F", "besselk,BesselBase,~F",
  "bessely,BesselBase,~F", "beta,~F", "betainc,~F", "betainc_regularized,~F",
  "binomial,CombinatorialFunction,~F", "BooleanFalse,BooleanAtom,Boolean,~B",
  "BooleanTrue,BooleanAtom,Boolean,~B", "Catalan,NumberSymbol,~A", "catalan,~F",
  "ceiling,RoundFunction,~F", "chebyshevt,OrthogonalPolynomial,~F",
  "chebyshevu,OrthogonalPolynomial,~F", "Chi,TrigonometricIntegral,~F",
  "Ci,TrigonometricIntegral,~F", "Complement,Set,Basic,Printable,EvalfMixin,object",
  "ComplexInfinity,~A", "conjugate,~F", "cos,TrigonometricFunction,~F",
  "cosh,HyperbolicFunction,~F", "cot,TrigonometricFunction,~F", "coth,HyperbolicFunction,~F",
  "csc,ReciprocalTrigonometricFunction,TrigonometricFunction,~F",
  "csch,ReciprocalHyperbolicFunction,HyperbolicFunction,~F", "Derivative,~E", "Determinant,~E",
  "DiracDelta,~F", "dirichlet_eta,~F", "divisor_sigma,~F",
  "Dummy,Symbol,AtomicExpr,Atom,Expr,Boolean,Basic,Printable,EvalfMixin,object", "Ei,~F",
  "elliptic_e,~F", "elliptic_f,~F", "elliptic_k,~F", "elliptic_pi,~F",
  "Equality,Relational,Boolean,Basic,Printable,EvalfMixin,object",
  "Equivalent,BooleanFunction,Application,Boolean,~B", "erf,~F", "erf2inv,~F", "erfc,~F",
  "erfcinv,~F", "erfi,~F", "erfinv,~F", "euler,~F", "EulerGamma,NumberSymbol,~A",
  "exp,ExpBase,~F", "Exp1,NumberSymbol,~A", "exp_polar,ExpBase,~F", "expint,~F",
  "ExprCondPair,Tuple,~B", "factorial,CombinatorialFunction,~F",
  "factorial2,CombinatorialFunction,~F", "FallingFactorial,CombinatorialFunction,~F",
  "fibonacci,~F", "FiniteSet,Set,Basic,Printable,EvalfMixin,object", "Float,Number,~A",
  "floor,RoundFunction,~F", "frac,~F", "fresnelc,FresnelIntegral,~F",
  "fresnels,FresnelIntegral,~F", "gamma,~F", "gegenbauer,OrthogonalPolynomial,~F", "genocchi,~F",
  "GoldenRatio,NumberSymbol,~A",
  "GreaterThan,_Greater,_Inequality,Relational,Boolean,Basic,Printable,EvalfMixin,object",
  "Half,RationalConstant,Rational,Number,~A", "hankel1,BesselBase,~F", "hankel2,BesselBase,~F",
  "harmonic,~F", "Heaviside,~F", "hermite,OrthogonalPolynomial,~F",
  "hermite_prob,OrthogonalPolynomial,~F",
  "hn1,SphericalHankelBase,SphericalBesselBase,BesselBase,~F",
  "hn2,SphericalHankelBase,SphericalBesselBase,BesselBase,~F", "hyper,TupleParametersBase,~F",
  "Idx,~E", "im,~F", "ImaginaryUnit,~A", "Implies,BooleanFunction,Application,Boolean,~B",
  "Indexed,~E", "IndexedBase,Expr,Basic,Printable,EvalfMixin,NotIterable,object",
  "Infinity,Number,~A", "Integer,Rational,Number,~A", "Integral,AddWithLimits,ExprWithLimits,~E",
  "Intersection,Set,LatticeOp,AssocOp,Basic,Printable,EvalfMixin,object",
  "Interval,Set,Basic,Printable,EvalfMixin,object", "Inverse,MatPow,MatrixExpr,~E",
  "ITE,BooleanFunction,Application,Boolean,~B", "jacobi,OrthogonalPolynomial,~F",
  "jacobi_symbol,~F", "jn,SphericalBesselBase,BesselBase,~F", "jtheta,~F", "kronecker_symbol,~F",
  "KroneckerDelta,~F", "laguerre,OrthogonalPolynomial,~F", "Lambda,~E", "LambertW,~F",
  "legendre,OrthogonalPolynomial,~F", "legendre_symbol,~F", "lerchphi,~F",
  "LessThan,_Less,_Inequality,Relational,Boolean,Basic,Printable,EvalfMixin,object", "li,~F",
  "Limit,~E", "log,~F", "loggamma,~F", "lowergamma,~F", "lucas,~F", "marcumq,~F",
  "MatAdd,MatrixExpr,Add,Expr,AssocOp,Basic,Printable,EvalfMixin,object",
  "mathieuc,MathieuBase,~F", "mathieucprime,MathieuBase,~F", "mathieus,MathieuBase,~F",
  "mathieusprime,MathieuBase,~F",
  "MatMul,MatrixExpr,Mul,Expr,AssocOp,Basic,Printable,EvalfMixin,object", "MatPow,MatrixExpr,~E",
  "MatrixSymbol,MatrixExpr,~E",
  "Max,MinMaxBase,Expr,LatticeOp,AssocOp,Application,Basic,Printable,EvalfMixin,object",
  "meijerg,TupleParametersBase,~F",
  "Min,MinMaxBase,Expr,LatticeOp,AssocOp,Application,Basic,Printable,EvalfMixin,object",
  "mobius,~F", "Mod,~F", "motzkin,~F", "Mul,Expr,AssocOp,Basic,Printable,EvalfMixin,object",
  "multigamma,~F", "NaN,Number,~A", "NegativeInfinity,Number,~A",
  "NegativeOne,IntegerConstant,Integer,Rational,Number,~A",
  "Not,BooleanFunction,Application,Boolean,~B", "One,IntegerConstant,Integer,Rational,Number,~A",
  "Or,LatticeOp,AssocOp,BooleanFunction,Application,Boolean,~B", "Order,~E", "owens_t,~F",
  "partition,~F", "periodic_argument,~F", "Pi,NumberSymbol,~A", "Piecewise,~F", "polar_lift,~F",
  "Poly,~B", "polygamma,~F", "polylog,~F", "Pow,~E", "primenu,~F", "primeomega,~F", "primepi,~F",
  "principal_branch,~F", "Product,ExprWithIntLimits,ExprWithLimits,~E", "Rational,Number,~A",
  "re,~F", "reduced_totient,~F", "Rem,~F", "riemann_xi,~F",
  "RisingFactorial,CombinatorialFunction,~F",
  "sec,ReciprocalTrigonometricFunction,TrigonometricFunction,~F",
  "sech,ReciprocalHyperbolicFunction,HyperbolicFunction,~F", "Shi,TrigonometricIntegral,~F",
  "Si,TrigonometricIntegral,~F", "sign,~F", "sin,TrigonometricFunction,~F", "sinc,~F",
  "SingularityFunction,~F", "sinh,HyperbolicFunction,~F", "stieltjes,~F",
  "StrictGreaterThan,_Greater,_Inequality,Relational,Boolean,Basic,Printable,EvalfMixin,object",
  "StrictLessThan,_Less,_Inequality,Relational,Boolean,Basic,Printable,EvalfMixin,object",
  "subfactorial,CombinatorialFunction,~F", "Subs,~E",
  "Sum,AddWithLimits,ExprWithIntLimits,ExprWithLimits,~E",
  "Symbol,AtomicExpr,Atom,Expr,Boolean,Basic,Printable,EvalfMixin,object",
  "tan,TrigonometricFunction,~F", "tanh,HyperbolicFunction,~F", "totient,~F", "Trace,~E",
  "Transpose,MatrixExpr,~E", "tribonacci,~F", "TribonacciConstant,NumberSymbol,~A", "Tuple,~B",
  "Unequality,Relational,Boolean,Basic,Printable,EvalfMixin,object",
  "Union,Set,LatticeOp,AssocOp,Basic,Printable,EvalfMixin,object", "uppergamma,~F",
  "Wild,Symbol,AtomicExpr,Atom,Expr,Boolean,Basic,Printable,EvalfMixin,object",
  "Xor,BooleanFunction,Application,Boolean,~B", "yn,SphericalBesselBase,BesselBase,~F", "Ynm,~F",
  "Zero,IntegerConstant,Integer,Rational,Number,~A", "zeta,~F")

.epath_mro_env <- new.env(parent = emptyenv())

.epath_mro_of <- function(cls) {
  if (!length(ls(.epath_mro_env))) {
    tails <- c("~F" = "DefinedFunction,Function,Application,Expr,Basic,Printable,EvalfMixin,object",
               "~A" = "AtomicExpr,Atom,Expr,Basic,Printable,EvalfMixin,object",
               "~E" = "Expr,Basic,Printable,EvalfMixin,object",
               "~B" = "Basic,Printable,object")
    for (s in .EPATH_MRO_TABLE) {
      for (k in names(tails)) s <- sub(k, tails[[k]], s, fixed = TRUE)
      parts <- strsplit(s, ",", fixed = TRUE)[[1L]]
      assign(parts[1L], parts, envir = .epath_mro_env)
    }
  }
  get0(cls, envir = .epath_mro_env, inherits = FALSE)
}

# The class names a value is an instance of, as [c.__name__ for c in
# type(x).mro()]; R lists stand for Python lists and tuples.
.epath_types <- function(x) {
  if (is_sym(x)) {
    cls <- switch(x$type,
      Function = x$name,
      Constant = switch(x$name, pi = "Pi", E = "Exp1", I = "ImaginaryUnit", oo = "Infinity",
                        "-oo" = "NegativeInfinity", zoo = "ComplexInfinity", nan = "NaN", x$name),
      BooleanAtom = if (identical(x$name, "true")) "BooleanTrue" else "BooleanFalse",
      Integer = if (x$value == 0) "Zero" else if (x$value == 1) "One" else if (x$value == -1) "NegativeOne" else "Integer",
      Rational = if (x$value == gmp::as.bigq(1, 2)) "Half" else "Rational",
      AccumBounds = "AccumulationBounds",
      x$type)
    m <- .epath_mro_of(cls)
    if (!is.null(m)) return(m)
    if (x$type == "Function")
      return(c(cls, "AppliedUndef", "Function", "Application", "Expr", "Basic", "Printable", "EvalfMixin", "object"))
    return(c(cls, "Expr", "Basic", "Printable", "EvalfMixin", "object"))
  }
  if (is.list(x)) return(c("list", "tuple", "object"))
  if (is.logical(x)) return(c("bool", "int", "object"))
  if (is.numeric(x)) return(c(if (all(x == round(x))) "int" else "float", "object"))
  if (is.character(x)) return(c("str", "object"))
  "object"
}

# hasattr(): every Basic has the is_* flags and its core methods; Tuple
# and Python containers are iterable.
# ponytail: a fixed attribute list, not every SymPy attribute; extend it
# when a path needs another one.
.epath_hasattr <- function(x, a) {
  if (is_sym(x)) {
    if (startsWith(a, "is_") || a %in% c("args", "func", "free_symbols", "atoms", "subs", "xreplace",
                                         "evalf", "n", "doit", "has", "expand", "simplify", "diff"))
      return(TRUE)
    return(x$type == "Tuple" && a %in% c("__iter__", "__len__", "__getitem__"))
  }
  is.list(x) && a %in% c("__iter__", "__len__", "__getitem__")
}

.epath_has <- function(x, attrs, types) {
  if (!length(attrs) && !length(types)) return(TRUE)
  if (length(attrs) && all(vapply(attrs, function(a) .epath_hasattr(x, a), TRUE))) return(TRUE)
  if (length(types) && any(.epath_types(x) %in% types)) return(TRUE)
  FALSE
}

# EPath._get_ordered_args: Add and Mul arguments in printing order
.epath_ordered_args <- function(e) {
  if (e$type == "Add") return(.as_ordered_terms(e))
  if (e$type == "Mul") return(.sympy_args(e))
  e$args
}

.epath_int <- function(s) {
  if (!grepl("^\\s*[+-]?[0-9]+(_[0-9]+)*\\s*$", s))
    stop(.sympy_condition("ValueError", sprintf("invalid literal for int() with base 10: '%s'", s)))
  as.integer(gsub("_", "", trimws(s), fixed = TRUE))
}

# slice(start, stop, step).indices(n), as 0-based positions
.epath_slice_positions <- function(sl, n) {
  step <- sl[[3L]] %||% 1L
  if (step == 0L) stop(.sympy_condition("ValueError", "slice step cannot be zero"))
  lower <- if (step > 0L) 0L else -1L
  upper <- if (step > 0L) n else n - 1L
  bound <- function(v, dflt) {
    if (is.null(v)) return(dflt)
    if (v < 0L) max(v + n, lower) else min(v, upper)
  }
  start <- bound(sl[[1L]], if (step > 0L) lower else upper)
  stop_ <- bound(sl[[2L]], if (step > 0L) upper else lower)
  out <- integer()
  i <- start
  while (if (step > 0L) i < stop_ else i > stop_) { out <- c(out, i); i <- i + step }
  out
}

# a Python index into a sequence of length n, as an R position or NA
.epath_pos <- function(i, n) {
  p <- if (i < 0L) n + i + 1L else i + 1L
  if (p < 1L || p > n) NA_integer_ else p
}

.epath_parse <- function(path) {
  if (!nzchar(path)) stop(.sympy_condition("ValueError", "empty EPath"))
  if (substr(path, 1L, 1L) != "/") stop(.sympy_condition("NotImplementedError", "non-root EPath"))
  path <- substring(path, 2L)
  out <- list()
  for (selector in strsplit(path, "/", fixed = TRUE)[[1L]][seq_len(max(1L, lengths(regmatches(path, gregexpr("/", path, fixed = TRUE))) + 1L))]) {
    selector <- trimws(if (is.na(selector)) "" else selector)
    if (!nzchar(selector)) stop(.sympy_condition("ValueError", "empty selector"))
    chars <- strsplit(selector, "")[[1L]]
    index <- 0L
    for (ch in chars) {
      if (grepl("^[[:alnum:]_|?]$", ch)) index <- index + 1L else break
    }
    attrs <- character(); types <- character()
    if (index) {
      elements <- substr(selector, 1L, index)
      selector <- substring(selector, index + 1L)
      for (el in strsplit(elements, "|", fixed = TRUE)[[1L]][seq_len(lengths(regmatches(elements, gregexpr("|", elements, fixed = TRUE))) + 1L)]) {
        el <- trimws(if (is.na(el)) "" else el)
        if (!nzchar(el)) stop(.sympy_condition("ValueError", "empty element"))
        if (endsWith(el, "?")) attrs <- c(attrs, substr(el, 1L, nchar(el) - 1L))
        else types <- c(types, el)
      }
    }
    span <- NULL
    if (selector != "*") {
      if (startsWith(selector, "[")) {
        i <- regexpr("]", selector, fixed = TRUE)
        if (i < 0L) stop(.sympy_condition("ValueError", "expected ']', got EOL"))
        sp <- substr(selector, 2L, i - 1L)
        if (!grepl(":", sp, fixed = TRUE)) {
          span <- .epath_int(sp)
        } else {
          pieces <- strsplit(sp, ":", fixed = TRUE)[[1L]]
          npieces <- lengths(regmatches(sp, gregexpr(":", sp, fixed = TRUE))) + 1L
          pieces <- c(pieces, rep("", npieces - length(pieces)))
          if (npieces > 3L) {
            pieces <- c(pieces[1:3], paste(pieces[-(1:3)], collapse = ":"))
            stop(.sympy_condition("TypeError", sprintf("slice expected at most 3 arguments, got %d", length(pieces))))
          }
          span <- lapply(pieces, function(p) if (!nzchar(p)) NULL else .epath_int(p))
          length(span) <- 3L
          class(span) <- "epath_slice"
        }
        selector <- substring(selector, i + 1L)
      }
      if (nzchar(selector)) stop(.sympy_condition("ValueError", "trailing characters in selector"))
    }
    out[[length(out) + 1L]] <- list(attrs = attrs, types = types, span = span)
  }
  out
}

.epath_select <- function(ep, expr) {
  result <- list()
  walk <- function(path, x) {
    if (!length(path)) { result[[length(result) + 1L]] <<- x; return(invisible()) }
    sel <- path[[1L]]; path <- path[-1L]
    if (is_sym(x)) args <- .epath_ordered_args(x)
    else if (is.list(x)) args <- x
    else return(invisible())
    sp <- sel$span
    if (inherits(sp, "epath_slice")) {
      args <- args[.epath_slice_positions(sp, length(args)) + 1L]
    } else if (!is.null(sp)) {
      p <- .epath_pos(sp, length(args))
      if (is.na(p)) return(invisible())
      args <- args[p]
    }
    for (a in args) if (.epath_has(a, sel$attrs, sel$types)) walk(path, a)
  }
  walk(ep, expr)
  result
}

.epath_apply <- function(ep, expr, func) {
  walk <- function(path, x) {
    if (!length(path)) return(func(x))
    sel <- path[[1L]]; path <- path[-1L]
    if (is_sym(x)) {
      if (!length(x$args)) return(x)
      args <- .epath_ordered_args(x); basic <- TRUE
    } else if (is.list(x)) {
      args <- x; basic <- FALSE
    } else return(x)
    n <- length(args)
    sp <- sel$span
    idx <- if (inherits(sp, "epath_slice")) .epath_slice_positions(sp, n) + 1L
           else if (!is.null(sp)) .epath_pos(sp, n)
           else seq_len(n)
    for (p in idx) {
      if (is.na(p)) next
      if (.epath_has(args[[p]], sel$attrs, sel$types)) args[p] <- list(walk(path, args[[p]]))
    }
    if (basic) sym_rebuild(x, args) else args
  }
  walk(ep, expr)
}

#' Manipulate expressions using paths
#'
#' `EPath(path)` compiles a path and `epath(path, expr, func)` selects or
#' modifies the parts of an expression (or a nested list of expressions)
#' that the path names, in the manner of XPath. A path is a sequence of
#' selectors, each starting with `/`:
#'
#' * `/*` selects all arguments.
#' * `/[0]`, `/[1:5]`, `/[1:5:2]` select by Python index or slice
#'   (0-based, negative indices count from the end).
#' * `/Symbol` or `/list|tuple` select by SymPy class name, including
#'   parent classes such as `TrigonometricFunction`, `Function` or `Atom`.
#' * `/__iter__?` selects by attribute.
#'
#' The arguments of a sum and a product are visited in printing order.
#'
#' With no `expr`, `epath()` returns the compiled path. With no `func`, it
#' returns the list of selected parts; otherwise it returns the expression
#' with `func` applied to each selected part. An `EPath` object has the
#' methods `$select(expr)` and `$apply(expr, func, args, kwargs)`.
#'
#' @param path A path string or an `EPath`.
#' @param expr An expression or a list of expressions.
#' @param func A function applied to the selected parts.
#' @param args,kwargs Extra positional and named arguments for `func`, as
#'   lists.
#' @return An `EPath`, a list of selected parts, or the modified expression.
#' @examples
#' x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z"); t <- Symbol("t")
#' expr <- list(list(list(x, 1), 2), list(list(3, y), z))
#' epath("/*/[0]/Symbol", expr)
#' epath("/*/[0]/Symbol", expr, function(e) e^2)
#' e2 <- t + sin(x + 1) + cos(x + y + sym_E)
#' epath("/*/*/Symbol", e2)
#' epath("/*/*/Symbol", e2, function(e) 2 * e)
#' @export
EPath <- function(path) {
  if (inherits(path, "EPath")) return(path)
  if (!is.character(path) || length(path) != 1L || is.na(path))
    stop(.sympy_condition("TypeError", "path must be a string"))
  ep <- .epath_parse(path)
  obj <- list(path = path, epath = ep)
  obj$select <- function(expr) .epath_select(ep, expr)
  obj$apply <- function(expr, func, args = NULL, kwargs = NULL)
    .epath_apply(ep, expr, function(e) do.call(func, c(list(e), args, kwargs)))
  structure(obj, class = "EPath")
}

#' @rdname EPath
#' @export
epath <- function(path, expr = NULL, func = NULL, args = NULL, kwargs = NULL) {
  p <- EPath(path)
  if (is.null(expr)) return(p)
  if (is.null(func)) return(p$select(expr))
  p$apply(expr, func, args, kwargs)
}

#' @export
print.EPath <- function(x, ...) {
  cat(sprintf("EPath('%s')\n", x$path))
  invisible(x)
}
