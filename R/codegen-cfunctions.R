# C99 math functions as unevaluated SymPy functions, after
# sympy/codegen/cfunctions.py. The code printers print them as the C library
# calls (expm1, log1p, exp2, log2, fma, log10, sqrt, cbrt, hypot, isnan,
# isinf). R already has expm1, log1p, log2 and log10, so those carry the
# sym_ prefix.

# Rewrite rules keyed by target and function name, read by rewrite() in
# R/piecewise-rewrite.R after its own table.
.cg_rewrites <- new.env()
.cg_add_rewrite <- function(target, name, fn) {
  cur <- .cg_rewrites[[target]]
  if (is.null(cur)) cur <- list()
  cur[[name]] <- fn
  assign(target, cur, envir = .cg_rewrites)
}

# exp.eval(arg) and log.eval(arg): the evaluated form, or NULL when the call
# stays unevaluated
.cg_fun_eval <- function(name, a) {
  r <- fcall(name, list(a))
  if (r$type == "Function" && r$name == name && length(r$args) == 1L && identical(r$args[[1L]]$key, a$key)) NULL else r
}

.cg_fact <- function(f, n) if (n %in% names(f)) unname(f[[n]]) else NA

.cg_expm1 <- function(x) .sub(fcall("exp", list(x)), S_One)
.cg_log1p <- function(x) fcall("log", list(.add(list(x, S_One))))
.cg_exp2 <- function(x) .pow(S_Two, x)
.cg_log2 <- function(x) .div(fcall("log", list(x)), fcall("log", list(S_Two)))
.cg_log10 <- function(x) .div(fcall("log", list(x)), fcall("log", list(.int(10L))))
.cg_hypot <- function(x, y) .pow(.add(list(.pow(x, S_Two), .pow(y, S_Two))), S_Half)

def_function("expm1",
  eval = function(args) {
    a <- args[[1L]]
    # a nonzero Float is evaluated at its precision (mpmath's expm1)
    if (a$type == "Float" && a$value != 0) {
      p <- .float_prec(a); x <- .float_mp(a)
      w <- p + 20 + max(0, -.mp_mag(x))
      r <- .mp_sub(.mp_exp(x, w), .mpf(1, 0), w)
      return(.mp_float(.mp_round(r$m, r$e, p), p))
    }
    r <- .cg_fun_eval("exp", a)
    if (!is.null(r)) .sub(r, S_One) else NULL
  },
  diff = function(args, i) fcall("exp", args),
  numeric = function(x) if (is.complex(x)) exp(x) - 1 else expm1(x),
  facts = function(args, afs) c(real = .cg_fact(afs[[1L]], "real"), finite = .cg_fact(afs[[1L]], "finite")))
.expand_hooks$expm1 <- function(e) .cg_expm1(e$args[[1L]])
.cg_add_rewrite("exp", "expm1", function(a) .cg_expm1(a))
.cg_add_rewrite("tractable", "expm1", function(a) .cg_expm1(a))

def_function("log1p",
  eval = function(args) {
    a <- args[[1L]]
    if (a$type %in% c("Integer", "Rational")) return(fcall("log", list(.add(list(a, S_One)))))
    if (a$type != "Float") return(.cg_fun_eval("log", .add(list(a, S_One))))
    fcall("log", list(.add(list(make_rational(gmp::as.bigq(a$value)), S_One))))
  },
  diff = function(args, i) .div(S_One, .add(list(args[[1L]], S_One))),
  numeric = function(x) if (is.complex(x)) log(x + 1) else log1p(x),
  facts = function(args, afs) {
    a <- args[[1L]]
    f1 <- .facts(.add(list(a, S_One)))
    fin <- if (isTRUE(.cg_fact(f1, "zero"))) FALSE else .cg_fact(afs[[1L]], "finite")
    c(real = .cg_fact(f1, "nonnegative"), finite = fin, positive = .cg_fact(afs[[1L]], "positive"),
      zero = .cg_fact(afs[[1L]], "zero"), nonnegative = .cg_fact(afs[[1L]], "nonnegative"))
  })
.expand_hooks$log1p <- function(e) .cg_log1p(e$args[[1L]])
.cg_add_rewrite("log", "log1p", function(a) .cg_log1p(a))
.cg_add_rewrite("tractable", "log1p", function(a) .cg_log1p(a))

def_function("exp2",
  eval = function(args) if (.expr_is_number(args[[1L]])) .cg_exp2(args[[1L]]) else NULL,
  diff = function(args, i) .mul(list(fcall("exp2", args), fcall("log", list(S_Two)))),
  numeric = function(x) 2^x)
.expand_hooks$exp2 <- function(e) .cg_exp2(e$args[[1L]])
.cg_add_rewrite("Pow", "exp2", function(a) .cg_exp2(a))
.cg_add_rewrite("tractable", "exp2", function(a) .cg_exp2(a))

# log.eval(arg, base): the evaluated log(arg)/log(base) when it is an atom
.cg_log_base_eval <- function(a, base) {
  r <- tryCatch(fcall("log", list(a, base)), error = function(e) NULL)
  if (!is.null(r) && !length(r$args)) r else NULL
}

def_function("log2",
  eval = function(args) {
    a <- args[[1L]]
    if (.expr_is_number(a)) return(.cg_log_base_eval(a, S_Two))
    if (a$type == "Pow" && identical(a$args[[1L]]$key, S_Two$key)) return(a$args[[2L]])
    NULL
  },
  diff = function(args, i) .div(S_One, .mul(list(fcall("log", list(S_Two)), args[[1L]]))),
  numeric = function(x) log2(x))
.expand_hooks$log2 <- function(e) .cg_log2(e$args[[1L]])
.cg_add_rewrite("log", "log2", function(a) .cg_log2(a))
.cg_add_rewrite("tractable", "log2", function(a) .cg_log2(a))

def_function("fma", nargs = 3L,
  diff = function(args, i) if (i %in% c(1L, 2L)) args[[3L - i]] else S_One,
  numeric = function(x, y, z) x * y + z)
.expand_hooks$fma <- function(e) .add(list(.mul(list(e$args[[1L]], e$args[[2L]])), e$args[[3L]]))

def_function("log10",
  eval = function(args) {
    a <- args[[1L]]
    if (.expr_is_number(a)) return(.cg_log_base_eval(a, .int(10L)))
    if (a$type == "Pow" && identical(a$args[[1L]]$key, .int(10L)$key)) return(a$args[[2L]])
    NULL
  },
  diff = function(args, i) .div(S_One, .mul(list(fcall("log", list(.int(10L))), args[[1L]]))),
  numeric = function(x) log10(x))
.expand_hooks$log10 <- function(e) .cg_log10(e$args[[1L]])
.cg_add_rewrite("log", "log10", function(a) .cg_log10(a))
.cg_add_rewrite("tractable", "log10", function(a) .cg_log10(a))

def_function("Sqrt",
  diff = function(args, i) .div(.pow(args[[1L]], Rational(-1, 2)), S_Two))
.expand_hooks$Sqrt <- function(e) .pow(e$args[[1L]], S_Half)
.cg_add_rewrite("Pow", "Sqrt", function(a) .pow(a, S_Half))
.cg_add_rewrite("tractable", "Sqrt", function(a) .pow(a, S_Half))

def_function("Cbrt",
  diff = function(args, i) .div(.pow(args[[1L]], Rational(-2, 3)), .int(3L)))
.expand_hooks$Cbrt <- function(e) .pow(e$args[[1L]], Rational(1, 3))
.cg_add_rewrite("Pow", "Cbrt", function(a) .pow(a, Rational(1, 3)))
.cg_add_rewrite("tractable", "Cbrt", function(a) .pow(a, Rational(1, 3)))

def_function("hypot", nargs = 2L,
  diff = function(args, i) .div(.mul(list(S_Two, args[[i]])), .mul(list(S_Two, fcall("hypot", args)))),
  numeric = function(x, y) sqrt(x^2 + y^2))
.expand_hooks$hypot <- function(e) .cg_hypot(e$args[[1L]], e$args[[2L]])
.cg_add_rewrite("Pow", "hypot", function(a, b) .cg_hypot(a, b))
.cg_add_rewrite("tractable", "hypot", function(a, b) .cg_hypot(a, b))

def_function("isnan",
  eval = function(args) {
    a <- args[[1L]]
    if (is_nan(a)) return(S_true)
    if (.expr_is_number(a)) return(S_false)
    NULL
  })
def_function("isinf",
  eval = function(args) {
    f <- .facts(args[[1L]])
    if (isTRUE(f[["infinite"]])) return(S_true)
    if (isTRUE(f[["finite"]])) return(S_false)
    NULL
  })

#' C99 math functions
#'
#' Unevaluated functions for the C99 math library, after SymPy's
#' `sympy.codegen.cfunctions`: `sym_expm1(x)` is `exp(x) - 1`, `sym_log1p(x)`
#' is `log(x + 1)`, `exp2(x)` is `2^x`, `sym_log2(x)` and `sym_log10(x)` are
#' the base 2 and base 10 logarithms, `fma(x, y, z)` is `x*y + z`, `Sqrt(x)`
#' and `Cbrt(x)` are the square and cube roots, `hypot(x, y)` is
#' `sqrt(x^2 + y^2)`, and `isnan(x)` and `isinf(x)` test for a NaN and an
#' infinity. They differentiate, and `expand(e, func = TRUE)` or
#' [rewrite()] turns them into elementary functions; the code printers print
#' them as the matching C calls.
#'
#' @param x,y,z Expressions.
#' @return A `sym` expression.
#' @examples
#' x <- Symbol("x")
#' sym_expm1(x)
#' expand(sym_log1p(x), func = TRUE)
#' diff(hypot(x, 2), x)
#' @name cfunctions
NULL

#' @rdname cfunctions
#' @export
sym_expm1 <- function(x) fcall("expm1", list(as_sym(x)))
#' @rdname cfunctions
#' @export
sym_log1p <- function(x) fcall("log1p", list(as_sym(x)))
#' @rdname cfunctions
#' @export
exp2 <- function(x) fcall("exp2", list(as_sym(x)))
#' @rdname cfunctions
#' @export
sym_log2 <- function(x) fcall("log2", list(as_sym(x)))
#' @rdname cfunctions
#' @export
fma <- function(x, y, z) fcall("fma", list(as_sym(x), as_sym(y), as_sym(z)))
#' @rdname cfunctions
#' @export
sym_log10 <- function(x) fcall("log10", list(as_sym(x)))
#' @rdname cfunctions
#' @export
Sqrt <- function(x) fcall("Sqrt", list(as_sym(x)))
#' @rdname cfunctions
#' @export
Cbrt <- function(x) fcall("Cbrt", list(as_sym(x)))
#' @rdname cfunctions
#' @export
hypot <- function(x, y) fcall("hypot", list(as_sym(x), as_sym(y)))
#' @rdname cfunctions
#' @export
isnan <- function(x) fcall("isnan", list(as_sym(x)))
#' @rdname cfunctions
#' @export
isinf <- function(x) fcall("isinf", list(as_sym(x)))
