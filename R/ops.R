# Operator overloading and R group generics for sym objects.

#' @export
chooseOpsMethod.sym <- function(x, y, mx, my, cl, reverse)
  !inherits(y, c("sym_ordinal", "sym_omega_power", "DifferentialOperator", "HolonomicFunction", "RecurrenceOperator", "physics_obj", "sym_ring_element", "pyobj"))


#' Arithmetic, comparison and logic operators for symbolic expressions
#'
#' `+`, `-`, `*`, `/`, `^`, `%%` and `%/%` build canonical symbolic
#' expressions. `==` and `!=` test structural equality, as in SymPy; use
#' [Eq()] for a symbolic equation. `<`, `<=`, `>` and `>=` return `TRUE` or
#' `FALSE` when the answer is decidable and an unevaluated relational
#' otherwise. `&`, `|` and `!` build logical expressions.
#'
#' @param e1,e2 Operands; at least one is a `sym`.
#' @return A `sym` expression or a logical value.
#' @name sym-operators
#' @examples
#' x <- Symbol("x")
#' (x + 1)^2 - 1
#' x + x == 2 * x
#' Symbol("p", positive = TRUE) > 0
#' @export
Ops.sym <- function(e1, e2) {
  if (missing(e2)) {
    a <- as_sym(e1)
    if (.Generic == "-" && .is_seq(a)) return(.seq_coeff_mul(a, S_NegativeOne))
    if (.Generic == "-" && .is_fourier(a)) return(scale(a, -1))
    if (.Generic == "-" && .is_matexpr(a)) return(doit(MatMul(S_NegativeOne, a)))
    if (.Generic == "-" && .fps_is(a)) return(.fps_neg(a))
    # Dimension.__neg__ returns the dimension
    if (.Generic == "-" && a$type == "Dimension") return(a)
    return(switch(.Generic,
                  "-" = .neg(a),
                  "+" = a,
                  "!" = Not(a),
                  stop("unsupported unary operator ", .Generic)))
  }
  h <- .ops_hooks$dispatch
  if (!is.null(h)) {
    r <- h(.Generic, e1, e2)
    if (!is.null(r)) return(r)
  }
  a <- as_sym(e1)
  b <- as_sym(e2)
  # AccumBounds has the higher _op_priority: its __add__, __sub__, __mul__
  # and __truediv__ methods handle the operator
  if (.Generic %in% c("+", "-", "*", "/") && (.is_ab(a) || .is_ab(b))) {
    r <- switch(.Generic,
      "+" = if (.is_ab(a)) .ab_add(a, b) else .ab_add(b, a),
      "-" = if (.is_ab(a)) .ab_add(a, .neg(b)) else .ab_add(.neg(b), a),
      "*" = if (.is_ab(a)) .ab_mul(a, b) else .ab_mul(b, a),
      "/" = .ab_div(a, b))
    if (!is.null(r)) return(r)
  }
  # Set.__mul__: the Cartesian product
  if (.Generic == "*" && is_set(a) && is_set(b)) return(ProductSet(a, b))
  switch(.Generic,
         "+" = .add(list(a, b)),
         "-" = .sub(a, b),
         "*" = .mul(list(a, b)),
         "/" = .div(a, b),
         "^" = .pow(a, b),
         "%%" = fcall("Mod", list(a, b)),
         "%/%" = fcall("floor", list(.div(a, b))),
         "==" = identical(a$key, b$key),
         "!=" = !identical(a$key, b$key),
         "<" = .rel_to_r(Lt(a, b)),
         "<=" = .rel_to_r(Le(a, b)),
         ">" = .rel_to_r(Gt(a, b)),
         ">=" = .rel_to_r(Ge(a, b)),
         "&" = And(a, b),
         "|" = Or(a, b),
         stop("unsupported operator ", .Generic))
}

.rel_to_r <- function(r) {
  if (r$type == "BooleanAtom") return(r$name == "true")
  r
}

#' @rdname sym-operators
#' @param x A `sym` expression.
#' @param ... Further arguments, such as `base` for `log()`.
#' @export
Math.sym <- function(x, ...) {
  x <- as_sym(x)
  switch(.Generic,
         abs = fcall("Abs", list(x)),
         sign = fcall("sign", list(x)),
         sqrt = if (.is_matexpr(x) && !.is_explicit(x)) MatPow(x, S_Half) else .sqrt(x),
         floor = fcall("floor", list(x)),
         ceiling = fcall("ceiling", list(x)),
         trunc = .mul(list(fcall("sign", list(x)), fcall("floor", list(fcall("Abs", list(x)))))),
         exp = fcall("exp", list(x)),
         log = {
           dots <- list(...)
           if (length(dots)) fcall("log", list(x, as_sym(dots[[1L]]))) else fcall("log", list(x))
         },
         expm1 = .sub(fcall("exp", list(x)), S_One),
         log1p = fcall("log", list(.add(list(S_One, x)))),
         log2 = fcall("log", list(x, S_Two)),
         log10 = fcall("log", list(x, .int(10))),
         cos = , sin = , tan = , acos = , asin = , atan = , cosh = , sinh = ,
         tanh = , acosh = , asinh = , atanh = fcall(.Generic, list(x)),
         cospi = fcall("cos", list(.mul(list(S_Pi, x)))),
         sinpi = fcall("sin", list(.mul(list(S_Pi, x)))),
         tanpi = fcall("tan", list(.mul(list(S_Pi, x)))),
         gamma = fcall("gamma", list(x)),
         lgamma = fcall("loggamma", list(x)),
         digamma = fcall("polygamma", list(S_Zero, x)),
         trigamma = fcall("polygamma", list(S_One, x)),
         round = , signif = {
           v <- N(x)
           if (!is_number_type(v)) stop(.Generic, "() needs a number")
           Float(get(.Generic, baseenv())(v$value, ...))
         },
         stop(.Generic, "() is not supported for sym objects"))
}

#' @rdname sym-operators
#' @param na.rm Ignored.
#' @export
Summary.sym <- function(..., na.rm = FALSE) {
  args <- lapply(list(...), as_sym)
  switch(.Generic,
         sum = .add(args),
         prod = .mul(args),
         max = fcall("Max", args),
         min = fcall("Min", args),
         range = list(fcall("Min", args), fcall("Max", args)),
         stop(.Generic, "() is not supported for sym objects"))
}

#' @rdname sym-operators
#' @param z A `sym` expression.
#' @export
Complex.sym <- function(z) {
  z <- as_sym(z)
  switch(.Generic,
         Re = fcall("re", list(z)),
         Im = fcall("im", list(z)),
         Mod = fcall("Abs", list(z)),
         Arg = fcall("arg", list(z)),
         Conj = fcall("conjugate", list(z)))
}

#' @export
c.sym <- function(...) lapply(list(...), as_sym)

# ---------------------------------------------------------------------------
# Relationals

# A Float infinity acts as oo or -oo, and a finite Float next to an exact
# number compares through its exact binary value, as SymPy Float does.
.rel_float <- function(z, other) {
  if (z$type != "Float" || is.nan(z$value)) return(z)
  if (is.infinite(z$value)) return(if (z$value > 0) S_Infinity else S_NegativeInfinity)
  if (is_exact_number(other)) return(make_rational(gmp::as.bigq(z$value)))
  z
}

.rel_eval <- function(type, a, b) {
  # Eq/Ne of sets: is_eq with Interval._eval_Eq and the _eval_is_eq handlers
  if ((type == "Equality" || type == "Unequality") && (is_set(a) || is_set(b))) {
    r <- .set_is_eq(a, b)
    if (is.null(r)) return(NULL)
    if (type == "Unequality") return(if (.is_true(r)) S_false else if (.is_false(r)) S_true else NULL)
    return(r)
  }
  # _eval_is_eq(MatrixExpr, MatrixExpr)
  if ((type == "Equality" || type == "Unequality") && .is_matexpr(a) && .is_matexpr(b)) {
    r <- .me_is_eq(a, b)
    if (is.null(r)) return(NULL)
    return(if (xor(r, type == "Unequality")) S_true else S_false)
  }
  a0 <- a
  a <- .rel_float(a, b)
  b <- .rel_float(b, a0)
  if (type == "Equality" || type == "Unequality") {
    nanish <- function(z) is_nan(z) || (z$type == "Float" && is.nan(z$value))
    # is_eq forms lhs - rhs only when both sides are Expr: two distinct
    # BooleanAtoms differ, a Boolean never equals a non-Symbol Expr, and any
    # other pair involving a Boolean stays undecided.
    boolish <- function(z) z$type %in% c(.BOOL_TYPES, .REL_TYPES, "Nand", "Nor", "Exclusive", "Contains")
    eq <- if (nanish(a) || nanish(b)) FALSE else if (identical(a$key, b$key)) TRUE
    else if (boolish(a) || boolish(b)) {
      if (a$type == "BooleanAtom" && b$type == "BooleanAtom") FALSE
      else if (a$type != "Symbol" && b$type != "Symbol" && xor(boolish(a), boolish(b))) FALSE
      else NA
    } else {
      # is_eq: an infinite side is decided by sign, reality or arg
      ie <- .coremm_is_eq_inf(a, b)
      if (!is.null(ie)) {
        if (is.na(ie)) return(NULL)
        return(if (xor(ie, type == "Unequality")) S_true else S_false)
      }
      # is_eq: no difference that would cancel a term not known to be finite
      if (.cancels_nonfinite(a, b)) return(NULL)
      d <- .sub(a, b)
      if (is_zero_num(d)) TRUE
      else if (is_number_type(d)) (if (d$type == "Float") d$value == 0 else FALSE)
      else if (is_inf_const(d) || is_nan(d)) FALSE
      else {
        fd <- .facts(d)
        if (isTRUE(fd[["zero"]])) TRUE
        # is_eq: dif.is_zero False decides, as for a nonreal difference
        else if (isTRUE(fd[["nonzero"]]) || (identical(fd[["zero"]], FALSE) && !identical(unname(fd["commutative"]), FALSE))) FALSE
        else if (a$type == "BooleanAtom" && b$type == "BooleanAtom") FALSE
        # Add._eval_is_extended_real through _fuzzy_group: extended real
        # terms and one that is not make a difference that is not real,
        # so not zero
        else if (d$type == "Add" && {
          er <- vapply(d$args, function(t) { v <- ask_fact(t, "extended_real"); if (length(v) == 1L) v else NA }, NA)
          !anyNA(er) && sum(!er) == 1L
        }) FALSE
        # a number difference whose evalf is clearly nonzero: the evalf sign
        # makes is_extended_real or is_zero False, as in _n2
        # ponytail: double precision with a fixed 1e-8 cutoff; use mpmath
        # evalf with a precision check if near cancellations show up
        else if (!length(free_symbols(d)) && {
          v <- tryCatch(.neval(d, list()), error = function(e) NA)
          length(v) == 1L && !is.na(v) && is.finite(Mod(v)) && Mod(v) > 1e-8
        }) FALSE
        else NA
      }
    }
    if (is.na(eq)) return(NULL)
    if (type == "Unequality") eq <- !eq
    return(if (eq) S_true else S_false)
  }
  if (a$type == "AccumulationBounds" || b$type == "AccumulationBounds") return(.ab_rel(type, a, b))
  if (is_const(a, "I") || is_const(b, "I")) stop("Invalid comparison of non-real I")
  if (is_inf_const(a) && identical(a$key, b$key) && !is_const(a, "zoo"))
    return(if (type %in% c("LessThan", "GreaterThan")) S_true else S_false)
  if (is_nan(a) || is_nan(b)) stop("Invalid NaN comparison")
  # SymPy decides a symbolic inequality from the difference only when both
  # sides are extended real.
  ext <- isTRUE(.facts(a)[["extended_real"]]) && isTRUE(.facts(b)[["extended_real"]])
  if ((length(free_symbols(a)) || length(free_symbols(b))) && !ext) return(NULL)
  d <- .sub(a, b)
  if (is_nan(d)) stop("Invalid NaN comparison")
  fd <- .facts(d)
  if (is_const(d, "oo")) fd[c("positive", "negative", "zero", "nonnegative", "nonpositive")] <-
    c(TRUE, FALSE, FALSE, TRUE, FALSE)
  if (is_const(d, "-oo")) fd[c("positive", "negative", "zero", "nonnegative", "nonpositive")] <-
    c(FALSE, TRUE, FALSE, FALSE, TRUE)
  if (identical(fd[["real"]], FALSE) && identical(fd[["extended_real"]], FALSE) &&
      length(free_symbols(d)) == 0L)
    stop("Invalid comparison of non-real ", .str(d))
  # anything extended real lies in [-oo, oo] (SymPy is_ge)
  ge_inf <- function(l, r) is_const(l, "oo") || is_const(r, "-oo")
  r <- if (!ext) NULL else switch(type,
              GreaterThan = if (ge_inf(a, b)) TRUE, LessThan = if (ge_inf(b, a)) TRUE,
              StrictLessThan = if (ge_inf(a, b)) FALSE, StrictGreaterThan = if (ge_inf(b, a)) FALSE)
  if (!is.null(r)) return(if (r) S_true else S_false)
  # is_ge reads diff.is_extended_nonnegative (the sign facts of an
  # extended real difference)
  sg <- function(n) {
    v <- fd[paste0("extended_", n)][[1L]]
    isTRUE(v) || (is.na(v %||% NA) && isTRUE(fd[[n]]))
  }
  r <- switch(type,
              StrictLessThan = if (sg("negative")) TRUE else if (sg("nonnegative")) FALSE else NA,
              LessThan = if (sg("nonpositive")) TRUE else if (sg("positive")) FALSE else NA,
              StrictGreaterThan = if (sg("positive")) TRUE else if (sg("nonpositive")) FALSE else NA,
              GreaterThan = if (sg("nonnegative")) TRUE else if (sg("negative")) FALSE else NA)
  if (is.na(r) && !length(free_symbols(d)) && isTRUE(fd[["real"]]) && !is_number_type(d) && !is_inf_const(d)) {
    # a comparable real constant: its sign from certified double bounds
    bd <- tryCatch(.nt_bounds(d), error = function(err) NULL)
    if (!is.null(bd)) {
      sg <- if (bd[[1L]] > 0) 1L else if (bd[[2L]] < 0) -1L else 0L
      if (sg != 0L) r <- switch(type, StrictLessThan = , LessThan = sg < 0, StrictGreaterThan = , GreaterThan = sg > 0)
    }
  }
  if (is.na(r)) return(NULL)
  if (r) S_true else S_false
}

# is_eq's guard: lhs - rhs would remove a common term (same coefficient,
# from as_coefficients_dict) whose finiteness is unknown, as in x + 1 vs x.
.cancels_nonfinite <- function(a, b) {
  fin <- function(z) isTRUE(.facts(z)[["finite"]])
  if (fin(a) && fin(b)) return(FALSE)
  cd <- function(e) {
    out <- list()
    for (t in if (e$type == "Add") e$args else list(e)) {
      sp <- as_coeff_mul_split(t)
      k <- sp[[2L]]$key
      out[[k]] <- list(sp[[2L]], if (is.null(out[[k]])) sp[[1L]] else num_add(out[[k]][[2L]], sp[[1L]]))
    }
    out
  }
  L <- cd(a); R <- cd(b)
  for (k in intersect(names(L), names(R))) {
    if (k == S_One$key || !identical(L[[k]][[2L]]$key, R[[k]][[2L]]$key)) next
    if (!fin(L[[k]][[1L]])) return(TRUE)
  }
  FALSE
}

.make_rel <- function(type, a, b, evaluate = TRUE) {
  a <- if (inherits(a, "symmatrix")) .as_block(a) else as_sym(a)
  b <- if (inherits(b, "symmatrix")) .as_block(b) else as_sym(b)
  if (evaluate) {
    r <- .rel_eval(type, a, b)
    if (!is.null(r)) return(r)
  }
  new_sym(type, list(a, b))
}

for (.t in .REL_TYPES) register_rebuilder(.t, local({
  tt <- .t
  function(e, args) .make_rel(tt, args[[1L]], args[[2L]])
}))
rm(.t)

#' Symbolic equations and inequalities
#'
#' `Eq(a, b)` is the equation `a = b`, `Ne()` is `a != b`, and `Lt()`, `Le()`,
#' `Gt()`, `Ge()` are the strict and non strict inequalities. When the truth
#' value is decidable they return `S_true` or `S_false`. Use `lhs()` and
#' `rhs()` to take them apart.
#'
#' @param a,b Expressions.
#' @param evaluate If `FALSE`, keep the relation unevaluated.
#' @param e A relational expression.
#' @return A relational `sym`, or a boolean atom.
#' @examples
#' x <- Symbol("x")
#' Eq(x^2, 4)
#' Lt(1, 2)
#' @export
Eq <- function(a, b = 0, evaluate = TRUE) .make_rel("Equality", a, b, evaluate)
#' @rdname Eq
#' @export
Ne <- function(a, b, evaluate = TRUE) .make_rel("Unequality", a, b, evaluate)
#' @rdname Eq
#' @export
Lt <- function(a, b, evaluate = TRUE) .make_rel("StrictLessThan", a, b, evaluate)
#' @rdname Eq
#' @export
Le <- function(a, b, evaluate = TRUE) .make_rel("LessThan", a, b, evaluate)
#' @rdname Eq
#' @export
Gt <- function(a, b, evaluate = TRUE) .make_rel("StrictGreaterThan", a, b, evaluate)
#' @rdname Eq
#' @export
Ge <- function(a, b, evaluate = TRUE) .make_rel("GreaterThan", a, b, evaluate)
#' @rdname Eq
#' @export
lhs <- function(e) as_sym(e)$args[[1L]]
#' @rdname Eq
#' @export
rhs <- function(e) as_sym(e)$args[[2L]]

is_relational <- function(e) e$type %in% .REL_TYPES

# ---------------------------------------------------------------------------
# Wrappers for functions that are not R group generics.

.w1 <- function(name) {
  force(name)
  function(x) fcall(name, list(as_sym(x)))
}

#' Elementary and complex functions
#'
#' Symbolic versions of trigonometric, hyperbolic, exponential, complex and
#' piecewise elementary functions. The R generics `sin()`, `cos()`, `exp()`,
#' `log()`, `sqrt()`, `abs()`, `Re()`, `Im()`, `sym_Mod()`, `Arg()`, `Conj()`,
#' `floor()`, `ceiling()`, `max()` and `min()` also work on `sym` objects.
#'
#' @param x,y Expressions.
#' @param n Root index.
#' @param k Branch index or derivative order.
#' @param h0 Value of the Heaviside function at zero.
#' @param delta_range For `KroneckerDelta()`, an optional range `(lower,
#'   upper)` of the indices; an index outside it gives zero.
#' @param ... Expressions.
#' @return A `sym` expression.
#' @name elementary-functions
#' @examples
#' x <- Symbol("x")
#' cot(sym_pi / 4)
#' Abs(-3 * x)
#' Max(1, 2, x)
NULL

#' @rdname elementary-functions
#' @export
cot <- .w1("cot")
#' @rdname elementary-functions
#' @export
sec <- .w1("sec")
#' @rdname elementary-functions
#' @export
csc <- .w1("csc")
#' @rdname elementary-functions
#' @export
acot <- .w1("acot")
#' @rdname elementary-functions
#' @export
asec <- .w1("asec")
#' @rdname elementary-functions
#' @export
acsc <- .w1("acsc")
#' @rdname elementary-functions
#' @export
coth <- .w1("coth")
#' @rdname elementary-functions
#' @export
sech <- .w1("sech")
#' @rdname elementary-functions
#' @export
csch <- .w1("csch")
#' @rdname elementary-functions
#' @export
acoth <- .w1("acoth")
#' @rdname elementary-functions
#' @export
asech <- .w1("asech")
#' @rdname elementary-functions
#' @export
acsch <- .w1("acsch")
#' @rdname elementary-functions
#' @export
sinc <- .w1("sinc")
#' @rdname elementary-functions
#' @export
Abs <- .w1("Abs")
#' @rdname elementary-functions
#' @export
re <- function(x) if (inherits(x, "symmatrix") || .is_matexpr(x)) .me_real_imag(x)[[1L]] else fcall("re", list(as_sym(x)))
#' @rdname elementary-functions
#' @export
im <- function(x) if (inherits(x, "symmatrix") || .is_matexpr(x)) .me_real_imag(x)[[2L]] else fcall("im", list(as_sym(x)))
#' @rdname elementary-functions
#' @export
conjugate <- function(x) if (inherits(x, "symmatrix") || .is_matexpr(x)) .me_conjugate(x) else fcall("conjugate", list(as_sym(x)))
#' @rdname elementary-functions
#' @export
arg <- .w1("arg")
#' @rdname elementary-functions
#' @export
frac <- .w1("frac")
#' @rdname elementary-functions
#' @export
Heaviside <- function(x, h0 = NULL) {
  if (is.null(h0)) fcall("Heaviside", list(as_sym(x)))
  else fcall("Heaviside", list(as_sym(x), as_sym(h0)))
}
#' @rdname elementary-functions
#' @export
DiracDelta <- function(x, k = 0) {
  if (identical(k, 0)) fcall("DiracDelta", list(as_sym(x)))
  else fcall("DiracDelta", list(as_sym(x), as_sym(k)))
}
#' @rdname elementary-functions
#' @export
KroneckerDelta <- function(x, y, delta_range = NULL) {
  if (is.null(delta_range)) return(fcall("KroneckerDelta", list(as_sym(x), as_sym(y))))
  if (!is_sym(delta_range)) delta_range <- do.call(Tuple, lapply(delta_range, as_sym))
  fcall("KroneckerDelta", list(as_sym(x), as_sym(y), delta_range))
}
#' @rdname elementary-functions
#' @export
LambertW <- function(x, k = 0) {
  if (identical(k, 0)) fcall("LambertW", list(as_sym(x)))
  else fcall("LambertW", list(as_sym(x), as_sym(k)))
}
#' @rdname elementary-functions
#' @export
Max <- function(...) fcall("Max", lapply(list(...), as_sym))
#' @rdname elementary-functions
#' @export
Min <- function(...) fcall("Min", lapply(list(...), as_sym))
#' @noRd
atan2 <- function(y, x) {
  if (is_sym(y) || is_sym(x)) return(fcall("atan2", list(as_sym(y), as_sym(x))))
  base::atan2(y, x)
}
#' @rdname elementary-functions
#' @export
root <- function(x, n, k = 0) {
  r <- .pow(as_sym(x), .inv(as_sym(n)))
  if (!identical(k, 0)) r <- .mul(list(r, fcall("exp", list(.div(.mul(list(S_Two, S_I, S_Pi, as_sym(k))), as_sym(n))))))
  r
}
#' @rdname elementary-functions
#' @export
cbrt <- function(x) .pow(as_sym(x), .rat(1, 3))
#' @rdname elementary-functions
#' @export
real_root <- function(x, n) {
  x <- as_sym(x)
  n <- as_sym(n)
  if (is_number_type(x) && num_sign(x) < 0 && n$type == "Integer" && .int_mod(n, 2L) == 1L)
    return(.neg(.pow(.neg(x), .inv(n))))
  .pow(x, .inv(n))
}
#' @noRd
Mod <- function(x, y, ...) {
  if (missing(x) || missing(y) || ...length()) return(.base_call(quote(base::Mod)))
  if (is_sym(x) || is_sym(y)) return(fcall("Mod", list(as_sym(x), as_sym(y))))
  x - y * floor(x / y)
}

def_function("Mod", nargs = 2L,
             eval = function(args) .mod_eval(args[[1L]], args[[2L]]),
             numeric = function(a, b) Re(a) - Re(b) * floor(Re(a) / Re(b)),
             str = function(e) paste0("Mod(", .str(e$args[[1L]]), ", ", .str(e$args[[2L]]), ")"),
             facts = function(args, afs) {
               out <- c(commutative = TRUE)
               if (.all_true(afs, "integer")) out["integer"] <- TRUE
               if (isTRUE(afs[[2L]][["positive"]]) && .all_true(afs, "real")) out["nonnegative"] <- TRUE
               out
             })

.is_mod <- function(e) e$type == "Function" && identical(e$name, "Mod")
.fz <- function(type, a, b) {
  r <- tryCatch(.rel_eval(type, a, b), error = function(e) NULL)
  if (is.null(r)) NA else .is_true(r)
}
.equal_valued_one <- function(e) is_number_type(e) && isTRUE(num_value_dbl(e) == 1)
.count_mod <- function(e) {
  n <- if (.is_mod(e)) 1L else 0L
  for (a in e$args) if (is_sym(a)) n <- n + .count_mod(a)
  n
}

# Python's int() of a numeric expression (Expr.__int__), or NULL where it
# raises TypeError.
.py_int_expr <- function(r) {
  if (length(free_symbols(r))) return(NULL)
  if (is_exact_number(r)) return(Integer(gmp::as.bigz(as.character(.trunc_bigq(r$value)))))
  v <- tryCatch(N(r, 40), error = function(e) NULL)
  if (is.null(v) || v$type != "Float" || !is.finite(num_value_dbl(v))) return(NULL)
  q <- .mp_to_q(.float_mp(v))
  r2 <- gmp::as.bigq(round(q * 100)) / 100
  i <- .trunc_bigq(r2)
  if (i == 0) return(S_Zero)
  I0 <- Integer(gmp::as.bigz(as.character(i)))
  if (gmp::denominator(r2) == 1) {
    if (isTRUE(.fz("StrictGreaterThan", r, I0))) return(I0)
    if (isTRUE(.fz("StrictLessThan", r, I0))) return(.sub(I0, S_One))
    ok <- .fz("Equality", r, I0)
    if (is.na(ok)) return(NULL)
    if (ok) return(I0)
    return(.sub(I0, if (i > 0) S_One else S_NegativeOne))
  }
  I0
}
.trunc_bigq <- function(q) {
  f <- gmp::as.bigq(.bigq_floor(q))
  if (q < 0 && f != q) f + 1 else f
}

# Mod's _eval_Mod hooks (Pow, factorial, binomial, partition).
.mod_hook <- function(p, q) {
  if (p$type == "Pow") {
    b <- p$args[[1L]]; e <- p$args[[2L]]
    if (isTRUE(ask_fact(e, "integer")) && isTRUE(ask_fact(e, "positive"))) {
      if (isTRUE(ask_fact(q, "integer")) && is_zero_num(fcall("Mod", list(b, q)))) return(S_Zero)
      if (b$type == "Integer" && e$type == "Integer" && q$type == "Integer")
        return(Integer(gmp::powm(gmp::as.bigz(b$value), gmp::as.bigz(e$value), gmp::as.bigz(q$value))))
    }
    return(NULL)
  }
  if (p$type == "Function" && identical(p$name, "factorial")) {
    n <- p$args[[1L]]
    if (isTRUE(ask_fact(n, "integer")) && isTRUE(ask_fact(n, "nonnegative")) && isTRUE(ask_fact(q, "integer"))) {
      aq <- fcall("Abs", list(q))
      d <- .sub(aq, n)
      if (isTRUE(ask_fact(d, "nonpositive"))) return(S_Zero)
      isprime <- ask_fact(aq, "prime")
      if (is_one_num(d)) {
        if (isTRUE(isprime)) return(fcall("Mod", list(S_NegativeOne, q)))
        if (identical(isprime, FALSE) && isTRUE(ask_fact(.sub(aq, Integer(6)), "nonnegative"))) return(S_Zero)
      } else if (n$type == "Integer" && q$type == "Integer") {
        m <- gmp::as.bigz(aq$value); nn <- as.numeric(n$value)
        fc <- gmp::as.bigz(1)
        dd <- as.numeric(d$value)
        if (isTRUE(isprime) && dd - 1 < nn) {
          for (k in seq_len(dd - 1)) fc <- (fc * k) %% m
          fc <- gmp::powm(fc, m - 2, m)
          if (dd %% 2 == 1) fc <- -fc
        } else {
          for (k in seq_len(nn)) fc <- (fc * k) %% m
        }
        return(fcall("Mod", list(Integer(fc), q)))
      }
    }
  }
  if (p$type == "Function" && identical(p$name, "binomial")) {
    # binomial._eval_Mod; its all-Integer branch needs an unevaluated
    # binomial, which R never builds (binomial of integers evaluates)
    if (any(vapply(list(p$args[[1L]], p$args[[2L]], q), function(z) identical(ask_fact(z, "integer"), FALSE), TRUE)))
      stop(.sympy_condition("ValueError", "Integers expected for binomial Mod"))
    return(NULL)
  }
  if (p$type == "Function" && identical(p$name, "partition") && q$type == "Integer") {
    # partition._eval_Mod: Ramanujan's congruences
    n <- p$args[[1L]]
    for (pr in list(c(5, 4), c(7, 5), c(11, 6)))
      if (q$value == pr[1L] && identical(fcall("Mod", list(n, q))$key, Integer(pr[2L])$key)) return(S_Zero)
  }
  NULL
}

# Mod.eval's number_eval: p %% q when both are numbers, or when p is known
# to lie within two multiples of q.
.mod_number_eval <- function(p, q) {
  if (is_zero_num(q)) stop(.sympy_condition("ZeroDivisionError", "Modulo by zero"))
  if (is_nan(p) || is_nan(q) || identical(ask_fact(p, "finite"), FALSE) ||
      identical(ask_fact(q, "finite"), FALSE)) return(S_NaN)
  if (is_zero_num(p) || identical(p$key, q$key) || identical(p$key, .neg(q)$key) ||
      (isTRUE(ask_fact(p, "integer")) && is_one_num(q))) return(S_Zero)
  if (is_number_type(q)) {
    if (is_number_type(p)) {
      if (is_exact_number(p) && is_exact_number(q))
        return(make_rational(p$value - q$value * gmp::as.bigq(.bigq_floor(p$value / q$value))))
      a <- num_value_dbl(p); b <- num_value_dbl(q)
      return(Float(a - b * floor(a / b)))
    }
    if (q$type == "Integer" && q$value == 2) {
      if (isTRUE(ask_fact(p, "even"))) return(S_Zero)
      if (isTRUE(ask_fact(p, "odd"))) return(S_One)
    }
  }
  rv <- .mod_hook(p, q)
  if (!is.null(rv)) return(rv)
  r <- .div(p, q)
  if (isTRUE(ask_fact(r, "integer"))) return(S_Zero)
  d <- .py_int_expr(r)
  if (!is.null(d)) {
    rv <- .sub(p, .mul(list(d, q)))
    if (isTRUE(.fz("StrictLessThan", .mul(list(rv, q)), S_Zero))) rv <- .add(list(rv, q))
    return(rv)
  }
  if (isTRUE(ask_fact(q, "positive"))) { c1 <- "LessThan"; c2 <- "StrictLessThan" }
  else if (isTRUE(ask_fact(q, "negative"))) { c1 <- "GreaterThan"; c2 <- "StrictGreaterThan" }
  else return(NULL)
  ls <- .mul(list(Integer(-2), q))
  r <- .sub(p, q)
  for (i in 1:4) {
    if (!isTRUE(.fz(c1, ls, p))) return(NULL)
    if (isTRUE(.fz(c2, r, ls))) return(.sub(p, ls))
    ls <- .add(list(ls, q))
  }
  NULL
}

# SymPy's Mod.eval.
.mod_eval <- function(p0, q0) {
  rv <- .mod_number_eval(p0, q0)
  if (!is.null(rv)) return(rv)
  p <- p0; q <- q0
  np <- .neg(p)
  if (.is_mod(p)) {
    qi <- p$args[[2L]]
    if (is_zero_num(fcall("Mod", list(qi, q)))) return(fcall("Mod", list(p$args[[1L]], q)))
    if (isTRUE(ask_fact(.mul(list(qi, .sub(q, qi))), "nonnegative"))) return(p)
  } else if (.is_mod(np)) {
    qi <- np$args[[2L]]
    if (is_zero_num(fcall("Mod", list(qi, q)))) return(fcall("Mod", list(.neg(np$args[[1L]]), q)))
    if (isTRUE(ask_fact(.mul(list(qi, .add(list(q, qi)))), "nonpositive"))) return(p)
  } else if (p$type == "Add") {
    ism <- vapply(p$args, .is_mod, TRUE)
    if (any(ism) && all(vapply(p$args[ism], function(m) identical(m$args[[2L]]$key, q$key), TRUE))) {
      net <- .add(c(p$args[!ism], lapply(p$args[ism], function(m) m$args[[1L]])))
      return(fcall("Mod", list(net, q)))
    }
  } else if (p$type == "Mul") {
    ism <- vapply(p$args, .is_mod, TRUE)
    non_mod_l <- p$args[!ism]; mod_l <- p$args[ism]
    allint <- all(vapply(p$args, function(t) isTRUE(ask_fact(t, "integer")), TRUE))
    if (length(mod_l) && all(vapply(mod_l, function(m) identical(m$args[[2L]]$key, q$key), TRUE)) &&
        allint && isTRUE(ask_fact(q, "integer"))) {
      nm <- lapply(non_mod_l, function(x) fcall("Mod", list(x, q)))
      mod <- list(); non_mod <- list()
      for (j in nm) if (.is_mod(j)) mod[[length(mod) + 1L]] <- j$args[[1L]] else non_mod[[length(non_mod) + 1L]] <- j
      net <- .mul(list(.mul(lapply(mod_l, function(i) i$args[[1L]])), .mul(mod), .mul(non_mod)))
      return(fcall("Mod", list(net, q)))
    }
    if (q$type == "Integer" && !is_one_num(q) && allint) {
      non_mod_l <- lapply(non_mod_l, function(i) if (i$type == "Integer") fcall("Mod", list(i, q)) else i)
      if (any(vapply(non_mod_l, is_zero_num, TRUE))) return(S_Zero)
    }
    p <- .mul(c(non_mod_l, mod_l))
  }
  # extract gcd; any further simplification is left to the user
  G <- tryCatch(gcd(p, q), error = function(e) S_One)
  if (!.equal_valued_one(G)) {
    p <- gcd_terms(.div(p, G), clear = FALSE, fraction = FALSE)
    q <- gcd_terms(.div(q, G), clear = FALSE, fraction = FALSE)
  }
  pwas <- p; qwas <- q
  if (p$type == "Add") {
    args <- lapply(p$args, function(i) {
      a <- fcall("Mod", list(i, q))
      if (.count_mod(a) > .count_mod(i)) i else a
    })
    if (!identical(vapply(args, `[[`, "", "key"), vapply(p$args, `[[`, "", "key"))) p <- .add(args)
  } else {
    sp <- as_coeff_mul_split(p); cp <- sp[[1L]]; p <- sp[[2L]]
    sq <- as_coeff_mul_split(q); cq <- sq[[1L]]; q <- sq[[2L]]
    ok <- FALSE
    if (!is_exact_number(cp) || !is_exact_number(cq)) {
      r <- .mod_number_eval(cp, cq)
      if (!is.null(r) && is_number_type(r) && isTRUE(num_value_dbl(r) == 0)) {
        G <- .mul(list(G, cq))
        p <- .mul(list(p, Integer(trunc(num_value_dbl(cp) / num_value_dbl(cq)))))
        ok <- TRUE
      }
    }
    if (!ok) { p <- .mul(list(cp, p)); q <- .mul(list(cq, q)) }
  }
  if (could_extract_minus_sign(p) && could_extract_minus_sign(q)) {
    G <- .neg(G); p <- .neg(p); q <- .neg(q)
  }
  rv <- .mod_number_eval(p, q)
  if (!is.null(rv)) return(.mul(list(rv, G)))
  if (G$type == "Float" && .equal_valued_one(G))
    return(fnode("Mod", list(.mul(list(p, G)), q)))
  if (G$type == "Mul" && G$args[[1L]]$type == "Float" && .equal_valued_one(G$args[[1L]])) {
    p <- .mul(list(G$args[[1L]], p))
    G <- .mul(G$args[-1L])
  }
  changed <- !identical(p$key, pwas$key) || !identical(q$key, qwas$key)
  if (is_one_num(G) && !changed && identical(p$key, p0$key) && identical(q$key, q0$key)) return(NULL)
  .mul(list(G, if (changed) fcall("Mod", list(p, q)) else fnode("Mod", list(p, q))))
}

# gmp's divq.bigz already rounds toward minus infinity.
.bigq_floor <- function(q) gmp::divq.bigz(gmp::numerator(q), gmp::denominator(q))
