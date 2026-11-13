# Core expression representation.
#
# Every expression is an immutable S3 list of class "sym" with fields
# type, args, name, value, extra and key. The key is a canonical string that
# encodes the whole tree; structural equality and canonical argument order are
# both defined through it. Only new_sym() writes keys and only sym_key() reads
# them outside this file, so the key scheme can be swapped for fingerprints.
# ponytail: recursive string keys cost O(n^2) memory on very deep chains;
# switch to hashed fingerprints if deep trees become common.

# Integer string to bigz. gmp reads a leading 0 as octal ("016" is 14), so
# decimal digits drop their leading zeros first; 0x, 0o and 0b prefixes are
# read as in Python.
.str_bigz <- function(s) {
  s <- gsub("_", "", trimws(s), fixed = TRUE)
  s <- sub("^\\+", "", s)
  if (grepl("^-?0[oO]", s)) s <- sub("0[oO]", "0", s)
  else if (!grepl("^-?0[xXbB]", s)) s <- sub("^(-?)0+([0-9])", "\\1\\2", s)
  z <- gmp::as.bigz(s)
  if (is.na(z)) stop("invalid integer literal: ", s, call. = FALSE)
  z
}

# Re-run the caller's own call, arguments unevaluated, with `fn` (a quoted
# function such as quote(stats::var)) in place of the masking function.
# The eval only replays the user's own call in the user's frame, so it runs
# no code the user did not write. Callers must not force any argument first,
# or that argument is evaluated twice.
.base_call <- function(fn) {
  cl <- sys.call(-1L)
  cl[[1L]] <- fn
  eval(cl, parent.frame(2L))
}

new_sym <- function(type, args = list(), name = NULL, value = NULL,
                    extra = NULL, key = NULL) {
  if (is.null(key)) key <- .compound_key(type, args, name)
  obj <- list(type = type, args = args, name = name, value = value,
              extra = extra, key = key)
  class(obj) <- "sym"
  obj
}

.compound_key <- function(type, args, name = NULL) {
  ks <- if (length(args)) vapply(args, .k, "") else character()
  paste0(type, if (!is.null(name)) paste0("[", nchar(name), ":", name, "]"),
         "(", paste(ks, collapse = ","), ")")
}

.k <- function(e) e$key

sym_key <- function(e) e$key

#' @export
`$<-.sym` <- function(x, name, value) {
  stop("sym objects are immutable", call. = FALSE)
}

#' @export
`[[<-.sym` <- function(x, i, value) {
  stop("sym objects are immutable", call. = FALSE)
}

#' Test whether an object is a symbolic expression
#'
#' @param x Any R object.
#' @return `TRUE` when `x` is a `sym` expression.
#' @export
is_sym <- function(x) inherits(x, "sym")

#' Type of a symbolic expression
#'
#' Returns the head of the expression tree, such as `"Add"`, `"Mul"`,
#' `"Pow"`, `"Symbol"`, `"Integer"` or `"Function"`. This plays the role of
#' SymPy's `func` and class checks such as `expr.is_Add`.
#'
#' @param e A `sym` expression.
#' @return A character scalar.
#' @export
sym_type <- function(e) as_sym(e)$type

#' Arguments of a symbolic expression
#'
#' The SymPy `args` attribute.
#'
#' @param e A `sym` expression.
#' @return A list of `sym` expressions.
#' @export
sym_args <- function(e) {
  # a mutable Matrix has no .args in SymPy
  if (inherits(e, "symmatrix"))
    stop(.sympy_condition("AttributeError", "'MutableDenseMatrix' object has no attribute 'args'"))
  .args_sympy(as_sym(e))
}

# the args of a node as SymPy's .args gives them (Add and Mul sorted by
# Basic.compare, see .sympy_args)
.args_sympy <- function(e) if (e$type %in% c("Add", "Mul")) .sympy_args(e) else e$args

#' Name of a symbol or function head
#'
#' @param e A `sym` expression.
#' @return The name as a character scalar, or `NULL` for unnamed heads.
#' @export
sym_name <- function(e) as_sym(e)$name

# ---------------------------------------------------------------------------
# Numbers

.NUMBER_TYPES <- c("Integer", "Rational", "Float")

is_number_type <- function(e) e$type %in% .NUMBER_TYPES
is_exact_number <- function(e) e$type == "Integer" || e$type == "Rational"

make_rational <- function(q) {
  den <- gmp::denominator(q)
  if (den == 1) {
    s <- as.character(gmp::numerator(q))
    new_sym("Integer", value = q, key = paste0("Z", s))
  } else {
    new_sym("Rational", value = q, key = paste0("Q", as.character(q)))
  }
}

#' Exact integers, rationals and floating point numbers
#'
#' `Integer()` builds an exact integer of arbitrary size, `Rational()` an
#' exact fraction reduced to lowest terms and `Float()` an arbitrary
#' precision binary floating point number, as SymPy's `Float`: 15 significant
#' digits by default, or `dps` digits. Strings are the safe way to pass
#' integers larger than `2^53`, and decimals with more than 15 digits, since R
#' doubles round them.
#'
#' @param p Numerator: a number, string, or `gmp` big integer.
#' @param q Denominator.
#' @param x A number or numeric string.
#' @param dps Decimal digits of precision for `Float()`.
#' @param precision Binary precision in bits for `Float()`; overrides `dps`.
#' @return A `sym` number.
#' @examples
#' Integer("123456789012345678901234567890")
#' Rational(6, 4)
#' Float(0.1)
#' Float("3.14159265358979323846264338327950288", 30)
#' @export
Integer <- function(p) {
  if (is_sym(p)) {
    if (p$type == "Integer") return(p)
    if (p$type == "Float" && p$value == round(p$value))
      return(make_rational(gmp::as.bigq(round(p$value))))
    stop("not an integer: ", format(p))
  }
  if (is.character(p)) return(make_rational(gmp::as.bigq(.str_bigz(p))))
  if (is.numeric(p) && p != round(p)) stop("not an integer: ", p)
  make_rational(gmp::as.bigq(p))
}

#' @rdname Integer
#' @export
Rational <- function(p, q = 1) {
  pq <- .to_bigq(p)
  qq <- .to_bigq(q)
  if (qq == 0) {
    if (pq == 0) return(S_NaN)
    return(S_ComplexInfinity)
  }
  make_rational(pq / qq)
}

.to_bigq <- function(x) {
  if (is_sym(x)) {
    if (is_exact_number(x)) return(x$value)
    stop("expected an exact number, got ", format(x))
  }
  if (inherits(x, "bigq")) return(x)
  if (inherits(x, "bigz")) return(gmp::as.bigq(x))
  if (is.character(x)) {
    if (grepl("/", x, fixed = TRUE)) {
      parts <- strsplit(x, "/", fixed = TRUE)[[1]]
      return(gmp::as.bigq(.str_bigz(parts[1]), .str_bigz(parts[2])))
    }
    if (grepl("[.eE]", x)) return(.decimal_to_bigq(x))
    return(gmp::as.bigq(.str_bigz(x)))
  }
  if (is.numeric(x)) {
    if (x == round(x) && abs(x) < 2^53) return(gmp::as.bigq(x))
    return(gmp::as.bigq(x))
  }
  stop("cannot convert to rational: ", class(x)[1])
}

# Exact value of a decimal string such as "0.1" or "1.5e-3".
.decimal_to_bigq <- function(s) {
  s <- tolower(trimws(s))
  ex <- 0L
  if (grepl("e", s, fixed = TRUE)) {
    parts <- strsplit(s, "e", fixed = TRUE)[[1]]
    s <- parts[1]
    ex <- as.integer(parts[2])
  }
  neg <- startsWith(s, "-")
  s <- sub("^[-+]", "", s)
  frac <- 0L
  if (grepl(".", s, fixed = TRUE)) {
    parts <- strsplit(s, ".", fixed = TRUE)[[1]]
    if (length(parts) == 1) parts <- c(parts, "")
    frac <- nchar(parts[2])
    s <- paste0(parts[1], parts[2])
  }
  if (s == "") s <- "0"
  num <- gmp::as.bigq(.str_bigz(s))
  pw <- ex - frac
  ten <- gmp::as.bigq(10)
  r <- if (pw >= 0) num * ten^pw else num / ten^(-pw)
  if (neg) -r else r
}

#' @rdname Integer
#' @export
Float <- function(x, dps = NULL, precision = NULL) {
  # SymPy's Float(x, dps): a decimal string is read exactly, a double keeps
  # its binary value, both rounded to the precision
  if (!is.null(dps) || !is.null(precision) || (is.character(x) && nchar(gsub("[^0-9]", "", sub("[eE].*$", "", x))) > 15)) {
    prec <- if (!is.null(precision)) precision else if (!is.null(dps)) .dps_to_prec(dps)
            else .dps_to_prec(nchar(sub("^0+", "", gsub("[^0-9]", "", sub("[eE].*$", "", x)))))
    q <- if (is.character(x)) .decimal_to_bigq(x) else if (is_sym(x) && x$type == "Float") .mp_to_q(.float_mp(x))
         else if (is_sym(x)) gmp::as.bigq(x$value) else .mp_to_q(.mp_from_dbl(as.double(x)))
    return(.mp_float(.mp_from_q(q, prec), prec))
  }
  # a decimal string is rounded correctly to 53 bits (R's own strtod is not
  # correctly rounded, for example for "1e-200")
  if (is.character(x) && grepl("^\\s*[-+]?([0-9]+\\.?[0-9]*|\\.[0-9]+)([eE][-+]?[0-9]+)?\\s*$", x))
    return(.mp_float(.mp_from_q(.decimal_to_bigq(trimws(x)), 53), 53))
  if (is_sym(x)) {
    if (is_inf_const(x) || is_nan(x)) return(x)
    v <- num_value_dbl(x)
    if (is.na(v)) stop("not a number: ", format(x))
    x <- v
  }
  x <- as.double(x)
  # SymPy's Float of an infinity or NaN is oo, -oo or nan
  if (is.nan(x)) return(S_NaN)
  if (is.infinite(x)) return(if (x > 0) S_Infinity else S_NegativeInfinity)
  new_sym("Float", value = x, key = paste0("F", sprintf("%.17g", x)))
}

# A Float node as is: arithmetic on Floats keeps the Float type even when a
# double overflows, so coefficient accumulators stay numeric.
.float_raw <- function(x) new_sym("Float", value = x, key = paste0("F", sprintf("%.17g", x)))

num_value_dbl <- function(e) {
  switch(e$type,
         Integer = , Rational = .bigq_dbl(e$value),
         Float = e$value,
         NA_real_)
}

# A rational as the nearest double, as mpmath's from_rational does; gmp's
# as.double() truncates (mpq_get_d), so 1/1000 would be one ulp low.
.bigq_dbl <- function(q) {
  n <- gmp::numerator(q); d <- gmp::denominator(q)
  if (abs(n) < 2^53 && d < 2^53) return(as.double(n) / as.double(d))
  r <- .mp_to_dbl(.mp_from_q(q, 53))
  # .mp_to_dbl underflows below the normal range
  if (r == 0 || is.infinite(r)) as.double(q) else r
}

num_sign <- function(e) {
  if (e$type == "Float") {
    if (!is.null(e$extra$mp)) { m <- e$extra$mp$m; return(if (m > 0) 1 else if (m < 0) -1 else 0) }
    return(sign(e$value))
  }
  v <- e$value
  if (v > 0) 1 else if (v < 0) -1 else 0
}

num_is_zero <- function(e) num_sign(e) == 0
num_is_one <- function(e) {
  if (e$type == "Float") return(if (!is.null(e$extra$mp)) e$extra$mp$m == 1 && e$extra$mp$e == 0 else e$value == 1)
  e$type == "Integer" && e$value == 1
}

num_add <- function(a, b) {
  if ((a$type == "Float" && !is.null(a$extra$mp)) || (b$type == "Float" && !is.null(b$extra$mp))) return(.num_mp_op("add", a, b))
  if (a$type == "Float" || b$type == "Float")
    return(.float_raw(num_value_dbl(a) + num_value_dbl(b)))
  make_rational(a$value + b$value)
}

num_mul <- function(a, b) {
  if ((a$type == "Float" && !is.null(a$extra$mp)) || (b$type == "Float" && !is.null(b$extra$mp))) return(.num_mp_op("mul", a, b))
  if (a$type == "Float" || b$type == "Float")
    return(.float_raw(num_value_dbl(a) * num_value_dbl(b)))
  make_rational(a$value * b$value)
}

num_neg <- function(a) {
  if (a$type == "Float" && !is.null(a$extra$mp)) return(.mp_float(.mp_neg(a$extra$mp), a$extra$prec))
  if (a$type == "Float") return(.float_raw(-a$value))
  make_rational(-a$value)
}

num_inv <- function(a) {
  if (a$type == "Float" && !is.null(a$extra$mp)) return(.num_mp_op("div", S_One, a))
  if (a$type == "Float") return(.float_raw(1 / a$value))
  if (a$value == 0) return(S_ComplexInfinity)
  make_rational(1 / a$value)
}

num_cmp <- function(a, b) {
  if ((a$type == "Float" && !is.null(a$extra$mp)) || (b$type == "Float" && !is.null(b$extra$mp))) {
    ma <- if (a$type == "Float") .float_mp(a) else .mp_from_q(a$value, 4000)
    mb <- if (b$type == "Float") .float_mp(b) else .mp_from_q(b$value, 4000)
    return(.mp_cmp(ma, mb))
  }
  if (a$type == "Float" || b$type == "Float") {
    d <- num_value_dbl(a) - num_value_dbl(b)
    return(sign(d))
  }
  d <- a$value - b$value
  if (d > 0) 1 else if (d < 0) -1 else 0
}

num_numer <- function(e) gmp::numerator(e$value)
num_denom <- function(e) gmp::denominator(e$value)

.bigq_pow_int <- function(q, n) {
  n <- as.integer(n)
  if (n == 0) return(gmp::as.bigq(1))
  if (n > 0) return(q^n)
  (1 / q)^(-n)
}

# ---------------------------------------------------------------------------
# Constants and singletons

new_constant <- function(name) new_sym("Constant", name = name,
                                       key = paste0("C", name))

S_Zero <- make_rational(gmp::as.bigq(0))
S_One <- make_rational(gmp::as.bigq(1))
S_NegativeOne <- make_rational(gmp::as.bigq(-1))
S_Two <- make_rational(gmp::as.bigq(2))
S_Half <- make_rational(gmp::as.bigq(1, 2))
S_Pi <- new_constant("pi")
S_E <- new_constant("E")
S_I <- new_constant("I")
S_Infinity <- new_constant("oo")
S_NegativeInfinity <- new_constant("-oo")
S_ComplexInfinity <- new_constant("zoo")
S_NaN <- new_constant("nan")
S_EulerGamma <- new_constant("EulerGamma")
S_Catalan <- new_constant("Catalan")
S_GoldenRatio <- new_constant("GoldenRatio")
S_TribonacciConstant <- new_constant("TribonacciConstant")
S_true <- new_sym("BooleanAtom", name = "true", key = "Btrue")
S_false <- new_sym("BooleanAtom", name = "false", key = "Bfalse")
S_EmptySet <- new_sym("EmptySet", key = "EmptySet()")

.INF_NAMES <- c("oo", "-oo", "zoo")

is_inf_const <- function(e) e$type == "Constant" && e$name %in% .INF_NAMES
is_nan <- function(e) e$type == "Constant" && e$name == "nan"
is_zero_num <- function(e) is_exact_number(e) && e$value == 0
is_one_num <- function(e) e$type == "Integer" && e$value == 1
is_negone_num <- function(e) e$type == "Integer" && e$value == -1
is_int_num <- function(e) e$type == "Integer"
is_const <- function(e, nm) e$type == "Constant" && e$name == nm

#' Mathematical constants and special values
#'
#' `pi`, Euler's number `E`, the imaginary unit `I`, positive and negative
#' infinity `oo` and `-oo`, complex infinity `zoo`, not a number `nan`, and the
#' constants `EulerGamma`, `Catalan`, `GoldenRatio` and `TribonacciConstant`.
#' They are exported as objects of the same names prefixed with `S_` (for
#' example `S_Pi`) and, for convenience, as `sym_pi`, `sym_E`, `sym_I`, `oo`,
#' `zoo` and `nan`.
#'
#' @name constants
#' @format `sym` objects.
#' @examples
#' sym_pi + 1
#' sym_I^2
#' oo + 1
NULL

#' @rdname constants
#' @export
sym_pi <- S_Pi
#' @rdname constants
#' @export
sym_E <- S_E
#' @rdname constants
#' @export
sym_I <- S_I
#' @rdname constants
#' @export
oo <- S_Infinity
#' @rdname constants
#' @export
zoo <- S_ComplexInfinity
#' @rdname constants
#' @export
nan <- S_NaN
#' @rdname constants
#' @export
EulerGamma <- S_EulerGamma
#' @rdname constants
#' @export
Catalan <- S_Catalan
#' @rdname constants
#' @export
GoldenRatio <- S_GoldenRatio
#' @rdname constants
#' @export
TribonacciConstant <- S_TribonacciConstant
#' @rdname constants
#' @export
S_Pi <- S_Pi
#' @rdname constants
#' @export
S_E <- S_E
#' @rdname constants
#' @export
S_I <- S_I
#' @rdname constants
#' @export
S_Infinity <- S_Infinity
#' @rdname constants
#' @export
S_NegativeInfinity <- S_NegativeInfinity
#' @rdname constants
#' @export
S_ComplexInfinity <- S_ComplexInfinity
#' @rdname constants
#' @export
S_NaN <- S_NaN
#' @rdname constants
#' @export
S_Zero <- S_Zero
#' @rdname constants
#' @export
S_One <- S_One
#' @rdname constants
#' @export
S_Half <- S_Half
#' @rdname constants
#' @export
S_NegativeOne <- S_NegativeOne
#' @rdname constants
#' @export
S_true <- S_true
#' @rdname constants
#' @export
S_false <- S_false
#' @rdname constants
#' @export
S_EmptySet <- S_EmptySet

# ---------------------------------------------------------------------------
# Symbols

.ASSUMPTION_NAMES <- c(
  "commutative", "complex", "real", "rational", "irrational", "integer",
  "noninteger", "even", "odd", "prime", "composite", "positive", "negative",
  "zero", "nonzero", "nonnegative", "nonpositive", "finite", "infinite",
  "imaginary", "algebraic", "transcendental", "extended_real",
  "extended_positive", "extended_negative", "extended_nonnegative",
  "extended_nonpositive", "extended_nonzero", "hermitian", "antihermitian",
  "polar"
)

.assumption_suffix <- function(facts) {
  if (all(is.na(facts[names(facts) != "commutative"])) &&
      isTRUE(facts[["commutative"]])) return("")
  paste0("{", paste(ifelse(is.na(facts), ".", ifelse(facts, "T", "F")),
                    collapse = ""), "}")
}

#' Create symbols
#'
#' `Symbol()` creates one symbol. `sym_symbols()` splits a string on spaces or
#' commas and returns one symbol or a list of them; ranges such as `"x0:3"`
#' expand to `x0, x1, x2`. `sym_var()` does the same and assigns the symbols in
#' the calling environment. Assumptions such as `positive = TRUE` or
#' `integer = TRUE` are recorded on the symbol and drive simplification, as in
#' SymPy. Symbols with the same name and different assumptions are distinct.
#'
#' @param name Symbol name.
#' @param names A string of names separated by spaces or commas.
#' @param ... Assumptions given as `name = TRUE` or `name = FALSE`. Valid
#'   names are `commutative`, `complex`, `real`, `rational`, `irrational`,
#'   `integer`, `noninteger`, `even`, `odd`, `prime`, `composite`,
#'   `positive`, `negative`, `zero`, `nonzero`, `nonnegative`, `nonpositive`,
#'   `finite`, `infinite`, `imaginary`, `algebraic`, `transcendental`,
#'   `extended_real`, `hermitian`, `antihermitian` and `polar`, and the
#'   Fermi level assumptions `above_fermi` and `below_fermi` used by
#'   [wicks()].
#' @param cls For `sym_symbols()`, the constructor to use instead of
#'   `Symbol()`, such as `Dummy` or `"Dummy"`.
#' @param envir Environment in which `sym_var()` assigns.
#' @return A `sym` symbol, or a list of symbols.
#' @examples
#' x <- Symbol("x")
#' xyz <- sym_symbols("x y z")
#' n <- Symbol("n", integer = TRUE, positive = TRUE)
#' @export
Symbol <- function(name, ...) {
  sp <- .split_other_facts(list(...))
  given <- .normalize_facts(sp$rest)
  facts <- .symbol_facts(given)
  if (identical(facts[["commutative"]], FALSE)) .nc_state$any <- TRUE
  extra <- list(facts = facts, given = given)
  if (length(sp$other)) extra$other <- sp$other
  new_sym("Symbol", name = name, extra = extra,
          key = paste0("S", nchar(name), ":", name, .assumption_suffix(facts),
                       .other_facts_suffix(sp$other)))
}

# Assumptions outside the fact system, such as the Fermi level ones of
# sympy.physics.secondquant (or any other name: Symbol does not check them),
# are kept as given (assumptions0) and are part of the identity of the symbol
.split_other_facts <- function(a) {
  is_o <- !(names(a) %in% .ASSUMPTION_NAMES)
  o <- a[is_o]
  if (length(o)) {
    o <- vapply(o, function(z) isTRUE(as.logical(z)), TRUE)
    o <- o[order(names(o), method = "radix")]
  }
  list(rest = a[!is_o], other = o)
}

.other_facts_suffix <- function(o) {
  if (!length(o)) return("")
  paste0("[", paste0(names(o), "=", ifelse(o, "T", "F"), collapse = ","), "]")
}

.normalize_facts <- function(f) {
  if (!length(f)) return(logical())
  bad <- setdiff(names(f), .ASSUMPTION_NAMES)
  if (length(bad)) stop("unknown assumption: ", paste(bad, collapse = ", "))
  v <- vapply(f, function(z) isTRUE(as.logical(z)), TRUE)
  names(v) <- names(f)
  v
}

.expand_symbol_names <- function(names) {
  parts <- strsplit(trimws(names), "[, ]+")[[1]]
  parts <- parts[nzchar(parts)]
  out <- character()
  for (p in parts) {
    m <- regmatches(p, regexec("^([^:]*?)([0-9]*):([0-9]+)$", p))[[1]]
    if (length(m)) {
      start <- if (nzchar(m[3])) as.integer(m[3]) else 0L
      end <- as.integer(m[4]) - 1L
      if (end >= start) out <- c(out, paste0(m[2], start:end))
    } else {
      out <- c(out, p)
    }
  }
  out
}

#' @noRd
symbols <- function(names, ..., cls = NULL) {
  if (missing(names)) return(.base_call(quote(graphics::symbols)))
  if (!is.character(names)) return(graphics::symbols(names, ...))
  nms <- .expand_symbol_names(names)
  mk <- if (is.null(cls)) Symbol else if (is.character(cls)) get(cls, mode = "function") else cls
  res <- lapply(nms, function(n) mk(n, ...))
  if (length(res) == 1 && !grepl("[, :]", trimws(names))) return(res[[1]])
  names(res) <- nms
  res
}

#' @noRd
var <- function(names, ..., envir = parent.frame()) {
  if (missing(names)) return(.base_call(quote(stats::var)))
  if (!is.character(names)) return(stats::var(names, ...))
  nms <- .expand_symbol_names(names)
  res <- lapply(nms, function(n) Symbol(n, ...))
  for (i in seq_along(nms)) assign(nms[i], res[[i]], envir = envir)
  invisible(if (length(res) == 1) res[[1]] else stats::setNames(res, nms))
}

.dummy_counter <- new.env()
.dummy_counter$n <- 0L

#' Dummy and wild symbols
#'
#' `Dummy()` creates a symbol that is unique even when another symbol has the
#' same name. `Wild()` creates a pattern variable for [sym_match()] and
#' [sym_replace()]; `exclude` lists expressions the wild may not contain and
#' `properties` is a list of predicate functions the match must satisfy.
#'
#' @param name Display name; by default `Dummy_<n>`, numbered by a session
#'   counter as SymPy does.
#' @param ... Assumptions, as for [Symbol()].
#' @param exclude A list of expressions the match must not contain.
#' @param properties A list of functions returning `TRUE` for valid matches.
#' @return A `sym` object.
#' @export
Dummy <- function(name = NULL, ...) {
  # Dummy._count: every Dummy takes the next number; an unnamed one is
  # named after it
  id <- .dummy_counter$n
  .dummy_counter$n <- id + 1L
  if (is.null(name)) name <- paste0("Dummy_", id)
  sp <- .split_other_facts(list(...))
  facts <- .symbol_facts(.normalize_facts(sp$rest))
  if (identical(facts[["commutative"]], FALSE)) .nc_state$any <- TRUE
  extra <- list(facts = facts, id = id)
  if (length(sp$other)) extra$other <- sp$other
  new_sym("Dummy", name = name, extra = extra,
          key = paste0("D", id, ":", nchar(name), ":", name,
                       .assumption_suffix(facts), .other_facts_suffix(sp$other)))
}

#' @rdname Dummy
#' @export
Wild <- function(name, exclude = list(), properties = list()) {
  exclude <- lapply(exclude, as_sym)
  new_sym("Wild", name = name,
          extra = list(exclude = exclude, properties = properties,
                       facts = .symbol_facts(logical())),
          key = paste0("W", nchar(name), ":", name, "{",
                       paste(vapply(exclude, .k, ""), collapse = ","), "}"))
}

is_symbol_like <- function(e) e$type %in% c("Symbol", "Dummy", "Wild")

# ---------------------------------------------------------------------------
# Conversion from R values

#' Convert R values to symbolic expressions
#'
#' `as_sym()` (alias `S()`) converts numbers, strings, logical values,
#' `gmp` big numbers, R complex numbers, formulas and calls into `sym`
#' expressions. Integer valued doubles become exact integers and other
#' doubles become [Float()] values. Strings are parsed with [sympify()].
#'
#' @param x Value to convert.
#' @return A `sym` expression.
#' @examples
#' S(2) / 3
#' S("x^2 + 1")
#' @export
as_sym <- function(x) {
  if (inherits(x, "sym")) return(x)
  # domain elements (polys-domains.R)
  if (inherits(x, "sym_ring_element")) return(x$expr)
  if (inherits(x, "sym_modint") || inherits(x, "sym_gaussian")) return(to_sympy(NULL, x))
  if (inherits(x, "bigz") || inherits(x, "bigq")) {
    if (length(x) != 1) stop("as_sym() needs a scalar")
    return(make_rational(gmp::as.bigq(x)))
  }
  if (is.character(x)) return(sympify(x))
  if (is.logical(x)) {
    if (length(x) != 1 || is.na(x)) stop("as_sym() needs TRUE or FALSE")
    return(if (x) S_true else S_false)
  }
  if (is.complex(x)) {
    if (length(x) != 1) stop("as_sym() needs a scalar")
    re <- as_sym(Re(x))
    im <- as_sym(Im(x))
    return(.add(list(re, .mul(list(im, S_I)))))
  }
  if (is.numeric(x)) {
    if (length(x) != 1) stop("as_sym() needs a scalar; use lapply() for vectors")
    if (is.na(x)) return(S_NaN)
    if (is.infinite(x)) return(if (x > 0) S_Infinity else S_NegativeInfinity)
    if (is.integer(x) || (x == round(x) && abs(x) < 2^53))
      return(make_rational(gmp::as.bigq(x)))
    return(Float(x))
  }
  if (is.language(x) || is.symbol(x)) return(sym_from_call(x))
  if (inherits(x, "symmatrix")) return(x)
  if (inherits(x, "sym_predicate_fn")) return(attr(x, "predicate"))
  stop("cannot convert object of class '", class(x)[1], "' to sym")
}

#' @rdname as_sym
#' @export
S <- function(x) as_sym(x)

# ---------------------------------------------------------------------------
# Generic rebuild: reconstruct an expression of the same head with new args.


sym_rebuild <- function(e, args) {
  t <- e$type
  if (t == "Add") return(.add(args))
  if (t == "Mul") return(.mul(args))
  if (t == "Pow") return(.pow(args[[1]], args[[2]]))
  # an undefined function marked by dynamicsymbols() keeps its mark; an
  # implemented function keeps its implementation
  if (t == "Function") {
    if (isTRUE(e$extra$undefined)) return(.pv_undef(e$name, args, e$extra))
    return(if (is.null(e$extra$imp)) fcall(e$name, args) else .cg_imp_call(e$name, e$extra$imp, args))
  }
  f <- .rebuilders[[t]]
  if (!is.null(f)) return(f(e, args))
  if (!length(args)) return(e)
  new_sym(t, args, name = e$name, value = e$value, extra = e$extra)
}

is_atom <- function(e) length(e$args) == 0
