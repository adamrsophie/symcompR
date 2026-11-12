# Sets, after sympy/sets.


is_set <- function(e) e$type %in% .SET_TYPES

.special_set <- function(name) new_sym("SpecialSet", name = name, key = paste0("SS", name))

#' Sets
#'
#' Intervals, finite sets, the standard number sets and set operations, like
#' SymPy's sets module. `Interval(a, b)` is closed by default; use
#' `left_open`/`right_open` or `Interval_open()`, `Interval_Lopen()` and
#' `Interval_Ropen()`. Number sets are `Reals`, `Integers`, `Naturals`
#' (starting at 1), `Naturals0`, `Rationals`, `Complexes` and
#' `UniversalSet`. `Union()`, `Intersection()`, `Complement()`,
#' `SymmetricDifference()` and `ProductSet()` combine sets and simplify when
#' possible. `contains()` tests membership; `measure()`, `boundary()`,
#' `set_inf()`, `set_sup()` and `is_subset()` query sets.
#'
#' @param a,b Endpoints.
#' @param left_open,right_open Openness of the endpoints.
#' @param ... Elements or sets.
#' @param s,t Sets.
#' @param x An element.
#' @return A set expression, or a logical/boolean value.
#' @examples
#' Interval(0, 1)
#' Union(Interval(0, 1), Interval(2, 3), Interval(1, 2))
#' Intersection(Interval(0, 5), Interval(3, 10))
#' contains(Interval(0, 1), Rational(1, 2))
#' Complement(Reals, FiniteSet(0))
#' @name sets
NULL

#' @rdname sets
#' @export
Reals <- .special_set("Reals")
#' @rdname sets
#' @export
Integers <- .special_set("Integers")
#' @rdname sets
#' @export
Naturals <- .special_set("Naturals")
#' @rdname sets
#' @export
Naturals0 <- .special_set("Naturals0")
#' @rdname sets
#' @export
Rationals <- .special_set("Rationals")
#' @rdname sets
#' @export
Complexes <- .special_set("Complexes")
#' @rdname sets
#' @export
UniversalSet <- .special_set("UniversalSet")

.lt_true <- function(a, b) isTRUE(.rel_to_r(tryCatch(Lt(a, b), error = function(e) NULL) %||% S_false))
.cmp <- function(a, b) {
  # -1, 0, 1 or NA
  if (identical(a$key, b$key)) return(0L)
  d <- .sub(a, b)
  if (is_const(a, "-oo") || is_const(b, "oo")) return(-1L)
  if (is_const(a, "oo") || is_const(b, "-oo")) return(1L)
  fd <- .facts(d)
  if (isTRUE(fd[["zero"]])) return(0L)
  if (isTRUE(fd[["positive"]])) return(1L)
  if (isTRUE(fd[["negative"]])) return(-1L)
  NA_integer_
}


#' @rdname sets
#' @export
Interval_open <- function(a, b) Interval(a, b, TRUE, TRUE)
#' @rdname sets
#' @export
Interval_Lopen <- function(a, b) Interval(a, b, TRUE, FALSE)
#' @rdname sets
#' @export
Interval_Ropen <- function(a, b) Interval(a, b, FALSE, TRUE)

#' @rdname sets
#' @export
FiniteSet <- function(...) {
  els <- list(...)
  if (length(els) == 1L && is.list(els[[1L]]) && !inherits(els[[1L]], "sym")) els <- els[[1L]]
  els <- lapply(els, as_sym)
  if (!length(els)) return(S_EmptySet)
  els <- els[!duplicated(vapply(els, .k, ""))]
  new_sym("FiniteSet", .sort_set_elems(els))
}

.sort_set_elems <- function(els) {
  if (length(els) < 2L) return(els)
  vals <- vapply(els, function(e) {
    if (length(free_symbols(e))) return(NA_real_)
    v <- tryCatch(.neval(e, list()), error = function(err) NA)
    if (length(v) != 1L || is.complex(v) || is.na(v)) NA_real_ else as.double(v)
  }, 0)
  o <- order(is.na(vals), vals, vapply(els, .k, ""), method = "radix")
  # Distinct numbers can share one double; order such runs exactly.
  v <- vals[o]
  runs <- rle(ifelse(is.na(v), NA_character_, sprintf("%.17g", v)))
  end <- cumsum(runs$lengths)
  for (r in which(runs$lengths > 1L & !is.na(runs$values))) {
    idx <- (end[r] - runs$lengths[r] + 1L):end[r]
    run <- els[o[idx]]
    if (!all(vapply(run, is_number_type, TRUE))) next
    q <- lapply(run, function(e) gmp::as.bigq(e$value))
    ro <- seq_along(run)
    for (i in seq_along(run)[-1L]) {
      j <- i - 1L
      cur <- ro[i]
      while (j > 0L && q[[ro[j]]] > q[[cur]]) {
        ro[j + 1L] <- ro[j]
        j <- j - 1L
      }
      ro[j + 1L] <- cur
    }
    o[idx] <- o[idx][ro]
  }
  els[o]
}

.finite_elems <- function(s) if (s$type == "FiniteSet") s$args else list()

#' @rdname sets
#' @export
EmptySet <- S_EmptySet

# SymPy default_sort_key (core/expr.py, Expr.sort_key), which the SymPy
# printers use to order FiniteSet elements. Keys are nested lists whose
# leaves are numbers or strings; .key_cmp compares them like Python tuples.
.FUNC_CLASS_INDEX <- c(exp = 10, log = 11, sin = 20, cos = 21, tan = 22, cot = 23, sinh = 30,
                       cosh = 31, tanh = 32, coth = 33, conjugate = 40, re = 41, im = 42, arg = 43)
# SymPy function classes whose nargs is Naturals0 (class_key index 0)
.NARGS0_FUNS <- c("carmichael", "exp_polar", "lerchphi", "meijerg", "LeviCivita")
# Max and Min are not Function subclasses: Basic.class_key
.BASIC_KEY_FUNS <- c("Max", "Min")
.NARGS0_TYPES <- c("Piecewise",
                   "CosineTransform", "FourierTransform", "HankelTransform", "InverseCosineTransform",
                   "InverseFourierTransform", "InverseHankelTransform", "InverseLaplaceTransform",
                   "InverseMellinTransform", "InverseSineTransform", "LaplaceTransform", "MellinTransform", "SineTransform")
.CONST_CLASS <- c(I = "ImaginaryUnit", pi = "Pi", E = "Exp1", zoo = "ComplexInfinity")

.class_key <- function(e) {
  t <- e$type
  # the number sets sort by their SymPy class name (Basic.class_key)
  if (t == "SpecialSet") return(list(5, 0, e$name))
  if (t %in% .NUMBER_TYPES || (t == "Constant" && e$name %in% c("oo", "-oo"))) return(list(1, 0, "Number"))
  if (t == "Constant") return(list(2, 0, if (is.na(.CONST_CLASS[e$name])) e$name else .CONST_CLASS[[e$name]]))
  if (t %in% c("Symbol", "Dummy", "Wild", "Object", "Str")) return(list(2, 0, if (!is.null(e$extra$sympy_class)) e$extra$sympy_class else t))
  # sympy.vector: base vectors and dyadics are atoms; their sums and
  # products subclass Add and Mul
  if (t %in% c("BaseVector", "BaseDyadic")) return(list(2, 0, t))
  if (t %in% c("VectorMul", "DyadicMul")) return(list(3, 0, t))
  if (t %in% c("VectorAdd", "DyadicAdd")) return(list(3, 1, t))
  # a physical quantity is an Atom
  if (t == "Quantity") return(list(2, 0, .sympy_class(e)))
  if (t == "Add") return(list(3, 1, "Add"))
  if (t == "Mul") return(list(3, 0, "Mul"))
  if (t == "Pow") return(list(3, 2, "Pow"))
  if (t == "Function" && e$name %in% .BASIC_KEY_FUNS) return(list(5, 0, e$name))
  if (t == "Function") {
    i <- .FUNC_CLASS_INDEX[e$name]
    return(list(4, if (!is.na(i)) i[[1L]] else if (.is_undefined_application(e) || e$name %in% .NARGS0_FUNS) 0 else 10000, e$name))
  }
  # Applications that are separate node types here
  if (t %in% .NARGS0_TYPES) return(list(4, 0, t))
  list(5, 0, t)
}

.ONE_SORT_KEY <- list(list(1, 0, "Number"), list(0, list()), list(), 1)

.sort_key <- function(e) {
  if (e$type == "Float" && is.nan(e$value)) e <- S_NaN
  if (e$type %in% .NUMBER_TYPES || (e$type == "Constant" && e$name %in% c("oo", "-oo")))
    return(list(.class_key(e), list(0, list()), list(),
                if (e$type == "Constant") (if (e$name == "oo") Inf else -Inf) else .sort_num(e)))
  if (e$type == "Constant" && e$name == "nan") return(list(list(1, 0, "NaN"), list(0, list()), list(), -Inf))
  # atoms, and Indexed (an atom in SymPy), sort by their string
  if (!length(e$args) || e$type == "Indexed") return(list(.class_key(e), list(1, list(.str(e))), .ONE_SORT_KEY, 1))
  # SymPy's MatrixSymbol args start with its name, a Str
  if (e$type == "MatrixSymbol")
    return(list(.class_key(e), list(3, list(list(list(5, 0, "Str"), list(1, list(e$name)), .ONE_SORT_KEY, 1),
                                           .sort_key(e$args[[1L]]), .sort_key(e$args[[2L]]))), .ONE_SORT_KEY, 1))
  sp <- as_coeff_mul_split(e)
  ex <- sp[[2L]]
  pw <- S_One
  if (ex$type == "Pow") {
    pw <- ex$args[[2L]]
    ex <- ex$args[[1L]]
  }
  args <- if (!length(ex$args)) list(.str(ex)) else {
    a <- if (ex$type == "Add") .add_print_order(ex$args) else
      if (ex$type %in% c("And", "Or")) .bool_print_order(ex$args) else ex$args
    k <- lapply(a, .sort_key)
    if (ex$type == "Mul") k <- .key_sort(k)
    k
  }
  list(.class_key(ex), list(length(args), args), .sort_key(pw), .sort_num(sp[[1L]]))
}

# exact rationals compare exactly (2^53 and 2^53 + 1 differ); floats as doubles
.sort_num <- function(e) {
  if (e$type %in% c("Integer", "Rational")) return(gmp::as.bigq(e$value))
  if (e$type == "Float" && (!is.null(e$extra$mp) || is.finite(e$value))) return(.mp_to_q(.float_mp(e)))
  num_value_dbl(e)
}

.key_cmp <- function(a, b) {
  if (is.list(a) && is.list(b)) {
    for (i in seq_len(min(length(a), length(b)))) {
      r <- .key_cmp(a[[i]], b[[i]])
      if (r != 0) return(r)
    }
    return(sign(length(a) - length(b)))
  }
  if (is.character(a) != is.character(b) || is.list(a) != is.list(b)) return(if (is.numeric(a) || inherits(a, "bigq")) -1 else 1)
  if (inherits(a, "bigq") || inherits(b, "bigq")) {
    if (!inherits(a, "bigq") && is.infinite(a)) return(if (a > 0) 1 else -1)
    if (!inherits(b, "bigq") && is.infinite(b)) return(if (b > 0) -1 else 1)
    return(if (a == b) 0 else if (a < b) -1 else 1)
  }
  if (a == b) 0 else if (is.character(a)) (if (order(c(a, b), method = "radix")[1L] == 1L) -1 else 1) else if (a < b) -1 else 1
}

# ponytail: insertion sort, O(n^2) key comparisons; fine for printed sets.
.key_sort <- function(keys, items = keys) {
  o <- integer()
  for (i in seq_along(keys)) {
    j <- length(o)
    while (j > 0L && .key_cmp(keys[[o[j]]], keys[[i]]) > 0) j <- j - 1L
    o <- append(o, i, after = j)
  }
  items[o]
}

register_str_printer("SpecialSet", function(e) e$name)
register_latex_printer("FiniteSet", function(e) {
  els <- .key_sort(lapply(e$args, .sort_key), e$args)
  paste0("\\left\\{", paste(vapply(els, .latex, ""), collapse = ", "), "\\right\\}")
})
register_latex_printer("EmptySet", function(e) "\\emptyset")
register_latex_printer("SpecialSet", function(e) {
  m <- c(Reals = "\\mathbb{R}", Complexes = "\\mathbb{C}", Integers = "\\mathbb{Z}",
         Naturals = "\\mathbb{N}", Naturals0 = "\\mathbb{N}_0", Rationals = "\\mathbb{Q}")
  if (is.na(m[e$name])) paste0("\\text{", e$name, "}") else m[[e$name]]
})
register_latex_printer("Interval", function(e) {
  i <- .iv(e)
  if (identical(i$a$key, i$b$key)) return(paste0("\\left\\{", .latex(i$a), "\\right\\}"))
  paste0("\\left", if (i$lo) "(" else "[", .latex(i$a), ", ", .latex(i$b), "\\right", if (i$ro) ")" else "]")
})
# SymPy parenthesizes a set operation nested in another one, since all of
# them share one precedence level.
.latex_setop_args <- function(e) vapply(if (e$type %in% c("Union", "Intersection")) .set_print_order(e$args) else e$args, function(a) {
  s <- .latex(a)
  if (a$type %in% c("Union", "Intersection", "Complement", "SymmetricDifference", "ProductSet"))
    paste0("\\left(", s, "\\right)") else s
}, "")
register_latex_printer("Union", function(e) paste(.latex_setop_args(e), collapse = " \\cup "))
register_latex_printer("Intersection", function(e) paste(.latex_setop_args(e), collapse = " \\cap "))
register_latex_printer("Complement", function(e) paste(.latex_setop_args(e), collapse = " \\setminus "))
register_str_printer("EmptySet", function(e) "EmptySet")
# SymPy Set._infimum_key: the numeric infimum when comparable, else oo.
.set_inf_dbl <- function(s) {
  v <- switch(s$type,
    Interval = { a <- s$args[[1L]]; if (is_const(a, "-oo")) -Inf else if (length(free_symbols(a))) NA_real_ else tryCatch(num_value_dbl(N(a)), error = function(e) NA_real_) },
    FiniteSet = {
      vals <- vapply(s$args, function(a) if (length(free_symbols(a))) NA_real_ else
        tryCatch({ z <- .neval(a, list()); if (is.complex(z) && Im(z) != 0) NA_real_ else as.double(Re(z)) }, error = function(e) NA_real_), 0)
      if (anyNA(vals)) NA_real_ else min(vals)
    },
    SpecialSet = switch(s$name, Reals = , Integers = , Rationals = -Inf, Naturals = 1, Naturals0 = 0, NA_real_),
    Range = { a <- .range_inf(s); if (is_const(a, "-oo")) -Inf else num_value_dbl(N(a)) },
    Union = min(vapply(s$args, .set_inf_dbl, 0)),
    NA_real_)
  if (is.na(v)) Inf else v
}
# Union and Intersection arguments are ordered(args, Set._infimum_key):
# infimum, then node count and default_sort_key. LaTeX prints them in that
# order; the str printer sorts them by default_sort_key, as LatticeOp does.
.set_print_order <- function(args) .key_sort(lapply(args, function(a) list(.set_inf_dbl(a), .nodes(a), .sort_key(a))), args)
register_str_printer("Union", function(e) paste0("Union(", paste(vapply(.key_sort(lapply(e$args, .sort_key), e$args), .str, ""), collapse = ", "), ")"))
register_str_printer("Intersection", function(e) paste0("Intersection(", paste(vapply(.key_sort(lapply(e$args, .sort_key), e$args), .str, ""), collapse = ", "), ")"))
register_str_printer("Complement", function(e) paste0("Complement(", .str(e$args[[1L]]), ", ", .str(e$args[[2L]]), ")"))

register_rebuilder("Interval", function(e, args) Interval(args[[1L]], args[[2L]], .is_true(args[[3L]]), .is_true(args[[4L]])))
register_rebuilder("FiniteSet", function(e, args) FiniteSet(args))
register_rebuilder("Union", function(e, args) do.call(Union, args))
register_rebuilder("Intersection", function(e, args) do.call(Intersection, args))
register_rebuilder("Complement", function(e, args) Complement(args[[1L]], args[[2L]]))






.max_sym <- function(a, b) { c <- .cmp(a, b); if (is.na(c)) NULL else if (c >= 0) a else b }
.min_sym <- function(a, b) { c <- .cmp(a, b); if (is.na(c)) NULL else if (c <= 0) a else b }





#' @rdname sets
#' @export
ProductSet <- function(...) {
  sets <- lapply(list(...), as_sym)
  # ProductSet.__new__: the nullary product is {()}, not the empty set, and
  # a product of finite sets stays unevaluated
  if (!length(sets)) return(FiniteSet(Tuple()))
  if (any(vapply(sets, function(s) s$type == "EmptySet", TRUE))) return(S_EmptySet)
  new_sym("ProductSet", sets)
}

register_str_printer("ProductSet", function(e) paste0("ProductSet(", paste(vapply(e$args, .str, ""), collapse = ", "), ")"))


# Union._measure: inclusion-exclusion over the intersections with nonzero measure
.union_measure <- function(args) {
  keys <- vapply(args, .k, "")
  sets <- lapply(seq_along(args), function(i) list(sos = i, inter = args[[i]]))
  total <- S_Zero; parity <- 1L
  while (length(sets)) {
    total <- .add(list(total, .mul(list(.int(parity), .add(lapply(sets, function(p) measure(p$inter)))))))
    nxt <- list(); seen <- character()
    for (p in sets) for (j in seq_along(args)) if (!(keys[j] %in% keys[p$sos])) {
      sos <- sort(c(p$sos, j))
      id <- paste(sos, collapse = ",")
      inter <- Intersection(args[[j]], p$inter)
      if (id %in% seen || isTRUE(is_zero_num(measure(inter)))) next
      seen <- c(seen, id); nxt[[length(nxt) + 1L]] <- list(sos = sos, inter = inter)
    }
    sets <- nxt; parity <- -parity
  }
  total
}





# ---------------------------------------------------------------------------
# ImageSet and ConditionSet

#' Lambda functions
#'
#' `Lambda(vars, expr)` is an anonymous symbolic function; apply it with
#' `apply_lambda()` or with `(...)` via `as.function()`.
#'
#' @param vars A symbol or list of symbols.
#' @param expr The body.
#' @param f A Lambda.
#' @param ... Arguments.
#' @return A `sym` Lambda.
#' @export
Lambda <- function(vars, expr) {
  vars <- .as_gens(list(vars))
  new_sym("Lambda", list(do.call(Tuple, vars), as_sym(expr)))
}

#' @rdname Lambda
#' @export
apply_lambda <- function(f, ...) {
  vals <- lapply(list(...), as_sym)
  vars <- f$args[[1L]]$args
  subs(f$args[[2L]], lapply(seq_along(vars), function(i) list(vars[[i]], vals[[i]])), simultaneous = TRUE)
}

register_bound_vars("Lambda", function(e) vapply(e$args[[1L]]$args, .k, ""))
register_latex_printer("Lambda", function(e) {
  v <- e$args[[1L]]$args
  vs <- if (length(v) == 1L) .latex(v[[1L]]) else paste0("\\left( ", paste(vapply(v, .latex, ""), collapse = ", \\  "), "\\right)")
  paste0("\\left( ", vs, " \\mapsto ", .latex(e$args[[2L]]), " \\right)")
})
register_str_printer("Lambda", function(e) {
  v <- e$args[[1L]]$args
  vs <- if (length(v) == 1L) .str(v[[1L]]) else paste0("(", paste(vapply(v, .str, ""), collapse = ", "), ")")
  paste0("Lambda(", vs, ", ", .str(e$args[[2L]]), ")")
})

register_str_printer("ImageSet", function(e) paste0("ImageSet(", .str(e$args[[1L]]), ", ", .str(e$args[[2L]]), ")"))

#' @rdname sets
#' @param sym A symbol.
#' @param condition A boolean condition.
#' @export
ConditionSet <- function(sym, condition, base = UniversalSet) {
  new_sym("ConditionSet", list(as_sym(sym), as_sym(condition), as_sym(base)))
}
register_bound_vars("ConditionSet", function(e) e$args[[1L]]$key)
register_str_printer("Contains", function(e) paste0("Contains(", .str(e$args[[1L]]), ", ", .str(e$args[[2L]]), ")"))

# ---------------------------------------------------------------------------
# Inequalities and solveset

#' Solve univariate inequalities
#'
#' Solves polynomial and rational inequalities in one real variable by
#' locating the real roots of numerator and denominator and testing the sign
#' on each piece, like SymPy's `solve_univariate_inequality()`. Returns a
#' relational (default) or a set with `relational = FALSE`.
#'
#' @param expr An inequality such as `x^2 > 4`.
#' @param x The variable.
#' @param relational Return a relational instead of a set.
#' @param domain The domain, by default the real line.
#' @return A boolean expression or a set.
#' @examples
#' x <- Symbol("x")
#' solve_univariate_inequality(x^2 > 4, x)
#' solve_univariate_inequality((x - 1) / (x + 2) <= 0, x, relational = FALSE)
#' @export
solve_univariate_inequality <- function(expr, x, relational = TRUE, domain = Reals) {
  expr <- as_sym(expr)
  x <- as_sym(x)
  if (expr$type == "BooleanAtom") return(if (relational) expr else if (.is_true(expr)) domain else S_EmptySet)
  op <- expr$type
  f <- .sub(expr$args[[1L]], expr$args[[2L]])
  nd <- fraction(together(f))
  num <- expand(nd[[1L]])
  den <- expand(nd[[2L]])
  crit <- list()
  zeros <- list()
  poles <- list()
  for (part in list(list(num, "z"), list(den, "p"))) {
    if (!.has1(part[[1L]], x)) next
    rs <- tryCatch(.solve_univariate(part[[1L]], x, check = FALSE), error = function(e) NULL)
    if (is.null(rs)) stop("could not solve the inequality")
    rs <- Filter(function(r) isTRUE(ask_fact(r, "real")), rs)
    rs <- Filter(function(r) .residual_ok(part[[1L]], x, r), rs)
    if (part[[2L]] == "z") zeros <- c(zeros, rs) else poles <- c(poles, rs)
  }
  # boundaries of the real domain: zeros of logarithm and even root arguments
  dom_pts <- list()
  walk_dom <- function(z) {
    if (z$type == "Function" && z$name == "log" && .has1(z$args[[1L]], x)) poles <<- c(poles, .real_zeros(z$args[[1L]], x))
    if (z$type == "Pow" && z$args[[2L]]$type == "Rational" && gmp::mod.bigz(gmp::denominator(z$args[[2L]]$value), 2) == 0 &&
        .has1(z$args[[1L]], x)) dom_pts <<- c(dom_pts, .real_zeros(z$args[[1L]], x))
    for (a in z$args) walk_dom(a)
  }
  walk_dom(f)
  crit <- .sort_set_elems(c(zeros, poles, dom_pts))
  crit <- crit[!duplicated(vapply(crit, .k, ""))]
  # test points are exact rationals when the bounds are, so close roots
  # keep their interval
  test_pt <- function(v) {
    val <- if (is_sym(v)) {
      ev <- tryCatch(subs(f, x, v), error = function(e) NULL)
      fe <- if (!is.null(ev) && !is_exact_number(ev) && !length(free_symbols(ev))) .facts(ev) else NULL
      if (!is.null(ev) && is_exact_number(ev)) num_sign(ev)
      else if (!is.null(fe) && isTRUE(fe[["positive"]])) 1
      else if (!is.null(fe) && isTRUE(fe[["negative"]])) -1
      else if (!is.null(fe) && isTRUE(fe[["zero"]])) 0
      else tryCatch(.neval(if (is.null(ev)) subs(f, x, Float(num_value_dbl(v))) else ev, list()), error = function(e) NA)
    } else tryCatch(.neval(subs(f, x, Float(v)), list()), error = function(e) NA)
    if (length(val) != 1L || is.na(val) || is.complex(val) && Im(val) != 0) return(NA)
    val <- Re(val)
    switch(op, StrictGreaterThan = val > 0, GreaterThan = val >= 0,
           StrictLessThan = val < 0, LessThan = val <= 0, Unequality = val != 0, Equality = val == 0)
  }
  exact_mid <- function(a, b) {
    ea <- is_exact_number(a); eb <- is_exact_number(b)
    if (is_inf_const(a) && is_inf_const(b)) return(S_Zero)
    if (is_inf_const(a) && eb) return(.sub(.sub(b, S_One), fcall("Abs", list(b))))
    if (is_inf_const(b) && ea) return(.add(list(a, S_One, fcall("Abs", list(a)))))
    if (ea && eb) return(.div(.add(list(a, b)), S_Two))
    NULL
  }
  cv <- vapply(crit, function(c) as.double(.neval(c, list())), 0)
  # periodic expressions are solved over one period, as SymPy does
  period <- .ineq_period(f, x)
  lo_sym <- S_NegativeInfinity
  hi_sym <- S_Infinity
  if (!is.null(period)) {
    pc <- .periodic_critical(f, x, period)
    crit <- pc$crit
    poles <- pc$poles
    zeros <- pc$zeros
    dom_pts <- c(dom_pts, list(S_Zero))
    cv <- vapply(crit, function(c) as.double(.neval(c, list())), 0)
    lo_sym <- S_Zero
    hi_sym <- period
  }
  pts <- c(as.double(.neval(lo_sym, list())), cv, as.double(.neval(hi_sym, list())))
  pieces <- list()
  for (k in seq_len(length(pts) - 1L)) {
    lo <- pts[k]; hi <- pts[k + 1L]
    a <- if (k == 1L) lo_sym else crit[[k - 1L]]
    b <- if (k == length(pts) - 1L) hi_sym else crit[[k]]
    mid <- exact_mid(a, b)
    if (is.null(mid))
      mid <- if (is.infinite(lo) && is.infinite(hi)) 0 else if (is.infinite(lo)) hi - 1 - abs(hi)
             else if (is.infinite(hi)) lo + 1 + abs(lo) else (lo + hi) / 2
    if (isTRUE(test_pt(mid))) pieces[[length(pieces) + 1L]] <- Interval(a, b, TRUE, TRUE)
  }
  nonstrict <- op %in% c("GreaterThan", "LessThan", "Equality")
  pk <- vapply(poles, .k, "")
  if (nonstrict) {
    for (z in zeros) if (!(z$key %in% pk)) pieces[[length(pieces) + 1L]] <- FiniteSet(z)
  }
  # points of the domain boundary where the inequality holds exactly
  zk <- vapply(zeros, .k, "")
  for (p in dom_pts) {
    if (p$key %in% pk || p$key %in% zk) next
    v <- tryCatch(subs(f, x, p), error = function(e) NULL)
    if (is.null(v) || !isTRUE(ask_fact(v, "real"))) next
    ok <- tryCatch(.rel_to_r(.make_rel(op, v, S_Zero)), error = function(e) NA)
    if (isTRUE(ok)) pieces[[length(pieces) + 1L]] <- FiniteSet(p)
  }
  # the domain itself when the inequality holds on the whole real line, as
  # SymPy (Interval(-oo, oo) is no longer turned into Reals)
  u <- do.call(Union, pieces)
  res <- if (.reals_like(u)) domain else Intersection(u, domain)
  if (!relational) return(res)
  as_relational(res, x)
}

# Period of an expression with trigonometric functions of linear arguments,
# or NULL when there are none. Errors if periods cannot be combined.
.ineq_period <- function(f, x) {
  tf <- Filter(function(z) z$name %in% c("sin", "cos", "tan", "cot", "sec", "csc") && .has1(z, x), atoms(f, "Function"))
  if (!length(tf)) return(NULL)
  ps <- lapply(tf, function(z) {
    la <- .linear_arg(z$args[[1L]], x)
    if (is.null(la)) stop("inequality solving needs linear arguments in trigonometric functions")
    base <- if (z$name %in% c("tan", "cot")) S_One else S_Two
    fcall("Abs", list(.div(.mul(list(base, S_Pi)), la$a)))
  })
  Reduce(function(a, b) {
    r <- simplify(.div(a, b))
    if (!is_exact_number(r)) stop("inequality solving needs commensurate periods")
    .mul(list(b, make_rational(gmp::as.bigq(gmp::numerator(r$value)))))
  }, ps)
}

# Critical points of a single-generator trigonometric expression in [0, T).
.periodic_critical <- function(f, x, period) {
  tf <- Filter(function(z) z$name %in% c("sin", "cos", "tan", "cot", "sec", "csc") && .has1(z, x), atoms(f, "Function"))
  if (length(unique(vapply(tf, function(z) z$args[[1L]]$key, ""))) != 1L || length(tf) != 1L)
    stop("inequality solving supports one trigonometric function of x")
  z <- tf[[1L]]
  la <- .linear_arg(z$args[[1L]], x)
  Tv <- as.double(.neval(period, list()))
  in_range <- function(v) v >= -1e-12 && v < Tv - 1e-12
  sols <- function(ang_base, step) {
    # x = (ang_base + k*step*pi - b)/a for integer k, inside [0, T)
    out <- list()
    for (k in -60:60) {
      v <- expand(.div(.sub(.add(list(ang_base, .mul(list(.int(k), .int(step), S_Pi)))), la$b), la$a))
      nv <- as.double(.neval(v, list()))
      if (in_range(nv)) out[[length(out) + 1L]] <- v
    }
    out
  }
  u <- Dummy("u")
  nd <- fraction(together(f))
  zeros <- list()
  if (.has1(nd[[1L]], x)) {
    pu <- subs(expand(nd[[1L]]), z, u)
    if (.has1(pu, x) || !is_polynomial(pu, u)) stop("inequality solving needs a polynomial in one trigonometric function")
    for (v in Filter(function(r) isTRUE(ask_fact(r, "real")), as.list(roots(pu, u)))) {
      zeros <- c(zeros, switch(z$name,
        sin = c(sols(fcall("asin", list(v)), 2L), sols(.sub(S_Pi, fcall("asin", list(v))), 2L)),
        cos = c(sols(fcall("acos", list(v)), 2L), sols(.neg(fcall("acos", list(v))), 2L)),
        tan = sols(fcall("atan", list(v)), 1L),
        cot = sols(fcall("acot", list(v)), 1L),
        list()))
    }
  }
  poles <- switch(z$name,
                  tan = , sec = sols(.mul(list(S_Half, S_Pi)), 1L),
                  cot = , csc = sols(S_Zero, 1L),
                  list())
  if (.has1(nd[[2L]], x)) {
    du <- subs(expand(nd[[2L]]), z, u)
    if (!.has1(du, x) && is_polynomial(du, u)) {
      for (v in Filter(function(r) isTRUE(ask_fact(r, "real")), as.list(roots(du, u)))) {
        poles <- c(poles, switch(z$name,
          sin = c(sols(fcall("asin", list(v)), 2L), sols(.sub(S_Pi, fcall("asin", list(v))), 2L)),
          cos = c(sols(fcall("acos", list(v)), 2L), sols(.neg(fcall("acos", list(v))), 2L)),
          list()))
      }
    }
  }
  all <- .sort_set_elems(c(zeros, poles))
  all <- all[!duplicated(vapply(all, .k, ""))]
  nz <- vapply(all, function(v) as.double(.neval(v, list())) > 1e-12, TRUE)
  list(crit = all[nz], zeros = zeros[!duplicated(vapply(zeros, .k, ""))], poles = poles)
}

.real_zeros <- function(g, x) {
  rs <- tryCatch(.solve_univariate(g, x, check = FALSE), error = function(e) list())
  Filter(function(r) isTRUE(ask_fact(r, "real")), rs)
}

#' Convert a set to a relational condition
#'
#' @param s A set.
#' @param x The variable.
#' @param ... For a `ProductSet()`, one more variable per further factor.
#' @return A boolean expression.
#' @export
as_relational <- function(s, x, ...) {
  s <- as_sym(s)
  # ProductSet.as_relational: one symbol per factor set
  if (s$type == "ProductSet") {
    xs <- lapply(c(list(x), list(...)), as_sym)
    if (length(xs) != length(s$args) || !all(vapply(xs, is_symbol_like, TRUE)))
      stop("number of symbols must match the number of sets", call. = FALSE)
    return(do.call(And, Map(as_relational, s$args, xs)))
  }
  x <- as_sym(x)
  switch(s$type,
         EmptySet = S_false,
         Interval = {
           i <- .iv(s)
           And(if (i$lo) Lt(i$a, x) else Le(i$a, x), if (i$ro) Lt(x, i$b) else Le(x, i$b))
         },
         FiniteSet = do.call(Or, lapply(s$args, function(p) Eq(x, p))),
         Range = {
           if (.range_empty(s)) return(S_false)
           p <- .range_parts(s)
           a <- if (is_const(p$start, "-oo")) .range_last(s) else p$start
           n <- .range_size(s)
           if (identical(n$key, S_One$key)) return(Eq(x, a))
           lo <- .range_inf(s); hi <- .range_sup(s)
           conds <- list(Eq(fcall("Mod", list(.sub(x, a), .zabs(.Z(p$step)) |> .zs())), S_Zero))
           if (!is_inf_const(lo)) conds <- c(conds, list(Ge(x, lo)))
           if (!is_inf_const(hi)) conds <- c(conds, list(Le(x, hi)))
           do.call(And, conds)
         },
         Union = {
           # Union.as_relational: two intervals touching at an excluded point
           # give And(Ne(x, p), x < sup, x > inf)
           if (length(s$args) == 2L && all(vapply(s$args, function(i) i$type == "Interval", NA))) {
             A <- .iv(s$args[[1L]]); B <- .iv(s$args[[2L]])
             if (identical(A$b$key, B$a$key) && A$ro && B$lo)
               return(And(Ne(x, A$b), Lt(x, B$b), Gt(x, A$a)))
           }
           do.call(Or, lapply(s$args, as_relational, x = x))
         },
         SpecialSet = if (s$name == "Reals") Lt(S_NegativeInfinity, x) & Lt(x, S_Infinity) else .contains(s, x),
         .contains(s, x))
}

#' Solve equations as sets
#'
#' `solveset(f, x, domain)` returns the solution set of `f = 0` (or of an
#' inequality) in `domain` (`Complexes` by default), like SymPy's
#' `solveset()`. Periodic trigonometric equations over the reals give
#' `ImageSet` unions; unsolvable equations give a `ConditionSet`.
#'
#' @param f An equation, expression or inequality.
#' @param x The variable.
#' @param domain `Complexes`, `Reals` or another set.
#' @return A set.
#' @examples
#' x <- Symbol("x")
#' solveset(x^2 - 1, x, Reals)
#' solveset(x^2 > 1, x, Reals)
#' solveset(sin(x), x, Reals)
#' @export
solveset <- function(f, x = NULL, domain = Complexes) {
  f <- as_sym(f)
  if (is.null(x)) x <- free_symbols(f)[[1L]]
  x <- as_sym(x)
  domain <- as_sym(domain)
  if (is_relational(f) && !(f$type %in% c("Equality", "Unequality"))) {
    return(solve_univariate_inequality(f, x, relational = FALSE,
                                       domain = if (domain$type == "SpecialSet" && domain$name == "Complexes") Reals else domain))
  }
  if (f$type == "Unequality") return(Complement(domain, solveset(Eq(f$args[[1L]], f$args[[2L]]), x, domain)))
  e <- .eq_to_expr(f)
  if (!.has1(e, x)) {
    z <- simplify(e)
    if (is_zero_num(z)) return(domain)
    if (isTRUE(ask_fact(z, "nonzero"))) return(S_EmptySet)
    return(ConditionSet(x, Eq(e, 0), domain))
  }
  is_real_dom <- domain$type == "SpecialSet" && domain$name == "Reals" || domain$type == "Interval"
  if (is_real_dom) {
    per <- .periodic_solutions(e, x)
    # the families are real already
    if (!is.null(per)) return(if (.reals_like(domain)) per else Intersection(per, domain))
  } else if (domain$type == "SpecialSet" && domain$name == "Complexes") {
    # SymPy returns the periodic families over the complexes too
    per <- .periodic_solutions(e, x, real = FALSE)
    if (is.null(per)) per <- .exp_solutions(e, x)
    if (!is.null(per)) return(per)
  }
  sols <- tryCatch(.solve_univariate(e, x), error = function(err) NULL)
  if (is.null(sols)) return(ConditionSet(x, Eq(e, 0), domain))
  if (is_real_dom) sols <- Filter(function(s) !identical(ask_fact(s, "real"), FALSE), sols)
  if (domain$type == "SpecialSet" && domain$name == "Complexes") return(if (length(sols)) FiniteSet(sols) else S_EmptySet)
  Intersection(FiniteSet(sols), domain)
}

# General solutions of sin(g(x)) = c style equations over the reals.
# exp(a*x + b) = v over the complexes: x = (log(v) + 2*pi*I*n - b)/a.
.exp_solutions <- function(e, x) {
  fs <- Filter(function(f) f$name == "exp" && .has1(f, x), atoms(e, "Function"))
  if (length(fs) != 1L) return(NULL)
  f <- fs[[1L]]
  if (!is_polynomial(f$args[[1L]], x) || length(.coeff_list(f$args[[1L]], x)) != 2L) return(NULL)
  u <- Dummy("u")
  pu <- subs(e, f, u)
  if (.has1(pu, x) || !is_polynomial(pu, u)) return(NULL)
  vals <- Filter(function(v) !is_zero_num(v), roots(pu, u))
  if (!length(vals)) return(S_EmptySet)
  cl <- .coeff_list(f$args[[1L]], x)
  a <- cl[[1L]]; b <- cl[[2L]]
  n <- Dummy("n", integer = TRUE)
  sets <- lapply(vals, function(v) {
    g <- .add(list(.mul(list(S_Two, S_Pi, S_I, n)), fcall("log", list(v))))
    ImageSet(Lambda(n, expand(.div(.sub(g, b), a))), Integers)
  })
  do.call(Union, sets[!duplicated(vapply(sets, .k, ""))])
}

.periodic_solutions <- function(e, x, real = TRUE) {
  fs <- Filter(function(f) f$name %in% c("sin", "cos", "tan") && .has1(f, x), atoms(e, "Function"))
  if (length(fs) != 1L) return(NULL)
  f <- fs[[1L]]
  if (!is_polynomial(f$args[[1L]], x) || length(.coeff_list(f$args[[1L]], x)) != 2L) return(NULL)
  u <- Dummy("u")
  pu <- subs(e, f, u)
  if (.has1(pu, x)) return(NULL)
  if (!is_polynomial(pu, u)) {
    # rational in u: numerator families minus the poles of the denominator
    nd <- fraction(together(pu))
    if (!is_polynomial(nd[[1L]], u) || !is_polynomial(nd[[2L]], u) || !.has1(nd[[2L]], u)) return(NULL)
    num <- .periodic_solutions(subs(nd[[1L]], u, f), x, real)
    den <- .periodic_solutions(subs(nd[[2L]], u, f), x, real)
    if (is.null(num) || is.null(den)) return(NULL)
    return(if (den$type == "EmptySet") num else Complement(num, den))
  }
  vals <- roots(pu, u, rootof = TRUE)
  # every root is needed, else the families would be incomplete
  if (sum(attr(vals, "multiplicity")) < degree(expand(pu), u)) return(NULL)
  n <- Dummy("n", integer = TRUE)
  cl <- .coeff_list(f$args[[1L]], x)
  a <- cl[[1L]]; b <- cl[[2L]]
  # each root gives its own families, wrapped in that root's condition only
  out <- list()
  for (v in vals) {
    cond <- NULL
    if (real) {
      isr <- ask_fact(v, "real")
      if (isFALSE(isr)) next
      if (f$name %in% c("sin", "cos")) {
        # exact range check; an undecidable one becomes a condition
        inr <- And(Ge(v, S_NegativeOne), Le(v, S_One))
        if (.is_false(inr)) next
        if (!.is_true(inr)) cond <- inr
      } else if (!isTRUE(isr)) cond <- And(Gt(v, S_NegativeInfinity), Lt(v, S_Infinity))
    }
    bs <- switch(f$name,
      sin = list(list(fcall("asin", list(v)), 2), list(.sub(S_Pi, fcall("asin", list(v))), 2)),
      cos = list(list(fcall("acos", list(v)), 2), list(.neg(fcall("acos", list(v))), 2)),
      tan = list(list(fcall("atan", list(v)), 1)))
    sets <- lapply(bs, function(z) {
      theta <- .add(list(.mul(list(.int(z[[2L]]), S_Pi, n)), z[[1L]]))
      ImageSet(Lambda(n, expand(.div(.sub(theta, b), a))), Integers)
    })
    sets <- sets[!duplicated(vapply(sets, .k, ""))]
    fam <- do.call(Union, sets)
    out[[length(out) + 1L]] <- if (is.null(cond)) fam else ConditionSet(x, cond, fam)
  }
  if (!length(out)) return(S_EmptySet)
  do.call(Union, out)
}
