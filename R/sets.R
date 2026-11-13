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
#' @param ... Elements or sets; for `ImageSet()`, further sets, one for each
#'   argument of a `Lambda` with several arguments.
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
  els <- els[!duplicated(vapply(els, .canon_key, ""))]
  new_sym("FiniteSet", .sort_set_elems(els))
}

# Lambda equality uses canonical_variables: elements that differ only in
# the bound symbols of a Lambda are one element of a set
.canon_key <- function(e) {
  if (!grepl("Lambda", e$key, fixed = TRUE)) return(e$key)
  lams <- Filter(function(n) is_sym(n) && n$type == "Lambda", preorder_traversal(e))
  vs <- unique(do.call(c, lapply(lams, function(l) .lambda_vars(l$args[[1L]]$args))))
  if (!length(vs)) return(e$key)
  reps <- lapply(seq_along(vs), function(i) list(vs[[i]], Symbol(paste0("_lambda_var", i))))
  .k(xreplace(e, reps))
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
  # a node may carry its own class key (the reduced Planck constant, an Atom)
  if (!is.null(e$extra$class_key)) return(e$extra$class_key)
  # the number sets sort by their SymPy class name (Basic.class_key)
  if (t == "SpecialSet") return(list(5, 0, e$name))
  if (t %in% .NUMBER_TYPES || (t == "Constant" && e$name %in% c("oo", "-oo"))) return(list(1, 0, "Number"))
  if (t == "Constant") return(list(2, 0, if (is.na(.CONST_CLASS[e$name])) e$name else .CONST_CLASS[[e$name]]))
  # a random symbol of sympy.stats is an Expr with arguments, not an Atom
  if (t == "Symbol" && !is.null(e$extra$rv_class)) return(list(5, 0, e$extra$sympy_class))
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
  # MatMul and MatAdd subclass Mul and Add
  if (t == "MatMul") return(list(3, 0, "MatMul"))
  if (t == "MatAdd") return(list(3, 1, "MatAdd"))
  if (t == "Function" && e$name %in% .BASIC_KEY_FUNS) return(list(5, 0, e$name))
  if (t == "Function") {
    i <- .FUNC_CLASS_INDEX[e$name]
    return(list(4, if (!is.na(i)) i[[1L]] else if (.is_undefined_application(e) || e$name %in% .NARGS0_FUNS) 0 else 10000, e$name))
  }
  # Applications that are separate node types here
  if (t %in% .NARGS0_TYPES) return(list(4, 0, t))
  if (!is.null(.qnode_def(t))) return(list(if (isTRUE(.qfind(t, "is_Atom")$val)) 2 else 5, 0, .qnode_def(t)$pyname))
  list(5, 0, if (!is.null(e$extra$sympy_class)) e$extra$sympy_class else t)
}

.ONE_SORT_KEY <- list(list(1, 0, "Number"), list(0, list()), list(), 1)

.sort_key <- function(e) {
  if (e$type == "Float" && is.nan(e$value)) e <- S_NaN
  if (e$type %in% .NUMBER_TYPES || (e$type == "Constant" && e$name %in% c("oo", "-oo")))
    return(list(.class_key(e), list(0, list()), list(),
                if (e$type == "Constant") (if (e$name == "oo") Inf else -Inf) else .sort_num(e)))
  if (e$type == "Constant" && e$name == "nan") return(list(list(1, 0, "NaN"), list(0, list()), list(), -Inf))
  # a random symbol sorts by its arguments, the symbol and the space
  if (e$type == "Symbol" && !is.null(e$extra$rv_class))
    return(list(.class_key(e), list(2, list(.sort_key(e$extra$symbol), .sort_key(e$extra$pspace))), .ONE_SORT_KEY, 1))
  # Dummy.sort_key: its name, then its dummy_index
  if (e$type == "Dummy" && !is.null(e$extra$id)) return(list(.class_key(e), list(2, list(e$name, e$extra$id)), .ONE_SORT_KEY, 1))
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
  args <- if (ex$type == "Symbol" && !is.null(ex$extra$rv_class)) list(.sort_key(ex$extra$symbol), .sort_key(ex$extra$pspace)) else
    if (!length(ex$args) || ex$type == "Indexed") list(.str(ex)) else {
    # Mul.as_ordered_factors: the commutative factors sorted, then the
    # noncommutative ones in their order
    a <- if (ex$type == "Add") .add_print_order(ex$args) else
      if (ex$type %in% c("And", "Or")) .bool_print_order(ex$args) else
      if (ex$type == "Mul") .mul_print_order(ex$args) else ex$args
    lapply(a, .sort_key)
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
        tryCatch({ z <- .neval(a, list()); if (length(z) != 1L || (is.complex(z) && Im(z) != 0)) NA_real_ else as.double(Re(z)) }, error = function(e) NA_real_), 0)
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
  # a distribution or finite density is called like a function (sympy.stats)
  if (is.function(f)) return(f(...))
  if (is_sym(f) && f$type != "Lambda") return(.rv_call(f, list(...)))
  vals <- lapply(list(...), as_sym)
  vars <- f$args[[1L]]$args
  # Lambda._match_signature: a tuple in the signature takes a tuple apart
  pairs <- list()
  match_sig <- function(vs, xs) for (i in seq_along(vs)) {
    if (vs[[i]]$type == "Tuple" && is_sym(xs[[i]]) && xs[[i]]$type == "Tuple") match_sig(vs[[i]]$args, xs[[i]]$args)
    else pairs[[length(pairs) + 1L]] <<- list(vs[[i]], xs[[i]])
  }
  match_sig(vars, vals)
  subs(f$args[[2L]], pairs, simultaneous = TRUE)
}

.lambda_vars <- function(l) do.call(c, lapply(l, function(v) if (v$type == "Tuple") .lambda_vars(v$args) else list(v)))
register_bound_vars("Lambda", function(e) vapply(.lambda_vars(e$args[[1L]]$args), .k, ""))
register_latex_printer("Lambda", function(e) {
  v <- e$args[[1L]]$args
  vs <- if (length(v) == 1L) .latex(v[[1L]]) else paste0("\\left( ", paste(vapply(v, .latex, ""), collapse = ", \\  "), "\\right)")
  paste0("\\left( ", vs, " \\mapsto ", .latex(e$args[[2L]]), " \\right)")
})
register_str_printer("Lambda", function(e) {
  v <- e$args[[1L]]$args
  # StrPrinter._print_Lambda: a one symbol signature prints as the symbol
  vs <- if (length(v) == 1L && v[[1L]]$type != "Tuple") .str(v[[1L]])
        else paste0("(", paste(vapply(v, .str, ""), collapse = ", "), if (length(v) == 1L) ",", ")")
  paste0("Lambda(", vs, ", ", .str(e$args[[2L]]), ")")
})

register_str_printer("ImageSet", function(e) paste0("ImageSet(", paste(vapply(e$args, .str, ""), collapse = ", "), ")"))

#' @rdname sets
#' @param sym A symbol.
#' @param condition A boolean condition.
#' @export
ConditionSet <- function(sym, condition, base = UniversalSet) {
  # ConditionSet.__new__
  err <- function(cls, msg) stop(.sympy_condition(cls, msg))
  islist <- function(v) is.list(v) && !is_sym(v)
  # sympify turns a tuple of symbols into a Tuple and fails on a list
  tup <- function(v) if (islist(v)) do.call(Tuple, lapply(v, tup)) else as_sym(v)
  if (islist(condition))
    err("SympifyError", paste0("SympifyError: [", paste(vapply(condition, function(b) .str(as_sym(b)), ""), collapse = ", "), "]"))
  if (islist(base))
    err("SympifyError", paste0("SympifyError: [", paste(vapply(base, function(b) .str(as_sym(b)), ""), collapse = ", "), "]"))
  sym <- tup(sym); condition <- as_sym(condition); base <- as_sym(base)
  flat <- if (sym$type == "Tuple") .lambda_vars(sym$args) else list(sym)
  if (anyDuplicated(vapply(flat, .k, ""))) err("BadSignatureError", "Duplicate symbols detected")
  if (!is_set(base)) err("TypeError", paste0("base set should be a Set object, not ", .str(base)))
  if (condition$type == "FiniteSet") condition <- do.call(And, lapply(condition$args, function(l) Eq(l, 0)))
  # as_Boolean: an expression other than a symbol is not a condition
  if (.is_expr_s(condition) && !(condition$type %in% c("Symbol", "Dummy", "Wild", .REL_TYPES, .BOOL_TYPES)))
    err("TypeError", paste0("expecting bool or Boolean, not `", .str(condition), "`."))
  if (identical(condition$key, S_true$key)) return(base)
  if (identical(condition$key, S_false$key)) return(S_EmptySet)
  if (base$type == "EmptySet") return(S_EmptySet)
  for (i in flat)
    # _diff_wrt: symbols, applied undefined functions, Indexed, Derivative
    if (!(is_symbol_like(i) || .is_undefined_application(i) || i$type %in% c("Indexed", "Derivative")))
      err("ValueError", paste0("`", .str(i), "` is not symbol-like"))
  inb <- tryCatch(.contains(base, sym), error = function(e) NULL)
  if (!is.null(inb) && identical(inb$key, S_false$key))
    err("TypeError", paste0("sym `", .str(sym), "` is not in base_set `", .str(base), "`"))
  know <- NULL
  if (base$type == "FiniteSet") {
    v <- vapply(base$args, function(el) .fb(subs(condition, sym, el)), NA)
    if (!any(is.na(v))) return(FiniteSet(base$args[v]))
    know <- FiniteSet(base$args[!is.na(v) & v])
    base <- FiniteSet(base$args[is.na(v)])
  }
  if (base$type == "ConditionSet") {
    s <- base$args[[1L]]; c <- base$args[[2L]]; b <- base$args[[3L]]
    sflat <- if (s$type == "Tuple") s$args else list(s)
    if (length(sflat) != length(flat)) err("BadSignatureError", "sym does not match sym of base set")
    fk <- function(e) vapply(free_symbols(e), .k, "")
    if (identical(s$key, sym$key)) {
      condition <- And(condition, c)
      base <- b
    } else if (!length(intersect(fk(c), fk(sym)))) {
      condition <- And(condition, xreplace(c, Map(list, sflat, flat)))
      base <- b
    } else if (!length(intersect(fk(condition), fk(s)))) {
      sym <- xreplace(sym, Map(list, flat, sflat))
      condition <- And(xreplace(condition, Map(list, flat, sflat)), c)
      base <- b
    }
  }
  if (condition$type == "Contains" && identical(condition$args[[1L]]$key, sym$key) && is_set(condition$args[[2L]]))
    return(Intersection(condition$args[[2L]], base))
  rv <- new_sym("ConditionSet", list(sym, condition, base))
  if (is.null(know)) rv else Union(know, rv)
}
# xreplace and subs rebuild a ConditionSet through ConditionSet.__new__
register_rebuilder("ConditionSet", function(e, args) ConditionSet(args[[1L]], args[[2L]], args[[3L]]))
register_bound_vars("ConditionSet", function(e) {
  s <- e$args[[1L]]
  if (s$type == "Tuple") vapply(s$args, .k, "") else s$key
})
register_str_printer("Contains", function(e) paste0("Contains(", .str(e$args[[1L]]), ", ", .str(e$args[[2L]]), ")"))

# ---------------------------------------------------------------------------
# Inequalities and solveset

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
