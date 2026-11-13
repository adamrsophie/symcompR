# solve() after sympy/solvers/solvers.py: solve, _solve, _solve_system,
# solve_linear, checksol, _invert, _tsolve and their helpers, keeping
# SymPy's control flow so that solutions, their order and their form match.

.fl <- function(flags, nm, default = NULL) {
  v <- flags[[nm]]
  if (is.null(v)) default else v
}

.sp_stop <- function(cls, msg) stop(.sympy_condition(cls, msg))

# sympy.utilities.misc.filldedent: a newline, then the text wrapped at 70
.filldedent <- function(s) paste0("\n", paste(strwrap(s, width = 71), collapse = "\n"))

# ---------------------------------------------------------------------------
# Solution dictionaries: ordered (key, value) pairs keyed by sym keys.

.sd_new <- function() list(k = list(), v = list())
.sd1 <- function(k, v) list(k = list(k), v = list(v))
.sd_idx <- function(d, k) if (!length(d$k)) NA_integer_ else match(k$key, vapply(d$k, .k, ""))
.sd_get <- function(d, k) {
  i <- .sd_idx(d, k)
  if (is.na(i)) NULL else d$v[[i]]
}
.sd_set <- function(d, k, v) {
  i <- .sd_idx(d, k)
  if (is.na(i)) {
    d$k[[length(d$k) + 1L]] <- k
    d$v[[length(d$v) + 1L]] <- v
  } else d$v[[i]] <- v
  d
}
.sd_del <- function(d, k) {
  i <- .sd_idx(d, k)
  if (!is.na(i)) {
    d$k <- d$k[-i]
    d$v <- d$v[-i]
  }
  d
}
.sd_map <- function(d, fn) {
  d$v <- lapply(d$v, fn)
  d
}
.sd_ident <- function(d) paste(sort(paste(vapply(d$k, .k, ""), vapply(d$v, .k, ""), sep = "=")), collapse = ";")
.sd_has_key <- function(d, k) !is.na(.sd_idx(d, k))

# Basic.subs with a dict: the larger keys first, ties by default_sort_key.
.sd_subs <- function(e, d) {
  if (!length(d$k)) return(e)
  o <- .key_sort(lapply(d$k, function(k) list(-.nodes(k), .sort_key(k))), seq_along(d$k))
  for (i in o) e <- subs(e, d$k[[i]], d$v[[i]])
  e
}

# _remove_duplicate_solutions
.sd_unique <- function(sols) {
  if (!length(sols)) return(sols)
  sols[!duplicated(vapply(sols, .sd_ident, ""))]
}

# default_sort_key of a Python dict
.sd_sort_key <- function(d) {
  items <- lapply(seq_along(d$k), function(i) list(.sort_key(d$k[[i]]), .sort_key(d$v[[i]])))
  list(length(items), .key_sort(items))
}

.sd_order_keys <- function(d) {
  if (length(d$k) < 2L) return(d)
  o <- .key_sort(lapply(d$k, function(a) list(.nodes(a), .sort_key(a))), seq_along(d$k))
  list(k = d$k[o], v = d$v[o])
}

# R form of a solution dict: a named list keyed by the str() of each key, or
# by the name of a Symbol or Dummy key, as the package's other code looks it up
.sd_to_r <- function(d) stats::setNames(d$v, vapply(d$k, function(k) if (k$type %in% c("Symbol", "Dummy")) k$name else .str(k), ""))

.sps_pyset_key <- function(xs) vapply(xs, .k, "")
.sp_ordered <- function(xs) if (length(xs) < 2L) xs else .key_sort(lapply(xs, function(a) list(.nodes(a), .sort_key(a))), xs)
.sp_sorted <- function(xs) if (length(xs) < 2L) xs else .key_sort(lapply(xs, .sort_key), xs)
.sp_uniq <- function(xs) if (!length(xs)) xs else xs[!duplicated(vapply(xs, .k, ""))]
.sp_in <- function(x, xs) length(xs) > 0L && x$key %in% vapply(xs, .k, "")

.has_free <- function(e, syms) {
  fs <- .free_sym_list(e)
  for (s in syms) {
    if (is_symbol_like(s)) {
      if (!is.null(fs[[s$key]])) return(TRUE)
    } else if (.contains_key(e, s$key)) return(TRUE)
  }
  FALSE
}

# Expr.as_independent(*syms, as_Add=)
.sp_indep <- function(e, syms, as_add = NULL) {
  if (is_zero_num(e)) return(list(e, e))
  if (is.null(as_add)) as_add <- e$type == "Add"
  want <- if (as_add) "Add" else "Mul"
  ident <- if (as_add) S_Zero else S_One
  if (e$type != want) {
    if (.has_free(e, syms)) return(list(ident, e))
    return(list(e, ident))
  }
  args <- .sympy_args(e)
  dep <- vapply(args, function(a) .has_free(a, syms), TRUE)
  if (as_add) return(list(.add(args[!dep]), .add(args[dep])))
  list(.mul(args[!dep]), .mul(args[dep]))
}

.sp_bool <- function(e) !(is_zero_num(e) || (e$type == "Float" && e$value == 0))

.sp_nan <- function(e) e$type == "Constant" && e$name == "nan"
.sp_illegal <- function(e) e$type == "Constant" && e$name %in% c("nan", "oo", "-oo", "zoo")
.sp_has_illegal <- function(e) {
  if (.sp_illegal(e)) return(TRUE)
  for (a in e$args) if (.sp_has_illegal(a)) return(TRUE)
  FALSE
}

.sp_is_number <- function(e) .expr_is_number(e)
.sp_fact <- function(e, f) {
  v <- tryCatch(ask_fact(e, f), error = function(err) NA)
  if (length(v) != 1L) NA else v
}

# ---------------------------------------------------------------------------
# Denominators

.sp_simple_dens <- function(f, symbols) {
  out <- list()
  for (d in do.call(denoms, c(list(f), symbols))) {
    if (d$type == "Pow" && is_number_type(d$args[[2L]])) {
      if (is_zero_num(d$args[[2L]])) next
      d <- d$args[[1L]]
    }
    out[[d$key]] <- d
  }
  unname(out)
}

.sp_reject_by_denominator <- function(d, sol, flags) {
  if (isTRUE(.sp_checksol(d, sol, flags))) return(TRUE)
  val <- .sd_subs(d, sol)
  .sp_is_number(val) && (isFALSE(.sp_fact(val, "finite")) || .sp_illegal(val))
}

.sp_reject_by_factor_finiteness <- function(fac, sol, flags) {
  if (is_polynomial(fac, sol$k)) return(FALSE)
  val <- .sd_subs(fac, sol)
  isFALSE(.sp_fact(val, "finite")) || .sp_illegal(val)
}

# ---------------------------------------------------------------------------
# checksol

.fuzzy_and_l <- function(vals) {
  rv <- TRUE
  for (v in vals) {
    if (isFALSE(v)) return(FALSE)
    if (is.na(v) || is.null(v)) rv <- NA
  }
  rv
}

# Expr.is_constant(*wrt, simplify=)
.sp_is_constant <- function(e, wrt = list(), simplify_ = TRUE) {
  if (.sp_is_number(e)) return(TRUE)
  free <- free_symbols(e)
  if (!length(free)) return(TRUE)
  if (length(wrt) && !any(vapply(wrt, function(w) .sp_in(w, free), TRUE))) return(TRUE)
  if (!length(wrt)) wrt <- free
  expr <- e
  if (simplify_) expr <- tryCatch(simplify(expr), error = function(err) expr)
  if (isTRUE(.sp_fact(expr, "zero"))) return(TRUE)
  wrt_number <- Filter(is_symbol_like, wrt)
  if (length(wrt_number) == length(free) && all(vapply(wrt_number, function(w) .sp_in(w, free), TRUE))) {
    num_at <- function(v) tryCatch({
      z <- subs(expr, lapply(free, function(s) list(s, v)), simultaneous = TRUE)
      if (.sp_nan(z)) NULL else z
    }, error = function(err) NULL)
    a <- num_at(S_Zero)
    if (!is.null(a)) {
      for (pt in list(S_One, Rational(-3L, 7L), .add(list(Rational(2L, 3L), .mul(list(Rational(5L, 11L), S_I)))))) {
        b <- num_at(pt)
        if (!is.null(b) && .sp_is_number(a) && .sp_is_number(b)) {
          va <- tryCatch(as.complex(.neval(a, list())), error = function(err) NA)
          vb <- tryCatch(as.complex(.neval(b, list())), error = function(err) NA)
          if (!is.na(va) && !is.na(vb) && base::Mod(va - vb) > 1e-10 * max(1, base::Mod(va))) return(FALSE)
        }
      }
    }
  }
  for (w in wrt_number) {
    deriv <- .diff1(expr, w)
    if (simplify_) deriv <- tryCatch(simplify(deriv), error = function(err) deriv)
    if (!is_zero_num(deriv)) return(FALSE)
  }
  cd <- .sp_fact(e, "zero")
  if (!is.na(cd)) return(cd)
  NA
}

.sp_numer <- function(e) .as_numer_denom(e)[[1L]]


.sp_checksol <- function(f, sol, flags = list()) {
  minimal <- isTRUE(.fl(flags, "minimal", FALSE))
  if (is.list(f) && !is_sym(f)) {
    if (!length(f)) .sp_stop("ValueError", "no functions to check")
    return(.fuzzy_and_l(lapply(f, function(fi) .sp_checksol(fi, sol, flags))))
  }
  if (inherits(f, "Poly")) f <- as_expr(f)
  f <- as_sym(f)
  if (.sp_is_number(f) && .is_expr_node(f)) return(.sp_fact(f, "zero"))
  if (f$type %in% c("Equality", "Unequality")) {
    if (f$args[[2L]]$type == "BooleanAtom") f <- sym_rebuild(f, rev(f$args))
    B <- f$args[[1L]]
    if (B$type == "BooleanAtom") {
      f <- .sd_subs(f, sol)
      if (f$type != "BooleanAtom") return(NA)
    } else if (f$type == "Equality") {
      f <- .sub(f$args[[1L]], f$args[[2L]])
    }
  }
  if (f$type == "BooleanAtom") return(.is_true(f))
  if (!is_relational(f) && !.sp_bool(f)) return(TRUE)
  if (any(vapply(sol$v, .sp_has_illegal, TRUE))) return(FALSE)
  numerical <- !isFALSE(.fl(flags, "numerical", TRUE))
  attempt <- -1L
  val <- NULL
  reps <- list()
  repeat {
    attempt <- attempt + 1L
    if (attempt == 0L) {
      val <- .sd_subs(f, sol)
      if (val$type == "Mul") {
        qs <- vapply(val$args, function(a) .has_type(a, "Quantity"), TRUE)
        if (any(qs)) val <- .mul(val$args[!qs])
      }
      if (.sp_has_illegal(val)) return(FALSE)
    } else if (attempt == 1L) {
      if (!.sp_is_number(val)) {
        ic <- .sp_is_constant(val, sol$k, simplify_ = !minimal)
        if (isFALSE(ic)) return(FALSE)
        val <- as_content_primitive(val)[[2L]]
        val <- .mexpand(.sp_numer(val))
      }
    } else if (attempt == 2L) {
      if (minimal) return(NA)
      if (!isFALSE(.fl(flags, "simplify", TRUE))) sol <- .sd_map(sol, function(v) simplify(v))
      val <- simplify(.sd_subs(f, sol))
      if (!isFALSE(.fl(flags, "force", TRUE))) {
        ps <- posify(val)
        val <- ps$expr
        reps <- ps$reps
        exval <- .mexpand(val)
        if (.sp_is_number(exval)) val <- exval
      }
    } else {
      saw <- FALSE
      walk <- function(p) {
        if (saw) return(invisible())
        if (p$type == "Pow" && p$args[[2L]]$type != "Integer") saw <<- TRUE
        else if (p$type == "Function") saw <<- TRUE
        if (!saw) for (a in p$args) walk(a)
      }
      walk(expand_mul(val))
      if (!saw) return(FALSE)
      if (!isFALSE(.fl(flags, "force", TRUE)) && length(reps)) val <- subs(val, reps)
      z <- .sp_fact(val, "zero")
      if (!is.na(z)) {
        if (.sp_is_number(val) && .has_fn(val, "LambertW")) return(NA)
        return(z)
      }
      break
    }
    if (is_exact_number(val)) return(is_zero_num(val))
    if (numerical && .sp_is_number(val)) {
      if (.sp_nan(val)) return(FALSE)
      v <- tryCatch(as.complex(.neval(val, list())), error = function(err) NA)
      if (is.na(v)) {
        vn <- tryCatch(N(val, 18), error = function(err) NULL)
        if (is.null(vn) || .sp_nan(vn)) return(FALSE)
        v <- tryCatch(as.complex(.neval(vn, list())), error = function(err) NA)
        if (is.na(v)) return(FALSE)
      }
      re <- Re(v); im <- Im(v)
      if (abs(re) < 1e-12) re <- 0
      if (abs(im) < 1e-12) im <- 0
      return(base::Mod(complex(real = re, imaginary = im)) < 1e-9)
    }
  }
  NA
}

.has_fn <- function(e, nm) {
  if (e$type == "Function" && e$name %in% nm) return(TRUE)
  for (a in e$args) if (.has_fn(a, nm)) return(TRUE)
  FALSE
}

# ---------------------------------------------------------------------------
# solve_linear

.sp_solve_linear <- function(lhs, rhs = S_Zero, symbols = list(), exclude = list()) {
  lhs <- as_sym(lhs)
  rhs <- as_sym(rhs)
  if (lhs$type == "Equality") {
    if (.sp_bool(rhs)) .sp_stop("ValueError", .filldedent(paste0("If lhs is an Equality, rhs must be 0 but was ", .sp_py_str(rhs))))
    rhs <- lhs$args[[2L]]
    lhs <- lhs$args[[1L]]
  }
  dens <- NULL
  eq <- if (is_zero_num(rhs)) lhs else .sub(lhs, rhs)
  nd <- .as_numer_denom(eq)
  n <- nd[[1L]]
  d <- nd[[2L]]
  if (!.sp_bool(n)) return(list(S_Zero, S_One))
  free <- free_symbols(n)
  if (!length(symbols)) {
    symbols <- free
  } else {
    bad <- Filter(function(s) !is_symbol_like(s), symbols)
    if (length(bad)) {
      eg <- if (length(symbols) == 1L) paste0("solve(", .str(eq), ", ", .str(symbols[[1L]]), ")") else
        paste0("solve(", .str(eq), ", *[", paste(vapply(symbols, .str, ""), collapse = ", "), "])")
      badstr <- if (length(bad) == 1L) .str(bad[[1L]]) else paste0("[", paste(vapply(bad, .str, ""), collapse = ", "), "]")
      .sp_stop("ValueError", .filldedent(paste0("solve_linear only handles symbols, not ", badstr,
                                    ". To isolate non-symbols use solve, e.g. >>> ", eg, " <<<.")))
    }
    symbols <- Filter(function(s) .sp_in(s, symbols), free)
  }
  symbols <- Filter(function(s) !.sp_in(s, exclude), symbols)
  if (!length(symbols)) return(list(S_Zero, S_One))
  derivs <- atoms(n, "Derivative")
  all_zero <- TRUE
  for (xi in .sp_sorted(symbols)) {
    dx <- Filter(function(der) .has1(der, xi), derivs)
    newn <- if (length(dx)) subs(n, lapply(dx, function(der) list(der, doit(der)))) else n
    dnewn_dxi <- .diff1(newn, xi)
    fr <- free_symbols(dnewn_dxi)
    same_set <- length(fr) == length(symbols) && all(vapply(fr, function(s) .sp_in(s, symbols), TRUE))
    if (.sp_bool(dnewn_dxi) && (!length(fr) || any(vapply(fr, function(s) .sp_bool(.diff1(dnewn_dxi, s)), TRUE)) || same_set)) {
      all_zero <- FALSE
      if (.sp_nan(dnewn_dxi)) break
      if (!.has1(dnewn_dxi, xi)) {
        vi <- .mul(list(S_NegativeOne, .inv(dnewn_dxi), subs(newn, xi, S_Zero)))
        if (is.null(dens)) dens <- .sp_simple_dens(eq, symbols)
        if (!any(vapply(dens, function(di) isTRUE(.sp_checksol(di, .sd1(xi, vi), list(minimal = TRUE))), TRUE))) {
          irep <- Filter(function(i) .sp_is_number(i$args[[1L]]), atoms(vi, "Integral"))
          if (length(irep)) vi <- subs(vi, lapply(irep, function(i) list(i, doit(i))))
          return(list(xi, expand_mul(vi)))
        }
      }
    }
  }
  if (all_zero) return(list(S_Zero, S_One))
  if (is_symbol_like(n)) return(list(S_Zero, S_Zero))
  list(n, d)
}

# ---------------------------------------------------------------------------
# _invert

.sp_ispow <- function(e) e$type == "Pow" || (e$type == "Function" && e$name == "exp")

.sp_inverse_fn <- function(nm) {
  inv <- c(exp = "log", log = "exp", tan = "atan", cot = "acot", asin = "sin", acos = "cos", atan = "tan",
           acot = "cot", asec = "sec", acsc = "csc", sinh = "asinh", tanh = "atanh", coth = "acoth",
           asinh = "sinh", acosh = "cosh", atanh = "tanh", acoth = "coth", asech = "sech", acsch = "csch",
           erf = "erfinv", erfc = "erfcinv", erfinv = "erf", erfcinv = "erfc", conjugate = "conjugate")
  if (!(nm %in% names(inv)) || is.na(inv[[nm]])) return(NULL)
  inv[[nm]]
}

.sp_invert <- function(eq, symbols, integer_power = FALSE) {
  eq <- as_sym(eq)
  # eq.func(*eq.args): a flat, evaluated Add or Mul
  if (eq$type == "Add") eq <- .add(eq$args) else if (eq$type == "Mul") eq <- .mul(eq$args)
  free <- free_symbols(eq)
  if (!length(symbols)) symbols <- free
  if (!any(vapply(symbols, function(s) .sp_in(s, free), TRUE))) return(list(eq, S_Zero))
  lhs <- eq
  rhs <- S_Zero
  repeat {
    was <- lhs
    repeat {
      id <- .sp_indep(lhs, symbols)
      indep <- id[[1L]]; dep <- id[[2L]]
      if (lhs$type == "Add") {
        if (is_zero_num(indep)) break
        lhs <- dep
        rhs <- .sub(rhs, indep)
      } else {
        if (is_one_num(indep)) break
        lhs <- dep
        rhs <- .div(rhs, indep)
      }
    }
    if (lhs$type == "Add") {
      tk <- character(); tv <- list(); td <- list()
      for (a in .sympy_args(lhs)) {
        id <- .sp_indep(a, symbols)
        k <- id[[2L]]$key
        if (k %in% tk) tv[[k]] <- c(tv[[k]], list(id[[1L]])) else {
          tk <- c(tk, k); tv[[k]] <- list(id[[1L]]); td[[k]] <- id[[2L]]
        }
      }
      if (any(lengths(tv) > 1L)) {
        lhs <- .add(lapply(tk, function(k) .mul(list(.add(tv[[k]]), td[[k]]))))
      }
    }
    if (lhs$type == "Add" && !.sp_bool(rhs) && length(lhs$args) == 2L && !is_polynomial(lhs, symbols)) {
      ab <- .sp_ordered(.sympy_args(lhs))
      a <- ab[[1L]]; b <- ab[[2L]]
      ia <- .sp_indep(a, symbols); ib <- .sp_indep(b, symbols)
      ai <- ia[[1L]]; ad <- ia[[2L]]; bi <- ib[[1L]]; bd <- ib[[2L]]
      if (.sp_ispow(ad) || .sp_ispow(bd)) {
        abe <- .as_base_exp(ad); bbe <- .as_base_exp(bd)
        if (identical(abe[[1L]]$key, bbe[[1L]]$key) && is.null(.extract_additively(abe[[2L]], bbe[[2L]]))) {
          lhs <- powsimp(powdenest(.div(ad, bd)))
          rhs <- .div(.neg(bi), ai)
        } else {
          rat <- .div(ad, bd)
          lhs2 <- powsimp(rat)
          if (!identical(lhs2$key, rat$key)) {
            lhs <- lhs2
            rhs <- .div(.neg(bi), ai)
          }
        }
      } else if (identical(ai$key, .neg(bi)$key)) {
        if (ad$type == "Function" && bd$type == "Function" && ad$name == bd$name) {
          if (length(ad$args) == 1L && length(bd$args) == 1L) lhs <- .sub(ad$args[[1L]], bd$args[[1L]])
          else if (length(ad$args) == length(bd$args)) .sp_stop("NotImplementedError", "equal function with more than 1 argument")
          else .sp_stop("ValueError", "function with different numbers of args")
        }
      }
    } else if (lhs$type == "Mul" && any(vapply(lhs$args, .sp_ispow, TRUE))) {
      lhs <- powsimp(powdenest(lhs))
    }
    if (lhs$type == "Function") {
      invn <- .sp_inverse_fn(lhs$name)
      if (!is.null(invn) && length(lhs$args) == 1L) {
        rhs <- fcall(invn, list(rhs))
        lhs <- lhs$args[[1L]]
      } else if (lhs$name == "atan2") {
        y <- lhs$args[[1L]]; x <- lhs$args[[2L]]
        lhs <- .mul(list(S_Two, fcall("atan", list(.div(y, .add(list(.sqrt(.add(list(.pow(x, S_Two), .pow(y, S_Two)))), x)))))))
      } else if (rhs$type == "Function" && rhs$name == lhs$name) {
        if (length(lhs$args) == 1L && length(rhs$args) == 1L) {
          lhs <- lhs$args[[1L]]
          rhs <- rhs$args[[1L]]
        } else if (length(lhs$args) == length(rhs$args)) .sp_stop("NotImplementedError", "equal function with more than 1 argument")
        else .sp_stop("ValueError", "function with different numbers of args")
      }
    }
    if (.sp_bool(rhs) && lhs$type == "Pow" && lhs$args[[2L]]$type == "Integer" && lhs$args[[2L]]$value < 0) {
      lhs <- .inv(lhs)
      rhs <- .inv(rhs)
    }
    if (lhs$type == "Pow") {
      ex <- lhs$args[[2L]]
      if ((ex$type == "Integer" && integer_power) ||
          (ex$type != "Integer" && length(symbols) > 1L &&
           sum(vapply(symbols, function(s) .sp_in(s, free_symbols(lhs$args[[1L]])), TRUE)) > 1L)) {
        rhs <- .pow(rhs, .inv(ex))
        lhs <- lhs$args[[1L]]
      }
    }
    if (identical(lhs$key, was$key)) break
  }
  list(rhs, lhs)
}


# Add(a, b, evaluate=False): kept unevaluated as SymPy does for Eq(L, R)
# Eq(L, R).rewrite(Add, evaluate=False): Add._from_args of the terms of L
# and of -R, without cancellation; one side zero gives the other side.
.sp_add_uneval <- function(L, R) {
  if (is_zero_num(L)) return(R)
  if (is_zero_num(R)) return(L)
  args <- c(.aargs(L), .aargs(.neg(R)))
  if (length(args) == 1L) args[[1L]] else new_sym("Add", args)
}

# Expr.extract_additively(c): NULL when it cannot be done
.extract_additively <- function(e, c) {
  if (is_zero_num(c)) return(e)
  if (identical(e$key, c$key)) return(S_Zero)
  if (is_zero_num(e)) return(NULL)
  if (is_number_type(e)) {
    if (!is_number_type(c)) return(NULL)
    co <- e
    diff <- .sub(co, c)
    if (isTRUE(.sp_fact(co, "positive")) && isTRUE(.sp_fact(diff, "positive")) ||
        isTRUE(.sp_fact(co, "negative")) && isTRUE(.sp_fact(diff, "negative"))) return(diff)
    return(NULL)
  }
  cc <- as_coeff_add_split(c)
  c0 <- cc[[1L]]
  xc <- cc[[2L]]
  ee <- as_coeff_add_split(e)
  s0 <- ee[[1L]]
  rest <- ee[[2L]]
  if (!is_zero_num(c0)) {
    xa <- .extract_additively(s0, c0)
    if (is.null(xa)) return(NULL)
    e <- .add(list(xa, rest))
    if (is_zero_num(xc)) return(e)
  }
  dd <- .sub(e, xc)
  if (.cops(dd) <= .cops(e)) {
    # only accept when the term of c actually appears in e
    ea <- vapply(.aargs(e), .k, "")
    if (all(vapply(.aargs(xc), function(t) t$key %in% ea, TRUE))) return(dd)
  }
  NULL
}

# ---------------------------------------------------------------------------
# solve() front end

.sp_iterable <- function(x) (is.list(x) && !is_sym(x) && !inherits(x, "Poly")) || (is_sym(x) && x$type == "Tuple") ||
  inherits(x, "symmatrix")

.sp_as_list <- function(w) {
  if (inherits(w, "symmatrix")) return(lapply(.mat_entries(w), as_sym))
  if (is_sym(w) && w$type == "Tuple") return(w$args)
  if (is.list(w) && !is_sym(w) && !inherits(w, "Poly")) return(lapply(w, function(z) if (inherits(z, "Poly") || inherits(z, "symmatrix")) z else as_sym(z)))
  list(if (inherits(w, "Poly") || inherits(w, "symmatrix")) w else as_sym(w))
}

.mat_entries <- function(M) unlist(.mrows(M), recursive = FALSE)

.sp_is_undef <- function(p) .is_undefined_application(p)

.sp_preorder <- function(e) {
  out <- list()
  walk <- function(x) {
    out[[length(out) + 1L]] <<- x
    for (a in .sympy_args(x)) walk(a)
  }
  walk(e)
  out
}

.sp_solve <- function(f, symbols, flags) {
  default <- TRUE
  for (k in c("cubics", "quartics", "quintics")) {
    v <- flags[[k]]
    default <- if (is.null(v)) default else isTRUE(as.logical(v))
    flags[[k]] <- default
  }
  implicit <- isTRUE(.fl(flags, "implicit", FALSE))
  as_set <- isTRUE(.fl(flags, "set", FALSE))
  bare_f <- !.sp_iterable(f)
  if (!is.null(flags$quick) && !isTRUE(flags$particular))
    .sp_stop("ValueError", "when using `quick`, `particular` should be True")
  if (isTRUE(flags$particular) && bare_f)
    .sp_stop("ValueError", .filldedent(paste("The 'particular/quick' flag is usually used with systems of equations.",
                                 "Either pass your equation in a list or consider using a solver like",
                                 "`diophantine` if you are looking for a solution in integers.")))
  f <- .sp_as_list(f)
  ordered_symbols <- NA
  if (!length(symbols)) {
    syms <- unique_syms(lapply(f, function(fi) if (is_sym(fi)) free_symbols(fi) else if (inherits(fi, "Poly")) free_symbols(as_expr(fi)) else list()))
    if (length(syms) < length(f)) {
      for (fi in f) {
        if (!is_sym(fi)) next
        stack <- list(fi)
        while (length(stack)) {
          p <- stack[[1L]]
          stack <- stack[-1L]
          if (.sp_is_undef(p)) {
            if (!as_set) flags$dict <- TRUE
            if (!.sp_in(p, syms)) syms[[length(syms) + 1L]] <- p
            next
          }
          stack <- c(.sympy_args(p), stack)
        }
      }
    }
    ordered_symbols <- FALSE
    symbols <- .sp_ordered(syms)
  } else {
    seq_given <- TRUE
    if (length(symbols) == 1L && .sp_iterable(symbols[[1L]])) {
      seq_given <- !isTRUE(attr(symbols[[1L]], "pyset"))
      symbols <- .sp_as_list(symbols[[1L]])
    }
    symbols <- lapply(symbols, as_sym)
    ordered_symbols <- length(symbols) > 0L && seq_given
    us <- .sp_uniq(symbols)
    if (length(us) != length(symbols)) {
      ordered_symbols <- FALSE
      symbols <- .sp_ordered(symbols)
    } else symbols <- us
  }
  if (length(symbols) != length(.sp_uniq(symbols))) .sp_stop("ValueError", "duplicate symbols given")
  exclude <- flags$exclude
  flags$exclude <- NULL
  if (length(exclude)) {
    ex <- .sp_as_list(exclude)
    exs <- unique_syms(lapply(ex, free_symbols))
    symbols <- Filter(function(s) !.sp_in(s, exs), symbols)
  }
  f <- Filter(function(s) !(is_sym(s) && .is_true(s)), f)
  i <- 0L
  while (i < length(f)) {
    i <- i + 1L
    fi <- f[[i]]
    if (is_sym(fi) && fi$type %in% c("Equality", "Unequality")) {
      L <- fi$args[[1L]]; R <- fi$args[[2L]]
      if (inherits(L, "symmatrix") || inherits(R, "symmatrix")) {
        fi <- .sub(L, R)
      } else {
        if (R$type == "BooleanAtom") { tmp <- L; L <- R; R <- tmp }
        if (L$type == "BooleanAtom") {
          if (fi$type == "Unequality") L <- Not(L)
          if (is_relational(R)) fi <- if (.is_false(L)) Not(R) else R
          else if (is_symbol_like(R)) return(L)
          else if (R$type %in% .BOOL_TYPES && is_symbol_like(Not(R))) return(Not(L))
          else .sp_stop("NotImplementedError", .filldedent("Unanticipated argument of Eq when other arg is True or False."))
        } else if (fi$type == "Equality") {
          fi <- .sp_add_uneval(L, R)
        }
      }
      f[[i]] <- fi
    }
    if (is_sym(fi) && is_relational(fi)) {
      if (length(symbols) != 1L) .sp_stop("ValueError", "can only solve for one symbol at a time")
      return(reduce_inequalities(f, symbols))
    }
    if (inherits(fi, "Poly")) f[[i]] <- as_expr(fi)
    if (inherits(f[[i]], "symmatrix")) {
      bare_f <- FALSE
      f <- c(f, .mat_entries(f[[i]]))
      f[[i]] <- S_Zero
      next
    }
    f[[i]] <- .map_bottom_up(f[[i]], function(w) {
      if (w$type == "Function" && w$name %in% c("sinh", "cosh", "tanh", "coth", "sech", "csch") && .has_free(w, symbols))
        rewrite(w, "exp") else w
    })
    freei <- free_symbols(f[[i]])
    if (length(freei) && all(vapply(freei, function(s) isTRUE(.sp_fact(s, "extended_real")) || isTRUE(.sp_fact(s, "imaginary")), TRUE))) {
      ri <- .as_real_imag(f[[i]])
      fr <- ri[[1L]]; fim <- ri[[2L]]
      had <- .reim_atoms(f[[i]])
      if (.sp_bool(fr) && .sp_bool(fim) && !identical(fr$key, fim$key) &&
          !any(vapply(list(fr, fim), function(z) length(setdiff(.reim_atoms(z), had)) > 0L, TRUE))) {
        bare_f <- FALSE
        f <- append(f[-i], list(fr, fim), after = i - 1L)
      }
    }
  }
  if (any(vapply(f, function(fi) is_sym(fi) && fi$type == "BooleanAtom", TRUE))) {
    if (as_set) return(.pytuple(list(list(), .sps_pyset(list()))))
    return(list())
  }
  for (i in seq_along(f)) {
    fi <- f[[i]]
    repeat {
      was <- fi
      fi <- .map_bottom_up(fi, function(ex) {
        if (ex$type == "Function" && ex$name %in% c("Abs", "Min", "Max") &&
            any(vapply(ex$args, function(a) .has_free(a, symbols), TRUE)))
          .sps_pw_rewrite(separatevars(ex)) else ex
      })
      if (identical(was$key, fi$key)) break
    }
    for (e in .sp_find(fi, function(ex) ex$type == "Function" && ex$name %in% c("Abs", "Min", "Max"))) {
      if (.has_free(e, symbols))
        .sp_stop("NotImplementedError", paste0("solving ", .sp_py_str(e), " when the argument is not real or imaginary."))
    }
    fi <- .map_bottom_up(fi, function(a) if (a$type == "Function" && a$name == "arg") rewrite(rewrite(a, "atan2"), "atan") else a)
    f[[i]] <- fi
  }
  freim <- Filter(function(fi) .has_fn(fi, c("re", "im")), f)
  if (length(freim)) {
    irf <- list()
    for (s in symbols) {
      if (isTRUE(.sp_fact(s, "real")) || isTRUE(.sp_fact(s, "imaginary"))) next
      rs <- fcall("re", list(s)); is <- fcall("im", list(s))
      if (any(vapply(freim, function(fi) .contains_key(fi, rs$key) || .contains_key(fi, is$key), TRUE)))
        irf[[length(irf) + 1L]] <- list(s, .add(list(rs, .mul(list(S_I, is)))))
    }
    if (length(irf)) {
      for (pr in irf) {
        f <- c(lapply(f, function(fi) xreplace(fi, list(pr))), list(.sub(pr[[1L]], pr[[2L]])))
        symbols <- c(symbols, list(fcall("re", list(pr[[1L]])), fcall("im", list(pr[[1L]]))))
      }
      bare_f <- FALSE
      flags$dict <- TRUE
    }
  }
  rc <- .sp_recast(f, symbols)
  f <- rc$f
  symbols <- rc$symbols
  swap_sym <- rc$restore
  newf <- list()
  for (fi in f) {
    ok <- FALSE
    if (.has_free(fi, symbols)) {
      ok <- TRUE
    } else if (.sp_is_number(fi)) {
      if (is_number_type(fi)) {
        if (!.sp_bool(fi)) next
        return(if (as_set) .pytuple(list(symbols, .sps_pyset(list()))) else list())
      }
      ok <- TRUE
    } else if (isTRUE(.sp_is_constant(fi))) ok <- TRUE
    if (ok) newf[[length(newf) + 1L]] <- fi
  }
  if (!length(newf)) return(if (as_set) .pytuple(list(symbols, .sps_pyset(list()))) else list())
  f <- newf
  # mask off non-invertible objects
  non_inverts <- list()
  seen <- character()
  for (fi in f) {
    walk <- function(p) {
      isexpr <- .is_expr_node(p) && p$type != "Piecewise"
      if (isexpr) {
        nonexempt <- p$type == "Function" && p$name %in% c("re", "im")
        if (!nonexempt && (!length(p$args) || .sp_in(p, symbols) || p$type %in% c("Add", "Mul") ||
                           (p$type == "Pow" && !implicit) || (p$type == "Function" && !implicit))) {
          for (a in .sympy_args(p)) walk(a)
          return(invisible())
        }
        if (!(p$key %in% seen)) {
          seen <<- c(seen, p$key)
          if (.has_free(p, symbols)) non_inverts[[length(non_inverts) + 1L]] <<- p
          else {
            for (a in .sympy_args(p)) walk(a)
            return(invisible())
          }
        }
      }
      invisible()
    }
    walk(fi)
  }
  ni_pairs <- lapply(non_inverts, function(p) list(p, Dummy()))
  if (length(ni_pairs)) f <- lapply(f, function(fi) subs(fi, ni_pairs))
  ni_restore <- lapply(ni_pairs, function(pr) list(pr[[2L]], subs(xreplace(pr[[1L]], swap_sym), swap_sym)))
  floats <- FALSE
  if (!isFALSE(.fl(flags, "rational", TRUE))) {
    for (i in seq_along(f)) {
      if (.has_type(f[[i]], "Float")) {
        floats <- TRUE
        f[[i]] <- nsimplify(f[[i]], rational = TRUE)
      }
    }
  }
  flags[["_denominators"]] <- .sp_simple_dens(f[[1L]], symbols)
  has_pw <- function(e) {
    if (e$type == "Piecewise") return(.has_free(e, symbols))
    any(vapply(e$args, has_pw, TRUE))
  }
  for (i in seq_along(f)) if (has_pw(f[[i]])) f[[i]] <- piecewise_fold(f[[i]])
  targs <- list()
  for (fi in f) targs <- c(targs, Filter(function(z) z$name %in% .TRIG_NAMES, atoms(fi, "Function")))
  targs <- .sp_uniq(targs)
  if (length(targs) > 1L) {
    isadd <- vapply(targs, function(t) t$args[[1L]]$type == "Add", TRUE)
    add <- Filter(function(t) .has_free(t, symbols), targs[isadd])
    other <- Filter(function(t) .has_free(t, symbols), targs[!isadd])
    trep <- list()
    for (t in add) {
      id <- .sp_indep(t$args[[1L]], symbols)
      ind <- id[[1L]]; dep <- id[[2L]]
      if (.sp_in(dep, symbols) || .sp_in(.neg(dep), symbols)) {
        n <- if (!is_number_type(ind)) Dummy() else ind
        trep[[length(trep) + 1L]] <- list(t, xreplace(TR10(fcall(t$name, list(.add(list(dep, n))))), list(list(n, ind))))
      }
    }
    if (length(other) && length(other) <= 2L) {
      base <- if (length(other) > 1L) .sps_gcd_terms(lapply(other, function(i) i$args[[1L]])) else other[[1L]]$args[[1L]]
      for (t in other) trep[[length(trep) + 1L]] <- list(t, TR11(t, base))
    }
    if (length(trep)) f <- lapply(f, function(fi) xreplace(fi, trep))
  }
  linear <- NA
  if (bare_f) {
    solution <- NULL
    if (length(symbols) != 1L) solution <- .sp_solve_undetermined(f[[1L]], symbols, flags)
    if (!length(solution)) solution <- .sp__solve(f[[1L]], symbols, flags)
  } else {
    ls <- .sp_solve_system(f, symbols, flags)
    linear <- ls[[1L]]
    solution <- ls[[2L]]
  }
  as_dict <- isTRUE(.fl(flags, "dict", FALSE))
  tuple_format <- function(s) lapply(s, function(d) .pytuple(lapply(symbols, function(x) { v <- .sd_get(d, x); if (is.null(v)) x else v })))
  unpack <- NULL
  if (as_dict || as_set) {
    unpack <- NULL
  } else if (bare_f) {
    if (length(symbols) == 1L) unpack <- function(s) lapply(s, function(d) .sd_get(d, symbols[[1L]]))
    else if (length(solution) == 1L && length(solution[[1L]]$k) == length(symbols)) unpack <- function(s) .sd_to_r(s[[1L]])
    else if (isTRUE(ordered_symbols)) unpack <- tuple_format
    else unpack <- function(s) lapply(s, .sd_to_r)
  } else {
    if (length(solution)) {
      if (isTRUE(linear) && length(solution) == 1L) unpack <- function(s) .sd_to_r(s[[1L]])
      else if (isTRUE(ordered_symbols)) unpack <- tuple_format
      else unpack <- function(s) lapply(s, .sd_to_r)
    }
  }
  if (length(ni_restore)) solution <- lapply(solution, function(s) .sd_map(s, function(v) subs(v, ni_restore)))
  if (length(swap_sym)) {
    symbols <- lapply(symbols, function(k) .swap_get(swap_sym, k))
    solution <- lapply(solution, function(s) {
      list(k = lapply(s$k, function(k) .swap_get(swap_sym, k)), v = lapply(s$v, function(v) subs(v, swap_sym)))
    })
  }
  check <- !isFALSE(.fl(flags, "check", TRUE))
  if (floats && length(solution) && is.null(flags$rational)) {
    solution <- lapply(solution, function(s) .sd_map(s, function(v) nfloat(v, exponent = FALSE)))
    solution <- .sd_unique(solution)
  }
  if (check && length(solution)) {
    solution <- Filter(function(sol) {
      v <- .fuzzy_and_l(lapply(seq_along(sol$k), function(j) .sp_check_assumptions(sol$v[[j]], sol$k[[j]])))
      !isFALSE(v)
    }, solution)
  }
  if (!length(solution)) return(if (as_set) .pytuple(list(symbols, .sps_pyset(list()))) else list())
  if (!as_set) {
    solution <- lapply(solution, .sd_order_keys)
    solution <- .key_sort(lapply(solution, .sd_sort_key), solution)
  }
  if (!(as_set || as_dict)) return(unpack(solution))
  if (as_dict) return(if (isTRUE(flags[["_raw"]])) solution else lapply(solution, .sd_to_r))
  k <- if (isTRUE(ordered_symbols)) symbols else .sp_ordered(.sp_uniq(do.call(c, lapply(solution, function(s) s$k))))
  tuples <- lapply(solution, function(s) .pytuple(lapply(k, function(ki) { v <- .sd_get(s, ki); if (is.null(v)) ki else v })))
  keys <- vapply(tuples, function(t) paste(vapply(t, .k, ""), collapse = "|"), "")
  tuples <- tuples[!duplicated(keys)]
  .pytuple(list(k, .sps_pyset(tuples)))
}

.sps_pyset <- function(x) structure(x, pyset = TRUE)

.swap_get <- function(swap, k) {
  for (p in swap) if (identical(p[[1L]]$key, k$key)) return(p[[2L]])
  k
}

.TRIG_NAMES <- c("sin", "cos", "tan", "cot", "sec", "csc")

.reim_atoms <- function(e) {
  out <- character()
  walk <- function(z) {
    if (z$type == "Function" && z$name %in% c("re", "im", "arg", "atan2")) out <<- c(out, z$key)
    for (a in z$args) walk(a)
  }
  walk(e)
  unique(out)
}

.sp_find <- function(e, pred) {
  out <- list()
  walk <- function(z) {
    if (pred(z)) out[[z$key]] <<- z
    for (a in z$args) walk(a)
  }
  walk(e)
  unname(out)
}

# check_assumptions(val, **symb.assumptions0)
.sp_check_assumptions <- function(val, symb) {
  given <- if (is_symbol_like(symb)) symb$extra$facts else NULL
  if (is.null(given)) return(TRUE)
  given <- given[!is.na(given)]
  given <- given[names(given) != "commutative"]
  if (!length(given)) return(TRUE)
  f <- .facts(val)
  res <- TRUE
  for (nm in names(given)) {
    v <- f[[nm]]
    if (is.null(v) || is.na(v)) { res <- NA; next }
    if (v != isTRUE(given[[nm]])) return(FALSE)
  }
  res
}

.sp_recast <- function(eqs, symbols) {
  orig <- symbols
  srt <- .sp_ordered(symbols)
  swap <- list()
  i <- 0L
  for (s in srt) {
    # a Dummy or Wild is a Symbol subclass
    if (!(s$type %in% c("Symbol", "Dummy", "Wild")) && !any(vapply(swap, function(p) identical(p[[1L]]$key, s$key), TRUE))) {
      swap[[length(swap) + 1L]] <- list(s, Dummy(paste0("X", i)))
      i <- i + 1L
    }
  }
  newf <- lapply(eqs, function(e) if (length(swap)) subs(e, swap) else e)
  list(f = newf, symbols = lapply(orig, function(s) .swap_get(swap, s)),
       restore = lapply(swap, function(p) list(p[[2L]], p[[1L]])))
}

.sps_gcd_terms <- function(xs) Reduce(function(a, b) gcd(a, b), xs)

# ---------------------------------------------------------------------------
# _solve: one equation

.sp_vsolve <- function(e, s, flags) lapply(.sp__solve(e, list(s), flags), function(d) .sd_get(d, s))

.sp__solve <- function(f, symbols, flags) {
  not_impl <- function(f) paste0("No algorithms are implemented to solve equation ", .sp_py_str(f))
  if (length(symbols) != 1L) {
    nonlin_s <- list()
    got_s <- list()
    result <- list()
    for (s in symbols) {
      sl <- .sp_solve_linear(f, symbols = list(s))
      xi <- sl[[1L]]; v <- sl[[2L]]
      if (identical(xi$key, s$key)) {
        if (!isFALSE(.fl(flags, "simplify", TRUE))) v <- simplify(v)
        vfree <- free_symbols(v)
        if (any(vapply(vfree, function(z) .sp_in(z, got_s), TRUE))) next
        got_s[[length(got_s) + 1L]] <- xi
        result[[length(result) + 1L]] <- .sd1(xi, v)
      } else if (.sp_bool(xi)) nonlin_s[[length(nonlin_s) + 1L]] <- s
    }
    if (!length(nonlin_s)) return(result)
    for (s in nonlin_s) {
      soln <- tryCatch(.sp__solve(f, list(s), flags), NotImplementedError = function(e) NULL)
      if (is.null(soln)) next
      for (sol in soln) {
        if (any(vapply(free_symbols(.sd_get(sol, s)), function(z) .sp_in(z, got_s), TRUE))) next
        if (!.sp_in(s, got_s)) got_s[[length(got_s) + 1L]] <- s
        result[[length(result) + 1L]] <- sol
      }
    }
    if (length(got_s)) return(result)
    .sp_stop("NotImplementedError", not_impl(f))
  }
  symbol <- symbols[[1L]]
  f <- .map_bottom_up(f, function(e) if (e$type == "Function" && e$name == "binomial" && .has1(e, symbol)) expand_func(e) else e)
  check <- !isFALSE(.fl(flags, "check", TRUE))
  flags$check <- check
  checkdens <- check
  result <- NULL
  msg <- ""
  f_num <- NULL
  if (f$type == "Mul") {
    vals <- list()
    for (m in .sympy_args(f)) {
      if (m$type == "Constant" && m$name %in% c("-oo", "zoo", "oo")) {
        vals <- list()
        break
      }
      vals <- c(vals, .sp_vsolve(m, symbol, flags))
    }
    vals <- .sp_uniq(vals)
    result <- lapply(vals, function(v) .sd1(symbol, v))
    if (check) {
      dens <- .fl(flags, "_denominators", .sp_simple_dens(f, symbols))
      result <- Filter(function(s) {
        !any(vapply(dens, function(den) .sp_reject_by_denominator(den, s, flags), TRUE)) &&
          !any(vapply(.sympy_args(f), function(m) .sp_reject_by_factor_finiteness(m, s, flags), TRUE))
      }, result)
    }
    check <- FALSE
    flags$simplify <- FALSE
  } else if (f$type == "Piecewise") {
    res <- list()
    if (any(vapply(f$args, function(p) isTRUE(.sp_fact(p[[1L]], "zero")), TRUE))) f <- simplify(f)
    cond <- S_true
    neg <- S_true
    for (pr in .sps_pw_pairs(f)) {
      expr <- pr[[1L]]; cnd <- pr[[2L]]
      cond <- And(neg, cnd)
      neg <- And(neg, Not(cond))
      if (isTRUE(.sp_fact(expr, "zero")) && !.is_false(simplify(cond)))
        .sp_stop("NotImplementedError", .filldedent(paste0("An expression is already zero when ", .sp_py_str(cond),
          ". This means that in this *region* the solution is zero but solve can only represent discrete, not interval, solutions. If this is a spurious interval it might be resolved with simplification of the Piecewise conditions.")))
      candidates <- .sp_vsolve(expr, symbol, flags)
      for (candidate in candidates) {
        if (.sp_in(candidate, res)) next
        v <- tryCatch(subs(cond, symbol, candidate), error = function(e) NULL)
        if (is.null(v)) next
        v <- tryCatch(simplify(v, ratio = 2), error = function(e) v)
        if (.is_false(v)) next
        if (.is_true(v)) res[[length(res) + 1L]] <- candidate
        else res[[length(res) + 1L]] <- Piecewise(list(candidate, v), list(S_NaN, S_true))
      }
    }
    return(lapply(res, function(r) .sd1(symbol, r)))
  } else {
    sl <- .sp_solve_linear(f, symbols = symbols)
    f_num <- sl[[1L]]; sol <- sl[[2L]]
    if (!.sp_bool(f_num) || .sp_nan(sol)) return(list())
    if (is_symbol_like(f_num)) {
      if (!isFALSE(.fl(flags, "simplify", TRUE))) sol <- simplify(sol)
      return(list(.sd1(f_num, sol)))
    }
    poly_gen <- NULL
    if (f_num$type != "Add") {
      add_args <- Filter(function(i) .has1(i, symbol), atoms(f_num, "Add"))
      if (length(add_args) == 1L) {
        gen <- add_args[[1L]]
        spart <- .as_base_exp(.sp_indep(gen, list(symbol))[[2L]])[[1L]]
        if (identical(spart$key, symbol$key)) poly_gen <- list(spart)
      }
    }
    pgens <- if (!is.null(poly_gen)) poly_gen else .sp_poly_gens(f_num)
    if (is.null(pgens)) {
      if (is.null(poly_gen) && !length(free_symbols(f_num))) {
        sf <- simplify(f_num)
        if (!identical(sf$key, f_num$key)) return(.sp__solve(sf, symbols, flags))
        .sp_stop("ValueError", "expression appears to be a constant")
      }
    }
    gens <- Filter(function(g) .has1(g, symbol), pgens)
    as_base_q <- function(x) {
      be <- .as_base_exp(x)
      b <- be[[1L]]; e <- be[[2L]]
      if (is_exact_number(e)) return(list(b, gmp::denominator(e$value)))
      if (e$type != "Mul") return(list(x, gmp::as.bigz(1)))
      cm <- as_coeff_mul_split(e)
      c <- cm[[1L]]; ee <- cm[[2L]]
      if (is_exact_number(c) && !is_one_num(c)) return(list(.pow(b, ee), gmp::denominator(c$value)))
      list(x, gmp::as.bigz(1))
    }
    if (length(gens) > 1L) {
      bq <- lapply(gens, as_base_q)
      bases <- .sp_uniq(lapply(bq, `[[`, 1L))
      qs <- vapply(bq, function(z) as.numeric(z[[2L]]), 0)
      if (length(bases) > 1L || !all(qs == 1)) {
        funcs <- Filter(function(b) b$type == "Function", bases)
        trig <- Filter(function(b) b$name %in% .TRIG_NAMES, funcs)
        other <- Filter(function(b) !(b$name %in% .TRIG_NAMES), funcs)
        if (!length(other) && length(trig) > 1L) {
          newf <- NULL
          if (f_num$type == "Add" && length(f_num$args) == 2L) {
            args_ <- .sympy_args(f_num)
            t2 <- lapply(args_, function(i) Filter(function(z) z$name %in% .TRIG_NAMES, atoms(i, "Function")))
            if (all(lengths(t2) == 1L)) {
              a <- t2[[1L]][[1L]]; b <- t2[[2L]][[1L]]
              if (a$name == "cos") { tmp <- a; a <- b; b <- tmp; args_ <- rev(args_) }
              if (a$name == "sin" && b$name == "cos" && identical(a$args[[1L]]$key, b$args[[1L]]$key)) {
                nd <- .as_numer_denom(.add(list(TR2i(.div(args_[[1L]], args_[[2L]])), S_One)))
                newf <- nd[[1L]]
                if (!is_number_type(nd[[2L]])) newf <- NULL
              }
            }
          }
          if (is.null(newf)) newf <- rewrite(TR1(f_num), "tan")
          if (!identical(newf$key, f_num$key)) {
            fl2 <- flags
            fl2$check <- FALSE
            result <- .sp__solve(newf, list(symbol), fl2)
          }
        }
        if (is.null(result) && length(funcs)) {
          funcs <- .sp_ordered(funcs)
          f1 <- funcs[[1L]]
          t <- Dummy("t")
          ftry <- subs(f_num, f1, t)
          if (!.has1(ftry, symbol)) {
            cv_sols <- .sp__solve(ftry, list(t), flags)
            cv_inv <- .sp_ordered(.sp_vsolve(.sub(t, f1), symbol, flags))[[1L]]
            result <- lapply(cv_sols, function(sol) .sd1(symbol, .sd_subs(cv_inv, sol)))
          }
        }
        if (is.null(result)) msg <- paste0("multiple generators [", paste(vapply(gens, .sp_py_str, ""), collapse = ", "), "]")
      } else {
        u <- bases[[1L]]
        t <- Dummy("t")
        inv <- .sp_vsolve(.sub(u, t), symbol, flags)
        if (u$type == "Pow" || (u$type == "Function" && u$name == "exp")) {
          ftry <- .map_bottom_up(f_num, function(w) {
            if (w$type == "Pow" || (w$type == "Function" && w$name == "exp")) {
              be <- .as_base_exp(w)
              expand_power_exp(.pow(be[[1L]], expand_mul(be[[2L]])))
            } else w
          })
          ftry <- subs(ftry, u, t)
          if (!.has1(ftry, symbol)) {
            soln <- .sp__solve(ftry, list(t), flags)
            result <- list()
            for (i in inv) for (s in soln) result[[length(result) + 1L]] <- .sd1(symbol, .sd_subs(i, s))
          }
        }
      }
    } else if (length(gens) == 1L) {
      g <- gens[[1L]]
      tsolve_pass <- isTRUE(flags$tsolve)
      flags$tsolve <- NULL
      if (!tsolve_pass) {
        u <- Dummy("u")
        pe <- if (identical(g$key, symbol$key)) f_num else subs(f_num, g, u)
        pv <- if (identical(g$key, symbol$key)) symbol else u
        if (!.has1(pe, symbol) || identical(g$key, symbol$key)) {
          deg <- degree(expand(pe), pv)
          flags$tsolve <- TRUE
          rr <- .roots_sympy(pe, pv, cubics = flags$cubics, quartics = flags$quartics, quintics = flags$quintics)
          soln <- unclass(rr)
          soln <- soln[seq_along(soln)]
          mult <- attr(rr, "multiplicity")
          if (sum(mult) < deg) {
            soln <- tryCatch(.sp_all_roots(pe, pv), NotImplementedError = function(e) {
              if (isFALSE(.fl(flags, "incomplete", TRUE)))
                .sp_stop("NotImplementedError", .filldedent("Neither high-order multivariate polynomials nor sorting of EX-domain polynomials is supported. If you want to see any results, pass keyword incomplete=True to solve; to see numerical values of roots for univariate expressions, use nroots."))
              soln
            })
          }
          if (!is.null(soln)) {
            if (!identical(g$key, symbol$key)) {
              t <- Dummy("t")
              inv <- tryCatch(.sp_vsolve(.sub(g, t), symbol, flags), NotImplementedError = function(e) NULL)
              if (is.null(inv)) soln <- NULL else {
                ss <- list()
                for (i in inv) for (s in soln) ss[[length(ss) + 1L]] <- subs(i, t, s)
                soln <- .sp_uniq(ss)
              }
            } else {
              check <- FALSE
            }
            if (!is.null(soln) && length(soln) > 2L && is.null(flags$simplify)) flags$simplify <- FALSE
          }
          if (!is.null(soln)) result <- lapply(soln, function(v) .sd1(symbol, v))
        }
      }
    }
  }
  if (is.null(result)) {
    unrad_flag <- .fl(flags, "_unrad", TRUE)
    flags[["_unrad"]] <- NULL
    if (!isFALSE(unrad_flag)) {
      u <- tryCatch(.sp_unrad(f_num, list(symbol)), ValueError = function(e) NULL, NotImplementedError = function(e) NULL)
      if (!is.null(u)) {
        eq <- u[[1L]]; cov <- u[[2L]]
        rv <- NULL
        if (length(cov)) {
          isym <- cov[[1L]]; ieq <- cov[[2L]]
          inv <- .sp_vsolve(ieq, symbol, flags)[[1L]]
          rv <- .sp_uniq(lapply(.sp__solve(eq, list(isym), flags), function(xi) .sd_subs(inv, xi)))
        } else {
          rv <- tryCatch(.sp_uniq(.sp_vsolve(eq, symbol, flags)), NotImplementedError = function(e) NULL)
        }
        if (!is.null(rv)) {
          result <- lapply(rv, function(v) .sd1(symbol, v))
          if (is.null(flags$simplify)) flags$simplify <- FALSE
        }
      }
    }
  }
  if (is.null(result)) {
    flags$tsolve <- NULL
    soln <- tryCatch(.sp_tsolve(f_num, symbol, flags), PolynomialError = function(e) NULL)
    if (!is.null(soln)) result <- lapply(soln, function(v) .sd1(symbol, v))
  }
  if (is.null(result)) .sp_stop("NotImplementedError", paste0(msg, "\n", not_impl(f)))
  result <- .sd_unique(result)
  if (!isFALSE(.fl(flags, "simplify", TRUE))) {
    result <- lapply(result, function(d) .sd_map(d, simplify))
    result <- .sd_unique(result)
    flags$simplify <- FALSE
  }
  if (checkdens) {
    dens <- .sp_simple_dens(f, symbols)
    result <- Filter(function(r) !any(vapply(dens, function(d) .sp_reject_by_denominator(d, r, flags), TRUE)), result)
  }
  if (check) result <- Filter(function(r) !isFALSE(.sp_checksol(f_num, r, flags)), result)
  result
}

# rewrite(Piecewise) as SymPy: Abs only for a real or imaginary argument
.sps_pw_rewrite <- function(e) .map_bottom_up(e, function(z) {
  if (z$type != "Function") return(z)
  if (z$name == "Abs") {
    a <- z$args[[1L]]
    if (isTRUE(.sp_fact(a, "extended_real"))) return(Piecewise(list(a, Ge(a, 0)), list(.neg(a), TRUE)))
    if (isTRUE(.sp_fact(a, "imaginary"))) {
      ia <- .mul(list(S_I, a))
      return(Piecewise(list(ia, Ge(ia, 0)), list(.neg(ia), TRUE)))
    }
    return(z)
  }
  if (z$name %in% c("Min", "Max")) return(rewrite(z, "Piecewise"))
  z
})

# str() as SymPy writes it in messages: powers with **
.sp_py_str <- function(e) gsub("^", "**", .str(e), fixed = TRUE)

.sps_pw_pairs <- function(pw) lapply(pw$args, function(p) if (is_sym(p) && p$type == "ExprCondPair") p$args else p)

# Poly(f).gens; NULL when f has no generators
.sp_poly_gens <- function(f) {
  g <- .find_gens(list(expand(f)))
  if (!length(g)) NULL else g
}

.sp_all_roots <- function(pe, x) {
  cf <- .coeff_list(expand(pe), x)
  if (!all(vapply(cf, is_exact_number, TRUE)))
    .sp_stop("NotImplementedError", "sorting of these roots is not supported")
  roots(pe, x, multiple = TRUE, rootof = TRUE)
}

# ---------------------------------------------------------------------------
# _tsolve

.sp_tsolve <- function(eq, sym, flags) {
  if (is.null(flags$tsolve_saw)) flags$tsolve_saw <- new.env()
  if (exists(eq$key, envir = flags$tsolve_saw, inherits = FALSE)) return(NULL)
  assign(eq$key, TRUE, envir = flags$tsolve_saw)
  inv <- .sp_invert(eq, list(sym))
  rhs <- inv[[1L]]; lhs <- inv[[2L]]
  if (identical(lhs$key, sym$key)) return(list(rhs))
  if (!.has_free(lhs, list(sym)) && !.has_free(rhs, list(sym))) return(NULL)
  r <- tryCatch(.sp_tsolve_main(eq, sym, lhs, rhs, flags), NotImplementedError = function(e) "notimpl")
  if (!identical(r, "notimpl") && !identical(r, "fallthrough")) return(r)
  # flags.pop(): the flag is removed whether or not the branch is taken
  biv <- .fl(flags, "bivariate", TRUE)
  flags$bivariate <- NULL
  if (!isFALSE(biv)) {
    r <- .sp_tsolve_bivariate(eq, sym, lhs, rhs, flags)
    if (!is.character(r)) return(r)
    if (!is.null(attr(r, "lhs"))) {
      lhs <- attr(r, "lhs")
      rhs <- attr(r, "rhs")
    }
  }
  frc <- .fl(flags, "force", TRUE)
  flags$force <- NULL
  if (!isFALSE(frc)) {
    flags$force <- FALSE
    ps <- posify(.sub(lhs, rhs))
    pos <- ps$expr
    reps <- ps$reps
    if (rhs$type == "Constant" && rhs$name == "zoo") return(list())
    u <- sym
    for (pr in reps) if (identical(pr[[2L]]$key, sym$key)) { u <- pr[[1L]]; break }
    if (.has1(pos, u)) {
      soln <- tryCatch(.sp_vsolve(pos, u, flags), NotImplementedError = function(e) NULL)
      if (!is.null(soln)) return(lapply(soln, function(s) subs(s, reps)))
    }
  }
  NULL
}

.sp_multi_inverses <- function(nm, rhs) {
  switch(nm,
    sin = list(fcall("asin", list(rhs)), .sub(S_Pi, fcall("asin", list(rhs)))),
    cos = list(fcall("acos", list(rhs)), .sub(.mul(list(S_Two, S_Pi)), fcall("acos", list(rhs)))),
    NULL)
}

.sp_tsolve_main <- function(eq, sym, lhs, rhs, flags) {
  if (lhs$type == "Mul" && isTRUE(.sp_fact(rhs, "positive")) && !.has_free(rhs, list(sym))) {
    factors <- lapply(.sympy_args(lhs), .as_base_exp)
    powers <- Filter(function(be) .has_free(be[[2L]], list(sym)), factors)
    if (length(powers) > 1L && all(vapply(factors, function(be) isTRUE(.sp_fact(be[[1L]], "positive")) &&
                                            !.has_free(be[[1L]], list(sym)) && isTRUE(.sp_fact(be[[2L]], "extended_real")), TRUE))) {
      lhs <- .add(lapply(factors, function(be) .mul(list(be[[2L]], fcall("log", list(be[[1L]]))))))
      rhs <- fcall("log", list(rhs))
    }
  }
  if (lhs$type == "Add") {
    f <- sym_factor(powdenest(.sub(lhs, rhs)))
    if (f$type == "Mul") return(.sp_vsolve(f, sym, flags))
    if (.sp_bool(rhs)) {
      f <- logcombine(lhs, force = !isFALSE(.fl(flags, "force", TRUE)))
      if (.count_fn(f, "log") != .count_fn(lhs, "log")) {
        if (f$type == "Function" && f$name == "log") return(.sp_vsolve(.sub(f$args[[1L]], fcall("exp", list(rhs))), sym, flags))
        return(.sp_tsolve(.sub(f, rhs), sym, flags))
      }
    }
  } else if (lhs$type == "Pow") {
    b <- lhs$args[[1L]]; ex <- lhs$args[[2L]]
    if (ex$type == "Integer") {
      if (!identical(.sub(lhs, rhs)$key, eq$key)) return(.sp_vsolve(.sub(lhs, rhs), sym, flags))
    }
    if (!.has_free(ex, list(sym))) return(.sp_vsolve(.sub(b, .pow(rhs, .inv(ex))), sym, flags))
    if (any(vapply(free_symbols(rhs), function(t) t$type == "Dummy", TRUE))) .sp_stop("NotImplementedError", "")
    if (!.sp_bool(rhs)) {
      sol_base <- .sp_vsolve(b, sym, flags)
      return(Filter(function(s) !is_zero_num(subs(ex, sym, s)), sol_base))
    }
    if (!.has1(b, sym)) {
      if (is_zero_num(b)) return(if (.sp_bool(rhs)) .sp_vsolve(ex, sym, flags) else list())
      rbe <- .as_base_exp(rhs)
      if (identical(b$key, rbe[[1L]]$key)) {
        sol <- .sp_vsolve(.sub(ex, rbe[[2L]]), sym, flags)
      } else {
        f <- .sub(fcall("exp", list(.mul(list(fcall("log", list(b)), ex)))), fcall("exp", list(fcall("log", list(rhs)))))
        sol <- .sp_vsolve(f, sym, flags)
      }
      equal <- function(e1, e2) {
        d <- Dummy()
        r <- .sp_checksol(.sub(e1, d), .sd1(d, e2))
        if (is.na(r)) {
          if (!identical(.sps_nsimplify(e1)$key, .sps_nsimplify(e2)$key)) return(FALSE)
          r <- .expr_equals(e1, e2)
        }
        r
      }
      e_rat <- .sps_nsimplify(.div(fcall("log", list(fcall("Abs", list(rhs)))), fcall("log", list(fcall("Abs", list(b))))))
      e_rat <- simplify(posify(e_rat)$expr)
      nd <- fraction(e_rat)
      if (is_zero_num(expand(.sub(.pow(b, nd[[1L]]), .pow(rhs, nd[[2L]]))))) {
        sol <- Filter(function(s) !isTRUE(equal(subs(ex, sym, s), e_rat)), sol)
        sol <- c(sol, .sp_vsolve(.sub(ex, e_rat), sym, flags))
      }
      return(.sp_uniq(sol))
    } else {
      sol <- list()
      logform <- .sub(.mul(list(ex, fcall("log", list(b)))), fcall("log", list(rhs)))
      if (!identical(logform$key, .sub(lhs, rhs)$key)) {
        s2 <- tryCatch(.sp_vsolve(logform, sym, flags), NotImplementedError = function(e) list())
        sol <- c(sol, s2)
      }
      check <- list()
      if (is_one_num(rhs)) {
        check <- c(check, .sp_vsolve(ex, sym, flags), .sp_vsolve(.sub(b, S_One), sym, flags), .sp_vsolve(.add(list(b, S_One)), sym, flags))
      } else if (is_exact_number(rhs)) {
        p <- gmp::numerator(rhs$value); q <- gmp::denominator(rhs$value)
        for (d in divisors(abs(p))) {
          if (d == 1) next
          il <- integer_log(p, d)
          e <- il[[1L]]; t <- il[[2L]]
          if (!isTRUE(t)) next
          for (s in divisors(abs(q))) {
            if (gmp::as.bigz(s)^e == q) {
              r <- make_rational(gmp::as.bigq(d, s))
              check <- c(check, .sp_vsolve(.sub(b, r), sym, flags), .sp_vsolve(.add(list(b, r)), sym, flags),
                         .sp_vsolve(.sub(ex, .int(e)), sym, flags))
            }
          }
        }
      } else if (isTRUE(.sp_fact(rhs, "irrational"))) {
        bl <- .as_base_exp(b)
        nd <- .as_numer_denom(.mul(list(bl[[2L]], ex)))
        rb <- .as_base_exp(sqrtdenest(rhs))
        check <- lapply(.sp_vsolve(.sub(b, rb[[1L]]), sym, flags), sqrtdenest)
        check <- c(check, lapply(.sp_vsolve(.sub(ex, rb[[2L]]), sym, flags), sqrtdenest))
        if (!is_one_num(.mul(list(bl[[2L]], nd[[2L]]))))
          check <- c(check, .sp_vsolve(.sub(.pow(bl[[1L]], nd[[1L]]), .pow(rhs, .mul(list(bl[[2L]], nd[[2L]])))), sym, flags))
      }
      for (s in check) {
        ok <- .sp_checksol(eq, .sd1(sym, s))
        if (is.na(ok)) ok <- isTRUE(.expr_equals(subs(eq, sym, s), S_Zero))
        if (isTRUE(ok)) sol[[length(sol) + 1L]] <- s
      }
      return(.sp_uniq(sol))
    }
  } else if (lhs$type == "Function" && length(lhs$args) == 1L) {
    mi <- .sp_multi_inverses(lhs$name, rhs)
    if (!is.null(mi)) {
      soln <- list()
      for (i in mi) soln <- c(soln, .sp_vsolve(.sub(lhs$args[[1L]], i), sym, flags))
      return(.sp_uniq(soln))
    } else if (lhs$name == "LambertW") {
      return(.sp_vsolve(.sub(lhs$args[[1L]], .mul(list(rhs, fcall("exp", list(rhs))))), sym, flags))
    }
  }
  rw <- .sps_rewrite_exp(lhs)
  if (!identical(rw$key, lhs$key)) return(.sp_vsolve(.sub(rw, rhs), sym, flags))
  "fallthrough"
}

# rebuild(e.rewrite(exp)) as SymPy: a power with a symbolic base and exponent
# becomes exp(log(base)*exp); one with a constant exponent re-evaluates to
# itself, so it is left alone
.sps_rewrite_exp <- function(e) {
  e <- rewrite(e, "exp")
  .map_bottom_up(e, function(z) {
    if (z$type != "Pow") return(z)
    b <- z$args[[1L]]; p <- z$args[[2L]]
    if (is_zero_num(b) || .has_fn(b, "exp") || .has_fn(p, "exp")) return(z)
    if (length(free_symbols(b)) && length(free_symbols(p)))
      return(fcall("exp", list(.mul(list(fcall("log", list(b)), p)))))
    z
  })
}

.count_fn <- function(e, nm) {
  n <- 0L
  walk <- function(z) {
    if (z$type == "Function" && z$name == nm) n <<- n + 1L
    for (a in z$args) walk(a)
  }
  walk(e)
  n
}


# ---------------------------------------------------------------------------
# solve_undetermined_coeffs

.sp_solve_undetermined <- function(g, symbols, flags) {
  fs <- free_symbols(g)
  if (any(vapply(fs, function(s) !.sp_in(s, symbols), TRUE))) {
    fl <- flags
    fl$dict <- TRUE
    fl$set <- NULL
    sol <- .sp_solve_undetermined_coeffs(g, symbols, list(), fl, raw = TRUE)
    if (is.list(sol) && length(sol) == 1L) return(sol)
  }
  NULL
}

.sp_solve_undetermined_coeffs <- function(equ, coeffs, syms, flags, raw = FALSE) {
  if (!length(coeffs) || !all(vapply(coeffs, is_symbol_like, TRUE))) .sp_stop("ValueError", "must provide symbols for coeffs")
  eq <- if (equ$type == "Equality") .sub(equ$args[[1L]], equ$args[[2L]]) else equ
  ceq <- cancel(eq)
  xeq <- .mexpand(.as_numer_denom(ceq)[[1L]])
  free <- free_symbols(xeq)
  coeffs <- Filter(function(s) .sp_in(s, coeffs), free)
  if (!length(coeffs)) return(if (isTRUE(flags$set)) .pytuple(list(list(), .sps_pyset(list()))) else list())
  if (!length(syms)) {
    id <- .sp_indep(xeq, coeffs, as_add = TRUE)
    dfree <- free_symbols(id[[2L]])
    syms <- Filter(function(s) .sp_in(s, free_symbols(id[[1L]])), dfree)
    if (!length(syms)) syms <- Filter(function(s) !.sp_in(s, coeffs), dfree)
    if (!length(syms)) syms <- list(Dummy())
  } else {
    if (length(syms) == 1L && .sp_iterable(syms[[1L]])) syms <- .sp_as_list(syms[[1L]])
    rc <- .sp_recast(list(xeq), syms)
    xeq <- rc$f[[1L]]
    syms <- rc$symbols
  }
  cd <- .sp_coefficients_dict(xeq, syms)
  gens <- Filter(function(g) !is_one_num(g), cd$k)
  if (any(vapply(gens, function(g) .has_free(g, coeffs), TRUE))) return(NULL)
  rc <- .sp_recast(list(xeq), gens)
  xeq <- rc$f[[1L]]
  gens <- rc$symbols
  system <- .sp_collect_values(xeq, gens)
  fl <- flags
  if (raw) {
    fl$dict <- TRUE
    return(.sp_solve_raw(system, coeffs, fl))
  }
  soln <- .sp_solve(system, list(.sps_pyset(coeffs)), fl)
  settings <- isTRUE(flags$dict) || isTRUE(flags$set)
  if (!is.list(soln) || !is.null(names(soln)) || settings || length(soln) != 1L) return(soln)
  soln[[1L]]
}

# solve() returning the internal list of solution dicts (for dict=True callers)
.sp_solve_raw <- function(system, coeffs, flags) {
  flags$dict <- TRUE
  flags[["_raw"]] <- TRUE
  .sp_solve(system, coeffs, flags)
}

.sp_find_key <- function(syms, n) {
  for (s in syms) if (identical(.str(s), n)) return(s)
  Symbol(n)
}

# as_coefficients_dict(*syms): term generators in syms with their coefficients
.sp_coefficients_dict <- function(e, syms) {
  out <- .sd_new()
  for (t in .aargs(expand(e))) {
    id <- .sp_indep(t, syms, as_add = FALSE)
    g <- id[[2L]]
    old <- .sd_get(out, g)
    out <- .sd_set(out, g, if (is.null(old)) id[[1L]] else .add(list(old, id[[1L]])))
  }
  out
}

# list(collect(e, gens, evaluate=False).values())
.sp_collect_values <- function(e, gens) {
  out <- .sd_new()
  for (t in .aargs(expand(e))) {
    id <- .sp_indep(t, gens, as_add = FALSE)
    g <- id[[2L]]
    old <- .sd_get(out, g)
    out <- .sd_set(out, g, if (is.null(old)) id[[1L]] else .add(list(old, id[[1L]])))
  }
  out$v
}

# ---------------------------------------------------------------------------
# Polynomial view of an expression in the given symbols (g.as_poly(*symbols)):
# list of terms, each list(monom = integer vector, coeff = sym), or NULL when
# g is not a polynomial in the symbols.

.sp_as_poly <- function(g, symbols) {
  e <- expand(g)
  n <- length(symbols)
  gk <- vapply(symbols, .k, "")
  terms <- list()
  keys <- character()
  for (t in .aargs(e)) {
    if (is_zero_num(t)) next
    coeff <- list()
    mon <- integer(n)
    for (fac in .margs(t)) {
      if (is_number_type(fac)) {
        coeff[[length(coeff) + 1L]] <- fac
        next
      }
      dp <- .decompose_power(fac)
      base <- dp[[1L]]
      ex <- dp[[2L]]
      if (ex < 0) {
        ex <- -ex
        base <- .pow(base, S_NegativeOne)
      }
      j <- match(base$key, gk)
      if (!is.na(j)) {
        mon[j] <- as.integer(ex)
      } else if (!.has_free(fac, symbols)) {
        coeff[[length(coeff) + 1L]] <- fac
      } else return(NULL)
    }
    k <- paste(mon, collapse = ",")
    i <- match(k, keys)
    cf <- .mul(coeff)
    if (is.na(i)) {
      keys <- c(keys, k)
      terms[[length(terms) + 1L]] <- list(monom = mon, coeff = cf)
    } else terms[[i]]$coeff <- .add(list(terms[[i]]$coeff, cf))
  }
  terms <- Filter(function(t) !is_zero_num(t$coeff), terms)
  structure(list(terms = terms, gens = symbols), class = "sp_poly")
}

.sp_poly_degree <- function(p, j = 1L) if (!length(p$terms)) -Inf else max(vapply(p$terms, function(t) t$monom[j], 0L))

.sp_poly_is_linear <- function(p) all(vapply(p$terms, function(t) sum(t$monom) <= 1L, TRUE))
.sp_poly_expr <- function(p) .add(lapply(p$terms, function(t)
  .mul(c(list(t$coeff), Map(function(g, k) .pow(g, .int(k)), p$gens[t$monom > 0], t$monom[t$monom > 0])))))
.sp_poly_free <- function(p) unique_syms(c(lapply(p$terms, function(t) free_symbols(t$coeff)),
                                           list(p$gens[vapply(seq_along(p$gens), function(j) any(vapply(p$terms, function(t) t$monom[j] > 0, TRUE)), TRUE)])))

# ---------------------------------------------------------------------------
# Linear systems: sympy.polys.solvers.solve_lin_sys over the coefficient field

.sps_fe <- function(e) if (is_exact_number(e)) e else cancel(e)
.sps_fe_add <- function(a, b) .sps_fe(.add(list(a, b)))
.sps_fe_sub <- function(a, b) .sps_fe(.sub(a, b))
.sps_fe_mul <- function(a, b) .sps_fe(.mul(list(a, b)))
.sps_fe_div <- function(a, b) .sps_fe(.div(a, b))

# reduced row echelon form of a list of rows (lists of field elements);
# returns list(rows, pivots) with 1-based pivot columns
.sps_rref <- function(M, ncol) {
  nrow <- length(M)
  piv <- integer()
  r <- 1L
  for (j in seq_len(ncol)) {
    if (r > nrow) break
    pr <- NA_integer_
    for (i in r:nrow) if (!is_zero_num(M[[i]][[j]])) { pr <- i; break }
    if (is.na(pr)) next
    if (pr != r) { tmp <- M[[r]]; M[[r]] <- M[[pr]]; M[[pr]] <- tmp }
    pv <- M[[r]][[j]]
    if (!is_one_num(pv)) M[[r]] <- lapply(M[[r]], function(z) if (is_zero_num(z)) z else .sps_fe_div(z, pv))
    for (i in seq_len(nrow)) {
      if (i == r) next
      c <- M[[i]][[j]]
      if (is_zero_num(c)) next
      M[[i]] <- lapply(seq_len(ncol), function(k) {
        rk <- M[[r]][[k]]
        if (is_zero_num(rk)) M[[i]][[k]] else .sps_fe_sub(M[[i]][[k]], .sps_fe_mul(c, rk))
      })
    }
    piv <- c(piv, j)
    r <- r + 1L
  }
  list(M, piv)
}

# eqs: list of sp_poly (linear). Returns a .sd dict or NULL (inconsistent).
.sp_solve_lin_sys <- function(polys, symbols) {
  n <- length(symbols)
  coeffs <- list()
  rhs <- list()
  for (p in polys) {
    co <- vector("list", n)
    r <- S_Zero
    any_c <- FALSE
    for (t in p$terms) {
      if (sum(t$monom) == 0L) r <- .sps_fe(.add(list(r, t$coeff)))
      else {
        j <- which(t$monom == 1L)
        co[[j]] <- t$coeff
        any_c <- TRUE
      }
    }
    if (!any_c) {
      if (!.sp_bool(r)) next
      return(NULL)
    }
    coeffs[[length(coeffs) + 1L]] <- co
    rhs[[length(rhs) + 1L]] <- r
  }
  # connected components of the symbols linked by equations
  comp <- seq_len(n)
  findc <- function(i) { while (comp[i] != i) i <- comp[i]; i }
  for (co in coeffs) {
    js <- which(!vapply(co, is.null, TRUE))
    if (length(js) > 1L) for (j in js[-1L]) {
      a <- findc(js[1L]); b <- findc(j)
      if (a != b) comp[max(a, b)] <- min(a, b)
    }
  }
  roots_ <- vapply(seq_len(n), findc, 0L)
  sol <- .sd_new()
  for (cr in unique(roots_)) {
    idx <- which(vapply(coeffs, function(co) { j <- which(!vapply(co, is.null, TRUE))[1L]; roots_[j] == cr }, TRUE))
    if (!length(idx)) next
    M <- lapply(idx, function(i) c(lapply(coeffs[[i]], function(z) if (is.null(z)) S_Zero else .sps_fe(z)), list(.sps_fe(.neg(rhs[[i]])))))
    rr <- .sps_rref(M, n + 1L)
    E <- rr[[1L]]; piv <- rr[[2L]]
    if (length(piv) && piv[length(piv)] == n + 1L) return(NULL)
    if (length(piv) == n) {
      for (i in seq_len(n)) sol <- .sd_set(sol, symbols[[i]], E[[i]][[n + 1L]])
    } else {
      for (i in seq_along(piv)) {
        p <- piv[i]
        terms <- list(E[[i]][[n + 1L]])
        if (p < n) for (j in (p + 1L):n) {
          if (!is_zero_num(E[[i]][[j]])) terms[[length(terms) + 1L]] <- .neg(.mul(list(E[[i]][[j]], symbols[[j]])))
        }
        sol <- .sd_set(sol, symbols[[p]], .add(terms))
      }
    }
  }
  keep <- vapply(seq_along(sol$k), function(i) !identical(sol$k[[i]]$key, sol$v[[i]]$key), TRUE)
  list(k = sol$k[keep], v = sol$v[keep])
}

# ---------------------------------------------------------------------------
# Polynomial systems: sympy/solvers/polysys.py

.sp_poly_rational <- function(p) all(vapply(p$terms, function(t) is_exact_number(t$coeff), TRUE))

.sp_solve_poly_system <- function(polys, gens, strict = FALSE) {
  exprs <- lapply(polys, function(p) .sp_poly_expr(structure(list(terms = p$terms, gens = gens), class = "sp_poly")))
  degs_ok <- function(e) all(vapply(gens, function(g) degree(e, g) <= 2, TRUE))
  if (length(exprs) == 2L && length(gens) == 2L && all(vapply(exprs, degs_ok, TRUE))) {
    r <- tryCatch(.sp_solve_biquadratic(exprs[[1L]], exprs[[2L]], gens), SolveFailed = function(e) "failed")
    if (!identical(r, "failed")) return(r)
  }
  .sp_solve_generic(exprs, gens, strict)
}

.sp_groebner <- function(F, gens) {
  ps <- lapply(F, function(f) .sp_as_poly(f, gens))
  if (any(vapply(ps, is.null, TRUE))) .sp_stop("NotImplementedError", "")
  if (all(vapply(ps, .sp_poly_rational, TRUE))) {
    G <- groebner(F, gens, order = "lex")
    return(lapply(G, identity))
  }
  .sp_groebner_field(ps, gens)
}

# Reduced lex Groebner basis over the field of rational functions in the
# parameters (sympy.polys.groebnertools over QQ(params)), each element then
# cleared of denominators as GroebnerBasis does for a ZZ[params] domain.
.sgb_cmp <- function(a, b) {
  d <- a - b
  nz <- which(d != 0L)
  if (!length(nz)) 0L else if (d[nz[1L]] > 0L) 1L else -1L
}
.sgb_norm <- function(terms) {
  keys <- vapply(terms, function(t) paste(t$monom, collapse = ","), "")
  out <- list()
  for (k in unique(keys)) {
    c <- .sps_fe(.add(lapply(terms[keys == k], `[[`, "coeff")))
    if (!is_zero_num(c)) out[[length(out) + 1L]] <- list(monom = terms[[match(k, keys)]]$monom, coeff = c)
  }
  if (length(out) > 1L) {
    o <- seq_along(out)
    for (i in 2:length(out)) {
      j <- i
      while (j > 1L && .sgb_cmp(out[[o[j]]]$monom, out[[o[j - 1L]]]$monom) > 0L) {
        tmp <- o[j]; o[j] <- o[j - 1L]; o[j - 1L] <- tmp
        j <- j - 1L
      }
    }
    out <- out[o]
  }
  out
}
.sgb_monic <- function(p) {
  lc <- p[[1L]]$coeff
  if (is_one_num(lc)) return(p)
  lapply(p, function(t) list(monom = t$monom, coeff = .sps_fe_div(t$coeff, lc)))
}
.sgb_mul_term <- function(p, m, c) lapply(p, function(t) list(monom = t$monom + m, coeff = .sps_fe_mul(t$coeff, c)))
.sgb_sub <- function(p, q) .sgb_norm(c(p, lapply(q, function(t) list(monom = t$monom, coeff = .neg(t$coeff)))))
.sgb_reduce <- function(p, G) {
  r <- list()
  while (length(p)) {
    lt <- p[[1L]]
    hit <- FALSE
    for (g in G) {
      if (all(g[[1L]]$monom <= lt$monom)) {
        p <- .sgb_sub(p, .sgb_mul_term(g, lt$monom - g[[1L]]$monom, .sps_fe_div(lt$coeff, g[[1L]]$coeff)))
        hit <- TRUE
        break
      }
    }
    if (!hit) {
      r[[length(r) + 1L]] <- lt
      p <- p[-1L]
    }
  }
  r
}
.sp_groebner_field <- function(ps, gens) {
  n <- length(gens)
  F <- Filter(length, lapply(ps, function(p) .sgb_norm(p$terms)))
  F <- lapply(F, .sgb_monic)
  if (!length(F)) return(list())
  if (any(vapply(F, function(p) all(p[[1L]]$monom == 0L), TRUE))) return(list(S_One))
  G <- F
  pairs <- list()
  for (k in seq_along(G)) for (i in seq_len(k - 1L)) pairs[[length(pairs) + 1L]] <- c(i, k)
  while (length(pairs)) {
    pr <- pairs[[1L]]
    pairs <- pairs[-1L]
    f <- G[[pr[1L]]]; g <- G[[pr[2L]]]
    a <- f[[1L]]$monom; b <- g[[1L]]$monom
    l <- pmax(a, b)
    if (all(l == a + b)) next
    sp <- .sgb_sub(.sgb_mul_term(f, l - a, S_One), .sgb_mul_term(g, l - b, S_One))
    h <- .sgb_reduce(sp, G)
    if (!length(h)) next
    h <- .sgb_monic(h)
    if (all(h[[1L]]$monom == 0L)) return(list(S_One))
    G[[length(G) + 1L]] <- h
    k <- length(G)
    for (i in seq_len(k - 1L)) pairs[[length(pairs) + 1L]] <- c(i, k)
  }
  keep <- rep(TRUE, length(G))
  for (i in seq_along(G)) for (j in seq_along(G)) {
    if (i != j && keep[j] && keep[i] && all(G[[j]][[1L]]$monom <= G[[i]][[1L]]$monom) &&
        (any(G[[j]][[1L]]$monom != G[[i]][[1L]]$monom) || j < i)) keep[i] <- FALSE
  }
  G <- G[keep]
  G <- lapply(seq_along(G), function(i) .sgb_monic(c(G[[i]][1L], .sgb_reduce(G[[i]][-1L], G[-i]))))
  o <- seq_along(G)
  if (length(G) > 1L) for (i in 2:length(G)) {
    j <- i
    while (j > 1L && .sgb_cmp(G[[o[j]]][[1L]]$monom, G[[o[j - 1L]]][[1L]]$monom) > 0L) {
      tmp <- o[j]; o[j] <- o[j - 1L]; o[j - 1L] <- tmp
      j <- j - 1L
    }
  }
  lapply(G[o], function(p) .sgb_clear_denoms(p, gens))
}
.sgb_clear_denoms <- function(p, gens) {
  dens <- lapply(p, function(t) .as_numer_denom(t$coeff)[[2L]])
  L <- S_One
  for (d in dens) if (!is_one_num(d)) L <- .sgb_lcm(L, d)
  .add(lapply(p, function(t) .mul(c(list(cancel(.mul(list(t$coeff, L)))),
                                    Map(function(g, k) .pow(g, .int(k)), gens[t$monom > 0L], t$monom[t$monom > 0L])))))
}
# lcm in ZZ[params] keeping the integer content, positive leading coefficient
.sgb_lcm <- function(a, b) {
  ca <- as_content_primitive(a); cb <- as_content_primitive(b)
  ic <- ilcm(.int_val(ca[[1L]]), .int_val(cb[[1L]]))
  pl <- lcm(ca[[2L]], cb[[2L]])
  pc <- as_content_primitive(pl)
  pl <- pc[[2L]]
  if (could_extract_minus_sign(pl)) pl <- .neg(pl)
  expand(.mul(list(as_sym(ic), pl)))
}

.sp_is_ground <- function(e, gens) !any(vapply(gens, function(g) .has1(e, g), TRUE))

.sp_solve_biquadratic <- function(f, g, gens) {
  G <- .sp_groebner(list(f, g), gens)
  if (length(G) == 1L && .sp_is_ground(G[[1L]], gens)) return(NULL)
  if (length(G) != 2L) stop(.sympy_condition("SolveFailed", ""))
  x <- gens[[1L]]; y <- gens[[2L]]
  p <- G[[1L]]; q <- G[[2L]]
  if (!.sp_is_ground(gcd(p, q), gens)) stop(.sympy_condition("SolveFailed", ""))
  pr <- .roots_sympy(p, x)
  p_roots <- lapply(pr[seq_along(pr)], function(e) rcollect(e, y))
  qr <- .roots_sympy(q, y)
  q_roots <- qr[seq_along(qr)]
  sols <- list()
  for (qv in q_roots) for (pv in p_roots) sols[[length(sols) + 1L]] <- list(subs(pv, y, qv), qv)
  .key_sort(lapply(sols, function(s) list(length(s), lapply(s, .sort_key))), sols)
}

.sp_solve_generic <- function(polys, gens, strict = FALSE) {
  is_univariate <- function(f, gs) {
    if (length(gs) <= 1L) return(TRUE)
    !any(vapply(gs[-length(gs)], function(g) .has1(f, g), TRUE))
  }
  subs_root <- function(f, gen, zero) {
    p <- subs(f, gen, zero)
    if (degree(f, gen) >= 2) p <- expand(p)
    p
  }
  rec <- function(system, gs, entry = FALSE) {
    if (length(system) == 1L && length(gs) == 1L) {
      zr <- .roots_sympy(system[[1L]], gs[[1L]], strict = strict)
      return(lapply(zr[seq_along(zr)], function(z) list(z)))
    }
    basis <- .sp_groebner(system, gs)
    if (length(basis) == 1L && .sp_is_ground(basis[[1L]], gs)) return(if (entry) NULL else list())
    univariate <- Filter(function(f) is_univariate(f, gs), basis)
    if (length(basis) < length(gs) || length(univariate) != 1L)
      .sp_stop("NotImplementedError", .filldedent("only zero-dimensional systems supported (finite number of solutions)"))
    f <- univariate[[1L]]
    gen <- gs[[length(gs)]]
    zr <- .roots_sympy(f, gen, strict = strict)
    zeros <- zr[seq_along(zr)]
    if (!length(zeros)) return(list())
    if (length(basis) == 1L) return(lapply(zeros, function(z) list(z)))
    out <- list()
    for (zero in zeros) {
      new_system <- list()
      new_gens <- gs[-length(gs)]
      for (b in basis[-length(basis)]) {
        eq <- subs_root(b, gen, zero)
        if (!is_zero_num(eq)) new_system[[length(new_system) + 1L]] <- eq
      }
      for (s in rec(new_system, new_gens)) out[[length(out) + 1L]] <- c(s, list(zero))
    }
    if (length(out) && length(out[[1L]]) != length(gs))
      .sp_stop("NotImplementedError", .filldedent("only zero-dimensional systems supported (finite number of solutions)"))
    out
  }
  result <- rec(polys, gens, entry = TRUE)
  if (is.null(result)) return(NULL)
  .key_sort(lapply(result, function(s) list(length(s), lapply(s, .sort_key))), result)
}

# ---------------------------------------------------------------------------
# _solve_system

.sp_solve_system <- function(exprs, symbols, flags) {
  if (!length(exprs)) return(list(FALSE, list()))
  if (!length(symbols)) {
    u <- Dummy()
    if (all(vapply(exprs, function(e) !.sp_bool(e) || isTRUE(.sp_checksol(u, .sd1(u, e))), TRUE)))
      return(list(TRUE, list(.sd_new())))
    return(list(TRUE, list()))
  }
  split_ <- .fl(flags, "_split", TRUE)
  flags[["_split"]] <- NULL
  if (!isFALSE(split_)) {
    exprsyms <- lapply(exprs, function(e) Filter(function(s) .has_free(e, list(s)), symbols))
    n <- length(exprs)
    comp <- seq_len(n)
    findc <- function(i) { while (comp[i] != i) i <- comp[i]; i }
    for (a in seq_len(n)) for (b in seq_len(a - 1L)) {
      if (any(vapply(exprsyms[[a]], function(s) .sp_in(s, exprsyms[[b]]), TRUE))) {
        ra <- findc(a); rb <- findc(b)
        if (ra != rb) comp[max(ra, rb)] <- min(ra, rb)
      }
    }
    roots_ <- vapply(seq_len(n), findc, 0L)
    if (length(unique(roots_)) > 1L) {
      subsols <- list()
      linear <- TRUE
      for (cr in unique(roots_)) {
        idx <- which(roots_ == cr)
        subsyms <- unique_syms(exprsyms[idx])
        subsyms <- subsyms[order(vapply(subsyms, function(s) match(s$key, vapply(symbols, .k, "")), 0L))]
        fl <- flags
        fl[["_split"]] <- FALSE
        r <- .sp_solve_system(exprs[idx], subsyms, fl)
        if (linear) linear <- isTRUE(r[[1L]])
        subsols[[length(subsols) + 1L]] <- r[[2L]]
      }
      sols <- list(.sd_new())
      for (ss in subsols) {
        new <- list()
        for (a in sols) for (b in ss) {
          d <- a
          for (j in seq_along(b$k)) d <- .sd_set(d, b$k[[j]], b$v[[j]])
          new[[length(new) + 1L]] <- d
        }
        sols <- new
      }
      return(list(linear, sols))
    }
  }
  polys <- list()
  dens <- list()
  failed <- list()
  result <- list()
  solved_syms <- list()
  linear <- TRUE
  manual <- isTRUE(.fl(flags, "manual", FALSE))
  check <- !isFALSE(.fl(flags, "check", TRUE))
  checkdens <- check
  for (g in exprs) {
    dens <- .sp_uniq(c(dens, .sp_simple_dens(g, symbols)))
    iv <- .sp_invert(g, symbols)
    i <- iv[[1L]]; d <- iv[[2L]]
    if (.sp_in(d, symbols)) {
      if (linear) linear <- identical(.sp_solve_linear(g, S_Zero, list(d))[[1L]]$key, d$key)
    }
    g <- .sub(d, i)
    g <- .as_numer_denom(g)[[1L]]
    if (manual) {
      failed[[length(failed) + 1L]] <- g
      next
    }
    p <- .sp_as_poly(g, symbols)
    if (!is.null(p)) polys[[length(polys) + 1L]] <- p else failed[[length(failed) + 1L]] <- g
  }
  if (length(polys)) {
    if (all(vapply(polys, .sp_poly_is_linear, TRUE))) {
      if (isTRUE(.fl(flags, "particular", FALSE))) {
        flags$particular <- NULL
        res <- .sp_minsolve_linear(polys, symbols, flags)
      } else res <- .sp_solve_lin_sys(polys, symbols)
      result <- if (!is.null(res) && length(res$k)) list(res) else list()
      if (length(failed)) solved_syms <- if (length(result)) result[[1L]]$k else list()
    } else {
      linear <- FALSE
      if (length(symbols) > length(polys)) {
        free <- unique_syms(lapply(polys, .sp_poly_free))
        free <- .sp_ordered(Filter(function(s) .sp_in(s, symbols), free))
        got_s <- list()
        result <- list()
        for (syms in combinations(free, min(length(free), length(polys)))) {
          res <- tryCatch(.sp_solve_poly_system(lapply(polys, function(p) .sp_as_poly(.sp_poly_expr(p), syms)), syms),
                          NotImplementedError = function(e) NULL)
          if (length(res)) {
            for (r in .sp_unique_tuples(res)) {
              skip <- FALSE
              for (r1 in r) if (length(got_s) && any(vapply(got_s, function(ss) .has_free(r1, list(ss)), TRUE))) skip <- TRUE
              if (!skip) {
                for (s in syms) if (!.sp_in(s, got_s)) got_s[[length(got_s) + 1L]] <- s
                result[[length(result) + 1L]] <- list(k = syms, v = r)
              }
            }
          }
        }
        if (length(got_s)) solved_syms <- got_s else failed <- c(failed, lapply(polys, .sp_poly_expr))
      } else {
        res <- tryCatch(.sp_solve_poly_system(polys, symbols), NotImplementedError = function(e) "failed")
        if (identical(res, "failed")) {
          failed <- c(failed, lapply(polys, .sp_poly_expr))
          solved_syms <- list()
        } else if (length(res)) {
          solved_syms <- symbols
          result <- lapply(.sp_unique_tuples(res), function(r) list(k = solved_syms, v = r))
        }
      }
    }
  }
  if (!length(result)) result <- list(.sd_new())
  if (length(failed)) {
    linear <- FALSE
    legal <- symbols
    u <- Dummy()
    ok_syms <- function(e, sort = FALSE) {
      rv <- Filter(function(s) .sp_in(s, legal), free_symbols(e))
      if (sort && length(rv) > 1L) {
        keyf <- function(sym) {
          ep <- .sp_as_poly(e, list(sym))
          cx <- if (is.null(ep)) list(Inf, Inf, Inf) else {
            lc <- .sp_poly_lc(ep)
            cs <- free_symbols(lc)
            list(max(vapply(ep$terms, function(t) t$monom[1L], 0L)), sum(vapply(cs, function(z) .sp_in(z, rv), TRUE)), length(cs))
          }
          c(cx, list(.sort_key(sym)))
        }
        rv <- .key_sort(lapply(rv, keyf), rv)
      }
      rv
    }
    ord <- order(vapply(failed, function(e) length(ok_syms(e)), 0L))
    # ordered(failed, key) : ties broken by (_nodes, default_sort_key)
    failed <- .key_sort(lapply(failed, function(e) list(length(ok_syms(e)), .nodes(e), .sort_key(e))), failed)
    for (eq in failed) {
      newresult <- list()
      bad_results <- list()
      hit <- FALSE
      broke <- FALSE
      for (r in result) {
        got_s <- list()
        eq2 <- .sd_subs(eq, r)
        if (check && length(r$k)) {
          b <- .sp_checksol(u, .sd1(u, eq2), list(minimal = TRUE))
          if (!is.na(b)) {
            if (b) newresult[[length(newresult) + 1L]] <- r else bad_results[[length(bad_results) + 1L]] <- r
            next
          }
        }
        oks <- ok_syms(eq2, sort = TRUE)
        if (!length(oks)) {
          if (length(r$k)) newresult[[length(newresult) + 1L]] <- r
          broke <- TRUE
          break
        }
        for (s in oks) {
          soln <- tryCatch(.sp_vsolve(eq2, s, flags), NotImplementedError = function(e) NULL)
          if (is.null(soln)) next
          for (sol in soln) {
            if (length(got_s) && any(vapply(got_s, function(ss) .has_free(sol, list(ss)), TRUE))) next
            rnew <- r
            for (j in seq_along(r$k)) rnew$v[[j]] <- subs(r$v[[j]], s, sol)
            rnew <- .sd_set(rnew, s, sol)
            iset <- paste(vapply(rnew$k, .k, ""), vapply(rnew$v, .k, ""))
            dup <- FALSE
            for (i in newresult) {
              if (length(i$k) < length(rnew$k)) {
                upd <- paste(vapply(i$k, .k, ""), vapply(i$v, function(v) xreplace(v, Map(list, rnew$k, rnew$v))$key, ""))
                if (all(upd %in% iset)) { dup <- TRUE; break }
              }
            }
            if (!dup) newresult[[length(newresult) + 1L]] <- rnew
          }
          hit <- TRUE
          got_s[[length(got_s) + 1L]] <- s
        }
        if (!hit) .sp_stop("NotImplementedError", paste0("could not solve ", .sp_py_str(eq2)))
      }
      if (!broke) {
        result <- newresult
        bk <- vapply(bad_results, .sd_ident, "")
        if (length(bk)) result <- Filter(function(r) !(.sd_ident(r) %in% bk), result)
      }
    }
  }
  if (!length(result)) return(list(FALSE, list()))
  default_simplify <- length(failed) > 0L
  if (isTRUE(.fl(flags, "simplify", default_simplify))) {
    result <- lapply(result, function(r) .sd_map(r, simplify))
    flags$simplify <- FALSE
  }
  if (checkdens) result <- Filter(function(r) !any(vapply(dens, function(d) isTRUE(.sp_checksol(d, r, flags)), TRUE)), result)
  if (check && !linear) result <- Filter(function(r) !any(vapply(exprs, function(e) isFALSE(.sp_checksol(e, r, flags)), TRUE)), result)
  result <- Filter(function(r) length(r$k) > 0L, result)
  list(linear, result)
}

.sp_unique_tuples <- function(res) res[!duplicated(vapply(res, function(r) paste(vapply(r, .k, ""), collapse = "|"), ""))]

.sp_poly_lc <- function(p) {
  # leading coefficient in lex order of the generators
  best <- NULL
  for (t in p$terms) {
    if (is.null(best)) { best <- t; next }
    d <- t$monom - best$monom
    nz <- which(d != 0L)
    if (length(nz) && d[nz[1L]] > 0L) best <- t
  }
  if (is.null(best)) S_Zero else best$coeff
}

.sp_minsolve_linear <- function(polys, symbols, flags) {
  .sp_stop("NotImplementedError", "particular=True is not supported")
}

# ---------------------------------------------------------------------------
# unrad: remove radicals by symbolic powering (solvers.py unrad)

.sp_mexpand_rec <- function(e) {
  repeat {
    was <- e
    e <- .mexpand(e)
    if (identical(was$key, e$key)) return(e)
  }
}

.sp_sqrt_depth <- function(p) {
  if (identical(p$key, S_I$key)) return(1L)
  if (!length(p$args)) return(0L)
  if (p$type %in% c("Add", "Mul")) return(max(vapply(p$args, .sp_sqrt_depth, 0L)))
  if (p$type == "Pow" && is_exact_number(p$args[[2L]]) && abs(as.numeric(p$args[[2L]]$value)) == 0.5)
    return(.sp_sqrt_depth(p$args[[1L]]) + 1L)
  0L
}

.sp_unrad_Q <- function(pw) {
  ex <- .as_base_exp(pw)[[2L]]
  c <- as_coeff_mul_split(ex)[[1L]]
  if (!is_exact_number(c)) return(gmp::as.bigz(1))
  gmp::denominator(gmp::as.bigq(c$value))
}

.sp_unrad_take <- function(d, syms) {
  for (pw in .margs(d)) {
    if (pw$type != "Pow") next
    if (.sp_unrad_Q(pw) == 1) next
    if (.has_free(pw, syms)) return(TRUE)
  }
  FALSE
}

.sp_unrad_rbl <- function(gens) {
  lcm <- gmp::as.bigz(1)
  rads <- list()
  bases <- list()
  for (g in gens) {
    q <- .sp_unrad_Q(g)
    if (q != 1) {
      rads[[g$key]] <- g
      lcm <- gmp::lcm.bigz(lcm, q)
      b <- .as_base_exp(g)[[1L]]
      bases[[b$key]] <- b
    }
  }
  list(rads = unname(rads), bases = unname(bases), lcm = as.integer(as.character(lcm)))
}

.sp_unrad_canonical <- function(eq, cov) {
  if (length(cov)) {
    p <- cov[[1L]]
    d <- Dummy(p$name)
    rep <- list(list(p, d))
    eq <- xreplace(eq, rep)
    cov <- list(d, xreplace(cov[[2L]], rep))
  }
  eq <- factor_terms(.sp_mexpand_rec(.as_numer_denom(eq)[[1L]]), clear = TRUE)
  if (eq$type == "Mul") {
    args <- list()
    for (f in eq$args) {
      if (.sp_is_number(f)) next
      args[[length(args) + 1L]] <- if (f$type == "Pow") f$args[[1L]] else f
    }
    eq <- .mul(args)
  }
  margs <- .margs(eq)
  changed <- FALSE
  for (i in seq_along(margs)) {
    if (could_extract_minus_sign(margs[[i]])) {
      margs[[i]] <- .neg(margs[[i]])
      changed <- TRUE
    }
  }
  if (changed) eq <- .mul(margs)
  list(eq, cov)
}

.sp_unrad <- function(eq, syms = list(), flags = list()) {
  uflags <- list(check = FALSE, simplify = FALSE)
  if (is_sym(eq) && eq$type == "Equality") eq <- .sub(eq$args[[1L]], eq$args[[2L]])
  else if (!.is_expr_node(eq)) return(NULL)
  cov <- .fl(flags, "cov", list())
  nwas <- flags$n
  rpt <- .fl(flags, "rpt", 0L)
  covf <- function(p, e) {
    if (length(cov)) {
      oldp <- cov[[1L]]; olde <- cov[[2L]]
      if (degree(e, p) %in% c(1, 2)) {
        cov <<- list(p, subs(olde, oldp, .sp_vsolve(e, p, uflags)[[1L]]))
      } else .sp_stop("NotImplementedError", "")
    } else cov <<- list(p, e)
  }
  eq <- powdenest(factor_terms(eq, radical = TRUE, clear = TRUE))
  eq <- .as_numer_denom(eq)[[1L]]
  eq <- .sp_mexpand_rec(eq)
  if (.sp_is_number(eq)) return(NULL)
  syms <- if (length(syms)) syms else free_symbols(eq)
  pg <- .sp_poly_gens(eq)
  gens <- Filter(function(g) .sp_unrad_take(g, syms), pg)
  if (!length(gens)) return(NULL)
  poly <- .sp_as_poly(eq, gens)
  if (is.null(poly)) return(NULL)
  if (any(vapply(gens, function(g) .has_free(.as_base_exp(g)[[2L]], syms), TRUE))) return(NULL)
  rbl <- .sp_unrad_rbl(gens)
  rads <- rbl$rads; bases <- rbl$bases; lcm <- rbl$lcm
  covsym <- Dummy("p", nonnegative = TRUE)
  newsyms <- list()
  for (r in rads) for (s in syms) if (.has_free(r, list(s)) && !.sp_in(s, newsyms)) newsyms[[length(newsyms) + 1L]] <- s
  if (length(newsyms) != length(syms)) syms <- newsyms
  rk <- vapply(rads, .k, "")
  groups <- list(list(key = "", terms = list()))
  gkeys <- ""
  for (t in .aargs(.sp_poly_expr(poly))) {
    key <- ""
    if (.sp_unrad_take(t, syms)) {
      tg <- .sp_poly_gens(t)
      common <- sort(match(intersect(vapply(tg, .k, ""), rk), rk))
      key <- paste(common, collapse = ",")
    }
    i <- match(key, gkeys)
    if (is.na(i)) {
      gkeys <- c(gkeys, key)
      groups[[length(groups) + 1L]] <- list(key = key, terms = list(t))
    } else groups[[i]]$terms[[length(groups[[i]]$terms) + 1L]] <- t
  }
  others <- .add(groups[[1L]]$terms)
  rterms <- lapply(groups[-1L], function(g) .add(g$terms))
  rterms <- rev(.sp_ordered(rterms))
  ok <- FALSE
  depth <- .sp_sqrt_depth(eq)
  is_pow2 <- function(n) n > 0 && bitwAnd(n, n - 1L) == 0L
  if (length(rterms) == 1L && !(rterms[[1L]]$type == "Add" && lcm > 2)) {
    eq <- .sub(.pow(rterms[[1L]], .int(lcm)), .pow(.neg(others), .int(lcm)))
    ok <- TRUE
  } else {
    if (length(rterms) == 1L && rterms[[1L]]$type == "Add") rterms <- .sympy_args(rterms[[1L]])
    if (length(bases) == 1L) {
      b <- bases[[1L]]
      x <- if (length(syms) > 1L) free_symbols(b) else syms
      x <- .sp_ordered(x)[[1L]]
      r <- tryCatch({
        inv <- .sp_vsolve(.sub(.pow(covsym, .int(lcm)), b), x, uflags)
        if (!length(inv)) .sp_stop("NotImplementedError", "")
        eq2 <- subs(subs(.sp_poly_expr(poly), b, .pow(covsym, .int(lcm))), x, inv[[1L]])
        covf(covsym, .sub(.pow(covsym, .int(lcm)), b))
        .sp_unrad_canonical(eq2, cov)
      }, NotImplementedError = function(e) NULL)
      if (!is.null(r)) return(r)
    }
    if (length(rterms) == 2L) {
      if (!.sp_bool(others)) {
        eq <- .sub(.pow(rterms[[1L]], .int(lcm)), .pow(.neg(rterms[[2L]]), .int(lcm)))
        ok <- TRUE
      } else if (!is_pow2(lcm)) {
        r0 <- rterms[[1L]]; r1 <- rterms[[2L]]
        if (isTRUE(flags[["_reverse"]])) { tmp <- r0; r0 <- r1; r1 <- tmp }
        i0 <- .sp_unrad_rbl(.sp_poly_gens(r0))
        i1 <- .sp_unrad_rbl(.sp_poly_gens(r1))
        for (reverse in 0:1) {
          if (reverse) {
            tmp <- i0; i0 <- i1; i1 <- tmp
            tmp <- r0; r0 <- r1; r1 <- tmp
          }
          rads1 <- .mul(i1$rads)
          lcm1 <- i1$lcm
          t1 <- .pow(rads1, .int(lcm1))
          c <- .sub(.pow(covsym, .int(lcm1)), t1)
          for (x in syms) {
            done <- tryCatch({
              sol <- .sp_vsolve(c, x, uflags)
              if (!length(sol)) .sp_stop("NotImplementedError", "")
              neweq <- .add(list(subs(r0, x, sol[[1L]]), .div(.mul(list(covsym, r1)), rads1), others))
              tmp <- .sp_unrad(neweq, list(covsym))
              if (!is.null(tmp)) {
                eq <- tmp[[1L]]
                newcov <- tmp[[2L]]
                if (length(newcov)) {
                  covf(newcov[[1L]], subs(c, covsym, .sp_vsolve(newcov[[2L]], covsym, uflags)[[1L]]))
                } else covf(covsym, c)
              } else {
                eq <- neweq
                covf(covsym, c)
              }
              TRUE
            }, NotImplementedError = function(e) {
              if (reverse) .sp_stop("NotImplementedError", "no successful change of variable found")
              FALSE
            })
            if (done) { ok <- TRUE; break }
          }
          if (ok) break
        }
      }
    } else if (length(rterms) == 3L) {
      info <- lapply(rterms, function(i) .sp_unrad_rbl(.sp_poly_gens(i)))
      if (info[[1L]]$lcm != 3) {
        info <- c(info[-1L], info[1L]); rterms <- c(rterms[-1L], rterms[1L])
      } else if (info[[2L]]$lcm != 3) {
        info <- c(info[-2L], info[2L]); rterms <- c(rterms[-2L], rterms[2L])
      }
      bkey <- function(i) paste(sort(vapply(i$bases, .k, "")), collapse = "|")
      if (info[[1L]]$lcm == 3 && info[[2L]]$lcm == 3) {
        if (bkey(info[[2L]]) != bkey(info[[3L]])) {
          tmp <- info[[1L]]; info[[1L]] <- info[[2L]]; info[[2L]] <- tmp
          tmp <- rterms[[1L]]; rterms[[1L]] <- rterms[[2L]]; rterms[[2L]] <- tmp
        }
        if (bkey(info[[2L]]) == bkey(info[[3L]])) {
          eq <- .add(list(.pow(rterms[[1L]], .int(3L)), .pow(.add(list(rterms[[2L]], rterms[[3L]], others)), .int(3L))))
          ok <- TRUE
        } else if (info[[3L]]$lcm < 5) {
          tt <- function(i) list(cancel(.div(rterms[[i]], .mul(info[[i]]$rads))), .mul(info[[i]]$bases))
          A1 <- tt(1L); B1 <- tt(2L)
          eq <- .sp_unrad_zz(A1[[1L]], A1[[2L]], B1[[1L]], B1[[2L]], .neg(rterms[[3L]]), others)
          ok <- TRUE
        }
      }
    }
    if (!ok) {
      if (is_pow2(lcm) && ((!.sp_bool(others) && length(rterms) == 4L) || length(rterms) < 4L)) {
        norm2 <- function(a, b) .add(list(.pow(a, S_Two), .pow(b, S_Two), .mul(list(S_Two, a, b))))
        if (length(rterms) == 4L) {
          eq <- .sub(norm2(rterms[[1L]], rterms[[2L]]), norm2(rterms[[3L]], rterms[[4L]]))
          ok <- TRUE
        } else if (length(rterms) == 3L) {
          eq <- .sub(norm2(rterms[[2L]], rterms[[3L]]), norm2(rterms[[1L]], others))
          ok <- TRUE
        } else if (length(rterms) == 2L) {
          eq <- .sub(.pow(rterms[[1L]], S_Two), norm2(rterms[[2L]], others))
          ok <- TRUE
        }
      }
    }
  }
  new_depth <- if (ok) .sp_sqrt_depth(eq) else depth
  rpt <- rpt + 1L
  if (!ok || (!is.null(nwas) && length(rterms) == nwas && new_depth == depth && rpt > 3L))
    .sp_stop("NotImplementedError", "Cannot remove all radicals")
  fl2 <- flags
  fl2$cov <- cov
  fl2$n <- length(rterms)
  fl2$rpt <- rpt
  neq <- .sp_unrad(eq, syms, fl2)
  if (!is.null(neq)) {
    eq <- neq[[1L]]
    cov <- neq[[2L]]
  }
  .sp_unrad_canonical(eq, cov)
}

.sp_unrad_zz <- function(aa, AA, bb, BB, cc, dd) {
  zz <- sympify(paste0("(c - d)*(A**3*a**9 + 3*A**2*B*a**6*b**3 - 3*A**2*a**6*c**3 + 9*A**2*a**6*c**2*d - 9*A**2*a**6*c*d**2 + ",
    "3*A**2*a**6*d**3 + 3*A*B**2*a**3*b**6 + 21*A*B*a**3*b**3*c**3 - 63*A*B*a**3*b**3*c**2*d + 63*A*B*a**3*b**3*c*d**2 - ",
    "21*A*B*a**3*b**3*d**3 + 3*A*a**3*c**6 - 18*A*a**3*c**5*d + 45*A*a**3*c**4*d**2 - 60*A*a**3*c**3*d**3 + 45*A*a**3*c**2*d**4 - ",
    "18*A*a**3*c*d**5 + 3*A*a**3*d**6 + B**3*b**9 - 3*B**2*b**6*c**3 + 9*B**2*b**6*c**2*d - 9*B**2*b**6*c*d**2 + 3*B**2*b**6*d**3 + ",
    "3*B*b**3*c**6 - 18*B*b**3*c**5*d + 45*B*b**3*c**4*d**2 - 60*B*b**3*c**3*d**3 + 45*B*b**3*c**2*d**4 - 18*B*b**3*c*d**5 + ",
    "3*B*b**3*d**6 - c**9 + 9*c**8*d - 36*c**7*d**2 + 84*c**6*d**3 - 126*c**5*d**4 + 126*c**4*d**5 - 84*c**3*d**6 + 36*c**2*d**7 - ",
    "9*c*d**8 + d**9)"))
  xreplace(zz, unname(Map(function(n, v) list(Symbol(n), v), c("a", "A", "b", "B", "c", "d"), list(aa, AA, bb, BB, cc, dd))))
}

# ---------------------------------------------------------------------------
# bivariate.py: _solve_lambert and bivariate_type, used by _tsolve

.sp_filtered_gens <- function(gens, symbol) {
  gs <- Filter(function(g) .has_free(g, list(symbol)), gens)
  out <- gs
  for (g in gs) {
    ag <- .inv(g)
    if (.sp_in(g, out) && .sp_in(ag, out)) {
      if (!is_one_num(.as_numer_denom(ag)[[2L]])) g <- ag
      out <- Filter(function(z) !identical(z$key, g$key), out)
    }
  }
  out
}

.sp_count_type <- function(e, pred) {
  n <- 0L
  walk <- function(z) {
    if (pred(z)) n <<- n + 1L
    for (a in z$args) walk(a)
  }
  walk(e)
  n
}

.sp_fpred <- function(func) switch(func,
  log = function(z) z$type == "Function" && z$name == "log",
  exp = function(z) z$type == "Function" && z$name == "exp",
  Pow = function(z) z$type == "Pow")

.sp_mostfunc <- function(lhs, func, X = NULL) {
  pred <- .sp_fpred(func)
  fterms <- .sp_find(lhs, pred)
  if (!is.null(X)) fterms <- Filter(function(t) .has_free(t, list(X)), fterms)
  if (length(fterms) == 1L) return(fterms[[1L]])
  if (length(fterms)) {
    o <- .sp_ordered(fterms)
    cnt <- vapply(o, function(t) .sp_count_type(t, pred), 0L)
    return(o[[which.max(cnt)]])
  }
  NULL
}

.sp_linab <- function(arg, symbol) {
  arg <- factor_terms(expand(arg))
  id <- .sp_indep(arg, list(symbol))
  ind <- id[[1L]]; dep <- id[[2L]]
  if (arg$type == "Mul" && dep$type == "Add") {
    r <- .sp_linab(dep, symbol)
    return(list(.mul(list(ind, r[[1L]])), .mul(list(ind, r[[2L]])), r[[3L]]))
  }
  if (arg$type != "Add") {
    b <- S_Zero
    a <- ind; x <- dep
  } else {
    b <- ind
    id2 <- .sp_indep(separatevars(dep), list(symbol), as_add = FALSE)
    a <- id2[[1L]]; x <- id2[[2L]]
  }
  if (could_extract_minus_sign(x)) {
    a <- .neg(a); x <- .neg(x)
  }
  list(a, b, x)
}

.sp_lambert <- function(eq, x) {
  eq <- .mexpand(expand_log(eq))
  mainlog <- .sp_mostfunc(eq, "log", x)
  if (is.null(mainlog)) return(list())
  other <- subs(eq, mainlog, S_Zero)
  if (.neg(other)$type == "Function" && .neg(other)$name == "log") {
    eq <- subs(.sub(eq, other), mainlog, mainlog$args[[1L]])
    mainlog <- mainlog$args[[1L]]
    if (!(mainlog$type == "Function" && mainlog$name == "log")) return(list())
    other <- .neg(.neg(other)$args[[1L]])
    eq <- .add(list(eq, other))
  }
  if (!.sp_in(x, free_symbols(other))) return(list())
  dfx <- .sp_linab(other, x)
  d <- dfx[[1L]]; f <- dfx[[2L]]; X2 <- dfx[[3L]]
  logterm <- collect(.sub(eq, other), list(mainlog))
  a <- as_coefficient(logterm, mainlog)
  if (is.null(a) || .sp_in(x, free_symbols(a))) return(list())
  logarg <- mainlog$args[[1L]]
  bcx <- .sp_linab(logarg, x)
  b <- bcx[[1L]]; c <- bcx[[2L]]; X1 <- bcx[[3L]]
  if (!identical(X1$key, X2$key)) return(list())
  u <- Dummy("rhs")
  xusolns <- solve(.sub(X1, u), x)
  sol <- list()
  nd <- .as_numer_denom(.div(.div(.sub(.mul(list(c, d)), .mul(list(b, f))), a), b))
  num <- nd[[1L]]
  pd <- as_coeff_mul_split(nd[[2L]])
  p <- pd[[1L]]; den <- pd[[2L]]
  e <- fcall("exp", list(.div(num, den)))
  t <- Dummy("t")
  rr <- .roots_sympy(.sub(.pow(t, p), e), t)
  args <- lapply(rr[seq_along(rr)], function(r) .mul(list(.div(d, .mul(list(a, b))), r)))
  for (arg in args) {
    for (k in c(-1L, 0L)) {
      w <- LambertW(arg, k)
      if (k != 0L && !isTRUE(.sps_lambertw_real(arg, k))) next
      rhs <- .add(list(.div(.neg(c), b), .mul(list(.div(a, d), w))))
      for (xu in xusolns) sol[[length(sol) + 1L]] <- subs(xu, u, rhs)
    }
  }
  sol
}

# LambertW(x, k).is_extended_real as SymPy decides it
.sps_lambertw_real <- function(x, k) {
  xe <- .add(list(x, .inv(S_E)))
  if (k == 0L) {
    if (isTRUE(.sp_fact(xe, "positive"))) return(TRUE)
    if (isTRUE(.sp_fact(xe, "nonpositive"))) return(FALSE)
  } else if (k == -1L) {
    if (isTRUE(.sp_fact(x, "negative")) && isTRUE(.sp_fact(xe, "positive"))) return(TRUE)
    if (isTRUE(.sp_fact(x, "nonpositive")) || isTRUE(.sp_fact(xe, "nonnegative"))) return(FALSE)
  } else if (isTRUE(.sp_fact(x, "extended_real"))) return(FALSE)
  NA
}

# nsimplify(e) as SymPy for an exact number expression: its numeric value
# identified as a simple exact number, else e unchanged
.sps_nsimplify <- function(e) {
  if (length(free_symbols(e)) || !length(e$args)) return(nsimplify(e))
  v <- tryCatch(as.complex(.neval(e, list())), error = function(err) NA_complex_)
  if (length(v) != 1L || is.na(v) || !is.finite(Re(v)) || !is.finite(Im(v))) return(e)
  ident <- function(w) {
    if (abs(w) < 1e-14) return(S_Zero)
    r <- .nsimplify_num(w, list(), 1e-10, FALSE)
    if (r$type == "Float") NULL else r
  }
  re <- ident(Re(v))
  im <- if (abs(Im(v)) > 1e-14 * max(1, abs(Re(v)))) ident(Im(v)) else S_Zero
  if (is.null(re) || is.null(im)) return(e)
  .add(list(re, .mul(list(im, S_I))))
}

# factor(e, deep=True); ponytail: the package's factor() stops on some
# expressions mixing x and log(x) that SymPy factors, so e is then kept as is
.sps_factor_deep <- function(e) {
  tryCatch(factor_sym(e, deep = TRUE), error = function(err) {
    if (grepl("time limit", conditionMessage(err))) stop(err)
    e
  })
}

.sp_is_func <- function(z, nm) z$type == "Function" && z$name %in% nm

.sp_solve_lambert <- function(f, symbol, gens) {
  even_deg <- function(expr, t) {
    nlhs <- xreplace(expr, list(list(t, .neg(symbol))))
    plhs <- xreplace(expr, list(list(t, symbol)))
    sols <- .sp_solve_lambert(nlhs, symbol, gens)
    if (!identical(plhs$key, nlhs$key)) sols <- c(sols, .sp_solve_lambert(plhs, symbol, gens))
    .sp_uniq(sols)
  }
  id <- .sp_indep(f, list(symbol), as_add = TRUE)
  rhs <- .neg(id[[1L]]); lhs <- id[[2L]]
  lamcheck <- Filter(function(tmp) .sp_is_func(tmp, c("exp", "log")) ||
                       (tmp$type == "Pow" && .sp_in(symbol, free_symbols(tmp$args[[2L]]))), gens)
  if (!length(lamcheck)) .sp_stop("NotImplementedError", "")
  if (lhs$type %in% c("Add", "Mul")) {
    t <- .sp_dummy_like(symbol, "t")
    lhs <- .map_bottom_up(lhs, function(i) {
      if (i$type == "Pow" && identical(i$args[[1L]]$key, symbol$key) && isTRUE(.sp_fact(i$args[[2L]], "even")))
        .pow(t, i$args[[2L]]) else i
    })
    if (lhs$type == "Add" && .has1(lhs, t)) {
      t_indep <- subs(lhs, t, S_Zero)
      t_term <- .sub(lhs, t_indep)
      rhs_ <- .sub(rhs, t_indep)
      if (t_term$type != "Add" && .sp_bool(rhs_) && !.sp_has_illegal(t_term)) {
        eq <- expand_log(.sub(fcall("log", list(t_term)), fcall("log", list(rhs_))))
        return(even_deg(eq, t))
      }
    } else if (lhs$type == "Mul" && .sp_bool(rhs)) {
      lhs <- expand_log(fcall("log", list(lhs)), force = TRUE)
      rhs <- fcall("log", list(rhs))
      if (.has1(lhs, t) && lhs$type == "Add") return(even_deg(.sub(lhs, rhs), t))
    }
    lhs <- xreplace(lhs, list(list(t, symbol)))
  }
  lhs <- powsimp(.sps_factor_deep(lhs))
  r <- Dummy()
  iv <- .sp_invert(.sub(lhs, r), list(symbol))
  lhs <- iv[[2L]]
  rhs <- xreplace(iv[[1L]], list(list(r, rhs)))
  soln <- list()
  mainlog <- .sp_mostfunc(lhs, "log", symbol)
  if (!is.null(mainlog)) {
    if (lhs$type == "Mul" && .sp_bool(rhs)) {
      soln <- .sp_lambert(.sub(fcall("log", list(lhs)), fcall("log", list(rhs))), symbol)
    } else if (lhs$type == "Add") {
      other <- subs(lhs, mainlog, S_Zero)
      if (.sp_bool(other) && other$type != "Add" &&
          length(Filter(function(tmp) .sp_in(symbol, free_symbols(tmp)), atoms(other, "Pow")))) {
        diff <- if (!.sp_bool(rhs)) .sub(fcall("log", list(other)), fcall("log", list(.sub(other, lhs)))) else
          .sub(fcall("log", list(.sub(lhs, other))), fcall("log", list(.sub(rhs, other))))
        soln <- .sp_lambert(expand_log(diff), symbol)
      } else soln <- .sp_lambert(.sub(lhs, rhs), symbol)
    }
  }
  if (!length(soln)) {
    mainexp <- .sp_mostfunc(lhs, "exp", symbol)
    if (!is.null(mainexp)) {
      lhs <- collect(lhs, list(mainexp))
      if (lhs$type == "Mul" && .sp_bool(rhs)) {
        soln <- .sp_lambert(expand_log(.sub(fcall("log", list(lhs)), fcall("log", list(rhs)))), symbol)
      } else if (lhs$type == "Add") {
        other <- subs(lhs, mainexp, S_Zero)
        mainterm <- .sub(lhs, other)
        rhs <- .sub(rhs, other)
        if (could_extract_minus_sign(mainterm) && could_extract_minus_sign(rhs)) {
          mainterm <- .neg(mainterm); rhs <- .neg(rhs)
        }
        soln <- .sp_lambert(expand_log(.sub(fcall("log", list(mainterm)), fcall("log", list(rhs)))), symbol)
      }
    }
  }
  if (!length(soln)) {
    mainpow <- .sp_mostfunc(lhs, "Pow", symbol)
    if (!is.null(mainpow) && .sp_in(symbol, free_symbols(mainpow$args[[2L]]))) {
      lhs <- collect(lhs, list(mainpow))
      if (lhs$type == "Mul" && .sp_bool(rhs)) {
        soln <- .sp_lambert(expand_log(.sub(fcall("log", list(lhs)), fcall("log", list(rhs)))), symbol)
      } else if (lhs$type == "Add") {
        other <- subs(lhs, mainpow, S_Zero)
        mainterm <- .sub(lhs, other)
        rhs <- .sub(rhs, other)
        soln <- .sp_lambert(expand_log(.sub(fcall("log", list(mainterm)), fcall("log", list(rhs)))), symbol)
      }
    }
  }
  if (!length(soln))
    .sp_stop("NotImplementedError", paste0(.sp_py_str(f), " does not appear to have a solution in terms of LambertW"))
  .sp_ordered(soln)
}

# Dummy carrying the assumptions of a symbol
.sp_dummy_like <- function(s, name) {
  facts <- if (is_symbol_like(s)) s$extra$facts else NULL
  if (is.null(facts) || !length(facts)) return(Dummy(name))
  facts <- facts[!is.na(facts)]
  do.call(Dummy, c(list(name), as.list(facts)))
}

.sp_bivariate_type <- function(f, x, y, first = TRUE) {
  u <- Dummy("u", positive = TRUE)
  if (first) {
    p <- .sp_as_poly(f, list(x, y))
    if (is.null(p)) .sp_stop("PolynomialError", "")
    f <- .sp_poly_expr(p)
    x_ <- Dummy(); y_ <- Dummy()
    # subs with a dict replaces the item with more nodes first
    nodes <- function(e) 1L + sum(vapply(e$args, nodes, 0L))
    reps <- if (nodes(y) > nodes(x)) list(list(y, y_), list(x, x_)) else list(list(x, x_), list(y, y_))
    rv <- .sp_bivariate_type(subs(f, reps), x_, y_, first = FALSE)
    if (!is.null(rv)) {
      reps <- list(list(x_, x), list(y_, y))
      return(list(xreplace(rv[[1L]], reps), xreplace(rv[[2L]], reps), rv[[3L]]))
    }
    return(NULL)
  }
  p <- .sp_as_poly(f, list(x, y))
  new <- list()
  okall <- TRUE
  for (a in .aargs(f)) {
    a <- .mexpand(subs(a, x, .div(u, y)))
    if (.sp_in(x, free_symbols(a)) || .sp_in(y, free_symbols(a))) { okall <- FALSE; break }
    new[[length(new) + 1L]] <- a
  }
  if (okall) return(list(.mul(list(x, y)), .add(new), u))
  okf <- function(f, v, c) {
    nw <- .mexpand(subs(f, v, c))
    fr <- free_symbols(nw)
    if (.sp_in(x, fr) || .sp_in(y, fr)) NULL else nw
  }
  coeff_mon <- function(p, mx, my) {
    for (t in p$terms) if (t$monom[1L] == mx && t$monom[2L] == my) return(t$coeff)
    S_Zero
  }
  d <- .sp_poly_degree(p, 1L)
  if (.sp_poly_degree(p, 2L) == d) {
    a <- root(coeff_mon(p, d, 0L), d)
    b <- root(coeff_mon(p, 0L, d), d)
    nw <- okf(f, x, .div(.sub(u, .mul(list(b, y))), a))
    if (!is.null(nw)) return(list(.add(list(.mul(list(a, x)), .mul(list(b, y)))), nw, u))
  }
  d <- .sp_poly_degree(p, 1L)
  if (.sp_poly_degree(p, 2L) == d) {
    for (itry in 1:2) {
      a <- root(coeff_mon(p, d, d), d)
      b <- root(coeff_mon(p, 0L, d), d)
      nw <- okf(f, x, .div(.div(.sub(u, .mul(list(b, y))), a), y))
      if (!is.null(nw)) return(list(.add(list(.mul(list(a, x, y)), .mul(list(b, y)))), nw, u))
      tmp <- x; x <- y; y <- tmp
      p <- .sp_as_poly(f, list(x, y))
    }
  }
  NULL
}

.sp_tsolve_bivariate <- function(eq, sym, lhs, rhs, flags) {
  logs <- .sp_find(eq, function(z) .sp_is_func(z, "log"))
  spows <- list()
  for (j in logs) for (i in atoms(j, "Pow")) if (identical(i$args[[1L]]$key, sym$key)) spows[[length(spows) + 1L]] <- i$args[[2L]]
  spow <- S_One
  if (length(spows)) {
    vals <- vapply(spows, function(z) tryCatch(as.numeric(.neval(z, list())), error = function(e) NA_real_), 0)
    if (!anyNA(vals)) spow <- spows[[which.min(vals)]]
  }
  if (!is_one_num(spow)) {
    p <- .pow(sym, spow)
    u <- Dummy("bivariate-cov")
    ueq <- subs(eq, p, u)
    if (!.has_free(ueq, list(sym))) {
      sol <- .sp_vsolve(ueq, u, flags)
      inv <- .sp_vsolve(.sub(p, u), sym, list())
      out <- list()
      for (i in inv) for (s in sol) out[[length(out) + 1L]] <- subs(i, u, s)
      return(out)
    }
  }
  pg <- .sp_poly_gens(eq)
  g <- .sp_filtered_gens(if (is.null(pg)) list() else pg, sym)
  up_or_log <- list()
  for (gi in g) {
    if (.sp_is_func(gi, c("exp", "log")) || (gi$type == "Pow" && identical(gi$args[[1L]]$key, S_E$key))) {
      up_or_log[[length(up_or_log) + 1L]] <- gi
    } else if (gi$type == "Pow") {
      gisimp <- powdenest(expand_power_exp(gi))
      if (gisimp$type == "Pow" && .sp_in(sym, free_symbols(gisimp$args[[2L]]))) up_or_log[[length(up_or_log) + 1L]] <- gi
    }
  }
  eq_down <- expand_log(expand_power_exp(eq))
  if (length(up_or_log)) eq_down <- subs(eq_down, lapply(up_or_log, function(z) list(z, S_Zero)))
  eq <- expand_power_exp(.add(list(.sps_factor_deep(eq_down), .sub(eq, eq_down))))
  iv <- .sp_invert(eq, list(sym))
  rhs <- iv[[1L]]; lhs <- iv[[2L]]
  if (.has1(lhs, sym)) {
    pg <- .sp_poly_gens(lhs)
    g <- .sp_filtered_gens(if (is.null(pg)) list() else pg, sym)
    eq_ <- .sub(lhs, rhs)
    sols <- tryCatch(.sp_solve_lambert(eq_, sym, g), NotImplementedError = function(e) NULL)
    if (!is.null(sols)) {
      for (n in seq_along(sols)) {
        s <- sols[[n]]
        ns <- .sps_nsimplify(s)
        if (!identical(ns$key, s$key) && count_ops(ns) <= count_ops(s)) {
          ok <- .sp_checksol(eq_, .sd1(sym, ns))
          if (is.na(ok)) ok <- isTRUE(.expr_equals(subs(eq_, sym, ns), S_Zero))
          if (isTRUE(ok)) sols[[n]] <- ns
        }
      }
      return(sols)
    }
    if (length(g) == 2L) {
      r <- tryCatch({
        gpu <- .sp_bivariate_type(.sub(lhs, rhs), g[[1L]], g[[2L]])
        if (is.null(gpu)) .sp_stop("NotImplementedError", "")
        fl <- flags
        fl$bivariate <- FALSE
        inversion <- .sp_tsolve(.sub(gpu[[1L]], gpu[[3L]]), sym, fl)
        if (length(inversion)) {
          sol <- .sp_vsolve(gpu[[2L]], gpu[[3L]], fl)
          out <- list()
          for (i in inversion) for (s in sol) out[[length(out) + 1L]] <- subs(i, gpu[[3L]], s)
          .sp_uniq(out)
        } else "fallthrough"
      }, NotImplementedError = function(e) "fallthrough", PolynomialError = function(e) "fallthrough")
      if (!is.character(r)) return(r)
    }
  }
  structure("fallthrough", lhs = lhs, rhs = rhs)
}

#' Remove radicals from an equation
#'
#' `unrad(eq, ...)` returns an equation free of the radicals in the symbols
#' `...` (all free symbols when none are given) whose roots include those of
#' `eq`, as SymPy's `unrad()`. The result is the pair `(eq, cov)`, where
#' `cov` is empty or holds a change of variables `p` and the equation
#' relating `p` to the original symbol; it is `NULL` when there are no
#' radicals to remove. An error of class `NotImplementedError` is raised
#' when the radicals cannot all be removed.
#'
#' @param eq An expression equal to zero, or an [Eq()].
#' @param ... Symbols of interest.
#' @return A two element `pytuple` or `NULL`.
#' @examples
#' x <- Symbol("x")
#' unrad(sqrt(x) + sqrt(1 - x))
#' unrad(x^Rational(5, 6) + 2)
#' @export
unrad <- function(eq, ...) {
  syms <- lapply(list(...), as_sym)
  r <- .sp_unrad(as_sym(eq), syms)
  if (is.null(r)) return(NULL)
  .pytuple(list(r[[1L]], r[[2L]]))
}
