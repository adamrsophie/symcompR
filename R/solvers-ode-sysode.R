# The systems part of sympy/solvers/ode/ode.py: dsolve() for a list of
# equations, classify_sysode(), and the solvers of the special two and three
# equation systems that dsolve_system() does not handle.

.ode_terms <- function(e) if (is_sym(e) && e$type == "Add") .sympy_args(e) else list(e)
.ode_factors <- function(e) if (is_sym(e) && e$type == "Mul") .sympy_args(e) else list(e)

# the independent variable of a system: a Symbol of a Derivative of the
# first equation
.ode_sys_t <- function(eq1) {
  d <- atoms(eq1, "Derivative")
  if (!length(d)) .ode_stop("IndexError", "list index out of range")
  Filter(function(s) s$type == "Symbol", atoms(d[[1L]], "Symbol"))[[1L]]
}

.ode_fc_key <- function(j, func, k) paste(j, func$key, k, sep = "|")

# dsolve() for an iterable of equations
.ode_dsolve_iterable <- function(eq, func, ics) {
  sol <- tryCatch(.ode_dsolve_system(eq, funcs = func, ics = ics, doit = TRUE), NotImplementedError = function(e) e)
  if (!inherits(sol, "condition")) return(if (length(sol) == 1L) sol[[1L]] else sol)
  match <- .ode_cls_sysode(eq, func)
  eq <- match$eq
  order <- match$order
  func <- match$func
  t <- .ode_sys_t(eq[[1L]])
  for (i in seq_along(eq)) for (func_ in func) {
    if (is.list(func_) && !is_sym(func_)) next
    c0 <- coeff(eq[[i]], .ode_d(func[[i]], t, .ode_order(eq[[i]], func[[i]])))
    if (isTRUE(.sp_fact(c0, "negative"))) eq[[i]] <- .neg(eq[[i]])
  }
  match$eq <- eq
  if (length(unique(order)) != 1L) .ode_stop("ValueError", "It solves only those systems of equations whose orders are equal")
  match$order <- order[[1L]]
  recur_len <- function(l) sum(vapply(l, function(it) if (is.list(it) && !is_sym(it)) recur_len(it) else 1L, 0L))
  if (recur_len(func) != length(eq))
    .ode_stop("ValueError", "dsolve() and classify_sysode() work with number of functions being equal to number of equations")
  if (is.null(match$type_of_equation)) .ode_nie("")
  nm <- sprintf("sysode_%s_%deq_order%d", if (isTRUE(match$is_linear)) "linear" else "nonlinear", match$no_of_equation, match$order)
  solvefunc <- switch(nm, sysode_linear_2eq_order1 = .ode_sysode_linear_2eq_order1,
                      sysode_nonlinear_2eq_order1 = .ode_sysode_nonlinear_2eq_order1,
                      sysode_nonlinear_3eq_order1 = .ode_sysode_nonlinear_3eq_order1,
                      .ode_stop("KeyError", paste0("'", nm, "'")))
  sols <- solvefunc(match)
  if (length(ics)) {
    fs <- list()
    for (s in sols) fs <- c(fs, free_symbols(s))
    ek <- character()
    for (e in eq) ek <- c(ek, vapply(free_symbols(e), .k, ""))
    constants <- Filter(function(s) !(s$key %in% ek), fs)
    constants <- constants[!duplicated(vapply(constants, .k, ""))]
    sc <- .ode_solve_ics(sols, func, constants, ics)
    return(lapply(sols, function(s) .ode_subs_dict(s, sc)))
  }
  sols
}

# classify_sysode as a named list; func_coeff is keyed by .ode_fc_key and
# its insertion order kept in fc_keys
.ode_cls_sysode <- function(eq, funcs = NULL) {
  eq <- if (is.list(eq) && !is_sym(eq)) lapply(eq, as_sym) else list(as_sym(eq))
  eq <- lapply(eq, function(fi) if (fi$type == "Equality") .sub(fi$args[[1L]], fi$args[[2L]]) else fi)
  t <- .ode_sys_t(eq[[1L]])
  n <- length(eq)
  if (n == 1L) .ode_stop("ValueError", "classify_sysode() works for systems of ODEs. For scalar ODEs, classify_ode should be used")
  funcs <- if (is.null(funcs)) .ode_extract_funcs(eq) else if (is.list(funcs) && !is_sym(funcs)) lapply(funcs, as_sym) else list(as_sym(funcs))
  funcs <- funcs[!duplicated(vapply(funcs, .k, ""))]
  if (length(funcs) != n)
    .ode_stop("ValueError", paste0("Number of functions given is not equal to the number of equations ", .ode_list_str(funcs)))
  order <- list()
  func_dict <- list()
  eq_no <- NULL
  for (func in funcs) {
    if (is.null(order[[func$key]]) || order[[func$key]] == 0L) {
      max_order <- 0L
      for (i in seq_along(eq)) {
        o <- .ode_order(eq[[i]], func)
        if (max_order < o) {
          max_order <- o
          eq_no <- i - 1L
        }
      }
      if (is.null(eq_no)) .ode_stop("UnboundLocalError", "cannot access local variable 'eq_no' where it is not associated with a value")
      k <- as.character(eq_no)
      func_dict[[k]] <- if (!is.null(func_dict[[k]])) list(func_dict[[k]], func) else func
      order[[func$key]] <- max_order
    }
  }
  funcs <- lapply(seq_along(func_dict) - 1L, function(i) {
    v <- func_dict[[as.character(i)]]
    if (is.null(v)) .ode_stop("KeyError", as.character(i))
    v
  })
  flat <- unlist(lapply(funcs, function(f) if (is.list(f) && !is_sym(f)) f else list(f)), recursive = FALSE)
  for (f in flat) if (length(f$args) != 1L)
    .ode_stop("ValueError", paste0("dsolve() and classify_sysode() work with functions of one variable only, not ", .str(f)))
  ofuncs <- lapply(names(order), function(k) flat[[match(k, vapply(flat, .k, ""))]])
  fck <- character()
  fcv <- list()
  fct <- list()
  linearity_check <- function(eqs, j, func, is_linear_) {
    for (k in 0:order[[func$key]]) {
      dk <- .ode_d(func, t, k)
      cf <- coeff(collect(expand(eqs), list(dk)), dk)
      key <- .ode_fc_key(j, func, k)
      if (!(key %in% fck)) {
        fck <<- c(fck, key)
        fct[[length(fct) + 1L]] <<- .pytuple(list(j, func, k))
      }
      fcv[[key]] <<- cf
      if (isTRUE(is_linear_)) {
        if (is_zero_num(cf)) {
          if (k == 0L) {
            coef <- .sp_indep(eqs, list(func), as_add = TRUE)[[2L]]
            ox <- .ode_order(eqs, func)
            if (ox >= 1L) for (xr in seq_len(ox))
              coef <- .sub(coef, .sp_indep(eqs, list(.ode_d(func, t, xr)), as_add = TRUE)[[2L]])
            if (!is_zero_num(coef)) is_linear_ <- FALSE
          } else if (!is_zero_num(.sp_indep(eqs, list(dk), as_add = TRUE)[[2L]])) is_linear_ <- FALSE
        } else {
          for (f2 in flat) {
            dep <- .sp_indep(cf, list(f2), as_add = TRUE)[[2L]]
            if (!is_zero_num(dep)) is_linear_ <- FALSE
          }
        }
      }
    }
    is_linear_
  }
  is_linear <- TRUE
  for (j in seq_along(eq)) for (func in funcs) {
    if (is.list(func) && !is_sym(func)) for (fe in func) is_linear <- linearity_check(eq[[j]], j - 1L, fe, is_linear)
    else is_linear <- linearity_check(eq[[j]], j - 1L, func, is_linear)
  }
  fc <- list(keys = fck, values = fcv[fck], tuples = fct)
  type <- NULL
  ov <- unlist(order)
  if (length(unique(ov)) == 1L) {
    order_eq <- ov[[1L]]
    if (is_linear) {
      if (n == 2L && order_eq == 1L) type <- .ode_check_linear_2eq_order1(eq, funcs, fc)
    } else if (n == 2L && order_eq == 1L) {
      r <- .ode_check_nonlinear_2eq_order1(eq, funcs, fc)
      type <- r$type
      eq <- r$eq
    } else if (n == 3L && order_eq == 1L) {
      r <- .ode_check_nonlinear_3eq_order1(eq, funcs, fc)
      type <- r$type
      eq <- r$eq
    }
  }
  m <- list(no_of_equation = n, eq = eq, func = funcs, order = ov, ofuncs = ofuncs, func_coeff = fc, is_linear = is_linear)
  m["type_of_equation"] <- list(type)
  m
}

#' Classify a system of ODEs
#'
#' `classify_sysode(eq, funcs)` returns the properties of the system of ODEs
#' `eq` that dsolve() uses to pick a solver, as SymPy's `classify_sysode()`:
#' the number of equations, the equations as expressions, the functions,
#' their orders, the coefficients of the functions and their derivatives,
#' the linearity, and the type of the system.
#'
#' @param eq A list of ODEs.
#' @param funcs A list of the unknown functions, or `NULL` to find them.
#' @return A dictionary keyed by the names of the properties.
#' @examples
#' t <- Symbol("t"); x <- Function("x"); y <- Function("y")
#' classify_sysode(list(Eq(Derivative(x(t), t), y(t)^3 * x(t)), Eq(Derivative(y(t), t), y(t)^5)))
#' @export
classify_sysode <- function(eq, funcs = NULL) {
  m <- .ode_cls_sysode(eq, funcs)
  fc <- .pydict(m$func_coeff$tuples, unname(m$func_coeff$values))
  od <- .pydict(m$ofuncs, lapply(unname(m$order), function(o) o))
  .pydict(list("no_of_equation", "eq", "func", "order", "func_coeff", "is_linear", "type_of_equation"),
          list(m$no_of_equation, m$eq, m$func, od, fc, m$is_linear, m$type_of_equation))
}

.ode_fc <- function(fc, j, func, k) {
  v <- fc$values[[.ode_fc_key(j, func, k)]]
  if (is.null(v)) .ode_stop("KeyError", paste0("(", j, ", ", .str(func), ", ", k, ")"))
  v
}

# ode.check_linear_2eq_order1
.ode_check_linear_2eq_order1 <- function(eq, func, fc) {
  x <- func[[1L]]
  y <- func[[2L]]
  t <- .ode_sys_t(eq[[1L]])
  r <- list()
  r$a1 <- .ode_fc(fc, 0L, x, 1L); r$a2 <- .ode_fc(fc, 1L, y, 1L)
  r$b1 <- .neg(.div(.ode_fc(fc, 0L, x, 0L), .ode_fc(fc, 0L, x, 1L)))
  r$b2 <- .neg(.div(.ode_fc(fc, 1L, x, 0L), .ode_fc(fc, 1L, y, 1L)))
  r$c1 <- .neg(.div(.ode_fc(fc, 0L, y, 0L), .ode_fc(fc, 0L, x, 1L)))
  r$c2 <- .neg(.div(.ode_fc(fc, 1L, y, 0L), .ode_fc(fc, 1L, y, 1L)))
  forcing <- list(S_Zero, S_Zero)
  for (i in 1:2) for (j in .ode_terms(eq[[i]])) if (!.has1(j, x) && !.has1(j, y)) forcing[[i]] <- .add(list(forcing[[i]], j))
  if (.has1(forcing[[1L]], t) || .has1(forcing[[2L]], t)) return(NULL)
  r$d1 <- forcing[[1L]]
  r$d2 <- forcing[[2L]]
  p <- 0L
  q <- S_Zero
  p1 <- cancel(.div(r$b2, .as_numer_denom(cancel(.div(r$b2, r$c2)))[[1L]]))
  p2 <- cancel(.div(r$b1, .as_numer_denom(cancel(.div(r$b1, r$c1)))[[1L]]))
  ps <- list(p1, p2)
  for (n in 0:1) for (j in .ode_factors(collect_const(ps[[n + 1L]]))) {
    if (!.has1(j, t)) q <- j
    if (!is_zero_num(q) && n == 0L) {
      if (identical(.div(.sub(.div(r$b2, j), r$b1), .sub(r$c1, .div(r$c2, j)))$key, j$key)) p <- 1L
    } else if (!is_zero_num(q) && n == 1L) {
      if (identical(.div(.sub(.div(r$b1, j), r$b2), .sub(r$c2, .div(r$c1, j)))$key, j$key)) p <- 2L
    }
  }
  if (!is_zero_num(r$d1) || !is_zero_num(r$d2)) return(NULL)
  if (!any(vapply(r[c("a1", "a2", "b1", "b2", "c1", "c2")], function(z) .has1(z, t), TRUE))) return(NULL)
  if (p) "type6" else "type7"
}

# ode.check_nonlinear_2eq_order1: the type and the equations, which it
# divides by the coefficients of the derivatives as SymPy does in place
.ode_check_nonlinear_2eq_order1 <- function(eq, func, fc) {
  t <- .ode_sys_t(eq[[1L]])
  f <- Wild("f")
  g <- Wild("g")
  u <- Dummy("u")
  v <- Dummy("v")
  for (func_ in func) {
    if (is.list(func_) && !is_sym(func_)) {
      x <- func[[1L]][[1L]]
      y <- func[[1L]][[2L]]
      tp <- .ode_type5_check(eq, x, y, t, f, g)
      ok <- function(tp) !is.null(tp[[1L]]) && !is.null(tp[[2L]]) && {
        sv <- function(e) subs(subs(e, .ode_d(x, t), u), .ode_d(y, t), v)
        !(.has1(sv(tp[[1L]]$f), t) || .has1(sv(tp[[2L]]$g), t))
      }
      type <- if (ok(tp)) "type5" else NULL
      if (is.null(type)) {
        tp <- .ode_type5_check(eq, y, x, t, f, g)
        type <- if (ok(tp)) "type5" else NULL
      }
      return(list(type = type, eq = eq))
    }
  }
  x <- func[[1L]]
  y <- func[[2L]]
  n <- Wild("n", exclude = list(x, y))
  for (i in 1:2) eq[[i]] <- .add(lapply(.ode_terms(eq[[i]]), function(z) .div(z, .ode_fc(fc, i - 1L, func[[i]], 1L))))
  dx <- .ode_d(x, t)
  dy <- .ode_d(y, t)
  sub_uv <- function(e) subs(subs(e, x, u), y, v)
  r <- sym_match(eq[[1L]], .sub(dx, .mul(list(.pow(x, n), f))))
  if (!is.null(r)) {
    gg <- .div(.sub(dy, eq[[2L]]), r$f)
    if (!(.has1(gg, x) || .has1(subs(gg, y, v), t) || .has1(sub_uv(r$f), t))) return(list(type = "type1", eq = eq))
  }
  r <- sym_match(eq[[1L]], .sub(dx, .mul(list(fcall("exp", list(.mul(list(n, x)))), f))))
  if (!is.null(r)) {
    gg <- .div(.sub(dy, eq[[2L]]), r$f)
    if (!(.has1(gg, x) || .has1(subs(gg, y, v), t) || .has1(sub_uv(r$f), t))) return(list(type = "type2", eq = eq))
  }
  r1 <- sym_match(eq[[1L]], .sub(dx, f))
  r2 <- sym_match(eq[[2L]], .sub(dy, g))
  if (!is.null(r1) && !is.null(r2) && !(.has1(sub_uv(r1$f), t) || .has1(sub_uv(r2$g), t))) return(list(type = "type3", eq = eq))
  if (is.null(r1) || is.null(r2)) .ode_stop("TypeError", "'NoneType' object is not subscriptable")
  nd <- .as_numer_denom(.div(sub_uv(r1$f), sub_uv(r2$g)))
  f1 <- Wild("f1", exclude = list(v, t)); f2 <- Wild("f2", exclude = list(v, t))
  g1 <- Wild("g1", exclude = list(u, t)); g2 <- Wild("g2", exclude = list(u, t))
  R1 <- sym_match(nd[[1L]], .mul(list(f1, g1)))
  R2 <- sym_match(nd[[2L]], .mul(list(f2, g2)))
  list(type = if (length(R1) && length(R2)) "type4" else NULL, eq = eq)
}

# the matches of the type5 system for the functions x and y
.ode_type5_check <- function(eq, x, y, t, f, g) {
  dx <- .ode_d(x, t)
  dy <- .ode_d(y, t)
  p1 <- function(fn, d, w) .add(list(.mul(list(t, d)), .neg(fn), w))
  p2 <- function(fn, d, w) .add(list(d, .neg(.div(fn, t)), .div(w, t)))
  both <- function(r) length(r[[1L]]) && length(r[[2L]])
  r <- list(sym_match(eq[[1L]], p1(x, dx, f)), sym_match(eq[[2L]], p1(y, dy, g)))
  if (!both(r)) r <- list(sym_match(eq[[1L]], p2(x, dx, f)), sym_match(eq[[2L]], p2(y, dy, g)))
  if (!both(r)) r <- list(sym_match(.neg(eq[[1L]]), p1(x, dx, f)), sym_match(.neg(eq[[2L]]), p1(y, dy, g)))
  if (!both(r)) r <- list(sym_match(.neg(eq[[1L]]), p2(x, dx, f)), sym_match(.neg(eq[[2L]]), p2(y, dy, g)))
  r
}

# ode.check_nonlinear_3eq_order1
.ode_check_nonlinear_3eq_order1 <- function(eq, func, fc) {
  x <- func[[1L]]; y <- func[[2L]]; z <- func[[3L]]
  t <- .ode_sys_t(eq[[1L]])
  u <- Dummy("u"); v <- Dummy("v"); w <- Dummy("w")
  ex <- list(x, y, z, t)
  a <- Wild("a", exclude = ex); b <- Wild("b", exclude = ex); c <- Wild("c", exclude = ex)
  f <- Wild("f")
  F1 <- Wild("F1"); F2 <- Wild("F2"); F3 <- Wild("F3")
  for (i in 1:3) eq[[i]] <- .add(lapply(.ode_terms(eq[[i]]), function(e) .div(e, .ode_fc(fc, i - 1L, func[[i]], 1L))))
  dx <- .ode_d(x, t); dy <- .ode_d(y, t); dz <- .ode_d(z, t)
  m <- function(e, p) sym_match(e, p)
  M <- function(...) .mul(list(...))
  out <- function(type) list(type = type, eq = eq)
  lin_ok <- function(ra, rb, rc, s3) {
    n1 <- .as_numer_denom(ra); n2 <- .as_numer_denom(rb); n3 <- .as_numer_denom(rc)
    eqs <- list(.sub(M(n1[[1L]], u), M(n1[[2L]], .sub(v, w))), .sub(M(n2[[1L]], v), M(n2[[2L]], .sub(w, u))),
                .sub(M(n3[[1L]], w), M(n3[[2L]], .sub(u, v))))
    length(.sp_solve(eqs, list(u, v), list())) > 0L
  }
  r1 <- m(eq[[1L]], .sub(dx, M(a, y, z)))
  r2 <- m(eq[[2L]], .sub(dy, M(b, z, x)))
  r3 <- m(eq[[3L]], .sub(dz, M(c, x, y)))
  if (length(r1) && length(r2) && length(r3) && lin_ok(r1$a, r2$b, r3$c)) return(out("type1"))
  r <- m(eq[[1L]], .sub(dx, M(y, z, f)))
  if (length(r)) {
    r1 <- m(collect_const(r$f), M(a, f))
    r2 <- m(.div(.sub(dy, eq[[2L]]), r1$f), M(b, z, x))
    r3 <- m(.div(.sub(dz, eq[[3L]]), r1$f), M(c, x, y))
  }
  if (length(r1) && length(r2) && length(r3) && lin_ok(r1$a, r2$b, r3$c)) return(out("type2"))
  swap <- function(r1) {
    if (.has1(eq[[2L]], r1$F2) && !.has1(eq[[2L]], r1$F3)) {
      tmp <- r1$F2; r1$F2 <- r1$F3; r1$F3 <- tmp
      tmp <- r1$c; r1$c <- .neg(r1$b); r1$b <- .neg(tmp)
    }
    r1
  }
  upd <- function(a1, a2) {
    if (is.null(a1)) .ode_stop("AttributeError", "'NoneType' object has no attribute 'update'")
    for (nm in names(a2)) a1[[nm]] <- a2[[nm]]
    a1
  }
  r <- m(eq[[1L]], .sub(dx, .sub(F2, F3)))
  if (length(r)) {
    r1 <- upd(m(collect_const(r$F2), M(c, F2)), m(collect_const(r$F3), M(b, F3)))
    r2 <- NULL
    r3 <- FALSE
    if (length(r1)) {
      r1 <- swap(r1)
      r2 <- m(eq[[2L]], .add(list(dy, .neg(M(a, r1$F3)), M(r1$c, F1))))
    }
    if (length(r2)) r3 <- identical(eq[[3L]]$key, .add(list(dz, .neg(M(r1$b, r2$F1)), M(r2$a, r1$F2)))$key)
    if (length(r1) && length(r2) && r3) return(out("type3"))
  }
  r <- m(eq[[1L]], .add(list(dx, .neg(M(z, F2)), M(y, F3))))
  if (length(r)) {
    r1 <- upd(m(collect_const(r$F2), M(c, F2)), m(collect_const(r$F3), M(b, F3)))
    r2 <- NULL
    r3 <- FALSE
    if (length(r1)) {
      r1 <- swap(r1)
      r2 <- m(.sub(dy, eq[[2L]]), .sub(M(a, x, r1$F3), M(r1$c, z, F1)))
    }
    if (length(r2)) r3 <- identical(.sub(dz, eq[[3L]])$key, .sub(M(r1$b, y, r2$F1), M(r2$a, x, r1$F2))$key)
    if (length(r1) && length(r2) && r3) return(out("type4"))
  }
  r <- m(.sub(dx, eq[[1L]]), M(x, .sub(F2, F3)))
  if (length(r)) {
    r1 <- upd(m(collect_const(r$F2), M(c, F2)), m(collect_const(r$F3), M(b, F3)))
    r2 <- NULL
    r3 <- FALSE
    if (length(r1)) {
      r1 <- swap(r1)
      r2 <- m(.sub(dy, eq[[2L]]), M(y, .sub(M(a, r1$F3), M(r1$c, F1))))
    }
    if (length(r2)) r3 <- identical(.sub(dz, eq[[3L]])$key, M(z, .sub(M(r1$b, r2$F1), M(r2$a, r1$F2)))$key)
    if (length(r1) && length(r2) && r3) return(out("type5"))
  }
  out(NULL)
}

# ode.sysode_linear_2eq_order1
.ode_sysode_linear_2eq_order1 <- function(match_) {
  func <- match_$func
  x <- func[[1L]]; y <- func[[2L]]
  fc <- match_$func_coeff
  eq <- match_$eq
  t <- .ode_sys_t(eq[[1L]])
  for (i in 1:2) eq[[i]] <- .add(lapply(.ode_terms(eq[[i]]), function(z) .div(z, .ode_fc(fc, i - 1L, func[[i]], 1L))))
  r <- list()
  r$a <- .neg(.div(.ode_fc(fc, 0L, x, 0L), .ode_fc(fc, 0L, x, 1L)))
  r$c <- .neg(.div(.ode_fc(fc, 1L, x, 0L), .ode_fc(fc, 1L, y, 1L)))
  r$b <- .neg(.div(.ode_fc(fc, 0L, y, 0L), .ode_fc(fc, 0L, x, 1L)))
  r$d <- .neg(.div(.ode_fc(fc, 1L, y, 0L), .ode_fc(fc, 1L, y, 1L)))
  forcing <- list(S_Zero, S_Zero)
  for (i in 1:2) for (j in .ode_terms(eq[[i]])) if (!.has1(j, x) && !.has1(j, y)) forcing[[i]] <- .add(list(forcing[[i]], j))
  if (.has1(forcing[[1L]], t) || .has1(forcing[[2L]], t))
    .ode_nie("Only homogeneous problems are supported (and constant inhomogeneity)")
  r$k1 <- forcing[[1L]]
  r$k2 <- forcing[[2L]]
  if (match_$type_of_equation == "type6") return(.ode_linear_2eq_order1_type6(x, y, t, r, eq))
  .ode_linear_2eq_order1_type7(x, y, t, r, eq)
}

.ode_linear_2eq_order1_type6 <- function(x, y, t, r, eq) {
  C <- .ode_get_numbered_constants(eq, num = 4L)
  C1 <- C[[1L]]
  p <- 0L
  q <- S_Zero
  s <- NULL
  p1 <- cancel(.div(r$c, .as_numer_denom(cancel(.div(r$c, r$d)))[[1L]]))
  p2 <- cancel(.div(r$a, .as_numer_denom(cancel(.div(r$a, r$b)))[[1L]]))
  ps <- list(p1, p2)
  for (n in 0:1) for (j in .ode_factors(collect_const(ps[[n + 1L]]))) {
    if (!.has1(j, t)) q <- j
    if (!is_zero_num(q) && n == 0L &&
        identical(.div(.sub(.div(r$c, j), r$a), .sub(r$b, .div(r$d, j)))$key, j$key)) {
      p <- 1L; s <- j
      break
    }
    if (!is_zero_num(q) && n == 1L &&
        identical(.div(.sub(.div(r$a, j), r$c), .sub(r$d, .div(r$b, j)))$key, j$key)) {
      p <- 2L; s <- j
      break
    }
  }
  if (p == 1L) {
    ex <- .mul(list(C1, fcall("exp", list(.mul(list(.neg(s), Integral(.sub(r$b, .div(r$d, s)), t)))))))
    equ <- .add(list(.ode_d(x, t), .neg(.mul(list(r$a, x))), .neg(.mul(list(r$b, .add(list(.mul(list(s, x)), ex)))))))
    hint1 <- .ode_classify(equ, NULL, kw = list())[[2L]]
    sol1 <- dsolve(equ, hint = paste0(hint1, "_Integral"))$args[[2L]]
    sol2 <- .add(list(.mul(list(s, sol1)), ex))
  } else if (p == 2L) {
    ex <- .mul(list(C1, fcall("exp", list(.mul(list(.neg(s), Integral(.sub(r$d, .div(r$b, s)), t)))))))
    equ <- .add(list(.ode_d(y, t), .neg(.mul(list(r$c, y))), .neg(.mul(list(r$d, s, y))), ex))
    hint1 <- .ode_classify(equ, NULL, kw = list())[[2L]]
    sol2 <- dsolve(equ, hint = paste0(hint1, "_Integral"))$args[[2L]]
    sol1 <- .add(list(.mul(list(s, sol2)), ex))
  } else .ode_stop("UnboundLocalError", "cannot access local variable 'sol1' where it is not associated with a value")
  list(Eq(x, sol1), Eq(y, sol2))
}

.ode_linear_2eq_order1_type7 <- function(x, y, t, r, eq) {
  C <- .ode_get_numbered_constants(eq, num = 4L)
  C1 <- C[[1L]]; C2 <- C[[2L]]
  a <- r$a; b <- r$b; c <- r$c; d <- r$d
  da <- diff(a, t); db <- diff(b, t); dc <- diff(c, t); dd <- diff(d, t)
  M <- function(...) .mul(list(...))
  e1 <- .add(list(M(a, b, c), .neg(M(.pow(b, S_Two), c)), M(a, db), .neg(M(da, b))))
  e2 <- .add(list(M(a, c, d), .neg(M(b, .pow(c, S_Two))), M(dc, d), .neg(M(c, dd))))
  m1 <- .add(list(M(a, b), M(b, d), db))
  m2 <- .add(list(M(a, c), M(c, d), dc))
  rhs <- function(e) dsolve(e)$args[[2L]]
  dx <- .ode_d(x, t); dy <- .ode_d(y, t)
  if (is_zero_num(e1)) {
    sol1 <- rhs(.sub(M(b, .ode_d(x, t, 2L)), M(m1, dx)))
    sol2 <- rhs(.add(list(dy, .neg(M(c, sol1)), .neg(M(d, y)))))
  } else if (is_zero_num(e2)) {
    sol2 <- rhs(.sub(M(c, .ode_d(y, t, 2L)), M(m2, dy)))
    sol1 <- rhs(.add(list(dx, .neg(M(a, x)), .neg(M(b, sol2)))))
  } else if (!.has1(.div(e1, b), t) && !.has1(.div(m1, b), t)) {
    sol1 <- rhs(.add(list(.ode_d(x, t, 2L), .neg(M(.div(m1, b), dx)), .neg(M(.div(e1, b), x)))))
    sol2 <- rhs(.add(list(dy, .neg(M(c, sol1)), .neg(M(d, y)))))
  } else if (!.has1(.div(e2, c), t) && !.has1(.div(m2, c), t)) {
    sol2 <- rhs(.add(list(.ode_d(y, t, 2L), .neg(M(.div(m2, c), dy)), .neg(M(.div(e2, c), y)))))
    sol1 <- rhs(.add(list(dx, .neg(M(a, x)), .neg(M(b, sol2)))))
  } else {
    x0 <- fcall("x0", list(t))
    y0 <- fcall("y0", list(t))
    F <- fcall("exp", list(Integral(a, t)))
    P <- fcall("exp", list(Integral(d, t)))
    I1 <- Integral(.div(M(b, F, P), .pow(x0, S_Two)), t)
    sol1 <- .add(list(M(C1, x0), M(C2, x0, I1)))
    sol2 <- .add(list(M(C1, y0), M(C2, .add(list(.div(M(F, P), x0), M(y0, I1))))))
  }
  list(Eq(x, sol1), Eq(y, sol2))
}

# ode.sysode_nonlinear_2eq_order1
.ode_sysode_nonlinear_2eq_order1 <- function(match_) {
  func <- match_$func
  eq <- match_$eq
  fc <- match_$func_coeff
  t <- .ode_sys_t(eq[[1L]])
  if (match_$type_of_equation == "type5") return(.ode_nonlinear_2eq_order1_type5(func, t, eq))
  x <- func[[1L]]; y <- func[[2L]]
  for (i in 1:2) eq[[i]] <- .add(lapply(.ode_terms(eq[[i]]), function(z) .div(z, .ode_fc(fc, i - 1L, func[[i]], 1L))))
  switch(match_$type_of_equation,
         type1 = .ode_nonlinear_2eq_order1_type12(x, y, t, eq, FALSE),
         type2 = .ode_nonlinear_2eq_order1_type12(x, y, t, eq, TRUE),
         type3 = .ode_nonlinear_2eq_order1_type3(x, y, t, eq),
         type4 = .ode_nonlinear_2eq_order1_type4(x, y, t, eq))
}

# types 1 (x' = x^n F(x, y)) and 2 (x' = exp(n x) F(x, y))
.ode_nonlinear_2eq_order1_type12 <- function(x, y, t, eq, type2) {
  C <- .ode_get_numbered_constants(eq, num = 2L)
  C1 <- C[[1L]]; C2 <- C[[2L]]
  n <- Wild("n", exclude = list(x, y))
  f <- Wild("f")
  u <- Symbol("u"); v <- Symbol("v")
  pat <- if (type2) fcall("exp", list(.mul(list(n, x)))) else .pow(x, n)
  r <- sym_match(eq[[1L]], .sub(.ode_d(x, t), .mul(list(pat, f))))
  g <- subs(.div(.sub(.ode_d(y, t), eq[[2L]]), r$f), y, v)
  F <- subs(subs(r$f, x, u), y, v)
  n <- r$n
  ig <- Integral(.inv(g), v)
  if (type2) {
    phi <- if (!is_zero_num(n)) .mul(list(.neg(.inv(n)), fcall("log", list(.sub(C1, .mul(list(n, ig))))))) else .add(list(C1, ig))
  } else {
    phi <- if (!identical(n$key, S_One$key)) .pow(.add(list(C1, .mul(list(.sub(S_One, n), ig)))), .inv(.sub(S_One, n)))
           else .mul(list(C1, fcall("exp", list(ig))))
  }
  phi <- doit(phi)
  ex <- .sub(.sub(doit(Integral(.inv(.mul(list(g, subs(F, u, phi)))), v)), t), C2)
  sol2 <- .sp_solve(ex, list(v), list())
  sol <- list()
  for (s in sol2) {
    sol[[length(sol) + 1L]] <- Eq(x, subs(phi, v, s))
    sol[[length(sol) + 1L]] <- Eq(y, s)
  }
  sol
}

.ode_nonlinear_2eq_order1_type3 <- function(x, y, t, eq) {
  C <- .ode_get_numbered_constants(eq, num = 4L)
  C2 <- C[[2L]]
  u <- Symbol("u")
  vu <- fcall("v", list(u))
  f <- Wild("f"); g <- Wild("g")
  r1 <- sym_match(eq[[1L]], .sub(.ode_d(x, t), f))
  r2 <- sym_match(eq[[2L]], .sub(.ode_d(y, t), g))
  F <- subs(subs(r1$f, x, u), y, vu)
  G <- subs(subs(r2$g, x, u), y, vu)
  sol2r <- dsolve(Eq(.ode_d(vu, u), .div(G, F)))
  if (is_sym(sol2r)) sol2r <- list(sol2r)
  sol1 <- list()
  sol2s <- NULL
  for (sol2s in sol2r)
    sol1 <- .sp_solve(.sub(.sub(doit(Integral(.inv(subs(F, vu, sol2s$args[[2L]])), u)), t), C2), list(u), list())
  sol <- list()
  for (s in sol1) {
    sol[[length(sol) + 1L]] <- Eq(x, s)
    sol[[length(sol) + 1L]] <- Eq(y, subs(sol2s$args[[2L]], u, s))
  }
  sol
}

.ode_nonlinear_2eq_order1_type4 <- function(x, y, t, eq) {
  C <- .ode_get_numbered_constants(eq, num = 2L)
  C1 <- C[[1L]]
  u <- Symbol("u"); v <- Symbol("v")
  Ut <- fcall("U", list(t)); Vt <- fcall("V", list(t))
  f <- Wild("f"); g <- Wild("g")
  f1 <- Wild("f1", exclude = list(v, t)); f2 <- Wild("f2", exclude = list(v, t))
  g1 <- Wild("g1", exclude = list(u, t)); g2 <- Wild("g2", exclude = list(u, t))
  r1 <- sym_match(eq[[1L]], .sub(.ode_d(x, t), f))
  r2 <- sym_match(eq[[2L]], .sub(.ode_d(y, t), g))
  suv <- function(e) subs(subs(e, x, u), y, v)
  nd <- .as_numer_denom(.div(suv(r1$f), suv(r2$g)))
  num <- nd[[1L]]
  R1 <- sym_match(num, .mul(list(f1, g1)))
  R2 <- sym_match(nd[[2L]], .mul(list(f2, g2)))
  phi <- .div(suv(r1$f), num)
  F1 <- R1$f1; F2 <- R2$f2; G1 <- R1$g1; G2 <- R2$g2
  ex <- .sub(.sub(doit(Integral(.div(F2, F1), u)), doit(Integral(.div(G1, G2), v))), C1)
  sol1r <- .sp_solve(ex, list(u), list())
  sol2r <- .sp_solve(ex, list(v), list())
  rhs <- function(e) dsolve(e)$args[[2L]]
  sol <- list()
  for (s in sol1r)
    sol[[length(sol) + 1L]] <- Eq(y, rhs(.sub(.ode_d(Vt, t), .mul(list(subs(subs(F2, u, s), v, Vt), subs(G2, v, Vt), subs(subs(phi, u, s), v, Vt))))))
  for (s in sol2r)
    sol[[length(sol) + 1L]] <- Eq(x, rhs(.sub(.ode_d(Ut, t), .mul(list(subs(F1, u, Ut), subs(subs(G1, v, s), u, Ut), subs(subs(phi, v, s), u, Ut))))))
  .pyset(sol)
}

.ode_nonlinear_2eq_order1_type5 <- function(func, t, eq) {
  C <- .ode_get_numbered_constants(eq, num = 2L)
  C1 <- C[[1L]]; C2 <- C[[2L]]
  f <- Wild("f"); g <- Wild("g")
  x <- NULL; y <- NULL; r <- NULL
  for (func_ in func) {
    if (is.list(func_) && !is_sym(func_)) {
      x <- func[[1L]][[1L]]
      y <- func[[1L]][[2L]]
      r <- .ode_type5_check(eq, x, y, t, f, g)
      if (!(length(r[[1L]]) && length(r[[2L]]))) {
        r <- .ode_type5_check(eq, y, x, t, f, g)
        tmp <- x; x <- y; y <- tmp
      }
    }
  }
  x1 <- .ode_d(x, t); y1 <- .ode_d(y, t)
  sb <- function(e) subs(subs(e, x1, C1), y1, C2)
  .pyset(list(Eq(x, .add(list(.mul(list(C1, t)), sb(r[[1L]]$f)))), Eq(y, .add(list(.mul(list(C2, t)), sb(r[[2L]]$g))))))
}

# ode.sysode_nonlinear_3eq_order1
.ode_sysode_nonlinear_3eq_order1 <- function(match_) {
  func <- match_$func
  x <- func[[1L]]; y <- func[[2L]]; z <- func[[3L]]
  eq <- match_$eq
  t <- .ode_sys_t(eq[[1L]])
  switch(match_$type_of_equation,
         type1 = .ode_nonlinear_3eq_order1_type12(x, y, z, t, eq, FALSE),
         type2 = .ode_nonlinear_3eq_order1_type12(x, y, z, t, eq, TRUE),
         type3 = .ode_nonlinear_3eq_order1_type35(x, y, z, t, eq, FALSE),
         type4 = .ode_nonlinear_3eq_order1_type4(x, y, z, t, eq),
         type5 = .ode_nonlinear_3eq_order1_type35(x, y, z, t, eq, TRUE))
}

.ode_nonlinear_3eq_order1_type12 <- function(x, y, z, t, eq, type2) {
  C <- .ode_get_numbered_constants(eq, num = 2L)
  C1 <- C[[1L]]; C2 <- C[[2L]]
  u <- Symbol("u"); v <- Symbol("v"); w <- Symbol("w")
  ex <- list(x, y, z, t)
  p <- Wild("p", exclude = ex); q <- Wild("q", exclude = ex); s <- Wild("s", exclude = ex)
  f <- Wild("f")
  M <- function(...) .mul(list(...))
  dx <- .ode_d(x, t); dy <- .ode_d(y, t); dz <- .ode_d(z, t)
  if (type2) {
    r1 <- sym_match(.sub(dx, eq[[1L]]), M(y, z, f))
    r <- sym_match(collect_const(r1$f), M(p, f))
    r <- c(r, sym_match(.div(.sub(dy, eq[[2L]]), r$f), M(q, z, x)))
    r <- c(r, sym_match(.div(.sub(dz, eq[[3L]]), r$f), M(s, x, y)))
  } else {
    r <- sym_match(.sub(dx, eq[[1L]]), M(p, y, z))
    r <- c(r, sym_match(.sub(dy, eq[[2L]]), M(q, z, x)))
    r <- c(r, sym_match(.sub(dz, eq[[3L]]), M(s, x, y)))
  }
  n1 <- .as_numer_denom(r$p); n2 <- .as_numer_denom(r$q); n3 <- .as_numer_denom(r$s)
  e3 <- if (type2) .add(list(.neg(M(n3[[2L]], u)), M(n3[[2L]], v), M(n3[[1L]], w)))
        else .add(list(M(n3[[2L]], u), .neg(M(n3[[2L]], v)), .neg(M(n3[[1L]], w))))
  val <- .sp_solve(list(.add(list(M(n1[[1L]], u), .neg(M(n1[[2L]], v)), M(n1[[2L]], w))),
                        .add(list(M(n2[[2L]], u), M(n2[[1L]], v), .neg(M(n2[[2L]], w)))), e3),
                   list(u, v), list())
  vals <- list(val[["v"]], val[["u"]])
  cc <- lcm(.as_numer_denom(vals[[1L]])[[2L]], .as_numer_denom(vals[[2L]])[[2L]])
  if (type2) {
    a <- subs(vals[[1L]], w, cc); b <- subs(vals[[2L]], w, cc)
  } else {
    b <- subs(vals[[1L]], w, cc); a <- subs(vals[[2L]], w, cc)
  }
  c <- cc
  y_x <- .sqrt(.div(.sub(.sub(M(c, C1), C2), M(a, .sub(c, a), .pow(x, S_Two))), M(b, .sub(c, b))))
  z_x <- .sqrt(.div(.sub(.sub(M(b, C1), C2), M(a, .sub(b, a), .pow(x, S_Two))), M(c, .sub(b, c))))
  z_y <- .sqrt(.div(.sub(.sub(M(a, C1), C2), M(b, .sub(a, b), .pow(y, S_Two))), M(c, .sub(a, c))))
  x_y <- .sqrt(.div(.sub(.sub(M(c, C1), C2), M(b, .sub(c, b), .pow(y, S_Two))), M(a, .sub(c, a))))
  x_z <- .sqrt(.div(.sub(.sub(M(b, C1), C2), M(c, .sub(b, c), .pow(z, S_Two))), M(a, .sub(b, a))))
  y_z <- .sqrt(.div(.sub(.sub(M(a, C1), C2), M(c, .sub(a, c), .pow(z, S_Two))), M(b, .sub(a, b))))
  rf <- if (type2) r$f else S_One
  sol1 <- dsolve(.sub(M(a, dx), M(.sub(b, c), y_x, z_x, rf)))
  sol2 <- dsolve(.sub(M(b, dy), M(.sub(c, a), z_y, x_y, rf)))
  sol3 <- dsolve(.sub(M(c, dz), M(.sub(a, b), x_z, y_z, rf)))
  list(sol1, sol2, sol3)
}

# the shared matching of types 3, 4 and 5
.ode_3eq_match <- function(x, y, z, t, eq, r1, pat2, F1, F2, F3, p, q, s) {
  r <- sym_match(collect_const(r1$F2), .mul(list(s, F2)))
  r2 <- sym_match(collect_const(r1$F3), .mul(list(q, F3)))
  for (nm in names(r2)) r[[nm]] <- r2[[nm]]
  if (.has1(eq[[2L]], r$F2) && !.has1(eq[[2L]], r$F3)) {
    tmp <- r$F2; r$F2 <- r$F3; r$F3 <- tmp
    tmp <- r$s; r$s <- .neg(r$q); r$q <- .neg(tmp)
  }
  r3 <- sym_match(.sub(.ode_d(y, t), eq[[2L]]), pat2(r))
  for (nm in names(r3)) r[[nm]] <- r3[[nm]]
  r
}

.ode_nonlinear_3eq_order1_type35 <- function(x, y, z, t, eq, type5) {
  C1 <- .ode_get_numbered_constants(eq, num = 1L)
  u <- Symbol("u"); v <- Symbol("v"); w <- Symbol("w")
  ex <- list(x, y, z, t)
  p <- Wild("p", exclude = ex); q <- Wild("q", exclude = ex); s <- Wild("s", exclude = ex)
  F1 <- Wild("F1"); F2 <- Wild("F2"); F3 <- Wild("F3")
  M <- function(...) .mul(list(...))
  if (type5) {
    r1 <- sym_match(eq[[1L]], .add(list(.ode_d(x, t), .neg(M(x, F2)), M(x, F3))))
    pat2 <- function(r) M(y, .sub(M(p, r$F3), M(r$s, F1)))
  } else {
    r1 <- sym_match(.sub(.ode_d(x, t), eq[[1L]]), .sub(F2, F3))
    pat2 <- function(r) .sub(M(p, r$F3), M(r$s, F1))
  }
  r <- .ode_3eq_match(x, y, z, t, eq, r1, pat2, F1, F2, F3, p, q, s)
  a <- r$p; b <- r$q; c <- r$s
  sx <- function(e) subs(subs(subs(e, x, u), y, v), z, w)
  F1 <- sx(r$F1); F2 <- sx(r$F2); F3 <- sx(r$F3)
  fu <- function(a1) fcall("u", list(a1)); fv <- function(a1) fcall("v", list(a1)); fw <- function(a1) fcall("w", list(a1))
  if (type5) {
    x_yz <- .pow(M(C1, .pow(v, .neg(b)), .pow(w, .neg(c))), .neg(a))
    y_zx <- .pow(M(C1, .pow(w, .neg(c)), .pow(u, .neg(a))), .neg(b))
    z_xy <- .pow(M(C1, .pow(u, .neg(a)), .pow(v, .neg(b))), .neg(c))
    A <- M(v, .sub(M(a, F3), M(c, F1))); B <- M(w, .sub(M(b, F1), M(a, F2))); D <- M(u, .sub(M(c, F2), M(b, F3)))
  } else {
    z_xy <- .div(.sub(.sub(C1, M(a, u)), M(b, v)), c)
    y_zx <- .div(.sub(.sub(C1, M(a, u)), M(c, w)), b)
    x_yz <- .div(.sub(.sub(C1, M(b, v)), M(c, w)), a)
    A <- .sub(M(a, F3), M(c, F1)); B <- .sub(M(b, F1), M(a, F2)); D <- .sub(M(c, F2), M(b, F3))
  }
  rhs <- function(e) dsolve(e)$args[[2L]]
  sd <- function(e, a1, b1) subs(e, a1, b1)
  y_x <- rhs(.sub(.ode_d(fv(u), u), subs(sd(.div(A, D), w, z_xy), v, fv(u))))
  z_x <- rhs(.sub(.ode_d(fw(u), u), subs(sd(.div(B, D), v, y_zx), w, fw(u))))
  z_y <- rhs(.sub(.ode_d(fw(v), v), subs(sd(.div(B, A), u, x_yz), w, fw(v))))
  x_y <- rhs(.sub(.ode_d(fu(v), v), subs(sd(.div(D, A), w, z_xy), u, fu(v))))
  y_z <- rhs(.sub(.ode_d(fv(w), w), subs(sd(.div(A, B), u, x_yz), v, fv(w))))
  x_z <- rhs(.sub(.ode_d(fu(w), w), subs(sd(.div(D, B), v, y_zx), u, fu(w))))
  sol1 <- rhs(.sub(.ode_d(fu(t), t), subs(subs(subs(D, v, y_x), w, z_x), u, fu(t))))
  sol2 <- rhs(.sub(.ode_d(fv(t), t), subs(subs(subs(A, u, x_y), w, z_y), v, fv(t))))
  sol3 <- rhs(.sub(.ode_d(fw(t), t), subs(subs(subs(B, u, x_z), v, y_z), w, fw(t))))
  list(sol1, sol2, sol3)
}

# SymPy's type4 solver calls a Symbol, which fails after the matching
.ode_nonlinear_3eq_order1_type4 <- function(x, y, z, t, eq) {
  .ode_get_numbered_constants(eq, num = 1L)
  .ode_stop("TypeError", "'Symbol' object is not callable")
}
