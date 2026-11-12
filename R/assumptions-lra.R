# Linear real arithmetic for ask(), after sympy/assumptions/lra_satask.py and
# the parts of sympy/logic/algorithms/lra_theory.py that decide whether
# linear constraints over the rationals can hold together. Where SymPy runs
# a simplex inside its DPLL search, this runs a DPLL search that checks the
# asserted constraints by Fourier and Motzkin elimination; both decide the
# same satisfiability.

.unhandled_input <- function(msg) stop(.sympy_condition("UnhandledInput", msg))

# ALLOWED_PRED and the WHITE_LIST of lra_satask.py (which names
# extended_negative twice and so leaves extended_nonnegative out)
.LRA_ALLOWED <- c("eq", "gt", "lt", "le", "ge")
.LRA_WHITE_LIST <- c(.LRA_ALLOWED, "positive", "negative", "zero", "nonzero", "nonpositive", "nonnegative",
                     "extended_positive", "extended_negative", "extended_nonpositive",
                     "extended_nonzero", "negative_infinite", "positive_infinite")
.LRA_POS_NEG_ZERO <- list(extended_positive = "positive", extended_negative = "negative",
                          extended_nonpositive = "nonpositive", extended_nonzero = "nonzero",
                          negative_infinite = FALSE, positive_infinite = FALSE)
.LRA_TO_BINREL <- c(positive = "gt", negative = "lt", zero = "eq", nonpositive = "le",
                    nonnegative = "ge", nonzero = "ne")

.is_builtin_applied <- function(p) {
  f <- .ap_fun(p)
  f$type == "Predicate" && !isTRUE(f$extra$custom)
}

.binrel <- function(name, a, b) new_sym("AppliedBinaryRelation", list(.qp(name), a, b))

# _pred_to_binrel: FALSE for the infinite predicates
.lra_pred_to_binrel <- function(pred) {
  if (!.is_applied(pred) || !.is_builtin_applied(pred)) return(pred)
  nm <- .ap_name(pred)
  arg <- .ap_args(pred)[[1L]]
  if (nm %in% names(.LRA_POS_NEG_ZERO)) {
    f <- .LRA_POS_NEG_ZERO[[nm]]
    if (isFALSE(f)) return(FALSE)
    nm <- f
  }
  if (nm %in% names(.LRA_TO_BINREL)) return(.binrel(.LRA_TO_BINREL[[nm]], arg, S_Zero))
  pred
}

# extract_pred_from_old_assum
.lra_old_assumptions <- function(all_exprs) {
  ret <- list()
  for (expr in all_exprs) {
    if (!length(free_symbols(expr))) next
    if (!isTRUE(.fact(expr, "real"))) .unhandled_input(paste0("LRASolver: ", .str(expr), " must be real"))
    if (expr$type == "Mul" && any(vapply(expr$args, function(a) !isTRUE(.fact(a, "real")), TRUE)))
      .unhandled_input(paste0("LRASolver: ", .str(expr), " must be real"))
    if (isTRUE(.fact(expr, "integer")) && !isTRUE(.fact(expr, "zero")))
      .unhandled_input(paste0("LRASolver: ", .str(expr), " is an integer"))
    if (identical(.fact(expr, "integer"), FALSE))
      .unhandled_input(paste0("LRASolver: ", .str(expr), " can't be an integer"))
    if (identical(.fact(expr, "rational"), FALSE))
      .unhandled_input(paste0("LRASolver: ", .str(expr), " is irational"))
    for (f in c("zero", "positive", "negative", "nonzero", "nonpositive", "nonnegative")) {
      if (isTRUE(.fact(expr, f))) {
        ret[[length(ret) + 1L]] <- .QA(f, expr)
        break
      }
    }
  }
  ret
}

# _preprocess: clauses over Q.eq, Q.gt, Q.lt, Q.ge and Q.le atoms, with
# x != y as x > y | x < y and a negated x != y as x = y. Returns the
# clauses and the atom of each new variable.
.lra_preprocess <- function(enc) {
  keys <- character()
  props <- list()
  id <- function(prop) {
    k <- if (isFALSE(prop)) "\rFalse" else prop$key
    i <- match(k, keys)
    if (is.na(i)) {
      keys <<- c(keys, k)
      props[[length(props) + 1L]] <<- prop
      i <- length(keys)
    }
    i
  }
  data <- lapply(enc$data, function(clause) {
    out <- integer()
    for (lit in clause) {
      if (lit == 0L) {
        out <- c(out, 0L)
        next
      }
      prop <- .lra_pred_to_binrel(enc$symbols[[abs(lit)]])
      negated <- lit < 0L
      sgn <- if (negated) -1L else 1L
      if (isFALSE(prop) || !.is_applied(prop)) {
        out <- c(out, sgn * id(prop))
        next
      }
      nm <- .ap_name(prop)
      a <- .ap_args(prop)
      if (negated && nm == "eq") {
        negated <- FALSE
        nm <- "ne"
      }
      if (nm == "ne") {
        if (negated) out <- c(out, id(.binrel("eq", a[[1L]], a[[2L]])))
        else out <- c(out, id(.binrel("gt", a[[1L]], a[[2L]])), id(.binrel("lt", a[[1L]], a[[2L]])))
        next
      }
      out <- c(out, sgn * id(prop))
    }
    unique(out)
  })
  list(data = data, props = props)
}

.is_rational_number <- function(e) is_sym(e) && e$type %in% c("Integer", "Rational")

# _sep_const_terms and _sep_const_coeff
.lra_sep_terms <- function(e) {
  terms <- if (e$type == "Add") .sympy_args(e) else list(e)
  const <- vapply(terms, function(t) !length(free_symbols(t)), TRUE)
  list(vars = .add(terms[!const]), const = .add(terms[const]))
}
.lra_sep_coeff <- function(e) {
  if (e$type == "Add") return(list(var = e, coeff = S_One))
  fs <- if (e$type == "Mul") .sympy_args(e) else list(e)
  const <- vapply(fs, function(f) !length(free_symbols(f)), TRUE)
  list(var = .mul(fs[!const]), coeff = .mul(fs[const]))
}

# LRASolver.from_encoded_cnf: the linear constraint of each atom, and the
# unit clauses of atoms that are true or false outright. A constraint is
# list(coef = named bigq coefficients over the atom terms, bound, upper,
# strict, eq); it reads sum(coef * term) <= bound (or >=, < or >).
.lra_atoms <- function(props) {
  cons <- vector("list", length(props))
  conflicts <- list()
  terms <- list()
  for (atom_id in seq_along(props)) {
    prop <- props[[atom_id]]
    if (isFALSE(prop) || (is_sym(prop) && .is_false(prop))) {
      conflicts[[length(conflicts) + 1L]] <- -atom_id
      next
    }
    if (!.is_applied(prop)) {
      if (.is_true(prop)) {
        conflicts[[length(conflicts) + 1L]] <- atom_id
        next
      }
      .ask_value_error(paste0("Unhandled Predicate: ", .str(prop)))
    }
    nm <- .ap_name(prop)
    a <- .ap_args(prop)
    lhs <- a[[1L]]; rhs <- a[[2L]]
    if (.is_nan_sym(lhs) || .is_nan_sym(rhs)) .ask_value_error(paste0(.str(prop), " contains nan"))
    if (isTRUE(.fact(lhs, "imaginary")) || isTRUE(.fact(rhs, "imaginary")))
      .unhandled_input(paste0(.str(prop), " contains an imaginary component"))
    if (.keq(lhs, S_Infinity) || .keq(rhs, S_Infinity))
      .unhandled_input(paste0(.str(prop), " contains infinity"))
    expr <- .sub(lhs, rhs)
    rel <- switch(nm, eq = Eq(expr, S_Zero), gt = Gt(expr, S_Zero), lt = Lt(expr, S_Zero),
                  le = Le(expr, S_Zero), ge = Ge(expr, S_Zero))
    if (.is_true(rel)) {
      conflicts[[length(conflicts) + 1L]] <- atom_id
      next
    }
    if (.is_false(rel)) {
      conflicts[[length(conflicts) + 1L]] <- -atom_id
      next
    }
    if (!length(free_symbols(expr))) .unhandled_input(paste0(.str(prop), " could not be simplified"))
    if (nm %in% c("ge", "gt")) expr <- .neg(expr)
    st <- .lra_sep_terms(expr)
    sc <- .lra_sep_coeff(st$vars)
    vars <- sc$var
    var_coeff <- sc$coeff
    const <- .div(st$const, var_coeff)
    tl <- if (vars$type == "Add") .sympy_args(vars) else list(vars)
    coef <- list()
    for (t in tl) {
      tc <- .lra_sep_coeff(t)
      k <- tc$var$key
      terms[[k]] <- tc$var
      coef[[k]] <- tc$coeff
    }
    if (length(tl) == 1L) coef[[1L]] <- S_One
    bound <- .neg(const)
    cons[[atom_id]] <- list(coef = coef, bound = bound, eq = nm == "eq",
                            upper = isTRUE(.is_true(Gt(var_coeff, S_Zero))), strict = nm %in% c("gt", "lt"))
  }
  # the atom terms must not share symbols
  fs <- lapply(terms, function(t) .keys_of(free_symbols(t)))
  if (length(fs) && length(unique(unlist(fs))) < length(unlist(fs)))
    .unhandled_input("Nonlinearity is not handled")
  for (cn in cons) {
    if (is.null(cn)) next
    if (length(cn$coef) > 1L && !all(vapply(cn$coef, .is_rational_number, TRUE)))
      .unhandled_input("Non-rational numbers are not handled")
  }
  for (cn in cons) {
    if (!is.null(cn) && !.is_rational_number(cn$bound)) .unhandled_input("Non-rational numbers are not handled")
  }
  list(cons = cons, conflicts = conflicts, terms = names(terms))
}

# Fourier and Motzkin elimination: whether a %*% x <= b (< where strict)
# has a rational solution. a is a list of bigq coefficient vectors.
.fm_feasible <- function(a, b, strict) {
  if (!length(a)) return(TRUE)
  nv <- length(a[[1L]])
  for (j in seq_len(nv)) {
    cj <- do.call(c, lapply(a, function(r) r[j]))
    pos <- which(cj > 0); neg <- which(cj < 0); zero <- which(cj == 0)
    na <- a[zero]; nb <- b[zero]; ns <- strict[zero]
    for (p in pos) for (n in neg) {
      na[[length(na) + 1L]] <- a[[p]] / cj[p] - a[[n]] / cj[n]
      nb <- c(nb, b[p] / cj[p] - b[n] / cj[n])
      ns <- c(ns, strict[p] || strict[n])
    }
    a <- na; b <- nb; strict <- ns
    if (!length(a)) return(TRUE)
  }
  !any(b < 0 | (strict & b == 0))
}

# The constraints asserted by an assignment: a true atom asserts its
# constraint, a false one its negation, except a false equality.
.lra_check <- function(info, val) {
  terms <- info$terms
  nv <- length(terms)
  a <- list(); b <- gmp::as.bigq(integer()); strict <- logical()
  add <- function(coef, bound, upper, st) {
    v <- gmp::as.bigq(rep(0L, nv))
    for (k in names(coef)) v[match(k, terms)] <- gmp::as.bigq(coef[[k]]$value)
    bq <- gmp::as.bigq(bound$value)
    if (!upper) { v <- -v; bq <- -bq }
    a[[length(a) + 1L]] <<- v
    b <<- c(b, bq)
    strict <<- c(strict, st)
  }
  for (i in which(val != 0L)) {
    cn <- info$cons[[i]]
    if (is.null(cn)) next
    if (cn$eq) {
      if (val[i] > 0L) {
        add(cn$coef, cn$bound, TRUE, FALSE)
        add(cn$coef, cn$bound, FALSE, FALSE)
      }
      next
    }
    neg <- val[i] < 0L
    add(cn$coef, cn$bound, cn$upper != neg, cn$strict != neg)
  }
  .fm_feasible(a, b, strict)
}

# satisfiable(..., use_lra_theory=True)
.lra_satisfiable <- function(pre) {
  if (any(vapply(pre$data, function(cl) length(cl) == 1L && cl == 0L, TRUE))) return(FALSE)
  info <- .lra_atoms(pre$props)
  clauses <- c(lapply(pre$data, function(cl) cl[cl != 0L]), info$conflicts)
  if (any(lengths(clauses) == 0L)) return(FALSE)
  n <- length(pre$props)
  val <- integer(n)
  if (!length(clauses)) return(.lra_check(info, val))
  f <- .sat_flatten(clauses)
  search <- function(val) {
    val <- .sat_propagate(f$lits, f$cid, f$ncl, val)
    if (is.null(val) || !.lra_check(info, val)) return(FALSE)
    lv <- sign(f$lits) * val[abs(f$lits)]
    sat <- tabulate(f$cid[lv == 1L], nbins = f$ncl) > 0L
    if (all(sat)) return(TRUE)
    cand <- f$lits[lv == 0L & !sat[f$cid]]
    tab <- table(cand)
    l <- as.integer(names(tab)[which.max(tab)])
    v1 <- val
    v1[abs(l)] <- sign(l)
    if (search(v1)) return(TRUE)
    v1[abs(l)] <- -sign(l)
    search(v1)
  }
  search(val)
}

#' Decide a proposition with linear real arithmetic
#'
#' `lra_satask()` decides propositions made of the relations `Q$eq`,
#' `Q$ne`, `Q$gt`, `Q$ge`, `Q$lt` and `Q$le` (and the sign predicates such
#' as `Q$positive`) between linear expressions of real symbols, under
#' assumptions of the same kind. It is the last step of [ask()]. Input it
#' does not handle (nonlinear terms, symbols that are not known to be real
#' or that are integers, other predicates) raises an `UnhandledInput`
#' condition.
#'
#' @param proposition A boolean expression of relations.
#' @param assumptions Assumptions, or `TRUE` for none.
#' @return `TRUE`, `FALSE` or `NA`; inconsistent assumptions raise a
#'   `ValueError`.
#' @examples
#' x <- Symbol("x", real = TRUE)
#' lra_satask(Q$gt(x, 0), Q$gt(x, 1))
#' lra_satask(Q$eq(x, 1), !Q$ne(x, 0))
#' @export
lra_satask <- function(proposition, assumptions = TRUE) {
  proposition <- .ask_bool(proposition)
  assumptions <- .ask_bool(assumptions)
  props <- .CNF_from_prop(proposition)
  nprops <- .CNF_from_prop(Not(proposition))
  enc <- .enc_cnf()
  .enc_from_cnf(enc, .CNF_from_prop(assumptions))
  .lra_check_satisfiability(props, nprops, enc)
}

# check_satisfiability
.lra_check_satisfiability <- function(prop, nprop, factbase) {
  sat_true <- .enc_copy(factbase)
  sat_false <- .enc_copy(factbase)
  .enc_add_from_cnf(sat_true, prop)
  .enc_add_from_cnf(sat_false, nprop)
  all_pred <- Filter(.is_applied, sat_true$symbols)
  all_exprs <- .pyset(unlist(lapply(all_pred, .ap_args), recursive = FALSE))
  for (pred in all_pred) {
    if (!.is_builtin_applied(pred) || !(.ap_name(pred) %in% c(.LRA_WHITE_LIST, "ne")))
      .unhandled_input(paste0("LRASolver: ", .str(pred), " is an unhandled predicate"))
  }
  for (expr in all_exprs) {
    if (.is_matrix_kind(expr)) .unhandled_input(paste0("LRASolver: ", .str(expr), " is of MatrixKind"))
    if (.is_nan_sym(expr)) .unhandled_input("LRASolver: nan")
  }
  for (assm in .lra_old_assumptions(all_exprs)) {
    for (enc in list(sat_true, sat_false)) {
      v <- .enc_arg(enc, list(lit = assm, neg = FALSE))
      enc$data[[length(enc$data) + 1L]] <- v
    }
  }
  can_be_true <- .lra_satisfiable(.lra_preprocess(sat_true))
  can_be_false <- .lra_satisfiable(.lra_preprocess(sat_false))
  if (can_be_true && can_be_false) return(NA)
  if (can_be_true) return(TRUE)
  if (can_be_false) return(FALSE)
  .ask_value_error("Inconsistent assumptions")
}
