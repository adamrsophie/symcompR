# The binary relations Q.eq, Q.ne, Q.gt, Q.ge, Q.lt and Q.le, after
# sympy/assumptions/relation/binrel.py and equality.py, with the functions
# is_eq() and is_ge() of sympy/core/relational.py taking assumptions as
# sympy/assumptions/wrapper.py does: without assumptions an expression's
# own assumptions answer, with them ask() does.

.BINREL_REVERSED <- c(eq = "eq", ne = "ne", gt = "lt", lt = "gt", ge = "le", le = "ge")
.BINREL_NEGATED <- c(eq = "ne", ne = "eq", gt = "le", ge = "lt", lt = "ge", le = "gt")

.binrel_reversed <- function(e) {
  a <- .ap_args(e)
  new_sym("AppliedBinaryRelation", list(.qp(.BINREL_REVERSED[[.ap_name(e)]]), a[[2L]], a[[1L]]))
}
.binrel_negated <- function(e) new_sym("AppliedBinaryRelation", c(list(.qp(.BINREL_NEGATED[[.ap_name(e)]])), .ap_args(e)))

# AppliedBinaryRelation._eval_ask
.binrel_eval_ask <- function(e, assumptions) {
  conj <- .ask_conjuncts(assumptions)
  keys <- vapply(conj, function(a) if (is_relational(a)) .rel_to_binrel(a)$key else a$key, "")
  rev <- .binrel_reversed(e)
  if (e$key %in% keys || rev$key %in% keys) return(TRUE)
  negs <- list(.binrel_negated(e), .binrel_negated(rev), new_sym("Not", list(e)), new_sym("Not", list(rev)))
  if (any(vapply(negs, .k, "") %in% keys)) return(FALSE)
  fn <- .ap_fun(e)
  args <- .ap_args(e)
  ret <- .pred_eval(fn, args, assumptions)
  if (!is.na(ret)) return(ret)
  args <- lapply(args, function(a) tryCatch(simplify(a), error = function(err) a))
  .pred_eval(fn, args, assumptions)
}

# EqualityPredicate.eval and the others: is_eq, is_neq, is_gt, ...
.binrel_eval <- function(name, args, assumptions) {
  a <- if (.is_true(.ask_bool(assumptions))) NULL else assumptions
  lhs <- args[[1L]]; rhs <- args[[2L]]
  switch(name,
    eq = .is_eq_a(lhs, rhs, a),
    ne = .fuzzy_not(.is_eq_a(lhs, rhs, a)),
    ge = .is_ge_a(lhs, rhs, a),
    le = .is_ge_a(rhs, lhs, a),
    gt = .fuzzy_not(.is_ge_a(rhs, lhs, a)),
    lt = .fuzzy_not(.is_ge_a(lhs, rhs, a)))
}

# AssumptionsWrapper(e, assumptions).is_<fact>
.awf <- function(e, fact, assumptions) {
  if (is.null(assumptions)) return(.fact(e, fact))
  ask(.QA(fact, e), assumptions)
}

.is_expr <- function(e) is_sym(e) && !.is_boolean_kind(e) && !.is_predicate(e) && !is_set(e) &&
  !(e$type %in% c("Tuple", "Str"))

.fuzzy_xor <- function(vals) {
  t <- 0L
  for (v in vals) {
    v <- .ask_fz(v)
    if (is.na(v)) return(NA)
    if (v) t <- t + 1L
  }
  t %% 2L == 1L
}

# Expr.is_comparable
.ask_is_comparable <- function(e) {
  if (!is_sym(e) || !.is_number(e)) return(FALSE)
  if (identical(.fact(e, "extended_real"), FALSE)) return(FALSE)
  if (is_exact_number(e) || e$type == "Float") return(TRUE)
  if (e$type == "Constant") return(!(e$name %in% c("nan", "zoo", "I")))
  ri <- .as_real_imag(e)
  if (!.py_falsy(ri[[2L]])) return(FALSE)
  n <- ri[[1L]]
  if (is_exact_number(n) || n$type == "Float") return(TRUE)
  !is.na(.ev2_sign(n))
}

# _n2: the sign of (a - b).evalf(2) when a and b are comparable, else NULL
.n2 <- function(a, b) {
  if (!(.ask_is_comparable(a) && .ask_is_comparable(b))) return(NULL)
  dif <- .sub(a, b)
  if (.is_nan_sym(dif)) return(NULL)
  s <- .ev2_sign(dif)
  if (is.na(s)) NULL else s
}

# is_ge(lhs, rhs, assumptions)
.is_ge_a <- function(lhs, rhs, assumptions) {
  if (!(.is_expr(lhs) && .is_expr(rhs))) .ask_type_error("Can only compare inequalities with Expr")
  n2 <- .n2(lhs, rhs)
  if (!is.null(n2)) return(n2 >= 0)
  W <- function(e, f) .awf(e, f, assumptions)
  if (isTRUE(W(lhs, "extended_real")) && isTRUE(W(rhs, "extended_real"))) {
    if ((isTRUE(W(lhs, "infinite")) && isTRUE(W(lhs, "extended_positive"))) ||
        (isTRUE(W(rhs, "infinite")) && isTRUE(W(rhs, "extended_negative")))) return(TRUE)
    diff <- .sub(lhs, rhs)
    if (!.is_nan_sym(diff)) {
      rv <- W(diff, "extended_nonnegative")
      if (!is.na(rv)) return(rv)
    }
  }
  NA
}

.add_terms <- function(e) if (is_sym(e) && e$type == "Add") .sympy_args(e) else list(e)

.as_coeff_Add <- function(e) {
  if (e$type %in% c(.NUMBER_TYPES)) return(list(e, S_Zero))
  if (e$type == "Add") {
    num <- vapply(e$args, function(a) a$type %in% .NUMBER_TYPES, TRUE)
    if (any(num)) return(list(e$args[[which(num)[1L]]], .add(e$args[-which(num)[1L]])))
  }
  list(S_Zero, e)
}

# is_eq(lhs, rhs, assumptions)
.is_eq_a <- function(lhs, rhs, assumptions) {
  if (is.null(assumptions) || is_set(lhs) || is_set(rhs)) {
    # without assumptions this is the evaluation of Eq(lhs, rhs)
    r <- tryCatch(.rel_eval("Equality", lhs, rhs), error = function(err) NULL)
    if (is.null(r)) return(NA)
    return(.is_true(r))
  }
  if (.keq(lhs, rhs)) return(TRUE)
  if (lhs$type == "BooleanAtom" && rhs$type == "BooleanAtom") return(FALSE)
  if (!(.is_symbol_node(lhs) || .is_symbol_node(rhs)) && .is_boolean_kind(lhs) != .is_boolean_kind(rhs)) return(FALSE)
  W <- function(e, f) .awf(e, f, assumptions)
  linf <- W(lhs, "infinite"); rinf <- W(rhs, "infinite")
  if (isTRUE(linf) || isTRUE(rinf)) {
    if (isTRUE(.fuzzy_xor(list(linf, rinf)))) return(FALSE)
    lr <- W(lhs, "extended_real"); rr <- W(rhs, "extended_real")
    if (isTRUE(.fuzzy_xor(list(lr, rr)))) return(FALSE)
    if (isTRUE(.fuzzy_and_vals(list(lr, rr)))) {
      if (isTRUE(.fuzzy_and_vals(list(linf, rinf))))
        return(.fuzzy_xor(list(W(lhs, "extended_positive"), .fuzzy_not(W(rhs, "extended_positive")))))
      if ((isTRUE(W(lhs, "extended_positive")) && isTRUE(W(rhs, "extended_nonpositive"))) ||
          (isTRUE(W(lhs, "extended_negative")) && isTRUE(W(rhs, "extended_nonnegative"))) ||
          (isTRUE(W(lhs, "extended_nonpositive")) && isTRUE(W(rhs, "extended_positive"))) ||
          (isTRUE(W(lhs, "extended_nonnegative")) && isTRUE(W(rhs, "extended_negative")))) return(FALSE)
      return(NA)
    }
    split_ri <- function(e) {
      out <- list(real = list(), imag = list(), none = list())
      for (t in .add_terms(e)) {
        g <- if (isTRUE(W(t, "extended_real"))) "real"
             else if (isTRUE(W(.mul(list(S_I, t)), "extended_real"))) "imag" else "none"
        out[[g]][[length(out[[g]]) + 1L]] <- t
      }
      out
    }
    lri <- split_ri(lhs)
    if (!length(lri$none)) {
      rri <- split_ri(rhs)
      if (!length(rri$none)) {
        eq_real <- .is_eq_a(.add(lri$real), .add(rri$real), assumptions)
        eq_imag <- .is_eq_a(.mul(list(S_I, .add(lri$imag))), .mul(list(S_I, .add(rri$imag))), assumptions)
        return(.fuzzy_and_vals(list(eq_real, eq_imag)))
      }
    }
    al <- fcall("arg", list(lhs)); ar <- fcall("arg", list(rhs))
    if (!(.is_nan_sym(al) && .is_nan_sym(ar))) return(.is_eq_a(al, ar, assumptions))
  }
  if (.is_expr(lhs) && .is_expr(rhs)) {
    if (!(isTRUE(W(lhs, "finite")) && isTRUE(W(rhs, "finite")))) {
      cd <- function(e) {
        out <- list()
        for (t in .add_terms(e)) {
          sp <- as_coeff_mul_split(t)
          k <- sp[[2L]]$key
          out[[k]] <- list(sp[[2L]], if (is.null(out[[k]])) sp[[1L]] else num_add(out[[k]][[2L]], sp[[1L]]))
        }
        out
      }
      L <- cd(lhs); R <- cd(rhs)
      for (k in intersect(names(L), names(R))) {
        if (k == S_One$key || !.keq(L[[k]][[2L]], R[[k]][[2L]])) next
        if (!isTRUE(W(L[[k]][[1L]], "finite"))) return(NA)
      }
    }
    dif <- .sub(lhs, rhs)
    z <- W(dif, "zero")
    if (!is.na(z)) {
      if (!z && isTRUE(W(dif, "commutative"))) return(FALSE)
      if (z) return(TRUE)
    }
    ct <- .as_coeff_Add(dif)
    cf <- ct[[1L]]; tt <- ct[[2L]]
    if (cf$type == "Float") {
      q <- .mp_to_q(.float_mp(cf))
      if (gmp::denominator(q) == 1) {
        if (identical(.fact(tt, "integer"), FALSE)) return(FALSE)
      } else if (identical(.fact(tt, "rational"), FALSE)) return(FALSE)
    }
    n2 <- .n2(lhs, rhs)
    if (!is.null(n2)) return(n2 == 0)
    nd <- .as_numer_denom(dif)
    n <- nd[[1L]]; d <- nd[[2L]]
    rv <- NA
    if (isTRUE(W(n, "zero"))) {
      rv <- W(d, "nonzero")
    } else if (isTRUE(W(n, "finite"))) {
      if (isTRUE(W(d, "infinite"))) rv <- TRUE
      else if (identical(W(n, "zero"), FALSE)) rv <- W(d, "infinite")
    } else if (.py_any(length(.add_terms(n)), function(i) W(.add_terms(n)[[i]], "infinite"))) {
      rv <- FALSE
    }
    if (!is.na(rv)) return(rv)
  }
  NA
}
