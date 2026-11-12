# Piecewise functions (sympy/functions/elementary/piecewise.py) and
# rewriting between function families (Expr.rewrite).

#' Piecewise defined functions
#'
#' `Piecewise(list(e1, c1), list(e2, c2), ...)` takes the value of the first
#' expression whose condition holds. Conditions that are `FALSE` are dropped
#' and a `TRUE` condition ends the list. `piecewise_fold()` moves arithmetic
#' inside the pieces.
#'
#' @param ... Pairs `list(expr, condition)`.
#' @param e An expression.
#' @param evaluate Evaluate the folded `Piecewise`.
#' @return A `sym` expression.
#' @examples
#' x <- Symbol("x")
#' p <- Piecewise(list(-x, x < 0), list(x, TRUE))
#' subs(p, x, -3)
#' diff(p, x)
#' @export
Piecewise <- function(...) {
  pairs <- list(...)
  pairs <- lapply(pairs, function(p) list(as_sym(p[[1L]]), .as_bool(p[[2L]])))
  .piecewise(pairs)
}

# Piecewise.__new__ and Piecewise.eval
.piecewise <- function(pairs) {
  args <- list()
  for (p in pairs) {
    if (.is_false(p[[2L]])) next
    args[[length(args) + 1L]] <- p
    if (.is_true(p[[2L]])) break
  }
  if (!length(args)) return(S_NaN)
  if (length(args) == 1L && .is_true(args[[1L]][[2L]])) return(args[[1L]][[1L]])
  newargs <- .pw_collapse(args)
  same <- length(newargs) == length(args) && all(vapply(seq_along(args), function(i)
    identical(args[[i]][[1L]]$key, newargs[[i]][[1L]]$key) && identical(args[[i]][[2L]]$key, newargs[[i]][[2L]]$key), TRUE))
  if (!same) return(.piecewise(newargs))
  new_sym("Piecewise", lapply(args, function(p) new_sym("ExprCondPair", list(p[[1L]], p[[2L]]))))
}

.REL_NEG_OP <- c("<" = ">=", "<=" = ">", ">" = "<=", ">=" = "<", "==" = "!=", "!=" = "==")
.rel_negated_canonical <- function(r) .rel_canonical(.rel_make(.REL_NEG_OP[[.rel_op(r)]], r$args[[1L]], r$args[[2L]]))

# _piecewise_collapse_arguments
.pw_collapse <- function(args) {
  newargs <- list(); current <- character()
  has_key <- function(c) c$key %in% current
  for (p in args) {
    expr <- p[[1L]]
    cond <- .map_bottom_up(p[[2L]], function(z) if (is_relational(z)) .canonical_coeff(z) else z)
    if (expr$type == "Piecewise") {
      unmatching <- list()
      for (ec in expr$args) {
        e <- ec$args[[1L]]; c <- ec$args[[2L]]
        if (has_key(c)) next
        if (identical(c$key, cond$key)) {
          if (!.is_true(c)) expr <- if (length(unmatching)) .piecewise(c(unmatching, list(list(e, c)))) else e
          break
        }
        unmatching[[length(unmatching) + 1L]] <- list(e, c)
      }
    }
    cands <- c(list(cond), if (cond$type == "And") cond$args else list())
    if (any(vapply(cands, has_key, TRUE))) next
    if (cond$type == "And") {
      nonredundant <- list(); falsified <- FALSE
      for (c in cond$args) {
        if (is_relational(c)) {
          if (has_key(.rel_negated_canonical(c))) next
          if (c$type %in% c("StrictLessThan", "StrictGreaterThan") &&
              has_key(.rel_make(if (c$type == "StrictLessThan") "<=" else ">=", c$args[[1L]], c$args[[2L]]))) {
            cond <- S_false; falsified <- TRUE; break
          }
        }
        nonredundant[[length(nonredundant) + 1L]] <- c
      }
      if (!falsified) cond <- do.call(And, nonredundant)
    } else if (is_relational(cond)) {
      if (has_key(.rel_negated_canonical(cond))) cond <- S_true
    }
    current <- c(current, cond$key)
    if (length(newargs)) {
      last <- newargs[[length(newargs)]]
      if (identical(last[[1L]]$key, expr$key)) {
        orcond <- Or(cond, last[[2L]])
        if (orcond$type %in% c("And", "Or")) orcond <- .distribute(orcond, "And", "Or")
        newargs[[length(newargs)]] <- list(expr, orcond)
        next
      }
      if (identical(last[[2L]]$key, cond$key)) next
    }
    newargs[[length(newargs) + 1L]] <- list(expr, cond)
  }
  newargs
}

register_rebuilder("Piecewise", function(e, args) .piecewise(lapply(args, function(t) list(t$args[[1L]], t$args[[2L]]))))
register_rebuilder("ExprCondPair", function(e, args) new_sym("ExprCondPair", list(args[[1L]], .as_bool(args[[2L]]))))
register_str_printer("Piecewise", function(e) {
  paste0("Piecewise(", paste(vapply(.pw_pairs(e), function(p) paste0("(", .str(p[[1L]]), ", ", .str(p[[2L]]), ")"), ""),
                            collapse = ", "), ")")
})
register_neval("Piecewise", function(e, env) {
  for (p in .pw_pairs(e)) if (.neval(p[[2L]], env) != 0) return(.neval(p[[1L]], env))
  NaN
})
register_neval("BooleanAtom", function(e, env) as.numeric(e$name == "true"))
register_diff_hook("Piecewise", function(e, x) {
  .piecewise(lapply(.pw_pairs(e), function(p) list(.diff1(p[[1L]], x), p[[2L]])))
})
.facts_hooks$Piecewise <- function(e) {
  fs <- lapply(.pw_pairs(e), function(p) .facts(p[[1L]]))
  out <- c(commutative = .all_true(fs, "commutative"))
  for (nm in c("real", "positive", "nonnegative", "integer", "finite", "negative", "nonpositive"))
    if (.all_true(fs, nm)) out[nm] <- TRUE
  out
}

#' @rdname Piecewise
#' @export
piecewise_fold <- function(e, evaluate = TRUE) {
  e <- as_sym(e)
  if (!.has_type(e, "Piecewise")) return(e)
  new_args <- list()
  if (e$type %in% c("Piecewise", "ExprCondPair")) {
    pairs <- if (e$type == "Piecewise") .pw_pairs(e) else list(list(e$args[[1L]], e$args[[2L]]))
    for (p in pairs) {
      ex <- p[[1L]]; c <- p[[2L]]
      if (ex$type != "Piecewise") ex <- piecewise_fold(ex)
      if (ex$type == "Piecewise") {
        for (q in .pw_pairs(ex)) new_args[[length(new_args) + 1L]] <- list(piecewise_fold(q[[1L]]), And(q[[2L]], c))
      } else new_args[[length(new_args) + 1L]] <- list(ex, c)
    }
  } else {
    args <- if (e$type %in% c("Add", "Mul")) .sympy_args(e) else e$args
    if (e$type == "Add" || (e$type == "Mul" && .is_comm(e))) {
      isp <- vapply(args, function(a) a$type == "Piecewise", NA)
      p <- args[isp]; args <- args[!isp]
      ck <- vapply(p, function(pw) paste(vapply(.pw_pairs(pw), function(q) q[[2L]]$key, ""), collapse = "|"), "")
      for (k in unique(ck)) {
        grp <- p[ck == k]
        if (length(grp) > 1L) {
          pargs <- lapply(grp, .pw_pairs)
          conds <- lapply(pargs, function(a) vapply(a, function(q) q[[2L]]$key, ""))
          n <- 0L
          while (n < min(lengths(conds)) && length(unique(vapply(conds, `[`, "", n + 1L))) == 1L) n <- n + 1L
          collected <- list()
          for (i in seq_len(n)) collected[[i]] <- list(sym_rebuild(e, lapply(pargs, function(a) a[[i]][[1L]])), pargs[[1L]][[i]][[2L]])
          remains <- list()
          for (a in pargs) {
            if (n == length(a)) next
            if (.is_true(a[[n + 1L]][[2L]])) remains[[length(remains) + 1L]] <- a[[n + 1L]][[1L]]
            else remains[[length(remains) + 1L]] <- .piecewise_noeval(a[(n + 1L):length(a)])
          }
          if (length(remains)) collected[[length(collected) + 1L]] <- list(sym_rebuild(e, remains), S_true)
          args[[length(args) + 1L]] <- .piecewise_noeval(collected)
        } else args <- c(args, grp)
      }
    }
    folded <- lapply(args, piecewise_fold)
    opts <- lapply(folded, function(a) if (a$type == "Piecewise") .pw_pairs(a) else list(list(a, S_true)))
    idx <- as.matrix(rev(expand.grid(lapply(rev(opts), seq_along))))
    for (r in seq_len(nrow(idx))) {
      ec <- Map(function(o, j) o[[j]], opts, idx[r, ])
      new_args[[length(new_args) + 1L]] <- list(sym_rebuild(e, lapply(ec, `[[`, 1L)), do.call(And, lapply(ec, `[[`, 2L)))
    }
  }
  if (is.null(evaluate)) {
    ks <- vapply(new_args, function(p) p[[2L]]$key, "")
    u <- unique(ks)
    last <- vapply(u, function(k) max(which(ks == k)), 0)
    new_args <- lapply(u[order(last)], function(k) new_args[[which(ks == k)[1L]]])
    rv <- .piecewise_noeval(new_args)
    return(rv)
  }
  if (isTRUE(evaluate)) .piecewise(new_args) else .piecewise_noeval(new_args)
}

# Piecewise(*args, evaluate=False)
.piecewise_noeval <- function(pairs) {
  args <- list()
  for (p in pairs) {
    if (.is_false(p[[2L]])) next
    args[[length(args) + 1L]] <- p
    if (.is_true(p[[2L]])) break
  }
  if (length(args) == 1L && .is_true(args[[1L]][[2L]])) return(args[[1L]][[1L]])
  new_sym("Piecewise", lapply(args, function(p) new_sym("ExprCondPair", list(p[[1L]], p[[2L]]))))
}

# ---------------------------------------------------------------------------
# rewrite

.rewrite_rules <- list(
  exp = list(
    sin = function(a) .div(.sub(fcall("exp", list(.mul(list(S_I, a)))), fcall("exp", list(.mul(list(.neg(S_I), a))))), .mul(list(S_Two, S_I))),
    cos = function(a) .div(.add(list(fcall("exp", list(.mul(list(S_I, a)))), fcall("exp", list(.mul(list(.neg(S_I), a)))))), S_Two),
    tan = function(a) {
      e1 <- fcall("exp", list(.mul(list(S_I, a)))); e2 <- fcall("exp", list(.mul(list(.neg(S_I), a))))
      .div(.mul(list(.neg(S_I), .sub(e1, e2))), .add(list(e1, e2)))
    },
    cot = function(a) {
      e1 <- fcall("exp", list(.mul(list(S_I, a)))); e2 <- fcall("exp", list(.mul(list(.neg(S_I), a))))
      .div(.mul(list(S_I, .add(list(e1, e2)))), .sub(e1, e2))
    },
    sinh = function(a) .div(.sub(fcall("exp", list(a)), fcall("exp", list(.neg(a)))), S_Two),
    cosh = function(a) .div(.add(list(fcall("exp", list(a)), fcall("exp", list(.neg(a))))), S_Two),
    tanh = function(a) .div(.sub(fcall("exp", list(a)), fcall("exp", list(.neg(a)))),
                            .add(list(fcall("exp", list(a)), fcall("exp", list(.neg(a)))))),
    coth = function(a) .div(.add(list(fcall("exp", list(a)), fcall("exp", list(.neg(a))))),
                            .sub(fcall("exp", list(a)), fcall("exp", list(.neg(a))))),
    sech = function(a) .div(S_Two, .add(list(fcall("exp", list(a)), fcall("exp", list(.neg(a)))))),
    csch = function(a) .div(S_Two, .sub(fcall("exp", list(a)), fcall("exp", list(.neg(a))))),
    sec = function(a) .div(S_Two, .add(list(fcall("exp", list(.mul(list(S_I, a)))), fcall("exp", list(.mul(list(.neg(S_I), a))))))),
    csc = function(a) .div(.mul(list(S_Two, S_I)), .sub(fcall("exp", list(.mul(list(S_I, a)))), fcall("exp", list(.mul(list(.neg(S_I), a))))))
  ),
  cos = list(
    sin = function(a) fnode("cos", list(.sub(a, .mul(list(S_Half, S_Pi))))),
    tan = function(a) .div(fnode("cos", list(.sub(a, .mul(list(S_Half, S_Pi))))), fcall("cos", list(a))),
    cot = function(a) .div(fcall("cos", list(a)), fnode("cos", list(.sub(a, .mul(list(S_Half, S_Pi)))))),
    sec = function(a) .inv(fcall("cos", list(a))),
    csc = function(a) .inv(fnode("cos", list(.sub(a, .mul(list(S_Half, S_Pi)))))),
    exp = function(a) {
      ia <- .mul(list(S_I, a))
      .add(list(fcall("cos", list(ia)), .mul(list(S_I, fcall("cos", list(.add(list(ia, .mul(list(S_Half, S_Pi))))))))))
    },
    cosh = function(a) fnode("cos", list(.mul(list(S_I, a))))
  ),
  sin = list(
    cos = function(a) fnode("sin", list(.add(list(a, .mul(list(S_Half, S_Pi)))))),
    tan = function(a) .div(.mul(list(S_Two, .sq(fcall("sin", list(a))))), fcall("sin", list(.mul(list(S_Two, a))))),
    cot = function(a) .div(fcall("sin", list(.mul(list(S_Two, a)))), .mul(list(S_Two, .sq(fcall("sin", list(a)))))),
    csc = function(a) .inv(fcall("sin", list(a))),
    sec = function(a) .inv(fnode("sin", list(.add(list(a, .mul(list(S_Half, S_Pi))))))),
    exp = function(a) {
      ia <- .mul(list(S_I, a))
      .sub(fcall("sin", list(.add(list(ia, .mul(list(S_Half, S_Pi)))))), .mul(list(S_I, fcall("sin", list(ia)))))
    }
  ),
  tan = list(
    sin = function(a) { t <- fcall("tan", list(.div(a, S_Two))); .div(.mul(list(S_Two, t)), .add(list(S_One, .sq(t)))) },
    cos = function(a) { t <- fcall("tan", list(.div(a, S_Two))); .div(.sub(S_One, .sq(t)), .add(list(S_One, .sq(t)))) },
    cot = function(a) .inv(fcall("tan", list(a)))
  ),
  log = list(
    asin = function(a) .mul(list(.neg(S_I), fcall("log", list(.add(list(.mul(list(S_I, a)), .sqrt(.sub(S_One, .sq(a))))))))),
    acos = function(a) .add(list(.mul(list(S_Half, S_Pi)), .mul(list(S_I, fcall("log", list(.add(list(.mul(list(S_I, a)), .sqrt(.sub(S_One, .sq(a))))))))))),
    atan = function(a) .mul(list(.div(S_I, S_Two), .sub(fcall("log", list(.sub(S_One, .mul(list(S_I, a))))), fcall("log", list(.add(list(S_One, .mul(list(S_I, a))))))))),
    asinh = function(a) fcall("log", list(.add(list(a, .sqrt(.add(list(.sq(a), S_One))))))),
    acosh = function(a) fcall("log", list(.add(list(a, .mul(list(.sqrt(.add(list(a, S_One))), .sqrt(.sub(a, S_One)))))))),
    atanh = function(a) .mul(list(S_Half, .sub(fcall("log", list(.add(list(S_One, a)))), fcall("log", list(.sub(S_One, a)))))),
    acot = function(a) {
      ia <- .div(S_I, a)
      .mul(list(.div(S_I, S_Two), .sub(fcall("log", list(.sub(S_One, ia))), fcall("log", list(.add(list(S_One, ia)))))))
    },
    asec = function(a) {
      l <- fcall("log", list(.add(list(.div(S_I, a), .sqrt(.sub(S_One, .inv(.sq(a))))))))
      .add(list(.mul(list(S_Half, S_Pi)), .mul(list(S_I, l))))
    },
    acsc = function(a) {
      l <- fcall("log", list(.add(list(.div(S_I, a), .sqrt(.sub(S_One, .inv(.sq(a))))))))
      .mul(list(.neg(S_I), l))
    },
    acoth = function(a) .mul(list(S_Half, .sub(fcall("log", list(.add(list(S_One, .inv(a))))), fcall("log", list(.sub(S_One, .inv(a))))))),
    asech = function(a) {
      r <- .inv(a)
      fcall("log", list(.add(list(r, .mul(list(.sqrt(.sub(r, S_One)), .sqrt(.add(list(r, S_One)))))))))
    },
    acsch = function(a) {
      r <- .inv(a)
      fcall("log", list(.add(list(r, .sqrt(.add(list(.sq(r), S_One)))))))
    }
  ),
  gamma = list(
    factorial = function(a) fcall("gamma", list(.add(list(a, S_One)))),
    binomial = function(n, k) .div(fcall("gamma", list(.add(list(n, S_One)))),
                                   .mul(list(fcall("gamma", list(.add(list(k, S_One)))), fcall("gamma", list(.add(list(.sub(n, k), S_One))))))),
    RisingFactorial = function(x, k) .div(fcall("gamma", list(.add(list(x, k)))), fcall("gamma", list(x))),
    FallingFactorial = function(x, k) .div(fcall("gamma", list(.add(list(x, S_One)))), fcall("gamma", list(.add(list(.sub(x, k), S_One))))),
    beta = function(a, b) .div(.mul(list(fcall("gamma", list(a)), fcall("gamma", list(b)))), fcall("gamma", list(.add(list(a, b)))))
  ),
  factorial = list(
    gamma = function(a) fcall("factorial", list(.sub(a, S_One)))
  ),
  Piecewise = list(
    Abs = function(a) Piecewise(list(a, Ge(a, 0)), list(.neg(a), TRUE)),
    Heaviside = function(a, h0 = S_Half) {
      if (is_zero_num(h0)) return(Piecewise(list(S_Zero, Le(a, 0)), list(S_One, TRUE)))
      if (is_one_num(h0)) return(Piecewise(list(S_Zero, Lt(a, 0)), list(S_One, TRUE)))
      Piecewise(list(S_Zero, Lt(a, 0)), list(h0, Eq(a, 0)), list(S_One, TRUE))
    },
    sign = function(a) Piecewise(list(S_One, Gt(a, 0)), list(S_NegativeOne, Lt(a, 0)), list(S_Zero, TRUE)),
    Max = function(a, b) Piecewise(list(a, Ge(a, b)), list(b, TRUE)),
    Min = function(a, b) Piecewise(list(a, Le(a, b)), list(b, TRUE))
  ),
  erf = list(
    erfc = function(a) .sub(S_One, fcall("erf", list(a))),
    erfi = function(a) .mul(list(.neg(S_I), fcall("erf", list(.mul(list(S_I, a))))))
  ),
  sqrt = list(),
  Heaviside = list()
)

#' Rewrite an expression in terms of other functions
#'
#' `rewrite(e, "exp")` expresses trigonometric and hyperbolic functions with
#' exponentials; `"cos"`, `"sin"` and `"tan"` convert between trigonometric
#' functions; `"log"` writes inverse trigonometric functions as logarithms;
#' `"gamma"` writes factorials and binomials with gamma functions;
#' `"factorial"` does the reverse; `"Piecewise"` rewrites `Abs`, `sign`,
#' `Heaviside`, `Max` and `Min`; `"erf"` rewrites `erfc` and `erfi`.
#'
#' @param e An expression.
#' @param target The target family as a string.
#' @return An expression.
#' @examples
#' x <- Symbol("x")
#' rewrite(sin(x), "exp")
#' rewrite(sym_factorial(Symbol("n")), "gamma")
#' @export
rewrite <- function(e, target) {
  if (inherits(e, "symmatrix")) return(.mat_map(e, function(z) rewrite(z, target)))
  e <- as_sym(e)
  rules <- .rewrite_rules[[target]]
  if (is.null(rules)) stop("unknown rewrite target: ", target)
  .map_bottom_up(e, function(z) {
    if (z$type == "Function" && !is.null(rules[[z$name]])) {
      r <- do.call(rules[[z$name]], z$args)
      if (!is.null(r)) return(r)
    }
    if (target == "exp" && z$type == "Pow" && !is_number_type(z$args[[2L]]) && !.has1(z$args[[1L]], S_E) &&
        is_number_type(z$args[[1L]]))
      return(fcall("exp", list(.mul(list(z$args[[2L]], fcall("log", list(z$args[[1L]])))))))
    z
  })
}
