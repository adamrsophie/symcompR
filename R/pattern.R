# Pattern matching with Wild symbols, after Basic.match and Basic.replace.

.wild_ok <- function(w, value) {
  for (ex in w$extra$exclude) if (.has1(value, ex)) return(FALSE)
  for (p in w$extra$properties) if (!isTRUE(p(value))) return(FALSE)
  TRUE
}

# The matcher of sympy/core/basic.py, operations.py, mul.py, add.py,
# power.py and symbol.py. A replacement dictionary is a named list from the
# wild's key to list(w = wild, v = value).
.rd_xrep <- function(p, d) if (length(d)) xreplace(p, unname(lapply(d, function(z) list(z$w, z$v)))) else p
.has_wildp <- function(z) .has_type(z, "Wild")

.matches <- function(p, e, d = list(), old = FALSE) {
  if (is.null(d)) d <- list()
  switch(p$type,
    Wild = .wild_matches(p, e, d),
    Add = .matches_commutative(p, e, d, old),
    Mul = .mul_matches(p, e, d, old),
    Pow = .pow_matches(p, e, d, old),
    if (!length(p$args)) (if (identical(p$key, e$key)) d else NULL) else .basic_matches(p, e, d, old))
}

# Wild.matches (WildFunction carries its own property)
.wild_matches <- function(w, e, d) {
  if (!.wild_ok(w, e)) return(NULL)
  d[[w$key]] <- list(w = w, v = e)
  d
}

# Basic.matches: the same class and the arguments in order
.basic_matches <- function(p, e, d, old) {
  if (!identical(p$type, e$type) || !identical(p$name, e$name)) return(NULL)
  if (identical(p$key, e$key)) return(d)
  if (length(p$args) != length(e$args)) return(NULL)
  for (i in seq_along(p$args)) {
    a <- p$args[[i]]; o <- e$args[[i]]
    if (identical(a$key, o$key)) next
    d <- if (is_relational(a)) tryCatch(.matches(.rd_xrep(a, d), o, d, old), error = function(err) NULL)
         else .matches(.rd_xrep(a, d), o, d, old)
    if (is.null(d)) return(NULL)
  }
  d
}

.op_identity <- function(t) if (t == "Add") S_Zero else S_One
.op_raw <- function(t, args) {
  if (!length(args)) return(.op_identity(t))
  if (length(args) == 1L) return(args[[1L]])
  new_sym(t, args)
}

# Mul._combine_inverse and Add._combine_inverse
.combine_inverse <- function(t, lhs, rhs) {
  if (t == "Mul") {
    if (identical(lhs$key, rhs$key)) return(S_One)
    if (any(vapply(list(lhs, rhs), function(i) i$type %in% c("Pow", "Mul"), TRUE))) {
      dI <- Dummy("I")
      pd <- function(z) { z <- xreplace(z, list(list(S_I, dI))); r <- list(); for (f in .margs(z)) { be <- .abe(f); k <- be[[1L]]$key; r[[k]] <- list(be[[1L]], if (is.null(r[[k]])) be[[2L]] else .add(list(r[[k]][[2L]], be[[2L]]))) }; r }
      A <- pd(lhs); B <- pd(rhs); blen <- length(B)
      for (k in names(B)) if (!is.null(A[[k]])) {
        v <- .sub(A[[k]][[2L]], B[[k]][[2L]]); B[[k]] <- NULL
        if (is_zero_num(v)) A[[k]] <- NULL else A[[k]][[2L]] <- v
      }
      if (length(B) != blen) {
        back <- list(list(dI, S_I))
        lhs <- xreplace(.mul(lapply(A, function(q) .pw(q[[1L]], q[[2L]]))), back)
        rhs <- xreplace(.mul(lapply(B, function(q) .pw(q[[1L]], q[[2L]]))), back)
      }
    }
    rv <- .div(lhs, rhs)
  } else {
    if (.contains_key(lhs, S_Infinity$key) || .contains_key(rhs, S_Infinity$key) ||
        .contains_key(lhs, S_NegativeInfinity$key) || .contains_key(rhs, S_NegativeInfinity$key)) {
      oo <- Dummy("oo")
      reps <- list(list(S_Infinity, oo), list(S_NegativeInfinity, .neg(oo)))
      eq <- .sub(xreplace(lhs, reps), xreplace(rhs, reps))
      if (.has1(eq, oo)) eq <- .map_bottom_up(eq, function(z) if (z$type == "Pow" && identical(z$args[[1L]]$key, oo$key)) z$args[[1L]] else z)
      rv <- xreplace(eq, list(list(oo, S_Infinity)))
    } else rv <- .sub(lhs, rhs)
  }
  srv <- signsimp(rv)
  if (is_number_type(srv)) srv else rv
}

# AssocOp._matches_simple of Mul and Add
.matches_simple <- function(p, e, d) {
  if (p$type == "Mul") {
    sp <- as_coeff_mul_split(p)
    terms <- .margs(sp[[2L]])
    if (length(terms) == 1L) return(.matches(terms[[1L]], .combine_inverse("Mul", e, sp[[1L]]), d))
  } else {
    ca <- .coeff_add_split(p)
    terms <- .aargs(ca[[2L]])
    if (length(terms) == 1L) return(.matches(terms[[1L]], .sub(e, ca[[1L]]), d))
  }
  NULL
}

.add_coeff_key <- function(xs) {
  # sorted(..., key = lambda x: x.args[0] if x.is_Mul and x.args[0].is_Number else 0), stable
  v <- vapply(xs, function(x) if (x$type == "Mul" && is_number_type(x$args[[1L]])) as.numeric(.sort_num(x$args[[1L]])) else 0, 0)
  xs[order(v, method = "radix")]
}

# AssocOp._matches_commutative
.matches_commutative <- function(p, e, d, old) {
  t <- p$type
  if (.is_expr_node(p) && !.is_expr_node(e)) return(NULL)
  if (identical(p$key, e$key)) return(d)
  r <- .matches_simple(p, e, d)
  if (!is.null(r)) return(r)
  isw <- vapply(p$args, function(a) .has_wildp(a) && !.contains_key(e, a$key), TRUE)
  wild_part <- p$args[isw]; exact_part <- p$args[!isw]
  if (!length(exact_part)) {
    wild_part <- ordered(wild_part)
    if (t == "Add") wild_part <- .add_coeff_key(wild_part)
  } else {
    exact <- .op_raw(t, exact_part)
    free <- vapply(free_symbols(e), .k, "")
    if (length(free) && length(setdiff(vapply(free_symbols(exact), .k, ""), free))) return(NULL)
    newexpr <- .combine_inverse(t, e, exact)
    if (!old && e$type %in% c("Add", "Mul")) {
      check <- newexpr
      if (.coeff_isneg(check)) check <- .neg(check)
      if (count_ops(check) > count_ops(e)) return(NULL)
    }
    return(.matches(.op_raw(t, wild_part), newexpr, d, old))
  }
  i <- 0L; saw <- character()
  while (!(e$key %in% saw)) {
    saw <- c(saw, e$key)
    args <- ordered(if (e$type == t) e$args else list(e))
    if (t == "Add" && e$type == "Add") args <- .add_coeff_key(args)
    expr_list <- c(list(.op_identity(t)), args)
    for (last_op in rev(expr_list)) for (w in rev(wild_part)) {
      d1 <- .matches(w, last_op, d)
      if (!is.null(d1)) {
        d2 <- .matches(.rd_xrep(p, d1), e, d1)
        if (!is.null(d2)) return(d2)
      }
    }
    if (i == 0L) {
      if (t == "Mul") {
        if (e$type == "Pow" && e$args[[2L]]$type == "Integer") {
          b <- e$args[[1L]]; n <- e$args[[2L]]
          e <- if (num_sign(n) > 0) new_sym("Mul", list(b, .pow(b, .sub(n, S_One)))) else new_sym("Mul", list(.inv(b), .pow(b, .add(list(n, S_One)))))
          i <- i + 1L
          next
        }
      } else {
        sp <- as_coeff_mul_split(e); c <- sp[[1L]]; ee <- sp[[2L]]
        if (is_number_type(c) && isTRUE(.sort_num(.abs_num(c)) > 1)) {
          e <- if (num_sign(c) > 0) new_sym("Add", list(ee, .mul(list(.sub(c, S_One), ee)))) else new_sym("Add", list(.neg(ee), .mul(list(.add(list(c, S_One)), ee))))
          i <- i + 1L
          next
        }
        was <- e; did <- character()
        for (w in rev(wild_part)) {
          cw <- .mul(Filter(function(f) !.has_wildp(f), .margs(w)))
          fr <- Filter(function(s) !(s$key %in% did), free_symbols(cw))
          if (length(fr)) { did <- c(did, vapply(fr, .k, "")); e <- collect(e, fr) }
        }
        if (!identical(e$key, was$key)) next
      }
    }
    break
  }
  NULL
}

# Mul.matches with the noncommutative matcher
.mul_matches <- function(p, e, d, old) {
  pc <- .is_comm(p); ec <- .is_comm(e)
  if (pc && ec) return(.matches_commutative(p, e, d, old))
  if (pc != ec) return(NULL)
  split <- function(z) { fs <- .margs(z); cm <- vapply(fs, .is_comm, TRUE); list(fs[cm], fs[!cm]) }
  s1 <- split(p); s2 <- split(e)
  c1 <- if (length(s1[[1L]])) s1[[1L]] else list(S_One); c2 <- if (length(s2[[1L]])) s2[[1L]] else list(S_One)
  d <- .matches(.mul(c1), .mul(c2), d, old)
  if (!length(d) && !identical(vapply(c1, .k, ""), vapply(c2, .k, ""))) return(NULL)
  expand_pows <- function(xs) do.call(c, lapply(xs, function(a) if (a$type == "Pow" && a$args[[2L]]$type == "Integer" && num_sign(a$args[[2L]]) > 0) rep(list(a$args[[1L]]), as.integer(a$args[[2L]]$value)) else list(a)))
  r <- .matches_noncomm(expand_pows(s1[[2L]]), expand_pows(s2[[2L]]), d)
  if (!length(r)) NULL else r
}

.matches_noncomm <- function(nodes, targets, d) {
  agenda <- list(); state <- c(0L, 0L); wd <- list()
  nN <- length(nodes); nT <- length(targets)
  while (state[2L] < nT && state[1L] < nN) {
    ni <- state[1L]; ti <- state[2L]
    node <- nodes[[ni + 1L]]
    if (node$type == "Wild") {
      k <- as.character(ni)
      wd[[k]] <- if (is.null(wd[[k]])) c(ti, ti) else c(wd[[k]][1L], ti)
    }
    sm <- .matches_new_states(wd, state, nodes, targets)
    if (!is.null(sm)) {
      agenda <- c(agenda, sm[[1L]])
      if (length(sm[[2L]])) for (k in names(sm[[2L]])) d[[k]] <- sm[[2L]][[k]]
    }
    if (!length(agenda)) return(NULL)
    state <- agenda[[length(agenda)]]; agenda <- agenda[-length(agenda)]
  }
  d
}

.matches_new_states <- function(wd, state, nodes, targets) {
  ni <- state[1L]; ti <- state[2L]
  nN <- length(nodes); nT <- length(targets)
  node <- nodes[[ni + 1L]]; target <- targets[[ti + 1L]]
  if (ti >= nT - 1L && ni < nN - 1L) return(NULL)
  if (node$type == "Wild") {
    be <- wd[[as.character(ni)]]
    terms <- targets[(be[1L] + 1L):(be[2L] + 1L)]
    mult <- if (length(terms) > 1L) .mul(terms) else terms[[1L]]
    m <- .wild_matches(node, mult, list())
    if (length(m)) {
      others <- Filter(function(k) identical(nodes[[as.integer(k) + 1L]]$key, node$key), names(wd))
      for (k in others) {
        ob <- wd[[k]]; cb <- wd[[as.character(ni)]]
        ot <- targets[(ob[1L] + 1L):(ob[2L] + 1L)]; ct <- targets[(cb[1L] + 1L):(cb[2L] + 1L)]
        for (j in seq_len(min(length(ot), length(ct)))) if (!identical(ct[[j]]$key, ot[[j]]$key)) return(NULL)
      }
      ns <- list(c(ni, ti + 1L))
      if (ni < nN - 1L) ns[[2L]] <- c(ni + 1L, ti + 1L)
      return(list(ns, m))
    }
    return(NULL)
  }
  if (ni >= nN - 1L && ti < nT - 1L) return(NULL)
  m <- .matches(node, target, list())
  if (length(m)) return(list(list(c(ni + 1L, ti + 1L)), m))
  if (identical(node$key, target$key)) return(list(list(c(ni + 1L, ti + 1L)), list()))
  NULL
}

# Pow.matches
.pow_matches <- function(p, e, d, old) {
  if (is_one_num(e) && e$type == "Integer") {
    r <- .matches(p$args[[2L]], S_Zero, d)
    if (!is.null(r)) return(r)
  }
  if (!.is_expr_node(e)) return(NULL)
  be <- .abe(e); b <- be[[1L]]; ee <- be[[2L]]
  sbe <- .abe(p); sb <- sbe[[1L]]; se <- sbe[[2L]]
  if (sb$type %in% c("Symbol", "Dummy", "Wild") && se$type == "Integer" && !is_zero_num(e)) {
    if (isTRUE(ask_fact(ee, "rational"))) return(.matches(sb, .pw(b, .div(ee, se)), d))
    return(.matches(sb, .pow(e, .inv(se)), d))
  }
  d2 <- .matches(p$args[[1L]], b, d)
  if (is.null(d2)) return(NULL)
  d2 <- .matches(.rd_xrep(p$args[[2L]], d2), ee, d2)
  if (is.null(d2)) return(.basic_matches(p, e, d, old))
  d2
}

#' Pattern matching
#'
#' `sym_match(e, pattern)` matches an expression against a pattern containing
#' [Wild()] symbols and returns a named list of bindings, or `NULL`. Sums and
#' products match up to commutativity, and a wild may take the identity.
#' `sym_replace(e, pattern, value)` rewrites every subexpression that matches,
#' bottom up; `value` is an expression in the same wild symbols or a function
#' of the named bindings. When `pattern` is an R predicate, `sym_replace()`
#' replaces every node it accepts by `value` (an expression, or an R function
#' of the node), as SymPy's `replace()` with a callable query.
#'
#' @param e An expression.
#' @param pattern A pattern containing wild symbols, a function class made
#'   by [Function()], which matches every application of that function, or
#'   for `sym_replace()` an R predicate of a node.
#' @param value The replacement; for a function class pattern, another
#'   function class or an R function called with the application's arguments.
#' @param old Use the old matching that skips the operation count check.
#' @param exact `sym_replace()` replaces only when no wild matched zero; by
#'   default when the pattern has more than one wild.
#' @return A named list or `NULL`; `sym_replace()` returns an expression.
#' @examples
#' x <- Symbol("x")
#' a <- Wild("a", exclude = list(x)); b <- Wild("b", exclude = list(x))
#' sym_match(3 * x + 2, a * x + b)
#' sym_replace(sin(x)^2 + 1, sin(Wild("u")), cos(Wild("u")))
#' @export
sym_match <- function(e, pattern, old = FALSE) {
  d <- .matches(as_sym(pattern), as_sym(e), list(), old)
  if (is.null(d)) return(NULL)
  out <- list()
  for (z in d) out[[z$w$name]] <- z$v
  out
}

#' @rdname sym_match
#' @export
sym_replace <- function(e, pattern, value, exact = NULL) {
  # a function class: replace every application, value(*expr.args)
  if (inherits(pattern, "sym_function")) {
    if (!is.function(value)) stop(.sympy_condition("TypeError", "given a type, replace() expects another type or a callable"))
    nm <- attr(pattern, "sym_name")
    return(.map_bottom_up(as_sym(e), function(z) {
      if (z$type == "Function" && identical(z$name, nm)) as_sym(do.call(value, z$args)) else z
    }))
  }
  # a predicate query: value(node) replaces each node it accepts
  if (is.function(pattern)) return(.cg_replace(e, pattern, value))
  pattern <- as_sym(pattern)
  if (is.null(exact)) exact <- length(unique(vapply(atoms(pattern, "Wild"), .k, ""))) > 1L
  .map_bottom_up(as_sym(e), function(z) {
    d <- .matches(pattern, z, list())
    if (is.null(d)) return(z)
    if (exact && !all(vapply(d, function(q) !is_zero_num(q$v), TRUE))) return(z)
    if (is.function(value)) {
      b <- list()
      for (q in d) b[[q$w$name]] <- q$v
      return(as_sym(do.call(value, b)))
    }
    if (!length(d)) return(as_sym(value))
    subs(as_sym(value), unname(lapply(d, function(q) list(q$w, q$v))))
  })
}

#' Wild functions
#'
#' `WildFunction(name, nargs)` is a wild symbol that matches any function
#' application, as SymPy's `WildFunction`; with `nargs` only applications
#' with that many arguments match.
#'
#' @param name The name.
#' @param nargs `NULL` (any number of arguments) or allowed argument counts.
#' @return A wild symbol.
#' @examples
#' x <- Symbol("x"); f <- Function("f")
#' F <- WildFunction("F")
#' sym_match(f(x), F)
#' sym_match(x + cos(x), x + F)
#' sym_match(f(x), WildFunction("G", nargs = 2))
#' @export
WildFunction <- function(name, nargs = NULL) {
  n <- if (is.null(nargs)) NULL else as.integer(nargs)
  new_sym("Wild", name = name,
          extra = list(exclude = list(), facts = .symbol_facts(logical()), wildfunction = TRUE,
                       properties = list(function(v) (v$type == "Function" || isTRUE(v$extra$wildfunction)) &&
                                           (is.null(n) || length(v$args) %in% n))),
          key = paste0("WF", nchar(name), ":", name, "[", paste(n, collapse = ","), "]"))
}
