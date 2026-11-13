# Algebraic normal form, Xnor, Exclusive, distribution and gate counts,
# after sympy/logic/boolalg.py.

# boolalg.is_literal
.is_literal_sp <- function(e) {
  if (e$type == "Not") return(.is_literal_sp(e$args[[1L]]))
  if (e$type == "BooleanAtom" || !length(e$args)) return(TRUE)
  if (!(e$type %in% .BOOL_TYPES) && all(vapply(e$args, function(a) !length(a$args), TRUE))) return(TRUE)
  FALSE
}

# the class constructor called with remove_true=False
.bool_cls <- function(type) switch(type,
  And = And, Or = Or,
  Xor = function(...) Xor(..., remove_true = FALSE),
  stop(.sympy_condition("TypeError", sprintf("%s does not take remove_true", type))))

#' Xnor and Exclusive
#'
#' `Xnor(...)` is the negation of [Xor()] and `Exclusive(...)` is true when
#' at most one argument is true, as in SymPy, where both evaluate to `Not`
#' and `And` expressions.
#'
#' @param ... Boolean expressions.
#' @return A boolean expression.
#' @examples
#' A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C")
#' Xnor(A, B)
#' Exclusive(A, B, C)
#' @export
Xnor <- function(...) Not(Xor(...))

#' @rdname Xnor
#' @export
Exclusive <- function(...) {
  args <- lapply(list(...), .as_bool)
  out <- list()
  n <- length(args)
  if (n >= 2L) for (i in 1:(n - 1L)) for (j in (i + 1L):n) out[[length(out) + 1L]] <- Not(And(args[[i]], args[[j]]))
  do.call(And, out)
}

# sympy.logic.boolalg._distribute
.distribute_sp <- function(e, outer, inner) {
  ocls <- .bool_cls(outer); icls <- .bool_cls(inner)
  if (e$type == inner) {
    args <- .bool_args(e)
    hit <- which(vapply(args, function(a) a$type == outer, TRUE))
    if (!length(hit)) return(e)
    conj <- args[[hit[1L]]]
    rest <- do.call(icls, args[-hit[1L]])
    return(do.call(ocls, lapply(.bool_args(conj), function(c) .distribute_sp(icls(c, rest), outer, inner))))
  }
  if (e$type == outer) return(do.call(ocls, lapply(.bool_args(e), .distribute_sp, outer = outer, inner = inner)))
  e
}

#' Distribution laws
#'
#' `distribute_and_over_or(expr)` gives a conjunction of disjunctions and
#' `distribute_xor_over_and(expr)` an exclusive or of conjunctions, as the
#' SymPy functions of the same names ([distribute_or_over_and()] gives a
#' disjunction of conjunctions).
#'
#' @param expr A boolean expression.
#' @return A boolean expression.
#' @examples
#' A <- Symbol("A"); B <- Symbol("B"); C <- Symbol("C")
#' distribute_and_over_or(Or(And(A, B), C))
#' distribute_xor_over_and(And(A, Xor(B, C)))
#' @export
distribute_and_over_or <- function(expr) .distribute_sp(.as_bool(expr), "And", "Or")

#' @rdname distribute_and_over_or
#' @export
distribute_xor_over_and <- function(expr) .distribute_sp(.as_bool(expr), "Xor", "And")

#' Algebraic normal form
#'
#' `to_anf(expr, deep)` rewrites a boolean expression as an exclusive or of
#' conjunctions of symbols (with `True` kept as a term), `is_anf(expr)`
#' tests that form, `anf_coeffs(truthvalues)` gives the coefficients of the
#' form from a truth table, and `ANFform(variables, truthvalues)` builds
#' the expression, as in `sympy.logic.boolalg`. `bool_minterm(k,
#' variables)`, `bool_maxterm()` and `bool_monomial()` give the `k`th
#' minterm, maxterm and monomial, where `k` is an integer or a vector of
#' bits.
#'
#' @param expr A boolean expression.
#' @param deep Convert the arguments too.
#' @param truthvalues A vector of 0 and 1 of length a power of two.
#' @param variables A list of symbols.
#' @param k An integer or a vector of bits.
#' @return A boolean expression, a logical value or an integer vector.
#' @examples
#' x <- Symbol("x"); y <- Symbol("y")
#' to_anf(Or(x, y))
#' is_anf(Xor(x, And(x, y)))
#' anf_coeffs(c(1, 1, 1, 0))
#' ANFform(list(x, y), c(1, 1, 1, 0))
#' bool_minterm(3, list(x, y))
#' @export
to_anf <- function(expr, deep = TRUE) {
  expr <- .as_bool(expr)
  if (is_anf(expr)) return(expr)
  .to_anf_m(expr, deep)
}

# BooleanFunction._to_anf
.bf_to_anf <- function(type, args, deep) {
  new <- lapply(args, function(a) {
    if (deep && (!.is_literal_sp(a) || a$type == "Not")) .to_anf_m(a, deep) else a
  })
  do.call(.bool_cls(type), new)
}

# the to_anf() methods
.to_anf_m <- function(e, deep) {
  switch(e$type,
    And = if (deep) distribute_xor_over_and(.bf_to_anf("And", .bool_args(e), deep)) else e,
    Or = {
      args <- .bool_args(e)
      terms <- list()
      for (j in seq_along(args)) for (cmb in utils::combn(length(args), j, simplify = FALSE)) {
        t <- do.call(And, args[cmb])
        terms[[length(terms) + 1L]] <- if (deep) to_anf(t, deep = deep) else t
      }
      do.call(Xor, c(terms, list(remove_true = FALSE)))
    },
    Not = .bf_to_anf("Xor", list(S_true, e$args[[1L]]), deep),
    Xor = .bf_to_anf("Xor", .bool_args(e), deep),
    Implies = .bf_to_anf("Xor", list(S_true, e$args[[1L]], And(e$args[[1L]], e$args[[2L]])), deep),
    Equivalent = {
      a <- do.call(And, e$args)
      b <- distribute_xor_over_and(do.call(And, lapply(e$args, function(x) to_anf(Not(x), deep = FALSE))))
      .bf_to_anf("Xor", list(a, b), deep)
    },
    ITE = do.call(ITE, lapply(e$args, function(a) if (deep && (!.is_literal_sp(a) || a$type == "Not")) .to_anf_m(a, deep) else a)),
    e)
}

#' @rdname to_anf
#' @export
is_anf <- function(expr) {
  expr <- .as_bool(expr)
  if (.is_literal_sp(expr) && expr$type != "Not") return(TRUE)
  if (expr$type == "And") return(all(vapply(expr$args, is_symbol_like, TRUE)))
  if (expr$type == "Xor") {
    for (a in expr$args) {
      if (a$type == "And") {
        if (!all(vapply(a$args, is_symbol_like, TRUE))) return(FALSE)
      } else if (.is_literal_sp(a)) {
        if (a$type == "Not") return(FALSE)
      } else return(FALSE)
    }
    return(TRUE)
  }
  FALSE
}

#' @rdname to_anf
#' @export
anf_coeffs <- function(truthvalues) {
  tv <- as.integer(unlist(truthvalues))
  n <- nchar(.int_to_bin(length(tv))) - 1L
  if (length(tv) != 2^n)
    stop(.sympy_condition("ValueError", sprintf("The number of truth values must be a power of two, got %d", length(tv))))
  coeffs <- lapply(tv, function(v) v)
  for (i in seq_len(n) - 1L) {
    tmp <- list()
    for (j in seq_len(2^(n - i - 1L)) - 1L) {
      a <- coeffs[[2L * j + 1L]]; b <- coeffs[[2L * j + 2L]]
      tmp[[length(tmp) + 1L]] <- c(a, bitwXor(a, b))
    }
    coeffs <- tmp
  }
  coeffs[[1L]]
}

.int_to_bin <- function(n) {
  if (n == 0) return("0")
  b <- character()
  while (n > 0) { b <- c(as.character(n %% 2), b); n <- n %/% 2 }
  paste(b, collapse = "")
}

# sympy.utilities.iterables.ibin
.anf_ibin <- function(k, n) {
  if (k < 0) stop(.sympy_condition("ValueError", "negative numbers are not allowed"))
  b <- as.integer(strsplit(.int_to_bin(k), "")[[1L]])
  if (length(b) < n) b <- c(rep(0L, n - length(b)), b)
  if (length(b) > n) stop(.sympy_condition("ValueError", sprintf("`bits` must be >= %d", length(b))))
  b
}

.anf_term <- function(term, variables) {
  temp <- variables[which(term == 1L)]
  if (!length(temp)) return(S_true)
  do.call(And, temp)
}

#' @rdname to_anf
#' @export
ANFform <- function(variables, truthvalues) {
  variables <- lapply(variables, as_sym)
  nv <- length(variables); n <- length(unlist(truthvalues))
  if (n != 2^nv)
    stop(.sympy_condition("ValueError", sprintf("The number of truth values must be equal to 2^%d, got %d", nv, n)))
  coeffs <- anf_coeffs(truthvalues)
  terms <- list()
  for (i in seq_len(2^nv) - 1L) if (coeffs[i + 1L] == 1L)
    terms[[length(terms) + 1L]] <- .anf_term(.anf_ibin(i, nv), variables)
  do.call(Xor, c(terms, list(remove_true = FALSE)))
}

#' @rdname to_anf
#' @export
bool_minterm <- function(k, variables) {
  variables <- lapply(variables, as_sym)
  k <- if (length(k) == 1L) .anf_ibin(k, length(variables)) else as.integer(k)
  do.call(And, lapply(which(k != 3L), function(n) if (k[n] == 1L) variables[[n]] else Not(variables[[n]])))
}

#' @rdname to_anf
#' @export
bool_maxterm <- function(k, variables) {
  variables <- lapply(variables, as_sym)
  k <- if (length(k) == 1L) .anf_ibin(k, length(variables)) else as.integer(k)
  do.call(Or, lapply(which(k != 3L), function(n) if (k[n] == 0L) variables[[n]] else Not(variables[[n]])))
}

#' @rdname to_anf
#' @export
bool_monomial <- function(k, variables) {
  variables <- lapply(variables, as_sym)
  k <- if (length(k) == 1L) .anf_ibin(k, length(variables)) else as.integer(k)
  .anf_term(k, variables)
}

#' Number of gate inputs
#'
#' `gateinputcount(expr)` counts the inputs of the logic gates (`And`, `Or`,
#' `Xor`, `Not` and `ITE`) in a boolean expression, as SymPy's
#' `gateinputcount()`.
#'
#' @param expr A boolean expression.
#' @return An integer.
#' @examples
#' a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c")
#' gateinputcount(Or(a, And(b, c)))
#' @export
gateinputcount <- function(expr) {
  if (!is_sym(expr)) expr <- as_sym(expr)
  ok <- expr$type %in% c(.BOOL_TYPES, .REL_TYPES, "BooleanAtom") || is_symbol_like(expr)
  if (!ok) stop(.sympy_condition("TypeError", "Expression must be Boolean"))
  if (expr$type %in% c("And", "Or", "Xor", "Not", "ITE"))
    return(length(expr$args) + sum(vapply(expr$args, gateinputcount, 0L)))
  0L
}
