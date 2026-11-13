# Boolean expression constructors, after sympy/logic/boolalg.py.

.as_bool <- function(e) {
  if (is.logical(e)) return(if (isTRUE(e)) S_true else S_false)
  as_sym(e)
}

.negate_rel <- function(e) {
  a <- e$args[[1L]]
  b <- e$args[[2L]]
  switch(e$type,
         Equality = Ne(a, b), Unequality = Eq(a, b),
         StrictLessThan = Ge(a, b), LessThan = Gt(a, b),
         StrictGreaterThan = Le(a, b), GreaterThan = Lt(a, b))
}

.and_or <- function(type, args) {
  absorb <- if (type == "And") .is_false else .is_true
  neutral <- if (type == "And") .is_true else .is_false
  flat <- list()
  for (a in args) {
    a <- .as_bool(a)
    if (a$type == type) flat <- c(flat, a$args) else flat[[length(flat) + 1L]] <- a
  }
  out <- list()
  for (a in flat) {
    if (absorb(a)) return(if (type == "And") S_false else S_true)
    if (neutral(a)) next
    out[[length(out) + 1L]] <- a
  }
  if (!length(out)) return(if (type == "And") S_true else S_false)
  out <- out[!duplicated(vapply(out, .k, ""))]
  keys <- vapply(out, .k, "")
  # And/Or._new_args_filter: only a relational and its negation collapse;
  # p & ~p, p | ~p and Q.real(x) & ~Q.real(x) stay unevaluated as in SymPy
  for (a in out) {
    if (is_relational(a)) {
      na <- .negate_rel(a)
      if (na$key %in% keys) return(if (type == "And") S_false else S_true)
    }
  }
  if (length(out) == 1L) return(out[[1L]])
  new_sym(type, .sort_by_key(out))
}

#' Boolean operators
#'
#' `And()`, `Or()`, `Not()`, `Xor()`, `Nand()`, `Nor()`, `Implies()`,
#' `Equivalent()` and `ITE()` build logical expressions. Symbols serve as
#' boolean variables, and `&`, `|` and `!` work on `sym` objects.
#'
#' @param ... Boolean expressions.
#' @param remove_true For `Xor()`: with `FALSE` a true argument is kept as
#'   an argument instead of negating the result, as SymPy's `remove_true`.
#' @param a,b,c Boolean expressions.
#' @return A boolean `sym` expression.
#' @examples
#' x <- Symbol("x"); y <- Symbol("y")
#' And(x, Or(y, Not(x)))
#' x & !x
#' @export
And <- function(...) .and_or("And", list(...))
#' @rdname And
#' @export
Or <- function(...) .and_or("Or", list(...))
#' @rdname And
#' @export
Not <- function(a) {
  a <- .as_bool(a)
  if (.is_true(a)) return(S_false)
  if (.is_false(a)) return(S_true)
  if (a$type == "Not") return(a$args[[1L]])
  if (is_relational(a)) return(.negate_rel(a))
  new_sym("Not", list(a))
}
#' @rdname And
#' @export
Xor <- function(..., remove_true = TRUE) {
  args <- lapply(list(...), .as_bool)
  flat <- list()
  for (a in args) if (a$type == "Xor") flat <- c(flat, a$args) else flat[[length(flat) + 1L]] <- a
  ntrue <- 0L
  rest <- list()
  for (a in flat) {
    if (.is_true(a)) ntrue <- ntrue + 1L
    else if (!.is_false(a)) rest[[length(rest) + 1L]] <- a
  }
  if (length(rest)) {
    k <- vapply(rest, .k, "")
    tab <- table(k)
    keepk <- names(tab)[tab %% 2L == 1L]
    rest <- rest[match(keepk, k)]
  }
  if (!length(rest)) return(if (ntrue %% 2L == 1L) S_true else S_false)
  # remove_true = FALSE keeps True as an argument, as SymPy's Xor does
  if (!remove_true && ntrue %% 2L == 1L) return(new_sym("Xor", .sort_by_key(c(rest, list(S_true)))))
  base <- if (length(rest) == 1L) rest[[1L]] else new_sym("Xor", .sort_by_key(rest))
  if (ntrue %% 2L == 1L) Not(base) else base
}
#' @rdname And
#' @export
Nand <- function(...) Not(And(...))
#' @rdname And
#' @export
Nor <- function(...) Not(Or(...))
#' @rdname And
#' @export
Implies <- function(a, b) {
  a <- .as_bool(a)
  b <- .as_bool(b)
  if (.is_false(a) || .is_true(b)) return(S_true)
  if (.is_true(a)) return(b)
  if (.is_false(b)) return(Not(a))
  if (identical(a$key, b$key)) return(S_true)
  new_sym("Implies", list(a, b))
}
#' @rdname And
#' @export
Equivalent <- function(...) {
  args <- lapply(list(...), .as_bool)
  args <- args[!duplicated(vapply(args, .k, ""))]
  if (any(vapply(args, .is_true, TRUE)) && any(vapply(args, .is_false, TRUE))) return(S_false)
  if (any(vapply(args, .is_true, TRUE))) return(do.call(And, args))
  if (any(vapply(args, .is_false, TRUE))) return(do.call(And, lapply(args, Not)))
  if (length(args) == 1L) return(S_true)
  new_sym("Equivalent", .sort_by_key(args))
}
#' @rdname And
#' @export
ITE <- function(a, b, c) {
  a <- .as_bool(a)
  b <- .as_bool(b)
  c <- .as_bool(c)
  if (.is_true(a)) return(b)
  if (.is_false(a)) return(c)
  if (identical(b$key, c$key)) return(b)
  if (.is_true(b) && .is_false(c)) return(a)
  if (.is_false(b) && .is_true(c)) return(Not(a))
  new_sym("ITE", list(a, b, c))
}

register_rebuilder("And", function(e, args) .and_or("And", args))
register_rebuilder("Or", function(e, args) .and_or("Or", args))
register_rebuilder("Not", function(e, args) Not(args[[1L]]))
register_rebuilder("Xor", function(e, args) do.call(Xor, args))
register_rebuilder("Implies", function(e, args) Implies(args[[1L]], args[[2L]]))
register_rebuilder("Equivalent", function(e, args) do.call(Equivalent, args))
register_rebuilder("ITE", function(e, args) ITE(args[[1L]], args[[2L]], args[[3L]]))

# SymPy stores And/Or arguments by ordered(): node count, then
# default_sort_key.
# _node_count: a Float counts one half
# SymPy's Heaviside stores its H0 (1/2) as a second argument
.nodes <- function(e) if (e$type == "Float") 0.5 else if (e$type == "Quantity") 3 else
  1 + sum(vapply(e$args, .nodes, 0)) + (e$type == "Function" && e$name == "Heaviside" && length(e$args) == 1L)
.bool_print_order <- function(args) .key_sort(lapply(args, function(a) list(.nodes(a), .sort_key(a))), args)
# an And, Or or Xor built with evaluate=False keeps its argument order
.bool_args <- function(e) if (isTRUE(e$extra$uneval)) e$args else .bool_print_order(e$args)

# SymPy prints the operands of & and | at the bitwise precedence levels, so
# relationals inside get parentheses, and moves a bound against -oo first.
register_str_printer("And", function(e) {
  args <- .bool_args(e)
  for (j in seq_along(args)) {
    a <- args[[j]]
    if (a$type %in% c("StrictLessThan", "LessThan", "StrictGreaterThan", "GreaterThan") &&
        (is_const(a$args[[1L]], "-oo") || is_const(a$args[[2L]], "-oo")))
      args <- c(list(a), args[-j])
  }
  paste(vapply(args, function(a) .paren(a, PREC$BitwiseAnd), ""), collapse = " & ")
})
register_str_printer("Or", function(e)
  paste(vapply(.bool_args(e), function(a) .paren(a, PREC$BitwiseOr), ""), collapse = " | "))
# a function application has the precedence of Func, below that of Not
register_str_printer("Not", function(e) {
  a <- e$args[[1L]]
  paste0("~", if (a$type == "Function" && !identical(a$name, "Mod")) paste0("(", .str(a), ")") else .paren(a, PREC$Not))
})

register_neval("And", function(e, env) as.numeric(all(vapply(e$args, function(a) .neval(a, env) != 0, TRUE))))
register_neval("Or", function(e, env) as.numeric(any(vapply(e$args, function(a) .neval(a, env) != 0, TRUE))))
register_neval("Not", function(e, env) as.numeric(.neval(e$args[[1L]], env) == 0))
for (.t in c("Equality", "Unequality", "StrictLessThan", "LessThan",
             "StrictGreaterThan", "GreaterThan")) {
  register_neval(.t, local({
    tt <- .t
    function(e, env) {
      a <- .neval(e$args[[1L]], env)
      b <- .neval(e$args[[2L]], env)
      as.numeric(switch(tt, Equality = a == b, Unequality = a != b,
                        StrictLessThan = Re(a) < Re(b), LessThan = Re(a) <= Re(b),
                        StrictGreaterThan = Re(a) > Re(b), GreaterThan = Re(a) >= Re(b)))
    }
  }))
}
rm(.t)
