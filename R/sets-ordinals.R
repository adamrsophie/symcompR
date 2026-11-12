# Ordinal arithmetic in Cantor normal form, after sympy/sets/ordinals.py.

#' Ordinals
#'
#' `Ordinal(...)` is an ordinal number in Cantor normal form, the sum of
#' terms `OmegaPower(exp, mult)` (that is `w^exp * mult`) with decreasing
#' exponents, as SymPy's `Ordinal`. `ord0` is the ordinal zero and
#' `ord_omega` is omega (SymPy's `omega`). `+`, `*` and comparisons follow
#' ordinal arithmetic, so `1 + ord_omega == ord_omega`; `ord_omega^a` raises
#' omega to an ordinal or integer. `is_successor_ordinal()`,
#' `is_limit_ordinal()` and `ordinal_degree()` are SymPy's properties.
#'
#' @param ... `OmegaPower` terms, with decreasing exponents.
#' @param exp An ordinal or a nonnegative integer.
#' @param mult A positive integer.
#' @param o An ordinal.
#' @param x An ordinal or omega power.
#' @param e1,e2 Ordinals or integers.
#' @param y,mx,my,cl,reverse Arguments of the `chooseOpsMethod()` generic.
#' @return An ordinal, an omega power, or a logical.
#' @examples
#' w <- ord_omega
#' (w + 3) * (w + 2)
#' 1 + w == w
#' Ordinal(OmegaPower(3, 2), OmegaPower(1, 4), OmegaPower(0, 7))
#' @name ordinals
NULL

.ord_new <- function(terms) structure(list(terms = terms), class = "sym_ordinal")
.ord_convert <- function(v) {
  if (inherits(v, "sym_ordinal")) return(v)
  z <- gmp::as.bigz(if (is_sym(v)) .Z(v) else v)
  if (z == 0) return(ord0)
  .ord_new(list(OmegaPower(0, z)))
}
.ord_eq <- function(a, b) {
  a <- .ord_convert(a); b <- .ord_convert(b)
  length(a$terms) == length(b$terms) &&
    all(vapply(seq_along(a$terms), function(i) .op_eq(a$terms[[i]], b$terms[[i]]), TRUE))
}
.op_eq <- function(p, q) .ord_eq(p$exp, q$exp) && p$mult == q$mult
.op_lt <- function(p, q) if (.ord_eq(p$exp, q$exp)) p$mult < q$mult else .ord_lt(p$exp, q$exp)
.ord_lt <- function(a, b) {
  a <- .ord_convert(a); b <- .ord_convert(b)
  for (i in seq_len(min(length(a$terms), length(b$terms))))
    if (!.op_eq(a$terms[[i]], b$terms[[i]])) return(.op_lt(a$terms[[i]], b$terms[[i]]))
  length(a$terms) < length(b$terms)
}

#' @rdname ordinals
#' @export
OmegaPower <- function(exp, mult) {
  m <- tryCatch(gmp::as.bigz(if (is_sym(mult)) .Z(mult) else mult), error = function(e) NA)
  if (length(m) != 1L || is.na(m) || m <= 0) stop("multiplicity must be a positive integer", call. = FALSE)
  structure(list(exp = .ord_convert(exp), mult = m), class = "sym_omega_power")
}

#' @rdname ordinals
#' @export
Ordinal <- function(...) {
  terms <- list(...)
  pw <- lapply(terms, `[[`, "exp")
  for (i in seq_len(max(0L, length(pw) - 1L))) if (.ord_lt(pw[[i]], pw[[i + 1L]])) stop("powers must be in decreasing order", call. = FALSE)
  .ord_new(terms)
}

#' @rdname ordinals
#' @export
ord0 <- structure(list(terms = list()), class = "sym_ordinal")

#' @rdname ordinals
#' @export
ord_omega <- structure(list(terms = list(structure(list(exp = structure(list(terms = list(structure(
  list(exp = ord0, mult = gmp::as.bigz(1)), class = "sym_omega_power"))), class = "sym_ordinal"),
  mult = gmp::as.bigz(1)), class = "sym_omega_power"))), class = "sym_ordinal")

.ord_lead <- function(o) { if (!length(o$terms)) stop("ordinal zero has no leading term", call. = FALSE); o$terms[[1L]] }
.ord_trail <- function(o) { if (!length(o$terms)) stop("ordinal zero has no trailing term", call. = FALSE); o$terms[[length(o$terms)]] }

#' @rdname ordinals
#' @export
is_successor_ordinal <- function(o) length(o$terms) > 0L && .ord_eq(.ord_trail(o)$exp, ord0)
#' @rdname ordinals
#' @export
is_limit_ordinal <- function(o) length(o$terms) > 0L && !.ord_eq(.ord_trail(o)$exp, ord0)
#' @rdname ordinals
#' @export
ordinal_degree <- function(o) .ord_lead(o)$exp

.ord_add <- function(a, b) {
  a <- .ord_convert(a); b <- .ord_convert(b)
  if (.ord_eq(b, ord0)) return(a)
  at <- a$terms; bt <- b$terms
  be <- ordinal_degree(b)
  r <- length(at)
  while (r >= 1L && .ord_lt(at[[r]]$exp, be)) r <- r - 1L
  terms <- if (r < 1L) bt
           else if (.ord_eq(at[[r]]$exp, be)) c(at[seq_len(r - 1L)], list(OmegaPower(be, at[[r]]$mult + bt[[1L]]$mult)), bt[-1L])
           else c(at[seq_len(r)], bt)
  do.call(Ordinal, terms)
}

.ord_mul <- function(a, b) {
  a <- .ord_convert(a); b <- .ord_convert(b)
  if (.ord_eq(a, ord0) || .ord_eq(b, ord0)) return(ord0)
  ae <- ordinal_degree(a); am <- .ord_lead(a)$mult
  s <- list()
  if (is_limit_ordinal(b)) {
    for (t in b$terms) s[[length(s) + 1L]] <- OmegaPower(.ord_add(ae, t$exp), t$mult)
  } else {
    for (t in b$terms[-length(b$terms)]) s[[length(s) + 1L]] <- OmegaPower(.ord_add(ae, t$exp), t$mult)
    s[[length(s) + 1L]] <- OmegaPower(ae, am * .ord_trail(b)$mult)
    s <- c(s, a$terms[-1L])
  }
  do.call(Ordinal, s)
}

#' @rdname ordinals
#' @export
Ops.sym_ordinal <- function(e1, e2) {
  switch(.Generic,
    `+` = .ord_add(e1, e2),
    `*` = .ord_mul(e1, e2),
    `^` = {
      if (!.ord_eq(e1, ord_omega)) stop("only omega can be raised to a power", call. = FALSE)
      Ordinal(OmegaPower(e2, 1))
    },
    `==` = .ord_eq(e1, e2), `!=` = !.ord_eq(e1, e2),
    `<` = .ord_lt(e1, e2), `<=` = .ord_eq(e1, e2) || .ord_lt(e1, e2),
    `>` = !(.ord_eq(e1, e2) || .ord_lt(e1, e2)), `>=` = !.ord_lt(e1, e2),
    stop("operator ", .Generic, " is not defined for ordinals", call. = FALSE))
}

#' @rdname ordinals
#' @export
chooseOpsMethod.sym_ordinal <- function(x, y, mx, my, cl, reverse) TRUE

#' @rdname ordinals
#' @export
Ops.sym_omega_power <- function(e1, e2) {
  if (!inherits(e1, "sym_omega_power") || !inherits(e2, "sym_omega_power"))
    stop("OmegaPower compares only with OmegaPower", call. = FALSE)
  switch(.Generic,
    `==` = .op_eq(e1, e2), `!=` = !.op_eq(e1, e2),
    `<` = .op_lt(e1, e2), `<=` = .op_eq(e1, e2) || .op_lt(e1, e2),
    `>` = .op_lt(e2, e1), `>=` = .op_eq(e1, e2) || .op_lt(e2, e1),
    stop("operator ", .Generic, " is not defined for OmegaPower", call. = FALSE))
}

#' @rdname ordinals
#' @export
format.sym_ordinal <- function(x, ...) {
  if (!length(x$terms)) return("ord0")
  parts <- vapply(x$terms, function(t) {
    e <- t$exp
    s <- if (.ord_eq(e, ord0)) as.character(t$mult)
         else if (.ord_eq(e, 1)) "w"
         else if (length(e$terms) > 1L || is_limit_ordinal(e)) paste0("w^(", format(e), ")")
         else paste0("w^", format(e))
    if (t$mult != 1 && !.ord_eq(e, ord0)) s <- paste0(s, "*", as.character(t$mult))
    s
  }, "")
  paste(parts, collapse = " + ")
}

#' @rdname ordinals
#' @export
print.sym_ordinal <- function(x, ...) { cat(format(x), "\n", sep = ""); invisible(x) }

#' @rdname ordinals
#' @export
format.sym_omega_power <- function(x, ...) paste0("OmegaPower(", format(x$exp), ", ", as.character(x$mult), ")")

#' @rdname ordinals
#' @export
print.sym_omega_power <- function(x, ...) { cat(format(x), "\n", sep = ""); invisible(x) }
