# Unification, after sympy/unify (core.py, usympy.py, rewrite.py).
# Generators become lists of every substitution SymPy would yield, in order.

#' Unification
#'
#' `unify(x, y, s, variables)` lists the substitutions (named lists from the
#' variables to expressions) that make `x` and `y` equal, treating the
#' symbols in `variables` as wildcards; sums, products, unions,
#' intersections and finite sets match associatively, and commutatively
#' where SymPy's do. `rewriterule(source, target, variables, condition,
#' assume)` returns a function giving the rewrites of an expression.
#' `unify_core()` is the generic algorithm on `unify_compound()`,
#' `unify_variable()` and `unify_condvariable()` trees.
#'
#' @param x,y Expressions or trees.
#' @param s A starting substitution.
#' @param variables Symbols that act as wildcards.
#' @param source,target Patterns.
#' @param condition A predicate on the matched values.
#' @param assume A boolean assumption the match must satisfy.
#' @param op,args,arg,valid Parts of a tree.
#' @param fns Named list with `is_associative` and `is_commutative`.
#' @param ... Unused.
#' @return A list of substitutions, or a rewrite function.
#' @examples
#' x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
#' unify(Symbol("a") + Symbol("b"), x + y, variables = list(x, y))
#' @name unify
NULL

#' @rdname unify
#' @export
unify_compound <- function(op, args) structure(list(op = op, args = args), class = "unify_compound")
#' @rdname unify
#' @export
unify_variable <- function(arg) structure(list(arg = arg), class = "unify_variable")
#' @rdname unify
#' @export
unify_condvariable <- function(arg, valid) structure(list(arg = arg, valid = valid), class = "unify_condvariable")

.u_isvar <- function(x) inherits(x, c("unify_variable", "unify_condvariable"))
.u_iscomp <- function(x) inherits(x, "unify_compound")
.u_isargs <- function(x) is.list(x) && !is_sym(x) && !.u_isvar(x) && !.u_iscomp(x) && !inherits(x, "unify_op")
.u_key <- function(x) {
  if (is_sym(x)) return(paste0("S", x$key))
  if (inherits(x, "unify_variable")) return(paste0("V(", .u_key(x$arg), ")"))
  if (inherits(x, "unify_condvariable")) return(paste0("C(", .u_key(x$arg), ")"))
  if (inherits(x, "unify_op")) return(paste0("O", x$type, ":", x$name))
  if (.u_iscomp(x)) return(paste0("[", .u_key(x$op), "|", paste(vapply(x$args, .u_key, ""), collapse = ","), "]"))
  if (is.list(x)) return(paste0("(", paste(vapply(x, .u_key, ""), collapse = ","), ")"))
  paste0("L", paste(format(x), collapse = ","))
}
.u_eq <- function(a, b) identical(.u_key(a), .u_key(b))

.u_occur <- function(var, x) {
  if (.u_eq(var, x)) return(TRUE)
  if (.u_iscomp(x)) return(.u_occur(var, x$args))
  if (.u_isargs(x)) return(any(vapply(x, function(xi) .u_occur(var, xi), NA)))
  FALSE
}

.u_assoc <- function(s, key, val) { s[[.u_key(key)]] <- list(k = key, v = val); s }

.u_unpack <- function(x) if (.u_iscomp(x) && length(x$args) == 1L) x$args[[1L]] else x

.u_allcomb <- function(A, B, ordered) {
  ordered <- if (ordered == "commutative") 11 else NULL
  if (length(A) < length(B)) { sm <- A; bg <- B; bgB <- TRUE } else { sm <- B; bg <- A; bgB <- .u_eq(A, B) }
  part_of <- function(it, part) lapply(part, function(ind) it[ind + 1L])
  lapply(kbins(seq_along(bg) - 1L, length(sm), ordered = ordered), function(part) {
    if (bgB) list(lapply(A, function(a) list(a)), part_of(B, part))
    else list(part_of(A, part), lapply(B, function(b) list(b)))
  })
}

#' @rdname unify
#' @export
unify_core <- function(x, y, s = list(), fns = list()) {
  is_comm <- if (is.null(fns$is_commutative)) function(z) FALSE else fns$is_commutative
  is_assoc <- if (is.null(fns$is_associative)) function(z) FALSE else fns$is_associative
  if (.u_eq(x, y)) return(list(s))
  if (.u_isvar(x)) return(.u_unify_var(x, y, s, fns))
  if (.u_isvar(y)) return(.u_unify_var(y, x, s, fns))
  out <- list()
  if (.u_iscomp(x) && .u_iscomp(y)) {
    for (sop in unify_core(x$op, y$op, s, fns)) {
      if (isTRUE(is_assoc(x)) && isTRUE(is_assoc(y))) {
        if (length(x$args) < length(y$args)) { a <- x; b <- y } else { a <- y; b <- x }
        combs <- .u_allcomb(a$args, b$args, if (isTRUE(is_comm(x)) && isTRUE(is_comm(y))) "commutative" else "associative")
        for (cb in combs) {
          aa <- lapply(cb[[1L]], function(arg) .u_unpack(unify_compound(a$op, arg)))
          bb <- lapply(cb[[2L]], function(arg) .u_unpack(unify_compound(b$op, arg)))
          out <- c(out, unify_core(aa, bb, sop, fns))
        }
      } else if (length(x$args) == length(y$args)) {
        out <- c(out, unify_core(x$args, y$args, sop, fns))
      }
    }
    return(out)
  }
  if (.u_isargs(x) && .u_isargs(y) && length(x) == length(y)) {
    if (!length(x)) return(list(s))
    for (sh in unify_core(x[[1L]], y[[1L]], s, fns)) out <- c(out, unify_core(x[-1L], y[-1L], sh, fns))
  }
  out
}

.u_unify_var <- function(var, x, s, fns) {
  k <- .u_key(var)
  if (!is.null(s[[k]])) return(unify_core(s[[k]]$v, x, s, fns))
  if (.u_occur(var, x)) return(list())
  if (inherits(var, "unify_condvariable")) return(if (isTRUE(var$valid(x))) list(.u_assoc(s, var, x)) else list())
  list(.u_assoc(s, var, x))
}

# ---------------------------------------------------------------------------
# usympy: SymPy trees

.u_op <- function(e) structure(list(type = e$type, name = if (e$type == "Function") e$name else ""), class = "unify_op")
.U_ASSOC <- c("Add", "Mul", "MatAdd", "MatMul", "Union", "Intersection", "FiniteSet", "Min", "Max", "And", "Or")
.U_COMM <- c("Add", "MatAdd", "Union", "Intersection", "FiniteSet")
.u_opname <- function(op) if (op$type == "Function") op$name else op$type

.u_is_associative <- function(x) .u_iscomp(x) && inherits(x$op, "unify_op") && .u_opname(x$op) %in% .U_ASSOC
.u_is_commutative <- function(x) {
  if (!.u_iscomp(x) || !inherits(x$op, "unify_op")) return(FALSE)
  if (.u_opname(x$op) %in% .U_COMM) return(TRUE)
  if (.u_opname(x$op) == "Mul") return(all(vapply(x$args, function(a) isTRUE(.is_comm(.u_construct(a))), NA)))
  FALSE
}

.u_deconstruct <- function(s, variables = list()) {
  if (is_sym(s) && any(vapply(variables, function(v) identical(v$key, s$key), NA))) return(unify_variable(s))
  if (.u_isvar(s)) return(s)
  if (!is_sym(s) || !length(s$args)) return(s)
  unify_compound(.u_op(s), lapply(s$args, .u_deconstruct, variables = variables))
}

.u_construct <- function(t) {
  if (.u_isvar(t)) return(t$arg)
  if (!.u_iscomp(t)) return(t)
  args <- lapply(t$args, .u_construct)
  op <- t$op
  # AssocOp, Pow and FiniteSet are built unevaluated (evaluate=False)
  if (op$type %in% c("Add", "Mul", "Pow", "FiniteSet") || op$name %in% c("Min", "Max") || op$type %in% c("And", "Or"))
    return(if (op$type == "Function") fcall(op$name, args) else new_sym(op$type, args))
  if (op$type == "Function") return(fcall(op$name, args))
  sym_rebuild(new_sym(op$type, args), args)
}

#' @rdname unify
#' @export
unify <- function(x, y, s = list(), variables = list(), ...) {
  dec <- function(z) .u_deconstruct(as_sym(z), variables)
  s0 <- list()
  if (length(s)) for (p in s) s0 <- .u_assoc(s0, dec(p[[1L]]), dec(p[[2L]]))
  fns <- list(is_associative = .u_is_associative, is_commutative = .u_is_commutative)
  lapply(unify_core(dec(x), dec(y), s0, fns), function(d) {
    out <- lapply(d, function(p) .u_construct(p$v))
    names(out) <- vapply(d, function(p) .str(.u_construct(p$k)), "")
    out
  })
}

#' @rdname unify
#' @export
rewriterule <- function(source, target, variables = list(), condition = NULL, assume = NULL) {
  source <- as_sym(source); target <- as_sym(target)
  function(expr, assumptions = TRUE) {
    out <- list()
    for (m in unify(source, expr, list(), variables = variables)) {
      vals <- lapply(variables, function(v) { k <- .str(v); if (!is.null(m[[k]])) m[[k]] else v })
      if (!is.null(condition) && !isTRUE(do.call(condition, vals))) next
      rules <- lapply(variables, function(v) { k <- .str(v); if (!is.null(m[[k]])) list(v, m[[k]]) else NULL })
      rules <- Filter(Negate(is.null), rules)
      if (!is.null(assume) && !isTRUE(ask(xreplace(assume, rules), assumptions))) next
      # strategies.tools.subs: a top down exact replacement whose parents are
      # rebuilt with Basic.__new__, so nothing evaluates
      e2 <- .u_subs_noeval(target, rules)
      # rebuild: construct(deconstruct(e)) evaluates the result
      out[[length(out) + 1L]] <- .u_construct(.u_deconstruct(e2))
    }
    out
  }
}

.u_subs_noeval <- function(e, rules) {
  for (r in rules) if (identical(e$key, r[[1L]]$key)) return(r[[2L]])
  if (!length(e$args)) return(e)
  args <- lapply(e$args, .u_subs_noeval, rules = rules)
  if (identical(vapply(args, `[[`, "", "key"), vapply(e$args, `[[`, "", "key"))) return(e)
  new_sym(e$type, args, name = e$name, value = e$value, extra = e$extra)
}
