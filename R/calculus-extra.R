# not_empty_in, after sympy/calculus/util.py.

#' Where an intersection with a finite set is not empty
#'
#' `not_empty_in(s, x)` takes the intersection of a finite set of
#' expressions in `x` with an interval or a union of intervals, and returns
#' the set of real `x` for which the intersection is not empty, as SymPy's
#' `not_empty_in()`. Points where an element is undefined are excluded.
#'
#' @param finset_intersection A `FiniteSet`, or its intersection with an
#'   interval or a union of intervals, or a union of such sets.
#' @param ... One symbol.
#' @return A set.
#' @examples
#' x <- Symbol("x")
#' not_empty_in(Intersection(FiniteSet(x/2), Interval(0, 1)), x)
#' not_empty_in(Intersection(FiniteSet(x^2/(x + 2)), Interval(1, Inf)), x)
#' @export
not_empty_in <- function(finset_intersection, ...) {
  syms <- list(...)
  if (!length(syms)) stop("One or more symbols must be given in syms.", call. = FALSE)
  if (length(syms) > 1L) stop("more than one variable is not handled", call. = FALSE)
  x <- as_sym(syms[[1L]])
  s <- as_sym(finset_intersection)
  if (s$type == "EmptySet") return(S_EmptySet)
  if (s$type == "Union") return(Union(not_empty_in(do.call(Union, s$args[-1L]), x), s$args[[1L]]))
  if (s$type == "FiniteSet") { fin <- s; dom <- Reals }
  else {
    k <- which(vapply(s$args, function(a) a$type == "FiniteSet", TRUE))
    if (s$type != "Intersection" || length(k) != 1L || length(s$args) != 2L)
      stop("A FiniteSet must be given, not ", .str(s), call. = FALSE)
    fin <- s$args[[k]]; dom <- s$args[[3L - k]]
  }
  elm_domain <- function(e, iv) {
    i <- if (iv$type == "SpecialSet" && iv$name == "Reals") list(a = S_NegativeInfinity, b = S_Infinity, lo = TRUE, ro = TRUE) else .iv(iv)
    sing <- solveset(fraction(together(e))[[2L]], x, Reals)
    d1 <- if (i$ro) { if (is_const(i$b, "oo")) Reals else solveset(Lt(e, i$b), x, Reals) } else solveset(Le(e, i$b), x, Reals)
    d2 <- if (i$lo) { if (is_const(i$a, "-oo")) Reals else solveset(Gt(e, i$a), x, Reals) } else solveset(Ge(e, i$a), x, Reals)
    Complement(Intersection(d1, d2), sing)
  }
  ivs <- if (dom$type == "Union") dom$args else list(dom)
  out <- lapply(ivs, function(iv) do.call(Union, lapply(fin$args, elm_domain, iv = iv)))
  do.call(Union, out)
}
