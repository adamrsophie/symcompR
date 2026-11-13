# Memoized recurrences, after sympy/utilities/memoization.py.

#' Memoized recurrences
#'
#' `recurrence_memo(initial)` returns a decorator: given `f(n, prev)`, which
#' computes term `n` from the list `prev` of the terms before it, it returns
#' a function of `n` that caches every term computed, starting from the
#' terms in `initial`. `assoc_recurrence_memo(base_seq)` does the same for a
#' two index sequence `F(n, m)` with `F(i, 0) = base_seq(i)` and
#' `f(i, j, prev)` giving `F(i, j)` for `1 <= j <= i`, where `prev[[i + 1]]`
#' holds the terms `F(i, 0), ..., F(i, j - 1)` (indices start at 0 as in
#' SymPy, so the list positions are shifted by one). As in SymPy, only
#' `m <= n` is cached. The function made by `recurrence_memo()` carries the
#' attributes `cache_length` and `fetch_item`.
#'
#' @param initial A list of the first terms.
#' @param base_seq A function of `i` giving `F(i, 0)`.
#' @return A decorator: a function taking `f` and returning the memoized
#'   function.
#' @examples
#' fib <- recurrence_memo(list(0, 1))(function(n, prev) prev[[n]] + prev[[n - 1]])
#' fib(30)
#' @export
recurrence_memo <- function(initial) {
  cache <- as.list(initial)
  function(f) {
    g <- function(n) {
      L <- length(cache)
      if (n < L) return(cache[[n + 1L]])
      for (i in L:n) cache[[length(cache) + 1L]] <<- f(i, cache)
      cache[[length(cache)]]
    }
    attr(g, "cache_length") <- function() length(cache)
    attr(g, "fetch_item") <- function(x) cache[[x + 1L]]
    g
  }
}

#' @rdname recurrence_memo
#' @export
assoc_recurrence_memo <- function(base_seq) {
  cache <- list()
  function(f) {
    function(n, m) {
      L <- length(cache)
      if (n < L) return(cache[[n + 1L]][[m + 1L]])
      for (i in L:n) {
        cache[[i + 1L]] <<- list(base_seq(i))
        for (j in seq_len(i)) cache[[i + 1L]][[j + 1L]] <<- f(i, j, cache)
      }
      cache[[n + 1L]][[m + 1L]]
    }
  }
}
