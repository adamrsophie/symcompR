# Permutation generators, after sympy/combinatorics/generators.py.

#' Permutation group elements
#'
#' `symmetric(n)`, `cyclic(n)`, `alternating(n)` and `dihedral(n)` list the
#' elements of the symmetric, cyclic, alternating and dihedral groups of
#' degree `n` as permutations, in SymPy's order (SymPy's generators of
#' `sympy.combinatorics.generators`).
#'
#' @param n The degree.
#' @return A list of `Permutation` objects.
#' @examples
#' cyclic(4)
#' length(alternating(4))
#' dihedral(3)
#' @export
symmetric <- function(n) {
  n <- as.integer(n)
  lapply(permutations(seq_len(n) - 1L), Permutation)
}

#' @rdname symmetric
#' @export
cyclic <- function(n) {
  n <- as.integer(n)
  gen <- seq_len(n) - 1L
  out <- list()
  for (i in seq_len(n)) {
    out[[i]] <- Permutation(gen)
    gen <- c(gen[-1L], gen[1L])
  }
  out
}

#' @rdname symmetric
#' @export
alternating <- function(n) Filter(perm_is_even, symmetric(n))

#' @rdname symmetric
#' @export
dihedral <- function(n) {
  n <- as.integer(n)
  if (n == 1L) return(list(Permutation(c(0L, 1L)), Permutation(c(1L, 0L))))
  if (n == 2L) return(list(Permutation(0:3), Permutation(c(1L, 0L, 3L, 2L)),
                           Permutation(c(2L, 3L, 0L, 1L)), Permutation(3:0)))
  gen <- seq_len(n) - 1L
  out <- list()
  for (i in seq_len(n)) {
    out[[length(out) + 1L]] <- Permutation(gen)
    out[[length(out) + 1L]] <- Permutation(rev(gen))
    gen <- c(gen[-1L], gen[1L])
  }
  out
}
