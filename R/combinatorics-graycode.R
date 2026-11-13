# Gray codes, after sympy/combinatorics/graycode.py. A GrayCode keeps the
# mutable state of SymPy's object (current code, cached rank, skip flag) in
# an environment.

# int(s, 2) for a bit string, as an integer while it fits
.gc_bin2int <- function(s) {
  z <- gmp::as.bigz(0)
  for (ch in strsplit(s, "", fixed = TRUE)[[1L]]) z <- 2 * z + (ch == "1")
  .gc_small(z)
}

.gc_small <- function(z) if (abs(z) < 2^31) as.integer(as.character(z)) else z

# f"{n:b}".rjust(width, "0")
.gc_int2bin <- function(n, width = 0L) {
  n <- gmp::as.bigz(n)
  s <- if (n == 0) "0" else as.character(n, b = 2)
  if (nchar(s) < width) s <- paste0(strrep("0", width - nchar(s)), s)
  s
}

#' Gray codes
#'
#' `GrayCode(n, start, rank)` is SymPy's `GrayCode`: the reflected binary code
#' on `n` bits, positioned at the bit string `start` or at the code of a given
#' `rank`. The fields `$n`, `$selections` (`2^n`), `$current` and `$rank` mirror
#' SymPy's. `graycode_next(g, delta)` moves `delta` codes on;
#' `generate_gray(g, start, rank)` returns the codes from the current one to
#' the end, and with `lazy = TRUE` an iterator (a function returning the next
#' code, then `NULL`) during which `graycode_skip(g)` skips the next code, as
#' SymPy's generator and `skip()` do. `graycode_unrank(n, rank)` gives the
#' code of a rank.
#'
#' `bin_to_gray()` and `gray_to_bin()` convert bit strings,
#' `random_bitstring(n)` draws one, `get_subset_from_bitstring(super_set,
#' bitstring)` selects elements and `graycode_subsets(set)` lists the
#' subsets in Gray code order.
#'
#' @param n Number of bits.
#' @param start A starting bit string.
#' @param rank A rank.
#' @param g A `GrayCode`.
#' @param delta Steps to move.
#' @param lazy Return an iterator instead of the codes.
#' @param bin_list A bit string.
#' @param super_set A vector, list or string of elements.
#' @param bitstring A bit string.
#' @param gray_code_set A vector or list of elements.
#' @return A `GrayCode`, a bit string, or a vector or list.
#' @examples
#' g <- GrayCode(3)
#' generate_gray(g)
#' GrayCode(5, start = "10010")$rank
#' graycode_next(GrayCode(6), -1)$current
#' bin_to_gray("111")
#' graycode_subsets(c("a", "b", "c"))
#' @export
GrayCode <- function(n, start = NULL, rank = NULL) {
  if (n < 1 || as.integer(n) != n)
    stop(.sympy_condition("ValueError", sprintf("Gray code dimension must be a positive integer, not %d", as.integer(n))))
  g <- new.env(parent = emptyenv())
  g$n <- as.integer(n)
  g$current_ <- 0L
  g$rank_ <- NULL
  g$skip_ <- FALSE
  if (!is.null(start)) {
    if (length(start) > 1L || !is.character(start)) {
      if (length(start) > n)
        stop(.sympy_condition("ValueError", sprintf("Gray code start has length %d but should not be greater than %d", length(start), as.integer(n))))
      start <- paste(start, collapse = "")
    }
    if (nchar(start) > n)
      stop(.sympy_condition("ValueError", sprintf("Gray code start has length %d but should not be greater than %d", nchar(start), as.integer(n))))
    g$current_ <- start
  } else if (!is.null(rank)) {
    if (as.integer(rank) != rank)
      stop(.sympy_condition("ValueError", sprintf("Gray code rank must be a positive integer, not %d", as.integer(rank))))
    g$rank_ <- .gc_small(gmp::as.bigz(rank) %% (gmp::as.bigz(2)^g$n))
    g$current_ <- graycode_unrank(g$n, g$rank_)
  }
  class(g) <- "GrayCode"
  g
}

.gc_current <- function(g) {
  rv <- g$current_
  if (!is.character(rv)) rv <- if (identical(rv, 0L) || rv == 0) "0" else .gc_int2bin(rv)
  if (nchar(rv) < g$n) rv <- paste0(strrep("0", g$n - nchar(rv)), rv)
  rv
}

#' @export
`$.GrayCode` <- function(x, name) {
  switch(name,
         current = .gc_current(x),
         rank = {
           if (is.null(get("rank_", envir = x))) assign("rank_", .gc_bin2int(gray_to_bin(.gc_current(x))), envir = x)
           get("rank_", envir = x)
         },
         selections = .gc_small(gmp::as.bigz(2)^get("n", envir = x)),
         get(name, envir = x))
}

#' @export
format.GrayCode <- function(x, ...) paste0("GrayCode(", get("n", envir = x), ")")

#' @export
print.GrayCode <- function(x, ...) {
  cat(format(x), "\n", sep = "")
  invisible(x)
}

#' @export
Ops.GrayCode <- function(e1, e2) {
  eq <- inherits(e2, "GrayCode") && get("n", envir = e1) == get("n", envir = e2)
  switch(.Generic, "==" = eq, "!=" = !eq,
         stop(.sympy_condition("TypeError", paste0("unsupported operand for ", .Generic))))
}

#' @rdname GrayCode
#' @export
graycode_next <- function(g, delta = 1) {
  sel <- gmp::as.bigz(2)^g$n
  GrayCode(g$n, rank = .gc_small((gmp::as.bigz(g$rank) + delta) %% sel))
}

#' @rdname GrayCode
#' @export
graycode_skip <- function(g) {
  assign("skip_", TRUE, envir = g)
  invisible(g)
}

#' @rdname GrayCode
#' @export
generate_gray <- function(g, start = NULL, rank = NULL, lazy = FALSE) {
  bits <- g$n
  if (is.null(start) && !is.null(rank)) start <- graycode_unrank(g$n, rank)
  if (!is.null(start)) assign("current_", start, envir = g)
  current <- g$current
  graycode_bin <- gray_to_bin(current)
  if (nchar(graycode_bin) > bits)
    stop(.sympy_condition("ValueError", sprintf("Gray code start has length %d but should not be greater than %d", nchar(graycode_bin), bits)))
  assign("current_", gmp::as.bigz(.gc_bin2int(current)), envir = g)
  i <- gmp::as.bigz(.gc_bin2int(graycode_bin))
  last <- gmp::as.bigz(2)^bits
  nxt <- function() {
    while (i < last) {
      out <- NULL
      if (get("skip_", envir = g)) assign("skip_", FALSE, envir = g) else out <- g$current
      bbtc <- .gc_xor(i, i + 1)
      gbtc <- .gc_xor(bbtc, bbtc %/% 2)
      assign("current_", .gc_xor(get("current_", envir = g), gbtc), envir = g)
      i <<- i + 1
      if (i >= last) assign("current_", 0L, envir = g)
      if (!is.null(out)) return(out)
    }
    NULL
  }
  if (lazy) return(nxt)
  out <- character()
  while (!is.null(code <- nxt())) out <- c(out, code)
  out
}

# bitwise xor of nonnegative big integers
.gc_xor <- function(a, b) {
  a <- gmp::as.bigz(a); b <- gmp::as.bigz(b)
  r <- gmp::as.bigz(0); p <- gmp::as.bigz(1)
  while (a > 0 || b > 0) {
    if ((a %% 2) != (b %% 2)) r <- r + p
    a <- a %/% 2; b <- b %/% 2; p <- 2 * p
  }
  r
}

#' @rdname GrayCode
#' @export
graycode_unrank <- function(n, rank) {
  un <- function(k, n) {
    if (n == 1) return(as.character(k %% 2))
    m <- gmp::as.bigz(2)^(n - 1)
    if (k < m) return(paste0("0", un(k, n - 1)))
    paste0("1", un(m - k %% m - 1, n - 1))
  }
  un(gmp::as.bigz(rank), as.integer(n))
}

#' @rdname GrayCode
#' @export
random_bitstring <- function(n) paste(sample(c("0", "1"), n, replace = TRUE), collapse = "")

#' @rdname GrayCode
#' @export
gray_to_bin <- function(bin_list) {
  s <- strsplit(bin_list, "", fixed = TRUE)[[1L]]
  b <- s[1L]
  for (i in seq_along(s)[-1L]) b[i] <- as.character(as.integer(b[i - 1L] != s[i]))
  paste(b, collapse = "")
}

#' @rdname GrayCode
#' @export
bin_to_gray <- function(bin_list) {
  s <- strsplit(bin_list, "", fixed = TRUE)[[1L]]
  b <- s[1L]
  for (i in seq_along(s)[-1L]) b[i] <- as.character(bitwXor(as.integer(s[i]), as.integer(s[i - 1L])))
  paste(b, collapse = "")
}

#' @rdname GrayCode
#' @export
get_subset_from_bitstring <- function(super_set, bitstring) {
  if (is.character(super_set) && length(super_set) == 1L) super_set <- strsplit(super_set, "", fixed = TRUE)[[1L]]
  bits <- strsplit(bitstring, "", fixed = TRUE)[[1L]]
  if (length(super_set) != length(bits))
    stop(.sympy_condition("ValueError", "The sizes of the lists are not equal"))
  super_set[bits == "1"]
}

#' @rdname GrayCode
#' @export
graycode_subsets <- function(gray_code_set) {
  lapply(generate_gray(GrayCode(length(gray_code_set))), function(b) get_subset_from_bitstring(gray_code_set, b))
}
