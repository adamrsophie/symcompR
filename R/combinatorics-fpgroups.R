# Free groups, finitely presented groups, coset enumeration, rewriting
# systems, group homomorphisms, polycyclic collectors and polyhedra, after
# sympy/combinatorics (free_groups.py, fp_groups.py, coset_table.py,
# rewritingsystem.py, rewritingsystem_fsm.py, homomorphisms.py,
# pc_groups.py and polyhedron.py).
#
# A word of a free group is a list with the group (an environment) and its
# array form split in two integer vectors: `s` holds the generator indices
# (1 based) and `e` the exponents, so SymPy's ((x, 2), (y, -1)) is
# s = c(1, 2), e = c(2, -1). A letter is a signed generator index (x is 1,
# x^-1 is -1). Coset numbers are 0 based, as in SymPy, and `NA` stands for
# None in coset tables.

# ---------------------------------------------------------------------------
# Small helpers

.fg_stop <- function(cls, msg = "") stop(.sympy_condition(cls, msg))

.py_tuple <- function(x) structure(as.list(x), class = "py_tuple")
.py_set <- function(x) structure(as.list(x), class = "py_set")

# Ordered dictionaries keyed by a string, with Python's insertion order:
# updating a key keeps its place, deleting and adding moves it to the end.
.od <- function() structure(list(k = character(), key = list(), v = list()), class = "py_od")
.od_has <- function(d, k) k %in% d$k
.od_get <- function(d, k) {
  i <- match(k, d$k)
  if (is.na(i)) NULL else d$v[[i]]
}
.od_set <- function(d, k, key, value) {
  i <- match(k, d$k)
  if (is.na(i)) {
    i <- length(d$k) + 1L
    d$k[i] <- k
    d$key[i] <- list(key)
  }
  d$v[i] <- list(value)
  d
}
.od_del <- function(d, k) {
  i <- match(k, d$k)
  if (!is.na(i)) {
    d$k <- d$k[-i]
    d$key <- d$key[-i]
    d$v <- d$v[-i]
  }
  d
}

# ---------------------------------------------------------------------------
# Python style strings for the objects of this file

.py_float_str <- function(x) {
  if (is.na(x)) return("nan")
  if (is.infinite(x)) return(if (x > 0) "inf" else "-inf")
  if (x == round(x) && abs(x) < 1e16) return(paste0(format(x, scientific = FALSE), ".0"))
  s <- format(x, digits = 17)
  for (d in 1:17) {
    t <- formatC(x, digits = d, format = "g")
    if (as.numeric(t) == x) {
      s <- t
      break
    }
  }
  s
}

.py_str <- function(x, top = TRUE) {
  if (is.null(x)) return("None")
  if (inherits(x, "FreeGroupElement")) return(.fw_str(x))
  if (inherits(x, "FreeGroup")) return(.fg_str(x))
  if (inherits(x, "FpGroup")) return(.fpg_str(x))
  if (inherits(x, "FpSubgroup")) return(.fps_str(x))
  if (inherits(x, "CosetTable")) return(.ct_str(x))
  if (inherits(x, "Permutation")) return(if (top) .fpp_str(x) else .fpp_repr(x))
  if (inherits(x, "PermutationGroup")) return(.group_str(x))
  if (inherits(x, "Polyhedron")) return(.poly_str(x))
  if (inherits(x, "sym")) return(sstr(x))
  if (inherits(x, "bigz") || inherits(x, "bigq")) {
    if (length(x) == 1L) return(as.character(x))
    return(paste0("[", paste(as.character(x), collapse = ", "), "]"))
  }
  if (is.environment(x)) return(format(x))
  if (inherits(x, "py_tuple")) {
    parts <- vapply(x, .py_str, "", top = FALSE)
    if (length(parts) == 1L) return(paste0("(", parts, ",)"))
    return(paste0("(", paste(parts, collapse = ", "), ")"))
  }
  if (inherits(x, "py_raw")) return(unclass(x))
  if (inherits(x, "py_dictlike")) return(paste0("{", paste(unlist(x), collapse = ", "), "}"))
  if (inherits(x, "py_od")) {
    items <- vapply(seq_along(x$k), function(i) paste0(.py_str(x$key[[i]], top = FALSE), ": ",
                                                       .py_str(x$v[[i]], top = FALSE)), "")
    return(paste0("{", paste(items, collapse = ", "), "}"))
  }
  if (inherits(x, "py_set")) {
    if (!length(x)) return("set()")
    return(paste0("{", paste(vapply(x, .py_str, "", top = FALSE), collapse = ", "), "}"))
  }
  if (is.list(x)) return(paste0("[", paste(vapply(x, .py_str, "", top = FALSE), collapse = ", "), "]"))
  scalar <- length(x) == 1L && !inherits(x, "AsIs")
  one <- function(v) {
    if (is.na(v)) return("None")
    if (is.logical(v)) return(if (v) "True" else "False")
    if (is.integer(v)) return(as.character(v))
    if (is.numeric(v)) return(.py_float_str(v))
    if (is.character(v)) return(if (top && scalar) v else paste0("'", v, "'"))
    as.character(v)
  }
  if (scalar) return(one(x))
  paste0("[", paste(vapply(seq_along(x), function(i) one(x[[i]]), ""), collapse = ", "), "]")
}

#' Python style string of a group theory object
#'
#' `fp_str()` returns the string SymPy's `str()` gives for the objects of
#' the free group, finitely presented group, coset table, homomorphism and
#' polyhedron functions, with `**` shown as `^`. R lists print as Python
#' lists, lists made by these functions as tuples where SymPy returns
#' tuples, `NA` as `None` and `TRUE` as `True`.
#'
#' @param x An object.
#' @return A character scalar.
#' @examples
#' G <- free_group("x, y")
#' x <- G[[2]]; y <- G[[3]]
#' fp_str(x^2 * y^-1)
#' fp_str(word_array_form(x^2 * y^-1))
#' @export
fp_str <- function(x) .py_str(x)

# ---------------------------------------------------------------------------
# Free groups

.fg_cache <- new.env(parent = emptyenv())

.fg_parse_symbols <- function(symbols) {
  if (is.null(symbols) || !length(symbols)) return(character())
  if (inherits(symbols, "sym")) return(sstr(symbols))
  if (is.character(symbols)) {
    if (length(symbols) == 1L && !nzchar(trimws(symbols))) return(character())
    return(unlist(lapply(symbols, .expand_symbol_names)))
  }
  if (is.list(symbols)) {
    if (all(vapply(symbols, is.character, NA))) return(unlist(lapply(symbols, .expand_symbol_names)))
    if (all(vapply(symbols, function(s) inherits(s, "sym"), NA))) return(vapply(symbols, sstr, ""))
  }
  .fg_stop("ValueError", paste0("The type of `symbols` must be one of the following: ",
                                "a str, Symbol/Expr or a sequence of one of these types"))
}

#' Free groups
#'
#' `FreeGroup(symbols)` is the free group on the given symbols; groups on the
#' same symbols are the same object, as in SymPy. `free_group(symbols)`
#' returns a list whose first element is the group and whose remaining
#' elements are its generators. `symbols` is a string such as `"x, y"` or
#' `"x:3"`, a character vector or a list of `Symbol`s. Words multiply with
#' `*`, are raised to integer powers with `^` and divided with `/`; `<`
#' compares them in the length lexicographic order of SymPy. `fp_order()`,
#' `fp_is_abelian()`, `fp_center()`, `fp_elements()` and `fp_contains()`
#' work on free groups too, and `fg_slice()` is SymPy's `F[i:j]`.
#'
#' @param symbols The generator names.
#' @param G A free group.
#' @param from,to Zero based slice limits, as in Python.
#' @param gen A generator of `G`.
#' @return A `FreeGroup` (an environment), a list, or an integer.
#' @examples
#' G <- free_group("x, y, z")
#' F <- G[[1]]; x <- G[[2]]; y <- G[[3]]
#' fp_str(F)
#' fp_str(x^2 * y * x^-1 * x)
#' fp_str(fg_slice(F, 1))
#' @export
FreeGroup <- function(symbols) {
  syms <- .fg_parse_symbols(symbols)
  key <- paste0("FreeGroup|", paste(syms, collapse = "|"), "|", length(syms))
  obj <- .fg_cache[[key]]
  if (!is.null(obj)) return(obj)
  obj <- new.env(parent = emptyenv())
  obj$symbols <- syms
  obj$rank <- length(syms)
  obj$generators <- .py_tuple(lapply(seq_along(syms), function(i) .fw(obj, i, 1L)))
  obj$identity <- .fw(obj, integer(), integer())
  obj$relators <- list()
  class(obj) <- "FreeGroup"
  assign(key, obj, envir = .fg_cache)
  obj
}

#' @rdname FreeGroup
#' @export
free_group <- function(symbols) {
  G <- FreeGroup(symbols)
  out <- c(list(G), as.list(G$generators))
  names(out) <- c("", G$symbols)
  out
}

#' @rdname FreeGroup
#' @export
xfree_group <- function(symbols) {
  G <- FreeGroup(symbols)
  list(G, G$generators)
}

.fg_str <- function(G) {
  if (G$rank > 30L) return(sprintf("<free group with %s generators>", G$rank))
  paste0("<free group on the generators ", .py_str(G$generators), ">")
}

#' @export
print.FreeGroup <- function(x, ...) {
  cat(.fg_str(x), "\n")
  invisible(x)
}

#' @export
format.FreeGroup <- function(x, ...) .fg_str(x)

#' @rdname FreeGroup
#' @export
fg_slice <- function(G, from = 0L, to = G$rank) {
  from <- max(0L, from)
  to <- min(G$rank, to)
  FreeGroup(if (to > from) G$symbols[(from + 1L):to] else character())
}

#' @rdname FreeGroup
#' @export
fg_index <- function(G, gen) {
  if (inherits(gen, "FreeGroupElement") && identical(gen$group, G)) {
    for (i in seq_along(G$generators)) if (.fw_eq(G$generators[[i]], gen)) return(i - 1L)
    .fg_stop("ValueError", "tuple.index(x): x not in tuple")
  }
  .fg_stop("ValueError", sprintf("expected a generator of Free Group %s, got %s", .fg_str(G), .py_str(gen)))
}

.is_fg <- function(G) inherits(G, "FreeGroup")
.is_fpg <- function(G) inherits(G, "FpGroup")
.is_pg <- function(G) inherits(G, "PermutationGroup")
.is_oo <- function(v) inherits(v, "sym") && identical(v$key, S_Infinity$key)

# ---------------------------------------------------------------------------
# Words

.fw <- function(G, s, e) structure(list(group = G, s = as.integer(s), e = as.integer(e)),
                                   class = "FreeGroupElement")

# _reduce_array_form
.fw_reduce <- function(s, e) {
  n <- length(s)
  rs <- integer(n)
  re <- integer(n)
  m <- 0L
  for (i in seq_len(n)) {
    if (e[i] == 0L) next
    if (m > 0L && rs[m] == s[i]) {
      x <- re[m] + e[i]
      if (x == 0L) m <- m - 1L else re[m] <- x
    } else {
      m <- m + 1L
      rs[m] <- s[i]
      re[m] <- e[i]
    }
  }
  list(s = rs[seq_len(m)], e = re[seq_len(m)])
}

.fw_new <- function(G, s, e) {
  r <- .fw_reduce(as.integer(s), as.integer(e))
  .fw(G, r$s, r$e)
}

.fw_letters <- function(w) rep(w$s * as.integer(sign(w$e)), abs(w$e))

.fw_from_letters <- function(G, L) .fw_new(G, abs(L), sign(L))

.fw_letter <- function(w, i) {
  L <- .fw_letters(w)
  if (i < 0L) i <- length(L) + i
  if (i < 0L || i >= length(L)) .fg_stop("IndexError", "tuple index out of range")
  l <- L[i + 1L]
  .fw(w$group, abs(l), sign(l))
}

.fw_len <- function(w) sum(abs(w$e))
.fw_is_id <- function(w) !length(w$s)
.fw_eq <- function(a, b) {
  inherits(b, "FreeGroupElement") && identical(a$group, b$group) &&
    identical(a$s, b$s) && identical(a$e, b$e)
}
.fw_key <- function(w) paste0(w$s, "^", w$e, collapse = " ")

.fw_str <- function(w) {
  if (.fw_is_id(w)) return("<identity>")
  nm <- w$group$symbols[w$s]
  paste0(ifelse(w$e == 1L, nm, paste0(nm, "^", w$e)), collapse = "*")
}

#' @export
format.FreeGroupElement <- function(x, ...) .fw_str(x)

#' @export
print.FreeGroupElement <- function(x, ...) {
  cat(.fw_str(x), "\n")
  invisible(x)
}

.fw_inv <- function(w) .fw(w$group, rev(w$s), -rev(w$e))

.fw_mul <- function(a, b) {
  if (!inherits(b, "FreeGroupElement") || !identical(a$group, b$group))
    .fg_stop("TypeError", "only FreeGroup elements of same FreeGroup can be multiplied")
  if (.fw_is_id(a)) return(b)
  if (.fw_is_id(b)) return(a)
  s <- c(a$s, b$s)
  e <- c(a$e, b$e)
  idx <- length(a$s)
  while (idx >= 1L && idx < length(s) && s[idx] == s[idx + 1L]) {
    e[idx] <- e[idx] + e[idx + 1L]
    s <- s[-(idx + 1L)]
    e <- e[-(idx + 1L)]
    if (e[idx] == 0L) {
      s <- s[-idx]
      e <- e[-idx]
      idx <- idx - 1L
    }
  }
  .fw(a$group, s, e)
}

# _concat_reduced_array_forms over a list of words of one group
.fw_prod <- function(G, words) {
  os <- integer()
  oe <- integer()
  for (w in words) {
    if (!inherits(w, "FreeGroupElement") || !identical(w$group, G))
      .fg_stop("TypeError", "only FreeGroup elements of same FreeGroup can be multiplied")
    ps <- w$s
    pe <- w$e
    if (!length(ps)) next
    n <- length(os)
    if (!n || os[n] != ps[1L]) {
      os <- c(os, ps)
      oe <- c(oe, pe)
      next
    }
    j <- 1L
    while (j <= length(ps) && length(os) && os[length(os)] == ps[j]) {
      x <- oe[length(os)] + pe[j]
      if (x != 0L) {
        oe[length(os)] <- x
        j <- j + 1L
        break
      }
      os <- os[-length(os)]
      oe <- oe[-length(oe)]
      j <- j + 1L
    }
    if (j <= length(ps)) {
      os <- c(os, ps[j:length(ps)])
      oe <- c(oe, pe[j:length(pe)])
    }
  }
  .fw(G, os, oe)
}

.fw_pow <- function(w, n) {
  if (length(n) != 1L || is.na(n) || n != round(n)) .fg_stop("ValueError", paste0(n, " is not an integer"))
  n <- as.integer(n)
  G <- w$group
  if (n == 0L) return(G$identity)
  if (n == 1L) return(w)
  if (.fw_is_id(w)) return(w)
  if (n == -1L) return(.fw_inv(w))
  if (length(w$s) > 1L && !.fw_is_cyc_red(w)) {
    cr <- .fw_cyclic_reduction(w)
    power <- .fw_pow(cr$word, n)
    if (.fw_is_id(cr$removed)) return(power)
    return(.fw_prod(G, list(cr$removed, power, .fw_inv(cr$removed))))
  }
  x <- if (n < 0L) .fw_inv(w) else w
  n <- abs(n)
  s <- x$s
  e <- x$e
  k <- length(s)
  if (k == 1L) return(.fw(G, s, e * n))
  if (s[1L] != s[k]) return(.fw(G, rep(s, n), rep(e, n)))
  mid <- if (k > 2L) 2:(k - 1L) else integer()
  bs <- c(s[mid], s[1L])
  be <- c(e[mid], e[1L] + e[k])
  .fw(G, c(s[1L], rep(bs, n - 1L), s[mid], s[k]), c(e[1L], rep(be, n - 1L), e[mid], e[k]))
}

# lenlex order of SymPy's FreeGroupElement.__lt__
.fw_lt <- function(a, b) {
  if (!inherits(b, "FreeGroupElement") || !identical(a$group, b$group))
    .fg_stop("TypeError", "only FreeGroup elements of same FreeGroup can be compared")
  l <- .fw_len(a)
  m <- .fw_len(b)
  if (l < m) return(TRUE)
  if (l > m) return(FALSE)
  la <- .fw_letters(a)
  lb <- .fw_letters(b)
  for (i in seq_len(l)) {
    p <- abs(la[i])
    q <- abs(lb[i])
    if (p < q) return(TRUE)
    if (p > q) return(FALSE)
    if (sign(la[i]) < sign(lb[i])) return(TRUE)
    if (sign(la[i]) > sign(lb[i])) return(FALSE)
  }
  FALSE
}

#' @export
Ops.FreeGroupElement <- function(e1, e2) {
  w1 <- inherits(e1, "FreeGroupElement")
  w2 <- !missing(e2) && inherits(e2, "FreeGroupElement")
  switch(.Generic,
         "*" = {
           if (!w1) .fg_stop("TypeError", sprintf("unsupported operand type(s) for *: '%s' and 'FreeGroupElement'", .py_type(e1)))
           .fw_mul(e1, e2)
         },
         "^" = {
           if (!w1) .fg_stop("TypeError", "unsupported operand type(s) for ** or pow()")
           .fw_pow(e1, e2)
         },
         "/" = {
           if (!w1 || !w2 || !identical(e1$group, e2$group))
             .fg_stop("TypeError", "only FreeGroup elements of same FreeGroup can be multiplied")
           .fw_mul(e1, .fw_inv(e2))
         },
         "==" = if (w1) .fw_eq(e1, e2) else .fw_eq(e2, e1),
         "!=" = !(if (w1) .fw_eq(e1, e2) else .fw_eq(e2, e1)),
         "<" = .fw_lt(e1, e2),
         "<=" = .fw_eq(e1, e2) || .fw_lt(e1, e2),
         ">" = {
           if (!w2 || !identical(e1$group, e2$group))
             .fg_stop("TypeError", "only FreeGroup elements of same FreeGroup can be compared")
           !(.fw_eq(e1, e2) || .fw_lt(e1, e2))
         },
         ">=" = !.fw_lt(e1, e2),
         "+" = .fg_stop("TypeError", sprintf("unsupported operand type(s) for +: '%s' and '%s'",
                                             .py_type(e1), .py_type(e2))),
         .fg_stop("TypeError", paste0("unsupported operand for FreeGroupElement: ", .Generic)))
}

.py_type <- function(v) {
  if (inherits(v, "FreeGroupElement")) return("FreeGroupElement")
  if (is.integer(v) || (is.numeric(v) && v == round(v))) return("int")
  if (is.numeric(v)) return("float")
  class(v)[1L]
}

.fw_is_cyc_red <- function(w) {
  if (.fw_is_id(w)) return(TRUE)
  L <- .fw_letters(w)
  L[1L] != -L[length(L)]
}

.fw_cyclic_reduction <- function(w) {
  word <- w
  G <- w$group
  g <- G$identity
  while (!.fw_is_cyc_red(word)) {
    k <- length(word$s)
    ex <- min(abs(word$e[1L]), abs(word$e[k]))
    start <- .fw_pow(.fw_letter(word, 0L), ex)
    end <- .fw_pow(.fw_letter(word, -1L), ex)
    word <- .fw_mul(.fw_mul(.fw_inv(start), word), .fw_inv(end))
    g <- .fw_mul(g, start)
  }
  list(word = word, removed = g)
}

.fw_check <- function(w) {
  if (!inherits(w, "FreeGroupElement")) .fg_stop("TypeError", "expected a FreeGroupElement")
  w
}

#' Words in free groups
#'
#' Methods of SymPy's `FreeGroupElement`. A word `w` is made from the
#' generators returned by [free_group()]. `word_array_form()` gives the pairs
#' (generator, exponent), `word_letter_form()` the letters (an inverse letter is
#' a negated symbol), `word_ext_rep()` the flattened array form and
#' `word_letter_form_elm()` the letters as words. `word_len()` is the length,
#' `word_inverse()` the inverse, `word_order()` the order (1 or `oo`) and
#' `word_commutator(w, v)` is `w^-1*v^-1*w*v`. `is_identity_word()` tests
#' for the identity. The other functions are named as in SymPy; positions
#' are zero based and ranges exclude their end, as in Python.
#'
#' @param w,v,word Words of the same free group.
#' @param gen A word (for most functions a generator or its inverse).
#' @param by The replacement word (default the identity).
#' @param words A list of words, or for `eliminate_words()` a two column
#'   specification `list(keys, values)` made with [word_dict()].
#' @param all Replace all occurrences (SymPy's `_all`).
#' @param inverse Also replace occurrences of the inverse of `gen`.
#' @param from_i,to_j Zero based positions.
#' @param strict When `FALSE`, clip the positions to the word.
#' @param start First position to search.
#' @param i A zero based syllable position.
#' @param removed Also return the removed conjugating word.
#' @param G A free group.
#' @return A word, a list, an integer or a logical.
#' @examples
#' G <- free_group("x, y, z")
#' x <- G[[2]]; y <- G[[3]]; z <- G[[4]]
#' w <- x^5 * y * x^2 * y^-4 * x
#' fp_str(eliminate_word(w, x, x^2))
#' fp_str(subword(w, 2, 6))
#' number_syllables(w)
#' exponent_sum(x^2 * y^3, x)
#' is_cyclic_conjugate(x * y, y * x)
#' @export
word_array_form <- function(w) {
  .fw_check(w)
  syms <- lapply(w$group$symbols, Symbol)
  .py_tuple(lapply(seq_along(w$s), function(i) .py_tuple(list(syms[[w$s[i]]], w$e[i]))))
}

#' @rdname word_array_form
#' @export
word_letter_form <- function(w) {
  .fw_check(w)
  syms <- lapply(w$group$symbols, Symbol)
  .py_tuple(lapply(.fw_letters(w), function(l) if (l > 0L) syms[[l]] else .neg(syms[[-l]])))
}

#' @rdname word_array_form
#' @export
word_letter_form_elm <- function(w) {
  .fw_check(w)
  lapply(.fw_letters(w), function(l) .fw(w$group, abs(l), sign(l)))
}

#' @rdname word_array_form
#' @export
word_ext_rep <- function(w) {
  .fw_check(w)
  syms <- lapply(w$group$symbols, Symbol)
  out <- list()
  for (i in seq_along(w$s)) out <- c(out, list(syms[[w$s[i]]], w$e[i]))
  .py_tuple(out)
}

#' @rdname word_array_form
#' @export
word_len <- function(w) .fw_len(.fw_check(w))

#' @rdname word_array_form
#' @export
is_identity_word <- function(w) .fw_is_id(.fw_check(w))

#' @rdname word_array_form
#' @export
word_inverse <- function(w) .fw_inv(.fw_check(w))

#' @rdname word_array_form
#' @export
word_order <- function(w) if (.fw_is_id(.fw_check(w))) S_One else S_Infinity

#' @rdname word_array_form
#' @export
word_commutator <- function(w, v) {
  if (!inherits(v, "FreeGroupElement") || !identical(w$group, v$group))
    .fg_stop("ValueError", "commutator of only FreeGroupElement of the same FreeGroup exists")
  .fw_mul(.fw_mul(.fw_mul(.fw_inv(w), .fw_inv(v)), w), v)
}

#' @rdname word_array_form
#' @export
word_prod <- function(G, words) .fw_prod(G, words)

#' @rdname word_array_form
#' @export
word_dict <- function(words, by) list(keys = words, values = by)

#' @rdname word_array_form
#' @export
eliminate_words <- function(w, words, all = FALSE, inverse = TRUE) {
  again <- TRUE
  new <- w
  is_dict <- !is.null(words$keys) && !inherits(words, "FreeGroupElement")
  keys <- if (is_dict) words$keys else words
  while (again) {
    again <- FALSE
    for (i in seq_along(keys)) {
      prev <- new
      new <- eliminate_word(new, keys[[i]], if (is_dict) words$values[[i]] else NULL,
                            all = all, inverse = inverse)
      if (!.fw_eq(new, prev)) again <- TRUE
    }
  }
  new
}

#' @rdname word_array_form
#' @export
eliminate_word <- function(w, gen, by = NULL, all = FALSE, inverse = TRUE) {
  if (is.null(by)) by <- w$group$identity
  if (.fw_is_id(gen)) .fg_stop("ValueError", "gen must be a non-empty word")
  if (is_independent(w, gen) || .fw_eq(gen, by)) return(w)
  if (.fw_eq(gen, w)) return(by)
  if (.fw_eq(.fw_inv(gen), by)) all <- FALSE
  word <- w
  l <- .fw_len(gen)
  i <- .fw_subword_index(word, gen)
  k <- 1L
  if (is.na(i)) {
    if (!inverse) return(word)
    i <- .fw_subword_index(word, .fw_inv(gen))
    k <- -1L
    if (is.na(i)) return(word)
  }
  word <- .fw_mul(.fw_mul(subword(word, 0L, i), .fw_pow(by, k)),
                  eliminate_word(subword(word, i + l, .fw_len(word)), gen, by))
  if (all) eliminate_word(word, gen, by, all = TRUE, inverse = inverse) else word
}

# zero based index of `word` in `w`, NA when absent
.fw_subword_index <- function(w, word, start = 0L) {
  l <- .fw_len(word)
  a <- .fw_letters(w)
  b <- .fw_letters(word)
  n <- length(a)
  if (start > n - l) return(NA_integer_)
  for (i in start:(n - l)) {
    if (!l || identical(a[(i + 1L):(i + l)], b)) return(i)
  }
  NA_integer_
}

#' @rdname word_array_form
#' @export
subword_index <- function(w, word, start = 0L) {
  i <- .fw_subword_index(w, word, start)
  if (is.na(i)) .fg_stop("ValueError", "The given word is not a subword of self")
  i
}

#' @rdname word_array_form
#' @export
is_dependent <- function(w, word) {
  !is.na(.fw_subword_index(w, word)) || !is.na(.fw_subword_index(w, .fw_inv(word)))
}

#' @rdname word_array_form
#' @export
is_independent <- function(w, word) !is_dependent(w, word)

#' @rdname word_array_form
#' @export
subword <- function(w, from_i, to_j, strict = TRUE) {
  G <- w$group
  if (!strict) {
    from_i <- max(from_i, 0L)
    to_j <- min(.fw_len(w), to_j)
  }
  if (from_i < 0L || to_j > .fw_len(w))
    .fg_stop("ValueError", "`from_i`, `to_j` must be positive and no greater than the length of associative word")
  if (to_j <= from_i) return(G$identity)
  L <- .fw_letters(w)[(from_i + 1L):to_j]
  .fw_from_letters(G, L)
}

#' @rdname word_array_form
#' @export
contains_generators <- function(w) {
  G <- w$group
  out <- list()
  seen <- integer()
  for (s in w$s) if (!(s %in% seen)) {
    seen <- c(seen, s)
    out[[length(out) + 1L]] <- .fw(G, s, 1L)
  }
  .py_set(out)
}

#' @rdname word_array_form
#' @export
cyclic_subword <- function(w, from_i, to_j) {
  G <- w$group
  l <- .fw_len(w)
  if (l == 0L) return(G$identity)
  L <- .fw_letters(w)
  period1 <- trunc(from_i / l)
  if (from_i >= l) {
    from_i <- from_i - l * period1
    to_j <- to_j - l * period1
  }
  diff <- to_j - from_i
  pyslice <- function(v, a, b) {
    n <- length(v)
    if (a < 0) a <- max(0, n + a)
    if (b < 0) b <- max(0, n + b)
    a <- min(a, n)
    b <- min(b, n)
    if (b <= a) v[0] else v[(a + 1L):b]
  }
  word <- pyslice(L, from_i, to_j)
  period2 <- trunc(to_j / l) - 1L
  word <- c(word, rep(L, max(0L, period2)), pyslice(L, 0L, diff - l + from_i - l * period2))
  .fw_from_letters(G, word)
}

#' @rdname word_array_form
#' @export
cyclic_conjugates <- function(w) {
  out <- list(w)
  keys <- .fw_key(w)
  l <- .fw_len(w)
  for (i in seq_len(max(0L, l - 1L))) {
    c <- cyclic_subword(w, i, i + l)
    k <- .fw_key(c)
    if (!(k %in% keys)) {
      keys <- c(keys, k)
      out[[length(out) + 1L]] <- c
    }
  }
  .py_set(out)
}

#' @rdname word_array_form
#' @export
is_cyclic_conjugate <- function(w, v) {
  if (.fw_len(w) != .fw_len(v)) return(FALSE)
  w1 <- identity_cyclic_reduction(w)
  w2 <- identity_cyclic_reduction(v)
  f <- function(u) {
    L <- .fw_letters(u)
    nm <- u$group$symbols[abs(L)]
    paste(ifelse(L > 0L, nm, paste0("-", nm)), collapse = " ")
  }
  str1 <- f(w1)
  str2 <- f(w2)
  if (nchar(str1) != nchar(str2)) return(FALSE)
  grepl(str1, paste(str2, str2), fixed = TRUE)
}

#' @rdname word_array_form
#' @export
number_syllables <- function(w) length(w$s)

.py_idx <- function(n, i) {
  j <- if (i < 0L) n + i else i
  if (j < 0L || j >= n) .fg_stop("IndexError", "tuple index out of range")
  j + 1L
}

#' @rdname word_array_form
#' @export
exponent_syllable <- function(w, i) w$e[.py_idx(length(w$e), i)]

#' @rdname word_array_form
#' @export
generator_syllable <- function(w, i) Symbol(w$group$symbols[w$s[.py_idx(length(w$s), i)]])

#' @rdname word_array_form
#' @export
sub_syllables <- function(w, from_i, to_j) {
  if (!is.numeric(from_i) || !is.numeric(to_j) || from_i != round(from_i) || to_j != round(to_j))
    .fg_stop("ValueError", "both arguments should be integers")
  G <- w$group
  if (to_j <= from_i) return(G$identity)
  n <- length(w$s)
  a <- min(max(from_i, 0L), n)
  b <- min(max(to_j, 0L), n)
  idx <- if (b > a) (a + 1L):b else integer()
  .fw_new(G, w$s[idx], w$e[idx])
}

#' @rdname word_array_form
#' @export
substituted_word <- function(w, from_i, to_j, by) {
  lw <- .fw_len(w)
  if (from_i >= to_j || from_i > lw || to_j > lw) .fg_stop("ValueError", "values should be within bounds")
  if (from_i == 0L && to_j == lw) return(by)
  if (from_i == 0L) return(.fw_mul(by, subword(w, to_j, lw)))
  if (to_j == lw) return(.fw_mul(subword(w, 0L, from_i), by))
  .fw_mul(.fw_mul(subword(w, 0L, from_i), by), subword(w, to_j, lw))
}

#' @rdname word_array_form
#' @export
is_cyclically_reduced <- function(w) .fw_is_cyc_red(w)

#' @rdname word_array_form
#' @export
identity_cyclic_reduction <- function(w) {
  word <- w
  G <- w$group
  while (!.fw_is_cyc_red(word)) {
    k <- length(word$s)
    exp1 <- word$e[1L]
    exp2 <- word$e[k]
    mid <- if (k > 2L) 2:(k - 1L) else integer()
    if (exp1 + exp2 == 0L) word <- .fw_new(G, word$s[mid], word$e[mid])
    else word <- .fw_new(G, c(word$s[1L], word$s[mid]), c(exp1 + exp2, word$e[mid]))
  }
  word
}

#' @rdname word_array_form
#' @export
cyclic_reduction <- function(w, removed = FALSE) {
  r <- .fw_cyclic_reduction(w)
  if (removed) return(.py_tuple(list(r$word, r$removed)))
  r$word
}

#' @rdname word_array_form
#' @export
power_of <- function(w, v) {
  if (.fw_is_id(w)) return(TRUE)
  l <- .fw_len(v)
  if (l == 1L) {
    gens <- contains_generators(w)
    s <- any(vapply(gens, function(g) .fw_eq(g, v) || .fw_eq(g, .fw_inv(v)), NA))
    return(length(gens) == 1L && s)
  }
  r <- .fw_cyclic_reduction(w)
  if (!.fw_is_id(r$removed)) {
    r2 <- .fw_cyclic_reduction(v)
    if (.fw_eq(r$removed, r2$removed)) return(power_of(r$word, r2$word))
    return(FALSE)
  }
  if (.fw_len(w) < l || .fw_len(w) %% l) return(FALSE)
  prefix <- subword(w, 0L, l)
  if (.fw_eq(prefix, v) || .fw_eq(.fw_inv(prefix), v)) return(power_of(subword(w, l, .fw_len(w)), v))
  FALSE
}

#' @rdname word_array_form
#' @export
exponent_sum <- function(w, gen) {
  if (.fw_len(gen) != 1L) .fg_stop("ValueError", "gen must be a generator or inverse of a generator")
  s <- gen$s[1L]
  gen$e[1L] * sum(w$e[w$s == s])
}

#' @rdname word_array_form
#' @export
generator_count <- function(w, gen) {
  if (.fw_len(gen) != 1L || gen$e[1L] < 0L) .fg_stop("ValueError", "gen must be a generator")
  s <- gen$s[1L]
  gen$e[1L] * sum(abs(w$e[w$s == s]))
}

#' @rdname word_array_form
#' @export
word_index <- function(w, gen) {
  if (.fw_len(gen) != 1L) .fg_stop("ValueError", "")
  i <- match(.fw_letters(gen)[1L], .fw_letters(w))
  if (is.na(i)) .fg_stop("ValueError", "tuple.index(x): x not in tuple")
  i - 1L
}

#' @rdname word_array_form
#' @export
word_contains <- function(w, gen) gen$s[1L] %in% w$s

# group level methods shared by free groups and fp groups

#' @rdname FreeGroup
#' @export
fp_contains <- function(G, gen) {
  if (.is_fpg(G)) G <- G$free_group
  inherits(gen, "FreeGroupElement") && identical(gen$group, G)
}

# ---------------------------------------------------------------------------
# Coset tables (coset_table.py)

.ct_opts <- new.env(parent = emptyenv())
.ct_opts$coset_table_max_limit <- 4096000L
.ct_opts$max_stack_size <- 100L

#' Coset table options
#'
#' Sets SymPy's class wide coset table limits: `coset_table_max_limit`, the
#' default largest number of cosets of an enumeration, and `max_stack_size`,
#' the deduction stack size at which the Felsch strategy runs a look ahead.
#'
#' @param coset_table_max_limit,max_stack_size New values, or `NULL` to keep
#'   the current one.
#' @return The previous values, invisibly.
#' @examples
#' old <- coset_table_options(max_stack_size = 10)
#' coset_table_options(max_stack_size = old$max_stack_size)
#' @export
coset_table_options <- function(coset_table_max_limit = NULL, max_stack_size = NULL) {
  old <- list(coset_table_max_limit = .ct_opts$coset_table_max_limit, max_stack_size = .ct_opts$max_stack_size)
  if (!is.null(coset_table_max_limit)) .ct_opts$coset_table_max_limit <- as.integer(coset_table_max_limit)
  if (!is.null(max_stack_size)) .ct_opts$max_stack_size <- as.integer(max_stack_size)
  invisible(old)
}

# column of a letter (1 based): generator i is column 2*i-1, its inverse 2*i
.ct_col <- function(l) ifelse(l > 0L, 2L * l - 1L, -2L * l)

.ct_letter <- function(x) {
  if (inherits(x, "FreeGroupElement")) {
    if (.fw_len(x) != 1L) .fg_stop("KeyError", .fw_str(x))
    return(x$s[1L] * x$e[1L])
  }
  as.integer(x)
}

#' Coset tables and coset enumeration
#'
#' `CosetTable(G, H)` starts a coset table of the subgroup generated by the
#' words `H` in the finitely presented group `G`. `coset_enumeration_r()`
#' runs the relator based (HLT) strategy, `coset_enumeration_c()` the coset
#' table based (Felsch) strategy and `modified_coset_enumeration_r()` the
#' HLT strategy that also records the subgroup words. The table is
#' `C$table`, a list with one integer vector per coset (zero based, `NA` for
#' an undefined entry); `C$p` is the coincidence vector. The `ct_*`
#' functions are the methods of SymPy's `CosetTable` (`ct_scan_and_fill`,
#' `ct_compress`, `ct_standardize`, ...); `x` is a generator or inverse of a
#' generator, cosets are zero based.
#'
#' @param fp_grp An `FpGroup`.
#' @param subgroup,Y A list of words generating the subgroup.
#' @param max_cosets Largest number of cosets.
#' @param draft A coset table to start from.
#' @param incomplete Return a partial table when `max_cosets` is reached.
#' @param C A coset table.
#' @param alpha,beta,gamma,k,lamda,coset Cosets.
#' @param x A generator or inverse generator (a word of length one).
#' @param word,w A word.
#' @param y A word of the subgroup free group (modified methods).
#' @param fill Define new cosets while scanning.
#' @param modified Use the modified methods.
#' @param q The coincidence queue (an integer vector).
#' @param R A list of relators.
#' @param R_c_x,R_c_x_inv Lists of relator conjugates.
#' @return A `CosetTable` (an environment), or what the SymPy method returns.
#' @examples
#' g <- free_group("x, y")
#' x <- g[[2]]; y <- g[[3]]
#' f <- FpGroup(g[[1]], list(x^3, y^3, x^-1 * y^-1 * x * y))
#' C <- coset_enumeration_r(f, list(x))
#' ct_compress(C); ct_standardize(C)
#' fp_str(C$table)
#' @export
CosetTable <- function(fp_grp, subgroup, max_cosets = NULL) {
  if (is.null(max_cosets) || !max_cosets) max_cosets <- .ct_opts$coset_table_max_limit
  C <- new.env(parent = emptyenv())
  C$fp_group <- fp_grp
  C$subgroup <- subgroup
  C$coset_table_limit <- max_cosets
  C$p <- 0L
  k <- length(fp_grp$generators)
  C$A <- as.integer(rbind(seq_len(k), -seq_len(k)))
  C$nA <- 2L * k
  C$table <- list(rep(NA_integer_, 2L * k))
  C$deduction_stack <- list()
  C$grp <- FreeGroup(if (length(subgroup)) paste0("a_", seq_along(subgroup) - 1L) else character())
  C$P <- list(vector("list", 2L * k))
  C$p_p <- list(C$grp$identity)
  class(C) <- "CosetTable"
  C
}

.ct_str <- function(C) {
  sprintf("Coset Table on %s with %s as subgroup generators", .py_str(C$fp_group), .py_str(C$subgroup))
}

#' @export
print.CosetTable <- function(x, ...) {
  cat(.ct_str(x), "\n")
  invisible(x)
}

#' @rdname CosetTable
#' @export
ct_omega <- function(C) which(C$p == seq_along(C$p) - 1L) - 1L

#' @rdname CosetTable
#' @export
ct_n <- function(C) if (!length(C$table)) 0L else max(ct_omega(C)) + 1L

#' @rdname CosetTable
#' @export
ct_A <- function(C) lapply(C$A, function(l) .fw(C$fp_group$identity$group, abs(l), sign(l)))

#' @rdname CosetTable
#' @export
ct_copy <- function(C) {
  # SymPy builds a fresh table on the same group and subgroup (default
  # limit, empty P) and copies the table, p and the deduction stack
  D <- new.env(parent = emptyenv())
  for (nm in c("fp_group", "subgroup", "A", "nA", "grp")) assign(nm, get(nm, envir = C), envir = D)
  D$coset_table_limit <- .ct_opts$coset_table_max_limit
  D$P <- list(vector("list", C$nA))
  D$p_p <- list(C$grp$identity)
  D$table <- C$table
  D$p <- C$p
  D$deduction_stack <- C$deduction_stack
  class(D) <- "CosetTable"
  D
}

#' @rdname CosetTable
#' @export
ct_is_complete <- function(C) {
  for (coset in ct_omega(C)) if (anyNA(C$table[[coset + 1L]])) return(FALSE)
  TRUE
}

.ct_limit_error <- function(C) {
  .fg_stop("ValueError", sprintf(paste0("the coset enumeration has defined more than %s cosets. ",
                                        "Try with a greater value max number of cosets "), C$coset_table_limit))
}

#' @rdname CosetTable
#' @export
ct_define <- function(C, alpha, x, modified = FALSE) {
  alpha <- as.integer(alpha)
  x <- .ct_letter(x)
  len_table <- length(C$table)
  if (len_table >= C$coset_table_limit) .ct_limit_error(C)
  C$table[[len_table + 1L]] <- rep(NA_integer_, C$nA)
  C$P[[len_table + 1L]] <- vector("list", C$nA)
  beta <- len_table
  C$p[len_table + 1L] <- beta
  C$table[[alpha + 1L]][.ct_col(x)] <- beta
  C$table[[beta + 1L]][.ct_col(-x)] <- alpha
  if (modified) {
    C$P[[alpha + 1L]][.ct_col(x)] <- list(C$grp$identity)
    C$P[[beta + 1L]][.ct_col(-x)] <- list(C$grp$identity)
    C$p_p[[beta + 1L]] <- C$grp$identity
  }
  invisible(NULL)
}

#' @rdname CosetTable
#' @export
ct_define_c <- function(C, alpha, x) {
  alpha <- as.integer(alpha)
  x <- .ct_letter(x)
  len_table <- length(C$table)
  if (len_table >= C$coset_table_limit) .ct_limit_error(C)
  C$table[[len_table + 1L]] <- rep(NA_integer_, C$nA)
  beta <- len_table
  C$p[len_table + 1L] <- beta
  C$table[[alpha + 1L]][.ct_col(x)] <- beta
  C$table[[beta + 1L]][.ct_col(-x)] <- alpha
  C$deduction_stack[[length(C$deduction_stack) + 1L]] <- c(alpha, x)
  invisible(NULL)
}

#' @rdname CosetTable
#' @export
ct_scan_c <- function(C, alpha, word) {
  alpha <- as.integer(alpha)
  L <- .fw_letters(word)
  cols <- .ct_col(L)
  icols <- .ct_col(-L)
  f <- alpha
  i <- 1L
  b <- alpha
  j <- length(L)
  while (i <= j && !is.na(C$table[[f + 1L]][cols[i]])) {
    f <- C$table[[f + 1L]][cols[i]]
    i <- i + 1L
  }
  if (i > j) {
    if (f != b) ct_coincidence_c(C, f, b)
    return(invisible(NULL))
  }
  while (j >= i && !is.na(C$table[[b + 1L]][icols[j]])) {
    b <- C$table[[b + 1L]][icols[j]]
    j <- j - 1L
  }
  if (j < i) {
    ct_coincidence_c(C, f, b)
  } else if (j == i) {
    C$table[[f + 1L]][cols[i]] <- b
    C$table[[b + 1L]][icols[i]] <- f
    C$deduction_stack[[length(C$deduction_stack) + 1L]] <- c(f, L[i])
  }
  invisible(NULL)
}

#' @rdname CosetTable
#' @export
ct_coincidence_c <- function(C, alpha, beta) {
  alpha <- as.integer(alpha)
  beta <- as.integer(beta)
  q <- ct_merge(C, alpha, beta, integer())
  while (length(q)) {
    gamma <- q[1L]
    q <- q[-1L]
    for (x in C$A) {
      cx <- .ct_col(x)
      cxi <- .ct_col(-x)
      delta <- C$table[[gamma + 1L]][cx]
      if (!is.na(delta)) {
        C$table[[delta + 1L]][cxi] <- NA_integer_
        C$deduction_stack[[length(C$deduction_stack) + 1L]] <- c(delta, -x)
        mu <- ct_rep(C, gamma)
        nu <- ct_rep(C, delta)
        if (!is.na(C$table[[mu + 1L]][cx])) {
          q <- ct_merge(C, nu, C$table[[mu + 1L]][cx], q)
        } else if (!is.na(C$table[[nu + 1L]][cxi])) {
          q <- ct_merge(C, mu, C$table[[nu + 1L]][cxi], q)
        } else {
          C$table[[mu + 1L]][cx] <- nu
          C$table[[nu + 1L]][cxi] <- mu
        }
      }
    }
  }
  invisible(NULL)
}

#' @rdname CosetTable
#' @export
ct_scan <- function(C, alpha, word, y = NULL, fill = FALSE, modified = FALSE) {
  alpha <- as.integer(alpha)
  L <- .fw_letters(word)
  cols <- .ct_col(L)
  icols <- .ct_col(-L)
  f <- alpha
  i <- 1L
  b <- alpha
  j <- length(L)
  b_p <- y
  if (modified) f_p <- C$grp$identity
  flag <- 0L
  while (fill || flag == 0L) {
    flag <- 1L
    while (i <= j && !is.na(C$table[[f + 1L]][cols[i]])) {
      if (modified) f_p <- .fw_mul(f_p, C$P[[f + 1L]][[cols[i]]])
      f <- C$table[[f + 1L]][cols[i]]
      i <- i + 1L
    }
    if (i > j) {
      if (f != b) {
        if (modified) ct_coincidence(C, f, b, .fw_mul(.fw_inv(f_p), y), modified = TRUE)
        else ct_coincidence(C, f, b)
      }
      return(invisible(NULL))
    }
    while (j >= i && !is.na(C$table[[b + 1L]][icols[j]])) {
      if (modified) b_p <- .fw_mul(b_p, C$P[[b + 1L]][[icols[j]]])
      b <- C$table[[b + 1L]][icols[j]]
      j <- j - 1L
    }
    if (j < i) {
      if (modified) ct_coincidence(C, f, b, .fw_mul(.fw_inv(f_p), b_p), modified = TRUE)
      else ct_coincidence(C, f, b)
    } else if (j == i) {
      C$table[[f + 1L]][cols[i]] <- b
      C$table[[b + 1L]][icols[i]] <- f
      if (modified) {
        C$P[[f + 1L]][[cols[i]]] <- .fw_mul(.fw_inv(f_p), b_p)
        C$P[[b + 1L]][[icols[i]]] <- .fw_mul(.fw_inv(b_p), f_p)
      }
      return(invisible(NULL))
    } else if (fill) {
      ct_define(C, f, L[i], modified = modified)
    }
  }
  invisible(NULL)
}

#' @rdname CosetTable
#' @export
ct_scan_check <- function(C, alpha, word) {
  alpha <- as.integer(alpha)
  L <- .fw_letters(word)
  cols <- .ct_col(L)
  icols <- .ct_col(-L)
  f <- alpha
  i <- 1L
  b <- alpha
  j <- length(L)
  while (i <= j && !is.na(C$table[[f + 1L]][cols[i]])) {
    f <- C$table[[f + 1L]][cols[i]]
    i <- i + 1L
  }
  if (i > j) return(f == b)
  while (j >= i && !is.na(C$table[[b + 1L]][icols[j]])) {
    b <- C$table[[b + 1L]][icols[j]]
    j <- j - 1L
  }
  if (j < i) return(FALSE)
  if (j == i) {
    C$table[[f + 1L]][cols[i]] <- b
    C$table[[b + 1L]][icols[i]] <- f
  }
  TRUE
}

#' @rdname CosetTable
#' @export
ct_merge <- function(C, k, lamda, q, w = NULL, modified = FALSE) {
  k <- as.integer(k)
  lamda <- as.integer(lamda)
  phi <- ct_rep(C, k, modified = modified)
  psi <- ct_rep(C, lamda, modified = modified)
  if (phi != psi) {
    mu <- min(phi, psi)
    v <- max(phi, psi)
    C$p[v + 1L] <- mu
    if (modified) {
      if (v == phi) {
        C$p_p[[phi + 1L]] <- .fw_mul(.fw_mul(.fw_inv(C$p_p[[k + 1L]]), w), C$p_p[[lamda + 1L]])
      } else {
        C$p_p[[psi + 1L]] <- .fw_mul(.fw_mul(.fw_inv(C$p_p[[lamda + 1L]]), .fw_inv(w)), C$p_p[[k + 1L]])
      }
    }
    q <- c(q, v)
  }
  q
}

#' @rdname CosetTable
#' @export
ct_rep <- function(C, k, modified = FALSE) {
  k <- as.integer(k)
  lamda <- k
  rho <- C$p[lamda + 1L]
  if (modified) s <- C$p
  while (rho != lamda) {
    if (modified) s[rho + 1L] <- lamda
    lamda <- rho
    rho <- C$p[lamda + 1L]
  }
  if (modified) {
    rho <- s[lamda + 1L]
    while (rho != k) {
      mu <- rho
      rho <- s[mu + 1L]
      C$p[rho + 1L] <- lamda
      C$p_p[[rho + 1L]] <- .fw_mul(C$p_p[[rho + 1L]], C$p_p[[mu + 1L]])
    }
  } else {
    mu <- k
    rho <- C$p[mu + 1L]
    while (rho != lamda) {
      C$p[mu + 1L] <- lamda
      mu <- rho
      rho <- C$p[mu + 1L]
    }
  }
  lamda
}

#' @rdname CosetTable
#' @export
ct_coincidence <- function(C, alpha, beta, w = NULL, modified = FALSE) {
  alpha <- as.integer(alpha)
  beta <- as.integer(beta)
  q <- ct_merge(C, alpha, beta, integer(), w = w, modified = modified)
  while (length(q)) {
    gamma <- q[1L]
    q <- q[-1L]
    for (x in C$A) {
      cx <- .ct_col(x)
      cxi <- .ct_col(-x)
      delta <- C$table[[gamma + 1L]][cx]
      if (!is.na(delta)) {
        C$table[[delta + 1L]][cxi] <- NA_integer_
        mu <- ct_rep(C, gamma, modified = modified)
        nu <- ct_rep(C, delta, modified = modified)
        if (!is.na(C$table[[mu + 1L]][cx])) {
          if (modified) {
            v <- .fw_mul(.fw_inv(C$p_p[[delta + 1L]]), .fw_inv(C$P[[gamma + 1L]][[cx]]))
            v <- .fw_mul(.fw_mul(v, C$p_p[[gamma + 1L]]), C$P[[mu + 1L]][[cx]])
            q <- ct_merge(C, nu, C$table[[mu + 1L]][cx], q, w = v, modified = TRUE)
          } else {
            q <- ct_merge(C, nu, C$table[[mu + 1L]][cx], q)
          }
        } else if (!is.na(C$table[[nu + 1L]][cxi])) {
          if (modified) {
            v <- .fw_mul(.fw_inv(C$p_p[[gamma + 1L]]), C$P[[gamma + 1L]][[cx]])
            v <- .fw_mul(.fw_mul(v, C$p_p[[delta + 1L]]), C$P[[mu + 1L]][[cxi]])
            q <- ct_merge(C, mu, C$table[[nu + 1L]][cxi], q, w = v, modified = TRUE)
          } else {
            q <- ct_merge(C, mu, C$table[[nu + 1L]][cxi], q)
          }
        } else {
          C$table[[mu + 1L]][cx] <- nu
          C$table[[nu + 1L]][cxi] <- mu
          if (modified) {
            v <- .fw_mul(.fw_mul(.fw_inv(C$p_p[[gamma + 1L]]), C$P[[gamma + 1L]][[cx]]), C$p_p[[delta + 1L]])
            C$P[[mu + 1L]][[cx]] <- v
            C$P[[nu + 1L]][[cxi]] <- .fw_inv(v)
          }
        }
      }
    }
  }
  invisible(NULL)
}

#' @rdname CosetTable
#' @export
ct_scan_and_fill <- function(C, alpha, word) ct_scan(C, alpha, word, fill = TRUE)

#' @rdname CosetTable
#' @export
ct_scan_and_fill_c <- function(C, alpha, word) {
  alpha <- as.integer(alpha)
  L <- .fw_letters(word)
  cols <- .ct_col(L)
  icols <- .ct_col(-L)
  f <- alpha
  i <- 1L
  b <- alpha
  j <- length(L)
  repeat {
    while (i <= j && !is.na(C$table[[f + 1L]][cols[i]])) {
      f <- C$table[[f + 1L]][cols[i]]
      i <- i + 1L
    }
    if (i > j) {
      if (f != b) ct_coincidence_c(C, f, b)
      return(invisible(NULL))
    }
    while (j >= i && !is.na(C$table[[b + 1L]][icols[j]])) {
      b <- C$table[[b + 1L]][icols[j]]
      j <- j - 1L
    }
    if (j < i) {
      ct_coincidence_c(C, f, b)
    } else if (j == i) {
      C$table[[f + 1L]][cols[i]] <- b
      C$table[[b + 1L]][icols[i]] <- f
      C$deduction_stack[[length(C$deduction_stack) + 1L]] <- c(f, L[i])
    } else {
      ct_define_c(C, f, L[i])
    }
  }
}

#' @rdname CosetTable
#' @export
ct_look_ahead <- function(C) {
  R <- C$fp_group$relators
  for (beta in ct_omega(C)) {
    for (w in R) {
      ct_scan(C, beta, w)
      if (C$p[beta + 1L] < beta) break
    }
  }
  invisible(NULL)
}

#' @rdname CosetTable
#' @export
ct_process_deductions <- function(C, R_c_x, R_c_x_inv) {
  while (length(C$deduction_stack)) {
    if (length(C$deduction_stack) >= .ct_opts$max_stack_size) {
      ct_look_ahead(C)
      C$deduction_stack <- list()
      next
    } else {
      d <- C$deduction_stack[[length(C$deduction_stack)]]
      C$deduction_stack[[length(C$deduction_stack)]] <- NULL
      alpha <- d[1L]
      x <- d[2L]
      if (C$p[alpha + 1L] == alpha) {
        for (w in R_c_x) {
          ct_scan_c(C, alpha, w)
          if (C$p[alpha + 1L] < alpha) break
        }
      }
    }
    beta <- C$table[[alpha + 1L]][.ct_col(x)]
    if (!is.na(beta) && C$p[beta + 1L] == beta) {
      for (w in R_c_x_inv) {
        ct_scan_c(C, beta, w)
        if (C$p[beta + 1L] < beta) break
      }
    }
  }
  invisible(NULL)
}

#' @rdname CosetTable
#' @export
ct_process_deductions_check <- function(C, R_c_x, R_c_x_inv) {
  while (length(C$deduction_stack)) {
    d <- C$deduction_stack[[length(C$deduction_stack)]]
    C$deduction_stack[[length(C$deduction_stack)]] <- NULL
    alpha <- d[1L]
    x <- d[2L]
    for (w in R_c_x) if (!ct_scan_check(C, alpha, w)) return(FALSE)
    beta <- C$table[[alpha + 1L]][.ct_col(x)]
    if (!is.na(beta)) {
      for (w in R_c_x_inv) if (!ct_scan_check(C, beta, w)) return(FALSE)
    }
  }
  TRUE
}

#' @rdname CosetTable
#' @export
ct_switch <- function(C, beta, gamma) {
  beta <- as.integer(beta)
  gamma <- as.integer(gamma)
  for (x in C$A) {
    cx <- .ct_col(x)
    z <- C$table[[gamma + 1L]][cx]
    C$table[[gamma + 1L]][cx] <- C$table[[beta + 1L]][cx]
    C$table[[beta + 1L]][cx] <- z
    for (alpha in seq_along(C$p) - 1L) {
      if (C$p[alpha + 1L] == alpha) {
        t <- C$table[[alpha + 1L]][cx]
        if (!is.na(t) && t == beta) C$table[[alpha + 1L]][cx] <- gamma
        else if (!is.na(t) && t == gamma) C$table[[alpha + 1L]][cx] <- beta
      }
    }
  }
  invisible(NULL)
}

#' @rdname CosetTable
#' @export
ct_standardize <- function(C) {
  gamma <- 1L
  n <- ct_n(C)
  for (alpha in seq_len(n) - 1L) {
    for (x in C$A) {
      beta <- C$table[[alpha + 1L]][.ct_col(x)]
      if (is.na(beta)) .fg_stop("TypeError", "'>=' not supported between instances of 'NoneType' and 'int'")
      if (beta >= gamma) {
        if (beta > gamma) ct_switch(C, gamma, beta)
        gamma <- gamma + 1L
        if (gamma == ct_n(C)) return(invisible(NULL))
      }
    }
  }
  invisible(NULL)
}

#' @rdname CosetTable
#' @export
ct_compress <- function(C) {
  gamma <- -1L
  chi <- which(C$p != seq_along(C$p) - 1L) - 1L
  for (alpha in ct_omega(C)) {
    gamma <- gamma + 1L
    # SymPy's loop also evaluates table[beta][A_dict_inv[x]] == gamma, a
    # comparison with no effect (sympy issue 27633)
    if (gamma != alpha) C$table[[gamma + 1L]] <- C$table[[alpha + 1L]]
  }
  C$p <- seq_len(gamma + 1L) - 1L
  C$table <- C$table[seq_len(gamma + 1L)]
  C$table <- lapply(C$table, function(row) row - findInterval(row, chi, left.open = TRUE))
  invisible(NULL)
}

#' @rdname CosetTable
#' @export
ct_conjugates <- function(C, R) {
  R_set <- list()
  keys <- character()
  for (rel in R) {
    for (cc in c(cyclic_conjugates(rel), cyclic_conjugates(.fw_inv(rel)))) {
      k <- .fw_key(cc)
      if (!(k %in% keys)) {
        keys <- c(keys, k)
        R_set[[length(R_set) + 1L]] <- cc
      }
    }
  }
  .ct_split_by_first(C, R_set)
}

.ct_split_by_first <- function(C, R_set) {
  if (any(vapply(R_set, .fw_is_id, NA))) .fg_stop("IndexError", "tuple index out of range")
  first <- vapply(R_set, function(w) .fw_letters(w)[1L], 0L)
  lapply(C$A, function(x) .py_set(R_set[first == x]))
}

#' @rdname CosetTable
#' @export
ct_coset_representative <- function(C, coset) {
  coset <- as.integer(coset)
  for (x in C$A) {
    gamma <- C$table[[coset + 1L]][.ct_col(x)]
    if (coset == 0L) return(C$fp_group$identity)
    if (!is.na(gamma) && gamma < coset) {
      return(.fw_mul(ct_coset_representative(C, gamma), .fw(C$fp_group$identity$group, abs(x), -sign(x))))
    }
  }
  NULL
}

#' @rdname CosetTable
#' @export
ct_modified_define <- function(C, alpha, x) ct_define(C, alpha, x, modified = TRUE)

#' @rdname CosetTable
#' @export
ct_modified_scan <- function(C, alpha, w, y, fill = FALSE) ct_scan(C, alpha, w, y = y, fill = fill, modified = TRUE)

#' @rdname CosetTable
#' @export
ct_modified_scan_and_fill <- function(C, alpha, w, y) ct_modified_scan(C, alpha, w, y, fill = TRUE)

#' @rdname CosetTable
#' @export
ct_modified_merge <- function(C, k, lamda, w, q) ct_merge(C, k, lamda, q, w = w, modified = TRUE)

#' @rdname CosetTable
#' @export
ct_modified_rep <- function(C, k) invisible(ct_rep(C, k, modified = TRUE))

#' @rdname CosetTable
#' @export
ct_modified_coincidence <- function(C, alpha, beta, w) ct_coincidence(C, alpha, beta, w = w, modified = TRUE)

#' @rdname CosetTable
#' @export
coset_enumeration_r <- function(fp_grp, Y, max_cosets = NULL, draft = NULL, incomplete = FALSE,
                                modified = FALSE) {
  C <- CosetTable(fp_grp, Y, max_cosets = max_cosets)
  if (!is.null(draft)) {
    C$table <- draft$table
    C$p <- draft$p
  }
  R <- fp_grp$relators
  for (i in seq_along(Y)) {
    if (modified) ct_modified_scan_and_fill(C, 0L, Y[[i]], C$grp$generators[[i]])
    else ct_scan_and_fill(C, 0L, Y[[i]])
  }
  alpha <- 0L
  while (alpha < ct_n(C)) {
    if (C$p[alpha + 1L] == alpha) {
      res <- tryCatch({
        for (w in R) {
          if (modified) ct_modified_scan_and_fill(C, alpha, w, C$grp$identity)
          else ct_scan_and_fill(C, alpha, w)
          if (C$p[alpha + 1L] < alpha) break
        }
        if (C$p[alpha + 1L] == alpha) {
          for (x in C$A) {
            if (is.na(C$table[[alpha + 1L]][.ct_col(x)])) ct_define(C, alpha, x, modified = modified)
          }
        }
        NULL
      }, ValueError = function(e) if (incomplete) C else stop(e))
      if (!is.null(res)) return(res)
    }
    alpha <- alpha + 1L
  }
  C
}

#' @rdname CosetTable
#' @export
modified_coset_enumeration_r <- function(fp_grp, Y, max_cosets = NULL, draft = NULL, incomplete = FALSE) {
  coset_enumeration_r(fp_grp, Y, max_cosets = max_cosets, draft = draft, incomplete = incomplete,
                      modified = TRUE)
}

#' @rdname CosetTable
#' @export
coset_enumeration_c <- function(fp_grp, Y, max_cosets = NULL, draft = NULL, incomplete = FALSE) {
  X <- fp_grp$generators
  R <- fp_grp$relators
  C <- CosetTable(fp_grp, Y, max_cosets = max_cosets)
  if (!is.null(draft)) {
    C$table <- draft$table
    C$p <- draft$p
    C$deduction_stack <- draft$deduction_stack
    for (alpha in seq_along(C$table) - 1L) for (xi in seq_along(X)) {
      if (!is.na(C$table[[alpha + 1L]][.ct_col(xi)]))
        C$deduction_stack[[length(C$deduction_stack) + 1L]] <- c(alpha, xi)
    }
  }
  R_cyc_red <- lapply(R, identity_cyclic_reduction)
  R_c_list <- ct_conjugates(C, R_cyc_red)
  for (w in Y) ct_scan_and_fill_c(C, 0L, w)
  for (x in C$A) ct_process_deductions(C, R_c_list[[.ct_col(x)]], R_c_list[[.ct_col(-x)]])
  alpha <- 0L
  while (alpha < length(C$table)) {
    if (C$p[alpha + 1L] == alpha) {
      res <- tryCatch({
        for (x in C$A) {
          if (C$p[alpha + 1L] != alpha) break
          if (is.na(C$table[[alpha + 1L]][.ct_col(x)])) {
            ct_define_c(C, alpha, x)
            ct_process_deductions(C, R_c_list[[.ct_col(x)]], R_c_list[[.ct_col(-x)]])
          }
        }
        NULL
      }, ValueError = function(e) if (incomplete) C else stop(e))
      if (!is.null(res)) return(res)
    }
    alpha <- alpha + 1L
  }
  C
}

# ---------------------------------------------------------------------------
# Finitely presented groups (fp_groups.py)

# stable sort of words in SymPy's length lexicographic order
.fw_sort_key <- function(w) {
  L <- .fw_letters(w)
  paste0(sprintf("%08d", length(L)), paste0(sprintf("%06d", 2L * abs(L) - (L < 0L)), collapse = ""))
}
.fw_sort <- function(words, decreasing = FALSE) {
  if (length(words) < 2L) return(unclass(as.list(words)))
  keys <- vapply(words, .fw_sort_key, "")
  unclass(as.list(words))[order(keys, decreasing = decreasing, method = "radix")]
}
# sorted() of lists of words (Python compares lists element by element)
.fw_sort_list <- function(lists) {
  if (length(lists) < 2L) return(lists)
  keys <- vapply(lists, function(l) paste(vapply(l, .fw_sort_key, ""), collapse = "|~"), "")
  lists[order(keys, method = "radix")]
}
# Python's list(set(words)): unique words, first occurrence order
.fw_unique <- function(words) {
  if (!length(words)) return(list())
  as.list(words)[!duplicated(vapply(words, .fw_key, ""))]
}
.fw_in <- function(w, words) {
  if (!length(words)) return(FALSE)
  .fw_key(w) %in% vapply(words, .fw_key, "")
}
.fw_keys <- function(words) vapply(words, .fw_key, "")

# Python int / SymPy Integer arithmetic for group orders
.py_mul <- function(a, b) {
  if (inherits(a, "sym") || inherits(b, "sym")) return(as_sym(a) * as_sym(b))
  a * b
}
.py_div <- function(a, b) {
  if (inherits(a, "sym") || inherits(b, "sym")) return(as_sym(a) / as_sym(b))
  as.double(a) / as.double(b)
}

#' Finitely presented groups
#'
#' `FpGroup(F, relators)` is the group with the generators of the free group
#' `F` and the given relators; `fp_group()` also returns its generators. The
#' `fp_*` functions are the methods of SymPy's `FpGroup`:
#' `fp_order()` (by coset enumeration), `fp_index()`, `fp_coset_enumeration()`,
#' `fp_coset_table()`, `fp_subgroup()` (a presentation of a subgroup by the
#' Reidemeister Schreier algorithm, optionally with the homomorphism into the
#' group), `fp_quotient()`, `fp_reduce()`, `fp_equals()` and
#' `fp_make_confluent()` (rewriting systems), `fp_most_frequent_generator()`
#' and `fp_random()`.
#'
#' @param fr_grp A free group.
#' @param relators A list of words.
#' @param G An `FpGroup`.
#' @param H,gens A list of words generating a subgroup.
#' @param N A list of words, or a normal `FpSubgroup`.
#' @param strategy `"relator_based"` (HLT) or `"coset_table_based"` (Felsch).
#' @param max_cosets,draft,incomplete Passed to the coset enumeration.
#' @param C A coset table to use.
#' @param homomorphism Also return the homomorphism from the subgroup.
#' @param word,word1,word2 Words.
#' @return An `FpGroup` (an environment), a list, a word, an integer, `oo`
#'   or a logical.
#' @examples
#' g <- free_group("x, y")
#' x <- g[[2]]; y <- g[[3]]
#' f <- FpGroup(g[[1]], list(x^4, y^2, x * y * x^-1 * y))
#' fp_order(f)
#' fp_str(fp_subgroup(f, list(x * y)))
#' @export
FpGroup <- function(fr_grp, relators) {
  G <- new.env(parent = emptyenv())
  G$free_group <- fr_grp
  G$relators <- relators
  G$generators <- fr_grp$generators
  G$identity <- fr_grp$identity
  G$ctab <- NULL
  G$is_standardized <- FALSE
  G$order_ <- NULL
  G$perm_iso <- NULL
  # ponytail: SymPy builds the rewriting system in __init__; it depends only
  # on the relators, so it is built on first use with the same result
  G$rws_ <- NULL
  class(G) <- "FpGroup"
  G
}

#' @rdname FpGroup
#' @export
fp_group <- function(fr_grp, relators = list()) {
  G <- FpGroup(fr_grp, relators)
  c(list(G), as.list(G$generators))
}

#' @rdname FpGroup
#' @export
xfp_group <- function(fr_grp, relators = list()) {
  G <- FpGroup(fr_grp, relators)
  list(G, G$generators)
}

.fpg_rws <- function(G) {
  if (is.null(G$rws_)) G$rws_ <- RewritingSystem(G)
  G$rws_
}

.fpg_str <- function(G) {
  if (G$free_group$rank > 30L) return(sprintf("<fp group with %s generators>", G$free_group$rank))
  sprintf("<fp group on the generators %s>", .py_str(G$generators))
}

#' @export
print.FpGroup <- function(x, ...) {
  cat(.fpg_str(x), "\n")
  invisible(x)
}

#' @rdname FpGroup
#' @export
fp_make_confluent <- function(G) {
  rws_make_confluent(.fpg_rws(G))
  invisible(NULL)
}

#' @rdname FpGroup
#' @export
fp_reduce <- function(G, word) rws_reduce(.fpg_rws(G), word)

#' @rdname FpGroup
#' @export
fp_equals <- function(G, word1, word2) {
  if (.fw_eq(fp_reduce(G, .fw_mul(word1, .fw_inv(word2))), G$identity)) return(TRUE)
  if (rws_is_confluent(.fpg_rws(G))) return(FALSE)
  NULL
}

.fpg_check_words <- function(G, gens) {
  if (!all(vapply(gens, function(g) inherits(g, "FreeGroupElement"), NA)))
    .fg_stop("ValueError", "Generators must be `FreeGroupElement`s")
  if (!all(vapply(gens, function(g) identical(g$group, G$free_group), NA)))
    .fg_stop("ValueError", "Given generators are not members of the group")
}

#' @rdname FpGroup
#' @export
fp_subgroup <- function(G, gens, C = NULL, homomorphism = FALSE) {
  .fpg_check_words(G, gens)
  r <- reidemeister_presentation(G, gens, C = C, homomorphism = homomorphism)
  g <- r[[1L]]
  rels <- r[[2L]]
  K <- if (length(g)) FpGroup(g[[1L]]$group, rels) else FpGroup(FreeGroup(character()), list())
  if (homomorphism) return(list(K, homomorphism(K, G, K$generators, r[[3L]], check = FALSE)))
  K
}

#' @rdname FpGroup
#' @export
fp_quotient <- function(G, N) {
  if (inherits(N, "FpSubgroup")) {
    if (!N$normal) .fg_stop("ValueError", "FpSubgroup must be normal")
    if (!(identical(N$parent, G) || identical(N$parent, G$free_group)))
      .fg_stop("ValueError", "Subgroup is not a subgroup of the group")
    gens <- N$generators
  } else {
    gens <- as.list(N)
  }
  .fpg_check_words(G, gens)
  FpGroup(G$free_group, c(as.list(G$relators), gens))
}

#' @rdname FpGroup
#' @export
fp_coset_enumeration <- function(G, H, strategy = "relator_based", max_cosets = NULL, draft = NULL,
                                 incomplete = FALSE) {
  if (is.null(max_cosets) || !max_cosets) max_cosets <- .ct_opts$coset_table_max_limit
  C <- if (strategy == "relator_based") {
    coset_enumeration_r(G, H, max_cosets = max_cosets, draft = draft, incomplete = incomplete)
  } else {
    coset_enumeration_c(G, H, max_cosets = max_cosets, draft = draft, incomplete = incomplete)
  }
  if (ct_is_complete(C)) ct_compress(C)
  C
}

.fpg_standardize_coset_table <- function(G) {
  ct_standardize(G$ctab)
  G$is_standardized <- TRUE
}

#' @rdname FpGroup
#' @export
fp_coset_table <- function(G, H, strategy = "relator_based", max_cosets = NULL, draft = NULL,
                           incomplete = FALSE) {
  if (!length(H)) {
    if (!is.null(G$ctab)) {
      if (!G$is_standardized) .fpg_standardize_coset_table(G)
    } else {
      G$ctab <- fp_coset_enumeration(G, list(), strategy, max_cosets = max_cosets, draft = draft,
                                     incomplete = incomplete)
      .fpg_standardize_coset_table(G)
    }
    return(G$ctab$table)
  }
  C <- fp_coset_enumeration(G, H, strategy, max_cosets = max_cosets, draft = draft, incomplete = incomplete)
  ct_standardize(C)
  C$table
}

.int_gcd <- function(v) Reduce(function(a, b) {
  a <- abs(a)
  b <- abs(b)
  while (b) {
    t <- a %% b
    a <- b
    b <- t
  }
  a
}, v)

#' @rdname FpGroup
#' @export
fp_order <- function(G, strategy = "relator_based") {
  if (.is_fg(G)) return(if (G$rank == 0L) S_One else S_Infinity)
  if (.is_pg(G)) return(.ppg_order(G))
  if (!is.null(G$order_)) return(G$order_)
  if (!is.null(G$ctab)) {
    G$order_ <- length(G$ctab$table)
  } else if (!length(G$relators)) {
    G$order_ <- fp_order(G$free_group)
  } else if (length(G$generators) == 1L) {
    G$order_ <- .int(abs(.int_gcd(vapply(G$relators, function(r) {
      if (!length(r$e)) .fg_stop("IndexError", "tuple index out of range")
      r$e[1L]
    }, 0L))))
  } else if (isTRUE(.fpg_is_infinite(G))) {
    G$order_ <- S_Infinity
  } else {
    r <- .fpg_finite_index_subgroup(G)
    if (!is.null(r$C)) {
      ind <- length(r$C$table)
      G$order_ <- .py_mul(ind, fp_order(fp_subgroup(G, r$gens, C = r$C)))
    } else {
      G$order_ <- fp_index(G, list())
    }
  }
  G$order_
}

# rank of an integer matrix, by exact elimination over QQ
.int_matrix_rank <- function(rows) {
  if (!length(rows)) return(0L)
  M <- do.call(rbind, lapply(rows, function(r) gmp::as.bigq(r)))
  M <- gmp::as.bigq(M)
  m <- nrow(M)
  n <- ncol(M)
  r <- 0L
  for (c in seq_len(n)) {
    if (r == m) break
    piv <- NA_integer_
    for (i in (r + 1L):m) if (M[i, c] != 0) {
      piv <- i
      break
    }
    if (is.na(piv)) next
    r <- r + 1L
    if (piv != r) {
      tmp <- M[r, ]
      M[r, ] <- M[piv, ]
      M[piv, ] <- tmp
    }
    if (r < m) for (i in (r + 1L):m) if (M[i, c] != 0) M[i, ] <- M[i, ] - M[i, c] / M[r, c] * M[r, ]
  }
  r
}

.fpg_is_infinite <- function(G) {
  used <- unique(unlist(lapply(G$relators, function(r) r$s)))
  if (!all(seq_along(G$generators) %in% used)) return(TRUE)
  rows <- lapply(G$relators, function(rel) vapply(G$generators, function(g) exponent_sum(rel, g), 0L))
  # 0 is an invariant factor exactly when the rank is below min(m, n)
  if (.int_matrix_rank(rows) < min(length(rows), length(G$generators))) return(TRUE)
  NULL
}

.fpg_finite_index_subgroup <- function(G, s = NULL) {
  gen <- fp_most_frequent_generator(G)
  rels <- c(as.list(G$generators), as.list(G$relators))
  others <- Filter(function(g) !.fw_eq(g, gen), G$generators)
  if (!length(s)) {
    if (length(G$generators) == 2L) {
      s <- c(list(gen), others)
    } else {
      rand <- G$free_group$identity
      i <- 0L
      while ((.fw_in(rand, rels) || .fw_in(.fw_inv(rand), rels) || .fw_is_id(rand)) && i < 10L) {
        rand <- fp_random(G)
        i <- i + 1L
      }
      s <- c(list(gen, rand), others)
    }
  }
  mid <- (length(s) + 1L) %/% 2L
  half1 <- s[seq_len(mid)]
  half2 <- s[-seq_len(mid)]
  draft1 <- NULL
  draft2 <- NULL
  m <- 200L
  C <- NULL
  half <- NULL
  while (is.null(C) && (m / 2 < .ct_opts$coset_table_max_limit)) {
    m <- min(m, .ct_opts$coset_table_max_limit)
    draft1 <- fp_coset_enumeration(G, half1, max_cosets = m, draft = draft1, incomplete = TRUE)
    if (ct_is_complete(draft1)) {
      C <- draft1
      half <- half1
    } else {
      draft2 <- fp_coset_enumeration(G, half2, max_cosets = m, draft = draft2, incomplete = TRUE)
      if (ct_is_complete(draft2)) {
        C <- draft2
        half <- half2
      }
    }
    if (is.null(C)) m <- m * 2L
  }
  if (is.null(C)) return(list(gens = NULL, C = NULL))
  ct_compress(C)
  list(gens = half, C = C)
}

#' @rdname FpGroup
#' @export
fp_most_frequent_generator <- function(G) {
  gens <- G$generators
  freqs <- vapply(gens, function(g) sum(vapply(G$relators, function(r) generator_count(r, g), 0L)), 0L)
  gens[[which.max(freqs)]]
}

#' @rdname FpGroup
#' @export
fp_random <- function(G) {
  gens <- G$generators
  k <- sample(2:3, 1L)
  .fw_prod(G$free_group, lapply(seq_len(k), function(i) .fw_pow(gens[[sample.int(length(gens), 1L)]], sample(c(1L, -1L), 1L))))
}

#' @rdname FpGroup
#' @export
fp_index <- function(G, H, strategy = "relator_based") {
  if (!length(H)) return(fp_order(G))
  length(fp_coset_enumeration(G, H, strategy)$table)
}

#' Subgroups of finitely presented and free groups
#'
#' `FpSubgroup(G, gens, normal)` is the subgroup of the `FpGroup` or free
#' group `G` generated by `gens` (or the normal closure, if `normal`).
#' `fps_contains()` tests membership (by Stallings foldings in a free group
#' and by coset enumeration or rewriting in an `FpGroup`), `fps_order()`
#' gives the order and `fps_to_FpGroup()` a presentation.
#'
#' @param G An `FpGroup` or a free group.
#' @param gens A list of words.
#' @param normal Take the normal closure.
#' @param S An `FpSubgroup`.
#' @param g A word.
#' @return An `FpSubgroup` (an environment), a logical, an order or a group.
#' @examples
#' g <- free_group("x, y")
#' x <- g[[2]]; y <- g[[3]]
#' f <- FpGroup(g[[1]], list(x^4, y^2, x * y * x^-1 * y))
#' S <- FpSubgroup(f, list(x * y))
#' fps_contains(S, (x * y)^-3)
#' fps_order(S)
#' @export
FpSubgroup <- function(G, gens, normal = FALSE) {
  S <- new.env(parent = emptyenv())
  S$parent <- G
  S$generators <- .fw_unique(Filter(function(g) !.fw_eq(g, G$identity), gens))
  S$C <- NULL
  S$normal <- normal
  class(S) <- "FpSubgroup"
  S
}

.fps_str <- function(S) {
  if (length(S$generators) > 30L) return(sprintf("<fp subgroup with %s generators>", length(S$generators)))
  sprintf("<fp subgroup on the generators %s>", .py_str(S$generators))
}

#' @export
print.FpSubgroup <- function(x, ...) {
  cat(.fps_str(x), "\n")
  invisible(x)
}

#' @rdname FpSubgroup
#' @export
fps_contains <- function(S, g) {
  if (S$normal) {
    quotient <- if (.is_fg(S$parent)) FpGroup(S$parent, S$generators) else fp_quotient(S$parent, S)
    return(isTRUE(fp_equals(quotient, g, quotient$identity)))
  }
  if (.is_fg(S$parent)) return(.stallings_contains(S$generators, g))
  if (is.null(S$C)) S$C <- fp_coset_enumeration(S$parent, S$generators)
  i <- 0L
  for (l in .fw_letters(g)) i <- S$C$table[[i + 1L]][.ct_col(l)]
  identical(i, 0L)
}

# Stallings folding of the subgroup graph of a free group (FpSubgroup)
.stallings_contains <- function(gens, g) {
  graph <- list()  # vertex key -> list(label key -> integer targets)
  vkey <- function(u) as.character(u)
  add_edge <- function(u, v, label) {
    ku <- vkey(u)
    kv <- vkey(v)
    lk <- as.character(label)
    if (is.null(graph[[ku]])) graph[[ku]] <<- list()
    graph[[ku]][[lk]] <<- unique(c(graph[[ku]][[lk]], v))
    li <- as.character(-label)
    if (is.null(graph[[kv]])) graph[[kv]] <<- list()
    graph[[kv]][[li]] <<- unique(c(graph[[kv]][[li]], u))
  }
  next_vertex_id <- 1L
  for (gen in gens) {
    current <- 0L
    letters <- .fw_letters(gen)
    for (i in seq_along(letters)) {
      if (i == length(letters)) {
        nxt <- 0L
      } else {
        nxt <- next_vertex_id
        next_vertex_id <- next_vertex_id + 1L
      }
      add_edge(current, nxt, letters[i])
      current <- nxt
    }
  }
  merge_vertices <- function(u, v) {
    if (u == v) return(invisible())
    for (node in names(graph)) {
      if (node == vkey(v)) {
        for (lk in names(graph[[node]])) for (target in graph[[node]][[lk]]) {
          if (target == v) target <- u
          add_edge(u, target, as.integer(lk))
        }
        graph[[vkey(v)]] <<- NULL
      } else {
        for (lk in names(graph[[node]])) {
          tg <- graph[[node]][[lk]]
          if (v %in% tg) {
            graph[[node]][[lk]] <<- unique(c(tg[tg != v], u))
            add_edge(as.integer(node), u, as.integer(lk))
          }
        }
      }
    }
  }
  repeat {
    folded <- TRUE
    for (u in names(graph)) {
      for (lk in names(graph[[u]])) {
        tg <- graph[[u]][[lk]]
        if (length(tg) > 1L) {
          st <- sort(tg)[1:2]
          merge_vertices(min(st), max(st))
          folded <- FALSE
          break
        }
      }
      if (!folded) break
    }
    if (folded) break
  }
  current <- 0L
  for (l in .fw_letters(g)) {
    tg <- graph[[vkey(current)]][[as.character(l)]]
    if (!length(tg)) return(FALSE)
    current <- tg[1L]
  }
  current == 0L
}

#' @rdname FpSubgroup
#' @export
fps_order <- function(S) {
  if (!length(S$generators)) return(S_One)
  if (.is_fg(S$parent)) return(S_Infinity)
  if (is.null(S$C)) S$C <- fp_coset_enumeration(S$parent, S$generators)
  .py_div(fp_order(S$parent), length(S$C$table))
}

#' @rdname FpSubgroup
#' @export
fps_to_FpGroup <- function(S) {
  if (.is_fg(S$parent)) return(FreeGroup(paste0("x_", seq_along(S$generators) - 1L)))
  fp_subgroup(S$parent, S$generators, C = S$C)
}

#' Low index subgroups
#'
#' Finds the subgroups of index at most `N` of the `FpGroup` `G` that contain
#' the words `Y` (each a generator or inverse of a generator), returning
#' their coset tables, as SymPy's `low_index_subgroups()`.
#'
#' @param G An `FpGroup`.
#' @param N The largest index.
#' @param Y A list of words.
#' @return A list of coset tables.
#' @examples
#' g <- free_group("x, y")
#' x <- g[[2]]; y <- g[[3]]
#' f <- FpGroup(g[[1]], list(x^2, y^3, (x * y)^4))
#' L <- low_index_subgroups(f, 4)
#' fp_str(lapply(L, function(C) C$table))
#' @export
low_index_subgroups <- function(G, N, Y = list()) {
  C <- CosetTable(G, list())
  R <- G$relators
  len_short_rel <- 5L
  R2 <- .fw_unique(Filter(function(rel) .fw_len(rel) > len_short_rel, R))
  R1 <- .fw_unique(lapply(Filter(function(rel) .fw_len(rel) <= len_short_rel, .fw_unique(R)),
                          identity_cyclic_reduction))
  R1_c_list <- ct_conjugates(C, R1)
  env <- new.env(parent = emptyenv())
  env$S <- list()
  .descendant_subgroups(env, C, R1_c_list, C$A[1L], R2, N, Y)
  env$S
}

.descendant_subgroups <- function(env, C, R1_c_list, x, R2, N, Y) {
  if (ct_is_complete(C)) {
    for (w in R2) for (alpha in ct_omega(C)) if (!ct_scan_check(C, alpha, w)) return(invisible())
    env$S[[length(env$S) + 1L]] <- C
  } else {
    found <- FALSE
    for (alpha in seq_along(C$table) - 1L) {
      for (xx in C$A) if (is.na(C$table[[alpha + 1L]][.ct_col(xx)])) {
        undefined_coset <- alpha
        undefined_gen <- xx
        found <- TRUE
        break
      }
      if (found) break
    }
    n <- ct_n(C)
    reach <- c(ct_omega(C), n)
    for (beta in reach) {
      if (beta < N) {
        if (beta == n || is.na(C$table[[beta + 1L]][.ct_col(-undefined_gen)]))
          .try_descendant(env, C, R1_c_list, R2, N, undefined_coset, undefined_gen, beta, Y)
      }
    }
  }
  invisible()
}

.try_descendant <- function(env, C, R1_c_list, R2, N, alpha, x, beta, Y) {
  D <- ct_copy(C)
  if (beta == ct_n(D) && beta < N) {
    D$table[[length(D$table) + 1L]] <- rep(NA_integer_, D$nA)
    D$p[length(D$p) + 1L] <- beta
  }
  D$table[[alpha + 1L]][.ct_col(x)] <- beta
  D$table[[beta + 1L]][.ct_col(-x)] <- alpha
  D$deduction_stack[[length(D$deduction_stack) + 1L]] <- c(alpha, x)
  if (!ct_process_deductions_check(D, R1_c_list[[.ct_col(x)]], R1_c_list[[.ct_col(-x)]])) return(invisible())
  for (w in Y) if (!ct_scan_check(D, 0L, w)) return(invisible())
  if (first_in_class(D, Y)) .descendant_subgroups(env, D, R1_c_list, x, R2, N, Y)
  invisible()
}

#' @rdname low_index_subgroups
#' @param C A coset table.
#' @export
first_in_class <- function(C, Y = list()) {
  n <- ct_n(C)
  lamda <- -1L
  nu <- rep(NA_integer_, n)
  mu <- rep(NA_integer_, n)
  next_alpha <- FALSE
  ylet <- vapply(Y, .ct_letter, 0L)
  for (alpha in seq_len(n - 1L)) {
    for (beta in seq_len(lamda + 1L) - 1L) nu[mu[beta + 1L] + 1L] <- NA_integer_
    for (w in ylet) {
      t <- C$table[[alpha + 1L]][.ct_col(w)]
      if (is.na(t) || t != alpha) {
        next_alpha <- TRUE
        break
      }
    }
    if (next_alpha) {
      next_alpha <- FALSE
      next
    }
    mu[1L] <- alpha
    nu[alpha + 1L] <- 0L
    lamda <- 0L
    for (beta in seq_len(n) - 1L) {
      if (is.na(mu[beta + 1L])) .fg_stop("TypeError", "list indices must be integers or slices, not NoneType")
      for (x in C$A) {
        gamma <- C$table[[beta + 1L]][.ct_col(x)]
        delta <- C$table[[mu[beta + 1L] + 1L]][.ct_col(x)]
        if (is.na(gamma) || is.na(delta)) {
          next_alpha <- TRUE
          break
        }
        if (is.na(nu[delta + 1L])) {
          lamda <- lamda + 1L
          nu[delta + 1L] <- lamda
          mu[lamda + 1L] <- delta
        }
        if (nu[delta + 1L] < gamma) return(FALSE)
        if (nu[delta + 1L] > gamma) {
          next_alpha <- TRUE
          break
        }
      }
      if (next_alpha) {
        next_alpha <- FALSE
        break
      }
    }
  }
  TRUE
}

#' Simplifying presentations
#'
#' `simplify_presentation(G)` returns a simplified isomorphic copy of the
#' `FpGroup` `G`; `simplify_presentation(gens, rels)` simplifies the lists of
#' generators and relators and returns both. With `change_gens = TRUE`
#' redundant generators are removed. `elimination_technique_1()` and
#' `elimination_technique_2()` are the Tietze transformations SymPy uses.
#'
#' @param ... An `FpGroup`, or the generators and the relators.
#' @param change_gens Remove redundant generators.
#' @param gens,rels Lists of words.
#' @param identity The identity of the free group.
#' @param C A coset table with Reidemeister data.
#' @return An `FpGroup` or a list of generators and relators.
#' @examples
#' g <- free_group("x, y")
#' x <- g[[2]]; y <- g[[3]]
#' fp_str(simplify_presentation(FpGroup(g[[1]], list(x^2, x^5, y^3)))$relators)
#' @export
simplify_presentation <- function(..., change_gens = FALSE) {
  args <- list(...)
  if (length(args) == 1L) {
    if (!.is_fpg(args[[1L]])) .fg_stop("TypeError", "The argument must be an instance of FpGroup")
    G <- args[[1L]]
    r <- simplify_presentation(G$generators, G$relators, change_gens = change_gens)
    if (length(r[[1L]])) return(FpGroup(r[[1L]][[1L]]$group, r[[2L]]))
    return(FpGroup(FreeGroup(character()), list()))
  } else if (length(args) == 2L) {
    gens <- args[[1L]]
    rels <- as.list(args[[2L]])
    if (!length(gens)) return(.py_tuple(list(gens, rels)))
    identity <- gens[[1L]]$group$identity
  } else {
    .fg_stop("RuntimeError", if (!length(args)) "Not enough arguments" else "Too many arguments")
  }
  prev_gens <- list()
  prev_rels <- list()
  sameset <- function(a, b) setequal(if (length(a)) .fw_keys(a) else character(),
                                     if (length(b)) .fw_keys(b) else character())
  while (!sameset(prev_rels, rels)) {
    prev_rels <- rels
    while (change_gens && !sameset(prev_gens, gens)) {
      prev_gens <- gens
      r <- elimination_technique_1(gens, rels, identity)
      gens <- r[[1L]]
      rels <- r[[2L]]
    }
    rels <- .simplify_relators(rels)
  }
  if (change_gens) {
    syms <- vapply(gens, function(g) g$group$symbols[g$s[1L]], "")
    Fn <- FreeGroup(syms)
    gens <- Fn$generators
    for (j in seq_along(rels)) {
      a <- rels[[j]]
      oldsyms <- a$group$symbols[a$s]
      idx <- match(oldsyms, syms)
      if (anyNA(idx)) .fg_stop("KeyError", oldsyms[is.na(idx)][1L])
      rels[[j]] <- .fw_prod(Fn, lapply(seq_along(idx), function(i) .fw_pow(gens[[idx[i]]], a$e[i])))
    }
  }
  .py_tuple(list(gens, rels))
}

.simplify_relators <- function(rels) {
  rels <- as.list(rels)
  if (!length(rels)) return(list())
  identity <- rels[[1L]]$group$identity
  exps <- .od()
  for (i in seq_along(rels)) {
    rel <- rels[[i]]
    if (length(rel$s) == 1L) {
      g <- .fw_letter(rel, 0L)
      ex <- abs(rel$e[1L])
      if (rel$e[1L] < 0L) {
        rels[[i]] <- .fw_inv(rels[[i]])
        g <- .fw_inv(g)
      }
      k <- .fw_key(g)
      if (.od_has(exps, k)) ex <- .int_gcd(c(ex, .od_get(exps, k)$e[1L]))
      exps <- .od_set(exps, k, g, .fw_pow(g, ex))
    }
  }
  one_syllables_words <- exps$v
  for (i in seq_along(rels)) {
    rel <- rels[[i]]
    if (.fw_in(rel, one_syllables_words)) next
    rel <- eliminate_words(rel, one_syllables_words, all = TRUE)
    for (g in contains_generators(rel)) {
      k <- .fw_key(g)
      if (.od_has(exps, k)) {
        ex <- .od_get(exps, k)$e[1L]
        max_exp <- (ex + 1L) %/% 2L
        rel <- eliminate_word(rel, .fw_pow(g, max_exp), .fw_pow(g, max_exp - ex), all = TRUE)
        rel <- eliminate_word(rel, .fw_pow(g, -max_exp), .fw_pow(g, -(max_exp - ex)), all = TRUE)
      }
    }
    rels[[i]] <- rel
  }
  rels <- lapply(rels, identity_cyclic_reduction)
  rels <- c(rels, one_syllables_words)
  rels <- .fw_sort(.fw_unique(rels))
  for (i in seq_along(rels)) if (.fw_eq(rels[[i]], identity)) {
    rels <- rels[-i]
    break
  }
  rels
}

#' @rdname simplify_presentation
#' @export
elimination_technique_1 <- function(gens, rels, identity) {
  rels <- .fw_sort(rels)
  gens <- as.list(gens)
  redundant <- .od()
  redundant_rels <- list()
  used_gens <- character()
  for (rel in rels) {
    contained_gens <- contains_generators(rel)
    if (any(.fw_keys(contained_gens) %in% redundant$k)) next
    contained_gens <- .fw_sort(contained_gens, decreasing = TRUE)
    for (gen in contained_gens) {
      if (generator_count(rel, gen) == 1L && !(.fw_key(gen) %in% used_gens)) {
        k <- exponent_sum(rel, gen)
        gen_index <- word_index(rel, .fw_pow(gen, k))
        bk <- subword(rel, gen_index + 1L, .fw_len(rel))
        fw <- subword(rel, 0L, gen_index)
        chi <- .fw_mul(bk, fw)
        redundant <- .od_set(redundant, .fw_key(gen), gen, .fw_pow(chi, -1L * k))
        cg <- contains_generators(chi)
        if (length(cg)) used_gens <- union(used_gens, .fw_keys(cg))
        redundant_rels[[length(redundant_rels) + 1L]] <- rel
        break
      }
    }
  }
  rels <- Filter(function(r) !.fw_in(r, redundant_rels), rels)
  rels <- lapply(rels, function(r) identity_cyclic_reduction(
    eliminate_words(r, list(keys = redundant$key, values = redundant$v), all = TRUE)))
  rels <- .fw_unique(rels)
  for (i in seq_along(rels)) if (.fw_eq(rels[[i]], identity)) {
    rels <- rels[-i]
    break
  }
  gens <- Filter(function(g) !(.fw_key(g) %in% redundant$k), gens)
  .py_tuple(list(gens, rels))
}

# Pg 175 [1]
.define_schreier_generators <- function(C, homomorphism = FALSE) {
  y <- character()
  gamma <- 1L
  f <- C$fp_group
  n <- ct_n(C)
  if (homomorphism) {
    gens_ <- .od()
    tau <- vector("list", n)
    tau[[1L]] <- f$identity
  }
  C$P <- lapply(seq_len(n), function(i) vector("list", C$nA))
  Fg <- f$identity$group
  for (alpha in ct_omega(C)) for (x in C$A) {
    cx <- .ct_col(x)
    beta <- C$table[[alpha + 1L]][cx]
    if (!is.na(beta) && beta == gamma) {
      C$P[[alpha + 1L]][[cx]] <- "<identity>"
      C$P[[beta + 1L]][[.ct_col(-x)]] <- "<identity>"
      gamma <- gamma + 1L
      if (homomorphism) tau[[beta + 1L]] <- .fw_mul(tau[[alpha + 1L]], .fw(Fg, abs(x), sign(x)))
    } else if (x > 0L && is.null(C$P[[alpha + 1L]][[cx]])) {
      y_alpha_x <- paste0(Fg$symbols[x], "_", alpha)
      y <- c(y, y_alpha_x)
      C$P[[alpha + 1L]][[cx]] <- y_alpha_x
      if (homomorphism) {
        gens_ <- .od_set(gens_, y_alpha_x, y_alpha_x,
                         .fw_mul(.fw_mul(tau[[alpha + 1L]], .fw(Fg, x, 1L)), .fw_inv(tau[[beta + 1L]])))
      }
    }
  }
  SF <- FreeGroup(y)
  C$schreier_free_group <- SF
  C$schreier_generators_ <- as.list(SF$generators)
  if (homomorphism) C$schreier_gen_elem <- gens_
  for (i in seq_along(C$P)) for (j in seq_len(C$nA)) {
    v <- C$P[[i]][[j]]
    if (is.character(v)) {
      if (v == "<identity>") {
        C$P[[i]][[j]] <- SF$identity
      } else {
        r <- C$schreier_generators_[[match(v, y)]]
        C$P[[i]][[j]] <- r
        beta <- C$table[[i]][j]
        C$P[[beta + 1L]][[j + 1L]] <- .fw_inv(r)
      }
    }
  }
  invisible(NULL)
}

.rewrite <- function(C, alpha, w) {
  factors <- list()
  for (l in .fw_letters(w)) {
    j <- .ct_col(l)
    factors[[length(factors) + 1L]] <- C$P[[alpha + 1L]][[j]]
    alpha <- C$table[[alpha + 1L]][j]
  }
  .fw_prod(C$schreier_free_group, factors)
}

.reidemeister_relators <- function(C) {
  R <- C$fp_group$relators
  n <- ct_n(C)
  rels <- list()
  for (word in R) for (coset in seq_len(n) - 1L) rels[[length(rels) + 1L]] <- .rewrite(C, coset, word)
  order_1_gens <- .fw_unique(Filter(function(i) .fw_len(i) == 1L, rels))
  rels <- Filter(function(rel) !.fw_in(rel, order_1_gens), rels)
  rels <- lapply(rels, function(w) eliminate_words(w, order_1_gens, all = TRUE))
  C$schreier_generators_ <- Filter(function(i) !(.fw_in(i, order_1_gens) || .fw_in(.fw_inv(i), order_1_gens)),
                                   C$schreier_generators_)
  i <- 1L
  while (i <= length(rels)) {
    w <- rels[[i]]
    j <- i + 1L
    while (j <= length(rels)) {
      if (is_cyclic_conjugate(w, rels[[j]])) rels <- rels[-j] else j <- j + 1L
    }
    i <- i + 1L
  }
  C$reidemeister_relators_ <- rels
  invisible(NULL)
}

#' @rdname simplify_presentation
#' @export
elimination_technique_2 <- function(C) {
  rels <- .fw_sort(C$reidemeister_relators_, decreasing = TRUE)
  gens <- C$schreier_generators_
  for (i in rev(seq_along(gens))) {
    if (i > length(rels)) .fg_stop("IndexError", "list index out of range")
    rel <- rels[[i]]
    for (j in rev(seq_along(gens))) {
      gen <- gens[[j]]
      if (generator_count(rel, gen) == 1L) {
        k <- exponent_sum(rel, gen)
        gen_index <- word_index(rel, .fw_pow(gen, k))
        bk <- subword(rel, gen_index + 1L, .fw_len(rel))
        fw <- subword(rel, 0L, gen_index)
        rep_by <- .fw_pow(.fw_mul(bk, fw), -1L * k)
        rels <- rels[-i]
        gens <- gens[-j]
        rels <- lapply(rels, function(r) eliminate_word(r, gen, rep_by))
        break
      }
    }
  }
  C$reidemeister_relators_ <- rels
  C$schreier_generators_ <- gens
  .py_tuple(list(gens, rels))
}

#' Reidemeister Schreier presentations
#'
#' Returns the generators and relators of a presentation of the subgroup of
#' the `FpGroup` `fp_grp` generated by the words `H`, simplified as in
#' SymPy's `reidemeister_presentation()`. With `homomorphism = TRUE` the
#' images of the new generators in `fp_grp` are returned as a third element.
#'
#' @param fp_grp An `FpGroup`.
#' @param H A list of words.
#' @param C A coset table of `H` to use.
#' @param homomorphism Also return the images of the generators.
#' @return A list of two (or three) tuples.
#' @examples
#' g <- free_group("x, y")
#' x <- g[[2]]; y <- g[[3]]
#' f <- FpGroup(g[[1]], list(x^3, y^5, (x * y)^2))
#' fp_str(reidemeister_presentation(f, list(x * y, x^-1 * y^-1 * x * y * x)))
#' @export
reidemeister_presentation <- function(fp_grp, H, C = NULL, homomorphism = FALSE) {
  if (is.null(C)) C <- coset_enumeration_r(fp_grp, H)
  ct_compress(C)
  ct_standardize(C)
  .define_schreier_generators(C, homomorphism = homomorphism)
  .reidemeister_relators(C)
  r <- simplify_presentation(C$schreier_generators_, C$reidemeister_relators_, change_gens = TRUE)
  C$schreier_generators <- .py_tuple(r[[1L]])
  C$reidemeister_relators <- .py_tuple(r[[2L]])
  if (homomorphism) {
    gens_ <- lapply(r[[1L]], function(g) .od_get(C$schreier_gen_elem, .fw_str(g)))
    return(list(C$schreier_generators, C$reidemeister_relators, gens_))
  }
  .py_tuple(list(C$schreier_generators, C$reidemeister_relators))
}

# ---------------------------------------------------------------------------
# Rewriting systems (rewritingsystem.py, rewritingsystem_fsm.py)

#' Rewriting systems for finitely presented groups
#'
#' `RewritingSystem(G)` builds the rewriting system SymPy attaches to an
#' `FpGroup` (also reachable with `fp_rewriting_system(G)`), with its rules
#' (`R$rules`, an ordered dictionary of left and right hand sides) and its
#' reduction automaton. `rws_reduce()` reduces a word with the rules,
#' `rws_reduce_using_automaton()` with the automaton,
#' `rws_make_confluent()` runs the Knuth Bendix completion,
#' `rws_is_confluent()` tests confluence, `rws_add_rule()` adds a rule (and
#' the rules its overlaps imply) and `rws_set_max()` sets the largest number
#' of rules. `rws_rules()` lists the rules as pairs.
#'
#' @param group,G An `FpGroup`.
#' @param R A rewriting system.
#' @param word,w1,w2 Words.
#' @param exclude A rule left hand side to skip.
#' @param check Only check, without adding rules.
#' @param n The largest number of rules.
#' @return A `RewritingSystem` (an environment), a word, a logical or a list.
#' @examples
#' g <- free_group("a, b")
#' a <- g[[2]]; b <- g[[3]]
#' G <- FpGroup(g[[1]], list(a * b * a^-1 * b^-1))
#' R <- fp_rewriting_system(G)
#' rws_is_confluent(R)
#' fp_str(fp_reduce(G, b^-1 * a))
#' @export
RewritingSystem <- function(group) {
  R <- new.env(parent = emptyenv())
  R$group <- group
  R$alphabet <- group$generators
  R$is_confluent_ <- NULL
  R$maxeqns <- 32767L
  R$tidyint <- 100L
  R$max_exceeded <- FALSE
  R$automaton <- NULL
  R$new_rules <- .od()
  R$rules <- .od()
  R$rules_cache <- character()
  class(R) <- "RewritingSystem"
  identity <- group$free_group$identity
  for (r in group$relators) rws_add_rule(R, r, identity)
  .rws_remove_redundancies(R)
  gens <- as.list(R$alphabet)
  gens <- c(gens, lapply(gens, .fw_inv))
  R$automaton <- .sm_new(paste0("Reduction automaton for ", .py_str(group)), gens)
  .rws_add_to_automaton(R, "rules")
  R
}

#' @rdname RewritingSystem
#' @export
fp_rewriting_system <- function(G) .fpg_rws(G)

#' @export
print.RewritingSystem <- function(x, ...) {
  cat("RewritingSystem with", length(x$rules$k), "rules\n")
  invisible(x)
}

#' @rdname RewritingSystem
#' @export
rws_rules <- function(R) lapply(seq_along(R$rules$k), function(i) .py_tuple(list(R$rules$key[[i]], R$rules$v[[i]])))

#' @rdname RewritingSystem
#' @export
rws_set_max <- function(R, n) {
  if (n > R$maxeqns) R$max_exceeded <- FALSE
  R$maxeqns <- n
  invisible(NULL)
}

#' @rdname RewritingSystem
#' @export
rws_is_confluent <- function(R) {
  if (is.null(R$is_confluent_)) R$is_confluent_ <- .rws_check_confluence(R)
  R$is_confluent_
}

.rws_add_rule_raw <- function(R, r1, r2) {
  if (length(R$rules$k) + 1L > R$maxeqns) {
    R$is_confluent_ <- .rws_check_confluence(R)
    R$max_exceeded <- TRUE
    .fg_stop("RuntimeError", "Too many rules were defined.")
  }
  k <- .fw_key(r1)
  R$rules <- .od_set(R$rules, k, r1, r2)
  if (!is.null(R$automaton)) R$new_rules <- .od_set(R$new_rules, k, r1, r2)
}

# a Python set of words: unique by key, insertion order
.wset_add <- function(s, w) {
  k <- .fw_key(w)
  if (!(k %in% s$k)) {
    s$k <- c(s$k, k)
    s$w[[length(s$w) + 1L]] <- w
  }
  s
}
.wset <- function() list(k = character(), w = list())

#' @rdname RewritingSystem
#' @export
rws_add_rule <- function(R, w1, w2, check = FALSE) .py_set(.rws_add_rule(R, w1, w2, check)$w)

.rws_add_rule <- function(R, w1, w2, check = FALSE) {
  new_keys <- .wset()
  if (.fw_eq(w1, w2)) return(new_keys)
  if (.fw_lt(w1, w2)) {
    t <- w1
    w1 <- w2
    w2 <- t
  }
  ck <- paste(.fw_key(w1), .fw_key(w2), sep = "|")
  if (ck %in% R$rules_cache) return(new_keys)
  R$rules_cache <- utils::tail(c(R$rules_cache, ck), 50L)
  s1 <- w1
  s2 <- w2
  if (.fw_len(s1) - .fw_len(s2) < 3L) {
    if (!.od_has(R$rules, .fw_key(s1))) {
      new_keys <- .wset_add(new_keys, s1)
      if (!check) .rws_add_rule_raw(R, s1, s2)
    }
    if (.fw_lt(.fw_inv(s1), .fw_inv(s2)) && !.od_has(R$rules, .fw_key(.fw_inv(s2)))) {
      new_keys <- .wset_add(new_keys, .fw_inv(s2))
      if (!check) .rws_add_rule_raw(R, .fw_inv(s2), .fw_inv(s1))
    }
  }
  while (.fw_len(s1) - .fw_len(s2) > -1L) {
    g <- .fw_letter(s1, .fw_len(s1) - 1L)
    s1 <- subword(s1, 0L, .fw_len(s1) - 1L)
    s2 <- .fw_mul(s2, .fw_inv(g))
    if (.fw_len(s1) - .fw_len(s2) < 0L) {
      if (!.od_has(R$rules, .fw_key(s2))) {
        if (!check) .rws_add_rule_raw(R, s2, s1)
        new_keys <- .wset_add(new_keys, s2)
      }
    } else if (.fw_len(s1) - .fw_len(s2) < 3L) {
      new <- .rws_add_rule(R, s1, s2, check)
      for (w in new$w) new_keys <- .wset_add(new_keys, w)
    }
  }
  while (.fw_len(w1) - .fw_len(w2) > -1L) {
    g <- .fw_letter(w1, 0L)
    w1 <- subword(w1, 1L, .fw_len(w1))
    w2 <- .fw_mul(.fw_inv(g), w2)
    if (.fw_len(w1) - .fw_len(w2) < 0L) {
      if (!.od_has(R$rules, .fw_key(w2))) {
        if (!check) .rws_add_rule_raw(R, w2, w1)
        new_keys <- .wset_add(new_keys, w2)
      }
    } else if (.fw_len(w1) - .fw_len(w2) < 3L) {
      new <- .rws_add_rule(R, w1, w2, check)
      for (w in new$w) new_keys <- .wset_add(new_keys, w)
    }
  }
  new_keys
}

.rws_remove_redundancies <- function(R, changes = FALSE) {
  removed <- .wset()
  added <- .wset()
  rules <- R$rules
  for (i in seq_along(rules$k)) {
    r <- rules$key[[i]]
    v <- rws_reduce(R, r, exclude = r)
    w <- rws_reduce(R, rules$v[[i]])
    if (!.fw_eq(v, r)) {
      R$rules <- .od_del(R$rules, rules$k[i])
      removed <- .wset_add(removed, r)
      if (.fw_lt(w, v)) {
        added <- .wset_add(added, v)
        R$rules <- .od_set(R$rules, .fw_key(v), v, w)
      } else if (.fw_lt(v, w)) {
        added <- .wset_add(added, w)
        R$rules <- .od_set(R$rules, .fw_key(w), w, v)
      }
    } else {
      R$rules <- .od_set(R$rules, .fw_key(v), v, w)
    }
  }
  if (changes) return(list(removed = removed, added = added))
  invisible(NULL)
}

#' @rdname RewritingSystem
#' @export
rws_make_confluent <- function(R, check = FALSE) {
  if (R$max_exceeded) return(R$is_confluent_)
  lhs <- R$rules$key
  overlaps_of <- function(r1, r2) {
    len1 <- .fw_len(r1)
    len2 <- .fw_len(r2)
    result <- list()
    for (j in seq_len(len1 + len2 - 1L)) {
      if (.fw_eq(subword(r1, len1 - j, len1 + len2 - j, strict = FALSE),
                 subword(r2, j - len1, j, strict = FALSE))) {
        a <- .fw_mul(subword(r1, 0L, len1 - j, strict = FALSE), subword(r2, 0L, j - len1, strict = FALSE))
        b <- subword(r2, j - len1, j, strict = FALSE)
        c <- .fw_mul(subword(r2, j, len2, strict = FALSE), subword(r1, len1 + len2 - j, len1, strict = FALSE))
        result[[length(result) + 1L]] <- .fw_mul(.fw_mul(a, b), c)
      }
    }
    result
  }
  # returns NULL (no change), "check", a word set, or FALSE
  process_overlap <- function(w, r1, r2) {
    s <- rws_reduce(R, eliminate_word(w, r1, .od_get(R$rules, .fw_key(r1))))
    t <- rws_reduce(R, eliminate_word(w, r2, .od_get(R$rules, .fw_key(r2))))
    if (!.fw_eq(s, t)) {
      if (check) return("check")
      return(tryCatch(.rws_add_rule(R, t, s, check), RuntimeError = function(e) FALSE))
    }
    NULL
  }
  added <- 0L
  i <- 1L
  while (i <= length(lhs)) {
    r1 <- lhs[[i]]
    i <- i + 1L
    j <- 1L
    while (j <= length(lhs)) {
      r2 <- lhs[[j]]
      j <- j + 1L
      if (.fw_eq(r1, r2)) next
      overlaps <- c(overlaps_of(r1, r2), overlaps_of(.fw_inv(r1), r2))
      if (!length(overlaps)) next
      for (w in overlaps) {
        new_keys <- process_overlap(w, r1, r2)
        if (identical(new_keys, "check")) return(FALSE)
        if (is.list(new_keys) && length(new_keys$w)) {
          lhs <- c(lhs, new_keys$w)
          added <- added + length(new_keys$w)
        } else if (isFALSE(new_keys)) {
          return(R$is_confluent_)
        }
      }
      if (added > R$tidyint && !check) {
        ra <- .rws_remove_redundancies(R, changes = TRUE)
        added <- 0L
        lk <- .fw_keys(lhs)
        if (length(ra$removed$k)) i <- min(match(ra$removed$k, lk))
        lhs <- lhs[!(lk %in% ra$removed$k)]
        lhs <- c(lhs, ra$added$w)
        if (.fw_key(r1) %in% ra$removed$k) break
      }
    }
  }
  R$is_confluent_ <- TRUE
  if (!check) .rws_remove_redundancies(R)
  TRUE
}

.rws_check_confluence <- function(R) rws_make_confluent(R, check = TRUE)

#' @rdname RewritingSystem
#' @export
rws_reduce <- function(R, word, exclude = NULL) {
  keep <- if (is.null(exclude)) rep(TRUE, length(R$rules$k)) else R$rules$k != .fw_key(exclude)
  lhs <- R$rules$key[keep]
  rhs <- R$rules$v[keep]
  inv_ok <- vapply(seq_along(lhs), function(i) !.fw_lt(.fw_inv(lhs[[i]]), .fw_inv(rhs[[i]])), NA)
  again <- TRUE
  new <- word
  while (again) {
    again <- FALSE
    for (i in seq_along(lhs)) {
      prev <- new
      new <- eliminate_word(new, lhs[[i]], rhs[[i]], all = TRUE, inverse = inv_ok[i])
      if (!.fw_eq(new, prev)) again <- TRUE
    }
  }
  new
}

.rws_compute_inverse_rules <- function(rules) {
  inv <- .od()
  for (i in seq_along(rules$k)) {
    rk <- .fw_inv(rules$key[[i]])
    rv <- .fw_inv(rules$v[[i]])
    if (.fw_lt(rv, rk)) inv <- .od_set(inv, .fw_key(rk), rk, rv)
    else inv <- .od_set(inv, .fw_key(rv), rv, rk)
  }
  inv
}

# state machine: states live in an environment keyed by the word key (or
# "start"), `order` keeps the insertion order of the keys
.sm_new <- function(name, alphabet) {
  M <- new.env(parent = emptyenv())
  M$name <- name
  M$automaton_alphabet <- alphabet
  M$states <- new.env(parent = emptyenv())
  M$order <- character()
  .sm_add_state(M, "start", NULL, "s")
  M
}

.sm_add_state <- function(M, key, word, type, rh_rule = NULL) {
  if (!(key %in% M$order)) M$order <- c(M$order, key)
  assign(key, list(name = word, type = type, rh_rule = rh_rule, transitions = character()), envir = M$states)
}

.sm_has <- function(M, key) exists(key, envir = M$states, inherits = FALSE)

# _add_to_automaton; `field` names the rule dictionary of R that SymPy
# passes (and updates in place with the inverse rules)
.rws_add_to_automaton <- function(R, field) {
  M <- R$automaton
  all_rules <- R[[field]]
  inv <- .rws_compute_inverse_rules(all_rules)
  for (i in seq_along(inv$k)) all_rules <- .od_set(all_rules, inv$k[i], inv$key[[i]], inv$v[[i]])
  assign(field, all_rules, envir = R)
  alphabet <- .wset()
  accept_states <- character()
  accept_words <- list()
  for (i in seq_along(all_rules$k)) {
    rule <- all_rules$key[[i]]
    rk <- all_rules$k[i]
    letters <- word_letter_form_elm(rule)
    for (l in letters) alphabet <- .wset_add(alphabet, l)
    arr <- letters
    for (t in seq_len(max(0L, length(arr) - 1L))) {
      arr[[t + 1L]] <- .fw_mul(arr[[t]], arr[[t + 1L]])
      elem <- arr[[t]]
      ek <- .fw_key(elem)
      if (!.sm_has(M, ek)) {
        .sm_add_state(M, ek, elem, "a")
        accept_states <- c(accept_states, ek)
        accept_words[[ek]] <- elem
      }
    }
    if (rk %in% accept_states) {
      st <- get(rk, envir = M$states)
      st$type <- "d"
      st$rh_rule <- all_rules$v[[i]]
      assign(rk, st, envir = M$states)
      accept_states <- accept_states[-match(rk, accept_states)]
    }
    if (!.sm_has(M, rk)) .sm_add_state(M, rk, rule, "d", all_rules$v[[i]])
  }
  next_state <- function(word, letter) {
    nxt <- .fw_mul(word, letter)
    while (.fw_len(nxt) && !.sm_has(M, .fw_key(nxt))) nxt <- subword(nxt, 1L, .fw_len(nxt))
    if (!.fw_len(nxt)) "start" else .fw_key(nxt)
  }
  for (key in M$order) {
    st <- get(key, envir = M$states)
    if (st$type == "s") {
      for (l in alphabet$w) {
        lk <- .fw_key(l)
        st$transitions[lk] <- if (.sm_has(M, lk)) lk else key
      }
      assign(key, st, envir = M$states)
    } else if (st$type == "a") {
      for (l in alphabet$w) st$transitions[.fw_key(l)] <- next_state(st$name, l)
      assign(key, st, envir = M$states)
    }
  }
  if (length(M$automaton_alphabet) != length(alphabet$w)) {
    for (key in accept_states) {
      st <- get(key, envir = M$states)
      for (l in M$automaton_alphabet) st$transitions[.fw_key(l)] <- next_state(st$name, l)
      assign(key, st, envir = M$states)
    }
  }
  invisible(NULL)
}

#' @rdname RewritingSystem
#' @export
rws_reduce_using_automaton <- function(R, word) {
  if (length(R$new_rules$k)) {
    .rws_add_to_automaton(R, "new_rules")
    R$new_rules <- .od()
  }
  M <- R$automaton
  flag <- TRUE
  while (flag) {
    flag <- FALSE
    current <- get("start", envir = M$states)
    letters <- word_letter_form_elm(word)
    for (i in seq_along(letters)) {
      lk <- .fw_key(letters[[i]])
      nk <- current$transitions[lk]
      if (is.na(nk)) .fg_stop("KeyError", .fw_str(letters[[i]]))
      nxt <- get(nk, envir = M$states)
      if (nxt$type == "d") {
        word <- substituted_word(word, i - 1L - .fw_len(nxt$name) + 1L, i, nxt$rh_rule)
        flag <- TRUE
        break
      }
      current <- nxt
    }
  }
  word
}

# ---------------------------------------------------------------------------
# Permutation groups: R/combinatorics-permgroups.R (SymPy's perm_groups.py)
# does the work; these are the two perm group methods that need this file
# (polycyclic_group and presentation) and a few small adapters.

# group orders as Python ints (R integers when they fit)
.ppg_order <- function(G) {
  o <- group_order(G)
  if (is.numeric(o) && o <= .Machine$integer.max) as.integer(o) else o
}
.ppg_slp_get <- function(G, g) G$strong_gens_slp[[.pkey(g)]]
.fpp_str <- function(p) .perm_str(p)
.fpp_repr <- function(p) .perm_str(p, repr = TRUE)
.fpp_mul <- function(p, q) p * q
.fpp_inv <- function(p) p^-1L
.fpp_pow <- function(p, k) p^as.integer(k)
.fpp_order <- function(p) as.integer(perm_order(p))
.fpp_eq <- function(p, q) inherits(p, "Permutation") && inherits(q, "Permutation") && identical(p$array, q$array)

.ppg_polycyclic_group <- function(G) {
  if (!is_solvable(G)) .fg_stop("ValueError", "The group must be solvable")
  der <- derived_series(G)
  pc_series <- list(der[[length(der)]])
  pc_sequence <- list()
  relative_order <- integer()
  der <- rev(der)
  for (i in seq_len(length(der) - 1L)) {
    H <- der[[i]]
    for (g in der[[i + 1L]]$gens) {
      if (!group_contains(H, g)) {
        H <- PermutationGroup(c(list(g), H$gens))
        pc_series <- c(list(H), pc_series)
        pc_sequence <- c(list(g), pc_sequence)
        relative_order <- c(as.integer(.ppg_order(pc_series[[1L]]) %/% .ppg_order(pc_series[[2L]])), relative_order)
      }
    }
  }
  PolycyclicGroup(pc_sequence, pc_series, relative_order)
}

# presentation (perm_groups.py, Chapter 6.1 of the Handbook)
.ppg_presentation <- function(G, eliminate_gens = TRUE) {
  if (!is.null(G$fp_presentation)) return(G$fp_presentation)
  gens <- G$gens
  len_g <- length(gens)
  if (len_g == 1L) {
    ord <- as.integer(perm_order(gens[[1L]]))
    if (ord == 1L) return(FreeGroup(character()))
    fx <- free_group("x")
    return(FpGroup(fx[[1L]], list(.fw_pow(fx[[2L]], ord))))
  }
  half_gens <- if (.ppg_order(G) > 20) gens[seq_len((len_g + 1L) %/% 2L)] else list()
  H <- PermutationGroup(half_gens)
  H_p <- .ppg_presentation(H)
  len_h <- length(H_p$generators)
  C <- coset_table(G, H)
  n <- length(C)
  Fg <- FreeGroup(paste0("x_", seq_len(len_g) - 1L))
  images <- lapply(seq_len(len_h), function(i) Fg$generators[[i]])
  R <- homomorphism(H_p, Fg, H_p$generators, images, check = FALSE)
  rels <- .hk_apply(R, if (.is_fpg(H_p)) H_p$relators else list())
  G_p <- FpGroup(Fg, rels)
  T <- homomorphism(G_p, G, G_p$generators, gens)
  C_p <- CosetTable(G_p, list())
  C_p$table <- lapply(seq_len(n), function(i) rep(NA_integer_, 2L * len_g))
  transversal <- vector("list", n)
  transversal[[1L]] <- G_p$identity
  for (i in seq_len(2L * len_h)) C_p$table[[1L]][i] <- 0L
  gamma <- 1L
  done <- FALSE
  for (alpha in seq_len(n) - 1L) {
    for (x in seq_len(2L * len_g) - 1L) {
      beta <- C[[alpha + 1L]][x + 1L]
      if (beta == gamma) {
        gen <- .fw_pow(G_p$generators[[x %/% 2L + 1L]], (-1L)^(x %% 2L))
        transversal[[beta + 1L]] <- .fw_mul(transversal[[alpha + 1L]], gen)
        C_p$table[[alpha + 1L]][x + 1L] <- beta
        C_p$table[[beta + 1L]][x + (-1L)^(x %% 2L) + 1L] <- alpha
        gamma <- gamma + 1L
        if (gamma == n) {
          done <- TRUE
          break
        }
      }
    }
    if (done) break
  }
  C_p$p <- seq_len(n) - 1L
  beta <- 0L
  x <- 0L
  while (!ct_is_complete(C_p)) {
    while (identical(C_p$table[[beta + 1L]][x + 1L], C[[beta + 1L]][x + 1L])) {
      x <- (x + 1L) %% (2L * len_g)
      if (x == 0L) beta <- (beta + 1L) %% n
    }
    gen <- .fw_pow(G_p$generators[[x %/% 2L + 1L]], (-1L)^(x %% 2L))
    new_rel <- .fw_mul(.fw_mul(transversal[[beta + 1L]], gen), .fw_inv(transversal[[C[[beta + 1L]][x + 1L] + 1L]]))
    perm <- .hk_apply(T, new_rel)
    nxt <- .fw_prod(Fg, lapply(generator_product(H, perm, original = TRUE), function(s) .fw_inv(.hk_invert(T, s))))
    new_rel <- .fw_mul(new_rel, nxt)
    G_p <- FpGroup(G_p$free_group, .fw_unique(c(list(new_rel), as.list(G_p$relators))))
    ct_scan_and_fill(C_p, 0L, new_rel)
    C_p <- fp_coset_enumeration(G_p, list(), strategy = "coset_table", draft = C_p, max_cosets = n, incomplete = TRUE)
  }
  G$fp_presentation <- simplify_presentation(G_p, change_gens = eliminate_gens)
  G$fp_presentation
}

# ---------------------------------------------------------------------------
# Elements of any of the three group kinds

.is_word <- function(e) inherits(e, "FreeGroupElement")
.is_perm <- function(e) inherits(e, "Permutation")
.g_key <- function(e) if (.is_word(e)) paste0("w", .fw_key(e)) else paste0("p", .pkey(e))
.g_mul <- function(a, b) if (.is_word(a)) .fw_mul(a, b) else .fpp_mul(a, b)
.g_inv <- function(a) if (.is_word(a)) .fw_inv(a) else .fpp_inv(a)
.g_pow <- function(a, k) if (.is_word(a)) .fw_pow(a, k) else .fpp_pow(a, k)
.g_is_id <- function(a) if (.is_word(a)) .fw_is_id(a) else .is_identity(a)
.g_eq <- function(a, b) if (.is_word(a)) .fw_eq(a, b) else .fpp_eq(a, b)
.g_in <- function(e, lst) any(vapply(lst, function(x) .g_eq(e, x), NA))

.grp_identity <- function(G) if (.is_pg(G)) group_identity(G) else G$identity
.grp_generators <- function(G) if (.is_pg(G)) G$gens else as.list(G$generators)
.grp_prod <- function(G, elems) {
  if (.is_pg(G) || (length(elems) && .is_perm(elems[[1L]]))) return(perm_prod(elems))
  F <- if (.is_fpg(G)) G$free_group else G
  .fw_prod(F, elems)
}
# Python's `g in G`
.grp_contains <- function(G, g) {
  if (.is_pg(G)) {
    if (!.is_perm(g)) .fg_stop("TypeError", sprintf(
      "A PermutationGroup contains only Permutations as elements, not elements of type %s", .py_class(g)))
    return(group_contains(G, g))
  }
  fp_contains(G, g)
}
.py_class <- function(g) {
  if (is.list(g) && !is.object(g)) return("<class 'list'>")
  if (.is_word(g)) return("<class 'sympy.combinatorics.free_groups.FreeGroupElement'>")
  paste0("<class '", class(g)[1L], "'>")
}
.grp_order <- function(G) if (.is_pg(G)) .ppg_order(G) else fp_order(G)
.is_one <- function(v) {
  if (inherits(v, "sym")) return(identical(v$key, S_One$key))
  isTRUE(v == 1)
}

# ---------------------------------------------------------------------------
# Group homomorphisms (homomorphisms.py). Homomorphisms between permutation
# groups are the package's `GroupHomomorphism` objects
# (R/combinatorics-permgroups.R); a homomorphism with a finitely presented
# or free domain or codomain is an `FpGroupHomomorphism` made here. The
# `.hk_*` helpers work on both kinds.

.hk_mine <- function(T) inherits(T, "FpGroupHomomorphism")
.hk_apply <- function(T, g) if (.hk_mine(T)) .fhom_apply(T, g) else .hom_apply(T, g)
.hk_invert <- function(T, g) if (.hk_mine(T)) .fhom_invert(T, g) else .hom_invert(T, g)
.hk_kernel <- function(T) if (.hk_mine(T)) .fhom_kernel(T) else .hom_kernel(T)
.hk_image <- function(T) if (.hk_mine(T)) .fhom_image(T) else .hom_image(T)
.hk_pairs <- function(T) if (.hk_mine(T)) list(keys = T$images$key, values = T$images$v) else
  list(keys = T$from, values = T$to)
.hk_img <- function(T, g) {
  pr <- .hk_pairs(T)
  i <- match(.g_key(g), vapply(pr$keys, .g_key, ""))
  if (is.na(i)) .fg_stop("KeyError", .py_str(g))
  pr$values[[i]]
}
.hk_order <- function(K) if (inherits(K, "FpSubgroup")) fps_order(K) else .grp_order(K)

#' Group homomorphisms
#'
#' `homomorphism(domain, codomain, gens, images)` is the homomorphism that
#' maps the generators `gens` of `domain` to `images` (the other generators
#' go to the identity); domain and codomain are permutation groups,
#' finitely presented groups or free groups. With `check = TRUE` the images
#' are checked against the relators of the domain (the presentation of a
#' permutation group). Between permutation groups the result is the
#' package's `GroupHomomorphism`; otherwise it is an `FpGroupHomomorphism`.
#' `hom_is_injective()`, `hom_is_surjective()`, `hom_is_isomorphism()` and
#' `hom_is_trivial()` test properties of either kind, `hom_images()` gives
#' the images of the generators (SymPy's `images` dictionary),
#' `hom_compose(T, S)` is
#' `T` after `S`, `hom_factor()` gives a surjection and an injection whose
#' composite is `T`, `hom_restrict_to()` restricts to a subgroup and
#' `hom_invert_subgroup()` gives a preimage subgroup. `group_isomorphism()`
#' searches an isomorphism and `is_isomorphic()` tests for one. Applying,
#' inverting, kernels and images are [hom_apply()], [hom_invert()],
#' [hom_kernel()] and [hom_image()].
#'
#' @param domain,codomain,G,H Groups.
#' @param gens A list of generators of the domain.
#' @param images A list of elements of the codomain; for
#'   `GroupHomomorphism()` an ordered dictionary as built by
#'   `homomorphism()`.
#' @param check Check that the images define a homomorphism.
#' @param T,other Homomorphisms.
#' @param isomorphism Also return the isomorphism.
#' @return A homomorphism (an environment), a group or a logical.
#' @examples
#' g <- free_group("a, b")
#' a <- g[[2]]; b <- g[[3]]
#' G <- FpGroup(g[[1]], list(a^3, b^3, (a * b)^2))
#' A <- AlternatingGroup(4)
#' T <- homomorphism(G, A, list(a, b), list(Permutation(0, 1, 2, size = 4), Permutation(1, 2, 3)))
#' hom_is_isomorphism(T)
#' is_isomorphic(G, A)
#' @export
GroupHomomorphism <- function(domain, codomain, images) {
  if (.is_pg(domain) && .is_pg(codomain)) return(.hom_new(domain, codomain, images$key, images$v))
  T <- new.env(parent = emptyenv())
  T$domain <- domain
  T$codomain <- codomain
  T$images <- images
  T$inverses <- NULL
  T$kernel <- NULL
  T$image <- NULL
  T$surjective_cert <- NULL
  T$surjective_cert_computed <- FALSE
  T$factorization <- NULL
  class(T) <- "FpGroupHomomorphism"
  T
}

#' @export
format.FpGroupHomomorphism <- function(x, ...) "GroupHomomorphism"

#' @export
print.FpGroupHomomorphism <- function(x, ...) {
  cat("GroupHomomorphism from", .py_str(x$domain), "to", .py_str(x$codomain), "\n")
  invisible(x)
}

.fhom_img <- function(T, g) {
  v <- .od_get(T$images, .g_key(g))
  if (is.null(v)) .fg_stop("KeyError", .py_str(g))
  v
}
.fhom_has <- function(T, g) .od_has(T$images, .g_key(g))

.fhom_invs <- function(T) {
  image <- .fhom_image(T)
  inverses <- .od()
  for (i in seq_along(T$images$k)) {
    k <- T$images$key[[i]]
    v <- T$images$v[[i]]
    if (!(.od_has(inverses, .g_key(v)) || .g_is_id(v))) inverses <- .od_set(inverses, .g_key(v), v, k)
  }
  perm_cod <- .is_pg(T$codomain)
  gens <- if (perm_cod) strong_gens(image) else image$generators
  for (g in gens) {
    if (.od_has(inverses, .g_key(g)) || .g_is_id(g)) next
    if (!perm_cod) .fg_stop("TypeError", "unsupported operand type(s) for ** or pow(): 'tuple' and 'int'")
    parts <- rev(.ppg_slp_get(image, g))
    w <- .grp_prod(T$domain, lapply(parts, function(s) {
      v <- .od_get(inverses, .g_key(s))
      if (!is.null(v)) v else .g_inv(.od_get(inverses, .g_key(.g_inv(s))))
    }))
    inverses <- .od_set(inverses, .g_key(g), g, w)
  }
  inverses
}

.fhom_invert <- function(T, g) {
  if (.is_perm(g) || .is_word(g)) {
    if (.is_fpg(T$codomain)) g <- fp_reduce(T$codomain, g)
    if (is.null(T$inverses)) T$inverses <- .fhom_invs(T)
    image <- .fhom_image(T)
    w <- .grp_identity(T$domain)
    if (.is_pg(T$codomain)) {
      if (!.grp_contains(image, g)) .fg_stop("ValueError", "Given element is not in the image of the homomorphism.")
      gens <- rev(generator_product(image, g))
      parts <- lapply(Filter(function(s) !.is_identity(s), gens), function(s) {
        v <- .od_get(T$inverses, .g_key(s))
        if (!is.null(v)) v else .g_inv(.od_get(T$inverses, .g_key(.g_inv(s))))
      })
      return(.grp_prod(T$domain, c(list(w), parts)))
    }
    current_coset <- 0L
    C <- modified_coset_enumeration_r(T$codomain, T$inverses$key)
    preimages <- T$inverses$v
    factors <- list()
    for (i in seq_along(g$s)) {
      s <- g$s[i]
      ex <- g$e[i]
      j <- .ct_col(if (ex < 0L) -s else s)
      for (r in seq_len(abs(ex))) {
        wi <- C$P[[current_coset + 1L]][[j]]
        if (!is.null(wi)) {
          for (t in seq_along(wi$s)) factors[[length(factors) + 1L]] <- .g_pow(preimages[[wi$s[t]]], wi$e[t])
        }
        current_coset <- C$table[[current_coset + 1L]][j]
      }
    }
    if (!identical(current_coset, 0L)) .fg_stop("ValueError", "Given element is not in the image of the homomorphism.")
    return(.grp_prod(T$domain, factors))
  }
  if (is.list(g) && !is.object(g)) return(lapply(g, function(e) .fhom_invert(T, e)))
  NULL
}

.fhom_kernel <- function(T) {
  if (is.null(T$kernel)) T$kernel <- .fhom_compute_kernel(T)
  T$kernel
}

.fhom_compute_kernel <- function(T) {
  G <- T$domain
  H <- T$codomain
  Ginf <- .is_oo(.grp_order(G))
  if (Ginf && (.is_fpg(H) || .is_fg(H) || .is_pg(H))) {
    perm_cod <- .is_pg(H)
    if (.is_fg(H)) H <- FpGroup(H, list())
    preimages <- .fhom_surjective_cert(T)
    if (!is.null(preimages)) {
      if (perm_cod) {
        perm_gens <- H$gens
        fpg <- .ppg_presentation(H, eliminate_gens = FALSE)
        if (.is_fg(fpg)) return(G)
        sym_pre <- lapply(seq_len(min(length(fpg$generators), length(perm_gens))),
                          function(i) .od_get(preimages, .g_key(perm_gens[[i]])))
        names(sym_pre) <- vapply(fpg$generators, function(g) g$group$symbols[g$s[1L]], "")[seq_along(sym_pre)]
      } else {
        sym_pre <- lapply(H$generators, function(g) .od_get(preimages, .g_key(g)))
        names(sym_pre) <- vapply(H$generators, function(g) g$group$symbols[g$s[1L]], "")
      }
      lift <- function(word) {
        if (perm_cod && .is_perm(word)) {
          factors <- rev(generator_product(H, word, original = TRUE))
          return(.grp_prod(G, lapply(Filter(function(s) !.is_identity(s), factors), function(s) {
            v <- .od_get(preimages, .g_key(s))
            if (!is.null(v)) v else .g_inv(.od_get(preimages, .g_key(.g_inv(s))))
          })))
        }
        syms <- word$group$symbols[word$s]
        .grp_prod(G, lapply(seq_along(syms), function(i) .g_pow(sym_pre[[syms[i]]], word$e[i])))
      }
      kernel_gens <- list()
      relators <- if (perm_cod) fpg$relators else H$relators
      for (relator in relators) {
        word <- lift(relator)
        if (!.g_is_id(word)) kernel_gens[[length(kernel_gens) + 1L]] <- word
      }
      for (generator in .grp_generators(G)) {
        word <- .g_mul(generator, .g_inv(lift(.fhom_img(T, generator))))
        if (!.g_is_id(word)) kernel_gens[[length(kernel_gens) + 1L]] <- word
      }
      return(FpSubgroup(G, kernel_gens, normal = TRUE))
    }
    return(.hk_kernel(hom_factor(T)[[1L]]))
  }
  if (Ginf) .fg_stop("NotImplementedError", "Kernel computation is not implemented for this group type")
  if (.is_fpg(G)) {
    r <- fp_to_perm_group(G)
    P <- r[[1L]]
    TT <- r[[2L]]
    perm_images <- .od()
    for (g in P$gens) perm_images <- .od_set(perm_images, .g_key(g), g, .fhom_apply(T, .fhom_invert(TT, g)))
    perm_kernel <- .hk_kernel(GroupHomomorphism(P, H, perm_images))
    return(FpSubgroup(G, .fhom_invert(TT, perm_kernel$gens), normal = TRUE))
  }
  if (.is_pg(G)) return(subgroup_search(G, function(g) .g_is_id(.fhom_apply(T, g))))
  .fg_stop("NotImplementedError", "Kernel computation is not implemented for this group type")
}

.fhom_image <- function(T) {
  if (is.null(T$image)) {
    values <- T$images$v[!duplicated(vapply(T$images$v, .g_key, ""))]
    T$image <- if (.is_pg(T$codomain)) subgroup(T$codomain, values) else FpSubgroup(T$codomain, values)
  }
  T$image
}

.fhom_apply <- function(T, g) {
  dom <- T$domain
  inside <- if (.is_pg(dom)) .grp_contains(dom, g) else fp_contains(dom, g)
  if (!inside) {
    if (is.list(g) && !is.object(g) || inherits(g, "py_tuple")) return(lapply(g, function(e) .fhom_apply(T, e)))
    .fg_stop("ValueError", "The supplied element does not belong to the domain")
  }
  if (.g_is_id(g)) return(.grp_identity(T$codomain))
  value <- .grp_identity(T$codomain)
  if (.is_pg(dom)) {
    for (h in generator_product(dom, g, original = TRUE)) {
      if (.fhom_has(T, h)) value <- .g_mul(.fhom_img(T, h), value)
      else value <- .g_mul(.g_inv(.fhom_img(T, .g_inv(h))), value)
    }
  } else {
    for (i in seq_along(g$s)) {
      gen <- .fw(g$group, g$s[i], 1L)
      value <- .g_mul(value, .g_pow(.fhom_img(T, gen), g$e[i]))
    }
  }
  value
}

#' @rdname GroupHomomorphism
#' @export
hom_images <- function(T) {
  pr <- .hk_pairs(T)
  d <- .od()
  for (i in seq_along(pr$keys)) d <- .od_set(d, .g_key(pr$keys[[i]]), pr$keys[[i]], pr$values[[i]])
  d
}

#' @rdname GroupHomomorphism
#' @export
hom_is_injective <- function(T) .is_one(.hk_order(.hk_kernel(T)))

.fhom_surjective_cert <- function(T) {
  if (!isTRUE(T$surjective_cert_computed)) {
    pre <- tryCatch({
      d <- .od()
      for (g in .grp_generators(T$codomain)) d <- .od_set(d, .g_key(g), g, .hk_invert(T, g))
      d
    }, ValueError = function(e) NULL)
    T$surjective_cert <- pre
    T$surjective_cert_computed <- TRUE
  }
  T$surjective_cert
}

#' @rdname GroupHomomorphism
#' @export
hom_is_surjective <- function(T) !is.null(.fhom_surjective_cert(T))

#' @rdname GroupHomomorphism
#' @export
hom_is_isomorphism <- function(T) hom_is_injective(T) && hom_is_surjective(T)

#' @rdname GroupHomomorphism
#' @export
hom_is_trivial <- function(T) .is_one(.hk_order(.hk_image(T)))

#' @rdname GroupHomomorphism
#' @export
hom_factor <- function(T) {
  if (!is.null(T$factorization)) return(T$factorization)
  dgens <- .grp_generators(T$domain)
  if (.is_pg(T$codomain)) {
    image <- .hk_image(T)
    surj <- homomorphism(T$domain, image, dgens, lapply(dgens, function(g) .hk_img(T, g)), check = FALSE)
    inj <- homomorphism(image, T$codomain, image$gens, image$gens)
    T$factorization <- list(surj, inj)
    return(T$factorization)
  }
  if (!(.is_fpg(T$codomain) || .is_fg(T$codomain)))
    .fg_stop("NotImplementedError", "Factorization is implemented only for permutation, finitely presented, and free codomains")
  codomain <- if (.is_fpg(T$codomain)) T$codomain else FpGroup(T$codomain, list())
  image_gens <- T$images$v[!duplicated(vapply(T$images$v, .g_key, ""))]
  r <- fp_subgroup(codomain, image_gens, homomorphism = TRUE)
  image_group <- r[[1L]]
  image_inclusion <- r[[2L]]
  surj_images <- lapply(dgens, function(g) .hk_invert(image_inclusion, .fhom_img(T, g)))
  surj <- homomorphism(T$domain, image_group, dgens, surj_images, check = FALSE)
  inj <- homomorphism(image_group, T$codomain, image_group$generators,
                      .hk_apply(image_inclusion, as.list(image_group$generators)), check = FALSE)
  T$factorization <- list(surj, inj)
  T$factorization
}

#' @rdname GroupHomomorphism
#' @export
hom_compose <- function(T, other) {
  im <- .hk_image(other)
  if (!.is_pg(im)) .fg_stop("AttributeError", "'FpSubgroup' object has no attribute 'is_subgroup'")
  if (!is_subgroup(im, T$domain))
    .fg_stop("ValueError", "The image of `other` must be a subgroup of the domain of `self`")
  images <- .od()
  for (g in .hk_pairs(other)$keys) images <- .od_set(images, .g_key(g), g, .hk_apply(T, .hk_apply(other, g)))
  GroupHomomorphism(other$domain, T$codomain, images)
}

#' @rdname GroupHomomorphism
#' @export
hom_restrict_to <- function(T, H) {
  if (!.is_pg(H) || !is_subgroup(H, T$domain)) .fg_stop("ValueError", "Given H is not a subgroup of the domain")
  images <- .od()
  for (g in H$gens) images <- .od_set(images, .g_key(g), g, .hk_apply(T, g))
  GroupHomomorphism(H, T$codomain, images)
}

#' @rdname GroupHomomorphism
#' @export
hom_invert_subgroup <- function(T, H) {
  if (!is_subgroup(H, .hk_image(T))) .fg_stop("ValueError", "Given H is not a subgroup of the image")
  gens <- list()
  P <- PermutationGroup(list(group_identity(.hk_image(T))))
  for (h in H$gens) {
    h_i <- .hk_invert(T, h)
    if (!group_contains(P, h_i)) {
      gens[[length(gens) + 1L]] <- h_i
      P <- PermutationGroup(gens)
    }
    for (k in .hk_kernel(T)$gens) {
      if (!group_contains(P, k * h_i)) {
        gens[[length(gens) + 1L]] <- k * h_i
        P <- PermutationGroup(gens)
      }
    }
  }
  P
}

#' @rdname GroupHomomorphism
#' @export
homomorphism <- function(domain, codomain, gens, images = list(), check = TRUE) {
  isgrp <- function(G) .is_pg(G) || .is_fpg(G) || .is_fg(G)
  if (!isgrp(domain)) .fg_stop("TypeError", "The domain must be a group")
  if (!isgrp(codomain)) .fg_stop("TypeError", "The codomain must be a group")
  generators <- .grp_generators(domain)
  gens <- as.list(gens)
  images <- as.list(images)
  if (!all(vapply(gens, function(g) .g_in(g, generators), NA)))
    .fg_stop("ValueError", "The supplied generators must be a subset of the domain's generators")
  if (!all(vapply(images, function(g) .grp_contains(codomain, g), NA)))
    .fg_stop("ValueError", "The images must be elements of the codomain")
  if (length(images) && length(images) != length(gens))
    .fg_stop("ValueError", "The number of images must be equal to the number of generators")
  images <- c(images, rep(list(.grp_identity(codomain)), length(generators) - length(images)))
  gens <- c(gens, Filter(function(g) !.g_in(g, gens), generators))
  d <- .od()
  for (i in seq_along(gens)) d <- .od_set(d, .g_key(gens[[i]]), gens[[i]], images[[i]])
  if (check && !.check_homomorphism(domain, codomain, d))
    .fg_stop("ValueError", "The given images do not define a homomorphism")
  GroupHomomorphism(domain, codomain, d)
}

.check_homomorphism <- function(domain, codomain, images) {
  pres <- if (.is_pg(domain)) .ppg_presentation(domain) else domain
  rels <- pres$relators
  gens <- pres$generators
  dgens <- .grp_generators(domain)
  syms <- vapply(gens, function(g) g$group$symbols[g$s[1L]], "")
  m <- min(length(syms), length(dgens))
  identity <- .grp_identity(codomain)
  image_of <- function(r) {
    rs <- r$group$symbols[r$s]
    .grp_prod(codomain, lapply(seq_along(rs), function(i) {
      k <- match(rs[i], syms[seq_len(m)])
      if (is.na(k)) .fg_stop("KeyError", rs[i])
      .g_pow(.od_get(images, .g_key(dgens[[k]])), r$e[i])
    }))
  }
  for (r in rels) {
    if (.is_fpg(codomain)) {
      s <- fp_equals(codomain, image_of(r), identity)
      if (is.null(s)) {
        fp_make_confluent(codomain)
        s <- fp_equals(codomain, image_of(r), identity)
        if (is.null(s)) .fg_stop("RuntimeError", paste0(
          "Can't determine if the images define a homomorphism. Try increasing the maximum number of rewriting rules ",
          "(group._rewriting_system.set_max(new_value); the current value is stored in group._rewriting_system.maxeqns)"))
      }
    } else {
      s <- .g_is_id(image_of(r))
    }
    if (!s) return(FALSE)
  }
  TRUE
}
#' @rdname GroupHomomorphism
#' @export
group_isomorphism <- function(G, H, isomorphism = TRUE) {
  if (!(.is_pg(G) || .is_fpg(G))) .fg_stop("TypeError", "The group must be a PermutationGroup or an FpGroup")
  if (!(.is_pg(H) || .is_fpg(H))) .fg_stop("TypeError", "The group must be a PermutationGroup or an FpGroup")
  if (.is_fpg(G) && .is_fpg(H)) {
    G <- simplify_presentation(G)
    H <- simplify_presentation(H)
    # SymPy compares (G.relators).sort() == (H.relators).sort(), which is
    # None == None, after sorting both relator lists in place
    G$relators <- .fw_sort(G$relators)
    H$relators <- .fw_sort(H$relators)
    same <- length(G$generators) == length(H$generators) &&
      all(vapply(seq_along(G$generators), function(i) .fw_eq(G$generators[[i]], H$generators[[i]]), NA))
    if (same) {
      if (!isomorphism) return(TRUE)
      return(.py_tuple(list(TRUE, homomorphism(G, H, G$generators, H$generators))))
    }
  }
  .H <- H
  g_order <- .grp_order(G)
  h_order <- .grp_order(H)
  if (.is_oo(g_order)) .fg_stop("NotImplementedError", "Isomorphism methods are not implemented for infinite groups.")
  if (.is_fpg(H)) {
    if (.is_oo(h_order)) .fg_stop("NotImplementedError", "Isomorphism methods are not implemented for infinite groups.")
    r <- fp_to_perm_group(H)
    .H <- r[[1L]]
    h_isomorphism <- r[[2L]]
  }
  num <- function(v) if (inherits(v, "sym")) as.numeric(as.character(v$value)) else as.numeric(v)
  if (num(g_order) != num(h_order) || fp_is_abelian(G) != fp_is_abelian(H)) {
    if (!isomorphism) return(FALSE)
    return(.py_tuple(list(FALSE, NULL)))
  }
  if (!isomorphism) {
    n <- num(g_order)
    if (igcd(n, as.numeric(as.character(totient(n)$value))) == 1) return(TRUE)
  }
  gens <- .grp_generators(G)
  try_images <- function(images) {
    d <- .od()
    for (i in seq_along(gens)) d <- .od_set(d, .g_key(gens[[i]]), gens[[i]], images[[i]])
    if (.check_homomorphism(G, .H, d)) {
      if (.is_fpg(H)) images <- .hk_invert(h_isomorphism, images)
      T <- homomorphism(G, H, .grp_generators(G), images, check = FALSE)
      if (hom_is_isomorphism(T)) return(T)
    }
    NULL
  }
  if (.is_pg(G)) {
    ords <- vapply(gens, .fpp_order, 0L)
    keys <- unique(ords)
    diz_gen <- lapply(keys, function(o) gens[ords == o])
    els <- group_generate(.H)
    eord <- vapply(els, .fpp_order, 0L)
    diz_imm <- lapply(keys, function(o) els[eord == o])
    lista <- lapply(seq_along(keys), function(i) .py_permutations(diz_imm[[i]], length(diz_gen[[i]])))
    for (m in .py_product(lista)) {
      images <- unlist(m, recursive = FALSE)
      T <- try_images(images)
      if (!is.null(T)) {
        if (!isomorphism) return(TRUE)
        return(.py_tuple(list(TRUE, T)))
      }
    }
  } else {
    els <- group_generate(.H)
    for (subset in .py_permutations(els, length(gens))) {
      images <- c(subset, rep(list(group_identity(.H)), length(gens) - length(subset)))
      T <- try_images(images)
      if (!is.null(T)) {
        if (!isomorphism) return(TRUE)
        return(.py_tuple(list(TRUE, T)))
      }
    }
  }
  if (!isomorphism) return(FALSE)
  .py_tuple(list(FALSE, NULL))
}

# itertools.permutations(x, k) and itertools.product(*lists), as lists
.py_permutations <- function(x, k) {
  n <- length(x)
  if (k > n) return(list())
  out <- list()
  rec <- function(chosen) {
    if (length(chosen) == k) {
      out[[length(out) + 1L]] <<- x[chosen]
      return(invisible())
    }
    for (i in seq_len(n)) if (!(i %in% chosen)) rec(c(chosen, i))
  }
  rec(integer())
  out
}
.py_product <- function(lists) {
  out <- list(list())
  for (l in lists) {
    nxt <- list()
    for (a in out) for (b in l) nxt[[length(nxt) + 1L]] <- c(a, list(b))
    out <- nxt
  }
  out
}

#' @rdname GroupHomomorphism
#' @export
is_isomorphic <- function(G, H) group_isomorphism(G, H, isomorphism = FALSE)

# ---------------------------------------------------------------------------
# FpGroup methods that go through an isomorphic permutation group

#' Permutation group methods of finitely presented groups
#'
#' SymPy answers these questions about a finite `FpGroup` through an
#' isomorphic permutation group from coset enumeration:
#' `fp_to_perm_group()` returns that group and the isomorphism,
#' `fp_is_abelian()`, `fp_is_cyclic()`, `fp_is_solvable()` and
#' `fp_is_nilpotent()` test properties, `fp_center()`,
#' `fp_derived_subgroup()`, `fp_derived_series()`,
#' `fp_lower_central_series()`, `fp_centralizer()` and
#' `fp_normal_closure()` give generators of subgroups as words, and
#' `fp_elements()`, `fp_abelian_invariants()` and `fp_composition_series()`
#' are as in SymPy. `fp_is_abelian()`, `fp_center()` and `fp_elements()`
#' also accept free groups.
#'
#' @param G An `FpGroup` (or a free group where noted).
#' @param other A list of words.
#' @return A list, a logical or a list of permutation groups.
#' @examples
#' g <- free_group("x, y")
#' x <- g[[2]]; y <- g[[3]]
#' G <- FpGroup(g[[1]], list(x^2, y^8, x * y * x^-1 * y))
#' fp_str(fp_center(G))
#' fp_is_abelian(G)
#' @export
fp_to_perm_group <- function(G) {
  if (.is_oo(fp_order(G))) .fg_stop("NotImplementedError", "Permutation presentation of infinite groups is not implemented")
  if (!is.null(G$perm_iso)) {
    T <- G$perm_iso
    P <- .hk_image(T)
  } else {
    C <- fp_coset_table(G, list())
    gens <- as.list(G$generators)
    images <- lapply(seq_along(gens), function(i) .perm_af(vapply(C, function(row) row[2L * i - 1L], 0L)))
    P <- PermutationGroup(images)
    T <- homomorphism(G, P, gens, images, check = FALSE)
    G$perm_iso <- T
  }
  list(P, T)
}

.fpg_perm_group_list <- function(G, f, ...) {
  r <- fp_to_perm_group(G)
  res <- f(r[[1L]], ...)
  if (.is_pg(res)) return(.hk_invert(r[[2L]], res$gens))
  lapply(res, function(H) .hk_invert(r[[2L]], H$gens))
}

#' @rdname fp_to_perm_group
#' @export
fp_derived_series <- function(G) .fpg_perm_group_list(G, derived_series)

#' @rdname fp_to_perm_group
#' @export
fp_lower_central_series <- function(G) .fpg_perm_group_list(G, lower_central_series)

#' @rdname fp_to_perm_group
#' @export
fp_center <- function(G) {
  if (.is_fg(G)) return(.py_set(list(G$identity)))
  .fpg_perm_group_list(G, group_center)
}

#' @rdname fp_to_perm_group
#' @export
fp_derived_subgroup <- function(G) .fpg_perm_group_list(G, derived_subgroup)

#' @rdname fp_to_perm_group
#' @export
fp_centralizer <- function(G, other) {
  T <- fp_to_perm_group(G)[[2L]]
  .fpg_perm_group_list(G, centralizer, .hk_apply(T, other))
}

#' @rdname fp_to_perm_group
#' @export
fp_normal_closure <- function(G, other) {
  T <- fp_to_perm_group(G)[[2L]]
  .fpg_perm_group_list(G, normal_closure, .hk_apply(T, other))
}

#' @rdname fp_to_perm_group
#' @export
fp_is_abelian <- function(G) {
  if (.is_fg(G)) return(G$rank %in% c(0L, 1L))
  if (.is_pg(G)) return(is_abelian(G))
  is_abelian(fp_to_perm_group(G)[[1L]])
}

#' @rdname fp_to_perm_group
#' @export
fp_is_nilpotent <- function(G) is_nilpotent(fp_to_perm_group(G)[[1L]])

#' @rdname fp_to_perm_group
#' @export
fp_is_solvable <- function(G) is_solvable(fp_to_perm_group(G)[[1L]])

#' @rdname fp_to_perm_group
#' @export
fp_elements <- function(G) {
  if (.is_fg(G)) {
    if (G$rank == 0L) return(.py_set(list(G$identity)))
    .fg_stop("ValueError", "Group contains infinitely many elements, hence cannot be represented")
  }
  r <- fp_to_perm_group(G)
  .hk_invert(r[[2L]], group_elements(r[[1L]]))
}

#' @rdname fp_to_perm_group
#' @export
fp_is_cyclic <- function(G) {
  if (length(G$generators) <= 1L) return(TRUE)
  r <- tryCatch(fp_to_perm_group(G), NotImplementedError = function(e)
    .fg_stop("NotImplementedError", "Check for infinite Cyclic group is not implemented"))
  is_cyclic(r[[1L]])
}

#' @rdname fp_to_perm_group
#' @export
fp_abelian_invariants <- function(G) {
  r <- tryCatch(fp_to_perm_group(G), NotImplementedError = function(e)
    .fg_stop("NotImplementedError", "abelian invariants is not implementedfor infinite group"))
  as.list(as.integer(abelian_invariants(r[[1L]])))
}

#' @rdname fp_to_perm_group
#' @export
fp_composition_series <- function(G) {
  r <- tryCatch(fp_to_perm_group(G), NotImplementedError = function(e)
    .fg_stop("NotImplementedError", "composition series is not implementedfor infinite group"))
  composition_series(r[[1L]])
}

#' Presentation of a permutation group
#'
#' Returns an `FpGroup` isomorphic to the permutation group `G`, built by
#' SymPy's algorithm (Handbook of Computational Group Theory, 6.1). A
#' cyclic group gives the presentation on one generator; the trivial
#' group gives the free group of rank 0.
#'
#' @param G A permutation group.
#' @param eliminate_gens Remove redundant generators.
#' @return An `FpGroup` or a free group.
#' @examples
#' P <- perm_group_presentation(DihedralGroup(4))
#' fp_order(P)
#' @export
perm_group_presentation <- function(G, eliminate_gens = TRUE) .ppg_presentation(G, eliminate_gens)

# ---------------------------------------------------------------------------
# Polyhedra (polyhedron.py)

.poly_key <- function(v) if (inherits(v, "sym")) sstr(v) else v
.poly_norm <- function(v) if (is.numeric(v) && length(v) == 1L && v == round(v)) as.integer(v) else v
# SymPy's minlex(f, directed=False): the least rotation of f or of its
# reverse, comparing elements by name (symbols) or value (integers)
.poly_minlex <- function(f) {
  f <- as.list(f)
  n <- length(f)
  if (n < 2L) return(f)
  keys <- vapply(f, function(v) as.character(.poly_key(v)), "")
  numeric_keys <- all(vapply(f, is.numeric, NA))
  kv <- if (numeric_keys) as.numeric(keys) else keys
  less <- function(a, b) {
    for (i in seq_along(a)) if (a[i] != b[i]) return(a[i] < b[i])
    FALSE
  }
  best_rot <- function(idx) {
    best <- idx
    for (r in seq_len(n - 1L)) {
      cand <- idx[c((r + 1L):n, seq_len(r))]
      if (less(kv[cand], kv[best])) best <- cand
    }
    best
  }
  b1 <- best_rot(seq_len(n))
  b2 <- best_rot(rev(seq_len(n)))
  f[if (less(kv[b2], kv[b1])) b2 else b1]
}

.poly_sort <- function(tuples) {
  if (!length(tuples)) return(list())
  numeric_keys <- all(vapply(tuples, function(t) all(vapply(t, is.numeric, NA)), NA))
  if (numeric_keys) {
    m <- max(vapply(tuples, length, 0L))
    key <- vapply(tuples, function(t) paste(c(sprintf("%04d", length(t)), sprintf("%012.0f", as.numeric(unlist(t)))), collapse = " "), "")
  } else {
    key <- vapply(tuples, function(t) paste(c(sprintf("%04d", length(t)), vapply(t, function(v) as.character(.poly_key(v)), "")), collapse = "\r"), "")
  }
  tuples <- tuples[!duplicated(key)]
  key <- key[!duplicated(key)]
  tuples[order(key, method = "radix")]
}

#' Polyhedra
#'
#' `Polyhedron(corners, faces, pgroup)` is SymPy's polyhedron: corners, faces
#' (stored in canonical rotation order, as a sorted set) and a group of
#' permutations of the corners. `poly_rotate()` applies a permutation (or
#' the generator of `pgroup` with the given zero based index) to the
#' corners in place, `poly_reset()` restores them and `poly_copy()` copies.
#' `poly_corners()`, `poly_faces()`, `poly_edges()`, `poly_pgroup()`,
#' `poly_size()`, `poly_array_form()` and `poly_cyclic_form()` are the
#' properties. `tetrahedron()`, `cube()`, `octahedron()`, `dodecahedron()`
#' and `icosahedron()` return new copies of the standard polyhedra, and
#' `polyhedron_faces(name)` gives their faces.
#'
#' @param corners A vector or list of corners.
#' @param faces A list of faces (vectors of corners).
#' @param pgroup A list of permutations.
#' @param P A polyhedron.
#' @param perm A permutation, or a zero based index into the group.
#' @param name One of `"tetrahedron"`, `"cube"`, `"octahedron"`,
#'   `"dodecahedron"`, `"icosahedron"`.
#' @return A `Polyhedron` (an environment) or one of its properties.
#' @examples
#' C <- cube()
#' poly_rotate(C, 0)
#' fp_str(poly_corners(C))
#' fp_str(poly_edges(tetrahedron()))
#' @export
Polyhedron <- function(corners, faces = list(), pgroup = list()) {
  corners <- lapply(as.list(corners), .poly_norm)
  faces <- lapply(faces, function(f) .poly_minlex(lapply(as.list(f), .poly_norm)))
  P <- new.env(parent = emptyenv())
  P$args_corners <- corners
  P$args_faces <- faces
  P$args_pgroup <- as.list(pgroup)
  P$corners <- corners
  P$faces <- .poly_sort(faces)
  if (length(pgroup) && length(pgroup[[1L]]$array) != length(corners))
    .fg_stop("ValueError", "Permutation size unequal to number of corners.")
  P$pgroup <- PermutationGroup(if (length(pgroup)) pgroup else list(.perm_af(seq_along(corners) - 1L)))
  P$edges <- NULL
  class(P) <- "Polyhedron"
  P
}

.poly_str <- function(P) {
  tup <- function(x) .py_str(.py_tuple(x))
  paste0("Polyhedron(", tup(P$args_corners), ", ",
         .py_str(.py_tuple(lapply(P$args_faces, .py_tuple))), ", ",
         .py_str(.py_tuple(lapply(P$args_pgroup, function(p) structure(.fpp_str(p), class = "py_raw")))), ")")
}

#' @export
print.Polyhedron <- function(x, ...) {
  cat(.poly_str(x), "\n")
  invisible(x)
}

#' @rdname Polyhedron
#' @export
poly_corners <- function(P) .py_tuple(P$corners)

#' @rdname Polyhedron
#' @export
poly_faces <- function(P) .py_set(lapply(P$faces, .py_tuple))

#' @rdname Polyhedron
#' @export
poly_pgroup <- function(P) P$pgroup

#' @rdname Polyhedron
#' @export
poly_size <- function(P) length(P$corners)

#' @rdname Polyhedron
#' @export
poly_array_form <- function(P) {
  keys <- vapply(P$args_corners, function(v) as.character(.poly_key(v)), "")
  match(vapply(P$corners, function(v) as.character(.poly_key(v)), ""), keys) - 1L
}

#' @rdname Polyhedron
#' @export
poly_cyclic_form <- function(P) lapply(cyclic_form(.perm_af(poly_array_form(P))), as.integer)

#' @rdname Polyhedron
#' @export
poly_edges <- function(P) {
  if (is.null(P$edges)) {
    out <- list()
    for (face in P$faces) {
      n <- length(face)
      for (i in seq_len(n)) {
        prev <- if (i == 1L) n else i - 1L
        pair <- list(face[[i]], face[[prev]])
        if (as.character(.poly_key(pair[[2L]])) < as.character(.poly_key(pair[[1L]])) &&
            !(is.numeric(pair[[1L]]) && is.numeric(pair[[2L]]))) pair <- pair[2:1]
        if (is.numeric(pair[[1L]]) && is.numeric(pair[[2L]]) && pair[[2L]] < pair[[1L]]) pair <- pair[2:1]
        out[[length(out) + 1L]] <- pair
      }
    }
    P$edges <- .poly_sort(out)
  }
  .py_set(lapply(P$edges, .py_tuple))
}

#' @rdname Polyhedron
#' @export
poly_rotate <- function(P, perm) {
  if (!inherits(perm, "Permutation")) {
    perm <- P$pgroup$gens[[as.integer(perm) + 1L]]
  } else if (length(perm$array) != length(P$corners)) {
    .fg_stop("ValueError", "Polyhedron and Permutation sizes differ.")
  }
  a <- perm$array
  P$corners <- P$corners[a + 1L]
  invisible(P)
}

#' @rdname Polyhedron
#' @export
poly_reset <- function(P) {
  P$corners <- P$args_corners
  invisible(P)
}

#' @rdname Polyhedron
#' @export
poly_copy <- function(P) Polyhedron(P$args_corners, P$args_faces, P$args_pgroup)

.pgroup_of_double <- function(polyh, ordered_faces, pgroup) {
  n <- length(ordered_faces[[1L]])
  fkeys <- vapply(ordered_faces, function(f) paste(unlist(f), collapse = ","), "")
  flat <- unlist(ordered_faces)
  lapply(pgroup, function(p) {
    h <- poly_copy(polyh)
    poly_rotate(h, p)
    c <- unlist(h$corners)
    vals <- c[flat + 1L]
    reorder <- split(vals, rep(seq_along(ordered_faces), each = n))
    reorder <- lapply(reorder, function(f) paste(unlist(.poly_minlex(f)), collapse = ","))
    .perm_af(match(unlist(reorder), fkeys) - 1L)
  })
}

.poly_cache <- new.env(parent = emptyenv())

.poly_calcs <- function() {
  if (!is.null(.poly_cache$data)) return(.poly_cache$data)
  cyc <- function(...) Permutation(list(...))
  tetrahedron_faces <- list(c(0, 1, 2), c(0, 2, 3), c(0, 3, 1), c(1, 2, 3))
  t_pgroup <- list(cyc(c(1, 2, 3), 0), cyc(c(0, 1, 2), 3), cyc(c(0, 3, 2), 1), cyc(c(0, 3, 1), 2),
                   cyc(c(0, 1), c(2, 3)), cyc(c(0, 2), c(1, 3)), cyc(c(0, 3), c(1, 2)))
  tetra <- Polyhedron(0:3, tetrahedron_faces, t_pgroup)
  cube_faces <- list(c(0, 1, 2, 3), c(0, 1, 5, 4), c(1, 2, 6, 5), c(2, 3, 7, 6), c(0, 3, 7, 4), c(4, 5, 6, 7))
  c_pgroup <- lapply(list(
    c(1, 2, 3, 0, 5, 6, 7, 4), c(4, 0, 3, 7, 5, 1, 2, 6), c(4, 5, 1, 0, 7, 6, 2, 3),
    c(1, 0, 4, 5, 2, 3, 7, 6), c(6, 2, 1, 5, 7, 3, 0, 4), c(6, 7, 3, 2, 5, 4, 0, 1),
    c(3, 7, 4, 0, 2, 6, 5, 1), c(4, 7, 6, 5, 0, 3, 2, 1), c(6, 5, 4, 7, 2, 1, 0, 3),
    c(0, 3, 7, 4, 1, 2, 6, 5), c(5, 1, 0, 4, 6, 2, 3, 7), c(5, 6, 2, 1, 4, 7, 3, 0),
    c(7, 4, 0, 3, 6, 5, 1, 2)), Permutation)
  cub <- Polyhedron(0:7, cube_faces, c_pgroup)
  octahedron_faces <- list(c(0, 1, 2), c(0, 2, 3), c(0, 3, 4), c(0, 1, 4), c(1, 2, 5), c(2, 3, 5), c(3, 4, 5), c(1, 4, 5))
  octa <- Polyhedron(0:5, octahedron_faces, .pgroup_of_double(cub, cube_faces, c_pgroup))
  dodecahedron_faces <- list(c(0, 1, 2, 3, 4), c(0, 1, 6, 10, 5), c(1, 2, 7, 11, 6), c(2, 3, 8, 12, 7),
                             c(3, 4, 9, 13, 8), c(0, 4, 9, 14, 5), c(5, 10, 16, 15, 14), c(6, 10, 16, 17, 11),
                             c(7, 11, 17, 18, 12), c(8, 12, 18, 19, 13), c(9, 13, 19, 15, 14), c(15, 16, 17, 18, 19))
  f0 <- Permutation(c(1, 2, 3, 4, 0, 6, 7, 8, 9, 5, 11, 12, 13, 14, 10, 16, 17, 18, 19, 15))
  f1 <- Permutation(c(5, 0, 4, 9, 14, 10, 1, 3, 13, 15, 6, 2, 8, 19, 16, 17, 11, 7, 12, 18))
  string_to_perm <- function(s) {
    rv <- list(Permutation(0:19))
    p <- NULL
    for (si in strsplit(s, "")[[1L]]) {
      if (!(si %in% c("0", "1"))) {
        count <- as.integer(si) - 1L
      } else {
        count <- 1L
        p <- if (si == "0") f0 else f1
      }
      rv <- c(rv, rep(list(p), count))
    }
    perm_prod(rev(rv))
  }
  strs <- strsplit("0104 140 014 0410 010 1403 03104 04103 102 120 1304 01303 021302 03130 0412041 041204103 04120410 041204104 041204102 10 01 1402 0140 04102 0412 1204 1302 0130 03120", " ")[[1L]]
  dodeca_pgroup <- c(list(f0, f1), lapply(strs, string_to_perm))
  dodeca <- Polyhedron(0:19, dodecahedron_faces, dodeca_pgroup)
  icosahedron_faces <- list(c(0, 1, 2), c(0, 2, 3), c(0, 3, 4), c(0, 4, 5), c(0, 1, 5), c(1, 6, 7), c(1, 2, 7),
                            c(2, 7, 8), c(2, 3, 8), c(3, 8, 9), c(3, 4, 9), c(4, 9, 10), c(4, 5, 10), c(5, 6, 10),
                            c(1, 5, 6), c(6, 7, 11), c(7, 8, 11), c(8, 9, 11), c(9, 10, 11), c(6, 10, 11))
  icosa <- Polyhedron(0:11, icosahedron_faces, .pgroup_of_double(dodeca, dodecahedron_faces, dodeca_pgroup))
  .poly_cache$data <- list(tetrahedron = tetra, cube = cub, octahedron = octa, dodecahedron = dodeca,
                           icosahedron = icosa)
  .poly_cache$data
}

.poly_standard <- function(name) poly_copy(.poly_calcs()[[name]])

#' @rdname Polyhedron
#' @export
tetrahedron <- function() .poly_standard("tetrahedron")

#' @rdname Polyhedron
#' @export
cube <- function() .poly_standard("cube")

#' @rdname Polyhedron
#' @export
octahedron <- function() .poly_standard("octahedron")

#' @rdname Polyhedron
#' @export
dodecahedron <- function() .poly_standard("dodecahedron")

#' @rdname Polyhedron
#' @export
icosahedron <- function() .poly_standard("icosahedron")

#' @rdname Polyhedron
#' @export
polyhedron_faces <- function(name) lapply(.poly_calcs()[[name]]$faces, function(f) .py_tuple(as.list(as.integer(unlist(f)))))

# ---------------------------------------------------------------------------
# Polycyclic groups (pc_groups.py)

#' Polycyclic groups and collectors
#'
#' `perm_polycyclic_group(G)` returns the polycyclic group of a solvable
#' permutation group `G` (SymPy's `G.polycyclic_group()`), with its
#' polycyclic generating sequence `$pcgs`, series `$pc_series`, relative
#' orders `$relative_order` and `Collector` `$collector`. The collector's
#' `$pc_presentation` is an ordered dictionary of power and conjugate
#' relations; `pc_presentation()` lists it as pairs. `collected_word()`,
#' `exponent_vector()`, `pc_depth()`, `leading_exponent()`,
#' `induced_pcgs()`, `constructive_membership_test()`,
#' `minimal_uncollected_subword()` and `pc_relations()` are the methods of
#' SymPy's `Collector`.
#'
#' @param G A solvable permutation group.
#' @param pc_sequence,pcgs A list of permutations.
#' @param pc_series A list of permutation groups.
#' @param relative_order An integer vector.
#' @param collector A collector (made when `NULL`).
#' @param free_group_,pc_presentation_ Optional free group and presentation.
#' @param C A collector.
#' @param PC A polycyclic group.
#' @param word A word of the collector's free group.
#' @param element A permutation.
#' @param gens,ipcgs Lists of permutations.
#' @param g A permutation.
#' @return A `PolycyclicGroup`, a `Collector`, a word, a list or an integer.
#' @examples
#' PC <- perm_polycyclic_group(SymmetricGroup(4))
#' PC$relative_order
#' exponent_vector(PC$collector, SymmetricGroup(4)$gens[[1]])
#' @export
PolycyclicGroup <- function(pc_sequence, pc_series, relative_order, collector = NULL) {
  PC <- structure(list(pcgs = pc_sequence, pc_series = pc_series, relative_order = as.integer(relative_order),
                       collector = if (is.null(collector)) Collector(pc_sequence, pc_series, relative_order) else collector),
                  class = "PolycyclicGroup")
  PC
}

#' @rdname PolycyclicGroup
#' @export
perm_polycyclic_group <- function(G) .ppg_polycyclic_group(G)

#' @rdname PolycyclicGroup
#' @export
pc_is_prime_order <- function(PC) all(vapply(PC$relative_order, function(o) isTRUE(isprime(o)), NA))

#' @rdname PolycyclicGroup
#' @export
pc_length <- function(PC) length(PC$pcgs)

#' @rdname PolycyclicGroup
#' @export
Collector <- function(pcgs, pc_series, relative_order, free_group_ = NULL, pc_presentation_ = NULL) {
  C <- new.env(parent = emptyenv())
  C$pcgs <- pcgs
  C$pc_series <- pc_series
  C$relative_order <- as.integer(relative_order)
  C$free_group <- if (is.null(free_group_)) FreeGroup(paste0("x:", length(pcgs))) else free_group_
  C$pc_presentation <- .od()
  class(C) <- "Collector"
  C$pc_presentation <- .pc_relators(C)
  C
}

#' @export
print.Collector <- function(x, ...) {
  cat("Collector with", length(x$pcgs), "generators\n")
  invisible(x)
}

#' @rdname PolycyclicGroup
#' @export
pc_presentation <- function(C) lapply(seq_along(C$pc_presentation$k), function(i)
  .py_tuple(list(C$pc_presentation$key[[i]], C$pc_presentation$v[[i]])))

#' @rdname PolycyclicGroup
#' @export
minimal_uncollected_subword <- function(C, word) {
  if (.fw_is_id(word)) return(NULL)
  re <- C$relative_order
  s <- word$s
  e <- word$e
  for (i in seq_along(s)) {
    r <- re[s[i]]
    if (!is.na(r) && r != 0L && (e[i] < 0L || e[i] > r - 1L)) return(list(c(s[i], e[i])))
  }
  for (i in seq_len(length(s) - 1L)) {
    if (s[i] > s[i + 1L]) return(list(c(s[i], e[i]), c(s[i + 1L], if (e[i + 1L] > 0L) 1L else -1L)))
  }
  NULL
}

#' @rdname PolycyclicGroup
#' @export
pc_relations <- function(C) {
  power <- .od()
  conj <- .od()
  P <- C$pc_presentation
  for (i in seq_along(P$k)) {
    if (length(P$key[[i]]$s) == 1L) power <- .od_set(power, P$k[i], P$key[[i]], P$v[[i]])
    else conj <- .od_set(conj, P$k[i], P$key[[i]], P$v[[i]])
  }
  list(power, conj)
}

.pc_subword_index <- function(word, w) {
  lw <- .fw_len(w)
  for (i in seq_len(max(0L, .fw_len(word) - lw + 1L)) - 1L) {
    if (.fw_eq(subword(word, i, i + lw), w)) return(c(i, i + lw))
  }
  c(-1L, -1L)
}

.pc_map_relation <- function(C, w) {
  s1 <- w$s[1L]
  s2 <- w$s[2L]
  key <- .fw_new(C$free_group, c(s2, s1, s2), c(-1L, 1L, 1L))
  .od_get(C$pc_presentation, .fw_key(key))
}

.pc_word <- function(C, sw) .fw_new(C$free_group, vapply(sw, `[`, 0L, 1L), vapply(sw, `[`, 0L, 2L))

#' @rdname PolycyclicGroup
#' @export
collected_word <- function(C, word) {
  F <- C$free_group
  repeat {
    w <- minimal_uncollected_subword(C, word)
    if (is.null(w)) break
    lh <- .pc_subword_index(word, .pc_word(C, w))
    low <- lh[1L]
    high <- lh[2L]
    if (low == -1L) next
    s1 <- w[[1L]][1L]
    e1 <- w[[1L]][2L]
    if (length(w) == 1L) {
      re <- C$relative_order[s1]
      q <- e1 %/% re
      r <- e1 - q * re
      key <- .fw(F, s1, re)
      pres <- .od_get(C$pc_presentation, .fw_key(key))
      if (inherits(pres, "FreeGroupElement")) {
        word_ <- .fw_new(F, c(s1, pres$s[1L]), c(r, q * pres$e[1L]))
      } else {
        word_ <- if (r != 0L) .fw_new(F, s1, r) else NULL
      }
      word <- eliminate_word(word, .pc_word(C, w), word_)
    }
    if (length(w) == 2L && w[[2L]][2L] > 0L) {
      s2 <- .fw(F, w[[2L]][1L], 1L)
      word_ <- .pc_map_relation(C, .pc_word(C, w))
      word_ <- .fw_mul(s2, .pc_pow(word_, e1, F))
      word <- substituted_word(word, low, high, word_)
    } else if (length(w) == 2L && w[[2L]][2L] < 0L) {
      s2 <- .fw(F, w[[2L]][1L], 1L)
      word_ <- .pc_map_relation(C, .pc_word(C, w))
      word_ <- .fw_mul(.fw_inv(s2), .pc_pow(word_, e1, F))
      word <- substituted_word(word, low, high, word_)
    }
  }
  word
}

# a presentation value is a word or the empty tuple () for the identity
.pc_pow <- function(v, e, F) if (inherits(v, "FreeGroupElement")) .fw_pow(v, e) else
  .fg_stop("TypeError", "unsupported operand type(s) for ** or pow(): 'tuple' and 'int'")

.pc_relators <- function(C) {
  F <- C$free_group
  rel_order <- C$relative_order
  pc_relators <- .od()
  pcgs <- C$pcgs
  perm_to_free <- .od()
  for (i in seq_along(pcgs)) {
    gen <- pcgs[[i]]
    s <- F$generators[[i]]
    perm_to_free <- .od_set(perm_to_free, .g_key(.fpp_inv(gen)), .fpp_inv(gen), .fw_inv(s))
    perm_to_free <- .od_set(perm_to_free, .g_key(gen), gen, s)
  }
  p2f <- function(g) {
    v <- .od_get(perm_to_free, .g_key(g))
    if (is.null(v)) .fg_stop("KeyError", .fpp_str(g))
    v
  }
  pcgs <- rev(pcgs)
  series <- rev(C$pc_series)
  rel_order <- rev(rel_order)
  collected_gens <- list()
  for (i in seq_along(pcgs)) {
    gen <- pcgs[[i]]
    re <- rel_order[i]
    relation <- .fw_pow(p2f(gen), re)
    G <- series[[i]]
    l <- rev(generator_product(G, .fpp_pow(gen, re), original = TRUE))
    word <- .fw_prod(F, lapply(l, p2f))
    word <- collected_word(C, word)
    pc_relators <- .od_set(pc_relators, .fw_key(relation), relation, if (.fw_is_id(word)) .py_tuple(list()) else word)
    C$pc_presentation <- pc_relators
    collected_gens[[length(collected_gens) + 1L]] <- gen
    if (length(collected_gens) > 1L) {
      conj <- collected_gens[[length(collected_gens)]]
      conjugator <- p2f(conj)
      for (j in seq_len(length(collected_gens) - 1L)) {
        conjugated <- p2f(collected_gens[[j]])
        relation <- .fw_mul(.fw_mul(.fw_inv(conjugator), conjugated), conjugator)
        gens <- .fpp_mul(.fpp_mul(.fpp_inv(conj), collected_gens[[j]]), conj)
        l <- rev(generator_product(G, gens, original = TRUE))
        word <- .fw_prod(F, lapply(l, p2f))
        word <- collected_word(C, word)
        pc_relators <- .od_set(pc_relators, .fw_key(relation), relation, if (.fw_is_id(word)) .py_tuple(list()) else word)
        C$pc_presentation <- pc_relators
      }
    }
  }
  pc_relators
}

#' @rdname PolycyclicGroup
#' @export
exponent_vector <- function(C, element) {
  F <- C$free_group
  # ponytail: SymPy rebuilds this group (and its Schreier Sims data) on every
  # call; it depends only on the pcgs, so it is built once per collector
  G <- C$pcgs_group
  if (is.null(G)) {
    G <- PermutationGroup(list())
    for (g in C$pcgs) G <- PermutationGroup(c(list(g), G$gens))
    C$pcgs_group <- G
  }
  gens <- rev(generator_product(G, element, original = TRUE))
  perm_to_free <- .od()
  for (i in seq_along(C$pcgs)) {
    g <- C$pcgs[[i]]
    sym <- F$generators[[i]]
    perm_to_free <- .od_set(perm_to_free, .g_key(.fpp_inv(g)), .fpp_inv(g), .fw_inv(sym))
    perm_to_free <- .od_set(perm_to_free, .g_key(g), g, sym)
  }
  w <- .fw_prod(F, lapply(gens, function(g) {
    v <- .od_get(perm_to_free, .g_key(g))
    if (is.null(v)) .fg_stop("KeyError", .fpp_str(g))
    v
  }))
  word <- collected_word(C, w)
  ev <- integer(F$rank)
  ev[word$s] <- word$e
  ev
}

#' @rdname PolycyclicGroup
#' @export
pc_depth <- function(C, element) {
  ev <- exponent_vector(C, element)
  nz <- which(ev != 0L)
  if (length(nz)) nz[1L] else length(C$pcgs) + 1L
}

#' @rdname PolycyclicGroup
#' @export
leading_exponent <- function(C, element) {
  ev <- exponent_vector(C, element)
  d <- pc_depth(C, element)
  if (d != length(C$pcgs) + 1L) return(ev[d])
  NULL
}

.pc_sift <- function(C, z, g) {
  h <- g
  d <- pc_depth(C, h)
  while (d < length(C$pcgs) && !identical(z[[d]], 1L)) {
    k <- z[[d]]
    e <- leading_exponent(C, h) * (1 / leading_exponent(C, k))
    e <- e %% C$relative_order[d]
    h <- .fpp_mul(.fpp_pow(k, -e), h)
    d <- pc_depth(C, h)
  }
  h
}

#' @rdname PolycyclicGroup
#' @export
induced_pcgs <- function(C, gens) {
  z <- rep(list(1L), length(C$pcgs))
  G <- as.list(gens)
  while (length(G)) {
    g <- G[[1L]]
    G <- G[-1L]
    h <- .pc_sift(C, z, g)
    d <- pc_depth(C, h)
    if (d <= length(C$pcgs) && identical(z[[d]], 1L)) {
      for (gen in z) {
        if (!identical(gen, 1L)) G[[length(G) + 1L]] <- .fpp_mul(.fpp_mul(.fpp_mul(.fpp_inv(h), .fpp_inv(gen)), h), gen)
      }
      z[[d]] <- h
    }
  }
  Filter(function(gen) !identical(gen, 1L), z)
}

#' @rdname PolycyclicGroup
#' @export
constructive_membership_test <- function(C, ipcgs, g) {
  e <- integer(length(ipcgs))
  h <- g
  d <- pc_depth(C, h)
  for (i in seq_along(ipcgs)) {
    gen <- ipcgs[[i]]
    while (pc_depth(C, gen) == d) {
      f <- leading_exponent(C, h) * leading_exponent(C, gen)
      f <- f %% C$relative_order[d]
      h <- .fpp_mul(.fpp_pow(gen, -f), h)
      e[i] <- f
      d <- pc_depth(C, h)
    }
  }
  if (.is_identity(h)) return(e)
  FALSE
}
