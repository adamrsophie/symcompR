# Iterators of sympy/utilities/iterables.py. Python generators become R
# functions returning a list of all items in the order SymPy yields them
# (each item a copy, so items do not alias as SymPy's reused lists do).

# itertools.permutations(v, r): positions in lexicographic order
.it_perms <- function(v, r = length(v)) {
  n <- length(v)
  if (r > n) return(list())
  if (r == 0L) return(list(v[0L]))
  out <- list()
  rec <- function(pre, rest) {
    if (length(pre) == r) { out[[length(out) + 1L]] <<- v[pre]; return(invisible()) }
    for (i in rest) rec(c(pre, i), setdiff(rest, i))
  }
  rec(integer(), seq_len(n))
  out
}
# itertools.combinations(v, r)
.it_combs <- function(v, r) {
  n <- length(v)
  if (r > n) return(list())
  if (r == 0L) return(list(v[0L]))
  cm <- utils::combn(n, r, simplify = FALSE)
  lapply(cm, function(i) v[i])
}
# sympy.core.sorting.ordered for plain R values and sym expressions
.it_ordered <- function(x) {
  if (!length(x)) return(x)
  if (is.list(x) && all(vapply(x, is_sym, NA))) return(.bool_print_order(x))
  if (is.list(x)) return(x)
  sort(x, method = "radix")
}
.it_same <- function(a, b) if (is_sym(a) && is_sym(b)) identical(a$key, b$key) else identical(a, b)
.it_key <- function(a) if (is_sym(a)) a$key else if (is.character(a)) paste0("s:", a) else paste0("n:", format(a, digits = 17))

#' Iterables
#'
#' Ports of the iterators of `sympy.utilities.iterables`, each returning a
#' list of everything the SymPy generator yields, in SymPy's order.
#' Sequences are R vectors or lists.
#'
#' @param seq,x,l,s,t A sequence (vector or list); for `take()` also a generator.
#' @param g Internal: the grouped form of `m`, as in SymPy.
#' @param m A sequence, or a named vector of counts, for the multiset
#'   functions; for `partitions()` the largest number of parts.
#' @param how A nested list of counts.
#' @param n,y,size Integers.
#' @param k An integer; for `partitions()` the largest part, for
#'   `multiset_partitions()` the number of parts.
#' @param bits Number of bits, or `"all"` for all bit strings of length `n`.
#' @param str Return strings.
#' @param G A graph `list(V, E)` with `E` a list of pairs.
#' @param key A key function for comparisons.
#' @param directed Directed rotation only.
#' @param op A binary predicate for `runs()`.
#' @param ordered `NULL`, `0`, `1`, `10` or `11`, as SymPy's `kbins()`.
#' @param sort Sort the partitions.
#' @param choice A function choosing an element of a vector.
#' @param strict Fail when no derangement exists.
#' @param free Count rotations and reflections as equal.
#' @param ... Sequences.
#' @return A list (or a vector or string where SymPy returns one).
#' @examples
#' ordered_partitions(5, 3)
#' kbins(1:4, 2)
#' multiset_permutations(c(1, 1, 2))
#' binary_partitions(5)
#' @name iterables
NULL

#' @rdname iterables
#' @export
reshape_seq <- function(seq, how) {
  sizes <- unlist(how)
  m <- sum(sizes)
  if (any(sizes < 0) || m <= 0 || length(seq) %% m)
    stop("template must contain nonnegative counts with a positive sum that divides the sequence length", call. = FALSE)
  i <- 0L
  make_group <- function(template) {
    res <- list()
    for (item in template) {
      if (is.numeric(item) && length(item) == 1L && !is.list(item)) {
        if (item > 0) res <- c(res, as.list(seq[(i + 1L):(i + item)]))
        i <<- i + item
      } else res[[length(res) + 1L]] <- make_group(item)
    }
    # a group of plain elements is a vector, as type(template)(result)
    if (is.atomic(seq) && all(vapply(res, function(r) is.atomic(r) && length(r) == 1L, NA))) unlist(res) else res
  }
  lapply(seq_len(length(seq) %/% m), function(r) make_group(how))
}

#' @rdname iterables
#' @export
multiset <- function(seq) {
  ks <- vapply(seq, .it_key, "")
  u <- !duplicated(ks)
  stats::setNames(as.integer(table(factor(ks, levels = ks[u]))), if (is.atomic(seq)) as.character(seq[u]) else ks[u])
}

#' @rdname iterables
#' @export
ibin <- function(n, bits = NULL, str = FALSE) {
  n <- as.numeric(n)
  if (n < 0) stop("negative numbers are not allowed", call. = FALSE)
  allbits <- is.character(bits)
  b <- if (is.null(bits)) 0 else if (allbits) -1 else as.numeric(bits)
  bin <- function(v, w) {
    s <- if (v == 0) "0" else paste(rev(as.integer(intToBits(v))[seq_len(floor(log2(v)) + 1)]), collapse = "")
    if (nchar(s) < w) s <- paste0(strrep("0", w - nchar(s)), s)
    s
  }
  if (b >= 0 && n > 0 && floor(log2(n)) + 1 > b && !is.null(bits))
    stop("`bits` must be >= ", floor(log2(n)) + 1, call. = FALSE)
  if (!str) {
    if (b >= 0) return(as.integer(strsplit(bin(n, b), "")[[1L]]))
    return(variations(0:1, n, repetition = TRUE))
  }
  if (b >= 0) return(bin(n, b))
  vapply(seq_len(2^n) - 1, function(i) bin(i, n), "")
}

#' @rdname iterables
#' @export
common_prefix <- function(...) {
  seqs <- list(...)
  if (!all(lengths(seqs) > 0)) return(seqs[[1L]][0L])
  if (length(seqs) == 1L) return(seqs[[1L]])
  m <- min(lengths(seqs)); i <- 0L
  for (j in seq_len(m)) {
    if (!all(vapply(seqs, function(s) .it_same(s[[j]], seqs[[1L]][[j]]), NA))) break
    i <- j
  }
  seqs[[1L]][seq_len(i)]
}

#' @rdname iterables
#' @export
common_suffix <- function(...) {
  seqs <- list(...)
  if (!all(lengths(seqs) > 0)) return(seqs[[1L]][0L])
  if (length(seqs) == 1L) return(seqs[[1L]])
  m <- min(lengths(seqs)); i <- 0L
  for (j in seq_len(m)) {
    if (!all(vapply(seqs, function(s) .it_same(s[[length(s) - j + 1L]], seqs[[1L]][[length(seqs[[1L]]) - j + 1L]]), NA))) break
    i <- j
  }
  s1 <- seqs[[1L]]
  if (!i) s1[0L] else s1[(length(s1) - i + 1L):length(s1)]
}

#' @rdname iterables
#' @export
strongly_connected_components <- function(G) {
  V <- G[[1L]]; E <- G[[2L]]
  vk <- vapply(V, .it_key, "")
  Gmap <- stats::setNames(rep(list(integer()), length(V)), vk)
  for (e in E) { a <- .it_key(e[[1L]]); Gmap[[a]] <- c(Gmap[[a]], match(.it_key(e[[2L]]), vk)) }
  .scc(V, vk, Gmap)
}

# Tarjan's algorithm, iterative, as _strongly_connected_components
.scc <- function(V, vk, Gmap) {
  lowlink <- indices <- rep(NA_integer_, length(V))
  stack <- integer(); onst <- logical(length(V))
  callstack <- list(); comps <- list()
  start <- function(v) {
    idx <- length(stack)
    indices[v] <<- lowlink[v] <<- idx
    stack <<- c(stack, v); onst[v] <<- TRUE
    callstack[[length(callstack) + 1L]] <<- list(v = v, it = Gmap[[vk[v]]], pos = 0L)
  }
  finish <- function(v1) {
    if (lowlink[v1] == indices[v1]) {
      comp <- integer()
      repeat {
        top <- stack[length(stack)]; stack <<- stack[-length(stack)]; onst[top] <<- FALSE
        comp <- c(comp, top)
        if (top == v1) break
      }
      comps[[length(comps) + 1L]] <<- V[rev(comp)]
    }
    v2 <- callstack[[length(callstack)]]$v
    callstack[[length(callstack)]] <<- NULL
    if (length(callstack)) {
      w <- callstack[[length(callstack)]]$v
      lowlink[w] <<- min(lowlink[w], lowlink[v2])
    }
  }
  for (v in seq_along(V)) {
    if (!is.na(indices[v])) next
    start(v)
    while (length(callstack)) {
      top <- length(callstack)
      cs <- callstack[[top]]
      if (cs$pos >= length(cs$it)) { finish(cs$v); next }
      callstack[[top]]$pos <- cs$pos + 1L
      v2 <- cs$it[cs$pos + 1L]
      if (is.na(indices[v2])) start(v2)
      else if (onst[v2]) lowlink[cs$v] <- min(lowlink[cs$v], indices[v2])
    }
  }
  comps
}

#' @rdname iterables
#' @export
connected_components <- function(G) {
  E <- unlist(lapply(G[[2L]], function(e) list(list(e[[1L]], e[[2L]]), list(e[[2L]], e[[1L]]))), recursive = FALSE)
  strongly_connected_components(list(G[[1L]], E))
}

#' @rdname iterables
#' @export
rotate_left <- function(x, y) {
  if (!length(x)) return(x)
  y <- y %% length(x)
  c(x[seq_along(x) > y], x[seq_len(y)])
}
#' @rdname iterables
#' @export
rotate_right <- function(x, y) {
  if (!length(x)) return(x)
  rotate_left(x, length(x) - y %% length(x))
}

#' @rdname iterables
#' @export
least_rotation <- function(x, key = identity) {
  S <- c(x, x)
  f <- rep(-1L, length(S))
  k <- 0L
  lt <- function(a, b) { ka <- key(a); kb <- key(b); if (is_sym(ka)) .is_true(Lt(ka, kb)) else ka < kb }
  ne <- function(a, b) !.it_same(a, b)
  for (j in seq_len(length(S) - 1L)) {
    sj <- S[[j + 1L]]
    i <- f[j - k]
    while (i != -1L && ne(sj, S[[k + i + 2L]])) {
      if (lt(sj, S[[k + i + 2L]])) k <- j - i - 1L
      i <- f[i + 1L]
    }
    if (ne(sj, S[[k + i + 2L]])) {
      if (lt(sj, S[[k + 1L]])) k <- j
      f[j - k + 1L] <- -1L
    } else f[j - k + 1L] <- i + 1L
  }
  k
}

#' @rdname iterables
#' @export
minlex <- function(seq, directed = TRUE, key = identity) {
  isstr <- is.character(seq) && length(seq) == 1L
  s <- if (isstr) strsplit(seq, "")[[1L]] else seq
  best <- rotate_left(s, least_rotation(s, key))
  if (!directed) {
    rs <- rev(s)
    rb <- rotate_left(rs, least_rotation(rs, key))
    # min(best, rbest, key=key): tuple comparison of the keys
    cmp <- 0L
    for (i in seq_along(best)) {
      a <- key(best[[i]]); b <- key(rb[[i]])
      if (!identical(a, b)) { cmp <- if (a < b) -1L else 1L; break }
    }
    if (cmp > 0L) best <- rb
  }
  if (isstr) paste(best, collapse = "") else best
}

.ms_groups <- function(m) {
  if (!is.null(names(m)) && is.numeric(m) && !is.list(m) && all(nzchar(names(m)))) {
    if (any(m < 0)) stop("counts cannot be negative", call. = FALSE)
    ks <- .it_ordered(names(m))
    return(lapply(ks, function(k) list(k, as.integer(m[[k]]))))
  }
  v <- .it_ordered(m)
  out <- list()
  for (i in seq_along(v)) {
    if (length(out) && .it_same(out[[length(out)]][[1L]], v[[i]])) out[[length(out)]][[2L]] <- out[[length(out)]][[2L]] + 1L
    else out[[length(out) + 1L]] <- list(v[[i]], 1L)
  }
  out
}
.rep_el <- function(k, v, proto) if (is.list(proto)) rep(list(k), v) else rep(k, v)

#' @rdname iterables
#' @export
multiset_combinations <- function(m, n, g = NULL) {
  proto <- m
  if (is.null(g)) {
    g <- .ms_groups(m)
    N <- sum(vapply(g, `[[`, 0L, 2L))
    if (n > N) return(list())
  } else N <- sum(vapply(g, `[[`, 0L, 2L))
  empty <- if (is.list(proto)) list() else proto[0L]
  if (n > N || !n) return(list(empty))
  out <- list()
  for (i in seq_along(g)) {
    k <- g[[i]][[1L]]; v <- g[[i]][[2L]]
    if (v >= n) { out[[length(out) + 1L]] <- .rep_el(k, n, proto); v <- n - 1L }
    for (vv in rev(seq_len(min(n, v)))) {
      for (j in multiset_combinations(proto, n - vv, if (i < length(g)) g[(i + 1L):length(g)] else list())) {
        rv <- c(.rep_el(k, vv, proto), j)
        if (length(rv) == n) out[[length(out) + 1L]] <- rv
      }
    }
  }
  out
}

#' @rdname iterables
#' @export
multiset_permutations <- function(m, size = NULL, g = NULL) {
  proto <- if (is.null(g)) m else g$proto
  do <- if (is.null(g)) .ms_groups(m) else g$do
  do <- Filter(function(d) d[[2L]] > 0, do)
  SUM <- sum(vapply(do, `[[`, 0L, 2L))
  empty <- if (is.list(proto)) list() else proto[0L]
  if (!length(do) || (!is.null(size) && (size > SUM || size < 1))) {
    if ((!length(do) && is.null(size)) || identical(size, 0) || identical(size, 0L)) return(list(empty))
    return(list())
  }
  if (identical(as.numeric(size), 1)) return(lapply(do, function(d) .rep_el(d[[1L]], 1L, proto)))
  if (length(do) == 1L) {
    k <- do[[1L]][[1L]]; v <- do[[1L]][[2L]]
    v <- if (is.null(size)) v else if (size <= v) size else 0L
    return(list(.rep_el(k, v, proto)))
  }
  if (all(vapply(do, `[[`, 0L, 2L) == 1L)) {
    ks <- if (is.list(proto)) lapply(do, `[[`, 1L) else unlist(lapply(do, `[[`, 1L))
    return(.it_perms(ks, if (is.null(size)) length(ks) else size))
  }
  size <- if (is.null(size)) SUM else size
  out <- list()
  for (i in seq_along(do)) {
    k <- do[[i]][[1L]]
    do[[i]][[2L]] <- do[[i]][[2L]] - 1L
    for (j in multiset_permutations(NULL, size - 1L, list(proto = proto, do = do)))
      if (length(j)) out[[length(out) + 1L]] <- c(.rep_el(k, 1L, proto), j)
    do[[i]][[2L]] <- do[[i]][[2L]] + 1L
  }
  out
}

#' @rdname iterables
#' @export
ordered_partitions <- function(n, m = NULL, sort = TRUE) {
  if (n < 1 || (!is.null(m) && m < 1)) return(list(integer()))
  out <- list()
  if (is.null(m)) {
    a <- rep(1L, n); y <- -1L; v <- n
    while (v > 0) {
      v <- v - 1L
      x <- a[v + 1L] + 1L
      while (y >= 2L * x) { a[v + 1L] <- x; y <- y - x; v <- v + 1L }
      w <- v + 1L
      while (x <= y) {
        a[v + 1L] <- x; a[w + 1L] <- y
        out[[length(out) + 1L]] <- a[seq_len(w + 1L)]
        x <- x + 1L; y <- y - 1L
      }
      a[v + 1L] <- x + y; y <- a[v + 1L] - 1L
      out[[length(out) + 1L]] <- a[seq_len(w)]
    }
    return(out)
  }
  if (m == 1) return(list(as.integer(n)))
  if (n == m) return(list(rep(1L, n)))
  for (b in seq_len(n %/% m)) {
    a <- rep(as.integer(b), m); x <- n - b * m
    if (!x) {
      if (sort) out[[length(out) + 1L]] <- a
    } else if (!sort && x <= m) {
      for (ax in ordered_partitions(x, sort = FALSE)) {
        mi <- length(ax); a2 <- a
        a2[(m - mi + 1L):m] <- ax + b
        out[[length(out) + 1L]] <- a2
      }
    } else {
      for (mi in seq_len(m - 1L)) for (ax in ordered_partitions(x, mi, sort = TRUE)) {
        a2 <- a; a2[(m - mi + 1L):m] <- ax + b
        out[[length(out) + 1L]] <- a2
      }
    }
  }
  out
}

#' @rdname iterables
#' @export
binary_partitions <- function(n) {
  power <- 2^ceiling(log2(n)); acc <- 0; p <- numeric()
  while (power >= 1) {
    if (acc + power <= n) { p <- c(p, power); acc <- acc + power }
    power <- floor(power / 2)
  }
  out <- list()
  last <- length(p) - 1L - (n %% 2)
  while (last >= 0) {
    out[[length(out) + 1L]] <- as.integer(p)
    if (p[last + 1L] == 2) { p[last + 1L] <- 1; p <- c(p, 1); last <- last - 1L; next }
    p <- c(p, 1)
    p[last + 1L] <- p[last + 1L] %/% 2
    x <- p[last + 1L]; p[last + 2L] <- x
    last <- last + 1L
    while (x > 1) {
      if (x <= length(p) - last - 1L) {
        p <- p[seq_len(length(p) - x + 1L)]
        last <- last + 1L
        p[last + 1L] <- x
      } else x <- x %/% 2
    }
  }
  c(out, list(rep(1L, n)))
}

#' @rdname iterables
#' @export
uniq <- function(seq) seq[!duplicated(vapply(seq, .it_key, ""))]

#' @rdname iterables
#' @export
generate_bell <- function(n) {
  n <- as.integer(n)
  if (n < 1L) stop("n must be a positive integer", call. = FALSE)
  if (n == 1L) return(list(0L))
  if (n == 2L) return(list(c(0L, 1L), c(1L, 0L)))
  if (n == 3L) return(list(c(0L, 1L, 2L), c(0L, 2L, 1L), c(2L, 0L, 1L), c(2L, 1L, 0L), c(1L, 2L, 0L), c(1L, 0L, 2L)))
  m <- n - 1L
  op <- c(0L, rep(-1L, m)); l <- 0:(n - 1L)
  out <- list()
  repeat {
    out[[length(out) + 1L]] <- l
    bi <- NA_integer_; bv <- -1L
    for (i in seq_len(n)) if (op[i] != 0L && l[i] > bv) { bi <- i; bv <- l[i] }
    if (is.na(bi)) break
    j <- bi + op[bi]
    tmp <- l[bi]; l[bi] <- l[j]; l[j] <- tmp
    tmp <- op[bi]; op[bi] <- op[j]; op[j] <- tmp
    if (j == 1L || j == m + 1L || l[j + op[j]] > l[j]) op[j] <- 0L
    for (i in seq_len(j - 1L)) if (l[i] > l[j]) op[i] <- 1L
    if (j < n) for (i in (j + 1L):n) if (l[i] > l[j]) op[i] <- -1L
  }
  out
}

#' @rdname iterables
#' @export
generate_involutions <- function(n) {
  Filter(function(p) all(p[p + 1L] == 0:(n - 1L)), .it_perms(0:(n - 1L)))
}

.set_derangements <- function(s) {
  if (length(s) < 2L) return(list())
  if (length(s) == 2L) return(list(s[c(2L, 1L)]))
  if (length(s) == 3L) return(list(s[c(2L, 3L, 1L)], s[c(3L, 1L, 2L)]))
  Filter(function(p) !any(vapply(seq_along(p), function(i) .it_same(p[[i]], s[[i]]), NA)), .it_perms(s))
}

#' @rdname iterables
#' @export
multiset_derangements <- function(s) {
  keys <- vapply(s, .it_key, "")
  ms <- table(factor(keys, levels = unique(keys)))
  mx <- max(ms); n <- length(s)
  if (mx * 2 > n) return(list())
  if (length(ms) == n) return(.set_derangements(s))
  Mk <- names(ms)[which(ms == mx)[1L]]
  M <- s[[match(Mk, keys)]]
  inonM <- which(keys != Mk); iM <- which(keys == Mk)
  out <- list()
  na <- if (is.list(s)) list(NULL) else NA
  blank <- function() if (is.list(s)) vector("list", n) else rep(s[1L], n)
  if (2 * mx == n) {
    rv <- blank()
    for (i in inonM) rv[[i]] <- M
    for (p in multiset_permutations(s[inonM])) {
      for (t in seq_along(iM)) rv[[iM[t]]] <- p[[t]]
      out[[length(out) + 1L]] <- rv
    }
    return(out)
  }
  if (n - 2 * mx == 1 && length(ms) == n - mx + 1) {
    rv <- blank()
    for (ii in seq_along(inonM)) {
      i1 <- inonM[ii]; ifill <- inonM[-ii]
      for (j in ifill) rv[[j]] <- M
      for (p in .it_perms(s[ifill])) {
        rv[[i1]] <- s[[i1]]
        for (t in seq_along(iM)) rv[[iM[t]]] <- p[[t]]
        k <- i1
        for (j in iM) {
          tmp <- rv[[j]]; rv[[j]] <- rv[[k]]; rv[[k]] <- tmp
          out[[length(out) + 1L]] <- rv
          k <- j
        }
      }
    }
    return(out)
  }
  # take = sorted(ms.items(), key=(count, element))
  ukeys <- names(ms)
  vals <- lapply(ukeys, function(k) s[[match(k, keys)]])
  ord <- order(as.integer(ms), .it_rank(vals))
  take <- lapply(ord, function(i) list(vals[[i]], ukeys[i], as.integer(ms[i])))
  filled <- rep(FALSE, n); rv <- blank()
  finish <- function() {
    a <- take[[2L]]; b <- take[[1L]]
    fa <- integer(); fb <- integer(); of <- integer()
    for (i in seq_len(n)) if (!filled[i]) {
      if (keys[i] == a[[2L]]) fb <- c(fb, i) else if (keys[i] == b[[2L]]) fa <- c(fa, i) else of <- c(of, i)
    }
    if (length(fa) > a[[3L]] || length(fb) > b[[3L]]) return(invisible())
    r2 <- rv
    for (i in fa) r2[[i]] <- a[[1L]]
    for (i in fb) r2[[i]] <- b[[1L]]
    for (pl in .it_combs(of, a[[3L]] - length(fa))) {
      r3 <- r2
      for (pp in pl) r3[[pp]] <- a[[1L]]
      for (i in setdiff(of, pl)) r3[[i]] <- b[[1L]]
      out[[length(out) + 1L]] <<- r3
    }
  }
  do <- function(j) {
    if (j == 1L) return(finish())
    Mj <- take[[j + 1L]]
    iopen <- which(!filled & keys != Mj[[2L]])
    for (ic in .it_combs(iopen, Mj[[3L]])) {
      for (ii in ic) { rv[[ii]] <<- Mj[[1L]]; filled[ii] <<- TRUE }
      do(j - 1L)
      for (ii in ic) filled[ii] <<- FALSE
    }
  }
  do(length(take) - 1L)
  out
}
.it_rank <- function(vals) {
  o <- .it_ordered(vals)
  ko <- vapply(o, .it_key, "")
  match(vapply(vals, .it_key, ""), ko)
}

#' @rdname iterables
#' @export
generate_derangements <- function(s) {
  if (!has_dups(s)) .set_derangements(s) else multiset_derangements(s)
}

#' @rdname iterables
#' @export
random_derangement <- function(t, choice = NULL, strict = TRUE) {
  isstr <- is.character(t) && length(t) == 1L
  tt <- if (isstr) strsplit(t, "")[[1L]] else t
  if (is.null(choice)) choice <- function(v) v[[sample.int(length(v), 1L)]]
  shuffle <- function(rv) {
    if (length(rv) > 1L) for (i in rev(seq_len(length(rv) - 1L))) {
      x <- choice(rv[seq_len(i + 1L)])
      j <- which(vapply(rv, function(r) .it_same(r, x), NA))[1L]
      tmp <- rv[[i + 1L]]; rv[[i + 1L]] <- rv[[j]]; rv[[j]] <- tmp
    }
    rv
  }
  keys <- vapply(tt, .it_key, "")
  ms <- table(factor(keys, levels = unique(keys)))
  ord <- order(as.integer(ms))
  ms <- ms[ord]
  Mk <- names(ms)[length(ms)]; mx <- as.integer(ms[length(ms)])
  n <- length(tt); xs <- 2 * mx - n
  if (xs > 0) {
    if (strict) stop("no derangement possible", call. = FALSE)
    opts <- shuffle(which(keys == Mk))
    stay <- sort(opts[seq_len(xs)])
    Mv <- tt[[match(Mk, keys)]]
    rv <- random_derangement(tt[-stay], choice)
    for (i in stay) rv <- append(rv, list(Mv), after = i - 1L)
    if (!is.list(tt)) rv <- unlist(rv)
    return(if (isstr) paste(rv, collapse = "") else rv)
  }
  if (n == length(ms)) {
    rv <- tt
    repeat {
      rv <- shuffle(rv)
      if (all(!vapply(seq_len(n), function(i) .it_same(rv[[i]], tt[[i]]), NA))) break
    }
    return(if (isstr) paste(rv, collapse = "") else rv)
  }
  repeat {
    rv <- vector("list", n); ok <- TRUE
    for (jj in rev(seq_along(ms))) {
      e <- names(ms)[jj]; c <- as.integer(ms[jj])
      opts <- which(vapply(rv, is.null, NA) & keys != e)
      if (length(opts) < c) { ok <- FALSE; break }
      opts <- shuffle(opts)
      for (i in seq_len(c)) rv[[opts[i]]] <- tt[[match(e, keys)]]
    }
    if (ok) break
  }
  if (!is.list(tt)) rv <- unlist(rv)
  if (isstr) paste(rv, collapse = "") else rv
}

#' @rdname iterables
#' @export
necklaces <- function(n, k, free = FALSE) {
  if (k == 0 && n > 0) return(list())
  a <- rep(0L, n)
  out <- list(a)
  if (n == 0) return(out)
  repeat {
    i <- n - 1L
    while (a[i + 1L] == k - 1L) { i <- i - 1L; if (i == -1L) return(out) }
    a[i + 1L] <- a[i + 1L] + 1L
    for (j in seq_len(n - i - 1L) - 1L) a[j + i + 2L] <- a[j + 1L]
    if (n %% (i + 1L) == 0L) {
      ok <- TRUE
      if (free) for (j in seq_len(n - 1L) - 1L) {
        # a <= a[j::-1] + a[-1:j:-1] as list comparison
        r <- c(a[(j + 1L):1L], if (j + 2L <= n) a[n:(j + 2L)])
        d <- which(a != r)
        if (length(d) && a[d[1L]] > r[d[1L]]) { ok <- FALSE; break }
      }
      if (ok) out[[length(out) + 1L]] <- a
    }
  }
}
#' @rdname iterables
#' @export
bracelets <- function(n, k) necklaces(n, k, free = TRUE)

#' @rdname iterables
#' @export
generate_oriented_forest <- function(n) {
  P <- (-1L):(n - 1L)
  out <- list()
  repeat {
    out[[length(out) + 1L]] <- P[-1L]
    if (P[n + 1L] > 0L) P[n + 1L] <- P[P[n + 1L] + 1L]
    else {
      found <- FALSE
      for (p in rev(seq_len(n - 1L))) {
        if (P[p + 1L] != 0L) {
          target <- P[p + 1L] - 1L
          q <- 0L
          for (qq in rev(seq_len(p - 1L))) if (P[qq + 1L] == target) { q <- qq; break }
          if (q == 0L) q <- 1L
          offset <- p - q
          for (i in p:n) P[i + 1L] <- P[i - offset + 1L]
          found <- TRUE
          break
        }
      }
      if (!found) break
    }
  }
  out
}

#' @rdname iterables
#' @export
runs <- function(seq, op = function(a, b) a > b) {
  if (!length(seq)) return(list())
  cycles <- list(); run <- seq[1L]
  for (ei in seq[-1L]) {
    if (isTRUE(op(ei, run[[length(run)]]))) run <- c(run, ei)
    else { cycles[[length(cycles) + 1L]] <- run; run <- ei }
  }
  c(cycles, list(run))
}

#' @rdname iterables
#' @export
sequence_partitions <- function(l, n) {
  if (n == 1 && length(l)) return(list(list(l)))
  out <- list()
  for (i in seq_len(length(l) - 1L))
    for (part in sequence_partitions(l[-seq_len(i)], n - 1)) out[[length(out) + 1L]] <- c(list(l[seq_len(i)]), part)
  out
}

#' @rdname iterables
#' @export
sequence_partitions_empty <- function(l, n) {
  if (n < 1) return(list())
  if (n == 1) return(list(list(l)))
  out <- list()
  for (i in 0:length(l))
    for (part in sequence_partitions_empty(l[seq_along(l) > i], n - 1))
      out[[length(out) + 1L]] <- c(list(l[seq_len(i)]), part)
  out
}

#' @rdname iterables
#' @export
kbins <- function(l, k, ordered = NULL) {
  if (is.null(ordered)) return(sequence_partitions(l, k))
  if (ordered == 11) return(unlist(lapply(multiset_permutations(l), sequence_partitions, n = k), recursive = FALSE))
  if (ordered == 0) return(multiset_partitions(l, k))
  if (ordered == 10) return(unlist(lapply(multiset_partitions(l, k), .it_perms), recursive = FALSE))
  if (ordered == 1) {
    out <- list()
    for (pp in partitions(length(l), k, size = TRUE)) {
      kgot <- pp[[1L]]; p <- pp[[2L]]
      if (kgot != k) next
      sizes <- sort(as.integer(names(p)))
      for (li in multiset_permutations(l)) {
        rv <- list(); i <- 0L
        for (sz in sizes) for (mm in seq_len(p[[as.character(sz)]])) {
          rv[[length(rv) + 1L]] <- li[(i + 1L):(i + sz)]; i <- i + sz
        }
        out[[length(out) + 1L]] <- rv
      }
    }
    return(out)
  }
  stop("ordered must be one of 00, 01, 10 or 11, not ", ordered, call. = FALSE)
}

#' @rdname iterables
#' @export
permute_signs <- function(t) {
  nz <- sum(t != 0)
  if (!nz) return(list(t))
  grid <- as.matrix(expand.grid(rep(list(c(1, -1)), nz)))[, nz:1, drop = FALSE]
  lapply(seq_len(nrow(grid)), function(r) {
    signs <- grid[r, ]; out <- t; k <- nz
    # signs.pop() takes from the end
    for (i in seq_along(t)) if (t[i] != 0) { out[i] <- t[i] * signs[k]; k <- k - 1L }
    out
  })
}

#' @rdname iterables
#' @export
signed_permutations <- function(t) unlist(lapply(multiset_permutations(t), permute_signs), recursive = FALSE)

#' @rdname iterables
#' @export
roundrobin <- function(...) {
  args <- list(...)
  its <- lapply(args, as.list)
  out <- list(); i <- 1L
  while (any(lengths(its) >= i)) {
    for (it in its) if (length(it) >= i) out[[length(out) + 1L]] <- it[[i]]
    i <- i + 1L
  }
  if (all(vapply(args, is.atomic, NA))) unlist(out) else out
}

#' @rdname iterables
#' @export
partitions <- function(n, m = NULL, k = NULL, size = FALSE) {
  mm <- m; kk <- k
  emit <- function(ms) if (size) list(sum(ms), ms) else ms
  if (n <= 0 || (!is.null(mm) && mm < 1) || (!is.null(kk) && kk < 1) ||
      (!is.null(mm) && !is.null(kk) && mm * kk < n))
    return(list(emit(stats::setNames(integer(), character()))))
  m <- if (is.null(mm)) n else min(mm, n)
  k <- min(if (is.null(kk)) n else kk, n)
  n <- as.integer(n); m <- as.integer(m); k <- as.integer(k)
  # ms as a named vector in dict insertion order; keys largest first
  ms <- stats::setNames(n %/% k, as.character(k)); r <- n %% k
  keys <- k
  if (r) { ms[as.character(r)] <- 1L; keys <- c(keys, r) }
  room <- m - n %/% k - as.integer(r > 0)
  out <- list(emit(ms))
  while (!identical(keys, 1L)) {
    if (keys[length(keys)] == 1L) {
      keys <- keys[-length(keys)]
      reuse <- ms[["1"]]; ms <- ms[names(ms) != "1"]
      room <- room + reuse
    } else reuse <- 0L
    repeat {
      i <- keys[length(keys)]
      ms[[as.character(i)]] <- ms[[as.character(i)]] - 1L
      newcount <- ms[[as.character(i)]]
      reuse <- reuse + i
      if (newcount == 0L) { keys <- keys[-length(keys)]; ms <- ms[names(ms) != as.character(i)] }
      room <- room + 1L
      i <- i - 1L
      q <- reuse %/% i; r <- reuse %% i
      need <- q + as.integer(r > 0)
      if (need > room) {
        if (!length(keys)) return(out)
        next
      }
      ms[as.character(i)] <- q
      keys <- c(keys, i)
      if (r) { ms[as.character(r)] <- 1L; keys <- c(keys, r) }
      break
    }
    room <- room - need
    out[[length(out) + 1L]] <- emit(ms)
  }
  out
}

# _set_partitions(n): (number of blocks, block of each element)
.set_partitions_q <- function(n) {
  p <- integer(n); q <- integer(n); nc <- 1L
  out <- list(list(nc, q))
  while (nc != n) {
    m <- n
    repeat {
      m <- m - 1L
      i <- q[m + 1L]
      if (p[i + 1L] != 1L) break
      q[m + 1L] <- 0L
    }
    i <- i + 1L
    q[m + 1L] <- i
    m <- m + 1L
    nc <- nc + m - n
    p[1L] <- p[1L] + n - m
    if (i == nc) { p[nc + 1L] <- 0L; nc <- nc + 1L }
    p[i] <- p[i] - 1L
    p[i + 1L] <- p[i + 1L] + 1L
    out[[length(out) + 1L]] <- list(nc, q)
  }
  out
}

#' @rdname iterables
#' @export
multiset_partitions <- function(m, k = NULL) {
  if (is.numeric(m) && length(m) == 1L && !is.list(m) && m == floor(m) && is.null(names(m))) {
    n <- as.integer(m)
    if (!is.null(k) && k > n) return(list())
    ms <- 0:(n - 1L)
    if (identical(as.numeric(k), 1)) return(list(list(ms)))
    out <- list()
    for (sq in .set_partitions_q(n)) {
      nc <- sq[[1L]]; q <- sq[[2L]]
      if (is.null(k) || nc == k) out[[length(out) + 1L]] <- lapply(seq_len(nc) - 1L, function(b) ms[q == b])
    }
    return(out)
  }
  if (is.character(m) && length(m) == 1L) m <- m
  if (!has_variety(m)) {
    n <- length(m)
    if (!is.null(k) && k > n) return(list())
    if (identical(as.numeric(k), 1)) return(list(list(m)))
    x <- m[1L]
    out <- list()
    for (sp in partitions(n, k, size = TRUE)) {
      if (is.null(k) || sp[[1L]] == k) {
        p <- sp[[2L]]; rv <- list()
        for (kk in sort(as.integer(names(p)))) rv <- c(rv, rep(list(rep(x, kk)), p[[as.character(kk)]]))
        out[[length(out) + 1L]] <- rv
      }
    }
    return(out)
  }
  m <- .it_ordered(m)
  n <- length(m)
  if (!is.null(k) && k > n) return(list())
  if (identical(as.numeric(k), 1)) return(list(list(m)))
  keys <- vapply(m, .it_key, "")
  u <- !duplicated(keys)
  elements <- m[u]
  mult <- as.integer(table(factor(keys, levels = keys[u])))
  if (length(elements) < n) {
    if (!is.null(k)) return(enum_range(mult, k - 1L, k, components = elements))
    return(multiset_partitions_taocp(mult, components = elements))
  }
  out <- list()
  for (sq in .set_partitions_q(n)) {
    nc <- sq[[1L]]; q <- sq[[2L]]
    if (is.null(k) || nc == k) out[[length(out) + 1L]] <- lapply(seq_len(nc) - 1L, function(b) m[q == b])
  }
  out
}
