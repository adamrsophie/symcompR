# Two dimensional pretty printing, after sympy/printing/pretty (stringPict
# boxes): powers as raised exponents, fractions over a bar, root signs and
# tall parentheses. Types without a 2D rule print on one line.

.pp_env <- new.env()
.pp_env$unicode <- NULL

.bx <- function(lines, base = 1L) {
  w <- max(0L, nchar(lines, type = "width"))
  pad <- w - nchar(lines, type = "width")
  list(lines = paste0(lines, strrep(" ", pad)), base = base, w = w, h = length(lines))
}
.bx_pad <- function(b, above, below) {
  sp <- strrep(" ", b$w)
  .bx(c(rep(sp, above), b$lines, rep(sp, below)), b$base + above)
}
.bx_hcat <- function(...) {
  bs <- Filter(Negate(is.null), list(...))
  if (length(bs) == 1L) return(bs[[1L]])
  up <- max(vapply(bs, function(b) b$base - 1L, 0L))
  dn <- max(vapply(bs, function(b) b$h - b$base, 0L))
  bs <- lapply(bs, function(b) .bx_pad(b, up - (b$base - 1L), dn - (b$h - b$base)))
  .bx(do.call(paste0, lapply(bs, `[[`, "lines")), up + 1L)
}
.bx_center <- function(b, w) {
  left <- (w - b$w) %/% 2L
  .bx(paste0(strrep(" ", left), b$lines, strrep(" ", w - b$w - left)), b$base)
}
.bx_frac <- function(num, den, u) {
  w <- max(num$w, den$w)
  bar <- strrep(if (u) "\u2500" else "-", w)
  .bx(c(.bx_center(num, w)$lines, bar, .bx_center(den, w)$lines), num$h + 1L)
}
.bx_sup <- function(b, e) {
  rows <- c(paste0(strrep(" ", b$w), e$lines), paste0(b$lines, strrep(" ", e$w)))
  .bx(rows, e$h + b$base)
}
.bx_paren <- function(b, u) {
  if (b$h == 1L) return(.bx(paste0("(", b$lines, ")"), 1L))
  h <- b$h
  if (u) {
    l <- c("\u239b", rep("\u239c", h - 2L), "\u239d"); r <- c("\u239e", rep("\u239f", h - 2L), "\u23a0")
  } else {
    l <- c("/", rep("|", h - 2L), "\\"); r <- c("\\", rep("|", h - 2L), "/")
  }
  .bx(paste0(l, b$lines, r), b$base)
}
.bx_root <- function(b, idx, u) {
  h <- b$h
  if (u && h == 1L && is.null(idx) && b$w == 1L) return(.bx(paste0("\u221a", b$lines), 1L))
  inner <- paste0(" ", b$lines, " ")
  pw <- h + 1L
  top <- paste0(strrep(" ", pw), strrep("_", b$w + 2L))
  rows <- vapply(seq_len(h), function(k) {
    if (k == h) paste0(if (u) "\u2572\u2571" else "\\/", strrep(" ", pw - 2L))
    else { col <- pw - (k - 1L); paste0(strrep(" ", col - 1L), if (u) "\u2571" else "/", strrep(" ", pw - col)) }
  }, "")
  out <- c(top, paste0(rows, inner))
  if (!is.null(idx)) {
    if (nchar(idx) <= pw) substr(out[1L], 1L, nchar(idx)) <- idx
    else {
      # a long index widens the margin before the root sign
      sh <- strrep(" ", nchar(idx) - pw)
      out <- c(paste0(idx, substring(out[1L], pw + 1L)), paste0(sh, out[-1L]))
    }
  }
  .bx(out, b$base + 1L)
}

# stringPict.stack: pictures from top to bottom, each line centered; the
# baseline is that of the second picture
.bx_stack <- function(...) {
  bs <- lapply(list(...), function(b) if (is.character(b)) .bx(b) else b)
  w <- max(vapply(bs, function(b) b$w, 0L))
  .bx(unlist(lapply(bs, function(b) .bx_center(b, w)$lines)), bs[[1L]]$h + bs[[2L]]$base)
}
# stringPict.below keeps the baseline of the top picture, above that of the
# bottom one
.bx_below <- function(b, ...) {
  s <- .bx_stack(b, ...)
  s$base <- b$base
  s
}
.bx_above <- function(b, ...) {
  s <- .bx_stack(..., b)
  s$base <- s$h - b$h + b$base
  s
}

# pretty_symbology.vobj for the bracket and integral signs (no middle piece)
.vobj <- function(symb, n, u) {
  v <- if (u) switch(symb, `[` = c("\u23a2", "\u23a1", "\u23a3", "["), `]` = c("\u23a5", "\u23a4", "\u23a6", "]"),
                     int = c("\u23ae", "\u2320", "\u2321", "\u222b"))
       else switch(symb, `[` = c("[", "[", "[", "["), `]` = c("]", "]", "]", "]"), int = c(" | ", "  /", "/  ", " | "))
  if (n == 1L) return(v[4L])
  c(v[2L], rep(v[1L], n - 2L), v[3L])
}

# PrettyPrinter._print_MatrixBase: the cells on a grid (columns two apart,
# rows one apart), centered in their columns, inside square brackets
.pp_matrix <- function(M, u) {
  D <- .bx("")
  if (M$nrow && M$ncol) {
    cells <- lapply(seq_len(M$nrow), function(i) lapply(seq_len(M$ncol), function(j) .pp(.mget(M, i, j), u)))
    maxw <- vapply(seq_len(M$ncol), function(j) max(vapply(cells, function(r) r[[j]]$w, 0L)), 0L)
    D <- NULL
    for (i in seq_len(M$nrow)) {
      row <- NULL
      for (j in seq_len(M$ncol)) {
        cl <- cells[[i]][[j]]
        d <- maxw[j] - cl$w
        cl <- .bx_hcat(.bx(strrep(" ", d %/% 2L)), cl, .bx(strrep(" ", d - d %/% 2L)))
        row <- if (is.null(row)) cl else .bx_hcat(row, .bx("  "), cl)
      }
      D <- if (is.null(D)) row else .bx_below(.bx_below(D, " "), row)
    }
  }
  D$base <- D$h %/% 2L + 1L
  .bx_hcat(.bx(.vobj("[", D$h, u), D$base), D, .bx(.vobj("]", D$h, u), D$base))
}

# PrettyPrinter._print_Integral
.pp_integral <- function(e, u) {
  f <- e$args[[1L]]
  lims <- e$args[-1L]
  arg <- .pp(f, u)
  if (f$type == "Add") arg <- .bx_paren(arg, u)
  for (l in lims) {
    pa <- .pp(l$args[[1L]], u)
    if (pa$w > 1L) pa <- .bx_paren(pa, u)
    arg <- .bx_hcat(arg, .bx(" d"), pa)
  }
  s <- NULL
  for (l in lims) {
    h <- arg$h
    H <- h + 2L + if (u) 0L else 2L
    pform <- .bx(.vobj("int", H, u), arg$base + (H - h) %/% 2L)
    la <- l$args
    if (length(la) > 1L) {
      if (length(la) == 2L) {
        pA <- .bx("")
        pB <- .pp(la[[2L]], u)
      } else {
        pA <- .pp(la[[2L]], u)
        pB <- .pp(la[[3L]], u)
      }
      if (!u) {
        pB <- .bx_hcat(.bx(strrep(" ", max(1L, 3L - pB$w))), pB)
        pA <- .bx_hcat(pA, .bx(strrep(" ", max(1L, 4L - pA$w))))
      }
      pform <- .bx_below(.bx_above(pform, pB), pA)
    }
    if (u) pform <- .bx_hcat(pform, .bx(" "))
    s <- if (is.null(s)) pform else .bx_hcat(pform, s)
  }
  .bx_hcat(s, arg)
}

# PrettyPrinter._print_Sum and its asum sign
.pp_asum <- function(hreq, ascii) {
  h <- max(hreq, 2L)
  d <- h %/% 2L
  w <- d + 1L
  more <- hreq %% 2L
  if (ascii) {
    lines <- c(paste0(strrep("_", w), " "), paste0("\\", strrep(" ", w - 1L), "`"))
    for (i in seq_len(d - 1L)) lines <- c(lines, paste0(strrep(" ", i), "\\", strrep(" ", w - i)))
    if (more) lines <- c(lines, paste0(strrep(" ", d), ")", strrep(" ", w - d)))
    for (i in rev(seq_len(d - 1L))) lines <- c(lines, paste0(strrep(" ", i), "/", strrep(" ", w - i)))
    lines <- c(lines, paste0("/", strrep("_", w - 1L), ","))
    return(list(d = d, h = h + more, lines = lines, more = more))
  }
  w <- w + more
  d <- d + more
  lines <- strrep("_", w)
  for (i in 0:(d - 1L)) lines <- c(lines, paste0(strrep(" ", i), "\u2572", strrep(" ", w - i - 1L)))
  for (i in (d - 1L):0) lines <- c(lines, paste0(strrep(" ", i), "\u2571", strrep(" ", w - i - 1L)))
  list(d = d, h = h + 2L * more, lines = c(lines, strrep("\u203e", w)), more = more)
}

.pp_sum <- function(e, u) {
  f <- e$args[[1L]]
  pf <- .pp(f, u)
  if (f$type == "Add") pf <- .bx_paren(pf, u)
  H <- pf$h + 2L
  first <- TRUE
  max_upper <- 0L
  sign_height <- 0L
  more <- 0L
  for (l in e$args[-1L]) {
    la <- l$args
    pU <- .pp(la[[3L]], u)
    pL <- .bx_hcat(.pp(la[[1L]], u), .bx(" = "), .pp(la[[2L]], u))
    max_upper <- max(max_upper, pU$h)
    a <- .pp_asum(H, !u)
    sg <- do.call(.bx_stack, as.list(a$lines))
    if (first) sign_height <- sg$h
    sg <- .bx_below(.bx_above(sg, pU), pL)
    if (first) {
      pf$base <- pf$base - (a$d - (pf$h %/% 2L - (pf$base - 1L)))
      first <- FALSE
    }
    sg <- .bx_hcat(sg, do.call(.bx_stack, as.list(rep(" ", a$h))))
    pf <- .bx_hcat(sg, pf)
    more <- a$more
  }
  pf$base <- max_upper + sign_height %/% 2L + (if (!u && !more) 1L else 0L) + 1L
  pf
}

# prettyForm.__pow__ wraps an exponent that is itself drawn as a power
.pp_is_pow <- function(z) {
  if (z$type == "Function") return(z$name == "exp" && length(z$args) == 1L)
  if (z$type != "Pow") return(FALSE)
  b <- z$args[[1L]]
  x <- z$args[[2L]]
  if (x$type == "Rational" && gmp::numerator(x$value) == 1) return(FALSE)
  if (is_number_type(x) && num_sign(x) < 0) return(FALSE)
  !(b$type == "Function" && length(b$args) && x$type == "Integer")
}
.pp_exp <- function(x, u) {
  b <- .pp(x, u)
  if (.pp_is_pow(x)) .bx_paren(b, u) else b
}

.pp <- function(e, u) {
  if (!is.null(.pp_env$hook) && !is.null(r <- .pp_env$hook(e, u))) return(r)
  txt <- function(s) .bx(s)
  par_if <- function(z, lvl) { b <- .pp(z, u); if (.prec(z) < lvl) .bx_paren(b, u) else b }
  switch(e$type,
    Integer = , Float = txt(.str(e)),
    Rational = {
      b <- .bx_frac(txt(as.character(abs(gmp::numerator(e$value)))), txt(as.character(gmp::denominator(e$value))), u)
      if (e$value < 0) .bx_hcat(txt("- "), b) else b
    },
    Symbol = , Dummy = txt(.str(e)),
    Constant = txt(switch(EXPR = e$name, pi = if (u) "\u03c0" else "pi", E = if (u) "\u212f" else "E", I = if (u) "\u2148" else "I",
                          oo = if (u) "\u221e" else "oo", `-oo` = if (u) "-\u221e" else "-oo", .str(e))),
    Add = {
      terms <- .add_print_order(e$args)
      out <- NULL
      for (i in seq_along(terms)) {
        z <- terms[[i]]
        neg <- could_extract_minus_sign(z)
        b <- if (neg) .pp(if (z$type == "Mul") .neg_mul_noflat(z) else .neg(z), u) else .pp(z, u)
        if (neg && z$type == "Add") b <- .bx_paren(b, u)
        sign <- if (i == 1L) { if (neg) (if (b$h > 1L) "- " else "-") else "" } else if (neg) " - " else " + "
        out <- if (is.null(out)) .bx_hcat(if (nzchar(sign)) txt(sign), b) else .bx_hcat(out, txt(sign), b)
      }
      out
    },
    Mul = {
      if (could_extract_minus_sign(e)) {
        b <- .pp(.neg_mul_noflat(e), u)
        return(.bx_hcat(txt(if (b$h > 1L) "- " else "-"), b))
      }
      nd <- fraction(e, exact = TRUE)
      if (!is_one_num(nd[[2L]])) return(.bx_frac(.pp(nd[[1L]], u), .pp(nd[[2L]], u), u))
      fs <- .mul_print_order(e$args)
      bs <- lapply(fs, function(z) par_if(z, PREC$Mul))
      dot <- txt(if (u) "\u22c5" else "*")
      Reduce(function(a, b) .bx_hcat(a, dot, b), bs)
    },
    Pow = {
      b <- e$args[[1L]]; x <- e$args[[2L]]
      if (identical(x$key, S_Half$key)) return(.bx_root(.pp(b, u), NULL, u))
      if (x$type == "Rational" && gmp::numerator(x$value) == 1)
        return(.bx_root(.pp(b, u), as.character(gmp::denominator(x$value)), u))
      if (is_number_type(x) && num_sign(x) < 0) return(.bx_frac(txt("1"), .pp(.pow(b, .neg(x)), u), u))
      if (is_const(b, "E")) return(.bx_sup(txt(if (u) "\u212f" else "e"), .pp_exp(x, u)))
      if (b$type == "Function" && length(b$args) && is_exact_number(x) && x$type == "Integer") {
        nm <- .bx_sup(txt(b$name), .pp(x, u))
        args <- Reduce(function(a, c) .bx_hcat(a, txt(", "), c), lapply(b$args, .pp, u = u))
        return(.bx_hcat(nm, .bx_paren(args, u)))
      }
      .bx_sup(par_if(b, PREC$Pow + 1), .pp_exp(x, u))
    },
    Function = {
      if (e$name == "exp" && length(e$args) == 1L) return(.bx_sup(txt(if (u) "\u212f" else "e"), .pp_exp(e$args[[1L]], u)))
      if (!length(e$args)) return(txt(paste0(e$name, "()")))
      args <- Reduce(function(a, c) .bx_hcat(a, txt(", "), c), lapply(e$args, .pp, u = u))
      .bx_hcat(txt(e$name), .bx_paren(args, u))
    },
    Integral = .pp_integral(e, u),
    Sum = .pp_sum(e, u),
    txt(.str(e)))
}

#' Pretty printing
#'
#' `pretty(e)` lays an expression out in two dimensions, as SymPy's
#' `pretty()`: raised exponents, fractions over a bar, root signs and tall
#' parentheses, in Unicode or ASCII. `pprint()` and `pretty_print()` print
#' it. `pprint_use_unicode(flag)` sets the default and returns the previous
#' value; `pprint_try_use_unicode()` turns Unicode on when the session
#' encoding is UTF-8. `pager_print()` shows the result in the pager.
#'
#' @param e,x An expression.
#' @param ... Unused.
#' @param use_unicode Unicode or ASCII; the session default when `NULL`.
#' @param flag `TRUE`, `FALSE`, or `NULL` to only query.
#' @return A character scalar; the printing functions return `NULL`
#'   invisibly.
#' @examples
#' x <- Symbol("x"); y <- Symbol("y")
#' cat(pretty((x + 1)/(y - 2), use_unicode = FALSE), "\n")
#' pprint(sqrt(x^2 + 1), use_unicode = FALSE)
#' @name pretty
NULL

.pp_unicode <- function(use_unicode) {
  if (!is.null(use_unicode)) return(isTRUE(use_unicode))
  if (!is.null(.pp_env$unicode)) return(.pp_env$unicode)
  isTRUE(l10n_info()[["UTF-8"]])
}

#' @rdname pretty
#' @export
pretty.sym <- function(x, use_unicode = NULL, ...) {
  b <- .pp(as_sym(x), .pp_unicode(use_unicode))
  paste(b$lines, collapse = "\n")
}

#' @rdname pretty
#' @export
pretty.symmatrix <- function(x, use_unicode = NULL, ...) {
  b <- .pp_matrix(x, .pp_unicode(use_unicode))
  paste(b$lines, collapse = "\n")
}

#' @rdname pretty
#' @export
pprint <- function(e, use_unicode = NULL) {
  e <- as_sym(e)
  cat(if (inherits(e, "symmatrix")) pretty.symmatrix(e, use_unicode) else pretty.sym(e, use_unicode), "\n", sep = "")
  invisible(NULL)
}

#' @rdname pretty
#' @export
pretty_print <- pprint

#' @rdname pretty
#' @export
pprint_use_unicode <- function(flag = NULL) {
  old <- .pp_unicode(NULL)
  if (!is.null(flag)) .pp_env$unicode <- isTRUE(flag)
  old
}

#' @rdname pretty
#' @export
pprint_try_use_unicode <- function() {
  .pp_env$unicode <- isTRUE(l10n_info()[["UTF-8"]])
  invisible(NULL)
}

#' @rdname pretty
#' @export
pager_print <- function(e, use_unicode = NULL) {
  f <- tempfile(fileext = ".txt")
  writeLines(pretty(as_sym(e), use_unicode = use_unicode), f)
  file.show(f, delete.file = TRUE)
  invisible(NULL)
}
