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
      if (is_const(b, "E")) return(.bx_sup(txt(if (u) "\u212f" else "e"), .pp(x, u)))
      if (b$type == "Function" && length(b$args) && is_exact_number(x) && x$type == "Integer") {
        nm <- .bx_sup(txt(b$name), .pp(x, u))
        args <- Reduce(function(a, c) .bx_hcat(a, txt(", "), c), lapply(b$args, .pp, u = u))
        return(.bx_hcat(nm, .bx_paren(args, u)))
      }
      .bx_sup(par_if(b, PREC$Pow + 1), .pp(x, u))
    },
    Function = {
      if (e$name == "exp" && length(e$args) == 1L) return(.bx_sup(txt(if (u) "\u212f" else "e"), .pp(e$args[[1L]], u)))
      if (!length(e$args)) return(txt(paste0(e$name, "()")))
      args <- Reduce(function(a, c) .bx_hcat(a, txt(", "), c), lapply(e$args, .pp, u = u))
      .bx_hcat(txt(e$name), .bx_paren(args, u))
    },
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
pprint <- function(e, use_unicode = NULL) {
  cat(pretty(as_sym(e), use_unicode = use_unicode), "\n", sep = "")
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
