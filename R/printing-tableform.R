# TableForm, after sympy/printing/tableform.py: a table of expressions with
# optional row and column headings, printed as text or as a LaTeX tabular.

#' Tables of expressions
#'
#' `TableForm(data, ...)` arranges expressions in a table, as SymPy's
#' `TableForm`. `data` is a list of rows, an R matrix, or an explicit
#' symbolic matrix; short rows are padded. `format()` and `print()` give the
#' text table, [latex()] a `tabular` environment, and `as.matrix()` the
#' entries as a symbolic matrix.
#'
#' @param data A list of rows, a matrix, or a symbolic matrix.
#' @param headings `NULL`, `"automatic"`, or a list of row and column
#'   headings, either of which may be `NULL` or `"automatic"`.
#' @param alignments `"l"`, `"r"`, `"c"` (or `"left"`, `"right"`,
#'   `"center"`, `"<"`, `">"`, `"^"`) for every column, or one per column,
#'   with an optional first one for the row headings.
#' @param formats `NULL`, or one entry per column for LaTeX output: `NULL`,
#'   a `sprintf()` format, or a function of the entry and its zero based row
#'   and column that returns a string (or `NULL` for the default).
#' @param pad The entry used for missing values and short rows.
#' @param wipe_zeros Print zero entries as blanks.
#' @param x A `TableForm`.
#' @param ... Unused.
#' @return A `TableForm` object; `format()` returns a string.
#' @examples
#' x <- Symbol("x")
#' TableForm(list(list(5, 7), list(4, 2), list(10, 3)), headings = "automatic")
#' latex(TableForm(list(list(x^2, 1/x)), headings = list(NULL, c("a", "b"))))
#' @export
TableForm <- function(data, headings = NULL, alignments = "l", formats = NULL, pad = NULL, wipe_zeros = TRUE) {
  if (inherits(data, "symmatrix")) data <- .mrows(data)
  else if (is.matrix(data)) data <- lapply(seq_len(nrow(data)), function(i) as.list(data[i, ]))
  data <- lapply(data, as.list)
  ok_none <- is.null(pad)
  padsym <- Symbol(if (ok_none) " " else as.character(pad))
  h <- length(data)
  w <- max(vapply(data, length, 0L))
  entry <- function(v) {
    if (is.null(v) || (length(v) == 1L && !is_sym(v) && is.na(v))) return(if (ok_none) NULL else padsym)
    if (is.character(v)) return(tryCatch(sympify(v), error = function(e) Symbol(v)))
    as_sym(v)
  }
  lines <- lapply(data, function(r) {
    out <- vector("list", w)
    for (j in seq_len(w)) {
      e <- if (j <= length(r)) entry(r[[j]]) else padsym
      if (!is.null(e)) out[[j]] <- e
    }
    out
  })
  if (identical(headings, "automatic")) headings <- list("automatic", "automatic")
  if (is.null(headings)) headings <- list(NULL, NULL)
  hd <- lapply(1:2, function(k) {
    v <- headings[[k]]
    if (identical(v, "automatic")) as.character(seq_len(c(h, w)[k])) else if (is.null(v)) NULL else vapply(as.list(v), function(z) if (is_sym(z)) .str(z) else as.character(z), "")
  })
  std <- function(a) {
    a <- tolower(trimws(a))
    if (nchar(a) > 1L) c(left = "l", right = "r", center = "c")[a] %|NA|% a
    else c(`<` = "l", `>` = "r", `^` = "c")[a] %|NA|% a
  }
  al <- vapply(alignments, std, "", USE.NAMES = FALSE)
  # a compact string such as "rcl" gives one alignment per character
  if (length(al) == 1L && !(al %in% c("l", "r", "c")) && nchar(alignments[1L]) > 1L)
    al <- vapply(strsplit(alignments[1L], "")[[1L]], std, "", USE.NAMES = FALSE)
  if (length(al) == 1L && al %in% c("l", "r", "c")) al <- rep(al, w)
  if (!all(al %in% c("l", "r", "c"))) stop("alignment \"", paste(alignments, collapse = ", "), "\" unrecognized", call. = FALSE)
  head_align <- "r"
  if (!is.null(hd[[1L]]) && length(al) == w + 1L) { head_align <- al[1L]; al <- al[-1L] }
  if (length(al) != w) stop("wrong number of alignments: expected ", w, " but got ", length(al), call. = FALSE)
  if (is.null(formats)) formats <- vector("list", w)
  structure(list(lines = lines, w = w, h = h, headings = hd, head_align = head_align, alignments = al,
                 formats = as.list(formats), wipe_zeros = wipe_zeros), class = "TableForm")
}

`%|NA|%` <- function(a, b) if (is.na(a)) b else unname(a)

# Python's str.center: the odd space goes left when the width is odd.
.py_center <- function(s, width) {
  marg <- width - nchar(s, type = "width")
  if (marg <= 0L) return(s)
  left <- marg %/% 2L + bitwAnd(bitwAnd(marg, width), 1L)
  paste0(strrep(" ", left), s, strrep(" ", marg - left))
}
.py_pad <- function(s, width, align) {
  p <- strrep(" ", max(0L, width - nchar(s, type = "width")))
  if (align == "l") paste0(s, p) else paste0(p, s)
}

#' @rdname TableForm
#' @export
format.TableForm <- function(x, ...) {
  widths <- integer(x$w)
  cells <- lapply(x$lines, function(r) vapply(seq_len(x$w), function(j) {
    s <- if (is.null(r[[j]])) "None" else .str(r[[j]])
    if (x$wipe_zeros && s == "0") s <- " "
    s
  }, ""))
  for (r in cells) widths <- pmax(widths, nchar(r, type = "width"))
  rh <- x$headings[[1L]]; ch <- x$headings[[2L]]
  if (!is.null(ch)) widths <- pmax(widths, nchar(ch, type = "width"))
  hw <- if (!is.null(rh)) max(nchar(rh, type = "width")) else 0L
  row <- function(d, rlab) {
    parts <- vapply(seq_len(x$w), function(j) .py_pad(d[j], widths[j], x$alignments[j]), "")
    if (!is.null(rh)) parts <- c(.py_pad(rlab, hw, x$head_align), "|", parts)
    paste(parts, collapse = " ")
  }
  out <- character()
  if (!is.null(ch)) {
    first <- row(ch, "")
    out <- c(out, first, strrep("-", nchar(first, type = "width")))
  }
  for (i in seq_along(cells)) {
    d <- vapply(seq_len(x$w), function(j) if (x$alignments[j] == "c") .py_center(cells[[i]][j], widths[j]) else cells[[i]][j], "")
    rl <- if (!is.null(rh)) { l <- rh[i]; if (x$head_align == "c") .py_center(l, hw) else l } else ""
    out <- c(out, row(d, rl))
  }
  paste(out, collapse = "\n")
}

#' @rdname TableForm
#' @export
print.TableForm <- function(x, ...) {
  cat(format(x), "\n", sep = "")
  invisible(x)
}

#' @rdname TableForm
#' @export
as.matrix.TableForm <- function(x, ...) Matrix(lapply(x$lines, function(r) lapply(r, function(e) if (is.null(e)) Symbol("None") else e)))

.tableform_latex <- function(x) {
  rh <- x$headings[[1L]]; ch <- x$headings[[2L]]
  al <- c(if (!is.null(rh)) x$head_align, x$alignments)
  s <- paste0("\\begin{tabular}{", paste(al, collapse = " "), "}\n")
  if (!is.null(ch)) s <- paste0(s, paste(c(if (!is.null(rh)) "", ch), collapse = " & "), " \\\\\n\\hline\n")
  for (i in seq_along(x$lines)) {
    d <- vapply(seq_len(x$w), function(j) {
      e <- x$lines[[i]][[j]]
      if (!is.null(e) && x$wipe_zeros && is_zero_num(e)) return(" ")
      f <- x$formats[[j]]
      if (is.function(f)) {
        v <- f(e, i - 1L, j - 1L)
        if (!is.null(v)) return(v)
        return(if (is.null(e)) "\\text{None}" else .latex(e))
      }
      if (is.null(e)) return("$\\text{None}$")
      if (is.character(f)) {
        # integers are formatted exactly; other numbers from their exact value
        if (grepl("%[-+ #0]*[0-9]*(\\.[0-9]*)?[di]", f) && e$type == "Integer")
          return(.py_format_int(f, gmp::as.bigz(e$value)))
        num <- grepl("%[-+ #0-9.]*[dieEfgG]", f)
        dbl <- if (is_exact_number(e)) as.numeric(gmp::as.bigq(e$value)) else num_value_dbl(N(e))
        val <- if (!num) .str(e) else if (grepl("%[-+ #0-9.]*[di]", f)) as.integer(dbl) else dbl
        return(sprintf(f, val))
      }
      paste0("$", .latex(e), "$")
    }, "")
    if (!is.null(rh)) d <- c(rh[i], d)
    s <- paste0(s, paste(d, collapse = " & "), " \\\\\n")
  }
  paste0(s, "\\end{tabular}")
}

# Python's "%d" formatting of an exact integer: flags, width and precision
# (minimum digits) applied to the decimal digits.
.py_format_int <- function(f, v) {
  m <- regmatches(f, regexec("%([-+ #0]*)([0-9]*)(\\.([0-9]*))?[di]", f))[[1L]]
  flags <- m[2L]; width <- if (nzchar(m[3L])) as.integer(m[3L]) else 0L
  prec <- if (nzchar(m[4L])) as.integer(if (nzchar(m[5L])) m[5L] else "0") else NA_integer_
  digits <- as.character(abs(v))
  if (!is.na(prec) && nchar(digits) < prec) digits <- paste0(strrep("0", prec - nchar(digits)), digits)
  sgn <- if (v < 0) "-" else if (grepl("+", flags, fixed = TRUE)) "+" else if (grepl(" ", flags, fixed = TRUE)) " " else ""
  body <- paste0(sgn, digits)
  pad <- width - nchar(body)
  if (pad > 0L) {
    body <- if (grepl("-", flags, fixed = TRUE)) paste0(body, strrep(" ", pad))
            else if (grepl("0", flags, fixed = TRUE) && is.na(prec)) paste0(sgn, strrep("0", pad), digits)
            else paste0(strrep(" ", pad), body)
  }
  pos <- regexpr(m[1L], f, fixed = TRUE)
  out <- paste0(substr(f, 1L, pos - 1L), body, substr(f, pos + nchar(m[1L]), nchar(f)))
  gsub("%%", "%", out, fixed = TRUE)
}
