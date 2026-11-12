# Vector printing, after sympy/physics/vector/printing.py: functions of
# time print without their argument and time derivatives with primes (str),
# dots (LaTeX) or dots above (pretty).

.latex_mode <- new.env()

.pv_undef_t <- function(e) .is_undefined_application(e) && length(e$args) == 1L && identical(e$args[[1L]]$key, .pv_t()$key)

.pv_deriv_vars <- function(e) {
  out <- list()
  for (vc in e$args[-1L]) {
    n <- as.integer(as.numeric(vc$args[[2L]]$value))
    out <- c(out, rep(list(vc$args[[1L]]), n))
  }
  out
}

.pv_plain_str <- function(e) {
  old <- .str_mode$hook
  .str_mode$hook <- NULL
  on.exit(.str_mode$hook <- old)
  .str(e)
}

# VectorStrPrinter._print_Function and _print_Derivative
.pv_vstr_hook <- function(e) {
  if (e$type == "Function" && .is_undefined_application(e)) {
    return(gsub(paste0("(", .str(.pv_t()), ")"), "", .pv_plain_str(e), fixed = TRUE))
  }
  if (e$type == "Derivative") {
    vars <- .pv_deriv_vars(e)
    tk <- .pv_t()$key
    if (any(vapply(vars, function(v) identical(v$key, tk), TRUE)) && .is_undefined_application(e$args[[1L]]))
      return(paste0(e$args[[1L]]$name, strrep(.pv_state$dstr, length(vars))))
    return(.pv_plain_str(e))
  }
  NULL
}

.pv_with_hook <- function(hook, expr) {
  old <- .str_mode$hook
  .str_mode$hook <- hook
  on.exit(.str_mode$hook <- old)
  expr
}

.pv_any_str <- function(expr, repr = FALSE) {
  if (is.null(expr)) return("None")
  if (inherits(expr, "physics_obj")) return(format(expr))
  if (is_sym(expr)) return(.str_top(expr))
  if (inherits(expr, "symmatrix")) return(sstr(expr))
  if (inherits(expr, "pytuple")) {
    s <- vapply(expr, .pv_any_str, "", repr = TRUE)
    return(paste0("(", paste(s, collapse = ", "), if (length(s) == 1L) ",", ")"))
  }
  if (inherits(expr, "pydict")) {
    if (!length(expr$keys)) return("{}")
    return(paste0("{", paste(vapply(seq_along(expr$keys), function(i)
      paste0(.pv_any_str(expr$keys[[i]], TRUE), ": ", .pv_any_str(expr$values[[i]], TRUE)), ""), collapse = ", "), "}"))
  }
  if (is.character(expr) && length(expr) == 1L) return(if (repr) paste0("'", expr, "'") else expr)
  if (is.list(expr)) return(paste0("[", paste(vapply(expr, .pv_any_str, "", repr = TRUE), collapse = ", "), "]"))
  .pv_str(expr)
}

#' Printing for classical mechanics
#'
#' `vsprint()` is SymPy's `vsprint()`: the string form of `expr` in which a
#' function of time `q(t)` prints as `q` and its derivatives as `q'`,
#' `q''`. `vsstrrepr()` is the same with strings quoted, `vprint()` prints
#' it, `vlatex()` gives LaTeX with dots for time derivatives and `vpprint()`
#' pretty prints. `init_vprinting()` makes these the default printers for
#' the session, as SymPy's `init_vprinting()`. `VectorLatexPrinter()` is
#' the printer behind `vlatex()`, with SymPy's `doprint()` method.
#'
#' @param expr An expression, vector, dyadic or matrix.
#' @param ... Ignored printer settings, accepted for compatibility.
#' @return A character scalar; `vpprint()` returns it invisibly after
#'   printing, `vprint()` and `init_vprinting()` return `NULL` invisibly and
#'   `VectorLatexPrinter()` returns a printer object.
#' @examples
#' q <- dynamicsymbols("q")
#' N <- ReferenceFrame("N")
#' vsprint(dynamicsymbols("q", 1) * N$x)
#' vlatex(dynamicsymbols("q", 2))
#' VectorLatexPrinter()$doprint(dynamicsymbols("q", 1))
#' @export
vsprint <- function(expr, ...) .pv_with_hook(.pv_vstr_hook, .pv_any_str(expr))

#' @rdname vsprint
#' @export
vsstrrepr <- function(expr, ...) .pv_with_hook(.pv_vstr_hook, .pv_any_str(expr, repr = TRUE))

#' @rdname vsprint
#' @export
vprint <- function(expr, ...) {
  s <- vsprint(expr, ...)
  if (s != "None") cat(s, "\n", sep = "")
  invisible(NULL)
}

# --- LaTeX

.pv_plain_latex <- function(e) {
  old <- .latex_mode$hook
  .latex_mode$hook <- NULL
  on.exit(.latex_mode$hook <- old)
  .latex(e)
}

.pv_latex_base_symbol <- function(e, exp = NULL) {
  base <- .latex(Symbol(e$name))
  if (is.null(exp)) return(base)
  if (grepl("^", base, fixed = TRUE)) base <- paste0("{", base, "}")
  sprintf("%s^{%s}", base, exp)
}

# VectorLatexPrinter._print_Function and _print_Derivative
.pv_vlatex_hook <- function(e) {
  if (e$type == "Function" && .pv_undef_t(e)) return(.pv_latex_base_symbol(e))
  if (e$type == "Pow" && .pv_undef_t(e$args[[1L]])) {
    x <- e$args[[2L]]
    if (!(is_number_type(x) && num_sign(x) < 0) && !identical(x$key, S_Half$key))
      return(.pv_latex_base_symbol(e$args[[1L]], .latex(x)))
  }
  if (e$type == "Derivative") {
    d <- doit(e)
    if (d$type != "Derivative") return(paste0("\\left(", .latex(d), "\\right)"))
    vars <- .pv_deriv_vars(d)
    tk <- .pv_t()$key
    if (!all(vapply(vars, function(v) identical(v$key, tk), TRUE))) return(NULL)
    f <- d$args[[1L]]
    if (!.is_undefined_application(f)) return(NULL)
    dots <- length(vars)
    base <- if (.pv_undef_t(f)) .pv_latex_base_symbol(f) else .pv_plain_latex(f)
    sp <- regmatches(base, regexpr("_", base, fixed = TRUE), invert = TRUE)[[1L]]
    b <- sp[1L]
    acc <- switch(as.character(dots), "1" = "\\dot", "2" = "\\ddot", "3" = "\\dddot", "4" = "\\ddddot", NULL)
    if (is.null(acc)) return(NULL)
    b <- sprintf("%s{%s}", acc, b)
    if (length(sp) != 1L) b <- paste0(b, "_", sp[2L])
    return(b)
  }
  NULL
}

.pv_vec_latex <- function(v) {
  ar <- .pv_vargs(v)
  if (!length(ar)) return("0")
  ol <- character()
  for (a in ar) {
    lv <- get("latex_vecs", envir = a[[2L]])
    for (j in 1:3) {
      c <- .pv_e(a[[1L]], j)
      if (identical(c$key, S_One$key)) ol <- c(ol, paste0(" + ", lv[j]))
      else if (identical(c$key, S_NegativeOne$key)) ol <- c(ol, paste0(" - ", lv[j]))
      else if (!identical(c$key, S_Zero$key)) {
        s <- .latex(c)
        if (c$type == "Add") s <- paste0("(", s, ")")
        if (startsWith(s, "-")) ol <- c(ol, paste0(" - ", substring(s, 2L), lv[j]))
        else ol <- c(ol, paste0(" + ", s, lv[j]))
      }
    }
  }
  out <- paste(ol, collapse = "")
  if (startsWith(out, " + ")) substring(out, 4L) else if (startsWith(out, " ")) substring(out, 2L) else out
}

.pv_dyad_latex <- function(d) {
  ar <- .pv_dargs(d)
  if (!length(ar)) return("0")
  ol <- character()
  for (v in ar) {
    pv <- paste0(.pv_vec_latex(v[[2L]]), "\\otimes ", .pv_vec_latex(v[[3L]]))
    if (identical(v[[1L]]$key, S_One$key)) ol <- c(ol, paste0(" + ", pv))
    else if (identical(v[[1L]]$key, S_NegativeOne$key)) ol <- c(ol, paste0(" - ", pv))
    else if (!identical(v[[1L]]$key, S_Zero$key)) {
      s <- .latex(v[[1L]])
      if (v[[1L]]$type == "Add") s <- paste0("(", s, ")")
      if (startsWith(s, "-")) ol <- c(ol, paste0(" - ", substring(s, 2L), pv))
      else ol <- c(ol, paste0(" + ", s, pv))
    }
  }
  out <- paste(ol, collapse = "")
  if (startsWith(out, " + ")) substring(out, 4L) else if (startsWith(out, " ")) substring(out, 2L) else out
}

.pv_any_latex <- function(expr) {
  if (.pv_is_vec(expr)) return(.pv_vec_latex(expr))
  if (.pv_is_dyad(expr)) return(.pv_dyad_latex(expr))
  if (inherits(expr, "physics_obj")) return(format(expr))
  latex(expr)
}

#' @rdname vsprint
#' @export
vlatex <- function(expr, ...) {
  old <- .latex_mode$hook
  .latex_mode$hook <- .pv_vlatex_hook
  on.exit(.latex_mode$hook <- old)
  .pv_any_latex(expr)
}

#' @rdname vsprint
#' @export
VectorLatexPrinter <- function(...) .pm_new_obj("VectorLatexPrinter")
.pvM_VectorLatexPrinter_doprint <- function(self, expr) vlatex(expr)

# --- pretty printing

.pv_vpp_hook <- function(e, u) {
  if (e$type == "Function" && .pv_undef_t(e)) return(.bx(e$name))
  if (e$type == "Derivative") {
    vars <- .pv_deriv_vars(e)
    tk <- .pv_t()$key
    f <- e$args[[1L]]
    if (!all(vapply(vars, function(v) identical(v$key, tk), TRUE)) || !.pv_undef_t(f) || length(vars) >= 5L) return(NULL)
    if (!u) return(.bx(paste0(f$name, strrep("'", length(vars)), "(t)")))
    dots <- c("\u0307", "\u0308", "\u20db", "\u20dc")[length(vars)]
    nm <- f$name
    # center_accent: the accent goes after the middle character
    mid <- (nchar(nm) + 1L) %/% 2L
    return(.bx(paste0(substr(nm, 1L, mid), dots, substring(nm, mid + 1L))))
  }
  NULL
}

.pv_vec_pretty <- function(v, u) {
  terms <- list()
  for (a in .pv_vargs(v)) {
    pv <- get("pretty_vecs", envir = a[[2L]])
    for (i in 1:3) {
      c <- .pv_e(a[[1L]], i)
      if (identical(c$key, S_Zero$key)) next
      if (identical(c$key, S_One$key)) terms[[length(terms) + 1L]] <- .bx(pv[i])
      else if (identical(c$key, S_NegativeOne$key)) terms[[length(terms) + 1L]] <- .bx(paste0("-", pv[i]))
      else {
        b <- .pp(c, u)
        if (c$type == "Add") b <- .bx_paren(b, u)
        terms[[length(terms) + 1L]] <- .bx_hcat(b, .bx(" "), .bx(pv[i]))
      }
    }
  }
  if (!length(terms)) return(.bx("0"))
  Reduce(function(a, b) .bx_hcat(a, .bx(" + "), b), terms)
}

#' @rdname vsprint
#' @param use_unicode Unicode or ASCII output for `vpprint()`.
#' @export
vpprint <- function(expr, use_unicode = NULL, ...) {
  u <- .pp_unicode(use_unicode)
  old <- .pp_env$hook
  .pp_env$hook <- .pv_vpp_hook
  on.exit(.pp_env$hook <- old)
  b <- if (.pv_is_vec(expr)) .pv_vec_pretty(expr, u) else .pp(as_sym(expr), u)
  s <- paste(b$lines, collapse = "\n")
  cat(s, "\n", sep = "")
  invisible(s)
}

#' @rdname vsprint
#' @export
init_vprinting <- function(...) {
  .str_mode$hook <- .pv_vstr_hook
  .latex_mode$hook <- .pv_vlatex_hook
  .pp_env$hook <- .pv_vpp_hook
  invisible(NULL)
}
