# The LaTeX printer, after sympy/printing/latex.py (LatexPrinter) and the
# traditional precedence of sympy/printing/precedence.py.

.lp_env <- new.env()

.LP_DEFAULTS <- list(
  full_prec = FALSE, fold_frac_powers = FALSE, fold_func_brackets = FALSE,
  fold_short_frac = NULL, inv_trig_style = "abbreviated", itex = FALSE,
  ln_notation = FALSE, long_frac_ratio = NULL, mat_delim = "[", mat_str = NULL,
  mode = "plain", mul_symbol = NULL, order = NULL, symbol_names = list(),
  root_notation = TRUE, mat_symbol_style = "plain", imaginary_unit = "i",
  gothic_re_im = FALSE, decimal_separator = "period", perm_cyclic = TRUE,
  parenthesize_super = TRUE, min = NULL, max = NULL, diff_operator = "d",
  adjoint_style = "dagger", disable_split_super_sub = FALSE,
  open_delim = c("(", ")"), closed_delim = c("[", "]"))

.lp_s <- function(name) .lp_env$s[[name]]

# LatexPrinter.__init__: the settings and the derived separators
.lp_make_settings <- function(user) {
  bad <- setdiff(names(user), names(.LP_DEFAULTS))
  if (length(bad)) .perr("TypeError", sprintf("Unknown setting '%s'.", bad[1L]))
  s <- .LP_DEFAULTS
  for (nm in names(user)) s[nm] <- list(user[[nm]])
  if (!s$mode %in% c("inline", "plain", "equation", "equation*"))
    .perr("ValueError", "'mode' must be one of 'inline', 'plain', 'equation' or 'equation*'")
  if (is.null(s$fold_short_frac) && s$mode == "inline") s$fold_short_frac <- TRUE
  tab <- c(ldot = " \\,.\\, ", dot = " \\cdot ", times = " \\times ")
  ms <- s$mul_symbol
  s$mul_symbol_latex <- if (is.null(ms)) " " else if (ms %in% names(tab)) tab[[ms]] else ms
  msn <- if (is.null(ms)) "dot" else ms
  s$mul_symbol_latex_numbers <- if (msn %in% names(tab)) tab[[msn]]
    else if (trimws(msn) %in% c("", " ", "\\", "\\,", "\\:", "\\;", "\\quad")) tab[["dot"]]
    else msn
  iu <- c(i = "i", ri = "\\mathrm{i}", ti = "\\text{i}", j = "j", rj = "\\mathrm{j}", tj = "\\text{j}")
  s$imaginary_unit_latex <- if (is.null(s$imaginary_unit)) "i" else if (s$imaginary_unit %in% names(iu)) iu[[s$imaginary_unit]] else s$imaginary_unit
  dop <- c(d = "d", rd = "\\mathrm{d}", td = "\\text{d}")
  s$diff_operator_latex <- if (is.null(s$diff_operator)) "d" else if (s$diff_operator %in% names(dop)) dop[[s$diff_operator]] else s$diff_operator
  s
}

#' LaTeX output
#'
#' Converts an expression to its LaTeX form, as SymPy's `latex()`.
#'
#' @param e A `sym` expression, matrix, `Poly` or a plain R value.
#' @param ... Printer settings, as in SymPy: `mul_symbol`, `fold_short_frac`,
#'   `long_frac_ratio`, `fold_frac_powers`, `fold_func_brackets`,
#'   `root_notation`, `inv_trig_style`, `ln_notation`, `mode`, `itex`,
#'   `mat_str`, `mat_delim`, `symbol_names`, `order`, `full_prec`, `min`,
#'   `max`, `imaginary_unit`, `diff_operator`, `gothic_re_im`,
#'   `decimal_separator`, `parenthesize_super`, `mat_symbol_style`,
#'   `adjoint_style`, `disable_split_super_sub`, `perm_cyclic`, `open_delim`
#'   and `closed_delim`.
#' @return A character scalar.
#' @examples
#' x <- Symbol("x")
#' latex(x^2 / 2 + sin(x))
#' latex(sin(x)^2)
#' latex(2 * x / 3, fold_short_frac = TRUE)
#' @export
latex <- function(e, ...) {
  if (inherits(e, "TableForm")) return(.tableform_latex(e))
  old <- .lp_env$s
  on.exit(.lp_env$s <- old)
  .lp_env$s <- .lp_make_settings(list(...))
  tex <- .latex(e)
  mode <- .lp_s("mode")
  if (mode == "plain") tex
  else if (mode == "inline") paste0("$", tex, "$")
  else if (isTRUE(.lp_s("itex"))) paste0("$$", tex, "$$")
  else paste0("\\begin{", mode, "}", tex, "\\end{", mode, "}")
}

# Printer._print: the hook of the physics vector printer, then the printer
# methods of this file, then the printers registered by other modules
# (SymPy's _latex methods), then _print_Basic.
.latex <- function(e, exp = NULL) {
  if (is.null(.lp_env$s)) {
    .lp_env$s <- .lp_make_settings(list())
    on.exit(.lp_env$s <- NULL)
  }
  if (!is.null(.latex_mode$hook) && is_sym(e) && !is.null(r <- .latex_mode$hook(e))) return(r)
  if (!is_sym(e)) return(.lp_nonsym(e))
  t <- e$type
  if (t == "Function") return(.lp_Function(e, exp))
  f <- .lp_printers[[t]]
  if (!is.null(f)) return(f(e, exp))
  f <- .latex_printers[[t]]
  if (!is.null(f)) return(f(e))
  .lp_Basic(e)
}

.lp_nonsym <- function(e) {
  if (is.null(e)) return("\\text{None}")
  if (inherits(e, "symmatrix")) return(.lp_matrix(e))
  if (inherits(e, "Poly")) return(.lp_Poly(e))
  if (inherits(e, "sym_function")) return(.lp_hprint_Function(attr(e, "sym_name")))
  if (inherits(e, "Permutation")) return(.lp_Permutation(e))
  if (inherits(e, "Cycle")) return(.lp_Cycle(e))
  if (inherits(e, "FreeGroupElement")) return(.lp_FreeGroupElement(e))
  if (inherits(e, "Quaternion")) return(.lp_Quaternion(e))
  if (inherits(e, "sym_ordinal")) return(.lp_ordinal(e))
  if (inherits(e, "sym_omega_power")) return(.lp_omega_power(e))
  if (is.function(e)) return(.lp_function_object(e))
  if (is.logical(e) && length(e) == 1L) return(if (e) "\\text{True}" else "\\text{False}")
  if (is.numeric(e) && length(e) == 1L && !is.object(e)) return(.latex(as_sym(e)))
  if (is.list(e) && !is.object(e)) {
    sep <- if (.lp_s("decimal_separator") == "comma") "; \\  " else ", \\  "
    return(paste0("\\left[ ", paste(vapply(e, .latex, ""), collapse = sep), "\\right]"))
  }
  if (is.character(e) && length(e) == 1L) return(paste0("\\mathtt{\\text{", latex_escape(e), "}}"))
  r <- tryCatch(as_sym(e), error = function(err) NULL)
  if (!is.null(r)) return(.latex(r))
  paste0("\\mathtt{\\text{", latex_escape(paste(format(e), collapse = " ")), "}}")
}

#' Escape a string for LaTeX
#'
#' Escapes the characters that LaTeX treats specially, so the string prints
#' as plain text, as SymPy's `latex_escape()`.
#'
#' @param s A character scalar.
#' @return The escaped string.
#' @examples
#' latex_escape("50% of $x_1$")
#' @export
latex_escape <- function(s) {
  s <- gsub("\\", "\\textbackslash", s, fixed = TRUE)
  for (ch in c("&", "%", "$", "#", "_", "{", "}")) s <- gsub(ch, paste0("\\", ch), s, fixed = TRUE)
  s <- gsub("~", "\\textasciitilde", s, fixed = TRUE)
  gsub("^", "\\textasciicircum", s, fixed = TRUE)
}

# ---------------------------------------------------------------------------
# Precedence

.LP_PREC_TRAD <- c(Lambda = 1, Xor = 10, Or = 20, And = 30, Add = 40, Mul = 50, Pow = 60, Not = 100,
                   Integral = 50, Sum = 50, Product = 50, Limit = 50, Derivative = 50, TensorProduct = 50,
                   Transpose = 60, Adjoint = 60, Dot = 49, Cross = 49, Gradient = 49, Divergence = 49,
                   Curl = 49, Laplacian = 49, DotProduct = 49, Union = 10, Intersection = 10,
                   Complement = 10, SymmetricDifference = 10, ProductSet = 10)

# precedence_traditional()
.lp_prec <- function(e) {
  if (!is_sym(e)) return(1000)
  t <- e$type
  if (t == "UnevaluatedExpr") return(.lp_prec(e$args[[1L]]))
  v <- .LP_PREC_TRAD[t]
  if (!is.na(v)) return(v[[1L]])
  if (t == "Function") {
    if (identical(e$name, "Mod")) return(50)
    if (e$name %in% c("Max", "Min")) return(1000)
    return(70)
  }
  if (t == "Piecewise") return(70)
  if (t %in% c("StrictLessThan", "LessThan", "StrictGreaterThan", "GreaterThan")) return(35)
  .prec(e)
}

.lp_parens <- function(s) paste0("\\left(", s, "\\right)")

# LatexPrinter.parenthesize
.lp_paren <- function(item, level, is_neg = FALSE, strict = FALSE) {
  p <- .lp_prec(item)
  if (is_neg && strict) return(.lp_parens(.latex(item)))
  if (p < level || (!strict && p <= level)) .lp_parens(.latex(item)) else .latex(item)
}

.lp_paren_super <- function(s) {
  if (grepl("^", s, fixed = TRUE)) {
    if (isTRUE(.lp_s("parenthesize_super"))) .lp_parens(s) else paste0("{", s, "}")
  } else s
}

.lp_do_exponent <- function(tex, exp) if (!is.null(exp)) sprintf("\\left(%s\\right)^{%s}", tex, exp) else tex
.lp_sup <- function(tex, exp) if (!is.null(exp)) sprintf("%s^{%s}", tex, exp) else tex

# ---------------------------------------------------------------------------
# Class tests

.LP_REL <- c("Equality", "Unequality", "StrictLessThan", "LessThan", "StrictGreaterThan", "GreaterThan")
.LP_BOOL <- c("And", "Or", "Not", "Xor", "Implies", "Equivalent", "ITE", "BooleanAtom", "Nand", "Nor", "Xnor")
.lp_is_rel <- function(e) e$type %in% .LP_REL
.lp_is_number <- function(e) is_sym(e) && (e$type %in% c("Integer", "Rational", "Float") ||
  (e$type == "Constant" && e$name %in% c("oo", "-oo", "zoo", "nan")))
.lp_is_function <- function(e) is_sym(e) && e$type %in% c("Function", "Piecewise", "And", "Or", "Not", "Xor")
.lp_is_atom <- function(e) is_sym(e) && (length(e$args) == 0L || e$type %in% c("Integer", "Rational", "Float"))
.lp_is_mul <- function(e) e$type %in% c("Mul", "MatMul")
.lp_is_unit <- function(e) {
  is_sym(e) && (e$type %in% c("Quantity", "Prefix", "PhysicalConstant") ||
                  (e$type == "Pow" && e$args[[1L]]$type == "PhysicalConstant"))
}
.lp_has <- function(e, pred) {
  if (!is_sym(e)) return(FALSE)
  if (pred(e)) return(TRUE)
  for (a in e$args) if (.lp_has(a, pred)) return(TRUE)
  FALSE
}
.lp_has_mod <- function(e) .lp_has(e, function(z) z$type == "Function" && identical(z$name, "Mod"))
.lp_cems <- function(e) {
  if (is_sym(e) && e$type == "MatMul") {
    a1 <- e$args[[1L]]
    return(is_exact_number(a1) && num_sign(a1) < 0)
  }
  isTRUE(tryCatch(could_extract_minus_sign(e), error = function(err) FALSE))
}
.lp_dsk <- function(xs) if (length(xs) < 2L) xs else .key_sort(lapply(xs, .sort_key), xs)

.lp_needs_brackets <- function(e) {
  !((e$type == "Integer" && e$value >= 0) ||
      (.lp_is_atom(e) && !identical(e$key, S_NegativeOne$key) && !(e$type %in% c("Integer", "Rational"))))
}

.lp_needs_function_brackets <- function(e) {
  if (!.lp_needs_brackets(e)) return(FALSE)
  if (e$type == "Mul" && !all(!vapply(e$args, .lp_is_function, NA))) return(TRUE)
  if (e$type == "Pow" && .lp_needs_brackets(e$args[[1L]])) return(TRUE)
  e$type == "Add" || .lp_is_function(e)
}

.lp_needs_mul_brackets <- function(e, first = FALSE, last = FALSE) {
  if (e$type == "UnevaluatedExpr") {
    a <- e$args[[1L]]
    if (a$type == "Mul" && .lp_cems(a)) return(TRUE)
  }
  if (.lp_is_mul(e)) {
    if (!first && .lp_cems(e)) return(TRUE)
  } else if (.lp_prec(e) < 50) {
    return(TRUE)
  } else if (.lp_is_rel(e)) {
    return(TRUE)
  }
  if (e$type == "Piecewise") return(TRUE)
  if (.lp_has_mod(e)) return(TRUE)
  if (!last && .lp_has(e, function(z) z$type %in% c("Integral", "Product", "Sum"))) return(TRUE)
  FALSE
}

.lp_needs_add_brackets <- function(e) {
  if (e$type == "UnevaluatedExpr" && e$args[[1L]]$type == "Add") return(TRUE)
  if (.lp_is_rel(e)) return(TRUE)
  if (.lp_has_mod(e)) return(TRUE)
  e$type == "Add"
}

# ---------------------------------------------------------------------------
# Atoms

.lp_Integer <- function(e, exp = NULL) as.character(gmp::numerator(e$value))

.lp_Rational <- function(e, exp = NULL) {
  p <- gmp::numerator(e$value); q <- gmp::denominator(e$value)
  if (q == 1) return(as.character(p))
  sign <- if (p < 0) "- " else ""
  p <- abs(p)
  if (isTRUE(.lp_s("fold_short_frac"))) sprintf("%s%s / %s", sign, as.character(p), as.character(q))
  else sprintf("%s\\frac{%s}{%s}", sign, as.character(p), as.character(q))
}

.lp_Float <- function(e, exp = NULL) {
  v <- e$value
  if (is.nan(v)) return("\\text{NaN}")
  if (is.null(e$extra$mp) && is.infinite(v)) return(if (v > 0) "\\infty" else "- \\infty")
  prec <- .float_prec(e)
  dps <- .prec_to_dps(prec)
  s <- .mp_to_str(.float_mp(e), dps, strip_zeros = !isTRUE(.lp_s("full_prec")),
                  min_fixed = .lp_s("min"), max_fixed = .lp_s("max"))
  sep <- .lp_s("mul_symbol_latex_numbers")
  if (grepl("e", s, fixed = TRUE)) {
    p <- strsplit(s, "e", fixed = TRUE)[[1L]]
    mant <- p[1L]; ex <- sub("^\\+", "", p[2L])
    if (.lp_s("decimal_separator") == "comma") mant <- gsub(".", "{,}", mant, fixed = TRUE)
    return(sprintf("%s%s10^{%s}", mant, sep, ex))
  }
  if (.lp_s("decimal_separator") == "comma") s <- gsub(".", "{,}", s, fixed = TRUE)
  s
}

.lp_Constant <- function(e, exp = NULL) {
  switch(EXPR = e$name, pi = "\\pi", E = "e", I = .lp_s("imaginary_unit_latex"), oo = "\\infty",
         `-oo` = "-\\infty", zoo = "\\tilde{\\infty}", nan = "\\text{NaN}", EulerGamma = "\\gamma",
         GoldenRatio = "\\phi", TribonacciConstant = "\\text{TribonacciConstant}", Catalan = "G",
         .lp_symbol_name(e$name))
}

.lp_BooleanAtom <- function(e, exp = NULL) if (e$name == "true") "\\text{True}" else "\\text{False}"

# ---------------------------------------------------------------------------
# Symbols: split_super_sub, translate and the modifiers

.LP_GREEKS <- c("alpha", "beta", "gamma", "delta", "epsilon", "zeta", "eta", "theta", "iota", "kappa",
                "lambda", "mu", "nu", "xi", "omicron", "pi", "rho", "sigma", "tau", "upsilon", "phi",
                "chi", "psi", "omega")

.LP_TEX_GREEK <- c(Alpha = "\\mathrm{A}", Beta = "\\mathrm{B}", Gamma = "\\Gamma", Delta = "\\Delta",
                   Epsilon = "\\mathrm{E}", Zeta = "\\mathrm{Z}", Eta = "\\mathrm{H}", Theta = "\\Theta",
                   Iota = "\\mathrm{I}", Kappa = "\\mathrm{K}", Lambda = "\\Lambda", Mu = "\\mathrm{M}",
                   Nu = "\\mathrm{N}", Xi = "\\Xi", omicron = "o", Omicron = "\\mathrm{O}", Pi = "\\Pi",
                   Rho = "\\mathrm{P}", Sigma = "\\Sigma", Tau = "\\mathrm{T}", Upsilon = "\\Upsilon",
                   Phi = "\\Phi", Chi = "\\mathrm{X}", Psi = "\\Psi", Omega = "\\Omega",
                   lamda = "\\lambda", Lamda = "\\Lambda", khi = "\\chi", Khi = "\\mathrm{X}",
                   varepsilon = "\\varepsilon", varkappa = "\\varkappa", varphi = "\\varphi",
                   varpi = "\\varpi", varrho = "\\varrho", varsigma = "\\varsigma", vartheta = "\\vartheta")

.LP_OTHER_SYMBOLS <- c("aleph", "beth", "daleth", "gimel", "ell", "eth", "hbar", "hslash", "mho", "wp")

.LP_MODIFIERS <- list(
  mathring = function(s) paste0("\\mathring{", s, "}"), ddddot = function(s) paste0("\\ddddot{", s, "}"),
  dddot = function(s) paste0("\\dddot{", s, "}"), ddot = function(s) paste0("\\ddot{", s, "}"),
  dot = function(s) paste0("\\dot{", s, "}"), check = function(s) paste0("\\check{", s, "}"),
  breve = function(s) paste0("\\breve{", s, "}"), acute = function(s) paste0("\\acute{", s, "}"),
  grave = function(s) paste0("\\grave{", s, "}"), tilde = function(s) paste0("\\tilde{", s, "}"),
  hat = function(s) paste0("\\hat{", s, "}"), bar = function(s) paste0("\\bar{", s, "}"),
  vec = function(s) paste0("\\vec{", s, "}"), prime = function(s) paste0("{", s, "}'"),
  prm = function(s) paste0("{", s, "}'"), bold = function(s) paste0("\\boldsymbol{", s, "}"),
  bm = function(s) paste0("\\boldsymbol{", s, "}"), cal = function(s) paste0("\\mathcal{", s, "}"),
  scr = function(s) paste0("\\mathscr{", s, "}"), frak = function(s) paste0("\\mathfrak{", s, "}"),
  norm = function(s) paste0("\\left\\|{", s, "}\\right\\|"),
  avg = function(s) paste0("\\left\\langle{", s, "}\\right\\rangle"),
  abs = function(s) paste0("\\left|{", s, "}\\right|"), mag = function(s) paste0("\\left|{", s, "}\\right|"))

# the modifier keys longest first, ties in insertion order (a stable sort)
.LP_MOD_KEYS <- names(.LP_MODIFIERS)[order(-nchar(names(.LP_MODIFIERS)))]

#' Translate a symbol name to LaTeX
#'
#' Converts Greek letter names, other special names and trailing modifiers
#' such as `hat` or `dot` to LaTeX, as SymPy's `translate()`.
#'
#' @param s A character scalar.
#' @return The LaTeX string.
#' @examples
#' latex_translate("alphahatdotprime")
#' @export
latex_translate <- function(s) {
  tex <- .LP_TEX_GREEK[s]
  if (!is.na(tex) && nzchar(tex)) return(tex[[1L]])
  lo <- tolower(s)
  if (lo %in% .LP_GREEKS) return(paste0("\\", lo))
  if (s %in% .LP_OTHER_SYMBOLS) return(paste0("\\", s))
  for (key in .LP_MOD_KEYS) {
    if (endsWith(lo, key) && nchar(s) > nchar(key))
      return(.LP_MODIFIERS[[key]](latex_translate(substr(s, 1L, nchar(s) - nchar(key)))))
  }
  s
}

.lp_split_super_sub <- function(name) {
  if (is.null(name) || grepl("{", name, fixed = TRUE)) return(list(name = name, supers = character(), subs = character()))
  if (isTRUE(.lp_s("disable_split_super_sub"))) {
    p <- list(name = gsub("^", "\\^", gsub("_", "\\_", name, fixed = TRUE), fixed = TRUE), supers = character(), subs = character())
  } else p <- .split_super_sub(name)
  list(name = latex_translate(p$name), supers = vapply(p$supers, latex_translate, "", USE.NAMES = FALSE),
       subs = vapply(p$subs, latex_translate, "", USE.NAMES = FALSE))
}

.lp_deal_with_super_sub <- function(string, style = "plain") {
  p <- .lp_split_super_sub(string)
  name <- p$name
  if (length(p$name) == 0L) return(string)
  if (style == "bold") name <- paste0("\\mathbf{", name, "}")
  if (length(p$supers)) name <- paste0(name, "^{", paste(p$supers, collapse = " "), "}")
  if (length(p$subs)) name <- paste0(name, "_{", paste(p$subs, collapse = " "), "}")
  name
}

.lp_symbol_name <- function(nm) .lp_deal_with_super_sub(nm)

.lp_Symbol <- function(e, exp = NULL, style = "plain") {
  sn <- .lp_s("symbol_names")
  if (length(sn)) {
    nm <- sn[[e$name]]
    if (!is.null(nm)) return(nm)
  }
  if (!is.null(e$extra$idx)) return(.lp_Idx(e))
  .lp_deal_with_super_sub(e$name, style = style)
}

.lp_Idx <- function(e) {
  label <- .lp_deal_with_super_sub(e$name)
  rg <- e$extra$range
  if (is.null(rg) || length(rg) < 2L || is.null(rg[[2L]])) return(label)
  lower <- if (is.null(rg[[1L]])) "0" else .latex(as_sym(rg[[1L]]))
  upper <- .latex(as_sym(rg[[2L]]))
  sprintf("{%s}_{%s\\mathrel{..}%s}", label, lower, upper)
}

# ---------------------------------------------------------------------------
# Add, Mul, Pow

.lp_ordered_terms <- function(e) {
  if (identical(.lp_s("order"), "none")) return(.sympy_args(e))
  is_o <- vapply(e$args, function(t) t$type == "Order", TRUE)
  if (any(is_o)) return(c(.add_print_order(e$args[!is_o], rev = TRUE), .add_print_order(e$args[is_o], rev = TRUE)))
  if (any(vapply(e$args, .me_has_nc_factor, NA))) return(.me_ordered_terms(e$args))
  if (identical(.lp_s("order"), "old")) return(rev(.add_print_order(e$args)))
  .add_print_order(e$args)
}

.lp_neg <- function(t) {
  if (t$type == "MatMul") {
    a <- t$args
    if (identical(a[[1L]]$key, S_NegativeOne$key)) a <- a[-1L] else a[[1L]] <- .neg(a[[1L]])
    if (length(a) == 1L) return(a[[1L]])
    return(new_sym("MatMul", a))
  }
  if (t$type == "Mul") .neg_mul_noflat(t) else .neg(t)
}

.lp_Add <- function(e, exp = NULL) {
  terms <- .lp_ordered_terms(e)
  tex <- ""
  for (i in seq_along(terms)) {
    term <- terms[[i]]
    if (i > 1L) {
      if (.lp_cems(term)) {
        tex <- paste0(tex, " - ")
        term <- .lp_neg(term)
      } else tex <- paste0(tex, " + ")
    }
    tt <- .latex(term)
    if (.lp_needs_add_brackets(term)) tt <- .lp_parens(tt)
    tex <- paste0(tex, tt)
  }
  tex
}

.lp_mul_factors <- function(e) {
  if (!(.lp_is_mul(e) || (e$type == "UnevaluatedExpr" && e$args[[1L]]$type == "Mul"))) return(NULL)
  if (e$type != "Mul") return(list(e))
  if (identical(.lp_s("order"), "old") || identical(.lp_s("order"), "none")) .sympy_args(e)
  else .mul_print_order(e$args)
}

.lp_convert_args <- function(args) {
  sep <- .lp_s("mul_symbol_latex"); numsep <- .lp_s("mul_symbol_latex_numbers")
  tex <- last <- ""
  n <- length(args)
  for (i in seq_len(n)) {
    term <- args[[i]]
    tt <- .latex(term)
    if (!.lp_is_unit(term)) {
      if (.lp_needs_mul_brackets(term, first = i == 1L, last = i == n)) tt <- .lp_parens(tt)
      if (grepl("[0-9][} ]*$", last) && grepl("^([0-9]|\\\\frac\\{[0-9]+\\}\\{[0-9]+\\})", tt))
        tex <- paste0(tex, numsep)
      else if (nzchar(tex)) tex <- paste0(tex, sep)
    } else if (nzchar(tex)) tex <- paste0(tex, sep)
    tex <- paste0(tex, tt)
    last <- tt
  }
  tex
}

.lp_convert <- function(e) {
  fs <- .lp_mul_factors(e)
  if (is.null(fs)) return(.latex(e))
  units <- vapply(fs, .lp_is_unit, NA)
  pre <- vapply(fs, function(z) z$type == "Prefix", NA)
  .lp_convert_args(c(fs[!units], fs[units & pre], fs[units & !pre]))
}

.lp_mul2 <- function(a, b) .mul(list(a, b))

.lp_Mul <- function(e, exp = NULL) {
  sep <- .lp_s("mul_symbol_latex")
  is_add <- function(z) z$type == "Add" || (z$type == "UnevaluatedExpr" && z$args[[1L]]$type == "Add")
  if (e$type == "Mul") {
    args <- e$args
    if (is_one_num(args[[1L]]) || any(vapply(args[-1L], .lp_is_number, NA))) return(.lp_convert_args(args))
  }
  include_parens <- FALSE
  expr <- e
  if (.lp_cems(expr)) {
    expr <- .lp_neg(expr)
    tex <- "- "
    if (is_add(expr) || (expr$type == "UnevaluatedExpr" && .lp_is_number(expr$args[[1L]]) && num_sign(expr$args[[1L]]) < 0)) {
      tex <- paste0(tex, "\\left(")
      include_parens <- TRUE
    }
  } else tex <- ""
  nd <- fraction(expr, exact = TRUE)
  numer <- nd[[1L]]; denom <- nd[[2L]]
  one_inv <- expr$type == "Mul" && any(vapply(expr$args, function(a) a$type == "Pow" && is_one_num(a$args[[1L]]) &&
                                                identical(a$args[[2L]]$key, S_NegativeOne$key), NA))
  if (is_one_num(denom) && !one_inv) {
    tex <- paste0(tex, .lp_convert(expr))
  } else {
    snumer <- .lp_convert(numer)
    sdenom <- .lp_convert(denom)
    ntok <- function(s) length(strsplit(trimws(s), "[[:space:]]+")[[1L]])
    ldenom <- ntok(sdenom)
    ratio <- .lp_s("long_frac_ratio")
    if (isTRUE(.lp_s("fold_short_frac")) && ldenom <= 2 && !grepl("^", sdenom, fixed = TRUE)) {
      if (.lp_needs_mul_brackets(numer, last = FALSE)) tex <- paste0(tex, sprintf("\\left(%s\\right) / %s", snumer, sdenom))
      else tex <- paste0(tex, sprintf("%s / %s", snumer, sdenom))
    } else if (!is.null(ratio) && ntok(snumer) > ratio * ldenom) {
      if (.lp_needs_mul_brackets(numer, last = TRUE)) {
        tex <- paste0(tex, sprintf("\\frac{1}{%s}%s\\left(%s\\right)", sdenom, sep, snumer))
      } else if (numer$type == "Mul") {
        a <- S_One; b <- S_One
        for (xx in .sympy_args(numer)) {
          if (.lp_needs_mul_brackets(xx, last = FALSE) || ntok(.lp_convert(.lp_mul2(a, xx))) > ratio * ldenom ||
              (isFALSE(.is_comm(b)) && isFALSE(.is_comm(xx)))) b <- .lp_mul2(b, xx)
          else a <- .lp_mul2(a, xx)
        }
        if (.lp_needs_mul_brackets(b, last = TRUE))
          tex <- paste0(tex, sprintf("\\frac{%s}{%s}%s\\left(%s\\right)", .lp_convert(a), sdenom, sep, .lp_convert(b)))
        else tex <- paste0(tex, sprintf("\\frac{%s}{%s}%s%s", .lp_convert(a), sdenom, sep, .lp_convert(b)))
      } else {
        tex <- paste0(tex, sprintf("\\frac{1}{%s}%s%s", sdenom, sep, snumer))
      }
    } else {
      tex <- paste0(tex, sprintf("\\frac{%s}{%s}", snumer, sdenom))
    }
  }
  if (include_parens) tex <- paste0(tex, "\\right)")
  tex
}

.lp_helper_std_power <- function(base_e, exp_e, template) {
  ex <- .latex(exp_e)
  if (exp_e$type == "Pow") ex <- .lp_paren_super(ex)
  base <- .lp_paren(base_e, 60)
  if (base_e$type %in% c("Symbol", "Dummy", "Wild")) base <- .lp_paren_super(base)
  else if (base_e$type == "Float") base <- paste0("{", base, "}")
  else if (base_e$type == "Derivative" && startsWith(base, "\\left(") && grepl("^\\\\left\\(\\\\d?d?dot", base) &&
           endsWith(base, "\\right)")) base <- substr(base, 7L, nchar(base) - 7L)
  sprintf(template, base, ex)
}

.lp_Pow <- function(e, exp = NULL) {
  b <- e$args[[1L]]; x <- e$args[[2L]]
  if (x$type %in% c("Integer", "Rational")) {
    p <- gmp::numerator(x$value); q <- gmp::denominator(x$value)
    if (abs(p) == 1 && q != 1 && isTRUE(.lp_s("root_notation"))) {
      base <- .latex(b)
      tex <- if (q == 2) sprintf("\\sqrt{%s}", base)
             else if (isTRUE(.lp_s("itex"))) sprintf("\\root{%s}{%s}", as.character(q), base)
             else sprintf("\\sqrt[%s]{%s}", as.character(q), base)
      return(if (p < 0) sprintf("\\frac{1}{%s}", tex) else tex)
    } else if (isTRUE(.lp_s("fold_frac_powers")) && q != 1) {
      base <- .lp_paren(b, 60)
      if (b$type %in% c("Symbol", "Dummy")) base <- .lp_paren_super(base)
      if (.lp_is_function(b)) return(.latex(b, exp = sprintf("%s/%s", as.character(p), as.character(q))))
      return(sprintf("%s^{%s/%s}", base, as.character(p), as.character(q)))
    } else if (p < 0 && .is_comm(b)) {
      if (is_one_num(b)) {
        if (identical(x$key, S_NegativeOne$key)) return("\\frac{1}{1}")
        return(sprintf("%s^{%s}", .str(b), .str(x)))
      }
      if (b$type %in% c("Integer", "Rational")) {
        bp <- gmp::numerator(b$value); bq <- gmp::denominator(b$value)
        if (bp * bq == abs(bq)) {
          if (identical(x$key, S_NegativeOne$key)) return(sprintf("\\frac{1}{\\frac{%s}{%s}}", as.character(bp), as.character(bq)))
          return(sprintf("\\frac{1}{(\\frac{%s}{%s})^{%s}}", as.character(bp), as.character(bq), .str(.neg(x))))
        }
      }
      if (identical(x$key, S_NegativeOne$key)) {
        if (isTRUE(.lp_s("fold_short_frac"))) return(sprintf("1 / %s", .lp_paren(b, 60)))
        return(sprintf("\\frac{1}{%s}", .latex(b)))
      }
      tex <- .lp_helper_std_power(b, .neg(x), "%s^{%s}")
      return(sprintf("\\frac{1}{%s}", tex))
    }
  }
  if (.lp_is_function(b)) return(.latex(b, exp = .latex(x)))
  .lp_helper_std_power(b, x, "%s^{%s}")
}

# ---------------------------------------------------------------------------
# Functions

.LP_ACCEPTED <- c("arcsin", "arccos", "arctan", "sin", "cos", "tan", "sinh", "cosh", "tanh", "sqrt", "ln",
                  "log", "sec", "csc", "cot", "coth", "re", "im", "frac", "root", "arg", "exp")

.LP_INV_TRIG <- c("asin", "acos", "atan", "acsc", "asec", "acot", "asinh", "acosh", "atanh", "acsch", "asech", "acoth")

# _hprint_Function: a LaTeX command, a single letter, or \operatorname
.lp_hprint_Function <- function(func) {
  func <- .lp_deal_with_super_sub(func)
  sup <- regexpr("^", func, fixed = TRUE) - 1L
  sub <- regexpr("_", func, fixed = TRUE) - 1L
  if (func %in% .LP_ACCEPTED) return(paste0("\\", func))
  if (nchar(func) == 1L || startsWith(func, "\\") || sub == 1L || sup == 1L) return(func)
  cut <- if (sup > 0L && sub > 0L) min(sub, sup) else if (sup > 0L) sup else if (sub > 0L) sub else -1L
  if (cut > 0L) sprintf("\\operatorname{%s}%s", substr(func, 1L, cut), substring(func, cut + 1L))
  else sprintf("\\operatorname{%s}", func)
}

.lp_Function <- function(e, exp = NULL) {
  nm <- e$name
  d <- .FUNS[[nm]]
  f <- .lp_funcs[[nm]]
  if (!is.null(f)) return(f(e, exp))
  if (!is.null(d) && !is.null(d$latex)) {
    tex <- d$latex(vapply(e$args, .latex, ""))
    return(.lp_do_exponent(tex, exp))
  }
  .lp_generic_function(e, nm, exp)
}

# LatexPrinter._print_Function for a function without its own method
.lp_generic_function <- function(e, func, exp = NULL) {
  args <- vapply(e$args, .latex, "")
  style <- .lp_s("inv_trig_style")
  power_case <- FALSE
  can_fold <- isTRUE(.lp_s("fold_func_brackets")) && length(args) == 1L && !.lp_needs_function_brackets(e$args[[1L]])
  if (func %in% .LP_INV_TRIG) {
    if (style == "full") func <- paste0(if (endsWith(func, "h")) "ar" else "arc", substring(func, 2L))
    else if (style == "power") {
      func <- substring(func, 2L)
      power_case <- TRUE
      if (!is.null(exp)) can_fold <- FALSE
    }
  }
  if (power_case) {
    name <- if (func %in% .LP_ACCEPTED) sprintf("\\%s^{-1}", func) else sprintf("\\operatorname{%s}^{-1}", func)
  } else if (!is.null(exp)) {
    name <- sprintf("%s^{%s}", .lp_paren_super(.lp_hprint_Function(func)), exp)
  } else name <- .lp_hprint_Function(func)
  a <- paste(args, collapse = ",")
  name <- if (can_fold) {
    if (func %in% .LP_ACCEPTED) paste0(name, " {", a, "}") else paste0(name, a)
  } else paste0(name, "{\\left(", a, " \\right)}")
  if (power_case && !is.null(exp)) name <- paste0(name, "^{", exp, "}")
  name
}

.lp_arg <- function(e, i = 1L) .latex(e$args[[i]])

# printers of the form NAME\left(a\right) with the exponent after NAME
.lp_one_arg <- function(sym) function(e, exp = NULL) {
  tex <- sprintf("\\left(%s\\right)", .lp_arg(e))
  if (!is.null(exp)) sprintf("%s^{%s}%s", sym, exp, tex) else paste0(sym, tex)
}

# printers whose power wraps the whole form in parentheses
.lp_wrap_exp <- function(fmt) function(e, exp = NULL) {
  tex <- do.call(sprintf, c(list(fmt), as.list(vapply(e$args, .latex, ""))))
  if (!is.null(exp)) tex <- paste0("\\left(", tex, "\\right)^{", exp, "}")
  tex
}

.lp_bessel <- function(sym) function(e, exp = NULL) {
  tex <- sym
  need_exp <- FALSE
  if (!is.null(exp)) {
    if (!grepl("^", tex, fixed = TRUE)) tex <- sprintf("%s^{%s}", tex, exp) else need_exp <- TRUE
  }
  tex <- sprintf("%s_{%s}\\left(%s\\right)", tex, .lp_arg(e, 1L), .lp_arg(e, 2L))
  if (need_exp) tex <- .lp_do_exponent(tex, exp)
  tex
}

.lp_airy <- function(notation, prime = FALSE) function(e, exp = NULL) {
  tex <- sprintf("\\left(%s\\right)", .lp_arg(e))
  if (prime) {
    if (!is.null(exp)) sprintf("{%s^\\prime}^{%s}%s", notation, exp, tex) else sprintf("%s^\\prime%s", notation, tex)
  } else {
    if (!is.null(exp)) sprintf("%s^{%s}%s", notation, exp, tex) else paste0(notation, tex)
  }
}

.lp_hprint_vec <- function(v) paste(vapply(v, .latex, ""), collapse = ", ")

.lp_number_poly <- function(letter) function(e, exp = NULL) {
  if (length(e$args) == 2L) {
    if (!is.null(exp)) return(sprintf("%s_{%s}^{%s}\\left(%s\\right)", letter, .lp_arg(e, 1L), exp, .lp_arg(e, 2L)))
    return(sprintf("%s_{%s}\\left(%s\\right)", letter, .lp_arg(e, 1L), .lp_arg(e, 2L)))
  }
  .lp_sup(sprintf("%s_{%s}", letter, .lp_arg(e, 1L)), exp)
}

.lp_mathieu <- function(ch, prime = FALSE) function(e, exp = NULL) {
  a <- vapply(e$args, .latex, "")
  sprintf("%s%s\\left(%s, %s, %s\\right)%s", ch, if (prime) "^{\\prime}" else "", a[1L], a[2L], a[3L],
          if (is.null(exp) || !nzchar(exp)) "" else sprintf("^{%s}", exp))
}

.lp_ntheory1 <- function(sym) function(e, exp = NULL) {
  if (!is.null(exp)) sprintf("\\left(%s\\left(%s\\right)\\right)^{%s}", sym, .lp_arg(e), exp)
  else sprintf("%s\\left(%s\\right)", sym, .lp_arg(e))
}

.lp_divisor_sigma <- function(star) function(e, exp = NULL) {
  tex <- if (length(e$args) == 2L) sprintf("_%s\\left(%s\\right)", .lp_arg(e, 2L), .lp_arg(e, 1L))
         else sprintf("\\left(%s\\right)", .lp_arg(e))
  s <- if (star) "\\sigma^*" else "\\sigma"
  if (!is.null(exp)) sprintf("%s^{%s}%s", s, exp, tex) else paste0(s, tex)
}

.lp_transform <- function(s, inverse = FALSE) function(e, exp = NULL)
  sprintf("\\mathcal{%s}%s_{%s}\\left[%s\\right]\\left(%s\\right)", s, if (inverse) "^{-1}" else "",
          .lp_arg(e, 2L), .lp_arg(e, 1L), .lp_arg(e, 3L))

.lp_tuple_args <- function(a) if (is_sym(a) && a$type == "Tuple") a$args else if (is.list(a) && !is_sym(a)) a else list(a)

.lp_funcs <- list(
  Abs = function(e, exp = NULL) .lp_sup(sprintf("\\left|{%s}\\right|", .lp_arg(e)), exp),
  re = function(e, exp = NULL) {
    a <- .lp_paren(e$args[[1L]], 1000)
    .lp_do_exponent(if (isTRUE(.lp_s("gothic_re_im"))) paste0("\\Re{", a, "}") else paste0("\\operatorname{re}{", a, "}"), exp)
  },
  im = function(e, exp = NULL) {
    a <- .lp_paren(e$args[[1L]], 1000)
    .lp_do_exponent(if (isTRUE(.lp_s("gothic_re_im"))) paste0("\\Im{", a, "}") else paste0("\\operatorname{im}{", a, "}"), exp)
  },
  conjugate = function(e, exp = NULL) .lp_sup(sprintf("\\overline{%s}", .lp_arg(e)), exp),
  polar_lift = function(e, exp = NULL) {
    arg <- sprintf("{\\left(%s \\right)}", .lp_arg(e))
    if (!is.null(exp)) sprintf("\\operatorname{polar\\_lift}^{%s}%s", exp, arg) else paste0("\\operatorname{polar\\_lift}", arg)
  },
  exp = function(e, exp = NULL) .lp_do_exponent(sprintf("e^{%s}", .lp_arg(e)), exp),
  exp_polar = function(e, exp = NULL) .lp_do_exponent(sprintf("e^{%s}", .lp_arg(e)), exp),
  log = function(e, exp = NULL) {
    tex <- if (length(e$args) == 2L) {
      base <- .lp_arg(e, 2L)
      if (nchar(base) == 1L) sprintf("\\log_%s{\\left(%s \\right)}", base, .lp_arg(e, 1L))
      else sprintf("\\log_{%s}{\\left(%s \\right)}", base, .lp_arg(e, 1L))
    } else if (!isTRUE(.lp_s("ln_notation"))) sprintf("\\log{\\left(%s \\right)}", .lp_arg(e))
    else sprintf("\\ln{\\left(%s \\right)}", .lp_arg(e))
    .lp_sup(tex, exp)
  },
  floor = function(e, exp = NULL) .lp_sup(sprintf("\\left\\lfloor{%s}\\right\\rfloor", .lp_arg(e)), exp),
  ceiling = function(e, exp = NULL) .lp_sup(sprintf("\\left\\lceil{%s}\\right\\rceil", .lp_arg(e)), exp),
  frac = function(e, exp = NULL) {
    if (is.null(exp)) sprintf("\\operatorname{frac}{\\left(%s\\right)}", .lp_arg(e))
    else sprintf("\\operatorname{frac}{\\left(%s\\right)}^{%s}", .lp_arg(e), exp)
  },
  Max = function(e, exp = NULL) .lp_variadic(e, "max", exp),
  Min = function(e, exp = NULL) .lp_variadic(e, "min", exp),
  Mod = function(e, exp = NULL) {
    a <- .lp_paren(e$args[[1L]], 50, strict = TRUE); b <- .lp_paren(e$args[[2L]], 50, strict = TRUE)
    if (!is.null(exp)) sprintf("\\left(%s \\bmod %s\\right)^{%s}", a, b, exp) else sprintf("%s \\bmod %s", a, b)
  },
  factorial = function(e, exp = NULL) .lp_sup(sprintf("%s!", .lp_paren(e$args[[1L]], 70)), exp),
  factorial2 = function(e, exp = NULL) .lp_sup(sprintf("%s!!", .lp_paren(e$args[[1L]], 70)), exp),
  subfactorial = function(e, exp = NULL) {
    tex <- sprintf("!%s", .lp_paren(e$args[[1L]], 70))
    if (!is.null(exp)) sprintf("\\left(%s\\right)^{%s}", tex, exp) else tex
  },
  binomial = function(e, exp = NULL) .lp_sup(sprintf("{\\binom{%s}{%s}}", .lp_arg(e, 1L), .lp_arg(e, 2L)), exp),
  RisingFactorial = function(e, exp = NULL)
    .lp_do_exponent(sprintf("{%s}^{\\left(%s\\right)}", .lp_paren(e$args[[1L]], 70), .lp_arg(e, 2L)), exp),
  FallingFactorial = function(e, exp = NULL)
    .lp_do_exponent(sprintf("{\\left(%s\\right)}_{%s}", .lp_arg(e, 1L), .lp_paren(e$args[[2L]], 70)), exp),
  gamma = .lp_one_arg("\\Gamma"),
  Chi = .lp_one_arg("\\operatorname{Chi}"),
  fresnels = .lp_one_arg("S"), fresnelc = .lp_one_arg("C"),
  dirichlet_eta = .lp_one_arg("\\eta"),
  airyai = .lp_airy("Ai"), airybi = .lp_airy("Bi"),
  airyaiprime = .lp_airy("Ai", TRUE), airybiprime = .lp_airy("Bi", TRUE),
  lowergamma = function(e, exp = NULL) {
    tex <- sprintf("\\left(%s, %s\\right)", .lp_arg(e, 1L), .lp_arg(e, 2L))
    if (!is.null(exp)) sprintf("\\gamma^{%s}%s", exp, tex) else paste0("\\gamma", tex)
  },
  uppergamma = function(e, exp = NULL) {
    tex <- sprintf("\\left(%s, %s\\right)", .lp_arg(e, 1L), .lp_arg(e, 2L))
    if (!is.null(exp)) sprintf("\\Gamma^{%s}%s", exp, tex) else paste0("\\Gamma", tex)
  },
  beta = function(e, exp = NULL) {
    y <- if (length(e$args) == 1L) 1L else 2L
    tex <- sprintf("\\left(%s, %s\\right)", .lp_arg(e, 1L), .lp_arg(e, y))
    if (!is.null(exp)) sprintf("\\operatorname{B}^{%s}%s", exp, tex) else paste0("\\operatorname{B}", tex)
  },
  betainc = function(e, exp = NULL) .lp_betainc(e, exp, "B"),
  betainc_regularized = function(e, exp = NULL) .lp_betainc(e, exp, "I"),
  expint = function(e, exp = NULL) {
    tex <- sprintf("\\left(%s\\right)", .lp_arg(e, 2L)); nu <- .lp_arg(e, 1L)
    if (!is.null(exp)) sprintf("\\operatorname{E}_{%s}^{%s}%s", nu, exp, tex) else sprintf("\\operatorname{E}_{%s}%s", nu, tex)
  },
  owens_t = function(e, exp = NULL) {
    tex <- sprintf("\\left(%s, %s\\right)", .lp_arg(e, 1L), .lp_arg(e, 2L))
    if (!is.null(exp)) sprintf("T^{%s}%s", exp, tex) else paste0("T", tex)
  },
  besselj = .lp_bessel("J"), besseli = .lp_bessel("I"), besselk = .lp_bessel("K"), bessely = .lp_bessel("Y"),
  yn = .lp_bessel("y"), jn = .lp_bessel("j"), hankel1 = .lp_bessel("H^{(1)}"), hankel2 = .lp_bessel("H^{(2)}"),
  hn1 = .lp_bessel("h^{(1)}"), hn2 = .lp_bessel("h^{(2)}"),
  hyper = function(e, exp = NULL) {
    ap <- .lp_tuple_args(e$args[[1L]]); bq <- .lp_tuple_args(e$args[[2L]])
    tex <- sprintf("{{}_{%d}F_{%d}\\left(\\begin{matrix} %s \\\\ %s \\end{matrix}\\middle| {%s} \\right)}",
                   length(ap), length(bq), .lp_hprint_vec(ap), .lp_hprint_vec(bq), .lp_arg(e, 3L))
    if (!is.null(exp)) tex <- sprintf("{%s}^{%s}", tex, exp)
    tex
  },
  meijerg = function(e, exp = NULL) {
    a <- .lp_tuple_args(e$args[[1L]]); b <- .lp_tuple_args(e$args[[2L]])
    an <- .lp_tuple_args(a[[1L]]); aother <- .lp_tuple_args(a[[2L]])
    bm <- .lp_tuple_args(b[[1L]]); bother <- .lp_tuple_args(b[[2L]])
    tex <- sprintf("{G_{%d, %d}^{%d, %d}\\left(\\begin{matrix} %s & %s \\\\%s & %s \\end{matrix} \\middle| {%s} \\right)}",
                   length(an) + length(aother), length(bm) + length(bother), length(bm), length(an),
                   .lp_hprint_vec(an), .lp_hprint_vec(aother), .lp_hprint_vec(bm), .lp_hprint_vec(bother), .lp_arg(e, 3L))
    if (!is.null(exp)) tex <- sprintf("{%s}^{%s}", tex, exp)
    tex
  },
  zeta = function(e, exp = NULL) {
    tex <- if (length(e$args) == 2L) sprintf("\\left(%s, %s\\right)", .lp_arg(e, 1L), .lp_arg(e, 2L))
           else sprintf("\\left(%s\\right)", .lp_arg(e))
    if (!is.null(exp)) sprintf("\\zeta^{%s}%s", exp, tex) else paste0("\\zeta", tex)
  },
  stieltjes = function(e, exp = NULL) {
    tex <- if (length(e$args) == 2L) sprintf("_{%s}\\left(%s\\right)", .lp_arg(e, 1L), .lp_arg(e, 2L))
           else sprintf("_{%s}", .lp_arg(e))
    if (!is.null(exp)) sprintf("\\gamma%s^{%s}", tex, exp) else paste0("\\gamma", tex)
  },
  lerchphi = function(e, exp = NULL) {
    tex <- sprintf("\\left(%s, %s, %s\\right)", .lp_arg(e, 1L), .lp_arg(e, 2L), .lp_arg(e, 3L))
    if (is.null(exp)) paste0("\\Phi", tex) else sprintf("\\Phi^{%s}%s", exp, tex)
  },
  polylog = function(e, exp = NULL) {
    tex <- sprintf("\\left(%s\\right)", .lp_arg(e, 2L))
    if (is.null(exp)) sprintf("\\operatorname{Li}_{%s}%s", .lp_arg(e, 1L), tex)
    else sprintf("\\operatorname{Li}_{%s}^{%s}%s", .lp_arg(e, 1L), exp, tex)
  },
  jacobi = .lp_wrap_exp("P_{%s}^{\\left(%s,%s\\right)}\\left(%s\\right)"),
  gegenbauer = .lp_wrap_exp("C_{%s}^{\\left(%s\\right)}\\left(%s\\right)"),
  chebyshevt = .lp_wrap_exp("T_{%s}\\left(%s\\right)"),
  chebyshevu = .lp_wrap_exp("U_{%s}\\left(%s\\right)"),
  legendre = .lp_wrap_exp("P_{%s}\\left(%s\\right)"),
  assoc_legendre = .lp_wrap_exp("P_{%s}^{\\left(%s\\right)}\\left(%s\\right)"),
  hermite = .lp_wrap_exp("H_{%s}\\left(%s\\right)"),
  laguerre = .lp_wrap_exp("L_{%s}\\left(%s\\right)"),
  assoc_laguerre = .lp_wrap_exp("L_{%s}^{\\left(%s\\right)}\\left(%s\\right)"),
  Ynm = .lp_wrap_exp("Y_{%s}^{%s}\\left(%s,%s\\right)"),
  Znm = .lp_wrap_exp("Z_{%s}^{%s}\\left(%s,%s\\right)"),
  mathieuc = .lp_mathieu("C"), mathieus = .lp_mathieu("S"),
  mathieucprime = .lp_mathieu("C", TRUE), mathieusprime = .lp_mathieu("S", TRUE),
  elliptic_k = function(e, exp = NULL) {
    tex <- sprintf("\\left(%s\\right)", .lp_arg(e))
    if (!is.null(exp)) sprintf("K^{%s}%s", exp, tex) else paste0("K", tex)
  },
  elliptic_f = function(e, exp = NULL) {
    tex <- sprintf("\\left(%s\\middle| %s\\right)", .lp_arg(e, 1L), .lp_arg(e, 2L))
    if (!is.null(exp)) sprintf("F^{%s}%s", exp, tex) else paste0("F", tex)
  },
  elliptic_e = function(e, exp = NULL) {
    tex <- if (length(e$args) == 2L) sprintf("\\left(%s\\middle| %s\\right)", .lp_arg(e, 1L), .lp_arg(e, 2L))
           else sprintf("\\left(%s\\right)", .lp_arg(e))
    if (!is.null(exp)) sprintf("E^{%s}%s", exp, tex) else paste0("E", tex)
  },
  elliptic_pi = function(e, exp = NULL) {
    tex <- if (length(e$args) == 3L) sprintf("\\left(%s; %s\\middle| %s\\right)", .lp_arg(e, 1L), .lp_arg(e, 2L), .lp_arg(e, 3L))
           else sprintf("\\left(%s\\middle| %s\\right)", .lp_arg(e, 1L), .lp_arg(e, 2L))
    if (!is.null(exp)) sprintf("\\Pi^{%s}%s", exp, tex) else paste0("\\Pi", tex)
  },
  jtheta = function(e, exp = NULL) {
    tex <- sprintf("\\vartheta_{%s}", .lp_arg(e, 1L))
    if (length(e$args) == 4L) tex <- paste0(tex, sprintf("^{(%s)}", .lp_arg(e, 4L)))
    tex <- paste0(tex, sprintf("\\left(%s, %s\\right)", .lp_arg(e, 2L), .lp_arg(e, 3L)))
    if (!is.null(exp)) sprintf("\\left(%s\\right)^{%s}", tex, exp) else tex
  },
  DiracDelta = function(e, exp = NULL) {
    tex <- if (length(e$args) == 1L || identical(e$args[[2L]]$key, S_Zero$key)) sprintf("\\delta\\left(%s\\right)", .lp_arg(e))
           else sprintf("\\delta^{\\left( %s \\right)}\\left( %s \\right)", .lp_arg(e, 2L), .lp_arg(e, 1L))
    if (!is.null(exp) && nzchar(exp)) tex <- sprintf("\\left(%s\\right)^{%s}", tex, exp)
    tex
  },
  SingularityFunction = function(e, exp = NULL) {
    shift <- .latex(.sub(e$args[[1L]], e$args[[2L]])); power <- .lp_arg(e, 3L)
    if (!is.null(exp)) sprintf("{\\left({\\langle %s \\rangle}^{%s}\\right)}^{%s}", shift, power, exp)
    else sprintf("{\\left\\langle %s \\right\\rangle}^{%s}", shift, power)
  },
  Heaviside = function(e, exp = NULL) {
    # pargs: the arguments without a default H(0) of 1/2
    pa <- e$args
    if (length(pa) == 2L && identical(pa[[2L]]$key, S_Half$key)) pa <- pa[1L]
    tex <- sprintf("\\theta\\left(%s\\right)", paste(vapply(pa, .latex, ""), collapse = ", "))
    if (!is.null(exp) && nzchar(exp)) tex <- sprintf("\\left(%s\\right)^{%s}", tex, exp)
    tex
  },
  KroneckerDelta = function(e, exp = NULL) {
    i <- .lp_arg(e, 1L); j <- .lp_arg(e, 2L)
    tex <- if (.lp_is_atom(e$args[[1L]]) && .lp_is_atom(e$args[[2L]])) sprintf("\\delta_{%s %s}", i, j)
           else sprintf("\\delta_{%s, %s}", i, j)
    if (!is.null(exp)) tex <- sprintf("\\left(%s\\right)^{%s}", tex, exp)
    tex
  },
  LeviCivita = function(e, exp = NULL) {
    ind <- vapply(e$args, .latex, "")
    tex <- if (all(vapply(e$args, .lp_is_atom, NA))) sprintf("\\varepsilon_{%s}", paste(ind, collapse = " "))
           else sprintf("\\varepsilon_{%s}", paste(ind, collapse = ", "))
    if (!is.null(exp) && nzchar(exp)) tex <- sprintf("\\left(%s\\right)^{%s}", tex, exp)
    tex
  },
  bernoulli = .lp_number_poly("B"), genocchi = .lp_number_poly("G"), fibonacci = .lp_number_poly("F"),
  tribonacci = .lp_number_poly("T"),
  bell = function(e, exp = NULL) {
    if (length(e$args) == 3L) {
      t1 <- sprintf("B_{%s, %s}", .lp_arg(e, 1L), .lp_arg(e, 2L))
      t2 <- sprintf("\\left(%s\\right)", paste(vapply(.lp_tuple_args(e$args[[3L]]), .latex, ""), collapse = ", "))
      return(if (!is.null(exp)) sprintf("%s^{%s}%s", t1, exp, t2) else paste0(t1, t2))
    }
    .lp_number_poly("B")(e, exp)
  },
  lucas = function(e, exp = NULL) .lp_sup(sprintf("L_{%s}", .lp_arg(e)), exp),
  mobius = function(e, exp = NULL) {
    if (is.null(exp)) sprintf("\\mu\\left(%s\\right)", .lp_arg(e)) else sprintf("\\mu^{%s}\\left(%s\\right)", exp, .lp_arg(e))
  },
  euler = function(e, exp = NULL) {
    tex <- .lp_sup(sprintf("E_{%s}", .lp_arg(e, 1L)), exp)
    if (length(e$args) == 2L) tex <- sprintf("%s\\left(%s\\right)", tex, .lp_arg(e, 2L))
    tex
  },
  catalan = function(e, exp = NULL) .lp_sup(sprintf("C_{%s}", .lp_arg(e)), exp),
  LambertW = function(e, exp = NULL) {
    ex <- if (!is.null(exp)) sprintf("^{%s}", exp) else ""
    if (length(e$args) == 1L) sprintf("W%s\\left(%s\\right)", ex, .lp_arg(e))
    else sprintf("W%s_{%s}\\left(%s\\right)", ex, .lp_arg(e, 2L), .lp_arg(e, 1L))
  },
  totient = .lp_ntheory1("\\phi"), reduced_totient = .lp_ntheory1("\\lambda"),
  primenu = .lp_ntheory1("\\nu"), primeomega = .lp_ntheory1("\\Omega"),
  divisor_sigma = .lp_divisor_sigma(FALSE), udivisor_sigma = .lp_divisor_sigma(TRUE)
)

.lp_variadic <- function(e, name, exp = NULL) {
  tex <- sprintf("\\%s\\left(%s\\right)", name, paste(vapply(.lp_dsk(e$args), .latex, ""), collapse = ", "))
  .lp_sup(tex, exp)
}

.lp_betainc <- function(e, exp, op) {
  a <- vapply(e$args, .latex, "")
  tex <- sprintf("\\left(%s, %s\\right)", a[1L], a[2L])
  if (!is.null(exp)) sprintf("\\operatorname{%s}_{(%s, %s)}^{%s}%s", op, a[3L], a[4L], exp, tex)
  else sprintf("\\operatorname{%s}_{(%s, %s)}%s", op, a[3L], a[4L], tex)
}

# ---------------------------------------------------------------------------
# Calculus

.lp_requires_partial <- function(e) {
  if (e$type == "Derivative") return(.lp_requires_partial(e$args[[1L]]))
  fs <- .lp_free(e)
  sum(vapply(fs, function(s) !isTRUE(.facts(s)[["integer"]]), NA)) > 1L
}

.lp_has_type <- function(e, ty) is_sym(e) && (e$type == ty || any(vapply(e$args, .lp_has_type, NA, ty = ty)))

# Free symbols as SymPy counts them: a Derivative contributes its expression
# and its derivative counts, not its variables.
.lp_free <- function(e) {
  if (!is_sym(e)) return(list())
  if (e$type == "Symbol") return(list(e))
  if (e$type == "Derivative") {
    l <- c(.lp_free(e$args[[1L]]), do.call(c, lapply(e$args[-1L], function(vc) .lp_free(vc$args[[2L]]))))
  } else if (!length(e$args) || !.lp_has_type(e, "Derivative")) return(free_symbols(e))
  else l <- do.call(c, lapply(e$args, .lp_free))
  if (!length(l)) return(list())
  l[!duplicated(vapply(l, function(s) s$key, ""))]
}

.lp_Derivative <- function(e, exp = NULL) {
  ex <- e$args[[1L]]
  diff_symbol <- if (.lp_requires_partial(ex)) "\\partial" else .lp_s("diff_operator_latex")
  tex <- ""
  dim <- 0
  for (vc in rev(e$args[-1L])) {
    v <- vc$args[[1L]]; num <- vc$args[[2L]]
    if (num$type == "Integer") dim <- dim + as.numeric(num$value) else dim <- .add(list(as_sym(dim), num))
    if (is_one_num(num)) tex <- paste0(tex, sprintf("%s %s", diff_symbol, .latex(v)))
    else tex <- paste0(tex, sprintf("%s %s^{%s}", diff_symbol, .lp_paren_super(.latex(v)), .latex(num)))
  }
  if (is.numeric(dim) && dim == 1) tex <- sprintf("\\frac{%s}{%s}", diff_symbol, tex)
  else tex <- sprintf("\\frac{%s^{%s}}{%s}", diff_symbol, .latex(as_sym(dim)), tex)
  prec <- 50
  if (!is.null(.lp_s("mul_symbol"))) prec <- prec + 1
  neg <- any(vapply(e$args, .lp_cems, NA))
  sprintf("%s %s", tex, .lp_paren(ex, prec, is_neg = neg, strict = TRUE))
}

.lp_Subs <- function(e, exp = NULL) {
  old <- .lp_tuple_args(e$args[[2L]]); new <- .lp_tuple_args(e$args[[3L]])
  subs <- paste(vapply(seq_along(old), function(i) paste0(.latex(old[[i]]), "=", .latex(new[[i]])), ""), collapse = "\\\\ ")
  sprintf("\\left. %s \\right|_{\\substack{ %s }}", .latex(e$args[[1L]]), subs)
}

.lp_Integral <- function(e, exp = NULL) {
  lims <- e$args[-1L]
  ds <- .lp_s("diff_operator_latex")
  tex <- ""
  symbols <- character()
  if (length(lims) <= 4L && all(vapply(lims, function(l) length(l$args) == 1L, NA))) {
    tex <- paste0("\\i", strrep("i", length(lims) - 1L), "nt")
    symbols <- vapply(lims, function(l) sprintf("\\, %s%s", ds, .latex(l$args[[1L]])), "")
  } else {
    for (lim in rev(lims)) {
      tex <- paste0(tex, "\\int")
      if (length(lim$args) > 1L) {
        if (.lp_s("mode") != "inline" && !isTRUE(.lp_s("itex"))) tex <- paste0(tex, "\\limits")
        if (length(lim$args) == 3L) tex <- paste0(tex, sprintf("_{%s}^{%s}", .latex(lim$args[[2L]]), .latex(lim$args[[3L]])))
        if (length(lim$args) == 2L) tex <- paste0(tex, sprintf("^{%s}", .latex(lim$args[[2L]])))
      }
      symbols <- c(sprintf("\\, %s%s", ds, .latex(lim$args[[1L]])), symbols)
    }
  }
  neg <- any(vapply(e$args, .lp_cems, NA))
  sprintf("%s %s%s", tex, .lp_paren(e$args[[1L]], 50, is_neg = neg, strict = TRUE), paste(symbols, collapse = ""))
}

.lp_SumProd <- function(op) function(e, exp = NULL) {
  lims <- e$args[-1L]
  if (length(lims) == 1L) {
    l <- lims[[1L]]$args
    tex <- sprintf("\\%s_{%s=%s}^{%s} ", op, .latex(l[[1L]]), .latex(l[[2L]]), .latex(l[[3L]]))
  } else {
    ineq <- vapply(lims, function(l) sprintf("%s \\leq %s \\leq %s", .latex(l$args[[2L]]), .latex(l$args[[1L]]), .latex(l$args[[3L]])), "")
    tex <- sprintf("\\%s_{\\substack{%s}} ", op, paste(ineq, collapse = "\\\\"))
  }
  f <- e$args[[1L]]
  paste0(tex, if (f$type == "Add") .lp_parens(.latex(f)) else .latex(f))
}

.lp_Limit <- function(e, exp = NULL) {
  ex <- e$args[[1L]]; z <- e$args[[2L]]; z0 <- e$args[[3L]]
  dir <- e$args[[4L]]$name
  tex <- sprintf("\\lim_{%s \\to ", .latex(z))
  if (dir == "+-" || is_const(z0, "oo") || is_const(z0, "-oo")) tex <- paste0(tex, .latex(z0), "}")
  else tex <- paste0(tex, .latex(z0), "^", dir, "}")
  if (ex$type %in% c("Add", "Mul", "And", "Or", "Xor") || (ex$type == "Function" && ex$name %in% c("Max", "Min")))
    sprintf("%s\\left(%s\\right)", tex, .latex(ex))
  else sprintf("%s %s", tex, .latex(ex))
}

.lp_Order <- function(e, exp = NULL) {
  s <- .latex(e$args[[1L]])
  pts <- lapply(e$args[-1L], function(t) t$args[[2L]])
  vars <- lapply(e$args[-1L], function(t) t$args[[1L]])
  if ((length(pts) && any(vapply(pts, function(p) !identical(p$key, S_Zero$key), NA))) || length(vars) > 1L) {
    s <- paste0(s, "; ")
    if (length(vars) > 1L) s <- paste0(s, .latex(do.call(Tuple, vars)))
    else if (length(vars)) s <- paste0(s, .latex(vars[[1L]]))
    s <- paste0(s, "\\rightarrow ")
    if (length(pts) > 1L) s <- paste0(s, .latex(do.call(Tuple, pts))) else s <- paste0(s, .latex(pts[[1L]]))
  }
  sprintf("O\\left(%s\\right)", s)
}

.lp_Lambda <- function(e, exp = NULL) {
  v <- e$args[[1L]]$args
  vs <- if (length(v) == 1L) .latex(v[[1L]]) else .lp_tuple(v)
  sprintf("\\left( %s \\mapsto %s \\right)", vs, .latex(e$args[[2L]]))
}

# ---------------------------------------------------------------------------
# Relations and logic

.lp_Relational <- function(e, exp = NULL) {
  op <- switch(e$type, Equality = "=", Unequality = "\\neq",
               StrictLessThan = if (isTRUE(.lp_s("itex"))) "\\lt" else "<", LessThan = "\\leq",
               StrictGreaterThan = if (isTRUE(.lp_s("itex"))) "\\gt" else ">", GreaterThan = "\\geq")
  sprintf("%s %s %s", .latex(e$args[[1L]]), op, .latex(e$args[[2L]]))
}

.lp_is_boolean <- function(e) is_sym(e) && e$type %in% .LP_BOOL
.lp_logop <- function(args, ch) {
  pr <- function(a) if (.lp_is_boolean(a) && a$type != "Not") .lp_parens(.latex(a)) else .latex(a)
  tex <- pr(args[[1L]])
  for (a in args[-1L]) tex <- paste0(tex, " ", ch, " ", pr(a))
  tex
}

.lp_bool_sorted <- function(e) .lp_dsk(e$args)

.lp_Not <- function(e, exp = NULL) {
  a <- e$args[[1L]]
  if (a$type == "Equivalent") return(.lp_logop(.lp_bool_sorted(a), "\\not\\Leftrightarrow"))
  if (a$type == "Implies") return(.lp_logop(a$args, "\\not\\Rightarrow"))
  if (.lp_is_boolean(a)) sprintf("\\neg \\left(%s\\right)", .latex(a)) else sprintf("\\neg %s", .latex(a))
}

.lp_Piecewise <- function(e, exp = NULL) {
  pr <- .pw_pairs(e)
  n <- length(pr)
  ec <- vapply(seq_len(n), function(i) {
    p <- pr[[i]]
    if (i == n && .is_true(p[[2L]])) sprintf("%s & \\text{otherwise}", .latex(p[[1L]]))
    else sprintf("%s & \\text{for}\\: %s", .latex(p[[1L]]), .latex(p[[2L]]))
  }, "")
  sprintf("\\begin{cases} %s \\end{cases}", paste(ec, collapse = " \\\\"))
}

# ---------------------------------------------------------------------------
# Containers

.lp_tuple <- function(items) {
  sep <- if (.lp_s("decimal_separator") == "comma") ";" else ","
  if (length(items) == 1L) return(paste0("\\left( ", .latex(items[[1L]]), sep, "\\right)"))
  paste0("\\left( ", paste(vapply(items, .latex, ""), collapse = paste0(sep, " \\  ")), "\\right)")
}

.lp_Dict <- function(e, exp = NULL) {
  ps <- e$args
  keys <- lapply(ps, function(p) p$args[[1L]])
  ord <- .key_sort(lapply(keys, .sort_key), seq_along(keys))
  items <- vapply(ord, function(i) sprintf("%s : %s", .latex(ps[[i]]$args[[1L]]), .latex(ps[[i]]$args[[2L]])), "")
  sprintf("\\left\\{ %s\\right\\}", paste(items, collapse = ", \\  "))
}

.lp_matrix_contents <- function(rows, ncol) {
  lines <- vapply(rows, function(r) paste(vapply(r, .latex, ""), collapse = " & "), "")
  ms <- .lp_s("mat_str")
  if (is.null(ms)) ms <- if (.lp_s("mode") == "inline") "smallmatrix" else if (ncol <= 10) "matrix" else "array"
  body <- paste(lines, collapse = "\\\\")
  if (ms == "array") body <- paste0("{", strrep("c", ncol), "}", body)
  sprintf("\\begin{%s}%s\\end{%s}", ms, body, ms)
}

.lp_matrix <- function(M) {
  out <- .lp_matrix_contents(.mrows(M), M$ncol)
  d <- .lp_s("mat_delim")
  if (!is.null(d) && nzchar(d)) out <- paste0("\\left", d, out, "\\right", c(`(` = ")", `[` = "]")[[d]])
  out
}

.lp_Poly <- function(p) {
  ms <- monoms(p); cs <- coeffs(p)
  ng <- .poly_ngens(p)
  gens <- p$gens[seq_len(ng)]
  terms <- character()
  for (i in seq_along(ms)) {
    mon <- ms[[i]]; coeff <- cs[[i]]
    s_monom <- ""
    for (j in seq_along(mon)) {
      if (mon[j] > 0) s_monom <- paste0(s_monom, if (mon[j] == 1) .latex(gens[[j]]) else .latex(.pow(gens[[j]], .int(mon[j]))))
    }
    if (coeff$type == "Add") {
      s_coeff <- if (nzchar(s_monom)) .lp_parens(.latex(coeff)) else .latex(coeff)
    } else {
      if (nzchar(s_monom)) {
        if (is_one_num(coeff)) { terms <- c(terms, "+", s_monom); next }
        if (identical(coeff$key, S_NegativeOne$key)) { terms <- c(terms, "-", s_monom); next }
      }
      s_coeff <- .latex(coeff)
    }
    s_term <- if (!nzchar(s_monom)) s_coeff else paste(s_coeff, s_monom)
    terms <- if (startsWith(s_term, "-")) c(terms, "-", substring(s_term, 2L)) else c(terms, "+", s_term)
  }
  if (!length(terms)) terms <- c("+", "0")
  if (terms[1L] %in% c("-", "+")) {
    mod <- terms[1L]; terms <- terms[-1L]
    if (mod == "-") terms[1L] <- paste0("-", terms[1L])
  }
  ex <- paste(terms, collapse = " ")
  dom <- if (!is.null(p$dom_str)) p$dom_str else if (!is.null(p$rr_prec)) "RR" else if (ng < length(p$gens))
    paste0(.domain_of(p$rep), "[", paste(vapply(p$gens[-seq_len(ng)], .str, ""), collapse = ","), "]") else .domain_of(p$rep)
  args <- paste(c(ex, vapply(gens, .latex, ""), paste0("domain=", .lp_domain(dom))), collapse = ", ")
  sprintf("\\operatorname{Poly}{\\left( %s \\right)}", args)
}

.lp_domain <- function(d) {
  base <- c(ZZ = "\\mathbb{Z}", QQ = "\\mathbb{Q}", RR = "\\mathbb{R}", CC = "\\mathbb{C}")
  if (d %in% names(base)) return(base[[d]])
  m <- regmatches(d, regexec("^(ZZ|QQ|RR|CC)([[(])(.*)[])]$", d))[[1L]]
  if (length(m)) {
    syms <- paste(vapply(strsplit(m[4L], ",")[[1L]], function(s) .latex(sympify(trimws(s))), ""), collapse = ", ")
    return(if (m[3L] == "[") sprintf("%s\\left[%s\\right]", base[[m[2L]]], syms) else sprintf("%s\\left(%s\\right)", base[[m[2L]]], syms))
  }
  sprintf("\\mathtt{\\text{%s}}", latex_escape(d))
}

# ---------------------------------------------------------------------------
# Sets

.lp_setop <- function(sym) function(e, exp = NULL) {
  args <- if (e$type %in% c("Union", "Intersection")) .set_print_order(e$args) else e$args
  pr <- .lp_prec(e)
  paste(vapply(args, function(a) .lp_paren(a, pr), ""), collapse = sym)
}

.lp_Interval <- function(e, exp = NULL) {
  i <- .iv(e)
  if (identical(i$a$key, i$b$key)) return(sprintf("\\left\\{%s\\right\\}", .latex(i$a)))
  l <- if (i$lo) .lp_s("open_delim")[1L] else .lp_s("closed_delim")[1L]
  r <- if (i$ro) .lp_s("open_delim")[2L] else .lp_s("closed_delim")[2L]
  sprintf("\\left%s%s, %s\\right%s", l, .latex(i$a), .latex(i$b), r)
}

.lp_FiniteSet <- function(e, exp = NULL) {
  items <- .lp_dsk(e$args)
  sep <- if (.lp_s("decimal_separator") == "comma") "; " else ", "
  sprintf("\\left\\{%s\\right\\}", paste(vapply(items, .latex, ""), collapse = sep))
}

.lp_ProductSet <- function(e, exp = NULL) {
  pr <- .lp_prec(e)
  sets <- e$args
  if (length(sets) >= 1L && length(unique(vapply(sets, `[[`, "", "key"))) == 1L)
    return(paste0(.lp_paren(sets[[1L]], pr), sprintf("^{%d}", length(sets))))
  paste(vapply(sets, function(s) .lp_paren(s, pr), ""), collapse = " \\times ")
}

.lp_SpecialSet <- function(e, exp = NULL) {
  m <- c(Reals = "\\mathbb{R}", Complexes = "\\mathbb{C}", Integers = "\\mathbb{Z}", Naturals = "\\mathbb{N}",
         Naturals0 = "\\mathbb{N}_0", Rationals = "\\mathbb{Q}", UniversalSet = "\\mathbb{U}")
  if (is.na(m[e$name])) .lp_Basic(e) else m[[e$name]]
}

.lp_ImageSet <- function(e, exp = NULL) {
  lam <- e$args[[1L]]
  sig <- lam$args[[1L]]$args
  bases <- e$args[-1L]
  xin <- paste(vapply(seq_along(sig), function(i) sprintf("%s \\in %s", .latex(sig[[i]]), .latex(bases[[i]])), ""), collapse = ", ")
  sprintf("\\left\\{%s\\; \\middle|\\; %s\\right\\}", .latex(lam$args[[2L]]), xin)
}

.lp_ConditionSet <- function(e, exp = NULL) {
  v <- e$args[[1L]]
  vars <- .latex(v)
  base <- e$args[[3L]]
  if (base$type == "SpecialSet" && base$name == "UniversalSet")
    return(sprintf("\\left\\{%s\\; \\middle|\\; %s \\right\\}", vars, .latex(e$args[[2L]])))
  sprintf("\\left\\{%s\\; \\middle|\\; %s \\in %s \\wedge %s \\right\\}", vars, vars, .latex(base), .latex(e$args[[2L]]))
}

# ---------------------------------------------------------------------------
# Matrix expressions

.lp_is_matsym <- function(m) m$type %in% c("MatrixSymbol", "BlockMatrix")
.lp_is_matexpr <- function(m) !(m$type %in% c("ImmutableMatrix", "ImmutableDenseMatrix")) && isTRUE(tryCatch(.is_matexpr(m), error = function(err) FALSE))

.lp_transpose_like <- function(sym) function(e, exp = NULL) {
  m <- e$args[[1L]]
  if (!.lp_is_matsym(m) && .lp_is_matexpr(m)) return(sprintf("\\left(%s\\right)^{%s}", .latex(m), sym))
  s <- .lp_paren(m, .lp_prec(e), strict = TRUE)
  if (grepl("^", s, fixed = TRUE)) sprintf("\\left(%s\\right)^{%s}", s, sym) else sprintf("%s^{%s}", s, sym)
}

.lp_MatMul <- function(e, exp = NULL) {
  pr <- .lp_prec(e)
  parens <- function(z) if (z$type == "Mul") .latex(z) else .lp_paren(z, pr, strict = FALSE)
  args <- e$args
  if (.lp_cems(e)) {
    if (identical(args[[1L]]$key, S_NegativeOne$key)) args <- args[-1L] else args[[1L]] <- .neg(args[[1L]])
    return(paste0("- ", paste(vapply(args, parens, ""), collapse = " ")))
  }
  paste(vapply(args, parens, ""), collapse = " ")
}

.lp_MatPow <- function(e, exp = NULL) {
  b <- e$args[[1L]]; x <- e$args[[2L]]
  if (b$type != "MatrixSymbol" && .lp_is_matexpr(b)) return(sprintf("\\left(%s\\right)^{%s}", .latex(b), .latex(x)))
  bs <- .latex(b)
  if (grepl("^", bs, fixed = TRUE)) sprintf("\\left(%s\\right)^{%s}", bs, .latex(x)) else sprintf("%s^{%s}", bs, .latex(x))
}

.lp_MatrixSymbol <- function(e, exp = NULL) .lp_Symbol(e, style = .lp_s("mat_symbol_style"))

.lp_MatrixElement <- function(e, exp = NULL)
  sprintf("{%s}_{%s,%s}", .lp_paren(e$args[[1L]], 1000, strict = TRUE), .latex(e$args[[2L]]), .latex(e$args[[3L]]))

# ---------------------------------------------------------------------------
# Fallback

.lp_Basic <- function(e, exp = NULL) {
  name <- .lp_deal_with_super_sub(e$type)
  if (length(e$args)) sprintf("\\operatorname{%s}\\left(%s\\right)", name, paste(vapply(e$args, .latex, ""), collapse = ", "))
  else sprintf("\\text{%s}", name)
}

.lp_printers <- list(
  Integer = .lp_Integer, Rational = .lp_Rational, Float = .lp_Float, Constant = .lp_Constant,
  BooleanAtom = .lp_BooleanAtom,
  Symbol = function(e, exp = NULL) .lp_Symbol(e), Dummy = function(e, exp = NULL) .lp_Symbol(e),
  Wild = function(e, exp = NULL) {
    if (isTRUE(e$extra$wildfunction)) paste0("\\operatorname{WildFunction}{\\left(", .lp_deal_with_super_sub(e$name), " \\right)}")
    else .lp_Symbol(e)
  },
  Add = .lp_Add, Mul = .lp_Mul, Pow = .lp_Pow,
  UnevaluatedExpr = function(e, exp = NULL) .latex(e$args[[1L]]),
  Derivative = .lp_Derivative, Integral = .lp_Integral, Subs = .lp_Subs,
  Sum = .lp_SumProd("sum"), Product = .lp_SumProd("prod"), Limit = .lp_Limit, Order = .lp_Order,
  Lambda = .lp_Lambda,
  Equality = .lp_Relational, Unequality = .lp_Relational, StrictLessThan = .lp_Relational,
  LessThan = .lp_Relational, StrictGreaterThan = .lp_Relational, GreaterThan = .lp_Relational,
  And = function(e, exp = NULL) .lp_logop(.lp_bool_sorted(e), "\\wedge"),
  Or = function(e, exp = NULL) .lp_logop(.lp_bool_sorted(e), "\\vee"),
  Xor = function(e, exp = NULL) .lp_logop(.lp_bool_sorted(e), "\\veebar"),
  Implies = function(e, exp = NULL) .lp_logop(e$args, "\\Rightarrow"),
  Equivalent = function(e, exp = NULL) .lp_logop(.lp_bool_sorted(e), "\\Leftrightarrow"),
  Not = .lp_Not, Piecewise = .lp_Piecewise,
  Tuple = function(e, exp = NULL) .lp_tuple(e$args),
  Dict = .lp_Dict,
  Interval = .lp_Interval, FiniteSet = .lp_FiniteSet, EmptySet = function(e, exp = NULL) "\\emptyset",
  SpecialSet = .lp_SpecialSet, UniversalSet = function(e, exp = NULL) "\\mathbb{U}",
  Union = .lp_setop(" \\cup "), Intersection = .lp_setop(" \\cap "), Complement = .lp_setop(" \\setminus "),
  SymmetricDifference = .lp_setop(" \\triangle "), ProductSet = .lp_ProductSet,
  ImageSet = .lp_ImageSet, ConditionSet = .lp_ConditionSet,
  PowerSet = function(e, exp = NULL) sprintf("\\mathcal{P}\\left(%s\\right)", .latex(e$args[[1L]])),
  Contains = function(e, exp = NULL) sprintf("%s \\in %s", .latex(e$args[[1L]]), .latex(e$args[[2L]])),
  CRootOf = function(e, exp = NULL) sprintf("\\operatorname{CRootOf} {\\left(%s, %s\\right)}", .latex(e$args[[1L]]), .lp_Integer(e$args[[2L]])),
  ComplexRootOf = function(e, exp = NULL) sprintf("\\operatorname{CRootOf} {\\left(%s, %s\\right)}", .latex(e$args[[1L]]), .lp_Integer(e$args[[2L]])),
  Indexed = function(e, exp = NULL) paste0("{", .latex(e$args[[1L]]), "}_{", paste(vapply(e$args[-1L], .latex, ""), collapse = ","), "}"),
  IndexedBase = function(e, exp = NULL) .latex(e$args[[1L]]),
  MatrixSymbol = .lp_MatrixSymbol, MatMul = .lp_MatMul, MatPow = .lp_MatPow, MatrixElement = .lp_MatrixElement,
  Transpose = .lp_transpose_like("T"),
  Adjoint = function(e, exp = NULL) .lp_transpose_like(switch(.lp_s("adjoint_style"), star = "\\ast", hermitian = "\\mathsf{H}", "\\dagger"))(e),
  Identity = function(e, exp = NULL) if (.lp_s("mat_symbol_style") == "plain") "\\mathbb{I}" else "\\mathbf{I}",
  ZeroMatrix = function(e, exp = NULL) if (.lp_s("mat_symbol_style") == "plain") "0" else "\\mathbf{0}",
  OneMatrix = function(e, exp = NULL) if (.lp_s("mat_symbol_style") == "plain") "1" else "\\mathbf{1}",
  HadamardProduct = function(e, exp = NULL) paste(vapply(e$args, function(a) .lp_paren(a, 60, strict = TRUE), ""), collapse = " \\circ "),
  KroneckerProduct = function(e, exp = NULL) paste(vapply(e$args, function(a) .lp_paren(a, 60, strict = TRUE), ""), collapse = " \\otimes "),
  HadamardPower = function(e, exp = NULL) {
    tpl <- if (.lp_prec(e$args[[2L]]) < 50) "%s^{\\circ \\left({%s}\\right)}" else "%s^{\\circ {%s}}"
    .lp_helper_std_power(e$args[[1L]], e$args[[2L]], tpl)
  },
  Trace = function(e, exp = NULL) sprintf("\\operatorname{tr}\\left(%s \\right)", .latex(e$args[[1L]])),
  Expectation = function(e, exp = NULL) sprintf("\\operatorname{E}\\left[%s\\right]", .latex(e$args[[1L]])),
  Variance = function(e, exp = NULL) sprintf("\\operatorname{Var}\\left(%s\\right)", .latex(e$args[[1L]])),
  Covariance = function(e, exp = NULL) sprintf("\\operatorname{Cov}\\left(%s\\right)", paste(vapply(e$args, .latex, ""), collapse = ", ")),
  Probability = function(e, exp = NULL) sprintf("\\operatorname{P}\\left(%s\\right)", .latex(e$args[[1L]])),
  LaplaceTransform = .lp_transform("L"), InverseLaplaceTransform = .lp_transform("L", TRUE),
  FourierTransform = .lp_transform("F"), InverseFourierTransform = .lp_transform("F", TRUE),
  MellinTransform = .lp_transform("M"), InverseMellinTransform = .lp_transform("M", TRUE),
  SineTransform = .lp_transform("SIN"), InverseSineTransform = .lp_transform("SIN", TRUE),
  CosineTransform = .lp_transform("COS"), InverseCosineTransform = .lp_transform("COS", TRUE),
  Str = function(e, exp = NULL) e$name
)

# ---------------------------------------------------------------------------
# Further printers of LatexPrinter

.lp_Range <- function(e, exp = NULL) {
  a <- e$args
  symbolic <- function() {
    cont <- if (identical(a[[1L]]$key, S_Zero$key)) {
      if (is_one_num(a[[3L]])) .latex(a[[2L]]) else paste(vapply(a, .latex, ""), collapse = ", ")
    } else if (is_one_num(a[[3L]])) paste(vapply(a[1:2], .latex, ""), collapse = ", ")
    else paste(vapply(a, .latex, ""), collapse = ", ")
    sprintf("\\text{Range}\\left(%s\\right)", cont)
  }
  p <- .range_parts(e)
  dots <- "\\ldots"
  if (is_inf_const(p$start) && is_inf_const(p$stop)) {
    ps <- if (num_sign(p$step) > 0) c(dots, "-1", "0", "1", dots) else c(dots, "1", "0", "-1", dots)
  } else if (is_inf_const(p$start)) {
    last <- .range_last(e)
    ps <- c(dots, .latex(.sub(last, p$step)), .latex(last))
  } else if (is_inf_const(p$stop)) {
    ps <- c(.latex(p$start), .latex(.add(list(p$start, p$step))), dots)
  } else {
    syms <- do.call(c, lapply(a, free_symbols))
    if (length(syms) && !all(vapply(syms, function(s) isTRUE(.facts(s)[["integer"]]), NA))) return(symbolic())
    n <- .range_size(e)
    if (!is_exact_number(n)) return(symbolic())
    k <- .int_val(n)
    if (k < 4) ps <- vapply(seq_len(k) - 1L, function(i) .latex(.range_item(e, i)), "")
    else ps <- c(.latex(.range_item(e, 0L)), .latex(.range_item(e, 1L)), dots, .latex(.range_last(e)))
    if (k <= 0) ps <- character()
  }
  paste0("\\left\\{", paste(ps, collapse = ", "), "\\right\\}")
}

# Permutation and Cycle (S3 objects)
.lp_Cycle_perm <- function(p) {
  a <- p$array
  if (!length(a)) return("\\left( \\right)")
  cyc <- cyclic_form(p)
  if (a[length(a)] == length(a) - 1L) cyc <- c(cyc, list(length(a) - 1L))
  paste0(vapply(cyc, function(cy) paste0("\\left( ", paste(cy, collapse = "\\; "), "\\right)"), ""), collapse = "")
}

.lp_Permutation <- function(p) {
  if (!identical(.lp_s("perm_cyclic"), FALSE)) return(.lp_Cycle_perm(p))
  a <- p$array
  if (!length(a)) return("\\left( \\right)")
  sprintf("\\begin{pmatrix} %s \\\\ %s \\end{pmatrix}", paste(seq_along(a) - 1L, collapse = " & "), paste(a, collapse = " & "))
}

.lp_Cycle <- function(cy) {
  if (!length(cy$keys)) return("\\left( \\right)")
  .lp_Cycle_perm(Permutation(cy))
}

.lp_omega_power <- function(op) {
  ex <- op$exp; mul <- as.character(op$mult)
  ex0 <- inherits(ex, "sym_ordinal") && !length(ex$terms)
  ex1 <- inherits(ex, "sym_ordinal") && length(ex$terms) == 1L && .lp_ord_is_one(ex)
  if (ex0) return(mul)
  ext <- .lp_ordinal(ex)
  if (mul != "1") { if (!ex1) sprintf("\\omega^{%s} %s", ext, mul) else sprintf("\\omega %s", mul) }
  else if (!ex1) sprintf("\\omega^{%s}", ext) else "\\omega"
}
.lp_ord_is_one <- function(o) {
  t <- o$terms[[1L]]
  inherits(t$exp, "sym_ordinal") && !length(t$exp$terms) && as.character(t$mult) == "1"
}
.lp_ordinal <- function(o) {
  if (!length(o$terms)) return("0")
  paste(vapply(o$terms, .lp_omega_power, ""), collapse = " + ")
}

.lp_FreeGroupElement <- function(w) {
  if (.fw_is_id(w)) return("1")
  nm <- w$group$symbols[w$s]
  parts <- ifelse(w$e == 1L, nm, paste0(nm, "^{", w$e, "}"))
  paste(parts, collapse = .lp_s("mul_symbol_latex"))
}

.lp_Quaternion <- function(q) {
  s <- vapply(q$c, function(i) .lp_paren(i, 50, strict = TRUE), "")
  paste(c(s[1L], paste(s[-1L], c("i", "j", "k"))), collapse = " + ")
}

# FunctionClass objects: the special function classes, else _hprint_Function
.LP_SPECIAL_CLASSES <- c(KroneckerDelta = "\\delta", gamma = "\\Gamma", lowergamma = "\\gamma",
                         beta = "\\operatorname{B}", DiracDelta = "\\delta", Chi = "\\operatorname{Chi}")

.lp_function_object <- function(f) {
  nm <- if (is.primitive(f)) {
    m <- regmatches(deparse(f), regexec("\\.Primitive\\(\"([^\"]+)\"\\)", deparse(f)))[[1L]]
    if (length(m)) m[2L] else NULL
  } else if (!is.null(environment(f)$name) && is.character(environment(f)$name)) environment(f)$name
  if (is.null(nm)) {
    ns <- asNamespace("symcompR")
    for (cand in getNamespaceExports(ns)) {
      g <- get0(cand, envir = ns, inherits = FALSE)
      if (is.function(g) && identical(g, f)) { nm <- cand; break }
    }
  }
  if (is.null(nm)) return(paste0("\\mathtt{\\text{", latex_escape(paste(deparse(f), collapse = " ")), "}}"))
  nm <- sub("^sym_", "", nm)
  if (nm %in% names(.LP_SPECIAL_CLASSES)) return(.LP_SPECIAL_CLASSES[[nm]])
  .lp_hprint_Function(nm)
}

.lp_immutable_rows <- function(e) {
  r <- .int_val(e$args[[1L]]); cc <- .int_val(e$args[[2L]])
  ent <- e$args[[3L]]$args
  list(rows = lapply(seq_len(r), function(i) ent[(i - 1L) * cc + seq_len(cc)]), ncol = cc)
}

.lp_ImmutableMatrix <- function(e, exp = NULL) {
  m <- .lp_immutable_rows(e)
  out <- .lp_matrix_contents(m$rows, m$ncol)
  d <- .lp_s("mat_delim")
  if (!is.null(d) && nzchar(d)) out <- paste0("\\left", d, out, "\\right", c(`(` = ")", `[` = "]")[[d]])
  out
}

.lp_block_rows <- function(e) {
  rows <- lapply(e$args[[1L]]$args, function(r) r$args)
  list(rows = rows, ncol = length(rows[[1L]]))
}

.lp_BlockMatrix <- function(e, exp = NULL) {
  m <- .lp_block_rows(e)
  out <- .lp_matrix_contents(m$rows, m$ncol)
  d <- .lp_s("mat_delim")
  if (!is.null(d) && nzchar(d)) out <- paste0("\\left", d, out, "\\right", c(`(` = ")", `[` = "]")[[d]])
  out
}

.lp_Determinant <- function(e, exp = NULL) {
  m <- e$args[[1L]]
  if (m$type == "BlockMatrix") { b <- .lp_block_rows(m); return(sprintf("\\left|{%s}\\right|", .lp_matrix_contents(b$rows, b$ncol))) }
  if (m$type == "ImmutableMatrix") { b <- .lp_immutable_rows(m); return(sprintf("\\left|{%s}\\right|", .lp_matrix_contents(b$rows, b$ncol))) }
  sprintf("\\left|{%s}\\right|", .latex(m))
}

# TransferFunction.to_expr: Mul(num, Pow(den, -1, evaluate=False), evaluate=False)
.lp_tf_expr <- function(tf) {
  num <- tf$args[[1L]]; den <- tf$args[[2L]]
  inv <- new_sym("Pow", list(den, S_NegativeOne))
  # AssocOp drops the identity even when evaluate=False
  if (identical(num$key, S_One$key)) return(inv)
  new_sym("Mul", list(num, inv))
}

.lp_mimo_matrix <- function(e) {
  rows <- lapply(e$args[[1L]]$args, function(r) lapply(r$args, .lp_tf_expr))
  out <- .lp_matrix_contents(rows, length(rows[[1L]]))
  d <- .lp_s("mat_delim")
  if (!is.null(d) && nzchar(d)) out <- paste0("\\left", d, out, "\\right", c(`(` = ")", `[` = "]")[[d]])
  out
}

.lp_tf_frac <- function(e, exp = NULL) sprintf("\\frac{%s}{%s}", .latex(e$args[[1L]]), .latex(e$args[[2L]]))

.lp_Feedback <- function(e, exp = NULL) {
  sys1 <- e$args[[1L]]; sys2 <- e$args[[2L]]
  sgn <- if (length(e$args) >= 3L) e$args[[3L]] else S_NegativeOne
  var <- if (sys1$type == "Series") sys1$args[[1L]]$args[[3L]] else sys1$args[[3L]]
  tf <- TransferFunction(1L, 1L, var)
  num_l <- if (sys1$type == "Series") sys1$args else list(sys1)
  den_l <- if (sys2$type == "Series") sys2$args else list(sys2)
  same <- function(a, b) identical(a$key, b$key)
  ser <- function(l) do.call(Series, l)
  d2 <- if (sys1$type == "Series" && sys2$type == "Series") ser(c(num_l, den_l))
        else if (sys1$type == "Series" && sys2$type == "TransferFunction") { if (same(sys2, tf)) ser(num_l) else ser(c(num_l, list(sys2))) }
        else if (sys1$type == "TransferFunction" && sys2$type == "Series") { if (same(sys1, tf)) ser(den_l) else ser(c(list(sys1), den_l)) }
        else if (same(sys1, tf)) ser(den_l) else if (same(sys2, tf)) ser(num_l) else ser(c(num_l, den_l))
  s <- if (identical(sgn$key, S_NegativeOne$key)) "+" else "-"
  sprintf("\\frac{%s}{%s %s %s}", .latex(sys1), .latex(tf), s, .latex(d2))
}

.lp_printers2 <- list(
  Range = .lp_Range,
  ImmutableMatrix = .lp_ImmutableMatrix, ImmutableDenseMatrix = .lp_ImmutableMatrix, BlockMatrix = .lp_BlockMatrix,
  Determinant = .lp_Determinant,
  MatAdd = .lp_Add,
  DotProduct = function(e, exp = NULL) {
    lv <- .lp_prec(e)
    sprintf("%s \\cdot %s", .lp_paren(e$args[[1L]], lv), .lp_paren(e$args[[2L]], lv))
  },
  ArraySymbol = function(e, exp = NULL) .latex(e$args[[1L]]),
  ArrayElement = function(e, exp = NULL)
    sprintf("{{%s}_{%s}}", .lp_paren(e$args[[1L]], 70, strict = TRUE), paste(vapply(.lp_tuple_args(e$args[[2L]]), .latex, ""), collapse = ", ")),
  PermutationMatrix = function(e, exp = NULL) {
    a <- vapply(.lp_tuple_args(e$args[[2L]]), function(z) .int_val(z), 0)
    sprintf("P_{%s}", .lp_Permutation(Permutation(as.integer(a))))
  },
  AppliedPermutation = function(e, exp = NULL) sprintf("\\sigma_{%s}(%s)", .latex(e$args[[1L]]), .latex(e$args[[2L]])),
  TransferFunction = .lp_tf_frac,
  DiscreteTransferFunction = function(e, exp = NULL)
    sprintf("\\frac{%s}{%s} \\text{ [st: } {%s} \\text{]}", .latex(e$args[[1L]]), .latex(e$args[[2L]]), .latex(e$args[[4L]])),
  Series = function(e, exp = NULL) { pr <- .lp_prec(e); paste(vapply(e$args, function(a) .lp_paren(a, pr), ""), collapse = " ") },
  Parallel = function(e, exp = NULL) paste(vapply(e$args, .latex, ""), collapse = " + "),
  MIMOParallel = function(e, exp = NULL) paste(vapply(e$args, .latex, ""), collapse = " + "),
  MIMOSeries = function(e, exp = NULL) {
    pr <- .lp_prec(e)
    paste(vapply(rev(e$args), function(a) if (a$type == "MIMOParallel") .lp_paren(a, pr) else .latex(a), ""), collapse = "\\cdot")
  },
  Feedback = .lp_Feedback,
  MIMOFeedback = function(e, exp = NULL) {
    sgn <- if (length(e$args) >= 3L) e$args[[3L]] else S_NegativeOne
    inv <- .latex(MIMOSeries(e$args[[2L]], e$args[[1L]]))
    sprintf("\\left(I_{\\tau} %s %s\\right)^{-1} \\cdot %s", if (identical(sgn$key, S_NegativeOne$key)) "+" else "-", inv, .latex(e$args[[1L]]))
  },
  TransferFunctionMatrix = function(e, exp = NULL) {
    mat <- .lp_mimo_matrix(e)
    el <- e$args[[1L]]$args[[1L]]$args[[1L]]
    if (el$type == "DiscreteTransferFunction") sprintf("\\underset{[st:\\ {%s}]}{%s_k}", sstr(el$args[[4L]]), mat)
    else sprintf("%s_\\tau", mat)
  },
  Manifold = function(e, exp = NULL) {
    p <- .lp_split_super_sub(.lp_strname(e$args[[1L]]))
    name <- sprintf("\\text{%s}", p$name)
    if (length(p$supers)) name <- paste0(name, "^{", paste(p$supers, collapse = " "), "}")
    if (length(p$subs)) name <- paste0(name, "_{", paste(p$subs, collapse = " "), "}")
    name
  },
  Patch = function(e, exp = NULL) sprintf("\\text{%s}_{%s}", .lp_strname(e$args[[1L]]), .latex(e$args[[2L]])),
  CoordSystem = function(e, exp = NULL)
    sprintf("\\text{%s}^{\\text{%s}}_{%s}", .lp_strname(e$args[[1L]]), .lp_strname(e$args[[2L]]$args[[1L]]), .latex(e$args[[2L]]$args[[2L]])),
  OrdinalOmega = function(e, exp = NULL) "\\omega",
  OrdinalZero = function(e, exp = NULL) "0"
)
for (.n in names(.lp_printers2)) .lp_printers[[.n]] <- .lp_printers2[[.n]]
rm(.n)

.lp_strname <- function(a) if (is_sym(a) && !is.null(a$name)) a$name else if (is.character(a)) a else .str(a)

# ReprPrinter._print_MatrixBase for a mutable matrix
.srepr_matrix <- function(m) {
  if (xor(m$nrow == 0L, m$ncol == 0L)) return(sprintf("%s(%d, %d, [])", .mclass(m), m$nrow, m$ncol))
  rows <- vapply(.mrows(m), function(r) paste0("[", paste(vapply(r, srepr, ""), collapse = ", "), "]"), "")
  paste0(.mclass(m), "([", paste(rows, collapse = ", "), "])")
}

# GLSLPrinter._print_MatrixBase with the default settings (mat_nested
# False, mat_separator ",\n", mat_transpose False, array_type "float",
# glsl_types True), with MatrixBase.table() padding each column to its
# widest entry
.glsl_matrix <- function(m, consts) {
  A <- if (m$ncol == 1L) t(m) else m
  rows <- .mrows(A)
  r <- A$nrow; cc <- A$ncol
  if (cc == 1L) return(.glsl(rows[[1L]][[1L]], consts))
  cells <- lapply(rows, function(row) vapply(row, .glsl, "", consts = consts))
  tab <- function(rowsep, rowstart, rowend) {
    if (!r || !cc) return("[]")
    w <- vapply(seq_len(cc), function(j) max(nchar(vapply(cells, `[`, "", j))), 0)
    paste(vapply(cells, function(row) paste0(rowstart, paste(mapply(function(z, k) formatC(z, width = k), row, w), collapse = ", "), rowend), ""),
          collapse = rowsep)
  }
  if (r <= 4L && cc <= 4L) {
    if (r == 1L) return(sprintf("vec%d%s", cc, tab("\n", "(", ")")))
    if (r == cc) return(sprintf("mat%d(%s)", r, tab(", ", "", "")))
    return(sprintf("mat%dx%d(%s)", cc, r, tab(", ", "", "")))
  }
  if (r == 1L || cc == 1L) return(sprintf("float[%d](%s)", r * cc, tab(",\n", "", "")))
  # GLSLPrinter.indent_code indents the lines inside the parentheses by 3
  sprintf("float[%d](\n%s\n) /* a %dx%d matrix */", r * cc, tab(",\n", "   ", ""), r, cc)
}

# MathMLContentPrinter._print_Integral (also used for Sum): nested applies,
# the last limit outermost
.mml_integral <- function(e, tag) {
  rec <- function(lims) {
    l <- lims[[1L]]$args
    out <- paste0("<apply><", tag, "/><bvar>", .mml_content(l[[1L]]), "</bvar>")
    if (length(l) == 3L) out <- paste0(out, "<lowlimit>", .mml_content(l[[2L]]), "</lowlimit><uplimit>", .mml_content(l[[3L]]), "</uplimit>")
    if (length(l) == 2L) out <- paste0(out, "<uplimit>", .mml_content(l[[2L]]), "</uplimit>")
    paste0(out, if (length(lims) == 1L) .mml_content(e$args[[1L]]) else rec(lims[-1L]), "</apply>")
  }
  rec(rev(e$args[-1L]))
}

# MathMLPresentationPrinter._print_MatrixBase with mat_delim "["
.mml_pres_matrix <- function(m) {
  rows <- vapply(.mrows(m), function(r) paste0("<mtr>", paste0("<mtd>", vapply(r, .mml_pres, ""), "</mtd>", collapse = ""), "</mtr>"), "")
  paste0("<mrow><mo>[</mo><mtable>", paste(rows, collapse = ""), "</mtable><mo>]</mo></mrow>")
}

# _get_latex_main of sympy/printing/preview.py
.latex_main <- function(expr, preamble = NULL, packages = character(), extra_preamble = NULL, euler = TRUE,
                        fontsize = NULL, ...) {
  if (is.null(preamble)) {
    pk <- c(packages, "amsmath", "amsfonts", if (euler) "euler")
    inc <- paste0("\n", paste0("\\usepackage{", pk, "}", collapse = "\n"))
    if (!is.null(extra_preamble)) inc <- paste0(inc, extra_preamble)
    fontsize <- if (is.null(fontsize)) "12pt" else if (is.numeric(fontsize)) paste0(fontsize, "pt") else fontsize
    preamble <- sprintf("\\documentclass[varwidth,%s]{standalone}\n%s\n\\begin{document}\n", fontsize, inc)
  } else if (length(packages) || !is.null(extra_preamble)) {
    stop("ValueError: The \"packages\" or \"extra_preamble\" keywordsmust not be set if a custom LaTeX preamble was specified",
         call. = FALSE)
  }
  body <- if (is.character(expr)) expr else paste0("$\\displaystyle ", latex(expr, mode = "plain", ...), "$")
  paste0(preamble, "\n", body, "\n\n\\end{document}")
}

#' Render an expression with LaTeX
#'
#' SymPy's `preview()`: typeset `expr` (an expression, or a LaTeX string)
#' in a standalone LaTeX document and compile it with
#' [tools::texi2pdf()]. `output = "pdf"` keeps the PDF; `"png"` converts it
#' with `pdftoppm` when that program is available; `"tex"` only writes the
#' source. With `viewer = "file"` the result goes to `filename`; otherwise
#' it is opened with the system viewer in an interactive session. A LaTeX
#' installation is needed for `"pdf"` and `"png"`.
#'
#' @param expr An expression or a LaTeX string.
#' @param output `"pdf"`, `"png"` or `"tex"`.
#' @param viewer `NULL`, or `"file"`.
#' @param euler Whether to use the Euler math fonts.
#' @param packages Extra LaTeX packages.
#' @param filename The output file for `viewer = "file"`.
#' @param preamble A full custom preamble.
#' @param extra_preamble Text added to the default preamble.
#' @param fontsize The font size, a number of points or a string.
#' @param ... Settings for [latex()].
#' @return The path of the output file, invisibly.
#' @examples
#' x <- Symbol("x")
#' f <- preview(x^2 + 1, output = "tex", viewer = "file", filename = tempfile(fileext = ".tex"))
#' cat(readLines(f), sep = "\n")
#' @export
preview <- function(expr, output = "pdf", viewer = NULL, euler = TRUE, packages = character(), filename = NULL,
                    preamble = NULL, extra_preamble = NULL, fontsize = NULL, ...) {
  output <- match.arg(output, c("pdf", "png", "tex"))
  src <- .latex_main(expr, preamble = preamble, packages = packages, extra_preamble = extra_preamble,
                     euler = euler, fontsize = fontsize, ...)
  if (identical(viewer, "file") && is.null(filename)) stop("ValueError: filename has to be specified if viewer=\"file\"", call. = FALSE)
  wd <- tempfile("preview"); dir.create(wd)
  tex <- file.path(wd, "texput.tex")
  writeLines(src, tex, useBytes = TRUE)
  out <- tex
  if (output != "tex") {
    old <- setwd(wd); on.exit(setwd(old), add = TRUE)
    tools::texi2pdf("texput.tex", clean = TRUE)
    out <- file.path(wd, "texput.pdf")
    if (output == "png") {
      if (!nzchar(Sys.which("pdftoppm"))) stop("RuntimeError: pdftoppm is needed for png output", call. = FALSE)
      system2("pdftoppm", c("-png", "-singlefile", "-r", "150", shQuote(out), shQuote(file.path(wd, "texput"))))
      out <- file.path(wd, "texput.png")
    }
  }
  if (identical(viewer, "file")) {
    file.copy(out, filename, overwrite = TRUE)
    return(invisible(filename))
  }
  if (interactive()) utils::browseURL(out)
  invisible(out)
}
