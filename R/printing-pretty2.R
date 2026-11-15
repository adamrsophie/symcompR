# Two dimensional pretty printing after sympy/printing/pretty: the
# stringPict and prettyForm boxes of stringpict.py, the symbology of
# pretty_symbology.py and the PrettyPrinter methods of pretty.py. A picture
# is a list of equal width lines with a baseline (0-based, as in SymPy) and
# a binding power.

.pp2_env <- new.env()
.pp2_env$u <- FALSE
.pp2_env$level <- 0L
.pp2_env$settings <- list()

.PP2 <- c(ATOM = 0L, FUNC = 1L, DIV = 2L, POW = 3L, MUL = 4L, ADD = 5L, NEG = 6L, OPEN = 7L, LINE = 8L)

.pp2_u <- function() isTRUE(.pp2_env$u)
.pp2_set <- function(name) .pp2_env$settings[[name]]

# line_width(): code points, combining marks not counted
.pp2_lw <- function(lines) {
  vapply(lines, function(l) {
    cp <- utf8ToInt(enc2utf8(l))
    if (anyNA(cp)) return(nchar(l))
    sum(!(cp %in% .pp2_combining))
  }, 0L, USE.NAMES = FALSE)
}

.pp2_center <- function(lines, w) {
  d <- w - .pp2_lw(lines)
  l <- d %/% 2L
  paste0(strrep(" ", l), lines, strrep(" ", d - l))
}

# a picture from lines that already have one width; a picture of empty
# lines loses one line, as joining and splitlines() do in SymPy
.pp2_mk <- function(lines, baseline = 0L, binding = 0L) {
  if (length(lines) > 1L && all(!nzchar(lines))) lines <- lines[-1L]
  if (!length(lines)) lines <- ""
  structure(list(picture = lines, baseline = as.integer(baseline), binding = as.integer(binding)), class = "pp2")
}

# stringPict(s) and prettyForm(s): split the string, center the lines
.pp2_pf <- function(s = "", baseline = 0L, binding = 0L) {
  lines <- if (length(s) == 1L) strsplit(s, "\n", fixed = TRUE)[[1L]] else s
  if (!length(lines)) lines <- ""
  w <- max(.pp2_lw(lines))
  structure(list(picture = .pp2_center(lines, w), baseline = as.integer(baseline), binding = as.integer(binding)),
            class = "pp2")
}

.pp2_as <- function(x) if (inherits(x, "pp2")) x else .pp2_pf(x)
.pp2_h <- function(p) length(p$picture)
.pp2_w <- function(p) .pp2_lw(p$picture[1L])
.pp2_str <- function(p) paste(p$picture, collapse = "\n")
.pp2_bind <- function(p, binding) { p$binding <- as.integer(binding); p }
.pp2_base <- function(p, baseline) { p$baseline <- as.integer(baseline); p }

# stringPict.next: side by side, on a common baseline
.pp2_next <- function(...) {
  objs <- lapply(list(...), .pp2_as)
  nb <- max(vapply(objs, function(o) o$baseline, 0L))
  below <- max(vapply(objs, function(o) .pp2_h(o) - o$baseline, 0L))
  nh <- nb + below
  pics <- lapply(objs, function(o) {
    e <- strrep(" ", .pp2_w(o))
    bp <- nb - o$baseline
    c(rep(e, bp), o$picture, rep(e, nh - .pp2_h(o) - bp))
  })
  .pp2_mk(do.call(paste0, pics), nb)
}
.pp2_right <- function(p, ...) .pp2_next(p, ...)
.pp2_left <- function(p, ...) .pp2_next(..., p)

.pp2_LINE <- structure(list(), class = "pp2_line")

# stringPict.stack: top to bottom, centered; the baseline of the second
.pp2_stack <- function(...) {
  objs <- lapply(list(...), function(a) if (inherits(a, "pp2_line")) a else .pp2_as(a))
  real <- Filter(function(o) !inherits(o, "pp2_line"), objs)
  nw <- max(vapply(real, .pp2_w, 0L))
  if (any(vapply(objs, inherits, NA, "pp2_line"))) {
    line <- .pp2_pf(.pp2_hobj("-", nw))
    objs <- lapply(objs, function(o) if (inherits(o, "pp2_line")) line else o)
  }
  lines <- unlist(lapply(objs, function(o) .pp2_center(o$picture, nw)))
  base <- if (length(objs) > 1L) .pp2_h(objs[[1L]]) + objs[[2L]]$baseline else 0L
  .pp2_mk(lines, base)
}
.pp2_below <- function(p, ...) .pp2_base(.pp2_stack(p, ...), p$baseline)
.pp2_above <- function(p, ...) {
  s <- .pp2_stack(..., p)
  .pp2_base(s, .pp2_h(s) - .pp2_h(p) + p$baseline)
}

.pp2_parens <- function(p, left = "(", right = ")", ifascii_nougly = FALSE) {
  h <- .pp2_h(p); b <- p$baseline
  if (ifascii_nougly && !.pp2_u()) { h <- 1L; b <- 0L }
  res <- p
  if (!is.null(left) && nzchar(left)) res <- .pp2_next(.pp2_pf(.pp2_vobj(left, h), b), p)
  if (!is.null(right) && nzchar(right)) res <- .pp2_next(res, .pp2_pf(.pp2_vobj(right, h), b))
  res
}

# pretty_symbology.xobj, vobj (a one piece object is split into its
# characters, as '\n'.join does) and hobj
.pp2_xobj <- function(symb, length) {
  if (length <= 0) stop("Length should be greater than 0", call. = FALSE)
  v <- (if (.pp2_u()) .pp2_xobj_u else .pp2_xobj_a)[[symb]]
  mid <- v$mid
  if (!is.null(mid)) { if (length %% 2L == 0L) length <- length + 1L } else mid <- v$ext
  if (length == 1L) return(v$c1)
  nxt <- (length - 2L) %/% 2L
  nmid <- (length - 2L) - nxt * 2L
  c(v$top, rep(v$ext, nxt), rep(mid, nmid), rep(v$ext, nxt), v$bot)
}
.pp2_vobj <- function(symb, h) {
  x <- .pp2_xobj(symb, h)
  if (length(x) == 1L) strsplit(x, "")[[1L]] else x
}
.pp2_hobj <- function(symb, w) paste(.pp2_xobj(symb, w), collapse = "")
.pp2_xsym <- function(s) .pp2_xsym_tab[[s]][if (.pp2_u()) 2L else 1L]

# pretty_atom: NULL stands for the KeyError of an ASCII atom
.pp2_atom <- function(name, default = NULL) {
  if (.pp2_u()) {
    if (name == "ImaginaryUnit" && identical(.pp2_set("imaginary_unit"), "j")) return("\U0001d55b")
    return(.pp2_atoms[[name]])
  }
  default
}

.pp2_center_accent <- function(s, accent) {
  mid <- nchar(s) %/% 2L + 1L
  paste0(substr(s, 1L, mid), accent, substring(s, mid + 1L))
}
.PP2_MODIFIERS <- list(
  mathring = function(s) .pp2_center_accent(s, "\u030a"), ddddot = function(s) .pp2_center_accent(s, "\u20dc"),
  dddot = function(s) .pp2_center_accent(s, "\u20db"), ddot = function(s) .pp2_center_accent(s, "\u0308"),
  dot = function(s) .pp2_center_accent(s, "\u0307"), check = function(s) .pp2_center_accent(s, "\u030c"),
  breve = function(s) .pp2_center_accent(s, "\u0306"), acute = function(s) .pp2_center_accent(s, "\u0301"),
  grave = function(s) .pp2_center_accent(s, "\u0300"), tilde = function(s) .pp2_center_accent(s, "\u0303"),
  hat = function(s) .pp2_center_accent(s, "\u0302"), bar = function(s) .pp2_center_accent(s, "\u0305"),
  vec = function(s) .pp2_center_accent(s, "\u20d7"), prime = function(s) paste0(s, "\u2032"),
  prm = function(s) paste0(s, "\u2032"), norm = function(s) paste0("\u2016", s, "\u2016"),
  avg = function(s) paste0("\u27e8", s, "\u27e9"), abs = function(s) paste0("|", s, "|"),
  mag = function(s) paste0("|", s, "|"))
.PP2_MOD_KEYS <- names(.PP2_MODIFIERS)[order(-nchar(names(.PP2_MODIFIERS)))]

.pp2_translate <- function(s, bold) {
  g <- if (bold) .pp2_greek_bold[[s]] else .pp2_greek[[s]]
  if (!is.null(g)) return(g)
  lo <- tolower(s)
  for (key in .PP2_MOD_KEYS)
    if (endsWith(lo, key) && nchar(s) > nchar(key))
      return(.PP2_MODIFIERS[[key]](.pp2_translate(substr(s, 1L, nchar(s) - nchar(key)), bold)))
  if (bold) return(paste(vapply(strsplit(s, "")[[1L]], function(ch) .pp2_bold[[ch]], ""), collapse = ""))
  s
}

# pretty_symbol: Greek names, modifiers, and Unicode sub and superscripts
.pp2_symbol <- function(name, bold = FALSE) {
  if (!.pp2_u()) return(name)
  p <- .split_super_sub(name)
  nm <- .pp2_translate(p$name, bold)
  plist <- function(l, mapping) {
    out <- character()
    for (s in l) {
      v <- mapping[[s]]
      if (is.null(v)) {
        chs <- strsplit(s, "")[[1L]]
        vs <- lapply(chs, function(ch) mapping[[ch]])
        if (!length(chs) || any(vapply(vs, is.null, NA))) return(NULL)
        v <- paste(unlist(vs), collapse = "")
      }
      out <- c(out, v)
    }
    out
  }
  psup <- plist(p$supers, .pp2_sup)
  psub <- if (!is.null(psup)) plist(p$subs, .pp2_sub)
  if (is.null(psub)) {
    if (length(p$subs)) nm <- paste0(nm, "_", paste(vapply(p$subs, .pp2_translate, "", bold = bold), collapse = "_"))
    if (length(p$supers)) nm <- paste0(nm, "__", paste(vapply(p$supers, .pp2_translate, "", bold = bold), collapse = "__"))
    return(nm)
  }
  paste0(nm, paste(psup, collapse = " "), paste(psub, collapse = " "))
}

# ---------------------------------------------------------------------------
# prettyForm arithmetic

.pp2_add <- function(...) {
  args <- list(...)
  res <- list()
  for (i in seq_along(args)) {
    a <- args[[i]]
    if (a$binding > .PP2[["NEG"]]) a <- .pp2_parens(a)
    if (i > 1L && a$binding != .PP2[["NEG"]]) res <- c(res, list(" + "))
    res <- c(res, list(a))
  }
  .pp2_bind(do.call(.pp2_next, res), .PP2[["ADD"]])
}

.pp2_div <- function(num, den) {
  if (num$binding == .PP2[["DIV"]]) num <- .pp2_parens(num)
  if (den$binding == .PP2[["DIV"]]) den <- .pp2_parens(den)
  if (num$binding == .PP2[["NEG"]]) num <- .pp2_right(num, " ")
  .pp2_bind(.pp2_stack(num, .pp2_LINE, den), .PP2[["DIV"]])
}

.pp2_mul <- function(...) {
  args <- list(...)
  if (length(args) == 1L) return(args[[1L]])
  wrap <- function(a) if (a$binding > .PP2[["MUL"]] && a$binding != .PP2[["NEG"]]) .pp2_parens(a) else a
  res <- list(wrap(args[[1L]]))
  star <- .pp2_xsym("*")
  for (a in args[-1L]) {
    if (a$picture[1L] != "\u00b0") res <- c(res, list(star))
    res <- c(res, list(wrap(a)))
  }
  is_str <- function(x, s) if (is.character(x)) identical(x, s) else identical(.pp2_str(x), s)
  i <- 1L
  while (i < length(res)) {
    if (is_str(res[[i]], "-1") && is_str(res[[i + 1L]], star)) {
      res <- c(res[seq_len(i - 1L)], list("-"), res[-seq_len(i + 1L)])
    }
    i <- i + 1L
  }
  first <- res[[1L]]
  first_line <- if (is.character(first)) first else first$picture[1L]
  if (startsWith(first_line, "-")) {
    bin <- .PP2[["NEG"]]
    if (is.character(first) && first == "-") {
      r <- res[[2L]]
      if (startsWith(r$picture[r$baseline + 1L], "-")) res[[1L]] <- "- "
    }
  } else bin <- .PP2[["MUL"]]
  .pp2_bind(do.call(.pp2_next, res), bin)
}

.pp2_pow <- function(a, b) {
  inline <- FALSE
  if (b$binding == .PP2[["POW"]]) b <- .pp2_parens(b)
  if (a$binding > .PP2[["FUNC"]]) a <- .pp2_parens(a)
  else if (a$binding == .PP2[["FUNC"]]) {
    if (.pp2_h(b) > 1L) a <- .pp2_parens(a) else inline <- TRUE
  }
  if (inline) {
    b$baseline <- a$prettyFunc$baseline + .pp2_h(b)
    func <- .pp2_right(a$prettyFunc, b)
    return(.pp2_right(func, a$prettyArgs))
  }
  top <- .pp2_left(b, strrep(" ", .pp2_w(a)))
  bot <- .pp2_right(a, strrep(" ", .pp2_w(b)))
  .pp2_bind(.pp2_above(bot, top), .PP2[["POW"]])
}

# ---------------------------------------------------------------------------
# PrettyPrinter

.pp2_print <- function(e) {
  .pp2_env$level <- .pp2_env$level + 1L
  on.exit(.pp2_env$level <- .pp2_env$level - 1L)
  if (!is_sym(e)) return(.pp2_nonsym(e))
  f <- .pp2_printers[[e$type]]
  if (!is.null(f)) return(f(e))
  if (e$type == "Function") return(.pp2_Function(e))
  .pp2_empty(e)
}

.pp2_empty <- function(e) .pp2_pf(if (is_sym(e)) .str(e) else paste(format(e), collapse = " "))

.pp2_nonsym <- function(e) {
  if (inherits(e, "pp2")) return(e)
  if (is.null(e)) return(.pp2_pf("None"))
  if (inherits(e, "symmatrix")) return(.pp2_MatrixBase(e))
  if (is.character(e) && length(e) == 1L) return(.pp2_pf(e))
  if (is.logical(e) && length(e) == 1L) return(.pp2_pf(if (e) "True" else "False"))
  if (is.numeric(e) && length(e) == 1L && !is.object(e)) return(.pp2_print(as_sym(e)))
  if (is.list(e) && !is.object(e)) return(.pp2_seq(e, "[", "]"))
  r <- tryCatch(as_sym(e), error = function(err) NULL)
  if (!is.null(r)) return(.pp2_print(r))
  .pp2_pf(paste(format(e), collapse = " "))
}

.pp2_seq <- function(seq, left = NULL, right = NULL, delimiter = ", ", parenthesize = function(x) FALSE,
                     ifascii_nougly = TRUE) {
  pf <- list()
  for (item in seq) {
    p <- .pp2_print(item)
    if (parenthesize(item)) p <- .pp2_parens(p)
    if (length(pf)) pf <- c(pf, list(delimiter))
    pf <- c(pf, list(p))
  }
  s <- if (!length(pf)) .pp2_pf("") else do.call(.pp2_next, pf)
  .pp2_parens(s, left, right, ifascii_nougly = ifascii_nougly)
}

.pp2_Symbol <- function(e, bold = FALSE) .pp2_pf(.pp2_symbol(e$name, bold))

.pp2_Float <- function(e) {
  fp <- .pp2_set("full_prec")
  if (identical(fp, "auto")) fp <- .pp2_env$level == 1L
  .pp2_pf(.fmt_float_sym(e, strip = !isTRUE(fp)))
}

.pp2_numer_denom <- function(p, q) {
  if (q == 1) return(if (p < 0) .pp2_pf(as.character(p), binding = .PP2[["NEG"]]) else .pp2_pf(as.character(p)))
  if (abs(p) >= 10 && abs(q) >= 10) {
    num <- if (p < 0) .pp2_pf(as.character(p), binding = .PP2[["NEG"]]) else .pp2_pf(as.character(p))
    return(.pp2_div(num, .pp2_pf(as.character(q))))
  }
  NULL
}

.pp2_Rational <- function(e) {
  r <- .pp2_numer_denom(gmp::numerator(e$value), gmp::denominator(e$value))
  if (is.null(r)) .pp2_empty(e) else r
}

.pp2_Constant <- function(e) {
  atom <- function(name, default = .str(e)) {
    a <- .pp2_atom(name)
    .pp2_pf(if (is.null(a)) default else a)
  }
  switch(EXPR = e$name,
    pi = atom("Pi"), E = atom("Exp1", "e"), I = atom("ImaginaryUnit"), oo = atom("Infinity"),
    `-oo` = atom("NegativeInfinity"),
    EulerGamma = if (.pp2_u()) .pp2_pf(.pp2_symbol("gamma")) else .pp2_pf("EulerGamma"),
    GoldenRatio = if (.pp2_u()) .pp2_pf(.pp2_symbol("phi")) else .pp2_pf("GoldenRatio"),
    Catalan = .pp2_pf("G"),
    .pp2_empty(e))
}

.pp2_SpecialSet <- function(e) {
  nm <- e$name
  if (nm == "Reals") {
    if (.pp2_u()) return(.pp2_pf(.pp2_atom("Reals")))
    return(.pp2_seq(list("-oo", "oo"), "(", ")"))
  }
  if (nm == "UniversalSet") return(.pp2_pf(if (.pp2_u()) .pp2_atom("Universe") else "UniversalSet"))
  a <- .pp2_atom(nm)
  .pp2_pf(if (is.null(a)) .str(e) else a)
}

# Add: negative terms become " - " and rational terms large fractions when
# another term is tall
.pp2_Add <- function(e) {
  ord <- .pp2_set("order")
  terms <- if (identical(ord, "rev-lex")) .add_print_order(e$args, rev = TRUE)
           else if (identical(ord, "none") || identical(ord, "old")) e$args else .add_print_order(e$args)
  pforms <- vector("list", length(terms))
  ind <- integer()
  negative <- function(p, i) {
    neg <- if (i == 1L) { if (.pp2_h(p) > 1L) "- " else "-" } else " - "
    q <- if (p$binding > .PP2[["NEG"]] || p$binding == .PP2[["ADD"]]) .pp2_parens(p) else p
    .pp2_bind(.pp2_next(neg, q), .PP2[["NEG"]])
  }
  for (i in seq_along(terms)) {
    t <- terms[[i]]
    if (t$type == "Mul" && could_extract_minus_sign(t)) {
      pforms[[i]] <- negative(.pp2_print(.neg_mul_noflat(t)), i)
    } else if (t$type == "Rational" && gmp::denominator(t$value) > 1) {
      ind <- c(ind, i)
    } else if (is_const(t, "-oo")) {
      pforms[[i]] <- negative(.pp2_print(S_Infinity), i)
    } else if (is_number_type(t) && .lp_is_number(t) && isTRUE(num_sign(t) < 0)) {
      pforms[[i]] <- negative(.pp2_print(.neg(t)), i)
    } else if (t$type %in% .PP2_REL || (t$type == "UnevaluatedExpr" && t$args[[1L]]$type == "Add")) {
      pforms[[i]] <- .pp2_parens(.pp2_print(t))
    } else pforms[[i]] <- .pp2_print(t)
  }
  if (length(ind)) {
    large <- any(vapply(pforms, function(p) !is.null(p) && .pp2_h(p) > 1L, NA))
    for (i in ind) {
      t <- terms[[i]]
      neg <- t$value < 0
      if (neg) t <- .neg(t)
      p <- if (large) .pp2_div(.pp2_pf(as.character(gmp::numerator(t$value))), .pp2_pf(as.character(gmp::denominator(t$value))))
           else .pp2_print(t)
      if (neg) p <- negative(p, i)
      pforms[[i]] <- p
    }
  }
  do.call(.pp2_add, pforms)
}

.pp2_Mul <- function(e) {
  args <- e$args
  if (is_one_num(args[[1L]]) || any(vapply(args[-1L], .lp_is_number, NA))) {
    strargs <- lapply(args, .pp2_print)
    negone <- identical(.pp2_str(strargs[[1L]]), "-1")
    if (negone) strargs[[1L]] <- .pp2_pf("1")
    obj <- do.call(.pp2_mul, strargs)
    if (negone) obj <- .pp2_pf(paste0("-", .pp2_str(obj)), obj$baseline, obj$binding)
    return(obj)
  }
  fs <- .mul_print_order(e$args)
  a <- list(); b <- list()
  for (it in fs) {
    if (it$type == "Pow" && .is_comm(it) && it$args[[2L]]$type %in% c("Integer", "Rational") && num_sign(it$args[[2L]]) < 0) {
      ex <- .neg(it$args[[2L]])
      b <- c(b, list(if (!is_one_num(ex)) new_sym("Pow", list(it$args[[1L]], ex)) else it$args[[1L]]))
    } else if (it$type %in% c("Integer", "Rational")) {
      p <- gmp::numerator(it$value); q <- gmp::denominator(it$value)
      if (p != 1) a <- c(a, list(Integer(as.character(p))))
      if (q != 1) b <- c(b, list(Integer(as.character(q))))
    } else a <- c(a, list(it))
  }
  addp <- function(z) z$type == "UnevaluatedExpr" && isTRUE(could_extract_minus_sign(z$args[[1L]]))
  pr <- function(z) if (addp(z)) .pp2_parens(.pp2_print(z)) else .pp2_print(z)
  a <- lapply(a, pr); b <- lapply(b, pr)
  if (!length(b)) return(do.call(.pp2_mul, a))
  if (!length(a)) a <- list(.pp2_print(S_One))
  .pp2_div(do.call(.pp2_mul, a), do.call(.pp2_mul, b))
}

.pp2_nth_root <- function(base, root) {
  bp <- .pp2_print(base)
  if (!identical(.pp2_set("use_unicode_sqrt_char"), FALSE) && .pp2_u() && is_one_num(.sub(root, Integer(1L))) &&
      .pp2_h(bp) == 1L && (.pp2_w(bp) == 1L || (base$type == "Integer" && base$value >= 0)))
    return(.pp2_left(bp, "\u221a"))
  zz <- .pp2_xobj("/", 1L)
  rootsign <- paste0(.pp2_xobj("\\", 1L), zz)
  rp <- .pp2_print(root)
  if (.pp2_h(rp) != 1L) return(.pp2_pow(.pp2_print(base), .pp2_print(.pow(root, S_NegativeOne))))
  is2 <- root$type == "Integer" && root$value == 2
  ex <- if (is2) "" else formatC(.pp2_str(rp), width = -2)
  if (nchar(ex) > 2L) rootsign <- paste0(strrep(" ", nchar(ex) - 2L), rootsign)
  rs <- .pp2_pf(paste0(ex, "\n", rootsign))
  rs$baseline <- 0L
  ll <- .pp2_h(bp) - 1L
  diag <- .pp2_pf(if (ll > 0L) vapply(seq_len(ll) - 1L, function(i) paste0(strrep(" ", ll - i - 1L), zz, strrep(" ", i)), "") else "")
  diag$baseline <- ll - 1L
  rs <- .pp2_right(rs, diag)
  rs$baseline <- max(1L, bp$baseline)
  s <- .pp2_pf(.pp2_hobj("_", 2L + .pp2_w(bp)))
  s <- .pp2_above(bp, s)
  .pp2_left(s, rs)
}

.pp2_Pow <- function(e) {
  b <- e$args[[1L]]; x <- e$args[[2L]]
  if (.is_comm(e)) {
    if (identical(x$key, S_NegativeOne$key)) return(.pp2_div(.pp2_pf("1"), .pp2_print(b)))
    nd <- fraction(x, exact = TRUE)
    n <- nd[[1L]]; d <- nd[[2L]]
    if (is_one_num(n) && .lp_is_atom(d) && x$type != "Integer" && (x$type == "Rational" || d$type == "Symbol") &&
        !identical(.pp2_set("root_notation"), FALSE))
      return(.pp2_nth_root(b, d))
    if (x$type %in% c("Integer", "Rational") && num_sign(x) < 0) return(.pp2_div(.pp2_pf("1"), .pp2_print(new_sym("Pow", list(b, .neg(x))))))
  }
  if (b$type %in% .PP2_REL) return(.pp2_pow(.pp2_parens(.pp2_print(b)), .pp2_print(x)))
  .pp2_pow(.pp2_print(b), .pp2_print(x))
}

.PP2_REL <- c("Equality", "Unequality", "StrictLessThan", "LessThan", "StrictGreaterThan", "GreaterThan")
.PP2_RELOP <- c(Equality = "==", Unequality = "!=", StrictLessThan = "<", LessThan = "<=", StrictGreaterThan = ">",
                GreaterThan = ">=")

.pp2_Relational <- function(e) {
  op <- .pp2_pf(paste0(" ", .pp2_xsym(.PP2_RELOP[[e$type]]), " "))
  .pp2_bind(.pp2_next(.pp2_print(e$args[[1L]]), op, .pp2_print(e$args[[2L]])), .PP2[["OPEN"]])
}

# _helper_print_function: the name as a Symbol, the arguments in parentheses
.pp2_helper_function <- function(name, args, sort = FALSE, func_name = NULL, delimiter = ", ", left = "(", right = ")",
                                 elementwise = FALSE) {
  if (sort) args <- .lp_dsk(args)
  if (is.null(func_name)) func_name <- name
  pfunc <- if (is.character(func_name)) .pp2_print(Symbol(func_name)) else .pp2_parens(.pp2_print(func_name))
  if (elementwise) pfunc <- .pp2_bind(.pp2_next(pfunc, .pp2_print(if (.pp2_u()) .pp2_atom("Modifier Letter Low Ring") else ".")), .PP2[["LINE"]])
  pargs <- .pp2_parens(.pp2_seq(args, delimiter = delimiter), left, right)
  p <- .pp2_bind(.pp2_next(pfunc, pargs), .PP2[["FUNC"]])
  p$prettyFunc <- pfunc
  p$prettyArgs <- pargs
  p
}

.pp2_fn <- function(e, sort = FALSE, func_name = NULL, left = "(", right = ")")
  .pp2_helper_function(if (is.null(e$name)) e$type else e$name, e$args, sort = sort, func_name = func_name, left = left, right = right)

.pp2_Function <- function(e) {
  f <- .pp2_funcs[[e$name]]
  if (!is.null(f)) return(f(e))
  .pp2_fn(e)
}

.pp2_fact <- function(suffix, left = FALSE) function(e) {
  x <- e$args[[1L]]
  p <- .pp2_print(x)
  if (!((x$type == "Integer" && x$value >= 0) || x$type == "Symbol")) p <- .pp2_parens(p)
  if (left) .pp2_left(p, suffix) else .pp2_right(p, suffix)
}

.pp2_number_function <- function(name) function(e) {
  p <- .pp2_pf(name)
  arg <- .pp2_print(e$args[[1L]])
  pa <- .pp2_below(.pp2_pf(strrep(" ", .pp2_w(arg))), arg)
  p <- .pp2_right(p, pa)
  if (length(e$args) == 1L) return(p)
  pargs <- .pp2_parens(.pp2_seq(list(e$args[[2L]])))
  out <- .pp2_bind(.pp2_next(p, pargs), .PP2[["FUNC"]])
  out$prettyFunc <- p
  out$prettyArgs <- pargs
  out
}

.pp2_vsep <- function(p1, p2, left = NULL, right = NULL, delimiter = "", ifascii_nougly = FALSE) {
  if (ifascii_nougly && !.pp2_u())
    return(.pp2_seq(list(p1, "|", p2), left = left, right = right, delimiter = delimiter, ifascii_nougly = TRUE))
  tmp <- .pp2_seq(list(p1, p2), left = left, right = right, delimiter = delimiter)
  sep <- .pp2_pf(.pp2_vobj("|", .pp2_h(tmp)), tmp$baseline)
  .pp2_seq(list(p1, sep, p2), left = left, right = right, delimiter = delimiter)
}

.pp2_greekname <- function(g, ascii) function(e) .pp2_fn(e, func_name = if (.pp2_u()) .pp2_greek[[g]] else ascii)

.pp2_funcs <- list(
  factorial = .pp2_fact("!"), factorial2 = .pp2_fact("!!"), subfactorial = .pp2_fact("!", left = TRUE),
  binomial = function(e) {
    n <- .pp2_print(e$args[[1L]]); k <- .pp2_print(e$args[[2L]])
    bar <- strrep(" ", max(.pp2_w(n), .pp2_w(k)))
    p <- .pp2_above(k, bar)
    p <- .pp2_above(p, n)
    p <- .pp2_parens(p, "(", ")")
    .pp2_base(p, (p$baseline + 1L) %/% 2L)
  },
  conjugate = function(e) { p <- .pp2_print(e$args[[1L]]); .pp2_above(p, .pp2_hobj("_", .pp2_w(p))) },
  Abs = function(e) .pp2_parens(.pp2_print(e$args[[1L]]), "|", "|"),
  floor = function(e) if (.pp2_u()) .pp2_parens(.pp2_print(e$args[[1L]]), "lfloor", "rfloor") else .pp2_fn(e),
  ceiling = function(e) if (.pp2_u()) .pp2_parens(.pp2_print(e$args[[1L]]), "lceil", "rceil") else .pp2_fn(e),
  exp = function(e) .pp2_pow(.pp2_pf(.pp2_atom("Exp1", "e")), .pp2_print(e$args[[1L]])),
  mathieuc = function(e) .pp2_fn(e, func_name = "C"), mathieus = function(e) .pp2_fn(e, func_name = "S"),
  mathieucprime = function(e) .pp2_fn(e, func_name = "C'"), mathieusprime = function(e) .pp2_fn(e, func_name = "S'"),
  lerchphi = .pp2_greekname("Phi", "lerchphi"), dirichlet_eta = .pp2_greekname("eta", "dirichlet_eta"),
  Heaviside = function(e) {
    fn <- if (.pp2_u()) .pp2_greek[["theta"]] else "Heaviside"
    if (length(e$args) > 1L && identical(e$args[[2L]]$key, S_Half$key))
      return(.pp2_left(.pp2_parens(.pp2_print(e$args[[1L]])), fn))
    .pp2_fn(e, func_name = fn)
  },
  fresnels = function(e) .pp2_fn(e, func_name = "S"), fresnelc = function(e) .pp2_fn(e, func_name = "C"),
  owens_t = function(e) .pp2_fn(e, func_name = "T"), airyai = function(e) .pp2_fn(e, func_name = "Ai"),
  airybi = function(e) .pp2_fn(e, func_name = "Bi"), airyaiprime = function(e) .pp2_fn(e, func_name = "Ai'"),
  airybiprime = function(e) .pp2_fn(e, func_name = "Bi'"), LambertW = function(e) .pp2_fn(e, func_name = "W"),
  beta = .pp2_greekname("Beta", "B"), betainc = function(e) .pp2_fn(e, func_name = "B'"),
  betainc_regularized = function(e) .pp2_fn(e, func_name = "I"),
  gamma = .pp2_greekname("Gamma", "Gamma"), uppergamma = .pp2_greekname("Gamma", "Gamma"),
  lowergamma = .pp2_greekname("gamma", "lowergamma"),
  DiracDelta = function(e) {
    if (!.pp2_u()) return(.pp2_fn(e))
    if (length(e$args) == 2L) {
      p <- .pp2_pow(.pp2_pf(.pp2_greek[["delta"]]), .pp2_parens(.pp2_print(e$args[[2L]])))
      p <- .pp2_right(p, " ")
      return(.pp2_right(p, .pp2_parens(.pp2_print(e$args[[1L]]))))
    }
    .pp2_left(.pp2_parens(.pp2_print(e$args[[1L]])), .pp2_greek[["delta"]])
  },
  polylog = function(e) {
    sub <- .pp2_str(.pp2_print(e$args[[1L]]))
    if (.pp2_u() && all(strsplit(sub, "")[[1L]] %in% names(.pp2_sub)))
      return(.pp2_fn(fcall(paste0("Li_", sub), list(e$args[[2L]]))))
    .pp2_fn(e)
  },
  expint = function(e) {
    sub <- .pp2_str(.pp2_print(e$args[[1L]]))
    if (.pp2_u() && all(strsplit(sub, "")[[1L]] %in% names(.pp2_sub)))
      return(.pp2_fn(fcall(paste0("E_", sub), list(e$args[[2L]]))))
    .pp2_fn(e)
  },
  Chi = function(e) {
    pfunc <- .pp2_pf("Chi")
    pargs <- .pp2_parens(.pp2_seq(e$args))
    p <- .pp2_bind(.pp2_next(pfunc, pargs), .PP2[["FUNC"]])
    p$prettyFunc <- pfunc
    p$prettyArgs <- pargs
    p
  },
  elliptic_e = function(e) {
    a0 <- .pp2_print(e$args[[1L]])
    p <- if (length(e$args) == 1L) a0 else .pp2_vsep(a0, .pp2_print(e$args[[2L]]))
    .pp2_left(.pp2_parens(p), "E")
  },
  elliptic_k = function(e) .pp2_left(.pp2_parens(.pp2_print(e$args[[1L]])), "K"),
  elliptic_f = function(e) .pp2_left(.pp2_parens(.pp2_vsep(.pp2_print(e$args[[1L]]), .pp2_print(e$args[[2L]]))), "F"),
  elliptic_pi = function(e) {
    nm <- if (.pp2_u()) .pp2_greek[["Pi"]] else "Pi"
    a0 <- .pp2_print(e$args[[1L]]); a1 <- .pp2_print(e$args[[2L]])
    p <- if (length(e$args) == 2L) .pp2_vsep(a0, a1)
         else .pp2_left(.pp2_left(.pp2_vsep(a1, .pp2_print(e$args[[3L]])), "; "), a0)
    .pp2_left(.pp2_parens(p), nm)
  },
  Mod = function(e) {
    p <- .pp2_print(e$args[[1L]])
    if (p$binding > .PP2[["MUL"]]) p <- .pp2_parens(p)
    p <- .pp2_right(.pp2_right(p, " mod "), .pp2_print(e$args[[2L]]))
    .pp2_bind(p, .PP2[["OPEN"]])
  },
  KroneckerDelta = function(e) {
    p <- .pp2_right(.pp2_right(.pp2_print(e$args[[1L]]), .pp2_pf(",")), .pp2_print(e$args[[2L]]))
    a <- .pp2_pf(if (.pp2_u()) .pp2_symbol("delta") else "d")
    top <- .pp2_left(p, strrep(" ", .pp2_w(a)))
    bot <- .pp2_right(a, strrep(" ", .pp2_w(p)))
    .pp2_bind(.pp2_below(bot, top), .PP2[["POW"]])
  },
  euler = .pp2_number_function("E"), catalan = .pp2_number_function("C"), bernoulli = .pp2_number_function("B"),
  bell = .pp2_number_function("B"), lucas = .pp2_number_function("L"), fibonacci = .pp2_number_function("F"),
  tribonacci = .pp2_number_function("T"),
  stieltjes = function(e) .pp2_number_function(if (.pp2_u()) .pp2_greek[["gamma"]] else "stieltjes")(e),
  primenu = function(e) .pp2_left(.pp2_parens(.pp2_print(e$args[[1L]])), if (.pp2_u()) .pp2_greek[["nu"]] else "nu"),
  primeomega = function(e) .pp2_left(.pp2_parens(.pp2_print(e$args[[1L]])), if (.pp2_u()) .pp2_greek[["Omega"]] else "Omega"),
  atan2 = function(e) .pp2_left(.pp2_parens(.pp2_seq(e$args)), "atan2"),
  SingularityFunction = function(e) {
    n <- .pp2_print(e$args[[3L]])
    shift <- .pp2_print(.sub(e$args[[1L]], e$args[[2L]]))
    if (.pp2_u()) return(.pp2_pow(.pp2_right(.pp2_right(.pp2_pf("<"), shift), ">"), n))
    .pp2_pow(.pp2_seq(list(shift), "<", ">", " "), n)
  }
)

# ---------------------------------------------------------------------------
# Calculus

.pp2_Derivative <- function(e) {
  ds_sym <- if (.lp_requires_partial(e$args[[1L]]) && .pp2_u()) "\u2202" else "d"
  x <- NULL
  total <- S_Zero
  for (vc in rev(e$args[-1L])) {
    v <- vc$args[[1L]]; num <- vc$args[[2L]]
    ds <- .pp2_left(.pp2_print(v), ds_sym)
    total <- .add(list(total, num))
    if (num$type != "Integer" || num$value > 1) ds <- .pp2_pow(ds, .pp2_pf(.str(num)))
    x <- if (is.null(x)) ds else .pp2_right(.pp2_right(x, " "), ds)
  }
  f <- .pp2_bind(.pp2_parens(.pp2_print(e$args[[1L]])), .PP2[["FUNC"]])
  p <- .pp2_pf(ds_sym)
  gt1 <- if (total$type == "Integer") total$value > 1 else TRUE
  if (gt1) p <- .pp2_pow(p, .pp2_pf(.str(total)))
  p <- .pp2_below(p, .pp2_LINE, x)
  p$baseline <- p$baseline + 1L
  p <- .pp2_next(p, f)
  .pp2_bind(p, .PP2[["MUL"]])
}

.pp2_Integral <- function(e) {
  f <- e$args[[1L]]
  lims <- e$args[-1L]
  pf <- .pp2_print(f)
  if (f$type == "Add") pf <- .pp2_parens(pf)
  arg <- pf
  for (l in lims) {
    pa <- .pp2_print(l$args[[1L]])
    if (.pp2_w(pa) > 1L) pa <- .pp2_parens(pa)
    arg <- .pp2_right(arg, " d", pa)
  }
  s <- NULL
  ascii <- !.pp2_u()
  for (l in lims) {
    h <- .pp2_h(arg)
    H <- h + 2L
    if (ascii) H <- H + 2L
    p <- .pp2_pf(.pp2_vobj("int", H))
    p$baseline <- arg$baseline + (H - h) %/% 2L
    la <- l$args
    if (length(la) > 1L) {
      if (length(la) == 2L) { pA <- .pp2_pf(""); pB <- .pp2_print(la[[2L]]) }
      else { pA <- .pp2_print(la[[2L]]); pB <- .pp2_print(la[[3L]]) }
      if (ascii) {
        pB <- .pp2_left(pB, strrep(" ", max(1L, 3L - .pp2_w(pB))))
        pA <- .pp2_right(pA, strrep(" ", max(1L, 4L - .pp2_w(pA))))
      }
      p <- .pp2_above(p, pB)
      p <- .pp2_below(p, pA)
    }
    if (!ascii) p <- .pp2_right(p, " ")
    s <- if (is.null(s)) p else .pp2_left(s, p)
  }
  .pp2_bind(.pp2_left(arg, s), .PP2[["MUL"]])
}

.pp2_sp_limits <- function(lim) {
  la <- lim$args
  lower <- .pp2_next(.pp2_print(la[[1L]]), .pp2_pf(paste0(" ", .pp2_xsym("=="), " ")), .pp2_print(la[[2L]]))
  list(lower = lower, upper = .pp2_print(la[[3L]]))
}

.pp2_Product <- function(e) {
  pfunc <- .pp2_print(e$args[[1L]])
  hc <- .pp2_xobj("_", 1L); cc <- hc; vc <- .pp2_xobj("|", 1L)
  if (.pp2_u()) { hc <- .pp2_xobj("-", 1L); cc <- .pp2_xobj("UpTack", 1L) }
  fh <- .pp2_h(pfunc)
  first <- TRUE; max_upper <- 0L; sign_height <- 0L
  for (lim in e$args[-1L]) {
    pl <- .pp2_sp_limits(lim)
    width <- ((fh + 2L) * 5L) %/% 3L - 2L
    lines <- c(paste0(hc, cc, strrep(hc, width - 2L), cc, hc),
               rep(paste0(" ", vc, strrep(" ", width - 2L), vc, " "), fh + 1L))
    sg <- do.call(.pp2_stack, as.list(lines))
    max_upper <- max(max_upper, .pp2_h(pl$upper))
    if (first) sign_height <- .pp2_h(sg)
    sg <- .pp2_below(.pp2_above(sg, pl$upper), pl$lower)
    if (first) { pfunc$baseline <- 0L; first <- FALSE }
    hh <- .pp2_h(sg)
    pad <- do.call(.pp2_stack, as.list(rep(" ", hh - 1L)))
    sg <- .pp2_right(sg, pad)
    pfunc <- .pp2_right(sg, pfunc)
  }
  pfunc$baseline <- max_upper + sign_height %/% 2L
  .pp2_bind(pfunc, .PP2[["MUL"]])
}

.pp2_asum <- function(hreq, ascii) {
  h <- max(hreq, 2L)
  d <- h %/% 2L
  w <- d + 1L
  more <- hreq %% 2L
  lines <- character()
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
  for (i in seq_len(d) - 1L) lines <- c(lines, paste0(strrep(" ", i), "\u2572", strrep(" ", w - i - 1L)))
  for (i in rev(seq_len(d) - 1L)) lines <- c(lines, paste0(strrep(" ", i), "\u2571", strrep(" ", w - i - 1L)))
  lines <- c(lines, strrep("\u203e", w))
  list(d = d, h = h + 2L * more, lines = lines, more = more)
}

.pp2_Sum <- function(e) {
  ascii <- !.pp2_u()
  f <- e$args[[1L]]
  pf <- .pp2_print(f)
  if (f$type == "Add") pf <- .pp2_parens(pf)
  H <- .pp2_h(pf) + 2L
  first <- TRUE; max_upper <- 0L; sign_height <- 0L; adj <- 0L
  for (lim in e$args[-1L]) {
    pl <- .pp2_sp_limits(lim)
    max_upper <- max(max_upper, .pp2_h(pl$upper))
    a <- .pp2_asum(H, ascii)
    sg <- do.call(.pp2_stack, as.list(a$lines))
    if (first) sign_height <- .pp2_h(sg)
    sg <- .pp2_below(.pp2_above(sg, pl$upper), pl$lower)
    if (first) {
      pf$baseline <- pf$baseline - (a$d - (.pp2_h(pf) %/% 2L - pf$baseline))
      first <- FALSE
    }
    pad <- do.call(.pp2_stack, as.list(rep(" ", a$h)))
    sg <- .pp2_right(sg, pad)
    pf <- .pp2_right(sg, pf)
    adj <- a$more
  }
  asc <- if (!adj) as.integer(ascii) else 0L
  pf$baseline <- max_upper + sign_height %/% 2L + asc
  .pp2_bind(pf, .PP2[["MUL"]])
}

.pp2_Limit <- function(e) {
  ex <- e$args[[1L]]; z <- e$args[[2L]]; z0 <- e$args[[3L]]
  dir <- if (length(e$args) >= 4L) e$args[[4L]]$name else "+"
  E <- .pp2_print(ex)
  if (.prec(ex) <= PREC$Mul) E <- .pp2_parens(E, "(", ")")
  Lim <- .pp2_pf("lim")
  la <- .pp2_print(z)
  la <- .pp2_right(la, if (.pp2_u()) paste0(.pp2_xobj("-", 1L), .pp2_atom("Arrow")) else "->")
  la <- .pp2_right(la, .pp2_print(z0))
  if (dir == "+-" || is_inf_const(z0)) dir <- ""
  else if (.pp2_u()) dir <- if (dir == "+") .pp2_atom("SuperscriptPlus") else .pp2_atom("SuperscriptMinus")
  la <- .pp2_right(la, .pp2_print(dir))
  Lim <- .pp2_below(Lim, la)
  .pp2_bind(.pp2_right(Lim, E), .PP2[["MUL"]])
}

.pp2_Order <- function(e) {
  p <- .pp2_print(e$args[[1L]])
  vp <- e$args[-1L]
  vars <- lapply(vp, function(t) t$args[[1L]])
  pts <- lapply(vp, function(t) t$args[[2L]])
  if ((length(pts) && any(vapply(pts, function(q) !is_zero_num(q), NA))) || length(vars) > 1L) {
    p <- .pp2_right(p, "; ")
    if (length(vars) > 1L) p <- .pp2_right(p, .pp2_print(do.call(Tuple, vars)))
    else if (length(vars)) p <- .pp2_right(p, .pp2_print(vars[[1L]]))
    p <- .pp2_right(p, if (.pp2_u()) paste0(" ", .pp2_atom("Arrow"), " ") else " -> ")
    p <- .pp2_right(p, if (length(pts) > 1L) .pp2_print(do.call(Tuple, pts)) else .pp2_print(pts[[1L]]))
  }
  .pp2_left(.pp2_parens(p), "O")
}

.pp2_Lambda <- function(e) {
  sig <- e$args[[1L]]
  arrow <- if (.pp2_u()) paste0(" ", .pp2_atom("ArrowFromBar"), " ") else " -> "
  if (sig$type == "Tuple" && length(sig$args) == 1L && sig$args[[1L]]$type == "Symbol") sig <- sig$args[[1L]]
  .pp2_bind(.pp2_next(.pp2_print(sig), arrow, .pp2_print(e$args[[2L]])), 8L)
}

.pp2_Subs <- function(e) {
  p <- .pp2_parens(.pp2_print(e$args[[1L]]))
  h <- if (.pp2_h(p) > 1L) .pp2_h(p) else 2L
  p <- .pp2_right(p, .pp2_pf(.pp2_vobj("|", h), p$baseline))
  b <- p$baseline
  p$baseline <- .pp2_h(p) - 1L
  tup <- function(a) if (a$type == "Tuple") a$args else list(a)
  vars <- tup(e$args[[2L]]); pts <- tup(e$args[[3L]])
  items <- lapply(seq_along(vars), function(i)
    .pp2_seq(list(.pp2_print(vars[[i]]), .pp2_xsym("=="), .pp2_print(pts[[i]])), delimiter = ""))
  p <- .pp2_right(p, .pp2_seq(items))
  p$baseline <- b
  p
}

# ---------------------------------------------------------------------------
# Logic

.pp2_boolean <- function(args, char) {
  pr <- function(a) {
    p <- .pp2_print(a)
    if (a$type %in% c(.LP_BOOL, .PP2_REL) && a$type != "Not") .pp2_parens(p) else p
  }
  p <- pr(args[[1L]])
  for (a in args[-1L]) p <- .pp2_right(.pp2_right(p, paste0(" ", char, " ")), pr(a))
  p
}
.pp2_logop <- function(atom, sort = TRUE) function(e, altchar = NULL) {
  if (!.pp2_u()) return(.pp2_fn(e, sort = sort))
  .pp2_boolean(if (sort) .lp_dsk(e$args) else e$args, if (is.null(altchar)) .pp2_atom(atom) else altchar)
}
.pp2_Not <- function(e) {
  if (!.pp2_u()) return(.pp2_helper_function("Not", e$args))
  a <- e$args[[1L]]
  if (a$type == "Equivalent") return(.pp2_printers$Equivalent(a, altchar = .pp2_atom("NotEquiv")))
  if (a$type == "Implies") return(.pp2_printers$Implies(a, altchar = .pp2_atom("NotArrow")))
  p <- .pp2_print(a)
  if (a$type %in% .LP_BOOL && a$type != "Not") p <- .pp2_parens(p)
  .pp2_left(p, .pp2_atom("Not"))
}

.pp2_Piecewise <- function(e) {
  n <- length(e$args)
  P <- lapply(e$args, function(ec) {
    cond <- ec$args[[2L]]
    list(.pp2_print(ec$args[[1L]]),
         if (identical(cond$key, S_true$key)) .pp2_pf("otherwise") else .pp2_right(.pp2_pf("for "), .pp2_print(cond)))
  })
  maxw <- vapply(1:2, function(j) max(vapply(P, function(r) .pp2_w(r[[j]]), 0L)), 0L)
  D <- NULL
  for (i in seq_len(n)) {
    row <- NULL
    for (j in 1:2) {
      p <- P[[i]][[j]]
      d <- maxw[j] - .pp2_w(p); wl <- d %/% 2L
      p <- .pp2_left(.pp2_right(p, strrep(" ", d - wl)), strrep(" ", wl))
      row <- if (is.null(row)) p else .pp2_right(.pp2_right(row, "  "), p)
    }
    D <- if (is.null(D)) row else .pp2_below(.pp2_below(D, " "), row)
  }
  D <- .pp2_parens(D, "{", "")
  D$baseline <- .pp2_h(D) %/% 2L
  .pp2_bind(D, .PP2[["OPEN"]])
}

# ---------------------------------------------------------------------------
# Containers and sets

.pp2_Tuple <- function(e) {
  a <- if (is_sym(e)) e$args else e
  if (length(a) == 1L) return(.pp2_parens(.pp2_next(.pp2_print(a[[1L]]), ","), "(", ")", ifascii_nougly = TRUE))
  .pp2_seq(a, "(", ")")
}

.pp2_Dict <- function(e) {
  keys <- lapply(e$args, function(t) t$args[[1L]]); vals <- lapply(e$args, function(t) t$args[[2L]])
  o <- if (length(keys) > 1L) .key_sort(lapply(keys, .sort_key), seq_along(keys)) else seq_along(keys)
  items <- lapply(unlist(o), function(i) .pp2_next(.pp2_print(keys[[i]]), ": ", .pp2_print(vals[[i]])))
  .pp2_seq(items, "{", "}")
}

.pp2_set_is <- function(s, types) is_sym(s) && s$type %in% types

.pp2_printers <- list(
  Symbol = function(e) .pp2_Symbol(e), Dummy = function(e) .pp2_Symbol(e), Wild = function(e) .pp2_Symbol(e),
  MatrixSymbol = function(e) .pp2_Symbol(e, identical(.pp2_set("mat_symbol_style"), "bold")),
  Float = .pp2_Float, Rational = .pp2_Rational, Integer = .pp2_Rational,
  Constant = .pp2_Constant, SpecialSet = .pp2_SpecialSet,
  EmptySet = function(e) .pp2_pf(if (.pp2_u()) .pp2_atom("EmptySet") else "EmptySet"),
  UniversalSet = function(e) .pp2_pf(if (.pp2_u()) .pp2_atom("Universe") else "UniversalSet"),
  Add = .pp2_Add, Mul = .pp2_Mul, Pow = .pp2_Pow,
  UnevaluatedExpr = function(e) .pp2_print(e$args[[1L]]),
  Equality = .pp2_Relational, Unequality = .pp2_Relational, StrictLessThan = .pp2_Relational,
  LessThan = .pp2_Relational, StrictGreaterThan = .pp2_Relational, GreaterThan = .pp2_Relational,
  Not = .pp2_Not, And = .pp2_logop("And"), Or = .pp2_logop("Or"), Xor = .pp2_logop("Xor"),
  Nand = .pp2_logop("Nand"), Nor = .pp2_logop("Nor"), Implies = .pp2_logop("Arrow", sort = FALSE),
  Equivalent = .pp2_logop("Equiv"),
  Piecewise = .pp2_Piecewise,
  Derivative = .pp2_Derivative, Integral = .pp2_Integral, Sum = .pp2_Sum, Product = .pp2_Product,
  Limit = .pp2_Limit, Order = .pp2_Order, Lambda = .pp2_Lambda, Subs = .pp2_Subs,
  Tuple = .pp2_Tuple, Dict = .pp2_Dict,
  FiniteSet = function(e) .pp2_seq(.lp_dsk(e$args), "{", "}", ", "),
  Interval = function(e) {
    a <- e$args
    if (identical(a[[1L]]$key, a[[2L]]$key)) return(.pp2_seq(a[1L], "{", "}"))
    .pp2_seq(a[1:2], if (identical(a[[3L]]$key, S_true$key)) "(" else "[", if (identical(a[[4L]]$key, S_true$key)) ")" else "]")
  },
  Union = function(e) .pp2_seq(e$args, NULL, NULL, paste0(" ", .pp2_atom("Union", "U"), " "),
                               parenthesize = function(s) .pp2_set_is(s, c("ProductSet", "Intersection", "Complement"))),
  Intersection = function(e) .pp2_seq(e$args, NULL, NULL, paste0(" ", .pp2_atom("Intersection", "n"), " "),
                                      parenthesize = function(s) .pp2_set_is(s, c("ProductSet", "Union", "Complement"))),
  Complement = function(e) .pp2_seq(e$args, NULL, NULL, " \\ ",
                                    parenthesize = function(s) .pp2_set_is(s, c("ProductSet", "Intersection", "Union"))),
  SymmetricDifference = function(e) {
    if (!.pp2_u()) stop(.sympy_condition("NotImplementedError", "ASCII pretty printing of SymmetricDifference is not implemented"))
    .pp2_seq(e$args, NULL, NULL, paste0(" ", .pp2_atom("SymmetricDifference"), " "))
  },
  ProductSet = function(e) {
    s <- e$args
    if (length(s) >= 1L && length(unique(vapply(s, function(z) z$key, ""))) == 1L)
      return(.pp2_pow(.pp2_print(s[[1L]]), .pp2_print(length(s))))
    ch <- if (.pp2_u()) .pp2_atom("Multiplication") else "x"
    .pp2_seq(s, NULL, NULL, paste0(" ", ch, " "),
             parenthesize = function(z) .pp2_set_is(z, c("Union", "Intersection", "ProductSet")))
  },
  ImageSet = function(e) {
    inn <- if (.pp2_u()) .pp2_atom("SmallElementOf") else "in"
    lam <- e$args[[1L]]
    sets <- e$args[-1L]
    sig <- lam$args[[1L]]
    sig <- if (sig$type == "Tuple") sig$args else list(sig)
    ex <- .pp2_print(lam$args[[2L]])
    if (length(sig) == 1L) {
      S <- .pp2_seq(list(sig[[1L]], inn, sets[[1L]]), delimiter = " ")
    } else {
      pargs <- unlist(lapply(seq_along(sig), function(i) list(sig[[i]], " ", inn, " ", sets[[i]], ", ")), recursive = FALSE)
      S <- .pp2_seq(pargs[-length(pargs)], delimiter = "")
    }
    .pp2_vsep(ex, S, left = "{", right = "}", ifascii_nougly = TRUE, delimiter = " ")
  },
  ConditionSet = function(e) {
    inn <- if (.pp2_u()) .pp2_atom("SmallElementOf") else "in"
    and <- if (.pp2_u()) .pp2_atom("And") else "and"
    sym <- e$args[[1L]]
    variables <- .pp2_seq(if (sym$type == "Tuple") list(sym) else list(sym))
    cond <- .pp2_print(e$args[[2L]])
    if (.pp2_u()) cond <- .pp2_parens(cond)
    base <- e$args[[3L]]
    if (base$type %in% c("UniversalSet") || (base$type == "SpecialSet" && base$name == "UniversalSet"))
      return(.pp2_vsep(variables, cond, left = "{", right = "}", ifascii_nougly = TRUE, delimiter = " "))
    C <- .pp2_seq(list(variables, inn, .pp2_print(base), and, cond), delimiter = " ")
    .pp2_vsep(variables, C, left = "{", right = "}", ifascii_nougly = TRUE, delimiter = " ")
  },
  Range = function(e) {
    dots <- if (.pp2_u()) .pp2_atom("Dots") else "..."
    p <- .range_parts(e)
    ps <- if (is_inf_const(p$start) && is_inf_const(p$stop)) {
      if (num_sign(p$step) > 0) list(dots, -1L, 0L, 1L, dots) else list(dots, 1L, 0L, -1L, dots)
    } else if (is_inf_const(p$start)) list(dots, .sub(.range_last(e), p$step), .range_last(e))
    else if (is_inf_const(p$stop)) list(p$start, .add(list(p$start, p$step)), dots)
    else {
      n <- .int_val(.range_size(e))
      if (n > 4) list(.range_item(e, 0L), .range_item(e, 1L), dots, .range_last(e))
      else lapply(seq_len(n) - 1L, function(i) .range_item(e, i))
    }
    .pp2_seq(ps, "{", "}", ", ")
  },
  Contains = function(e) {
    if (.pp2_u()) return(.pp2_bind(.pp2_next(.pp2_print(e$args[[1L]]), paste0(" ", .pp2_atom("ElementOf"), " "),
                                             .pp2_print(e$args[[2L]])), 8L))
    .pp2_pf(.str(e))
  },
  # matrices
  MatAdd = function(e) {
    s <- NULL
    for (it in e$args) {
      p <- .pp2_print(it)
      if (is.null(s)) { s <- p; next }
      c1 <- if (it$type == "MatMul") it$args[[1L]] else S_One
      neg <- (is_exact_number(c1) && isTRUE(num_sign(c1) < 0)) || (c1$type == "Mul" && isTRUE(could_extract_minus_sign(c1)))
      s <- .pp2_next(s, if (neg) " " else " + ")
      s <- .pp2_next(s, p)
    }
    s
  },
  MatMul = function(e) {
    n <- length(e$args)
    a <- lapply(e$args, function(z) if (z$type %in% c("Add", "MatAdd", "HadamardProduct", "KroneckerProduct") && n > 1L)
      .pp2_parens(.pp2_print(z)) else .pp2_print(z))
    do.call(.pp2_mul, a)
  },
  MatPow = function(e) {
    b <- e$args[[1L]]
    p <- .pp2_print(b)
    if (b$type != "MatrixSymbol" && .lp_is_matexpr(b)) p <- .pp2_parens(p)
    .pp2_pow(p, .pp2_print(e$args[[2L]]))
  },
  Transpose = function(e) {
    m <- e$args[[1L]]
    p <- .pp2_print(m)
    if (!(m$type %in% c("MatrixSymbol", "BlockMatrix")) && .lp_is_matexpr(m)) p <- .pp2_parens(p)
    .pp2_pow(p, .pp2_pf("T"))
  },
  Adjoint = function(e) {
    m <- e$args[[1L]]
    p <- .pp2_print(m)
    dag <- .pp2_pf(if (.pp2_u()) .pp2_atom("Dagger") else "+")
    if (!(m$type %in% c("MatrixSymbol", "BlockMatrix")) && .lp_is_matexpr(m)) p <- .pp2_parens(p)
    .pp2_pow(p, dag)
  },
  Identity = function(e) .pp2_pf(if (.pp2_u()) .pp2_atom("IdentityMatrix") else "I"),
  ZeroMatrix = function(e) .pp2_pf(if (.pp2_u()) .pp2_atom("ZeroMatrix") else "0"),
  OneMatrix = function(e) .pp2_pf(if (.pp2_u()) .pp2_atom("OneMatrix") else "1"),
  Trace = function(e) {
    D <- .pp2_parens(.pp2_print(e$args[[1L]]), "(", ")")
    D$baseline <- .pp2_h(D) %/% 2L
    .pp2_left(D, "tr")
  },
  Determinant = function(e) {
    m <- e$args[[1L]]
    if (m$type == "BlockMatrix") {
      b <- .lp_block_rows(m)
      D <- .pp2_matrix_contents(b$rows, b$ncol)
    } else D <- .pp2_print(m)
    D$baseline <- .pp2_h(D) %/% 2L
    .pp2_parens(D, "|", "|")
  },
  MatrixElement = function(e) {
    par <- e$args[[1L]]
    i <- e$args[[2L]]; j <- e$args[[3L]]
    if (par$type == "MatrixSymbol" && i$type == "Integer" && j$type == "Integer")
      return(.pp2_print(Symbol(paste0(par$name, "_", .str(i), .str(j)))))
    pfunc <- .pp2_parens(.pp2_print(par))
    pidx <- .pp2_parens(.pp2_seq(list(i, j), delimiter = ", "), "[", "]")
    p <- .pp2_bind(.pp2_next(pfunc, pidx), .PP2[["FUNC"]])
    p$prettyFunc <- pfunc
    p$prettyArgs <- pidx
    p
  },
  HadamardProduct = function(e) .pp2_seq(e$args, NULL, NULL, if (.pp2_u()) .pp2_atom("Ring") else ".*",
                                         parenthesize = function(z) z$type %in% c("MatAdd", "MatMul", "HadamardProduct")),
  KroneckerProduct = function(e) .pp2_seq(e$args, NULL, NULL, if (.pp2_u()) paste0(" ", .pp2_atom("TensorProduct"), " ") else " x ",
                                          parenthesize = function(z) z$type %in% c("MatAdd", "MatMul")),
  DotProduct = function(e) do.call(.pp2_mul, lapply(e$args, .pp2_print)),
  Covariance = function(e) .pp2_fn(e, func_name = "Cov"), Variance = function(e) .pp2_fn(e, func_name = "Var"),
  Probability = function(e) .pp2_fn(e, func_name = "P"),
  Expectation = function(e) .pp2_fn(e, func_name = "E", left = "[", right = "]"),
  CRootOf = function(e) .pp2_left(.pp2_parens(.pp2_seq(list(.pp2_Add_lex(e$args[[1L]]), e$args[[2L]]))), "CRootOf"),
  ComplexRootOf = function(e) .pp2_left(.pp2_parens(.pp2_seq(list(.pp2_Add_lex(e$args[[1L]]), e$args[[2L]]))), "CRootOf"),
  Str = function(e) .pp2_pf(e$name),
  Tr = function(e) {
    D <- .pp2_parens(.pp2_print(e$args[[1L]]), "(", ")")
    D$baseline <- .pp2_h(D) %/% 2L
    .pp2_left(D, "Tr")
  },
  RootSum = function(e) {
    a <- list(.pp2_Add_lex(e$args[[1L]]))
    if (length(e$args) > 1L) {
      fn <- e$args[[2L]]
      ident <- fn$type == "Lambda" && length(fn$args[[1L]]$args) == 1L && identical(fn$args[[1L]]$args[[1L]]$key, fn$args[[2L]]$key)
      if (!ident) a <- c(a, list(.pp2_print(fn)))
    }
    .pp2_left(.pp2_parens(.pp2_seq(a)), "RootSum")
  },
  SeqFormula = function(e) .pp2_SeqFormula(e), SeqPer = function(e) .pp2_SeqFormula(e),
  SeqAdd = function(e) .pp2_SeqFormula(e), SeqMul = function(e) .pp2_SeqFormula(e),
  Inverse = function(e) {
    m <- e$args[[1L]]
    p <- .pp2_print(m)
    if (m$type != "MatrixSymbol" && .lp_is_matexpr(m)) p <- .pp2_parens(p)
    .pp2_pow(p, .pp2_print(S_NegativeOne))
  },
  BlockMatrix = function(e) {
    b <- .lp_block_rows(e)
    if (length(b$rows) == 1L && b$ncol == 1L) return(.pp2_print(b$rows[[1L]][[1L]]))
    D <- .pp2_matrix_contents(b$rows, b$ncol)
    D$baseline <- .pp2_h(D) %/% 2L
    .pp2_parens(D, "[", "]")
  },
  ITE = function(e) .pp2_print(Piecewise(list(e$args[[2L]], e$args[[1L]]), list(e$args[[3L]], S_true))),
  Assignment = function(e) .pp2_assign(e, ":="), AddAugmentedAssignment = function(e) .pp2_assign(e, "+="),
  SubAugmentedAssignment = function(e) .pp2_assign(e, "-="), MulAugmentedAssignment = function(e) .pp2_assign(e, "*="),
  DivAugmentedAssignment = function(e) .pp2_assign(e, "/="), ModAugmentedAssignment = function(e) .pp2_assign(e, "%="),
  BooleanAtom = function(e) .pp2_pf(if (identical(e$name, "true")) "True" else "False")
)

.pp2_SeqFormula <- function(s) {
  dots <- if (.pp2_u()) .pp2_atom("Dots") else "..."
  b <- .seq_bounds(s)
  if (length(free_symbols(b[[1L]])) || length(free_symbols(b[[2L]]))) return(.pp2_fn(s))
  ps <- if (is_const(b[[1L]], "-oo")) {
    c(list(dots), lapply(3:0, function(k) .seq_coeff(s, .sub(b[[2L]], Integer(k)))))
  } else if (is_inf_const(b[[2L]]) || isTRUE(.int_val(.seq_length(s)) > 4)) c(sequence_terms(s, 4L), list(dots))
  else sequence_terms(s, 1e6)
  .pp2_seq(ps, "[", "]")
}

.pp2_assign <- function(e, op)
  .pp2_next(.pp2_print(e$args[[1L]]), .pp2_pf(paste0(" ", .pp2_xsym(op), " ")), .pp2_print(e$args[[2L]]))

.pp2_Add_lex <- function(e) if (e$type == "Add") .pp2_Add(e) else .pp2_print(e)

# _print_matrix_contents and _print_MatrixBase
.pp2_matrix_contents <- function(rows, ncol) {
  if (!length(rows) || !ncol) return(.pp2_pf(""))
  Ms <- lapply(rows, function(r) lapply(r, .pp2_print))
  maxw <- vapply(seq_len(ncol), function(j) max(vapply(Ms, function(r) .pp2_w(r[[j]]), 0L)), 0L)
  D <- NULL
  for (r in Ms) {
    row <- NULL
    for (j in seq_len(ncol)) {
      s <- r[[j]]
      d <- maxw[j] - .pp2_w(s); l <- d %/% 2L
      s <- .pp2_left(.pp2_right(s, strrep(" ", d - l)), strrep(" ", l))
      row <- if (is.null(row)) s else .pp2_right(.pp2_right(row, "  "), s)
    }
    D <- if (is.null(D)) row else .pp2_below(.pp2_below(D, " "), row)
  }
  D
}

.pp2_MatrixBase <- function(M, lp = "[", rp = "]") {
  D <- .pp2_matrix_contents(if (M$nrow) .mrows(M) else list(), M$ncol)
  D$baseline <- .pp2_h(D) %/% 2L
  .pp2_parens(D, lp, rp)
}

.pp2_printers$ImmutableMatrix <- function(e) {
  m <- .lp_immutable_rows(e)
  D <- .pp2_matrix_contents(m$rows, m$ncol)
  D$baseline <- .pp2_h(D) %/% 2L
  .pp2_parens(D, "[", "]")
}

# PrettyPrinter.doprint with wrap_line False
.pp2_doprint <- function(e, use_unicode, settings = list()) {
  old <- list(u = .pp2_env$u, level = .pp2_env$level, settings = .pp2_env$settings)
  on.exit({ .pp2_env$u <- old$u; .pp2_env$level <- old$level; .pp2_env$settings <- old$settings })
  .pp2_env$u <- isTRUE(use_unicode)
  .pp2_env$level <- 0L
  defaults <- list(order = NULL, full_prec = "auto", use_unicode_sqrt_char = TRUE, root_notation = TRUE,
                   mat_symbol_style = "plain", imaginary_unit = "i", perm_cyclic = TRUE)
  .pp2_env$settings <- utils::modifyList(defaults, settings)
  p <- .pp2_print(e)
  .pp2_str(p)
}

# ---------------------------------------------------------------------------
# Hypergeometric and Meijer G functions

# str.splitlines(): a final newline adds no empty line
.pp2_splitlines <- function(s) {
  if (endsWith(s, "\n")) s <- substr(s, 1L, nchar(s) - 1L)
  l <- strsplit(paste0(s, "\n."), "\n", fixed = TRUE)[[1L]]
  l[-length(l)]
}

.pp2_annotated <- function(letter) {
  if (.pp2_u()) {
    if (letter == "F") return(list(sz = 2L, t = 0L, b = 2L, add = 0L,
                                   img = "\u250c\u2500\n\u251c\u2500\n\u2575"))
    return(list(sz = 3L, t = 0L, b = 3L, add = 1L,
                img = "\u256d\u2500\u256e\n\u2502\u2576\u2510\n\u2570\u2500\u256f"))
  }
  if (letter == "F") return(list(sz = 3L, t = 0L, b = 3L, add = 0L, img = " _\n|_\n|\n"))
  list(sz = 3L, t = 0L, b = 3L, add = 1L, img = " __\n/__\n\\_|")
}

.pp2_hvec <- function(v) {
  if (!length(v)) return(.pp2_pf(" "))
  D <- v[[1L]]
  for (p in v[-1L]) D <- .pp2_right(.pp2_right(D, ", "), p)
  D
}

.pp2_letter <- function(letter, D) {
  above <- .pp2_h(D) %/% 2L - 1L
  below <- .pp2_h(D) - above - 1L
  an <- .pp2_annotated(letter)
  s <- paste0(strrep("\n", max(0L, above - an$t)), an$img, strrep("\n", max(0L, below - an$b)))
  lines <- .pp2_splitlines(s)
  w <- max(.pp2_lw(lines))
  F <- .pp2_mk(paste0(lines, strrep(" ", w - .pp2_lw(lines))), above + an$sz)
  list(F = F, above = above, add = (an$sz + 1L) %/% 2L, annot_add = an$add)
}

.pp2_hyper_body <- function(D, z) {
  P <- .pp2_print(z)
  P$baseline <- .pp2_h(P) %/% 2L
  D$baseline <- .pp2_h(D) %/% 2L
  P <- .pp2_left(P, " ")
  D <- .pp2_right(D, " ")
  D <- .pp2_vsep(D, P)
  .pp2_parens(D, "(", ")")
}

.pp2_hyper <- function(e) {
  ap <- lapply(.lp_tuple_args(e$args[[1L]]), .pp2_print)
  bq <- lapply(.lp_tuple_args(e$args[[2L]]), .pp2_print)
  D <- .pp2_below(.pp2_below(.pp2_hvec(ap), " "), .pp2_hvec(bq))
  D <- .pp2_hyper_body(D, e$args[[3L]])
  L <- .pp2_letter("F", D)
  F <- .pp2_left(L$F, .pp2_print(length(ap)))
  F <- .pp2_right(F, .pp2_print(length(bq)))
  F$baseline <- L$above + L$add
  .pp2_right(F, " ", D)
}

.pp2_meijerg <- function(e) {
  a <- .lp_tuple_args(e$args[[1L]]); b <- .lp_tuple_args(e$args[[2L]])
  pr <- function(t) lapply(.lp_tuple_args(t), .pp2_print)
  v <- list(pr(a[[1L]]), pr(a[[2L]]), pr(b[[1L]]), pr(b[[2L]]))
  vp <- lapply(v, .pp2_hvec)
  for (i in 1:2) {
    idx <- c(i, i + 2L)
    maxw <- max(.pp2_w(vp[[idx[1L]]]), .pp2_w(vp[[idx[2L]]]))
    for (k in idx) {
      l <- (maxw - .pp2_w(vp[[k]])) %/% 2L
      r <- maxw - l - .pp2_w(vp[[k]])
      vp[[k]] <- .pp2_right(.pp2_left(vp[[k]], strrep(" ", l)), strrep(" ", r))
    }
  }
  D1 <- .pp2_below(.pp2_right(vp[[1L]], "  ", vp[[2L]]), " ")
  D2 <- .pp2_right(vp[[3L]], "  ", vp[[4L]])
  D <- .pp2_hyper_body(.pp2_below(D1, D2), e$args[[3L]])
  an <- .pp2_annotated("G")
  above <- .pp2_h(D) %/% 2L - 1L
  below <- .pp2_h(D) - above - 1L
  lines <- .pp2_splitlines(paste0(strrep("\n", max(0L, above - an$t)), an$img, strrep("\n", max(0L, below - an$b))))
  w <- max(.pp2_lw(lines))
  F <- .pp2_mk(paste0(lines, strrep(" ", w - .pp2_lw(lines))), above + an$sz)
  pp <- .pp2_print(length(a[[1L]]$args) + length(a[[2L]]$args))
  pq <- .pp2_print(length(b[[1L]]$args) + length(b[[2L]]$args))
  pm <- .pp2_print(length(v[[3L]])); pn <- .pp2_print(length(v[[1L]]))
  adjust <- function(p1, p2) {
    d <- .pp2_w(p1) - .pp2_w(p2)
    if (d > 0) list(p1, .pp2_left(p2, strrep(" ", d))) else if (d < 0) list(.pp2_left(p1, strrep(" ", -d)), p2) else list(p1, p2)
  }
  t1 <- adjust(pp, pm); pp <- t1[[1L]]; pm <- t1[[2L]]
  t2 <- adjust(pq, pn); pq <- t2[[1L]]; pn <- t2[[2L]]
  pu <- .pp2_right(pm, ", ", pn)
  pl <- .pp2_right(pp, ", ", pq)
  ht <- F$baseline - above - 2L
  if (ht > 0) pu <- .pp2_below(pu, strrep("\n", ht))
  p <- .pp2_below(pu, pl)
  F$baseline <- above
  F <- .pp2_right(F, p)
  F$baseline <- above + an$add
  .pp2_right(F, " ", D)
}
.pp2_funcs$hyper <- .pp2_hyper
.pp2_funcs$meijerg <- .pp2_meijerg

# ---------------------------------------------------------------------------
# The helpers of sympy.printing.pretty.pretty_symbology

.pp2_with_u <- function(u, expr) {
  old <- .pp2_env$u
  .pp2_env$u <- isTRUE(u)
  on.exit(.pp2_env$u <- old)
  expr
}

#' Pretty printing symbology
#'
#' The helpers of `sympy.printing.pretty.pretty_symbology` that the pretty
#' printer draws with. `pretty_symbol(name)` renders a symbol name with
#' Greek letters and Unicode sub- and superscripts; `xsym(s)` gives the
#' character for an operator such as `"=="` or `"*"`; `hobj(symb, width)`
#' and `vobj(symb, height)` draw a horizontal or vertical object such as
#' `"-"` or `"("`; `center_pad(wstring, wtarget)` the padding that centers
#' a string, `center(string, width)` the centered string; `annotated(letter)`
#' the stylized `"F"` or `"G"` of hypergeometric and Meijer G functions as
#' `list(size, top, bottom, add, picture)`. `greek_unicode` maps Greek
#' letter names to their characters and `sup` characters to their
#' superscripts.
#'
#' @param name A symbol name.
#' @param bold_name Whether to use bold letters.
#' @param s An operator name.
#' @param symb An object name.
#' @param width,height The size.
#' @param wstring,wtarget The string width and the target width.
#' @param string A string.
#' @param fillchar The padding character.
#' @param letter `"F"` or `"G"`.
#' @param use_unicode Whether to use Unicode characters.
#' @return A string, a pair of strings, or a list.
#' @examples
#' pretty_symbol("alpha_1")
#' xsym("==")
#' hobj("-", 5)
#' cat(vobj("(", 3), "\n")
#' center("ab", 6)
#' @name pretty_symbology
#' @export
pretty_symbol <- function(name, bold_name = FALSE, use_unicode = TRUE) .pp2_with_u(use_unicode, .pp2_symbol(name, bold_name))
#' @rdname pretty_symbology
#' @export
xsym <- function(s, use_unicode = TRUE) .pp2_with_u(use_unicode, .pp2_xsym(s))
#' @rdname pretty_symbology
#' @export
hobj <- function(symb, width, use_unicode = TRUE) .pp2_with_u(use_unicode, .pp2_hobj(symb, as.integer(width)))
#' @rdname pretty_symbology
#' @export
vobj <- function(symb, height, use_unicode = TRUE) .pp2_with_u(use_unicode, paste(.pp2_vobj(symb, as.integer(height)), collapse = "\n"))
#' @rdname pretty_symbology
#' @export
center_pad <- function(wstring, wtarget, fillchar = " ") {
  if (wstring > wtarget) stop("ValueError: not enough space for string", call. = FALSE)
  d <- wtarget - wstring; l <- d %/% 2
  c(strrep(fillchar, l), strrep(fillchar, d - l))
}
#' @rdname pretty_symbology
#' @export
center <- function(string, width, fillchar = " ") {
  p <- center_pad(.pp2_lw(string), width, fillchar)
  paste0(p[1L], string, p[2L])
}
#' @rdname pretty_symbology
#' @export
annotated <- function(letter, use_unicode = TRUE) {
  a <- .pp2_with_u(use_unicode, .pp2_annotated(letter))
  list(a$sz, a$t, a$b, a$add, a$img)
}
#' @rdname pretty_symbology
#' @format `greek_unicode` and `sup` are named lists of strings.
#' @export
greek_unicode <- .pp2_greek
#' @rdname pretty_symbology
#' @export
sup <- .pp2_sup
