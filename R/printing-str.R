# String printer, after sympy/printing/str.py. Output uses R syntax (^ for
# powers) so it can be read back with sympify().

PREC <- list(Lambda = 1, Xor = 10, Or = 20, And = 30, Relational = 35,
             BitwiseOr = 36, BitwiseXor = 37, BitwiseAnd = 38,
             Add = 40, Mul = 50, Pow = 60, Func = 70, Not = 100, Atom = 1000)


.prec <- function(e) {
  t <- e$type
  switch(t,
         Add = PREC$Add,
         Mul = if (any(vapply(e$args, function(a) a$type == "UnevaluatedExpr", TRUE))) PREC$Mul
               else if (could_extract_minus_sign(e)) PREC$Add else PREC$Mul,
         Pow = PREC$Pow,
         Rational = if (e$value < 0) PREC$Add else PREC$Mul,
         Integer = if (e$value < 0) PREC$Add else PREC$Atom,
         Float = if (e$value < 0) PREC$Add else PREC$Atom,
         Constant = if (e$name == "-oo") PREC$Add else PREC$Atom,
         StrictLessThan = , LessThan = , StrictGreaterThan = , GreaterThan = PREC$Relational,
         And = PREC$And, Or = PREC$Or, Xor = PREC$Xor, Not = PREC$Not,
         Implies = PREC$Xor, Equivalent = PREC$Xor,
         Equality = , Unequality = PREC$Mul,
         # PRECEDENCE_VALUES gives Mod the precedence of Mul
         Function = if (identical(e$name, "Mod")) PREC$Mul else PREC$Atom,
         UnevaluatedExpr = .prec(e$args[[1L]]) - 0.5,
         # PRECEDENCE_VALUES and the class MRO for matrix expressions
         MatAdd = PREC$Add, MatPow = , Inverse = , HadamardPower = PREC$Pow,
         HadamardProduct = , KroneckerProduct = PREC$Mul,
         MatMul = if (is_number_type(e$args[[1L]]) && could_extract_minus_sign(e$args[[1L]])) PREC$Add else PREC$Mul,
         PREC$Atom)
}

.paren <- function(e, level, strict = FALSE) {
  s <- .str(e)
  p <- .prec(e)
  if (p < level || (!strict && p <= level)) paste0("(", s, ")") else s
}

.fmt_float <- function(v) {
  if (is.nan(v)) return("nan")
  if (is.infinite(v)) return(if (v > 0) "oo" else "-oo")
  s <- format(v, digits = 15)
  if (!grepl("[.eE]", s)) s <- paste0(s, ".0")
  s
}

# Printing order of Add terms: lexicographic in the generators, highest
# degree first, number-like terms last.
.add_print_order <- function(args, rev = FALSE) {
  n <- length(args)
  if (n < 2L) return(args)
  # SymPy as_ordered_terms special case: a positive number plus a negative
  # number times one factor prints as 5 - 2*x
  if (n == 2L && !rev) {
    isnum <- vapply(args, function(a) is_number_type(a) || (a$type == "Constant" && a$name %in% c("pi", "E", "EulerGamma", "Catalan", "GoldenRatio")), TRUE)
    if (sum(isnum) == 1L) {
      a0 <- args[[which(isnum)]]
      a1 <- args[[which(!isnum)]]
      if (a1$type == "Mul" && length(a1$args) == 2L && is_number_type(a1$args[[1L]]) &&
          num_sign(a1$args[[1L]]) < 0 && isTRUE(.facts(a0)[["positive"]]))
        return(list(a0, a1))
    }
  }
  info <- lapply(args, .term_monomial)
  gens <- unique(unlist(lapply(info, function(i) names(i$mon))))
  if (length(gens) > 1L) {
    # generators in SymPy default_sort_key order
    gsym <- lapply(gens, function(g) {
      for (i in info) if (!is.null(i$gsyms[[g]])) return(i$gsyms[[g]])
    })
    gens <- .key_sort(lapply(gsym, .sort_key), gens)
  }
  mat <- matrix(0, nrow = n, ncol = max(1L, length(gens)))
  for (i in seq_len(n)) {
    m <- info[[i]]$mon
    if (length(m)) mat[i, match(names(m), gens)] <- m
  }
  cols <- lapply(seq_len(ncol(mat)), function(j) -mat[, j])
  coefstr <- vapply(info, function(i) i$cstr, "")
  # within a monomial, SymPy sorts by the complex value of the numeric
  # part: real before imaginary, then by imaginary and real parts
  cv <- vapply(args, function(a) {
    sp <- as_coeff_mul_split(a)
    facs <- if (sp[[2L]]$type == "Mul") sp[[2L]]$args else list(sp[[2L]])
    num <- c(list(sp[[1L]]), Filter(.expr_is_number, facs))
    v <- tryCatch(as.complex(.neval(.mul(num), list())), error = function(e) NA_complex_)
    if (length(v) != 1L) NA_complex_ else v
  }, 0i)
  re <- Re(cv); im <- Im(cv)
  ord <- do.call(order, c(cols, list(.nc_rank(lapply(info, `[[`, "nc")), im != 0, im, re, coefstr), list(method = "radix", decreasing = rev)))
  args[ord]
}

# Dense ranks of the noncommutative parts, compared as tuples of sort keys.
.nc_rank <- function(nc) {
  n <- length(nc)
  if (!any(lengths(nc) > 0L)) return(integer(n))
  keys <- lapply(nc, function(l) lapply(l, .sort_key))
  cmp <- function(a, b) {
    for (i in seq_len(min(length(a), length(b)))) { c <- .key_cmp(a[[i]], b[[i]]); if (c != 0) return(c) }
    length(a) - length(b)
  }
  o <- integer()
  for (i in seq_len(n)) {
    j <- length(o)
    while (j > 0L && cmp(keys[[o[j]]], keys[[i]]) > 0) j <- j - 1L
    o <- append(o, i, after = j)
  }
  rk <- integer(n); r <- 0L
  for (k in seq_along(o)) {
    if (k == 1L || cmp(keys[[o[k - 1L]]], keys[[o[k]]]) != 0) r <- r + 1L
    rk[o[k]] <- r
  }
  rk
}

# Generator and integer degree of a factor (sympy/core/exprtools.py,
# decompose_power): x^(3/2) is sqrt(x) cubed, exp(-2*y) is exp(y)^-2.
.decompose_power <- function(f) {
  if (f$type == "Pow") {
    b <- f$args[[1L]]; e <- f$args[[2L]]
  } else if (f$type == "Function" && f$name == "exp") {
    b <- NULL; e <- f$args[[1L]]
  } else {
    return(list(f, 1))
  }
  mk <- function(t) if (is.null(b)) fcall("exp", list(t)) else .pow(b, t)
  if (is_number_type(e)) {
    if (e$type == "Float") return(list(f, 1))
    q <- gmp::denominator(e$value)
    # exp(k).as_base_exp() is (E, k)
    base <- if (q == 1) (if (is.null(b)) S_E else b) else mk(make_rational(gmp::as.bigq(1, q)))
    return(list(base, as.numeric(gmp::numerator(e$value))))
  }
  sp <- as_coeff_mul_split(e)
  c1 <- sp[[1L]]
  if (!(c1$type %in% c("Integer", "Rational")) || is_one_num(c1)) return(list(f, 1))
  q <- gmp::denominator(c1$value)
  list(mk(.div(sp[[2L]], make_rational(gmp::as.bigq(q)))), as.numeric(gmp::numerator(c1$value)))
}

.term_monomial <- function(t) {
  sp <- as_coeff_mul_split(t)
  rest <- sp[[2L]]
  facs <- if (rest$type == "Mul") rest$args else if (is_one_num(rest)) list() else list(rest)
  mon <- numeric()
  gsyms <- list()
  numlike <- TRUE
  cparts <- character()
  nc <- list()
  for (f in facs) {
    # Expr.is_number: no free symbols and no undefined function
    if (.expr_is_number(f)) {
      cparts <- c(cparts, .str(f))
      next
    }
    numlike <- FALSE
    # noncommutative factors form the ordered ncpart of SymPy's term key
    if (!.is_comm(f)) { nc[[length(nc) + 1L]] <- f; next }
    d <- .decompose_power(f)
    b <- d[[1L]]
    nm <- .str(b)
    gsyms[[nm]] <- b
    mon[nm] <- (if (is.na(mon[nm])) 0 else mon[nm]) + d[[2L]]
  }
  list(mon = mon, numlike = numlike, gsyms = gsyms, nc = nc,
       cstr = paste(cparts, collapse = "*"))
}

# Printing order of Mul factors: SymPy's as_ordered_factors(), which sorts
# the commutative factors by default_sort_key and keeps noncommutative ones
# last in their order.
.mul_print_order <- function(fs) {
  if (length(fs) < 2L) return(fs)
  comm <- vapply(fs, .is_comm, TRUE)
  if (!all(comm)) return(c(.mul_print_order(fs[comm]), fs[!comm]))
  .key_sort(lapply(fs, .sort_key), fs)
}

.str_add <- function(e) {
  is_o <- vapply(e$args, function(t) t$type == "Order", TRUE)
  if (any(is_o)) {
    # Expr.as_ordered_terms with Order terms: both groups sorted in reverse
    terms <- c(.add_print_order(e$args[!is_o], rev = TRUE), .add_print_order(e$args[is_o], rev = TRUE))
  } else if (any(vapply(e$args, .me_has_nc_factor, NA))) {
    terms <- .me_ordered_terms(e$args)
  } else {
    terms <- .add_print_order(e$args)
  }
  out <- character(length(terms))
  for (i in seq_along(terms)) {
    t <- terms[[i]]
    if (i == 1L) {
      out[i] <- .str_add_term(t)
    } else if (could_extract_minus_sign(t)) {
      nt <- if (t$type == "Mul") .neg_mul_noflat(t) else .neg(t)
      # -1/y^2 prints as a fraction, as SymPy's Mul printer does
      out[i] <- if (nt$type == "Pow" && .is_neg_exp(nt$args[[2L]]) && .is_comm(nt))
        paste0(" - 1/", .paren(.pow(nt$args[[1L]], .neg(nt$args[[2L]])), PREC$Mul, strict = TRUE))
      else if (nt$type == "Add") paste0(" - (", .str(nt), ")")
      # the negative term is printed at the precedence of Add, as SymPy does
      else if (nt$type == "Mul") paste0(" - ", .str_mul(nt, level = PREC$Add))
      else paste0(" - ", .str_add_term(nt))
    } else {
      out[i] <- paste0(" + ", .str_add_term(t))
    }
  }
  paste(out, collapse = "")
}

.str_add_term <- function(t) {
  # an Add term (a MatAdd inside Add) is parenthesized, as in StrPrinter._print_Add
  if (.prec(t) < PREC$Add || t$type == "MatAdd") paste0("(", .str(t), ")") else .str(t)
}

# -e for a Mul with a negative leading number, without redistributing an
# unevaluated coefficient such as -2*(x + 2) over its Add.
.neg_mul_noflat <- function(e) {
  c1 <- e$args[[1L]]
  if (!is_number_type(c1) || num_sign(c1) >= 0) return(.neg(e))
  rest <- e$args[-1L]
  if (is_one_num(.neg(c1))) return(if (length(rest) == 1L) rest[[1L]] else new_sym("Mul", rest))
  new_sym("Mul", c(list(.neg(c1)), rest))
}

.is_neg_exp <- function(ex) {
  if (is_number_type(ex)) return(num_sign(ex) < 0)
  could_extract_minus_sign(ex) && ex$type == "Mul"
}

# StrPrinter._print_Mul for an unevaluated Mul (a leading 1, or a number or
# an integer power after the first factor): the factors in their stored
# order, the negative powers as a denominator. NULL for any other Mul.
.str_mul_uneval <- function(e) {
  args <- e$args
  numlike <- function(a) is_number_type(a) || is_inf_const(a) || is_nan(a)
  if (!((args[[1L]]$type == "Integer" && is_one_num(args[[1L]])) ||
        any(vapply(args[-1L], function(a) numlike(a) ||
                     (a$type == "Pow" && all(vapply(a$args, function(z) z$type == "Integer", TRUE))), TRUE))))
    return(NULL)
  exp_neg <- function(a) {
    ex <- a$args[[2L]]
    c1 <- if (is_number_type(ex)) ex else if (ex$type == "Mul" && is_number_type(ex$args[[1L]])) ex$args[[1L]] else S_One
    num_sign(c1) < 0
  }
  isd <- vapply(args, function(a) a$type == "Pow" && exp_neg(a), TRUE)
  d <- lapply(args[isd], function(di) {
    ex <- di$args[[2L]]
    e2 <- if (is_number_type(ex)) .neg(ex) else {
      dargs <- ex$args
      dargs[[1L]] <- .neg(dargs[[1L]])
      if (length(dargs) == 1L) dargs[[1L]] else new_sym("Mul", dargs)
    }
    if (is_one_num(e2)) di$args[[1L]] else new_sym("Pow", list(di$args[[1L]], e2))
  })
  n <- args[!isd]
  pre <- character()
  if (length(n) && n[[1L]]$type != "Add" && could_extract_minus_sign(n[[1L]])) {
    pre <- .str(n[[1L]])
    n <- n[-1L]
  }
  prec <- .prec(e)
  nf <- c(pre, vapply(n, function(a) .paren(a, prec), ""))
  if (!length(nf)) nf <- "1"
  pre <- character()
  if (length(d) > 1L && could_extract_minus_sign(d[[1L]])) {
    pre <- .str(d[[1L]])
    d <- d[-1L]
  }
  df <- c(pre, vapply(d, function(a) .paren(a, prec), ""))
  ns <- paste(nf, collapse = "*")
  if (length(df) > 1L) return(paste0(ns, "/(", paste(df, collapse = "*"), ")"))
  if (length(df)) return(paste0(ns, "/", df))
  ns
}

# level: the precedence the factors are parenthesized against, that of the
# printed product (a negative product has the precedence of Add)
.str_mul <- function(e, pow_paren = TRUE, level = PREC$Mul) {
  s <- .str_mul_uneval(e)
  if (!is.null(s)) return(s)
  if (could_extract_minus_sign(e)) {
    ne <- .neg_mul_noflat(e)
    # -x^(-2*a) prints as -1/x^(2*a), as SymPy's Mul printer does
    if (ne$type == "Pow" && .is_neg_exp(ne$args[[2L]]) && .is_comm(ne))
      return(paste0("-1/", .paren(.pow(ne$args[[1L]], .neg(ne$args[[2L]])), PREC$Mul, strict = TRUE)))
    s <- if (ne$type == "Mul") .str_mul(ne, pow_paren = FALSE, level = PREC$Add) else .str(ne)
    if (.prec(ne) < PREC$Mul) s <- paste0("(", s, ")")
    return(paste0("-", s))
  }
  fs <- e$args
  num <- list()
  den <- list()
  coef_num <- NULL
  coef_den <- NULL
  lead <- NULL
  paren_keys <- character()
  if (fs[[1L]]$type %in% c("Integer", "Float")) {
    lead <- .str(fs[[1L]])
    fs <- fs[-1L]
  }
  for (f in fs) {
    if (f$type == "Rational") {
      n <- gmp::numerator(f$value)
      d <- gmp::denominator(f$value)
      if (n != 1) coef_num <- as.character(n)
      coef_den <- as.character(d)
    } else if (f$type == "Pow" && .is_neg_exp(f$args[[2L]]) && .is_comm(f)) {
      den[[length(den) + 1L]] <- .pow(f$args[[1L]], .neg(f$args[[2L]]))
      # an unevaluated 1/(1/y) keeps its parentheses (issue 21537)
      if (pow_paren && is_negone_num(f$args[[2L]]) && f$args[[1L]]$type %in% c("Mul", "Pow"))
        paren_keys <- c(paren_keys, f$args[[1L]]$key)
    } else {
      num[[length(num) + 1L]] <- f
    }
  }
  num <- .mul_print_order(num)
  den <- .mul_print_order(den)
  nstr <- vapply(num, function(f) .paren(f, level, strict = FALSE), "")
  if (!is.null(coef_num)) nstr <- c(coef_num, nstr)
  if (!is.null(lead)) nstr <- c(lead, nstr)
  if (!length(nstr)) nstr <- "1"
  s <- paste(nstr, collapse = "*")
  dstr <- vapply(den, function(f) { p <- .paren(f, level, strict = FALSE)
    if (f$key %in% paren_keys) paste0("(", p, ")") else p }, "")
  if (!is.null(coef_den)) dstr <- c(coef_den, dstr)
  if (!length(dstr)) return(s)
  if (length(dstr) == 1L) return(paste0(s, "/", dstr))
  paste0(s, "/(", paste(dstr, collapse = "*"), ")")
}

.str_pow <- function(e) {
  b <- e$args[[1L]]
  x <- e$args[[2L]]
  if (x$type == "Rational" && x$value == gmp::as.bigq(1, 2))
    return(paste0("sqrt(", .str(b), ")"))
  nc <- !.is_comm(b)
  if (!nc && x$type == "Rational" && x$value == gmp::as.bigq(-1, 2))
    return(paste0("1/sqrt(", .str(b), ")"))
  # noncommutative inverses keep their place: A^(-1)*B, not B/A
  if (!nc && is_negone_num(x)) return(paste0("1/", .paren(b, PREC$Mul, strict = FALSE)))
  bs <- .paren(b, PREC$Pow, strict = FALSE)
  if (is_number_type(b) && b$type == "Rational") bs <- paste0("(", .str(b), ")")
  xs <- .str(x)
  py <- isTRUE(.str_mode$python)
  # Python's Rational(p, q) needs no parentheses as an exponent
  if (!(py && x$type == "Rational") && (.prec(x) < PREC$Atom || x$type == "Rational")) xs <- paste0("(", xs, ")")
  if (py && is_number_type(b) && b$type == "Rational") bs <- .str(b)
  paste0(bs, if (py) "**" else "^", xs)
}

.REL_OPS <- c(Equality = "==", Unequality = "!=", StrictLessThan = "<",
              LessThan = "<=", StrictGreaterThan = ">", GreaterThan = ">=")


.str_mode <- new.env()

.str <- function(e) {
  if (!is.null(.str_mode$hook) && !is.null(r <- .str_mode$hook(e))) return(r)
  t <- e$type
  switch(t,
         Integer = as.character(gmp::numerator(e$value)),
         Rational = if (isTRUE(.str_mode$python)) paste0("Rational(", gmp::numerator(e$value), ", ", gmp::denominator(e$value), ")")
                    else as.character(e$value),
         Float = .fmt_float_sym(e, strip = TRUE),
         Symbol = e$name,
         Dummy = paste0("_", e$name),
         Wild = paste0(e$name, "_"),
         Constant = e$name,
         BooleanAtom = if (e$name == "true") "True" else "False",
         Add = .str_add(e),
         Mul = .str_mul(e),
         Pow = .str_pow(e),
         Function = {
           d <- .FUNS[[e$name]]
           if (!is.null(d) && !is.null(d$str)) d$str(e)
           else paste0(e$name, "(", paste(vapply(e$args, .str, ""), collapse = ", "), ")")
         },
         Equality = paste0("Eq(", .str(e$args[[1L]]), ", ", .str(e$args[[2L]]), ")"),
         Unequality = paste0("Ne(", .str(e$args[[1L]]), ", ", .str(e$args[[2L]]), ")"),
         StrictLessThan = , LessThan = , StrictGreaterThan = , GreaterThan = {
           paste(.paren(e$args[[1L]], PREC$Relational), .REL_OPS[[t]],
                 .paren(e$args[[2L]], PREC$Relational))
         },
         {
           f <- .str_printers[[t]]
           if (!is.null(f)) f(e)
           else paste0(t, "(", paste(vapply(e$args, .str, ""), collapse = ", "), ")")
         })
}

#' Convert an expression to a string
#'
#' `sstr()` gives the SymPy style string form with R operators, so the result
#' can be parsed again with [sympify()]. `srepr()` shows the expression tree.
#'
#' @param e A `sym` expression.
#' @return A character scalar.
#' @examples
#' x <- Symbol("x")
#' sstr(x^2 + 2*x + 1)
#' srepr(x^2 + 1)
#' @export
sstr <- function(e) {
  if (inherits(e, "physics_obj")) return(format(e))
  if (inherits(e, "TensExpr")) return(.tn_str(e))
  # free group words, fp groups, coset tables and polyhedra print as their str()
  if (inherits(e, c("FreeGroupElement", "FreeGroup", "FpGroup", "FpSubgroup", "CosetTable", "Polyhedron")))
    return(fp_str(e))
  # permutations, cycles, groups and cosets print as their own str()
  if (inherits(e, c("Permutation", "Cycle", "PermutationGroup", "Coset", "SymmetricPermutationGroup",
                    "TransitiveSubgroupName", "galois_group_result")))
    return(format(e))
  # MatrixBase.__str__: Matrix(<rows as a Python list>), or the shape when empty
  if (inherits(e, "symmatrix")) {
    if (e$nrow == 0L || e$ncol == 0L) return(sprintf("Matrix(%d, %d, [])", e$nrow, e$ncol))
    # str(self.tolist()): each entry prints as its own repr, so Floats keep their zeros
    rows <- vapply(.mrows(e), function(r) paste0("[", paste(vapply(r, .str_top, ""), collapse = ", "), "]"), "")
    return(paste0("Matrix([", paste(rows, collapse = ", "), "])"))
  }
  if (inherits(e, "SymArray")) return(.ae_str_ndarray(e$data, e$shape))
  .str_top(as_sym(e))
}

# SymPy prints a Float at top level with all its digits, and strips trailing
# zeros when it is part of a larger expression; Prefix.__str__ is the
# abbreviation, while inside an expression a prefix prints as its name
.str_top <- function(e) if (e$type == "Float") .fmt_float_sym(e, strip = FALSE) else
  if (e$type == "Prefix") .str(e$args[[2L]]) else .str(e)

#' @rdname sstr
#' @export
srepr <- function(e) {
  e <- as_sym(e)
  t <- e$type
  # codegen nodes print through their own repr (R/codegen-ast.R)
  f <- .srepr_printers[[t]]
  if (!is.null(f)) return(f(e))
  switch(t,
         Integer = paste0("Integer(", .str(e), ")"),
         Rational = paste0("Rational(", gmp::numerator(e$value), ", ",
                           gmp::denominator(e$value), ")"),
         Float = paste0("Float('", .fmt_float_sym(e, strip = TRUE), "', precision=", .float_prec(e), ")"),
         Symbol = if (isTRUE(e$extra$idx)) {
           r <- e$extra$range
           paste0("Idx(Symbol('", e$name, "', integer=True)",
                  if (!is.null(r)) paste0(", Tuple(", srepr(r[[1L]]), ", ", srepr(r[[2L]]), ")") else "", ")")
         } else {
           # the assumptions as given, sorted, as SymPy's srepr
           given <- if (!is.null(e$extra$given)) e$extra$given else { f <- e$extra$facts; f[!is.na(f) & names(f) != "commutative"] }
           if (length(given)) given <- given[order(names(given), method = "radix")]
           if (length(given)) paste0("Symbol('", e$name, "', ",
                                     paste0(names(given), "=", ifelse(given, "True", "False"), collapse = ", "), ")")
           else paste0("Symbol('", e$name, "')")
         },
         Dummy = paste0("Dummy('", e$name, "')"),
         Wild = paste0(if (isTRUE(e$extra$wildfunction)) "WildFunction('" else "Wild('", e$name, "')"),
         Constant = e$name,
         BooleanAtom = .str(e),
         Predicate = , UndefinedPredicate = .str(e),
         Function = paste0(e$name, "(", paste(vapply(e$args, srepr, ""), collapse = ", "), ")"),
         paste0(t, "(", paste(vapply(e$args, srepr, ""), collapse = ", "), ")"))
}

#' @export
format.sym <- function(x, ...) .str_top(x)

#' @export
as.character.sym <- function(x, ...) .str_top(x)

#' @export
print.sym <- function(x, ...) {
  # init_printing() chooses the console form
  mode <- getOption("symcompR.print", "str")
  if (identical(mode, "pretty")) { pprint(x); return(invisible(x)) }
  if (identical(mode, "latex")) { cat(latex(x), "\n", sep = ""); return(invisible(x)) }
  cat(.str_top(x), "\n", sep = "")
  invisible(x)
}

#' @export
toString.sym <- function(x, ...) .str_top(x)

# Expr.is_number: no free symbols, no undefined function, and no argument
# that is a Basic but not an Expr (a coordinate system, for one)
.NON_EXPR_TYPES <- c("CoordSystem", "Patch", "Manifold", "Point", "Object", "Str",
                     # Quantity.is_number and Dimension.is_number are False
                     "Quantity", "Dimension",
                     # MatrixSymbol.free_symbols is the symbol itself
                     "MatrixSymbol")
.expr_is_number <- function(f) !length(free_symbols(f)) && !.has_undef(f) && !.has_type(f, .NON_EXPR_TYPES)
