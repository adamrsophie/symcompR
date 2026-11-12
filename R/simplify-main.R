# simplify() after sympy/simplify/simplify.py, with collect_abs,
# inversecombine, ratsimp, sum_simplify and product_simplify of the same
# module.

.BOOL_TYPES <- c("And", "Or", "Not", "Xor", "Implies", "Equivalent", "ITE", "BooleanAtom")
.NONEXPR_TYPES <- c(.BOOL_TYPES, "Tuple", "ExprCondPair", "FiniteSet", "Interval", "Union", "Intersection",
                    "Complement", "EmptySet", "ProductSet", "ImageSet", "ConditionSet")
.is_expr_node <- function(e) !(is_relational(e) || e$type %in% .NONEXPR_TYPES || inherits(e, "symmatrix"))

.BESSEL_FUNS <- c("besselj", "bessely", "besseli", "besselk", "hankel1", "hankel2", "jn", "yn", "hn1", "hn2")
.COMB_FUNS <- c("factorial", "factorial2", "subfactorial", "RisingFactorial", "FallingFactorial", "rf", "ff",
                "binomial", "MultiFactorial")
.INVERSE_FUNS <- c(exp = "log", log = "exp", sin = "asin", cos = "acos", tan = "atan", cot = "acot",
                   sec = "asec", csc = "acsc", asin = "sin", acos = "cos", atan = "tan", acot = "cot",
                   asec = "sec", acsc = "csc", sinh = "asinh", cosh = "acosh", tanh = "atanh",
                   coth = "acoth", sech = "asech", csch = "acsch", asinh = "sinh", acosh = "cosh",
                   atanh = "tanh", acoth = "coth", asech = "sech", acsch = "csch", erf = "erfinv",
                   erfinv = "erf", erfc = "erfcinv", erfcinv = "erfc")

#' Simplify an expression
#'
#' Runs SymPy's `simplify()` heuristics: rational simplification, `powsimp`,
#' cancellation and `together`, `factor_terms`, then the specialized
#' simplifications (Piecewise, hypergeometric, Bessel, trigonometric,
#' logarithmic, combinatorial, sums, integrals and products) and keeps the
#' shortest result under `measure`, returning the input when the result is
#' more than `ratio` times as complex.
#'
#' @param e An expression, matrix or list.
#' @param ratio Accept a result only if it is not more than `ratio` times as
#'   complex as the input (`Inf` keeps the shortest).
#' @param measure The complexity measure, `count_ops()` by default.
#' @param rational `FALSE` keeps Floats, `TRUE` converts them to rationals,
#'   `NULL` converts them and back at the end.
#' @param inverse Also combine functions with their inverses
#'   (`inversecombine()`).
#' @param doit Evaluate unevaluated objects at the end.
#' @param ... Ignored.
#' @return A `sym` expression.
#' @examples
#' x <- Symbol("x")
#' simplify(sin(x)^2 + cos(x)^2)
#' simplify((x^3 + x^2 - x - 1) / (x^2 + 2*x + 1))
#' @export
simplify <- function(e, ratio = 1.7, measure = count_ops, rational = FALSE, inverse = FALSE, doit = TRUE, ...) {
  if (inherits(e, "symmatrix")) return(.mat_map(e, function(z) simplify(z, ratio, measure, rational, inverse, doit)))
  if (is.list(e) && !inherits(e, "sym")) return(lapply(e, simplify, ratio = ratio, measure = measure, rational = rational, inverse = inverse, doit = doit))
  if (.is_bd(e)) return(.bd_simplify(e, ratio = ratio, measure = measure, rational = rational, inverse = inverse, doit = doit))
  e <- as_sym(e)
  if (isTRUE(rational)) e <- nsimplify(e, rational = TRUE)
  .simplify(e, list(ratio = ratio, measure = measure, rational = rational, inverse = inverse, doit = doit))
}

.simp_measure <- function(z, kw) as.numeric(kw$measure(z))

.simplify <- function(expr, kw) {
  measure <- function(z) .simp_measure(z, kw)
  ratio <- kw$ratio
  sub <- function(z) .simplify(z, kw)
  shorter <- function(...) {
    ch <- list(...)
    if (all(vapply(ch, function(z) identical(z$key, ch[[1L]]$key), TRUE))) return(ch[[1L]])
    ch[[which.min(vapply(ch, measure, 0))]]
  }
  done <- function(z) {
    rv <- if (isTRUE(kw$doit)) doit(z) else z
    shorter(rv, collect_abs(rv))
  }
  if (.is_expr_node(expr) && isTRUE(ask_fact(expr, "zero"))) return(if (is_number_type(expr)) expr else S_Zero)
  # _eval_simplify of the classes that define one
  r <- .eval_simplify(expr, kw)
  if (!is.null(r)) return(r)
  original_expr <- expr <- collect_abs(signsimp(expr))
  if (!length(expr$args)) return(expr)
  if (isTRUE(kw$inverse) && .has_type(expr, "Function")) {
    expr <- inversecombine(expr)
    if (!length(expr$args)) return(expr)
  }
  handled <- function(z) z$type %in% c("Add", "Mul", "Pow") || (z$type == "Function" && z$name %in% c("exp", "exp_polar"))
  expr <- .map_bottom_up(expr, function(z) {
    if (length(z$args) && .is_expr_node(z) && !handled(z)) sym_rebuild(z, lapply(z$args, sub)) else z
  })
  if (!handled(expr)) return(done(expr))
  if (!.is_comm(expr)) expr <- nc_simplify(expr)
  floats <- FALSE
  if (!isFALSE(kw$rational) && .has_type(expr, "Float")) {
    floats <- TRUE
    expr <- nsimplify(expr, rational = TRUE)
  }
  expr <- .map_bottom_up(expr, function(w) if (length(w$args) && .is_expr_node(w)) .expr_normal(w) else w)
  safe_ratsimp <- function(arg) {
    if (!.is_comm(arg) || !length(arg$args)) return(arg)
    nd <- .as_numer_denom(arg)
    if (is_number_type(nd[[2L]]) || !.has_type(nd[[1L]], "Add")) return(arg)
    ratsimp(arg)
  }
  expr <- .map_bottom_up(expr, function(z) if (z$type == "Function" && z$name == "exp") fcall("exp", list(safe_ratsimp(z$args[[1L]]))) else z)
  expr <- .map_bottom_up(expr, function(z) if (z$type == "Pow") .pow(safe_ratsimp(z$args[[1L]]), safe_ratsimp(z$args[[2L]])) else z)
  expr <- .mul(as_content_primitive(powsimp(expr)))
  .e <- cancel(expr)
  expr1 <- shorter(.e, cancel(.mexpand(.e)))
  expr2 <- shorter(together(expr, deep = TRUE), together(expr1, deep = TRUE))
  expr <- if (is.infinite(ratio)) expr2 else shorter(expr2, expr1, expr)
  expr <- factor_terms(expr, sign = FALSE)
  if (.has_function(expr, "sign")) expr <- rewrite(expr, "Abs")
  if (.has_type(expr, "Piecewise")) {
    expr <- piecewise_fold(expr)
    expr <- done(expr)
    if (.has_type(expr, "Piecewise")) {
      expr <- piecewise_fold(expr)
      if (.has_function(expr, "KroneckerDelta")) expr <- kroneckersimp(expr)
      if (.has_type(expr, "Piecewise")) {
        kw2 <- kw; kw2$doit <- FALSE
        expr <- piecewise_simplify(expr, deep = TRUE, kw = kw2)
        if (.has_type(expr, "Piecewise")) return(shorter(expr, factor_terms(expr)))
      }
    }
  }
  if (.has_function(expr, c("hyper", "meijerg"))) expr <- hyperexpand(expr)
  if (.has_function(expr, "KroneckerDelta")) expr <- kroneckersimp(expr)
  if (.has_function(expr, .BESSEL_FUNS)) expr <- besselsimp(expr)
  if (.has_function(expr, .TRIGS_ALL)) expr <- trigsimp(expr, deep = TRUE)
  if (.has_function(expr, "log")) expr <- shorter(expand_log(expr, deep = TRUE), logcombine(expr))
  if (.has_function(expr, c(.COMB_FUNS, "gamma"))) expr <- combsimp(expr)
  if (.has_type(expr, "Sum")) expr <- sum_simplify(expr, kw)
  if (.has_type(expr, "Integral")) {
    ints <- unique(atoms(expr, "Integral"))
    expr <- xreplace(expr, lapply(ints, function(i) list(i, factor_terms(i))))
  }
  if (.has_type(expr, "Product")) expr <- product_simplify(expr, kw)
  if (.has_type(expr, "Quantity")) expr <- quantity_simplify(expr)
  short <- shorter(powsimp(expr, combine = "exp", deep = TRUE), powsimp(expr), expr)
  short <- shorter(short, cancel(short))
  short <- shorter(short, factor_terms(short), expand_power_exp(expand_mul(short)))
  if (.has_function(short, c(.TRIGS_ALL, "exp"))) short <- exptrigsimp(short)
  # hollow 2-arg Muls 2*(x + 1) are multiplied out
  expr <- .xreplace_td(short, function(z) z$type == "Mul" && length(z$args) == 2L && is_number_type(z$args[[1L]]) &&
                         z$args[[2L]]$type == "Add" && .is_comm(z), function(z) .mul(z$args))
  nd <- .as_numer_denom(expr)
  if (nd[[2L]]$type == "Add") {
    fr <- fraction(radsimp(.inv(nd[[2L]]), symbolic = FALSE, max_terms = 1L))
    if (!(is_one_num(fr[[1L]]) && fr[[1L]]$type == "Integer")) expr <- .div(expand(.mul(list(nd[[1L]], fr[[1L]]))), fr[[2L]])
  }
  if (could_extract_minus_sign(expr)) {
    fr <- fraction(expr)
    if (!is_zero_num(fr[[2L]])) expr <- signsimp(.div(.neg(fr[[1L]]), .neg(fr[[2L]])))
  }
  if (measure(expr) > ratio * measure(original_expr)) expr <- original_expr
  if (floats && is.null(kw$rational)) expr <- nfloat(expr, exponent = FALSE)
  done(expr)
}

# The _eval_simplify methods; NULL when the class has none
.eval_simplify <- function(e, kw) {
  measure <- function(z) .simp_measure(z, kw)
  sub <- function(z) .simplify(z, kw)
  t <- e$type
  if (!length(e$args)) return(e)
  if (t == "Function") {
    nm <- e$name
    if (nm == "hyper") return(hyperexpand(e))
    if (nm %in% .BESSEL_FUNS) return(besselsimp(e))
    if (nm %in% .COMB_FUNS) {
      r <- combsimp(e)
      return(if (measure(r) <= kw$ratio * measure(e)) r else e)
    }
    if (nm == "sign") return(fcall("sign", list(factor_terms(e$args[[1L]]))))
    if (nm == "log") {
      r <- fcall("log", list(sub(e$args[[1L]])))
      if (isTRUE(kw$inverse)) r <- inversecombine(r)
      r <- expand_log(r, deep = TRUE)
      return(if (measure(r) < measure(e)) r else e)
    }
    return(NULL)
  }
  if (t == "Piecewise") return(piecewise_simplify(e, kw = kw))
  if (t == "ExprCondPair") return(sym_rebuild(e, lapply(e$args, sub)))
  if (t == "Integral") {
    r <- factor_terms(e)
    if (r$type == "Integral") return(sym_rebuild(r, lapply(r$args, sub)))
    return(sub(r))
  }
  if (t == "Sum") return(.sum_eval_simplify(e, kw))
  if (t == "Product") {
    r <- product_simplify(e, kw)
    return(if (isTRUE(kw$doit)) doit(r) else r)
  }
  if (t == "FiniteSet") return(do.call(FiniteSet, lapply(e$args, sub)))
  if (t == "Tuple") return(sym_rebuild(e, lapply(e$args, sub)))
  if (is_relational(e)) return(.relational_simplify(e, kw))
  if (t %in% .BOOL_TYPES) return(.boolean_simplify(e, kw))
  NULL
}

#' Simplification helpers
#'
#' `collect_abs()` joins products of absolute values, `|x|*|y| = |x*y|`.
#' `inversecombine()` simplifies functions of their inverses such as
#' `log(exp(x))` and `asin(sin(x))`. `ratsimp()` writes a rational function
#' as polynomial plus a proper part. `sum_simplify()` and
#' `product_simplify()` combine sums and products with the same limits, as
#' the SymPy functions of the same names.
#'
#' @param e,expr An expression.
#' @param kw Options of `simplify()` as a list (internal).
#' @return An expression.
#' @examples
#' x <- Symbol("x", real = TRUE); y <- Symbol("y", real = TRUE)
#' collect_abs(abs(x) * abs(y))
#' inversecombine(log(exp(x)))
#' @name simplify-helpers
NULL

#' @rdname simplify-helpers
#' @export
collect_abs <- function(expr) {
  expr <- as_sym(expr)
  ab <- function(m) {
    fs <- .margs(m)
    comm <- vapply(fs, .is_comm, TRUE)
    c <- fs[comm]; nc <- fs[!comm]
    a <- list(); o <- list()
    for (i in c) {
      if (i$type == "Function" && i$name == "Abs") a[[length(a) + 1L]] <- i$args[[1L]]
      else if (i$type == "Pow" && i$args[[1L]]$type == "Function" && i$args[[1L]]$name == "Abs" && isTRUE(ask_fact(i$args[[2L]], "real")))
        a[[length(a) + 1L]] <- .pow(i$args[[1L]]$args[[1L]], i$args[[2L]])
      else o[[length(o) + 1L]] <- i
    }
    if (length(a) < 2L && !any(vapply(a, function(i) i$type == "Pow" && isTRUE(ask_fact(i$args[[2L]], "negative")), TRUE))) return(m)
    absarg <- .mul(a)
    A <- fcall("Abs", list(absarg))
    if (!.has_function(A, "Abs")) return(.mul(c(list(A), o, nc)))
    if (!(A$type == "Function" && A$name == "Abs")) A <- fnode("Abs", list(absarg))
    args <- c(list(A), o)
    args <- args[order(vapply(args, function(z) !is_number_type(z), TRUE))]
    new_sym("Mul", c(args, nc))
  }
  expr <- .map_bottom_up(expr, function(z) if (z$type == "Mul") ab(z) else z)
  .map_bottom_up(expr, function(z) if (z$type == "Pow") ab(z) else z)
}

#' @rdname simplify-helpers
#' @export
inversecombine <- function(expr) {
  .bottom_up_args(as_sym(expr), function(rv) {
    if (rv$type == "Function" && rv$name == "log") {
      a <- rv$args[[1L]]
      if (a$type == "Function" && a$name == "exp") rv <- a$args[[1L]]
    } else if (rv$type == "Function" && rv$name %in% names(.INVERSE_FUNS)) {
      a <- rv$args[[1L]]
      if (length(rv$args) == 1L && length(a$args) == 1L && a$type == "Function" && a$name == .INVERSE_FUNS[[rv$name]])
        rv <- a$args[[1L]]
    }
    if (rv$type == "Function" && rv$name == "exp" && rv$args[[1L]]$type == "Function" && rv$args[[1L]]$name == "log")
      rv <- rv$args[[1L]]$args[[1L]]
    rv
  })
}

#' @rdname simplify-helpers
#' @export
ratsimp <- function(expr) {
  expr <- as_sym(expr)
  fg <- .as_numer_denom(cancel(expr)); f <- fg[[1L]]; g <- fg[[2L]]
  r <- tryCatch(reduced(f, list(g), field = TRUE, expand = FALSE), error = function(err) NULL)
  if (is.null(r)) return(.div(f, g))
  .add(list(.add(r$q), cancel(.div(r$r, g))))
}

#' @rdname simplify-helpers
#' @export
sum_simplify <- function(e, kw = .simp_kw()) {
  s <- as_sym(e)
  if (s$type != "Add") {
    adds <- Filter(function(a) .has_type(a, "Sum"), unique(atoms(s, "Add")))
    if (length(adds)) s <- xreplace(s, lapply(adds, function(a) list(a, sum_simplify(a, kw))))
  }
  s <- expand(s)
  if (s$type != "Add") return(s)
  s_t <- list(); o_t <- list()
  for (term in s$args) {
    fs <- .margs(term)
    isum <- vapply(fs, function(i) i$type == "Sum", TRUE)
    if (!any(isum)) { o_t[[length(o_t) + 1L]] <- term; next }
    s_t[[length(s_t) + 1L]] <- .mul(c(list(.mul(fs[!isum])), lapply(fs[isum], .sum_eval_simplify, kw = kw)))
  }
  .add(c(list(.sum_combine(s_t, kw)), o_t))
}

.simp_kw <- function() list(ratio = 1.7, measure = count_ops, rational = FALSE, inverse = FALSE, doit = TRUE)

# Sum._eval_simplify
.sum_eval_simplify <- function(e, kw) {
  f <- e$args[[1L]]; limits <- e$args[-1L]
  if (!isFALSE(kw$deep)) f <- .simplify(f, kw)
  s_t <- list(); o_t <- list()
  for (term in .aargs(expand(f))) {
    if (.has_type(term, "Sum")) {
      out <- lapply(.margs(expand(term)), function(st) if (st$type == "Sum") .sum_eval_simplify(st, kw) else st)
      s_t[[length(s_t) + 1L]] <- .mul(out)
    } else o_t[[length(o_t) + 1L]] <- term
  }
  result <- .add(c(list(.sum_combine(s_t, kw)), o_t))
  factor_terms(sym_rebuild(e, c(list(result), limits)))
}

.sum_combine <- function(s_t, kw) {
  used <- rep(FALSE, length(s_t))
  for (method in 0:1) {
    for (i in seq_along(s_t)) {
      if (used[i]) next
      for (j in seq_along(s_t)) {
        if (used[j] || i == j) next
        temp <- .sum_add(s_t[[i]], s_t[[j]], method, kw)
        if (temp$type %in% c("Sum", "Mul")) { s_t[[i]] <- temp; used[j] <- TRUE }
      }
    }
  }
  .add(s_t[!used])
}

.lim3 <- function(l) l$args
.sum_add <- function(self, other, method, kw) {
  refactor <- function(v) {
    fs <- .margs(v)
    k <- which(vapply(fs, function(x) x$type == "Sum", TRUE))[1L]
    sv <- fs[[k]]
    sym_rebuild(sv, c(list(.mul(c(fs[-k], list(sv$args[[1L]])))), sv$args[-1L]))
  }
  rs <- if (self$type == "Mul") refactor(self) else self
  ro <- if (other$type == "Mul") refactor(other) else other
  if (identical(rs$type, ro$type) && rs$type == "Sum") {
    ls <- rs$args[-1L]; lo <- ro$args[-1L]
    if (method == 0L) {
      if (identical(vapply(ls, .k, ""), vapply(lo, .k, "")))
        return(factor_terms(sym_rebuild(rs, c(list(.add(list(rs$args[[1L]], ro$args[[1L]]))), ls))))
    } else if (is_zero_num(.simplify(.sub(rs$args[[1L]], ro$args[[1L]]), kw)) && length(ls) == 1L && length(lo) == 1L) {
      a <- .lim3(ls[[1L]]); b <- .lim3(lo[[1L]])
      if (identical(a[[1L]]$key, b[[1L]]$key)) {
        if (identical(b[[2L]]$key, .add(list(a[[3L]], S_One))$key))
          return(factor_terms(sym_rebuild(rs, list(rs$args[[1L]], Tuple(a[[1L]], a[[2L]], b[[3L]])))))
        if (identical(a[[2L]]$key, .add(list(b[[3L]], S_One))$key))
          return(factor_terms(sym_rebuild(rs, list(rs$args[[1L]], Tuple(a[[1L]], b[[2L]], a[[3L]])))))
      }
    }
  }
  .add(list(self, other))
}

#' @rdname simplify-helpers
#' @export
product_simplify <- function(e, kw = .simp_kw()) {
  s <- as_sym(e)
  terms <- .margs(s)
  p_t <- list(); o_t <- list()
  deep <- !isFALSE(kw$deep)
  for (term in terms) {
    if (term$type == "Product") {
      p_t[[length(p_t) + 1L]] <- if (deep) sym_rebuild(term, c(list(.simplify(term$args[[1L]], kw)), term$args[-1L])) else term
    } else o_t[[length(o_t) + 1L]] <- term
  }
  used <- rep(FALSE, length(p_t))
  for (method in 0:1) {
    for (i in seq_along(p_t)) {
      if (used[i]) next
      for (j in seq_along(p_t)) {
        if (used[j] || i == j) next
        tmp <- .product_mul(p_t[[i]], p_t[[j]], method, kw)
        if (tmp$type == "Product") { p_t[[i]] <- tmp; used[j] <- TRUE }
      }
    }
  }
  result <- .mul(o_t)
  for (i in seq_along(p_t)) if (!used[i]) result <- .mul(list(result, p_t[[i]]))
  result
}

.product_mul <- function(self, other, method, kw) {
  if (identical(self$type, other$type)) {
    ls <- self$args[-1L]; lo <- other$args[-1L]
    if (method == 0L) {
      if (identical(vapply(ls, .k, ""), vapply(lo, .k, "")))
        return(sym_rebuild(self, c(list(.mul(list(self$args[[1L]], other$args[[1L]]))), ls)))
    } else if (is_zero_num(.simplify(.sub(self$args[[1L]], other$args[[1L]]), kw)) && length(ls) == 1L && length(lo) == 1L) {
      a <- .lim3(ls[[1L]]); b <- .lim3(lo[[1L]])
      if (identical(a[[1L]]$key, b[[1L]]$key)) {
        if (identical(b[[2L]]$key, .add(list(a[[3L]], S_One))$key))
          return(sym_rebuild(self, list(self$args[[1L]], Tuple(a[[1L]], a[[2L]], b[[3L]]))))
        if (identical(a[[2L]]$key, .add(list(b[[3L]], S_One))$key))
          return(sym_rebuild(self, list(self$args[[1L]], Tuple(a[[1L]], b[[2L]], a[[3L]]))))
      }
    }
  }
  .mul(list(self, other))
}

#' Simplify a noncommutative expression
#'
#' `nc_simplify(expr)` simplifies a noncommutative product by collecting
#' repeated subterms into powers and inverses, `a*b*a*b*c` becoming
#' `(a*b)^2*c`, as SymPy's `nc_simplify()`. It works on products of
#' noncommutative symbols and on matrix expressions.
#'
#' @param expr An expression.
#' @param deep Also simplify the arguments of the products found.
#' @return An expression.
#' @examples
#' a <- Symbol("a", commutative = FALSE); b <- Symbol("b", commutative = FALSE)
#' nc_simplify(a * b * a * b * a * b)
#' @export
nc_simplify <- function(expr, deep = TRUE) {
  expr <- as_sym(expr)
  ismat <- expr$type %in% c("MatAdd", "MatMul", "MatPow", "MatrixSymbol", "Inverse", "Transpose")
  if (ismat) {
    expr <- doit(expr)
    tA <- "MatAdd"; tM <- "MatMul"; tP <- "MatPow"; tS <- "MatrixSymbol"
    fA <- function(xs) do.call(.madd, xs); fM <- function(xs) if (length(xs) == 1L) xs[[1L]] else do.call(.mm, xs)
    fP <- function(b, e) MatPow(b, e)
  } else {
    tA <- "Add"; tM <- "Mul"; tP <- "Pow"; tS <- "Symbol"
    fA <- .add; fM <- .mul; fP <- function(b, e) .pow(b, as_sym(e))
  }
  eqk <- function(a, b) identical(a$key, b$key)
  seqeq <- function(xs, ys) length(xs) == length(ys) && all(vapply(seq_along(xs), function(i) eqk(xs[[i]], ys[[i]]), TRUE))
  overlaps <- function(args) {
    n <- length(args)
    m <- list(lapply(args[-1L], function(a) if (eqk(a, args[[1L]])) c(1L, 0L) else 0L))
    if (n > 1L) for (i in seq_len(n - 1L)) {
      ov <- list()
      if (n - i - 1L > 0L) for (j in 0:(n - i - 2L)) {
        o <- integer()
        for (v in m[[i]][[j + 2L]]) {
          idx <- j + i + 1L + v
          if (idx < n && eqk(args[[i + 1L]], args[[idx + 1L]])) o <- c(o, v + 1L)
        }
        ov[[length(ov) + 1L]] <- c(o, 0L)
      }
      m[[i + 1L]] <- ov
    }
    m
  }
  reduce_inverses <- function(a0) {
    inv_tot <- 0; inverses <- list(); args <- list()
    for (arg in a0) {
      if (arg$type == tP && isTRUE(ask_fact(arg$args[[2L]], "extended_negative"))) {
        inverses <- c(list(fP(arg, -1L)), inverses); inv_tot <- inv_tot + 1
      } else {
        if (length(inverses) == 1L) args[[length(args) + 1L]] <- fP(inverses[[1L]], -1L)
        else if (length(inverses) > 1L) { args[[length(args) + 1L]] <- fP(fM(inverses), -1L); inv_tot <- inv_tot - (length(inverses) - 1) }
        inverses <- list()
        args[[length(args) + 1L]] <- arg
      }
    }
    if (length(inverses)) { args[[length(args) + 1L]] <- fP(fM(inverses), -1L); inv_tot <- inv_tot - (length(inverses) - 1) }
    list(inv_tot, args)
  }
  get_score <- function(s) {
    if (s$type == tP) return(get_score(s$args[[1L]]))
    if (s$type %in% c(tA, tM)) return(sum(vapply(s$args, get_score, 0)))
    1
  }
  compare <- function(s, alt) if (!eqk(s, alt) && get_score(alt) < get_score(s)) alt else s
  if (!(expr$type %in% c(tA, tM, tP)) || (!ismat && .is_comm(expr))) return(expr)
  args <- expr$args
  if (expr$type == tP) return(if (deep) doit(fP(nc_simplify(args[[1L]]), args[[2L]])) else expr)
  if (expr$type == tA) return(doit(fA(lapply(args, nc_simplify, deep = deep))))
  comm <- vapply(args, .is_comm, TRUE)
  if (!ismat) {
    com_coeff <- .mul(args[comm])
    if (!(is_one_num(com_coeff))) return(.mul(list(com_coeff, nc_simplify(.div(expr, com_coeff), deep = deep))))
    args <- args[!comm]
  }
  ri <- reduce_inverses(args); inv_tot <- ri[[1L]]; args <- ri[[2L]]
  invert <- FALSE
  if (inv_tot > length(args) / 2) { invert <- TRUE; args <- lapply(rev(args), function(a) fP(a, -1L)) }
  if (deep) args <- lapply(args, nc_simplify)
  m <- overlaps(args)
  n <- length(args)
  simps <- list(); post <- S_One; pre <- S_One; max_simp_coeff <- 0; simp <- NULL
  skey <- function(st) paste(vapply(st, .k, ""), collapse = "\u0001")
  sub_args <- function(a, b) if (b >= a) args[a:b] else list()
  rep <- NULL
  if (n > 1L) for (i in 1:(n - 1L)) {
    simp_coeff <- 0; l <- 0L; p <- 0
    if (i < n - 1L) rep <- m[[i + 1L]][[1L]]
    start <- i; end <- i + 1L
    subterm <- NULL; exp <- 1
    if (i == n - 1L || identical(rep, 0L)) {
      ai <- args[[i + 1L]]
      if (ai$type == tP && ai$args[[1L]]$type != tS) {
        subterm <- ai$args[[1L]]$args
        l <- length(subterm)
        if (i - l >= 0 && seqeq(sub_args(i - l + 1L, i), subterm) && eqk(ai$args[[1L]], fM(subterm))) { p <- p + 1; start <- start - l }
        if (seqeq(sub_args(i + 2L, min(n, i + 1L + l)), subterm) && eqk(ai$args[[1L]], fM(subterm))) { p <- p + 1; end <- end + l }
      }
      if (p) p <- .add(list(.int(p), ai$args[[2L]])) else next
    } else {
      l <- rep[1L]
      start <- start - (l - 1L)
      subterm <- sub_args(start + 1L, end)
      p <- 2
      end <- end + l
    }
    sk <- skey(subterm)
    if (!is.null(simps[[sk]]) && simps[[sk]] >= start) next
    while (end < n) {
      if (l %in% m[[end]][[1L]]) { p <- if (is_sym(p)) .add(list(p, S_One)) else p + 1; end <- end + l }
      else if (args[[end + 1L]]$type == tP && seqeq(args[[end + 1L]]$args[[1L]]$args, subterm)) {
        p <- .add(list(if (is_sym(p)) p else .int(p), args[[end + 1L]]$args[[2L]])); end <- end + 1L
      } else break
    }
    pre_exp <- S_Zero; pre_arg <- S_One
    if (start - l >= 0 && seqeq(sub_args(start - l + 2L, start), subterm[-1L])) {
      if (subterm[[1L]]$type == tP) { pre_arg <- subterm[[1L]]$args[[1L]]; exp <- subterm[[1L]]$args[[2L]] }
      else { pre_arg <- subterm[[1L]]; exp <- S_One }
      asl <- args[[start - l + 1L]]
      if (asl$type == tP && eqk(asl$args[[1L]], pre_arg)) { pre_exp <- .sub(asl$args[[2L]], exp); start <- start - l; p <- if (is_sym(p)) .add(list(p, S_One)) else p + 1 }
      else if (eqk(asl, pre_arg)) { pre_exp <- .sub(S_One, exp); start <- start - l; p <- if (is_sym(p)) .add(list(p, S_One)) else p + 1 }
    }
    post_exp <- S_Zero; post_arg <- S_One
    if (end + l - 1L < n && seqeq(sub_args(end + 1L, end + l - 1L), subterm[-length(subterm)])) {
      if (subterm[[length(subterm)]]$type == tP) { post_arg <- subterm[[length(subterm)]]$args[[1L]]; exp <- subterm[[length(subterm)]]$args[[2L]] }
      else { post_arg <- subterm[[length(subterm)]]; exp <- S_One }
      ael <- args[[end + l]]
      if (ael$type == tP && eqk(ael$args[[1L]], post_arg)) { post_exp <- .sub(ael$args[[2L]], exp); end <- end + l; p <- if (is_sym(p)) .add(list(p, S_One)) else p + 1 }
      else if (eqk(ael, post_arg)) { post_exp <- .sub(S_One, exp); end <- end + l; p <- if (is_sym(p)) .add(list(p, S_One)) else p + 1 }
    }
    exp <- as_sym(exp)
    if (!is_zero_num(post_exp) && exp$type == "Integer" && .int_mod(exp, 2L) == 0L && start > 0) {
      exp <- .div(exp, S_Two)
      pe_ <- S_One; po_ <- S_One
      asm <- args[[start]]
      if (asm$type == tP && eqk(asm$args[[1L]], post_arg)) { po_ <- .add(list(post_exp, exp)); pe_ <- .sub(asm$args[[2L]], exp) }
      else if (eqk(asm, post_arg)) { po_ <- .add(list(post_exp, exp)); pe_ <- .sub(S_One, exp) }
      if (is_zero_num(pe_) || is_zero_num(po_)) {
        if (is_zero_num(pre_exp)) start <- start - 1L
        post_exp <- po_; pre_exp <- pe_; pre_arg <- post_arg
        subterm <- c(list(fP(post_arg, exp)), subterm[-length(subterm)], list(fP(post_arg, exp)))
      }
    }
    simp_coeff <- simp_coeff + end - start
    if (!is_zero_num(post_exp)) simp_coeff <- simp_coeff - 1
    if (!is_zero_num(pre_exp)) simp_coeff <- simp_coeff - 1
    simps[[skey(subterm)]] <- end
    if (simp_coeff > max_simp_coeff) {
      max_simp_coeff <- simp_coeff
      simp <- list(start, fM(subterm), if (is_sym(p)) p else .int(p), end, l)
      pre <- if (is_zero_num(pre_exp)) S_One else fP(pre_arg, pre_exp)
      post <- if (is_zero_num(post_exp)) S_One else fP(post_arg, post_exp)
    }
  }
  mulp <- function(a, b) if (ismat) { if (is_one_num(a)) b else if (is_one_num(b)) a else .mm(a, b) } else .mul(list(a, b))
  if (!is.null(simp)) {
    subterm <- fP(nc_simplify(simp[[2L]], deep = deep), simp[[3L]])
    pre <- nc_simplify(mulp(if (simp[[1L]] > 0) fM(args[seq_len(simp[[1L]])]) else S_One, pre), deep = deep)
    post <- mulp(post, nc_simplify(if (simp[[4L]] < n) fM(args[(simp[[4L]] + 1L):n]) else S_One, deep = deep))
    res <- mulp(mulp(pre, subterm), post)
    if (!is_one_num(pre) || !is_one_num(post)) res <- nc_simplify(res, deep = FALSE)
  } else res <- fM(args)
  if (invert) res <- fP(res, -1L)
  if (!ismat) {
    f_expr <- factor_nc(expr)
    if (!eqk(f_expr, expr)) res <- compare(res, nc_simplify(f_expr, deep = deep))
  } else res <- doit(res)
  res
}
