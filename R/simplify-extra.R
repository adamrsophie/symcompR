# differentiate_finite (sympy/calculus/finite_diff.py), collect_const
# (sympy/simplify/radsimp.py) and kroneckersimp (sympy/simplify/simplify.py).

#' Finite difference approximation of a derivative
#'
#' `differentiate_finite(expr, ...)` replaces the derivative of `expr` with
#' respect to the variables in `...` by a finite difference of `expr`
#' itself, so the product rule does not hold for it (as SymPy's
#' `differentiate_finite()`). `points` is the step, or a list of points.
#'
#' @param expr An expression.
#' @param ... Variables, each optionally followed by an order.
#' @param points A step or a list of points.
#' @param x0 The point at which the derivative is approximated.
#' @param wrt Unused; kept for SymPy compatibility.
#' @return An expression.
#' @examples
#' x <- Symbol("x"); h <- Symbol("h"); f <- Function("f"); g <- Function("g")
#' differentiate_finite(f(x)*g(x), x, points = list(x - h, x + h))
#' @export
differentiate_finite <- function(expr, ..., points = 1, x0 = NULL, wrt = NULL) {
  e <- as_sym(expr)
  vars <- list(...)
  if (!length(vars)) {
    fs <- free_symbols(e)
    if (length(fs) != 1L) stop("specify the differentiation variable", call. = FALSE)
    vars <- list(fs[[1L]])
  }
  spec <- .diff_spec(vars)
  # an outer derivative merges with the requested one (SymPy's unevaluated diff)
  if (e$type == "Derivative") {
    spec <- c(lapply(e$args[-1L], function(t) list(t$args[[1L]], .int_val(t$args[[2L]]))), spec)
    e <- e$args[[1L]]
  }
  # symbols commute: all occurrences of a variable merge into one order
  merged <- list()
  for (sp in spec) {
    k <- sp[[1L]]$key
    if (is.null(merged[[k]])) merged[[k]] <- sp else merged[[k]][[2L]] <- merged[[k]][[2L]] + sp[[2L]]
  }
  merged <- unname(merged)
  for (sp in merged) {
    v <- sp[[1L]]
    order <- sp[[2L]]
    xc <- if (is.null(x0)) v else as_sym(x0)
    pts <- .fd_points(points, order, v, xc, !is.null(x0))
    e <- apply_finite_diff(order, pts, lapply(pts, function(p) .fd_subs(e, v, p, points)), xc)
  }
  e
}

# The stencil points: a list as given, or steps of h centered at xc (a step
# that depends on v is taken at xc when the center is given).
.fd_points <- function(points, order, v, xc, at_center) {
  seqp <- !is_sym(points) && (is.list(points) || length(points) > 1L)
  pts <- if (seqp) lapply(points, as_sym) else {
    h <- as_sym(points)
    if (at_center) h <- subs(h, v, xc)
    ks <- if (order %% 2L == 0L) seq(-order %/% 2L, order %/% 2L) else seq(-order, order, by = 2L)
    lapply(ks, function(k) .add(list(xc, .mul(list(h, if (order %% 2L == 0L) .int(k) else Rational(k, 2L))))))
  }
  if (length(pts) < order + 1L) stop("Too few points for order ", order, call. = FALSE)
  pts
}

# e at v = p; an embedded derivative in v becomes SymPy's Subs, a finite
# difference over the same points centered at p.
.fd_subs <- function(e, v, p, points) {
  if (identical(p$key, v$key)) return(e)
  ds <- list()
  e <- .map_bottom_up(e, function(z) {
    if (z$type != "Derivative" || !any(vapply(z$args[-1L], function(t) identical(t$args[[1L]]$key, v$key), TRUE))) return(z)
    d <- Dummy("fd")
    ds[[length(ds) + 1L]] <<- list(d, z)
    d
  })
  e <- subs(e, v, p)
  if (!length(ds)) return(e)
  xreplace(e, lapply(ds, function(dz) {
    z <- dz[[2L]]
    # only v is discretized; derivatives in other variables stay
    tv <- Filter(function(t) identical(t$args[[1L]]$key, v$key), z$args[-1L])
    other <- Filter(function(t) !identical(t$args[[1L]]$key, v$key), z$args[-1L])
    f <- if (length(other)) new_sym("Derivative", c(list(z$args[[1L]]), other)) else z$args[[1L]]
    cnt <- as.integer(sum(vapply(tv, function(t) as.numeric(.int_val(t$args[[2L]])), 0)))
    pts <- .fd_points(points, cnt, v, p, TRUE)
    list(dz[[1L]], apply_finite_diff(cnt, pts, lapply(pts, function(pt) subs(f, v, pt)), p))
  }))
}

# SymPy's Factors of a product: list of factor keys -> list(base, exponent),
# with a number split into its sign, numerator and denominator.
.factors_map <- function(m) {
  out <- list()
  put <- function(b, e) {
    k <- b$key
    out[[k]] <<- if (is.null(out[[k]])) list(b, e) else list(b, .add(list(out[[k]][[2L]], e)))
  }
  facs <- if (m$type == "Mul") m$args else list(m)
  # exprtools.Factors keeps the noncommutative part as one ordered factor
  nc <- Filter(function(f) !.is_comm(f), facs)
  if (length(nc)) {
    facs <- Filter(.is_comm, facs)
    put(.mul(nc), S_One)
  }
  for (f in facs) {
    if (is_exact_number(f)) {
      # numbers stay whole: sign, numerator and denominator (exprtools.Factors)
      q <- gmp::as.bigq(f$value)
      if (q < 0) put(S_NegativeOne, S_One)
      p <- abs(gmp::as.bigz(gmp::numerator(q))); d <- gmp::as.bigz(gmp::denominator(q))
      if (p != 1) put(.zs(p), S_One)
      if (d != 1) put(.zs(d), S_NegativeOne)
      next
    }
    be <- .base_exp(f)
    put(be[[1L]], be[[2L]])
  }
  out
}

# m / v when v divides m without turning an integer power fractional.
.factors_div <- function(m, v) {
  fm <- .factors_map(m)
  fv <- .factors_map(v)
  for (k in names(fv)) {
    if (is.null(fm[[k]])) return(NULL)
    d <- .sub(fm[[k]][[2L]], fv[[k]][[2L]])
    if (isTRUE(.facts(d)[["negative"]])) return(NULL)
    if (fm[[k]][[2L]]$type == "Integer" && d$type != "Integer") return(NULL)
    fm[[k]][[2L]] <- d
  }
  .mul(lapply(fm, function(p) .pow(p[[1L]], p[[2L]])))
}

#' Collect common numeric factors
#'
#' `collect_const(expr, ...)` groups the terms of a sum that share the
#' numeric factors in `...` (by default every numeric factor that occurs),
#' as SymPy's `collect_const()`. With `Numbers = FALSE` only non-rational
#' factors such as `sqrt(3)` are collected.
#'
#' @param expr An expression.
#' @param ... Factors to collect.
#' @param Numbers Also collect rational numbers.
#' @return An expression.
#' @examples
#' s <- Symbol("s")
#' collect_const(sqrt(3)*s + sqrt(7)*s + sqrt(3) + sqrt(7))
#' @export
collect_const <- function(expr, ..., Numbers = TRUE) {
  expr <- as_sym(expr)
  if (expr$type != "Add") return(expr)
  vars <- lapply(list(...), as_sym)
  recurse <- !length(vars)
  if (recurse) {
    for (a in expr$args) for (m in if (a$type == "Mul") a$args else list(a))
      if (!length(free_symbols(m))) vars[[m$key]] <- m
    vars <- unname(vars)
  }
  if (!Numbers) vars <- Filter(function(v) !is_number_type(v), vars)
  vars <- .key_sort(lapply(vars, .sort_key), vars)
  i <- 1L
  while (i <= length(vars)) {
    v <- vars[[i]]
    i <- i + 1L
    grouped <- list(); rest <- list()
    for (m in expr$args) {
      q <- .factors_div(m, v)
      if (is.null(q)) rest[[length(rest) + 1L]] <- m else grouped[[length(grouped) + 1L]] <- q
    }
    if (length(grouped) < 2L) next
    inner <- .add(grouped)
    if (recurse && !identical(inner$key, expr$key)) vars[[length(vars) + 1L]] <- inner
    term <- if (Numbers && is_number_type(v) && inner$type == "Add") new_sym("Mul", list(v, inner)) else .mul(list(v, inner))
    expr <- if (length(rest)) .add(c(rest, list(term))) else term
    if (expr$type != "Add") break
  }
  expr
}

#' Simplify products of Kronecker deltas
#'
#' `kroneckersimp(expr)` sets to zero the products of `KroneckerDelta`s that
#' cannot both be one, such as `KroneckerDelta(0, j)*KroneckerDelta(1, j)`.
#'
#' @param expr An expression.
#' @return An expression.
#' @examples
#' j <- Symbol("j")
#' kroneckersimp(1 + KroneckerDelta(0, j)*KroneckerDelta(1, j))
#' @export
kroneckersimp <- function(expr) {
  expr <- as_sym(expr)
  isd <- function(a) a$type == "Function" && a$name == "KroneckerDelta"
  cancels <- function(a1, a2) {
    for (i1 in 1:2) for (i2 in 1:2) {
      if (.is_true(Eq(a1[[i1]], a2[[i2]])) && .is_false(Eq(a1[[3L - i1]], a2[[3L - i2]]))) return(TRUE)
    }
    FALSE
  }
  .map_bottom_up(expr, function(m) {
    if (m$type != "Mul") return(m)
    ds <- Filter(isd, m$args)
    if (length(ds) < 2L) return(m)
    for (p in utils::combn(length(ds), 2L, simplify = FALSE))
      if (cancels(ds[[p[1L]]]$args, ds[[p[2L]]]$args)) return(.mul(list(S_Zero, m)))
    m
  })
}

#' Simplify Bessel functions
#'
#' `besselsimp(expr)` simplifies Bessel functions as SymPy's `besselsimp()`:
#' factors of `I` and `-1` (and their polar forms) are taken out of the
#' arguments of `besselj` and `besseli`, functions of half-integer order are
#' rewritten with trigonometric functions, functions of integer order above
#' one with functions of lower order, the recurrence
#' `J(a, x) = 2*(a + 1)/x*J(a + 1, x) - J(a + 2, x)` removes redundant terms,
#' and a changed expression is factored.
#'
#' @param expr An expression.
#' @return An expression.
#' @examples
#' z <- Symbol("z"); nu <- Symbol("nu")
#' besselsimp(besselj(nu, z * polar_lift(-1)))
#' besselsimp(besseli(Rational(-1, 2), z))
#' besselsimp(z * besseli(0, z) + z * besseli(2, z) / 2 + besseli(1, z))
#' @export
besselsimp <- function(expr) {
  expr <- as_sym(expr)
  orig <- expr
  repl_fun <- function(e, name, fn) .map_bottom_up(e, function(z) {
    if (z$type == "Function" && z$name == name && length(z$args) == 2L) fn(z$args[[1L]], z$args[[2L]]) else z
  })
  mk <- function(name) function(nu, z) fcall(name, list(nu, z))
  replacer <- function(fro, to, factors) function(nu, z) {
    zs <- if (z$type == "Mul") z$args else list(z)
    fk <- vapply(factors, .k, "")
    if (any(vapply(zs, function(a) a$key %in% fk, TRUE))) to(nu, z) else fro(nu, z)
  }
  torewrite <- function(fro, to) function(nu, z) rewrite(fro(nu, z), to)
  tominus <- function(fro) function(nu, z) .mul(list(fcall("exp", list(.mul(list(S_I, S_Pi, nu)))),
                                                    fro(nu, .mul(list(fcall("exp_polar", list(.mul(list(.neg(S_I), S_Pi)))), z)))))
  ifactors <- list(S_I, fcall("exp_polar", list(.mul(list(S_I, S_Pi, S_Half)))), fcall("exp_polar", list(.mul(list(.neg(S_I), S_Pi, S_Half)))))
  expr <- repl_fun(expr, "besselj", replacer(mk("besselj"), torewrite(mk("besselj"), "besseli"), ifactors))
  expr <- repl_fun(expr, "besseli", replacer(mk("besseli"), torewrite(mk("besseli"), "besselj"), ifactors))
  minusfactors <- list(S_NegativeOne, fcall("exp_polar", list(.mul(list(S_I, S_Pi)))))
  expr <- repl_fun(expr, "besselj", replacer(mk("besselj"), tominus(mk("besselj")), minusfactors))
  expr <- repl_fun(expr, "besseli", replacer(mk("besseli"), tominus(mk("besseli")), minusfactors))
  z0 <- Dummy("z")
  expander <- function(name) function(nu, z) {
    m <- fcall("Mod", list(nu, S_One))
    if (identical(m$key, S_Half$key)) {
      r <- expand_func(rewrite(rewrite(fcall(name, list(nu, z0)), "besselj"), "jn"))
      return(simplify(subs(trigsimp(unpolarify(r)), z0, z)))
    }
    if (nu$type == "Integer" && nu$value > 1) return(expand_func(fcall(name, list(nu, z))))
    fcall(name, list(nu, z))
  }
  for (nm in c("besselj", "bessely", "besseli", "besselk")) expr <- repl_fun(expr, nm, expander(nm))
  # the three-term recurrence J(a, x) = 2*(a + 1)/x*J(a + 1, x) - J(a + 2, x)
  use_rec <- function(name, e) {
    repeat {
      bs <- Filter(function(z) z$type == "Function" && z$name == name, unique(atoms(e, "Function")))
      if (!length(bs)) return(e)
      keys <- lapply(bs, function(b) {
        ca <- fcall("re", list(b$args[[1L]]))
        c0 <- if (ca$type == "Add" && is_number_type(ca$args[[1L]])) ca$args[[1L]] else if (is_number_type(ca)) ca else S_Zero
        list(.sort_key(.sub(ca, c0)), as.numeric(num_value_dbl(c0)))
      })
      bs <- .key_sort(keys, bs)
      done <- TRUE
      for (ba in bs) {
        a <- ba$args[[1L]]; xx <- ba$args[[2L]]
        b1 <- fcall(name, list(.add(list(a, S_One)), xx)); b2 <- fcall(name, list(.add(list(a, S_Two)), xx))
        if (.contains_key(e, b1$key) && .contains_key(e, b2$key)) {
          e <- subs(e, ba, .sub(.mul(list(S_Two, .add(list(a, S_One)), .inv(xx), b1)), b2))
          done <- FALSE
          break
        }
      }
      if (done) return(e)
    }
  }
  if (.has_func(expr, "besselj")) expr <- use_rec("besselj", expr)
  if (.has_func(expr, "bessely")) expr <- use_rec("bessely", expr)
  if (!identical(expr$key, orig$key)) expr <- factor_sym(expr)
  expr
}

.has_func <- function(e, name) any(vapply(atoms(e, "Function"), function(z) z$name == name, TRUE))
