# Hypergeometric and Meijer G functions, a port of
# sympy/functions/special/hyper.py: the hyper and meijerg constructors, their
# derivatives, convergence data, and the HyperRep_* functions that
# hyperexpand() returns before rewriting them to elementary functions.

# Parameters as a list of sym: a Tuple, an R list or an atomic vector.
.tuple_args <- function(v) {
  if (is_sym(v)) return(if (v$type == "Tuple") v$args else list(v))
  lapply(as.list(v), as_sym)
}

# _prep_tuple: hypergeometric functions are unbranched in their parameters
.prep_tuple <- function(v) do.call(Tuple, lapply(v, unpolarify))

.ordered_list <- function(v) if (length(v) > 1L) ordered(v) else v

# hyper.__new__: parameters common to ap and bq cancel (Counter difference),
# and each list is put in ordered() order
.hyper_canon <- function(ap, bq) {
  ka <- vapply(ap, `[[`, "", "key"); kb <- vapply(bq, `[[`, "", "key")
  for (k in intersect(ka, kb)) {
    m <- min(sum(ka == k), sum(kb == k))
    ia <- which(ka == k)[seq_len(m)]; ib <- which(kb == k)[seq_len(m)]
    ap <- ap[-ia]; ka <- ka[-ia]; bq <- bq[-ib]; kb <- kb[-ib]
  }
  srt <- function(v, keys) {
    u <- .ordered_list(v[!duplicated(keys)])
    unlist(lapply(u, function(x) rep(list(x), sum(keys == x$key))), recursive = FALSE)
  }
  list(srt(ap, ka), srt(bq, kb))
}

.hyper_eval <- function(args) {
  cn <- .hyper_canon(.tuple_args(args[[1L]]), .tuple_args(args[[2L]]))
  z <- args[[3L]]
  p <- length(cn[[1L]]); q <- length(cn[[2L]])
  # hyper.eval: the argument is unbranched where the series converges
  if (p <= q || (p == q + 1L && .is_true(Le(fcall("Abs", list(z)), S_One)))) z <- unpolarify(z)
  fnode("hyper", list(.prep_tuple(cn[[1L]]), .prep_tuple(cn[[2L]]), z))
}

.hy_ap <- function(e) e$args[[1L]]$args
.hy_bq <- function(e) e$args[[2L]]$args

# TupleParametersBase._eval_derivative: hyper has no parameter derivative
.hyper_diffx <- function(e, x) {
  if (.has1(e$args[[1L]], x) || .has1(e$args[[2L]], x)) return(.derivative_node(e, list(list(x, 1L))))
  ap <- .hy_ap(e); bq <- .hy_bq(e)
  fac <- .div(.mul(ap), .mul(bq))
  d <- .mul(list(fac, hyper(lapply(ap, function(a) .add(list(a, S_One))), lapply(bq, function(b) .add(list(b, S_One))), e$args[[3L]])))
  .mul(list(d, .diff1(e$args[[3L]], x)))
}

.pfq_numeric <- function(ap, bq, z) {
  for (v in ap) {
    j <- match(v, bq)
    if (!is.na(j)) { ap <- ap[-match(v, ap)]; bq <- bq[-j] }
  }
  if (Mod(z) >= 1 && !any(ap <= 0 & ap == trunc(ap))) stop("hyper() series needs |z| < 1")
  .pfq(ap, bq, z)
}

def_function("hyper", nargs = 3L, eval = .hyper_eval, diffx = .hyper_diffx, numeric = .pfq_numeric)

#' @rdname special-functions-extra
#' @export
hyper <- function(ap, bq, z) {
  tup <- function(v) do.call(Tuple, .tuple_args(v))
  fcall("hyper", list(tup(ap), tup(bq), as_sym(z)))
}

# hyper.eta, radius_of_convergence and convergence_statement
.hyper_eta <- function(e) .sub(.add(.hy_ap(e)), .add(.hy_bq(e)))

.hyper_radius <- function(e) {
  ap <- .hy_ap(e); bq <- .hy_bq(e)
  ap_nonpos <- Filter(function(a) isTRUE(ask_fact(a, "integer")) && isTRUE(ask_fact(a, "nonpositive")), ap)
  aints <- Filter(function(a) a$type == "Integer", ap_nonpos)
  bints <- Filter(function(b) b$type == "Integer" && num_sign(b) <= 0, bq)
  if (length(ap_nonpos) && all(vapply(bq, function(b) identical(ask_fact(b, "integer"), FALSE) || identical(ask_fact(b, "nonpositive"), FALSE), TRUE)))
    return(S_Infinity)
  if (length(aints) || length(bints)) {
    if (length(aints) < length(bints)) return(S_Zero)
    popped <- FALSE
    for (b in bints) {
      cancelled <- FALSE
      while (length(aints)) {
        a <- aints[[length(aints)]]; aints <- aints[-length(aints)]
        if (a$value >= b$value) { cancelled <- TRUE; break }
        popped <- TRUE
      }
      if (!cancelled) return(S_Zero)
    }
    if (length(aints) || popped) return(S_Infinity)
  }
  if (length(ap) == length(bq) + 1L) return(S_One)
  if (length(ap) <= length(bq)) return(S_Infinity)
  S_Zero
}

.hyper_convergence <- function(e) {
  R <- .hyper_radius(e)
  if (is_zero_num(R)) return(S_false)
  if (identical(R$key, S_Infinity$key)) return(S_true)
  et <- .hyper_eta(e); z <- e$args[[3L]]
  re_e <- fcall("re", list(et)); az <- fcall("Abs", list(z))
  c1 <- And(Lt(re_e, 0), Le(az, 1))
  c2 <- And(Le(0, re_e), Lt(re_e, 1), Le(az, 1), Ne(z, 1))
  c3 <- And(Ge(re_e, 1), Lt(az, 1))
  Or(c1, c2, c3)
}

# hyper._eval_expand_func
.expand_hooks$hyper <- function(e) {
  ap <- .hy_ap(e); bq <- .hy_bq(e)
  if (length(ap) == 2L && length(bq) == 1L && is_one_num(e$args[[3L]])) {
    a <- ap[[1L]]; b <- ap[[2L]]; c <- bq[[1L]]
    g <- function(v) fcall("gamma", list(v))
    return(.div(.div(.mul(list(g(c), g(.sub(.sub(c, a), b)))), g(.sub(c, a))), g(.sub(c, b))))
  }
  hyperexpand(e)
}

# hyper._eval_rewrite_as_Sum
.hyper_as_sum <- function(ap, bq, z) {
  e <- fnode("hyper", list(ap, bq, z))
  n <- Dummy("n", integer = TRUE)
  rf <- function(a) fcall("RisingFactorial", list(a, n))
  coeff <- .div(.mul(lapply(ap$args, rf)), .mul(lapply(bq$args, rf)))
  term <- .div(.mul(list(coeff, .pow(z, n))), fcall("factorial", list(n)))
  Piecewise(list(Sum(term, list(n, S_Zero, S_Infinity)), .hyper_convergence(e)), list(e, S_true))
}

# ---------------------------------------------------------------------------
# meijerg

.mg_an <- function(e) e$args[[1L]]$args[[1L]]$args
.mg_aother <- function(e) e$args[[1L]]$args[[2L]]$args
.mg_bm <- function(e) e$args[[2L]]$args[[1L]]$args
.mg_bother <- function(e) e$args[[2L]]$args[[2L]]$args
.mg_ap <- function(e) c(.mg_an(e), .mg_aother(e))
.mg_bq <- function(e) c(.mg_bm(e), .mg_bother(e))
.mg_nu <- function(e) .sub(.add(.mg_bq(e)), .add(.mg_ap(e)))
.mg_delta <- function(e) .sub(.int(length(.mg_bm(e)) + length(.mg_an(e))), .rat(length(.mg_ap(e)) + length(.mg_bq(e)), 2L))

# meijerg.__new__
.meijerg_eval <- function(args) {
  if (length(args) == 5L) args <- list(Tuple(args[[1L]], args[[2L]]), Tuple(args[[3L]], args[[4L]]), args[[5L]])
  tr <- function(p) {
    ps <- .tuple_args(p)
    if (length(ps) != 2L) stop("wrong argument", call. = FALSE)
    do.call(Tuple, lapply(ps, function(i) .prep_tuple(.ordered_list(.tuple_args(i)))))
  }
  a0 <- tr(args[[1L]]); a1 <- tr(args[[2L]])
  if (.has_any(Tuple(a0, a1), list(S_Infinity, S_ComplexInfinity, S_NegativeInfinity)))
    stop("G-function parameters must be finite", call. = FALSE)
  for (a in a0$args[[1L]]$args) for (b in a1$args[[1L]]$args) {
    d <- .sub(a, b)
    if (d$type == "Integer" && num_sign(d) > 0)
      stop("no parameter a1, ..., an may differ from any b1, ..., bm by a positive integer", call. = FALSE)
  }
  fnode("meijerg", list(a0, a1, args[[3L]]))
}

.mg_build <- function(an, aother, bm, bother, z)
  fcall("meijerg", list(Tuple(do.call(Tuple, an), do.call(Tuple, aother)), Tuple(do.call(Tuple, bm), do.call(Tuple, bother)), z))

# meijerg.fdiff(3)
.mg_fdiff_z <- function(e) {
  an <- .mg_an(e); bm <- .mg_bm(e); z <- e$args[[3L]]
  if (length(an)) {
    a <- an; a[[1L]] <- .sub(a[[1L]], S_One)
    G <- .mg_build(a, .mg_aother(e), bm, .mg_bother(e), z)
    return(.mul(list(.inv(z), .add(list(.mul(list(.sub(an[[1L]], S_One), e)), G)))))
  }
  if (length(bm)) {
    b <- bm; b[[1L]] <- .add(list(b[[1L]], S_One))
    G <- .mg_build(an, .mg_aother(e), b, .mg_bother(e), z)
    return(.mul(list(.inv(z), .sub(.mul(list(bm[[1L]], e)), G))))
  }
  S_Zero
}

.mod1_zero <- function(d) is_zero_num(fcall("Mod", list(simplify(d), S_One)))

# meijerg._diff_wrt_parameter (idx is 0 based, over ap + bq)
.mg_diff_param <- function(e, idx) {
  an <- .mg_an(e); ap <- .mg_aother(e); bm <- .mg_bm(e); bq <- .mg_bother(e); z <- e$args[[3L]]
  if (idx < length(an)) an <- an[-(idx + 1L)]
  else {
    idx <- idx - length(an)
    if (idx < length(ap)) ap <- ap[-(idx + 1L)]
    else {
      idx <- idx - length(ap)
      if (idx < length(bm)) bm <- bm[-(idx + 1L)] else bq <- bq[-(idx - length(bm) + 1L)]
    }
  }
  pairs_of <- function(l1, l2) {
    pairs <- list()
    while (length(l1)) {
      x <- l1[[length(l1)]]; l1 <- l1[-length(l1)]
      found <- NULL
      for (i in seq_along(l2)) if (.mod1_zero(.sub(x, l2[[i]]))) { found <- i; break }
      if (is.null(found)) stop("Derivative not expressible as G-function?", call. = FALSE)
      pairs[[length(pairs) + 1L]] <- list(x, l2[[found]])
      l2 <- l2[-found]
    }
    pairs
  }
  pairs1 <- pairs_of(an, bq); pairs2 <- pairs_of(ap, bm)
  an0 <- .mg_an(e); ap0 <- .mg_aother(e); bm0 <- .mg_bm(e); bq0 <- .mg_bother(e)
  res <- .mul(list(fcall("log", list(z)), e))
  for (pr in pairs1) {
    sgn <- 1L; n <- .sub(pr[[1L]], pr[[2L]]); base <- pr[[2L]]
    if (.is_true(Lt(n, 0))) { sgn <- -1L; n <- .sub(pr[[2L]], pr[[1L]]); base <- pr[[1L]] }
    for (k in seq_len(as.integer(.int_val(n))) - 1L)
      res <- .sub(res, .mul(list(.int(sgn), .mg_build(c(an0, list(.add(list(base, .int(k + 1L))))), ap0, bm0, c(bq0, list(.add(list(base, .int(k))))), z))))
  }
  for (pr in pairs2) {
    sgn <- 1L; n <- .sub(pr[[2L]], pr[[1L]]); base <- pr[[1L]]
    if (.is_true(Lt(n, 0))) { sgn <- -1L; n <- .sub(pr[[1L]], pr[[2L]]); base <- pr[[2L]] }
    for (k in seq_len(as.integer(.int_val(n))) - 1L)
      res <- .sub(res, .mul(list(.int(sgn), .mg_build(an0, c(ap0, list(.add(list(base, .int(k + 1L))))), c(bm0, list(.add(list(base, .int(k))))), bq0, z))))
  }
  res
}

# TupleParametersBase._eval_derivative for meijerg
.meijerg_diffx <- function(e, x) {
  res <- S_Zero
  if (.has1(e$args[[1L]], x) || .has1(e$args[[2L]], x)) {
    ps <- c(.mg_ap(e), .mg_bq(e))
    for (i in seq_along(ps)) {
      m <- .diff1(ps[[i]], x)
      if (!is_zero_num(m)) {
        r <- tryCatch(.mg_diff_param(e, i - 1L), error = function(err) NULL)
        if (is.null(r)) return(.derivative_node(e, list(list(x, 1L))))
        res <- .add(list(res, .mul(list(r, m))))
      }
    }
  }
  .add(list(res, .mul(list(.mg_fdiff_z(e), .diff1(e$args[[3L]], x)))))
}

def_function("meijerg", nargs = c(3L, 5L), eval = .meijerg_eval, diffx = .meijerg_diffx)

#' @rdname special-functions-extra
#' @export
meijerg <- function(...) {
  args <- list(...)
  if (length(args) == 5L) args <- list(list(args[[1L]], args[[2L]]), list(args[[3L]], args[[4L]]), args[[5L]])
  if (length(args) != 3L) stop("args must be either as, as', bs, bs', z or as, bs, z", call. = FALSE)
  pair <- function(p) {
    ps <- if (is_sym(p) && p$type == "Tuple") p$args else as.list(p)
    if (length(ps) != 2L) stop("wrong argument", call. = FALSE)
    Tuple(do.call(Tuple, .tuple_args(ps[[1L]])), do.call(Tuple, .tuple_args(ps[[2L]])))
  }
  fcall("meijerg", list(pair(args[[1L]]), pair(args[[2L]]), as_sym(args[[3L]])))
}

.expand_hooks$meijerg <- function(e) hyperexpand(e)

# meijerg.get_period
.mg_get_period <- function(e) {
  compute <- function(l) {
    for (i in seq_along(l)) {
      if (!.is_rat(l[[i]])) return(S_Infinity)
      for (j in seq_along(l)[seq_along(l) > i]) if (.mod1_zero(.sub(l[[i]], l[[j]]))) return(S_Infinity)
    }
    if (!length(l)) return(S_One)
    .int(Reduce(function(a, b) gmp::lcm.bigz(a, b), lapply(l, function(x) gmp::denominator(x$value))))
  }
  beta <- compute(.mg_bm(e)); alpha <- compute(.mg_an(e))
  p <- length(.mg_ap(e)); q <- length(.mg_bq(e))
  twopi <- .mul(list(S_Two, S_Pi))
  if (p == q) {
    if (identical(alpha$key, S_Infinity$key) || identical(beta$key, S_Infinity$key)) return(S_Infinity)
    return(.mul(list(twopi, .int(gmp::lcm.bigz(gmp::numerator(alpha$value), gmp::numerator(beta$value))))))
  }
  .mul(list(twopi, if (p < q) beta else alpha))
}

# meijerg.integrand
.mg_integrand <- function(e, s) {
  g <- function(v) fcall("gamma", list(v))
  num <- c(list(.pow(e$args[[3L]], s)), lapply(.mg_bm(e), function(b) g(.sub(b, s))),
           lapply(.mg_an(e), function(a) g(.add(list(S_One, .neg(a), s)))))
  den <- c(lapply(.mg_bother(e), function(b) g(.add(list(S_One, .neg(b), s)))), lapply(.mg_aother(e), function(a) g(.sub(a, s))))
  .div(.mul(num), .mul(den))
}

# ---------------------------------------------------------------------------
# HyperRep functions: branched representations used by hyperexpand(). Each has
# _expr_small, _expr_small_minus, _expr_big and _expr_big_minus; rewriting to
# "nonrep" or "nonrepsmall" (rules in rewrite-hyper.R) picks the branch from
# the polar argument.

.HYPERREP <- list()

.def_hyperrep <- function(name, nargs, small, small_minus, big, big_minus) {
  .HYPERREP[[name]] <<- list(small = small, small_minus = small_minus, big = big, big_minus = big_minus)
  def_function(name, nargs = nargs, eval = function(args) {
    # HyperRep.eval: unbranched in all but the last argument
    k <- length(args)
    newargs <- c(lapply(args[-k], unpolarify), args[k])
    if (any(vapply(seq_len(k - 1L), function(i) !identical(newargs[[i]]$key, args[[i]]$key), TRUE)))
      return(fnode(name, newargs))
    NULL
  })
  invisible(name)
}

.hyperrep_nonrep <- function(name, args) {
  f <- .HYPERREP[[name]]
  k <- length(args)
  bf <- extract_branch_factor(args[[k]], allow_half = TRUE)
  x <- bf[[1L]]; n <- bf[[2L]]
  newargs <- c(args[-k], list(x))
  minus <- FALSE
  if (n$type != "Integer") { minus <- TRUE; n <- .sub(n, S_Half) }
  if (minus) {
    small <- do.call(f$small_minus, newargs); big <- do.call(f$big_minus, c(newargs, list(n)))
  } else {
    small <- do.call(f$small, newargs); big <- do.call(f$big, c(newargs, list(n)))
  }
  if (identical(big$key, small$key)) return(small)
  Piecewise(list(big, Gt(fcall("Abs", list(x)), S_One)), list(small, S_true))
}

.hyperrep_nonrepsmall <- function(name, args) {
  f <- .HYPERREP[[name]]
  k <- length(args)
  bf <- extract_branch_factor(args[[k]], allow_half = TRUE)
  newargs <- c(args[-k], list(bf[[1L]]))
  if (bf[[2L]]$type != "Integer") return(do.call(f$small_minus, newargs))
  do.call(f$small, newargs)
}

.hr_is <- function(n, fact) isTRUE(ask_fact(n, fact))
.ipi <- function(...) .mul(c(list(S_Pi, S_I), list(...)))

.hr_power1_small <- function(a, x) .pow(.sub(S_One, x), a)
.hr_power1_small_minus <- function(a, x) .pow(.add(list(S_One, x)), a)
.hr_power1_big <- function(a, x, n) {
  if (.hr_is(a, "integer")) return(.hr_power1_small(a, x))
  .mul(list(.pow(.sub(x, S_One), a), fcall("exp", list(.ipi(.sub(.mul(list(S_Two, n)), S_One), a)))))
}
.hr_power1_big_minus <- function(a, x, n) {
  if (.hr_is(a, "integer")) return(.hr_power1_small_minus(a, x))
  .mul(list(.pow(.add(list(S_One, x)), a), fcall("exp", list(.ipi(S_Two, n, a)))))
}
.def_hyperrep("HyperRep_power1", 2L, .hr_power1_small, .hr_power1_small_minus, .hr_power1_big, .hr_power1_big_minus)

.p2 <- function(a) .pow(S_Two, .sub(.mul(list(S_Two, a)), S_One))
.def_hyperrep("HyperRep_power2", 2L,
  function(a, x) .mul(list(.p2(a), .pow(.add(list(S_One, sqrt(.sub(S_One, x)))), .sub(S_One, .mul(list(S_Two, a)))))),
  function(a, x) .mul(list(.p2(a), .pow(.add(list(S_One, sqrt(.add(list(S_One, x))))), .sub(S_One, .mul(list(S_Two, a)))))),
  function(a, x, n) {
    sgn <- S_NegativeOne
    if (.hr_is(n, "odd")) { sgn <- S_One; n <- .sub(n, S_One) }
    .mul(list(.p2(a), .pow(.add(list(S_One, .mul(list(sgn, S_I, sqrt(.sub(x, S_One)))))), .sub(S_One, .mul(list(S_Two, a)))),
              fcall("exp", list(.neg(.ipi(S_Two, n, a))))))
  },
  function(a, x, n) {
    sgn <- if (.hr_is(n, "odd")) S_NegativeOne else S_One
    .mul(list(sgn, .p2(a), .pow(.add(list(sqrt(.add(list(S_One, x))), sgn)), .sub(S_One, .mul(list(S_Two, a)))),
              fcall("exp", list(.neg(.ipi(S_Two, a, n))))))
  })

.def_hyperrep("HyperRep_log1", 1L,
  function(x) fcall("log", list(.sub(S_One, x))),
  function(x) fcall("log", list(.add(list(S_One, x)))),
  function(x, n) .add(list(fcall("log", list(.sub(x, S_One))), .ipi(.sub(.mul(list(S_Two, n)), S_One)))),
  function(x, n) .add(list(fcall("log", list(.add(list(S_One, x)))), .ipi(S_Two, n))))

.def_hyperrep("HyperRep_atanh", 1L,
  function(x) .div(fcall("atanh", list(sqrt(x))), sqrt(x)),
  function(x) .div(fcall("atan", list(sqrt(x))), sqrt(x)),
  function(x, n) {
    h <- .div(.ipi(), S_Two)
    .div(if (.hr_is(n, "even")) .add(list(fcall("acoth", list(sqrt(x))), h)) else .sub(fcall("acoth", list(sqrt(x))), h), sqrt(x))
  },
  function(x, n) {
    if (.hr_is(n, "even")) return(.div(fcall("atan", list(sqrt(x))), sqrt(x)))
    .div(.sub(fcall("atan", list(sqrt(x))), S_Pi), sqrt(x))
  })

.hr_asin1 <- list(
  small = function(z) .div(fcall("asin", list(sqrt(z))), sqrt(z)),
  small_minus = function(z) .div(fcall("asinh", list(sqrt(z))), sqrt(z)),
  big = function(z, n) .mul(list(.pow(S_NegativeOne, n), .add(list(.div(.mul(list(.sub(S_Half, n), S_Pi)), sqrt(z)),
                                                                     .div(.mul(list(S_I, fcall("acosh", list(sqrt(z))))), sqrt(z)))))),
  big_minus = function(z, n) .mul(list(.pow(S_NegativeOne, n), .add(list(.div(fcall("asinh", list(sqrt(z))), sqrt(z)),
                                                                           .div(.mul(list(n, S_Pi, S_I)), sqrt(z)))))))
.def_hyperrep("HyperRep_asin1", 1L, .hr_asin1$small, .hr_asin1$small_minus, .hr_asin1$big, .hr_asin1$big_minus)

.def_hyperrep("HyperRep_asin2", 1L,
  function(z) .div(.hr_asin1$small(z), .hr_power1_small(S_Half, z)),
  function(z) .div(.hr_asin1$small_minus(z), .hr_power1_small_minus(S_Half, z)),
  function(z, n) .div(.hr_asin1$big(z, n), .hr_power1_big(S_Half, z, n)),
  function(z, n) .div(.hr_asin1$big_minus(z, n), .hr_power1_big_minus(S_Half, z, n)))

.e2pi <- function(...) fcall("exp", list(.mul(c(list(S_Two, S_Pi, S_I), list(...)))))
.def_hyperrep("HyperRep_sqrts1", 2L,
  function(a, z) .div(.add(list(.pow(.sub(S_One, sqrt(z)), .mul(list(S_Two, a))), .pow(.add(list(S_One, sqrt(z))), .mul(list(S_Two, a))))), S_Two),
  function(a, z) .mul(list(.pow(.add(list(S_One, z)), a), fcall("cos", list(.mul(list(.mul(list(S_Two, a)), fcall("atan", list(sqrt(z))))))))),
  function(a, z, n) {
    if (.hr_is(n, "even"))
      return(.div(.add(list(.mul(list(.pow(.add(list(sqrt(z), S_One)), .mul(list(S_Two, a))), .e2pi(n, a))),
                            .mul(list(.pow(.sub(sqrt(z), S_One), .mul(list(S_Two, a))), .e2pi(.sub(n, S_One), a))))), S_Two))
    n <- .sub(n, S_One)
    .div(.add(list(.mul(list(.pow(.sub(sqrt(z), S_One), .mul(list(S_Two, a))), .e2pi(a, .add(list(n, S_One))))),
                   .mul(list(.pow(.add(list(sqrt(z), S_One)), .mul(list(S_Two, a))), .e2pi(a, n))))), S_Two)
  },
  function(a, z, n) {
    at <- .mul(list(.mul(list(S_Two, a)), fcall("atan", list(sqrt(z)))))
    if (.hr_is(n, "even")) return(.mul(list(.pow(.add(list(S_One, z)), a), .e2pi(n, a), fcall("cos", list(at)))))
    .mul(list(.pow(.add(list(S_One, z)), a), .e2pi(n, a), fcall("cos", list(.sub(at, .mul(list(S_Two, S_Pi, a)))))))
  })

.def_hyperrep("HyperRep_sqrts2", 2L,
  function(a, z) .div(.mul(list(sqrt(z), .sub(.pow(.sub(S_One, sqrt(z)), .mul(list(S_Two, a))), .pow(.add(list(S_One, sqrt(z))), .mul(list(S_Two, a)))))), S_Two),
  function(a, z) .mul(list(sqrt(z), .pow(.add(list(S_One, z)), a), fcall("sin", list(.mul(list(.mul(list(S_Two, a)), fcall("atan", list(sqrt(z))))))))),
  function(a, z, n) {
    if (.hr_is(n, "even"))
      return(.mul(list(.div(sqrt(z), S_Two), .sub(.mul(list(.pow(.sub(sqrt(z), S_One), .mul(list(S_Two, a))), .e2pi(a, .sub(n, S_One)))),
                                                   .mul(list(.pow(.add(list(sqrt(z), S_One)), .mul(list(S_Two, a))), .e2pi(a, n)))))))
    n <- .sub(n, S_One)
    .mul(list(.div(sqrt(z), S_Two), .sub(.mul(list(.pow(.sub(sqrt(z), S_One), .mul(list(S_Two, a))), .e2pi(a, .add(list(n, S_One))))),
                                         .mul(list(.pow(.add(list(sqrt(z), S_One)), .mul(list(S_Two, a))), .e2pi(a, n))))))
  },
  function(a, z, n) {
    at <- .mul(list(.mul(list(S_Two, a)), fcall("atan", list(sqrt(z)))))
    if (.hr_is(n, "even")) return(.mul(list(.pow(.add(list(S_One, z)), a), .e2pi(n, a), sqrt(z), fcall("sin", list(at)))))
    .mul(list(.pow(.add(list(S_One, z)), a), .e2pi(n, a), sqrt(z), fcall("sin", list(.sub(at, .mul(list(S_Two, S_Pi, a)))))))
  })

.def_hyperrep("HyperRep_log2", 1L,
  function(z) fcall("log", list(.add(list(S_Half, .div(sqrt(.sub(S_One, z)), S_Two))))),
  function(z) fcall("log", list(.add(list(S_Half, .div(sqrt(.add(list(S_One, z))), S_Two))))),
  function(z, n) {
    base <- .add(list(.ipi(.sub(n, S_Half)), fcall("log", list(.div(sqrt(z), S_Two)))))
    as <- .mul(list(S_I, fcall("asin", list(.inv(sqrt(z))))))
    if (.hr_is(n, "even")) .add(list(base, as)) else .sub(base, as)
  },
  function(z, n) {
    if (.hr_is(n, "even")) return(.add(list(.ipi(n), fcall("log", list(.add(list(S_Half, .div(sqrt(.add(list(S_One, z))), S_Two))))))))
    .add(list(.ipi(n), fcall("log", list(.sub(.div(sqrt(.add(list(S_One, z))), S_Two), S_Half)))))
  })

.def_hyperrep("HyperRep_cosasin", 2L,
  function(a, z) fcall("cos", list(.mul(list(.mul(list(S_Two, a)), fcall("asin", list(sqrt(z))))))),
  function(a, z) fcall("cosh", list(.mul(list(.mul(list(S_Two, a)), fcall("asinh", list(sqrt(z))))))),
  function(a, z, n) fcall("cosh", list(.add(list(.mul(list(.mul(list(S_Two, a)), fcall("acosh", list(sqrt(z))))), .ipi(a, .sub(.mul(list(S_Two, n)), S_One)))))),
  function(a, z, n) fcall("cosh", list(.add(list(.mul(list(.mul(list(S_Two, a)), fcall("asinh", list(sqrt(z))))), .ipi(S_Two, a, n))))))

.def_hyperrep("HyperRep_sinasin", 2L,
  function(a, z) .mul(list(.div(sqrt(z), sqrt(.sub(S_One, z))), fcall("sin", list(.mul(list(.mul(list(S_Two, a)), fcall("asin", list(sqrt(z))))))))),
  function(a, z) .mul(list(.neg(.div(sqrt(z), sqrt(.add(list(S_One, z))))), fcall("sinh", list(.mul(list(.mul(list(S_Two, a)), fcall("asinh", list(sqrt(z))))))))),
  function(a, z, n) .mul(list(.neg(.inv(sqrt(.sub(S_One, .inv(z))))),
                              fcall("sinh", list(.add(list(.mul(list(.mul(list(S_Two, a)), fcall("acosh", list(sqrt(z))))), .ipi(a, .sub(.mul(list(S_Two, n)), S_One)))))))),
  function(a, z, n) .mul(list(.neg(.inv(sqrt(.add(list(S_One, .inv(z)))))),
                              fcall("sinh", list(.add(list(.mul(list(.mul(list(S_Two, a)), fcall("asinh", list(sqrt(z))))), .ipi(S_Two, a, n))))))))
