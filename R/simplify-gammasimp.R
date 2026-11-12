# gammasimp and combsimp after sympy/simplify/gammasimp.py and
# sympy/simplify/combsimp.py, with the gamma and factorial rewrites of the
# combinatorial functions.

.fn <- function(nm, ...) fcall(nm, list(...))

# expr.rewrite(gamma, piecewise = piecewise)
.rewrite_gamma <- function(e, piecewise = TRUE) .map_bottom_up(e, function(z) {
  if (z$type != "Function") return(z)
  a <- z$args
  M <- function(...) .mul(list(...)); G <- function(x) .fn("gamma", x)
  switch(z$name,
    factorial = G(.add(list(a[[1L]], S_One))),
    binomial = .div(G(.add(list(a[[1L]], S_One))), M(G(.add(list(a[[2L]], S_One))), G(.add(list(.sub(a[[1L]], a[[2L]]), S_One))))),
    RisingFactorial = {
      x <- a[[1L]]; k <- a[[2L]]
      p1 <- .div(G(.add(list(x, k))), G(x))
      p2 <- .div(M(.pow(S_NegativeOne, k), G(.sub(S_One, x))), G(.add(list(.neg(k), .neg(x), S_One))))
      if (!piecewise) { if (.is_true(Le(x, S_Zero))) p2 else p1 }
      else Piecewise(list(p1, Gt(x, S_Zero)), list(p2, S_true))
    },
    FallingFactorial = {
      x <- a[[1L]]; k <- a[[2L]]
      p1 <- .div(G(.add(list(x, S_One))), G(.add(list(.sub(x, k), S_One))))
      p2 <- .div(M(.pow(S_NegativeOne, k), G(.sub(k, x))), G(.neg(x)))
      if (!piecewise) { if (.is_true(Lt(x, S_Zero))) p2 else p1 }
      else Piecewise(list(p1, Ge(x, S_Zero)), list(p2, S_true))
    },
    factorial2 = {
      n <- a[[1L]]
      M(.pow(S_Two, .div(n, S_Two)), G(.add(list(.div(n, S_Two), S_One))),
        Piecewise(list(S_One, Eq(.fn("Mod", n, S_Two), S_Zero)), list(.sqrt(.div(S_Two, S_Pi)), Eq(.fn("Mod", n, S_Two), S_One))))
    },
    subfactorial = {
      n <- a[[1L]]
      M(.add(list(M(.pow(S_NegativeOne, .add(list(n, S_One))), .fn("exp", .neg(M(S_I, S_Pi, n))), .fn("lowergamma", .add(list(n, S_One)), S_NegativeOne)),
                  G(.add(list(n, S_One))))), .fn("exp", S_NegativeOne))
    },
    catalan = { n <- a[[1L]]; .div(M(.pow(.int(4), n), G(.add(list(n, S_Half)))), M(G(S_Half), G(.add(list(n, S_Two))))) },
    beta = if (length(a) == 2L) .div(M(G(a[[1L]]), G(a[[2L]])), G(.add(list(a[[1L]], a[[2L]])))) else z,
    z)
})

# expr.rewrite(factorial)
.rewrite_factorial <- function(e) .map_bottom_up(e, function(z) {
  if (z$type != "Function") return(z)
  a <- z$args
  F <- function(x) .fn("factorial", x)
  switch(z$name,
    gamma = F(.sub(a[[1L]], S_One)),
    binomial = .div(F(a[[1L]]), .mul(list(F(a[[2L]]), F(.sub(a[[1L]], a[[2L]]))))),
    catalan = .div(F(.mul(list(S_Two, a[[1L]]))), .mul(list(F(.add(list(a[[1L]], S_One))), F(a[[1L]])))),
    RisingFactorial = {
      x <- a[[1L]]; k <- a[[2L]]
      if (isTRUE(ask_fact(x, "integer")) && isTRUE(ask_fact(k, "integer")))
        Piecewise(list(.div(F(.sub(.add(list(k, x)), S_One)), F(.sub(x, S_One))), Gt(x, S_Zero)),
                  list(.div(.mul(list(.pow(S_NegativeOne, k), F(.neg(x)))), F(.sub(.neg(k), x))), S_true))
      else z
    },
    FallingFactorial = {
      x <- a[[1L]]; k <- a[[2L]]
      if (isTRUE(ask_fact(x, "integer")) && isTRUE(ask_fact(k, "integer")))
        Piecewise(list(.div(F(x), F(.sub(x, k))), Ge(x, S_Zero)),
                  list(.div(.mul(list(.pow(S_NegativeOne, k), F(.sub(.sub(k, x), S_One)))), F(.sub(.neg(x), S_One))), S_true))
      else z
    },
    z)
})

# The local _rf of gammasimp: a rising factorial that splits off integer
# offsets
.gs_rf <- function(a, b) {
  if (b$type == "Integer") {
    if (is_zero_num(b)) return(S_One)
    n <- as.integer(b$value)
    if (n > 0L) return(.mul(lapply(0:(n - 1L), function(i) .add(list(a, .int(i))))))
    return(.inv(.mul(lapply(1:(-n), function(i) .sub(a, .int(i))))))
  }
  if (b$type == "Add") {
    cb <- .coeff_add_split(b); c <- cb[[1L]]; b_ <- cb[[2L]]
    if (c$type == "Integer" && !is_zero_num(c)) {
      if (num_sign(c) > 0) return(.mul(list(.gs_rf(a, b_), .gs_rf(.add(list(a, b_)), c))))
      return(.div(.gs_rf(a, b_), .gs_rf(.add(list(a, b_, c)), .neg(c))))
    }
  }
  if (a$type == "Add") {
    ca <- .coeff_add_split(a); c <- ca[[1L]]; a_ <- ca[[2L]]
    if (c$type == "Integer" && !is_zero_num(c)) {
      if (num_sign(c) > 0) return(.div(.mul(list(.gs_rf(a_, b), .gs_rf(.add(list(a_, b)), c))), .gs_rf(a_, c)))
      return(.div(.mul(list(.gs_rf(a_, b), .gs_rf(.add(list(a_, c)), .neg(c)))), .gs_rf(.add(list(a_, b, c)), .neg(c))))
    }
  }
  fnode(".gs_rf", list(a, b))
}

.replace_fun <- function(e, name, f) .map_bottom_up(e, function(z) if (z$type == "Function" && z$name == name) do.call(f, z$args) else z)
.count_fun <- function(e, name) { n <- 0L; walk <- function(z) { if (z$type == "Function" && z$name == name) n <<- n + 1L; for (a in z$args) walk(a) }; walk(e); n }

#' Simplify gamma functions and combinatorial expressions
#'
#' `gammasimp(expr)` simplifies products and quotients of gamma functions
#' using the functional equation, the reflection and duplication formulas
#' and the multiplication theorem, as SymPy's `gammasimp()`. `combsimp()`
#' does the same for factorials and binomials, keeping the result in terms
#' of them when the arguments are integers.
#'
#' @param expr An expression.
#' @return An expression.
#' @examples
#' x <- Symbol("x"); n <- Symbol("n", integer = TRUE)
#' gammasimp(gamma(x + 1) / gamma(x))
#' combsimp(factorial(n) / factorial(n - 3))
#' @export
gammasimp <- function(expr) {
  expr <- .rewrite_gamma(as_sym(expr), piecewise = TRUE)
  fs <- unique(c(atoms(expr, "Function"), atoms(expr, "Piecewise")))
  isg <- vapply(fs, function(f) f$type == "Function" && f$name == "gamma", TRUE)
  if (!any(isg)) return(expr)
  f <- fs[!isg]
  if (length(f)) {
    f <- ordered(f)
    dum <- lapply(f, function(i) Dummy("g"))
    simp <- lapply(f, function(fi) sym_rebuild(fi, lapply(fi$args, function(a) .gammasimp(a, as_comb = FALSE))))
    d <- xreplace(expr, Map(list, f, dum))
    return(xreplace(.gammasimp(d, as_comb = FALSE), Map(list, dum, simp)))
  }
  .gammasimp(expr, as_comb = FALSE)
}

#' @rdname gammasimp
#' @export
combsimp <- function(expr) {
  expr <- .rewrite_gamma(as_sym(expr), piecewise = FALSE)
  nonint <- FALSE
  walk <- function(z) {
    if (z$type == "Function" && z$name == "gamma" && !isTRUE(ask_fact(z$args[[1L]], "integer"))) nonint <<- TRUE
    for (a in z$args) walk(a)
  }
  walk(expr)
  if (nonint) return(gammasimp(expr))
  expr <- .gammasimp(expr, as_comb = TRUE)
  .gamma_as_comb(expr)
}

.gamma_as_comb <- function(expr) {
  expr <- .rewrite_factorial(expr)
  f <- function(rv) {
    if (rv$type != "Mul") return(rv)
    bk <- list(); be <- list()
    for (t in rv$args) {
      p <- .abe(t); k <- p[[1L]]$key
      if (is.null(bk[[k]])) { bk[[k]] <- p[[1L]]; be[[k]] <- p[[2L]] } else be[[k]] <- .add(list(be[[k]], p[[2L]]))
    }
    nd <- list(list(), list())
    for (k in names(bk)) {
      b <- bk[[k]]; ex <- be[[k]]
      if (b$type == "Function" && b$name == "factorial" && ex$type == "Integer") {
        m <- as.integer(ex$value)
        if (m > 0L) nd[[1L]] <- c(nd[[1L]], rep(list(b$args[[1L]]), m))
        else nd[[2L]] <- c(nd[[2L]], rep(list(b$args[[1L]]), -m))
        be[[k]] <- S_Zero
      }
    }
    if (!length(nd[[1L]]) || !length(nd[[2L]])) return(rv)
    hit <- FALSE
    addk <- function(b, v) { k <- b$key; if (is.null(bk[[k]])) { bk[[k]] <<- b; be[[k]] <<- S_Zero }; be[[k]] <<- .add(list(be[[k]], v)) }
    for (m in 1:2) {
      i <- 1L
      while (i <= length(nd[[m]])) {
        ai <- nd[[m]][[i]]; found <- FALSE
        if (i < length(nd[[m]])) for (j in (i + 1L):length(nd[[m]])) {
          aj <- nd[[m]][[j]]
          s <- .add(list(ai, aj))
          w <- which(vapply(nd[[3L - m]], function(z) identical(z$key, s$key), TRUE))
          if (length(w)) {
            hit <- TRUE
            nd[[3L - m]] <- nd[[3L - m]][-w[1L]]
            nd[[m]] <- nd[[m]][-c(i, j)]
            addk(.fn("binomial", s, if (count_ops(ai) < count_ops(aj)) ai else aj), if (m == 1L) S_NegativeOne else S_One)
            found <- TRUE
            break
          }
        }
        if (!found) i <- i + 1L
      }
    }
    if (hit) {
      num <- .mul(c(lapply(names(bk), function(k) .pw(bk[[k]], be[[k]])), lapply(nd[[1L]], function(a) .fn("factorial", a))))
      return(.div(num, .mul(lapply(nd[[2L]], function(a) .fn("factorial", a)))))
    }
    rv
  }
  .bottom_up_args(expr, f)
}

.gammasimp <- function(expr, as_comb) {
  expr <- .replace_fun(expr, "gamma", function(n) .gs_rf(S_One, expand(.sub(n, S_One))))
  expr <- if (as_comb) .replace_fun(expr, ".gs_rf", function(a, b) .fn("gamma", .add(list(b, S_One))))
          else .replace_fun(expr, ".gs_rf", function(a, b) .div(.fn("gamma", .add(list(a, b))), .fn("gamma", a)))
  was <- factor_sym(expr)
  expr <- .rule_gamma(was, 0L, as_comb)
  if (!identical(expr$key, was$key)) expr <- factor_sym(expr)
  .replace_fun(expr, "gamma", function(n) if (.is_rat(n)) expand_func(.fn("gamma", n)) else .fn("gamma", n))
}

.gamma_rat <- function(x) {
  was <- .count_fun(x, "gamma")
  xx <- .replace_fun(x, "gamma", function(n) .replace_fun(.gs_rf(S_One, expand(.sub(n, S_One))), ".gs_rf",
                                                           function(a, b) .div(.fn("gamma", .add(list(a, b))), .fn("gamma", a))))
  if (.count_fun(xx, "gamma") < was) xx else x
}

.gamma_factor <- function(x) {
  if (x$type == "Function" && x$name == "gamma") return(TRUE)
  if (x$type %in% c("Add", "Mul")) return(any(vapply(x$args, .gamma_factor, TRUE)))
  if (x$type == "Pow" && (isTRUE(ask_fact(x$args[[2L]], "integer")) || isTRUE(ask_fact(x$args[[1L]], "positive")))) return(.gamma_factor(x$args[[1L]]))
  FALSE
}

.rule_gamma <- function(expr, level, as_comb) {
  rg <- function(z, l) .rule_gamma(z, l, as_comb)
  if (!length(expr$args)) return(expr)
  if (level == 0L) {
    expr <- sym_rebuild(expr, lapply(expr$args, rg, l = 1L))
    level <- 1L
  }
  if (expr$type != "Mul") return(expr)
  if (level == 1L) {
    comm <- vapply(expr$args, .is_comm, TRUE)
    if (!any(comm)) return(expr)
    if (!all(comm)) return(.mul(list(rg(.mul(expr$args[comm]), 2L), .mul(expr$args[!comm]))))
    level <- 2L
  }
  if (level == 2L) {
    gf <- vapply(expr$args, .gamma_factor, TRUE)
    gamma_ind <- .mul(expr$args[!gf]); d <- .mul(expr$args[gf])
    ndd <- .as_numer_denom(d); nd <- ndd[[1L]]; dd <- ndd[[2L]]
    for (ipass in 0:1) {
      args <- ordered(.margs(nd))
      for (i in seq_along(args)) {
        ni <- args[[i]]
        if (ni$type == "Add") {
          r <- .as_numer_denom(.add(lapply(ni$args, function(a) rg(.gamma_rat(.div(a, dd)), 3L))))
          args[[i]] <- r[[1L]]; dd <- r[[2L]]
          if (!.has_function(dd, "gamma")) break
        }
      }
      nd <- .mul(args)
      if (ipass == 0L && !.gamma_factor(nd)) break
      t <- nd; nd <- dd; dd <- t
    }
    expr <- .div(.mul(list(gamma_ind, nd)), dd)
    if (!(expr$type == "Mul" && (.gamma_factor(dd) || .gamma_factor(nd)))) return(expr)
    level <- 3L
  }
  if (level == 3L) {
    repeat {
      was <- expr
      expr <- rg(expr, 4L)
      if (identical(expr$key, was$key)) return(expr)
    }
  }
  ng <- list(); dg <- list(); no <- list(); do <- list()
  explicate <- function(p) {
    if (is_one_num(p) && p$type == "Integer") return(list(NA, list()))
    be <- .abe(p); b <- be[[1L]]; e <- be[[2L]]
    if (e$type == "Integer") {
      m <- as.integer(e$value)
      if (b$type == "Function" && b$name == "gamma") return(list(TRUE, rep(list(b$args[[1L]]), m)))
      return(list(FALSE, rep(list(b), m)))
    }
    list(FALSE, list(p))
  }
  newargs <- ordered(expr$args)
  while (length(newargs)) {
    t <- newargs[[length(newargs)]]; newargs <- newargs[-length(newargs)]
    nd <- .as_numer_denom(t)
    r <- explicate(nd[[1L]])
    if (isTRUE(r[[1L]])) ng <- c(ng, r[[2L]]) else if (isFALSE(r[[1L]])) no <- c(no, r[[2L]])
    r <- explicate(nd[[2L]])
    if (isTRUE(r[[1L]])) dg <- c(dg, r[[2L]]) else if (isFALSE(r[[1L]])) do <- c(do, r[[2L]])
  }
  # the groups are updated in place as SymPy's lists are
  G <- list(ng = ng, dg = dg, no = no, do = do)
  isint <- function(z) z$type == "Integer"
  if (!as_comb) {
    # reflection: gamma(x)*gamma(1 - x) = pi/sin(pi*x)
    for (tr in list(c("ng", "no", "do"), c("dg", "do", "no"))) {
      gam <- G[[tr[1L]]]; new <- list()
      while (length(gam)) {
        g1 <- gam[[length(gam)]]; gam <- gam[-length(gam)]
        if (isTRUE(ask_fact(g1, "integer"))) { new[[length(new) + 1L]] <- g1; next }
        found <- FALSE
        for (i in seq_along(gam)) {
          n <- .sub(.add(list(g1, gam[[i]])), S_One)
          if (!isint(n)) next
          G[[tr[2L]]] <- c(G[[tr[2L]]], list(S_Pi))
          G[[tr[3L]]] <- c(G[[tr[3L]]], list(.fn("sin", .mul(list(S_Pi, g1)))))
          gam <- gam[-i]
          nn <- as.integer(n$value)
          if (nn > 0L) G[[tr[2L]]] <- c(G[[tr[2L]]], lapply(0:(nn - 1L), function(k) .add(list(S_One, .neg(g1), .int(k)))))
          else if (nn < 0L) G[[tr[3L]]] <- c(G[[tr[3L]]], lapply(0:(-nn - 1L), function(k) .sub(.neg(g1), .int(k))))
          found <- TRUE
          break
        }
        if (!found) new[[length(new) + 1L]] <- g1
      }
      G[[tr[1L]]] <- new
    }
    # duplication: gamma(2*y + n) against gamma(y)
    for (tr in list(c("ng", "dg", "no", "do"), c("dg", "ng", "do", "no"))) {
      repeat {
        hitxy <- NULL
        for (xi in seq_along(G[[tr[1L]]])) {
          for (yi in seq_along(G[[tr[2L]]])) {
            n <- .sub(G[[tr[1L]]][[xi]], .mul(list(S_Two, G[[tr[2L]]][[yi]])))
            if (isint(n)) { hitxy <- c(xi, yi); break }
          }
          if (!is.null(hitxy)) break
        }
        if (is.null(hitxy)) break
        y <- G[[tr[2L]]][[hitxy[2L]]]
        n <- .sub(G[[tr[1L]]][[hitxy[1L]]], .mul(list(S_Two, y)))
        G[[tr[1L]]] <- G[[tr[1L]]][-hitxy[1L]]; G[[tr[2L]]] <- G[[tr[2L]]][-hitxy[2L]]
        nn <- as.integer(n$value)
        if (nn > 0L) G[[tr[3L]]] <- c(G[[tr[3L]]], lapply(0:(nn - 1L), function(k) .add(list(.mul(list(S_Two, y)), .int(k)))))
        else if (nn < 0L) G[[tr[4L]]] <- c(G[[tr[4L]]], lapply(0:(-nn - 1L), function(k) .sub(.sub(.mul(list(S_Two, y)), S_One), .int(k))))
        G[[tr[1L]]] <- c(G[[tr[1L]]], list(.add(list(y, S_Half))))
        G[[tr[3L]]] <- c(G[[tr[3L]]], list(.pow(S_Two, .sub(.mul(list(S_Two, y)), S_One))))
        G[[tr[4L]]] <- c(G[[tr[4L]]], list(.sqrt(S_Pi)))
      }
    }
    for (tr in list(c("ng", "no", "do"), c("dg", "do", "no"))) {
      r <- .mult_thm(G[[tr[1L]]], G[[tr[2L]]])
      G[[tr[1L]]] <- r[[1L]]; G[[tr[2L]]] <- r[[2L]]
    }
  }
  # fuzzy matching of gamma arguments against the other factors
  inv <- list()
  compute_ST <- function(z) {
    if (!is.null(inv[[z$key]])) return(inv[[z$key]])
    fun <- vapply(unique(atoms(z, "Function")), .k, "")
    pw <- vapply(atoms(z, "Pow"), function(p) p$args[[2L]]$key, "")
    list(vapply(free_symbols(z), .k, ""), sort(unique(c(fun, pw))))
  }
  update_ST <- function(z) inv[[z$key]] <<- compute_ST(z)
  for (z in c(G$ng, G$dg, G$no, G$do)) update_ST(z)
  find_fuzzy <- function(l, x) {
    if (!length(l)) return(NULL)
    st1 <- compute_ST(x)
    for (y in l) {
      st2 <- inv[[y$key]]
      if (is.null(st2)) st2 <- compute_ST(y)
      if (!identical(st1[[2L]], st2[[2L]]) || (!length(intersect(st1[[1L]], st2[[1L]])) && (length(st1[[1L]]) || length(st2[[1L]])))) next
      a <- length(free_symbols(cancel(.div(x, y)))); b <- length(st1[[1L]]); c <- length(st2[[1L]])
      if (a == 0L && (b > 0L || c > 0L)) return(y)
    }
    NULL
  }
  rm1 <- function(l, y) { w <- which(vapply(l, function(z) identical(z$key, y$key), TRUE))[1L]; l[-w] }
  for (tr in list(c("ng", "no", "do"), c("dg", "do", "no"))) {
    gam <- G[[tr[1L]]]; new <- list()
    while (length(gam)) {
      g <- gam[[length(gam)]]; gam <- gam[-length(gam)]
      cont <- TRUE
      while (cont) {
        cont <- FALSE
        y <- find_fuzzy(G[[tr[2L]]], g)
        if (!is.null(y)) {
          G[[tr[2L]]] <- rm1(G[[tr[2L]]], y)
          if (!identical(y$key, g$key)) { q <- .div(y, g); G[[tr[2L]]] <- c(G[[tr[2L]]], list(q)); update_ST(q) }
          g <- .add(list(g, S_One)); cont <- TRUE
        }
        gm1 <- .sub(g, S_One)
        y <- find_fuzzy(G[[tr[3L]]], gm1)
        if (!is.null(y)) {
          G[[tr[3L]]] <- rm1(G[[tr[3L]]], y)
          if (!identical(y$key, gm1$key)) { q <- .div(gm1, y); G[[tr[2L]]] <- c(G[[tr[2L]]], list(q)); update_ST(q) }
          g <- .sub(g, S_One); cont <- TRUE
        }
      }
      new[[length(new) + 1L]] <- g
    }
    G[[tr[1L]]] <- new
  }
  .mul(list(.div(.mul(lapply(G$ng, function(g) .fn("gamma", g))), .mul(lapply(G$dg, function(g) .fn("gamma", g)))),
            .div(.mul(G$no), .mul(G$do))))
}

# The multiplication theorem of gamma for arguments in arithmetic
# progression with rational step 1/n
.mult_thm <- function(gammas, numer) {
  rk <- list(); rv <- list()
  for (g in gammas) {
    ca <- .coeff_add_split(g); c <- ca[[1L]]; resid <- ca[[2L]]
    k <- resid$key
    if (is.null(rk[[k]])) { rk[[k]] <- resid; rv[[k]] <- list() }
    rv[[k]][[length(rv[[k]]) + 1L]] <- c
  }
  keys <- .key_sort(lapply(rk, .sort_key), names(rk))
  for (k in keys) {
    resid <- rk[[k]]
    coeffs <- rv[[k]]
    coeffs <- coeffs[order(vapply(coeffs, function(c) as.numeric(.sort_num(c)), 0))]
    new <- list()
    repeat {
      run <- .gs_run(coeffs)
      if (is.null(run)) break
      coeffs <- run$coeffs
      n <- run$n; ui <- run$first
      for (u in run$other) {
        con <- .sub(.add(list(resid, u)), S_One)
        m <- as.integer(floor(as.numeric(.sort_num(.sub(u, ui)))))
        if (m > 0L) for (kk in 0:(m - 1L)) numer[[length(numer) + 1L]] <- .sub(con, .int(kk))
      }
      con <- .mul(list(.int(n), .add(list(resid, ui))))
      numer[[length(numer) + 1L]] <- .mul(list(.pow(.mul(list(S_Two, S_Pi)), .div(.int(n - 1L), S_Two)), .pow(.int(n), .sub(S_Half, con))))
      new[[length(new) + 1L]] <- con
    }
    rv[[k]] <- c(lapply(coeffs, function(c) .add(list(resid, c))), new)
  }
  g <- list()
  for (k in keys) g <- c(g, rv[[k]])
  list(g, numer)
}

.gs_run <- function(coeffs) {
  u <- coeffs[!duplicated(vapply(coeffs, .k, ""))]
  frac1 <- function(q) { v <- q$value; v - floor(v) }
  for (i in seq_along(u)) {
    js <- if (i < length(u)) (i + 1L):length(u) else integer()
    dj <- lapply(js, function(j) list(frac1(.sub(u[[j]], u[[i]])), j))
    for (dd in dj) {
      one <- dd[[1L]]
      if (gmp::numerator(one) == 1 && gmp::denominator(one) != 1) {
        n <- as.integer(gmp::denominator(one))
        got <- i; get <- seq_len(n - 1L); done <- FALSE
        for (d2 in dj) {
          m <- n * d2[[1L]]
          if (gmp::denominator(m) == 1 && as.integer(m) %in% get) {
            get <- setdiff(get, as.integer(m)); got <- c(got, d2[[2L]])
            if (!length(get)) { done <- TRUE; break }
          }
        }
        if (!done) next
        vals <- lapply(got, function(j) u[[j]])
        for (c in vals) { w <- which(vapply(coeffs, function(z) identical(z$key, c$key), TRUE))[1L]; coeffs <- coeffs[-w] }
        return(list(n = n, first = vals[[1L]], other = vals[-1L], coeffs = coeffs))
      }
    }
  }
  NULL
}
