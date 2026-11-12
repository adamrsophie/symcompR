# Factors, Term, gcd_terms, factor_terms and as_content_primitive after
# sympy/core/exprtools.py and the as_content_primitive methods of
# sympy/core/{add,mul,power,numbers,expr}.py.

.isnumber <- function(z) is_sym(z) && is_number_type(z)
.is_rat <- function(z) z$type %in% c("Integer", "Rational")
.num_lt <- function(a, b) {
  d <- .sub(a, b)
  if (is_number_type(d)) return(num_sign(d) < 0)
  v <- ask_fact(d, "negative")
  if (is.na(v)) stop("cannot determine truth value of relational")
  v
}

# SymPy's _keep_coeff
.keep_coeff <- function(coeff, factors, clear = TRUE, sign = FALSE) {
  if (!is_number_type(coeff)) {
    if (is_number_type(factors)) { t <- factors; factors <- coeff; coeff <- t }
    else return(.mul(list(coeff, factors)))
  }
  if (is_one_num(factors) && factors$type == "Integer") return(coeff)
  if (is_one_num(coeff) && coeff$type == "Integer") return(factors)
  if (identical(coeff$key, S_NegativeOne$key) && !sign) return(.neg(factors))
  if (factors$type == "Add") {
    if (!clear && coeff$type == "Rational") {
      args <- lapply(factors$args, function(i) { sp <- as_coeff_mul_split(i); list(.keep_coeff(sp[[1L]], coeff), sp[[2L]]) })
      if (any(vapply(args, function(a) a[[1L]]$type == "Integer", TRUE)))
        return(.add(lapply(args, function(a) if (is_one_num(a[[1L]])) a[[2L]] else .mul(list(a[[1L]], a[[2L]])))))
    }
    return(new_sym("Mul", list(coeff, factors)))
  }
  if (factors$type == "Mul") {
    margs <- factors$args
    if (is_number_type(margs[[1L]])) {
      margs[[1L]] <- .mul(list(margs[[1L]], coeff))
      if (is_one_num(margs[[1L]])) margs <- margs[-1L]
    } else margs <- c(list(coeff), margs)
    if (length(margs) == 1L) return(margs[[1L]])
    if (any(vapply(margs, function(m) m$type == "Add", TRUE)) && is_number_type(margs[[1L]])) return(new_sym("Mul", margs))
    return(.mul(margs))
  }
  m <- .mul(list(coeff, factors))
  if (is_number_type(m) && !is_number_type(factors)) m <- new_sym("Mul", list(coeff, factors))
  m
}

# as_coeff_Mul(rational = TRUE)
.coeff_mul_rat <- function(e) {
  sp <- as_coeff_mul_split(e)
  if (.is_rat(sp[[1L]])) sp else list(S_One, e)
}
.margs <- function(e) if (e$type == "Mul") e$args else list(e)
.aargs <- function(e) if (e$type == "Add") e$args else list(e)

# Add.primitive
.add_primitive <- function(e) add_primitive(e)

#' Content and primitive part
#'
#' `as_content_primitive(e)` returns `list(c, p)` with `e = c*p`, `c` a
#' positive rational and `p` with integer coefficients of gcd one, as
#' SymPy's `Expr.as_content_primitive()`. `radical = TRUE` also takes out
#' common radicals of a sum; `clear = FALSE` keeps a denominator inside a
#' sum when some term would keep an integer coefficient.
#'
#' @param e An expression.
#' @param radical,clear See above.
#' @return A list of the content and the primitive part.
#' @examples
#' x <- Symbol("x")
#' as_content_primitive((2*x + 4)^2)
#' as_content_primitive(2*sqrt(2) + 4*sqrt(10), radical = TRUE)
#' @export
as_content_primitive <- function(e, radical = FALSE, clear = TRUE) .acp(as_sym(e), radical, clear)

.acp <- function(e, radical = FALSE, clear = TRUE) {
  t <- e$type
  if (t %in% c("Integer", "Rational")) {
    if (is_zero_num(e)) return(list(S_One, e))
    if (num_sign(e) > 0) return(list(e, S_One))
    return(list(.neg(e), S_NegativeOne))
  }
  if (t == "Add") return(.acp_add(e, radical, clear))
  if (t == "Mul") {
    coef <- S_One; args <- list()
    for (a in e$args) {
      r <- .acp(a, radical, clear)
      coef <- .mul(list(coef, r[[1L]]))
      if (!(is_one_num(r[[2L]]) && r[[2L]]$type == "Integer")) args[[length(args) + 1L]] <- r[[2L]]
    }
    return(list(coef, .mul(args)))
  }
  if (t == "Pow" || (t == "Function" && e$name == "exp")) return(.acp_pow(e, radical, clear))
  list(S_One, e)
}

.acp_add <- function(e, radical, clear) {
  terms <- lapply(e$args, function(a) { r <- .acp(a, radical, clear); .keep_coeff(r[[1L]], r[[2L]]) })
  pr <- .add_primitive(.add(terms))
  con <- pr[[1L]]; prim <- pr[[2L]]
  if (!clear && con$type != "Integer" && prim$type == "Add") {
    nd <- .as_numer_denom(con); d <- nd[[2L]]
    pp <- .div(prim, d)
    if (any(vapply(.aargs(pp), function(a) as_coeff_mul_split(a)[[1L]]$type == "Integer", TRUE))) {
      prim <- pp; con <- nd[[1L]]
    } else con <- .div(nd[[1L]], d)
  }
  if (radical && prim$type == "Add") {
    args <- prim$args
    rads <- list(); common_q <- NULL; ok <- TRUE
    for (m in args) {
      term_rads <- list()
      for (ai in .margs(m)) {
        if (ai$type == "Pow") {
          be <- .abe(ai); b <- be[[1L]]; ex <- be[[2L]]
          if (.is_rat(ex) && b$type == "Integer") {
            q <- as.character(gmp::denominator(ex$value))
            term_rads[[q]] <- c(term_rads[[q]], list(abs(gmp::as.bigz(b$value))^as.integer(gmp::numerator(ex$value))))
          }
        }
      }
      if (!length(term_rads)) { ok <- FALSE; break }
      common_q <- if (is.null(common_q)) names(term_rads) else intersect(common_q, names(term_rads))
      if (!length(common_q)) { ok <- FALSE; break }
      rads[[length(rads) + 1L]] <- term_rads
    }
    if (ok) {
      G <- list()
      for (q in common_q) {
        vals <- lapply(rads, function(r) Reduce(`*`, r[[q]]))
        g <- Reduce(gmp::gcd.bigz, vals, gmp::as.bigz(0))
        if (g != 1) G[[length(G) + 1L]] <- .pow(.int(g), make_rational(gmp::as.bigq(1, gmp::as.bigz(q))))
      }
      if (length(G)) {
        Gm <- .mul(G)
        prim <- .mul(list(Gm, .add(lapply(args, function(ai) .div(ai, Gm)))))
      }
    }
  }
  list(con, prim)
}

.acp_pow <- function(e, radical, clear) {
  be <- .abe(e); b <- be[[1L]]; ex <- be[[2L]]
  rb <- .acp(b, radical, clear); b <- .keep_coeff(rb[[1L]], rb[[2L]])
  re <- .acp(ex, radical, clear); ce <- re[[1L]]; pe <- re[[2L]]
  if (.is_rat(b)) {
    ht <- .coeff_add_split(pe); h <- ht[[1L]]; tt <- ht[[2L]]
    if (.is_rat(h) && !is_zero_num(b)) {
      ceh <- .mul(list(ce, h))
      c <- .pw(b, ceh)
      r <- S_Zero
      if (!.is_rat(c)) {
        p <- gmp::numerator(ceh$value); q <- gmp::denominator(ceh$value)
        iceh <- p %/% q; r <- .int(p %% q)
        c <- .pw(b, .int(iceh))
      }
      cq <- .int(gmp::denominator(ceh$value))
      return(list(c, .pw(b, .keep_coeff(ce, .add(list(tt, .div(.div(r, ce), cq)))))))
    }
  }
  ex <- .keep_coeff(ce, pe)
  if (.is_rat(ex) && b$type == "Mul") {
    ht <- .acp(b, radical, clear); h <- ht[[1L]]; tt <- ht[[2L]]
    cm <- as_coeff_mul_split(.pw(h, ex)); c <- cm[[1L]]; m <- cm[[2L]]
    mme <- .abe(m); m <- mme[[1L]]; me <- mme[[2L]]
    if ((is_one_num(m) && m$type == "Integer") || identical(me$key, ex$key))
      return(list(c, .pw(.keep_coeff(m, tt), ex)))
  }
  list(S_One, .pw(b, ex))
}

# as_coeff_Add: the numeric term and the rest
.coeff_add_split <- function(e) {
  if (e$type == "Add" && is_number_type(e$args[[1L]])) return(list(e$args[[1L]], .add(e$args[-1L])))
  if (is_number_type(e)) return(list(e, S_Zero))
  list(S_Zero, e)
}

# ---------------------------------------------------------------------------
# Factors: an ordered map from bases to exponents

.F <- function(b = list(), e = list()) structure(list(b = b, e = e), class = "symfactors")
.F_zero <- function() .F(list2(S_Zero), list2(S_One))
list2 <- function(z) stats::setNames(list(z), z$key)
.F_is_zero <- function(f) length(f$b) == 1L && identical(names(f$b), S_Zero$key)
.F_put <- function(f, b, e) { f$b[[b$key]] <- b; f$e[[b$key]] <- e; f }
.F_del <- function(f, k) { f$b[[k]] <- NULL; f$e[[k]] <- NULL; f }

# Factors(expr)
.F_new <- function(x = NULL) {
  if (is.null(x) || inherits(x, "symfactors")) return(if (is.null(x)) .F() else x)
  if (is_one_num(x) && x$type == "Integer") return(.F())
  if (is_zero_num(x) && is_number_type(x)) return(.F_zero())
  if (is_number_type(x)) {
    f <- .F(); n <- x
    if (num_sign(n) < 0) { f <- .F_put(f, S_NegativeOne, S_One); n <- .neg(n) }
    if (!(is_one_num(n) && n$type == "Integer")) {
      if (n$type %in% c("Float", "Integer") || is_const(n, "oo")) f <- .F_put(f, n, S_One)
      else if (n$type == "Rational") {
        if (gmp::numerator(n$value) != 1) f <- .F_put(f, .int(gmp::numerator(n$value)), S_One)
        f <- .F_put(f, .int(gmp::denominator(n$value)), S_NegativeOne)
      } else stop("Expected Float|Rational|Integer")
    }
    return(f)
  }
  if (!length(x$args)) return(.F_put(.F(), x, S_One))
  fs <- .margs(x)
  comm <- vapply(fs, .is_comm, TRUE)
  c <- fs[comm]; nc <- fs[!comm]
  isI <- vapply(c, function(z) is_const(z, "I"), TRUE)
  ni <- sum(isI); c <- c[!isI]
  f <- .F()
  for (t in c) {
    be <- if (.is_rat(t)) list(t, S_One) else .abe(t)
    k <- be[[1L]]$key
    f <- .F_put(f, be[[1L]], if (is.null(f$e[[k]])) be[[2L]] else .add(list(f$e[[k]], be[[2L]])))
  }
  for (k in names(f$b)) {
    r <- f$b[[k]]
    if (r$type == "Rational") {
      ex <- f$e[[k]]; f <- .F_del(f, k)
      p <- .int(gmp::numerator(r$value)); q <- .int(gmp::denominator(r$value))
      f <- .F_put(f, p, .add(list(if (is.null(f$e[[p$key]])) S_Zero else f$e[[p$key]], ex)))
      f <- .F_put(f, q, .sub(if (is.null(f$e[[q$key]])) S_Zero else f$e[[q$key]], ex))
    }
  }
  if (ni) f <- .F_put(f, S_I, .add(list(if (is.null(f$e[[S_I$key]])) S_Zero else f$e[[S_I$key]], .int(ni))))
  if (length(nc)) f <- .F_put(f, new_sym("Mul", nc), S_One)
  f
}

# Factors(dict): numeric powers of I, -1 and 1 normalized
.F_dict <- function(f) {
  hk <- names(f$b)[vapply(f$b, function(k) is_const(k, "I") || (k$type == "Integer" && abs(k$value) == 1), TRUE)]
  if (!length(hk)) return(f)
  i1 <- S_One
  for (k in hk) {
    if (!.isnumber(f$e[[k]])) next
    i1 <- .mul(list(i1, .pow(f$b[[k]], f$e[[k]])))
    f <- .F_del(f, k)
  }
  if (!(is_one_num(i1) && i1$type == "Integer")) {
    for (a in .margs(i1)) {
      if (identical(a$key, S_NegativeOne$key)) f <- .F_put(f, a, S_One)
      else if (is_const(a, "I")) f <- .F_put(f, S_I, S_One)
      else if (a$type == "Pow") f <- .F_put(f, a$args[[1L]], .add(list(if (is.null(f$e[[a$args[[1L]]$key]])) S_Zero else f$e[[a$args[[1L]]$key]], a$args[[2L]])))
      else if (is_number_type(a) && is_one_num(a)) f <- .F_put(f, a, S_One)
      else if (is_number_type(a) && is_one_num(.neg(a))) { f <- .F_put(f, .neg(a), S_One); f <- .F_put(f, S_NegativeOne, S_One) }
      else stop("unexpected factor in i1")
    }
  }
  f
}

.F_as_expr <- function(f) {
  args <- list()
  for (k in names(f$b)) {
    b <- f$b[[k]]; ex <- f$e[[k]]
    if (!(is_one_num(ex) && ex$type == "Integer")) {
      if (ex$type == "Integer") {
        be <- .abe(b)
        args[[length(args) + 1L]] <- .pw(be[[1L]], .keep_coeff(ex, be[[2L]]))
      } else args[[length(args) + 1L]] <- .pw(b, ex)
    } else args[[length(args) + 1L]] <- b
  }
  .mul(args)
}

.F_mul <- function(a, b) {
  b <- .F_new(b)
  if (.F_is_zero(a) || .F_is_zero(b)) return(.F_zero())
  f <- a
  for (k in names(b$b)) {
    ex <- b$e[[k]]
    if (!is.null(f$b[[k]])) {
      ex <- .add(list(f$e[[k]], ex))
      if (is_zero_num(ex)) { f <- .F_del(f, k); next }
    }
    f <- .F_put(f, b$b[[k]], ex)
  }
  .F_dict(f)
}

.F_normal <- function(a, b) {
  if (!inherits(b, "symfactors")) {
    b <- .F_new(b)
    if (.F_is_zero(b)) return(list(.F(), .F_zero()))
    if (.F_is_zero(a)) return(list(.F_zero(), .F()))
  }
  sf <- a; of <- b
  for (k in names(a$b)) {
    if (is.null(b$b[[k]])) next
    se <- a$e[[k]]; oe <- b$e[[k]]
    ex <- .sub(se, oe)
    if (is_zero_num(ex)) { sf <- .F_del(sf, k); of <- .F_del(of, k) }
    else if (.isnumber(ex)) {
      if (num_sign(ex) > 0) { sf$e[[k]] <- ex; of <- .F_del(of, k) }
      else { sf <- .F_del(sf, k); of$e[[k]] <- .neg(ex) }
    } else {
      r <- extract_additively(se, oe)
      if (!is.null(r)) {
        if (!is_zero_num(r)) { sf$e[[k]] <- r; of <- .F_del(of, k) }
        else { sf <- .F_del(sf, k); of <- .F_del(of, k) }
      } else {
        sa <- .coeff_add_split(se); sc <- sa[[1L]]
        if (!is_zero_num(sc)) {
          oa <- .coeff_add_split(oe); oc <- oa[[1L]]
          diff <- .sub(sc, oc)
          if (num_sign(diff) > 0) { sf$e[[k]] <- .sub(sf$e[[k]], oc); oe <- oa[[2L]] }
          else if (num_sign(diff) < 0) { sf$e[[k]] <- .sub(sf$e[[k]], sc); of$e[[k]] <- .sub(of$e[[k]], sc); oe <- .sub(oa[[2L]], diff) }
          else { sf$e[[k]] <- sa[[2L]]; oe <- oa[[2L]] }
        }
        if (!is_zero_num(oe)) of$e[[k]] <- oe else of <- .F_del(of, k)
      }
    }
  }
  list(.F_dict(sf), .F_dict(of))
}

.F_div <- function(a, b) {
  quo <- a; rem <- .F()
  if (!inherits(b, "symfactors")) {
    b <- .F_new(b)
    if (.F_is_zero(b)) stop("division by zero")
    if (.F_is_zero(a)) return(list(.F_zero(), .F()))
  }
  for (k in names(b$b)) {
    ex <- b$e[[k]]
    if (!is.null(quo$b[[k]])) {
      d <- .sub(quo$e[[k]], ex)
      if (.isnumber(d)) {
        if (num_sign(d) <= 0) quo <- .F_del(quo, k)
        if (num_sign(d) >= 0) { if (!is_zero_num(d)) quo <- .F_put(quo, b$b[[k]], d); next }
        ex <- .neg(d)
      } else {
        r <- extract_additively(quo$e[[k]], ex)
        if (!is.null(r)) {
          if (!is_zero_num(r)) quo$e[[k]] <- r else quo <- .F_del(quo, k)
        } else {
          oe <- ex
          sa <- .coeff_add_split(quo$e[[k]]); sc <- sa[[1L]]
          if (!is_zero_num(sc)) {
            oa <- .coeff_add_split(oe); oc <- oa[[1L]]
            diff <- .sub(sc, oc)
            if (num_sign(diff) > 0) { quo$e[[k]] <- .sub(quo$e[[k]], oc); oe <- oa[[2L]] }
            else if (num_sign(diff) < 0) { quo$e[[k]] <- .sub(quo$e[[k]], sc); oe <- .sub(oa[[2L]], diff) }
            else { quo$e[[k]] <- sa[[2L]]; oe <- oa[[2L]] }
          }
          if (!is_zero_num(oe)) rem <- .F_put(rem, b$b[[k]], oe)
        }
        next
      }
    }
    rem <- .F_put(rem, b$b[[k]], ex)
  }
  list(.F_dict(quo), .F_dict(rem))
}

.F_pow <- function(a, n) {
  if (n < 0) stop("expected non-negative integer")
  f <- .F()
  if (n != 0) for (k in names(a$b)) f <- .F_put(f, a$b[[k]], .mul(list(a$e[[k]], .int(n))))
  .F_dict(f)
}

.F_gcd <- function(a, b) {
  if (!inherits(b, "symfactors")) {
    b <- .F_new(b)
    if (.F_is_zero(b)) return(a)
  }
  f <- .F()
  for (k in names(a$b)) {
    if (is.null(b$b[[k]])) next
    lt <- ask_fact(.sub(a$e[[k]], b$e[[k]]), "negative")
    if (isTRUE(lt)) f <- .F_put(f, a$b[[k]], a$e[[k]])
    else if (isFALSE(lt)) f <- .F_put(f, a$b[[k]], b$e[[k]])
  }
  .F_dict(f)
}

.F_lcm <- function(a, b) {
  if (!inherits(b, "symfactors")) {
    b <- .F_new(b)
    if (.F_is_zero(a) || .F_is_zero(b)) return(.F_zero())
  }
  f <- a
  for (k in names(b$b)) {
    ex <- b$e[[k]]
    if (!is.null(f$b[[k]])) { cur <- f$e[[k]]; ex <- if (.num_lt(cur, ex)) ex else cur }
    f <- .F_put(f, b$b[[k]], ex)
  }
  .F_dict(f)
}

# ---------------------------------------------------------------------------
# Term: coefficient, numerator and denominator factors

.T <- function(coeff, numer = .F(), denom = .F()) list(coeff = coeff, numer = numer, denom = denom)

.T_new <- function(term) {
  if (!.is_comm(term)) stop("commutative expression expected")
  sp <- as_coeff_mul_split(term); coeff <- sp[[1L]]
  fs <- if (is_one_num(sp[[2L]]) && sp[[2L]]$type == "Integer") list() else .margs(sp[[2L]])
  nb <- list(); ne <- list(); db <- list(); de <- list()
  for (fct in fs) {
    dp <- .decompose_power(fct); base <- dp[[1L]]; ex <- dp[[2L]]
    if (base$type == "Add") {
      pr <- .add_primitive(base); base <- pr[[2L]]
      coeff <- .mul(list(coeff, .pow(pr[[1L]], .int(ex))))
    }
    k <- base$key
    if (ex > 0) { nb[[k]] <- base; ne[[k]] <- (if (is.null(ne[[k]])) 0 else ne[[k]]) + as.numeric(ex) }
    else { db[[k]] <- base; de[[k]] <- (if (is.null(de[[k]])) 0 else de[[k]]) - as.numeric(ex) }
  }
  mk <- function(b, e) .F_dict(.F(b, lapply(e, .int)))
  .T(coeff, mk(nb, ne), mk(db, de))
}

.T_as_expr <- function(t) .mul(list(t$coeff, .div(.F_as_expr(t$numer), .F_as_expr(t$denom))))
.T_mul <- function(a, b) {
  nd <- .F_normal(.F_mul(a$numer, b$numer), .F_mul(a$denom, b$denom))
  .T(.mul(list(a$coeff, b$coeff)), nd[[1L]], nd[[2L]])
}
.T_inv <- function(a) .T(.inv(a$coeff), a$denom, a$numer)
.T_quo <- function(a, b) .T_mul(a, .T_inv(b))
.T_gcd <- function(a, b) .T(.num_gcd(a$coeff, b$coeff), .F_gcd(a$numer, b$numer), .F_gcd(a$denom, b$denom))

# ---------------------------------------------------------------------------
# gcd_terms and factor_terms

.gcd_terms_ <- function(terms, isprimitive = FALSE, fraction = TRUE) {
  if (is_sym(terms)) terms <- .aargs(terms)
  terms <- Filter(function(t) !is_zero_num(t), terms)
  terms <- lapply(terms, .T_new)
  if (!length(terms)) return(list(S_Zero, S_Zero, S_One))
  if (length(terms) == 1L) {
    cont <- terms[[1L]]$coeff
    numer <- .F_as_expr(terms[[1L]]$numer); denom <- .F_as_expr(terms[[1L]]$denom)
  } else {
    cont <- terms[[1L]]
    for (t in terms[-1L]) cont <- .T_gcd(cont, t)
    terms <- lapply(terms, function(t) .T_quo(t, cont))
    if (fraction) {
      denom <- terms[[1L]]$denom
      for (t in terms[-1L]) denom <- .F_lcm(denom, t$denom)
      numers <- lapply(terms, function(t) .mul(list(t$coeff, .F_as_expr(.F_mul(t$numer, .F_div(denom, t$denom)[[1L]])))))
    } else {
      numers <- lapply(terms, .T_as_expr)
      denom <- .F()
    }
    cont <- .T_as_expr(cont)
    numer <- .add(numers)
    denom <- .F_as_expr(denom)
  }
  if (!isprimitive && numer$type == "Add") {
    pr <- .add_primitive(numer); cont <- .mul(list(cont, pr[[1L]])); numer <- pr[[2L]]
  }
  list(cont, numer, denom)
}

#' Common factors of terms
#'
#' `factor_terms(e)` pulls out the factors shared by all terms of every sum
#' in `e`, without expanding: `3*x*y*(x + 2*y)`. `gcd_terms(e)` does this for
#' the top level sum and also puts it over a common denominator, as SymPy's
#' functions of the same names.
#'
#' @param e An expression.
#' @param isprimitive The terms are already primitive.
#' @param clear Clear a rational coefficient from a sum (`FALSE` keeps a
#'   denominator inside when some term keeps an integer coefficient).
#' @param fraction Put the result over a common denominator.
#' @param radical Also take out common radicals.
#' @param sign Take out a common negative sign.
#' @return An expression.
#' @examples
#' x <- Symbol("x"); y <- Symbol("y")
#' factor_terms(3*x^2*y + 6*x*y^2)
#' gcd_terms(x / 2 + 1 / x)
#' @name factor_terms
NULL

#' @rdname factor_terms
#' @export
gcd_terms <- function(e, isprimitive = FALSE, clear = TRUE, fraction = TRUE) {
  if (inherits(e, "symmatrix")) return(.mat_map(e, function(z) gcd_terms(z, isprimitive, clear, fraction)))
  addlike <- is.list(e) && !is_sym(e)
  if (!addlike) e <- as_sym(e)
  if (addlike || e$type == "Add") {
    terms <- if (addlike) lapply(e, as_sym) else e$args
    reps <- list()
    terms <- lapply(terms, function(a) {
      if (.is_comm(a)) return(a)
      fs <- .margs(a); comm <- vapply(fs, .is_comm, TRUE)
      d <- Dummy("d")
      reps[[length(reps) + 1L]] <<- list(d, .mul(fs[!comm]))
      .mul(c(fs[comm], list(d)))
    })
    r <- .gcd_terms_(terms, isprimitive, fraction)
    cont <- r[[1L]]; numer <- r[[2L]]; denom <- r[[3L]]
    if (length(reps)) numer <- xreplace(numer, reps)
    cf <- as_coeff_mul_split(cont); coeff <- cf[[1L]]; factors <- cf[[2L]]
    if (!clear) {
      c2 <- as_coeff_mul_split(coeff); c <- c2[[1L]]; c_ <- c2[[2L]]
      if (c$type != "Integer" && numer$type == "Add") {
        nd <- .as_numer_denom(c)
        nn <- .div(numer, nd[[2L]])
        if (any(vapply(.aargs(nn), function(a) as_coeff_mul_split(a)[[1L]]$type == "Integer", TRUE))) {
          numer <- nn; coeff <- .mul(list(nd[[1L]], c_))
        }
      }
    }
    return(.keep_coeff(coeff, .div(.mul(list(factors, numer)), denom), clear = clear))
  }
  if (!length(e$args)) return(e)
  if (e$type == "Mul") {
    sp <- as_coeff_mul_split(e)
    rest <- if (is_one_num(sp[[2L]]) && sp[[2L]]$type == "Integer") list() else .margs(sp[[2L]])
    return(.keep_coeff(sp[[1L]], .mul(lapply(rest, gcd_terms, isprimitive = isprimitive, clear = clear, fraction = fraction)), clear = clear))
  }
  sym_rebuild(e, lapply(e$args, function(a) if (length(a$args) || is_sym(a)) gcd_terms(a, isprimitive, clear, fraction) else a))
}

#' @rdname factor_terms
#' @export
factor_terms <- function(e, radical = FALSE, clear = FALSE, fraction = FALSE, sign = TRUE) {
  if (inherits(e, "symmatrix")) return(.mat_map(e, function(z) factor_terms(z, radical, clear, fraction, sign)))
  .factor_terms(as_sym(e), radical, clear, fraction, sign)
}

.factor_terms <- function(e, radical, clear, fraction, sign) {
  do <- function(z) .factor_terms(z, radical, clear, fraction, sign)
  if (!length(e$args)) return(e)
  if (e$type %in% c("Pow", "Function", "Tuple") || is_relational(e) ||
      e$type %in% c("And", "Or", "Not", "Xor", "Implies", "Equivalent", "Piecewise", "ExprCondPair")) {
    args <- lapply(e$args, do)
    if (identical(vapply(args, .k, ""), vapply(e$args, .k, ""))) return(e)
    return(sym_rebuild(e, args))
  }
  if (e$type %in% c("Sum", "Integral")) return(.factor_sum_int(e, radical, clear, fraction, sign))
  cp <- .acp(e, radical, clear); cont <- cp[[1L]]; p <- cp[[2L]]
  if (p$type == "Add") {
    la <- lapply(p$args, do)
    if (all(vapply(la, function(a) !is.null(extract_multiplicatively(as_coeff_mul_split(a)[[1L]], S_NegativeOne)), TRUE))) {
      cont <- .neg(cont); la <- lapply(la, .neg)
    }
    special <- list()
    for (i in seq_along(la)) {
      ex <- .abe(la[[i]])[[2L]]
      if (ex$type == "Mul" && !identical(ex$key, .mul(ex$args)$key)) {
        d <- Dummy("d"); special[[length(special) + 1L]] <- list(d, la[[i]]); la[[i]] <- d
      }
    }
    p <- gcd_terms(.add(la), isprimitive = TRUE, clear = clear, fraction = fraction)
    if (length(special)) p <- xreplace(p, special)
  } else if (length(p$args)) {
    p <- sym_rebuild(p, lapply(p$args, do))
  }
  .keep_coeff(cont, p, clear = clear, sign = sign)
}

.factor_sum_int <- function(e, radical, clear, fraction, sign) {
  result <- e$args[[1L]]
  if (is_zero_num(result)) return(S_Zero)
  limits <- e$args[-1L]
  wrt <- lapply(limits, function(l) if (l$type == "Tuple") l$args[[1L]] else l)
  f <- .factor_terms(result, radical, clear, fraction, sign)
  ts <- .margs(f)
  if (f$type == "Add") ts <- f$args
  dep <- function(t) any(vapply(wrt, function(w) .has1(t, w), TRUE))
  if (f$type == "Add") {
    dd <- vapply(f$args, dep, TRUE)
    i <- .add(f$args[!dd]); d <- .add(f$args[dd])
    return(.add(list(.mul(list(i, sym_rebuild(e, c(list(S_One), limits)))), sym_rebuild(e, c(list(d), limits)))))
  }
  dd <- vapply(.margs(f), dep, TRUE)
  i <- .mul(.margs(f)[!dd]); d <- .mul(.margs(f)[dd])
  .mul(list(i, sym_rebuild(e, c(list(d), limits))))
}
