# solve_univariate_inequality, ported from sympy/solvers/inequalities.py
# together with its helpers _pt and polytools._nsort. It finds the critical
# points with solvify() and tests a point in each interval between them.

# Python truth value of a relational: TypeError when it does not evaluate
.sset_truth <- function(r) {
  if (.is_true(r)) return(TRUE)
  if (.is_false(r)) return(FALSE)
  .sset_err("TypeError", paste0("cannot determine truth value of Relational: ", .str(r)))
}

# _pt: a point between start and end
.sset_pt <- function(start, end) {
  sinf <- isTRUE(.sset_fact(start, "infinite")); einf <- isTRUE(.sset_fact(end, "infinite"))
  if (!sinf && !einf) return(.div(.add(list(start, end)), S_Two))
  if (sinf && einf) return(S_Zero)
  if ((sinf && is.na(.sset_fact(start, "extended_positive"))) || (einf && is.na(.sset_fact(end, "extended_positive"))))
    .sset_err("ValueError", "cannot proceed with unsigned infinite values")
  if ((einf && isTRUE(.sset_fact(end, "extended_negative"))) || (sinf && isTRUE(.sset_fact(start, "extended_positive")))) {
    tmp <- start; start <- end; end <- tmp
    tmp <- sinf; sinf <- einf; einf <- tmp
  }
  if (einf) {
    if (isTRUE(.sset_fact(start, "extended_positive"))) return(.mul(list(start, S_Two)))
    if (isTRUE(.sset_fact(start, "extended_negative"))) return(.mul(list(start, S_Half)))
    return(.add(list(start, S_One)))
  }
  if (isTRUE(.sset_fact(end, "extended_positive"))) return(.mul(list(end, S_Half)))
  if (isTRUE(.sset_fact(end, "extended_negative"))) return(.mul(list(end, S_Two)))
  .sub(end, S_One)
}

# polytools._nsort(roots, separated = TRUE): list(real roots, complex roots)
.sset_nsort <- function(roots) {
  if (!length(roots)) return(list(list(), list()))
  if (!all(vapply(roots, function(r) isTRUE(.is_number(r)), NA))) .sset_err("NotImplementedError", "")
  key <- lapply(roots, function(r) {
    ri <- .as_real_imag(r)
    vapply(ri, function(p) as.numeric(.as_real_imag(N(p, 2))[[1L]]), 0)
  })
  im <- vapply(key, function(k) as.numeric(k[2L] != 0), 0)
  o <- order(im, vapply(key, `[`, 0, 1L), vapply(key, `[`, 0, 2L), method = "radix")
  roots <- roots[o]; im <- im[o]
  list(roots[im == 0], roots[im == 1])
}

# a TypeError, including the comparisons ops.R rejects with a plain error
.sset_type_error <- function(err)
  inherits(err, "TypeError") || grepl("^Invalid (comparison|NaN comparison)", conditionMessage(err))

.sset_in_list <- function(e, lst) any(vapply(lst, function(z) .sset_eq(z, e), NA))

#' Solve a univariate inequality
#'
#' `solve_univariate_inequality(expr, x)` solves a real inequality in one
#' variable, as SymPy's function of the same name: the critical points come
#' from [solvify()], the singularities and the discontinuities of the domain,
#' and a test point decides each interval between them.
#'
#' @param expr A relational such as `x^2 > 4`.
#' @param x The variable.
#' @param relational Return a relational (default) instead of a set.
#' @param domain A subset of the real line, `Reals` by default.
#' @param continuous Whether `expr` is known to be continuous on `domain`.
#' @return A boolean expression, or a set when `relational = FALSE`.
#' @examples
#' x <- Symbol("x")
#' solve_univariate_inequality(x^2 > 4, x)
#' solve_univariate_inequality((x - 1) / (x + 2) <= 0, x, relational = FALSE)
#' @export
solve_univariate_inequality <- function(expr, x, relational = TRUE, domain = Reals, continuous = FALSE) {
  expr <- as_sym(expr); gen <- as_sym(x); domain <- as_sym(domain)
  if (identical(.sset_subset(domain, Reals), FALSE))
    .sset_err("NotImplementedError", .filldedent("Inequalities in the complex domain are not supported. Try the real domain by setting domain=S.Reals"))
  if (!.sset_set_is(domain, "Reals")) {
    rv <- Intersection(solve_univariate_inequality(expr, gen, relational = FALSE, continuous = continuous), domain)
    if (relational) rv <- as_relational(rv, gen)
    return(rv)
  }
  .gen <- gen
  .domain <- domain
  if (identical(.sset_fact(gen, "extended_real"), FALSE)) {
    rv <- S_EmptySet
    return(if (!relational) rv else as_relational(rv, .gen))
  }
  if (is.na(.sset_fact(gen, "extended_real")) || !(gen$type %in% c("Symbol", "Dummy"))) {
    gen <- Dummy("gen", extended_real = TRUE)
    expr <- xreplace(expr, list(list(.gen, gen)))
  }
  rv <- NULL
  if (.is_true(expr)) {
    rv <- domain
  } else if (.is_false(expr)) {
    rv <- S_EmptySet
  } else {
    func <- function(a, b) .make_rel(expr$type, a, b)
    e <- .sub(expr$args[[1L]], expr$args[[2L]])
    period <- periodicity(e, gen)
    if (!is.null(period) && is_zero_num(period)) {
      e <- expand_mul(e)
      const <- func(e, S_Zero)
      if (.is_true(const)) rv <- domain
      else if (.is_false(const)) rv <- S_EmptySet
    } else if (!is.null(period)) {
      frange <- function_range(e, gen, domain)
      rel <- .REL_OPS[[expr$type]]
      if (rel %in% c("<", "<=")) {
        if (.sset_truth(func(set_sup(frange), S_Zero))) rv <- domain
        else if (!.sset_truth(func(set_inf(frange), S_Zero))) rv <- S_EmptySet
      } else if (rel %in% c(">", ">=")) {
        if (.sset_truth(func(set_inf(frange), S_Zero))) rv <- domain
        else if (!.sset_truth(func(set_sup(frange), S_Zero))) rv <- S_EmptySet
      }
      if (identical(.sub(set_sup(domain), set_inf(domain))$key, S_Infinity$key)) {
        domain <- Intersection(Interval(S_Zero, period, FALSE, TRUE), .domain)
        .domain <- domain
      }
    }
    if (is.null(rv)) {
      nd <- .as_numer_denom(e)
      fail <- function()
        .sset_err("NotImplementedError", .filldedent(paste0("The inequality, ", .str(subs(expr, gen, Symbol("x"))),
                                                            ", cannot be solved using solve_univariate_inequality.")))
      solns <- tryCatch({
        if (!.has1(nd[[1L]], gen) && length(free_symbols(e)) > 1L) .sset_err("ValueError", "")
        s <- solvify(e, gen, domain)
        if (is.null(s)) .sset_err("ValueError", "")
        s
      }, ValueError = function(err) fail(), NotImplementedError = function(err) fail())
      expanded_e <- expand_mul(e)
      valid <- function(xv) {
        v <- subs(expanded_e, gen, expand_mul(xv))
        r <- tryCatch(func(v, S_Zero), error = function(err) if (.sset_type_error(err)) S_false else stop(err))
        if (.is_true(r) || .is_false(r)) return(r)
        if (identical(.sset_fact(v, "extended_real"), FALSE)) return(S_false)
        v <- N(v, 2)
        if (isTRUE(.is_comparable(v))) return(func(v, S_Zero))
        .sset_err("NotImplementedError", paste0("relationship did not evaluate: ", .str(r)))
      }
      singularities <- list()
      for (d in denoms(expr, gen)) singularities <- c(singularities, solvify(d, gen, domain))
      if (!continuous) domain <- continuous_domain(expanded_e, gen, domain)
      include_x <- expr$type %in% c("LessThan", "GreaterThan", "Equality")
      sorted <- tryCatch({
        disc <- Complement(boundary(domain), FiniteSet(set_inf(domain), set_sup(domain)))
        discontinuities <- .sset_elems(disc)
        critical_points <- Intersection(FiniteSet(c(solns, singularities, discontinuities)),
                                        Interval(set_inf(domain), set_sup(domain),
                                                 !.py_in(domain, set_inf(domain)), !.py_in(domain, set_sup(domain))))
        cps <- .sset_elems(critical_points)
        if (all(vapply(cps, function(r) isTRUE(.is_number(r)), NA))) {
          reals <- .sset_nsort(cps)[[1L]]
        } else {
          rl <- vapply(cps, function(r) .sset_fact(r, "extended_real"), NA)
          if (any(is.na(rl))) .sset_err("NotImplementedError", "")
          reals <- cps[rl]
          if (length(reals) > 1L) reals <- .sort_set_elems(reals)
        }
        list(discontinuities, critical_points, reals)
      }, NotImplementedError = function(err) .sset_err("NotImplementedError", "sorting of these roots is not supported"))
      discontinuities <- sorted[[1L]]; critical_points <- sorted[[2L]]; reals <- sorted[[3L]]
      make_real <- Reals
      coeffI <- coeff(expanded_e, S_I)
      check <- TRUE
      if (!is_zero_num(coeffI)) {
        im_sol <- FiniteSet(list())
        im_sol <- tryCatch({
          a <- solveset(coeffI, gen, domain)
          if (a$type != "Interval") {
            for (z in .sset_elems(a))
              if (!.sset_in_list(z, singularities) && .is_true(valid(z)) && isTRUE(.sset_fact(z, "extended_real")))
                im_sol <- Union(im_sol, FiniteSet(z))
          } else {
            start <- set_inf(a); end <- set_sup(a)
            for (z in .sset_nsort(.sset_elems(Union(critical_points, FiniteSet(end))))[[1L]]) {
              valid_start <- .is_true(valid(start))
              if (!.sset_eq(start, end)) {
                valid_z <- .is_true(valid(z))
                pt <- .sset_pt(start, z)
                if (!.sset_in_list(pt, singularities) && isTRUE(.sset_fact(pt, "extended_real")) && .is_true(valid(pt))) {
                  im_sol <- Union(im_sol, Interval(start, z, !valid_start, !valid_z))
                }
              }
              start <- z
            }
            for (s in singularities) im_sol <- Complement(im_sol, FiniteSet(s))
          }
          im_sol
        }, error = function(err) { if (!.sset_type_error(err)) stop(err); check <<- FALSE; Reals })
        if (.sset_is_empty(im_sol))
          .sset_err("ValueError", .filldedent(paste0(.str(subs(expr, gen, .gen)), " contains imaginary parts which cannot be made 0 for any value of ",
                                                     .str(.gen), " satisfying the inequality, leading to relations like I < 0.")))
        make_real <- Intersection(make_real, im_sol)
      }
      sol_sets <- list(S_EmptySet)
      start <- set_inf(domain)
      if (.py_in(domain, start) && .is_true(valid(start)) && isTRUE(.sset_fact(start, "finite")))
        sol_sets[[length(sol_sets) + 1L]] <- FiniteSet(start)
      for (xv in reals) {
        end <- xv
        if (.is_true(valid(.sset_pt(start, end)))) sol_sets[[length(sol_sets) + 1L]] <- Interval(start, end, TRUE, TRUE)
        si <- which(vapply(singularities, function(z) .sset_eq(z, xv), NA))
        if (length(si)) {
          singularities <- singularities[-si[1L]]
        } else {
          di <- which(vapply(discontinuities, function(z) .sset_eq(z, xv), NA))
          if (length(di)) {
            discontinuities <- discontinuities[-di[1L]]
            .valid <- .is_true(valid(xv))
          } else .valid <- include_x
          if (.valid) sol_sets[[length(sol_sets) + 1L]] <- FiniteSet(xv)
        }
        start <- end
      }
      end <- set_sup(domain)
      if (.py_in(domain, end) && .is_true(valid(end)) && isTRUE(.sset_fact(end, "finite")))
        sol_sets[[length(sol_sets) + 1L]] <- FiniteSet(end)
      if (.is_true(valid(.sset_pt(start, end)))) sol_sets[[length(sol_sets) + 1L]] <- Interval(start, end, TRUE, TRUE)
      if (!is_zero_num(coeffI) && check) {
        rv <- Intersection(make_real, .domain)
      } else {
        rv <- subs(Intersection(do.call(Union, sol_sets), make_real, .domain), gen, .gen)
      }
    }
  }
  if (!relational) rv else as_relational(rv, .gen)
}

# Boolean.as_set() for a condition in one real variable
.sset_as_set <- function(cond) {
  if (.is_true(cond)) return(UniversalSet)
  if (.is_false(cond)) return(S_EmptySet)
  free <- free_symbols(cond)
  if (length(free) != 1L)
    .sset_err("NotImplementedError", "Sorry, as_set has not yet been implemented for multivariate expressions")
  x <- free[[1L]]
  rel_set <- function(r) tryCatch(solve_univariate_inequality(r, x, relational = FALSE),
                                  NotImplementedError = function(err) ConditionSet(x, r, Reals))
  rels <- Filter(is_relational, .sset_subexprs(cond))
  reps <- list()
  for (r in rels) {
    per <- periodicity(.sub(r$args[[1L]], r$args[[2L]]), x)
    if (is.null(per) || is_zero_num(per)) next
    s <- rel_set(r)
    if (s$type %in% c("EmptySet", "UniversalSet") || .sset_set_is(s, "Reals")) {
      reps[[length(reps) + 1L]] <- list(r, as_relational(s, x))
      next
    }
    .sset_err("NotImplementedError", .filldedent("as_set is not implemented for relationals with periodic solutions"))
  }
  if (length(reps)) cond <- subs(cond, reps)
  eval_set <- function(c) {
    if (.is_true(c)) return(UniversalSet)
    if (.is_false(c)) return(S_EmptySet)
    switch(c$type,
           And = do.call(Intersection, lapply(c$args, eval_set)),
           Or = do.call(Union, lapply(c$args, eval_set)),
           Not = Complement(Reals, eval_set(c$args[[1L]])),
           if (is_relational(c)) rel_set(c) else .sset_err("NotImplementedError", paste0("as_set is not implemented for ", .str(c))))
  }
  eval_set(cond)
}
