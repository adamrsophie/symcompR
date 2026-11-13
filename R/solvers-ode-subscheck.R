# Checking solutions of ODEs by substitution, after
# sympy/solvers/ode/subscheck.py: sub_func_doit, checkodesol and
# checksysodesol.

# subscheck.sub_func_doit: func replaced by new, its derivatives evaluated
.ode_sub_func_doit <- function(eq, func, new) {
  reps <- list(list(func, new))
  for (d in atoms(eq, "Derivative")) {
    if (identical(d$args[[1L]]$key, func$key)) {
      reps[[length(reps) + 1L]] <- list(d, .diff_spec_eval(new, .ode_variable_count(d)))
    } else {
      reps[[length(reps) + 1L]] <- list(d, doit(xreplace(d, list(list(func, new))), deep = FALSE))
    }
  }
  xreplace(eq, reps)
}

.ode_is_seq <- function(x) is.list(x) && !is_sym(x)

#' Check solutions of ODEs
#'
#' `checkodesol(ode, sol, func)` substitutes the solution `sol` into the
#' ODE `ode` and returns `(TRUE, 0)` when it is a solution, or `(FALSE,
#' residual)`, as SymPy's `checkodesol()`. A list of solutions gives a list
#' of results, and a list of equations is checked as a system with
#' [checksysodesol()].
#'
#' @param ode The ODE, an expression or an equation.
#' @param sol A solution, as an equation or the value of the function, or a
#'   list of them.
#' @param func The unknown function; detected when `NULL`.
#' @param order The order of the ODE, or `"auto"`.
#' @param solve_for_func Whether to try solving an implicit solution for
#'   the function first.
#' @return A tuple `(TRUE, 0)` or `(FALSE, residual)`, or a list of them.
#' @examples
#' x <- Symbol("x"); f <- Function("f")
#' checkodesol(Derivative(f(x), x) - f(x), Eq(f(x), exp(x)), f(x))
#' @export
checkodesol <- function(ode, sol, func = NULL, order = "auto", solve_for_func = TRUE) {
  if (.ode_is_seq(ode)) return(checksysodesol(ode, sol, func = func))
  ode <- as_sym(ode)
  if (!is.null(func)) func <- as_sym(func)
  if (ode$type != "Equality") ode <- Eq(ode, S_Zero)
  if (!.ode_is_seq(sol)) sol <- as_sym(sol)
  if (is.null(func)) {
    func <- tryCatch(.ode_preprocess(ode$args[[1L]])[[2L]], ValueError = function(e) NULL)
    if (is.null(func)) {
      fs <- list()
      for (s in (if (.ode_is_seq(sol)) sol else list(sol))) fs <- c(fs, .ode_undef_atoms(s))
      if (length(fs)) fs <- fs[!duplicated(vapply(fs, .k, ""))]
      if (length(fs) != 1L) .ode_stop("ValueError", "must pass func arg to checkodesol for this case.")
      func <- fs[[1L]]
    }
  }
  if (!.is_undefined_application(func) || length(func$args) != 1L)
    .ode_stop("ValueError", paste0("func must be a function of one variable, not ", .str(func)))
  if (.ode_is_seq(sol)) return(lapply(sol, function(i) checkodesol(ode, i, order = order, solve_for_func = solve_for_func)))
  if (sol$type != "Equality") sol <- Eq(func, sol)
  else if (identical(sol$args[[2L]]$key, func$key)) sol <- Eq(sol$args[[2L]], sol$args[[1L]])
  if (identical(order, "auto")) order <- .ode_order(ode, func)
  solved <- identical(sol$args[[1L]]$key, func$key) && !.has1(sol$args[[2L]], func)
  if (solve_for_func && !solved) {
    rhs <- .sp_solve(sol, list(func), list())
    if (length(rhs)) {
      eqs <- lapply(rhs, function(t) Eq(func, t))
      if (length(rhs) == 1L) eqs <- eqs[[1L]]
      return(checkodesol(ode, eqs, order = order, solve_for_func = FALSE))
    }
  }
  x <- func$args[[1L]]
  if (.has_type(sol, "Order")) {
    Oterm <- getO(sol$args[[2L]])
    solrhs <- removeO(sol$args[[2L]])
    sorder <- Oterm$args[[1L]]$args[[2L]]
    odesubs <- expand(doit(subs(.sub(ode$args[[1L]], ode$args[[2L]]), func, solrhs)))
    neworder <- O(.pow(x, .sub(sorder, .int(order))))
    odesubs <- .add(list(odesubs, neworder))
    residual <- removeO(odesubs)
    return(.pytuple(list(is_zero_num(residual), residual)))
  }
  s <- TRUE
  testnum <- 0L
  repeat {
    if (!.ode_true(s)) break
    if (testnum == 0L) {
      ode_diff <- .sub(ode$args[[1L]], ode$args[[2L]])
      if (identical(sol$args[[1L]]$key, func$key)) {
        s <- .ode_sub_func_doit(ode_diff, func, sol$args[[2L]])
        s <- besselsimp(s)
      } else {
        testnum <- testnum + 1L
        next
      }
      ss <- simplify(rewrite(s, "exp"))
      s <- if (.ode_true(ss)) expand(ss, force = TRUE) else S_Zero
      testnum <- testnum + 1L
    } else if (testnum == 1L) {
      s <- simplify(.add(list(trigsimp(.sub(.ode_d(sol$args[[1L]], x, order), .ode_d(sol$args[[2L]], x, order))),
                              .neg(trigsimp(ode$args[[1L]])), trigsimp(ode$args[[2L]]))))
      testnum <- testnum + 1L
    } else if (testnum == 2L) {
      if (identical(sol$args[[1L]]$key, func$key) && !.has1(sol$args[[2L]], func)) diffsols <- list(`0` = sol$args[[2L]])
      else if (identical(sol$args[[2L]]$key, func$key) && !.has1(sol$args[[1L]], func)) diffsols <- list(`0` = sol$args[[1L]])
      else diffsols <- list()
      solx <- .sub(sol$args[[1L]], sol$args[[2L]])
      failed <- FALSE
      for (i in seq_len(order)) {
        if (i == 1L) {
          ds <- .diff1(solx, x)
          sdf <- tryCatch(.sp_solve(ds, list(.ode_d(func, x, i)), list()), NotImplementedError = function(e) list())
          if (!length(sdf)) {
            testnum <- testnum + 1L
            failed <- TRUE
            break
          }
          diffsols[["1"]] <- sdf[[1L]]
        } else diffsols[[as.character(i)]] <- .diff1(diffsols[[as.character(i - 1L)]], x)
      }
      if (failed) next
      lhs <- ode$args[[1L]]
      rhs <- ode$args[[2L]]
      for (i in seq(order, 0L)) {
        if (i == 0L && is.null(diffsols[["0"]])) break
        lhs <- .ode_sub_func_doit(lhs, .ode_d(func, x, i), diffsols[[as.character(i)]])
        rhs <- .ode_sub_func_doit(rhs, .ode_d(func, x, i), diffsols[[as.character(i)]])
        ob <- simplify(Eq(lhs, rhs))
        if (ob$type == "BooleanAtom") {
          if (ob$name == "true") {
            lhs <- S_Zero
            rhs <- S_Zero
          }
        } else {
          lhs <- ob$args[[1L]]
          rhs <- ob$args[[2L]]
        }
      }
      num <- trigsimp(.as_numer_denom(.sub(lhs, rhs))[[1L]])
      fd <- Dummy("func")
      num <- subs(num, func, fd)
      pe <- posify(num)
      s <- xreplace(xreplace(simplify(pe$expr), pe$reps), list(list(fd, func)))
      testnum <- testnum + 1L
    } else break
  }
  if (!.ode_true(s)) return(.pytuple(list(TRUE, if (isTRUE(s) || isFALSE(s)) S_Zero else s)))
  if (isTRUE(s)) .ode_nie(paste0("Unable to test if ", .str(sol), " is a solution to ", .str(ode), "."))
  .pytuple(list(FALSE, s))
}

#' Check solutions of systems of ODEs
#'
#' `checksysodesol(eqs, sols)` substitutes the solutions `sols` (a list of
#' equations, one per function) into the system `eqs`, as SymPy's
#' `checksysodesol()`, and returns `(TRUE, residuals)` when all residuals
#' vanish, `(FALSE, residuals)` otherwise.
#'
#' @param eqs A list of ODEs.
#' @param sols A list of solutions.
#' @param func The unknown functions; detected when `NULL`.
#' @return A tuple of a logical and the list of residuals.
#' @examples
#' t <- Symbol("t"); f <- Function("f"); g <- Function("g")
#' checksysodesol(list(Eq(Derivative(f(t), t), g(t)), Eq(Derivative(g(t), t), -f(t))),
#'                list(Eq(f(t), sin(t)), Eq(g(t), cos(t))))
#' @export
checksysodesol <- function(eqs, sols, func = NULL) {
  eqs <- if (.ode_is_seq(eqs)) lapply(eqs, as_sym) else list(as_sym(eqs))
  sols <- if (.ode_is_seq(sols)) lapply(sols, as_sym) else list(as_sym(sols))
  eqs <- lapply(eqs, .eq_to_expr)
  funcs <- list()
  if (is.null(func)) {
    for (eq in eqs) for (d in atoms(eq, "Derivative")) funcs <- c(funcs, .ode_undef_atoms(d))
    if (length(funcs)) funcs <- .sp_sorted(funcs[!duplicated(vapply(funcs, .k, ""))])
  } else funcs <- if (.ode_is_seq(func)) lapply(func, as_sym) else list(as_sym(func))
  if (!all(vapply(funcs, function(fn) .is_undefined_application(fn) && length(fn$args) == 1L, TRUE)) &&
      length(unique(vapply(funcs, function(fn) paste(vapply(fn$args, .k, ""), collapse = ","), ""))) != 1L)
    .ode_stop("ValueError", paste0("func must be a function of one variable, not ", .str(funcs[[length(funcs)]])))
  for (sol in sols) if (length(.ode_undef_atoms(sol)) != 1L) .ode_stop("ValueError", "solutions should have one function only")
  if (length(funcs) != length(unique(vapply(sols, function(s) s$args[[1L]]$key, ""))))
    .ode_stop("ValueError", "number of solutions provided does not match the number of equations")
  dk <- character()
  dv <- list()
  for (sol in sols) {
    fn <- .ode_undef_atoms(sol)[[1L]]
    if (identical(sol$args[[2L]]$key, fn$key)) sol <- Eq(sol$args[[2L]], sol$args[[1L]])
    solved <- identical(sol$args[[1L]]$key, fn$key) && !.has1(sol$args[[2L]], fn)
    if (!solved) {
      rhs <- .sp_solve(sol, list(fn), list())
      if (!length(rhs)) .ode_nie("")
    } else rhs <- sol$args[[2L]]
    dk <- c(dk, fn$key)
    dv[[length(dv) + 1L]] <- rhs
  }
  checkeq <- list()
  for (eq in eqs) {
    for (fn in funcs) eq <- .ode_sub_func_doit(eq, fn, dv[[match(fn$key, dk)]])
    ss <- simplify(eq)
    if (!is_zero_num(ss)) {
      eq <- expand(ss, force = TRUE)
      if (!is_zero_num(eq)) eq <- simplify(sqrtdenest(eq))
    } else eq <- S_Zero
    checkeq[[length(checkeq) + 1L]] <- eq
  }
  ok <- length(unique(vapply(checkeq, .k, ""))) == 1L && is_zero_num(checkeq[[1L]])
  .pytuple(list(ok, checkeq))
}
