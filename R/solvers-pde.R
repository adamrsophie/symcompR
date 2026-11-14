# First order linear partial differential equations in two variables and
# separation of variables, after sympy/solvers/pde.py.

.PDE_HINTS <- c("1st_linear_constant_coeff_homogeneous", "1st_linear_constant_coeff",
                "1st_linear_constant_coeff_Integral", "1st_linear_variable_coeff")

# The unknown function: the undefined application inside the derivatives.
.pde_func <- function(e) {
  out <- list()
  walk <- function(z) {
    if (z$type == "Derivative" && .is_undefined_application(z$args[[1L]])) out[[z$args[[1L]]$key]] <<- z$args[[1L]]
    for (a in z$args) walk(a)
  }
  walk(e)
  if (length(out) != 1L) stop("The function cannot be automatically detected for ", .str(e), ".", call. = FALSE)
  out[[1L]]
}

# Total order of the highest derivative of fx.
.pde_order <- function(e, fx) {
  m <- 0L
  walk <- function(z) {
    if (z$type == "Derivative" && identical(z$args[[1L]]$key, fx$key))
      m <<- max(m, sum(vapply(z$args[-1L], function(t) as.integer(as.numeric(t$args[[2L]]$value)), 0L)))
    for (a in z$args) walk(a)
  }
  walk(e)
  m
}

.as_independent <- function(e, v) {
  if (e$type == "Add") {
    dep <- vapply(e$args, function(t) .has1(t, v), TRUE)
    return(list(.add(e$args[!dep]), .add(e$args[dep])))
  }
  facs <- if (e$type == "Mul") e$args else list(e)
  dep <- vapply(facs, function(t) .has1(t, v), TRUE)
  list(.mul(facs[!dep]), .mul(facs[dep]))
}

.has_any <- function(e, vs) any(vapply(vs, function(v) .has1(e, v), TRUE))

# b*fx + c*fy + d*f + e with b, c, d, e free of f and its derivatives.
.pde_linear_parts <- function(e, func, fx, fy) {
  D <- list(Dummy("p"), Dummy("q"), Dummy("u"))
  z <- expand(xreplace(xreplace(e, list(list(fx, D[[1L]]), list(fy, D[[2L]]))), list(list(func, D[[3L]]))))
  if (.has1(z, func)) return(NULL)
  cs <- lapply(D, function(d) .linear_coeff(z, d))
  rest <- expand(.sub(z, .add(Map(function(c, d) .mul(list(c, d)), cs, D))))
  if (any(vapply(c(cs, list(rest)), .has_any, TRUE, vs = D))) return(NULL)
  list(b = cs[[1L]], c = cs[[2L]], d = cs[[3L]], e = rest)
}

# SymPy _preprocess: f^p = 0 with p > 0 means f = 0.
.pde_prep <- function(e) {
  e <- doit(e)
  if (e$type == "Pow" && isTRUE(.facts(e$args[[2L]])[["positive"]])) e <- e$args[[1L]]
  e
}

#' First order partial differential equations
#'
#' `classify_pde(eq, func)` lists the solution methods that apply to `eq`,
#' a PDE in an unknown function of two variables such as `f(x, y)`.
#' `pdsolve(eq, func)` solves first order linear equations with constant or
#' variable coefficients; the arbitrary function is `F` unless `solvefun`
#' is given. `checkpdesol(eq, sol)` returns `list(ok, residual)`.
#' `pde_separate(eq, fun, sep, strategy)` separates variables with the
#' product (`"mul"`) or sum (`"add"`) ansatz `sep`, returning the two
#' separated sides or `NULL`.
#'
#' @param eq A PDE, as an equation or an expression equal to zero.
#' @param func The unknown function application, found when `NULL`.
#' @param dict Return a named list with the order, the default hint and the
#'   matches.
#' @param hint A method name or `"default"`.
#' @param solvefun The arbitrary function.
#' @param sol A solution, as an equation or the right side.
#' @param solve_for_func Solve `sol` for `func` first when needed.
#' @param fun The unknown function application, such as `u(x, t)`.
#' @param sep A list of function applications, such as `list(X(x), T(t))`.
#' @param strategy `"mul"` or `"add"`.
#' @return A character vector, a list or an equation.
#' @examples
#' x <- Symbol("x"); y <- Symbol("y"); f <- Function("f")
#' u <- f(x, y)
#' pdsolve(diff(u, x) + diff(u, y) + u)
#' classify_pde(x*diff(u, x) + diff(u, y))
#' @name pdsolve
NULL

#' @rdname pdsolve
#' @export
classify_pde <- function(eq, func = NULL, dict = FALSE) {
  e <- .pde_prep(.eq_to_expr(eq))
  func <- if (is.null(func)) .pde_func(e) else as_sym(func)
  if (length(func$args) != 2L)
    stop("Right now only partial differential equations of two variables are supported", call. = FALSE)
  x <- func$args[[1L]]
  y <- func$args[[2L]]
  fx <- .diff1(func, x)
  fy <- .diff1(func, y)
  order <- .pde_order(e, func)
  hints <- list(order = order)
  if (!order) return(if (dict) c(hints, list(default = NULL)) else character())
  e <- expand(e)
  reduced <- e
  if (e$type == "Add") {
    power <- NULL
    vars <- lapply(0:order, function(k) c(rep(list(x), order - k), rep(list(y), k)))
    for (vs in vars) {
      cf <- coeff(e, Reduce(function(a, v) .diff1(a, v), vs, func))
      if (identical(cf$key, S_One$key) || is_zero_num(cf)) next
      facs <- if (cf$type == "Mul") cf$args else list(cf)
      dep <- vapply(facs, function(t) .has1(t, func), TRUE)
      pw <- if (!any(dep)) 0 else if (sum(dep) == 1L) {
        t <- facs[dep][[1L]]
        if (identical(t$key, func$key)) 1
        else if (t$type == "Pow" && identical(t$args[[1L]]$key, func$key) && !.has_any(t$args[[2L]], list(x, y)))
          num_value_dbl(t$args[[2L]]) else NA
      } else NA
      if (!is.na(pw) && (is.null(power) || pw < power)) power <- pw
    }
    if (!is.null(power) && power != 0) {
      den <- .pow(func, as_sym(power))
      reduced <- .add(lapply(e$args, function(a) .div(a, den)))
    }
  }
  if (order == 1L) {
    reduced <- collect(reduced, func)
    m <- .pde_linear_parts(reduced, func, fx, fy)
    if (!is.null(m)) {
      const <- !any(vapply(m[c("b", "c", "d")], .has_any, TRUE, vs = list(x, y)))
      if (const) {
        if (is_zero_num(m$e)) hints[["1st_linear_constant_coeff_homogeneous"]] <- m
        else if (!is_zero_num(simplify(.add(list(.pow(m$b, S_Two), .pow(m$c, S_Two)))))) {
          hints[["1st_linear_constant_coeff"]] <- m
          hints[["1st_linear_constant_coeff_Integral"]] <- m
        }
      } else hints[["1st_linear_variable_coeff"]] <- m
    }
  }
  found <- .PDE_HINTS[.PDE_HINTS %in% names(hints)]
  if (dict) return(c(hints, list(default = if (length(found)) found[[1L]], ordered_hints = found)))
  found
}

#' @rdname pdsolve
#' @export
pdsolve <- function(eq, func = NULL, hint = "default", solvefun = NULL) {
  if (is.null(solvefun)) solvefun <- Function("F")
  e <- .pde_prep(.eq_to_expr(eq))
  func <- if (is.null(func)) .pde_func(e) else as_sym(func)
  cl <- classify_pde(e, func, dict = TRUE)
  if (identical(hint, "default")) hint <- cl$default
  if (is.null(hint) || is.null(cl[[hint]])) stop("psolve: Cannot solve ", .str(e), call. = FALSE)
  m <- cl[[hint]]
  x <- func$args[[1L]]
  y <- func$args[[2L]]
  b <- m$b; c <- m$c; d <- m$d; ee <- .neg(m$e)
  bc <- .add(list(.pow(b, S_Two), .pow(c, S_Two)))
  if (hint == "1st_linear_constant_coeff_homogeneous") {
    rhs <- .mul(list(fcall("exp", list(.mul(list(.neg(.div(d, bc)), .add(list(.mul(list(b, x)), .mul(list(c, y)))))))),
                     solvefun(.sub(.mul(list(c, x)), .mul(list(b, y))))))
    return(Eq(func, rhs, evaluate = FALSE))
  }
  if (hint %in% c("1st_linear_constant_coeff", "1st_linear_constant_coeff_Integral")) {
    xi <- Symbol("xi"); eta <- Symbol("eta")
    xs <- .div(.add(list(.mul(list(b, xi)), .mul(list(c, eta)))), bc)
    ys <- .div(.sub(.mul(list(c, xi)), .mul(list(b, eta))), bc)
    expterm <- fcall("exp", list(.mul(list(.neg(.div(d, bc)), xi))))
    # SymPy's Mul folds exp(a)*exp(b) into exp(a + b); powsimp does that here
    integrand <- powsimp(xreplace(.div(ee, expterm), list(list(x, xs), list(y, ys))))
    u <- .add(list(.mul(list(b, x)), .mul(list(c, y))))
    v <- .sub(.mul(list(c, x)), .mul(list(b, y)))
    if (hint == "1st_linear_constant_coeff_Integral") {
      genterm <- .div(Integral(integrand, list(xi, u)), bc)
      return(Eq(func, Subs(.mul(list(expterm, .add(list(solvefun(eta), genterm)))),
                           Tuple(eta, xi), Tuple(v, u)), evaluate = FALSE))
    }
    A <- integrate(integrand, xi)
    G <- .div(xreplace(A, list(list(xi, u), list(eta, v))), bc)
    rhs <- .mul(list(xreplace(expterm, list(list(xi, u))), .add(list(solvefun(v), G))))
    return(Eq(func, simplify(rhs), evaluate = FALSE))
  }
  .pde_variable_coeff(func, x, y, b, c, d, ee, solvefun)
}

# Generated constants: symbols of sol absent from the ODE and from x, y.
.pde_constants <- function(sol, ode, xy) {
  xk <- vapply(xy, .k, "")
  Filter(function(s) !.has1(ode, s) && !(s$key %in% xk), free_symbols(sol))
}
.first_sol <- function(sol) if (is.list(sol) && !is_sym(sol)) sol[[1L]] else sol

.pde_simplify_coeff <- function(sol, syms, fun, arg) {
  for (s in syms) sol <- subs(sol, s, fun(arg))
  simplify(subs(sol, Symbol("eta"), arg))
}

.pde_variable_coeff <- function(func, x, y, b, c, d, e, solvefun) {
  eta <- Symbol("eta")
  F1 <- Function(func$name)
  if (is_zero_num(d) && (is_zero_num(b) || is_zero_num(c))) {
    if (!is_zero_num(c)) return(Eq(func, .add(list(solvefun(x), integrate(.div(e, c), y))), evaluate = FALSE))
    if (!is_zero_num(b)) return(Eq(func, .add(list(solvefun(y), integrate(.div(e, b), x))), evaluate = FALSE))
  }
  one_var <- function(v, coef, other) {
    fv <- F1(v)
    ode <- .sub(.add(list(.mul(list(coef, .diff1(fv, v))), .mul(list(d, fv)))), e)
    sol <- .first_sol(dsolve(ode, fv))
    rhs <- sol$args[[2L]]
    Eq(func, .pde_simplify_coeff(rhs, .pde_constants(rhs, ode, list(x, y)), solvefun, other), evaluate = FALSE)
  }
  if (is_zero_num(c)) return(one_var(x, b, y))
  if (is_zero_num(b)) return(one_var(y, c, x))
  dm <- Function("d")
  dx <- dm(x)
  h <- subs(.div(c, b), y, dx)
  sol <- .first_sol(dsolve(.sub(.diff1(dx, x), h), dx))
  cs <- .pde_constants(sol$args[[2L]], h, list(x, y))
  if (length(cs) != 1L)
    stop("Cannot solve the partial differential equation due to inability of constantsimp", call. = FALSE)
  etat <- subs(solve(.sub(sol$args[[1L]], sol$args[[2L]]), cs[[1L]])[[1L]], dx, y)
  ysub <- solve(.sub(eta, etat), y)[[1L]]
  fx1 <- F1(x)
  deq <- subs(.sub(.add(list(.mul(list(b, .diff1(fx1, x))), .mul(list(d, fx1)))), e), y, ysub)
  final <- .first_sol(dsolve(deq, fx1))$args[[2L]]
  Eq(func, .pde_simplify_coeff(final, .pde_constants(final, deq, list(x, y)), solvefun, etat), evaluate = FALSE)
}

#' @rdname pdsolve
#' @export
checkpdesol <- function(eq, sol, func = NULL, solve_for_func = TRUE) {
  e <- .eq_to_expr(eq)
  if (is.null(func)) func <- .pde_func(doit(e))
  func <- as_sym(func)
  if (is.list(sol) && !is_sym(sol)) return(lapply(sol, function(s) checkpdesol(e, s, func, solve_for_func)))
  sol <- as_sym(sol)
  if (sol$type != "Equality") sol <- Eq(func, sol, evaluate = FALSE)
  else if (identical(sol$args[[2L]]$key, func$key)) sol <- Eq(sol$args[[2L]], sol$args[[1L]], evaluate = FALSE)
  solved <- identical(sol$args[[1L]]$key, func$key) && !.has1(sol$args[[2L]], func)
  if (solve_for_func && !solved) {
    r <- solve(.sub(sol$args[[1L]], sol$args[[2L]]), func)
    if (length(r) == 1L) return(checkpdesol(e, Eq(func, r[[1L]], evaluate = FALSE), func, FALSE))
    if (length(r)) return(lapply(r, function(t) checkpdesol(e, Eq(func, t, evaluate = FALSE), func, FALSE)))
  }
  if (identical(sol$args[[1L]]$key, func$key)) {
    s <- simplify(doit(subs(e, func, sol$args[[2L]])))
    return(list(is_zero_num(s), s))
  }
  stop("Unable to test if ", .str(sol), " is a solution to ", .str(e), ".", call. = FALSE)
}

#' @rdname pdsolve
#' @export
pde_separate <- function(eq, fun, sep, strategy = "mul") {
  if (!strategy %in% c("add", "mul")) stop(.sympy_condition("ValueError", paste0("Unknown strategy: ", strategy)))
  e <- .eq_to_expr(eq)
  fun <- as_sym(fun)
  sep <- lapply(sep, as_sym)
  orig <- vapply(fun$args, .k, "")
  sargs <- unlist(lapply(sep, `[[`, "args"), recursive = FALSE)
  sk <- vapply(sargs, .k, "")
  if (length(sk) != length(orig)) stop(.sympy_condition("ValueError", "Variable counts do not match"))
  if (anyDuplicated(sk)) stop(.sympy_condition("ValueError", "Duplicate substitution arguments detected"))
  if (!setequal(sk, orig)) stop(.sympy_condition("ValueError", "Arguments do not match"))
  functions <- if (strategy == "add") .add(sep) else .mul(sep)
  result <- doit(subs(e, fun, functions))
  if (strategy == "mul") {
    terms <- if (result$type == "Add") result$args else list(result)
    result <- .add(lapply(terms, function(t) .div(t, functions)))
  }
  .pde_split(result, sargs[[1L]], sargs[-1L])
}

#' @rdname pdsolve
#' @export
pde_separate_add <- function(eq, fun, sep) pde_separate(eq, fun, sep, strategy = "add")

#' @rdname pdsolve
#' @export
pde_separate_mul <- function(eq, fun, sep) pde_separate(eq, fun, sep, strategy = "mul")

# SymPy _separate: move the terms in dep to one side and the rest to the
# other, dividing out common factors.
.pde_split <- function(eq, dep, others) {
  terms <- if (eq$type == "Add") eq$args else list(eq)
  marked <- list()
  for (t in terms) {
    hit <- if (t$type == "Mul") any(vapply(t$args, function(i) i$type == "Derivative" && !.has_any(i, others), TRUE))
           else t$type == "Derivative" && !.has_any(t, others)
    if (hit) marked[[t$key]] <- t
  }
  div <- list()
  for (t in marked) {
    ai <- .as_independent(expand(t), dep)
    if (.has_any(ai[[2L]], others)) return(NULL)
    div[[ai[[1L]]$key]] <- ai[[1L]]
  }
  if (length(div)) terms <- lapply(terms, function(t) simplify(.add(lapply(div, function(i) .div(t, i)))))
  div <- list()
  lhs <- S_Zero
  rhs <- S_Zero
  for (t in terms) {
    if (!.has_any(t, others)) {
      lhs <- .add(list(lhs, t))
      next
    }
    ai <- .as_independent(expand(t), dep)
    if (.has_any(ai[[2L]], others)) return(NULL)
    div[[ai[[2L]]$key]] <- ai[[2L]]
    rhs <- .sub(rhs, expand(t))
  }
  fulldiv <- .add(unname(div))
  lhs <- expand(simplify(.div(lhs, fulldiv)))
  rhs <- expand(simplify(.div(rhs, fulldiv)))
  if (.has_any(lhs, others) || .has1(rhs, dep)) return(NULL)
  list(lhs, rhs)
}
