# Rewriting of expressions for numerical evaluation, after
# sympy/codegen/rewriting.py.
#
# An optimization is an R function of an expression (class
# `c(<class>, "Optimization", "function")`) carrying the attributes
# `cost_function` and `priority`, as SymPy's `Optimization` objects. The
# predefined optimizations build their patterns when they are called, as the
# expression constructors are defined in files collated after this one.

.cg_optim_obj <- function(fn, cost_function = NULL, priority = 1, cls = "ReplaceOptim", ...) {
  structure(fn, class = c(cls, "Optimization", "function"),
            cost_function = cost_function, priority = priority, ...)
}

# The class of a node as SymPy's `expr.func`: the name of a function
# application, otherwise the node type.
.cg_head <- function(e) if (is_sym(e)) { if (e$type == "Function") e$name else e$type } else ""

# The head of the nodes an R function of symbols builds (`exp` gives "exp",
# `Pow` gives "Pow").
.cg_func_head <- function(f) {
  if (is.character(f)) return(f)
  d <- list(Dummy(), Dummy())
  for (n in 1:2) {
    r <- tryCatch(do.call(f, d[seq_len(n)]), error = function(err) NULL)
    if (is_sym(r)) return(.cg_head(r))
  }
  .cg_err("TypeError", "the function does not build a symbolic expression")
}

# Basic.count(query): the number of nodes in preorder satisfying `pred`.
.cg_count <- function(e, pred) {
  n <- as.integer(isTRUE(pred(e)))
  for (a in e$args) if (is_sym(a)) n <- n + .cg_count(a, pred)
  n
}

# Basic.replace(query, value): `query` is a pattern (a `sym` with wild
# symbols), a function name or node type, or an R predicate; `value` is an
# expression or an R function, applied bottom up to every node.
.cg_replace <- function(e, query, value, exact = NULL) {
  e <- as_sym(e)
  if (is_sym(query)) return(sym_replace(e, query, value, exact))
  if (is.character(query)) {
    q <- function(z) .cg_head(z) == query
    vf <- if (is.function(value)) function(z) do.call(value, z$args) else function(z) value
  } else {
    q <- query
    vf <- if (is.function(value)) value else function(z) value
  }
  .map_bottom_up(e, function(z) {
    if (!isTRUE(q(z))) return(z)
    v <- vf(z)
    if (is.null(v)) z else as_sym(v)
  })
}

#' Rewriting expressions for numerical evaluation
#'
#' `ReplaceOptim(query, value)` is an optimization that replaces the
#' subexpressions matching `query` by `value`, as SymPy's `ReplaceOptim`.
#' `sym_optimize(expr, optimizations)` applies optimizations in order of
#' decreasing priority, keeping the cheaper expression when an
#' optimization has a cost function (SymPy's `optimize()`).
#' `FuncMinusOneOptim(func, func_m_1)` rewrites sums such as `exp(x) - 1` as
#' `expm1(x)`, and `create_expand_pow_optimization(limit)` writes small
#' integer powers as repeated products.
#'
#' The predefined optimizations are `exp2_opt`, `sinc_opt1`, `sinc_opt2`,
#' `log2_opt`, `log2const_opt`, `logsumexp_2terms_opt`, `expm1_opt`,
#' `cosm1_opt`, `powm1_opt`, `log1p_opt`, `matinv_opt`, `logaddexp_opt` and
#' `logaddexp2_opt`, with the collections `sinc_opts`, `optims_c99`,
#' `optims_numpy` and `optims_scipy`. `matinv_opt` rewrites
#' the product of `Inverse(A)` and a column `x` as `MatrixSolve(A, x)` when
#' the assumptions context holds `Q$fullrank(A)` for the matrix symbol `A`
#' (see [assuming()] and [global_assumptions]).
#'
#' @param query A pattern with wild symbols, the name of a function or node
#'   type, or an R predicate of a node.
#' @param value The replacement: an expression (using the wild symbols of a
#'   pattern), or an R function of the matched wild symbols (for a pattern),
#'   of the arguments (for a name) or of the node (for a predicate).
#' @param cost_function `NULL`, or an R function giving the cost of an
#'   expression.
#' @param priority Optimizations of higher priority are applied first.
#' @param expr An expression.
#' @param optimizations A list of optimizations.
#' @param func,func_m_1 R functions building `f(x)` and `f(x) - 1` (such as
#'   `exp` and `sym_expm1`), or their names.
#' @param opportunistic Whether to rewrite whenever the constant term gets
#'   smaller, rather than only when it cancels.
#' @param limit The largest exponent to expand.
#' @param base_req An R predicate on the base of a power; the default accepts
#'   symbols.
#' @return An optimization (an R function of an expression), or for
#'   `sym_optimize()` an expression.
#' @examples
#' x <- Symbol("x")
#' sym_optimize(3*exp(x) - 3, list(expm1_opt))
#' sym_optimize(log(x)/log(Integer(2)) + log(x + 1), optims_c99)
#' ccode(sym_optimize(x^4 + x^2, list(create_expand_pow_optimization(4))))
#' @name codegen-rewriting
NULL

#' @rdname codegen-rewriting
#' @export
ReplaceOptim <- function(query, value, cost_function = NULL, priority = 1) {
  force(query); force(value)
  .cg_optim_obj(function(expr) .cg_replace(expr, query, value), cost_function, priority,
                query = query, value = value)
}

# Optimization.cheapest: the first argument of least cost
.cg_cheapest <- function(opt, ...) {
  args <- list(...)
  cf <- attr(opt, "cost_function")
  costs <- vapply(args, function(a) as.numeric(cf(a)), 0)
  args[[which.min(costs)]]
}

#' @rdname codegen-rewriting
#' @export
sym_optimize <- function(expr, optimizations) {
  if (inherits(optimizations, "Optimization")) optimizations <- list(optimizations)
  pr <- vapply(optimizations, function(o) as.numeric(attr(o, "priority")), 0)
  expr <- as_sym(expr)
  for (optim in optimizations[order(-pr)]) {
    new_expr <- optim(expr)
    expr <- if (is.null(attr(optim, "cost_function"))) new_expr else .cg_cheapest(optim, expr, new_expr)
  }
  expr
}

# Wild symbols of the predefined optimizations
.cg_rw_wild <- function(n) switch(n,
  u = Wild("u", properties = list(function(x) !.expr_is_number(x) && x$type != "Add")),
  v = Wild("v"),
  w = Wild("w"),
  n = Wild("n", properties = list(function(x) .expr_is_number(x))))

# A ReplaceOptim whose query and value `build()` returns when called
.cg_lazy_optim <- function(build, cost_function = NULL) {
  force(build)
  .cg_optim_obj(function(expr) {
    b <- build()
    .cg_replace(expr, b[[1L]], b[[2L]])
  }, cost_function, 1)
}

#' @rdname codegen-rewriting
#' @export
exp2_opt <- .cg_optim_obj(function(expr) .cg_replace(expr,
  function(p) p$type == "Pow" && identical(p$args[[1L]]$key, S_Two$key),
  function(p) exp2(p$args[[2L]])))

#' @rdname codegen-rewriting
#' @export
sinc_opt1 <- .cg_lazy_optim(function() {
  w <- .cg_rw_wild("w")
  list(sin(w) / w, sinc(w))
})

#' @rdname codegen-rewriting
#' @export
sinc_opt2 <- .cg_lazy_optim(function() {
  w <- .cg_rw_wild("w"); n <- .cg_rw_wild("n")
  list(sin(n * w) / w, n * sinc(n * w))
})

#' @rdname codegen-rewriting
#' @export
sinc_opts <- list(sinc_opt1, sinc_opt2)

# division and evaluation of transcendental functions are the expensive
# floating point operations
.cg_log2_cost <- function(expr) .cg_count(expr, function(e)
  (e$type == "Pow" && isTRUE(ask_fact(e$args[[2L]], "negative"))) ||
    (.cg_head(e) %in% c("log", "log2") && !.expr_is_number(e$args[[1L]])))

#' @rdname codegen-rewriting
#' @export
log2_opt <- .cg_lazy_optim(function() {
  v <- .cg_rw_wild("v"); w <- .cg_rw_wild("w")
  list(v * log(w) / log(S_Two), v * sym_log2(w))
}, cost_function = .cg_log2_cost)

#' @rdname codegen-rewriting
#' @export
log2const_opt <- .cg_lazy_optim(function() {
  w <- .cg_rw_wild("w")
  list(log(S_Two) * sym_log2(w), log(w))
})

#' @rdname codegen-rewriting
#' @export
logsumexp_2terms_opt <- .cg_optim_obj(function(expr) .cg_replace(expr,
  function(l) .cg_head(l) == "log" && l$args[[1L]]$type == "Add" && length(l$args[[1L]]$args) == 2L &&
    all(vapply(l$args[[1L]]$args, function(t) .cg_head(t) == "exp", TRUE)),
  function(l) {
    a <- lapply(l$args[[1L]]$args, function(e) e$args[[1L]])
    do.call(Max, a) + sym_log1p(exp(do.call(Min, a)))
  }))

# sift(args, pred, binary=True)
.cg_sift2 <- function(args, pred) {
  k <- vapply(args, function(a) isTRUE(pred(a)), TRUE)
  list(args[k], args[!k])
}
.cg_num_double <- function(z) tryCatch(as.numeric(N(z)), error = function(err) NA_real_)

#' @rdname codegen-rewriting
#' @export
FuncMinusOneOptim <- function(func, func_m_1, opportunistic = TRUE) {
  weight <- 10
  fh <- .cg_func_head(func)
  fm1h <- .cg_func_head(func_m_1)
  fm1 <- if (is.character(func_m_1)) function(...) fcall(func_m_1, lapply(list(...), as_sym)) else func_m_1
  is_f <- function(a) .cg_head(a) == fh
  replace_in_add <- function(e) {
    g <- .cg_sift2(e$args, .expr_is_number)
    numsum <- if (length(g[[1L]])) .add(g[[1L]]) else S_Zero
    g2 <- .cg_sift2(g[[2L]], function(a) .cg_count(a, is_f) > 0L)
    if (identical(numsum$key, S_Zero$key)) return(e)
    substituted <- list(); untouched <- list()
    for (with_func in g2[[1L]]) {
      coeff <- NULL
      if (with_func$type == "Mul") {
        fc <- .cg_sift2(with_func$args, is_f)
        if (length(fc[[1L]]) == 1L && length(fc[[2L]]) == 1L) { fn <- fc[[1L]][[1L]]; coeff <- fc[[2L]][[1L]] }
      } else if (is_f(with_func)) {
        fn <- with_func; coeff <- S_One
      }
      if (!is.null(coeff) && .expr_is_number(coeff) &&
          isTRUE(sign(.cg_num_double(coeff)) == -sign(.cg_num_double(numsum)))) {
        cs <- .add(list(coeff, numsum))
        do_substitute <- if (opportunistic) isTRUE(abs(.cg_num_double(cs)) < abs(.cg_num_double(numsum)))
                         else identical(cs$key, S_Zero$key)
        if (do_substitute) {
          numsum <- cs
          substituted[[length(substituted) + 1L]] <- .mul(list(coeff, as_sym(do.call(fm1, fn$args))))
          next
        }
      }
      untouched[[length(untouched) + 1L]] <- with_func
    }
    .add(c(list(numsum), substituted, untouched, g2[[2L]]))
  }
  cost <- function(expr) count_ops(expr) - weight * .cg_count(expr, function(e) .cg_head(e) == fm1h)
  query <- function(e) e$type == "Add"
  self <- NULL
  call <- function(expr) {
    expr <- as_sym(expr)
    alt1 <- .cg_replace(expr, query, replace_in_add)
    alt2 <- .cg_replace(sym_factor(expr), query, replace_in_add)
    .cg_cheapest(self, alt1, alt2)
  }
  self <- .cg_optim_obj(call, cost, 1, cls = c("FuncMinusOneOptim", "ReplaceOptim"),
                        query = query, value = replace_in_add, func = func, func_m_1 = func_m_1,
                        opportunistic = opportunistic)
  self
}

#' @rdname codegen-rewriting
#' @export
expm1_opt <- FuncMinusOneOptim("exp", "expm1")

#' @rdname codegen-rewriting
#' @export
cosm1_opt <- FuncMinusOneOptim("cos", "cosm1")

#' @rdname codegen-rewriting
#' @export
powm1_opt <- FuncMinusOneOptim("Pow", "powm1")

#' @rdname codegen-rewriting
#' @export
log1p_opt <- .cg_optim_obj(function(expr) .cg_replace(expr,
  function(e) .cg_head(e) == "log",
  function(l) {
    u <- .cg_rw_wild("u")
    l <- expand_log(.cg_replace(l, "log", function(arg, ...) log(sym_factor(arg))))
    .cg_replace(l, log(u + 1L), sym_log1p(u))
  }))

.cg_is_symbol <- function(b) is_sym(b) && b$type %in% c("Symbol", "Dummy", "Wild", "Indexed", "IndexedBase", "MatrixSymbol")

#' @rdname codegen-rewriting
#' @export
create_expand_pow_optimization <- function(limit, base_req = NULL) {
  if (is.null(base_req)) base_req <- .cg_is_symbol
  force(limit)
  ReplaceOptim(
    function(e) e$type == "Pow" && isTRUE(base_req(e$args[[1L]])) && e$args[[2L]]$type == "Integer" &&
      abs(.int_val(e$args[[2L]])) <= limit,
    function(p) {
      n <- .int_val(p$args[[2L]])
      prod <- function(k) UnevaluatedExpr(new_sym("Mul", rep(list(p$args[[1L]]), k)))
      if (n > 0) prod(n) else .pow(prod(-n), S_NegativeOne)
    })
}

# _matinv_predicate: an inverse of a full rank matrix symbol times a vector

.cg_matinv_predicate <- function(expr) {
  if (expr$type == "MatMul" && length(expr$args) == 2L) {
    left <- expr$args[[1L]]; right <- expr$args[[2L]]
    if (left$type == "Inverse" && identical(as_sym(.mshape(right)[[2L]])$key, S_One$key)) {
      inv_arg <- left$args[[1L]]
      if (inv_arg$type == "MatrixSymbol") return(isTRUE(ask(Q$fullrank(inv_arg))))
    }
  }
  FALSE
}

#' @rdname codegen-rewriting
#' @export
matinv_opt <- .cg_optim_obj(function(expr) .cg_replace(expr, .cg_matinv_predicate,
  function(e) MatrixSolve(e$args[[1L]]$args[[1L]], e$args[[2L]])))

#' @rdname codegen-rewriting
#' @export
logaddexp_opt <- .cg_lazy_optim(function() {
  v <- .cg_rw_wild("v"); w <- .cg_rw_wild("w")
  list(log(exp(v) + exp(w)), logaddexp(v, w))
})

#' @rdname codegen-rewriting
#' @export
logaddexp2_opt <- .cg_lazy_optim(function() {
  v <- .cg_rw_wild("v"); w <- .cg_rw_wild("w")
  list(log(.pow(S_Two, v) + .pow(S_Two, w)), logaddexp2(v, w) * log(S_Two))
})

#' @rdname codegen-rewriting
#' @export
optims_c99 <- list(expm1_opt, log1p_opt, exp2_opt, log2_opt, log2const_opt)

#' @rdname codegen-rewriting
#' @export
optims_numpy <- c(optims_c99, list(logaddexp_opt, logaddexp2_opt), sinc_opts)

#' @rdname codegen-rewriting
#' @export
optims_scipy <- list(cosm1_opt, powm1_opt)
