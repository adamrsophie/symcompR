# Rewriting expressions for numerical evaluation

`ReplaceOptim(query, value)` is an optimization that replaces the
subexpressions matching `query` by `value`, as SymPy's `ReplaceOptim`.
`sym_optimize(expr, optimizations)` applies optimizations in order of
decreasing priority, keeping the cheaper expression when an optimization
has a cost function (SymPy's
[`optimize()`](https://rdrr.io/r/stats/optimize.html)).
`FuncMinusOneOptim(func, func_m_1)` rewrites sums such as `exp(x) - 1`
as `expm1(x)`, and `create_expand_pow_optimization(limit)` writes small
integer powers as repeated products.

## Usage

``` r
ReplaceOptim(query, value, cost_function = NULL, priority = 1)

sym_optimize(expr, optimizations)

exp2_opt(expr)

sinc_opt1(expr)

sinc_opt2(expr)

sinc_opts

log2_opt(expr)

log2const_opt(expr)

logsumexp_2terms_opt(expr)

FuncMinusOneOptim(func, func_m_1, opportunistic = TRUE)

expm1_opt(expr)

cosm1_opt(expr)

powm1_opt(expr)

log1p_opt(expr)

create_expand_pow_optimization(limit, base_req = NULL)

matinv_opt(expr)

logaddexp_opt(expr)

logaddexp2_opt(expr)

optims_c99

optims_numpy

optims_scipy
```

## Arguments

- query:

  A pattern with wild symbols, the name of a function or node type, or
  an R predicate of a node.

- value:

  The replacement: an expression (using the wild symbols of a pattern),
  or an R function of the matched wild symbols (for a pattern), of the
  arguments (for a name) or of the node (for a predicate).

- cost_function:

  `NULL`, or an R function giving the cost of an expression.

- priority:

  Optimizations of higher priority are applied first.

- expr:

  An expression.

- optimizations:

  A list of optimizations.

- func, func_m_1:

  R functions building `f(x)` and `f(x) - 1` (such as `exp` and
  `sym_expm1`), or their names.

- opportunistic:

  Whether to rewrite whenever the constant term gets smaller, rather
  than only when it cancels.

- limit:

  The largest exponent to expand.

- base_req:

  An R predicate on the base of a power; the default accepts symbols.

## Value

An optimization (an R function of an expression), or for
`sym_optimize()` an expression.

## Details

The predefined optimizations are `exp2_opt`, `sinc_opt1`, `sinc_opt2`,
`log2_opt`, `log2const_opt`, `logsumexp_2terms_opt`, `expm1_opt`,
`cosm1_opt`, `powm1_opt`, `log1p_opt`, `matinv_opt`, `logaddexp_opt` and
`logaddexp2_opt`, with the collections `sinc_opts`, `optims_c99`,
`optims_numpy` and `optims_scipy`. `matinv_opt` rewrites the product of
`Inverse(A)` and a column `x` as `MatrixSolve(A, x)` when the
assumptions context holds `Q$fullrank(A)` for the matrix symbol `A` (see
[`assuming()`](https://adamrsophie.github.io/symcompR/reference/assumptions-context.md)
and
[global_assumptions](https://adamrsophie.github.io/symcompR/reference/assumptions-context.md)).

## Examples

``` r
x <- Symbol("x")
sym_optimize(3*exp(x) - 3, list(expm1_opt))
#> 3*expm1(x)
sym_optimize(log(x)/log(Integer(2)) + log(x + 1), optims_c99)
#> log1p(x) + log2(x)
ccode(sym_optimize(x^4 + x^2, list(create_expand_pow_optimization(4))))
#> [1] "x*x + x*x*x*x"
```
