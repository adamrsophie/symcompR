# Expand an expression

Distributes products over sums and expands powers of sums, like SymPy's
`expand()`. Hints switch individual rules: `mul` (distribute),
`multinomial` (powers of sums), `power_exp` (`x^(a+b) = x^a*x^b`),
`power_base` (`(x*y)^a = x^a*y^a` when valid), `log`
(`log(x*y) = log(x) + log(y)` when valid), `trig` (angle addition),
`complex` (split into real and imaginary parts) and `func` (special
function identities). `force = TRUE` ignores the validity assumptions of
`power_base` and `log`.

## Usage

``` r
expand(
  e,
  deep = TRUE,
  mul = TRUE,
  multinomial = TRUE,
  power_exp = TRUE,
  power_base = TRUE,
  log = TRUE,
  trig = FALSE,
  complex = FALSE,
  func = FALSE,
  force = FALSE,
  basic = TRUE,
  factor = FALSE,
  ...
)

expand_mul(e, deep = TRUE, force = FALSE)

expand_multinomial(e, deep = TRUE, force = FALSE)

expand_log(e, deep = TRUE, force = FALSE, factor = FALSE)

expand_power_exp(e, deep = TRUE, force = FALSE)

expand_power_base(e, deep = TRUE, force = FALSE)

expand_trig(e, deep = TRUE, force = FALSE)

expand_complex(e, deep = TRUE, force = FALSE)

expand_func(e, deep = TRUE, force = FALSE)
```

## Arguments

- e:

  A `sym` expression, matrix or list.

- deep:

  Expand inside function arguments too.

- mul, multinomial, power_exp, power_base, log, trig, complex, func:

  Hints.

- force:

  Apply `log` and `power_base` regardless of assumptions.

- basic:

  The `basic` hint, as in SymPy.

- factor:

  For the `log` hint: factor integer arguments of logarithms.

- ...:

  Ignored.

## Value

A `sym` expression.

## Details

The helpers `expand_mul()`, `expand_multinomial()`, `expand_log()`,
`expand_power_exp()`, `expand_power_base()`, `expand_trig()`,
`expand_complex()` and `expand_func()` apply one rule.

## Examples

``` r
x <- Symbol("x"); y <- Symbol("y")
expand((x + y)^3)
#> x^3 + 3*x^2*y + 3*x*y^2 + y^3
expand((x + 1) * (x - 1))
#> x^2 - 1
expand_trig(sin(x + y))
#> sin(x)*cos(y) + sin(y)*cos(x)
```
