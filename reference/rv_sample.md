# Draw random samples

`rv_sample(expr, condition, size)` is SymPy's
[`sample()`](https://rdrr.io/r/base/sample.html): one realization of a
random expression, or a vector of `size` of them. The draws follow the
exact law of each distribution: R's own generators for the families
SymPy hands to SciPy (Normal, Exponential, Gamma, Cauchy, Poisson and
others), and inverse transform sampling of the cdf over the whole
support for the rest. Samples use R's random number generator, so they
differ from SciPy's draws for the same seed; call
[`set.seed()`](https://rdrr.io/r/base/Random.html) (or pass `seed`) for
reproducible draws.

## Usage

``` r
rv_sample(expr, condition = NULL, size = 1, seed = NULL, ...)
```

## Arguments

- expr:

  A random expression.

- condition:

  An optional condition the draws must satisfy; a number here is taken
  as `size`.

- size:

  Number of draws.

- seed:

  An optional seed passed to
  [`set.seed()`](https://rdrr.io/r/base/Random.html).

- ...:

  Ignored.

## Value

A numeric vector.

## Examples

``` r
set.seed(1)
rv_sample(Die("D", 6), size = 5)
#> [1] 3 4 5 1 3
```
