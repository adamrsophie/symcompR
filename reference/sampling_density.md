# Sampling density

`sampling_density(expr)` of `sympy.stats`: the frequencies of sampled
values. SymPy samples through scipy, numpy or pymc; this function
follows SymPy without those libraries and raises SymPy's error. Use
[`rv_sample()`](https://adamrsophie.github.io/symcompR/reference/rv_sample.md)
for samples drawn with R's random number generator.

## Usage

``` r
sampling_density(
  expr,
  given_condition = NULL,
  library = "scipy",
  numsamples = 1,
  seed = NULL,
  ...
)
```

## Arguments

- expr:

  A random expression.

- given_condition:

  An optional condition.

- library:

  The sampling library SymPy would use.

- numsamples:

  Number of samples.

- seed:

  Ignored.

- ...:

  Ignored.

## Value

Does not return; raises a `ValueError` condition.

## Examples

``` r
r <- tryCatch(sampling_density(Normal("X", 0, 1)), error = conditionMessage)
r
#> [1] "Failed to import scipy"
```
