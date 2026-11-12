# Expand hypergeometric functions

`hyperexpand(f)` rewrites the generalized hypergeometric functions
([`hyper()`](https://adamrsophie.github.io/symcompR/reference/special-functions-extra.md))
and Meijer G functions
([`meijerg()`](https://adamrsophie.github.io/symcompR/reference/special-functions-extra.md))
in `f` in terms of named special functions, as SymPy's `hyperexpand()`.
Functions it cannot express are left unchanged.

## Usage

``` r
hyperexpand(f, allow_hyper = FALSE, rewrite = "default", place = NULL)
```

## Arguments

- f:

  An expression.

- allow_hyper:

  Allow the result to contain simpler
  [`hyper()`](https://adamrsophie.github.io/symcompR/reference/special-functions-extra.md)
  terms when a Meijer G function has no closed form.

- rewrite:

  The branch rewriting of the result: `"default"`, `"nonrepsmall"`,
  `"nonrep"`, or `NULL` for none.

- place:

  Where the expansion should be valid for a Meijer G function with two
  representations: `0` or `zoo`.

## Value

An expression.

## Examples

``` r
z <- Symbol("z")
hyperexpand(hyper(list(1, 1), list(2), -z) * z)
#> log(z + 1)
hyperexpand(hyper(list(), list(Rational(1, 2)), -z^2 / 4))
#> cos(z)
```
