# Laplace helpers and singularity function integration

`laplace_correspondence(f, fdict)` replaces
`LaplaceTransform(f(t), t, s)` by `F(s)` (and the inverse transform
back) for the pairs `list(list(f, F))`.
`laplace_initial_conds(f, t, fdict)` substitutes initial values
`list(list(f, c(f0, f1, ...)))` for `f(0)` and the derivatives of `f` at
0. `singularityintegrate(f, x)` integrates singularity functions, or
returns `NULL`.

## Usage

``` r
laplace_correspondence(f, fdict)

laplace_initial_conds(f, t, fdict)

singularityintegrate(f, x)
```

## Arguments

- f:

  An expression.

- fdict:

  Pairs of functions, or of a function and its initial values.

- t, x:

  Variables.

## Value

An expression or `NULL`.

## Examples

``` r
x <- Symbol("x"); a <- Symbol("a")
singularityintegrate(SingularityFunction(x, a, 3), x)
#> SingularityFunction(x, a, 4)/4
```
