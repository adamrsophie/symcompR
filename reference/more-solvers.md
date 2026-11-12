# More equation solvers

`decompogen(f, x)` splits `f` into a list of functions whose composition
is `f`. `homogeneous_order(eq, ...)` is the order of homogeneity of `eq`
in the given symbols or function applications, or `NULL`.
`solve_linear_system_LU(M, syms)` solves an augmented matrix with LU
decomposition. `solve_triangulated(polys, gens)` solves a zero
dimensional polynomial system through a lex Groebner basis and returns
its rational solutions as lists ordered like `gens`.
`factor_system(eqs, gens)` splits a polynomial system into irreducible
subsystems: [`list()`](https://rdrr.io/r/base/list.html) means no
solution and `list(list())` means every value solves it.

## Usage

``` r
decompogen(f, x)

homogeneous_order(eq, ...)

solve_linear_system_LU(M, syms)

solve_triangulated(polys, gens)

factor_system(eqs, gens = NULL)
```

## Arguments

- f, eq:

  An expression.

- x:

  A symbol.

- ...:

  Symbols or function applications.

- M:

  An augmented matrix.

- syms:

  A list of symbols.

- polys, eqs:

  Lists of polynomial expressions, each equal to zero.

- gens:

  Generators; all free symbols by default for `factor_system()`.

## Value

A list, a number or `NULL`.

## Examples

``` r
x <- Symbol("x"); y <- Symbol("y")
decompogen(sin(x)^2 + sin(x) + 1, x)
#> [[1]]
#> x^2 + x + 1
#> 
#> [[2]]
#> sin(x)
#> 
homogeneous_order(x^2 + x*y, x, y)
#> 2
solve_triangulated(list(x^2 + y^2 - 1, x - y + 1), list(x, y))
#> [[1]]
#> [[1]][[1]]
#> -1
#> 
#> [[1]][[2]]
#> 0
#> 
#> 
#> [[2]]
#> [[2]][[1]]
#> 0
#> 
#> [[2]][[2]]
#> 1
#> 
#> 
factor_system(list(x^2 - 1, y - 1))
#> [[1]]
#> [[1]][[1]]
#> x + 1
#> 
#> [[1]][[2]]
#> y - 1
#> 
#> 
#> [[2]]
#> [[2]][[1]]
#> x - 1
#> 
#> [[2]][[2]]
#> y - 1
#> 
#> 
```
