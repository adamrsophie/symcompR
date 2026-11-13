# Text plot

`textplot(expr, a, b, W, H)` prints an ASCII plot of a function of one
variable on `[a, b]`, with the layout of SymPy's `textplot()`.

## Usage

``` r
textplot(expr, a, b, W = 55, H = 21)
```

## Arguments

- expr:

  An expression in one variable.

- a, b:

  Range; numbers or numeric expressions such as `2*pi`.

- W, H:

  Width and height in characters.

## Value

The character rows, invisibly.

## Examples

``` r
textplot(Symbol("x")^2, -2, 2, W = 20, H = 6)
#>       4 |\                  /
#>         | \                / 
#>     2.4 |--\--------------/--
#>         |   ..          ..   
#>         |     ..      ..     
#>       0 |____________________
#>          -2        0.105263 2
```
