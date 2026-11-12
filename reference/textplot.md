# Text plot

`textplot(expr, a, b, W, H)` prints an ASCII plot of a function of one
variable on `[a, b]`, like SymPy's `textplot()`.

## Usage

``` r
textplot(expr, a, b, W = 55, H = 18)
```

## Arguments

- expr:

  An expression in one variable.

- a, b:

  Range.

- W, H:

  Width and height in characters.

## Value

The character rows, invisibly.
