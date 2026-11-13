# Series, parallel and feedback interconnections

`Series(...)` is the series and `Parallel(...)` the parallel connection
of SISO systems; `Feedback(sys1, sys2, sign)` is the closed loop with
`sys1` in the feed forward path and `sys2` in the feedback path (`sign`
-1 for negative feedback). `$doit()` evaluates an interconnection to a
transfer function (or a state space model), `$to_expr()` gives it as an
expression and `$rewrite("TransferFunction")` rewrites it. After
`sympy.physics.control.lti`.

## Usage

``` r
Series(..., evaluate = FALSE)

Parallel(..., evaluate = FALSE)

Feedback(sys1, sys2 = NULL, sign = -1L)
```

## Arguments

- ...:

  The systems.

- evaluate:

  Whether to evaluate the interconnection.

- sys1, sys2:

  The systems of the feed forward and feedback paths.

- sign:

  The sign of the feedback.

## Value

A `Series`, `Parallel` or `Feedback` object.

## Examples

``` r
s <- Symbol("s")
G1 <- TransferFunction(1, s + 1, s); G2 <- TransferFunction(s, s + 2, s)
Series(G1, G2)$doit()
#> TransferFunction(s, (s + 1)*(s + 2), s)
Parallel(G1, G2)$doit()
#> TransferFunction(s*(s + 1) + s + 2, (s + 1)*(s + 2), s)
Feedback(G1, G2)$doit()
#> TransferFunction((s + 1)*(s + 2), (s + 1)*(s + (s + 1)*(s + 2)), s)
```
