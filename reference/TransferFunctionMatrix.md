# Transfer function matrices and MIMO interconnections

`TransferFunctionMatrix(arg)` is the transfer function matrix of a MIMO
system, built from a list of rows of SISO systems.
`TransferFunctionMatrix$from_Matrix(matrix, var)` builds one from a
matrix of rational expressions. `MIMOSeries()`, `MIMOParallel()` and
`MIMOFeedback()` are the interconnections of MIMO systems. A transfer
function matrix has the properties `$var`, `$num_inputs`, `$num_outputs`
and `$shape` and the methods `$transpose()`, `$elem_poles()`,
`$elem_zeros()`, `$eval_frequency()`, `$simplify()` and `$expand()`;
`tfm[i, j]` is an element or a submatrix. After
`sympy.physics.control.lti`.

## Usage

``` r
TransferFunctionMatrix(arg)

# S3 method for class 'TransferFunctionMatrix'
x[i, j, drop = TRUE]

MIMOSeries(..., evaluate = FALSE)

MIMOParallel(..., evaluate = FALSE)

MIMOFeedback(sys1, sys2, sign = -1L)
```

## Arguments

- arg:

  A list of rows of SISO systems.

- x:

  A transfer function matrix.

- i, j:

  The row and column indices.

- drop:

  Whether a single element is returned as a SISO system.

- ...:

  The systems.

- evaluate:

  Whether to evaluate the interconnection.

- sys1, sys2:

  The systems of the feed forward and feedback paths.

- sign:

  The sign of the feedback.

## Value

A `TransferFunctionMatrix`, `MIMOSeries`, `MIMOParallel` or
`MIMOFeedback` object.

## Examples

``` r
s <- Symbol("s")
G1 <- TransferFunction(1, s + 1, s); G2 <- TransferFunction(s, s + 2, s)
M <- TransferFunctionMatrix(list(list(G1, G2), list(G2, G1)))
M$shape
#> (2, 2) 
MIMOSeries(M, M)$doit()
#> TransferFunctionMatrix(((TransferFunction(s^2*(s + 1)^2 + (s + 2)^2, (s + 1)^2*(s + 2)^2, s), TransferFunction(2*s*(s + 1)*(s + 2), (s + 1)^2*(s + 2)^2, s)), (TransferFunction(2*s*(s + 1)*(s + 2), (s + 1)^2*(s + 2)^2, s), TransferFunction(s^2*(s + 1)^2 + (s + 2)^2, (s + 1)^2*(s + 2)^2, s))))
```
