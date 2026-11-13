# Motion parameters from kinematic conditions

Acceleration, velocity and position of a particle in `frame` from one of
them and boundary conditions, as SymPy's `get_motion_params()`.

## Usage

``` r
get_motion_params(frame, ...)
```

## Arguments

- frame:

  A
  [`ReferenceFrame()`](https://adamrsophie.github.io/symcompR/reference/ReferenceFrame.md).

- ...:

  Any of `acceleration`, `velocity`, `position`, `timevalue`,
  `timevalue1` and `timevalue2`.

## Value

A list of three vectors: acceleration, velocity and position.

## Examples

``` r
N <- ReferenceFrame("N")
get_motion_params(N, position = Symbol("t")^2 * N$x)
#> (2*N.x, 2*t*N.x, t^2*N.x) 
```
