# Partial velocities

The partial velocities of each velocity vector with respect to the
generalized speeds, as SymPy's `partial_velocity()`.

## Usage

``` r
partial_velocity(vel_vecs, gen_speeds, frame)
```

## Arguments

- vel_vecs:

  A list of velocity vectors.

- gen_speeds:

  A list of generalized speeds.

- frame:

  The
  [`ReferenceFrame()`](https://adamrsophie.github.io/symcompR/reference/ReferenceFrame.md)
  the velocities are in.

## Value

A list (one per velocity) of lists of vectors.

## Examples

``` r
N <- ReferenceFrame("N")
u <- dynamicsymbols("u1 u2")
partial_velocity(list(u[[1]] * N$x + u[[2]] * N$y), u, N)
#> [[1]]
#> [[1]][[1]]
#> N.x
#> 
#> [[1]][[2]]
#> N.y
#> 
#> 
```
