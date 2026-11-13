# Attribute access for classical mechanics objects

Frames, vectors, dyadics, points and the mechanics classes expose the
attributes, properties and methods of their SymPy counterparts through
`$`, so `N$x` is the unit vector of frame `N` and `v$dt(N)` is the time
derivative of `v` in `N`. Assigning with `$<-` runs a property setter
where SymPy has one.

## Usage

``` r
# S3 method for class 'physics_obj'
x$name

# S3 method for class 'physics_obj'
x$name <- value
```

## Arguments

- x:

  A physics object.

- name:

  Attribute name.

- value:

  New value.

## Value

The attribute value, or a function for a method.

## Examples

``` r
N <- ReferenceFrame("N")
N$x
#> N.x
N$name
#> [1] "N"
```
