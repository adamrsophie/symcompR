# Dimension systems

`DimensionSystem(base_dims, derived_dims, dimensional_dependencies)` is
SymPy's coherent set of dimensions: the base dimensions, the derived
ones and how the derived ones depend on the base ones. Dimensions may be
given as `Dimension` objects or names, and `dimensional_dependencies` as
a named list such as `list(velocity = list(length = 1, time = -1))`, a
list of `list(dimension, dependencies)` pairs, or a `units_dict`.
`get_dimensional_dependencies()` expresses a dimension through the base
dimensions as a `units_dict` (fields `k` and `v` hold the keys and the
powers, in Python's dict order); `equivalent_dims()` compares two
dimensions that way and `is_dimensionless()` tests for no dependencies.
`extend()` adds base and derived dimensions. The other functions are the
legacy helpers of SymPy (`list_can_dims`, `dim_can_vector`, ...).

## Usage

``` r
DimensionSystem(
  base_dims,
  derived_dims = list(),
  dimensional_dependencies = list()
)

base_dims(dim_sys)

derived_dims(dim_sys)

dimensional_dependencies(dim_sys)

get_dimensional_dependencies(dim_sys, name, mark_dimensionless = FALSE)

equivalent_dims(dim_sys, dim1, dim2)

extend(x, ...)

# S3 method for class 'DimensionSystem'
extend(x, new_base_dims, new_derived_dims = list(), new_dim_deps = NULL, ...)

is_dimensionless(dim_sys, dimension)

list_can_dims(dim_sys)

dim_can_vector(dim_sys, dim)

inv_can_transf_matrix(dim_sys)

can_transf_matrix(dim_sys)

dim_vector(dim_sys, dim)

print_dim_base(dim_sys, dim)

system_dim(x)

is_consistent(x)
```

## Arguments

- base_dims, new_base_dims:

  Base dimensions.

- derived_dims, new_derived_dims:

  Derived dimensions.

- dimensional_dependencies, new_dim_deps:

  Dependencies of the derived dimensions.

- dim_sys, x:

  A `DimensionSystem` (or, for `extend()`, `system_dim()` and
  `is_consistent()`, a `UnitSystem`).

- name, dim, dim1, dim2, dimension:

  Dimensions or dimension expressions.

- mark_dimensionless:

  Return `{Dimension(1): 1}` for a dimensionless input.

- ...:

  Further arguments of the `UnitSystem` method.

## Value

A `DimensionSystem`, `units_dict`, list, matrix, expression, integer or
`TRUE`/`FALSE`.

## Examples

``` r
u <- physics_units()
get_dimensional_dependencies(u$dimsys_SI, u$force)
#> {Dimension(mass): 1, Dimension(length): 1, Dimension(time): -2}
equivalent_dims(u$dimsys_SI, u$energy, u$force * u$length)
#> [1] TRUE
ds <- DimensionSystem(list(u$length, u$time), list(u$velocity))
base_dims(ds)
#> [[1]]
#> Dimension(length, L)
#> 
#> [[2]]
#> Dimension(time, T)
#> 
```
