# Integration over polytopes

`polytope_integrate(poly, expr)` integrates a polynomial in `x` and `y`
over a polygon, or in `x`, `y` and `z` over a polyhedron, as
`sympy.integrals.intpoly.polytope_integrate()`. A polygon is a
[`Polygon()`](https://adamrsophie.github.io/symcompR/reference/Point.md)
(the vertex representation, whose float coordinates are made rational)
or a list of hyperplanes `list(a, b)` meaning `a . (x, y) <= b`; a
polyhedron is a list whose first entry is the list of vertices and whose
other entries are the faces, each a vector of zero based vertex indices,
as in SymPy. The vertices are expected in clockwise order unless
`clockwise = TRUE` sorts them. With `max_degree`, `expr` is a list of
polynomials and the result is a named list of the integrals of those of
degree at most `max_degree` (of every monomial up to that degree when
`expr` is `NULL`), named by the polynomials.

## Usage

``` r
polytope_integrate(poly, expr = NULL, clockwise = FALSE, max_degree = NULL)

intpoly_strip(monom)

main_integrate3d(expr, facets, vertices, hp_params, max_degree = NULL)

main_integrate(expr, facets, hp_params, max_degree = NULL)

polygon_integrate(facet, hp_param, index, facets, vertices, expr, degree)

distance_to_side(point, line_seg, A)

lineseg_integrate(polygon, index, line_seg, expr, degree)

integration_reduction(facets, index, a, b, expr, dims, degree)

left_integral2D(m, index, facets, x0, expr, gens)

integration_reduction_dynamic(
  facets,
  index,
  a,
  b,
  expr,
  degree,
  dims,
  x_index,
  y_index,
  max_index,
  x0,
  monomial_values,
  monom_index,
  vertices = NULL,
  hp_param = NULL
)

left_integral3D(facets, index, expr, vertices, hp_param, degree)

gradient_terms(binomial_power = 0, no_of_gens = 2)

hyperplane_parameters(poly, vertices = NULL)

intpoly_cross_product(v1, v2, v3)

best_origin(a, b, lineseg, expr)

intpoly_decompose(expr, separate = FALSE)

point_sort(poly, normal = NULL, clockwise = TRUE)

intpoly_norm(point)

intpoly_intersection(geom_1, geom_2, intersection_type)

is_vertex(ent)
```

## Arguments

- poly:

  A polygon, a list of hyperplanes or a polyhedron.

- expr:

  A polynomial, or a list of them.

- clockwise:

  Sort the vertices of a polygon clockwise first.

- max_degree:

  The maximum degree of the polynomials.

- monom:

  A monomial.

- facets:

  The facets: segments of a polygon, or faces of a polyhedron as zero
  based vertex indices.

- vertices:

  The vertices of a polyhedron.

- hp_params:

  The hyperplane parameters of the facets.

- facet:

  A face, as zero based vertex indices.

- hp_param:

  The hyperplane parameters of the face.

- index:

  The zero based index of the facet.

- degree:

  The degree of the polynomial.

- point:

  A point.

- line_seg:

  A segment, as a list of two points.

- A:

  The normal of the face.

- polygon:

  The face the segment belongs to (unused, as in SymPy).

- a, b:

  The hyperplane `a . x = b` of the facet.

- dims:

  The symbols of the dimensions.

- m:

  The number of facets.

- x0:

  The origin of the facet.

- gens:

  The symbols.

- x_index, y_index, max_index:

  Indices of the monomial, as in SymPy.

- monomial_values:

  The table of monomials with their values.

- monom_index:

  The zero based index of the monomial in that table.

- binomial_power:

  The maximum degree.

- no_of_gens:

  The number of generators, 2 or 3.

- v1, v2, v3:

  Three points.

- lineseg:

  A segment.

- separate:

  Return the set of terms instead.

- normal:

  The normal of the plane of three dimensional points.

- geom_1, geom_2:

  Two segments or two hyperplanes.

- intersection_type:

  `"segment2D"`, `"segment3D"` or `"plane2D"`.

- ent:

  An object.

## Value

An expression or a named list.

## Details

The helpers follow the SymPy module: `intpoly_decompose(expr, separate)`
splits a polynomial by degree (into its terms with `separate = TRUE`),
`hyperplane_parameters()` gives the facet hyperplanes, `best_origin()`,
`point_sort()`, `intpoly_norm()`, `intpoly_intersection()`,
`is_vertex()`, `gradient_terms()`, `intpoly_cross_product()` and
`intpoly_strip()` are the geometric and algebraic steps, and
`main_integrate()`, `main_integrate3d()`, `polygon_integrate()`,
`lineseg_integrate()`, `distance_to_side()`, `integration_reduction()`,
`integration_reduction_dynamic()`, `left_integral2D()` and
`left_integral3D()` carry out the reduction to lower dimensions.

## Examples

``` r
x <- Symbol("x"); y <- Symbol("y")
polytope_integrate(Polygon(Point(0, 0), Point(0, 1), Point(1, 1), Point(1, 0)), x * y)
#> 1/4
polytope_integrate(list(list(c(-1, 0), 0), list(c(1, 2), 4), list(c(0, -1), 0)), 1)
#> 4
cube <- list(list(c(0, 0, 0), c(0, 0, 5), c(0, 5, 0), c(0, 5, 5), c(5, 0, 0), c(5, 0, 5),
                  c(5, 5, 0), c(5, 5, 5)),
             c(3, 7, 6, 2), c(1, 5, 7, 3), c(5, 4, 6, 7), c(0, 4, 5, 1),
             c(2, 0, 1, 3), c(2, 6, 4, 0))
polytope_integrate(cube, x^2 + y^2 + x * y + Symbol("z")^2)
#> 15625/4
intpoly_decompose(x^2 + 3 * y + 7)
#> $`0`
#> 7
#> 
#> $`2`
#> x^2
#> 
#> $`1`
#> 3*y
#> 
```
