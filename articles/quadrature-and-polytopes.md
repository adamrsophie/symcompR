# Quadrature and integration over polytopes

This article covers two numerical and geometric corners of SymPy’s
integration module: the nodes and weights of Gaussian quadrature rules,
and the exact integration of polynomials over polygons and polyhedra.

``` r

library(symcompR)
x <- Symbol("x"); y <- Symbol("y")
```

## Gaussian quadrature

The quadrature functions return the nodes and weights of Gaussian
quadrature rules of order `n`, to `n_digits` significant digits. The
helper `nw()` prints the nodes and the weights as character vectors:

``` r

nw <- function(r) lapply(setNames(r, c("nodes", "weights")),
                         function(v) vapply(v, sstr, ""))
nw(gauss_legendre(3, 5))
#> $nodes
#> [1] "-0.77460" "0"        "0.77460" 
#> 
#> $weights
#> [1] "0.55556" "0.88889" "0.55556"
nw(gauss_laguerre(3, 5))
#> $nodes
#> [1] "0.41577" "2.2943"  "6.2899" 
#> 
#> $weights
#> [1] "0.71109"  "0.27852"  "0.010389"
nw(gauss_hermite(3, 5))
#> $nodes
#> [1] "-1.2247" "0"       "1.2247" 
#> 
#> $weights
#> [1] "0.29541" "1.1816"  "0.29541"
nw(gauss_gen_laguerre(3, -Rational(1, 2), 5))
#> $nodes
#> [1] "0.19016" "1.7845"  "5.5253" 
#> 
#> $weights
#> [1] "1.4493"    "0.31413"   "0.0090600"
nw(gauss_chebyshev_t(3, 5))
#> $nodes
#> [1] "0.86602"  "0"        "-0.86602"
#> 
#> $weights
#> [1] "1.0472" "1.0472" "1.0472"
nw(gauss_chebyshev_u(3, 5))
#> $nodes
#> [1] "0.70711"  "0"        "-0.70711"
#> 
#> $weights
#> [1] "0.39270" "0.78540" "0.39270"
nw(gauss_jacobi(3, Rational(1, 2), -Rational(1, 2), 5))
#> $nodes
#> [1] "-0.90097" "-0.22252" "0.62349" 
#> 
#> $weights
#> [1] "1.7063"  "1.0973"  "0.33795"
nw(gauss_lobatto(3, 5))
#> $nodes
#> [1] "-1" "0"  "1" 
#> 
#> $weights
#> [1] "0.33333" "1.3333"  "0.33333"
```

## Integration over polytopes

[`polytope_integrate()`](https://adamrsophie.github.io/symcompR/reference/polytope_integrate.md)
integrates a polynomial over a polygon or a polyhedron. A polyhedron is
a list whose first element holds the vertices and whose other elements
are the faces, as zero based vertex indices.

``` r

polytope_integrate(Polygon(Point(0, 0), Point(0, 1), Point(1, 1), Point(1, 0)), x * y)
#> 1/4
cube <- list(list(c(0, 0, 0), c(0, 0, 5), c(0, 5, 0), c(0, 5, 5), c(5, 0, 0),
                  c(5, 0, 5), c(5, 5, 0), c(5, 5, 5)),
             c(2, 6, 7, 3), c(3, 7, 5, 1), c(7, 6, 4, 5), c(1, 5, 4, 0),
             c(3, 1, 0, 2), c(0, 4, 6, 2))
vertices <- cube[[1]]; faces <- cube[-1]
```

The algorithm reduces the integral over the region to integrals over its
facets with the generalized Stokes theorem. Its building blocks are
exported:
[`hyperplane_parameters()`](https://adamrsophie.github.io/symcompR/reference/polytope_integrate.md)
gives the hyperplane of each facet,
[`main_integrate()`](https://adamrsophie.github.io/symcompR/reference/polytope_integrate.md)
and
[`main_integrate3d()`](https://adamrsophie.github.io/symcompR/reference/polytope_integrate.md)
integrate over a polygon and a polyhedron, and the other functions are
the steps on single facets and line segments.

``` r

hp3 <- hyperplane_parameters(faces, vertices)
main_integrate3d(1, faces, vertices, hp3)
#> -125
polygon_integrate(faces[[1]], list(c(0, 1, 0), 5), 0, faces, vertices, 1, 0)
#> -25
left_integral3D(faces, 3, 1, vertices, list(c(0, -1, 0), -5), 0)
#> -50
tri <- Polygon(Point(0, 3), Point(5, 3), Point(1, 1))
facets <- list(Segment(Point(0, 3), Point(5, 3)), Segment(Point(5, 3), Point(1, 1)),
               Segment(Point(1, 1), Point(0, 3)))
hp <- hyperplane_parameters(tri)
sstr(do.call(Tuple, lapply(hp, function(h) Tuple(do.call(Tuple, h[[1]]), h[[2]]))))
#> [1] "(((0, 1), 3), ((1, -2), -1), ((-2, -1), -3))"
main_integrate(x^2 + y^2, facets, hp)
#> 325/6
x0 <- Point(0, 3)
integration_reduction(facets, 0, hp[[1]][[1]], hp[[1]][[2]], 1, list(x, y), 0)
#> 5
left_integral2D(3, 0, facets, x0, 1, list(x, y))
#> 5
mv <- list(list(0, 0, 0, 0), list(1, 0, 0, 5), list(y, 0, 1, 15), list(x, 1, 0, NULL))
integration_reduction_dynamic(facets, 0, hp[[1]][[1]], hp[[1]][[2]], x, 1, list(x, y),
                              1, 0, 1, x0, mv, 3)
#> 25/2
lineseg_integrate(list(c(0, 5, 0), c(5, 5, 0), c(5, 5, 5), c(0, 5, 5)), 0,
                  list(c(0, 5, 0), c(5, 5, 0)), 1, 0)
#> 5
distance_to_side(c(0, 0, 0), list(c(0, 0, 1), c(0, 1, 0)), c(1, 0, 0))
#> -sqrt(2)/2
```

Geometric and algebraic helpers:

``` r

sstr(do.call(Tuple, lapply(gradient_terms(2), function(r) do.call(Tuple, r))))
#> [1] "((1, 0, 0, 0), (y, 0, 1, 0), (y^2, 0, 2, 0), (x, 1, 0, 0), (x*y, 1, 1, 0), (x^2, 2, 0, 0))"
unlist(lapply(best_origin(c(2, 1), 3, Segment2D(Point(0, 3), Point(1, 1)), x^3 * y^7), sstr))
#> [1] "0" "3"
intpoly_norm(Point(2, 7))
#> sqrt(53)
unlist(lapply(intpoly_intersection(Segment2D(Point(1, 1), Point(3, 5)),
                                   Segment2D(Point(2, 0), Point(2, 5)), "segment2D"), sstr))
#> [1] "2" "3"
unlist(lapply(intpoly_intersection(list(c(-1, 0), 0), list(c(0, 1), 1), "plane2D"), sstr))
#> [1] "0" "1"
is_vertex(c(2, 3))
#> [1] TRUE
is_vertex(Point(2, 3))
#> [1] TRUE
point_sort(list(Point(0, 0), Point(1, 0), Point(1, 1)))
#> [[1]]
#> Point2D(1, 1) 
#> 
#> [[2]]
#> Point2D(1, 0) 
#> 
#> [[3]]
#> Point2D(0, 0)
intpoly_strip(3 * x^2)
#> [[1]]
#> 3
#> 
#> [[2]]
#> x^2
unlist(lapply(intpoly_cross_product(c(0, 0, 0), c(1, 0, 0), c(0, 1, 0)), sstr))
#> [1] "0" "0" "1"
intpoly_decompose(x^2 + x * y + x + y + x^3 * y^2 + y^5)
#> $`5`
#> x^3*y^2 + y^5
#> 
#> $`2`
#> x^2 + x*y
#> 
#> $`1`
#> x + y
intpoly_decompose(x^2 + x * y + x + y + x^3 * y^2 + y^5, TRUE)
#> [[1]]
#> x^3*y^2
#> 
#> [[2]]
#> x*y
#> 
#> [[3]]
#> x^2
#> 
#> [[4]]
#> y^5
#> 
#> [[5]]
#> x
#> 
#> [[6]]
#> y
```
