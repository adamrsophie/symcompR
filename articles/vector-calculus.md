# Vector calculus

``` r

library(symcompR)
x <- Symbol("x"); y <- Symbol("y"); t <- Symbol("t"); theta <- Symbol("theta")
a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c")
```

This vignette covers `sympy.vector`: coordinate systems, vectors and
dyadics, the differential operators, orientation of coordinate systems,
and integration over parametric and implicit regions.

## Coordinate systems and vectors

[`CoordSys3D()`](https://adamrsophie.github.io/symcompR/reference/CoordSys3D.md)
creates a Cartesian coordinate system with base scalars `N$x`, `N$y`,
`N$z` and base vectors `N$i`, `N$j`, `N$k`. Vectors are built from them
with ordinary arithmetic.

``` r

N <- CoordSys3D("N")
v <- 3 * N$i + 4 * N$j
v
#> 3*N.i + 4*N.j
magnitude(v)
#> 5
normalize(v)
#> 3/5*N.i + 4/5*N.j
dot(N$i + 2 * N$j, 3 * N$i - N$k)
#> 3
cross(N$i, N$j)
#> N.k
cross(N$i + N$j, N$j + N$k)
#> N.i + (-1)*N.j + N.k
to_matrix(N$x * N$i + 2 * N$j, N)
#> Matrix(3, 1)
#> [N.x]
#> [2  ]
#> [0  ]
matrix_to_vector(Matrix(list(1, 2, 3)), N)
#> N.i + 2*N.j + 3*N.k
```

The classes behind vectors can be used directly.
[`BaseScalar()`](https://adamrsophie.github.io/symcompR/reference/CoordSys3D.md)
and
[`BaseVector()`](https://adamrsophie.github.io/symcompR/reference/VectorAdd.md)
are the coordinate variables and unit vectors,
[`VectorAdd()`](https://adamrsophie.github.io/symcompR/reference/VectorAdd.md),
[`VectorMul()`](https://adamrsophie.github.io/symcompR/reference/VectorAdd.md),
and
[`VectorZero()`](https://adamrsophie.github.io/symcompR/reference/VectorAdd.md)
build vectors,
[`measure_number()`](https://adamrsophie.github.io/symcompR/reference/VectorAdd.md)
is the coefficient of a base vector, and
[`vector_equals()`](https://adamrsophie.github.io/symcompR/reference/magnitude.md)
tests equality.
[`Dot()`](https://adamrsophie.github.io/symcompR/reference/Cross.md) and
[`Cross()`](https://adamrsophie.github.io/symcompR/reference/Cross.md)
are the unevaluated products.

``` r

BaseScalar(0, N)
#> N.x
BaseVector(1, N)
#> N.j
VectorAdd(a * N$i, b * N$j)
#> a*N.i + b*N.j
VectorMul(3, N$k)
#> 3*N.k
VectorZero()
#> 0
measure_number(a * N$i)
#> a
vector_equals(N$i + N$j, N$j + N$i)
#> [1] TRUE
Dot(N$i, N$j)
#> Dot(N.i, N.j)
doit(Dot(N$i + N$j, N$j))
#> 1
Cross(N$i, N$j)
#> Cross(N.i, N.j)
doit(Cross(N$i, N$j))
#> N.k
components(a * N$i + b * N$k)
#> {N.i: a, N.k: b}
base_scalars(N)
#> (N.x, N.y, N.z)
base_vectors(N)
#> (N.i, N.j, N.k)
orthogonalize(N$i + N$j, N$i, orthonormal = TRUE)
#> [[1]]
#> (sqrt(2)/2)*N.i + (sqrt(2)/2)*N.j
#> 
#> [[2]]
#> (sqrt(2)/2)*N.i + (-sqrt(2)/2)*N.j
```

## Differential operators

``` r

grad(N$x^2 * N$y + N$z)
#> 2*N.x*N.y*N.i + N.x^2*N.j + N.k
gradient(N$x^2 * N$y)
#> 2*N.x*N.y*N.i + N.x^2*N.j
divergence(N$x * N$i + N$y * N$j + N$z * N$k)
#> 3
curl(-N$y * N$i + N$x * N$j)
#> 2*N.k
laplacian(N$x^2 + N$y^2 + N$z^2)
#> 6
directional_derivative(N$x * N$y * N$z, N$i + N$j)
#> N.x*N.z + N.y*N.z
```

A field is conservative when it is the gradient of a scalar potential,
and solenoidal when its divergence vanishes.

``` r

F <- N$y * N$z * N$i + N$x * N$z * N$j + N$x * N$y * N$k
is_conservative(F)
#> [1] TRUE
scalar_potential(F, N)
#> N.x*N.y*N.z
is_solenoidal(-N$y * N$i + N$x * N$j)
#> [1] TRUE
```

[`separate()`](https://adamrsophie.github.io/symcompR/reference/magnitude.md)
splits a vector into its components in each coordinate system it
involves.

``` r

separate(N$x * N$i + N$y * N$j)
#> {N: N.x*N.i + N.y*N.j}
```

The operators also exist as unevaluated objects,
[`Gradient()`](https://adamrsophie.github.io/symcompR/reference/Gradient.md),
[`Divergence()`](https://adamrsophie.github.io/symcompR/reference/Gradient.md),
[`Curl()`](https://adamrsophie.github.io/symcompR/reference/Gradient.md),
and
[`Laplacian()`](https://adamrsophie.github.io/symcompR/reference/Gradient.md),
which
[`doit()`](https://adamrsophie.github.io/symcompR/reference/doit.md)
evaluates, and as methods of the del operator
[`Del()`](https://adamrsophie.github.io/symcompR/reference/Del.md).

``` r

Gradient(N$x * N$y * N$z)
#> Gradient(N.x*N.y*N.z)
doit(Divergence(N$x * N$i + N$z^2 * N$j))
#> 1
doit(Curl(N$x * N$i + N$y * N$j + N$z * N$k))
#> 0
doit(Laplacian(N$x^2 + N$y^2 + N$z^2))
#> 6
delop <- Del()
del_gradient(delop, N$x^2 * N$y, doit = TRUE)
#> 2*N.x*N.y*N.i + N.x^2*N.j
del_dot(delop, N$x * N$i + N$y * N$j, doit = TRUE)
#> 2
del_cross(delop, -N$y * N$i + N$x * N$j, doit = TRUE)
#> 2*N.k
scalar_potential_difference(F, N, N$origin, locate_new(N, "P1", N$i + N$j + N$k)$origin)
#> 1
```

## Rotated, translated, and curvilinear systems

[`orient_new_axis()`](https://adamrsophie.github.io/symcompR/reference/CoordSys3D.md)
creates a coordinate system rotated about an axis, and
[`express()`](https://adamrsophie.github.io/symcompR/reference/express.md)
rewrites vectors in another system.

``` r

A <- orient_new_axis(N, "A", theta, N$k)
express(A$i, N)
#> (cos(theta))*N.i + (sin(theta))*N.j
rotation_matrix(A, N)
#> Matrix(3, 3)
#> [cos(theta) , sin(theta), 0]
#> [-sin(theta), cos(theta), 0]
#> [0          , 0         , 1]
```

[`create_new()`](https://adamrsophie.github.io/symcompR/reference/CoordSys3D.md)
with a `transformation` defines curvilinear coordinates; the Lamé
coefficients are the scale factors of the coordinate lines.

``` r

Sph <- create_new(N, "S", transformation = "spherical")
lame_coefficients(Sph)
#> (1, S.r, S.r*sin(S.theta))
```

[`orient_new()`](https://adamrsophie.github.io/symcompR/reference/CoordSys3D.md)
takes a list of orienters, which describe the rotation:
[`AxisOrienter()`](https://adamrsophie.github.io/symcompR/reference/AxisOrienter.md)
(about an axis),
[`BodyOrienter()`](https://adamrsophie.github.io/symcompR/reference/AxisOrienter.md)
(three successive rotations about the body axes),
[`SpaceOrienter()`](https://adamrsophie.github.io/symcompR/reference/AxisOrienter.md)
(about the space axes), and
[`QuaternionOrienter()`](https://adamrsophie.github.io/symcompR/reference/AxisOrienter.md).
The `orient_new_*()` functions are shortcuts.

``` r

q1 <- Symbol("q1"); q2 <- Symbol("q2"); q3 <- Symbol("q3"); q0 <- Symbol("q0")
orienter_rotation_matrix(AxisOrienter(a, N$k), N)
#> Matrix(3, 3)
#> [cos(a) , sin(a), 0]
#> [-sin(a), cos(a), 0]
#> [0      , 0     , 1]
B <- orient_new_body(N, "B", q1, 0, 0, "123")
rotation_matrix(B, N)
#> Matrix(3, 3)
#> [1, 0       , 0      ]
#> [0, cos(q1) , sin(q1)]
#> [0, -sin(q1), cos(q1)]
Sp <- orient_new_space(N, "Sp", q1, 0, 0, "123")
rotation_matrix(Sp, N)
#> Matrix(3, 3)
#> [1, 0       , 0      ]
#> [0, cos(q1) , sin(q1)]
#> [0, -sin(q1), cos(q1)]
D <- orient_new(N, "D", list(AxisOrienter(q1, N$j)))
rotation_matrix(D, N)
#> Matrix(3, 3)
#> [cos(q1), 0, -sin(q1)]
#> [0      , 1, 0       ]
#> [sin(q1), 0, cos(q1) ]
vector_equals(express(D$k, N), express(orient_new_axis(N, "D2", q1, N$j)$k, N))
#> [1] TRUE
Dq <- orient_new_quaternion(N, "Q", 1, 0, 0, 0)
rotation_matrix(Dq, N)
#> Matrix(3, 3)
#> [1, 0, 0]
#> [0, 1, 0]
#> [0, 0, 1]
O1 <- orient_new(N, "O1", list(BodyOrienter(q1, q2, q3, "123")))
O2 <- orient_new_body(N, "O2", q1, q2, q3, "123")
simplify(rotation_matrix(O1, N) - rotation_matrix(O2, N))
#> Matrix(3, 3)
#> [0, 0, 0]
#> [0, 0, 0]
#> [0, 0, 0]
O3 <- orient_new(N, "O3", list(SpaceOrienter(q1, 0, 0, "123")))
rotation_matrix(O3, N)
#> Matrix(3, 3)
#> [1, 0       , 0      ]
#> [0, cos(q1) , sin(q1)]
#> [0, -sin(q1), cos(q1)]
O4 <- orient_new(N, "O4", list(QuaternionOrienter(1, 0, 0, 0)))
rotation_matrix(O4, N)
#> Matrix(3, 3)
#> [1, 0, 0]
#> [0, 1, 0]
#> [0, 0, 1]
```

[`locate_new()`](https://adamrsophie.github.io/symcompR/reference/CoordSys3D.md)
translates the origin.
[`position_wrt()`](https://adamrsophie.github.io/symcompR/reference/lame_coefficients.md)
gives the position of one origin relative to another,
[`express_coordinates()`](https://adamrsophie.github.io/symcompR/reference/lame_coefficients.md)
the coordinates of a point, and
[`scalar_map()`](https://adamrsophie.github.io/symcompR/reference/lame_coefficients.md)
how the coordinate variables of two systems are related.

``` r

C <- locate_new(N, "C", a * N$i + b * N$j + c * N$k)
position_wrt(C, N)
#> a*N.i + b*N.j + c*N.k
express_coordinates(N$origin, C)
#> (-a, -b, -c)
scalar_map(N, C)
#> {N.x: C.x + a, N.y: C.y + b, N.z: C.z + c}
```

A coordinate system defined by a transformation keeps the map from its
parent, and
[`change_of_basis_matrix_from()`](https://adamrsophie.github.io/symcompR/reference/lame_coefficients.md)
relates the base vectors.

``` r

Cyl <- create_new(N, "Cyl", transformation = "cylindrical")
transformation_to_parent(Cyl)
#> (Cyl.r*cos(Cyl.theta), Cyl.r*sin(Cyl.theta), Cyl.z)
vapply(transformation_from_parent(create_new(N, "Sp2", transformation = "spherical")), sstr, "")
#> [1] "sqrt(N.x^2 + N.y^2 + N.z^2)"          
#> [2] "acos(N.z/sqrt(N.x^2 + N.y^2 + N.z^2))"
#> [3] "atan2(N.y, N.x)"
transformation_from_parent_function(CoordSys3D("Ca", "cartesian"))(x, y, t)
#> (x, y, t)
Sw <- create_new(N, "Sw", transformation = function(x, y, z) list(y, x, z))
change_of_basis_matrix_from(Sw, N)
#> Matrix(3, 3)
#> [0, 1, 0]
#> [1, 0, 0]
#> [0, 0, 1]
lame_coefficients(Cyl)
#> (1, Cyl.r, 1)
```

## Dyadics

The outer product of two vectors is a dyadic; `&` is the dot product
with a vector.

``` r

outer_sym(N$i, N$j)
#> (N.i|N.j)
(N$i | N$j) & N$j
#> N.i
```

[`BaseDyadic()`](https://adamrsophie.github.io/symcompR/reference/VectorAdd.md),
[`DyadicAdd()`](https://adamrsophie.github.io/symcompR/reference/VectorAdd.md),
[`DyadicMul()`](https://adamrsophie.github.io/symcompR/reference/VectorAdd.md),
and
[`DyadicZero()`](https://adamrsophie.github.io/symcompR/reference/VectorAdd.md)
are the dyadic counterparts of the vector classes, and
[`matrix_to_dyadic()`](https://adamrsophie.github.io/symcompR/reference/express.md)
builds a dyadic from its matrix of components.

``` r

BaseDyadic(N$i, N$j)
#> (N.i|N.j)
DyadicAdd(BaseDyadic(N$i, N$j), BaseDyadic(N$j, N$i))
#> (N.i|N.j) + (N.j|N.i)
DyadicMul(2, BaseDyadic(N$k, N$k))
#> 2*(N.k|N.k)
DyadicZero()
#> 0
to_matrix(matrix_to_dyadic(Matrix(list(list(1, 2, 0), list(0, 1, 0), list(0, 0, 3))), N), N)
#> Matrix(3, 3)
#> [1, 2, 0]
#> [0, 1, 0]
#> [0, 0, 3]
```

## Integrals over regions

[`ParametricRegion()`](https://adamrsophie.github.io/symcompR/reference/ParametricRegion.md)
describes a curve, surface, or volume by a parametrization, and
[`vector_integrate()`](https://adamrsophie.github.io/symcompR/reference/ParametricIntegral.md)
integrates scalar or vector fields over it. Geometry objects such as
segments are accepted as well.

``` r

quarter <- ParametricRegion(list(cos(t), sin(t)), list(t, 0, sym_pi / 2))
vector_integrate(N$x * N$y, quarter)
#> 1/2
disk <- ParametricRegion(list(t * cos(theta), t * sin(theta)),
                         list(t, 0, 2), list(theta, 0, 2 * sym_pi))
vector_integrate(1, disk)
#> 4*pi
vector_integrate(N$x^2, Segment(Point(0, 0), Point(1, 1)))
#> sqrt(2)/3
ParametricIntegral(N$x * N$i + N$y * N$j, ParametricRegion(list(t, t), list(t, 0, 1)))
#> 1
```

[`ImplicitRegion()`](https://adamrsophie.github.io/symcompR/reference/ImplicitRegion.md)
describes a region by an equation.
[`regular_point()`](https://adamrsophie.github.io/symcompR/reference/ImplicitRegion.md)
finds a point on it,
[`singular_points()`](https://adamrsophie.github.io/symcompR/reference/ImplicitRegion.md)
its singular points,
[`implicit_degree()`](https://adamrsophie.github.io/symcompR/reference/ImplicitRegion.md)
the degree of the equation, and
[`rational_parametrization()`](https://adamrsophie.github.io/symcompR/reference/ImplicitRegion.md)
a rational parametrization of a conic.
[`conic_coeff()`](https://adamrsophie.github.io/symcompR/reference/ImplicitRegion.md)
returns the coefficients of a conic equation and
[`parametric_region_list()`](https://adamrsophie.github.io/symcompR/reference/ParametricRegion.md)
converts geometry objects to parametric regions.

``` r

circ <- ImplicitRegion(list(x, y), x^2 + y^2 - 4)
regular_point(circ)
#> (0, -2)
singular_points(ImplicitRegion(list(x, y), y^2 - x^3))
#> {(0, 0)}
implicit_degree(circ)
#> 2
rational_parametrization(circ)
#> (4*t/(t^2 + 1), 4*t^2/(t^2 + 1) - 2)
conic_coeff(list(x, y), x^2 - 4 * x * y - 3 * y^2 + 4 * x + 8 * y - 5)
#> (1, -4, -3, 4, 8, -5)
parametric_region_list(Point(-5, 12))
#> [[1]]
#> ParametricRegion((-5, 12))
parametric_region_list(Segment(Point(0, 0), Point(1, 3)))
#> [[1]]
#> ParametricRegion((t, 3*t), (t, 0, 1))
```
