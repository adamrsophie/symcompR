# Differential geometry

``` r

library(symcompR)
t <- Symbol("t")
```

This vignette covers `sympy.diffgeom`: manifolds, patches, coordinate
systems, scalar and vector fields, differential forms, integral curves,
and the curvature of a metric.

## Fields and forms on \mathbb{R}^2

[`diffgeom_rn()`](https://adamrsophie.github.io/symcompR/reference/diffgeom.md)
returns the predefined manifold \mathbb{R}^2 with rectangular and polar
coordinate systems. Vector fields act on scalar fields as derivations,
and one-forms act on vector fields;
[`rcall()`](https://adamrsophie.github.io/symcompR/reference/diffgeom.md)
applies them.

``` r

rn <- diffgeom_rn()
fx <- base_scalar(rn$R2_r, 0); fy <- base_scalar(rn$R2_r, 1)
e_x <- base_vector(rn$R2_r, 0); e_y <- base_vector(rn$R2_r, 1)
dx <- base_oneform(rn$R2_r, 0); dy <- base_oneform(rn$R2_r, 1)
rcall(e_x, fx^2 * fy)
#> 2*x*y
rcall(dx, e_x)
#> 1
rcall(Differential(fx^2 * fy), e_x)
#> 2*x*y
rcall(Commutator(e_x, fx * e_y), fy)
#> 1
rcall(WedgeProduct(dx, dy), e_x, e_y)
#> 1
LieDerivative(e_x, fx * fy)
#> y
```

Coordinate changes and their Jacobians:

``` r

cs_transform(rn$R2_p, rn$R2_r)
#> Matrix(2, 1)
#> [rho*cos(theta)]
#> [rho*sin(theta)]
transformation(rn$R2_r, rn$R2_p)
#> Lambda((x, y), Matrix([
#> [sqrt(x^2 + y^2)],
#> [    atan2(y, x)]]))
jacobian_dg(rn$R2_p, rn$R2_r)
#> Matrix(2, 2)
#> [cos(theta), -rho*sin(theta)]
#> [sin(theta), rho*cos(theta) ]
simplify(jacobian_determinant(rn$R2_p, rn$R2_r))
#> rho
```

## Metrics and curvature

The Christoffel symbols of the Euclidean metric dr^2 + r^2 d\theta^2 in
polar coordinates:

``` r

dr <- base_oneform(rn$R2_p, 0); dth <- base_oneform(rn$R2_p, 1)
r <- base_scalar(rn$R2_p, 0)
Gam <- metric_to_Christoffel_2nd(TensorProduct(dr, dr) + r^2 * TensorProduct(dth, dth))
shape(Gam)
#> [1] 2 2 2
arr_get(Gam, c(1, 2, 2))
#> -rho
arr_get(Gam, c(2, 1, 2))
#> 1/rho
```

These are \Gamma^r\_{\theta\theta} = -r and \Gamma^\theta\_{r\theta} =
1/r. The only other nonzero symbol is \Gamma^\theta\_{\theta r} = 1/r,
by symmetry in the lower indices.

The first kind symbols, the Riemann tensor, and the Ricci tensor:

``` r

g_flat <- TensorProduct(dx, dx) + TensorProduct(dy, dy)
shape(metric_to_Christoffel_1st(g_flat))
#> [1] 2 2 2
Rie <- metric_to_Riemann_components(TensorProduct(dr, dr) + r^2 * TensorProduct(dth, dth))
shape(Rie)
#> [1] 2 2 2 2
unique(vapply(seq_len(16), function(k) sstr(arr_get(Rie, arrayInd(k, rep(2, 4))[1, ])), ""))
#> [1] "0"
Ric <- metric_to_Ricci_components(g_flat)
arr_get(Ric, c(1, 1))
#> 0
```

The polar metric is flat, so every component of its Riemann tensor is
zero.

The covariant derivative along a vector field uses the Christoffel
symbols:

``` r

ch <- metric_to_Christoffel_2nd(g_flat)
rcall(CovarDerivativeOp(fx * e_x, ch), fx)
#> x
rcall(BaseCovarDerivativeOp(rn$R2_r, 0, ch), fx)
#> 1
```

## Manifolds, patches, and coordinate systems

[`Manifold()`](https://adamrsophie.github.io/symcompR/reference/diffgeom.md),
[`Patch()`](https://adamrsophie.github.io/symcompR/reference/diffgeom.md),
and
[`CoordSystem()`](https://adamrsophie.github.io/symcompR/reference/diffgeom.md)
build new spaces. The relations argument lists the transformation
between two coordinate systems.

``` r

M <- Manifold("M", 2)
U <- Patch("U", M)
M
#> M
U
#> U
a <- Symbol("a"); b <- Symbol("b"); c <- Symbol("c"); d <- Symbol("d")
rel <- list(list(c("C1", "C2"), list(a, b), list(2 * a, 3 * b)))
C1 <- CoordSystem("C1", U, list(a, b), rel)
C2 <- CoordSystem("C2", U, list(c, d), rel)
coord_symbols(C1)
#> [[1]]
#> a
#> 
#> [[2]]
#> b
cs_transform(C1, C2)
#> Matrix(2, 1)
#> [2*a]
#> [3*b]
base_scalars(rn$R2_r)
#> [[1]]
#> x
#> 
#> [[2]]
#> y
base_vectors(rn$R2_r)
#> [[1]]
#> e_x
#> 
#> [[2]]
#> e_y
base_oneforms(rn$R2_p)
#> [[1]]
#> drho
#> 
#> [[2]]
#> dtheta
```

Points and coordinate changes:

``` r

x0 <- Symbol("x0"); y0 <- Symbol("y0")
pt <- point(rn$R2_r, list(x0, y0))
point_to_coords(rn$R2_p, pt)
#> Matrix(2, 1)
#> [sqrt(x0^2 + y0^2)]
#> [atan2(y0, x0)    ]
e_r <- base_vector(rn$R2_p, 0)
vectors_in_basis(e_r, rn$R2_r)
#> sin(theta)*e_y + cos(theta)*e_x
```

## Orders of tensor fields

``` r

contravariant_order(TensorProduct(e_x, e_y, e_x))
#> [1] 3
covariant_order(dx + dr)
#> [1] 1
covariant_order(WedgeProduct(dx, dy))
#> [1] 2
twoform_to_matrix(WedgeProduct(dx, dy))
#> Matrix(2, 2)
#> [0, -1]
#> [1, 0 ]
```

## Integral curves

[`intcurve_series()`](https://adamrsophie.github.io/symcompR/reference/diffgeom.md)
expands the integral curve of a vector field through a point as a power
series in the parameter, and
[`intcurve_diffequ()`](https://adamrsophie.github.io/symcompR/reference/diffgeom.md)
returns the differential equations that define it.

``` r

intcurve_series(e_x, t, point(rn$R2_r, list(1, 2)), n = 3)
#> Matrix(2, 1)
#> [t + 1]
#> [2    ]
intcurve_series(-fy * e_x + fx * e_y, t, point(rn$R2_r, list(1, 0)), n = 4)
#> Matrix(2, 1)
#> [1 - t^2/2 ]
#> [-t^3/6 + t]
intcurve_diffequ(-fy * e_x + fx * e_y, t, point(rn$R2_r, list(1, 0)))
#> [[1]]
#> [[1]][[1]]
#> f_1(t) + Derivative(f_0(t), t)
#> 
#> [[1]][[2]]
#> -f_0(t) + Derivative(f_1(t), t)
#> 
#> 
#> [[2]]
#> [[2]][[1]]
#> f_0(0) - 1
#> 
#> [[2]][[2]]
#> f_1(0)
```

The rotation field -y\\\partial_x + x\\\partial_y has the circle (\cos
t, \sin t) as its integral curve through (1, 0); the series above are
the Taylor polynomials of \cos t and \sin t.
