# Kinematic differential equations

The equations relating the time derivatives of orientation coordinates
to body fixed angular speeds, as SymPy's `kinematic_equations()`.

## Usage

``` r
kinematic_equations(speeds, coords, rot_type, rot_order = "")
```

## Arguments

- speeds:

  A list of three speeds.

- coords:

  A list of three (or four, for quaternions) coordinates.

- rot_type:

  `"body"`, `"space"` or `"quaternion"`.

- rot_order:

  Rotation order such as `"123"` or `"XYZ"`.

## Value

A list of expressions equal to zero.

## Examples

``` r
u <- dynamicsymbols("u1 u2 u3")
q <- dynamicsymbols("q1 q2 q3")
kinematic_equations(u, q, "body", "313")
#> [[1]]
#> -(u1(t)*sin(q3(t)) + u2(t)*cos(q3(t)))/sin(q2(t)) + Derivative(q1(t), t)
#> 
#> [[2]]
#> -u1(t)*cos(q3(t)) + u2(t)*sin(q3(t)) + Derivative(q2(t), t)
#> 
#> [[3]]
#> (u1(t)*sin(q3(t)) + u2(t)*cos(q3(t)))*cos(q2(t))/sin(q2(t)) - u3(t) + Derivative(q3(t), t)
#> 
```
