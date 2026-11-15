# Abstract base classes of the physics modules

`ActuatorBase()`, `PathwayBase()`, `WrappingGeometryBase()`,
`BodyBase()` and `MethodBase()` are SymPy's abstract base classes of
actuators, pathways, wrapping geometries, bodies and equation of motion
methods: calling one raises a `TypeError`, and their names are the
classes that the concrete objects inherit, so
`inherits(x, "PathwayBase")` tests for a pathway.
`LoadBase(location, vector)` is the named pair that
[`Force()`](https://adamrsophie.github.io/symcompR/reference/Force.md)
and
[`Torque()`](https://adamrsophie.github.io/symcompR/reference/Force.md)
extend. `BosonicBasis()` is the base class of
[`VarBosonicBasis()`](https://adamrsophie.github.io/symcompR/reference/secondquant.md)
and
[`FixedBosonicBasis()`](https://adamrsophie.github.io/symcompR/reference/secondquant.md).
`Vector()` and `Dyadic()` are the vectors and dyadics of
`sympy.physics.vector`,
[`PhysVector()`](https://adamrsophie.github.io/symcompR/reference/PhysVector.md)
and
[`PhysDyadic()`](https://adamrsophie.github.io/symcompR/reference/PhysVector.md).

## Usage

``` r
ActuatorBase(...)

PathwayBase(...)

WrappingGeometryBase(...)

BodyBase(...)

LoadBase(location, vector)

MethodBase(...)

BosonicBasis()

Vector(inlist)

Dyadic(inlist)
```

## Arguments

- ...:

  Ignored; abstract classes cannot be built.

- location, vector:

  The point or frame and the vector of a load.

- inlist:

  See
  [`PhysVector()`](https://adamrsophie.github.io/symcompR/reference/PhysVector.md)
  and
  [`PhysDyadic()`](https://adamrsophie.github.io/symcompR/reference/PhysVector.md).

## Value

`BosonicBasis()` an empty basis; `LoadBase()` a load; `Vector()` and
`Dyadic()` a vector and a dyadic; the others signal an error.

## Examples

``` r
tryCatch(PathwayBase(), error = conditionMessage)
#> [1] "Can't instantiate abstract class PathwayBase without an implementation for abstract methods 'extension_velocity', 'length', 'to_loads'"
inherits(LinearPathway(PhysPoint("pA"), PhysPoint("pB")), "PathwayBase")
#> [1] TRUE
inherits(VarBosonicBasis(3), "BosonicBasis")
#> [1] TRUE
N <- ReferenceFrame("N")
Vector(list(list(Matrix(list(1, 2, 0)), N)))
#> N.x + 2*N.y
```
