# Base classes of sympy.physics: the mechanics bases (but LoadBase, a
# namedtuple) are Python ABCs, so constructing one raises TypeError as in SymPy; their names are the
# S3 classes the concrete objects inherit. BosonicBasis is the plain base
# class of the second quantization bases.

.abstract_base <- function(name, methods)
  function(...) stop(.sympy_condition("TypeError", sprintf(
    "Can't instantiate abstract class %s without an implementation for abstract method%s %s", name,
    if (length(methods) > 1L) "s" else "", paste0("'", methods, "'", collapse = ", "))))

#' Abstract base classes of the physics modules
#'
#' `ActuatorBase()`, `PathwayBase()`, `WrappingGeometryBase()`,
#' `BodyBase()` and `MethodBase()` are SymPy's abstract base classes of
#' actuators, pathways, wrapping geometries, bodies and equation of motion
#' methods: calling one raises a `TypeError`, and their
#' names are the classes that the concrete objects inherit, so
#' `inherits(x, "PathwayBase")` tests for a pathway. `LoadBase(location,
#' vector)` is the named pair that [Force()] and [Torque()] extend. `BosonicBasis()` is the
#' base class of [VarBosonicBasis()] and [FixedBosonicBasis()]. `Vector()`
#' and `Dyadic()` are the vectors and dyadics of `sympy.physics.vector`,
#' [PhysVector()] and [PhysDyadic()].
#'
#' @param ... Ignored; abstract classes cannot be built.
#' @param location,vector The point or frame and the vector of a load.
#' @param inlist See [PhysVector()] and [PhysDyadic()].
#' @return `BosonicBasis()` an empty basis; `LoadBase()` a load; `Vector()` and `Dyadic()` a
#'   vector and a dyadic; the others signal an error.
#' @examples
#' tryCatch(PathwayBase(), error = conditionMessage)
#' inherits(LinearPathway(PhysPoint("pA"), PhysPoint("pB")), "PathwayBase")
#' inherits(VarBosonicBasis(3), "BosonicBasis")
#' N <- ReferenceFrame("N")
#' Vector(list(list(Matrix(list(1, 2, 0)), N)))
#' @name physics-abstract-bases
NULL

#' @rdname physics-abstract-bases
#' @export
ActuatorBase <- .abstract_base("ActuatorBase", "to_loads")
#' @rdname physics-abstract-bases
#' @export
PathwayBase <- .abstract_base("PathwayBase", c("extension_velocity", "length", "to_loads"))
#' @rdname physics-abstract-bases
#' @export
WrappingGeometryBase <- .abstract_base("WrappingGeometryBase",
  c("geodesic_end_vectors", "geodesic_length", "point"))
#' @rdname physics-abstract-bases
#' @export
BodyBase <- .abstract_base("BodyBase", c("angular_momentum", "kinetic_energy", "linear_momentum", "parallel_axis"))
#' @rdname physics-abstract-bases
#' @export
LoadBase <- function(location, vector) .pm_namedtuple("LoadBase", list(location = location, vector = vector))
#' @rdname physics-abstract-bases
#' @export
MethodBase <- .abstract_base("MethodBase", c("_form_eoms", "acceleration_constraints", "bodies", "forcing", "forcing_full", "frame",
  "holonomic_constraints", "loads", "mass_matrix", "mass_matrix_full", "nonholonomic_constraints", "q", "u",
  "velocity_constraints"))
#' @rdname physics-abstract-bases
#' @export
BosonicBasis <- function() structure(list(), class = "BosonicBasis")
#' @rdname physics-abstract-bases
#' @export
Vector <- function(inlist) PhysVector(inlist)
#' @rdname physics-abstract-bases
#' @export
Dyadic <- function(inlist) PhysDyadic(inlist)
.pvR_LoadBase <- function(x) sprintf("LoadBase(location=%s, vector=%s)", .pv_str(x$location), .pv_str(x$vector))
