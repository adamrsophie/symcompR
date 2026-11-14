# Vector calculus in Cartesian coordinates, after sympy/vector. The geometry
# entities are in matgeo-geometry-*.R.

# ---------------------------------------------------------------------------
# Vector calculus

#' Vector calculus in Cartesian coordinates
#'
#' Vectors are column matrices of three components. `grad()`, `divergence()`,
#' `curl()` and `laplacian()` differentiate with respect to `vars` (default
#' `x`, `y`, `z`). `dot()` and `cross()` come from the matrix module.
#' `is_conservative()` and `is_solenoidal()` test fields, and
#' `scalar_potential()` integrates a conservative field.
#'
#' @param f A scalar expression.
#' @param v A vector field: a list of three expressions or a column matrix.
#' @param vars A list of three symbols. When `vars` is not given and the
#'   field is a [CoordSys3D()] vector or refers to base scalars, the
#'   functions of `sympy.vector` are used instead; `scalar_potential(v,
#'   coord_sys)` then takes the coordinate system.
#' @param doit Evaluate the derivatives (coordinate system fields only).
#' @return An expression or matrix.
#' @examples
#' x <- Symbol("x"); y <- Symbol("y"); z <- Symbol("z")
#' curl(list(-y, x, 0))
#' divergence(list(x, y, z))
#' @export
grad <- function(f, vars = .xyz(), doit = TRUE) {
  if (missing(vars) && .is_vector_field(f)) return(.vec_gradient(f, doit = doit))
  gradient(f, vars)
}

.xyz <- function() list(Symbol("x"), Symbol("y"), Symbol("z"))

.vcomp <- function(v) if (inherits(v, "symmatrix")) v$data else lapply(v, as_sym)

#' @rdname grad
#' @export
divergence <- function(v, vars = .xyz(), doit = TRUE) {
  if (.pv_dispatch(v, vars)) return(.pv_divergence(v, vars))
  if (missing(vars) && .is_vector_field(v)) return(.vec_divergence(v, doit = doit))
  v <- .vcomp(v)
  .add(Map(function(c, s) .diff1(c, as_sym(s)), v, vars))
}

#' @rdname grad
#' @export
curl <- function(v, vars = .xyz(), doit = TRUE) {
  if (.pv_dispatch(v, vars)) return(.pv_curl(v, vars))
  if (missing(vars) && .is_vector_field(v)) return(.vec_curl(v, doit = doit))
  v <- .vcomp(v)
  x <- as_sym(vars[[1L]]); y <- as_sym(vars[[2L]]); z <- as_sym(vars[[3L]])
  Matrix(list(.sub(.diff1(v[[3L]], y), .diff1(v[[2L]], z)),
              .sub(.diff1(v[[1L]], z), .diff1(v[[3L]], x)),
              .sub(.diff1(v[[2L]], x), .diff1(v[[1L]], y))))
}

#' @rdname grad
#' @export
laplacian <- function(f, vars = .xyz()) {
  if (missing(vars) && .is_vector_field(f)) return(.vec_laplacian(f))
  f <- as_sym(f)
  .add(lapply(vars, function(s) .diff1(.diff1(f, as_sym(s)), as_sym(s))))
}

#' @rdname grad
#' @export
is_conservative <- function(v, vars = .xyz()) {
  if (.pv_dispatch(v) || .pv_is0(v)) return(.pv_is_conservative(v))
  if (is_sym(v)) return(.vec_is_conservative(v))
  all(vapply(Matrix(curl(v, vars))$data, function(z) is_zero_num(simplify(z)), TRUE))
}

#' @rdname grad
#' @export
is_solenoidal <- function(v, vars = .xyz()) {
  if (.pv_dispatch(v) || .pv_is0(v)) return(.pv_is_solenoidal(v))
  if (is_sym(v)) return(.vec_is_solenoidal(v))
  is_zero_num(simplify(divergence(v, vars)))
}

#' @rdname grad
#' @export
scalar_potential <- function(v, vars = .xyz()) {
  if (.pv_dispatch(v, vars)) return(.pv_scalar_potential(v, vars))
  if (is_sym(v)) return(.vec_scalar_potential(v, vars))
  if (!is_conservative(v, vars)) stop("field is not conservative")
  v <- .vcomp(v)
  out <- S_Zero
  for (i in seq_along(vars)) {
    s <- as_sym(vars[[i]])
    rem <- simplify(.sub(v[[i]], .diff1(out, s)))
    out <- .add(list(out, integrate(rem, s)))
  }
  out
}
