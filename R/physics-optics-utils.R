# Ray optics formulas, after sympy/physics/optics/utils.py.

.refractive_index_of_medium <- function(medium) {
  if (inherits(medium, "Medium")) medium$refractive_index else as_sym(medium)
}

# float(x) in Python: NULL when x is not a real number
.optics_float <- function(x) {
  if (is.numeric(x) && !is.object(x) && length(x) == 1L) return(as.numeric(x))
  if (is_sym(x) && !inherits(x, "pyobj") && x$type %in% c("Integer", "Rational", "Float")) return(as.numeric(N(x)))
  NULL
}

.optics_is_seq <- function(x) (is.list(x) && !is.object(x)) || inherits(x, "pytuple") ||
  (is.atomic(x) && !is.object(x) && length(x) > 1L)

.optics_col <- function(x) Matrix(lapply(as.list(unclass(x)), function(v) list(as_sym(v))))

.optics_dir <- function(r) .optics_col(lapply(.ldir(r), simplify))

# Point.__new__ turns Floats into Rationals
.optics_point <- function(coords) {
  do.call(Point, lapply(coords, function(c) {
    c <- as_sym(c)
    if (.has_type(c, "Float")) c <- simplify(nsimplify(c, rational = TRUE))
    c
  }))
}

.optics_ray <- function(p, drs) {
  d <- .optics_point(drs$data)
  Ray3D(p, .optics_point(Map(function(a, b) .add(list(a, b)), p$coords, d$coords)))
}

.optics_sumsq <- function(M) sqrt(Reduce(function(s, v) s + v^2, M$data, 0L))

#' Refraction and reflection at an interface
#'
#' `refraction_angle(incident, medium1, medium2, normal, plane)` is the
#' angle of refraction for an angle of incidence, or the direction (a matrix
#' or a `Ray3D`) of the refracted ray for an incident direction given as a
#' matrix, a sequence or a `Ray3D` with the normal or the plane of the
#' interface. `deviation()` is the angle of deviation of the refracted ray.
#' `fresnel_coefficients(angle, medium1, medium2)` is the list of reflection
#' and transmission coefficients `(R_p, R_s, T_p, T_s)`, or `(R_p, R_s)` past
#' the critical angle. `brewster_angle()` and `critical_angle()` are the
#' Brewster and critical angles. The media are [Medium()] objects or
#' refractive indices. After `sympy.physics.optics.utils`.
#'
#' @param incident The angle of incidence or the incident direction.
#' @param medium1,medium2 The media on the two sides.
#' @param normal The normal of the interface.
#' @param plane The `Plane` of the interface.
#' @param angle_of_incidence The angle of incidence.
#' @return A `sym` expression, a `symmatrix`, a `Ray3D` or a list.
#' @examples
#' refraction_angle(0.5, 1, 2)
#' refraction_angle(list(1, 1, 1), 1, 1, list(0, 0, 1))
#' critical_angle(1.33, 1)
#' @export
refraction_angle <- function(incident, medium1, medium2, normal = NULL, plane = NULL) {
  n1 <- .refractive_index_of_medium(medium1)
  n2 <- .refractive_index_of_medium(medium2)
  angle_of_incidence <- .optics_float(incident)
  critical_angle_ <- tryCatch(critical_angle(medium1, medium2), ValueError = function(e) NULL,
                              TypeError = function(e) NULL)
  if (!is.null(angle_of_incidence)) {
    if (!is.null(normal) || !is.null(plane)) .physics_error("ValueError", "Normal/plane not allowed if incident is an angle")
    if (!(0 <= angle_of_incidence && angle_of_incidence < pi * 0.5))
      .physics_error("ValueError", "Angle of incidence not in range [0:pi/2)")
    if (!is.null(critical_angle_) && !is_zero_num(critical_angle_) && isTRUE(angle_of_incidence > as.numeric(N(critical_angle_))))
      .physics_error("ValueError", "Ray undergoes total internal reflection")
    return(asin(n1 * sin(Float(angle_of_incidence)) / n2))
  }
  return_ray <- FALSE
  if (!is.null(plane) && !is.null(normal)) .physics_error("ValueError", "Either plane or normal is acceptable.")
  if (!inherits(incident, "symmatrix")) {
    if (.optics_is_seq(incident)) .incident <- .optics_col(incident)
    else if (inherits(incident, "GeoRay")) .incident <- .optics_dir(incident)
    else .physics_error("TypeError", "incident should be a Matrix, Ray3D, or sequence")
  } else {
    .incident <- incident
  }
  if (!is.null(plane)) {
    if (!inherits(plane, "GeoPlane")) .physics_error("TypeError", "plane should be an instance of geometry.plane.Plane")
    if (inherits(incident, "GeoRay")) {
      return_ray <- TRUE
      intersection_pt <- intersection(plane, incident)[[1L]]
    }
    .normal <- .optics_col(plane$normal)
  } else {
    if (!inherits(normal, "symmatrix")) {
      if (.optics_is_seq(normal)) {
        .normal <- .optics_col(normal)
      } else if (inherits(normal, "GeoRay")) {
        .normal <- .optics_dir(normal)
        if (inherits(incident, "GeoRay")) {
          intersection_pt <- intersection(incident, normal)
          if (length(intersection_pt) == 0L) .physics_error("ValueError", "Normal isn't concurrent with the incident ray.")
          return_ray <- TRUE
          intersection_pt <- intersection_pt[[1L]]
        }
      } else {
        .physics_error("TypeError", "Normal should be a Matrix, Ray3D, or sequence")
      }
    } else {
      .normal <- normal
    }
  }
  eta <- n1 / n2
  mag_incident <- .optics_sumsq(.incident)
  mag_normal <- .optics_sumsq(.normal)
  .incident <- .incident / mag_incident
  .normal <- .normal / mag_normal
  c1 <- -.add(Map(function(a, b) .mul(list(a, b)), .incident$data, .normal$data))
  cs2 <- 1 - eta^2 * (1 - c1^2)
  if (isTRUE(ask_fact(cs2, "negative"))) return(S_Zero)
  drs <- .incident * eta + .normal * (eta * c1 - sqrt(cs2))
  drs <- drs * mag_incident
  if (!return_ray) return(drs)
  .optics_ray(intersection_pt, drs)
}

#' @rdname refraction_angle
#' @export
fresnel_coefficients <- function(angle_of_incidence, medium1, medium2) {
  angle_of_incidence <- as_sym(angle_of_incidence)
  if (!(isTRUE(0 <= 2 * angle_of_incidence) && isTRUE(2 * angle_of_incidence < S_Pi)))
    .physics_error("ValueError", "Angle of incidence not in range [0:pi/2)")
  n1 <- .refractive_index_of_medium(medium1)
  n2 <- .refractive_index_of_medium(medium2)
  angle_of_refraction <- asin(n1 * sin(angle_of_incidence) / n2)
  onset <- tryCatch(critical_angle(n1, n2), ValueError = function(e) NULL)
  if (is.null(onset) || isTRUE(onset > angle_of_incidence)) {
    a <- angle_of_incidence; r <- angle_of_refraction
    R_s <- -sin(a - r) / sin(a + r)
    R_p <- tan(a - r) / tan(a + r)
    T_s <- 2 * sin(r) * cos(a) / sin(a + r)
    T_p <- 2 * sin(r) * cos(a) / (sin(a + r) * cos(a - r))
    return(list(R_p, R_s, T_p, T_s))
  }
  n <- n2 / n1
  a <- angle_of_incidence
  R_s <- cancel((cos(a) - S_I * sqrt(sin(a)^2 - n^2)) / (cos(a) + S_I * sqrt(sin(a)^2 - n^2)))
  R_p <- cancel((n^2 * cos(a) - S_I * sqrt(sin(a)^2 - n^2)) / (n^2 * cos(a) + S_I * sqrt(sin(a)^2 - n^2)))
  list(R_p, R_s)
}

#' @rdname refraction_angle
#' @export
deviation <- function(incident, medium1, medium2, normal = NULL, plane = NULL) {
  refracted <- refraction_angle(incident, medium1, medium2, normal = normal, plane = plane)
  angle_of_incidence <- .optics_float(incident)
  if (!is.null(angle_of_incidence)) return(Float(as.numeric(N(refracted))) - Float(angle_of_incidence))
  if (is_sym(refracted) && is_zero_num(refracted)) return(NULL)
  if (inherits(refracted, "GeoRay")) refracted <- .optics_dir(refracted)
  if (!inherits(incident, "symmatrix")) {
    if (.optics_is_seq(incident)) .incident <- .optics_col(incident)
    else if (inherits(incident, "GeoRay")) .incident <- .optics_dir(incident)
    else .physics_error("TypeError", "incident should be a Matrix, Ray3D, or sequence")
  } else {
    .incident <- incident
  }
  if (is.null(plane)) {
    if (!inherits(normal, "symmatrix")) {
      if (.optics_is_seq(normal)) .normal <- .optics_col(normal)
      else if (inherits(normal, "GeoRay")) .normal <- .optics_dir(normal)
      else .physics_error("TypeError", "normal should be a Matrix, Ray3D, or sequence")
    } else {
      .normal <- normal
    }
  } else {
    .normal <- .optics_col(plane$normal)
  }
  mag_incident <- .optics_sumsq(.incident)
  mag_normal <- .optics_sumsq(.normal)
  mag_refracted <- .optics_sumsq(refracted)
  .incident <- .incident / mag_incident
  .normal <- .normal / mag_normal
  refracted <- refracted / mag_refracted
  dot <- function(a, b) .add(Map(function(u, v) .mul(list(u, v)), a$data, b$data))
  i <- acos(dot(.incident, .normal))
  r <- acos(dot(refracted, .normal))
  i - r
}

#' @rdname refraction_angle
#' @export
brewster_angle <- function(medium1, medium2) {
  n1 <- .refractive_index_of_medium(medium1)
  n2 <- .refractive_index_of_medium(medium2)
  sym_atan2(n2, n1)
}

#' @rdname refraction_angle
#' @export
critical_angle <- function(medium1, medium2) {
  n1 <- .refractive_index_of_medium(medium1)
  n2 <- .refractive_index_of_medium(medium2)
  gt <- n2 > n1
  if (!is.logical(gt)) .physics_error("TypeError", "cannot determine truth value of Relational")
  if (isTRUE(gt)) .physics_error("ValueError", "Total internal reflection impossible for n1 < n2")
  asin(n2 / n1)
}

#' Lens and mirror formulas
#'
#' `lens_makers_formula(n_lens, n_surr, r1, r2, d)` is the focal length of a
#' lens with radii of curvature `r1` and `r2` and thickness `d`.
#' `mirror_formula()` and `lens_formula()` give the third of the focal length
#' `focal_length`, the object distance `u` and the image distance `v` from
#' the other two (any may be infinite). `hyperfocal_distance(f, N, c)` is the
#' hyperfocal distance and `transverse_magnification(si, so)` the transverse
#' magnification. After `sympy.physics.optics.utils`.
#'
#' @param n_lens,n_surr The refractive indices of the lens and its
#'   surroundings, or [Medium()] objects.
#' @param r1,r2 The radii of curvature of the two surfaces.
#' @param d The thickness of the lens.
#' @param focal_length,u,v The focal length and the object and image
#'   distances.
#' @param f The focal length.
#' @param N The f-number.
#' @param c The circle of confusion.
#' @param si,so The image and object distances.
#' @return A `sym` expression.
#' @examples
#' f <- Symbol("f"); u <- Symbol("u")
#' lens_formula(focal_length = f, u = u)
#' mirror_formula(u = oo, v = Symbol("v"))
#' transverse_magnification(30, 15)
#' @export
lens_makers_formula <- function(n_lens, n_surr, r1, r2, d = 0L) {
  n_lens <- .refractive_index_of_medium(n_lens)
  n_surr <- .refractive_index_of_medium(n_surr)
  d <- as_sym(d)
  # 1/r1 of a plain number is a Python float division
  inv <- function(r) {
    if (is_sym(r)) return(1 / r)
    if (r == 0) .physics_error("ZeroDivisionError", "division by zero")
    Float(1 / r)
  }
  focal_length <- 1 / ((n_lens - n_surr) / n_surr * (inv(r1) - inv(r2) + (((n_lens - n_surr) * d) / (n_lens * r1 * r2))))
  if (identical(as_sym(focal_length)$key, zoo$key)) return(oo)
  focal_length
}

.optics_is_oo <- function(x) !is.null(x) && identical(x$key, oo$key)

.optics_truthy <- function(x) !is.null(x) && !(is_sym(x) && is_zero_num(x)) && !(is.numeric(x) && x == 0)

.optics_conj <- function(focal_length, u, v, forms) {
  if (.optics_truthy(focal_length) && .optics_truthy(u) && .optics_truthy(v))
    .physics_error("ValueError", "Please provide only two parameters")
  if (!is.null(focal_length)) focal_length <- as_sym(focal_length)
  if (!is.null(u)) u <- as_sym(u)
  if (!is.null(v)) v <- as_sym(v)
  .u <- Symbol("u"); .v <- Symbol("v"); .f <- Symbol("f")
  lim <- function(e, x) limit(e, x, oo)
  if (is.null(focal_length)) {
    g <- forms$f
    if (.optics_is_oo(u) && .optics_is_oo(v)) return(lim(lim(g(.u, .v), .u), .v))
    if (.optics_is_oo(u)) return(lim(g(.u, v), .u))
    if (.optics_is_oo(v)) return(lim(g(u, .v), .v))
    return(g(u, v))
  }
  if (is.null(u)) {
    g <- forms$u
    if (.optics_is_oo(v) && .optics_is_oo(focal_length)) return(lim(lim(g(.v, .f), .v), .f))
    if (.optics_is_oo(v)) return(lim(g(.v, focal_length), .v))
    if (.optics_is_oo(focal_length)) return(lim(g(v, .f), .f))
    return(g(v, focal_length))
  }
  if (is.null(v)) {
    g <- forms$v
    if (.optics_is_oo(u) && .optics_is_oo(focal_length)) return(lim(lim(g(.u, .f), .u), .f))
    if (.optics_is_oo(u)) return(lim(g(.u, focal_length), .u))
    if (.optics_is_oo(focal_length)) return(lim(g(u, .f), .f))
    return(g(u, focal_length))
  }
  NULL
}

#' @rdname lens_makers_formula
#' @export
mirror_formula <- function(focal_length = NULL, u = NULL, v = NULL) {
  .optics_conj(focal_length, u, v, list(
    f = function(u, v) v * u / (v + u),
    u = function(v, f) v * f / (v - f),
    v = function(u, f) u * f / (u - f)))
}

#' @rdname lens_makers_formula
#' @export
lens_formula <- function(focal_length = NULL, u = NULL, v = NULL) {
  .optics_conj(focal_length, u, v, list(
    f = function(u, v) v * u / (u - v),
    u = function(v, f) v * f / (f - v),
    v = function(u, f) u * f / (u + f)))
}

#' @rdname lens_makers_formula
#' @export
hyperfocal_distance <- function(f, N, c) {
  f <- as_sym(f); N <- as_sym(N); c <- as_sym(c)
  (1 / (N * c)) * (f^2)
}

#' @rdname lens_makers_formula
#' @export
transverse_magnification <- function(si, so) {
  si <- as_sym(si); so <- as_sym(so)
  -(si / so)
}
