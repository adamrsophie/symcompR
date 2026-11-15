# Planes, curves, parabolas and the geometry utilities, after
# sympy/geometry/plane.py, curve.py, parabola.py and util.py, with the
# functions that dispatch SymPy's entity methods.

.mg_plane_raw <- function(p1, n) {
  .mg_gnew("Plane", "GeoPlane", list(p1 = p1, normal_vector = n, p = p1, normal = n))
}

#' Planes
#'
#' `Plane()` builds a plane through three points or through a point with a
#' `normal_vector`, as SymPy does. The point and the normal vector are in
#' `$p1` and `$normal_vector`. The functions are SymPy's methods:
#' `is_coplanar()`, `parallel_plane()`, `perpendicular_plane()`,
#' `projection_line()` and the shared [distance()], [intersection()],
#' [projection()], [angle_between()], [is_parallel()],
#' [is_perpendicular()], [perpendicular_line()] and [arbitrary_point()].
#'
#' @param p1 A point of the plane.
#' @param a,b Two more points, when `normal_vector` is not given.
#' @param normal_vector The normal vector, a list or vector of length 3.
#' @param evaluate Make a point of the normal vector (rationalizing floats).
#' @param e,pl A plane.
#' @param o,pt,line An entity, a point or a line.
#' @param ... Points (for `perpendicular_plane()`) or planes.
#' @return An entity, an expression, a logical value or a list.
#' @examples
#' p <- Plane(Point3D(1, 1, 1), normal_vector = c(1, 4, 7))
#' equation(p)
#' distance(p, Point3D(1, 2, 3))
#' intersection(Plane(Point3D(0, 0, 0), normal_vector = c(1, 0, 0)),
#'              Plane(Point3D(0, 0, 0), normal_vector = c(0, 1, 0)))
#' @export
Plane <- function(p1, a = NULL, b = NULL, normal_vector = NULL, evaluate = TRUE) {
  p1 <- .mg_as_pt(p1, 3L)
  if (!is.null(a) && !is.null(b)) {
    p2 <- .mg_as_pt(a, 3L); p3 <- .mg_as_pt(b, 3L)
    if (is_collinear(p1, p2, p3)) .mg_err("ValueError", "Enter three non-collinear points")
    n <- cross(Matrix(direction_ratio(p1, p2)), Matrix(direction_ratio(p1, p3)))$data
  } else {
    if (is.null(normal_vector)) normal_vector <- a
    sq <- .mg_seq(normal_vector)
    if (is.null(sq) || length(sq) != 3L)
      .mg_errf("ValueError", "Either provide 3 3D points or a point with a normal vector expressed as a sequence of length 3")
    n <- if (isTRUE(evaluate)) Point3D(sq)$coords else lapply(sq, as_sym)
    if (all(vapply(n, function(z) isTRUE(ask_fact(z, "zero")), TRUE))) .mg_err("ValueError", "Normal vector cannot be zero vector")
  }
  .mg_plane_raw(p1, n)
}

.mg_plane_xyz <- function() lapply(1:3, function(i) Dummy("plane_dummy"))

.mg_plane_equation <- function(e, x = NULL, y = NULL, z = NULL) {
  given <- Filter(Negate(is.null), list(x, y, z))
  fk <- vapply(.mg_geo_free(e), .k, "")
  if (any(vapply(given, function(s) is_sym(s) && s$key %in% fk, TRUE)))
    .mg_err("ValueError", "Provided symbols clash with symbols used to define the plane")
  v <- Map(function(s, nm) if (is.null(s)) Symbol(nm, real = TRUE) else .mg_param(s), list(x, y, z), c("x", "y", "z"))
  d <- direction_ratio(e$p1, Point3D(v))
  .add(Map(function(i, j) .mul(list(i, j)), d, e$normal_vector))
}

.mg_plane_has <- function(e, o) {
  xyz <- .mg_plane_xyz()
  k <- .mg_plane_equation(e, xyz[[1L]], xyz[[2L]], xyz[[3L]])
  if (inherits(o, "GeoLine")) {
    d <- Point3D(.mg_line_arbitrary(o, Dummy("plane_dummy"))$coords)
    return(isTRUE(.mg_eq0(subs(k, Map(list, xyz, d$coords)))))
  }
  o <- tryCatch(Point(o, dim = 3L), error = function(err) NULL)
  if (is.null(o)) return(FALSE)
  isTRUE(.mg_eq0(xreplace(k, Map(list, xyz, o$coords))))
}

# Expr.equals(0): exact when the difference simplifies to a number
.mg_eq0 <- function(e) {
  d <- simplify(as_sym(e))
  if (is_exact_number(d)) return(is_zero_num(d))
  .expr_equals(d, S_Zero)
}

.mg_plane_angle <- function(e, o) {
  n <- e$normal_vector
  nrm <- function(v) .sqrt(.add(lapply(v, function(i) .pow(i, S_Two))))
  dt <- function(a, b) .add(Map(function(u, w) .mul(list(u, w)), a, b))
  if (inherits(o, "GeoLinearEntity3D")) {
    dr <- direction_ratio(o)
    return(asin(.div(dt(n, dr), .mul(list(nrm(n), nrm(dr))))))
  }
  if (inherits(o, "GeoPlane")) {
    m <- o$normal_vector
    return(acos(.div(dt(n, m), .mul(list(nrm(n), nrm(m))))))
  }
  NULL
}

.mg_plane_arbitrary <- function(e, u = NULL, v = NULL) {
  circle <- is.null(v)
  if (circle) u <- .mg_param(u, "t")
  else { u <- .mg_param(u, "u"); v <- .mg_param(v, "v") }
  n <- e$normal_vector
  x <- n[[1L]]; y <- n[[2L]]; z <- n[[3L]]
  pc <- e$p1$coords
  w1 <- if (isTRUE(ask_fact(x, "zero")) && isTRUE(ask_fact(y, "zero"))) list(S_One, S_Zero, S_Zero) else list(.neg(y), x, S_Zero)
  w2 <- cross(Matrix(list(x, y, z)), Matrix(w1))$data
  if (circle) {
    nr <- function(w) { s <- .sqrt(.add(lapply(w, function(i) .pow(i, S_Two)))); lapply(w, function(i) .div(i, s)) }
    w1 <- nr(w1); w2 <- nr(w2)
    return(Point3D(lapply(1:3, function(i) .add(list(pc[[i]], .mul(list(w1[[i]], cos(u))), .mul(list(w2[[i]], sin(u))))))))
  }
  Point3D(lapply(1:3, function(i) .add(list(pc[[i]], .mul(list(w1[[i]], u)), .mul(list(w2[[i]], v))))))
}

.mg_planes_concurrent <- function(planes) {
  for (i in planes) if (!inherits(i, "GeoPlane"))
    .mg_err("ValueError", "All objects should be Planes but got ", if (inherits(i, "GeoEntity")) sprintf("<class 'sympy.geometry.%s.%s'>", .mg_pymod(i), i$gtype) else class(i)[1L])
  planes <- planes[!duplicated(vapply(planes, .mg_geo_key, ""))]
  if (length(planes) < 2L) return(FALSE)
  first <- planes[[1L]]; planes <- planes[-1L]
  sol <- .mg_plane_intersection(first, planes[[1L]])
  if (!length(sol)) return(FALSE)
  line <- sol[[1L]]
  for (i in planes[-1L]) {
    l <- .mg_plane_intersection(first, i)
    if (!length(l) || !.geo_contains(line, l[[1L]])) return(FALSE)
  }
  TRUE
}

.mg_plane_distance <- function(e, o) {
  if (!inherits(o, "GeoEntity")) o <- Point(o, dim = 3L)
  n <- Point3D(e$normal_vector)
  N <- .sqrt(.mg_pdot(n, n))
  if (inherits(o, "GeoPoint3D")) return(.div(Abs(.mg_pdot(n, .mg_pt_sub(o, e$p1))), N))
  if (inherits(o, "GeoLinearEntity3D")) {
    A <- o$p1; B <- o$p2; P <- e$p1
    a <- .mg_pdot(n, .mg_pt_sub(A, P)); c <- .mg_pdot(n, .mg_pt_sub(B, A))
    if (inherits(o, "GeoSegment")) {
      b <- .mg_pdot(n, .mg_pt_sub(B, P))
      return(.div(.sub(.add(list(Abs(a), Abs(b))), Abs(c)), .mul(list(S_Two, N))))
    }
    if (inherits(o, "GeoRay")) return(.div(.mul(list(Abs(a), .add(list(S_One, Min(S_Zero, sign(.mul(list(a, c)))))))), N))
    return(.div(.mul(list(Abs(a), .sub(S_One, Abs(sign(c))))), N))
  }
  if (inherits(o, "GeoPlane")) {
    m <- Point3D(o$normal_vector)
    a <- .mg_pdot(n, .mg_pt_sub(o$p1, e$p1))
    delta <- .sub(.mul(list(.mg_pdot(n, n), .mg_pdot(m, m))), .pow(.mg_pdot(n, m), S_Two))
    return(.div(.mul(list(Abs(a), .sub(S_One, sign(delta)))), N))
  }
  .mg_err("NotImplementedError", "distance between Plane and ", o$gtype, " is not implemented")
}

.mg_plane_equals <- function(e, o) {
  if (!inherits(o, "GeoPlane")) return(FALSE)
  r <- cancel(.div(.mg_plane_equation(e), .mg_plane_equation(o)))
  !length(free_symbols(r))
}

.mg_plane_parallel <- function(e, l) {
  if (inherits(l, "GeoLinearEntity3D"))
    return(is_zero_num(.add(Map(function(i, j) .mul(list(i, j)), direction_ratio(l), e$normal_vector))))
  if (inherits(l, "GeoPlane")) return(all(vapply(cross(Matrix(l$normal_vector), Matrix(e$normal_vector))$data, is_zero_num, TRUE)))
  NULL
}

.mg_plane_perpendicular <- function(e, l) {
  if (inherits(l, "GeoLinearEntity3D"))
    return(all(vapply(cross(Matrix(direction_ratio(l)), Matrix(e$normal_vector))$data, is_zero_num, TRUE)))
  if (inherits(l, "GeoPlane"))
    return(is_zero_num(.add(Map(function(i, j) .mul(list(i, j)), l$normal_vector, e$normal_vector))))
  FALSE
}

.mg_plane_intersection <- function(e, o) {
  if (!inherits(o, "GeoEntity")) o <- Point(o, dim = 3L)
  if (inherits(o, "GeoPoint")) return(if (.mg_plane_has(e, o)) list(o) else list())
  if (inherits(o, "GeoLine")) {
    o <- switch(.mg_kind(o), Segment = Segment3D(o$p1, o$p2), Ray = Ray3D(o$p1, o$p2), Line3D(o$p1, o$p2))
    if (.mg_plane_has(e, o)) return(list(o))
    tt <- Dummy("plane_dummy")
    a <- .mg_line_arbitrary(o, tt)
    n <- Point3D(e$normal_vector)
    c <- solve(.mg_pdot(.mg_pt_sub(a, e$p1), n), tt)
    if (!length(c)) return(list())
    c <- Filter(function(i) !identical(ask_fact(i, "real"), FALSE), c)
    if (length(c) > 1L) c <- Filter(function(i) isTRUE(ask_fact(i, "real")), c)
    if (length(c) != 1L) .mg_err("Undecidable", "not sure which point is real")
    p <- .mg_geo_subs(a, tt, c[[1L]])
    if (!.geo_contains(o, p)) return(list())
    return(list(p))
  }
  if (inherits(o, "GeoPlane")) {
    if (.mg_plane_equals(e, o)) return(list(e))
    if (isTRUE(.mg_plane_parallel(e, o))) return(list())
    xyz <- lapply(c("x", "y", "z"), function(nm) Dummy(nm))
    c <- cross(Matrix(list(e$normal_vector)), Matrix(list(o$normal_vector)))$data
    d <- .mg_plane_equation(e, xyz[[1L]], xyz[[2L]], xyz[[3L]])
    f <- .mg_plane_equation(o, xyz[[1L]], xyz[[2L]], xyz[[3L]])
    res <- linsolve(list(d, f), xyz[[1L]], xyz[[2L]], xyz[[3L]])
    r <- if (is_sym(res)) res$args[[1L]] else res[[1L]]
    r <- if (is_sym(r)) r$args else r
    for (s in xyz) r <- lapply(r, function(z) subs(z, s, S_Zero))
    return(list(Line3D(Point3D(r), direction_ratio = c)))
  }
  .mg_isect1(o, e)
}

#' @rdname Plane
#' @export
is_coplanar <- function(e, o) {
  if (inherits(o, "GeoPlane")) {
    xyz <- .mg_plane_xyz()
    r <- cancel(.div(do.call(.mg_plane_equation, c(list(e), xyz)), do.call(.mg_plane_equation, c(list(o), xyz))))
    return(!any(vapply(xyz, function(s) .has1(r, s), TRUE)))
  }
  if (inherits(o, "GeoPoint3D")) return(.mg_plane_has(e, o))
  if (inherits(o, "GeoLinearEntity3D")) return(.mg_plane_has(e, o$p1) && .mg_plane_has(e, o$p2))
  if (inherits(o, "GeoEntity")) return(all(vapply(e$normal_vector[1:2], is_zero_num, TRUE)))
  NULL
}

#' @rdname Plane
#' @export
parallel_plane <- function(pl, pt) Plane(pt, normal_vector = pl$normal_vector)

#' @rdname Plane
#' @export
perpendicular_plane <- function(pl, ...) {
  pts <- list(...)
  if (length(pts) > 2L) .mg_err("ValueError", "No more than 2 pts should be provided.")
  if (!length(pts)) pts <- list(pl$p1)
  n <- pl$normal_vector
  xy0 <- is_zero_num(n[[1L]]) && is_zero_num(n[[2L]])
  if (length(pts) == 1L) {
    dir <- if (xy0) Point3D(0, 1, 0) else Point3D(0, 0, 1)
    pts[[2L]] <- .mg_pt_add(.mg_as_pt(pts[[1L]], 3L), dir)
  }
  p1 <- .mg_as_pt(pts[[1L]], 3L); p2 <- .mg_as_pt(pts[[2L]], 3L)
  l <- Line3D(p1, p2)
  nl <- Line3D(p1, direction_ratio = n)
  if (.mg_line_contains(nl, l)) {
    p3 <- if (xy0) Point3D(0, 1, 0) else Point3D(0, 0, 1)
    if (.mg_line_contains(l, p3)) p3 <- .mg_pt_mul(p3, 2)
  } else p3 <- .mg_pt_add(p1, Point3D(n))
  Plane(p1, p2, p3)
}

#' @rdname Plane
#' @export
projection_line <- function(pl, line) {
  if (!inherits(line, "GeoLine")) .mg_err("NotImplementedError", "Enter a linear entity only")
  a <- .mg_plane_projection(pl, line$p1); b <- .mg_plane_projection(pl, line$p2)
  if (.mg_geo_eq(a, b)) return(a)
  switch(.mg_kind(line), Line = Line3D(a, b), Ray = Ray3D(a, b), Segment3D(a, b))
}

.mg_plane_projection <- function(e, pt) {
  rv <- Point(pt, dim = 3L)
  if (.mg_plane_has(e, rv)) return(rv)
  n <- Point3D(e$normal_vector)
  d <- .div(.mg_pdot(.mg_pt_sub(rv, e$p1), n), .mg_pdot(n, n))
  .mg_pt_sub(rv, .mg_pt_mul(n, d))
}

# ---------------------------------------------------------------------------
# Curves

#' Curves and parabolas
#'
#' `Curve()` is a parametric curve `list(x(t), y(t))` over
#' `list(t, a, b)`; `Parabola()` the parabola with a focus and a directrix.
#' The fields are `$functions`, `$limits` and `$parameter`, and `$focus`
#' and `$directrix`. The functions are SymPy's methods: `focal_length()`,
#' `p_parameter()`, `vertex()`, `axis_of_symmetry()`, and the shared
#' [equation()], [intersection()], [arbitrary_point()] and [geo_length()].
#'
#' @param functions A list of the coordinate functions.
#' @param limits A list `list(t, a, b)`.
#' @param focus The focus.
#' @param directrix The directrix, a line.
#' @param e A curve or a parabola.
#' @return An entity, an expression or a list.
#' @examples
#' t <- Symbol("t")
#' C <- Curve(list(t, t^2), list(t, 0, 1))
#' arbitrary_point(C)
#' P <- Parabola(Point(0, 0), Line(Point(5, 8), Point(7, 8)))
#' vertex(P)
#' equation(P)
#' @export
Curve <- function(functions, limits) {
  f <- .mg_seq(functions); l <- .mg_seq(limits)
  if (is.null(f) || length(f) != 2L)
    .mg_err("ValueError", "Function argument should be (x(t), y(t)) but got ", .mg_tstr(functions))
  if (is.null(l) || length(l) != 3L)
    .mg_err("ValueError", "Limit argument should be (t, tmin, tmax) but got ", .mg_tstr(limits))
  f <- lapply(f, as_sym); l <- lapply(l, as_sym)
  .mg_gnew("Curve", "GeoCurve", list(functions = f, limits = l, parameter = l[[1L]]))
}

# a sequence shown as a Python tuple
.mg_tstr <- function(a) {
  sq <- .mg_seq(a)
  if (is.null(sq)) return(.mg_gstr(a))
  paste0("(", paste(vapply(sq, .mg_gstr, ""), collapse = ", "), if (length(sq) == 1L) ",", ")")
}

.mg_curve_subs <- function(e, a) {
  if (length(a) >= 2L && is_sym(a[[1L]]) && identical(a[[1L]]$key, e$parameter$key))
    return(Point(lapply(e$functions, function(f) subs(f, a[[1L]], a[[2L]]))))
  NULL
}

.mg_curve_arbitrary <- function(e, parameter = "t") {
  if (is.null(parameter)) return(Point(e$functions))
  tnew <- if (is_sym(parameter)) parameter else if (identical(parameter, e$parameter$name)) e$parameter else Symbol(parameter, real = TRUE)
  t <- e$parameter
  if (!identical(tnew$name, t$name) && tnew$name %in% vapply(.mg_geo_free(e), function(s) s$name, ""))
    .mg_err("ValueError", "Symbol ", tnew$name, " already appears in object and cannot be used as a parameter.")
  Point(lapply(e$functions, function(w) subs(w, t, tnew)))
}

.mg_curve_rotate <- function(e, angle = 0, pt = NULL) {
  pt <- if (!is.null(pt)) .mg_pt_neg(.mg_as_pt(pt, 2L)) else Point(0, 0)
  rv <- .mg_entity_translate(e, pt$coords[[1L]], pt$coords[[2L]])
  f <- .new_mat(c(rv$functions, list(S_Zero)), 1L, 3L)
  f <- .matmul(f, rot_axis3(angle))
  rv <- Curve(f$data[1:2], e$limits)
  pt <- .mg_pt_neg(pt)
  .mg_entity_translate(rv, pt$coords[[1L]], pt$coords[[2L]])
}

# ---------------------------------------------------------------------------
# Parabolas

#' @rdname Curve
#' @export
Parabola <- function(focus = NULL, directrix = NULL) {
  focus <- if (!is.null(focus)) .mg_as_pt(focus, 2L) else Point(0, 0)
  directrix <- Line(directrix)
  if (.mg_line_contains(directrix, focus)) .mg_err("ValueError", "The focus must not be a point of directrix")
  .mg_gnew("Parabola", "GeoParabola", list(focus = focus, directrix = directrix))
}

#' @rdname Curve
#' @export
axis_of_symmetry <- function(e) perpendicular_line(e$directrix, e$focus)

#' @rdname Curve
#' @export
focal_length <- function(e) .div(.mg_line_distance(e$directrix, e$focus), S_Two)

#' @rdname Curve
#' @export
p_parameter <- function(e) {
  m <- slope(e$directrix)
  if (identical(m$key, S_Infinity$key)) p <- sign(.add(list(e$focus$coords[[1L]], geo_coefficients(e$directrix)[[3L]])))
  else if (is_zero_num(m)) p <- sign(.add(list(e$focus$coords[[2L]], geo_coefficients(e$directrix)[[3L]])))
  else {
    d <- .mg_line_projection(e$directrix, e$focus)
    p <- sign(.sub(e$focus$coords[[1L]], d$coords[[1L]]))
  }
  .mul(list(p, focal_length(e)))
}

#' @rdname Curve
#' @export
vertex <- function(e) {
  f <- e$focus$coords
  m <- slope(e$directrix)
  if (identical(m$key, S_Infinity$key)) return(Point(.sub(f[[1L]], p_parameter(e)), f[[2L]]))
  if (is_zero_num(m)) return(Point(f[[1L]], .sub(f[[2L]], p_parameter(e))))
  .mg_isect1(axis_of_symmetry(e), e)[[1L]]
}

.mg_parabola_equation <- function(e, x = "x", y = "y") {
  x <- .mg_param(x, "x"); y <- .mg_param(y, "y")
  m <- slope(e$directrix)
  if (identical(m$key, S_Infinity$key)) {
    v <- vertex(e)
    t1 <- .mul(list(.int(4), p_parameter(e), .sub(x, v$coords[[1L]])))
    t2 <- .pow(.sub(y, v$coords[[2L]]), S_Two)
  } else if (is_zero_num(m)) {
    v <- vertex(e)
    t1 <- .mul(list(.int(4), p_parameter(e), .sub(y, v$coords[[2L]])))
    t2 <- .pow(.sub(x, v$coords[[1L]]), S_Two)
  } else {
    a <- e$focus$coords[[1L]]; b <- e$focus$coords[[2L]]
    cf <- geo_coefficients(e$directrix)
    t1 <- .add(list(.pow(.sub(x, a), S_Two), .pow(.sub(y, b), S_Two)))
    t2 <- .div(.pow(.mg_line_equation(e$directrix, x, y), S_Two), .add(list(.pow(cf[[1L]], S_Two), .pow(cf[[2L]], S_Two))))
  }
  .sub(t1, t2)
}

.mg_parabola_intersection <- function(e, o) {
  x <- Symbol("x", real = TRUE); y <- Symbol("y", real = TRUE)
  pe <- .mg_parabola_equation(e)
  if (inherits(o, "GeoParabola")) {
    if (.mg_geo_eq(o, e)) return(list(o))
    return(.mg_sort_ents(.mg_solve_points(list(pe, .mg_parabola_equation(o)), list(x, y))))
  }
  if (inherits(o, "GeoPoint2D")) {
    if (is_zero_num(simplify(subs(pe, list(list(x, o$coords[[1L]]), list(y, o$coords[[2L]])))))) return(list(o))
    return(list())
  }
  if (inherits(o, c("GeoSegment2D", "GeoRay2D"))) {
    pts <- .mg_solve_points(list(pe, .mg_line_equation(Line2D(o$p1, o$p2))), list(x, y))
    return(.mg_sort_ents(Filter(function(p) .mg_line_contains(o, p), pts)))
  }
  if (inherits(o, c("GeoLine2D", "GeoEllipse"))) {
    oe <- if (inherits(o, "GeoEllipse")) .mg_ellipse_equation(o) else .mg_line_equation(o)
    return(.mg_sort_ents(.mg_solve_points(list(pe, oe), list(x, y))))
  }
  if (inherits(o, "GeoLinearEntity3D")) .mg_err("TypeError", "Entity must be two dimensional, not three dimensional")
  .mg_err("TypeError", "Wrong type of argument were put")
}

# ---------------------------------------------------------------------------
# Transformations

.mg_scale_about <- function(e, x, y, pt) {
  pt <- .mg_as_pt(pt, 2L)
  r <- .mg_entity_translate(e, .neg(pt$coords[[1L]]), .neg(pt$coords[[2L]]))
  r <- .mg_entity_scale(r, x, y)
  .mg_entity_translate(r, pt$coords[[1L]], pt$coords[[2L]])
}

.mg_entity_rotate <- function(e, angle, pt = NULL) {
  angle <- as_sym(angle)
  if (inherits(e, "GeoPoint")) {
    if (!inherits(e, "GeoPoint2D")) .mg_err("AttributeError", "'", e$gtype, "' object has no attribute 'rotate'")
    return(.mg_pt_rotate(e, angle, pt))
  }
  if (inherits(e, "GeoEllipse")) return(.mg_ellipse_rotate(e, angle, pt))
  if (inherits(e, "GeoCurve")) return(.mg_curve_rotate(e, angle, pt))
  if (inherits(e, "GeoRegularPolygon"))
    return(RegularPolygon(.mg_pt_rotate(e$center, angle, pt), e$radius, e$n, .add(list(e$rot, angle))))
  .mg_geo_map(e, function(p) .mg_entity_rotate(p, angle, pt))
}

.mg_entity_scale <- function(e, x = 1, y = 1, pt = NULL, z = 1) {
  if (inherits(e, "GeoPoint")) return(.mg_pt_scale(e, x, y, z, pt))
  if (inherits(e, "GeoEllipse")) return(.mg_ellipse_scale(e, x, y, pt))
  if (!is.null(pt)) return(.mg_scale_about(e, x, y, pt))
  x <- as_sym(x); y <- as_sym(y)
  if (inherits(e, "GeoCurve")) return(Curve(list(.mul(list(e$functions[[1L]], x)), .mul(list(e$functions[[2L]], y))), e$limits))
  if (inherits(e, "GeoRegularPolygon")) {
    if (!identical(x$key, y$key)) return(.mg_entity_scale(do.call(Polygon, e$vertices), x, y))
    return(RegularPolygon(e$center, .mul(list(e$radius, x)), e$n, e$rot))
  }
  if (inherits(e, "GeoPlane")) .mg_err("AttributeError", "'Tuple' object has no attribute 'scale'")
  .mg_geo_map(e, function(p) .mg_pt_scale(p, x, y))
}

.mg_entity_translate <- function(e, x = 0, y = 0, z = 0) {
  if (inherits(e, "GeoPoint")) return(.mg_pt_translate(e, x, y, z))
  if (inherits(e, "GeoCurve"))
    return(Curve(list(.add(list(e$functions[[1L]], as_sym(x))), .add(list(e$functions[[2L]], as_sym(y)))), e$limits))
  .mg_geo_map(e, function(p) .mg_pt_translate(p, x, y, z))
}

# GeometryEntity.reflect
.mg_reflect <- function(g, line) {
  if (inherits(g, "GeoEllipse")) return(.mg_ellipse_reflect(g, line))
  if (inherits(g, "GeoRegularPolygon")) {
    v <- g$vertices[[1L]]
    d <- .mg_pt_sub(v, g$center)
    cc <- .mg_reflect(g$center, line); vv <- .mg_reflect(v, line)
    dd <- .mg_pt_sub(vv, cc)
    ang <- closing_angle(Ray(list(0, 0), dd), Ray(list(0, 0), d))
    return(RegularPolygon(cc, .neg(g$radius), g$n, .add(list(g$rot, ang))))
  }
  s <- slope(line)
  pts <- .mg_geo_points(g)
  if (isTRUE(ask_fact(s, "zero"))) {
    v <- line$p1$coords[[2L]]
    if (!.mg_truthy(v)) return(.mg_entity_scale(g, 1, -1))
    reps <- lapply(pts, function(p) .mg_pt_translate(p, 0, .mul(list(S_Two, .sub(v, p$coords[[2L]])))))
  } else if (identical(s$key, S_Infinity$key)) {
    v <- line$p1$coords[[1L]]
    if (!.mg_truthy(v)) return(.mg_entity_scale(g, -1, 1))
    reps <- lapply(pts, function(p) .mg_pt_translate(p, .mul(list(S_Two, .sub(v, p$coords[[1L]])))))
  } else {
    a <- atan(s)
    cf <- geo_coefficients(line)
    d <- .div(.neg(cf[[3L]]), cf[[2L]])
    xx <- Dummy("entity_dummy"); yy <- Dummy("entity_dummy")
    o <- Point(0, 0)
    xf <- .mg_pt_translate(Point(xx, yy), 0, .neg(d))
    xf <- .mg_pt_rotate(xf, .neg(a), o)
    xf <- .mg_pt_scale(xf, 1, -1)
    xf <- .mg_pt_rotate(xf, a, o)
    xf <- .mg_pt_translate(xf, 0, d)
    reps <- lapply(pts, function(p) Point(lapply(xf$coords, function(z) xreplace(z, list(list(xx, p$coords[[1L]]), list(yy, p$coords[[2L]]))))))
  }
  keys <- vapply(pts, .mg_pkey, "")
  .mg_geo_map(g, function(p) { i <- match(.mg_pkey(p), keys); if (is.na(i)) p else reps[[i]] })
}

#' Transformations of geometric entities
#'
#' `geo_rotate()`, `geo_scale()`, `geo_reflect()` and `translate()` are
#' SymPy's `rotate`, `scale`, `reflect` and `translate` methods of points,
#' lines, polygons, ellipses and curves; `geo_transform()` applies a
#' transformation matrix to a point. `geo_rotation_matrix()`,
#' `geo_scale_matrix()` and `geo_translation_matrix()` are the 3 by 3
#' matrices of `sympy.geometry.entity.rotate`, `scale` and `translate`.
#'
#' @param e An entity.
#' @param angle The rotation angle.
#' @param pt The center of the rotation or scaling.
#' @param x,y,z The scale factors or translation offsets.
#' @param line The mirror line.
#' @return The transformed entity.
#' @examples
#' geo_rotate(Point(1, 0), pi / 2)
#' geo_scale(Circle(Point(0, 0), 1), 2, 3)
#' geo_reflect(Point(1, 2), Line(Point(0, 0), Point(1, 1)))
#' translate(Point(1, 2), 3, 4)
#' geo_transform(Point(1, 0), geo_rotation_matrix(pi / 2))
#' @name geometry-transformations
NULL

#' @rdname geometry-transformations
#' @param th The angle of the rotation matrix.
#' @export
geo_rotation_matrix <- function(th) {
  th <- as_sym(th)
  s <- sin(th); c <- cos(th)
  Matrix(list(list(c, s, 0), list(.neg(s), c, 0), list(0, 0, 1)))
}

#' @rdname geometry-transformations
#' @export
geo_translation_matrix <- function(x, y) Matrix(list(list(1, 0, 0), list(0, 1, 0), list(x, y, 1)))

#' @rdname geometry-transformations
#' @export
geo_scale_matrix <- function(x, y, pt = NULL) {
  rv <- Matrix(list(list(x, 0, 0), list(0, y, 0), list(0, 0, 1)))
  if (is.null(pt)) return(rv)
  pt <- Point(pt, dim = 2L)
  geo_translation_matrix(.neg(pt$coords[[1L]]), .neg(pt$coords[[2L]])) %*% rv %*%
    geo_translation_matrix(pt$coords[[1L]], pt$coords[[2L]])
}

#' @rdname geometry-transformations
#' @export
geo_rotate <- function(e, angle, pt = NULL) .mg_entity_rotate(e, angle, pt)

#' @rdname geometry-transformations
#' @export
geo_scale <- function(e, x = 1, y = 1, pt = NULL, z = 1) {
  if (ambient_dimension(e) == 3L && !is.null(pt) && is.null(.mg_seq(pt))) {
    zz <- pt; pt <- if (identical(z, 1)) NULL else z; z <- zz
  }
  .mg_entity_scale(e, x, y, pt, z)
}

#' @rdname geometry-transformations
#' @export
geo_reflect <- function(e, line) .mg_reflect(e, line)

#' @rdname geometry-transformations
#' @export
translate <- function(e, x = 0, y = 0, z = 0) .mg_entity_translate(e, x, y, z)

# ---------------------------------------------------------------------------
# Methods shared by several entities

#' Measurements and relations of geometric entities
#'
#' These functions are SymPy's entity methods and properties that several
#' classes share: `distance()`, `intersection()` (of any number of entities,
#' as `sympy.geometry.intersection`), `area()`, `centroid()` (of points,
#' segments or polygons, as `sympy.geometry.centroid`), `geo_length()`
#' (SymPy's `length`), `equation()`, `projection()`, `angle_between()`,
#' `is_parallel()`, `is_perpendicular()`, `normal_vector()`,
#' `arbitrary_point()`, `parameter_value()`, `plot_interval()`,
#' `random_point()`, `encloses()`, `equals()`, `are_similar()`,
#' `are_coplanar()`, `convex_hull()`, `closest_points()`,
#' `farthest_points()`, `line_equation()` and `circle_equation()`.
#'
#' @param a,b,e,e1,e2,p,q,o Entities.
#' @param ... Entities or points.
#' @param pairwise Intersect every pair of the entities.
#' @param x,y,z Names or symbols of the variables.
#' @param pt A point.
#' @param slope The slope of the major axis of a rotated ellipse.
#' @param scalar For a vector projection, return the scalar component.
#' @param parameter Parameter name or symbol.
#' @param v The second parameter of a plane.
#' @param other A point.
#' @param t The parameter symbol.
#' @param seed A seed for the random point.
#' @param polygon Return a polygon (or segment or point) rather than the
#'   upper and lower hulls.
#' @return An entity, an expression, a logical value or a list.
#' @examples
#' intersection(Circle(Point(0, 0), 1), Line(Point(0, 0), Point(1, 0)))
#' area(Circle(Point(0, 0), 2))
#' centroid(Point(0, 0), Point(2, 0), Point(1, 3))
#' convex_hull(Point(0, 0), Point(1, 1), Point(2, 0), Point(1, 0.5))
#' @name geometry-shared
NULL

#' @rdname geometry-shared
#' @export
distance <- function(p, q) {
  if (inherits(p, "GeoPoint")) {
    if (!inherits(q, "GeoEntity"))
      q <- tryCatch(Point(q, dim = length(p$coords)), TypeError = function(err)
        .mg_err("TypeError", "not recognized as a GeometricEntity: <class '", .mg_pytype(q), "'>"))
    if (inherits(q, "GeoPoint")) return(.mg_pt_distance(p, q))
    return(distance(q, p))
  }
  if (inherits(p, "GeoLine")) return(.mg_line_distance(p, q))
  if (inherits(p, "GeoPlane")) return(.mg_plane_distance(p, q))
  if (inherits(p, "GeoPolygon")) return(.mg_poly_distance(p, q))
  if (!inherits(p, "GeoEntity")) return(distance(Point(p), q))
  .mg_err("TypeError", "distance between ", p$gtype, " and ", if (inherits(q, "GeoEntity")) q$gtype else class(q)[1L], " is not defined")
}

.mg_isect1 <- function(a, b) {
  if (inherits(a, "GeoPoint")) return(.mg_pt_intersection(a, b))
  if (inherits(a, "GeoLine")) return(.mg_line_intersection(a, b))
  if (inherits(a, "GeoEllipse")) return(.mg_ellipse_intersection(a, b))
  if (inherits(a, "GeoPolygon")) return(.mg_poly_intersection(a, b))
  if (inherits(a, "GeoPlane")) return(.mg_plane_intersection(a, b))
  if (inherits(a, "GeoParabola")) return(.mg_parabola_intersection(a, b))
  .mg_err("NotImplementedError", "")
}

#' @rdname geometry-shared
#' @export
intersection <- function(..., pairwise = FALSE) {
  ents <- list(...)
  if (length(ents) >= 1L && .is_ab(ents[[1L]])) return(.ab_intersection(ents[[1L]], as_sym(ents[[2L]])))
  if (length(ents) <= 1L) return(list())
  prec <- NULL
  ents <- lapply(ents, function(e) {
    if (!inherits(e, "GeoEntity")) e <- Point(e)
    fl <- .mg_float_atoms(e)
    if (length(fl)) {
      prec <<- min(c(prec, vapply(fl, .float_prec, 0)))
      e <- .mg_geo_map(e, NULL, .mg_rationalize)
    }
    e
  })
  if (!isTRUE(pairwise)) {
    res <- .mg_isect1(ents[[1L]], ents[[2L]])
    if (inherits(res, "GeoEntity")) res <- list(res)
    for (en in ents[-(1:2)]) {
      nr <- list()
      for (r in res) { s <- .mg_isect1(r, en); nr <- c(nr, if (inherits(s, "GeoEntity")) list(s) else s) }
      res <- nr
    }
  } else {
    ans <- list()
    for (j in seq_along(ents)) for (k in seq_along(ents)) if (k > j) ans <- c(ans, intersection(ents[[j]], ents[[k]]))
    res <- .mg_sort_ents(ans[!duplicated(vapply(ans, .mg_geo_key, ""))])
  }
  if (!is.null(prec)) res <- lapply(res, .mg_geo_evalf, digits = .prec_to_dps(prec))
  res
}

.mg_float_atoms <- function(e) {
  if (inherits(e, "GeoEntity")) return(do.call(c, c(list(list()), lapply(.mg_geo_args(e), .mg_float_atoms))))
  if (is_sym(e)) return(atoms(e, "Float"))
  list()
}

#' @rdname geometry-shared
#' @export
area <- function(e) {
  if (inherits(e, "GeoEllipse")) return(.mg_ellipse_area(e))
  if (inherits(e, "GeoPolygon")) return(.mg_poly_area(e))
  .mg_err("AttributeError", "'", if (inherits(e, "GeoEntity")) e$gtype else class(e)[1L], "' object has no attribute 'area'")
}

#' @rdname geometry-shared
#' @export
centroid <- function(...) {
  args <- list(...)
  if (length(args) == 1L && inherits(args[[1L]], "GeoPolygon")) return(.mg_poly_centroid(args[[1L]]))
  if (!length(args)) return(NULL)
  if (all(vapply(args, inherits, TRUE, "GeoPoint"))) {
    c <- Point(0, 0)
    for (g in args) c <- .mg_pt_add(c, g)
    den <- .int(length(args))
  } else if (all(vapply(args, inherits, TRUE, "GeoSegment"))) {
    c <- Point(0, 0); L <- list()
    for (g in args) {
      l <- geo_length(g)
      c <- .mg_pt_add(c, .mg_pt_mul(midpoint(g$p1, g$p2), l))
      L[[length(L) + 1L]] <- l
    }
    den <- .add(L)
  } else if (all(vapply(args, inherits, TRUE, "GeoPolygon"))) {
    c <- Point(0, 0); A <- list()
    for (g in args) {
      a <- .mg_poly_area(g)
      c <- .mg_pt_add(c, .mg_pt_mul(.mg_poly_centroid(g), a))
      A[[length(A) + 1L]] <- a
    }
    den <- .add(A)
  } else return(NULL)
  c <- .mg_pt_div(c, den)
  Point(lapply(c$coords, simplify))
}

#' @rdname geometry-shared
#' @export
geo_length <- function(e) {
  if (inherits(e, "GeoPoint")) return(S_Zero)
  if (inherits(e, "GeoSegment")) return(.mg_pt_distance(e$p1, e$p2))
  if (inherits(e, "GeoLine")) return(S_Infinity)
  if (inherits(e, "GeoRegularPolygon")) return(.mg_regpoly_side(e))
  if (inherits(e, "GeoCurve")) {
    v <- e$limits[[1L]]
    sp <- .sqrt(.add(lapply(e$functions, function(f) .pow(.diff1(f, v), S_Two))))
    return(sym_integrate(sp, list(v, e$limits[[2L]], e$limits[[3L]])))
  }
  .mg_err("AttributeError", "'", if (inherits(e, "GeoEntity")) e$gtype else class(e)[1L], "' object has no attribute 'length'")
}

#' @rdname geometry-shared
#' @export
equation <- function(e, x = NULL, y = NULL, z = NULL, slope = NULL) {
  if (inherits(e, "GeoPlane")) return(.mg_plane_equation(e, x, y, z))
  dx <- function(v, d) if (is.null(v)) d else v
  if (inherits(e, "GeoLine")) return(.mg_line_equation(e, dx(x, "x"), dx(y, "y"), dx(z, "z")))
  if (inherits(e, "GeoEllipse")) return(.mg_ellipse_equation(e, dx(x, "x"), dx(y, "y"), slope))
  if (inherits(e, "GeoParabola")) return(.mg_parabola_equation(e, dx(x, "x"), dx(y, "y")))
  .mg_err("AttributeError", "this entity has no equation")
}

#' @rdname geometry-shared
#' @export
line_equation <- function(e, x = "x", y = "y") .mg_line_equation(e, x, y)

#' @rdname geometry-shared
#' @export
circle_equation <- function(e, x = "x", y = "y") .mg_ellipse_equation(e, x, y)

#' @rdname geometry-shared
#' @export
projection <- function(e, pt, scalar = FALSE) {
  if (.is_vec(e)) return(.vec_projection(e, pt, scalar))
  if (inherits(e, "GeoPlane")) return(.mg_plane_projection(e, pt))
  .mg_line_projection(e, pt)
}

#' @rdname geometry-shared
#' @export
angle_between <- function(e1, e2) {
  if (inherits(e1, "GeoPlane")) return(.mg_plane_angle(e1, e2))
  .mg_line_angle(e1, e2)
}

#' @rdname geometry-shared
#' @export
is_parallel <- function(e1, e2) {
  if (inherits(e1, "GeoPlane")) return(.mg_plane_parallel(e1, e2))
  .mg_is_parallel(e1, e2)
}

#' @rdname geometry-shared
#' @export
is_perpendicular <- function(e1, e2) {
  if (inherits(e1, "GeoPlane")) return(.mg_plane_perpendicular(e1, e2))
  .mg_is_perpendicular(e1, e2)
}

#' @rdname geometry-shared
#' @export
normal_vector <- function(e) new_sym("Tuple", e$normal_vector)

#' @rdname geometry-shared
#' @export
arbitrary_point <- function(e, parameter = "t", v = NULL) {
  if (inherits(e, "GeoLine")) return(.mg_line_arbitrary(e, parameter))
  if (inherits(e, "GeoEllipse")) return(.mg_ellipse_arbitrary(e, parameter))
  if (inherits(e, "GeoPolygon")) {
    pc <- .mg_poly_arbitrary(e, parameter)
    return(new_sym("Piecewise", lapply(pc, function(p) new_sym("ExprCondPair", list(new_sym("Tuple", p[[1L]]$coords), p[[2L]])))))
  }
  if (inherits(e, "GeoPlane")) return(.mg_plane_arbitrary(e, if (identical(parameter, "t")) NULL else parameter, v))
  if (inherits(e, "GeoCurve")) return(.mg_curve_arbitrary(e, parameter))
  .mg_err("AttributeError", "this entity has no arbitrary_point")
}

#' @rdname geometry-shared
#' @export
parameter_value <- function(e, other, t, v = NULL) {
  if (!inherits(other, "GeoEntity")) other <- Point(other, dim = ambient_dimension(e))
  if (!inherits(other, "GeoPoint")) .mg_err("ValueError", "other must be a point")
  if (inherits(e, "GeoPlane")) {
    if (.mg_geo_eq(other, e$p1)) return(other)
    t <- if (is.character(t)) Symbol(t) else t
    if (!is.null(v) && is.character(v)) v <- Symbol(v)
    if (!(is_sym(t) && t$type == "Symbol") || (!is.null(v) && !(is_sym(v) && v$type == "Symbol")))
      .mg_err("ValueError", "expecting 1 or 2 symbols")
    if (is.null(v)) {
      delta <- .mg_pt_sub(.mg_plane_arbitrary(e, t), e$p1)
      eq <- .mg_pt_sub(delta, unit(.mg_pt_sub(other, e$p1)))
      sol <- solve(eq$coords, t, dict = TRUE)
    } else {
      pt <- .mg_plane_arbitrary(e, t, v)
      sol <- solve(.mg_pt_sub(pt, other)$coords, list(t, v), dict = TRUE)
    }
    if (!length(sol)) .mg_err("ValueError", "Given point is not on Plane")
    out <- sol[[1L]]
    attr(out, "dict") <- TRUE
    return(out)
  }
  if (inherits(e, "GeoPolygon")) {
    if (length(.mg_geo_free(other))) .mg_err("NotImplementedError", "non-numeric coordinates")
    T <- Dummy("polygon_dummy", real = TRUE)
    unknown <- FALSE
    for (pc in .mg_poly_arbitrary(e, T)) {
      sol <- solve(.mg_pt_sub(pc[[1L]], other)$coords, T, dict = TRUE)
      if (!length(sol)) next
      val <- sol[[1L]][[1L]]
      if (.is_true(simplify(subs(pc[[2L]], T, val)))) return(.mg_named(t, val))
      unknown <- TRUE
    }
    if (unknown) .mg_err("ValueError", "Given point may not be on ", e$gtype)
    .mg_err("ValueError", "Given point is not on ", e$gtype)
  }
  T <- Dummy("entity_dummy", real = TRUE)
  sol <- solve(.mg_pt_sub(arbitrary_point(e, T), other)$coords, T, dict = TRUE)
  if (!length(sol)) .mg_err("ValueError", "Given point is not on ", e$gtype)
  .mg_named(t, sol[[1L]][[1L]])
}

.mg_named <- function(k, v) {
  out <- list(v)
  names(out) <- .str(as_sym(k))
  attr(out, "dict") <- TRUE
  out
}

#' @rdname geometry-shared
#' @export
plot_interval <- function(e, parameter = "t") {
  t <- .mg_param(parameter)
  if (inherits(e, "GeoInfLine")) return(list(t, .int(-5), .int(5)))
  if (inherits(e, "GeoRay")) return(list(t, S_Zero, .int(10)))
  if (inherits(e, "GeoSegment") || inherits(e, "GeoPolygon")) return(list(t, S_Zero, S_One))
  if (inherits(e, "GeoEllipse")) return(list(t, .neg(S_Pi), S_Pi))
  if (inherits(e, "GeoCurve")) {
    t <- if (is_sym(parameter)) parameter else if (identical(parameter, e$parameter$name)) e$parameter else Symbol(parameter, real = TRUE)
    return(c(list(t), e$limits[-1L]))
  }
  .mg_err("AttributeError", "this entity has no plot_interval")
}

#' @rdname geometry-shared
#' @export
random_point <- function(e, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (inherits(e, "GeoLine")) {
    tt <- Dummy("line_dummy")
    v <- if (inherits(e, "GeoRay")) abs(stats::rnorm(1)) else if (inherits(e, "GeoSegment")) stats::runif(1) else stats::rnorm(1)
    return(.mg_geo_subs(.mg_line_arbitrary(e, tt), tt, nsimplify(Float(v), rational = TRUE)))
  }
  if (inherits(e, "GeoEllipse")) {
    r <- nsimplify(Float(stats::runif(1)), rational = TRUE)
    c <- .sub(.mul(list(S_Two, r)), S_One)
    s <- .sqrt(.sub(S_One, .pow(c, S_Two)))
    return(Point(.add(list(e$center$coords[[1L]], .mul(list(e$hradius, c)))), .add(list(e$center$coords[[2L]], .mul(list(e$vradius, s))))))
  }
  if (inherits(e, "GeoPlane")) {
    u <- nsimplify(Float(2 * stats::rnorm(1) - 1), rational = TRUE)
    w <- nsimplify(Float(2 * stats::rnorm(1) - 1), rational = TRUE)
    uu <- Dummy("u"); vv <- Dummy("v")
    p <- .mg_plane_arbitrary(e, uu, vv)
    return(Point3D(lapply(p$coords, function(z) subs(z, list(list(uu, u), list(vv, w))))))
  }
  .mg_err("AttributeError", "this entity has no random_point")
}

#' @rdname geometry-shared
#' @export
encloses <- function(e, o) {
  enc <- function(p) {
    r <- encloses_point(e, p)
    isTRUE(r)
  }
  if (inherits(o, "GeoPoint")) return(encloses_point(e, o))
  if (inherits(o, "GeoSegment")) return(all(vapply(list(o$p1, o$p2), enc, TRUE)))
  if (inherits(o, c("GeoRay", "GeoInfLine"))) return(FALSE)
  if (inherits(o, "GeoEllipse"))
    return(enc(o$center) && enc(Point(.add(list(o$center$coords[[1L]], o$hradius)), o$center$coords[[2L]])) &&
             !length(.mg_isect1(e, o)))
  if (inherits(o, "GeoPolygon")) {
    if (inherits(o, "GeoRegularPolygon") && !enc(o$center)) return(FALSE)
    return(all(vapply(o$vertices, enc, TRUE)))
  }
  .mg_err("NotImplementedError", "")
}

#' @rdname geometry-shared
#' @export
equals <- function(a, b) {
  if (inherits(a, "GeoPoint")) {
    if (!inherits(b, "GeoPoint") || length(a$coords) != length(b$coords)) return(FALSE)
    return(all(vapply(Map(.expr_equals, a$coords, b$coords), isTRUE, TRUE)))
  }
  if (inherits(a, "GeoLine")) {
    if (inherits(a, "GeoInfLine")) return(inherits(b, "GeoInfLine") && is_collinear(a$p1, b$p1, a$p2, b$p2))
    if (inherits(a, "GeoRay")) return(inherits(b, "GeoRay") && .mg_geo_eq(a$source, b$source) && .mg_line_contains(a, b$p2))
    if (!identical(class(a)[[1L]], class(b)[[1L]])) return(FALSE)
    ka <- sort(c(.mg_pkey(a$p1), .mg_pkey(a$p2))); kb <- sort(c(.mg_pkey(b$p1), .mg_pkey(b$p2)))
    return(identical(ka, kb))
  }
  if (inherits(a, "GeoPlane")) return(.mg_plane_equals(a, b))
  if (inherits(a, "GeoEntity")) return(.mg_geo_eq(a, b))
  if (inherits(a, "symmatrix")) {
    b <- .as_mat(b)
    if (a$nrow != b$nrow || a$ncol != b$ncol) return(FALSE)
    r <- Map(.expr_equals, a$data, b$data)
    if (any(vapply(r, isFALSE, TRUE))) return(FALSE)
    if (all(vapply(r, isTRUE, TRUE))) return(TRUE)
    return(NULL)
  }
  .mg_eq0(.sub(as_sym(a), as_sym(b)))
}

#' @rdname geometry-shared
#' @export
are_similar <- function(e1, e2) {
  if (inherits(e1, "GeoEntity") && inherits(e2, "GeoEntity") && .mg_geo_eq(e1, e2)) return(TRUE)
  r <- tryCatch(is_similar(e1, e2), NotImplementedError = function(err) NULL)
  if (!is.null(r)) return(r)
  r <- tryCatch(is_similar(e2, e1), NotImplementedError = function(err) NULL)
  if (!is.null(r)) return(r)
  .mg_err("GeometryError", "Cannot test similarity between ", e1$gtype, " and ", e2$gtype)
}

#' @rdname geometry-shared
#' @export
are_coplanar <- function(...) {
  e <- list(...)
  e <- e[!duplicated(vapply(e, function(z) if (inherits(z, "GeoEntity")) .mg_geo_key(z) else paste(z), ""))]
  for (i in seq_along(e)) if (inherits(e[[i]], "GeoPlane")) {
    pl <- e[[i]]
    return(all(vapply(e[-i], function(p) isTRUE(is_coplanar(p, pl)), TRUE)))
  }
  if (length(e) && all(vapply(e, inherits, TRUE, "GeoPoint3D"))) {
    if (length(e) < 3L) return(FALSE)
    a <- e[[length(e)]]; b <- e[[length(e) - 1L]]
    e <- e[seq_len(length(e) - 2L)]
    e <- Filter(function(i) !is_collinear(a, b, i), e)
    if (!length(e)) return(FALSE)
    p <- Plane(a, b, e[[length(e)]])
    for (i in e[-length(e)]) if (!.mg_plane_has(p, i)) return(FALSE)
    return(TRUE)
  }
  pt3d <- list()
  for (i in e) {
    if (inherits(i, "GeoPoint3D")) pt3d <- c(pt3d, list(i))
    else if (inherits(i, "GeoLinearEntity3D")) pt3d <- c(pt3d, list(i$p1, i$p2))
    else if (inherits(i, "GeoEntity")) for (p in .mg_geo_args(i)) if (inherits(p, "GeoPoint")) pt3d <- c(pt3d, list(Point3D(c(p$coords, list(S_Zero)))))
  }
  do.call(are_coplanar, pt3d)
}

# Plane.are_concurrent
.mg_plane_are_concurrent <- function(...) .mg_planes_concurrent(list(...))

.mg_pymod <- function(e) {
  if (inherits(e, "GeoPoint")) return("point")
  if (inherits(e, "GeoLine")) return("line")
  if (inherits(e, "GeoEllipse")) return("ellipse")
  if (inherits(e, "GeoPolygon")) return("polygon")
  if (inherits(e, "GeoPlane")) return("plane")
  if (inherits(e, "GeoCurve")) return("curve")
  "parabola"
}

.mg_find <- function(x, eq) {
  fs <- free_symbols(eq)
  xs <- Filter(function(s) if (is.character(x)) identical(s$name, x) else identical(s$key, as_sym(x)$key), fs)
  if (!length(xs)) .mg_err("ValueError", "could not find ", if (is.character(x)) x else .str(x))
  if (length(xs) != 1L) .mg_err("ValueError", "ambiguous ", if (is.character(x)) x else .str(x))
  xs[[1L]]
}

.mg_pt_sort <- function(p) p[order(vapply(p, .mg_ent_sort_key, ""), method = "radix")]

#' @rdname geometry-shared
#' @export
convex_hull <- function(..., polygon = TRUE) {
  p <- list()
  for (e in list(...)) {
    if (!inherits(e, "GeoEntity")) e <- Point(e)
    if (inherits(e, "GeoPoint")) p <- c(p, list(e))
    else if (inherits(e, "GeoSegment")) p <- c(p, list(e$p1, e$p2))
    else if (inherits(e, "GeoPolygon")) p <- c(p, e$vertices)
    else .mg_err("NotImplementedError", "Convex hull for ", e$gtype, " not implemented.")
  }
  p <- .mg_uniq_pts(p)
  if (any(vapply(p, function(z) length(z$coords) != 2L, TRUE))) .mg_err("ValueError", "Can only compute the convex hull in two dimensions")
  if (length(p) == 1L) return(if (polygon) p[[1L]] else list(p[[1L]], NULL))
  if (length(p) == 2L) { s <- Segment(p[[1L]], p[[2L]]); return(if (polygon) s else list(s, NULL)) }
  orient <- function(a, b, c) .sub(.mul(list(.sub(b$coords[[2L]], a$coords[[2L]]), .sub(c$coords[[1L]], a$coords[[1L]]))),
                                    .mul(list(.sub(b$coords[[1L]], a$coords[[1L]]), .sub(c$coords[[2L]], a$coords[[2L]]))))
  U <- list(); L <- list()
  for (pi in .mg_pt_sort(p)) {
    while (length(U) > 1L && .mg_rbool(Le(orient(U[[length(U) - 1L]], U[[length(U)]], pi), S_Zero))) U[[length(U)]] <- NULL
    while (length(L) > 1L && .mg_rbool(Ge(orient(L[[length(L) - 1L]], L[[length(L)]], pi), S_Zero))) L[[length(L)]] <- NULL
    U <- c(U, list(pi)); L <- c(L, list(pi))
  }
  U <- rev(U)
  hull <- c(L, if (length(U) > 2L) U[2:(length(U) - 1L)] else list())
  if (length(hull) == 2L) { s <- Segment(hull[[1L]], hull[[2L]]); return(if (polygon) s else list(s, NULL)) }
  if (polygon) return(do.call(Polygon, hull))
  list(rev(U), L)
}

.mg_hypot <- function(dx, dy) .sqrt(.add(list(.pow(dx, S_Two), .pow(dy, S_Two))))

#' @rdname geometry-shared
#' @export
closest_points <- function(...) {
  p <- .mg_pt_sort(.mg_uniq_pts(lapply(list(...), function(a) .mg_as_pt(a, 2L))))
  if (length(p) < 2L) .mg_err("ValueError", "At least 2 distinct points must be given.")
  dd <- function(i, j) .mg_hypot(.sub(p[[i]]$coords[[1L]], p[[j]]$coords[[1L]]), .sub(p[[i]]$coords[[2L]], p[[j]]$coords[[2L]]))
  rv <- list(c(1L, 2L))
  best <- dd(1L, 2L)
  left <- 1L
  box <- c(1L, 2L)
  if (length(p) > 2L) for (i in 3:length(p)) {
    while (left < i && .mg_rbool(Gt(.sub(p[[i]]$coords[[1L]], p[[left]]$coords[[1L]]), best))) {
      box <- box[-1L]; left <- left + 1L
    }
    for (j in box) {
      d <- dd(i, j)
      if (.mg_rbool(Lt(d, best))) rv <- list(c(j, i))
      else if (.mg_rbool(Eq(d, best))) rv <- c(rv, list(c(j, i)))
      else next
      best <- d
    }
    box <- c(box, i)
  }
  lapply(rv, function(pr) list(p[[pr[1L]]], p[[pr[2L]]]))
}

#' @rdname geometry-shared
#' @export
farthest_points <- function(...) {
  args <- lapply(list(...), function(a) .mg_as_pt(a, 2L))
  hl <- do.call(convex_hull, c(args, list(polygon = FALSE)))
  U <- hl[[1L]]; L <- hl[[2L]]
  pairs <- list()
  if (is.null(L)) {
    if (inherits(U, "GeoPoint")) .mg_err("ValueError", "At least two distinct points must be given.")
    pairs <- list(list(U$p1, U$p2))
  } else {
    i <- 1L; j <- length(L)
    while (i < length(U) || j > 1L) {
      pairs[[length(pairs) + 1L]] <- list(U[[i]], L[[j]])
      if (i == length(U)) j <- j - 1L
      else if (j == 1L) i <- i + 1L
      else if (.mg_rbool(Gt(.mul(list(.sub(U[[i + 1L]]$coords[[2L]], U[[i]]$coords[[2L]]), .sub(L[[j]]$coords[[1L]], L[[j - 1L]]$coords[[1L]]))),
                           .mul(list(.sub(L[[j]]$coords[[2L]], L[[j - 1L]]$coords[[2L]]), .sub(U[[i + 1L]]$coords[[1L]], U[[i]]$coords[[1L]])))))) i <- i + 1L
      else j <- j - 1L
    }
  }
  rv <- list(); diam <- S_Zero
  for (pr in pairs) {
    pr <- .mg_pt_sort(pr)
    d <- .mg_hypot(.sub(pr[[1L]]$coords[[1L]], pr[[2L]]$coords[[1L]]), .sub(pr[[1L]]$coords[[2L]], pr[[2L]]$coords[[2L]]))
    if (.mg_rbool(Gt(d, diam))) rv <- list(pr)
    else if (.mg_rbool(Eq(d, diam))) rv <- c(rv, list(pr))
    else next
    diam <- d
  }
  rv[!duplicated(vapply(rv, function(pr) paste(.mg_pkey(pr[[1L]]), .mg_pkey(pr[[2L]])), ""))]
}

# __contains__ of an entity (the `in` operator)
.geo_contains <- function(s, x) {
  if (inherits(s, "GeoPoint")) {
    if (length(s$coords) %in% 2:3) return(inherits(x, "GeoPoint") && .mg_geo_eq(s, x))
    return(any(vapply(s$coords, function(z) identical(z$key, as_sym(x)$key), TRUE)))
  }
  if (inherits(s, "GeoLine")) return(.mg_line_contains(s, x))
  if (inherits(s, "GeoEllipse")) return(.mg_ellipse_has(s, x))
  if (inherits(s, "GeoPolygon")) return(.mg_poly_has(s, x))
  if (inherits(s, "GeoPlane")) return(.mg_plane_has(s, x))
  if (inherits(s, "GeoParabola")) return(inherits(x, "GeoParabola") && .mg_geo_eq(s, x))
  if (identical(s$gtype, x$gtype)) return(.mg_geo_eq(s, x))
  .mg_err("NotImplementedError", "")
}
