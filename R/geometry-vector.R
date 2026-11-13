# Geometry and vector calculus, after sympy/geometry and sympy/vector.

.pt_coords <- function(p) if (inherits(p, "GeoPoint")) p$coords else lapply(as.list(p), as_sym)

#' Plane geometry
#'
#' `Point()`, `Line()`, `Segment()`, `Circle()`, `Polygon()` and `Triangle()`
#' build geometric entities with exact coordinates. `distance()` measures
#' between points, `midpoint()`, `area()`, `perimeter()`, `centroid()`,
#' `slope()`, `is_parallel()`, `is_perpendicular()`, `is_collinear()`,
#' `intersection()`, `line_equation()` and `circle_equation()` follow SymPy.
#'
#' @param ... Coordinates for `Point()`; vertices for `Polygon()`.
#' @param p,q,a,b,P Points or entities.
#' @param center,radius Circle data.
#' @param l,m Lines.
#' @param x,y Symbols for equations.
#' @return A geometric entity, expression or list.
#' @examples
#' A <- Point(0, 0); B <- Point(3, 4)
#' distance(A, B)
#' area(Triangle(Point(0, 0), Point(4, 0), Point(0, 3)))
#' intersection(Line(Point(0, 0), Point(1, 1)), Line(Point(0, 2), Point(2, 0)))
#' @export
Point <- function(...) {
  a <- list(...)
  if (length(a) == 1L && (is.list(a[[1L]]) || length(a[[1L]]) > 1L)) a <- as.list(a[[1L]])
  structure(list(coords = lapply(a, as_sym)), class = c("GeoPoint", "GeoEntity"))
}

#' @export
format.GeoPoint <- function(x, ...) {
  n <- length(x$coords)
  nm <- if (n == 2L) "Point2D" else if (n == 3L) "Point3D" else "Point"
  paste0(nm, "(", paste(vapply(x$coords, .str, ""), collapse = ", "), ")")
}

#' @export
print.GeoEntity <- function(x, ...) {
  cat(format(x), "\n")
  invisible(x)
}

#' @export
Ops.GeoPoint <- function(e1, e2) {
  if (.Generic == "==") return(all(Map(function(a, b) is_zero_num(simplify(.sub(a, b))), e1$coords, e2$coords) == TRUE))
  if (.Generic == "!=") return(!(e1 == e2))
  stop("unsupported operator for points")
}

#' @rdname Point
#' @export
distance <- function(p, q) {
  if (inherits(p, "GeoInfLine") && inherits(q, "GeoInfLine") && length(p$p1$coords) == 3L) {
    # parallel lines: a point's distance; skew lines: along the common normal
    cr <- .vcross(.ldir(p), .ldir(q))
    if (all(vapply(cr, .geo_zero, TRUE))) return(distance(q$p1, p))
    return(simplify(.div(fcall("Abs", list(.vdot(.vsub(q$p1$coords, p$p1$coords), cr))), .vnorm(cr))))
  }
  if (inherits(p, c("GeoLine", "GeoPlane")) && inherits(q, "GeoPoint")) return(distance(q, p))
  if (inherits(p, "GeoLine") && inherits(q, "GeoPlane")) return(distance(q, p))
  if (inherits(p, "GeoPlane")) return(.dist_plane(p, q))
  if (inherits(q, "GeoPlane")) return(.dist_plane(q, p))
  if (inherits(q, "GeoSegment")) return(.dist_point_segment(p, q))
  if (inherits(q, "GeoRay")) return(.dist_point_ray(p, q))
  if (inherits(q, "GeoLine")) return(.dist_point_infline(p, q))
  simplify(.sqrt(.add(Map(function(a, b) .sq(.sub(a, b)), p$coords, q$coords))))
}

# Plane.distance in SymPy, for points, segments, rays, lines and planes.
.dist_plane <- function(pl, o) {
  n <- pl$normal
  N <- .sqrt(.vdot(n, n))
  sgn <- function(v) fcall("sign", list(simplify(v)))
  if (inherits(o, "GeoPoint")) return(simplify(.div(fcall("Abs", list(.vdot(n, .vsub(o$coords, pl$p$coords)))), N)))
  if (inherits(o, "GeoPlane")) {
    m <- o$normal
    a <- .vdot(n, .vsub(o$p$coords, pl$p$coords))
    delta <- .sub(.mul(list(.vdot(n, n), .vdot(m, m))), .sq(.vdot(n, m)))
    return(simplify(.div(.mul(list(fcall("Abs", list(a)), .sub(S_One, sgn(delta)))), N)))
  }
  A <- o$p1$coords; B <- o$p2$coords
  a <- .vdot(n, .vsub(A, pl$p$coords))
  c <- .vdot(n, .vsub(B, A))
  if (inherits(o, "GeoSegment")) {
    b <- .vdot(n, .vsub(B, pl$p$coords))
    return(simplify(.div(.sub(.add(list(fcall("Abs", list(a)), fcall("Abs", list(b)))), fcall("Abs", list(c))), .mul(list(S_Two, N)))))
  }
  if (inherits(o, "GeoRay"))
    return(simplify(.div(.mul(list(fcall("Abs", list(a)), .add(list(S_One, Min(S_Zero, sgn(.mul(list(a, c)))))))), N)))
  simplify(.div(.mul(list(fcall("Abs", list(a)), .sub(S_One, fcall("Abs", list(sgn(c)))))), N))
}

# Distance to the infinite line through l, in any dimension.
.dist_point_infline <- function(p, l) {
  d <- length(l$p1$coords)
  if (d == 2L) return(.dist_point_line(p, l))
  if (d == 3L) return(.dist_point_line3d(p, l))
  distance(p, .line_projection(l, p))
}

.line_projection <- function(l, p) {
  d <- .ldir(l)
  t <- simplify(.div(.vdot(.vsub(p$coords, l$p1$coords), d), .vdot(d, d)))
  do.call(Point, Map(function(a, u) simplify(.add(list(a, .mul(list(t, u))))), l$p1$coords, d))
}

# Ray.distance in SymPy: the projection when it lies on the ray, else the source.
.dist_point_ray <- function(p, r) {
  if (.geo_on_line(r, p) && .on_segment(r, p)) return(S_Zero)
  pr <- .line_projection(r, p)
  t <- .line_param(r, pr)
  ok <- Ge(t, S_Zero)
  if (ok$type != "BooleanAtom") stop("cannot decide where the point projects onto the ray", call. = FALSE)
  distance(p, if (.is_true(ok)) pr else r$p1)
}

.geo_on_line <- function(l, p) !is.null(.line_param(l, p))

#' @rdname Point
#' @export
midpoint <- function(p, q = NULL) {
  if (is.null(q) && inherits(p, "GeoSegment")) return(midpoint(p$p1, p$p2))
  do.call(Point, Map(function(a, b) .div(.add(list(a, b)), S_Two), p$coords, q$coords))
}

.mk_line <- function(p, q, cls) {
  if (isTRUE(p == q)) stop("a line needs two distinct points")
  structure(list(p1 = p, p2 = q), class = c(cls, "GeoLine", "GeoEntity"))
}

#' @rdname Point
#' @export
Line <- function(p, q) .mk_line(p, q, "GeoInfLine")

#' @rdname Point
#' @export
Segment <- function(p, q) .mk_line(p, q, "GeoSegment")

#' @export
format.GeoLine <- function(x, ...) {
  nm <- if (inherits(x, "GeoSegment")) "Segment" else if (inherits(x, "GeoRay")) "Ray" else "Line"
  d <- length(x$p1$coords)
  paste0(nm, if (d %in% 2:3) paste0(d, "D"), "(", format(x$p1), ", ", format(x$p2), ")")
}

# Coefficients a, b, c of a*x + b*y + c = 0.
.line_abc <- function(l) {
  p <- l$p1$coords; q <- l$p2$coords
  a <- .sub(p[[2L]], q[[2L]])
  b <- .sub(q[[1L]], p[[1L]])
  c <- .sub(.mul(list(p[[1L]], q[[2L]])), .mul(list(q[[1L]], p[[2L]])))
  list(a = a, b = b, c = c)
}

.dist_point_segment <- function(p, s) {
  d <- Map(.sub, s$p2$coords, s$p1$coords)
  dot <- function(u) simplify(.add(Map(function(a, b) .mul(list(a, b)), d, u)))
  t1 <- Ge(dot(Map(.sub, p$coords, s$p1$coords)), S_Zero)
  t2 <- Le(dot(Map(.sub, p$coords, s$p2$coords)), S_Zero)
  if (t1$type != "BooleanAtom" || t2$type != "BooleanAtom")
    stop("cannot decide where the point projects onto the segment", call. = FALSE)
  s1 <- .is_true(t1)
  s2 <- .is_true(t2)
  if (s1 && s2) return(.dist_point_infline(p, s))
  distance(p, if (s1) s$p2 else s$p1)
}

.dist_point_line <- function(p, l) {
  k <- .line_abc(l)
  simplify(.div(fcall("Abs", list(.add(list(.mul(list(k$a, p$coords[[1L]])), .mul(list(k$b, p$coords[[2L]])), k$c)))),
                .sqrt(.add(list(.sq(k$a), .sq(k$b))))))
}

#' @rdname Point
#' @export
slope <- function(l) {
  p <- l$p1$coords; q <- l$p2$coords
  dx <- simplify(.sub(q[[1L]], p[[1L]]))
  if (is_zero_num(dx)) return(S_Infinity)
  simplify(.div(.sub(q[[2L]], p[[2L]]), dx))
}

#' @rdname Point
#' @export
line_equation <- function(l, x = Symbol("x"), y = Symbol("y")) {
  k <- .line_abc(l)
  expand(.add(list(.mul(list(k$a, x)), .mul(list(k$b, y)), k$c)))
}

#' @rdname Point
#' @export
is_parallel <- function(l, m) {
  if (inherits(l, "GeoPlane") || inherits(m, "GeoPlane") || length(l$p1$coords) == 3L) {
    u <- if (inherits(l, "GeoPlane")) l$normal else .ldir(l)
    v <- if (inherits(m, "GeoPlane")) m$normal else .ldir(m)
    if (xor(inherits(l, "GeoPlane"), inherits(m, "GeoPlane"))) return(.geo_zero(.vdot(u, v)))
    return(all(vapply(.vcross(u, v), .geo_zero, TRUE)))
  }
  a <- .line_abc(l); b <- .line_abc(m)
  is_zero_num(simplify(.sub(.mul(list(a$a, b$b)), .mul(list(a$b, b$a)))))
}

#' @rdname Point
#' @export
is_perpendicular <- function(l, m) {
  a <- .line_abc(l); b <- .line_abc(m)
  is_zero_num(simplify(.add(list(.mul(list(a$a, b$a)), .mul(list(a$b, b$b))))))
}

#' @rdname Point
#' @export
is_collinear <- function(...) {
  pts <- list(...)
  if (length(pts) < 3L) return(TRUE)
  p <- pts[[1L]]$coords; q <- pts[[2L]]$coords
  all(vapply(pts[-(1:2)], function(r) {
    d <- .sub(.mul(list(.sub(q[[1L]], p[[1L]]), .sub(r$coords[[2L]], p[[2L]]))),
              .mul(list(.sub(q[[2L]], p[[2L]]), .sub(r$coords[[1L]], p[[1L]]))))
    is_zero_num(simplify(d))
  }, TRUE))
}

#' @rdname Point
#' @export
Circle <- function(center, radius) {
  structure(list(center = center, radius = as_sym(radius)), class = c("GeoCircle", "GeoEntity"))
}

#' @export
format.GeoCircle <- function(x, ...) paste0("Circle(", format(x$center), ", ", .str(x$radius), ")")

#' @rdname Point
#' @export
circle_equation <- function(l, x = Symbol("x"), y = Symbol("y")) {
  expand(.sub(.add(list(.sq(.sub(x, l$center$coords[[1L]])), .sq(.sub(y, l$center$coords[[2L]])))), .sq(l$radius)))
}

#' @rdname Point
#' @export
Polygon <- function(...) {
  v <- list(...)
  cls <- if (length(v) == 3L) c("GeoTriangle", "GeoPolygon", "GeoEntity") else c("GeoPolygon", "GeoEntity")
  structure(list(vertices = v), class = cls)
}

#' @rdname Point
#' @export
Triangle <- function(a, b, P) Polygon(a, b, P)

#' @export
format.GeoPolygon <- function(x, ...) {
  nm <- if (inherits(x, "GeoTriangle")) "Triangle" else "Polygon"
  paste0(nm, "(", paste(vapply(x$vertices, format, ""), collapse = ", "), ")")
}

.signed_area <- function(v) {
  n <- length(v)
  .div(.add(lapply(seq_len(n), function(i) {
    p <- v[[i]]$coords; q <- v[[if (i == n) 1L else i + 1L]]$coords
    .sub(.mul(list(p[[1L]], q[[2L]])), .mul(list(q[[1L]], p[[2L]])))
  })), S_Two)
}

#' @rdname Point
#' @export
area <- function(p) {
  if (inherits(p, "GeoCircle")) return(.mul(list(S_Pi, .sq(p$radius))))
  if (inherits(p, "GeoEllipse")) return(.mul(list(S_Pi, p$hradius, p$vradius)))
  simplify(.signed_area(p$vertices))
}

#' @rdname Point
#' @export
perimeter <- function(p) {
  if (inherits(p, c("GeoCircle", "GeoEllipse"))) return(circumference(p))
  v <- p$vertices
  n <- length(v)
  simplify(.add(lapply(seq_len(n), function(i) distance(v[[i]], v[[if (i == n) 1L else i + 1L]]))))
}

#' @rdname Point
#' @export
centroid <- function(p) {
  v <- p$vertices
  n <- length(v)
  A <- .signed_area(v)
  coord <- function(k) simplify(.div(.add(lapply(seq_len(n), function(i) {
    a <- v[[i]]$coords; b <- v[[if (i == n) 1L else i + 1L]]$coords
    .mul(list(.add(list(a[[k]], b[[k]])), .sub(.mul(list(a[[1L]], b[[2L]])), .mul(list(b[[1L]], a[[2L]])))))
  })), .mul(list(.int(6), A))))
  Point(coord(1L), coord(2L))
}

.on_segment <- function(s, pt) {
  if (inherits(s, "GeoRay")) {
    t <- .line_param(s, pt)
    return(is.null(t) || !.is_false(Ge(t, 0)))
  }
  if (!inherits(s, "GeoSegment")) return(TRUE)
  in_range <- function(k) {
    a <- s$p1$coords[[k]]; b <- s$p2$coords[[k]]; v <- pt$coords[[k]]
    d1 <- .cmp(.sub(v, a), S_Zero); d2 <- .cmp(.sub(v, b), S_Zero)
    is.na(d1) || is.na(d2) || d1 * d2 <= 0L
  }
  all(vapply(seq_along(s$p1$coords), in_range, TRUE))
}

# SymPy's _span_test: 0 at the start point, 1 ahead of it, -1 behind.
.span_test <- function(l, pt) {
  if (isTRUE(pt == l$p1)) return(0L)
  if (.is_true(Gt(simplify(.vdot(.ldir(l), .vsub(pt$coords, l$p1$coords))), S_Zero))) 1L else -1L
}

# Intersection of two linear entities on one line (LinearEntity.intersection).
.isect_collinear <- function(a, b) {
  # an overlap that shrinks to one point is that point
  Segment <- function(p, q) if (isTRUE(p == q)) p else .mk_line(p, q, "GeoSegment")
  if (inherits(a, "GeoInfLine")) return(list(b))
  if (inherits(b, "GeoInfLine")) return(list(a))
  ray_seg <- function(r, sg) {
    st1 <- .span_test(r, sg$p1); st2 <- .span_test(r, sg$p2)
    if (st1 < 0L && st2 < 0L) return(list())
    if (st1 >= 0L && st2 >= 0L) return(list(sg))
    list(Segment(r$p1, if (st1 >= 0L) sg$p1 else sg$p2))
  }
  if (inherits(a, "GeoRay") && inherits(b, "GeoRay")) {
    if (.is_true(Gt(simplify(.vdot(.ldir(a), .ldir(b))), S_Zero)))
      return(if (.span_test(a, b$p1) >= 0L) list(b) else list(a))
    st <- .span_test(a, b$p1)
    if (st < 0L) return(list())
    if (st == 0L) return(list(b$p1))
    return(list(Segment(a$p1, b$p1)))
  }
  if (inherits(a, "GeoRay")) return(ray_seg(a, b))
  if (inherits(b, "GeoRay")) return(ray_seg(b, a))
  inside <- function(s1, s2) .on_segment(s1, s2$p1) && .on_segment(s1, s2$p2)
  if (inside(a, b)) return(list(b))
  if (inside(b, a)) return(list(a))
  s1 <- a; s2 <- b
  if (.is_true(Lt(simplify(.vdot(.ldir(s1), .ldir(s2))), S_Zero))) s2 <- Segment(s2$p2, s2$p1)
  if (.span_test(s1, s2$p1) < 0L) { tmp <- s1; s1 <- s2; s2 <- tmp }
  if (.span_test(s2, s1$p2) < 0L) return(list())
  list(Segment(s2$p1, s1$p2))
}

#' @rdname Point
#' @export
intersection <- function(a, b) {
  if (.is_ab(a)) return(.ab_intersection(a, as_sym(b)))
  if (inherits(a, "GeoLine") && inherits(b, "GeoPlane")) return(intersection(b, a))
  if (inherits(a, "GeoPlane") && inherits(b, "GeoLine")) return(.isect_plane_line(a, b))
  if (inherits(a, "GeoLine") && inherits(b, c("GeoEllipse", "GeoParabola"))) return(.isect_conic_line(b, a))
  if (inherits(a, c("GeoEllipse", "GeoParabola")) && inherits(b, "GeoLine")) return(.isect_conic_line(a, b))
  if (inherits(a, "GeoLine") && inherits(b, "GeoLine") && length(a$p1$coords) == 3L) return(.isect_lines3d(a, b))
  if (inherits(a, "GeoCircle") && inherits(b, "GeoLine")) return(.isect_circle_line(a, b))
  if (inherits(a, "GeoLine") && inherits(b, "GeoCircle")) return(.isect_circle_line(b, a))
  if (inherits(a, "GeoCircle") && inherits(b, "GeoCircle")) return(.isect_circles(a, b))
  if (inherits(a, "GeoLine") && inherits(b, "GeoLine")) {
    k1 <- .line_abc(a); k2 <- .line_abc(b)
    det <- simplify(.sub(.mul(list(k1$a, k2$b)), .mul(list(k2$a, k1$b))))
    if (is_zero_num(det)) {
      on <- is_zero_num(simplify(.add(list(.mul(list(k2$a, a$p1$coords[[1L]])), .mul(list(k2$b, a$p1$coords[[2L]])), k2$c))))
      return(if (on) .isect_collinear(a, b) else list())
    }
    x <- simplify(.div(.sub(.mul(list(k1$b, k2$c)), .mul(list(k2$b, k1$c))), det))
    y <- simplify(.div(.sub(.mul(list(k2$a, k1$c)), .mul(list(k1$a, k2$c))), det))
    pt <- Point(x, y)
    if (.on_segment(a, pt) && .on_segment(b, pt)) return(list(pt))
    return(list())
  }
  stop("intersection not implemented for these entities")
}

.isect_circle_line <- function(c, l) {
  t <- Dummy("t")
  p <- l$p1$coords; q <- l$p2$coords
  px <- .add(list(p[[1L]], .mul(list(t, .sub(q[[1L]], p[[1L]])))))
  py <- .add(list(p[[2L]], .mul(list(t, .sub(q[[2L]], p[[2L]])))))
  eq <- .sub(.add(list(.sq(.sub(px, c$center$coords[[1L]])), .sq(.sub(py, c$center$coords[[2L]])))), .sq(c$radius))
  ts <- Filter(function(v) !identical(ask_fact(v, "real"), FALSE), roots(expand(eq), t))
  pts <- lapply(ts, function(v) Point(simplify(subs(px, t, v)), simplify(subs(py, t, v))))
  .ordered_points(Filter(function(pt) .on_segment(l, pt), pts))
}

.isect_circles <- function(a, b) {
  x <- Dummy("x"); y <- Dummy("y")
  ea <- circle_equation(a, x, y); eb <- circle_equation(b, x, y)
  sols <- solve(list(ea, eb), list(x, y), dict = TRUE)
  .ordered_points(lapply(sols, function(s) Point(simplify(s[[x$name]]), simplify(s[[y$name]]))))
}

# list(ordered(points)) as Ellipse.intersection returns them: node count,
# then default_sort_key of the coordinates
.ordered_points <- function(pts) {
  .key_sort(lapply(pts, function(p) list(1 + sum(vapply(p$coords, .nodes, 0)), lapply(p$coords, .sort_key))), pts)
}

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
  is_zero_matrix(curl(v, vars))
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
