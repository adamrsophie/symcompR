# 3D points and lines, rays, planes, ellipses, parabolas and curves, after
# sympy/geometry (point.py, line.py, plane.py, ellipse.py, parabola.py,
# curve.py, util.py).

.vsub <- function(a, b) Map(.sub, a, b)
.vdot <- function(a, b) .add(Map(function(u, v) .mul(list(u, v)), a, b))
.vcross <- function(a, b) list(.sub(.mul(list(a[[2L]], b[[3L]])), .mul(list(a[[3L]], b[[2L]]))),
                               .sub(.mul(list(a[[3L]], b[[1L]])), .mul(list(a[[1L]], b[[3L]]))),
                               .sub(.mul(list(a[[1L]], b[[2L]])), .mul(list(a[[2L]], b[[1L]]))))
.vnorm <- function(a) simplify(.sqrt(.add(lapply(a, .sq))))
.geo_zero <- function(v) is_zero_num(simplify(v))
.ldir <- function(l) .vsub(l$p2$coords, l$p1$coords)
.ldim <- function(l) length(l$p1$coords)

#' Geometry in three dimensions, rays, planes and conics
#'
#' `Point3D()`, `Line3D()`, `Segment3D()`, `Ray()` (`Ray2D()`, `Ray3D()`),
#' `Plane()`, `Ellipse()`, `Parabola()` and `Curve()` extend the plane
#' geometry of [Point()]. [distance()], [intersection()], [is_parallel()],
#' [area()] and [contains()] accept them. The accessors follow SymPy's
#' properties: `direction_ratio()`, `direction_cosine()`,
#' `arbitrary_point()`, `equation()`, `normal_vector()`, `projection()`,
#' `angle_between()`, `circumference()`, `eccentricity()`, `foci()`,
#' `major()`, `minor()` (an ellipse's minor radius; for a matrix the
#' [sym_det()] minor), `focus_distance()`, `encloses_point()`,
#' `focal_length()`, `vertex()`, `p_parameter()`, `axis_of_symmetry()`,
#' `translate()` and `geo_length()` (SymPy's `length`). `are_similar()`
#' tests similarity.
#'
#' @param ... Coordinates, or points.
#' @param p,q Points.
#' @param direction_ratio A direction instead of the second point.
#' @param p1,p2,p3 Points of a plane.
#' @param normal_vector The normal of a plane.
#' @param center,hradius,vradius Ellipse data.
#' @param focus,directrix Parabola data (an axis parallel directrix).
#' @param functions,limits Curve coordinates and `list(t, a, b)`.
#' @param e,e1,e2 Geometric entities.
#' @param t A parameter symbol.
#' @param x,y,z Symbols for equations, or offsets for `translate()`.
#' @param pt A point, or for `projection()` of a [CoordSys3D()] vector the
#'   vector to project (`Vector.projection`).
#' @param scalar Return the scalar projection factor instead.
#' @return An entity, an expression, a list or a logical.
#' @examples
#' P <- Point3D(1, 2, 3); Q <- Point3D(4, 6, 8)
#' L <- Line3D(P, Q)
#' direction_ratio(L)
#' distance(L, Point3D(0, 0, 0))
#' pl <- Plane(Point3D(0, 0, 0), normal_vector = c(1, 1, 1))
#' intersection(pl, Line3D(Point3D(0, 0, 5), Point3D(1, 1, 5)))
#' E <- Ellipse(Point(0, 0), 3, 2)
#' foci(E)
#' @name geometry-3d
NULL

#' @rdname geometry-3d
#' @export
Point3D <- function(...) {
  p <- Point(...)
  if (length(p$coords) != 3L) stop("Point3D needs three coordinates", call. = FALSE)
  p
}

#' @rdname geometry-3d
#' @export
Line3D <- function(p, q = NULL, direction_ratio = NULL) {
  if (is.null(q)) q <- do.call(Point, Map(function(a, b) .add(list(a, as_sym(b))), p$coords, as.list(direction_ratio)))
  .mk_line(p, q, "GeoInfLine")
}
#' @rdname geometry-3d
#' @export
Segment3D <- function(p, q) .mk_line(p, q, "GeoSegment")
#' @rdname geometry-3d
#' @export
Ray <- function(p, q) .mk_line(p, q, "GeoRay")
#' @rdname geometry-3d
#' @export
Ray2D <- Ray
#' @rdname geometry-3d
#' @export
Ray3D <- Ray

#' @rdname geometry-3d
#' @export
direction_ratio <- function(e) lapply(.ldir(e), simplify)
#' @rdname geometry-3d
#' @export
direction_cosine <- function(e) {
  d <- .ldir(e); n <- .vnorm(d)
  lapply(d, function(v) radsimp(simplify(.div(v, n))))
}

#' @rdname geometry-3d
#' @export
arbitrary_point <- function(e, t = Symbol("t")) {
  t <- as_sym(t)
  if (inherits(e, "GeoCurve")) {
    v <- e$limits[[1L]]
    return(do.call(Point, lapply(e$functions, function(f) subs(f, v, t))))
  }
  d <- .ldir(e)
  do.call(Point, Map(function(a, b) .add(list(a, .mul(list(t, b)))), e$p1$coords, d))
}

# Parameter of pt on line l (p1 + t*d), or NULL when pt is off the line.
.line_param <- function(l, pt) {
  d <- .ldir(l); w <- .vsub(pt$coords, l$p1$coords)
  k <- which(vapply(d, function(v) !.geo_zero(v), TRUE))[1L]
  t <- simplify(.div(w[[k]], d[[k]]))
  if (!all(Map(function(wi, di) .geo_zero(.sub(wi, .mul(list(t, di)))), w, d) == TRUE)) return(NULL)
  t
}

.geo_contains <- function(e, pt) {
  if (inherits(e, "GeoLine")) {
    t <- .line_param(e, pt)
    if (is.null(t)) return(FALSE)
    if (inherits(e, "GeoSegment")) return(isTRUE(.is_true(Ge(t, 0))) && isTRUE(.is_true(Le(t, 1))))
    if (inherits(e, "GeoRay")) return(isTRUE(.is_true(Ge(t, 0))))
    return(TRUE)
  }
  if (inherits(e, "GeoPlane")) return(.geo_zero(.vdot(e$normal, .vsub(pt$coords, e$p$coords))))
  if (inherits(e, "GeoCircle")) return(.geo_zero(circle_equation(e, pt$coords[[1L]], pt$coords[[2L]])))
  if (inherits(e, "GeoEllipse")) return(.geo_zero(equation(e, pt$coords[[1L]], pt$coords[[2L]])))
  stop("contains() is not implemented for this entity", call. = FALSE)
}

.dist_point_line3d <- function(pt, l) {
  d <- .ldir(l)
  simplify(.div(.vnorm(.vcross(.vsub(pt$coords, l$p1$coords), d)), .vnorm(d)))
}

.isect_lines3d <- function(a, b) {
  d1 <- .ldir(a); d2 <- .ldir(b)
  if (all(vapply(.vcross(d1, d2), .geo_zero, TRUE)))
    return(if (!is.null(.line_param(a, b$p1))) .isect_collinear(a, b) else list())
  t <- Dummy("t"); s <- Dummy("s")
  eqs <- Map(function(p, u, q, v) .sub(.add(list(p, .mul(list(t, u)))), .add(list(q, .mul(list(s, v))))),
             a$p1$coords, d1, b$p1$coords, d2)
  sol <- tryCatch(solve(eqs, list(t, s)), error = function(e) list())
  if (!length(sol)) return(list())
  tv <- sol[[1L]][[t$name]]
  pt <- do.call(Point, Map(function(p, u) simplify(.add(list(p, .mul(list(tv, u))))), a$p1$coords, d1))
  if (.on_segment(a, pt) && .on_segment(b, pt)) list(pt) else list()
}

#' @rdname geometry-3d
#' @export
Plane <- function(p1, p2 = NULL, p3 = NULL, normal_vector = NULL) {
  n <- if (!is.null(normal_vector)) lapply(as.list(normal_vector), as_sym)
       else lapply(.vcross(.vsub(p2$coords, p1$coords), .vsub(p3$coords, p1$coords)), simplify)
  if (all(vapply(n, is_zero_num, TRUE))) stop("the points of a plane must not be collinear", call. = FALSE)
  structure(list(p = p1, normal = n), class = c("GeoPlane", "GeoEntity"))
}
#' @export
format.GeoPlane <- function(x, ...) paste0("Plane(", format(x$p), ", (", paste(vapply(x$normal, .str, ""), collapse = ", "), "))")

#' @rdname geometry-3d
#' @export
normal_vector <- function(e) e$normal

#' @rdname geometry-3d
#' @export
equation <- function(e, x = Symbol("x"), y = Symbol("y"), z = Symbol("z")) {
  x <- as_sym(x); y <- as_sym(y); z <- as_sym(z)
  if (inherits(e, "GeoPlane")) return(expand(.vdot(e$normal, .vsub(list(x, y, z), e$p$coords))))
  if (inherits(e, "GeoCircle")) return(circle_equation(e, x, y))
  if (inherits(e, "GeoLine")) return(line_equation(e, x, y))
  if (inherits(e, "GeoEllipse")) {
    c <- e$center$coords
    return(.sub(.add(list(.pow(.div(.sub(x, c[[1L]]), e$hradius), S_Two), .pow(.div(.sub(y, c[[2L]]), e$vradius), S_Two))), S_One))
  }
  if (inherits(e, "GeoParabola")) {
    v <- vertex(e)$coords; p <- p_parameter(e)
    if (.parabola_horizontal(e))
      return(expand(.add(list(.neg(.sq(.sub(x, v[[1L]]))), .mul(list(.int(4), p, .sub(y, v[[2L]])))))))
    return(expand(.add(list(.neg(.sq(.sub(y, v[[2L]]))), .mul(list(.int(4), p, .sub(x, v[[1L]])))))))
  }
  stop("equation() is not implemented for this entity", call. = FALSE)
}

#' @rdname geometry-3d
#' @export
projection <- function(e, pt, scalar = FALSE) {
  if (.is_vec(e)) return(.vec_projection(e, pt, scalar))
  n <- e$normal
  k <- simplify(.div(.vdot(n, .vsub(pt$coords, e$p$coords)), .vdot(n, n)))
  do.call(Point, Map(function(a, b) simplify(.sub(a, .mul(list(k, b)))), pt$coords, n))
}

#' @rdname geometry-3d
#' @export
angle_between <- function(e1, e2) {
  a <- if (inherits(e1, "GeoPlane")) e1$normal else .ldir(e1)
  b <- if (inherits(e2, "GeoPlane")) e2$normal else .ldir(e2)
  # a plane and a line meet at the complement of the normal's angle: asin
  f <- if (xor(inherits(e1, "GeoPlane"), inherits(e2, "GeoPlane"))) "asin" else "acos"
  fcall(f, list(radsimp(simplify(.div(.vdot(a, b), .mul(list(.vnorm(a), .vnorm(b))))))))
}

.isect_plane_line <- function(pl, l) {
  d <- .ldir(l); nd <- simplify(.vdot(pl$normal, d))
  if (.geo_zero(nd)) return(if (.geo_contains(pl, l$p1)) list(l) else list())
  t <- simplify(.div(.neg(.vdot(pl$normal, .vsub(l$p1$coords, pl$p$coords))), nd))
  pt <- do.call(Point, Map(function(p, u) simplify(.add(list(p, .mul(list(t, u))))), l$p1$coords, d))
  if (.on_segment(l, pt)) list(pt) else list()
}

#' @rdname geometry-3d
#' @export
Ellipse <- function(center, hradius, vradius) {
  hradius <- as_sym(hradius); vradius <- as_sym(vradius)
  if (identical(simplify(.sub(hradius, vradius))$key, S_Zero$key)) return(Circle(center, hradius))
  structure(list(center = center, hradius = hradius, vradius = vradius), class = c("GeoEllipse", "GeoEntity"))
}
#' @export
format.GeoEllipse <- function(x, ...) paste0("Ellipse(", format(x$center), ", ", .str(x$hradius), ", ", .str(x$vradius), ")")

#' @rdname geometry-3d
#' @export
major <- function(e) {
  if (inherits(e, "GeoCircle")) return(e$radius)
  if (.is_true(Lt(e$hradius, e$vradius))) e$vradius else e$hradius
}
.ellipse_minor <- function(e) {
  if (inherits(e, "GeoCircle")) return(e$radius)
  if (identical(major(e)$key, e$hradius$key)) e$vradius else e$hradius
}
#' @rdname geometry-3d
#' @export
focus_distance <- function(e) simplify(.sqrt(.sub(.sq(major(e)), .sq(.ellipse_minor(e)))))
#' @rdname geometry-3d
#' @export
eccentricity <- function(e) {
  if (inherits(e, "GeoParabola")) return(S_One)
  if (inherits(e, "GeoCircle")) return(S_Zero)
  simplify(.div(focus_distance(e), major(e)))
}
#' @rdname geometry-3d
#' @export
foci <- function(e) {
  if (inherits(e, "GeoCircle")) return(list(e$center, e$center))
  c <- e$center$coords; fd <- focus_distance(e)
  if (identical(major(e)$key, e$hradius$key))
    list(Point(.sub(c[[1L]], fd), c[[2L]]), Point(.add(list(c[[1L]], fd)), c[[2L]]))
  else list(Point(c[[1L]], .sub(c[[2L]], fd)), Point(c[[1L]], .add(list(c[[2L]], fd))))
}
#' @rdname geometry-3d
#' @export
circumference <- function(e) {
  if (inherits(e, "GeoCircle")) return(.mul(list(S_Two, S_Pi, e$radius)))
  ecc <- eccentricity(e)
  if (identical(ecc$key, S_One$key)) return(.mul(list(.int(4), major(e))))
  if (is_zero_num(ecc)) return(.mul(list(S_Two, S_Pi, e$hradius)))
  .mul(list(.int(4), major(e), fcall("elliptic_e", list(.sq(ecc)))))
}
#' @rdname geometry-3d
#' @export
encloses_point <- function(e, pt) {
  if (inherits(e, "GeoPolygon")) return(.polygon_encloses(e, pt))
  v <- simplify(equation(e, pt$coords[[1L]], pt$coords[[2L]]))
  .is_true(Lt(v, 0))
}

# Polygon.encloses_point: the orientation test for a convex polygon, else
# the crossings of a ray from the point
.polygon_encloses <- function(poly, p) {
  v <- poly$vertices
  n <- length(v)
  side <- function(i) Segment(v[[i]], v[[i %% n + 1L]])
  on_side <- function(sg) {
    a <- sg$p1$coords; b <- sg$p2$coords; q <- p$coords
    cr <- simplify(.sub(.mul(list(.sub(b[[1L]], a[[1L]]), .sub(q[[2L]], a[[2L]]))),
                        .mul(list(.sub(b[[2L]], a[[2L]]), .sub(q[[1L]], a[[1L]])))))
    is_zero_num(cr) && .on_segment(sg, p)
  }
  if (any(vapply(v, function(w) isTRUE(w == p), TRUE)) || any(vapply(seq_len(n), function(i) on_side(side(i)), TRUE)))
    return(FALSE)
  lit <- lapply(v, function(w) lapply(1:2, function(k) simplify(.sub(w$coords[[k]], p$coords[[k]]))))
  if (any(vapply(lit, function(w) length(free_symbols(w[[1L]])) + length(free_symbols(w[[2L]])) > 0L, TRUE))) return(NA)
  at <- function(i) lit[[i %% n + 1L]]
  if (.polygon_is_convex(v)) {
    orientation <- NULL
    for (i in 0:n) {
      a <- at(i); b <- at(i + 1L)
      test <- ask_fact(.sub(.mul(list(.neg(a[[2L]]), .sub(b[[1L]], a[[1L]]))), .mul(list(.neg(a[[1L]]), .sub(b[[2L]], a[[2L]])))), "negative")
      if (is.null(orientation)) orientation <- test
      else if (!identical(test, orientation)) return(FALSE)
    }
    return(TRUE)
  }
  lt <- function(a, b) .is_true(Lt(a, b))
  hit_odd <- FALSE
  p1 <- at(0L)
  for (i in 1:n) {
    p2 <- at(i)
    if (lt(Min(p1[[2L]], p2[[2L]]), S_Zero)) {
      if (!lt(Max(p1[[2L]], p2[[2L]]), S_Zero) && !lt(Max(p1[[1L]], p2[[1L]]), S_Zero) &&
          !identical(p1[[2L]]$key, p2[[2L]]$key)) {
        xi <- .add(list(.div(.mul(list(.neg(p1[[2L]]), .sub(p2[[1L]], p1[[1L]]))), .sub(p2[[2L]], p1[[2L]])), p1[[1L]]))
        if (identical(p1[[1L]]$key, p2[[1L]]$key) || !lt(xi, S_Zero)) hit_odd <- !hit_odd
      }
    }
    p1 <- p2
  }
  hit_odd
}

# Polygon.is_convex: one orientation at every vertex and no crossing sides
.polygon_is_convex <- function(v) {
  n <- length(v)
  at <- function(i) v[[i %% n + 1L]]
  cw <- function(a, b, c) {
    ba <- lapply(1:2, function(k) .sub(b$coords[[k]], a$coords[[k]]))
    ca <- lapply(1:2, function(k) .sub(c$coords[[k]], a$coords[[k]]))
    r <- ask_fact(simplify(.sub(.mul(list(ba[[1L]], ca[[2L]])), .mul(list(ca[[1L]], ba[[2L]])))), "nonpositive")
    if (is.na(r)) stop(.sympy_condition("ValueError", "Can't determine orientation"))
    r
  }
  c0 <- cw(at(-2L), at(-1L), at(0L))
  for (i in seq_len(n - 1L)) if (xor(c0, cw(at(i - 2L), at(i - 1L), at(i)))) return(FALSE)
  sides <- lapply(seq_len(n), function(i) Segment(at(i - 1L), at(i)))
  same <- function(a, b) isTRUE(a == b)
  for (i in seq_len(n)) {
    si <- sides[[i]]
    js <- seq_len(max(0L, i - 2L))
    if (i == n) js <- js[-1L]
    for (j in js) {
      sj <- sides[[j]]
      if (!any(vapply(list(si$p1, si$p2), function(q) same(sj$p1, q) || same(sj$p2, q), TRUE)) &&
          length(intersection(si, sj))) return(FALSE)
    }
  }
  TRUE
}

.isect_conic_line <- function(e, l) {
  t <- Dummy("t")
  p <- l$p1$coords; d <- .ldir(l)
  px <- .add(list(p[[1L]], .mul(list(t, d[[1L]])))); py <- .add(list(p[[2L]], .mul(list(t, d[[2L]]))))
  eq <- fraction(together(expand(equation(e, px, py))))[[1L]]
  ts <- Filter(function(v) !identical(ask_fact(v, "real"), FALSE), roots(expand(eq), t))
  pts <- lapply(ts, function(v) Point(simplify(subs(px, t, v)), simplify(subs(py, t, v))))
  .ordered_points(Filter(function(pt) .on_segment(l, pt), pts))
}

#' @rdname geometry-3d
#' @export
Parabola <- function(focus, directrix) structure(list(focus = focus, directrix = directrix), class = c("GeoParabola", "GeoEntity"))
#' @export
format.GeoParabola <- function(x, ...) paste0("Parabola(", format(x$focus), ", ", format(x$directrix), ")")

.parabola_horizontal <- function(e) is_zero_num(simplify(.sub(e$directrix$p1$coords[[2L]], e$directrix$p2$coords[[2L]])))
#' @rdname geometry-3d
#' @export
focal_length <- function(e) simplify(.div(distance(e$focus, e$directrix), S_Two))
#' @rdname geometry-3d
#' @export
p_parameter <- function(e) {
  f <- e$focus$coords
  if (.parabola_horizontal(e)) simplify(.div(.sub(f[[2L]], e$directrix$p1$coords[[2L]]), S_Two))
  else simplify(.div(.sub(f[[1L]], e$directrix$p1$coords[[1L]]), S_Two))
}
#' @rdname geometry-3d
#' @export
vertex <- function(e) {
  f <- e$focus$coords; p <- p_parameter(e)
  if (.parabola_horizontal(e)) Point(f[[1L]], simplify(.sub(f[[2L]], p))) else Point(simplify(.sub(f[[1L]], p)), f[[2L]])
}
#' @rdname geometry-3d
#' @export
axis_of_symmetry <- function(e) {
  f <- e$focus$coords
  if (.parabola_horizontal(e)) Line(e$focus, Point(f[[1L]], .add(list(f[[2L]], S_One))))
  else Line(e$focus, Point(.add(list(f[[1L]], S_One)), f[[2L]]))
}

#' @rdname geometry-3d
#' @export
Curve <- function(functions, limits) {
  structure(list(functions = lapply(functions, as_sym), limits = lapply(limits, as_sym)), class = c("GeoCurve", "GeoEntity"))
}
#' @export
format.GeoCurve <- function(x, ...) paste0("Curve((", paste(vapply(x$functions, .str, ""), collapse = ", "), "), (",
                                           paste(vapply(x$limits, .str, ""), collapse = ", "), "))")

#' @rdname geometry-3d
#' @export
translate <- function(e, x = 0, y = 0, z = 0) {
  dim <- if (inherits(e, "GeoCurve")) length(e$functions) else length(e$coords)
  off <- list(as_sym(x), as_sym(y), as_sym(z))
  if (dim < 3L && !is_zero_num(off[[3L]])) stop("a z offset needs a 3D entity", call. = FALSE)
  off <- off[seq_len(dim)]
  if (inherits(e, "GeoCurve")) return(Curve(Map(function(f, o) .add(list(f, o)), e$functions, off), e$limits))
  if (inherits(e, "GeoPoint")) return(do.call(Point, Map(function(a, o) .add(list(a, o)), e$coords, off)))
  stop("translate() is not implemented for this entity", call. = FALSE)
}

#' @rdname geometry-3d
#' @export
geo_length <- function(e) {
  if (inherits(e, "GeoSegment")) return(distance(e$p1, e$p2))
  if (inherits(e, "GeoCurve")) {
    v <- e$limits[[1L]]
    sp <- .sqrt(.add(lapply(e$functions, function(f) .sq(.diff1(f, v)))))
    return(simplify(sym_integrate(sp, list(v, e$limits[[2L]], e$limits[[3L]]))))
  }
  if (inherits(e, "GeoLine")) return(S_Infinity)
  stop("geo_length() is not implemented for this entity", call. = FALSE)
}

#' @rdname geometry-3d
#' @export
are_similar <- function(e1, e2) {
  if (inherits(e1, "GeoLine") && inherits(e2, "GeoLine"))
    return(is_parallel(e1, e2) && !is.null(.line_param(e1, e2$p1)))
  if (inherits(e1, "GeoPolygon") && inherits(e2, "GeoPolygon") && length(e1$vertices) == 3L && length(e2$vertices) == 3L) {
    sides <- function(T) {
      v <- T$vertices
      s <- lapply(1:3, function(i) distance(v[[i]], v[[i %% 3L + 1L]]))
      s[order(vapply(s, function(z) num_value_dbl(N(z)), 0))]
    }
    a <- sides(e1); b <- sides(e2)
    r <- simplify(.div(a[[1L]], b[[1L]]))
    return(all(vapply(2:3, function(i) .geo_zero(.sub(.div(a[[i]], b[[i]]), r)), TRUE)))
  }
  stop("are_similar() is not implemented for these entities", call. = FALSE)
}
