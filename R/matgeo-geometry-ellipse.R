# Ellipses and circles, after sympy/geometry/ellipse.py.

.mg_ellipse_raw <- function(center, h, v) .mg_gnew("Ellipse", "GeoEllipse", list(center = center, hradius = h, vradius = v))

.mg_circle_raw <- function(center, r) {
  .mg_gnew("Circle", c("GeoCircle", "GeoEllipse"), list(center = center, radius = r, hradius = r, vradius = Abs(r)))
}

#' Ellipses and circles
#'
#' `Ellipse()` builds an ellipse from its center and two of `hradius`,
#' `vradius` and `eccentricity`; equal radii give a `Circle`, a zero radius
#' a segment. `Circle()` takes a center and a radius, three points, or an
#' equation in `x` and `y`. The fields `$center`, `$hradius`, `$vradius`
#' and, for a circle, `$radius` hold the defining data.
#'
#' The functions are SymPy's methods and properties: `major()`, `minor()`,
#' `foci()`, `focus_distance()`, `eccentricity()`, `circumference()`,
#' `periapsis()`, `apoapsis()`, `semilatus_rectum()`, `encloses_point()`,
#' `is_tangent()`, `tangent_lines()`, `normal_lines()`, `evolute()`,
#' `auxiliary_circle()`, `director_circle()`, `second_moment_of_area()`,
#' `polar_second_moment_of_area()` and `section_modulus()`.
#'
#' @param center The center.
#' @param hradius,vradius The horizontal and vertical radii.
#' @param eccentricity The eccentricity.
#' @param ... For `Circle()`: a center and a radius, three points, or an
#'   equation; for `normal_lines()`, `prec`.
#' @param x,y Names or symbols of the variables of an equation.
#' @param evaluate Replace floats in the points by rationals.
#' @param e An ellipse or circle.
#' @param o,p An entity or a point.
#' @param point A reference point.
#' @return An entity, an expression, a logical value or a list.
#' @examples
#' e <- Ellipse(Point(0, 0), 3, 2)
#' foci(e)
#' eccentricity(e)
#' equation(Circle(Point(1, 2), 3))
#' Circle(Point(0, 0), Point(1, 0), Point(0, 1))
#' tangent_lines(Circle(Point(0, 0), 1), Point(2, 0))
#' @export
Ellipse <- function(center = NULL, hradius = NULL, vradius = NULL, eccentricity = NULL) {
  if (!is.null(hradius)) hradius <- as_sym(hradius)
  if (!is.null(vradius)) vradius <- as_sym(vradius)
  if (is.null(center)) center <- Point(0, 0)
  else {
    if (length(.mg_seq(center)) != 2L) .mg_err("ValueError", "The center of \"Ellipse\" must be a two dimensional point")
    center <- .mg_as_pt(center, 2L)
  }
  if (sum(!vapply(list(hradius, vradius, eccentricity), is.null, TRUE)) != 2L)
    .mg_errf("ValueError", "Exactly two arguments of \"hradius\", \"vradius\", and \"eccentricity\" must not be None.")
  if (!is.null(eccentricity)) {
    eccentricity <- as_sym(eccentricity)
    if (isTRUE(ask_fact(eccentricity, "negative")))
      .mg_err("GeometryError", "Eccentricity of ellipse/circle should lie between [0, 1)")
    if (is.null(hradius)) hradius <- .div(vradius, .sqrt(.sub(S_One, .pow(eccentricity, S_Two))))
    else if (is.null(vradius)) vradius <- .mul(list(hradius, .sqrt(.sub(S_One, .pow(eccentricity, S_Two)))))
  }
  if (identical(hradius$key, vradius$key)) return(Circle(center, hradius))
  if (is_zero_num(hradius) || is_zero_num(vradius))
    return(Segment(Point(.sub(center$coords[[1L]], hradius), .sub(center$coords[[2L]], vradius)),
                   Point(.add(list(center$coords[[1L]], hradius)), .add(list(center$coords[[2L]], vradius)))))
  if (identical(ask_fact(hradius, "real"), FALSE) || identical(ask_fact(vradius, "real"), FALSE))
    .mg_err("GeometryError", "Invalid value encountered when computing hradius / vradius.")
  .mg_ellipse_raw(center, hradius, vradius)
}

#' @rdname Ellipse
#' @export
Circle <- function(..., x = "x", y = "y", evaluate = TRUE) {
  a <- list(...)
  if (length(a) == 1L && .mg_is_eqn(a[[1L]])) {
    eq <- expand(a[[1L]])
    if (eq$type == "Equality") eq <- .sub(eq$args[[1L]], eq$args[[2L]])
    xs <- .mg_find(x, eq); ys <- .mg_find(y, eq)
    cf <- tryCatch(linear_coeffs(eq, .pow(xs, S_Two), .pow(ys, S_Two), xs, ys),
                   error = function(e) .mg_err("GeometryError", "The given equation is not that of a circle."))
    A <- cf[[1L]]; B <- cf[[2L]]; C <- cf[[3L]]; D <- cf[[4L]]; E <- cf[[5L]]
    if (is_zero_num(A) || is_zero_num(B) || !identical(A$key, B$key))
      .mg_err("GeometryError", "The given equation is not that of a circle.")
    cx <- .div(.div(.neg(C), A), S_Two)
    cy <- .div(.div(.neg(D), B), S_Two)
    r2 <- .sub(.add(list(.pow(cx, S_Two), .pow(cy, S_Two))), .div(E, A))
    return(Circle(list(cx, cy), .sqrt(r2), evaluate = evaluate))
  }
  if (length(a) == 3L) {
    pts <- lapply(a, function(p) Point(p, dim = 2L, evaluate = evaluate))
    t <- do.call(Triangle, pts)
    if (!inherits(t, "GeoTriangle")) return(t)
    return(.mg_circle_new(circumcenter(t), circumradius(t)))
  }
  if (length(a) == 2L) {
    c <- Point(a[[1L]], dim = 2L, evaluate = evaluate)
    r <- tryCatch(Point(a[[2L]], 0, evaluate = evaluate)$coords[[1L]],
                  ValueError = function(e) .mg_err("GeometryError", "Circle with imaginary radius is not permitted"))
    return(.mg_circle_new(c, r))
  }
  .mg_err("GeometryError", "Circle.__new__ received unknown arguments")
}

.mg_circle_new <- function(c, r) {
  if (is_zero_num(r)) return(c)
  .mg_circle_raw(c, r)
}

.mg_need_ellipse <- function(e) if (!inherits(e, "GeoEllipse")) .mg_err("TypeError", "expected an Ellipse")

#' @rdname Ellipse
#' @export
major <- function(e) {
  if (inherits(e, "GeoCircle")) return(e$radius)
  a <- e$hradius; b <- e$vradius
  o <- .mg_fz(Lt(.sub(b, a), S_Zero))
  if (isTRUE(o)) return(a)
  if (identical(o, FALSE)) return(b)
  e$hradius
}

.ellipse_minor <- function(e) {
  if (inherits(e, "GeoCircle")) return(e$radius)
  a <- e$hradius; b <- e$vradius
  o <- .mg_fz(Lt(.sub(a, b), S_Zero))
  if (isTRUE(o)) return(a)
  if (identical(o, FALSE)) return(b)
  e$vradius
}

#' @rdname Ellipse
#' @export
foci <- function(e) {
  c <- e$center
  if (identical(e$hradius$key, e$vradius$key)) return(list(c, c))
  fd <- .sqrt(.sub(.pow(major(e), S_Two), .pow(.ellipse_minor(e), S_Two)))
  if (identical(e$hradius$key, .ellipse_minor(e)$key)) return(list(.mg_pt_add(c, Point(0, .neg(fd))), .mg_pt_add(c, Point(0, fd))))
  if (identical(e$hradius$key, major(e)$key)) return(list(.mg_pt_add(c, Point(.neg(fd), 0)), .mg_pt_add(c, Point(fd, 0))))
  NULL
}

#' @rdname Ellipse
#' @export
focus_distance <- function(e) .mg_pt_distance(e$center, foci(e)[[1L]])

#' @rdname Ellipse
#' @export
eccentricity <- function(e) {
  if (inherits(e, "GeoParabola")) return(S_One)
  .div(focus_distance(e), major(e))
}

#' @rdname Ellipse
#' @export
circumference <- function(e) {
  if (inherits(e, "GeoCircle")) return(.mul(list(S_Two, S_Pi, e$radius)))
  ec <- eccentricity(e)
  if (identical(ec$key, S_One$key)) return(.mul(list(.int(4), major(e))))
  if (is_zero_num(ec)) return(.mul(list(S_Two, S_Pi, e$hradius)))
  .mul(list(.int(4), major(e), elliptic_e(.pow(ec, S_Two))))
}

#' @rdname Ellipse
#' @export
periapsis <- function(e) .mul(list(major(e), .sub(S_One, eccentricity(e))))
#' @rdname Ellipse
#' @export
apoapsis <- function(e) .mul(list(major(e), .add(list(S_One, eccentricity(e)))))
#' @rdname Ellipse
#' @export
semilatus_rectum <- function(e) .mul(list(major(e), .sub(S_One, .pow(eccentricity(e), S_Two))))

.mg_ellipse_equation <- function(e, x = "x", y = "y", slope = NULL) {
  x <- .mg_param(x, "x"); y <- .mg_param(y, "y")
  dx <- .sub(x, e$center$coords[[1L]]); dy <- .sub(y, e$center$coords[[2L]])
  if (inherits(e, "GeoCircle")) return(.sub(.add(list(.pow(dx, S_Two), .pow(dy, S_Two))), .pow(major(e), S_Two)))
  if (!is.null(slope)) {
    slope <- as_sym(slope)
    L <- .pow(.sub(dy, .mul(list(slope, dx))), S_Two)
    l <- .pow(.add(list(.mul(list(slope, dy)), dx)), S_Two)
    h <- .add(list(S_One, .pow(slope, S_Two)))
    b <- .mul(list(h, .pow(major(e), S_Two)))
    a <- .mul(list(h, .pow(.ellipse_minor(e), S_Two)))
    return(.sub(.add(list(.div(l, b), .div(L, a))), S_One))
  }
  .sub(.add(list(.pow(.div(dx, e$hradius), S_Two), .pow(.div(dy, e$vradius), S_Two))), S_One)
}

.mg_ellipse_has <- function(e, o) {
  if (inherits(o, "GeoPoint")) {
    xx <- Dummy("ellipse_dummy", real = TRUE); yy <- Dummy("ellipse_dummy", real = TRUE)
    res <- subs(.mg_ellipse_equation(e, xx, yy), list(list(xx, o$coords[[1L]]), list(yy, o$coords[[2L]])))
    return(is_zero_num(trigsimp(simplify(res))))
  }
  if (inherits(o, "GeoEllipse")) return(.mg_geo_eq(e, o))
  FALSE
}

#' @rdname Ellipse
#' @export
encloses_point <- function(e, p) {
  if (inherits(e, "GeoRegularPolygon")) return(.mg_regpoly_encloses(e, p))
  if (inherits(e, "GeoPolygon")) return(.mg_poly_encloses(e, p))
  p <- .mg_as_pt(p, 2L)
  if (.mg_ellipse_has(e, p)) return(FALSE)
  f <- foci(e)
  if (length(f) == 2L) {
    h1 <- .mg_pt_distance(f[[1L]], p); h2 <- .mg_pt_distance(f[[2L]], p)
    test <- .sub(.mul(list(S_Two, major(e))), .add(list(h1, h2)))
  } else test <- .sub(e$radius, .mg_pt_distance(e$center, p))
  ask_fact(test, "positive")
}

.mg_ellipse_xy <- function() list(Dummy("ellipse_dummy", real = TRUE), Dummy("ellipse_dummy", real = TRUE))

# solve(..., set=True)[1] as a list of points
.mg_solve_points <- function(eqs, xs) {
  sol <- solve(eqs, xs, set = TRUE)
  s <- if (is.list(sol) && length(sol) == 2L && !is_sym(sol)) sol[[2L]] else sol
  if (is_sym(s)) s <- s$args
  lapply(s, function(tp) Point(if (is_sym(tp)) tp$args else tp))
}

.mg_sort_ents <- function(lst) {
  if (length(lst) < 2L) return(lst)
  keys <- vapply(lst, .mg_ent_sort_key, "")
  lst[order(keys, method = "radix")]
}

# a sort key for points that agrees with SymPy's ordered() on numbers
.mg_ent_sort_key <- function(p) {
  if (!inherits(p, "GeoPoint")) {
    a <- .mg_geo_args(p)
    if (!length(a) || !all(vapply(a, inherits, TRUE, "GeoEntity"))) return(paste0("~", format(p)))
    return(paste0("~", p$gtype, "(", paste(vapply(a, .mg_ent_sort_key, ""), collapse = ","), ")"))
  }
  paste(vapply(p$coords, function(z) {
    v <- tryCatch(as.numeric(.neval(z, list())), error = function(e) NA_real_)
    if (length(v) != 1L || is.na(v) || !(is_exact_number(z) || z$type == "Float")) return(paste0("~", default_sort_key(z)))
    sprintf("%s%025.12f", if (v < 0) "0" else "1", if (v < 0) 1e12 + v else v)
  }, ""), collapse = "|")
}

.mg_ellipse_intersection <- function(e, o) {
  xy <- .mg_ellipse_xy(); xx <- xy[[1L]]; yy <- xy[[2L]]
  if (inherits(o, "GeoPoint")) return(if (.mg_ellipse_has(e, o)) list(o) else list())
  if (inherits(o, c("GeoSegment2D", "GeoRay2D"))) {
    pts <- .mg_solve_points(list(.mg_ellipse_equation(e, xx, yy), .mg_line_equation(Line(o$p1, o$p2), xx, yy)), list(xx, yy))
    return(.mg_sort_ents(Filter(function(p) .mg_line_contains(o, p), pts)))
  }
  if (inherits(o, "GeoPolygon")) return(.mg_isect1(o, e))
  if (inherits(o, c("GeoEllipse", "GeoLine2D"))) {
    if (.mg_geo_eq(o, e)) return(e)
    oe <- if (inherits(o, "GeoEllipse")) .mg_ellipse_equation(o, xx, yy) else .mg_line_equation(o, xx, yy)
    return(.mg_sort_ents(.mg_solve_points(list(.mg_ellipse_equation(e, xx, yy), oe), list(xx, yy))))
  }
  if (inherits(o, "GeoLinearEntity3D")) .mg_err("TypeError", "Entity must be two dimensional, not three dimensional")
  .mg_err("TypeError", "Intersection not handled for ", o$gtype)
}

#' @rdname Ellipse
#' @export
is_tangent <- function(e, o) {
  if (inherits(o, "GeoPoint2D")) return(FALSE)
  if (inherits(o, "GeoEllipse")) {
    it <- .mg_ellipse_intersection(e, o)
    if (inherits(it, "GeoEllipse")) return(TRUE)
    if (length(it)) return(all(vapply(it, function(i) equals(tangent_lines(e, i)[[1L]], tangent_lines(o, i)[[1L]]), TRUE)))
    return(FALSE)
  }
  if (inherits(o, "GeoLine2D")) {
    hit <- .mg_ellipse_intersection(e, o)
    if (!length(hit)) return(FALSE)
    if (length(hit) == 1L) return(TRUE)
    return(equals(hit[[1L]], hit[[2L]]))
  }
  if (inherits(o, c("GeoSegment2D", "GeoRay2D"))) {
    it <- .mg_ellipse_intersection(e, o)
    if (length(it) == 1L) return(.mg_line_contains(tangent_lines(e, it[[1L]])[[1L]], o))
    return(FALSE)
  }
  if (inherits(o, "GeoPolygon")) return(all(vapply(sides(o), function(s) is_tangent(e, s), TRUE)))
  if (inherits(o, c("GeoLinearEntity3D", "GeoPoint3D"))) .mg_err("TypeError", "Entity must be two dimensional, not three dimensional")
  .mg_err("TypeError", "Is_tangent not handled for ", .mg_pytype(o))
}

#' @rdname Ellipse
#' @export
tangent_lines <- function(e, p) {
  p <- .mg_as_pt(p, 2L)
  if (isTRUE(encloses_point(e, p))) return(list())
  if (.mg_ellipse_has(e, p)) {
    delta <- .mg_pt_sub(e$center, p)
    rise <- .mul(list(.pow(e$vradius, S_Two), delta$coords[[1L]]))
    run <- .neg(.mul(list(.pow(e$hradius, S_Two), delta$coords[[2L]])))
    p2 <- Point(simplify(.add(list(p$coords[[1L]], run))), simplify(.add(list(p$coords[[2L]], rise))))
    return(list(Line(p, p2)))
  }
  f <- foci(e)
  if (length(f) == 2L) {
    test <- .sub(.sub(.mul(list(S_Two, e$hradius)), .mg_pt_distance(f[[1L]], p)), .mg_pt_distance(f[[2L]], p))
  } else test <- .sub(e$radius, .mg_pt_distance(e$center, p))
  if (.is_number(test) && isTRUE(ask_fact(test, "positive"))) return(list())
  xy <- .mg_ellipse_xy(); xx <- xy[[1L]]; yy <- xy[[2L]]
  eq <- .mg_ellipse_equation(e, xx, yy)
  dydx <- idiff(eq, yy, xx)
  sl <- slope(Line(p, Point(xx, yy)))
  tp <- solve(list(.sub(sl, dydx), eq), list(xx, yy))
  tp <- lapply(tp, function(z) if (is_sym(z)) z$args else z)
  if (length(tp) == 1L) {
    if (identical(tp[[1L]][[1L]]$key, p$coords[[1L]]$key) || identical(tp[[1L]][[2L]]$key, p$coords[[2L]]$key))
      return(list(Line(p, .mg_pt_add(p, Point(1, 0))), Line(p, .mg_pt_add(p, Point(0, 1)))))
    return(list(Line(p, .mg_pt_add(p, Point(0, 1))), Line(p, Point(tp[[1L]]))))
  }
  list(Line(p, Point(tp[[1L]])), Line(p, Point(tp[[2L]])))
}

#' @rdname Ellipse
#' @export
normal_lines <- function(e, p, ...) {
  prec <- list(...)$prec
  p <- .mg_as_pt(p, 2L)
  rv <- list()
  if (identical(p$coords[[1L]]$key, e$center$coords[[1L]]$key)) rv <- c(rv, list(Line(e$center, slope = S_Infinity)))
  if (identical(p$coords[[2L]]$key, e$center$coords[[2L]]$key)) rv <- c(rv, list(Line(e$center, slope = 0)))
  if (length(rv)) return(rv)
  xy <- .mg_ellipse_xy(); xx <- xy[[1L]]; yy <- xy[[2L]]
  eq <- .mg_ellipse_equation(e, xx, yy)
  dydx <- idiff(eq, yy, xx)
  norm <- .div(S_NegativeOne, dydx)
  sl <- slope(Line(p, list(xx, yy)))
  yis <- solve(.sub(sl, norm), yy)[[1L]]
  xeq <- expand(numer(together(subs(eq, yy, yis))))
  if (length(free_symbols(xeq)) != 1L) .mg_err("NotImplementedError", "intersections for the general ellipse are not supported")
  xsol <- tryCatch(real_roots(Poly(xeq, xx)), error = function(err) solve(xeq, xx))
  pts <- lapply(xsol, function(i) Point(i, solve(subs(eq, xx, i), yy)[[1L]]))
  slopes <- lapply(pts, function(pt) subs(norm, list(list(xx, pt$coords[[1L]]), list(yy, pt$coords[[2L]]))))
  if (!is.null(prec)) {
    pts <- lapply(pts, function(pt) .mg_geo_evalf(pt, prec))
    slopes <- lapply(slopes, function(s) if (.is_number(s) && !(s$type == "Constant")) N(s, prec) else s)
  }
  Map(function(pt, s) Line(pt, slope = s), pts, slopes)
}

#' @rdname Ellipse
#' @export
evolute <- function(e, x = "x", y = "y") {
  if (inherits(e, "GeoCircle")) .mg_err("NotImplementedError", "Evolute of arbitrary Ellipse is not supported.")
  x <- .mg_param(x, "x"); y <- .mg_param(y, "y")
  t1 <- .pow(.mul(list(e$hradius, .sub(x, e$center$coords[[1L]]))), .rat(2L, 3L))
  t2 <- .pow(.mul(list(e$vradius, .sub(y, e$center$coords[[2L]]))), .rat(2L, 3L))
  .sub(.add(list(t1, t2)), .pow(.sub(.pow(e$hradius, S_Two), .pow(e$vradius, S_Two)), .rat(2L, 3L)))
}

#' @rdname Ellipse
#' @export
auxiliary_circle <- function(e) Circle(e$center, Max(e$hradius, e$vradius))

#' @rdname Ellipse
#' @export
director_circle <- function(e) Circle(e$center, .sqrt(.add(list(.pow(e$hradius, S_Two), .pow(e$vradius, S_Two)))))

.mg_ellipse_area <- function(e) simplify(.mul(list(S_Pi, e$hradius, e$vradius)))

.mg_ellipse_arbitrary <- function(e, parameter = "t") {
  t <- .mg_param(parameter)
  .mg_check_param(e, t)
  Point(.add(list(e$center$coords[[1L]], .mul(list(e$hradius, cos(t))))),
        .add(list(e$center$coords[[2L]], .mul(list(e$vradius, sin(t))))))
}

.mg_ellipse_rotate <- function(e, angle = 0, pt = NULL) {
  angle <- as_sym(angle)
  if (identical(e$hradius$key, e$vradius$key)) return(Circle(.mg_pt_rotate(e$center, angle, pt), e$hradius))
  if (isTRUE(ask_fact(.div(angle, S_Pi), "integer"))) return(Ellipse(.mg_pt_rotate(e$center, angle, pt), e$hradius, e$vradius))
  if (isTRUE(ask_fact(.div(.mul(list(S_Two, angle)), S_Pi), "integer")))
    return(Ellipse(.mg_pt_rotate(e$center, angle, pt), e$vradius, e$hradius))
  .mg_err("NotImplementedError", "Only rotations of pi/2 are currently supported for Ellipse.")
}

.mg_ellipse_scale <- function(e, x = 1, y = 1, pt = NULL) {
  if (!is.null(pt)) return(.mg_scale_about(e, x, y, pt))
  x <- as_sym(x); y <- as_sym(y)
  c <- .mg_pt_scale(e$center, x, y)
  if (inherits(e, "GeoCircle")) {
    ax <- Abs(x); ay <- Abs(y)
    if (identical(ax$key, ay$key)) return(Circle(c, .mul(list(ax, e$radius))))
    return(Ellipse(c, hradius = .mul(list(e$radius, ax)), vradius = .mul(list(e$radius, ay))))
  }
  Ellipse(c, hradius = .mul(list(e$hradius, x)), vradius = .mul(list(e$vradius, y)))
}

.mg_ellipse_reflect <- function(e, line) {
  if (inherits(e, "GeoCircle")) return(Circle(.mg_reflect(e$center, line), .neg(e$radius)))
  s <- slope(line)
  if (is_zero_num(s) || identical(s$key, S_Infinity$key))
    return(Ellipse(.mg_reflect(e$center, line), .neg(e$hradius), e$vradius))
  xx <- Symbol("_x", real = TRUE); yy <- Symbol("_y", real = TRUE)
  ex <- .mg_ellipse_equation(e, xx, yy)
  p <- .mg_reflect(Point(xx, yy), line)
  res <- subs(ex, list(list(xx, p$coords[[1L]]), list(yy, p$coords[[2L]])), simultaneous = TRUE)
  .mg_errf("NotImplementedError", "General Ellipse is not supported but the equation of the reflected Ellipse is given by the zeros of: f(_x, _y) = ", .str(res))
}

#' @rdname Ellipse
#' @export
second_moment_of_area <- function(e, point = NULL) {
  if (inherits(e, "GeoPolygon")) return(.mg_poly_second_moment(e, point))
  h <- e$hradius; v <- e$vradius
  Ixx <- .div(.mul(list(S_Pi, h, .pow(v, .int(3)))), .int(4))
  Iyy <- .div(.mul(list(S_Pi, .pow(h, .int(3)), v)), .int(4))
  Ixy <- S_Zero
  if (is.null(point)) return(list(Ixx, Iyy, Ixy))
  pc <- lapply(.mg_seq(point), as_sym)
  A <- .mg_ellipse_area(e)
  dx <- .sub(pc[[1L]], e$center$coords[[1L]]); dy <- .sub(pc[[2L]], e$center$coords[[2L]])
  list(.add(list(Ixx, .mul(list(A, .pow(dy, S_Two))))), .add(list(Iyy, .mul(list(A, .pow(dx, S_Two))))),
       .add(list(Ixy, .mul(list(A, dx, dy)))))
}

#' @rdname Ellipse
#' @export
polar_second_moment_of_area <- function(e) {
  sm <- second_moment_of_area(e)
  .add(list(sm[[1L]], sm[[2L]]))
}

#' @rdname Ellipse
#' @export
section_modulus <- function(e, point = NULL) {
  cc <- if (inherits(e, "GeoPolygon")) centroid(e) else e$center
  xc <- cc$coords[[1L]]; yc <- cc$coords[[2L]]
  if (is.null(point)) {
    b <- bounds(e)
    yv <- .mg_max(list(.sub(yc, b[[2L]]), .sub(b[[4L]], yc)))
    xv <- .mg_max(list(.sub(xc, b[[1L]]), .sub(b[[3L]], xc)))
  } else {
    point <- .mg_as_pt(point, 2L)
    yv <- .sub(point$coords[[2L]], yc); xv <- .sub(point$coords[[1L]], xc)
  }
  sm <- second_moment_of_area(e)
  list(.div(sm[[1L]], yv), .div(sm[[2L]], xv))
}
