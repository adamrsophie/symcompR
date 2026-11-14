# Lines, rays and segments, after sympy/geometry/line.py.

.MG_LINE_CLASSES <- list(
  Line2D = c("GeoLine2D", "GeoInfLine", "GeoLinearEntity2D", "GeoLine"),
  Line3D = c("GeoLine3D", "GeoInfLine", "GeoLinearEntity3D", "GeoLine"),
  Line = c("GeoInfLine", "GeoLine"),
  Segment2D = c("GeoSegment2D", "GeoSegment", "GeoLinearEntity2D", "GeoLine"),
  Segment3D = c("GeoSegment3D", "GeoSegment", "GeoLinearEntity3D", "GeoLine"),
  Segment = c("GeoSegment", "GeoLine"),
  Ray2D = c("GeoRay2D", "GeoRay", "GeoLinearEntity2D", "GeoLine"),
  Ray3D = c("GeoRay3D", "GeoRay", "GeoLinearEntity3D", "GeoLine"),
  Ray = c("GeoRay", "GeoLine"))

.mg_line_raw <- function(gtype, p1, p2) {
  f <- list(p1 = p1, p2 = p2)
  if (startsWith(gtype, "Ray")) f$source <- p1
  .mg_gnew(gtype, .MG_LINE_CLASSES[[gtype]], f)
}

# LinearEntity.__new__
.mg_linear_new <- function(kind, p1, p2) {
  pr <- .mg_normdim(list(p1, p2))
  p1 <- pr[[1L]]; p2 <- pr[[2L]]
  d <- length(p1$coords)
  gtype <- if (d == 2L) paste0(kind, "2D") else if (d == 3L) paste0(kind, "3D") else kind
  if (.mg_geo_eq(p1, p2)) .mg_err("ValueError", gtype, ".__new__ requires two unique Points.")
  .mg_line_raw(gtype, p1, p2)
}

# rebuild a linear entity of the same kind from two points
.mg_line_new <- function(gtype, p1, p2) {
  if (startsWith(gtype, "Line")) return(Line(p1, p2))
  if (startsWith(gtype, "Ray")) return(Ray(p1, p2))
  Segment(p1, p2)
}

.mg_kind <- function(l) if (inherits(l, "GeoInfLine")) "Line" else if (inherits(l, "GeoRay")) "Ray" else "Segment"

.mg_is_eqn <- function(x) is_sym(x) && x$type != "Tuple"

#' Lines, rays and segments
#'
#' `Line()`, `Ray()` and `Segment()` build linear entities through two
#' points, as SymPy does; the dimension of the points gives `Line2D`,
#' `Line3D` and so on. A 2D line may also be given by a point and `slope`,
#' or by a linear equation in `x` and `y`; a ray by a point and `angle`; a
#' 3D line or ray by a point and `direction_ratio`. A segment between equal
#' points is that point.
#'
#' The other functions are the SymPy methods: `direction()`, `slope()`,
#' `geo_coefficients()` (SymPy's `coefficients`), `perpendicular_line()`,
#' `parallel_line()`, `perpendicular_segment()`,
#' `perpendicular_bisector()`, `smallest_angle_between()`,
#' `are_concurrent()`, `bisectors()`, `xdirection()`, `ydirection()`,
#' `zdirection()` and `closing_angle()`. The fields `$p1`, `$p2` and, for a
#' ray, `$source` hold the defining points.
#'
#' @param p1 A point, a linear entity, or (for `Line()`) an equation.
#' @param p2,pt A second point.
#' @param slope The slope of a 2D line.
#' @param angle The angle of a 2D ray.
#' @param direction_ratio The direction of a 3D line or ray.
#' @param x,y Names or symbols of the variables of an equation.
#' @param l,l1,l2,r1,r2 Linear entities.
#' @param p A point.
#' @param ... Lines.
#' @return An entity, an expression, a logical value or a list.
#' @examples
#' l <- Line(Point(0, 0), Point(1, 1))
#' slope(l)
#' equation(l)
#' intersection(l, Line(Point(0, 2), Point(2, 0)))
#' Segment(Point(1, 1), Point(1, 1))
#' Line(Point(0, 0), slope = 2)
#' @export
Line <- function(p1, p2 = NULL, slope = NULL, direction_ratio = NULL, x = NULL, y = NULL) {
  if (.mg_is_eqn(p1) && is.null(p2)) {
    eq <- p1
    if (eq$type == "Equality") eq <- .sub(eq$args[[1L]], eq$args[[2L]])
    miss <- Symbol("?")
    kw <- !is.null(x) || !is.null(y)
    fm <- function(v, dflt) {
      if (is.null(v)) { if (kw) return(miss); v <- dflt }
      tryCatch(.mg_find(v, eq), ValueError = function(e) miss)
    }
    xs <- fm(x, "x"); ys <- fm(y, "y")
    cf <- linear_coeffs(eq, xs, ys)
    a <- cf[[1L]]; b <- cf[[2L]]; c <- cf[[3L]]
    if (.mg_truthy(b)) return(Line(list(S_Zero, .div(.neg(c), b)), slope = .div(.neg(a), b)))
    if (.mg_truthy(a)) return(Line(list(.div(.neg(c), a), S_Zero), slope = S_Infinity))
    .mg_err("ValueError", "not found in equation: ", paste(setdiff(c("x", "y"), c(xs$name, ys$name)), collapse = ", "))
  }
  if (inherits(p1, "GeoLine")) {
    if (!is.null(p2)) .mg_err("ValueError", "If p1 is a LinearEntity, p2 must be None.")
    d <- length(p1$p1$coords)
  } else {
    p1 <- .mg_as_pt(p1)
    d <- length(p1$coords)
    if (!is.null(p2)) p2 <- .mg_as_pt(p2)
  }
  if (d == 2L) return(Line2D(p1, p2, slope = slope))
  if (d == 3L) return(Line3D(p1, p2, direction_ratio = direction_ratio))
  .mg_linear_new("Line", p1, p2)
}

#' @rdname Line
#' @export
Line2D <- function(p1, pt = NULL, slope = NULL) {
  if (inherits(p1, "GeoLine")) {
    if (!is.null(pt)) .mg_err("ValueError", "When p1 is a LinearEntity, pt should be None")
    pr <- .mg_normdim(list(p1$p1, p1$p2), dim = 2L)
    p1 <- pr[[1L]]; pt <- pr[[2L]]
  } else p1 <- .mg_as_pt(p1, 2L)
  if (!is.null(pt) && is.null(slope)) {
    p2 <- tryCatch(.mg_as_pt(pt, 2L), error = function(e)
      .mg_errf("ValueError", "The 2nd argument was not a valid Point. If it was a slope, enter it with keyword \"slope\"."))
  } else if (!is.null(slope) && is.null(pt)) {
    slope <- as_sym(slope)
    if (identical(ask_fact(slope, "finite"), FALSE)) { dx <- S_Zero; dy <- S_One }
    else { dx <- S_One; dy <- slope }
    p2 <- Point(.add(list(p1$coords[[1L]], dx)), .add(list(p1$coords[[2L]], dy)), evaluate = FALSE)
  } else .mg_err("ValueError", "A 2nd Point or keyword \"slope\" must be used.")
  .mg_linear_new("Line", p1, p2)
}

#' @rdname Line
#' @export
Line3D <- function(p1, pt = NULL, direction_ratio = NULL) {
  if (inherits(p1, "GeoLinearEntity3D")) {
    if (!is.null(pt)) .mg_err("ValueError", "if p1 is a LinearEntity, pt must be None.")
    pt <- p1$p2; p1 <- p1$p1
  } else p1 <- .mg_as_pt(p1, 3L)
  dr <- if (is.null(direction_ratio)) list() else lapply(.mg_seq(direction_ratio), as_sym)
  if (!is.null(pt) && !length(dr)) pt <- .mg_as_pt(pt, 3L)
  else if (length(dr) == 3L && is.null(pt)) pt <- Point3D(Map(function(a, b) .add(list(a, b)), p1$coords, dr))
  else .mg_err("ValueError", "A 2nd Point or keyword \"direction_ratio\" must be used.")
  .mg_l3_new("Line", p1, pt)
}

.mg_l3_new <- function(kind, p1, p2) {
  p1 <- .mg_as_pt(p1, 3L); p2 <- .mg_as_pt(p2, 3L)
  if (.mg_geo_eq(p1, p2)) .mg_err("ValueError", kind, "3D.__new__ requires two unique Points.")
  .mg_line_raw(paste0(kind, "3D"), p1, p2)
}

#' @rdname Line
#' @export
Ray <- function(p1, p2 = NULL, angle = NULL, direction_ratio = NULL) {
  p1 <- .mg_as_pt(p1)
  if (!is.null(p2)) {
    pr <- .mg_normdim(list(p1, .mg_as_pt(p2)))
    p1 <- pr[[1L]]; p2 <- pr[[2L]]
  }
  d <- length(p1$coords)
  if (d == 2L) return(Ray2D(p1, p2, angle = angle))
  if (d == 3L) return(Ray3D(p1, p2, direction_ratio = direction_ratio))
  .mg_linear_new("Ray", p1, p2)
}

#' @rdname Line
#' @export
Ray2D <- function(p1, pt = NULL, angle = NULL) {
  p1 <- .mg_as_pt(p1, 2L)
  if (!is.null(pt) && is.null(angle)) {
    p2 <- tryCatch(.mg_as_pt(pt, 2L), error = function(e)
      .mg_errf("ValueError", "The 2nd argument was not a valid Point; if it was meant to be an angle it should be given with keyword \"angle\"."))
    if (.mg_geo_eq(p1, p2)) .mg_err("ValueError", "A Ray requires two distinct points.")
  } else if (!is.null(angle) && is.null(pt)) {
    angle <- as_sym(angle)
    cf <- .pi_coeff_sym(angle)
    p2 <- NULL
    if (!is.null(cf)) {
      if (cf$type %in% c("Integer", "Rational")) {
        q <- gmp::denominator(gmp::as.bigq(cf$value)); pn <- gmp::numerator(gmp::as.bigq(cf$value))
        if (q == 2 && pn == 1) p2 <- .mg_pt_add(p1, Point(0, 1))
        else if (q == 2 && pn == 3) p2 <- .mg_pt_add(p1, Point(0, -1))
        else if (q == 1 && pn == 0) p2 <- .mg_pt_add(p1, Point(1, 0))
        else if (q == 1 && pn == 1) p2 <- .mg_pt_add(p1, Point(-1, 0))
      }
      c <- if (is.null(p2)) .mul(list(cf, S_Pi)) else NULL
    } else c <- Mod(angle, .mul(list(S_Two, S_Pi)))
    if (is.null(p2)) {
      m <- .div(.mul(list(S_Two, c)), S_Pi)
      left <- And(Lt(S_One, m), Lt(m, .int(3)))
      xx <- Piecewise(list(.int(-1), left), list(Piecewise(list(S_Zero, Eq(Mod(m, S_One), S_Zero)), list(S_One, S_true)), S_true))
      yy <- Piecewise(list(.neg(tan(c)), left),
                      list(Piecewise(list(S_One, Eq(m, S_One)), list(.int(-1), Eq(m, .int(3))), list(tan(c), S_true)), S_true))
      p2 <- .mg_pt_add(p1, Point(xx, yy))
    }
  } else .mg_err("ValueError", "A 2nd point or keyword \"angle\" must be used.")
  .mg_linear_new("Ray", p1, p2)
}

#' @rdname Line
#' @export
Ray3D <- function(p1, pt = NULL, direction_ratio = NULL) {
  if (inherits(p1, "GeoLinearEntity3D")) {
    if (!is.null(pt)) .mg_err("ValueError", "If p1 is a LinearEntity, pt must be None")
    pt <- p1$p2; p1 <- p1$p1
  } else p1 <- .mg_as_pt(p1, 3L)
  dr <- if (is.null(direction_ratio)) list() else lapply(.mg_seq(direction_ratio), as_sym)
  if (!is.null(pt) && !length(dr)) pt <- .mg_as_pt(pt, 3L)
  else if (length(dr) == 3L && is.null(pt)) pt <- Point3D(Map(function(a, b) .add(list(a, b)), p1$coords, dr))
  else .mg_errf("ValueError", "A 2nd Point or keyword \"direction_ratio\" must be used.")
  .mg_l3_new("Ray", p1, pt)
}

#' @rdname Line
#' @export
Segment <- function(p1, p2) {
  pr <- .mg_normdim(list(.mg_as_pt(p1), .mg_as_pt(p2)))
  p1 <- pr[[1L]]; p2 <- pr[[2L]]
  d <- length(p1$coords)
  if (d %in% 2:3) {
    if (.mg_geo_eq(p1, p2)) return(p1)
    return(.mg_line_raw(paste0("Segment", d, "D"), p1, p2))
  }
  .mg_linear_new("Segment", p1, p2)
}

#' @rdname Line
#' @export
Segment2D <- function(p1, p2) Segment(.mg_as_pt(p1, 2L), .mg_as_pt(p2, 2L))

#' @rdname Line
#' @export
Segment3D <- function(p1, p2) Segment(.mg_as_pt(p1, 3L), .mg_as_pt(p2, 3L))

# ---------------------------------------------------------------------------
# LinearEntity methods

.mg_ldir <- function(l) .mg_pt_sub(l$p2, l$p1)

#' @rdname Line
#' @export
direction <- function(l) .mg_ldir(l)

.mg_need_linear <- function(l1, l2) {
  if (!inherits(l1, "GeoLine") || !inherits(l2, "GeoLine"))
    .mg_err("TypeError", "Must pass only LinearEntity objects")
}

.mg_abs_pt <- function(v) .mg_pt_distance(origin(v), v)

.mg_line_angle <- function(l1, l2) {
  .mg_need_linear(l1, l2)
  v1 <- .mg_ldir(l1); v2 <- .mg_ldir(l2)
  acos(.div(.mg_pdot(v1, v2), .mul(list(.mg_abs_pt(v1), .mg_abs_pt(v2)))))
}

#' @rdname Line
#' @export
smallest_angle_between <- function(l1, l2) {
  .mg_need_linear(l1, l2)
  v1 <- .mg_ldir(l1); v2 <- .mg_ldir(l2)
  acos(.div(Abs(.mg_pdot(v1, v2)), .mul(list(.mg_abs_pt(v1), .mg_abs_pt(v2)))))
}

# _symbol(parameter, real=True)
.mg_param <- function(parameter, dflt = "t") {
  if (is.null(parameter)) parameter <- dflt
  if (is_sym(parameter)) return(parameter)
  Symbol(as.character(parameter), real = TRUE)
}

.mg_check_param <- function(e, t) {
  if (t$name %in% vapply(.mg_geo_free(e), function(s) if (is.null(s$name)) "" else s$name, ""))
    .mg_errf("ValueError", "Symbol ", t$name, " already appears in object and cannot be used as a parameter.")
}

.mg_line_arbitrary <- function(l, parameter = "t") {
  t <- .mg_param(parameter)
  .mg_check_param(l, t)
  .mg_pt_add(l$p1, .mg_pt_mul(.mg_pt_sub(l$p2, l$p1), t))
}

#' @rdname Line
#' @export
are_concurrent <- function(...) {
  ents <- list(...)
  if (length(ents) && all(vapply(ents, inherits, TRUE, "GeoPlane"))) return(.mg_planes_concurrent(ents))
  if (length(ents) < 2L) return(FALSE)
  res <- do.call(intersection, ents)
  length(res) == 1L && inherits(res[[1L]], "GeoPoint")
}

.mg_line_contains <- function(l, o) {
  if (!inherits(o, "GeoEntity")) o <- Point(o, dim = length(l$p1$coords))
  if (inherits(l, "GeoInfLine")) {
    if (inherits(o, "GeoPoint")) return(is_collinear(o, l$p1, l$p2))
    if (inherits(o, "GeoLine")) return(is_collinear(l$p1, l$p2, o$p1, o$p2))
    return(FALSE)
  }
  if (inherits(l, "GeoRay")) {
    if (inherits(o, "GeoPoint")) {
      if (is_collinear(l$p1, l$p2, o))
        return(.mg_rbool(Ge(.mg_pdot(.mg_pt_sub(l$p2, l$p1), .mg_pt_sub(o, l$p1)), S_Zero)))
      return(FALSE)
    }
    if (inherits(o, "GeoRay")) {
      if (is_collinear(l$p1, l$p2, o$p1, o$p2))
        return(.mg_rbool(Gt(.mg_pdot(.mg_pt_sub(l$p2, l$p1), .mg_pt_sub(o$p2, o$p1)), S_Zero)))
      return(FALSE)
    }
    if (inherits(o, "GeoSegment")) return(.mg_line_contains(l, o$p1) && .mg_line_contains(l, o$p2))
    return(FALSE)
  }
  # Segment
  if (inherits(o, "GeoPoint")) {
    if (is_collinear(o, l$p1, l$p2)) {
      if (inherits(l, "GeoSegment2D")) {
        vert <- .expr_equals(.div(S_One, slope(l)), S_Zero)
        if (identical(vert, FALSE)) {
          isin <- Le(.mul(list(.sub(l$p1$coords[[1L]], o$coords[[1L]]), .sub(l$p2$coords[[1L]], o$coords[[1L]]))), S_Zero)
          if (is_sym(isin) && isin$type == "BooleanAtom") return(isin$name == "true")
        }
        if (identical(vert, TRUE)) {
          isin <- Le(.mul(list(.sub(l$p1$coords[[2L]], o$coords[[2L]]), .sub(l$p2$coords[[2L]], o$coords[[2L]]))), S_Zero)
          if (is_sym(isin) && isin$type == "BooleanAtom") return(isin$name == "true")
        }
      }
      d1 <- .mg_pt_sub(o, l$p1); d2 <- .mg_pt_sub(o, l$p2); d <- .mg_pt_sub(l$p2, l$p1)
      r <- simplify(Eq(.sub(.add(list(.mg_abs_pt(d1), .mg_abs_pt(d2))), .mg_abs_pt(d)), S_Zero))
      if (is_sym(r) && r$type == "BooleanAtom") return(r$name == "true")
      .mg_err("Undecidable", "Cannot determine if ", format(o), " is in ", format(l))
    }
  }
  if (inherits(o, "GeoSegment")) return(.mg_line_contains(l, o$p1) && .mg_line_contains(l, o$p2))
  FALSE
}

.mg_span_test <- function(l, o) {
  if (.mg_geo_eq(l$p1, o)) return(0L)
  if (.mg_rbool(Gt(.mg_pdot(.mg_ldir(l), .mg_pt_sub(o, l$p1)), S_Zero))) 1L else -1L
}

.mg_has_float <- function(e) {
  for (p in .mg_geo_points(e)) for (z in p$coords) if (length(atoms(z, "Float"))) return(TRUE)
  FALSE
}

.mg_line_intersection <- function(self, other) {
  par_rays <- function(r1, r2) {
    if (.mg_rbool(Gt(.mg_pdot(.mg_ldir(r1), .mg_ldir(r2)), S_Zero)))
      return(if (.mg_span_test(r1, r2$p1) >= 0) list(r2) else list(r1))
    st <- .mg_span_test(r1, r2$p1)
    if (st < 0) return(list())
    if (st == 0) return(list(r2$p1))
    list(Segment(r1$p1, r2$p1))
  }
  par_ray_seg <- function(ray, seg) {
    s1 <- .mg_span_test(ray, seg$p1); s2 <- .mg_span_test(ray, seg$p2)
    if (s1 < 0 && s2 < 0) return(list())
    if (s1 >= 0 && s2 >= 0) return(list(seg))
    if (s1 >= 0) return(list(Segment(ray$p1, seg$p1)))
    list(Segment(ray$p1, seg$p2))
  }
  par_segs <- function(s1, s2) {
    if (.mg_line_contains(s1, s2)) return(list(s2))
    if (.mg_line_contains(s2, s1)) return(list(s1))
    if (.mg_rbool(Lt(.mg_pdot(.mg_ldir(s1), .mg_ldir(s2)), S_Zero))) s2 <- Segment(s2$p2, s2$p1)
    if (.mg_span_test(s1, s2$p1) < 0) { tmp <- s1; s1 <- s2; s2 <- tmp }
    if (.mg_span_test(s2, s1$p2) < 0) return(list())
    list(Segment(s2$p1, s1$p2))
  }
  if (!inherits(other, "GeoEntity")) other <- Point(other, dim = length(self$p1$coords))
  if (inherits(other, "GeoPoint")) return(if (.mg_line_contains(self, other)) list(other) else list())
  if (!inherits(other, "GeoLine")) return(.mg_isect1(other, self))
  pts <- .mg_normdim(list(self$p1, self$p2, other$p1, other$p2))
  rank <- do.call(affine_rank, pts)
  if (rank == 1L) {
    if (inherits(self, "GeoInfLine")) return(list(other))
    if (inherits(other, "GeoInfLine")) return(list(self))
    if (inherits(self, "GeoRay") && inherits(other, "GeoRay")) return(par_rays(self, other))
    if (inherits(self, "GeoRay") && inherits(other, "GeoSegment")) return(par_ray_seg(self, other))
    if (inherits(self, "GeoSegment") && inherits(other, "GeoRay")) return(par_ray_seg(other, self))
    return(par_segs(self, other))
  }
  if (rank == 2L) {
    l1 <- Line(pts[[1L]], pts[[2L]]); l2 <- Line(pts[[3L]], pts[[4L]])
    d1 <- .mg_ldir(l1); d2 <- .mg_ldir(l2)
    if (is_scalar_multiple(d1, d2)) return(list())
    v <- .mg_pt_sub(l2$p1, l1$p1)
    n <- length(d1$coords)
    rows <- lapply(seq_len(n), function(i) list(d1$coords[[i]], .neg(d2$coords[[i]]), v$coords[[i]]))
    rr <- rref(Matrix(rows), simplify = TRUE)
    if (length(rr$pivots) != 2L)
      .mg_err("GeometryError", "Failed when solving Mx=b when M=", sstr(Matrix(lapply(rows, `[`, 1:2))), " and b=", sstr(Matrix(lapply(rows, `[[`, 3L))))
    coeff <- .mget(rr$matrix, 1L, 3L)
    pt <- .mg_pt_add(.mg_pt_mul(d1, coeff), self$p1)
    if (inherits(self, "GeoInfLine") && inherits(other, "GeoInfLine")) return(list(pt))
    if ((inherits(self, "GeoInfLine") || .mg_line_contains(self, pt)) && .mg_line_contains(other, pt)) return(list(pt))
    if (!.mg_has_float(self) && !.mg_has_float(other)) return(list())
    tt <- Dummy("line_dummy"); uu <- Dummy("line_dummy")
    eqs <- Map(.sub, .mg_line_arbitrary(self, tt)$coords, .mg_line_arbitrary(other, uu)$coords)
    sol <- solve(eqs, list(tt, uu), dict = TRUE)[[1L]]
    getv <- function(s) { for (nm in names(sol)) if (identical(nm, s$key) || identical(nm, .str(s))) return(sol[[nm]]); NULL }
    ok <- function(p, l) {
      if (inherits(l, "GeoInfLine")) return(TRUE)
      if (inherits(l, "GeoRay")) return(isTRUE(ask_fact(p, "nonnegative")))
      isTRUE(ask_fact(p, "nonnegative")) && isTRUE(ask_fact(.sub(S_One, p), "nonnegative"))
    }
    if (ok(getv(tt), self) && ok(getv(uu), other)) return(list(pt))
    return(list())
  }
  list()
}

#' @rdname Line
#' @export
perpendicular_line <- function(l, p) {
  if (inherits(l, "GeoPlane")) return(Line3D(p, direction_ratio = l$normal_vector))
  p <- .mg_as_pt(p, length(l$p1$coords))
  if (inherits(l, "GeoLinearEntity2D")) return(Line(p, .mg_pt_add(p, orthogonal_direction(.mg_ldir(l)))))
  if (.mg_line_contains(l, p)) p <- .mg_pt_add(p, orthogonal_direction(.mg_ldir(l)))
  Line(p, .mg_line_projection(l, p))
}

#' @rdname Line
#' @export
parallel_line <- function(l, p) {
  p <- .mg_as_pt(p, length(l$p1$coords))
  Line(p, .mg_pt_add(p, .mg_ldir(l)))
}

#' @rdname Line
#' @export
perpendicular_segment <- function(l, p) {
  p <- .mg_as_pt(p, length(l$p1$coords))
  if (.mg_line_contains(l, p)) return(p)
  pl <- perpendicular_line(l, p)
  p2 <- .mg_line_intersection(Line(l$p1, l$p2), pl)[[1L]]
  Segment(p, p2)
}

#' @rdname Line
#' @export
perpendicular_bisector <- function(l, p = NULL) {
  m <- midpoint(l$p1, l$p2)
  pl <- perpendicular_line(l, m)
  if (!is.null(p)) {
    p2 <- .mg_as_pt(p, length(l$p1$coords))
    if (.mg_line_contains(pl, p2)) return(Segment(p2, m))
  }
  pl
}

.mg_line_projection <- function(l, other) {
  if (!inherits(other, "GeoEntity")) other <- Point(other, dim = length(l$p1$coords))
  pp <- function(p) .mg_pt_add(project(.mg_pt_sub(p, l$p1), .mg_ldir(l)), l$p1)
  if (inherits(other, "GeoPoint")) return(pp(other))
  if (inherits(other, "GeoLine")) {
    a <- pp(other$p1); b <- pp(other$p2)
    if (.mg_geo_eq(a, b)) return(a)
    pr <- .mg_line_new(other$gtype, a, b)
    res <- .mg_line_intersection(l, pr)
    if (!length(res)) return(S_EmptySet)
    if (length(res) == 1L && inherits(res[[1L]], "GeoPoint")) return(res[[1L]])
    pr <- res[[1L]]
    if (.mg_rbool(Lt(.mg_pdot(.mg_ldir(l), .mg_ldir(pr)), S_Zero))) pr <- .mg_line_new(pr$gtype, pr$p2, pr$p1)
    return(pr)
  }
  .mg_err("GeometryError", "Do not know how to project ", format(other), " onto ", format(l))
}

.mg_line_distance <- function(l, o) {
  if (inherits(l, "GeoLine3D") && !inherits(o, "GeoPoint")) {
    if (is.list(o) && !inherits(o, "GeoEntity")) o <- Point3D(o)
    if (inherits(o, "GeoLine3D")) {
      if (.mg_geo_eq(l, o)) return(S_Zero)
      if (.mg_is_parallel(l, o)) return(.mg_line_distance(l, o$p1))
      n <- cross(Matrix(direction_ratio(l)), Matrix(direction_ratio(o)))
      return(.mg_plane_distance(Plane(l$p1, normal_vector = n$data), o$p1))
    }
    if (inherits(o, "GeoPlane")) return(.mg_plane_distance(o, l))
    if (!inherits(o, "GeoPoint")) .mg_err("NotImplementedError", format(o), " has type ", o$gtype, ", which is unsupported")
  }
  if (!inherits(o, "GeoEntity")) o <- Point(o, dim = length(l$p1$coords))
  if (inherits(l, "GeoInfLine")) {
    if (.mg_line_contains(l, o)) return(S_Zero)
    return(geo_length(perpendicular_segment(l, o)))
  }
  if (inherits(l, "GeoRay")) {
    if (.mg_line_contains(l, o)) return(S_Zero)
    pr <- .mg_line_projection(Line(l$p1, l$p2), o)
    if (.mg_line_contains(l, pr)) return(.mg_abs_pt(.mg_pt_sub(o, pr)))
    return(.mg_abs_pt(.mg_pt_sub(o, l$source)))
  }
  if (inherits(o, "GeoPoint")) {
    vp1 <- .mg_pt_sub(o, l$p1); vp2 <- .mg_pt_sub(o, l$p2)
    s1 <- .mg_rbool(Ge(.mg_pdot(.mg_ldir(l), vp1), S_Zero))
    s2 <- .mg_rbool(Le(.mg_pdot(.mg_ldir(l), vp2), S_Zero))
    if (s1 && s2) return(.mg_line_distance(Line(l$p1, l$p2), o))
    if (s1 && !s2) return(.mg_abs_pt(vp2))
    if (!s1 && s2) return(.mg_abs_pt(vp1))
  }
  .mg_err("NotImplementedError", "")
}

.mg_is_parallel <- function(l1, l2) {
  .mg_need_linear(l1, l2)
  is_scalar_multiple(.mg_ldir(l1), .mg_ldir(l2))
}

.mg_is_perpendicular <- function(l1, l2) {
  .mg_need_linear(l1, l2)
  r <- .expr_equals(S_Zero, .mg_pdot(.mg_ldir(l1), .mg_ldir(l2)))
  if (is.null(r)) NA else r
}

#' @rdname Line
#' @export
slope <- function(l) {
  d <- .mg_pt_sub(l$p1, l$p2)$coords
  if (identical(d[[1L]]$key, S_Zero$key)) return(S_Infinity)
  simplify(.div(d[[2L]], d[[1L]]))
}

#' @rdname Line
#' @export
geo_coefficients <- function(l) {
  p1 <- l$p1; p2 <- l$p2
  if (identical(p1$coords[[1L]]$key, p2$coords[[1L]]$key)) return(list(S_One, S_Zero, .neg(p1$coords[[1L]])))
  if (identical(p1$coords[[2L]]$key, p2$coords[[2L]]$key)) return(list(S_Zero, S_One, .neg(p1$coords[[2L]])))
  x1 <- p1$coords[[1L]]; y1 <- p1$coords[[2L]]; x2 <- p2$coords[[1L]]; y2 <- p2$coords[[2L]]
  lapply(list(.sub(y1, y2), .sub(x2, x1), .sub(.mul(list(x1, y2)), .mul(list(y1, x2)))), simplify)
}

.mg_line_equation <- function(l, x = "x", y = "y", z = "z") {
  if (inherits(l, "GeoLinearEntity3D")) {
    vs <- lapply(list(x, y, z, "k"), .mg_param)
    k <- vs[[4L]]
    d <- direction_ratio(l$p1, l$p2)
    c1 <- l$p1$coords
    eqs <- lapply(1:3, function(i) .add(list(.neg(.mul(list(d[[i]], k))), vs[[i]], .neg(c1[[i]]))))
    kk <- NULL
    for (i in seq_along(eqs)) if (.has1(eqs[[i]], k)) {
      kk <- solve(eqs[[i]], k)[[1L]]
      eqs <- eqs[-i]
      break
    }
    return(new_sym("Tuple", lapply(eqs, function(e) .as_numer_denom(cancel(subs(e, k, kk)))[[1L]])))
  }
  x <- .mg_param(x, "x"); y <- .mg_param(y, "y")
  p1 <- l$p1; p2 <- l$p2
  if (identical(p1$coords[[1L]]$key, p2$coords[[1L]]$key)) return(.sub(x, p1$coords[[1L]]))
  if (identical(p1$coords[[2L]]$key, p2$coords[[2L]]$key)) return(.sub(y, p1$coords[[2L]]))
  cf <- geo_coefficients(l)
  .add(list(.mul(list(cf[[1L]], x)), .mul(list(cf[[2L]], y)), cf[[3L]]))
}

.mg_dirn <- function(l, i) {
  a <- l$p1$coords[[i]]; b <- l$p2$coords[[i]]
  if (.mg_rbool(Lt(a, b))) return(S_Infinity)
  if (identical(a$key, b$key)) return(S_Zero)
  S_NegativeInfinity
}

#' @rdname Line
#' @export
xdirection <- function(l) .mg_dirn(l, 1L)
#' @rdname Line
#' @export
ydirection <- function(l) .mg_dirn(l, 2L)
#' @rdname Line
#' @export
zdirection <- function(l) .mg_dirn(l, 3L)

#' @rdname Line
#' @export
closing_angle <- function(r1, r2) {
  if (!inherits(r1, "GeoRay2D") || !inherits(r2, "GeoRay2D")) .mg_err("TypeError", "Both arguments must be Ray2D objects.")
  d1 <- .mg_ldir(r1)$coords; d2 <- .mg_ldir(r2)$coords
  a1 <- atan2(d1[[2L]], d1[[1L]]); a2 <- atan2(d2[[2L]], d2[[1L]])
  if (.mg_rbool(Lt(.mul(list(a1, a2)), S_Zero))) {
    tp <- .mul(list(S_Two, S_Pi))
    if (.mg_rbool(Lt(a1, S_Zero))) a1 <- .add(list(tp, a1))
    if (.mg_rbool(Lt(a2, S_Zero))) a2 <- .add(list(tp, a2))
  }
  .sub(a1, a2)
}

#' @rdname Line
#' @export
bisectors <- function(l, l2 = NULL) {
  if (inherits(l, "GeoTriangle")) return(.mg_tri_bisectors(l))
  if (inherits(l, "GeoPolygon")) return(.mg_poly_bisectors(l, l2))
  if (!inherits(l2, "GeoLine")) .mg_err("GeometryError", "Expecting LinearEntity, not ", .mg_gstr(l2))
  l1 <- l
  if (length(l1$p1$coords) != length(l2$p1$coords)) {
    if (inherits(l1, "GeoLine2D")) { tmp <- l1; l1 <- l2; l2 <- tmp }
    p1 <- .mg_normdim(list(l1$p1, l2$p1), on_morph = "ignore")[[2L]]
    p2 <- .mg_normdim(list(l1$p2, l2$p2), on_morph = "ignore")[[2L]]
    l2 <- Line(p1, p2)
  }
  pt <- intersection(l1, l2)
  if (!length(pt)) .mg_err("GeometryError", "The lines do not intersect")
  pt <- pt[[1L]]
  if (inherits(pt, "GeoInfLine")) return(list(l))
  d1 <- unit(.mg_ldir(l1)); d2 <- unit(.mg_ldir(l2))
  list(Line(pt, .mg_pt_add(.mg_pt_add(pt, d1), d2)), Line(pt, .mg_pt_sub(.mg_pt_add(pt, d1), d2)))
}

# min and max of SymPy numbers as Python does them
.mg_min <- function(v) { r <- v[[1L]]; for (z in v[-1L]) if (.mg_rbool(Lt(z, r))) r <- z; r }
.mg_max <- function(v) { r <- v[[1L]]; for (z in v[-1L]) if (.mg_rbool(Gt(z, r))) r <- z; r }

#' @rdname Line
#' @export
bounds <- function(l) {
  if (inherits(l, "GeoPoint2D")) return(list(l$coords[[1L]], l$coords[[2L]], l$coords[[1L]], l$coords[[2L]]))
  if (inherits(l, "GeoEllipse")) {
    c <- l$center$coords; h <- l$hradius; v <- l$vradius
    return(list(.sub(c[[1L]], h), .sub(c[[2L]], v), .add(list(c[[1L]], h)), .add(list(c[[2L]], v))))
  }
  verts <- if (inherits(l, "GeoPolygon")) l$vertices else list(l$p1, l$p2)
  xs <- lapply(verts, point_x); ys <- lapply(verts, point_y)
  list(.mg_min(xs), .mg_min(ys), .mg_max(xs), .mg_max(ys))
}
