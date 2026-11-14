# Polygons, triangles and regular polygons, after sympy/geometry/polygon.py.

.mg_poly_raw <- function(gtype, verts) {
  cls <- switch(gtype, Triangle = c("GeoTriangle", "GeoPolygon"), "GeoPolygon")
  .mg_gnew(gtype, cls, list(vertices = verts))
}

# Python list helpers with negative indices
.mg_pyi <- function(i, n) if (i < 0) n + i + 1L else i + 1L

.mg_dedupe <- function(vs) {
  nod <- list()
  for (p in vs) {
    if (length(nod) && .mg_geo_eq(p, nod[[length(nod)]])) next
    nod[[length(nod) + 1L]] <- p
  }
  if (length(nod) > 1L && .mg_geo_eq(nod[[length(nod)]], nod[[1L]])) nod[[length(nod)]] <- NULL
  nod
}

.mg_dict <- function(keys, vals) {
  out <- vals
  names(out) <- vapply(keys, .mg_gstr, "")
  attr(out, "dict") <- TRUE
  out
}

#' Polygons and triangles
#'
#' `Polygon()` builds a polygon from its vertices as SymPy does: repeated
#' and collinear vertices are removed, so that three vertices give a
#' `Triangle`, two a segment and one a point; `n` gives a regular polygon.
#' `Triangle()` also accepts `sss`, `asa` or `sas` (side lengths and angles
#' in degrees). `RegularPolygon()` takes a center, a radius, the number of
#' sides and a rotation. The vertices are in `$vertices`.
#'
#' The functions are SymPy's methods and properties: `sides()`, `angles()`,
#' `perimeter()`, `is_convex()`, `encloses()`, `cut_section()`,
#' `first_moment_of_area()`, the triangle centers `circumcenter()`,
#' `incenter()`, `orthocenter()`, `excenters()` and the related
#' `circumradius()`, `circumcircle()`, `inradius()`, `incircle()`,
#' `exradii()`, `altitudes()`, `medians()`, `medial()`,
#' `nine_point_circle()`, `eulerline()`, `is_right()`, `is_equilateral()`,
#' `is_isosceles()`, `is_scalene()`, `is_similar()`, and for regular
#' polygons `apothem()`, `interior_angle()`, `exterior_angle()`,
#' `rotation()` and `spin()` (which returns the spun polygon). Results
#' keyed by vertices or sides are named lists.
#'
#' @param ... Vertices; for `Triangle()` also `sss`, `asa` or `sas`.
#' @param n The number of sides of a regular polygon.
#' @param c,r The center and radius of a regular polygon.
#' @param rot,angle Rotation angles.
#' @param p,t,t1,t2 Polygons, triangles or other entities.
#' @param line A line.
#' @param point A reference point.
#' @return An entity, an expression, a logical value or a list.
#' @examples
#' p <- Polygon(Point(0, 0), Point(1, 0), Point(1, 1), Point(0, 1))
#' area(p)
#' perimeter(p)
#' t <- Triangle(Point(0, 0), Point(4, 0), Point(0, 3))
#' circumcenter(t)
#' incenter(t)
#' Polygon(Point(0, 0), Point(1, 0), Point(2, 0), Point(1, 1))
#' @export
Polygon <- function(..., n = 0) {
  args <- list(...)
  if (!identical(n, 0) && !identical(n, 0L)) {
    if (length(args) == 2L) args <- c(args, list(n))
    else if (length(args) == 3L) args <- c(args[1:2], list(n), args[3L])
    return(do.call(RegularPolygon, args))
  }
  nod <- .mg_dedupe(lapply(args, function(a) .mg_as_pt(a, 2L)))
  i <- -3L
  while (i < length(nod) - 3L && length(nod) > 2L) {
    m <- length(nod)
    a <- nod[[.mg_pyi(i, m)]]; b <- nod[[.mg_pyi(i + 1L, m)]]; c <- nod[[.mg_pyi(i + 2L, m)]]
    if (is_collinear(a, b, c)) {
      nod[[.mg_pyi(i + 1L, m)]] <- NULL
      if (.mg_geo_eq(a, c)) nod[[.mg_pyi(i, length(nod))]] <- NULL
    } else i <- i + 1L
  }
  if (length(nod) > 3L) return(.mg_poly_raw("Polygon", nod))
  if (length(nod) == 3L) return(do.call(Triangle, nod))
  if (length(nod) == 2L) return(Segment(nod[[1L]], nod[[2L]]))
  nod[[1L]]
}

#' @rdname Polygon
#' @export
Triangle <- function(...) {
  args <- list(...)
  kw <- names(args)
  if (length(args) == 1L && !is.null(kw) && kw %in% c("sss", "asa", "sas")) {
    v <- lapply(lapply(args[[1L]], as_sym), simplify)
    return(switch(kw, sss = .mg_sss(v[[1L]], v[[2L]], v[[3L]]), asa = .mg_asa(v[[1L]], v[[2L]], v[[3L]]),
                  sas = .mg_sas(v[[1L]], v[[2L]], v[[3L]])))
  }
  if (length(args) != 3L) .mg_err("GeometryError", "Triangle instantiates with three points or a valid keyword.")
  nod <- .mg_dedupe(lapply(args, function(a) .mg_as_pt(a, 2L)))
  i <- -3L
  while (i < length(nod) - 3L && length(nod) > 2L) {
    m <- length(nod)
    tri <- list(nod[[.mg_pyi(i, m)]], nod[[.mg_pyi(i + 1L, m)]], nod[[.mg_pyi(i + 2L, m)]])
    tri <- tri[order(vapply(tri, .mg_ent_sort_key, ""), method = "radix")]
    if (is_collinear(tri[[1L]], tri[[2L]], tri[[3L]])) {
      nod[[.mg_pyi(i, m)]] <- tri[[1L]]
      nod[[.mg_pyi(i + 1L, m)]] <- NULL
    }
    i <- i + 1L
  }
  if (length(nod) == 3L) return(.mg_poly_raw("Triangle", nod))
  if (length(nod) == 2L) return(Segment(nod[[1L]], nod[[2L]]))
  nod[[1L]]
}

.mg_tan_deg <- function(d) tan(rad(d))

.mg_asa <- function(d1, l, d2) {
  xy <- .mg_line_intersection(Line(list(0, 0), slope = .mg_tan_deg(d1)), Line(list(l, 0), slope = .mg_tan_deg(.sub(.int(180), d2))))[[1L]]
  Triangle(list(0, 0), list(l, 0), xy)
}

.mg_sss <- function(l1, l2, l3) {
  c1 <- Circle(list(0, 0), l3)
  c2 <- Circle(list(l1, 0), l2)
  inter <- Filter(function(a) isTRUE(ask_fact(a$coords[[2L]], "nonnegative")), .mg_ellipse_intersection(c1, c2))
  if (!length(inter)) return(NULL)
  Triangle(list(0, 0), list(l1, 0), inter[[1L]])
}

.mg_sas <- function(l1, d, l2) {
  Triangle(Point(0, 0), Point(l2, 0), Point(.mul(list(cos(rad(d)), l1)), .mul(list(sin(rad(d)), l1))))
}

#' @rdname Polygon
#' @export
RegularPolygon <- function(c, r, n, rot = 0) {
  r <- as_sym(r); n <- as_sym(n); rot <- as_sym(rot)
  c <- .mg_as_pt(c, 2L)
  if (n$type %in% c("Integer", "Rational", "Float")) {
    if (!(n$type == "Integer" || n$type == "Float" && .mg_float_is_int(n))) .mg_err("ValueError", .str(n), " is not an integer")
    if (.mg_num(n) < 3) .mg_err("GeometryError", "n must be a >= 3, not ", .str(n))
  }
  if (.is_number(rot)) rot <- Mod(rot, .div(.mul(list(S_Two, S_Pi)), n))
  .mg_regpoly_raw(c, r, n, rot)
}

.mg_regpoly_raw <- function(c, r, n, rot) {
  e <- .mg_gnew("RegularPolygon", c("GeoRegularPolygon", "GeoPolygon"), list(center = c, radius = r, n = n, rot = rot))
  e$vertices <- .mg_regpoly_vertices(e)
  e
}

.mg_regpoly_vertices <- function(e) {
  if (e$n$type != "Integer") return(list())
  c <- e$center$coords; r <- Abs(e$radius); v <- .div(.mul(list(S_Two, S_Pi)), e$n)
  lapply(seq_len(as.integer(.mg_num(e$n))) - 1L, function(k) {
    ang <- .add(list(.mul(list(.int(k), v)), e$rot))
    Point(.add(list(c[[1L]], .mul(list(r, cos(ang))))), .add(list(c[[2L]], .mul(list(r, sin(ang))))))
  })
}

# ---------------------------------------------------------------------------

.mg_poly_area <- function(p) {
  if (inherits(p, "GeoRegularPolygon"))
    return(.div(.mul(list(sign(p$radius), p$n, .pow(.mg_regpoly_side(p), S_Two))), .mul(list(.int(4), tan(.div(S_Pi, p$n))))))
  v <- p$vertices; m <- length(v)
  terms <- lapply(seq_len(m), function(i) {
    a <- v[[if (i == 1L) m else i - 1L]]$coords; b <- v[[i]]$coords
    .sub(.mul(list(a[[1L]], b[[2L]])), .mul(list(b[[1L]], a[[2L]])))
  })
  .div(simplify(.add(terms)), S_Two)
}

.mg_regpoly_side <- function(p) .mul(list(p$radius, S_Two, sin(.div(S_Pi, p$n))))

.mg_is_clockwise <- function(a, b, c) {
  ba <- .mg_pt_sub(b, a); ca <- .mg_pt_sub(c, a)
  ta <- simplify(.sub(.mul(list(ba$coords[[1L]], ca$coords[[2L]])), .mul(list(ca$coords[[1L]], ba$coords[[2L]]))))
  res <- ask_fact(ta, "nonpositive")
  if (is.na(res)) .mg_err("ValueError", "Can't determine orientation")
  res
}

#' @rdname Polygon
#' @export
angles <- function(p) {
  v <- p$vertices; m <- length(v)
  if (inherits(p, "GeoRegularPolygon")) return(.mg_dict(v, rep(list(interior_angle(p)), m)))
  keys <- list(); vals <- list()
  for (i in seq_len(m) - 1L) {
    a <- v[[.mg_pyi(i - 2L, m)]]; b <- v[[.mg_pyi(i - 1L, m)]]; c <- v[[.mg_pyi(i, m)]]
    ra <- .mg_line_angle(Ray(b, a), Ray(b, c))
    keys[[length(keys) + 1L]] <- b
    vals[[length(vals) + 1L]] <- if (.mg_is_clockwise(a, b, c)) .sub(.mul(list(S_Two, S_Pi)), ra) else ra
  }
  wrong <- ask_fact(.sub(.div(.sub(.div(.add(vals), S_Pi), S_One), .int(m - 2L)), S_One), "positive")
  if (isTRUE(wrong)) vals <- lapply(vals, function(z) .sub(.mul(list(S_Two, S_Pi)), z))
  else if (is.na(wrong)) .mg_err("ValueError", "could not determine Polygon orientation.")
  .mg_dict(keys, vals)
}

#' @rdname Polygon
#' @export
perimeter <- function(p) {
  v <- p$vertices; m <- length(v)
  simplify(.add(lapply(seq_len(m), function(i) .mg_pt_distance(v[[if (i == 1L) m else i - 1L]], v[[i]]))))
}

.mg_poly_centroid <- function(p) {
  if (inherits(p, "GeoRegularPolygon")) return(p$center)
  A <- .mg_poly_area(p)
  if (is_zero_num(A)) .mg_err("GeometryError", "Centroid is undefined for a polygon with zero area.")
  Af <- .div(S_One, .mul(list(.int(6), A)))
  v <- p$vertices; m <- length(v)
  cx <- list(); cy <- list()
  for (i in seq_len(m)) {
    a <- v[[if (i == 1L) m else i - 1L]]$coords; b <- v[[i]]$coords
    w <- .sub(.mul(list(a[[1L]], b[[2L]])), .mul(list(b[[1L]], a[[2L]])))
    cx[[i]] <- .mul(list(w, .add(list(a[[1L]], b[[1L]]))))
    cy[[i]] <- .mul(list(w, .add(list(a[[2L]], b[[2L]]))))
  }
  Point(simplify(.mul(list(Af, .add(cx)))), simplify(.mul(list(Af, .add(cy)))))
}

.mg_poly_second_moment <- function(p, point = NULL) {
  v <- p$vertices; m <- length(v)
  Ixx <- list(); Iyy <- list(); Ixy <- list()
  for (i in seq_len(m)) {
    a <- v[[if (i == 1L) m else i - 1L]]$coords; b <- v[[i]]$coords
    x1 <- a[[1L]]; y1 <- a[[2L]]; x2 <- b[[1L]]; y2 <- b[[2L]]
    w <- .sub(.mul(list(x1, y2)), .mul(list(x2, y1)))
    Ixx[[i]] <- .mul(list(.add(list(.pow(y1, S_Two), .mul(list(y1, y2)), .pow(y2, S_Two))), w))
    Iyy[[i]] <- .mul(list(.add(list(.pow(x1, S_Two), .mul(list(x1, x2)), .pow(x2, S_Two))), w))
    Ixy[[i]] <- .mul(list(.add(list(.mul(list(x1, y2)), .mul(list(S_Two, x1, y1)), .mul(list(S_Two, x2, y2)), .mul(list(x2, y1)))), w))
  }
  A <- .mg_poly_area(p)
  cc <- .mg_poly_centroid(p)$coords
  cx <- cc[[1L]]; cy <- cc[[2L]]
  Ixxc <- .sub(.div(.add(Ixx), .int(12)), .mul(list(A, .pow(cy, S_Two))))
  Iyyc <- .sub(.div(.add(Iyy), .int(12)), .mul(list(A, .pow(cx, S_Two))))
  Ixyc <- .sub(.div(.add(Ixy), .int(24)), .mul(list(A, cx, cy)))
  if (is.null(point)) return(list(Ixxc, Iyyc, Ixyc))
  pc <- lapply(.mg_seq(point), as_sym)
  list(.add(list(Ixxc, .mul(list(A, .pow(.sub(pc[[2L]], cy), S_Two))))),
       .add(list(Iyyc, .mul(list(A, .pow(.sub(pc[[1L]], cx), S_Two))))),
       .add(list(Ixyc, .mul(list(A, .mul(list(.sub(pc[[1L]], cx), .sub(pc[[2L]], cy))))))))
}

#' @rdname Polygon
#' @export
first_moment_of_area <- function(p, point = NULL) {
  if (!is.null(point)) {
    cc <- .mg_poly_centroid(p)
    point <- .mg_as_pt(point, 2L)
  } else {
    point <- .mg_poly_centroid(p)
    cc <- point
  }
  xc <- cc$coords[[1L]]; yc <- cc$coords[[2L]]
  hp <- cut_section(p, Line(point, slope = 0))
  vp <- cut_section(p, Line(point, slope = S_Infinity))
  le <- function(a, b) .mg_rbool(Le(.mg_poly_area(a), .mg_poly_area(b)))
  p1 <- if (le(hp[[1L]], hp[[2L]])) hp[[1L]] else hp[[2L]]
  p2 <- if (le(vp[[1L]], vp[[2L]])) vp[[1L]] else vp[[2L]]
  list(.mul(list(.sub(.mg_poly_centroid(p1)$coords[[2L]], yc), .mg_poly_area(p1))),
       .mul(list(.sub(.mg_poly_centroid(p2)$coords[[1L]], xc), .mg_poly_area(p2))))
}

#' @rdname Polygon
#' @export
sides <- function(p) {
  v <- p$vertices; m <- length(v)
  lapply(seq_len(m), function(i) Segment(v[[i]], v[[if (i == m) 1L else i + 1L]]))
}

.mg_poly_is_convex <- function(p) {
  v <- p$vertices; m <- length(v)
  cw <- .mg_is_clockwise(v[[m - 1L]], v[[m]], v[[1L]])
  for (i in seq_len(m - 1L)) {
    if (xor(cw, .mg_is_clockwise(v[[.mg_pyi(i - 2L, m)]], v[[.mg_pyi(i - 1L, m)]], v[[i + 1L]]))) return(FALSE)
  }
  s <- sides(p); ns <- length(s)
  for (i in seq_len(ns) - 1L) {
    si <- s[[i + 1L]]
    pts <- c(.mg_pkey(si$p1), .mg_pkey(si$p2))
    lo <- if (i == ns - 1L) 1L else 0L
    if (lo < i - 1L) for (j in lo:(i - 2L)) {
      sj <- s[[j + 1L]]
      if (!(.mg_pkey(sj$p1) %in% pts) && !(.mg_pkey(sj$p2) %in% pts)) {
        if (length(.mg_isect1(si, sj))) return(FALSE)
      }
    }
  }
  TRUE
}

.mg_poly_encloses <- function(p, pt) {
  pt <- .mg_as_pt(pt, 2L)
  if (.mg_pkey(pt) %in% vapply(p$vertices, .mg_pkey, "") || any(vapply(sides(p), function(s) .mg_line_contains(s, pt), TRUE)))
    return(FALSE)
  lit <- list()
  for (v in p$vertices) {
    lit[[length(lit) + 1L]] <- .mg_pt_sub(v, pt)
    if (length(.mg_geo_free(lit[[length(lit)]]))) return(NA)
  }
  poly <- do.call(Polygon, lit)
  args <- poly$vertices; m <- length(args)
  idx <- seq(-m, 0L)
  if (.mg_poly_is_convex(poly)) {
    orient <- NULL
    for (i in idx) {
      a <- args[[.mg_pyi(i, m)]]; b <- args[[.mg_pyi(i + 1L, m)]]
      ax <- a$coords[[1L]]; ay <- a$coords[[2L]]; bx <- b$coords[[1L]]; by <- b$coords[[2L]]
      test <- ask_fact(.sub(.mul(list(.neg(ay), .sub(bx, ax))), .mul(list(.neg(ax), .sub(by, ay)))), "negative")
      if (is.null(orient)) orient <- test
      else if (!identical(test, orient)) return(FALSE)
    }
    return(TRUE)
  }
  hit <- FALSE
  p1x <- args[[1L]]$coords[[1L]]; p1y <- args[[1L]]$coords[[2L]]
  for (i in idx[-1L]) {
    p2 <- args[[.mg_pyi(i, m)]]
    p2x <- p2$coords[[1L]]; p2y <- p2$coords[[2L]]
    if (.mg_rbool(Gt(S_Zero, .mg_min(list(p1y, p2y)))) && .mg_rbool(Le(S_Zero, .mg_max(list(p1y, p2y)))) &&
        .mg_rbool(Le(S_Zero, .mg_max(list(p1x, p2x)))) && !identical(p1y$key, p2y$key)) {
      xin <- .add(list(.div(.mul(list(.neg(p1y), .sub(p2x, p1x))), .sub(p2y, p1y)), p1x))
      if (identical(p1x$key, p2x$key) || .mg_rbool(Le(S_Zero, xin))) hit <- !hit
    }
    p1x <- p2x; p1y <- p2y
  }
  hit
}

.mg_regpoly_encloses <- function(e, p) {
  p <- .mg_as_pt(p, 2L)
  d <- geo_length(Segment(e$center, p))
  if (.mg_rbool(Ge(d, e$radius))) return(FALSE)
  if (.mg_rbool(Lt(d, apothem(e)))) return(TRUE)
  .mg_poly_encloses(e, p)
}

.mg_poly_arbitrary <- function(p, parameter = "t") {
  t <- .mg_param(parameter)
  .mg_check_param(p, t)
  per <- perimeter(p)
  start <- S_Zero
  pieces <- list()
  for (s in sides(p)) {
    frac <- .div(geo_length(s), per)
    end <- .add(list(start, frac))
    pt <- .mg_geo_subs(.mg_line_arbitrary(s, t), t, .div(.sub(t, start), frac))
    pieces[[length(pieces) + 1L]] <- list(pt, And(Le(start, t), Lt(t, end)))
    start <- end
  }
  pieces
}

#' @rdname Polygon
#' @export
cut_section <- function(p, line) {
  ip <- .mg_isect1(p, line)
  if (!length(ip)) .mg_err("ValueError", "This line does not intersect the polygon")
  pts <- c(p$vertices, list(p$vertices[[1L]]))
  xx <- Dummy("polygon_dummy", real = TRUE); yy <- Dummy("polygon_dummy", real = TRUE)
  eq <- .mg_line_equation(line, xx, yy)
  a <- coeff(eq, xx); b <- coeff(eq, yy)
  up <- list(); lo <- list()
  prev <- TRUE; prevp <- NULL
  for (pt in pts) {
    cmp <- if (.mg_truthy(b)) .div(subs(eq, list(list(xx, pt$coords[[1L]]), list(yy, pt$coords[[2L]]))), b)
           else .div(subs(eq, xx, pt$coords[[1L]]), a)
    if (.mg_rbool(Gt(cmp, S_Zero))) {
      if (!prev) {
        np <- .mg_line_intersection(Line(pt, prevp), line)
        up <- c(up, list(np[[1L]])); lo <- c(lo, list(np[[1L]]))
      }
      up <- c(up, list(pt))
      prev <- TRUE
    } else {
      if (prev && !is.null(prevp)) {
        np <- .mg_line_intersection(Line(pt, prevp), line)
        up <- c(up, list(np[[1L]])); lo <- c(lo, list(np[[1L]]))
      }
      lo <- c(lo, list(pt))
      prev <- FALSE
    }
    prevp <- pt
  }
  upp <- if (length(up)) do.call(Polygon, up) else NULL
  lop <- if (length(lo)) do.call(Polygon, lo) else NULL
  list(if (inherits(upp, "GeoPolygon")) upp else NULL, if (inherits(lop, "GeoPolygon")) lop else NULL)
}

.mg_poly_intersection <- function(p, o) {
  k <- if (inherits(o, "GeoPolygon")) sides(o) else list(o)
  res <- list()
  for (s in sides(p)) for (s1 in k) res <- c(res, .mg_isect1(s, s1))
  res <- res[!duplicated(vapply(res, .mg_geo_key, ""))]
  pts <- Filter(function(e) inherits(e, "GeoPoint"), res)
  segs <- Filter(function(e) inherits(e, "GeoSegment"), res)
  if (length(pts) && length(segs)) {
    inseg <- vapply(pts, function(q) any(vapply(segs, function(s) .mg_line_contains(s, q), TRUE)), TRUE)
    if (any(inseg)) return(.mg_sort_ents(c(segs, pts[!inseg])))
    return(.mg_sort_ents(c(segs, pts)))
  }
  .mg_sort_ents(res)
}

.mg_poly_has <- function(p, o) {
  if (inherits(o, "GeoPolygon")) return(.mg_geo_eq(p, o))
  if (inherits(o, "GeoSegment")) return(any(vapply(sides(p), function(s) .mg_line_contains(s, o), TRUE)))
  if (inherits(o, "GeoPoint")) {
    if (.mg_pkey(o) %in% vapply(p$vertices, .mg_pkey, "")) return(TRUE)
    for (s in sides(p)) if (.mg_line_contains(s, o)) return(TRUE)
  }
  FALSE
}

.mg_poly_distance <- function(p, o) {
  if (inherits(o, "GeoPoint")) {
    d <- S_Infinity
    for (s in sides(p)) {
      cur <- .mg_line_distance(s, o)
      if (is_zero_num(cur)) return(S_Zero)
      if (.mg_rbool(Lt(cur, d))) d <- cur
    }
    return(d)
  }
  if (inherits(o, "GeoPolygon") && .mg_poly_is_convex(p) && .mg_poly_is_convex(o)) return(.mg_poly_poly_distance(p, o))
  .mg_err("NotImplementedError", "")
}

.mg_lt_n <- function(a, b) as.numeric(.neval(a, list())) < as.numeric(.neval(b, list()))

# Polygon._do_poly_distance (rotating calipers)
.mg_poly_poly_distance <- function(e1, e2) {
  c1 <- .mg_poly_centroid(e1); c2 <- .mg_poly_centroid(e2)
  r1 <- S_Zero; r2 <- S_Zero
  for (v in e1$vertices) { r <- .mg_pt_distance(c1, v); if (.mg_rbool(Lt(r1, r))) r1 <- r }
  for (v in e2$vertices) { r <- .mg_pt_distance(c2, v); if (.mg_rbool(Lt(r2, r))) r2 <- r }
  if (.mg_rbool(Le(.mg_pt_distance(c1, c2), .add(list(r1, r2))))) warning("Polygons may intersect producing erroneous output", call. = FALSE)
  ymax <- Point(0, S_NegativeInfinity); ymin <- Point(0, S_Infinity)
  for (v in e1$vertices) if (.mg_rbool(Gt(point_y(v), point_y(ymax))) || (identical(point_y(v)$key, point_y(ymax)$key) && .mg_rbool(Gt(point_x(v), point_x(ymax))))) ymax <- v
  for (v in e2$vertices) if (.mg_rbool(Lt(point_y(v), point_y(ymin))) || (identical(point_y(v)$key, point_y(ymin)$key) && .mg_rbool(Lt(point_x(v), point_x(ymin))))) ymin <- v
  mind <- .mg_pt_distance(ymax, ymin)
  conn <- function(e) {
    cn <- list()
    for (s in sides(e)) {
      k1 <- .mg_pkey(s$p1); k2 <- .mg_pkey(s$p2)
      cn[[k1]] <- c(cn[[k1]], list(s$p2)); cn[[k2]] <- c(cn[[k2]], list(s$p1))
    }
    cn
  }
  cn1 <- conn(e1); cn2 <- conn(e2)
  sup <- Line(Point(0, 0), Point(1, 0))
  pick <- function(cn, v, cmp) {
    a <- cn[[.mg_pkey(v)]][[1L]]; b <- cn[[.mg_pkey(v)]][[2L]]
    g1 <- .mg_line_angle(sup, Line(v, a)); g2 <- .mg_line_angle(sup, Line(v, b))
    if (.mg_rbool(cmp(g1, g2))) return(a)
    if (.mg_rbool(cmp(g2, g1))) return(b)
    if (.mg_rbool(Gt(.mg_pt_distance(v, a), .mg_pt_distance(v, b)))) b else a
  }
  n1 <- pick(cn1, ymax, Lt); n2 <- pick(cn2, ymin, Gt)
  cur1 <- ymax; cur2 <- ymin
  step <- function(cn, cur, nxt) {
    l <- cn[[.mg_pkey(nxt)]]
    list(nxt, if (!.mg_geo_eq(l[[1L]], cur)) l[[1L]] else l[[2L]])
  }
  repeat {
    a1 <- .mg_line_angle(sup, Line(cur1, n1))
    a2 <- .sub(S_Pi, .mg_line_angle(sup, Line(cur2, n2)))
    if (isTRUE(.mg_fz(Lt(a1, a2)))) {
      sup <- Line(cur1, n1)
      d <- .mg_line_distance(Segment(cur1, n1), cur2)
      if (.mg_lt_n(d, mind)) mind <- d
      s <- step(cn1, cur1, n1); cur1 <- s[[1L]]; n1 <- s[[2L]]
    } else if (isTRUE(.mg_fz(Gt(a1, a2)))) {
      sup <- Line(n2, cur2)
      d <- .mg_line_distance(Segment(cur2, n2), cur1)
      if (.mg_lt_n(d, mind)) mind <- d
      s <- step(cn2, cur2, n2); cur2 <- s[[1L]]; n2 <- s[[2L]]
    } else {
      sup <- Line(cur1, n1)
      m1 <- .mg_line_distance(Segment(cur1, n1), n2)
      m2 <- .mg_line_distance(Segment(cur2, n2), n1)
      d <- if (.mg_rbool(Lt(m2, m1))) m2 else m1
      if (.mg_lt_n(d, mind)) mind <- d
      s <- step(cn1, cur1, n1); cur1 <- s[[1L]]; n1 <- s[[2L]]
      s <- step(cn2, cur2, n2); cur2 <- s[[1L]]; n2 <- s[[2L]]
    }
    if (.mg_geo_eq(cur1, ymax) && .mg_geo_eq(cur2, ymin)) break
  }
  mind
}

# Polygon._hashable_content: the vertices from the least rotation, in the
# smaller of the two directions
.mg_poly_canonical <- function(verts) {
  ref <- function(pl) {
    keys <- vapply(pl, .mg_pkey, "")
    u <- !duplicated(keys)
    up <- pl[u]
    o <- order(vapply(up, .mg_ent_sort_key, ""), method = "radix")
    rk <- match(keys, keys[u][o]) - 1L
    list(r = rk, d = up[o])
  }
  least_rot <- function(r) {
    n <- length(r)
    best <- 0L
    for (k in seq_len(n - 1L)) {
      a <- r[c((k + 1L):n, seq_len(k))]
      b <- r[c((best + 1L):n, seq_len(best))[seq_len(n)]]
      cmp <- which(a != b)
      if (length(cmp) && a[cmp[1L]] < b[cmp[1L]]) best <- k
    }
    if (best) r[c((best + 1L):n, seq_len(best))] else r
  }
  s1 <- ref(verts); s2 <- ref(rev(verts))
  rn <- least_rot(s1$r); rr <- least_rot(s2$r)
  dif <- which(rn != rr)
  r <- if (length(dif) && rn[dif[1L]] < rr[dif[1L]]) rn else rr
  lapply(r, function(i) s1$d[[i + 1L]])
}

.mg_poly_bisectors <- function(p, prec = NULL) {
  pts <- c(p$vertices, list(p$vertices[[1L]]))
  if (.mg_is_clockwise(pts[[1L]], pts[[2L]], pts[[3L]])) pts <- rev(pts)
  an <- angles(p)
  keys <- list(); vals <- list()
  for (v in p$vertices) {
    a <- an[[format(v)]]
    i <- match(.mg_pkey(v), vapply(pts, .mg_pkey, ""))
    ray <- .mg_entity_rotate(Ray(pts[[i]], pts[[i + 1L]]), .div(a, S_Two), v)
    d <- .mg_ldir(ray)
    ray <- Ray(ray$p1, .mg_pt_add(ray$p1, .mg_pt_div(d, .mg_pt_distance(d, Point(0, 0)))))
    if (!is.null(prec)) ray <- Ray(ray$p1, .mg_geo_evalf(ray$p2, prec))
    keys[[length(keys) + 1L]] <- v; vals[[length(vals) + 1L]] <- ray
  }
  .mg_dict(keys, vals)
}

# ---------------------------------------------------------------------------
# Regular polygons

#' @rdname Polygon
#' @export
apothem <- function(p) .mul(list(p$radius, cos(.div(S_Pi, p$n))))
#' @rdname Polygon
#' @export
interior_angle <- function(p) .div(.mul(list(.sub(p$n, S_Two), S_Pi)), p$n)
#' @rdname Polygon
#' @export
exterior_angle <- function(p) .div(.mul(list(S_Two, S_Pi)), p$n)
#' @rdname Polygon
#' @export
rotation <- function(p) p$rot
#' @rdname Polygon
#' @export
spin <- function(p, angle) .mg_regpoly_raw(p$center, p$radius, p$n, .add(list(p$rot, as_sym(angle))))

# ---------------------------------------------------------------------------
# Triangles

.mg_side_lengths <- function(t) lapply(sides(t), geo_length)

#' @rdname Polygon
#' @export
is_similar <- function(t1, t2) {
  if (inherits(t1, "GeoTriangle")) {
    if (!inherits(t2, "GeoPolygon")) return(FALSE)
    s1 <- .mg_side_lengths(t1); s2 <- .mg_side_lengths(t2)
    sim <- function(u, v) {
      e <- Map(function(a, b) simplify(.div(a, b)), u, v)
      identical(e[[1L]]$key, e[[2L]]$key) && identical(e[[2L]]$key, e[[3L]]$key)
    }
    perms <- list(c(1, 2, 3), c(1, 3, 2), c(2, 1, 3), c(2, 3, 1), c(3, 1, 2), c(3, 2, 1))
    return(any(vapply(perms, function(pm) sim(s1[pm], s2), TRUE)))
  }
  if (inherits(t1, "GeoLine")) return(.mg_line_contains(Line(t1$p1, t1$p2), t2))
  if (inherits(t1, "GeoCircle") && inherits(t2, "GeoCircle")) return(TRUE)
  .mg_err("NotImplementedError", "")
}

.mg_lkeys <- function(t) vapply(.mg_side_lengths(t), function(z) .k(z), "")

#' @rdname Polygon
#' @export
is_equilateral <- function(t) length(unique(.mg_lkeys(t))) <= 1L
#' @rdname Polygon
#' @export
is_isosceles <- function(t) anyDuplicated(.mg_lkeys(t)) > 0L
#' @rdname Polygon
#' @export
is_scalene <- function(t) anyDuplicated(.mg_lkeys(t)) == 0L

#' @rdname Polygon
#' @export
is_right <- function(t) {
  s <- sides(t)
  isTRUE(.mg_is_perpendicular(s[[1L]], s[[2L]])) || isTRUE(.mg_is_perpendicular(s[[2L]], s[[3L]])) ||
    isTRUE(.mg_is_perpendicular(s[[1L]], s[[3L]]))
}

#' @rdname Polygon
#' @export
altitudes <- function(t) {
  s <- sides(t); v <- t$vertices
  .mg_dict(v, list(perpendicular_segment(s[[2L]], v[[1L]]), perpendicular_segment(s[[3L]], v[[2L]]),
                   perpendicular_segment(s[[1L]], v[[3L]])))
}

#' @rdname Polygon
#' @export
orthocenter <- function(t) {
  a <- altitudes(t)
  .mg_line_intersection(Line(a[[1L]]), Line(a[[2L]]))[[1L]]
}

#' @rdname Polygon
#' @export
circumcenter <- function(t) {
  if (inherits(t, "GeoRegularPolygon")) return(t$center)
  b <- lapply(sides(t), perpendicular_bisector)
  .mg_line_intersection(b[[1L]], b[[2L]])[[1L]]
}

#' @rdname Polygon
#' @export
circumradius <- function(t) {
  if (inherits(t, "GeoRegularPolygon")) return(t$radius)
  .mg_pt_distance(circumcenter(t), t$vertices[[1L]])
}

#' @rdname Polygon
#' @export
circumcircle <- function(t) {
  if (inherits(t, "GeoRegularPolygon")) return(Circle(t$center, t$radius))
  Circle(circumcenter(t), circumradius(t))
}

.mg_tri_bisectors <- function(t) {
  s <- lapply(sides(t), Line); v <- t$vertices; c <- incenter(t)
  seg <- function(i, j) Segment(v[[i]], .mg_line_intersection(Line(v[[i]], c), s[[j]])[[1L]])
  .mg_dict(v, list(seg(1L, 2L), seg(2L, 3L), seg(3L, 1L)))
}

#' @rdname Polygon
#' @export
incenter <- function(t) {
  s <- sides(t)
  l <- lapply(c(2L, 3L, 1L), function(i) geo_length(s[[i]]))
  p <- .add(l)
  v <- t$vertices
  x <- simplify(.div(.add(Map(function(a, b) .mul(list(a, point_x(b))), l, v)), p))
  y <- simplify(.div(.add(Map(function(a, b) .mul(list(a, point_y(b))), l, v)), p))
  Point(x, y)
}

#' @rdname Polygon
#' @export
inradius <- function(t) {
  if (inherits(t, "GeoRegularPolygon")) return(apothem(t))
  simplify(.div(.mul(list(S_Two, .mg_poly_area(t))), perimeter(t)))
}

#' @rdname Polygon
#' @export
incircle <- function(t) {
  if (inherits(t, "GeoRegularPolygon")) return(Circle(t$center, apothem(t)))
  Circle(incenter(t), inradius(t))
}

#' @rdname Polygon
#' @export
exradii <- function(t) {
  s <- sides(t)
  a <- geo_length(s[[1L]]); b <- geo_length(s[[2L]]); c <- geo_length(s[[3L]])
  hs <- .div(.add(list(a, b, c)), S_Two)
  A <- .mg_poly_area(t)
  .mg_dict(s, lapply(list(a, b, c), function(z) simplify(.div(A, .sub(hs, z)))))
}

#' @rdname Polygon
#' @export
excenters <- function(t) {
  s <- sides(t); v <- t$vertices
  a <- geo_length(s[[2L]]); b <- geo_length(s[[3L]]); c <- geo_length(s[[1L]])
  X <- lapply(v, point_x); Y <- lapply(v, point_y)
  f <- function(w, ca, cb, cc) simplify(.div(.add(list(.mul(list(ca, a, w[[1L]])), .mul(list(cb, b, w[[2L]])), .mul(list(cc, c, w[[3L]])))),
                                              .add(list(.mul(list(ca, a)), .mul(list(cb, b)), .mul(list(cc, c))))))
  m1 <- S_NegativeOne; p1 <- S_One
  .mg_dict(s, list(Point(f(X, m1, p1, p1), f(Y, m1, p1, p1)), Point(f(X, p1, m1, p1), f(Y, p1, m1, p1)),
                   Point(f(X, p1, p1, m1), f(Y, p1, p1, m1))))
}

#' @rdname Polygon
#' @export
medians <- function(t) {
  s <- sides(t); v <- t$vertices
  .mg_dict(v, list(Segment(v[[1L]], midpoint(s[[2L]])), Segment(v[[2L]], midpoint(s[[3L]])), Segment(v[[3L]], midpoint(s[[1L]]))))
}

#' @rdname Polygon
#' @export
medial <- function(t) {
  s <- sides(t)
  Triangle(midpoint(s[[1L]]), midpoint(s[[2L]]), midpoint(s[[3L]]))
}

#' @rdname Polygon
#' @export
nine_point_circle <- function(t) do.call(Circle, medial(t)$vertices)

#' @rdname Polygon
#' @export
eulerline <- function(t) {
  if (is_equilateral(t)) return(orthocenter(t))
  Line(orthocenter(t), circumcenter(t))
}
