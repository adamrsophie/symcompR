# Geometric entities and points, after sympy/geometry/entity.py and
# sympy/geometry/point.py.
#
# Entities are S3 lists. Every entity has the field `gtype` (the SymPy class
# name used for printing) and the classes below, ending in "GeoEntity":
#   Point2D, Point3D, Point            GeoPoint (coords)
#   Line2D, Segment3D, Ray2D, ...      GeoLine (p1, p2) with GeoInfLine,
#                                      GeoSegment or GeoRay
#   Ellipse, Circle                    GeoEllipse (center, hradius, vradius),
#                                      Circle also radius
#   Polygon, Triangle, RegularPolygon  GeoPolygon (vertices)
#   Plane                              GeoPlane (p1, normal_vector)
#   Curve                              GeoCurve (functions, limits)
#   Parabola                           GeoParabola (focus, directrix)

.mg_err <- function(cls, ...) stop(.sympy_condition(cls, paste0(...)))
# filldedent(): a leading newline, then the text filled to 70 columns
.mg_errf <- function(cls, ...) .mg_err(cls, "\n", paste(strwrap(paste0(...), width = 71), collapse = "\n"))

# The Python class name SymPy shows for an object
.mg_pytype <- function(x) {
  if (inherits(x, "GeoEntity")) return(x$gtype)
  if (is_sym(x)) return(x$type)
  if (is.numeric(x) && length(x) == 1L) return(if (x == round(x)) "int" else "float")
  class(x)[1L]
}

.mg_is_ent <- function(x) inherits(x, "GeoEntity")
.mg_is_pt <- function(x) inherits(x, "GeoPoint")

# Python truthiness of a SymPy expression: only a zero number is false
.mg_truthy <- function(e) {
  e <- as_sym(e)
  !(e$type %in% c("Integer", "Rational", "Float") && .mg_num(e) == 0) &&
    !(e$type == "BooleanAtom" && e$name == "false")
}

# bool() of a SymPy relational: TypeError when undecided
.mg_rbool <- function(r) {
  if (is.logical(r)) return(isTRUE(r))
  if (is_sym(r) && r$type == "BooleanAtom") return(r$name == "true")
  .mg_err("TypeError", "cannot determine truth value of Relational")
}

# fuzzy value of a relational: TRUE, FALSE or NA
.mg_fz <- function(r) {
  if (is.logical(r)) return(r[1L])
  if (is_sym(r) && r$type == "BooleanAtom") return(r$name == "true")
  NA
}

.mg_seq <- function(x) {
  if (inherits(x, "GeoPoint")) return(x$coords)
  if (is_sym(x)) {
    if (x$type == "Tuple") return(x$args)
    return(NULL)
  }
  if (inherits(x, "symmatrix")) return(x$data)
  if (is.list(x) && !inherits(x, "GeoEntity")) return(x)
  if (is.atomic(x) && length(x) > 1L) return(as.list(x))
  NULL
}

.mg_pkey <- function(p) paste(vapply(p$coords, .k, ""), collapse = "|")

.mg_gnew <- function(gtype, cls, fields) structure(c(list(gtype = gtype), fields), class = c(cls, "GeoEntity"))

# ---------------------------------------------------------------------------
# Points

#' Points
#'
#' `Point()` builds a point of any dimension, as SymPy's `Point`: two
#' coordinates give a `Point2D`, three a `Point3D`. Float coordinates are
#' replaced by rationals unless `evaluate = FALSE`; `dim` pads with zeros or
#' drops trailing zero coordinates. Points support `+`, `-`, `*` and `/`
#' (each coordinate is simplified, as in SymPy) and `==`, which compares the
#' coordinates structurally.
#'
#' The functions below are SymPy's `Point` methods and properties:
#' `point_x()`, `point_y()` and `point_z()` (the `x`, `y` and `z`
#' properties), `ambient_dimension()`, `affine_rank()`, `is_collinear()`,
#' `are_collinear()`, `is_concyclic()`, `is_scalar_multiple()`,
#' `midpoint()`, `origin()`, `orthogonal_direction()`, `project()`,
#' `taxicab_distance()`, `canberra_distance()`, `unit()`,
#' `direction_ratio()`, `direction_cosine()` and `geo_transform()` (SymPy's
#' `transform`).
#'
#' @param ... Coordinates, a list of coordinates, or a point; for
#'   `affine_rank()`, `is_collinear()`, `are_collinear()` and
#'   `is_concyclic()`, points.
#' @param evaluate Replace floats by rationals.
#' @param dim The dimension of the point.
#' @param on_morph What to do when `dim` changes the dimension: `"ignore"`,
#'   `"warn"` or `"error"`.
#' @param p,q,a,b Points.
#' @param matrix A 3 by 3 (or, in 3D, 4 by 4) transformation matrix.
#' @return A point, an expression, a logical value or a list.
#' @examples
#' A <- Point(0, 0); B <- Point(3, 4)
#' distance(A, B)
#' A + B
#' Point(0.5, 0.25)
#' is_collinear(Point(0, 0), Point(1, 1), Point(2, 2))
#' midpoint(A, B)
#' @export
Point <- function(..., evaluate = TRUE, dim = NULL, on_morph = "ignore") {
  a <- list(...)
  coords <- if (length(a) == 1L) a[[1L]] else a
  if (inherits(coords, "GeoPoint")) {
    if (is.null(dim) || length(coords$coords) == dim) return(coords)
    evaluate <- FALSE
  }
  sq <- .mg_seq(coords)
  if (is.null(sq) && !(is.list(coords) && !length(coords))) {
    .mg_err("TypeError", "\nExpecting sequence of coordinates, not `", .mg_pytype(coords), "`")
  }
  if (is.null(sq)) sq <- list()
  if (!length(sq) && !is.null(dim)) sq <- rep(list(S_Zero), dim)
  coords <- lapply(sq, function(z) if (is_sym(z)) z else as_sym(z))
  if (is.null(dim)) dim <- length(coords)
  if (length(coords) < 2L)
    .mg_err("ValueError", "\nPoint requires 2 or more coordinates or keyword `dim` > 1.")
  if (length(coords) != dim) {
    msg <- sprintf("Dimension of (%s) needs to be changed from %d to %d.",
                   paste(vapply(coords, .str, ""), collapse = ", "), length(coords), dim)
    if (on_morph == "error") .mg_err("ValueError", msg)
    else if (on_morph == "warn") warning(msg, call. = FALSE)
    else if (on_morph != "ignore") .mg_err("ValueError", "\non_morph value should be 'error', 'warn' or 'ignore'.")
  }
  if (length(coords) > dim && any(vapply(coords[(dim + 1L):length(coords)], .mg_truthy, TRUE)))
    .mg_err("ValueError", "Nonzero coordinates cannot be removed.")
  for (z in coords) {
    if (.is_number(z) && identical(ask_fact(im(z), "zero"), FALSE))
      .mg_err("ValueError", "Imaginary coordinates are not permitted.")
  }
  coords <- c(coords, rep(list(S_Zero), max(0L, dim - length(coords))))[seq_len(dim)]
  if (isTRUE(evaluate)) coords <- lapply(coords, .mg_rationalize)
  .mg_point(coords)
}

# replace the floats of an expression by simplify(nsimplify(f, rational = TRUE))
.mg_rationalize <- function(z) {
  fl <- atoms(z, "Float")
  if (!length(fl)) return(z)
  xreplace(z, lapply(fl, function(f) list(f, simplify(nsimplify(f, rational = TRUE)))))
}

.mg_point <- function(coords) {
  n <- length(coords)
  if (n == 2L) return(.mg_gnew("Point2D", c("GeoPoint2D", "GeoPoint"), list(coords = coords)))
  if (n == 3L) return(.mg_gnew("Point3D", c("GeoPoint3D", "GeoPoint"), list(coords = coords)))
  .mg_gnew("Point", "GeoPoint", list(coords = coords))
}

#' @rdname Point
#' @export
Point2D <- function(..., evaluate = TRUE) Point(..., evaluate = evaluate, dim = 2L)

#' @rdname Point
#' @export
Point3D <- function(..., evaluate = TRUE) Point(..., evaluate = evaluate, dim = 3L)

.mg_as_pt <- function(p, dim = NULL, evaluate = TRUE) {
  if (inherits(p, "GeoPoint") && (is.null(dim) || length(p$coords) == dim)) return(p)
  Point(p, dim = dim, evaluate = evaluate)
}

#' @rdname Point
#' @export
point_x <- function(p) p$coords[[1L]]
#' @rdname Point
#' @export
point_y <- function(p) p$coords[[2L]]
#' @rdname Point
#' @export
point_z <- function(p) p$coords[[3L]]

# Point._normalize_dimension
.mg_normdim <- function(pts, dim = NULL, on_morph = "warn") {
  if (is.null(dim)) dim <- max(vapply(pts, function(p) length(p$coords), 1L))
  if (all(vapply(pts, function(p) length(p$coords) == dim, TRUE))) return(pts)
  lapply(pts, function(p) Point(p, dim = dim, on_morph = on_morph))
}

.mg_pt_add <- function(a, b) {
  o <- tryCatch(Point(b, evaluate = FALSE), TypeError = function(e) NULL, ValueError = function(e) NULL)
  if (is.null(o)) .mg_err("GeometryError", "Don't know how to add ", .mg_gstr(b), " and a Point object")
  pr <- .mg_normdim(list(a, o))
  Point(Map(function(u, v) simplify(.add(list(u, v))), pr[[1L]]$coords, pr[[2L]]$coords), evaluate = FALSE)
}
.mg_pt_neg <- function(a) Point(lapply(a$coords, .neg), evaluate = FALSE)
.mg_pt_sub <- function(a, b) {
  sq <- .mg_seq(b)
  if (is.null(sq)) .mg_err("TypeError", "cannot subtract a non sequence from a Point")
  .mg_pt_add(a, lapply(sq, function(z) .neg(as_sym(z))))
}
.mg_pt_mul <- function(a, f) {
  f <- as_sym(f)
  Point(lapply(a$coords, function(z) simplify(.mul(list(z, f)))), evaluate = FALSE)
}
.mg_pt_div <- function(a, f) {
  f <- as_sym(f)
  Point(lapply(a$coords, function(z) simplify(.div(z, f))), evaluate = FALSE)
}

#' @export
chooseOpsMethod.GeoEntity <- function(x, y, mx, my, cl, reverse) TRUE

#' Operators on geometric entities
#'
#' `==` and `!=` compare entities as SymPy does: points, lines and circles
#' structurally, polygons up to rotation and reversal of their vertices.
#' Points also support `+`, `-`, `*`, `/` and unary minus.
#'
#' @param e1,e2 Entities, or a point and a scalar or a coordinate list.
#' @return An entity or a logical value.
#' @examples
#' Point(1, 2) + Point(3, 4)
#' Point(1, 2) * 3
#' Line(Point(0, 0), Point(1, 1)) == Line(Point(0, 0), Point(2, 2))
#' @name geometry-operators
#' @export
Ops.GeoEntity <- function(e1, e2) {
  if (missing(e2)) {
    if (.Generic == "-" && inherits(e1, "GeoPoint")) return(.mg_pt_neg(e1))
    if (.Generic == "+") return(e1)
    .mg_err("TypeError", "bad operand type for unary ", .Generic, ": '", e1$gtype, "'")
  }
  if (.Generic == "==") return(.mg_geo_eq(e1, e2))
  if (.Generic == "!=") return(!.mg_geo_eq(e1, e2))
  p1 <- inherits(e1, "GeoPoint"); p2 <- inherits(e2, "GeoPoint")
  if (p1 || p2) {
    switch(.Generic,
           "+" = return(if (p1) .mg_pt_add(e1, e2) else .mg_pt_add(e2, e1)),
           "-" = return(if (p1) .mg_pt_sub(e1, e2) else .mg_pt_add(.mg_pt_neg(e2), e1)),
           "*" = {
             if (p1 && p2) .mg_err("TypeError", "cannot multiply two Points")
             return(if (p1) .mg_pt_mul(e1, e2) else .mg_pt_mul(e2, e1))
           },
           "/" = if (p1 && !p2) return(.mg_pt_div(e1, e2)))
  }
  .mg_err("TypeError", "unsupported operand type(s) for ", .Generic)
}

#' @export
Math.GeoEntity <- function(x, ...) {
  if (.Generic == "abs" && inherits(x, "GeoPoint")) return(distance(origin(x), x))
  .mg_err("TypeError", "bad operand type for ", .Generic, "(): '", x$gtype, "'")
}

# structural equality (Basic.__eq__), with Polygon and Ellipse overrides
.mg_geo_eq <- function(a, b) {
  if (!inherits(a, "GeoEntity") || !inherits(b, "GeoEntity")) {
    if (inherits(a, "GeoPoint") || inherits(b, "GeoPoint")) return(FALSE)
    return(FALSE)
  }
  if (inherits(a, "GeoPoint") || inherits(b, "GeoPoint")) {
    if (!(inherits(a, "GeoPoint") && inherits(b, "GeoPoint")) || length(a$coords) != length(b$coords)) return(FALSE)
    return(identical(.mg_pkey(a), .mg_pkey(b)))
  }
  if (inherits(a, "GeoEllipse") && !inherits(a, "GeoCircle") || inherits(b, "GeoEllipse") && !inherits(b, "GeoCircle")) {
    return(inherits(a, "GeoEllipse") && inherits(b, "GeoEllipse") && .mg_geo_eq(a$center, b$center) &&
             identical(a$hradius$key, b$hradius$key) && identical(a$vradius$key, b$vradius$key))
  }
  if (inherits(a, "GeoRegularPolygon") && inherits(b, "GeoRegularPolygon")) return(identical(.mg_geo_key(a), .mg_geo_key(b)))
  if (inherits(a, "GeoPolygon") && inherits(b, "GeoPolygon")) {
    if (!identical(a$gtype, b$gtype) || inherits(a, "GeoRegularPolygon") || inherits(b, "GeoRegularPolygon")) return(FALSE)
    ka <- vapply(.mg_poly_canonical(a$vertices), .mg_pkey, "")
    kb <- vapply(.mg_poly_canonical(b$vertices), .mg_pkey, "")
    return(identical(ka, kb))
  }
  if (inherits(a, "GeoEllipse") && inherits(b, "GeoEllipse"))
    return(.mg_geo_eq(a$center, b$center) && identical(a$hradius$key, b$hradius$key) && identical(a$vradius$key, b$vradius$key))
  identical(.mg_geo_key(a), .mg_geo_key(b))
}

# the SymPy args of an entity
.mg_geo_args <- function(e) {
  switch(class(e)[[1L]],
         GeoPoint2D = , GeoPoint3D = , GeoPoint = e$coords,
         GeoCircle = list(e$center, e$radius),
         GeoEllipse = list(e$center, e$hradius, e$vradius),
         GeoRegularPolygon = list(e$center, e$radius, e$n, e$rot),
         GeoTriangle = , GeoPolygon = e$vertices,
         GeoPlane = list(e$p1, new_sym("Tuple", e$normal_vector)),
         GeoCurve = list(new_sym("Tuple", e$functions), new_sym("Tuple", e$limits)),
         GeoParabola = list(e$focus, e$directrix),
         list(e$p1, e$p2))
}

.mg_geo_key <- function(e) {
  paste0(e$gtype, "(", paste(vapply(.mg_geo_args(e), function(a) if (inherits(a, "GeoEntity")) .mg_geo_key(a) else .k(as_sym(a)), ""), collapse = ","), ")")
}

.mg_gstr <- function(a) {
  if (inherits(a, "GeoEntity")) return(format(a))
  if (is_sym(a)) {
    if (a$type == "Tuple") return(paste0("(", paste(vapply(a$args, .mg_gstr, ""), collapse = ", "), if (length(a$args) == 1L) ",", ")"))
    return(.str(a))
  }
  if (is.list(a)) return(paste0("[", paste(vapply(a, .mg_gstr, ""), collapse = ", "), "]"))
  as.character(a)
}

#' @export
format.GeoEntity <- function(x, ...) {
  if (inherits(x, "GeoRegularPolygon"))
    return(sprintf("RegularPolygon(%s, %s, %s, %s)", format(x$center), .str(x$radius), .str(x$n), .str(x$rot)))
  paste0(x$gtype, "(", paste(vapply(.mg_geo_args(x), .mg_gstr, ""), collapse = ", "), ")")
}

#' @export
print.GeoEntity <- function(x, ...) {
  cat(format(x), "\n")
  invisible(x)
}

#' @export
as.character.GeoEntity <- function(x, ...) format(x)

# all points of an entity, in order (atoms(Point))
.mg_geo_points <- function(e) {
  if (inherits(e, "GeoPoint")) return(list(e))
  out <- list()
  for (a in .mg_geo_args(e)) if (inherits(a, "GeoEntity")) out <- c(out, .mg_geo_points(a))
  out[!duplicated(vapply(out, .mg_pkey, ""))]
}

# rebuild an entity, mapping its points with fp and its other expressions
# with fe (as Basic.func(*args) does)
.mg_geo_map <- function(e, fp = NULL, fe = identity) {
  if (inherits(e, "GeoPoint")) {
    if (!is.null(fp)) return(fp(e))
    return(Point(lapply(e$coords, fe)))
  }
  mp <- function(p) .mg_geo_map(p, fp, fe)
  switch(class(e)[[1L]],
         GeoCircle = Circle(mp(e$center), fe(e$radius)),
         GeoEllipse = Ellipse(mp(e$center), fe(e$hradius), fe(e$vradius)),
         GeoRegularPolygon = RegularPolygon(mp(e$center), fe(e$radius), fe(e$n), fe(e$rot)),
         GeoTriangle = do.call(Triangle, lapply(e$vertices, mp)),
         GeoPolygon = do.call(Polygon, lapply(e$vertices, mp)),
         GeoPlane = Plane(mp(e$p1), normal_vector = lapply(e$normal_vector, fe)),
         GeoCurve = Curve(lapply(e$functions, fe), c(list(e$limits[[1L]]), lapply(e$limits[-1L], fe))),
         GeoParabola = Parabola(mp(e$focus), mp(e$directrix)),
         .mg_line_new(e$gtype, mp(e$p1), mp(e$p2)))
}

.mg_geo_subs <- function(e, ...) {
  a <- list(...)
  if (length(a) > 2L) .mg_err("TypeError", "Basic.subs() takes from 2 to 3 positional arguments but ", length(a) + 1L, " were given")
  if (length(a) == 1L && !is.list(a[[1L]]))
    .mg_errf("ValueError", "When a single argument is passed to subs it should be a dictionary of old: new pairs or an iterable of (old, new) tuples.")
  # subs of a point by a point replaces that point
  if (length(a) >= 2L && (inherits(a[[1L]], "GeoEntity") || is.list(a[[1L]]) && !is_sym(a[[1L]]) && length(a[[1L]]) %in% 2:3 && !is.list(a[[1L]][[1L]]))) {
    old <- .mg_as_pt(a[[1L]]); new <- .mg_as_pt(a[[2L]])
    k <- .mg_pkey(old)
    return(.mg_geo_map(e, function(p) if (identical(.mg_pkey(p), k)) new else p))
  }
  if (inherits(e, "GeoCurve")) {
    rv <- .mg_curve_subs(e, a)
    if (!is.null(rv)) return(rv)
  }
  .mg_geo_map(e, NULL, function(z) do.call(subs, c(list(z), a)))
}

.mg_geo_evalf <- function(e, digits = 15) {
  f <- function(z) N(z, digits)
  if (inherits(e, "GeoPoint")) return(Point(lapply(e$coords, f), evaluate = FALSE))
  if (inherits(e, "GeoCircle")) return(.mg_circle_raw(.mg_geo_evalf(e$center, digits), f(e$radius)))
  if (inherits(e, "GeoPlane")) return(.mg_plane_raw(.mg_geo_evalf(e$p1, digits), lapply(e$normal_vector, f)))
  if (inherits(e, "GeoPolygon") && !inherits(e, "GeoRegularPolygon")) {
    v <- lapply(e$vertices, .mg_geo_evalf, digits = digits)
    return(.mg_poly_raw(e$gtype, v))
  }
  if (inherits(e, "GeoLine")) return(.mg_line_raw(e$gtype, .mg_geo_evalf(e$p1, digits), .mg_geo_evalf(e$p2, digits)))
  .mg_geo_map(e, function(p) .mg_geo_evalf(p, digits), f)
}

.mg_geo_free <- function(e) {
  out <- list()
  for (a in .mg_geo_args(e)) out <- c(out, if (inherits(a, "GeoEntity")) .mg_geo_free(a) else free_symbols(as_sym(a)))
  if (inherits(e, "GeoCurve")) out <- Filter(function(s) !identical(s$key, e$parameter$key), out)
  out[!duplicated(vapply(out, .k, ""))]
}

#' @rdname Point
#' @export
ambient_dimension <- function(p) {
  if (inherits(p, "GeoPoint")) return(length(p$coords))
  if (inherits(p, "GeoLine")) return(length(p$p1$coords))
  if (inherits(p, "GeoPlane")) return(3L)
  if (inherits(p, "GeoCurve")) return(length(p$functions))
  2L
}

.mg_pt_list <- function(args) lapply(args, function(a) .mg_as_pt(a))

.mg_uniq_pts <- function(pts) pts[!duplicated(vapply(pts, .mg_pkey, ""))]

#' @rdname Point
#' @export
affine_rank <- function(...) {
  args <- list(...)
  if (!length(args)) return(-1L)
  pts <- .mg_normdim(.mg_pt_list(args))
  o <- pts[[1L]]
  pts <- lapply(pts[-1L], function(p) .mg_pt_sub(p, o))
  if (!length(pts)) return(0L)
  m <- Matrix(lapply(pts, function(p) p$coords))
  isz <- function(x) {
    if (.is_number(x)) return(abs(as.complex(.neval(x, list()))) < 1e-12)
    ask_fact(x, "zero")
  }
  length(.matint_rref(m, isz, FALSE, TRUE)$pivots)
}

# Point.are_coplanar
.mg_point_are_coplanar <- function(...) {
  pts <- list(...)
  if (length(pts) <= 1L) return(TRUE)
  pts <- .mg_normdim(.mg_pt_list(pts))
  if (length(pts[[1L]]$coords) == 2L) return(TRUE)
  pts <- .mg_uniq_pts(pts)
  do.call(affine_rank, pts) <= 2L
}

#' @rdname Point
#' @export
is_collinear <- function(...) {
  pts <- .mg_uniq_pts(.mg_normdim(.mg_pt_list(list(...))))
  do.call(affine_rank, pts) <= 1L
}

#' @rdname Point
#' @export
are_collinear <- function(...) is_collinear(...)

#' @rdname Point
#' @export
is_concyclic <- function(...) {
  pts <- .mg_uniq_pts(.mg_normdim(.mg_pt_list(list(...))))
  if (!(do.call(affine_rank, pts) <= 2L)) return(FALSE)
  o <- pts[[1L]]
  pts <- lapply(pts, function(p) .mg_pt_sub(p, o))
  m <- Matrix(lapply(pts, function(p) c(p$coords, list(.mg_pdot(p, p)))))
  piv <- .rref(m)$pivots
  !((length(o$coords) + 1L) %in% piv)
}

.mg_pdot <- function(a, b) {
  bc <- if (inherits(b, "GeoPoint")) b$coords else lapply(.mg_seq(b), as_sym)
  .add(Map(function(u, v) .mul(list(u, v)), a$coords, bc))
}

#' @rdname Point
#' @export
is_scalar_multiple <- function(p, q) {
  pr <- .mg_normdim(list(.mg_as_pt(p), .mg_as_pt(q)))
  s <- pr[[1L]]; o <- pr[[2L]]
  if (length(s$coords) == 2L) {
    rv <- .expr_equals(.sub(.mul(list(s$coords[[1L]], o$coords[[2L]])), .mul(list(o$coords[[1L]], s$coords[[2L]]))), S_Zero)
    if (is.null(rv)) .mg_err("Undecidable", "\nCannot determine if ", format(s), " is a scalar\nmultiple of                     ", format(o))
  }
  rank_matrix(Matrix(list(s$coords, o$coords))) < 2L
}

.mg_pt_is_zero <- function(p) {
  nz <- vapply(p$coords, function(z) ask_fact(z, "nonzero"), NA)
  if (any(nz %in% TRUE)) return(FALSE)
  if (anyNA(nz)) return(NA)
  TRUE
}

#' @rdname Point
#' @export
midpoint <- function(p, q = NULL) {
  if (is.null(q) && inherits(p, "GeoSegment")) return(midpoint(p$p1, p$p2))
  pr <- .mg_normdim(list(.mg_as_pt(p), .mg_as_pt(q)))
  Point(Map(function(a, b) simplify(.mul(list(.add(list(a, b)), S_Half))), pr[[1L]]$coords, pr[[2L]]$coords))
}

#' @rdname Point
#' @export
origin <- function(p) Point(rep(list(S_Zero), length(p$coords)), evaluate = FALSE)

#' @rdname Point
#' @export
orthogonal_direction <- function(p) {
  d <- length(p$coords)
  z <- rep(list(S_Zero), d)
  if (isTRUE(ask_fact(p$coords[[1L]], "zero"))) return(Point(c(list(S_One), z[-1L])))
  if (isTRUE(ask_fact(p$coords[[2L]], "zero"))) return(Point(c(list(S_Zero, S_One), z[-(1:2)])))
  Point(c(list(.neg(p$coords[[2L]]), p$coords[[1L]]), z[-(1:2)]))
}

#' @rdname Point
#' @export
project <- function(a, b) {
  if (inherits(a, "symmatrix")) {
    b <- Matrix(b)
    return(applyfunc(b, function(z) .mul(list(z, .div(dot(a, b), dot(b, b))))))
  }
  pr <- .mg_normdim(list(.mg_as_pt(a), .mg_as_pt(b)))
  a <- pr[[1L]]; b <- pr[[2L]]
  if (isTRUE(.mg_pt_is_zero(b))) .mg_err("ValueError", "Cannot project to the zero vector.")
  .mg_pt_mul(b, .div(.mg_pdot(a, b), .mg_pdot(b, b)))
}

#' @rdname Point
#' @export
taxicab_distance <- function(p, q) {
  pr <- .mg_normdim(list(.mg_as_pt(p), .mg_as_pt(q)))
  .add(Map(function(a, b) Abs(.sub(a, b)), pr[[1L]]$coords, pr[[2L]]$coords))
}

#' @rdname Point
#' @export
canberra_distance <- function(p, q) {
  pr <- .mg_normdim(list(.mg_as_pt(p), .mg_as_pt(q)))
  if (isTRUE(.mg_pt_is_zero(pr[[1L]])) && isTRUE(.mg_pt_is_zero(pr[[2L]])))
    .mg_err("ValueError", "Cannot project to the zero vector.")
  .add(Map(function(a, b) .div(Abs(.sub(a, b)), .add(list(Abs(a), Abs(b)))), pr[[1L]]$coords, pr[[2L]]$coords))
}

#' @rdname Point
#' @export
unit <- function(p) .mg_pt_div(p, distance(origin(p), p))

.mg_pt_distance <- function(p, q) {
  pr <- .mg_normdim(list(p, .mg_as_pt(q)))
  .sqrt(.add(Map(function(a, b) .pow(.sub(a, b), S_Two), pr[[1L]]$coords, pr[[2L]]$coords)))
}

.mg_pt_rotate <- function(p, angle, pt = NULL) {
  angle <- as_sym(angle)
  c <- cos(angle); s <- sin(angle)
  rv <- p
  if (!is.null(pt)) {
    pt <- .mg_as_pt(pt, 2L)
    rv <- .mg_pt_sub(rv, pt)
  }
  x <- rv$coords[[1L]]; y <- rv$coords[[2L]]
  rv <- Point(.sub(.mul(list(c, x)), .mul(list(s, y))), .add(list(.mul(list(s, x)), .mul(list(c, y)))))
  if (!is.null(pt)) rv <- .mg_pt_add(rv, pt)
  rv
}

.mg_pt_scale <- function(p, x = 1, y = 1, z = 1, pt = NULL) {
  if (!is.null(pt)) {
    pt <- .mg_as_pt(pt, length(p$coords))
    neg <- lapply(pt$coords, .neg)
    r <- do.call(.mg_pt_translate, c(list(p), neg))
    r <- .mg_pt_scale(r, x, y, z)
    return(do.call(.mg_pt_translate, c(list(r), pt$coords)))
  }
  f <- list(as_sym(x), as_sym(y), as_sym(z))
  Point(Map(function(a, b) .mul(list(a, b)), p$coords, f[seq_along(p$coords)]))
}

.mg_pt_translate <- function(p, x = 0, y = 0, z = 0) {
  f <- list(as_sym(x), as_sym(y), as_sym(z))
  Point(Map(function(a, b) .add(list(a, b)), p$coords, f[seq_along(p$coords)]))
}

#' @rdname Point
#' @export
geo_transform <- function(p, matrix) {
  d <- length(p$coords)
  if (!inherits(matrix, "symmatrix") || matrix$nrow != d + 1L || matrix$ncol != d + 1L)
    .mg_err("ValueError", sprintf("matrix must be a %dx%d matrix", d + 1L, d + 1L))
  m <- if (d == 3L) t(matrix) else matrix
  r <- .matmul(.new_mat(c(p$coords, list(S_One)), 1L, d + 1L), m)
  Point(r$data[seq_len(d)])
}

#' @rdname Point
#' @export
direction_ratio <- function(p, q = NULL) {
  if (inherits(p, "GeoLine")) return(direction_ratio(p$p1, p$p2))
  q <- .mg_as_pt(q)
  Map(.sub, q$coords, p$coords)
}

#' @rdname Point
#' @export
direction_cosine <- function(p, q = NULL) {
  if (inherits(p, "GeoLine")) return(direction_cosine(p$p1, p$p2))
  a <- direction_ratio(p, q)
  b <- .sqrt(.add(lapply(a, function(z) .pow(z, S_Two))))
  lapply(a, function(z) .div(z, b))
}

# Point.intersection
.mg_pt_intersection <- function(p, o) {
  if (!inherits(o, "GeoEntity")) o <- Point(o, dim = if (length(p$coords) == 3L) 3L else NULL)
  if (inherits(o, "GeoPoint")) {
    if (.mg_geo_eq(p, o)) return(list(p))
    if (length(p$coords) == 3L) return(list())
    pr <- .mg_normdim(list(p, o))
    if (.mg_geo_eq(pr[[1L]], p) && .mg_geo_eq(pr[[1L]], pr[[2L]])) return(list(p))
    return(list())
  }
  intersection(o, p)
}
