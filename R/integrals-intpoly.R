# Integration of polynomials over polygons and polyhedra, after
# sympy/integrals/intpoly.py (Chin, Lasserre and Sukumar, 2015). Points are
# handled as lists of coordinates; a polyhedron is a list whose first entry
# holds the vertices and whose other entries are the faces, given by
# zero based vertex indices as in SymPy.

.ip_x <- function() Symbol("x")
.ip_y <- function() Symbol("y")
.ip_z <- function() Symbol("z")

# the coordinates of a point; a geometric point has its floats made
# rational, as SymPy's Point does
.ip_pt <- function(p) {
  if (inherits(p, "GeoPoint")) return(lapply(p$coords, .ip_rat))
  if (is_sym(p)) return(list(p))
  lapply(as.list(p), as_sym)
}

.ip_rat <- function(e) {
  e <- as_sym(e)
  fl <- Filter(function(a) a$type == "Float", .atoms_all(e))
  if (!length(fl)) return(e)
  xreplace(e, lapply(fl, function(f) list(f, nsimplify(f, rational = TRUE))))
}

.atoms_all <- function(e) {
  if (!length(e$args)) return(list(e))
  do.call(c, lapply(e$args, .atoms_all))
}

.ip_same <- function(p, q) length(p) == length(q) && all(vapply(seq_along(p), function(i) identical(p[[i]]$key, q[[i]]$key), TRUE))

# the vertices of a polygon after SymPy's Polygon(): consecutive duplicates
# and collinear points are removed
.ip_polygon <- function(poly) {
  vs <- lapply(poly$vertices, function(v) {
    v <- if (inherits(v, "GeoPoint")) .ip_pt(v) else lapply(lapply(as.list(v), as_sym), .ip_rat)
    v
  })
  nodup <- list()
  for (p in vs) {
    if (length(nodup) && .ip_same(p, nodup[[length(nodup)]])) next
    nodup[[length(nodup) + 1L]] <- p
  }
  if (length(nodup) > 1L && .ip_same(nodup[[length(nodup)]], nodup[[1L]])) nodup[[length(nodup)]] <- NULL
  at <- function(i) nodup[[(if (i < 0L) length(nodup) + i else i) + 1L]]
  pop <- function(i) nodup[[(if (i < 0L) length(nodup) + i else i) + 1L]] <<- NULL
  i <- -3L
  while (i < length(nodup) - 3L && length(nodup) > 2L) {
    a <- at(i); b <- at(i + 1L); c3 <- at(i + 2L)
    cr <- .sub(.mul(list(.sub(b[[1L]], a[[1L]]), .sub(c3[[2L]], a[[2L]]))),
               .mul(list(.sub(b[[2L]], a[[2L]]), .sub(c3[[1L]], a[[1L]]))))
    if (.geo_zero(cr)) {
      pop(i + 1L)
      if (.ip_same(a, c3)) pop(i)
    } else i <- i + 1L
  }
  nodup
}

.ip_sides <- function(vs) {
  n <- length(vs)
  lapply(seq_len(n), function(i) list(vs[[i]], vs[[if (i == n) 1L else i + 1L]]))
}

.ip_seg <- function(s) {
  if (inherits(s, "GeoLine")) return(list(.ip_pt(s$p1), .ip_pt(s$p2)))
  lapply(s, .ip_pt)
}

.ip_is_polygon <- function(p) inherits(p, "GeoPolygon")

.ip_truth <- function(r) isTRUE(.is_true(r))

#' Integration over polytopes
#'
#' `polytope_integrate(poly, expr)` integrates a polynomial in `x` and `y`
#' over a polygon, or in `x`, `y` and `z` over a polyhedron, as
#' `sympy.integrals.intpoly.polytope_integrate()`. A polygon is a
#' [Polygon()] (the vertex representation, whose float coordinates are made
#' rational) or a list of hyperplanes `list(a, b)` meaning `a . (x, y) <= b`;
#' a polyhedron is a list whose first entry is the list of vertices and
#' whose other entries are the faces, each a vector of zero based vertex
#' indices, as in SymPy. The vertices are expected in clockwise order unless
#' `clockwise = TRUE` sorts them. With `max_degree`, `expr` is a list of
#' polynomials and the result is a named list of the integrals of those of
#' degree at most `max_degree` (of every monomial up to that degree when
#' `expr` is `NULL`), named by the polynomials.
#'
#' The helpers follow the SymPy module: `intpoly_decompose(expr, separate)`
#' splits a polynomial by degree (into its terms with `separate = TRUE`),
#' `hyperplane_parameters()` gives the facet hyperplanes, `best_origin()`,
#' `point_sort()`, `intpoly_norm()`, `intpoly_intersection()`,
#' `is_vertex()`, `gradient_terms()`, `intpoly_cross_product()` and
#' `intpoly_strip()` are the geometric and algebraic steps, and
#' `main_integrate()`, `main_integrate3d()`, `polygon_integrate()`,
#' `lineseg_integrate()`, `distance_to_side()`, `integration_reduction()`,
#' `integration_reduction_dynamic()`, `left_integral2D()` and
#' `left_integral3D()` carry out the reduction to lower dimensions.
#'
#' @param poly A polygon, a list of hyperplanes or a polyhedron.
#' @param expr A polynomial, or a list of them.
#' @param clockwise Sort the vertices of a polygon clockwise first.
#' @param max_degree The maximum degree of the polynomials.
#' @return An expression or a named list.
#' @examples
#' x <- Symbol("x"); y <- Symbol("y")
#' polytope_integrate(Polygon(Point(0, 0), Point(0, 1), Point(1, 1), Point(1, 0)), x * y)
#' polytope_integrate(list(list(c(-1, 0), 0), list(c(1, 2), 4), list(c(0, -1), 0)), 1)
#' cube <- list(list(c(0, 0, 0), c(0, 0, 5), c(0, 5, 0), c(0, 5, 5), c(5, 0, 0), c(5, 0, 5),
#'                   c(5, 5, 0), c(5, 5, 5)),
#'              c(3, 7, 6, 2), c(1, 5, 7, 3), c(5, 4, 6, 7), c(0, 4, 5, 1), c(2, 0, 1, 3), c(2, 6, 4, 0))
#' polytope_integrate(cube, x^2 + y^2 + x * y + Symbol("z")^2)
#' intpoly_decompose(x^2 + 3 * y + 7)
#' @export
polytope_integrate <- function(poly, expr = NULL, clockwise = FALSE, max_degree = NULL) {
  if (clockwise) {
    if (.ip_is_polygon(poly)) {
      poly <- structure(list(vertices = point_sort(poly)), class = c("GeoPolygon", "GeoEntity"))
    } else stop(.sympy_condition("TypeError", "clockwise=True works for only 2-PolytopeV-representation input"))
  }
  vertices <- NULL
  if (.ip_is_polygon(poly)) {
    hp_params <- hyperplane_parameters(poly)
    facets <- .ip_sides(.ip_polygon(poly))
  } else if (length(poly[[1L]]) == 2L) {
    plen <- length(poly)
    if (length(poly[[1L]][[1L]]) == 2L) {
      hp_params <- lapply(poly, function(h) list(lapply(as.list(h[[1L]]), as_sym), as_sym(h[[2L]])))
      ints <- lapply(seq_len(plen), function(i)
        intpoly_intersection(hp_params[[if (i == 1L) plen else i - 1L]], hp_params[[i]], "plane2D"))
      ints <- lapply(ints, function(p) lapply(p, .ip_rat))
      facets <- .ip_sides(ints)
    } else stop(.sympy_condition("NotImplementedError", "Integration for H-representation 3Dcase not implemented yet."))
  } else {
    vertices <- lapply(poly[[1L]], .ip_pt)
    facets <- poly[-1L]
    hp_params <- hyperplane_parameters(facets, vertices)
    if (is.null(max_degree)) {
      if (is.null(expr)) stop(.sympy_condition("TypeError", "Input expression must be a valid SymPy expression"))
      return(main_integrate3d(expr, facets, vertices, hp_params))
    }
  }
  if (!is.null(max_degree)) {
    if (!is.null(expr)) {
      if (is_sym(expr) || !is.list(expr))
        stop(.sympy_condition("TypeError", sprintf("'%s' object is not iterable", if (is_sym(expr)) expr$type else "int")))
      f_expr <- list()
      for (e in expr) {
        d <- intpoly_decompose(e)
        if ((length(d) == 1L && names(d) == "0") || total_degree(Poly(as_sym(e))) <= max_degree)
          f_expr[[length(f_expr) + 1L]] <- e
      }
      expr <- f_expr
    }
    result_dict <- if (length(hp_params[[1L]][[1L]]) == 3L)
      main_integrate3d(0, facets, vertices, hp_params, max_degree)
    else main_integrate(0, facets, hp_params, max_degree)
    if (is.null(expr)) return(result_dict)
    result <- list()
    for (p in expr) {
      p <- as_sym(p)
      k <- sstr(p)
      if (k %in% names(result)) next
      if (is_zero_num(p)) { result[["0"]] <- S_Zero; next }
      v <- S_Zero
      for (monom in intpoly_decompose(p, separate = TRUE)) {
        cm <- intpoly_strip(nsimplify(monom))
        v <- .add(list(v, .mul(list(result_dict[[sstr(cm[[2L]])]], cm[[1L]]))))
      }
      result[[k]] <- v
    }
    return(result)
  }
  if (is.null(expr)) stop(.sympy_condition("TypeError", "Input expression must be a valid SymPy expression"))
  main_integrate(expr, facets, hp_params)
}

#' @rdname polytope_integrate
#' @param monom A monomial.
#' @export
intpoly_strip <- function(monom) {
  monom <- as_sym(monom)
  if (is_zero_num(monom)) return(list(S_Zero, S_Zero))
  if (.expr_is_number(monom)) return(list(monom, S_One))
  cf <- LC(monom)
  list(cf, .div(monom, cf))
}

.ip_polynomial_integrate <- function(polynomials, facets, hp_params) {
  dims <- list(.ip_x(), .ip_y())
  out <- S_Zero
  for (deg in names(polynomials)) {
    pc <- S_Zero
    for (i in seq_along(hp_params)) {
      hp <- hp_params[[i]]
      vb <- integration_reduction(facets, i - 1L, hp[[1L]], hp[[2L]], polynomials[[deg]], dims, as_sym(deg))
      pc <- .add(list(pc, .mul(list(vb, .div(hp[[2L]], intpoly_norm(hp[[1L]]))))))
    }
    pc <- .div(pc, .add(list(.int(length(dims)), as_sym(deg))))
    out <- .add(list(out, pc))
  }
  out
}

#' @rdname polytope_integrate
#' @param facets The facets: segments of a polygon, or faces of a polyhedron
#'   as zero based vertex indices.
#' @param vertices The vertices of a polyhedron.
#' @param hp_params The hyperplane parameters of the facets.
#' @export
main_integrate3d <- function(expr, facets, vertices, hp_params, max_degree = NULL) {
  vertices <- lapply(vertices, .ip_pt)
  dims <- list(.ip_x(), .ip_y(), .ip_z())
  dl <- length(dims)
  if (!is.null(max_degree) && max_degree != 0) {
    gt <- gradient_terms(max_degree, 3)
    pos <- list()
    for (zi in seq_along(gt)) for (xi in seq_along(gt[[zi]])) for (yi in seq_along(gt[[zi]][[xi]]))
      pos[[length(pos) + 1L]] <- c(zi, xi, yi)
    result <- list()
    for (p in pos) result[[sstr(gt[[p[1L]]][[p[2L]]][[p[3L]]][[1L]])]] <- S_Zero
    for (fc in seq_along(hp_params)) {
      a <- hp_params[[fc]][[1L]]; b <- hp_params[[fc]][[2L]]
      x0 <- vertices[[facets[[fc]][1L] + 1L]]
      for (i in seq_along(pos)) {
        p <- pos[[i]]
        monom <- gt[[p[1L]]][[p[2L]]][[p[3L]]]
        degree <- monom[[2L]] + monom[[3L]] + monom[[4L]]
        vf <- if (is_zero_num(b)) S_Zero else
          integration_reduction_dynamic(facets, fc - 1L, a, b, monom[[1L]], degree, dims, monom[[7L]], monom[[6L]],
                                        monom[[5L]], x0, gt, i - 1L, vertices, hp_params[[fc]])
        gt[[p[1L]]][[p[2L]]][[p[3L]]][[8L]] <- vf
        k <- sstr(monom[[1L]])
        result[[k]] <- .add(list(result[[k]], .div(.mul(list(vf, .div(b, intpoly_norm(a)))), .int(dl + degree))))
      }
    }
    return(result)
  }
  out <- S_Zero
  polys <- intpoly_decompose(expr)
  for (deg in names(polys)) {
    pc <- S_Zero
    for (i in seq_along(facets)) {
      hp <- hp_params[[i]]
      if (is_zero_num(hp[[2L]])) next
      pi_ <- polygon_integrate(facets[[i]], hp, i - 1L, facets, vertices, expr, as_sym(deg))
      pc <- .add(list(pc, .mul(list(pi_, .div(hp[[2L]], intpoly_norm(hp[[1L]]))))))
    }
    pc <- .div(pc, .add(list(.int(dl), as_sym(deg))))
    out <- .add(list(out, pc))
  }
  out
}

#' @rdname polytope_integrate
#' @export
main_integrate <- function(expr, facets, hp_params, max_degree = NULL) {
  facets <- lapply(facets, .ip_seg)
  hp_params <- lapply(hp_params, function(h) list(lapply(as.list(h[[1L]]), as_sym), as_sym(h[[2L]])))
  dims <- list(.ip_x(), .ip_y())
  dl <- length(dims)
  if (!is.null(max_degree) && max_degree != 0) {
    gt <- c(list(list(0, 0, 0, 0)), gradient_terms(max_degree))
    result <- list()
    for (fc in seq_along(hp_params)) {
      a <- hp_params[[fc]][[1L]]; b <- hp_params[[fc]][[2L]]
      x0 <- facets[[fc]][[1L]]
      for (i in seq_along(gt)) {
        monom <- gt[[i]]
        m <- monom[[1L]]
        k <- if (is_sym(m)) sstr(m) else format(m)
        degree <- S_Zero
        if (is_zero_num(b)) vb <- S_Zero
        else {
          degree <- as_sym(monom[[2L]] + monom[[3L]])
          vb <- integration_reduction_dynamic(facets, fc - 1L, a, b, m, degree, dims, monom[[2L]], monom[[3L]],
                                              max_degree, x0, gt, i - 1L)
        }
        gt[[i]][[4L]] <- vb
        add <- .div(.mul(list(as_sym(vb), .div(b, intpoly_norm(a)))), .add(list(.int(dl), degree)))
        result[[k]] <- if (is.null(result[[k]])) add else .add(list(result[[k]], add))
      }
    }
    return(result)
  }
  if (is_sym(expr) || !is.list(expr))
    return(.ip_polynomial_integrate(intpoly_decompose(expr), facets, hp_params))
  out <- list()
  for (e in expr) out[[sstr(as_sym(e))]] <- .ip_polynomial_integrate(intpoly_decompose(e), facets, hp_params)
  out
}

#' @rdname polytope_integrate
#' @param facet A face, as zero based vertex indices.
#' @param hp_param The hyperplane parameters of the face.
#' @param index The zero based index of the facet.
#' @param degree The degree of the polynomial.
#' @export
polygon_integrate <- function(facet, hp_param, index, facets, vertices, expr, degree) {
  vertices <- lapply(vertices, .ip_pt)
  expr <- as_sym(expr); degree <- as_sym(degree)
  if (is_zero_num(expr)) return(S_Zero)
  out <- S_Zero
  x0 <- vertices[[facet[1L] + 1L]]
  n <- length(facet)
  for (i in seq_len(n)) {
    side <- list(vertices[[facet[i] + 1L]], vertices[[facet[if (i == n) 1L else i + 1L] + 1L]])
    out <- .add(list(out, .mul(list(distance_to_side(x0, side, hp_param[[1L]]),
                                     lineseg_integrate(facet, i - 1L, side, expr, degree)))))
  }
  if (!.expr_is_number(expr)) {
    expr <- .add(list(.mul(list(diff(expr, .ip_x()), x0[[1L]])), .mul(list(diff(expr, .ip_y()), x0[[2L]])),
                      .mul(list(diff(expr, .ip_z()), x0[[3L]]))))
    out <- .add(list(out, polygon_integrate(facet, hp_param, index, facets, vertices, expr, .sub(degree, S_One))))
  }
  .div(out, .add(list(degree, S_Two)))
}

#' @rdname polytope_integrate
#' @param point A point.
#' @param line_seg A segment, as a list of two points.
#' @param A The normal of the face.
#' @export
distance_to_side <- function(point, line_seg, A) {
  point <- .ip_pt(point); line_seg <- lapply(line_seg, .ip_pt); A <- lapply(as.list(A), as_sym)
  nA <- intpoly_norm(A)
  rev_normal <- lapply(A, function(i) .div(.neg(i), nA))
  vec <- lapply(1:3, function(i) .sub(line_seg[[2L]][[i]], line_seg[[1L]][[i]]))
  nv <- intpoly_norm(vec)
  vec <- lapply(vec, function(v) .div(v, nv))
  n_side <- intpoly_cross_product(list(0, 0, 0), rev_normal, vec)
  vx0 <- lapply(1:3, function(i) .sub(line_seg[[1L]][[i]], point[[i]]))
  .add(c(list(S_Zero), lapply(1:3, function(i) .mul(list(vx0[[i]], n_side[[i]])))))
}

#' @rdname polytope_integrate
#' @param polygon The face the segment belongs to (unused, as in SymPy).
#' @export
lineseg_integrate <- function(polygon, index, line_seg, expr, degree) {
  line_seg <- lapply(line_seg, .ip_pt)
  expr <- as_sym(expr); degree <- as_sym(degree)
  if (is_zero_num(expr)) return(S_Zero)
  x0 <- line_seg[[1L]]
  dist <- intpoly_norm(lapply(1:3, function(i) .sub(line_seg[[2L]][[i]], line_seg[[1L]][[i]])))
  out <- .mul(list(dist, subs(expr, list(list(.ip_x(), line_seg[[2L]][[1L]]), list(.ip_y(), line_seg[[2L]][[2L]]),
                                       list(.ip_z(), line_seg[[2L]][[3L]])))))
  expr <- .add(list(.mul(list(diff(expr, .ip_x()), x0[[1L]])), .mul(list(diff(expr, .ip_y()), x0[[2L]])),
                    .mul(list(diff(expr, .ip_z()), x0[[3L]]))))
  out <- .add(list(out, lineseg_integrate(polygon, index, line_seg, expr, .sub(degree, S_One))))
  .div(out, .add(list(degree, S_One)))
}

#' @rdname polytope_integrate
#' @param a,b The hyperplane `a . x = b` of the facet.
#' @param dims The symbols of the dimensions.
#' @export
integration_reduction <- function(facets, index, a, b, expr, dims, degree) {
  facets <- lapply(facets, .ip_seg)
  expr <- as_sym(expr); degree <- as_sym(degree)
  if (is_zero_num(expr)) return(expr)
  value <- S_Zero
  x0 <- facets[[index + 1L]][[1L]]
  m <- length(facets)
  gens <- list(.ip_x(), .ip_y())
  inner <- .add(list(.mul(list(diff(expr, gens[[1L]]), x0[[1L]])), .mul(list(diff(expr, gens[[2L]]), x0[[2L]]))))
  if (!is_zero_num(inner))
    value <- .add(list(value, integration_reduction(facets, index, a, b, inner, dims, .sub(degree, S_One))))
  value <- .add(list(value, left_integral2D(m, index, facets, x0, expr, gens)))
  .div(value, .add(list(.int(length(dims)), degree, S_NegativeOne)))
}

#' @rdname polytope_integrate
#' @param m The number of facets.
#' @param x0 The origin of the facet.
#' @param gens The symbols.
#' @export
left_integral2D <- function(m, index, facets, x0, expr, gens) {
  facets <- lapply(facets, .ip_seg)
  x0 <- .ip_pt(x0)
  value <- S_Zero
  for (j in seq_len(m) - 1L) {
    isect <- NULL
    if (j %in% c((index - 1L) %% m, (index + 1L) %% m))
      isect <- intpoly_intersection(facets[[index + 1L]], facets[[j + 1L]], "segment2D")
    if (length(isect)) {
      d0 <- intpoly_norm(lapply(seq_along(isect), function(i) .sub(isect[[i]], x0[[i]])))
      if (is_vertex(isect)) {
        rule <- lapply(seq_along(gens), function(i) list(gens[[i]], isect[[i]]))
        value <- .add(list(value, .mul(list(d0, subs(as_sym(expr), rule)))))
      }
    }
  }
  value
}

#' @rdname polytope_integrate
#' @param x_index,y_index,max_index Indices of the monomial, as in SymPy.
#' @param monomial_values The table of monomials with their values.
#' @param monom_index The zero based index of the monomial in that table.
#' @export
integration_reduction_dynamic <- function(facets, index, a, b, expr, degree, dims, x_index, y_index, max_index,
                                          x0, monomial_values, monom_index, vertices = NULL, hp_param = NULL) {
  value <- S_Zero
  m <- length(facets)
  if (!is_sym(expr) && expr == 0) return(expr)
  expr <- as_sym(expr); degree <- as_sym(degree)
  if (is_zero_num(expr)) return(expr)
  num <- function(v) if (is.null(v)) S_Zero else as_sym(v)
  if (length(dims) == 2L) {
    facets <- lapply(facets, .ip_seg)
    x0 <- .ip_pt(x0)
    if (!.expr_is_number(expr)) {
      mv <- monomial_values[[monom_index + 1L]]
      xd <- mv[[2L]]; yd <- mv[[3L]]
      xi <- if (xd > 0) monom_index - max_index + x_index - 2L else 0L
      yi <- if (yd > 0) monom_index - 1L else 0L
      pyi <- function(i) (if (i < 0L) length(monomial_values) + i else i) + 1L
      xv <- num(monomial_values[[pyi(xi)]][[4L]]); yv <- num(monomial_values[[pyi(yi)]][[4L]])
      value <- .add(list(value, .mul(list(as_sym(xd), xv, x0[[1L]])), .mul(list(as_sym(yd), yv, x0[[2L]]))))
    }
    value <- .add(list(value, left_integral2D(m, index, facets, x0, expr, dims)))
  } else {
    vertices <- lapply(vertices, .ip_pt)
    x0 <- .ip_pt(x0)
    zi <- max_index
    if (!.expr_is_number(expr)) {
      xd <- y_index; yd <- zi - x_index - y_index; zd <- x_index
      g <- function(t, u, w) num(monomial_values[[t + 1L]][[u + 1L]][[w + 1L]][[8L]])
      xv <- if (xd > 0) g(zi - 1L, y_index - 1L, x_index) else S_Zero
      yv <- if (yd > 0) g(zi - 1L, y_index, x_index) else S_Zero
      zv <- if (zd > 0) g(zi - 1L, y_index, x_index - 1L) else S_Zero
      value <- .add(list(value, .mul(list(as_sym(xd), xv, x0[[1L]])), .mul(list(as_sym(yd), yv, x0[[2L]])),
                         .mul(list(as_sym(zd), zv, x0[[3L]]))))
    }
    value <- .add(list(value, left_integral3D(facets, index, expr, vertices, hp_param, degree)))
  }
  .div(value, .add(list(.int(length(dims)), degree, S_NegativeOne)))
}

#' @rdname polytope_integrate
#' @export
left_integral3D <- function(facets, index, expr, vertices, hp_param, degree) {
  vertices <- lapply(vertices, .ip_pt)
  value <- S_Zero
  facet <- facets[[index + 1L]]
  x0 <- vertices[[facet[1L] + 1L]]
  n <- length(facet)
  for (i in seq_len(n)) {
    side <- list(vertices[[facet[i] + 1L]], vertices[[facet[if (i == n) 1L else i + 1L] + 1L]])
    value <- .add(list(value, .mul(list(distance_to_side(x0, side, hp_param[[1L]]),
                                         lineseg_integrate(facet, i - 1L, side, expr, degree)))))
  }
  value
}

#' @rdname polytope_integrate
#' @param binomial_power The maximum degree.
#' @param no_of_gens The number of generators, 2 or 3.
#' @export
gradient_terms <- function(binomial_power = 0, no_of_gens = 2) {
  x <- .ip_x(); y <- .ip_y(); z <- .ip_z()
  binomial_power <- as.integer(binomial_power)
  mono <- function(i, j, k) .mul(list(.pow(x, .int(i)), .pow(y, .int(j)), .pow(z, .int(k))))
  if (no_of_gens == 2) {
    out <- list()
    for (xc in 0:binomial_power) for (yc in 0:(binomial_power - xc))
      out[[length(out) + 1L]] <- list(mono(xc, yc, 0L), xc, yc, 0L)
    return(out)
  }
  lapply(0:binomial_power, function(zc) lapply(0:zc, function(xc) lapply((zc - xc):0, function(yc)
    list(mono(xc, yc, zc - yc - xc), xc, yc, zc - yc - xc, zc, xc, zc - yc - xc, 0L))))
}

#' @rdname polytope_integrate
#' @export
hyperplane_parameters <- function(poly, vertices = NULL) {
  if (.ip_is_polygon(poly)) {
    vs <- .ip_polygon(poly)
    vs <- c(vs, list(vs[[1L]]))
    return(lapply(seq_len(length(vs) - 1L), function(i) {
      v1 <- vs[[i]]; v2 <- vs[[i + 1L]]
      a1 <- .sub(v1[[2L]], v2[[2L]]); a2 <- .sub(v2[[1L]], v1[[1L]])
      b <- .sub(.mul(list(v2[[1L]], v1[[2L]])), .mul(list(v2[[2L]], v1[[1L]])))
      f <- gcd_list(list(a1, a2, b))
      list(list(.div(a1, f), .div(a2, f)), .div(b, f))
    }))
  }
  vertices <- lapply(vertices, .ip_pt)
  lapply(poly, function(face) {
    v <- lapply(face[1:3], function(k) vertices[[k + 1L]])
    normal <- intpoly_cross_product(v[[1L]], v[[2L]], v[[3L]])
    b <- .add(c(list(S_Zero), lapply(1:3, function(j) .mul(list(normal[[j]], v[[1L]][[j]])))))
    f <- gcd_list(normal)
    if (is_zero_num(f)) f <- S_One
    list(lapply(normal, function(j) .div(j, f)), .div(b, f))
  })
}

#' @rdname polytope_integrate
#' @param v1,v2,v3 Three points.
#' @export
intpoly_cross_product <- function(v1, v2, v3) {
  v1 <- .ip_pt(v1); v2 <- .ip_pt(v2); v3 <- .ip_pt(v3)
  v2 <- lapply(1:3, function(j) .sub(v2[[j]], v1[[j]]))
  v3 <- lapply(1:3, function(j) .sub(v3[[j]], v1[[j]]))
  list(.sub(.mul(list(v3[[3L]], v2[[2L]])), .mul(list(v3[[2L]], v2[[3L]]))),
       .sub(.mul(list(v3[[1L]], v2[[3L]])), .mul(list(v3[[3L]], v2[[1L]]))),
       .sub(.mul(list(v3[[2L]], v2[[1L]])), .mul(list(v3[[1L]], v2[[2L]]))))
}

#' @rdname polytope_integrate
#' @param lineseg A segment.
#' @export
best_origin <- function(a, b, lineseg, expr) {
  pts <- .ip_seg(lineseg)
  a1 <- pts[[1L]][[1L]]; b1 <- pts[[1L]][[2L]]
  zero <- function(v) is_zero_num(as_sym(v))
  x_axis_cut <- function() {
    p <- pts[[1L]]; q <- pts[[2L]]
    if (zero(p[[2L]])) return(p)
    if (zero(q[[2L]])) return(q)
    if (.ip_truth(Lt(.div(p[[2L]], q[[2L]]), S_Zero)))
      return(list(.add(list(.div(.mul(list(p[[2L]], .sub(p[[1L]], q[[1L]]))), .sub(q[[2L]], p[[2L]])), p[[1L]])), S_Zero))
    list()
  }
  y_axis_cut <- function() {
    p <- pts[[1L]]; q <- pts[[2L]]
    if (zero(p[[1L]])) return(p)
    if (zero(q[[1L]])) return(q)
    if (.ip_truth(Lt(.div(p[[1L]], q[[1L]]), S_Zero)))
      return(list(S_Zero, .add(list(.div(.mul(list(p[[1L]], .sub(p[[2L]], q[[2L]]))), .sub(q[[1L]], p[[1L]])), p[[2L]]))))
    list()
  }
  a0 <- a[[1L]]; a2 <- a[[2L]]
  if (zero(a0)) return(if (length(y_axis_cut())) list(S_Zero, b / a2) else list(a1, b1))
  if (zero(a2)) return(if (length(x_axis_cut())) list(b / a0, S_Zero) else list(a1, b1))
  expr <- as_sym(expr)
  gens <- list(.ip_x(), .ip_y())
  gk <- vapply(gens, .k, "")
  pw <- list(S_Zero, S_Zero)
  addp <- function(g, e) { i <- match(g$key, gk); if (!is.na(i)) pw[[i]] <<- .add(list(pw[[i]], e)) }
  if (.expr_is_number(expr)) return(list(a1, b1))
  if (expr$type == "Add") {
    for (mo in expr$args) {
      if (mo$type == "Pow") addp(mo$args[[1L]], mo$args[[2L]])
      else for (u in mo$args) {
        if (!length(u$args)) addp(u, S_One)
        else if (length(u$args) == 2L) addp(u$args[[1L]], u$args[[2L]])
      }
    }
  } else if (expr$type == "Mul") {
    for (u in expr$args) {
      if (!length(u$args)) addp(u, S_One)
      else if (length(u$args) == 2L) addp(u$args[[1L]], u$args[[2L]])
    }
  } else if (expr$type == "Pow") {
    i <- match(expr$args[[1L]]$key, gk)
    if (!is.na(i)) pw[[i]] <- expr$args[[2L]]
  } else if (expr$type == "Symbol") addp(expr, S_One)
  if (.ip_truth(Ge(pw[[1L]], pw[[2L]]))) {
    if (length(y_axis_cut())) list(S_Zero, b / a2)
    else if (length(x_axis_cut())) list(b / a0, S_Zero)
    else list(a1, b1)
  } else {
    if (length(x_axis_cut())) list(b / a0, S_Zero)
    else if (length(y_axis_cut())) list(S_Zero, b / a2)
    else list(a1, b1)
  }
}

#' @rdname polytope_integrate
#' @param separate Return the set of terms instead.
#' @export
intpoly_decompose <- function(expr, separate = FALSE) {
  d <- list()
  e <- if (is_sym(expr)) expr else as_sym(expr)
  put <- function(k, v) d[[sstr(as_sym(k))]] <<- v
  if (is_sym(expr) && !.expr_is_number(e)) {
    if (e$type == "Symbol") put(1, e)
    else if (e$type == "Add") {
      syms <- free_symbols(e)
      degs <- lapply(e$args, function(mo) {
        dl <- if (length(syms)) do.call(degree_list, c(list(mo), syms)) else 0
        list(sum(dl), mo)
      })
      if (separate) return(.ip_unique(lapply(degs, `[[`, 2L)))
      for (dg in degs) {
        k <- as.character(dg[[1L]])
        d[[k]] <- if (!is.null(d[[k]])) .add(list(d[[k]], dg[[2L]])) else dg[[2L]]
      }
    } else if (e$type == "Pow") put(e$args[[2L]], e)
    else {
      deg <- S_Zero
      for (t in e$args) {
        if (!length(t$args) && t$type == "Symbol") deg <- .add(list(deg, S_One))
        else if (length(t$args) == 2L) deg <- .add(list(deg, t$args[[2L]]))
      }
      put(deg, e)
    }
  } else d[["0"]] <- expr
  if (separate) return(.ip_unique(unname(d)))
  d
}

.ip_unique <- function(xs) {
  k <- vapply(xs, function(v) .k(as_sym(v)), "")
  xs[!duplicated(k)]
}

# CPython's list.sort for short lists: the leading run, then a binary
# insertion sort, with lt(a, b) as the order
.ip_pysort <- function(xs, lt) {
  n <- length(xs)
  if (n < 2L) return(xs)
  run <- 2L
  if (lt(xs[[2L]], xs[[1L]])) {
    while (run < n && lt(xs[[run + 1L]], xs[[run]])) run <- run + 1L
    xs[seq_len(run)] <- rev(xs[seq_len(run)])
  } else {
    while (run < n && !lt(xs[[run + 1L]], xs[[run]])) run <- run + 1L
  }
  if (run < n) for (i in (run + 1L):n) {
    pivot <- xs[[i]]
    l <- 1L; r <- i
    while (l < r) {
      p <- (l + r) %/% 2L
      if (lt(pivot, xs[[p]])) r <- p else l <- p + 1L
    }
    if (l < i) xs[(l + 1L):i] <- xs[l:(i - 1L)]
    xs[[l]] <- pivot
  }
  xs
}

#' @rdname polytope_integrate
#' @param normal The normal of the plane of three dimensional points.
#' @export
point_sort <- function(poly, normal = NULL, clockwise = TRUE) {
  pts <- if (.ip_is_polygon(poly)) .ip_polygon(poly) else lapply(poly, .ip_pt)
  n <- length(pts)
  if (n < 2L) return(lapply(pts, function(p) Point(p)))
  order <- if (clockwise) 1L else -1L
  dim <- length(pts[[1L]])
  center <- lapply(seq_len(min(dim, 3L)), function(k) .div(.add(lapply(pts, `[[`, k)), .int(n)))
  T_ <- function(r) .ip_truth(r)
  cmp2 <- function(a, b) {
    ax <- .sub(a[[1L]], center[[1L]]); bx <- .sub(b[[1L]], center[[1L]])
    ay <- .sub(a[[2L]], center[[2L]]); by <- .sub(b[[2L]], center[[2L]])
    if (T_(Ge(ax, S_Zero)) && T_(Lt(bx, S_Zero))) return(-order)
    if (T_(Lt(ax, S_Zero)) && T_(Ge(bx, S_Zero))) return(order)
    if (is_zero_num(ax) && is_zero_num(bx)) {
      if (T_(Ge(ay, S_Zero)) || T_(Ge(by, S_Zero))) return(if (T_(Gt(a[[2L]], b[[2L]]))) -order else order)
      return(if (T_(Gt(b[[2L]], a[[2L]]))) -order else order)
    }
    det <- .sub(.mul(list(ax, by)), .mul(list(bx, ay)))
    if (T_(Lt(det, S_Zero))) return(-order)
    if (T_(Gt(det, S_Zero))) return(order)
    first <- .add(list(.mul(list(ax, ax)), .mul(list(ay, ay))))
    second <- .add(list(.mul(list(bx, bx)), .mul(list(by, by))))
    if (T_(Gt(first, second))) -order else order
  }
  cmp3 <- function(a, b) {
    det <- intpoly_cross_product(center, a, b)
    dp <- .add(lapply(1:3, function(i) .mul(list(det[[i]], as_sym(normal[[i]])))))
    if (T_(Lt(dp, S_Zero))) return(-order)
    if (T_(Gt(dp, S_Zero))) return(order)
    NULL
  }
  cmp <- if (dim == 2L) cmp2 else cmp3
  lt <- function(a, b) {
    r <- cmp(a, b)
    if (is.null(r)) stop(.sympy_condition("TypeError", "'<' not supported between instances of 'NoneType' and 'int'"))
    r < 0
  }
  lapply(.ip_pysort(pts, lt), function(p) Point(p))
}

#' @rdname polytope_integrate
#' @export
intpoly_norm <- function(point) {
  if (inherits(point, "GeoPoint")) {
    c0 <- .ip_pt(point)
    if (length(c0) == 2L) return(.pow(.add(list(.pow(c0[[1L]], S_Two), .pow(c0[[2L]], S_Two))), S_Half))
    return(.pow(.add(list(.pow(c0[[1L]], S_Two), .pow(c0[[2L]], S_Two), c0[[3L]])), S_Half))
  }
  cs <- lapply(as.list(point), as_sym)
  .pow(.add(c(list(S_Zero), lapply(cs, function(v) .pow(v, S_Two)))), S_Half)
}

#' @rdname polytope_integrate
#' @param geom_1,geom_2 Two segments or two hyperplanes.
#' @param intersection_type `"segment2D"`, `"segment3D"` or `"plane2D"`.
#' @export
intpoly_intersection <- function(geom_1, geom_2, intersection_type) {
  if (substr(intersection_type, 1L, nchar(intersection_type) - 2L) == "segment") {
    g1 <- .ip_seg(geom_1); g2 <- .ip_seg(geom_2)
    x1 <- g1[[1L]][[1L]]; y1 <- g1[[1L]][[2L]]; x2 <- g1[[2L]][[1L]]; y2 <- g1[[2L]][[2L]]
    x3 <- g2[[1L]][[1L]]; y3 <- g2[[1L]][[2L]]; x4 <- g2[[2L]][[1L]]; y4 <- g2[[2L]][[2L]]
    den <- .sub(.mul(list(.sub(x1, x2), .sub(y3, y4))), .mul(list(.sub(y1, y2), .sub(x3, x4))))
    if (!is_zero_num(den)) {
      t1 <- .sub(.mul(list(x1, y2)), .mul(list(y1, x2)))
      t2 <- .sub(.mul(list(x3, y4)), .mul(list(x4, y3)))
      return(list(.div(.sub(.mul(list(t1, .sub(x3, x4))), .mul(list(t2, .sub(x1, x2)))), den),
                  .div(.sub(.mul(list(t1, .sub(y3, y4))), .mul(list(t2, .sub(y1, y2)))), den)))
    }
  }
  if (substr(intersection_type, 1L, nchar(intersection_type) - 2L) == "plane" && intersection_type == "plane2D") {
    a1 <- lapply(as.list(geom_1[[1L]]), as_sym); a2 <- lapply(as.list(geom_2[[1L]]), as_sym)
    b1 <- as_sym(geom_1[[2L]]); b2 <- as_sym(geom_2[[2L]])
    den <- .sub(.mul(list(a1[[1L]], a2[[2L]])), .mul(list(a2[[1L]], a1[[2L]])))
    if (!is_zero_num(den))
      return(list(.div(.sub(.mul(list(b1, a2[[2L]])), .mul(list(b2, a1[[2L]]))), den),
                  .div(.sub(.mul(list(b2, a1[[1L]])), .mul(list(b1, a2[[1L]]))), den)))
  }
  NULL
}

#' @rdname polytope_integrate
#' @param ent An object.
#' @export
is_vertex <- function(ent) {
  if (inherits(ent, "GeoPoint")) return(TRUE)
  if (is.list(ent) && !is_sym(ent) && !inherits(ent, "GeoEntity")) return(length(ent) %in% c(2L, 3L))
  if (is.atomic(ent) && length(ent) > 1L) return(length(ent) %in% c(2L, 3L))
  FALSE
}
