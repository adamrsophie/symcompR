# Differential operators and vector functions, after
# sympy/vector/operators.py, deloperator.py and functions.py.

# The coordinate systems an expression refers to, through its base scalars
# and base vectors (sympy.vector.operators._get_coord_systems).
.get_coord_systems <- function(e) {
  out <- list()
  walk <- function(z) {
    if (!is_sym(z)) return(invisible())
    sys <- if (z$type %in% c("Symbol", "BaseVector") && !is.null(z$extra$system)) z$extra$system
           else if (.is_coordsys(z)) z else NULL
    if (!is.null(sys)) {
      if (is.null(out[[sys$key]])) out[[sys$key]] <<- sys
      return(invisible())
    }
    for (a in z$args) walk(a)
  }
  walk(e)
  unname(out)
}

.split_mul_args_wrt_coordsys <- function(e) {
  keys <- character(); vals <- list()
  for (i in e$args) {
    k <- paste(sort(vapply(.get_coord_systems(i), .k, ""), method = "radix"), collapse = "|")
    j <- match(k, keys)
    if (is.na(j)) { keys <- c(keys, k); vals[[length(vals) + 1L]] <- i }
    else vals[[j]] <- .mul_any(vals[[j]], i)
  }
  vals
}

#' Unevaluated differential operators
#'
#' `Gradient(f)`, `Divergence(v)`, `Curl(v)` and `Laplacian(f)` represent the
#' operators without computing them; [doit()] evaluates.
#'
#' @param expr A scalar or vector field.
#' @return A `sym` expression.
#' @examples
#' R <- CoordSys3D("R")
#' g <- Gradient(R$x * R$y * R$z)
#' g
#' doit(g)
#' @export
Gradient <- function(expr) new_sym("Gradient", list(as_sym(expr)))

#' @rdname Gradient
#' @export
Divergence <- function(expr) new_sym("Divergence", list(as_sym(expr)))

#' @rdname Gradient
#' @export
Curl <- function(expr) new_sym("Curl", list(as_sym(expr)))

#' @rdname Gradient
#' @export
Laplacian <- function(expr) new_sym("Laplacian", list(as_sym(expr)))

register_doit("Gradient", function(e) .vec_gradient(e$args[[1L]], doit = TRUE))
register_doit("Divergence", function(e) .vec_divergence(e$args[[1L]], doit = TRUE))
register_doit("Curl", function(e) .vec_curl(e$args[[1L]], doit = TRUE))
register_doit("Laplacian", function(e) .vec_laplacian(e$args[[1L]]))


.D <- function(f, x) .derivative_node(as_sym(f), list(list(x, 1L)))

.is_vecish <- function(e) .is_vec(e) || e$type %in% c("Cross", "Gradient")

.vec_curl <- function(vect, doit = TRUE) {
  vect <- as_sym(vect)
  cs <- .get_coord_systems(vect)
  if (!length(cs)) return(.VZERO)
  if (length(cs) == 1L) {
    sys <- cs[[1L]]
    bv <- .cs(sys)$base_vectors; bs <- .cs(sys)$base_scalars; h <- .cs(sys)$lame
    i <- bv[[1L]]; j <- bv[[2L]]; k <- bv[[3L]]
    x <- bs[[1L]]; y <- bs[[2L]]; z <- bs[[3L]]
    vx <- .vec_dot_method(vect, i); vy <- .vec_dot_method(vect, j); vz <- .vec_dot_method(vect, k)
    out <- .VZERO
    t1 <- .bd_div(.bd_mul(list(.sub(.D(.mul(list(vz, h[[3L]])), y), .D(.mul(list(vy, h[[2L]])), z)), i), "V"), .mul(list(h[[2L]], h[[3L]])))
    out <- .bd_add2(out, t1)
    t2 <- .bd_div(.bd_mul(list(.sub(.D(.mul(list(vx, h[[1L]])), z), .D(.mul(list(vz, h[[3L]])), x)), j), "V"), .mul(list(h[[1L]], h[[3L]])))
    out <- .bd_add2(out, t2)
    t3 <- .bd_div(.bd_mul(list(.sub(.D(.mul(list(vy, h[[2L]])), x), .D(.mul(list(vx, h[[1L]])), y)), k), "V"), .mul(list(h[[2L]], h[[1L]])))
    out <- .bd_add2(out, t3)
    if (doit) return(doit(out))
    return(out)
  }
  if (vect$type %in% c("Add", "VectorAdd")) {
    args <- tryCatch({
      c1 <- cs[[1L]]
      lapply(vect$args, function(a) express(a, c1, variables = TRUE))
    }, ValueError = function(err) vect$args)
    return(.bd_add(lapply(args, function(a) .vec_curl(a, doit = doit)), "V"))
  }
  if (vect$type %in% c("Mul", "VectorMul")) {
    vector <- Filter(.is_vecish, vect$args)[[1L]]
    scalar <- .mul(Filter(function(a) !.is_vecish(a), vect$args))
    res <- .add_any(doit(Cross(.vec_gradient(scalar), vector)), .mul_any(scalar, .vec_curl(vector, doit = doit)))
    if (doit) return(doit(res))
    return(res)
  }
  if (vect$type %in% c("Cross", "Curl", "Gradient")) return(Curl(vect))
  .vec_stop("ValueError", "Invalid argument for curl")
}

.vec_divergence <- function(vect, doit = TRUE) {
  vect <- as_sym(vect)
  cs <- .get_coord_systems(vect)
  if (!length(cs)) return(S_Zero)
  if (length(cs) == 1L) {
    if (vect$type %in% c("Cross", "Curl", "Gradient")) return(Divergence(vect))
    sys <- cs[[1L]]
    bv <- .cs(sys)$base_vectors; bs <- .cs(sys)$base_scalars; h <- .cs(sys)$lame
    hh <- .mul(list(h[[1L]], h[[2L]], h[[3L]]))
    vx <- .div(.diff_conditional(.vec_dot_method(vect, bv[[1L]]), bs[[1L]], h[[2L]], h[[3L]]), hh)
    vy <- .div(.diff_conditional(.vec_dot_method(vect, bv[[2L]]), bs[[2L]], h[[3L]], h[[1L]]), hh)
    vz <- .div(.diff_conditional(.vec_dot_method(vect, bv[[3L]]), bs[[3L]], h[[1L]], h[[2L]]), hh)
    res <- .add(list(.add(list(vx, vy)), vz))
    if (doit) return(doit(res))
    return(res)
  }
  if (vect$type %in% c("Add", "VectorAdd")) return(.add(lapply(vect$args, function(a) .vec_divergence(a, doit = doit))))
  if (vect$type %in% c("Mul", "VectorMul")) {
    vector <- Filter(.is_vecish, vect$args)[[1L]]
    scalar <- .mul(Filter(function(a) !.is_vecish(a), vect$args))
    res <- .add(list(Dot(vector, .vec_gradient(scalar)), .mul(list(scalar, .vec_divergence(vector, doit = doit)))))
    if (doit) return(doit(res))
    return(res)
  }
  if (vect$type %in% c("Cross", "Curl", "Gradient")) return(Divergence(vect))
  .vec_stop("ValueError", "Invalid argument for divergence")
}

.christoffel_symbol_2nd_kind <- function(h, q, i, j, k) {
  if (i == j && i == k) return(.div(.D(h[[i]], q[[i]]), h[[i]]))
  if (i == k && i != j) return(.div(.D(h[[i]], q[[j]]), h[[i]]))
  if (j == k && i != j) return(.div(.D(h[[j]], q[[i]]), h[[j]]))
  if (i == j && i != k) return(.mul(list(S_NegativeOne, .D(h[[i]], q[[k]]), h[[i]], .pow(h[[k]], .int(-2)))))
  S_Zero
}

.vec_gradient <- function(field, doit = TRUE) {
  field <- as_sym(field)
  cs <- .get_coord_systems(field)
  if (!length(cs)) return(.VZERO)
  if (!.is_bd(field)) {
    if (length(cs) > 1L) {
      if (field$type == "Add") return(.bd_add(lapply(field$args, .vec_gradient), "V"))
      if (field$type == "Mul") {
        s <- .split_mul_args_wrt_coordsys(field)
        return(.bd_add(lapply(s, function(i) .mul_any(.div(field, i), .vec_gradient(i))), "V"))
      }
      return(Gradient(field))
    }
    sys <- cs[[1L]]
    h <- .cs(sys)$lame; bv <- .cs(sys)$base_vectors; bs <- .cs(sys)$base_scalars
    terms <- lapply(1:3, function(n) .bd_mul(list(.div(.D(field, bs[[n]]), h[[n]]), bv[[n]]), "V"))
    out <- .bd_add2(.bd_add2(terms[[1L]], terms[[2L]]), terms[[3L]])
    if (doit) return(doit(out))
    return(out)
  }
  if (.is_vec(field)) {
    if (length(cs) > 1L) {
      if (field$type == "VectorAdd") return(.bd_add(lapply(field$args, .vec_gradient), "D"))
      if (field$type == "VectorMul") {
        s <- .split_mul_args_wrt_coordsys(field)
        return(.bd_add(lapply(s, function(i) .mul_any(.bd_div(field, i), .vec_gradient(i))), "D"))
      }
      return(Gradient(field))
    }
    sys <- cs[[1L]]
    e <- .cs(sys)$base_vectors; q <- .cs(sys)$base_scalars; h <- .cs(sys)$lame
    A <- lapply(e, function(u) .vec_dot_method(field, u))
    s <- .DZERO
    for (i in 1:3) for (j in 1:3) {
      component <- .sub(.div(.D(A[[j]], q[[i]]), h[[i]]),
                        .mul(list(.div(A[[j]], .mul(list(h[[i]], h[[j]]))), .D(h[[j]], q[[i]]))))
      for (k in 1:3) {
        G <- .christoffel_symbol_2nd_kind(h, q, i, k, j)
        component <- .add(list(component, .mul(list(.div(h[[j]], .mul(list(h[[i]], h[[k]]))), A[[k]], G))))
      }
      s <- .bd_add2(s, .bd_mul(list(component, .vec_outer(e[[i]], e[[j]])), "D"))
    }
    if (doit) return(doit(s))
    return(s)
  }
  Gradient(field)
}

.diff_conditional <- function(expr, base_scalar, coeff_1, coeff_2) {
  new_expr <- express(expr, base_scalar$extra$system, variables = TRUE)
  arg <- .mul(list(coeff_1, coeff_2, new_expr))
  if (is_zero_num(arg)) S_Zero else .D(arg, base_scalar)
}

.vec_laplacian <- function(expr) {
  expr <- as_sym(expr)
  if (.is_vec(expr))
    return(doit(.bd_add2(.vec_gradient(.vec_divergence(expr)), .vec_curl(.vec_curl(expr)), sub = TRUE)))
  doit(.vec_divergence(.vec_gradient(expr, doit = FALSE), doit = FALSE))
}

# ---------------------------------------------------------------------------
# The del operator

#' The del operator
#'
#' `Del()` returns the vector differential operator. Called on a scalar
#' field it gives the gradient; `delop & v` gives the divergence,
#' `cross(delop, v)` the curl and `v & delop` a function computing the
#' directional derivative along `v`. Results stay unevaluated unless
#' `doit = TRUE`.
#'
#' @return A function of class `"Del"`.
#' @examples
#' C <- CoordSys3D("C")
#' delop <- Del()
#' delop(C$x * C$y * C$z)
#' doit(delop & (C$x * C$i))
#' @export
Del <- function() {
  structure(function(scalar_field, doit = FALSE) .vec_gradient(scalar_field, doit = doit), class = c("Del", "function"))
}

#' @export
print.Del <- function(x, ...) {
  cat("delop\n")
  invisible(x)
}

#' @export
format.Del <- function(x, ...) "delop"

#' @export
Ops.Del <- function(e1, e2) {
  r <- .del_ops(.Generic, e1, e2)
  if (is.null(r)) .vec_stop("TypeError", paste0("unsupported operand for the del operator: ", .Generic))
  r
}

.del_ops <- function(gen, e1, e2) {
  if (gen == "&") {
    if (inherits(e1, "Del")) return(.vec_divergence(e2, doit = FALSE))
    return(.vec_dot_method(as_sym(e1), e2))
  }
  if (gen == "==") return(inherits(e1, "Del") && inherits(e2, "Del"))
  if (gen == "!=") return(!(inherits(e1, "Del") && inherits(e2, "Del")))
  NULL
}

#' @rdname Del
#' @param delop The del operator.
#' @param field A scalar or vector field.
#' @param doit Evaluate the derivatives.
#' @export
del_gradient <- function(delop, field, doit = FALSE) .vec_gradient(field, doit = doit)

#' @rdname Del
#' @export
del_dot <- function(delop, field, doit = FALSE) .vec_divergence(field, doit = doit)

#' @rdname Del
#' @export
del_cross <- function(delop, field, doit = FALSE) .vec_curl(field, doit = doit)

# ---------------------------------------------------------------------------
# sympy/vector/functions.py

#' Express vectors, dyadics and scalars in a coordinate system
#'
#' `express(expr, system)` re-expresses a vector, dyadic or scalar field in
#' `system`; with `variables = TRUE` the base scalars of other systems are
#' substituted too. `matrix_to_vector()` and `matrix_to_dyadic()` turn
#' component matrices into vectors and dyadics, `orthogonalize()` applies
#' Gram-Schmidt, `directional_derivative()` differentiates a field along a
#' vector, and `scalar_potential_difference()` evaluates the potential
#' difference of a field between two points.
#'
#' @param expr A vector, dyadic or scalar expression.
#' @param system,system2 Coordinate systems ([CoordSys3D()]).
#' @param variables Substitute the coordinate variables.
#' @param matrix A 3 by 1 (or 3 by 3) matrix.
#' @param ... Vectors.
#' @param orthonormal Normalize the result.
#' @param field A scalar or vector field.
#' @param direction_vector A vector.
#' @param coord_sys A coordinate system.
#' @param point1,point2 Points.
#' @param origin For physics reference frames, the origin point used by
#'   [scalar_potential_difference()]; ignored for `CoordSys3D` systems.
#' @return A vector, dyadic, expression or list.
#' @examples
#' N <- CoordSys3D("N")
#' A <- orient_new_axis(N, "A", Symbol("q"), N$k)
#' express(A$i, N)
#' matrix_to_vector(Matrix(list(1, 2, 3)), N)
#' orthogonalize(N$i + 2 * N$j, 2 * N$i + 3 * N$j)
#' @export
express <- function(expr, system, system2 = NULL, variables = FALSE) {
  if (.pv_dispatch(system)) return(.pv_express(expr, system, system2, variables))
  expr <- as_sym(expr)
  if (is_zero_num(expr) || expr$type == "VectorZero") return(expr)
  if (!.is_coordsys(system)) .vec_stop("TypeError", "system should be a CoordSys3D                         instance")
  if (.is_vec(expr)) {
    if (!is.null(system2)) .vec_stop("ValueError", "system2 should not be provided for                                 Vectors")
    sd <- list()
    if (variables) {
      sl <- Filter(function(s) !identical(s$key, system$key), .atom_systems(expr))
      for (f in sl) {
        m <- scalar_map(f, system)
        ks <- attr(m, "keys")
        for (i in seq_along(ks)) sd[[ks[[i]]$key]] <- list(ks[[i]], m[[i]])
      }
      expr <- .subs_dict(expr, sd)
    }
    outvec <- .VZERO
    parts <- .vec_separate(expr)
    for (i in seq_along(parts$keys)) {
      x <- parts$keys[[i]]
      if (!identical(x$key, system$key)) {
        temp <- .matmul(change_of_basis_matrix_from(system, x), to_matrix(parts$values[[i]], x))
        outvec <- .bd_add2(outvec, matrix_to_vector(temp, system))
      } else {
        outvec <- .bd_add2(outvec, parts$values[[i]])
      }
    }
    return(.subs_dict(outvec, sd))
  }
  if (.is_dya(expr)) {
    if (is.null(system2)) system2 <- system
    if (!.is_coordsys(system2)) .vec_stop("TypeError", "system2 should be a CoordSys3D                             instance")
    out <- .DZERO
    for (cm in .bd_components(expr)) {
      o <- .vec_outer(express(cm[[1L]]$args[[1L]], system, variables = variables),
                      express(cm[[1L]]$args[[2L]], system2, variables = variables))
      out <- .bd_add2(out, .mul_any(express(cm[[2L]], system, variables = variables), o))
    }
    return(out)
  }
  if (!is.null(system2)) .vec_stop("ValueError", "system2 should not be provided for                                 Vectors")
  if (variables) {
    sl <- list()
    for (s in .get_scalar_atoms(expr)) {
      sy <- s$extra$system
      if (!identical(sy$key, system$key) && is.null(sl[[sy$key]])) sl[[sy$key]] <- sy
    }
    sd <- list()
    for (f in sl) {
      m <- scalar_map(f, system)
      ks <- attr(m, "keys")
      for (i in seq_along(ks)) sd[[ks[[i]]$key]] <- list(ks[[i]], m[[i]])
    }
    return(.subs_dict(expr, sd))
  }
  expr
}

# the base scalars of an expression
.get_scalar_atoms <- function(e) {
  out <- list()
  walk <- function(z) {
    if (z$type == "Symbol" && isTRUE(z$extra$base_scalar)) { out[[z$key]] <<- z; return(invisible()) }
    for (a in z$args) walk(a)
  }
  walk(e)
  unname(out)
}

# the systems of the BaseScalar and BaseVector atoms
.atom_systems <- function(e) {
  out <- list()
  walk <- function(z) {
    if ((z$type == "Symbol" && isTRUE(z$extra$base_scalar)) || z$type == "BaseVector") {
      s <- z$extra$system
      out[[s$key]] <<- s
      return(invisible())
    }
    for (a in z$args) walk(a)
  }
  walk(e)
  unname(out)
}

# Basic.subs with a dictionary: atoms of equal size are substituted in
# default_sort_key order, one after the other
.subs_dict <- function(e, sd) {
  if (!length(sd)) return(e)
  sd <- unname(sd)
  o <- .key_sort(lapply(sd, function(p) .sort_key(p[[1L]])), seq_along(sd))
  for (i in o) e <- subs(e, sd[[i]][[1L]], sd[[i]][[2L]])
  e
}

#' @rdname express
#' @export
matrix_to_vector <- function(matrix, system) {
  out <- .VZERO
  vects <- .cs(system)$base_vectors
  for (i in seq_along(matrix$data)) out <- .bd_add2(out, .mul_any(matrix$data[[i]], vects[[i]]))
  out
}

#' @rdname express
#' @export
matrix_to_dyadic <- function(matrix, system) {
  if (matrix$nrow != 3L || matrix$ncol != 3L) stop(ShapeError("A 3 x 3 matrix is required."))
  vects <- .cs(system)$base_vectors
  args <- list()
  for (i in 1:3) for (j in 1:3)
    args[[length(args) + 1L]] <- .mul_any(.mget(matrix, i, j), .vec_outer(vects[[i]], vects[[j]]))
  .bd_add(args, "D")
}

#' @rdname express
#' @export
orthogonalize <- function(..., orthonormal = FALSE) {
  vlist <- list(...)
  if (!all(vapply(vlist, .is_vec, TRUE))) .vec_stop("TypeError", "Each element must be of Type Vector")
  ortho <- list()
  for (i in seq_along(vlist)) {
    term <- vlist[[i]]
    for (j in seq_len(i - 1L)) term <- .bd_add2(term, .vec_projection(ortho[[j]], vlist[[i]]), sub = TRUE)
    if (isTRUE(vector_equals(term, .VZERO))) .vec_stop("ValueError", "Vector set not linearly independent")
    ortho[[length(ortho) + 1L]] <- term
  }
  if (orthonormal) ortho <- lapply(ortho, normalize)
  ortho
}

#' @rdname express
#' @export
directional_derivative <- function(field, direction_vector) {
  field <- as_sym(field)
  cs <- .get_coord_systems(field)
  if (!length(cs)) return(if (.is_vec(field)) .VZERO else S_Zero)
  sys <- cs[[1L]]
  field <- express(field, sys, variables = TRUE)
  g <- .vec_gradient(field)
  if (g$type == "Gradient") return(Dot(direction_vector, g))
  .bd_dot_method(direction_vector, g)
}

#' @rdname express
#' @export
scalar_potential_difference <- function(field, coord_sys, point1, point2, origin = NULL) {
  if (.pv_dispatch(coord_sys)) return(.pv_scalar_potential_difference(field, coord_sys, point1, point2, origin))
  if (!.is_coordsys(coord_sys)) .vec_stop("TypeError", "coord_sys must be a CoordSys3D")
  field <- as_sym(field)
  scalar_fn <- if (.is_vec(field)) .vec_scalar_potential(field, coord_sys) else field
  origin <- coord_sys$origin
  position1 <- express(.vpoint_position_wrt(point1, origin), coord_sys, variables = TRUE)
  position2 <- express(.vpoint_position_wrt(point2, origin), coord_sys, variables = TRUE)
  sc <- .cs(coord_sys)$base_scalars
  bv <- .cs(coord_sys)$base_vectors
  sd1 <- lapply(1:3, function(i) list(sc[[i]], .vec_dot_method(bv[[i]], position1)))
  sd2 <- lapply(1:3, function(i) list(sc[[i]], .vec_dot_method(bv[[i]], position2)))
  .sub(.subs_dict(scalar_fn, sd2), .subs_dict(scalar_fn, sd1))
}

.vec_is_conservative <- function(field) {
  if (!.is_vec(field)) .vec_stop("TypeError", "field should be a Vector")
  if (field$type == "VectorZero") return(TRUE)
  identical(.bd_simplify(.vec_curl(field))$key, .VZERO$key)
}

.vec_is_solenoidal <- function(field) {
  if (!.is_vec(field)) .vec_stop("TypeError", "field should be a Vector")
  if (field$type == "VectorZero") return(TRUE)
  is_zero_num(simplify(.vec_divergence(field)))
}

.vec_scalar_potential <- function(field, coord_sys) {
  if (!.vec_is_conservative(field)) .vec_stop("ValueError", "Field is not conservative")
  if (field$type == "VectorZero") return(S_Zero)
  if (!.is_coordsys(coord_sys)) .vec_stop("TypeError", "coord_sys must be a CoordSys3D")
  field <- express(field, coord_sys, variables = TRUE)
  dims <- .cs(coord_sys)$base_vectors
  scalars <- .cs(coord_sys)$base_scalars
  tf <- integrate(.vec_dot_method(field, dims[[1L]]), scalars[[1L]])
  for (i in 2:3) {
    pd <- .diff1(tf, scalars[[i]])
    pd <- .sub(.vec_dot_method(field, dims[[i]]), pd)
    tf <- .add(list(tf, integrate(pd, scalars[[i]])))
  }
  tf
}

# Route the Cartesian helpers of geometry-vector.R: a field holds a
# coordinate system when it is a vector or refers to base scalars.
.is_vector_field <- function(v) is_sym(v) && (.is_bd(v) || length(.get_coord_systems(v)) > 0L)
