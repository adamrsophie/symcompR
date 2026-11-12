# Coordinate systems, base scalars, orienters and points, after
# sympy/vector/coordsysrect.py, scalar.py, orienters.py and point.py.
#
# A CoordSys3D is an atomic sym node (type "CoordSys3D", class
# c("CoordSys3D", "sym")) whose key encodes its name, transformation and
# parent, so equal systems compare equal as in SymPy. Everything else lives in
# an environment in `extra$env`; `C$x`, `C$i` and `C$origin` read it. A base
# scalar is a Symbol named "C.x" that carries its system and index, so the
# rest of the package (diff, integrate, solve, simplify) treats it as the
# symbol it is in SymPy.

.is_coordsys <- function(e) is_sym(e) && identical(e$type, "CoordSys3D")
.is_vpoint <- function(e) is_sym(e) && identical(e$type, "VectorPoint")
.cs <- function(sys) sys$extra$env

#' @export
`$.CoordSys3D` <- function(x, name) {
  if (name %in% c("type", "args", "name", "value", "extra", "key")) return(.subset2(x, name))
  env <- .subset2(x, "extra")$env
  if (name == "origin") return(env$origin)
  i <- match(name, env$variable_names)
  if (!is.na(i)) return(env$base_scalars[[i]])
  i <- match(name, env$vector_names)
  if (!is.na(i)) return(env$base_vectors[[i]])
  .vec_stop("AttributeError", paste0("'CoordSys3D' object has no attribute '", name, "'"))
}

#' @export
`$.VectorPoint` <- function(x, name) {
  if (name %in% c("type", "args", "name", "value", "extra", "key")) return(.subset2(x, name))
  .vec_stop("AttributeError", paste0("'Point' object has no attribute '", name, "'"))
}

register_str_printer("CoordSys3D", function(e) e$name)
register_str_printer("VectorPoint", function(e) e$name)

.check_strings <- function(arg_name, arg) {
  msg <- paste0(arg_name, " must be an iterable of 3 string-types")
  if (length(arg) != 3L) .vec_stop("ValueError", msg)
  for (s in arg) if (!is.character(s)) .vec_stop("TypeError", msg)
}

# Lambda transformations: list(vars, exprs), called by simultaneous xreplace
.vec_lam_call <- function(lam, args) {
  lapply(lam$exprs, function(e) xreplace(e, lapply(seq_along(lam$vars), function(i) list(lam$vars[[i]], args[[i]]))))
}

.vec_lam_key <- function(lam) paste0("Lambda(", paste(vapply(lam$vars, .k, ""), collapse = ","), ";",
                                 paste(vapply(lam$exprs, .k, ""), collapse = ","), ")")

#' Three dimensional coordinate systems
#'
#' `CoordSys3D(name)` creates a Cartesian coordinate system. Its base
#' scalars and base vectors are reached with `$`: `C$x`, `C$y`, `C$z`,
#' `C$i`, `C$j`, `C$k`, and its origin with `C$origin`. A `transformation`
#' of `"spherical"` or `"cylindrical"` (or `"cartesian"`) gives curvilinear
#' systems, whose scalars default to `r, theta, phi` and `r, theta, z`; a
#' list of two lists `list(list(x, y, z), list(f1, f2, f3))` or an R function
#' of three arguments gives a user defined orthogonal system, and
#' `list(rotation_matrix, location)` a rotated and translated child of
#' `parent`.
#'
#' `locate_new()`, `orient_new()`, `orient_new_axis()`, `orient_new_body()`,
#' `orient_new_space()`, `orient_new_quaternion()` and `create_new()` create
#' child systems, as the methods of SymPy's `CoordSys3D` do.
#'
#' @param name The name of the system.
#' @param transformation See Details.
#' @param parent The parent system, or `NULL`.
#' @param location The position of the origin, a vector in `parent`.
#' @param rotation_matrix The rotation matrix with respect to `parent`.
#' @param vector_names,variable_names Three names for the base vectors and
#'   base scalars.
#' @param sys A coordinate system, or (for `locate_new()`) a point.
#' @param position The position of the new origin or point.
#' @param orienters An orienter or a list of orienters.
#' @param angle,angle1,angle2,angle3 Rotation angles.
#' @param axis The rotation axis, a vector.
#' @param rotation_order A rotation order string such as `"123"` or `"XYZ"`.
#' @param q0,q1,q2,q3 Quaternion parameters.
#' @return A `CoordSys3D`, or a point for `locate_new()` on a point.
#' @examples
#' N <- CoordSys3D("N")
#' N$x * N$i + N$y * N$j
#' q <- Symbol("q")
#' A <- orient_new_axis(N, "A", q, N$k)
#' express(A$i, N)
#' S <- create_new(N, "S", transformation = "spherical")
#' lame_coefficients(S)
#' @export
CoordSys3D <- function(name, transformation = NULL, parent = NULL, location = NULL,
                       rotation_matrix = NULL, vector_names = NULL, variable_names = NULL) {
  name <- as.character(name)
  ttype <- "tuple"
  lam <- NULL
  if (!is.null(transformation)) {
    if (!is.null(location) || !is.null(rotation_matrix))
      .vec_stop("ValueError", "specify either `transformation` or `location`/`rotation_matrix`")
    if (is.list(transformation) && !is_sym(transformation) && !inherits(transformation, "symmatrix") && is.null(names(transformation))) {
      if (inherits(transformation[[1L]], "symmatrix")) {
        rotation_matrix <- transformation[[1L]]
        location <- transformation[[2L]]
        transformation <- NULL
      } else {
        ttype <- "lambda"
        lam <- list(vars = lapply(transformation[[1L]], as_sym), exprs = lapply(transformation[[2L]], as_sym))
      }
    } else if (is.function(transformation)) {
      ttype <- "lambda"
      ds <- list(Dummy("x1"), Dummy("x2"), Dummy("x3"))
      lam <- list(vars = ds, exprs = lapply(transformation(ds[[1L]], ds[[2L]], ds[[3L]]), as_sym))
    } else if (is.character(transformation)) {
      ttype <- "str"
    } else if (is_sym(transformation) && transformation$type == "Str") {
      ttype <- "str"
      transformation <- transformation$name
    } else {
      .vec_stop("TypeError", paste0("transformation: wrong type <class '", if (is.list(transformation)) "dict" else class(transformation)[1L], "'>"))
    }
  }
  if (is.null(rotation_matrix)) {
    rotation_matrix <- eye(3)
  } else if (!inherits(rotation_matrix, "symmatrix")) {
    .vec_stop("TypeError", "rotation_matrix should be an ImmutableMatrix instance")
  }
  if (!is.null(parent)) {
    if (!.is_coordsys(parent)) .vec_stop("TypeError", "parent should be a CoordSys3D/None")
    if (is.null(location)) {
      location <- .VZERO
    } else {
      if (!.is_vec(location)) .vec_stop("TypeError", "location should be a Vector")
      for (s in free_symbols(location)) if (isTRUE(s$extra$base_scalar))
        .vec_stop("ValueError", "location should not contain BaseScalars")
    }
    origin <- .vpoint_locate_new(parent$origin, paste0(name, ".origin"), location)
  } else {
    location <- .VZERO
    origin <- .vpoint(paste0(name, ".origin"))
  }
  env <- new.env(parent = emptyenv())
  if (ttype == "tuple") {
    tkey <- paste0("T(", paste(vapply(rotation_matrix$data, .k, ""), collapse = ","), ";", location$key, ")")
    rot <- rotation_matrix
    lambda_transformation <- .compose_rotation_and_translation(t(rot), location, parent)
    lambda_lame <- .lame_lambda("cartesian")
    l <- .vec_projections(location)
    lambda_inverse <- function(x, y, z) .matmul(rot, Matrix(list(.sub(x, l[[1L]]), .sub(y, l[[2L]]), .sub(z, l[[3L]]))))
  } else if (ttype == "str") {
    trname <- transformation
    tkey <- paste0("Str(", trname, ")")
    lambda_transformation <- .transformation_lambdas(trname)
    if (!is.null(parent)) {
      lc <- lame_coefficients(parent)
      if (!all(vapply(lc, is_one_num, TRUE)))
        .vec_stop("ValueError", "Parent for pre-defined coordinate system should be Cartesian.")
    }
    lambda_lame <- .lame_lambda(trname)
    lambda_inverse <- .inv_trans_equations(trname)
  } else {
    if (!.check_orthogonality(lam))
      .vec_stop("ValueError", "The transformation equation does not create orthogonal coordinate system")
    tkey <- .vec_lam_key(lam)
    lambda_transformation <- function(x, y, z) .vec_lam_call(lam, list(x, y, z))
    lambda_lame <- .calculate_lame_coeff(lam)
    lambda_inverse <- NULL
  }
  if (is.null(variable_names)) {
    variable_names <- if (ttype == "lambda") c("x1", "x2", "x3")
      else if (ttype == "str" && transformation == "spherical") c("r", "theta", "phi")
      else if (ttype == "str" && transformation == "cylindrical") c("r", "theta", "z")
      else c("x", "y", "z")
  }
  if (is.null(vector_names)) vector_names <- c("i", "j", "k")
  key <- paste0("CoordSys3D[", nchar(name), ":", name, "|", tkey, "|", if (!is.null(parent)) parent$key else "", "]")
  sys <- new_sym("CoordSys3D", name = name, extra = list(env = env), key = key)
  class(sys) <- c("CoordSys3D", "sym")
  .check_strings("vector_names", vector_names)
  vector_names <- as.character(unlist(vector_names))
  variable_names <- as.character(unlist(variable_names))
  env$name <- name
  env$vector_names <- vector_names
  env$variable_names <- variable_names
  env$base_vectors <- lapply(0:2, function(i) .new_base_vector(i, sys, vector_names[i + 1L]))
  bs <- lapply(0:2, function(i) .new_base_scalar(i, sys, variable_names[i + 1L]))
  env$base_scalars <- bs
  env$ttype <- ttype
  env$transformation <- if (ttype == "str") transformation else if (ttype == "lambda") lam else list(rotation_matrix, location)
  env$lambda_transformation <- lambda_transformation
  env$lame <- lambda_lame(bs[[1L]], bs[[2L]], bs[[3L]])
  env$lambda_inverse <- lambda_inverse
  env$parent <- parent
  env$parent_rotation_matrix <- rotation_matrix
  env$origin <- origin
  if (!is.null(parent)) {
    X <- Matrix(lambda_transformation(bs[[1L]], bs[[2L]], bs[[3L]]))
    h <- env$lame
    J <- jacobian(X, bs)
    is_curvilinear <- any(vapply(J$data, function(z) any(vapply(free_symbols(z), function(s) isTRUE(s$extra$base_scalar), TRUE)), TRUE))
    if (is_curvilinear) {
      d <- J$data
      for (r in 1:3) for (cc in 1:3) d[[(r - 1L) * 3L + cc]] <- .div(.mget(J, r, cc), h[[cc]])
      J <- .new_mat(d, 3L, 3L)
    }
    T_to_parent <- J
  } else {
    is_curvilinear <- FALSE
    T_to_parent <- eye(3)
  }
  env$T_to_parent <- T_to_parent
  env$is_cartesian <- all(vapply(env$lame, is_one_num, TRUE))
  env$is_curvilinear <- is_curvilinear
  sys
}

.new_base_scalar <- function(index, system, vname) {
  nm <- paste0(system$name, ".", vname)
  new_sym("Symbol", name = nm,
          extra = list(facts = .symbol_facts(logical()), given = logical(), base_scalar = TRUE,
                       sympy_class = "BaseScalar", system = system, index = index),
          key = paste0("S", nchar(nm), ":", nm, "{BS", index, "|", system$key, "}"))
}

#' @rdname CoordSys3D
#' @param index An integer 0, 1 or 2.
#' @param system A coordinate system.
#' @export
BaseScalar <- function(index, system) {
  if (!.is_coordsys(system)) .vec_stop("TypeError", "system should be a CoordSys3D")
  if (!(as.numeric(index) %in% 0:2)) .vec_stop("ValueError", "Invalid index specified.")
  .cs(system)$base_scalars[[as.integer(index) + 1L]]
}

.compose_rotation_and_translation <- function(rot, translation, parent) {
  r <- function(x, y, z) .matmul(rot, Matrix(list(x, y, z)))$data
  if (is.null(parent)) return(r)
  dd <- lapply(.cs(parent)$base_vectors, function(i) .vec_dot_method(translation, i))
  function(x, y, z) {
    p <- r(x, y, z)
    list(.add(list(p[[1L]], dd[[1L]])), .add(list(p[[2L]], dd[[2L]])), .add(list(p[[3L]], dd[[3L]])))
  }
}

.transformation_lambdas <- function(nm) {
  if (nm == "cartesian") return(function(x, y, z) list(x, y, z))
  if (nm == "spherical") return(function(r, theta, phi) list(
    .mul(list(r, fcall("sin", list(theta)), fcall("cos", list(phi)))),
    .mul(list(r, fcall("sin", list(theta)), fcall("sin", list(phi)))),
    .mul(list(r, fcall("cos", list(theta))))))
  if (nm == "cylindrical") return(function(r, theta, h) list(
    .mul(list(r, fcall("cos", list(theta)))), .mul(list(r, fcall("sin", list(theta)))), h))
  .vec_stop("ValueError", "Wrong set of parameters.Type of coordinate system is defined")
}

.lame_lambda <- function(nm) {
  if (nm == "cartesian") return(function(x, y, z) list(S_One, S_One, S_One))
  if (nm == "spherical") return(function(r, theta, phi) list(S_One, r, .mul(list(r, fcall("sin", list(theta))))))
  if (nm == "cylindrical") return(function(r, theta, h) list(S_One, r, S_One))
  .vec_stop("ValueError", "Wrong set of parameters. Type of coordinate system is not defined")
}

.inv_trans_equations <- function(nm) {
  if (nm == "cartesian") return(function(x, y, z) list(x, y, z))
  if (nm == "spherical") return(function(x, y, z) {
    rr <- .sqrt(.add(list(.sq(x), .sq(y), .sq(z))))
    list(rr, fcall("acos", list(.div(z, rr))), fcall("atan2", list(y, x)))
  })
  if (nm == "cylindrical") return(function(x, y, z) list(.sqrt(.add(list(.sq(x), .sq(y)))), fcall("atan2", list(y, x)), z))
  .vec_stop("ValueError", "Wrong set of parameters.Type of coordinate system is defined")
}

.check_orthogonality <- function(lam) {
  ds <- list(Dummy("x1"), Dummy("x2"), Dummy("x3"))
  eqs <- .vec_lam_call(lam, ds)
  v <- lapply(ds, function(d) lapply(eqs, function(eq) .diff1(eq, d)))
  if (any(vapply(v, function(i) is_zero_num(simplify(.add(i))), TRUE))) return(FALSE)
  dotp <- function(a, b) simplify(.add(Map(function(p, q) .mul(list(p, q)), a, b)))
  is_zero_num(dotp(v[[1L]], v[[2L]])) && is_zero_num(dotp(v[[2L]], v[[3L]])) && is_zero_num(dotp(v[[3L]], v[[1L]]))
}

.calculate_lame_coeff <- function(lam) {
  function(x1, x2, x3) {
    eqs <- .vec_lam_call(lam, list(x1, x2, x3))
    lapply(list(x1, x2, x3), function(v) .sqrt(.add(lapply(eqs, function(eq) .sq(.diff1(eq, v))))))
  }
}

# ---------------------------------------------------------------------------
# Accessors

#' Properties of coordinate systems
#'
#' `base_vectors()` and `base_scalars()` give the basis of a system (they
#' also accept a diffgeom coordinate system), `lame_coefficients()` its Lame
#' coefficients, `transformation_to_parent()`,
#' `transformation_from_parent()` and `transformation_from_parent_function()`
#' its transformation equations, `rotation_matrix(sys, other)` the direction
#' cosine matrix, `change_of_basis_matrix_from(sys, other)` the change of
#' basis matrix, `scalar_map(sys, other)` the base scalars of `sys` in terms
#' of those of `other`, `position_wrt()` the position of a system or point
#' relative to another, and `express_coordinates()` the coordinates of a
#' point.
#'
#' @param sys A coordinate system (or a point for `position_wrt()` and
#'   `express_coordinates()`).
#' @param other Another coordinate system or point.
#' @param system A coordinate system.
#' @return A list, matrix, function or vector.
#' @examples
#' A <- CoordSys3D("A")
#' B <- orient_new_axis(A, "B", Symbol("q"), A$k)
#' rotation_matrix(A, B)
#' scalar_map(B, A)
#' @export
lame_coefficients <- function(sys) .pytuple(.cs(sys)$lame)

#' @rdname lame_coefficients
#' @export
transformation_to_parent <- function(sys) {
  bs <- .cs(sys)$base_scalars
  .pytuple(.cs(sys)$lambda_transformation(bs[[1L]], bs[[2L]], bs[[3L]]))
}

#' @rdname lame_coefficients
#' @export
transformation_from_parent <- function(sys) {
  env <- .cs(sys)
  if (is.null(env$parent)) .vec_stop("ValueError", "no parent coordinate system, use `transformation_from_parent_function()`")
  pb <- .cs(env$parent)$base_scalars
  r <- env$lambda_inverse(pb[[1L]], pb[[2L]], pb[[3L]])
  if (inherits(r, "symmatrix")) r else .pytuple(r)
}

#' @rdname lame_coefficients
#' @export
transformation_from_parent_function <- function(sys) {
  f <- .cs(sys)$lambda_inverse
  if (is.null(f)) return(NULL)
  function(x, y, z) {
    r <- f(as_sym(x), as_sym(y), as_sym(z))
    if (inherits(r, "symmatrix")) r else .pytuple(r)
  }
}

# CoordSys3D._calculate_inv_trans_equations: solve the transformation for
# the coordinates of the parent.
.calculate_inv_trans_equations <- function(sys) {
  env <- .cs(sys)
  lam <- env$transformation
  x1 <- Dummy("x1", real = TRUE); x2 <- Dummy("x2", real = TRUE); x3 <- Dummy("x3", real = TRUE)
  x <- Dummy("x"); y <- Dummy("y"); z <- Dummy("z")
  eqs <- .vec_lam_call(lam, list(x1, x2, x3))
  sol <- solve(list(.sub(eqs[[1L]], x), .sub(eqs[[2L]], y), .sub(eqs[[3L]], z)), list(x1, x2, x3), dict = TRUE)[[1L]]
  solved <- list(sol[[x1$name]], sol[[x2$name]], sol[[x3$name]])
  env$lambda_inverse <- function(a1, a2, a3) lapply(solved, function(i) subs(i, list(list(x, a1), list(y, a2), list(z, a3))))
  invisible(sys)
}

#' @rdname lame_coefficients
#' @export
rotation_matrix <- function(sys, other) {
  if (!.is_coordsys(other)) .vec_stop("TypeError", paste0(.pystr(other), " is not a CoordSys3D"))
  if (identical(other$key, sys$key)) return(eye(3))
  sp <- .cs(sys)$parent
  op <- .cs(other)$parent
  if (!is.null(sp) && identical(other$key, sp$key)) return(.cs(sys)$parent_rotation_matrix)
  if (!is.null(op) && identical(op$key, sys$key)) return(t(.cs(other)$parent_rotation_matrix))
  p <- .vpath(sys, other)
  result <- eye(3)
  for (i in seq_len(p$index)) result <- .matmul(result, .cs(p$path[[i]])$parent_rotation_matrix)
  if (p$index + 2L <= length(p$path))
    for (i in (p$index + 2L):length(p$path)) result <- .matmul(result, t(.cs(p$path[[i]])$parent_rotation_matrix))
  result
}

#' @rdname lame_coefficients
#' @export
change_of_basis_matrix_from <- function(sys, other) {
  if (!.is_coordsys(other)) .vec_stop("TypeError", paste0(.pystr(other), " is not a CoordSys3D"))
  if (identical(sys$key, other$key)) return(eye(3))
  p <- .vpath(sys, other)
  result <- eye(3)
  mult <- FALSE
  for (i in seq_len(p$index)) {
    e <- .cs(p$path[[i]])
    result <- .matmul(result, if (e$is_cartesian || e$is_curvilinear) t(e$T_to_parent) else inverse(e$T_to_parent))
    mult <- TRUE
  }
  if (p$index + 2L <= length(p$path))
    for (i in (p$index + 2L):length(p$path)) {
      result <- .matmul(result, .cs(p$path[[i]])$T_to_parent)
      mult <- TRUE
    }
  # TR5 rewrites an immutable product; the unmultiplied identity is mutable
  # in SymPy and comes back unchanged
  if (mult) result <- .mat_map(result, TR5)
  result
}

#' @rdname lame_coefficients
#' @export
position_wrt <- function(sys, other) {
  if (.is_coordsys(sys)) sys <- sys$origin
  .vpoint_position_wrt(sys, other)
}

#' @rdname lame_coefficients
#' @export
express_coordinates <- function(sys, system) {
  pv <- .vpoint_position_wrt(sys, system$origin)
  .pytuple(to_matrix(pv, system)$data)
}

.VEC_TRIG_FUNS <- c("sin", "cos", "tan", "cot", "sec", "csc")

#' @rdname lame_coefficients
#' @export
scalar_map <- function(sys, other) {
  p <- .vpath(other, sys)
  systems <- p$path
  root_idx <- p$index
  fb <- .cs(sys)$base_scalars
  current <- list(keys = fb, vals = fb)
  apply_cur <- function(e) subs(e, lapply(seq_along(current$keys), function(i) list(current$keys[[i]], current$vals[[i]])))
  for (s in systems[seq_len(root_idx)]) {
    pt <- transformation_to_parent(s)
    par <- .cs(s)$parent
    vals <- lapply(pt, apply_cur)
    current <- list(keys = .cs(par)$base_scalars, vals = vals)
  }
  if (root_idx + 2L <= length(systems)) for (s in systems[(root_idx + 2L):length(systems)]) {
    env <- .cs(s)
    if (is.null(env$lambda_inverse)) {
      pbs <- .cs(env$parent)$base_scalars
      tp <- transformation_to_parent(s)
      eqs <- Map(.sub, pbs, tp)
      sol <- solve(eqs, env$base_scalars, dict = TRUE)[[1L]]
      xs <- lapply(env$base_scalars, function(k) sol[[k$name]])
    } else {
      xs <- transformation_from_parent(s)
      if (inherits(xs, "symmatrix")) xs <- xs$data
    }
    xs <- lapply(xs, apply_cur)
    current <- list(keys = env$base_scalars, vals = xs)
  }
  post <- function(expr) {
    pattern <- function(t) {
      t$type == "Mul" &&
        any(vapply(t$args, function(a) a$type == "Function" && a$name %in% .VEC_TRIG_FUNS, TRUE)) &&
        any(vapply(t$args, function(a) a$type == "Add" && .has_trig(a), TRUE))
    }
    found <- FALSE
    walk <- function(z) {
      if (found) return(invisible())
      if (pattern(z)) { found <<- TRUE; return(invisible()) }
      for (a in z$args) walk(a)
    }
    walk(expr)
    if (!found) return(expr)
    TR5(factor_terms(expand(expr)))
  }
  vals <- lapply(current$vals, post)
  structure(stats::setNames(vals, vapply(current$keys, function(k) k$name, "")), class = c("scalar_map", "list"),
            keys = current$keys)
}

.has_trig <- function(e) {
  if (e$type == "Function" && e$name %in% .VEC_TRIG_FUNS) return(TRUE)
  any(vapply(e$args, .has_trig, TRUE))
}

#' @export
format.scalar_map <- function(x, ...) format(.pydict(attr(x, "keys"), unname(unclass(x))))

#' @export
print.scalar_map <- function(x, ...) {
  cat(format(x), "\n")
  invisible(x)
}

# ---------------------------------------------------------------------------
# New systems

#' @rdname CoordSys3D
#' @export
locate_new <- function(sys, name, position, vector_names = NULL, variable_names = NULL) {
  if (.is_vpoint(sys)) return(.vpoint_locate_new(sys, name, position))
  env <- .cs(sys)
  if (is.null(variable_names)) variable_names <- env$variable_names
  if (is.null(vector_names)) vector_names <- env$vector_names
  CoordSys3D(name, location = position, vector_names = vector_names, variable_names = variable_names, parent = sys)
}

#' @rdname CoordSys3D
#' @export
orient_new <- function(sys, name, orienters, location = NULL, vector_names = NULL, variable_names = NULL) {
  env <- .cs(sys)
  if (is.null(variable_names)) variable_names <- env$variable_names
  if (is.null(vector_names)) vector_names <- env$vector_names
  if (inherits(orienters, "Orienter")) {
    final_matrix <- .orienter_matrix(orienters, sys)
    final_matrix <- trigsimp(final_matrix)
  } else {
    final_matrix <- eye(3)
    for (o in orienters) final_matrix <- .matmul(final_matrix, .orienter_matrix(o, sys))
  }
  CoordSys3D(name, rotation_matrix = final_matrix, vector_names = vector_names,
             variable_names = variable_names, location = location, parent = sys)
}

#' @rdname CoordSys3D
#' @export
orient_new_axis <- function(sys, name, angle, axis, location = NULL, vector_names = NULL, variable_names = NULL) {
  orient_new(sys, name, AxisOrienter(angle, axis), location = location,
             vector_names = vector_names, variable_names = variable_names)
}

#' @rdname CoordSys3D
#' @export
orient_new_body <- function(sys, name, angle1, angle2, angle3, rotation_order, location = NULL,
                            vector_names = NULL, variable_names = NULL) {
  orient_new(sys, name, BodyOrienter(angle1, angle2, angle3, rotation_order), location = location,
             vector_names = vector_names, variable_names = variable_names)
}

#' @rdname CoordSys3D
#' @export
orient_new_space <- function(sys, name, angle1, angle2, angle3, rotation_order, location = NULL,
                             vector_names = NULL, variable_names = NULL) {
  orient_new(sys, name, SpaceOrienter(angle1, angle2, angle3, rotation_order), location = location,
             vector_names = vector_names, variable_names = variable_names)
}

#' @rdname CoordSys3D
#' @export
orient_new_quaternion <- function(sys, name, q0, q1, q2, q3, location = NULL,
                                  vector_names = NULL, variable_names = NULL) {
  orient_new(sys, name, QuaternionOrienter(q0, q1, q2, q3), location = location,
             vector_names = vector_names, variable_names = variable_names)
}

#' @rdname CoordSys3D
#' @export
create_new <- function(sys, name, transformation, variable_names = NULL, vector_names = NULL) {
  CoordSys3D(name, parent = sys, transformation = transformation,
             variable_names = variable_names, vector_names = vector_names)
}

# ---------------------------------------------------------------------------
# Orienters

#' Orienters
#'
#' `AxisOrienter(angle, axis)`, `BodyOrienter()`, `SpaceOrienter()` and
#' `QuaternionOrienter()` describe rotations for [orient_new()].
#' `orienter_rotation_matrix()` gives the rotation matrix (an axis orienter
#' needs the system the axis is expressed in).
#'
#' @param angle,angle1,angle2,angle3 Angles.
#' @param axis A vector.
#' @param rot_order A rotation order such as `"123"` or `"XYZ"`.
#' @param q0,q1,q2,q3 Quaternion parameters.
#' @param orienter An orienter.
#' @param system A coordinate system (axis orienters only).
#' @return An orienter, or a matrix.
#' @examples
#' A <- CoordSys3D("A")
#' a <- Symbol("a")
#' orienter_rotation_matrix(AxisOrienter(a, A$k), A)
#' orienter_rotation_matrix(BodyOrienter(a, 0, 0, "123"))
#' @export
AxisOrienter <- function(angle, axis) {
  if (!.is_vec(axis)) .vec_stop("TypeError", "axis should be a Vector")
  structure(list(angle = as_sym(angle), axis = axis), class = c("AxisOrienter", "Orienter"))
}

.three_angle <- function(cls, in_order, angle1, angle2, angle3, rot_order) {
  if (is_sym(rot_order) && rot_order$type == "Str") rot_order <- rot_order$name
  approved <- c("123", "231", "312", "132", "213", "321", "121", "131", "212", "232", "313", "323", "")
  original <- rot_order
  ro <- toupper(as.character(rot_order))
  if (nchar(ro) != 3L) .vec_stop("TypeError", "rot_order should be a str of length 3")
  ro <- chartr("XYZ", "123", ro)
  if (!(ro %in% approved)) .vec_stop("TypeError", "Invalid rot_type parameter")
  a <- as.integer(strsplit(ro, "")[[1L]])
  angles <- list(as_sym(angle1), as_sym(angle2), as_sym(angle3))
  m <- if (in_order) .matmul(.matmul(.vrot(a[1L], angles[[1L]]), .vrot(a[2L], angles[[2L]])), .vrot(a[3L], angles[[3L]]))
       else .matmul(.matmul(.vrot(a[3L], angles[[3L]]), .vrot(a[2L], angles[[2L]])), .vrot(a[1L], angles[[1L]]))
  structure(list(angles = angles, rot_order = original, parent_orient = t(m)), class = c(cls, "Orienter"))
}

.vrot <- function(axis, angle) {
  switch(axis, t(rot_axis1(angle)), t(rot_axis2(angle)), t(rot_axis3(angle)))
}

#' @rdname AxisOrienter
#' @export
BodyOrienter <- function(angle1, angle2, angle3, rot_order) .three_angle("BodyOrienter", TRUE, angle1, angle2, angle3, rot_order)

#' @rdname AxisOrienter
#' @export
SpaceOrienter <- function(angle1, angle2, angle3, rot_order) .three_angle("SpaceOrienter", FALSE, angle1, angle2, angle3, rot_order)

#' @rdname AxisOrienter
#' @export
QuaternionOrienter <- function(q0, q1, q2, q3) {
  q0 <- as_sym(q0); q1 <- as_sym(q1); q2 <- as_sym(q2); q3 <- as_sym(q3)
  two <- S_Two
  m <- Matrix(list(
    list(.add(list(.sq(q0), .sq(q1), .neg(.sq(q2)), .neg(.sq(q3)))),
         .mul(list(two, .sub(.mul(list(q1, q2)), .mul(list(q0, q3))))),
         .mul(list(two, .add(list(.mul(list(q0, q2)), .mul(list(q1, q3))))))),
    list(.mul(list(two, .add(list(.mul(list(q1, q2)), .mul(list(q0, q3)))))),
         .add(list(.sq(q0), .neg(.sq(q1)), .sq(q2), .neg(.sq(q3)))),
         .mul(list(two, .sub(.mul(list(q2, q3)), .mul(list(q0, q1)))))),
    list(.mul(list(two, .sub(.mul(list(q1, q3)), .mul(list(q0, q2))))),
         .mul(list(two, .add(list(.mul(list(q0, q1)), .mul(list(q2, q3)))))),
         .add(list(.sq(q0), .neg(.sq(q1)), .neg(.sq(q2)), .sq(q3))))))
  structure(list(q = list(q0, q1, q2, q3), parent_orient = t(m)), class = c("QuaternionOrienter", "Orienter"))
}

.orienter_matrix <- function(o, system = NULL) {
  if (inherits(o, "AxisOrienter")) {
    axis <- normalize(express(o$axis, system))
    ax <- to_matrix(axis, system)$data
    theta <- o$angle
    A <- Matrix(ax)
    AAt <- .matmul(A, t(A))
    K <- Matrix(list(list(S_Zero, .neg(ax[[3L]]), ax[[2L]]), list(ax[[3L]], S_Zero, .neg(ax[[1L]])),
                     list(.neg(ax[[2L]]), ax[[1L]], S_Zero)))
    po <- (eye(3) - AAt) * fcall("cos", list(theta)) + K * fcall("sin", list(theta)) + AAt
    return(t(po))
  }
  o$parent_orient
}

#' @rdname AxisOrienter
#' @export
orienter_rotation_matrix <- function(orienter, system = NULL) .orienter_matrix(orienter, system)

# ---------------------------------------------------------------------------
# Points (sympy/vector/point.py)

.vpoint <- function(name, position = .VZERO, parent_point = NULL) {
  name <- as.character(name)
  if (!.is_vec(position)) .vec_stop("TypeError", paste0("position should be an instance of Vector, not ", .vec_py_class(position)))
  if (!is.null(parent_point) && !.is_vpoint(parent_point))
    .vec_stop("TypeError", paste0("parent_point should be an instance of Point, not ", .vec_py_class(parent_point)))
  env <- new.env(parent = emptyenv())
  env$pos <- position
  env$parent <- parent_point
  p <- new_sym("VectorPoint", name = name, extra = list(env = env),
               key = paste0("VectorPoint[", nchar(name), ":", name, "|", position$key, "|",
                            if (!is.null(parent_point)) parent_point$key else "", "]"))
  class(p) <- c("VectorPoint", "sym")
  p
}

.vpoint_locate_new <- function(p, name, position) .vpoint(name, position, p)

.vparent <- function(o) o$extra$env$parent
.vroot <- function(o) {
  while (!is.null(.vparent(o))) o <- .vparent(o)
  o
}

# sympy.vector.functions._path
.vpath <- function(from_object, to_object) {
  if (!identical(.vroot(from_object)$key, .vroot(to_object)$key))
    .vec_stop("ValueError", paste0("No connecting path found between ", .str(from_object), " and ", .str(to_object)))
  other_path <- list()
  obj <- to_object
  while (!is.null(.vparent(obj))) {
    other_path[[length(other_path) + 1L]] <- obj
    obj <- .vparent(obj)
  }
  other_path[[length(other_path) + 1L]] <- obj
  okeys <- vapply(other_path, .k, "")
  from_path <- list()
  obj <- from_object
  while (!(obj$key %in% okeys)) {
    from_path[[length(from_path) + 1L]] <- obj
    obj <- .vparent(obj)
  }
  index <- length(from_path)
  j <- match(obj$key, okeys)
  list(index = index, path = c(from_path, rev(other_path[seq_len(j)])))
}

.vpoint_position_wrt <- function(self, other) {
  if (!.is_vpoint(other) && !.is_coordsys(other))
    .vec_stop("TypeError", paste0(.pystr(other), "is not a Point or CoordSys3D"))
  if (.is_coordsys(other)) other <- other$origin
  if (identical(other$key, self$key)) return(.VZERO)
  sp <- .vparent(self); op <- .vparent(other)
  if (!is.null(sp) && identical(other$key, sp$key)) return(self$extra$env$pos)
  if (!is.null(op) && identical(op$key, self$key)) return(.bd_mul(list(S_NegativeOne, other$extra$env$pos), "V"))
  p <- .vpath(self, other)
  result <- .VZERO
  for (i in seq_len(p$index)) result <- .bd_add2(result, p$path[[i]]$extra$env$pos)
  if (p$index + 2L <= length(p$path))
    for (i in (p$index + 2L):length(p$path)) result <- .bd_add2(result, p$path[[i]]$extra$env$pos, sub = TRUE)
  result
}
