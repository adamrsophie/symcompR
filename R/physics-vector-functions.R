# Functions of sympy/physics/vector/functions.py and fieldfunctions.py.

.pv_express <- function(expr, frame, frame2 = NULL, variables = FALSE) {
  .pv_check_frame(frame)
  if (.pv_is0(expr)) return(expr)
  if (.pv_is_vec(expr)) {
    if (variables) {
      frame_list <- lapply(.pv_vargs(expr), `[[`, 2L)
      keys <- list(); vals <- list()
      for (f in frame_list) {
        m <- .pv_variable_map(f, frame)
        for (i in seq_along(m$keys)) {
          j <- which(vapply(keys, function(k) identical(k$key, m$keys[[i]]$key), TRUE))
          if (length(j)) vals[[j]] <- m$values[[i]]
          else { keys[[length(keys) + 1L]] <- m$keys[[i]]; vals[[length(vals) + 1L]] <- m$values[[i]] }
        }
      }
      expr <- .pvM_PhysVector_subs(expr, .pydict(keys, vals))
    }
    outvec <- .pv_vec_new(list())
    for (v in .pv_vargs(expr)) {
      if (!identical(v[[2L]], frame)) {
        temp <- .pv_matmul(.pv_dcm(frame, v[[2L]]), v[[1L]])
        if (.pv_state$simp) temp <- .pv_mmap(temp, function(x) trigsimp(x, method = "fu"))
        outvec <- .pv_vadd(outvec, .pv_vec_new(list(list(temp, frame))))
      } else {
        outvec <- .pv_vadd(outvec, .pv_vec_new(list(v)))
      }
    }
    return(outvec)
  }
  if (.pv_is_dyad(expr)) {
    if (is.null(frame2)) frame2 <- frame
    .pv_check_frame(frame2)
    ol <- .pv_dyad0()
    for (v in .pv_dargs(expr)) {
      ol <- .pv_dadd(ol, .pv_dmul(.pv_outer(.pv_express(v[[2L]], frame, variables = variables),
                                            .pv_express(v[[3L]], frame2, variables = variables)),
                                  .pv_express(v[[1L]], frame, variables = variables)))
    }
    return(ol)
  }
  if (variables) {
    ismat <- inherits(expr, "symmatrix")
    if (!ismat) expr <- as_sym(expr)
    frames <- list()
    fs <- if (ismat) .pv_unique_syms(expr) else free_symbols(expr)
    for (x in fs) {
      if (.pv_is_coordsym(x) && !identical(x$extra$frame, frame) && !.pv_in(x$extra$frame, frames))
        frames[[length(frames) + 1L]] <- x$extra$frame
    }
    keys <- list(); vals <- list()
    for (f in frames) {
      m <- .pv_variable_map(f, frame)
      keys <- c(keys, m$keys); vals <- c(vals, m$values)
    }
    if (!length(keys)) return(expr)
    d <- .pydict(keys, vals)
    if (ismat) return(.pv_mmap(expr, function(z) .pv_subs(z, d)))
    return(.pv_subs(expr, d))
  }
  expr
}

.pv_unique_syms <- function(M) {
  out <- list()
  for (z in M$data) for (s in free_symbols(z)) if (!any(vapply(out, function(o) identical(o$key, s$key), TRUE))) out[[length(out) + 1L]] <- s
  out
}

.pv_time_derivative <- function(expr, frame, order = 1) {
  t <- .pv_t()
  .pv_check_frame(frame)
  if (order == 0) return(expr)
  if (order %% 1 != 0 || order < 0) .pv_err("ValueError", "Unsupported value of order entered")
  if (.pv_is_vec(expr)) {
    outlist <- list()
    for (v in .pv_vargs(expr)) {
      if (identical(v[[2L]], frame)) {
        outlist[[length(outlist) + 1L]] <- list(.pv_mdiff(.pv_express(v[[1L]], frame, variables = TRUE), t), frame)
      } else {
        vv <- .pv_vec_new(list(v))
        outlist <- c(outlist, .pv_vargs(.pv_vadd(.pv_time_derivative(vv, v[[2L]]),
                                                 .pv_vcross(.pv_ang_vel_in(v[[2L]], frame), vv))))
      }
    }
    return(.pv_time_derivative(.pv_vec_new(outlist), frame, order - 1))
  }
  if (.pv_is_dyad(expr)) {
    ol <- .pv_dyad0()
    for (v in .pv_dargs(expr)) {
      ol <- .pv_dadd(ol, .pv_dmul(.pv_outer(v[[2L]], v[[3L]]), .pv_diff(v[[1L]], t)))
      ol <- .pv_dadd(ol, .pv_dmul(.pv_outer(.pv_time_derivative(v[[2L]], frame), v[[3L]]), v[[1L]]))
      ol <- .pv_dadd(ol, .pv_dmul(.pv_outer(v[[2L]], .pv_time_derivative(v[[3L]], frame)), v[[1L]]))
    }
    return(.pv_time_derivative(ol, frame, order - 1))
  }
  e <- .pv_express(expr, frame, variables = TRUE)
  if (inherits(e, "symmatrix")) {
    for (i in seq_len(order)) e <- .pv_mdiff(e, t)
    return(e)
  }
  diff(as_sym(e), t, order)
}

#' Time derivative in a reference frame
#'
#' The derivative with respect to time of a vector, dyadic or scalar as seen
#' from `frame`, as SymPy's `time_derivative()`. The time symbol is `t`.
#'
#' @param expr A `PhysVector`, `PhysDyadic` or scalar expression.
#' @param frame A [ReferenceFrame()].
#' @param order Order of the derivative.
#' @return An object of the same kind as `expr`.
#' @examples
#' N <- ReferenceFrame("N")
#' q <- dynamicsymbols("q")
#' A <- N$orientnew("A", "Axis", list(q, N$x))
#' time_derivative(A$y, N)
#' @export
time_derivative <- function(expr, frame, order = 1) .pv_time_derivative(expr, frame, order)

.pv_dot_dispatch <- function(vec1, vec2) {
  if (!(.pv_is_vec(vec1) || .pv_is_dyad(vec1))) .pv_err("TypeError", "Dot product is between two vectors")
  .pv_ops("&", vec1, vec2)
}

.pv_cross_dispatch <- function(vec1, vec2) {
  if (!(.pv_is_vec(vec1) || .pv_is_dyad(vec1))) .pv_err("TypeError", "Cross product is between two vectors")
  .pv_ops("^", vec1, vec2)
}

.pv_outer_dispatch <- function(vec1, vec2) {
  if (!.pv_is_vec(vec1)) .pv_err("TypeError", "Outer product is between two Vectors")
  .pv_outer(vec1, vec2)
}

# true when a call to a function shared with sympy.vector belongs here
.pv_dispatch <- function(...) any(vapply(list(...), function(x) inherits(x, "physics_obj"), TRUE))

#' Kinematic differential equations
#'
#' The equations relating the time derivatives of orientation coordinates
#' to body fixed angular speeds, as SymPy's `kinematic_equations()`.
#'
#' @param speeds A list of three speeds.
#' @param coords A list of three (or four, for quaternions) coordinates.
#' @param rot_type `"body"`, `"space"` or `"quaternion"`.
#' @param rot_order Rotation order such as `"123"` or `"XYZ"`.
#' @return A list of expressions equal to zero.
#' @examples
#' u <- dynamicsymbols("u1 u2 u3")
#' q <- dynamicsymbols("q1 q2 q3")
#' kinematic_equations(u, q, "body", "313")
#' @export
kinematic_equations <- function(speeds, coords, rot_type, rot_order = "") {
  rot_order <- .pv_rot_order(rot_order)
  rot_type <- tolower(rot_type)
  if (!is.list(speeds) || is_sym(speeds)) .pv_err("TypeError", "Need to supply speeds in a list")
  if (length(speeds) != 3L) .pv_err("TypeError", "Need to supply 3 body-fixed speeds")
  if (!is.list(coords) || is_sym(coords)) .pv_err("TypeError", "Need to supply coordinates in a list")
  speeds <- lapply(speeds, as_sym)
  coords <- lapply(coords, as_sym)
  t <- .pv_t()
  if (rot_type %in% c("body", "space")) {
    if (!(rot_order %in% c(.PV_APPROVED_ORDERS, "1", "2", "3"))) .pv_err("ValueError", "Not an acceptable rotation order")
    if (length(coords) != 3L) .pv_err("ValueError", "Need 3 coordinates for body or space")
    w1 <- speeds[[1L]]; w2 <- speeds[[2L]]; w3 <- speeds[[3L]]
    if (is_zero_num(w1) && is_zero_num(w2) && is_zero_num(w3)) return(rep(list(S_Zero), 3L))
    q1 <- coords[[1L]]; q2 <- coords[[2L]]; q3 <- coords[[3L]]
    q1d <- .pv_diff(q1, t); q2d <- .pv_diff(q2, t); q3d <- .pv_diff(q3, t)
    s1 <- sin(q1); s2 <- sin(q2); s3 <- sin(q3)
    c1 <- cos(q1); c2 <- cos(q2); c3 <- cos(q3)
    if (rot_type == "body") {
      r <- switch(rot_order,
        "123" = list(q1d - (w1 * c3 - w2 * s3) / c2, q2d - w1 * s3 - w2 * c3, q3d - (-w1 * c3 + w2 * s3) * s2 / c2 - w3),
        "231" = list(q1d - (w2 * c3 - w3 * s3) / c2, q2d - w2 * s3 - w3 * c3, q3d - w1 - (-w2 * c3 + w3 * s3) * s2 / c2),
        "312" = list(q1d - (-w1 * s3 + w3 * c3) / c2, q2d - w1 * c3 - w3 * s3, q3d - (w1 * s3 - w3 * c3) * s2 / c2 - w2),
        "132" = list(q1d - (w1 * c3 + w3 * s3) / c2, q2d + w1 * s3 - w3 * c3, q3d - (w1 * c3 + w3 * s3) * s2 / c2 - w2),
        "213" = list(q1d - (w1 * s3 + w2 * c3) / c2, q2d - w1 * c3 + w2 * s3, q3d - (w1 * s3 + w2 * c3) * s2 / c2 - w3),
        "321" = list(q1d - (w2 * s3 + w3 * c3) / c2, q2d - w2 * c3 + w3 * s3, q3d - w1 - (w2 * s3 + w3 * c3) * s2 / c2),
        "121" = list(q1d - (w2 * s3 + w3 * c3) / s2, q2d - w2 * c3 + w3 * s3, q3d - w1 + (w2 * s3 + w3 * c3) * c2 / s2),
        "131" = list(q1d - (-w2 * c3 + w3 * s3) / s2, q2d - w2 * s3 - w3 * c3, q3d - w1 - (w2 * c3 - w3 * s3) * c2 / s2),
        "212" = list(q1d - (w1 * s3 - w3 * c3) / s2, q2d - w1 * c3 - w3 * s3, q3d - (-w1 * s3 + w3 * c3) * c2 / s2 - w2),
        "232" = list(q1d - (w1 * c3 + w3 * s3) / s2, q2d + w1 * s3 - w3 * c3, q3d + (w1 * c3 + w3 * s3) * c2 / s2 - w2),
        "313" = list(q1d - (w1 * s3 + w2 * c3) / s2, q2d - w1 * c3 + w2 * s3, q3d + (w1 * s3 + w2 * c3) * c2 / s2 - w3),
        "323" = list(q1d - (-w1 * c3 + w2 * s3) / s2, q2d - w1 * s3 - w2 * c3, q3d - (w1 * c3 - w2 * s3) * c2 / s2 - w3),
        NULL)
      return(r)
    }
    r <- switch(rot_order,
      "123" = list(q1d - w1 - (w2 * s1 + w3 * c1) * s2 / c2, q2d - w2 * c1 + w3 * s1, q3d - (w2 * s1 + w3 * c1) / c2),
      "231" = list(q1d - (w1 * c1 + w3 * s1) * s2 / c2 - w2, q2d + w1 * s1 - w3 * c1, q3d - (w1 * c1 + w3 * s1) / c2),
      "312" = list(q1d - (w1 * s1 + w2 * c1) * s2 / c2 - w3, q2d - w1 * c1 + w2 * s1, q3d - (w1 * s1 + w2 * c1) / c2),
      "132" = list(q1d - w1 - (-w2 * c1 + w3 * s1) * s2 / c2, q2d - w2 * s1 - w3 * c1, q3d - (w2 * c1 - w3 * s1) / c2),
      "213" = list(q1d - (w1 * s1 - w3 * c1) * s2 / c2 - w2, q2d - w1 * c1 - w3 * s1, q3d - (-w1 * s1 + w3 * c1) / c2),
      "321" = list(q1d - (-w1 * c1 + w2 * s1) * s2 / c2 - w3, q2d - w1 * s1 - w2 * c1, q3d - (w1 * c1 - w2 * s1) / c2),
      "121" = list(q1d - w1 + (w2 * s1 + w3 * c1) * c2 / s2, q2d - w2 * c1 + w3 * s1, q3d - (w2 * s1 + w3 * c1) / s2),
      "131" = list(q1d - w1 - (w2 * c1 - w3 * s1) * c2 / s2, q2d - w2 * s1 - w3 * c1, q3d - (-w2 * c1 + w3 * s1) / s2),
      "212" = list(q1d - (-w1 * s1 + w3 * c1) * c2 / s2 - w2, q2d - w1 * c1 - w3 * s1, q3d - (w1 * s1 - w3 * c1) / s2),
      "232" = list(q1d + (w1 * c1 + w3 * s1) * c2 / s2 - w2, q2d + w1 * s1 - w3 * c1, q3d - (w1 * c1 + w3 * s1) / s2),
      "313" = list(q1d + (w1 * s1 + w2 * c1) * c2 / s2 - w3, q2d - w1 * c1 + w2 * s1, q3d - (w1 * s1 + w2 * c1) / s2),
      "323" = list(q1d - (w1 * c1 - w2 * s1) * c2 / s2 - w3, q2d - w1 * s1 - w2 * c1, q3d - (-w1 * c1 + w2 * s1) / s2),
      NULL)
    return(r)
  }
  if (rot_type == "quaternion") {
    if (rot_order != "") .pv_err("ValueError", "Cannot have rotation order for quaternion")
    if (length(coords) != 4L) .pv_err("ValueError", "Need 4 coordinates for quaternion")
    e0 <- coords[[1L]]; e1 <- coords[[2L]]; e2 <- coords[[3L]]; e3 <- coords[[4L]]
    w <- .pv_colv(c(speeds, list(S_Zero)))
    E <- .pv_rows(list(list(e0, -e3, e2, e1), list(e3, e0, -e1, e2), list(-e2, e1, e0, e3), list(-e1, -e2, -e3, e0)))
    edots <- .pv_colv(lapply(list(e1, e2, e3, e0), function(i) .pv_diff(i, t)))
    rhs <- .pv_matmul(.pv_mscale(.pv_T(w), Float(0.5)), .pv_T(E))
    return(.pv_msub(.pv_T(edots), rhs)$data)
  }
  .pv_err("ValueError", "Not an approved rotation type for this function")
}

#' Motion parameters from kinematic conditions
#'
#' Acceleration, velocity and position of a particle in `frame` from one of
#' them and boundary conditions, as SymPy's `get_motion_params()`.
#'
#' @param frame A [ReferenceFrame()].
#' @param ... Any of `acceleration`, `velocity`, `position`, `timevalue`,
#'   `timevalue1` and `timevalue2`.
#' @return A list of three vectors: acceleration, velocity and position.
#' @examples
#' N <- ReferenceFrame("N")
#' get_motion_params(N, position = Symbol("t")^2 * N$x)
#' @export
get_motion_params <- function(frame, ...) {
  kwargs <- list(...)
  process <- function(vectdiff, condition, variable, ordinate, frame) {
    if (!.pv_is0(condition)) condition <- .pv_express(condition, frame, variables = TRUE)
    if (.pv_vec_eq(vectdiff, .pv_vec0())) return(list(S_Zero, S_Zero, condition))
    vectdiff1 <- .pv_express(vectdiff, frame)
    vectdiff2 <- .pv_time_derivative(vectdiff, frame)
    vectdiff0 <- .pv_vec0()
    lims <- list(variable, ordinate, variable)
    for (dim in .pv_frame_iter(frame)) {
      function1 <- .pv_vdot(vectdiff1, dim)
      abscissa <- .pv_subs(.pv_vdot(dim, condition), .pydict(list(variable), list(ordinate)))
      vectdiff0 <- .pv_vadd(vectdiff0, .pv_vmul(dim, .add(list(integrate(function1, lims), abscissa))))
    }
    list(vectdiff2, vectdiff, vectdiff0)
  }
  .pv_check_frame(frame)
  mode <- if (!is.null(kwargs$acceleration)) 2L else if (!is.null(kwargs$velocity)) 1L else 0L
  conds <- c("acceleration", "velocity", "position", "timevalue", "timevalue1", "timevalue2")
  for (i in seq_along(conds)) {
    x <- conds[i]
    if (is.null(kwargs[[x]])) kwargs[[x]] <- if (i <= 3L) .pv_vec0() else S_Zero
    else if (i <= 3L) .pv_check_vector(kwargs[[x]])
    else kwargs[[x]] <- as_sym(kwargs[[x]])
  }
  t <- .pv_t()
  if (mode == 2L) {
    vel <- process(kwargs$acceleration, kwargs$velocity, t, kwargs$timevalue2, frame)[[3L]]
    pos <- process(vel, kwargs$position, t, kwargs$timevalue1, frame)[[3L]]
    return(.pytuple(list(kwargs$acceleration, vel, pos)))
  }
  if (mode == 1L) return(.pytuple(process(kwargs$velocity, kwargs$position, t, kwargs$timevalue1, frame)))
  vel <- .pv_time_derivative(kwargs$position, frame)
  acc <- .pv_time_derivative(vel, frame)
  .pytuple(list(acc, vel, kwargs$position))
}

# sympy.polys.matrices.linsolve._lin_eq2dict
.pv_lin_eq2dict <- function(a, syms) {
  keys <- vapply(syms, .k, "")
  if (a$key %in% keys) return(list(S_Zero, .pydict(list(a), list(S_One))))
  if (a$type == "Add") {
    tk <- list(); tv <- list()
    coeff_list <- list()
    for (ai in .sympy_args(a)) {
      r <- .pv_lin_eq2dict(ai, syms)
      coeff_list[[length(coeff_list) + 1L]] <- r[[1L]]
      for (i in seq_along(r[[2L]]$keys)) {
        k <- r[[2L]]$keys[[i]]
        j <- which(vapply(tk, function(z) identical(z$key, k$key), TRUE))
        if (length(j)) tv[[j]] <- c(tv[[j]], list(r[[2L]]$values[[i]]))
        else { tk[[length(tk) + 1L]] <- k; tv[[length(tv) + 1L]] <- list(r[[2L]]$values[[i]]) }
      }
    }
    return(list(.add(coeff_list), .pydict(tk, lapply(tv, .add))))
  }
  if (a$type == "Mul") {
    terms <- NULL
    terms_coeff <- NULL
    coeff_list <- list()
    for (ai in .sympy_args(a)) {
      r <- .pv_lin_eq2dict(ai, syms)
      if (!length(r[[2L]]$keys)) coeff_list[[length(coeff_list) + 1L]] <- r[[1L]]
      else if (is.null(terms)) { terms <- r[[2L]]; terms_coeff <- r[[1L]] }
      else .pv_err("NonlinearError", .pv_fill(paste0("nonlinear cross-term: ", .str(a))))
    }
    coeff <- .mul(coeff_list)
    if (is.null(terms)) return(list(coeff, .pydict(list(), list())))
    terms <- .pydict(terms$keys, lapply(terms$values, function(c) .mul(list(coeff, c))))
    if (is_zero_num(terms_coeff)) return(list(S_Zero, terms))
    return(list(.mul(list(coeff, terms_coeff)), terms))
  }
  if (!any(vapply(syms, function(s) .has1(a, s), TRUE))) return(list(a, .pydict(list(), list())))
  .pv_err("NonlinearError", paste0("nonlinear term: ", .str(a)))
}

# sympy.solvers.solveset.linear_eq_to_matrix for expressions
.pv_linear_eq_to_matrix <- function(eqs, syms) {
  syms <- lapply(syms, as_sym)
  keys <- vapply(syms, .k, "")
  if (anyDuplicated(keys)) .pv_err("ValueError", "Symbols must be unique")
  A <- rep(list(S_Zero), length(eqs) * length(syms))
  cs <- list()
  for (r in seq_along(eqs)) {
    e <- as_sym(eqs[[r]])
    res <- if (e$type == "Equality") {
      l <- .pv_lin_eq2dict(e$args[[1L]], syms)
      rr <- .pv_lin_eq2dict(e$args[[2L]], syms)
      tk <- l[[2L]]$keys; tv <- l[[2L]]$values
      for (i in seq_along(rr[[2L]]$keys)) {
        k <- rr[[2L]]$keys[[i]]
        j <- which(vapply(tk, function(z) identical(z$key, k$key), TRUE))
        if (length(j)) tv[[j]] <- .sub(tv[[j]], rr[[2L]]$values[[i]])
        else { tk[[length(tk) + 1L]] <- k; tv[[length(tv) + 1L]] <- .neg(rr[[2L]]$values[[i]]) }
      }
      keep <- !vapply(tv, is_zero_num, TRUE)
      list(.sub(l[[1L]], rr[[1L]]), .pydict(tk[keep], tv[keep]))
    } else .pv_lin_eq2dict(e, syms)
    cs[[r]] <- res[[1L]]
    d <- res[[2L]]
    for (i in seq_along(d$keys)) {
      col <- match(d$keys[[i]]$key, keys)
      A[[(r - 1L) * length(syms) + col]] <- d$values[[i]]
    }
  }
  list(.new_mat(A, length(eqs), length(syms)), .new_mat(lapply(cs, .neg), length(eqs), 1L))
}

.pv_partial_velocity <- function(vel_vecs, gen_speeds, frame) {
  if (!is.list(vel_vecs) || inherits(vel_vecs, "physics_obj")) .pv_err("TypeError", "Velocity vectors must be contained in an iterable.")
  if (!(is.list(gen_speeds) || inherits(gen_speeds, "symmatrix")) || is_sym(gen_speeds))
    .pv_err("TypeError", "Generalized speeds must be contained in an iterable")
  gen_speeds <- .pv_seq(gen_speeds)
  out <- list()
  for (vel in vel_vecs) {
    partials <- rep(list(.pv_vec0()), length(gen_speeds))
    for (cr in .pv_vargs(vel)) {
      mat <- .pv_linear_eq_to_matrix(cr[[1L]]$data, gen_speeds)[[1L]]
      dirs <- .pv_frame_iter(cr[[2L]])
      for (i in seq_along(gen_speeds)) for (dim in 1:3) {
        m <- mat$data[[(dim - 1L) * length(gen_speeds) + i]]
        if (!is_zero_num(m)) partials[[i]] <- .pv_vadd(partials[[i]], .pv_vmul(dirs[[dim]], m))
      }
    }
    out[[length(out) + 1L]] <- partials
  }
  out
}

#' Partial velocities
#'
#' The partial velocities of each velocity vector with respect to the
#' generalized speeds, as SymPy's `partial_velocity()`.
#'
#' @param vel_vecs A list of velocity vectors.
#' @param gen_speeds A list of generalized speeds.
#' @param frame The [ReferenceFrame()] the velocities are in.
#' @return A list (one per velocity) of lists of vectors.
#' @examples
#' N <- ReferenceFrame("N")
#' u <- dynamicsymbols("u1 u2")
#' partial_velocity(list(u[[1]] * N$x + u[[2]] * N$y), u, N)
#' @export
partial_velocity <- function(vel_vecs, gen_speeds, frame) .pv_partial_velocity(vel_vecs, gen_speeds, frame)

# symbols()' name expansion: ranges such as `q0:3`, `a:c` and `x1:3_J`,
# several ranges in one name giving their product.
.pv_expand_names <- function(names) {
  parts <- strsplit(trimws(names), "[,[:space:]]+")[[1L]]
  parts <- parts[nzchar(parts)]
  out <- character()
  rx <- "[0-9]*:[0-9]+|[a-zA-Z]?:[a-zA-Z]"
  for (p in parts) {
    if (!grepl(":", p, fixed = TRUE)) {
      out <- c(out, p)
      next
    }
    m <- gregexpr(rx, p)[[1L]]
    ranges <- regmatches(p, list(m))[[1L]]
    lits <- regmatches(p, list(m), invert = TRUE)[[1L]]
    seqs <- list()
    for (i in seq_along(lits)) {
      if (nzchar(lits[i])) seqs[[length(seqs) + 1L]] <- lits[i]
      if (i <= length(ranges)) {
        r <- ranges[i]
        if (endsWith(r, ":")) .pv_err("ValueError", "missing end range")
        ab <- strsplit(r, ":", fixed = TRUE)[[1L]]
        a <- ab[1L]; b <- ab[2L]
        if (grepl("[0-9]$", b)) {
          ai <- if (nzchar(a)) as.integer(a) else 0L
          bi <- as.integer(b)
          v <- if (bi > ai) as.character(seq.int(ai, bi - 1L)) else character()
        } else {
          letters_all <- c(letters, LETTERS)
          if (!nzchar(a)) a <- "a"
          ia <- match(a, letters_all); ib <- match(b, letters_all)
          v <- if (ib >= ia) letters_all[ia:ib] else character()
        }
        seqs[[length(seqs) + 1L]] <- v
        if (!length(v)) break
      }
    }
    if (any(lengths(seqs) == 0L)) next
    combos <- ""
    for (sq in seqs) combos <- as.vector(t(outer(combos, sq, paste0)))
    out <- c(out, combos)
  }
  out
}

#' Dynamic symbols
#'
#' Functions of time `t`, as SymPy's `dynamicsymbols()`: `q(t)`, or its
#' `level`-th time derivative. Names follow [sym_symbols()]: `"q1 q2"`,
#' `"q1, q2"` or `"q1:3"` give a list.
#'
#' @param names Names separated by spaces or commas, or a range.
#' @param level Order of the time derivative.
#' @param ... Assumptions such as `positive = TRUE`, given to the functions.
#' @return A `sym` expression or a list of them.
#' @examples
#' q <- dynamicsymbols("q")
#' dynamicsymbols("q", 1)
#' dynamicsymbols("q1:4")
#' @export
dynamicsymbols <- function(names, level = 0, ...) {
  t <- .pv_t()
  nms <- .pv_expand_names(names)
  # assumptions go on the undefined function, as Function(name, **assumptions)
  given <- .normalize_facts(list(...))
  facts <- .symbol_facts(given)
  if (identical(facts[["commutative"]], FALSE)) .nc_state$any <- TRUE
  make <- function(n) {
    # an undefined function of t, as Function(name)(t); a name this package
    # also uses for a function (beta, gamma) is marked so it stays undefined
    e <- if (length(given)) .pv_undef(n, list(t), list(undefined = TRUE, facts = facts, given = given))
         else if (is.null(.FUNS[[n]])) fcall(n, list(t))
         else new_sym("Function", list(t), name = n, extra = list(undefined = TRUE))
    for (i in seq_len(level)) e <- diff(e, t)
    e
  }
  res <- lapply(nms, make)
  if (length(res) == 1L && !endsWith(trimws(names), ",") && !grepl(":", names, fixed = TRUE)) return(res[[1L]])
  res
}

# an applied undefined function that carries assumptions; they are part of
# its key, as the assumptions are part of SymPy's UndefinedFunction
.pv_undef <- function(name, args, extra) {
  key <- if (is.null(extra$facts)) NULL
         else paste0(.compound_key("Function", args, name), .assumption_suffix(extra$facts))
  new_sym("Function", args, name = name, extra = extra, key = key)
}

# --- fieldfunctions.py

.pv_curl <- function(vect, frame) {
  .pv_check_vector(vect)
  if (.pv_is0(vect)) return(.pv_vec0())
  vect <- .pv_express(vect, frame, variables = TRUE)
  fx <- .pv_frame_iter(frame)
  vx <- .pv_vdot(vect, fx[[1L]]); vy <- .pv_vdot(vect, fx[[2L]]); vz <- .pv_vdot(vect, fx[[3L]])
  c1 <- frame[1]; c2 <- frame[2]; c3 <- frame[3]
  out <- .pv_vec0()
  out <- .pv_vadd(out, .pv_vmul(fx[[1L]], .sub(.pv_diff(vz, c2), .pv_diff(vy, c3))))
  out <- .pv_vadd(out, .pv_vmul(fx[[2L]], .sub(.pv_diff(vx, c3), .pv_diff(vz, c1))))
  out <- .pv_vadd(out, .pv_vmul(fx[[3L]], .sub(.pv_diff(vy, c1), .pv_diff(vx, c2))))
  out
}

.pv_divergence <- function(vect, frame) {
  .pv_check_vector(vect)
  if (.pv_is0(vect)) return(S_Zero)
  vect <- .pv_express(vect, frame, variables = TRUE)
  fx <- .pv_frame_iter(frame)
  out <- S_Zero
  for (i in 1:3) out <- .add(list(out, .pv_diff(.pv_vdot(vect, fx[[i]]), frame[i])))
  out
}

.pv_gradient <- function(scalar, frame) {
  .pv_check_frame(frame)
  out <- .pv_vec0()
  scalar <- .pv_express(scalar, frame, variables = TRUE)
  fx <- .pv_frame_iter(frame)
  for (i in 1:3) out <- .pv_vadd(out, .pv_vmul(fx[[i]], .pv_diff(scalar, frame[i])))
  out
}

.pv_is_conservative <- function(field) {
  if (.pv_is0(field)) return(TRUE)
  frame <- .pv_vargs(field)[[1L]][[2L]]
  .pv_vec_eq(.pvM_PhysVector_simplify(.pv_curl(field, frame)), .pv_vec0())
}

.pv_is_solenoidal <- function(field) {
  if (.pv_is0(field)) return(TRUE)
  frame <- .pv_vargs(field)[[1L]][[2L]]
  identical(simplify(.pv_divergence(field, frame))$key, S_Zero$key)
}

.pv_scalar_potential <- function(field, frame) {
  if (!.pv_is_conservative(field)) .pv_err("ValueError", "Field is not conservative")
  if (.pv_is0(field)) return(S_Zero)
  .pv_check_frame(frame)
  field <- .pv_express(field, frame, variables = TRUE)
  dims <- .pv_frame_iter(frame)
  tf <- integrate(.pv_vdot(field, dims[[1L]]), frame[1])
  for (i in 2:3) {
    pd <- .pv_diff(tf, frame[i])
    pd <- .sub(.pv_vdot(field, dims[[i]]), pd)
    tf <- .add(list(tf, integrate(pd, frame[i])))
  }
  tf
}

.pv_scalar_potential_difference <- function(field, frame, point1, point2, origin) {
  .pv_check_frame(frame)
  scalar_fn <- if (.pv_is_vec(field)) .pv_scalar_potential(field, frame) else as_sym(field)
  position1 <- .pv_express(.pv_pos_from(point1, origin), frame, variables = TRUE)
  position2 <- .pv_express(.pv_pos_from(point2, origin), frame, variables = TRUE)
  fx <- .pv_frame_iter(frame)
  k <- list(frame[1], frame[2], frame[3])
  d1 <- .pydict(k, lapply(fx, function(x) .pv_vdot(x, position1)))
  d2 <- .pydict(k, lapply(fx, function(x) .pv_vdot(x, position2)))
  .sub(.pv_subs(scalar_fn, d2), .pv_subs(scalar_fn, d1))
}
