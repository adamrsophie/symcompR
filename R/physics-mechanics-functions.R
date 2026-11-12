# Functions of classical mechanics, after sympy/physics/mechanics/functions.py.

.pm_check_bodies <- function(bodies) {
  for (e in bodies) if (!inherits(e, c("RigidBody", "Particle")))
    .pv_err("TypeError", "*body must have only Particle or RigidBody")
}

#' Momentum, energy and the Lagrangian of a system of bodies
#'
#' `linear_momentum()`, `angular_momentum()` (about `point`),
#' `kinetic_energy()` and `potential_energy()` sum the contributions of the
#' given particles and rigid bodies, as the SymPy functions of the same
#' names. `Lagrangian()` is the kinetic minus the potential energy and
#' `center_of_mass()` the position of the mass center from `point`.
#'
#' @param frame The [ReferenceFrame()] in which the quantity is measured.
#' @param point A [PhysPoint()].
#' @param ... Particles and rigid bodies.
#' @return A `PhysVector` or a scalar expression.
#' @examples
#' N <- ReferenceFrame("N")
#' P <- PhysPoint("P")
#' P$set_vel(N, 10 * N$x)
#' p <- Particle("p", P, 2)
#' linear_momentum(N, p)
#' kinetic_energy(N, p)
#' @export
linear_momentum <- function(frame, ...) {
  if (!.pv_is_frame(frame)) .pv_err("TypeError", "Please specify a valid ReferenceFrame")
  linear_momentum_sys <- .pv_vec0()
  for (e in list(...)) {
    if (!inherits(e, c("RigidBody", "Particle"))) .pv_err("TypeError", "*body must have only Particle or RigidBody")
    linear_momentum_sys <- linear_momentum_sys + e$linear_momentum(frame)
  }
  linear_momentum_sys
}

#' @rdname linear_momentum
#' @export
angular_momentum <- function(point, frame, ...) {
  if (!.pv_is_frame(frame)) .pv_err("TypeError", "Please enter a valid ReferenceFrame")
  if (!.pv_is_point(point)) .pv_err("TypeError", "Please specify a valid Point")
  angular_momentum_sys <- .pv_vec0()
  for (e in list(...)) {
    if (!inherits(e, c("RigidBody", "Particle"))) .pv_err("TypeError", "*body must have only Particle or RigidBody")
    angular_momentum_sys <- angular_momentum_sys + e$angular_momentum(point, frame)
  }
  angular_momentum_sys
}

#' @rdname linear_momentum
#' @export
kinetic_energy <- function(frame, ...) {
  if (!.pv_is_frame(frame)) .pv_err("TypeError", "Please enter a valid ReferenceFrame")
  ke_sys <- S_Zero
  for (e in list(...)) {
    if (!inherits(e, c("RigidBody", "Particle"))) .pv_err("TypeError", "*body must have only Particle or RigidBody")
    ke_sys <- ke_sys + e$kinetic_energy(frame)
  }
  ke_sys
}

#' @rdname linear_momentum
#' @export
potential_energy <- function(...) {
  pe_sys <- S_Zero
  for (e in list(...)) {
    if (!inherits(e, c("RigidBody", "Particle"))) .pv_err("TypeError", "*body must have only Particle or RigidBody")
    pe_sys <- pe_sys + e[["_potential_energy"]]
  }
  pe_sys
}

#' @rdname linear_momentum
#' @export
center_of_mass <- function(point, ...) {
  bodies <- list(...)
  if (!length(bodies)) .pv_err("TypeError", "No bodies(instances of Particle or Rigidbody) were passed.")
  total_mass <- S_Zero
  vec <- .pv_vec0()
  for (i in bodies) {
    total_mass <- total_mass + i$mass
    masscenter <- i$masscenter
    vec <- vec + i$mass * .pv_pos_from(masscenter, point)
  }
  vec / total_mass
}

#' @rdname linear_momentum
#' @export
Lagrangian <- function(frame, ...) {
  if (!.pv_is_frame(frame)) .pv_err("TypeError", "Please supply a valid ReferenceFrame")
  .pm_check_bodies(list(...))
  kinetic_energy(frame, ...) - potential_energy(...)
}

# ---------------------------------------------------------------------------
# find_dynamicsymbols

.pm_is_dynsym <- function(e, t = .pv_t()) {
  if (!(e$type == "Derivative" || .is_undefined_application(e))) return(FALSE)
  fs <- free_symbols(e)
  length(fs) == 1L && identical(fs[[1L]]$key, t$key)
}

# expr.atoms(AppliedUndef, Derivative), in preorder
.pm_undef_deriv_atoms <- function(e) {
  out <- list()
  seen <- character()
  walk <- function(x) {
    if (!is_sym(x)) return(invisible(NULL))
    if ((x$type == "Derivative" || .is_undefined_application(x)) && !(x$key %in% seen)) {
      seen <<- c(seen, x$key)
      out[[length(out) + 1L]] <<- x
    }
    for (a in x$args) walk(a)
  }
  if (inherits(e, "symmatrix")) for (z in e$data) walk(z) else walk(as_sym(e))
  out
}

#' Dynamic symbols of an expression
#'
#' `find_dynamicsymbols()` returns the functions of time and their
#' derivatives found in `expression`, as SymPy's `find_dynamicsymbols()`.
#' A vector needs the frame in which to express it.
#'
#' @param expression An expression, matrix or vector.
#' @param exclude Dynamic symbols to leave out.
#' @param reference_frame The [ReferenceFrame()] for a vector expression.
#' @return A list of expressions (a set in SymPy).
#' @examples
#' q <- dynamicsymbols("q")
#' find_dynamicsymbols(q^2 + dynamicsymbols("q", 1))
#' @export
find_dynamicsymbols <- function(expression, exclude = NULL, reference_frame = NULL) {
  t <- .pv_t()
  exclude_set <- list()
  if (length(exclude)) {
    if (is_sym(exclude) || !(is.list(exclude) || inherits(exclude, "symmatrix")))
      .pv_err("TypeError", "exclude kwarg must be iterable")
    exclude_set <- .pv_seq(exclude)
  }
  if (.pv_is_vec(expression)) {
    if (is.null(reference_frame))
      .pv_err("ValueError", "You must provide reference_frame when passing a vector expression, got None.")
    expression <- .pv_to_matrix(expression, reference_frame)
  }
  ex <- vapply(exclude_set, function(z) as_sym(z)$key, "")
  res <- Filter(function(i) .pm_is_dynsym(i, t) && !(i$key %in% ex), .pm_undef_deriv_atoms(expression))
  structure(res, pyset = TRUE)
}

# ---------------------------------------------------------------------------
# msubs

# dict_merge of the substitution dictionaries: later keys win
.pm_merge_dicts <- function(dicts) {
  keys <- list()
  vals <- list()
  for (d in dicts) {
    pairs <- if (inherits(d, "pydict")) Map(list, d$keys, d$values) else .subs_pairs(d)
    for (p in pairs) {
      k <- if (is.character(p[[1L]])) Symbol(p[[1L]]) else as_sym(p[[1L]])
      i <- match(k$key, vapply(keys, .k, ""))
      if (is.na(i)) {
        keys[[length(keys) + 1L]] <- k
        vals[[length(vals) + 1L]] <- as_sym(p[[2L]])
      } else vals[[i]] <- as_sym(p[[2L]])
    }
  }
  list(keys = keys, vals = vals)
}

.pm_crawl <- function(expr, func) {
  val <- func(expr)
  if (!is.null(val)) return(val)
  args <- lapply(expr$args, .pm_crawl, func = func)
  if (identical(vapply(args, .k, ""), vapply(expr$args, .k, ""))) expr else sym_rebuild(expr, args)
}

.pm_sub_func <- function(expr, sd) {
  i <- match(expr$key, vapply(sd$keys, .k, ""))
  if (!is.na(i)) return(sd$vals[[i]])
  if (!length(expr$args) || expr$type == "Derivative") return(expr)
  NULL
}

.pm_tan_repl_func <- function(expr) {
  if (expr$type == "Function" && identical(expr$name, "tan"))
    return(fcall("sin", expr$args) / fcall("cos", expr$args))
  if (!length(expr$args) || expr$type == "Derivative") return(expr)
  NULL
}

# `a < 0` of an exponent, as the truth value SymPy needs
.pm_lt0 <- function(a) {
  if (is_number_type(a)) return(num_sign(a) < 0)
  neg <- ask_fact(a, "extended_negative")
  if (isTRUE(neg)) return(TRUE)
  if (isTRUE(ask_fact(a, "extended_nonnegative"))) return(FALSE)
  .pv_err("TypeError", sprintf("cannot determine truth value of Relational: %s < 0", .str(a)))
}

.pm_fraction_decomp <- function(expr) {
  if (expr$type != "Mul") return(list(expr, S_One))
  num <- list()
  den <- list()
  for (a in .sympy_args(expr)) {
    if (a$type == "Pow" && .pm_lt0(a$args[[2L]])) den[[length(den) + 1L]] <- .div(S_One, a)
    else num[[length(num) + 1L]] <- a
  }
  if (!length(den)) return(list(expr, S_One))
  list(.mul(num), .mul(den))
}

.pm_smart_subs <- function(expr, sd) {
  expr <- .pm_crawl(expr, .pm_tan_repl_func)
  recurser <- function(expr) {
    nd <- .pm_fraction_decomp(expr)
    if (!identical(nd[[2L]]$key, S_One$key)) {
      denom_subbed <- recurser(nd[[2L]])
      if (identical(N(denom_subbed)$key, S_Zero$key)) {
        expr <- simplify(expr)
      } else {
        num_subbed <- recurser(nd[[1L]])
        return(num_subbed / denom_subbed)
      }
    }
    val <- .pm_sub_func(expr, sd)
    if (!is.null(val)) return(val)
    args <- lapply(expr$args, recurser)
    if (identical(vapply(args, .k, ""), vapply(expr$args, .k, ""))) expr else sym_rebuild(expr, args)
  }
  recurser(expr)
}

#' Substitution for mechanics expressions
#'
#' `msubs()` substitutes into expressions, matrices, vectors and dyadics as
#' SymPy's `msubs()`: a structural replacement that leaves derivatives
#' whole unless they appear as keys. With `smart = TRUE` tangents are
#' rewritten and denominators that vanish are simplified first.
#'
#' @param expr An expression, matrix, vector or dyadic.
#' @param ... Substitution dictionaries: lists of `list(old, new)` pairs.
#' @param smart Use the smart substitution.
#' @return The substituted object.
#' @examples
#' q <- dynamicsymbols("q")
#' msubs(q^2 + dynamicsymbols("q", 1), list(list(q, 2)))
#' @export
msubs <- function(expr, ..., smart = FALSE) {
  sd <- .pm_merge_dicts(list(...))
  func <- if (smart) function(x) .pm_smart_subs(x, sd) else function(x) .pm_crawl(x, function(z) .pm_sub_func(z, sd))
  if (inherits(expr, "symmatrix")) return(.pv_mmap(expr, func))
  if (inherits(expr, c("PhysVector", "PhysDyadic"))) return(expr$applyfunc(func))
  func(as_sym(expr))
}

#' @rdname vsprint
#' @export
mprint <- function(expr, ...) vprint(expr, ...)

#' @rdname vsprint
#' @export
msprint <- function(expr, ...) vsprint(expr, ...)

#' @rdname vsprint
#' @export
mpprint <- function(expr, use_unicode = NULL, ...) vpprint(expr, use_unicode = use_unicode, ...)

#' @rdname vsprint
#' @export
mlatex <- function(expr, ...) vlatex(expr, ...)

#' @rdname vsprint
#' @export
mechanics_printing <- function(...) init_vprinting(...)

# ---------------------------------------------------------------------------
# Private helpers of the equation of motion classes

.pm_f_list_parser <- function(fl, ref_frame) {
  if (!length(fl)) return(list(list(), list()))
  vel_list <- list()
  f_list <- list()
  for (pair in fl) {
    obj <- pair[[1L]]
    force <- pair[[2L]]
    if (.pv_is_frame(obj)) v <- .pv_ang_vel_in(obj, ref_frame)
    else if (.pv_is_point(obj)) v <- .pv_point_vel(obj, ref_frame)
    else .pv_err("TypeError", "First entry in each forcelist pair must be a point or frame.")
    vel_list[[length(vel_list) + 1L]] <- v
    f_list[[length(f_list) + 1L]] <- force
  }
  list(vel_list, f_list)
}

.pm_iterable <- function(x) {
  if (is.null(x)) return(list())
  if (inherits(x, "symmatrix")) return(x$data)
  if (is.list(x) && !is_sym(x) && !inherits(x, "physics_obj")) return(unclass(x))
  list(x)
}

.pm_set_str <- function(xs) paste0("{", paste(vapply(xs, .pv_str, ""), collapse = ", "), "}")

.pm_dups <- function(xs) {
  seen <- character()
  out <- list()
  for (x in xs) {
    k <- as_sym(x)$key
    if (k %in% seen) {
      if (!(k %in% vapply(out, .k, ""))) out[[length(out) + 1L]] <- as_sym(x)
    } else seen <- c(seen, k)
  }
  out
}

.pm_intersect <- function(a, b) {
  bk <- vapply(b, function(z) as_sym(z)$key, "")
  out <- list()
  for (x in a) {
    k <- as_sym(x)$key
    if (k %in% bk && !(k %in% vapply(out, .k, ""))) out[[length(out) + 1L]] <- as_sym(x)
  }
  out
}

.pm_validate_coordinates <- function(coordinates = NULL, speeds = NULL, check_duplicates = TRUE,
                                     is_dynamicsymbols = TRUE, u_auxiliary = NULL) {
  t <- .pv_t()
  coordinates <- .pm_iterable(coordinates)
  speeds <- .pm_iterable(speeds)
  u_auxiliary <- .pm_iterable(u_auxiliary)
  msgs <- character()
  if (check_duplicates) {
    coord_duplicates <- .pm_dups(coordinates)
    speed_duplicates <- .pm_dups(speeds)
    aux_duplicates <- .pm_dups(u_auxiliary)
    overlap_coords <- .pm_intersect(coordinates, speeds)
    overlap_aux <- .pm_intersect(c(coordinates, speeds), u_auxiliary)
    if (length(coord_duplicates))
      msgs <- c(msgs, sprintf("The generalized coordinates %s are duplicated, all generalized coordinates should be unique.",
                              .pm_set_str(coord_duplicates)))
    if (length(speed_duplicates))
      msgs <- c(msgs, sprintf("The generalized speeds %s are duplicated, all generalized speeds should be unique.",
                              .pm_set_str(speed_duplicates)))
    if (length(aux_duplicates))
      msgs <- c(msgs, sprintf("The auxiliary speeds %s are duplicated, all auxiliary speeds should be unique.",
                              .pm_set_str(aux_duplicates)))
    if (length(overlap_coords))
      msgs <- c(msgs, sprintf("%s are defined as both generalized coordinates and generalized speeds.",
                              .pm_set_str(overlap_coords)))
    if (length(overlap_aux))
      msgs <- c(msgs, sprintf("The auxiliary speeds %s are also defined as generalized coordinates or generalized speeds.",
                              .pm_set_str(overlap_aux)))
  }
  if (is_dynamicsymbols) {
    for (coordinate in coordinates) if (!.pm_is_dynsym(as_sym(coordinate), t))
      msgs <- c(msgs, sprintf("Generalized coordinate \"%s\" is not a dynamicsymbol.", .pv_str(coordinate)))
    for (speed in speeds) if (!.pm_is_dynsym(as_sym(speed), t))
      msgs <- c(msgs, sprintf("Generalized speed \"%s\" is not a dynamicsymbol.", .pv_str(speed)))
    for (aux in u_auxiliary) if (!.pm_is_dynsym(as_sym(aux), t))
      msgs <- c(msgs, sprintf("Auxiliary speed \"%s\" is not a dynamicsymbol.", .pv_str(aux)))
  }
  if (length(msgs)) .pv_err("ValueError", paste(msgs, collapse = "\n"))
  invisible(NULL)
}

.pm_parse_linear_solver <- function(linear_solver) {
  if (is.function(linear_solver)) return(linear_solver)
  function(A, b) .pm_matrix_solve(A, b, linear_solver)
}

.pm_matrix_solve <- function(A, b, method = "GJ") {
  if (identical(method, "LU")) return(.pm_LUsolve(A, b))
  solve(A, b)
}

# Matrix.LUdecomposition_Simple(rankcheck = TRUE), as in
# sympy/matrices/decompositions.py
.pm_LU_simple <- function(M) {
  rows <- M$nrow
  cols <- M$ncol
  lu <- .mrows(M)
  row_swaps <- list()
  if (rows == 0L || cols == 0L) return(list(lu = lu, swaps = row_swaps))
  rank_err <- function() .pv_err("ValueError", "Rank of matrix is strictly less than number of rows or columns. Pass keyword argument rankcheck=False to compute the LU decomposition of this matrix.")
  pivot_col <- 0L
  for (pivot_row in seq_len(rows - 1L) - 1L) {
    iszeropivot <- TRUE
    pivot_row_offset <- NULL
    pairs <- list()
    while (pivot_col != cols && iszeropivot) {
      pv <- NULL
      pairs <- list()
      ind <- list()
      for (i in seq.int(0L, rows - pivot_row - 1L)) {
        val <- lu[[pivot_row + i + 1L]][[pivot_col + 1L]]
        z <- is_zero(val)
        if (identical(z, FALSE)) { pv <- list(i, val); break }
        if (is.na(z)) ind[[length(ind) + 1L]] <- list(i, val)
      }
      if (is.null(pv) && length(ind)) {
        # simpfunc is None in LUdecomposition_Simple: the first undecided entry
        pv <- ind[[1L]]
      }
      iszeropivot <- is.null(pv)
      pivot_row_offset <- if (is.null(pv)) NULL else pv[[1L]]
      if (iszeropivot) pivot_col <- pivot_col + 1L
    }
    if (pivot_col != pivot_row) rank_err()
    candidate <- if (is.null(pivot_row_offset)) NULL else pivot_row + pivot_row_offset
    if (is.null(candidate) && iszeropivot) return(list(lu = lu, swaps = row_swaps))
    if (pivot_row != candidate) {
      row_swaps[[length(row_swaps) + 1L]] <- c(pivot_row, candidate)
      a <- lu[[pivot_row + 1L]]
      b <- lu[[candidate + 1L]]
      lu[[pivot_row + 1L]] <- b
      lu[[candidate + 1L]] <- a
    }
    start_col <- pivot_col + 1L
    if (pivot_row + 1L < rows) for (row in seq.int(pivot_row + 1L, rows - 1L)) {
      lu[[row + 1L]][[pivot_row + 1L]] <- .div(lu[[row + 1L]][[pivot_col + 1L]], lu[[pivot_row + 1L]][[pivot_col + 1L]])
      if (start_col < cols) for (c in seq.int(start_col, cols - 1L))
        lu[[row + 1L]][[c + 1L]] <- .sub(lu[[row + 1L]][[c + 1L]], .mul(list(lu[[row + 1L]][[pivot_row + 1L]], lu[[pivot_row + 1L]][[c + 1L]])))
    }
    pivot_col <- pivot_col + 1L
    if (pivot_col == cols) return(list(lu = lu, swaps = row_swaps))
  }
  k <- min(rows, cols)
  if (isTRUE(is_zero(lu[[k]][[k]]))) rank_err()
  list(lu = lu, swaps = row_swaps)
}

# Matrix.LUsolve, as in sympy/matrices/solvers.py
.pm_LUsolve <- function(M, rhs) {
  if (rhs$nrow != M$nrow) .pv_err("ShapeError", "``M`` and ``rhs`` must have the same number of rows.")
  m <- M$nrow
  n <- M$ncol
  if (m < n) .pv_err("NotImplementedError", "Underdetermined systems not supported.")
  r <- tryCatch(.pm_LU_simple(M), ValueError = function(e) .pv_err("NonInvertibleMatrixError", "Matrix det == 0; not invertible."))
  A <- r$lu
  b <- .mrows(rhs)
  for (sw in r$swaps) {
    tmp <- b[[sw[1L] + 1L]]
    b[[sw[1L] + 1L]] <- b[[sw[2L] + 1L]]
    b[[sw[2L] + 1L]] <- tmp
  }
  nc <- rhs$ncol
  op <- function(i, j, scale) b[[i]] <<- lapply(seq_len(nc), function(k) .sub(b[[i]][[k]], .mul(list(scale, b[[j]][[k]]))))
  for (i in seq_len(m)) for (j in seq_len(min(i - 1L, n))) op(i, j, A[[i]][[j]])
  if (m > n) {
    for (i in seq.int(n + 1L, m)) for (k in seq_len(nc))
      if (!isTRUE(is_zero(b[[i]][[k]]))) .pv_err("ValueError", "The system is inconsistent.")
    b <- b[seq_len(n)]
  }
  for (i in rev(seq_len(n))) {
    if (i < n) for (j in seq.int(i + 1L, n)) op(i, j, A[[i]][[j]])
    scale <- A[[i]][[i]]
    b[[i]] <- lapply(b[[i]], function(x) .mul(list(.pow(scale, S_NegativeOne), x)))
  }
  if (!length(b)) return(.new_mat(list(), 0L, nc))
  .from_rows(b)
}
