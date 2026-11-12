# Classical mechanics vectors, after sympy/physics/vector/frame.py.
#
# SymPy's physics.vector and sympy.vector share several names. In this
# package the physics.vector classes are ReferenceFrame, PhysVector,
# PhysDyadic and PhysPoint; dot(), cross(), express(), outer_sym(), curl(),
# divergence(), gradient(), is_conservative(), is_solenoidal(),
# scalar_potential() and scalar_potential_difference() dispatch on the class
# of their first argument, so the same names serve both modules.
#
# Frames, points and the mechanics objects are mutable in SymPy and compared
# by identity; here they are environments. Vectors and dyadics are immutable
# lists. Python methods and properties are reached with `$`: a property is a
# function .pvP_<class>_<name>(x), a method .pvM_<class>_<name>(x, ...) and a
# setter .pvS_<class>_<name>(x, value), looked up along the S3 class vector.

.pv_ns <- environment()

.pv_state <- new.env()
.pv_state$frame_count <- 0L
.pv_state$uid <- 0L
.pv_state$simp <- FALSE
.pv_state$dstr <- "'"

# dynamicsymbols._t: the time symbol, a class attribute SymPy lets callers reset
.pv_t <- function() if (is.null(.pv_state$t)) Symbol("t") else .pv_state$t

.pv_next_uid <- function() {
  .pv_state$uid <- .pv_state$uid + 1L
  .pv_state$uid
}

.PV_EXC_PARENTS <- list(VectorTypeError = "TypeError", ArgumentIndexError = "ValueError",
                        FrozenInstanceError = "AttributeError")

.pv_err <- function(cls, msg) {
  parents <- .PV_EXC_PARENTS[[cls]]
  if (is.null(parents)) parents <- .SYMPY_EXCEPTIONS[[cls]]
  stop(structure(list(message = paste(msg, collapse = ""), call = NULL),
                 class = c(cls, parents, "sympy_error", "error", "condition")))
}

# textwrap.fill(width = 70) behind SymPy's filldedent()
.pv_fill <- function(s) {
  words <- strsplit(trimws(gsub("\\s+", " ", s)), " ", fixed = TRUE)[[1L]]
  lines <- character()
  cur <- ""
  for (w in words) {
    if (!nzchar(cur)) cur <- w
    else if (nchar(cur) + 1L + nchar(w) <= 70L) cur <- paste(cur, w)
    else {
      lines <- c(lines, cur)
      cur <- w
    }
  }
  if (nzchar(cur)) lines <- c(lines, cur)
  paste0("\n", paste(lines, collapse = "\n"))
}

# Python's type() of a value, as printed in error messages
.pv_pyclass <- function(v) {
  if (is.null(v)) return("<class 'NoneType'>")
  if (is.character(v)) return("<class 'str'>")
  if (is.logical(v)) return("<class 'bool'>")
  if (is.numeric(v)) return(if (v == round(v)) "<class 'int'>" else "<class 'float'>")
  if (is_sym(v)) {
    special <- if (identical(v$key, S_Zero$key)) "Zero" else if (identical(v$key, S_One$key)) "One"
               else if (identical(v$key, S_NegativeOne$key)) "NegativeOne" else if (identical(v$key, S_Half$key)) "Half"
    if (!is.null(special)) return(paste0("<class 'sympy.core.numbers.", special, "'>"))
    mod <- switch(v$type,
                  Symbol = "sympy.core.symbol.Symbol", Dummy = "sympy.core.symbol.Dummy",
                  Integer = , Rational = , Float = paste0("sympy.core.numbers.", v$type),
                  Add = "sympy.core.add.Add", Mul = "sympy.core.mul.Mul",
                  Pow = "sympy.core.power.Pow", NULL)
    if (is.null(mod) && v$type == "Constant") {
      mod <- switch(v$name, pi = "sympy.core.numbers.Pi", E = "sympy.core.numbers.Exp1",
                    paste0("sympy.core.numbers.", v$name))
    }
    if (is.null(mod) && v$type == "Function") {
      if (.is_undefined_application(v)) return(v$name)
      mod <- paste0("sympy.functions.", v$name)
    }
    if (is.null(mod)) mod <- paste0("sympy.core.", v$type)
    return(paste0("<class '", mod, "'>"))
  }
  if (inherits(v, "symmatrix")) return("<class 'sympy.matrices.immutable.ImmutableDenseMatrix'>")
  cls <- class(v)[1L]
  mod <- .PV_PYMODULE[[cls]]
  if (!is.null(mod)) return(paste0("<class '", mod, "'>"))
  if (is.list(v)) return("<class 'list'>")
  paste0("<class '", cls, "'>")
}

.PV_PYMODULE <- list(
  ReferenceFrame = "sympy.physics.vector.frame.ReferenceFrame",
  PhysVector = "sympy.physics.vector.vector.Vector",
  PhysDyadic = "sympy.physics.vector.dyadic.Dyadic",
  PhysPoint = "sympy.physics.vector.point.Point",
  Particle = "sympy.physics.mechanics.particle.Particle",
  RigidBody = "sympy.physics.mechanics.rigidbody.RigidBody",
  Inertia = "sympy.physics.mechanics.inertia.Inertia",
  Force = "sympy.physics.mechanics.loads.Force",
  Torque = "sympy.physics.mechanics.loads.Torque",
  Body = "sympy.physics.mechanics.body.Body",
  PinJoint = "sympy.physics.mechanics.joint.PinJoint",
  PrismaticJoint = "sympy.physics.mechanics.joint.PrismaticJoint",
  CylindricalJoint = "sympy.physics.mechanics.joint.CylindricalJoint",
  PlanarJoint = "sympy.physics.mechanics.joint.PlanarJoint",
  SphericalJoint = "sympy.physics.mechanics.joint.SphericalJoint",
  WeldJoint = "sympy.physics.mechanics.joint.WeldJoint",
  KanesMethod = "sympy.physics.mechanics.kane.KanesMethod",
  LagrangesMethod = "sympy.physics.mechanics.lagrange.LagrangesMethod",
  System = "sympy.physics.mechanics.system.System",
  ZerothOrderActivation = "sympy.physics.biomechanics.activation.ZerothOrderActivation",
  FirstOrderActivationDeGroote2016 = "sympy.physics.biomechanics.activation.FirstOrderActivationDeGroote2016",
  MusculotendonDeGroote2016 = "sympy.physics.biomechanics.musculotendon.MusculotendonDeGroote2016",
  CharacteristicCurveCollection = "sympy.physics.biomechanics.curve.CharacteristicCurveCollection",
  pytuple = "tuple"
)

# Python's str() of a value
.pv_str <- function(v) {
  if (is.null(v)) return("None")
  if (inherits(v, "physics_obj")) return(format(v))
  if (is_sym(v) || inherits(v, "symmatrix")) return(.pv_sstr(v))
  if (is.logical(v) && length(v) == 1L) return(if (v) "True" else "False")
  if (is.character(v) && length(v) == 1L) return(v)
  if (is.numeric(v) && length(v) == 1L) return(format(v))
  if (is.list(v)) return(paste0("[", paste(vapply(v, .pv_repr, ""), collapse = ", "), "]"))
  paste(format(v), collapse = " ")
}

.pv_repr <- function(v) if (is.character(v) && length(v) == 1L) paste0("'", v, "'") else .pv_str(v)

# sstr() honoring the active vector printer
.pv_sstr <- function(e) {
  if (inherits(e, "physics_obj")) return(format(e))
  if (inherits(e, "symmatrix")) return(sstr(e))
  .str_top(as_sym(e))
}

# ---------------------------------------------------------------------------
# Attribute access: properties, fields and methods.

.pv_missing <- structure(list(), class = "pv_missing")

.pv_lookup <- function(x, prefix, name) {
  for (cl in class(x)) {
    f <- get0(paste0(prefix, cl, "_", name), envir = .pv_ns, inherits = FALSE)
    if (!is.null(f)) return(f)
  }
  NULL
}

.pv_getattr <- function(x, name) {
  p <- .pv_lookup(x, ".pvP_", name)
  if (!is.null(p)) return(p(x))
  v <- if (is.environment(x)) get0(name, envir = x, inherits = FALSE, ifnotfound = .pv_missing)
       else if (name %in% names(unclass(x))) .subset2(x, name) else .pv_missing
  if (!inherits(v, "pv_missing")) return(v)
  m <- .pv_lookup(x, ".pvM_", name)
  if (!is.null(m)) return(function(...) m(x, ...))
  .pv_err("AttributeError", sprintf("'%s' object has no attribute '%s'", .pv_pyname(x), name))
}

.pv_pyname <- function(x) {
  cl <- class(x)[1L]
  switch(cl, PhysVector = "Vector", PhysDyadic = "Dyadic", PhysPoint = "Point", cl)
}

#' Attribute access for classical mechanics objects
#'
#' Frames, vectors, dyadics, points and the mechanics classes expose the
#' attributes, properties and methods of their SymPy counterparts through
#' `$`, so `N$x` is the unit vector of frame `N` and `v$dt(N)` is the time
#' derivative of `v` in `N`. Assigning with `$<-` runs a property setter
#' where SymPy has one.
#'
#' @param x A physics object.
#' @param name Attribute name.
#' @param value New value.
#' @return The attribute value, or a function for a method.
#' @examples
#' N <- ReferenceFrame("N")
#' N$x
#' N$name
#' @name physics-attributes
#' @export
`$.physics_obj` <- function(x, name) .pv_getattr(x, name)

#' @rdname physics-attributes
#' @export
`$<-.physics_obj` <- function(x, name, value) {
  s <- .pv_lookup(x, ".pvS_", name)
  if (!is.null(s)) {
    s(x, value)
    return(invisible(x))
  }
  if (!is.environment(x)) .pv_err("AttributeError", sprintf("can't set attribute '%s'", name))
  if (!is.null(.pv_lookup(x, ".pvP_", name)))
    .pv_err("AttributeError", sprintf("property '%s' of '%s' object has no setter", name, .pv_pyname(x)))
  assign(name, value, envir = x)
  invisible(x)
}

#' @export
print.physics_obj <- function(x, ...) {
  cat(format(x), "\n", sep = "")
  invisible(x)
}

# ---------------------------------------------------------------------------
# Insertion ordered dictionaries keyed by frames, points and other objects:
# a named list of list(key, val) entries, named by the key's id.

.pv_id <- function(o) {
  if (is.environment(o)) return(paste0("#", get("_uid", envir = o)))
  if (is_sym(o)) return(o$key)
  if (is.character(o)) return(paste0("'", o))
  stop("unhashable key")
}

.pvod_has <- function(d, k) !is.null(d[[.pv_id(k)]])
.pvod_get <- function(d, k) {
  e <- d[[.pv_id(k)]]
  if (is.null(e)) NULL else e$val
}
.pvod_set <- function(d, k, v) {
  d[[.pv_id(k)]] <- list(key = k, val = v)
  d
}
.pvod_del <- function(d, k) {
  d[[.pv_id(k)]] <- NULL
  d
}
.od_keys <- function(d) lapply(unname(d), `[[`, "key")
.od_vals <- function(d) lapply(unname(d), `[[`, "val")

.od_pydict <- function(d) .pydict(.od_keys(d), .od_vals(d))

# ---------------------------------------------------------------------------
# Small 3 x 3 matrix helpers that build entries as SymPy's dense matrix
# products and sums do.

.pv_col <- function(...) .new_mat(lapply(list(...), as_sym), 3L, 1L)
.pv_colv <- function(xs) .new_mat(lapply(xs, as_sym), length(xs), 1L)
.pv_rows <- function(rows) .from_rows(lapply(rows, function(r) lapply(r, as_sym)))
.pv_e <- function(M, i) M$data[[i]]
.pv_T <- function(M) {
  d <- vector("list", length(M$data))
  for (i in seq_len(M$nrow)) for (j in seq_len(M$ncol))
    d[[(j - 1L) * M$nrow + i]] <- M$data[[(i - 1L) * M$ncol + j]]
  .new_mat(d, M$ncol, M$nrow)
}
.pv_matmul <- function(A, B) {
  if (A$ncol != B$nrow) .pv_err("ShapeError", sprintf("Matrix size mismatch: (%d, %d) * (%d, %d).", A$nrow, A$ncol, B$nrow, B$ncol))
  d <- vector("list", A$nrow * B$ncol)
  for (i in seq_len(A$nrow)) for (j in seq_len(B$ncol)) {
    terms <- lapply(seq_len(A$ncol), function(k) .mul(list(A$data[[(i - 1L) * A$ncol + k]], B$data[[(k - 1L) * B$ncol + j]])))
    d[[(i - 1L) * B$ncol + j]] <- .add(terms)
  }
  .new_mat(d, A$nrow, B$ncol)
}
.pv_madd <- function(A, B) .new_mat(Map(function(a, b) .add(list(a, b)), A$data, B$data), A$nrow, A$ncol)
.pv_msub <- function(A, B) .new_mat(Map(.sub, A$data, B$data), A$nrow, A$ncol)
.pv_mscale <- function(M, s) .new_mat(lapply(M$data, function(a) .mul(list(s, a))), M$nrow, M$ncol)
.pv_mmap <- function(M, f) .new_mat(lapply(M$data, f), M$nrow, M$ncol)
.pv_meq <- function(A, B) A$nrow == B$nrow && A$ncol == B$ncol &&
  all(vapply(seq_along(A$data), function(k) identical(A$data[[k]]$key, B$data[[k]]$key), TRUE))
.pv_mzero <- function(M) all(vapply(M$data, function(a) identical(a$key, S_Zero$key), TRUE))
.pv_eye3 <- function() .pv_rows(list(list(1, 0, 0), list(0, 1, 0), list(0, 0, 1)))
.pv_zeros <- function(n, m) .new_mat(rep(list(S_Zero), n * m), n, m)
.pv_mdiff <- function(M, x) .pv_mmap(M, function(a) .pv_diff(a, x))

.pv_diff <- function(e, x) diff(as_sym(e), as_sym(x))

# Python's `x == 0` for the values found in this module
.pv_is0 <- function(x) {
  if (is.null(x)) return(FALSE)
  if (inherits(x, c("PhysVector", "PhysDyadic"))) return(!length(.subset2(x, "args")))
  if (is.numeric(x) && length(x) == 1L) return(x == 0)
  if (is_sym(x)) return(identical(x$key, S_Zero$key))
  FALSE
}

.pv_is_frame <- function(x) inherits(x, "ReferenceFrame")
.pv_is_vec <- function(x) inherits(x, "PhysVector")
.pv_is_dyad <- function(x) inherits(x, "PhysDyadic")
.pv_is_point <- function(x) inherits(x, "PhysPoint")

.pv_check_frame <- function(other) {
  if (!.pv_is_frame(other)) {
    ReferenceFrame("A")
    .pv_err("VectorTypeError", .pv_fill(sprintf("Expected an instance of %s, but received object '%s' of %s.",
                                               .PV_FRAME_CLASS, .pv_str(other), .pv_pyclass(other))))
  }
  invisible(other)
}
.PV_FRAME_CLASS <- "<class 'sympy.physics.vector.frame.ReferenceFrame'>"

.pv_check_vector <- function(other) {
  if (!.pv_is_vec(other)) .pv_err("TypeError", "A Vector must be supplied")
  other
}

.pv_check_dyadic <- function(other) {
  if (!.pv_is_dyad(other)) .pv_err("TypeError", "A Dyadic must be supplied")
  other
}

# Python's translate(str(order), 'XYZxyz', '123123')
.pv_rot_order <- function(o) chartr("XYZxyz", "123123", if (is.numeric(o)) format(o, scientific = FALSE) else as.character(o))

.PV_APPROVED_ORDERS <- c("123", "231", "312", "132", "213", "321", "121",
                         "131", "212", "232", "313", "323", "")

.pv_seq <- function(x) {
  if (is.null(x)) return(list())
  if (is_sym(x) || inherits(x, "physics_obj")) return(list(x))
  if (inherits(x, "symmatrix")) return(x$data)
  as.list(x)
}

# Python's list(x) and x[i] (i zero based), with Python's errors for scalars.
.pv_tname <- function(v) sub("^<class '(.*[.])?([^.']*)'>$", "\\2", .pv_pyclass(v))

.pv_pylist <- function(x) {
  if (is_sym(x) || (is.numeric(x) && length(x) == 1L))
    .pv_err("TypeError", sprintf("'%s' object is not iterable", .pv_tname(x)))
  .pv_seq(x)
}

.pv_pyindex <- function(x, i) {
  if (is_sym(x) || (is.numeric(x) && length(x) == 1L))
    .pv_err("TypeError", sprintf("'%s' object is not subscriptable", .pv_tname(x)))
  s <- .pv_seq(x)
  if (i >= length(s)) .pv_err("IndexError", sprintf("%s index out of range", if (inherits(x, "pytuple")) "tuple" else "list"))
  s[[i + 1L]]
}

# ---------------------------------------------------------------------------
# CoordinateSym

#' Coordinate variables of a reference frame
#'
#' A `CoordinateSym` is a symbol tied to a [ReferenceFrame()] and one of its
#' three axes, as SymPy's `sympy.physics.vector.CoordinateSym`. Each frame
#' creates its own three; `N[1]` returns the first. Two coordinate symbols
#' are equal only when frame and index agree.
#'
#' @param name Display name.
#' @param frame The `ReferenceFrame` it belongs to.
#' @param index The axis, 0, 1 or 2 as in SymPy.
#' @return A `sym` symbol.
#' @examples
#' A <- ReferenceFrame("A")
#' A[1]
#' CoordinateSym("A_x", A, 0)
#' @export
CoordinateSym <- function(name, frame, index) {
  .pv_check_frame(frame)
  if (!(index %in% 0:2)) .pv_err("ValueError", "Invalid index specified")
  uid <- get("_uid", envir = frame)
  facts <- .symbol_facts(logical())
  new_sym("Symbol", name = name, extra = list(facts = facts, given = logical(), frame = frame, index = as.integer(index)),
          key = paste0("S", nchar(name), ":", name, .assumption_suffix(facts), "@F", uid, ".", index))
}

.pv_is_coordsym <- function(e) is_sym(e) && e$type == "Symbol" && !is.null(e$extra$frame)

# ---------------------------------------------------------------------------
# ReferenceFrame

.pv_check_str3 <- function(v, what, msg_type) {
  if (!(is.character(v) || is.list(v))) .pv_err("TypeError", what)
  if (length(v) != 3L) .pv_err("ValueError", "Supply 3 indices")
  if (!all(vapply(as.list(v), is.character, TRUE))) .pv_err("TypeError", msg_type)
  vapply(as.list(v), identity, "")
}

#' Reference frames for classical mechanics
#'
#' `ReferenceFrame()` creates a frame with three orthonormal unit vectors, as
#' SymPy's `sympy.physics.vector.ReferenceFrame`. Frames are mutable and
#' compared by identity. Unit vectors are `N$x`, `N$y` and `N$z` (or
#' `N["x"]`), and `N[1]` to `N[3]` are the coordinate variables. Orient a
#' frame with the methods `orient_axis`, `orient_dcm`, `orient_explicit`,
#' `orient_body_fixed`, `orient_space_fixed`, `orient_quaternion`, `orient`
#' and `orientnew`, and query it with `dcm`, `ang_vel_in`, `ang_acc_in`,
#' `variable_map` and `partial_velocity`, all reached with `$`.
#'
#' @param name Frame name.
#' @param indices Optional three strings naming the unit vectors.
#' @param latexs Optional three LaTeX strings for the unit vectors.
#' @param variables Optional three names for the coordinate variables.
#' @return A `ReferenceFrame`.
#' @examples
#' N <- ReferenceFrame("N")
#' q <- dynamicsymbols("q")
#' A <- N$orientnew("A", "Axis", list(q, N$z))
#' A$dcm(N)
#' A$ang_vel_in(N)
#' @export
ReferenceFrame <- function(name, indices = NULL, latexs = NULL, variables = NULL) {
  if (!is.character(name) || length(name) != 1L) .pv_err("TypeError", "Need to supply a valid name")
  fr <- new.env(parent = emptyenv())
  assign("_uid", .pv_next_uid(), envir = fr)
  if (!is.null(indices)) {
    if (!(is.character(indices) || is.list(indices))) .pv_err("TypeError", "Supply the indices as a list")
    if (length(indices) != 3L) .pv_err("ValueError", "Supply 3 indices")
    if (!all(vapply(as.list(indices), is.character, TRUE))) .pv_err("TypeError", "Indices must be strings")
    indices <- vapply(as.list(indices), identity, "")
    fr$str_vecs <- paste0(name, "['", indices, "']")
    fr$pretty_vecs <- paste0(tolower(name), "_", indices)
    fr$latex_vecs <- sprintf("\\mathbf{\\hat{%s}_{%s}}", tolower(name), indices)
    fr$indices <- indices
  } else {
    fr$str_vecs <- paste0(name, c(".x", ".y", ".z"))
    fr$pretty_vecs <- paste0(tolower(name), c("_x", "_y", "_z"))
    fr$latex_vecs <- sprintf("\\mathbf{\\hat{%s}_%s}", tolower(name), c("x", "y", "z"))
    fr$indices <- c("x", "y", "z")
  }
  if (!is.null(latexs)) {
    if (!(is.character(latexs) || is.list(latexs))) .pv_err("TypeError", "Supply the indices as a list")
    if (length(latexs) != 3L) .pv_err("ValueError", "Supply 3 indices")
    if (!all(vapply(as.list(latexs), is.character, TRUE))) .pv_err("TypeError", "Latex entries must be strings")
    fr$latex_vecs <- vapply(as.list(latexs), identity, "")
  }
  fr$name <- name
  fr[["_var_dict"]] <- list()
  fr[["_dcm_dict"]] <- list()
  fr[["_dcm_cache"]] <- list()
  fr[["_ang_vel_dict"]] <- list()
  fr[["_ang_acc_dict"]] <- list()
  fr[["_cur"]] <- 0L
  class(fr) <- c("ReferenceFrame", "physics_obj")
  fr[["_x"]] <- .pv_vec(list(list(.pv_col(1, 0, 0), fr)))
  fr[["_y"]] <- .pv_vec(list(list(.pv_col(0, 1, 0), fr)))
  fr[["_z"]] <- .pv_vec(list(list(.pv_col(0, 0, 1), fr)))
  if (!is.null(variables)) {
    if (!(is.character(variables) || is.list(variables))) .pv_err("TypeError", "Supply the variable names as a list/tuple")
    if (length(variables) != 3L) .pv_err("ValueError", "Supply 3 variable names")
    if (!all(vapply(as.list(variables), is.character, TRUE))) .pv_err("TypeError", "Variable names must be strings")
    variables <- vapply(as.list(variables), identity, "")
  } else {
    variables <- paste0(name, c("_x", "_y", "_z"))
  }
  fr$varlist <- .pytuple(list(CoordinateSym(variables[1L], fr, 0L),
                              CoordinateSym(variables[2L], fr, 1L),
                              CoordinateSym(variables[3L], fr, 2L)))
  .pv_state$frame_count <- .pv_state$frame_count + 1L
  fr$index <- .pv_state$frame_count
  fr
}

.pv_frame_dict <- function(fr, num) get(c("_dcm_dict", "_ang_vel_dict", "_ang_acc_dict")[num + 1L], envir = fr)

#' @export
format.ReferenceFrame <- function(x, ...) get("name", envir = x)

#' @export
`[.ReferenceFrame` <- function(x, i) {
  if (!is.character(i)) {
    if (i >= 1 && i <= 3) return(get("varlist", envir = x)[[i]])
    .pv_err("ValueError", "Invalid index provided")
  }
  ind <- get("indices", envir = x)
  if (ind[1L] == i) return(x[["_x"]])
  if (ind[2L] == i) return(x[["_y"]])
  if (ind[3L] == i) return(x[["_z"]])
  .pv_err("ValueError", "Not a defined index")
}

.pv_frame_iter <- function(fr) list(fr[["_x"]], fr[["_y"]], fr[["_z"]])

# the string attributes are Python lists
.pvP_ReferenceFrame_indices <- function(self) as.list(self[["indices"]])
.pvP_ReferenceFrame_latex_vecs <- function(self) as.list(self[["latex_vecs"]])
.pvP_ReferenceFrame_str_vecs <- function(self) as.list(self[["str_vecs"]])
.pvP_ReferenceFrame_pretty_vecs <- function(self) as.list(self[["pretty_vecs"]])
.pvP_ReferenceFrame_x <- function(self) self[["_x"]]
.pvP_ReferenceFrame_y <- function(self) self[["_y"]]
.pvP_ReferenceFrame_z <- function(self) self[["_z"]]
.pvP_ReferenceFrame_xx <- function(self) .pv_outer(self[["_x"]], self[["_x"]])
.pvP_ReferenceFrame_xy <- function(self) .pv_outer(self[["_x"]], self[["_y"]])
.pvP_ReferenceFrame_xz <- function(self) .pv_outer(self[["_x"]], self[["_z"]])
.pvP_ReferenceFrame_yx <- function(self) .pv_outer(self[["_y"]], self[["_x"]])
.pvP_ReferenceFrame_yy <- function(self) .pv_outer(self[["_y"]], self[["_y"]])
.pvP_ReferenceFrame_yz <- function(self) .pv_outer(self[["_y"]], self[["_z"]])
.pvP_ReferenceFrame_zx <- function(self) .pv_outer(self[["_z"]], self[["_x"]])
.pvP_ReferenceFrame_zy <- function(self) .pv_outer(self[["_z"]], self[["_y"]])
.pvP_ReferenceFrame_zz <- function(self) .pv_outer(self[["_z"]], self[["_z"]])
.pvP_ReferenceFrame_u <- function(self) .pv_dadd(.pv_dadd(.pvP_ReferenceFrame_xx(self), .pvP_ReferenceFrame_yy(self)), .pvP_ReferenceFrame_zz(self))
.pvP_ReferenceFrame__dlist <- function(self) list(self[["_dcm_dict"]], self[["_ang_vel_dict"]], self[["_ang_acc_dict"]])

# Identity based membership in a list of environments
.pv_in <- function(x, xs) any(vapply(xs, function(y) identical(x, y), TRUE))
.pv_path_eq <- function(a, b) length(a) == length(b) && all(vapply(seq_along(a), function(i) identical(a[[i]], b[[i]]), TRUE))
.pv_path_in <- function(p, ps) any(vapply(ps, function(q) .pv_path_eq(p, q), TRUE))

# The shortest connecting path, as ReferenceFrame._dict_list and
# Point._pdict_list: every simple path from self is grown while the list is
# iterated, then the paths not ending at other are dropped.
.pv_connect_path <- function(self, other, dict_of) {
  paths <- list(list(self))
  old <- list(list())
  same <- function(a, b) length(a) == length(b) && all(vapply(seq_along(a), function(i) .pv_path_eq(a[[i]], b[[i]]), TRUE))
  while (!same(paths, old)) {
    old <- paths
    k <- 1L
    while (k <= length(paths)) {
      fl <- paths[[k]]
      for (adj in .od_keys(dict_of(fl[[length(fl)]]))) {
        if (!.pv_in(adj, fl)) {
          cp <- c(fl, list(adj))
          if (!.pv_path_in(cp, paths)) paths[[length(paths) + 1L]] <- cp
        }
      }
      k <- k + 1L
    }
  }
  for (cp in old) {
    if (!identical(cp[[length(cp)]], other)) {
      for (j in seq_along(paths)) if (.pv_path_eq(paths[[j]], cp)) { paths[[j]] <- NULL; break }
    }
  }
  if (!length(paths)) return(NULL)
  lens <- vapply(paths, length, 0L)
  paths[[order(lens, method = "radix")[1L]]]
}

.pv_dict_list <- function(self, other, num) {
  if (!(num %in% 0:2)) .pv_err("ValueError", "Valid values for num are 0, 1, or 2.")
  p <- .pv_connect_path(self, other, function(f) .pv_frame_dict(f, num))
  if (!is.null(p)) return(p)
  ct <- c("orientation", "angular velocity", "angular acceleration")[num + 1L]
  .pv_err("ValueError", sprintf("No connecting %s path found between %s and %s.", ct, self$name, other$name))
}
.pvM_ReferenceFrame__dict_list <- .pv_dict_list

.pv_w_diff_dcm <- function(self, otherframe) {
  dcm2diff <- .pv_dcm(otherframe, self)
  diffed <- .pv_mdiff(dcm2diff, .pv_t())
  angvelmat <- .pv_matmul(diffed, .pv_T(dcm2diff))
  w1 <- trigsimp(expand(.pv_e(angvelmat, 8L)))
  w2 <- trigsimp(expand(.pv_e(angvelmat, 3L)))
  w3 <- trigsimp(expand(.pv_e(angvelmat, 4L)))
  .pv_vec_new(list(list(.pv_col(w1, w2, w3), otherframe)))
}
.pvM_ReferenceFrame__w_diff_dcm <- .pv_w_diff_dcm

.pv_variable_map <- function(self, otherframe) {
  .pv_check_frame(otherframe)
  key <- paste0(.pv_id(otherframe), "|", .pv_state$simp)
  vd <- self[["_var_dict"]]
  if (!is.null(vd[[key]])) return(vd[[key]])
  vars_matrix <- .pv_matmul(.pv_dcm(self, otherframe), .pv_colv(as.list(get("varlist", envir = otherframe))))
  vl <- get("varlist", envir = self)
  keys <- list(); vals <- list()
  for (i in 1:3) {
    keys[[i]] <- vl[[i]]
    vals[[i]] <- if (.pv_state$simp) trigsimp(.pv_e(vars_matrix, i), method = "fu") else .pv_e(vars_matrix, i)
  }
  m <- .pydict(keys, vals)
  vd[[key]] <- m
  self[["_var_dict"]] <- vd
  m
}
.pvM_ReferenceFrame_variable_map <- .pv_variable_map

.pv_ang_acc_in <- function(self, otherframe) {
  .pv_check_frame(otherframe)
  d <- self[["_ang_acc_dict"]]
  if (.pvod_has(d, otherframe)) return(.pvod_get(d, otherframe))
  .pv_time_derivative(.pv_ang_vel_in(self, otherframe), otherframe)
}
.pvM_ReferenceFrame_ang_acc_in <- .pv_ang_acc_in

.pv_ang_vel_in <- function(self, otherframe) {
  .pv_check_frame(otherframe)
  flist <- .pv_dict_list(self, otherframe, 1L)
  outvec <- .pv_vec0()
  for (i in seq_len(length(flist) - 1L))
    outvec <- .pv_vadd(outvec, .pvod_get(flist[[i]][["_ang_vel_dict"]], flist[[i + 1L]]))
  outvec
}
.pvM_ReferenceFrame_ang_vel_in <- .pv_ang_vel_in

.pv_dcm <- function(self, otherframe) {
  .pv_check_frame(otherframe)
  cache <- self[["_dcm_cache"]]
  if (.pvod_has(cache, otherframe)) return(.pvod_get(cache, otherframe))
  flist <- .pv_dict_list(self, otherframe, 0L)
  outdcm <- .pv_eye3()
  for (i in seq_len(length(flist) - 1L))
    outdcm <- .pv_matmul(outdcm, .pvod_get(flist[[i]][["_dcm_dict"]], flist[[i + 1L]]))
  self[["_dcm_cache"]] <- .pvod_set(self[["_dcm_cache"]], otherframe, outdcm)
  otherframe[["_dcm_cache"]] <- .pvod_set(otherframe[["_dcm_cache"]], self, .pv_T(outdcm))
  outdcm
}
.pvM_ReferenceFrame_dcm <- .pv_dcm

.pv_set_dcm <- function(self, parent, parent_orient) {
  frames <- .od_keys(self[["_dcm_cache"]])
  if (.pv_in(parent, frames)) {
    dcm_dict_del <- list()
    dcm_cache_del <- list()
    for (frame in frames) {
      if (.pvod_has(self[["_dcm_dict"]], frame)) dcm_dict_del[[length(dcm_dict_del) + 1L]] <- frame
      dcm_cache_del[[length(dcm_cache_del) + 1L]] <- frame
    }
    for (frame in dcm_dict_del) frame[["_dcm_dict"]] <- .pvod_del(frame[["_dcm_dict"]], self)
    for (frame in dcm_cache_del) frame[["_dcm_cache"]] <- .pvod_del(frame[["_dcm_cache"]], self)
    self[["_dcm_dict"]] <- list()
    self[["_dcm_cache"]] <- list()
  } else {
    visited <- list()
    queue <- frames
    cont <- TRUE
    while (length(queue) && cont) {
      node <- queue[[1L]]
      queue <- queue[-1L]
      if (!.pv_in(node, visited)) {
        visited[[length(visited) + 1L]] <- node
        for (neighbor in .od_keys(node[["_dcm_dict"]])) {
          if (identical(neighbor, parent)) {
            warning("Loops are defined among the orientation of frames. This is likely not desired and may cause errors in your calculations.", call. = FALSE)
            cont <- FALSE
            break
          }
          queue[[length(queue) + 1L]] <- neighbor
        }
      }
    }
  }
  pT <- .pv_T(parent_orient)
  self[["_dcm_dict"]] <- .pvod_set(self[["_dcm_dict"]], parent, pT)
  parent[["_dcm_dict"]] <- .pvod_set(parent[["_dcm_dict"]], self, parent_orient)
  self[["_dcm_cache"]] <- .pvod_set(self[["_dcm_cache"]], parent, pT)
  parent[["_dcm_cache"]] <- .pvod_set(parent[["_dcm_cache"]], self, parent_orient)
  invisible(NULL)
}
.pvM_ReferenceFrame__dcm <- .pv_set_dcm

.pv_set_wvec <- function(self, parent, wvec) {
  self[["_ang_vel_dict"]] <- .pvod_set(self[["_ang_vel_dict"]], parent, wvec)
  parent[["_ang_vel_dict"]] <- .pvod_set(parent[["_ang_vel_dict"]], self, .pv_vneg(wvec))
  self[["_var_dict"]] <- list()
  invisible(NULL)
}

.pv_orient_axis <- function(self, parent, axis, angle) {
  .pv_check_frame(parent)
  if (!.pv_is_vec(axis) && .pv_is_vec(angle)) {
    tmp <- axis; axis <- angle; angle <- tmp
  }
  axis <- .pv_check_vector(axis)
  theta <- as_sym(angle)
  if (!.pv_vec_eq(.pv_time_derivative(axis, parent), 0)) .pv_err("ValueError", "Axis cannot be time-varying.")
  unit_axis <- .pv_normalize(.pv_express(axis, parent))
  # axis.args[0][0]
  if (!length(.subset2(unit_axis, "args"))) .pv_err("IndexError", "list index out of range")
  u <- .subset2(unit_axis, "args")[[1L]][[1L]]
  uuT <- .pv_matmul(u, .pv_T(u))
  K <- .pv_rows(list(list(0, .neg(.pv_e(u, 3L)), .pv_e(u, 2L)),
                     list(.pv_e(u, 3L), 0, .neg(.pv_e(u, 1L))),
                     list(.neg(.pv_e(u, 2L)), .pv_e(u, 1L), 0)))
  ct <- fcall("cos", list(theta))
  st <- fcall("sin", list(theta))
  pori <- .pv_madd(.pv_madd(.pv_mmap(.pv_msub(.pv_eye3(), uuT), function(a) .mul(list(a, ct))),
                            .pv_mmap(K, function(a) .mul(list(a, st)))), uuT)
  .pv_set_dcm(self, parent, pori)
  thetad <- .pv_diff(theta, .pv_t())
  wvec <- .pv_vmul(.pv_normalize(.pv_express(axis, parent)), thetad)
  .pv_set_wvec(self, parent, wvec)
}
.pvM_ReferenceFrame_orient_axis <- function(self, parent, axis, angle) invisible(.pv_orient_axis(self, parent, axis, angle))

.pv_is_matrix <- function(m) inherits(m, "symmatrix")

.pvM_ReferenceFrame_orient_explicit <- function(self, parent, dcm) {
  .pv_check_frame(parent)
  if (!.pv_is_matrix(dcm)) .pv_err("TypeError", "Amounts must be a SymPy Matrix type object.")
  .pvM_ReferenceFrame_orient_dcm(self, parent, .pv_T(dcm))
}

.pvM_ReferenceFrame_orient_dcm <- function(self, parent, dcm) {
  .pv_check_frame(parent)
  if (!.pv_is_matrix(dcm)) .pv_err("TypeError", "Amounts must be a SymPy Matrix type object.")
  .pv_set_dcm(self, parent, .pv_T(dcm))
  wvec <- .pv_w_diff_dcm(self, parent)
  .pv_set_wvec(self, parent, wvec)
}

.pv_rot <- function(axis, angle) {
  c <- fcall("cos", list(angle)); s <- fcall("sin", list(angle))
  switch(axis,
         .pv_rows(list(list(1, 0, 0), list(0, c, .neg(s)), list(0, s, c))),
         .pv_rows(list(list(c, 0, s), list(0, 1, 0), list(.neg(s), 0, c))),
         .pv_rows(list(list(c, .neg(s), 0), list(s, c, 0), list(0, 0, 1))))
}
.pvM_ReferenceFrame__rot <- function(self, axis, angle) .pv_rot(axis, as_sym(angle))

.pv_parse_rotations <- function(angles, rotation_order) {
  amounts <- .pv_seq(angles)
  for (i in seq_along(amounts)) if (!.pv_is_vec(amounts[[i]])) amounts[[i]] <- as_sym(amounts[[i]])
  rot_order <- .pv_rot_order(rotation_order)
  if (!(rot_order %in% .PV_APPROVED_ORDERS)) .pv_err("TypeError", "The rotation order is not a valid order.")
  ro <- as.integer(strsplit(rot_order, "")[[1L]])
  # Python's `len(amounts) == 3 & len(rot_order) == 3` binds as a chain
  # len(amounts) == (3 & len(rot_order)) == 3
  if (!(length(amounts) == bitwAnd(3L, length(ro)) && bitwAnd(3L, length(ro)) == 3L))
    .pv_err("TypeError", "Body orientation takes 3 values & 3 orders")
  mats <- Map(function(o, a) .pv_rot(o, a), ro, amounts)
  list(amounts = amounts, rot_order = ro, rot_matrices = mats)
}

.pv_rot_vecs <- function(amounts, ro) {
  lapply(seq_along(ro), function(i) {
    v <- rep(list(S_Zero), 3L)
    v[[ro[i]]] <- .pv_diff(amounts[[i]], .pv_t())
    .new_mat(v, 3L, 1L)
  })
}

.pvM_ReferenceFrame_orient_body_fixed <- function(self, parent, angles, rotation_order) {
  .pv_check_frame(parent)
  p <- .pv_parse_rotations(angles, rotation_order)
  R <- p$rot_matrices
  .pv_set_dcm(self, parent, .pv_matmul(.pv_matmul(R[[1L]], R[[2L]]), R[[3L]]))
  rv <- .pv_rot_vecs(p$amounts, p$rot_order)
  u <- .pv_madd(rv[[3L]], .pv_matmul(.pv_T(R[[3L]]), .pv_madd(rv[[2L]], .pv_matmul(.pv_T(R[[2L]]), rv[[1L]]))))
  wvec <- .pv_vadd(.pv_vadd(.pv_vmul(self[["_x"]], .pv_e(u, 1L)), .pv_vmul(self[["_y"]], .pv_e(u, 2L))),
                   .pv_vmul(self[["_z"]], .pv_e(u, 3L)))
  .pv_set_wvec(self, parent, wvec)
}

.pvM_ReferenceFrame_orient_space_fixed <- function(self, parent, angles, rotation_order) {
  .pv_check_frame(parent)
  p <- .pv_parse_rotations(angles, rotation_order)
  R <- p$rot_matrices
  .pv_set_dcm(self, parent, .pv_matmul(.pv_matmul(R[[3L]], R[[2L]]), R[[1L]]))
  rv <- .pv_rot_vecs(p$amounts, p$rot_order)
  u <- .pv_madd(rv[[1L]], .pv_matmul(.pv_T(R[[1L]]), .pv_madd(rv[[2L]], .pv_matmul(.pv_T(R[[2L]]), rv[[3L]]))))
  wvec <- .pv_vadd(.pv_vadd(.pv_vmul(self[["_x"]], .pv_e(u, 1L)), .pv_vmul(self[["_y"]], .pv_e(u, 2L))),
                   .pv_vmul(self[["_z"]], .pv_e(u, 3L)))
  .pv_set_wvec(self, parent, wvec)
}

.pvM_ReferenceFrame_orient_quaternion <- function(self, parent, numbers) {
  .pv_check_frame(parent)
  numbers <- .pv_pylist(numbers)
  for (i in seq_along(numbers)) if (!.pv_is_vec(numbers[[i]])) numbers[[i]] <- as_sym(numbers[[i]])
  if (length(numbers) != 4L) .pv_err("TypeError", "Amounts are a list or tuple of length 4")
  q0 <- numbers[[1L]]; q1 <- numbers[[2L]]; q2 <- numbers[[3L]]; q3 <- numbers[[4L]]
  M <- .pv_rows(list(
    list(q0^2 + q1^2 - q2^2 - q3^2, 2 * (q1 * q2 - q0 * q3), 2 * (q0 * q2 + q1 * q3)),
    list(2 * (q1 * q2 + q0 * q3), q0^2 - q1^2 + q2^2 - q3^2, 2 * (q2 * q3 - q0 * q1)),
    list(2 * (q1 * q3 - q0 * q2), 2 * (q0 * q1 + q2 * q3), q0^2 - q1^2 - q2^2 + q3^2)))
  .pv_set_dcm(self, parent, M)
  t <- .pv_t()
  q0d <- .pv_diff(q0, t); q1d <- .pv_diff(q1, t); q2d <- .pv_diff(q2, t); q3d <- .pv_diff(q3, t)
  w1 <- 2 * (q1d * q0 + q2d * q3 - q3d * q2 - q0d * q1)
  w2 <- 2 * (q2d * q0 + q3d * q1 - q1d * q3 - q0d * q2)
  w3 <- 2 * (q3d * q0 + q1d * q2 - q2d * q1 - q0d * q3)
  wvec <- .pv_vec_new(list(list(.pv_col(w1, w2, w3), self)))
  .pv_set_wvec(self, parent, wvec)
}

.pv_orient_dispatch <- function(child, parent, rot_type, amounts, rot_order) {
  rot_order <- .pv_rot_order(rot_order)
  rot_type <- toupper(rot_type)
  if (!(rot_order %in% .PV_APPROVED_ORDERS)) .pv_err("TypeError", "The supplied order is not an approved type")
  if (rot_type == "AXIS") {
    ax <- .pv_pyindex(amounts, 1L)
    .pv_orient_axis(child, parent, ax, .pv_pyindex(amounts, 0L))
  } else if (rot_type == "DCM") .pvM_ReferenceFrame_orient_explicit(child, parent, amounts)
  else if (rot_type == "BODY") .pvM_ReferenceFrame_orient_body_fixed(child, parent, amounts, rot_order)
  else if (rot_type == "SPACE") .pvM_ReferenceFrame_orient_space_fixed(child, parent, amounts, rot_order)
  else if (rot_type == "QUATERNION") .pvM_ReferenceFrame_orient_quaternion(child, parent, amounts)
  else .pv_err("NotImplementedError", "That is not an implemented rotation")
}

.pvM_ReferenceFrame_orient <- function(self, parent, rot_type, amounts, rot_order = "") {
  .pv_check_frame(parent)
  invisible(.pv_orient_dispatch(self, parent, rot_type, amounts, rot_order))
}

.pvM_ReferenceFrame_orientnew <- function(self, newname, rot_type, amounts, rot_order = "",
                                          variables = NULL, indices = NULL, latexs = NULL) {
  newframe <- ReferenceFrame(newname, variables = variables, indices = indices, latexs = latexs)
  .pv_orient_dispatch(newframe, self, rot_type, amounts, rot_order)
  newframe
}

.pvM_ReferenceFrame_set_ang_acc <- function(self, otherframe, value) {
  if (.pv_is0(value)) value <- .pv_vec0()
  value <- .pv_check_vector(value)
  .pv_check_frame(otherframe)
  self[["_ang_acc_dict"]] <- .pvod_set(self[["_ang_acc_dict"]], otherframe, value)
  otherframe[["_ang_acc_dict"]] <- .pvod_set(otherframe[["_ang_acc_dict"]], self, .pv_vneg(value))
  invisible(NULL)
}

.pvM_ReferenceFrame_set_ang_vel <- function(self, otherframe, value) {
  if (.pv_is0(value)) value <- .pv_vec0()
  value <- .pv_check_vector(value)
  .pv_check_frame(otherframe)
  self[["_ang_vel_dict"]] <- .pvod_set(self[["_ang_vel_dict"]], otherframe, value)
  otherframe[["_ang_vel_dict"]] <- .pvod_set(otherframe[["_ang_vel_dict"]], self, .pv_vneg(value))
  invisible(NULL)
}

.pvM_ReferenceFrame_partial_velocity <- function(self, frame, ...) {
  vel <- .pv_ang_vel_in(self, frame)
  partials <- .pv_partial_velocity(list(vel), list(...), frame)[[1L]]
  if (length(partials) == 1L) partials[[1L]] else .pytuple(partials)
}

# Python-visible dictionaries of a frame, for inspection
.pvP_ReferenceFrame__dcm_dict <- function(self) .od_pydict(self[["_dcm_dict"]])
.pvP_ReferenceFrame__dcm_cache <- function(self) .od_pydict(self[["_dcm_cache"]])
.pvP_ReferenceFrame__ang_vel_dict <- function(self) .od_pydict(self[["_ang_vel_dict"]])
.pvP_ReferenceFrame__ang_acc_dict <- function(self) .od_pydict(self[["_ang_acc_dict"]])
.pvP_ReferenceFrame__var_dict <- function(self) {
  vd <- self[["_var_dict"]]
  .pydict(lapply(names(vd), function(k) k), unname(vd))
}
