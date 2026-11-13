# Gaussian optics, after sympy/physics/optics/gaussopt.py.
#
# The ray transfer matrices and the geometric ray are explicit matrices whose
# class also names the optical element; `%*%` follows their __mul__, so a
# matrix times a ray is a ray and a matrix times a beam parameter is the
# transformed beam parameter. A BeamParameter is an immutable object.

.optics_mat <- function(M, cls) {
  class(M) <- c(cls, if (!identical(cls, "RayTransferMatrix") && !identical(cls, "GeometricRay")) "RayTransferMatrix", "symmatrix")
  attr(M, "pyclass") <- cls[1L]
  M
}

.optics_args_str <- function(args) {
  s <- vapply(args, function(a) if (inherits(a, "symmatrix")) sstr(a) else .str(as_sym(a)), "")
  paste0("(", paste(s, collapse = ", "), if (length(args) == 1L) "," else "", ")")
}

#' Ray transfer matrices of Gaussian optics
#'
#' `RayTransferMatrix(A, B, C, D)` (or a 2 by 2 matrix) is the ABCD matrix of
#' an optical element; `FreeSpace(d)`, `FlatRefraction(n1, n2)`,
#' `CurvedRefraction(R, n1, n2)`, `FlatMirror()`, `CurvedMirror(R)` and
#' `ThinLens(f)` are the matrices of common elements. `GeometricRay(h, angle)`
#' (or a 2 by 1 matrix) is a ray. The product `%*%` of a matrix with another
#' matrix, a ray or a [BeamParameter()] gives the transformed object.
#' `M$A`, `M$B`, `M$C` and `M$D` are the entries; `ray$height` and
#' `ray$angle` the parts of a ray. After `sympy.physics.optics.gaussopt`.
#'
#' @param ... The four entries, or one 2 by 2 matrix.
#' @param d A distance.
#' @param n1,n2 The refractive indices before and after the interface.
#' @param R The radius of curvature.
#' @param f The focal distance.
#' @return A `symmatrix` of class `RayTransferMatrix` (or `GeometricRay`).
#' @examples
#' d <- Symbol("d"); f <- Symbol("f")
#' ThinLens(f) %*% FreeSpace(d)
#' FreeSpace(d) %*% GeometricRay(Symbol("h"), Symbol("angle"))
#' @export
RayTransferMatrix <- function(...) .rtm_new(list(...), "RayTransferMatrix")

.rtm_new <- function(args, cls) {
  if (length(args) == 4L) {
    temp <- Matrix(list(list(args[[1L]], args[[2L]]), list(args[[3L]], args[[4L]])))
  } else if (length(args) == 1L && inherits(args[[1L]], "symmatrix") && args[[1L]]$nrow == 2L && args[[1L]]$ncol == 2L) {
    temp <- .new_mat(args[[1L]]$data, 2L, 2L)
  } else {
    .physics_error("ValueError", filldedent(paste0("\n    Expecting 2x2 Matrix or the 4 elements of\n    the Matrix but got ", .optics_args_str(args))))
  }
  .optics_mat(temp, cls)
}

#' @rdname RayTransferMatrix
#' @export
FreeSpace <- function(d) .rtm_new(list(1L, d, 0L, 1L), "FreeSpace")

#' @rdname RayTransferMatrix
#' @export
FlatRefraction <- function(n1, n2) {
  n1 <- as_sym(n1); n2 <- as_sym(n2)
  .rtm_new(list(1L, 0L, 0L, n1 / n2), "FlatRefraction")
}

#' @rdname RayTransferMatrix
#' @export
CurvedRefraction <- function(R, n1, n2) {
  R <- as_sym(R); n1 <- as_sym(n1); n2 <- as_sym(n2)
  .rtm_new(list(1L, 0L, (n1 - n2) / R / n2, n1 / n2), "CurvedRefraction")
}

#' @rdname RayTransferMatrix
#' @export
FlatMirror <- function() .rtm_new(list(1L, 0L, 0L, 1L), "FlatMirror")

#' @rdname RayTransferMatrix
#' @export
CurvedMirror <- function(R) {
  R <- as_sym(R)
  .rtm_new(list(1L, 0L, -2 / R, 1L), "CurvedMirror")
}

#' @rdname RayTransferMatrix
#' @export
ThinLens <- function(f) {
  f <- as_sym(f)
  .rtm_new(list(1L, 0L, -1 / f, 1L), "ThinLens")
}

#' @rdname RayTransferMatrix
#' @export
GeometricRay <- function(...) {
  args <- list(...)
  if (length(args) == 1L && inherits(args[[1L]], "symmatrix") && args[[1L]]$nrow == 2L && args[[1L]]$ncol == 1L) {
    temp <- .new_mat(args[[1L]]$data, 2L, 1L)
  } else if (length(args) == 2L) {
    temp <- Matrix(list(list(args[[1L]]), list(args[[2L]])))
  } else {
    .physics_error("ValueError", filldedent(paste0("\n    Expecting 2x1 Matrix or the 2 elements of\n    the Matrix but got ", .optics_args_str(args))))
  }
  .optics_mat(temp, "GeometricRay")
}

.plain_mat <- function(M) .new_mat(M$data, M$nrow, M$ncol)

#' @export
`$.RayTransferMatrix` <- function(x, name) {
  switch(name,
         A = .mget(x, 1L, 1L), B = .mget(x, 1L, 2L), C = .mget(x, 2L, 1L), D = .mget(x, 2L, 2L),
         .subset2(x, name))
}

#' @export
`$.GeometricRay` <- function(x, name) {
  switch(name, height = .subset2(x, "data")[[1L]], angle = .subset2(x, "data")[[2L]], .subset2(x, name))
}

#' @export
matrixOps.RayTransferMatrix <- function(x, y) {
  if (.Generic == "%*%" && inherits(x, "RayTransferMatrix")) {
    if (inherits(y, "RayTransferMatrix")) return(RayTransferMatrix(.matmul(.plain_mat(x), .plain_mat(y))))
    if (inherits(y, "GeometricRay")) return(GeometricRay(.matmul(.plain_mat(x), .plain_mat(y))))
    if (inherits(y, "BeamParameter")) {
      temp <- .matmul(.plain_mat(x), Matrix(list(list(y$q), list(S_One))))
      q <- expand(.div(temp$data[[1L]], temp$data[[2L]]), complex = TRUE)
      return(BeamParameter(y$wavelen, together(re(q)), z_r = together(im(q))))
    }
    return(matrixOps.symmatrix(.plain_mat(x), if (inherits(y, "symmatrix")) .plain_mat(y) else y))
  }
  matrixOps.symmatrix(.plain_mat(x), if (inherits(y, "symmatrix")) .plain_mat(y) else y)
}

#' The parameter of a Gaussian beam
#'
#' `BeamParameter(wavelen, z, z_r = NULL, w = NULL, n = 1)` is the complex
#' beam parameter `q = z + I*z_r` of a Gaussian beam of wavelength `wavelen`
#' at distance `z` from the waist, given by the Rayleigh range `z_r` or the
#' waist `w`, in a medium of refractive index `n`. Its properties are
#' `$wavelen`, `$z`, `$z_r`, `$n`, `$q`, `$radius`, `$w`, `$w_0`,
#' `$divergence`, `$gouy` and `$waist_approximation_limit`.
#' `waist2rayleigh()` and `rayleigh2waist()` convert between the waist and
#' the Rayleigh range. After `sympy.physics.optics.gaussopt`.
#'
#' @param wavelen The wavelength.
#' @param z The distance to the waist.
#' @param z_r The Rayleigh range.
#' @param w The waist.
#' @param n The refractive index of the medium.
#' @return A `BeamParameter` object, or a `sym` expression.
#' @examples
#' p <- BeamParameter(530e-9, 1, w = 1e-3)
#' p$q
#' N(p$w_0)
#' waist2rayleigh(Symbol("w"), Symbol("wavelen"))
#' @export
BeamParameter <- function(wavelen, z, z_r = NULL, w = NULL, n = 1L) {
  wavelen <- as_sym(wavelen); z <- as_sym(z); n <- as_sym(n)
  if (!is.null(z_r) && is.null(w)) {
    z_r <- as_sym(z_r)
  } else if (!is.null(w) && is.null(z_r)) {
    z_r <- waist2rayleigh(as_sym(w), wavelen, n)
  } else if (is.null(z_r) && is.null(w)) {
    .physics_error("ValueError", "Must specify one of w and z_r.")
  }
  .pyobj_new("BeamParameter", list(wavelen, z, z_r, n))
}

local({
  bp <- function(i) function(self) self$args[[i]]
  .pyclass_def("BeamParameter", expr = TRUE,
    ctor = function(wavelen, z, z_r, n) BeamParameter(wavelen, z, z_r = z_r, n = n),
    props = list(
      wavelen = bp(1L), z = bp(2L), z_r = bp(3L), n = bp(4L),
      q = function(self) self$z + S_I * self$z_r,
      radius = function(self) self$z * (1 + (self$z_r / self$z)^2),
      w = function(self) self$w_0 * sqrt(1 + (self$z / self$z_r)^2),
      w_0 = function(self) sqrt(self$z_r / (S_Pi * self$n) * self$wavelen),
      divergence = function(self) self$wavelen / S_Pi / self$w_0,
      gouy = function(self) sym_atan2(self$z, self$z_r),
      waist_approximation_limit = function(self) 2 * self$wavelen / S_Pi))
})

#' @rdname BeamParameter
#' @export
waist2rayleigh <- function(w, wavelen, n = 1L) {
  w <- as_sym(w); wavelen <- as_sym(wavelen)
  w^2 * n * S_Pi / wavelen
}

#' @rdname BeamParameter
#' @export
rayleigh2waist <- function(z_r, wavelen) {
  z_r <- as_sym(z_r); wavelen <- as_sym(wavelen)
  sqrt(z_r / S_Pi * wavelen)
}

#' Conjugation relations of geometric and Gaussian optics
#'
#' `geometric_conj_ab(a, b)` is the focal distance that images a point at
#' distance `a` to distance `b`; `geometric_conj_af(a, f)` and
#' `geometric_conj_bf(b, f)` give the other distance from one distance and
#' the focal distance. `gaussian_conj(s_in, z_r_in, f)` gives the list
#' `(s_out, z_r_out, m)` for a Gaussian beam through a thin lens, and
#' `conjugate_gauss_beams(wavelen, waist_in, waist_out, f = f)` the list
#' `(s_in, s_out, f)` of the optical setup that conjugates two beams. After
#' `sympy.physics.optics.gaussopt`.
#'
#' @param a,b The distances of the object and the image.
#' @param f The focal distance.
#' @param s_in The distance of the waist before the lens.
#' @param z_r_in The Rayleigh range of the incoming beam.
#' @param wavelen The wavelength.
#' @param waist_in,waist_out The waists of the two beams.
#' @param ... The focal distance `f`, the only setup supported.
#' @return A `sym` expression or a list.
#' @examples
#' a <- Symbol("a"); b <- Symbol("b"); f <- Symbol("f")
#' geometric_conj_ab(a, b)
#' gaussian_conj(Symbol("s_in"), Symbol("z_r_in"), f)[[1]]
#' @export
geometric_conj_ab <- function(a, b) {
  a <- as_sym(a); b <- as_sym(b)
  if (isTRUE(ask_fact(a, "infinite")) || isTRUE(ask_fact(b, "infinite")))
    return(if (isTRUE(ask_fact(b, "infinite"))) a else b)
  a * b / (a + b)
}

#' @rdname geometric_conj_ab
#' @export
geometric_conj_af <- function(a, f) {
  a <- as_sym(a); f <- as_sym(f)
  -geometric_conj_ab(a, -f)
}

#' @rdname geometric_conj_ab
#' @export
geometric_conj_bf <- function(a, f) geometric_conj_af(a, f)

#' @rdname geometric_conj_ab
#' @export
gaussian_conj <- function(s_in, z_r_in, f) {
  s_in <- as_sym(s_in); z_r_in <- as_sym(z_r_in); f <- as_sym(f)
  s_out <- 1 / (-1 / (s_in + z_r_in^2 / (s_in - f)) + 1 / f)
  m <- 1 / sqrt((1 - (s_in / f)^2) + (z_r_in / f)^2)
  z_r_out <- z_r_in / ((1 - (s_in / f)^2) + (z_r_in / f)^2)
  .pytuple(list(s_out, z_r_out, m))
}

#' @rdname geometric_conj_ab
#' @export
conjugate_gauss_beams <- function(wavelen, waist_in, waist_out, ...) {
  kwargs <- list(...)
  wavelen <- as_sym(wavelen); waist_in <- as_sym(waist_in); waist_out <- as_sym(waist_out)
  m <- waist_out / waist_in
  z <- waist2rayleigh(waist_in, wavelen)
  if (length(kwargs) != 1L) .physics_error("ValueError", "The function expects only one named argument")
  if ("dist" %in% names(kwargs))
    .physics_error("NotImplementedError", filldedent("\n    Currently only focal length is supported as a parameter"))
  if ("f" %in% names(kwargs)) {
    f <- as_sym(kwargs$f)
    s_in <- f * (1 - sqrt(1 / m^2 - z^2 / f^2))
    s_out <- gaussian_conj(s_in, z, f)[[1L]]
  } else if ("s_in" %in% names(kwargs)) {
    .physics_error("NotImplementedError", filldedent("\n    Currently only focal length is supported as a parameter"))
  } else {
    .physics_error("ValueError", filldedent("\n    The functions expects the focal length as a named argument"))
  }
  .pytuple(list(s_in, s_out, f))
}
