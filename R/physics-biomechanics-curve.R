# Musculotendon characteristic curves, after
# sympy/physics/biomechanics/curve.py.
#
# Each curve is a function of the expression tree, registered with
# def_function() so that it prints, differentiates (fdiff), evaluates
# numerically and is released by doit() as SymPy's classes are. With
# `evaluate = FALSE`, doit() wraps the inner expressions in
# UnevaluatedExpr() as SymPy does, so that the result keeps its form.

.BM_CURVE_NAMES <- c("TendonForceLengthDeGroote2016", "TendonForceLengthInverseDeGroote2016",
                     "FiberForceLengthPassiveDeGroote2016", "FiberForceLengthPassiveInverseDeGroote2016",
                     "FiberForceLengthActiveDeGroote2016", "FiberForceVelocityDeGroote2016",
                     "FiberForceVelocityInverseDeGroote2016")

.bm_ns <- environment()

.bm_ue <- function(e) UnevaluatedExpr(e)

# ArgumentIndexError of sympy/core/function.py
.bm_argindex_error <- function(e, argindex)
  .pv_err("ArgumentIndexError", sprintf("Invalid operation with argument number %d for Function %s",
                                        as.integer(argindex), .str_top(e)))

# Function.doit with the curve's own deep and evaluate handling
.bm_curve_doit <- function(formula) {
  force(formula)
  function(e, deep = TRUE, evaluate = TRUE, ...) {
    args <- e$args
    if (deep) args <- lapply(args, function(a) doit(a, deep = TRUE, evaluate = evaluate, ...))
    formula(args, evaluate)
  }
}

.bm_def_curve <- function(name, nconst, formula, fdiff, numeric, latex, inverse = NULL) {
  doit_fn <- .bm_curve_doit(formula)
  def_function(name, nargs = nconst + 1L,
               diff = function(args, i) fdiff(args, i, fnode(name, args)),
               numeric = numeric,
               latex = function(args) sprintf(latex, args[[1L]]),
               doit = doit_fn,
               code = function(e) doit_fn(e, deep = FALSE, evaluate = FALSE),
               inverse = if (!is.null(inverse)) function() get(inverse, envir = .bm_ns))
}

.bm_curve <- function(name, args) fcall(name, lapply(args, as_sym))

# curve.with_defaults(x) for a curve given by its constructor
.bm_with_defaults <- function(curve, x) {
  for (nm in .BM_CURVE_NAMES)
    if (identical(curve, get(nm, envir = .bm_ns))) return(get(paste0(nm, "_with_defaults"), envir = .bm_ns)(x))
  .pv_err("TypeError", "not a characteristic curve")
}

#' Musculotendon characteristic curves
#'
#' The characteristic curves of De Groote et al. (2016), as SymPy's
#' `sympy.physics.biomechanics` curve classes: the tendon force length
#' curve, the passive and active fiber force length curves, the fiber force
#' velocity curve and the inverses of three of them. Each call builds an
#' unevaluated function of its argument and constants; [doit()] gives the
#' formula, with `doit(f, evaluate = FALSE)` keeping its arrangement, and
#' [diff()] differentiates through it. `inverse(f)` returns the constructor
#' of the inverse curve. The `*_with_defaults()` variants use the published
#' constants, as SymPy's `with_defaults()` class methods.
#'
#' `CharacteristicCurveFunction()` is the abstract base class and raises an
#' error, as SymPy does.
#'
#' @param l_T_tilde Normalized tendon length.
#' @param fl_T Normalized tendon force.
#' @param l_M_tilde Normalized fiber length.
#' @param fl_M_pas Normalized passive fiber force.
#' @param v_M_tilde Normalized fiber velocity.
#' @param fv_M Fiber force velocity multiplier.
#' @param c0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11 Curve constants.
#' @param ... Arguments, none accepted.
#' @return A `sym` expression.
#' @examples
#' l <- Symbol("l_T_tilde")
#' f <- TendonForceLengthDeGroote2016_with_defaults(l)
#' f
#' doit(f)
#' diff(f, l)
#' inverse(f)
#' @name characteristic-curves
NULL

#' @rdname characteristic-curves
#' @export
CharacteristicCurveFunction <- function(...) {
  n <- length(list(...))
  if (n) .pv_err("TypeError", sprintf("CharacteristicCurveFunction takes exactly 0 arguments (%d given)", n))
  .pv_err("TypeError", paste0("Cannot directly instantiate 'CharacteristicCurveFunction', instances of ",
                              "characteristic curves must be of a concrete subclass."))
}

# ---------------------------------------------------------------------------
# TendonForceLengthDeGroote2016

#' @rdname characteristic-curves
#' @export
TendonForceLengthDeGroote2016 <- function(l_T_tilde, c0, c1, c2, c3)
  .bm_curve("TendonForceLengthDeGroote2016", list(l_T_tilde, c0, c1, c2, c3))

#' @rdname characteristic-curves
#' @export
TendonForceLengthDeGroote2016_with_defaults <- function(l_T_tilde)
  TendonForceLengthDeGroote2016(l_T_tilde, Float("0.2"), Float("0.995"), Float("0.25"), Float("33.93669377311689"))

.bm_def_curve(
  "TendonForceLengthDeGroote2016", 4L,
  formula = function(a, evaluate) {
    l_T_tilde <- a[[1L]]; c0 <- a[[2L]]; c1 <- a[[3L]]; c2 <- a[[4L]]; c3 <- a[[5L]]
    if (evaluate) return(c0 * exp(c3 * (l_T_tilde - c1)) - c2)
    c0 * exp(c3 * .bm_ue(l_T_tilde - c1)) - c2
  },
  fdiff = function(a, argindex, self) {
    l_T_tilde <- a[[1L]]; c0 <- a[[2L]]; c1 <- a[[3L]]; c2 <- a[[4L]]; c3 <- a[[5L]]
    switch(as.character(argindex),
           `1` = c0 * c3 * exp(c3 * .bm_ue(l_T_tilde - c1)),
           `2` = exp(c3 * .bm_ue(l_T_tilde - c1)),
           `3` = -c0 * c3 * exp(c3 * .bm_ue(l_T_tilde - c1)),
           `4` = as_sym(-1L),
           `5` = c0 * (l_T_tilde - c1) * exp(c3 * .bm_ue(l_T_tilde - c1)),
           .bm_argindex_error(self, argindex))
  },
  numeric = function(l_T_tilde, c0, c1, c2, c3) c0 * exp(c3 * (l_T_tilde - c1)) - c2,
  latex = "\\operatorname{fl}^T \\left( %s \\right)",
  inverse = "TendonForceLengthInverseDeGroote2016")

# ---------------------------------------------------------------------------
# TendonForceLengthInverseDeGroote2016

#' @rdname characteristic-curves
#' @export
TendonForceLengthInverseDeGroote2016 <- function(fl_T, c0, c1, c2, c3)
  .bm_curve("TendonForceLengthInverseDeGroote2016", list(fl_T, c0, c1, c2, c3))

#' @rdname characteristic-curves
#' @export
TendonForceLengthInverseDeGroote2016_with_defaults <- function(fl_T)
  TendonForceLengthInverseDeGroote2016(fl_T, Float("0.2"), Float("0.995"), Float("0.25"), Float("33.93669377311689"))

.bm_def_curve(
  "TendonForceLengthInverseDeGroote2016", 4L,
  formula = function(a, evaluate) {
    fl_T <- a[[1L]]; c0 <- a[[2L]]; c1 <- a[[3L]]; c2 <- a[[4L]]; c3 <- a[[5L]]
    if (evaluate) return(log((fl_T + c2) / c0) / c3 + c1)
    log(.bm_ue((fl_T + c2) / c0)) / c3 + c1
  },
  fdiff = function(a, argindex, self) {
    fl_T <- a[[1L]]; c0 <- a[[2L]]; c1 <- a[[3L]]; c2 <- a[[4L]]; c3 <- a[[5L]]
    switch(as.character(argindex),
           `1` = 1 / (c3 * (fl_T + c2)),
           `2` = -1 / (c0 * c3),
           `3` = S_One,
           `4` = 1 / (c3 * (fl_T + c2)),
           `5` = -log(.bm_ue((fl_T + c2) / c0)) / c3^2,
           .bm_argindex_error(self, argindex))
  },
  numeric = function(fl_T, c0, c1, c2, c3) log((fl_T + c2) / c0) / c3 + c1,
  latex = "\\left( \\operatorname{fl}^T \\right)^{-1} \\left( %s \\right)",
  inverse = "TendonForceLengthDeGroote2016")

# ---------------------------------------------------------------------------
# FiberForceLengthPassiveDeGroote2016

#' @rdname characteristic-curves
#' @export
FiberForceLengthPassiveDeGroote2016 <- function(l_M_tilde, c0, c1)
  .bm_curve("FiberForceLengthPassiveDeGroote2016", list(l_M_tilde, c0, c1))

#' @rdname characteristic-curves
#' @export
FiberForceLengthPassiveDeGroote2016_with_defaults <- function(l_M_tilde)
  FiberForceLengthPassiveDeGroote2016(l_M_tilde, Float("0.6"), Float("4.0"))

.bm_def_curve(
  "FiberForceLengthPassiveDeGroote2016", 2L,
  formula = function(a, evaluate) {
    l_M_tilde <- a[[1L]]; c0 <- a[[2L]]; c1 <- a[[3L]]
    if (evaluate) return((exp(c1 * (l_M_tilde - 1) / c0) - 1) / (exp(c1) - 1))
    (exp(c1 * .bm_ue(l_M_tilde - 1) / c0) - 1) / (exp(c1) - 1)
  },
  fdiff = function(a, argindex, self) {
    l_M_tilde <- a[[1L]]; c0 <- a[[2L]]; c1 <- a[[3L]]
    switch(as.character(argindex),
           `1` = c1 * exp(c1 * .bm_ue(l_M_tilde - 1) / c0) / (c0 * (exp(c1) - 1)),
           `2` = -c1 * exp(c1 * .bm_ue(l_M_tilde - 1) / c0) * .bm_ue(l_M_tilde - 1) / (c0^2 * (exp(c1) - 1)),
           `3` = -exp(c1) * (-1 + exp(c1 * .bm_ue(l_M_tilde - 1) / c0)) / (exp(c1) - 1)^2 +
             exp(c1 * .bm_ue(l_M_tilde - 1) / c0) * (l_M_tilde - 1) / (c0 * (exp(c1) - 1)),
           .bm_argindex_error(self, argindex))
  },
  numeric = function(l_M_tilde, c0, c1) (exp(c1 * (l_M_tilde - 1) / c0) - 1) / (exp(c1) - 1),
  latex = "\\operatorname{fl}^M_{pas} \\left( %s \\right)",
  inverse = "FiberForceLengthPassiveInverseDeGroote2016")

# ---------------------------------------------------------------------------
# FiberForceLengthPassiveInverseDeGroote2016

#' @rdname characteristic-curves
#' @export
FiberForceLengthPassiveInverseDeGroote2016 <- function(fl_M_pas, c0, c1)
  .bm_curve("FiberForceLengthPassiveInverseDeGroote2016", list(fl_M_pas, c0, c1))

#' @rdname characteristic-curves
#' @export
FiberForceLengthPassiveInverseDeGroote2016_with_defaults <- function(fl_M_pas)
  FiberForceLengthPassiveInverseDeGroote2016(fl_M_pas, Float("0.6"), Float("4.0"))

.bm_def_curve(
  "FiberForceLengthPassiveInverseDeGroote2016", 2L,
  formula = function(a, evaluate) {
    fl_M_pas <- a[[1L]]; c0 <- a[[2L]]; c1 <- a[[3L]]
    if (evaluate) return(c0 * log(fl_M_pas * (exp(c1) - 1) + 1) / c1 + 1)
    c0 * log(.bm_ue(fl_M_pas * (exp(c1) - 1)) + 1) / c1 + 1
  },
  fdiff = function(a, argindex, self) {
    fl_M_pas <- a[[1L]]; c0 <- a[[2L]]; c1 <- a[[3L]]
    switch(as.character(argindex),
           `1` = c0 * (exp(c1) - 1) / (c1 * (fl_M_pas * (exp(c1) - 1) + 1)),
           `2` = log(fl_M_pas * (exp(c1) - 1) + 1) / c1,
           `3` = c0 * fl_M_pas * exp(c1) / (c1 * (fl_M_pas * (exp(c1) - 1) + 1)) -
             c0 * log(fl_M_pas * (exp(c1) - 1) + 1) / c1^2,
           .bm_argindex_error(self, argindex))
  },
  numeric = function(fl_M_pas, c0, c1) c0 * log(fl_M_pas * (exp(c1) - 1) + 1) / c1 + 1,
  latex = "\\left( \\operatorname{fl}^M_{pas} \\right)^{-1} \\left( %s \\right)",
  inverse = "FiberForceLengthPassiveDeGroote2016")

# ---------------------------------------------------------------------------
# FiberForceLengthActiveDeGroote2016

#' @rdname characteristic-curves
#' @export
FiberForceLengthActiveDeGroote2016 <- function(l_M_tilde, c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11)
  .bm_curve("FiberForceLengthActiveDeGroote2016", list(l_M_tilde, c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11))

#' @rdname characteristic-curves
#' @export
FiberForceLengthActiveDeGroote2016_with_defaults <- function(l_M_tilde)
  FiberForceLengthActiveDeGroote2016(l_M_tilde, Float("0.814"), Float("1.06"), Float("0.162"), Float("0.0633"),
                                     Float("0.433"), Float("0.717"), Float("-0.0299"), Float("0.2"), Float("0.1"),
                                     Float("1.0"), Float("0.354"), Float("0.0"))

# one Gaussian term c_a*exp(-((l - c_b)/(c_c + c_d*l))^2/2), its derivative
# with respect to l and its partials with respect to the four constants
.bm_gauss <- function(l, ca, cb, cc, cd, ue = TRUE) {
  num <- if (ue) .bm_ue(l - cb) else l - cb
  ca * exp(-(num / (cc + cd * l))^2 / 2)
}
.bm_gauss_e <- function(l, cb, cc, cd) exp(-(l - cb)^2 / (2 * (cc + cd * l)^2))

.bm_def_curve(
  "FiberForceLengthActiveDeGroote2016", 12L,
  formula = function(a, evaluate) {
    l <- a[[1L]]; k <- a[-1L]
    .bm_gauss(l, k[[1L]], k[[2L]], k[[3L]], k[[4L]], !evaluate) +
      .bm_gauss(l, k[[5L]], k[[6L]], k[[7L]], k[[8L]], !evaluate) +
      .bm_gauss(l, k[[9L]], k[[10L]], k[[11L]], k[[12L]], !evaluate)
  },
  fdiff = function(a, argindex, self) {
    l <- a[[1L]]; k <- a[-1L]
    dl <- function(ca, cb, cc, cd)
      ca * (cd * (l - cb)^2 / (cc + cd * l)^3 + (cb - l) / (cc + cd * l)^2) * .bm_gauss_e(l, cb, cc, cd)
    if (argindex == 1L)
      return(dl(k[[1L]], k[[2L]], k[[3L]], k[[4L]]) + dl(k[[5L]], k[[6L]], k[[7L]], k[[8L]]) +
               dl(k[[9L]], k[[10L]], k[[11L]], k[[12L]]))
    if (argindex < 1L || argindex > 13L) .bm_argindex_error(self, argindex)
    g <- (argindex - 2L) %/% 4L
    j <- (argindex - 2L) %% 4L
    ca <- k[[4L * g + 1L]]; cb <- k[[4L * g + 2L]]; cc <- k[[4L * g + 3L]]; cd <- k[[4L * g + 4L]]
    switch(j + 1L,
           .bm_gauss_e(l, cb, cc, cd),
           ca * (l - cb) / (cc + cd * l)^2 * .bm_gauss_e(l, cb, cc, cd),
           ca * (l - cb)^2 / (cc + cd * l)^3 * .bm_gauss_e(l, cb, cc, cd),
           ca * l * (l - cb)^2 / (cc + cd * l)^3 * .bm_gauss_e(l, cb, cc, cd))
  },
  numeric = function(l, c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11) {
    g <- function(ca, cb, cc, cd) ca * exp(-((l - cb) / (cc + cd * l))^2 / 2)
    g(c0, c1, c2, c3) + g(c4, c5, c6, c7) + g(c8, c9, c10, c11)
  },
  latex = "\\operatorname{fl}^M_{act} \\left( %s \\right)")

# ---------------------------------------------------------------------------
# FiberForceVelocityDeGroote2016

#' @rdname characteristic-curves
#' @export
FiberForceVelocityDeGroote2016 <- function(v_M_tilde, c0, c1, c2, c3)
  .bm_curve("FiberForceVelocityDeGroote2016", list(v_M_tilde, c0, c1, c2, c3))

#' @rdname characteristic-curves
#' @export
FiberForceVelocityDeGroote2016_with_defaults <- function(v_M_tilde)
  FiberForceVelocityDeGroote2016(v_M_tilde, Float("-0.318"), Float("-8.149"), Float("-0.374"), Float("0.886"))

.bm_def_curve(
  "FiberForceVelocityDeGroote2016", 4L,
  formula = function(a, evaluate) {
    v_M_tilde <- a[[1L]]; c0 <- a[[2L]]; c1 <- a[[3L]]; c2 <- a[[4L]]; c3 <- a[[5L]]
    if (evaluate) return(c0 * log(c1 * v_M_tilde + c2 + sqrt((c1 * v_M_tilde + c2)^2 + 1)) + c3)
    c0 * log(c1 * v_M_tilde + c2 + sqrt(.bm_ue(c1 * v_M_tilde + c2)^2 + 1)) + c3
  },
  fdiff = function(a, argindex, self) {
    v_M_tilde <- a[[1L]]; c0 <- a[[2L]]; c1 <- a[[3L]]; c2 <- a[[4L]]; c3 <- a[[5L]]
    root <- function() sqrt(.bm_ue(c1 * v_M_tilde + c2)^2 + 1)
    switch(as.character(argindex),
           `1` = c0 * c1 / root(),
           `2` = log(c1 * v_M_tilde + c2 + root()),
           `3` = c0 * v_M_tilde / root(),
           `4` = c0 / root(),
           `5` = S_One,
           .bm_argindex_error(self, argindex))
  },
  numeric = function(v_M_tilde, c0, c1, c2, c3) c0 * log(c1 * v_M_tilde + c2 + sqrt((c1 * v_M_tilde + c2)^2 + 1)) + c3,
  latex = "\\operatorname{fv}^M \\left( %s \\right)",
  inverse = "FiberForceVelocityInverseDeGroote2016")

# ---------------------------------------------------------------------------
# FiberForceVelocityInverseDeGroote2016

#' @rdname characteristic-curves
#' @export
FiberForceVelocityInverseDeGroote2016 <- function(fv_M, c0, c1, c2, c3)
  .bm_curve("FiberForceVelocityInverseDeGroote2016", list(fv_M, c0, c1, c2, c3))

#' @rdname characteristic-curves
#' @export
FiberForceVelocityInverseDeGroote2016_with_defaults <- function(fv_M)
  FiberForceVelocityInverseDeGroote2016(fv_M, Float("-0.318"), Float("-8.149"), Float("-0.374"), Float("0.886"))

.bm_def_curve(
  "FiberForceVelocityInverseDeGroote2016", 4L,
  formula = function(a, evaluate) {
    fv_M <- a[[1L]]; c0 <- a[[2L]]; c1 <- a[[3L]]; c2 <- a[[4L]]; c3 <- a[[5L]]
    if (evaluate) return((sinh((fv_M - c3) / c0) - c2) / c1)
    (sinh(.bm_ue(fv_M - c3) / c0) - c2) / c1
  },
  fdiff = function(a, argindex, self) {
    fv_M <- a[[1L]]; c0 <- a[[2L]]; c1 <- a[[3L]]; c2 <- a[[4L]]; c3 <- a[[5L]]
    switch(as.character(argindex),
           `1` = cosh((fv_M - c3) / c0) / (c0 * c1),
           `2` = (c3 - fv_M) * cosh((fv_M - c3) / c0) / (c0^2 * c1),
           `3` = (c2 - sinh((fv_M - c3) / c0)) / c1^2,
           `4` = -1 / c1,
           `5` = -cosh((fv_M - c3) / c0) / (c0 * c1),
           .bm_argindex_error(self, argindex))
  },
  numeric = function(fv_M, c0, c1, c2, c3) (sinh((fv_M - c3) / c0) - c2) / c1,
  latex = "\\left( \\operatorname{fv}^M \\right)^{-1} \\left( %s \\right)",
  inverse = "FiberForceVelocityDeGroote2016")

# ---------------------------------------------------------------------------
# CharacteristicCurveCollection

.BM_COLLECTION_FIELDS <- c("tendon_force_length", "tendon_force_length_inverse", "fiber_force_length_passive",
                           "fiber_force_length_passive_inverse", "fiber_force_length_active",
                           "fiber_force_velocity", "fiber_force_velocity_inverse")

#' Collection of characteristic curves
#'
#' A frozen record of the seven characteristic curve constructors a
#' musculotendon model uses, as SymPy's `CharacteristicCurveCollection`
#' dataclass. Fields are read with `$` and cannot be reassigned.
#'
#' @param tendon_force_length,tendon_force_length_inverse,fiber_force_length_passive,fiber_force_length_passive_inverse,fiber_force_length_active,fiber_force_velocity,fiber_force_velocity_inverse
#'   Curve constructors, such as [TendonForceLengthDeGroote2016()].
#' @return A `CharacteristicCurveCollection` object.
#' @examples
#' curves <- CharacteristicCurveCollection(
#'   TendonForceLengthDeGroote2016, TendonForceLengthInverseDeGroote2016,
#'   FiberForceLengthPassiveDeGroote2016, FiberForceLengthPassiveInverseDeGroote2016,
#'   FiberForceLengthActiveDeGroote2016, FiberForceVelocityDeGroote2016,
#'   FiberForceVelocityInverseDeGroote2016)
#' curves$tendon_force_length(Symbol("l"), 1, 2, 3, 4)
#' @export
CharacteristicCurveCollection <- function(tendon_force_length, tendon_force_length_inverse, fiber_force_length_passive,
                                          fiber_force_length_passive_inverse, fiber_force_length_active,
                                          fiber_force_velocity, fiber_force_velocity_inverse) {
  fields <- list(tendon_force_length, tendon_force_length_inverse, fiber_force_length_passive,
                 fiber_force_length_passive_inverse, fiber_force_length_active,
                 fiber_force_velocity, fiber_force_velocity_inverse)
  names(fields) <- .BM_COLLECTION_FIELDS
  structure(fields, class = c("CharacteristicCurveCollection", "physics_obj"))
}

.bm_curve_label <- function(f) {
  for (nm in .BM_CURVE_NAMES) if (identical(f, get(nm, envir = .bm_ns))) return(nm)
  .pv_objrepr(f)
}

.pvR_CharacteristicCurveCollection <- function(x)
  sprintf("CharacteristicCurveCollection(%s)",
          paste(sprintf("%s=%s", .BM_COLLECTION_FIELDS, vapply(unclass(x), .bm_curve_label, "")), collapse = ", "))

#' @export
format.CharacteristicCurveCollection <- function(x, ...) .pv_objrepr(x)

for (.f in .BM_COLLECTION_FIELDS) {
  assign(paste0(".pvS_CharacteristicCurveCollection_", .f), local({
    fld <- .f
    function(self, value) .pv_err("FrozenInstanceError", sprintf("cannot assign to field '%s'", fld))
  }))
}
rm(.f)

.pvM_CharacteristicCurveCollection___iter__ <- function(self) unname(unclass(self))

.bm_degroote_curves <- function()
  CharacteristicCurveCollection(TendonForceLengthDeGroote2016, TendonForceLengthInverseDeGroote2016,
                                FiberForceLengthPassiveDeGroote2016, FiberForceLengthPassiveInverseDeGroote2016,
                                FiberForceLengthActiveDeGroote2016, FiberForceVelocityDeGroote2016,
                                FiberForceVelocityInverseDeGroote2016)
