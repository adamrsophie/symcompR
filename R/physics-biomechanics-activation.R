# Activation dynamics, after sympy/physics/biomechanics/activation.py and
# the name handling of sympy/physics/biomechanics/_mixin.py.

# _NamedMixin.name setter: a nonempty string, set once
.bm_set_name <- function(self, name) {
  if (.pm_has(self, "_name"))
    .pv_err("AttributeError", sprintf("Can't set attribute `name` to %s as it is immutable.", .pv_objrepr(name)))
  if (!(is.character(name) && length(name) == 1L))
    .pv_err("TypeError", sprintf("Name %s passed to `name` was of type %s, must be <class 'str'>.",
                                 .pv_objrepr(name), .pv_pyclass(name)))
  if (!nzchar(name))
    .pv_err("ValueError", sprintf("Name %s is invalid, must be a nonzero length <class 'type'>.", .pv_objrepr(name)))
  self[["_name"]] <- name
}

.bm_named_classes <- c("ActivationBase", "MusculotendonBase")
for (.c in .bm_named_classes) {
  assign(paste0(".pvP_", .c, "_name"), function(self) self[["_name"]])
  assign(paste0(".pvS_", .c, "_name"), function(self, value) .bm_set_name(self, value))
}
rm(.c)

.bm_abstract <- function(cls, methods)
  .pv_err("TypeError", sprintf("Can't instantiate abstract class %s without an implementation for abstract method%s %s",
                               cls, if (length(methods) > 1L) "s" else "", paste0("'", methods, "'", collapse = ", ")))

.bm_colmat <- function(xs) if (length(xs)) .pv_colv(xs) else zeros(0L, 1L)

#' Activation dynamics
#'
#' Models of how a muscle's activation follows its neural excitation, as
#' SymPy's `sympy.physics.biomechanics.activation` classes.
#' `ZerothOrderActivation()` makes the activation equal to the excitation;
#' `FirstOrderActivationDeGroote2016()` uses the first order differential
#' equation of De Groote et al. (2016), with time constants and a smoothing
#' rate that default to symbols named after the model. The
#' `*_with_defaults()` variants use the published constants.
#' `ActivationBase()` is the abstract base class and raises an error, as
#' SymPy does.
#'
#' The objects expose SymPy's attributes through `$`: `name`, `excitation`
#' (`e`), `activation` (`a`), `order`, `state_vars` (`x`), `input_vars`
#' (`r`), `constants` (`p`), `M`, `F` and the method `rhs()`.
#'
#' @param name Name of the model, used in the symbols it creates.
#' @param activation_time_constant,deactivation_time_constant,smoothing_rate
#'   Constants of the first order model; `NULL` gives symbols.
#' @return An activation dynamics object.
#' @examples
#' act <- FirstOrderActivationDeGroote2016("bicep")
#' act$x
#' act$rhs()
#' ZerothOrderActivation("bicep")$r
#' @name activation-dynamics
NULL

#' @rdname activation-dynamics
#' @export
ActivationBase <- function(name)
  .bm_abstract("ActivationBase", c("F", "M", "constants", "input_vars", "order", "p", "r", "rhs",
                                   "state_vars", "with_defaults", "x"))

.bm_activation_init <- function(self, name) {
  .bm_set_name(self, .pv_str(name))
  self[["_e"]] <- dynamicsymbols(paste0("e_", .pv_str(name)))
  self[["_a"]] <- dynamicsymbols(paste0("a_", .pv_str(name)))
}

.pvP_ActivationBase_excitation <- function(self) self[["_e"]]
.pvP_ActivationBase_e <- function(self) self[["_e"]]
.pvP_ActivationBase_activation <- function(self) self[["_a"]]
.pvP_ActivationBase_a <- function(self) self[["_a"]]
.pvM_ActivationBase___eq__ <- function(self, other) {
  if (!identical(class(self)[1L], class(other)[1L])) return(FALSE)
  identical(self$name, other$name)
}
.pvR_ActivationBase <- function(x) sprintf("%s(%s)", class(x)[1L], .pv_objrepr(x$name))

#' @export
format.ActivationBase <- function(x, ...) .pv_objrepr(x)

# ---------------------------------------------------------------------------
# ZerothOrderActivation

#' @rdname activation-dynamics
#' @export
ZerothOrderActivation <- function(name) {
  self <- .pm_new_obj(c("ZerothOrderActivation", "ActivationBase"))
  .bm_activation_init(self, name)
  self[["_a"]] <- self[["_e"]]
  self
}

#' @rdname activation-dynamics
#' @export
ZerothOrderActivation_with_defaults <- function(name) ZerothOrderActivation(name)

.pvP_ZerothOrderActivation_order <- function(self) 0L
.pvP_ZerothOrderActivation_state_vars <- function(self) zeros(0L, 1L)
.pvP_ZerothOrderActivation_x <- function(self) zeros(0L, 1L)
.pvP_ZerothOrderActivation_input_vars <- function(self) .pv_colv(list(self[["_e"]]))
.pvP_ZerothOrderActivation_r <- function(self) .pv_colv(list(self[["_e"]]))
.pvP_ZerothOrderActivation_constants <- function(self) zeros(0L, 1L)
.pvP_ZerothOrderActivation_p <- function(self) zeros(0L, 1L)
.pvP_ZerothOrderActivation_M <- function(self) zeros(0L, 0L)
.pvP_ZerothOrderActivation_F <- function(self) zeros(0L, 1L)
.pvM_ZerothOrderActivation_rhs <- function(self) zeros(0L, 1L)

# ---------------------------------------------------------------------------
# FirstOrderActivationDeGroote2016

#' @rdname activation-dynamics
#' @export
FirstOrderActivationDeGroote2016 <- function(name, activation_time_constant = NULL, deactivation_time_constant = NULL,
                                             smoothing_rate = NULL) {
  self <- .pm_new_obj(c("FirstOrderActivationDeGroote2016", "ActivationBase"))
  .bm_activation_init(self, name)
  .pvS_FirstOrderActivationDeGroote2016_activation_time_constant(self, activation_time_constant)
  .pvS_FirstOrderActivationDeGroote2016_deactivation_time_constant(self, deactivation_time_constant)
  .pvS_FirstOrderActivationDeGroote2016_smoothing_rate(self, smoothing_rate)
  self
}

#' @rdname activation-dynamics
#' @export
FirstOrderActivationDeGroote2016_with_defaults <- function(name)
  FirstOrderActivationDeGroote2016(name, Float("0.015"), Float("0.060"), Float("10.0"))

.bm_constant_setter <- function(field, attr, prefix, sep_repr = FALSE) {
  force(field); force(attr); force(prefix); force(sep_repr)
  function(self, value) {
    if (.pm_has(self, field)) {
      old <- self[[field]]
      .pv_err("AttributeError", sprintf("Can't set attribute `%s` to %s as it is immutable and already has value %s.",
                                        attr, .pv_objrepr(value), if (sep_repr) .pv_objrepr(old) else .pv_str(old)))
    }
    self[[field]] <- if (is.null(value)) Symbol(paste0(prefix, self$name)) else value
  }
}

.pvS_FirstOrderActivationDeGroote2016_activation_time_constant <- .bm_constant_setter("_tau_a", "activation_time_constant", "tau_a_")
.pvS_FirstOrderActivationDeGroote2016_deactivation_time_constant <- .bm_constant_setter("_tau_d", "deactivation_time_constant", "tau_d_")
.pvS_FirstOrderActivationDeGroote2016_smoothing_rate <- .bm_constant_setter("_b", "smoothing_rate", "b_", TRUE)

.pvP_FirstOrderActivationDeGroote2016_activation_time_constant <- function(self) self[["_tau_a"]]
.pvP_FirstOrderActivationDeGroote2016_tau_a <- function(self) self[["_tau_a"]]
.pvP_FirstOrderActivationDeGroote2016_deactivation_time_constant <- function(self) self[["_tau_d"]]
.pvP_FirstOrderActivationDeGroote2016_tau_d <- function(self) self[["_tau_d"]]
.pvP_FirstOrderActivationDeGroote2016_smoothing_rate <- function(self) self[["_b"]]
.pvP_FirstOrderActivationDeGroote2016_b <- function(self) self[["_b"]]
.pvP_FirstOrderActivationDeGroote2016_order <- function(self) 1L
.pvP_FirstOrderActivationDeGroote2016_state_vars <- function(self) .pv_colv(list(self[["_a"]]))
.pvP_FirstOrderActivationDeGroote2016_x <- function(self) .pv_colv(list(self[["_a"]]))
.pvP_FirstOrderActivationDeGroote2016_input_vars <- function(self) .pv_colv(list(self[["_e"]]))
.pvP_FirstOrderActivationDeGroote2016_r <- function(self) .pv_colv(list(self[["_e"]]))

.bm_symbolic_constants <- function(cs) .bm_colmat(Filter(function(c) !.expr_is_number(as_sym(c)), cs))

.pvP_FirstOrderActivationDeGroote2016_constants <- function(self)
  .bm_symbolic_constants(list(self[["_tau_a"]], self[["_tau_d"]], self[["_b"]]))
.pvP_FirstOrderActivationDeGroote2016_p <- .pvP_FirstOrderActivationDeGroote2016_constants
.pvP_FirstOrderActivationDeGroote2016_M <- function(self) .pv_colv(list(S_One))
.pvP_FirstOrderActivationDeGroote2016_F <- function(self) .pv_colv(list(.bm_da_eqn(self)))
.pvM_FirstOrderActivationDeGroote2016_rhs <- function(self) .pv_colv(list(.bm_da_eqn(self)))

# the cached property _da_eqn
.bm_da_eqn <- function(self) {
  if (.pm_has(self, "_da_eqn")) return(self[["_da_eqn"]])
  HALF <- Rational(1L, 2L)
  b <- as_sym(self[["_b"]]); e <- self[["_e"]]; a <- self[["_a"]]
  a0 <- HALF * tanh(b * (e - a))
  a1 <- HALF + Rational(3L, 2L) * a
  a2 <- (HALF + a0) / (as_sym(self[["_tau_a"]]) * a1)
  a3 <- a1 * (HALF - a0) / as_sym(self[["_tau_d"]])
  r <- (a2 + a3) * (e - a)
  self[["_da_eqn"]] <- r
  r
}

.pvM_FirstOrderActivationDeGroote2016___eq__ <- function(self, other) {
  if (!identical(class(self)[1L], class(other)[1L])) return(FALSE)
  same <- function(p, q) if (is.character(p) || is.character(q)) identical(p, q) else identical(as_sym(p)$key, as_sym(q)$key)
  same(self$name, other$name) && same(self$tau_a, other$tau_a) && same(self$tau_d, other$tau_d) && same(self$b, other$b)
}

.pvR_FirstOrderActivationDeGroote2016 <- function(x)
  sprintf("%s(%s, activation_time_constant=%s, deactivation_time_constant=%s, smoothing_rate=%s)",
          class(x)[1L], .pv_objrepr(x$name), .pv_objrepr(x$tau_a), .pv_objrepr(x$tau_d), .pv_objrepr(x$b))
