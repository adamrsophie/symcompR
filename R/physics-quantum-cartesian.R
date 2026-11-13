# Operators and states of one and three dimensional Cartesian position and
# momentum, after sympy/physics/quantum/cartesian.py.

.qL2_real <- function(cls, args) L2(Interval(S_NegativeInfinity, S_Infinity))

.qclass("XOp", "HermitianOperator",
  default_args = .cm(function(cls) list("X")),
  `_eval_hilbert_space` = .cm(.qL2_real),
  `_eval_commutator_PxOp` = function(self, other, ...) .mul(list(sym_I, hbar)),
  `_apply_operator_XKet` = function(self, ket, ...) .qtimes(.qa(ket, "position"), ket),
  `_apply_operator_PositionKet3D` = function(self, ket, ...) .qtimes(.qa(ket, "position_x"), ket),
  `_represent_PxKet` = function(self, basis, index = 1, ...) {
    states <- .qmeth(basis, "_enumerate_state", 2, start_index = index)
    coord1 <- .qa(states[[1L]], "momentum")
    coord2 <- .qa(states[[2L]], "momentum")
    d <- DifferentialOperator(coord1)
    delta <- fcall("DiracDelta", list(.sub(coord1, coord2)))
    .mul(list(sym_I, hbar, .qtimes(d, delta)))
  })

.qclass("YOp", "HermitianOperator",
  default_args = .cm(function(cls) list("Y")),
  `_eval_hilbert_space` = .cm(.qL2_real),
  `_apply_operator_PositionKet3D` = function(self, ket, ...) .qtimes(.qa(ket, "position_y"), ket))

.qclass("ZOp", "HermitianOperator",
  default_args = .cm(function(cls) list("Z")),
  `_eval_hilbert_space` = .cm(.qL2_real),
  `_apply_operator_PositionKet3D` = function(self, ket, ...) .qtimes(.qa(ket, "position_z"), ket))

.qclass("PxOp", "HermitianOperator",
  default_args = .cm(function(cls) list("Px")),
  `_eval_hilbert_space` = .cm(.qL2_real),
  `_apply_operator_PxKet` = function(self, ket, ...) .qtimes(.qa(ket, "momentum"), ket),
  `_represent_XKet` = function(self, basis, index = 1, ...) {
    states <- .qmeth(basis, "_enumerate_state", 2, start_index = index)
    coord1 <- .qa(states[[1L]], "position")
    coord2 <- .qa(states[[2L]], "position")
    d <- DifferentialOperator(coord1)
    delta <- fcall("DiracDelta", list(.sub(coord1, coord2)))
    .mul(list(S_NegativeOne, sym_I, hbar, .qtimes(d, delta)))
  })

.qlabel_str <- function(x) if (is_sym(x)) .str_top(x) else as.character(x)

.qlowercase_labels <- function(ops) {
  if (!(is.list(ops) && !is_sym(ops))) ops <- list(ops)
  lapply(ops, function(arg) tolower(.qlabel_str(.qa(arg, "label")[[1L]])))
}

.quppercase_labels <- function(ops) {
  if (!(is.list(ops) && !is_sym(ops))) ops <- list(ops)
  lapply(ops, function(arg) {
    s <- .qlabel_str(.qa(arg, "label")[[1L]])
    paste0(toupper(substr(s, 1L, 1L)), substring(s, 2L))
  })
}

.qops_to_lower_state <- .cm(function(cls, op, ...) .qfunc(cls, .qlowercase_labels(op)))
.qstate_to_upper_ops <- function(self, op_class, ...) {
  if (is.list(op_class)) op_class <- op_class[[1L]]
  .qfunc(op_class, .quppercase_labels(self))
}

.qenumerate_continuous_1D <- function(state, num_states, index_list = list(), start_index = 1, ...) {
  if (!length(index_list)) index_list <- as.list(seq(.qint(start_index), length.out = .qint(num_states)))
  label <- state$args[[1L]]
  lapply(index_list, function(ind) .qfunc(state$type, list(paste0(.qlabel_str(label), "_", .qint(ind)))))
}

.qclass("XKet", "Ket",
  `_operators_to_state` = .qops_to_lower_state,
  `_state_to_operators` = .qstate_to_upper_ops,
  default_args = .cm(function(cls) list("x")),
  dual_class = .cm(function(cls) "XBra"),
  position = .prop(function(self) .qa(self, "label")[[1L]]),
  `_enumerate_state` = function(self, num_states, ...) .qenumerate_continuous_1D(self, num_states, ...),
  `_eval_innerproduct_XBra` = function(self, bra, ...) fcall("DiracDelta", list(.sub(.qa(self, "position"), .qa(bra, "position")))),
  `_eval_innerproduct_PxBra` = function(self, bra, ...) {
    .mul(list(exp(.mul(list(S_NegativeOne, sym_I, .qa(self, "position"), .qa(bra, "momentum"), .pow(hbar, S_NegativeOne)))),
              .pow(sqrt(.mul(list(S_Two, sym_pi, hbar))), S_NegativeOne)))
  })

.qclass("XBra", "Bra",
  default_args = .cm(function(cls) list("x")),
  dual_class = .cm(function(cls) "XKet"),
  position = .prop(function(self) .qa(self, "label")[[1L]]))

.qclass("PositionState3D", "State",
  `_operators_to_state` = .qops_to_lower_state,
  `_state_to_operators` = .qstate_to_upper_ops,
  default_args = .cm(function(cls) list("x", "y", "z")),
  position_x = .prop(function(self) .qa(self, "label")[[1L]]),
  position_y = .prop(function(self) .qa(self, "label")[[2L]]),
  position_z = .prop(function(self) .qa(self, "label")[[3L]]))

.qclass("PositionKet3D", c("Ket", "PositionState3D"),
  `_eval_innerproduct_PositionBra3D` = function(self, bra, ...) {
    dd <- function(nm) fcall("DiracDelta", list(.sub(.qa(self, nm), .qa(bra, nm))))
    .mul(list(dd("position_x"), dd("position_y"), dd("position_z")))
  },
  dual_class = .cm(function(cls) "PositionBra3D"))

.qclass("PositionBra3D", c("Bra", "PositionState3D"),
  dual_class = .cm(function(cls) "PositionKet3D"))

.qclass("PxKet", "Ket",
  `_operators_to_state` = .qops_to_lower_state,
  `_state_to_operators` = .qstate_to_upper_ops,
  default_args = .cm(function(cls) list("px")),
  dual_class = .cm(function(cls) "PxBra"),
  momentum = .prop(function(self) .qa(self, "label")[[1L]]),
  `_enumerate_state` = function(self, num_states, ...) .qenumerate_continuous_1D(self, num_states, ...),
  `_eval_innerproduct_XBra` = function(self, bra, ...) {
    .mul(list(exp(.mul(list(sym_I, .qa(self, "momentum"), .qa(bra, "position"), .pow(hbar, S_NegativeOne)))),
              .pow(sqrt(.mul(list(S_Two, sym_pi, hbar))), S_NegativeOne)))
  },
  `_eval_innerproduct_PxBra` = function(self, bra, ...) fcall("DiracDelta", list(.sub(.qa(self, "momentum"), .qa(bra, "momentum")))))

.qclass("PxBra", "Bra",
  default_args = .cm(function(cls) list("px")),
  dual_class = .cm(function(cls) "PxKet"),
  momentum = .prop(function(self) .qa(self, "label")[[1L]]))

#' Cartesian position and momentum
#'
#' The operators and states of `sympy.physics.quantum.cartesian`: the
#' position operators `XOp()`, `YOp()`, `ZOp()` and the momentum operator
#' `PxOp()`, with their eigenstates `XKet()`, `XBra()`, `PxKet()`, `PxBra()`
#' and the three dimensional `PositionKet3D()` and `PositionBra3D()`.
#' Inner products give delta functions and plane waves; [represent()] gives
#' the differential operators.
#'
#' @param ... Labels.
#' @return A quantum operator or state.
#' @examples
#' x <- Symbol("x"); px <- Symbol("px")
#' qapply(XOp() * XKet(x))
#' XBra(x) * PxKet(px)
#' @name cartesian-quantum
NULL

#' @rdname cartesian-quantum
#' @export
XOp <- function(...) .qfunc("XOp", list(...))
#' @rdname cartesian-quantum
#' @export
YOp <- function(...) .qfunc("YOp", list(...))
#' @rdname cartesian-quantum
#' @export
ZOp <- function(...) .qfunc("ZOp", list(...))
#' @rdname cartesian-quantum
#' @export
PxOp <- function(...) .qfunc("PxOp", list(...))
#' @rdname cartesian-quantum
#' @export
XKet <- function(...) .qfunc("XKet", list(...))
#' @rdname cartesian-quantum
#' @export
XBra <- function(...) .qfunc("XBra", list(...))
#' @rdname cartesian-quantum
#' @export
PositionState3D <- function(...) .qfunc("PositionState3D", list(...))
#' @rdname cartesian-quantum
#' @export
PositionKet3D <- function(...) .qfunc("PositionKet3D", list(...))
#' @rdname cartesian-quantum
#' @export
PositionBra3D <- function(...) .qfunc("PositionBra3D", list(...))
#' @rdname cartesian-quantum
#' @export
PxKet <- function(...) .qfunc("PxKet", list(...))
#' @rdname cartesian-quantum
#' @export
PxBra <- function(...) .qfunc("PxBra", list(...))
