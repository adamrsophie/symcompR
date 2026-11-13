# Bosonic and fermionic operators and Fock states, after
# sympy/physics/quantum/boson.py and fermion.py.

.qladder_args <- function(args) {
  if (!(length(args) %in% c(1L, 2L)))
    .qstop("ValueError", "1 or 2 parameters expected, got ", .qrepr(args))
  if (length(args) == 1L) return(list(args[[1L]], S_One))
  b <- args[[2L]]
  list(args[[1L]], if (is.logical(b)) as_sym(as.integer(b)) else Integer(b))
}

.qindependent <- function(...) isTRUE(list(...)$independent)

.qladder_printers <- list(
  name = .prop(function(self) self$args[[1L]]),
  is_annihilation = .prop(function(self) .qtruthy(self$args[[2L]])),
  `_print_contents_latex` = function(self, p) {
    if (.qa(self, "is_annihilation")) sprintf("{%s}", .str(.qa(self, "name")))
    else sprintf("{{%s}^\\dagger}", .str(.qa(self, "name")))
  },
  `_print_contents` = function(self, p) {
    if (.qa(self, "is_annihilation")) .str(.qa(self, "name")) else sprintf("Dagger(%s)", .str(.qa(self, "name")))
  })

do.call(.qclass, c(list("BosonOp", "Operator",
  default_args = .cm(function(cls) list("a", TRUE)),
  `__new__` = .cm(function(cls, args, ...) .qsuper("Operator", "__new__", cls, .qladder_args(args))),
  `_eval_commutator_BosonOp` = function(self, other, ...) {
    if (.qeq(.qa(self, "name"), .qa(other, "name"))) {
      if (!.qa(self, "is_annihilation") && .qa(other, "is_annihilation")) return(S_NegativeOne)
    } else if (.qindependent(...)) {
      return(S_Zero)
    }
    NULL
  },
  `_eval_commutator_FermionOp` = function(self, other, ...) S_Zero,
  `_eval_anticommutator_BosonOp` = function(self, other, ...) {
    if (.qindependent(...)) return(.mul(list(S_Two, self, other)))
    NULL
  },
  `_eval_adjoint` = function(self) BosonOp(.str(.qa(self, "name")), !.qa(self, "is_annihilation"))),
  .qladder_printers))

.qclass("BosonFockKet", "Ket",
  `__new__` = .cm(function(cls, args, ...) {
    if (length(args) != 1L) .qstop("TypeError", "BosonFockKet.__new__() takes 2 positional arguments but ", length(args) + 1L, " were given")
    .qsuper("Ket", "__new__", cls, args)
  }),
  n = .prop(function(self) .qa(self, "label")[[1L]]),
  dual_class = .cm(function(cls) "BosonFockBra"),
  `_eval_hilbert_space` = .cm(function(cls, label) FockSpace()),
  `_eval_innerproduct_BosonFockBra` = function(self, bra, ...) KroneckerDelta(.qa(self, "n"), .qa(bra, "n")),
  `_apply_from_right_to_BosonOp` = function(self, op, ...) {
    n <- .qa(self, "n")
    if (.qa(op, "is_annihilation")) .mul(list(.sqrt(n), BosonFockKet(.sub(n, S_One))))
    else .mul(list(.sqrt(.add(list(n, S_One))), BosonFockKet(.add(list(n, S_One)))))
  })

.qclass("BosonFockBra", "Bra",
  `__new__` = .cm(function(cls, args, ...) {
    if (length(args) != 1L) .qstop("TypeError", "BosonFockBra.__new__() takes 2 positional arguments but ", length(args) + 1L, " were given")
    .qsuper("Bra", "__new__", cls, args)
  }),
  n = .prop(function(self) .qa(self, "label")[[1L]]),
  dual_class = .cm(function(cls) "BosonFockKet"),
  `_eval_hilbert_space` = .cm(function(cls, label) FockSpace()))

.qclass("BosonCoherentKet", "Ket",
  alpha = .prop(function(self) .qa(self, "label")[[1L]]),
  dual_class = .cm(function(cls) "BosonCoherentBra"),
  `_eval_hilbert_space` = .cm(function(cls, label) HilbertSpace()),
  `_eval_innerproduct_BosonCoherentBra` = function(self, bra, ...) {
    a <- .qa(self, "alpha"); b <- .qa(bra, "alpha")
    if (.qeq(a, b)) return(S_One)
    exp(.mul(list(S_NegativeOne, .add(list(.pow(Abs(a), S_Two), .pow(Abs(b), S_Two),
                                            .mul(list(S_NegativeOne, S_Two, conjugate(b), a)))), S_Half)))
  },
  `_apply_from_right_to_BosonOp` = function(self, op, ...) {
    if (.qa(op, "is_annihilation")) return(.mul(list(.qa(self, "alpha"), self)))
    NULL
  })

.qclass("BosonCoherentBra", "Bra",
  alpha = .prop(function(self) .qa(self, "label")[[1L]]),
  dual_class = .cm(function(cls) "BosonCoherentKet"),
  `_apply_operator_BosonOp` = function(self, op, ...) {
    if (!.qa(op, "is_annihilation")) return(.mul(list(.qa(self, "alpha"), self)))
    NULL
  })

#' Bosonic operators and states
#'
#' `BosonOp(name, annihilation = TRUE)` is the bosonic ladder operator of
#' `sympy.physics.quantum.boson` (its adjoint is the creation operator),
#' `BosonFockKet(n)` and `BosonFockBra(n)` are number states and
#' `BosonCoherentKet(alpha)` and `BosonCoherentBra(alpha)` coherent states.
#'
#' @param ... The name and whether the operator annihilates.
#' @param n The occupation number.
#' @param alpha The eigenvalue of the annihilation operator.
#' @return A quantum operator or state.
#' @examples
#' a <- BosonOp("a")
#' Dagger(a)
#' qapply(a * BosonFockKet(2))
#' doit(Commutator(a, Dagger(a)))
#' @name boson
NULL

#' @rdname boson
#' @export
BosonOp <- function(...) .qfunc("BosonOp", list(...))
#' @rdname boson
#' @export
BosonFockKet <- function(n) .qfunc("BosonFockKet", list(n))
#' @rdname boson
#' @export
BosonFockBra <- function(n) .qfunc("BosonFockBra", list(n))
#' @rdname boson
#' @export
BosonCoherentKet <- function(alpha) .qfunc("BosonCoherentKet", list(alpha))
#' @rdname boson
#' @export
BosonCoherentBra <- function(alpha) .qfunc("BosonCoherentBra", list(alpha))

# ---------------------------------------------------------------------------
# Fermions

do.call(.qclass, c(list("FermionOp", "Operator",
  default_args = .cm(function(cls) list("c", TRUE)),
  `__new__` = .cm(function(cls, args, ...) .qsuper("Operator", "__new__", cls, .qladder_args(args))),
  `_eval_commutator_FermionOp` = function(self, other, ...) {
    if (.qindependent(...)) return(S_Zero)
    NULL
  },
  `_eval_anticommutator_FermionOp` = function(self, other, ...) {
    if (.qeq(.qa(self, "name"), .qa(other, "name"))) {
      if (!.qa(self, "is_annihilation") && .qa(other, "is_annihilation")) return(S_One)
    } else if (.qindependent(...)) {
      return(.mul(list(S_Two, self, other)))
    }
    NULL
  },
  `_eval_anticommutator_BosonOp` = function(self, other, ...) .mul(list(S_Two, self, other)),
  `_eval_commutator_BosonOp` = function(self, other, ...) S_Zero,
  `_eval_adjoint` = function(self) FermionOp(.str(.qa(self, "name")), !.qa(self, "is_annihilation")),
  `_eval_power` = function(self, ex) {
    if (.qeq(ex, S_Zero)) return(S_One)
    if (.qeq(ex, S_One)) return(self)
    if (isTRUE(.qgt(ex, S_One)) && isTRUE(ask_fact(ex, "integer"))) return(S_Zero)
    if (isTRUE(.qgt(S_Zero, ex)) || identical(ask_fact(ex, "integer"), FALSE))
      .qstop("ValueError", "Fermionic operators can only be raised to a positive integer power")
    NULL
  }), .qladder_printers))

.qfermion_n <- function(args, cls) {
  if (length(args) != 1L) .qstop("TypeError", cls, ".__new__() takes 2 positional arguments but ", length(args) + 1L, " were given")
  n <- args[[1L]]
  ok <- (is.numeric(n) || is.logical(n)) && as.numeric(n) %in% c(0, 1) ||
    is_sym(n) && n$type == "Integer" && .qint(n) %in% c(0, 1)
  if (!ok) .qstop("ValueError", "n must be 0 or 1")
}

.qclass("FermionFockKet", "Ket",
  `__new__` = .cm(function(cls, args, ...) {
    .qfermion_n(args, "FermionFockKet")
    .qsuper("Ket", "__new__", cls, args)
  }),
  n = .prop(function(self) .qa(self, "label")[[1L]]),
  dual_class = .cm(function(cls) "FermionFockBra"),
  `_eval_hilbert_space` = .cm(function(cls, label) HilbertSpace()),
  `_eval_innerproduct_FermionFockBra` = function(self, bra, ...) KroneckerDelta(.qa(self, "n"), .qa(bra, "n")),
  `_apply_from_right_to_FermionOp` = function(self, op, ...) {
    n <- .qint(.qa(self, "n"))
    if (.qa(op, "is_annihilation")) {
      if (n == 1) FermionFockKet(0) else S_Zero
    } else if (n == 0) FermionFockKet(1) else S_Zero
  })

.qclass("FermionFockBra", "Bra",
  `__new__` = .cm(function(cls, args, ...) {
    .qfermion_n(args, "FermionFockBra")
    .qsuper("Bra", "__new__", cls, args)
  }),
  n = .prop(function(self) .qa(self, "label")[[1L]]),
  dual_class = .cm(function(cls) "FermionFockKet"))

#' Fermionic operators and states
#'
#' `FermionOp(name, annihilation = TRUE)` is the fermionic ladder operator
#' of `sympy.physics.quantum.fermion` and `FermionFockKet(n)` and
#' `FermionFockBra(n)` its number states, with `n` 0 or 1.
#'
#' @param ... The name and whether the operator annihilates.
#' @param n The occupation number, 0 or 1.
#' @return A quantum operator or state.
#' @examples
#' c1 <- FermionOp("c")
#' doit(AntiCommutator(c1, Dagger(c1)))
#' qapply(Dagger(c1) * FermionFockKet(0))
#' @name fermion
NULL

#' @rdname fermion
#' @export
FermionOp <- function(...) .qfunc("FermionOp", list(...))
#' @rdname fermion
#' @export
FermionFockKet <- function(n) .qfunc("FermionFockKet", list(n))
#' @rdname fermion
#' @export
FermionFockBra <- function(n) .qfunc("FermionFockBra", list(n))

# ---------------------------------------------------------------------------
# Operator ordering (sympy/physics/quantum/operatorordering.py)

.qexpand_powers <- function(factors) {
  out <- list()
  for (f in .qargs(factors)) {
    if (f$type == "Pow" && f$args[[2L]]$type == "Integer" && .qint(f$args[[2L]]) > 0) {
      for (i in seq_len(.qint(f$args[[2L]]))) out[[length(out) + 1L]] <- f$args[[1L]]
    } else {
      out[[length(out) + 1L]] <- f
    }
  }
  out
}

.qsame_list <- function(a, b) length(a) == length(b) && identical(vapply(a, .k, ""), vapply(b, .k, ""))

.qnormal_ordered_form_factor <- function(product, independent, recursive_limit, depth) {
  factors <- .qexpand_powers(product)
  new_factors <- list()
  add <- function(x) new_factors[[length(new_factors) + 1L]] <<- x
  n <- 1L
  nf <- length(factors)
  while (n < nf) {
    current <- factors[[n]]; nxt <- factors[[n + 1L]]
    if (!all(vapply(list(current, nxt), function(f) .qisa(f, c("FermionOp", "BosonOp")), NA))) {
      add(current); n <- n + 1L; next
    }
    k1 <- list(as.integer(.qa(current, "is_annihilation")), .str(.qa(current, "name")))
    k2 <- list(as.integer(.qa(nxt, "is_annihilation")), .str(.qa(nxt, "name")))
    if (.key_cmp(k1, k2) <= 0) {
      add(current); n <- n + 1L; next
    }
    n <- n + 2L
    if (.qa(current, "is_annihilation") && !.qa(nxt, "is_annihilation")) {
      if (.qisa(current, "BosonOp") && .qisa(nxt, "BosonOp")) {
        if (!.qeq(current$args[[1L]], nxt$args[[1L]])) {
          c <- if (independent) S_Zero else Commutator(current, nxt)
          add(.add(list(.mul(list(nxt, current)), c)))
        } else {
          add(.add(list(.mul(list(nxt, current)), S_One)))
        }
      } else if (.qisa(current, "FermionOp") && .qisa(nxt, "FermionOp")) {
        if (!.qeq(current$args[[1L]], nxt$args[[1L]])) {
          c <- if (independent) S_Zero else AntiCommutator(current, nxt)
          add(.add(list(.mul(list(S_NegativeOne, nxt, current)), c)))
        } else {
          add(.add(list(.mul(list(S_NegativeOne, nxt, current)), S_One)))
        }
      }
    } else if (.qa(current, "is_annihilation") == .qa(nxt, "is_annihilation") &&
               .qisa(current, "FermionOp") && .qisa(nxt, "FermionOp")) {
      add(.mul(list(S_NegativeOne, nxt, current)))
    } else {
      add(.mul(list(nxt, current)))
    }
  }
  if (n == nf) add(factors[[nf]])
  if (.qsame_list(new_factors, factors)) return(product)
  normal_ordered_form(expand(.mul(new_factors)), recursive_limit = recursive_limit,
                      .recursive_depth = depth + 1L, independent = independent)
}

#' Normal ordering of bosonic and fermionic operators
#'
#' `normal_ordered_form()` moves creation operators to the left of
#' annihilation operators using the commutation relations, and
#' `normal_order()` does the same without the commutator terms, as in
#' `sympy.physics.quantum.operatorordering`.
#'
#' @param expr An expression of [BosonOp()] and [FermionOp()] operators.
#' @param independent Treat operators with different names as independent.
#' @param recursive_limit The maximum number of recursions.
#' @param .recursive_depth The current recursion depth.
#' @return An expression.
#' @examples
#' a <- BosonOp("a")
#' normal_ordered_form(a * Dagger(a))
#' normal_order(a * Dagger(a))
#' @name operator-ordering
NULL

#' @rdname operator-ordering
#' @export
normal_ordered_form <- function(expr, independent = FALSE, recursive_limit = 10, .recursive_depth = 0L) {
  if (.recursive_depth > recursive_limit) {
    warning("Too many recursions, aborting", call. = FALSE)
    return(expr)
  }
  if (expr$type == "Add")
    return(.add(lapply(.qargs(expr), function(t)
      if (t$type == "Mul") .qnormal_ordered_form_factor(t, independent, recursive_limit, .recursive_depth) else t)))
  if (expr$type == "Mul") return(.qnormal_ordered_form_factor(expr, independent, recursive_limit, .recursive_depth))
  expr
}

.qnormal_order_factor <- function(product, recursive_limit, depth) {
  factors <- .qexpand_powers(product)
  new_factors <- list()
  add <- function(x) new_factors[[length(new_factors) + 1L]] <<- x
  n <- 1L
  nf <- length(factors)
  while (n < nf) {
    f <- factors[[n]]; g <- factors[[n + 1L]]
    if (.qisa(f, "BosonOp") && .qa(f, "is_annihilation")) {
      if (!.qisa(g, "BosonOp")) add(f)
      else if (.qa(g, "is_annihilation")) add(f)
      else { add(.mul(list(g, f))); n <- n + 1L }
    } else if (.qisa(f, "FermionOp") && .qa(f, "is_annihilation")) {
      if (!.qisa(g, "FermionOp")) add(f)
      else if (.qa(g, "is_annihilation")) add(f)
      else { add(.mul(list(S_NegativeOne, g, f))); n <- n + 1L }
    } else {
      add(f)
    }
    n <- n + 1L
  }
  if (n == nf) add(factors[[nf]])
  if (.qsame_list(new_factors, factors)) return(product)
  normal_order(expand(.mul(new_factors)), recursive_limit = recursive_limit, .recursive_depth = depth + 1L)
}

#' @rdname operator-ordering
#' @export
normal_order <- function(expr, recursive_limit = 10, .recursive_depth = 0L) {
  if (.recursive_depth > recursive_limit) {
    warning("Too many recursions, aborting", call. = FALSE)
    return(expr)
  }
  if (expr$type == "Add")
    return(.add(lapply(.qargs(expr), function(t)
      if (t$type == "Mul") .qnormal_order_factor(t, recursive_limit, .recursive_depth) else t)))
  if (expr$type == "Mul") return(.qnormal_order_factor(expr, recursive_limit, .recursive_depth))
  expr
}

# ---------------------------------------------------------------------------
# Particle in a box (sympy/physics/quantum/piab.py)

.qpiab_hs <- function() L2(Interval(S_NegativeInfinity, oo))

.qclass("PIABHamiltonian", "HermitianOperator",
  `_eval_hilbert_space` = .cm(function(cls, label) .qpiab_hs()),
  `_apply_operator_PIABKet` = function(self, ket, ...) {
    n <- .qa(ket, "label")[[1L]]
    m <- Symbol("m"); L <- Symbol("L")
    .mul(list(.pow(n, S_Two), .pow(sym_pi, S_Two), .pow(hbar, S_Two),
              .pow(.mul(list(S_Two, m, .pow(L, S_Two))), S_NegativeOne), ket))
  })

.qclass("PIABKet", "Ket",
  `_eval_hilbert_space` = .cm(function(cls, args) .qpiab_hs()),
  dual_class = .cm(function(cls) "PIABBra"),
  `_represent_default_basis` = function(self, ...) .qmeth(self, "_represent_XOp", NULL, ...),
  `_represent_XOp` = function(self, basis, subs = NULL, ...) {
    x <- Symbol("x"); n <- Symbol("n"); L <- Symbol("L")
    s <- sin(.mul(list(n, sym_pi, x, .pow(L, S_NegativeOne))))
    if (length(subs)) s <- subs(s, subs)
    .mul(list(.sqrt(.mul(list(S_Two, .pow(L, S_NegativeOne)))), s))
  },
  `_eval_innerproduct_PIABBra` = function(self, bra, ...) KroneckerDelta(.qa(bra, "label")[[1L]], .qa(self, "label")[[1L]]))

.qclass("PIABBra", "Bra",
  `_eval_hilbert_space` = .cm(function(cls, label) .qpiab_hs()),
  dual_class = .cm(function(cls) "PIABKet"))

#' Particle in a box
#'
#' The Hamiltonian `PIABHamiltonian()` and the eigenstates `PIABKet(n)`
#' and `PIABBra(n)` of a particle of mass `m` in a box of length `L`, as in
#' `sympy.physics.quantum.piab`.
#'
#' @param ... The label.
#' @return A quantum operator or state.
#' @examples
#' n <- Symbol("n")
#' qapply(PIABHamiltonian("H") * PIABKet(n))
#' @name piab
NULL

#' @rdname piab
#' @export
PIABHamiltonian <- function(...) .qfunc("PIABHamiltonian", list(...))
#' @rdname piab
#' @export
PIABKet <- function(...) .qfunc("PIABKet", list(...))
#' @rdname piab
#' @export
PIABBra <- function(...) .qfunc("PIABBra", list(...))
