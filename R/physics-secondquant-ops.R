# Second quantization: operators, Fock states and their algebra, after
# sympy/physics/secondquant.py. The classes live in the quantum class table
# .QC (see physics-quantum-aaa-qexpr.R). The secondquant Dagger, Commutator
# and InnerProduct are sq_Dagger(), sq_Commutator() and sq_InnerProduct(),
# since the names Dagger and Commutator belong to sympy.physics.quantum.

.SQ_ERR <- "SecondQuantizationError"
.sq_stop <- function(name, ...) stop(.qerror(name, paste0(...), bases = .SQ_ERR))

# index.assumptions0.get(nm)
.sq_af <- function(e, nm) is_sym(e) && isTRUE(e$extra$other[nm][[1L]])

# a tie breaker standing in for Python's hash()
.sq_tb <- function(e) e$key

.sq_seq <- function(x) {
  if (is_sym(x) && x$type == "Tuple") return(x$args)
  if (is_sym(x)) return(list(x))
  lapply(as.list(x), .qS)
}

.sq_raw <- function(cls, args) .qraw(cls, NULL, args)
.sq_craw <- function(cls, args) {
  .qm_state$on <- TRUE
  new_sym(cls, args)
}

.sq_basic_str <- function(self, p) paste0(.qpyname(self$type), "(", paste(vapply(self$args, .str, ""), collapse = ", "), ")")

# ---------------------------------------------------------------------------
# Dagger

.sq_dagger_eval <- function(arg) {
  if (!is.null(.qnode_def(arg$type)) && .qhas(arg, "_dagger_")) return(.qmeth(arg, "_dagger_"))
  if (arg$type %in% c("Symbol", "Dummy") && .is_comm(arg)) return(conjugate(arg))
  if (arg$type == "Add") return(.add(lapply(.sympy_args(arg), sq_Dagger)))
  if (arg$type == "Mul") return(.mul(lapply(rev(.sympy_args(arg)), sq_Dagger)))
  if (is_number_type(arg)) return(arg)
  if (arg$type == "Pow") return(.pow(sq_Dagger(arg$args[[1L]]), arg$args[[2L]]))
  if (identical(arg$key, sym_I$key)) return(.neg(arg))
  if (.qisa(arg, "Function") && all(vapply(arg$args, .is_comm, NA)))
    return(sym_rebuild(arg, lapply(arg$args, sq_Dagger)))
  NULL
}

.qclass("sqDagger", "Expr", pyname = "Dagger", is_commutative = FALSE,
  `__new__` = .cm(function(cls, args, ...) {
    arg <- .qS(args[[1L]])
    r <- .sq_dagger_eval(arg)
    if (!is.null(r)) return(r)
    .sq_raw(cls, list(arg))
  }),
  `_dagger_` = function(self) self$args[[1L]],
  `_sympystr` = .sq_basic_str,
  `_sympyrepr` = .qbasic_printers[[2L]],
  `_latex` = function(self, p) paste0("\\operatorname{Dagger}{\\left(", .latex(self$args[[1L]]), " \\right)}"))

#' Second quantization
#'
#' Creation and annihilation operators of bosons and fermions, occupation
#' number (Fock) states and the algebra of `sympy.physics.secondquant`.
#' `B(k)` and `Bd(k)` annihilate and create a boson in state `k`; `sym_F(k)`
#' and `Fd(k)` do the same for a fermion. A symbolic index may carry the
#' assumption `above_fermi = TRUE` (a particle state) or
#' `below_fermi = TRUE` (a hole state). `BKet()`, `BBra()`, `FKet()` and
#' `FBra()` are occupation number states. `sq_Dagger()` is the Hermitian
#' conjugate, `sq_Commutator()` the commutator, `sq_InnerProduct()` the
#' inner product of a bra and a ket and `apply_operators()` applies the
#' operators of an expression to the states on their right.
#'
#' @param k The state (index) acted upon: an integer or a symbol.
#' @param arg,a,b,e Expressions.
#' @param occupations The occupation numbers (bosons) or the occupied
#'   states (fermions), as a list.
#' @param fermi_level The number of states below the Fermi level.
#' @param symbol,upper,lower The name, upper and lower indices of an
#'   antisymmetric tensor.
#' @param bra,ket A bra and a ket.
#' @param op An operator.
#' @param basis A basis from `VarBosonicBasis()` or `FixedBosonicBasis()`.
#' @param n_max The number of states of a variable particle number basis.
#' @param n_particles,n_levels The number of particles and of levels of a
#'   fixed particle number basis.
#' @return An expression, a list or, for the bases, an object of class
#'   `VarBosonicBasis` or `FixedBosonicBasis`.
#' @examples
#' p <- Symbol("p"); q <- Symbol("q")
#' Fd(p) * sym_F(q)
#' sq_Dagger(B(0))
#' apply_operators(Bd(0) * BKet(list(1, 2)))
#' sq_Commutator(B(0), Bd(0))
#' @name secondquant
NULL

#' @rdname secondquant
#' @export
sq_Dagger <- function(arg) .qfunc("sqDagger", list(arg))

# ---------------------------------------------------------------------------
# Antisymmetric tensors

.qclass("TensorSymbol", "Expr", node = FALSE, is_commutative = TRUE)

.qclass("AntiSymmetricTensor", "TensorSymbol",
  `__new__` = .cm(function(cls, args, ...) {
    su <- tryCatch(.sq_sort_anticommuting_fermions(.sq_seq(args[[2L]]), .sq_key_index),
                   ViolationOfPauliPrinciple = function(e) NULL)
    if (is.null(su)) return(S_Zero)
    sl <- tryCatch(.sq_sort_anticommuting_fermions(.sq_seq(args[[3L]]), .sq_key_index),
                   ViolationOfPauliPrinciple = function(e) NULL)
    if (is.null(sl)) return(S_Zero)
    obj <- .sq_craw(cls, list(.qS(args[[1L]]), do.call(Tuple, su[[1L]]), do.call(Tuple, sl[[1L]])))
    if ((su[[2L]] + sl[[2L]]) %% 2) .neg(obj) else obj
  }),
  symbol = .prop(function(self) self$args[[1L]]),
  upper = .prop(function(self) self$args[[2L]]),
  lower = .prop(function(self) self$args[[3L]]),
  `__str__` = function(self) sprintf("%s(%s,%s)", .str(self$args[[1L]]), .str(self$args[[2L]]), .str(self$args[[3L]])),
  `_sympystr` = .sq_basic_str,
  `_sympyrepr` = .qbasic_printers[[2L]],
  `_latex` = function(self, p) sprintf("{%s^{%s}_{%s}}", .str(self$args[[1L]]),
                                       paste(vapply(self$args[[2L]]$args, .latex, ""), collapse = ""),
                                       paste(vapply(self$args[[3L]]$args, .latex, ""), collapse = "")))

#' @rdname secondquant
#' @export
AntiSymmetricTensor <- function(symbol, upper, lower) .qfunc("AntiSymmetricTensor", list(symbol, upper, lower))

# ---------------------------------------------------------------------------
# Operators

.sq_op_str <- function(self) sprintf("%s(%s)", .qa(self, "op_symbol"), .str(self$args[[1L]]))
.sq_op_repr <- function(self) sprintf("%s(%s)", .qpyname(self$type), .str(self$args[[1L]]))

.qclass("SqOperator", "Expr", op_symbol = "sq", is_commutative = FALSE,
  `__new__` = .cm(function(cls, args, ...) .sq_raw(cls, list(.qS(args[[1L]])))),
  state = .prop(function(self) self$args[[1L]]),
  is_symbolic = .prop(function(self) self$args[[1L]]$type != "Integer"),
  `__str__` = .sq_op_str,
  `__repr__` = function(self) NULL,
  apply_operator = function(self, state) .qni("implement apply_operator in a subclass"),
  `_sympystr` = .sq_basic_str,
  `_sympyrepr` = .qbasic_printers[[2L]],
  `_latex` = .sq_basic_str)
.qclass("BosonicOperator", "SqOperator")
.qclass("Annihilator", "SqOperator")
.qclass("Creator", "SqOperator")

.sq_lat <- function(fmt0, fmt) function(self, p) {
  st <- self$args[[1L]]
  if (st$type == "Integer" && num_is_zero(st)) fmt0 else sprintf(fmt, .latex(st))
}

.qclass("AnnihilateBoson", c("BosonicOperator", "Annihilator"), op_symbol = "b",
  `_dagger_` = function(self) CreateBoson(self$args[[1L]]),
  apply_operator = function(self, state) {
    if (!.qa(self, "is_symbolic") && .qisa(state, "FockStateKet")) {
      el <- self$args[[1L]]
      return(.mul(list(sqrt(.qmeth(state, "__getitem__", el)), .qmeth(state, "down", el))))
    }
    .mul(list(self, state))
  },
  `__repr__` = .sq_op_repr,
  `_latex` = .sq_lat("b_{0}", "b_{%s}"))

.qclass("CreateBoson", c("BosonicOperator", "Creator"), op_symbol = "b+",
  `_dagger_` = function(self) AnnihilateBoson(self$args[[1L]]),
  apply_operator = function(self, state) {
    if (!.qa(self, "is_symbolic") && .qisa(state, "FockStateKet")) {
      el <- self$args[[1L]]
      return(.mul(list(sqrt(.add(list(.qmeth(state, "__getitem__", el), S_One))), .qmeth(state, "up", el))))
    }
    .mul(list(self, state))
  },
  `__repr__` = .sq_op_repr,
  `_latex` = .sq_lat("{b^\\dagger_{0}}", "{b^\\dagger_{%s}}"))

#' @rdname secondquant
#' @export
SqOperator <- function(k) .qfunc("SqOperator", list(k))
#' @rdname secondquant
#' @export
BosonicOperator <- function(k) .qfunc("BosonicOperator", list(k))
#' @rdname secondquant
#' @export
AnnihilateBoson <- function(k) .qfunc("AnnihilateBoson", list(k))
#' @rdname secondquant
#' @export
CreateBoson <- function(k) .qfunc("CreateBoson", list(k))
#' @rdname secondquant
#' @export
B <- AnnihilateBoson
#' @rdname secondquant
#' @export
Bd <- CreateBoson

.sq_fermi_attrs <- list(
  is_restricted = .prop(function(self) {
    s <- self$args[[1L]]
    if (.sq_af(s, "below_fermi")) return(-1L)
    if (.sq_af(s, "above_fermi")) return(1L)
    0L
  }),
  is_above_fermi = .prop(function(self) !.sq_af(self$args[[1L]], "below_fermi")),
  is_below_fermi = .prop(function(self) !.sq_af(self$args[[1L]], "above_fermi")),
  is_only_below_fermi = .prop(function(self) .qa(self, "is_below_fermi") && !.qa(self, "is_above_fermi")),
  is_only_above_fermi = .prop(function(self) .qa(self, "is_above_fermi") && !.qa(self, "is_below_fermi")),
  `_sortkey` = function(self) {
    h <- .sq_tb(self)
    label <- .str(self$args[[1L]])
    if (.qa(self, "is_only_q_creator")) return(list(1, label, h))
    if (.qa(self, "is_only_q_annihilator")) return(list(4, label, h))
    if (.qisa(self, "Annihilator")) return(list(3, label, h))
    if (.qisa(self, "Creator")) return(list(2, label, h))
    NULL
  })
do.call(.qclass, c(list("FermionicOperator", "SqOperator"), .sq_fermi_attrs))

# the operator applied to a fermion ket, or to a product led by one
.sq_fermion_apply <- function(self, state, move) {
  if (.qisa(state, "FockStateFermionKet")) return(.qmeth(state, move, self$args[[1L]]))
  if (state$type == "Mul") {
    cn <- .args_cnc(state)
    if (length(cn$nc) && .qisa(cn$nc[[1L]], "FockStateFermionKet"))
      return(.mul(c(cn$c, list(.qmeth(cn$nc[[1L]], move, self$args[[1L]])), cn$nc[-1L])))
  }
  .mul(list(self, state))
}

.qclass("AnnihilateFermion", c("FermionicOperator", "Annihilator"), op_symbol = "f",
  `_dagger_` = function(self) CreateFermion(self$args[[1L]]),
  apply_operator = function(self, state) .sq_fermion_apply(self, state, "down"),
  is_q_creator = .prop(function(self) if (.qa(self, "is_below_fermi")) -1L else 0L),
  is_q_annihilator = .prop(function(self) if (.qa(self, "is_above_fermi")) 1L else 0L),
  is_only_q_creator = .prop(function(self) .qa(self, "is_only_below_fermi")),
  is_only_q_annihilator = .prop(function(self) .qa(self, "is_only_above_fermi")),
  `__repr__` = .sq_op_repr,
  `_latex` = .sq_lat("a_{0}", "a_{%s}"))

.qclass("CreateFermion", c("FermionicOperator", "Creator"), op_symbol = "f+",
  `_dagger_` = function(self) AnnihilateFermion(self$args[[1L]]),
  apply_operator = function(self, state) .sq_fermion_apply(self, state, "up"),
  is_q_creator = .prop(function(self) if (.qa(self, "is_above_fermi")) 1L else 0L),
  is_q_annihilator = .prop(function(self) if (.qa(self, "is_below_fermi")) -1L else 0L),
  is_only_q_creator = .prop(function(self) .qa(self, "is_only_above_fermi")),
  is_only_q_annihilator = .prop(function(self) .qa(self, "is_only_below_fermi")),
  `__repr__` = .sq_op_repr,
  `_latex` = .sq_lat("{a^\\dagger_{0}}", "{a^\\dagger_{%s}}"))

#' @rdname secondquant
#' @export
FermionicOperator <- function(k) .qfunc("FermionicOperator", list(k))
#' @rdname secondquant
#' @export
AnnihilateFermion <- function(k) .qfunc("AnnihilateFermion", list(k))
#' @rdname secondquant
#' @export
CreateFermion <- function(k) .qfunc("CreateFermion", list(k))
#' @rdname secondquant
#' @export
sym_F <- AnnihilateFermion
#' @rdname secondquant
#' @export
Fd <- CreateFermion

# ---------------------------------------------------------------------------
# Fock states

# FermionState.fermi_level: a class attribute that FermionState.__new__
# overwrites on the class it builds
.sq_fermi <- new.env()
.sq_fl <- function(cls) {
  for (c in .QC[[cls]]$mro) { v <- .sq_fermi[[c]]; if (!is.null(v)) return(v) }
  0
}

.sq_labels_str <- function(l) .str(do.call(Tuple, l))

.qclass("FockState", "Expr", is_commutative = FALSE, lbracket = "", rbracket = "",
  lbracket_latex = "", rbracket_latex = "",
  `__new__` = .cm(function(cls, args, ...) .sq_raw(cls, list(do.call(Tuple, .sq_seq(args[[1L]]))))),
  `__getitem__` = function(self, i) self$args[[1L]]$args[[.qint(i) + 1L]],
  `__len__` = function(self) length(self$args[[1L]]$args),
  `_labels` = function(self) self$args[[1L]]$args,
  `__str__` = function(self) paste0(.qa(self, "lbracket"), .sq_labels_str(.qmeth(self, "_labels")), .qa(self, "rbracket")),
  `__repr__` = function(self) sprintf("FockState(%s)", .str(self$args[[1L]])),
  `_sympystr` = .sq_basic_str,
  `_sympyrepr` = .qbasic_printers[[2L]],
  `_latex` = function(self, p) paste0(.qa(self, "lbracket_latex"), .latex(do.call(Tuple, .qmeth(self, "_labels"))),
                                      .qa(self, "rbracket_latex")))

.qclass("BosonState", "FockState",
  up = function(self, i) {
    i <- .qint(i) + 1L
    occ <- self$args[[1L]]$args
    occ[[i]] <- .add(list(occ[[i]], S_One))
    .qfunc(self$type, list(occ))
  },
  down = function(self, i) {
    i <- .qint(i) + 1L
    occ <- self$args[[1L]]$args
    if (.qeq(occ[[i]], S_Zero)) return(S_Zero)
    occ[[i]] <- .sub(occ[[i]], S_One)
    .qfunc(self$type, list(occ))
  })

.sq_in <- function(x, lst) any(vapply(lst, function(y) .qeq(x, y), NA))

.qclass("FermionState", "FockState",
  `__new__` = .cm(function(cls, args, fermi_level = 0, ...) {
    occ <- .sq_seq(args[[1L]])
    if (length(args) > 1L) fermi_level <- args[[2L]]
    if (is_sym(fermi_level)) fermi_level <- .qint(fermi_level)
    sign <- 0L
    if (length(occ) > 1L) {
      r <- tryCatch(.sq_sort_anticommuting_fermions(occ, .sq_key_index),
                    ViolationOfPauliPrinciple = function(e) NULL)
      if (is.null(r)) return(S_Zero)
      occ <- r[[1L]]; sign <- r[[2L]]
    }
    assign(cls, fermi_level, envir = .sq_fermi)
    if (.qmeth(cls, "_count_holes", occ) > fermi_level) return(S_Zero)
    obj <- .sq_raw(cls, list(do.call(Tuple, occ)))
    if (sign %% 2) .neg(obj) else obj
  }),
  fermi_level = .prop(function(self) .sq_fl(self$type)),
  up = function(self, i) {
    present <- .sq_in(i, self$args[[1L]]$args)
    if (.qmeth(self, "_only_above_fermi", i)) {
      if (present) return(S_Zero)
      return(.qmeth(self, "_add_orbit", i))
    }
    if (.qmeth(self, "_only_below_fermi", i)) {
      if (present) return(.qmeth(self, "_remove_orbit", i))
      return(S_Zero)
    }
    if (present) {
      hole <- Dummy("i", below_fermi = TRUE)
      return(.mul(list(KroneckerDelta(i, hole), .qmeth(self, "_remove_orbit", i))))
    }
    particle <- Dummy("a", above_fermi = TRUE)
    .mul(list(KroneckerDelta(i, particle), .qmeth(self, "_add_orbit", i)))
  },
  down = function(self, i) {
    present <- .sq_in(i, self$args[[1L]]$args)
    if (.qmeth(self, "_only_above_fermi", i)) {
      if (present) return(.qmeth(self, "_remove_orbit", i))
      return(S_Zero)
    }
    if (.qmeth(self, "_only_below_fermi", i)) {
      if (present) return(S_Zero)
      return(.qmeth(self, "_add_orbit", i))
    }
    if (present) {
      hole <- Dummy("i", below_fermi = TRUE)
      return(.mul(list(KroneckerDelta(i, hole), .qmeth(self, "_add_orbit", i))))
    }
    particle <- Dummy("a", above_fermi = TRUE)
    .mul(list(KroneckerDelta(i, particle), .qmeth(self, "_remove_orbit", i)))
  },
  `_only_below_fermi` = .cm(function(cls, i) {
    if (is_number_type(i)) return(.qtruthy(.le_num(i, .sq_fl(cls))))
    .sq_af(i, "below_fermi")
  }),
  `_only_above_fermi` = .cm(function(cls, i) {
    if (is_number_type(i)) return(!.qtruthy(.le_num(i, .sq_fl(cls))))
    if (.sq_af(i, "above_fermi")) return(TRUE)
    .sq_fl(cls) == 0
  }),
  `_remove_orbit` = function(self, i) {
    occ <- self$args[[1L]]$args
    pos <- which(vapply(occ, function(y) .qeq(i, y), NA))[1L] - 1L
    occ <- occ[-(pos + 1L)]
    r <- .qfunc(self$type, list(occ), fermi_level = .qa(self, "fermi_level"))
    if (pos %% 2) .mul(list(S_NegativeOne, r)) else r
  },
  `_add_orbit` = function(self, i) .qfunc(self$type, list(c(list(i), self$args[[1L]]$args)), fermi_level = .qa(self, "fermi_level")),
  `_count_holes` = .cm(function(cls, occupations) sum(vapply(occupations, function(i) isTRUE(.qmeth(cls, "_only_below_fermi", i)), NA))),
  `_negate_holes` = function(self, occupations) lapply(occupations, function(i)
    if (isTRUE(.qmeth(self, "_only_below_fermi", i)) && is_number_type(i)) .neg(i) else i),
  `__repr__` = function(self) {
    fl <- .qa(self, "fermi_level")
    if (fl != 0) sprintf("FockStateKet(%s, fermi_level=%s)", .str(self$args[[1L]]), fl)
    else sprintf("FockStateKet(%s)", .str(self$args[[1L]]))
  },
  `_labels` = function(self) .qmeth(self, "_negate_holes", self$args[[1L]]$args))

# i <= fermi_level for a number i
.le_num <- function(i, fl) if (num_cmp(i, as_sym(fl)) <= 0) S_true else S_false

.qclass("FockStateKet", "FockState", lbracket = "|", rbracket = ">",
  lbracket_latex = "\\left|", rbracket_latex = "\\right\\rangle")
.qclass("FockStateBra", "FockState", lbracket = "<", rbracket = "|",
  lbracket_latex = "\\left\\langle", rbracket_latex = "\\right|")
.qclass("FockStateBosonKet", c("BosonState", "FockStateKet"),
  `_dagger_` = function(self) .qfunc("FockStateBosonBra", self$args))
.qclass("FockStateBosonBra", c("BosonState", "FockStateBra"),
  `_dagger_` = function(self) .qfunc("FockStateBosonKet", self$args))
.qclass("FockStateFermionKet", c("FermionState", "FockStateKet"),
  `_dagger_` = function(self) .qfunc("FockStateFermionBra", self$args))
.qclass("FockStateFermionBra", c("FermionState", "FockStateBra"),
  `_dagger_` = function(self) .qfunc("FockStateFermionKet", self$args))

#' @rdname secondquant
#' @export
FockState <- function(occupations) .qfunc("FockState", list(occupations))
#' @rdname secondquant
#' @export
BosonState <- function(occupations) .qfunc("BosonState", list(occupations))
#' @rdname secondquant
#' @export
FermionState <- function(occupations, fermi_level = 0) .qfunc("FermionState", list(occupations), fermi_level = fermi_level)
#' @rdname secondquant
#' @export
FockStateKet <- function(occupations) .qfunc("FockStateKet", list(occupations))
#' @rdname secondquant
#' @export
FockStateBra <- function(occupations) .qfunc("FockStateBra", list(occupations))
#' @rdname secondquant
#' @export
FockStateBosonKet <- function(occupations) .qfunc("FockStateBosonKet", list(occupations))
#' @rdname secondquant
#' @export
FockStateBosonBra <- function(occupations) .qfunc("FockStateBosonBra", list(occupations))
#' @rdname secondquant
#' @export
FockStateFermionKet <- function(occupations, fermi_level = 0) .qfunc("FockStateFermionKet", list(occupations), fermi_level = fermi_level)
#' @rdname secondquant
#' @export
FockStateFermionBra <- function(occupations, fermi_level = 0) .qfunc("FockStateFermionBra", list(occupations), fermi_level = fermi_level)
#' @rdname secondquant
#' @export
BKet <- FockStateBosonKet
#' @rdname secondquant
#' @export
BBra <- FockStateBosonBra
#' @rdname secondquant
#' @export
FKet <- FockStateFermionKet
#' @rdname secondquant
#' @export
FBra <- FockStateFermionBra

# FockStateBra.__mul__: a bra times a ket is their inner product
local({
  old <- .ops_hooks$dispatch
  .ops_hooks$dispatch <- function(gen, e1, e2) {
    if (gen == "*" && is_sym(e1) && is_sym(e2) && .qisa(e1, "FockStateBra") && .qisa(e2, "FockStateKet"))
      return(sq_InnerProduct(e1, e2))
    if (is.null(old)) NULL else old(gen, e1, e2)
  }
})

# ---------------------------------------------------------------------------
# Applying operators

.sq_apply_mul <- function(m) {
  if (m$type != "Mul") return(m)
  cn <- .args_cnc(m)
  c_part <- cn$c; nc <- cn$nc
  n <- length(nc)
  if (n < 2L) return(m)
  last <- nc[[n]]; ntl <- nc[[n - 1L]]
  if (!.qisa(last, "FockStateKet")) return(m)
  rest <- c(c_part, nc[seq_len(n - 2L)])
  if (.qisa(ntl, "SqOperator")) {
    if (.qa(ntl, "is_symbolic")) return(m)
    result <- .qmeth(ntl, "apply_operator", last)
    if (.qeq(result, S_Zero)) return(S_Zero)
    return(.sq_apply_mul(.mul(c(rest, list(result)))))
  }
  if (ntl$type == "Pow") {
    base <- ntl$args[[1L]]; ex <- ntl$args[[2L]]
    if (.qisa(base, "SqOperator") && ex$type == "Integer") {
      if (.qa(base, "is_symbolic")) return(m)
      result <- last
      for (k in seq_len(.qint(ex))) {
        result <- .qmeth(base, "apply_operator", result)
        if (.qeq(result, S_Zero)) break
      }
      if (.qeq(result, S_Zero)) return(S_Zero)
      return(.sq_apply_mul(.mul(c(rest, list(result)))))
    }
    return(m)
  }
  if (.qisa(ntl, "FockStateBra")) {
    result <- sq_InnerProduct(ntl, last)
    if (.qeq(result, S_Zero)) return(S_Zero)
    return(.sq_apply_mul(.mul(c(rest, list(result)))))
  }
  m
}

#' @rdname secondquant
#' @export
apply_operators <- function(e) {
  e <- expand(as_sym(e))
  muls <- atoms(e, "Mul")
  subs(e, lapply(muls, function(m) list(m, .sq_apply_mul(m))))
}

#' @rdname secondquant
#' @export
sq_InnerProduct <- function(bra, ket) {
  if (!.qisa(bra, "FockStateBra")) .qstop("TypeError", "must be a bra")
  if (!.qisa(ket, "FockStateKet")) .qstop("TypeError", "must be a ket")
  result <- S_One
  bs <- bra$args[[1L]]$args; ks <- ket$args[[1L]]$args
  for (k in seq_len(min(length(bs), length(ks)))) {
    result <- .mul(list(result, KroneckerDelta(bs[[k]], ks[[k]])))
    if (.qeq(result, S_Zero)) break
  }
  result
}

#' @rdname secondquant
#' @export
matrix_rep <- function(op, basis) {
  n <- .qlen_basis(basis)
  a <- zeros(n)
  for (i in seq_len(n)) for (j in seq_len(n))
    a[i, j] <- apply_operators(.mul(list(sq_Dagger(basis$basis[[i]]), op, basis$basis[[j]])))
  a
}

.qlen_basis <- function(b) length(b$basis)

# the basis with its methods index(state) (0 based), state(i) and b[i]
.sq_basis <- function(cls, basis, ...) {
  b <- list(basis = basis, n_basis = length(basis), ...)
  b$index <- function(state) {
    k <- which(vapply(basis, function(s) .qeq(s, state), NA))
    if (!length(k)) .qstop("ValueError", .qmeth(state, "__repr__"), " is not in list")
    k[1L] - 1L
  }
  b$state <- function(i) basis[[.qint(i) + 1L]]
  b$`__getitem__` <- b$state
  class(b) <- cls
  b
}

#' @rdname secondquant
#' @export
VarBosonicBasis <- function(n_max) {
  .sq_basis("VarBosonicBasis", lapply(seq_len(.qint(n_max)) - 1L, function(i) FockStateBosonKet(list(i))), n_max = n_max)
}

#' @rdname secondquant
#' @export
FixedBosonicBasis <- function(n_particles, n_levels) {
  # the nested loops for i0 in range(n_levels) for i1 in range(i0 + 1) ...
  np <- .qint(n_particles); nl <- .qint(n_levels)
  locs <- list()
  rec <- function(prefix, k) {
    if (k > np) { locs[[length(locs) + 1L]] <<- prefix; return(invisible()) }
    top <- if (k == 1L) nl else prefix[k - 1L] + 1L
    for (v in seq_len(top) - 1L) rec(c(prefix, v), k + 1L)
  }
  rec(integer(), 1L)
  basis <- lapply(locs, function(t) {
    occ <- integer(nl)
    for (lv in t) occ[lv + 1L] <- occ[lv + 1L] + 1L
    FockStateBosonKet(as.list(occ))
  })
  .sq_basis("FixedBosonicBasis", basis, n_particles = n_particles, n_levels = n_levels,
            particle_locations = locs)
}

.sq_basis_repr <- function(x) paste0("[", paste(vapply(x$basis, function(s) .qmeth(s, "__repr__"), ""), collapse = ", "), "]")

#' @export
format.VarBosonicBasis <- function(x, ...) .sq_basis_repr(x)
#' @export
format.FixedBosonicBasis <- function(x, ...) .sq_basis_repr(x)
#' @export
print.VarBosonicBasis <- function(x, ...) { cat(format(x), "\n", sep = ""); invisible(x) }
#' @export
print.FixedBosonicBasis <- function(x, ...) { cat(format(x), "\n", sep = ""); invisible(x) }
#' @export
length.VarBosonicBasis <- function(x) length(unclass(x)$basis)
#' @export
length.FixedBosonicBasis <- function(x) length(unclass(x)$basis)
