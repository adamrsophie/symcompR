# Second quantization: normal ordering, contractions, Wick's theorem and the
# handling of dummy indices, after sympy/physics/secondquant.py.

# ---------------------------------------------------------------------------
# Sorting strings of anticommuting operators

.sq_key_operator <- function(op) .qmeth(op, "_sortkey")

.sq_key_index <- function(index) {
  h <- .sq_tb(index)
  label <- .str(index)
  base <- if (index$type == "Dummy") 20 else 10
  if (.sq_af(index, "above_fermi")) return(list(base, label, h))
  if (.sq_af(index, "below_fermi")) return(list(base + 1, label, h))
  list(base + 2, label, h)
}

.sq_key_str <- function(k) paste0("(", k[[1L]], ", '", k[[2L]], "', ", k[[3L]], ")")

# _sort_anticommuting_fermions: list(sorted, number of swaps)
.sq_sort_anticommuting_fermions <- function(string1, key = .sq_key_operator) {
  verified <- FALSE
  sign <- 0L
  n <- length(string1)
  rng <- seq_len(max(n - 1L, 0L))
  rev <- if (n >= 3L) seq.int(n - 2L, 1L) else integer()
  keys <- lapply(string1, key)
  pos <- seq_len(n)
  pass <- function(ix) {
    for (i in ix) {
      cmp <- .key_cmp(keys[[i]], keys[[i + 1L]])
      if (cmp == 0) .sq_stop("ViolationOfPauliPrinciple", "[", .sq_key_str(keys[[i]]), ", ", .sq_key_str(keys[[i + 1L]]), "]")
      if (cmp > 0) {
        verified <<- FALSE
        keys[c(i, i + 1L)] <<- keys[c(i + 1L, i)]
        pos[c(i, i + 1L)] <<- pos[c(i + 1L, i)]
        sign <<- sign + 1L
      }
    }
  }
  while (!verified) {
    verified <- TRUE
    pass(rng)
    if (verified) break
    pass(rev)
  }
  list(string1[pos], sign)
}

# ---------------------------------------------------------------------------
# Commutator

.qclass("sqCommutator", "Function", pyname = "Commutator", is_commutative = FALSE,
  eval = .cm(function(cls, a, b) {
    if (!(.qtruthy(a) && .qtruthy(b))) return(S_Zero)
    if (.qeq(a, b)) return(S_Zero)
    if (.is_comm(a) || .is_comm(b)) return(S_Zero)
    a <- expand(a)
    if (a$type == "Add") return(.add(lapply(.sympy_args(a), function(t) sq_Commutator(t, b))))
    b <- expand(b)
    if (b$type == "Add") return(.add(lapply(.sympy_args(b), function(t) sq_Commutator(a, t))))
    ca <- .args_cnc(a); cb <- .args_cnc(b)
    c_part <- c(ca$c, cb$c)
    if (length(c_part))
      return(.mul(list(.mul(c_part), sq_Commutator(.qmul_raw(ca$nc), .qmul_raw(cb$nc)))))
    if (.qisa(a, "BosonicOperator") && .qisa(b, "BosonicOperator")) {
      if (.qisa(b, "CreateBoson") && .qisa(a, "AnnihilateBoson"))
        return(KroneckerDelta(a$args[[1L]], b$args[[1L]]))
      if (.qisa(a, "CreateBoson") && .qisa(b, "AnnihilateBoson"))
        return(.mul(list(S_NegativeOne, KroneckerDelta(a$args[[1L]], b$args[[1L]]))))
      return(S_Zero)
    }
    if (.qisa(a, "FermionicOperator") && .qisa(b, "FermionicOperator"))
      return(.sub(wicks(.mul(list(a, b))), wicks(.mul(list(b, a)))))
    if (.key_cmp(.sort_key(a), .sort_key(b)) > 0) return(.mul(list(S_NegativeOne, sq_Commutator(b, a))))
    NULL
  }),
  `__new__` = .cm(function(cls, args, ...) {
    a <- .qS(args[[1L]]); b <- .qS(args[[2L]])
    r <- .qmeth(cls, "eval", a, b)
    if (!is.null(r)) return(r)
    .sq_raw(cls, list(a, b))
  }),
  doit = function(self, wicks = FALSE, ...) {
    a <- self$args[[1L]]; b <- self$args[[2L]]
    if (isTRUE(wicks)) {
      a <- .sq_doit(a, wicks = wicks, ...); b <- .sq_doit(b, wicks = wicks, ...)
      r <- tryCatch(.sub(wicks(.mul(list(a, b))), wicks(.mul(list(b, a)))),
                    ContractionAppliesOnlyToFermions = function(e) NULL,
                    WicksTheoremDoesNotApply = function(e) NULL)
      if (!is.null(r)) return(r)
    }
    .sq_doit(.sub(.mul(list(a, b)), .mul(list(b, a))), wicks = wicks, ...)
  },
  `__repr__` = function(self) sprintf("Commutator(%s,%s)", .str_top(self$args[[1L]]), .str_top(self$args[[2L]])),
  `__str__` = function(self) sprintf("[%s,%s]", .str_top(self$args[[1L]]), .str_top(self$args[[2L]])),
  `_sympystr` = .sq_basic_str,
  `_sympyrepr` = .qbasic_printers[[2L]],
  `_latex` = function(self, p) sprintf("\\left[%s,%s\\right]", .latex(self$args[[1L]]), .latex(self$args[[2L]])))

#' @rdname secondquant
#' @export
sq_Commutator <- function(a, b) .qfunc("sqCommutator", list(a, b))

# Commutator.doit and NO.doit override Basic.doit: the arguments are not
# evaluated first
register_doit("sqCommutator", structure(function(e, deep = TRUE, ...) .qmeth(e, "doit", ...), shallow = TRUE))

# expr.doit(**hints) passing the hints to the secondquant objects
.sq_doit <- function(e, ...) {
  if (!is_sym(e) || !length(e$args)) return(e)
  if (!is.null(.qnode_def(e$type)) && .qhas(e, "doit")) return(.qmeth(e, "doit", ...))
  args <- lapply(e$args, .sq_doit, ...)
  if (identical(vapply(args, .k, ""), vapply(e$args, .k, ""))) {
    if (e$type %in% c("Add", "Mul", "Pow", "Tuple") || !is.null(.qnode_def(e$type))) return(e)
    return(doit(e))
  }
  r <- sym_rebuild(e, args)
  if (is_sym(r) && !(r$type %in% c("Add", "Mul", "Pow", "Tuple")) && is.null(.qnode_def(r$type))) doit(r) else r
}

# ---------------------------------------------------------------------------
# Normal ordering

.qclass("NO", "Expr", is_commutative = FALSE,
  `__new__` = .cm(function(cls, args, ...) {
    arg <- expand(.qS(args[[1L]]))
    if (arg$type == "Add") return(.add(lapply(.sympy_args(arg), NO)))
    if (arg$type == "Mul") {
      cn <- .args_cnc(arg)
      seq <- cn$nc
      if (length(cn$c)) {
        coeff <- .mul(cn$c)
        if (!length(seq)) return(coeff)
      } else coeff <- S_One
      newseq <- list(); foundit <- FALSE
      for (fac in seq) {
        if (.qisa(fac, "NO")) { newseq <- c(newseq, fac$args); foundit <- TRUE }
        else newseq <- c(newseq, list(fac))
      }
      if (foundit) return(.mul(list(coeff, NO(.mul(newseq)))))
      if (.qisa(seq[[1L]], "BosonicOperator")) .qni("")
      r <- tryCatch(.sq_sort_anticommuting_fermions(seq), ViolationOfPauliPrinciple = function(e) NULL)
      if (is.null(r)) return(S_Zero)
      newseq <- r[[1L]]; sign <- r[[2L]]
      if (sign %% 2) return(.mul(list(S_NegativeOne, coeff, NO(.mul(newseq)))))
      if (sign) return(.mul(list(coeff, NO(.mul(newseq)))))
      if (!.qeq(coeff, S_One)) return(.mul(list(coeff, NO(.mul(newseq)))))
      return(.sq_raw(cls, list(.mul(newseq))))
    }
    arg
  }),
  has_q_creators = .prop(function(self) .qa(.sympy_args(self$args[[1L]])[[1L]], "is_q_creator")),
  has_q_annihilators = .prop(function(self) { a <- .sympy_args(self$args[[1L]]); .qa(a[[length(a)]], "is_q_annihilator") }),
  doit = function(self, remove_brackets = TRUE, ...) {
    if (isTRUE(remove_brackets)) return(.qmeth(self, "_remove_brackets"))
    NO(.sq_doit(self$args[[1L]], remove_brackets = remove_brackets, ...))
  },
  `_remove_brackets` = function(self) {
    subslist <- list()
    for (i in .qmeth(self, "iter_q_creators")) {
      op <- .qmeth(self, "__getitem__", i)
      if (.qtruthy(.qa(op, "is_q_annihilator"))) {
        st <- op$args[[1L]]
        if (st$type == "Dummy") {
          assume <- st$extra$other
          assume <- assume[setdiff(names(assume), "above_fermi")]
          assume[["below_fermi"]] <- TRUE
          below <- do.call(Dummy, c(list("i"), as.list(assume)))
          assume <- assume[setdiff(names(assume), "below_fermi")]
          assume[["above_fermi"]] <- TRUE
          above <- do.call(Dummy, c(list("a"), as.list(assume)))
          split <- .add(list(.mul(list(.qfunc(op$type, list(below)), KroneckerDelta(below, st))),
                             .mul(list(.qfunc(op$type, list(above)), KroneckerDelta(above, st)))))
          subslist[[length(subslist) + 1L]] <- list(op, split)
        } else {
          .sq_stop("SubstitutionOfAmbigousOperatorFailed", .str(op))
        }
      }
    }
    if (length(subslist)) {
      result <- NO(subs(self, subslist))
      if (result$type == "Add") return(.add(lapply(.sympy_args(result), function(t) .sq_doit(t))))
      return(NULL)
    }
    self$args[[1L]]
  },
  `__getitem__` = function(self, i) .sympy_args(self$args[[1L]])[[.qint(i) + 1L]],
  `__len__` = function(self) length(self$args[[1L]]$args),
  iter_q_annihilators = function(self) {
    ops <- .sympy_args(self$args[[1L]])
    out <- integer()
    for (i in rev(seq_along(ops))) {
      if (.qtruthy(.qa(ops[[i]], "is_q_annihilator"))) out <- c(out, i - 1L) else break
    }
    out
  },
  iter_q_creators = function(self) {
    ops <- .sympy_args(self$args[[1L]])
    out <- integer()
    for (i in seq_along(ops)) {
      if (.qtruthy(.qa(ops[[i]], "is_q_creator"))) out <- c(out, i - 1L) else break
    }
    out
  },
  get_subNO = function(self, i) {
    a <- .sympy_args(self$args[[1L]])
    NO(.qmul_raw(a[-(.qint(i) + 1L)]))
  },
  `_latex` = function(self, p) sprintf("\\left\\{%s\\right\\}", .latex(self$args[[1L]])),
  `__repr__` = function(self) sprintf("NO(%s)", .str_top(self$args[[1L]])),
  `__str__` = function(self) sprintf(":%s:", .str_top(self$args[[1L]])),
  `_sympystr` = .sq_basic_str,
  `_sympyrepr` = .qbasic_printers[[2L]])

#' Normal ordering, contractions and Wick's theorem
#'
#' `NO()` is the normal ordering of a product of fermion operators.
#' `contraction()` is the contraction of two fermion operators, `wicks()`
#' expands a product of fermion operators by Wick's theorem,
#' `evaluate_deltas()` removes the Kronecker deltas of a product by
#' substitution, `substitute_dummies()` renames the dummy indices of an
#' expression into a canonical form, `PermutationOperator()` is the operator
#' exchanging two indices and `simplify_index_permutations()` gathers the
#' terms related by such permutations, as `sympy.physics.secondquant`.
#'
#' @param arg,expr,e,a,b Expressions or operators.
#' @param i,j Indices.
#' @param simplify_kronecker_deltas,expand,simplify_dummies,keep_only_fully_contracted
#'   Options of `wicks()`.
#' @param new_indices Whether `substitute_dummies()` makes new dummies.
#' @param pretty_indices A list with the letters `above`, `below` and
#'   `general` for the new dummies.
#' @param permutation_operators A list of `PermutationOperator()`s.
#' @return An expression.
#' @examples
#' i <- Symbol("i", below_fermi = TRUE); a <- Symbol("a", above_fermi = TRUE)
#' p <- Symbol("p"); q <- Symbol("q")
#' NO(sym_F(i) * Fd(a))
#' contraction(Fd(i), sym_F(i))
#' wicks(sym_F(p) * Fd(q))
#' @name secondquant-wicks
NULL

#' @rdname secondquant-wicks
#' @export
NO <- function(arg) .qfunc("NO", list(arg))
register_doit("NO", structure(function(e, deep = TRUE, ...) .qmeth(e, "doit", ...), shallow = TRUE))

#' @rdname secondquant-wicks
#' @export
contraction <- function(a, b) {
  if (.qisa(b, "FermionicOperator") && .qisa(a, "FermionicOperator")) {
    sa <- a$args[[1L]]; sb <- b$args[[1L]]
    if (.qisa(a, "AnnihilateFermion") && .qisa(b, "CreateFermion")) {
      if (.sq_af(sb, "below_fermi") || .sq_af(sa, "below_fermi")) return(S_Zero)
      if (.sq_af(sb, "above_fermi") || .sq_af(sa, "above_fermi")) return(KroneckerDelta(sa, sb))
      return(.mul(list(KroneckerDelta(sa, sb), KroneckerDelta(sb, Dummy("a", above_fermi = TRUE)))))
    }
    if (.qisa(b, "AnnihilateFermion") && .qisa(a, "CreateFermion")) {
      if (.sq_af(sb, "above_fermi") || .sq_af(sa, "above_fermi")) return(S_Zero)
      if (.sq_af(sb, "below_fermi") || .sq_af(sa, "below_fermi")) return(KroneckerDelta(sa, sb))
      return(.mul(list(KroneckerDelta(sa, sb), KroneckerDelta(sb, Dummy("i", below_fermi = TRUE)))))
    }
    return(S_Zero)
  }
  t <- vapply(list(a, b), function(x) .qisa(x, "FermionicOperator"), NA)
  .sq_stop("ContractionAppliesOnlyToFermions", "(", paste(ifelse(t, "True", "False"), collapse = ", "), ")")
}

# ---------------------------------------------------------------------------
# Kronecker deltas and the Fermi level

.kd_af <- function(d, k, nm) .sq_af(d$args[[k]], nm)
.kd_is_above <- function(d) !(.kd_af(d, 1L, "below_fermi") || .kd_af(d, 2L, "below_fermi"))
.kd_is_below <- function(d) !(.kd_af(d, 1L, "above_fermi") || .kd_af(d, 2L, "above_fermi"))
.kd_equal_info <- function(d) {
  if (.kd_af(d, 1L, "below_fermi") && .kd_af(d, 2L, "below_fermi")) return(TRUE)
  if (.kd_af(d, 1L, "above_fermi") && .kd_af(d, 2L, "above_fermi")) return(TRUE)
  .kd_is_below(d) && .kd_is_above(d)
}
.kd_pref <- function(d) {
  if (!.kd_is_above(d)) return(if (.kd_af(d, 1L, "below_fermi")) 0L else 1L)
  if (!.kd_is_below(d)) return(if (.kd_af(d, 1L, "above_fermi")) 0L else 1L)
  0L
}
.kd_preferred <- function(d) if (.kd_pref(d)) d$args[[2L]] else d$args[[1L]]
.kd_killable <- function(d) if (.kd_pref(d)) d$args[[1L]] else d$args[[2L]]
.is_kd <- function(e) is_sym(e) && e$type == "Function" && e$name == "KroneckerDelta"

# the Fermi level properties of KroneckerDelta
.kd_attr <- function(d, name) switch(name,
  is_above_fermi = .kd_is_above(d), is_below_fermi = .kd_is_below(d),
  is_only_above_fermi = .kd_af(d, 1L, "above_fermi") || .kd_af(d, 2L, "above_fermi"),
  is_only_below_fermi = .kd_af(d, 1L, "below_fermi") || .kd_af(d, 2L, "below_fermi"),
  indices_contain_equal_information = .kd_equal_info(d),
  preferred_index = .kd_preferred(d), killable_index = .kd_killable(d), NULL)

#' @rdname secondquant-wicks
#' @export
evaluate_deltas <- function(e) {
  e <- as_sym(e)
  if (e$type == "Add") return(.add(lapply(.sympy_args(e), evaluate_deltas)))
  if (e$type != "Mul") return(e)
  deltas <- list()
  cnt <- list()
  for (f in .sympy_args(e)) {
    for (s in free_symbols(f)) cnt[[s$key]] <- if (is.null(cnt[[s$key]])) 0L else cnt[[s$key]] + 1L
    if (.is_kd(f)) deltas[[length(deltas) + 1L]] <- f
  }
  issym <- function(x) x$type %in% c("Symbol", "Dummy", "Wild")
  for (d in deltas) {
    kil <- .kd_killable(d); pre <- .kd_preferred(d)
    if (issym(kil) && isTRUE(cnt[[kil$key]] > 0L)) {
      e <- subs(e, kil, pre)
      if (length(deltas) > 1L) return(evaluate_deltas(e))
    } else if (issym(pre) && isTRUE(cnt[[pre$key]] > 0L) && .kd_equal_info(d)) {
      e <- subs(e, pre, kil)
      if (length(deltas) > 1L) return(evaluate_deltas(e))
    }
  }
  e
}

# ---------------------------------------------------------------------------
# Dummy indices

.sq_atoms_dummy <- function(e) atoms(e, "Dummy")
.sq_make_args <- function(e, type) if (e$type == type) .sympy_args(e) else list(e)
.sq_has_dups <- function(x) anyDuplicated(x) > 0L
.sq_sort_by <- function(items, keys) .key_sort(keys, items)

#' @rdname secondquant-wicks
#' @export
substitute_dummies <- function(expr, new_indices = FALSE, pretty_indices = list()) {
  expr <- as_sym(expr)
  if (new_indices) {
    la <- if (is.null(pretty_indices$above)) "" else pretty_indices$above
    lb <- if (is.null(pretty_indices$below)) "" else pretty_indices$below
    lg <- if (is.null(pretty_indices$general)) "" else pretty_indices$general
    pick <- function(letters, pre, n) if (n < nchar(letters)) substr(letters, n + 1L, n + 1L) else paste0(pre, "_", n - nchar(letters))
  }
  aboves <- list(); belows <- list(); generals <- list()
  dummies <- .sq_atoms_dummy(expr)
  dummies <- .sq_sort_by(dummies, lapply(dummies, .sort_key))
  na <- ni <- np <- 0L
  for (d in dummies) {
    assum <- d$extra$other
    if (.sq_af(d, "above_fermi")) {
      if (new_indices) { sym <- pick(la, "a", na); na <- na + 1L }
      which <- "a"
    } else if (.sq_af(d, "below_fermi")) {
      if (new_indices) { sym <- pick(lb, "i", ni); ni <- ni + 1L }
      which <- "i"
    } else {
      if (new_indices) { sym <- pick(lg, "p", np); np <- np + 1L }
      which <- "p"
    }
    nd <- if (new_indices) do.call(Dummy, c(list(sym), as.list(assum))) else d
    if (which == "a") aboves[[length(aboves) + 1L]] <- nd
    else if (which == "i") belows[[length(belows) + 1L]] <- nd
    else generals[[length(generals) + 1L]] <- nd
  }
  expr <- expand(expr)
  new_terms <- list()
  for (term in .sq_make_args(expr, "Add")) {
    ki <- 0L; ka <- 0L; kp <- 0L
    nxt <- function(l, k) {
      if (k >= length(l)) .qstop("StopIteration", "")
      l[[k + 1L]]
    }
    ordered <- .sq_get_ordered_dummies(term)
    sk <- character(); sv <- list(); skey <- list()
    for (d in ordered) {
      if (.sq_af(d, "below_fermi")) { v <- nxt(belows, ki); ki <- ki + 1L }
      else if (.sq_af(d, "above_fermi")) { v <- nxt(aboves, ka); ka <- ka + 1L }
      else { v <- nxt(generals, kp); kp <- kp + 1L }
      if (!(d$key %in% sk)) { sk <- c(sk, d$key); skey[[length(skey) + 1L]] <- d }
      sv[[d$key]] <- v
    }
    subslist <- list(); final_subs <- list()
    for (n in seq_along(sk)) {
      k <- skey[[n]]; v <- sv[[sk[n]]]
      if (identical(k$key, v$key)) next
      if (v$key %in% sk) {
        if (sv[[v$key]]$key %in% sk) {
          x <- Dummy("x")
          subslist[[length(subslist) + 1L]] <- list(k, x)
          final_subs[[length(final_subs) + 1L]] <- list(x, v)
        } else {
          final_subs <- c(list(list(k, v)), final_subs)
        }
      } else {
        subslist[[length(subslist) + 1L]] <- list(k, v)
      }
    }
    subslist <- c(subslist, final_subs)
    new_terms[[length(new_terms) + 1L]] <- if (length(subslist)) subs(term, subslist) else term
  }
  .add(new_terms)
}

# the KeyPrinter of SymPy: a dummy prints as (name_index)
.sq_kprint <- function(e) {
  old <- .str_mode$dummy_key
  .str_mode$dummy_key <- TRUE
  on.exit(.str_mode$dummy_key <- old)
  .str(e)
}

.sq_find_all <- function(s, pat) {
  m <- gregexpr(pat, s, fixed = TRUE)[[1L]]
  if (m[1L] == -1L) integer() else as.integer(m) - 1L
}

.sq_get_ordered_dummies <- function(mul, verbose = FALSE) {
  args <- .sq_make_args(mul, "Mul")
  fac_dum <- lapply(args, .sq_atoms_dummy)
  fac_repr <- vapply(args, .sq_kprint, "")
  all_dums <- list()
  for (fd in fac_dum) for (d in fd) all_dums[[d$key]] <- d
  mask <- vapply(all_dums, function(d) if (.sq_af(d, "below_fermi")) "0" else if (.sq_af(d, "above_fermi")) "1" else "2", "")
  dum_repr <- vapply(all_dums, .sq_kprint, "")
  has_d <- function(k, d) any(vapply(fac_dum[[k]], function(x) identical(x$key, d$key), NA))
  key_of <- function(d) {
    ks <- which(vapply(seq_along(args), function(k) has_d(k, d), NA))
    other <- unique(unlist(lapply(ks, function(k) vapply(fac_dum[[k]], .k, ""))))
    other <- setdiff(other, d$key)
    masked <- fac_repr[ks]
    for (d2 in other) masked <- gsub(dum_repr[[d2]], mask[[d2]], masked, fixed = TRUE)
    all_masked <- gsub(dum_repr[[d$key]], mask[[d$key]], masked, fixed = TRUE)
    if (.sq_has_dups(all_masked)) return(list(mask[[d$key]], as.list(sort(all_masked, method = "radix"))))
    o <- order(all_masked, method = "radix")
    ks <- ks[o]; masked <- masked[o]
    all_masked <- all_masked[o]
    pos_val <- list()
    for (n in seq_along(ks)) {
      fac <- args[[ks[n]]]
      if (.qisa(fac, "AntiSymmetricTensor")) {
        if (.sq_in(d, fac$args[[2L]]$args)) pos_val[[length(pos_val) + 1L]] <- "u"
        if (.sq_in(d, fac$args[[3L]]$args)) pos_val[[length(pos_val) + 1L]] <- "l"
      } else if (.qisa(fac, "Creator")) {
        pos_val[[length(pos_val) + 1L]] <- "u"
      } else if (.qisa(fac, "Annihilator")) {
        pos_val[[length(pos_val) + 1L]] <- "l"
      } else if (.qisa(fac, "NO")) {
        for (op in .sympy_args(fac$args[[1L]])) if (has(op, d))
          pos_val[[length(pos_val) + 1L]] <- if (.qisa(op, "Creator")) "u" else "l"
      } else {
        for (fp in .sq_find_all(masked[n], dum_repr[[d$key]])) pos_val[[length(pos_val) + 1L]] <- fp
      }
    }
    list(mask[[d$key]], as.list(all_masked), pos_val[[1L]], pos_val[[length(pos_val)]])
  }
  dumkey <- lapply(all_dums, key_of)
  result <- unname(.key_sort(dumkey, all_dums))
  kstr <- vapply(dumkey, function(k) paste(deparse(k), collapse = ""), "")
  if (.sq_has_dups(kstr)) {
    grp <- split(names(all_dums), factor(kstr, levels = unique(kstr)))
    grp <- grp[lengths(grp) >= 2L]
    gkeys <- lapply(grp, function(g) dumkey[[g[1L]]])
    grp <- .key_sort(gkeys, unname(grp))
    result <- .sq_determine_ambiguous(mul, result, grp)
  }
  result
}

.sq_factory <- new.env()
.sq_factory$counter <- 0L

.sq_determine_ambiguous <- function(term, ordered, ambiguous_groups) {
  all_amb <- unlist(ambiguous_groups)
  okeys <- vapply(ordered, .k, "")
  all_ordered <- setdiff(okeys, all_amb)
  if (!length(all_ordered)) {
    g1 <- okeys[okeys %in% ambiguous_groups[[1L]]]
    all_ordered <- g1[1L]
    ambiguous_groups[[1L]] <- setdiff(ambiguous_groups[[1L]], g1[1L])
  }
  stored <- .sq_factory$counter
  subslist <- list()
  for (d in ordered) if (d$key %in% all_ordered) {
    nondum <- Symbol(sprintf("_]\"]_%d", .sq_factory$counter))
    .sq_factory$counter <- .sq_factory$counter + 1L
    subslist[[length(subslist) + 1L]] <- list(d, nondum)
  }
  newterm <- subs(term, subslist)
  neworder <- .sq_get_ordered_dummies(newterm)
  .sq_factory$counter <- stored
  for (group in ambiguous_groups) {
    og <- Filter(function(d) d$key %in% group, neworder)
    og <- rev(og)
    result <- list()
    for (d in ordered) {
      if (d$key %in% group) {
        result[[length(result) + 1L]] <- og[[length(og)]]
        og <- og[-length(og)]
      } else result[[length(result) + 1L]] <- d
    }
    ordered <- result
  }
  ordered
}

# ---------------------------------------------------------------------------
# Wick's theorem

.sq_get_contractions <- function(string1, keep_only_fully_contracted = FALSE) {
  result <- if (keep_only_fully_contracted && length(string1)) list() else list(NO(.mul(string1)))
  n <- length(string1)
  for (i in seq_len(max(n - 1L, 0L))) {
    for (j in seq.int(i + 1L, n)) {
      c <- contraction(string1[[i]], string1[[j]])
      if (.qtruthy(c)) {
        coeff <- if ((j - i + 1L) %% 2L) .mul(list(S_NegativeOne, c)) else c
        oplist <- c(string1[seq.int(i + 1L, length.out = j - i - 1L)], string1[seq.int(j + 1L, length.out = n - j)])
        head <- .mul(string1[seq_len(i - 1L)])
        if (length(oplist)) {
          result[[length(result) + 1L]] <- .mul(list(coeff, NO(.mul(list(head,
            .sq_get_contractions(oplist, keep_only_fully_contracted))))))
        } else {
          result[[length(result) + 1L]] <- .mul(list(coeff, NO(head)))
        }
      }
    }
    if (keep_only_fully_contracted) break
  }
  .add(result)
}

#' @rdname secondquant-wicks
#' @export
wicks <- function(e, simplify_kronecker_deltas = FALSE, expand = TRUE, simplify_dummies = FALSE,
                  keep_only_fully_contracted = FALSE) {
  e <- as_sym(e)
  if (!.qtruthy(e)) return(S_Zero)
  kw <- list(simplify_kronecker_deltas = simplify_kronecker_deltas, expand = expand,
             simplify_dummies = simplify_dummies, keep_only_fully_contracted = keep_only_fully_contracted)
  if (.qisa(e, c("NO", "FermionicOperator"))) return(if (keep_only_fully_contracted) S_Zero else e)
  e <- .sq_doit(e, wicks = TRUE)
  e <- symcompR::expand(e)
  if (e$type == "Add") {
    terms <- lapply(.sympy_args(e), function(t) do.call(wicks, c(list(t), kw)))
    return(if (simplify_dummies) substitute_dummies(.add(terms)) else .add(terms))
  }
  if (e$type == "Mul") {
    c_part <- list(); string1 <- list()
    for (f in .sympy_args(e)) {
      if (.is_comm(f)) c_part[[length(c_part) + 1L]] <- f else string1[[length(string1) + 1L]] <- f
    }
    n <- length(string1)
    if (n == 0L) result <- e
    else if (n == 1L) {
      if (keep_only_fully_contracted) return(S_Zero)
      result <- e
    } else {
      if (.qisa(string1[[1L]], "BosonicOperator")) .qni("")
      result <- .mul(list(.mul(c_part), .sq_get_contractions(string1, keep_only_fully_contracted)))
    }
    if (expand) result <- symcompR::expand(result)
    if (simplify_kronecker_deltas) result <- evaluate_deltas(result)
    return(result)
  }
  e
}

# ---------------------------------------------------------------------------
# Index permutations

.qclass("PermutationOperator", "Expr", is_commutative = TRUE,
  `__new__` = .cm(function(cls, args, ...) {
    ij <- lapply(args, .qS)
    ij <- .key_sort(lapply(ij, .sort_key), ij)
    .sq_craw(cls, ij)
  }),
  get_permuted = function(self, expr) {
    i <- self$args[[1L]]; j <- self$args[[2L]]
    if (has(expr, i) && has(expr, j)) {
      tmp <- Dummy()
      expr <- subs(expr, i, tmp)
      expr <- subs(expr, j, i)
      expr <- subs(expr, tmp, j)
      return(.mul(list(S_NegativeOne, expr)))
    }
    expr
  },
  `_sympystr` = .sq_basic_str,
  `_sympyrepr` = .qbasic_printers[[2L]],
  `_latex` = function(self, p) sprintf("P(%s%s)", .latex(self$args[[1L]]), .latex(self$args[[2L]])))

#' @rdname secondquant-wicks
#' @export
PermutationOperator <- function(i, j) .qfunc("PermutationOperator", list(i, j))

#' @rdname secondquant-wicks
#' @export
simplify_index_permutations <- function(expr, permutation_operators) {
  get_indices <- function(e, ind) {
    out <- list()
    for (a in e$args) {
      if (.sq_in(a, ind)) out[[length(out) + 1L]] <- a
      else if (length(a$args)) out <- c(out, get_indices(a, ind))
    }
    out
  }
  keep_one <- function(a, b, ind) {
    ka <- .sort_key(do.call(Tuple, get_indices(a, ind)))
    kb <- .sort_key(do.call(Tuple, get_indices(b, ind)))
    if (.key_cmp(kb, ka) < 0) b else a
  }
  # a set of expressions, as a list keyed by expression key
  setof <- function(l) { names(l) <- vapply(l, .k, ""); l[!duplicated(names(l))] }
  expr <- expand(as_sym(expr))
  if (expr$type != "Add") return(expr)
  terms <- setof(.sympy_args(expr))
  for (P in permutation_operators) {
    new_terms <- list(); on_hold <- list()
    while (length(terms)) {
      term <- terms[[1L]]; terms <- terms[-1L]
      permuted <- .qmeth(P, "get_permuted", term)
      pk <- permuted$key
      if (pk %in% c(names(terms), names(on_hold))) {
        if (pk %in% names(terms)) terms[[pk]] <- NULL else on_hold[[pk]] <- NULL
        keep <- keep_one(term, permuted, P$args)
        new_terms <- setof(c(new_terms, list(.mul(list(P, keep)))))
      } else {
        permuted1 <- permuted
        permuted <- substitute_dummies(permuted)
        pk <- permuted$key
        if (identical(permuted1$key, pk)) {
          on_hold <- setof(c(on_hold, list(term)))
        } else if (pk %in% c(names(terms), names(on_hold))) {
          if (pk %in% names(terms)) terms[[pk]] <- NULL else on_hold[[pk]] <- NULL
          keep <- keep_one(term, permuted, P$args)
          new_terms <- setof(c(new_terms, list(.mul(list(P, keep)))))
        } else {
          new_terms <- setof(c(new_terms, list(term)))
        }
      }
    }
    terms <- setof(c(new_terms, on_hold))
  }
  .add(unname(terms))
}
