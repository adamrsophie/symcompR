# Vectors and dyadics of classical mechanics, after
# sympy/physics/vector/vector.py and dyadic.py.

.pv_vec <- function(args) structure(list(args = args), class = c("PhysVector", "physics_obj"))
.pv_vec0 <- function() .pv_vec(list())
.pv_vargs <- function(v) .subset2(v, "args")

# Vector.__init__: a list of (matrix, frame) pairs is summed per frame in
# first appearance order; a dict (a pydict keyed by frames) is used as is.
.pv_vec_new <- function(inlist) {
  if (.pv_is0(inlist)) inlist <- list()
  d <- list()
  if (inherits(inlist, "pydict")) {
    for (i in seq_along(inlist$keys)) d <- .pvod_set(d, inlist$keys[[i]], inlist$values[[i]])
  } else {
    for (inp in inlist) {
      inp <- as.list(inp)
      fr <- inp[[2L]]
      if (.pvod_has(d, fr)) d <- .pvod_set(d, fr, .pv_madd(.pvod_get(d, fr), inp[[1L]]))
      else d <- .pvod_set(d, fr, Matrix(inp[[1L]]))
    }
  }
  args <- list()
  for (e in d) if (!.pv_mzero(e$val)) args[[length(args) + 1L]] <- list(e$val, e$key)
  .pv_vec(args)
}

#' Vectors and dyadics of classical mechanics
#'
#' `PhysVector()` builds a vector of `sympy.physics.vector` from a list of
#' `list(matrix, frame)` pairs, and `PhysVector(0)` is the zero vector;
#' vectors are usually made from the unit vectors of a [ReferenceFrame()],
#' as `x * N$x + y * N$y`. `PhysDyadic()` builds a dyadic from a list of
#' `list(coefficient, vector1, vector2)` triples; `N$x | N$y` is the outer
#' product. The names differ from SymPy's `Vector` and `Dyadic` because
#' `sympy.vector` has classes of the same names.
#'
#' Operators follow SymPy: `+`, `-`, scalar `*` and `/`, `==`, `&` (dot
#' product), `^` (cross product) and `|` (outer product). Methods such as
#' `dot`, `cross`, `outer`, `express`, `dt`, `diff`, `magnitude`,
#' `normalize`, `subs`, `simplify`, `to_matrix`, `doit`, `applyfunc`,
#' `xreplace`, `separate` and `free_symbols` are reached with `$`.
#'
#' @param inlist `0`, or a list of pairs (vectors) or triples (dyadics).
#' @return A `PhysVector` or `PhysDyadic`.
#' @examples
#' N <- ReferenceFrame("N")
#' x <- Symbol("x")
#' v <- x * N$x + 2 * N$y
#' v
#' v$magnitude()
#' (N$x | N$y) & N$y
#' PhysVector(0)
#' @export
PhysVector <- function(inlist) .pv_vec_new(inlist)

.pv_vadd <- function(self, other) {
  if (.pv_is0(other)) return(self)
  other <- .pv_check_vector(other)
  .pv_vec_new(c(.pv_vargs(self), .pv_vargs(other)))
}

.pv_vmul <- function(self, other) {
  other <- as_sym(other)
  .pv_vec_new(lapply(.pv_vargs(self), function(a) list(.pv_mscale(a[[1L]], other), a[[2L]])))
}

.pv_vneg <- function(self) .pv_vmul(self, S_NegativeOne)

.pv_vsub <- function(self, other) {
  o <- if (.pv_is_vec(other)) .pv_vmul(other, S_NegativeOne)
       else if (.pv_is_dyad(other)) .pv_dmul(other, S_NegativeOne)
       else .mul(list(as_sym(other), S_NegativeOne))
  .pv_vadd(self, o)
}

.pv_vdiv <- function(self, other) .pv_vmul(self, .div(S_One, as_sym(other)))

.pv_vdot <- function(self, other) {
  if (.pv_is_dyad(other)) {
    ol <- .pv_vec0()
    for (v in .pv_dargs(other))
      ol <- .pv_vadd(ol, .pv_vmul(.pv_vmul(v[[3L]], v[[1L]]), .pv_vdot(v[[2L]], self)))
    return(ol)
  }
  other <- .pv_check_vector(other)
  out <- S_Zero
  for (v1 in .pv_vargs(self)) for (v2 in .pv_vargs(other)) {
    m <- .pv_matmul(.pv_matmul(.pv_T(v2[[1L]]), .pv_dcm(v2[[2L]], v1[[2L]])), v1[[1L]])
    out <- .add(list(out, m$data[[1L]]))
  }
  if (.pv_state$simp) trigsimp(out) else out
}

.pv_vec_eq <- function(self, other) {
  if (.pv_is0(other)) other <- .pv_vec0()
  if (!.pv_is_vec(other)) return(FALSE)
  a <- .pv_vargs(self); b <- .pv_vargs(other)
  if (!length(a) && !length(b)) return(TRUE)
  if (!length(a) || !length(b)) return(FALSE)
  frame <- a[[1L]][[2L]]
  diffv <- .pv_vsub(self, other)
  for (v in .pv_frame_iter(frame)) if (!is_zero_num(expand(.pv_vdot(diffv, v)))) return(FALSE)
  TRUE
}

.pv_outer <- function(self, other) {
  other <- .pv_check_vector(other)
  ol <- .pv_dyad0()
  for (v in .pv_vargs(self)) for (v2 in .pv_vargs(other)) {
    fv <- .pv_frame_iter(v[[2L]]); fv2 <- .pv_frame_iter(v2[[2L]])
    for (i in 1:3) for (j in 1:3)
      ol <- .pv_dadd(ol, .pv_dyad_new(list(list(.mul(list(.pv_e(v[[1L]], i), .pv_e(v2[[1L]], j))), fv[[i]], fv2[[j]]))))
  }
  ol
}

.pv_det3 <- function(m) {
  # m[0] holds vectors, rows 1 and 2 scalars
  t1 <- .pv_vmul(m[[1L]][[1L]], .sub(.mul(list(m[[2L]][[2L]], m[[3L]][[3L]])), .mul(list(m[[2L]][[3L]], m[[3L]][[2L]]))))
  t2 <- .pv_vmul(m[[1L]][[2L]], .sub(.mul(list(m[[2L]][[3L]], m[[3L]][[1L]])), .mul(list(m[[2L]][[1L]], m[[3L]][[3L]]))))
  t3 <- .pv_vmul(m[[1L]][[3L]], .sub(.mul(list(m[[2L]][[1L]], m[[3L]][[2L]])), .mul(list(m[[2L]][[2L]], m[[3L]][[1L]]))))
  .pv_vadd(.pv_vadd(t1, t2), t3)
}

.pv_vcross <- function(self, other) {
  if (.pv_is_dyad(other)) {
    ol <- .pv_dyad0()
    for (v in .pv_dargs(other))
      ol <- .pv_dadd(ol, .pv_dmul(.pv_outer(.pv_vcross(self, v[[2L]]), v[[3L]]), v[[1L]]))
    return(ol)
  }
  other <- .pv_check_vector(other)
  if (!length(.pv_vargs(other))) return(.pv_vec0())
  outlist <- list()
  for (v in .pv_vargs(other)) {
    fx <- .pv_frame_iter(v[[2L]])
    vv <- .pv_vec_new(list(v))
    m <- list(fx, lapply(fx, function(b) .pv_vdot(self, b)), lapply(fx, function(b) .pv_vdot(vv, b)))
    outlist <- c(outlist, .pv_vargs(.pv_det3(m)))
  }
  .pv_vec_new(outlist)
}

.pv_separate <- function(self) {
  keys <- list(); vals <- list()
  for (x in .pv_vargs(self)) {
    keys[[length(keys) + 1L]] <- x[[2L]]
    vals[[length(vals) + 1L]] <- .pv_vec_new(list(x))
  }
  .pydict(keys, vals)
}

.pv_vdiff <- function(self, var, frame, var_in_dcm = TRUE) {
  .pv_check_frame(frame)
  var <- as_sym(var)
  inlist <- list()
  for (vc in .pv_vargs(self)) {
    mn <- vc[[1L]]
    cf <- vc[[2L]]
    if (identical(cf, frame)) {
      inlist[[length(inlist) + 1L]] <- list(.pv_mdiff(mn, var), frame)
    } else if (!var_in_dcm || .pv_mzero(.pv_mdiff(.pv_dcm(frame, cf), var))) {
      inlist[[length(inlist) + 1L]] <- list(.pv_mdiff(mn, var), cf)
    } else {
      re <- .pv_express(.pv_vec_new(list(vc)), frame)
      deriv <- .pv_mdiff(.pv_vargs(re)[[1L]][[1L]], var)
      inlist <- c(inlist, .pv_vargs(.pv_vec_new(list(list(deriv, frame)))))
    }
  }
  .pv_vec_new(inlist)
}

.pv_to_matrix <- function(self, reference_frame) {
  .pv_colv(lapply(.pv_frame_iter(reference_frame), function(u) .pv_vdot(self, u)))
}

.pv_vdict <- function(self, f) {
  keys <- list(); vals <- list()
  for (v in .pv_vargs(self)) {
    keys[[length(keys) + 1L]] <- v[[2L]]
    vals[[length(vals) + 1L]] <- f(v[[1L]])
  }
  .pv_vec_new(.pydict(keys, vals))
}

.pv_magnitude <- function(self) .sqrt(.pv_vdot(self, self))
.pv_normalize <- function(self) .pv_vdiv(.pv_vec(.pv_vargs(self)), .pv_magnitude(self))

# expr.subs(*args) as SymPy reads its arguments: a dict (an R list of pairs
# or a named list) is applied with the most complex keys first.
.pv_subs <- function(e, ...) {
  a <- list(...)
  if (length(a) == 1L && (inherits(a[[1L]], "pydict") || (is.list(a[[1L]]) && !is_sym(a[[1L]])))) {
    d <- a[[1L]]
    pairs <- if (inherits(d, "pydict")) Map(list, d$keys, d$values) else .subs_pairs(d)
    pairs <- lapply(pairs, function(p) {
      if (is.character(p[[1L]])) list(.resolve_name(as_sym(e), p[[1L]]), as_sym(p[[2L]]))
      else list(as_sym(p[[1L]]), as_sym(p[[2L]]))
    })
    if (length(pairs) > 1L) {
      keys <- lapply(pairs, function(p) list(-.pv_nodes(p[[1L]]), .sort_key(p[[1L]])))
      pairs <- .key_sort(keys, pairs)
    }
    return(subs(e, pairs))
  }
  do.call(subs, c(list(e), a))
}

# sympy.core.sorting._nodes
.pv_nodes <- function(e) {
  if (e$type == "Derivative") {
    n <- .pv_nodes(e$args[[1L]])
    for (vc in e$args[-1L]) {
      c <- vc$args[[2L]]
      n <- n + if (is_number_type(c)) as.numeric(c$value) else .pv_nodes(c)
    }
    return(n)
  }
  .pv_node_count(e)
}
.pv_node_count <- function(e) if (e$type == "Float") 0.5 else 1 + sum(vapply(e$args, .pv_node_count, 0))

.pv_msubs <- function(M, ...) .pv_mmap(M, function(z) .pv_subs(z, ...))

.pv_angle_between <- function(self, vec) {
  v1 <- .pv_normalize(self)
  v2 <- .pv_normalize(vec)
  fcall("acos", list(.pv_vdot(v1, v2)))
}

.pv_free_syms_set <- function(e) structure(free_symbols(e), pyset = TRUE)

# --- Vector methods reached with `$`

.pvM_PhysVector_dot <- function(self, other) .pv_vdot(self, other)
.pvM_PhysVector_cross <- function(self, other) .pv_vcross(self, other)
.pvM_PhysVector_outer <- function(self, other) .pv_outer(self, other)
.pvM_PhysVector_separate <- function(self) .pv_separate(self)
.pvM_PhysVector_diff <- function(self, var, frame, var_in_dcm = TRUE) .pv_vdiff(self, var, frame, var_in_dcm)
.pvM_PhysVector_express <- function(self, otherframe, variables = FALSE) .pv_express(self, otherframe, variables = variables)
.pvM_PhysVector_to_matrix <- function(self, reference_frame) .pv_to_matrix(self, reference_frame)
.pvM_PhysVector_doit <- function(self, ...) .pv_vdict(self, function(m) .pv_mmap(m, function(z) doit(z)))
.pvM_PhysVector_dt <- function(self, otherframe) .pv_time_derivative(self, otherframe)
.pvM_PhysVector_simplify <- function(self) .pv_vdict(self, function(m) .pv_mmap(m, simplify))
.pvM_PhysVector_subs <- function(self, ...) {
  a <- list(...)
  .pv_vdict(self, function(m) .pv_mmap(m, function(z) do.call(.pv_subs, c(list(z), a))))
}
.pvM_PhysVector_magnitude <- function(self) .pv_magnitude(self)
.pvM_PhysVector_normalize <- function(self) .pv_normalize(self)
.pvM_PhysVector_applyfunc <- function(self, f) {
  if (!is.function(f)) .pv_err("TypeError", "`f` must be callable.")
  .pv_vdict(self, function(m) .pv_mmap(m, function(z) as_sym(f(z))))
}
.pvM_PhysVector_angle_between <- function(self, vec) .pv_angle_between(self, vec)
.pvM_PhysVector_free_symbols <- function(self, reference_frame) .pv_free_syms_set(.pv_to_matrix(self, reference_frame))
.pvM_PhysVector_free_dynamicsymbols <- function(self, reference_frame) find_dynamicsymbols(self, reference_frame = reference_frame)
.pv_evalf_subs <- function(subs) {
  if (is.null(subs)) return(list())
  if (inherits(subs, "pydict")) return(stats::setNames(subs$values, vapply(subs$keys, function(k) as_sym(k)$name, "")))
  subs
}
.pvM_PhysVector_evalf <- function(self, n = 15, subs = NULL, ...) {
  if (!length(.pv_vargs(self))) return(self)
  s <- .pv_evalf_subs(subs)
  .pv_vec_new(lapply(.pv_vargs(self), function(a) list(.pv_mmap(a[[1L]], function(z) N(z, n, subs = s)), a[[2L]])))
}
.pvM_PhysVector_n <- .pvM_PhysVector_evalf
.pvM_PhysVector_xreplace <- function(self, rule) {
  if (missing(rule)) .pv_err("TypeError", "Vector.xreplace() missing 1 required positional argument: 'rule'")
  .pv_vec_new(lapply(.pv_vargs(self), function(a) list(.pv_mmap(a[[1L]], function(z) .pv_xreplace(z, rule)), a[[2L]])))
}
.pvP_PhysVector_func <- function(self) PhysVector

.pv_xreplace <- function(e, rule) {
  if (inherits(rule, "pydict")) rule <- Map(list, rule$keys, rule$values)
  else if (is.list(rule) && length(rule) && all(vapply(rule, is_sym, TRUE))) {
    # Basic.xreplace with a list: `node in rule` holds, then rule[node] fails
    keys <- vapply(rule, .k, "")
    hit <- Filter(function(z) z$key %in% keys, c(list(e), .all_subexprs(e)))
    if (length(hit)) .pv_err("TypeError", sprintf("list indices must be integers or slices, not %s",
                                                  sub("^<class '(.*\\.)?([^.']*)'>$", "\\2", .pv_pyclass(hit[[1L]]))))
    return(e)
  }
  xreplace(e, rule)
}

# --- printing

# Vector._sympystr: frames sorted by creation index when there are several
.pv_coef_str <- function(c, prefix_sep) {
  s <- .str(c)
  if (c$type == "Add") s <- paste0("(", s, ")")
  if (startsWith(s, "-")) list(" - ", substring(s, 2L)) else list(" + ", s)
}

.pv_vec_str <- function(v, order = TRUE) {
  ar <- .pv_vargs(v)
  if (!length(ar)) return(.str(S_Zero))
  if (order && length(ar) > 1L) {
    idx <- vapply(ar, function(a) get("index", envir = a[[2L]]), 0L)
    ar <- ar[order(idx, method = "radix")]
  }
  ol <- character()
  for (a in ar) {
    sv <- get("str_vecs", envir = a[[2L]])
    for (j in 1:3) {
      c <- .pv_e(a[[1L]], j)
      if (identical(c$key, S_One$key)) ol <- c(ol, paste0(" + ", sv[j]))
      else if (identical(c$key, S_NegativeOne$key)) ol <- c(ol, paste0(" - ", sv[j]))
      else if (!identical(c$key, S_Zero$key)) {
        p <- .pv_coef_str(c)
        ol <- c(ol, paste0(p[[1L]], p[[2L]], "*", sv[j]))
      }
    }
  }
  out <- paste(ol, collapse = "")
  if (startsWith(out, " + ")) substring(out, 4L) else if (startsWith(out, " ")) substring(out, 2L) else out
}

#' @export
format.PhysVector <- function(x, ...) .pv_vec_str(x)

# --- operators

.pv_ops <- function(gen, e1, e2) {
  if (missing(e2)) {
    if (gen == "-") return(if (.pv_is_vec(e1)) .pv_vneg(e1) else .pv_dmul(e1, S_NegativeOne))
    if (gen == "+") return(e1)
    .pv_err("TypeError", sprintf("bad operand type for unary %s", gen))
  }
  v1 <- .pv_is_vec(e1); v2 <- .pv_is_vec(e2)
  d1 <- .pv_is_dyad(e1); d2 <- .pv_is_dyad(e2)
  switch(gen,
    "+" = {
      if (v1) return(.pv_vadd(e1, e2))
      if (d1) return(.pv_dadd(e1, e2))
      if (v2) return(.pv_vadd(e2, e1))
      return(.pv_dadd(e2, e1))
    },
    "-" = {
      if (v1) return(.pv_vsub(e1, e2))
      if (d1) return(.pv_dadd(e1, if (d2) .pv_dmul(e2, S_NegativeOne) else .pv_scalar_neg(e2)))
      # __rsub__: (-1 * self) + other
      if (v2) return(.pv_vadd(.pv_vneg(e2), e1))
      return(.pv_dadd(.pv_dmul(e2, S_NegativeOne), e1))
    },
    "*" = {
      if (v1 || d1) {
        if (v2 || d2) .pv_err("TypeError", "unsupported operand type(s) for *")
        return(if (v1) .pv_vmul(e1, e2) else .pv_dmul(e1, e2))
      }
      return(if (v2) .pv_vmul(e2, e1) else .pv_dmul(e2, e1))
    },
    "/" = {
      if (v1) return(.pv_vdiv(e1, e2))
      if (d1) return(.pv_dmul(e1, .div(S_One, as_sym(e2))))
      .pv_err("TypeError", "unsupported operand type(s) for /")
    },
    "==" = return(.pv_eq_any(e1, e2)),
    "!=" = return(!.pv_eq_any(e1, e2)),
    "&" = {
      if (v1) return(.pv_vdot(e1, e2))
      if (d1) return(.pv_ddot(e1, e2))
      # Vector.__rand__
      if (v2) return(.pv_vdot(e2, e1))
      .pv_unsupported(gen, e1, e2)
    },
    "^" = {
      if (v1) return(.pv_vcross(e1, e2))
      if (d1) return(.pv_dcross(e1, e2))
      .pv_unsupported(gen, e1, e2)
    },
    "|" = {
      if (v1) return(.pv_outer(e1, e2))
      .pv_unsupported(gen, e1, e2)
    },
    .pv_unsupported(gen, e1, e2))
}

.pv_unsupported <- function(gen, e1, e2)
  .pv_err("TypeError", sprintf("unsupported operand type(s) for %s: '%s' and '%s'", gen, .pv_tname(e1), .pv_tname(e2)))

.pv_scalar_neg <- function(x) if (is.numeric(x)) -x else .neg(as_sym(x))

.pv_eq_any <- function(e1, e2) {
  if (.pv_is_vec(e1)) return(.pv_vec_eq(e1, e2))
  if (.pv_is_dyad(e1)) return(.pv_dyad_eq(e1, e2))
  if (.pv_is_vec(e2)) return(.pv_vec_eq(e2, e1))
  .pv_dyad_eq(e2, e1)
}

#' @export
Ops.PhysVector <- function(e1, e2) .pv_ops(.Generic, e1, e2)

#' @export
Ops.PhysDyadic <- function(e1, e2) .pv_ops(.Generic, e1, e2)

#' @export
chooseOpsMethod.physics_obj <- function(x, y, mx, my, cl, reverse) TRUE

#' @export
diff.PhysVector <- function(x, var, frame, var_in_dcm = TRUE, ...) .pv_vdiff(x, var, frame, var_in_dcm)

# ---------------------------------------------------------------------------
# Dyadic

.pv_dyad <- function(args) structure(list(args = args), class = c("PhysDyadic", "physics_obj"))
.pv_dyad0 <- function() .pv_dyad(list())
.pv_dargs <- function(d) .subset2(d, "args")

# Dyadic.__init__: terms with the same printed unit vectors are merged, and
# terms with a zero coefficient or zero vector are dropped.
.pv_dyad_new <- function(inlist) {
  if (.pv_is0(inlist)) inlist <- list()
  args <- list()
  strs <- list()
  for (t in inlist) {
    t <- as.list(t)
    t[[1L]] <- as_sym(t[[1L]])
    s1 <- .pv_str(t[[2L]]); s2 <- .pv_str(t[[3L]])
    hit <- 0L
    for (i in seq_along(args)) {
      if (strs[[i]][1L] == s1 && strs[[i]][2L] == s2) { hit <- i; break }
    }
    if (hit) args[[hit]][[1L]] <- .add(list(args[[hit]][[1L]], t[[1L]]))
    else {
      args[[length(args) + 1L]] <- t
      strs[[length(strs) + 1L]] <- c(s1, s2)
    }
  }
  keep <- vapply(args, function(a) !(.pv_is0(a[[1L]]) || .pv_is0(a[[2L]]) || .pv_is0(a[[3L]])), TRUE)
  .pv_dyad(args[keep])
}

#' @rdname PhysVector
#' @export
PhysDyadic <- function(inlist) .pv_dyad_new(inlist)

.pv_dadd <- function(self, other) {
  other <- .pv_check_dyadic(other)
  .pv_dyad_new(c(.pv_dargs(self), .pv_dargs(other)))
}

.pv_dmul <- function(self, other) {
  other <- as_sym(other)
  .pv_dyad_new(lapply(.pv_dargs(self), function(a) list(.mul(list(other, a[[1L]])), a[[2L]], a[[3L]])))
}

.pv_ddot <- function(self, other) {
  if (.pv_is_dyad(other)) {
    ol <- .pv_dyad0()
    for (v in .pv_dargs(self)) for (v2 in .pv_dargs(other))
      ol <- .pv_dadd(ol, .pv_dmul(.pv_outer(v[[2L]], v2[[3L]]), .mul(list(v[[1L]], v2[[1L]], .pv_vdot(v[[3L]], v2[[2L]])))))
    return(ol)
  }
  other <- .pv_check_vector(other)
  ol <- .pv_vec0()
  for (v in .pv_dargs(self))
    ol <- .pv_vadd(ol, .pv_vmul(.pv_vmul(v[[2L]], v[[1L]]), .pv_vdot(v[[3L]], other)))
  ol
}

.pv_dcross <- function(self, other) {
  other <- .pv_check_vector(other)
  ol <- .pv_dyad0()
  for (v in .pv_dargs(self))
    ol <- .pv_dadd(ol, .pv_dmul(.pv_outer(v[[2L]], .pv_vcross(v[[3L]], other)), v[[1L]]))
  ol
}

.pv_dyad_eq <- function(self, other) {
  if (.pv_is0(other)) other <- .pv_dyad0()
  other <- .pv_check_dyadic(other)
  a <- .pv_dargs(self); b <- .pv_dargs(other)
  if (!length(a) && !length(b)) return(TRUE)
  if (!length(a) || !length(b)) return(FALSE)
  # set(self.args) == set(other.args)
  member <- function(t, ts) any(vapply(ts, function(u) identical(t[[1L]]$key, u[[1L]]$key) &&
                                          .pv_vec_eq(t[[2L]], u[[2L]]) && .pv_vec_eq(t[[3L]], u[[3L]]), TRUE))
  all(vapply(a, member, TRUE, ts = b)) && all(vapply(b, member, TRUE, ts = a))
}

.pv_dyad_str <- function(d) {
  ar <- .pv_dargs(d)
  if (!length(ar)) return(.str(S_Zero))
  ol <- character()
  for (v in ar) {
    pv <- paste0("(", format(v[[2L]]), "|", format(v[[3L]]), ")")
    if (identical(v[[1L]]$key, S_One$key)) ol <- c(ol, paste0(" + ", pv))
    else if (identical(v[[1L]]$key, S_NegativeOne$key)) ol <- c(ol, paste0(" - ", pv))
    else if (!identical(v[[1L]]$key, S_Zero$key)) {
      p <- .pv_coef_str(v[[1L]])
      ol <- c(ol, paste0(p[[1L]], p[[2L]], "*", pv))
    }
  }
  out <- paste(ol, collapse = "")
  if (startsWith(out, " + ")) substring(out, 4L) else if (startsWith(out, " ")) substring(out, 2L) else out
}

#' @export
format.PhysDyadic <- function(x, ...) .pv_dyad_str(x)

.pv_dsum <- function(ds) Reduce(.pv_dadd, ds, .pv_dyad0())

.pvM_PhysDyadic_dot <- function(self, other) .pv_ddot(self, other)
.pvM_PhysDyadic_cross <- function(self, other) .pv_dcross(self, other)
.pvM_PhysDyadic_express <- function(self, frame1, frame2 = NULL) .pv_express(self, frame1, frame2)
.pvM_PhysDyadic_to_matrix <- function(self, reference_frame, second_reference_frame = NULL) {
  if (is.null(second_reference_frame)) second_reference_frame <- reference_frame
  vals <- list()
  for (i in .pv_frame_iter(reference_frame)) for (j in .pv_frame_iter(second_reference_frame))
    vals[[length(vals) + 1L]] <- .pv_vdot(.pv_vdot(i, self), j)
  .new_mat(vals, 3L, 3L)
}
.pvM_PhysDyadic_doit <- function(self, ...) .pv_dsum(lapply(.pv_dargs(self), function(v) .pv_dyad_new(list(list(doit(v[[1L]]), v[[2L]], v[[3L]])))))
.pvM_PhysDyadic_dt <- function(self, frame) .pv_time_derivative(self, frame)
.pvM_PhysDyadic_simplify <- function(self) {
  out <- .pv_dyad0()
  for (v in .pv_dargs(self)) out <- .pv_dadd(out, .pv_dyad_new(list(list(simplify(v[[1L]]), v[[2L]], v[[3L]]))))
  out
}
.pvM_PhysDyadic_subs <- function(self, ...) {
  a <- list(...)
  .pv_dsum(lapply(.pv_dargs(self), function(v) .pv_dyad_new(list(list(do.call(.pv_subs, c(list(v[[1L]]), a)), v[[2L]], v[[3L]])))))
}
.pvM_PhysDyadic_applyfunc <- function(self, f) {
  if (!is.function(f)) .pv_err("TypeError", "`f` must be callable.")
  out <- .pv_dyad0()
  for (v in .pv_dargs(self)) out <- .pv_dadd(out, .pv_dmul(.pv_outer(v[[2L]], v[[3L]]), as_sym(f(v[[1L]]))))
  out
}
.pvM_PhysDyadic_evalf <- function(self, n = 15, subs = NULL, ...) {
  if (!length(.pv_dargs(self))) return(self)
  s <- .pv_evalf_subs(subs)
  .pv_dyad_new(lapply(.pv_dargs(self), function(v) list(N(v[[1L]], n, subs = s), v[[2L]], v[[3L]])))
}
.pvM_PhysDyadic_n <- .pvM_PhysDyadic_evalf
.pvM_PhysDyadic_xreplace <- function(self, rule) {
  if (missing(rule)) .pv_err("TypeError", "Dyadic.xreplace() missing 1 required positional argument: 'rule'")
  .pv_dyad_new(lapply(.pv_dargs(self), function(v) list(.pv_xreplace(v[[1L]], rule), v[[2L]], v[[3L]])))
}
.pvP_PhysDyadic_func <- function(self) PhysDyadic
