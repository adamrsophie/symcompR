# Linearization of equations of motion, after
# sympy/physics/mechanics/linearize.py.

#' Linearizer of equations of motion
#'
#' `Linearizer()` holds the general form of a system's equations of motion
#' and linearizes them about an operating point, as SymPy's `Linearizer`.
#' It is usually made by the `to_linearizer()` method of [KanesMethod()]
#' or [LagrangesMethod()]; `linearize(op_point, A_and_B, simplify)`
#' (reached with `$`) returns `M`, `A` and `B`, or `A` and `B`.
#'
#' @param f_0,f_1,f_2,f_3,f_4,f_c,f_v,f_a The matrices of the general form.
#' @param q,u Generalized coordinates and speeds.
#' @param q_i,q_d,u_i,u_d Independent and dependent coordinates and speeds.
#' @param r Input variables.
#' @param lams Lagrange multipliers.
#' @param linear_solver `"LU"` or a function of `A` and `b`.
#' @return A `Linearizer` object.
#' @examples
#' q <- dynamicsymbols("q"); u <- dynamicsymbols("u")
#' qd <- dynamicsymbols("q", 1); ud <- dynamicsymbols("u", 1)
#' k <- Symbol("k")
#' lin <- Linearizer(list(qd), list(-u), list(ud), list(k * q), list(0),
#'                   list(), list(), list(), list(q), list(u),
#'                   q_i = list(q), u_i = list(u))
#' lin$linearize(A_and_B = TRUE)
#' @export
Linearizer <- function(f_0, f_1, f_2, f_3, f_4, f_c, f_v, f_a, q, u, q_i = NULL,
                       q_d = NULL, u_i = NULL, u_d = NULL, r = NULL, lams = NULL,
                       linear_solver = "LU") {
  self <- new.env(parent = emptyenv())
  class(self) <- c("Linearizer", "physics_obj")
  assign("_uid", .pv_next_uid(), envir = self)
  t <- .pv_t()
  self[["linear_solver"]] <- .pm_parse_linear_solver(linear_solver)
  self[["f_0"]] <- .pm_colmat(f_0)
  self[["f_1"]] <- .pm_colmat(f_1)
  self[["f_2"]] <- .pm_colmat(f_2)
  self[["f_3"]] <- .pm_colmat(f_3)
  self[["f_4"]] <- .pm_colmat(f_4)
  self[["f_c"]] <- .pm_colmat(f_c)
  self[["f_v"]] <- .pm_colmat(f_v)
  self[["f_a"]] <- .pm_colmat(f_a)
  self[["q"]] <- .pm_colmat(q)
  self[["u"]] <- .pm_colmat(u)
  self[["q_i"]] <- .pm_none_handler(q_i)
  self[["q_d"]] <- .pm_none_handler(q_d)
  self[["u_i"]] <- .pm_none_handler(u_i)
  self[["u_d"]] <- .pm_none_handler(u_d)
  self[["r"]] <- .pm_none_handler(r)
  self[["lams"]] <- .pm_none_handler(lams)
  self[["_qd"]] <- .pm_mdiff_t(self[["q"]])
  self[["_ud"]] <- .pm_mdiff_t(self[["u"]])
  uk <- vapply(self[["u"]]$data, .k, "")
  self[["_qd_dup"]] <- .pm_colmat(lapply(self[["_qd"]]$data, function(v) if (v$key %in% uk) Dummy() else v))
  self[["_dims"]] <- list(l = .pm_len(self[["f_c"]]), m = .pm_len(self[["f_v"]]), n = .pm_len(self[["q"]]),
                          o = .pm_len(self[["u"]]), s = .pm_len(self[["r"]]), k = .pm_len(self[["lams"]]))
  for (nm in c("_Pq", "_Pqi", "_Pqd", "_Pu", "_Pui", "_Pud", "_C_0", "_C_1", "_C_2", "perm_mat"))
    self[[nm]] <- NULL
  self[["_setup_done"]] <- FALSE
  self
}

.pm_lin_setup <- function(self) {
  .pm_lin_perm(self)
  .pm_lin_blocks(self)
  .pm_lin_coef(self)
  self[["_setup_done"]] <- TRUE
}

.pvM_Linearizer__setup <- function(self) invisible(.pm_lin_setup(self))

.pm_lin_perm <- function(self) {
  d <- self[["_dims"]]
  l <- d$l; m <- d$m; n <- d$n; o <- d$o; k <- d$k
  if (n != 0L) {
    self[["_Pq"]] <- .pm_permutation_matrix(self[["q"]], .pm_vstack(list(self[["q_i"]], self[["q_d"]])))
    if (l > 0L) {
      self[["_Pqi"]] <- .pm_cols(self[["_Pq"]], 0L, n - l)
      self[["_Pqd"]] <- .pm_cols(self[["_Pq"]], n - l, n)
    } else {
      self[["_Pqi"]] <- self[["_Pq"]]
      self[["_Pqd"]] <- .pm_empty()
    }
  }
  if (o != 0L) {
    self[["_Pu"]] <- .pm_permutation_matrix(self[["u"]], .pm_vstack(list(self[["u_i"]], self[["u_d"]])))
    if (m > 0L) {
      self[["_Pui"]] <- .pm_cols(self[["_Pu"]], 0L, o - m)
      self[["_Pud"]] <- .pm_cols(self[["_Pu"]], o - m, o)
    } else {
      self[["_Pui"]] <- self[["_Pu"]]
      self[["_Pud"]] <- .pm_empty()
    }
  }
  Pqi <- if (is.null(self[["_Pqi"]])) .pm_empty() else self[["_Pqi"]]
  Pui <- if (is.null(self[["_Pui"]])) .pm_empty() else self[["_Pui"]]
  P_col1 <- .pm_vstack(list(Pqi, .mech_zeros(o + k, n - l)))
  P_col2 <- .pm_vstack(list(.mech_zeros(n, o - m), Pui, .mech_zeros(k, o - m)))
  self[["perm_mat"]] <- if (.pm_truthy(P_col1)) {
    if (.pm_truthy(P_col2)) .pm_rjoin(P_col1, P_col2) else P_col1
  } else P_col2
}

.pm_eye_minus <- function(n, A) .pv_msub(.mech_eye(n), A)

.pm_lin_coef <- function(self) {
  d <- self[["_dims"]]
  l <- d$l; m <- d$m; n <- d$n; o <- d$o
  ls <- self[["linear_solver"]]
  if (l > 0L) {
    f_c_jac_q <- .pm_jacobian(self[["f_c"]], self[["q"]])
    self[["_C_0"]] <- .pv_matmul(.pm_eye_minus(n, .pv_matmul(self[["_Pqd"]], ls(.pv_matmul(f_c_jac_q, self[["_Pqd"]]), f_c_jac_q))),
                                 self[["_Pqi"]])
  } else self[["_C_0"]] <- .mech_eye(n)
  if (m > 0L) {
    f_v_jac_u <- .pm_jacobian(self[["f_v"]], self[["u"]])
    temp <- .pv_matmul(f_v_jac_u, self[["_Pud"]])
    if (n != 0L) {
      f_v_jac_q <- .pm_jacobian(self[["f_v"]], self[["q"]])
      self[["_C_1"]] <- .pv_matmul(.pm_mneg(self[["_Pud"]]), ls(temp, f_v_jac_q))
    } else self[["_C_1"]] <- .mech_zeros(o, n)
    self[["_C_2"]] <- .pv_matmul(.pm_eye_minus(o, .pv_matmul(self[["_Pud"]], ls(temp, f_v_jac_u))), self[["_Pui"]])
  } else {
    self[["_C_1"]] <- .mech_zeros(o, n)
    self[["_C_2"]] <- .mech_eye(o)
  }
}

.pm_lin_blocks <- function(self) {
  d <- self[["_dims"]]
  l <- d$l; m <- d$m; n <- d$n; o <- d$o; s <- d$s; k <- d$k
  J <- .pm_jacobian
  E <- .pm_empty()
  f0 <- self[["f_0"]]; f1 <- self[["f_1"]]; f2 <- self[["f_2"]]; f3 <- self[["f_3"]]
  f4 <- self[["f_4"]]; fa <- self[["f_a"]]; q <- self[["q"]]; u <- self[["u"]]
  if (n != 0L) {
    self[["_M_qq"]] <- J(f0, self[["_qd"]])
    self[["_A_qq"]] <- .pm_mneg(J(.pv_madd(f0, f1), q))
  } else self[["_M_qq"]] <- self[["_A_qq"]] <- E
  if (n != 0L && m != 0L) {
    self[["_M_uqc"]] <- J(fa, self[["_qd_dup"]])
    self[["_A_uqc"]] <- .pm_mneg(J(fa, q))
  } else self[["_M_uqc"]] <- self[["_A_uqc"]] <- E
  if (n != 0L && o - m + k != 0L) {
    self[["_M_uqd"]] <- J(f3, self[["_qd_dup"]])
    self[["_A_uqd"]] <- .pm_mneg(J(.pv_madd(.pv_madd(f2, f3), f4), q))
  } else self[["_M_uqd"]] <- self[["_A_uqd"]] <- E
  if (o != 0L && m != 0L) {
    self[["_M_uuc"]] <- J(fa, self[["_ud"]])
    self[["_A_uuc"]] <- .pm_mneg(J(fa, u))
  } else self[["_M_uuc"]] <- self[["_A_uuc"]] <- E
  if (o != 0L && o - m + k != 0L) {
    self[["_M_uud"]] <- J(f2, self[["_ud"]])
    self[["_A_uud"]] <- .pm_mneg(J(.pv_madd(f2, f3), u))
  } else self[["_M_uud"]] <- self[["_A_uud"]] <- E
  self[["_A_qu"]] <- if (o != 0L && n != 0L) .pm_mneg(J(f1, u)) else E
  self[["_M_uld"]] <- if (k != 0L && o - m + k != 0L) J(f4, self[["lams"]]) else E
  self[["_B_u"]] <- if (s != 0L && o - m + k != 0L) .pm_mneg(J(f3, self[["r"]])) else E
}

# op_point: a dict, a list of dicts, or a list of pairs
.pm_op_point_dict <- function(op_point) {
  if (is.null(op_point)) return(.pydict(list(), list()))
  if (inherits(op_point, "pydict")) return(op_point)
  if (is.list(op_point) && length(op_point) && all(vapply(op_point, function(z) inherits(z, "pydict"), TRUE))) {
    m <- .pm_merge_dicts(op_point)
  } else if (is.list(op_point)) {
    m <- .pm_merge_dicts(list(op_point))
  } else return(.pydict(list(), list()))
  .pydict(m$keys, m$vals)
}

.pvM_Linearizer_linearize <- function(self, op_point = NULL, A_and_B = FALSE, simplify = FALSE) {
  if (!isTRUE(self[["_setup_done"]])) .pm_lin_setup(self)
  opd <- .pm_op_point_dict(op_point)
  d <- self[["_dims"]]
  l <- d$l; m <- d$m; n <- d$n; o <- d$o; s <- d$s; k <- d$k
  g <- function(nm) self[[nm]]
  mm <- .pv_matmul
  ma <- .pv_madd
  if (o != 0L) col2 <- .pm_vstack(list(.mech_zeros(n, o), g("_M_uuc"), g("_M_uud")))
  if (k != 0L) col3 <- .pm_vstack(list(.mech_zeros(n + m, k), g("_M_uld")))
  if (n != 0L) {
    col1 <- .pm_vstack(list(g("_M_qq"), g("_M_uqc"), g("_M_uqd")))
    M <- if (o != 0L && k != 0L) .pm_rjoin(.pm_rjoin(col1, col2), col3)
         else if (o != 0L) .pm_rjoin(col1, col2) else col1
  } else if (k != 0L) M <- .pm_rjoin(col2, col3) else M <- col2
  M_eq <- msubs(M, opd)
  C_0 <- g("_C_0"); C_1 <- g("_C_1"); C_2 <- g("_C_2")
  if (n != 0L) {
    r1c1 <- g("_A_qq")
    if (o != 0L) r1c1 <- ma(r1c1, mm(g("_A_qu"), C_1))
    r1c1 <- mm(r1c1, C_0)
    if (m != 0L) {
      r2c1 <- g("_A_uqc")
      if (o != 0L) r2c1 <- ma(r2c1, mm(g("_A_uuc"), C_1))
      r2c1 <- mm(r2c1, C_0)
    } else r2c1 <- .pm_empty()
    if (o - m + k != 0L) {
      r3c1 <- g("_A_uqd")
      if (o != 0L) r3c1 <- ma(r3c1, mm(g("_A_uud"), C_1))
      r3c1 <- mm(r3c1, C_0)
    } else r3c1 <- .pm_empty()
    col1 <- .pm_vstack(list(r1c1, r2c1, r3c1))
  } else col1 <- .pm_empty()
  if (o != 0L) {
    r1c2 <- if (n != 0L) mm(g("_A_qu"), C_2) else .pm_empty()
    r2c2 <- if (m != 0L) mm(g("_A_uuc"), C_2) else .pm_empty()
    r3c2 <- if (o - m + k != 0L) mm(g("_A_uud"), C_2) else .pm_empty()
    col2 <- .pm_vstack(list(r1c2, r2c2, r3c2))
  } else col2 <- .pm_empty()
  Amat <- if (.pm_truthy(col1)) { if (.pm_truthy(col2)) .pm_rjoin(col1, col2) else col1 } else col2
  Amat_eq <- msubs(Amat, opd)
  if (s != 0L && o - m + k != 0L) {
    Bmat <- .pm_cjoin(.mech_zeros(n + m, s), g("_B_u"))
    Bmat_eq <- msubs(Bmat, opd)
  } else Bmat_eq <- .pm_empty()
  simp <- function(M) .pv_mmap(M, base::get("simplify", envir = .pv_ns))
  if (isTRUE(A_and_B)) {
    ls <- g("linear_solver")
    PT <- .pv_T(g("perm_mat"))
    A_cont <- mm(PT, ls(M_eq, Amat_eq))
    B_cont <- if (.pm_truthy(Bmat_eq)) mm(PT, ls(M_eq, Bmat_eq)) else Bmat_eq
    if (isTRUE(simplify)) {
      A_cont <- simp(A_cont)
      B_cont <- simp(B_cont)
    }
    return(.pytuple(list(A_cont, B_cont)))
  }
  if (isTRUE(simplify)) {
    M_eq <- simp(M_eq)
    Amat_eq <- simp(Amat_eq)
    Bmat_eq <- simp(Bmat_eq)
  }
  .pytuple(list(M_eq, Amat_eq, Bmat_eq))
}

.pm_flatten <- function(x) if (inherits(x, "symmatrix")) x$data else lapply(unclass(x), as_sym)

.pm_permutation_matrix <- function(orig_vec, per_vec) {
  orig_vec <- .pm_flatten(orig_vec)
  per_vec <- .pm_flatten(per_vec)
  ok <- vapply(orig_vec, .k, "")
  pk <- vapply(per_vec, .k, "")
  if (!setequal(ok, pk)) .pv_err("ValueError", "orig_vec and per_vec must be the same length, and contain the same symbols.")
  ind_list <- match(pk, ok)
  n <- length(orig_vec)
  p <- .mech_zeros(n, n)
  for (i in seq_along(ind_list)) p$data[[(i - 1L) * n + ind_list[i]]] <- S_One
  p
}
