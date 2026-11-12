# Differential geometry, after sympy/diffgeom (diffgeom.py and rn.py).
# A coordinate system is a keyed leaf holding its patch, symbols and
# relations, so the fields built on it have no free symbols. Fields are
# applied with rcall(), SymPy's Basic.rcall.

#' Differential geometry
#'
#' `Manifold(name, dim)`, `Patch(name, manifold)` and
#' `CoordSystem(name, patch, symbols, relations)` build coordinate systems;
#' `relations` is a list of `list(c(from, to), syms, exprs)` entries.
#' `transformation()`, `cs_transform()` (SymPy's `CoordSystem.transform`), `jacobian_dg()` (SymPy's
#' `CoordSystem.jacobian`), `jacobian_determinant()`, `point()`,
#' `point_to_coords()`, `base_scalar()`, `base_scalars()`, `base_vector()`,
#' `base_vectors()`, `base_oneform()` and `base_oneforms()` are SymPy's
#' methods. `Differential()`, `TensorProduct()`, `WedgeProduct()`,
#' `Commutator()`, `LieDerivative()`, `BaseCovarDerivativeOp()` and
#' `CovarDerivativeOp()` build fields and operators, applied with
#' `rcall(field, ...)`. `contravariant_order()`, `covariant_order()`,
#' `vectors_in_basis()`, `twoform_to_matrix()`, `metric_to_Christoffel_1st()`,
#' `metric_to_Christoffel_2nd()`, `metric_to_Riemann_components()`,
#' `metric_to_Ricci_components()`, `intcurve_series()` and
#' `intcurve_diffequ()` are SymPy's functions. `R2`, `R2_r`, `R2_p`, `R3`,
#' `R3_r`, `R3_c` and `R3_s` (from `sympy.diffgeom.rn`) are built with
#' `diffgeom_rn()`.
#'
#' @param name A name.
#' @param dim The dimension.
#' @param manifold,patch,cs,sys,coord_sys A manifold, patch or coordinate system.
#' @param symbols Coordinate symbols.
#' @param relations Transformations between systems.
#' @param coordinates,coords Coordinates.
#' @param coord_index An index from 0.
#' @param pt A point.
#' @param form_field,expr,v_field,v1,v2,wrt,vector_field Fields.
#' @param ... Fields for `rcall()`, `TensorProduct()` and `WedgeProduct()`.
#' @param index An index.
#' @param christoffel Christoffel symbols, an array.
#' @param param The curve parameter.
#' @param start_point A point.
#' @param n The number of terms.
#' @param coeffs Return the coefficients.
#' @param to_sys A coordinate system.
#' @param .strict Return -1 for expressions unrelated to the module.
#' @return A sym object, matrix, array or list.
#' @examples
#' rn <- diffgeom_rn()
#' x <- base_scalar(rn$R2_r, 0); y <- base_scalar(rn$R2_r, 1)
#' e_x <- base_vector(rn$R2_r, 0)
#' rcall(e_x, x^2 * y)
#' @name diffgeom
NULL

#' @rdname diffgeom
#' @export
Manifold <- function(name, dim) new_sym("Manifold", list(.Str(name), as_sym(dim)))
register_str_printer("Manifold", function(e) e$args[[1L]]$name)
#' @rdname diffgeom
#' @export
Patch <- function(name, manifold) new_sym("Patch", list(.Str(name), manifold))
register_str_printer("Patch", function(e) e$args[[1L]]$name)

.dg_dim <- function(x) {
  switch(x$type,
    Manifold = .int_val(x$args[[2L]]),
    Patch = .dg_dim(x$args[[2L]]),
    CoordSystem = .dg_dim(x$extra$patch),
    BaseScalarField = , BaseVectorField = .dg_dim(x$args[[1L]]),
    Point = .dg_dim(x$args[[1L]]))
}

#' @rdname diffgeom
#' @export
CoordSystem <- function(name, patch, symbols = NULL, relations = list()) {
  dim <- .dg_dim(patch)
  if (is.null(symbols)) symbols <- lapply(seq_len(dim) - 1L, function(i) Symbol(paste0(name, "_", i), real = TRUE))
  symbols <- lapply(symbols, as_sym)
  rels <- lapply(relations, function(r) {
    if (is_sym(r[[2L]]) && r[[2L]]$type == "Lambda") list(from = r[[1L]][1L], to = r[[1L]][2L], syms = r[[2L]]$args[[1L]]$args, exprs = r[[2L]]$args[[2L]]$args)
    else list(from = r[[1L]][1L], to = r[[1L]][2L], syms = lapply(r[[2L]], as_sym), exprs = lapply(r[[3L]], as_sym))
  })
  rk <- paste(vapply(rels, function(r) paste(r$from, r$to, paste(vapply(r$exprs, `[[`, "", "key"), collapse = ","), sep = ">"), ""), collapse = ";")
  key <- paste0("CS:", name, "|", patch$key, "|", paste(vapply(symbols, `[[`, "", "key"), collapse = ","), "|", rk)
  new_sym("CoordSystem", list(), name = name, extra = list(patch = patch, symbols = symbols, relations = rels,
                                                         dummy = Dummy("cs")), key = key)
}
register_str_printer("CoordSystem", function(e) e$name)

#' @rdname diffgeom
#' @export
coord_symbols <- function(cs) cs$extra$symbols

.dg_rel <- function(cs, from, to) {
  for (r in cs$extra$relations) if (r$from == from && r$to == to) return(r)
  NULL
}

.dg_solve_inverse <- function(sym1, sym2, exprs, n1, n2) {
  eqs <- Map(function(a, b) .sub(a, b), sym2, exprs)
  ret <- solve(eqs, sym1, dict = TRUE)
  if (!length(ret)) stop("Cannot solve inverse relation from ", n1, " to ", n2, ".", call. = FALSE)
  if (length(ret) > 1L) stop("Obtained multiple inverse relation from ", n1, " to ", n2, ".", call. = FALSE)
  ret[[1L]]
}

.dg_dijkstra <- function(cs, to) {
  rels <- cs$extra$relations
  graph <- list()
  for (r in rels) {
    graph[[r$from]] <- unique(c(graph[[r$from]], r$to))
    graph[[r$to]] <- unique(c(graph[[r$to]], r$from))
  }
  pd <- lapply(graph, function(g) list(d = 0, path = character(), v = 0))
  visit <- function(s) {
    pd[[s]]$v <<- 1
    for (ns in graph[[s]]) {
      dist <- pd[[s]]$d + 1
      if (pd[[ns]]$d >= dist || !length(pd[[ns]]$path)) {
        pd[[ns]]$d <<- dist
        pd[[ns]]$path <<- c(pd[[s]]$path, s)
      }
    }
  }
  visit(cs$name)
  repeat {
    md <- max(vapply(pd, `[[`, 0, "d"))
    ns <- NULL
    for (s in names(pd)) if (0 < pd[[s]]$d && pd[[s]]$d <= md && !pd[[s]]$v) { md <- pd[[s]]$d; ns <- s }
    if (is.null(ns)) break
    visit(ns)
  }
  res <- c(pd[[to]]$path, to)
  if (identical(res, to)) stop("Two coordinate systems are not connected.", call. = FALSE)
  res
}

.dg_transf_exprs <- function(cs, sys) {
  if (identical(cs$key, sys$key)) return(cs$extra$symbols)
  r <- .dg_rel(cs, cs$name, sys$name)
  if (!is.null(r)) return(r$exprs)
  if (!is.null(.dg_rel(cs, sys$name, cs$name))) {
    fwd <- cs_transform(sys, cs)$data
    inv <- .dg_solve_inverse(sys$extra$symbols, cs$extra$symbols, fwd, sys$name, cs$name)
    return(lapply(sys$extra$symbols, function(s) inv[[s$name]]))
  }
  path <- .dg_dijkstra(cs, sys$name)
  transforms <- list()
  for (k in seq_len(length(path) - 1L)) {
    s1 <- path[k]; s2 <- path[k + 1L]
    r <- .dg_rel(cs, s1, s2)
    if (!is.null(r)) transforms[[length(transforms) + 1L]] <- list(r$syms, r$exprs)
    else {
      r <- .dg_rel(cs, s2, s1)
      d1 <- lapply(r$syms, function(z) Dummy())
      ret <- .dg_solve_inverse(r$syms, d1, r$exprs, s2, s1)
      transforms[[length(transforms) + 1L]] <- list(d1, lapply(r$syms, function(s) ret[[s$name]]))
    }
  }
  exprs <- cs$extra$symbols
  for (t in transforms) exprs <- lapply(t[[2L]], function(e) subs(e, Map(list, t[[1L]], exprs)))
  exprs
}

#' @rdname diffgeom
#' @export
transformation <- function(cs, sys) Lambda(cs$extra$symbols, .mat_imm(Matrix(.dg_transf_exprs(cs, sys))))
.mat_imm <- function(M) new_sym("ImmutableMatrix", list(.int(M$nrow), .int(M$ncol), do.call(Tuple, M$data)))

#' @rdname diffgeom
#' @export
cs_transform <- function(cs, sys, coordinates = NULL) {
  if (is.null(coordinates)) coordinates <- cs$extra$symbols
  coordinates <- lapply(coordinates, as_sym)
  if (identical(cs$key, sys$key)) return(Matrix(coordinates))
  ex <- .dg_transf_exprs(cs, sys)
  Matrix(lapply(ex, function(e) subs(e, Map(list, cs$extra$symbols, coordinates), simultaneous = TRUE)))
}

#' @rdname diffgeom
#' @export
jacobian_dg <- function(cs, sys, coordinates = NULL) {
  res <- jacobian(cs_transform(cs, sys)$data, cs$extra$symbols)
  if (!is.null(coordinates)) res <- .mat_map(res, function(z) subs(z, Map(list, cs$extra$symbols, lapply(coordinates, as_sym))))
  res
}
#' @rdname diffgeom
#' @export
jacobian_determinant <- function(cs, sys, coordinates = NULL) det(jacobian_dg(cs, sys, coordinates))

#' @rdname diffgeom
#' @export
point <- function(cs, coords) new_sym("Point", list(cs, do.call(Tuple, lapply(coords, as_sym))))
register_str_printer("Point", function(e) paste0("Point(", e$args[[1L]]$name, ", Matrix([", paste(paste0("[", vapply(e$args[[2L]]$args, .str, ""), "]"), collapse = ", "), "]))"))
#' @rdname diffgeom
#' @export
point_to_coords <- function(cs, pt) .dg_point_coords(pt, cs)
.dg_point_coords <- function(pt, sys = NULL) {
  cs <- pt$args[[1L]]; co <- pt$args[[2L]]$args
  if (is.null(sys)) Matrix(co) else cs_transform(cs, sys, co)
}

#' @rdname diffgeom
#' @export
base_scalar <- function(cs, coord_index) new_sym("BaseScalarField", list(cs, .int(coord_index)))
#' @rdname diffgeom
#' @export
base_scalars <- function(cs) if (.is_coordsys(cs)) .pytuple(.cs(cs)$base_scalars) else lapply(seq_len(.dg_dim(cs)) - 1L, function(i) base_scalar(cs, i))
#' @rdname diffgeom
#' @export
base_vector <- function(cs, coord_index) { if (.is_bd(cs)) return(.bd_base(cs)); .nc_state$any <- TRUE; new_sym("BaseVectorField", list(cs, .int(coord_index))) }
#' @rdname diffgeom
#' @export
base_vectors <- function(cs) if (.is_coordsys(cs)) .pytuple(.cs(cs)$base_vectors) else lapply(seq_len(.dg_dim(cs)) - 1L, function(i) base_vector(cs, i))
#' @rdname diffgeom
#' @export
base_oneform <- function(cs, coord_index) Differential(base_scalar(cs, coord_index))
#' @rdname diffgeom
#' @export
base_oneforms <- function(cs) lapply(seq_len(.dg_dim(cs)) - 1L, function(i) base_oneform(cs, i))

.dg_idx <- function(f) .int_val(f$args[[2L]])
.dg_cs <- function(f) f$args[[1L]]
.dg_symname <- function(f) .dg_cs(f)$extra$symbols[[.dg_idx(f) + 1L]]$name
register_str_printer("BaseScalarField", function(e) .dg_symname(e))
register_str_printer("BaseVectorField", function(e) paste0("e_", .dg_symname(e)))
register_latex_printer("BaseScalarField", function(e) .latex(Symbol(.dg_symname(e))))
register_latex_printer("BaseVectorField", function(e) paste0("\\partial_{", .latex(Symbol(.dg_symname(e))), "}"))

.facts_hooks$BaseScalarField <- function(e) c(commutative = TRUE)
for (.t in c("BaseVectorField", "Differential", "TensorProduct", "WedgeProduct", "Commutator", "LieDerivative",
             "BaseCovarDerivativeOp", "CovarDerivativeOp"))
  assign(.t, function(e) c(commutative = FALSE), envir = .facts_hooks)
rm(.t)

#' @rdname diffgeom
#' @export
Differential <- function(form_field) {
  form_field <- as_sym(form_field)
  if (contravariant_order(form_field)) stop("A vector field was supplied as an argument to Differential.", call. = FALSE)
  if (form_field$type == "Differential") return(S_Zero)
  .nc_state$any <- TRUE
  new_sym("Differential", list(form_field))
}
register_str_printer("Differential", function(e) {
  f <- e$args[[1L]]
  if (f$type == "BaseScalarField") paste0("d", .dg_symname(f)) else paste0("d(", .str(f), ")")
})

.DG_NONE <- new_sym("DGNone", list(), name = "None", key = "DGNone")
.dg_isnone <- function(z) is_sym(z) && identical(z$key, "DGNone")

#' @rdname diffgeom
#' @export
TensorProduct <- function(...) {
  args <- lapply(list(...), as_sym)
  ords <- vapply(args, function(m) covariant_order(m) + contravariant_order(m), 0)
  scalar <- .mul(args[ords == 0])
  multi <- args[ords != 0]
  if (!length(multi)) return(scalar)
  if (length(multi) == 1L) return(.mul(list(scalar, multi[[1L]])))
  .nc_state$any <- TRUE
  .mul(list(scalar, new_sym("TensorProduct", multi)))
}
#' @rdname diffgeom
#' @export
WedgeProduct <- function(...) {
  args <- lapply(list(...), as_sym)
  ords <- vapply(args, function(m) covariant_order(m) + contravariant_order(m), 0)
  scalar <- .mul(args[ords == 0])
  multi <- args[ords != 0]
  if (!length(multi)) return(scalar)
  if (length(multi) == 1L) return(.mul(list(scalar, multi[[1L]])))
  .nc_state$any <- TRUE
  .mul(list(scalar, new_sym("WedgeProduct", multi)))
}

.dg_vecs_of <- function(v) Filter(function(z) z$type == "BaseVectorField", atoms(v, "BaseVectorField"))
.dg_find_coords <- function(e) {
  fs <- atoms(e, c("BaseScalarField", "BaseVectorField"))
  cs <- lapply(fs, .dg_cs)
  cs[!duplicated(vapply(cs, `[[`, "", "key"))]
}
# expr.coeff(b) for a linear combination of base vectors
.dg_coeff <- function(e, b) {
  terms <- if (e$type == "Add") e$args else list(e)
  out <- list()
  for (t in terms) {
    fs <- if (t$type == "Mul") t$args else list(t)
    hit <- vapply(fs, function(z) identical(z$key, b$key), NA)
    if (sum(hit) == 1L) out[[length(out) + 1L]] <- .mul(fs[!hit])
  }
  .add(out)
}

#' @rdname diffgeom
#' @export
Commutator <- function(v1, v2) {
  v1 <- as_sym(v1); v2 <- as_sym(v2)
  if (covariant_order(v1) || contravariant_order(v1) != 1 || covariant_order(v2) || contravariant_order(v2) != 1)
    stop("Only commutators of vector fields are supported.", call. = FALSE)
  if (identical(v1$key, v2$key)) return(S_Zero)
  cs <- .dg_find_coords(.add(list(v1, v2)))
  cs <- c(.dg_find_coords(v1), .dg_find_coords(v2)); cs <- cs[!duplicated(vapply(cs, `[[`, "", "key"))]
  if (length(cs) == 1L) {
    if (all(vapply(list(v1, v2), function(v) v$type == "BaseVectorField", NA))) return(S_Zero)
    b1 <- .dg_vecs_of(v1); b2 <- .dg_vecs_of(v2)
    c1 <- lapply(b1, function(b) .dg_coeff(expand(v1), b)); c2 <- lapply(b2, function(b) .dg_coeff(expand(v2), b))
    res <- list()
    for (i in seq_along(b1)) for (j in seq_along(b2))
      res[[length(res) + 1L]] <- .sub(.mul(list(c1[[i]], rcall(b1[[i]], c2[[j]]), b2[[j]])),
                                      .mul(list(c2[[j]], rcall(b2[[j]], c1[[i]]), b1[[i]])))
    return(.add(res))
  }
  .nc_state$any <- TRUE
  new_sym("Commutator", list(v1, v2))
}

#' @rdname diffgeom
#' @export
LieDerivative <- function(v_field, expr) {
  v_field <- as_sym(v_field); expr <- as_sym(expr)
  if (contravariant_order(v_field) != 1 || covariant_order(v_field))
    stop("Lie derivatives are defined only with respect to vector fields. The supplied argument was not a vector field.", call. = FALSE)
  if (covariant_order(expr) > 0) { .nc_state$any <- TRUE; return(new_sym("LieDerivative", list(v_field, expr))) }
  if (length(atoms(expr, "BaseVectorField"))) return(Commutator(v_field, expr))
  rcall(v_field, expr)
}

#' @rdname diffgeom
#' @export
BaseCovarDerivativeOp <- function(coord_sys, index, christoffel) {
  .nc_state$any <- TRUE
  new_sym("BaseCovarDerivativeOp", list(coord_sys, .int(index)), extra = list(ch = christoffel))
}
#' @rdname diffgeom
#' @export
CovarDerivativeOp <- function(wrt, christoffel) {
  wrt <- as_sym(wrt)
  cs <- lapply(.dg_vecs_of(wrt), .dg_cs)
  if (length(unique(vapply(cs, `[[`, "", "key"))) > 1L) stop("NotImplementedError", call. = FALSE)
  if (contravariant_order(wrt) != 1 || covariant_order(wrt))
    stop("Covariant derivatives are defined only with respect to vector fields. The supplied argument was not a vector field.", call. = FALSE)
  .nc_state$any <- TRUE
  new_sym("CovarDerivativeOp", list(wrt), extra = list(ch = christoffel))
}

# ---------------------------------------------------------------------------
# calls

#' @rdname diffgeom
#' @export
rcall <- function(expr, ...) .dg_rcall(as_sym(expr), lapply(list(...), function(z) if (is.null(z)) .DG_NONE else z))

.DG_CALLABLE <- c("BaseScalarField", "BaseVectorField", "Differential", "TensorProduct", "WedgeProduct", "Commutator",
                  "LieDerivative", "BaseCovarDerivativeOp", "CovarDerivativeOp", "Lambda")
.dg_rcall <- function(e, args) {
  if (e$type %in% .DG_CALLABLE) return(.dg_call(e, args))
  if (!length(e$args) || e$type %in% c("Symbol", "Dummy")) return(e)
  sym_rebuild(e, lapply(e$args, .dg_rcall, args = args))
}

.dg_call <- function(e, args) {
  switch(e$type,
    Lambda = do.call(apply_lambda, c(list(e), args)),
    BaseScalarField = {
      pt <- args[[1L]]
      if (length(args) != 1L || !is_sym(pt) || pt$type != "Point") return(e)
      co <- .dg_point_coords(pt, .dg_cs(e))
      doit(simplify(co$data[[.dg_idx(e) + 1L]]))
    },
    BaseVectorField = .dg_vector_call(e, args[[1L]]),
    Differential = .dg_differential_call(e, args),
    TensorProduct = .dg_tensor_call(e, args),
    WedgeProduct = .dg_wedge_call(e, args),
    Commutator = .sub(rcall(e$args[[1L]], rcall(e$args[[2L]], args[[1L]])), rcall(e$args[[2L]], rcall(e$args[[1L]], args[[1L]]))),
    LieDerivative = {
      v <- e$args[[1L]]; ex <- e$args[[2L]]
      lead <- rcall(v, do.call(rcall, c(list(ex), args)))
      # SymPy multiplies the arguments with the commutator in place
      rest <- .add(lapply(seq_along(args), function(i) {
        a <- args; a[[i]] <- Commutator(v, args[[i]])
        .mul(a)
      }))
      .sub(lead, rest)
    },
    BaseCovarDerivativeOp = .dg_bcov_call(e, args[[1L]]),
    CovarDerivativeOp = {
      vecs <- .dg_vecs_of(e$args[[1L]])
      ops <- lapply(vecs, function(v) BaseCovarDerivativeOp(.dg_cs(v), .dg_idx(v), e$extra$ch))
      rcall(subs(e$args[[1L]], Map(list, vecs, ops)), args[[1L]])
    })
}

.dg_vector_call <- function(e, sf) {
  if (.dg_isnone(sf)) return(e)
  sf <- as_sym(sf)
  if (covariant_order(sf) || contravariant_order(sf))
    stop("Only scalar fields can be supplied as arguments to vector fields.", call. = FALSE)
  cs <- .dg_cs(e)
  bs <- atoms(sf, "BaseScalarField")
  dv <- cs$extra$dummy
  dfs <- lapply(seq_along(bs) - 1L, function(i) fcall(paste0("_#_", i), list(dv)))
  dres <- .diff1(subs(sf, Map(list, bs, dfs), simultaneous = TRUE), dv)
  coords <- cs$extra$symbols
  dder <- lapply(dfs, function(f) .diff1(f, dv))
  dsub <- lapply(bs, function(b) .mget(jacobian_dg(cs, .dg_cs(b), coords), .dg_idx(b) + 1L, .dg_idx(e) + 1L))
  dres <- subs(dres, Map(list, dder, dsub), simultaneous = TRUE)
  res <- subs(dres, Map(list, dfs, bs), simultaneous = TRUE)
  res <- subs(res, Map(list, coords, base_scalars(cs)), simultaneous = TRUE)
  doit(res)
}

.dg_differential_call <- function(e, vf) {
  for (a in vf) if (!.dg_isnone(a) && (contravariant_order(a) != 1 || covariant_order(a)))
    stop("The arguments supplied to Differential should be vector fields or Nones.", call. = FALSE)
  k <- length(vf); f <- e$args[[1L]]
  if (k == 1L) return(if (!.dg_isnone(vf[[1L]]) && !is_zero_num(vf[[1L]])) rcall(vf[[1L]], f) else e)
  ret <- list()
  for (i in seq_len(k)) {
    t <- rcall(vf[[i]], do.call(rcall, c(list(f), vf[-i])))
    ret[[length(ret) + 1L]] <- .mul(list(.int((-1)^(i - 1L)), t))
    if (i < k) for (j in (i + 1L):k) {
      cm <- Commutator(vf[[i]], vf[[j]])
      if (!is_zero_num(cm)) {
        t2 <- do.call(rcall, c(list(f), list(cm), vf[-c(i, j)]))
        ret[[length(ret) + 1L]] <- .mul(list(.int((-1)^(i + j - 2L)), t2))
      }
    }
  }
  .add(ret)
}

.dg_order <- function(z) covariant_order(z) + contravariant_order(z)
.dg_tensor_call <- function(e, fields) {
  tot <- .dg_order(e)
  if (length(fields) != tot) fields <- c(fields, rep(list(.DG_NONE), tot - length(fields)))
  ords <- vapply(e$args, .dg_order, 0)
  ends <- cumsum(ords); starts <- c(1L, utils::head(ends, -1L) + 1L)
  mult <- lapply(seq_along(e$args), function(i) {
    fs <- if (ends[i] >= starts[i]) fields[starts[i]:ends[i]] else list()
    do.call(rcall, c(list(e$args[[i]]), fs))
  })
  do.call(TensorProduct, mult)
}

.perm_sign <- function(p) {
  s <- 1L; p <- unlist(p)
  for (i in seq_along(p)) for (j in seq_along(p)) if (i < j && p[i] > p[j]) s <- -s
  s
}
.dg_wedge_call <- function(e, fields) {
  ords <- vapply(e$args, .dg_order, 0)
  mul <- .inv(.mul(lapply(ords, function(o) fcall("factorial", list(.int(o))))))
  perms <- .it_perms(fields)
  signs <- vapply(.it_perms(seq_along(fields) - 1L), .perm_sign, 0L)
  tp <- do.call(TensorProduct, e$args)
  .mul(list(mul, .add(Map(function(p, s) .mul(list(.dg_tensor_call(tp, p), .int(s))), perms, signs))))
}

.dg_bcov_call <- function(e, field) {
  field <- as_sym(field)
  if (covariant_order(field) != 0) stop("NotImplementedError", call. = FALSE)
  cs <- .dg_cs(e); idx <- .dg_idx(e); ch <- e$extra$ch
  field <- vectors_in_basis(field, cs)
  wv <- base_vector(cs, idx); ws <- base_scalar(cs, idx)
  vecs <- atoms(field, "BaseVectorField")
  dfs <- lapply(seq_along(vecs) - 1L, function(i) fcall(paste0("_#_", i), list(ws)))
  dres <- subs(field, Map(list, vecs, dfs), simultaneous = TRUE)
  dres <- rcall(wv, dres)
  dres <- subs(dres, Map(list, dfs, vecs), simultaneous = TRUE)
  derivs <- lapply(vecs, function(v) .add(lapply(seq_len(.dg_dim(.dg_cs(v))) - 1L, function(k)
    .mul(list(arr_get(ch, c(k + 1L, .dg_idx(wv) + 1L, .dg_idx(v) + 1L)), base_vector(.dg_cs(v), k))))))
  tosubs <- lapply(dfs, function(d) rcall(wv, d))
  res <- subs(dres, Map(list, tosubs, derivs), simultaneous = TRUE)
  res <- subs(res, Map(list, dfs, vecs), simultaneous = TRUE)
  doit(res)
}

# ---------------------------------------------------------------------------
# orders and metric helpers

#' @rdname diffgeom
#' @export
contravariant_order <- function(expr, .strict = FALSE) {
  expr <- as_sym(expr)
  switch(expr$type,
    Add = { o <- vapply(expr$args, contravariant_order, 0)
            if (length(unique(o)) != 1L) stop("Misformed expression containing contravariant fields of varying order.", call. = FALSE)
            o[1L] },
    Mul = { o <- vapply(expr$args, contravariant_order, 0); nz <- o[o != 0]
            if (length(nz) > 1L) stop("Misformed expression containing multiplication between vectors.", call. = FALSE)
            if (!length(nz)) 0 else nz[1L] },
    Pow = { if (covariant_order(expr$args[[1L]]) || covariant_order(expr$args[[2L]]))
              stop("Misformed expression containing a power of a vector.", call. = FALSE)
            0 },
    BaseVectorField = 1,
    TensorProduct = , WedgeProduct = sum(vapply(expr$args, contravariant_order, 0)),
    if (!.strict || length(atoms(expr, "BaseScalarField"))) 0 else -1)
}

#' @rdname diffgeom
#' @export
covariant_order <- function(expr, .strict = FALSE) {
  expr <- as_sym(expr)
  switch(expr$type,
    Add = { o <- vapply(expr$args, covariant_order, 0)
            if (length(unique(o)) != 1L) stop("Misformed expression containing form fields of varying order.", call. = FALSE)
            o[1L] },
    Mul = { o <- vapply(expr$args, covariant_order, 0); nz <- o[o != 0]
            if (length(nz) > 1L) stop("Misformed expression containing multiplication between forms.", call. = FALSE)
            if (!length(nz)) 0 else nz[1L] },
    Pow = { if (covariant_order(expr$args[[1L]]) || covariant_order(expr$args[[2L]]))
              stop("Misformed expression containing a power of a form.", call. = FALSE)
            0 },
    Differential = covariant_order(expr$args[[1L]]) + 1,
    TensorProduct = , WedgeProduct = sum(vapply(expr$args, covariant_order, 0)),
    if (!.strict || length(atoms(expr, "BaseScalarField"))) 0 else -1)
}

#' @rdname diffgeom
#' @export
vectors_in_basis <- function(expr, to_sys) {
  expr <- as_sym(expr)
  vecs <- atoms(expr, "BaseVectorField")
  nv <- lapply(vecs, function(v) {
    cs <- .dg_cs(v)
    jac <- jacobian_dg(cs, to_sys, base_scalars(cs))
    col <- t(jac) %*% Matrix(base_vectors(to_sys))
    col$data[[.dg_idx(v) + 1L]]
  })
  subs(expr, Map(list, vecs, nv), simultaneous = TRUE)
}

#' @rdname diffgeom
#' @export
twoform_to_matrix <- function(expr) {
  expr <- as_sym(expr)
  if (covariant_order(expr) != 2 || contravariant_order(expr)) stop("The input expression is not a two-form.", call. = FALSE)
  cs <- .dg_find_coords(expr)
  if (length(cs) != 1L) stop("The input expression concerns more than one coordinate systems, hence there is no unambiguous way to choose a coordinate system for the matrix.", call. = FALSE)
  vs <- base_vectors(cs[[1L]])
  ex <- expand(expr)
  Matrix(lapply(vs, function(v2) lapply(vs, function(v1) rcall(ex, v1, v2))))
}

#' @rdname diffgeom
#' @export
metric_to_Christoffel_1st <- function(expr) {
  m <- twoform_to_matrix(expr)
  if (!isTRUE(is_symmetric(m))) stop("The two-form representing the metric is not symmetric.", call. = FALSE)
  cs <- .dg_find_coords(expr)[[1L]]
  dm <- lapply(base_vectors(cs), function(d) .mat_map(m, function(z) rcall(d, z)))
  n <- .dg_dim(cs); ix <- seq_len(n)
  Array(lapply(ix, function(i) lapply(ix, function(j) lapply(ix, function(k)
    .div(.sub(.add(list(.mget(dm[[k]], i, j), .mget(dm[[j]], i, k))), .mget(dm[[i]], j, k)), S_Two)))))
}

#' @rdname diffgeom
#' @export
metric_to_Christoffel_2nd <- function(expr) {
  ch1 <- metric_to_Christoffel_1st(expr)
  cs <- .dg_find_coords(expr)[[1L]]
  n <- .dg_dim(cs); ix <- seq_len(n)
  m <- twoform_to_matrix(expr)
  sf <- unique(unlist(lapply(m$data, function(z) atoms(z, "BaseScalarField")), recursive = FALSE))
  sf <- sf[!duplicated(vapply(sf, `[[`, "", "key"))]
  dums <- cs$extra$symbols[as.integer(vapply(sf, .dg_idx, 0)) + 1L]
  mi <- .mat_map(inverse(.mat_map(m, function(z) subs(z, Map(list, sf, dums), simultaneous = TRUE))),
                 function(z) subs(z, Map(list, dums, sf), simultaneous = TRUE))
  Array(lapply(ix, function(i) lapply(ix, function(j) lapply(ix, function(k)
    .add(lapply(ix, function(l) .mul(list(.mget(mi, i, l), arr_get(ch1, c(l, j, k))))))))))
}

#' @rdname diffgeom
#' @export
metric_to_Riemann_components <- function(expr) {
  ch2 <- metric_to_Christoffel_2nd(expr)
  cs <- .dg_find_coords(expr)[[1L]]
  n <- .dg_dim(cs); ix <- seq_len(n); bv <- base_vectors(cs)
  dch <- function(i, j, k, d) rcall(bv[[d]], arr_get(ch2, c(i, j, k)))
  Array(lapply(ix, function(rho) lapply(ix, function(sig) lapply(ix, function(mu) lapply(ix, function(nu) {
    a <- .sub(dch(rho, sig, nu, mu), dch(rho, sig, mu, nu))
    b <- .add(lapply(ix, function(l) .sub(.mul(list(arr_get(ch2, c(rho, l, mu)), arr_get(ch2, c(l, sig, nu)))),
                                          .mul(list(arr_get(ch2, c(rho, l, nu)), arr_get(ch2, c(l, sig, mu)))))))
    .add(list(a, b))
  })))))
}

#' @rdname diffgeom
#' @export
metric_to_Ricci_components <- function(expr) {
  rm4 <- metric_to_Riemann_components(expr)
  cs <- .dg_find_coords(expr)[[1L]]
  ix <- seq_len(.dg_dim(cs))
  Array(lapply(ix, function(i) lapply(ix, function(j) .add(lapply(ix, function(k) arr_get(rm4, c(k, i, k, j)))))))
}

#' @rdname diffgeom
#' @export
intcurve_series <- function(vector_field, param, start_point, n = 6, coord_sys = NULL, coeffs = FALSE) {
  vector_field <- as_sym(vector_field)
  if (contravariant_order(vector_field) != 1 || covariant_order(vector_field)) stop("The supplied field was not a vector field.", call. = FALSE)
  iter <- function(sf, i) { for (k in seq_len(i)) sf <- rcall(vector_field, sf); sf }
  cs <- if (is.null(coord_sys)) start_point$args[[1L]] else coord_sys
  terms <- lapply(base_scalars(cs), function(cf) lapply(seq_len(n) - 1L, function(i)
    .div(.mul(list(.pow(param, .int(i)), rcall(iter(cf, i), start_point))), fcall("factorial", list(.int(i))))))
  if (coeffs) return(lapply(seq_len(n), function(i) Matrix(lapply(terms, `[[`, i))))
  Matrix(lapply(terms, .add))
}

#' @rdname diffgeom
#' @export
intcurve_diffequ <- function(vector_field, param, start_point, coord_sys = NULL) {
  vector_field <- as_sym(vector_field)
  if (contravariant_order(vector_field) != 1 || covariant_order(vector_field)) stop("The supplied field was not a vector field.", call. = FALSE)
  cs <- if (is.null(coord_sys)) start_point$args[[1L]] else coord_sys
  gam <- lapply(seq_len(.dg_dim(start_point$args[[1L]])) - 1L, function(i) fcall(paste0("f_", i), list(param)))
  ap <- point(cs, gam)
  cfs <- base_scalars(cs)
  eqs <- lapply(cfs, function(cf) simplify(.sub(.diff1(rcall(cf, ap), param), rcall(rcall(vector_field, cf), ap))))
  ini <- lapply(cfs, function(cf) simplify(.sub(subs(rcall(cf, ap), param, S_Zero), rcall(cf, start_point))))
  list(eqs, ini)
}

# ---------------------------------------------------------------------------
# sympy.diffgeom.rn

#' @rdname diffgeom
#' @export
diffgeom_rn <- function() {
  x <- Symbol("x", real = TRUE); y <- Symbol("y", real = TRUE)
  r <- Symbol("rho", nonnegative = TRUE); th <- Symbol("theta", nonnegative = TRUE)
  R2 <- Manifold("R^2", 2); R2o <- Patch("origin", R2)
  rel2 <- list(list(c("rectangular", "polar"), list(x, y), list(sqrt(x^2 + y^2), atan2(y, x))),
               list(c("polar", "rectangular"), list(r, th), list(r * cos(th), r * sin(th))))
  R2_r <- CoordSystem("rectangular", R2o, list(x, y), rel2)
  R2_p <- CoordSystem("polar", R2o, list(r, th), rel2)
  z <- Symbol("z", real = TRUE)
  rho <- Symbol("rho", nonnegative = TRUE); psi <- Symbol("psi", nonnegative = TRUE)
  rr <- Symbol("r", nonnegative = TRUE); phi <- Symbol("phi", nonnegative = TRUE)
  R3 <- Manifold("R^3", 3); R3o <- Patch("origin", R3)
  rel3 <- list(
    list(c("rectangular", "cylindrical"), list(x, y, z), list(sqrt(x^2 + y^2), atan2(y, x), z)),
    list(c("cylindrical", "rectangular"), list(rho, psi, z), list(rho * cos(psi), rho * sin(psi), z)),
    list(c("rectangular", "spherical"), list(x, y, z), list(sqrt(x^2 + y^2 + z^2), acos(z / sqrt(x^2 + y^2 + z^2)), atan2(y, x))),
    list(c("spherical", "rectangular"), list(rr, th, phi), list(rr * sin(th) * cos(phi), rr * sin(th) * sin(phi), rr * cos(th))),
    list(c("cylindrical", "spherical"), list(rho, psi, z), list(sqrt(rho^2 + z^2), acos(z / sqrt(rho^2 + z^2)), psi)),
    list(c("spherical", "cylindrical"), list(rr, th, phi), list(rr * sin(th), phi, rr * cos(th))))
  list(R2 = R2, R2_origin = R2o, R2_r = R2_r, R2_p = R2_p, R3 = R3, R3_origin = R3o,
       R3_r = CoordSystem("rectangular", R3o, list(x, y, z), rel3),
       R3_c = CoordSystem("cylindrical", R3o, list(rho, psi, z), rel3),
       R3_s = CoordSystem("spherical", R3o, list(rr, th, phi), rel3))
}
