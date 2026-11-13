# DomainMatrix, after sympy/polys/matrices/domainmatrix.py, rref.py,
# domainscalar.py, eigen.py and normalforms.py. A DomainMatrix wraps a DDM
# (dense, R/polys-ddm.R) or an SDM (sparse, R/polys-sdm.R) of elements of a
# domain (R/polys-domelem.R). Methods are called with `$` as in SymPy, with
# 0-based indices; `[` indexes from 1 as R does.
# Limitation: DFM (dense matrices over python-flint) is not available; the
# dense format is always DDM, as SymPy uses when flint is not installed.

# ---------------------------------------------------------------------------
# Python str() and repr() of the values the methods return

.dmpy_str <- function(v) {
  if (is.null(v)) return("None")
  if (inherits(v, c("DDM", "SDM"))) return(.dm_str(v))
  if (inherits(v, c("DomainMatrix", "DomainScalar", "pydict", "pytuple", "sym_domain", "PowerSeriesElement",
                    "PowerSeriesRing")))
    return(.dmpy_repr(v))
  if (inherits(v, c("bigz", "bigq", "PolyElement", "FracElement", "EXElement", "sym_modint", "sym_gaussian", "ANP")))
    return(.de_str(v))
  if (is_sym(v)) return(gsub("^", "**", sstr(v), fixed = TRUE))
  if (inherits(v, "symmatrix")) return(gsub("^", "**", sstr(v), fixed = TRUE))
  if (is.character(v) && length(v) == 1L) return(v)
  .dmpy_repr(v)
}
.dmpy_repr <- function(v) {
  if (is.null(v)) return("None")
  if (inherits(v, "DomainMatrix")) return(.dm_dm_repr(v))
  if (inherits(v, c("DDM", "SDM"))) return(.dm_repr(v))
  if (inherits(v, "DomainScalar")) return(.dmpy_repr(v$element))
  if (inherits(v, "sym_domain")) return(format(v))
  if (inherits(v, "Poly")) return(paste(utils::capture.output(print(v)), collapse = "\n"))
  if (inherits(v, "PowerSeriesElement") || inherits(v, "PowerSeriesRing")) return(format(v))
  if (inherits(v, "pytuple")) {
    s <- vapply(v, .dmpy_repr, "")
    return(paste0("(", paste(s, collapse = ", "), if (length(s) == 1L) ",", ")"))
  }
  if (inherits(v, "pydict")) {
    if (!length(v$keys)) return("{}")
    return(paste0("{", paste(vapply(seq_along(v$keys), function(i) paste0(.dmpy_repr(v$keys[[i]]), ": ", .dmpy_repr(v$values[[i]])), ""),
                             collapse = ", "), "}"))
  }
  if (inherits(v, c("bigz", "bigq", "PolyElement", "FracElement", "EXElement", "sym_modint", "sym_gaussian", "ANP")))
    return(.de_repr(v))
  if (is_sym(v) || inherits(v, "symmatrix")) return(gsub("^", "**", sstr(v), fixed = TRUE))
  if (is.logical(v) && length(v) == 1L) return(if (is.na(v)) "None" else if (v) "True" else "False")
  if (is.character(v) && length(v) == 1L) return(paste0("'", v, "'"))
  if (is.integer(v) && length(v) == 1L) return(as.character(v))
  if (is.double(v) && length(v) == 1L) return(if (v == round(v) && abs(v) < 1e15) paste0(format(v, scientific = FALSE), ".0") else as.character(v))
  if (is.list(v) && !is.object(v)) return(paste0("[", paste(vapply(v, .dmpy_repr, ""), collapse = ", "), "]"))
  paste(format(v), collapse = " ")
}
# str() of a DDM or SDM
.dm_str <- function(a) {
  if (inherits(a, "DDM")) {
    rows <- .subset2(a, "rows")
    return(paste0("[", paste(vapply(rows, function(r) paste0("[", paste(vapply(r, .dmpy_str, ""), collapse = ", "), "]"), ""),
                             collapse = ", "), "]"))
  }
  d <- .subset2(a, "d")
  paste0("{", paste(vapply(names(d), function(i) paste0(i, ": {", paste(vapply(names(d[[i]]), function(j)
    paste0(j, ": ", .dmpy_str(d[[i]][[j]])), ""), collapse = ", "), "}"), ""), collapse = ", "), "}")
}
# repr() of a DDM or SDM
.dm_repr <- function(a) {
  s <- .subset2(a, "shape"); K <- .subset2(a, "domain")
  body <- if (inherits(a, "DDM")) {
    rows <- .subset2(a, "rows")
    paste0("[", paste(vapply(rows, function(r) paste0("[", paste(vapply(r, .dmpy_repr, ""), collapse = ", "), "]"), ""), collapse = ", "), "]")
  } else {
    d <- .subset2(a, "d")
    paste0("{", paste(vapply(names(d), function(i) paste0(i, ": {", paste(vapply(names(d[[i]]), function(j)
      paste0(j, ": ", .dmpy_repr(d[[i]][[j]])), ""), collapse = ", "), "}"), ""), collapse = ", "), "}")
  }
  gsub("**", "^", paste0(if (inherits(a, "DDM")) "DDM(" else "SDM(", body, ", ", .shape_str(s), ", ", format(K), ")"), fixed = TRUE)
}
.dm_dm_repr <- function(A) {
  rep <- .subset2(A, "rep")
  gsub("**", "^", paste0("DomainMatrix(", .dm_str(rep), ", ", .shape_str(.subset2(rep, "shape")), ", ", format(.subset2(rep, "domain")), ")"),
       fixed = TRUE)
}

# Ops of DDM and SDM: `+` and `-` with the same class, `*` with an element
.dm_ops <- function(op, e1, e2, cls) {
  if (missing(e2)) {
    if (op == "-") return(e1$neg())
    if (op == "+") return(e1)
  }
  is_m <- function(z) inherits(z, cls)
  switch(op,
    `+` = { if (!(is_m(e1) && is_m(e2))) stop(.py_err("TypeError", "unsupported operand type(s) for +")); e1$add(e2) },
    `-` = { if (!(is_m(e1) && is_m(e2))) stop(.py_err("TypeError", "unsupported operand type(s) for -")); e1$sub(e2) },
    `*` = if (is_m(e1) && is_m(e2)) { if (cls == "SDM") e1$matmul(e2) else stop(.py_err("TypeError", "unsupported operand type(s) for *")) }
          else if (is_m(e1)) e1$mul(.de_convert(.subset2(e1, "domain"), e2)) else e2$rmul(.de_convert(.subset2(e2, "domain"), e1)),
    `==` = is_m(e1) && is_m(e2) && .dm_rep_eq(e1, e2),
    `!=` = !(is_m(e1) && is_m(e2) && .dm_rep_eq(e1, e2)),
    stop("operator ", op, " is not defined for ", cls, call. = FALSE))
}
.dm_rep_eq <- function(a, b) {
  if (!identical(class(a), class(b)) || !.de_same(.subset2(a, "domain"), .subset2(b, "domain"))) return(FALSE)
  if (!identical(.subset2(a, "shape"), .subset2(b, "shape"))) return(FALSE)
  la <- unlist(.dm_rows(a), recursive = FALSE); lb <- unlist(.dm_rows(b), recursive = FALSE)
  all(mapply(.de_eq, la, lb))
}
.dm_rows <- function(rep) if (inherits(rep, "DDM")) .subset2(rep, "rows") else .sdm_to_list(rep)

# ---------------------------------------------------------------------------
# DomainMatrix

.dm_from_rep <- function(rep) structure(list(rep = rep), class = "DomainMatrix")

# an entry of DM() or from_list(): a pair c(p, q) is domain(p, q)
.dm_elem <- function(K, e) {
  if (is.list(e) && !is.object(e) && length(e) == 2L) e <- unlist(e)
  if (is.numeric(e) && length(e) == 2L) {
    k <- .de_kind(K)
    if (k %in% c("ZZ_I", "QQ_I")) return(.dom_convert(.dom(K), e[1L], e[2L]))
    return(.de_convert(K, gmp::as.bigq(e[1L]) / gmp::as.bigq(e[2L])))
  }
  if ((inherits(e, "bigz") || inherits(e, "bigq")) && length(e) == 2L) return(.de_convert(K, gmp::as.bigq(e[1L]) / gmp::as.bigq(e[2L])))
  .de_convert(K, e)
}
.dm_rowslist <- function(rows) {
  if (is.matrix(rows)) return(lapply(seq_len(nrow(rows)), function(i) as.list(rows[i, ])))
  lapply(rows, function(r) if (is.list(r) && !is.object(r)) r else as.list(r))
}

#' Matrices over polynomial domains
#'
#' `DomainMatrix(rows, shape, domain)` is SymPy's
#' `sympy.polys.matrices.DomainMatrix`: a matrix whose entries are elements
#' of a domain such as `ZZ`, `QQ`, `GF(p)`, `ZZ_I`, `QQ[x]`, `QQ(x)`,
#' `QQ<sqrt(2)>`, `RR` or `EX`. `rows` is a list of rows (each a list of
#' elements) for a dense matrix, or a nested named list
#' `list(i = list(j = element))` with 0-based indices for a sparse one. `DM(rows,
#' domain)` converts plain values with the domain (`c(1, 2)` is the pair
#' `QQ(1, 2)` or `ZZ_I(1, 2)`). The matrix methods of SymPy are called with
#' `$`: `A$rref()`, `A$rref_den()`, `A$det()`, `A$inv()`, `A$charpoly()`,
#' `A$nullspace()`, `A$lu()`, `A$fflu()`, `A$convert_to(QQ)`, `A$to_Matrix()`
#' and so on, with 0-based indices and Python tuples and lists as results.
#' `A + B`, `A - B`, `A * B` (matrix product or scalar multiple), `A^n` and
#' `A / x` work as in SymPy; `A[i, j]` indexes from 1 and returns a
#' [DomainScalar()] or a submatrix. `DomainMatrix_eye()`, `DomainMatrix_zeros()`,
#' `DomainMatrix_ones()`, `DomainMatrix_diag()`, `DomainMatrix_from_list()`,
#' `DomainMatrix_from_list_sympy()`, `DomainMatrix_from_dict_sympy()`,
#' `DomainMatrix_from_Matrix()`, `DomainMatrix_from_rep()`,
#' `DomainMatrix_from_list_flat()`, `DomainMatrix_from_dod()` and
#' `DomainMatrix_from_dok()` are the class methods. Matrices over flint
#' (`DFM`) are not available.
#'
#' @param rows A list of rows, or a nested named list of nonzero entries.
#' @param shape The number of rows and columns.
#' @param domain The domain of the entries.
#' @param fmt `"dense"` or `"sparse"` to force a format.
#' @param x A `DomainMatrix`.
#' @param name A method name.
#' @param i,j Row and column indices (from 1).
#' @param drop Unused.
#' @param e1,e2 Matrices or scalars.
#' @param rep A `DDM` or `SDM`.
#' @param nrows,ncols The shape.
#' @param M A [Matrix()].
#' @param diagonal A list of diagonal elements.
#' @param elements,data A flat list of elements and its layout.
#' @param dod,dok A dict of dicts or a dict of keys.
#' @param elemsdict A nested named list of entries.
#' @param ... Options passed to [construct_domain()], or unused.
#' @return A `DomainMatrix`, or the value of the method called.
#' @examples
#' A <- DM(list(list(1, 2), list(3, 4)), ZZ)
#' A
#' A$det()
#' A$rref()
#' A$charpoly()
#' A$convert_to(QQ)$inv()
#' A * A
#' A[1, 2]
#' @export
DomainMatrix <- function(rows, shape, domain, fmt = NULL) {
  if (inherits(rows, c("DDM", "SDM"))) stop(.py_err("TypeError", "Use from_rep to initialize from SDM/DDM"))
  if (!is.list(rows) && !is.matrix(rows)) stop(.py_err("TypeError", "Input should be list-of-lists or dict-of-dicts"))
  num <- function(e) if (is.numeric(e) && !is.object(e) && length(e) == 1L) .de_convert(domain, e) else e
  if (is.list(rows)) rows <- lapply(rows, function(r) if (is.list(r)) lapply(r, num) else r)
  rep <- if (!is.null(names(rows))) SDM(rows, shape, domain) else DDM(.dm_rowslist(rows), shape, domain)
  if (!is.null(fmt)) {
    rep <- if (fmt == "sparse") rep$to_sdm() else if (fmt == "dense") rep$to_ddm()
           else stop(.py_err("ValueError", "fmt should be 'sparse' or 'dense'"))
  }
  .dm_from_rep(rep)
}
#' @rdname DomainMatrix
#' @export
DM <- function(rows, domain) DomainMatrix_from_list(rows, domain)
#' @rdname DomainMatrix
#' @export
DomainMatrix_from_rep <- function(rep) {
  if (!inherits(rep, c("DDM", "SDM"))) stop(.py_err("TypeError", "rep should be of type DDM or SDM"))
  .dm_from_rep(rep)
}
#' @rdname DomainMatrix
#' @export
DomainMatrix_from_list <- function(rows, domain) {
  rows <- .dm_rowslist(rows)
  nr <- length(rows); nc <- if (nr) length(rows[[1L]]) else 0L
  DomainMatrix(lapply(rows, function(r) lapply(r, function(e) .dm_elem(domain, e))), c(nr, nc), domain)
}
#' @rdname DomainMatrix
#' @export
DomainMatrix_from_list_sympy <- function(nrows, ncols, rows, ...) {
  rows <- .dm_rowslist(rows)
  items <- lapply(unlist(rows, recursive = FALSE), as_sym)
  r <- .dm_get_domain(items, ...)
  K <- r[[1L]]; el <- r[[2L]]
  DomainMatrix(lapply(seq_len(nrows), function(i) el[(i - 1L) * ncols + seq_len(ncols)]), c(nrows, ncols), K)
}
#' @rdname DomainMatrix
#' @export
DomainMatrix_from_dict_sympy <- function(nrows, ncols, elemsdict, ...) {
  d <- .dm_as_dod(elemsdict)
  .sdm_check_bounds(d, nrows, ncols)
  items <- lapply(unlist(lapply(d, unname), recursive = FALSE), as_sym)
  r <- .dm_get_domain(items, ...)
  idx <- 0L
  for (i in names(d)) for (j in names(d[[i]])) { idx <- idx + 1L; d[[i]][[j]] <- r[[2L]][[idx]] }
  DomainMatrix(d, c(nrows, ncols), r[[1L]])
}
#' @rdname DomainMatrix
#' @export
DomainMatrix_from_Matrix <- function(M, fmt = "sparse", ...) {
  nr <- M$nrow; nc <- M$ncol
  rows <- lapply(seq_len(nr), function(i) lapply(seq_len(nc), function(j) .mget(M, i, j)))
  if (fmt == "dense") return(DomainMatrix_from_list_sympy(nr, nc, rows, ...))
  dod <- list()
  for (i in seq_len(nr)) {
    r <- list()
    for (j in seq_len(nc)) if (!is_zero_num(rows[[i]][[j]])) r[[.sk(j - 1L)]] <- rows[[i]][[j]]
    if (length(r)) dod[[.sk(i - 1L)]] <- r
  }
  if (!length(dod)) dod <- stats::setNames(list(), character())
  DomainMatrix_from_dict_sympy(nr, nc, dod, ...)
}
# construct_domain(items) with the elements as domain elements
.dm_get_domain <- function(items, ...) {
  K <- construct_domain(items, ...)[[1L]]
  if (isTRUE(list(...)$field) && .de_has_assoc_field(K)) K <- .de_get_field(K)
  list(K, lapply(items, function(e) .de_convert(K, e)))
}
#' @rdname DomainMatrix
#' @export
DomainMatrix_eye <- function(shape, domain) {
  shape <- as.integer(unlist(shape)); if (length(shape) == 1L) shape <- c(shape, shape)
  .dm_from_rep(.sdm_eye(shape, domain))
}
#' @rdname DomainMatrix
#' @export
DomainMatrix_diag <- function(diagonal, domain, shape = NULL) {
  if (is.null(shape)) shape <- c(length(diagonal), length(diagonal))
  .dm_from_rep(.sdm_diag(diagonal, domain, .dm_shape(shape)))
}
#' @rdname DomainMatrix
#' @export
DomainMatrix_zeros <- function(shape, domain, fmt = "sparse") .dm_from_rep(.sdm_new(list(), .dm_shape(shape), domain))
#' @rdname DomainMatrix
#' @export
DomainMatrix_ones <- function(shape, domain) {
  sh <- .dm_shape(shape)
  .dm_from_rep(.ddm_new(rep(list(rep(list(.de_one(domain)), sh[2L])), sh[1L]), sh, domain))
}
#' @rdname DomainMatrix
#' @export
DomainMatrix_from_list_flat <- function(elements, shape, domain) {
  sh <- .dm_shape(shape)
  if (length(elements) != sh[1L] * sh[2L]) stop(.dm_err("DMBadInputError", "Inconsistent flat-list shape"))
  .dm_from_rep(.ddm_new(lapply(seq_len(sh[1L]), function(i) elements[(i - 1L) * sh[2L] + seq_len(sh[2L])]), sh, domain))
}
#' @rdname DomainMatrix
#' @export
DomainMatrix_from_dod <- function(dod, shape, domain) {
  d <- .dm_as_dod(dod)
  d <- lapply(d, function(r) Filter(.de_nz, r)); d <- Filter(length, d)
  .dm_from_rep(.sdm_new(d, .dm_shape(shape), domain))
}
#' @rdname DomainMatrix
#' @export
DomainMatrix_from_dok <- function(dok, shape, domain) {
  d <- list()
  keys <- if (inherits(dok, "pydict")) dok$keys else lapply(names(dok), function(k) as.integer(strsplit(k, ",")[[1L]]))
  vals <- if (inherits(dok, "pydict")) dok$values else unname(dok)
  for (t in seq_along(keys)) {
    if (!.de_nz(vals[[t]])) next
    ij <- as.integer(unlist(keys[[t]]))
    i <- .sk(ij[1L]); if (is.null(d[[i]])) d[[i]] <- list()
    d[[i]][[.sk(ij[2L])]] <- vals[[t]]
  }
  .dm_from_rep(.sdm_new(d, .dm_shape(shape), domain))
}

.dm_shape_of <- function(A) .subset2(.subset2(A, "rep"), "shape")
.dm_dom <- function(A) .subset2(.subset2(A, "rep"), "domain")
.dm_fmt <- function(A) if (inherits(.subset2(A, "rep"), "DDM")) "dense" else "sparse"
.dm_to_sparse <- function(A) if (.dm_fmt(A) == "sparse") A else .dm_from_rep(.ddm_to_sdm(.subset2(A, "rep")))
.dm_to_dense <- function(A) if (.dm_fmt(A) == "dense") A else .dm_from_rep(.sdm_to_ddm(.subset2(A, "rep")))
.dm_convert_to <- function(A, K) {
  if (.de_same(K, .dm_dom(A))) return(A)
  .dm_from_rep(.subset2(A, "rep")$convert_to(K))
}
.de_unify <- function(K0, K1) {
  if (.de_same(K0, K1)) return(K0)
  k0 <- .de_kind(K0); k1 <- .de_kind(K1)
  if (k0 == "EXRAW") return(K0)
  if (k1 == "EXRAW") return(K1)
  if (k0 == "Algebraic" && k1 %in% c("ZZ", "QQ")) return(K0)
  if (k1 == "Algebraic" && k0 %in% c("ZZ", "QQ")) return(K1)
  if (k0 == "GF" || k1 == "GF") {
    if (k0 == "GF" && k1 %in% c("ZZ", "QQ")) return(K0)
    if (k1 == "GF" && k0 %in% c("ZZ", "QQ")) return(K1)
  }
  .dom_unify(K0, K1)
}
.dm_unify_domain <- function(mats) {
  doms <- lapply(mats, .dm_dom)
  K <- Reduce(.de_unify, doms)
  lapply(mats, function(M) .dm_convert_to(M, K))
}
.dm_unify <- function(mats, fmt = NULL) {
  mats <- .dm_unify_domain(mats)
  if (!is.null(fmt)) {
    fmts <- unique(vapply(mats, .dm_fmt, ""))
    if (length(fmts) > 1L) mats <- lapply(mats, if (fmt == "sparse") .dm_to_sparse else .dm_to_dense)
  }
  mats
}
.dm_check <- function(A, op, B, ashape, bshape) {
  if (!.de_same(.dm_dom(A), .dm_dom(B)))
    stop(.dm_err("DMDomainError", paste0("Domain mismatch: ", format(.dm_dom(A)), " ", op, " ", format(.dm_dom(B)))))
  if (!identical(as.integer(ashape), as.integer(bshape)))
    stop(.dm_err("DMShapeError", paste0("Shape mismatch: ", .shape_str(.dm_shape_of(A)), " ", op, " ", .shape_str(.dm_shape_of(B)))))
  if (.dm_fmt(A) != .dm_fmt(B))
    stop(.dm_err("DMFormatError", paste0("Format mismatch: ", .dm_fmt(A), " ", op, " ", .dm_fmt(B))))
}
.dm_to_field <- function(A) .dm_convert_to(A, .de_get_field(.dm_dom(A)))
.dm_has_assoc_field <- function(K) .de_has_assoc_field(K)

# rref.py
.dm_row_density <- function(A) {
  d <- .subset2(.subset2(A, "rep")$to_sdm(), "d"); ncols <- .dm_shape_of(A)[2L]
  if (!length(d)) return(list(0, 0L, ncols))
  list(sum(vapply(d, length, 0L)) / length(d), length(d), ncols)
}
.dm_flat_nz <- function(A) unname(unlist(lapply(unname(.subset2(.subset2(A, "rep")$to_sdm(), "d")), unname), recursive = FALSE))
.dm_bitlen <- function(z) { z <- abs(gmp::as.bigz(z)); if (z == 0) 0L else nchar(as.character(z, b = 2)) }
.dm_choose_ZZ <- function(A) {
  PARAM <- 10000
  rd <- .dm_row_density(A); density <- rd[[1L]]; nrows_nz <- rd[[2L]]; ncols <- rd[[3L]]
  if (nrows_nz < 10L) return(if (density < ncols / 2) "GJ" else "FF")
  if (density < 5) return("GJ")
  if (density > 5 + PARAM / nrows_nz) return("FF")
  bits <- max(vapply(.dm_flat_nz(A), .dm_bitlen, 0L), 1L)
  wideness <- max(1, 2 / 3 * ncols / nrows_nz)
  max_density <- (5 + PARAM / (nrows_nz * bits^2)) * wideness
  if (density < max_density) "GJ" else "FF"
}
.dm_choose_QQ <- function(A) {
  rd <- .dm_row_density(A); density <- rd[[1L]]; ncols <- rd[[3L]]
  if (density < min(5, ncols / 2)) return("GJ")
  el <- .dm_flat_nz(A)
  numer_bits <- max(vapply(el, function(e) .dm_bitlen(gmp::numerator(e)), 0L), 1L)
  denom_lcm <- gmp::as.bigz(1)
  for (e in el) {
    denom_lcm <- gmp::lcm.bigz(denom_lcm, gmp::denominator(e))
    if (.dm_bitlen(denom_lcm) > 5 * numer_bits) return("GJ")
  }
  if (.dm_bitlen(denom_lcm) < 50L) "CD" else "FF"
}
.dm_rref_choose <- function(A, method, denominator) {
  if (method != "auto") {
    if (endsWith(method, "_dense")) return(list(sub("_dense$", "", method), "dense"))
    return(list(method, "sparse"))
  }
  use_fmt <- "sparse"
  k <- .de_kind(.dm_dom(A))
  method <- if (k == "ZZ") .dm_choose_ZZ(A)
    else if (k == "QQ") .dm_choose_QQ(A)
    else if (k %in% c("RR", "CC")) { use_fmt <- "dense"; "GJ" }
    else if (k == "EX" && .dm_fmt(A) == "dense" && !denominator) { use_fmt <- "dense"; "GJ" }
    else if (denominator) "FF" else "GJ"
  list(method, use_fmt)
}
.dm_to_fmt <- function(A, fmt) if (fmt == "dense") .dm_to_dense(A) else .dm_to_sparse(A)
.dm_to_field_if <- function(A) if (.de_has_assoc_field(.dm_dom(A))) .dm_to_field(A) else A
.dm_rref_GJ <- function(A) {
  rep <- .subset2(A, "rep")
  if (inherits(rep, "SDM")) {
    r <- .sdm_irref(.subset2(rep, "d"))
    return(list(.dm_from_rep(.sdm_new(r$rref, .subset2(rep, "shape"), .subset2(rep, "domain"))), .pytuple(as.list(r$pivots))))
  }
  K <- .subset2(rep, "domain")
  r <- .ddm_irref(.subset2(rep, "rows"), .de_kind(K) %in% c("RR", "CC"))
  list(.dm_from_rep(.ddm_new(r$rows, .subset2(rep, "shape"), K)), .pytuple(as.list(r$pivots)))
}
.dm_rref_den_FF <- function(A) {
  rep <- .subset2(A, "rep"); K <- .subset2(rep, "domain")
  if (inherits(rep, "SDM")) {
    r <- .sdm_rref_den(.subset2(rep, "d"), K)
    return(list(.dm_from_rep(.sdm_new(r$rref, .subset2(rep, "shape"), K)), r$den, .pytuple(as.list(r$pivots))))
  }
  r <- .ddm_irref_den(.subset2(rep, "rows"), K)
  list(.dm_from_rep(.ddm_new(r$rows, .subset2(rep, "shape"), K)), r$den, .pytuple(as.list(r$pivots)))
}
.dm_rref <- function(A, method = "auto") {
  ch <- .dm_rref_choose(A, method, FALSE)
  old <- .dm_fmt(A); A <- .dm_to_fmt(A, ch[[2L]])
  if (ch[[1L]] == "GJ") { r <- .dm_rref_GJ(.dm_to_field_if(A)); M <- r[[1L]]; piv <- r[[2L]] }
  else if (ch[[1L]] == "FF") { r <- .dm_rref_den_FF(A); M <- .dm_div(.dm_to_field_if(r[[1L]]), r[[2L]], .dm_dom(r[[1L]])); piv <- r[[3L]] }
  else if (ch[[1L]] == "CD") {
    Mr <- .dm_clear_denoms_rowwise(A, TRUE)[[2L]]
    r <- .dm_rref_den_FF(Mr); M <- .dm_div(.dm_to_field_if(r[[1L]]), r[[2L]], .dm_dom(r[[1L]])); piv <- r[[3L]]
  } else stop(.py_err("ValueError", paste0("Unknown method for rref: ", ch[[1L]])))
  .pytuple(list(.dm_to_fmt(M, old), piv))
}
.dm_rref_den <- function(A, method = "auto", keep_domain = TRUE) {
  ch <- .dm_rref_choose(A, method, TRUE)
  old <- .dm_fmt(A); A <- .dm_to_fmt(A, ch[[2L]])
  if (ch[[1L]] == "FF") { r <- .dm_rref_den_FF(A); M <- r[[1L]]; den <- r[[2L]]; piv <- r[[3L]] }
  else if (ch[[1L]] == "GJ") {
    r <- .dm_rref_GJ(.dm_to_field_if(A)); Mf <- r[[1L]]; piv <- r[[2L]]
    if (keep_domain && !.de_same(.dm_dom(Mf), .dm_dom(A))) {
      M <- .dm_clear_denoms(Mf, TRUE)[[2L]]
      den <- if (length(piv)) M$rep$getitem(0L, piv[[1L]]) else .de_one(.dm_dom(M))
    } else { M <- Mf; den <- .de_one(.dm_dom(M)) }
  } else if (ch[[1L]] == "CD") {
    Mr <- .dm_clear_denoms_rowwise(A, TRUE)[[2L]]
    r <- .dm_rref_den_FF(Mr); Mrr <- r[[1L]]; den <- r[[2L]]; piv <- r[[3L]]
    if (keep_domain && !.de_same(.dm_dom(Mrr), .dm_dom(A))) {
      M <- .dm_div(.dm_to_field_if(Mrr), den, .dm_dom(Mrr)); den <- .de_one(.dm_dom(A))
    } else {
      M <- Mrr
      den <- if (length(piv)) M$rep$getitem(0L, piv[[1L]]) else .de_one(.dm_dom(M))
    }
  } else stop(.py_err("ValueError", paste0("Unknown method for rref: ", ch[[1L]])))
  .pytuple(list(.dm_to_fmt(M, old), den, piv))
}

# DomainMatrix.__truediv__ with an element of the domain K0
.dm_div <- function(A, lamda, K0) {
  K <- .dm_dom(A)
  lam <- if (inherits(lamda, "DomainScalar")) lamda else DomainScalar(.de_convert(K0, lamda), K0)
  A <- .dm_to_field(A)
  Kf <- .de_unify(.dm_dom(A), lam$domain)
  A <- .dm_convert_to(A, Kf); el <- .de_convert(Kf, lam$element, lam$domain)
  if (.de_is_zero(el)) stop(.py_err("ZeroDivisionError", ""))
  if (.de_is_one(Kf, el)) return(A)
  .dm_from_rep(.subset2(A, "rep")$mul(.de_recip(el)))
}

# dup_clear_denoms on a list of elements of K0 with ring K1
.de_clear_denoms <- function(elems, K0, K1, convert) {
  common <- .de_one(K1)
  for (c in elems) common <- .de_lcm(K1, common, .de_denom(K0, c))
  f <- if (.de_is_one(K1, common)) elems else lapply(elems, function(c) c * .de_convert(K0, common, K1))
  if (convert) f <- lapply(f, function(c) .de_convert(K1, c, K0))
  list(common, f)
}
.dm_clear_denoms <- function(A, convert = FALSE) {
  K0 <- .dm_dom(A)
  K1 <- if (.de_has_assoc_ring(K0)) .de_get_ring(K0) else K0
  rep <- .subset2(A, "rep"); s <- .subset2(rep, "shape")
  if (inherits(rep, "DDM")) {
    el <- unlist(.subset2(rep, "rows"), recursive = FALSE)
    nz <- vapply(el, .de_nz, NA)
    r <- .de_clear_denoms(el[nz], K0, K1, convert)
    Kn <- if (convert) K1 else K0
    out <- rep(list(.de_zero(Kn)), length(el)); out[nz] <- r[[2L]]
    M <- .dm_from_rep(.ddm_new(lapply(seq_len(s[1L]), function(i) out[(i - 1L) * s[2L] + seq_len(s[2L])]), s, Kn))
  } else {
    d <- .subset2(rep, "d")
    el <- unlist(lapply(unname(d), unname), recursive = FALSE)
    r <- .de_clear_denoms(el, K0, K1, convert)
    Kn <- if (convert) K1 else K0
    idx <- 0L
    for (i in names(d)) for (j in names(d[[i]])) { idx <- idx + 1L; d[[i]][[j]] <- r[[2L]][[idx]] }
    M <- .dm_from_rep(.sdm_new(d, s, Kn))
  }
  .pytuple(list(DomainScalar(r[[1L]], K1), M))
}
.dm_clear_denoms_rowwise <- function(A, convert = FALSE) {
  K0 <- .dm_dom(A)
  K1 <- if (.de_has_assoc_ring(K0)) .de_get_ring(K0) else K0
  d <- .subset2(.subset2(A, "rep")$to_sdm(), "d")
  m <- .dm_shape_of(A)[1L]
  diagonals <- rep(list(.de_one(K0)), m)
  dnum <- list()
  for (i in names(d)) {
    r <- .de_clear_denoms(unname(d[[i]]), K0, K1, convert)
    diagonals[[as.integer(i) + 1L]] <- r[[1L]]
    dnum[[i]] <- stats::setNames(r[[2L]], names(d[[i]]))
  }
  Kn <- if (convert) K1 else K0
  den <- DomainMatrix_diag(lapply(diagonals, function(z) .de_convert(K1, z)), K1)
  rep <- .subset2(A, "rep")
  num <- if (inherits(rep, "DDM")) .dm_from_rep(.sdm_to_ddm(.sdm_new(dnum, .subset2(rep, "shape"), Kn)))
         else .dm_from_rep(.sdm_new(dnum, .subset2(rep, "shape"), Kn))
  .pytuple(list(den, num))
}

# dup_content and dup_primitive over K
.de_content <- function(K, elems) {
  k <- .de_kind(K)
  if (!length(elems)) return(.de_zero(K))
  if (.de_is_field(K) && !(k %in% c("QQ", "FracField"))) return(.de_one(K))
  cont <- .de_zero(K)
  for (c in elems) { cont <- .de_gcd(K, cont, c); if (k == "ZZ" && cont == 1) break }
  cont
}
.dm_flat_elems <- function(A) {
  rep <- .subset2(A, "rep")
  if (inherits(rep, "SDM")) return(unlist(lapply(unname(.subset2(rep, "d")), unname), recursive = FALSE))
  Filter(.de_nz, unlist(.subset2(rep, "rows"), recursive = FALSE))
}
.dm_map_nz <- function(A, f, K = .dm_dom(A)) {
  rep <- .subset2(A, "rep")
  if (inherits(rep, "SDM")) {
    d <- lapply(.subset2(rep, "d"), function(r) lapply(r, f))
    return(.dm_from_rep(.sdm_new(d, .subset2(rep, "shape"), K)))
  }
  z <- .de_zero(K)
  .dm_from_rep(.ddm_new(lapply(.subset2(rep, "rows"), function(r) lapply(r, function(e) if (.de_nz(e)) f(e) else z)), .subset2(rep, "shape"), K))
}
.de_canonical_unit <- function(K, a) {
  if (inherits(a, "sym_gaussian")) {
    u <- .gi_canonical_unit(c(as.numeric(a$re), as.numeric(a$im)))
    return(structure(list(re = gmp::as.bigq(u[1L]), im = gmp::as.bigq(u[2L]), kind = a$kind), class = "sym_gaussian"))
  }
  if (.de_is_negative(K, a)) return(-.de_one(K))
  if (.de_kind(K) %in% c("QQ", "RR", "CC", "GF", "Algebraic", "EX", "EXRAW") && .de_nz(a) && .de_kind(K) != "QQ") return(.de_one(K))
  .de_one(K)
}
.dm_cancel_denom <- function(A, denom) {
  K <- .dm_dom(A)
  if (.de_is_zero(denom)) stop(.py_err("ZeroDivisionError", "denominator is zero"))
  if (.de_is_one(K, denom)) return(.pytuple(list(A, denom)))
  el <- .dm_flat_elems(A)
  u <- if (.de_is_negative(K, denom)) -.de_one(K) else .de_canonical_unit(K, denom)
  content <- .de_content(K, el)
  common <- .de_gcd(K, content, denom)
  if (!.de_is_one(K, content)) {
    common <- .de_gcd(K, content, denom)
    if (!.de_is_one(K, common)) {
      A <- .dm_map_nz(A, function(e) .de_exquo(K, e, common))
      denom <- .de_quo(K, denom, common)
    }
  }
  if (!.de_is_one(K, u)) {
    A <- .dm_map_nz(A, function(e) e * u)
    denom <- u * denom
  } else if (.de_is_one(K, common)) return(.pytuple(list(A, denom)))
  .pytuple(list(A, denom))
}
.dm_cancel_denom_elementwise <- function(A, denom) {
  K <- .dm_dom(A)
  if (.de_is_zero(denom)) stop(.py_err("ZeroDivisionError", "denominator is zero"))
  if (.de_is_one(K, denom)) return(.pytuple(list(A, DomainMatrix_ones(.dm_shape_of(A), K))))
  .pytuple(list(.dm_map_nz(A, function(e) .de_cofactors(K, e, denom)[[2L]]),
                .dm_map_nz(A, function(e) .de_cofactors(K, e, denom)[[3L]])))
}

# dup helpers on lists of elements, leading coefficient first
.dl_strip <- function(f, K) { while (length(f) && .de_is_zero(f[[1L]])) f <- f[-1L]; f }
.dl_mul <- function(f, g, K) {
  if (!length(f) || !length(g)) return(list())
  n <- length(f) + length(g) - 1L
  out <- rep(list(.de_zero(K)), n)
  for (i in seq_along(f)) for (j in seq_along(g)) out[[i + j - 1L]] <- out[[i + j - 1L]] + f[[i]] * g[[j]]
  .dl_strip(out, K)
}

.dm_charpoly_base <- function(M) {
  K <- .dm_dom(M)
  s <- .dm_shape_of(M)
  density <- if (s[1L]) length(.dm_flat_elems(M)) / s[1L]^2 else 0
  M <- if (density < 0.5) .dm_to_sparse(M) else .dm_to_dense(M)
  cd <- .de_is_field(K) && .de_has_assoc_ring(K)
  if (cd) {
    r <- .dm_clear_denoms(M, TRUE); d <- r[[1L]]$element; M <- r[[2L]]; Kr <- .dm_dom(M)
  }
  cp <- .subset2(M, "rep")$charpoly()
  if (cd) {
    cp <- lapply(cp, function(c) .de_convert(K, c, Kr))
    # dup_transform(cp, [1, 0], [1/d]): coefficient i scaled by (1/d)^i
    dinv <- .de_recip(.de_convert(K, d, Kr))
    sc <- .de_one(K)
    for (i in seq_along(cp)) { cp[[i]] <- cp[[i]] * sc; sc <- sc * dinv }
  }
  cp
}
.dm_collect_factors <- function(fl) {
  keys <- character(); facs <- list(); mult <- integer()
  for (p in fl) {
    k <- paste(vapply(p[[1L]], .de_repr, ""), collapse = "|")
    w <- match(k, keys)
    if (is.na(w)) { keys <- c(keys, k); facs[[length(facs) + 1L]] <- p[[1L]]; mult <- c(mult, p[[2L]]) }
    else mult[w] <- mult[w] + p[[2L]]
  }
  o <- order(vapply(facs, length, 0L), method = "radix")
  lapply(o, function(i) .pytuple(list(facs[[i]], mult[i])))
}
.dm_charpoly_factor_blocks <- function(M) {
  s <- .dm_shape_of(M)
  if (s[1L] != s[2L]) stop(.dm_err("DMNonSquareMatrixError", "not square"))
  comps <- .subset2(M, "rep")$scc()
  bf <- lapply(comps, function(ind) list(.dm_charpoly_base(.dm_from_rep(.subset2(M, "rep")$extract(as.list(ind), as.list(ind)))), 1L))
  .dm_collect_factors_keep(bf)
}
# _collect_factors: Counter of factors, then sorted by (len(f), e, f)
.dm_collect_factors_keep <- function(fl) {
  keys <- character(); facs <- list(); mult <- integer()
  for (p in fl) {
    k <- paste(vapply(p[[1L]], .de_repr, ""), collapse = "|")
    w <- match(k, keys)
    if (is.na(w)) { keys <- c(keys, k); facs[[length(facs) + 1L]] <- p[[1L]]; mult <- c(mult, p[[2L]]) }
    else mult[w] <- mult[w] + p[[2L]]
  }
  if (!length(facs)) return(list())
  # _sort_factors(multiple=True): key (len(f), e, f)
  o <- .dm_sort_factor_order(facs, mult)
  lapply(o, function(i) .pytuple(list(facs[[i]], mult[i])))
}
.dm_sort_factor_order <- function(facs, mult) {
  n <- length(facs)
  idx <- seq_len(n)
  less <- function(a, b) {
    la <- length(facs[[a]]); lb <- length(facs[[b]])
    if (la != lb) return(la < lb)
    if (mult[a] != mult[b]) return(mult[a] < mult[b])
    for (t in seq_len(la)) {
      x <- facs[[a]][[t]]; y <- facs[[b]][[t]]
      if (!.de_eq(x, y)) return(isTRUE(.de_lt(x, y)))
    }
    FALSE
  }
  for (i in seq_len(n)) for (j in seq_len(n - i)) if (less(idx[j + 1L], idx[j])) { t <- idx[j]; idx[j] <- idx[j + 1L]; idx[j + 1L] <- t }
  idx
}
.de_lt <- function(x, y) {
  if (inherits(x, c("bigz", "bigq"))) return(x < y)
  if (inherits(x, "sym_modint")) return(x$val < y$val)
  if (is.numeric(x)) return(x < y)
  .de_repr(x) < .de_repr(y)
}
.dm_charpoly <- function(M) {
  K <- .dm_dom(M)
  cp <- list(.de_one(K))
  for (fm in .dm_charpoly_factor_blocks(M)) for (t in seq_len(fm[[2L]])) cp <- .dl_mul(cp, fm[[1L]], K)
  cp
}

.dm_nullspace <- function(A, divide_last = FALSE) {
  K <- .dm_dom(A)
  if (divide_last && !.de_is_field(K)) stop(.dm_err("DMNotAField", "Cannot normalize vectors over a non-field"))
  if (divide_last) { r <- .dm_rref(A); Ar <- r[[1L]]; piv <- r[[2L]] }
  else {
    r <- .dm_rref_den(A); Ar <- r[[1L]]; den <- r[[2L]]; piv <- r[[3L]]
    u <- .de_canonical_unit(K, den)
    if (!.de_is_one(K, u)) Ar <- .dm_from_rep(.subset2(Ar, "rep")$mul(u))
  }
  res <- .subset2(Ar, "rep")$nullspace_from_rref(as.list(unlist(piv)))
  .dm_from_rep(res[[1L]])
}
.dm_inv <- function(A) {
  K <- .dm_dom(A)
  if (!.de_is_field(K)) stop(.dm_err("DMNotAField", "Not a field"))
  s <- .dm_shape_of(A)
  if (s[1L] != s[2L]) stop(.dm_err("DMNonSquareMatrixError"))
  .dm_from_rep(.subset2(A, "rep")$inv())
}
.dm_det <- function(A) {
  s <- .dm_shape_of(A)
  if (s[1L] != s[2L]) stop(.dm_err("DMNonSquareMatrixError"))
  .subset2(A, "rep")$det()
}
.dm_matmul <- function(A, B) {
  .dm_check(A, "*", B, .dm_shape_of(A)[2L], .dm_shape_of(B)[1L])
  .dm_from_rep(.subset2(A, "rep")$matmul(.subset2(B, "rep")))
}
.dm_scalarmul <- function(A, lamda, reverse) {
  K <- .dm_dom(A)
  if (.de_is_zero(lamda) && .de_kind(K) != "EXRAW") return(DomainMatrix_zeros(.dm_shape_of(A), K))
  if (.de_is_one(K, lamda)) return(A)
  .dm_from_rep(if (reverse) .subset2(A, "rep")$rmul(lamda) else .subset2(A, "rep")$mul(lamda))
}
.dmx_pow <- function(A, n) {
  s <- .dm_shape_of(A)
  if (s[1L] != s[2L]) stop(.dm_err("DMNonSquareMatrixError", "Power of a nonsquare matrix"))
  n <- as.integer(as.character(n))
  if (n < 0L) stop(.py_err("NotImplementedError", "Negative powers"))
  if (n == 0L) return(DomainMatrix_eye(s[1L], .dm_dom(A)))
  if (n == 1L) return(A)
  if (n %% 2L == 1L) return(.dmx_mul(A, .dmx_pow(A, n - 1L)))
  r <- .dmx_pow(A, n %/% 2L)
  .dmx_mul(r, r)
}
.dmx_mul <- function(A, B) {
  if (inherits(B, "DomainMatrix")) {
    ab <- .dm_unify(list(A, B), "dense")
    return(.dm_matmul(ab[[1L]], ab[[2L]]))
  }
  if (inherits(B, "DomainScalar")) {
    K <- .de_unify(.dm_dom(A), B$domain)
    return(.dm_scalarmul(.dm_convert_to(A, K), .de_convert(K, B$element, B$domain), FALSE))
  }
  .dm_scalarmul(A, .de_convert(.dm_dom(A), B), FALSE)
}
.dm_rmul <- function(B, A) {
  if (inherits(B, "DomainScalar")) {
    K <- .de_unify(.dm_dom(A), B$domain)
    return(.dm_scalarmul(.dm_convert_to(A, K), .de_convert(K, B$element, B$domain), TRUE))
  }
  .dm_scalarmul(A, .de_convert(.dm_dom(A), B), TRUE)
}
.dm_addsub <- function(A, B, op) {
  if (!inherits(B, "DomainMatrix") || !inherits(A, "DomainMatrix"))
    stop(.py_err("TypeError", paste0("unsupported operand type(s) for ", op)))
  ab <- .dm_unify(list(A, B), "dense"); A <- ab[[1L]]; B <- ab[[2L]]
  .dm_check(A, op, B, .dm_shape_of(A), .dm_shape_of(B))
  .dm_from_rep(if (op == "+") .subset2(A, "rep")$add(.subset2(B, "rep")) else .subset2(A, "rep")$sub(.subset2(B, "rep")))
}
.dm_eq <- function(A, B) {
  inherits(A, "DomainMatrix") && inherits(B, "DomainMatrix") && .de_same(.dm_dom(A), .dm_dom(B)) &&
    .dm_rep_eq(.subset2(A, "rep"), .subset2(B, "rep"))
}
.dm_hstack <- function(A, ...) {
  ms <- .dm_unify(c(list(A), list(...)), .dm_fmt(A))
  .dm_from_rep(do.call(.subset2(ms[[1L]], "rep")$hstack, lapply(ms[-1L], function(M) .subset2(M, "rep"))))
}
.dm_vstack <- function(A, ...) {
  ms <- .dm_unify(c(list(A), list(...)), "dense")
  .dm_from_rep(do.call(.subset2(ms[[1L]], "rep")$vstack, lapply(ms[-1L], function(M) .subset2(M, "rep"))))
}
.dm_extract <- function(A, rws, cls) .dm_from_rep(.subset2(A, "rep")$extract(as.list(as.integer(unlist(rws))), as.list(as.integer(unlist(cls)))))
# A[i, j] with 0-based index vectors (slices already expanded)
.dm_extract_idx <- function(A, ri, ci) {
  rep <- .subset2(A, "rep")
  if (inherits(rep, "SDM")) return(.dm_from_rep(.sdm_extract_slice(rep, ri, ci)))
  rows <- .subset2(rep, "rows")
  .dm_from_rep(.ddm_new(lapply(ri + 1L, function(i) rows[[i]][ci + 1L]), c(length(ri), length(ci)), .subset2(rep, "domain")))
}
.dm_solve_den_rref <- function(A, b) {
  s <- .dm_shape_of(A); bs <- .dm_shape_of(b); m <- s[1L]; n <- s[2L]
  if (m != bs[1L]) stop(.dm_err("DMShapeError", "Matrix equation shape mismatch."))
  if (m < n) stop(.dm_err("DMShapeError", "Underdetermined matrix equation."))
  Aaug <- .dm_hstack(A, b)
  r <- .dm_rref_den(Aaug); piv <- unlist(r[[3L]])
  if (length(piv) != n || (length(piv) && piv[length(piv)] >= n)) stop(.dm_err("DMNonInvertibleMatrixError", "Non-unique solution."))
  .pytuple(list(.dm_extract_idx(r[[1L]], seq_len(n) - 1L, seq_len(bs[2L]) + n - 1L), r[[2L]]))
}
.dm_adj_poly_det <- function(A, cp = NULL) {
  s <- .dm_shape_of(A)
  if (s[1L] != s[2L]) stop(.dm_err("DMNonSquareMatrixError", "Matrix must be square"))
  if (is.null(cp)) cp <- .dm_charpoly(A)
  n <- length(cp)
  if (n %% 2L) { detA <- cp[[n]]; f <- lapply(cp[-n], function(c) -c) }
  else { detA <- -cp[[n]]; f <- cp[-n] }
  list(f, detA)
}
.dm_eval_poly_mul <- function(A, p, B) {
  s <- .dm_shape_of(A); bs <- .dm_shape_of(B)
  if (s[1L] != s[2L]) stop(.dm_err("DMNonSquareMatrixError", "Matrix must be square"))
  if (bs[1L] != s[2L]) stop(.dm_err("DMShapeError", "Matrices are not aligned"))
  if (!.de_same(.dm_dom(A), .dm_dom(B))) stop(.dm_err("DMDomainError", "Matrices must have the same domain"))
  if (!length(p)) return(DomainMatrix_zeros(bs, .dm_dom(B)))
  pAB <- .dm_rmul(p[[1L]], B)
  for (pi in p[-1L]) pAB <- .dm_addsub(.dmx_mul(A, pAB), .dm_rmul(pi, B), "+")
  pAB
}
.dm_eval_poly <- function(A, p) {
  s <- .dm_shape_of(A); K <- .dm_dom(A)
  if (s[1L] != s[2L]) stop(.dm_err("DMNonSquareMatrixError", "Matrix must be square"))
  if (!length(p)) return(DomainMatrix_zeros(s, K))
  I <- DomainMatrix_eye(s, K)
  if (length(p) == 1L) return(.dm_rmul(p[[1L]], I))
  pA <- .dm_rmul(p[[1L]], I)
  for (pi in p[-1L]) pA <- .dm_addsub(.dmx_mul(A, pA), .dm_rmul(pi, I), "+")
  pA
}
.dm_solve_den_charpoly <- function(A, b, cp = NULL, check = TRUE) {
  ab <- .dm_unify(list(A, b)); A <- ab[[1L]]; b <- ab[[2L]]
  s <- .dm_shape_of(A); bs <- .dm_shape_of(b)
  if (s[1L] != s[2L]) stop(.dm_err("DMNonSquareMatrixError", "Matrix must be square"))
  if (bs[1L] != s[1L]) stop(.dm_err("DMShapeError", "Matrix and vector must have the same number of rows"))
  fd <- .dm_adj_poly_det(A, cp)
  if (check && .de_is_zero(fd[[2L]])) stop(.dm_err("DMNonInvertibleMatrixError", "Matrix is not invertible"))
  .pytuple(list(.dm_eval_poly_mul(A, fd[[1L]], b), fd[[2L]]))
}
.dm_adj_det <- function(A) {
  m <- .dm_shape_of(A)[1L]
  r <- .dm_solve_den_charpoly(A, DomainMatrix_eye(c(m, m), .dm_dom(A)), check = FALSE)
  adjA <- r[[1L]]
  if (.dm_fmt(A) == "dense") adjA <- .dm_to_dense(adjA)
  .pytuple(list(adjA, r[[2L]]))
}
.dm_solve_den <- function(A, b, method = NULL) {
  s <- .dm_shape_of(A); bs <- .dm_shape_of(b)
  if (s[1L] != bs[1L]) stop(.dm_err("DMShapeError", "Matrix equation shape mismatch."))
  if (is.null(method)) method <- "rref"
  else if (method == "charpoly" && s[1L] != s[2L]) stop(.dm_err("DMNonSquareMatrixError", "method='charpoly' requires a square matrix."))
  if (method == "charpoly") .dm_solve_den_charpoly(A, b)
  else if (method == "rref") .dm_solve_den_rref(A, b)
  else stop(.dm_err("DMBadInputError", "method should be 'rref' or 'charpoly'"))
}
.dm_lu <- function(A) {
  if (!.de_is_field(.dm_dom(A))) stop(.dm_err("DMNotAField", "Not a field"))
  r <- .subset2(A, "rep")$lu()
  .pytuple(list(.dm_from_rep(r[[1L]]), .dm_from_rep(r[[2L]]), r[[3L]]))
}
.dm_lu_solve <- function(A, rhs) {
  if (.dm_shape_of(A)[1L] != .dm_shape_of(rhs)[1L]) stop(.dm_err("DMShapeError", "Shape"))
  if (!.de_is_field(.dm_dom(A))) stop(.dm_err("DMNotAField", "Not a field"))
  .dm_from_rep(.subset2(A, "rep")$lu_solve(.subset2(rhs, "rep")))
}
.dm_ground_eigenvals <- function(A) {
  K <- .dm_dom(A)
  p <- .dm_charpoly(A)
  fl <- .dl_factor_list(p, K)
  keys <- list(); vals <- list()
  for (bf in fl) {
    base <- bf[[1L]]
    if (length(base) == 2L) { keys[[length(keys) + 1L]] <- -base[[2L]] / base[[1L]]; vals[[length(vals) + 1L]] <- bf[[2L]] }
  }
  .pydict(keys, vals)
}
# dup_factor_list over ZZ or QQ through Poly factor_list
.dl_factor_list <- function(p, K) {
  X <- Dummy("x")
  e <- .add(lapply(seq_along(p), function(i) .mul(list(.de_to_sympy(K, p[[i]]), .pow(X, .int(length(p) - i))))))
  fl <- factor_list(e, X)
  lapply(fl[[2L]]$args, function(t) {
    cl <- .coeff_list(expand(t$args[[1L]]), X)
    list(lapply(cl, function(z) .de_convert(K, z)), as.integer(t$args[[2L]]$value))
  })
}

.dm_method <- function(A, name) {
  rep <- .subset2(A, "rep"); K <- .subset2(rep, "domain"); s <- .subset2(rep, "shape")
  switch(name,
    rep = rep, shape = .pytuple(as.list(s)), domain = K,
    is_zero_matrix = rep$is_zero_matrix(), is_upper = rep$is_upper(), is_lower = rep$is_lower(),
    is_diagonal = rep$is_diagonal(), is_square = s[1L] == s[2L],
    getitem = function(i, j) DomainScalar(rep$getitem(i, j), K),
    getitem_sympy = function(i, j) .de_to_sympy(K, rep$getitem(i, j)),
    setitem = function(i, j, value) .dm_from_rep(rep$setitem(i, j, value)),
    extract = function(rowslist, colslist) .dm_extract(A, rowslist, colslist),
    copy = function() A,
    convert_to = function(K2) .dm_convert_to(A, K2),
    to_sympy = function() .dm_convert_to(A, EXRAW),
    to_field = function() .dm_to_field(A),
    to_sparse = function() .dm_to_sparse(A), to_dense = function() .dm_to_dense(A),
    to_ddm = function() rep$to_ddm(), to_sdm = function() rep$to_sdm(),
    to_dfm = function() stop(.py_err("NotImplementedError", "DFM requires python-flint")),
    to_dfm_or_ddm = function() rep$to_ddm(),
    unify = function(..., fmt = NULL) .pytuple(.dm_unify(c(list(A), list(...)), fmt)),
    unify_eq = function(B) {
      if (!identical(s, .dm_shape_of(B))) return(FALSE)
      ab <- .dm_unify(list(A, B)); .dm_eq(ab[[1L]], ab[[2L]])
    },
    to_Matrix = function() {
      rows <- rep$to_list()
      if (!length(rows) || !length(rows[[1L]])) return(zeros(s[1L], s[2L]))
      Matrix(lapply(rows, function(r) lapply(r, function(e) .de_to_sympy(K, e))))
    },
    to_list = function() rep$to_list(), to_list_flat = function() rep$to_list_flat(),
    to_flat_nz = function() rep$to_flat_nz(),
    to_dod = function() rep$to_dod(), to_dok = function() rep$to_dok(),
    iter_values = function() rep$iter_values(),
    nnz = function() rep$nnz(),
    transpose = function() .dm_from_rep(rep$transpose()),
    flat = function() rep$to_list_flat(),
    diagonal = function() rep$diagonal(),
    rank = function() length(.dm_rref(A)[[2L]]),
    hstack = function(...) .dm_hstack(A, ...),
    vstack = function(...) .dm_vstack(A, ...),
    applyfunc = function(func, domain = NULL) .dm_from_rep(rep$applyfunc(func, if (is.null(domain)) K else domain)),
    add = function(B) { .dm_check(A, "+", B, s, .dm_shape_of(B)); .dm_from_rep(rep$add(.subset2(B, "rep"))) },
    sub = function(B) { .dm_check(A, "-", B, s, .dm_shape_of(B)); .dm_from_rep(rep$sub(.subset2(B, "rep"))) },
    neg = function() .dm_from_rep(rep$neg()),
    mul = function(b) .dm_from_rep(rep$mul(b)), rmul = function(b) .dm_from_rep(rep$rmul(b)),
    matmul = function(B) .dm_matmul(A, B),
    scalarmul = function(lamda) .dm_scalarmul(A, lamda, FALSE),
    rscalarmul = function(lamda) .dm_scalarmul(A, lamda, TRUE),
    mul_elementwise = function(B) .dm_from_rep(rep$mul_elementwise(.subset2(B, "rep"))),
    pow = function(n) .dmx_pow(A, n),
    scc = function() {
      if (s[1L] != s[2L]) stop(.dm_err("DMNonSquareMatrixError", "Matrix must be square for scc"))
      lapply(rep$scc(), as.list)
    },
    clear_denoms = function(convert = FALSE) .dm_clear_denoms(A, convert),
    clear_denoms_rowwise = function(convert = FALSE) .dm_clear_denoms_rowwise(A, convert),
    cancel_denom = function(denom) .dm_cancel_denom(A, denom),
    cancel_denom_elementwise = function(denom) .dm_cancel_denom_elementwise(A, denom),
    content = function() .de_content(K, .dm_flat_elems(A)),
    primitive = function() {
      cont <- .de_content(K, .dm_flat_elems(A))
      prim <- if (.de_is_zero(cont) || .de_is_one(K, cont)) A else .dm_map_nz(A, function(e) .de_exquo(K, e, cont))
      .pytuple(list(cont, prim))
    },
    rref = function(method = "auto") .dm_rref(A, method),
    rref_den = function(method = "auto", keep_domain = TRUE) .dm_rref_den(A, method, keep_domain),
    columnspace = function() {
      if (!.de_is_field(K)) stop(.dm_err("DMNotAField", "Not a field"))
      piv <- unlist(.dm_rref(A)[[2L]])
      .dm_extract(A, seq_len(s[1L]) - 1L, piv)
    },
    rowspace = function() {
      if (!.de_is_field(K)) stop(.dm_err("DMNotAField", "Not a field"))
      piv <- unlist(.dm_rref(A)[[2L]])
      .dm_extract(A, seq_along(piv) - 1L, seq_len(s[2L]) - 1L)
    },
    nullspace = function(divide_last = FALSE) .dm_nullspace(A, divide_last),
    nullspace_from_rref = function(pivots = NULL) .dm_from_rep(rep$nullspace_from_rref(pivots)[[1L]]),
    inv = function() .dm_inv(A),
    det = function() .dm_det(A),
    adj_det = function() .dm_adj_det(A),
    adjugate = function() .dm_adj_det(A)[[1L]],
    inv_den = function(method = NULL) .dm_solve_den(A, DomainMatrix_eye(s, K), method),
    solve_den = function(b, method = NULL) .dm_solve_den(A, b, method),
    solve_den_rref = function(b) .dm_solve_den_rref(A, b),
    solve_den_charpoly = function(b, cp = NULL, check = TRUE) .dm_solve_den_charpoly(A, b, cp, check),
    adj_poly_det = function(cp = NULL) .pytuple(.dm_adj_poly_det(A, cp)),
    eval_poly = function(p) .dm_eval_poly(A, p),
    eval_poly_mul = function(p, B) .dm_eval_poly_mul(A, p, B),
    lu = function() .dm_lu(A),
    qr = function() { r <- rep$qr(); .pytuple(list(.dm_from_rep(r[[1L]]), .dm_from_rep(r[[2L]]))) },
    lu_solve = function(rhs) .dm_lu_solve(A, rhs),
    fflu = function() .pytuple(lapply(rep$fflu(), .dm_from_rep)),
    charpoly = function() .dm_charpoly(A),
    charpoly_base = function() .dm_charpoly_base(A),
    charpoly_berk = function() rep$charpoly(),
    charpoly_factor_blocks = function() .dm_charpoly_factor_blocks(A),
    charpoly_factor_list = function() {
      out <- list()
      for (fm in .dm_charpoly_factor_blocks(A)) for (fj in .dl_factor_list(fm[[1L]], K)) out[[length(out) + 1L]] <- list(fj[[1L]], fm[[2L]] * fj[[2L]])
      .dm_collect_factors_keep(out)
    },
    ground_eigenvals = function() .dm_ground_eigenvals(A),
    lll = function(delta = gmp::as.bigq(3, 4)) .dm_from_rep(rep$lll(delta)),
    lll_transform = function(delta = gmp::as.bigq(3, 4)) { r <- rep$lll_transform(delta); .pytuple(list(.dm_from_rep(r[[1L]]), .dm_from_rep(r[[2L]]))) },
    eye = function(shape, domain) DomainMatrix_eye(shape, domain),
    zeros = function(shape, domain, fmt = "sparse") DomainMatrix_zeros(shape, domain),
    ones = function(shape, domain) DomainMatrix_ones(shape, domain),
    diag = function(diagonal, domain, shape = NULL) DomainMatrix_diag(diagonal, domain, shape),
    from_rep = function(rep2) DomainMatrix_from_rep(rep2),
    from_list_flat = function(elements, shape, domain) DomainMatrix_from_list_flat(elements, shape, domain),
    from_flat_nz = function(elements, data, domain) {
      keys <- data[[1L]]; sh <- data[[2L]]
      DomainMatrix_from_dok(.pydict(as.list(keys), elements), sh, domain)
    },
    from_dod_like = function(dod, domain = NULL) {
      M <- DomainMatrix_from_dod(dod, s, if (is.null(domain)) K else domain)
      if (.dm_fmt(A) == "dense") .dm_to_dense(M) else M
    },
    choose_domain = function(...) {
      fz <- rep$to_sdm()
      d <- .subset2(fz, "d")
      items <- lapply(unlist(lapply(unname(d), unname), recursive = FALSE), function(e) .de_to_sympy(K, e))
      r <- .dm_get_domain(items, ...)
      idx <- 0L
      for (i in names(d)) for (j in names(d[[i]])) { idx <- idx + 1L; d[[i]][[j]] <- r[[2L]][[idx]] }
      M <- .dm_from_rep(.sdm_new(d, s, r[[1L]]))
      if (.dm_fmt(A) == "dense") .dm_to_dense(M) else M
    },
    stop(.py_err("AttributeError", paste0("'DomainMatrix' object has no attribute '", name, "'"))))
}

#' @rdname DomainMatrix
#' @export
`$.DomainMatrix` <- function(x, name) .dm_method(x, name)
#' @rdname DomainMatrix
#' @export
format.DomainMatrix <- function(x, ...) .dm_dm_repr(x)
#' @rdname DomainMatrix
#' @export
print.DomainMatrix <- function(x, ...) { cat(format(x), "\n", sep = ""); invisible(x) }
#' @rdname DomainMatrix
#' @export
`[.DomainMatrix` <- function(x, i, j, drop = FALSE) {
  s <- .dm_shape_of(x)
  ri <- if (missing(i)) seq_len(s[1L]) - 1L else as.integer(i) - 1L
  ci <- if (missing(j)) seq_len(s[2L]) - 1L else as.integer(j) - 1L
  if (!missing(i) && !missing(j) && length(ri) == 1L && length(ci) == 1L) return(.dm_method(x, "getitem")(ri, ci))
  if (any(ri < 0L | ri >= s[1L])) stop(.py_err("IndexError", "Row index out of range"))
  if (any(ci < 0L | ci >= s[2L])) stop(.py_err("IndexError", "Column index out of range"))
  .dm_extract_idx(x, ri, ci)
}
#' @rdname DomainMatrix
#' @export
Ops.DomainMatrix <- function(e1, e2) {
  if (missing(e2)) {
    if (.Generic == "-") return(.dm_from_rep(.subset2(e1, "rep")$neg()))
    if (.Generic == "+") return(e1)
  }
  switch(.Generic,
    `+` = .dm_addsub(e1, e2, "+"),
    `-` = .dm_addsub(e1, e2, "-"),
    `*` = if (inherits(e1, "DomainMatrix")) .dmx_mul(e1, e2) else .dm_rmul(e1, e2),
    `^` = .dmx_pow(e1, e2),
    `/` = .dm_div(e1, if (inherits(e2, "DomainScalar")) e2 else if (.de_is_field(.dm_dom(e1)) && !(is.numeric(e2) && e2 == round(e2))) DomainScalar(.de_convert(.dm_dom(e1), e2), .dm_dom(e1)) else DomainScalar(.de_convert(ZZ, e2), ZZ), ZZ),
    `==` = .dm_eq(e1, e2),
    `!=` = !.dm_eq(e1, e2),
    stop("operator ", .Generic, " is not defined for DomainMatrix", call. = FALSE))
}
#' @rdname DomainMatrix
#' @export
chooseOpsMethod.DomainMatrix <- function(x, y, mx, my, cl, reverse) TRUE

# ---------------------------------------------------------------------------
# DomainScalar

#' Elements of a domain with their domain
#'
#' `DomainScalar(element, domain)` is SymPy's
#' `sympy.polys.matrices.domainscalar.DomainScalar`, as returned by indexing
#' a [DomainMatrix()]. Arithmetic unifies the domains; `$to_sympy()`,
#' `$to_domain(K)`, `$convert_to(K)`, `$unify(other)`, `$is_zero()` and
#' `$is_one()` are its methods, and `DomainScalar_from_sympy(expr)` builds one
#' from an expression.
#'
#' @param element An element of `domain`.
#' @param domain A domain.
#' @param expr An expression.
#' @param x A `DomainScalar`.
#' @param name A method name.
#' @param e1,e2 Scalars or integers.
#' @param ... Unused.
#' @return A `DomainScalar`, or the value of the method.
#' @examples
#' a <- DomainScalar(ZZ(3), ZZ)
#' b <- DomainScalar(QQ(1, 2), QQ)
#' a + b
#' a$to_sympy()
#' @export
DomainScalar <- function(element, domain) {
  if (!inherits(domain, "sym_domain")) stop(.py_err("TypeError", "domain should be of type Domain"))
  structure(list(element = element, domain = domain), class = "DomainScalar")
}
#' @rdname DomainScalar
#' @export
DomainScalar_from_sympy <- function(expr) {
  r <- .dm_get_domain(list(as_sym(expr)))
  DomainScalar(r[[2L]][[1L]], r[[1L]])
}
.ds_to_domain <- function(a, K) DomainScalar(.de_convert(K, .subset2(a, "element"), .subset2(a, "domain")), K)
.ds_unify <- function(a, b) {
  K <- .de_unify(.subset2(a, "domain"), .subset2(b, "domain"))
  list(.ds_to_domain(a, K), .ds_to_domain(b, K))
}
#' @rdname DomainScalar
#' @export
`$.DomainScalar` <- function(x, name) {
  el <- .subset2(x, "element"); K <- .subset2(x, "domain")
  switch(name,
    element = el, domain = K,
    to_sympy = function() .de_to_sympy(K, el),
    to_domain = function(domain) .ds_to_domain(x, domain),
    convert_to = function(domain) .ds_to_domain(x, domain),
    unify = function(other) .pytuple(.ds_unify(x, other)),
    is_zero = function() .de_is_zero(el),
    is_one = function() .de_is_one(K, el),
    stop(.py_err("AttributeError", paste0("'DomainScalar' object has no attribute '", name, "'"))))
}
#' @rdname DomainScalar
#' @export
format.DomainScalar <- function(x, ...) gsub("**", "^", .dmpy_repr(.subset2(x, "element")), fixed = TRUE)
#' @rdname DomainScalar
#' @export
print.DomainScalar <- function(x, ...) { cat(format(x), "\n", sep = ""); invisible(x) }
#' @rdname DomainScalar
#' @export
Ops.DomainScalar <- function(e1, e2) {
  if (missing(e2)) {
    if (.Generic == "-") return(DomainScalar(-e1$element, e1$domain))
    return(e1)
  }
  if (.Generic == "^") return(DomainScalar(e1$element^as.integer(as.character(e2)), e1$domain))
  if (.Generic == "*" && !inherits(e2, "DomainScalar") && is.numeric(e2)) e2 <- DomainScalar(.de_convert(ZZ, e2), ZZ)
  if (.Generic == "*" && !inherits(e1, "DomainScalar") && is.numeric(e1)) e1 <- DomainScalar(.de_convert(ZZ, e1), ZZ)
  if (!inherits(e1, "DomainScalar") || !inherits(e2, "DomainScalar")) {
    if (.Generic == "==") return(FALSE)
    if (.Generic == "!=") return(TRUE)
    stop(.py_err("TypeError", paste0("unsupported operand type(s) for ", .Generic)))
  }
  u <- .ds_unify(e1, e2); a <- u[[1L]]$element; b <- u[[2L]]$element; K <- u[[1L]]$domain
  switch(.Generic,
    `+` = DomainScalar(a + b, K), `-` = DomainScalar(a - b, K), `*` = DomainScalar(a * b, K),
    `%/%` = DomainScalar(.de_quo(K, a, b), K), `%%` = DomainScalar(.de_rem(K, a, b), K),
    `==` = .de_eq(e1$element, e2$element) && .de_same(e1$domain, e2$domain),
    `!=` = !(.de_eq(e1$element, e2$element) && .de_same(e1$domain, e2$domain)),
    stop("operator ", .Generic, " is not defined for DomainScalar", call. = FALSE))
}
#' @rdname DomainScalar
#' @export
chooseOpsMethod.DomainScalar <- function(x, y, mx, my, cl, reverse) TRUE

# ---------------------------------------------------------------------------
# normalforms.py over a principal ideal domain

.de_is_pid <- function(K) .de_kind(K) %in% c("ZZ", "QQ", "RR", "CC", "QQ_I", "ZZ_I", "GF", "Algebraic", "FracField", "EX", "EXRAW") ||
  (.de_kind(K) == "PolyRing" && length(.dk_gens(K)) == 1L && .de_is_field(.dk_ground(K)))
# K.div(a, b): list(q, r)
.de_div <- function(K, a, b) {
  k <- .de_kind(K)
  if (k == "ZZ") return(list(gmp::divq.bigz(a, b), gmp::mod.bigz(a, b)))
  if (k == "PolyRing") {
    pr <- qp_div_multi(a$p, list(b$p))
    return(list(.pe(pr$q[[1L]], K), .pe(pr$r, K)))
  }
  list(a / b, .de_zero(K))
}
# K.gcdex(a, b): list(s, t, h) with s*a + t*b = h
.de_gcdex <- function(K, a, b) {
  k <- .de_kind(K)
  if (k == "ZZ") return(.igcdex(a, b))
  if (k == "PolyRing") {
    X <- .dk_gens(K)[[1L]]
    r <- gcdex(.de_pe_sym(a), .de_pe_sym(b), X)
    return(lapply(r$args, function(z) .pe_from_sym(K, z)))
  }
  if (.de_is_zero(a) && .de_is_zero(b)) return(list(.de_one(K), .de_zero(K), .de_zero(K)))
  if (.de_nz(a)) list(.de_recip(a), .de_zero(K), .de_one(K)) else list(.de_zero(K), .de_recip(b), .de_one(K))
}
.de_canon_unit_pid <- function(K, a) {
  k <- .de_kind(K)
  if (k == "ZZ") return(if (a < 0) gmp::as.bigz(-1) else gmp::as.bigz(1))
  if (k == "PolyRing") {
    lc <- qp_lc(a$p)
    return(.pe_const(K, if (.dk(.dk_ground(K)) == "ZZ") (if (lc < 0) -1 else 1) else 1 / lc))
  }
  if (.de_is_field(K)) return(.de_recip(a))
  .de_one(K)
}
.nf_snd <- function(m, K, rows, cols, full) {
  if (!.de_is_pid(K)) stop(.py_err("ValueError", paste0("The matrix entries must be over a principal ideal domain, but got ", format(K))))
  zero <- .de_zero(K); one <- .de_one(K)
  eye <- function(n) lapply(seq_len(n), function(j) lapply(seq_len(n), function(i) if (i == j) one else zero))
  if (rows == 0L || cols == 0L) return(if (full) list(list(), eye(rows), eye(cols)) else list())
  s <- if (full) eye(rows); t <- if (full) eye(cols)
  cv <- function(v) .de_convert(K, v)
  add_rows <- function(M, i, j, a, b, c, d) {
    a <- cv(a); b <- cv(b); c <- cv(c); d <- cv(d)
    for (k in seq_along(M[[1L]])) { e <- M[[i]][[k]]; M[[i]][[k]] <- a * e + b * M[[j]][[k]]; M[[j]][[k]] <- c * e + d * M[[j]][[k]] }
    M
  }
  add_cols <- function(M, i, j, a, b, c, d) {
    a <- cv(a); b <- cv(b); c <- cv(c); d <- cv(d)
    for (k in seq_along(M)) { e <- M[[k]][[i]]; M[[k]][[i]] <- a * e + b * M[[k]][[j]]; M[[k]][[j]] <- c * e + d * M[[k]][[j]] }
    M
  }
  clear_column <- function() {
    pivot <- m[[1L]][[1L]]
    if (rows >= 2L) for (j in 2:rows) {
      if (.de_is_zero(m[[j]][[1L]])) next
      dr <- .de_div(K, m[[j]][[1L]], pivot)
      if (.de_is_zero(dr[[2L]])) {
        m <<- add_rows(m, 1L, j, 1, 0, -dr[[1L]], 1); if (full) s <<- add_rows(s, 1L, j, 1, 0, -dr[[1L]], 1)
      } else {
        g <- .de_gcdex(K, pivot, m[[j]][[1L]])
        d0 <- .de_exquo(K, m[[j]][[1L]], g[[3L]]); dj <- .de_exquo(K, pivot, g[[3L]])
        m <<- add_rows(m, 1L, j, g[[1L]], g[[2L]], d0, -dj); if (full) s <<- add_rows(s, 1L, j, g[[1L]], g[[2L]], d0, -dj)
        pivot <- g[[3L]]
      }
    }
  }
  clear_row <- function() {
    pivot <- m[[1L]][[1L]]
    if (cols >= 2L) for (j in 2:cols) {
      if (.de_is_zero(m[[1L]][[j]])) next
      dr <- .de_div(K, m[[1L]][[j]], pivot)
      if (.de_is_zero(dr[[2L]])) {
        m <<- add_cols(m, 1L, j, 1, 0, -dr[[1L]], 1); if (full) t <<- add_cols(t, 1L, j, 1, 0, -dr[[1L]], 1)
      } else {
        g <- .de_gcdex(K, pivot, m[[1L]][[j]])
        d0 <- .de_exquo(K, m[[1L]][[j]], g[[3L]]); dj <- .de_exquo(K, pivot, g[[3L]])
        m <<- add_cols(m, 1L, j, g[[1L]], g[[2L]], d0, -dj); if (full) t <<- add_cols(t, 1L, j, g[[1L]], g[[2L]], d0, -dj)
        pivot <- g[[3L]]
      }
    }
  }
  ind <- which(vapply(m, function(r) .de_nz(r[[1L]]), NA))
  if (length(ind) && ind[1L] != 1L) {
    tt <- m[[1L]]; m[[1L]] <- m[[ind[1L]]]; m[[ind[1L]]] <- tt
    if (full) { tt <- s[[1L]]; s[[1L]] <- s[[ind[1L]]]; s[[ind[1L]]] <- tt }
  } else {
    ind <- which(vapply(m[[1L]], .de_nz, NA))
    if (length(ind) && ind[1L] != 1L) {
      sw <- function(r) { tt <- r[[1L]]; r[[1L]] <- r[[ind[1L]]]; r[[ind[1L]]] <- tt; r }
      m <- lapply(m, sw); if (full) t <- lapply(t, sw)
    }
  }
  more <- function() (cols >= 2L && any(vapply(m[[1L]][-1L], .de_nz, NA))) || (rows >= 2L && any(vapply(m[-1L], function(r) .de_nz(r[[1L]]), NA)))
  while (more()) { clear_column(); clear_row() }
  if (.de_nz(m[[1L]][[1L]])) {
    c <- .de_canon_unit_pid(K, m[[1L]][[1L]])
    if (!.de_is_one(K, c)) {
      m[[1L]][[1L]] <- m[[1L]][[1L]] * c
      if (full) s[[1L]] <- lapply(s[[1L]], function(e) e * c)
    }
  }
  if (rows == 1L || cols == 1L) invs <- list()
  else {
    ret <- .nf_snd(lapply(m[-1L], function(r) r[-1L]), K, rows - 1L, cols - 1L, full)
    if (full) {
      invs <- ret[[1L]]
      s2 <- c(list(c(list(one), rep(list(zero), rows - 1L))), lapply(ret[[2L]], function(r) c(list(zero), r)))
      t2 <- c(list(c(list(one), rep(list(zero), cols - 1L))), lapply(ret[[3L]], function(r) c(list(zero), r)))
      s <- .ddm_imatmul(.ddm_zeros(c(rows, rows), K)[["rows"]], s2, s)
      t <- .ddm_imatmul(.ddm_zeros(c(cols, cols), K)[["rows"]], t, t2)
    } else invs <- ret
  }
  if (.de_nz(m[[1L]][[1L]])) {
    result <- c(list(m[[1L]][[1L]]), invs)
    for (i in seq_len(length(result) - 1L)) {
      a <- result[[i]]; b <- result[[i + 1L]]
      if (.de_nz(b) && .de_nz(.de_div(K, b, a)[[2L]])) {
        if (full) { g <- .de_gcdex(K, a, b); x <- g[[1L]]; y <- g[[2L]]; d <- g[[3L]] } else d <- .de_gcd(K, a, b)
        alpha <- .de_div(K, a, d)[[1L]]
        if (full) {
          beta <- .de_div(K, b, d)[[1L]]
          s <- add_rows(s, i, i + 1L, 1, 0, x, 1)
          t <- add_cols(t, i, i + 1L, 1, y, 0, 1)
          s <- add_rows(s, i, i + 1L, 1, -alpha, 0, 1)
          t <- add_cols(t, i, i + 1L, 1, 0, -beta, 1)
          s <- add_rows(s, i, i + 1L, 0, 1, -1, 0)
        }
        result[[i + 1L]] <- b * alpha
        result[[i]] <- d
      } else break
    }
  } else {
    if (full) {
      if (rows > 1L) s <- c(s[-1L], s[1L])
      if (cols > 1L) t <- lapply(t, function(r) c(r[-1L], r[1L]))
    }
    result <- c(invs, list(m[[1L]][[1L]]))
  }
  if (full) list(result, s, t) else result
}
.dm_invariant_factors <- function(A) {
  .pytuple(.nf_snd(.subset2(A, "rep")$to_list(), .dm_dom(A), .dm_shape_of(A)[1L], .dm_shape_of(A)[2L], FALSE))
}
.dm_smith_normal_form <- function(A) {
  invs <- .nf_snd(.subset2(A, "rep")$to_list(), .dm_dom(A), .dm_shape_of(A)[1L], .dm_shape_of(A)[2L], FALSE)
  DomainMatrix_diag(invs, .dm_dom(A), .dm_shape_of(A))
}

#' Smith normal form decomposition and test
#'
#' `smith_normal_decomp(m)` returns the Smith normal form `smf` of a
#' [DomainMatrix()] over a principal ideal domain with invertible matrices
#' `s` and `t` such that `smf == s * m * t`; `is_smith_normal_form(m)` tests
#' whether a matrix is in Smith normal form. [smith_normal_form()],
#' [invariant_factors()] and [hermite_normal_form()] also accept a
#' `DomainMatrix`.
#'
#' @param m A `DomainMatrix`.
#' @return A tuple of three matrices, or `TRUE` or `FALSE`.
#' @examples
#' m <- DM(list(list(12, 6, 4), list(3, 9, 6), list(2, 16, 14)), ZZ)
#' smith_normal_decomp(m)
#' is_smith_normal_form(DM(list(list(1, 0), list(0, 2)), ZZ))
#' @export
smith_normal_decomp <- function(m) {
  K <- .dm_dom(m); sh <- .dm_shape_of(m)
  r <- .nf_snd(.subset2(m, "rep")$to_list(), K, sh[1L], sh[2L], TRUE)
  smf <- .dm_to_dense(DomainMatrix_diag(r[[1L]], K, sh))
  .pytuple(list(smf, DomainMatrix(r[[2L]], c(sh[1L], sh[1L]), K), DomainMatrix(r[[3L]], c(sh[2L], sh[2L]), K)))
}
#' @rdname smith_normal_decomp
#' @export
is_smith_normal_form <- function(m) {
  K <- .dm_dom(m); sh <- .dm_shape_of(m); L <- .subset2(m, "rep")$to_list()
  for (i in seq_len(sh[1L])) for (j in seq_len(sh[2L])) if (i != j && .de_nz(L[[i]][[j]])) return(FALSE)
  up <- min(sh)
  if (up >= 2L) for (i in 2:up) {
    if (.de_is_zero(L[[i - 1L]][[i - 1L]])) { if (.de_nz(L[[i]][[i]])) return(FALSE) }
    else if (.de_nz(.de_div(K, L[[i]][[i]], L[[i - 1L]][[i - 1L]])[[2L]])) return(FALSE)
  }
  TRUE
}
.dm_hnf <- function(A, D = NULL, check_rank = FALSE) {
  if (.de_kind(.dm_dom(A)) != "ZZ") stop(.dm_err("DMDomainError", "Matrix must be over domain ZZ."))
  sh <- .dm_shape_of(A)
  rows <- lapply(.subset2(A, "rep")$to_list(), function(r) do.call(c, r))
  if (!is.null(D) && (!check_rank || length(.dm_rref(.dm_convert_to(A, QQ))[[2L]]) == sh[1L])) {
    if (inherits(D, "bigq") || (is_sym(D) && D$type != "Integer")) stop(.dm_err("DMDomainError", "Modulus D must be positive element of domain ZZ."))
    if (sh[2L] < sh[1L]) stop(.dm_err("DMShapeError", "Matrix must have at least as many columns as rows."))
    W <- .hnf_mod_D(rows, sh[2L], D)
    return(DomainMatrix(lapply(W, function(r) lapply(seq_along(r), function(k) r[k])), c(sh[1L], sh[1L]), ZZ))
  }
  h <- .hnf(rows, sh[2L])
  nc <- if (length(h)) length(h[[1L]]) else 0L
  DomainMatrix(lapply(h, function(r) lapply(seq_len(nc), function(k) r[k])), c(sh[1L], nc), ZZ)
}
