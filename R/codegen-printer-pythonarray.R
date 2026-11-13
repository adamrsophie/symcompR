# Array printing for the Python code printers, after the ArrayPrinter class
# of sympy/printing/pycode.py, and the CuPy and JAX printers of
# sympy/printing/numpy.py. ArrayPrinter is a mixin that comes first in the
# bases of NumPyPrinter, TensorflowPrinter and TorchPrinter; it prints array
# expressions as einsum, transpose and fold calls of the array module named by
# the `_module` attribute of the printer.

# StrPrinter._print_ArraySymbol and _print_ArrayElement
.cgp_set("StrPrinter",
         `_print_List` = function(p, e) .cgp_call(p, "_print_list", e$args),
         `_print_ArraySymbol` = function(p, e) .cgpp(p, .prt_array_name(e)),
         `_print_ArrayElement` = function(p, e)
           sprintf("%s[%s]", .cgp_paren(p, e$args[[1L]], .CG_PRECEDENCE[["Func"]], strict = TRUE),
                   paste(vapply(e$args[[2L]]$args, function(i) .cgpp(p, i), ""), collapse = ", ")))
# expr.name: a Symbol for ArraySymbol, the label string for IndexedBase
.prt_array_name <- function(e) if (e$type == "IndexedBase") sstr(e$args[[1L]]) else e$args[[1L]]

# AbstractPythonCodePrinter._expand_fold_binary_op and _expand_reduce_binary_op,
# _print_ITE and _print_Sum
.prt_fold <- function(p, op, args) {
  n <- length(args)
  if (n == 1L) return(.cgpp(p, args[[1L]]))
  sprintf("%s(%s, %s)", .cgp_module_format(p, op), .prt_fold(p, op, args[-n]), .cgpp(p, args[[n]]))
}
.prt_reduce <- function(p, op, args) {
  n <- length(args)
  if (n == 1L) return(.cgpp(p, args[[1L]]))
  h <- n %/% 2L
  sprintf("%s(%s, %s)", .cgp_module_format(p, op), .prt_reduce(p, op, args[seq_len(h)]),
          .prt_reduce(p, op, args[(h + 1L):n]))
}
.cgp_set("AbstractPythonCodePrinter",
         `_expand_fold_binary_op` = .prt_fold,
         `_expand_reduce_binary_op` = .prt_reduce,
         `_print_ITE` = function(p, e) .cgpp(p, .cg_ite_piecewise(e)),
         `_print_Sum` = function(p, e) {
           lims <- rev(e$args[-1L])
           loops <- vapply(lims, function(l) sprintf("for %s in range(%s, %s+1)", .cgpp(p, l$args[[1L]]),
                                                     .cgpp(p, l$args[[2L]]), .cgpp(p, l$args[[3L]])), "")
           sprintf("(builtins.sum(%s %s))", .cgpp(p, e$args[[1L]]), paste(loops, collapse = " "))
         })

# The letters of an einsum string: a to z, then A to Z
.prt_letters <- function() {
  pool <- c(letters, LETTERS)
  k <- 0L
  function() {
    k <<- k + 1L
    if (k > length(pool)) .cg_err("ValueError", "out of letters")
    pool[k]
  }
}

# ArrayPrinter._get_einsum_string: the einsum subscripts of operands with
# ranks `subranks` under the index groups `cind` (0-based), with the free and
# the repeated letters
.prt_einsum_string <- function(subranks, cind) {
  nxt <- .prt_letters()
  dkeys <- integer(); dvals <- integer()
  for (g in cind) for (j in g) {
    if (j %in% dkeys) dvals[match(j, dkeys)] <- min(g) else { dkeys <- c(dkeys, j); dvals <- c(dvals, min(g)) }
  }
  counter <- 0L
  indices <- list()
  for (r in subranks) {
    li <- integer()
    for (k in seq_len(r)) {
      li <- c(li, if (counter %in% dkeys) dvals[match(counter, dkeys)] else counter)
      counter <- counter + 1L
    }
    indices[[length(indices) + 1L]] <- li
  }
  mkeys <- integer(); mvals <- character()
  s <- ""; free <- character(); dum <- character()
  for (li in indices) {
    for (j in li) {
      if (!(j %in% mkeys)) { l <- nxt(); mkeys <- c(mkeys, j); mvals <- c(mvals, l) } else l <- mvals[match(j, mkeys)]
      s <- paste0(s, l)
      if (j %in% dkeys) { if (!(l %in% dum)) dum <- c(dum, l) } else free <- c(free, l)
    }
    s <- paste0(s, ",")
  }
  list(string = substr(s, 1L, nchar(s) - 1L), free = free, dum = dum)
}

.prt_amod <- function(p, attr) .cgp_module_format(p, paste0(.cgp_a(p, "_module"), ".", .cgp_a(p, attr)))
.prt_join <- function(p, args, sep) paste(vapply(args, function(a) .cgpp(p, a), ""), collapse = sep)
.prt_is_tp <- function(e) is_sym(e) && e$type == "ArrayTensorProduct"

# _ArrayExpr.shape of an operand, as a list
.prt_shape <- function(e) .ae_shape(e)

# ArrayPrinter._arrayify: convert_indexed_to_array, or the expression itself
# when that fails
.prt_arrayify <- function(e) tryCatch(convert_indexed_to_array(e), error = function(err) e)

# _is_atomic_code: an identifier or a single call
.prt_is_atomic_code <- function(code) {
  if (grepl("^[A-Za-z_][A-Za-z0-9_]*$", code)) return(TRUE)
  k <- regexpr("(", code, fixed = TRUE)
  if (k < 0L) return(FALSE)
  head <- substr(code, 1L, k - 1L); rest <- substring(code, k + 1L)
  if (!grepl("^[A-Za-z_][A-Za-z0-9_]*$", gsub(".", "", head, fixed = TRUE)) || !endsWith(rest, ")")) return(FALSE)
  ch <- strsplit(rest, "", fixed = TRUE)[[1L]]
  depth <- 1L
  for (i in seq_along(ch)) {
    depth <- depth + (ch[i] == "(") - (ch[i] == ")")
    if (depth == 0L) return(i == length(ch))
  }
  FALSE
}

# _ElementwiseOperand: a Dummy named after the printed operand, printed as
# its name, with the precedence of the operand
.prt_operand <- function(code, prec) {
  d <- unclass(Dummy(code))
  d$extra$prt_operand <- TRUE
  d$extra$prt_prec <- prec
  structure(d, class = "sym")
}

# Lambda.variables[0]
.prt_lambda_var <- function(lam) { v <- lam$args[[1L]]; if (v$type == "Tuple") v$args[[1L]] else v }

# ArrayPrinter._elementwise_body
.prt_elementwise_body <- function(p, e) {
  operand <- e$args[[2L]]
  if (is_sym(operand) && operand$type == "ArrayElementwiseApplyFunc") operand <- .prt_elementwise_body(p, operand)
  code <- .cgpp(p, operand)
  prec <- if (.prt_is_atomic_code(code)) .CG_PRECEDENCE[["Atom"]] else .cgp_prec(operand)
  ph <- .prt_operand(code, prec)
  lam <- e$args[[1L]]
  xreplace(lam$args[[2L]], list(list(.prt_lambda_var(lam), ph)))
}

.cgp_def("ArrayPrinter",
         methods = list(
           `_arrayify` = function(p, e) .prt_arrayify(e),
           `_get_einsum_string` = function(p, subranks, cind) .prt_einsum_string(subranks, cind),
           `_print_ArrayTensorProduct` = function(p, e) {
             nxt <- .prt_letters()
             cs <- paste(vapply(.ae_sub_ndim_list(e), function(r) paste(vapply(seq_len(r), function(k) nxt(), ""), collapse = ""), ""),
                         collapse = ",")
             sprintf("%s(\"%s\", %s)", .prt_amod(p, "_einsum"), cs, .prt_join(p, e$args, ", "))
           },
           `_print_ArrayContraction` = function(p, e) {
             base <- e$args[[1L]]
             cind <- .ae_contr_ind(e)
             if (.prt_is_tp(base)) {
               elems <- .prt_join(p, base$args, ",")
               ndim <- .ae_sub_ndim_list(base)
             } else {
               elems <- .cgpp(p, base)
               ndim <- length(.prt_shape(base))
             }
             es <- .prt_einsum_string(ndim, cind)
             if (!length(cind)) return(.cgpp(p, base))
             elems <- if (.prt_is_tp(base)) .prt_join(p, base$args, ",") else .cgpp(p, base)
             sprintf("%s(\"%s->%s\", %s)", .prt_amod(p, "_einsum"), es$string,
                     paste(sort(es$free, method = "radix"), collapse = ""), elems)
           },
           `_print_ArrayDiagonal` = function(p, e) {
             dind <- .ae_diag_ind(e)
             if (.prt_is_tp(e$args[[1L]])) {
               subranks <- .ae_sub_ndim_list(e$args[[1L]])
               elems <- e$args[[1L]]$args
             } else {
               subranks <- .ae_sub_ndim_list(e)
               elems <- list(e$args[[1L]])
             }
             es <- .prt_einsum_string(subranks, dind)
             sprintf("%s(\"%s->%s\", %s)", .prt_amod(p, "_einsum"), es$string,
                     paste(c(es$free, es$dum), collapse = ""), .prt_join(p, elems, ", "))
           },
           `_print_PermuteDims` = function(p, e)
             sprintf("%s(%s, %s)", .prt_amod(p, "_transpose"), .cgpp(p, e$args[[1L]]), .cgpp(p, as.list(.ae_paf(e$args[[2L]])))),
           `_print_ArrayAdd` = function(p, e)
             .prt_fold(p, paste0(.cgp_a(p, "_module"), ".", .cgp_a(p, "_add")), e$args),
           `_print_ArrayElementwiseApplyFunc` = function(p, e) .cgpp(p, .prt_elementwise_body(p, e)),
           `_print_Dummy` = function(p, e) {
             if (isTRUE(e$extra$prt_operand)) return(e$name)
             .cgp_super(p, "ArrayPrinter", "_print_Dummy", e)
           },
           `_print_OneArray` = function(p, e)
             sprintf("%s((%s,))", .prt_amod(p, "_ones"), .prt_join(p, e$args, ",")),
           `_print_ZeroArray` = function(p, e)
             sprintf("%s((%s,))", .prt_amod(p, "_zeros"), .prt_join(p, e$args, ",")),
           `_print_Assignment` = function(p, e)
             sprintf("%s = %s", .cgpp(p, .prt_arrayify(e$args[[1L]])), .cgpp(p, .prt_arrayify(e$args[[2L]]))),
           `_print_IndexedBase` = function(p, e) .cgp_call(p, "_print_ArraySymbol", e)
         ))

# NumPyPrinter methods of numpy.py that the base port lacks, and the array
# function names
.cgp_classes$NumPyPrinter$parents <- c("ArrayPrinter", "PythonCodePrinter")
.cgp_classes$NumPyPrinter$attrs[c("_add", "_einsum", "_transpose", "_ones", "_zeros")] <-
  list("add", "einsum", "transpose", "ones", "zeros")
rm(list = ls(.cgp_mro_cache), envir = .cgp_mro_cache)

.prt_matrix_tolist <- function(e) {
  sh <- .cgp_mshape(e)
  lapply(seq_len(sh[1L]) - 1L, function(i) lapply(seq_len(sh[2L]) - 1L, function(j) .cgp_mget(e, i, j)))
}

# NDimArray.tolist(): nested lists of the entries
.prt_ndarray_tolist <- function(data, shape) {
  if (length(shape) == 1L) return(data)
  n <- shape[1L]; step <- length(data) %/% n
  lapply(seq_len(n), function(k) .prt_ndarray_tolist(data[(k - 1L) * step + seq_len(step)], shape[-1L]))
}
.prt_ndarray_parts <- function(e) {
  if (inherits(e, "SymArray")) return(list(data = e$data, shape = as.integer(e$shape)))
  list(data = e$args, shape = as.integer(unlist(e$extra$shape)))
}

.cgp_set("NumPyPrinter",
         `_print_Contains` = function(p, e) {
           item <- e$args[[1L]]; s <- e$args[[2L]]
           if (identical(s$key, Integers$key)) return(sprintf("equal(mod(%s, 1), 0)", .cgpp(p, item)))
           .cg_err("NotImplementedError", sprintf("NumPy printing for Contains(%s, %s) not implemented", sstr(item), sstr(s)))
         },
         `_print_DotProduct` = function(p, e) {
           a1 <- e$args[[1L]]; a2 <- e$args[[2L]]
           if (!identical(.mshape(a1)[[1L]]$key, S_One$key)) a1 <- Transpose(a1)
           if (!identical(.mshape(a2)[[2L]]$key, S_One$key)) a2 <- Transpose(a2)
           sprintf("%s(%s, %s)", .cgp_np(p, "dot"), .cgpp(p, a1), .cgpp(p, a2))
         },
         `_print_Identity` = function(p, e) {
           sh <- .mshape(e)
           if (all(vapply(sh, function(d) d$type == "Integer", NA)))
             return(sprintf("%s(%s)", .cgp_np(p, "eye"), .cgpp(p, sh[[1L]])))
           .cg_err("NotImplementedError", "Symbolic matrix dimensions are not yet supported for identity matrices")
         },
         `_print_BlockMatrix` = function(p, e)
{
           # the blocks are a Tuple of row Tuples here
           rows <- lapply(e$args[[1L]]$args, function(r) r$args)
           sprintf("%s(%s)", .cgp_np(p, "block"), .cgpp(p, rows))
         },
         `_print_NDimArray` = function(p, e) {
           pr <- .prt_ndarray_parts(e)
           if (!length(pr$shape))
             return(sprintf("%s(%s)", .cgp_np(p, "array"), .cgpp(p, pr$data[[1L]])))
           if (0L %in% pr$shape)
             return(sprintf("%s(%s)", .prt_amod(p, "_zeros"), .cgpp(p, structure(as.list(pr$shape), class = "pytuple"))))
           sprintf("%s(%s)", .cgp_np(p, "array"), .cgpp(p, .prt_ndarray_tolist(pr$data, pr$shape)))
         },
         `_print_isinf` = .cgp_print_known_func,
         `_print_isnan` = .cgp_print_known_func)

# SciPyPrinter methods of numpy.py that the base port lacks
.prt_sp <- function(p, name) .cgp_module_format(p, paste0("scipy.special.", name))
.prt_sp_index <- function(fun, k) function(p, e) sprintf("%s(%s)[%d]", .prt_sp(p, fun), .cgpp(p, e$args[[1L]]), k)
.cgp_set("SciPyPrinter",
         `_print_assoc_legendre` = function(p, e)
           sprintf("%s(%s, %s, %s)", .prt_sp(p, "lpmv"), .cgpp(p, e$args[[2L]]), .cgpp(p, e$args[[1L]]), .cgpp(p, e$args[[3L]])),
         `_print_lowergamma` = function(p, e)
           sprintf("%s(%s)*%s(%s, %s)", .prt_sp(p, "gamma"), .cgpp(p, e$args[[1L]]), .prt_sp(p, "gammainc"),
                   .cgpp(p, e$args[[1L]]), .cgpp(p, e$args[[2L]])),
         `_print_uppergamma` = function(p, e)
           sprintf("%s(%s)*%s(%s, %s)", .prt_sp(p, "gamma"), .cgpp(p, e$args[[1L]]), .prt_sp(p, "gammaincc"),
                   .cgpp(p, e$args[[1L]]), .cgpp(p, e$args[[2L]])),
         `_print_betainc_regularized` = function(p, e) {
           f <- .prt_sp(p, "betainc")
           a <- vapply(e$args, function(z) .cgpp(p, z), "")
           sprintf("%s(%s, %s, %s) - %s(%s, %s, %s)", f, a[1L], a[2L], a[4L], f, a[1L], a[2L], a[3L])
         },
         `_print_fresnels` = .prt_sp_index("fresnel", 0L),
         `_print_fresnelc` = .prt_sp_index("fresnel", 1L),
         `_print_owens_t` = function(p, e) sprintf("%s(%s, %s)", .prt_sp(p, "owens_t"), .cgpp(p, e$args[[1L]]), .cgpp(p, e$args[[2L]])),
         `_print_airyai` = .prt_sp_index("airy", 0L),
         `_print_airyaiprime` = .prt_sp_index("airy", 1L),
         `_print_airybi` = .prt_sp_index("airy", 2L),
         `_print_airybiprime` = .prt_sp_index("airy", 3L),
         `_print_Si` = .prt_sp_index("sici", 0L),
         `_print_Ci` = .prt_sp_index("sici", 1L))

# CuPy and JAX
.cgp_def("CuPyPrinter", "NumPyPrinter",
         attrs = list(`_module` = "cupy", `_kf` = .cg_prefixed(.cg_kf_numpy, "cupy."), `_kc` = .cg_prefixed(.cg_kc_numpy, "cupy.")))
.cgp_add_known("CuPyPrinter", names(.cg_kf_numpy), names(.cg_kc_numpy))
.cgp_def("JaxPrinter", "NumPyPrinter",
         attrs = list(`_module` = "jax.numpy", `_kf` = .cg_prefixed(.cg_kf_numpy, "jax.numpy."),
                      `_kc` = .cg_prefixed(.cg_kc_numpy, "jax.numpy.")),
         init = function(p, settings, nxt) {
           nxt(settings)
           p$printmethod <- "_jaxcode"
         },
         methods = list(
           `_print_And` = function(p, e)
             sprintf("%s(%s.asarray([%s]), axis=0)", .cgp_np(p, "all"), .cgp_module_format(p, .cgp_a(p, "_module")), .prt_join(p, e$args, ",")),
           `_print_Or` = function(p, e)
             sprintf("%s(%s.asarray([%s]), axis=0)", .cgp_np(p, "any"), .cgp_module_format(p, .cgp_a(p, "_module")), .prt_join(p, e$args, ","))
         ))
.cgp_add_known("JaxPrinter", names(.cg_kf_numpy), names(.cg_kc_numpy))

#' Array code printers
#'
#' Printer objects after SymPy's `CuPyPrinter` and `JaxPrinter`
#' (`sympy.printing.numpy`), which print Python code for the CuPy and JAX
#' array libraries as `NumPyPrinter()` does for NumPy. The NumPy, SciPy,
#' CuPy, JAX, TensorFlow and PyTorch printers print array expressions
#' (`ArrayTensorProduct()`, `ArrayContraction()`, `ArrayDiagonal()`,
#' `PermuteDims()`, `ArrayAdd()`, `ZeroArray()`, `OneArray()`) as einsum,
#' transpose and add calls.
#'
#' @param settings A named list of printer settings.
#' @return A printer, for [printer_doprint()].
#' @examples
#' x <- Symbol("x")
#' printer_doprint(CuPyPrinter(), sin(x) + sqrt(x))
#' printer_doprint(JaxPrinter(), And(x > 1, x < 2))
#' M <- MatrixSymbol("M", 2, 2)
#' N <- MatrixSymbol("N", 2, 2)
#' printer_doprint(NumPyPrinter(), ArrayContraction(ArrayTensorProduct(M, N), c(1, 2)))
#' @name array-code-printers
NULL

#' @rdname array-code-printers
#' @export
CuPyPrinter <- function(settings = NULL) .cgp_new("CuPyPrinter", settings)
#' @rdname array-code-printers
#' @export
JaxPrinter <- function(settings = NULL) .cgp_new("JaxPrinter", settings)
