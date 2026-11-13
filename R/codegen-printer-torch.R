# PyTorch code printer, after sympy/printing/pytorch.py. The mapping and the
# shared helpers are in R/codegen-printer-tensorflow.R.

.prt_torch_mapping <- c(
  Abs = "torch.abs", sign = "torch.sign", ceiling = "torch.ceil", floor = "torch.floor", log = "torch.log",
  exp = "torch.exp", Sqrt = "torch.sqrt", cos = "torch.cos", acos = "torch.acos", sin = "torch.sin",
  asin = "torch.asin", tan = "torch.tan", atan = "torch.atan", atan2 = "torch.atan2", cosh = "torch.cosh",
  acosh = "torch.acosh", sinh = "torch.sinh", asinh = "torch.asinh", tanh = "torch.tanh", atanh = "torch.atanh",
  Pow = "torch.pow", re = "torch.real", im = "torch.imag", arg = "torch.angle", erf = "torch.erf",
  loggamma = "torch.lgamma", Equality = "torch.eq", Unequality = "torch.ne", StrictGreaterThan = "torch.gt",
  StrictLessThan = "torch.lt", LessThan = "torch.le", GreaterThan = "torch.ge", And = "torch.logical_and",
  Or = "torch.logical_or", Not = "torch.logical_not", Max = "torch.maximum", Min = "torch.minimum",
  MatAdd = "torch.add", HadamardProduct = "torch.mul", Trace = "torch.trace", Determinant = "torch.det")

.prt_grad <- function(res, var) sprintf("torch.autograd.grad(%s, %s, create_graph=True)[0]", res, var)

.prt_torch_fun <- .prt_print_mapped("TorchPrinter")
.cgp_def("TorchPrinter", c("ArrayPrinter", "AbstractPythonCodePrinter"),
         attrs = list(printmethod = "_torchcode", module = "sympy.printing.pytorch", mapping = .prt_torch_mapping,
                      `_default_settings` = c(.cgp_classes$AbstractPythonCodePrinter$attrs[["_default_settings"]],
                                              list(torch_version = NULL, requires_grad = FALSE, dtype = "torch.float64")),
                      `_module` = "torch", `_einsum` = "einsum", `_add` = "add", `_transpose` = "t",
                      `_ones` = "ones", `_zeros` = "zeros"),
         init = function(p, settings, nxt) {
           nxt(settings)
           p$requires_grad <- .cgp_s(p, "requires_grad")
           p$dtype <- .cgp_s(p, "dtype")
           p$torch_version <- .cgp_s(p, "torch_version")
           invisible(NULL)
         },
         methods = list(
           `_print_Function` = .prt_torch_fun, `_print_Expr` = .prt_torch_fun, `_print_Application` = .prt_torch_fun,
           `_print_MatrixExpr` = .prt_torch_fun, `_print_Relational` = .prt_torch_fun, `_print_Not` = .prt_torch_fun,
           `_print_And` = .prt_torch_fun, `_print_Or` = .prt_torch_fun, `_print_HadamardProduct` = .prt_torch_fun,
           `_print_Trace` = .prt_torch_fun, `_print_Determinant` = .prt_torch_fun,
           `_print_Inverse` = function(p, e) sprintf("%s(%s)", .cgp_module_format(p, "torch.linalg.inv"), .cgpp(p, e$args[[1L]])),
           `_print_Transpose` = function(p, e) {
             a <- e$args[[1L]]
             sh <- .prt_shape(a)
             if (.cgp_is_matrix(a) && identical(sh[[1L]]$key, sh[[2L]]$key))
               return(sprintf("torch.transpose(%s).t()", .cgpp(p, a)))
             idx <- paste(seq_along(sh) - 1L, collapse = ", ")
             sprintf("%s.permute(%s)", .cgpp(p, a), paste(rev(strsplit(idx, "", fixed = TRUE)[[1L]]), collapse = ""))
           },
           `_print_PermuteDims` = function(p, e)
             sprintf("%s.permute(%s)", .cgpp(p, e$args[[1L]]), paste(.ae_paf(e$args[[2L]]), collapse = ", ")),
           `_print_Derivative` = function(p, e) {
             vars <- .prt_deriv_vars(e)
             arg <- e$args[[1L]]
             if (length(vars) > 1L) {
               res <- .cgpp(p, arg)
               keys <- character(); groups <- list(); counts <- integer()
               for (v in vars) {
                 k <- match(v$key, keys)
                 if (is.na(k)) { keys <- c(keys, v$key); groups[[length(groups) + 1L]] <- v; counts <- c(counts, 1L) }
                 else counts[k] <- counts[k] + 1L
               }
               for (g in seq_along(groups)) for (r in seq_len(counts[g])) res <- .prt_grad(res, .cgpp(p, groups[[g]]))
               return(res)
             }
             if (length(vars) == 1L) return(sprintf("torch.autograd.grad(%s, %s)[0]", .cgpp(p, arg), .cgpp(p, vars[[1L]])))
             .cgpp(p, arg)
           },
           `_print_Piecewise` = function(p, e) .prt_where(p, e, "torch.where"),
           `_print_Pow` = function(p, e, ...) {
             if (identical(e$args[[2L]]$key, S_Half$key))
               return(sprintf("%s(%s)", .cgp_module_format(p, "torch.sqrt"), .cgpp(p, e$args[[1L]])))
             sprintf("%s(%s, %s)", .cgp_module_format(p, "torch.pow"), .cgpp(p, e$args[[1L]]), .cgpp(p, e$args[[2L]]))
           },
           `_print_MatMul` = .prt_matmul("torch.matmul"),
           `_print_MatPow` = .prt_matpow("torch.mm"),
           `_print_MatrixBase` = function(p, e) {
             params <- c(.prt_matrix_data(p, e), paste0("dtype=", p$dtype))
             if (isTRUE(p$requires_grad)) params <- c(params, "requires_grad=True")
             sprintf("%s(%s)", .cgp_module_format(p, "torch.tensor"), paste(params, collapse = ", "))
           },
           `_print_isnan` = function(p, e) sprintf("torch.isnan(%s)", .cgpp(p, e$args[[1L]])),
           `_print_isinf` = function(p, e) sprintf("torch.isinf(%s)", .cgpp(p, e$args[[1L]])),
           `_print_Identity` = function(p, e) {
             sh <- .mshape(e)
             if (all(vapply(sh, function(d) d$type == "Integer", NA)))
               return(sprintf("%s(%s)", .cgp_module_format(p, "torch.eye"), .cgpp(p, sh[[1L]])))
             sprintf("%s(%s, %s)", .cgp_module_format(p, "torch.eye"), .cgpp(p, sh[[1L]]), .cgpp(p, sh[[2L]]))
           },
           `_print_ZeroMatrix` = function(p, e) sprintf("%s(%s)", .cgp_module_format(p, "torch.zeros"), .cgpp(p, do.call(Tuple, .mshape(e)))),
           `_print_OneMatrix` = function(p, e) sprintf("%s(%s)", .cgp_module_format(p, "torch.ones"), .cgpp(p, do.call(Tuple, .mshape(e)))),
           `_print_conjugate` = function(p, e) sprintf("%s(%s)", .cgp_module_format(p, "torch.conj"), .cgpp(p, e$args[[1L]])),
           `_print_ImaginaryUnit` = function(p, e) "1j",
           `_print_Heaviside` = function(p, e) {
             # SymPy keeps the default H0 = 1/2 as the second argument
             a <- c(.cgpp(p, e$args[[1L]]), .cgpp(p, if (length(e$args) > 1L) e$args[[2L]] else S_Half))
             sprintf("%s(%s, %s)", .cgp_module_format(p, "torch.heaviside"), a[1L], a[2L])
           },
           `_print_gamma` = function(p, e) sprintf("%s(%s)", .cgp_module_format(p, "torch.special.gamma"), .cgpp(p, e$args[[1L]])),
           `_print_polygamma` = function(p, e) {
             if (identical(e$args[[1L]]$key, S_Zero$key))
               return(sprintf("%s(%s)", .cgp_module_format(p, "torch.special.digamma"), .cgpp(p, e$args[[2L]])))
             .cg_err("NotImplementedError", "PyTorch only supports digamma (0th order polygamma)")
           }
         ))

#' @rdname tensorflow-torch-code
#' @export
TorchPrinter <- function(settings = NULL) .cgp_new("TorchPrinter", settings)

#' @rdname tensorflow-torch-code
#' @export
torch_code <- function(expr, requires_grad = FALSE, dtype = "torch.float64", assign_to = NULL)
  printer_doprint(TorchPrinter(list(requires_grad = requires_grad, dtype = dtype)), expr, assign_to = assign_to)
