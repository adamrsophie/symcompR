# Rendering of abstract syntax trees as source files, after
# sympy/codegen/cutils.py, futils.py and pyutils.py.

.cg_printer_obj <- function(printer, settings) {
  if (inherits(printer, "cg_printer")) return(printer)
  if (is.character(printer)) return(.cgp_new(printer, settings))
  printer(settings)
}

#' Source files and modules from abstract syntax trees
#'
#' `render_as_source_file(content)` prints `content` as C source preceded
#' by the `#include` lines it needs (SymPy's `cutils.render_as_source_file`).
#' `py_render_as_module(content)` prints it as a Python module preceded by
#' its imports (SymPy's `pyutils.render_as_module`), and
#' `f_render_as_module(definitions, name)` prints Fortran definitions as a
#' module (SymPy's `futils.render_as_module`).
#'
#' @param content An expression or abstract syntax tree.
#' @param Printer A printer constructor such as [C99CodePrinter()], or a
#'   printer class name.
#' @param settings `NULL`, or a list of printer settings.
#' @param standard The Python standard.
#' @param definitions A list of definitions, such as [FunctionDefinition()]
#'   nodes.
#' @param name The name of the module.
#' @param declarations A list of declarations.
#' @param printer_settings `NULL` (Fortran 2003 in free form), or a list of
#'   printer settings.
#' @return A character scalar.
#' @examples
#' x <- Symbol("x", real = TRUE)
#' cat(render_as_source_file(Assignment(x, sin(x))), "\n")
#' cat(py_render_as_module(Print(list(x), "%f")), "\n")
#' @export
render_as_source_file <- function(content, Printer = C99CodePrinter, settings = NULL) {
  printer <- .cg_printer_obj(Printer, if (is.null(settings)) list() else settings)
  code_str <- .cgp_doprint(printer, content, NULL)
  includes <- paste(sprintf("#include <%s>", unclass(printer$headers)), collapse = "\n")
  paste0(includes, "\n\n", code_str)
}

#' @rdname render_as_source_file
#' @export
py_render_as_module <- function(content, standard = "python3") {
  printer <- .cgp_new("PythonCodePrinter", list(standard = standard))
  pystr <- .cgp_doprint(printer, content, NULL)
  mi <- unclass(printer$module_imports)
  imports <- if (isTRUE(.cgp_s(printer, "fully_qualified_modules")))
    sprintf("import %s", names(mi))
  else vapply(names(mi), function(k) sprintf("from %s import %s", k, paste(mi[[k]], collapse = ", ")), "")
  paste0(paste(imports, collapse = "\n"), "\n\n", pystr)
}

#' @rdname render_as_source_file
#' @export
f_render_as_module <- function(definitions, name, declarations = list(), printer_settings = NULL) {
  if (is.null(printer_settings)) printer_settings <- list(standard = 2003, source_format = "free")
  printer <- .cgp_new("FCodePrinter", printer_settings)
  dummy <- Dummy()
  if (is_sym(definitions) && definitions$type == "Module")
    .cg_err("ValueError", "This function expects to construct a module on its own.")
  mod <- Module(name, c(as.list(declarations), list(dummy)), definitions)
  fstr <- .cgp_doprint(printer, mod, NULL)
  mu <- printer$module_uses
  uses <- vapply(names(mu), function(k) sprintf("use %s, only: %s", k, paste(mu[[k]], collapse = ", ")), "")
  module_use_str <- sprintf("   %s\n", paste(uses, collapse = "   \n"))
  module_use_str <- paste0(module_use_str, "   implicit none\n", "   private\n")
  defs <- if (is_sym(definitions)) definitions$args else definitions
  nms <- unlist(lapply(defs, function(d) {
    n <- tryCatch(ast_attr(d, "name"), error = function(err) NULL)
    if (is.null(n) || .is_none(n)) NULL else .cg_text(n)
  }))
  module_use_str <- paste0(module_use_str, sprintf("   public %s\n", paste(nms, collapse = ", ")))
  gsub(.cgp_doprint(printer, dummy, NULL), module_use_str, fstr, fixed = TRUE)
}
