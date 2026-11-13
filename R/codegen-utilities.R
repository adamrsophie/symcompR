# Routines and source files, after sympy/utilities/codegen.py.

.cgu_obj <- function(cls, ...) structure(list(...), class = c(cls, "cg_cgobj"))
.cgu_is <- function(x, cls) inherits(x, cls)
.cgu_cls <- c(InputArgument = "Argument Variable", OutputArgument = "Argument ResultBase Variable",
              InOutArgument = "Argument ResultBase Variable", Result = "Variable ResultBase")

#' Data types and arguments of code generation routines
#'
#' `DataType()` holds the names of a data type in the target languages, and
#' `get_default_datatype()` chooses one for an expression, as in
#' `sympy.utilities.codegen`. `InputArgument()`, `OutputArgument()`,
#' `InOutArgument()` and `Result()` describe the arguments and results of a
#' routine made by [make_routine()]; their fields are read with
#' [ast_attr()].
#'
#' @param cname,fname,pyname,jlname,octname,rsname The type names in C,
#'   Fortran, Python, Julia, Octave and Rust.
#' @param expr An expression.
#' @param complex_allowed Whether a complex type is allowed. `NULL` reads the
#'   option `symcompR.complex_allowed` (default `FALSE`), which plays the
#'   part of SymPy's module flag `COMPLEX_ALLOWED`.
#' @param name A symbol (or, for `Result()`, a name).
#' @param result_var The variable or element the result is stored in.
#' @param datatype A `DataType()`.
#' @param dimensions A list of `list(lower, upper)` bounds.
#' @param precision The precision, or `NULL`.
#' @return An object of the corresponding class.
#' @examples
#' x <- Symbol("x")
#' get_default_datatype(x)
#' InputArgument(x)
#' @name codegen-datatypes
NULL

#' @rdname codegen-datatypes
#' @export
DataType <- function(cname, fname, pyname, jlname, octname, rsname)
  .cgu_obj("DataType", cname = cname, fname = fname, pyname = pyname, jlname = jlname, octname = octname, rsname = rsname)

.cgu_default_datatypes <- function() list(int = DataType("int", "INTEGER*4", "int", "", "", "i32"),
                                         float = DataType("double", "REAL*8", "float", "", "", "f64"),
                                         complex = DataType("double", "COMPLEX*16", "complex", "", "", "float"))

#' @rdname codegen-datatypes
#' @export
get_default_datatype <- function(expr, complex_allowed = NULL) {
  if (is.null(complex_allowed)) complex_allowed <- isTRUE(getOption("symcompR.complex_allowed", FALSE))
  dts <- .cgu_default_datatypes()
  final <- if (isTRUE(complex_allowed)) "complex" else "float"
  if (inherits(expr, "symmatrix")) {
    dt <- "int"
    for (el in expr$data) {
      if (dt == "int" && !isTRUE(ask_fact(el, "integer"))) dt <- "float"
      if (dt == "float" && !isTRUE(ask_fact(el, "real"))) return(dts[[final]])
    }
    return(dts[[dt]])
  }
  if (is_sym(expr) && !(.is_matexpr(expr))) {
    if (isTRUE(ask_fact(expr, "integer"))) return(dts$int)
    if (isTRUE(ask_fact(expr, "real"))) return(dts$float)
  }
  dts[[final]]
}

.cgu_var <- function(cls, name, datatype = NULL, dimensions = NULL, precision = NULL, ...) {
  if (!(is_sym(name) && name$type %in% c("Symbol", "Dummy", "MatrixSymbol")))
    .cg_err("TypeError", "The first argument must be a SymPy symbol.")
  if (is.null(datatype)) datatype <- get_default_datatype(name)
  else if (!inherits(datatype, "DataType"))
    .cg_err("TypeError", "The (optional) `datatype' argument must be an instance of the DataType class.")
  if (length(dimensions) && !is.list(dimensions)) .cg_err("TypeError", "The dimensions argument must be a sequence of tuples")
  dt <- list(C = datatype$cname, FORTRAN = datatype$fname, JULIA = datatype$jlname, OCTAVE = datatype$octname,
             PYTHON = datatype$pyname, RUST = datatype$rsname)
  structure(c(list(name = name, datatype = dt, dimensions = dimensions, precision = precision), list(...)),
            class = c(cls, strsplit(.cgu_cls[[cls]], " ")[[1L]], "cg_cgobj"))
}

#' @rdname codegen-datatypes
#' @export
InputArgument <- function(name, datatype = NULL, dimensions = NULL, precision = NULL)
  .cgu_var("InputArgument", name, datatype, dimensions, precision)
#' @rdname codegen-datatypes
#' @export
OutputArgument <- function(name, result_var, expr, datatype = NULL, dimensions = NULL, precision = NULL)
  .cgu_var("OutputArgument", name, datatype, dimensions, precision, expr = expr, result_var = result_var)
#' @rdname codegen-datatypes
#' @export
InOutArgument <- function(name, result_var, expr, datatype = NULL, dimensions = NULL, precision = NULL) {
  if (is.null(datatype)) datatype <- get_default_datatype(expr)
  .cgu_var("InOutArgument", name, datatype, dimensions, precision, expr = expr, result_var = result_var)
}
#' @rdname codegen-datatypes
#' @export
Result <- function(expr, name = NULL, result_var = NULL, datatype = NULL, dimensions = NULL, precision = NULL) {
  if (!(is_sym(expr) || inherits(expr, "symmatrix"))) .cg_err("TypeError", "The first argument must be a SymPy expression.")
  if (is.null(name)) name <- sprintf("result_%s", .cgu_hash(expr))
  if (is.null(datatype)) datatype <- get_default_datatype(expr)
  if (is.character(name)) {
    name <- if (inherits(expr, "symmatrix") || .is_matexpr(expr)) {
      sh <- .mshape(expr)
      MatrixSymbol(name, sh[[1L]], sh[[2L]])
    } else Symbol(name)
  }
  if (is.null(result_var)) result_var <- name
  .cgu_var("Result", name, datatype, dimensions, precision, expr = expr, result_var = result_var)
}
# abs(hash(expr)): any stable number derived from the expression
.cgu_hash <- function(expr) {
  k <- if (is_sym(expr)) expr$key else paste(vapply(expr$data, .k, ""), collapse = ",")
  v <- utf8ToInt(k)
  h <- 0
  for (c in v) h <- (h * 31 + c) %% 2147483647
  format(h, scientific = FALSE)
}

.cgu_get_datatype <- function(v, language) {
  dt <- v$datatype[[toupper(language)]]
  if (is.null(dt)) .cg_err("CodeGenError", sprintf("Has datatypes for languages: %s", paste(names(v$datatype), collapse = ", ")))
  dt
}

#' @export
format.cg_cgobj <- function(x, ...) {
  r <- function(v) if (is.character(v)) .py_repr_str(v) else if (is.null(v)) "None" else if (is.list(v) && !is_sym(v)) .cgu_fmt_list(v) else .py_str_of(v)
  cls <- class(x)[1L]
  switch(cls,
         DataType = sprintf("<DataType %s>", x$cname),
         InputArgument = sprintf("%s(%s)", cls, r(x$name)),
         OutputArgument = , InOutArgument = sprintf("%s(%s, %s, %s)", cls, r(x$name), r(x$result_var), r(x$expr)),
         Result = sprintf("%s(%s, %s, %s)", cls, r(x$expr), r(x$name), r(x$result_var)),
         Routine = sprintf("Routine(%s, %s, %s, %s, %s)", .py_repr_str(x$name), .cgu_fmt_list(x$arguments),
                           .cgu_fmt_list(x$results), .cgu_fmt_set(x$local_vars), .cgu_fmt_set(x$global_vars)),
         cls)
}
.cgu_fmt_list <- function(l) paste0("[", paste(vapply(l, function(z) if (inherits(z, "cg_cgobj")) format(z) else .py_str_of(z), ""), collapse = ", "), "]")
.cgu_fmt_set <- function(l) if (!length(l)) "set()" else paste0("{", paste(vapply(l, function(z) if (inherits(z, "cg_cgobj")) format(z) else .py_str_of(z), ""), collapse = ", "), "}")
#' @export
print.cg_cgobj <- function(x, ...) {
  cat(format(x), "\n")
  invisible(x)
}

.cgu_free <- function(e) {
  fs <- .cgu_freeall(list(e))
  ix <- vapply(.cgu_atoms(list(e), "Indexed"), .k, "")
  Filter(function(s) !(s$key %in% ix), fs)
}
.cgu_key <- function(z) if (is_sym(z)) z$key else if (inherits(z, "cg_cgobj")) paste0("\001", format(z)) else paste(deparse(z), collapse = "")
.cgu_uniq <- function(l) l[!duplicated(vapply(l, .cgu_key, ""))]
.cgu_in <- function(s, l) .cgu_key(s) %in% vapply(l, .cgu_key, "")
.cgu_minus <- function(a, b) Filter(function(s) !.cgu_in(s, b), a)
.cgu_label <- function(s) if (.is_idx(s)) .cg_idx_label(s) else s

#' Routines for code generation
#'
#' `Routine()` holds the arguments, results and variables of one routine, as
#' SymPy's `sympy.utilities.codegen.Routine`. `make_routine()` builds one
#' from an expression, as SymPy's `make_routine()`.
#'
#' @param name The routine name.
#' @param arguments,results Lists of arguments and results.
#' @param local_vars,global_vars Lists of local and global variables.
#' @param expr An expression, an equation, or a list of them.
#' @param argument_sequence The order of the arguments, or `NULL`.
#' @param language The language whose conventions to follow.
#' @return A routine (class `Routine`); its fields are read with
#'   [ast_attr()].
#' @examples
#' x <- Symbol("x"); y <- Symbol("y")
#' r <- make_routine("f", x^2 + y)
#' ast_attr(r, "arguments")
#' @export
Routine <- function(name, arguments, results, local_vars, global_vars) {
  input_symbols <- list(); symbols <- list()
  for (arg in arguments) {
    if (inherits(arg, "OutputArgument")) symbols <- c(symbols, .cgu_free(arg$expr))
    else if (inherits(arg, "InputArgument")) input_symbols <- c(input_symbols, list(arg$name))
    else if (inherits(arg, "InOutArgument")) {
      input_symbols <- c(input_symbols, list(arg$name)); symbols <- c(symbols, .cgu_free(arg$expr))
    } else .cg_err("ValueError", sprintf("Unknown Routine argument: %s", format(arg)))
  }
  for (r in results) {
    if (!inherits(r, "Result")) .cg_err("ValueError", sprintf("Unknown Routine result: %s", format(r)))
    symbols <- c(symbols, .cgu_free(r$expr))
  }
  local_symbols <- list()
  for (r in local_vars) {
    if (inherits(r, "Result")) { symbols <- c(symbols, .cgu_free(r$expr)); local_symbols <- c(local_symbols, list(r$name)) }
    else local_symbols <- c(local_symbols, list(r))
  }
  symbols <- .cgu_uniq(lapply(symbols, .cgu_label))
  notcovered <- .cgu_minus(symbols, c(input_symbols, local_symbols, global_vars))
  if (length(notcovered))
    .cg_err("ValueError", paste0("Symbols needed for output are not in input ", paste(vapply(notcovered, sstr, ""), collapse = ", ")))
  .cgu_obj("Routine", name = name, arguments = arguments, results = results, local_vars = local_vars, global_vars = global_vars)
}
.cgu_variables <- function(r) .cgu_uniq(c(Filter(is_sym, r$local_vars), lapply(r$arguments, `[[`, "name"),
                                          lapply(r$results, `[[`, "result_var")))
.cgu_result_variables <- function(r) c(Filter(function(a) inherits(a, c("OutputArgument", "InOutArgument")), r$arguments), r$results)

.cgu_attr <- function(x, name) {
  if (inherits(x, "cg_codegen")) return(.cgu_cg_attr(x, name))
  if (inherits(x, "Routine") && name == "variables") return(.cgu_variables(x))
  if (inherits(x, "Routine") && name == "result_variables") return(.cgu_result_variables(x))
  if (name == "get_datatype") return(function(language) .cgu_get_datatype(x, language))
  # Variable keeps its datatypes in the attribute _datatype, a dict by language
  if (name == "_datatype" && !is.null(x[["datatype"]])) return(x[["datatype"]])
  if (name %in% setdiff(names(x), "datatype")) return(x[[name]])
  .cg_err("AttributeError", sprintf("'%s' object has no attribute '%s'", class(x)[1L], name))
}

.cgu_seq <- function(expr) is.list(expr) && !is_sym(expr) && !inherits(expr, "symmatrix")
.cgu_atoms <- function(exprs, type) .cgu_uniq(do.call(c, c(list(list()), lapply(exprs, function(e)
  if (inherits(e, "symmatrix")) do.call(c, c(list(list()), lapply(e$data, atoms, types = type))) else atoms(e, type)))))
# free symbols, with the symbols in the shapes of IndexedBase objects
.cgu_freeall <- function(exprs) {
  fs <- do.call(c, c(list(list()), lapply(exprs, function(e)
    if (inherits(e, "symmatrix")) do.call(c, c(list(list()), lapply(e$data, free_symbols))) else free_symbols(e))))
  shp <- do.call(c, c(list(list()), lapply(.cgu_atoms(exprs, "Indexed"), function(a) a$args[[1L]]$extra$shape)))
  .cgu_uniq(c(fs, do.call(c, c(list(list()), lapply(shp, free_symbols)))))
}
.cgu_sort_str <- function(l) l[order(vapply(l, function(z) .py_str_of(z), ""), method = "radix")]
.cgu_dims0 <- function(s) lapply(.cgu_shape(s), function(d) list(S_Zero, .sub(as_sym(d), S_One)))
.cgu_shape <- function(s) {
  if (inherits(s, "symmatrix")) return(list(Integer(s$nrow), Integer(s$ncol)))
  if (s$type == "Indexed") {
    if (!is.null(s$args[[1L]]$extra$shape)) return(s$args[[1L]]$extra$shape)
    return(lapply(s$args[-1L], function(i) .add(list(.cgp_idx_upper(i), S_One, .neg(.cgp_idx_lower(i))))))
  }
  if (s$type == "IndexedBase") return(ast_attr(s, "shape"))
  .mshape(s)
}
.cgu_has <- function(e, s) if (inherits(e, "symmatrix")) any(vapply(e$data, function(z) has(z, s), TRUE)) else has(e, s)
.cgu_eq_parts <- function(e) list(e$args[[1L]], e$args[[2L]])

.cgu_argseq <- function(arg_list, argument_sequence, make_input) {
  if (is.null(argument_sequence)) return(arg_list)
  argument_sequence <- lapply(argument_sequence, function(a) {
    if (is.function(a) && !is.null(attr(a, "base"))) a <- attr(a, "base")
    if (is_sym(a) && a$type == "IndexedBase") a$args[[1L]] else a
  })
  names_ <- lapply(arg_list, `[[`, "name")
  missing <- Filter(function(x) !.cgu_in(x$name, argument_sequence), arg_list)
  if (length(missing)) {
    msg <- sprintf("Argument list didn't specify: %s ", paste(vapply(missing, function(m) sstr(m$name), ""), collapse = ", "))
    .cg_err("CodeGenArgumentListError", sprintf("(%s, %s)", .py_repr_str(msg), .cgu_fmt_list(missing)))
  }
  lapply(argument_sequence, function(s) {
    hit <- which(vapply(names_, function(n) identical(n$key, s$key), TRUE))
    if (length(hit)) arg_list[[hit[1L]]] else make_input(s)
  })
}

# cse() of a list, with explicit matrices flattened and rebuilt as SymPy's
# cse() does
.cgu_cse <- function(exprs) {
  unwrap <- function(e) if (is_sym(e) && e$type == "ImmutableMatrix") as_explicit(e) else e
  exprs <- lapply(exprs, unwrap)
  shp <- lapply(exprs, function(e) if (inherits(e, "symmatrix")) c(e$nrow, e$ncol) else NULL)
  flat <- lapply(exprs, function(e) if (inherits(e, "symmatrix")) do.call(Tuple, e$data) else e)
  r <- cse(flat)
  r$reduced <- Map(function(z, s) if (is.null(s)) z else .new_mat(z$args, s[1L], s[2L]), r$reduced, shp)
  r
}
.cgu_routine_generic <- function(cg, name, expr, argument_sequence = NULL, global_vars = NULL) {
  local_vars <- list(); local_symbols <- list()
  local_expressions <- list()
  if (isTRUE(cg$cse)) {
    if (.cgu_seq(expr)) {
      if (!length(expr)) .cg_err("ValueError", "No expression given")
      for (e in expr) if (!(is_sym(e) && e$type == "Equality"))
        .cg_err("CodeGenError", sprintf("Lists of expressions must all be Equalities. %s is not.", .py_str_of(e)))
      r <- .cgu_cse(lapply(expr, function(e) e$args[[2L]]))
      expr <- Map(function(e, rhs) Eq(e$args[[1L]], rhs), expr, r$reduced)
    } else if (is_sym(expr) && expr$type == "Equality") {
      r <- .cgu_cse(list(expr$args[[2L]])); expr <- Eq(expr$args[[1L]], r$reduced[[1L]])
    } else { r <- .cgu_cse(list(expr)); expr <- r$reduced[[1L]] }
    local_vars <- lapply(r$replacements, function(pr) Result(pr[[2L]], pr[[1L]]))
    local_symbols <- lapply(r$replacements, `[[`, 1L)
    local_expressions <- lapply(r$replacements, `[[`, 2L)
  }
  if (.cgu_seq(expr)) {
    if (!length(expr)) .cg_err("ValueError", "No expression given")
    expressions <- expr
  } else expressions <- list(expr)
  if (isTRUE(cg$cse)) {
    if (length(.cgu_atoms(expressions, "Idx"))) .cg_err("CodeGenError", "CSE and Indexed expressions do not play well together yet")
  } else {
    local_vars <- .cgu_uniq(lapply(Filter(.is_idx, .cgu_atoms(expressions, c("Idx", "Symbol"))), .cgu_label))
    local_symbols <- local_vars
  }
  global_vars <- if (is.null(global_vars)) list() else global_vars
  symbols <- .cgu_minus(.cgu_freeall(c(expressions, local_expressions)), c(local_symbols, global_vars))
  new_symbols <- list()
  for (s in symbols) {
    if (.is_idx(s)) new_symbols <- c(new_symbols, free_symbols(.cgp_idx_upper(s)), free_symbols(.cgp_idx_lower(s)))
    else if (s$type != "Indexed") new_symbols <- c(new_symbols, list(s))
  }
  symbols <- .cgu_uniq(new_symbols)
  return_val <- list(); output_args <- list()
  for (e in expressions) {
    if (is_sym(e) && e$type == "Equality") {
      out_arg <- e$args[[1L]]; e <- e$args[[2L]]
      if (out_arg$type == "Indexed") { dims <- .cgu_dims0(out_arg); symbol <- out_arg$args[[1L]]$args[[1L]] }
      else if (out_arg$type %in% c("Symbol", "Dummy")) { dims <- list(); symbol <- out_arg }
      else if (out_arg$type == "MatrixSymbol") { dims <- .cgu_dims0(out_arg); symbol <- out_arg }
      else .cg_err("CodeGenError", "Only Indexed, Symbol, or MatrixSymbol can define output arguments.")
      if (.cgu_has(e, symbol)) output_args[[length(output_args) + 1L]] <- InOutArgument(symbol, out_arg, e, dimensions = dims)
      else output_args[[length(output_args) + 1L]] <- OutputArgument(symbol, out_arg, e, dimensions = dims)
      if (!.cgu_in(symbol, local_vars)) symbols <- .cgu_minus(symbols, list(symbol))
    } else if (inherits(e, "symmatrix") || (is_sym(e) && e$type %in% c("ImmutableMatrix", "MatrixSlice"))) {
      sh <- .cgu_shape(e)
      out_arg <- MatrixSymbol(sprintf("out_%s", .cgu_hash(e)), sh[[1L]], sh[[2L]])
      output_args[[length(output_args) + 1L]] <- OutputArgument(out_arg, out_arg, e, dimensions = .cgu_dims0(out_arg))
    } else return_val[[length(return_val) + 1L]] <- Result(e)
  }
  array_symbols <- list()
  for (a in .cgu_atoms(c(expressions, local_expressions), "Indexed")) array_symbols[[a$args[[1L]]$args[[1L]]$key]] <- a
  for (a in .cgu_atoms(c(expressions, local_expressions), "MatrixSymbol")) array_symbols[[a$key]] <- a
  arg_list <- lapply(.cgu_sort_str(symbols), function(s) {
    if (s$key %in% names(array_symbols)) InputArgument(s, dimensions = .cgu_dims0(array_symbols[[s$key]])) else InputArgument(s)
  })
  output_args <- output_args[order(vapply(output_args, function(a) sstr(a$name), ""), method = "radix")]
  arg_list <- c(arg_list, output_args)
  arg_list <- .cgu_argseq(arg_list, argument_sequence, function(s) {
    if (s$type %in% c("IndexedBase", "MatrixSymbol")) InputArgument(s, dimensions = .cgu_dims0(s)) else InputArgument(s)
  })
  Routine(name, arg_list, return_val, local_vars, global_vars)
}

# JuliaCodeGen, OctaveCodeGen and RustCodeGen routines
.cgu_routine_results <- function(cg, name, expr, argument_sequence, global_vars, lang) {
  expressions <- if (.cgu_seq(expr)) { if (!length(expr)) .cg_err("ValueError", "No expression given"); expr } else list(expr)
  local_vars <- .cgu_uniq(lapply(Filter(.is_idx, .cgu_atoms(expressions, c("Idx", "Symbol"))), .cgu_label))
  global_vars <- if (is.null(global_vars)) list() else global_vars
  old <- .cgu_minus(.cgu_freeall(expressions), c(local_vars, global_vars))
  symbols <- list()
  if (lang == "rust") symbols <- Filter(function(s) s$type != "Indexed", old)
  else for (s in old) {
    if (.is_idx(s)) symbols <- c(symbols, free_symbols(.cgp_idx_upper(s)), free_symbols(.cgp_idx_lower(s)))
    else if (s$type != "Indexed") symbols <- c(symbols, list(s))
  }
  symbols <- .cgu_uniq(symbols)
  return_vals <- list(); output_args <- list()
  for (i in seq_along(expressions)) {
    e <- expressions[[i]]
    if (is_sym(e) && e$type == "Equality") {
      out_arg <- e$args[[1L]]; e <- e$args[[2L]]; symbol <- out_arg
      if (out_arg$type == "Indexed") {
        symbol <- out_arg$args[[1L]]$args[[1L]]
        if (lang != "octave") {
          dims <- lapply(.cgu_shape(out_arg), function(d) list(S_One, d))
          output_args[[length(output_args) + 1L]] <- InOutArgument(symbol, out_arg, e, dimensions = dims)
        }
      }
      if (!(out_arg$type %in% c("Indexed", "Symbol", "Dummy", "MatrixSymbol")))
        .cg_err("CodeGenError", "Only Indexed, Symbol, or MatrixSymbol can define output arguments.")
      return_vals[[length(return_vals) + 1L]] <- Result(e, name = symbol, result_var = out_arg)
      if (!.cgu_has(e, symbol)) {
        if (!.cgu_in(symbol, symbols)) .cg_err("KeyError", sstr(symbol))
        symbols <- .cgu_minus(symbols, list(symbol))
      }
    } else return_vals[[length(return_vals) + 1L]] <- Result(e, name = sprintf("out%d", i))
  }
  output_args <- output_args[order(vapply(output_args, function(a) sstr(a$name), ""), method = "radix")]
  arg_list <- c(output_args, lapply(.cgu_sort_str(symbols), InputArgument))
  arg_list <- .cgu_argseq(arg_list, argument_sequence, InputArgument)
  Routine(name, arg_list, return_vals, local_vars, global_vars)
}

# ---------------------------------------------------------------------------
# Code generators: environments of class "cg_codegen"

.cgu_header_comment <- function(project)
  sprintf("Code generated with symcompR %s\n\nSee http://www.sympy.org/ for more information.\n\nThis file is part of '%s'\n",
          as.character(utils::packageVersion("symcompR")), project)
.cgu_new <- function(cls, project = "project", printer = NULL, cse = FALSE, preprocessor_statements = NULL) {
  cg <- new.env()
  cg$cls <- cls; cg$project <- project; cg$cse <- cse
  cg$printer <- if (!is.null(printer)) printer else switch(cls,
    CCodeGen = , C99CodeGen = .cgp_new("C99CodePrinter"), C89CodeGen = .cgp_new("C89CodePrinter"),
    FCodeGen = .cgp_new("FCodePrinter"), JuliaCodeGen = .cgp_new("JuliaCodePrinter"),
    OctaveCodeGen = .cgp_new("OctaveCodePrinter"), RustCodeGen = .cgp_new("RustCodePrinter"))
  cg$preprocessor_statements <- if (is.null(preprocessor_statements)) "#include <math.h>" else preprocessor_statements
  class(cg) <- c(cls, "cg_codegen")
  cg
}
.cgu_with_settings <- function(cg, settings, f) {
  p <- cg$printer
  ori <- p$settings[names(settings)]
  had <- names(settings) %in% names(p$settings)
  for (k in names(settings)) p$settings[k] <- list(settings[[k]])
  on.exit({ for (i in seq_along(settings)) { k <- names(settings)[i]; if (had[i]) p$settings[k] <- list(ori[[k]]) else p$settings[[k]] <- NULL } })
  f(p)
}
.cgu_doprint_tuple <- function(cg, settings, expr, assign_to) .cgu_with_settings(cg, settings, function(p) .cgp_doprint(p, expr, assign_to))
.cgu_get_symbol <- function(cg, s) {
  p <- cg$printer
  if (isTRUE(.cgp_s(p, "human", TRUE))) return(trimws(.cgp_doprint(p, s)))
  r <- .cgp_doprint(p, s)
  if (length(r[[1L]]) || length(r[[2L]])) .cg_err("ValueError", sprintf("Failed to print %s", .py_str_of(s)))
  trimws(r[[3L]])
}
.cgu_lang <- function(cg) switch(cg$cls, CCodeGen = , C89CodeGen = , C99CodeGen = "c", FCodeGen = "f", JuliaCodeGen = "julia",
                                 OctaveCodeGen = "octave", RustCodeGen = "rust")
.cgu_indent_str <- function(p, s, settings = list()) {
  old <- p$settings[names(settings)]
  had <- names(settings) %in% names(p$settings)
  for (k in names(settings)) p$settings[k] <- list(settings[[k]])
  on.exit({ for (i in seq_along(settings)) { k <- names(settings)[i]; if (had[i]) p$settings[k] <- list(old[[k]]) else p$settings[[k]] <- NULL } })
  lines <- regmatches(s, gregexpr("[^\n]*\n|[^\n]+$", s))[[1L]]
  paste(.cgp_call(p, "indent_code", lines), collapse = "")
}
.cgu_write_out <- function(f, s) if (nzchar(s)) cat(s, file = f, sep = "")
.cgu_sorted_pairs <- function(l) l[order(vapply(l, function(q) paste0("(", .py_str_of(q[[1L]]), ", ", .py_str_of(q[[2L]]), ")"), ""), method = "radix")]
.cgu_sorted <- function(l) l[order(vapply(l, .py_str_of, ""), method = "radix")]

.cgu_header <- function(cg) {
  tmp <- .py_splitlines(.cgu_header_comment(cg$project))
  switch(.cgu_lang(cg),
    c = c(paste0("/", strrep("*", 78), "\n"), sprintf(" *%s*\n", vapply(tmp, .py_center, "", width = 76L)), paste0(" ", strrep("*", 78), "/\n")),
    f = c(paste0("!", strrep("*", 78), "\n"), sprintf("!*%s*\n", vapply(tmp, .py_center, "", width = 76L)), paste0("!", strrep("*", 78), "\n")),
    julia = vapply(tmp, function(l) if (l == "") "#\n" else sprintf("#   %s\n", l), "", USE.NAMES = FALSE),
    octave = vapply(tmp, function(l) if (l == "") "%\n" else sprintf("%%   %s\n", l), "", USE.NAMES = FALSE),
    rust = c("/*\n", vapply(tmp, function(l) paste0(sub("\\s+$", "", sprintf(" *%s", .py_center(l, 76L))), "\n"), "", USE.NAMES = FALSE), " */\n"))
}

.cgu_c_prototype <- function(cg, routine) {
  if (length(routine$results) > 1L) .cg_err("CodeGenError", "C only supports a single or no return value.")
  ctype <- if (length(routine$results) == 1L) .cgu_get_datatype(routine$results[[1L]], "C") else "void"
  ta <- vapply(routine$arguments, function(arg) {
    nm <- .cgp_doprint(cg$printer, arg$name)
    if (length(arg$dimensions) || inherits(arg, "ResultBase")) sprintf("%s *%s", .cgu_get_datatype(arg, "C"), nm)
    else sprintf("%s %s", .cgu_get_datatype(arg, "C"), nm)
  }, "")
  sprintf("%s %s(%s)", ctype, routine$name, paste(ta, collapse = ", "))
}
.cgu_rust_prototype <- function(cg, routine) {
  results <- vapply(routine$results, .cgu_get_datatype, "", language = "Rust")
  rstype <- if (length(results) == 1L) paste0(" -> ", results) else if (length(results) > 1L) paste0(" -> (", paste(results, collapse = ", "), ")") else ""
  ta <- vapply(routine$arguments, function(arg) {
    nm <- .cgp_doprint(cg$printer, arg$name)
    if (length(arg$dimensions) || inherits(arg, "ResultBase")) sprintf("*%s: %s", nm, .cgu_get_datatype(arg, "Rust"))
    else sprintf("%s: %s", nm, .cgu_get_datatype(arg, "Rust"))
  }, "")
  sprintf("fn %s(%s)%s", routine$name, paste(ta, collapse = ", "), rstype)
}
.cgu_deref <- function(routine) lapply(Filter(function(a) inherits(a, "ResultBase") && !length(a$dimensions), routine$arguments), `[[`, "name")

.cgu_opening <- function(cg, routine) {
  switch(.cgu_lang(cg),
    c = sprintf("%s {\n", .cgu_c_prototype(cg, routine)),
    rust = sprintf("%s {\n", .cgu_rust_prototype(cg, routine)),
    f = {
      cl <- character()
      if (length(routine$results) > 1L) .cg_err("CodeGenError", "Fortran only supports a single or no return value.")
      if (length(routine$results) == 1L) cl <- c(.cgu_get_datatype(routine$results[[1L]], "fortran"), "function")
      else cl <- "subroutine"
      args <- paste(vapply(routine$arguments, function(a) .cgu_get_symbol(cg, a$name), ""), collapse = ", ")
      call_sig <- sprintf("%s(%s)\n", routine$name, args)
      call_sig <- paste0(paste(.py_textwrap(call_sig, 60L), collapse = " &\n"), "\n")
      c(paste(c(cl, call_sig), collapse = " "), "implicit none\n")
    },
    julia = {
      args <- character()
      for (arg in routine$arguments) {
        if (inherits(arg, "OutputArgument")) .cg_err("CodeGenError", sprintf("Julia: invalid argument of type %s", .cgu_type_str(arg)))
        args <- c(args, .cgu_get_symbol(cg, arg$name))
      }
      sprintf("function %s(%s)\n", routine$name, paste(args, collapse = ", "))
    },
    octave = {
      outs <- vapply(routine$results, function(r) .cgu_get_symbol(cg, r$name), "")
      o <- if (length(outs) > 1L) paste0("[", paste(outs, collapse = ", "), "]") else paste(outs, collapse = "")
      args <- character()
      for (arg in routine$arguments) {
        if (inherits(arg, c("OutputArgument", "InOutArgument")))
          .cg_err("CodeGenError", sprintf("Octave: invalid argument of type %s", .cgu_type_str(arg)))
        args <- c(args, .cgu_get_symbol(cg, arg$name))
      }
      sprintf("function %s = %s(%s)\n", o, routine$name, paste(args, collapse = ", "))
    })
}
.cgu_type_str <- function(x) sprintf("<class 'sympy.utilities.codegen.%s'>", class(x)[1L])

# textwrap.wrap(text, width, break_long_words=False)
.py_textwrap <- function(text, width) {
  words <- strsplit(gsub("\\s+", " ", trimws(text)), " ")[[1L]]
  lines <- character(); cur <- ""
  for (w in words) {
    if (!nzchar(cur)) cur <- w
    else if (nchar(cur) + 1L + nchar(w) <= width) cur <- paste(cur, w)
    else { lines <- c(lines, cur); cur <- w }
  }
  if (nzchar(cur)) lines <- c(lines, cur)
  lines
}

.cgu_declare_arguments <- function(cg, routine) {
  if (.cgu_lang(cg) != "f") return(character())
  arrays <- character(); scalars <- character()
  fprint <- function(s) .cgu_get_symbol(cg, s)
  for (arg in routine$arguments) {
    if (inherits(arg, "InputArgument")) ti <- sprintf("%s, intent(in)", .cgu_get_datatype(arg, "fortran"))
    else if (inherits(arg, "InOutArgument")) ti <- sprintf("%s, intent(inout)", .cgu_get_datatype(arg, "fortran"))
    else if (inherits(arg, "OutputArgument")) ti <- sprintf("%s, intent(out)", .cgu_get_datatype(arg, "fortran"))
    else .cg_err("CodeGenError", sprintf("Unknown Argument type: %s", .cgu_type_str(arg)))
    if (length(arg$dimensions)) {
      dimstr <- paste(vapply(arg$dimensions, function(d) sprintf("%s:%s", fprint(.add(list(d[[1L]], S_One))), fprint(.add(list(d[[2L]], S_One)))), ""), collapse = ", ")
      arrays <- c(arrays, sprintf("%s, dimension(%s) :: %s\n", ti, dimstr, fprint(arg$name)))
    } else scalars <- c(scalars, sprintf("%s :: %s\n", ti, fprint(arg$name)))
  }
  c(scalars, arrays)
}
.cgu_declare_globals <- function(cg, routine) {
  if (.cgu_lang(cg) != "octave" || !length(routine$global_vars)) return(character())
  s <- paste(sort(vapply(routine$global_vars, function(g) .cgu_get_symbol(cg, g), ""), method = "radix"), collapse = " ")
  paste0("global ", s, "\n")
}
.cgu_declare_locals <- function(cg, routine) {
  lang <- .cgu_lang(cg)
  if (lang == "f") {
    return(vapply(.cgu_sorted(routine$local_vars), function(v)
      sprintf("%s :: %s\n", get_default_datatype(v)$fname, .cgu_get_symbol(cg, v)), "", USE.NAMES = FALSE))
  }
  if (lang != "c") return(character())
  deref <- .cgu_deref(routine)
  out <- character()
  for (result in routine$local_vars) {
    if (!inherits(result, "Result")) next
    if (!identical(result$name$key, result$result_var$key))
      .cg_err("CodeGenError", sprintf("Result variable and name should match: %s", format(result)))
    assign_to <- result$name
    t <- .cgu_get_datatype(result, "c")
    if (inherits(result$expr, "symmatrix") || .is_matexpr(result$expr)) {
      dims <- .cgu_shape(result$expr)
      out <- c(out, sprintf("%s %s[%s];\n", t, sstr(assign_to), sstr(.mul(list(dims[[1L]], dims[[2L]])))))
      prefix <- ""
    } else prefix <- sprintf("const %s ", t)
    r <- .cgu_doprint_tuple(cg, list(human = FALSE, dereference = deref, strict = FALSE), result$expr, assign_to)
    for (q in .cgu_sorted_pairs(r[[1L]])) out <- c(out, sprintf("double const %s = %s;\n", .py_str_of(q[[1L]]), .py_str_of(q[[2L]])))
    out <- c(out, sprintf("%s%s\n", prefix, r[[3L]]))
  }
  out
}
.cgu_unsupported_name <- function(obj) if (.cgp_is_function(obj)) .cgp_fname(obj) else .py_str_of(obj)

.cgu_call_printer <- function(cg, routine) {
  lang <- .cgu_lang(cg)
  decl <- character(); code <- character()
  if (lang == "c") {
    deref <- .cgu_deref(routine)
    return_val <- NULL
    st <- list(human = FALSE, dereference = deref, strict = FALSE)
    for (result in .cgu_result_variables(routine)) {
      if (inherits(result, "Result")) {
        assign_to <- Symbol(paste0(routine$name, "_result"))
        code <- c(code, sprintf("%s %s;\n", .cgu_get_datatype(result, "c"), sstr(assign_to)))
        return_val <- sstr(assign_to)
      } else assign_to <- result$result_var
      r <- tryCatch(.cgu_doprint_tuple(cg, st, result$expr, assign_to), AssignmentError = function(err) NULL)
      if (is.null(r)) {
        assign_to <- result$result_var
        code <- c(code, sprintf("%s %s;\n", .cgu_get_datatype(result, "c"), sstr(assign_to)))
        r <- .cgu_doprint_tuple(cg, st, result$expr, assign_to)
      }
      for (q in .cgu_sorted_pairs(r[[1L]])) code <- c(code, sprintf("double const %s = %s;\n", .py_str_of(q[[1L]]), .py_str_of(q[[2L]])))
      code <- c(code, sprintf("%s\n", r[[3L]]))
    }
    if (!is.null(return_val)) code <- c(code, sprintf("   return %s;\n", return_val))
    return(code)
  }
  if (lang == "f") {
    for (result in .cgu_result_variables(routine)) {
      assign_to <- if (inherits(result, "Result")) Symbol(routine$name) else result$result_var
      r <- .cgu_doprint_tuple(cg, list(human = FALSE, source_format = "free", standard = 95, strict = FALSE), result$expr, assign_to)
      for (q in .cgu_sorted_pairs(r[[1L]]))
        decl <- c(decl, sprintf("%s, parameter :: %s = %s\n", get_default_datatype(q[[1L]])$fname, .py_str_of(q[[1L]]), .py_str_of(q[[2L]])))
      for (obj in .cgu_sorted(r[[2L]]))
        decl <- c(decl, sprintf("%s :: %s\n", get_default_datatype(obj)$fname, .cgu_unsupported_name(obj)))
      code <- c(code, sprintf("%s\n", r[[3L]]))
    }
    return(c(decl, code))
  }
  returns <- character()
  deref <- .cgu_deref(routine)
  for (result in routine$results) {
    if (!inherits(result, "Result")) .cg_err("CodeGenError", "unexpected object in Routine results")
    assign_to <- result$result_var
    if (lang == "rust") returns <- c(returns, sstr(result$result_var))
    r <- .cgu_doprint_tuple(cg, list(human = FALSE, strict = FALSE), result$expr, assign_to)
    for (q in .cgu_sorted_pairs(r[[1L]])) decl <- c(decl, switch(lang,
      julia = sprintf("%s = %s\n", .py_str_of(q[[1L]]), .py_str_of(q[[2L]])),
      octave = sprintf("  %s = %s;  %% constant\n", .py_str_of(q[[1L]]), .py_str_of(q[[2L]])),
      rust = sprintf("const %s: f64 = %s;\n", .py_str_of(q[[1L]]), .py_str_of(q[[2L]]))))
    for (obj in .cgu_sorted(r[[2L]])) decl <- c(decl, switch(lang,
      julia = sprintf("# unsupported: %s\n", .cgu_unsupported_name(obj)),
      octave = sprintf("  %% unsupported: %s\n", .cgu_unsupported_name(obj)),
      rust = sprintf("// unsupported: %s\n", .cgu_unsupported_name(obj))))
    code <- c(code, sprintf(if (lang == "rust") "let %s\n" else "%s\n", r[[3L]]))
  }
  if (lang == "rust") {
    if (length(returns) > 1L) returns <- paste0("(", paste(returns, collapse = ", "), ")")
    return(c(decl, code, returns, "\n"))
  }
  c(decl, code)
}
.cgu_ending <- function(cg, routine) {
  switch(.cgu_lang(cg), c = , rust = "}\n",
         f = if (length(routine$results) == 1L) "end function\n" else "end subroutine\n",
         julia = paste0("return ", paste(vapply(routine$results, function(r) .cgu_get_symbol(cg, r$name), ""), collapse = ", "), "\nend\n"),
         octave = "end\n")
}
.cgu_indent <- function(cg, s) {
  switch(.cgu_lang(cg),
         c = , rust = { x <- s; attr(x, "str") <- TRUE; .cgp_call(cg$printer, "indent_code", x) },
         f = .cgu_indent_str(cg$printer, s, list(human = FALSE, source_format = "free", strict = FALSE)),
         julia = .cgu_indent_str(.cgp_new("JuliaCodePrinter", list(human = FALSE, strict = FALSE)), s),
         octave = .cgu_indent_str(cg$printer, s, list(human = FALSE, strict = FALSE)))
}
.cgu_preprocessor <- function(cg, prefix) {
  if (.cgu_lang(cg) != "c") return(character())
  paste0(c(sprintf("#include \"%s.h\"", basename(prefix)), cg$preprocessor_statements), "\n")
}

.cgu_dump_code <- function(cg, routines, f, prefix, header = TRUE, empty = TRUE) {
  lines <- .cgu_preprocessor(cg, prefix)
  for (routine in routines) {
    if (empty) lines <- c(lines, "\n")
    lines <- c(lines, .cgu_opening(cg, routine), .cgu_declare_arguments(cg, routine), .cgu_declare_globals(cg, routine),
               .cgu_declare_locals(cg, routine))
    if (empty) lines <- c(lines, "\n")
    lines <- c(lines, .cgu_call_printer(cg, routine))
    if (empty) lines <- c(lines, "\n")
    lines <- c(lines, .cgu_ending(cg, routine))
  }
  code <- .cgu_indent(cg, paste(lines, collapse = ""))
  if (header) code <- paste0(paste(.cgu_header(cg), collapse = ""), code)
  .cgu_write_out(f, code)
}
.cgu_dump_f95 <- function(cg, routines, f, prefix, header = TRUE, empty = TRUE) {
  for (r in routines) {
    v <- vapply(.cgu_variables(r), .py_str_of, "")
    if (length(unique(tolower(v))) < length(unique(v)))
      .cg_err("CodeGenError", sprintf("Fortran ignores case. Got symbols: %s", paste(v, collapse = ", ")))
  }
  .cgu_dump_code(cg, routines, f, prefix, header, empty)
}
.cgu_dump_c_h <- function(cg, routines, f, prefix, header = TRUE, empty = TRUE) {
  out <- character()
  pr <- function(s = "") out <<- c(out, paste0(s, "\n"))
  if (header) pr(paste(.cgu_header(cg), collapse = ""))
  guard <- sprintf("%s__%s__H", toupper(gsub(" ", "_", cg$project)), toupper(gsub("/", "_", prefix)))
  if (empty) pr()
  pr(sprintf("#ifndef %s", guard)); pr(sprintf("#define %s", guard))
  if (empty) pr()
  for (routine in routines) pr(sprintf("%s;", .cgu_c_prototype(cg, routine)))
  if (empty) pr()
  pr("#endif")
  if (empty) pr()
  cat(paste(out, collapse = ""), file = f, sep = "")
}
.cgu_dump_f_h <- function(cg, routines, f, prefix, header = TRUE, empty = TRUE) {
  out <- character()
  if (header) out <- c(out, paste0(paste(.cgu_header(cg), collapse = ""), "\n"))
  if (empty) out <- c(out, "\n")
  for (routine in routines) out <- c(out, .cgu_get_interface(cg, routine))
  if (empty) out <- c(out, "\n")
  cat(paste(out, collapse = ""), file = f, sep = "")
}
.cgu_get_interface <- function(cg, routine)
  paste(c("interface\n", .cgu_opening(cg, routine), .cgu_declare_arguments(cg, routine), .cgu_ending(cg, routine), "end interface\n"), collapse = "")
.cgu_dump_m <- function(cg, routines, f, prefix, header = TRUE, empty = TRUE, inline = TRUE) {
  lines <- character()
  for (i in seq_along(routines)) {
    routine <- routines[[i]]
    if (i > 1L && empty) lines <- c(lines, "\n")
    lines <- c(lines, .cgu_opening(cg, routine))
    if (i == 1L) {
      if (routine$name != prefix) .cg_err("ValueError", "Octave function name should match prefix")
      if (header) lines <- c(lines, paste0("%", toupper(prefix), "  Autogenerated by symcompR\n"), paste(.cgu_header(cg), collapse = ""))
    }
    lines <- c(lines, .cgu_declare_arguments(cg, routine), .cgu_declare_globals(cg, routine), .cgu_declare_locals(cg, routine))
    if (empty) lines <- c(lines, "\n")
    lines <- c(lines, .cgu_call_printer(cg, routine))
    if (empty) lines <- c(lines, "\n")
    lines <- c(lines, .cgu_ending(cg, routine))
  }
  .cgu_write_out(f, .cgu_indent(cg, paste(lines, collapse = "")))
}
.cgu_dump_fns <- function(cg) switch(.cgu_lang(cg),
  c = list(list(.cgu_dump_code, "c"), list(.cgu_dump_c_h, "h")),
  f = list(list(.cgu_dump_f95, "f90"), list(.cgu_dump_f_h, "h")),
  julia = list(list(.cgu_dump_code, "jl")), octave = list(list(.cgu_dump_m, "m")), rust = list(list(.cgu_dump_code, "rs")))

.cgu_routine <- function(cg, name, expr, argument_sequence = NULL, global_vars = NULL) {
  num <- function(e) if (is.numeric(e) && length(e) == 1L) as_sym(e) else e
  expr <- if (.cgu_seq(expr)) lapply(expr, num) else num(expr)
  lang <- .cgu_lang(cg)
  if (lang %in% c("c", "f")) .cgu_routine_generic(cg, name, expr, argument_sequence, global_vars)
  else .cgu_routine_results(cg, name, expr, argument_sequence, global_vars, lang)
}
.cgu_write <- function(cg, routines, prefix, to_files = FALSE, header = TRUE, empty = TRUE) {
  res <- list()
  for (df in .cgu_dump_fns(cg)) {
    filename <- sprintf("%s.%s", prefix, df[[2L]])
    if (to_files) {
      con <- file(filename, "w")
      df[[1L]](cg, routines, con, prefix, header, empty)
      close(con)
    } else {
      con <- rawConnection(raw(0), "w")
      df[[1L]](cg, routines, con, prefix, header, empty)
      s <- rawToChar(rawConnectionValue(con)); close(con)
      Encoding(s) <- "UTF-8"
      res[[length(res) + 1L]] <- structure(list(filename, s), class = "pytuple")
    }
  }
  if (to_files) invisible(NULL) else res
}
.cgu_cg_attr <- function(cg, name) {
  switch(name,
         routine = function(name, expr, argument_sequence = NULL, global_vars = NULL) .cgu_routine(cg, name, expr, argument_sequence, global_vars),
         write = function(routines, prefix, to_files = FALSE, header = TRUE, empty = TRUE) .cgu_write(cg, routines, prefix, to_files, header, empty),
         dump_code = , dump_c = , dump_jl = , dump_rs = function(routines, f, prefix, header = TRUE, empty = TRUE) .cgu_dump_code(cg, routines, f, prefix, header, empty),
         dump_f95 = function(routines, f, prefix, header = TRUE, empty = TRUE) .cgu_dump_f95(cg, routines, f, prefix, header, empty),
         dump_m = function(routines, f, prefix, header = TRUE, empty = TRUE, inline = TRUE) .cgu_dump_m(cg, routines, f, prefix, header, empty, inline),
         dump_h = function(routines, f, prefix, header = TRUE, empty = TRUE)
           if (.cgu_lang(cg) == "c") .cgu_dump_c_h(cg, routines, f, prefix, header, empty) else .cgu_dump_f_h(cg, routines, f, prefix, header, empty),
         get_prototype = function(routine) if (.cgu_lang(cg) == "rust") .cgu_rust_prototype(cg, routine) else .cgu_c_prototype(cg, routine),
         get_interface = function(routine) .cgu_get_interface(cg, routine),
         printer = cg$printer, project = cg$project, cse = cg$cse, preprocessor_statements = cg$preprocessor_statements,
         .cg_err("AttributeError", sprintf("'%s' object has no attribute '%s'", cg$cls, name)))
}
#' @export
print.cg_codegen <- function(x, ...) {
  cat("<", x$cls, ">\n", sep = "")
  invisible(x)
}

#' Code generators
#'
#' Code generators write routines as source files, as the classes of
#' `sympy.utilities.codegen`: `CCodeGen()`, `C89CodeGen()`, `C99CodeGen()`,
#' `FCodeGen()`, `JuliaCodeGen()`, `OctaveCodeGen()` and `RustCodeGen()`.
#' Their methods (`routine`, `write`, `dump_c`, `dump_h`, `dump_f95`,
#' `dump_jl`, `dump_m`, `dump_rs`, `get_prototype`, ...) are read with
#' [ast_attr()]; the dump methods write to a connection.
#' `get_code_generator()` returns the generator for a language, and
#' `codegen()` writes the source files of named expressions.
#'
#' @param project The project name.
#' @param printer A code printer object, or `NULL` for the default.
#' @param preprocessor_statements C preprocessor lines.
#' @param cse Whether to eliminate common subexpressions.
#' @param language "C", "C89", "C99", "F95", "JULIA", "OCTAVE" or "RUST".
#' @param standard The C standard, for language "C".
#' @param name_expr A `list(name, expr)` pair, or a list of them.
#' @param prefix The file name prefix.
#' @param to_files Whether to write files instead of returning them.
#' @param header,empty Whether to write the header comment and empty lines.
#' @param argument_sequence,global_vars As in [make_routine()].
#' @param code_gen A code generator, instead of `language`.
#' @return A code generator; `codegen()` returns a list of
#'   `(file name, contents)` pairs.
#' @examples
#' x <- Symbol("x"); y <- Symbol("y")
#' res <- codegen(list("f", x * y), "C99", header = FALSE, empty = FALSE)
#' cat(res[[1]][[2]])
#' @name code-generators
NULL
#' @rdname code-generators
#' @export
CCodeGen <- function(project = "project", printer = NULL, preprocessor_statements = NULL, cse = FALSE)
  .cgu_new("CCodeGen", project, printer, cse, preprocessor_statements)
#' @rdname code-generators
#' @export
C89CodeGen <- function(project = "project", printer = NULL, preprocessor_statements = NULL, cse = FALSE)
  .cgu_new("C89CodeGen", project, printer, cse, preprocessor_statements)
#' @rdname code-generators
#' @export
C99CodeGen <- function(project = "project", printer = NULL, preprocessor_statements = NULL, cse = FALSE)
  .cgu_new("C99CodeGen", project, printer, cse, preprocessor_statements)
#' @rdname code-generators
#' @export
FCodeGen <- function(project = "project", printer = NULL) .cgu_new("FCodeGen", project, printer)
#' @rdname code-generators
#' @export
JuliaCodeGen <- function(project = "project", printer = NULL) .cgu_new("JuliaCodeGen", project, printer)
#' @rdname code-generators
#' @export
OctaveCodeGen <- function(project = "project", printer = NULL) .cgu_new("OctaveCodeGen", project, printer)
#' @rdname code-generators
#' @export
RustCodeGen <- function(project = "project", printer = NULL) .cgu_new("RustCodeGen", project, printer)

#' @rdname code-generators
#' @export
get_code_generator <- function(language, project = NULL, standard = NULL, printer = NULL) {
  if (language == "C" && !is.null(standard)) {
    if (tolower(standard) == "c89") language <- "C89" else if (tolower(standard) == "c99") language <- "C99"
  }
  cls <- c(C = "CCodeGen", C89 = "C89CodeGen", C99 = "C99CodeGen", F95 = "FCodeGen", JULIA = "JuliaCodeGen",
           OCTAVE = "OctaveCodeGen", RUST = "RustCodeGen")[toupper(language)]
  if (is.na(cls)) .cg_err("ValueError", sprintf("Language '%s' is not supported.", language))
  .cgu_new(cls, if (is.null(project)) "project" else project, printer)
}

#' @rdname code-generators
#' @export
codegen <- function(name_expr, language = NULL, prefix = NULL, project = "project", to_files = FALSE, header = TRUE,
                    empty = TRUE, argument_sequence = NULL, global_vars = NULL, standard = NULL, code_gen = NULL,
                    printer = NULL) {
  if (is.null(language)) {
    if (is.null(code_gen)) .cg_err("ValueError", "Need either language or code_gen")
  } else {
    if (!is.null(code_gen)) .cg_err("ValueError", "You cannot specify both language and code_gen.")
    code_gen <- get_code_generator(language, project, standard, printer)
  }
  if (is.character(name_expr[[1L]])) name_expr <- list(name_expr)
  if (is.null(prefix)) prefix <- name_expr[[1L]][[1L]]
  routines <- lapply(name_expr, function(ne) .cgu_routine(code_gen, ne[[1L]], ne[[2L]], argument_sequence, global_vars))
  .cgu_write(code_gen, routines, prefix, to_files, header, empty)
}

#' @rdname Routine
#' @export
make_routine <- function(name, expr, argument_sequence = NULL, global_vars = NULL, language = "F95")
  .cgu_routine(get_code_generator(language), name, expr, argument_sequence, global_vars)
