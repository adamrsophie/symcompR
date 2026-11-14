# Function objects used as SymPy classes: has(e, sin) and
# sym_replace(e, sin, cos) take a function type as Basic.has and
# Basic.replace do.

.fs_base_names <- c(abs = "Abs", sign = "sign", floor = "floor", ceiling = "ceiling",
                    exp = "exp", log = "log", cos = "cos", sin = "sin", tan = "tan",
                    acos = "acos", asin = "asin", atan = "atan", cosh = "cosh",
                    sinh = "sinh", tanh = "tanh", acosh = "acosh", asinh = "asinh",
                    atanh = "atanh", gamma = "gamma", lgamma = "loggamma",
                    digamma = "digamma", trigamma = "trigamma", factorial = "factorial")

.fs_core_types <- c("Add", "Mul", "Pow", "Symbol", "Integer", "Rational", "Float",
                    "Piecewise", "Derivative", "Integral", "Sum", "Product",
                    "And", "Or", "Not", "Xor", "Implies", "Equivalent", "Eq", "Ne",
                    "Lambda", "Tuple", "Interval", "FiniteSet", "Union", "Intersection")

# The SymPy class name an R function object stands for, or NULL when it is
# an ordinary function (a predicate for replace).
.fs_fn_class <- function(f) {
  if (!is.function(f)) return(NULL)
  if (inherits(f, "sym_function")) return(attr(f, "sym_name"))
  # only base functions and the package's own exports stand for classes
  env <- environment(f)
  if (!is.null(env)) {
    top <- topenv(env)
    if (!identical(top, .BaseNamespaceEnv) && !identical(top, baseenv()) && !identical(top, asNamespace("symcompR"))) return(NULL)
  }
  for (nm in names(.fs_base_names)) {
    if (identical(f, get(nm, envir = baseenv()))) return(.fs_base_names[[nm]])
  }
  ns <- asNamespace("symcompR")
  for (nm in getNamespaceExports(ns)) {
    if (!exists(nm, envir = ns, inherits = FALSE)) next
    if (!identical(f, get(nm, envir = ns))) next
    cls <- sub("^sym_", "", nm)
    if (cls == "Eq") cls <- "Equality"
    if (cls == "Ne") cls <- "Unequality"
    if (exists(cls, envir = .FUNS, inherits = FALSE) || cls %in% c(.fs_core_types, "Equality", "Unequality")) return(cls)
  }
  NULL
}

# Basic.has with a class: any node of that class.
.fs_has_head <- function(e, cls) {
  if (!is_sym(e)) return(FALSE)
  if (identical(.cg_head(e), cls)) return(TRUE)
  for (a in e$args) if (.fs_has_head(a, cls)) return(TRUE)
  FALSE
}

# sympify(str) for input outside R's grammar: parse_expr(s, local_dict =
# locals), NULL when that fails too
.fs_sympify_python <- function(s, locals, rational) {
  if (rational) return(NULL)
  ld <- if (length(locals) && !is.null(names(locals))) locals else NULL
  r <- tryCatch(parse_expr(s, local_dict = ld), error = function(e) NULL)
  if (is.null(r)) return(NULL)
  r
}

# Symbol._hashable_content: the name, then sorted(assumptions0.items()) as
# (name, bool) pairs; a pair is the string name + "\001" + 0/1, which sorts
# as the Python tuple does
.fs_assumption_items <- function(e) {
  hit <- .fs_hc_memo[[e$key]]
  if (!is.null(hit)) return(hit)
  fx <- e$extra$facts
  fx <- if (is.null(fx)) c(commutative = TRUE) else fx[!is.na(fx)]
  # hermitian and antihermitian are deduced here from real and imaginary,
  # which SymPy's fact rules do not do
  hg <- c("hermitian", "antihermitian")
  fx <- fx[!(names(fx) %in% setdiff(hg, names(e$extra$given)))]
  o <- e$extra$other
  if (length(o)) fx <- c(fx, o)
  if (!length(fx)) return(list())
  fx <- fx[order(names(fx), method = "radix")]
  r <- as.list(paste0(names(fx), "\001", ifelse(fx, "1", "0")))
  assign(e$key, r, envir = .fs_hc_memo)
  r
}
.fs_hc_memo <- new.env(parent = emptyenv())

# obj.method(*args, **kw) in parse_expr for the other methods: the exported
# function of the same name (or sym_<name>) called with obj first
.fs_pe_method <- function(obj, attr) {
  if (!(is_sym(obj) || inherits(obj, "symmatrix"))) return(NULL)
  if (!grepl("^[a-z][A-Za-z0-9_]*$", attr)) return(NULL)
  ns <- asNamespace("symcompR")
  ex <- getNamespaceExports(ns)
  if (attr == "inv") attr <- "inverse"
  cls <- if (is_sym(obj)) "sym" else "symmatrix"
  fn <- if (paste0("sym_", attr) %in% ex) get(paste0("sym_", attr), envir = ns)
        else if (attr %in% ex) get(attr, envir = ns)
        else if (!is.null(utils::getS3method(attr, cls, optional = TRUE, envir = ns))) utils::getS3method(attr, cls, envir = ns)
        else return(NULL)
  if (!is.function(fn)) return(NULL)
  .pe_callable(function(args, kw) do.call(fn, c(list(obj), lapply(args, .pe_to_r), lapply(kw, .pe_to_r))), attr)
}
