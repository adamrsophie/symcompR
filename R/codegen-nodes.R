# Language specific codegen nodes, after sympy/codegen/cnodes.py,
# cxxnodes.py, fnodes.py, pynodes.py, abstract_nodes.py, matrix_nodes.py,
# numpy_nodes.py and scipy_nodes.py.

# Nodes that derive from Basic only (not from Token)
.cg_basic_new <- function(type, args) new_sym(type, lapply(args, .cg_sympify))

# ---------------------------------------------------------------------------
# abstract_nodes.py and pynodes.py

.cg_def("List", token = FALSE, tail = c("Tuple", "Basic"), module = "sympy.codegen.abstract_nodes",
        str = function(e, repr, ...) paste0("[", paste(vapply(e$args, .cg_p, "", repr = repr), collapse = ", "), "]"))
cg_List <- function(...) .cg_basic_new("List", list(...))

#' Codegen nodes for C, C++, Python and lists
#'
#' Nodes after SymPy's `sympy.codegen.cnodes`, `cxxnodes`, `pynodes` and
#' `abstract_nodes`. `alignof()` and `sizeof()` are calls to the C
#' operators; `CommaOperator()` is the C comma operator; `Label()` and
#' `goto()` are a labeled statement and a jump; `PreIncrement()`,
#' `PostIncrement()`, `PreDecrement()` and `PostDecrement()` are the C
#' increment and decrement operators; `struct()` and `c_union()` declare a
#' struct and a union (SymPy's `union`); `void` is the void type and
#' `restrict`, `volatile` and `static` are attributes. `using()` is a C++
#' using declaration and `constexpr` an attribute. `List()` is a list, which
#' prints with square brackets, and `NumExprEvaluate()` a call to
#' `numexpr.evaluate`.
#'
#' @param arg A type name or an expression.
#' @param ... Arguments or fields of the node.
#' @return A node.
#' @examples
#' x <- Symbol("x")
#' ccode(sizeof("unsigned int"))
#' ccode(CommaOperator(PreIncrement(x), 2 * x))
#' List(1, 2, 3)
#' @name codegen-cnodes
NULL

#' @rdname codegen-cnodes
#' @export
List <- function(...) .cg_basic_new("List", list(...))

.cg_def("NumExprEvaluate", "Token", fields = "expr", module = "sympy.codegen.pynodes")
#' @rdname codegen-cnodes
#' @export
NumExprEvaluate <- function(...) .cg_new("NumExprEvaluate", list(...))

# ---------------------------------------------------------------------------
# cnodes.py

#' @rdname codegen-cnodes
#' @export
void <- Type("void")
#' @rdname codegen-cnodes
#' @export
restrict <- Attribute("restrict")
#' @rdname codegen-cnodes
#' @export
volatile <- Attribute("volatile")
#' @rdname codegen-cnodes
#' @export
static <- Attribute("static")

#' @rdname codegen-cnodes
#' @export
alignof <- function(arg) FunctionCall("alignof", list(if (is.character(arg)) String(arg) else arg))
#' @rdname codegen-cnodes
#' @export
sizeof <- function(arg) FunctionCall("sizeof", list(if (is.character(arg)) String(arg) else arg))

.cg_def("CommaOperator", token = FALSE, tail = "Basic", module = "sympy.codegen.cnodes")
#' @rdname codegen-cnodes
#' @export
CommaOperator <- function(...) .cg_basic_new("CommaOperator", list(...))

.cg_body_construct <- function(itr) if (is_sym(itr) && itr$type == "CodeBlock") itr else do.call(CodeBlock, .cg_iter(itr))

.cg_def("Label", "Node", fields = c("name", "body"), defaults = list(body = none), module = "sympy.codegen.cnodes",
        construct = list(name = function(a) String(a), body = .cg_body_construct))
#' @rdname codegen-cnodes
#' @export
Label <- function(...) .cg_new("Label", list(...))

.cg_def("goto", "Token", fields = "label", module = "sympy.codegen.cnodes", construct = list(label = function(a) Label(a)))
#' @rdname codegen-cnodes
#' @export
goto <- function(...) .cg_new("goto", list(...))

for (.cg_n in c("PreDecrement", "PostDecrement", "PreIncrement", "PostIncrement"))
  .cg_def(.cg_n, token = FALSE, tail = "Basic", module = "sympy.codegen.cnodes")
rm(.cg_n)
#' @rdname codegen-cnodes
#' @export
PreDecrement <- function(...) .cg_basic_new("PreDecrement", list(...))
#' @rdname codegen-cnodes
#' @export
PostDecrement <- function(...) .cg_basic_new("PostDecrement", list(...))
#' @rdname codegen-cnodes
#' @export
PreIncrement <- function(...) .cg_basic_new("PreIncrement", list(...))
#' @rdname codegen-cnodes
#' @export
PostIncrement <- function(...) .cg_basic_new("PostIncrement", list(...))

.cg_def("struct", "Node", fields = c("name", "declarations"), defaults = list(name = none), module = "sympy.codegen.cnodes",
        construct = list(name = function(a) String(a),
                         declarations = function(a) .cg_tuple(lapply(.cg_iter(a), function(d) Declaration(d)))))
.cg_def("union", "struct", ctor = "c_union")
#' @rdname codegen-cnodes
#' @export
struct <- function(...) .cg_new("struct", list(...))
#' @rdname codegen-cnodes
#' @export
c_union <- function(...) .cg_new("union", list(...))

# ---------------------------------------------------------------------------
# cxxnodes.py

.cg_def("using", "Token", fields = c("type", "alias"), defaults = list(alias = none), module = "sympy.codegen.cxxnodes",
        construct = list(type = function(a) Type(a), alias = function(a) String(a)))
#' @rdname codegen-cnodes
#' @export
using <- function(...) .cg_new("using", list(...))
#' @rdname codegen-cnodes
#' @export
constexpr <- Attribute("constexpr")

# ---------------------------------------------------------------------------
# matrix_nodes.py

.cg_def("MatrixSolve", "Token", fields = c("matrix", "vector"), module = "sympy.codegen.matrix_nodes",
        tail = c("CodegenAST", "MatrixExpr", "Expr", "Basic"),
        construct = list(matrix = function(a) .cg_sympify(a), vector = function(a) .cg_sympify(a)))

#' Solve a linear system in generated code
#'
#' `MatrixSolve(A, b)` stands for the solution `x` of `A x = b`, as SymPy's
#' `sympy.codegen.matrix_nodes.MatrixSolve`; the NumPy printer prints it as
#' `numpy.linalg.solve(A, b)`. It differentiates as
#' `MatrixSolve(A, b' - A' MatrixSolve(A, b))`.
#'
#' @param ... The matrix and the vector.
#' @return A node.
#' @examples
#' A <- MatrixSymbol("A", 3, 3); b <- MatrixSymbol("b", 3, 1)
#' printer_doprint(NumPyPrinter(), MatrixSolve(A, b))
#' @export
MatrixSolve <- function(...) .cg_new("MatrixSolve", list(...))

# MatrixSolve._eval_derivative
register_diff_hook("MatrixSolve", function(e, x) {
  v <- .cg_vals(e)
  A <- .unwrap_explicit(v$matrix); b <- .unwrap_explicit(v$vector)
  MatrixSolve(v$matrix, MatAdd(diff(b, x), MatMul(S_NegativeOne, diff(A, x), MatrixSolve(v$matrix, v$vector),
                                                   evaluate = TRUE)))
})

# ---------------------------------------------------------------------------
# fnodes.py

#' Fortran codegen nodes
#'
#' Nodes after SymPy's `sympy.codegen.fnodes`: `Program()`, `Module()`,
#' `Subroutine()` and `SubroutineCall()` build program units; `Do()` is a do
#' loop, `ArrayConstructor()` an array constructor and `ImpliedDoLoop()` an
#' implied do loop; `GoTo()` is a computed goto and `FortranReturn()` a
#' return statement; `fortran_use()` and `use_rename()` are use statements
#' (SymPy's `use`); `Extent()` is an array extent, `fortran_dimension()` the
#' dimension attribute (SymPy's `dimension`) and `fortran_array()` declares
#' an array (SymPy's `array`); `allocated()`, `lbound()`, `ubound()`,
#' `fortran_shape()`, `size()` and `fortran_reshape()` call the intrinsic functions
#' (SymPy's `shape`); `bind_C()` is the bind(C) attribute; `isign()`,
#' `dsign()`, `cmplx()`, `kind()` and `fortran_merge()` (SymPy's `merge`)
#' are intrinsic functions; `literal_sp()` and `literal_dp()` are single and
#' double precision literals; `sum_()` and `product_()` are the array
#' reductions.
#'
#' @param ... Arguments or fields of the node.
#' @param symbol,source,array A symbol, an array or an expression.
#' @param dim,kind,shape,pad,order,intent,attrs,value,type,name Optional
#'   arguments of the Fortran functions and attributes.
#' @param x,y,z Expressions.
#' @return A node or expression.
#' @examples
#' x <- Symbol("x"); i <- Symbol("i", integer = TRUE); n <- Symbol("n", integer = TRUE)
#' fcode(Do(list(AddAugmentedAssignment(x, i)), i, 1, n), source_format = "free")
#' fcode(isign(1, x), source_format = "free")
#' @name codegen-fnodes
NULL

.cg_fmod <- "sympy.codegen.fnodes"

#' @rdname codegen-fnodes
#' @export
pure <- Attribute("pure")
#' @rdname codegen-fnodes
#' @export
elemental <- Attribute("elemental")
#' @rdname codegen-fnodes
#' @export
intent_in <- Attribute("intent_in")
#' @rdname codegen-fnodes
#' @export
intent_out <- Attribute("intent_out")
#' @rdname codegen-fnodes
#' @export
intent_inout <- Attribute("intent_inout")
#' @rdname codegen-fnodes
#' @export
allocatable <- Attribute("allocatable")

.cg_def("Program", "Token", fields = c("name", "body"), module = .cg_fmod,
        construct = list(name = function(a) String(a), body = function(b) do.call(CodeBlock, .cg_iter(b))))
#' @rdname codegen-fnodes
#' @export
Program <- function(...) .cg_new("Program", list(...))

.cg_def("use_rename", "Token", fields = c("local", "original"), module = .cg_fmod,
        construct = list(local = function(a) String(a), original = function(a) String(a)))
#' @rdname codegen-fnodes
#' @export
use_rename <- function(...) .cg_new("use_rename", list(...))

.cg_fname_arg <- function(arg) {
  if (is_sym(arg) && (.cg_isa(arg, "Variable") || .cg_isa(arg, "FunctionPrototype") || arg$type %in% c("Symbol", "Dummy")))
    return(if (arg$type %in% c("Symbol", "Dummy")) arg$name
           else if (!is.null(.cg_vals(arg)$name)) .cg_vals(arg)$name else .cg_vals(arg)$symbol)
  if (is_sym(arg) && !is.null(.cg_cls(arg)) && "name" %in% .cg_cls(arg)$fields) return(.cg_vals(arg)$name)
  String(arg)
}

.cg_def("use", "Token", fields = c("namespace", "rename", "only"), defaults = list(rename = none, only = none),
        module = .cg_fmod, ctor = "fortran_use",
        construct = list(namespace = .cg_fname_arg,
                         rename = function(args) .cg_tuple(lapply(.cg_iter(args), function(a)
                           if (.cg_isa(a, "use_rename")) a else do.call(use_rename, .cg_iter(a)))),
                         only = function(args) .cg_tuple(lapply(.cg_iter(args), function(a)
                           if (.cg_isa(a, "use_rename")) a else .cg_fname_arg(a)))))
#' @rdname codegen-fnodes
#' @export
fortran_use <- function(...) .cg_new("use", list(...))

.cg_def("Module", "Token", fields = c("name", "declarations", "definitions"),
        defaults = list(declarations = new_sym("Tuple")), module = .cg_fmod,
        construct = list(name = function(a) String(a),
                         declarations = function(args) do.call(CodeBlock, lapply(.cg_iter(args), function(a) if (is.character(a)) String(a) else a)),
                         definitions = function(a) do.call(CodeBlock, .cg_iter(a))))
#' @rdname codegen-fnodes
#' @export
Module <- function(...) .cg_new("Module", list(...))

.cg_def("Subroutine", "Node", fields = c("name", "parameters", "body", "attrs"), module = .cg_fmod,
        construct = list(name = function(a) String(a),
                         parameters = function(params) .cg_tuple(lapply(.cg_iter(params), function(q) variable_deduced(q))),
                         body = .cg_body_construct))
#' @rdname codegen-fnodes
#' @export
Subroutine <- function(...) .cg_new("Subroutine", list(...))

.cg_def("SubroutineCall", "Token", fields = c("name", "subroutine_args"), module = .cg_fmod,
        construct = list(name = .cg_fname_arg, subroutine_args = function(a) .cg_mk_tuple(a)))
#' @rdname codegen-fnodes
#' @export
SubroutineCall <- function(...) .cg_new("SubroutineCall", list(...))

.cg_def("Do", "Token", fields = c("body", "counter", "first", "last", "step", "concurrent"),
        defaults = list(step = Integer(1), concurrent = S_false), module = .cg_fmod,
        construct = list(body = function(b) do.call(CodeBlock, .cg_iter(b)), counter = function(a) .cg_sympify(a),
                         first = function(a) .cg_sympify(a), last = function(a) .cg_sympify(a),
                         step = function(a) .cg_sympify(a),
                         concurrent = function(a) if (.cgp_truth(a)) S_true else S_false))
#' @rdname codegen-fnodes
#' @export
Do <- function(...) .cg_new("Do", list(...))

.cg_def("ArrayConstructor", "Token", fields = "elements", module = .cg_fmod,
        construct = list(elements = function(a) .cg_mk_tuple(a)))
#' @rdname codegen-fnodes
#' @export
ArrayConstructor <- function(...) .cg_new("ArrayConstructor", list(...))

.cg_def("ImpliedDoLoop", "Token", fields = c("expr", "counter", "first", "last", "step"),
        defaults = list(step = Integer(1)), module = .cg_fmod,
        construct = list(expr = function(a) .cg_sympify(a), counter = function(a) .cg_sympify(a),
                         first = function(a) .cg_sympify(a), last = function(a) .cg_sympify(a),
                         step = function(a) .cg_sympify(a)))
#' @rdname codegen-fnodes
#' @export
ImpliedDoLoop <- function(...) .cg_new("ImpliedDoLoop", list(...))

.cg_def("Extent", token = FALSE, tail = "Basic", module = .cg_fmod,
        str = function(e, repr, ...) {
          if (repr) return(.cg_basic_print(e, TRUE))
          if (!length(e$args)) return(":")
          paste(vapply(e$args, sstr, ""), collapse = ":")
        })
#' @rdname codegen-fnodes
#' @export
Extent <- function(...) {
  args <- list(...)
  if (length(args) == 2L) return(new_sym("Extent", lapply(args, .cg_sympify)))
  if (!length(args) || (length(args) == 1L && (is.null(args[[1L]]) || identical(args[[1L]], ":"))))
    return(new_sym("Extent"))
  .cg_err("ValueError", "Expected 0 or 2 args (or one argument == None or ':')")
}
#' @rdname codegen-fnodes
#' @export
assumed_extent <- new_sym("Extent")

#' @rdname codegen-fnodes
#' @export
fortran_dimension <- function(...) {
  args <- list(...)
  if (length(args) > 7L) .cg_err("ValueError", "Fortran only supports up to 7 dimensional arrays")
  params <- lapply(args, function(a) {
    if (is_sym(a) && a$type == "Extent") return(a)
    if (is.character(a)) return(if (a == ":") Extent() else String(a))
    if (.cg_is_iterable(a)) return(do.call(Extent, .cg_iter(a)))
    .cg_sympify(a)
  })
  if (!length(args)) .cg_err("ValueError", "Need at least one dimension")
  Attribute("dimension", params)
}
#' @rdname codegen-fnodes
#' @export
assumed_size <- fortran_dimension("*")

#' @rdname codegen-fnodes
#' @export
fortran_array <- function(symbol, dim, intent = NULL, attrs = list(), value = NULL, type = NULL) {
  if (is_sym(dim) && .cg_isa(dim, "Attribute")) {
    if (.cg_text(dim$args[[1L]]) != "dimension")
      .cg_err("ValueError", sprintf("Got an unexpected Attribute argument as dim: %s", sstr(dim)))
  } else dim <- do.call(fortran_dimension, .cg_iter(dim))
  attrs <- c(.cg_iter(attrs), list(dim))
  if (!is.null(intent)) {
    known <- list(intent_in, intent_out, intent_inout)
    if (!any(vapply(known, function(k) is_sym(intent) && identical(k$key, intent$key), TRUE)))
      intent <- list(`in` = intent_in, out = intent_out, inout = intent_inout)[[intent]]
    attrs <- c(attrs, list(intent))
  }
  if (is.null(type)) return(variable_deduced(symbol, value = value, attrs = attrs))
  Variable(symbol, type, value = value, attrs = attrs)
}

.cg_printable <- function(arg) if (is.character(arg)) String(arg) else .cg_sympify(arg)

#' @rdname codegen-fnodes
#' @export
allocated <- function(array) FunctionCall("allocated", list(.cg_printable(array)))
#' @rdname codegen-fnodes
#' @export
lbound <- function(array, dim = NULL, kind = NULL)
  FunctionCall("lbound", c(list(.cg_printable(array)), if (!is.null(dim)) list(.cg_printable(dim)),
                           if (!is.null(kind)) list(.cg_printable(kind))))
#' @rdname codegen-fnodes
#' @export
ubound <- function(array, dim = NULL, kind = NULL)
  FunctionCall("ubound", c(list(.cg_printable(array)), if (!is.null(dim)) list(.cg_printable(dim)),
                           if (!is.null(kind)) list(.cg_printable(kind))))
#' @rdname codegen-fnodes
#' @export
fortran_shape <- function(source, kind = NULL)
  FunctionCall("shape", c(list(.cg_printable(source)), if (!is.null(kind)) list(.cg_printable(kind))))
#' @rdname codegen-fnodes
#' @export
size <- function(array, dim = NULL, kind = NULL)
  FunctionCall("size", c(list(.cg_printable(array)), if (!is.null(dim)) list(.cg_printable(dim)),
                         if (!is.null(kind)) list(.cg_printable(kind))))
#' @rdname codegen-fnodes
#' @export
fortran_reshape <- function(source, shape, pad = NULL, order = NULL) {
  nested <- function(arg) {
    if (is.list(arg) && !is_sym(arg)) return(ArrayConstructor(lapply(arg, nested)))
    .cg_printable(arg)
  }
  arr <- function(arg) {
    if (is.list(arg) && !is_sym(arg)) return(ArrayConstructor(arg))
    .cg_printable(arg)
  }
  kw <- list()
  if (!is.null(pad)) kw$pad <- nested(pad)
  if (!is.null(order)) kw$order <- arr(order)
  KeywordFunctionCall("reshape", list(nested(source), arr(shape)), kw)
}
#' @rdname codegen-fnodes
#' @export
bind_C <- function(name = NULL) Attribute("bind_C", if (!is.null(name)) list(String(name)) else list())

.cg_def("GoTo", "Token", fields = c("labels", "expr"), defaults = list(expr = none), module = .cg_fmod,
        construct = list(labels = function(a) .cg_mk_tuple(a), expr = function(a) .cg_sympify(a)))
#' @rdname codegen-fnodes
#' @export
GoTo <- function(...) .cg_new("GoTo", list(...))

.cg_def("FortranReturn", "Token", fields = "return_value", defaults = list(return_value = none), module = .cg_fmod,
        construct = list(return_value = function(a) .cg_sympify(a)))
#' @rdname codegen-fnodes
#' @export
FortranReturn <- function(...) .cg_new("FortranReturn", list(...))

# FFunction and F95Function: functions printed by name, with a required
# Fortran standard
.cg_ffunctions <- c(isign = 77, dsign = 77, cmplx = 77, kind = 77, merge = 95)
for (.cg_n in names(.cg_ffunctions))
  def_function(.cg_n, nargs = c(isign = 2L, dsign = 2L, cmplx = 2L, kind = 1L, merge = 3L)[[.cg_n]])
rm(.cg_n)
#' @rdname codegen-fnodes
#' @export
isign <- function(x, y) fcall("isign", list(as_sym(x), as_sym(y)))
#' @rdname codegen-fnodes
#' @export
dsign <- function(x, y) fcall("dsign", list(as_sym(x), as_sym(y)))
#' @rdname codegen-fnodes
#' @export
cmplx <- function(x, y) fcall("cmplx", list(as_sym(x), as_sym(y)))
#' @rdname codegen-fnodes
#' @export
kind <- function(x) fcall("kind", list(as_sym(x)))
#' @rdname codegen-fnodes
#' @export
fortran_merge <- function(x, y, z) fcall("merge", list(as_sym(x), as_sym(y), as_sym(z)))

# _literal(Float): a Float that prints with a Fortran exponent letter
.cg_literal <- function(cls, x) {
  f <- Float(x)
  new_sym("Float", value = f$value, extra = c(f$extra, list(literal = cls)), key = paste0(f$key, "#", cls))
}
#' @rdname codegen-fnodes
#' @export
literal_sp <- function(x) .cg_literal("literal_sp", x)
#' @rdname codegen-fnodes
#' @export
literal_dp <- function(x) .cg_literal("literal_dp", x)

.cg_def("sum_", "Token", fields = c("array", "dim", "mask"), defaults = list(dim = none, mask = none),
        module = .cg_fmod, expr = TRUE,
        construct = list(array = function(a) .cg_sympify(a), dim = function(a) .cg_sympify(a)))
.cg_def("product_", "sum_")
#' @rdname codegen-fnodes
#' @export
sum_ <- function(...) .cg_new("sum_", list(...))
#' @rdname codegen-fnodes
#' @export
product_ <- function(...) .cg_new("product_", list(...))

# ---------------------------------------------------------------------------
# numpy_nodes.py and scipy_nodes.py

.cg_two_log <- function(x) .div(fcall("log", list(x)), fcall("log", list(S_Two)))

# _eval_simplify of logaddexp, cosm1 and powm1: the evaluated form unless it
# equals the unevaluated one (`uneval`, built with new_sym). For logaddexp and
# logaddexp2 the sums inside the logarithms are compared, so that no
# evaluation runs on the unevaluated nodes.
def_function("logaddexp", nargs = 2L,
             diff = function(args, i) {
               wrt <- args[[i]]; other <- args[[3L - i]]
               .div(S_One, .add(list(S_One, fcall("exp", list(.sub(other, wrt))))))
             },
             simplify = function(e, sub) {
               ab <- lapply(e$args, sub)
               inner <- .add(list(fcall("exp", list(ab[[1L]])), fcall("exp", list(ab[[2L]]))))
               ex <- lapply(ab, function(z) new_sym("Function", list(z), name = "exp"))
               if (!identical(inner$key, new_sym("Add", ex)$key)) fcall("log", list(inner)) else logaddexp(ab[[1L]], ab[[2L]])
             },
             evalf = function(e) rewrite(e, "log"))
def_function("logaddexp2", nargs = 2L,
             diff = function(args, i) {
               wrt <- args[[i]]; other <- args[[3L - i]]
               .div(S_One, .add(list(S_One, .pow(S_Two, .sub(other, wrt)))))
             },
             simplify = function(e, sub) {
               ab <- lapply(e$args, function(z) factor(sub(z)))
               inner <- .add(list(.pow(S_Two, ab[[1L]]), .pow(S_Two, ab[[2L]])))
               p2 <- lapply(ab, function(z) new_sym("Pow", list(S_Two, z)))
               if (!identical(inner$key, new_sym("Add", p2)$key)) .cg_two_log(inner) else logaddexp2(ab[[1L]], ab[[2L]])
             },
             evalf = function(e) rewrite(e, "log"))
.cg_add_rewrite("log", "logaddexp", function(a, b) fcall("log", list(.add(list(fcall("exp", list(a)), fcall("exp", list(b)))))))
.cg_add_rewrite("log", "logaddexp2", function(a, b) .cg_two_log(.add(list(.pow(S_Two, a), .pow(S_Two, b)))))

.cg_sorted_new <- function(name, args) {
  args <- lapply(args, as_sym)
  fcall(name, .key_sort(lapply(args, .sort_key), args))
}

#' NumPy and SciPy codegen functions
#'
#' Functions after SymPy's `sympy.codegen.numpy_nodes` and `scipy_nodes`:
#' `logaddexp(a, b)` is `log(exp(a) + exp(b))`, `logaddexp2(a, b)` is
#' `log2(2^a + 2^b)`, `amin()` and `amax()` are array reductions,
#' `np_maximum()` and `np_minimum()` are the elementwise maximum and minimum
#' (SymPy's `maximum` and `minimum`), `cosm1(x)` is `cos(x) - 1` and
#' `powm1(x, y)` is `x^y - 1`. The Python printers print them as the NumPy
#' and SciPy calls; [rewrite()] turns them into elementary functions.
#'
#' @param a,b,x,y Expressions.
#' @param ... Arguments.
#' @return An expression or node.
#' @examples
#' x <- Symbol("x"); y <- Symbol("y")
#' logaddexp(x, y)
#' diff(logaddexp(x, y), x)
#' rewrite(cosm1(x), "cos")
#' @name codegen-numpy-nodes
NULL

#' @rdname codegen-numpy-nodes
#' @export
logaddexp <- function(a, b) .cg_sorted_new("logaddexp", list(a, b))
#' @rdname codegen-numpy-nodes
#' @export
logaddexp2 <- function(a, b) .cg_sorted_new("logaddexp2", list(a, b))

.cg_def("amin", "Token", fields = c("array", "axis"), defaults = list(axis = none), module = "sympy.codegen.numpy_nodes",
        construct = list(axis = function(a) .cg_sympify(a)))
.cg_def("amax", "amin")
#' @rdname codegen-numpy-nodes
#' @export
amin <- function(...) .cg_new("amin", list(...))
#' @rdname codegen-numpy-nodes
#' @export
amax <- function(...) .cg_new("amax", list(...))

def_function("maximum", nargs = NULL)
def_function("minimum", nargs = NULL)
.cg_add_rewrite("Max", "maximum", function(...) do.call(Max, list(...)))
.cg_add_rewrite("Min", "minimum", function(...) do.call(Min, list(...)))
#' @rdname codegen-numpy-nodes
#' @export
np_maximum <- function(...) fcall("maximum", lapply(list(...), as_sym))
#' @rdname codegen-numpy-nodes
#' @export
np_minimum <- function(...) fcall("minimum", lapply(list(...), as_sym))

def_function("cosm1", nargs = 1L, diff = function(args, i) .neg(fcall("sin", args)),
             simplify = function(e, sub) {
               x <- e$args[[1L]]
               candidate <- .add(list(fcall("cos", list(sub(x))), S_NegativeOne))
               uneval <- new_sym("Add", list(new_sym("Function", list(x), name = "cos"), S_NegativeOne))
               if (!identical(candidate$key, uneval$key)) candidate else cosm1(x)
             })
.cg_add_rewrite("cos", "cosm1", function(a) .add(list(fcall("cos", list(a)), S_NegativeOne)))
def_function("powm1", nargs = 2L,
             diff = function(args, i) {
               if (i == 1L) return(.div(.mul(list(.pow(args[[1L]], args[[2L]]), args[[2L]])), args[[1L]]))
               .mul(list(fcall("log", list(args[[1L]])), .pow(args[[1L]], args[[2L]])))
             },
             simplify = function(e, sub) {
               x <- e$args[[1L]]; y <- e$args[[2L]]
               candidate <- .add(list(.pow(sub(x), sub(y)), S_NegativeOne))
               uneval <- new_sym("Add", list(new_sym("Pow", list(x, y)), S_NegativeOne))
               if (!identical(candidate$key, uneval$key)) candidate else powm1(x, y)
             })
.cg_add_rewrite("Pow", "powm1", function(a, b) .add(list(.pow(a, b), S_NegativeOne)))
#' @rdname codegen-numpy-nodes
#' @export
cosm1 <- function(x) fcall("cosm1", list(as_sym(x)))
#' @rdname codegen-numpy-nodes
#' @export
powm1 <- function(x, y) fcall("powm1", list(as_sym(x), as_sym(y)))
