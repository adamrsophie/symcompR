# Parsing strings and R calls into expressions, after
# sympy/parsing/sympy_parser.py. Input is never evaluated as R code: the
# parse tree is walked and only whitelisted operators and registered
# function names are interpreted.


.CONST_NAMES <- list(pi = function() S_Pi, E = function() S_E, I = function() S_I,
                     oo = function() S_Infinity, `Inf` = function() S_Infinity,
                     zoo = function() S_ComplexInfinity, nan = function() S_NaN,
                     `NaN` = function() S_NaN,
                     EulerGamma = function() S_EulerGamma,
                     Catalan = function() S_Catalan,
                     GoldenRatio = function() S_GoldenRatio,
                     TribonacciConstant = function() S_TribonacciConstant,
                     True = function() S_true, False = function() S_false)

.NUM_PATTERN <- "(?<![A-Za-z0-9_.])((?:[0-9]+\\.?[0-9]*|\\.[0-9]+)(?:[eE][-+]?[0-9]+)?)(?![A-Za-z_])"

#' Parse a string into a symbolic expression
#'
#' `sympify()` reads a mathematical expression written with R or Python
#' operators (`^` or `**` for powers) and returns a `sym` expression. Number
#' literals keep their exact digits, so `"123456789012345678901234567890"`
#' and `"1/3"` stay exact, while decimals become [Float()] values. Names
#' become symbols except for known constants (`pi`, `E`, `I`, `oo`) and
#' functions; unknown function names become undefined functions. The input
#' is parsed, never evaluated as R code. `==` in a string means [Eq()].
#'
#' @param s A character string, or a value accepted by [as_sym()].
#' @param locals Optional named list mapping names to `sym` objects or R
#'   functions, used before the built-in names.
#' @param rational If `TRUE`, decimal literals become exact rationals.
#' @return A `sym` expression.
#' @examples
#' sympify("x**2 + 2*x + 1")
#' sympify("sin(x)^2 + cos(x)^2")
#' sympify("1/3 + 0.5")
#' @export
sympify <- function(s, locals = list(), rational = FALSE) {
  if (!is.character(s)) return(as_sym(s))
  if (length(s) != 1L) return(lapply(s, sympify, locals = locals, rational = rational))
  txt <- gsub("**", "^", s, fixed = TRUE)
  txt <- .mark_tuples(txt)
  txt <- gsub(perl = TRUE, .NUM_PATTERN, ".symnum('\\1')", txt)
  ex <- tryCatch(parse(text = txt, keep.source = FALSE),
                 error = function(e) stop(SympifyError(paste0("could not parse '", s, "': ", conditionMessage(e)))))
  if (length(ex) != 1L) stop(SympifyError(paste0("expected a single expression in '", s, "'")))
  sym_from_call(ex[[1L]], locals = locals, rational = rational, from_string = TRUE)
}

#' @rdname sympify
#' @param x An R call, formula, or name.
#' @param from_string Internal flag; `TRUE` when called from `sympify()`.
#' @export
sym_from_call <- function(x, locals = list(), rational = FALSE, from_string = FALSE) {
  rec <- function(z) sym_from_call(z, locals, rational, from_string)
  if (inherits(x, "formula")) x <- x[[length(x)]]
  if (is.numeric(x) || is.logical(x) || is.complex(x)) return(as_sym(x))
  if (is.character(x)) return(sympify(x, locals, rational))
  if (is.symbol(x)) {
    nm <- as.character(x)
    if (!is.null(locals[[nm]])) return(as_sym(locals[[nm]]))
    if (!is.null(.CONST_NAMES[[nm]])) return(.CONST_NAMES[[nm]]())
    return(Symbol(nm))
  }
  if (!is.call(x)) stop("unsupported input of class ", class(x)[1])
  op <- x[[1L]]
  if (!is.symbol(op)) stop("unsupported call head in ", deparse(x)[1])
  f <- as.character(op)
  args <- as.list(x)[-1L]
  if (f == ".symnum") {
    lit <- args[[1L]]
    if (grepl("[.eE]", lit)) {
      if (rational) return(make_rational(.decimal_to_bigq(lit)))
      return(Float(as.numeric(lit)))
    }
    return(make_rational(gmp::as.bigq(.str_bigz(lit))))
  }
  if (f == "(") return(rec(args[[1L]]))
  if (f == ".tuple") return(do.call(Tuple, lapply(args, rec)))
  if (f %in% c("+", "-") && length(args) == 1L) {
    a <- rec(args[[1L]])
    return(if (f == "-") .neg(a) else a)
  }
  bin <- c("+", "-", "*", "/", "^", "%%", "%/%", "==", "!=", "<", "<=", ">", ">=",
           "&", "|", "&&", "||")
  if (f %in% bin && length(args) == 2L) {
    a <- rec(args[[1L]])
    b <- rec(args[[2L]])
    return(switch(f,
                  "+" = .add(list(a, b)), "-" = .sub(a, b),
                  "*" = .mul(list(a, b)), "/" = .div(a, b), "^" = .pow(a, b),
                  "%%" = fcall("Mod", list(a, b)),
                  "%/%" = fcall("floor", list(.div(a, b))),
                  "==" = Eq(a, b), "!=" = Ne(a, b), "<" = Lt(a, b), "<=" = Le(a, b),
                  ">" = Gt(a, b), ">=" = Ge(a, b),
                  "&" = , "&&" = And(a, b), "|" = , "||" = Or(a, b)))
  }
  if (f == "!" && length(args) == 1L) return(Not(rec(args[[1L]])))
  if (!is.null(locals[[f]])) {
    fn <- locals[[f]]
    return(as_sym(do.call(fn, lapply(args, rec))))
  }
  fn <- .parse_functions[[f]]
  # a Python tuple argument of an R constructor is an R list
  untuple <- function(a) if (is_sym(a) && a$type == "Tuple") a$args else a
  if (!is.null(fn)) return(do.call(fn, lapply(lapply(args, rec), untuple)))
  rfn <- if (f %in% c("Piecewise", "Sum", "Product", "Integral", "integrate", "summation", "product", "Limit", "limit"))
    get0(f, envir = asNamespace("symcompR"), mode = "function") else NULL
  if (!is.null(rfn)) return(do.call(rfn, lapply(lapply(args, rec), untuple)))
  if (!is.null(.FUNS[[f]])) return(fcall(f, lapply(args, rec)))
  if (grepl("^[A-Za-z][A-Za-z0-9_.]*$", f)) return(fcall(f, lapply(args, rec)))
  stop("unsupported syntax '", f, "' in ", paste(deparse(x), collapse = " "))
}

.reg_parse_defaults <- function() {
  p <- function(nm, fn) register_parse_function(nm, fn)
  p("sqrt", function(x) .sqrt(x))
  p("abs", function(x) fcall("Abs", list(x)))
  p("Rational", function(a, b = S_One) .div(a, b))
  p("Integer", function(a) a)
  p("Float", function(a) Float(a))
  p("S", function(a) a)
  p("Symbol", function(a) a)
  p("ln", function(x) fcall("log", list(x)))
  p("rf", function(x, k) fcall("RisingFactorial", list(x, k)))
  p("ff", function(x, k) fcall("FallingFactorial", list(x, k)))
  p("Pow", function(a, b) .pow(a, b))
  p("Add", function(...) .add(list(...)))
  p("Mul", function(...) .mul(list(...)))
  p("Eq", function(a, b = S_Zero) Eq(a, b))
  p("Ne", function(a, b) Ne(a, b))
  p("Lt", function(a, b) Lt(a, b))
  p("Le", function(a, b) Le(a, b))
  p("Gt", function(a, b) Gt(a, b))
  p("Ge", function(a, b) Ge(a, b))
  p("exp", function(x) fcall("exp", list(x)))
  p("cbrt", function(x) .pow(x, .rat(1, 3)))
  p("root", function(x, n) .pow(x, .inv(n)))
}
.reg_parse_defaults()

# Python tuples, a parenthesized group with a top-level comma that is not a
# call, become .tuple(...) calls R can parse.
.mark_tuples <- function(txt) {
  ch <- strsplit(txt, "", fixed = TRUE)[[1L]]
  n <- length(ch)
  if (!n || !any(ch %in% c(",", "[", "("))) return(txt)
  # Python tuples "(a, b)", "(a,)", "()" and lists "[a, b]" become .tuple(...)
  open_tuple <- logical(n); drop <- logical(n); close_paren <- logical(n)
  stack <- integer(); commas <- integer()
  prev_nonspace <- function(i) { j <- i - 1L; while (j >= 1L && ch[j] == " ") j <- j - 1L; j }
  for (i in seq_len(n)) {
    if (ch[i] %in% c("(", "[")) {
      j <- prev_nonspace(i)
      iscall <- j >= 1L && grepl("[]A-Za-z0-9_.)]", ch[j])
      stack <- c(stack, if (iscall) -i else i); commas <- c(commas, 0L)
    } else if (ch[i] %in% c(")", "]") && length(stack)) {
      k <- stack[length(stack)]; top <- length(stack)
      j <- prev_nonspace(i)
      trailing <- j >= 1L && ch[j] == ","
      if (k > 0L && (ch[k] == "[" || commas[top] > 0L || j == k)) {
        open_tuple[k] <- TRUE
        if (ch[k] == "[") close_paren[i] <- TRUE
        if (trailing) drop[j] <- TRUE
      }
      stack <- stack[-top]; commas <- commas[-top]
    } else if (ch[i] == "," && length(stack)) commas[length(commas)] <- commas[length(commas)] + 1L
  }
  if (!any(open_tuple)) return(txt)
  ch[open_tuple] <- ".tuple("
  ch[close_paren] <- ")"
  ch[drop] <- ""
  paste(ch, collapse = "")
}
