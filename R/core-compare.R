# Basic.compare: SymPy's internal order of the arguments of Add and Mul
# (`_addsort`, `_mulsort`). Algorithms that walk `expr.args` in SymPy see the
# arguments in this order, so ports of such algorithms sort with it.

.ORDERING_OF_CLASSES <- c("Zero", "One", "Half", "Infinity", "NaN", "NegativeOne", "NegativeInfinity",
  "Integer", "Rational", "Float", "Exp1", "Pi", "ImaginaryUnit", "Symbol", "Wild", "Pow", "Mul", "Add",
  "Derivative", "Integral", "Abs", "Sign", "Sqrt", "Floor", "Ceiling", "Re", "Im", "Arg", "Conjugate",
  "Exp", "Log", "Sin", "Cos", "Tan", "Cot", "ASin", "ACos", "ATan", "ACot", "Sinh", "Cosh", "Tanh",
  "Coth", "ASinh", "ACosh", "ATanh", "ACoth", "RisingFactorial", "FallingFactorial", "factorial",
  "binomial", "Gamma", "LowerGamma", "UpperGamma", "PolyGamma", "Erf", "Chebyshev", "Chebyshev2",
  "Function", "WildFunction", "Lambda", "Order", "Equality", "Unequality", "StrictGreaterThan",
  "StrictLessThan", "GreaterThan", "LessThan")

# the SymPy class name of a node
.sympy_class <- function(e) {
  t <- e$type
  if (t == "Symbol" && !is.null(e$extra$sympy_class)) return(e$extra$sympy_class)
  if (t == "Integer") {
    if (e$value == 0) return("Zero")
    if (e$value == 1) return("One")
    if (e$value == -1) return("NegativeOne")
    return("Integer")
  }
  if (t == "Rational") return(if (e$value == gmp::as.bigq(1, 2)) "Half" else "Rational")
  if (t == "Constant") return(switch(EXPR = e$name, oo = "Infinity", "-oo" = "NegativeInfinity", zoo = "ComplexInfinity",
                                     nan = "NaN", E = "Exp1", pi = "Pi", I = "ImaginaryUnit", e$name))
  if (t == "Function") return(e$name)
  if (t == "BooleanAtom") return(if (e$name == "true") "BooleanTrue" else "BooleanFalse")
  if (t == "Quantity") return(if (isTRUE(e$extra$pc)) "PhysicalConstant" else "Quantity")
  t
}

# Python compares strings by code point; R's > follows the locale collation
.cmp_codepoints <- function(a, b) {
  if (identical(a, b)) return(0L)
  if (identical(sort(c(a, b), method = "radix")[1L], a)) -1L else 1L
}

.cmp_name <- function(n1, n2) {
  if (n1 == n2) return(0L)
  i1 <- match(n1, .ORDERING_OF_CLASSES); i2 <- match(n2, .ORDERING_OF_CLASSES)
  if (is.na(i1) && is.na(i2)) return(.cmp_codepoints(n1, n2))
  i1 <- if (is.na(i1)) length(.ORDERING_OF_CLASSES) + 1L else i1
  i2 <- if (is.na(i2)) length(.ORDERING_OF_CLASSES) + 1L else i2
  sign(i1 - i2)
}

# _hashable_content as a list of sym nodes and plain values
.hashable <- function(e) {
  t <- e$type
  if (t == "Integer") return(list(as.double(e$value)))
  if (t == "Rational") return(list(as.double(gmp::numerator(e$value)), as.double(gmp::denominator(e$value))))
  if (t == "Float") return(list(as.double(e$value)))
  # a BaseScalar hashes its index and system; a CoordSys3D its name first
  if (t %in% c("Symbol", "Wild")) return(if (isTRUE(e$extra$base_scalar)) list(.int(e$extra$index), e$extra$system) else list(e$name))
  if (t == "CoordSys3D") return(list(e$name))
  if (t == "Dummy") return(list(e$name, e$key))
  # Quantity args are its name and abbreviation
  if (t == "Quantity") return(list(Symbol(e$name), e$extra$abbrev))
  if (t %in% c("Add", "Mul")) return(.sympy_args(e))
  if (!length(e$args)) return(list())
  e$args
}

.basic_compare <- function(a, b) {
  if (identical(a$key, b$key)) return(0L)
  c <- .cmp_name(.sympy_class(a), .sympy_class(b))
  if (c != 0L) return(c)
  st <- .hashable(a); ot <- .hashable(b)
  if (length(st) != length(ot)) return(if (length(st) > length(ot)) 1L else -1L)
  for (i in seq_along(st)) {
    l <- st[[i]]; r <- ot[[i]]
    c <- if (is_sym(l) && is_sym(r)) .basic_compare(l, r) else if (identical(l, r)) 0L
         else if (is.character(l) && is.character(r)) .cmp_codepoints(l, r) else if (l > r) 1L else -1L
    if (c != 0L) return(c)
  }
  0L
}

# the args of an Add or Mul in SymPy's internal order (insertion sort, stable)
.sympy_args <- function(e) {
  args <- e$args
  if (length(args) < 2L) return(args)
  # Mul: the commutative part sorted, the noncommutative part in order
  if (e$type == "Mul") {
    comm <- vapply(args, .is_comm, TRUE)
    if (!all(comm)) {
      cp <- args[comm]
      if (length(cp) > 1L) cp <- .sympy_args(new_sym("Mul", cp))
      return(c(cp, args[!comm]))
    }
  }
  o <- integer()
  for (i in seq_along(args)) {
    j <- length(o)
    while (j > 0L && .basic_compare(args[[o[j]]], args[[i]]) > 0L) j <- j - 1L
    o <- append(o, i, after = j)
  }
  args[o]
}
