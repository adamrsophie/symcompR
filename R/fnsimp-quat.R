# Quaternion.add, _generic_mul with a scalar and Quaternion.pow of
# sympy/algebras/quaternion.py.

# A scalar z with a quaternion q: over the real field a complex z is the
# quaternion re(z) + im(z)*i; any other commutative z scales or shifts q.
.fs_quat_scalar <- function(op, e1, e2) {
  qleft <- inherits(e1, "Quaternion")
  q <- if (qleft) e1 else e2
  z <- as_sym(if (qleft) e2 else e1)
  rf <- !identical(q$real_field, FALSE)
  if (rf && isTRUE(ask_fact(z, "complex"))) {
    zq <- Quaternion(fcall("re", list(z)), fcall("im", list(z)), 0, 0)
    if (op == "+") return(q + zq)
    return(if (qleft) q * zq else zq * q)
  }
  if (op == "+") {
    cs <- q$c
    cs[[1L]] <- .add(list(cs[[1L]], z))
    return(Quaternion(cs[[1L]], cs[[2L]], cs[[3L]], cs[[4L]], real_field = q$real_field))
  }
  cs <- lapply(q$c, function(u) .mul(list(z, u)))
  Quaternion(cs[[1L]], cs[[2L]], cs[[3L]], cs[[4L]], real_field = q$real_field)
}

# Quaternion.pow: binary powering of an integer power, the inverse first
# for a negative one
.fs_quat_pow <- function(q, p) {
  pz <- if (is_sym(p)) p else tryCatch(as_sym(p), error = function(e) NULL)
  if (!inherits(q, "Quaternion") || is.null(pz) || pz$type != "Integer") {
    cls <- if (is.null(pz)) class(p)[1L] else .fs_py_class(pz)
    stop(.sympy_condition("TypeError", sprintf("unsupported operand type(s) for ** or pow(): 'Quaternion' and '%s'", cls)))
  }
  k <- .int_val(pz)
  if (k < 0) {
    q <- quat_inverse(q)
    k <- -k
  }
  if (k == 1) return(q)
  res <- Quaternion(1, 0, 0, 0)
  while (k > 0) {
    if (k %% 2 == 1) res <- res * q
    q <- q * q
    k <- k %/% 2
  }
  res
}

# The Python class name SymPy shows for a non integer exponent
.fs_py_class <- function(e) {
  if (e$type == "Rational" && identical(e$key, S_Half$key)) return("Half")
  if (e$type == "Constant") return(switch(EXPR = e$name, pi = "Pi", E = "Exp1", I = "ImaginaryUnit", e$name))
  e$type
}
