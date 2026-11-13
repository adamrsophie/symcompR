# The algebra of Pauli matrices as noncommutative symbols, after
# sympy/physics/paulialgebra.py.

.pa_delta <- function(i, j) if (i == j) 1L else 0L
.pa_epsilon <- function(i, j, k) {
  t <- paste(i, j, k)
  if (t %in% c("1 2 3", "2 3 1", "3 1 2")) return(1L)
  if (t %in% c("1 3 2", "3 2 1", "2 1 3")) return(-1L)
  0L
}

.is_pauli <- function(x) is_sym(x) && identical(x$extra$sympy_class, "Pauli")

#' Pauli algebra
#'
#' `Pauli(i, label)` is the Pauli matrix `sigma_i` (`i` in 1, 2, 3) of
#' `sympy.physics.paulialgebra` as a noncommutative Hermitian symbol named
#' `label` followed by `i`. The product of two of them with the same label
#' is evaluated, and so are their positive integer powers.
#' `evaluate_pauli_product()` evaluates the products of Pauli symbols in an
#' expression.
#'
#' @param i 1, 2 or 3.
#' @param label The name of the symbol.
#' @param arg An expression.
#' @return A symbol or an expression.
#' @examples
#' s1 <- Pauli(1); s2 <- Pauli(2); s3 <- Pauli(3)
#' s1 * s2
#' evaluate_pauli_product(sym_I * s2 * s3)
#' @name pauli-algebra
NULL

#' @rdname pauli-algebra
#' @export
Pauli <- function(i, label = "sigma") {
  ii <- if (is_sym(i)) .qint(i) else i
  if (!(length(ii) == 1L && ii %in% 1:3)) .qstop("IndexError", "Invalid Pauli index")
  ii <- as.integer(ii)
  given <- .normalize_facts(list(commutative = FALSE, hermitian = TRUE))
  .nc_state$any <- TRUE
  new_sym("Symbol", name = sprintf("%s%d", label, ii),
          extra = list(facts = .symbol_facts(given), given = given, sympy_class = "Pauli",
                       pauli_i = ii, pauli_label = label),
          key = paste0("Pauli:", ii, ":", nchar(label), ":", label))
}

# Pauli.__mul__: the product of two Pauli symbols with the same label
.pauli_mul <- function(a, b) {
  if (!(.is_pauli(a) && .is_pauli(b) && identical(a$extra$pauli_label, b$extra$pauli_label))) return(NULL)
  j <- a$extra$pauli_i; k <- b$extra$pauli_i; lab <- a$extra$pauli_label
  .add(list(as_sym(.pa_delta(j, k)),
            .mul(list(sym_I, as_sym(.pa_epsilon(j, k, 1L)), Pauli(1L, lab))),
            .mul(list(sym_I, as_sym(.pa_epsilon(j, k, 2L)), Pauli(2L, lab))),
            .mul(list(sym_I, as_sym(.pa_epsilon(j, k, 3L)), Pauli(3L, lab)))))
}

# the binary product a*b with Python's operator dispatch
.pauli_times <- function(a, b) {
  r <- .pauli_mul(a, b)
  if (!is.null(r)) return(r)
  .mul(list(a, b))
}

local({
  old <- .ops_hooks$dispatch
  .ops_hooks$dispatch <- function(gen, e1, e2) {
    if (gen == "*" && .is_pauli(e1) && .is_pauli(e2)) {
      r <- .pauli_mul(e1, e2)
      if (!is.null(r)) return(r)
    }
    if (is.null(old)) NULL else old(gen, e1, e2)
  }
  oldpow <- .arith_hooks$pow_pre
  .arith_hooks$pow_pre <- function(b, e) {
    if (.is_pauli(b) && is_sym(e) && e$type == "Integer" && .qint(e) > 0)
      return(if (.qint(e) %% 2 == 0) S_One else b)
    if (is.null(oldpow)) NULL else oldpow(b, e)
  }
})

.pa_as_coeff_mul <- function(e) {
  if (e$type == "Mul") {
    a <- .sympy_args(e)
    if (is_number_type(a[[1L]])) return(list(a[[1L]], a[-1L]))
    return(list(S_One, a))
  }
  list(S_One, list(e))
}

#' @rdname pauli-algebra
#' @export
evaluate_pauli_product <- function(arg) {
  arg <- as_sym(arg)
  start <- arg; end <- arg
  if (arg$type == "Pow" && .is_pauli(arg$args[[1L]])) {
    if (isTRUE(ask_fact(arg$args[[2L]], "odd"))) return(arg$args[[1L]])
    return(S_One)
  }
  if (arg$type == "Add") return(.add(lapply(.sympy_args(arg), evaluate_pauli_product)))
  if (arg$type == "QTensorProduct") return(.qtp(lapply(arg$args, evaluate_pauli_product)))
  if (arg$type != "Mul") return(arg)
  first <- TRUE
  while (!.qeq(start, end) || (first && .qeq(start, arg) && .qeq(end, arg))) {
    first <- FALSE
    start <- end
    tmp <- .pa_as_coeff_mul(start)
    sigma_product <- S_One; com_product <- S_One; keeper <- S_One
    for (el in tmp[[2L]]) {
      if (.is_pauli(el)) {
        sigma_product <- .pauli_times(sigma_product, el)
      } else if (!isTRUE(ask_fact(el, "commutative"))) {
        if (el$type == "Pow" && .is_pauli(el$args[[1L]])) {
          if (isTRUE(ask_fact(el$args[[2L]], "odd"))) sigma_product <- .pauli_times(sigma_product, el$args[[1L]])
        } else if (el$type == "QTensorProduct") {
          keeper <- .pauli_times(.pauli_times(keeper, sigma_product), .qtp(lapply(el$args, evaluate_pauli_product)))
          sigma_product <- S_One
        } else {
          keeper <- .pauli_times(.pauli_times(keeper, sigma_product), el)
          sigma_product <- S_One
        }
      } else {
        com_product <- .pauli_times(com_product, el)
      }
    }
    end <- .pauli_times(.pauli_times(.pauli_times(tmp[[1L]], keeper), sigma_product), com_product)
    if (.qeq(end, arg)) break
  }
  end
}
