# Rewrite rules of the hypergeometric functions (R/functions-hyper.R); this
# file loads after .rewrite_rules is created in R/piecewise-rewrite.R.

.rewrite_rules$Sum$hyper <- .hyper_as_sum
for (.nm in names(.HYPERREP)) local({
  name <- .nm
  .rewrite_rules$nonrep[[name]] <<- function(...) .hyperrep_nonrep(name, list(...))
  .rewrite_rules$nonrepsmall[[name]] <<- function(...) .hyperrep_nonrepsmall(name, list(...))
})
rm(.nm)

# MinMaxBase._eval_rewrite_as_Abs: (a + M)/2 +- |a - M|/2 with M the Max or
# Min of the remaining arguments
.minmax_as_abs <- function(name) function(...) {
  a <- list(...)
  if (length(a) == 1L) return(a[[1L]])
  # LatticeOp.args: ordered by default_sort_key
  a <- .key_sort(lapply(a, .sort_key), a)
  rest <- if (length(a) == 2L) a[[2L]] else fcall(name, a[-1L])
  s <- .div(.add(list(a[[1L]], rest)), S_Two)
  d <- .div(fcall("Abs", list(.sub(a[[1L]], rest))), S_Two)
  rewrite(if (name == "Max") .add(list(s, d)) else .sub(s, d), "Abs")
}
.rewrite_rules$Abs$Max <- .minmax_as_abs("Max")
.rewrite_rules$Abs$Min <- .minmax_as_abs("Min")

# _eval_rewrite_as_tanh of the hyperbolic functions (tanh of half the
# argument) and sin._eval_rewrite_as_sinc
.th_half <- function(a) fcall("tanh", list(.mul(list(S_Half, a))))
.rewrite_rules$tanh <- list(
  sinh = function(a) { th <- .th_half(a); .div(.mul(list(S_Two, th)), .sub(S_One, .sq(th))) },
  cosh = function(a) { th <- .th_half(a); .div(.add(list(S_One, .sq(th))), .sub(S_One, .sq(th))) },
  coth = function(a) .inv(fcall("tanh", list(a))),
  csch = function(a) { th <- .th_half(a); .div(.sub(S_One, .sq(th)), .mul(list(S_Two, th))) },
  sech = function(a) { th <- .th_half(a); .div(.sub(S_One, .sq(th)), .add(list(S_One, .sq(th)))) })
.rewrite_rules$sinc <- list(sin = function(a) .mul(list(a, fcall("sinc", list(a)))))

# SingularityFunction._eval_rewrite_as_Heaviside
.rewrite_rules$Heaviside$SingularityFunction <- function(x, a, n) {
  H <- fcall("Heaviside", list(.sub(x, a)))
  if (identical(n$key, .int(-2L)$key)) return(diff(H, free_symbols(x)[[1L]], 2))
  if (identical(n$key, S_NegativeOne$key)) return(diff(H, free_symbols(x)[[1L]], 1))
  if (isTRUE(ask_fact(n, "nonnegative"))) return(.mul(list(.pow(.sub(x, a), n), fcall("Heaviside", list(.sub(x, a), S_One)))))
  NULL
}
