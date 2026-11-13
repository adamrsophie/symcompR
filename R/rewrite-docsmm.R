# Rewrite rules of subfactorial, after
# sympy/functions/combinatorial/factorials.py. This file loads after
# .rewrite_rules is created in R/piecewise-rewrite.R.

# subfactorial._eval_rewrite_as_uppergamma
.rewrite_rules$uppergamma$subfactorial <- function(a)
  .div(fcall("uppergamma", list(.add(list(a, S_One)), S_NegativeOne)), S_E)

# subfactorial._eval_rewrite_as_gamma
.rewrite_rules$gamma$subfactorial <- function(a) {
  a1 <- .add(list(a, S_One))
  .mul(list(.add(list(.mul(list(.pow(S_NegativeOne, a1), fcall("exp", list(.neg(.mul(list(S_I, S_Pi, a))))),
                                fcall("lowergamma", list(a1, S_NegativeOne)))),
                      fcall("gamma", list(a1)))),
            fcall("exp", list(S_NegativeOne))))
}
