# The Factorization handlers of sympy/assumptions/handlers/matrices.py and
# the predicates of sympy/matrices/expressions/factorizations.py.

# Factorization.predicates of each factor class
.docsmm_factor_predicates <- list(
  LofLU = "lower_triangular", UofLU = "upper_triangular",
  LofCholesky = "lower_triangular", UofCholesky = "upper_triangular",
  QofQR = "orthogonal", RofQR = "upper_triangular",
  EigenVectors = "orthogonal", EigenValues = "diagonal",
  UofSVD = "orthogonal", SofSVD = "diagonal", VofSVD = "orthogonal")

# _Factorization: True when the predicate is one of expr.predicates
.docsmm_factorization <- function(key) {
  force(key)
  function(expr, assumptions) if (key %in% .docsmm_factor_predicates[[expr$type]]) TRUE else NA
}

invisible(lapply(c("orthogonal", "unitary", "upper_triangular", "lower_triangular", "diagonal"),
                 function(k) .ask_register(k, "Factorization", .docsmm_factorization(k))))
.ask_register("real_elements", "Factorization", .closed_elements("real_elements"))
.ask_register("complex_elements", "Factorization", .closed_elements("complex_elements"))
