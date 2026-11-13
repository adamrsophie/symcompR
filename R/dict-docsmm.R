# Item access on a Python dictionary: d$name is d["name"] for a string key,
# as the dictionaries of classify_ode() and classify_sysode() are read in
# SymPy. The fields keys and values hold the dictionary itself.

#' @export
`$.pydict` <- function(x, name) {
  if (name %in% c("keys", "values")) return(.subset2(x, name))
  ks <- .subset2(x, "keys")
  for (i in seq_along(ks))
    if (is.character(ks[[i]]) && identical(ks[[i]], name)) return(.subset2(x, "values")[[i]])
  NULL
}
