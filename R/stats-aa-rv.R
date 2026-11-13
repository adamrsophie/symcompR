# Random variables, after sympy/stats/rv.py: random symbols, probability
# spaces and domains, and the user functions pspace(), random_symbols(),
# given(), E(), P(), rv_density(), rv_cdf(), characteristic_function(),
# moment_generating_function(), rv_quantile(), where(), rv_sample(),
# independent() and dependent().
#
# Probability spaces, domains and distributions are sym nodes whose type is
# the SymPy class name, so they print and compare as SymPy's Basic objects.
# Their methods live in a small class table (.defm) searched along the
# Python method resolution order (.MRO).

# ---------------------------------------------------------------------------
# Class table

.STM <- new.env(hash = TRUE)
.MRO <- new.env(hash = TRUE)
.PARENTS <- new.env(hash = TRUE)

.defclass <- function(cls, parents = character()) {
  assign(cls, parents, envir = .PARENTS)
  invisible(cls)
}

# C3 linearization, as Python computes a class's method resolution order
.c3 <- function(cls) {
  hit <- .MRO[[cls]]
  if (!is.null(hit)) return(hit)
  parents <- .PARENTS[[cls]]
  seqs <- c(lapply(parents, .c3), list(parents))
  seqs <- Filter(length, seqs)
  out <- cls
  while (length(seqs)) {
    cand <- NULL
    for (s in seqs) {
      h <- s[1L]
      if (!any(vapply(seqs, function(t) h %in% t[-1L], NA))) { cand <- h; break }
    }
    if (is.null(cand)) stop("inconsistent class hierarchy for ", cls)
    out <- c(out, cand)
    seqs <- Filter(length, lapply(seqs, function(t) t[t != cand]))
  }
  assign(cls, out, envir = .MRO)
  out
}

.defm <- function(cls, name, fn) {
  tab <- .STM[[cls]]
  if (is.null(tab)) tab <- list()
  tab[[name]] <- fn
  assign(cls, tab, envir = .STM)
  invisible(NULL)
}

.mro <- function(obj) {
  if (is.null(.PARENTS[[obj$type]])) return(obj$type)
  .c3(obj$type)
}

.findm <- function(obj, name) {
  for (c in .mro(obj)) {
    f <- .STM[[c]][[name]]
    if (!is.null(f)) return(f)
  }
  NULL
}

.callm <- function(obj, name, ...) {
  f <- .findm(obj, name)
  if (is.null(f)) .serr("AttributeError", paste0("'", obj$type, "' object has no attribute '", name, "'"))
  f(obj, ...)
}

.isa <- function(obj, cls) is_sym(obj) && cls %in% .mro(obj)

.serr <- function(cls, msg) stop(structure(class = c(cls, "error", "condition"), list(message = msg, call = NULL)))

# ---------------------------------------------------------------------------
# Helpers mirroring Python idioms

# _value_check: raise ValueError when a condition is False; TRUE when all
# conditions are True
.fuzzy_bool <- function(c) {
  if (is.null(c)) return(NA)
  if (is.logical(c)) return(if (length(c) == 1L) c else NA)
  if (is_sym(c) && c$type == "BooleanAtom") return(c$name == "true")
  NA
}

.value_check <- function(condition, message) {
  if (!is.list(condition) || is_sym(condition)) condition <- list(condition)
  rv <- TRUE
  for (c in condition) {
    a <- .fuzzy_bool(c)
    if (identical(a, FALSE)) .serr("ValueError", message)
    if (isTRUE(rv)) rv <- a
  }
  isTRUE(rv)
}

# Python's `a and b`: a when a is falsy (False or None), else b
.pyand <- function(a, b) if (isTRUE(.fuzzy_bool(a))) b else a

.symbol_converter <- function(sym) {
  if (is.character(sym)) sym <- Symbol(sym)
  if (!(is_sym(sym) && sym$type %in% c("Symbol", "Dummy")))
    .serr("TypeError", paste0(.rv_pystr(sym), " is neither a Symbol nor a string"))
  sym
}

.rv_pystr <- function(x) {
  if (is.character(x)) return(x)
  if (is.logical(x)) return(if (isTRUE(x)) "True" else "False")
  if (is.numeric(x)) return(format(x))
  if (is.list(x) && !is_sym(x)) return(paste0("[", paste(vapply(x, .rv_pystr, ""), collapse = ", "), "]"))
  sstr(x)
}

.st_fact <- function(e, f) {
  e <- as_sym(e)
  # SymPy's Float.is_integer: True for zero, False with a fractional part,
  # None otherwise
  if (f == "integer" && e$type == "Float") {
    v <- as.numeric(e$value)
    return(if (v == 0) TRUE else if (v != round(v)) FALSE else NA)
  }
  ask_fact(e, f)
}

# Python truth of a sympy comparison (Relational.__bool__ raises)
.truth <- function(x) {
  if (is.logical(x)) return(isTRUE(x))
  if (is_sym(x) && x$type == "BooleanAtom") return(x$name == "true")
  .serr("TypeError", "cannot determine truth value of Relational")
}

.as_bool_sym <- function(x) if (is.logical(x)) as_sym(x) else x

.eq_struct <- function(a, b) is_sym(a) && is_sym(b) && identical(a$key, b$key)

.uniq <- function(l) l[!duplicated(vapply(l, .k, ""))]

.in_list <- function(x, l) length(l) > 0L && x$key %in% vapply(l, .k, "")

.sorted_basic <- function(l) if (length(l) < 2L) l else .key_sort(lapply(l, .sort_key), l)

# a FiniteSet of Basic objects keeps SymPy's order: infimum first, then the
# default sort key; objects without an infimum sort by the key alone
.fset_basic <- function(l) .sorted_basic(.uniq(l))

# ---------------------------------------------------------------------------
# Random symbols

.PSpace0 <- function() new_sym("PSpace", list())

.RandomSymbol <- function(symbol, pspace = NULL, cls = "RandomSymbol") {
  if (is.null(pspace)) pspace <- .PSpace0()
  if (cls != "RandomIndexedSymbol") symbol <- .symbol_converter(symbol)
  if (!.isa(pspace, "PSpace")) .serr("TypeError", "pspace variable should be of type PSpace")
  if (cls == "JointRandomSymbol" && .isa(pspace, "SinglePSpace")) cls <- "RandomSymbol"
  sf <- if (symbol$type %in% c("Symbol", "Dummy")) symbol$extra$facts else .facts(symbol)
  ps_real <- isTRUE(.ps_is_real(pspace))
  given <- c(finite = TRUE)
  if (isTRUE(sf[["real"]]) || ps_real) given["real"] <- TRUE
  if (!is.na(sf[["positive"]])) given["positive"] <- sf[["positive"]]
  if (!is.na(sf[["integer"]])) given["integer"] <- sf[["integer"]]
  if (identical(sf[["commutative"]], FALSE)) given["commutative"] <- FALSE
  facts <- .symbol_facts(given)
  name <- if (cls == "RandomIndexedSymbol") .str(symbol) else symbol$name
  new_sym("Symbol", name = name,
          extra = list(facts = facts, sympy_class = cls, rv_class = cls, symbol = symbol, pspace = pspace),
          key = paste0("RS[", cls, "](", symbol$key, ";", pspace$key, ")"))
}

.is_rs <- function(e) is_sym(e) && e$type %in% c("Symbol", "RandomIndexedSymbol") && !is.null(e$extra$rv_class)
.rs_symbol <- function(rv) if (rv$type == "RandomIndexedSymbol") rv$args[[1L]] else rv$extra$symbol
.rs_pspace <- function(rv) if (rv$type == "RandomIndexedSymbol") rv$args[[2L]] else rv$extra$pspace

register_str_printer("RandomSymbolNode", function(e) e$name)

#' Random symbols in an expression
#'
#' `random_symbols(expr)` lists the random symbols of `expr`, sorted by
#' name; `is_random(expr)` tells whether there are any; `pspace(expr)` is
#' the probability space they live on. These follow `sympy.stats`.
#'
#' @param expr An expression.
#' @return A list of random symbols, a logical, or a probability space.
#' @examples
#' X <- Normal("X", 0, 1); Y <- Normal("Y", 0, 1)
#' random_symbols(2 * X + Y)
#' is_random(X^2)
#' pspace(X)
#' @export
random_symbols <- function(expr) {
  if (is.null(expr) || is.logical(expr) || is.numeric(expr) || is.character(expr)) return(list())
  if (inherits(expr, "symmatrix")) return(.rv_sort(.uniq(unlist(lapply(expr$data, .rv_atoms), recursive = FALSE))))
  if (!is_sym(expr)) return(list())
  .rv_sort(.rv_atoms(expr))
}

.rv_atoms <- function(e) {
  out <- list()
  walk <- function(x) {
    if (.is_rs(x)) { out[[x$key]] <<- x; return(invisible()) }
    if (!is.null(x$extra$rv_matrix)) { out[[x$key]] <<- x; return(invisible()) }
    for (a in x$args) if (is_sym(a)) walk(a)
  }
  walk(e)
  unname(out)
}

.rv_sort <- function(l) {
  if (length(l) < 2L) return(l)
  l[order(vapply(l, function(r) .rv_sym_name(r), ""), method = "radix")]
}

.rv_sym_name <- function(r) {
  s <- .rs_symbol(r)
  if (is.null(s)) return(r$name)
  if (s$type %in% c("Symbol", "Dummy")) s$name else .str(s)
}

#' @rdname random_symbols
#' @export
is_random <- function(expr) {
  if (is.null(expr) || is.logical(expr) || is.numeric(expr) || is.character(expr)) return(FALSE)
  if (inherits(expr, "symmatrix")) return(any(vapply(expr$data, is_random, NA)))
  if (inherits(expr, "rv_process")) return(TRUE)
  if (!is_sym(expr)) return(FALSE)
  if (.is_rs(expr)) return(TRUE)
  if (!is.null(expr$extra$rv_matrix)) return(TRUE)
  length(.rv_atoms(expr)) > 0L || any(vapply(free_symbols(expr), .is_rs, NA))
}

#' @rdname random_symbols
#' @export
pspace <- function(expr) {
  if (inherits(expr, "rv_process")) return(.process_pspace(expr))
  expr <- if (inherits(expr, "symmatrix")) expr else as_sym(expr)
  if (.is_rs(expr) && !is.null(.rs_pspace(expr))) return(.rs_pspace(expr))
  if (is_sym(expr)) {
    rm <- Filter(function(z) !is.null(z$extra$rv_matrix), .all_subexprs(expr))
    if (length(rm)) return(rm[[1L]]$extra$pspace)
  }
  rvs <- random_symbols(expr)
  if (!length(rvs)) .serr("ValueError", paste0("Expression containing Random Variable expected, not ", .rv_pystr(expr)))
  ps <- lapply(rvs, .rs_pspace)
  if (all(vapply(ps, function(p) identical(p$key, ps[[1L]]$key), NA))) return(ps[[1L]])
  for (p in ps) if (.isa(p, "CompoundPSpace") || .isa(p, "StochasticPSpace")) return(p)
  do.call(.IndependentProductPSpace, ps)
}

# ---------------------------------------------------------------------------
# Probability spaces

.defclass("PSpace")
.defclass("SinglePSpace", "PSpace")
.defclass("ProductPSpace", "PSpace")
.defclass("IndependentProductPSpace", "ProductPSpace")

.ps_is_real <- function(ps) {
  f <- .findm(ps, "is_real")
  if (is.null(f)) NA else f(ps)
}
.ps_flag <- function(ps, flag) {
  f <- .findm(ps, flag)
  if (is.null(f)) NA else f(ps)
}
.ps_cont <- function(ps) isTRUE(.ps_flag(ps, "is_Continuous"))
.ps_disc <- function(ps) isTRUE(.ps_flag(ps, "is_Discrete"))
.ps_fin <- function(ps) isTRUE(.ps_flag(ps, "is_Finite"))

.ps_domain <- function(ps) .callm(ps, "domain")
.ps_values <- function(ps) .callm(ps, "values")
.ps_symbols <- function(ps) .callm(ps, "symbols")
.ps_pdf <- function(ps) .callm(ps, "pdf")

.defm("PSpace", "domain", function(ps) ps$args[[1L]])
.defm("PSpace", "density", function(ps) ps$args[[2L]])
.defm("PSpace", "symbols", function(ps) .dom_symbols(.ps_domain(ps)))
.defm("PSpace", "values", function(ps) .fset_basic(lapply(.ps_symbols(ps), function(s) .RandomSymbol(s, ps))))
.defm("PSpace", "where", function(ps, condition) .serr("NotImplementedError", ""))
.defm("PSpace", "compute_density", function(ps, expr, ...) .serr("NotImplementedError", ""))
.defm("PSpace", "probability", function(ps, condition, ...) .serr("NotImplementedError", ""))
.defm("PSpace", "compute_expectation", function(ps, expr, ...) .serr("NotImplementedError", ""))

.defm("SinglePSpace", "symbol", function(ps) ps$args[[1L]])
.defm("SinglePSpace", "distribution", function(ps) ps$args[[2L]])
.defm("SinglePSpace", "value", function(ps) .RandomSymbol(ps$args[[1L]], ps))
.defm("SinglePSpace", "pdf", function(ps) .dist_pdf(ps$args[[2L]], ps$args[[1L]]))

.ps_value <- function(ps) .callm(ps, "value")

# IndependentProductPSpace(*spaces)
.IndependentProductPSpace <- function(...) {
  spaces <- list(...)
  rs <- list()
  for (sp in spaces) for (v in .ps_values(sp)) rs[[v$key]] <- v
  syms <- .uniq(lapply(unname(rs), .rs_symbol))
  nsym <- sum(vapply(spaces, function(sp) {
    d <- .findm(sp, "distribution")
    dist <- if (is.null(d)) NULL else d(sp)
    if (!is.null(dist) && (dist$type %in% c("CompoundDistribution", "MarginalDistribution"))) 0L else length(.ps_symbols(sp))
  }, 0L))
  if (length(syms) < nsym) .serr("ValueError", "Overlapping Random Variables")
  cls <- if (all(vapply(spaces, .ps_fin, NA))) "ProductFinitePSpace" else "IndependentProductPSpace"
  new_sym(cls, .fset_basic(spaces))
}

.defm("IndependentProductPSpace", "spaces", function(ps) ps$args)
.defm("IndependentProductPSpace", "rs_space_dict", function(ps) {
  d <- list()
  for (sp in ps$args) for (v in .ps_values(sp)) d[[v$key]] <- list(v, sp)
  d
})
.defm("IndependentProductPSpace", "symbols", function(ps)
  .fset_basic(lapply(.callm(ps, "rs_space_dict"), function(p) .rs_symbol(p[[1L]]))))
.defm("IndependentProductPSpace", "values", function(ps) .uniq(unlist(lapply(ps$args, .ps_values), recursive = FALSE)))
.defm("IndependentProductPSpace", "pdf", function(ps) {
  p <- .mul(lapply(ps$args, .ps_pdf))
  vals <- .ps_values(ps)
  subs(p, lapply(vals, function(rv) list(rv, .rs_symbol(rv))))
})
.defm("IndependentProductPSpace", "domain", function(ps) do.call(.ProductDomain, lapply(ps$args, .ps_domain)))
.defm("IndependentProductPSpace", "density", function(ps) .serr("NotImplementedError", "Density not available for ProductSpaces"))
.defm("IndependentProductPSpace", "compute_expectation", function(ps, expr, rvs = NULL, evaluate = FALSE, ...) {
  if (is.null(rvs) || !length(rvs)) rvs <- .ps_values(ps)
  for (space in ps$args) {
    sv <- .ps_values(space)
    sub <- Filter(function(r) .in_list(r, sv), rvs)
    expr <- .callm(space, "compute_expectation", expr, sub, evaluate = FALSE, ...)
  }
  if (evaluate && is_sym(expr)) return(doit(expr))
  expr
})
.defm("IndependentProductPSpace", "probability", function(ps, condition, ...) {
  cond_inv <- FALSE
  if (condition$type == "Unequality") {
    condition <- Eq(condition$args[[1L]], condition$args[[2L]])
    cond_inv <- TRUE
  } else if (condition$type == "And") {
    return(.mul(lapply(.sympy_args(condition), function(a) .callm(ps, "probability", a))))
  } else if (condition$type == "Or") {
    return(.add(lapply(.sympy_args(condition), function(a) .callm(ps, "probability", a))))
  }
  expr <- .sub(condition$args[[1L]], condition$args[[2L]])
  rvs <- random_symbols(expr)
  dens <- .callm(ps, "compute_density", expr)
  if (any(vapply(rvs, function(rv) .ps_cont(pspace(rv)), NA))) {
    if (.in_list(expr, .ps_values(ps))) {
      others <- Filter(function(r) !identical(r$key, expr$key), .ps_values(ps))
      syms <- lapply(others, .rs_symbol)
      pdf <- .callm(.ps_domain(ps), "integrate", .ps_pdf(ps), syms)
      return(Lambda(.rs_symbol(expr), pdf))
    }
    dens <- .ContinuousDistributionHandmade(dens)
    z <- Dummy("z", real = TRUE)
    space <- .SingleContinuousPSpace(z, dens)
  } else {
    dens <- .DiscreteDistributionHandmade(dens)
    z <- Dummy("z", integer = TRUE)
    space <- .SingleDiscretePSpace(z, dens)
  }
  result <- .callm(space, "probability", .rel_same(condition, .ps_value(space), S_Zero))
  if (cond_inv) .sub(S_One, result) else result
})
.defm("IndependentProductPSpace", "compute_density", function(ps, expr, ...) {
  rvs <- random_symbols(expr)
  if (any(vapply(rvs, function(rv) .ps_cont(pspace(rv)), NA))) {
    z <- Dummy("z", real = TRUE)
    e <- .callm(ps, "compute_expectation", DiracDelta(.sub(expr, z)), ...)
  } else {
    z <- Dummy("z", integer = TRUE)
    e <- .callm(ps, "compute_expectation", KroneckerDelta(expr, z), ...)
  }
  Lambda(z, e)
})
.defm("IndependentProductPSpace", "compute_cdf", function(ps, expr, ...)
  .serr("ValueError", "CDF not well defined on multivariate expressions"))
.defm("IndependentProductPSpace", "conditional_space", function(ps, condition, normalize = TRUE, ...) {
  rvs <- random_symbols(condition)
  condition <- xreplace(condition, unname(lapply(.ps_values(ps), function(rv) list(rv, .rs_symbol(rv)))))
  pss <- lapply(rvs, pspace)
  if (any(vapply(pss, .ps_cont, NA))) {
    domain <- .ConditionalContinuousDomain(.ps_domain(ps), condition)
    cls <- "ContinuousPSpace"
  } else if (any(vapply(pss, .ps_disc, NA))) {
    domain <- .ConditionalDiscreteDomain(.ps_domain(ps), condition)
    cls <- "DiscretePSpace"
  } else if (all(vapply(pss, .ps_fin, NA))) {
    return(.STM[["FinitePSpace"]]$conditional_space(ps, condition))
  }
  density <- NULL
  if (normalize) {
    syms <- .dom_symbols(domain)
    replacement <- unname(lapply(.ps_symbols(ps), function(s) list(s, Dummy(.str(s)))))
    norm <- .callm(domain, "compute_expectation", .ps_pdf(ps), ...)
    pdf <- .div(.ps_pdf(ps), xreplace(norm, replacement))
    density <- Lambda(syms, pdf)
  }
  new_sym(cls, list(domain, density))
})

.rel_same <- function(rel, a, b) {
  switch(rel$type, Equality = Eq(a, b), Unequality = Ne(a, b), StrictLessThan = Lt(a, b),
         LessThan = Le(a, b), StrictGreaterThan = Gt(a, b), GreaterThan = Ge(a, b))
}

# ---------------------------------------------------------------------------
# Random domains

.defclass("RandomDomain")
.defclass("SingleDomain", "RandomDomain")
.defclass("ConditionalDomain", "RandomDomain")
.defclass("ProductDomain", "RandomDomain")

.dom_symbols <- function(d) .callm(d, "symbols")
.dom_set <- function(d) .callm(d, "set")
.dom_as_boolean <- function(d) .callm(d, "as_boolean")

.defm("RandomDomain", "symbols", function(d) .finite_elems(d$args[[1L]]))
.defm("RandomDomain", "set", function(d) d$args[[2L]])
.defm("SingleDomain", "symbol", function(d) d$args[[1L]])
.defm("SingleDomain", "symbols", function(d) list(d$args[[1L]]))
.defm("SingleDomain", "set", function(d) d$args[[2L]])
.defm("ConditionalDomain", "symbols", function(d) .dom_symbols(d$args[[1L]]))
.defm("ConditionalDomain", "fulldomain", function(d) d$args[[1L]])
.defm("ConditionalDomain", "condition", function(d) d$args[[2L]])
.defm("ConditionalDomain", "set", function(d) .serr("NotImplementedError", "Set of Conditional Domain not Implemented"))
.defm("ConditionalDomain", "as_boolean", function(d) And(.dom_as_boolean(d$args[[1L]]), d$args[[2L]]))
.defm("ProductDomain", "domains", function(d) d$args)
.defm("ProductDomain", "symbols", function(d) .fset_basic(unlist(lapply(d$args, .dom_symbols), recursive = FALSE)))
.defm("ProductDomain", "set", function(d) do.call(ProductSet, lapply(d$args, .dom_set)))
.defm("ProductDomain", "as_boolean", function(d) do.call(And, lapply(d$args, .dom_as_boolean)))
.defm("ProductDomain", "sym_domain_dict", function(d) {
  out <- list()
  for (dom in d$args) for (s in .dom_symbols(dom)) out[[s$key]] <- dom
  out
})

.ConditionalDomain_new <- function(cls, fulldomain, condition) {
  condition <- xreplace(.as_bool_sym(condition), lapply(random_symbols(condition), function(rs) list(rs, .rs_symbol(rs))))
  new_sym(cls, list(fulldomain, condition))
}

.ProductDomain <- function(...) {
  doms <- list()
  for (d in list(...)) if (.isa(d, "ProductDomain")) doms <- c(doms, d$args) else doms[[length(doms) + 1L]] <- d
  doms <- .fset_basic(doms)
  cls <- "ProductDomain"
  if (all(vapply(doms, function(d) .isa(d, "FiniteDomain"), NA))) cls <- "ProductFiniteDomain"
  if (all(vapply(doms, function(d) .isa(d, "ContinuousDomain"), NA))) cls <- "ProductContinuousDomain"
  if (all(vapply(doms, function(d) .isa(d, "DiscreteDomain"), NA))) cls <- "ProductDiscreteDomain"
  new_sym(cls, doms)
}

# StrPrinter._print_RandomDomain
.print_domain <- function(d) {
  if (!is.null(.findm(d, "as_boolean"))) {
    b <- tryCatch(.dom_as_boolean(d), error = function(e) NULL)
    if (!is.null(b)) return(paste0("Domain: ", .str(.as_bool_sym(b))))
  }
  if (!is.null(.findm(d, "set")))
    return(paste0("Domain: ", .str(do.call(FiniteSet, .dom_symbols(d))), " in ", .str(.dom_set(d))))
  paste0("Domain on ", .str(do.call(FiniteSet, .dom_symbols(d))))
}

# ---------------------------------------------------------------------------
# Distributions

.defclass("Distribution")

.ARGN <- new.env(hash = TRUE)

# A distribution class: its argument names and its methods.
# the formals start with a dot so that a method named p, a or c in ...
# is not partially matched to them
.defdist <- function(.cls, .parent, .argnames = character(), ...) {
  .defclass(.cls, .parent)
  assign(.cls, .argnames, envir = .ARGN)
  m <- list(...)
  for (n in names(m)) .defm(.cls, n, m[[n]])
  invisible(.cls)
}

# named arguments of a distribution (NamedArgsMixin)
.da <- function(d, name) {
  i <- match(name, .ARGN[[d$type]])
  if (is.na(i)) .serr("AttributeError", paste0("'", d$type, "' object has no attribute '", name, "'"))
  d$args[[i]]
}
.dl <- function(d) stats::setNames(d$args[seq_along(.ARGN[[d$type]])], .ARGN[[d$type]])

.dist_pdf <- function(d, x) .callm(d, "pdf", as_sym(x))
.dist_set <- function(d) .callm(d, "set")

# ---------------------------------------------------------------------------
# rv(): build a random symbol from a distribution class

.rv_make <- function(name, cls, args, check = TRUE, single = "SingleContinuousPSpace") {
  args <- lapply(args, function(a) if (inherits(a, "symmatrix") || is.list(a) && !is_sym(a)) a else as_sym(a))
  dist <- new_sym(cls, args)
  if (check) {
    ck <- .findm(dist, "check")
    if (!is.null(ck)) do.call(ck, args)
  }
  ps <- new_sym(single, list(.symbol_converter(name), dist))
  if (any(vapply(args, is_random, NA))) ps <- .CompoundPSpace(.symbol_converter(name), .CompoundDistribution(dist))
  .ps_value(ps)
}

# ---------------------------------------------------------------------------
# given()

#' Conditional random expressions and events
#'
#' `given(expr, condition)` re-expresses `expr` on the probability space
#' restricted by `condition`. `where(condition)` is the domain on which
#' `condition` holds. `independent()` and `dependent()` test whether two
#' random expressions depend on each other. `rs_swap(a, b)` pairs random
#' symbols with the same underlying symbol.
#'
#' @param expr A random expression.
#' @param condition A relational or Boolean in random symbols, or a random
#'   symbol (treated as fixed).
#' @param given_condition An optional further condition.
#' @param a,b Random expressions or lists of random symbols.
#' @param ... Ignored.
#' @return An expression, a domain, a logical, or a list of pairs.
#' @examples
#' X <- Normal("X", 0, 1)
#' Y <- given(X, X > 0)
#' E(Y)
#' where(X^2 <= 1)
#' D <- Die("D", 6)
#' rv_density(given(D, D > 3))
#' @export
given <- function(expr, condition = NULL, ...) {
  if (!is_random(condition) || isTRUE(.pspace_independent(expr, condition))) return(expr)
  if (.is_rs(condition)) condition <- Eq(condition, .rs_symbol(condition))
  condsymbols <- random_symbols(condition)
  if (condition$type == "Equality" && length(condsymbols) == 1L &&
      !.isa(.ps_domain(pspace(expr)), "ConditionalDomain")) {
    rv <- condsymbols[[1L]]
    results <- solveset(condition, rv)
    if (results$type == "Intersection" && any(vapply(results$args, function(a) identical(a$key, Reals$key), NA)))
      results <- results$args[[2L]]
    sums <- 0
    for (res in .finite_elems(results)) {
      temp <- subs(expr, rv, res)
      if (.is_true(.as_bool_sym(temp))) return(TRUE)
      if (!.is_false(.as_bool_sym(temp))) {
        if (identical(sums, 0) && is_relational(expr)) sums <- subs(expr, rv, res)
        else sums <- .add(list(as_sym(sums), subs(expr, rv, res)))
      }
    }
    if (identical(sums, 0) || (is_sym(sums) && is_zero_num(sums))) return(FALSE)
    return(sums)
  }
  fullspace <- pspace(Tuple(expr, condition))
  space <- .callm(fullspace, "conditional_space", condition, ...)
  swap <- rs_swap(.ps_values(fullspace), .ps_values(space))
  if (inherits(expr, "symmatrix")) return(.mat_map(expr, function(z) xreplace(z, swap)))
  xreplace(expr, swap)
}

#' @rdname given
#' @export
rs_swap <- function(a, b) {
  lapply(a, function(rsa) {
    m <- Filter(function(rsb) identical(.rs_symbol(rsa)$key, .rs_symbol(rsb)$key), b)
    list(rsa, m[[1L]])
  })
}

rv_subs <- function(expr, symbols = NULL) {
  if (is.null(symbols)) symbols <- random_symbols(expr)
  if (!length(symbols)) return(expr)
  subs(expr, lapply(symbols, function(rv) list(rv, .rs_symbol(rv))))
}

.pspace_independent <- function(a, b) {
  a_symbols <- .ps_symbols(pspace(b))
  b_symbols <- .ps_symbols(pspace(a))
  ra <- random_symbols(a); rb <- random_symbols(b)
  if (any(vapply(ra, function(r) .in_list(r, rb), NA))) return(FALSE)
  if (!any(vapply(a_symbols, function(s) .in_list(s, b_symbols), NA))) return(TRUE)
  NA
}

#' @rdname given
#' @export
dependent <- function(a, b) {
  if (isTRUE(.pspace_independent(a, b))) return(FALSE)
  z <- Symbol("z", real = TRUE)
  !.same_obj(rv_density(a, Eq(b, z)), rv_density(a)) || !.same_obj(rv_density(b, Eq(a, z)), rv_density(b))
}

.same_obj <- function(a, b) {
  if (is_sym(a) && is_sym(b)) return(identical(a$key, b$key))
  identical(a, b)
}

#' @rdname given
#' @export
independent <- function(a, b) !dependent(a, b)

#' @rdname given
#' @export
where <- function(condition, given_condition = NULL, ...) {
  if (!is.null(given_condition)) return(where(given(condition, given_condition, ...), ...))
  .callm(pspace(condition), "where", .as_bool_sym(condition))
}

# ---------------------------------------------------------------------------
# E, P, density, cdf and friends

#' Expectation, probability and distribution functions
#'
#' `E(expr, condition)` (alias `expectation()`) is the expected value of a
#' random expression; `P(condition, given_condition)` (alias
#' `probability()`) the probability of an event. `rv_density()` is SymPy's
#' `density()`: a distribution, a `Lambda`, or a finite density mapping;
#' `rv_cdf()` is `cdf()`, and `characteristic_function()`,
#' `moment_generating_function()` and `rv_quantile()` (SymPy's
#' `quantile()`) return `Lambda` functions. Apply the result to a value with
#' [apply_lambda()]. For compatibility, a non-Boolean second argument to
#' `rv_density()` and `rv_cdf()` is a point at which to evaluate. `H()`
#' (alias `entropy()`) is the entropy.
#'
#' @param expr A random expression.
#' @param condition An optional condition (Relational or Boolean).
#' @param given_condition An optional condition for `P()`.
#' @param numsamples Sampling is not supported; must be `NULL`.
#' @param evaluate `FALSE` keeps integrals and sums unevaluated.
#' @param ... Ignored, except by `Density()`: further `list(state, probability)`
#'   pairs of a quantum density operator ([density-operators]).
#' @return An expression, a distribution, a `Lambda` or a mapping.
#' @examples
#' X <- Normal("X", 0, 1)
#' E(X^2)
#' P(X > 0)
#' rv_density(X)
#' apply_lambda(rv_cdf(X), 0)
#' D <- Die("D", 6)
#' E(D)
#' P(D > 4)
#' @export
E <- function(expr, condition = NULL, numsamples = NULL, evaluate = TRUE, ...) {
  if (!is_random(expr)) return(if (is.numeric(expr) || is.logical(expr)) as_sym(expr) else expr)
  if (!is.null(numsamples)) .serr("NotImplementedError", "sampling is not supported")
  e <- Expectation(expr, condition)
  if (evaluate) return(.expectation_doit(e))
  e
}

#' @rdname E
#' @export
expectation <- E

#' @rdname E
#' @export
P <- function(condition, given_condition = NULL, numsamples = NULL, evaluate = TRUE, ...) {
  if (!is.null(numsamples)) .serr("NotImplementedError", "sampling is not supported")
  p <- Probability(condition, given_condition)
  if (evaluate) return(.probability_doit(p, evaluate = evaluate))
  p
}

#' @rdname E
#' @export
probability <- P

.is_cond <- function(x) {
  if (is.logical(x)) return(TRUE)
  is_sym(x) && (is_relational(x) || x$type %in% c("And", "Or", "Not", "BooleanAtom", "Xor", "Implies", "Equivalent") || .is_rs(x))
}

#' @rdname E
#' @export
rv_density <- function(expr, condition = NULL, evaluate = TRUE, numsamples = NULL, ...) {
  if (!is.null(numsamples)) .serr("NotImplementedError", "sampling is not supported")
  if (!is.null(condition) && !.is_cond(condition)) return(apply_lambda(rv_density(expr), condition))
  .density_doit(expr, condition, evaluate = evaluate, ...)
}

#' @rdname E
#' @export
Density <- function(expr, condition = NULL, ...) {
  # sympy.physics.quantum.Density: the arguments are list(state, probability)
  if ((is.list(expr) && !is_sym(expr) && !inherits(expr, "symmatrix")) ||
      (is_sym(expr) && expr$type == "Tuple" && length(expr$args) == 2L && .qhas_node(expr$args[[1L]])))
    return(do.call(.qdensity, c(list(expr), if (!is.null(condition)) list(condition), list(...))))
  args <- list(if (inherits(expr, "symmatrix")) expr else as_sym(expr))
  if (!is.null(condition)) args[[2L]] <- .as_bool_sym(as_sym(condition))
  new_sym("Density", args)
}

register_doit("Density", function(e, deep = TRUE) .density_doit(e$args[[1L]], if (length(e$args) > 1L) e$args[[2L]] else NULL))

.density_doit <- function(expr, condition = NULL, evaluate = TRUE, ...) {
  if (is_sym(expr) && .isa(expr, "SingleFiniteDistribution")) return(.fdist_dict(expr))
  if (!is.null(condition)) expr <- given(expr, condition, ...)
  if (!length(random_symbols(expr)) && !(is_sym(expr) && !is.null(expr$extra$rv_matrix))) {
    x <- Symbol("x")
    return(Lambda(x, DiracDelta(.sub(x, as_sym(expr)))))
  }
  if (.is_rs(expr)) {
    ps <- .rs_pspace(expr)
    if ((.isa(ps, "SinglePSpace") || .isa(ps, "JointPSpace") || .isa(ps, "MatrixPSpace")) && !is.null(.findm(ps, "distribution")))
      return(.callm(ps, "distribution"))
  }
  if (is_sym(expr) && !is.null(expr$extra$rv_matrix) && .isa(expr$extra$pspace, "MatrixPSpace"))
    return(expr$extra$pspace$args[[2L]])
  if (is_sym(expr) && !is.null(expr$extra$rv_matrix) && .isa(expr$extra$pspace, "RandomMatrixPSpace"))
    return(expr$extra$pspace$args[[2L]])
  ps <- pspace(expr)
  result <- if (.isa(ps, "CompoundPSpace")) .callm(ps, "compute_density", expr, compound_evaluate = evaluate)
            else .callm(ps, "compute_density", expr)
  if (evaluate && is_sym(result) && !.isa(result, "Distribution")) return(doit(result))
  result
}

#' @rdname E
#' @export
rv_cdf <- function(expr, condition = NULL, evaluate = TRUE, ...) {
  if (!is.null(condition) && !.is_cond(condition)) return(apply_lambda(rv_cdf(expr), condition))
  if (!is.null(condition)) return(rv_cdf(given(expr, condition, ...), ...))
  result <- .callm(pspace(expr), "compute_cdf", expr)
  if (evaluate && is_sym(result)) return(doit(result))
  result
}

#' @rdname E
#' @export
characteristic_function <- function(expr, condition = NULL, evaluate = TRUE, ...) {
  if (!is.null(condition)) return(characteristic_function(given(expr, condition, ...), ...))
  result <- .callm(pspace(expr), "compute_characteristic_function", expr)
  if (evaluate && is_sym(result)) return(doit(result))
  result
}

#' @rdname E
#' @export
moment_generating_function <- function(expr, condition = NULL, evaluate = TRUE, ...) {
  if (!is.null(condition)) return(moment_generating_function(given(expr, condition, ...), ...))
  result <- .callm(pspace(expr), "compute_moment_generating_function", expr)
  if (evaluate && is_sym(result)) return(doit(result))
  result
}

#' @rdname E
#' @export
rv_quantile <- function(expr, evaluate = TRUE, ...) {
  result <- .callm(pspace(expr), "compute_quantile", expr)
  if (evaluate && is_sym(result)) return(doit(result))
  result
}

# ---------------------------------------------------------------------------
# Sampling. SymPy samples through scipy, numpy or pymc, whose generators R
# cannot reproduce; samples here use R's own generator (set.seed()).

#' Draw random samples
#'
#' `rv_sample(expr, condition, size)` is SymPy's `sample()`: one realization
#' of a random expression, or a vector of `size` of them. Samples use R's
#' random number generator; call `set.seed()` for reproducible draws.
#'
#' @param expr A random expression.
#' @param condition An optional condition the draws must satisfy; a number here is taken as `size`.
#' @param size Number of draws.
#' @param seed An optional seed passed to `set.seed()`.
#' @param ... Ignored.
#' @return A numeric vector.
#' @examples
#' set.seed(1)
#' rv_sample(Die("D", 6), size = 5)
#' @export
rv_sample <- function(expr, condition = NULL, size = 1, seed = NULL, ...) {
  # earlier versions took the size second: rv_sample(X, 10)
  if (is.numeric(condition)) { size <- condition; condition <- NULL }
  if (!is.null(seed)) set.seed(seed)
  rvs <- if (is.null(condition)) .ps_values(pspace(expr)) else .ps_values(pspace(Tuple(expr, condition)))
  rvs <- Filter(function(r) .has1(as_sym(expr), r) || (!is.null(condition) && .has1(as_sym(condition), r)), rvs)
  out <- numeric(0)
  while (length(out) < size) {
    draws <- lapply(rvs, function(rv) .rv_draw(rv, size))
    for (i in seq_len(size)) {
      vals <- lapply(draws, function(d) as_sym(d[i]))
      rep <- Map(function(r, v) list(r, v), rvs, vals)
      if (!is.null(condition)) {
        ok <- subs(as_sym(condition), rep)
        if (!(is_sym(ok) && ok$type == "BooleanAtom")) .serr("ValueError", "Conditions must not contain free symbols")
        if (!.is_true(ok)) next
      }
      out <- c(out, as.numeric(N(subs(as_sym(expr), rep))))
      if (length(out) >= size) break
    }
  }
  out[seq_len(size)]
}

.rv_draw <- function(rv, size) {
  ps <- .rs_pspace(rv)
  d <- .callm(ps, "distribution")
  f <- .findm(d, "rsample")
  if (!is.null(f)) return(f(d, size))
  if (.isa(d, "SingleFiniteDistribution")) {
    dd <- .fdist_items(d)
    vals <- vapply(dd, function(p) as.numeric(N(p[[1L]])), 0)
    probs <- vapply(dd, function(p) as.numeric(N(p[[2L]])), 0)
    return(vals[sample.int(length(vals), size, replace = TRUE, prob = probs)])
  }
  # inverse transform on a numeric grid of the cdf
  z <- Dummy("z", real = TRUE)
  cdf <- tryCatch(.dist_cdf(d, z), error = function(e) NULL)
  pdf <- .dist_pdf(d, z)
  s <- .dist_set(d)
  lo <- tryCatch(as.numeric(N(set_inf(s))), error = function(e) -Inf)
  hi <- tryCatch(as.numeric(N(set_sup(s))), error = function(e) Inf)
  if (.isa(d, "SingleDiscreteDistribution")) {
    a <- if (is.finite(lo)) lo else -200; b <- if (is.finite(hi)) hi else a + 2000
    k <- seq(a, b)
    f <- lambdify(z, pdf)
    w <- vapply(k, function(v) { r <- tryCatch(as.numeric(f(v)), error = function(e) 0); if (is.finite(r)) r else 0 }, 0)
    return(k[sample.int(length(k), size, replace = TRUE, prob = pmax(w, 0))])
  }
  # a grid 40 standard deviations wide around the mean, inside the support
  m <- tryCatch(as.numeric(N(E(rv))), error = function(e) NA)
  sd <- tryCatch(as.numeric(N(std(rv))), error = function(e) NA)
  if (length(m) != 1L || !is.finite(m)) m <- 0
  if (length(sd) != 1L || !is.finite(sd) || sd <= 0) sd <- 3
  a <- max(lo, m - 20 * sd); b <- min(hi, m + 20 * sd)
  grid <- seq(a, b, length.out = 20001L)
  f <- lambdify(z, pdf)
  w <- vapply(grid, function(v) { r <- tryCatch(as.numeric(f(v)), error = function(e) 0); if (length(r) == 1L && is.finite(r)) r else 0 }, 0)
  grid[sample.int(length(grid), size, replace = TRUE, prob = pmax(w, 0))]
}

# ---------------------------------------------------------------------------
# Attribute access and indexing for the battery and for users


#' Attributes and items of stats objects
#'
#' `rv_attr(x, name, ...)` reads a SymPy attribute of a random symbol,
#' probability space, domain or distribution (such as `"pspace"`,
#' `"domain"`, `"set"`, `"symbol"`, `"distribution"`, `"as_boolean"`), or
#' calls a method with further arguments. `rv_getitem(x, i)` indexes a
#' finite density, a joint random symbol or a stochastic process.
#'
#' @param x A stats object.
#' @param name The attribute or method name.
#' @param i An index or key.
#' @param ... Method arguments.
#' @return The attribute value.
#' @examples
#' X <- Normal("X", 0, 1)
#' rv_attr(rv_attr(X, "pspace"), "domain")
#' rv_attr(where(X^2 <= 1), "set")
#' @export
rv_attr <- function(x, name, ...) {
  if (inherits(x, "rv_process")) return(.process_attr(x, name, ...))
  if (name == "T" && (inherits(x, "symmatrix") || (is_sym(x) && x$type %in% .MATEXPR))) return(if (inherits(x, "symmatrix")) t(x) else Transpose(x))
  if (inherits(x, "symmatrix")) {
    if (name == "is_Matrix") return(TRUE)
    .serr("AttributeError", paste0("'Matrix' object has no attribute '", name, "'"))
  }
  if (is_sym(x) && !is.null(x$extra$rv_matrix)) {
    if (name == "pspace") return(x$extra$pspace)
    if (name == "symbol") return(x$extra$symbol)
    if (name == "args") return(Tuple(x$extra$symbol, x$args[[1L]], x$args[[2L]], x$extra$pspace))
  }
  if (.is_rs(x) && x$type == "Symbol" && name == "args") return(Tuple(x$extra$symbol, x$extra$pspace))
  if (.is_rs(x)) {
    switch(name, pspace = return(.rs_pspace(x)), symbol = return(.rs_symbol(x)), name = return(x$name),
           is_real = return(ask_fact(x, "real")), is_finite = return(ask_fact(x, "finite")),
           is_positive = return(ask_fact(x, "positive")), is_integer = return(ask_fact(x, "integer")),
           is_commutative = return(ask_fact(x, "commutative")), key = return(.ris_key(x)),
           NULL)
    if (name == "diff") return(diff(x, ...))
  }
  if (is_sym(x) && x$type == "FiniteDensity") {
    if (name == "dict") return(x)
    if (name == "items") return(x)
    if (name == "keys") return(lapply(x$args, function(t) t$args[[1L]]))
    if (name == "values") return(lapply(x$args, function(t) t$args[[2L]]))
  }
  # the pdf of a handmade distribution is its first argument, a Lambda
  if (is_sym(x) && name == "pdf" && !...length() && grepl("DistributionHandmade$", x$type)) return(x$args[[1L]])
  # the set of a finite distribution is a set of its values
  if (is_sym(x) && name == "set" && !...length() && .isa(x, "SingleFiniteDistribution")) {
    s <- .callm(x, "set")
    return(if (is_sym(s)) s else do.call(FiniteSet, s))
  }
  f <- if (is_sym(x)) .findm(x, name) else NULL
  if (!is.null(f)) return(f(x, ...))
  if (is_sym(x) && name == "equals") return(.expr_equals(x, ...))
  if (is_sym(x) && name %in% c("expr")) return(x$args[[1L]])
  if (is_sym(x) && name == "condition" && length(x$args) > 1L) return(x$args[[2L]])
  if (is_sym(x) && !is.null(.ARGN[[x$type]]) && name %in% .ARGN[[x$type]]) return(.da(x, name))
  if (is_sym(x) && x$type %in% c("Function", "Probability", "Expectation", "Variance", "Covariance", "Moment", "CentralMoment")) {
    if (name == "doit") return(doit(x))
    if (name == "evaluate_integral" || name == "evaluate_sum") return(doit(rewrite(x, "Integral")))
  }
  if (is_sym(x) && grepl("^is_", name)) {
    fl <- substring(name, 4)
    if (fl %in% c("Continuous", "Discrete", "Finite")) return(.ps_flag(x, name))
    return(ask_fact(x, fl))
  }
  .serr("AttributeError", paste0("'", if (is_sym(x)) .sympy_class(x) else class(x)[1], "' object has no attribute '", name, "'"))
}

#' @rdname rv_attr
#' @export
rv_getitem <- function(x, i, ...) {
  if (inherits(x, "rv_process")) return(.process_getitem(x, i))
  if (inherits(x, "symmatrix")) {
    idx <- c(list(i), list(...))
    if (length(idx) == 2L) return(.mget(x, as.integer(idx[[1L]]) + 1L, as.integer(idx[[2L]]) + 1L))
    return(x$data[[as.integer(i) + 1L]])
  }
  if (is.list(x) && !is_sym(x)) return(x[[as.integer(i) + 1L]])
  if (.is_rs(x) && .isa(.rs_pspace(x), "JointPSpace")) return(.joint_getitem(x, i))
  if (is_sym(x) && x$type %in% c("FiniteDensity", "Dict")) {
    k <- as_sym(i)
    for (t in x$args) if (identical(t$args[[1L]]$key, k$key)) return(t$args[[2L]])
    .serr("KeyError", .str(k))
  }
  if (is_sym(x) && .isa(x, "SingleFiniteDistribution")) return(rv_getitem(.fdist_dict(x), i))
  if (is_sym(x) && x$type %in% c("Tuple", "FiniteSet")) return(x$args[[as.integer(i) + 1L]])
  .serr("TypeError", paste0("'", if (is_sym(x)) .sympy_class(x) else class(x)[1], "' object is not subscriptable"))
}

# Calling a distribution, a finite density or a Lambda: SymPy's f(x)
.rv_call <- function(f, args) {
  if (f$type == "FiniteDensity") {
    k <- as_sym(args[[1L]])
    for (t in f$args) if (identical(t$args[[1L]]$key, k$key)) return(t$args[[2L]])
    return(S_Zero)
  }
  conv <- function(a) {
    # an IndexedBase is an R function carrying its node
    if (is.function(a) && !is.null(attr(a, "base"))) return(attr(a, "base"))
    if (inherits(a, "symmatrix") || (is.list(a) && !is_sym(a))) a else as_sym(a)
  }
  if (!is.null(.findm(f, "call"))) return(do.call(.callm, c(list(f, "call"), lapply(args, conv))))
  .serr("TypeError", paste0("'", .sympy_class(f), "' object is not callable"))
}

.defm("Distribution", "call", function(d, ...) .callm(d, "pdf", ...))
