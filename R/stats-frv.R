# Finite probability spaces, after sympy/stats/frv.py and frv_types.py.

.defclass("FiniteDomain", "RandomDomain")
.defclass("SingleFiniteDomain", "FiniteDomain")
.defclass("ProductFiniteDomain", c("ProductDomain", "FiniteDomain"))
.defclass("ConditionalFiniteDomain", c("ConditionalDomain", "ProductFiniteDomain"))
.defclass("SingleFiniteDistribution", "Distribution")
.defclass("FinitePSpace", "PSpace")
.defclass("SingleFinitePSpace", c("SinglePSpace", "FinitePSpace"))
.defclass("ProductFinitePSpace", c("IndependentProductPSpace", "FinitePSpace"))

.defm("FiniteDomain", "is_Finite", function(d) TRUE)
.defm("FinitePSpace", "is_Finite", function(ps) TRUE)

# ---------------------------------------------------------------------------
# Python containers emulated here
#
# A finite density is a dict in SymPy; FiniteDensity keeps its insertion
# order and prints like a Python dict.

.FD <- function(keys = list(), vals = list()) new_sym("FiniteDensity", Map(function(k, v) Tuple(k, v), keys, vals))
.fd_keys <- function(d) lapply(d$args, function(t) t$args[[1L]])
.fd_vals <- function(d) lapply(d$args, function(t) t$args[[2L]])

register_str_printer("FiniteDensity", function(e)
  paste0("{", paste(vapply(e$args, function(t) paste0(.str(t$args[[1L]]), ": ", .str(t$args[[2L]])), ""), collapse = ", "), "}"))

# An accumulating dict: d[key] = d.get(key, 0) + prob
.acc_new <- function() { e <- new.env(); e$keys <- list(); e$vals <- list(); e$idx <- character(); e }
.acc_add <- function(acc, key, val) {
  i <- match(key$key, acc$idx)
  if (is.na(i)) {
    acc$idx <- c(acc$idx, key$key)
    acc$keys[[length(acc$keys) + 1L]] <- key
    acc$vals[[length(acc$vals) + 1L]] <- val
  } else acc$vals[[i]] <- .add(list(acc$vals[[i]], val))
  invisible(acc)
}
.acc_set <- function(acc, key, val) {
  i <- match(key$key, acc$idx)
  if (is.na(i)) .acc_add(acc, key, val) else acc$vals[[i]] <- val
  invisible(acc)
}
.acc_fd <- function(acc) .FD(acc$keys, acc$vals)

# Iteration order of a Python set of integers: CPython's open addressing
# table with hash(n) = n (and hash(-1) = -2). Other elements, whose hashes
# vary between Python processes, keep SymPy's sort order.
.pyset_order <- function(items) {
  items <- .uniq(items)
  if (!length(items) || !all(vapply(items, function(z) z$type == "Integer", NA))) return(.sorted_basic(items))
  vals <- vapply(items, function(z) as.numeric(z$value), 0)
  if (any(abs(vals) > 1e9)) return(.sorted_basic(items))
  h <- ifelse(vals == -1, -2, vals)
  mask <- 7; table <- rep(NA_integer_, 8L); fill <- 0
  insert <- function(table, mask, j) {
    hh <- h[j]
    # perturb is the hash as an unsigned 64 bit integer, kept as two 32 bit
    # halves; only the low half reaches the index
    lo <- hh %% 2^32; hi <- if (hh < 0) 2^32 - 1 else 0
    i <- hh %% (mask + 1)
    repeat {
      probes <- if (i + 9 <= mask) 9 else 0
      for (q in 0:probes) {
        if (is.na(table[i + q + 1L])) { table[i + q + 1L] <- j; return(table) }
      }
      lo <- floor(lo / 32) + (hi %% 32) * 2^27; hi <- floor(hi / 32)
      i <- (i * 5 + 1 + lo) %% (mask + 1)
    }
  }
  for (j in seq_along(items)) {
    table <- insert(table, mask, j)
    fill <- fill + 1
    if (fill * 5 >= mask * 3) {
      used <- fill
      newsize <- 8
      while (newsize <= used * 4) newsize <- newsize * 2
      old <- table[!is.na(table)]
      mask <- newsize - 1; table <- rep(NA_integer_, newsize)
      for (k in old) table <- insert(table, mask, k)
    }
  }
  items[table[!is.na(table)]]
}

# ---------------------------------------------------------------------------
# Elements of finite domains: frozensets of (symbol, value) pairs, kept as
# lists of pairs sorted by symbol

.elem_sort <- function(pairs) pairs[order(vapply(pairs, function(p) p[[1L]]$name, ""), method = "radix")]
.elem_union <- function(elems) .elem_sort(unlist(elems, recursive = FALSE))
.elem_key <- function(el) paste(vapply(el, function(p) paste0(p[[1L]]$key, "=", p[[2L]]$key), ""), collapse = "&")
.elem_sym <- function(el) new_sym("FiniteSet", .sorted_basic(lapply(el, function(p) Tuple(p[[1L]], p[[2L]]))))

.fdom_iter <- function(d) .callm(d, "iter")

.defm("FiniteDomain", "elements", function(d) d$args[[1L]])
.defm("FiniteDomain", "symbols", function(d) .fset_basic(lapply(.fdom_iter(d)[[1L]], function(p) p[[1L]])))
# a FiniteSet of Dicts, in SymPy's order
.defm("FiniteDomain", "dict", function(d) new_sym("FiniteSet", .fset_basic(lapply(.fdom_iter(d), function(el) Dict(lapply(el, `[[`, 1L), lapply(el, `[[`, 2L))))))
.defm("FiniteDomain", "as_boolean", function(d)
  do.call(Or, lapply(.fdom_iter(d), function(el) do.call(And, lapply(el, function(p) Eq(p[[1L]], p[[2L]]))))))
.defm("FiniteDomain", "contains", function(d, el) .elem_key(el) %in% vapply(.fdom_iter(d), .elem_key, ""))

.SingleFiniteDomain <- function(symbol, set) {
  if (!(is_sym(set) && set$type %in% c("FiniteSet", "Intersection"))) set <- if (length(set)) do.call(FiniteSet, set) else S_EmptySet
  new_sym("SingleFiniteDomain", list(symbol, set))
}
.defm("SingleFiniteDomain", "symbol", function(d) d$args[[1L]])
.defm("SingleFiniteDomain", "symbols", function(d) list(d$args[[1L]]))
.defm("SingleFiniteDomain", "set", function(d) d$args[[2L]])
.defm("SingleFiniteDomain", "iter", function(d) {
  s <- d$args[[2L]]
  if (s$type != "FiniteSet") .serr("TypeError", "'Intersection' object is not iterable")
  lapply(s$args, function(v) list(list(d$args[[1L]], v)))
})
.defm("SingleFiniteDomain", "elements", function(d) new_sym("FiniteSet", .sorted_basic(lapply(.fdom_iter(d), .elem_sym))))

.defm("ProductFiniteDomain", "iter", function(d) {
  out <- list(list())
  for (dom in d$args) {
    els <- .fdom_iter(dom)
    nxt <- list()
    for (a in out) for (b in els) nxt[[length(nxt) + 1L]] <- c(a, b)
    out <- nxt
  }
  lapply(out, .elem_sort)
})
.defm("ProductFiniteDomain", "elements", function(d) new_sym("FiniteSet", .fset_basic(lapply(.fdom_iter(d), .elem_sym))))
.defm("ProductFiniteDomain", "as_boolean", function(d)
  do.call(Or, lapply(.fdom_iter(d), function(el) do.call(And, lapply(el, function(p) Eq(p[[1L]], p[[2L]]))))))

.ConditionalFiniteDomain <- function(domain, condition) {
  if (isTRUE(condition) || .is_true(as_sym(condition))) return(domain)
  new_sym("ConditionalFiniteDomain", list(domain, rv_subs(.as_bool_sym(as_sym(condition)))))
}
.cfd_test <- function(d, el) {
  val <- xreplace(d$args[[2L]], el)
  if (is.logical(val)) return(val)
  if (val$type == "BooleanAtom") return(val$name == "true")
  if (val$type == "Equality") return(identical(val$args[[1L]]$key, val$args[[2L]]$key))
  .serr("ValueError", paste0("Undecidable if ", .str(val)))
}
.defm("ConditionalFiniteDomain", "iter", function(d) Filter(function(el) .cfd_test(d, el), .fdom_iter(d$args[[1L]])))
.defm("ConditionalFiniteDomain", "test", function(d, el) .cfd_test(d, el))
.defm("ConditionalFiniteDomain", "symbols", function(d) .dom_symbols(d$args[[1L]]))
.defm("ConditionalFiniteDomain", "set", function(d) {
  fd <- d$args[[1L]]
  if (fd$type == "SingleFiniteDomain") {
    s <- fd$args[[1L]]
    keep <- Filter(function(v) .cfd_test(d, list(list(s, v))) , .finite_elems(fd$args[[2L]]))
    return(do.call(FiniteSet, keep))
  }
  .serr("NotImplementedError", "Not implemented on multi-dimensional conditional domain")
})
.defm("ConditionalFiniteDomain", "as_boolean", function(d) .STM[["FiniteDomain"]]$as_boolean(d))
.defm("ConditionalFiniteDomain", "elements", function(d) new_sym("FiniteSet", .fset_basic(lapply(.fdom_iter(d), .elem_sym))))

# free symbols of a finite domain (Basic.free_symbols over its args)
.fdom_free <- function(d) {
  switch(d$type,
         SingleFiniteDomain = .uniq(c(list(d$args[[1L]]), free_symbols(d$args[[2L]]))),
         ConditionalFiniteDomain = .uniq(c(.fdom_free(d$args[[1L]]), free_symbols(d$args[[2L]]))),
         .uniq(unlist(lapply(d$args, .fdom_free), recursive = FALSE)))
}

# ---------------------------------------------------------------------------
# Finite distributions

.fdist_symbolic <- function(d) {
  f <- .findm(d, "is_symbolic")
  if (is.null(f)) FALSE else isTRUE(f(d))
}

# SingleFiniteDistribution.dict as a FiniteDensity (or Density when symbolic)
.fdist_dict <- function(d) {
  f <- .findm(d, "dict")
  f(d)
}
.fdist_items <- function(d) {
  dd <- .fdist_dict(d)
  lapply(dd$args, function(t) list(t$args[[1L]], t$args[[2L]]))
}

.defm("SingleFiniteDistribution", "dict", function(d) {
  if (.fdist_symbolic(d)) return(Density(d))
  s <- .callm(d, "set")
  ks <- if (is.list(s) && !is_sym(s)) s else .finite_elems(s)
  .FD(ks, lapply(ks, function(k) .callm(d, "pmf", k)))
})
.defm("SingleFiniteDistribution", "pmf", function(d, ...) .serr("NotImplementedError", ""))
.defm("SingleFiniteDistribution", "call", function(d, ...) {
  pm <- .callm(d, "pmf", ...)
  pm
})
.defm("SingleFiniteDistribution", "values", function(d) .fd_vals(.fdist_dict(d)))
.defm("SingleFiniteDistribution", "items", function(d) .fdist_dict(d))
.defm("SingleFiniteDistribution", "is_symbolic", function(d) FALSE)

.defdfin <- function(.cls, .argnames, ...) .defdist(.cls, "SingleFiniteDistribution", .argnames, ...)

.rvf <- function(name, cls, args, check = TRUE) .rv_make(name, cls, args, check = check, single = "SingleFinitePSpace")

# ---------------------------------------------------------------------------
# FinitePSpace

.FinitePSpace <- function(domain, keys, vals) {
  new_sym("FinitePSpace", list(domain, Dict(keys, vals)), extra = list(dk = keys, dv = vals))
}

# (keys, values) of the private density
.fps_density <- function(ps) .callm(ps, "_density")
.defm("FinitePSpace", "_density", function(ps) list(keys = ps$extra$dk, vals = ps$extra$dv))

.fps_prob_of <- function(ps, el, dens = .fps_density(ps)) {
  if (!length(dens$keys)) return(S_Zero)
  if (is_sym(dens$keys[[1L]]) && dens$keys[[1L]]$type == "FiniteSet") {
    k <- .elem_sym(el)$key
  } else k <- el[[1L]][[2L]]$key
  i <- match(k, vapply(dens$keys, .k, ""))
  if (is.na(i)) S_Zero else dens$vals[[i]]
}

.defm("FinitePSpace", "where", function(ps, condition) .ConditionalFiniteDomain(.ps_domain(ps), condition))

.defm("FinitePSpace", "compute_density", function(ps, expr, ...) {
  expr <- rv_subs(as_sym(expr), .ps_values(ps))
  acc <- .acc_new()
  dens <- .fps_density(ps)
  for (el in .fdom_iter(.ps_domain(ps))) {
    val <- xreplace(expr, el)
    .acc_add(acc, as_sym(val), .fps_prob_of(ps, el, dens))
  }
  .acc_fd(acc)
})

# Python's sorted() on SymPy keys: numbers compare, anything else raises
.py_sorted <- function(keys) {
  if (length(keys) < 2L) return(keys)
  num <- vapply(keys, function(k) is_sym(k) && isTRUE(ask_fact(k, "extended_real")) && .rv_is_number(k), NA)
  if (!all(num)) {
    if (all(vapply(keys, function(k) k$type == "BooleanAtom", NA)))
      return(keys[order(vapply(keys, function(k) k$name == "true", NA))])
    .serr("TypeError", "cannot determine truth value of Relational")
  }
  keys[order(vapply(keys, function(k) as.numeric(N(k)), 0), method = "radix")]
}

.defm("FinitePSpace", "compute_cdf", function(ps, expr, ...) {
  d <- .callm(ps, "compute_density", expr)
  keys <- .fd_keys(d); vals <- .fd_vals(d)
  sk <- .py_sorted(keys)
  cum <- S_Zero; outv <- list()
  for (k in sk) {
    cum <- .add(list(cum, vals[[match(k$key, vapply(keys, .k, ""))]]))
    outv[[length(outv) + 1L]] <- cum
  }
  .FD(sk, outv)
})

.defm("FinitePSpace", "compute_characteristic_function", function(ps, expr, ...) {
  d <- .callm(ps, "compute_density", expr)
  t <- Dummy("t", real = TRUE)
  Lambda(t, .add(c(list(S_Zero), lapply(d$args, function(kv) .mul(list(exp(.mul(list(S_I, kv$args[[1L]], t))), kv$args[[2L]]))))))
})
.defm("FinitePSpace", "compute_moment_generating_function", function(ps, expr, ...) {
  d <- .callm(ps, "compute_density", expr)
  t <- Dummy("t", real = TRUE)
  Lambda(t, .add(c(list(S_Zero), lapply(d$args, function(kv) .mul(list(exp(.mul(list(kv$args[[1L]], t))), kv$args[[2L]]))))))
})

.is_logic_rel <- function(e) is.logical(e) || (is_sym(e) && (is_relational(e) || e$type %in% c("And", "Or", "Not", "BooleanAtom", "Xor", "Implies", "Equivalent", "ITE")))

.defm("FinitePSpace", "compute_expectation", function(ps, expr, rvs = NULL, ...) {
  if (is.null(rvs) || !length(rvs)) rvs <- .ps_values(ps)
  expr <- rv_subs(as_sym(expr), rvs)
  els <- .fdom_iter(.ps_domain(ps))
  dens <- .fps_density(ps)
  probs <- lapply(els, function(el) .fps_prob_of(ps, el, dens))
  if (.is_logic_rel(expr)) {
    parse <- lapply(els, function(el) el[[1L]][[2L]])
    bools <- lapply(els, function(el) xreplace(expr, el))
  } else {
    parse <- lapply(els, function(el) xreplace(expr, el))
    bools <- rep(list(S_true), length(els))
  }
  terms <- Map(function(p, e, b) Piecewise(list(.mul(list(p, e)), b), list(S_Zero, S_true)), probs, parse, bools)
  .add(c(list(S_Zero), terms))
})

.defm("FinitePSpace", "compute_quantile", function(ps, expr, ...) {
  cdf <- .callm(ps, "compute_cdf", expr)
  p <- Dummy("p", real = TRUE)
  pieces <- list(list(S_NaN, Or(Lt(p, 0), Gt(p, 1))))
  for (kv in cdf$args) pieces[[length(pieces) + 1L]] <- list(kv$args[[1L]], Le(p, kv$args[[2L]]))
  Lambda(p, do.call(Piecewise, pieces))
})

.defm("FinitePSpace", "probability", function(ps, condition, ...) {
  condition <- .as_bool_sym(as_sym(condition))
  cond_symbols <- lapply(random_symbols(condition), .rs_symbol)
  cond <- rv_subs(condition)
  pss <- .ps_symbols(ps)
  bad <- Filter(function(s) !.in_list(s, pss), cond_symbols)
  if (length(bad)) .serr("ValueError", paste0("Cannot compare foreign random symbols, frozenset({", paste(vapply(bad, .str, ""), collapse = ", "), "})"))
  if (is_relational(condition) && !all(vapply(free_symbols(cond), function(s) .in_list(s, .fdom_free(.ps_domain(ps))), NA))) {
    rv <- if (condition$args[[2L]]$type == "Symbol" && is.null(condition$args[[2L]]$extra$rv_class)) condition$args[[1L]] else condition$args[[2L]]
    dens <- .fps_density(ps)
    terms <- lapply(.fdom_iter(.ps_domain(ps)), function(el)
      Piecewise(list(.fps_prob_of(ps, el, dens), subs(condition, rv, el[[1L]][[2L]])), list(S_Zero, S_true)))
    return(.add(c(list(S_Zero), terms)))
  }
  dens <- .fps_density(ps)
  .add(c(list(S_Zero), lapply(.fdom_iter(.callm(ps, "where", condition)), function(el) .fps_prob_of(ps, el, dens))))
})

.fps_conditional <- function(ps, condition) {
  domain <- .callm(ps, "where", condition)
  prob <- .callm(ps, "probability", condition)
  dens <- .fps_density(ps)
  keys <- list(); vals <- list()
  for (i in seq_along(dens$keys)) {
    k <- dens$keys[[i]]
    el <- lapply(k$args, function(t) list(t$args[[1L]], t$args[[2L]]))
    if (isTRUE(.callm(domain, "test", el))) {
      keys[[length(keys) + 1L]] <- k
      vals[[length(vals) + 1L]] <- .div(dens$vals[[i]], prob)
    }
  }
  .FinitePSpace(domain, keys, vals)
}
.defm("FinitePSpace", "conditional_space", function(ps, condition, ...) .fps_conditional(ps, condition))
.defm("FinitePSpace", "test", function(d, el) TRUE)
.defm("FiniteDomain", "test", function(d, el) TRUE)

# ---------------------------------------------------------------------------
# SingleFinitePSpace

.sfps_dist <- function(ps) ps$args[[2L]]
.sfps_symbolic <- function(ps) .fdist_symbolic(ps$args[[2L]])
.sfps_as_fps <- function(ps) {
  items <- .fdist_items(ps$args[[2L]])
  .FinitePSpace(.ps_domain(ps), lapply(items, `[[`, 1L), lapply(items, `[[`, 2L))
}

.defm("SingleFinitePSpace", "domain", function(ps) .SingleFiniteDomain(ps$args[[1L]], .callm(ps$args[[2L]], "set")))
.defm("SingleFinitePSpace", "density", function(ps) ps$args[[2L]])
.defm("SingleFinitePSpace", "pmf", function(ps, expr) .callm(ps$args[[2L]], "pmf", expr))
.defm("SingleFinitePSpace", "_density", function(ps) {
  items <- .fdist_items(ps$args[[2L]])
  list(keys = lapply(items, function(kv) .elem_sym(list(list(ps$args[[1L]], kv[[1L]])))), vals = lapply(items, `[[`, 2L))
})
.defm("SingleFinitePSpace", "compute_characteristic_function", function(ps, expr, ...) {
  if (.sfps_symbolic(ps)) {
    d <- .callm(ps, "compute_density", expr)
    t <- Dummy("t", real = TRUE); ki <- Dummy("ki")
    dist <- ps$args[[2L]]
    return(Lambda(t, Sum(.rv_lcall(d, ki) * exp(S_I * ki * t), list(ki, .callm(dist, "low"), .callm(dist, "high")))))
  }
  expr <- rv_subs(as_sym(expr), .ps_values(ps))
  .STM[["FinitePSpace"]]$compute_characteristic_function(.sfps_as_fps(ps), expr)
})
.defm("SingleFinitePSpace", "compute_moment_generating_function", function(ps, expr, ...) {
  if (.sfps_symbolic(ps)) {
    d <- .callm(ps, "compute_density", expr)
    t <- Dummy("t", real = TRUE); ki <- Dummy("ki")
    dist <- ps$args[[2L]]
    return(Lambda(t, Sum(.rv_lcall(d, ki) * exp(ki * t), list(ki, .callm(dist, "low"), .callm(dist, "high")))))
  }
  expr <- rv_subs(as_sym(expr), .ps_values(ps))
  .STM[["FinitePSpace"]]$compute_moment_generating_function(.sfps_as_fps(ps), expr)
})
.defm("SingleFinitePSpace", "compute_quantile", function(ps, expr, ...) {
  if (.sfps_symbolic(ps))
    .serr("NotImplementedError", "Computing quantile for random variables with symbolic dimension because the bounds of searching the required value is undetermined.")
  expr <- rv_subs(as_sym(expr), .ps_values(ps))
  .STM[["FinitePSpace"]]$compute_quantile(.sfps_as_fps(ps), expr)
})
.defm("SingleFinitePSpace", "compute_density", function(ps, expr, ...) {
  if (.sfps_symbolic(ps)) {
    rv <- random_symbols(expr)[[1L]]
    k <- Dummy("k", integer = TRUE)
    cond <- if (!.is_logic_rel(expr)) S_true else subs(expr, rv, k)
    dist <- ps$args[[2L]]
    return(Lambda(k, Piecewise(list(.callm(dist, "pmf", k), And(Ge(k, .callm(dist, "low")), Le(k, .callm(dist, "high")), cond)), list(S_Zero, S_true))))
  }
  expr <- rv_subs(as_sym(expr), .ps_values(ps))
  .STM[["FinitePSpace"]]$compute_density(.sfps_as_fps(ps), expr)
})
.defm("SingleFinitePSpace", "compute_cdf", function(ps, expr, ...) {
  if (.sfps_symbolic(ps)) {
    d <- .callm(ps, "compute_density", expr)
    k <- Dummy("k"); ki <- Dummy("ki")
    return(Lambda(k, Sum(.rv_lcall(d, ki), list(ki, .callm(ps$args[[2L]], "low"), k))))
  }
  expr <- rv_subs(as_sym(expr), .ps_values(ps))
  .STM[["FinitePSpace"]]$compute_cdf(.sfps_as_fps(ps), expr)
})
.defm("SingleFinitePSpace", "compute_expectation", function(ps, expr, rvs = NULL, ...) {
  if (.sfps_symbolic(ps)) {
    rv <- random_symbols(expr)[[1L]]
    k <- Dummy("k", integer = TRUE)
    expr <- subs(as_sym(expr), rv, k)
    cond <- if (!.is_logic_rel(expr)) S_true else expr
    dist <- ps$args[[2L]]
    func <- if (!.is_true(.as_bool_sym(cond))) .mul(list(.callm(dist, "pmf", k), k)) else .mul(list(.callm(dist, "pmf", k), expr))
    return(doit(Sum(Piecewise(list(func, cond), list(S_Zero, S_true)), list(k, .callm(dist, "low"), .callm(dist, "high")))))
  }
  expr <- rv_subs(as_sym(expr), rvs)
  .STM[["FinitePSpace"]]$compute_expectation(.sfps_as_fps(ps), expr, rvs, ...)
})
.defm("SingleFinitePSpace", "probability", function(ps, condition, ...) {
  if (.sfps_symbolic(ps))
    .serr("NotImplementedError", "Currently, probability queries are not supported for random variables with symbolic sized distributions.")
  condition <- rv_subs(.as_bool_sym(as_sym(condition)))
  .STM[["FinitePSpace"]]$probability(.sfps_as_fps(ps), condition)
})
.defm("SingleFinitePSpace", "conditional_space", function(ps, condition, ...) {
  if (.sfps_symbolic(ps))
    .serr("NotImplementedError", "Conditional space is not supported for random variables with symbolic sized distributions.")
  .fps_conditional(ps, condition)
})

# ---------------------------------------------------------------------------
# ProductFinitePSpace

.defm("ProductFinitePSpace", "domain", function(ps) do.call(.ProductDomain, lapply(ps$args, .ps_domain)))
.defm("ProductFinitePSpace", "_density", function(ps) {
  acc <- .acc_new()
  combos <- list(list(elems = list(), probs = list()))
  for (sp in ps$args) {
    dd <- .fps_density(sp)
    nxt <- list()
    for (cb in combos) for (i in seq_along(dd$keys))
      nxt[[length(nxt) + 1L]] <- list(elems = c(cb$elems, list(dd$keys[[i]])), probs = c(cb$probs, list(dd$vals[[i]])))
    combos <- nxt
  }
  for (cb in combos) {
    el <- .elem_union(lapply(cb$elems, function(k) lapply(k$args, function(t) list(t$args[[1L]], t$args[[2L]]))))
    .acc_add(acc, .elem_sym(el), .mul(cb$probs))
  }
  list(keys = acc$keys, vals = acc$vals)
})
.defm("ProductFinitePSpace", "density", function(ps) { d <- .fps_density(ps); Dict(d$keys, d$vals) })
.defm("ProductFinitePSpace", "probability", function(ps, condition, ...) .STM[["FinitePSpace"]]$probability(ps, condition))
.defm("ProductFinitePSpace", "compute_density", function(ps, expr, ...) .STM[["FinitePSpace"]]$compute_density(ps, expr))
.defm("ProductFinitePSpace", "compute_expectation", function(ps, expr, rvs = NULL, ...) .STM[["IndependentProductPSpace"]]$compute_expectation(ps, expr, rvs, ...))
.defm("ProductFinitePSpace", "conditional_space", function(ps, condition, ...) .fps_conditional(ps, condition))
.defm("ProductFinitePSpace", "where", function(ps, condition) .ConditionalFiniteDomain(.ps_domain(ps), condition))

# ---------------------------------------------------------------------------
# Finite distribution types

#' Finite random variables
#'
#' Constructors for the finite distributions of `sympy.stats`: `FiniteRV()`
#' from a mapping of values to probabilities, `DiscreteUniform()`, `Die()`,
#' `Bernoulli()`, `Coin()`, `Binomial()`, `BetaBinomial()`,
#' `Hypergeometric()`, `Rademacher()`, `IdealSoliton()` and
#' `RobustSoliton()`.
#'
#' @param name The name of the random variable.
#' @param density A `Dict()` (or a list of values with names) mapping
#'   values to probabilities.
#' @param items A list of equally likely values.
#' @param sides Number of sides of the die.
#' @param p,succ,fail Success probability and the success and failure
#'   values.
#' @param n,alpha,beta,N,m,k,delta,c Distribution parameters, as in SymPy.
#' @param check Whether to check the arguments.
#' @return A random symbol.
#' @examples
#' D <- Die("D", 6)
#' E(D)
#' P(D > 3)
#' rv_density(Coin("C"))
#' @name finite_rv
NULL

.defdfin("FiniteDistributionHandmade", character(),
  dict = function(d) d$args[[1L]],
  pmf = function(d, x) {
    x <- Symbol("x")
    dd <- d$args[[1L]]
    Lambda(x, do.call(Piecewise, c(lapply(dd$args, function(t) list(t$args[[2L]], Eq(t$args[[1L]], x))), list(list(S_Zero, S_true)))))
  },
  set = function(d) .pyset_order(.fd_keys(d$args[[1L]])),
  check = function(density) {
    for (p in .fd_vals(density)) .value_check(list(p >= 0, p <= 1), "Probability at a point must be between 0 and 1.")
    val <- .add(c(list(S_Zero), .fd_vals(density)))
    .value_check(!identical(.fuzzy_bool(Eq(val, 1)), FALSE), "Total Probability must be 1.")
  })

.as_fd <- function(density) {
  if (is_sym(density) && density$type == "FiniteDensity") return(density)
  if (is_sym(density) && density$type == "Dict") return(.FD(lapply(density$args, function(t) t$args[[1L]]), lapply(density$args, function(t) t$args[[2L]])))
  if (is.list(density) && !is.null(names(density))) return(.FD(lapply(names(density), function(n) as_sym(type.convert(n, as.is = TRUE))), lapply(density, as_sym)))
  stop("density must be a Dict or a named list", call. = FALSE)
}

#' @rdname finite_rv
#' @export
FiniteRV <- function(name, density, check = FALSE) .rvf(name, "FiniteDistributionHandmade", list(.as_fd(density)), check = check)

.defdfin("DiscreteUniformDistribution", character(),
  check = function(...) {
    a <- list(...)
    if (length(.uniq(a)) != length(a)) {
      ks <- .uniq(a)
      n <- length(a)
      w <- vapply(ks, function(k) sum(vapply(a, function(z) identical(z$key, k$key), NA)), 0)
      ws <- paste0("{", paste(vapply(seq_along(ks), function(i) paste0(.str(ks[[i]]), ": ", .str(Rational(w[i], n))), ""), collapse = ", "), "}")
      .serr("ValueError", paste0("\nRepeated args detected but set expected. For a distribution having\ndifferent weights for each item use the following:\nS(\"FiniteRV('X', ", ws, ")\")"))
    }
  },
  p = function(d) Rational(1, length(d$args)),
  dict = function(d) { ks <- .callm(d, "set"); .FD(ks, rep(list(Rational(1, length(d$args))), length(ks))) },
  set = function(d) .pyset_order(d$args),
  pmf = function(d, x) if (.in_list(as_sym(x), d$args)) Rational(1, length(d$args)) else S_Zero)

#' @rdname finite_rv
#' @export
DiscreteUniform <- function(name, items) .rvf(name, "DiscreteUniformDistribution", lapply(items, as_sym))

.pmf_arg_check <- function(x) {
  # explicit matrices and undefined functions are not expressions here
  cls <- if (inherits(x, "symmatrix")) "sympy.matrices.immutable.ImmutableDenseMatrix" else
    if (inherits(x, "sym_function")) "sympy.core.function.UndefinedFunction" else NULL
  if (is.null(cls)) {
    x <- as_sym(x)
    if (isTRUE(.rv_is_number(x)) || x$type %in% c("Symbol", "Dummy") || is_random(x)) return(x)
    cls <- paste0("sympy.core.", tolower(x$type), ".", x$type)
  }
  .serr("ValueError", paste0("'x' expected as an argument of type 'number', 'Symbol', or 'RandomSymbol' not <class '", cls, "'>"))
}
.rv_is_number <- function(x) length(free_symbols(x)) == 0L && !is_random(x)

.defdfin("DieDistribution", c("sides"),
  check = function(sides) .value_check(list(.st_fact(sides, "positive"), .st_fact(sides, "integer")), "number of sides must be a positive integer."),
  is_symbolic = function(d) !.rv_is_number(.da(d, "sides")),
  high = function(d) .da(d, "sides"),
  low = function(d) S_One,
  set = function(d) {
    if (.fdist_symbolic(d)) return(Intersection(Naturals0, Interval(0, .da(d, "sides"))))
    .pyset_order(lapply(seq_len(as.integer(.da(d, "sides")$value)), as_sym))
  },
  pmf = function(d, x) {
    x <- .pmf_arg_check(x)
    cond <- And(Ge(x, 1), Le(x, .da(d, "sides")), Contains(x, Integers))
    Piecewise(list(S_One / .da(d, "sides"), cond), list(S_Zero, S_true))
  })

#' @rdname finite_rv
#' @export
Die <- function(name, sides = 6) .rvf(name, "DieDistribution", list(sides))

.defdfin("BernoulliDistribution", c("p", "succ", "fail"),
  check = function(p, succ, fail) .value_check(list(p >= 0, p <= 1), "p should be in range [0, 1]."),
  set = function(d) .pyset_order(list(.da(d, "succ"), .da(d, "fail"))),
  pmf = function(d, x) {
    p <- .dl(d)
    x <- as_sym(x)
    if (p$succ$type == "Symbol" && p$fail$type == "Symbol")
      return(Piecewise(list(p$p, as_sym(identical(x$key, p$succ$key))), list(1 - p$p, as_sym(identical(x$key, p$fail$key))), list(S_Zero, S_true)))
    Piecewise(list(p$p, Eq(x, p$succ)), list(1 - p$p, Eq(x, p$fail)), list(S_Zero, S_true))
  })

#' @rdname finite_rv
#' @export
Bernoulli <- function(name, p, succ = 1, fail = 0) .rvf(name, "BernoulliDistribution", list(p, succ, fail))

#' @rdname finite_rv
#' @export
Coin <- function(name, p = S_Half) .rvf(name, "BernoulliDistribution", list(p, Symbol("H"), Symbol("T")))

.defdfin("BinomialDistribution", c("n", "p", "succ", "fail"),
  check = function(n, p, succ, fail) {
    .value_check(list(.st_fact(n, "integer"), .st_fact(n, "nonnegative")), "'n' must be nonnegative integer.")
    .value_check(list(p <= 1, p >= 0), "p should be in range [0, 1].")
  },
  high = function(d) .da(d, "n"),
  low = function(d) S_Zero,
  is_symbolic = function(d) !.rv_is_number(.da(d, "n")),
  set = function(d) {
    if (.fdist_symbolic(d)) return(Intersection(Naturals0, Interval(0, .da(d, "n"))))
    .pyset_order(.fd_keys(.fdist_dict(d)))
  },
  pmf = function(d, x) {
    p <- .dl(d)
    x <- .pmf_arg_check(x)
    cond <- And(Ge(x, 0), Le(x, p$n), Contains(x, Integers))
    Piecewise(list(sym_binomial(p$n, x) * p$p^x * (1 - p$p)^(p$n - x), cond), list(S_Zero, S_true))
  },
  dict = function(d) {
    if (.fdist_symbolic(d)) return(Density(d))
    p <- .dl(d)
    acc <- .acc_new()
    for (k in 0:as.integer(p$n$value)) .acc_set(acc, k * p$succ + (p$n - k) * p$fail, .callm(d, "pmf", k))
    .acc_fd(acc)
  })

#' @rdname finite_rv
#' @export
Binomial <- function(name, n, p, succ = 1, fail = 0) .rvf(name, "BinomialDistribution", list(n, p, succ, fail))

.defdfin("BetaBinomialDistribution", c("n", "alpha", "beta"),
  check = function(n, alpha, beta) {
    .value_check(list(.st_fact(n, "integer"), .st_fact(n, "nonnegative")), paste0("'n' must be nonnegative integer. n = ", .rv_pystr(n), "."))
    .value_check(alpha > 0, paste0("'alpha' must be: alpha > 0 . alpha = ", .str(alpha)))
    .value_check(beta > 0, paste0("'beta' must be: beta > 0 . beta = ", .str(beta)))
  },
  high = function(d) .da(d, "n"),
  low = function(d) S_Zero,
  is_symbolic = function(d) !.rv_is_number(.da(d, "n")),
  set = function(d) {
    if (.fdist_symbolic(d)) return(Intersection(Naturals0, Interval(0, .da(d, "n"))))
    .pyset_order(lapply(0:as.integer(.da(d, "n")$value), as_sym))
  },
  pmf = function(d, k) {
    p <- .dl(d)
    sym_binomial(p$n, k) * sym_beta(k + p$alpha, p$n - k + p$beta) / sym_beta(p$alpha, p$beta)
  })

#' @rdname finite_rv
#' @export
BetaBinomial <- function(name, n, alpha, beta) .rvf(name, "BetaBinomialDistribution", list(n, alpha, beta))

.defdfin("HypergeometricDistribution", c("N", "m", "n"),
  check = function(n, N, m) {
    .value_check(list(.st_fact(N, "integer"), .st_fact(N, "nonnegative")), paste0("'N' must be nonnegative integer. N = ", .rv_pystr(N), "."))
    .value_check(list(.st_fact(n, "integer"), .st_fact(n, "nonnegative")), paste0("'n' must be nonnegative integer. n = ", .rv_pystr(n), "."))
    .value_check(list(.st_fact(m, "integer"), .st_fact(m, "nonnegative")), paste0("'m' must be nonnegative integer. m = ", .rv_pystr(m), "."))
  },
  is_symbolic = function(d) !all(vapply(d$args, .rv_is_number, NA)),
  high = function(d) { p <- .dl(d); Piecewise(list(p$n, as_sym(!identical(.fuzzy_bool(Lt(p$n, p$m)), FALSE))), list(p$m, S_true)) },
  low = function(d) { p <- .dl(d); Piecewise(list(0, as_sym(!identical(.fuzzy_bool(Gt(0, p$n + p$m - p$N)), FALSE))), list(p$n + p$m - p$N, S_true)) },
  set = function(d) {
    if (.fdist_symbolic(d)) return(Intersection(Naturals0, Interval(.callm(d, "low"), .callm(d, "high"))))
    p <- vapply(d$args, function(z) as.numeric(z$value), 0)
    lo <- max(0, p[3] + p[2] - p[1]); hi <- min(p[3], p[2])
    if (hi < lo) return(list())
    .pyset_order(lapply(lo:hi, as_sym))
  },
  pmf = function(d, k) {
    p <- .dl(d)
    sym_binomial(p$m, k) * sym_binomial(p$N - p$m, p$n - k) / sym_binomial(p$N, p$n)
  })

#' @rdname finite_rv
#' @export
Hypergeometric <- function(name, N, m, n) .rvf(name, "HypergeometricDistribution", list(N, m, n))

.defdfin("RademacherDistribution", character(),
  set = function(d) .pyset_order(list(as_sym(-1), as_sym(1))),
  pmf = function(d, ...) {
    k <- Dummy("k")
    f <- Lambda(k, Piecewise(list(S_Half, Or(Eq(k, -1), Eq(k, 1))), list(S_Zero, S_true)))
    if (...length()) apply_lambda(f, ...) else f
  })

#' @rdname finite_rv
#' @export
Rademacher <- function(name) .rvf(name, "RademacherDistribution", list())

.defdfin("IdealSolitonDistribution", c("k"),
  check = function(k) .value_check(.pyand(.st_fact(k, "integer"), .st_fact(k, "positive")), "'k' must be a positive integer."),
  low = function(d) S_One,
  high = function(d) .da(d, "k"),
  set = function(d) .pyset_order(lapply(seq_len(as.integer(.da(d, "k")$value)), as_sym)),
  dict = function(d) {
    k <- .da(d, "k")
    if (k$type == "Symbol") return(Density(d))
    kk <- as.integer(k$value)
    keys <- list(as_sym(1)); vals <- list(Rational(1, kk))
    if (kk >= 2) for (i in 2:kk) { keys[[length(keys) + 1L]] <- as_sym(i); vals[[length(vals) + 1L]] <- Rational(1, i * (i - 1)) }
    .FD(keys, vals)
  },
  pmf = function(d, x) {
    x <- .pmf_arg_check(x)
    k <- .da(d, "k")
    xi <- as_sym(ask_fact(x, "integer"))
    cond1 <- And(Eq(x, 1), xi)
    cond2 <- And(Ge(x, 1), Le(x, k), xi)
    Piecewise(list(1 / k, cond1), list(1 / (x * (x - 1)), cond2), list(S_Zero, S_true))
  })

#' @rdname finite_rv
#' @export
IdealSoliton <- function(name, k) .rvf(name, "IdealSolitonDistribution", list(k))

.rs_R <- function(d) { p <- .dl(d); p$c * log(p$k / p$delta) * p$k^Float(0.5) }
.py_round <- function(x) {
  v <- as.numeric(N(x))
  r <- round(v)
  if (abs(v - trunc(v)) == 0.5) r <- 2 * round(v / 2)
  as.integer(r)
}
.rs_Z <- function(d) {
  p <- .dl(d); R <- .rs_R(d)
  z <- S_Zero
  hi <- .py_round(p$k / R)
  if (hi > 1) for (i in 1:(hi - 1)) z <- z + 1 / as_sym(i)
  z <- z + log(R / p$delta)
  1 + z * R / p$k
}

.defdfin("RobustSolitonDistribution", c("k", "delta", "c"),
  check = function(k, delta, c) {
    .value_check(.pyand(.st_fact(k, "integer"), .st_fact(k, "positive")), "'k' must be a positive integer")
    .value_check(.pyand(Gt(delta, 0), Le(delta, 1)), "'delta' must be a real number in the interval (0,1)")
    .value_check(.st_fact(c, "positive"), "'c' must be a positive real number.")
  },
  R = function(d) .rs_R(d),
  Z = function(d) .rs_Z(d),
  low = function(d) S_One,
  high = function(d) .da(d, "k"),
  set = function(d) .pyset_order(lapply(seq_len(as.integer(.da(d, "k")$value)), as_sym)),
  is_symbolic = function(d) !(all(vapply(d$args, .rv_is_number, NA))),
  pmf = function(d, x) {
    x <- .pmf_arg_check(x)
    p <- .dl(d); R <- .rs_R(d)
    xi <- as_sym(ask_fact(x, "integer"))
    cond1 <- And(Eq(x, 1), xi)
    cond2 <- And(Ge(x, 1), Le(x, p$k), xi)
    rho <- Piecewise(list(1 / p$k, cond1), list(1 / (x * (x - 1)), cond2), list(S_Zero, S_true))
    rk <- .py_round(p$k / R)
    cond1 <- And(Ge(x, 1), Le(x, rk - 1))
    cond2 <- Eq(x, rk)
    tau <- Piecewise(list(R / (p$k * x), cond1), list(R * log(R / p$delta) / p$k, cond2), list(S_Zero, S_true))
    (rho + tau) / .rs_Z(d)
  })

#' @rdname finite_rv
#' @export
RobustSoliton <- function(name, k, delta, c) .rvf(name, "RobustSolitonDistribution", list(k, delta, c))
