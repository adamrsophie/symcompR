# Markov chains of sympy/stats/stochastic_process_types.py: the
# MarkovProcess mixin, DiscreteMarkovChain, ContinuousMarkovChain and the
# assumption objects TransitionMatrixOf, GeneratorMatrixOf and
# StochasticStateSpaceOf, with sample_stochastic_process(). A chain is a
# process node (see stats-stochastic.R) whose arguments are the symbol, the
# state space (a Tuple or a Range) and the transition or generator matrix
# (an ImmutableMatrix node or a MatrixSymbol).

.defclass("MarkovProcess", "StochasticProcess")
.defclass("DiscreteMarkovChain", c("DiscreteTimeStochasticProcess", "MarkovProcess"))
.defclass("ContinuousMarkovChain", c("ContinuousTimeStochasticProcess", "MarkovProcess"))
.defclass("TransitionMatrixOf")
.defclass("GeneratorMatrixOf", "TransitionMatrixOf")
.defclass("StochasticStateSpaceOf")

# ---------------------------------------------------------------------------
# Converters (_set_converter, _state_converter, _matrix_checks)

# Python's sympify of a state: strings are parsed, processes are their nodes
.stp_state_item <- function(i) {
  if (is.character(i)) return(sympify(i))
  if (inherits(i, "rv_process")) return(attr(i, "node"))
  as_sym(i)
}

# the hash key of a raw state, for the uniqueness check of a list
.stp_raw_key <- function(i) {
  if (is.character(i)) return(paste0("str:", i))
  .stp_state_item(i)$key
}

.stp_state_converter <- function(itr) {
  if (is_sym(itr) && itr$type %in% c("Tuple", "FiniteSet"))
    return(do.call(Tuple, lapply(itr$args, .stp_state_item)))
  if ((is.list(itr) && !is_sym(itr) && !inherits(itr, "symmatrix") && !is.function(itr)) ||
      (is.atomic(itr) && !is.null(itr) && !is.object(itr))) {
    itr <- as.list(itr)
    keys <- vapply(itr, .stp_raw_key, "")
    if (anyDuplicated(keys)) .serr("ValueError", "The state space must have unique elements.")
    return(do.call(Tuple, lapply(itr, .stp_state_item)))
  }
  if (is_sym(itr) && itr$type == "Range") {
    n <- .range_size(itr)
    if (n$type == "Integer") {
      k <- .int_val(n)
      return(do.call(Tuple, lapply(seq_len(k) - 1L, function(i) .range_item(itr, i))))
    }
    return(itr)
  }
  .serr("TypeError", paste0(.rv_pystr(itr), " is not an instance of list/tuple/set/Range/Tuple/FiniteSet."))
}

.stp_matrix_checks <- function(m) {
  if (inherits(m, "symmatrix")) m <- .as_block(m)
  if (!(is_sym(m) && m$type %in% c("ImmutableMatrix", "MatrixSymbol")))
    .serr("TypeError", "Transition probabilities either should be a Matrix or a MatrixSymbol.")
  sh <- .mshape(m)
  if (!identical(sh[[1L]]$key, sh[[2L]]$key)) .serr("NonSquareMatrixError", paste0(.str(m), " is not a square matrix"))
  m
}

# the size of a converted state space (Range objects give a symbolic size)
.stp_ss_size <- function(ss) {
  if (ss$type == "Range") {
    p <- .range_parts(ss)
    return(ceiling(.div(.sub(p$stop, p$start), p$step)))
  }
  .int(length(ss$args))
}

.stp_explicit <- function(m) if (is_sym(m) && m$type == "ImmutableMatrix") as_explicit(m) else m
.stp_is_msym <- function(m) is_sym(m) && m$type == "MatrixSymbol"
.stp_shape0 <- function(m) if (inherits(m, "symmatrix")) .int(m$nrow) else .mshape(m)[[1L]]

# MarkovProcess._sanity_checks
.stp_sanity_checks <- function(state_space, trans_probs) {
  if (is.null(state_space) && is.null(trans_probs)) {
    n <- Dummy("n", integer = TRUE, nonnegative = TRUE)
    state_space <- .stp_state_converter(Range(n))
    trans_probs <- .stp_matrix_checks(MatrixSymbol("_T", n, n))
  } else if (is.null(state_space)) {
    trans_probs <- .stp_matrix_checks(trans_probs)
    state_space <- .stp_state_converter(Range(.stp_shape0(trans_probs)))
  } else if (is.null(trans_probs)) {
    state_space <- .stp_state_converter(state_space)
    n <- .stp_ss_size(state_space)
    trans_probs <- MatrixSymbol("_T", n, n)
  } else {
    state_space <- .stp_state_converter(state_space)
    trans_probs <- .stp_matrix_checks(trans_probs)
    if (!identical(.stp_ss_size(state_space)$key, .stp_shape0(trans_probs)$key))
      .serr("ValueError", "The size of the state space and the number of rows of the transition matrix must be the same.")
  }
  list(state_space, trans_probs)
}

# ---------------------------------------------------------------------------
# MarkovProcess

.defm("MarkovProcess", "number_of_states", function(p) .stp_shape0(p$args[[3L]]))
.defm("MarkovProcess", "_state_index", function(p) p$args[[2L]])
.defm("MarkovProcess", "state_space", function(p) {
  ss <- p$args[[2L]]
  if (ss$type %in% c("FiniteSet", "Range")) ss else do.call(FiniteSet, ss$args)
})

# the elements of a finite set, Tuple or finite Range, in iteration order
.stp_elems <- function(s) {
  if (s$type %in% c("FiniteSet", "Tuple")) return(s$args)
  if (s$type == "Range") {
    n <- .range_size(s)
    if (n$type != "Integer") .serr("TypeError", paste0("Cannot iterate over symbolic Range"))
    return(lapply(seq_len(.int_val(n)) - 1L, function(i) .range_item(s, i)))
  }
  if (s$type == "EmptySet") return(list())
  if (s$type == "Union" && all(vapply(s$args, function(a) a$type == "FiniteSet", NA))) return(unlist(lapply(s$args, .stp_elems), recursive = FALSE))
  .serr("TypeError", paste0("'", .sympy_class(s), "' object is not iterable"))
}
.stp_len <- function(s) length(.stp_elems(s))

# index_of: the position of each state, for an integer number of states
.stp_index_of <- function(p) {
  if (p$type == "DiscreteMarkovChain") {
    if (.callm(p, "number_of_states")$type != "Integer") return(integer())
    st <- .stp_elems(p$args[[2L]])
  } else {
    if (.callm(p, "number_of_states")$type != "Integer") return(integer())
    st <- .stp_elems(.callm(p, "state_space"))
  }
  out <- integer()
  for (i in seq_along(st)) out[[st[[i]]$key]] <- i - 1L
  out
}
.stp_index_get <- function(idx, z) {
  if (!is_sym(z)) return(z)
  i <- match(z$key, names(idx))
  if (is.na(i)) z else .int(idx[[i]])
}

.stp_is_ris <- function(e) is_sym(e) && e$type == "RandomIndexedSymbol"
.stp_ris_atoms <- function(e) {
  if (is.null(e) || !is_sym(e)) return(list())
  .uniq(Filter(.stp_is_ris, .all_subexprs(e)))
}
.stp_is_symbolic_key <- function(k) is_sym(k) && k$type %in% c("Symbol", "Dummy", "Indexed")
# Python truth of a numeric comparison
.stp_true <- function(r) isTRUE(.fuzzy_bool(.as_bool_sym(r)))
.stp_bool <- function(x) if (isTRUE(x)) S_true else S_false

.stp_replace_with_index <- function(p, condition) {
  if (is_sym(condition) && is_relational(condition)) {
    lhs <- condition$args[[1L]]; rhs <- condition$args[[2L]]
    if (!.stp_is_ris(lhs)) { tmp <- lhs; lhs <- rhs; rhs <- tmp }
    idx <- .stp_index_of(p)
    condition <- .rel_same(condition, .stp_index_get(idx, lhs), .stp_index_get(idx, rhs))
  }
  condition
}

.stp_is_tmo <- function(e) is_sym(e) && e$type %in% c("TransitionMatrixOf", "GeneratorMatrixOf")

# MarkovProcess._extract_information
.stp_extract_information <- function(p, given_condition) {
  trans_probs <- p$args[[3L]]
  state_index <- p$args[[2L]]
  if (is_sym(given_condition) && given_condition$type == "And") {
    gcs <- .sympy_args(given_condition)
    given_condition <- S_true
    for (gc in gcs) {
      if (.stp_is_tmo(gc)) trans_probs <- gc$args[[2L]]
      if (is_sym(gc) && gc$type == "StochasticStateSpaceOf") state_index <- gc$args[[2L]]
      if (is_sym(gc) && is_relational(gc)) given_condition <- And(given_condition, gc)
    }
  }
  if (.stp_is_tmo(given_condition)) {
    trans_probs <- given_condition$args[[2L]]
    given_condition <- S_true
  }
  if (is_sym(given_condition) && given_condition$type == "StochasticStateSpaceOf") {
    state_index <- given_condition$args[[2L]]
    given_condition <- S_true
  }
  list(trans_probs, state_index, given_condition)
}

# Python's sum(): 0 plus each item, left to right
.stp_pysum <- function(items, start = S_Zero) {
  s <- start
  for (it in items) s <- .add(list(s, it))
  s
}

.stp_check_trans_probs <- function(trans_probs, row_sum = 1) {
  if (.stp_is_msym(trans_probs)) return(invisible(NULL))
  M <- .stp_explicit(trans_probs)
  for (row in .mrows(M)) {
    d <- .sub(.stp_pysum(row), as_sym(row_sum))
    if (!.stp_true(Eq(d, S_Zero)))
      .serr("ValueError", paste0("Values in a row must sum to ", row_sum, ". If you are using Float or floats then please use Rational."))
  }
  invisible(NULL)
}

.stp_work_out_state_index <- function(p, state_index, given_condition, trans_probs) {
  if (!is.null(given_condition)) {
    rand_var <- Filter(function(r) !.stp_is_ris(r), random_symbols(given_condition))
    if (length(rand_var) == 1L) state_index <- .callm(.rs_pspace(rand_var[[1L]]), "set")
  }
  sym_cond <- .callm(p, "number_of_states")$type != "Integer"
  n0 <- .stp_shape0(trans_probs)
  if (!sym_cond && !identical(.int(.stp_len(state_index))$key, n0$key))
    .serr("ValueError", "state space is not compatible with the transition probabilities.")
  if (!(n0$type %in% c("Symbol", "Dummy")))
    state_index <- do.call(FiniteSet, lapply(seq_len(.int_val(n0)) - 1L, .int))
  state_index
}

# MarkovProcess._preprocess: list(is_insufficient, matrix, state index, given)
.stp_preprocess <- function(p, given_condition, evaluate) {
  if (!evaluate) return(list(TRUE, NULL, NULL, NULL))
  info <- .stp_extract_information(p, given_condition)
  trans_probs <- info[[1L]]; state_index <- info[[2L]]; given_condition <- info[[3L]]
  is_insufficient <- FALSE
  if (is.null(trans_probs) || is.null(given_condition)) {
    is_insufficient <- TRUE
  } else {
    .stp_check_trans_probs(trans_probs, if (p$type == "DiscreteMarkovChain") 1 else 0)
    state_index <- .stp_work_out_state_index(p, state_index, given_condition, trans_probs)
  }
  list(is_insufficient, trans_probs, state_index, given_condition)
}

# Relational.as_set and Boolean.as_set in the random indexed symbol
.stp_as_set <- function(cond) {
  rvs <- .stp_ris_atoms(cond)
  if (length(rvs) != 1L) return(.sset_as_set(cond))
  x <- Dummy("x", real = TRUE)
  c2 <- xreplace(cond, list(list(rvs[[1L]], x)))
  .sset_as_set(c2)
}

# the matrix power of an explicit or symbolic transition matrix
.stp_mpow <- function(m, k) {
  if (.stp_is_msym(m)) return(MatPow(m, k))
  M <- .stp_explicit(m)
  k <- as_sym(k)
  if (k$type == "Integer" && .int_val(k) >= 1) return(.stp_pow_rec(M, .int_val(k)))
  if (k$type != "Integer") {
    # _matrix_pow_by_jordan_blocks, unsimplified, for a diagonal Jordan form
    jf <- tryCatch(jordan_form(M), error = function(e) NULL)
    if (!is.null(jf)) {
      J <- jf$J; n <- J$nrow
      if (all(vapply(seq_len(n * n), function(q) (q - 1L) %/% n == (q - 1L) %% n || is_zero_num(J$data[[q]]), NA))) {
        D <- diag_matrix(lapply(seq_len(n), function(i) .pow(.mget(J, i, i), k)))
        return(.matmul(.matmul(jf$P, D), inverse(jf$P)))
      }
    }
  }
  M %^% k
}
# MatrixBase._eval_pow_by_recursion
.stp_pow_rec <- function(M, num) {
  if (num == 1) return(M)
  if (num %% 2 == 1) { a <- M; b <- .stp_pow_rec(M, num - 1) }
  else { a <- b <- .stp_pow_rec(M, num %/% 2) }
  .matmul(a, b)
}
.stp_getitem <- function(M, i, j) {
  if (inherits(M, "symmatrix")) return(MatrixElement(.as_block(M), i, j))
  .me_getitem(M, as_sym(i), as_sym(j))
}

.stp_mat_of <- function(p, mat) {
  node_cls <- if (p$type == "DiscreteMarkovChain") "TransitionMatrixOf" else "GeneratorMatrixOf"
  new_sym(node_cls, list(p, .stp_matrix_checks(mat)))
}
.stp_ss_of <- function(p, state_space) {
  ss <- .stp_state_converter(state_space)
  new_sym("StochasticStateSpaceOf", list(p, Range(.stp_ss_size(ss))))
}

.stp_probability <- function(p, condition, given_condition = NULL, evaluate = TRUE, ...) {
  pre <- .stp_preprocess(p, given_condition, evaluate)
  check <- pre[[1L]]; mat <- pre[[2L]]; state_index <- pre[[3L]]; new_given_condition <- pre[[4L]]
  rv <- .stp_ris_atoms(condition)
  symbolic <- any(vapply(rv, function(s) .stp_is_symbolic_key(.ris_key(s)), NA))
  if (check) return(Probability(condition, new_given_condition))
  trans_probs <- if (p$type == "ContinuousMarkovChain") .stp_cmc_transition(p, mat) else mat
  condition <- .stp_replace_with_index(p, condition)
  given_condition <- .stp_replace_with_index(p, given_condition)
  new_given_condition <- .stp_replace_with_index(p, new_given_condition)
  if (is_relational(condition)) {
    gcs <- if (new_given_condition$type == "And") .sympy_args(new_given_condition) else list(new_given_condition)
    min_key_rv <- .stp_ris_atoms(new_given_condition)
    if (length(min_key_rv)) {
      min_key_rv <- min_key_rv[[1L]]
      for (r in rv) {
        if (.stp_is_symbolic_key(.ris_key(min_key_rv)) || .stp_is_symbolic_key(.ris_key(r))) next
        if (.stp_true(Gt(.ris_key(min_key_rv), .ris_key(r)))) return(Probability(condition))
      }
    } else {
      return(Probability(condition))
    }
    if (symbolic) return(.stp_symbolic_probability(p, condition, new_given_condition, rv, min_key_rv))
    if (length(rv) > 1L) {
      rv0 <- condition$args[[1L]]; rv1 <- condition$args[[2L]]
      if (.stp_true(Lt(.ris_key(rv0), .ris_key(rv1)))) {
        tmp <- rv0; rv0 <- rv1; rv1 <- tmp
        lhs <- condition$args[[1L]]; rhs <- condition$args[[2L]]
        condition <- switch(condition$type, StrictGreaterThan = Lt(lhs, rhs), StrictLessThan = Gt(lhs, rhs),
                            GreaterThan = Le(lhs, rhs), LessThan = Ge(lhs, rhs), condition)
      }
      s <- S_Zero
      n <- .stp_len(.callm(p, "state_space"))
      prob <- function(c, g) .stp_probability(p, c, g)
      if (condition$type %in% c("Equality", "Unequality")) {
        for (i in seq_len(n) - 1L)
          s <- .add(list(s, .mul(list(prob(Eq(rv0, i), Eq(rv1, i)), prob(Eq(rv1, i), new_given_condition)))))
        return(if (condition$type == "Equality") s else .sub(S_One, s))
      }
      upper <- if (condition$type %in% c("GreaterThan", "StrictLessThan")) 1L else 0L
      greater <- condition$type %in% c("GreaterThan", "StrictGreaterThan")
      for (i in seq_len(n) - 1L) {
        if (i <= n %/% 2L) {
          for (j in seq_len(i + upper) - 1L)
            s <- .add(list(s, .mul(list(prob(Eq(rv0, i), Eq(rv1, j)), prob(Eq(rv1, j), new_given_condition)))))
        } else {
          s <- .add(list(s, prob(Eq(rv0, i), new_given_condition)))
          if (i + upper <= n - 1L) for (j in (i + upper):(n - 1L))
            s <- .sub(s, .mul(list(prob(Eq(rv0, i), Eq(rv1, j)), prob(Eq(rv1, j), new_given_condition))))
        }
      }
      return(if (greater) s else .sub(S_One, s))
    }
    rv <- rv[[1L]]
    states <- .stp_as_set(condition)
    prob <- list(); gstate <- NULL
    for (gc in gcs) {
      if (.has1(gc, min_key_rv)) {
        if (.has_type(gc, "Probability")) {
          if (gc$args[[1L]]$type == "Probability") { pv <- gc$args[[2L]]; gp <- gc$args[[1L]] }
          else { pv <- gc$args[[1L]]; gp <- gc$args[[2L]] }
          gr <- gp$args[[1L]]
          gset <- Intersection(.stp_as_set(gr), state_index)
          gstate <- .stp_elems(gset)[[1L]]
          prob[[gset$key]] <- list(gset, pv)
        } else {
          gstate <- if (.stp_is_ris(gc$args[[1L]])) gc$args[[2L]] else gc$args[[1L]]
        }
      }
    }
    iset <- .callm(p, "index_set")
    for (k in list(.ris_key(rv), .ris_key(min_key_rv))) {
      inside <- .fuzzy_bool(Contains(k, iset))
      if (identical(inside, FALSE)) .serr("IndexError", "The timestamps of the process are not in it's index set.")
    }
    nos <- .callm(p, "number_of_states")
    if (!(nos$type %in% c("Symbol", "Dummy"))) states <- Intersection(states, state_index)
    n0 <- .stp_shape0(mat)
    for (state in .stp_elems(Union(states, FiniteSet(gstate)))) {
      if (state$type != "Integer" || .stp_true(Ge(state, n0)))
        .serr("IndexError", paste0("No information is available for (", .str(gstate), ", ", .str(state),
                                   ") in transition probabilities of shape, (", .str(n0), ", ", .str(.mshape(mat)[[2L]]),
                                   "). State space is zero indexed."))
    }
    if (length(prob)) {
      gstates <- do.call(Union, unname(lapply(prob, `[[`, 1L)))
      ng <- .stp_len(gstates)
      if (ng == 1L) {
        gstate <- .stp_elems(gstates)[[1L]]
        gprob <- prob[[1L]][[2L]]
        prob[[gstates$key]] <- list(gstates, gprob)
      } else if (ng == .stp_len(state_index) - 1L) {
        rest <- Complement(state_index, gstates)
        gstate <- .stp_elems(rest)[[1L]]
        gprob <- .sub(S_One, .stp_pysum(lapply(prob, `[[`, 2L)))
        prob[[rest$key]] <- list(rest, gprob)
      } else {
        .serr("ValueError", "Conflicting information.")
      }
    } else {
      gprob <- S_One
    }
    sts <- .stp_elems(states)
    if (identical(min_key_rv$key, rv$key)) {
      return(.stp_pysum(lapply(sts, function(state) {
        fs <- FiniteSet(state)
        hit <- prob[[fs$key]]
        if (is.null(hit)) .serr("KeyError", .str(fs))
        hit[[2L]]
      })))
    }
    dk <- .sub(.ris_key(rv), .ris_key(min_key_rv))
    M <- if (p$type == "ContinuousMarkovChain") .stp_explicit(apply_lambda(trans_probs, dk)) else .stp_mpow(trans_probs, dk)
    return(.mul(list(gprob, .stp_pysum(lapply(sts, function(state) .stp_getitem(M, gstate, state))))))
  }
  if (condition$type == "Not")
    return(.sub(S_One, .stp_probability(p, condition$args[[1L]], given_condition, evaluate, ...)))
  if (condition$type == "And") {
    compute_later <- list(); state2cond <- list(); ris_of <- list()
    for (expr in .sympy_args(condition)) {
      if (is_relational(expr)) {
        ris <- .stp_ris_atoms(expr)[[1L]]
        if (is.null(state2cond[[ris$key]])) { state2cond[[ris$key]] <- S_true; ris_of[[ris$key]] <- ris }
        state2cond[[ris$key]] <- And(state2cond[[ris$key]], expr)
      } else {
        compute_later[[length(compute_later) + 1L]] <- expr
      }
    }
    for (k in names(state2cond)) {
      cset <- Intersection(.stp_as_set(state2cond[[k]]), state_index)
      if (.stp_len(cset) == 0L) return(S_Zero)
      state2cond[[k]] <- as_relational(cset, ris_of[[k]])
    }
    keys <- names(state2cond)
    ord <- .stp_sort_by_key(lapply(keys, function(k) .ris_key(ris_of[[k]])))
    sorted_ris <- keys[ord]
    prod <- .stp_probability(p, state2cond[[sorted_ris[1L]]], given_condition, evaluate, ...)
    for (i in seq_along(sorted_ris)[-1L]) {
      ri <- sorted_ris[i]; prev_ri <- sorted_ris[i - 1L]
      if (state2cond[[ri]]$type != "Equality")
        .serr("ValueError", paste0("The process is in multiple states at ", .str(ris_of[[ri]]), ", unable to determine the probability."))
      mat_of <- .stp_mat_of(p, mat)
      prod <- .mul(list(prod, .stp_probability(p, state2cond[[ri]],
                                               And(state2cond[[prev_ri]], mat_of, .stp_ss_of(p, state_index)),
                                               evaluate, ...)))
    }
    for (expr in compute_later) prod <- .mul(list(prod, .stp_probability(p, expr, given_condition, evaluate, ...)))
    return(prod)
  }
  if (condition$type == "Or")
    return(.stp_pysum(lapply(.sympy_args(condition), function(expr) .stp_probability(p, expr, given_condition, evaluate, ...))))
  .serr("NotImplementedError", paste0("Mechanism for handling (", .str(condition), ", ", .stp_none(given_condition),
                                      ") queries hasn't been implemented yet."))
}

.stp_none <- function(x) if (is.null(x)) "None" else .str(x)

# sorted() with numeric keys, stable
.stp_sort_by_key <- function(keys) {
  n <- length(keys)
  ord <- seq_len(n)
  if (n < 2L) return(ord)
  for (i in 2:n) {
    j <- i
    while (j > 1L && .stp_true(Lt(keys[[ord[j]]], keys[[ord[j - 1L]]]))) {
      tmp <- ord[j]; ord[j] <- ord[j - 1L]; ord[j - 1L] <- tmp
      j <- j - 1L
    }
  }
  ord
}

.stp_symbolic_probability <- function(p, condition, new_given_condition, rv, min_key_rv) {
  if (!is_relational(condition)) return(Probability(condition, new_given_condition))
  curr_state <- if (.stp_is_ris(new_given_condition$args[[1L]])) new_given_condition$args[[2L]] else new_given_condition$args[[1L]]
  next_state <- if (.stp_is_ris(condition$args[[1L]])) condition$args[[2L]] else condition$args[[1L]]
  dk <- .sub(.ris_key(rv[[1L]]), .ris_key(min_key_rv))
  if (condition$type %in% c("Equality", "Unequality")) {
    if (p$type == "DiscreteMarkovChain") {
      P <- .stp_mpow(p$args[[3L]], dk)
    } else {
      P <- matrix_exp(.mat_map(.stp_explicit(p$args[[3L]]), function(z) .mul(list(z, dk))))
    }
    el <- .stp_getitem(P, curr_state, next_state)
    prob <- if (condition$type == "Equality") el else .sub(S_One, el)
    return(Piecewise(list(prob, Gt(.ris_key(rv[[1L]]), .ris_key(min_key_rv))), list(Probability(condition), S_true)))
  }
  upper <- 1L
  greater <- FALSE
  if (condition$type %in% c("GreaterThan", "StrictLessThan")) upper <- 0L
  if (condition$type %in% c("GreaterThan", "StrictGreaterThan")) greater <- TRUE
  k <- Dummy("k")
  condition <- if (.stp_is_ris(condition$args[[1L]])) Eq(condition$args[[1L]], k) else Eq(condition$args[[2L]], k)
  total <- Sum(.stp_probability(p, condition, new_given_condition),
               list(k, .add(list(next_state, .int(upper))), set_sup(.callm(p, "state_space"))))
  kc <- Gt(.ris_key(rv[[1L]]), .ris_key(min_key_rv))
  if (greater) Piecewise(list(total, kc), list(Probability(condition), S_true))
  else Piecewise(list(.sub(S_One, total), kc), list(Probability(condition), S_true))
}

.stp_expectation <- function(p, expr, condition = NULL, evaluate = TRUE, ...) {
  pre <- .stp_preprocess(p, condition, evaluate)
  check <- pre[[1L]]; mat <- pre[[2L]]; state_index <- pre[[3L]]; condition <- pre[[4L]]
  if (check) return(Expectation(expr, condition))
  rvs <- random_symbols(expr)
  if (is_sym(expr) && !(is_relational(expr) || expr$type %in% c("And", "Or", "Not", "BooleanAtom")) && is_sym(condition) && condition$type == "Equality" && length(rvs) == 1L) {
    condition <- .stp_replace_with_index(p, condition)
    rv <- rvs[[1L]]
    lhsg <- condition$args[[1L]]; rhsg <- condition$args[[2L]]
    if (!.stp_is_ris(lhsg)) { tmp <- lhsg; lhsg <- rhsg; rhsg <- tmp }
    if (!.stp_true(Contains(rhsg, state_index)))
      .serr("ValueError", paste0(.str(rhsg), " state is not in the state space."))
    if (.stp_true(Lt(.ris_key(rv), .ris_key(lhsg))))
      .serr("ValueError", paste0("Incorrect given condition is given, expectation time ", .str(.ris_key(rv)),
                                 " < time ", .str(.ris_key(rv))))
    cond <- And(condition, .stp_mat_of(p, mat), .stp_ss_of(p, state_index))
    sidx <- p$args[[2L]]
    return(.stp_pysum(lapply(.stp_elems(state_index), function(s) {
      st <- if (sidx$type == "Range") .range_item(sidx, .int_val(s)) else sidx$args[[.int_val(s) + 1L]]
      .mul(list(.stp_probability(p, Eq(rv, s), cond), subs(expr, rv, st)))
    })))
  }
  .serr("NotImplementedError", paste0("Mechanism for handling (", .str(expr), ", ", .stp_none(condition),
                                      ") queries hasn't been implemented yet."))
}

.defm("MarkovProcess", "probability", function(p, condition, given_condition = NULL, evaluate = TRUE, ...)
  .stp_probability(p, condition, given_condition, evaluate, ...))
.defm("MarkovProcess", "expectation", function(p, expr, condition = NULL, evaluate = TRUE, ...)
  .stp_expectation(p, expr, condition, evaluate, ...))
.defm("MarkovProcess", "replace_with_index", function(p, condition) .stp_replace_with_index(p, condition))
.defm("MarkovProcess", "index_of", function(p) {
  idx <- .stp_index_of(p)
  if (!length(idx)) return(Dict(list(), list()))
  st <- lapply(names(idx), function(k) NULL)
  ss <- if (p$type == "DiscreteMarkovChain") .stp_elems(p$args[[2L]]) else .stp_elems(.callm(p, "state_space"))
  keys <- list(); vals <- list()
  for (s in ss) {
    if (!.in_list(s, keys)) keys[[length(keys) + 1L]] <- s
  }
  Dict(keys, lapply(keys, function(s) .int(idx[[s$key]])))
})

# ---------------------------------------------------------------------------
# DiscreteMarkovChain

.defm("DiscreteMarkovChain", "index_set", function(p) Naturals0)
.defm("DiscreteMarkovChain", "transition_probabilities", function(p) .stp_explicit(p$args[[3L]]))

.stp_require_explicit <- function(T) {
  if (.stp_is_msym(T)) .serr("NotImplementedError", "Cannot perform the operation with a symbolic matrix.")
  .stp_explicit(T)
}

.stp_communication_classes <- function(p) {
  n <- .callm(p, "number_of_states")
  T <- .stp_require_explicit(p$args[[3L]])
  nn <- .int_val(n)
  V <- seq_len(nn) - 1L
  nz <- function(i, j) !is_zero_num(.mget(T, i + 1L, j + 1L))
  E <- list()
  for (i in V) for (j in V) if (nz(i, j)) E[[length(E) + 1L]] <- list(.int(i), .int(j))
  classes <- strongly_connected_components(list(lapply(V, .int), E))
  classes <- lapply(classes, function(cl) vapply(cl, function(z) as.integer(.int_val(as_sym(z))), 0L))
  recurrence <- list(); periods <- integer()
  for (class_ in classes) {
    is_recurrent <- S_true
    for (i in class_) {
      row <- lapply(class_, function(j) .mget(T, i + 1L, j + 1L))
      if (!.stp_true(Eq(.sub(.stp_pysum(row), S_One), S_Zero))) { is_recurrent <- S_false; break }
    }
    recurrence[[length(recurrence) + 1L]] <- is_recurrent
    # breadth first search for the period (Python sets of small integers
    # iterate in increasing order)
    non_tree <- integer()
    visited <- class_[1L]; newly_visited <- class_[1L]
    level <- integer(); level[as.character(class_[1L])] <- 0L
    current_level <- 0L
    done <- FALSE
    while (!done) {
      done <- length(visited) == length(class_)
      current_level <- current_level + 1L
      for (i in sort(newly_visited)) {
        newly_visited <- class_[vapply(class_, function(j) nz(i, j), NA)]
        new_tree <- setdiff(newly_visited, visited)
        for (j in new_tree) level[as.character(j)] <- current_level
        new_non_tree <- intersect(newly_visited, visited)
        vals <- vapply(new_non_tree, function(j) level[[as.character(i)]] - level[[as.character(j)]] + 1L, 0L)
        non_tree <- union(non_tree, vals)
        visited <- union(visited, new_tree)
      }
    }
    pos <- non_tree[non_tree > 0L]
    periods <- c(periods, if (!length(pos)) length(class_) else if (length(pos) == 1L) pos else Reduce(function(a, b) { while (b) { t <- a %% b; a <- b; b <- t }; a }, pos))
  }
  sidx <- p$args[[2L]]
  state_of <- function(i) if (sidx$type == "Range") .range_item(sidx, i) else sidx$args[[i + 1L]]
  lapply(seq_along(classes), function(k)
    .pytuple(list(lapply(classes[[k]], state_of), recurrence[[k]], .int(periods[k]))))
}

.stp_decompose <- function(p) {
  trans_probs <- .stp_require_explicit(p$args[[3L]])
  classes <- .stp_communication_classes(p)
  r_states <- list(); t_states <- list()
  for (cl in classes) {
    if (.is_true(cl[[2L]])) r_states <- c(r_states, cl[[1L]]) else t_states <- c(t_states, cl[[1L]])
  }
  states <- c(r_states, t_states)
  idx <- .stp_index_of(p)
  indexes <- vapply(states, function(s) idx[[s$key]], 0L) + 1L
  nr <- length(r_states); nt <- length(t_states)
  sub <- function(rows, cols) .new_mat(unlist(lapply(rows, function(i) lapply(cols, function(j) .mget(trans_probs, i, j))), recursive = FALSE),
                                      length(rows), length(cols))
  A <- sub(indexes[seq_len(nr)], indexes[seq_len(nr)])
  B <- sub(indexes[nr + seq_len(nt)], indexes[seq_len(nr)])
  C <- sub(indexes[nr + seq_len(nt)], indexes[nr + seq_len(nt)])
  .pytuple(list(states, A, B, C))
}

.stp_canonical_form <- function(p) {
  d <- .stp_decompose(p)
  A <- d[[2L]]; B <- d[[3L]]; C <- d[[4L]]
  O <- zeros(A$nrow, C$ncol)
  n <- A$nrow + C$nrow
  rows <- c(lapply(seq_len(A$nrow), function(i) c(.mrows(A)[[i]], .mrows(O)[[i]])),
            lapply(seq_len(B$nrow), function(i) c(.mrows(B)[[i]], .mrows(C)[[i]])))
  M <- .new_mat(unlist(rows, recursive = FALSE), n, n)
  .pytuple(list(d[[1L]], M))
}

.stp_stationary_distribution <- function(p, condition_set = FALSE) {
  trans_probs <- p$args[[3L]]
  n <- .callm(p, "number_of_states")
  if (is_zero_num(n)) return(.new_mat(list(), 1L, 0L))
  if (.stp_is_msym(trans_probs)) {
    wm <- MatrixSymbol("wm", 1, n)
    if (condition_set) return(ConditionSet(wm, Eq(MatMul(wm, trans_probs), wm)))
    return(Lambda(list(wm, trans_probs), Eq(MatMul(wm, trans_probs), wm)))
  }
  nn <- .int_val(n)
  P <- .stp_explicit(trans_probs)
  a <- t(P - eye(nn))
  for (j in seq_len(nn)) a$data[[j]] <- S_One
  b <- zeros(nn, 1L)
  b$data[[1L]] <- S_One
  soln <- .stp_elems(linsolve(list(a, b)))[[1L]]
  .new_mat(soln$args, 1L, nn)
}

.stp_fundamental_matrix <- function(p) {
  Q <- .stp_decompose(p)[[4L]]
  if (Q$nrow > 0L) {
    I <- eye(Q$nrow)
    if (is_zero_num(det(I - Q))) .serr("ValueError", "The fundamental matrix doesn't exist.")
    return(inverse(I - Q))
  }
  P <- .stp_explicit(p$args[[3L]])
  I <- eye(P$nrow)
  w <- .stp_stationary_distribution(p)
  W <- .new_mat(rep(w$data, P$nrow), P$nrow, P$nrow)
  if (is_zero_num(det(I - P + W))) .serr("ValueError", "The fundamental matrix doesn't exist.")
  inverse(I - P + W)
}

.stp_absorbing_probabilities <- function(p) {
  R <- .stp_decompose(p)[[3L]]
  N <- .stp_fundamental_matrix(p)
  if (is.null(R) || is.null(N)) return(NULL)
  .matmul(N, R)
}

.defm("DiscreteMarkovChain", "communication_classes", function(p) .stp_communication_classes(p))
.defm("DiscreteMarkovChain", "decompose", function(p) .stp_decompose(p))
.defm("DiscreteMarkovChain", "canonical_form", function(p) .stp_canonical_form(p))
.defm("DiscreteMarkovChain", "fundamental_matrix", function(p) .stp_fundamental_matrix(p))
.defm("DiscreteMarkovChain", "absorbing_probabilities", function(p) .stp_absorbing_probabilities(p))
.defm("DiscreteMarkovChain", "absorbing_probabilites", function(p) .stp_absorbing_probabilities(p))
.defm("DiscreteMarkovChain", "is_regular", function(p) {
  tuples <- .stp_communication_classes(p)
  if (!length(tuples)) return(S_false)
  And(.stp_bool(length(tuples) == 1L), .stp_bool(.int_val(tuples[[1L]][[3L]]) == 1))
})
.defm("DiscreteMarkovChain", "is_ergodic", function(p) {
  tuples <- .stp_communication_classes(p)
  if (!length(tuples)) return(S_false)
  .stp_bool(length(tuples) == 1L)
})
.defm("DiscreteMarkovChain", "is_absorbing_state", function(p, state) {
  tp <- p$args[[3L]]
  if (tp$type == "ImmutableMatrix" && .stp_true(Lt(as_sym(state), .stp_shape0(tp)))) {
    s <- .int_val(as_sym(state))
    return(identical(.mget(as_explicit(tp), s + 1L, s + 1L)$key, S_One$key))
  }
  NULL
})
.defm("DiscreteMarkovChain", "is_absorbing_chain", function(p) {
  A <- .stp_decompose(p)[[2L]]
  r <- A$nrow
  And(.stp_bool(r > 0L), .stp_bool(r > 0L && isTRUE(A == eye(r))))
})
.defm("DiscreteMarkovChain", "stationary_distribution", function(p, condition_set = FALSE) .stp_stationary_distribution(p, condition_set))
.defm("DiscreteMarkovChain", "fixed_row_vector", function(p) .stp_stationary_distribution(p))
# a property in SymPy: further arguments call the result (a Lambda)
.defm("DiscreteMarkovChain", "limiting_distribution", function(p, ...) {
  r <- .stp_stationary_distribution(p)
  if (!...length()) return(r)
  # Lambda.__call__ substitutes the matrices for the matrix symbols
  vals <- lapply(list(...), function(v) .as_mexpr(as_matexpr(v)))
  xreplace(r$args[[2L]], Map(list, r$args[[1L]]$args, vals))
})

# the first state is drawn uniformly; SymPy draws the later ones with
# sample_iter(FiniteRV(...)), which needs scipy and fails without it
.defm("DiscreteMarkovChain", "sample", function(p) {
  started <- FALSE
  function() {
    tp <- p$args[[3L]]
    if (tp$type != "ImmutableMatrix") .serr("ValueError", "Transition Matrix must be provided for sampling")
    if (started) .serr("ValueError", "Failed to import scipy")
    states <- .stp_elems(.callm(p, "state_space"))
    started <<- TRUE
    states[[sample.int(length(states), 1L)]]
  }
})

#' Markov chains
#'
#' `DiscreteMarkovChain()` and `ContinuousMarkovChain()` of `sympy.stats`:
#' finite, time homogeneous Markov chains given by a state space and a
#' transition matrix (or a generator matrix for continuous time). A
#' discrete chain is indexed with [rv_getitem()] (SymPy's `Y[t]`), a
#' continuous chain is called, `C(t)`. Queries go through [P()] and [E()];
#' pass a transition matrix, a generator matrix or a state space in the given
#' condition with `TransitionMatrixOf()`, `GeneratorMatrixOf()` and
#' `StochasticStateSpaceOf()`, combined with `&`. Methods such as
#' `"communication_classes"`, `"decompose"`, `"canonical_form"`,
#' `"fundamental_matrix"`, `"absorbing_probabilities"`, `"is_regular"`,
#' `"is_ergodic"`, `"is_absorbing_chain"`, `"stationary_distribution"`,
#' `"fixed_row_vector"`, `"limiting_distribution"` and
#' `"transition_probabilities"` are called with [rv_attr()]. Python tuples come
#' back as lists.
#'
#' @param sym The name of the chain.
#' @param state_space Optional state space: a list of states, a `Tuple`, a
#'   `FiniteSet` or a `Range`; by default `Range(n)`.
#' @param trans_probs Optional transition matrix (a `Matrix` or a
#'   `MatrixSymbol`); by default `MatrixSymbol("_T", n, n)`.
#' @param gen_mat Optional generator matrix of a continuous chain.
#' @param process A Markov chain.
#' @param matrix A transition or generator matrix.
#' @return A process (an R function), or an assumption to use in a given
#'   condition.
#' @examples
#' T <- Matrix(list(list(Rational(1, 2), Rational(1, 2)), list(Rational(1, 4), Rational(3, 4))))
#' Y <- DiscreteMarkovChain("Y", list(0, 1), T)
#' P(Eq(rv_getitem(Y, 3), 1), Eq(rv_getitem(Y, 1), 0))
#' rv_attr(Y, "stationary_distribution")
#' rv_attr(Y, "communication_classes")
#' G <- Matrix(list(list(-1, 1), list(1, -1)))
#' C <- ContinuousMarkovChain("C", list(0, 1), G)
#' rv_attr(C, "limiting_distribution")
#' @name markov_chains
NULL

#' @rdname markov_chains
#' @export
DiscreteMarkovChain <- function(sym, state_space = NULL, trans_probs = NULL) {
  sym <- .symbol_converter(sym)
  sc <- .stp_sanity_checks(state_space, trans_probs)
  .mkproc(new_sym("DiscreteMarkovChain", list(sym, sc[[1L]], sc[[2L]])))
}

# ---------------------------------------------------------------------------
# ContinuousMarkovChain

.defm("ContinuousMarkovChain", "index_set", function(p) Reals)
.defm("ContinuousMarkovChain", "generator_matrix", function(p) .stp_explicit(p$args[[3L]]))

# transition_probabilities(gen_mat): a Lambda of the time
.stp_cmc_transition <- function(p, gen_mat = NULL) {
  t <- Dummy("t")
  if (!is.null(gen_mat) && (inherits(gen_mat, "symmatrix") || (is_sym(gen_mat) && gen_mat$type == "ImmutableMatrix"))) {
    G <- .stp_explicit(gen_mat)
    if (is_diagonalizable(G)) {
      d <- diagonalize(G)
      Q <- d$P; D <- d$D
      eD <- diag_matrix(lapply(seq_len(D$nrow), function(i) exp(.mul(list(t, .mget(D, i, i))))))
      return(Lambda(t, .as_block(.matmul(.matmul(Q, eD), inverse(Q)))))
    }
  }
  if (!is.null(gen_mat)) {
    g <- if (inherits(gen_mat, "symmatrix")) .as_block(gen_mat) else gen_mat
    return(Lambda(t, exp(MatMul(t, g))))
  }
  NULL
}
.defm("ContinuousMarkovChain", "transition_probabilities", function(p, gen_mat = NULL) .stp_cmc_transition(p, gen_mat))

.defm("ContinuousMarkovChain", "limiting_distribution", function(p) {
  gen_mat <- p$args[[3L]]
  if (is.null(gen_mat)) return(NULL)
  if (.stp_is_msym(gen_mat)) {
    wm <- MatrixSymbol("wm", 1, .stp_shape0(gen_mat))
    return(Lambda(list(wm, gen_mat), Eq(MatMul(wm, gen_mat), wm)))
  }
  G <- as_explicit(gen_mat)
  w <- IndexedBase("w")
  wi <- lapply(seq_len(G$nrow) - 1L, function(i) Indexed("w", i))
  wm <- .new_mat(wi, 1L, length(wi))
  eqs <- .matmul(wm, G)$data
  eqs[[length(eqs) + 1L]] <- .sub(.stp_pysum(wi), S_One)
  soln <- .stp_elems(linsolve(eqs, wi))[[1L]]
  .new_mat(soln$args, 1L, length(wi))
})

#' @rdname markov_chains
#' @export
ContinuousMarkovChain <- function(sym, state_space = NULL, gen_mat = NULL) {
  sym <- .symbol_converter(sym)
  sc <- .stp_sanity_checks(state_space, gen_mat)
  .mkproc(new_sym("ContinuousMarkovChain", list(sym, sc[[1L]], sc[[2L]])))
}

# ---------------------------------------------------------------------------
# Assumptions for the given condition

#' @rdname markov_chains
#' @export
TransitionMatrixOf <- function(process, matrix) {
  process <- .pnode(process)
  if (!(is_sym(process) && process$type == "DiscreteMarkovChain"))
    .serr("ValueError", "Currently only DiscreteMarkovChain support TransitionMatrixOf.")
  new_sym("TransitionMatrixOf", list(process, .stp_matrix_checks(matrix)))
}

#' @rdname markov_chains
#' @export
GeneratorMatrixOf <- function(process, matrix) {
  process <- .pnode(process)
  if (!(is_sym(process) && process$type == "ContinuousMarkovChain"))
    .serr("ValueError", "Currently only ContinuousMarkovChain support GeneratorMatrixOf.")
  new_sym("GeneratorMatrixOf", list(process, .stp_matrix_checks(matrix)))
}

#' @rdname markov_chains
#' @export
StochasticStateSpaceOf <- function(process, state_space) {
  process <- .pnode(process)
  if (!(is_sym(process) && process$type %in% c("DiscreteMarkovChain", "ContinuousMarkovChain")))
    .serr("ValueError", "Currently only DiscreteMarkovChain and ContinuousMarkovChain support StochasticStateSpaceOf.")
  .stp_ss_of(process, state_space)
}

.defm("TransitionMatrixOf", "process", function(e) e$args[[1L]])
.defm("TransitionMatrixOf", "matrix", function(e) .stp_explicit(e$args[[2L]]))
.defm("StochasticStateSpaceOf", "process", function(e) e$args[[1L]])
.defm("StochasticStateSpaceOf", "state_index", function(e) e$args[[2L]])

#' Sample a stochastic process
#'
#' `sample_stochastic_process(process)` of `sympy.stats`: a generator of the
#' states of a discrete Markov chain. It is an R function; each call returns
#' the next state. The first state is drawn uniformly from the state space
#' with R's random number generator (use `set.seed()`). SymPy draws the later
#' states with scipy, so a second call raises SymPy's error when scipy is
#' missing.
#'
#' @param process A stochastic process.
#' @return A function returning the next sampled state at each call.
#' @examples
#' set.seed(1)
#' T <- Matrix(list(list(Rational(1, 2), Rational(1, 2)), list(Rational(1, 4), Rational(3, 4))))
#' g <- sample_stochastic_process(DiscreteMarkovChain("Y", list(0, 1), T))
#' g()
#' @export
sample_stochastic_process <- function(process) {
  if (!inherits(process, "rv_process"))
    .serr("ValueError", "Process must be an instance of Stochastic Process")
  .callm(attr(process, "node"), "sample")
}
