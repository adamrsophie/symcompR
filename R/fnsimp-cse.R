# Common subexpression elimination after sympy/simplify/cse_main.py: opt_cse
# (negations, negative powers and the common arguments of sums and products
# through match_common_args), tree_cse and cse.
#
# An Unevaluated node of SymPy is a list of class "fs_uneval" with func,
# args and an id; SymPy hashes it by identity, so tables key it by that id.

.fs_cse_env <- new.env(parent = emptyenv())
.fs_cse_env$uid <- 0L

.fs_uneval <- function(func, args) {
  .fs_cse_env$uid <- .fs_cse_env$uid + 1L
  structure(list(func = func, args = args, id = .fs_cse_env$uid), class = "fs_uneval")
}
.fs_is_uneval <- function(z) inherits(z, "fs_uneval")
# Environment names are limited to 10000 bytes: longer keys are interned.
.fs_intern <- utils::hashtab()
.fs_ek <- function(k) {
  if (nchar(k, "bytes") < 9000L) return(k)
  id <- utils::gethash(.fs_intern, k)
  if (is.null(id)) { id <- paste0("\001K", utils::numhash(.fs_intern) + 1L); utils::sethash(.fs_intern, k, id) }
  id
}
.fs_ukey <- function(z) if (.fs_is_uneval(z)) paste0("\001U", z$id) else .fs_ek(z$key)
.fs_uneval_basic <- function(z) {
  switch(z$func, Add = .add(z$args), Mul = .mul(z$args), Pow = .pow(z$args[[1L]], z$args[[2L]]))
}
.fs_cse_args <- function(z) {
  if (.fs_is_uneval(z)) return(z$args)
  if (z$type %in% c("Add", "Mul")) return(.sympy_args(z))
  z$args
}
.fs_cse_free <- function(z) {
  if (.fs_is_uneval(z)) return(unlist(lapply(z$args, .fs_cse_free), recursive = FALSE))
  free_symbols(z)
}

# ordered sets of integers: insertion order kept
.fs_os_diff <- function(a, b) a[!(a %in% b)]
.fs_os_inter <- function(a, b) a[a %in% b]
.fs_os_union <- function(a, b) c(a, b[!(b %in% a)])

# FuncArgTracker
.fs_tracker <- function(funcs) {
  tr <- new.env(parent = emptyenv())
  tr$value_numbers <- new.env(parent = emptyenv())
  tr$values <- list()
  tr$arg_to_funcset <- list()
  tr$func_to_argset <- list()
  for (fi in seq_along(funcs)) {
    s <- integer()
    for (a in .fs_cse_args(funcs[[fi]])) {
      n <- .fs_tr_number(tr, a)
      if (!(n %in% s)) s <- c(s, n)
      if (!(fi %in% tr$arg_to_funcset[[n]])) tr$arg_to_funcset[[n]] <- c(tr$arg_to_funcset[[n]], fi)
    }
    tr$func_to_argset[[fi]] <- s
  }
  tr
}
.fs_tr_number <- function(tr, value) {
  k <- .fs_ukey(value)
  n <- tr$value_numbers[[k]]
  if (!is.null(n)) return(n)
  n <- length(tr$values) + 1L
  assign(k, n, envir = tr$value_numbers)
  tr$values[[n]] <- value
  tr$arg_to_funcset[[n]] <- integer()
  n
}
.fs_tr_values <- function(tr, argset) tr$values[sort(argset)]
.fs_tr_stop <- function(tr, fi) {
  for (a in tr$func_to_argset[[fi]]) tr$arg_to_funcset[[a]] <- setdiff(tr$arg_to_funcset[[a]], fi)
}
.fs_tr_common <- function(tr, argset, min_func_i) {
  if (!length(argset)) return(integer())
  funcsets <- lapply(argset, function(a) tr$arg_to_funcset[[a]])
  li <- which.max(vapply(funcsets, length, 0L))
  cnt_k <- integer(); cnt_v <- integer()
  bump <- function(fi, by) {
    w <- match(fi, cnt_k)
    if (is.na(w)) { cnt_k <<- c(cnt_k, fi); cnt_v <<- c(cnt_v, by) } else cnt_v[w] <<- cnt_v[w] + by
  }
  for (j in seq_along(funcsets)) {
    if (j == li) next
    for (fi in funcsets[[j]]) if (fi >= min_func_i) bump(fi, 1L)
  }
  largest <- funcsets[[li]]
  if (length(largest) <= length(cnt_k)) {
    smaller <- largest
    for (fi in smaller) {
      w <- match(fi, cnt_k)
      if (is.na(w) || cnt_v[w] < 1L) next
      cnt_v[w] <- cnt_v[w] + 1L
    }
  } else {
    for (w in seq_along(cnt_k)) {
      if (cnt_v[w] < 1L) next
      if (cnt_k[w] %in% largest) cnt_v[w] <- cnt_v[w] + 1L
    }
  }
  keep <- cnt_v >= 2L
  k <- cnt_k[keep]; v <- cnt_v[keep]
  k[order(v, k)]
}
.fs_tr_subset <- function(tr, argset, restrict) {
  idx <- tr$arg_to_funcset[[argset[1L]]]
  if (!is.null(restrict)) idx <- .fs_os_inter(idx, restrict)
  for (a in argset[-1L]) idx <- .fs_os_inter(idx, tr$arg_to_funcset[[a]])
  idx
}
.fs_tr_update <- function(tr, fi, new_argset) {
  old <- tr$func_to_argset[[fi]]
  for (d in .fs_os_diff(old, new_argset)) tr$arg_to_funcset[[d]] <- setdiff(tr$arg_to_funcset[[d]], fi)
  for (a in .fs_os_diff(new_argset, old)) tr$arg_to_funcset[[a]] <- c(tr$arg_to_funcset[[a]], fi)
  tr$func_to_argset[[fi]] <- unique(new_argset)
}

.fs_match_common_args <- function(func_class, funcs, opt) {
  if (!length(funcs)) return(invisible())
  funcs <- funcs[order(vapply(funcs, function(f) length(f$args), 0L))]
  tr <- .fs_tracker(funcs)
  changed <- integer()
  for (i in seq_along(funcs)) {
    cands <- .fs_tr_common(tr, tr$func_to_argset[[i]], i + 1L)
    while (length(cands)) {
      j <- cands[1L]
      cands <- cands[-1L]
      com <- .fs_os_inter(tr$func_to_argset[[i]], tr$func_to_argset[[j]])
      if (length(com) <= 1L) next
      diff_i <- .fs_os_diff(tr$func_to_argset[[i]], com)
      if (length(diff_i)) {
        cf <- .fs_uneval(func_class, .fs_tr_values(tr, com))
        cn <- .fs_tr_number(tr, cf)
        .fs_tr_update(tr, i, .fs_os_union(diff_i, cn))
        changed <- union(changed, i)
      } else {
        cn <- .fs_tr_number(tr, funcs[[i]])
      }
      diff_j <- .fs_os_diff(tr$func_to_argset[[j]], com)
      .fs_tr_update(tr, j, .fs_os_union(diff_j, cn))
      changed <- union(changed, j)
      for (k in .fs_tr_subset(tr, com, cands)) {
        diff_k <- .fs_os_diff(tr$func_to_argset[[k]], com)
        .fs_tr_update(tr, k, .fs_os_union(diff_k, cn))
        changed <- union(changed, k)
      }
    }
    if (i %in% changed) assign(.fs_ek(funcs[[i]]$key), .fs_uneval(func_class, .fs_tr_values(tr, tr$func_to_argset[[i]])), envir = opt)
    .fs_tr_stop(tr, i)
  }
  invisible()
}

.fs_is_atomish <- function(z) !length(z$args) || z$type == "Order"
.fs_is_iter <- function(z) z$type == "Tuple"

.fs_opt_cse <- function(exprs) {
  opt <- new.env(parent = emptyenv())
  adds <- list(); muls <- list()
  seen <- new.env(parent = emptyenv())
  find <- function(z) {
    if (!is_sym(z) || .fs_is_atomish(z)) return(invisible())
    if (.fs_is_iter(z)) { lapply(z$args, find); return(invisible()) }
    if (exists(.fs_ek(z$key), envir = seen, inherits = FALSE)) return(invisible())
    assign(.fs_ek(z$key), TRUE, envir = seen)
    lapply(.fs_cse_args(z), find)
    if (!.fs_is_matrix_expr_node(z) && could_extract_minus_sign(z)) {
      ng <- if (z$type == "Add") .add(lapply(z$args, .neg)) else .neg(z)
      if (length(ng$args)) {
        assign(.fs_ek(z$key), .fs_uneval("Mul", list(S_NegativeOne, ng)), envir = opt)
        assign(.fs_ek(ng$key), TRUE, envir = seen)
        z <- ng
      }
    }
    if (z$type == "Mul") {
      if (is.null(muls[[z$key]])) muls[[z$key]] <<- z
    } else if (z$type == "Add") {
      if (is.null(adds[[z$key]])) adds[[z$key]] <<- z
    } else if (z$type == "Pow") {
      b <- z$args[[1L]]; ex <- z$args[[2L]]
      if (could_extract_minus_sign(ex)) assign(.fs_ek(z$key), .fs_uneval("Pow", list(.pow(b, .neg(ex)), S_NegativeOne)), envir = opt)
    }
    invisible()
  }
  for (e in exprs) if (is_sym(e)) find(e)
  cmuls <- list()
  for (m in muls) {
    comm <- vapply(m$args, .is_comm, TRUE)
    c <- m$args[comm]; nc <- m$args[!comm]
    if (length(c)) {
      c_mul <- .mul(c)
      if (length(nc)) {
        new_obj <- if (is_one_num(c_mul)) .mul(nc) else new_sym("Mul", list(c_mul, .mul(nc)))
        assign(.fs_ek(m$key), new_obj, envir = opt)
      }
      if (length(c) > 1L && is.null(cmuls[[c_mul$key]])) cmuls[[c_mul$key]] <- c_mul
    }
  }
  .fs_match_common_args("Add", unname(adds), opt)
  .fs_match_common_args("Mul", unname(cmuls), opt)
  opt
}

.fs_is_matrix_expr_node <- function(z) z$type %in% c("MatrixSymbol", "MatMul", "MatAdd", "MatPow", "Inverse", "Transpose", "MatrixElement")

.fs_tree_cse <- function(exprs, symbols, opt, order, ignore) {
  to_elim <- new.env(parent = emptyenv())
  seen <- new.env(parent = emptyenv())
  excluded <- character()
  ign_keys <- vapply(ignore, .k, "")
  optget <- function(z) if (.fs_is_uneval(z)) NULL else get0(.fs_ek(z$key), envir = opt, inherits = FALSE)
  find <- function(z) {
    if (!(is_sym(z) || .fs_is_uneval(z))) return(invisible())
    if (!.fs_is_uneval(z)) {
      if (z$type %in% c("RootOf", "CRootOf")) return(invisible())
      if (.fs_is_atomish(z) || z$type %in% c("MatrixSymbol", "MatrixElement")) {
        if (z$type == "Symbol") excluded <<- union(excluded, z$name)
        return(invisible())
      }
    }
    if (!.fs_is_uneval(z) && .fs_is_iter(z)) {
      args <- z$args
    } else {
      k <- .fs_ukey(z)
      if (exists(k, envir = seen, inherits = FALSE)) {
        if (!length(ign_keys) || !any(vapply(.fs_cse_free(z), .k, "") %in% ign_keys)) {
          assign(k, TRUE, envir = to_elim)
          return(invisible())
        }
      }
      assign(k, TRUE, envir = seen)
      o <- optget(z)
      if (!is.null(o)) z <- o
      args <- .fs_cse_args(z)
    }
    lapply(args, find)
    invisible()
  }
  for (e in exprs) if (is_sym(e)) find(e)
  sym_i <- 0L
  next_sym <- function() {
    repeat {
      sym_i <<- sym_i + 1L
      s <- if (is.function(symbols)) symbols(sym_i)
           else if (sym_i <= length(symbols)) symbols[[sym_i]]
           else stop(.sympy_condition("ValueError", "Symbols iterator ran out of symbols."))
      if (!(s$name %in% excluded)) return(s)
    }
  }
  subs_tab <- new.env(parent = emptyenv())
  replacements <- list()
  rebuild <- function(z) {
    if (!(is_sym(z) || .fs_is_uneval(z))) return(z)
    if (!length(z$args)) return(z)
    if (!.fs_is_uneval(z) && .fs_is_iter(z)) return(sym_rebuild(z, lapply(z$args, rebuild)))
    k <- .fs_ukey(z)
    hit <- get0(k, envir = subs_tab, inherits = FALSE)
    if (!is.null(hit)) return(hit)
    orig_k <- k
    o <- optget(z)
    if (!is.null(o)) z <- o
    if (order != "none" && !.fs_is_uneval(z) && z$type %in% c("Mul", "Add")) {
      if (z$type == "Mul") {
        comm <- vapply(z$args, .is_comm, TRUE)
        c <- z$args[comm]; nc <- z$args[!comm]
        args <- if (length(c) == 1L && is_one_num(c[[1L]])) nc else c(.fs_ordered(c), nc)
      } else args <- .fs_ordered(z$args)
    } else args <- .fs_cse_args(z)
    new_args <- lapply(args, rebuild)
    same <- identical(vapply(new_args, .fs_ukey, ""), vapply(args, .fs_ukey, ""))
    new_expr <- if (.fs_is_uneval(z)) .fs_uneval_basic(list(func = z$func, args = new_args))
                else if (!same) sym_rebuild(z, new_args) else z
    if (exists(orig_k, envir = to_elim, inherits = FALSE)) {
      s <- next_sym()
      assign(orig_k, s, envir = subs_tab)
      replacements[[length(replacements) + 1L]] <<- list(s, new_expr)
      return(s)
    }
    new_expr
  }
  reduced <- lapply(exprs, function(e) if (is_sym(e)) rebuild(e) else e)
  list(replacements = replacements, reduced = reduced)
}

# numbered_symbols(prefix, cls=Symbol): the i-th symbol of the stream
.fs_numbered <- function(prefix) {
  force(prefix)
  function(i) Symbol(paste0(prefix, i - 1L))
}

.fs_cse <- function(exprs, symbols = NULL, optimizations = NULL, postprocess = NULL,
                    order = "canonical", ignore = list(), list = TRUE, prefix = "x") {
  if (!isTRUE(list)) {
    single <- is_sym(exprs) || inherits(exprs, "symmatrix")
    r <- .fs_cse(if (single) base::list(exprs) else exprs, symbols, optimizations, postprocess = NULL,
                 order = order, ignore = ignore, prefix = prefix)
    if (single) r$reduced <- r$reduced[[1L]]
    return(if (is.null(postprocess)) r else postprocess(r$replacements, r$reduced))
  }
  if (is.numeric(exprs) || is_sym(exprs) || inherits(exprs, "symmatrix")) exprs <- base::list(exprs)
  exprs <- lapply(exprs, function(e) if (inherits(e, "symmatrix")) e else as_sym(e))
  if (is_sym(ignore)) ignore <- base::list(ignore)
  copy <- exprs
  flat <- lapply(exprs, function(e) if (inherits(e, "symmatrix")) do.call(Tuple, e$data) else e)
  if (is.null(optimizations)) optimizations <- base::list()
  if (identical(optimizations, "basic"))
    stop(.sympy_condition("NotImplementedError", "the basic optimizations of cse are not available"))
  for (op in optimizations) if (is.function(op[[1L]])) flat <- lapply(flat, op[[1L]])
  symbols <- if (is.null(symbols)) .fs_numbered(prefix) else if (is_sym(symbols)) base::list(symbols) else as.list(symbols)
  opt <- .fs_opt_cse(flat)
  r <- .fs_tree_cse(flat, symbols, opt, order, ignore)
  post <- function(z) {
    for (op in rev(optimizations)) if (length(op) > 1L && is.function(op[[2L]])) z <- op[[2L]](z)
    z
  }
  r$replacements <- lapply(r$replacements, function(p) base::list(p[[1L]], post(p[[2L]])))
  r$reduced <- lapply(r$reduced, post)
  for (i in seq_along(copy)) {
    if (inherits(copy[[i]], "symmatrix")) {
      m <- copy[[i]]
      m$data <- r$reduced[[i]]$args
      r$reduced[[i]] <- m
    }
  }
  if (is.null(postprocess)) return(r)
  postprocess(r$replacements, r$reduced)
}

.fs_ordered <- function(l) if (length(l) < 2L) l else .bool_print_order(l)
