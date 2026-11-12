# The Earley parser of lark (lark/parsers/xearley.py and earley.py, with the
# SPPF of earley_forest.py and the tree shaping of parse_tree_builder.py), as
# the LaTeX parser of sympy/parsing/latex/lark runs it: a dynamic lexer that
# matches the terminals expected at each character, explicit ambiguity, and
# keep_all_tokens. The grammar is the compiled one in parsing-latex-grammar.R.
# Item sets keep insertion order, as lark's OrderedSet does, so the readings of
# an ambiguous input come out in lark's order.

.lark_env <- new.env(parent = emptyenv())

.lark_grammar <- function() {
  if (!is.null(.lark_env$g)) return(.lark_env$g)
  R <- .LARK_RULES
  g <- new.env(parent = emptyenv())
  g$origin <- vapply(R, `[[`, "", "origin")
  g$exp <- lapply(R, `[[`, "expansion")
  g$isterm <- lapply(R, `[[`, "is_term")
  g$len <- lengths(g$exp)
  g$order <- vapply(R, `[[`, 0L, "order")
  g$expand1 <- vapply(R, `[[`, FALSE, "expand1")
  alias <- vapply(R, `[[`, "", "alias")
  g$name <- ifelse(is.na(alias), g$origin, alias)
  g$has_alias <- !is.na(alias)
  # children that are inlined into their parent: nonterminals named _rule
  g$inline <- lapply(seq_along(R), function(r) !g$isterm[[r]] & startsWith(g$exp[[r]], "_"))
  terms <- unique(unlist(Map(function(e, t) e[t], g$exp, g$isterm)))
  g$terminals <- terms
  g$origins <- unique(g$origin)
  g$pred <- .LARK_PREDICTIONS
  g$pred_term <- lapply(g$pred, function(rr) vapply(rr, function(r) g$len[r] > 0L && g$isterm[[r]][1L], TRUE))
  g$patterns <- .LARK_TERMINALS
  .lark_env$g <- g
  g
}

# lark's term matcher: re.compile(pattern).match(text, pos), so lookbehinds
# see the text before pos. Each terminal is matched at every position in one
# vectorized call; one character of left context is enough for the
# lookbehinds of the grammar.
.lark_matcher <- function(text) {
  cache <- new.env(parent = emptyenv())
  pats <- .lark_grammar()$patterns
  n <- nchar(text)
  subj <- if (n > 1L) substring(text, seq_len(n - 1L), n) else character()
  function(term, pos) {
    v <- cache[[term]]
    if (is.null(v)) {
      v <- rep(NA_integer_, n)
      if (n) {
        m0 <- regexpr(paste0("^(?:", pats[[term]], ")"), text, perl = TRUE)
        if (m0 > 0L && attr(m0, "match.length") > 0L) v[1L] <- attr(m0, "match.length")
        if (n > 1L) {
          m <- regexpr(paste0("^(?s:.)\\K(?:", pats[[term]], ")"), subj, perl = TRUE)
          ok <- m > 0L & attr(m, "match.length") > 0L
          v[-1L][ok] <- attr(m, "match.length")[ok]
        }
      }
      assign(term, v, envir = cache)
    }
    len <- v[pos + 1L]
    if (is.na(len)) NA_integer_ else pos + len
  }
}

.lark_set <- function() {
  s <- new.env(parent = emptyenv())
  s$ids <- integer()
  s$keys <- character()
  s
}

.lark_unexpected <- function(cls, message) {
  structure(class = c(cls, "UnexpectedInput", "LarkError", "sympy_error", "error", "condition"),
            list(message = message, call = NULL))
}

.lark_context <- function(chars, pos) {
  n <- length(chars)
  before <- if (pos > 0L) paste(chars[max(pos - 40L, 0L) + seq_len(pos - max(pos - 40L, 0L))], collapse = "") else ""
  before <- sub("^.*\n", "", before)
  after <- if (pos < n) paste(chars[(pos + 1L):min(pos + 40L, n)], collapse = "") else ""
  after <- sub("\n.*$", "", after)
  paste0(before, after, "\n", strrep(" ", nchar(before)), "^\n")
}

.lark_expected <- function(names) paste0("Expected one of: \n\t* ", paste(names, collapse = "\n\t* "), "\n")

# Parse `text` and return the parse tree, with an "_ambig" node at the root
# when the input has several readings.
.lark_parse <- function(text, start = "latex_string") {
  g <- .lark_grammar()
  chars <- strsplit(text, "", fixed = TRUE)[[1L]]
  nchar_ <- length(chars)
  match <- .lark_matcher(text)
  is_terminal <- function(sym) sym %in% g$terminals

  # items: rule, ptr, start, node
  it_rule <- integer(); it_ptr <- integer(); it_start <- integer(); it_node <- integer()
  new_item <- function(rule, ptr, st, node = 0L) {
    k <- length(it_rule) + 1L
    it_rule[k] <<- rule; it_ptr[k] <<- ptr; it_start[k] <<- st; it_node[k] <<- node
    k
  }
  item_key <- function(id) paste(it_rule[id], it_ptr[id], it_start[id], sep = ".")
  expect_of <- function(id) {
    r <- it_rule[id]; p <- it_ptr[id]
    if (p == g$len[r]) NA_character_ else g$exp[[r]][p + 1L]
  }
  set_add <- function(s, id) {
    k <- item_key(id)
    if (k %in% s$keys) return(FALSE)
    s$keys <- c(s$keys, k)
    s$ids <- c(s$ids, id)
    TRUE
  }
  # add the predicted items Item(rule, 0, i) for the rules rr, in order;
  # returns the ids added
  set_add_new <- function(s, rr, i) {
    if (!length(rr)) return(integer())
    keys <- paste(rr, 0L, i, sep = ".")
    new <- !(keys %in% s$keys)
    if (!any(new)) return(integer())
    rr <- rr[new]
    k0 <- length(it_rule)
    ids <- k0 + seq_along(rr)
    it_rule[ids] <<- rr; it_ptr[ids] <<- 0L; it_start[ids] <<- i; it_node[ids] <<- 0L
    s$keys <- c(s$keys, keys[new])
    s$ids <- c(s$ids, ids)
    ids
  }

  # SPPF: symbol nodes, packed nodes and token nodes
  sn_inter <- logical(); sn_children <- list(); sn_pkeys <- list()
  pk_rule <- integer(); pk_left <- integer(); pk_right <- integer(); pk_parent <- integer()
  tk_type <- character(); tk_value <- character(); tk_start <- integer()
  node_for <- function(cache, id, end) {
    r <- it_rule[id]; p <- it_ptr[id]
    complete <- p == g$len[r]
    key <- if (complete) paste0("o", g$origin[r], ":", it_start[id]) else paste0("r", r, ".", p, ":", it_start[id])
    nd <- cache[[key]]
    if (is.null(nd)) {
      nd <- length(sn_inter) + 1L
      sn_inter[nd] <<- !complete
      sn_children[nd] <<- list(integer())
      sn_pkeys[nd] <<- list(character())
      assign(key, nd, envir = cache)
    }
    nd
  }
  add_family <- function(nd, rule, left, right) {
    k <- paste(left, right)
    if (k %in% sn_pkeys[[nd]]) return(invisible())
    p <- length(pk_rule) + 1L
    pk_rule[p] <<- rule; pk_left[p] <<- left; pk_right[p] <<- right; pk_parent[p] <<- nd
    sn_pkeys[[nd]] <<- c(sn_pkeys[[nd]], k)
    sn_children[[nd]] <<- c(sn_children[[nd]], p)
    invisible()
  }

  columns <- list(.lark_set())
  to_scan <- .lark_set()
  orig_cache <- list()
  originators <- function(col, sym) {
    key <- paste0(col, ":", sym)
    v <- orig_cache[[key]]
    if (is.null(v)) {
      ids <- columns[[col + 1L]]$ids
      ex <- vapply(ids, expect_of, "")
      v <- ids[!is.na(ex) & ex == sym]
      orig_cache[[key]] <<- v
    }
    v
  }
  for (r in g$pred[[start]]) {
    id <- new_item(r, 0L, 0L)
    if (is_terminal(expect_of(id))) set_add(to_scan, id) else set_add(columns[[1L]], id)
  }

  predict_and_complete <- function(i, node_cache) {
    column <- columns[[i + 1L]]
    stack <- column$ids
    while (length(stack)) {
      id <- stack[length(stack)]
      stack <- stack[-length(stack)]
      r <- it_rule[id]
      if (it_ptr[id] == g$len[r]) {
        sym <- g$origin[r]
        for (o in originators(it_start[id], sym)) {
          ni <- new_item(it_rule[o], it_ptr[o] + 1L, it_start[o])
          nd <- node_for(node_cache, ni, i)
          it_node[ni] <<- nd
          add_family(nd, it_rule[o], it_node[o], it_node[id])
          if (is_terminal(expect_of(ni))) set_add(to_scan, ni)
          else if (set_add(column, ni)) stack <- c(stack, ni)
        }
      } else {
        ex <- expect_of(id)
        if (!is_terminal(ex)) {
          rr <- g$pred[[ex]]; ft <- g$pred_term[[ex]]
          set_add_new(to_scan, rr[ft], i)
          stack <- c(stack, set_add_new(column, rr[!ft], i))
        }
      }
    }
  }

  delayed <- new.env(parent = emptyenv())
  push_delayed <- function(end, ids, toks) {
    k <- as.character(end)
    old <- delayed[[k]]
    if (is.null(old)) old <- list(ids = integer(), toks = integer())
    assign(k, list(ids = c(old$ids, ids), toks = c(old$toks, toks)), envir = delayed)
  }
  line <- 1L; col <- 1L
  scan <- function(i) {
    node_cache <- new.env(parent = emptyenv())
    for (id in to_scan$ids) {
      term <- expect_of(id)
      e <- match(term, i)
      if (!is.na(e)) {
        tk <- length(tk_type) + 1L
        tk_type[tk] <<- term; tk_value[tk] <<- substr(text, i + 1L, e); tk_start[tk] <<- i
        push_delayed(e, id, tk)
      }
    }
    for (x in .LARK_IGNORE) {
      e <- match(x, i)
      if (!is.na(e)) {
        push_delayed(e, to_scan$ids, rep(0L, length(to_scan$ids)))
        cids <- columns[[i + 1L]]$ids
        done <- cids[vapply(cids, function(c) it_ptr[c] == g$len[it_rule[c]] && g$origin[it_rule[c]] == start, TRUE)]
        push_delayed(e, done, rep(0L, length(done)))
      }
    }
    next_to_scan <- .lark_set()
    next_set <- .lark_set()
    columns[[i + 2L]] <<- next_set
    k <- as.character(i + 1L)
    d <- delayed[[k]]
    if (!is.null(d)) {
      for (j in seq_along(d$ids)) {
        id <- d$ids[j]; tk <- d$toks[j]
        if (tk) {
          ni <- new_item(it_rule[id], it_ptr[id] + 1L, it_start[id])
          nd <- node_for(node_cache, ni, i + 1L)
          it_node[ni] <<- nd
          add_family(nd, it_rule[id], it_node[id], -tk)
        } else ni <- id
        if (is_terminal(expect_of(ni))) set_add(next_to_scan, ni) else set_add(next_set, ni)
      }
      rm(list = k, envir = delayed)
    }
    if (!length(next_set$ids) && !length(ls(delayed)) && !length(next_to_scan$ids)) {
      allowed <- sort(unique(vapply(to_scan$ids, expect_of, "")), method = "radix")
      msg <- paste0(sprintf("No terminal matches '%s' in the current parser context, at line %d col %d",
                            chars[i + 1L], line, col),
                    "\n\n", .lark_context(chars, i), .lark_expected(allowed))
      stop(.lark_unexpected("UnexpectedCharacters", msg))
    }
    list(next_to_scan, node_cache)
  }

  node_cache <- new.env(parent = emptyenv())
  i <- 0L
  for (ch in chars) {
    predict_and_complete(i, node_cache)
    res <- scan(i)
    to_scan <- res[[1L]]; node_cache <- res[[2L]]
    if (ch == "\n") { line <- line + 1L; col <- 1L } else col <- col + 1L
    i <- i + 1L
  }
  predict_and_complete(i, node_cache)

  last <- columns[[i + 1L]]$ids
  sol <- unique(it_node[last[vapply(last, function(c) it_ptr[c] == g$len[it_rule[c]] &&
                                      g$origin[it_rule[c]] == start && it_start[c] == 0L && it_node[c] > 0L, TRUE)]])
  if (!length(sol)) {
    expected <- vapply(to_scan$ids, expect_of, "")
    stop(.lark_unexpected("UnexpectedEOF", paste0("Unexpected end-of-input. ", .lark_expected(expected))))
  }
  forest <- list(sn_inter = sn_inter, sn_children = sn_children, pk_rule = pk_rule, pk_left = pk_left,
                 pk_right = pk_right, pk_parent = pk_parent, tk_type = tk_type, tk_value = tk_value,
                 tk_start = tk_start)
  .lark_forest_to_tree(forest, sol[1L], g)
}

# lark trees: list(data, children); tokens: list(type, value, start)
.lark_tree <- function(data, children) structure(list(data = data, children = children), class = "lark_tree")
.lark_token <- function(type, value, start = NA_integer_) structure(list(type = type, value = value, start = start), class = "lark_token")
.is_tree <- function(x, data = NULL) inherits(x, "lark_tree") && (is.null(data) || x$data %in% data)

# ForestToParseTree with resolve_ambiguity=False, then the callbacks that
# ParseTreeBuilder wraps with ChildFilter, AmbiguousExpander and
# AmbiguousIntermediateExpander. The grammar has no cycles and no priorities,
# so packed nodes sort by (is_empty, rule.order) and keep insertion order.
.lark_forest_to_tree <- function(f, root, g) {
  cache <- list()
  tr_node <- function(ref) {
    if (ref < 0L) return(.lark_token(f$tk_type[-ref], f$tk_value[-ref], f$tk_start[-ref]))
    pks <- f$sn_children[[ref]]
    ord <- order(f$pk_left[pks] == 0L & f$pk_right[pks] == 0L, g$order[f$pk_rule[pks]], seq_along(pks))
    data <- lapply(pks[ord], tr_packed)
    if (f$sn_inter[ref]) {
      if (length(data) > 1L)
        return(.lark_tree("_iambig", lapply(data, function(c) .lark_tree("_inter", c))))
      return(data[[1L]])
    }
    flat <- list()
    for (d in data) flat <- if (.is_tree(d, "_ambig")) c(flat, d$children) else c(flat, list(d))
    if (length(flat) > 1L) .lark_tree("_ambig", flat) else flat[[1L]]
  }
  tr_packed <- function(p) {
    key <- as.character(p)
    if (!is.null(cache[[key]])) return(cache[[key]])
    children <- list()
    l <- f$pk_left[p]; r <- f$pk_right[p]
    if (l != 0L) {
      lv <- tr_node(l)
      if (l > 0L && f$sn_inter[l] && !inherits(lv, c("lark_tree", "lark_token"))) children <- c(children, lv)
      else children <- c(children, list(lv))
    }
    if (r != 0L) children <- c(children, list(tr_node(r)))
    out <- if (f$sn_inter[f$pk_parent[p]]) children else .lark_callback(f$pk_rule[p], children, g)
    cache[[key]] <<- out
    out
  }
  tr_node(root)
}

.lark_callback <- function(rule, children, g) {
  # AmbiguousIntermediateExpander
  collapse <- function(ch) {
    if (length(ch) && .is_tree(ch[[1L]], "_iambig")) {
      res <- list()
      rest <- ch[-1L]
      for (gc in ch[[1L]]$children) {
        col <- collapse(gc$children)
        if (!is.null(col)) res <- c(res, lapply(col, function(t) .lark_tree("_inter", c(t$children, rest))))
        else res <- c(res, list(.lark_tree("_inter", c(gc$children, rest))))
      }
      return(res)
    }
    NULL
  }
  collapsed <- collapse(children)
  if (!is.null(collapsed))
    return(.lark_tree("_ambig", lapply(collapsed, function(t) .lark_ambig_expand(rule, t$children, g))))
  .lark_ambig_expand(rule, children, g)
}

.lark_ambig_expand <- function(rule, children, g) {
  amb <- integer()
  for (i in seq_along(children)) {
    ch <- children[[i]]
    if (.is_tree(ch, "_ambig")) {
      amb <- c(amb, i)
      kids <- list()
      for (k in ch$children) kids <- if (.is_tree(k, "_ambig")) c(kids, k$children) else c(kids, list(k))
      children[[i]] <- .lark_tree("_ambig", kids)
    }
  }
  if (!length(amb)) return(.lark_child_filter(rule, children, g))
  expand <- lapply(seq_along(children), function(i) if (i %in% amb) children[[i]]$children else list(children[[i]]))
  # itertools.product: the last position varies fastest
  combos <- list(list())
  for (e in expand) combos <- unlist(lapply(combos, function(cmb) lapply(e, function(v) c(cmb, list(v)))), recursive = FALSE)
  .lark_tree("_ambig", lapply(combos, function(cmb) .lark_child_filter(rule, cmb, g)))
}

.lark_child_filter <- function(rule, children, g) {
  inl <- g$inline[[rule]]
  out <- list()
  for (i in seq_along(children)) {
    out <- if (inl[i]) c(out, children[[i]]$children) else c(out, list(children[[i]]))
  }
  if (g$expand1[rule] && !g$has_alias[rule] && length(out) == 1L) return(out[[1L]])
  .lark_tree(g$name[rule], out)
}
