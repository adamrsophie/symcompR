# SMT-LIB parser, after sympy/parsing/smtlib (the lark grammar
# smtlib.lark, LarkSMTLibParser and SMTLibTransformer). The text is first
# parsed completely, so syntax errors come before any other error, and the
# commands are then evaluated in order.

.SMT_COMMANDS <- c("assert", "check-sat", "check-sat-assuming", "declare-const", "declare-fun",
  "declare-sort", "define-const", "define-fun", "define-fun-rec", "define-sort", "declare-datatype",
  "declare-datatypes", "echo", "exit", "get-assertions", "get-assignment", "get-info", "get-model",
  "get-option", "get-proof", "get-unsat-assumptions", "get-unsat-core", "get-value", "pop", "push",
  "reset", "reset-assertions", "set-info", "set-logic", "set-option")

.SMT_UNSUPPORTED <- c("concat", "bvnot", "bvand", "bvor", "bvneg", "bvadd", "bvmul",
  "bvudiv", "bvurem", "bvshl", "bvlshr", "bvsub", "bvult", "bvxor",
  "bvnand", "bvnor", "bvxnor", "bvcomp", "bvsdiv", "bvsrem",
  "bvsmod", "bvashr", "bvule", "bvugt", "bvuge", "bvslt", "bvsle",
  "bvsgt", "bvsge", "bv2nat", "bv2int", "ubv_to_int",
  "ext_rotate_left", "ext_rotate_right", "str.len", "str.++",
  "str.at", "str.contains", "str.indexof", "str.replace",
  "str.substr", "str.prefixof", "str.suffixof", "str.to.int",
  "int.to.str", "select", "store", "const")

.smt_syntax <- function(msg, st) {
  before <- substr(st$s, 1L, st$pos - 1L)
  line <- lengths(regmatches(before, gregexpr("\n", before))) + 1L
  col <- st$pos - max(c(0L, gregexpr("\n", before)[[1L]]))
  stop(.sympy_condition("SMTLibSyntaxError",
                        sprintf("Syntax error: %s at line %d, column %d.", msg, line, col)))
}

# the lexer: whitespace and comments are skipped; `num` tells whether a
# numeral or decimal is acceptable here, as the contextual lexer of lark
# only tries the terminals that the parser state accepts
.smt_lex <- function(st, num) {
  rest <- substring(st$s, st$pos)
  sk <- regexpr("^([ \t\f\r\n]+|;[^\n]*)+", rest, perl = TRUE)
  if (sk > 0L) {
    st$pos <- st$pos + attr(sk, "match.length")
    rest <- substring(rest, attr(sk, "match.length") + 1L)
  }
  if (!nzchar(rest)) return(list(type = "EOF", value = "", pos = st$pos))
  tok <- function(type, n) list(type = type, value = substr(rest, 1L, n), pos = st$pos, len = n)
  m <- function(p) { r <- regexpr(p, rest, perl = TRUE); if (r > 0L) attr(r, "match.length") else 0L }
  c1 <- substr(rest, 1L, 1L)
  if (c1 == "(") return(tok("LPAR", 1L))
  if (c1 == ")") return(tok("RPAR", 1L))
  if (num) {
    n <- m("^-?(0|[1-9][0-9]*)\\.[0-9]+")
    if (n) return(tok("DECIMAL", n))
    n <- m("^-?(0|[1-9][0-9]*)")
    if (n) return(tok("NUMERAL", n))
  }
  if (c1 == "\"") {
    # "([^"]|"")*" with the backtracking of a regular expression
    chars <- strsplit(rest, "")[[1L]]
    i <- 2L; last_pair <- 0L
    while (i <= length(chars)) {
      if (chars[i] == "\"") {
        if (i < length(chars) && chars[i + 1L] == "\"") { last_pair <- i; i <- i + 2L; next }
        return(tok("STRING", i))
      }
      i <- i + 1L
    }
    if (last_pair) return(tok("STRING", last_pair))
  }
  n <- m("^#x[0-9a-fA-F]+")
  if (n) return(tok("HEXADECIMAL", n))
  n <- m("^#b[01]+")
  if (n) return(tok("BINARY", n))
  n <- m("^\\|[^|\\\\]*\\|")
  if (n) return(tok("QUOTED_SYMBOL", n))
  n <- m("^:[a-zA-Z0-9_+*/%?!$&~^<>=@.-]+")
  if (n) return(tok("KEYWORD", n))
  n <- m("^[a-zA-Z_+*/%?!$&~^<>=@.-][a-zA-Z0-9_+*/%?!$&~^<>=@.-]*")
  if (n) return(tok("SYMBOL", n))
  .smt_syntax(sprintf("No terminal matches '%s'", c1), st)
}

.smt_peek <- function(st, num = TRUE) {
  p <- st$pos
  t <- .smt_lex(st, num)
  st$pos <- p
  t
}

.smt_next <- function(st, num = TRUE) {
  t <- .smt_lex(st, num)
  if (t$type != "EOF") st$pos <- t$pos + t$len
  t
}

.smt_unexpected <- function(st, t) {
  st$pos <- t$pos
  if (t$type == "EOF") .smt_syntax("Unexpected end of input", st)
  .smt_syntax(sprintf("Unexpected token '%s'", t$value), st)
}

.smt_expect <- function(st, type, num = FALSE) {
  t <- .smt_next(st, num)
  if (!(t$type %in% type)) .smt_unexpected(st, t)
  t
}

.smt_symbol <- function(st) {
  t <- .smt_expect(st, c("SYMBOL", "QUOTED_SYMBOL"))
  if (t$type == "QUOTED_SYMBOL") substr(t$value, 2L, nchar(t$value) - 1L) else t$value
}

.smt_numeral <- function(st) .smt_const(.smt_expect(st, "NUMERAL", num = TRUE))

# the token callbacks of SMTLibTransformer
.smt_const <- function(t) {
  v <- t$value
  switch(t$type,
    NUMERAL = make_rational(gmp::as.bigq(gmp::as.bigz(sub("^(-?)0+(?=[0-9])", "\\1", v, perl = TRUE)))),
    DECIMAL = Float(v),
    HEXADECIMAL = make_rational(gmp::as.bigq(gmp::as.bigz(paste0("0x", substring(v, 3L))))),
    BINARY = make_rational(gmp::as.bigq(gmp::as.bigz(paste0("0b", substring(v, 3L))))),
    STRING = Symbol(gsub("\"\"", "\"", substr(v, 2L, nchar(v) - 1L), fixed = TRUE)))
}

.SMT_SPEC <- c("NUMERAL", "DECIMAL", "HEXADECIMAL", "BINARY", "STRING")

.smt_rpar <- function(st) invisible(.smt_expect(st, "RPAR"))

.smt_sort <- function(st) {
  t <- .smt_next(st, FALSE)
  if (t$type == "SYMBOL") return(t$value)
  if (t$type == "QUOTED_SYMBOL") return(substr(t$value, 2L, nchar(t$value) - 1L))
  if (t$type != "LPAR") .smt_unexpected(st, t)
  h <- .smt_peek(st, FALSE)
  if (h$type == "SYMBOL" && h$value == "_") {
    .smt_next(st, FALSE)
    s <- .smt_symbol(st)
    nums <- list(.smt_numeral(st))
    while (.smt_peek(st)$type == "NUMERAL") nums[[length(nums) + 1L]] <- .smt_numeral(st)
    .smt_rpar(st)
    return(c(list("_", s), nums))
  }
  s <- .smt_symbol(st)
  out <- list(s, .smt_sort(st))
  while (.smt_peek(st, FALSE)$type != "RPAR") out[[length(out) + 1L]] <- .smt_sort(st)
  .smt_rpar(st)
  out
}

.smt_sexpr <- function(st) {
  t <- .smt_next(st, TRUE)
  if (t$type %in% .SMT_SPEC) return(.smt_const(t))
  if (t$type %in% c("SYMBOL", "QUOTED_SYMBOL", "KEYWORD")) return(t$value)
  if (t$type != "LPAR") .smt_unexpected(st, t)
  out <- list()
  while (.smt_peek(st)$type != "RPAR") out[[length(out) + 1L]] <- .smt_sexpr(st)
  .smt_rpar(st)
  out
}

# attribute or option: keyword | keyword attribute_value
.smt_attribute <- function(st) {
  k <- .smt_expect(st, "KEYWORD")$value
  t <- .smt_peek(st)
  if (t$type %in% c("RPAR", "KEYWORD")) return(list(k))
  if (t$type %in% c(.SMT_SPEC, "SYMBOL", "QUOTED_SYMBOL")) {
    .smt_next(st)
    return(list(k, if (t$type %in% .SMT_SPEC) .smt_const(t) else t$value))
  }
  if (t$type == "LPAR") return(list(k, .smt_sexpr(st)))
  .smt_unexpected(st, t)
}

.smt_sorted_var <- function(st) {
  .smt_expect(st, "LPAR")
  v <- .smt_symbol(st)
  s <- .smt_sort(st)
  .smt_rpar(st)
  list(v, s)
}

.smt_indexed <- function(st) {
  s <- .smt_symbol(st)
  .smt_numeral(st)
  while (.smt_peek(st)$type == "NUMERAL") .smt_numeral(st)
  .smt_rpar(st)
  s
}

# the qualified identifier after "(" "_" or "(" "as"
.smt_qual_paren <- function(st, head) {
  if (head == "_") return(list(k = "indexed", name = .smt_indexed(st)))
  t <- .smt_peek(st, FALSE)
  if (t$type == "LPAR") {
    .smt_next(st)
    u <- .smt_next(st, FALSE)
    if (!(u$type == "SYMBOL" && u$value == "_")) .smt_unexpected(st, u)
    name <- .smt_indexed(st)
    .smt_sort(st)
    .smt_rpar(st)
    return(list(k = "as_indexed", name = name))
  }
  name <- .smt_symbol(st)
  .smt_sort(st)
  .smt_rpar(st)
  list(k = "as", name = name)
}

.smt_terms_until_rpar <- function(st, min = 1L) {
  out <- list(.smt_term(st))
  while (.smt_peek(st)$type != "RPAR") out[[length(out) + 1L]] <- .smt_term(st)
  .smt_rpar(st)
  out
}

.smt_term <- function(st) {
  t <- .smt_next(st, TRUE)
  if (t$type %in% .SMT_SPEC) return(list(k = "const", v = .smt_const(t)))
  if (t$type %in% c("SYMBOL", "QUOTED_SYMBOL")) {
    name <- if (t$type == "QUOTED_SYMBOL") substr(t$value, 2L, nchar(t$value) - 1L) else t$value
    return(list(k = "qid", q = list(k = "simple", name = name)))
  }
  if (t$type != "LPAR") .smt_unexpected(st, t)
  h <- .smt_next(st, FALSE)
  if (h$type == "SYMBOL" && h$value %in% c("let", "forall", "exists", "match", "!", "_", "as")) {
    switch(h$value,
      let = {
        .smt_expect(st, "LPAR")
        bs <- list()
        repeat {
          .smt_expect(st, "LPAR")
          v <- .smt_symbol(st)
          bs[[length(bs) + 1L]] <- list(v, .smt_term(st))
          .smt_rpar(st)
          if (.smt_peek(st)$type == "RPAR") break
        }
        .smt_rpar(st)
        body <- .smt_term(st)
        .smt_rpar(st)
        return(list(k = "let", b = bs, body = body))
      },
      forall = , exists = {
        .smt_expect(st, "LPAR")
        repeat {
          .smt_sorted_var(st)
          if (.smt_peek(st)$type == "RPAR") break
        }
        .smt_rpar(st)
        body <- .smt_term(st)
        .smt_rpar(st)
        return(list(k = h$value, body = body))
      },
      match = {
        e <- .smt_term(st)
        .smt_expect(st, "LPAR")
        cases <- list()
        repeat {
          .smt_expect(st, "LPAR")
          p <- .smt_next(st, FALSE)
          if (p$type == "LPAR") {
            .smt_symbol(st); .smt_symbol(st)
            while (.smt_peek(st, FALSE)$type != "RPAR") .smt_symbol(st)
            .smt_rpar(st)
          } else if (!(p$type %in% c("SYMBOL", "QUOTED_SYMBOL"))) .smt_unexpected(st, p)
          cases[[length(cases) + 1L]] <- .smt_term(st)
          .smt_rpar(st)
          if (.smt_peek(st)$type == "RPAR") break
        }
        .smt_rpar(st)
        .smt_rpar(st)
        return(list(k = "match", e = e, cases = cases))
      },
      "!" = {
        body <- .smt_term(st)
        .smt_attribute(st)
        while (.smt_peek(st)$type != "RPAR") .smt_attribute(st)
        .smt_rpar(st)
        return(list(k = "annotate", body = body))
      },
      return(list(k = "qid", q = .smt_qual_paren(st, h$value))))
  }
  if (h$type == "LPAR") {
    u <- .smt_next(st, FALSE)
    if (!(u$type == "SYMBOL" && u$value %in% c("_", "as"))) .smt_unexpected(st, u)
    q <- .smt_qual_paren(st, u$value)
  } else if (h$type %in% c("SYMBOL", "QUOTED_SYMBOL")) {
    q <- list(k = "simple", name = if (h$type == "QUOTED_SYMBOL") substr(h$value, 2L, nchar(h$value) - 1L) else h$value)
  } else .smt_unexpected(st, h)
  list(k = "app", q = q, args = .smt_terms_until_rpar(st))
}

.smt_constructor <- function(st) {
  .smt_expect(st, "LPAR")
  .smt_symbol(st)
  while (.smt_peek(st, FALSE)$type == "LPAR") {
    .smt_next(st)
    .smt_symbol(st); .smt_sort(st)
    .smt_rpar(st)
  }
  .smt_rpar(st)
}

.smt_datatype_dec <- function(st) {
  .smt_expect(st, "LPAR")
  t <- .smt_peek(st, FALSE)
  if (t$type == "SYMBOL" && t$value == "par") {
    .smt_next(st)
    .smt_expect(st, "LPAR")
    .smt_symbol(st)
    while (.smt_peek(st, FALSE)$type != "RPAR") .smt_symbol(st)
    .smt_rpar(st)
    .smt_expect(st, "LPAR")
    .smt_constructor(st)
    while (.smt_peek(st, FALSE)$type != "RPAR") .smt_constructor(st)
    .smt_rpar(st)
    .smt_rpar(st)
    return(invisible())
  }
  .smt_constructor(st)
  while (.smt_peek(st, FALSE)$type != "RPAR") .smt_constructor(st)
  .smt_rpar(st)
}

.smt_command <- function(st) {
  .smt_expect(st, "LPAR")
  h <- .smt_next(st, FALSE)
  if (!(h$type %in% c("SYMBOL", "QUOTED_SYMBOL"))) .smt_unexpected(st, h)
  if (h$type == "QUOTED_SYMBOL" || !(h$value %in% .SMT_COMMANDS)) {
    name <- if (h$type == "QUOTED_SYMBOL") substr(h$value, 2L, nchar(h$value) - 1L) else h$value
    while (.smt_peek(st)$type != "RPAR") .smt_sexpr(st)
    .smt_rpar(st)
    return(list(cmd = "generic", name = name))
  }
  cmd <- h$value
  out <- list(cmd = cmd)
  switch(cmd,
    assert = { out$term <- .smt_term(st) },
    "check-sat-assuming" = {
      .smt_expect(st, "LPAR")
      while (.smt_peek(st, FALSE)$type != "RPAR") {
        t <- .smt_next(st, FALSE)
        if (t$type == "LPAR") {
          u <- .smt_next(st, FALSE)
          if (!(u$type == "SYMBOL" && u$value == "not")) .smt_unexpected(st, u)
          .smt_symbol(st)
          .smt_rpar(st)
        } else if (!(t$type %in% c("SYMBOL", "QUOTED_SYMBOL"))) .smt_unexpected(st, t)
      }
      .smt_rpar(st)
    },
    "declare-const" = { out$name <- .smt_symbol(st); out$sort <- .smt_sort(st) },
    "declare-fun" = {
      out$name <- .smt_symbol(st)
      .smt_expect(st, "LPAR")
      dom <- list()
      while (.smt_peek(st, FALSE)$type != "RPAR") dom[[length(dom) + 1L]] <- .smt_sort(st)
      .smt_rpar(st)
      out$domain <- dom
      out$sort <- .smt_sort(st)
    },
    "declare-sort" = { out$name <- .smt_symbol(st); out$arity <- .smt_numeral(st) },
    "define-const" = { out$name <- .smt_symbol(st); out$sort <- .smt_sort(st); out$term <- .smt_term(st) },
    "define-fun" = , "define-fun-rec" = {
      out$name <- .smt_symbol(st)
      .smt_expect(st, "LPAR")
      vars <- list()
      while (.smt_peek(st, FALSE)$type != "RPAR") vars[[length(vars) + 1L]] <- .smt_sorted_var(st)
      .smt_rpar(st)
      out$vars <- vars
      out$sort <- .smt_sort(st)
      out$term <- .smt_term(st)
    },
    "define-sort" = {
      out$name <- .smt_symbol(st)
      .smt_expect(st, "LPAR")
      ps <- character()
      while (.smt_peek(st, FALSE)$type != "RPAR") ps <- c(ps, .smt_symbol(st))
      .smt_rpar(st)
      out$params <- ps
      out$sort <- .smt_sort(st)
    },
    "declare-datatype" = { .smt_symbol(st); .smt_datatype_dec(st) },
    "declare-datatypes" = {
      .smt_expect(st, "LPAR")
      while (.smt_peek(st, FALSE)$type != "RPAR") {
        .smt_expect(st, "LPAR"); .smt_symbol(st); .smt_numeral(st); .smt_rpar(st)
      }
      .smt_rpar(st)
      .smt_expect(st, "LPAR")
      while (.smt_peek(st, FALSE)$type != "RPAR") .smt_datatype_dec(st)
      .smt_rpar(st)
    },
    echo = .smt_expect(st, "STRING"),
    "get-info" = , "get-option" = .smt_expect(st, "KEYWORD"),
    "get-value" = {
      .smt_expect(st, "LPAR")
      .smt_terms_until_rpar(st)
    },
    pop = , push = { out$n <- .smt_numeral(st) },
    "set-info" = , "set-option" = .smt_attribute(st),
    "set-logic" = { out$logic <- .smt_symbol(st) },
    NULL)
  .smt_rpar(st)
  out
}

.smt_parse <- function(s) {
  st <- new.env()
  st$s <- s; st$pos <- 1L
  cmds <- list()
  repeat {
    t <- .smt_peek(st, FALSE)
    if (t$type == "EOF") break
    cmds[[length(cmds) + 1L]] <- .smt_command(st)
  }
  cmds
}

# ---------------------------------------------------------------------------
# SMTLibTransformer

.smt_state <- function(tr = new.env()) {
  tr$symbols <- list(); tr$functions <- character(); tr$macros <- list()
  tr$assertions <- list(); tr$sorts <- list(); tr$sort_aliases <- list()
  tr$stack <- list()
  tr
}

.smt_resolve_sort <- function(tr, sort) {
  seen <- character()
  while (is.character(sort) && !is.null(tr$sort_aliases[[sort]]) && !(sort %in% seen)) {
    seen <- c(seen, sort)
    sort <- tr$sort_aliases[[sort]]
  }
  sort
}

.smt_indexed_identifier <- function(name) {
  if (startsWith(name, "bv") && grepl("^[0-9]+$", substring(name, 3L))) return(Integer(substring(name, 3L)))
  stop(.sympy_condition("NotImplementedError", sprintf("Indexed identifier (_ %s ...) is not supported", name)))
}

# the transformed qual_identifier: a name, or an expression for (_ bvN w)
.smt_qual <- function(q) {
  if (q$k %in% c("indexed", "as_indexed")) .smt_indexed_identifier(q$name) else q$name
}

.smt_xreplace <- function(e, keys, vals) {
  # a dict built from the pairs: the last value of a repeated key wins
  keep <- !duplicated(keys, fromLast = TRUE)
  if (!any(keep)) return(e)
  xreplace(e, unname(Map(function(k, v) list(Symbol(k), v), keys[keep], vals[keep])))
}

.smt_chain <- function(f, ps) {
  if (length(ps) == 2L) return(f(ps[[1L]], ps[[2L]]))
  do.call(And, lapply(seq_len(length(ps) - 1L), function(i) f(ps[[i]], ps[[i + 1L]])))
}

.smt_eval <- function(tr, t) {
  switch(t$k,
    const = t$v,
    qid = {
      atom <- .smt_qual(t$q)
      if (is_sym(atom)) return(atom)
      if (atom %in% names(tr$symbols)) return(tr$symbols[[atom]])
      if (atom == "true") return(S_true)
      if (atom == "false") return(S_false)
      Symbol(atom)
    },
    app = {
      opn <- .smt_qual(t$q)
      ps <- lapply(t$args, .smt_eval, tr = tr)
      .smt_apply(tr, opn, ps)
    },
    let = {
      vals <- lapply(t$b, function(b) .smt_eval(tr, b[[2L]]))
      body <- .smt_eval(tr, t$body)
      .smt_xreplace(body, vapply(t$b, `[[`, "", 1L), vals)
    },
    forall = {
      .smt_eval(tr, t$body)
      stop(.sympy_condition("NotImplementedError", "Quantifiers (forall) are not supported since SymPy has no representation for them"))
    },
    exists = {
      .smt_eval(tr, t$body)
      stop(.sympy_condition("NotImplementedError", "Quantifiers (exists) are not supported since SymPy has no representation for them"))
    },
    match = {
      .smt_eval(tr, t$e)
      for (cs in t$cases) .smt_eval(tr, cs)
      stop(.sympy_condition("NotImplementedError", "Match expressions are not fully supported yet."))
    },
    annotate = .smt_eval(tr, t$body))
}

.smt_apply <- function(tr, opn, ps) {
  if (is_sym(opn))
    stop(.sympy_condition("UnknownSMTLibOperatorError", sprintf("Cannot apply arguments to expression %s", sstr(opn))))
  op <- opn
  if (op %in% names(tr$macros)) {
    m <- tr$macros[[op]]
    n <- min(length(m$vars), length(ps))
    return(.smt_xreplace(m$body, m$vars[seq_len(n)], ps[seq_len(n)]))
  }
  if (op %in% tr$functions) return(fcall(op, ps))
  p1 <- function() ps[[1L]]
  switch(op,
    not = return(Not(p1())),
    and = return(do.call(And, ps)),
    or = return(do.call(Or, ps)),
    "=>" = {
      r <- ps[[length(ps)]]
      for (p in rev(ps[-length(ps)])) r <- Implies(p, r)
      return(r)
    },
    xor = return(do.call(Xor, ps)),
    "=" = return(.smt_chain(Eq, ps)),
    distinct = {
      if (length(ps) == 2L) return(Ne(ps[[1L]], ps[[2L]]))
      out <- list()
      for (i in seq_along(ps)) for (j in seq_along(ps)) if (j > i) out[[length(out) + 1L]] <- Ne(ps[[i]], ps[[j]])
      return(do.call(And, out))
    },
    "<" = return(.smt_chain(Lt, ps)),
    "<=" = return(.smt_chain(Le, ps)),
    ">" = return(.smt_chain(Gt, ps)),
    ">=" = return(.smt_chain(Ge, ps)),
    "+" = return(.add(ps)),
    "-" = return(if (length(ps) == 1L) .neg(p1()) else .sub(p1(), .add(ps[-1L]))),
    "*" = return(.mul(ps)),
    "/" = return(if (length(ps) == 2L) .div(ps[[1L]], ps[[2L]]) else .div(p1(), .mul(ps[-1L]))),
    pow = , "^" = return(.pow(ps[[1L]], ps[[2L]])),
    div = {
      r <- p1()
      for (p in ps[-1L]) r <- .div(.sub(r, Mod(r, Abs(p))), p)
      return(r)
    },
    mod = return(Mod(p1(), Abs(ps[[2L]]))),
    ite = return(Piecewise(list(ps[[2L]], ps[[1L]]), list(ps[[3L]], TRUE))),
    abs = return(Abs(p1())),
    to_real = return(p1()),
    to_int = return(floor(p1())),
    is_int = return(Eq(p1(), floor(p1()))))
  if (op %in% .SMT_UNSUPPORTED)
    stop(.sympy_condition("NotImplementedError", sprintf("The SMT-LIB operator '%s' has no SymPy equivalent and is not supported", op)))
  stop(.sympy_condition("UnknownSMTLibOperatorError", sprintf("Unknown operator or function: %s", op)))
}

.smt_add_assertion <- function(tr, a) tr$assertions[[length(tr$assertions) + 1L]] <- a

.smt_declare_const <- function(tr, name, sort) {
  sym <- Symbol(name)
  sort <- .smt_resolve_sort(tr, sort)
  if (identical(sort, "Real")) .smt_add_assertion(tr, Q$real(sym))
  else if (identical(sort, "Int")) .smt_add_assertion(tr, Q$integer(sym))
  tr$symbols[name] <- list(sym)
}

.smt_run <- function(tr, c) {
  switch(c$cmd,
    assert = .smt_add_assertion(tr, .smt_eval(tr, c$term)),
    "declare-const" = .smt_declare_const(tr, c$name, c$sort),
    "declare-fun" = {
      if (!length(c$domain)) .smt_declare_const(tr, c$name, c$sort)
      else tr$functions <- union(tr$functions, c$name)
    },
    "declare-sort" = { tr$sorts[c$name] <- list(c$arity) },
    "declare-datatype" = stop(.sympy_condition("NotImplementedError", "Algebraic datatypes (declare-datatype) are not supported")),
    "declare-datatypes" = stop(.sympy_condition("NotImplementedError", "Algebraic datatypes (declare-datatypes) are not supported")),
    "define-fun" = {
      body <- .smt_eval(tr, c$term)
      if (!length(c$vars)) tr$symbols[c$name] <- list(body)
      else tr$macros[c$name] <- list(list(vars = vapply(c$vars, `[[`, "", 1L), body = body))
    },
    "define-fun-rec" = {
      .smt_eval(tr, c$term)
      stop(.sympy_condition("NotImplementedError", "Recursive macros are not supported in SymPy because AST generation evaluates bottom-up."))
    },
    "define-const" = { tr$symbols[c$name] <- list(.smt_eval(tr, c$term)) },
    "define-sort" = {
      if (length(c$params))
        stop(.sympy_condition("NotImplementedError", "Parametric sort definitions (define-sort with parameters) are not supported"))
      tr$sort_aliases[c$name] <- list(.smt_resolve_sort(tr, c$sort))
    },
    generic = stop(.sympy_condition("UnknownSMTLibCommandError", sprintf("Unknown SMT-LIB command: %s", c$name))),
    pop = {
      n <- as.numeric(c$n$value)
      for (i in seq_len(max(0, n))) {
        if (!length(tr$stack)) break
        top <- tr$stack[[length(tr$stack)]]
        tr$stack[[length(tr$stack)]] <- NULL
        for (f in names(top)) assign(f, top[[f]], envir = tr)
      }
    },
    push = {
      n <- as.numeric(c$n$value)
      for (i in seq_len(max(0, n)))
        tr$stack[[length(tr$stack) + 1L]] <- list(assertions = tr$assertions, symbols = tr$symbols,
          macros = tr$macros, functions = tr$functions, sorts = tr$sorts, sort_aliases = tr$sort_aliases)
    },
    reset = .smt_state(tr),
    "reset-assertions" = { tr$assertions <- list() },
    NULL)
  invisible()
}

#' Parse SMT-LIB
#'
#' `parse_smtlib(s)` parses the text of an SMT-LIB script, as
#' `sympy.parsing.smtlib.parse_smtlib()`, and returns the declared symbols
#' and the assertions. Declaring a constant of sort `Int` or `Real` adds the
#' assumption `Q$integer()` or `Q$real()` to the assertions. Terms built
#' from the core, integer and real theories become expressions (`ite` gives
#' a `Piecewise`, `div` and `mod` are Euclidean), `define-fun` macros are
#' expanded and `push`, `pop`, `reset` and `reset-assertions` act on the
#' state. Invalid syntax signals an `SMTLibSyntaxError`, an unknown command
#' an `UnknownSMTLibCommandError`, an unknown operator an
#' `UnknownSMTLibOperatorError`, and quantifiers, datatypes, recursive or
#' parametric definitions and the bitvector, string and array operators a
#' `NotImplementedError`.
#'
#' @param s A character string with the SMT-LIB text.
#' @return A list with `symbols`, a named list of the symbols by name, and
#'   `assertions`, a list of expressions and predicates.
#' @examples
#' r <- parse_smtlib("(declare-const x Int) (declare-const y Int) (assert (= (+ x y) 10))")
#' r$symbols$x
#' r$assertions
#' @export
parse_smtlib <- function(s) {
  cmds <- .smt_parse(paste(s, collapse = "\n"))
  tr <- .smt_state()
  for (c in cmds) .smt_run(tr, c)
  list(symbols = tr$symbols, assertions = tr$assertions)
}

#' @export
format.sym_predicate <- function(x, ...) paste0("Q.", x$fact, "(", .str(x$expr), ")")
