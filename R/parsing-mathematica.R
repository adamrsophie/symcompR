# Wolfram Mathematica parser, after sympy/parsing/mathematica.py
# (MathematicaParser): named characters are resolved, the input is tokenized,
# the tokens are arranged into a FullForm list by the operator precedence
# table, and the FullForm list is turned into an expression. A FullForm list
# is an R list whose first element is the head; atoms are character strings.

.mma_syntax_error <- function(message) .pe_err("SyntaxError", message)
.mma_index_error <- function() .pe_err("IndexError", "list index out of range")

#' Parse a Mathematica expression
#'
#' `parse_mathematica()` translates a string of Wolfram Mathematica code into
#' an expression, as SymPy's `parse_mathematica` does: standard input form
#' and FullForm are both read, known heads such as `Sin`, `Log`, `Plus` or
#' `Power` become their SymPy counterparts, and any other head becomes an
#' undefined function. `parse_mathematica_to_fullformlist()` stops at the
#' FullForm, returned as nested lists of strings.
#'
#' @param s A character string of Mathematica code.
#' @return For `parse_mathematica()` an expression; for
#'   `parse_mathematica_to_fullformlist()` a string or a nested list.
#' @examples
#' parse_mathematica("Sin[x]^2 Tan[y]")
#' parse_mathematica("Log[2, 8] + 1/2")
#' parse_mathematica_to_fullformlist("x*(a + b)")
#' @export
parse_mathematica <- function(s) .mma_from_fullformlist_to_sympy(parse_mathematica_to_fullformlist(s))

#' @rdname parse_mathematica
#' @export
parse_mathematica_to_fullformlist <- function(s) {
  s2 <- .mma_named_characters_to_unicode(enc2utf8(s))
  s3 <- .mma_from_mathematica_to_tokens(s2)
  .mma_from_tokens_to_fullformlist(s3)
}

.mma_named_characters_to_unicode <- function(s) {
  if (grepl("\\[", s, fixed = TRUE)) {
    m <- gregexpr("\\\\\\[([a-zA-Z]+)\\]", s, perl = TRUE)[[1L]]
    if (m[1L] != -1L) {
      cs <- attr(m, "capture.start")[, 1L]; cl <- attr(m, "capture.length")[, 1L]
      names_ <- substring(s, cs, cs + cl - 1L)
      for (nm in names_) if (is.na(.MMA_NAMED_CHARACTERS[nm]))
        stop(.pe_value_error(sprintf("Unknown Mathematica named character: %s", nm)))
      regmatches(s, list(m)) <- list(unname(.MMA_NAMED_CHARACTERS[names_]))
    }
  }
  if (grepl("\\[", s, fixed = TRUE)) stop(.mma_syntax_error(sprintf("Unmatched '\\[' in '%s'", s)))
  s
}

# ---------------------------------------------------------------------------
# The operator precedence table, loosest first. Each level is
# list(type, strategy, ops) where ops maps a token to a head or a function.

.MMA_NUMBER <- "(?:[0-9]+(?:\\.[0-9]*)?|\\.[0-9]+)"
.MMA_OUT <- "%[0-9]+|%+"
.MMA_UNSET <- "=\\.(?![0-9])"
.MMA_ENCLOSURE_OPEN <- c("(", "[", "[[", "{", "\u2329")
.MMA_ENCLOSURE_CLOSE <- c(")", "]", "]]", "}", "\u232a")
.MMA_WHITESPACE_SENSITIVE <- c("_", "_.", "__", "___", "#", "##")
.MMA_BLANKS <- c("_", "_.", "__", "___")

.mma_get_neg <- function(x) {
  if (is.character(x) && grepl(paste0("^", .MMA_NUMBER), x, perl = TRUE)) return(paste0("-", x))
  list("Times", "-1", x)
}
.mma_get_inv <- function(x) list("Power", x, "-1")
.mma_get_derivative <- function(x) {
  if (is.list(x) && length(x) == 2L && is.list(x[[1L]]) && length(x[[1L]]) == 2L && identical(x[[1L]][[1L]], "Derivative"))
    return(list(list("Derivative", as.character(as.integer(x[[1L]][[2L]]) + 1L)), x[[2L]]))
  list(list("Derivative", "1"), x)
}

.mma_precedence <- function() {
  if (!is.null(.lark_env$mma_prec)) return(.lark_env$mma_prec)
  P <- function(type, strat, ops) list(type = type, strat = strat, ops = ops)
  prec <- list(
    P("Postfix", NA, list(";" = function(x) if (is.list(x) && length(x) && identical(x[[1L]], "CompoundExpression")) c(x, list("Null")) else list("CompoundExpression", x, "Null"))),
    P("Infix", "Flat", list(";" = "CompoundExpression")),
    P("Infix", "Right", list("=" = "Set", ":=" = "SetDelayed")),
    P("Postfix", NA, list("=." = "Unset")),
    P("Infix", "Right", list("+=" = "AddTo", "-=" = "SubtractFrom", "*=" = "TimesBy", "/=" = "DivideBy")),
    P("Infix", "Right", list("\u2234" = "Therefore")),
    P("Infix", "Left", list("//" = function(x, y) list(y, x))),
    P("Postfix", NA, list("&" = "Function")),
    P("Infix", "Left", list("/." = "ReplaceAll")),
    P("Infix", "Right", list("->" = "Rule", ":>" = "RuleDelayed", "\uf522" = "Rule", "\uf51f" = "RuleDelayed")),
    P("Infix", "Flat", list("\u2194" = "LeftRightArrow")),
    P("Infix", "Left", list("/;" = "Condition")),
    P("Infix", "Left", list(":" = function(x, y) if (.mma_is_symbol(x)) list("Pattern", x, y) else list("Optional", x, y))),
    P("Infix", "Flat", list("|" = "Alternatives")),
    P("Postfix", NA, list(".." = "Repeated", "..." = "RepeatedNull")),
    P("Infix", "Flat", list("||" = "Or")),
    P("Infix", "Flat", list("&&" = "And")),
    P("Prefix", NA, list("!" = "Not")),
    P("Infix", "Flat", list("===" = "SameQ", "=!=" = "UnsameQ")),
    P("Infix", "Flat", list("==" = "Equal", "!=" = "Unequal", "<=" = "LessEqual", "<" = "Less", ">=" = "GreaterEqual", ">" = "Greater")),
    P("Infix", "Flat", list("\u2248" = "TildeTilde")),
    P("Infix", "Flat", list(";;" = "Span")),
    P("Infix", "Flat", list("+" = "Plus", "-" = "Plus")),
    P("Infix", "Flat", list("\u2295" = "CirclePlus")),
    P("Infix", "Flat", list("\u22c6" = "Star")),
    P("Infix", "Flat", list("*" = "Times", "/" = "Times")),
    P("Infix", "Flat", list("\u2297" = "CircleTimes")),
    P("Infix", "Flat", list("." = "Dot")),
    P("Prefix", NA, list("-" = function(x) .mma_get_neg(x), "+" = function(x) list("Plus", x))),
    P("Prefix", NA, list("\u2207" = "Del", "\u25a1" = "Square")),
    P("Infix", "Right", list("^" = "Power")),
    P("Postfix", NA, list("'" = function(x) .mma_get_derivative(x))),
    P("Infix", "Right", list("@@" = "Apply", "/@" = "Map", "//@" = "MapAll",
                             "@@@" = function(x, y) list("Apply", x, y, list("List", "1")))),
    P("Infix", "Right", list("@" = function(x, y) list(x, y))),
    P("Postfix", NA, list("!" = "Factorial", "!!" = "Factorial2", "--" = "Decrement")),
    P("Infix", NA, list("[" = function(x, y) c(list(x), y), "[[" = function(x, y) c(list("Part", x), y))),
    P("Prefix", NA, list("{" = function(x) c(list("List"), x), "(" = function(x) x[[1L]],
                         "\u2329" = function(x) c(list("AngleBracket"), x))),
    P("Infix", NA, list("?" = "PatternTest")),
    P("Postfix", NA, list("_" = function(x) list("Pattern", x, list("Blank")),
                          "_." = function(x) list("Optional", list("Pattern", x, list("Blank"))),
                          "__" = function(x) list("Pattern", x, list("BlankSequence")),
                          "___" = function(x) list("Pattern", x, list("BlankNullSequence")))),
    P("Infix", NA, list("_" = "Blank", "__" = "BlankSequence", "___" = "BlankNullSequence")),
    P("Prefix", NA, list("#" = "Slot", "##" = "SlotSequence")),
    P("Infix", "Flat", list("::" = function(...) {
      a <- list(...)
      c(list("MessageName", a[[1L]]), lapply(a[-1L], function(n) if (is.character(n)) list("_Str", n) else n))
    })))
  .lark_env$mma_prec <- prec
  prec
}

.mma_missing_default <- function(token) {
  switch(token, "#" = list("Slot", "1"), "##" = list("SlotSequence", "1"), "_" = list("Blank"),
         "_." = list("Optional", list("Blank")), "__" = list("BlankSequence"), "___" = list("BlankNullSequence"))
}
.MMA_MISSING_DEFAULT_KEYS <- c("#", "##", "_", "_.", "__", "___")
.MMA_ENCLOSURE_PREFIXES <- c("(", "{", "\u2329")
.mma_enclosure_build <- function(tok, group) {
  switch(tok, "(" = group[[1L]], "{" = c(list("List"), group), c(list("AngleBracket"), group))
}

.mma_tighter_than_derivative <- function() {
  if (!is.null(.lark_env$mma_tighter)) return(.lark_env$mma_tighter)
  prec <- .mma_precedence()
  prime <- which(vapply(prec, function(l) "'" %in% names(l$ops), TRUE))[1L]
  out <- character()
  for (l in prec[seq.int(prime + 1L, length(prec))])
    if (l$type == "Infix") out <- c(out, setdiff(names(l$ops), c("[", "[[")))
  .lark_env$mma_tighter <- out
  out
}

# ---------------------------------------------------------------------------
# Token predicates

.mma_is_str <- function(t) is.character(t) && length(t) == 1L
.mma_tok_is <- function(t, values) .mma_is_str(t) && t %in% values
# memoized: the identifier pattern is large, and the same tokens are tested
# many times
.mma_memo <- new.env(parent = emptyenv())
.mma_re_match <- function(re, s) .mma_re_cached("m", re, s, paste0("(*UTF)^(?:", re, ")"))
.mma_re_full <- function(re, s) .mma_re_cached("f", re, s, paste0("(*UTF)^(?:", re, ")$"))
.mma_re_cached <- function(kind, re, s, pat) {
  key <- paste(kind, if (identical(re, .MMA_LITERAL)) "L" else re, s, sep = "\r")
  v <- .mma_memo[[key]]
  if (is.null(v)) {
    v <- grepl(pat, s, perl = TRUE)
    assign(key, v, envir = .mma_memo)
  }
  v
}
.mma_is_op <- function(t) {
  if (is.list(t)) return(FALSE)
  if (.mma_re_match(.MMA_LITERAL, t)) return(FALSE)
  if (.mma_re_match(paste0("-?", .MMA_NUMBER), t)) return(FALSE)
  TRUE
}
.mma_is_symbol <- function(t) .mma_is_str(t) && .mma_re_full(.MMA_LITERAL, t)
.mma_is_valid_star1 <- function(t) if (.mma_tok_is(t, c(")", "}"))) TRUE else !.mma_is_op(t)
.mma_is_valid_star2 <- function(t) if (.mma_tok_is(t, c("(", "{"))) TRUE else !.mma_is_op(t)
.mma_in_dict <- function(t, ops) .mma_is_str(t) && t %in% names(ops)
.mma_eq <- function(a, b) if (.mma_is_str(a) && .mma_is_str(b)) a == b else if (is.list(a) && is.list(b)) identical(a, b) else FALSE
.mma_last <- function(x) {
  if (is.character(x)) { if (!nchar(x)) stop(.pe_err("IndexError", "string index out of range")); return(substring(x, nchar(x))) }
  if (!length(x)) stop(.mma_index_error())
  x[[length(x)]]
}
.mma_len <- function(x) if (is.character(x)) nchar(x) else length(x)

# ---------------------------------------------------------------------------
# Tokenizer

.mma_escape <- function(s) gsub("([.^$*+?()\\[\\]{}|\\\\-])", "\\\\\\1", s, perl = TRUE)

.mma_tokenizer <- function() {
  if (!is.null(.lark_env$mma_tok)) return(.lark_env$mma_tok)
  esc <- c(.MMA_ENCLOSURE_OPEN, .MMA_ENCLOSURE_CLOSE)
  for (l in .mma_precedence()) for (k in names(l$ops)) if (k != "=.") esc <- c(esc, k)
  esc <- esc[order(-nchar(esc))]
  pat <- paste0("(*UTF)(", paste(c(.MMA_LITERAL, .MMA_NUMBER, .MMA_OUT, .MMA_UNSET, vapply(esc, .mma_escape, ""), ",", "\n"),
                           collapse = "|"), ")")
  .lark_env$mma_tok <- pat
  pat
}

.mma_from_mathematica_to_tokens <- function(code) {
  tokenizer <- .mma_tokenizer()
  splits <- list()
  repeat {
    st <- regexpr("\"", code, fixed = TRUE)
    if (st == -1L) {
      if (nchar(code) > 0L) splits[[length(splits) + 1L]] <- code
      break
    }
    rest <- substring(code, st + 1L)
    me <- regexpr("(?<!\\\\)\"", rest, perl = TRUE)
    if (me == -1L) stop(.mma_syntax_error("mismatch in string \"  \" expression"))
    en <- st + me
    if (st > 1L) splits[[length(splits) + 1L]] <- substr(code, 1L, st - 1L)
    splits[[length(splits) + 1L]] <- list("_Str", gsub("\\\"", "\"", substr(code, st + 1L, en - 1L), fixed = TRUE))
    code <- substring(code, en + 1L)
  }
  for (i in seq_along(splits)) {
    cs <- splits[[i]]
    if (is.list(cs)) next
    repeat {
      ps <- regexpr("(*", cs, fixed = TRUE)
      if (ps == -1L) break
      pe <- regexpr("*)", cs, fixed = TRUE)
      if (pe == -1L || pe < ps) stop(.mma_syntax_error("mismatch in comment (*  *) code"))
      cs <- paste0(substr(cs, 1L, ps - 1L), substring(cs, pe + 2L))
    }
    splits[[i]] <- cs
  }
  token_split <- function(code) {
    if (is.list(code)) return(list(code))
    m <- gregexpr(tokenizer, code, perl = TRUE)[[1L]]
    has <- m[1L] != -1L
    if (has || !grepl("[^\\x01-\\x7f]", code, perl = TRUE)) {
      out <- list()
      if (!has) return(out)
      starts <- as.integer(m); lens <- attr(m, "match.length")
      toks <- substring(code, starts, starts + lens - 1L)
      for (i in seq_along(toks)) {
        tok <- toks[i]
        if (i > 1L && starts[i - 1L] + lens[i - 1L] != starts[i]) {
          prev <- toks[i - 1L]
          if ((tok %in% .MMA_WHITESPACE_SENSITIVE && !.mma_is_op(prev)) ||
              (prev %in% .MMA_WHITESPACE_SENSITIVE && !.mma_is_op(tok))) out[[length(out) + 1L]] <- "*"
        }
        out[[length(out) + 1L]] <- tok
      }
      return(out)
    }
    list(code)
  }
  tokens <- unlist(lapply(splits, token_split), recursive = FALSE)
  if (is.null(tokens)) tokens <- list()
  while (length(tokens) && .mma_tok_is(tokens[[1L]], "\n")) tokens <- tokens[-1L]
  while (length(tokens) && .mma_tok_is(tokens[[length(tokens)]], "\n")) tokens <- tokens[-length(tokens)]
  for (i in seq_along(tokens)) {
    tk <- tokens[[i]]
    if (.mma_is_str(tk) && (startsWith(tk, "System`") || startsWith(tk, "Global`")) && .mma_is_symbol(tk)) {
      name <- sub("^[^`]*`", "", tk)
      if (nchar(name)) tokens[[i]] <- name
    }
  }
  for (i in seq_along(tokens)) {
    tk <- tokens[[i]]
    if (.mma_is_str(tk) && startsWith(tk, "%") && .mma_re_full(.MMA_OUT, tk)) {
      rest <- substring(tk, 2L)
      tokens[[i]] <- if (nchar(rest) && grepl("^[0-9]+$", rest)) list("Out", rest)
                     else if (nchar(tk) == 1L) list("Out") else list("Out", as.character(-nchar(tk)))
    }
  }
  tokens
}

# ---------------------------------------------------------------------------
# Tokens to FullForm. Python mutates the token lists in place; here the
# helpers return the updated lists.

.mma_from_tokens_to_fullformlist <- function(tokens) {
  stack <- list(list())
  open_seq <- character()
  pointer <- 1L
  top <- function() stack[[length(stack)]]
  while (pointer <= length(tokens)) {
    token <- tokens[[pointer]]
    if (.mma_tok_is(token, .MMA_ENCLOSURE_OPEN)) {
      stack[[length(stack)]] <- c(top(), list(token))
      open_seq <- c(open_seq, token)
      stack[[length(stack) + 1L]] <- list()
    } else if (.mma_tok_is(token, ",")) {
      if (.mma_len(top()) == 0L) {
        if (length(stack) < 2L || !length(open_seq)) stop(.mma_index_error())
        if (.mma_eq(.mma_last(stack[[length(stack) - 1L]]), open_seq[length(open_seq)]))
          stop(.mma_syntax_error(sprintf("%s cannot be followed by comma ,", open_seq[length(open_seq)])))
      }
      stack[[length(stack)]] <- .mma_parse_after_braces(top())
      stack[[length(stack) + 1L]] <- list()
    } else if (.mma_tok_is(token, .MMA_ENCLOSURE_CLOSE)) {
      ind <- match(token, .MMA_ENCLOSURE_CLOSE)
      if (!length(open_seq)) stop(.mma_index_error())
      if (.MMA_ENCLOSURE_OPEN[ind] != open_seq[length(open_seq)]) {
        if (token == "]]" && open_seq[length(open_seq)] == "[") {
          if (length(open_seq) < 2L) stop(.mma_index_error())
          if (open_seq[length(open_seq) - 1L] == "[") {
            tokens <- append(tokens, list("]"), after = pointer)
          } else if (open_seq[length(open_seq) - 1L] == "[[") {
            if (pointer + 1L > length(tokens)) stop(.mma_index_error())
            nx <- tokens[[pointer + 1L]]
            if (.mma_tok_is(nx, "]")) tokens[[pointer + 1L]] <- "]]"
            else if (.mma_tok_is(nx, "]]")) {
              tokens[[pointer + 1L]] <- "]]"
              tokens <- append(tokens, list("]"), after = pointer + 1L)
            } else stop(.mma_syntax_error("unmatched enclosure"))
          }
        } else stop(.mma_syntax_error("unmatched enclosure"))
      }
      if (.mma_len(top()) == 0L && length(stack) < 2L) stop(.mma_index_error())
      if (.mma_len(top()) == 0L && .mma_eq(.mma_last(stack[[length(stack) - 1L]]), "("))
        stop(.mma_syntax_error("( ) not valid syntax"))
      stack[[length(stack)]] <- .mma_parse_after_braces(top(), TRUE)
      new_el <- list()
      while (!.mma_eq(.mma_last(top()), open_seq[length(open_seq)])) {
        new_el[[length(new_el) + 1L]] <- top()
        stack[[length(stack)]] <- NULL
      }
      new_el <- rev(new_el)
      if (open_seq[length(open_seq)] == "(" && length(new_el) != 1L)
        stop(.mma_syntax_error(sprintf("( must be followed by one expression, %i detected", length(new_el))))
      stack[[length(stack)]] <- c(top(), list(new_el))
      open_seq <- open_seq[-length(open_seq)]
    } else {
      stack[[length(stack)]] <- c(top(), list(token))
    }
    pointer <- pointer + 1L
  }
  if (length(stack) != 1L) stop(.pe_err("RuntimeError", "Stack should have only one element"))
  .mma_parse_after_braces(stack[[1L]])
}

# returns list(lines, tokens)
.mma_remove_newlines <- function(lines, tokens, inside) {
  pointer <- 1L
  size <- length(tokens)
  while (pointer <= size) {
    token <- tokens[[pointer]]
    if (.mma_tok_is(token, "\n")) {
      if (inside) { tokens <- tokens[-pointer]; size <- size - 1L; next }
      if (pointer == 1L) { tokens <- tokens[-1L]; size <- size - 1L; next }
      if (pointer > 2L) {
        prev <- tryCatch(.mma_parse_after_braces(tokens[seq_len(pointer - 1L)], inside),
                         SyntaxError = function(e) NULL)
        if (is.null(prev)) { tokens <- tokens[-pointer]; size <- size - 1L; next }
      } else prev <- tokens[[1L]]
      if (.mma_len(prev) > 0L && .mma_eq(if (is.list(prev)) prev[[1L]] else substr(prev, 1L, 1L), "CompoundExpression"))
        lines <- c(lines, prev[-1L])
      else lines <- c(lines, list(prev))
      tokens <- tokens[-seq_len(pointer - 1L)]
      size <- size - (pointer - 1L)
      pointer <- 1L
      next
    }
    pointer <- pointer + 1L
  }
  list(lines = lines, tokens = tokens)
}

.mma_add_missing_asterisks <- function(tokens) {
  size <- length(tokens)
  pointer <- 1L
  while (pointer <= size) {
    if (pointer > 1L && .mma_is_valid_star1(tokens[[pointer - 1L]]) && .mma_is_valid_star2(tokens[[pointer]])) {
      if (.mma_tok_is(tokens[[pointer]], "(")) {
        tokens[[pointer]] <- "*"
        tokens[[pointer + 1L]] <- tokens[[pointer + 1L]][[1L]]
      } else {
        tokens <- append(tokens, list("*"), after = pointer - 1L)
        pointer <- pointer + 1L
        size <- size + 1L
      }
    }
    pointer <- pointer + 1L
  }
  tokens
}

# _complete_left_operand: returns list(tokens, removed)
.mma_complete_left_operand <- function(tokens, pointer) {
  removed <- 0L
  while (pointer > 1L) {
    if (.mma_tok_is(tokens[[pointer - 1L]], "'")) {
      start <- pointer - 1L
      while (start > 1L && .mma_tok_is(tokens[[start - 1L]], "'")) start <- start - 1L
      r <- .mma_complete_left_operand(tokens, start)
      tokens <- r$tokens
      if (r$removed) { removed <- removed + r$removed; pointer <- pointer - r$removed; next }
      if (start == 1L || .mma_is_op(tokens[[start - 1L]])) return(list(tokens = tokens, removed = removed))
      if (start > 2L && .mma_tok_is(tokens[[start - 2L]], .mma_tighter_than_derivative()))
        return(list(tokens = tokens, removed = removed))
      head <- tokens[[start - 1L]]
      for (k in seq_len(pointer - start)) head <- .mma_get_derivative(head)
      tokens <- c(tokens[seq_len(start - 2L)], list(head), tokens[seq.int(pointer, length.out = length(tokens) - pointer + 1L)])
      removed <- removed + (pointer - start)
      pointer <- start
    } else if (pointer > 2L && .mma_tok_is(tokens[[pointer - 2L]], c("[", "[["))) {
      r <- .mma_complete_left_operand(tokens, pointer - 2L)
      tokens <- r$tokens
      if (r$removed) { removed <- removed + r$removed; pointer <- pointer - r$removed; next }
      if (pointer < 4L || .mma_is_op(tokens[[pointer - 3L]])) return(list(tokens = tokens, removed = removed))
      head <- tokens[[pointer - 3L]]; args <- tokens[[pointer - 1L]]
      node <- if (.mma_tok_is(tokens[[pointer - 2L]], "[")) c(list(head), args) else c(list("Part", head), args)
      tokens <- c(tokens[seq_len(pointer - 4L)], list(node), tokens[seq.int(pointer, length.out = length(tokens) - pointer + 1L)])
      removed <- removed + 2L
      pointer <- pointer - 2L
    } else if (pointer > 2L && .mma_tok_is(tokens[[pointer - 2L]], .MMA_ENCLOSURE_PREFIXES)) {
      node <- .mma_enclosure_build(tokens[[pointer - 2L]], tokens[[pointer - 1L]])
      tokens <- c(tokens[seq_len(pointer - 3L)], list(node), tokens[seq.int(pointer, length.out = length(tokens) - pointer + 1L)])
      removed <- removed + 1L
      pointer <- pointer - 1L
    } else return(list(tokens = tokens, removed = removed))
  }
  list(tokens = tokens, removed = removed)
}

.mma_check_op_compatible <- function(op1, op2) {
  if (.mma_eq(op1, op2)) return(TRUE)
  if (!.mma_is_str(op1)) stop(.pe_type_error("unhashable type: 'list'"))
  if (op1 %in% c("*", "/") && op2 %in% c("*", "/")) return(TRUE)
  if (op1 %in% c("+", "-") && op2 %in% c("+", "-")) return(TRUE)
  FALSE
}

.mma_call_op <- function(op, node) if (is.function(op)) do.call(op, node) else node

.mma_parse_after_braces <- function(tokens, inside = FALSE) {
  changed <- FALSE
  rn <- .mma_remove_newlines(list(), tokens, inside)
  lines <- rn$lines; tokens <- rn$tokens
  prec <- .mma_precedence()
  for (lev in rev(prec)) {
    op_type <- lev$type; strat <- lev$strat; op_dict <- lev$ops
    if ("*" %in% names(op_dict)) tokens <- .mma_add_missing_asterisks(tokens)
    size <- length(tokens)
    pointer <- 1L
    while (pointer <= size) {
      token <- tokens[[pointer]]
      if (.mma_in_dict(token, op_dict)) {
        op_name <- op_dict[[token]]
        named <- is.character(op_name)
        node <- if (named) list(op_name) else list()
        first_index <- if (named) 2L else 1L
        if (token %in% c("+", "-") && op_type == "Prefix" && pointer > 1L && !.mma_is_op(tokens[[pointer - 1L]])) {
          pointer <- pointer + 1L
          next
        }
        if (token == "+" && op_type == "Prefix" &&
            ((pointer > 2L && .mma_tok_is(tokens[[pointer - 1L]], "+") && !.mma_is_op(tokens[[pointer - 2L]])) ||
             (pointer == 1L && size > 2L && .mma_tok_is(tokens[[pointer + 2L]], c("+", "-"))))) {
          tokens <- tokens[-pointer]
          size <- size - 1L
          changed <- TRUE
          next
        }
        if (op_type == "Infix" && token %in% .MMA_BLANKS) {
          if (pointer == size || .mma_is_op(tokens[[pointer + 1L]])) { pointer <- pointer + 1L; next }
          blank <- list(op_name, tokens[[pointer + 1L]])
          if (pointer > 1L && .mma_is_symbol(tokens[[pointer - 1L]])) {
            tokens <- c(tokens[seq_len(pointer - 2L)], list(list("Pattern", tokens[[pointer - 1L]], blank)),
                        tokens[seq.int(pointer + 2L, length.out = size - pointer - 1L)])
            size <- size - 2L
          } else {
            tokens <- c(tokens[seq_len(pointer - 1L)], list(blank), tokens[seq.int(pointer + 2L, length.out = size - pointer - 1L)])
            size <- size - 1L
            pointer <- pointer + 1L
          }
          changed <- TRUE
          next
        }
        if (op_type %in% c("Infix", "Postfix")) {
          r <- .mma_complete_left_operand(tokens, pointer)
          tokens <- r$tokens
          if (r$removed) {
            pointer <- pointer - r$removed
            size <- size - r$removed
            changed <- TRUE
          }
        }
        if (op_type == "Infix") {
          if (pointer == 1L || pointer == size || .mma_is_op(tokens[[pointer - 1L]]) || .mma_is_op(tokens[[pointer + 1L]])) {
            pointer <- pointer + 1L
            next
          }
        }
        if (!(token %in% .MMA_MISSING_DEFAULT_KEYS)) {
          if (op_type == "Postfix" && (pointer == 1L || .mma_is_op(tokens[[pointer - 1L]]))) { pointer <- pointer + 1L; next }
          if (op_type == "Prefix" && (pointer == size || .mma_is_op(tokens[[pointer + 1L]]))) { pointer <- pointer + 1L; next }
        }
        changed <- TRUE
        placed <- TRUE
        if (op_type == "Infix") {
          arg1 <- tokens[[pointer - 1L]]
          arg2 <- tokens[[pointer + 1L]]
          tokens <- tokens[-c(pointer - 1L, pointer + 1L)]
          if (token == "/") arg2 <- .mma_get_inv(arg2) else if (token == "-") arg2 <- .mma_get_neg(arg2)
          pointer <- pointer - 1L
          size <- size - 2L
          node <- c(node, list(arg1))
          if (identical(strat, "Flat")) {
            while (pointer + 2L <= size && .mma_check_op_compatible(tokens[[pointer + 1L]], token)) {
              node <- c(node, list(arg2))
              other_op <- tokens[[pointer + 1L]]
              arg2 <- tokens[[pointer + 2L]]
              tokens <- tokens[-c(pointer + 1L, pointer + 2L)]
              if (.mma_tok_is(other_op, "/")) arg2 <- .mma_get_inv(arg2)
              else if (.mma_tok_is(other_op, "-")) arg2 <- .mma_get_neg(arg2)
              size <- size - 2L
            }
            node <- c(node, list(arg2))
          } else if (identical(strat, "Right")) {
            operands <- list(arg1, arg2)
            chain <- token
            while (pointer + 2L <= size && .mma_in_dict(tokens[[pointer + 1L]], op_dict)) {
              chain <- c(chain, tokens[[pointer + 1L]])
              operands[[length(operands) + 1L]] <- tokens[[pointer + 2L]]
              tokens <- tokens[-c(pointer + 1L, pointer + 2L)]
              size <- size - 2L
            }
            rest <- operands[[length(operands)]]
            if (length(chain) > 1L) for (i in seq.int(length(chain), 2L)) {
              link <- op_dict[[chain[i]]]
              rest <- if (is.character(link)) list(link, operands[[i]], rest) else link(operands[[i]], rest)
            }
            node <- if (named) list(op_name, operands[[1L]], rest) else list(operands[[1L]], rest)
          } else if (identical(strat, "Left")) {
            while (pointer + 1L <= size && .mma_tok_is(tokens[[pointer + 1L]], token)) {
              node[[first_index]] <- if (named) list(op_name, node[[first_index]], arg2) else op_name(node[[first_index]], arg2)
              arg2 <- tokens[[pointer + 2L]]
              tokens <- tokens[-c(pointer + 1L, pointer + 2L)]
              size <- size - 2L
            }
            node <- c(node, list(arg2))
          } else {
            node <- c(node, list(arg2))
          }
        } else if (op_type == "Prefix") {
          if (token %in% c("#", "##")) {
            nxt <- if (pointer + 1L <= size) tokens[[pointer + 1L]] else NULL
            is_number <- .mma_is_str(nxt) && .mma_re_full(.MMA_NUMBER, nxt)
            if (is.null(nxt) || .mma_is_op(nxt) || (token == "##" && !is_number)) {
              tokens[[pointer]] <- .mma_missing_default(token)
              placed <- FALSE
              node <- NULL
            } else if (is_number) {
              node <- c(node, list(nxt)); tokens <- tokens[-(pointer + 1L)]; size <- size - 1L
            } else {
              node <- c(node, list(list("_Str", nxt))); tokens <- tokens[-(pointer + 1L)]; size <- size - 1L
            }
          } else if (pointer == size || .mma_is_op(tokens[[pointer + 1L]])) {
            tokens[[pointer]] <- .mma_missing_default(token)
            placed <- FALSE
            node <- NULL
          } else {
            node <- c(node, list(tokens[[pointer + 1L]]))
            tokens <- tokens[-(pointer + 1L)]
            size <- size - 1L
          }
        } else {
          if (pointer == 1L || .mma_is_op(tokens[[pointer - 1L]]) ||
              (token %in% .MMA_BLANKS && !.mma_is_symbol(tokens[[pointer - 1L]]))) {
            tokens[[pointer]] <- .mma_missing_default(token)
            pointer <- pointer + 1L
            next
          }
          node <- c(node, list(tokens[[pointer - 1L]]))
          tokens <- tokens[-(pointer - 1L)]
          pointer <- pointer - 1L
          size <- size - 1L
        }
        if (placed) {
          if (!named) {
            nn <- do.call(op_name, node)
            node <- nn
          }
          tokens[[pointer]] <- node
        } else if (!named && !is.null(node)) {
          tokens[[pointer]] <- do.call(op_name, node)
        }
      }
      pointer <- pointer + 1L
    }
  }
  if (length(tokens) > 1L || (!length(lines) && !length(tokens))) {
    if (changed) return(.mma_parse_after_braces(tokens, inside))
    stop(.mma_syntax_error("unable to create a single AST for the expression"))
  }
  if (length(lines)) {
    if (!length(tokens)) stop(.mma_index_error())
    if (.mma_len(tokens[[1L]]) && .mma_eq(if (is.list(tokens[[1L]])) tokens[[1L]][[1L]] else substr(tokens[[1L]], 1L, 1L), "CompoundExpression"))
      tokens <- tokens[[1L]][-1L]
    return(c(list("CompoundExpression"), lines, tokens))
  }
  tokens[[1L]]
}

# ---------------------------------------------------------------------------
# FullForm strings and lists to expressions

.mma_from_fullform_to_fullformlist <- function(wmexpr) {
  m <- gregexpr("[][,]", wmexpr, perl = TRUE)[[1L]]
  out <- list()
  path <- integer()
  get_at <- function(p) { x <- out; for (i in p) x <- x[[i]]; x }
  append_at <- function(p, v) {
    if (!length(p)) { out[[length(out) + 1L]] <<- v; return(length(out)) }
    x <- get_at(p)
    k <- length(x) + 1L
    out[[c(p, k)]] <<- v
    k
  }
  last <- 1L
  if (m[1L] != -1L) for (j in seq_along(m)) {
    position <- m[j]
    ch <- substr(wmexpr, position, position)
    last_expr <- trimws(gsub("[][,]", "", substr(wmexpr, last, position - 1L)))
    if (ch == ",") {
      if (last_expr != "") append_at(path, last_expr)
    } else if (ch == "]") {
      if (last_expr != "") append_at(path, last_expr)
      path <- path[-length(path)]
    } else {
      k <- append_at(path, list(last_expr))
      path <- c(path, k)
    }
    last <- position + 1L
  }
  out[[1L]]
}

.mma_converter_fullformsympy <- function(expr) {
  if (is.list(expr)) {
    if (!length(expr)) stop(.pe_value_error("Empty list of expressions"))
    return(fcall(expr[[1L]], lapply(expr[-1L], .mma_converter_fullformsympy)))
  }
  Symbol(expr)
}

.mma_chain3 <- function(s) {
  e <- .mma_converter_fullformsympy(.mma_from_fullform_to_fullformlist(s))
  for (head in names(.mma_node_conversions())) {
    conv <- .mma_node_conversions()[[head]]
    rep_ <- function(x) {
      if (!is_sym(x) || !length(x$args)) return(x)
      args <- lapply(x$args, rep_)
      if (x$type == "Function" && identical(x$name, head)) return(conv(args))
      sym_rebuild(x, args)
    }
    e <- rep_(e)
  }
  e
}

.mma_rational <- function(args) {
  for (a in args) if (!(a$type %in% c("Integer", "Rational", "Float")))
    stop(.pe_type_error(sprintf("invalid input: %s", sstr(a))))
  if (length(args) == 1L) return(args[[1L]])
  .div(args[[1L]], args[[2L]])
}

.mma_parse_function <- function(args) {
  if (length(args) == 1L) {
    arg <- args[[1L]]
    slots <- Filter(function(a) a$type == "Function" && identical(a$name, "Slot"), .mma_applied_functions(arg))
    if (!length(slots)) stop(.pe_value_error("max() iterable argument is empty"))
    nums <- lapply(slots, function(a) a$args[[1L]])
    best <- nums[[1L]]
    for (v in nums[-1L]) if (isTRUE(.pe_truthy(Gt(v, best)))) best <- v
    if (best$type == "Integer") {
      n <- as.integer(best$value)
      vars <- lapply(seq_len(n) - 1L, function(i) Dummy(paste0("dummy", i)))
      rule <- lapply(seq_len(n), function(i) list(fcall("Slot", list(.int(i))), vars[[i]]))
      return(.mma_lambda(vars, xreplace(arg, rule)))
    }
    return(.mma_lambda(list(), arg))
  }
  if (length(args) == 2L) {
    v <- args[[1L]]
    vars <- if (v$type == "Tuple") v$args else list(v)
    return(.mma_lambda(vars, args[[2L]]))
  }
  stop(.mma_syntax_error("Function node expects 1 or 2 arguments"))
}

# Lambda.__new__ returns the identity function for Lambda(x, x), which
# prints as Lambda(x, x)
.mma_lambda <- function(vars, body) {
  if (length(vars) == 1L && identical(vars[[1L]]$key, body$key)) return(Lambda(list(Symbol("x")), Symbol("x")))
  Lambda(vars, body)
}

.mma_applied_functions <- function(e) {
  out <- list(); seen <- character()
  walk <- function(x) {
    if (!is_sym(x)) return(invisible())
    if (x$type == "Function" && !(x$key %in% seen)) { seen <<- c(seen, x$key); out[[length(out) + 1L]] <<- x }
    for (a in x$args) walk(a)
  }
  walk(e)
  out
}

.mma_node_conversions <- function() {
  if (!is.null(.lark_env$mma_conv)) return(.lark_env$mma_conv)
  f1 <- function(name) { force(name); function(a) fcall(name, a) }
  conv <- list(
    Times = function(a) .mul(a), Plus = function(a) .add(a), Power = function(a) .pow(a[[1L]], a[[2L]]),
    Rational = .mma_rational,
    Log = function(a) fcall("log", rev(a)), Log2 = function(a) fcall("log", list(a[[1L]], .int(2L))),
    Log10 = function(a) fcall("log", list(a[[1L]], .int(10L))), Exp = f1("exp"),
    Sqrt = function(a) .pow(a[[1L]], S_Half),
    Sin = f1("sin"), Cos = f1("cos"), Tan = f1("tan"), Cot = f1("cot"), Sec = f1("sec"), Csc = f1("csc"),
    ArcSin = f1("asin"), ArcCos = f1("acos"),
    ArcTan = function(a) if (length(a) == 2L) fcall("atan2", rev(a)) else fcall("atan", a),
    ArcCot = f1("acot"), ArcSec = f1("asec"), ArcCsc = f1("acsc"),
    Sinh = f1("sinh"), Cosh = f1("cosh"), Tanh = f1("tanh"), Coth = f1("coth"), Sech = f1("sech"), Csch = f1("csch"),
    ArcSinh = f1("asinh"), ArcCosh = f1("acosh"), ArcTanh = f1("atanh"), ArcCoth = f1("acoth"),
    ArcSech = f1("asech"), ArcCsch = f1("acsch"),
    Expand = function(a) expand(a[[1L]]), Im = f1("im"), Re = f1("re"),
    Flatten = function(a) stop(.pe_err("NotImplementedError", "Flatten is not available")),
    Polylog = f1("polylog"), Cancel = function(a) cancel(a[[1L]]),
    TrigExpand = function(a) expand_trig(a[[1L]]), Sign = f1("sign"), Simplify = function(a) simplify(a[[1L]]),
    Defer = function(a) UnevaluatedExpr(a[[1L]]), Identity = function(a) a[[1L]],
    Null = function(a) S_Zero, Mod = function(a) fcall("Mod", a), Max = function(a) do.call(Max, a),
    Min = function(a) do.call(Min, a), Pochhammer = f1("RisingFactorial"), ExpIntegralEi = f1("Ei"),
    SinIntegral = f1("Si"), CosIntegral = f1("Ci"), AiryAi = f1("airyai"), AiryAiPrime = f1("airyaiprime"),
    AiryBi = f1("airybi"), AiryBiPrime = f1("airybiprime"), LogIntegral = f1("li"), PrimePi = f1("primepi"),
    Prime = function(a) prime(a[[1L]]), PrimeQ = function(a) isprime(a[[1L]]),
    List = function(a) do.call(Tuple, a), Greater = function(a) do.call(Gt, a), GreaterEqual = function(a) do.call(Ge, a),
    Less = function(a) do.call(Lt, a), LessEqual = function(a) do.call(Le, a), Equal = function(a) do.call(Eq, a),
    Or = function(a) do.call(Or, a), And = function(a) do.call(And, a), Not = function(a) Not(a[[1L]]),
    Function = .mma_parse_function, Factorial = f1("factorial"))
  .lark_env$mma_conv <- conv
  conv
}

.MMA_ATOMS <- c("I", "Pi", "ExponentialE", "ImaginaryI", "ImaginaryJ")

.mma_call_head <- function(head, args) {
  if (is.function(head)) return(head(args))
  if (is_sym(head) && head$type == "Lambda") return(do.call(apply_lambda, c(list(head), args)))
  stop(.pe_type_error(sprintf("'%s' object is not callable", .pe_class_name(head))))
}

.mma_derivative_to_sympy <- function(expr, recurse) {
  applied <- is.list(expr[[1L]]) && length(expr[[1L]]) == 2L && is.list(expr[[1L]][[1L]]) &&
    identical(expr[[1L]][[1L]][[1L]], "Derivative")
  bare <- is.list(expr[[1L]]) && identical(expr[[1L]][[1L]], "Derivative") && length(expr) == 2L
  if (!(applied || bare)) return(NULL)
  if (applied) { spec <- expr[[1L]][[1L]]; fname <- expr[[1L]][[2L]] } else { spec <- expr[[1L]]; fname <- expr[[2L]] }
  orders <- spec[-1L]
  if (!all(vapply(orders, function(o) .mma_is_str(o) && grepl("^[0-9]+$", o), TRUE))) return(NULL)
  orders <- as.integer(unlist(orders))
  conv <- .mma_node_conversions()
  func <- if (.mma_is_str(fname)) {
    if (!is.null(conv[[fname]])) conv[[fname]] else { fn <- fname; function(a) fcall(fn, a) }
  } else {
    h <- recurse(fname)
    function(a) .mma_call_head(h, a)
  }
  if (bare) {
    slots <- lapply(orders, function(o) Dummy())
    d <- do.call(Derivative, c(list(func(slots)), lapply(seq_along(slots), function(i) list(slots[[i]], orders[i]))))
    return(Lambda(slots, doit(d, deep = FALSE)))
  }
  args <- lapply(expr[-1L], recurse)
  if (length(orders) != length(args)) return(NULL)
  keys <- vapply(args, function(a) if (is_sym(a)) a$key else "", "")
  if (!anyDuplicated(keys) && all(vapply(args, function(a) is_sym(a) && a$type == "Symbol", TRUE))) {
    d <- do.call(Derivative, c(list(func(args)), lapply(seq_along(args), function(i) list(args[[i]], orders[i]))))
    return(doit(d, deep = FALSE))
  }
  slots <- lapply(args, function(a) Dummy())
  d <- do.call(Derivative, c(list(func(slots)), lapply(seq_along(slots), function(i) list(slots[[i]], orders[i]))))
  doit(Subs(d, do.call(Tuple, slots), do.call(Tuple, args)), deep = FALSE)
}

.mma_from_fullformlist_to_sympy <- function(ffl) {
  conv <- .mma_node_conversions()
  recurse <- function(expr) {
    if (is.list(expr)) {
      d <- .mma_derivative_to_sympy(expr, recurse)
      if (!is.null(d)) return(d)
      if (is.list(expr[[1L]])) head <- recurse(expr[[1L]])
      else {
        h <- expr[[1L]]
        head <- if (!is.null(conv[[h]])) conv[[h]] else function(a) fcall(h, lapply(a, .pe_to_sym))
      }
      return(.mma_call_head(head, lapply(expr[-1L], recurse)))
    }
    v <- .pe_sympify_string(expr)
    switch(expr, I = S_I, Pi = S_Pi, ExponentialE = S_E, ImaginaryI = S_I, ImaginaryJ = S_I, v)
  }
  recurse(ffl)
}
