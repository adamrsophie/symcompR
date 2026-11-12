# parse_expr and its token transformations, after sympy/parsing/sympy_parser.py.
# The input goes through a tokenizer for Python's expression tokens, the
# transformations of SymPy (each a function of the token list), and SymPy's
# untokenize; the resulting code is parsed by a small parser for Python
# expressions and evaluated with Python's semantics for SymPy objects and
# Python numbers, strings, tuples and lists. Names resolve as in Python's
# eval(code, global_dict, local_dict). Nothing is evaluated as R code.

.pe_err <- function(cls, message, bases = character()) {
  structure(class = c(cls, bases, "error", "condition"),
            list(message = paste(as.character(message), collapse = ""), call = NULL))
}
.token_error <- function(message = "") .pe_err("TokenError", message)
.pe_type_error <- function(message) .pe_err("TypeError", message)
.pe_value_error <- function(message) .pe_err("ValueError", message)
.pe_syntax_error <- function(message = "invalid syntax (<string>, line 1)") .pe_err("SyntaxError", message)

# ---------------------------------------------------------------------------
# Python values that are not SymPy objects

.pe_int <- function(v) structure(list(v = gmp::as.bigz(v)), class = "pe_int")
.pe_float <- function(v) structure(list(v = as.double(v)), class = "pe_float")
.pe_str_val <- function(v) structure(list(v = v), class = "pe_str")
.pe_none <- structure(list(), class = "pe_none")
.pe_tuple <- function(x) structure(x, class = "pe_tuple")
.pe_list <- function(x) structure(x, class = "pe_list")
.pe_callable <- function(f, name = "function") structure(f, class = c("pe_callable", "function"), pe_name = name)

.pe_float_repr <- function(v) {
  if (is.nan(v)) return("nan")
  if (is.infinite(v)) return(if (v > 0) "inf" else "-inf")
  s <- NULL
  for (d in 1:17) {
    s <- formatC(v, digits = d, format = "g")
    if (as.double(s) == v) break
  }
  s <- trimws(s)
  if (grepl("e", s)) {
    m <- regmatches(s, regexec("^(-?[0-9.]+)e([-+])0*([0-9]+)$", s))[[1L]]
    ex <- as.integer(paste0(m[3L], m[4L]))
    if (ex >= -5L && ex < 16L) s <- format(v, scientific = FALSE, digits = 17)
    else return(paste0(m[2L], "e", m[3L], sprintf("%02d", as.integer(m[4L]))))
  }
  if (!grepl("[.e]", s)) s <- paste0(s, ".0")
  s
}

.pe_repr_str <- function(s) {
  q <- if (grepl("'", s, fixed = TRUE) && !grepl("\"", s, fixed = TRUE)) "\"" else "'"
  out <- gsub("\\", "\\\\", s, fixed = TRUE)
  if (q == "'") out <- gsub("'", "\\'", out, fixed = TRUE)
  out <- gsub("\n", "\\n", out, fixed = TRUE)
  out <- gsub("\t", "\\t", out, fixed = TRUE)
  out <- gsub("\r", "\\r", out, fixed = TRUE)
  paste0(q, out, q)
}

# str() and repr() of a value as Python prints it
.pe_str <- function(v) .pe_fmt(v, FALSE)
.pe_repr <- function(v) .pe_fmt(v, TRUE)
.pe_fmt <- function(v, repr) {
  if (is.null(v) || inherits(v, "pe_none")) return("None")
  if (is_sym(v) || inherits(v, "symmatrix")) return(sstr(v))
  if (inherits(v, "lark_tree")) return(format(v))
  if (inherits(v, "pe_int")) return(as.character(v$v))
  if (inherits(v, "pe_float")) return(.pe_float_repr(v$v))
  if (inherits(v, "pe_str")) return(if (repr) .pe_repr_str(v$v) else v$v)
  if (is.character(v) && length(v) == 1L) return(if (repr) .pe_repr_str(v) else v)
  if (is.logical(v) && length(v) == 1L) return(if (v) "True" else "False")
  if (inherits(v, "sym_function")) return(attr(v, "sym_name"))
  if (inherits(v, "pe_callable")) return(attr(v, "pe_name"))
  if (inherits(v, "pe_tuple")) {
    el <- vapply(unclass(v), .pe_fmt, "", repr = TRUE)
    return(if (length(el) == 1L) paste0("(", el, ",)") else paste0("(", paste(el, collapse = ", "), ")"))
  }
  if (inherits(v, "pe_set")) return(if (length(v)) paste0("{", paste(vapply(unclass(v), .pe_fmt, "", repr = TRUE), collapse = ", "), "}") else "set()")
  if (is.list(v) && !is.object(v) || inherits(v, "pe_list"))
    return(paste0("[", paste(vapply(unclass(v), .pe_fmt, "", repr = TRUE), collapse = ", "), "]"))
  if (is.function(v)) return("<function>")
  paste(format(v), collapse = " ")
}

#' @export
print.pe_int <- function(x, ...) { cat(.pe_str(x), "\n"); invisible(x) }
#' @export
print.pe_float <- function(x, ...) { cat(.pe_str(x), "\n"); invisible(x) }
#' @export
print.pe_str <- function(x, ...) { cat(.pe_repr(x), "\n"); invisible(x) }
#' @export
print.pe_tuple <- function(x, ...) { cat(.pe_str(x), "\n"); invisible(x) }
#' @export
print.pe_list <- function(x, ...) { cat(.pe_str(x), "\n"); invisible(x) }

# ---------------------------------------------------------------------------
# Tokenizer: the tokens of Python 3.12's tokenize module for one logical line

.PY_OPS <- c("**=", "//=", ">>=", "<<=", "...", "->", ":=", "**", "//", "<<", ">>", "<=", ">=", "==", "!=",
             "+=", "-=", "*=", "/=", "%=", "&=", "|=", "^=", "@=")
.PY_NUMBER <- paste0(
  "^(?:(?:[0-9](?:_?[0-9])*[jJ]|(?:(?:[0-9](?:_?[0-9])*\\.(?:[0-9](?:_?[0-9])*)?|\\.[0-9](?:_?[0-9])*)",
  "(?:[eE][-+]?[0-9](?:_?[0-9])*)?|[0-9](?:_?[0-9])*[eE][-+]?[0-9](?:_?[0-9])*)[jJ])",
  "|(?:(?:[0-9](?:_?[0-9])*\\.(?:[0-9](?:_?[0-9])*)?|\\.[0-9](?:_?[0-9])*)(?:[eE][-+]?[0-9](?:_?[0-9])*)?",
  "|[0-9](?:_?[0-9])*[eE][-+]?[0-9](?:_?[0-9])*)",
  "|(?:0[xX](?:_?[0-9a-fA-F])+|0[bB](?:_?[01])+|0[oO](?:_?[0-7])+|0(?:_?0)*|[1-9](?:_?[0-9])*))")
.PY_NAME <- "^[\\p{L}\\p{Nl}_][\\p{L}\\p{Nl}\\p{Mn}\\p{Mc}\\p{Nd}\\p{Pc}]*"

.pe_tok <- function(type, value) c(type, value)

.pe_tokenize <- function(s) {
  n <- nchar(s)
  toks <- list()
  i <- 1L
  depth <- 0L
  add <- function(type, value) toks[[length(toks) + 1L]] <<- c(type, value)
  while (i <= n) {
    ch <- substr(s, i, i)
    if (ch %in% c(" ", "\t", "\f")) { i <- i + 1L; next }
    if (ch == "#") break
    if (ch == "\\") {
      if (i == n) stop(.token_error(sprintf("('unexpected EOF in multi-line statement', (1, %d))", i)))
      i <- i + 1L
      next
    }
    if (ch == "\n") {
      if (depth == 0L) add("NEWLINE", "\n")
      i <- i + 1L
      next
    }
    rest <- substring(s, i)
    if (grepl("^[0-9]|^\\.[0-9]", rest)) {
      m <- regmatches(rest, regexpr(.PY_NUMBER, rest, perl = TRUE))
      add("NUMBER", m)
      i <- i + nchar(m)
      next
    }
    m <- regmatches(rest, regexpr(.PY_NAME, rest, perl = TRUE))
    if (length(m) && nchar(m)) {
      nx <- substr(s, i + nchar(m), i + nchar(m))
      if (nx %in% c("'", "\"") && tolower(m) %in% c("r", "u", "b", "f", "br", "rb", "fr", "rf")) {
        j <- .pe_string_end(s, i + nchar(m), i)
        add("STRING", substr(s, i, j))
        i <- j + 1L
        next
      }
      add("NAME", m)
      i <- i + nchar(m)
      next
    }
    if (ch %in% c("'", "\"")) {
      j <- .pe_string_end(s, i, i)
      add("STRING", substr(s, i, j))
      i <- j + 1L
      next
    }
    op <- NULL
    for (o in .PY_OPS) if (startsWith(rest, o)) { op <- o; break }
    if (is.null(op)) op <- ch
    if (op %in% c("(", "[", "{")) depth <- depth + 1L
    if (op %in% c(")", "]", "}")) depth <- max(depth - 1L, 0L)
    add("OP", op)
    i <- i + nchar(op)
  }
  if (depth > 0L) stop(.token_error("('unexpected EOF in multi-line statement', (1, 0))"))
  add("NEWLINE", "")
  add("ENDMARKER", "")
  toks
}

.pe_string_end <- function(s, q, start) {
  n <- nchar(s)
  quote <- substr(s, q, q)
  triple <- substr(s, q, q + 2L) == strrep(quote, 3L)
  if (triple) {
    j <- q + 3L
    while (j <= n - 2L) {
      if (substr(s, j, j) == "\\") { j <- j + 2L; next }
      if (substr(s, j, j + 2L) == strrep(quote, 3L)) return(j + 2L)
      j <- j + 1L
    }
    stop(.token_error(sprintf("('EOF in multi-line string', (1, %d))", start)))
  }
  j <- q + 1L
  while (j <= n) {
    c1 <- substr(s, j, j)
    if (c1 == "\\") { j <- j + 2L; next }
    if (c1 == quote) return(j)
    if (c1 == "\n") break
    j <- j + 1L
  }
  stop(.token_error(sprintf("('unterminated string literal (detected at line 1)', (1, %d))", start)))
}

# untokenize in compat mode: NAME and NUMBER get a space after them, two
# adjacent strings a space between them
.pe_untokenize <- function(toks) {
  out <- character(length(toks))
  prevstring <- FALSE
  for (k in seq_along(toks)) {
    t <- toks[[k]]
    v <- t[2L]
    if (t[1L] %in% c("NAME", "NUMBER")) v <- paste0(v, " ")
    if (t[1L] == "STRING") {
      if (prevstring) v <- paste0(" ", v)
      prevstring <- TRUE
    } else prevstring <- FALSE
    out[k] <- v
  }
  paste(out, collapse = "")
}

# ---------------------------------------------------------------------------
# The namespaces: local_dict is an environment; the default global_dict is
# SymPy's (from sympy import *, the builtin functions, max = Max, min = Min).

.pe_ns <- function(local_dict, global_dict) {
  e <- new.env(parent = emptyenv())
  e$local <- local_dict
  e$global <- global_dict
  e
}
.pe_in_local <- function(ns, name) exists(name, envir = ns$local, inherits = FALSE)
.pe_in_global <- function(ns, name) {
  if (is.null(ns$global)) name %in% .PY_GLOBAL_OBJECTS || name %in% .PY_GLOBAL_OTHER
  else name %in% names(ns$global)
}
.pe_global_object <- function(ns, name) {
  if (is.null(ns$global)) return(name %in% .PY_GLOBAL_OBJECTS)
  v <- ns$global[[name]]
  is_sym(v) || is.function(v) || inherits(v, "symmatrix")
}
.pe_is_callable_value <- function(v) (is.function(v) || inherits(v, "pe_callable")) && !is_sym(v)

.pe_token_callable <- function(tok, ns) {
  name <- tok[2L]
  if (.pe_in_local(ns, name)) {
    v <- get(name, envir = ns$local)
    if (!(is.null(v) || (is.logical(v) && !isTRUE(v)))) return(.pe_is_callable_value(v))
  }
  if (is.null(ns$global)) return(name %in% .PY_GLOBAL_CALLABLE)
  .pe_is_callable_value(ns$global[[name]])
}

# ---------------------------------------------------------------------------
# Token helpers. A token is c(type, value); grouped tokens are a
# ParenthesisGroup (class pe_pg, a list of tokens) or an AppliedFunction
# (class pe_af with fields fn, args, exponent).

.tk_type <- function(t) if (inherits(t, "pe_af")) "<AppliedFunction>" else if (inherits(t, "pe_pg")) "<ParenthesisGroup>" else t[1L]
.tk_val <- function(t) if (is.character(t)) t[2L] else "<group>"
.tk_is <- function(t, type, value) is.character(t) && t[1L] == type && t[2L] == value

.pe_add_factorial_tokens <- function(name, result) {
  if (!length(result) || .tk_val(result[[length(result)]]) == "(") stop(.token_error())
  beginning <- list(c("NAME", name), c("OP", "("))
  end <- list(c("OP", ")"))
  diff <- 0L
  len <- length(result)
  for (index in 0:(len - 1L)) {
    i <- len - index
    tokval <- .tk_val(result[[i]])
    if (tokval == ")") diff <- diff + 1L
    else if (tokval == "(") diff <- diff - 1L
    if (diff == 0L) {
      if (i - 1L >= 1L && .tk_type(result[[i - 1L]]) == "NAME")
        return(c(result[seq_len(i - 2L)], beginning, result[(i - 1L):len], end))
      return(c(result[seq_len(i - 1L)], beginning, result[i:len], end))
    }
  }
  result
}

.pe_group_parentheses <- function(recursor) {
  function(tokens, local_dict, global_dict) {
    result <- list()
    stacks <- list()
    for (token in tokens) {
      if (.tk_type(token) == "OP") {
        if (.tk_val(token) == "(") {
          stacks[[length(stacks) + 1L]] <- list()
        } else if (.tk_val(token) == ")") {
          if (!length(stacks)) stop(.token_error("Mismatched parentheses"))
          top <- c(stacks[[length(stacks)]], list(token))
          stacks[[length(stacks)]] <- NULL
          if (length(stacks)) {
            stacks[[length(stacks)]] <- c(stacks[[length(stacks)]], top)
          } else {
            inner <- recursor(top[seq_len(length(top) - 2L) + 1L], local_dict, global_dict)
            result[[length(result) + 1L]] <- structure(c(top[1L], inner, top[length(top)]), class = "pe_pg")
          }
          next
        }
      }
      if (length(stacks)) stacks[[length(stacks)]] <- c(stacks[[length(stacks)]], list(token))
      else result[[length(result) + 1L]] <- token
    }
    if (length(stacks)) stop(.token_error("Mismatched parentheses"))
    result
  }
}

.pe_apply_functions <- function(tokens, local_dict, global_dict) {
  ns <- .pe_ns(local_dict, global_dict)
  result <- list()
  symbol <- NULL
  for (tok in tokens) {
    if (inherits(tok, "pe_pg")) {
      if (!is.null(symbol) && .pe_token_callable(symbol, ns)) {
        result[[length(result)]] <- structure(list(fn = symbol, args = tok, exponent = list()), class = "pe_af")
        symbol <- NULL
      } else result <- c(result, unclass(tok))
    } else if (.tk_type(tok) == "NAME") {
      symbol <- tok
      result[[length(result) + 1L]] <- tok
    } else {
      symbol <- NULL
      result[[length(result) + 1L]] <- tok
    }
  }
  result
}

.pe_flatten <- function(result) {
  out <- list()
  for (tok in result) {
    if (inherits(tok, "pe_af")) out <- c(out, list(tok$fn), unclass(tok$args))
    else out[[length(out) + 1L]] <- tok
  }
  out
}

.pe_implicit_multiplication <- function(tokens, local_dict, global_dict) {
  ns <- .pe_ns(local_dict, global_dict)
  result <- list()
  skip <- FALSE
  star <- c("OP", "*")
  n <- length(tokens)
  for (k in seq_len(max(n - 1L, 0L))) {
    tok <- tokens[[k]]; nextTok <- tokens[[k + 1L]]
    result[[length(result) + 1L]] <- tok
    if (skip) { skip <- FALSE; next }
    if (.tk_is(tok, "OP", ".") && .tk_type(nextTok) == "NAME") { skip <- TRUE; next }
    if (inherits(tok, "pe_af")) {
      if (inherits(nextTok, "pe_af")) result[[length(result) + 1L]] <- star
      else if (.tk_is(nextTok, "OP", "(")) {
        if (tok$fn[2L] == "Function") result[[length(result)]]$fn <- c(tok$fn[1L], "Symbol")
        result[[length(result) + 1L]] <- star
      } else if (.tk_type(nextTok) == "NAME") result[[length(result) + 1L]] <- star
    } else if (.tk_is(tok, "OP", ")")) {
      if (inherits(nextTok, "pe_af") || .tk_type(nextTok) == "NAME" || .tk_is(nextTok, "OP", "("))
        result[[length(result) + 1L]] <- star
    } else if (.tk_type(tok) == "NAME" && !.pe_token_callable(tok, ns)) {
      if (inherits(nextTok, "pe_af") || (.tk_type(nextTok) == "NAME" && .pe_token_callable(nextTok, ns)))
        result[[length(result) + 1L]] <- star
      else if (.tk_is(nextTok, "OP", "(")) result[[length(result) + 1L]] <- star
      else if (.tk_type(nextTok) == "NAME") result[[length(result) + 1L]] <- star
    }
  }
  if (n) result[[length(result) + 1L]] <- tokens[[n]]
  result
}

.pe_implicit_application <- function(tokens, local_dict, global_dict) {
  ns <- .pe_ns(local_dict, global_dict)
  result <- list()
  appendParen <- 0L
  skip <- 0L
  exponentSkip <- FALSE
  n <- length(tokens)
  for (k in seq_len(max(n - 1L, 0L))) {
    tok <- tokens[[k]]; nextTok <- tokens[[k + 1L]]
    result[[length(result) + 1L]] <- tok
    if (.tk_type(tok) == "NAME" && !(.tk_type(nextTok) %in% c("OP", "ENDMARKER", "NEWLINE"))) {
      if (.pe_token_callable(tok, ns)) {
        result[[length(result) + 1L]] <- c("OP", "(")
        appendParen <- appendParen + 1L
      }
    } else if (.tk_type(tok) == "NAME" && .tk_type(nextTok) == "OP" && .tk_val(nextTok) == "**") {
      if (.pe_token_callable(tok, ns)) exponentSkip <- TRUE
    } else if (exponentSkip) {
      if (inherits(tok, "pe_af") || .tk_is(tok, "OP", "*")) {
        if (!.tk_is(nextTok, "OP", "*")) {
          if (!.tk_is(nextTok, "OP", "(")) {
            result[[length(result) + 1L]] <- c("OP", "(")
            appendParen <- appendParen + 1L
          }
          exponentSkip <- FALSE
        }
      }
    } else if (appendParen) {
      if (.tk_type(nextTok) == "OP" && .tk_val(nextTok) %in% c("^", "**", "*")) {
        skip <- 1L
        next
      }
      if (skip) { skip <- skip - 1L; next }
      result[[length(result) + 1L]] <- c("OP", ")")
      appendParen <- appendParen - 1L
    }
  }
  if (n) result[[length(result) + 1L]] <- tokens[[n]]
  if (appendParen) result <- c(result, rep(list(c("OP", ")")), appendParen))
  result
}

#' Token transformations of parse_expr
#'
#' The transformations of SymPy's `parse_expr()`, each a function of the
#' token list, the local namespace and the global namespace that returns
#' the transformed token list. `standard_transformations` is the default
#' list: `lambda_notation`, `auto_symbol`, `repeated_decimals`,
#' `auto_number` and `factorial_notation`. `convert_xor` reads `^` as a
#' power, `implicit_multiplication` and `implicit_application` add the
#' omitted `*` and parentheses, `split_symbols` splits names like `xy`,
#' `function_exponentiation` reads `sin**2(x)`, `t_rationalize` (SymPy's
#' `rationalize`) makes floats
#' exact and `convert_equals_signs` reads `=` as [Eq()].
#' `split_symbols_custom(predicate)` builds a symbol splitter.
#'
#' @param tokens A list of tokens, each `c(type, value)`.
#' @param local_dict An environment holding the local names.
#' @param global_dict A named list of global names, or `NULL` for SymPy's.
#' @param predicate A function of a name returning `TRUE` when it is split.
#' @return The transformed token list; `split_symbols_custom()` returns a
#'   transformation.
#' @examples
#' parse_expr("2x", transformations = c(standard_transformations, implicit_multiplication))
#' parse_expr("x^2", transformations = c(standard_transformations, convert_xor))
#' @name parse_expr_transformations
NULL

#' @rdname parse_expr_transformations
#' @export
function_exponentiation <- function(tokens, local_dict, global_dict) {
  ns <- .pe_ns(local_dict, global_dict)
  result <- list()
  exponent <- list()
  consuming_exponent <- FALSE
  level <- 0L
  n <- length(tokens)
  for (k in seq_len(max(n - 1L, 0L))) {
    tok <- tokens[[k]]; nextTok <- tokens[[k + 1L]]
    if (.tk_type(tok) == "NAME" && .tk_type(nextTok) == "OP" && .tk_val(nextTok) == "**") {
      if (.pe_token_callable(tok, ns)) consuming_exponent <- TRUE
    } else if (consuming_exponent) {
      if (.tk_is(tok, "NAME", "Function")) tok <- c("NAME", "Symbol")
      exponent[[length(exponent) + 1L]] <- tok
      if (.tk_type(tok) == "OP" && .tk_type(nextTok) == "OP" && .tk_val(tok) == ")" && .tk_val(nextTok) == "(")
        consuming_exponent <- FALSE
      if (.tk_type(tok) == "OP" && .tk_type(nextTok) == "OP" && .tk_val(tok) == "*" && .tk_val(nextTok) == "(") {
        consuming_exponent <- FALSE
        exponent[[length(exponent)]] <- NULL
      }
      next
    } else if (length(exponent) && !consuming_exponent) {
      if (.tk_type(tok) == "OP") {
        if (.tk_val(tok) == "(") level <- level + 1L
        else if (.tk_val(tok) == ")") level <- level - 1L
      }
      if (level == 0L) {
        result <- c(result, list(tok), exponent)
        exponent <- list()
        next
      }
    }
    result[[length(result) + 1L]] <- tok
  }
  if (n) result[[length(result) + 1L]] <- tokens[[n]]
  if (length(exponent)) result <- c(result, exponent)
  result
}

.PY_GREEK_LETTERS <- c("alpha", "beta", "chi", "delta", "digamma", "epsilon", "eta", "gamma", "heta", "iota",
                       "kappa", "koppa", "lamda", "mu", "nu", "omega", "omicron", "phi", "pi", "psi", "rho",
                       "sampi", "san", "sho", "sigma", "stigma", "tau", "theta", "upsilon", "xi", "zeta")

.token_splittable <- function(token_name) {
  if (grepl("_", token_name, fixed = TRUE)) return(FALSE)
  if (tolower(token_name) %in% .PY_GREEK_LETTERS) return(FALSE)
  nchar(token_name) > 1L
}

#' @rdname parse_expr_transformations
#' @export
split_symbols_custom <- function(predicate) {
  force(predicate)
  function(tokens, local_dict, global_dict) {
    ns <- .pe_ns(local_dict, global_dict)
    result <- list()
    split <- FALSE
    split_previous <- FALSE
    for (tok in tokens) {
      if (split_previous) { split_previous <- FALSE; next }
      split_previous <- FALSE
      if (.tk_type(tok) == "NAME" && .tk_val(tok) %in% c("Symbol", "Function")) {
        split <- TRUE
      } else if (split && .tk_type(tok) == "NAME") {
        symbol <- substr(tok[2L], 2L, nchar(tok[2L]) - 1L)
        if (isTRUE(predicate(symbol))) {
          result <- result[seq_len(length(result) - 2L)]
          chars <- strsplit(symbol, "", fixed = TRUE)[[1L]]
          i <- 1L
          while (i <= length(chars)) {
            char <- chars[i]
            if (.pe_in_local(ns, char) || .pe_in_global(ns, char)) {
              result[[length(result) + 1L]] <- c("NAME", char)
            } else if (grepl("^[0-9]$", char)) {
              digits <- char
              while (i < length(chars) && grepl("^[0-9]$", chars[i + 1L])) {
                i <- i + 1L
                digits <- paste0(digits, chars[i])
              }
              result <- c(result, list(c("NAME", "Number"), c("OP", "("), c("NAME", paste0("'", digits, "'")), c("OP", ")")))
            } else {
              result <- c(result, list(c("NAME", "Symbol"), c("OP", "("), c("NAME", paste0("'", char, "'")), c("OP", ")")))
            }
            i <- i + 1L
          }
          split <- FALSE
          split_previous <- TRUE
          next
        } else split <- FALSE
      }
      result[[length(result) + 1L]] <- tok
    }
    result
  }
}

#' @rdname parse_expr_transformations
#' @export
split_symbols <- split_symbols_custom(.token_splittable)

#' @rdname parse_expr_transformations
#' @export
implicit_multiplication <- function(tokens, local_dict, global_dict) {
  res1 <- .pe_group_parentheses(implicit_multiplication)(tokens, local_dict, global_dict)
  res2 <- .pe_apply_functions(res1, local_dict, global_dict)
  res3 <- .pe_implicit_multiplication(res2, local_dict, global_dict)
  .pe_flatten(res3)
}

#' @rdname parse_expr_transformations
#' @export
implicit_application <- function(tokens, local_dict, global_dict) {
  res1 <- .pe_group_parentheses(implicit_application)(tokens, local_dict, global_dict)
  res2 <- .pe_apply_functions(res1, local_dict, global_dict)
  res3 <- .pe_implicit_application(res2, local_dict, global_dict)
  .pe_flatten(res3)
}

#' @rdname parse_expr_transformations
#' @export
implicit_multiplication_application <- function(tokens, local_dict, global_dict) {
  for (step in list(split_symbols, implicit_multiplication, implicit_application, function_exponentiation))
    tokens <- step(tokens, local_dict, global_dict)
  tokens
}

.PY_KEYWORDS <- c("False", "None", "True", "and", "as", "assert", "async", "await", "break", "class",
                  "continue", "def", "del", "elif", "else", "except", "finally", "for", "from", "global",
                  "if", "import", "in", "is", "lambda", "nonlocal", "not", "or", "pass", "raise", "return",
                  "try", "while", "with", "yield")

#' @rdname parse_expr_transformations
#' @export
auto_symbol <- function(tokens, local_dict, global_dict) {
  ns <- .pe_ns(local_dict, global_dict)
  result <- list()
  prevTok <- c("-1", "")
  tokens <- c(tokens, list(c("-1", "")))
  for (k in seq_len(length(tokens) - 1L)) {
    tok <- tokens[[k]]; nextTok <- tokens[[k + 1L]]
    if (.tk_type(tok) == "NAME") {
      name <- tok[2L]
      nextVal <- .tk_val(nextTok)
      if (name %in% .PY_KEYWORDS ||
          (prevTok[1L] == "OP" && prevTok[2L] == ".") ||
          (prevTok[1L] == "OP" && prevTok[2L] %in% c("(", ",") && .tk_type(nextTok) == "OP" && nextVal == "=") ||
          .pe_in_local(ns, name)) {
        result[[length(result) + 1L]] <- c("NAME", name)
        next
      } else if (.pe_in_global(ns, name) && .pe_global_object(ns, name)) {
        result[[length(result) + 1L]] <- c("NAME", name)
        next
      }
      result <- c(result, list(c("NAME", if (nextVal != "(") "Symbol" else "Function"), c("OP", "("),
                               c("NAME", .pe_repr_str(name)), c("OP", ")")))
    } else result[[length(result) + 1L]] <- tok
    prevTok <- if (is.character(tok)) tok else c("<group>", "")
  }
  result
}

#' @rdname parse_expr_transformations
#' @export
lambda_notation <- function(tokens, local_dict, global_dict) {
  first <- tokens[[1L]]
  if (.tk_is(first, "NAME", "lambda")) {
    n <- length(tokens)
    if (n == 2L || (n == 3L && .tk_type(tokens[[2L]]) == "NEWLINE")) return(tokens)
    if (n > 2L) {
      result <- list(c("NAME", "Lambda"), c("OP", "("), c("OP", "("), c("OP", ")"), c("OP", ")"))
      flag <- FALSE
      for (tok in tokens[-1L]) {
        if (.tk_is(tok, "OP", ":")) { tok <- c("OP", ","); flag <- TRUE }
        if (!flag && .tk_type(tok) == "OP" && .tk_val(tok) %in% c("*", "**"))
          stop(.token_error("Starred arguments in lambda not supported"))
        pos <- if (flag) length(result) - 1L else length(result) - 2L
        result <- append(result, list(tok), after = pos)
      }
      return(result)
    }
  }
  tokens
}

#' @rdname parse_expr_transformations
#' @export
factorial_notation <- function(tokens, local_dict, global_dict) {
  result <- list()
  nfactorial <- 0L
  for (tok in tokens) {
    if (.tk_is(tok, "OP", "!")) nfactorial <- nfactorial + 1L
    else if (.tk_type(tok) == "ERRORTOKEN") {
      if (tok[2L] == "!") nfactorial <- nfactorial + 1L
      else { nfactorial <- 0L; result[[length(result) + 1L]] <- c("OP", tok[2L]) }
    } else {
      if (nfactorial == 1L) result <- .pe_add_factorial_tokens("factorial", result)
      else if (nfactorial == 2L) result <- .pe_add_factorial_tokens("factorial2", result)
      else if (nfactorial > 2L) stop(.token_error())
      nfactorial <- 0L
      result[[length(result) + 1L]] <- tok
    }
  }
  result
}

#' @rdname parse_expr_transformations
#' @export
convert_xor <- function(tokens, local_dict, global_dict) {
  lapply(tokens, function(t) if (.tk_is(t, "OP", "^")) c("OP", "**") else t)
}

#' @rdname parse_expr_transformations
#' @export
repeated_decimals <- function(tokens, local_dict, global_dict) {
  result <- list()
  is_digit <- function(s) grepl("^[0-9_]*$", s)
  num <- list()
  for (tok in tokens) {
    toknum <- .tk_type(tok); tokval <- .tk_val(tok)
    if (toknum == "NUMBER") {
      if (!length(num) && grepl(".", tokval, fixed = TRUE) && !grepl("e", tolower(tokval), fixed = TRUE) &&
          !grepl("j", tolower(tokval), fixed = TRUE)) num[[length(num) + 1L]] <- tok
      else if (is_digit(tokval) && (length(num) == 2L || (length(num) == 3L && is_digit(num[[3L]][2L]))))
        num[[length(num) + 1L]] <- tok
      else num <- list()
    } else if (toknum == "OP") {
      if (tokval == "[" && length(num) == 1L) num[[length(num) + 1L]] <- c("OP", tokval)
      else if (tokval == "]" && length(num) >= 3L) num[[length(num) + 1L]] <- c("OP", tokval)
      else if (tokval == "." && !length(num)) num[[length(num) + 1L]] <- c("NUMBER", "0.")
      else num <- list()
    } else num <- list()
    result[[length(result) + 1L]] <- tok
    if (length(num) && num[[length(num)]][2L] == "]") {
      result <- result[seq_len(length(result) - length(num))]
      parts <- strsplit(paste0(num[[1L]][2L], " "), ".", fixed = TRUE)[[1L]]
      pre <- parts[1L]; post <- trimws(parts[2L])
      repetend <- num[[3L]][2L]
      if (length(num) == 5L) repetend <- paste0(repetend, num[[4L]][2L])
      pre <- gsub("_", "", pre, fixed = TRUE); post <- gsub("_", "", post, fixed = TRUE)
      repetend <- gsub("_", "", repetend, fixed = TRUE)
      zeros <- strrep("0", nchar(post))
      post2 <- sub("^0+", "", post); repetends <- sub("^0+", "", repetend)
      a <- if (nchar(pre)) pre else "0"
      b <- if (nchar(post2)) post2 else "0"; cc <- paste0("1", zeros)
      d <- repetends; e <- paste0(strrep("9", nchar(repetend)), zeros)
      result <- c(result, list(c("OP", "("), c("NAME", "Integer"), c("OP", "("), c("NUMBER", a), c("OP", ")"),
                               c("OP", "+"), c("NAME", "Rational"), c("OP", "("), c("NUMBER", b), c("OP", ","),
                               c("NUMBER", cc), c("OP", ")"), c("OP", "+"), c("NAME", "Rational"), c("OP", "("),
                               c("NUMBER", d), c("OP", ","), c("NUMBER", e), c("OP", ")"), c("OP", ")")))
      num <- list()
    }
  }
  result
}

#' @rdname parse_expr_transformations
#' @export
auto_number <- function(tokens, local_dict, global_dict) {
  result <- list()
  for (tok in tokens) {
    if (.tk_type(tok) == "NUMBER") {
      number <- tok[2L]
      postfix <- list()
      if (grepl("[jJ]$", number)) {
        number <- substr(number, 1L, nchar(number) - 1L)
        postfix <- list(c("OP", "*"), c("NAME", "I"))
      }
      if (grepl(".", number, fixed = TRUE) || (grepl("[eE]", number) && !grepl("^0[xX]", number)))
        seq_ <- list(c("NAME", "Float"), c("OP", "("), c("NUMBER", .pe_repr_str(number)), c("OP", ")"))
      else seq_ <- list(c("NAME", "Integer"), c("OP", "("), c("NUMBER", number), c("OP", ")"))
      result <- c(result, seq_, postfix)
    } else result[[length(result) + 1L]] <- tok
  }
  result
}

#' @rdname parse_expr_transformations
#' @export
t_rationalize <- function(tokens, local_dict, global_dict) {
  result <- list()
  passed_float <- FALSE
  for (tok in tokens) {
    if (.tk_type(tok) == "NAME") {
      if (tok[2L] == "Float") { passed_float <- TRUE; tok <- c("NAME", "Rational") }
      result[[length(result) + 1L]] <- tok
    } else if (passed_float && .tk_type(tok) == "NUMBER") {
      passed_float <- FALSE
      result[[length(result) + 1L]] <- c("STRING", tok[2L])
    } else result[[length(result) + 1L]] <- tok
  }
  result
}

.pe_transform_equals_sign <- function(tokens, local_dict, global_dict) {
  if (!any(vapply(tokens, .tk_is, TRUE, type = "OP", value = "="))) return(tokens)
  result <- list(c("NAME", "Eq"), c("OP", "("))
  for (tok in tokens) result[[length(result) + 1L]] <- if (.tk_is(tok, "OP", "=")) c("OP", ",") else tok
  c(result, list(c("OP", ")")))
}

#' @rdname parse_expr_transformations
#' @export
convert_equals_signs <- function(tokens, local_dict, global_dict) {
  res1 <- .pe_group_parentheses(convert_equals_signs)(tokens, local_dict, global_dict)
  res2 <- .pe_apply_functions(res1, local_dict, global_dict)
  res3 <- .pe_transform_equals_sign(res2, local_dict, global_dict)
  .pe_flatten(res3)
}

#' @rdname parse_expr_transformations
#' @export
standard_transformations <- list(lambda_notation, auto_symbol, repeated_decimals, auto_number, factorial_notation)

.pe_transformation_list <- function() {
  list(lambda_notation, auto_symbol, repeated_decimals, auto_number, factorial_notation,
       implicit_multiplication_application, convert_xor, implicit_application, implicit_multiplication,
       convert_equals_signs, function_exponentiation, t_rationalize)
}

#' Select parse_expr transformations by number
#'
#' The `T[...]` shortcut of SymPy's parser: the transformations at the given
#' 0-based positions in SymPy's numbered list (0 `lambda_notation`, 1
#' `auto_symbol`, 2 `repeated_decimals`, 3 `auto_number`, 4
#' `factorial_notation`, 5 `implicit_multiplication_application`, 6
#' `convert_xor`, 7 `implicit_application`, 8 `implicit_multiplication`, 9
#' `convert_equals_signs`, 10 `function_exponentiation`, 11 `t_rationalize`).
#' Negative positions count from the end, as in Python.
#'
#' @param ... Integer positions or vectors of positions.
#' @return A list of transformations.
#' @examples
#' parse_expr("2x", transformations = parser_T(0:4, 8))
#' @export
parser_T <- function(...) {
  t <- .pe_transformation_list()
  idx <- unlist(list(...))
  lapply(idx, function(i) {
    i <- as.integer(i)
    if (i < -length(t) || i >= length(t)) stop(.pe_err("IndexError", "range object index out of range"))
    t[[if (i < 0L) length(t) + i + 1L else i + 1L]]
  })
}

# ---------------------------------------------------------------------------
# A parser for Python expressions over the tokens of the transformed code

.pe_parse <- function(code) {
  toks <- Filter(function(t) !(t[1L] %in% c("NEWLINE", "ENDMARKER")), .pe_tokenize(trimws(code)))
  pos <- 1L
  n <- length(toks)
  peek <- function(k = 0L) if (pos + k <= n) toks[[pos + k]] else c("END", "")
  is_op <- function(v, k = 0L) { t <- peek(k); t[1L] == "OP" && t[2L] %in% v }
  is_kw <- function(v, k = 0L) { t <- peek(k); t[1L] == "NAME" && t[2L] %in% v }
  fail <- function() stop(.pe_syntax_error())
  expect_op <- function(v) { if (!is_op(v)) fail(); pos <<- pos + 1L }
  advance <- function() { t <- peek(); pos <<- pos + 1L; t }

  p_expressions <- function() {
    e <- p_expr()
    if (!is_op(",")) return(e)
    elts <- list(e)
    while (is_op(",")) {
      pos <<- pos + 1L
      if (pos > n || is_op(c(")", "]", "}"))) break
      elts[[length(elts) + 1L]] <- p_expr()
    }
    list(k = "Tuple", elts = elts)
  }
  p_expr <- function() {
    if (is_kw("lambda")) {
      pos <<- pos + 1L
      params <- character()
      while (!is_op(":")) {
        t <- advance()
        if (t[1L] != "NAME") fail()
        params <- c(params, t[2L])
        if (is_op(",")) pos <<- pos + 1L else if (!is_op(":")) fail()
      }
      pos <<- pos + 1L
      return(list(k = "Lambda", params = params, body = p_expr()))
    }
    e <- p_or()
    if (is_kw("if")) {
      pos <<- pos + 1L
      test <- p_or()
      if (!is_kw("else")) fail()
      pos <<- pos + 1L
      return(list(k = "IfExp", test = test, body = e, orelse = p_expr()))
    }
    e
  }
  p_or <- function() {
    e <- p_and()
    if (!is_kw("or")) return(e)
    vals <- list(e)
    while (is_kw("or")) { pos <<- pos + 1L; vals[[length(vals) + 1L]] <- p_and() }
    list(k = "BoolOp", op = "Or", values = vals)
  }
  p_and <- function() {
    e <- p_not()
    if (!is_kw("and")) return(e)
    vals <- list(e)
    while (is_kw("and")) { pos <<- pos + 1L; vals[[length(vals) + 1L]] <- p_not() }
    list(k = "BoolOp", op = "And", values = vals)
  }
  p_not <- function() {
    if (is_kw("not")) { pos <<- pos + 1L; return(list(k = "UnaryOp", op = "Not", x = p_not())) }
    p_compare()
  }
  p_compare <- function() {
    left <- p_bitor()
    ops <- character(); comps <- list()
    repeat {
      op <- NULL
      if (is_op(c("==", "!=", "<", "<=", ">", ">="))) { op <- peek()[2L]; pos <<- pos + 1L }
      else if (is_kw("in")) { op <- "in"; pos <<- pos + 1L }
      else if (is_kw("not") && is_kw("in", 1L)) { op <- "not in"; pos <<- pos + 2L }
      else if (is_kw("is")) {
        pos <<- pos + 1L
        if (is_kw("not")) { pos <<- pos + 1L; op <- "is not" } else op <- "is"
      }
      if (is.null(op)) break
      ops <- c(ops, op); comps[[length(comps) + 1L]] <- p_bitor()
    }
    if (!length(ops)) return(left)
    list(k = "Compare", left = left, ops = ops, comps = comps)
  }
  binlevel <- function(next_, opsv) {
    function() {
      e <- next_()
      while (is_op(opsv)) {
        op <- advance()[2L]
        e <- list(k = "BinOp", op = op, l = e, r = next_())
      }
      e
    }
  }
  p_factor <- function() {
    if (is_op(c("+", "-", "~"))) {
      op <- advance()[2L]
      return(list(k = "UnaryOp", op = switch(op, "+" = "UAdd", "-" = "USub", "~" = "Invert"), x = p_factor()))
    }
    p_power()
  }
  p_term <- binlevel(p_factor, c("*", "/", "//", "%", "@"))
  p_sum <- binlevel(p_term, c("+", "-"))
  p_shift <- binlevel(p_sum, c("<<", ">>"))
  p_bitand <- binlevel(p_shift, "&")
  p_bitxor <- binlevel(p_bitand, "^")
  p_bitor <- binlevel(p_bitxor, "|")
  p_power <- function() {
    e <- p_primary()
    if (is_op("**")) { pos <<- pos + 1L; return(list(k = "BinOp", op = "**", l = e, r = p_factor())) }
    e
  }
  p_args <- function() {
    args <- list(); kw <- list()
    while (!is_op(")")) {
      if (is_op("*")) { pos <<- pos + 1L; args[[length(args) + 1L]] <- list(k = "Starred", v = p_expr()) }
      else if (is_op("**")) { pos <<- pos + 1L; kw[[length(kw) + 1L]] <- list(name = NA_character_, v = p_expr()) }
      else if (peek()[1L] == "NAME" && is_op("=", 1L)) {
        nm <- advance()[2L]; pos <<- pos + 1L
        kw[[length(kw) + 1L]] <- list(name = nm, v = p_expr())
      } else {
        if (length(kw)) fail()
        args[[length(args) + 1L]] <- p_expr()
      }
      if (is_op(",")) pos <<- pos + 1L else if (!is_op(")")) fail()
    }
    pos <<- pos + 1L
    list(args = args, kw = kw)
  }
  p_primary <- function() {
    e <- p_atom()
    repeat {
      if (is_op("(")) {
        pos <<- pos + 1L
        a <- p_args()
        e <- list(k = "Call", f = e, args = a$args, kw = a$kw)
      } else if (is_op("[")) {
        pos <<- pos + 1L
        idx <- p_expressions()
        expect_op("]")
        e <- list(k = "Subscript", v = e, idx = idx)
      } else if (is_op(".") && peek(1L)[1L] == "NAME") {
        pos <<- pos + 1L
        e <- list(k = "Attr", v = e, attr = advance()[2L])
      } else break
    }
    e
  }
  p_atom <- function() {
    t <- peek()
    if (t[1L] == "NAME") {
      if (t[2L] %in% c("True", "False", "None")) {
        pos <<- pos + 1L
        return(list(k = "Const", v = switch(t[2L], True = TRUE, False = FALSE, None = .pe_none)))
      }
      if (t[2L] %in% .PY_KEYWORDS) fail()
      pos <<- pos + 1L
      return(list(k = "Name", id = t[2L]))
    }
    if (t[1L] == "NUMBER") { pos <<- pos + 1L; return(list(k = "Const", v = .pe_number_literal(t[2L]))) }
    if (t[1L] == "STRING") {
      s <- ""
      while (peek()[1L] == "STRING") s <- paste0(s, .pe_string_literal(advance()[2L]))
      return(list(k = "Const", v = .pe_str_val(s)))
    }
    if (is_op("(")) {
      pos <<- pos + 1L
      if (is_op(")")) { pos <<- pos + 1L; return(list(k = "Tuple", elts = list())) }
      e <- p_expr()
      if (is_op(",")) {
        elts <- list(e)
        while (is_op(",")) {
          pos <<- pos + 1L
          if (is_op(")")) break
          elts[[length(elts) + 1L]] <- p_expr()
        }
        e <- list(k = "Tuple", elts = elts)
      }
      expect_op(")")
      return(e)
    }
    if (is_op("[")) {
      pos <<- pos + 1L
      elts <- list()
      while (!is_op("]")) {
        elts[[length(elts) + 1L]] <- p_expr()
        if (is_op(",")) pos <<- pos + 1L else if (!is_op("]")) fail()
      }
      pos <<- pos + 1L
      return(list(k = "List", elts = elts))
    }
    if (is_op("{")) {
      pos <<- pos + 1L
      elts <- list(); keys <- list()
      while (!is_op("}")) {
        e <- p_expr()
        if (is_op(":")) { pos <<- pos + 1L; keys[[length(keys) + 1L]] <- e; e <- p_expr() }
        elts[[length(elts) + 1L]] <- e
        if (is_op(",")) pos <<- pos + 1L else if (!is_op("}")) fail()
      }
      pos <<- pos + 1L
      return(if (length(keys)) list(k = "Dict", keys = keys, values = elts) else list(k = "Set", elts = elts))
    }
    fail()
  }
  if (!n) fail()
  e <- p_expressions()
  if (pos <= n) fail()
  e
}

.pe_number_literal <- function(s) {
  s <- gsub("_", "", s, fixed = TRUE)
  if (grepl("[jJ]$", s)) return(structure(list(re = 0, im = as.double(sub("[jJ]$", "", s))), class = "pe_complex"))
  if (grepl("^0[xXoObB]", s) || !grepl("[.eE]", s)) return(.pe_int(.str_bigz(s)))
  .pe_float(as.double(s))
}

.pe_string_literal <- function(s) {
  prefix <- tolower(sub("^([A-Za-z]*).*$", "\\1", s))
  body <- substring(s, nchar(prefix) + 1L)
  q <- if (startsWith(body, "'''") || startsWith(body, "\"\"\"")) 3L else 1L
  body <- substr(body, q + 1L, nchar(body) - q)
  if (grepl("r", prefix, fixed = TRUE)) return(body)
  out <- character()
  i <- 1L; n <- nchar(body)
  while (i <= n) {
    ch <- substr(body, i, i)
    if (ch == "\\" && i < n) {
      nx <- substr(body, i + 1L, i + 1L)
      rep_ <- switch(nx, n = "\n", t = "\t", r = "\r", "\\" = "\\", "'" = "'", "\"" = "\"", NULL)
      if (!is.null(rep_)) { out <- c(out, rep_); i <- i + 2L; next }
      out <- c(out, ch); i <- i + 1L; next
    }
    out <- c(out, ch); i <- i + 1L
  }
  paste(out, collapse = "")
}
