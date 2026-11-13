# Helpers for the tests built from the SymPy comparison battery: .pyrun()
# evaluates a call and formats its value the way Python prints the SymPy
# result (str() at the top level and repr() inside containers, with ** shown
# as ^), or an error as "ERR <class>: <message>".

.pyfloat <- function(x) {
  if (is.na(x)) return("nan")
  if (is.infinite(x)) return(if (x > 0) "inf" else "-inf")
  for (d in 15:17) {
    s <- formatC(x, digits = d, format = "g")
    if (as.numeric(s) == x) break
  }
  s <- trimws(s)
  if (grepl("e", s)) {
    m <- regmatches(s, regexec("^(-?[0-9.]+)e([+-])([0-9]+)$", s))[[1]]
    ex <- as.integer(m[4])
    if (m[3] == "+" && ex < 16) s <- format(x, scientific = FALSE, digits = 17)
    else s <- sprintf("%se%s%02d", m[2], m[3], ex)
  }
  if (!grepl("[.e]|inf|nan", s)) s <- paste0(s, ".0")
  s
}

# Python sorts the strings with ** where R prints ^
.pyord <- function(s) order(gsub("^", "**", s, fixed = TRUE), method = "radix")

.pyfmt <- function(v, top = TRUE) {
  if (is.null(v)) return("None")
  if (inherits(v, "pytuple")) {
    s <- vapply(v, .pyfmt, "", top = FALSE)
    return(paste0("(", paste(s, collapse = ", "), if (length(s) == 1L) ",", ")"))
  }
  if (inherits(v, "pydict")) {
    if (!length(v$keys)) return("{}")
    return(paste0("{", paste(vapply(seq_along(v$keys), function(i)
      paste0(.pyfmt(v$keys[[i]], FALSE), ": ", .pyfmt(v$values[[i]], FALSE)), ""), collapse = ", "), "}"))
  }
  if (isTRUE(attr(v, "pyset"))) {
    if (!length(v)) return("set()")
    s <- vapply(v, .pyfmt, "", top = FALSE)
    return(paste0("{", paste(s[.pyord(s)], collapse = ", "), "}"))
  }
  if (inherits(v, "sym") || inherits(v, "symmatrix")) return(sstr(v))
  if (inherits(v, "Poly")) return(paste(format(v), collapse = " "))
  if (inherits(v, "bigz") || inherits(v, "bigq")) {
    if (length(v) == 1L) return(as.character(v))
    return(paste0("[", paste(as.character(v), collapse = ", "), "]"))
  }
  if (is.character(v) && length(v) == 1L && is.null(names(v))) return(if (top) v else paste0("'", v, "'"))
  if (is.logical(v) && length(v) == 1L && is.null(names(v)))
    return(if (is.na(v)) "None" else if (v) "True" else "False")
  if (is.integer(v) && length(v) == 1L && is.null(names(v))) return(as.character(v))
  if (is.double(v) && length(v) == 1L && is.null(names(v))) return(.pyfloat(v))
  if (is.object(v)) return(paste(format(v, repr = !top), collapse = "\n"))
  if (is.atomic(v)) v <- as.list(v)
  if (is.list(v)) {
    nm <- names(v)
    if (length(v) && !is.null(nm) && all(nzchar(nm)))
      return(paste0("{", paste(vapply(seq_along(v), function(i) paste0(nm[i], ": ", .pyfmt(v[[i]], FALSE)), ""),
                               collapse = ", "), "}"))
    return(paste0("[", paste(vapply(v, .pyfmt, "", top = FALSE), collapse = ", "), "]"))
  }
  paste(format(v), collapse = " ")
}

.pyL <- function(x) as.list(x)
.pyT <- function(x) structure(as.list(x), class = "pytuple")
.pytuple <- function(...) structure(list(...), class = "pytuple")
.pyset <- function(x) structure(as.list(x), pyset = TRUE)
.pyD <- function(x) structure(list(keys = as.list(names(x)), values = unname(as.list(x))), class = "pydict")
# a dict printed with its keys sorted
.pySD <- function(d) {
  if (!is.list(d)) return(d)
  if (!length(d)) return("{}")
  nm <- names(d)
  o <- .pyord(nm)
  paste0("{", paste0(nm[o], ": ", vapply(d[o], .pyfmt, "", top = FALSE), collapse = ", "), "}")
}
# a list of models (dicts) as sorted strings
.pySMS <- function(ms) {
  s <- vapply(ms, function(m) if (is.list(m)) .pySD(m) else .pyfmt(m), "")
  as.list(s[.pyord(s)])
}

.pyerr_class <- function(e) {
  cl <- setdiff(class(e), c("simpleError", "error", "condition", "sympy_error"))
  if (length(cl)) cl[1L] else "Error"
}
# the result, or only the class of the error raised
.pyEC <- function(expr) tryCatch(expr, error = function(e) paste0("ERR ", .pyerr_class(e)))

.pyrun <- function(expr) {
  res <- tryCatch(gsub("**", "^", .pyfmt(suppressWarnings(expr)), fixed = TRUE),
                  error = function(e) gsub("**", "^", paste0("ERR ", .pyerr_class(e), ": ", conditionMessage(e)), fixed = TRUE))
  gsub("\n", "\\n", res, fixed = TRUE)
}
