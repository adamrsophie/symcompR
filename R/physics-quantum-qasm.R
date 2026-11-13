# QASM circuit descriptions (sympy/physics/quantum/qasm.py), the mock
# measurement gates and gate factories of circuitplot.py, and Shor's
# algorithm (sympy/physics/quantum/shor.py).

.qclass("Mz", "OneQubitGate", measurement = TRUE, gate_name = "Mz", gate_name_latex = "M_z")
.qclass("Mx", "OneQubitGate", measurement = TRUE, gate_name = "Mx", gate_name_latex = "M_x")

#' Measurement gates and gate factories
#'
#' `Mz(target)` and `Mx(target)` are the mock measurement gates of
#' `sympy.physics.quantum.circuitplot`. `CreateOneQubitGate(name)` makes a
#' new one-qubit gate class and returns its constructor, and
#' `CreateCGate(name)` returns a function `(controls, target)` building a
#' controlled gate of a new one-qubit gate. `labeller(n)` gives wire labels.
#' Circuit plotting itself needs matplotlib and is not available.
#'
#' @param ... The target qubit.
#' @param name,latexname The gate name and its LaTeX form.
#' @param n The number of qubits.
#' @param symbol The label prefix.
#' @return A gate, a gate constructor or a character vector.
#' @examples
#' Mz(0)
#' Q <- CreateOneQubitGate("Q")
#' Q(1)
#' labeller(3)
#' @name circuit-gates
NULL

#' @rdname circuit-gates
#' @export
Mz <- function(...) .qfunc("Mz", list(...))
#' @rdname circuit-gates
#' @export
Mx <- function(...) .qfunc("Mx", list(...))

#' @rdname circuit-gates
#' @export
CreateOneQubitGate <- function(name, latexname = NULL) {
  if (is.null(latexname) || !nzchar(latexname)) latexname <- name
  cls <- paste0(name, "Gate")
  .qclass(cls, "OneQubitGate", gate_name = name, gate_name_latex = latexname)
  function(...) .qfunc(cls, list(...))
}

#' @rdname circuit-gates
#' @export
CreateCGate <- function(name, latexname = NULL) {
  if (is.null(latexname) || !nzchar(latexname)) latexname <- name
  onequbitgate <- CreateOneQubitGate(name, latexname)
  function(ctrls, target) CGate(do.call(Tuple, as.list(ctrls)), onequbitgate(target))
}

#' @rdname circuit-gates
#' @export
labeller <- function(n, symbol = "q") sprintf("%s_%d", symbol, n - seq_len(n))

# ---------------------------------------------------------------------------
# QASM

flip_index <- function(i, n) n - i - 1L

trim <- function(line) {
  if (!grepl("#", line, fixed = TRUE)) return(line)
  strsplit(line, "#", fixed = TRUE)[[1L]][1L]
}

get_index <- function(target, labels) flip_index(match(target, labels) - 1L, length(labels))

.qasm_get_indices <- function(targets, labels) vapply(unlist(targets), get_index, 0, labels = labels, USE.NAMES = FALSE)

.qisspace <- function(s) nzchar(s) && grepl("^[[:space:]]+$", s)

# like Python, a single string is iterated character by character
nonblank <- function(args) {
  if (is.character(args) && length(args) == 1L) args <- as.list(strsplit(args, "")[[1L]])
  out <- list()
  for (line in args) {
    line <- trim(line)
    if (.qisspace(line)) next
    out[[length(out) + 1L]] <- line
  }
  out
}

fixcommand <- function(c) {
  c <- gsub("-", "", tolower(c), fixed = TRUE)
  if (c == "def") return("qdef")
  c
}

fullsplit <- function(line) {
  words <- strsplit(trimws(line), "[[:space:]]+")[[1L]]
  rest <- paste(words[-1L], collapse = " ")
  parts <- if (nzchar(rest)) strsplit(rest, ",", fixed = TRUE)[[1L]] else ""
  if (endsWith(rest, ",")) parts <- c(parts, "")
  list(fixcommand(words[1L]), as.list(trimws(parts)))
}

stripquotes <- function(s) gsub("'", "", gsub("\"", "", s, fixed = TRUE), fixed = TRUE)

#' QASM circuits
#'
#' `Qasm(...)` reads lines of QASM (`qubit`, `h`, `x`, `z`, `s`, `t`,
#' `cnot`, `swap`, `cphase`, `toffoli`, `c-x`, `c-z`, `measure`, `def`) as
#' in `sympy.physics.quantum.qasm`. The result is an environment whose
#' methods are called with `$`, for example `q$get_circuit()` and
#' `q$get_labels()`. `read_qasm()` reads a string and `read_qasm_file()` a
#' file.
#'
#' @param ... Lines of QASM.
#' @param lines A string of QASM lines.
#' @param filename A file of QASM lines.
#' @return A `Qasm` object.
#' @examples
#' q <- Qasm("qubit q0", "qubit q1", "h q0", "cnot q0,q1")
#' q$get_circuit()
#' @name qasm
NULL

#' @rdname qasm
#' @export
Qasm <- function(...) {
  self <- new.env()
  self$defs <- list(); self$circuit <- list(); self$labels <- character(); self$inits <- list()
  push <- function(g) self$circuit[[length(self$circuit) + 1L]] <- g
  self$add <- function(...) {
    for (line in nonblank(as.list(c(...)))) {
      fs <- fullsplit(line)
      command <- fs[[1L]]; rest <- fs[[2L]]
      if (!is.null(self$defs[[command]])) {
        f <- self$defs[[command]]
        indices <- self$indices(rest)
        if (length(indices) == 1L) push(f(indices[1L])) else push(f(indices[-length(indices)], indices[length(indices)]))
      } else if (command %in% .qasm_methods) {
        do.call(self[[command]], as.list(rest))
      } else {
        message(sprintf("Function %s not defined. Skipping", command))
      }
    }
    invisible(NULL)
  }
  self$get_circuit <- function() if (!length(self$circuit)) S_One else .mul(rev(self$circuit))
  self$get_labels <- function() rev(self$labels)
  self$plot <- function() .qni("Circuit plotting needs matplotlib and is not available in R.")
  self$qubit <- function(arg, init = NULL) {
    self$labels <- c(self$labels, arg)
    if (length(init) && nzchar(init)) self$inits[[arg]] <- init
    invisible(NULL)
  }
  self$indices <- function(args) .qasm_get_indices(args, self$labels)
  self$index <- function(arg) get_index(arg, self$labels)
  self$nop <- function(...) invisible(NULL)
  self$x <- function(arg) push(XGate(self$index(arg)))
  self$z <- function(arg) push(ZGate(self$index(arg)))
  self$h <- function(arg) push(HadamardGate(self$index(arg)))
  self$s <- function(arg) push(PhaseGate(self$index(arg)))
  self$t <- function(arg) push(TGate(self$index(arg)))
  self$measure <- function(arg) push(Mz(self$index(arg)))
  self$cnot <- function(a1, a2) push(do.call(CNOT, as.list(self$indices(c(a1, a2)))))
  self$swap <- function(a1, a2) push(do.call(SWAP, as.list(self$indices(c(a1, a2)))))
  self$cphase <- function(a1, a2) push(do.call(CPHASE, as.list(self$indices(c(a1, a2)))))
  self$toffoli <- function(a1, a2, a3) {
    i <- self$indices(c(a1, a2, a3))
    push(CGateS(list(i[1L], i[2L]), XGate(i[3L])))
  }
  self$cx <- function(a1, a2) { i <- self$indices(c(a1, a2)); push(CGate(i[1L], XGate(i[2L]))) }
  self$cz <- function(a1, a2) { i <- self$indices(c(a1, a2)); push(CGate(i[1L], ZGate(i[2L]))) }
  self$defbox <- function(...) message("defbox not supported yet. Skipping: ", paste(c(...), collapse = ", "))
  self$qdef <- function(name, ncontrols, symbol) {
    ncontrols <- as.integer(ncontrols)
    command <- fixcommand(name)
    symbol <- stripquotes(symbol)
    self$defs[[command]] <- if (ncontrols > 0) CreateCGate(symbol) else CreateOneQubitGate(symbol)
    invisible(NULL)
  }
  self$add(...)
  class(self) <- "Qasm"
  self
}

.qasm_methods <- c("add", "get_circuit", "get_labels", "plot", "qubit", "indices", "index", "nop", "x", "z", "h", "s",
                   "t", "measure", "cnot", "swap", "cphase", "toffoli", "cx", "cz", "defbox", "qdef")

#' @export
print.Qasm <- function(x, ...) {
  cat("<Qasm circuit on", length(x$labels), "qubits>\n")
  invisible(x)
}

#' @rdname qasm
#' @export
read_qasm <- function(lines) do.call(Qasm, as.list(strsplit(lines, "\r\n|\n|\r")[[1L]]))

#' @rdname qasm
#' @export
read_qasm_file <- function(filename) do.call(Qasm, as.list(readLines(filename)))

# ---------------------------------------------------------------------------
# Shor's algorithm

.qclass("CMod", "Gate",
  `_eval_args` = .cm(function(cls, args, ...) .qni("The CMod gate has not been completed.")),
  t = .prop(function(self) .qa(self, "label")[[1L]]),
  a = .prop(function(self) .qa(self, "label")[[2L]]),
  N = .prop(function(self) .qa(self, "label")[[3L]]),
  `_apply_operator_Qubit` = function(self, qubits, ...) {
    t <- .qint(.qa(self, "t"))
    n <- 1; k <- 0
    for (i in seq_len(t) - 1L) {
      k <- k + n * .qint(.qmeth(qubits, "__getitem__", t + i))
      n <- n * 2
    }
    out <- .qint(.qmod_int(.pow(.qa(self, "a"), as_sym(k)), .qa(self, "N")))
    outarray <- qubits$args[seq_len(t)]
    for (i in rev(seq_len(t) - 1L)) outarray[[length(outarray) + 1L]] <- as_sym(bitwAnd(bitwShiftR(out, i), 1L))
    do.call(Qubit, outarray)
  })

.qmod_int <- function(a, b) as_sym(gmp::mod.bigz(gmp::as.bigz(.qint_big(a)), gmp::as.bigz(.qint_big(b))))
.qint_big <- function(x) gmp::numerator(as_sym(x)$value)

#' Shor's algorithm
#'
#' `shor(N)` factors `N` with the period finding circuit of
#' `sympy.physics.quantum.shor`; `period_find()` runs that circuit,
#' `getr()` and `ratioize()` turn a measurement into a period through its
#' continued fraction. The modular exponentiation gate `CMod` is not
#' completed in SymPy either, so `shor()` and `period_find()` raise a
#' `NotImplementedError`.
#'
#' @param N The number to factor, or the bound on the period.
#' @param a The base.
#' @param ... The arguments of the `CMod` gate.
#' @param x,y The numerator and denominator of the measured fraction.
#' @param list The continued fraction terms.
#' @return A factor, a period or an expression.
#' @examples
#' getr(513, 1024, 10)
#' @name shor
NULL

#' @rdname shor
#' @export
CMod <- function(...) .qfunc("CMod", list(...))

#' @rdname shor
#' @export
shor <- function(N) {
  N <- .qint(N)
  a <- sample.int(N - 2L, 1L) - 1L + 2L
  g <- .qint(igcd(N, a))
  if (g != 1) return(as_sym(g))
  r <- period_find(a, N)
  if (.qint(r) %% 2 == 1) shor(N)
  half <- .mul(list(r, S_Half))
  list(igcd(.sub(.pow(as_sym(a), half), S_One), N), igcd(.add(list(.pow(as_sym(a), half), S_One)), N))
}

#' @rdname shor
#' @export
getr <- function(x, y, N) ratioize(continued_fraction_periodic(x, y), N)

#' @rdname shor
#' @export
ratioize <- function(list, N) {
  first <- as_sym(list[[1L]])
  if (isTRUE(.qgt(first, as_sym(N)))) return(S_Zero)
  if (length(list) == 1L) return(first)
  .add(list(first, ratioize(list[-1L], N)))
}

#' @rdname shor
#' @export
period_find <- function(a, N) {
  epsilon <- 0.5
  t <- as.integer(2 * ceiling(log(.qint(N), 2)))
  start <- rep(list(0L), t)
  factor <- .pow(.pow(S_Two, as_sym(t)), .neg(S_Half))
  qubits <- S_Zero
  for (arr in variations(list(0L, 1L), t, repetition = TRUE)) qubits <- .add(list(qubits, do.call(Qubit, c(as.list(arr), start))))
  circuit <- expand(.mul(list(factor, qubits)))
  circuit <- .mul(list(.qfunc("CMod", list(t, a, N)), circuit))
  circuit <- qapply(circuit)
  for (i in seq_len(t) - 1L) circuit <- measure_partial_oneshot(circuit, i)
  circuit <- qapply(.mul(list(.qmeth(QFT(t, t * 2L), "decompose"), circuit)), floatingPoint = TRUE)
  for (i in seq_len(t) - 1L) circuit <- measure_partial_oneshot(circuit, i + t)
  register <- if (.qisa(circuit, "Qubit")) circuit
              else if (circuit$type == "Mul") .qargs(circuit)[[length(.qargs(circuit))]]
              else { l <- .qargs(circuit)[[length(.qargs(circuit))]]; .qargs(l)[[length(.qargs(l))]] }
  n <- 1; answer <- 0
  for (i in seq_len(length(register$args) %/% 2) - 1L) {
    answer <- answer + n * .qint(.qmeth(register, "__getitem__", i + t))
    n <- n * 2
  }
  if (answer == 0) .qstop("OrderFindingException", sprintf("Order finder returned 0. Happens with chance %f", epsilon))
  getr(answer, 2^t, N)
}
