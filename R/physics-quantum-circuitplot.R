# CircuitPlot of sympy/physics/quantum/circuitplot.py, drawn with base R
# graphics instead of matplotlib. The plot_gate methods of the gate classes
# (gate.py) are the branches of .cp_plot_gate().

#' Plot quantum circuits
#'
#' `circuit_plot(c, nqubits)` (also `CircuitPlot()`) draws a circuit, a
#' gate or a product of gates, on `nqubits` wires: one qubit gates as boxes,
#' controlled gates as a control line with control points (with
#' `simplify_cgate`, a controlled `X` is drawn as a not point and a
#' controlled `Z` as a control point), swaps as crosses and measured wires
#' doubled after the measurement. The gates act from right to left, so the
#' rightmost factor of a product is drawn first. Options are SymPy's
#' attributes: `scale`, `fontsize` (here a `cex` multiplier, default 1),
#' `linewidth`, `control_radius`, `not_radius`, `swap_delta`, `labels`
#' (wire labels, see [labeller()]), `inits` (a named list of initial
#' states) and `label_buffer`.
#'
#' @param c A gate or a product of gates.
#' @param nqubits The number of wires.
#' @param ... Options, see Details.
#' @return Invisibly, a `CircuitPlot` object with the grids and gates.
#' @examples
#' circuit_plot(CNOT(1, 0) * HadamardGate(1), 2, labels = labeller(2))
#' circuit_plot(SwapGate(0, 1) * Mz(0), 2)
#' @name circuit_plot
#' @export
circuit_plot <- function(c, nqubits, ...) {
  opts <- utils::modifyList(list(scale = 1, fontsize = 1, linewidth = 1, control_radius = 0.05,
                                 not_radius = 0.15, swap_delta = 0.05, labels = character(),
                                 inits = list(), label_buffer = 0.5), list(...))
  self <- c(list(circuit = c, nqubits = as.integer(nqubits)), opts)
  self$gates <- .cp_gates(c)
  self$ngates <- length(self$gates)
  if (!self$ngates) stop("ValueError: the circuit has no gates", call. = FALSE)
  self$wire_grid <- (seq_len(self$nqubits) - 1) * self$scale
  self$gate_grid <- (seq_len(self$ngates) - 1) * self$scale
  self <- structure(self, class = "CircuitPlot")
  off <- 0.5 * self$scale
  xl <- c(self$gate_grid[1L] - off - self$scale - if (length(self$labels)) self$label_buffer + 0.25 else 0,
          self$gate_grid[self$ngates] + off + self$scale)
  op <- graphics::par(mar = c(0, 0, 0, 0)); on.exit(graphics::par(op))
  graphics::plot.new()
  graphics::plot.window(xl, c(self$wire_grid[1L] - off, self$wire_grid[self$nqubits] + off), asp = 1)
  .cp_plot_wires(self)
  for (i in seq_len(self$ngates)) .cp_plot_gate(self, self$gates[[i]], i)
  invisible(self)
}

#' @rdname circuit_plot
#' @export
CircuitPlot <- circuit_plot

#' @export
print.CircuitPlot <- function(x, ...) {
  cat(sprintf("CircuitPlot of %s on %d qubits\n", sstr(x$circuit), x$nqubits)); invisible(x)
}

.cp_gates <- function(c) {
  if (is_sym(c) && c$type == "Mul") return(Filter(function(g) .qisa(g, "Gate"), rev(c$args)))
  if (.qisa(c, "Gate")) return(list(c))
  list()
}

.cp_int <- function(v) as.integer(.qint(v))

.cp_measurements <- function(self) {
  m <- list()
  for (i in seq_along(self$gates)) {
    g <- self$gates[[i]]
    if (isTRUE(.qtry_ni(tryCatch(.qa(g, "measurement"), error = function(e) FALSE), FALSE)))
      for (t in .qa(g, "targets")) {
        k <- as.character(.cp_int(t))
        if (is.null(m[[k]]) || m[[k]] > i) m[[k]] <- i
      }
  }
  m
}

.cp_line <- function(self, x, y) graphics::lines(x, y, lwd = self$linewidth)

.cp_plot_wires <- function(self) {
  xdata <- c(self$gate_grid[1L] - self$scale, self$gate_grid[self$ngates] + self$scale)
  for (i in seq_len(self$nqubits)) {
    y <- self$wire_grid[i]
    .cp_line(self, xdata, c(y, y))
    if (length(self$labels)) {
      lab <- self$labels[i]; init <- self$inits[[lab]]
      buf <- if (!is.null(init)) 0.25 else 0
      txt <- if (!is.null(init)) sprintf("|%s> = |%s>", lab, init) else sprintf("|%s>", lab)
      graphics::text(xdata[1L] - self$label_buffer - buf, y, txt, cex = self$fontsize)
    }
  }
  # doubled wires after a measurement, and doubled control lines off them
  ism <- .cp_measurements(self); dy <- 0.04
  xstop <- self$gate_grid[self$ngates]
  for (w in names(ism)) {
    y <- self$wire_grid[as.integer(w) + 1L] + dy
    .cp_line(self, c(self$gate_grid[ism[[w]]], xstop + self$scale), c(y, y))
  }
  for (i in seq_along(self$gates)) {
    g <- self$gates[[i]]
    if (.qisa(g, c("CGate", "CGateS"))) {
      wires <- vapply(c(.qa(g, "controls"), .qa(g, "targets")), .cp_int, 0L)
      for (w in wires) {
        k <- as.character(w)
        if (!is.null(ism[[k]]) && self$gate_grid[i] > self$gate_grid[ism[[k]]])
          .cp_line(self, rep(self$gate_grid[i] - dy, 2), self$wire_grid[range(wires) + 1L])
      }
    }
  }
}

.cp_box <- function(self, t, gi, wi) {
  x <- self$gate_grid[gi]; y <- self$wire_grid[wi + 1L]
  h <- 0.2 * self$scale
  w <- max(h, graphics::strwidth(t, cex = self$fontsize) / 2 + 0.05)
  graphics::rect(x - w, y - h, x + w, y + h, col = "white", border = "black", lwd = self$linewidth)
  graphics::text(x, y, t, cex = self$fontsize)
}
.cp_control_line <- function(self, gi, lo, hi) .cp_line(self, rep(self$gate_grid[gi], 2), self$wire_grid[c(lo, hi) + 1L])
.cp_control_point <- function(self, gi, wi)
  graphics::symbols(self$gate_grid[gi], self$wire_grid[wi + 1L], circles = self$control_radius * self$scale,
                    inches = FALSE, add = TRUE, bg = "black", fg = "black", lwd = self$linewidth)
.cp_not_point <- function(self, gi, wi) {
  x <- self$gate_grid[gi]; y <- self$wire_grid[wi + 1L]; r <- self$not_radius
  graphics::symbols(x, y, circles = r, inches = FALSE, add = TRUE, bg = "white", fg = "black", lwd = self$linewidth)
  .cp_line(self, c(x, x), c(y - r, y + r))
}
.cp_swap_point <- function(self, gi, wi) {
  x <- self$gate_grid[gi]; y <- self$wire_grid[wi + 1L]; d <- self$swap_delta
  .cp_line(self, c(x - d, x + d), c(y - d, y + d))
  .cp_line(self, c(x - d, x + d), c(y + d, y - d))
}

.cp_plot_gate <- function(self, g, gi) {
  if (.qisa(g, "CGate")) {
    ctl <- vapply(.qa(g, "controls"), .cp_int, 0L)
    tg <- vapply(.qa(g, "targets"), .cp_int, 0L)
    .cp_control_line(self, gi, min(ctl, tg), max(ctl, tg))
    for (c in ctl) .cp_control_point(self, gi, c)
    inner <- .qa(g, "gate")
    if (isTRUE(.qa(g, "simplify_cgate")) && .qa(inner, "gate_name") == "X") .cp_not_point(self, gi, tg[1L])
    else if (isTRUE(.qa(g, "simplify_cgate")) && .qa(inner, "gate_name") == "Z") .cp_control_point(self, gi, tg[1L])
    else .cp_plot_gate(self, inner, gi)
  } else if (.qisa(g, "SwapGate")) {
    tg <- vapply(.qa(g, "targets"), .cp_int, 0L)
    .cp_control_line(self, gi, min(tg), max(tg))
    .cp_swap_point(self, gi, min(tg)); .cp_swap_point(self, gi, max(tg))
  } else {
    .cp_box(self, .qa(g, "gate_name"), gi, .cp_int(.qa(g, "targets")[[1L]]))
  }
}
